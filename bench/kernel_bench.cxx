/*
  SeWaS
  Copyright (C) 2018  ANEO

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published
  by the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
==============================================================================*/

/* Isolated-kernel benchmark for the 4th-order staggered derivative (Phase 1 of
 * sewas_transformation_plan_en.md).
 *
 * Compares two data layouts on the same bundle of P synthetic tiles, P being the SIMD width
 * (SWS::packetSize):
 *
 *  - std: today's SpatialBlockField layout (one tile per array, z contiguous). Dx/Dy go
 *    through CentralFDOperator::apply<d> ("eigen"), Dz through applyZ<Legacy|Unaligned>.
 *    "raw" runs the generic DLI loop body below with scalar elements on each std tile, as a
 *    control separating the layout effect from Eigen's per-(i,j)-column overhead.
 *  - dli: Data Layout Interleaving, as in Legolas++ (github.com/LaurentPlagne/Legolas). The P
 *    tiles are interleaved cell by cell, so element (i,j,k) is a SIMD pack holding that cell
 *    for all P tiles. The kernel is plain scalar-looking code over packs, identical for the
 *    three directions up to the tap stride; every tap is one aligned vector load, so Dz loses
 *    its unaligned/shuffle penalty.
 *
 * Every variant writes its derivative into an output field (no reduction in the timed loop),
 * and the DLI results are checked lane by lane against the std results before timing. */

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CentralFDOperator.hxx"
#include "Config.hxx"

namespace {

using Real = SWS::RealType;
using Tile = SWS::SpatialBlockField<Real>;

constexpr int P = SWS::packetSize;

// GCC/Clang vector extension: arithmetic on Pack compiles to native SIMD (same mechanism as
// Legolas::NativeSimd). Scalar operands are broadcast.
typedef Real Pack __attribute__((vector_size(P * sizeof(Real))));

struct Args
{
  int size = 50;
  int reps = 100;
};

Args
parseArgs(int argc, char** argv)
{
  Args args;

  for (int i = 1; i < argc; i++) {
    const std::string arg = argv[i];
    if (arg == "--size" && i + 1 < argc) {
      args.size = std::atoi(argv[++i]);
    } else if (arg == "--reps" && i + 1 < argc) {
      args.reps = std::atoi(argv[++i]);
    } else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: kernel_bench [--size N] [--reps N]\n"
                << "  --size N  interior edge length of each synthetic cubic tile (default 50)\n"
                << "  --reps N  number of sweeps over the tile bundle per variant (default 100)\n";
      exit(0);
    }
  }

  return args;
}

/* P tiles interleaved cell by cell, Y_MAJOR order like SpatialBlockField. */
struct DliField
{
  int cx, cy, cz;
  std::vector<Pack> data; // C++17 aligned new honors alignof(Pack)

  DliField(const int nx, const int ny, const int nz)
    : cx(nx)
    , cy(ny)
    , cz(nz)
    , data(static_cast<size_t>(nx) * ny * nz, Pack{})
  {}

  inline size_t index(const int i, const int j, const int k) const
  {
    return (static_cast<size_t>(i) * cy + j) * cz + k;
  }
};

struct Halo
{
  int hx, hy, hz;
};

/* Generic tile kernel: one body for all directions, d only selects the tap stride (and DK
 * only shifts Z, as in CentralFDOperator::applyZ). Coefficients/offsets match
 * CentralFDOperator: c1*f[+d1] + c2*f[-d2] + c3*f[+d3] + c4*f[-d4], d = (1, 2, 0, 1).
 *
 * T = Pack runs it on a DliField (P tiles per call, explicit SIMD across tiles); T = Real
 * runs the very same loop on one std-layout tile ("raw", SIMD left to the compiler's
 * auto-vectorizer along k), which separates the layout effect from Eigen's per-column
 * overhead in the std/eigen variants. */
template<SWS::Directions d, typename T>
void
derivative(const T* __restrict in, T* __restrict out, const int cx, const int cy, const int cz, const Halo& h,
           const Real ds, const short DK)
{
  const Real c1 = (-1. / 24) / ds, c2 = (+1. / 24) / ds, c3 = (+9. / 8) / ds, c4 = (-9. / 8) / ds;

  const long s = (d == SWS::X) ? long(cy) * cz : (d == SWS::Y) ? long(cz) : 1L;
  const int dk = (d == SWS::Z) ? DK : 0;
  const long o1 = (1 + dk) * s, o2 = -(2 - dk) * s, o3 = (0 + dk) * s, o4 = -(1 - dk) * s;

  for (int i = h.hx; i < cx - h.hx; i++) {
    for (int j = h.hy; j < cy - h.hy; j++) {
      const T* __restrict src = in + (long(i) * cy + j) * cz;
      T* __restrict dst = out + (long(i) * cy + j) * cz;
      for (int k = h.hz; k < cz - h.hz; k++) {
        dst[k] = c1 * src[k + o1] + c2 * src[k + o2] + c3 * src[k + o3] + c4 * src[k + o4];
      }
    }
  }
}

template<SWS::Directions d>
void
dliDerivative(const DliField& f, DliField& out, const Halo& h, const Real ds, const short DK)
{
  derivative<d>(f.data.data(), out.data.data(), f.cx, f.cy, f.cz, h, ds, DK);
}

template<SWS::Directions d>
void
rawDerivative(const Tile& f, Tile& out, const Halo& h, const Real ds, const short DK)
{
  derivative<d>(f.data().data(), out.data().data(), int(f.dimension(SWS::X)), int(f.dimension(SWS::Y)),
                int(f.dimension(SWS::Z)), h, ds, DK);
}

double
timeIt(const int reps, const std::function<void()>& sweep)
{
  sweep(); // warm-up: page in outputs, fill caches

  const auto start = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; r++) {
    sweep();
  }
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

/* Output checksums are printed so the compiler can't prove the stores dead. */
Real
checksum(const std::vector<std::unique_ptr<Tile>>& tiles)
{
  Real s = 0;
  for (const auto& t : tiles) {
    s += t->data().sum();
  }
  return s;
}

Real
checksum(const DliField& f)
{
  Real s = 0;
  for (const auto& p : f.data) {
    for (int l = 0; l < P; l++) {
      s += p[l];
    }
  }
  return s;
}

/* Max |dli - std| / max |std| over the interior, all lanes. */
double
compare(const DliField& dli, const std::vector<std::unique_ptr<Tile>>& ref, const Halo& h)
{
  double maxDiff = 0, maxRef = 0;
  for (int p = 0; p < P; p++) {
    for (int i = h.hx; i < dli.cx - h.hx; i++) {
      for (int j = h.hy; j < dli.cy - h.hy; j++) {
        for (int k = h.hz; k < dli.cz - h.hz; k++) {
          const double r = (*ref[p])(i, j, k);
          maxDiff = std::max(maxDiff, std::abs(double(dli.data[dli.index(i, j, k)][p]) - r));
          maxRef = std::max(maxRef, std::abs(r));
        }
      }
    }
  }
  return maxRef > 0 ? maxDiff / maxRef : maxDiff;
}

void
report(const std::string& layout, const char dir, const std::string& variant, const double t,
       const double cells, const Real sum)
{
  std::cout << std::left << std::setw(5) << layout << std::setw(5) << dir << std::setw(11) << variant
            << std::right << std::fixed << std::setprecision(4) << std::setw(10) << t << std::setw(12)
            << std::setprecision(3) << cells / t * 1e-9 << std::defaultfloat << "\n";
  std::cerr << "  (" << layout << " D" << dir << " " << variant << " checksum: " << std::setprecision(17) << sum
            << std::setprecision(6) << ")\n";
}

} // namespace

int
main(int argc, char** argv)
{
  const auto args = parseArgs(argc, argv);

  const Halo h{ CentralFDOperator::hnx(), CentralFDOperator::hny(), CentralFDOperator::hnz() };

  const int cx = h.hx + args.size + h.hx;
  const int cy = h.hy + args.size + h.hy;
  const int cz = h.hz + args.size + h.hz;

  const Real ds = 0.1;
  CentralFDOperator fdo(ds, ds, ds);

  // P distinct tiles in the std layout, and the same data interleaved in a DliField.
  std::vector<std::unique_ptr<Tile>> in, out;
  DliField dliIn(cx, cy, cz), dliOut(cx, cy, cz);

  for (int p = 0; p < P; p++) {
    in.emplace_back(std::make_unique<Tile>(cx, cy, cz));
    out.emplace_back(std::make_unique<Tile>(cx, cy, cz));
    in[p]->setHaloSize(h.hx, h.hy, h.hz);
    out[p]->setHaloSize(h.hx, h.hy, h.hz);
    *out[p] = 0.0;

    for (int i = 0; i < cx; i++) {
      for (int j = 0; j < cy; j++) {
        for (int k = 0; k < cz; k++) {
          // Non-polynomial in every direction (the 4th-order stencil is exact up to degree
          // 3, so a low-degree field would hide wrong taps/coefficients) and bounded, so the
          // tap sum doesn't cancel catastrophically and the lane check below stays tight; p
          // shifts the phase so lanes differ.
          const auto v = static_cast<Real>(std::sin(0.37 * k + 0.11 * i + 0.5 * p) * std::cos(0.23 * j + 0.3 * p) +
                                           0.5 * std::sin(0.19 * i + 0.29 * j + 0.41 * k + 0.7 * p));
          (*in[p])(i, j, k) = v;
          dliIn.data[dliIn.index(i, j, k)][p] = v;
        }
      }
    }
  }

  const int nI = cx - 2 * h.hx, nJ = cy - 2 * h.hy, nK = cz - 2 * h.hz;
  const double cells = double(args.reps) * P * nI * nJ * nK;
  const double footprintMiB = 2.0 * P * cx * cy * cz * sizeof(Real) / (1024. * 1024.);

  std::cout << "kernel_bench: " << P << " tiles of " << cx << "x" << cy << "x" << cz << " (interior " << nI << "x"
            << nJ << "x" << nK << "), " << (sizeof(Real) == 4 ? "fp32" : "fp64") << ", " << args.reps
            << " reps, in+out footprint " << std::setprecision(3) << footprintMiB << " MiB\n";

  // Correctness: DLI vs std (Unaligned for Z, apply<d> for X/Y), all lanes, before timing.
  auto stdX = [&] {
    for (int p = 0; p < P; p++)
      for (int i = h.hx; i < cx - h.hx; i++)
        for (int j = h.hy; j < cy - h.hy; j++)
          (*out[p])(i, j) = fdo.apply<SWS::X>(*in[p], i, j, SWS::DK);
  };
  auto stdY = [&] {
    for (int p = 0; p < P; p++)
      for (int i = h.hx; i < cx - h.hx; i++)
        for (int j = h.hy; j < cy - h.hy; j++)
          (*out[p])(i, j) = fdo.apply<SWS::Y>(*in[p], i, j, SWS::DK);
  };
  auto stdZLegacy = [&] {
    for (int p = 0; p < P; p++)
      for (int i = h.hx; i < cx - h.hx; i++)
        for (int j = h.hy; j < cy - h.hy; j++)
          (*out[p])(i, j) = fdo.applyZ<SWS::KernelVariant::Legacy>(*in[p], i, j, SWS::DK);
  };
  auto stdZUnaligned = [&] {
    for (int p = 0; p < P; p++)
      for (int i = h.hx; i < cx - h.hx; i++)
        for (int j = h.hy; j < cy - h.hy; j++)
          (*out[p])(i, j) = fdo.applyZ<SWS::KernelVariant::Unaligned>(*in[p], i, j, SWS::DK);
  };
  auto rawX = [&] {
    for (int p = 0; p < P; p++)
      rawDerivative<SWS::X>(*in[p], *out[p], h, ds, SWS::DK);
  };
  auto rawY = [&] {
    for (int p = 0; p < P; p++)
      rawDerivative<SWS::Y>(*in[p], *out[p], h, ds, SWS::DK);
  };
  auto rawZ = [&] {
    for (int p = 0; p < P; p++)
      rawDerivative<SWS::Z>(*in[p], *out[p], h, ds, SWS::DK);
  };
  auto dliX = [&] { dliDerivative<SWS::X>(dliIn, dliOut, h, ds, SWS::DK); };
  auto dliY = [&] { dliDerivative<SWS::Y>(dliIn, dliOut, h, ds, SWS::DK); };
  auto dliZ = [&] { dliDerivative<SWS::Z>(dliIn, dliOut, h, ds, SWS::DK); };

  // For each direction: reference = std/eigen output, then check raw and dli against it.
  std::vector<std::unique_ptr<Tile>> ref;
  for (int p = 0; p < P; p++) {
    ref.emplace_back(std::make_unique<Tile>(cx, cy, cz));
    ref[p]->setHaloSize(h.hx, h.hy, h.hz);
  }
  auto snapshot = [&] {
    for (int p = 0; p < P; p++)
      ref[p]->data() = out[p]->data();
  };
  auto compareStd = [&] {
    double maxDiff = 0, maxRef = 0;
    for (int p = 0; p < P; p++)
      for (int i = h.hx; i < cx - h.hx; i++)
        for (int j = h.hy; j < cy - h.hy; j++)
          for (int k = h.hz; k < cz - h.hz; k++) {
            maxDiff = std::max(maxDiff, std::abs(double((*out[p])(i, j, k)) - (*ref[p])(i, j, k)));
            maxRef = std::max(maxRef, std::abs(double((*ref[p])(i, j, k))));
          }
    return maxRef > 0 ? maxDiff / maxRef : maxDiff;
  };

  struct Check
  {
    char dir;
    std::function<void()> eigen, raw, dli;
  };
  const Check checks[] = { { 'x', stdX, rawX, dliX }, { 'y', stdY, rawY, dliY }, { 'z', stdZUnaligned, rawZ, dliZ } };
  const double tol = sizeof(Real) == 4 ? 1e-5 : 1e-12;
  bool ok = true;
  std::cout << std::scientific << std::setprecision(2);
  for (const auto& c : checks) {
    c.eigen();
    snapshot();
    c.raw();
    const double errRaw = compareStd();
    c.dli();
    const double errDli = compare(dliOut, ref, h);
    std::cout << "check D" << c.dir << ": max rel diff vs std/eigen  raw " << errRaw << "  dli " << errDli << "\n";
    ok = ok && errRaw < tol && errDli < tol;
  }
  std::cout << std::defaultfloat;
  if (!ok) {
    std::cerr << "kernel_bench: raw/DLI results disagree with the std/eigen kernels\n";
    return 1;
  }

  std::cout << "\nlayout dir variant        time [s]  Gcells/s\n";
  report("std", 'x', "eigen", timeIt(args.reps, stdX), cells, checksum(out));
  report("std", 'y', "eigen", timeIt(args.reps, stdY), cells, checksum(out));
  report("std", 'z', "legacy", timeIt(args.reps, stdZLegacy), cells, checksum(out));
  report("std", 'z', "unaligned", timeIt(args.reps, stdZUnaligned), cells, checksum(out));
  report("std", 'x', "raw", timeIt(args.reps, rawX), cells, checksum(out));
  report("std", 'y', "raw", timeIt(args.reps, rawY), cells, checksum(out));
  report("std", 'z', "raw", timeIt(args.reps, rawZ), cells, checksum(out));
  report("dli", 'x', "dli", timeIt(args.reps, dliX), cells, checksum(dliOut));
  report("dli", 'y', "dli", timeIt(args.reps, dliY), cells, checksum(dliOut));
  report("dli", 'z', "dli", timeIt(args.reps, dliZ), cells, checksum(dliOut));

  return 0;
}
