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

/* Isolated-kernel benchmark for CentralFDOperator::applyZ<Variant> (Phase 1.1 of
 * sewas_transformation_plan_en.md): applies both the legacy scalar and unaligned Eigen Dz
 * kernels --reps times over a cached synthetic tile's (i,j) interior and reports elapsed
 * time / throughput for each, so the two can be A/B compared without touching the full
 * sewas binary. Complements bench/run_baseline.sh's integrated (full-run) benchmark. */

#include <chrono>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <string>

#include "CentralFDOperator.hxx"
#include "Config.hxx"

namespace {

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
                << "  --size N  interior edge length of the synthetic cubic tile (default 50)\n"
                << "  --reps N  number of times each variant sweeps the tile's (i,j) interior "
                   "(default 100)\n";
      exit(0);
    }
  }

  return args;
}

/* Times applyZ<Variant> over reps sweeps of fijk's (i,j) interior. Accumulates a checksum
 * from every returned tap and prints it, so the compiler can't prove the loop is dead and
 * elide it. */
template<SWS::KernelVariant Variant>
double
timeVariant(const CentralFDOperator& fdo, const SWS::SpatialBlockField<SWS::RealType>& fijk, const int reps)
{
  SWS::RealType checksum = 0.0;

  const auto iStart = fijk.iStart(), iEnd = fijk.iEnd();
  const auto jStart = fijk.jStart(), jEnd = fijk.jEnd();

  const auto start = std::chrono::steady_clock::now();

  for (int r = 0; r < reps; r++) {
    for (int i = iStart; i < iEnd; i++) {
      for (int j = jStart; j < jEnd; j++) {
        checksum += fdo.applyZ<Variant>(fijk, i, j, SWS::DK).sum();
      }
    }
  }

  const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

  std::cerr << "  (checksum: " << std::setprecision(17) << checksum << std::setprecision(6) << ")\n";

  return elapsed;
}

} // namespace

int
main(int argc, char** argv)
{
  const auto args = parseArgs(argc, argv);

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hny();
  const auto hnz = CentralFDOperator::hnz();

  const int cx = hnx + args.size + hnx;
  const int cy = hny + args.size + hny;
  const int cz = hnz + args.size + hnz;

  const SWS::RealType ds = 0.1;
  CentralFDOperator fdo(ds, ds, ds);

  SWS::SpatialBlockField<SWS::RealType> fijk(cx, cy, cz);
  fijk.setHaloSize(hnx, hny, hnz);

  for (int i = 0; i < cx; i++) {
    for (int j = 0; j < cy; j++) {
      for (int k = 0; k < cz; k++) {
        // Non-polynomial in k (and varying with i, j): the 4th-order stencil is exact on
        // polynomials up to degree 3, so a linear field would give a constant derivative
        // and a checksum blind to wrong tap offsets or coefficients.
        fijk(i, j, k) = static_cast<SWS::RealType>(std::sin(0.37 * k + 0.11 * i) * std::cos(0.23 * j) +
                                                   1e-3 * k * k * k * k);
      }
    }
  }

  const int interiorX = cx - 2 * hnx;
  const int interiorY = cy - 2 * hny;
  const long long taps = (long long)args.reps * interiorX * interiorY;

  std::cout << "kernel_bench: tile " << cx << "x" << cy << "x" << cz << " (interior " << interiorX << "x"
            << interiorY << "x" << (cz - 2 * hnz) << "), " << args.reps << " reps\n\n";

  std::cout << "legacy:\n";
  const auto legacyTime = timeVariant<SWS::KernelVariant::Legacy>(fdo, fijk, args.reps);
  std::cout << "  " << legacyTime << " s  (" << (taps / legacyTime) << " (i,j)-taps/s)\n\n";

  std::cout << "unaligned:\n";
  const auto unalignedTime = timeVariant<SWS::KernelVariant::Unaligned>(fdo, fijk, args.reps);
  std::cout << "  " << unalignedTime << " s  (" << (taps / unalignedTime) << " (i,j)-taps/s)\n";

  return 0;
}
