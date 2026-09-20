#!/usr/bin/env python3
"""Roofline figures for each run in a baseline CSV (see run_baseline.sh):
arithmetic intensity of the stencil kernel, and the incompressible
communication time implied by the halo volume at a given network bandwidth.

Arithmetic intensity (FLOP/byte) is derived from the exact operation count in
the velocity/stress update kernels (src/LinearSeismicWaveModel.cxx) and the
distinct fields each phase touches (include/DataSet.hxx, include/Config.hxx),
not a rough estimate -- see FLOPS_PER_CELL_PER_TIMESTEP and
BYTES_PER_CELL_PER_TIMESTEP below for the derivation.

Halo volume is the geometric inter-rank exchange volume implied by the
process grid: for each axis, every one of the (P-1)/(Q-1)/(R-1) inter-rank
boundaries along it exchanges a halo of depth HN in both directions, across
the full cross-section of the other two axes, for every field that gets
halo-exchanged.

Usage:
    roofline.py --csv bench/baseline.csv --bandwidth-gbps 12.5 [--out bench/roofline.csv]

--bandwidth-gbps is the injectable network bandwidth (GB/s, decimal: 1 GB/s =
1e9 bytes/s) to compute the communication floor against. Measure it for the
target network (e.g. an MPI ping-pong/bisection-bandwidth benchmark) --
this script doesn't measure it itself.
"""

import argparse
import csv
import json
import math
import sys

# CentralFDOperator's 4th-order stencil taps (include/CentralFDOperator.hxx):
# d1=1, d2=2, d3=0, d4=1 -> halo width hn = max(d1..d4) = 2, on every axis.
# Each directional derivative evaluation is 4 multiplies + 3 adds = 7 FLOPs
# (confirmed by the "7 flops" comment at its scalar fallback implementation).
HN = 2
FLOPS_PER_STENCIL_TAP = 7

# Exact operation count per cell for one full (velocity + stress) timestep,
# counted directly from src/LinearSeismicWaveModel.cxx's update expressions:
#
# Velocity (vx, vy, vz): each is `field += (tap+tap+tap) * dt * buoyancy`
#   -- 3 stencil taps (3*7=21) + 2 adds (summing 3 taps) + 1 mul (*dt)
#   + 1 mul (*buoyancy) + 1 add (+=) = 26 FLOPs/component, x3 components.
_VELOCITY_FLOPS = 3 * (3 * FLOPS_PER_STENCIL_TAP + 2 + 1 + 1 + 1)

# Stress, diagonal (xx, yy, zz): `field += ((tap+tap)*lambda + (lambda+2*mu)*tap) * dt`
#   -- 3 taps (21) + 1 add + 1 mul(*lambda) + 1 mul(2*mu) + 1 add(lambda+2mu)
#   + 1 mul(*tap) + 1 add(combine) + 1 mul(*dt) + 1 add(+=) = 29 FLOPs, x3.
# Stress, off-diagonal (xy, xz, yz): `field += (tap+tap) * mu * dt`
#   -- 2 taps (14) + 1 add + 1 mul(*mu) + 1 mul(*dt) + 1 add(+=) = 18 FLOPs, x3.
_STRESS_FLOPS = 3 * (3 * FLOPS_PER_STENCIL_TAP + 8) + 3 * (2 * FLOPS_PER_STENCIL_TAP + 4)

FLOPS_PER_CELL_PER_TIMESTEP = _VELOCITY_FLOPS + _STRESS_FLOPS  # 219

# Distinct fields each phase's working set touches, per cell, in itemsize units
# (1 = one field read or written once; read+write of the same field counts 2):
#
# Velocity phase: writes vx,vy,vz (3*2=6); reads all 6 stress components (the
# union of what vx/vy/vz's update expressions reference) and all 3 buoyancy
# components (include/DataSet.hxx: Buoyancy is DIM=3-sized, one used per
# velocity component) -- 6+3 read-only.
_VELOCITY_ITEMS = 3 * 2 + 6 + 3

# Stress phase: writes all 6 stress components (6*2=12); reads vx,vy,vz (3,
# the union across all 6 update expressions); reads elasticity coefficients
# -- lambda_ and mu_ are each their own NB_STRESS_FIELD_COMPONENTS=6-sized
# array (include/DataSet.hxx), but only lambda(sc)+mu(sc) for sc in
# {xx,yy,zz} (3+3) and mu(sc) for sc in {xy,xz,yz} (3) are ever read (the 3
# off-diagonal lambda_ entries are allocated but never used) -- 3+6=9 read-only.
_STRESS_ITEMS = 6 * 2 + 3 + 9

ITEMS_PER_CELL_PER_TIMESTEP = _VELOCITY_ITEMS + _STRESS_ITEMS  # 39

PRECISION_ITEMSIZE = {
    "float64": 8,
    "float32": 4,
}

# Fields halo-exchanged per phase (include/HaloManager.hxx): velocity phase
# exchanges the 3 velocity components' boundary cells, stress phase the 6
# stress components'.
HALO_FIELDS_PER_TIMESTEP = 3 + 6


def derive_domain(dfile, cx, cy, cz, P, Q, R):
    """Reproduce SEWASParameterManager::parseDataFile()'s padding so nx/ny/nz
    match what the run actually used."""
    with open(dfile) as f:
        d = json.load(f)
    nx = math.ceil(d["lx"] / d["ds"])
    ny = math.ceil(d["ly"] / d["ds"])
    nz = math.ceil(d["lz"] / d["ds"])

    nxx = math.ceil(nx / cx)
    nyy = math.ceil(ny / cy)
    nzz = math.ceil(nz / cz)

    lnxx = math.ceil(nxx / P)
    lnyy = math.ceil(nyy / Q)
    lnzz = math.ceil(nzz / R)

    nxx, nyy, nzz = lnxx * P, lnyy * Q, lnzz * R
    nx, ny, nz = nxx * cx, nyy * cy, nzz * cz

    return nx, ny, nz


def halo_bytes_per_step(nx, ny, nz, P, Q, R, itemsize):
    """Total inter-rank halo exchange volume for one full timestep, summed
    over every rank-to-rank boundary in the process grid."""
    bytes_x = 2 * max(P - 1, 0) * ny * nz * HN * itemsize
    bytes_y = 2 * max(Q - 1, 0) * nx * nz * HN * itemsize
    bytes_z = 2 * max(R - 1, 0) * nx * ny * HN * itemsize
    return (bytes_x + bytes_y + bytes_z) * HALO_FIELDS_PER_TIMESTEP


def process_row(row, bandwidth_bytes_per_s):
    cx, cy, cz = int(row["cx"]), int(row["cy"]), int(row["cz"])
    P, Q, R = int(row["P"]), int(row["Q"]), int(row["R"])
    precision = row["precision"]

    if precision not in PRECISION_ITEMSIZE:
        raise SystemExit(
            f"Unrecognized precision label {precision!r} (expected one of "
            f"{list(PRECISION_ITEMSIZE)}); pass a matching --precision-label to run_baseline.sh"
        )
    itemsize = PRECISION_ITEMSIZE[precision]

    flop_per_byte = FLOPS_PER_CELL_PER_TIMESTEP / (ITEMS_PER_CELL_PER_TIMESTEP * itemsize)

    nx, ny, nz = derive_domain(row["dfile"], cx, cy, cz, P, Q, R)
    halo_bytes = halo_bytes_per_step(nx, ny, nz, P, Q, R, itemsize)
    comm_floor_s = halo_bytes / bandwidth_bytes_per_s

    time_per_step_s = float(row["time_per_step_s"])
    steering_ratio = comm_floor_s / time_per_step_s if time_per_step_s > 0 else float("nan")

    return {
        **row,
        "arithmetic_intensity_flop_per_byte": f"{flop_per_byte:.6f}",
        "halo_bytes_per_step": str(halo_bytes),
        "comm_floor_s": f"{comm_floor_s:.6f}",
        "comm_floor_over_compute": f"{steering_ratio:.6f}",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True, help="baseline CSV from run_baseline.sh")
    parser.add_argument("--bandwidth-gbps", type=float, required=True, help="injectable network bandwidth, GB/s (decimal: 1 GB/s = 1e9 B/s)")
    parser.add_argument("--out", help="output CSV (default: print to stdout)")
    args = parser.parse_args()

    bandwidth_bytes_per_s = args.bandwidth_gbps * 1e9

    with open(args.csv, newline="") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise SystemExit(f"{args.csv}: no rows")

    out_rows = [process_row(row, bandwidth_bytes_per_s) for row in rows]

    out_file = open(args.out, "w", newline="") if args.out else sys.stdout
    writer = csv.DictWriter(out_file, fieldnames=list(out_rows[0].keys()))
    writer.writeheader()
    writer.writerows(out_rows)
    if args.out:
        out_file.close()
        print(f"wrote {len(out_rows)} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
