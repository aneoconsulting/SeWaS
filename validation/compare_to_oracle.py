#!/usr/bin/env python3
"""Compare a candidate SeWaS run's final velocity field against the sequential
oracle's.

Reads the raw per-(rank, component, local tile) binary dumps produced by
`sewas --dump-velocity=<dir>` (see src/main.cxx), reconstructs each run's
global velocity field, and reports the L2 and Linf norms of the difference
per component. Exits non-zero if any component's relative Linf exceeds
--linf-threshold.

Each dump file holds one tile's PADDED extent (cx+2*hn, cy+2*hn, cz+2*hn),
halo margin included, flattened in the SpatialBlockField's own order
(Y_MAJOR: k fastest) -- this script slices out the physical
[hn:hn+cx, hn:hn+cy, hn:hn+cz] interior before comparing.

The oracle is assumed to have been run with --P 1 --Q 1 --R 1 (its local
tile grid IS the global tile grid). The candidate's process grid (P, Q, R)
is derived from the ratio of the oracle's global tile counts to the
candidate's local tile counts, so it doesn't need to be passed explicitly.

Usage:
    compare_to_oracle.py --oracle-dir DIR --candidate-dir DIR \\
        --cx CX --cy CY --cz CZ [--hn 2] [--precision float64] \\
        [--linf-threshold 1e-10]
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

FILENAME_RE = re.compile(r"velocity_rank(\d+)_([xyz])_(\d+)_(\d+)_(\d+)\.bin$")

PRECISION_TO_DTYPE = {
    "float32": np.float32,
    "float64": np.float64,
}


def parse_dump_dir(dump_dir):
    """Return {(rank, component, ii, jj, kk): filepath} for a dump directory."""
    entries = {}
    for path in glob.glob(os.path.join(dump_dir, "velocity_rank*.bin")):
        m = FILENAME_RE.search(os.path.basename(path))
        if not m:
            continue
        rank, comp, ii, jj, kk = m.groups()
        entries[(int(rank), comp, int(ii), int(jj), int(kk))] = path
    if not entries:
        raise SystemExit(f"No velocity_rank*.bin files found in {dump_dir}")
    return entries


def load_tile_interior(path, cx, cy, cz, hn, dtype):
    ccx, ccy, ccz = cx + 2 * hn, cy + 2 * hn, cz + 2 * hn
    raw = np.fromfile(path, dtype=dtype)
    expected = ccx * ccy * ccz
    if raw.size != expected:
        raise SystemExit(
            f"{path}: expected {expected} elements ({ccx}x{ccy}x{ccz} padded tile), "
            f"got {raw.size}. Check --cx/--cy/--cz/--hn/--precision match the run."
        )
    # C order == Y_MAJOR (k fastest), matches Indexer<Y_MAJOR> / SpatialBlockField.
    padded = raw.reshape(ccx, ccy, ccz)
    return padded[hn : hn + cx, hn : hn + cy, hn : hn + cz]


def reconstruct_run(dump_dir, component, cx, cy, cz, hn, dtype):
    """Return (world, lnxx, lnyy, lnzz, {(rank, ii, jj, kk): interior array})."""
    entries = parse_dump_dir(dump_dir)
    comp_entries = {k: v for k, v in entries.items() if k[1] == component}
    if not comp_entries:
        raise SystemExit(f"No '{component}' velocity dumps found in {dump_dir}")

    ranks = sorted({r for (r, c, ii, jj, kk) in comp_entries})
    world = len(ranks)
    if ranks != list(range(world)):
        raise SystemExit(f"{dump_dir}: non-contiguous rank set {ranks}")

    lnxx = max(ii for (r, c, ii, jj, kk) in comp_entries) + 1
    lnyy = max(jj for (r, c, ii, jj, kk) in comp_entries) + 1
    lnzz = max(kk for (r, c, ii, jj, kk) in comp_entries) + 1

    expected_tiles = world * lnxx * lnyy * lnzz
    if len(comp_entries) != expected_tiles:
        raise SystemExit(
            f"{dump_dir}: found {len(comp_entries)} '{component}' tile dumps, "
            f"expected {expected_tiles} ({world} ranks x {lnxx}x{lnyy}x{lnzz} local "
            "tiles) -- a rank's dump is incomplete or the tiling isn't uniform."
        )

    tiles = {
        (r, ii, jj, kk): load_tile_interior(path, cx, cy, cz, hn, dtype)
        for (r, c, ii, jj, kk), path in comp_entries.items()
    }
    return world, lnxx, lnyy, lnzz, tiles


def assemble_global(world, lnxx, lnyy, lnzz, tiles, P, Q, R, cx, cy, cz, dtype):
    nxx, nyy, nzz = lnxx * P, lnyy * Q, lnzz * R
    nx, ny, nz = nxx * cx, nyy * cy, nzz * cz
    out = np.empty((nx, ny, nz), dtype=dtype)

    for rank in range(world):
        # Inverse of rank_of(ii,jj,kk) = P*Q*floor(kk/lnzz) + P*floor(jj/lnyy) + floor(ii/lnxx),
        # see ../../graviton_tests/parallelism.md.
        ii_block = rank % P
        jj_block = (rank // P) % Q
        kk_block = rank // (P * Q)

        for lii in range(lnxx):
            for ljj in range(lnyy):
                for lkk in range(lnzz):
                    gii = ii_block * lnxx + lii
                    gjj = jj_block * lnyy + ljj
                    gkk = kk_block * lnzz + lkk
                    block = tiles[(rank, lii, ljj, lkk)]
                    out[
                        gii * cx : (gii + 1) * cx,
                        gjj * cy : (gjj + 1) * cy,
                        gkk * cz : (gkk + 1) * cz,
                    ] = block
    return out


def reconstruct_and_assemble(dump_dir, component, cx, cy, cz, hn, dtype, P=None, Q=None, R=None):
    world, lnxx, lnyy, lnzz, tiles = reconstruct_run(dump_dir, component, cx, cy, cz, hn, dtype)

    if P is None:
        # Oracle case: single rank, its local grid IS the global grid.
        if world != 1:
            raise SystemExit(f"{dump_dir}: expected a single-rank oracle dump, found {world} ranks")
        P, Q, R = 1, 1, 1
    elif P * Q * R != world:
        raise SystemExit(
            f"{dump_dir}: derived process grid {P}x{Q}x{R} = {P * Q * R} ranks, "
            f"but the dump has {world} ranks"
        )

    global_field = assemble_global(world, lnxx, lnyy, lnzz, tiles, P, Q, R, cx, cy, cz, dtype)
    return global_field, (lnxx, lnyy, lnzz), (P, Q, R)


def compare(oracle_dir, candidate_dir, cx, cy, cz, hn, dtype, linf_threshold):
    results = {}
    ok = True

    for component in ("x", "y", "z"):
        oracle_field, (o_lnxx, o_lnyy, o_lnzz), _ = reconstruct_and_assemble(
            oracle_dir, component, cx, cy, cz, hn, dtype
        )
        # Oracle ran with P=Q=R=1, so its local tile grid is the global one.
        oracle_nxx, oracle_nyy, oracle_nzz = o_lnxx, o_lnyy, o_lnzz

        cand_world, cand_lnxx, cand_lnyy, cand_lnzz, cand_tiles = reconstruct_run(
            candidate_dir, component, cx, cy, cz, hn, dtype
        )
        if oracle_nxx % cand_lnxx or oracle_nyy % cand_lnyy or oracle_nzz % cand_lnzz:
            raise SystemExit(
                f"Candidate's local tile grid ({cand_lnxx}x{cand_lnyy}x{cand_lnzz}) doesn't "
                f"evenly divide the oracle's global tile grid ({oracle_nxx}x{oracle_nyy}x{oracle_nzz}) "
                "-- the two runs weren't given the same domain/tiling."
            )
        P = oracle_nxx // cand_lnxx
        Q = oracle_nyy // cand_lnyy
        R = oracle_nzz // cand_lnzz

        candidate_field = assemble_global(
            cand_world, cand_lnxx, cand_lnyy, cand_lnzz, cand_tiles, P, Q, R, cx, cy, cz, dtype
        )

        if candidate_field.shape != oracle_field.shape:
            raise SystemExit(
                f"Component {component}: oracle shape {oracle_field.shape} != "
                f"candidate shape {candidate_field.shape}"
            )

        diff = candidate_field.astype(np.float64) - oracle_field.astype(np.float64)
        oracle_f64 = oracle_field.astype(np.float64)

        l2 = np.sqrt(np.sum(diff * diff))
        linf = np.max(np.abs(diff))
        oracle_l2 = np.sqrt(np.sum(oracle_f64 * oracle_f64))
        oracle_linf = np.max(np.abs(oracle_f64))

        rel_l2 = l2 / oracle_l2 if oracle_l2 > 0 else l2
        rel_linf = linf / oracle_linf if oracle_linf > 0 else linf

        results[component] = {
            "l2": l2,
            "linf": linf,
            "rel_l2": rel_l2,
            "rel_linf": rel_linf,
            "candidate_grid": (P, Q, R),
        }

        if rel_linf > linf_threshold:
            ok = False

    return ok, results


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--oracle-dir", required=True, help="--dump-velocity dir from the P=Q=R=1 oracle run")
    parser.add_argument("--candidate-dir", required=True, help="--dump-velocity dir from the run being validated")
    parser.add_argument("--cx", type=int, required=True)
    parser.add_argument("--cy", type=int, required=True)
    parser.add_argument("--cz", type=int, required=True)
    parser.add_argument("--hn", type=int, default=2, help="halo width; must match CentralFDOperator::hnx() (default 2)")
    parser.add_argument("--precision", choices=PRECISION_TO_DTYPE.keys(), default="float64")
    parser.add_argument("--linf-threshold", type=float, default=1e-10, help="max acceptable relative Linf per component")
    args = parser.parse_args()

    dtype = PRECISION_TO_DTYPE[args.precision]

    ok, results = compare(
        args.oracle_dir, args.candidate_dir, args.cx, args.cy, args.cz, args.hn, dtype, args.linf_threshold
    )

    for component, r in results.items():
        P, Q, R = r["candidate_grid"]
        print(
            f"v{component}: L2={r['l2']:.6e} (rel {r['rel_l2']:.3e})  "
            f"Linf={r['linf']:.6e} (rel {r['rel_linf']:.3e})  candidate grid={P}x{Q}x{R}"
        )

    if ok:
        print(f"PASS: all components within relative Linf < {args.linf_threshold:.1e}")
    else:
        print(f"FAIL: at least one component exceeded relative Linf {args.linf_threshold:.1e}")

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
