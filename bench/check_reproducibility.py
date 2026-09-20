#!/usr/bin/env python3
"""Check run-to-run reproducibility of a baseline CSV (see run_baseline.sh):
groups rows by configuration (tile size, process grid, precision, scheduler),
and for every configuration with more than one row left after discarding
--warmup of its earliest repeats, reports the relative spread of
total_time_s and time_per_step_s across the rest.

A history-based StarPU scheduler (e.g. dmda, run_baseline.sh's default) has
no performance-model data for a given (tile, grid, precision) combination
until it's actually seen that combination run once, so its first run tends
to schedule worse -- and less consistently -- than later ones. --warmup
discards that cold-start run per configuration before comparing the rest, so
it isn't mistaken for genuine run-to-run noise. Set it to 0 to compare every
repeat as-is.

Usage:
    check_reproducibility.py --csv bench/baseline.csv [--threshold 0.02] [--warmup 1]

Exits non-zero if any configuration's relative spread (after discarding
warm-up runs) exceeds --threshold on either column.
"""

import argparse
import csv
import statistics
import sys

GROUP_KEYS = ["cx", "cy", "cz", "P", "Q", "R", "precision", "scheduler"]
METRICS = ["total_time_s", "time_per_step_s"]
ORDER_KEY = "timestamp"


def relative_spread(values):
    """(max - min) / mean -- zero for a single value, well-defined for >=1."""
    if len(values) == 1:
        return 0.0
    mean = statistics.mean(values)
    if mean == 0:
        return float("inf") if max(values) != min(values) else 0.0
    return (max(values) - min(values)) / mean


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True, help="baseline CSV from run_baseline.sh")
    parser.add_argument("--threshold", type=float, default=0.02, help="max acceptable relative spread (default 0.02 = 2%%)")
    parser.add_argument("--warmup", type=int, default=1, help="earliest repeats to discard per configuration before comparing (default 1)")
    args = parser.parse_args()

    with open(args.csv, newline="") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise SystemExit(f"{args.csv}: no rows")

    groups = {}
    for row in rows:
        key = tuple(row[k] for k in GROUP_KEYS)
        groups.setdefault(key, []).append(row)

    ok = True
    skipped_configs = 0

    for key, group_rows in sorted(groups.items()):
        label = "cx={} cy={} cz={} P={} Q={} R={} precision={} scheduler={}".format(*key)

        group_rows = sorted(group_rows, key=lambda r: r[ORDER_KEY])
        warm_rows = group_rows[args.warmup :]

        if len(warm_rows) < 2:
            skipped_configs += 1
            print(
                f"{label}: {len(group_rows)} run(s), only {len(warm_rows)} left after discarding "
                f"{min(args.warmup, len(group_rows))} warm-up -- can't check reproducibility, repeat it more"
            )
            continue

        spreads = {}
        for metric in METRICS:
            values = [float(r[metric]) for r in warm_rows]
            spreads[metric] = relative_spread(values)

        worst = max(spreads.values())
        status = "PASS" if worst <= args.threshold else "FAIL"
        if worst > args.threshold:
            ok = False

        detail = "  ".join(f"{m}: spread={spreads[m]:.4%}" for m in METRICS)
        print(f"{label}: {len(warm_rows)}/{len(group_rows)} runs (after warm-up), {detail}  [{status}]")

    if skipped_configs:
        print(f"\n{skipped_configs} configuration(s) had too few runs left after discarding warm-up -- repeat them more.")

    if ok:
        print(f"\nPASS: all configurations within {args.threshold:.1%} relative spread (after discarding {args.warmup} warm-up run(s) each)")
    else:
        print(f"\nFAIL: at least one configuration exceeded {args.threshold:.1%} relative spread (after discarding {args.warmup} warm-up run(s) each)")

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
