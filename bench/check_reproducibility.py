#!/usr/bin/env python3
"""Check run-to-run reproducibility of a baseline CSV (see run_baseline.sh):
groups rows by configuration (tile size, process grid, precision, scheduler),
and for every configuration with more than one row, reports the relative
spread of total_time_s and time_per_step_s across those repeats.

Usage:
    check_reproducibility.py --csv bench/baseline.csv [--threshold 0.02]

Exits non-zero if any configuration's relative spread exceeds --threshold on
either column.
"""

import argparse
import csv
import statistics
import sys

GROUP_KEYS = ["cx", "cy", "cz", "P", "Q", "R", "precision", "scheduler"]
METRICS = ["total_time_s", "time_per_step_s"]


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
    single_run_configs = 0

    for key, group_rows in sorted(groups.items()):
        label = "cx={} cy={} cz={} P={} Q={} R={} precision={} scheduler={}".format(*key)
        if len(group_rows) == 1:
            single_run_configs += 1
            print(f"{label}: only 1 run, can't check reproducibility")
            continue

        spreads = {}
        for metric in METRICS:
            values = [float(r[metric]) for r in group_rows]
            spreads[metric] = relative_spread(values)

        worst = max(spreads.values())
        status = "PASS" if worst <= args.threshold else "FAIL"
        if worst > args.threshold:
            ok = False

        detail = "  ".join(f"{m}: spread={spreads[m]:.4%}" for m in METRICS)
        print(f"{label}: {len(group_rows)} runs, {detail}  [{status}]")

    if single_run_configs:
        print(f"\n{single_run_configs} configuration(s) had only 1 run -- repeat them to actually check reproducibility.")

    if ok:
        print(f"\nPASS: all repeated configurations within {args.threshold:.1%} relative spread")
    else:
        print(f"\nFAIL: at least one configuration exceeded {args.threshold:.1%} relative spread")

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
