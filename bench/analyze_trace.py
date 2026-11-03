#!/usr/bin/env python3
"""Decompose a StarPU run into pure-compute / overlap / pure-comm / idle time.

Reads the `tasks.rec` and `comms.rec` files `starpu_fxt_tool` produces from a
raw FxT trace (`prof_file_*`) and merges each into a set of non-overlapping
intervals: task execution windows (StartTime..EndTime, across all workers and
ranks) and MPI point-to-point transfers (SendTime..RecvTime). The union of
each set, and their pairwise intersection, gives four disjoint buckets over
the trace's observed time span:

    pure compute = time executing, no MPI transfer in flight
    overlap      = time executing AND transferring
    pure comm    = time transferring, nothing executing
    idle         = neither

These are wall-clock durations over the trace's own span (earliest event to
latest event), which only covers the time StarPU itself was actively
submitting/running tasks -- it's a subset of a full run's elapsed time (it
excludes setup/teardown before StarPU starts and after it stops).

Usage:
    analyze_trace.py --tasks tasks.rec --comms comms.rec
"""

import argparse
import re
import sys

RECORD_KV_RE = re.compile(r"^([A-Za-z]+):\s*(.*)$")


def parse_records(path):
    """Yield one dict per blank-line-separated 'Key: value' record block."""
    record = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                if record:
                    yield record
                    record = {}
                continue
            m = RECORD_KV_RE.match(line)
            if m:
                record[m.group(1)] = m.group(2)
    if record:
        yield record


# starpu_fxt_tool writes StartTime/EndTime/SendTime/RecvTime in milliseconds
# (src/debug/traces/starpu_fxt.c divides the raw nanosecond event time by 1e6).
MS_PER_S = 1000.0


def task_intervals(tasks_rec_path):
    # Internal StarPU bookkeeping records (e.g. data-acquire callbacks) share
    # this file but carry no StartTime/EndTime -- they never occupy a worker,
    # so skip them rather than treat a missing field as an error.
    intervals = []
    for rec in parse_records(tasks_rec_path):
        if "StartTime" not in rec or "EndTime" not in rec:
            continue
        start = float(rec["StartTime"]) / MS_PER_S
        end = float(rec["EndTime"]) / MS_PER_S
        if end > start:
            intervals.append((start, end))
    return intervals


def comm_intervals(comms_rec_path):
    intervals = []
    for rec in parse_records(comms_rec_path):
        if "SendTime" not in rec or "RecvTime" not in rec:
            continue
        start = float(rec["SendTime"]) / MS_PER_S
        end = float(rec["RecvTime"]) / MS_PER_S
        if end > start:
            intervals.append((start, end))
    return intervals


def merge_intervals(intervals):
    """Sort and merge overlapping/touching intervals into a disjoint union."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def union_duration(merged):
    return sum(end - start for start, end in merged)


def intersection_duration(merged_a, merged_b):
    """Sum of overlap between two disjoint, sorted interval lists."""
    i, j = 0, 0
    total = 0.0
    while i < len(merged_a) and j < len(merged_b):
        a_start, a_end = merged_a[i]
        b_start, b_end = merged_b[j]
        lo = max(a_start, b_start)
        hi = min(a_end, b_end)
        if lo < hi:
            total += hi - lo
        if a_end < b_end:
            i += 1
        else:
            j += 1
    return total


def analyze(tasks_rec_path, comms_rec_path):
    compute = merge_intervals(task_intervals(tasks_rec_path))
    comm = merge_intervals(comm_intervals(comms_rec_path))

    all_events = compute + comm
    if not all_events:
        raise SystemExit("No task or communication records found -- empty trace?")
    span_start = min(s for s, _ in all_events)
    span_end = max(e for _, e in all_events)
    span = span_end - span_start

    compute_union = union_duration(compute)
    comm_union = union_duration(comm)
    overlap = intersection_duration(compute, comm)

    pure_compute = compute_union - overlap
    pure_comm = comm_union - overlap
    idle = span - (compute_union + comm_union - overlap)

    return {
        "span_s": span,
        "pure_compute_s": pure_compute,
        "overlap_s": overlap,
        "pure_comm_s": pure_comm,
        "idle_s": idle,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tasks", required=True, help="tasks.rec path (from starpu_fxt_tool)")
    parser.add_argument("--comms", required=True, help="comms.rec path (from starpu_fxt_tool)")
    args = parser.parse_args()

    result = analyze(args.tasks, args.comms)

    for key, value in result.items():
        print(f"{key}={value:.6f}")


if __name__ == "__main__":
    main()
