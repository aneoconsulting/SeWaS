#!/bin/bash
# Sweeps tile size x process grid x precision, running each combination and
# appending one CSV row per run: total/per-step wall time, and -- when the
# binary was built with FxT tracing (SEWAS_STARPU_ENABLE_FXT/FXT_ROOT, see
# cmake/resources/starpu/CMakeLists.txt) -- a decomposition of the run into
# pure-compute / overlap / pure-comm / idle time (bench/analyze_trace.py).
#
# Usage:
#   bench/run_baseline.sh --binary PATH --precision-label LABEL \
#       --tiles "CXxCYxCZ[,CXxCYxCZ...]" --grids "PxQxR[,PxQxR...]" \
#       --nthreads N [options...]
#
# To cover the full tile-size x node-count x precision matrix, invoke this
# script once per binary/precision (double vs float are separate compile-time
# builds -- see SEWAS_SINGLE_PRECISION in CMakeLists.txt); each invocation can
# already sweep --tiles and --grids in one pass, appending to the same --csv.
#
# Options:
#   --binary PATH           sewas binary to run (required)
#   --precision-label LABEL free-form label recorded in the CSV, e.g.
#                            "float64"/"float32" (required)
#   --tiles LIST             comma-separated CXxCYxCZ combinations (required)
#   --grids LIST             comma-separated PxQxR combinations (required)
#   --nthreads N             StarPU worker threads per rank (required)
#   --dfile PATH             domain/material JSON (default data/input/TestA.json)
#   --scheduler NAME          STARPU_SCHED value, fixed and recorded for every
#                            run in the sweep since it changes the balance
#                            being measured (default dmda)
#   --launcher-template CMD  launcher prefix with a "{np}" placeholder for the
#                            rank count, e.g. "mpirun -np {np} --host h1,h2,h3,h4";
#                            "{np}" is substituted with P*Q*R for each grid.
#                            Default: "mpirun -np {np} --oversubscribe" (single
#                            node, oversubscribing cores -- fine for P*Q*R <=
#                            the node's core count, not for a real multi-node run)
#                            For a real multi-node run, the template also needs
#                            to forward STARPU_FXT_TRACE, STARPU_GENERATE_TRACE
#                            and STARPU_FXT_PREFIX to every rank (e.g. OpenMPI's
#                            "-x STARPU_FXT_TRACE -x STARPU_GENERATE_TRACE -x
#                            STARPU_FXT_PREFIX"), the same way it already needs
#                            to forward LD_LIBRARY_PATH -- mpirun does not
#                            forward the launching shell's environment to
#                            remote ranks on its own. Without this, only the
#                            launching rank's trace gets recorded, so
#                            starpu_fxt_tool sees just one of several ranks and
#                            can't resolve any MPI transfers between them
#                            (comm/overlap silently come out as zero, not an
#                            error).
#   --csv PATH                output file (default bench/baseline.csv);
#                            appended to, header written only if the file is new
#   --trace-dir DIR           where to keep each run's raw FxT trace + the
#                            starpu_fxt_tool outputs derived from it, named
#                            <cx>x<cy>x<cz>_<P>x<Q>x<R>_<precision>/ (default
#                            bench/traces). For a multi-node run this MUST be
#                            on storage shared across every node the launcher
#                            can place a rank on (e.g. NFS-mounted $HOME) --
#                            every rank writes its raw trace under here
#                            (STARPU_FXT_PREFIX), and this script only looks
#                            for those files on the node it itself runs on.
#                            Skipped (with a warning) if the binary wasn't
#                            built with FxT support (no prof_file produced).
#   --starpu-fxt-tool PATH    starpu_fxt_tool binary (default: same directory
#                            as --binary's StarPU install, if found on PATH)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

BINARY=""
PRECISION_LABEL=""
TILES=""
GRIDS=""
NTHREADS=""
DFILE="$REPO_ROOT/data/input/TestA.json"
SCHEDULER="dmda"
LAUNCHER_TEMPLATE="mpirun -np {np} --oversubscribe"
CSV="$SCRIPT_DIR/baseline.csv"
TRACE_DIR="$SCRIPT_DIR/traces"
STARPU_FXT_TOOL="$(command -v starpu_fxt_tool || true)"

while [ $# -gt 0 ]; do
  case "$1" in
    --binary) BINARY="$2"; shift 2 ;;
    --precision-label) PRECISION_LABEL="$2"; shift 2 ;;
    --tiles) TILES="$2"; shift 2 ;;
    --grids) GRIDS="$2"; shift 2 ;;
    --nthreads) NTHREADS="$2"; shift 2 ;;
    --dfile) DFILE="$2"; shift 2 ;;
    --scheduler) SCHEDULER="$2"; shift 2 ;;
    --launcher-template) LAUNCHER_TEMPLATE="$2"; shift 2 ;;
    --csv) CSV="$2"; shift 2 ;;
    --trace-dir) TRACE_DIR="$2"; shift 2 ;;
    --starpu-fxt-tool) STARPU_FXT_TOOL="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

for required in BINARY PRECISION_LABEL TILES GRIDS NTHREADS; do
  if [ -z "${!required}" ]; then
    echo "Missing required option: --$(echo "$required" | tr '[:upper:]' '[:lower:]' | tr '_' '-')" >&2
    exit 1
  fi
done

if [ -n "$STARPU_FXT_TOOL" ]; then
  echo "[run_baseline] trace decomposition: enabled ($STARPU_FXT_TOOL)"
else
  echo "[run_baseline] trace decomposition: disabled (starpu_fxt_tool not found; pass --starpu-fxt-tool or put it on PATH)"
fi

if [ ! -f "$CSV" ]; then
  mkdir -p "$(dirname "$CSV")"
  echo "timestamp,cx,cy,cz,P,Q,R,nthreads,precision,scheduler,dfile,total_time_s,time_per_step_s,trace_span_s,pure_compute_s,overlap_s,pure_comm_s,idle_s" > "$CSV"
fi

parse_stat() {
  # $1: MetricsManager section name, $2: sewas stdout log
  # Prints that section's "Elapsed Time" AVG column (the first numeric field
  # after "Elapsed Time (s) :").
  awk -v section="$1" '
    $0 ~ "^" section "$" { in_section = 1; next }
    in_section && /^[A-Za-z]/ && $0 !~ "^" section "$" { in_section = 0 }
    in_section && /Elapsed Time/ { print $5; exit }
  ' "$2"
}

run_one() {
  local cx="$1" cy="$2" cz="$3" P="$4" Q="$5" R="$6"
  local np=$((P * Q * R))
  local run_id="${cx}x${cy}x${cz}_${P}x${Q}x${R}_${PRECISION_LABEL}"
  local run_dir
  run_dir="$(mktemp -d)"
  trap 'rm -rf "$run_dir"' RETURN

  local launcher="${LAUNCHER_TEMPLATE//\{np\}/$np}"

  echo "[run_baseline] $run_id: $launcher $BINARY --nthreads $NTHREADS"

  # Every rank writes its raw trace under STARPU_FXT_PREFIX, and ranks can run
  # on different physical nodes -- this has to be reachable from all of them,
  # not just wherever this script itself runs, so it lives under --trace-dir
  # (the caller's responsibility to put on shared storage for a real
  # multi-node run) rather than a node-local mktemp.
  local trace_out_dir="$TRACE_DIR/$run_id"
  local fxt_prefix="$trace_out_dir/raw/"
  mkdir -p "$fxt_prefix"
  local log="$run_dir/run.log"

  STARPU_SCHED="$SCHEDULER" \
  STARPU_FXT_TRACE=1 STARPU_GENERATE_TRACE=1 STARPU_FXT_PREFIX="$fxt_prefix" \
  $launcher "$BINARY" \
    --cx "$cx" --cy "$cy" --cz "$cz" \
    --P "$P" --Q "$Q" --R "$R" --nthreads "$NTHREADS" \
    --dfile="$DFILE" \
    > "$log" 2>&1

  local total_time
  local core_sim_time
  total_time="$(parse_stat Global "$log")"
  core_sim_time="$(parse_stat "Core simulation" "$log")"

  if [ -z "$total_time" ] || [ -z "$core_sim_time" ]; then
    echo "[run_baseline] $run_id: FAILED to parse timing from $log (see it for the run's own output)" >&2
    cat "$log" >&2
    return 1
  fi

  local nt
  nt="$(python3 -c "
import json
with open('$DFILE') as f:
    d = json.load(f)
print(int(-(-d['tmax'] // d['dt'])))  # ceil, matches SEWASParameterManager::parseDataFile
")"
  local time_per_step
  time_per_step="$(python3 -c "print($core_sim_time / $nt)")"

  local trace_span="" pure_compute="" overlap="" pure_comm="" idle=""
  if [ -n "$STARPU_FXT_TOOL" ]; then
    local prof_files=("$fxt_prefix"prof_file_*)
    if [ -e "${prof_files[0]}" ]; then
      local fxt_args=()
      for f in "${prof_files[@]}"; do
        fxt_args+=(-i "$f")
      done
      (cd "$trace_out_dir" && "$STARPU_FXT_TOOL" "${fxt_args[@]}" -o paje.trace > fxt_tool.log 2>&1) || {
        echo "[run_baseline] $run_id: starpu_fxt_tool failed, see $trace_out_dir/fxt_tool.log" >&2
      }
      if [ -f "$trace_out_dir/tasks.rec" ] && [ -f "$trace_out_dir/comms.rec" ]; then
        local analysis
        analysis="$(python3 "$SCRIPT_DIR/analyze_trace.py" --tasks "$trace_out_dir/tasks.rec" --comms "$trace_out_dir/comms.rec")"
        trace_span="$(echo "$analysis" | grep '^span_s=' | cut -d= -f2)"
        pure_compute="$(echo "$analysis" | grep '^pure_compute_s=' | cut -d= -f2)"
        overlap="$(echo "$analysis" | grep '^overlap_s=' | cut -d= -f2)"
        pure_comm="$(echo "$analysis" | grep '^pure_comm_s=' | cut -d= -f2)"
        idle="$(echo "$analysis" | grep '^idle_s=' | cut -d= -f2)"
      fi
    else
      echo "[run_baseline] $run_id: no prof_file produced -- binary likely wasn't built with FxT support" >&2
    fi
  fi

  echo "$(date -Iseconds),$cx,$cy,$cz,$P,$Q,$R,$NTHREADS,$PRECISION_LABEL,$SCHEDULER,$DFILE,$total_time,$time_per_step,$trace_span,$pure_compute,$overlap,$pure_comm,$idle" >> "$CSV"
}

IFS=',' read -ra TILE_LIST <<< "$TILES"
IFS=',' read -ra GRID_LIST <<< "$GRIDS"

for tile in "${TILE_LIST[@]}"; do
  IFS='x' read -r cx cy cz <<< "$tile"
  for grid in "${GRID_LIST[@]}"; do
    IFS='x' read -r P Q R <<< "$grid"
    run_one "$cx" "$cy" "$cz" "$P" "$Q" "$R"
  done
done

echo "[run_baseline] done, results in $CSV"
