#!/bin/bash
# make check: runs a small case (4 global tiles, TestA's 100 timesteps)
# through both the sequential oracle and the distributed engine under test,
# then compares their final velocity fields.
#
# Configurable via env vars (all have defaults matching a local dev build):
#   ORACLE_BIN      path to a SeWaS binary built with neither
#                    SEWAS_WITH_STARPU nor SEWAS_WITH_PARSEC (SEWASSequential)
#   CANDIDATE_BIN    path to the SeWaS binary under test (StarPU by default
#                    on this project)
#   CANDIDATE_LAUNCH launcher prefix for the candidate run, e.g.
#                    "mpirun -np 4" for a real distributed run; empty runs it
#                    directly (single rank)
#   CANDIDATE_P/Q/R  process grid for the candidate run (default 2 2 1,
#                    matching TestA's 2x2x1 global tile grid at cx=cy=cz=100)
#   DFILE            domain/material JSON (default data/input/TestA.json)
#   LINF_THRESHOLD   passed through to compare_to_oracle.py (default 1e-10)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

ORACLE_BIN="${ORACLE_BIN:-$REPO_ROOT/build/sewas}"
CANDIDATE_BIN="${CANDIDATE_BIN:-$REPO_ROOT/build-starpu/sewas}"
CANDIDATE_LAUNCH="${CANDIDATE_LAUNCH:-}"
CANDIDATE_P="${CANDIDATE_P:-2}"
CANDIDATE_Q="${CANDIDATE_Q:-2}"
CANDIDATE_R="${CANDIDATE_R:-1}"
DFILE="${DFILE:-$REPO_ROOT/data/input/TestA.json}"
LINF_THRESHOLD="${LINF_THRESHOLD:-1e-10}"

# TestA is nx=ny=200, nz=100; cx=cy=cz=100 gives a 2x2x1 = 4 global-tile case,
# matching the plan's "e.g. 4 tiles, 100 steps" example.
CX=100
CY=100
CZ=100

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

ORACLE_DIR="$WORKDIR/oracle"
CANDIDATE_DIR="$WORKDIR/candidate"
mkdir -p "$ORACLE_DIR" "$CANDIDATE_DIR"

echo "[run_check] oracle:    $ORACLE_BIN --P 1 --Q 1 --R 1"
"$ORACLE_BIN" \
  --cx "$CX" --cy "$CY" --cz "$CZ" \
  --P 1 --Q 1 --R 1 --nthreads 1 \
  --dfile="$DFILE" \
  --dump-velocity="$ORACLE_DIR" \
  > "$WORKDIR/oracle.log" 2>&1

echo "[run_check] candidate: $CANDIDATE_LAUNCH $CANDIDATE_BIN --P $CANDIDATE_P --Q $CANDIDATE_Q --R $CANDIDATE_R"
$CANDIDATE_LAUNCH "$CANDIDATE_BIN" \
  --cx "$CX" --cy "$CY" --cz "$CZ" \
  --P "$CANDIDATE_P" --Q "$CANDIDATE_Q" --R "$CANDIDATE_R" --nthreads 2 \
  --dfile="$DFILE" \
  --dump-velocity="$CANDIDATE_DIR" \
  > "$WORKDIR/candidate.log" 2>&1

python3 "$SCRIPT_DIR/compare_to_oracle.py" \
  --oracle-dir "$ORACLE_DIR" \
  --candidate-dir "$CANDIDATE_DIR" \
  --cx "$CX" --cy "$CY" --cz "$CZ" \
  --linf-threshold "$LINF_THRESHOLD"
