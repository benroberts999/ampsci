#!/usr/bin/env bash

set -euo pipefail

# Hot-function profiler: where the CPU time goes (perf sampling).
#
# Build once:   make MODE=profile
# Run:          ./src/tools/profile_hotspots.sh input.in
#
# Options:  --threads N   OMP threads          {default: all cores}
#           --blas N      OPENBLAS/MKL threads  {default: 1}
#           --freq F      sampling Hz           {default: 999}
#           --top N       size of hot list      {default: 40}
#           --exe CMD     executable            {default: ./ampsci}
#           --out DIR     output dir            {default: ./profile/hotspots}
#
# Writes (each file starts with a one-line "how to read this"):
#   SUMMARY.txt          headline + top 15 functions
#   1_hot_functions.txt  flat list, time spent IN each function (self time)
#   2_callgraph.txt      top-down tree, time UNDER each call path
#   3_hot_callers.txt    who calls the hottest function
#
# Reading it: malloc/free/lock near the top => allocation is the cost (see
# profile_alloc.sh). Your own compute on top => real work.

DEFAULT_OUT=profile/hotspots
source "$(dirname "$0")/profile_common.sh"

freq=999
topn=40
# Pull out the tool-specific flags, pass the rest to the common parser.
args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --freq) freq="$2"; shift 2 ;;
    --top)  topn="$2"; shift 2 ;;
    *)      args+=("$1"); shift ;;
  esac
done

profile_parse "${args[@]}" || { profile_usage "$0"; exit 1; }
profile_check_perf || exit 1

data="$OUTDIR/perf.data"
echo "Recording CPU samples (perf, ${freq} Hz) ..."
# Frame-pointer call graph: MODE=profile keeps frame pointers, so this is far
# smaller/faster than dwarf and accurate enough for hotspots.
perf record -F "$freq" -g --call-graph fp -o "$data" -- \
  "$EXE" "$INPUT" >"$OUTDIR/run.log" 2>&1 || true

# Report slicing is best-effort: `grep | head` makes head close the pipe early
# (SIGPIPE), which would trip set -e/pipefail. Turn both off here.
set +e +o pipefail

# Dump the two perf views to raw files ONCE, then slice with grep/sed.
flat="$OUTDIR/.flat.raw"
tree="$OUTDIR/.tree.raw"
perf report -i "$data" --stdio --no-children -g none 2>/dev/null > "$flat"
perf report -i "$data" --stdio -g graph,0.5,caller 2>/dev/null > "$tree"

{
  echo "# How to read: % = time spent INSIDE this function itself (not callees)."
  echo "# This is your hot list. malloc/free here => allocation is the cost."
  echo
  grep -vE '^#|^[[:space:]]*$' "$flat" | head -n "$topn"
} > "$OUTDIR/1_hot_functions.txt"

{
  echo "# How to read: % = time UNDER this call (function + everything it calls)."
  echo "# Follow the heaviest branch down to see which path dominates."
  echo
  grep -vE '^#' "$tree" | head -120
} > "$OUTDIR/2_callgraph.txt"

hot="$(grep -vE '^#|^[[:space:]]*$' "$flat" | head -1 \
      | sed -E 's/^[ 0-9.%]+[a-z]+ +[^ ]+ +\[.\] //')"
perf report -i "$data" --stdio -g graph,0,callee 2>/dev/null > "$OUTDIR/.callee.raw"
{
  echo "# How to read: who calls the hottest function:"
  echo "#   $hot"
  echo "# Tells you which code path to attack to remove that cost."
  echo
  grep -F -A 40 "$hot" "$OUTDIR/.callee.raw" | head -45
} > "$OUTDIR/3_hot_callers.txt"

{
  echo "ampsci hot-function profile"
  echo "  input   : $INPUT"
  echo "  threads : OMP=$OMP_NUM_THREADS BLAS=$OPENBLAS_NUM_THREADS"
  echo
  echo "Top 15 functions by self-time (full list: 1_hot_functions.txt):"
  grep -vE '^#|^[[:space:]]*$' "$flat" | head -15 \
    | sed -E 's/\[\.\] //; s/\(.*//; s/^/  /'
  echo
  echo "Other reports: 1_hot_functions.txt, 2_callgraph.txt, 3_hot_callers.txt"
  echo "Explore interactively: perf report -i $data"
} > "$OUTDIR/SUMMARY.txt"

rm -f "$flat" "$tree" "$OUTDIR/.callee.raw"
echo
cat "$OUTDIR/SUMMARY.txt"
echo
echo "All reports in: $OUTDIR/"
