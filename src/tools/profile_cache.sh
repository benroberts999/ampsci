#!/usr/bin/env bash

set -euo pipefail

# Cache profiler: cache-miss rates + false-sharing check. Use to find out why
# parallel scaling stalls (memory/cache bound vs threads fighting over a line).
#
# Build once:   make MODE=profile
# Run:          ./src/tools/profile_cache.sh input.in
#
# Options:  --threads N  --blas N  --exe CMD  --out DIR   {default: ./profile/cache}
#           --no-c2c     skip the (slower) false-sharing run
#
# Writes:
#   CACHE.txt          IPC + L1/LLC miss rates (perf stat -dd)
#   FALSE_SHARING.txt  cache lines bouncing between cores (perf c2c)
#
# Reading it:
#   IPC < ~1.0                 -> CPU stalling on memory.
#   high LLC-load-miss %       -> spilling to RAM: bandwidth-bound (more threads
#                                 won't help; this is usually the scaling limit).
#   many remote HITMs on a line your threads share -> false sharing (fixable by
#                                 padding/separating the data).

DEFAULT_OUT=profile/cache
source "$(dirname "$0")/profile_common.sh"

do_c2c=1
args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-c2c) do_c2c=0; shift ;;
    *)        args+=("$1"); shift ;;
  esac
done

profile_parse "${args[@]}" || { profile_usage "$0"; exit 1; }
profile_check_perf || exit 1

# 1) Cache miss rates + IPC (perf stat -dd: one run, negligible overhead).
echo "Measuring cache counters (perf stat -dd) ..."
statout="$OUTDIR/.stat.raw"
perf stat -dd -o "$statout" -- "$EXE" "$INPUT" >"$OUTDIR/run.log" 2>&1 || true

set +e +o pipefail
rep="$OUTDIR/CACHE.txt"
{
  echo "# How to read:"
  echo "#   insn per cycle (IPC) < ~1.0  -> CPU stalling on memory."
  echo "#   L1-dcache-load-misses %      -> poor cache locality."
  echo "#   LLC-load-misses %            -> spilling to RAM (bandwidth-bound)."
  echo "#   '<not supported>/<not counted>' for cache events on Intel P/E hybrid"
  echo "#   CPUs is normal; pin to P-cores to measure them:"
  echo "#     OMP_PROC_BIND=close OMP_PLACES=cores ./src/tools/profile_cache.sh ..."
  echo
  echo "== Key metrics =="
  grep -iE "insn per cycle|dcache-load-misses|LLC-load-misses|cache-misses" \
    "$statout" || echo "(counters not available on this CPU)"
  echo
  echo "== Full perf stat -dd =="
  cat "$statout"
} > "$rep"

# 2) False-sharing detector (perf c2c: separate run, slower).
fs="$OUTDIR/FALSE_SHARING.txt"
if [[ $do_c2c -eq 1 ]]; then
  echo "Checking false sharing (perf c2c) ..."
  c2cdata="$OUTDIR/c2c.data"
  if perf c2c record -o "$c2cdata" -- "$EXE" "$INPUT" >>"$OUTDIR/run.log" 2>&1; then
    perf c2c report -i "$c2cdata" --stdio 2>/dev/null > "$OUTDIR/.c2c.raw"
    {
      echo "# How to read: false sharing = different threads writing the SAME"
      echo "# 64-byte line. In the table below, lines with many remote HITMs that"
      echo "# your hot data sits on are the culprits (fix: pad/separate the data)."
      echo
      sed -n '/Shared Data Cache Line Table/,/Cacheline/p' "$OUTDIR/.c2c.raw" | head -40
      echo
      echo "(full report: perf c2c report -i $c2cdata)"
    } > "$fs"
    rm -f "$OUTDIR/.c2c.raw"
  else
    echo "perf c2c not available on this CPU/kernel (skipped)." > "$fs"
  fi
else
  echo "skipped (--no-c2c)" > "$fs"
fi

echo
cat "$rep"
echo
echo "False sharing: $fs"
echo "All reports in: $OUTDIR/"
