#!/usr/bin/env bash

set -euo pipefail

# Allocation profiler: how many heap allocations, how big, and from where
# (heaptrack). Use when you suspect memory allocation is a bottleneck.
#
# Build once:   make MODE=profile
# Needs:        heaptrack   (sudo apt install heaptrack)
# Run:          ./src/tools/profile_alloc.sh input.in
#
# Options:  --threads N  --blas N  --exe CMD  --out DIR   {default: ./profile/alloc}
#
# Writes:
#   ALLOCATIONS.txt   summary: total allocations, peak heap, temporaries, leaks,
#                     and the call sites doing the most allocating.
#   heaptrack.zst     raw data (open with: heaptrack_gui)
#
# Reading it: a huge allocation COUNT with a small peak heap means lots of tiny
# short-lived allocations. That only matters if profile_hotspots.sh shows
# malloc/free hot -- otherwise it's harmless.

DEFAULT_OUT=profile/alloc
source "$(dirname "$0")/profile_common.sh"

profile_parse "$@" || { profile_usage "$0"; exit 1; }

if ! command -v heaptrack >/dev/null; then
  echo "heaptrack not found. Install with: sudo apt install heaptrack" >&2
  exit 1
fi

ht="$OUTDIR/heaptrack.zst"
echo "Recording allocations (heaptrack) ..."
heaptrack -o "${ht%.zst}" "$EXE" "$INPUT" >"$OUTDIR/run.log" 2>&1 || true

set +e +o pipefail
rep="$OUTDIR/ALLOCATIONS.txt"
raw="$OUTDIR/.heap.raw"
heaptrack_print "$ht" 2>/dev/null > "$raw"   # parse the .zst once, then slice
{
  echo "# How to read: allocation counts and the call sites doing the most."
  echo "# Many cheap allocations are fine -- only a worry if profile_hotspots.sh"
  echo "# shows malloc/free hot. (heaptrack's own summary is at the end.)"
  echo
  echo "== Summary =="
  grep -E "total runtime|calls to allocation functions:|temporary memory allocations:|peak heap memory|total memory leaked:" "$raw"
  echo
  echo "== Most-allocating call paths (count + source line) =="
  sed -n '/MOST CALLS TO ALLOCATION FUNCTIONS/,/MOST TEMPORARY/p' "$raw" \
    | grep -E "calls to allocation functions with|at src/" | head -24
} > "$rep"
rm -f "$raw"

echo
cat "$rep"
echo
echo "Full detail (interactive): heaptrack_gui $ht"
echo "All reports in: $OUTDIR/"
