#!/usr/bin/env bash
#
# Shared helpers for the profile_*.sh tools. SOURCED, not run directly.
# Provides common argument parsing, run setup, and perf availability checks so
# each tool script stays short and only contains its own tool + formatting.

# Print the leading comment block of a script as its help text.
profile_usage() { grep '^#' "$1" | sed 's/^# \{0,1\}//'; }

# Parse args common to all tools. Sets globals: THREADS BLAS EXE INPUT OUTDIR.
# A script sets DEFAULT_OUT before calling. Returns 1 on bad/empty args.
#   --threads N  --blas N  --exe CMD  --out DIR  <input-file>
profile_parse() {
  THREADS="$(nproc)"
  BLAS=1
  EXE="./ampsci"
  OUTDIR=""
  INPUT=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --threads) THREADS="$2"; shift 2 ;;
      --blas)    BLAS="$2";    shift 2 ;;
      --exe)     EXE="$2";     shift 2 ;;
      --out)     OUTDIR="$2";  shift 2 ;;
      -h|--help) return 1 ;;
      *)         INPUT="$1";   shift ;;
    esac
  done

  if [[ -z "$INPUT" ]]; then
    echo "Error: no input file given." >&2
    return 1
  fi
  if [[ ! -x "$EXE" ]]; then
    echo "Error: '$EXE' not found. Build a profiling binary first:" >&2
    echo "    make MODE=profile" >&2
    return 1
  fi
  [[ -n "$OUTDIR" ]] || OUTDIR="./${DEFAULT_OUT:-profile_out}"

  export OMP_NUM_THREADS="$THREADS"
  export OPENBLAS_NUM_THREADS="$BLAS"
  export MKL_NUM_THREADS="$BLAS"

  mkdir -p "$OUTDIR"

  echo "  exe     : $EXE $INPUT"
  echo "  threads : OMP=$OMP_NUM_THREADS  BLAS=$OPENBLAS_NUM_THREADS"
  echo "  output  : $OUTDIR/"
  echo
  return 0
}

# Check perf can sample (installed + permissive paranoid level, or root).
# Prints the one-line fix and returns 1 if not.
profile_check_perf() {
  if ! command -v perf >/dev/null; then
    echo "perf not found. Install linux-tools for your kernel." >&2
    return 1
  fi
  local p
  p="$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo 4)"
  if [[ "$p" -gt 2 && "${EUID:-$(id -u)}" -ne 0 ]]; then
    echo "perf cannot sample: kernel.perf_event_paranoid = $p (need <= 2)." >&2
    echo "Run this one line once (resets at reboot), then re-run:" >&2
    echo "    sudo sysctl kernel.perf_event_paranoid=1" >&2
    return 1
  fi
  return 0
}
