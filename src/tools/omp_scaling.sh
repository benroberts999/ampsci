#!/usr/bin/env bash

set -euo pipefail

# This program runs a given command in parallel using OpenMP for different 
# thread counts. Times the runs to test the scaling with OpenMP.
# Outputs a list of times and efficiencies.
# 
# Usage:
#   ./omp_scaling.sh [--threads M] [--inc D] [--repeat N] --exe executable [args...]
#
# threads     set maximum number of threads it will use {default 16}
# inc         set the increment {default 1}
# repeat      set how many times to average each rim {default 1}
# exe         (and any following flags) set the command that is executed
#
# exe is required; rest are optional
#
# Examples:
#   ./omp_scaling.sh --threads 6 --repeat 1 --inc 1 --exe ./tests "MBPT: Feynman unit tests"
#   ./omp_scaling.sh --exe ./tests [unit]
#   ./omp_scaling.sh --threads 24 --repeat 5 --inc 2 --exe ./tests [unit]
#   ./omp_scaling.sh --threads 24 --repeat 1 --inc 2 --exe ./ampsci ampsci.in

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 [--threads M] [--inc D] [--repeat N] --exe executable [args...]"
  exit 1
fi

max_threads=16
inc=2
repeat=1

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads)
      threads_opt="$2"
      shift 2
      ;;
    --inc)
      inc="$2"
      shift 2
      ;;
    --repeat)
      repeat="$2"
      shift 2
      ;;
    --exe)
      shift
      exe=("$@")
      break
      ;;
    *)
      echo "Unknown argument: $1"
      echo "Usage: $0 [--threads M] [--inc D] [--repeat N] --exe executable [args...]"
      exit 1
      ;;
  esac
done

echo "Executable: ${exe[*]}"
echo "Max threads: $max_threads"
echo "Increment: $inc"
echo "Repeats: $repeat"
echo

# Build thread list
thread_list=()


for ((t=1; t<=max_threads; t+=inc)); do
  thread_list+=("$t")
done


results=()
single_thread=0
for threads in "${thread_list[@]}"; do
  export OMP_NUM_THREADS=$threads

  echo "Running with OMP_NUM_THREADS=$threads ..."
  echo $OMP_NUM_THREADS

  total=0

  for ((i=1; i<=repeat; ++i)); do
    start=$(date +%s.%N)
    "${exe[@]}" >/dev/null 2>&1
    end=$(date +%s.%N)

    elapsed=$(echo "$end - $start" | bc)
    total=$(echo "$total + $elapsed" | bc)

    echo "     run $i: $elapsed s"
  done

  avg=$(echo "scale=6; $total / $repeat" | bc)

  # Save single-thread time
  if [[ $threads -eq 1 ]]; then
    single_thread="$avg"
  fi
  efficiency=$(echo "scale=3; $single_thread  / ($avg * $threads)" | bc)

  echo "Average    : $avg s"
  echo "Efficiencey: $efficiency "
  echo

  results+=("$threads $avg $efficiency")
done

echo "Summary (threads, time_in_seconds, efficient):"
for r in "${results[@]}"; do
  t=$(echo "$r" | awk '{print $1}')
  e=$(echo "$r" | awk '{print $2}')
  s=$(echo "$r" | awk '{print $3}')
  echo "$t  $e  $s"
done