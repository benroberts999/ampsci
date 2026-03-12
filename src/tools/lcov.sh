#!/bin/bash

set -euo pipefail

DRY_RUN=false

for arg in "$@"; do
  case "$arg" in
    --dry-run)
      DRY_RUN=true
      ;;
    -h|--help)
      cat <<EOF
Usage: $0 [--dry-run]

Options:
  --dry-run   Compile with coverage flags, reset counters, but do not run tests
  -h, --help  Show this help message
EOF
      exit 0
      ;;
    *)
      echo "Unknown option: $arg" >&2
      exit 1
      ;;
  esac
done

# Build tests with coverage enabled
make CXX="g++-11 --coverage -g" MODE=dev tests

# Remove old report files
rm -f coverage.info cov-info.txt

# Optional: remove stale gcda files so old runs don't contaminate coverage
find . -name '*.gcda' -delete

# Reset counters
lcov --gcov-tool gcov-11 --zerocounters --directory .

# Run all unit tests, unless 'dry' run, then only v. fast qip tests
if [ "$DRY_RUN" = true ]; then
  echo "Dry run: running small [qip] test subset."
  ./tests "[qip]"
else
  ./tests "[unit]"
fi

# Capture coverage
lcov --gcov-tool gcov-11 \
  --capture \
  --directory . \
  --output-file coverage.info

# Remove files we don't want in the final report
lcov --gcov-tool gcov-11 \
  --remove coverage.info \
  '*/catch2/*' \
  '*/fmt/*' \
  '*/json/*' \
  '/*.tests*' \
  '*/version/*' \
  '/usr/*' \
  '*/Modules/*' \
  '*/Kionisation/*' \
  '*/Physics/periodicTable.*' \
  '*/DiracODE/ComplexDirac.*' \
  '*/MBPT/Ladder.*' \
  '*/qip/Widgets.*' \
  '*/main.cpp' \
  '*/tests.cpp' \
  --output-file coverage.info

# Print summary
lcov --list coverage.info | tee cov-info.txt

# View output locally:
# https://lcov-viewer.netlify.app/
# https://github.com/eugenezinovyev/lcov-viewer