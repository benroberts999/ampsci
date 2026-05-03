#!/bin/bash
# Format files changed since last commit (git diff HEAD).
# clang-format <16 omits space before '{' in constructor bodies; normalise with sed.
# clang-format <16 changes '} //' alignment: normalise to exactly one space.
# Usage: clang_format.sh [clang-format-binary]

CLANG_FORMAT=${1:-clang-format}

which "$CLANG_FORMAT" > /dev/null 2>&1 || { echo "$CLANG_FORMAT not found, skipping"; exit 0; }
"$CLANG_FORMAT" --version
echo "Running clang-format (changed files only)"
git clang-format --binary "$CLANG_FORMAT" --force HEAD || true
git diff --name-only HEAD \
  | grep -E '\.(cpp|hpp|ipp)$' \
  | xargs -r sed -i 's/){/) {/g'
git diff --name-only HEAD \
  | grep -E '\.(cpp|hpp|ipp)$' \
  | xargs -r sed -i 's/}[[:space:]]*\/\//} \/\//g'
