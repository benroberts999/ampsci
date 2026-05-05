#!/bin/bash
# Format files changed since last commit (git diff HEAD).
# clang-format <16 omits space before '{' in constructor bodies; normalise with sed.
# clang-format <16 changes '} //' alignment: normalise to exactly one space.
# Usage: clang_format.sh [clang-format-binary]

CLANG_FORMAT=${1:-clang-format}

which "$CLANG_FORMAT" > /dev/null 2>&1 || { echo "$CLANG_FORMAT not found, skipping"; exit 0; }
"$CLANG_FORMAT" --version

# BSD sed (macOS) requires an explicit empty-string argument with -i; GNU sed does not.
if sed --version 2>/dev/null | grep -q GNU; then
  SED_I=(sed -i)
else
  SED_I=(sed -i '')
fi

echo "Running clang-format (changed files only)"
# --force rewrites files even if they have unstaged changes
git clang-format --binary "$CLANG_FORMAT" --force HEAD || true

# Apply the same sed normalisations as clang_format_all.sh (see that script for rationale).
# Scope is limited to files changed since HEAD (same set git clang-format just touched).
git diff --name-only HEAD \
  | grep -E '\.(cpp|hpp|ipp)$' \
  | xargs -r "${SED_I[@]}" 's/){/) {/g'
git diff --name-only HEAD \
  | grep -E '\.(cpp|hpp|ipp)$' \
  | xargs -r "${SED_I[@]}" 's/}[[:space:]]*\/\//} \/\//g'
