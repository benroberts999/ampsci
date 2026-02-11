#!/usr/bin/env bash
set -euo pipefail

SRC_DIR="src"

# Process directories deepest first
mapfile -t dirs < <(
  find "$SRC_DIR" -type d -print |
  awk '{ print length, $0 }' |
  sort -nr |
  cut -d" " -f2-
)

for dir in "${dirs[@]}"; do
  out="$dir/include.hpp"

  # Headers in THIS directory only (not recursive), excluding include.hpp
  mapfile -t headers < <(
    find "$dir" -maxdepth 1 -type f -name '*.hpp' ! -name 'include.hpp' -print | sort
  )

  # Immediate child directories that already have include.hpp
  mapfile -t children < <(
    find "$dir" -mindepth 1 -maxdepth 1 -type d \
      -exec test -f "{}/include.hpp" \; -print | sort
  )

  # If no headers AND no child includes → skip
  if [ ${#headers[@]} -eq 0 ] && [ ${#children[@]} -eq 0 ]; then
    # Remove stale include.hpp if it exists
    [ -f "$out" ] && rm -f "$out"
    continue
  fi

  {
    echo "#pragma once"

    # Include headers in this directory
    for h in "${headers[@]}"; do
      rel="${h#"$SRC_DIR/"}"
      echo "#include \"${rel}\""
    done

    # Include child include.hpp files
    for c in "${children[@]}"; do
      rel="${c#"$SRC_DIR/"}/include.hpp"
      echo "#include \"${rel}\""
    done
  } > "$out"

  echo "Wrote $out"
done

echo "Done."
