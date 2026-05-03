#!/bin/bash
# Format the entire project in place.
# clang-format <16 omits space before '{' in constructor bodies; normalise with sed.
# clang-format <16 changes '} //' alignment: normalise to exactly one space.
# Usage: clang_format_all.sh [clang-format-binary] [src-dir]

CLANG_FORMAT=${1:-clang-format}
SRC=${2:-src}

which "$CLANG_FORMAT" > /dev/null 2>&1 || { echo "$CLANG_FORMAT not found, skipping"; exit 0; }
"$CLANG_FORMAT" --version
echo "This will run project-wide clang-format in place. Suggest to commit first!"
read -p "Format entire project? [y/N] " ans && [ "$ans" = "y" ] || { echo "Aborted."; exit 1; }
echo "Running clang-format (whole project)"
find "$SRC" \
  -type f \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print \
  | xargs -r "$CLANG_FORMAT" -i -verbose
echo 'Forcing single space between ) { for consistency across CF versions'
find "$SRC" \
  -type f \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print \
  | xargs -r sed -i 's/){/) {/g'
echo 'Forcing single space between } // for consistency across CF versions'
find "$SRC" \
  -type f \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print \
  | xargs -r sed -i 's/}[[:space:]]*\/\//} \/\//g'
