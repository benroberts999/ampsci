#!/bin/bash
# Format the entire project in place, then print which files were changed.
# Runs clang-format on every source file, then applies two sed normalisations to
# compensate for behaviour differences across clang-format versions:
#   - clang-format <16 omits the space in ){ -> ) {
#   - clang-format <16 collapses the space in } // -> } //
# Usage: clang_format_all.sh [clang-format-binary] [src-dir]

CLANG_FORMAT=${1:-clang-format}
SRC=${2:-src}

which "$CLANG_FORMAT" > /dev/null 2>&1 || { echo "$CLANG_FORMAT not found, skipping"; exit 0; }
"$CLANG_FORMAT" --version
echo "This will run project-wide clang-format in place. Suggest to commit first!"
read -p "Format entire project? [y/N] " ans && [ "$ans" = "y" ] || { echo "Aborted."; exit 1; }

# BSD sed (macOS) requires an explicit empty-string argument with -i; GNU sed does not.
if sed --version 2>/dev/null | grep -q GNU; then
  SED_I=(sed -i)
else
  SED_I=(sed -i '')
fi

# Snapshot checksums before any changes so we can report which files were modified,
# regardless of any pre-existing uncommitted changes in the working tree.
TMPFILE=$(mktemp)
find "$SRC" -type f \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print \
  | sort | xargs shasum > "$TMPFILE"

echo "Running clang-format (whole project)"
find "$SRC" \
  -type f \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print \
  | xargs "$CLANG_FORMAT" -i -verbose
echo 'Forcing single space between ) { for consistency across CF versions'
find "$SRC" \
  -type f \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print \
  | xargs "${SED_I[@]}" 's/){/) {/g'
echo 'Forcing single space between } // for consistency across CF versions'
find "$SRC" \
  -type f \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print \
  | xargs "${SED_I[@]}" 's/}[[:space:]]*\/\//} \/\//g'

# Recompute checksums and diff against the snapshot; '>' lines are files whose
# hash changed (i.e. were modified by clang-format or sed above).
echo ''
CHANGED=$(find "$SRC" -type f \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print \
  | sort | xargs shasum | diff "$TMPFILE" - | grep '^>' | awk '{print $NF}')
rm "$TMPFILE"
if [ -n "$CHANGED" ]; then
  echo "Files modified by clang-format:"
  echo "$CHANGED" | sed 's/^/  /'
else
  echo "No files were modified."
fi
