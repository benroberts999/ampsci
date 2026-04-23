#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

rm -rf "$ROOT/doc/html/"

doxygen "$SCRIPT_DIR/Doxyfile" 2>/dev/null \
    && cp "$SCRIPT_DIR/ampsci.html" "$ROOT/doc/ampsci-documentation.html" 2>/dev/null \
    || true

cp "$ROOT/doc/tex/ampsci.pdf" "$ROOT/doc/html/ampsci.pdf" 2>/dev/null || true
cp -r "$ROOT/doc/img" "$ROOT/doc/html/" 2>/dev/null || true
cp "$SCRIPT_DIR/search.js" "$ROOT/doc/html/search/search.js"

{
    echo '<script id="searchdata" type="text/xmldata">'
    cat "$ROOT/doc/searchdata.xml"
    echo '</script>'
} >> "$ROOT/doc/html/search.html"

python3 "$SCRIPT_DIR/merge_api.py" "$ROOT/doc/html" 2>/dev/null || true
