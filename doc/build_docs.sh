#!/usr/bin/env bash
# Builds all ampsci documentation:
#   1. PDF manual from LaTeX source (doc/tex/)
#   2. HTML API reference using Doxygen (doc/doxygen/)
#   3. Post-processes the HTML: copies assets, patches search, merges API pages
#
# Can be run directly or via: make docs
# Safe to run from any directory.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DOXYGEN_DIR="$SCRIPT_DIR/doxygen"

# Build PDF from LaTeX; non-fatal if pdflatex not available
( cd "$ROOT/doc/tex" && make 2>/dev/null || : )
cp "$ROOT/doc/tex/ampsci.pdf" "$ROOT/doc/ampsci.pdf" 2>/dev/null || :

# Rebuild HTML from scratch
rm -rf "$ROOT/doc/html/"

# Run Doxygen; copy the top-level redirect page if it was generated
doxygen "$DOXYGEN_DIR/Doxyfile" 2>/dev/null \
    && cp "$DOXYGEN_DIR/ampsci.html" "$ROOT/doc/ampsci-documentation.html" 2>/dev/null \
    || true

# Copy PDF and images into the HTML tree so they are accessible from the docs
cp "$ROOT/doc/tex/ampsci.pdf" "$ROOT/doc/html/ampsci.pdf" 2>/dev/null || true
cp -r "$ROOT/doc/img" "$ROOT/doc/html/" 2>/dev/null || true

# Replace the default Doxygen search.js with the custom one (fixes short-query behaviour)
cp "$DOXYGEN_DIR/search.js" "$ROOT/doc/html/search/search.js"

# Inject searchdata.xml into search.html so the client-side search works offline
{
    echo '<script id="searchdata" type="text/xmldata">'
    cat "$ROOT/doc/searchdata.xml"
    echo '</script>'
} >> "$ROOT/doc/html/search.html"

# Merge hand-written API pages into the Doxygen HTML output
python3 "$DOXYGEN_DIR/merge_api.py" "$ROOT/doc/html" 2>/dev/null || true
