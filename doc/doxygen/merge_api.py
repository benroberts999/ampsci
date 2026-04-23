#!/usr/bin/env python3
"""
Post-processes Doxygen output to merge namespaces.html and annotated.html
into a single page (annotated.html), then redirects namespaces.html there.

namespaces.html has all namespaces (including those with no classes).
annotated.html has all classes, including top-level ones not in any namespace.
The merged page shows everything.
"""

import re
import sys
import os


AUTO_EXPAND_NAVTREE = True  # expand sidebar page-nodes on click; toggle to disable
HIDE_UNDOC_NAMESPACES = True  # hide namespaces with no doc comment from Full API list


def extract_contents_div(html):
    """Return the inner HTML of <div class="contents">...</div>."""
    m = re.search(r'<div class="contents">(.*?)\n</div>', html, re.DOTALL)
    if not m:
        raise ValueError('Could not find <div class="contents"> block')
    return m.group(1)


def extract_directory_table(contents_html):
    """Return the <table class="directory">...</table> block."""
    m = re.search(r'(<table class="directory">.*?</table>)', contents_html, re.DOTALL)
    if not m:
        raise ValueError('Could not find <table class="directory">')
    return m.group(1)


def extract_levels_block(contents_html):
    """Return the [detail level ...] div, if present."""
    m = re.search(r'(<div class="levels">.*?</div>)', contents_html, re.DOTALL)
    return m.group(1) if m else ""


def filter_undoc_namespaces(table_html):
    """Remove top-level undocumented namespace rows and all their children."""
    rows = re.findall(r"<tr[^>]*>.*?</tr>", table_html, re.DOTALL)

    undoc_idx = set()
    for row in rows:
        id_m = re.search(r'id="row_(\d+)_"', row)  # top-level only: row_N_
        if not id_m:
            continue
        icon_m = re.search(r'<span class="icon">(.)</span>', row)
        if not icon_m or icon_m.group(1) != "N":
            continue
        desc_m = re.search(r'<td class="desc">(.*?)</td>', row, re.DOTALL)
        if desc_m and not desc_m.group(1).strip():
            undoc_idx.add(id_m.group(1))

    if not undoc_idx:
        return table_html

    kept = []
    for row in rows:
        id_m = re.search(r'id="row_(\d+)_', row)  # first index in row_N[_M...]_
        if id_m and id_m.group(1) in undoc_idx:
            continue
        kept.append(row)

    return re.sub(
        r'(<table class="directory">)(.*?)(</table>)',
        lambda m: m.group(1) + "\n" + "\n".join(kept) + "\n" + m.group(3),
        table_html,
        flags=re.DOTALL,
    )


def get_top_level_rows(table_html):
    """
    Return rows that are top-level entries (id="row_N_" with a single number)
    and are NOT namespaces (icon text != 'N').
    These are classes/structs not inside any namespace.
    """
    # Match complete <tr>...</tr> blocks
    rows = re.findall(r"<tr[^>]*>.*?</tr>", table_html, re.DOTALL)
    result = []
    for row in rows:
        id_m = re.search(r'id="(row_\d+_)"', row)
        if not id_m:
            continue
        # Top-level: id has exactly one number segment, e.g. row_3_
        rid = id_m.group(1)
        if re.fullmatch(r"row_\d+_", rid):
            icon_m = re.search(r'<span class="icon">(.)</span>', row)
            if icon_m and icon_m.group(1) != "N":
                result.append(row)
    return result


def max_top_level_index(table_html):
    """Return the highest top-level row index (the N in row_N_)."""
    indices = [
        int(m)
        for m in re.findall(r'id="row_(\d+)_"', table_html)
        if re.fullmatch(r"\d+", m)
    ]
    return max(indices) if indices else -1


def reindex_row(row_html, old_index, new_index):
    """Replace all occurrences of the old top-level index with the new one."""
    # Replace in id, onclick, and style attributes that reference this index
    row_html = row_html.replace(f'"row_{old_index}_"', f'"row_{new_index}_"')
    row_html = row_html.replace(f"'row_{old_index}_'", f"'row_{new_index}_'")
    row_html = row_html.replace(
        f"toggleFolder('{old_index}_')", f"toggleFolder('{new_index}_')"
    )
    row_html = row_html.replace(f'"arr_{old_index}_"', f'"arr_{new_index}_"')
    row_html = row_html.replace(f"'arr_{old_index}_'", f"'arr_{new_index}_'")
    return row_html


def build_redirect(target_url, title="Redirecting..."):
    return f"""<!DOCTYPE html>
<html>
<head>
<meta http-equiv="refresh" content="0; url={target_url}" />
<title>{title}</title>
</head>
<body>
<p>This page has moved. <a href="{target_url}">Click here</a> if not redirected.</p>
</body>
</html>
"""


def merge(html_dir):
    ns_path = os.path.join(html_dir, "namespaces.html")
    cl_path = os.path.join(html_dir, "annotated.html")

    with open(ns_path) as f:
        ns_html = f.read()
    with open(cl_path) as f:
        cl_html = f.read()

    ns_contents = extract_contents_div(ns_html)
    cl_contents = extract_contents_div(cl_html)

    ns_table = extract_directory_table(ns_contents)
    if HIDE_UNDOC_NAMESPACES:
        ns_table = filter_undoc_namespaces(ns_table)
    cl_table = extract_directory_table(cl_contents)

    # Find top-level class rows in annotated.html not already covered by namespaces
    extra_rows = get_top_level_rows(cl_table)

    if extra_rows:
        next_idx = max_top_level_index(ns_table) + 1
        renumbered = []
        for i, row in enumerate(extra_rows):
            old_m = re.search(r'id="row_(\d+)_"', row)
            if old_m:
                old_idx = int(old_m.group(1))
                row = reindex_row(row, old_idx, next_idx + i)
            renumbered.append(row)

        merged_rows = "\n".join(renumbered)
        merged_table = ns_table.replace("</table>", merged_rows + "\n</table>", 1)
    else:
        merged_table = ns_table

    levels = extract_levels_block(ns_contents)
    merged_contents = (
        '<div class="textblock">Complete API reference: '
        "Documents all ampsci namespaces, classes, and functions. "
        "Namespaces group related functionality; expand each to browse its members, or use search to find something specific.</div>"
        '<div class="directory">\n' + levels + "\n" + merged_table + "\n</div>"
    )

    # Replace the contents div in namespaces.html (use it as the merged page)
    replacement = '<div class="contents">\n' + merged_contents + "\n</div>"
    new_ns_html = re.sub(
        r'<div class="contents">.*?\n</div>',
        lambda _: replacement,
        ns_html,
        flags=re.DOTALL,
        count=1,
    )

    with open(ns_path, "w") as f:
        f.write(new_ns_html)

    with open(cl_path, "w") as f:
        f.write(build_redirect("namespaces.html", "Namespaces & Classes"))

    if AUTO_EXPAND_NAVTREE:
        navtree_js_path = os.path.join(html_dir, "navtree.js")
        if os.path.exists(navtree_js_path):
            with open(navtree_js_path) as f:
                navtree_js = f.read()
            # Auto-expand current node's children on page load (like index/pages do)
            navtree_js = navtree_js.replace(
                'if (rootBase=="index" || rootBase=="pages" || rootBase=="search") {',
                'if (rootBase=="index" || rootBase=="pages" || rootBase=="search" || n.childrenData) {',
            )
            # Toggle collapse when clicking a link to the already-current page
            navtree_js = navtree_js.replace(
                "a.onclick = function() { storeLink(link); }",
                "a.onclick = function() { storeLink(link); if (link!=='index.html' && window.location.pathname.endsWith(link) && node.expandToggle) { node.expandToggle.onclick(); return false; } }",
            )
            with open(navtree_js_path, "w") as f:
                f.write(navtree_js)

    # Remove the "Class List" (annotated.html) duplicate from the sidebar navtree
    navtree_path = os.path.join(html_dir, "navtreedata.js")
    if os.path.exists(navtree_path):
        with open(navtree_path) as f:
            navtree = f.read()
        # Strip any child entry linking to annotated.html inside the Full API block
        navtree = re.sub(
            r',?\s*\[\s*"[^"]*",\s*"annotated\.html",\s*"[^"]*"\s*\]',
            "",
            navtree,
        )
        # Remove children from the "Full API" navtree entry so it's a plain link
        # Pattern accounts for one level of bracket nesting (no DOTALL greedy trap)
        navtree = re.sub(
            r'(\[\s*"Full API",\s*"namespaces\.html",\s*)\[[^\[\]]*(?:\[[^\[\]]*\][^\[\]]*)*\]',
            r"\1null",
            navtree,
        )
        with open(navtree_path, "w") as f:
            f.write(navtree)
        print("Removed annotated.html entry from navtreedata.js")

    print(f"Merged: {len(extra_rows)} top-level class rows added to namespaces.html")
    print(f"Redirected: annotated.html -> namespaces.html")


if __name__ == "__main__":
    html_dir = sys.argv[1] if len(sys.argv) > 1 else "doc/html"
    merge(html_dir)
