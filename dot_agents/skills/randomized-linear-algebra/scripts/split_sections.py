#!/usr/bin/env python3
"""Split MT_merged.md into per-section files and regenerate the section manifest.

convert_latex.sh produces references/MT_merged.md (pandoc output + spliced
algorithm divs). This script splits it at top-level ('# ') headings into
references/sections/NN-slug.md files and writes references/sections/README.md
(the manifest that book_lookup.py loads).

Section filenames are fixed (by source order) so the addressing layer stays
stable across re-conversions; only the line ranges in the manifest shift.

Usage: split_sections.py <merged.md> <sections_dir>
"""
from __future__ import annotations

import os
import re
import sys

# Fixed filename order (matches the survey's section order). The Kind column is
# read from each heading's {#sec:...} attribute (or "intro"/"frontmatter").
FILENAMES = [
    "00-frontmatter.md",
    "01-introduction.md",
    "02-linear-algebra-preliminaries.md",
    "03-probability-preliminaries.md",
    "04-trace-estimation-by-sampling.md",
    "05-schatten-norm-estimation-by-sampling.md",
    "06-maximum-eigenvalues-and-trace-functions.md",
    "07-matrix-approximation-by-sampling.md",
    "08-randomized-embeddings.md",
    "09-structured-random-embeddings.md",
    "10-how-to-use-random-embeddings.md",
    "11-the-randomized-rangefinder.md",
    "12-error-estimation-and-adaptivity.md",
    "13-finding-natural-bases-qr-id-and-cur.md",
    "14-nystrom-approximation.md",
    "15-single-view-algorithms.md",
    "16-factoring-matrices-of-full-or-nearly-full-rank.md",
    "17-general-linear-solvers.md",
    "18-linear-solvers-for-graph-laplacians.md",
    "19-kernel-matrices-in-machine-learning.md",
    "20-high-accuracy-approximation-of-kernel-matrices.md",
]


def main() -> None:
    if len(sys.argv) != 3:
        print("usage: split_sections.py <merged.md> <sections_dir>", file=sys.stderr)
        sys.exit(64)
    merged_path, sections_dir = sys.argv[1], sys.argv[2]
    with open(merged_path) as f:
        lines = f.readlines()  # keep newlines; 1-indexed conceptually

    # Find top-level heading line numbers (1-indexed) and their {#sec:...} label.
    heading_re = re.compile(r"^# .+?(?:\{(#(sec:[^}]+))\})?\s*$")
    headings: list[tuple[int, str | None]] = []  # (lineno, label or None)
    for i, line in enumerate(lines, start=1):
        m = heading_re.match(line.rstrip("\n"))
        if m:
            headings.append((i, m.group(2)))  # group(2) is 'sec:...' or None

    if len(headings) + 1 != len(FILENAMES):
        print(
            f"error: found {len(headings)} headings but expected {len(FILENAMES) - 1}; "
            "section count changed -- update FILENAMES in split_sections.py",
            file=sys.stderr,
        )
        sys.exit(1)

    # Section i spans [start_i, start_{i+1} - 1]; frontmatter spans [1, first_heading-1].
    starts = [1] + [h[0] for h in headings]
    ends = [h[0] - 1 for h in headings] + [len(lines)]
    kinds: list[str] = ["frontmatter"]
    for _, label in headings:
        kinds.append(label if label else "intro")

    os.makedirs(sections_dir, exist_ok=True)
    rows: list[str] = []
    for idx, fname in enumerate(FILENAMES):
        s, e, kind = starts[idx], ends[idx], kinds[idx]
        with open(os.path.join(sections_dir, fname), "w") as f:
            f.writelines(lines[s - 1 : e])
        rows.append(f"| `{fname}` | {s}-{e} | {kind} |")

    # Manifest header (matches MANIFEST_ROW_RE in book_lookup.py).
    manifest = [
        "# Section split manifest\n",
        "\n",
        "Source: Martinsson & Tropp, *Randomized Numerical Linear Algebra* (arXiv:2002.01387v3).\n",
        "Converted via `scripts/convert_latex.sh` (pandoc 3.10 + algorithm extractor)\n",
        "and split by `scripts/split_sections.py`. Line ranges refer to `../MT_merged.md`.\n",
        "The `Kind` column holds the section label (`sec:...`) used by `book_lookup.py section`.\n",
        "\n",
        "| File | Lines | Kind |\n",
        "| --- | --- | --- |\n",
    ]
    for r in rows:
        manifest.append(r + "\n")
    with open(os.path.join(sections_dir, "README.md"), "w") as f:
        f.writelines(manifest)

    print(f"split {merged_path} -> {len(FILENAMES)} files in {sections_dir}")


if __name__ == "__main__":
    main()
