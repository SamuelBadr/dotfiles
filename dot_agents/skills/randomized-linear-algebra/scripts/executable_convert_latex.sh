#!/usr/bin/env bash
# Convert Martinsson & Tropp (arXiv:2002.01387) LaTeX source to markdown.
#
# Usage: convert_latex.sh <input.tex> <output.md>
#
# Pipeline:
#   1. sed preprocessor: rewrite old-style _{\rm X} / ^{\rm X} as _{\mathrm{X}}
#      (pandoc's math parser rejects \rm inside sub/superscripts; this is the
#      only class of conversion warning in the source).
#   2. pandoc 3.10: expand custom macros (\mtx -> \bm{\mathsf{}}, etc.),
#      render math as $...$ / $$...$$, emit markdown headings. The algorithm
#      (float) + algorithmic (pseudocode) packages are opaque to pandoc, so
#      \begin{algorithm} blocks vanish from the output.
#   3. convert_algorithms.py: extract the 18 dropped algorithm environments
#      from the LaTeX, parse the algorithmic body to indented pseudocode, and
#      emit ::: {#alg:label .algorithm} divs. These are spliced back into the
#      pandoc output immediately after the paragraph that first references each
#      algorithm, so the canonical pseudocode sits where a reader expects it.
#   4. split_sections.py: split the merged markdown at top-level headings into
#      references/sections/NN-slug.md and regenerate the section manifest.
#
# Pinned to pandoc 3.10 (recorded in MT_merged.md provenance header). A major
# pandoc version bump may shift output; re-verify the addressing-layer regexes
# in book_lookup.py after any upgrade.
set -euo pipefail

if [ $# -ne 2 ]; then
  echo "usage: $0 <input.tex> <output.md>" >&2
  exit 64
fi

in="$1"; out="$2"
here="$(cd "$(dirname "$0")" && pwd)"

# Step 1+2: preprocess \rm, then pandoc convert to a temp file.
tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT
sed -E 's/([_^])\{\\rm ([A-Za-z]+)\}/\1{\\mathrm{\2}}/g' "$in" \
  | pandoc -f latex -t markdown -o "$tmp"

# Step 3: splice extracted algorithm blocks back in after their first reference.
python3 "$here/convert_algorithms.py" "$in" > "$tmp.algs"

python3 - "$tmp" "$tmp.algs" "$out" <<'PY'
import sys, re

with open(sys.argv[1]) as f:
    md = f.read()
with open(sys.argv[2]) as f:
    algs = f.read()

# Each algorithm block starts with "::: {#alg:LABEL .algorithm}". Splice it
# into the markdown right after the line containing 'reference="alg:LABEL"}'
# (the pandoc-rendered cross-reference, which spans two lines ending in '}').
md_lines = md.splitlines(keepends=True)
inserted = set()
for block in re.split(r"(?m)^(?=^::: \{#alg:)", algs):
    m = re.match(r"::: \{#(alg:[^}\s]+)", block)
    if not m:
        continue
    label = m.group(1)
    if label in inserted:
        continue
    marker = f'reference="{label}"}}'
    splice_at = None
    for i, line in enumerate(md_lines):
        if marker in line:
            splice_at = i + 1
            break
    if splice_at is None:
        # No cross-reference found; append at end so the block is not lost.
        md_lines.append("\n" + block)
    else:
        md_lines.insert(splice_at, block + "\n")
    inserted.add(label)

with open(sys.argv[3], "w") as f:
    f.writelines(md_lines)
PY

echo "converted $in -> $out"

# Step 4: split into per-section files + regenerate the manifest.
python3 "$here/split_sections.py" "$out" "$(dirname "$out")/sections"
