---
name: marker-pdf
description: |
  Trigger when the user mentions a local filesystem path ending in .pdf.
  Do not trigger for URLs (http/https), arXiv links, or when no local path is provided.
---

Goal: Prefer LaTeX source (best math fidelity). If unavailable, fall back to Marker extraction.

Sandbox-safe dirs (always):
OUT_ROOT="/tmp/codex_papers"
mkdir -p "$OUT_ROOT"

Step 1 — Identify arXiv ID (fast):
1) Extract PDF path(s).
2) Verify: `test -f "<PDF>"`
3) Extract first ~2 pages of text quickly:
   - Prefer: `pdftotext -f 1 -l 2 -layout "<PDF>" "$OUT_ROOT/head.txt"`
   - If pdftotext missing, use python pypdf:
     python3 - <<'PY'
     from pypdf import PdfReader
     import sys
     r=PdfReader(sys.argv[1])
     txt=[]
     for i in range(min(2, len(r.pages))):
         txt.append(r.pages[i].extract_text() or "")
     print("\n".join(txt))
     PY "<PDF>" > "$OUT_ROOT/head.txt"

4) Grep for arXiv IDs:
   - New-style: \b(arXiv:)?\d{4}\.\d{4,5}(v\d+)?\b
   - Old-style: \b(arXiv:)?[a-z\-]+/\d{7}(v\d+)?\b
   Example: `rg -o '(arXiv:)?\d{4}\.\d{4,5}(v\d+)?|(?:arXiv:)?[a-z\-]+/\d{7}(v\d+)?' "$OUT_ROOT/head.txt" | head -n 1`

If an arXiv ID is found:
Step 2 — Download TeX source from arXiv (best math):
5) Normalize ID (strip "arXiv:" prefix; keep version if present).
6) Download source tarball (use /src which always gives tar.gz):
   SRC_DIR="$OUT_ROOT/arxiv_src/<ID>"
   mkdir -p "$SRC_DIR"
   curl -L --fail "https://arxiv.org/src/<ID>" -o "$SRC_DIR/src.tar.gz"
7) Extract:
   tar -xzf "$SRC_DIR/src.tar.gz" -C "$SRC_DIR"
8) Find main TeX:
   - List .tex files: `find "$SRC_DIR" -maxdepth 3 -name '*.tex'`
   - Heuristic for main file: contains \documentclass OR is the largest .tex.
9) Read main .tex + included files needed for the asked section.
10) Use TeX as source-of-truth (math is exact).

If no arXiv ID or source fetch fails:
Step 3 — Fallback to Marker (two-pass):
11) Run marker_single without --force_ocr (Pass 1).
12) If math is degraded, rerun with --force_ocr on a narrow --page_range only.

Marker execution must use sandbox-safe caches/fonts:
CACHE_DIR="$OUT_ROOT/marker_cache"
FONT_DIR="$CACHE_DIR/fonts"
mkdir -p "$CACHE_DIR" "$FONT_DIR"
ENV:
XDG_CACHE_HOME="$CACHE_DIR"
HF_HOME="$CACHE_DIR/hf"
TRANSFORMERS_CACHE="$CACHE_DIR/transformers"
FONT_DIR="$FONT_DIR"
FONT_PATH="$FONT_DIR/GoNotoCurrent-Regular.ttf"

Never claim success unless commands exit 0 and outputs exist.