#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BOOK = ROOT / "references" / "Matrix_Computations_merged.md"
SECTIONS_DIR = ROOT / "references" / "sections"
SECTIONS_MANIFEST = SECTIONS_DIR / "README.md"

HEADING_RE = re.compile(r"^(#{1,6})\s+(.*\S)\s*$")
NUMBERED_TITLE_RE = re.compile(r"^(\d+(?:\.\d+)*)\b")
MANIFEST_ROW_RE = re.compile(r"^\| `([^`]+)` \| (\d+)-(\d+) \| (.+) \|$")
SECTION_REF_RE = re.compile(r"^§?(\d+(?:\.\d+){1,2})$")
PAREN_NUMBER_RE = re.compile(r"^\((\d+(?:\.\d+){2})\)$")
ITEM_NUMBER_RE = re.compile(r"^(\d+(?:\.\d+){2})$")
ITEM_KIND_RE = re.compile(
    r"^(theorem|lemma|corollary|algorithm|equation|eq|alg)\s*[: ]\s*\(?\s*(\d+(?:\.\d+){2})\s*\)?$",
    re.IGNORECASE,
)

ITEM_ALIASES = {
    "eq": "equation",
    "alg": "algorithm",
}


def load_lines() -> list[str]:
    if not BOOK.exists():
        print(f"error: missing book file: {BOOK}", file=sys.stderr)
        sys.exit(1)
    return BOOK.read_text(errors="ignore").splitlines()


def load_manifest() -> list[dict[str, object]]:
    if not SECTIONS_MANIFEST.exists():
        return []
    rows: list[dict[str, object]] = []
    for line in SECTIONS_MANIFEST.read_text(errors="ignore").splitlines():
        match = MANIFEST_ROW_RE.match(line.strip())
        if not match:
            continue
        rows.append(
            {
                "file": match.group(1),
                "start": int(match.group(2)),
                "end": int(match.group(3)),
                "kind": match.group(4),
            }
        )
    return rows


def build_headings(lines: list[str]) -> list[dict[str, object]]:
    headings: list[dict[str, object]] = []
    for lineno, line in enumerate(lines, start=1):
        match = HEADING_RE.match(line)
        if not match:
            continue
        title = match.group(2)
        numbered = NUMBERED_TITLE_RE.match(title)
        level = len(numbered.group(1).split(".")) if numbered else len(match.group(1))
        headings.append({"level": level, "title": title, "start": lineno})

    for i, heading in enumerate(headings):
        level = int(heading["level"])
        end = len(lines)
        for later in headings[i + 1 :]:
            if int(later["level"]) <= level:
                end = int(later["start"]) - 1
                break
        heading["end"] = end
    return headings


def normalize(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip().lower())


def best_manifest_match(start: int, end: int, manifest: list[dict[str, object]]) -> dict[str, object] | None:
    candidates = []
    for row in manifest:
        row_start = int(row["start"])
        row_end = int(row["end"])
        if row_start <= start and end <= row_end:
            span = row_end - row_start
            kind = str(row["kind"])
            penalty = 1 if kind in {"frontmatter", "index", "chapter-intro"} else 0
            candidates.append((penalty, span, row))
    if candidates:
        candidates.sort(key=lambda item: (item[0], item[1]))
        return candidates[0][2]

    start_candidates = []
    for row in manifest:
        row_start = int(row["start"])
        row_end = int(row["end"])
        if row_start <= start <= row_end:
            span = row_end - row_start
            kind = str(row["kind"])
            penalty = 1 if kind in {"frontmatter", "index", "chapter-intro"} else 0
            start_candidates.append((penalty, span, row))
    if not start_candidates:
        return None
    start_candidates.sort(key=lambda item: (item[0], item[1]))
    return start_candidates[0][2]


def format_match(start: int, end: int, manifest: list[dict[str, object]]) -> str:
    match = best_manifest_match(start, end, manifest)
    if match is None:
        return ""
    return f"  -> references/sections/{match['file']}"


def print_heading_matches(matches: list[dict[str, object]], manifest: list[dict[str, object]]) -> None:
    for match in matches:
        level = int(match["level"])
        indent = "  " * (level - 1)
        title = str(match["title"])
        start = int(match["start"])
        end = int(match["end"])
        print(f"{start:>6}-{end:<6} {indent}{title}{format_match(start, end, manifest)}")


def print_manifest_rows(rows: list[dict[str, object]]) -> None:
    for row in rows:
        print(f"references/sections/{row['file']}  ({row['start']}-{row['end']}, {row['kind']})")


def normalize_section_ref(text: str) -> str:
    match = SECTION_REF_RE.match(text.strip())
    if not match:
        print(f"error: invalid section reference: {text}", file=sys.stderr)
        sys.exit(1)
    return match.group(1)


def normalize_item_kind(kind: str) -> str:
    return ITEM_ALIASES.get(kind.lower(), kind.lower())


def normalize_item_number(text: str) -> str:
    stripped = text.strip()
    paren = PAREN_NUMBER_RE.match(stripped)
    if paren:
        return paren.group(1)
    plain = ITEM_NUMBER_RE.match(stripped)
    if plain:
        return plain.group(1)
    print(f"error: invalid item number: {text}", file=sys.stderr)
    sys.exit(1)


def extract_excerpt(line: str, width: int = 120) -> str:
    cleaned = re.sub(r"\s+", " ", line.strip())
    if len(cleaned) <= width:
        return cleaned
    return cleaned[: width - 3] + "..."


def find_item_matches(kind: str, number: str, lines: list[str]) -> list[dict[str, object]]:
    escaped = re.escape(number)
    if kind == "equation":
        primary = re.compile(rf"\\tag\s*\{{\s*{escaped}\s*\}}")
        fallback = re.compile(rf"Equation\s*\({escaped}\)")
    else:
        label = kind.capitalize()
        primary = re.compile(rf"^\s*{label}\s+{escaped}(?:\b|\.)")
        fallback = re.compile(rf"\b{label}\s+{escaped}(?:\b|\.)")

    matches: list[dict[str, object]] = []
    for lineno, line in enumerate(lines, start=1):
        if primary.search(line):
            matches.append({"line": lineno, "text": line, "quality": 0})
    if matches:
        return matches

    for lineno, line in enumerate(lines, start=1):
        if fallback.search(line):
            matches.append({"line": lineno, "text": line, "quality": 1})
    return matches


def print_item_matches(kind: str, number: str, matches: list[dict[str, object]], manifest: list[dict[str, object]]) -> None:
    for match in matches:
        lineno = int(match["line"])
        excerpt = extract_excerpt(str(match["text"]))
        print(f"{kind} {number} @ {lineno}: {excerpt}{format_match(lineno, lineno, manifest)}")


def cmd_toc(args: argparse.Namespace, headings: list[dict[str, object]], manifest: list[dict[str, object]]) -> None:
    query = normalize(args.query or "")
    matches = []
    for heading in headings:
        if int(heading["level"]) > args.max_level:
            continue
        title = normalize(str(heading["title"]))
        if query and query not in title:
            continue
        matches.append(heading)
    print_heading_matches(matches[: args.limit], manifest)


def cmd_heading(args: argparse.Namespace, headings: list[dict[str, object]], manifest: list[dict[str, object]]) -> None:
    tokens = [token for token in normalize(" ".join(args.query)).split(" ") if token]
    matches = []
    for heading in headings:
        title = normalize(str(heading["title"]))
        if all(token in title for token in tokens):
            matches.append(heading)
    print_heading_matches(matches[: args.limit], manifest)


def cmd_text(args: argparse.Namespace, lines: list[str], manifest: list[dict[str, object]]) -> None:
    query = normalize(" ".join(args.query))
    if not query:
        return
    shown = 0
    for lineno, line in enumerate(lines, start=1):
        if query in normalize(line):
            print(f"{lineno:>6}: {line}{format_match(lineno, lineno, manifest)}")
            shown += 1
            if shown >= args.limit:
                break


def cmd_section(args: argparse.Namespace, manifest: list[dict[str, object]]) -> None:
    ref = normalize_section_ref(args.reference)
    exact = [row for row in manifest if str(row["kind"]) == ref]
    if exact:
        print_manifest_rows(exact)
        return

    prefix_rows = [
        row
        for row in manifest
        if str(row["kind"]) == f"{ref} intro" or str(row["kind"]).startswith(f"{ref}.")
    ]
    if prefix_rows:
        prefix_rows.sort(key=lambda row: (0 if str(row["kind"]) == f"{ref} intro" else 1, int(row["start"])))
        print_manifest_rows(prefix_rows)
        return

    print(f"no split file found for section reference {ref}")


def cmd_item(kind: str, reference: str, lines: list[str], manifest: list[dict[str, object]]) -> None:
    normalized_kind = normalize_item_kind(kind)
    number = normalize_item_number(reference)
    matches = find_item_matches(normalized_kind, number, lines)
    if not matches:
        print(f"no {normalized_kind} found for reference {number}")
        return
    print_item_matches(normalized_kind, number, matches, manifest)


def cmd_ref(args: argparse.Namespace, lines: list[str], manifest: list[dict[str, object]]) -> None:
    reference = " ".join(args.reference).strip()
    if not reference:
        print("error: missing reference", file=sys.stderr)
        sys.exit(1)

    if SECTION_REF_RE.match(reference):
        cmd_section(argparse.Namespace(reference=reference), manifest)
        return

    parenthesized = PAREN_NUMBER_RE.match(reference)
    if parenthesized:
        cmd_item("equation", parenthesized.group(1), lines, manifest)
        return

    item_match = ITEM_KIND_RE.match(reference)
    if item_match:
        cmd_item(item_match.group(1), item_match.group(2), lines, manifest)
        return

    print(
        "error: could not parse reference. Use forms like '6.5', '§7.7.6', 'Theorem 2.4.1', 'Algorithm 5.1.1', or '(5.1.2)'.",
        file=sys.stderr,
    )
    sys.exit(1)


def cmd_file_at(args: argparse.Namespace, manifest: list[dict[str, object]]) -> None:
    match = best_manifest_match(args.line, args.line, manifest)
    if match is None:
        print("no containing split file found")
        return
    print_manifest_rows([match])


def main() -> None:
    parser = argparse.ArgumentParser(description="Locate sections and labeled results in Matrix Computations markdown")
    sub = parser.add_subparsers(dest="command", required=True)

    toc = sub.add_parser("toc", help="list headings, optionally filtered by a keyword")
    toc.add_argument("query", nargs="?", default="")
    toc.add_argument("--max-level", type=int, default=3)
    toc.add_argument("--limit", type=int, default=80)

    heading = sub.add_parser("heading", help="find headings whose titles contain all query words")
    heading.add_argument("query", nargs="+", help="heading words to match")
    heading.add_argument("--limit", type=int, default=20)

    text = sub.add_parser("text", help="find literal text matches in the book")
    text.add_argument("query", nargs="+", help="text to search for")
    text.add_argument("--limit", type=int, default=20)

    section = sub.add_parser("section", help="resolve a section reference like 6.5, §6.5, or 7.7.6")
    section.add_argument("reference")

    theorem = sub.add_parser("theorem", help="resolve a theorem reference like 2.4.1")
    theorem.add_argument("reference")

    lemma = sub.add_parser("lemma", help="resolve a lemma reference like 2.6.1")
    lemma.add_argument("reference")

    corollary = sub.add_parser("corollary", help="resolve a corollary reference like 2.4.2")
    corollary.add_argument("reference")

    algorithm = sub.add_parser("algorithm", help="resolve an algorithm reference like 5.1.1")
    algorithm.add_argument("reference")

    equation = sub.add_parser("equation", help="resolve an equation reference like 5.1.2")
    equation.add_argument("reference")

    ref = sub.add_parser("ref", help="resolve a generic internal reference like '6.5', 'Theorem 2.4.1', or '(5.1.2)'")
    ref.add_argument("reference", nargs="+")

    file_at = sub.add_parser("file-at", help="show which split file contains a merged-file line number")
    file_at.add_argument("line", type=int)

    args = parser.parse_args()
    lines = load_lines()
    manifest = load_manifest()
    headings = build_headings(lines)

    if args.command == "toc":
        cmd_toc(args, headings, manifest)
    elif args.command == "heading":
        cmd_heading(args, headings, manifest)
    elif args.command == "text":
        cmd_text(args, lines, manifest)
    elif args.command == "section":
        cmd_section(args, manifest)
    elif args.command == "theorem":
        cmd_item("theorem", args.reference, lines, manifest)
    elif args.command == "lemma":
        cmd_item("lemma", args.reference, lines, manifest)
    elif args.command == "corollary":
        cmd_item("corollary", args.reference, lines, manifest)
    elif args.command == "algorithm":
        cmd_item("algorithm", args.reference, lines, manifest)
    elif args.command == "equation":
        cmd_item("equation", args.reference, lines, manifest)
    elif args.command == "ref":
        cmd_ref(args, lines, manifest)
    elif args.command == "file-at":
        cmd_file_at(args, manifest)


if __name__ == "__main__":
    main()
