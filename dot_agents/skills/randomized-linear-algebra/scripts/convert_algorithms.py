#!/usr/bin/env python3
"""Extract algorithm environments from the LaTeX source and render them as
markdown divs that pandoc otherwise drops.

The `algorithm` (float) + `algorithmic` (pseudocode) packages are opaque to
pandoc, so \\begin{algorithm}...\\end{algorithm} blocks vanish entirely from the
converted markdown. This script extracts each block, parses the algorithmic
body into indented pseudocode, and emits a `::: {#label .algorithm}` div with
the caption as a header and the body as a fenced code block.

Inline math in the pseudocode (e.g. \\mtx{A}, \\vct{\\omega}) is passed through
unchanged; pandoc expands the macros when it processes the spliced output. The
div anchor and label let `book_lookup.py algorithm <label>` resolve identically
to theorem/lemma lookups.

Usage: convert_algorithms.py <input.tex>
Writes extracted algorithm blocks to stdout, one per algorithm, in source order.
"""
from __future__ import annotations

import re
import sys


def parse_brace_arg(line: str, start: int) -> tuple[str | None, str | None]:
    """Find the first '{' at or after `start` and return its balanced-brace
    argument plus whatever follows the matching '}'. Tolerates the tab/space
    the source places between a command keyword and its argument (e.g.
    '\\For\\t{...}', '\\If\\t\\t{...}') and nested braces inside the argument
    (e.g. '\\abs{\\xi_{i-1}}'). Returns (None, None) if no '{' or unbalanced.
    """
    i = line.find("{", start)
    if i < 0:
        return None, None
    depth = 0
    for j in range(i, len(line)):
        if line[j] == "{":
            depth += 1
        elif line[j] == "}":
            depth -= 1
            if depth == 0:
                return line[i + 1 : j], line[j + 1 :].strip()
    return None, None  # unbalanced


def clean(text: str | None) -> str:
    """Strip the trailing LaTeX line-continuation marker (%) and whitespace."""
    if not text:
        return ""
    return re.sub(r"\s*%$", "", text).strip()


_UMLAUTS = {"o": "\u00f6", "a": "\u00e4", "u": "\u00fc", "O": "\u00d6", "A": "\u00c4", "U": "\u00dc"}


def unwrap_accents(text: str) -> str:
    """Expand German umlaut accents ('\"o', '\"{o}', '{\"o}') to precomposed
    chars. Only the forms that appear in the source; other accents pass through."""
    text = re.sub(r'\{\\"([aAuUoO])\}', lambda m: _UMLAUTS[m.group(1)], text)
    text = re.sub(r'\\"\\{([aAuUoO])\}', lambda m: _UMLAUTS[m.group(1)], text)
    return re.sub(r'\\"([aAuUoO])(?![a-zA-Z])', lambda m: _UMLAUTS[m.group(1)], text)


def unwrap_text(text: str) -> str:
    """Unwrap pure-formatting macros (\\texttt{X}, \\textbf{X}, \\textsc{X},
    \\textit{X}) to their content X, unescape \\_ to _, and expand umlauts.
    Math macros (\\mtx, \\vct, \\norm, ...) pass through verbatim -- they are
    readable shorthand documented in SKILL.md."""
    for cmd in ("texttt", "textbf", "textsc", "textit"):
        text = re.sub(r"\\" + cmd + r"\{([^{}]*)\}", r"\1", text)
    text = text.replace(r"\_", "_")
    return unwrap_accents(text)


def split_comment(text: str) -> tuple[str, str]:
    """Split text at the first inline \\Comment command into (main, comment).
    \\Comment never appears inside $...$ math in this source, so a plain search
    is safe. Returns (text, '') when there is no inline comment."""
    m = re.search(r"\\Comment\b", text)
    if m:
        return text[:m.start()], clean(text[m.end():])
    return text, ""


def extract_algorithms(tex: str) -> list[dict[str, str]]:
    """Return a list of {label, caption, body_md} for each algorithm environment."""
    alg_re = re.compile(r"\\begin\{algorithm\}.*?\\end\{algorithm\}", re.DOTALL)
    results: list[dict[str, str]] = []
    for alg in alg_re.findall(tex):
        # \caption{\textit{TITLE.} ...}: balanced-brace match the caption arg,
        # then the \textit{...} title inside it, so accented chars like {\"o}
        # (Nystrom) and trailing description text are not truncated by a naive
        # [^}]+ regex that stops at the first '}'.
        cap_i = alg.find(r"\caption{")
        title = None
        if cap_i >= 0:
            cap_arg, _ = parse_brace_arg(alg, cap_i + len(r"\caption"))
            if cap_arg is not None:
                ti_i = cap_arg.find(r"\textit{")
                if ti_i >= 0:
                    title, _ = parse_brace_arg(cap_arg, ti_i + len(r"\textit"))
        label_m = re.search(r"\\label\{(alg:[^}]+)\}", alg)
        body_m = re.search(r"\\begin\{algorithmic\}(?:\[\d\])?(.*?)\\end\{algorithmic\}", alg, re.DOTALL)
        if not (title and label_m and body_m):
            continue
        results.append(
            {
                "label": label_m.group(1),
                "caption": unwrap_accents(title.rstrip(".")),
                # caption + label sit at the top of the algorithmic body;
                # drop everything up to and including the \label line so only
                # the pseudocode reaches render_body.
                "body": re.split(r"\\label\{alg:[^}]+\}", body_m.group(1), maxsplit=1)[-1],
            }
        )
    return results


def render_body(body: str) -> str:
    """Convert algorithmic commands to indented pseudocode.

    Inline math (\\mtx{A}, \\F^n, etc.) passes through; pandoc expands macros
    when it processes the spliced output. The control commands we handle:
    Require/Ensure (I/O), Function, For/While/If/ElsIf/Else, State, Statex
    (blank line), Comment, Return, and the End* delimiters that close blocks.

    Block-header arguments are parsed with balanced-brace matching so that
    whitespace between the keyword and '{' (the source uses '\\For\\t{...}',
    '\\If\\t\\t{...}') and nested braces inside conditions (e.g.
    '\\abs{\\xi_{i-1}}') are handled correctly; a prior regex approach silently
    dropped any header it failed to match. Trailing content on a header line
    (e.g. '\\If{cond} \\textbf{break}') is emitted as the first statement inside
    the new block. An inline \\Comment on a \\State or header line is split off
    and rendered as a '# ...' comment.
    """
    indent = 0
    out: list[str] = []

    def emit(text: str) -> None:
        out.append("    " * indent + unwrap_text(text))

    def emit_rest(rest: str | None) -> None:
        """Emit trailing content after a block header: a statement (if any) and
        a trailing inline comment (as its own '# ...' line)."""
        rest = clean(rest)
        if not rest:
            return
        main, comment = split_comment(rest)
        if main.strip():
            emit(main.strip())
        if comment:
            emit(f"# {comment}")

    for raw in body.splitlines():
        line = raw.strip()
        if not line or line == "%":
            continue

        if line.startswith(r"\Require"):
            emit(f"**Input:** {clean(line[len(r'\Require'):])}")
        elif line.startswith(r"\Ensure"):
            emit(f"**Output:** {clean(line[len(r'\Ensure'):])}")
        elif line.startswith(r"\Function"):
            name, rest = parse_brace_arg(line, len(r"\Function"))
            if name is not None and rest is not None:
                args, _ = parse_brace_arg(rest, 0)
                emit(f"function {name}({args or ''}):")
                indent += 1
        elif line.startswith(r"\For"):
            cond, rest = parse_brace_arg(line, len(r"\For"))
            if cond is not None:
                emit(f"for {clean(cond)}:")
                indent += 1
                emit_rest(rest)
        elif line.startswith(r"\While"):
            cond, rest = parse_brace_arg(line, len(r"\While"))
            if cond is not None:
                emit(f"while {clean(cond)}:")
                indent += 1
                emit_rest(rest)
        elif line.startswith(r"\If"):
            cond, rest = parse_brace_arg(line, len(r"\If"))
            if cond is not None:
                emit(f"if {clean(cond)}:")
                indent += 1
                emit_rest(rest)
        elif line.startswith(r"\ElsIf"):
            cond, _ = parse_brace_arg(line, len(r"\ElsIf"))
            indent = max(0, indent - 1)
            if cond is not None:
                emit(f"elif {clean(cond)}:")
                indent += 1
        elif line.startswith(r"\Else"):
            indent = max(0, indent - 1)
            emit("else:")
            indent += 1
        elif (
            line.startswith(r"\EndFor")
            or line.startswith(r"\EndFunction")
            or line.startswith(r"\EndIf")
            or line.startswith(r"\EndWhile")
        ):
            indent = max(0, indent - 1)
        elif line.startswith(r"\label"):
            # \label{line:...} inside a body is a cross-ref anchor for a source
            # line, not pseudocode; drop it. (The \label{alg:...} caption label
            # is stripped earlier in extract_algorithms.)
            continue
        elif line.startswith(r"\Statex"):
            emit("")
        elif line.startswith(r"\State"):
            content = clean(line[len(r"\State"):])
            main, comment = split_comment(content)
            main = main.strip()
            if main.startswith(r"\Return"):
                text = f"return {clean(main[len(r'\Return'):])}"
            else:
                text = main
            if comment:
                emit(f"{text}    # {comment}")
            else:
                emit(text)
        elif line.startswith(r"\Comment"):
            emit(f"# {clean(line[len(r'\Comment'):])}")
        elif line.startswith(r"\Return"):
            emit(f"return {clean(line[len(r'\Return'):])}")
        else:
            # Unrecognized line; pass through verbatim rather than drop it.
            emit(clean(line))
    return "\n".join(out)


def main() -> None:
    if len(sys.argv) != 2:
        print("usage: convert_algorithms.py <input.tex>", file=sys.stderr)
        sys.exit(64)
    with open(sys.argv[1]) as f:
        tex = f.read()
    for alg in extract_algorithms(tex):
        body = render_body(alg["body"])
        print(f"::: {{#{alg['label']} .algorithm}}")
        print(f"**Algorithm. {alg['caption']}.**")
        print()
        print("```")
        print(body)
        print("```")
        print(":::")
        print()


if __name__ == "__main__":
    main()
