# linear-algebra skill

A Pi skill for implementing and validating numerical linear algebra algorithms in Julia using *Matrix Computations* by Golub and Van Loan as the primary reference.

## Files

- `SKILL.md` — skill frontmatter and agent instructions
- `scripts/book_lookup.py` — lookup helper for topics and internal references
- `references/sections/` — primary split references
- `references/chapters/` — chapter-level references
- `references/Matrix_Computations_merged.md` — full merged source
- `references/topic-map.md` — quick chapter/topic map
- `references/errata.md` — known book errata tracked locally

## Recommended usage

Inside Pi, either let the skill auto-trigger or force it with:

```bash
/skill:linear-algebra
```

Typical lookup flow:

```bash
python3 scripts/book_lookup.py heading "conjugate gradient"
python3 scripts/book_lookup.py section 11.3
python3 scripts/book_lookup.py theorem 2.4.1
python3 scripts/book_lookup.py algorithm 5.1.1
python3 scripts/book_lookup.py ref '(5.1.2)'
```

The helper prints the best matching file under `references/sections/`; open that file first.

## Retrieval model

- Topic search: `heading`, `text`, `toc`
- Section references: `section`, `ref '§6.5'`
- Labeled results: `theorem`, `lemma`, `corollary`, `algorithm`, `equation`
- Generic references: `ref 'Theorem 2.4.1'`, `ref '(5.1.2)'`
- Line mapping: `file-at`

## Design intent

- Use the book for algorithm choice, derivation, and stability guidance.
- Use Julia best practices separately via the `julia` skill when needed.
- Read the smallest relevant section first, then widen to chapter context only if necessary.

## Notes

- The source markdown is fixed and copied into the skill, so the split files do not need regeneration.
- The markdown transcription is expected to be reliable. If something looks wrong, check `references/errata.md` first, then consider a book erratum or notation issue before assuming extraction error.
