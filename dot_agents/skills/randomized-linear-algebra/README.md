# randomized-linear-algebra skill

A Pi skill for implementing and validating randomized numerical linear algebra algorithms using *Randomized Numerical Linear Algebra: Foundations & Algorithms* by Martinsson and Tropp (arXiv:2002.01387) as the primary reference.

## What it does

Provides a retrieval layer over the Martinsson & Tropp survey so the agent grounds randomized-NLA implementations in the canonical source rather than memory. Covers trace estimation, randomized SVD/eigensolvers, the randomized rangefinder, random embeddings and sketching, Nyström approximation, single-view algorithms, and sketch-and-solve linear systems.

## Structure

- `SKILL.md` — trigger, retrieval workflow, lookup commands, scope boundary
- `scripts/book_lookup.py` — local BM25 + Snowball-stemmed lookup; resolves theorems/algorithms/equations by sequential number or named label
- `~/.agents/nla-knowledge/` — shared OKF-inspired knowledge bundle and `scripts/nla_lookup.py` CLI (preferred cross-skill lookup layer)
- `scripts/convert_latex.sh` + `scripts/convert_algorithms.py` + `scripts/split_sections.py` — reproducible LaTeX→markdown conversion (pandoc + `\rm` preprocessor + algorithm extractor + section split)
- `references/MT_merged.md` — full converted survey (source of truth)
- `references/sections/` — 20 section files + manifest
- `references/topic-map.md` — quick orientation
- `references/errata.md` — known corrections (starts empty)
- `CONTEXT.md` — domain model (ubiquitous language)
- `docs/adr/` — architectural decisions

## Relationship to linear-algebra skill

This skill covers **randomized** matrix computations. The `linear-algebra` skill covers **deterministic** methods (LU, QR, deterministic SVD, Krylov methods without random sketching) via Golub & Van Loan. The two cross-link bidirectionally.

## Dependencies

- `pandoc` 3.10 (for conversion; not needed for lookup)
- `pystemmer` (optional; enables Snowball stemming for text/heading search; falls back to unstemmed)
