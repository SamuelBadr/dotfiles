---
name: linear-algebra
description: Guide implementation and validation of numerical linear algebra algorithms using Golub and Van Loan's Matrix Computations as the primary reference. Use when working with matrix factorizations, least squares, conditioning, eigenvalue and SVD methods, sparse Krylov or iterative methods, matrix functions, or related scientific computing code.
---

# linear-algebra skill

Use this skill for scientific-computing tasks where the agent should choose, justify, or implement a numerical linear algebra method with *Matrix Computations* by Golub and Van Loan as the primary source of truth.

## Primary references

- Section-level references: `references/sections/`
- Section manifest: `references/sections/README.md`
- Chapter-level references: `references/chapters/`
- Full merged source: `references/Matrix_Computations_merged.md`
- Quick topic map: `references/topic-map.md`
- Known errata: `references/errata.md`
- Lookup helper: `scripts/book_lookup.py`
- Shared knowledge bundle (preferred when available): `~/.agents/nla-knowledge/`

## Retrieval workflow

1. Do not answer from memory when the book is relevant.
2. Prefer the shared bundle CLI when available: `python3 ~/.agents/nla-knowledge/scripts/nla_lookup.py lookup <reference>` or `search <query> --source gvl`. It provides stable IDs (`gvl:theorem:2.4.1`) and generated indexes across both NLA skills.
3. Use local `scripts/book_lookup.py` as the fallback or for line/file debugging inside this skill.
4. For a labeled item (theorem/algorithm/equation), lookup returns canonical content directly. Open the section file only to widen context.
5. Pull in chapter context only when the local section is not enough.
6. Use the merged file only for exact line-range lookups or debugging a reference.

## Lookup commands

Shared bundle lookup/search (preferred):

```bash
python3 ~/.agents/nla-knowledge/scripts/nla_lookup.py lookup gvl:theorem:2.4.1
python3 ~/.agents/nla-knowledge/scripts/nla_lookup.py lookup 'theorem 2.4.1'
python3 ~/.agents/nla-knowledge/scripts/nla_lookup.py search cholesky positive definite --source gvl
```

Local topic search (fallback / local debugging):

```bash
python3 scripts/book_lookup.py heading "conjugate gradient"
python3 scripts/book_lookup.py text "backward error"
python3 scripts/book_lookup.py toc qr
```

Internal cross-references:

```bash
python3 scripts/book_lookup.py section 6.5
python3 scripts/book_lookup.py section 7.7.6
python3 scripts/book_lookup.py ref '§11.3'
```

Labeled results:

```bash
python3 scripts/book_lookup.py theorem 2.4.1
python3 scripts/book_lookup.py lemma 2.6.1
python3 scripts/book_lookup.py corollary 2.4.2
python3 scripts/book_lookup.py algorithm 5.1.1
python3 scripts/book_lookup.py equation 5.1.2
python3 scripts/book_lookup.py ref 'Theorem 2.4.1'
python3 scripts/book_lookup.py ref '(5.1.2)'
```

Line-to-file mapping:

```bash
python3 scripts/book_lookup.py file-at 30619
```

For `text`/`heading`, the helper prints file hints; open that file to widen. For labeled items (`theorem`/`algorithm`/`equation`), it prints the canonical content directly.

## Implementation stance

When writing code:

1. Identify the mathematically appropriate algorithm from the book.
2. Prefer numerically stable formulations discussed in the text over naive formulas.
3. Then implement the algorithm idiomatically in the target language.
4. State important assumptions clearly: dimensions, rank conditions, symmetry, definiteness, conditioning, stopping criteria, and error model.
5. When useful, cite the chapter/section or theorem that motivates the implementation.

For randomized NLA methods (randomized SVD, trace estimation, sketching, Nyström approximation), see the `randomized-linear-algebra` skill.

## Common entry points

- LU / Gaussian elimination / pivoting → Chapter 3
- SPD / banded / Toeplitz / Poisson systems → Chapter 4
- Householder / Givens / QR / least squares → Chapter 5
- Regularization / constrained least squares / total least squares → Chapter 6
- Unsymmetric eigenvalue problems / Hessenberg / QR / Schur / pseudospectra → Chapter 7
- Symmetric eigenvalue problems / Jacobi / SVD computation → Chapter 8
- Matrix functions / exponential / sign / square root / logarithm → Chapter 9
- Lanczos / sparse SVD / Arnoldi / Jacobi-Davidson → Chapter 10
- Conjugate gradient / GMRES / MINRES / LSQR / preconditioning / multigrid → Chapter 11
- Kronecker products / tensor contractions → Chapter 12

## Notes

- Prefer nearby definitions, theorems, examples, and algorithms over isolated search hits.
- The markdown transcription is expected to be reliable. If a formula or algorithm seems inconsistent, first check `references/errata.md`, then consider a book erratum or notation issue before blaming extraction.
- If the book does not fully answer the question, say so explicitly before using general knowledge.
