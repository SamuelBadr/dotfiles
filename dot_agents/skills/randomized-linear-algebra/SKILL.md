---
name: randomized-linear-algebra
description: Guide implementation and validation of randomized numerical linear algebra algorithms using Martinsson & Tropp's "Randomized Numerical Linear Algebra" survey as the primary reference. Use when working with randomized SVD or eigensolvers, the randomized rangefinder, trace or Schatten-norm estimation, random embeddings / sketching, Nyström approximation, single-view matrix approximation, or sketch-and-solve linear systems.
---

# randomized-linear-algebra skill

Use this skill for scientific-computing tasks where the agent should choose, justify, or implement a randomized numerical linear algebra method with *Randomized Numerical Linear Algebra: Foundations & Algorithms* by Martinsson & Tropp (arXiv:2002.01387) as the primary source of truth.

## Primary references

- Section-level references: `references/sections/`
- Section manifest: `references/sections/README.md`
- Full merged source: `references/MT_merged.md`
- Quick topic map: `references/topic-map.md`
- Known errata: `references/errata.md`
- Lookup helper: `scripts/book_lookup.py`
- Shared knowledge bundle (preferred when available): `~/.agents/nla-knowledge/`
- Conversion pipeline: `scripts/convert_latex.sh` + `scripts/convert_algorithms.py` + `scripts/split_sections.py`

## Scope boundary

This skill covers **randomized** matrix computations: methods that use random sampling, sketching, or probability to approximate matrix operations.

For **deterministic** numerical linear algebra (LU, QR, deterministic SVD, conditioning, perturbation theory, Krylov methods without random sketching), use the `linear-algebra` skill instead.

Randomized Krylov methods (where the randomization is load-bearing to the analysis) belong here. A plain deterministic Lanczos implementation belongs to `linear-algebra`.

## Retrieval workflow

1. Do not answer from memory when the survey is relevant.
2. Prefer the shared bundle CLI when available: `python3 ~/.agents/nla-knowledge/scripts/nla_lookup.py lookup <reference>` or `search <query> --source mt`. It provides stable IDs (`mt:algorithm:slq`) and generated indexes across both NLA skills.
3. Use local `scripts/book_lookup.py` as the fallback or for line/file debugging inside this skill.
4. For a labeled item (theorem/algorithm/equation), lookup returns canonical content directly. Open the section file only to widen context.
5. Pull in broader section context only when the canonical item is not enough.
6. Use the merged file only for exact line-range lookups or debugging a reference.

## Lookup commands

Shared bundle lookup/search (preferred):

```bash
python3 ~/.agents/nla-knowledge/scripts/nla_lookup.py lookup mt:algorithm:slq
python3 ~/.agents/nla-knowledge/scripts/nla_lookup.py lookup 'algorithm trace-est'
python3 ~/.agents/nla-knowledge/scripts/nla_lookup.py lookup 'theorem 8'
python3 ~/.agents/nla-knowledge/scripts/nla_lookup.py search stochastic lanczos log determinant --source mt
```

Local topic search (fallback / local debugging):

```bash
python3 scripts/book_lookup.py heading "randomized rangefinder"
python3 scripts/book_lookup.py text "trace estimation"
python3 scripts/book_lookup.py toc
```

Section resolution (by ordinal or label):

```bash
python3 scripts/book_lookup.py section 4
python3 scripts/book_lookup.py section sec:trace-est
```

Labeled items — by sequential number or named label:

```bash
python3 scripts/book_lookup.py theorem 8
python3 scripts/book_lookup.py theorem rand-power
python3 scripts/book_lookup.py algorithm trace-est
python3 scripts/book_lookup.py equation rand-power-nogap
python3 scripts/book_lookup.py ref 'theorem 8'
python3 scripts/book_lookup.py ref 'theorem rand-power'
```

Line-to-file mapping:

```bash
python3 scripts/book_lookup.py file-at 1742
```

For `text`/`heading`, the helper prints file hints; open that file to widen. For labeled items (`theorem`/`algorithm`/`equation`), it prints the canonical content directly.

## Addressing scheme

The survey numbers theorems/lemmas/corollaries sequentially (one shared counter) and cross-references them by name ("Theorem 8") and by named label (`thm:rand-power`). Algorithms and equations are referenced by named label only (`alg:trace-est`, `eqn:rand-power-nogap`); they are not sequentially numbered. The lookup accepts either form: `theorem 8` (sequential) or `theorem rand-power` (label).

## Implementation stance

When writing code:

1. Identify the mathematically appropriate algorithm from the survey.
2. Prefer numerically stable formulations discussed in the text over naive formulas.
3. Then implement the algorithm idiomatically in the target language.
4. State important assumptions clearly: dimensions, rank conditions, symmetry, definiteness, conditioning, stopping criteria, and error model.
5. When useful, cite the section or theorem that motivates the implementation.

## Common entry points

- Trace / Schatten-norm estimation → §4–5
- Maximum eigenvalues, trace functions, randomized power/Krylov → §6
- Matrix approximation by sampling (leverage scores, column sampling) → §7
- Random / structured random embeddings (Gaussian, SRT, Hadamard) → §8–9
- Sketch-and-solve, overdetermined least squares → §10
- Randomized rangefinder, fixed-precision approximation → §11
- Error estimation and adaptive rank selection → §12
- Interpolative decomposition, QR, CUR factorizations → §13
- Nyström approximation (PSD matrices) → §14
- Single-view / streaming algorithms → §15
- Randomized full-rank factorizations (LU/QR/SVD) → §16
- Sketch-and-solve linear systems → §17
- Graph Laplacian solvers, randomized preconditioners → §18
- Kernel matrices (ML), coordinate Nyström → §19
- High-accuracy kernel approximation, rank-structured methods → §20

## Notes

- Prefer nearby definitions, theorems, examples, and algorithms over isolated search hits.
- The markdown is converted from the LaTeX source via `scripts/convert_latex.sh`. Algorithm pseudocode is extracted from the LaTeX `algorithmic` environments by `scripts/convert_algorithms.py` (pandoc drops them entirely). Math macros (`\mtx{A}`, `\vct{\omega}`, `\norm{...}`) pass through unexpanded inside code fences as readable shorthand (matrix/vector); pure-formatting macros (`\texttt{randn}`→`randn`, `\textbf{break}`→`break`) and `\_`→`_` are unwrapped.
- If a formula or algorithm seems inconsistent, first check `references/errata.md`, then consider a source issue or notation discrepancy before blaming the conversion.
- If the survey does not fully answer the question, say so explicitly before using general knowledge.
- For deterministic NLA methods, see the `linear-algebra` skill.
