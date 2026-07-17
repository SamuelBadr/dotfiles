# Section split manifest

Source: Martinsson & Tropp, *Randomized Numerical Linear Algebra* (arXiv:2002.01387v3).
Converted via `scripts/convert_latex.sh` (pandoc 3.10 + algorithm extractor)
and split by `scripts/split_sections.py`. Line ranges refer to `../MT_merged.md`.
The `Kind` column holds the section label (`sec:...`) used by `book_lookup.py section`.

| File | Lines | Kind |
| --- | --- | --- |
| `00-frontmatter.md` | 1-26 | frontmatter |
| `01-introduction.md` | 27-413 | intro |
| `02-linear-algebra-preliminaries.md` | 414-718 | sec:lin-alg |
| `03-probability-preliminaries.md` | 719-820 | sec:prob |
| `04-trace-estimation-by-sampling.md` | 821-1386 | sec:trace-est |
| `05-schatten-norm-estimation-by-sampling.md` | 1387-1604 | sec:schatten-p |
| `06-maximum-eigenvalues-and-trace-functions.md` | 1605-2202 | sec:max-eig |
| `07-matrix-approximation-by-sampling.md` | 2203-2656 | sec:matrix-mc |
| `08-randomized-embeddings.md` | 2657-3126 | sec:gauss |
| `09-structured-random-embeddings.md` | 3127-3534 | sec:dimension-reduction |
| `10-how-to-use-random-embeddings.md` | 3535-3807 | sec:overdet-ls |
| `11-the-randomized-rangefinder.md` | 3808-4847 | sec:random-rangefinder |
| `12-error-estimation-and-adaptivity.md` | 4848-5255 | sec:error-est |
| `13-finding-natural-bases-qr-id-and-cur.md` | 5256-5786 | sec:natural |
| `14-nystrom-approximation.md` | 5787-6024 | sec:nystrom |
| `15-single-view-algorithms.md` | 6025-6418 | sec:singlepass |
| `16-factoring-matrices-of-full-or-nearly-full-rank.md` | 6419-7138 | sec:full |
| `17-general-linear-solvers.md` | 7139-7461 | sec:linear-solve |
| `18-linear-solvers-for-graph-laplacians.md` | 7462-7968 | sec:sparse-cholesky |
| `19-kernel-matrices-in-machine-learning.md` | 7969-8700 | sec:kernel |
| `20-high-accuracy-approximation-of-kernel-matrices.md` | 8701-9763 | sec:rankstructured |
