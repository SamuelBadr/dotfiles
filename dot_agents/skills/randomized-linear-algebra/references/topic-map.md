# Randomized NLA topic map

Quick orientation for Martinsson & Tropp, *Randomized Numerical Linear Algebra* (arXiv:2002.01387). Use before opening the full merged source.

- §1 Introduction — the case for randomness; what randomization accomplishes; design considerations
- §2 Linear algebra preliminaries — norms, SVD, pseudoinverse, projectors, spectral functions, restricted singular values
- §3 Probability preliminaries — concentration, matrix Gaussian series, subgaussian/subexponential norms, decoupling
- §4 Trace estimation by sampling — Hutchinson's estimator, randomized power/Krylov methods, bootstrapping
- §5 Schatten *p*-norm estimation — sampling methods for Schatten norms beyond trace
- §6 Maximum eigenvalues and trace functions — randomized power/Krylov for extremal eigenvalues, Lanczos for trace functions
- §7 Matrix approximation by sampling — column/row sampling, leverage scores, uniform sampling bounds
- §8 Randomized embeddings — Gaussian/Subsampled Randomized Transforms (SRTs), embeddings for least squares
- §9 Structured random embeddings — SRT/Hadamard-based maps, coherence reduction
- §10 How to use random embeddings — overdetermined least squares, sketch-and-solve, choosing embedding dimension
- §11 The randomized rangefinder — randomized range finder, power iteration, Krylov rangefinder; fixed-precision approximation
- §12 Error estimation and adaptivity — a posteriori error bounds, randomized QB, on-the-fly rank selection
- §13 Finding natural bases — QR, Interpolative Decomposition (ID), CUR factorization from the randomized range
- §14 Nyström approximation — PSD matrix approximation, single-pass Nyström, error analysis
- §15 Single-view algorithms — streaming/one-pass matrix approximation under update constraints
- §16 Factoring full/nearly-full-rank matrices — randomized LU/QR/SVD, column-pivoted QR, rank-revealing factorizations
- §17 General linear solvers — randomized least squares, normal equations sketching, iterative refinement
- §18 Linear solvers for graph Laplacians — sparse Cholesky, support theory, randomized preconditioners for graph Laplacians
- §19 Kernel matrices in machine learning — randomized approximation of kernel/dense PSD matrices, Nystöm for ML kernels
- §20 High-accuracy approximation of kernel matrices — fast multipole / rank-structured approximation, randomized methods for kernel interpolation
