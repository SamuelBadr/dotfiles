# Chapter 4

# Special Linear Systems

4.1 Diagonal Dominance and Symmetry   
4.2 Positive Definite Systems   
4.3 Banded Systems   
4.4 Symmetric Indefinite Systems   
4.5 Block Tridiagonal Systems   
4.6 Vandermonde Systems   
4.7 Classical Methods for Toeplitz Systems   
4.8 Circulant and Discrete Poisson Systems

It is a basic tenet of numerical analysis that solution procedures should exploit structure whenever it is present. In numerical linear algebra, this translates into an expectation that algorithms for general linear systems can be streamlined in the presence of such properties as symmetry, definiteness, and bandedness. Two themes prevail:

• There are important classes of matrices for which it is safe not to pivot when computing the LU or a related factorization.

• There are important classes of matrices with highly structured LU factorizations that can be computed quickly, sometimes, very quickly.

Challenges arise when a fast, but unstable, LU factorization is available.

Symmetry and diagonal dominance are prime examples of exploitable matrix structure and we use these properties to introduce some key ideas in §4.1. In §4.2 we examine the case when A is both symmetric and positive definite, deriving the stable Cholesky factorization. Unsymmetric positive definite systems are also investigated. In §4.3, banded versions of the LU and Cholesky factorizations are discussed and this is followed in §4.4 with a treatment of the symmetric indefinite problem. Block matrix ideas and sparse matrix ideas come together when the matrix of coefficients is block tridiagonal. This important class of systems receives a special treatment in §4.5.

Classical methods for Vandermonde and Toeplitz systems are considered in §4.6 and §4.7. In §4.8 we connect the fast transform discussion in §1.4 to the problem of solving circulant systems and systems that arise when the Poisson problem is discretized using finite differences.

Before we get started, we clarify some terminology associated with structured problems that pertains to this chapter and beyond. Banded matrices and block-banded matrices are examples of sparse matrices, meaning that the vast majority of their entries are zero. Linear equation methods that are appropriate when the zero-nonzero pattern is more arbitrary are discussed in Chapter 11. Toeplitz, Vandermonde, and circulant matrices are data sparse. A matrix $A \in \mathbb { R } ^ { m \times n }$ is data sparse if it can be parameterized with many fewer than O(mn) numbers. Cauchy-like systems and semiseparable systems are considered in §12.1 and §12.2.

# Reading Notes

Knowledge of Chapters 1, 2, and 3 is assumed. Within this chapter there are the following dependencies:

$$
\begin{array}{c c c c c c c c c c} \S 4. 1 & \to & \S 4. 2 & \to & \S 4. 3 & \to & \S 4. 4 \\ \downarrow & & & & \downarrow & & \\ \S 4. 6 & & & & \S 4. 5 & \to & \S 4. 7 & \to & \S 4. 8 \end{array}
$$

Global references include Stewart( MABD), Higham (ASNA), Watkins (FMC), Trefethen and Bau (NLA), Demmel (ANLA), and Ipsen (NMA).
