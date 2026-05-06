# Chapter 10

# Large Sparse Eigenvalue Problems

10.1 The Symmetric Lanczos Process   
10.2 Lanczos, Quadrature, and Approximation   
10.3 Practical Lanczos Procedures   
10.4 Large Sparse SVD Frameworks   
10.5 Krylov Methods for Unsymmetric Problems   
10.6 Jacobi-Davidson and Related Methods

The Lanczos process computes a sequence of partial tridiagonalizations that are orthogonally related to a given symmetric matrix A. It is of particular interest if A is large and sparse because, instead of updating A along the way as in the Householder method of §8.2, it simply relies on matrix-vector products. Equally important, information about A’s extremal eigenvalues tends to emerge fairly early during the iteration, making the method very useful in situations where just a few of A’s largest or smallest eigenvalues are desired, together with the corresponding eigenvectors.

The derivation and exact arithmetic attributes of the method are presented in §10.1, including its extraordinary convergence properties. Central to the discussion is the connection to an underlying Krylov subspace that is defined by the starting vector. In §10.2 we point out connections between Gauss quadrature and the Lanczos process that can be used to estimate expressions of the form $u ^ { T } f ( A ) u$ where f (A) is a function of a large, sparse symmetric positive definite matrix A. Unfortunately, a “math book” implementation of the Lanczos method is practically useless because of roundoff error. This makes it necessary to enlist the help of various “workarounds,” which we describe in §10.3. A sparse SVD framework based on Golub-Kahan bidiagonalization is detailed in §10.4. We also introduce the idea of a randomized SVD. The last two sections deal with the more difficult unsymmetric problem. The Arnoldi iteration is a Krylov subspace iteration like Lanczos. To make it effective, it is necessary to extract valuable “restart information” from the Hessenberg matrix sequence that it produces. This is discussed in §10.5 together with a brief presentation of the unsymmetric Lanczos framework. In the last section we derive the Jacobi-Davidson method, which combines Newton ideas with Rayleigh-Ritz refinement.

# Reading Notes

Familiarity with Chapters 5, 7, and 8 is recommended. Within this chapter there are the following dependencies:

$$
\begin{array}{c c c c c c c c} \S 1 0. 1 & \to & \S 1 0. 3 & \to & \S 1 0. 5 & \to & \S 1 0. 6 \\ \downarrow & & \downarrow & & \\ \S 1 0. 2 & & \S 1 0. 4 & & \end{array}
$$

General references for this chapter include Parlett (SEP), Stewart (MAE), Watkins (MEP), Chatelin (EOM), Cullum and Willoughby (LALSE), Meurant (LCG), Saad (NMLE), Kressner (NMSE), and EIG TEMPLATES.
