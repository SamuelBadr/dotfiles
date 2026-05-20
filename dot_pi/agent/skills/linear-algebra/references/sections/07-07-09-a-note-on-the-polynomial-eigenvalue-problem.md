# 7.7.9 A Note on the Polynomial Eigenvalue Problem

More general than the generalized eigenvalue problem is the polynomial eigenvalue problem. Here we are given matrices $A _ { 0 } , \dots , A _ { d } \in \mathbb { C } ^ { n \times n }$ and determine $\lambda \in \mathbb { C }$ and $0 \neq x \in \mathbb { C } ^ { n }$ so that

$$
P (\lambda) x = 0 \tag {7.7.3}
$$

where the λ-matrix P (λ) is defined by

$$
P (\lambda) = A _ {0} + \lambda A _ {1} + \dots + \lambda^ {d} A _ {d}. \tag {7.7.4}
$$

We assume $A _ { d } \neq 0$ and regard d as the degree of $P ( \lambda )$ . The theory behind the polynomial eigenvalue problem is nicely developed in Lancaster (1966).

It is possible to convert (7.7.3) into an equivalent linear eigenvalue problem with larger dimension. For example, suppose d = 3 and

$$
L (\lambda) = \left[ \begin{array}{c c c} 0 & 0 & A _ {0} \\ - I & 0 & A _ {1} \\ 0 & - I & A _ {2} \end{array} \right] + \lambda \left[ \begin{array}{c c c} I & 0 & 0 \\ 0 & I & 0 \\ 0 & 0 & A _ {3} \end{array} \right]. \tag {7.7.5}
$$

If

$$
L (\lambda) \left[ \begin{array}{l} u _ {1} \\ u _ {2} \\ x \end{array} \right] = \left[ \begin{array}{l} 0 \\ 0 \\ 0 \end{array} \right],
$$

then

$$
0 = A _ {0} x + \lambda u _ {1} = A _ {0} + \lambda \left(A _ {1} x + \lambda u _ {2}\right) = A _ {0} + \lambda \left(A _ {1} x + \lambda \left(A _ {2} + \lambda A _ {3}\right)\right) x = P (\lambda) x.
$$

In general, we say that L(λ) is a linearization of P (λ) if there are dn-by-dn λ-matrices S(λ) and T (λ), each with constant nonzero determinants, so that

$$
S (\lambda) \left[ \begin{array}{c c} P (\lambda) & 0 \\ 0 & I _ {(d - 1) n} \end{array} \right] T (\lambda) = L (\lambda) \tag {7.7.6}
$$

has unit degree. With this conversion, the A − λB methods just discussed can be applied to find the required eigenvalues and eigenvectors.

Recent work has focused on how to choose the λ-transformations S(λ) and $T ( \lambda )$ so that special structure in P (λ) is reflected in L(λ). See Mackey, Mackey, Mehl, and Mehrmann (2006). The idea is to think of (7.7.6) as a factorization and to identify the transformations that produce a properly structured L(λ). To appreciate this solution framework it is necessary to have a facility with λ-matrix manipulation and to that end we briefly examine the λ-matrix transformations behind the above linearization. If

$$
P _ {1} (\lambda) = A _ {1} + \lambda A _ {2} + \dots + \lambda^ {d - 1} A _ {d}
$$

then

$$
P (\lambda) = A _ {0} + \lambda P _ {1} (\lambda)
$$

and it is easy to verify that

$$
\left[ \begin{array}{c c} I _ {n} & - \lambda I _ {n} \\ 0 & I _ {n} \end{array} \right] \left[ \begin{array}{c c} A _ {0} + \lambda P _ {1} (\lambda) & 0 \\ 0 & I _ {n} \end{array} \right] \left[ \begin{array}{c c} 0 & I _ {n} \\ - I _ {n} & P _ {1} (\lambda) \end{array} \right] = \left[ \begin{array}{c c} \lambda I _ {n} & A _ {0} \\ - I _ {n} & P _ {1} (\lambda) \end{array} \right].
$$

Notice that the transformation matrices have unit determinant and that the λ-matrix on the right-hand side has degree d − 1. The process can be repeated. If

$$
P _ {2} (\lambda) = A _ {2} + \lambda A _ {3} + \dots + \lambda^ {d - 2} A _ {d}
$$

then

$$
P _ {1} (\lambda) = A _ {1} + \lambda P _ {2} (\lambda)
$$

and

$$
\left[ \begin{array}{c c c} I _ {n} & 0 & 0 \\ \hline 0 & I _ {n} & - \lambda I _ {n} \\ 0 & 0 & I _ {n} \end{array} \right] \left[ \begin{array}{c c c} \lambda I _ {n} & A _ {0} & 0 \\ \hline - I _ {n} & P _ {1} (\lambda) & 0 \\ 0 & 0 & I _ {n} \end{array} \right] \left[ \begin{array}{c c c} I _ {n} & 0 & 0 \\ \hline 0 & 0 & I _ {n} \\ 0 & - I _ {n} & P _ {2} (\lambda) \end{array} \right] =
$$

$$
\left[ \begin{array}{c c c} \lambda I _ {n} & 0 & A _ {0} \\ \hline - I _ {n} & \lambda I _ {n} & A _ {1} \\ 0 & - I _ {n} & P _ {2} (\lambda) \end{array} \right].
$$

Note that the matrix on the right has degree d − 2. A straightforward induction argument can be assembled to establish that if the dn-by-dn matrices $S ( \lambda )$ and $T ( \lambda )$ are defined by

$$
S (\lambda) = \left[ \begin{array}{c c c c c} I _ {n} & - \lambda I _ {n} & 0 & \dots & 0 \\ 0 & I _ {n} & - \lambda I _ {n} & & \vdots \\ 0 & & \ddots & \ddots & \\ \vdots & & & I _ {n} & - \lambda I _ {n} \\ 0 & 0 & \dots & 0 & I _ {n} \end{array} \right], \quad T (\lambda) = \left[ \begin{array}{c c c c c} 0 & 0 & 0 & \dots & I \\ - I _ {n} & 0 & & & P _ {1} (\lambda) \\ 0 & - I _ {n} & \ddots & & \vdots \\ \vdots & & \ddots & \ddots & P _ {d - 2} (\lambda) \\ 0 & 0 & \dots & - I _ {n} & P _ {d - 1} (\lambda) \end{array} \right]
$$

where

$$
P _ {k} (\lambda) = A _ {k} + \lambda A _ {k + 1} + \dots + \lambda^ {d - k} A _ {d},
$$

then

$$
S (\lambda) \left[ \begin{array}{c c} P (\lambda) & 0 \\ 0 & I _ {(d - 1) n} \end{array} \right] T (\lambda) = \left[ \begin{array}{c c c c c} \lambda I _ {n} & 0 & 0 & \dots & A _ {0} \\ - I _ {n} & \lambda I _ {n} & & & A _ {1} \\ 0 & - I _ {n} & \ddots & & \vdots \\ \vdots & & \ddots & \lambda I _ {n} & A _ {d - 2} \\ 0 & 0 & \dots & - I _ {n} & A _ {d - 1} + \lambda A _ {d} \end{array} \right].
$$

Note that, if we solve the linearized problem using the QZ algorithm, then $O ( ( d n ) ^ { 3 } )$ flops are required.

# Problems

P7.7.1 Suppose A and B are in $\mathbb { R } ^ { n \times n }$ and that

$$
U ^ {T} B V = \left[ \begin{array}{c c} D & 0 \\ 0 & 0 \\ r & n - r \end{array} \right] _ {n - r} ^ {r}, \quad U = \left[ \begin{array}{c c} U _ {1} & U _ {2} \\ r & n - r \end{array} \right], \quad V = \left[ \begin{array}{c c} V _ {1} & V _ {2} \\ r & n - r \end{array} \right],
$$

is the SVD of B, where D is r-by-r and $r = { \mathsf { r a n k } } ( B )$ . Show that if $\lambda ( A , B ) = \mathbb { C }$ then $U _ { 2 } ^ { T } A V _ { 2 }$ is singular.

P7.7.2 Suppose A and B are in $\mathbb { R } ^ { n \times n }$ . Give an algorithm for computing orthogonal Q and Z such that $Q ^ { T } A \dot { Z }$ is upper Hessenberg and $Z ^ { T } B Q$ is upper triangular.

# P7.7.3 Suppose

$$
A = \left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ 0 & A _ {2 2} \end{array} \right] \quad \text {and} \quad B = \left[ \begin{array}{c c} B _ {1 1} & B _ {1 2} \\ 0 & B _ {2 2} \end{array} \right]
$$

with $A _ { 1 1 } , B _ { 1 1 } \in \mathbb { R } ^ { k \times k }$ and $A _ { 2 2 } , B _ { 2 2 } \in \mathbb { R } ^ { j \times j }$ . Under what circumstances do there exist

$$
X = \left[ \begin{array}{c c} I _ {k} & X _ {1 2} \\ 0 & I _ {j} \end{array} \right] \quad \text {and} \quad Y = \left[ \begin{array}{c c} I _ {k} & Y _ {1 2} \\ 0 & I _ {j} \end{array} \right]
$$

so that $Y ^ { - 1 } A X$ and $Y ^ { - 1 } B X$ are both block diagonal? This is the generalized Sylvester equation problem. Specify an algorithm for the case when $A _ { 1 1 } , A _ { 2 2 } , B _ { 1 1 }$ , and $B _ { 2 2 }$ are upper triangular. See K˚agstr¨om (1994).

P7.7.4 Suppose $\mu \not \in \lambda ( A , B )$ . Relate the eigenvalues and eigenvectors of $A _ { 1 } = ( A - \mu B ) ^ { - 1 } A$ and $B _ { 1 } = ( A - \mu B ) ^ { - 1 } B$ to the generalized eigenvalues and eigenvectors of $A - \lambda B$ .

P7.7.5 What does the generalized Schur decomposition say about the pencil $A - \lambda A ^ { T } ?$ Hint: If $T \in \mathbb { R } ^ { n \times n }$ is upper triangular, then $\mathcal { E } _ { n } T \mathcal { E } _ { n }$ is lower triangular where $\mathcal { E } _ { n }$ is the exchange permutation defined in §1.2.11.

# P7.7.6 Prove that

$$
L _ {1} (\lambda) = \left[ \begin{array}{c c c c} A _ {3} + \lambda A _ {4} & A _ {2} & A _ {1} & A _ {0} \\ - I _ {n} & 0 & 0 & 0 \\ 0 & - I _ {n} & 0 & 0 \\ 0 & 0 & - I _ {n} & 0 \end{array} \right], L _ {2} (\lambda) = \left[ \begin{array}{c c c c} A _ {3} + \lambda A _ {4} & - I _ {n} & 0 & 0 \\ A _ {2} & 0 & - I _ {n} & 0 \\ A _ {1} & 0 & 0 & - I _ {n} \\ A _ {0} & 0 & 0 & 0 \end{array} \right]
$$

are linearizations of

$$
P (\lambda) = A _ {0} + \lambda A _ {1} + \lambda^ {2} A _ {2} + \lambda^ {3} A _ {3} + \lambda^ {4} A _ {4}.
$$

Specify the λ-matrix transformations that relate diag $( P ( \lambda ) , I _ { 3 n } )$ to both $L _ { 1 } ( \lambda )$ and $L _ { 2 } ( \lambda )$ .

# Notes and References for §7.7

For background to the generalized eigenvalue problem we recommend Stewart(IMC), Stewart and Sun (MPT), and Watkins (MEP) and:

B. K˚agstr¨om and A. Ruhe (1983). Matrix Pencils, Proceedings Pite Havsbad, 1982, Lecture Notes in Mathematics Vol. 973, Springer-Verlag, New York.

QZ-related papers include:

C.B. Moler and G.W. Stewart (1973). “An Algorithm for Generalized Matrix Eigenvalue Problems,” SIAM J. Numer. Anal. 10, 241–256.

L. Kaufman (1974). “The LZ Algorithm to Solve the Generalized Eigenvalue Problem,” SIAM J. Numer. Anal. 11, 997–1024.

R.C. Ward (1975). “The Combination Shift QZ Algorithm,” SIAM J. Numer. Anal. 12, 835–853.

C.F. Van Loan (1975). “A General Matrix Eigenvalue Algorithm,” SIAM J. Numer. Anal. 12, 819–834.

L. Kaufman (1977). “Some Thoughts on the QZ Algorithm for Solving the Generalized Eigenvalue Problem,” ACM Trans. Math. Softw. 3, 65–75.

R.C. Ward (1981). “Balancing the Generalized Eigenvalue Problem,” SIAM J. Sci. Stat. Comput. 2, 141–152.

P. Van Dooren (1982). “Algorithm 590: DSUBSP and EXCHQZ: Fortran Routines for Computing Deflating Subspaces with Specified Spectrum,” ACM Trans. Math. Softw. 8, 376–382.

K. Dackland and B. K˚agstr¨om (1999). “Blocked Algorithms and Software for Reduction of a Regular Matrix Pair to Generalized Schur Form,” ACM Trans. Math. Softw. 25, 425–454.

D.S. Watkins (2000). “Performance of the QZ Algorithm in the Presence of Infinite Eigenvalues,” SIAM J. Matrix Anal. Applic. 22, 364–375.

B. K˚agstr¨om, D. Kressner, E.S. Quintana-Orti, and G. Quintana-Orti (2008). “Blocked Algorithms for the Reduction to Hessenberg-Triangular Form Revisited,” BIT 48, 563–584.

Many algorithmic ideas associated with the A − λI problem extend to the A − λB problem:

A. Jennings and M.R. Osborne (1977). “Generalized Eigenvalue Problems for Certain Unsymmetric Band Matrices,” Lin. Alg. Applic. 29, 139–150.

V.N. Kublanovskaya (1984). “AB Algorithm and Its Modifications for the Spectral Problem of Linear Pencils of Matrices,” Numer. Math. 43, 329–342.   
Z. Bai, J. Demmel, and M. Gu (1997). “An Inverse Free Parallel Spectral Divide and Conquer Algorithm for Nonsymmetric Eigenproblems,” Numer. Math. 76, 279–308.   
G.H. Golub and Q. Ye (2000). “Inexact Inverse Iteration for Generalized Eigenvalue Problems,” BIT 40, 671–684.   
F. Tisseur (2001). “Newton’s Method in Floating Point Arithmetic and Iterative Refinement of Generalized Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 22, 1038–1057.   
D. Lemonnier and P. Van Dooren (2006). “Balancing Regular Matrix Pencils,” SIAM J. Matrix Anal. Applic. 28, 253–263.   
R. Granat, B. K˚agstr¨om, and D. Kressner (2007). “Computing Periodic Deflating Subspaces Associated with a Specified Set of Eigenvalues,” BIT 47, 763–791.   
The perturbation theory for the generalized eigenvalue problem is treated in:   
G.W. Stewart (1972). “On the Sensitivity of the Eigenvalue Problem Ax = λBx,” SIAM J. Numer. Anal. 9, 669–686.   
G.W. Stewart (1973). “Error and Perturbation Bounds for Subspaces Associated with Certain Eigenvalue Problems,” SIAM Review 15, 727–764.   
G.W. Stewart (1975). “Gershgorin Theory for the Generalized Eigenvalue Problem Ax = λBx,” Math. Comput. 29, 600–606.   
A. Pokrzywa (1986). “On Perturbations and the Equivalence Orbit of a Matrix Pencil,” Lin. Alg. Applic. 82, 99–121.   
J. Sun (1995). “Perturbation Bounds for the Generalized Schur Decomposition,” SIAM J. Matrix Anal. Applic. 16, 1328–1340.   
R. Bhatia and R.–C. Li (1996). “On Perturbations of Matrix Pencils with Real Spectra. II,” Math. Comput. 65, 637–645.   
J.-P. Dedieu (1997). “Condition Operators, Condition Numbers, and Condition Number Theorem for the Generalized Eigenvalue Problem,” Lin. Alg. Applic. 263, 1–24.   
D.J. Higham and N.J. Higham (1998). “Structured Backward Error and Condition of Generalized Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 20, 493–512.   
R. Byers, C. He, and V. Mehrmann (1998). “Where is the Nearest Non-Regular Pencil?,” Lin. Alg. Applic. 285, 81–105.   
V. Frayss and V. Toumazou (1998). “A Note on the Normwise Perturbation Theory for the Regular Generalized Eigenproblem,” Numer. Lin. Alg. 5, 1–10.   
R.–C. Li (2003). “On Perturbations of Matrix Pencils with Real Spectra, A Revisit,” Math. Comput. 72, 715–728.   
S. Bora and V. Mehrmann (2006). “Linear Perturbation Theory for Structured Matrix Pencils Arising in Control Theory,” SIAM J. Matrix Anal. Applic. 28, 148–169.   
X.S. Chen (2007). “On Perturbation Bounds of Generalized Eigenvalues for Diagonalizable Pairs,” Numer. Math. 107, 79–86.   
The Kronecker structure of the pencil A − λB is analogous to Jordan structure of A − λI and it can provide useful information about the underlying application. Papers concerned with this important decomposition include:   
J.H. Wilkinson (1978). “Linear Differential Equations and Kronecker’s Canonical Form,” in Recent Advances in Numerical Analysis, C. de Boor and G.H. Golub (eds.), Academic Press, New York, 231–265.   
J.H. Wilkinson (1979). “Kronecker’s Canonical Form and the QZ Algorithm,” Lin. Alg. Applic. 28, 285–303.   
P. Van Dooren (1979). “The Computation of Kronecker’s Canonical Form of a Singular Pencil,” Lin. Alg. Applic. 27, 103–140.   
J.W. Demmel (1983). “The Condition Number of Equivalence Transformations that Block Diagonalize Matrix Pencils,” SIAM J. Numer. Anal. 20, 599–610.   
J.W. Demmel and B. K˚agstr¨om (1987). “Computing Stable Eigendecompositions of Matrix Pencils,” Linear Alg. Applic. 88/89, 139–186.   
B. K˚agstr¨om (1985). “The Generalized Singular Value Decomposition and the General A − λB Problem,” BIT 24, 568–583.   
B. K˚agstr¨om (1986). “RGSVD: An Algorithm for Computing the Kronecker Structure and Reducing Subspaces of Singular A λB Pencils,” SIAM J. Sci. Stat. Comput. 7, 185–211.

J. Demmel and B. K˚agstr¨om (1986). “Stably Computing the Kronecker Structure and Reducing Subspaces of Singular Pencils A  λB for Uncertain Data,” in Large Scale Eigenvalue Problems, J. Cullum and R.A. Willoughby (eds.), North-Holland, Amsterdam.   
T. Beelen and P. Van Dooren (1988). “An Improved Algorithm for the Computation of Kronecker’s Canonical Form of a Singular Pencil,” Lin. Alg. Applic. 105, 9–65.   
E. Elmroth and B. K˚agstr¨om(1996). “The Set of 2-by-3 Matrix Pencils — Kronecker Structures and Their Transitions under Perturbations,” SIAM J. Matrix Anal. Applic. 17, 1–34.   
A. Edelman, E. Elmroth, and B. K˚agstr¨om (1997). “A Geometric Approach to Perturbation Theory of Matrices and Matrix Pencils Part I: Versal Defformations,” SIAM J. Matrix Anal. Applic. 18, 653–692.   
E. Elmroth, P. Johansson, and B. K˚agstr¨om (2001). “Computation and Presentation of Graphs Displaying Closure Hierarchies of Jordan and Kronecker Structures,” Num. Lin. Alg. 8, 381–399.   
Just as the Schur decomposition can be used to solve the Sylvester equation problem $A _ { 1 } X - X A _ { 2 } = B _ { 1 }$ , the generalized Schur decomposition can be used to solve the generalized Sylvester equation problem where matrices X and Y are sought so that $A _ { 1 } X - Y A _ { 2 } = B _ { 1 }$ and $A _ { 3 } X - Y A _ { 4 } = B _ { 2 }$ , see:   
W. Enright and S. Serbin (1978). “A Note on the Efficient Solution of Matrix Pencil Systems,” BIT 18, 276–81.   
B. K˚agstr¨om and L. Westin (1989). “Generalized Schur Methods with Condition Estimators for Solving the Generalized Sylvester Equation,” IEEE Trans. Autom. Contr. AC-34, 745–751.   
B. K˚agstr¨om (1994). “A Perturbation Analysis of the Generalized Sylvester Equation (AR−LB, DR− LE) = (C, F ),” SIAM J. Matrix Anal. Applic. 15, 1045–1060.   
J.-G. Sun (1996). “Perturbation Analysis of System Hessenberg and Hessenberg-Triangular Forms,” Lin. Alg. Applic. 241–3, 811–849.   
B. K˚agstr¨om and P. Poromaa (1996). “LAPACK-style Algorithms and Software for Solving the Generalized Sylvester Equation and Estimating the Separation Between Regular Matrix Pairs,” ACM Trans. Math. Softw. 22, 78–103.   
I. Jonsson and B. K˚agstr¨om (2002). “Recursive Blocked Algorithms for Solving Triangular Systems– Part II: Two-sided and Generalized Sylvester and Lyapunov Matrix Equations,” ACM Trans. Math. Softw. 28, 416–435.   
R. Granat and B. K˚agstr¨om (2010). “Parallel Solvers for Sylvester-Type Matrix Equations with Applications in Condition Estimation, Part I: Theory and Algorithms,” ACM Trans. Math. Softw. 37, Article 32.   
Rectangular generalized eigenvalue problems also arise. In this setting the goal is to reduce the rank of A − λB, see:   
G.W. Stewart (1994). “Perturbation Theory for Rectangular Matrix Pencils,” Lin. Alg. Applic. 208/209, 297–301.   
G. Boutry, M. Elad, G.H. Golub, and P. Milanfar (2005). “The Generalized Eigenvalue Problem for Nonsquare Pencils Using a Minimal Perturbation Approach,” SIAM J. Matrix Anal. Applic. 27, 582–601.   
D. Chu and G.H. Golub (2006). “On a Generalized Eigenvalue Problem for Nonsquare Pencils,” SIAM J. Matrix Anal. Applic. 28, 770–787.   
References for the polynomial eigenvalue problem include:   
P. Lancaster (1966). Lambda-Matrices and Vibrating Systems, Pergamon Press, Oxford, U.K.   
I. Gohberg, P. Lancaster, and L. Rodman (1982). Matrix Polynomials, Academic Press, New York.   
F. Tisseur (2000). “Backward Error and Condition of Polynomial Eigenvalue Problems,” Lin. Alg. Applic. 309, 339–361.   
J.-P. Dedieu and F. Tisseur (2003). “Perturbation Theory for Homogeneous Polynomial Eigenvalue Problems,” Lin. Alg. Applic. 358, 71–94.   
N.J. Higham, D.S. Mackey, and F. Tisseur (2006). “The Conditioning of Linearizations of Matrix Polynomials,” SIAM J. Matrix Anal. Applic. 28, 1005–1028.   
D.S. Mackey, N. Mackey, C. Mehl, V. Mehrmann (2006). “Vector Spaces of Linearizations for Matrix Polynomials,” SIAM J. Matrix Anal. Applic. 28, 971–1004.   
The structured quadratic eigenvalue problem is discussed briefly in §8.7.9.
