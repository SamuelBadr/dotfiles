# 7.8 Hamiltonian and Product Eigenvalue Problems

Two structured unsymmetric eigenvalue problems are considered. The Hamiltonian matrix eigenvalue problem comes with its own special Schur decomposition. Orthogonal symplectic similarity transformations are used to bring about the required reduction. The product eigenvalue problem involves computing the eigenvalues of a product like $A _ { 1 } A _ { 2 } ^ { - 1 } A _ { 3 }$ without actually forming the product or the designated inverses. For detailed background to these problems, see Kressner (NMGS) and Watkins (MEP).

# 7.8.1 Hamiltonian Matrix Eigenproblems

Hamiltonian and symplectic matrices are introduced in §1.3.10. Their 2-by-2 block structure provide a nice framework for practicing block matrix manipulation, see P1.3.2 and P2.5.4. We now describe some interesting eigenvalue problems that involve these matrices. For a given n, we define the matrix $\check { J } \in \mathbb { R } ^ { 2 n \times 2 n }$ by

$$
J = \left[ \begin{array}{c c} 0 & I _ {n} \\ - I _ {n} & 0 \end{array} \right]
$$

and proceed to work with the families of 2-by-2 block structured matrices that are displayed in Figure 7.8.1. We mention four important facts concerning these matrices.

<table><tr><td>Family</td><td>Definition</td><td colspan="2">What They Look Like</td></tr><tr><td>Hamiltonian</td><td> $JM = (JM)^T$ </td><td> $M = \begin{bmatrix} A & G \\ F & -A^T \end{bmatrix}$ </td><td> $G$  symmetric $F$  symmetric</td></tr><tr><td>SkewHamiltonian</td><td> $JN = -(JN)^T$ </td><td> $N = \begin{bmatrix} A & G \\ F & A^T \end{bmatrix}$ </td><td> $G$  skew-symmetric $F$  skew-symmetric</td></tr><tr><td>Symplectic</td><td> $JS = S^{-T}J$ </td><td> $S = \begin{bmatrix} S_{11} & S_{12} \\ S_{21} & S_{22} \end{bmatrix}$ </td><td> $S_{11}^T S_{21}$  symmetric $S_{22}^T S_{12}$  symmetric $S_{11}^T S_{22} = I + S_{21}^T S_{12}$ </td></tr><tr><td>OrthogonalSymplectic</td><td> $JQ = QJ$ </td><td> $Q = \begin{bmatrix} Q_1 & Q_2 \\ -Q_2 & Q_1 \end{bmatrix}$ </td><td> $Q_1^T Q_2$  symmetric $I = Q_1^T Q_1 + Q_2^T Q_2$ </td></tr></table>

Figure 7.8.1. Hamiltonian and symplectic structures

(1) Symplectic similarity transformations preserve Hamiltonian structure:

$$
J (S ^ {- 1} M S) = (J S ^ {- 1} J ^ {T}) (J M J ^ {T}) (J S) = - S ^ {T} M ^ {T} S ^ {- T} J = (J (S ^ {- 1} M S)) ^ {T}.
$$

(2) The square of a Hamiltonian matrix is skew-Hamiltonian:

$$
J M ^ {2} = (J M J ^ {T}) (J M) = - M ^ {T} (J M) ^ {T} = - M ^ {2 T} J ^ {T} = - (J M ^ {2}) ^ {T}.
$$

(3) If M is a Hamiltonian matrix and $\lambda \in \lambda ( M )$ , then $- \lambda \in \lambda ( M )$ :

$$
M \left[ \begin{array}{c} u \\ v \end{array} \right] = \lambda \left[ \begin{array}{c} u \\ v \end{array} \right] \qquad \Rightarrow \qquad M ^ {T} \left[ \begin{array}{c} v \\ - u \end{array} \right] = - \lambda \left[ \begin{array}{c} v \\ - u \end{array} \right].
$$

(4) If S is symplectic and $\lambda \in \lambda ( S )$ , then $1 / \lambda \in \lambda ( S )$ :

$$
S \left[ \begin{array}{l} u \\ v \end{array} \right] = \lambda \left[ \begin{array}{l} u \\ v \end{array} \right] \quad \Rightarrow \quad S ^ {T} \left[ \begin{array}{l} v \\ - u \end{array} \right] = \frac {1}{\lambda} \left[ \begin{array}{l} v \\ - u \end{array} \right].
$$

Symplectic versions of Householder and Givens transformations have a prominanent role to play in Hamiltonian matrix computations. If $P = I _ { n } - 2 v v ^ { T }$ is a Householder matrix, then diag(P, P ) is a symplectic orthogonal matrix. Likewise, if $G \in \mathbb { R } ^ { 2 n \times 2 n }$ is a Givens rotation that involves planes i and $i + n$ , then G is a symplectic orthogonal matrix. Combinations of these transformations can be used to introduce zeros. For example, a Householder-Givens-Householder sequence can do this:

$$
\left[ \begin{array}{c} \times \\ \times \\ \times \\ \times \\ \hline \times \\ \times \\ \times \\ \times \\ \end{array} \right] \quad \stackrel {{\mathrm{diag} (P _ {1}, P _ {1})}} {{\longrightarrow}} \quad \left[ \begin{array}{c} \times \\ \times \\ \times \\ \times \\ \hline \times \\ 0 \\ 0 \\ 0 \end{array} \right] \quad \stackrel {{G _ {1, 5}}} {{\longrightarrow}} \quad \left[ \begin{array}{c} \times \\ \times \\ \times \\ \times \\ \hline 0 \\ 0 \\ 0 \\ 0 \end{array} \right] \quad \stackrel {{\mathrm{diag} (P _ {2}, P _ {2})}} {{\longrightarrow}} \quad \left[ \begin{array}{c} \times \\ 0 \\ 0 \\ 0 \\ \hline 0 \\ 0 \\ 0 \\ 0 \end{array} \right].
$$

This kind of vector reduction can be sequenced to produce a constructive proof of a structured Schur decomposition for Hamiltonian matrices. Suppose λ is a real eigenvalue of a Hamiltonian matrix M and that $x \in \mathbb { R } ^ { 2 n }$ is a unit 2-norm vector with $M x = \lambda x$ . If $Q _ { 1 } \in \mathbb { R } ^ { 2 n \times 2 n }$ is an orthogonal symplectic matrix and $Q _ { 1 } ^ { T } x = e _ { 1 }$ , then it follows from $( Q _ { 1 } ^ { T } M Q _ { 1 } ) ( Q _ { 1 } ^ { T } x ) = \lambda ( Q _ { 1 } ^ { T } x )$ that

$$
Q _ {1} ^ {T} M Q _ {1} = \left[ \begin{array}{c c c c c c c c} \lambda & \times & \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times & \times & \times \\ \hline 0 & 0 & 0 & 0 & - \lambda & 0 & 0 & 0 \\ 0 & \times & \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times & \times & \times \end{array} \right].
$$

The “extra” zeros follow from the Hamiltonian structure of $Q _ { 1 } ^ { T } M Q _ { 1 }$ . The process can be repeated on the 6-by-6 Hamiltonian submatrix defined by rows and columns 2-3-4- 6-7-8. Together with the assumption that M has no purely imaginary eigenvalues, it is possible to show that an orthogonal symplectic matrix Q exists so that

$$
Q ^ {T} M Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ - Q _ {2} & Q _ {1} \end{array} \right] ^ {T} \left[ \begin{array}{c c} A & F \\ G & - A ^ {T} \end{array} \right] \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ - Q _ {2} & Q _ {1} \end{array} \right] = \left[ \begin{array}{c c} T & R \\ 0 & - T ^ {T} \end{array} \right] \tag {7.8.1}
$$

where $T \in \mathbb { R } ^ { n \times n }$ is upper quasi-triangular. This is the real Hamiltonian-Schur decomposition. See Paige and Van Loan (1981) and, for a more general version, Lin, Mehrmann, and Xu (1999).

One reason that the Hamiltonian eigenvalue problem is so important is its connection to the algebraic Ricatti equation

$$
G + X A + A ^ {T} X - X F X = 0. \tag {7.8.2}
$$

This quadratic matrix problem arises in optimal control and a symmetric solution is sought so that the eigenvalues of $A - F X$ are in the open left half plane. Modest assumptions typically ensure that M has no eigenvalues on the imaginary axis and that the matrix $Q _ { 1 }$ in (7.8.1) is nonsingular. If we compare (2,1) blocks in (7.8.1), then

$$
Q _ {2} ^ {T} A Q _ {1} - Q _ {2} ^ {T} F Q _ {2} + Q _ {1} ^ {T} G Q _ {1} + Q _ {1} ^ {T} A ^ {T} Q _ {2} = 0.
$$

It follows from $I _ { n } = Q _ { 1 } ^ { T } Q _ { 1 } + Q _ { 2 } ^ { T } Q _ { 2 }$ that $X = Q _ { 2 } Q _ { 1 } ^ { - 1 }$ is symmetric and that it satisfies (7.8.2). From (7.8.1) it is easy to show that $A { \bf - } ^ { \prime } F X = Q _ { 1 } T Q _ { 1 } ^ { - 1 }$ and so the eigenvalues of $A - F X$ are the eigenvalues of T . It follows that the desired solution to the algebraic Ricatti equation can be obtained by computing the real Hamiltonian-Schur decomposition and ordering the eigenvalues so that $\lambda ( T )$ is in the left half plane.

How might the real Hamiltonian-Schur form be computed? One idea is to reduce $M$ to some condensed Hamiltonian form and then devise a structure-preserving QRiteration. Regarding the former task, it is easy to compute an orthogonal symplectic $U _ { 0 }$ so that

$$
U _ {0} ^ {T} M U _ {0} = \left[ \begin{array}{c c} H & R \\ D & - H ^ {T} \end{array} \right] \tag {7.8.3}
$$

where $H \in \mathbb { R } ^ { n \times n }$ is upper Hessenberg and D is diagonal. Unfortunately, a structurepreserving QR iteration that maintains this condensed form has yet to be devised. This impasse prompts consideration of methods that involve the skew-Hamiltonian matrix $N = M ^ { 2 }$ . Because the (2,1) block of a skew-Hamiltonian matrix is skew-symmetric, it has a zero diagonal. Symplectic similarity transforms preserve skew-Hamiltonian structure, and it is straightforward to compute an orthogonal symplectic matrix $V _ { 0 }$ such that

$$
V _ {0} ^ {T} M ^ {2} V _ {0} = \left[ \begin{array}{c c} H & R \\ 0 & H ^ {T} \end{array} \right], \tag {7.8.4}
$$

where H is upper Hessenberg. If $U ^ { T } H U = T$ is the real Schur form of H and and $Q = V _ { 0 } \cdot \mathrm { d i a g } ( U , U )$ , then

$$
Q ^ {T} M ^ {2} Q = \left[ \begin{array}{c c} T & U ^ {T} R U \\ 0 & T ^ {T} \end{array} \right]
$$

is the real skew-Hamiltonian Schur form. See Van Loan (1984). It does not follow that $Q ^ { T } M Q$ is in Schur-Hamiltonian form. Moreover, the quality of the computed small eigenvalues is not good because of the explicit squaring of M. However, these shortfalls can be overcome in an efficient numerically sound way, see Chu, Lie, and Mehrmann (2007) and the references therein. Kressner (NMSE, p. 175–208) and Watkins (MEP, p. 319–341) have in-depth treatments of the Hamiltonian eigenvalue problem.

# 7.8.2 Product Eigenvalue Problems

Using SVD and QZ, we can compute the eigenvalues of $A ^ { T } A$ and $B ^ { - 1 } A$ without forming products or inverses. The intelligent computation of the Hamiltonian-Schur decomposition involves a correspondingly careful handling of the product M -times-M . In this subsection we further develop this theme by discussing various product decompositions. Here is an example that suggests how we might compute the Hessenberg decomposition of

$$
A = A _ {3} A _ {2} A _ {1}
$$

where $A _ { 1 } , A _ { 2 } , A _ { 3 } \in \mathbb { R } ^ { n \times n }$ . Instead of forming this product explicitly, we compute orthogonal $U _ { 1 } , U _ { 2 } , U _ { 3 } \in \mathbb { R } ^ { n \times n }$ such that

$$
U _ {1} ^ {T} A _ {3} U _ {3} = H _ {3} \quad \text {(upper Hessenberg)},
$$

$$
U _ {3} ^ {T} A _ {2} U _ {2} = T _ {2} \quad (\text { upper   triangular }), \tag {7.8.5}
$$

$$
U _ {2} ^ {T} A _ {1} U _ {1} = T _ {1} \quad \text {(upper triangular)}.
$$

It follows that

$$
U _ {1} ^ {T} A U _ {1} = (U _ {1} ^ {T} A _ {3} U _ {3}) (U _ {3} ^ {T} A _ {2} U _ {2}) (U _ {2} ^ {T} A _ {1} U _ {1}) = H _ {3} T _ {2} T _ {1}
$$

is upper Hessenberg. A procedure for doing this would start by computing the QR factorizations

$$
Q _ {2} ^ {T} A _ {1} = R _ {1}, \quad Q _ {3} ^ {T} (A _ {2} Q _ {2}) = R _ {2}.
$$

If $\tilde { A } _ { 3 } = A _ { 3 } Q _ { 3 }$ , then $A = \tilde { A } _ { 3 } R _ { 2 } R _ { 1 }$ . The next phase involves reducing ${ \tilde { A } } _ { 3 }$ to Hessenberg form with Givens transformations coupled with “bulge chasing” to preserve the triangular structures already obtained. The process is similar to the reduction of $A - \lambda B$ to Hessenberg-triangular form; see §7.7.4.

Now suppose we want to compute the real Schur form of A

$$
Q _ {1} ^ {T} A _ {3} Q _ {3} = T _ {3} \quad \text {(upper quasi - triangular)},
$$

$$
Q _ {3} ^ {T} A _ {2} Q _ {2} = T _ {2} \quad \text {(upper triangular)}, \tag {7.8.6}
$$

$$
Q _ {2} ^ {T} A _ {1} Q _ {1} = T _ {1} \quad \text {(upper triangular)},
$$

where $Q _ { 1 } , Q _ { 2 } , Q _ { 3 } \in \mathbb { R } ^ { n \times n }$ are orthogonal. Without loss of generality we may assume that $\{ A _ { 3 } , A _ { 2 } , A _ { 1 } \}$ is in Hessenberg-triangular-triangular form. Analogous to the QZ iteration, the next phase is to produce a sequence of converging triplets

$$
\{A _ {3} ^ {(k)}, A _ {2} ^ {(k)}, A _ {1} ^ {(k)} \} \rightarrow \{T _ {3}, T _ {2}, T _ {1} \} \tag {7.8.7}
$$

with the property that all the iterates are in Hessenberg-triangular-triangular form.

Product decompositions (7.8.5) and (7.8.6) can be framed as structured decompositions of block-cyclic 3-by-3 matrices. For example, if

$$
U = \left[ \begin{array}{c c c} U _ {1} & 0 & 0 \\ 0 & U _ {2} & 0 \\ 0 & 0 & U _ {3} \end{array} \right]
$$

then we have the following restatement of (7.8.5):

$$
U ^ {T} \left[ \begin{array}{c c c} 0 & 0 & A _ {3} \\ A _ {1} & 0 & 0 \\ 0 & A _ {2} & 0 \end{array} \right] U = \left[ \begin{array}{c c c} 0 & 0 & H _ {3} \\ T _ {1} & 0 & 0 \\ 0 & T _ {2} & 0 \end{array} \right] = \tilde {H}.
$$

Consider the zero-nonzero structure of this matrix for the case $n = 4$ :

$$
\tilde {H} = \left[ \begin{array}{c c c c c c c c c c c c} 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times \\ \hline \times & \times & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & \times & \times & \times & \times & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & \times & \times & \times & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & 0 & 0 \end{array} \right].
$$

Using the perfect shuffle $\mathcal { P } _ { 3 4 }$ (see §1.2.11) we also have

$$
\mathcal {P} _ {3 4} \tilde {H} \mathcal {P} _ {3 4} = \left[ \begin{array}{c c c c c c c c c c c c} 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times \\ \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 \\ 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 \\ \hline 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times \\ 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 \end{array} \right].
$$

Note that this is a highly structured 12-by-12 upper Hessenberg matrix. This connection makes it possible to regard the product-QR iteration as a structure-preserving

QR iteration. For a detailed discussion about this connection and its implications for both analysis and computation, see Kressner (NMSE, pp. 146–174) and Watkins(MEP, pp. 293–303). We mention that with the “technology” that has been developed, it is possible to solve product eigenvalue problems where the factor matrices that define A are rectangular. Square nonsingular factors can also participate through their inverses, $\mathrm { e . g . , } A = A _ { 3 } A _ { 2 } ^ { - 1 } A _ { 1 }$ .

# Problems

P7.8.1 What can you say about the eigenvalues and eigenvectors of a symplectic matrix?

P7.8.2 Suppose $S _ { 1 } , S _ { 2 } \in \mathbb { R } ^ { n \times n }$ are both skew-symmetric and let $A = S _ { 1 } S _ { 2 }$ . Show that the nonzero eigenvalues of A are not simple. How would you compute these eigenvalues?

P7.8.3 Relate the eigenvalues and eigenvectors of

$$
A = \left[ \begin{array}{c c c c} 0 & A _ {1} & 0 & 0 \\ 0 & 0 & A _ {2} & 0 \\ 0 & 0 & 0 & A _ {3} \\ A _ {4} & 0 & 0 & 0 \end{array} \right].
$$

to the eigenvalues and eigenvectors of $\tilde { A } = A _ { 1 } A _ { 2 } A _ { 3 } A _ { 4 }$ . Assume that the diagonal blocks are square.

# Notes and References for §7.8

The books by Kressner(NMSE) and Watkins (MEP) have chapters on product eigenvalue problems and Hamiltonian eigenvalue problems. The sometimes bewildering network of interconnections that exist among various structured classes of matrices is clarified in:

A. Bunse-Gerstner, R. Byers, and V. Mehrmann (1992). “A Chart of Numerical Methods for Structured Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 13, 419–453.

Papers concerned with the Hamiltonian Schur decomposition include:

A.J. Laub and K. Meyer (1974). “Canonical Forms for Symplectic and Hamiltonian Matrices,” J. Celestial Mechanics 9, 213–238.

C.C. Paige and C. Van Loan (1981). “A Schur Decomposition for Hamiltonian Matrices,” Lin. Alg. Applic. 41, 11–32.

V. Mehrmann (1991). Autonomous Linear Quadratic Control Problems, Theory and Numerical Solution, Lecture Notes in Control and Information Sciences No. 163, Springer-Verlag, Heidelberg.

W.-W. Lin, V. Mehrmann, and H. Xu (1999). “Canonical Forms for Hamiltonian and Symplectic Matrices and Pencils,” Lin. Alg. Applic. 302/303, 469–533.

Various methods for Hamiltonian eigenvalue problems have been devised that exploit the rich underlying structure, see:

C. Van Loan (1984). “A Symplectic Method for Approximating All the Eigenvalues of a Hamiltonian Matrix,” Lin. Alg. Applic. 61, 233–252.

R. Byers (1986) “A Hamiltonian QR Algorithm,” SIAM J. Sci. Stat. Comput. 7, 212–229.

P. Benner, R. Byers, and E. Barth (2000). “Algorithm 800: Fortran 77 Subroutines for Computing the Eigenvalues of Hamiltonian Matrices. I: the Square-Reduced Method,” ACM Trans. Math. Softw. 26, 49–77.

H. Fassbender, D.S. Mackey and N. Mackey (2001). “Hamilton and Jacobi Come Full Circle: Jacobi Algorithms for Structured Hamiltonian Eigenproblems,” Lin. Alg. Applic. 332-4, 37–80.

D.S. Watkins (2006). “On the Reduction of a Hamiltonian Matrix to Hamiltonian Schur Form,” ETNA 23, 141–157.

D.S. Watkins (2004). “On Hamiltonian and Symplectic Lanczos Processes,” Lin. Alg. Applic. 385, 23–45.

D. Chu, X. Liu, and V. Mehrmann (2007). “A Numerical Method for Computing the Hamiltonian Schur Form,” Numer. Math. 105, 375–412.

Generalized eigenvalue problems that involve Hamiltonian matrices also arise:

P. Benner, V. Mehrmann, and H. Xu (1998). “A Numerically Stable, Structure Preserving Method for Computing the Eigenvalues of Real Hamiltonian or Symplectic Pencils,” Numer. Math. 78, 329–358.


---

<!-- golub_450_499 -->

C. Mehl (2000). “Condensed Forms for Skew-Hamiltonian/Hamiltonian Pencils,” SIAM J. Matrix Anal. Applic. 21, 454–476.   
V. Mehrmann and D.S. Watkins (2001). “Structure–Preserving Methods for Computing Eigenpairs of Large Sparse Skew–Hamiltonian/Hamiltonian Pencils,” SIAM J. Sci. Comput. 22, 1905–1925.   
P. Benner and R. Byers, V. Mehrmann, and H. Xu (2002). “Numerical Computation of Deflating Subspaces of Skew-Hamiltonian/Hamiltonian Pencils,” SIAM J. Matrix Anal. Applic. 24, 165– 190.

Methods for symplectic eigenvalue problems are discussed in:

P. Benner, H. Fassbender and D.S. Watkins (1999). “SR and SZ Algorithms for the Symplectic (Butterfly) Eigenproblem,” Lin. Alg. Applic. 287, 41–76.

The Golub-Kahan SVD algorithm that we discuss in the next chapter does not form $A ^ { T } A$ or $A A ^ { T }$ despite the rich connection to the Schur decompositions of those matrices. From that point on there has been an appreciation for the numerical dangers associated with explicit products. Here is a sampling of the literature:

C. Van Loan (1975). “A General Matrix Eigenvalue Algorithm,” SIAM J. Numer. Anal. 12, 819–834.

M.T. Heath, A.J. Laub, C.C. Paige, and R.C. Ward (1986). “Computing the SVD of a Product of Two Matrices,” SIAM J. Sci. Stat. Comput. 7, 1147–1159.

R. Mathias (1998). “Analysis of Algorithms for Orthogonalizing Products of Unitary Matrices,” Num. Lin. Alg. 3, 125–145.

G. Golub, K. Solna, and P. Van Dooren (2000). “Computing the SVD of a General Matrix Product/Quotient,” SIAM J. Matrix Anal. Applic. 22, 1–19.

D.S. Watkins (2005). “Product Eigenvalue Problems,” SIAM Review 47, 3–40.

R. Granat and B. Kgstrom (2006). “Direct Eigenvalue Reordering in a Product of Matrices in Periodic Schur Form,” SIAM J. Matrix Anal. Applic. 28, 285–300.

Finally we mention that there is a substantial body of work concerned with structured error analysis and structured perturbation theory for structured matrix problems, see:

F. Tisseur (2003). “A Chart of Backward Errors for Singly and Doubly Structured Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 24, 877–897.

R. Byers and D. Kressner (2006). “Structured Condition Numbers for Invariant Subspaces,” SIAM J. Matrix Anal. Applic. 28, 326–347.

M. Karow, D. Kressner, and F. Tisseur (2006). “Structured Eigenvalue Condition Numbers,” SIAM J. Matrix Anal. Applic. 28, 1052–1068.
