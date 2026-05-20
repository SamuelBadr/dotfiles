# 12.2.10 Eigenvalues of an Orthogonal Upper Hessenberg Matrix

We close with an eigenvalue problem that has quasiseparable structure. Suppose $H \in \mathbb { R } ^ { n \times n }$ is an upper Hessenberg matrix that is also orthogonal. Our goal is to compute $\lambda ( H )$ . Note that each eigenvalue is on the unit circle. Without loss of generality we may assume that the subdiagonal entries are nonzero.

If n is odd, then it must have a real eigenvalue because the eigenvalues of a real matrix come in complex conjugate pairs. In this case it is possible to deflate the problem by carefully working with the eigenvector equation $H x = x { \mathrm { ~ ( o r ~ } } H x = - x { \mathrm { ) } }$ . Thus, we may assume that n is even.

For $1 \leq k \leq n - 1$ , define the reflection $G _ { k } \in \mathbb { R } ^ { n \times n }$ by

$$
G _ {k} = G (\phi_ {k}) = \mathrm{diag} (I _ {k - 1}, R (\phi_ {k}), I _ {n - k - 1})
$$

where

$$
R (\phi_ {k}) = \left[ \begin{array}{c c} - \cos (\phi_ {k}) & \sin (\phi_ {k}) \\ \sin (\phi_ {k}) & \cos (\phi_ {k}) \end{array} \right], \qquad 0 <   \phi_ {k} <   \pi .
$$

These transformations can be used to represent the QR factorization of H. Indeed, as for the Givens process described in §5.2.6, we can compute $G _ { 1 } , \ldots , G _ { n - 1 }$ so that

$$
G _ {n - 1} \dots G _ {1} H = G _ {n} \equiv \operatorname{diag} (1, \dots , 1, - c _ {n}).
$$

The matrix $G _ { n }$ is the $^ { 6 } R ^ { 5 }$ matrix. It is diagonal because an orthogonal upper triangular matrix must be diagonal. Since the determinant of a matrix is the product of its eigenvalues, the value of $c _ { n }$ is either +1 or −1. If $c _ { n } = - 1$ , then det $( H ) = - 1$ , which in turn implies that H has a real eigenvalue and we can deflate the problem. Thus, we may assume that

$$
H = G _ {1} \dots G _ {n}, \quad G _ {n} = \operatorname{diag} (1, \dots , 1, - 1), \quad n = 2 m \tag {12.2.18}
$$

and that our goal is to compute

$$
\lambda (H) = \{\cos (\theta_ {1}) \pm i \cdot \sin (\theta_ {1}), \dots , \cos (\theta_ {m}) \pm i \cdot \sin (\theta_ {m}) \}. \tag {12.2.19}
$$

Note that (12.2.4) and (12.2.18) tell us that H is quasiseparable.

Ammar, Gragg, and Reichel (1986) propose an interesting $O ( n ^ { 2 } )$ method that computes the required eigenvalues by setting up a pair of m-by-m bidiagonal SVD problems. Three facts are required:

Fact 1. H is similar to $\tilde { H } = H _ { o } H _ { e }$ where

$$
H _ {o} = G _ {1} G _ {3} \dots G _ {n - 1} = \operatorname{diag} (R (\phi_ {1}), R (\phi_ {3}), \dots , R (\phi_ {n - 1})),
$$

$$
H _ {e} = G _ {2} G _ {4} \dots G _ {n} = \operatorname{diag} (1, R (\phi_ {2}), R (\phi_ {4}), \dots , R (\phi_ {n - 2}), - 1).
$$

Fact 2. The matrices

$$
C = \frac {H _ {o} + H _ {e}}{2}, \qquad S = \frac {H _ {o} - H _ {e}}{2}
$$

are symmetric and tridiagonal. Moreover, their eigenvalues are given by

$$
\lambda (C) = \{\pm \cos (\theta_ {1} / 2), \dots , \pm \cos (\theta_ {m} / 2) \},
$$

$$
\lambda (S) = \{\pm \sin (\theta_ {1} / 2), \dots , \pm \sin (\theta_ {m} / 2) \}.
$$

Fact 3. If

$$
Q _ {o} = \operatorname{diag} (R (\phi_ {1} / 2), R (\phi_ {3} / 2), \dots , R (\phi_ {n - 1} / 2)),
$$

$$
Q _ {e} = \mathrm{diag} (1, R (\phi_ {2} / 2), R (\phi_ {4} / 2), \ldots , R (\phi_ {n - 2} / 2), - 1),
$$

then perfect shuffle permutations of the matrices

$$
C ^ {(1)} = Q _ {o} C Q _ {e}, \qquad S ^ {(1)} = Q _ {o} S Q _ {e}
$$

expose a pair of m-by-m bidiagonal matrices $B _ { c }$ and $B _ { s }$ with the property that

$$
\sigma (B _ {c}) = \{\cos (\theta_ {1} / 2), \ldots , \cos (\theta_ {m} / 2) \},
$$

$$
\sigma (B _ {s}) = \{\sin (\theta_ {1} / 2), \dots , \sin (\theta_ {m} / 2) \}.
$$

Once the bidiagonal matrices $B _ { c }$ and $B _ { s }$ are set up (which involves $O ( n )$ work), then their singular values can be computed via Golub-Kahan SVD algorithm. The angle $\theta _ { k }$ can be accurately determined from sin $\left( \theta _ { k } / 2 \right)$ if $0 < \theta _ { k } < \pi / 2$ and from $\cos ( \theta _ { k } / 2 )$ otherwise. See Ammar, Gragg, and Reichel (1986) for more details.

# Problems

P12.2.1 Rigorously prove that the matrix $\boldsymbol { B } ( \boldsymbol { r } ) ^ { - 1 }$ is semiseparable.

P12.2.2 Prove that A is quasiseparable if and only if $A = \mathbf { S } ( u , t , v , d , p , r , q )$ for appropriately chosen vectors u, v, t, d, p, r, and q.

P12.2.3 How many flops are required to execute the n-by-n matrix vector product y = Ax where $A = \mathbf { S } ( u , v , t , d , p , q , r )$ .

P12.2.4 Refer to (12.2.4). Determine u, v, t, d, p, q, and r so that $M = \mathbf { S } ( u , v , t , d , p , q , r )$

P12.2.5 Suppose $\mathbf { S } ( u , v , t , d , v , u , t )$ is symmetric positive definite and semiseparable. Show that its Cholesky factor is semiseparable and give an algorithm for computing its quasiseparable representation.

P12.2.6 Verify the three facts in §12.2.3.

P12.2.7 Develop a fast method for solving the upper triangular system $T x = y$ where $_ T$ is the matrix $T = \mathrm { d i a g } ( d ) + \mathrm { t r i u } ( p q ^ { T } , 1 ) \ . * B ( r ) ^ { - 1 }$ with $p , q , d , y \in \mathbb { R } ^ { n }$ and $r \in \mathbb { R } ^ { n - 1 }$ .

P12.2.8 Verify (12.2.7).

P12.2.9 Prove (12.2.14).

P12.2.10 Assume that A is an N-by-N block matrix that has the sequentially separable structure illustrated in $( 1 2 . 2 . 1 7 )$ . Assume that the blocks are each $m { \mathrm { - } } \mathrm { b y } - m$ . Give a fast algorithm for computing y = Ax where x ∈ IRNm. $y = A x$ $\boldsymbol { x } \in \mathbb { R } ^ { N m }$

P12.2.11 It can be shown that

$$
A = \left[ \begin{array}{c c c c} A _ {1} & B _ {1} ^ {T} & 0 & 0 \\ B _ {1} & A _ {2} & B _ {2} ^ {T} & 0 \\ 0 & B _ {2} & A _ {3} & B _ {3} ^ {T} \\ 0 & 0 & B _ {3} & A _ {4} \end{array} \right] \Rightarrow A ^ {- 1} = \left[ \begin{array}{c c c c} U _ {1} V _ {1} ^ {T} & V _ {1} U _ {2} ^ {T} & V _ {1} U _ {3} ^ {T} & V _ {1} U _ {4} ^ {T} \\ U _ {2} V _ {1} ^ {T} & U _ {2} V _ {2} ^ {T} & V _ {2} U _ {3} ^ {T} & V _ {2} U _ {4} ^ {T} \\ U _ {3} V _ {1} ^ {T} & U _ {3} V _ {2} ^ {T} & U _ {3} V _ {3} ^ {T} & V _ {3} U _ {4} ^ {T} \\ U _ {4} V _ {1} ^ {T} & U _ {4} V _ {2} ^ {T} & U _ {4} V _ {3} ^ {T} & U _ {4} V _ {4} ^ {T} \end{array} \right],
$$

assuming that A is symmetric positive definite and that the $B _ { i }$ are nonsingular. Give an algorithm that computes $U _ { 1 } , \dots , U _ { 4 }$ and $V _ { 1 } , \ldots , V _ { 4 }$ .

P12.2.12 Suppose $a , b , f , g \in \mathbb { R } ^ { n }$ and that $A = \mathsf { t r i u } ( a b ^ { T } + f g ^ { T } )$ is nonsingular. (a) Given $x \in \mathbb { R } ^ { n }$ , show how to compute efficiently $y = A x$ . (b) Given $\boldsymbol { y } \in \mathbb { R } ^ { n }$ , show how to compute $\boldsymbol { x } \in \mathbb { R } ^ { n }$ so that $A x = y . \ \mathrm { ( c ) }$ Given $y , d \in \mathbb { R } ^ { n }$ , show how to compute x so that $y = ( A + D ) x$ where it is assumed that $D = \mathrm { d i a g } ( d )$ and $A + D$ are nonsingular.

P12.2.13 Verify the three facts in 12.2.10 for the case $n = 8 .$

P12.2.14 Show how to compute the eigenvalues of an orthogonal matrix $A \in \mathbb { R } ^ { n \times n }$ by computing the Schur decompositions of $( A + A ^ { T } ) / 2$ and $( A - A ^ { T } ) / 2$ .

# Notes and References for §12.2

For all matters concerning structured rank matrix computations, see:

R. Vandebril, M. Van Barel, and N. Mastronardi (2008). Matrix Computations and Semiseparable Matrices, Vol. I Linear Systems, Johns Hopkins University Press, Baltimore, MD.

R. Vandebril, M. Van Barel, and N. Mastronardi (2008). Matrix Computations and Semiseparable Matrices, Vol. II Eigenvalue and Singular Value Methods, Johns Hopkins University Press, Baltimore, MD.

As we have seen, working with the “right” representation is critically important in order to realize an efficient implementation. For more details, see:

R. Vandebril, M. Van Barel, and N. Mastronardi (2005). “A Note on the Representation and Definition of Semiseparable Matrices,” Num. Lin. Alg. Applic. 12, 839–858.

References concerned with the fast solution of linear equations and least squares problems with structured rank include:

I. Gohberg, T. Kailath, and I Koltracht (1985) “Linear Complexity Algorithm for Semiseparable Matrices,” Integral Equations Operator Theory 8, 780–804.

Y. Eidelman and I. Gohberg (1997). “Inversion Formulas and Linear Complexity Algorithm for Diagonal-Plus-Semiseparable Matrices,” Comput. Math. Applic. 33, 69–79.

P. Dewilde and A.J. van der Veen (1998). Time-Varying Systems and Computations, Kluwer Academic, Boston, MA,

S. Chandrasekaran and M. Gu (2003). “Fast and Stable Algorithms for Banded-Plus-Semiseparable Systems of Linear Equations,” SIAM J. Matrix Anal. Applic. 25, 373–384.

S. Chandrasekaran, P. Dewilde, M. Gu, T. Pals, X. Sun, A.J. Van Der Veen, and D. White (2005). “Some Fast Algorithms for Sequentially Semiseparable Representations,” SIAM J. Matrix Anal. Applic. 27, 341–364.

E. Van Camp, N. Mastronardi, and M. Van Barel (2004). “Two Fast Algorithms for Solving Diagonal-Plus-Semiseparable Linear Systems,” J. Comput. Appl. Math. 164, 731–747.   
T. Bella, Y. Eidelman, I. Gohberg, V. Koltracht, and V. Olshevsky (2009). “A Fast Bjorck-Pereyra-Type Algorithm for Solving Hessenberg-Quasiseparable-Vandermonde Systems SIAM. J. Matrix Anal. Applic. 31, 790–815.   
J. Xia and M. Gu (2010). “Robust Approximate Cholesky Factorization of Rank-Structured Symmetric Positive Definite Matrices,” SIAM J. Matrix Anal. Applic. 31, 2899–2920.   
For discussion of methods that exploit hierarchical rank structure, see:   
S. B¨orm, L. Grasedyck, and W. Hackbusch (2003). “Introduction to Hierarchical Matrices with Applications,” Engin. Anal. Boundary Elements 27, 405–422.   
S. Chandrasekaran, M. Gu, and T. Pals (2006). “A Fast ULV Decomposition Solver for Hierarchically Semiseparable Representations,” SIAM J. Matrix Anal. Applic. 28, 603–622.   
S. Chandrasekaran, M. Gu, X. Sun, J. Xia, and J. Zhu (2007). “A Superfast Algorithm for Toeplitz Systems of Linear Equations,” SIAM J. Matrix Anal. Applic. 29, 1247–1266.   
S. Chandrasekaran, M. Gu, J. Xia, and J. Zhu (2007). “A Fast QR Algorithm for Companion Matrices,” Oper. Theory Adv. Applic. 179, 111–143.   
J. Xia, S. Chandrasekaran, M. Gu, and X.S. Li (2010). “Fast algorithms for Hierarchically Semiseparable Matrices,” Numer. Lin. Alg. Applic. 17, 953–976.   
S. Chandrasekaran, P. Dewilde, M. Gu, and N. Somasunderam (2010). “On the Numerical rank of the Off-Diagonal Blocks of Schur Complements of Discretized Elliptic PDEs,” SIAM J. Matrix Anal. Applic. 31, 2261–2290.   
P.G. Martinsson (2011). “A Fast Randomized Algorithm for Computing a Hierarchically Semi-Separable Representation of a Matrix,” SIAM J. Matrix Anal. Applic. 32, 1251–1274.   
J. Xia (2012). “On the Complexity of Some Hierarchical Structured Matrix Algorithms,” SIAM J. Matrix Anal. Applic. 33, 388–410.   
Reductions to tridiagonal, bidiagonal, and Hessenberg form are essential “front ends” for many eigenvalue and singular value procedures. There are ways to proceed when rank structure is present, see:   
N. Mastronardi, S. Chandrasekaran, and S. van Huffel (2001). “Fast and Stable Reduction of Diagonal Plus Semi-Separable Matrices to Tridiagonal and Bidiagonal Form,” BIT 41, 149–157.   
M. Van Barel, R. Vandebril, and N. Mastronardi (2005). “An Orthogonal Similarity Reduction of a Matrix into Semiseparable Form,” SIAM J. Matrix Anal. Applic. 27, 176–197.   
M. Van Barel, E. Van Camp, N. Mastronardi (2005). “Orthogonal Similarity Transformation into Block-Semiseparable Matrices of Semiseparability Rank,” Num. Lin. Alg. 12, 981–1000.   
R. Vandebril, E. Van Camp, M. Van Barel, and N. Mastronardi (2006). “Orthogonal Similarity Transformation of a Symmetric Matrix into a Diagonal-Plus-Semiseparable One with Free Choice of the Diagonal,” Numer. Math. 102, 709–726.   
Y. Eidelman, I. Gohberg, and L. Gemignani (2007). “On the Fast reduction of a Quasiseparable Matrix to Hessenberg and Tridiagonal Forms,” Lin. Alg. Applic. 420, 86–101.   
R. Vandebril, E. Van Camp, M. Van Barel, and N. Mastronardi (2006). “On the Convergence Properties of the Orthogonal Similarity Transformations to Tridiagonal and Semiseparable (Plus Diagonal) Form,” Numer. Math. 104, 205–239.   
Papers concerned with various structured rank eigenvalue iterations include:   
R. Vandebril, M. Van Barel, and N. Mastronardi (2004). “A QR Method for Computing the Singular Values via Semiseparable Matrices,” Numer. Math. 99, 163–195.   
R. Vandebril, M. Van Barel, N. Mastronardi (2005). “An Implicit QR algorithm for Symmetric Semiseparable Matrices,” Num. Lin. Alg. 12, 625–658.   
N. Mastronardi, E. Van Camp, and M. Van Barel (2005). “Divide and Conquer Algorithms for Computing the Eigendecomposition of Symmetric Diagonal-plus-Semiseparable Matrices,” Numer. Alg. 39, 379–398.   
Y. Eidelman, I. Gohberg, and V. Olshevsky (2005). “The QR Iteration Method for Hermitian Quasiseparable Matrices of an Arbitrary Order,” Lin. Alg. Applic. 404, 305–324.   
Y. Vanberghen, R. Vandebril, M. Van Barel (2008). “A QZ-Method Based on Semiseparable Matrices,” J. Comput. Appl. Math. 218, 482–491.   
M. Van Barel, Y. Vanberghen, and P. Van Dooren (2010). “Using Semiseparable Matrices to Compute the SVD of a General Matrix Product/Quotient,” J. Comput. Appl. Math. 234, 3175–3180.   
Our discussion of the orthogonal matrix eigenvalue problem is based on:

G.S. Ammar, W.B. Gragg, and L. Reichel (1985). “On the Eigenproblem for Orthogonal Matrices,” Proc. IEEE Conference on Decision and Control, 1963–1966.

There is an extensive literature concerned with unitary/orthogonal eigenvalue problem including:

P.J. Eberlein and C.P. Huang (1975). “Global Convergence of the QR Algorithm for Unitary Matrices with Some Results for Normal Matrices,” SIAM J. Numer. Anal. 12, 421–453.

A. Bunse-Gerstner and C. He (1995). “On a Sturm Sequence of Polynomials for Unitary Hessenberg Matrices,” SIAM J. Matrix Anal. Applic. 16, 1043–1055.

B. Bohnhorst, A. Bunse-Gerstner, and H. Fassbender (2000). “On the Perturbation Theory for Unitary Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 21, 809–824.

M. Gu, R. Guzzo, X.-B. Chi, and X.-O. Cao (2003). “A Stable Divide and Conquer Algorithm for the Unitary Eigenproblem,” SIAM J. Matrix Anal. Applic. 25, 385–404.

M. Stewart (2006). “An Error Analysis of a Unitary Hessenberg QR Algorithm,” SIAM J. Matrix Anal. Applic. 28, 40–67.

R.J.A. David and D.S. Watkins (2006). “Efficient Implementation of the Multishift QR Algorithm for the Unitary Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 28, 623–633.

For a nice introduction to this problem, see Watkins (MEP, pp. 341–346).
