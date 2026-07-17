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

# 4.1 Diagonal Dominance and Symmetry

Pivoting is a serious concern in the context of high-performance computing because the cost of moving data around rivals the cost of computation. Equally important, pivoting can destroy exploitable structure. For example, if A is symmetric, then it involves half the data of a general A. Our intuition (correctly) tells us that we should be able to solve a symmetric Ax = b problem with half the arithmetic. However, in the context of Gaussian elimination with pivoting, symmetry can be destroyed at the very start of the reduction, e.g.,

$$
\left[ \begin{array}{c c c} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{array} \right] \left[ \begin{array}{c c c} a & b & c \\ b & d & e \\ c & e & f \end{array} \right] = \left[ \begin{array}{c c c} c & e & f \\ b & d & e \\ a & b & c \end{array} \right].
$$

Taking advantage of symmetry and other patterns and identifying situations where pivoting is unnecessary are typical activities in the realm of structured Ax = b solving. The goal is to expose computational shortcuts and to justify their use through analysis.

# 4.1.1 Diagonal Dominance and the LU Factorization

If A’s diagonal entries are large compared to its off-diagonal entries, then we anticipate that it is safe to compute $A = L U$ without pivoting. Consider the $n = 2$ case:

$$
\left[ \begin{array}{c c} a & b \\ c & d \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ c / a & 1 \end{array} \right] \left[ \begin{array}{c c} a & b \\ 0 & d - (c / a) b \end{array} \right].
$$

If a and d “dominate” b and c in magnitude, then the elements of L and U will be nicely bounded. To quantify this we make a definition. We say that $A \in \mathbb { R } ^ { n \times n }$ is row diagonally dominant if

$$
\left|a_{ii}\right|\geq \sum_{\substack{j = 1\\ j\neq i}}^{n}\left|a_{ij}\right|,\qquad i = 1:n. \tag{4.1.1}
$$

Similarly, column diagonal dominance means that $| a _ { j j } |$ is larger than the sum of all off-diagonal element magnitudes in the same column. If these inequalities are strict, then A is strictly (row/column) diagonally dominant. A diagonally dominant matrix can be singular, e.g., the 2-by-2 matrix of 1’s. However, if a nonsingular matrix is diagonally dominant, then it has a “safe” LU factorization.

Theorem 4.1.1. If A is nonsingular and column diagonally dominant, then it has an LU factorization and the entries in $L = ( \ell _ { i j } )$ satisfy $| l _ { i j } | \le 1$ .

Proof. We proceed by induction. The theorem is obviously true if $n = 1$ . Assume that it is true for $( n - 1 ) – \mathrm { b y } – ( n - 1 )$ nonsingular matrices that are column diagonally dominant. Partition $A \in \mathbb { R } ^ { n \times n }$ as follows:

$$
A   =   \left[ \begin{array}{c c} \alpha & w ^ {T} \\ v & C \end{array} \right], \qquad \alpha \in \mathbb {R},   v, w \in \mathbb {R} ^ {n - 1},   C \in \mathbb {R} ^ {(n - 1) \times (n - 1)}.
$$

If $\alpha = 0$ , then $v = 0$ and A is singular. Thus, $\alpha \neq 0$ and we have the factorization

$$
\left[ \begin{array}{c c} \alpha & w ^ {T} \\ v & C \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & I \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ 0 & B \end{array} \right] \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & I \end{array} \right], \tag {4.1.2}
$$

where

$$
B = C - \frac {1}{\alpha} v w ^ {T}.
$$

Since det $( A ) = \alpha \cdot \mathsf { d e t } ( B )$ , it follows that B is nonsingular. It is also column diagonally dominant because

$$
\begin{array}{l} \sum_{\substack{i = 1\\ i\neq j}}^{n - 1}|b_{ij}| = \sum_{\substack{i = 1\\ i\neq j}}^{n - 1}|c_{ij} - v_{i}w_{j} / \alpha |\leq \sum_{\substack{i = 1\\ i\neq j}}^{n - 1}|c_{ij}| + \frac{|w_{j}|}{|\alpha|}\sum_{\substack{i = 1\\ i\neq j}}^{n - 1}|v_{i}| \\ <   (| c _ {j j} | - | w _ {j} |) + \frac {| w _ {j} |}{| \alpha |} (| \alpha | - | v _ {j} |) \leq \left| c _ {j j} - \frac {w _ {j} v _ {j}}{\alpha} \right| = | b _ {j j} |. \\ \end{array}
$$

By induction, B has an LU factorization $L _ { 1 } U _ { 1 }$ and so from (4.1.2) we have

$$
A   =   \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & L _ {1} \end{array} \right] \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & U _ {1} \end{array} \right]   \equiv   L U.
$$

The entries in $| v / \alpha |$ are bounded by 1 because A is column diagonally dominant. By induction, the same can be said about the entries in $| L _ { 1 } |$ . Thus, the entries in $| L |$ are all bounded by 1 completing the proof.

The theorem shows that Gaussian elimination without pivoting is a stable solution procedure for a column diagonally dominant matrix. If the diagonal elements strictly dominate the off-diagonal elements, then we can actually bound $\parallel A ^ { - 1 } \parallel$ .

Theorem 4.1.2. If $A \in \mathbb { R } ^ { n \times n }$ and

$$
\delta = \min_ {1 \leq j \leq n} \left(| a _ {j j} | - \sum_ {\substack {i = 1 \\ i \neq j}} ^ {n} | a _ {i j} |\right) > 0 \tag{4.1.3}
$$

then

$$
\| A ^ {- 1} \| _ {1} \leq 1 / \delta .
$$

Proof. Define $D = \mathrm { d i a g } ( a _ { 1 1 } , \ldots , a _ { n n } )$ and $E = A - D$ . If e is the column n-vector of 1’s, then

$$
e ^ {T} | E | \leq e ^ {T} | D | - \delta e ^ {T}.
$$

If $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , then $\begin{array} { r } { D x = A x - E x } \end{array}$ and

$$
| D | | x | \leq | A x | + | E | | x |.
$$

Thus,

$$
e ^ {T} | D | | x | \leq e ^ {T} | A x | + e ^ {T} | E | | x | \leq \| A x \| _ {1} + \left(e ^ {T} | D | - \delta e ^ {T}\right) | x |
$$

and so $\delta \parallel x \parallel _ { 1 } = \delta e ^ { T } | x | \leq \parallel A x \parallel _ { 1 }$ . The bound on $\left. \mathbf { A } ^ { - 1 } \right. _ { 1 }$ follows from the fact that for any $y \in \mathbb { R } ^ { n }$ ,

$$
\delta \| A ^ {- 1} y \| _ {1} \leq \| A (A ^ {- 1} y) \| _ {1} = \| y \| _ {1}. \quad \square
$$

The “dominance” factor δ defined in (4.1.3) is important because it has a bearing on the condition of the linear system. Moreover, if it is too small, then diagonal dominance may be lost during the elimination process because of roundoff. That is, the computed version of the B matrix in (4.1.2) may not be column diagonally dominant.

# 4.1.2 Symmetry and the $\mathbf { \mathbf { \mathbf { L } } } \mathbf { \mathbf { D } } \mathbf { \mathbf { \mathbf { L } } } ^ { T }$ Factorization

If A is symmetric and has an LU factorization $A = L U$ , then L and U have a connection. For example, if $n = 2$ we have

$$
\begin{array}{l} \left[ \begin{array}{c c} a & c \\ c & d \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ c / a & 1 \end{array} \right] \cdot \left[ \begin{array}{c c} a & c \\ 0 & d - (c / a) c \end{array} \right] \\ = \left[ \begin{array}{c c} 1 & 0 \\ c / a & 1 \end{array} \right] \cdot \left(\left[ \begin{array}{c c} a & 0 \\ 0 & d - (c / a) c \end{array} \right] \left[ \begin{array}{c c} 1 & c / a \\ 0 & 1 \end{array} \right]\right). \\ \end{array}
$$

It appears that U is a row scaling of $L ^ { T }$ . Here is a result that makes this precise.

Theorem 4.1.3. (LDLT Factorization) $I f A \in \mathbb { R } ^ { n \times n }$ is symmetric and the principal submatrix $A ( 1 { : } k , 1 { : } k )$ is nonsingular for $k = 1 { : } n - 1$ , then there exists a unit lower triangular matrix L and a diagonal matrix

$$
D = \operatorname{diag} (d _ {1}, \dots , d _ {n})
$$

such that $A = L D L ^ { T }$ . The factorization is unique.

Proof. By Theorem 3.2.1 we know that A has an LU factorization $A = L U$ . Since the matrix

$$
L ^ {- 1} A L ^ {- T} = U L ^ {- T}
$$

is both symmetric and upper triangular, it must be diagonal. The theorem follows by setting $D = U L ^ { - T }$ and the uniqueness of the LU factorization.

Note that once we have the $\mathrm { L D L } ^ { T }$ factorization, then solving $A x = b$ is a 3-step process:

$$
L z = b, \qquad D y = z, \qquad L ^ {T} x = y.
$$

This works because $A x = L ( D ( L ^ { T } x ) ) = L ( D y ) = L z = b .$

Because there is only one triangular matrix to compute, it is not surprising that the factorization $A = L D L ^ { T }$ requires half as many flops to compute as $A = L U$ . To see this we derive a Gaxpy-rich procedure that, for $j = 1 { : } n$ , computes $L ( j + 1 { : } n , j )$ and $d _ { j }$ in step j. Note that

$$
A (j: n, j) = L (j: n, 1: j) \cdot v (1: j)
$$

where

$$
v (1: j) = \left[ \begin{array}{c} d _ {1} \ell_ {j 1} \\ d _ {2} \ell_ {j 2} \\ \vdots \\ d _ {j - 1} \ell_ {j, j - 1} \\ d _ {j} \end{array} \right].
$$

From this we conclude that

$$
d _ {j} = a _ {j j} - \sum_ {k = 1} ^ {j - 1} d _ {k} \ell_ {j k} ^ {2}.
$$

With $d _ { j }$ available, we can rearrange the equation

$$
\begin{array}{l} A (j + 1: n, j) = L (j + 1: n, 1: j) \cdot v (1: j) \\ = L (j + 1: n, 1: j - 1) \cdot v (1: j - 1) + d _ {j} \cdot L (j + 1: n, j) \\ \end{array}
$$

to get a recipe for $L ( j + 1 { : } n , j )$ :

$$
L (j + 1: n, j) = \frac {1}{d _ {j}} \left(A (j + 1: n, j) - L (j + 1: n, 1: j - 1) \cdot v (1: j - 1)\right).
$$

Properly sequenced, we obtain the following overall procedure:

for $j = 1:n$ for $i = 1:j - 1$ $v(i) = L(j,i) \cdot d(i)$ end $d(j) = A(j,j) - L(j,1:j - 1) \cdot v(1:j - 1)$ $L(j + 1:n,j) = (A(j + 1:n,j) - L(j + 1:n,1:j - 1) \cdot v(1:j - 1)) / d(j)$ end

With overwriting we obtain the following procedure.

Algorithm 4.1.1 (LDLT) If $A \in \mathbb { R } ^ { n \times n }$ is symmetric and has an LU factorization, then this algorithm computes a unit lower triangular matrix L and a diagonal matrix $D =$ diag $( d _ { 1 } , \ldots , d _ { n } )$ so $\bar { \boldsymbol { A } } = \boldsymbol { L } \boldsymbol { D } \boldsymbol { L } ^ { T }$ . The entry $a _ { i j }$ is overwritten with $\ell _ { i j } { \mathrm { ~ i f ~ } } i > j$ and with $d _ { i }$ if $i = j$ .

for $j = 1:n$ for $i = 1:j - 1$ $v(i) = A(j,i)A(i,i)$ end $A(j,j) = A(j,j) - A(j,1:j - 1) \cdot v(1:j - 1)$ $A(j + 1:n,j) = (A(j + 1:n,j) - A(j + 1:n,1:j - 1) \cdot v(1:j - 1)) / A(j,j)$

end

This algorithm requires $n ^ { 3 } / 3$ flops, about half the number of flops involved in Gaussian elimination.

The computed solution $\hat { x }$ to Ax = b obtained via Algorithm 4.1.1 and the usual triangular system solvers of §3.1 can be shown to satisfy a perturbed system $( A + E ) { \hat { x } } =$ $b ,$ where

$$
| E | \leq n \mathbf {u} \left(2 | A | + 4 | \hat {L} | | \hat {D} | | \hat {L} ^ {T} |\right) + O (\mathbf {u} ^ {2}) \tag {4.1.4}
$$

and $\hat { L }$ and $\hat { D }$ are the computed versions of L and D, respectively.

As in the case of the LU factorization considered in the previous chapter, the upper bound in (4.1.4) is without limit unless A has some special property that guarantees stability. In the next section, we show that if A is symmetric and positive definite, then Algorithm 4.1.1 not only runs to completion, but is extremely stable. If $A$ is symmetric but not positive definite, then, as we discuss in §4.4, it is necessary to consider alternatives to the $\mathrm { L D L } ^ { T }$ factorization.

# Problems

P4.1.1 Show that if all the inequalities in (4.1.1) are strict inequalities, then A is nonsingular.

P4.1.2 State and prove a result similar to Theorem 4.1.2 that applies to a row diagonally dominant matrix. In particular, show that $\| \ b { A } ^ { - 1 } \| _ { \infty } \leq 1 / \delta$ where $\delta$ measures the strength of the row diagonal dominance as defined in Equation 4.1.3.

P4.1.3 Suppose A is column diagonally dominant, symmetric, and nonsingular and that $A = L D L ^ { T }$ .

What can you say about the size of entries in L and D? Give the smallest upper bound you can for $\parallel L \parallel _ { 1 }$ .

# Notes and References for §4.1

The unsymmetric analog of Algorithm 4.1.2 is related to the methods of Crout and Doolittle. See Stewart (IMC, pp. 131–149) and also:

G.E. Forsythe (1960). “Crout with Pivoting,” Commun. ACM 3, 507–508.

W.M. McKeeman (1962). “Crout with Equilibration and Iteration,” Commun. ACM 5, 553–555.

H.J. Bowdler, R.S. Martin, G. Peters, and J.H. Wilkinson (1966), “Solution of Real and Complex Systems of Linear Equations,” Numer. Math. 8, 217–234.

Just as algorithms can be tailored to exploit structure, so can error analysis and perturbation theory:

C. de Boor and A. Pinkus (1977). “A Backward Error Analysis for Totally Positive Linear Systems,” Numer. Math. 27, 485–490.

J.R. Bunch, J.W. Demmel, and C.F. Van Loan (1989). “The Strong Stability of Algorithms for Solving Symmetric Linear Systems,” SIAM J. Matrix Anal. Applic. 10, 494–499.

A. Barrlund (1991). “Perturbation Bounds for the LDLT and LU Decompositions,” BIT 31, 358–363.

D.J. Higham and N.J. Higham (1992). “Backward Error and Condition of Structured Linear Systems,” SIAM J. Matrix Anal. Applic. 13, 162–175.

J.M. Pe˜na (2004). “LDU Decompositions with L and U Well Conditioned,” ETNA 18, 198–208.

J-G. Sun (2004). “A Note on Backward Errors for Structured Linear Systems,” Numer. Lin. Alg. 12, 585–603.

R. Cant´o, P. Koev, B. Ricarte, and M. Urbano (2008). “LDU Factorization of Nonsingular Totally Positive Matrices,” SIAM J. Matrix Anal. Applic. 30, 777–782.

Numerical issues that associated with the factorization of a diagonaly dominant matrix are discussed in:

J.M. Pe˜na (1998). “Pivoting Strategies Leading to Diagonal Dominance by Rows,” Numer. Math. 81, 293–304.

M. Mendoza, M. Raydan, and P. Tarazaga (1999). “Computing the Nearest Diagonally Dominant Matrix,” Numer. Lin. Alg. 5, 461–474.

A. George and K.D. Ikramov (2005). “Gaussian Elimination Is Stable for the Inverse of a Diagonally Dominant Matrix,” Math. Comput. 73, 653–657.

J.M. Pe˜na (2007). “Strict Diagonal Dominance and Optimal Bounds for the Skeel Condition Number,” SIAM J. Numer. Anal. 45, 1107–1108.

F. Dopico and P. Koev (2011). “Perturbation Theory for the LDU Factorization and Accurate Computations for Diagonally Dominant Matrices,” Numer. Math. 119, 337–371.

# 4.2 Positive Definite Systems

A matrix $A \in \mathbb { R } ^ { n \times n }$ is positive definite if $x ^ { T } A x > 0$ for all nonzero $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , positive semidefinite if $x ^ { T } A x \geq 0$ for all $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , and indefinite if we can find $x , y \in \mathbb { R } ^ { n }$ so $\left( x ^ { T } A x \right) \left( y ^ { T } A y \right) < 0 . \mathrm { ~ S } _  \mathrm  ~ \scriptsize ~  ~$ ymmetric positive definite systems constitute one of the most important classes of special $A x = b$ problems. Consider the 2-by-2 symmetric case. If

$$
A = \left[ \begin{array}{l l} \alpha & \beta \\ \beta & \gamma \end{array} \right]
$$

is positive definite then

$$
x = [ 1, 0 ] ^ {T} \Rightarrow x ^ {T} A x = \alpha > 0,
$$

$$
x = [ 0, 1 ] ^ {T} \quad \Rightarrow x ^ {T} A x = \gamma > 0,
$$

$$
x = [ 1, 1 ] ^ {T} \Rightarrow x ^ {T} A x = \alpha + 2 \beta + \gamma > 0,
$$

$$
x = [ 1, - 1 ] ^ {T} \Rightarrow x ^ {T} A x = \alpha - 2 \beta + \gamma > 0.
$$

The last two equations imply $\left| \beta \right| \le ( \alpha { + } \gamma ) / 2$ . From these results we see that the largest entry in A is on the diagonal and that it is positive. This turns out to be true in general. (See Theorem 4.2.8 below.) A symmetric positive definite matrix has a diagonal that is sufficiently “weighty” to preclude the need for pivoting. A special factorization called the Cholesky factorization is available for such matrices. It exploits both symmetry and definiteness and its implementation is the main focus of this section. However, before those details are pursued we discuss unsymmetric positive definite matrices. This class of matrices is important in its own right and and presents interesting pivot-related issues.

# 4.2.1 Positive Definiteness

Suppose $A \in \mathbb { R } ^ { n \times n }$ is positive definite. It is obvious that a positive definite matrix is nonsingular for otherwise we could find a nonzero x so $x ^ { T } A x = 0$ . However, much more is implied by the positivity of the quadratic form $x ^ { T } A x$ as the following results show.

Theorem 4.2.1. If $A \in \mathbb { R } ^ { n \times n }$ is positive definite and $\boldsymbol { X } \in \mathbb { R } ^ { n \times k }$ has rank k, then $B = X ^ { T } A X \in \mathbb { R } ^ { k \times k }$ is also positive definite.

Proof. If $z \in \mathbb { R } ^ { k }$ satisfies $0 \geq z ^ { T } B z = ( X z ) ^ { T } A ( X z )$ , then $X z = 0$ . But since X has full column rank, this implies that $z = 0 . \qquad \bigtriangledown$

Corollary 4.2.2. If A is positive definite, then all its principal submatrices are positive definite. In particular, all the diagonal entries are positive.

Proof. If v is an integer length-k vector with $1 \leq v _ { 1 } < \cdot \cdot \cdot < v _ { k } \leq n$ , then $X = I _ { n } ( : , v )$ （ is a rank-k matrix made up of columns $v _ { 1 } , \ldots , v _ { k }$ of the identity. It follows from Theorem 4.2.1 that $A ( v , v ) = X ^ { T } A X$ is positive definite.

Theorem 4.2.3. The matrix $A \in \mathbb { R } ^ { n \times n }$ is positive definite if and only if the symmetric matrix

$$
T = \frac {A + A ^ {T}}{2}
$$

has positive eigenvalues.

Proof. Note that $x ^ { T } A x = x ^ { T } T x$ . If $T x = \lambda x$ then $x ^ { T } A x = \lambda \cdot x ^ { T } x$ . Thus, if A is positive definite then λ is positive. Conversely, suppose T has positive eigenvalues and $Q ^ { T } T Q = \mathrm { d i a g } ( \lambda _ { i } )$ is its Schur decomposition. (See §2.1.7.) It follows that if $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $y = Q ^ { T } x$ , then

$$
x ^ {T} A x = x ^ {T} T x = y ^ {T} (Q ^ {T} T Q) y = \sum_ {k = 1} ^ {n} \lambda_ {k} y _ {k} ^ {2} > 0,
$$

completing the proof of the theorem.

□

Corollary 4.2.4. If A is positive definite, then it has an LU factorization and the diagonal entries of U are positive.

Proof. From Corollary 4.2.2, it follows that the submatrices $A ( 1 { : } k , 1 { : } k )$ are nonsingular for $k = 1 { : } n$ and so from Theorem 3.2.1 the factorization $A = L U$ exists. If we apply Theorem 4.2.1 with $X = ( L ^ { - 1 } ) ^ { T } = L ^ { - T }$ , then $B = X ^ { T } A X = L ^ { - 1 } ( L U ) L ^ { - 1 } = U \bar { L } ^ { - \bar { T } }$ is positive definite and therefore has positive diagonal entries. The corollary follows because $L ^ { - T }$ is unit upper triangular and this implies $b _ { i i } = u _ { i i } , ~ i = 1 { : } n$ . □

The mere existence of an LU factorization does not mean that its computation is advisable because the resulting factors may have unacceptably large elements. For example, if $\epsilon > 0$ , then the matrix

$$
A = \left[ \begin{array}{c c} \epsilon & m \\ - m & \epsilon \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ - m / \epsilon & 1 \end{array} \right] \left[ \begin{array}{c c} \epsilon & m \\ 0 & 1 + m ^ {2} / \epsilon \end{array} \right]
$$

is positive definite. However, if $m / \epsilon \gg 1$ , then it appears that some kind of pivoting is in order. This prompts us to pose an interesting question. Are there conditions that guarantee when it is safe to compute the LU-without-pivoting factorization of a positive definite matrix?

# 4.2.2 Unsymmetric Positive Definite Systems

The positive definiteness of a general matrix A is inherited from its symmetric part:

$$
T = \frac {A + A ^ {T}}{2}.
$$

Note that for any square matrix we have $A = T + S$ where

$$
S = \frac {A - A ^ {T}}{2}
$$

is the skew-symmetric part of A. Recall that a matrix S is skew symmetric if $S ^ { T } = - S$ . If S is skew-symmetric, then $x ^ { T } S x = 0$ for all $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $s _ { i i } = 0 , i = 1 { : } n$ . It follows that A is positive definite if and only if its symmetric part is positive definite.

The derivation and analysis of methods for positive definite systems require an understanding about how the symmetric and skew-symmetric parts interact during the LU process.

Theorem 4.2.5. Suppose

$$
A = \left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right] + \left[ \begin{array}{c c} 0 & - w ^ {T} \\ w & C \end{array} \right]
$$

is positive definite and that $B \in \mathbb { R } ^ { ( n - 1 ) \times ( n - 1 ) }$ is symmetric and $C \in \mathbb { R } ^ { ( n - 1 ) \times ( n - 1 ) }$ is skew-symmetric. Then it follows that

$$
A = \left[ \begin{array}{c c} 1 & 0 \\ (v + w) / \alpha & I \end{array} \right] \left[ \begin{array}{c c} \alpha & (v - w) ^ {T} \\ 0 & B _ {1} + C _ {1} \end{array} \right] \tag {4.2.1}
$$

where

$$
B _ {1} = B - \frac {1}{\alpha} \left(v v ^ {T} - w w ^ {T}\right) \tag {4.2.2}
$$

is symmetric positive definite and

$$
C _ {1} = C - \frac {1}{\alpha} \left(w v ^ {T} - v w ^ {T}\right) \tag {4.2.3}
$$

is skew-symmetric.

Proof. Since $\alpha \neq 0$ it follows that (4.2.1) holds. It is obvious from their definitions that $B _ { 1 }$ is symmetric and that $C _ { 1 }$ is skew-symmetric. Thus, all we have to show is that $B _ { 1 }$ is positive definite i.e.,

$$
0 <   z ^ {T} B _ {1} z = z ^ {T} B z - \frac {1}{\alpha} (v ^ {T} z) ^ {2} + \frac {1}{\alpha} (w ^ {T} z) ^ {2} \tag {4.2.4}
$$

for all nonzero $z \in \mathbb { R } ^ { n - 1 }$ . For any $\mu \in \mathbb { R }$ and $0 \neq z \in \mathbb { R } ^ { n - 1 }$ we have

$$
\begin{array}{l} 0 <   \left[ \begin{array}{c} \mu \\ z \end{array} \right] ^ {T} A \left[ \begin{array}{c} \mu \\ z \end{array} \right] = \left[ \begin{array}{c} \mu \\ z \end{array} \right] ^ {T} \left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right] \left[ \begin{array}{c} \mu \\ z \end{array} \right] \\ = \mu^ {2} \alpha + 2 \mu v ^ {T} z + z ^ {T} B z. \\ \end{array}
$$

If $\mu = - ( v ^ { T } z ) / \alpha$ , then

$$
0 <   z ^ {T} B z - \frac {1}{\alpha} (v ^ {T} z) ^ {2},
$$

which establishes the inequality (4.2.4).

From (4.2.1) we see that if $B _ { 1 } + C _ { 1 } = L _ { 1 } U _ { 1 }$ is the LU factorization, then $A = L U$ where

$$
L   =   \left[ \begin{array}{c c} 1 & 0 \\ (v + w) / \alpha & L _ {1} \end{array} \right] \left[ \begin{array}{c c} \alpha & (v - w) ^ {T} \\ 0 & U _ {1} \end{array} \right].
$$

Thus, the theorem shows that triangular factors in $A = L U$ are nicely bounded if S is not too big compared to $T ^ { - 1 }$ . Here is a result that makes this precise:

Theorem 4.2.6. Let $A \in \mathbb { R } ^ { n \times n }$ be positive definite and set $T = ( A + A ^ { T } ) / 2$ and $S = ( A - A ^ { T } ) / 2$ . If A = LU is the LU factorization, then

$$
\left\| | L | | U | \right\| _ {F} \leq n \left(\left\| T \right\| _ {2} + \left\| S T ^ {- 1} S \right\| _ {2}\right). \tag {4.2.5}
$$

Proof. See Golub and Van Loan (1979).

The theorem suggests when it is safe not to pivot. Assume that the computed factors $\hat { L }$ and $\hat { U }$ satisfy

$$
\| | \hat {L} | | \hat {U} | \| _ {F} \leq c \| | L | | U | \| _ {F}, \tag {4.2.6}
$$

where c is a constant of modest size. It follows from (4.2.1) and the analysis in §3.3 that if these factors are used to compute a solution to $A x = b ,$ then the computed solution ˆx satisfies $( { \boldsymbol { A } } + { \boldsymbol { E } } ) { \hat { \boldsymbol { x } } } = { \boldsymbol { b } }$ with

$$
\left\| E \right\| _ {F} \leq \mathbf {u} \left(2 n \| A \| _ {F} + 4 c n ^ {2} \left(\| T \| _ {2} + \| S T ^ {- 1} S \| _ {2}\right)\right) + O \left(\mathbf {u} ^ {2}\right). \tag {4.2.7}
$$

It is easy to show that $\parallel T \parallel _ { 2 } \leq \parallel A \parallel _ { 2 }$ , and so it follows that if

$$
\Omega = \frac {\| S T ^ {- 1} S \| _ {2}}{\| A \| _ {2}} \tag {4.2.8}
$$

is not too large, then it is safe not to pivot. In other words, the norm of the skewsymmetric part S has to be modest relative to the condition of the symmetric part T . Sometimes it is possible to estimate Ω in an application. This is trivially the case when A is symmetric for then $\Omega = 0$ .

# 4.2.3 Symmetric Positive Definite Systems

If we apply the above results to a symmetric positive definite matrix we know that the factorization $A = L U$ exists and is stable to compute. The computation of the factorization $A = L D L ^ { T }$ via Algorithm 4.1.2 is also stable and exploits symmetry. However, for symmetric positive definite systems it is often handier to work with a variation of $L \dot { D } L ^ { T }$ .

Theorem 4.2.7 (Cholesky Factorization). If $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite, then there exists a unique lower triangular $G \in \mathbb { R } ^ { n \times n }$ with positive diagonal entries such that $A = G G ^ { T }$ .

Proof. From Theorem 4.1.3, there exists a unit lower triangular L and a diagonal

$$
D = \mathrm{diag} (d _ {1}, \ldots , d _ {n})
$$

such that $A = L D L ^ { T }$ . Theorem 4.2.1 tells us that $L ^ { - 1 } A L ^ { - T } = D$ is positive definite. Thus, the $d _ { k }$ are positive and the matrix $G = L \operatorname { d i a g } ( { \sqrt { d _ { 1 } } } , \ldots , { \sqrt { d _ { n } } } )$ is real and lower triangular with positive diagonal entries. It also satisfies $A = G G ^ { T }$ . Uniqueness follows from the uniqueness of the $\mathrm { \Delta L D L ^ { T } }$ factorization.

The factorization $A = G G ^ { T }$ is known as the Cholesky factorization and G is the Cholesky factor. Note that if we compute the Cholesky factorization and solve the triangular systems $G y = b$ and $G ^ { T } x = y .$ , then $b = G y = G ( G ^ { T } x ) = ( G G ^ { T } ) x = A x$ .

# 4.2.4 The Cholesky Factor is not a Square Root

A matrix $X \in \mathbb { R } ^ { n \times n }$ that satisfies $A = X ^ { 2 }$ is a square root of A. Note that if A symmetric, positive definite, and not diagonal, then its Cholesky factor is not a square root. However, if $A = G G ^ { T }$ and $\mathbf { \mathit { X } } = { \bar { U } } { \bar { \Sigma } } U ^ { T }$ where $G = U \Sigma V ^ { T }$ is the SVD, then

$$
X ^ {2} = (U \Sigma U ^ {T}) (U \Sigma U ^ {T}) = U \Sigma^ {2} U ^ {T} = (U \Sigma V ^ {T}) (U \Sigma V ^ {T}) ^ {T} = G G ^ {T} = A.
$$

Thus, a symmetric positive definite matrix A has a symmetric positive definite square root denoted by $A ^ { 1 \bar { / } 2 }$ . We have more to say about matrix square roots in §9.4.2.

# 4.2.5 A Gaxpy-Rich Cholesky Factorization

Our proof of the Cholesky factorization in Theorem 4.2.7 is constructive. However, we can develop a more effective procedure by comparing columns in $A = G G ^ { T }$ . If $A \in \mathbb { R } ^ { n \times n }$ and $1 \leq j \leq n$ , then

$$
A (:, j) = \sum_ {k = 1} ^ {j} G (j, k) \cdot G (:, k).
$$

This says that

$$
G (j, j) G (:, j) = A (:, j) - \sum_ {k = 1} ^ {j - 1} G (j, k) \cdot G (:, k) \equiv v. \tag {4.2.9}
$$

If the first $j - 1$ columns of G are known, then v is computable. It follows by equating components in (4.2.9) that

$$
G (j: n, j) = v (j: n) / \sqrt {v (j)}
$$

and so we obtain

for $j = 1:n$ $v(j:n) = A(j:n, j)$ for $k = 1:j - 1$ $v(j:n) = v(j:n) - G(j, k) \cdot G(j:n, k)$ end $G(j:n, j) = v(j:n) / \sqrt{v(j)}$ end

It is possible to arrange the computations so that G overwrites the lower triangle of A.

Algorithm 4.2.1 (Gaxpy Cholesky) Given a symmetric positive definite $A \in \mathbb { R } ^ { n \times n }$ , the following algorithm computes a lower triangular G such that $A = G G ^ { T }$ . For all $i \geq j , G ( i , j )$ overwrites $A ( i , j )$ .

for $j = 1:n$ if $j > 1$ $A(j:n,j) = A(j:n,j) - A(j:n,1:j - 1) \cdot A(j,1:j - 1)^T$ end $A(j:n,j) = A(j:n,j) / \sqrt{A(j,j)}$ end

This algorithm requires $n ^ { 3 } / 3$ flops.

# 4.2.6 Stability of the Cholesky Process

In exact arithmetic, we know that a symmetric positive definite matrix has a Cholesky factorization. Conversely, if the Cholesky process runs to completion with strictly positive square roots, then A is positive definite. Thus, to find out if a matrix A is positive definite, we merely try to compute its Cholesky factorization using any of the methods given above.

The situation in the context of roundoff error is more interesting. The numerical stability of the Cholesky algorithm roughly follows from the inequality

$$
g _ {i j} ^ {2} \leq \sum_ {k = 1} ^ {i} g _ {i k} ^ {2} = a _ {i i}.
$$

This shows that the entries in the Cholesky triangle are nicely bounded. The same conclusion can be reached from the equation $\parallel G \parallel _ { 2 } ^ { 2 } = \parallel A \parallel _ { 2 }$ .

The roundoff errors associated with the Cholesky factorization have been extensively studied in a classical paper by Wilkinson (1968). Using the results in this paper, it can be shown that if ˆx is the computed solution to $A x = b$ , obtained via the Cholesky process, then ˆx solves the perturbed system

$$
(A + E) \hat {x} = b \quad \| E \| _ {2} \leq c _ {n} \mathbf {u} \| A \| _ {2},
$$

where $c _ { n }$ is a small constant that depends upon n. Moreover, Wilkinson shows that if $q _ { n } \mathbf { u } \kappa _ { 2 } ( A ) \leq 1$ where $q _ { n }$ is another small constant, then the Cholesky process runs to completion, i.e., no square roots of negative numbers arise.

It is important to remember that symmetric positive definite linear systems can be ill-conditioned. Indeed, the eigenvalues and singular values of a symmetric positive definite matrix are the same. This follows from (2.4.1) and Theorem 4.2.3. Thus,

$$
\kappa_ {2} (A) = \frac {\lambda_ {\max} (A)}{\lambda_ {\min} (A)}.
$$

The eigenvalue $\lambda _ { \mathrm { m i n } } ( A )$ is the “distance to trouble” in the Cholesky setting. This prompts us to consider a permutation strategy that steers us away from using small diagonal elements that jeopardize the factorization process.

# 4.2.7 The LDLT Factorization with Symmetric Pivoting

With an eye towards handling ill-conditioned symmetric positive definite systems, we return to the $\mathrm { L D L } ^ { T }$ factorization and develop an outer product implementation with pivoting. We first observe that if A is symmetric and $P _ { 1 }$ is a permutation, then $P _ { 1 } A$ is not symmetric. On the other hand, $P _ { 1 } A P _ { 1 } ^ { T }$ is symmetric suggesting that we consider the following factorization:

$$
P _ {1} A P _ {1} ^ {T} = \left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{c c} \alpha & 0 \\ 0 & \tilde {A} \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & I _ {n - 1} \end{array} \right] ^ {T}
$$

where

$$
\tilde {A} = B - \frac {1}{\alpha} v v ^ {T}.
$$

Note that with this kind of symmetric pivoting, the new (1,1) entry α is some diagonal entry $a _ { i i }$ . Our plan is to choose $P _ { 1 }$ so that α is the largest of A’s diagonal entries. If we apply the same strategy recursively to A˜ and compute

$$
\tilde {P} \tilde {A} \tilde {P} ^ {T} = \tilde {L} \tilde {D} \tilde {L} ^ {T},
$$

then we emerge with the factorization

$$
P A P ^ {T} = L D L ^ {T} \tag {4.2.10}
$$

where

$$
P = \left[ \begin{array}{c c} 1 & 0 \\ 0 & \tilde {P} \end{array} \right] P _ {1}, \qquad L = \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & \tilde {L} \end{array} \right], \qquad D = \left[ \begin{array}{c c} \alpha & 0 \\ 0 & \tilde {D} \end{array} \right].
$$

By virtue of this pivot strategy,

$$
d _ {1} \geq d _ {2} \geq \dots \geq d _ {n} > 0.
$$

Here is a nonrecursive implementation of the overall algorithm:

Algorithm 4.2.2 (Outer Product $\mathbf { \mathbf { \mathbf { \mathbf { L } } } } \mathbf { \mathbf { \mathbf { D } } } \mathbf { \mathbf \mathbf { \mathbf { \mathbf { L } } } } ^ { T }$ with Pivoting) Given a symmetric positive semidefinite $A \in \mathbb { R } ^ { n \times n }$ , the following algorithm computes a permutation P , a unit lower triangular L, and a diagonal matrix $D = \mathrm { d i a g } ( d _ { 1 } , \ldots , d _ { n } )$ so $P A P ^ { T } = L D L ^ { T }$ with $d _ { 1 } \geq d _ { 2 } \geq \dots \geq d _ { n } > 0 .$ The matrix element $a _ { i j }$ is overwritten by $d _ { i }$ if $i = j$ and by $\ell _ { i j } { \mathrm { ~ i f ~ } } i > j . P = P _ { 1 } \cdot \cdot \cdot P _ { n } $ where $P _ { k }$ is the identity with rows k and $p i v ( k )$ interchanged.

$$
\begin{array}{l} \text { for } k = 1: n \\ p i v (k) = j \text { where } a _ {j j} = \max \{a _ {k k}, \dots , a _ {n n} \} \\ A (k,:) \leftrightarrow A (j,:) \\ A (:, k) \leftrightarrow A (:, j) \\ \alpha = A (k, k) \\ v = A (k + 1: n, k) \\ A (k + 1: n, k) = v / \alpha \\ A (k + 1: n, k + 1: n) = A (k + 1: n, k + 1: n) - v v ^ {T} / \alpha \\ \end{array}
$$

end

If symmetry is exploited in the outer product update, then $n ^ { 3 } / 3$ flops are required. To solve $A x = b$ given $P A P ^ { T } = L D L ^ { T }$ , we proceed as follows:

$$
L w = P b, \qquad D y = w, \qquad L ^ {T} z = y, \qquad x = P ^ {T} z.
$$

We mention that Algorithm 4.2.2 can be implemented in a way that only references the lower trianglar part of A.

It is reasonable to ask why we even bother with the $\mathrm { L D L } ^ { T }$ factorization given that it appears to offer no real advantage over the Cholesky factorization. There are two reasons. First, it is more efficient in narrow band situations because it avoids square roots; see §4.3.6. Second, it is a graceful way to introduce factorizations of the form

$$
P A P ^ {T} = \left( \begin{array}{c} \mathrm{lower} \\ \mathrm{triangular} \end{array} \right) \times \left( \begin{array}{c} \mathrm{simple} \\ \mathrm{matrix} \end{array} \right) \times \left( \begin{array}{c} \mathrm{lower} \\ \mathrm{triangular} \end{array} \right) ^ {T},
$$

where P is a permutation arising from a symmetry-exploiting pivot strategy. The symmetric indefinite factorizations that we develop in §4.4 fall under this heading as does the “rank revealing” factorization that we are about to discuss for semidefinite problems.

# 4.2.8 The Symmetric Semidefinite Case

A symmetric matrix $A \in \mathbb { R } ^ { n \times n }$ is positive semidefinite if

$$
x ^ {T} A x \geq 0
$$

for every $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . It is easy to show that if $A \in \mathbb { R } ^ { n \times n }$ is symmetric and positive semidefinite, then its eigenvalues satisfy

$$
0 = \lambda_ {n} (A) = \dots = \lambda_ {r + 1} (A) <   \lambda_ {r} (A) \leq \dots \leq \lambda_ {1} (A) \tag {4.2.11}
$$

where r is the rank of A. Our goal is to show that Algorithm 4.2.2 can be used to estimate r and produce a streamlined version of (4.2.10). But first we establish some useful properties.

Theorem 4.2.8. If $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive semidefinite, then

$$
\left| a _ {i j} \right| \leq \left(a _ {i i} + a _ {j j}\right) / 2, \tag {4.2.12}
$$

$$
\left| a _ {i j} \right| \leq \sqrt {a _ {i i} a _ {j j}}, \quad (i \neq j), \tag {4.2.13}
$$

$$
\max \left| a _ {i j} \right| = \max a _ {i i}, \tag {4.2.14}
$$

$$
a _ {i i} = 0 \Rightarrow A (i,:) = 0, A (:, i) = 0. \tag {4.2.15}
$$

Proof. Let $e _ { i }$ denote the ith column of $I _ { n }$ . Since

$$
x = e _ {i} + e _ {j} \Rightarrow 0 \leq x ^ {T} A x = a _ {i i} + 2 a _ {i j} + a _ {j j},
$$

$$
x = e _ {i} - e _ {j} \Rightarrow 0 \leq x ^ {T} A x = a _ {i i} - 2 a _ {i j} + a _ {j j},
$$

it follows that

$$
- 2 a _ {i j} \leq a _ {i i} + a _ {j j},
$$

$$
2 a _ {i j} \leq a _ {i i} + a _ {j j}.
$$

These two equations confirm (4.2.12), which in turn implies (4.2.14).

To prove (4.2.13), set $x = \tau e _ { i } + e _ { j }$ where $\tau \in \mathbb { R }$ . It follows that

$$
0 <   x ^ {T} A x = a _ {i i} \tau^ {2} + 2 a _ {i j} \tau + a _ {j j}
$$

must hold for all τ . This is a quadratic equation in τ and for the inequality to hold, the discriminant $4 a _ { i j } ^ { 2 } - 4 a _ { i i } a _ { j j }$ must be negative, i.e., $| a _ { i j } | \leq \sqrt { a _ { i i } a _ { j j } }$ . The implication in (4.2.15) follows immediately from (4.2.13).

Let us examine what happens when Algorithm 4.2.2 is applied to a rank-r positive semidefinite matrix. If $k \leq r$ , then after k steps we have the factorization

$$
\tilde {P} A \tilde {P} ^ {T} = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & I _ {n - k} \end{array} \right] \left[ \begin{array}{c c} D _ {k} & 0 \\ 0 & A _ {k} \end{array} \right] \left[ \begin{array}{c c} L _ {1 1} ^ {T} & L _ {2 1} ^ {T} \\ 0 & I _ {n - k} \end{array} \right] \tag {4.2.16}
$$

where $D _ { k } = \mathrm { d i a g } ( d _ { 1 } , \ldots , d _ { k } ) \in \mathbb { R } ^ { k \times k }$ and $d _ { 1 } \geq \dots \geq d _ { k } \geq 0$ . By virtue of the pivot strategy, if $d _ { k } = 0$ , then $A _ { k }$ has a zero diagonal. Since $A _ { k }$ is positive semidefinite, it follows from (4.2.15) that $A _ { k } = 0$ . This contradicts the assumption that A has rank r unless $k = r$ . Thus, if $k \leq r .$ , then $d _ { k } > 0$ . Moreover, we must have $A _ { r } = 0$ since A has the same rank as diag $( D _ { r } , A _ { r } )$ . It follows from (4.2.16) that

$$
P A P ^ {T} = \left[ \begin{array}{l} L _ {1 1} \\ L _ {2 1} \end{array} \right] D _ {r} \left[ \begin{array}{l} L _ {1 1} ^ {T} \mid L _ {2 1} ^ {T} \end{array} \right] \tag {4.2.17}
$$

where $D _ { r } = \operatorname { d i a g } ( d _ { 1 } , \ldots , d _ { r } )$ has positive diagonal entries, $L _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ is unit lower triangular, and $L _ { 2 1 } \in \mathbb { R } ^ { ( n - r ) \times r }$ . If $\ell _ { j }$ is the jth column of the L-matrix, then we can rewrite (4.2.17) as a sum of rank-1 matrices:

$$
P A P ^ {T} = \sum_ {j = 1} ^ {r} d _ {j} \ell_ {j} \ell_ {j} ^ {T}.
$$

This can be regarded as a relatively cheap alternative to the SVD rank-1 expansion.

It is important to note that our entire semidefinite discussion has been an exact arithmetic discussion. In practice, a threshold tolerance for small diagonal entries has to be built into Algorithm 4.2.2. If the diagonal of the computed $A _ { k }$ in (4.2.16) is sufficiently small, then the loop can be terminated and $\tilde { r }$ can be regarded as the numerical rank of A. For more details, see Higham (1989).

# 4.2.9 Block Cholesky

Just as there are block methods for computing the LU factorization, so are there are block methods for computing the Cholesky factorization. Paralleling the derivation of the block LU algorithm in §3.2.11, we start by blocking $A = G G ^ { T }$ as follows

$$
\left[ \begin{array}{l l} A _ {1 1} & A _ {2 1} ^ {T} \\ A _ {2 1} & A _ {2 2} \end{array} \right] = \left[ \begin{array}{c c} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right] \left[ \begin{array}{c c} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right] ^ {T}. \tag {4.2.18}
$$

Here, $A _ { 1 1 } \in \mathbb { R } ^ { r \times r } , A _ { 2 2 } \in \mathbb { R } ^ { ( n - r ) \times ( n - r ) }$ , r is a blocking parameter, and G is partitioned conformably. Comparing blocks in (4.2.18) we conclude that

$$
A _ {1 1} = G _ {1 1} G _ {1 1} ^ {T},
$$

$$
A _ {2 1} = G _ {2 1} G _ {1 1} ^ {T},
$$

$$
A _ {2 2} = G _ {2 1} G _ {2 1} ^ {T} + G _ {2 2} G _ {2 2} ^ {T},
$$

which suggests the following 3-step procedure:

Step 1: Compute the Cholesky factorization of $A _ { 1 1 }$ to get $G _ { 1 1 }$ .

Step 2: Solve a lower triangular multiple-right-hand-side system for $G _ { 2 1 }$ .

Step 3: Compute the Cholesky factor $G _ { 2 2 }$ of $A _ { 2 2 } - G _ { 2 1 } G _ { 2 1 } ^ { T } = A _ { 2 2 } - A _ { 2 1 } A _ { 1 1 } ^ { - 1 } A _ { 2 1 } ^ { T }$ . In recursive form we obtain the following algorithm.

Algorithm 4.2.3 (Recursive Block Cholesky) Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite and r is a positive integer. The following algorithm computes a lower triangular $G \in \mathbb { R } ^ { n \times n }$ so $\bar { A } = G G ^ { T }$ .

function $G = { \mathsf { B l o c k C h o l e s k y } } ( A , n , r )$

if $n \leq r$

Compute the Cholesky factorization $A = G G ^ { T }$ .

else

Compute the Cholesky factorization $A ( 1 { : } r , 1 { : } r ) = G _ { 1 1 } G _ { 1 1 } ^ { T }$

Solve $G _ { 2 1 } G _ { 1 1 } ^ { T } = A ( r + 1 { : } n , 1 { : } r )$ for $G _ { 2 1 }$

$$
\tilde {A} = A (r + 1: n, r + 1: n) - G _ {2 1} G _ {2 1} ^ {T}
$$

$$
G _ {2 2} = \text { BlockCholesky } (\tilde {A}, n - r, r)
$$

$$
G = \left[ \begin{array}{c c} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right]
$$

end

If symmetry is exploited in the computation of ${ \tilde { A } } ,$ then this algorithm requires $n ^ { 3 } / 3$ flops. A careful accounting of flops reveals that the level-3 fraction is about $1 - 1 / N ^ { 2 }$ where $N \approx n / r$ . The “small” Cholesky computation for $G _ { 1 1 }$ and the “thin” solution process for $G _ { 2 1 }$ are dominated by the “large” level-3 update for $\tilde { A }$ .

To develop a nonrecursive implementation, we assume for clarity that $n = N r$ where N is a positive integer and consider the partitioning

$$
\left[ \begin{array}{c c c} A _ {1 1} & \dots & A _ {1 N} \\ \vdots & \ddots & \vdots \\ A _ {N 1} & \dots & A _ {N N} \end{array} \right] = \left[ \begin{array}{c c c} G _ {1 1} & \dots & 0 \\ \vdots & \ddots & \vdots \\ G _ {N 1} & \dots & G _ {N N} \end{array} \right] \left[ \begin{array}{c c c} G _ {1 1} & \dots & 0 \\ \vdots & \ddots & \vdots \\ G _ {N 1} & \dots & G _ {N N} \end{array} \right] ^ {T} \tag {4.2.19}
$$

where all blocks are r-by-r. By equating $( i , j )$ blocks with $i \geq j$ it follows that

$$
A _ {i j} = \sum_ {k = 1} ^ {j} G _ {i k} G _ {j k} ^ {T}.
$$

Define

$$
S = A _ {i j} - \sum_ {k = 1} ^ {j - 1} G _ {i k} G _ {j k} ^ {T} = A _ {i j} - \left[ G _ {i 1} \left| \dots \right| G _ {i, j - 1} \right] \left[ \begin{array}{c} G _ {j 1} ^ {T} \\ \vdots \\ G _ {j, j - 1} ^ {T} \end{array} \right].
$$

If $i = j$ , then $G _ { j j }$ is the Cholesky factor of S. If $i > j$ , then $G _ { i j } G _ { j j } ^ { T } = S$ and $G _ { i j }$ is the solution to a triangular multiple right hand side problem. Properly sequenced, these equations can be arranged to compute all the G-blocks.

Algorithm 4.2.4 (Nonrecursive Block Cholesky) Given a symmetric positive definite $A \in \mathbb { R } ^ { n \times n }$ with $n = N r$ with blocking (4.2.19), the following algorithm computes a lower triangular $G \in \mathbb { R } ^ { n \times n }$ such that $A { \stackrel { \cdot } { = } } G G ^ { T }$ . The lower triangular part of A is overwritten by the lower triangular part of G.

for $j = 1:N$ for $i = j:N$ Compute $S = A_{ij} - \sum_{k=1}^{j-1} G_{ik} G_{jk}^T$ .
    if $i = j$ Compute Cholesky factorization $S = G_{jj} G_{jj}^T$ .
    else
    Solve $G_{ij} G_{jj}^T = S$ for $G_{ij}$ .
    end $A_{ij} = G_{ij}$ .
    end
end

The overall process involves $n ^ { 3 } / 3$ flops like the other Cholesky procedures that we have developed. The algorithm is rich in matrix multiplication with a level-3 fraction given by $1 - \left( 1 / N ^ { 2 } \right)$ . The algorithm can be easily modified to handle the case when r does not divide n.

# 4.2.10 Recursive Blocking

It is instructive to look a little more deeply into the implementation of a block Cholesky factorization as it is an occasion to stress the importance of designing data structures that are tailored to the problem at hand. High-performance matrix computations are filled with tensions and tradeoffs. For example, a successful pivot strategy might balance concerns about stability and memory traffic. Another tension is between performance and memory constraints. As an example of this, we consider how to achieve level-3 performance in a Cholesky implementation given that the matrix is represented in packed format. This data structure houses the lower (or upper) triangular portion of a matrix $A \in \mathbb { R } ^ { n \times n }$ in a vector of length $N = n ( n + 1 ) / 2$ . The symvec arrangement stacks the lower triangular subcolumns, e.g.,

$$
\operatorname{symvec} (A) = \left[ a _ {1 1} a _ {2 1} a _ {3 1} a _ {4 1} a _ {2 2} a _ {3 2} a _ {4 2} a _ {3 3} a _ {4 3} a _ {4 4} \right] ^ {T}. \tag {4.2.20}
$$

This layout is not very friendly when it comes to block Cholesky calculations because the assembly of an A-block (say $A ( i _ { 1 } { : } i _ { 2 } , j _ { 1 } { : } j _ { 2 } ) )$ involves irregular memory access patterns. To realize a high-performance matrix multiplication it is usually necessary to have the matrices laid out conventionally as full rectangular arrays that are contiguous in memory, e.g.,

$$
\operatorname{vec} (A) = \left[ a _ {1 1} a _ {2 1} a _ {3 1} a _ {4 1} a _ {1 2} a _ {2 2} a _ {3 2} a _ {4 2} a _ {1 3} a _ {2 3} a _ {3 3} a _ {4 3} a _ {1 4} a _ {2 4} a _ {3 4} a _ {4 4} \right] ^ {T}. \tag {4.2.21}
$$

(Recall that we introduced the vec operation in §1.3.7.) Thus, the challenge is to develop a high performance block algorithm that overwrites a symmetric positive definite A in packed format with its Cholesky factor G in packed format. Toward that end, we present the main ideas behind a recursive data structure that supports level-3 computation and is storage efficient. As memory hierarchies get deeper and more complex, recursive data structures are an interesting way to address the problem of blocking for performance.

The starting point is once again a 2-by-2 blocking of the equation $A = G G ^ { T }$ :

$$
{\left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right]} = {\left[ \begin{array}{l l} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right]} {\left[ \begin{array}{l l} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right]} ^ {T}.
$$

However, unlike in (4.2.18) where $A _ { 1 1 }$ has a chosen block size, we now assume that $A _ { 1 1 } \in \mathbb { R } ^ { m \times m }$ where $m = \operatorname { c e i l } ( n / 2 )$ . In other words, the four blocks are roughly the same size. As before, we equate entries and identify the key subcomputations:

<table><tr><td> $G_{11}G_{11}^{T} = A_{11}$ </td><td>half-sized Cholesky.</td></tr><tr><td> $G_{21}G_{11}^{T} = A_{21}$ </td><td>multiple-right-hand-side triangular solve.</td></tr><tr><td> $\tilde{A}_{22} = A_{22} - G_{21}G_{21}^{T}$ </td><td>symmetric matrix multiplication update.</td></tr><tr><td> $G_{22}G_{22}^{T} = \tilde{A}_{22}$ </td><td>half-sized Cholesky.</td></tr></table>

Our goal is to develop a symmetry-exploiting, level-3-rich procedure that overwrites A with its Cholesky factor G. To do this we introduce the mixed packed format. An n = 9 example with A11 ∈ IR5×5 $n = 9$ $A _ { 1 1 } \in \mathbb { R } ^ { 5 \times 5 }$ serves to distinguish this layout from the conventional packed format layout:

<table><tr><td>1</td><td></td><td></td><td></td><td></td></tr><tr><td>2 10</td><td></td><td></td><td></td><td></td></tr><tr><td>3 11 18</td><td></td><td></td><td></td><td></td></tr><tr><td>4 12 19 25</td><td></td><td></td><td></td><td></td></tr><tr><td>5 13 20 26 31</td><td></td><td></td><td></td><td></td></tr><tr><td>6 14 21 27 32</td><td></td><td></td><td></td><td>36</td></tr><tr><td>7 15 22 28 33</td><td></td><td></td><td></td><td>37 40</td></tr><tr><td>8 16 23 29 34</td><td></td><td></td><td></td><td>38 41 43</td></tr><tr><td>9 17 24 30 35</td><td></td><td></td><td></td><td>39 42 44 45</td></tr></table>

Packed format

<table><tr><td>1</td><td></td><td></td><td></td><td></td><td></td></tr><tr><td>2</td><td>6</td><td></td><td></td><td></td><td></td></tr><tr><td>3</td><td> $\overline{t}$ </td><td>10</td><td></td><td></td><td></td></tr><tr><td>4</td><td>8</td><td>11</td><td>13</td><td></td><td></td></tr><tr><td>5</td><td>9</td><td>12</td><td>14</td><td>15</td><td></td></tr><tr><td>16</td><td>20</td><td>24</td><td>28</td><td>32</td><td>36</td></tr><tr><td>17</td><td>21</td><td>25</td><td>29</td><td>33</td><td>37 40</td></tr><tr><td>18</td><td>22</td><td>26</td><td>30</td><td>34</td><td>38 41 43</td></tr><tr><td>19</td><td>23</td><td>27</td><td>31</td><td>35</td><td>39 42 44 45</td></tr></table>

Mixed packed format

Notice how the entries from $A _ { 1 1 }$ and $A _ { 2 1 }$ are shuffled with the conventional packed format layout. On the other hand, with the mixed packed format layout, the 15 entries that define $A _ { 1 1 }$ are followed by the 20 numbers that define $A _ { 2 1 }$ which in turn are followed by the 10 numbers that define $A _ { 2 2 }$ . The process can be repeated on $A _ { 1 1 }$ and

<table><tr><td colspan="3">12 43 5 6</td><td colspan="2"></td><td rowspan="2" colspan="2"></td></tr><tr><td colspan="3">7 9 118 10 12</td><td colspan="2">1314 15</td></tr><tr><td rowspan="3" colspan="3">16 20 2417 21 2518 22 2619 23 27</td><td rowspan="3" colspan="2">28 3229 3330 3431 35</td><td colspan="2">3637 38</td></tr><tr><td rowspan="2" colspan="2">39 4140 42</td></tr><tr></tr></table>

Thus, the key to this recursively defined data layout is the idea of representing square diagonal blocks in a mixed packed format. To be precise, recall the definition of vec and symvec in (4.2.20) and (4.2.21). If $C \in \mathbb { R } ^ { q \times q }$ is such a block, then

$$
\operatorname{mixvec} (C) = \left[ \begin{array}{c} \operatorname{symvec} \left(C _ {1 1}\right) \\ \operatorname{vec} \left(C _ {2 1}\right) \\ \operatorname{symvec} \left(C _ {2 2}\right) \end{array} \right] \tag {4.2.22}
$$

where $m = \mathrm { c e i l } ( q / 2 ) , C _ { 1 1 } = C ( 1 { : } m , 1 { : } m ) , C _ { 2 2 } = C ( m + 1 { : } n , m + 1 { : } n )$ , and $C _ { 2 1 } =$ $C ( m + 1 { : } n , 1 { : } m )$ . Notice that since $C _ { 2 1 }$ is conventionally stored, it is ready to be engaged in a high-performance matrix multiplication.

We now outline a recursive, divide-and-conquer block Cholesky procedure that works with A in packed format. To achieve high performance the incoming A is converted to mixed format at each level of the recursion. Assuming the existence of a triangular system solve procedure TriSol (for the system $G _ { 2 1 } G _ { 1 1 } ^ { T } = A _ { 2 1 } )$ and a symmetric update procedure SymUpdate (for $A _ { 2 2 }  A _ { 2 2 } - G _ { 2 1 } G _ { 2 1 } ^ { T } )$ we have the following framework:

function G = PackedBlockCholesky(A)

{A and G in packed format}

$$
n = \operatorname{size} (A)
$$

if $n \leq n$ min

G is obtained via any level-2, packed-format Cholesky method .

else

Set $m = { \mathsf { c e i l } } ( n / 2 )$ and overwrite A’s packed-format representation with its mixed-format representation.

$$
G _ {1 1} = \text { PackedBlockCholesky } (A _ {1 1})
$$

$$
G _ {2 1} = \operatorname{TriSol} \left(G _ {1 1}, A _ {2 1}\right)
$$

$$
A _ {2 2} = \text { SymUpdate } (A _ {2 2}, G _ {2 1})
$$

$$
G _ {2 2} = \text { PackedBlockCholesky } (A _ {2 2})
$$

end

Here, $n _ { \mathrm { m i n } }$ is a threshold dimension below which it is not possible to achieve level-3 performance. To take full advantage of the mixed format, the procedures TriSol and SymUpdate require a recursive design based on blockings that halve problem size. For example, TriSol should take the incoming packed format $A _ { 1 1 }$ , convert it to mixed format, and solve a 2-by-2 blocked system of the form

$$
\left[ \begin{array}{c c} X _ {1} & X _ {2} \end{array} \right] \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \end{array} \right] ^ {T} = \left[ \begin{array}{c c} B _ {1} & B _ {2} \end{array} \right].
$$

This sets up a recursive solution based on the half-sized problems

$$
\begin{array}{l} X _ {1} L _ {1 1} ^ {T} = B _ {1}, \\ X _ {2} L _ {2 2} ^ {T} = B _ {2} - X _ {1} L _ {2 1} ^ {T}. \\ \end{array}
$$

Likewise, SymUpdate should take the incoming packed format $A _ { 2 2 }$ , convert it to mixed format, and block the required update as follows:

$$
\left[ \begin{array}{c c} C _ {1 1} & C _ {2 1} ^ {T} \\ C _ {2 1} & C _ {2 2} \end{array} \right] = \left[ \begin{array}{c c} C _ {1 1} & C _ {2 1} ^ {T} \\ C _ {2 1} & C _ {2 2} \end{array} \right] - \left[ \begin{array}{c} Y _ {1} \\ Y _ {2} \end{array} \right] \left[ \begin{array}{c} Y _ {1} \\ Y _ {2} \end{array} \right] ^ {T}.
$$

The evaluation is recursive and based on the half-sized updates

$$
\begin{array}{l} C _ {1 1} = C _ {1 1} - Y _ {1} Y _ {1} ^ {T}, \\ C _ {2 1} = C _ {2 1} - Y _ {2} Y _ {1} ^ {T}, \\ C _ {2 2} = C _ {2 2} - Y _ {2} Y _ {2} ^ {T}. \\ \end{array}
$$

Of course, if the incoming matrices are small enough relative to $n _ { \mathrm { m i n } }$ , then TriSol and SymUpdate carry out their tasks conventionally without any further subdivisions.

Overall, it can be shown that PackedBlockCholesky has a level-3 fraction approximately equal to $1 - O ( n _ { \operatorname* { m i n } } / n )$ .

# Problems

P4.2.1 Suppose that $H = A + i B$ is Hermitian and positive definite with A, $B \in \mathbb { R } ^ { n \times n }$ . This means that $x ^ { H } H x > 0$ whenever x 
= 0. (a) Show that

$$
C = \left[ \begin{array}{c c} A & - B \\ B & A \end{array} \right]
$$

is symmetric and positive definite. (b) Formulate an algorithm for solving $( A + i B ) ( x + i y ) = ( b + i c )$ , where b, c, x, and y are in $\mathbb { R } ^ { n }$ . It should involve $8 n ^ { 3 } / 3$ flops. How much storage is required?

P4.2.2 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and positive definite. Give an algorithm for computing an upper triangular matrix $R \in \mathbb { R } ^ { n \times n }$ such that $\dot { A } = R R ^ { T }$ .

P4.2.3 Let $A \in \mathbb { R } ^ { n \times n }$ be positive definite and set $T = ( A + A ^ { T } ) / 2$ and $S = ( A - A ^ { T } ) / 2$ . (a) Show that $\parallel A ^ { - 1 } \parallel _ { 2 } \leq \parallel T ^ { - 1 } \parallel _ { 2 }$ and $x ^ { T } A ^ { - 1 } x \leq x ^ { T } T ^ { - 1 } x$ for all $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . (b) Show that if $A = L D M ^ { T }$ , then $d _ { k } \geq 1 / \Vert \ T ^ { - 1 } \Vert _ { 2 }$ for $k = 1 { : } n$ .

P4.2.4 Find a 2-by-2 real matrix A with the property that $x ^ { T } A x > 0$ for all real nonzero 2-vectors but which is not positive definite when regarded as a member of C2×2 . $\mathbb { C } ^ { 2 \times 2 }$

P4.2.5 Suppose $A \in \mathbb { R } ^ { n \times n }$ has a positive diagonal. Show that if both A and $A ^ { T }$ are strictly diagonally

dominant, then A is positive definite.

P4.2.6 Show that the function $f ( x ) = \sqrt { x ^ { T } A x } / 2$ is a vector norm on $\mathbb { R } ^ { n }$ if and only if A is positive definite.

P4.2.7 Modify Algorithm 4.2.1 so that if the square root of a negative number is encountered, then the algorithm finds a unit vector x so that $x ^ { T } A x < 0$ and terminates.

P4.2.8 Develop an outer product implementation of Algorithm 4.2.1 and a gaxpy implementation of Algorithm 4.2.2.

P4.2.9 Assume that $A \in \mathbb { C } ^ { n \times n }$ is Hermitian and positive definite. Show that if $a _ { 1 1 } = \cdots = a _ { n n } = 1$ and $| a _ { i j } | < 1$ for all $i \neq j$ , then $\mathrm { d i a g } ( A ^ { - 1 } ) \geq \mathrm { d i a g } ( ( \mathsf { R e } ( A ) ) ^ { - 1 } )$ .

P4.2.10 Suppose $A = I + u u ^ { T }$ where $A \in \mathbb { R } ^ { n \times n }$ and $\parallel u \parallel _ { 2 } = 1$ . Give explicit formulae for the diagonal and subdiagonal of A’s Cholesky factor.

P4.2.11 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite and that its Cholesky factor is available. Let $e _ { k } = I _ { n } ( : , k )$ . For $1 \leq i < j \leq n$ , let $\alpha _ { i j }$ be the smallest real that makes $\mathbf { \bar { \Phi } } A + \alpha ( e _ { i } e _ { j } ^ { T } + e _ { j } e _ { i } ^ { T } )$ singular. Likewise, let $\alpha _ { i i }$ be the smallest real that makes $( A + \alpha e _ { i } e _ { i } ^ { T } )$ singular. Show how to compute these quantities using the Sherman-Morrison-Woodbury formula. How many flops are required to find all the $\alpha _ { i j } ?$

P4.2.12 Show that if

$$
M = \left[ \begin{array}{c c} A & B \\ B ^ {T} & C \end{array} \right]
$$

is symmetric positive definite and A and C are square, then

$$
M ^ {- 1} = \left[ \begin{array}{c c} A ^ {- 1} + A ^ {- 1} B S ^ {- 1} B ^ {T} A ^ {- 1} & - A ^ {- 1} B S ^ {- 1} \\ S ^ {- 1} B ^ {T} A ^ {- 1} & S ^ {- 1} \end{array} \right], \qquad S = C - B ^ {T} A ^ {- 1} B.
$$

P4.2.13 Suppose $\sigma \in \mathbb { R }$ and $u \in \mathbb { R } ^ { n }$ . Under what conditions can we find a matrix $X \in \mathbb { R } ^ { n \times n }$ so that $X ( I + \sigma u u ^ { T } ) X = I _ { n } ?$ Give an efficient algorithm for computing X if it exists.

P4.2.14 Suppose $D = \operatorname { d i a g } ( d _ { 1 } , \dotsc , d _ { n } )$ with $d _ { i } > 0$ for all i. Give an efficient algorithm for computing the largest entry in the matrix $( D + C C ^ { T } ) ^ { - 1 }$ where $C \in \mathbb { R } ^ { n \times r }$ . Hint: Use the Sherman-Morrison-Woodbury formula.

P4.2.15 Suppose $A ( \lambda )$ has continuously differentiable entries and is always symmetric and positive definite. If $f ( \lambda ) = \log ( { \mathsf { d e t } } ( A ( \lambda ) ) )$ , then how would you compute $f ^ { \prime } ( 0 ) \smash { \ ? }$

P4.2.16 Suppose $A \in \mathbb { R } ^ { n \times n }$ is a rank-r symmetric positive semidefinite matrix. Assume that it costs one dollar to evaluate each $\boldsymbol { a } _ { i j }$ . Show how to compute the factorization (4.2.17) spending only $O ( n r )$ dollars on $a _ { i j }$ evaluation.

P4.2.17 The point of this problem is to show that from the complexity point of view, if you have a fast matrix multiplication algorithm, then you have an equally fast matrix inversion algorithm, and vice versa. (a) Suppose $F _ { n }$ is the number of flops required by some method to form the inverse of an n-by-n matrix. Assume that there exists a constant $c _ { 1 }$ and a real number α such that $F _ { n } \leq c _ { 1 } n ^ { \alpha }$ for all n. Show that there is a method that can compute the $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ matrix product $A B$ with fewer than $c _ { 2 } n ^ { \alpha }$ flops where $c _ { 2 }$ is a constant independent of n. Hint: Consider the inverse of

$$
C = \left[ \begin{array}{c c c} I _ {n} & A & 0 \\ 0 & I _ {n} & B \\ 0 & 0 & I _ {n} \end{array} \right].
$$

(b) Let $G _ { n }$ be the number of flops required by some method to form the $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ matrix product $A B$ . Assume that there exists a constant $c _ { 1 }$ and a real number α such that $G _ { n } \leq c _ { 1 } n ^ { \alpha }$ for all n. Show that there is a method that can invert a nonsingular n-by-n matrix A with fewer than $c _ { 2 } n ^ { \alpha }$ flops where $c _ { 2 }$ is a constant. Hint: First show that the result applies for triangular matrices by applying recursion to

$$
\left[ \begin{array}{c c} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right] ^ {- 1} = \left[ \begin{array}{c c} G _ {1 1} ^ {- 1} & 0 \\ - G _ {2 2} ^ {- 1} G _ {2 1} G _ {1 1} ^ {- 1} & G _ {2 2} ^ {- 1} \end{array} \right].
$$

Then observe that for general A, $A ^ { - 1 } = A ^ { T } ( A A ^ { T } ) ^ { - 1 } = A ^ { T } G ^ { - T } G ^ { - 1 }$ where $A A ^ { T } = G G ^ { T }$ is the Cholesky factorization.

# Notes and References for §4.2

For an in-depth theoretical treatment of positive definiteness, see:

R. Bhatia (2007). Positive Definite Matrices, Princeton University Press, Princeton, NJ.

The definiteness of the quadratic form $x ^ { T }$ Ax can frequently be established by considering the mathematics of the underlying problem. For example, the discretization of certain partial differential operators gives rise to provably positive definite matrices. Aspects of the unsymmetric positive definite problem are discussed in:

A. Buckley (1974). “A Note on Matrices A = I + H, H Skew-Symmetric,” Z. Angew. Math. Mech. 54, 125–126.   
A. Buckley (1977). “On the Solution of Certain Skew-Symmetric Linear Systems,” SIAM J. Numer. Anal. 14, 566–570.   
G.H. Golub and C. Van Loan (1979). “Unsymmetric Positive Definite Linear Systems,” Lin. Alg. Applic. 28, 85–98.   
R. Mathias (1992). “Matrices with Positive Definite Hermitian Part: Inequalities and Linear Systems,” SIAM J. Matrix Anal. Applic. 13, 640–654.   
K.D. Ikramov and A.B. Kucherov (2000). “Bounding the growth factor in Gaussian elimination for Buckley’s class of complex symmetric matrices,” Numer. Lin. Alg. 7, 269–274.

Complex symmetric matrices have the property that their real and imaginary parts are each symmetric. The following paper shows that if they are also positive definite, then the $\mathrm { \dot { L } D L ^ { T } }$ factorization is safe to compute without pivoting:

S. Serbin (1980). “On Factoring a Class of Complex Symmetric Matrices Without Pivoting,” Math. Comput. 35, 1231–1234.

Historically important Algol implementations of the Cholesky factorization include:

R.S. Martin, G. Peters, and J.H. Wilkinson (1965). “Symmetric Decomposition of a Positive Definite Matrix,” Numer. Math. 7, 362–83.   
R.S. Martin, G. Peters, and J.H. Wilkinson (1966). “Iterative Refinement of the Solution of a Positive Definite System of Equations,” Numer. Math. 8, 203–16.   
F.L. Bauer and C. Reinsch (1971). “Inversion of Positive Definite Matrices by the Gauss-Jordan Method,” in Handbook for Automatic Computation Vol. 2, Linear Algebra, J.H. Wilkinson and C. Reinsch (eds.), Springer-Verlag, New York, 45–49.

For roundoff error analysis of Cholesky, see:

J.H. Wilkinson (1968). “A Priori Error Analysis of Algebraic Processes,” Proceedings of the International Congress on Mathematics, Izdat. Mir, 1968, Moscow, 629–39.   
J. Meinguet (1983). “Refined Error Analyses of Cholesky Factorization,” SIAM J. Numer. Anal. 20, 1243–1250.   
A. Kielbasinski (1987). “A Note on Rounding Error Analysis of Cholesky Factorization,” Lin. Alg. Applic. 88/89, 487–494.   
N.J. Higham (1990). “Analysis of the Cholesky Decomposition of a Semidefinite Matrix,” in Reliable Numerical Computation, M.G. Cox and S.J. Hammarling (eds.), Oxford University Press, Oxford, U.K., 161–185.   
J-Guang Sun (1992). “Rounding Error and Perturbation Bounds for the Cholesky and $L D L ^ { T }$ Factorizations,” Lin. Alg. Applic. 173, 77–97.

The floating point determination of positive definiteness is an interesting problem, see:

S.M. Rump (2006). “Verification of Positive Definiteness,” BIT 46, 433–452.

The question of how the Cholesky triangle G changes when $A = G G ^ { T }$ is perturbed is analyzed in:

G.W. Stewart (1977). “Perturbation Bounds for the QR Factorization of a Matrix,” SIAM J. Num. Anal. 14, 509–18.

Z. Dram˘ac, M. Omladi˘c, and K. Veseli˘c (1994). “On the Perturbation of the Cholesky Factorization,” SIAM J. Matrix Anal. Applic. 15, 1319–1332.

X-W. Chang, C.C. Paige, and G.W. Stewart (1996). “New Perturbation Analyses for the Cholesky Factorization,” IMA J. Numer. Anal. 16, 457–484.


---

<!-- golub_200_249 -->

G.W. Stewart (1997) “On the Perturbation of LU and Cholesky Factors,” IMA J. Numer. Anal. 17, 1–6.   
Nearness/sensitivity issues associated with positive semidefiniteness are presented in:   
N.J. Higham (1988). “Computing a Nearest Symmetric Positive Semidefinite Matrix,” Lin. Alg. Applic. 103, 103–118.   
The numerical issues associated with semi-definite rank determination are covered in:   
P.C. Hansen and P.Y. Yalamov (2001). “Computing Symmetric Rank-Revealing Decompositions via Triangular Factorization,” SIAM J. Matrix Anal. Applic. 23, 443–458.   
M. Gu and L. Miranian (2004). “Strong Rank-Revealing Cholesky Factorization,” ETNA 17, 76–92.   
The issues that surround level-3 performance of packed-format Cholesky are discussed in:   
F.G. Gustavson (1997). “Recursion Leads to Automatic Variable Blocking for Dense Linear-Algebra Algorithms,” IBM J. Res. Dev. 41, 737–756.   
F.G. Gustavson, A. Henriksson, I. Jonsson, B. K˚agstr¨om, , and P. Ling (1998). “Recursive Blocked Data Formats and BLAS’s for Dense Linear Algebra Algorithms,” Applied Parallel Computing Large Scale Scientific and Industrial Problems, Lecture Notes in Computer Science, Springer-Verlag, 1541/1998, 195–206.   
F.G. Gustavson and I. Jonsson (2000). “Minimal Storage High-Performance Cholesky Factorization via Blocking and Recursion,” IBM J. Res. Dev. 44, 823–849.   
B.S. Andersen, J. Wasniewski, and F.G. Gustavson (2001). “A Recursive Formulation of Cholesky Factorization of a Matrix in Packed Storage,” ACM Trans. Math. Softw. 27, 214–244.   
E. Elmroth, F. Gustavson, I. Jonsson, and B. K˚agstr¨om, (2004). “Recursive Blocked Algorithms and Hybrid Data Structures for Dense Matrix Library Software,” SIAM Review 46, 3–45.   
F.G. Gustavson, J. Wasniewski, J.J. Dongarra, and J. Langou (2010). “Rectangular Full Packed Format for Cholesky’s Algorithm: Factorization, Solution, and Inversion,” ACM Trans. Math. Softw. 37, Article 19.   
Other high-performance Cholesky implementations include:   
F.G. Gustavson, L. Karlsson, and B. K˚agstr¨om, (2009). “Distributed SBP Cholesky Factorization Algorithms with Near-Optimal Scheduling,” ACM Trans. Math. Softw. 36, Article 11.   
G. Ballard, J. Demmel, O. Holtz, and O. Schwartz (2010). “Communication-Optimal Parallel and Sequential Cholesky,” SIAM J. Sci. Comput. 32, 3495–3523.   
P. Bientinesi, B. Gunter, and R.A. van de Geijn (2008). “Families of Algorithms Related to the Inversion of a Symmetric Positive Definite Matrix,” ACM Trans. Math. Softw. 35, Article 3.   
M.D. Petkovi´c and P.S. Stanimirovi´c (2009). “Generalized Matrix Inversion is not Harder than Matrix Multiplication,” J. Comput. Appl. Math. 230, 270–282.

# 4.3 Banded Systems

In many applications that involve linear systems, the matrix of coefficients is banded. This is the case whenever the equations can be ordered so that each unknown $x _ { i }$ appears in only a few equations in a “neighborhood”of the ith equation. Recall from §1.2.1 that $A = \left( a _ { i j } \right)$ has upper bandwidth q if $a _ { i j } = 0$ whenever $j > i + q$ and lower bandwidth p if $a _ { i j } = 0$ whenever $i > j + p$ . Substantial economies can be realized when solving banded systems because the triangular factors in LU, $\mathrm { G G } ^ { T }$ , and $\mathrm { L D L } ^ { T }$ are also banded.

# 4.3.1 Band LU Factorization

Our first result shows that if A is banded and $A = L U$ , then L inherits the lower bandwidth of A and U inherits the upper bandwidth of A.

Theorem 4.3.1. Suppose $A \in \mathbb { R } ^ { n \times n }$ has an LU factorization $A = L U$ . If A has upper bandwidth q and lower bandwidth p, then U has upper bandwidth q and L has lower bandwidth p.

Proof. The proof is by induction on n. Since

$$
A   =   \left[ \begin{array}{c c} \alpha & w ^ {T} \\ v & B \end{array} \right]   =   \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ 0 & B - v w ^ {T} / \alpha \end{array} \right] \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & I _ {n - 1} \end{array} \right]  .
$$

It is clear that $B - v w ^ { T } / \alpha$ has upper bandwidth q and lower bandwidth p because only the first q components of w and the first p components of v are nonzero. Let $L _ { 1 } U _ { 1 }$ be the LU factorization of this matrix. Using the induction hypothesis and the sparsity of w and v, it follows that

$$
L = \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & L _ {1} \end{array} \right], \qquad U = \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & U _ {1} \end{array} \right]
$$

have the desired bandwidth properties and satisfy $A = L U . \quad \ \bigtriangledown$

The specialization of Gaussian elimination to banded matrices having an LU factorization is straightforward.

Algorithm 4.3.1 (Band Gaussian Elimination) Given $A \in \mathbb { R } ^ { n \times n }$ with upper bandwidth q and lower bandwidth $p ,$ the following algorithm computes the factorization $A = L U$ , assuming it exists. $A ( i , j )$ is overwritten by $L ( i , j ) { \mathrm { ~ i f ~ } } i > j$ and by $U ( i , j )$ otherwise.

for k = 1:n - 1
    for i = k + 1:min{k + p, n}
    A(i, k) = A(i, k)/A(k, k)
    end
    for j = k + 1:min{k + q, n}
    for i = k + 1:min{k + p, n}
    A(i, j) = A(i, j) - A(i, k)·A(k, j)
    end
    end
end

If $n \gg p$ and $n \gg q ,$ , then this algorithm involves about $2 n p q$ flops. Effective implementations would involve band matrix data structures; see §1.2.5. A band version of Algorithm 4.1.1 (LDLT) is similar and we leave the details to the reader.

# 4.3.2 Band Triangular System Solving

Banded triangular system solving is also fast. Here are the banded analogues of Algorithms 3.1.3 and 3.1.4:

Algorithm 4.3.2 (Band Forward Substitution) Let $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ be a unit lower triangular matrix with lower bandwidth $p .$ Given $b \in \mathbb { R } ^ { n }$ , the following algorithm overwrites b with the solution to $L x = b$ .

for $j = 1:n$ for $i = j + 1:\min \{j + p,n\}$ $b(i) = b(i) - L(i,j)\cdot b(j)$ end   
end

If $n \gg p ,$ then this algorithm requires about 2np flops.

Algorithm 4.3.3 (Band Back Substitution) Let $U \in \mathbb { R } ^ { n \times n }$ be a nonsingular upper triangular matrix with upper bandwidth q. Given $b \in \mathbb { R } ^ { n }$ , the following algorithm overwrites b with the solution to $U x = b$ .

for $j = n: -1:1$ $b(j) = b(j)/U(j, j)$ for $i = \max\{1, j - q\}: j - 1$ $b(i) = b(i) - U(i, j) \cdot b(j)$ end
end

If $n \gg q$ , then this algorithm requires about 2nq flops.

# 4.3.3 Band Gaussian Elimination with Pivoting

Gaussian elimination with partial pivoting can also be specialized to exploit band structure in A. However, if $P A = L U$ , then the band properties of L and U are not quite so simple. For example, if A is tridiagonal and the first two rows are interchanged at the very first step of the algorithm, then $u _ { 1 3 }$ is nonzero. Consequently, row interchanges expand bandwidth. Precisely how the band enlarges is the subject of the following theorem.

Theorem 4.3.2. Suppose $A \in \mathbb { R } ^ { n \times n }$ is nonsingular and has upper and lower bandwidths q and $p ,$ respectively. If Gaussian elimination with partial pivoting is used to compute Gauss transformations

$$
M _ {j} = I - \alpha^ {(j)} e _ {j} ^ {T} \qquad j = 1: n - 1
$$

and permutations $P _ { 1 } , \ldots , P _ { n - 1 }$ such that $M _ { n - 1 } P _ { n - 1 } \cdot \cdot \cdot M _ { 1 } P _ { 1 } A = U$ is upper triangular, then U has upper bandwidth $p + q$ and $\alpha _ { i } ^ { ( j ) } = \theta$ whenever $i \leq j$ or $i > j + p$ .

Proof. Let $P A = L U$ be the factorization computed by Gaussian elimination with partial pivoting and recall that $P = P _ { n - 1 } \cdot \cdot \cdot P _ { 1 }$ . Write $P ^ { T } = \ [ e _ { s _ { 1 } } \ | \cdot \cdot \cdot | \ e _ { s _ { n } } \ ]$ , where $\{ s _ { 1 } , . . . , s _ { n } \}$ is a permutation of $\{ 1 , 2 , . . . , n \}$ . If $s _ { i } > i + p$ then it follows that the leading i-by-i principal submatrix of P A is singular, since $[ P A ] _ { i j } = a _ { s _ { i } , j }$ for $j = 1 \colon s _ { i } - p - 1$ and $s _ { i } - p - 1 \geq i$ . This implies that U and A are singular, a contradiction. Thus, $s _ { i } \leq i + p$ for $i = 1 { : } n$ and therefore, P A has upper bandwidth $p + q$ . It follows from Theorem 4.3.1 that U has upper bandwidth $p + q .$ . The assertion about the $\alpha ^ { ( j ) }$ can be verified by observing that $M _ { j }$ need only zero elements $( j + 1 , j ) , \dotsc , ( j + p , j )$ of the partially reduced matrix $P _ { j } M _ { j - 1 } P _ { j - 1 } \cdot \cdot \cdot _ { 1 } P _ { 1 } A$ .

Thus, pivoting destroys band structure in the sense that U becomes “wider” than A’s upper triangle, while nothing at all can be said about the bandwidth of L. However, since the jth column of L is a permutation of the jth Gauss vector $\alpha _ { j }$ , it follows that L has at most $p + 1$ nonzero elements per column.

# 4.3.4 Hessenberg LU

As an example of an unsymmetric band matrix computation, we show how Gaussian elimination with partial pivoting can be applied to factor an upper Hessenberg matrix H. (Recall that if H is upper Hessenberg then $h _ { i j } = 0 , i > j + 1 . )$ After $k - 1$ steps of Gaussian elimination with partial pivoting we are left with an upper Hessenberg matrix of the form

$$
\left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right], \qquad k = 3, n = 5.
$$

By virtue of the special structure of this matrix, we see that the next permutation, $P _ { 3 } .$ , is either the identity or the identity with rows 3 and 4 interchanged. Moreover, the next Gauss transformation $M _ { k }$ has a single nonzero multiplier in the $( k + 1 , k )$ position. This illustrates the kth step of the following algorithm.

Algorithm 4.3.4 (Hessenberg LU) Given an upper Hessenberg matrix $H \in \mathbb { R } ^ { n \times n }$ , the following algorithm computes the upper triangular matrix $M _ { n - 1 } P _ { n - 1 } \cdot \cdot \cdot M _ { 1 } P _ { 1 } H = U$ where each $P _ { k }$ is a permutation and each $M _ { k }$ is a Gauss transformation whose entries are bounded by unity. $H ( i , k )$ is overwritten with $U ( i , k )$ if $i \leq k$ and by $- [ M _ { k } ] _ { k + 1 , k }$ if $i = k + 1$ . An integer vector $p i v ( 1 { : } n - 1 )$ encodes the permutations. If $P _ { k } = I$ , then $p i v ( k ) = 0$ . If $P _ { k }$ interchanges rows k and $k + 1$ , then $p i v ( k ) = 1$ .

for k = 1:n - 1
    if $|H(k,k)| < |H(k+1,k)|$ $piv(k) = 1; H(k,k:n) \leftrightarrow H(k+1,k:n)$ else $piv(k) = 0$ end
    if $H(k,k) \neq 0$ $\tau = H(k+1,k)/H(k,k)$ $H(k+1,k+1:n) = H(k+1,k+1:n) - \tau \cdot H(k,k+1:n)$ $H(k+1,k) = \tau$ end
end

This algorithm requires $n ^ { 2 }$ flops.

# 4.3.5 Band Cholesky

The rest of this section is devoted to banded $A x = b$ problems where the matrix A is also symmetric positive definite. The fact that pivoting is unnecessary for such matrices leads to some very compact, elegant algorithms. In particular, it follows from Theorem 4.3.1 that if $A = G G ^ { T }$ is the Cholesky factorization of $A ,$ , then G has the same lower bandwidth as A. This leads to the following banded version of Algorithm 4.2.1.

Algorithm 4.3.5 (Band Cholesky) Given a symmetric positive definite $A \in \mathbb { R } ^ { n \times n }$ with bandwidth $p ,$ the following algorithm computes a lower triangular matrix $G$ with lower bandwidth $p$ such that $A = G G ^ { T }$ . For all $i \geq j , G ( i , j )$ overwrites $A ( i , j )$ .

$$
\begin{array}{l} \text { for } j = 1: n \\ \text { for } k = \max (1, j - p): j - 1 \\ \lambda = \min (k + p, n) \\ A (j: \lambda , j) = A (j: \lambda , j) - A (j, k) \cdot A (j: \lambda , k) \\ \lambda = \min (j + p, n) \\ A (j: \lambda , j) = A (j: \lambda , j) / \sqrt {A (j , j)} \\ \end{array}
$$

end

If $n \gg p .$ , then this algorithm requires about $n ( p ^ { 2 } + 3 p )$ flops and n square roots. Of course, in a serious implementation an appropriate data structure for A should be used. For example, if we just store the nonzero lower triangular part, then a $( p + 1 ) – \mathrm { b y } – n$ array would suffice.

If our band Cholesky procedure is coupled with appropriate band triangular solve routines, then approximately $n p ^ { 2 } + 7 n p + 2 n$ flops and n square roots are required to solve $A x = b$ . For small $p$ it follows that the square roots represent a significant portion of the computation and it is preferable to use the LDLT approach. Indeed, a careful flop count of the steps $A = L D L ^ { T } , L y = b , D z = y$ , and $L ^ { T } x \ = \ z$ reveals that $n p ^ { 2 } + 8 n p + n$ flops and no square roots are needed.

# 4.3.6 Tridiagonal System Solving

As a sample narrow band $\mathrm { L D L } ^ { T }$ solution procedure, we look at the case of symmetric positive definite tridiagonal systems. Setting

$$
L = \left[ \begin{array}{c c c c} 1 & 0 & \dots & 0 \\ \ell_ {1} & 1 & & \vdots \\ \vdots & \ddots & \ddots & 0 \\ 0 & \dots & \ell_ {n - 1} & 1 \end{array} \right]
$$

and $D = \operatorname { d i a g } ( d _ { 1 } , \dotsc , d _ { n } )$ , we deduce from the equation $A = L D L ^ { T }$ that

$$
\begin{array}{l} a _ {1 1} \quad = d _ {1}, \\ a _ {k, k - 1} = \ell_ {k - 1} d _ {k - 1}, \quad k = 2: n, \\ a _ {k k} \quad = d _ {k} + \ell_ {k - 1} ^ {2} d _ {k - 1} = d _ {k} + \ell_ {k - 1} a _ {k, k - 1}, \quad k = 2: n. \\ \end{array}
$$

Thus, the $d _ { i }$ and $\ell _ { i }$ can be resolved as follows:

$$
\begin{array}{l} d _ {1} = a _ {1 1} \\ \ell_ {k - 1} = a _ {k, k - 1} / d _ {k - 1} \\ d _ {k} = a _ {k k} - \ell_ {k - 1} a _ {k, k - 1} \\ \end{array}
$$

end

To obtain the solution to Ax = b we solve $L y = b , \ D z \ = \ y$ , and $L ^ { T } x = z$ . With overwriting we obtain

Algorithm 4.3.6 (Symmetric, Tridiagonal, Positive Definite System Solver) Given an n-by-n symmetric, tridiagonal, positive definite matrix A and $b \in \mathbb { R } ^ { n }$ , the following algorithm overwrites b with the solution to $A x = b $ . It is assumed that the diagonal of A is stored in α(1:n) and the superdiagonal in $\beta ( 1 : n - 1 )$ .

for k = 2:n

$$
t = \beta (k - 1), \beta (k - 1) = t / \alpha (k - 1), \alpha (k) = \alpha (k) - t \cdot \beta (k - 1)
$$

end

for k = 2:n

$$
b (k) = b (k) - \beta (k - 1) \cdot \beta (k - 1)
$$

end

$$
b (n) = b (n) / \alpha (n)
$$

for k = n − 1: − 1:1

$$
b (k) = b (k) / \alpha (k) - \beta (k) \cdot b (k + 1)
$$

end

This algorithm requires 8n flops.

# 4.3.7 Vectorization Issues

The tridiagonal example brings up a sore point: narrow band problems and vectorization do not mix. However, it is sometimes the case that large, independent sets of such problems must be solved at the same time. Let us examine how such a computation could be arranged in light of the issues raised in §1.5. For simplicity, assume that we must solve the n-by-n unit lower bidiagonal systems

$$
A ^ {(k)} x ^ {(k)} = b ^ {(k)}, \qquad k = 1: m,
$$

and that $m \gg n$ . Suppose we have arrays $E ( 1 { : } n - 1 , 1 { : } m )$ and $B ( 1 { : } n , 1 { : } m )$ with the property that $E ( 1 { : } n - 1 , k )$ houses the subdiagonal of $A ^ { ( k ) }$ and $B ( 1 { : } n , k )$ houses the kth right-hand side $b ^ { ( k ) }$ . We can overwrite $b ^ { ( k ) }$ with the solution $x ^ { ( k ) }$ as follows:

for k = 1:m

for i = 2:n

$$
B (i, k) = B (i, k) - E (i - 1, k) \cdot B (i - 1, k)
$$

end

end

This algorithm sequentially solves each bidiagonal system in turn. Note that the inner loop does not vectorize because of the dependence of $B ( i , k )$ on $B ( i - 1 , k )$ . However, if we interchange the order of the two loops, then the calculation does vectorize:

for $i = 2 { : } n$

$$
B (i,:) = B (i,:) - E (i - 1,:). * B (i - 1,:)
$$

end

A column-oriented version can be obtained simply by storing the matrix subdiagonals by row in E and the right-hand sides by row in B:

for $i = 2 { : } n$

$$
B (:, i) = B (:, i) - E (:, i - 1). * B (:, i - 1)
$$

end

Upon completion, the transpose of solution $x ^ { ( k ) }$ is housed on $B ( k , : )$ .

# 4.3.8 The Inverse of a Band Matrix

In general, the inverse of a nonsingular band matrix A is full. However, the off-diagonal blocks of $A ^ { - 1 }$ have low rank.

Theorem 4.3.3. Suppose

$$
A = \left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right]
$$

is nonsingular and has lower bandwidth p and upper bandwidth q. Assume that the diagonal blocks are square. If

$$
A ^ {- 1} = X = \left[ \begin{array}{l l} X _ {1 1} & X _ {1 2} \\ X _ {2 1} & X _ {2 2} \end{array} \right]
$$

is partitioned conformably, then

$$
\operatorname{rank} (X _ {2 1}) \leq p, \tag {4.3.1}
$$

$$
\operatorname{rank} (X _ {1 2}) \leq q. \tag {4.3.2}
$$

Proof. Assume that $A _ { 1 1 }$ and $A _ { 2 2 }$ are nonsingular. From the equation $A X = I$ we conclude that

$$
A _ {2 1} X _ {1 1} + A _ {2 2} X _ {2 1} = 0,
$$

$$
A _ {1 1} X _ {1 2} + A _ {1 2} X _ {2 2} = 0,
$$

and so

$$
\operatorname{rank} \left(X _ {2 1}\right) = \operatorname{rank} \left(A _ {2 2} ^ {- 1} A _ {2 1} X _ {1 1}\right) \leq \operatorname{rank} \left(A _ {2 1}\right)
$$

$$
\operatorname{rank} \left(X _ {1 2}\right) = \operatorname{rank} \left(A _ {1 1} ^ {- 1} A _ {1 2} X _ {2 2}\right) \leq \operatorname{rank} \left(A _ {1 2}\right).
$$

From the bandedness assumptions it follows that $A _ { 2 1 }$ has at most p nonzero rows and $A _ { 1 2 }$ has at most q nonzero rows. Thus, rank $( A _ { 2 1 } ) \le p$ and rank $( A _ { 1 2 } ) \le q$ which proves the theorem for the case when both $A _ { 1 1 }$ and $A _ { 2 2 }$ are nonsingular. A simple limit argument can be used to handle the situation when $A _ { 1 1 }$ and/or $A _ { 2 2 }$ are singular. See P4.3.11.

It can actually be shown that rank $( A _ { 2 1 } ) = \mathsf { r a n k } ( X _ { 2 1 } )$ and rank $( A _ { 1 2 } ) = \mathsf { r a n k } ( X _ { 1 2 } )$ . See Strang and Nguyen (2004). As we will see in §11.5.9 and §12.2, the low-rank, offdiagonal structure identified by the theorem has important algorithmic ramifications.

# 4.3.9 Band Matrices with Banded Inverse

If $A \in \mathbb { R } ^ { n \times n }$ is a product

$$
A = F _ {1} \dots F _ {N} \tag {4.3.3}
$$

and each $F _ { i } \in \mathbb { R } ^ { n \times n }$ is block diagonal with 1-by-1 and 2-by-2 diagonal blocks, then it follows that both A and

$$
A ^ {- 1} = F _ {N} ^ {- 1} \dots F _ {1} ^ {- 1}
$$

are banded, assuming that N is not too big. For example, if

$$
A = \left[ \begin{array}{c c c c c c c c c} \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & \times & \times & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & \times & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & \times & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & \times & \times & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times \end{array} \right] \left[ \begin{array}{c c c c c c c c c} \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & \times & \times & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & \times & \times & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & \times & \times & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & \times & \times & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times \end{array} \right]
$$

then

$$
A = \left[ \begin{array}{c c c c c c c c c} \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & \times & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & \times & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & \times & \times & \times & \times & 0 & 0 & 0 \\ 0 & 0 & \times & \times & \times & \times & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & \times & \times & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times & \times & \times & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times \end{array} \right], \quad A ^ {- 1} = \left[ \begin{array}{c c c c c c c c c} \times & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & \times & \times & \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & \times & \times & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & \times & \times & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & \times & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times \end{array} \right].
$$

Strang (2010a, 2010b) has pointed out a very important “reverse” fact. If A and $A ^ { - 1 }$ are banded, then there is a factorization of the form (4.3.3) with relatively small N . Indeed, he shows that N is very small for certain types of matrices that arise in signal processing. An important consequence of this is that both the forward transform Ax and the inverse transform $A ^ { - 1 } x$ can be computed very fast.

# Problems

P4.3.1 Develop a version of Algorithm 4.3.1 which assumes that the matrix A is stored in band format style. (See §1.2.5.)

P4.3.2 Show how the output of Algorithm 4.3.4 can be used to solve the upper Hessenberg system $H x = b $ .

P4.3.3 Show how Algorithm 4.3.4 could be used to solve a lower hessenberg system $H x = b$ .

P4.3.4 Give an algorithm for solving an unsymmetric tridiagonal system $A x = b$ that uses Gaussian elimination with partial pivoting. It should require only four n-vectors of floating point storage for the factorization.

P4.3.5 (a) For $C \in \mathbb { R } ^ { n \times n }$ define the profile indices $m ( C , i ) = \operatorname* { m i n } \{ j { : } c _ { i j } \neq 0 \}$ , where $i = 1 { : } n$ . Show that if $\mathbf { \Psi } _ { A } \stackrel { \mathbf { \tilde { = } } } { = } \mathbf { \Psi } _ { G } \mathbf { { \cal { G } } } ^ { T }$ is the Cholesky factorization of A, then $m ( A , i ) = m ( G , i )$ for $i = 1 { : } n$ (We say that G has the same profile as A.) (b) Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite with profile indices $m _ { i } = m ( A , i )$ where $i = 1 { : } n$ . Assume that A is stored in a one-dimensional array v as follows:

$$
v = (a _ {1 1}, a _ {2, m _ {2}}, \ldots , a _ {2 2}, a _ {3, m _ {3}}, \ldots , a _ {3 3}, \ldots , a _ {n, m _ {n}}, \ldots , a _ {n n}).
$$

Give an algorithm that overwrites v with the corresponding entries of the Cholesky factor G and then uses this factorization to solve $A x = b ,$ . How many flops are required? (c) For $C \in \mathbb { R } ^ { n \times n }$ define $p ( C , i )$ $= \operatorname* { m a x } \{ j { : } c _ { i j } \neq 0 \}$ . Suppose that $A \in \mathbb { R } ^ { n \times n }$ has an LU factorization $A = L U$ and that

$$
m (A, 1) \leq m (A, 2) \leq \dots \leq m (A, n),
$$

$$
p (A, 1) \leq p (A, 2) \leq \dots \leq p (A, n).
$$

Show that $m ( A , i ) = m ( L , i )$ and $p ( A , i ) = p ( U , i )$ for $i = 1 { : } n$ .

P4.3.6 Develop a gaxpy version of Algorithm 4.3.1.

P4.3.7 Develop a unit stride, vectorizable algorithm for solving the symmetric positive definite tridiagonal systems $A ^ { ( k ) } x ^ { ( k ) } = b ^ { ( \dot { k } ) }$ . Assume that the diagonals, superdiagonals, and right hand sides are stored by row in arrays $D , E ,$ and B and that $b ^ { ( k ) }$ is overwritten with $x ^ { ( k ) }$ .

P4.3.8 Give an example of a 3-by-3 symmetric positive definite matrix whose tridiagonal part is not positive definite.

P4.3.9 Suppose a symmetric positive definite matrix $A \in \mathbb { R } ^ { n \times n }$ has the “arrow structure”, e.g.,

$$
A = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & \times & 0 & 0 & 0 \\ \times & 0 & \times & 0 & 0 \\ \times & 0 & 0 & \times & 0 \\ \times & 0 & 0 & 0 & \times \end{array} \right].
$$

(a) Show how the linear system $A x = b$ can be solved with $O ( n )$ flops using the Sherman-Morrison-Woodbury formula. (b) Determine a permutation matrix P so that the Cholesky factorization

$$
P A P ^ {T} = G G ^ {T}
$$

can be computed with O(n) flops.

P4.3.10 Suppose $A \in \mathbb { R } ^ { n \times n }$ is tridiagonal, positive definite, but not symmetric. Give an efficient algorithm for computing the largest entry of $\vert S T ^ { - 1 } S \vert$ where $S = ( A - A ^ { T } ) / 2$ and $T = ( A + A ^ { T } ) / 2$ .

P4.3.11 Show that if $A \in \mathbb { R } ^ { n \times n }$ and $\epsilon > 0$ , then there is a $B \in \mathbb { R } ^ { n \times n }$ such that $\| A - B \| \leq \epsilon$ and B has the property that all its principal submatrices are nonsingular. Use this result to formally complete the proof of Theorem 4.3.3.

P4.3.12 Give an upper bound on the bandwidth of the matrix A in (4.3.3).

P4.3.13 Show that $A ^ { T }$ and $A ^ { - 1 }$ have the same upper and lower bandwidths in (4.3.3).

P4.3.14 For the $A = F _ { 1 } F _ { 2 }$ example in 4.3.9, show that A(2:3, :), A(4:5, :), A(6:7, :), . . . each consist of two singular 2-by-2 blocks.

# Notes and References for §4.3

Representative papers on the topic of banded systems include:

R.S. Martin and J.H. Wilkinson (1965). “Symmetric Decomposition of Positive Definite Band Matrices,” Numer. Math. 7, 355–61.   
R. S. Martin and J.H. Wilkinson (1967). “Solution of Symmetric and Unsymmetric Band Equations and the Calculation of Eigenvalues of Band Matrices,” Numer. Math. 9, 279–301.   
E.L. Allgower (1973). “Exact Inverses of Certain Band Matrices,” Numer. Math. 21, 279–284.   
Z. Bohte (1975). “Bounds for Rounding Errors in the Gaussian Elimination for Band Systems,” J. Inst. Math. Applic. 16, 133–142.   
L. Kaufman (2007). “The Retraction Algorithm for Factoring Banded Symmetric Matrices,” Numer. Lin. Alg. Applic. 14, 237–254.   
C. Vomel and J. Slemons (2009). “Twisted Factorization of a Banded Matrix,” BIT 49, 433–447.

Tridiagonal systems are particularly important, see:

C. Fischer and R.A. Usmani (1969). “Properties of Some Tridiagonal Matrices and Their Application to Boundary Value Problems,” SIAM J. Numer. Anal. 6, 127–142.   
D.J. Rose (1969). “An Algorithm for Solving a Special Class of Tridiagonal Systems of Linear Equations,” Commun. ACM 12, 234–236.   
M.A. Malcolm and J. Palmer (1974). “A Fast Method for Solving a Class of Tridiagonal Systems of Linear Equations,” Commun. ACM 17, 14–17.   
N.J. Higham (1986). “Efficient Algorithms for Computing the Condition Number of a Tridiagonal Matrix,” SIAM J. Sci. Stat. Comput. 7, 150–165.   
N.J. Higham (1990). “Bounding the Error in Gaussian Elimination for Tridiagonal Systems,” SIAM J. Matrix Anal. Applic. 11, 521–530.   
I.S. Dhillon (1998). “Reliable Computation of the Condition Number of a Tridiagonal Matrix in O(n) Time,” SIAM J. Matrix Anal. Applic. 19, 776–796.   
I. Bar-On and M. Leoncini (2000). “Reliable Solution of Tridiagonal Systems of Linear Equations,” SIAM J. Numer. Anal. 38, 1134–1153.   
M.I. Bueno and F.M. Dopico (2004). “Stability and Sensitivity of Tridiagonal LU Factorization without Pivoting,” BIT 44, 651–673.   
J.R. Bunch and R.F. Marcia (2006). “A Simplified Pivoting Strategy for Symmetric Tridiagonal Matrices,” Numer. Lin. Alg. 13, 865–867.

For a discussion of parallel methods for banded problems, see:

H.S. Stone (1975). “Parallel Tridiagonal Equation Solvers,” ACM Trans. Math. Softw. 1, 289–307.   
I. Bar-On, B. Codenotti and M. Leoncini (1997). “A Fast Parallel Cholesky Decomposition Algorithm for Tridiagonal Symmetric Matrices,” SIAM J. Matrix Anal. Applic. 18, 403–418.   
G.H. Golub, A.H. Sameh, and V. Sarin (2001). “A Parallel Balance Scheme for Banded Linear Systems,” Num. Lin. Alg. 8, 297–316.   
S. Rao and Sarita (2008). “Parallel Solution of Large Symmetric Tridiagonal Linear Systems,” Parallel Comput. 34, 177–197.

Papers that are concerned with the structure of the inverse of a band matrix include:

E. Asplund (1959). “Inverses of Matrices {aij } Which Satisfy aij = 0 for j > i + p,” Math. Scand. 7, 57–60.   
C.A. Micchelli (1992). “Banded Matrices with Banded Inverses,” J. Comput. Appl. Math. 41, 281-300.   
G. Strang and T. Nguyen (2004). “The Interplay of Ranks of Submatrices,” SIAM Review 46, 637–648.   
G. Strang (2010a). “Fast Transforms: Banded Matrices with Banded Inverses,” Proc. National Acad. Sciences 107, 12413-12416.   
G. Strang (2010b). “Banded Matrices with Banded Inverses and A = LP U,” Proceedings International Congress of Chinese Mathematicians, Beijing.

A pivotal result in this arena is the nullity theorem, a more general version of Theorem 4.3.3, see:   
R. Vandebril, M. Van Barel, and N. Mastronardi (2008). Matrix Computations and Semiseparable Matrices, Volume I Linear Systems, Johns Hopkins University Press, Baltimore, MD., 37–40.

# 4.4 Symmetric Indefinite Systems

Recall that a matrix whose quadratic form $x ^ { T }$ Ax takes on both positive and negative values is indefinite. In this section we are concerned with symmetric indefinite linear systems. The $\mathrm { L D L ^ { T } }$ factorization is not always advisable as the following 2-by-2 example illustrates:

$$
\left[ \begin{array}{c c} \epsilon & 1 \\ 1 & 0 \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ 1 / \epsilon & 1 \end{array} \right] \left[ \begin{array}{c c} \epsilon & 0 \\ 0 & - 1 / \epsilon \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ 1 / \epsilon & 1 \end{array} \right] ^ {T}.
$$

Of course, any of the pivot strategies in §3.4 could be invoked. However, they destroy symmetry and, with it, the chance for a “Cholesky speed” symmetric indefinite system solver. Symmetric pivoting, i.e., data reshufflings of the form $A  P A P ^ { T }$ , must be used as we discussed in §4.2.8. Unfortunately, symmetric pivoting does not always stabilize the $\mathrm { L D L ^ { T } }$ computation. If $\epsilon _ { 1 }$ and $\epsilon _ { 2 }$ are small, then regardless of $P _ { \cdot }$ the matrix

$$
\tilde {A} = P \left[ \begin{array}{c c} \epsilon_ {1} & 1 \\ 1 & \epsilon_ {2} \end{array} \right] P ^ {T}
$$

has small diagonal entries and large numbers surface in the factorization. With symmetric pivoting, the pivots are always selected from the diagonal and trouble results if these numbers are small relative to what must be zeroed off the diagonal. Thus, $\mathrm { L D L ^ { T } }$ with symmetric pivoting cannot be recommended as a reliable approach to symmetric indefinite system solving. It seems that the challenge is to involve the off-diagonal entries in the pivoting process while at the same time maintaining symmetry.

In this section we discuss two ways to do this. The first method is due to Aasen (1971) and it computes the factorization

$$
P A P ^ {T} = L T L ^ {T}, \tag {4.4.1}
$$

where $L = ( \ell _ { i j } )$ is unit lower triangular and $T$ is tridiagonal. P is a permutation chosen such that $| \ell _ { i j } | \le 1$ . In contrast, the diagonal pivoting method due to Bunch and Parlett (1971) computes a permutation P such that

$$
P A P ^ {T} = L D L ^ {T}, \tag {4.4.2}
$$

where D is a direct sum of 1-by-1 and 2-by-2 pivot blocks. Again, $P$ is chosen so that the entries in the unit lower triangular L satisfy $| \ell _ { i j } | \le 1$ . Both factorizations involve $n ^ { 3 } / 3$ flops and once computed, can be used to solve $A x = b$ with $O ( n ^ { 2 } )$ work:

$$
P A P ^ {T} = L T L ^ {T}, L z = P b, T w = z, L ^ {T} y = w, x = P ^ {T} y \Rightarrow A x = b,
$$

$$
P A P ^ {T} = L D L ^ {T}, L z = P b, D w = z, L ^ {T} y = w, x = P ^ {T} y \Rightarrow A x = b.
$$

A few comments need to be made about the $T w = z$ and $D w = z$ systems that arise when these methods are invoked.

In Aasen’s method, the symmetric indefinite tridiagonal system $T w = z$ is solved in $O ( n )$ time using band Gaussian elimination with pivoting. Note that there is no serious price to pay for the disregard of symmetry at this level since the overall process is $O ( n ^ { 3 } )$ .

In the diagonal pivoting approach, the $D w = z$ system amounts to a set of 1-by-1 and $2 \mathrm { - b y { - } 2 }$ symmetric indefinite systems. The 2-by-2 problems can be handled via Gaussian elimination with pivoting. Again, there is no harm in disregarding symmetry during this $O ( n )$ phase of the calculation. Thus, the central issue in this section is the efficient computation of the factorizations (4.4.1) and (4.4.2).

# 4.4.1 The Parlett-Reid Algorithm

Parlett and Reid (1970) show how to compute (4.4.1) using Gauss transforms. Their algorithm is sufficiently illustrated by displaying the k = 2 step for the case $n = 5$ . At the beginning of this step the matrix A has been transformed to

$$
A ^ {(1)} = M _ {1} P _ {1} A P _ {1} ^ {T} M _ {1} ^ {T} = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & 0 & 0 & 0 \\ \beta_ {1} & \alpha_ {2} & v _ {3} & v _ {4} & v _ {5} \\ 0 & v _ {3} & \times & \times & \times \\ 0 & v _ {4} & \times & \times & \times \\ 0 & v _ {5} & \times & \times & \times \end{array} \right],
$$

where $P _ { 1 }$ is a permutation chosen so that the entries in the Gauss transformation $M _ { 1 }$ are bounded by unity in modulus. Scanning the vector $[ v _ { 3 } \ v _ { 4 } \ v _ { 5 } ] ^ { T }$ for its largest entry, we now determine a 3-by-3 permutation ${ \tilde { P } } _ { 2 }$ such that

$$
\tilde {P} _ {2} \left[ \begin{array}{l} v _ {3} \\ v _ {4} \\ v _ {5} \end{array} \right] = \left[ \begin{array}{l} \tilde {v} _ {3} \\ \tilde {v} _ {4} \\ \tilde {v} _ {5} \end{array} \right] \qquad \Rightarrow \qquad | \tilde {v} _ {3} | = \max \{| \tilde {v} _ {3} |, | \tilde {v} _ {4} |, | \tilde {v} _ {5} | \}.
$$

If this maximal element is zero, we set $M _ { 2 } = P _ { 2 } = I$ and proceed to the next step. Otherwise, we set $P _ { 2 } = \mathrm { d i a g } ( I _ { 2 } , \tilde { P } _ { 2 } )$ and $M _ { 2 } = I \ - \ \alpha ^ { ( 2 ) } e _ { 3 } ^ { T }$ with

$$
\alpha^ {(2)} = \left[ \begin{array}{l l l l l} 0 & 0 & 0 & \tilde {v} _ {4} / \tilde {v} _ {3} & \tilde {v} _ {5} / \tilde {v} _ {3} \end{array} \right] ^ {T}.
$$

Observe that

$$
A ^ {(2)} = M _ {2} P _ {2} A ^ {(1)} P _ {2} ^ {T} M _ {2} ^ {T} = \left[ \begin{array}{l l l l l} \alpha_ {1} & \beta_ {1} & 0 & 0 & 0 \\ \beta_ {1} & \alpha_ {2} & \tilde {v} _ {3} & 0 & 0 \\ 0 & \tilde {v} _ {3} & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & \times & \times & \times \end{array} \right].
$$

In general, the process continues for $n - 2$ steps leaving us with a tridiagonal matrix

$$
T = A ^ {(n - 2)} = (M _ {n - 2} P _ {n - 2} \dots M _ {1} P _ {1}) A (M _ {n - 2} P _ {n - 2} \dots M _ {1} P _ {1}) ^ {T}.
$$

It can be shown that (4.4.1) holds with $P = P _ { n - 2 } \cdot \cdot \cdot P _ { 1 }$ and

$$
L = \left(M _ {n - 2} P _ {n - 2} \dots M _ {1} P _ {1} P ^ {T}\right) ^ {- 1}.
$$

Analysis of L reveals that its first column is $e _ { 1 }$ and that its subdiagonal entries in column k with $k > 1$ are “made up” of the multipliers in $M _ { k - 1 }$ .

The efficient implementation of the Parlett-Reid method requires care when computing the update

$$
A ^ {(k)} = M _ {k} (P _ {k} A ^ {(k - 1)} P _ {k} ^ {T}) M _ {k} ^ {T}. \tag {4.4.3}
$$

To see what is involved with a minimum of notation, suppose B = BT ∈ IR(n−k)×(n−k) $B = B ^ { T } \in \mathbb { R } ^ { ( n - k ) \times ( n - k ) }$ has and that we wish to form

$$
B _ {+} = (I - w e _ {1} ^ {T}) B (I - w e _ {1} ^ {T}) ^ {T},
$$

where $w \in \mathbb { R } ^ { n - k }$ and $e _ { 1 }$ is the first column of $I _ { n - k }$ . Such a calculation is at the heart of (4.4.3). If we set

$$
u = B e _ {1} - \frac {b _ {1 1}}{2} w,
$$

then $\begin{array} { r } { B _ { + } = B - w u ^ { T } - u w ^ { T } } \end{array}$ and its lower triangular portion can be formed in $2 ( n - k ) ^ { 2 }$ flops. Summing this quantity as k ranges from 1 to $n - 2$ indicates that the Parlett-Reid procedure requires $2 n ^ { 3 } / 3$ flops—twice the volume of work associated with Cholesky.

# 4.4.2 The Method of Aasen

An $n ^ { 3 } / 3$ approach to computing (4.4.1) due to Aasen (1971) can be derived by reconsidering some of the computations in the Parlett-Reid approach. We examine the no-pivoting case first where the goal is to compute a unit lower triangular matrix L with $L ( : , 1 ) = e _ { 1 }$ and a tridiagonal matrix

$$
T = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \beta_ {n - 1} \\ 0 & \dots & & \beta_ {n - 1} & \alpha_ {n} \end{array} \right].
$$

such that $A = L T L ^ { T }$ . The Aasen method is structured as follows:

for $j = 1:n$ $\{\alpha(1:j - 1), \beta(1:j - 1) \text{ and } L(:,1:j) \text{ are known}\}$ Compute $\alpha_j$ .  
    if $j \leq n - 1$ Compute $\beta_j$ .  
    end  
    if $j \leq n - 2$ Compute $L(j + 2:n, j + 1)$ .  
    end  
end

To develop recipes for $\alpha _ { j } , \beta _ { j }$ , and $L ( j + 2 ; n , j + 1 )$ , we compare the jth columns in the equation $A = L H$ where $\boldsymbol { H } ^ { \setminus } = \boldsymbol { T } \boldsymbol { L } ^ { T }$ . Noting that H is an upper Hessenberg matrix we obtain

$$
A (:, j) = L H (:, j) = \sum_ {k = 1} ^ {j + 1} L (:, k) \cdot h (k), \tag {4.4.5}
$$

where $h ( 1 { : } j + 1 ) = H ( 1 { : } j + 1 , j )$ and we assume that $j \le n - 1$ . It follows that

$$
h _ {j + 1} \cdot L (j + 1: n, j + 1) = v (j + 1: n), \tag {4.4.6}
$$

where

$$
v (j + 1: n) = A (j + 1: n, j) - L (j + 1: n, 1: j) \cdot h (1: j). \tag {4.4.7}
$$

Since L is unit lower triangular and $L ( : , 1 : j )$ is known, this gives us a working recipe for $L ( j + 2 { : } n , j + 1 )$ provided we know $h ( 1 { : } j )$ . Indeed, from (4.4.6) and (4.4.7) it is easy to show that

$$
L (j + 2: n, j + 1) = v (j + 2: n) / v (j + 1). \tag {4.4.8}
$$

To compute $h ( 1 { : } j )$ we turn to the equation $H = T L ^ { T }$ and examine its jth column. The case $j = 5$ amply displays what is going on:

$$
\left[ \begin{array}{l} h _ {1} \\ h _ {2} \\ h _ {3} \\ h _ {4} \\ h _ {5} \\ h _ {6} \end{array} \right] = \left[ \begin{array}{l l l l l} \alpha_ {1} & \beta_ {1} & 0 & 0 & 0 \\ \beta_ {1} & \alpha_ {2} & \beta_ {2} & 0 & 0 \\ 0 & \beta_ {2} & \alpha_ {3} & \beta_ {3} & 0 \\ 0 & 0 & \beta_ {3} & \alpha_ {4} & \beta_ {4} \\ 0 & 0 & 0 & \beta_ {4} & \alpha_ {5} \\ 0 & 0 & 0 & 0 & \beta_ {5} \end{array} \right] \left[ \begin{array}{c} 0 \\ \ell_ {5 2} \\ \ell_ {5 3} \\ \ell_ {5 4} \\ 1 \end{array} \right] = \left[ \begin{array}{c} \beta_ {1} \ell_ {5 2} \\ \alpha_ {2} \ell_ {5 2} + \beta_ {2} \ell_ {5 3} \\ \beta_ {2} \ell_ {5 2} + \alpha_ {3} \ell_ {5 3} + \beta_ {3} \ell_ {5 4} \\ \beta_ {3} \ell_ {5 3} + \alpha_ {4} \ell_ {5 4} + \beta_ {4} \\ \beta_ {4} \ell_ {5 4} + \alpha_ {5} \\ \beta_ {5} \end{array} \right] \tag {4.4.9}
$$

At the start of step $j ,$ we know $\alpha ( 1 { : } j - 1 ) , \beta ( 1 { : } j - 1 )$ and $L ( : , 1 : j )$ . Thus, we can determine $h ( 1 { : } j - 1 )$ as follows

$$
h _ {1} = \beta_ {1} \ell_ {j 2}
$$

for $k = 1 { : } j - 1$

$$
h _ {k} = \beta_ {k - 1} \ell_ {j, k - 1} + \alpha_ {k} \ell_ {j k} + \beta_ {k} \ell_ {j, k + 1} \tag {4.4.10}
$$

end

Equation (4.4.5) gives us a formula for $h _ { j }$ :

$$
h _ {j} = A (j, j) - \sum_ {k = 1} ^ {j - 1} L (j, k) h _ {k}. \tag {4.4.11}
$$

From (4.4.9) we infer that

$$
\alpha_ {j} = h _ {j} - \beta_ {j - 1} \ell_ {j, j - 1}, \tag {4.4.12}
$$

$$
\beta_ {j} = h _ {j + 1}. \tag {4.4.13}
$$

Combining these equations with (4.4.4), (4.4.7), (4.4.8), (4.4.10), and (4.4.11) we obtain the Aasen method without pivoting:

$$
L = I _ {n}
$$

for j = 1:n

$$
\text { if } j = 1
$$

$$
\alpha_ {1} = a _ {1 1}
$$

$$
v (2: n) = A (2: n, 1)
$$

else

$$
h _ {1} = \beta_ {1} \cdot \ell_ {j 2}
$$

for k = 2:j − 1

$$
h _ {k} = \beta_ {k - 1} \ell_ {j, k - 1} + \alpha_ {k} \ell_ {j k} + \beta_ {k} \ell_ {j, k + 1}
$$

end

$$
h _ {j} = a _ {j j} - L (j, 1: j - 1) \cdot h (1: j - 1)
$$

$$
\alpha_ {j} = h _ {j} - \beta_ {j - 1} \ell_ {j, j - 1} \tag {4.4.14}
$$

$$
v (j + 1: n) = A (j + 1: n, j) - L (j + 1: n, 1: j) \cdot h (1: j)
$$

end

$$
\text { if   } j <   = n - 1
$$

$$
\beta_ {j} = v (j + 1)
$$

end

$$
\text { if } j <   = n - 2
$$

$$
L (j + 2: n, j + 1) = v (j + 2: n) / v (j + 1)
$$

end

end

The dominant operation each pass through the j-loop is an (n−j)-by-j gaxpy operation. Accounting for the associated flops we see that the overall Aasen ccomputation involves $n ^ { 3 } / 3$ flops, the same as for the Cholesky factorization.

As it now stands, the columns of L are scalings of the v-vectors in (4.4.14). If any of these scalings are large, i.e., if any $v ( j + 1 )$ is small, then we are in trouble. To circumvent this problem, it is only necessary to permute the largest component of $v ( j + 1 { : } n )$ to the top position. Of course, this permutation must be suitably applied to the unreduced portion of A and the previously computed portion of L. With pivoting, Aasen’s method is stable in the same sense that Gaussian elimination with partial pivoting is stable.

In a practical implementation of the Aasen algorithm, the lower triangular portion of A would be overwritten with L and T , e.g.,

$$
A \leftarrow \left[ \begin{array}{c c c c c} \alpha_ {1} & & & & \\ \beta_ {1} & \alpha_ {2} & & & \\ \ell_ {3 2} & \beta_ {2} & \alpha_ {3} & & \\ \ell_ {4 2} & \ell_ {4 3} & \beta_ {3} & \alpha_ {4} & \\ \ell_ {5 2} & \ell_ {5 3} & \ell_ {5 4} & \beta_ {4} & \alpha_ {5} \end{array} \right].
$$

Notice that the columns of L are shifted left in this arrangement.

# 4.4.3 Diagonal Pivoting Methods

We next describe the computation of the block $L D L ^ { T }$ factorization (4.4.2). We follow the discussion in Bunch and Parlett (1971). Suppose

$$
P _ {1} A P _ {1} ^ {T} = \left[ \begin{array}{c c} E & C ^ {T} \\ C & B \end{array} \right] _ {n - s} ^ {s}
$$

where $P _ { 1 }$ is a permutation matrix and s = 1 or 2. If A is nonzero, then it is always possible to choose these quantities so that E is nonsingular, thereby enabling us to write

$$
P _ {1} A P _ {1} ^ {T} = \left[ \begin{array}{c c} I _ {s} & 0 \\ C E ^ {- 1} & I _ {n - s} \end{array} \right] \left[ \begin{array}{c c} E & 0 \\ 0 & B - C E ^ {- 1} C ^ {T} \end{array} \right] \left[ \begin{array}{c c} I _ {s} & E ^ {- 1} C ^ {T} \\ 0 & I _ {n - s} \end{array} \right].
$$

For the sake of stability, the s-by-s “pivot” E should be chosen so that the entries in

$$
\tilde {A} = \left(\tilde {a} _ {i j}\right) \equiv B - C E ^ {- 1} C ^ {T} \tag {4.4.15}
$$

are suitably bounded. To this end, let $\alpha \in ( 0 , 1 )$ be given and define the size measures

$$
\begin{array}{l} \mu_ {0} = \max _ {i, j} | a _ {i j} |, \\ \mu_ {1} = \max _ {i} | a _ {i i} |. \\ \end{array}
$$

The Bunch-Parlett pivot strategy is as follows:

$$
\begin{array}{l} \mathbf {i f} \mu_ {1} \geq \alpha \mu_ {0} \\ s = 1 \\ \text { Choose } P _ {1} \text { so } | e _ {1 1} | = \mu_ {1}. \\ s = 2 \\ \text { Choose } P _ {1} \text { so } | e _ {2 1} | = \mu_ {0}. \\ \end{array}
$$

It is easy to verify from (4.4.15) that if s = 1, then

$$
\left| \tilde {a} _ {i j} \right| \leq \left(1 + \alpha^ {- 1}\right) \mu_ {0}, \tag {4.4.16}
$$

while s = 2 implies

$$
\left| \tilde {a} _ {i j} \right| \leq \frac {3 - \alpha}{1 - \alpha} \mu_ {0}. \tag {4.4.17}
$$

By equating $( 1 + \alpha ^ { - 1 } ) ^ { 2 }$ , the growth factor that is associated with two $s = 1$ steps, and $( 3 - \alpha ) / ( 1 - \alpha )$ , the corresponding s = 2 factor, Bunch and Parlett conclude that $\alpha = ( 1 + { \sqrt { 1 7 } } ) / 8$ is optimum from the standpoint of minimizing the bound on element growth.

The reductions outlined above can be repeated on the order- $( n - s )$ symmetric matrix ${ \tilde { A } } .$ . A simple induction argument establishes that the factorization (4.4.2) exists and that $n ^ { 3 } / 3$ flops are required if the work associated with pivot determination is ignored.

# 4.4.4 Stability and Efficiency

Diagonal pivoting with the above strategy is shown by Bunch (1971) to be as stable as Gaussian elimination with complete pivoting. Unfortunately, the overall process requires between $n ^ { 3 } / 1 2$ and $n ^ { 3 } / 6$ comparisons, since $\mu _ { 0 }$ involves a two-dimensional search at each stage of the reduction. The actual number of comparisons depends on the total number of 2-by-2 pivots but in general the Bunch-Parlett method for computing (4.4.2) is considerably slower than the technique of Aasen. See Barwell and George (1976).

This is not the case with the diagonal pivoting method of Bunch and Kaufman (1977). In their scheme, it is only necessary to scan two columns at each stage of the reduction. The strategy is fully illustrated by considering the very first step in the reduction:

$\alpha = (1 + \sqrt{17}) / 8$ $\lambda = |a_{r1}| = \max \{|a_{21}|, \ldots, |a_{n1}|\}$ if $\lambda > 0$ if $|a_{11}| \geq \alpha \lambda$ Set $s = 1$ and $P_{1} = I$ .   
else $\sigma = |a_{pr}| = \max \{|a_{1r}, \ldots, |a_{r-1,r}|, |a_{r+1,r}|, \ldots, |a_{nr}|\}$ if $\sigma |a_{11}| \geq \alpha \lambda^2$ Set $s = 1$ and $P_{1} = I$ elseif $|a_{rr}| \geq \alpha \sigma$ Set $s = 1$ and choose $P_{1}$ so $(P_{1}^{T}AP_{1})_{11} = a_{rr}$ .   
else   
Set $s = 2$ and choose $P_{1}$ so $(P_{1}^{T}AP_{1})_{21} = a_{rp}$ .   
end   
end   
end

Overall, the Bunch-Kaufman algorithm requires $n ^ { 3 } / 3$ flops, $O ( n ^ { 2 } )$ comparisons, and, like all the methods of this section, $n ^ { 2 } / 2$ storage.

# 4.4.5 A Note on Equilibrium Systems

A very important class of symmetric indefinite matrices have the form

$$
A = \left[ \begin{array}{c c} C & B \\ B ^ {T} & 0 \\ n & p \end{array} \right] _ {p} ^ {n} \tag {4.4.18}
$$

where C is symmetric positive definite and B has full column rank. These conditions ensure that A is nonsingular.

Of course, the methods of this section apply to A. However, they do not exploit its structure because the pivot strategies “wipe out” the zero (2,2) block. On the other hand, here is a tempting approach that does exploit A’s block structure:

Step 1. Compute the Cholesky factorization $C = G G ^ { T }$ .

Step 2. Solve GK = B for $K \in \mathbb { R } ^ { n \times p }$ .

Step 3. Compute the Cholesky factorization $H H ^ { T } = K ^ { T } K = B ^ { T } C ^ { - 1 } B .$

From this it follows that

$$
A = \left[ \begin{array}{c c} G & 0 \\ K ^ {T} & H \end{array} \right] \left[ \begin{array}{c c} G ^ {T} & K \\ 0 & - H ^ {T} \end{array} \right].
$$

In principle, this triangular factorization can be used to solve the equilibrium system

$$
\left[ \begin{array}{c c} C & B \\ B ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} x \\ y \end{array} \right] = \left[ \begin{array}{l} f \\ g \end{array} \right]. \tag {4.4.19}
$$

However, it is clear by considering steps (b) and (c) above that the accuracy of the computed solution depends upon $\kappa ( C )$ and this quantity may be much greater than $\kappa ( A )$ . The situation has been carefully analyzed and various structure-exploiting algorithms have been proposed. A brief review of the literature is given at the end of the section.

It is interesting to consider a special case of (4.4.19) that clarifies what it means for an algorithm to be stable and illustrates how perturbation analysis can structure the search for better methods. In several important applications, $g = 0 , C$ is diagonal, and the solution subvector y is of primary importance. A manipulation shows that this vector is specified by

$$
y = (B ^ {T} C ^ {- 1} B) ^ {- 1} B ^ {T} C ^ {- 1} f. \tag {4.4.20}
$$

Looking at this we are again led to believe that $\kappa ( C )$ should have a bearing on the accuracy of the computed y. However, it can be shown that

$$
\left\| \left(B ^ {T} C ^ {- 1} B\right) ^ {- 1} B ^ {T} C ^ {- 1} \right\| \leq \psi_ {B} \tag {4.4.21}
$$

where the upper bound $\psi _ { B }$ is independent of C, a result that (correctly) suggests that y is not sensitive to perturbations in C. A stable method for computing this vector should respect this, meaning that the accuracy of the computed y should be independent of C. Vavasis (1994) has developed a method with this property. It involves the careful assembly of a matrix $V \in \mathbb { R } ^ { n \times ( n - p ) }$ whose columns are a basis for the nullspace of $B ^ { T } C ^ { - 1 }$ . The n-by-n linear system

$$
[ B \mid V ] \left[ \begin{array}{l} y \\ q \end{array} \right] = f
$$

is then solved implying $f = B y + V q$ . Thus, $B ^ { T } C ^ { - 1 } f = B ^ { T } C ^ { - 1 } B y$ and (4.4.20) holds.

# Problems

P4.4.1 Show that if all the 1-by-1 and 2-by-2 principal submatrices of an n-by-n symmetric matrix A are singular, then A is zero.

P4.4.2 Show that no 2-by-2 pivots can arise in the Bunch-Kaufman algorithm if A is positive definite.

P4.4.3 Arrange (4.4.14) so that only the lower triangular portion of A is referenced and so that $\alpha ( j )$ overwrites $A ( j , j )$ for $j = 1 { : } n , \beta ( j )$ overwrites $A ( j + 1 , j )$ for $j = 1 { : } n - 1$ , and $L ( i , j )$ overwrites $A ( i , j - 1 )$ for $j = 2 { : } n - 1$ and $i = j + 1 { : } n$ .

P4.4.4 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and strictly diagonally dominant. Give an algorithm that computes the factorization

$$
\Pi A \Pi^ {T} = \left[ \begin{array}{c c} R & 0 \\ S & - M \end{array} \right] \left[ \begin{array}{c c} R ^ {T} & S ^ {T} \\ 0 & M ^ {T} \end{array} \right]
$$

where Π is a permuation and the diagonal blocks R and M are lower triangular.

P4.4.5 A symmetric matrix A is quasidefinite if it has the form

$$
A = \left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & - A _ {2 2} \end{array} \right] _ {p} ^ {n}
$$

with $A _ { 1 1 }$ and $A _ { 2 2 }$ positive definite. (a) Show that such a matrix has an $\mathrm { L D L } ^ { T }$ factorization with the property that

$$
D = \left[ \begin{array}{c c} D _ {1} & 0 \\ 0 & - D _ {2} \end{array} \right]
$$

where $D _ { 1 } \in \mathbb { R } ^ { n \times n }$ and $D _ { 2 } \in \mathbb { R } ^ { p \times p }$ have positive diagonal entries. (b) Show that if A is quasidefinite then all its principal submatrices are nonsingular. This means that ${ \bf \vec { P } } { \bf \vec { { A } } } { \cal { P } } ^ { T }$ has an $\mathrm { L D L } ^ { T }$ factorization for any permutation matrix P .

P4.4.6 Prove (4.4.16) and (4.4.17).

P4.4.7 Show that $- ( B ^ { T } C ^ { - 1 } B ) ^ { - 1 }$ is the (2,2) block of $A ^ { - 1 }$ where A is given by equation (4.4.18).

P4.4.8 The point of this problem is to consider a special case of (4.4.21). Define the matrix

$$
M (\alpha) = (B ^ {T} C ^ {- 1} B) ^ {- 1} B ^ {T} C ^ {- 1}
$$

where $C = ( I _ { n } + \alpha e _ { k } e _ { k } ^ { T } )$ , α > −1, and $e _ { k } = I _ { n } ( : , k )$ . (Note that C is just the identity with α added to the (k, k) entry.) Assume that $B \in \mathbb { R } ^ { n \times p }$ has rank p and show that

$$
M (\alpha) = (B ^ {T} B) ^ {- 1} B ^ {T} \left(I _ {n} - \frac {\alpha}{1 + \alpha w ^ {T} w} e _ {k} w ^ {T}\right)
$$

where

$$
w = (I _ {n} - B (B ^ {T} B) ^ {- 1} B ^ {T}) e _ {k}.
$$

Show that if $\parallel w \parallel _ { 2 } = 0 \mathrm { o r } \parallel w \parallel _ { 2 } = 1$ , then $\parallel M ( \alpha ) \parallel _ { 2 } = 1 / \sigma _ { \mathrm { m i n } } ( B )$ . Show that if $0 < \parallel w \parallel _ { 2 } < 1$ , then

$$
\| M (\alpha) \| _ {2} \leq \max \left\{\frac {1}{1 - \| w \| _ {2}}, 1 + \frac {1}{\| w \| _ {2}} \right\} \Bigg / \sigma_ {\min} (B).
$$

Thus,  $M ( \alpha ) \parallel _ { 2 }$ has an α-independent upper bound.

# Notes and References for §4.4

The basic references for computing (4.4.1) are as follows:

J.O. Aasen (1971). “On the Reduction of a Symmetric Matrix to Tridiagonal Form,” BIT 11, 233–242.

B.N. Parlett and J.K. Reid (1970). “On the Solution of a System of Linear Equations Whose Matrix Is Symmetric but not Definite,”BIT 10, 386–397.

The diagonal pivoting literature includes:

J.R. Bunch and B.N. Parlett (1971). “Direct Methods for Solving Symmetric Indefinite Systems of Linear Equations,” SIAM J. Numer. Anal. 8, 639–655.

J.R. Bunch (1971). “Analysis of the Diagonal Pivoting Method,” SIAM J. Numer. Anal. 8, 656–680.

J.R. Bunch (1974). “Partial Pivoting Strategies for Symmetric Matrices,” SIAM J. Numer. Anal. 11, 521–528.

J.R. Bunch, L. Kaufman, and B.N. Parlett (1976). “Decomposition of a Symmetric Matrix,” Numer. Math. 27, 95–109.   
J.R. Bunch and L. Kaufman (1977). “Some Stable Methods for Calculating Inertia and Solving Symmetric Linear Systems,” Math. Comput. 31, 162–79.   
M.T. Jones and M.L. Patrick (1993). “Bunch-Kaufman Factorization for Real Symmetric Indefinite Banded Matrices,” SIAM J. Matrix Anal. Applic. 14, 553–559.   
Because “future” columns must be scanned in the pivoting process, it is awkward (but possible) to obtain a gaxpy-rich diagonal pivoting algorithm. On the other hand, Aasen’s method is naturally rich in gaxpys. Block versions of both procedures are possible. Various performance issues are discussed in:   
V. Barwell and J.A. George (1976). “A Comparison of Algorithms for Solving Symmetric Indefinite Systems of Linear Equations,” ACM Trans. Math. Softw. 2, 242–251.   
M.T. Jones and M.L. Patrick (1994). “Factoring Symmetric Indefinite Matrices on High-Performance Architectures,” SIAM J. Matrix Anal. Applic. 15, 273–283.   
Another idea for a cheap pivoting strategy utilizes error bounds based on more liberal interchange criteria, an idea borrowed from some work done in the area of sparse elimination methods, see:   
R. Fletcher (1976). “Factorizing Symmetric Indefinite Matrices,” Lin. Alg. Applic. 14, 257–272.   
Before using any symmetric Ax = b solver, it may be advisable to equilibrate A. An $O ( n ^ { 2 } )$ algorithm for accomplishing this task is given in:   
J.R. Bunch (1971). “Equilibration of Symmetric Matrices in the Max-Norm,”J. ACM 18, 566–572.   
N.J. Higham (1997). “Stability of the Diagonal Pivoting Method with Partial Pivoting,” SIAM J. Matrix Anal. Applic. 18, 52–65.   
Procedures for skew-symmetric systems similar to the methods that we have presented in this section also exist:   
J.R. Bunch (1982). “A Note on the Stable Decomposition of Skew Symmetric Matrices,” Math. Comput. 158, 475–480.   
J. Bunch (1982). “Stable Decomposition of Skew-Symmetric Matrices,” Math. Comput. 38, 475–479.   
P. Benner, R. Byers, H. Fassbender, V. Mehrmann, and D. Watkins (2000). “Cholesky-like Factorizations of Skew-Symmetric Matrices,” ETNA 11, 85–93.   
For a discussion of symmetric indefinite system solvers that are also banded or sparse, see:   
C. Ashcraft, R.G. Grimes, and J.G. Lewis (1998). “Accurate Symmetric Indefinite Linear Equation Solvers,” SIAM J. Matrix Anal. Applic. 20, 513–561.   
S.H. Cheng and N.J. Higham (1998). “A Modified Cholesky Algorithm Based on a Symmetric Indefinite Factorization,” SIAM J. Matrix Anal. Applic. 19, 1097–1110.   
J. Zhao, W. Wang, and W. Ren (2004). “Stability of the Matrix Factorization for Solving Block Tridiagonal Symmetric Indefinite Linear Systems,” BIT 44, 181–188.   
H. Fang and D.P. O’Leary (2006). “Stable Factorizations of Symmetric Tridiagonal and Triadic Matrices,” SIAM J. Matrix Anal. Applic. 28, 576–595.   
D. Irony and S. Toledo (2006). “The Snap-Back Pivoting Method for Symmetric Banded Indefinite Matrices,” SIAM J. Matrix Anal. Applic. 28, 398–424.

The equilibrium system literature is scattered among the several application areas where it has an important role to play. Nice overviews with pointers to this literature include:

G. Strang (1988). “A Framework for Equilibrium Equations,” SIAM Review 30, 283–297.

S.A. Vavasis (1994). “Stable Numerical Algorithms for Equilibrium Systems,” SIAM J. Matrix Anal. Applic. 15, 1108–1131.

P.E. Gill, M.A. Saunders, and J.R. Shinnerl (1996). “On the Stability of Cholesky Factorization for Symmetric Quasidefinite Systems,” SIAM J. Matrix Anal. Applic. 17, 35–46.

G.H. Golub and C. Greif (2003). “On Solving Block-Structured Indefinite Linear Systems,” SIAM J. Sci. Comput. 24, 2076–2092.

For a discussion of (4.4.21), see:

G.W. Stewart (1989). “On Scaled Projections and Pseudoinverses,” Lin. Alg. Applic. 112, 189–193.

D.P. O’Leary (1990). “On Bounds for Scaled Projections and Pseudoinverses,” Lin. Alg. Applic. 132, 115–117.

M.J. Todd (1990). “A Dantzig-Wolfe-like Variant of Karmarkar’s Interior-Point Linear Programming Algorithm,” Oper. Res. 38, 1006–1018.

An equilibrium system is a special case of a saddle point system. See §11.5.10.

# 4.5 Block Tridiagonal Systems

Block tridiagonal linear systems of the form

$$
\left[ \begin{array}{c c c c c} D _ {1} & F _ {1} & & \dots & 0 \\ E _ {1} & D _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & F _ {N - 1} \\ 0 & \dots & & E _ {N - 1} & D _ {N} \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \\ \vdots \\ \vdots \\ x _ {N} \end{array} \right] = \left[ \begin{array}{c} b _ {1} \\ b _ {2} \\ \vdots \\ \vdots \\ b _ {N} \end{array} \right]. \tag {4.5.1}
$$

frequently arise in practice. We assume for clarity that all blocks are $q { \mathrm { - b y - } } q$ . In this section we discuss both a block LU approach to this problem as well as a pair of divide-and-conquer schemes.

# 4.5.1 Block Tridiagonal LU Factorization

If

$$
A = \left[ \begin{array}{c c c c c} D _ {1} & F _ {1} & & \dots & 0 \\ E _ {1} & D _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & F _ {N - 1} \\ 0 & \dots & & E _ {N - 1} & D _ {N} \end{array} \right] \tag {4.5.2}
$$

then by comparing blocks in

$$
A = \left[ \begin{array}{c c c c c} I & & & \dots & 0 \\ L _ {1} & I & & & \vdots \\ & \ddots & \ddots & & \\ \vdots & & \ddots & & \\ 0 & \dots & & L _ {N - 1} & I \end{array} \right] \left[ \begin{array}{c c c c c} U _ {1} & F _ {1} & & \dots & 0 \\ 0 & U _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & F _ {N - 1} \\ 0 & \dots & & 0 & U _ {N} \end{array} \right] \tag {4.5.3}
$$

we formally obtain the following algorithm for computing the $L _ { i }$ and $U _ { i }$ :

$$
U _ {1} = D _ {1}
$$

for i = 2:N

$$
\text { Solve } L _ {i - 1} U _ {i - 1} = E _ {i - 1} \text { for } L _ {i - 1}. \tag {4.5.4}
$$

$$
U _ {i} = D _ {i} - L _ {i - 1} F _ {i - 1}
$$

end

The procedure is defined as long as the $U _ { i }$ are nonsingular.

Having computed the factorization (4.5.3), the vector x in (4.5.1) can be obtained via block forward elimination and block back substitution:

$$
y _ {1} = b _ {1}
$$

for i = 2:N

$$
y _ {i} = b _ {i} - L _ {i - 1} y _ {i - 1}
$$

end (4.5.5)

Solve UN xN = yN for $x _ { N }$

for i = N − 1: − 1:1

Solve Uixi = yi − Fixi+1 for xi

end

To carry out both (4.5.4) and (4.5.5), each $U _ { i }$ must be factored since linear systems involving these submatrices are solved. This could be done using Gaussian elimination with pivoting. However, this does not guarantee the stability of the overall process.

# 4.5.2 Block Diagonal Dominance

In order to obtain satisfactory bounds on the $L _ { i }$ and $U _ { i }$ it is necessary to make additional assumptions about the underlying block matrix. For example, if we have

$$
\left\| D _ {i} ^ {- 1} \right\| _ {1} \left(\left\| F _ {i - 1} \right\| _ {1} + \left\| E _ {i} \right\| _ {1}\right) <   1, \quad E _ {N} \equiv F _ {0} \equiv 0, \tag {4.5.6}
$$

for $i = 1 { : } N$ , then the factorization (4.5.3) exists and it is possible to show that the $L _ { i }$ and $U _ { i }$ satisfy the inequalities

$$
\left\| L _ {i} \right\| _ {1} \leq 1, \tag {4.5.7}
$$

$$
\left\| U _ {i} \right\| _ {1} \leq \left\| A _ {n} \right\| _ {1}. \tag {4.5.8}
$$

The conditions (4.5.6) define a type of block diagonal dominance.

# 4.5.3 Block-Cyclic Reduction

We next describe the method of block-cyclic reduction that can be used to solve some important special instances of the block tridiagonal system (4.5.1). For simplicity, we assume that A has the form

$$
A = \left[ \begin{array}{c c c c c} D & F & & \dots & 0 \\ F & D & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & F \\ 0 & \dots & & F & D \end{array} \right] \in \mathbb {R} ^ {N q \times N q} \tag {4.5.9}
$$

where $F$ and $D$ are $q { \mathrm { - b y - } } q$ matrices that satisfy ${ D F } = F { D }$ . We also assume that $N = 2 ^ { k } - 1$ . These conditions hold in certain important applications such as the discretization of Poisson’s equation on a rectangle. (See §4.8.4.)

The basic idea behind cyclic reduction is to halve repeatedly the dimension of the problem on hand repeatedly until we are left with a single q-by-q system for the unknown subvector $x _ { 2 } k - 1$ . This system is then solved by standard means. The previously eliminated $x _ { i }$ are found by a back-substitution process.

The general procedure is adequately illustrated by considering the case $N = 7 { : }$

$$
\begin{array}{l} b _ {1} = D x _ {1} + F x _ {2}, \\ b _ {2} = F x _ {1} + D x _ {2} + F x _ {3}, \\ b _ {3} = F x _ {2} + D x _ {3} + F x _ {4}, \\ b _ {4} = F x _ {3} + D x _ {4} + F x _ {5}, \\ b _ {5} = F x _ {4} + D x _ {5} + F x _ {6}, \\ b _ {6} = F x _ {5} + D x _ {6} + F x _ {7}, \\ b _ {7} = F x _ {6} + D x _ {7}. \\ \end{array}
$$

For i = 2, 4, and 6 we multiply equations $i - 1$ , i, and $i + 1$ by $F , \mathrm { ~ } - D _ { \mathrm { { i } } }$ , and $F _ { ; }$ , respectively, and add the resulting equations to obtain

$$
\begin{array}{l} (2 F ^ {2} - D ^ {2}) x _ {2} + F ^ {2} x _ {4} = F (b _ {1} + b _ {3}) - D b _ {2}, \\ F ^ {2} x _ {2} + (2 F ^ {2} - D ^ {2}) x _ {4} + F ^ {2} x _ {6} = F (b _ {3} + b _ {5}) - D b _ {4}, \\ F ^ {2} x _ {4} + (2 F ^ {2} - D ^ {2}) x _ {6} = F (b _ {5} + b _ {7}) - D b _ {6}. \\ \end{array}
$$

Thus, with this tactic we have removed the odd-indexed $x _ { i }$ and are left with a reduced block tridiagonal system of the form

$$
\begin{array}{l} D ^ {(1)} x _ {2} + F ^ {(1)} x _ {4} = b _ {2} ^ {(1)}, \\ F ^ {(1)} x _ {2} + D ^ {(1)} x _ {4} + F ^ {(1)} x _ {6} = b _ {4} ^ {(1)}, \\ F ^ {(1)} x _ {4} + D ^ {(1)} x _ {6} = b _ {6} ^ {(1)}, \\ \end{array}
$$

where $D ^ { ( 1 ) } = 2 F ^ { 2 } - D ^ { 2 }$ and $F ^ { ( 1 ) } = F ^ { 2 }$ commute. Applying the same elimination strategy as above, we multiply these three equations respectively by $F ^ { ( 1 ) } , - D ^ { ( 1 ) }$ , and $F ^ { ( 1 ) }$ . When these transformed equations are added together, we obtain the single equation

$$
\left(2 [ F ^ {(1)} ] ^ {2} - D ^ {(1) ^ {2}}\right) x _ {4} = F ^ {(1)} \left(b _ {2} ^ {(1)} + b _ {6} ^ {(1)}\right) - D ^ {(1)} b _ {4} ^ {(1)},
$$

which we write as

$$
D ^ {(2)} x _ {4} = b ^ {(2)}.
$$

This completes the cyclic reduction. We now solve this (small) q-by-q system for $x _ { 4 }$ . The vectors $x _ { 2 }$ and $x _ { 6 }$ are then found by solving the systems

$$
\begin{array}{l} D ^ {(1)} x _ {2} = b _ {2} ^ {(1)} - F ^ {(1)} x _ {4}, \\ D ^ {(1)} x _ {6} = b _ {6} ^ {(1)} - F ^ {(1)} x _ {4}. \\ \end{array}
$$

Finally, we use the first, third, fifth, and seventh equations in the original system to compute $x _ { 1 } , x _ { 3 } , x _ { 5 }$ , and $x _ { 7 }$ , respectively.

The amount of work required to perform these recursions for general N depends greatly upon the sparsity of the $D ^ { ( p ) }$ and $F ^ { ( p ) }$ . In the worst case when these matrices are full, the overall flop count has order $\mathrm { l o g } ( N ) q ^ { 3 }$ . Care must be exercised in order to ensure stability during the reduction. For further details, see Buneman (1969).

# 4.5.4 The SPIKE Framework

A bandwidth-p matrix $A \in \mathbb { R } ^ { N q \times N q }$ can also be regarded as a block tridiagonal matrix with banded diagonal blocks and low-rank off-diagonal blocks. Here is an example where $N = 4 , q = 7$ , and $p = 2$ :

$$
A = \left[ \begin{array}{c c c c c c c c c} \times & \times & \times \\ \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times
$$

Note that the diagonal blocks have bandwidth p and the blocks along the subdiagonal and superdiagonal have rank p. The low rank of the off-diagonal blocks makes it possible to formulate a divide-and-conquer procedure known as the “SPIKE” algorithm. The method is of interest because it parallelizes nicely. Our brief discussion is based on Polizzi and Sameh (2007).

Assume for clarity that the diagonal blocks $D _ { 1 } , \ldots , D _ { 4 }$ are sufficiently well conditioned. If we premultiply the above matrix by the inverse of $\mathrm { d i a g } ( D _ { 1 } , D _ { 2 } , D _ { 3 } , D _ { 4 } )$ , then we obtain

$$
\tilde {A} = \left[ \begin{array}c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
$$

With this maneuver, the original linear system

$$
\left[ \begin{array}{c c c c} D _ {1} & F _ {1} & 0 & 0 \\ E _ {1} & D _ {2} & F _ {2} & 0 \\ 0 & E _ {2} & D _ {3} & F _ {3} \\ 0 & 0 & E _ {3} & D _ {4} \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \\ x _ {3} \\ x _ {4} \end{array} \right] = \left[ \begin{array}{c} b _ {1} \\ b _ {2} \\ b _ {3} \\ b _ {4} \end{array} \right], \tag {4.5.13}
$$

which corresponds to (4.5.11), transforms to

$$
\left[ \begin{array}{c c c c} I _ {7} & \tilde {F} _ {1} & 0 & 0 \\ \tilde {E} _ {1} & I _ {7} & \tilde {F} _ {2} & 0 \\ 0 & \tilde {E} _ {2} & I _ {7} & \tilde {F} _ {3} \\ 0 & 0 & \tilde {E} _ {3} & I _ {7} \end{array} \right] \left[ \begin{array}{l} x _ {1} \\ x _ {2} \\ x _ {3} \\ x _ {4} \end{array} \right] = \left[ \begin{array}{l} \tilde {b} _ {1} \\ \tilde {b} _ {2} \\ \tilde {b} _ {3} \\ \tilde {b} _ {4} \end{array} \right], \tag {4.5.14}
$$

where $D _ { i } \tilde { b } _ { i } = b _ { i } , D _ { i } \tilde { F } _ { i } = F _ { i }$ , and $D _ { i + 1 } { \tilde { E } } _ { i } = E _ { i }$ . Next, we refine the blocking (4.5.14) by turning each submatrix into a 3-by-3 block matrix and each subvector into a 3-by-1 block vector as follows:

$$
\left[ \begin{array}{c c c c c c c c c c c c} I _ {2} & 0 & 0 & K _ {1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & I _ {3} & 0 & H _ {1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & I _ {2} & G _ {1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & R _ {1} & I _ {2} & 0 & 0 & K _ {2} & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & S _ {1} & 0 & I _ {3} & 0 & H _ {2} & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & T _ {1} & 0 & 0 & I _ {2} & G _ {2} & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & R _ {2} & I _ {2} & 0 & 0 & K _ {3} & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & S _ {2} & 0 & I _ {3} & 0 & H _ {3} & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & T _ {2} & 0 & 0 & I _ {2} & G _ {3} & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & R _ {3} & I _ {q} & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & S _ {3} & 0 & I _ {m} & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & T _ {3} & 0 & 0 & I _ {q} \end{array} \right] \left[ \begin{array}{c} w _ {1} \\ y _ {1} \\ z _ {1} \\ \hline w _ {2} \\ y _ {2} \\ z _ {2} \\ \hline w _ {3} \\ y _ {3} \\ z _ {3} \\ \hline w _ {4} \\ y _ {4} \\ z _ {4} \end{array} \right] = \left[ \begin{array}{c} c _ {1} \\ d _ {1} \\ f _ {1} \\ \hline c _ {2} \\ d _ {2} \\ f _ {2} \\ \hline c _ {3} \\ d _ {3} \\ f _ {3} \\ \hline c _ {4} \\ d _ {4} \\ f _ {4} \end{array} \right]. \quad (4. 5. 1 5)
$$

The block rows and columns in this equation can be reordered to produce the following equivalent system:

$$
\left[ \begin{array}{c c c c c c c c c c c c} I _ {2} & 0 & K _ {1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & I _ {2} & G _ {1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & R _ {1} & I _ {2} & 0 & K _ {2} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & T _ {1} & 0 & I _ {2} & G _ {2} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & R _ {2} & I _ {2} & 0 & K _ {3} & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & T _ {2} & 0 & I _ {2} & G _ {3} & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & R _ {3} & I _ {2} & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & T _ {3} & 0 & I _ {2} & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & H _ {1} & 0 & 0 & 0 & 0 & 0 & I _ {3} & 0 & 0 & 0 \\ 0 & S _ {1} & 0 & 0 & H _ {2} & 0 & 0 & 0 & 0 & I _ {3} & 0 & 0 \\ 0 & 0 & 0 & S _ {2} & 0 & 0 & H _ {3} & 0 & 0 & 0 & I _ {3} & 0 \\ 0 & 0 & 0 & 0 & 0 & S _ {3} & 0 & 0 & 0 & 0 & 0 & I _ {3} \end{array} \right] \left[ \begin{array}{c} w _ {1} \\ z _ {1} \\ w _ {2} \\ z _ {2} \\ w _ {3} \\ z _ {3} \\ w _ {4} \\ z _ {4} \\ \hline y _ {1} \\ y _ {2} \\ y _ {3} \\ y _ {4} \end{array} \right] = \left[ \begin{array}{c} c _ {1} \\ f _ {1} \\ c _ {2} \\ f _ {2} \\ c _ {3} \\ f _ {3} \\ c _ {4} \\ f _ {4} \\ \hline d _ {1} \\ d _ {2} \\ d _ {3} \\ d _ {4} \end{array} \right]. \tag {4.5.16}
$$

If we assume that $N \gg 1$ , then the (1,1) block is a relatively small banded matrix that define the $z _ { i }$ and $w _ { i }$ . Once these quantities are computed, then the remaining unknowns follow from a decoupled set of large matrix-vector multiplications, $\mathrm { e . g . } , y _ { 1 } = d _ { 1 } - H _ { 1 } w _ { 2 }$ , $y _ { 2 } = d _ { 2 } - S _ { 1 } z _ { 1 } - H _ { 2 } w _ { 3 } , y _ { 3 } = d _ { 3 } - S _ { 2 } z _ { 2 } - H _ { 3 } w _ { 4 }$ , and $y _ { 4 } = d _ { 4 } - S _ { 3 } z _ { 3 }$ . Thus, in a fourprocessor execution of this method, there are (short) communications that involves the $w _ { i }$ and $z _ { i }$ and a lot of large, local gaxpy computations.

# Problems

P4.5.1 (a) Show that a block diagonally dominant matrix is nonsingular. (b) Verify that (4.5.6) implies (4.5.7) and (4.5.8).

P4.5.2 Write a recursive function $\boldsymbol { x } = \mathsf { C R } ( D , F , N , b )$ that returns the solution to $A x = b$ where A is specified by (4.5.9). Assume that $N = 2 ^ { k } - 1$ for some positive integer k, D, $F \in \mathbb { R } ^ { q \times q } ,$ and $b \in \mathbb { R } ^ { N q }$ .

P4.5.3 How would you solve a system of the form

$$
\left[ \begin{array}{c c} D _ {1} & F _ {1} \\ E _ {1} & D _ {2} \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right] = \left[ \begin{array}{c} b _ {1} \\ b _ {2} \end{array} \right]
$$

where $D _ { 1 }$ and $D _ { 2 }$ are diagonal and $F _ { 1 }$ and $E _ { 1 }$ are tridiagonal? Hint: Use the perfect shuffle permutation.

P4.5.4 In the simplified SPIKE framework that we presented in §4.5.4, we treat A as an N-by-N block matrix with q-by-q blocks. It is assumed that $A \in \mathbb { R } ^ { N q \times N q }$ has bandwidth p and that $p \ll q$ . For this general case, describe the block sizes that result when the transition from (4.5.11) to (4.5.16) is carried out. Assuming that A’s band is dense, what fraction of flops are gaxpy flops?

# Notes and References for 4.5

The following papers provide insight into the various nuances of block matrix computations:

J.M. Varah (1972). “On the Solution of Block-Tridiagonal Systems Arising from Certain Finite-Difference Equations,” Math. Comput. 26, 859–868.   
R. Fourer (1984). “Staircase Matrices and Systems,” SIAM Review 26, 1–71.   
M.L. Merriam (1985). “On the Factorization of Block Tridiagonals With Storage Constraints,” SIAM J. Sci. Stat. Comput. 6, 182-192.   
The property of block diagonal dominance and its various implications is the central theme in:   
D.G. Feingold and R.S. Varga (1962). “Block Diagonally Dominant Matrices and Generalizations of the Gershgorin Circle Theorem,”Pacific J. Math. 12, 1241–1250.   
R.S. Varga (1976). “On Diagonal Dominance Arguments for Bounding $\| \ A ^ { - 1 } \ \| _ { \infty } , \stackrel {  } { _ {  } }$ Lin. Alg. Applic. 14, 211–217.   
Early methods that involve the idea of cyclic reduction are described in:   
R.W. Hockney (1965). “A Fast Direct Solution of Poisson’s Equation Using Fourier Analysis, ”J. ACM 12, 95–113.   
B.L. Buzbee, G.H. Golub, and C.W. Nielson (1970). “On Direct Methods for Solving Poisson’s Equations,” SIAM J. Numer. Anal. 7, 627–656.   
The accumulation of the right-hand side must be done with great care, for otherwise there would be a significant loss of accuracy. A stable way of doing this is described in:   
O. Buneman (1969). “A Compact Non-Iterative Poisson Solver,” Report 294, Stanford University Institute for Plasma Research, Stanford, CA.   
Other literature concerned with cyclic reduction includes:   
F.W. Dorr (1970). “The Direct Solution of the Discrete Poisson Equation on a Rectangle,” SIAM Review 12, 248–263.

B.L. Buzbee, F.W. Dorr, J.A. George, and G.H. Golub (1971). “The Direct Solution of the Discrete Poisson Equation on Irregular Regions,”SIAM J. Numer. Anal. 8, 722–736.   
F.W. Dorr (1973). “The Direct Solution of the Discrete Poisson Equation in $O ( n ^ { 2 } )$ Operations,” SIAM Review 15, 412–415.   
P. Concus and G.H. Golub (1973). “Use of Fast Direct Methods for the Efficient Numerical Solution of Nonseparable Elliptic Equations,” SIAM J. Numer. Anal. 10, 1103–1120.   
B.L. Buzbee and F.W. Dorr (1974). “The Direct Solution of the Biharmonic Equation on Rectangular Regions and the Poisson Equation on Irregular Regions,”SIAM J. Numer. Anal. 11, 753–763.   
D. Heller (1976). “Some Aspects of the Cyclic Reduction Algorithm for Block Tridiagonal Linear Systems,”SIAM J. Numer. Anal. 13, 484–496.

Various generalizations and extensions to cyclic reduction have been proposed:

P.N. Swarztrauber and R.A. Sweet (1973). “The Direct Solution of the Discrete Poisson Equation on a Disk,” SIAM J. Numer. Anal. 10, 900–907.   
R.A. Sweet (1974). “A Generalized Cyclic Reduction Algorithm,” SIAM J. Num. Anal. 11, 506–20.   
M.A. Diamond and D.L.V. Ferreira (1976). “On a Cyclic Reduction Method for the Solution of Poisson’s Equation,” SIAM J. Numer. Anal. 13, 54–70.   
R.A. Sweet (1977). “A Cyclic Reduction Algorithm for Solving Block Tridiagonal Systems of Arbitrary Dimension,” SIAM J. Numer. Anal. 14, 706–720.   
P.N. Swarztrauber and R. Sweet (1989). “Vector and Parallel Methods for the Direct Solution of Poisson’s Equation,” J. Comput. Appl. Math. 27, 241–263.   
S. Bondeli and W. Gander (1994). “Cyclic Reduction for Special Tridiagonal Systems,” SIAM J. Matrix Anal. Applic. 15, 321–330.

A 2-by-2 block system with very thin (1,2) and (2,1) blocks is referred to as a bordered linear system. Special techniques for problems with this structure are discussed in:

W. Govaerts and J.D. Pryce (1990). “Block Elimination with One Iterative Refinement Solves Bordered Linear Systems Accurately,” BIT 30, 490–507.   
W. Govaerts (1991). “Stable Solvers and Block Elimination for Bordered Systems,” SIAM J. Matrix Anal. Applic. 12, 469–483.   
W. Govaerts and J.D. Pryce (1993). “Mixed Block Elimination for Linear Systems with Wider Borders,” IMA J. Numer. Anal. 13, 161–180.

Systems that are block bidiagonal, block Hessenberg, and block triangular also occur, see:

G. Fairweather and I. Gladwell (2004). “Algorithms for Almost Block Diagonal Linear Systems,” SIAM Review 46, 49–58.   
U. von Matt and G. W. Stewart (1996). “Rounding Errors in Solving Block Hessenberg Systems,” Math. Comput. 65, 115–135.   
L. Gemignani and G. Lotti (2003). “Efficient and Stable Solution of M-Matrix Linear Systems of (Block) Hessenberg Form,” SIAM J. Matrix Anal. Applic. 24, 852–876.   
M. Hegland and M.R. Osborne (1998). “Wrap-Around Partitioning for Block Bidiagonal Linear Systems,” IMA J. Numer. Anal. 18, 373–383.   
T. Rossi and J. Toivanen (1999). “A Parallel Fast Direct Solver for Block Tridiagonal Systems with Separable Matrices of Arbitrary Dimension,” SIAM J. Sci. Comput. 20, 1778–1793.   
I.M. Spitkovsky and D. Yong (2000). “Almost Periodic Factorization of Certain Block Triangular Matrix Functions,” Math. Comput. 69, 1053–1070.

The SPIKE framework supports many different options according to whether the band is sparse or dense. Also, steps have to be taken if the diagonal blocks are ill-conditioned, see:

E. Polizzi and A. Sameh (2007). “SPIKE: A Parallel Environment for Solving Banded Linear Systems,” Comput. Fluids 36, 113–120.   
C.C.K. Mikkelsen and M. Manguoglu (2008). “Analysis of the Truncated SPIKE Algorithm,” SIAM J. Matrix Anal. Applic. 30, 1500–1519.

# 4.6 Vandermonde Systems

Suppose x $( 0 : n ) \in \mathbb { R } ^ { n + 1 }$ . A matrix $V \in \mathbb { R } ^ { ( n + 1 ) \times ( n + 1 ) }$ of the form

$$
V = V (x _ {0}, \ldots , x _ {n}) = \left[ \begin{array}{c c c c} 1 & 1 & \dots & 1 \\ x _ {0} & x _ {1} & \dots & x _ {n} \\ \vdots & \vdots & & \vdots \\ x _ {0} ^ {n} & x _ {1} ^ {n} & \dots & x _ {n} ^ {n} \end{array} \right]
$$

is said to be a Vandermonde matrix. Note that the discrete Fourier transform matrix (§1.4.1) is a very special complex Vandermonde matrix.

In this section, we show how the systems $V ^ { T } a = f = f ( 0 { : } n )$ and $V z = b = b ( 0 { : } n )$ can be solved in $O ( n ^ { 2 } )$ flops. For convenience, vectors and matrices are subscripted from 0 in this section.

# 4.6.1 Polynomial Interpolation: $V ^ { T } a = f$

Vandermonde systems arise in many approximation and interpolation problems. Indeed, the key to obtaining a fast Vandermonde solver is to recognize that solving $V ^ { T } a = f$ is equivalent to polynomial interpolation. This follows because if $V ^ { T } a = f$ and

$$
p (x) = \sum_ {j = 0} ^ {n} a _ {j} x ^ {j}, \tag {4.6.1}
$$

then $p ( x _ { i } ) = f _ { i }$ for $i = 0 { : } n$ .

Recall that if the $x _ { i }$ are distinct then there is a unique polynomial of degree n that interpolates $( x _ { 0 } , f _ { 0 } ) , \ldots , ( x _ { n } , f _ { n } )$ . Consequently, V is nonsingular as long as the $x _ { i }$ are distinct. We assume this throughout the section.

The first step in computing the $a _ { j }$ of (4.6.1) is to calculate the Newton representation of the interpolating polynomial $p \mathrm { : }$

$$
p (x) = \sum_ {k = 0} ^ {n} c _ {k} \left(\prod_ {i = 0} ^ {k - 1} (x - x _ {i})\right). \tag {4.6.2}
$$

The constants $c _ { k }$ are divided differences and may be determined as follows:

$$
c (0: n) = f (0: n)
$$

$$
\text { for } k = 0: n - 1
$$

$$
\text { for   } i = n: - 1: k + 1 \tag {4.6.3}
$$

$$
c _ {i} = \left(c _ {i} - c _ {i - 1}\right) / \left(x _ {i} - x _ {i - k - 1}\right)
$$

end

end

See Conte and deBoor (1980).

The next task is to generate the coefficients $a _ { 0 } , \ldots , a _ { n }$ in (4.6.1) from the Newton representation coefficients $c _ { 0 } , \ldots , c _ { n }$ . Define the polynomials $p _ { n } ( x ) , \ldots , p _ { 0 } ( x )$ by the iteration

$$
p _ {n} (x) = c _ {n}
$$

for k = n − 1 : −1 : 0

$$
p _ {k} (x) = c _ {k} + (x - x _ {k}) \cdot p _ {k + 1} (x)
$$

end

and observe that $p _ { 0 } ( x ) = p ( x )$ . Writing

$$
p _ {k} (x) = a _ {k} ^ {(k)} + a _ {k + 1} ^ {(k)} x + \dots + a _ {n} ^ {(k)} x ^ {n - k}
$$

and equating like powers of x in the equation $p _ { k } = c _ { k } + ( x - x _ { k } ) p _ { k + 1 }$ gives the following recursion for the coefficients $a _ { i } ^ { ( k ) }$ (k)

$$
a _ {n} ^ {(n)} = c _ {n}
$$

for k = n−1: −1 : 0

$$
a _ {k} ^ {(k)} = c _ {k} - x _ {k} a _ {k + 1} ^ {(k + 1)}
$$

for i = k + 1 : n − 1

$$
a _ {i} ^ {(k)} = a _ {i} ^ {(k + 1)} - x _ {k} a _ {i + 1} ^ {(k + 1)}
$$

end

$$
a _ {n} ^ {(k)} = a _ {n} ^ {(k + 1)}
$$

end

Consequently, the coefficients $a _ { i } = a _ { i } ^ { ( 0 ) }$ can be calculated as follows:

$$
a (0: n) = c (0: n)
$$

for k = n−1: −1 : 0

for i = k:n − 1 (4.6.4)

$$
a _ {i} = a _ {i} - x _ {k} a _ {i + 1}
$$

end

end

Combining this iteration with (4.6.3) gives the following algorithm.

Algorithm 4.6.1 Given $x ( 0 : n ) \in \mathbb { R } ^ { n + 1 }$ with distinct entries and $f = f ( 0 : n ) \in \mathbb { R } ^ { n + 1 }$ , the following algorithm overwrites f with the solution $a = a ( 0 : n )$ to the Vandermonde system $V ( x _ { 0 } , \ldots , x _ { n } ) ^ { T } a = f .$ .

for k = 0 : n − 1

for i = n: −1 :k + 1

$$
f (i) = (f (i) - f (i - 1)) / (x (i) - x (i - k - 1))
$$

end

end

for k = n − 1: −1 : 0

for i = k : n − 1

$$
f (i) = f (i) - f (i + 1) \cdot x (k)
$$

end

end

This algorithm requires $5 n ^ { 2 } / 2$ flops.

# 4.6.2 The System $V z = b$

Now consider the system $V z = b$ . To derive an efficient algorithm for this problem, we describe what Algorithm 4.6.1 does in matrix-vector language. Define the lower bidiagonal matrix $L _ { k } ( \alpha ) \in \mathbb { R } ^ { ( n + 1 ) \times ( n + 1 ) }$ by

$$
L _ {k} (\alpha) = \left[ \begin{array}{c c c c c c c} I _ {k} & & & 0 & & & \\ \hline & 1 & 0 & \dots & & & 0 \\ & - \alpha & 1 & & & & \\ & 0 & \ddots & \ddots & & & \\ 0 & \vdots & & \ddots & \ddots & & \vdots \\ & & & & \ddots & 1 & \\ & 0 & & \dots & & - \alpha & 1 \end{array} \right]
$$

and the diagonal matrix $D _ { k }$ by

$$
D _ {k} = \operatorname{diag} (\underbrace {1 , \ldots , 1} _ {k + 1}, x _ {k + 1} - x _ {0}, \ldots , x _ {n} - x _ {n - k - 1}).
$$

With these definitions it is easy to verify from (4.6.3) that, if $f = f ( 0 : n ) \operatorname { a n d } c = c ( 0 : n )$ is the vector of divided differences, then

$$
c = U ^ {T} f
$$

where U is the upper triangular matrix defined by

$$
U ^ {T} = D _ {n - 1} ^ {- 1} L _ {n - 1} (1) \dots D _ {0} ^ {- 1} L _ {0} (1).
$$

Similarly, from (4.6.4) we have

$$
a = L ^ {T} c,
$$

where L is the unit lower triangular matrix defined by

$$
L ^ {T} = L _ {0} (x _ {0}) ^ {T} \dots L _ {n - 1} (x _ {n - 1}) ^ {T}.
$$

It follows that $a = V ^ { - T } f$ is given by

$$
a = L ^ {T} U ^ {T} f.
$$

Thus,

$$
V ^ {- T} = L ^ {T} U ^ {T}
$$

which shows that Algorithm 4.6.1 solves $V ^ { T } a ~ = ~ f$ by tacitly computing the $^ { 6 6 } \mathrm { U L }$ factorization” of $V ^ { - 1 }$ . Consequently, the solution to the system $V z = b$ is given by

$$
\begin{array}{l} z = V ^ {- 1} b = U (L b) \\ = \left(L _ {0} (1) ^ {T} D _ {0} ^ {- 1} \dots L _ {n - 1} (1) ^ {T} D _ {n - 1} ^ {- 1}\right) \left(L _ {n - 1} (x _ {n - 1}) \dots L _ {0} (x _ {0}) b\right). \\ \end{array}
$$

This observation gives rise to the following algorithm:

Algorithm 4.6.2 Given $x ( 0 : n ) \in \mathbb { R } ^ { n + 1 }$ with distinct entries and $b = b ( 0 : n ) \in \mathbb { R } ^ { n + 1 }$ , the following algorithm overwrites b with the solution $z = z ( 0 : n )$ to the Vandermonde system $V ( x _ { 0 } , \ldots , x _ { n } ) z = b .$

for k = 0:n - 1
    for i = n: -1:k + 1
    b(i) = b(i) - x(k)b(i - 1)
    end
end
for k = n - 1: -1: 0
    for i = k + 1:n
    b(i) = b(i)/(x(i) - x(i - k - 1))
    end
    for i = k:n - 1
    b(i) = b(i) - b(i + 1)
    end
end

This algorithm requires $5 n ^ { 2 } / 2$ flops.

Algorithms 4.6.1 and 4.6.2 are discussed and analyzed by Bj¨orck and Pereyra (1970). Their experience is that these algorithms frequently produce surprisingly accurate solutions, even if V is ill-conditioned.

We mention that related techniques have been developed and analyzed for confluent Vandermonde systems, e.g., systems of the form

$$
\left[ \begin{array}{c c c c} 1 & 1 & 0 & 1 \\ x _ {0} & x _ {1} & 1 & x _ {3} \\ x _ {0} ^ {2} & x _ {1} ^ {2} & 2 x _ {1} & x _ {3} ^ {2} \\ x _ {0} ^ {3} & x _ {1} ^ {3} & 3 x _ {1} ^ {2} & x _ {3} ^ {3} \end{array} \right] ^ {T} \left[ \begin{array}{c} a _ {0} \\ a _ {1} \\ a _ {2} \\ a _ {3} \end{array} \right] = \left[ \begin{array}{c} f _ {0} \\ f _ {1} \\ f _ {2} \\ f _ {3} \end{array} \right].
$$

See Higham (1990).

# Problems

P4.6.1 Show that if $V = V ( x _ { 0 } , \ldots , x _ { n } )$ , then

$$
\det (V) = \prod_ {n \geq i > j \geq 0} (x _ {i} - x _ {j}).
$$

P4.6.2 (Gautschi 1975) Verify the following inequality for the n = 1 case above:

$$
\| V^{-1}\|_{\infty}\leq \max_{0\leq k\leq n}\prod_{\substack{i = 0\\ i\neq k}}^{n}\frac{1 + |x_{i}|}{|x_{k} - x_{i}|}  .
$$

Equality results if the $x _ { i }$ are all on the same ray in the complex plane.

# Notes and References for §4.6

Our discussion of Vandermonde linear systems is drawn from the following papers:

The divided difference computations we discussed are detailed in:   
Error analyses of Vandermonde system solvers include:   
Interesting theoretical results concerning the condition of Vandermonde systems may be found in:   
A. Bj¨orck and V. Pereyra (1970). “Solution of Vandermonde Systems of Equations,”Math. Comput. 24, 893–903.   
A. Bj¨orck and T. Elfving (1973). “Algorithms for Confluent Vandermonde Systems,” Numer. Math. 21, 130–37.   
S.D. Conte and C. de Boor (1980). Elementary Numerical Analysis: An Algorithmic Approach, Third Edition, McGraw-Hill, New York, Chapter 2.   
N.J. Higham (1987). “Error Analysis of the Bj¨orck-Pereyra Algorithms for Solving Vandermonde Systems,” Numer. Math. 50, 613–632.   
N.J. Higham (1988). “Fast Solution of Vandermonde-like Systems Involving Orthogonal Polynomials,” IMA J. Numer. Anal. 8, 473–486.   
N.J. Higham (1990). “Stability Analysis of Algorithms for Solving Confluent Vandermonde-like Systems,” SIAM J. Matrix Anal. Applic. 11, 23–41.   
S.G. Bartels and D.J. Higham (1992). “The Structured Sensitivity of Vandermonde-Like Systems,” Numer. Math. 62, 17–34.   
J.M. Varah (1993). “Errors and Perturbations in Vandermonde Systems,” IMA J. Numer. Anal. 13, 1–12.   
W. Gautschi (1975). “Norm Estimates for Inverses of Vandermonde Matrices,”Numer. Math. 23, 337–347.   
W. Gautschi (1975). “Optimally Conditioned Vandermonde Matrices,” Numer. Math. 24, 1–12.   
J-G. Sun (1998). “Bounds for the Structured Backward Errors of Vandermonde Systems,” SIAM J. Matrix Anal. Applic. 20, 45–59.   
B.K. Alpert (1996). “Condition Number of a Vandermonde Matrix,” SIAM Review 38, 314–314.   
B. Beckermann (2000). “The condition number of real Vandermonde, Krylov and positive definite Hankel matrices,” Numer. Math. 85, 553–577.   
The basic algorithms presented can be extended to cover confluent Vandermonde systems, block Vandermonde systems, and Vandermonde systems with other polynomial bases:   
G. Galimberti and V. Pereyra (1970). “Numerical Differentiation and the Solution of Multidimensional Vandermonde Systems,” Math. Comput. 24, 357–364.   
G. Galimberti and V. Pereyra (1971). “Solving Confluent Vandermonde Systems of Hermitian Type,” Numer. Math. 18, 44–60.   
H. Van de Vel (1977). “Numerical Treatment of Generalized Vandermonde Systems of Equations,” Lin. Alg. Applic. 17, 149–174.   
G.H. Golub and W.P Tang (1981). “The Block Decomposition of a Vandermonde Matrix and Its Applications,”BIT 21, 505–517.   
D. Calvetti and L. Reichel (1992). “A Chebychev-Vandermonde Solver,” Lin. Alg. Applic. 172, 219–229.   
D. Calvetti and L. Reichel (1993). “Fast Inversion of Vandermonde-Like Matrices Involving Orthogonal Polynomials,” BIT 33, 473–484.   
H. Lu (1994). “Fast Solution of Confluent Vandermonde Linear Systems,” SIAM J. Matrix Anal. Applic. 15, 1277–1289.   
H. Lu (1996). “Solution of Vandermonde-like Systems and Confluent Vandermonde-like Systems,” SIAM J. Matrix Anal. Applic. 17, 127–138.   
M.-R. Skrzipek (2004). “Inversion of Vandermonde-Like Matrices,” BIT 44, 291–306.   
J.W. Demmel and P. Koev (2005). “The Accurate and Efficient Solution of a Totally Positive Generalized Vandermonde Linear System,” SIAM J. Matrix Anal. Applic. 27, 142–152.   
The displacement rank idea that we discuss in §12.1 can also be used to develop fast methods for Vandermonde systems.

# 4.7 Classical Methods for Toeplitz Systems

Matrices whose entries are constant along each diagonal arise in many applications and are called Toeplitz matrices. Formally, $T \in \mathbb { R } ^ { n \times n }$ is Toeplitz if there exist scalars $r _ { - n + 1 } , \ldots , r _ { 0 } , \ldots , r _ { n - 1 }$ such that $a _ { i j } = r _ { j - i }$ for all i and j. Thus,

$$
T = \left[ \begin{array}{l l l l} r _ {0} & r _ {1} & r _ {2} & r _ {3} \\ r _ {- 1} & r _ {0} & r _ {1} & r _ {2} \\ r _ {- 2} & r _ {- 1} & r _ {0} & r _ {1} \\ r _ {- 3} & r _ {- 2} & r _ {- 1} & r _ {0} \end{array} \right] = \left[ \begin{array}{l l l l} 3 & 1 & 7 & 6 \\ 4 & 3 & 1 & 7 \\ 0 & 4 & 3 & 1 \\ 9 & 0 & 4 & 3 \end{array} \right]
$$

is Toeplitz. In this section we show that Toeplitz systems can be solved in $O ( n ^ { 2 } )$ flops The discussion focuses on the important case when T is also symmetric and positive definite, but we also include a few comments about general Toeplitz systems. An alternative approach to Toeplitz system solving based on displacement rank is given in §12.1.

# 4.7.1 Persymmetry

The key fact that makes it possible to solve a Toeplitz system $T x = b$ so fast has to do with the structure of $T ^ { - 1 }$ . Toeplitz matrices belong to the larger class of persymmetric matrices. We say that $B \in \mathbb { R } ^ { n \times n }$ is persymmetric if

$$
\mathcal {E} _ {n} B \mathcal {E} _ {n} = B ^ {T}
$$

where ${ \mathcal { E } } _ { n }$ is the n-by-n exchange matrix defined in §1.2.11, e.g.,

$$
\mathcal {E} _ {4} = \left[ \begin{array}{c c c c} 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \end{array} \right].
$$

If B is persymmetric, then ${ \mathcal { E } } _ { n } B$ is symmetric. This means that B is symmetric about its antidiagonal. Note that the inverse of a persymmetric matrix is also persymmetric:

$$
\mathcal {E} _ {n} B ^ {- 1} \mathcal {E} _ {n} = (\mathcal {E} _ {n} B \mathcal {E} _ {n}) ^ {- 1} = (B ^ {T}) ^ {- 1} = (B ^ {- 1}) ^ {T}.
$$

Thus, the inverse of a nonsingular Toeplitz matrix is persymmetric.

# 4.7.2 Three Problems

Assume that we have scalars $r _ { 1 } , \ldots , r _ { n }$ such that for $k = 1 { : } n$ the matrices

$$
T _ {k} = \left[ \begin{array}{c c c c c} 1 & r _ {1} & \dots & r _ {k - 2} & r _ {k - 1} \\ r _ {1} & 1 & \ddots & & r _ {k - 2} \\ \vdots & \ddots & \ddots & \ddots & \vdots \\ r _ {k - 2} & & \ddots & \ddots & r _ {1} \\ r _ {k - 1} & r _ {k - 2} & \dots & r _ {1} & 1 \end{array} \right]
$$

are positive definite. (There is no loss of generality in normalizing the diagonal.) We set out to describe three important algorithms:

• Durbin’s algorithm for the Yule-Walker problem $T _ { n } y = - [ r _ { 1 } , \ldots , r _ { n } ] ^ { T }$   
• Levinson’s algorithm for the general right-hand-side problem $T _ { n } x = b$   
• Trench’s algorithm for computing $B = T _ { n } ^ { - 1 }$

# 4.7.3 Solving the Yule-Walker Equations

We begin by presenting Durbin’s algorithm for the Yule-Walker equations which arise in conjunction with certain linear prediction problems. Suppose for some k that satisfies $1 \leq k \leq n - 1$ we have solved the kth order Yule-Walker system $T _ { k } y = - r =$ $- [ r _ { 1 } , \ldots , r _ { k } ] ^ { T }$ . We now show how the $( k + 1 ) \mathrm { s t }$ order Yule-Walker system

$$
\left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ r ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] \left[ \begin{array}{c} z \\ \alpha \end{array} \right] = - \left[ \begin{array}{c} r \\ r _ {k + 1} \end{array} \right]
$$

can be solved in O(k) flops. First observe that

$$
z = T _ {k} ^ {- 1} (- r - \alpha \mathcal {E} _ {k} r) = y - \alpha T _ {k} ^ {- 1} \mathcal {E} _ {k} r
$$

and

$$
\alpha = - r _ {k + 1} - r ^ {T} \mathcal {E} _ {k} z.
$$

Since $T _ { k } ^ { - 1 }$ is persymmetric, $T _ { k } ^ { - 1 } { \mathcal { E } } _ { k } = { \mathcal { E } } _ { k } T _ { k } ^ { - 1 }$ and thus

$$
z = y - \alpha \mathcal {E} _ {k} T _ {k} ^ {- 1} r = y + \alpha \mathcal {E} _ {k} y.
$$

By substituting this into the above expression for α we find

$$
\alpha = - r _ {k + 1} - r ^ {T} \mathcal {E} _ {k} (y + \alpha \mathcal {E} _ {k} y) = - (r _ {k + 1} + r ^ {T} \mathcal {E} _ {k} y) / (1 + r ^ {T} y).
$$

The denominator is positive because $T _ { k + 1 }$ is positive definite and because

$$
\left[ \begin{array}{c c} I & \mathcal {E} _ {k} y \\ 0 & 1 \end{array} \right] ^ {T} \left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ r ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] \left[ \begin{array}{c c} I & \mathcal {E} _ {k} y \\ 0 & 1 \end{array} \right] = \left[ \begin{array}{c c} T _ {k} & 0 \\ 0 & 1 + r ^ {T} y \end{array} \right].
$$

We have illustrated the kth step of an algorithm proposed by Durbin (1960). It proceeds by solving the Yule-Walker systems

$$
T _ {k} y ^ {(k)} = - r ^ {(k)} = - \left[ r _ {1}, \dots , r _ {k} \right] ^ {T}
$$

for k = 1:n as follows:

$$
y ^ {(1)} = - r _ {1}
$$

for $k = 1 { : } n - 1$

$$
\beta_ {k} = 1 + [ r ^ {(k)} ] ^ {T} y ^ {(k)}
$$

$$
\alpha_ {k} = - (r _ {k + 1} + r ^ {(k) ^ {T}} \mathcal {E} _ {k} y ^ {(k)}) / \beta_ {k} \tag {4.7.1}
$$

$$
z ^ {(k)} = y ^ {(k)} + \alpha_ {k} \mathcal {E} _ {k} y ^ {(k)}
$$

$$
y ^ {(k + 1)} = \left[ \begin{array}{c} z ^ {(k)} \\ \alpha_ {k} \end{array} \right]
$$

end

As it stands, this algorithm would require $3 n ^ { 2 }$ flops to generate $y = y ^ { ( n ) }$ . It is possible, however, to reduce the amount of work even further by exploiting some of the above expressions:

$$
\begin{array}{l} \beta_ {k} = 1 + \left[ r ^ {(k)} \right] ^ {T} y ^ {(k)} \\ = 1 + \left[ \begin{array}{c} r ^ {(k - 1)} \\ r _ {k} \end{array} \right] ^ {T} \left[ \begin{array}{c} y ^ {(k - 1)} + \alpha_ {k - 1} \mathcal {E} _ {k - 1} y ^ {(k - 1)} \\ \alpha_ {k - 1} \end{array} \right] \\ = \left(1 + \left[ r ^ {(k - 1)} \right] ^ {T} y ^ {(k - 1)}\right) + \alpha_ {k - 1} \left(\left[ r ^ {(k - 1)} \right] ^ {T} \mathcal {E} _ {k - 1} y ^ {(k - 1)} + r _ {k}\right) \\ = \beta_ {k - 1} + \alpha_ {k - 1} (- \beta_ {k - 1} \alpha_ {k - 1}) \\ = (1 - \alpha_ {k - 1} ^ {2}) \beta_ {k - 1}. \\ \end{array}
$$

Using this recursion we obtain the following algorithm:

Algorithm 4.7.1 (Durbin) Given real numbers $r _ { 0 } , r _ { 1 } , \ldots , r _ { n }$ with $r _ { 0 } = 1$ such that $T = ( r _ { | i - j | } ) \in \mathbb { R } ^ { n \times n }$ is positive definite, the following algorithm computes $\boldsymbol { y } \in \mathbb { R } ^ { n }$ such that $T \dot { y } = - [ r _ { 1 } , \ldots , r _ { n } ] ^ { T }$ .

$$
y (1) = - r (1); \beta = 1; \alpha = - r (1)
$$

for $k = 1 { : } n - 1$

$$
\beta = (1 - \alpha^ {2}) \beta
$$

$$
\alpha = - \left(r (k + 1) + r (k: - 1: 1) ^ {T} y (1: k)\right) / \beta
$$

$$
z (1: k) = y (1: k) + \alpha y (k: - 1: 1)
$$

$$
y (1: k + 1) = \left[ \begin{array}{c} z (1: k) \\ \alpha \end{array} \right]
$$

end

This algorithm requires $2 n ^ { 2 }$ flops. We have included an auxiliary vector z for clarity, but it can be avoided.

# 4.7.4 The General Right-Hand-Side Problem

With a little extra work, it is possible to solve a symmetric positive definite Toeplitz system that has an arbitrary right-hand side. Suppose that we have solved the system

$$
T _ {k} x = b = \left[ b _ {1}, \dots , b _ {k} \right] ^ {T} \tag {4.7.2}
$$

for some k satisfying $1 \leq k < n$ and that we now wish to solve

$$
\left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ r ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] \left[ \begin{array}{l} v \\ \mu \end{array} \right] = \left[ \begin{array}{c} b \\ b _ {k + 1} \end{array} \right]. \tag {4.7.3}
$$

Here, $\boldsymbol { r } = [ r _ { 1 } , \ldots , r _ { k } ] ^ { T }$ as above. Assume also that the solution to the order-k Yule-Walker system $T _ { k } y = - r$ is also available. From $T _ { k } v + \mu \mathcal { E } _ { k } r = b$ it follows that

$$
v = T _ {k} ^ {- 1} (b - \mu \mathcal {E} _ {k} r) = x - \mu T _ {k} ^ {- 1} \mathcal {E} _ {k} r = x + \mu \mathcal {E} _ {k} y
$$

and so

$$
\begin{array}{l} \mu = b _ {k + 1} - r ^ {T} \mathcal {E} _ {k} v \\ = b _ {k + 1} - r ^ {T} \mathcal {E} _ {k} x - \mu r ^ {T} y \\ = \left(b _ {k + 1} - r ^ {T} \mathcal {E} _ {k} x\right) / \left(1 + r ^ {T} y\right). \\ \end{array}
$$

Consequently, we can effect the transition from (4.7.2) to (4.7.3) in $O ( k )$ flops.

Overall, we can efficiently solve the system $T _ { n } x = b$ by solving the systems

$$
T _ {k} x ^ {(k)} = b ^ {(k)} = [ b _ {1}, \dots , b _ {k} ] ^ {T}
$$

and

$$
T _ {k} y ^ {(k)} = - r ^ {(k)} = - [ r _ {1}, \dots , r _ {k} ] ^ {T}
$$

“in parallel” for $k = 1 { : } n$ . This is the gist of the Levinson algorithm.

Algorithm 4.7.2 (Levinson) Given $b \in \mathbb { R } ^ { n }$ and real numbers $1 = r _ { 0 } , r _ { 1 } , . . . , r _ { n }$ such that $T = ( r _ { | i - j | } ) \in \mathbb { R } ^ { n \times n }$ is positive definite, the following algorithm computes $\boldsymbol { x } \in \mathbb { R } ^ { n }$ such that $T x = b$ .

$$
\begin{array}{l} y (1) = - r (1); x (1) = b (1); \beta = 1; \alpha = - r (1) \\ \beta = (1 - \alpha^ {2}) \beta \\ \mu = \left(b (k + 1) - r (1: k) ^ {T} x (k: - 1: 1)\right) / \beta \\ v (1: k) = x (1: k) + \mu \cdot y (k: - 1: 1) \\ x (1: k + 1) = \left[ \begin{array}{c} v (1: k) \\ \mu \end{array} \right] \\ z (1: k) = y (1: k) + \alpha \cdot y (k: - 1: 1) \\ y (1: k + 1) = \left[ \begin{array}{c} z (1: k) \\ \alpha \end{array} \right] \\ \end{array}
$$

$$
\begin{array}{l} \beta = (1 - \alpha^ {2}) \beta \\ \mu = \left(b (k + 1) - r (1: k) ^ {T} x (k: - 1: 1)\right) / \beta \\ v (1: k) = x (1: k) + \mu \cdot y (k: - 1: 1) \\ x (1: k + 1) = \left[ \begin{array}{c} v (1: k) \\ \mu \end{array} \right] \\ \end{array}
$$

$$
\begin{array}{l} \alpha = - \left(r (k + 1) + r (1: k) ^ {T} y (k: - 1: 1)\right) / \beta \\ z (1: k) = y (1: k) + \alpha \cdot y (k: - 1: 1) \\ y (1: k + 1) = \left[ \begin{array}{c} z (1: k) \\ \alpha \end{array} \right] \\ \end{array}
$$

end

This algorithm requires $4 n ^ { 2 }$ flops. The vectors z and v are for clarity and can be avoided in a detailed implementation.

# 4.7.5 Computing the Inverse

One of the most surprising properties of a symmetric positive definite Toeplitz matrix $T _ { n }$ is that its complete inverse can be calculated in $O ( n ^ { 2 } )$ flops. To derive the algorithm for doing this, partition $T _ { n } ^ { - 1 }$ as follows:

$$
T _ {n} ^ {- 1} = \left[ \begin{array}{c c} A & E r \\ r ^ {T} E & 1 \end{array} \right] ^ {- 1} = \left[ \begin{array}{c c} B & v \\ v ^ {T} & \gamma \end{array} \right] \tag {4.7.4}
$$

where $A = T _ { n - 1 } , E = \mathcal { E } _ { n - 1 }$ , and $r = [ r _ { 1 } , \ldots , r _ { n - 1 } ] ^ { T }$ . From the equation

$$
\left[ \begin{array}{c c} A & E r \\ r ^ {T} E & 1 \end{array} \right] \left[ \begin{array}{c} v \\ \gamma \end{array} \right] = \left[ \begin{array}{c} 0 \\ 1 \end{array} \right]
$$

it follows that $A v = - \gamma E r = - \gamma E ( r _ { 1 } , \ldots , r _ { n - 1 } ) ^ { T }$ and $\gamma = 1 - r ^ { T } E v$ . If y solves the order-(n−1) Yule-Walker system $A y = - r$ , then these expressions imply that

$$
\gamma = 1 / (1 + r ^ {T} y),
$$

$$
v = \gamma E y.
$$

Thus, the last row and column of $T _ { n } ^ { - 1 }$ are readily obtained.

It remains for us to develop working formulae for the entries of the submatrix B in (4.7.4). Since $A B + \mathcal { E } r v ^ { T } = I _ { n - 1 }$ , it follows that

$$
B = A ^ {- 1} - (A ^ {- 1} E r) v ^ {T} = A ^ {- 1} + \frac {v v ^ {T}}{\gamma}.
$$

Now since $A = T _ { n - 1 }$ is nonsingular and Toeplitz, its inverse is persymmetric. Thus,

$$
\begin{array}{l} b _ {i j} = (A ^ {- 1}) _ {i j} + \frac {v _ {i} v _ {j}}{\gamma} \\ = (A ^ {- 1}) _ {n - j, n - i} + \frac {v _ {i} v _ {j}}{\gamma} \tag {4.7.5} \\ = b _ {n - j, n - i} - \frac {v _ {n - j} v _ {n - i}}{\gamma} + \frac {v _ {i} v _ {j}}{\gamma} \\ = b _ {n - j, n - i} + \frac {1}{\gamma} \left(v _ {i} v _ {j} - v _ {n - j} v _ {n - i}\right). \\ \end{array}
$$

This indicates that although B is not persymmetric, we can readily compute an element $b _ { i j }$ from its reflection across the northeast-southwest axis. Coupling this with the fact that $A ^ { - 1 }$ is persymmetric enables us to determine B from its “edges” to its “interior.”

Because the order of operations is rather cumbersome to describe, we preview the formal specification of the algorithm pictorially. To this end, assume that we know the last column and row of $T _ { n } ^ { - 1 }$ :

$$
T _ {n} ^ {- 1} = \left[ \begin{array}{l l l l l l} u & u & u & u & u & k \\ u & u & u & u & u & k \\ u & u & u & u & u & k \\ u & u & u & u & u & k \\ u & u & u & u & u & k \\ k & k & k & k & k & k \end{array} \right].
$$

Here $^ { 6 } u ^ { \dag }$ and $^ { 6 } k ^ { 7 }$ denote the unknown and the known entries, respectively, and $n =$ 6. Alternately exploiting the persymmetry of $T _ { n } ^ { - 1 }$ and the recursion (4.7.5), we can compute B, the leading $( n - 1 ) – \mathrm { b y } – ( n - 1 )$ block of $T _ { n } ^ { - 1 }$ , as follows:

$$
\stackrel {\mathrm{persym}} {\longrightarrow} \left[ \begin{array}{l l l l l l} k & k & k & k & k & k \\ k & u & u & u & u & k \\ k & u & u & u & u & k \\ k & u & u & u & u & k \\ k & u & u & u & u & k \\ k & k & k & k & k & k \end{array} \right] \stackrel {(4. 7. 5)} {\longrightarrow} \left[ \begin{array}{l l l l l l} k & k & k & k & k & k \\ k & u & u & u & k & k \\ k & u & u & u & k & k \\ k & u & u & u & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \end{array} \right] \stackrel {\mathrm{persym}} {\longrightarrow} \left[ \begin{array}{l l l l l l} k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & u & u & k & k \\ k & k & u & u & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \end{array} \right]
$$

$$
\stackrel {(4. 7. 5)} {\longrightarrow} \left[ \begin{array}{c c c c c c} k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & u & k & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \end{array} \right] \stackrel {\text {persym}} {\longrightarrow} \left[ \begin{array}{c c c c c c} k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \end{array} \right].
$$

Of course, when computing a matrix that is both symmetric and persymmetric, such as $T _ { n } ^ { - 1 }$ , it is only necessary to compute the “upper wedge” of the matrix—e.g.,

$$
\begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ & \times & \times & \times & \times \\ & & \times & \times \end{array} \qquad (n = 6).
$$

With this last observation, we are ready to present the overall algorithm.

Algorithm 4.7.3 (Trench) Given real numbers $1 = r _ { 0 } , r _ { 1 } , . . . , r _ { n }$ such that $T =$ $( r _ { | i - j | } ) \in \mathbb { R } ^ { n \times n }$ is positive definite, the following algorithm computes $B = T _ { n } ^ { - 1 }$ . Only those $b _ { i j }$ for which $i \leq j$ and $i + j \le n + 1$ are computed.

Use Algorithm 4.7.1 to solve $T _ { n - 1 } y = - ( r _ { 1 } , \ldots , r _ { n - 1 } ) ^ { T }$

$$
\gamma = 1 / (1 + r (1: n - 1) ^ {T} y (1: n - 1))
$$

$$
v (1: n - 1) = \gamma y (n - 1: - 1: 1)
$$

$$
B (1, 1) = \gamma
$$

$$
B (1, 2: n) = v (n - 1: - 1: 1) ^ {T}
$$

for i = 2 : floor((n − 1)/2) + 1

for j = i:n − i + 1

$$
B (i, j) = B (i - 1, j - 1) + (v (n + 1 - j) v (n + 1 - i) - v (i - 1) v (j - 1)) / \gamma
$$

end

end

This algorithm requires $1 3 n ^ { 2 } / 4$ flops.

# 4.7.6 Stability Issues

Error analyses for the above algorithms have been performed by Cybenko (1978), and we briefly report on some of his findings.

The key quantities turn out to be the $\alpha _ { k }$ in (4.7.1). In exact arithmetic these scalars satisfy

$$
\left| \alpha_ {k} \right| <   1
$$

and can be used to bound $\| T ^ { - 1 } \| _ { 1 }$ :

$$
\max \left\{\frac {1}{\prod_ {j = 1} ^ {n - 1} (1 - \alpha_ {j} ^ {2})}, \frac {1}{\prod_ {j = 1} ^ {n - 1} (1 - \alpha_ {j})} \right\} \leq \| T _ {n} ^ {- 1} \| \leq \prod_ {j = 1} ^ {n - 1} \frac {1 + | \alpha_ {j} |}{1 - | \alpha_ {j} |}. \tag {4.7.6}
$$

Moreover, the solution to the Yule-Walker system $T _ { n } y = - r ( 1 { : } n )$ satisfies

$$
\| y \| _ {1} = \left(\prod_ {k = 1} ^ {n} (1 + \alpha_ {k})\right) - 1 \tag {4.7.7}
$$

provided all the $\alpha _ { k }$ are nonnegative.

Now if ˆx is the computed Durbin solution to the Yule-Walker equations, then the vector $\boldsymbol { r } _ { D } = T _ { n } \boldsymbol { \hat { x } } + \boldsymbol { r }$ can be bounded as follows

$$
\parallel r _ {D} \parallel \approx \mathbf {u} \prod_ {k = 1} ^ {n} (1 + | \hat {\alpha} _ {k} |),
$$

where $\hat { \alpha } _ { k }$ is the computed version of $\alpha _ { k }$ . By way of comparison, since each $| r _ { i } |$ is bounded by unity, it follows that $\parallel r _ { c } \parallel \approx \mathbf { u } \parallel y \parallel _ { 1 }$ where $r _ { C }$ is the residual associated with the computed solution obtained via the Cholesky factorization. Note that the two residuals are of comparable magnitude provided (4.7.7) holds. Experimental evidence suggests that this is the case even if some of the $\alpha _ { k }$ are negative. Similar comments apply to the numerical behavior of the Levinson algorithm.

For the Trench method, the computed inverse $\hat { B }$ of $T _ { n } ^ { - 1 }$ can be shown to satisfy

$$
\frac {\| T _ {n} ^ {- 1} - \hat {B} \| _ {1}}{\| T _ {n} ^ {- 1} \| _ {1}} \approx \mathbf {u} \prod_ {k = 1} ^ {n} \frac {1 + | \hat {\alpha} _ {k} |}{1 - | \hat {\alpha} _ {k} |}.
$$

In light of (4.7.7) we see that the right-hand side is an approximate upper bound for $\mathbf { u } \parallel T _ { n } ^ { - 1 } \parallel$ which is approximately the size of the relative error when $T _ { n } ^ { - 1 }$ is calculated using the Cholesky factorization.

# 4.7.7 A Toeplitz Eigenvalue Problem

Our discussion of the symmetric eigenvalue problem begins in Chapter 8. However, we are able to describe a solution procedure for an important Toeplitz eigenvalue problem that does not require the heavy machinery from that later chapter. Suppose

$$
T = \left[ \begin{array}{l l} 1 & r ^ {T} \\ r & B \end{array} \right]
$$

is symmetric, positive definite, and Toeplitz with $r \in \mathbb { R } ^ { n - 1 }$ . Cybenko and Van Loan (1986) show how to pair the Durbin algorithm with Newton’s method to compute $\lambda _ { \operatorname* { m i n } } ( T )$ assuming that

$$
\lambda_ {\min} (T) <   \lambda_ {\min} (B). \tag {4.7.8}
$$

This assumption is typically the case in practice. If

$$
\left[ \begin{array}{c c} 1 & r ^ {T} \\ r & B \end{array} \right] \left[ \begin{array}{c} \alpha \\ y \end{array} \right] = \lambda_ {\min} \left[ \begin{array}{c} \alpha \\ y \end{array} \right],
$$

then $y = - \alpha ( B - \lambda _ { \operatorname* { m i n } } I ) ^ { - 1 } r , \alpha \neq 0$ , and

$$
\alpha + r ^ {T} \left[ - \alpha (B - \lambda_ {\mathrm{min}} I) ^ {- 1} r \right] = \lambda_ {\mathrm{min}} \alpha .
$$

Thus, $\lambda _ { \mathrm { m i n } }$ is a zero of the rational function

$$
f (\lambda) = 1 - \lambda - r ^ {T} (B - \lambda I) ^ {- 1} r.
$$

Note that if $\lambda < \lambda _ { \operatorname* { m i n } } ( B )$ , then

$$
f ^ {\prime} (\lambda) = - 1 - \left\| (B - \lambda I) ^ {- 1} r \right\| _ {2} ^ {2} \leq - 1,
$$

$$
f ^ {\prime \prime} (\lambda) = - 2 r ^ {T} (B - \lambda I) ^ {- 3} r \leq 0.
$$

Using these facts it can be shown that if

$$
\lambda_ {\min} (T) \leq \lambda^ {(0)} <   \lambda_ {\min} (B), \tag {4.7.9}
$$

then the Newton iteration

$$
\lambda^ {(k + 1)} = \lambda^ {(k)} - \frac {f (\lambda^ {(k)})}{f ^ {\prime} (\lambda^ {(k)})} \tag {4.7.10}
$$

converges to $\lambda _ { \operatorname* { m i n } } ( T )$ monotonically from the right. The iteration has the form

$$
\lambda^ {(k + 1)} = \lambda^ {(k)} + \frac {1 + r ^ {T} w - \lambda^ {(k)}}{1 + w ^ {T} w},
$$

where w solves the “shifted” Yule-Walker system

$$
(B - \lambda^ {(k)} I) w = - r.
$$

Since $\lambda ^ { ( k ) } < \lambda _ { \operatorname* { m i n } } ( B )$ , this system is positive definite and the Durbin algorithm (Algorithm 4.7.1) can be applied to the normalized Toeplitz matrix $( B - \lambda ^ { ( k ) } I ) / ( 1 - \dot { \lambda } ^ { ( k ) } )$ .

The Durbin algorithm can also be used to determine a starting value $\lambda ^ { ( 0 ) }$ that satisfies (4.7.9). If that algorithm is applied to

$$
T _ {\lambda} = (T - \lambda I) / (1 - \lambda)
$$

then it runs to completion if $T _ { \lambda }$ is positive definite. In this case, the $\beta _ { k }$ defined in (4.7.1) are all positive. On the other hand, if $k \leq n - 1 , \beta _ { k } \leq 0$ and $\beta _ { 1 } , \ldots , \beta _ { k - 1 }$ are all positive, then it follows that $T _ { \lambda } ( 1 { : } k , 1 { : } k )$ is positive definite but that $T _ { \lambda } ( 1 { : } k + 1 , k + 1 )$ i s not. Let $m ( \lambda )$ be the index of the first nonpositive $\beta$ and observe that if $m ( \lambda ^ { ( 0 ) } ) = n { - } 1$ , then $B - \lambda ^ { ( 0 ) } I$ is positive definite and $T - \lambda ^ { ( 0 ) } I$ is not, thereby establishing (4.7.9). A bisection scheme can be formulated to compute $\lambda ^ { ( 0 ) }$ with this property:

$$
L = 0
$$

$$
R = 1 - \left| r _ {1} \right|
$$

$$
\mu = (L + R) / 2
$$

while $m ( \mu ) \neq n - 1$

$$
\text { if } m (\mu) <   n - 1
$$

$$
R = \mu
$$

else (4.7.11)

$$
L = \mu
$$

end

$$
\mu = (L + R) / 2
$$

end

$$
\lambda^ {(0)} = \mu
$$

At all times during the iteration we have $m ( L ) \leq n - 1 \leq m ( R )$ . The initial value for R follows from the inequality

$$
0 <   \lambda_ {\min} (T) <   \lambda_ {\min} (B) \leq \lambda_ {\min} \left(\left[ \begin{array}{c c} 1 & r _ {1} \\ r _ {1} & 1 \end{array} \right]\right) = 1 - | r _ {1} |.
$$

Note that the iterations in (4.7.10) and (4.7.11) involve at most $O ( n ^ { 2 } )$ flops per pass. A heuristic argument that O(log n) iterations are required is given by Cybenko and Van Loan (1986).

# 4.7.8 Unsymmetric Toeplitz System Solving

We close with some remarks about unsymmetric Toeplitz system-solving. Suppose we are given scalars $r _ { 1 } , \ldots , r _ { n - 1 } , p _ { 1 } , \ldots , p _ { n - 1 }$ , and $b _ { 1 } , \ldots , b _ { n }$ and that we want to solve a linear system $T x = b$ of the form

$$
\left[ \begin{array}{c c c c c} 1 & r _ {1} & r _ {2} & r _ {3} & r _ {4} \\ p _ {1} & 1 & r _ {1} & r _ {2} & r _ {3} \\ p _ {2} & p _ {1} & 1 & r _ {1} & r _ {2} \\ p _ {3} & p _ {2} & p _ {1} & 1 & r _ {1} \\ p _ {4} & p _ {3} & p _ {2} & p _ {1} & 1 \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \\ x _ {3} \\ x _ {4} \\ x _ {5} \end{array} \right] = \left[ \begin{array}{c} b _ {1} \\ b _ {2} \\ b _ {3} \\ b _ {4} \\ b _ {5} \end{array} \right] \qquad (n = 5).
$$

Assume that $T _ { k } = T ( 1 { : } k , 1 { : } k )$ is nonsingular for $k = 1 { : } n$ . It can shown that if we have the solutions to the k-by-k systems

$$
T _ {k} ^ {T} y = - r = - \left[ r _ {1} r _ {2} \dots r _ {k} \right] ^ {T},
$$

$$
T _ {k} w = - p = - \left[ p _ {1} p _ {2} \dots p _ {k} \right] ^ {T}, \tag {4.7.12}
$$

$$
T _ {k} x = b = \left[ b _ {1} b _ {2} \dots b _ {k} \right] ^ {T},
$$

then we can obtain solutions to

$$
\begin{array}{l} \left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ p ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] ^ {T} \left[ \begin{array}{c} z \\ \alpha \end{array} \right] = - \left[ \begin{array}{c} r \\ r _ {k + 1} \end{array} \right], \\ \left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ p ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] \left[ \begin{array}{l} u \\ \nu \end{array} \right] = - \left[ \begin{array}{c} p \\ p _ {k + 1} \end{array} \right], \tag {4.7.13} \\ \left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ p ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] \left[ \begin{array}{c} v \\ \mu \end{array} \right] = \left[ \begin{array}{c} b \\ b _ {k + 1} \end{array} \right] \\ \end{array}
$$

in $O ( k )$ flops. The update formula derivations are very similar to the Levinson algorithm derivations in §4.7.3. Thus, if the process is repeated for $k = 1 { : } n - 1$ , then we emerge with the solution to $T x = T _ { n } x = b$ . Care must be exercised if a $T _ { k }$ matrix is singular or ill-conditioned. One strategy involves a lookahead idea. In this framework, one might transition from the $T _ { k }$ problem directly to the $T _ { k + 2 }$ problem if it is deemed that the $T _ { k + 1 }$ problem is dangerously ill-conditioned. See Chan and Hansen (1992). An alternative approach based on displacement rank is given in §12.1.

# Problems

P4.7.1 For any $v \in \mathbb { R } ^ { n }$ define the vectors $v _ { + } = ( v + \mathcal { E } _ { n } v ) / 2$ and $v _ { - } = ( v - \mathcal { E } _ { n } v ) / 2$ . Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and persymmetric. Show that if $A x = b$ then $A x _ { + } = b _ { + }$ and $A x _ { - } = b _ { - }$ .

P4.7.2 Let $U \in \mathbb { R } ^ { n \times n }$ be the unit upper triangular matrix with the property that $U ( 1 { : } k - 1 , k ) =$ $\mathcal { E } _ { k - 1 } y ^ { ( k - 1 ) }$ where $y ^ { ( k ) }$ is defined by (4.7.1). Show that $U ^ { T } T _ { n } U = \operatorname { d i a g } ( 1 , \beta _ { 1 } , . . . , \beta _ { n - 1 } )$ .

P4.7.3 Suppose that $z \in \mathbb { R } ^ { n }$ and that $S \in \mathbb { R } ^ { n \times n }$ is orthogonal. Show that if $X = \left[ z , \ S z , \ . . . , S ^ { n - 1 } z \right]$ , then $X ^ { T } X$ is Toeplitz.

P4.7.4 Consider the $\mathrm { L D L ^ { T } }$ factorization of an n-by-n symmetric, tridiagonal, positive definite Toeplitz matrix. Show that $d _ { n }$ and $\ell _ { n , n - 1 }$ converge as $n \to \infty$ .

P4.7.5 Show that the product of two lower triangular Toeplitz matrices is Toeplitz.

P4.7.6 Give an algorithm for determining $\mu \in \mathbb { R }$ such that $T _ { n } + \mu \left( e _ { n } e _ { 1 } ^ { T } + e _ { 1 } e _ { n } ^ { T } \right)$ is singular. Assume $T _ { n } = ( r _ { | i - j | } )$ is positive definite, with $r _ { 0 } = 1$ .

P4.7.7 Suppose $T \in \mathbb { R } ^ { n \times n }$ is symmetric, positive definite, and Toeplitz with unit diagonal. What is the smallest perturbation of the the ith diagonal that makes T semidefinite?

P4.7.8 Rewrite Algorithm 4.7.2 so that it does not require the vectors z and v.

P4.7.9 Give an algorithm for computing $\kappa _ { \infty } ( T _ { k } )$ for $k = 1 { : } n$ .

P4.7.10 A p-by-p block matrix $A ~ = ~ ( A _ { i j } )$ with m-by-m blocks is block Toeplitz if there exist $A _ { - p + 1 } , \hdots , A _ { - 1 } , A _ { 0 } , A _ { 1 } , \hdots , A _ { p - 1 } \in \mathbb { R } ^ { m \times m }$ so that $\begin{array} { r } { A _ { i j } = A _ { i - j } , \mathrm { e . g . } } \end{array}$ ,

$$
A = \left[ \begin{array}{c c c c} A _ {0} & A _ {1} & A _ {2} & A _ {3} \\ A _ {- 1} & A _ {0} & A _ {1} & A _ {2} \\ A _ {- 2} & A _ {- 1} & A _ {0} & A _ {1} \\ A _ {- 3} & A _ {- 2} & A _ {- 1} & A _ {0} \end{array} \right].
$$

(a) Show that there is a permutation Π such that

$$
\Pi^ {T} A \Pi =: \left[ \begin{array}{c c c c} T _ {1 1} & T _ {1 2} & \dots & T _ {1 m} \\ T _ {2 1} & T _ {2 2} & & \vdots \\ \vdots & & \ddots & \vdots \\ T _ {m 1} & \dots & & T _ {m m} \end{array} \right]
$$

where each $T _ { i j }$ is p-by-p and Toeplitz. Each $T _ { i j }$ should be “made up” of $( i , j )$ entries selected from the $A _ { k }$ matrices. (b) What can you say about the $T _ { i j } \mathrm { ~ i f ~ } A _ { k } = A _ { - k } , \bar { k } = 1 : p - 1 ?$

P4.7.11 Show how to compute the solutions to the systems in (4.7.13) given that the solutions to the systems in (4.7.12) are available. Assume that all the matrices involved are nonsingular. Proceed to develop a fast unsymmetric Toeplitz solver for $T x = b$ assuming that $T \ ' _ { \mathrm { s } }$ leading principal submatrices are all nonsingular.

P4.7.12 Consider the order-k Yule-Walker system $T _ { k } y ^ { ( k ) } = - r ^ { ( k ) }$ that arises in (4.7.1). Show that if $y ^ { ( k ) } = [ y _ { k 1 } , \dots , y _ { k k } ] ^ { T }$ for k = 1:n − 1 and

$$
L = \left[ \begin{array}{c c c c c c} 1 & 0 & 0 & 0 & \dots & 0 \\ y _ {1 1} & 1 & 0 & 0 & \dots & 0 \\ y _ {2 2} & y _ {2 1} & 1 & 0 & \dots & 0 \\ \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\ y _ {n - 1, n - 1} & y _ {n - 1, n - 2} & y _ {n - 1, n - 3} & \dots & y _ {n - 1, 1} & 1 \end{array} \right],
$$

then $L ^ { T } T _ { n } L \ = \ \mathrm { d i a g } ( 1 , \beta _ { 1 } , . . . , \beta _ { n - 1 } )$ where $\beta _ { k } = 1 + r ^ { ( k ) ^ { T } } y ^ { ( k ) }$ . Thus, the Durbin algorithm can be thought of as a fast method for computing and $\mathrm { L D L } ^ { T }$ factorization of $T _ { n } ^ { - 1 }$ .

P4.7.13 Show how the Trench algorithm can be used to obtain an initial bracketing interval for the bisection scheme (4.7.11).

# Notes and References for §4.7

The original references for the three algorithms described in this section are as follows:

J. Durbin (1960). “The Fitting of Time Series Models,”Rev. Inst. Int. Stat. 28, 233–243.   
N. Levinson (1947). “The Weiner RMS Error Criterion in Filter Design and Prediction,” J. Math. Phys. 25, 261–278.   
W.F. Trench (1964). “An Algorithm for the Inversion of Finite Toeplitz Matrices,” J. SIAM 12, 515–522.   
As is true with the “fast algorithms” area in general, unstable Toeplitz techniques abound and caution must be exercised, see:   
G. Cybenko (1978). “Error Analysis of Some Signal Processing Algorithms,” PhD Thesis, Princeton University.   
G. Cybenko (1980). “The Numerical Stability of the Levinson-Durbin Algorithm for Toeplitz Systems of Equations,” SIAM J. Sci. Stat. Comput. 1, 303-319.   
J.R. Bunch (1985). “Stability of Methods for Solving Toeplitz Systems of Equations,” SIAM J. Sci. Stat. Comput. 6, 349–364.   
E. Linzer (1992). “On the Stability of Solution Methods for Band Toeplitz Systems,” Lin. Alg. Applic. 170, 1–32.   
J.M. Varah (1994). “Backward Error Estimates for Toeplitz Systems,” SIAM J. Matrix Anal. Applic. 15, 408–417.   
A.W. Bojanczyk, R.P. Brent, F.R. de Hoog, and D.R. Sweet (1995). “On the Stability of the Bareiss and Related Toeplitz Factorization Algorithms,” SIAM J. Matrix Anal. Applic. 16, 40–57.   
M.T. Chu, R.E. Funderlic, and R.J. Plemmons (2003). “Structured Low Rank Approximation,” Lin. Alg. Applic. 366, 157–172.   
A. Bottcher and S. M. Grudsky (2004). “Structured Condition Numbers of Large Toeplitz Matrices are Rarely Better than Usual Condition Numbers,” Num. Lin. Alg. 12, 95–102.   
J.-G. Sun (2005). “A Note on Backwards Errors for Structured Linear Systems,” Numer. Lin. Alg. Applic. 12, 585–603.   
P. Favati, G. Lotti, and O. Menchi (2010). “Stability of the Levinson Algorithm for Toeplitz-Like Systems,” SIAM J. Matrix Anal. Applic. 31, 2531–2552.   
Papers concerned with the lookahead idea include:   
T.F. Chan and P. Hansen (1992). “A Look-Ahead Levinson Algorithm for Indefinite Toeplitz Systems,” SIAM J. Matrix Anal. Applic. 13, 490–506.   
M. Gutknecht and M. Hochbruck (1995). “Lookahead Levinson and Schur Algorithms for Nonhermitian Toeplitz Systems,” Numer. Math. 70, 181–227.

M. Van Barel and A. Bultheel (1997). “A Lookahead Algorithm for the Solution of Block Toeplitz Systems,” Lin. Alg. Applic. 266, 291–335.

Various Toeplitz eigenvalue computations are presented in:

G. Cybenko and C. Van Loan (1986). “Computing the Minimum Eigenvalue of a Symmetric Positive Definite Toeplitz Matrix,” SIAM J. Sci. Stat. Comput. 7, 123–131.

W.F. Trench (1989). “Numerical Solution of the Eigenvalue Problem for Hermitian Toeplitz Matrices,” SIAM J. Matrix Anal. Appl. 10, 135–146.

H. Voss (1999). “Symmetric Schemes for Computing the Minimum Eigenvalue of a Symmetric Toeplitz Matrix,” Lin. Alg. Applic. 287, 359–371.

A. Melman (2004). “Computation of the Smallest Even and Odd Eigenvalues of a Symmetric Positive-Definite Toeplitz Matrix,” SIAM J. Matrix Anal. Applic. 25, 947–963.

# 4.8 Circulant and Discrete Poisson Systems

If $A \in \mathbb { C } ^ { n \times n }$ has a factorization of the form

$$
V ^ {- 1} A V = \Lambda = \mathrm{diag} (\lambda_ {1}, \dots , \lambda_ {n}), \tag {4.8.1}
$$

then the columns of V are eigenvectors and the $\lambda _ { i }$ are the corresponding eigenvalues2. In principle, such a decomposition can be used to solve a nonsingular Au = b problem:

$$
u = A ^ {- 1} b = (V \Lambda V ^ {- 1}) ^ {- 1} b = V (\Lambda^ {- 1} (V ^ {- 1} b)). \tag {4.8.2}
$$

However, if this solution framework is to rival the efficiency of Gaussian elimination or the Cholesky factorization, then V and Λ need to be very special. We say that A has a fast eigenvalue decomposition (4.8.1) if

(1) Matrix-vector products of the form $y = V x$ require O(n log n) flops to evaluate.   
(2) The eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { n }$ require O(n log n) flops to evaluate.   
(3) Matrix-vector products of the form $\tilde { b } = V ^ { - 1 } b$ require O(n log n) flops to evaluate.

If these three properties hold, then it follows from (4.8.2) that $O ( n \log n )$ flops are required to solve $A u = b$ .

Circulant systems and related discrete Poisson systems lend themselves to this strategy and are the main concern of this section. In these applications, the V -matrices are associated with the discrete Fourier transform and various sine and cosine transforms. (Now is the time to review §1.4.1 and §1.4.2 and to recall that we have n log n methods for the DFT, DST, DST2, and DCT.) It turns out that fast methods exist for the inverse of these transforms and that is important because of (3). We will not be concerned with precise flop counts because in the fast transform “business”, some n are friendlier than others from the efficiency point of view. While this issue may be important in practice, it is not something that we have to worry about in our brief, proof-of-concept introduction. Our discussion is modeled after §4.3–§4.5 in Van Loan (FFT) where the reader can find complete derivations and greater algorithmic detail. The interconnection between boundary conditions and fast transforms is a central theme and in that regard we also recommend Strang (1999).

# 4.8.1 The Inverse of the DFT Matrix

Recall from §1.4.1 that the DFT matrix $F _ { n } \in \mathbb { C } ^ { n \times n }$ is defined by

$$
[ F _ {n} ] _ {k j} = \omega_ {n} ^ {(k - 1) (j - 1)}, \qquad \omega_ {n} = \cos \left(\frac {2 \pi}{n}\right) - i \sin \left(\frac {2 \pi}{n}\right).
$$

It is easy to verify that

$$
F _ {n} ^ {H} = \bar {F} _ {n}
$$

and so for all p and q that satisfy $0 \leq p < n$ and $0 \leq q < n$ we have

$$
F _ {n} (:, p + 1) ^ {H} F _ {n} (:, q + 1) = \sum_ {k = 0} ^ {n - 1} \bar {\omega} _ {n} ^ {k p} \omega_ {n} ^ {k q} = \sum_ {k = 0} ^ {n - 1} \omega_ {n} ^ {k (q - p)}.
$$

If $q = p ,$ , then this sum equals n. Otherwise,

$$
\sum_ {k = 0} ^ {n - 1} \omega_ {n} ^ {k (q - p)} = \frac {1 - \omega_ {n} ^ {n (q - p)}}{1 - \omega_ {n} ^ {q - p}} = \frac {1 - 1}{1 - \omega_ {n} ^ {q - p}} = 0.
$$

It follows that

$$
n I _ {n} = F _ {n} ^ {H} F _ {n} = \bar {F} _ {n} F _ {n}.
$$

Thus, the DFT matrix is a scaled unitary matrix and

$$
F _ {n} ^ {- 1} = \frac {1}{n} \bar {F} _ {n}.
$$

A fast Fourier transform procedure for $F _ { n } x$ can be turned into a fast inverse Fourier transform procedure for $F _ { n } ^ { - 1 } x$ . Since

$$
y = F _ {n} ^ {- 1} x = \frac {1}{n} \bar {F} _ {n} x,
$$

simply replace each reference to $\omega _ { n }$ with a reference to $\bar { \omega } _ { n }$ and scale. See Algorithm 1.4.1.

# 4.8.2 Circulant Systems

A circulant matrix is a Toeplitz matrix with “wraparound”, e.g.,

$$
C (z) = \left[ \begin{array}{l l l l l} z _ {0} & z _ {4} & z _ {3} & z _ {2} & z _ {1} \\ z _ {1} & z _ {0} & z _ {4} & z _ {3} & z _ {2} \\ z _ {2} & z _ {1} & z _ {0} & z _ {4} & z _ {3} \\ z _ {3} & z _ {2} & z _ {1} & z _ {0} & z _ {4} \\ z _ {4} & z _ {3} & z _ {2} & z _ {1} & z _ {0} \end{array} \right].
$$

We assume that the vector z is complex. Any circulant $C ( z ) \in \mathbb { C } ^ { n \times n }$ is a linear combination of $I _ { n } , D _ { n } , \ldots , D _ { n } ^ { n - 1 }$ where $\mathcal { D } _ { n }$ is the downshift permutation defined in §1.2.11. For example, if $n = 5$ , then

$$
\mathcal {D} _ {5} = \left[ \begin{array}{l l l l l} 0 & 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \end{array} \right]
$$

and

$$
\mathcal {D} _ {5} ^ {2} = \left[ \begin{array}{l l l l l} 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \end{array} \right], \quad \mathcal {D} _ {5} ^ {3} = \left[ \begin{array}{l l l l l} 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 \end{array} \right], \quad \mathcal {D} _ {5} ^ {4} = \left[ \begin{array}{l l l l l} 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 & 0 \end{array} \right].
$$

Thus, the 5-by-5 circulant matrix displayed above is given by

$$
C (z) = z _ {0} I + z _ {1} \mathcal {D} _ {n} + z _ {2} \mathcal {D} _ {n} ^ {2} + z _ {3} \mathcal {D} _ {n} ^ {3} + z _ {4} \mathcal {D} _ {n} ^ {4}.
$$

Note that $\mathcal { D } _ { 5 } ^ { 5 } = I _ { 5 }$ . More generally,

$$
z = \left[ \begin{array}{c} z _ {0} \\ z _ {1} \\ \vdots \\ z _ {n - 1} \end{array} \right] \quad \Rightarrow \quad C (z) = \sum_ {k = 0} ^ {n - 1} z _ {k} \mathcal {D} _ {n} ^ {k}. \tag {4.8.3}
$$

Note that if $V ^ { - 1 } { \mathcal { D } } _ { n } V = \Lambda$ is diagonal, then

$$
V ^ {- 1} C (z) V = V ^ {- 1} \left(\sum_ {k = 0} ^ {n - 1} z _ {k} \mathcal {D} _ {n} ^ {k}\right) V = \sum_ {k = 0} ^ {n - 1} z _ {k} \left(V ^ {- 1} \mathcal {D} _ {n} V ^ {- 1}\right) ^ {k} = \sum_ {k = 0} ^ {n - 1} z _ {k} \Lambda^ {k} \tag {4.8.4}
$$

is diagonal. It turns out that the DFT matrix diagonalizes the downshift permutation.

Lemma 4.8.1. $I f V = F _ { n }$ , then $V ^ { - 1 } { \mathcal D } _ { n } V = \Lambda = \mathrm { d i a g } ( \lambda _ { 1 } , \dots , \lambda _ { n } )$ where

$$
\lambda_ {j + 1} = \bar {\omega} _ {n} ^ {j} = \cos \left(\frac {2 j \pi}{n}\right) + i \sin \left(\frac {2 j \pi}{n}\right)
$$

for $j = 0 { : } n - 1$ .

Proof. For $j = 0 { : } n - 1$ we have

$$
\mathcal {D} _ {n} F _ {n} (:, j + 1) = \mathcal {D} _ {n} \left[ \begin{array}{c} 1 \\ \omega_ {n} ^ {j} \\ \omega_ {n} ^ {2 j} \\ \vdots \\ \omega_ {n} ^ {(n - 1) j} \end{array} \right] = \left[ \begin{array}{c} \omega_ {n} ^ {(n - 1) j} \\ 1 \\ \omega_ {n} ^ {j} \\ \vdots \\ \omega_ {n} ^ {(n - 2) j} \end{array} \right] = \bar {\omega} _ {n} ^ {j} \left[ \begin{array}{c} 1 \\ \omega_ {n} ^ {j} \\ \omega_ {n} ^ {2 j} \\ \vdots \\ \omega_ {n} ^ {(n - 1) j} \end{array} \right].
$$

This vector is precisely $F _ { n } \Lambda ( : , j + 1 )$ . Thus, $\mathcal { D } _ { n } V = V \Lambda$ , i.e., $V ^ { - 1 } { \mathcal { D } } _ { n } V = \Lambda$ .

It follows from (4.8.4) that any circulant $C ( z )$ is diagonalized by $F _ { n }$ and the eigenvalues of $C ( z )$ can be computed fast.

Theorem 4.8.2. Suppose $z \in \mathbb { C } ^ { n }$ and C(z) are defined by (4.8.3). $I f V = F _ { n }$ and $\lambda \ : = \ : \bar { F } _ { n } z$ , then $V ^ { - 1 } C ( z ) V = \mathrm { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ .

Proof. Define

$$
f = \left[ \begin{array}{c} 1 \\ \bar {\omega} _ {n} \\ \vdots \\ \bar {\omega} _ {n} ^ {n - 1} \end{array} \right]
$$

and note that the columns of ${ \bar { F } } _ { n }$ are componentwise powers of this vector. In particular, $\bar { F } _ { n } ( : , k + 1 ) = f . \hat { \cdot } k$ where $[ f . \hat { \cdot } k ] _ { j } = f _ { j } ^ { k }$ . Since $\Lambda = \operatorname { d i a g } ( f )$ , it follows from Lemma 4.8.1 that

$$
\begin{array}{l} V ^ {- 1} C (z) V = \sum_ {k = 0} ^ {n - 1} z _ {k} \Lambda^ {k} = \sum_ {k = 0} ^ {n - 1} z _ {k} \operatorname{diag} (f) ^ {k} = \sum_ {k = 0} ^ {n - 1} z _ {k} \operatorname{diag} (f. \hat {} k) \\ = \operatorname{diag} \left(\sum_ {k = 0} ^ {n - 1} z _ {k} f. \hat {k}\right) = \operatorname{diag} \left(\bar {F} _ {n} z\right) \\ \end{array}
$$

completing the proof of the theorem

Thus, the eigenvalues of the circulant matrix $C ( z )$ are the components of the vector $\bar { F } _ { n } z$ . Using this result we obtain the following algorithm.

Algorithm 4.8.1 If $z \in \mathbb { C } ^ { n } , \quad y \in \mathbb { C } ^ { n }$ , and $C ( z )$ is nonsingular, then the following algorithm solves the linear system $C ( z ) x = y$ .

Use an FFT to compute $c = { \bar { F } } _ { n } y$ and $d = { \bar { F } } _ { n } z$ .

$$
w = c. / d
$$

Use an FFT to compute $u \ = \ F _ { n } w$

$$
x = u / n
$$

This algorithm requires $O ( n \log n )$ flops.

# 4.8.3 The Discretized Poisson Equation in One Dimension

We now turn our attention to a family of real matrices that have real, fast eigenvalue decompositions. The starting point in the discussion is the differential equation

$$
\frac {d ^ {2} u}{d x ^ {2}} = - f (x) \quad \alpha \leq u (x) \leq \beta , \tag {4.8.5}
$$

together with one of four possible specifications of $u ( x )$ on the boundary.

Dirichlet-Dirichlet (DD): $u ( \alpha ) = u _ { \alpha } , \qquad u ( \beta ) = u _ { \beta } ,$

Dirichlet-Neumann (DN): $u ( \alpha ) = u _ { \alpha } , \qquad u ^ { \prime } ( \beta ) = u _ { \beta } ^ { \prime } ,$

Neumann-Neumann (NN): $u ^ { \prime } ( \alpha ) = u _ { \alpha } ^ { \prime } , \qquad u ^ { \prime } ( \beta ) = u _ { \beta } ^ { \prime } .$

Periodic (P): $u ( \alpha ) = u ( \beta ) .$ .

By replacing the derivatives in (4.8.5) with divided differences, we obtain a system of linear equations. Indeed, if m is a positive integer and

$$
h = \frac {\beta - \alpha}{m},
$$

then for $i = 1 { : } m - 1$ we have

$$
\frac {\frac {u _ {i + 1} - u _ {i}}{h} - \frac {u _ {i} - u _ {i - 1}}{h}}{h} = \frac {u _ {i - 1} - 2 u _ {i} + u _ {i + 1}}{h ^ {2}} = - f _ {i} \tag {4.8.6}
$$

where $f _ { i } = f ( \alpha { + } i h )$ and $u _ { i } \approx u ( \alpha { + } i h )$ . To appreciate this discretization we display the linear equations that result when $m = 5$ for the various possible boundary conditions. The matrices $\mathcal { T } _ { n } ^ { ( { D D } ) } , \mathcal { T } _ { n } ^ { ( { D N } ) } , \mathcal { T } _ { n } ^ { ( { N N } ) }$ , and ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ are formally defined afterwards.

For the Dirichlet-Dirichlet problem, the system is 4-by-4 and tridiagonal:

$$
\mathcal {T} _ {4} ^ {(D D)} \cdot u (1: 4) \equiv \left[ \begin{array}{r r r r} 2 & - 1 & 0 & 0 \\ - 1 & 2 & - 1 & 0 \\ 0 & - 1 & 2 & - 1 \\ 0 & 0 & - 1 & 2 \end{array} \right] \left[ \begin{array}{l} u _ {1} \\ u _ {2} \\ u _ {3} \\ u _ {4} \end{array} \right] = \left[ \begin{array}{l} h ^ {2} f _ {1} + u _ {\alpha} \\ h ^ {2} f _ {2} \\ h ^ {2} f _ {3} \\ h ^ {2} f _ {4} + u _ {\beta} \end{array} \right].
$$

For the Dirichlet-Neumann problem the system is still tridiagonal, but $u _ { 5 }$ joins $u _ { 1 } , \ldots , u _ { 4 }$ as an unknown:

$$
\mathcal {T} _ {5} ^ {(D N)} \cdot u (1: 5) \equiv \left[ \begin{array}{r r r r r} 2 & - 1 & 0 & 0 & 0 \\ - 1 & 2 & - 1 & 0 & 0 \\ 0 & - 1 & 2 & - 1 & 0 \\ 0 & 0 & - 1 & 2 & - 1 \\ 0 & 0 & 0 & - 2 & 2 \end{array} \right] \left[ \begin{array}{l} u _ {1} \\ u _ {2} \\ u _ {3} \\ u _ {4} \\ u _ {5} \end{array} \right] = \left[ \begin{array}{c} h ^ {2} f _ {1} + u _ {\alpha} \\ h ^ {2} f _ {2} \\ h ^ {2} f _ {3} \\ h ^ {2} f _ {4} \\ 2 h u _ {\beta} ^ {\prime} \end{array} \right].
$$

The new equation on the bottom is derived from the approximation $u ^ { \prime } ( \beta ) \approx ( u _ { 5 } - u _ { 4 } ) / h$ . (The scaling of this equation by 2 simplifies some of the derivations below.) For the Neumann-Neumann problem, $u _ { 5 }$ and $u _ { 0 }$ need to be determined:

$$
\mathcal {T} _ {6} ^ {(N N)} \cdot u (0 {:} 5) \equiv \left[ \begin{array}{r r r r r r} 2 & - 2 & 0 & 0 & 0 & 0 \\ - 1 & 2 & - 1 & 0 & 0 & 0 \\ 0 & - 1 & 2 & - 1 & 0 & 0 \\ 0 & 0 & - 1 & 2 & - 1 & 0 \\ 0 & 0 & 0 & - 1 & 2 & - 1 \\ 0 & 0 & 0 & 0 & - 2 & 2 \end{array} \right] \left[ \begin{array}{l} u _ {0} \\ u _ {1} \\ u _ {2} \\ u _ {3} \\ u _ {4} \\ u _ {5} \end{array} \right] = \left[ \begin{array}{l} - 2 h u _ {\alpha} ^ {\prime} \\ h ^ {2} f _ {1} \\ h ^ {2} f _ {2} \\ h ^ {2} f _ {3} \\ h ^ {2} f _ {3} \\ 2 h u _ {\beta} ^ {\prime} \end{array} \right].
$$

Finally, for the periodic problem we have

$$
\mathcal {T} _ {5} ^ {(P)} \cdot u (1: 5) \equiv \left[ \begin{array}{r r r r r} 2 & - 1 & 0 & 0 & - 1 \\ - 1 & 2 & - 1 & 0 & 0 \\ 0 & - 1 & 2 & - 1 & 0 \\ 0 & 0 & - 1 & 2 & - 1 \\ - 1 & 0 & 0 & - 1 & 2 \end{array} \right] \left[ \begin{array}{l} u _ {1} \\ u _ {2} \\ u _ {3} \\ u _ {4} \\ u _ {5} \end{array} \right] = \left[ \begin{array}{l} h ^ {2} f _ {1} \\ h ^ {2} f _ {2} \\ h ^ {2} f _ {3} \\ h ^ {2} f _ {4} \\ h ^ {2} f _ {5} \end{array} \right].
$$

The first and last equations use the conditions $u _ { 0 } = u _ { 5 }$ and $u _ { 1 } = u _ { 6 }$ . These constraints follow from the assumption that u has period $\beta - \alpha$ .

As we show below, the n-by-n matrix

$$
\mathcal {T} _ {n} ^ {(D D)} = \left[ \begin{array}{c c c c} 2 & - 1 & \dots & 0 \\ - 1 & 2 & \ddots & \vdots \\ \vdots & \ddots & \ddots & - 1 \\ 0 & \dots & - 1 & 2 \end{array} \right] \tag {4.8.7}
$$

and its low-rank adjustments

$$
\mathcal {T} _ {n} ^ {(D N)} = \mathcal {T} _ {n} ^ {(D D)} - e _ {n} e _ {n - 1} ^ {T}, \tag {4.8.8}
$$

$$
\mathcal {T} _ {n} ^ {(N N)} = \mathcal {T} _ {n} ^ {(D D)} - e _ {n} e _ {n - 1} ^ {T} - e _ {1} e _ {2} ^ {T}, \tag {4.8.9}
$$

$$
\mathcal {T} _ {n} ^ {(P)} = \mathcal {T} _ {n} ^ {(D D)} - e _ {1} e _ {n} ^ {T} - e _ {n} e _ {1} ^ {T}. \tag {4.8.10}
$$

have fast eigenvalue decompositions. However, the existence of O(n log n) methods for these systems is not very interesting because algorithms based on Gaussian elimination are faster: $O ( n )$ versus O(n log n). Things get much more interesting when we discretize the 2-dimensional analogue of (4.8.5).

# 4.8.4 The Discretized Poisson Equation in Two Dimensions

To launch the 2D discussion, suppose $F ( x , y )$ is defined on the rectangle

$$
R = \{(x, y): \alpha_ {x} \leq x \leq \beta_ {x}, \alpha_ {y} \leq y \leq \beta_ {y} \}
$$

and that we wish to find a function u that satisfies

$$
\frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}} = - F (x, y) \tag {4.8.11}
$$

on R and has its value prescribed on the boundary of R. This is Poisson’s equation with Dirichlet boundary conditions. Our plan is to approximate u at the grid points $( \alpha _ { x } + i h _ { x } , \alpha _ { y } + j h _ { y } )$ where $i = 1 { : } m _ { 1 } - 1 , j = 1 { : } m _ { 2 } - 1$ , and

$$
h _ {x} = \frac {\beta_ {x} - \alpha_ {x}}{m _ {1}} \qquad h _ {y} = \frac {\beta_ {y} - \alpha_ {y}}{m _ {2}}.
$$

Refer to Figure 4.8.1, which displays the case when $m _ { 1 } = 6$ and $m _ { 2 } = 5$ . Notice that there are two kinds of grid points. The function u is known at the $\mathbf { \widetilde { \Gamma } } ^ { ( 6 } \bullet ^ { 9 }$ grid points on the boundary. The function u is to be determined at the “◦” grid points in the interior. The interior grid points have been indexed in a top-to-bottom, left-to-right order. The idea is to have $u _ { k }$ approximate the value of $u ( x , y )$ at grid point k.

As in the one-dimensional problem considered §4.8.3, we use divided differences to obtain a set of linear equations that define the unknowns. An interior grid point P has a north (N ), east (E), south (S), and west (W ) neighbor. Using this “compass point” notation we obtain the following approximation to (4.8.11) at $P 3$

$$
\begin{array}{r l r} \frac {u (E) - u (P)}{h _ {x}} - \frac {u (P) - u (W)}{h _ {x}} & + & \frac {u (N) - u (P)}{h _ {y}} - \frac {u (P) - u (S)}{h _ {y}} \\ \hline h _ {x} & + & \frac {h _ {y}}{h _ {y}} = - F (P) \end{array}
$$

![](images/golub_200_249__541d3ad665f3a4144be5b650e770d934fb9ac0e547bc33cba358850ba3e872af.jpg)

<details>
<summary>text_image</summary>

1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
</details>

Figure 4.8.1. A grid with $m _ { 1 } = 6$ and $m _ { 2 } = 5$ .

The x-partial and y-partial have been replaced by second-order divided differences. Assume for clarity that the horizontal and vertical grid spacings are equal, i.e., $h _ { x } =$ $h _ { y } = h$ . With this assumption, the linear equation at point $P$ has the form

$$
4 u (P) - u (N) - u (E) - u (S) - u (W) = h ^ {2} F (P).
$$

In our example, there are 20 such equations. It should be noted that some of $P ^ { * } { \mathrm { s } }$ neighbors may be on the boundary, in which case the corresponding linear equation involves fewer than 5 unknowns. For example, if P is the third grid point then we see from Figure 4.8.1 that the north neighbor N is on the boundary. It follows that the associated linear equation has the form

$$
4 u (P) - u (E) - u (S) - u (W) = h ^ {2} F (P) + u (N).
$$

Reasoning like this, we conclude that the matrix of coefficients has the following block tridiagonal form

$$
A = \left[ \begin{array}{c c c c} \mathcal {T} _ {5} ^ {(D D)} & 0 & 0 & 0 \\ 0 & \mathcal {T} _ {5} ^ {(D D)} & 0 & 0 \\ 0 & 0 & \mathcal {T} _ {5} ^ {(D D)} & 0 \\ 0 & 0 & 0 & \mathcal {T} _ {5} ^ {(D D)} \end{array} \right] + \left[ \begin{array}{c c c c} 2 I _ {5} & - I _ {5} & 0 & 0 \\ - I _ {5} & 2 I _ {5} & - I _ {5} & 0 \\ 0 & - I _ {5} & 2 I _ {5} & - I _ {5} \\ 0 & 0 & - I _ {5} & 2 I _ {5} \end{array} \right]
$$

i.e.,

$$
A = I _ {4} \otimes \mathcal {T} _ {5} ^ {(D D)} + \mathcal {T} _ {4} ^ {(D D)} \otimes I _ {5}.
$$

Notice that the first matrix is associated with the x-partials while the second matrix is associated with the y-partials. The right-hand side in $A u = b$ is made up of $F _ { - }$ evaluations and specified values of $u ( x , y )$ on the boundary.


---

<!-- golub_250_299 -->

Extrapolating from our example, we conclude that the matrix of coefficients is an $( m _ { 2 } - 1 ) – \mathrm { b y } – ( m _ { 2 } - 1 )$ block tridiagonal matrix with $( m _ { 1 } - 1 ) – \mathrm { b y } – ( m _ { 1 } - 1 )$ blocks:

$$
A = I _ {m _ {2} - 1} \otimes \mathcal {T} _ {m _ {1} - 1} ^ {(D D)} + \mathcal {T} _ {m _ {2} - 1} ^ {(D D)} \otimes I _ {m _ {1} - 1}.
$$

Alternative specifications along the boundary lead to systems with similar structure, $\mathrm { e . g . }$ ,

$$
A u \equiv \left(I _ {n _ {2}} \otimes A _ {1} + A _ {2} \otimes I _ {n _ {1}}\right) u = b. \tag {4.8.12}
$$

For example, if we impose Dirichet-Neumann, Neumann-Neumann, or periodic boundary coequal the, or and right edges of the rectangular domain accordingly. Likewise, if we impose Dir $R ,$ then et-Ne $A _ { 1 }$ willann, $\mathscr { T } _ { m _ { 1 } } ^ { ( D N ) } , \mathscr { T } _ { m _ { 1 } + 1 } ^ { ( N N ) }$ $\mathcal { T } _ { m _ { 1 } } ^ { \left( P \right) }$ Neumann-Neumann, or periodic boundary conditions along the bottom and top edges of R, then A2 will equal T (DNm2 $R ,$ $A _ { 2 }$ $\mathcal { T } _ { m _ { 2 } } ^ { ( D N ) } , \mathcal { T } _ { m _ { 2 } + 1 } ^ { ( N N ) }$ , or $\mathcal { T } _ { m _ { 2 } } ^ { \left( P \right) }$ . If the system (4.8.12) is nonsingular and $A _ { 1 }$ and $A _ { 2 }$ have fast eigenvalue decompositions, then it can be solved with just O(N log N) flops where $N = n _ { 1 } n _ { 2 }$ . To see why this is possible, assume that

$$
V ^ {- 1} A _ {1} V = D _ {1} = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n _ {1}}\right), \tag {4.8.13}
$$

$$
W ^ {- 1} A _ {2} W = D _ {2} = \operatorname{diag} \left(\mu_ {1}, \dots , \mu_ {n _ {2}}\right) \tag {4.8.14}
$$

are fast eigenvalue decompositions. Using facts about the Kronecker product that are set forth in §1.3.6–§1.3.8, we can reformulate (4.8.12) as a matrix equation

$$
A _ {1} U + U A _ {2} ^ {T} = B
$$

where $U = { \mathsf { r e s h a p e } } ( u , n _ { 1 } , n _ { 2 } )$ and $\boldsymbol { B } = \mathsf { r e s h a p e } ( b , n _ { 1 } , n _ { 2 } )$ . Substituting the above eigenvalue decompositions into this equation we obtain

$$
D _ {1} \tilde {U} + \tilde {U} D _ {2} = \tilde {B},
$$

where $\tilde { U } = ( \tilde { u } _ { i j } ) = V ^ { - 1 } U W ^ { - T }$ and $\begin{array} { r } { \tilde { B } = ( \tilde { b } _ { i j } ) = V ^ { - 1 } B W ^ { - T } } \end{array}$ . Note how easy it is to solve this transformed system because $D _ { 1 }$ and $D _ { 2 }$ are diagonal:

$$
\tilde {u} _ {i j} = \frac {\tilde {b} _ {i j}}{\lambda_ {i} + \mu_ {j}} \qquad i = 1: n _ {1}, j = 1: n _ {2}.
$$

For this to be well-defined, no eigenvalue of $A _ { 1 }$ can be the negative of an eigenvalue of $A _ { 2 }$ . In our example, all the $\lambda _ { i }$ and $\mu _ { i }$ are positive. Overall we obtain

Algorithm 4.8.2 (Fast Poisson Solver Framework) Assume that $A _ { 1 } \in \mathbb { R } ^ { n _ { 1 } \times n _ { 1 } }$ and $A _ { 2 } \in \mathbb { R } ^ { n _ { 2 } \times n _ { 2 } }$ have fast eigenvalue decompositions (4.8.13) and (4.8.14) and that the matrix $A = I _ { n _ { 2 } } \otimes A _ { 1 } + A _ { 2 } \otimes I _ { n _ { 1 } }$ is nonsingular. The following algorithm solves the linear system Au = b where b ∈ IRn1n2 . $A u = b$ $b \in \mathbb { R } ^ { n _ { 1 } n _ { 2 } }$

$$
\tilde {B} = (W ^ {- 1} (V ^ {- 1} B) ^ {T}) ^ {T} \text {   where   } B = \operatorname{reshape} (b, n _ {1}, n _ {2})
$$

for $i = 1 { : } n _ { 1 }$

$$
\text { for } j = 1: n _ {2}
$$

$$
\tilde {u} _ {i j} = \tilde {b} _ {i j} / (\lambda_ {i} + \mu_ {j})
$$

end

end

$$
u = \operatorname{reshape} (U, n _ {1} n _ {2}, 1) \text {   where   } U = (W (V \tilde {U}) ^ {T}) ^ {T}
$$

The following table accounts for the work involved:

<table><tr><td>Operation</td><td>How Many?</td><td>Work</td></tr><tr><td> $V^{-1}$  times  $n_1$ -vector</td><td> $n_2$ </td><td> $O(n_2 \cdot n_1 \cdot \log n_1)$ </td></tr><tr><td> $W^{-1}$  times  $n_2$ -vector</td><td> $n_1$ </td><td> $O(n_1 \cdot n_2 \cdot \log n_2)$ </td></tr><tr><td> $V$  times  $n_1$ -vector</td><td> $n_2$ </td><td> $O(n_2 \cdot n_1 \cdot \log n_1)$ </td></tr><tr><td> $W$  times  $n_2$ -vector</td><td> $n_1$ </td><td> $O(n_1 \cdot n_2 \cdot \log n_2)$ </td></tr></table>

Adding up the operation counts, we see that $O ( n _ { 1 } n _ { 2 } \log ( n _ { 1 } n _ { 2 } ) ) = O ( N \log N )$ flops are required where $N = n _ { 1 } n _ { 2 }$ is the size of the matrix A.

Below we show that the matrices $\mathcal { T } _ { n } ^ { ( { D D } ) } , \mathcal { T } _ { n } ^ { ( { D N } ) } , \mathcal { T } _ { n } ^ { ( { N N } ) }$ , and ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ have fast eigenvalue decompositions and this means that Algorithm 4.8.2 can be used to solve discrete Poisson systems. To appreciate the speedup over conventional methods, suppose $A _ { 1 } = \mathcal { T } _ { n _ { 1 } } ^ { \left( D D \right) }$ and $A _ { 2 } = \mathcal { T } _ { n _ { 2 } } ^ { \left( D D \right) }$ . It can be shown that A is symmetric positive definite with bandwidth $n _ { 1 } + 1$ . Solving $A u = b$ using Algorithm 4.3.5 (band Cholesky) would require $O ( n _ { 1 } ^ { 3 } n _ { 2 } ) = O ( N n _ { 1 } ^ { 2 } )$ flops.

# 4.8.5 The Inverse of the DST and DCT Matrices

The eigenvector matrices for $\mathcal { T } _ { n } ^ { ( { D D } ) } , \mathcal { T } _ { n } ^ { ( { D N } ) } , \mathcal { T } _ { n } ^ { ( { N N } ) }$ , and T (P ) n ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ are associated with the fast trigonometric transforms presented in §1.4.2. It is incumbent upon us to show that the inverse of these transforms can also be computed fast. We do this for the discrete sine transform (DST) and the discrete cosine transform (DCT) and leave similar fast inverse verifications to the exercises at the end of the section.

By considering the blocks of the DFT matrix $F _ { 2 m }$ , we can determine the inverses of the transform matrices $\mathrm { D S T } ( m - 1 )$ and $\mathrm { D C T } ( m + 1 )$ . Recall from §1.4.2 that if $C _ { r } \in \mathbb { R } ^ { r \times r }$ and $S _ { r } \in \mathbb { R } ^ { r \times r }$ are defined by

$$
[ C _ {r} ] _ {k j} = \cos \left(\frac {k j \pi}{r + 1}\right), [ S _ {r} ] _ {k j} = \sin \left(\frac {k j \pi}{r + 1}\right)
$$

then

$$
F _ {2 m} = \left[ \begin{array}{c c c c} 1 & e ^ {T} & 1 & e ^ {T} \\ e & C - i S & v & (C + i S) E \\ 1 & v ^ {T} & (- 1) ^ {m} & v ^ {T} E \\ e & E (C + i S) & E v & E (C - i S) E \end{array} \right]
$$

where $C = C _ { m - 1 } , S = S _ { m - 1 } , E = \mathcal { E } _ { m - 1 }$ , and

$$
e ^ {T} = (\underbrace {1 , 1 , \ldots , 1} _ {m - 1}) \qquad v ^ {T} = (\underbrace {- 1 , 1 , \ldots , (- 1) ^ {m - 1}} _ {m - 1}).
$$

By comparing the (2,1), (2,2), (2,3), and (2,4) blocks in the equation $2 m I = \bar { F } _ { 2 m } F _ { 2 m }$ we conclude that

$$
0 = 2 C e + e + v,
$$

$$
\begin{array}{l} 2 m I _ {m - 1} = 2 C ^ {2} + 2 S ^ {2} + e e ^ {T} + v v ^ {T}, \\ 0 = 2 C v + e + (- 1) ^ {m} v, \\ 0 = 2 C ^ {2} - 2 S ^ {2} + e e ^ {T} + v v ^ {T}. \\ \end{array}
$$

It follows that $2 S ^ { 2 } = m I _ { m - 1 }$ and $2 C ^ { 2 } = m I _ { m - 1 } - e e ^ { T } - v v ^ { T }$ . Using these equations it is easy to verify that

$$
S _ {m - 1} ^ {- 1} = \frac {2}{m} S _ {m - 1}
$$

and

$$
\left[ \begin{array}{c c c} 1 / 2 & e ^ {T} & 1 / 2 \\ e / 2 & C _ {m - 1} & v / 2 \\ 1 / 2 & v ^ {T} & (- 1) ^ {m} / 2 \end{array} \right] ^ {- 1} = \frac {2}{m} \left[ \begin{array}{c c c} 1 / 2 & e ^ {T} & 1 / 2 \\ e / 2 & C _ {m - 1} & v / 2 \\ 1 / 2 & v ^ {T} & (- 1) ^ {m} / 2 \end{array} \right].
$$

Thus, it follows from the definitions (1.4.8) and (1.4.10) that

$$
V = \mathrm{DST} (m - 1) \Rightarrow V ^ {- 1} = \frac {2}{m} \mathrm{DST} (m - 1),
$$

$$
V = \mathrm{DCT} (m + 1) \Rightarrow V ^ {- 1} = \frac {2}{m} \mathrm{DCT} (m + 1).
$$

In both cases, the inverse transform is a multiple of the “forward” transform and can be computed fast. See Algorithms 1.4.2 and 1.4.3.

# 4.8.6 Four Fast Eigenvalue Decompositions

The matrices $\mathcal { T } _ { n } ^ { ( { D D } ) } , \mathcal { T } _ { n } ^ { ( { D N } ) } , \mathcal { T } _ { n } ^ { ( { N N } ) }$ , and ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ do special things to vectors of sines and cosines.

Lemma 4.8.3. Define the real n-vectors s(θ) and $c ( \theta )$ by

$$
s (\theta) = \left[ \begin{array}{c} s _ {1} \\ \vdots \\ s _ {n} \end{array} \right], \quad c (\theta) = \left[ \begin{array}{c} c _ {0} \\ \vdots \\ c _ {n - 1} \end{array} \right], \tag {4.8.15}
$$

where $s _ { k } = \sin ( k \theta )$ and $c _ { k } = \cos ( k \theta )$ . $\textit { I f e } _ { k } = \textit { I } _ { n } ( : , k )$ and $\lambda = 4 \sin ^ { 2 } ( \theta / 2 )$ , then

$$
\mathcal {T} _ {n} ^ {(D D)} \cdot s (\theta) = \lambda \cdot s (\theta) + s _ {n + 1} e _ {n}, \tag {4.8.16}
$$

$$
\mathcal {T} _ {n} ^ {(D D)} \cdot c (\theta) = \lambda \cdot c (\theta) + c _ {1} e _ {1} + c _ {n} e _ {n}, \tag {4.8.17}
$$

$$
\mathcal {T} _ {n} ^ {(D N)} \cdot s (\theta) = \lambda \cdot s (\theta) + (s _ {n + 1} - s _ {n - 1}) e _ {n}, \tag {4.8.18}
$$

$$
\mathcal {T} _ {n} ^ {(N N)} \cdot c (\theta) = \lambda \cdot c (\theta) + (c _ {n} - c _ {n - 2}) e _ {n}, \tag {4.8.19}
$$

$$
\mathcal {T} _ {n} ^ {(P)} \cdot s (\theta) = \lambda \cdot s (\theta) - s _ {n} e _ {1} + (s _ {n + 1} - s _ {1}) e _ {n}, \tag {4.8.20}
$$

$$
\mathcal {T} _ {n} ^ {(P)} \cdot c (\theta) = \lambda \cdot c (\theta) + (c _ {1} - c _ {n - 1}) e _ {1} + (c _ {n} - 1) e _ {n}. \tag {4.8.21}
$$

Proof. The proof is mainly an exercise in using the trigonometric identities

$$
s _ {k - 1} = c _ {1} s _ {k} - s _ {1} c _ {k}, \quad c _ {k - 1} = c _ {1} c _ {k} + s _ {1} s _ {k},
$$

$$
s _ {k + 1} = c _ {1} s _ {k} + s _ {1} c _ {k}, \quad c _ {k + 1} = c _ {1} c _ {k} - s _ {1} s _ {k}.
$$

For example, if $y = { \mathcal { T } _ { n } ^ { \left( D D \right) } } s ( \theta )$ , then

$$
y _ {k} = \left\{ \begin{array}{l l} 2 s _ {1} - s _ {2} = 2 s _ {1} (1 - c _ {1}), & \text { if } k = 1, \\ - s _ {k - 1} + 2 s _ {k} - s _ {k + 1} = 2 s _ {k} (1 - c _ {1}), & \text { if } 2 \leq k \leq n - 1, \\ - s _ {n - 1} + 2 s _ {n} = 2 s _ {n} (1 - c _ {1}) + s _ {n + 1}, & \text { if } k = n. \end{array} \right.
$$

Equation (4.8.16) follows since $( 1 - c _ { 1 } ) = 1 - \cos ( \theta ) = 2 \sin ^ { 2 } ( \theta / 2 )$ . The proof of (4.8.17) is similar while the remaining equations follow from Equations (4.8.8)–(4.8.10).

Notice that (4.8.16)-(4.8.21) are eigenvector equations except for the $" e _ { 1 } \ "$ and $" e _ { n } "$ terms. By choosing the right value for θ, we can make these residuals disappear, thereby obtaining recipes for the eigensystems of $\mathcal { T } _ { n } ^ { ( { D D } ) } , \mathcal { T } _ { n } ^ { ( { D N } ) } , \mathcal { T } _ { n } ^ { ( { N N } ) }$ T n N) , and T (P ) n . ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$

The Dirichlet-Dirichlet Matrix

If j is an integer and $\theta = j \pi / ( n + 1 )$ , then $s _ { n + 1 } = \sin ( ( n + 1 ) \theta ) = 0$ . It follows from (4.8.16) that

$$
\mathcal {T} _ {n} ^ {(D D)} s (\theta_ {j}) = 4 \sin^ {2} (\theta_ {j} / 2) s (\theta_ {j}), \qquad \theta_ {j} = \frac {j \pi}{n + 1},
$$

for $j = 1 { : } n$ . Thus, the columns of the matrix $V _ { n } ^ { \left( D D \right) } \in \mathbb { R } ^ { n \times n }$ defined by

$$
[ V _ {n} ^ {(D D)} ] _ {k j} = \sin \left(\frac {k j \pi}{n + 1}\right)
$$

are eigenvectors for T (DD) n ${ \mathcal { T } } _ { n } ^ { ( D D ) }$ and the corresponding eigenvalues are given by

$$
\lambda_ {j} = 4 \sin^ {2} \left(\frac {j \pi}{2 (n + 1)}\right),
$$

for $j = 1 { : } n$ . Note that $V _ { n } ^ { \left( { D D } \right) } = \mathrm { D S T } ( n )$ . It follows that ${ \mathcal { T } } _ { n } ^ { ( D D ) }$ has a fast eigenvalue decomposition.

The Dirichlet-Neumann Matrix

If j is an integer and $\theta = ( 2 j - 1 ) \pi / ( 2 n )$ , then $s _ { n + 1 } - s _ { n - 1 } \ = \ 2 s _ { 1 } c _ { n } = 0$ . It follows from (4.8.18) that

$$
\mathcal {T} _ {n} ^ {(D N)} \cdot s (\theta_ {j}) = 4 \sin^ {2} (\theta_ {j} / 2) \cdot s (\theta_ {j}), \qquad \theta_ {j} = \frac {(2 j - 1) \pi}{2 n},
$$

for j = 1:n. Thus, the columns of the matrix V (DN) n $j = 1 { : } n$ $V _ { n } ^ { \left( D N \right) } \in \mathbb { R } ^ { n \times n }$ defined by

$$
[ V _ {n} ^ {(D N)} ] _ {k j} = \sin \left(\frac {k (2 j - 1) \pi}{2 n}\right)
$$

are eigenvectors of the matrix ${ \mathcal { T } } _ { n } ^ { ( D N ) }$ and the corresponding eigenvalues are given by

$$
\lambda_ {j} = 4 \sin^ {2} \left(\frac {(2 j - 1) \pi}{4 n}\right)
$$

for $j = 1 { : } n$ . Comparing with (1.4.13) we see that that $V _ { n } ^ { ( D N ) } = { \mathrm { D S T 2 } } ( n )$ . The inverse DST2 can be evaluated fast. See Van Loan (FFT, p. 242) for details, but also P4.8.11. It follows that $\mathcal { T } ^ { ( D N ) }$ has a fast eigenvalue decomposition.

The Neumann-Neumann Matrix

If j is an integer and $\theta = ( j - 1 ) \pi / ( n - 1 )$ , then $c _ { n } - c _ { n - 2 } = - 2 s _ { 1 } s _ { n - 1 } = 0$ . It follows from (4.8.19) that

$$
\mathcal {T} _ {n} ^ {(N N)} \cdot c (\theta_ {j}) = 4 \sin^ {2} \left(\frac {\theta_ {j}}{2}\right) \cdot c (\theta_ {j}), \qquad \theta_ {j} = \frac {(j - 1) \pi}{n - 1}.
$$

Thus, the columns of the matrix $V _ { n } ^ { \left( D N \right) } \in \mathbb { R } ^ { n \times n }$ defined by

$$
[ V _ {n} ^ {(N N)} ] _ {k j} = \cos \left(\frac {(k - 1) (j - 1) \pi}{n - 1}\right)
$$

are eigenvectors of the matrix ${ \mathcal { T } } _ { n } ^ { ( D N ) }$ and the corresponding eigenvalues are given by

$$
\lambda_ {j} = 4 \sin^ {2} \left(\frac {(j - 1) \pi}{2 (n - 1)}\right)
$$

for $j = 1 { : } n$ . Comparing with (1.4.10) we see that

$$
V _ {n} ^ {(N N)} = \mathrm{DCT} (n) \cdot \mathrm{diag} (2, I _ {n - 2}, 2)
$$

and therefore $\mathcal { T } ^ { ( N N ) }$ has a fast eigenvalue decomposition.

The Periodic Matrix

We can proceed to work out the eigenvalue decomposition for ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ as we did in the previous three cases, i.e., by zeroing the residuals in (4.8.20) and (4.8.21). However, $\mathcal { T } _ { n } ^ { \left( P \right) }$ is a circulant matrix and so we know from Theorem 4.8.2 that

$$
F _ {n} ^ {- 1} \mathcal {T} _ {n} ^ {(P)} F _ {n} = \mathrm{diag} (\lambda_ {1}, \ldots , \lambda_ {n})
$$

where

$$
\lambda = \bar {F} _ {n} \left[ \begin{array}{c} 2 \\ - 1 \\ 0 \\ \vdots \\ - 1 \end{array} \right] = 2 \bar {F} _ {n} (:, 1) - \bar {F} _ {n} (:, 2) - \bar {F} _ {n} (:, n).
$$

It can be shown that

$$
\lambda_ {j} = 4 \sin^ {2} \left(\frac {(j - 1) \pi}{n}\right)
$$

for $j = 1 { : } n$ . It follows that ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ has a fast eigenvalue decomposition. However, since this matrix is real it is preferable to have a real V -matrix. Using the facts that

$$
\lambda_ {j} = \lambda_ {n + 2 - j} \tag {4.8.22}
$$

and

$$
\bar {F} _ {n} (:, j) = F _ {n} (:, (n + 2 - j)) \tag {4.8.23}
$$

for $j = 2 { : } n$ , it can be shown that if $m = { \mathsf { c e i l } } ( ( n + 1 ) / 2 )$ and

$$
V _ {n} ^ {(P)} = \left[ \operatorname{Re} \left(F _ {n} (:, 1: m) \mid \operatorname{Im} \left(F _ {n} (:, m + 1: n)\right) \right] \right. \tag {4.8.24}
$$

then

$$
\mathcal {T} _ {n} ^ {(P)} V _ {n} ^ {(P)} (:, j) = \lambda_ {j} V _ {n} ^ {(P)} (:, j) \tag {4.8.25}
$$

for $j = 1 { : } n$ . Manipulations with this real matrix and its inverse can be carried out rapidly as discussed in Van Loan (FFT, Chap. 4).

# 4.8.7 A Note on Symmetry and Boundary Conditions

In our presentation, the matrices ${ \mathcal { T } } _ { n } ^ { ( D N ) }$ ) and T (NN) n a ${ \mathcal { T } } _ { n } ^ { ( N N ) }$ re not symmetric. However, a simple diagonal similarity transformation changes this. For example, if $D = \mathrm { d i a g } ( I _ { n - 1 } , \sqrt { 2 } )$ , then $D ^ { - 1 } \mathcal { T } _ { n } ^ { ( D N ) } D$ is symmetric. Working with symmetric second difference matrices has certain attractions, i.e., the automatic orthogonality of the eigenvector matrix. See Strang (1999).

# Problems

P4.8.1 Suppose $z \in \mathbb { R } ^ { n }$ has the property that $z ( 2 : n ) = \mathcal { E } _ { n - 1 } z ( 2 : n )$ . Show that $C ( z )$ is symmetric and $\bar { F } _ { n } z$ is real.

P4.8.2 As measured in the Frobenius norm, what is the nearest real circulant matrix to a given real Toeplitz matrix?

P4.8.3 Given $x , z \in \mathbb { C } ^ { n }$ , show how to compute $y = C ( z ) \cdot x$ in O(n log n) flops. In this case, y is the cyclic convolution of x and z.

P4.8.4 Suppose $\begin{array} { l } { { a } } \end{array} = \left[ \begin{array} { l } { { a } _ { - n + 1 } , \ldots , { a } _ { - 1 } , { a } _ { 0 } , { a } _ { 1 } , \ldots , { a } _ { n - 1 } } \end{array} \right]$ and let $T = ( t _ { k j } )$ be the n-by-n Toeplitz matrix defined by $t _ { k j } = a _ { k - j }$ . Thus, if $\boldsymbol { a } = \left[ a _ { - 2 } , a _ { - 1 } , a _ { 0 } , a _ { 1 } , a _ { 2 } \right]$ , then

$$
T = T (a) = \left[ \begin{array}{c c c} a _ {0} & a _ {- 1} & a _ {- 2} \\ a _ {1} & a _ {0} & a _ {- 1} \\ a _ {2} & a _ {1} & a _ {0} \end{array} \right].
$$

It is possible to “embed” T into a circulant, e.g.,

$$
C = \left[ \begin{array}{c c c c c c c c} a _ {0} & a _ {- 1} & a _ {- 2} & 0 & 0 & 0 & a _ {2} & a _ {1} \\ a _ {1} & a _ {0} & a _ {- 1} & a _ {- 2} & 0 & 0 & 0 & a _ {2} \\ a _ {2} & a _ {1} & a _ {0} & a _ {- 1} & a _ {- 2} & 0 & 0 & 0 \\ 0 & a _ {2} & a _ {1} & a _ {0} & a _ {- 1} & a _ {- 2} & 0 & 0 \\ 0 & 0 & a _ {2} & a _ {1} & a _ {0} & a _ {- 1} & a _ {- 2} & 0 \\ 0 & 0 & 0 & a _ {2} & a _ {1} & a _ {0} & a _ {- 1} & a _ {- 2} \\ a _ {- 2} & 0 & 0 & 0 & a _ {2} & a _ {1} & a _ {0} & a _ {- 1} \\ a _ {- 1} & a _ {- 2} & 0 & 0 & 0 & a _ {2} & a _ {1} & a _ {0} \end{array} \right].
$$

Given $a _ { - n + 1 } , \ldots , a _ { - 1 } , 1 _ { 0 } , a _ { 1 } , \ldots , a _ { n - 1 }$ and $m \geq 2 n - 1$ , show how to construct a vector $v \in \mathbb { C } ^ { m }$ so that if $C = C ( v )$ , then $C ( 1 { : } n , 1 { : } n ) = T$ . Note that v is not unique if $m > 2 n - 1$ .

P4.8.5 Complete the proof of Lemma 4.8.3.

P4.8.6 Show how to compute a Toeplitz-vector product $y = T u$ in n log n time using the embedding idea outlined in the previous problem and the fact that circulant matrices have a fast eigenvalue decomposition.

P4.8.7 Give a complete specification of the vector b in (4.8.12) if $A _ { 1 } = \mathcal { T } _ { n _ { 1 } } ^ { \left( D D \right) } , A _ { 2 } = \mathcal { T } _ { n _ { 2 } } ^ { \left( D D \right) }$ D) , A2 = T (DD) n2 , and $u ( x , y ) = 0$ on the boundary of the rectangular domain R. In terms of the underlying grid, $n _ { 1 } = m _ { 1 } - 1$ and $n _ { 2 } = m _ { 2 } - 1$ .

P4.8.8 Give a complete specification of the vector b in (4.8.12) if $A _ { 1 } \ = \ T _ { n _ { 1 } } ^ { ( D N ) } , \ A _ { 2 } \ = \ T _ { n _ { 2 } } ^ { ( D N ) }$ , $u ( x , y ) = 0$ on the bottom and left edge of $R , u _ { x } ( x , y ) = 0$ along the right edge of R, and $u _ { y } ( x , y ) = 0$ along the top edge of R. In terms of the underlying grid, $n _ { 1 } = m _ { 1 }$ and $n _ { 2 } = m _ { 2 }$ .

P4.8.9 Define a Neumann-Dirichlet matrix ${ \mathcal { T } } _ { n } ^ { ( N D ) }$ that would arise in conjunction with (4.8.5) if $u ^ { \prime } ( \alpha )$ and $u ( \beta )$ were specified. Show that ${ \mathcal { T } } _ { n } ^ { ( N D ) }$ has a fast eigenvalue decomposition.

${ \mathcal { T } } _ { n } ^ { ( N N ) }$ and  woul ${ \mathcal { T } } _ { n } ^ { ( P ) }$ are singular. (a) Assum solve the linear system b is in the range of subject to the const $A =$ In2 ⊗ T (P ) n1 + $I _ { n _ { 2 } } \otimes \mathcal { T } _ { n _ { 1 } } ^ { ( P ) } + \mathcal { T } _ { n _ { 2 } } ^ { ( P ) } \otimes I _ { n _ { 1 } }$ $A u = b$ Repeat part (a) replacing $\mathcal { T } _ { n _ { 1 } } ^ { \left( P \right) }$ with $\mathcal { T } _ { n _ { 1 } } ^ { ( N N ) }$ and $\mathcal { T } _ { n _ { 2 } } ^ { ( P ) }$ with $\mathcal { T } _ { n _ { 2 } } ^ { ( N N ) }$ .

P4.8.11 Let V be the matrix that defines the $\mathrm { D S T 2 } ( n )$ transformation in (1.4.12). (a) Show that

$$
V ^ {T} V = \frac {n}{2} I _ {n} + \frac {1}{2} v v ^ {T}
$$

where $v = [ 1 , - 1 , 1 , . . . , ( - 1 ) ^ { n } ] ^ { T }$ . (b) Verify that

$$
V ^ {- 1} = \frac {2}{n} \left(I - \frac {1}{2 n} v v ^ {T}\right) V ^ {T}.
$$

(c) Show how to compute $V ^ { - 1 } x$ rapidly.

P4.8.12 Verify (4.8.22), (4.8.23), and (4.8.25).

P4.8.13 Show that if $V = V _ { 2 m } ^ { \left( P \right) }$ V (P ) defined in (4.8.24), then

$$
V ^ {T} V = m \left(I _ {n} + e _ {1} e _ {1} ^ {T} + e _ {m + 1} e _ {m + 1} ^ {T}\right).
$$

What can you say about $V ^ { T } V { \mathrm { ~ i f ~ } } V = V _ { 2 m - 1 } ^ { ( P ) } ?$

# Notes and References for §4.8

As we mentioned, this section is based on Van Loan (FFT). For more details about fast Poisson solvers, see:

R.W. Hockney (1965). “A Fast Direct Solution of Poisson’s Equation Using Fourier Analysis,” J. Assoc. Comput. Mach. 12, 95–113.   
B. Buzbee, G. Golub, and C. Nielson (1970). “On Direct Methods for Solving Poisson’s Equation,” SIAM J. Numer. Anal. 7, 627–656.   
F. Dorr (1970). “The Direct Solution of the Discrete Poisson Equation on a Rectangle,” SIAM Review 12, 248–263.   
R. Sweet (1973). “Direct Methods for the Solution of Poisson’s Equation on a Staggered Grid,” J. Comput. Phys. 12, 422–428.   
P.N. Swarztrauber (1974). “A Direct Method for the Discrete Solution of Separable Elliptic Equations,” SIAM J. Numer. Anal. 11, 1136–1150.   
P.N. Swarztrauber (1977). “The Methods of Cyclic Reduction, Fourier Analysis and Cyclic Reduction-Fourier Analysis for the Discrete Solution of Poisson’s Equation on a Rectangle,” SIAM Review 19, 490–501.

There are actually eight variants of the discrete cosine transform each of which corresponds to the location of the Neumann conditions and how the divided difference approximations are set up. For a unified, matrix-based treatment, see:

G. Strang (1999). “The Discrete Cosine Transform,” SIAM Review 41, 135–147.
