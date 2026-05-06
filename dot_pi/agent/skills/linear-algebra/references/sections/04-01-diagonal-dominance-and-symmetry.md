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
