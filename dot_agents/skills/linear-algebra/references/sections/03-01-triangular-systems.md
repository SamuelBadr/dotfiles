# 3.1 Triangular Systems

Traditional factorization methods for linear systems involve the conversion of the given square system to a triangular system that has the same solution. This section is about the solution of triangular systems.

# 3.1.1 Forward Substitution

Consider the following 2-by-2 lower triangular system:

$$
\left[ \begin{array}{c c} \ell_ {1 1} & 0 \\ \ell_ {2 1} & \ell_ {2 2} \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right] = \left[ \begin{array}{c} b _ {1} \\ b _ {2} \end{array} \right].
$$

If $\ell _ { 1 1 } \ell _ { 2 2 } \neq 0$ , then the unknowns can be determined sequentially:

$$
x _ {1} = b _ {1} / \ell_ {1 1},
$$

$$
x _ {2} = (b _ {2} - \ell_ {2 1} x _ {1}) / \ell_ {2 2}.
$$

This is the 2-by-2 version of an algorithm known as forward substitution. The general procedure is obtained by solving the ith equation in $L x = b$ for $x _ { i } )$ :

$$
x _ {i} = \left. \left(b _ {i} - \sum_ {j = 1} ^ {i - 1} \ell_ {i j} x _ {j}\right) \right/ \ell_ {i i}.
$$

If this is evaluated for i = 1:n, then a complete specification of x is obtained. Note that at the ith stage the dot product of L(i, 1:i − 1) and x(1:i − 1) is required. Since $b _ { i }$ is involved only in the formula for $x _ { i }$ , the former may be overwritten by the latter.

Algorithm 3.1.1 (Row-Oriented Forward Substitution) If $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ is lower triangular and b ∈ IRn, then this algorithm overwrites b with the solution to $L x = b$ . L is assumed to be nonsingular.

$$
b (1) = b (1) / L (1, 1)
$$

for i = 2:n

$$
b (i) = (b (i) - L (i, 1: i - 1) \cdot b (1: i - 1)) / L (i, i)
$$

end

This algorithm requires $n ^ { 2 }$ flops. Note that L is accessed by row. The computed solution ˆx can be shown to satisfy

$$
(L + F) \hat {x} = b \quad | F | \leq n \mathbf {u} | L | + O \left(\mathbf {u} ^ {2}\right). \tag {3.1.1}
$$

For a proof, see Higham (ASNA, pp. 141-142). It says that the computed solution exactly satisfies a slightly perturbed system. Moreover, each entry in the perturbing matrix F is small relative to the corresponding element of L.

# 3.1.2 Back Substitution

The analogous algorithm for an upper triangular system $U x = b$ is called back substitution. The recipe for $x _ { i }$ is prescribed by

$$
x _ {i} = \left. \left(b _ {i} - \sum_ {j = i + 1} ^ {n} u _ {i j} x _ {j}\right) \right/ u _ {i i}
$$

and once again $b _ { i }$ can be overwritten by $x _ { i } .$ .

Algorithm 3.1.2 (Row-Oriented Back Substitution) If $U \in \mathbb { R } ^ { n \times n }$ is upper triangular and $b \in \mathbb { R } ^ { n }$ , then the following algorithm overwrites b with the solution to $U x = b$ . U is assumed to be nonsingular.

$$
b (n) = b (n) / U (n, n)
$$

for $i = n - 1 \colon - 1 \colon 1$

$$
b (i) = (b (i) - U (i, i + 1: n) \cdot b (i + 1: n)) / U (i, i)
$$

end

This algorithm requires $n ^ { 2 }$ flops and accesses $U$ by row. The computed solution ˆx obtained by the algorithm can be shown to satisfy

$$
(U + F) \hat {x} = b, \quad | F | \leq n \mathbf {u} | U | + O \left(\mathbf {u} ^ {2}\right). \tag {3.1.2}
$$

# 3.1.3 Column-Oriented Versions

Column-oriented versions of the above procedures can be obtained by reversing loop orders. To understand what this means from the algebraic point of view, consider forward substitution. Once $x _ { 1 }$ is resolved, it can be removed from equations 2 through n leaving us with the reduced system

$$
L (2: n, 2: n) x (2: n) = b (2: n) - x (1) \cdot L (2: n, 1).
$$

We next compute $x _ { 2 }$ and remove it from equations 3 through n, etc. Thus, if this approach is applied to

$$
\left[ \begin{array}{c c c} 2 & 0 & 0 \\ 1 & 5 & 0 \\ 7 & 9 & 8 \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \\ x _ {3} \end{array} \right] = \left[ \begin{array}{c} 6 \\ 2 \\ 5 \end{array} \right],
$$

we find $x _ { 1 } = 3$ and then deal with the 2-by-2 system

$$
\left[ \begin{array}{l l} 5 & 0 \\ 9 & 8 \end{array} \right] \left[ \begin{array}{l} x _ {2} \\ x _ {3} \end{array} \right] = \left[ \begin{array}{l} 2 \\ 5 \end{array} \right] - 3 \left[ \begin{array}{l} 1 \\ 7 \end{array} \right] = \left[ \begin{array}{l} - 1 \\ - 1 6 \end{array} \right].
$$

Here is the complete procedure with overwriting.

Algorithm 3.1.3 (Column-Oriented Forward Substitution) If the matrix $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ is lower triangular and $b \in \mathbb { R } ^ { n }$ , then this algorithm overwrites b with the solution to $L x = b $ . L is assumed to be nonsingular.

for $j = 1 \colon n - 1$

$$
b (j) = b (j) / L (j, j)
$$

$$
b (j + 1: n) = b (j + 1: n) - b (j) \cdot L (j + 1: n, j)
$$

end

$$
b (n) = b (n) / L (n, n)
$$

It is also possible to obtain a column-oriented saxpy procedure for back substitution.

Algorithm 3.1.4 (Column-Oriented Back Substitution) If $U \in \mathbb { R } ^ { n \times n }$ is upper triangular and $b \in \mathbb { R } ^ { n }$ , then this algorithm overwrites b with the solution to $U x = b$ . U is assumed to be nonsingular.

for $j = n \colon - 1 \colon 2$

$$
b (j) = b (j) / U (j, j)
$$

$$
b (1: j - 1) = b (1: j - 1) - b (j) \cdot U (1: j - 1, j)
$$

end

$$
b (1) = b (1) / U (1, 1)
$$

Note that the dominant operation in both Algorithms 3.1.3 and 3.1.4 is the saxpy operation. The roundoff behavior of these implementations is essentially the same as for the dot product versions.

# 3.1.4 Multiple Right-Hand Sides

Consider the problem of computing a solution $X \in \mathbb { R } ^ { n \times q }$ to $L X = B$ where $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ is lower triangular and $B \in \mathbb { R } ^ { n \times q }$ . This is the multiple-right-hand-side problem and it amounts to solving q separate triangular systems, i.e., $L X ( : , j ) = B ( : , j ) , j = 1 : q ,$ . Interestingly, the computation can be blocked in such a way that the resulting algorithm is rich in matrix multiplication, assuming that q and n are large enough. This turns out to be important in subsequent sections where various block factorization schemes are discussed.

It is sufficient to consider just the lower triangular case as the derivation of block back substitution is entirely analogous. We start by partitioning the equation $L X = B$ as follows:

$$
\left[ \begin{array}{c c c c} L _ {1 1} & 0 & \dots & 0 \\ L _ {2 1} & L _ {2 2} & \dots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ L _ {N 1} & L _ {N 2} & \dots & L _ {N N} \end{array} \right] \left[ \begin{array}{c} X _ {1} \\ X _ {2} \\ \vdots \\ X _ {N} \end{array} \right] = \left[ \begin{array}{c} B _ {1} \\ B _ {2} \\ \vdots \\ B _ {N} \end{array} \right]. \tag {3.1.3}
$$

Assume that the diagonal blocks are square. Paralleling the development of Algorithm 3.1.3, we solve the system $L _ { 1 1 } X _ { 1 } = B _ { 1 }$ for $X _ { 1 }$ and then remove $X _ { 1 }$ from block equations 2 through N:

$$
\left[ \begin{array}{c c c c} L _ {2 2} & 0 & \dots & 0 \\ L _ {3 2} & L _ {3 3} & \dots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ L _ {N 2} & L _ {N 3} & \dots & L _ {N N} \end{array} \right] \left[ \begin{array}{c} X _ {2} \\ X _ {3} \\ \vdots \\ X _ {N} \end{array} \right] = \left[ \begin{array}{c} B _ {2} \\ B _ {3} \\ \vdots \\ B _ {N} \end{array} \right] - \left[ \begin{array}{c} L _ {2 1} \\ L _ {3 1} \\ \vdots \\ L _ {N 1} \end{array} \right] X _ {1}.
$$

Continuing in this way we obtain the following block forward elimination scheme:

for $j = 1:N$ Solve $L_{jj}X_j = B_j$ for $i = j + 1:N$ $B_i = B_i - L_{ij}X_j$ (3.1.4) end end

Notice that the i-loop oversees a single block saxpy update of the form

$$
\left[ \begin{array}{c} B _ {j + 1} \\ \vdots \\ B _ {N} \end{array} \right] = \left[ \begin{array}{c} B _ {j + 1} \\ \vdots \\ B _ {N} \end{array} \right] - \left[ \begin{array}{c} L _ {j + 1, j} \\ \vdots \\ L _ {N, j} \end{array} \right] X _ {j}.
$$

To realize level-3 performance, the submatrices in (3.1.3) must be sufficiently large in dimension.

# 3.1.5 The Level-3 Fraction

It is handy to adopt a measure that quantifies the amount of matrix multiplication in a given algorithm. To this end we define the level-3 fraction of an algorithm to be the fraction of flops that occur in the context of matrix multiplication. We call such flops level-3 flops.

Let us determine the level-3 fraction for (3.1.4) with the simplifying assumption that $n = r N$ . (The same conclusions hold with the unequal blocking described above.) Because there are N applications of r-by-r forward elimination (the level-2 portion of the computation) and $\mathbf { \hat { \eta } } _ { n } \mathbf { \hat { \eta } } _ { n } ^ { 2 }$ flops overall, the level-3 fraction is approximately given by

$$
1 - \frac {N r ^ {2}}{n ^ {2}} = 1 - \frac {1}{N}.
$$

Thus, for large N almost all flops are level-3 flops. It makes sense to choose N as large as possible subject to the constraint that the underlying architecture can achieve a high level of performance when processing block saxpys that have width $r = n / N$ or greater.

# 3.1.6 Nonsquare Triangular System Solving

The problem of solving nonsquare, m-by-n triangular systems deserves some attention. Consider the lower triangular case when $m \geq n$ , i.e.,

$$
\left[ \begin{array}{l} L _ {1 1} \\ L _ {2 1} \end{array} \right] x \quad = \quad \left[ \begin{array}{l} b _ {1} \\ b _ {2} \end{array} \right] \qquad \begin{array}{l} L _ {1 1} \in \mathbb {R} ^ {n \times n}, \qquad b _ {1} \in \mathbb {R} ^ {n}, \\ L _ {2 1} \in \mathbb {R} ^ {(m - n) \times n}, \quad b _ {2} \in \mathbb {R} ^ {m - n}. \end{array}
$$

Assume that $L _ { 1 1 }$ is lower triangular and nonsingular. If we apply forward elimination to $L _ { 1 1 } x = b _ { 1 }$ , then x solves the system provided $L _ { 2 1 } ( L _ { 1 1 } ^ { - 1 } b _ { 1 } ) = b _ { 2 }$ . Otherwise, there is no solution to the overall system. In such a case least squares minimization may be appropriate. See Chapter 5.

Now consider the lower triangular system Lx = b when the number of columns n exceeds the number of rows m. We can apply forward substitution to the square system $L ( 1 { : } m , 1 { : } m ) x ( 1 { : } m , 1 { : } m ) = b$ and prescribe an arbitrary value for $x ( m + 1 { : } n )$ . See §5.6 for additional comments on systems that have more unknowns than equations. The handling of nonsquare upper triangular systems is similar. Details are left to the reader.

# 3.1.7 The Algebra of Triangular Matrices

A unit triangular matrix is a triangular matrix with 1’s on the diagonal. Many of the triangular matrix computations that follow have this added bit of structure. It clearly poses no difficulty in the above procedures.

For future reference we list a few properties about products and inverses of triangular and unit triangular matrices.

• The inverse of an upper (lower) triangular matrix is upper (lower) triangular.   
• The product of two upper (lower) triangular matrices is upper (lower) triangular.   
• The inverse of a unit upper (lower) triangular matrix is unit upper (lower) triangular.   
• The product of two unit upper (lower) triangular matrices is unit upper (lower) triangular.

# Problems

P3.1.1 Give an algorithm for computing a nonzero $z \in \mathbb { R } ^ { n }$ such that $U z = 0$ where $U \in \mathbb { R } ^ { n \times n }$ is upper triangular with $u _ { n n } = 0$ and $u _ { 1 1 } \cdot \cdot \cdot u _ { n - 1 , n - 1 } \neq 0$ .

P3.1.2 Suppose $L = I _ { n } - N$ is unit lower triangular where $N \in \mathbb { R } ^ { n \times n }$ . Show that

$$
L ^ {- 1} = I _ {n} + N + N ^ {2} + \dots + N ^ {n - 1}.
$$

What is the value of $\parallel L ^ { - 1 } \parallel _ { F } \mathrm { i f } \ N _ { i j } = 1$ for all $i > j ?$

P3.1.3 Write a detailed version of (3.1.4). Do not assume that N divides $n .$

P3.1.4 Prove all the facts about triangular matrices that are listed in 3.1.7.

P3.1.5 Suppose $S , T \in \mathbb { R } ^ { n \times n }$ are upper triangular and that $( S T - \lambda I ) x = b$ is a nonsingular system. Give an $O ( n ^ { 2 } )$ algorithm for computing x. Note that the explicit formation of $S T - \lambda I$ requires $O ( n ^ { 3 } )$ flops. Hint: Suppose

$$
S _ {+} = \left[ \begin{array}{l l} \sigma & u ^ {T} \\ 0 & S _ {c} \end{array} \right], \quad T _ {+} = \left[ \begin{array}{l l} \tau & v ^ {T} \\ 0 & T _ {c} \end{array} \right], \quad b _ {+} = \left[ \begin{array}{l} \beta \\ b _ {c} \end{array} \right],
$$

where $S _ { + } = S ( k - 1 { : } n , k - 1 { : } n ) , T _ { + } = T ( k - 1 { : } n , k - 1 { : } n ) , b _ { + } = b ( k - 1 { : } n )$ , and $\sigma , \tau , \beta \in \mathbb { R }$ . Show that if we have a vector $x _ { c }$ such that

$$
(S _ {c} T _ {c} - \lambda I) x _ {c} = b _ {c}
$$

and $w _ { c } = T _ { c } x _ { c }$ is available, then

$$
x _ {+} = \left[ \begin{array}{l} \gamma \\ x _ {c} \end{array} \right], \quad \gamma = \frac {\beta - \sigma v ^ {T} x _ {c} - u ^ {T} w _ {c}}{\sigma \tau - \lambda}
$$

solves $( S _ { + } T _ { + } - \lambda I ) x _ { + } = b _ { + }$ . Observe that $x _ { + }$ and $w _ { + } = T _ { + } x _ { + }$ each require $O ( n - k )$ flops.

P3.1.6 Suppose the matrices $R 1 , \ldots , R _ { p } \in \mathbb { R } ^ { n \times n }$ are all upper triangular. Give an $O ( p n ^ { 2 } )$ algorithm for solving the system $( R _ { 1 } \cdot \cdot \cdot R _ { p } - \lambda I ) x = b$ assuming that the matrix of coefficients is nonsingular. Hint. Generalize the solution to the previous problem.

P3.1.7 Suppose $L , K \in \mathbb { R } ^ { n \times n }$ are lower triangular and $B \in \mathbb { R } ^ { n \times n }$ . Give an algorithm for computing $X \in \mathbb { R } ^ { n \times n }$ so that $L X K = B$ .

# Notes and References for §3.1

The accuracy of a computed solution to a triangular system is often surprisingly good, see:

N.J. Higham (1989). “The Accuracy of Solutions to Triangular Systems,” SIAM J. Numer. Anal. 26, 1252–1265.

Solving systems of the form $( T _ { p } \cdot \cdot \cdot T _ { 1 } - \lambda I ) x = b$ where each $T _ { i }$ is triangular is considered in:

C.D. Martin and C.F. Van Loan (2002). “Product Triangular Systems with Shift,” SIAM J. Matrix Anal. Applic. 24, 292–301.

The trick to obtaining an $O ( p n ^ { 2 } )$ procedure that does not involve any matrix-matrix multiplications is to look carefully at the back-substitution recursions. See P3.1.6.

A survey of parallel triangular system solving techniques and their stabilty is given in:

N.J. Higham (1995). “Stability of Parallel Triangular System Solvers,” SIAM J. Sci. Comput. 16, 400–413.
