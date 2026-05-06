# Chapter 3

# General Linear Systems

3.1 Triangular Systems   
3.2 The LU Factorization   
3.3 Roundoff Error in Gaussian Elimination   
3.4 Pivoting   
3.5 Improving and Estimating Accuracy   
3.6 Parallel LU

The problem of solving a linear system Ax = b is central to scientific computation. In this chapter we focus on the method of Gaussian elimination, the algorithm of choice if A is square, dense, and unstructured. Other methods are applicable if A does not fall into this category, see Chapter 4, Chapter 11, §12.1, and §12.2. Solution procedures for triangular systems are discussed first. These are followed by a derivation of Gaussian elimination that makes use of Gauss transformations. The process of eliminating unknowns from equations is described in terms of the factorization $A = L U$ where L is lower triangular and U is upper triangular. Unfortunately, the derived method behaves poorly on a nontrivial class of problems. An error analysis pinpoints the difficulty and sets the stage for a discussion of pivoting, a permutation strategy that keeps the numbers “nice” during the elimination. Practical issues associated with scaling, iterative improvement, and condition estimation are covered. A framework for computing the LU factorization in parallel is developed in the final section.

# Reading Notes

Familiarity with Chapter 1, §§2.1–2.5, and §2.7 is assumed. The sections within this chapter depend upon each other as follows:

$$
\begin{array}{c c c c c c c c} & & & & & \S 3. 5 \\ & & & & & \uparrow \\ \S 3. 1 & \to & \S 3. 2 & \to & \S 3. 3 & \to & \S 3. 4 \\ & & & & & \downarrow \\ & & & & & \S 3. 6 \end{array}
$$

Useful global references include Forsythe and Moler (SLAS), Stewart( MABD), Higham (ASNA), Watkins (FMC), Trefethen and Bau (NLA), Demmel (ANLA), and Ipsen (NMA).

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

# 3.2 The LU Factorization

Triangular system solving is an easy $O ( n ^ { 2 } )$ computation. The idea behind Gaussian elimination is to convert a given system Ax = b to an equivalent triangular system. The conversion is achieved by taking appropriate linear combinations of the equations. For example, in the system

$$
3 x _ {1} + 5 x _ {2} = 9,
$$

$$
6 x _ {1} + 7 x _ {2} = 4,
$$

if we multiply the first equation by 2 and subtract it from the second we obtain

$$
3 x _ {1} + 5 x _ {2} = \quad 9,
$$

$$
- 3 x _ {2} = - 1 4.
$$

This is n = 2 Gaussian elimination. Our objective in this section is to describe the procedure in the language of matrix factorizations. This means showing that the algorithm computes a unit lower triangular matrix L and an upper triangular matrix U so that A = LU , e.g.,

$$
{\left[ \begin{array}{l l} 3 & 5 \\ 6 & 7 \end{array} \right]} = {\left[ \begin{array}{l l} 1 & 0 \\ 2 & 1 \end{array} \right]} {\left[ \begin{array}{l l} 3 & 5 \\ 0 & - 3 \end{array} \right]}.
$$

The solution to the original Ax = b problem is then found by a two-step triangular solve process:

$$
L y = b, \quad U x = y \quad \Longrightarrow \quad A x = L U x = L y = b. \tag {3.2.1}
$$

The LU factorization is a “high-level” algebraic description of Gaussian elimination. Linear equation solving is not about the matrix vector product $A ^ { - 1 } b$ but about computing LU and using it effectively; see §3.4.9. Expressing the outcome of a matrix algorithm in the “language” of matrix factorizations is a productive exercise, one that is repeated many times throughout this book. It facilitates generalization and highlights connections between algorithms that can appear very different at the scalar level.

# 3.2.1 Gauss Transformations

To obtain a factorization description of Gaussian elimination as it is traditionally presented, we need a matrix description of the zeroing process. At the $n = 2$ level, if $v _ { 1 } \neq 0$ and $\tau = v _ { 2 } / v _ { 1 }$ , then

$$
\left[ \begin{array}{c c} 1 & 0 \\ - \tau & 1 \end{array} \right] \left[ \begin{array}{c} v _ {1} \\ v _ {2} \end{array} \right] = \left[ \begin{array}{c} v _ {1} \\ 0 \end{array} \right].
$$

More generally, suppose $v \in \mathbb { R } ^ { n }$ with $v _ { k } \neq 0$ . If

$$
\tau^ {T} = [ \underbrace {0 , \ldots , 0} _ {k}, \tau_ {k + 1}, \ldots , \tau_ {n} ], \qquad \tau_ {i} = \frac {v _ {i}}{v _ {k}}, \quad i = k + 1: n,
$$

and we define

$$
M _ {k} = I _ {n} - \tau e _ {k} ^ {T}, \tag {3.2.2}
$$

then

$$
M _ {k} v = \left[ \begin{array}{c c c c c c} 1 & \dots & 0 & 0 & \dots & 0 \\ \vdots & \ddots & \vdots & \vdots & & \vdots \\ 0 & & 1 & 0 & & 0 \\ 0 & & - \tau_ {k + 1} & 1 & & 0 \\ \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & \dots & - \tau_ {n} & 0 & \dots & 1 \end{array} \right] \left[ \begin{array}{c} v _ {1} \\ \vdots \\ v _ {k} \\ v _ {k + 1} \\ \vdots \\ v _ {n} \end{array} \right] = \left[ \begin{array}{c} v _ {1} \\ \vdots \\ v _ {k} \\ 0 \\ \vdots \\ 0 \end{array} \right].
$$

A matrix of the form $M _ { k } = I _ { n } - \tau e _ { k } ^ { T } \in \mathbb { R } ^ { n \times n }$ is a Gauss transformation if the first k components of $\tau \in \mathbb { R } ^ { n }$ are zero. Such a matrix is unit lower triangular. The components of $\tau ( k + 1 { : } n )$ are called multipliers. The vector τ is called the Gauss vector.

# 3.2.2 Applying Gauss Transformations

Multiplication by a Gauss transformation is particularly simple. If $C \in \mathbb { R } ^ { n \times r }$ and $M _ { k } = I _ { n } - \tau e _ { k } ^ { T }$ is a Gauss transformation, then

$$
M _ {k} C = (I _ {n} - \tau e _ {k} ^ {T}) C = C - \tau (e _ {k} ^ {T} C) = C - \tau C (k,:) \nonumber
$$

is an outer product update. Since $\tau ( 1 { : } k ) = 0$ only $C ( k + 1 { : } n , : )$ is affected and the update $C = M _ { k } C$ can be computed row by row as follows:

for $i = k + 1 { : } n$

$$
C (i,:) = C (i,:) - \tau_ {i} \cdot C (k,:) \tag {1}
$$

end

This computation requires $2 ( n - k ) r$ flops. Here is an example:

$$
C   =   \left[ \begin{array}{c c c} 1 & 4 & 7 \\ 2 & 5 & 8 \\ 3 & 6 & 1 0 \end{array} \right],    \tau   =   \left[ \begin{array}{c} 0 \\ 1 \\ - 1 \end{array} \right] \qquad \Longrightarrow \qquad (I - \tau e _ {1} ^ {T}) C = \left[ \begin{array}{c c c} 1 & 4 & 7 \\ 1 & 1 & 1 \\ 4 & 1 0 & 1 7 \end{array} \right].
$$

# 3.2.3 Roundoff Properties of Gauss Transformations

If ˆτ is the computed version of an exact Gauss vector τ , then it is easy to verify that

$$
\hat {\tau} = \tau + e, \qquad | e | \leq \mathbf {u} | \tau |.
$$

If ˆτ is used in a Gauss transform update and $\mathrm { H } ( ( I _ { n } - \hat { \tau } e _ { k } ^ { T } ) C )$ denotes the computed result, then

$$
\mathrm{fl} \left((I _ {n} - \hat {\tau} e _ {k} ^ {T}) C\right) = (I - \tau e _ {k} ^ {T}) C + E  ,
$$

where

$$
| E | \leq 3 \mathbf {u} (| C | + | \tau | | C (k,:) |) + O (\mathbf {u} ^ {2}).
$$

Clearly, if τ has large components, then the errors in the update may be large in comparison to |C|. For this reason, care must be exercised when Gauss transformations are employed, a matter that is pursued in §3.4.

# 3.2.4 Upper Triangularizing

Assume that $A \in \mathbb { R } ^ { n \times n }$ . Gauss transformations $M _ { 1 } , \dots , M _ { n - 1 }$ can usually be found such that $M _ { n - 1 } \cdot \cdot \cdot M _ { 2 } M _ { 1 } A = U$ is upper triangular. To see this we first look at the $n = 3$ case. Suppose

$$
A = \left[ \begin{array}{l l l} 1 & 4 & 7 \\ 2 & 5 & 8 \\ 3 & 6 & 1 0 \end{array} \right]
$$

and note that

$$
M _ {1} = \left[ \begin{array}{r r r} 1 & 0 & 0 \\ - 2 & 1 & 0 \\ - 3 & 0 & 1 \end{array} \right] \quad \Rightarrow \quad M _ {1} A = \left[ \begin{array}{r r r} 1 & 4 & 7 \\ 0 & - 3 & - 6 \\ 0 & - 6 & - 1 1 \end{array} \right].
$$

Likewise, in the second step we have

$$
M _ {2} = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & - 2 & 1 \end{array} \right] \quad \Rightarrow \quad M _ {2} (M _ {1} A) = \left[ \begin{array}{c c c} 1 & 4 & 7 \\ 0 & - 3 & - 6 \\ 0 & 0 & 1 \end{array} \right].
$$

Extrapolating from this example to the general n case we conclude two things.

• At the start of the kth step we have a matrix $A ^ { ( k - 1 ) } = M _ { k - 1 } \cdot \cdot \cdot M _ { 1 } A$ that is upper triangular in columns 1 through k − 1.   
• The multipliers in the kth Gauss transform $M _ { k }$ are based on $A ^ { ( k - 1 ) } ( k + 1 { : } n , k )$ and a(k−kk $a _ { k k } ^ { ( k - 1 ) }$ must be nonzero in order to proceed.

Noting that complete upper triangularization is achieved after n − 1 steps, we obtain the following rough draft of the overall process:

$$
A ^ {(1)} = A
$$

for $k = 1 { : } n - 1$

For $i = k + 1 { : } n$ , determine the multipliers $\tau _ { i } ^ { ( k ) } = a _ { i k } ^ { ( k ) } / a _ { k k } ^ { ( k ) }$ = a ik /akk . (3.2.3)

Apply $M _ { k } = I - \tau ^ { ( k ) } e _ { k } ^ { T }$ to obtain $A ^ { ( k + 1 ) } = M _ { k } A ^ { ( k ) }$ .

end

For this process to be well-defined, the matrix entries $a _ { 1 1 } ^ { ( 1 ) } , a _ { 2 2 } ^ { ( 2 ) } , \dots , a _ { n - 1 , n - 1 } ^ { ( n - 1 ) }$ must be nonzero. These quantities are called pivots.

# 3.2.5 Existence

If no zero pivots are encountered in (3.2.3), then Gauss transformations $M _ { 1 } , \dots , M _ { n - 1 }$ are generated such that $M _ { n - 1 } \cdot \cdot \cdot M _ { 1 } A \ = \ U$ is upper triangular. It is easy to check that if $M _ { k } = I _ { n } - \tau ^ { ( k ) } e _ { k } ^ { T }$ − , then its inverse is prescribed by $\Breve { M _ { k } ^ { - 1 } } = I _ { n } + \tau ^ { ( k ) } \Breve { e _ { k } ^ { T } }$ and so

$$
A = L U \tag {3.2.4}
$$

where

$$
L = M _ {1} ^ {- 1} \dots M _ {n - 1} ^ {- 1}. \tag {3.2.5}
$$

It is clear that L is a unit lower triangular matrix because each $M _ { k } ^ { - 1 }$ is unit lower triangular. The factorization (3.2.4) is called the LU factorization.

The LU factorization may not exist. For example, it is impossible to find $l _ { i j }$ and $u _ { i j }$ so

$$
{\left[ \begin{array}{l l l} 1 & 2 & 3 \\ 2 & 4 & 7 \\ 3 & 5 & 3 \end{array} \right]} = {\left[ \begin{array}{l l l} 1 & 0 & 0 \\ \ell_ {2 1} & 1 & 0 \\ \ell_ {3 1} & \ell_ {3 2} & 1 \end{array} \right]} {\left[ \begin{array}{l l l} u _ {1 1} & u _ {1 2} & u _ {1 3} \\ 0 & u _ {2 2} & u _ {2 3} \\ 0 & 0 & u _ {3 3} \end{array} \right]}.
$$

To see this, equate entries and observe that we must have $u _ { 1 1 } = 1 , u _ { 1 2 } = 2 , \ell _ { 2 1 } = 2$ , $u _ { 2 2 } = 0$ , and $\ell _ { 3 1 } = 3$ . But then the (3,2) entry gives us the contradictory equation $5 = \ell _ { 3 1 } u _ { 1 2 } + \ell _ { 3 2 } u _ { 2 2 } = 6$ . For this example, the pivot $a _ { 2 2 } ^ { ( 1 ) } = a _ { 2 2 } - ( a _ { 2 1 } / a _ { 1 1 } ) a _ { 1 2 }$ is zero.

It turns out that the kth pivot in (3.2.3) is zero if $A ( 1 { : } k , 1 { : } k )$ is singular. A submatrix of the form A(1:k, 1:k) is called a leading principal submatrix.

Theorem 3.2.1. (LU Factorization). If $A \in \mathbb { R } ^ { n \times n }$ and det $( A ( 1 { : } k , 1 { : } k ) ) \neq 0$ for $k = 1 { : } n - 1$ , then there exists a unit lower triangular $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ and an upper triangular $U \in \mathbb { R } ^ { n \times n }$ such that $A \ = \ L U$ . If this is the case and A is nonsingular, then the factorization is unique and det $( A ) = u _ { 1 1 } \cdot \cdot \cdot u _ { n n }$ .

Proof. Suppose k − 1 steps in (3.2.3) have been executed. At the beginning of step k the matrix A has been overwritten by $M _ { k - 1 } \cdot \cdot \cdot M _ { 1 } A = A ^ { ( k - 1 ) }$ . Since Gauss transformations are unit lower triangular, it follows by looking at the leading k-by-k portion of this equation that

$$
\det (A (1: k, 1: k)) = a _ {1 1} ^ {(k - 1)} \dots a _ {k k} ^ {(k - 1)}. \tag {3.2.6}
$$

Thus, if $A ( 1 { : } k , 1 { : } k )$ is nonsingular, then the kth pivot akk $a _ { k k } ^ { ( k - 1 ) }$ is nonzero.

As for uniqueness, if $A = L _ { 1 } U _ { 1 }$ and $A = L _ { 2 } U _ { 2 }$ are two LU factorizations of a nonsingular A, then $L _ { 2 } ^ { - 1 } L _ { 1 } = U _ { 2 } U _ { 1 } ^ { - 1 }$ . Since $L _ { 2 } ^ { - 1 } L _ { 1 }$ is unit lower triangular and $U _ { 2 } U _ { 1 } ^ { - 1 }$ is upper triangular, it follows that both of these matrices must equal the identity. Hence, $L _ { 1 } = L _ { 2 }$ and $U _ { 1 } = U _ { 2 }$ . Finally, if $A = L U$ , then

$$
\det (A) = \det (L U) = \det (L) \det (U) = \det (U).
$$

It follows that det $( A ) \ = \ u _ { 1 1 } \cdot \cdot \cdot u _ { n n }$

□

# 3.2.6 L Is the Matrix of Multipliers

It turns out that the construction of L is not nearly so complicated as Equation (3.2.5) suggests. Indeed,

$$
\begin{array}{l} L = M _ {1} ^ {- 1} \dots M _ {n - 1} ^ {- 1} \\ = \left(I _ {n} - \tau^ {(1)} e _ {1} ^ {T}\right) ^ {- 1} \dots \left(I _ {n} - \tau^ {(n - 1)} e _ {n - 1} ^ {T}\right) ^ {- 1} \\ = \left(I _ {n} + \tau^ {(1)} e _ {1} ^ {T}\right) \dots \left(I _ {n} + \tau^ {(n - 1)} e _ {n - 1} ^ {T}\right) \\ = I _ {n} + \sum_ {k = 1} ^ {n - 1} \tau^ {(k)} e _ {k} ^ {T} \\ \end{array}
$$

showing that

$$
L (k + 1: n, k) = \tau^ {(k)} (k + 1: n) \quad k = 1: n - 1. \tag {3.2.7}
$$

In other words, the kth column of L is defined by the multipliers that arise in the k-th step of (3.2.3). Consider the example in §3.2.4:

$$
\tau^ {(1)} = \left[ \begin{array}{l} 0 \\ 2 \\ 3 \end{array} \right], \tau^ {(2)} = \left[ \begin{array}{l} 0 \\ 0 \\ 2 \end{array} \right] \quad \Rightarrow \quad \left[ \begin{array}{l l l} 1 & 4 & 7 \\ 2 & 5 & 8 \\ 3 & 6 & 1 0 \end{array} \right] = \left[ \begin{array}{l l l} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 3 & 2 & 1 \end{array} \right] \left[ \begin{array}{l l l} 1 & 4 & 7 \\ 0 & - 3 & - 6 \\ 0 & 0 & 1 \end{array} \right].
$$

# 3.2.7 The Outer Product Point of View

Since the application of a Gauss transformation to a matrix involves an outer product, we can regard (3.2.3) as a sequence of outer product updates. Indeed, if

$$
A = \left[ \begin{array}{c c} \alpha & w ^ {T} \\ v & B \\ 1 & n - 1 \end{array} \right] _ {n - 1} ^ {1}
$$

then the first step in Gaussian elimination results in the decomposition

$$
\left[ \begin{array}{c c} \alpha & w ^ {T} \\ z & B \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ z / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ 0 & B - z w ^ {T} / \alpha \end{array} \right] \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & I _ {n - 1} \end{array} \right].
$$

Steps 2 through n − 1 compute the LU factorization

$$
B - z w ^ {T} / \alpha = L _ {1} U _ {1}
$$

for then

$$
A = \left[ \begin{array}{l l} 1 & 0 \\ z / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{l l} 1 & 0 \\ 0 & L _ {1} U _ {1} \end{array} \right] \left[ \begin{array}{l l} \alpha & w ^ {T} \\ 0 & I _ {n - 1} \end{array} \right] = \left[ \begin{array}{l l} 1 & 0 \\ z / \alpha & L _ {1} \end{array} \right] \left[ \begin{array}{l l} \alpha & w ^ {T} \\ 0 & U _ {1} \end{array} \right] \equiv L U.
$$

# 3.2.8 Practical Implementation

Let us consider the efficient implementation of (3.2.3). First, because zeros have already been introduced in columns 1 through $k - 1$ , the Gauss transformation update need only be applied to columns k through n. Of course, we need not even apply the kth Gauss transform to $A ( : , k )$ since we know the result. So the efficient thing to do is simply to update $A ( k + 1 { : } n , k + 1 { : } n )$ . Also, the observation (3.2.7) suggests that we can overwrite $A ( k + 1 { : } n , k )$ with $L ( k + 1 { : } n , k )$ since the latter houses the multipliers that are used to zero the former. Overall we obtain:

Algorithm 3.2.1 (Outer Product LU) Suppose $A \in \mathbb { R } ^ { n \times n }$ has the property that $A ( 1 { : } k , 1 { : } k )$ is nonsingular for $k = 1 { : } n - 1$ . This algorithm computes the factorization $A = L U$ where L is unit lower triangular and U is upper triangular. For $i = 1 { : } n - 1$ , $A ( i , i { : } n )$ is overwritten by $U ( i , i ; n )$ while $A ( i + 1 { : } n , i )$ is overwritten by $L ( i + 1 { : } n , i )$ .

for $k = 1 { : } n - 1$

$$
\rho = k + 1: n
$$

$$
A (\rho , k) = A (\rho , k) / A (k, k)
$$

$$
A (\rho , \rho) = A (\rho , \rho) - A (\rho , k) \cdot A (k, \rho)
$$

end

This algorithm involves $2 n ^ { 3 } / 3$ flops and it is one of several formulations of Gaussian elimination. Note that the k-th step involves an $( n - k ) – \mathrm { b y } – ( n - k )$ outer product.

# 3.2.9 Other Versions

Similar to matrix-matrix multiplication, Gaussian elimination is a triple-loop procedure that can be arranged in several ways. Algorithm 3.2.1 corresponds to the $\ " k i j \ "$ version of Gaussian elimination if we compute the outer product update row by row:

for k = 1:n - 1
    A(k + 1:n, k) = A(k + 1:n, k)/A(k, k)
    for i = k + 1:n
    for j = k + 1:n
    A(i, j) = A(i, j) - A(i, k)·A(k, j)
    end
    end
end

There are five other versions: $k j i , i k j , i j k , j i k$ , and $j k i$ . The last of these results in an implementation that features a sequence of gaxpys and forward eliminations which we now derive at the vector level.

The plan is to compute the jth columns of L and $U$ in step j. If $j = 1$ , then by comparing the first columns in $A = L U$ we conclude that

$$
L (2: n, j) = A (2: n, 1) / A (1, 1)
$$

and $U ( 1 , 1 ) = A ( 1 , 1 )$ . Now assume that $L ( : , 1 : j - 1 )$ and $U ( 1 { : } j - 1 , 1 { : } j - 1 )$ are known. To get the jth columns of $L$ and U we equate the jth columns in the equation $A = L U$ and infer from the vector equation $A ( : , j ) = L U ( : , j )$ that

$$
A (1: j - 1, j) = L (1: j - 1, 1: j - 1) \cdot U (1: j - 1, j)
$$

and

$$
A (j: n, j) = \sum_ {k = 1} ^ {j} L (j: n, k) \cdot U (k, j).
$$

The first equation is a lower triangular linear system that can be solved for the vector $U ( 1 ; j - 1 , j )$ . Once this is accomplished, the second equation can be rearranged to produce recipes for $U ( j , j )$ and $L ( j + 1 { : } n , j )$ . Indeed, if we set

$$
\begin{array}{l} v (j: n) = A (j: n, j) - \sum_ {k = 1} ^ {j - 1} L (j: n, k) U (k, j) \\ = A (j: n, j) - L (j: n, 1: j - 1) \cdot U (1: j - 1, j), \\ \end{array}
$$

then $L ( j + 1 ; n , j ) = v ( j + 1 ; n ) / v ( j )$ and $U ( j , j ) = v ( j )$ . Thus, $L ( j + 1 { : } n , j )$ is a scaled gaxpy and we obtain the following alternative to Algorithm 3.2.1:

Algorithm 3.2.2 (Gaxpy LU) Suppose $A \in \mathbb { R } ^ { n \times n }$ has the property that $A ( 1 { : } k , 1 { : } k )$ i s nonsingular for $k = 1 { : } n - 1$ . This algorithm computes the factorization $A = L U$ where L is unit lower triangular and U is upper triangular.

Initialize L to the identity and U to the zero matrix.

for $j = 1:n$ if $j = 1$ $v = A(:,1)$ else $\tilde{a} = A(:,j)$ Solve $L(1:j-1,1:j-1) \cdot z = \tilde{a}(1:j-1)$ for $z \in \mathbb{R}^{j-1}$ . $U(1:j-1,j) = z$ $v(j:n) = \tilde{a}(j:n) - L(j:n,1:j-1) \cdot z$ end $U(j,j) = v(j)$ $L(j+1:n,j) = v(j+1:n)/v(j)$

end

(We chose to have separate arrays for L and U for clarity; it is not necessary in practice.) Algorithm 3.2.2 requires $2 n ^ { 3 } / 3$ flops, the same volume of floating point work required by Algorithm 3.2.1. However, from §1.5.2 there is less memory traffic associated with a gaxpy than with an outer product, so the two implementations could perform differently in practice. Note that in Algorithm 3.2.2, the original $A ( : , j )$ is untouched until step j.

The terms right-looking and left-looking are sometimes applied to Algorithms 3.2.1 and 3.2.2. In the outer-product implementation, after $L ( k { : } n , k )$ is determined, the columns to the right of $A ( : , k )$ are updated so it is a right-looking procedure. In contrast, subcolumns to the left of $A ( : , k )$ are accessed in gaxpy LU before $L ( k + 1 { : } n , k )$ is produced so that implementation left-looking.

# 3.2.10 The LU Factorization of a Rectangular Matrix

The LU factorization of a rectangular matrix $A \in \mathbb { R } ^ { n \times r }$ can also be performed. The $n > r$ case is illustrated by

$$
{\left[ \begin{array}{l l} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{array} \right]} = {\left[ \begin{array}{l l} 1 & 0 \\ 3 & 1 \\ 5 & 2 \end{array} \right]} {\left[ \begin{array}{l l} 1 & 2 \\ 0 & - 2 \end{array} \right]}
$$

while

$$
{\left[ \begin{array}{l l l} 1 & 2 & 3 \\ 4 & 5 & 6 \end{array} \right]} = {\left[ \begin{array}{l l} 1 & 0 \\ 4 & 1 \end{array} \right]} {\left[ \begin{array}{l l l} 1 & 2 & 3 \\ 0 & - 3 & - 6 \end{array} \right]}
$$

depicts the $n < r$ situation. The LU factorization of $A \in \mathbb { R } ^ { n \times r }$ is guaranteed to exist if $A ( 1 { : } k , 1 { : } k )$ is nonsingular for $k = 1 { : } \mathrm { m i n } \{ n , r \}$ .

The square LU factorization algorithms above needs only minor alterations to handle the rectangular case. For example, if $n > r$ , then Algorithm 3.2.1 modifies to the following:

for $k = 1 { : } r$

$$
\rho = k + 1: n
$$

$$
A (\rho , k) = A (\rho , k) / A (k, k)
$$

$\mathbf { i f } \ k < r$

$$
\mu = k + 1: r
$$

$$
A (\rho , \mu) = A (\rho , \mu) - A (\rho , k) \cdot A (k, \mu)
$$

end

end

This calculation requires $n r ^ { 2 } - r ^ { 3 } / 3$ flops. Upon completion, A is overwritten by the strictly lower triangular portion of $\ b { L } \in \mathbb { R } ^ { n \times r }$ and the upper triangular portion of $U \in \mathbb { R } ^ { r \times r }$ .

# 3.2.11 Block LU

It is possible to organize Gaussian elimination so that matrix multiplication becomes the dominant operation. Partition $A \in \mathbb { R } ^ { n \times n }$ as follows:

$$
\begin{array}{r l} {A =} & {\left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right] _ {n - r} ^ {r}} \\ & {r \quad n - r} \end{array}
$$

where $r$ is a blocking parameter. Suppose we compute the LU factorization

$$
\left[ \begin{array}{l} A _ {1 1} \\ A _ {2 1} \end{array} \right] = \left[ \begin{array}{l} L _ {1 1} \\ L _ {2 1} \end{array} \right] U _ {1 1}.
$$

Here, $L _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ is unit lower triangular and $U _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ is upper triangular and assumed to be nonsingular. If we solve $L _ { 1 1 } U _ { 1 2 } = A _ { 1 2 }$ for $U _ { 1 2 } \in \mathbb { R } ^ { r \times n - r }$ , then

$$
\left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right] = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & I _ {n - r} \end{array} \right] \left[ \begin{array}{c c} I _ {r} & 0 \\ 0 & \tilde {A} \end{array} \right] \left[ \begin{array}{c c} U _ {1 1} & U _ {1 2} \\ 0 & I _ {n - r} \end{array} \right],
$$

where

$$
\tilde {A} = A _ {2 2} - L _ {2 1} U _ {1 2} = A _ {2 2} - A _ {2 1} A _ {1 1} ^ {- 1} A _ {1 2} \tag {3.2.9}
$$

is the Schur complement of $A _ { 1 1 }$ in A. Note that if

$$
\tilde {A} = L _ {2 2} U _ {2 2}
$$

is the LU factorization of ${ \tilde { A } } .$ then

$$
A = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \end{array} \right] \left[ \begin{array}{c c} U _ {1 1} & U _ {1 2} \\ 0 & U _ {2 2} \end{array} \right]
$$

is the LU factorization of A. This lays the groundwork for a recursive implementation.

Algorithm 3.2.3 (Recursive Block LU) Suppose $A \in \mathbb { R } ^ { n \times n }$ has an LU factorization and r is a positive integer. The following algorithm computes unit lower triangular $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ and upper triangular $U \in \mathbb { R } ^ { n \times n }$ so $A = L U$ .

function $[ L , U ] = { \mathsf { B l o c k L U } } ( A , n , r )$

$n \leq r$

Compute the LU factorization A = LU using (say) Algorithm 3.2.1.

else

Use (3.2.8) to compute the LU factorization $A ( : , 1 : r ) = \left[ \begin{array} { l } { L _ { 1 1 } } \\ { L _ { 2 1 } } \end{array} \right] U _ { 1 1 }$ .

Solve $L _ { 1 1 } U _ { 1 2 } = A ( 1 { : } r , r + 1 { : } n )$ for $U _ { 1 2 }$ .

$$
\tilde {A} = A (r + 1: n, r + 1: n) - L _ {2 1} U _ {1 2}
$$

$$
\left[ L _ {2 2}, U _ {2 2} \right] = \operatorname{BlockLU} (\tilde {A}, n - r, r)
$$

$$
L = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \end{array} \right], U = \left[ \begin{array}{c c} U _ {1 1} & U _ {1 2} \\ 0 & U _ {2 2} \end{array} \right]
$$

end

The following table explains where the flops come from:

<table><tr><td>Activity</td><td>Flops</td></tr><tr><td> $L_{11}, L_{21}, U_{11}$ </td><td> $nr^{2} - r^{3}/3$ </td></tr><tr><td> $U_{12}$ </td><td> $(n - r)r^{2}$ </td></tr><tr><td> $\tilde{A}$ </td><td> $2(n - r)^{2}$ </td></tr></table>

If $n \gg r ,$ , then there are a total of about $2 n ^ { 3 } / 3$ flops, the same volume of atithmetic as Algorithms 3.2.1 and 3.2.2. The vast majority of these flops are the level-3 flops associated with the production of ${ \tilde { A } } .$ .

The actual level-3 fraction, a concept developed in §3.1.5, is more easily derived from a nonrecursive implementation. Assume for clarity that $n = N r$ where N is a positive integer and that we want to compute

$$
\left[ \begin{array}{c c c} A _ {1 1} & \dots & A _ {1 N} \\ \vdots & \ddots & \vdots \\ A _ {N 1} & \dots & A _ {N N} \end{array} \right] = \left[ \begin{array}{c c c} L _ {1 1} & \dots & 0 \\ \vdots & \ddots & \vdots \\ L _ {N 1} & \dots & L _ {N N} \end{array} \right] \left[ \begin{array}{c c c} U _ {1 1} & \dots & U _ {1 N} \\ \vdots & \ddots & \vdots \\ 0 & \dots & U _ {N N} \end{array} \right] \tag {3.2.10}
$$

where all blocks are $r { \mathrm { - } } \mathrm { b y } { \mathrm { - } } r$ . Analogously to Algorithm 3.2.3 we have the following.

Algorithm 3.2.4 (Nonrecursive Block LU) Suppose $A \in \mathbb { R } ^ { n \times n }$ has an $L U$ factorization and r is a positive integer. The following algorithm computes unit lower triangular $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ and upper triangular $U \in \mathbb { R } ^ { n \times n }$ so $A = L U$ .

for $k = 1 { : } N$

Rectangular Gaussian elimination:

$$
\left[ \begin{array}{c} A _ {k k} \\ \vdots \\ A _ {N k} \end{array} \right] = \left[ \begin{array}{c} L _ {k k} \\ \vdots \\ L _ {N k} \end{array} \right] U _ {k k}
$$

Multiple right hand side solve:

$$
L _ {k k} \left[ \begin{array}{c c c c} U _ {k, k + 1} & \ldots & U _ {k N} \end{array} \right] = \left[ \begin{array}{c c c c} A _ {k, k + 1} & \ldots & A _ {k N} \end{array} \right]
$$

Level-3 updates:

$$
A _ {i j} = A _ {i j} - L _ {i k} U _ {k j}, \quad i = k + 1: N, j = k + 1: N
$$

end

Here is the flop situation during the kth pass through the loop:

<table><tr><td>Activity</td><td>Flops</td></tr><tr><td>Gaussian elimination</td><td> $(N - k + 1)r^{3} - r^{3}/3$ </td></tr><tr><td>Multiple RHS solve</td><td> $(N - k)r^{3}$ </td></tr><tr><td>Level-3 updates</td><td> $2(N - k)^{2}r^{2}$ </td></tr></table>

Summing these quantities for $k = 1 { : } N$ we find that the level-3 fraction is approximately

$$
{\frac {2 n ^ {3} / 3}{2 n ^ {3} / 3 + n ^ {2} r}} = 1 - {\frac {3}{2 N}}.
$$

Thus, for large N almost all arithmetic takes place in the context of matrix multiplication. This ensures a favorable amount of data reuse as discussed in §1.5.4.

# Problems

P3.2.1 Verify Equation (3.2.6).

P3.2.2 Suppose the entries of $A ( \epsilon ) \in \mathbb { R } ^ { n \times n }$ are continuously differentiable functions of the scalar -. Assume that $A \equiv A ( 0 )$ and all its principal submatrices are nonsingular. Show that for sufficiently small -, the matrix $A ( \epsilon )$ has an LU factorization $A ( \epsilon ) = L ( \epsilon ) U ( \epsilon )$ and that $L ( \epsilon )$ and $U ( \epsilon )$ are both continuously differentiable.

P3.2.3 Suppose we partition $A \in \mathbb { R } ^ { n \times n }$

$$
A = \left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right]
$$

where $A _ { 1 1 }$ is $r { \mathrm { - b y } } - r$ and nonsingular. Let S be the Schur complement of $A _ { 1 1 }$ in A as defined in (3.2.9). Show that after r steps of Algorithm 3.2.1, $A ( r + 1 { : } n , r + 1 { : } n )$ houses S. How could S be obtained after r steps of Algorithm 3.2.2?

P3.2.4 Suppose $A \in \mathbb { R } ^ { n \times n }$ has an LU factorization. Show how Ax = b can be solved without storing the multipliers by computing the LU factorization of the n-by-(n + 1) matrix [A b].

P3.2.5 Describe a variant of Gaussian elimination that introduces zeros into the columns of A in the order, $n \colon - 1 { : } 2$ and which produces the factorization $A = U L$ where U is unit upper triangular and L is lower triangular.

P3.2.6 Matrices in $\mathbb { R } ^ { n \times n }$ of the form $N ( y , k ) = I - y e _ { k } ^ { T }$ where $\boldsymbol { y } \in \mathbb { R } ^ { n }$ are called Gauss-Jordan transformations. (a) Give a formula for $N ( y , k ) ^ { - 1 }$ assuming it exists. (b) Given $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , under what conditions can y be found so $N ( y , k ) x = e _ { k } ? \ ( \mathrm { c } )$ Give an algorithm using Gauss-Jordan transformations that overwrites A with $A ^ { - 1 }$ . What conditions on A ensure the success of your algorithm?

P3.2.7 Extend Algorithm 3.2.2 so that it can also handle the case when A has more rows than columns.

P3.2.8 Show how A can be overwritten with L and U in Algorithm 3.2.2. Give a 3-loop specification so that unit stride access prevails.

P3.2.9 Develop a version of Gaussian elimination in which the innermost of the three loops oversees a dot product.

# Notes and References for §3.2

The method of Gaussian elimination has a long and interesting history, see:

J.F. Grcar (2011). “How Ordinary Elimination Became Gaussian Elimination,” Historica Mathematica, 38, 163–218.   
J.F. Grcar (2011). “Mathematicians of Gaussian Elimination,” Notices of the AMS 58, 782–792.   
Schur complements (3.2.9) arise in many applications. For a survey of both practical and theoretical interest, see:   
R.W. Cottle (1974). “Manifestations of the Schur Complement,” Lin. Alg. Applic. 8, 189–211.   
Schur complements are known as “Gauss transforms” in some application areas. The use of Gauss-Jordan transformations (P3.2.6) is detailed in Fox (1964). See also:   
T. Dekker and W. Hoffman (1989). “Rehabilitation of the Gauss-Jordan Algorithm,” Numer. Math. 54, 591–599.   
As we mentioned, inner product versions of Gaussian elimination have been known and used for some time. The names of Crout and Doolittle are associated with these techniques, see:   
G.E. Forsythe (1960). “Crout with Pivoting,” Commun. ACM 3, 507–508.   
W.M. McKeeman (1962). “Crout with Equilibration and Iteration,” Commun. ACM. 5, 553–555.   
Loop orderings and block issues in LU computations are discussed in:   
J.J. Dongarra, F.G. Gustavson, and A. Karp (1984). “Implementing Linear Algebra Algorithms for Dense Matrices on a Vector Pipeline Machine,” SIAM Review 26, 91–112.   
J.M. Ortega (1988). “The ijk Forms of Factorization Methods I: Vector Computers,” Parallel Comput. 7, 135–147.   
D.H. Bailey, K.Lee, and H.D. Simon (1991). “Using Strassen’s Algorithm to Accelerate the Solution of Linear Systems,” J. Supercomput. 4, 357–371.   
J.W. Demmel, N.J. Higham, and R.S. Schreiber (1995). “Stability of Block LU Factorization,” Numer. Lin. Alg. Applic. 2, 173–190.   
Suppose $A = L U$ and $A + \Delta A = ( L + \Delta L ) ( U + \Delta U )$ are LU factorizations. Bounds on the perturbations ∆L and ∆U in terms of ∆A are given in:   
G.W. Stewart (1997). “On the Perturbation of LU and Cholesky Factors,” IMA J. Numer. Anal. 17, 1–6.   
X.-W. Chang and C.C. Paige (1998). “On the Sensitivity of the LU factorization,” BIT 38, 486–501.

In certain limited domains, it is possible to solve linear systems exactly using rational arithmetic. For a snapshot of the challenges, see:

P. Alfeld and D.J. Eyre (1991). “The Exact Analysis of Sparse Rectangular Linear Systems,” ACM Trans. Math. Softw. 17, 502–518.   
P. Alfeld (2000). “Bivariate Spline Spaces and Minimal Determining Sets,” J. Comput. Appl. Math. 119, 13–27.

# 3.3 Roundoff Error in Gaussian Elimination

We now assess the effect of rounding errors when the algorithms in the previous two sections are used to solve the linear system $A x = b$ . A much more detailed treatment of roundoff error in Gaussian elimination is given in Higham (ASNA).

# 3.3.1 Errors in the LU Factorization

Let us see how the error bounds for Gaussian elimination compare with the ideal bounds derived in §2.7.11. We work with the infinity norm for convenience and focus our attention on Algorithm 3.2.1, the outer product version. The error bounds that we derive also apply to the gaxpy formulation (Algorithm 3.2.2). Our first task is to quantify the roundoff errors associated with the computed triangular factors.

Theorem 3.3.1. Assume that A is an n-by-n matrix of floating point numbers. If no zero pivots are encountered during the execution of Algorithm 3.2.1, then the computed triangular matrices $\hat { L }$ and Uˆ satisfy

$$
\hat {L} \hat {U} = A + H, \tag {3.3.1}
$$

$$
| H | \leq 2 (n - 1) \mathbf {u} \left(| A | + | \hat {L} | | \hat {U} |\right) + O \left(\mathbf {u} ^ {2}\right). \tag {3.3.2}
$$

Proof. The proof is by induction on n. The theorem obviously holds for $n = 1$ . Assume that $n \geq 2$ and that the theorem holds for all $( n - 1 ) – \mathrm { b y } – ( n - 1 )$ floating point matrices. If A is partitioned as follows

$$
A = \left[ \begin{array}{c c} \alpha & w ^ {T} \\ v & B \end{array} \right] _ {n - 1} ^ {1}
$$

then the first step in Algorithm 3.2.1 is to compute

$$
\hat {z} = \mathsf {f l} (v / \alpha), \qquad \hat {C} = \mathsf {f l} (\hat {z} w ^ {T}), \qquad \hat {A} _ {1} = \mathsf {f l} (B - \hat {C}),
$$

from which we conclude that

$$
\hat {z} = v / \alpha + f, \tag {3.3.3}
$$

$$
| f | \leq \mathbf {u} | v / \alpha |, \tag {3.3.4}
$$

$$
\hat {C} = \hat {z} w ^ {T} + F _ {1}, \tag {3.3.5}
$$

$$
\left| F _ {1} \right| \leq \mathbf {u} \left| \hat {z} \right| \left| w ^ {T} \right|, \tag {3.3.6}
$$

$$
\hat {A} _ {1} = B - \left(\hat {z} w ^ {T} + F _ {1}\right) + F _ {2}, \tag {3.3.7}
$$

$$
\left| F _ {2} \right| \leq \mathbf {u} \left(\left| B \right| + \left| \hat {z} \right| \left| w ^ {T} \right|\right) + O \left(\mathbf {u} ^ {2}\right), \tag {3.3.8}
$$

$$
| \hat {A} _ {1} | \leq | B | + | \hat {z} | | w ^ {T} | + O (\mathbf {u}). \tag {3.3.9}
$$

The algorithm proceeds to compute the LU factorization of $\hat { A } _ { 1 }$ . By induction, the computed factors $\hat { L } _ { 1 }$ and $\hat { U } _ { 1 }$ satisfy

$$
\hat {L} _ {1} \hat {U} _ {1} = \hat {A} _ {1} + H _ {1} \tag {3.3.10}
$$

where

$$
\left| H _ {1} \right| \leq 2 (n - 2) \mathbf {u} \left(\left| \hat {A} _ {1} \right| + \left| \hat {L} _ {1} \right| \left| \hat {U} _ {1} \right|\right) + O \left(\mathbf {u} ^ {2}\right). \tag {3.3.11}
$$

If

$$
\hat {L} = \left[ \begin{array}{c c} 1 & 0 \\ \hat {z} & \hat {L} _ {1} \end{array} \right], \qquad \hat {U} = \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & \hat {U} _ {1} \end{array} \right],
$$

then it is easy to verify that

$$
\hat {L} \hat {U} = A + H
$$

where

$$
H = \left[ \begin{array}{c c} 0 & 0 \\ \alpha f & H _ {1} - F _ {1} + F _ {2} \end{array} \right]. \tag {3.3.12}
$$

To prove the theorem we must verify (3.3.2), i.e.,

$$
| H | \leq 2 (n - 1) \mathbf {u} \left[ \begin{array}{c c} 2 | \alpha | & 2 | w ^ {T} | \\ | v | + | \alpha | | f | & | B | + | \hat {L} _ {1} | | \hat {U} _ {1} | + | \hat {z} | | w ^ {T} | \end{array} \right] + O (\mathbf {u} ^ {2}).
$$

Considering (3.3.12), this is obviously the case if

$$
\left| H _ {1} \right| + \left| F _ {1} \right| + \left| F _ {2} \right| \leq 2 (n - 1) \mathbf {u} \left(\left| B \right| + \left| \hat {z} \right| \left| w ^ {T} \right| + \left| \hat {L} _ {1} \right| \left| \hat {U} _ {1} \right|\right) + O \left(\mathbf {u} ^ {2}\right). \tag {3.3.13}
$$

Using (3.3.9) and (3.3.11) we have

$$
| H _ {1} | \leq 2 (n - 2) \mathbf {u} \left(| B | + | \hat {z} | | w ^ {T} | + | \hat {L} _ {1} | | \hat {U} _ {1} |\right) + O (\mathbf {u} ^ {2}),
$$

while (3.3.6) and (3.3.8) imply

$$
| F _ {1} | + | F _ {2} | \leq \mathbf {u} (| B | + 2 | \hat {z} | | w |) + O (\mathbf {u} ^ {2}).
$$

These last two results establish (3.3.13) and therefore the theorem.

We mention that if A is m-by-n, then the theorem applies with n replaced by the smaller of n and m in Equation 3.3.2.

# 3.3.2 Triangular Solving with Inexact Triangles

We next examine the effect of roundoff error when $\hat { L }$ and $\hat { U }$ are used by the triangular system solvers of §3.1.

Theorem 3.3.2. Let $\hat { L }$ and $\hat { U }$ be the computed LU factors obtained by Algorithm 3.2.1 when it is applied to an n-by-n floating point matrix A. If the methods of §3.1 are used to produce the computed solution yˆ to Lyˆ = b and the computed solution xˆ to $\hat { U } x = \hat { y }$ , then $( { \boldsymbol { A } } + { \boldsymbol { E } } ) { \hat { \boldsymbol { x } } } = { \boldsymbol { b } }$ with

$$
| E | \leq n \mathbf {u} (2 | A | + 4 | \hat {L} | | \hat {U} |) + O \left(\mathbf {u} ^ {2}\right). \tag {3.3.14}
$$

Proof. From (3.1.1) and (3.1.2) we have

$$
(\hat {L} + F) \hat {y} = b, \quad | F | \leq n \mathbf {u} | \hat {L} | + O (\mathbf {u} ^ {2}),
$$

$$
(\hat {U} + G) \hat {x} = \hat {y}, \quad | G | \leq n \mathbf {u} | \hat {U} | + O (\mathbf {u} ^ {2}),
$$

and thus

$$
(\hat {L} + F) (\hat {U} + G) \hat {x} = (\hat {L} \hat {U} + F \hat {U} + \hat {L} G + F G) \hat {x} = b.
$$

If follows from Theorem 3.3.1 that $\hat { L } \hat { U } = A + H$ with

$$
| H | \leq 2 (n - 1) \mathbf {u} (| A | + | \hat {L} | | \hat {U} |) + O (\mathbf {u} ^ {2}),
$$

and so by defining

$$
E = H + F \hat {U} + \hat {L} G + F G
$$

we find (A + E)ˆx = b. Moreover,

$$
| E | \leq | H | + | F | | \hat {U} | + | \hat {L} | | G | + O (\mathbf {u} ^ {2})
$$

$$
\leq 2 n \mathbf {u} \left(| A | + | \hat {L} | | \hat {U} |\right) + 2 n \mathbf {u} \left(| \hat {L} | | \hat {U} |\right) + O (\mathbf {u} ^ {2}),
$$

completing the proof of the theorem.

If it were not for the possibility of a large $| \hat { L } | | \hat { U } |$ term, (3.3.14) would compare favorably with the ideal bound (2.7.21). (The factor n is of no consequence, cf. the Wilkinson quotation in §2.7.7.) Such a possibility exists, for there is nothing in Gaussian elimination to rule out the appearance of small pivots. If a small pivot is encountered, then we can expect large numbers to be present in $\hat { L }$ and $\hat { U }$ .

We stress that small pivots are not necessarily due to ill-conditioning as the example

$$
A = \left[ \begin{array}{c c} \epsilon & 1 \\ 1 & 0 \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ 1 / \epsilon & 1 \end{array} \right] \left[ \begin{array}{c c} \epsilon & 1 \\ 0 & - 1 / \epsilon \end{array} \right]
$$

shows. Thus, Gaussian elimination can give arbitrarily poor results, even for wellconditioned problems. The method is unstable. For example, suppose 3-digit floating point arithmetic is used to solve

$$
\left[ \begin{array}{l l}. 0 0 1 & 1. 0 0 \\ 1. 0 0 & 2. 0 0 \end{array} \right] \left[ \begin{array}{l} x _ {1} \\ x _ {2} \end{array} \right] = \left[ \begin{array}{l} 1. 0 0 \\ 3. 0 0 \end{array} \right].
$$

(See §2.7.1.) Applying Gaussian elimination we get

$$
\hat {L} = \left[ \begin{array}{c c} 1 & 0 \\ 1 0 0 0 & 1 \end{array} \right], \qquad \hat {U} = \left[ \begin{array}{c c}. 0 0 1 & 1 \\ 0 & - 1 0 0 0 \end{array} \right],
$$

and a calculation shows that

$$
\hat {L} \hat {U} = \left[ \begin{array}{c c}. 0 0 1 & 1 \\ 1 & 2 \end{array} \right] + \left[ \begin{array}{c c} 0 & 0 \\ 0 & - 2 \end{array} \right] \equiv A + H.
$$

If we go on to solve the problem using the triangular system solvers of §3.1, then using the same precision arithmetic we obtain a computed solution $\hat { x } = [ 0 \ , \ 1 ] ^ { T }$ . This is in contrast to the exact solution $\boldsymbol { x } = [ 1 . 0 0 2 \ldots , . 9 9 8 \ldots ] ^ { T }$ .

# Problems

P3.3.1 Show that if we drop the assumption that A is a floating point matrix in Theorem 3.3.1, then Equation 3.3.2 holds with the coefficient “2”replaced by “3.”

P3.3.2 Suppose A is an n-by-n matrix and that $\hat { L }$ and $\hat { U }$ are produced by Algorithm 3.2.1. (a) How many flops are required to compute $\| | \hat { L } | | \hat { U } | \| _ { \infty } ?$ (b) Show $\mathsf { f l } \big ( \big | \hat { L } \big | \big | \hat { U } \big | \big ) \leq ( 1 + 2 n \mathbf { u } ) \big | \hat { L } \big | \big | \hat { U } \big | + O \big ( \mathbf { u } ^ { 2 } \big )$ .

# Notes and References for 3.3

The original roundoff analysis of Gaussian elimination appears in:

J.H. Wilkinson (1961). “Error Analysis of Direct Methods of Matrix Inversion,” J. ACM 8, 281–330.

Various improvements and insights regarding the bounds and have been made over the years, see:

B.A. Chartres and J.C. Geuder (1967). “Computable Error Bounds for Direct Solution of Linear Equations,” J. ACM 14, 63–71.

J.K. Reid (1971). “A Note on the Stability of Gaussian Elimination,” J. Inst. Math. Applic. 8, 374–75.

C.C. Paige (1973). “An Error Analysis of a Method for Solving Matrix Equations,” Math. Comput. 27, 355–59.

H.H. Robertson (1977). “The Accuracy of Error Estimates for Systems of Linear Algebraic Equations,” J. Inst. Math. Applic. 20, 409–14.

J.J. Du Croz and N.J. Higham (1992). “Stability of Methods for Matrix Inversion,” IMA J. Numer. Anal. 12, 1–19.

J.M. Banoczi, N.C. Chiu, G.E. Cho, and I.C.F. Ipsen (1998). “The Lack of Influence of the Right–Hand Side on the Accuracy of Linear System Solution,” SIAM J. Sci. Comput. 20, 203–227.

P. Amodio and F. Mazzia (1999). “A New Approach to Backward Error Analysis of LU Factorization BIT 39, 385–402.

An interesting account of von Neuman’s contributions to the numerical analysis of Gaussian elimination is detailed in:

J.F. Grcar (2011). “John von Neuman’s Analysis of Gaussian Elimination and the Origins of Modern Numerical Analysis,” SIAM Review 53, 607–682.

# 3.4 Pivoting

The analysis in the previous section shows that we must take steps to ensure that no large entries appear in the computed triangular factors $\hat { L }$ and $\hat { U }$ . The example

$$
A = \left[ \begin{array}{c c}. 0 0 0 1 & 1 \\ 1 & 1 \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ 1 0 0 0 0 & 1 \end{array} \right] \left[ \begin{array}{c c}. 0 0 0 1 & 1 \\ 0 & - 9 9 9 9 \end{array} \right] = L U
$$


---

<!-- golub_150_199 -->

correctly identifies the source of the difficulty: relatively small pivots. A way out of this difficulty is to interchange rows. For example, if P is the permutation

$$
P = \left[ \begin{array}{l l} 0 & 1 \\ 1 & 0 \end{array} \right]
$$

then

$$
P A = \left[ \begin{array}{c c} 1 & 1 \\ . 0 0 0 1 & 1 \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ . 0 0 0 1 & 1 \end{array} \right] \left[ \begin{array}{c c} 1 & 1 \\ 0 & . 9 9 9 9 \end{array} \right] = L U.
$$

Observe that the triangular factors have modestly sized entries.

In this section we show how to determine a permuted version of A that has a reasonably stable LU factorization. There are several ways to do this and they each corresponds to a different pivoting strategy. Partial pivoting, complete pivoting, and rook pivoting are considered. The efficient implementation of these strategies and their properties are discussed. We begin with a few comments about permutation matrices that can be used to swap rows or columns.

# 3.4.1 Interchange Permutations

The stabilizations of Gaussian elimination that are developed in this section involve data movements such as the interchange of two matrix rows. In keeping with our desire to describe all computations in “matrix terms,” we use permutation matrices to describe this process. (Now is a good time to review §1.2.8–§1.2.11.) Interchange permutations are particularly important. These are permutations obtained by swapping two rows in the identity, e.g.,

$$
\Pi = \left[ \begin{array}{c c c c} 0 & 0 & 0 & 1 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 1 & 0 & 0 & 0 \end{array} \right].
$$

Interchange permutations can be used to describe row and column swapping. If $\ b { A } \in \mathbb { R } ^ { 4 \times 4 }$ , then Π·A is A with rows 1 and 4 interchanged while A·Π is A with columns 1 and 4 swapped.

If $P = \Pi _ { m } \cdot \cdot \cdot \Pi _ { 1 }$ and each $\Pi _ { k }$ is the identity with rows k and $p i v ( k )$ interchanged, then $p i v ( 1 { : } m )$ encodes P . Indeed, $\boldsymbol { x } \in \mathbb { R } ^ { n }$ can be overwritten by $P x$ as follows:

$$
\begin{array}{l} \text { for } k = 1: m \\ x (k) \leftrightarrow x (p i v (k)) \\ \end{array}
$$

Here, the $^ { 6 6 }  ^ { \mathfrak { N } }$ notation means “swap contents.” Since each $\Pi _ { k }$ is symmetric, we have $P ^ { T } = \Pi _ { 1 } \cdot \cdot \cdot \Pi _ { m }$ . Thus, the piv representation can also be used to overwrite x with $P ^ { T } x \mathrm { : }$ :

$$
\begin{array}{l} \text { for } k = m: - 1: 1 \\ x (k) \leftrightarrow x (p i v (k)) \\ \end{array}
$$

We remind the reader that although no floating point arithmetic is involved in a permutation operation, permutations move data and have a nontrivial effect upon performance.

# 3.4.2 Partial Pivoting

Interchange permutations can be used in LU computations to guarantee that no multiplier is greater than 1 in absolute value. Suppose

$$
A = \left[ \begin{array}{c c c} 3 & 1 7 & 1 0 \\ 2 & 4 & - 2 \\ 6 & 1 8 & - 1 2 \end{array} \right].
$$

To get the smallest possible multipliers in the first Gauss transformation, we need $a _ { 1 1 }$ to be the largest entry in the first column. Thus, if $\Pi _ { 1 }$ is the interchange permutation

$$
\Pi_ {1} = \left[ \begin{array}{l l l} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{array} \right]
$$

then

$$
\Pi_ {1} A = \left[ \begin{array}{c c c} 6 & 1 8 & - 1 2 \\ 2 & 4 & - 2 \\ 3 & 1 7 & 1 0 \end{array} \right].
$$

It follows that

$$
M _ {1} = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ - 1 / 3 & 1 & 0 \\ - 1 / 2 & 0 & 1 \end{array} \right] \qquad \Longrightarrow \qquad M _ {1} \Pi_ {1} A = \left[ \begin{array}{c c c} 6 & 1 8 & - 1 2 \\ 0 & - 2 & 2 \\ 0 & 8 & 1 6 \end{array} \right].
$$

To obtain the smallest possible multiplier in $M _ { 2 }$ , we need to swap rows 2 and 3. Thus, if

$$
\Pi_ {2} = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{array} \right] \qquad \text {and} \qquad M _ {2} = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 1 / 4 & 1 \end{array} \right],
$$

then

$$
M _ {2} \Pi_ {2} M _ {1} \Pi_ {1} A = \left[ \begin{array}{c c c} 6 & 1 8 & - 1 2 \\ 0 & 8 & 1 6 \\ 0 & 0 & 6 \end{array} \right].
$$

For general n we have

for $k = 1 { : } n - 1$

Find an interchange permutation $\Pi _ { k } \in \mathbb { R } ^ { n \times n }$ that swaps

A(k, k) with the largest element in $\left| A ( k { : } n , k ) \right|$ .

$$
A = \Pi_ {k} A \tag {3.4.1}
$$

Determine the Gauss transformation $M _ { k } = I _ { n } - \tau ^ { ( k ) } e _ { k } ^ { T }$ such that if

v is the kth column of $M _ { k } A$ , then $v ( k + 1 { : } n ) = 0$ .

$$
A = M _ {k} A
$$

end

This particular row interchange strategy is called partial pivoting and upon completion, we have

$$
M _ {n - 1} \Pi_ {n - 1} \dots M _ {1} \Pi_ {1} A = U \tag {3.4.2}
$$

where U is upper triangular. As a consequence of the partial pivoting, no multiplier is larger than one in absolute value.

# 3.4.3 Where is L?

It turns out that (3.4.1) computes the factorization

$$
P A = L U \tag {3.4.3}
$$

where $P = \Pi _ { n - 1 } \cdot \cdot \cdot \Pi _ { 1 }$ , U is upper triangular, and L is unit lower triangular with $| \ell _ { i j } | \le 1$ . We show that $L ( k + 1 { : } n , k )$ is a permuted version of $M _ { k }$ ’s multipliers. From (3.4.2) it can be shown that

$$
\tilde {M} _ {n - 1} \dots \tilde {M} _ {1} P A = U \tag {3.4.4}
$$

where

$$
\tilde {M} _ {k} = (\Pi_ {n - 1} \dots \Pi_ {k + 1}) M _ {k} (\Pi_ {k + 1} \dots \Pi_ {n - 1}) \tag {3.4.5}
$$

for $k = 1 { : } n - 1$ . For example, in the $n = 4$ case we have

$$
\tilde {M} _ {3} \tilde {M} _ {2} \tilde {M} _ {1} P A = M _ {3} \cdot (\Pi_ {3} M _ {2} \Pi_ {3}) \cdot (\Pi_ {3} \Pi_ {2} M _ {1} \Pi_ {2} \Pi_ {3}) \cdot (\Pi_ {3} \Pi_ {2} \Pi_ {1}) A
$$

since the $\Pi _ { i }$ are symmetric. Moreover,

$$
\tilde {M} _ {k} = (\Pi_ {n - 1} \dots \Pi_ {k + 1}) \cdot (I _ {n} - \tau^ {(k)} e _ {k} ^ {T}) \cdot (\Pi_ {k + 1} \dots \Pi_ {n - 1}) = I _ {n} - \tilde {\tau} ^ {(k)} e _ {k} ^ {T}
$$

with $\tilde { \tau } ^ { ( k ) } = \Pi _ { n - 1 } \cdot \cdot \cdot \Pi _ { k + 1 } \tau ^ { ( k ) }$ . This shows that $\tilde { M _ { k } }$ is a Gauss transformation. The transformation from $\tau ^ { ( k ) }$ to $\tilde { \tau } ^ { ( k ) }$ is easy to implement in practice.

Algorithm 3.4.1 (Outer Product LU with Partial Pivoting) This algorithm computes the factorization $P A = L U$ where P is a permutation matrix encoded by $p i v ( 1 { : } n - 1 )$ , $L$ is unit lower triangular with $| \ell _ { i j } | \le 1$ , and U is upper triangular. For $i = 1 { : } n$ , $A ( i , i { : } n )$ is overwritten by $U ( i , i ; n )$ and $A ( i + 1 { : } n , i )$ is overwritten by $L ( i + 1 { : } n , i )$ . The permutation $P$ is given by $P = \Pi _ { n - 1 } \cdot \cdot \cdot \Pi _ { 1 }$ where $\Pi _ { k }$ is an interchange permutation obtained by swapping rows k and $p i v ( k )$ of $I _ { n }$ .

for $k = 1 { : } n - 1$

Determine $\mu$ with $k \leq \mu \leq n$ so |A(µ, k)| =  A(k:n, k) ∞

$$
p i v (k) = \mu
$$

$$
A (k,:) \leftrightarrow A (\mu ,:)
$$

$ { \mathbf { i } } \mathbf { f } \ A ( k , k ) \neq 0$

$$
\rho = k + 1: n
$$

$$
A (\rho , k) = A (\rho , k) / A (k, k)
$$

$$
A (\rho , \rho) = A (\rho , \rho) - A (\rho , k) A (k, \rho)
$$

end

end

The floating point overhead associated with partial pivoting is minimal from the standpoint of arithmetic as there are only $O ( n ^ { 2 } )$ comparisons associated with the search for the pivots. The overall algorithm involves $2 n ^ { 3 } / 3$ flops.

If Algorithm 3.4.1 is applied to

$$
A = \left[ \begin{array}{c c c} 3 & 1 7 & 1 0 \\ 2 & 4 & - 2 \\ 6 & 1 8 & - 1 2 \end{array} \right],
$$

then upon completion

$$
A = \left[ \begin{array}{c c c} 6 & 1 8 & - 1 2 \\ 1 / 2 & 8 & 1 6 \\ 1 / 3 & - 1 / 4 & 6 \end{array} \right]
$$

and $p i v = [ 3 , 3 ]$ . These two quantities encode all the information associated with the reduction:

$$
\left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{array} \right] \left[ \begin{array}{c c c} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{array} \right] A = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 1 / 2 & 1 & 0 \\ 1 / 3 & - 1 / 4 & 1 \end{array} \right] \left[ \begin{array}{c c c} 6 & 1 8 & - 1 2 \\ 0 & 8 & 1 6 \\ 0 & 0 & 6 \end{array} \right].
$$

To compute the solution to $A x \ = \ b$ after invoking Algorithm 3.4.1, we solve $L y = P b$ for y and $U x = y$ for x. Note that b can be overwritten by P b as follows

for $k = 1 { : } n - 1$

$$
b (k) \leftrightarrow b (p i v (k))
$$

end

We mention that if Algorithm 3.4.1 is applied to the problem,

$$
\left[ \begin{array}{c c}. 0 0 1 & 1. 0 0 \\ 1. 0 0 & 2. 0 0 \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right] = \left[ \begin{array}{c} 1. 0 0 \\ 3. 0 0 \end{array} \right],
$$

using 3-digit floating point arithmetic, then

$$
P   =   \left[ \begin{array}{c c} 0 & 1 \\ 1 & 0 \end{array} \right], \quad \hat {L}   =   \left[ \begin{array}{c c} 1. 0 0 & 0 \\ . 0 0 1 & 1. 0 0 \end{array} \right], \quad \hat {U}   =   \left[ \begin{array}{c c} 1. 0 0 & 2. 0 0 \\ 0 & 1. 0 0 \end{array} \right],
$$

and $\boldsymbol { \hat { x } } = [ 1 . 0 0 , . 9 9 6 ] ^ { T }$ . Recall from §3.3.2 that if Gaussian elimination without pivoting is applied to this problem, then the computed solution has O(1) error.

We mention that Algorithm 3.4.1 always runs to completion. If $A ( k { : } n , k ) = 0$ i n step k, then $M _ { k } = I _ { n }$ .

# 3.4.4 The Gaxpy Version

In §3.2 we developed outer product and gaxpy schemes for computing the LU factorization. Having just incorporated pivoting in the outer product version, it is equally straight forward to do the same with the gaxpy approach. Referring to Algorithm 3.2.2, we simply search the vector $\left| v ( j { : } n ) \right|$ in that algorithm for its maximal element and proceed accordingly.

Algorithm 3.4.2 (Gaxpy LU with Partial Pivoting) This algorithm computes the factorization $P A = L U$ where $P$ is a permutation matrix encoded by $p i v ( 1 : n - 1 )$ , $L$ is unit lower triangular with $| \ell _ { i j } | \le 1$ , and $U$ is upper triangular. For $i = 1 { : } n$ , $A ( i , i { : } n )$ is overwritten by $U ( i , i ; n )$ and $A ( i + 1 { : } n , i )$ is overwritten by $L ( i + 1 { : } n , i )$ . The permutation $P$ is given by $P = \Pi _ { n - 1 } \cdot \cdot \cdot \Pi _ { 1 }$ where $\Pi _ { k }$ is an interchange permutation obtained by swapping rows k and $p i v ( k )$ of $I _ { n }$ .

Initialize $L$ to the identity and $U$ to the zero matrix.

for j = 1:n
    if j = 1
    v = A(:, 1)
    else $\tilde{a} = \Pi_{j-1} \cdots \Pi_1 A(:, j)$ Solve $L(1:j-1, 1:j-1)z = \tilde{a}(1:j-1)$ for $z \in \mathbb{R}^{j-1}$ $U(1:j-1, j) = z$ , $v(j:n) = \tilde{a}(j:n) - L(j:n, 1:j-1) \cdot z$ end
    Determine $\mu$ with $j \leq \mu \leq n$ so $|v(\mu)| = \| v(j:n) \|_{\infty}$ and set $piv(j) = \mu$ $v(j) \leftrightarrow v(\mu)$ , $L(j, 1:j-1) \leftrightarrow L(\mu, 1:j-1)$ , $U(j, j) = v(j)$ if $v(j) \neq 0$ $L(j+1:n, j) = v(j+1:n)/v(j)$ end
end

As with Algorithm 3.4.1, this procedure requires $2 n ^ { 3 } / 3$ flops and $O ( n ^ { 2 } )$ comparisons.

# 3.4.5 Error Analysis and the Growth Factor

We now examine the stability that is obtained with partial pivoting. This requires an accounting of the rounding errors that are sustained during elimination and during the triangular system solving. Bearing in mind that there are no rounding errors associated with permutation, it is not hard to show using Theorem 3.3.2 that the computed solution $\hat { x }$ satisfies $( { \boldsymbol { A } } + { \boldsymbol { E } } ) { \hat { \boldsymbol { x } } } \ = \ b$ where

$$
| E | \leq n \mathbf {u} \left(2 | A | + 4 \hat {P} ^ {T} | \hat {L} | | \hat {U} |\right) + O (\mathbf {u} ^ {2}). \tag {3.4.6}
$$

Here we are assuming that $\hat { P } , \hat { L } .$ and $\hat { U }$ are the computed analogs of $P , L ,$ , and $U$ as produced by the above algorithms. Pivoting implies that the elements of $\hat { L }$ are bounded by one. Thus $\| { \hat { L } } \| _ { \infty } \leq n$ and we obtain the bound

$$
\| E \| _ {\infty} \leq n \mathbf {u} \left(2 \| A \| _ {\infty} + 4 n \| \hat {U} \| _ {\infty}\right) + O \left(\mathbf {u} ^ {2}\right). \tag {3.4.7}
$$

The problem now is to bound $\| \hat { U } \| _ { \infty }$ . Define the growth factor $\rho$ by

$$
\rho = \max _ {i, j, k} \frac {\left| \hat {a} _ {i j} ^ {(k)} \right|}{\| A \| _ {\infty}} \tag {3.4.8}
$$

where $\hat { A } ^ { ( k ) }$ is the computed version of the matrix $A ^ { ( k ) } = M _ { k } \Pi _ { k } \cdot \cdot \cdot M _ { 1 } \Pi _ { 1 } A$ . It follows that

$$
\| E \| _ {\infty} \leq 6 n ^ {3} \rho \| A \| _ {\infty} \mathbf {u} + O \left(\mathbf {u} ^ {2}\right). \tag {3.4.9}
$$

Whether or not this compares favorably with the ideal bound (2.7.20) hinges upon the size of the growth factor of $\rho .$ (The factor $n ^ { 3 }$ is not an operating factor in practice and may be ignored in this discussion.)

The growth factor measures how large the A-entries become during the process of elimination. Whether or not we regard Gaussian elimination with partial pivoting is safe to use depends upon what we can say about this quantity. From an average-case point of view, experiments by Trefethen and Schreiber (1990) suggest that $\rho$ is usually in the vicinity of $n ^ { 2 / 3 }$ . However, from the worst-case point of view, $\rho$ can be as large as $2 ^ { n - 1 }$ . In particular, if $A \in \mathbb { R } ^ { n \times n }$ is defined by

$$
a _ {i j} = \left\{ \begin{array}{r l} 1 & \text { if } i = j \text { or } j = n, \\ - 1 & \text { if } i > j, \\ 0 & \text { otherwise }, \end{array} \right.
$$

then there is no swapping of rows during Gaussian elimination with partial pivoting. We emerge with $A = L U$ and it can be shown that $u _ { n n } = 2 ^ { n - 1 }$ . For example,

$$
\left[ \begin{array}{r r r r} 1 & 0 & 0 & 1 \\ - 1 & 1 & 0 & 1 \\ - 1 & - 1 & 1 & 1 \\ - 1 & - 1 & - 1 & 1 \end{array} \right] = \left[ \begin{array}{r r r r} 1 & 0 & 0 & 0 \\ - 1 & 1 & 0 & 0 \\ - 1 & - 1 & 1 & 0 \\ - 1 & - 1 & - 1 & 1 \end{array} \right] \left[ \begin{array}{r r r r} 1 & 0 & 0 & 1 \\ 0 & 1 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 8 \end{array} \right].
$$

Understanding the behavior of $\rho$ requires an intuition about what makes the $U _ { - }$ factor large. Since $P A = L U$ implies $U = L ^ { - 1 } P A$ it would appear that the size of $L ^ { - 1 }$ is relevant. However, Stewart (1997) discusses why one can expect the L-factor to be well conditioned.

Although there is still more to understand about $\rho ,$ the consensus is that serious element growth in Gaussian elimination with partial pivoting is extremely rare. The method can be used with confidence.

# 3.4.6 Complete Pivoting

Another pivot strategy called complete pivoting has the property that the associated growth factor bound is considerably smaller than $2 ^ { n - 1 }$ . Recall that in partial pivoting, the kth pivot is determined by scanning the current subcolumn $A ( k { : } n , k )$ . In complete pivoting, the largest entry in the current submatrix $A ( k { : } n , k { : } n )$ is permuted into the $( k , k )$ position. Thus, we compute the upper triangularization

$$
M _ {n - 1} \Pi_ {n - 1} \dots M _ {1} \Pi_ {1} A \Gamma_ {1} \dots \Gamma_ {n - 1} = U.
$$

In step k we are confronted with the matrix

$$
A ^ {(k - 1)} = M _ {k - 1} \Pi_ {k - 1} \dots M _ {1} \Pi_ {1} A \Gamma_ {1} \dots \Gamma_ {k - 1}
$$

and determine interchange permutations $\Pi _ { k }$ and $\Gamma _ { k }$ such that

$$
\left| \left(\Pi_ {k} A ^ {(k - 1)} \Gamma_ {k}\right) _ {k k} \right| = \max _ {k \leq i, j \leq n} \left| \left(\Pi_ {k} A ^ {(k - 1)} \Gamma_ {k}\right) _ {i j} \right|.
$$

Algorithm 3.4.3 (Outer Product LU with Complete Pivoting) This algorithm computes the factorization $P A Q ^ { T } = L U$ where P is a permutation matrix encoded by $p i v ( 1 { : } n - 1 )$ , Q is a permutation matrix encoded by $c o l p i v ( 1 { : } n - 1 )$ , L is unit lower triangular with $| \ell _ { i j } | \le 1$ , and $U$ is upper triangular. For $i = 1 { : } n , A ( i , i { : } n )$ is overwritten by $U ( i , i ; n )$ and $A ( i + 1 { : } n , i )$ is overwritten by $L ( i + 1 { : } n , i )$ . The permutation $P$ is given by $P = \Pi _ { n - 1 } \cdot \cdot \cdot \Pi _ { 1 }$ where $\Pi _ { k }$ is an interchange permutation obtained by swapping rows k and rowpiv(k) of $I _ { n }$ . The permutation $Q$ is given by $Q = \Gamma _ { n - 1 } \cdot \cdot \cdot \Gamma _ { 1 }$ where $\Gamma _ { k }$ is an interchange permutation obtained by swapping rows k and colpiv(k) of $I _ { n }$ .

for k = 1:n - 1
    Determine $\mu$ with $k \leq \mu \leq n$ and $\lambda$ with $k \leq \lambda \leq n$ so $|A(\mu, \lambda)| = \max\{|A(i, j)| : i = k:n, j = k:n\}$ rowpiv(k) = $\mu$ $A(k, 1:n) \leftrightarrow A(\mu, 1:n)$ colpiv(k) = $\lambda$ $A(1:n, k) \leftrightarrow A(1:n, \lambda)$ if $A(k, k) \neq 0$ $\rho = k + 1:n$ $A(\rho, k) = A(\rho, k)/A(k, k)$ $A(\rho, \rho) = A(\rho, \rho) - A(\rho, k)A(k, \rho)$ end
end

This algorithm requires $2 n ^ { 3 } / 3$ flops and $O ( n ^ { 3 } )$ comparisons. Unlike partial pivoting, complete pivoting involves a significant floating point arithmetic overhead because of the two-dimensional search at each stage.

With the factorization $P A Q ^ { T } = L U$ in hand the solution to $A x = b$ proceeds as follows:

Step 1. Solve $L z = P b$ for z.

Step 2. Solve $U y = z \ { \mathrm { f o r } } \ y$ .

Step 3. Set $x = Q ^ { T } y .$ .

The rowpiv and colpiv representations can be used to form P b and $Q y$ , respectively.

Wilkinson (1961) has shown that in exact arithmetic the elements of the matrix $A ^ { ( k ) } = M _ { k } \Pi _ { k } \cdot \cdot \cdot M _ { 1 } \Pi _ { 1 } A \Gamma _ { 1 } \cdot \cdot \cdot \Gamma _ { k }$ satisfy

$$
\left| a _ {i j} ^ {(k)} \right| \leq k ^ {1 / 2} \left(2 \cdot 3 ^ {1 / 2} \dots k ^ {1 / k - 1}\right) ^ {1 / 2} \max \left| a _ {i j} \right|. \tag {3.4.10}
$$

The upper bound is a rather slow-growing function of k. This fact coupled with vast empirical evidence suggesting that $\rho$ is always modestly sized $( \mathrm { e . g } , \rho = 1 0 )$ permit us to conclude that Gaussian elimination with complete pivoting is stable. The method solves a nearby linear system $( { \boldsymbol { A } } + { \boldsymbol { E } } ) { \hat { \boldsymbol { x } } } = { \boldsymbol { b } }$ in the sense of (2.7.21). However, in general there is little reason to choose complete pivoting over partial pivoting. A possible exception is when A is rank deficient. In principal, complete pivoting can be used to reveal the rank of a matrix. Suppose rank $( A ) = r < n$ . It follows that at the beginning of step $r + 1$ , $A ( r + 1 { : } n , r + 1 { : } n ) = 0$ . This implies that $\Pi _ { k } = \Gamma _ { k } = M _ { k } = I$ for $k = r + 1$ :n and so the algorithm can be terminated after step r with the following factorization in hand:

$$
P A Q ^ {T}   =   L U   =   \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & I _ {n - r} \end{array} \right] \left[ \begin{array}{c c} U _ {1 1} & U _ {1 2} \\ 0 & 0 \end{array} \right]  .
$$

Here, $L _ { 1 1 }$ and $U _ { 1 1 }$ are $r { \mathrm { - } } \mathrm { b y } { \mathrm { - } } r$ and $L _ { 2 1 }$ and $U _ { 1 2 } ^ { T }$ are $( n - r ) – \mathrm { b y } – r$ . Thus, Gaussian elimination with complete pivoting can in principle be used to determine the rank of a matrix. Nevertheless, roundoff errors make the probability of encountering an exactly zero pivot remote. In practice one would have to “declare” A to have rank k if the pivot element in step $k + 1$ was sufficiently small. The numerical rank determination problem is discussed in detail in §5.5.

# 3.4.7 Rook Pivoting

A third type of LU stablization strategy called rook pivoting provides an interesting alternative to partial pivoting and complete pivoting. As with complete pivoting, it computes the factorization $P A Q \ = \ L U$ . However, instead of choosing as pivot the largest value in $\left| A ( k { : } n , k { : } n ) \right|$ , it searches for an element of that submatrix that is maximal in both its row and column. Thus, if

$$
A (k: n, k: n) = \left[ \begin{array}{c c c c} 2 4 & 3 6 & 1 3 & 6 1 \\ 4 2 & 6 7 & 7 2 & 5 0 \\ 3 8 & 1 1 & 3 6 & 4 3 \\ 5 2 & 3 7 & 4 8 & 1 6 \end{array} \right],
$$

then $^ { 6 6 } 7 2 ^ { 5 }$ would be identified by complete pivoting while $^ { 6 4 } { \bf 5 } 2 , ^ { 5 } \quad ^ { 6 } { \bf 7 } 2 , ^ { 5 }$ or $^ { 6 6 1 }$ would be acceptable with the rook pivoting strategy. To implement rook pivoting, the scanand-swap portion of Algorithm 3.4.3 is changed to

$$
\mu = k, \lambda = k, \tau = | a _ {\mu \lambda} |, s = 0
$$

$$
\text { while } \tau <   \| (A (k: n, \lambda) \| _ {\infty} \vee \tau <   \| (A (\mu , k: n) \| _ {\infty}
$$

$$
\text { if } \mod (s, 2) = 0
$$

$$
\text { Update } \mu \text { so   that } | a _ {\mu \lambda} | = \| (A (k: n, \lambda) \| _ {\infty} \text { with } k \leq \mu \leq n.
$$

$$
\text { Update } \lambda \text { so that } | a _ {\mu \lambda} | = \| (A (\mu , k: n) \| _ {\infty} \text { with } k \leq \lambda \leq n.
$$

end

$$
s = s + 1
$$

end

$$
\operatorname{rowpiv} (k) = \mu , A (k,:) \leftrightarrow A (\mu ,:) \operatorname{colpiv} (k) = \lambda , A (:, k) \leftrightarrow A (:, \lambda)
$$

The search for a larger $| a _ { \mu \lambda } |$ involves alternate scans of $A ( k { : } n , \lambda )$ and $A ( \mu , k { : } n )$ . The value of $\tau$ is monotone increasing and that ensures termination of the while-loop. In theory, the exit value of s could be $O ( n - k ) ^ { 2 } )$ , but in practice its value is $O ( 1 )$ . See Chang (2002). The bottom line is that rook pivoting represents the same $O ( n ^ { 2 } )$ overhead as partial pivoting, but that it induces the same level of reliability as complete pivoting.

# 3.4.8 A Note on Underdetermined Systems

If $A \in \mathbb { R } ^ { m \times n }$ with $m < n$ , rank $( A ) = m$ , and $b \in \mathbb { R } ^ { m }$ , then the linear system $A x = b$ is said to be underdetermined. Note that in this case there are an infinite number of solutions. With either complete or rook pivoting, it is possible to compute an LU factorization of the form

$$
P A Q ^ {T} = L \left[ U _ {1} \mid U _ {2} \right] \tag {3.4.11}
$$

where P and $Q$ are permutations, $\boldsymbol { L } \in \mathbb { R } ^ { m \times m }$ is unit lower triangular, and $U _ { 1 } \in \mathbb { R } ^ { m }$ ×m is nonsingular and upper triangular. Note that

$$
A x = b \Leftrightarrow (P A Q ^ {T}) (Q x) = (P b) \Leftrightarrow L \left[ U _ {1} \mid U _ {2} \right] \left[ \begin{array}{l} z _ {1} \\ z _ {2} \end{array} \right] = L (U _ {1} z _ {1} + U _ {2} z _ {2}) = c
$$

where $c = P b$ and

$$
\left[ \begin{array}{l} z _ {1} \\ z _ {2} \end{array} \right] = Q x.
$$

This suggests the following solution procedure:

Step 1. Solve $L y = P b$ for $\boldsymbol { y } \in \mathbb { R } ^ { m }$ .

Step 2. Choose $z _ { 2 } \in \mathbb { R } ^ { n - m }$ and solve $U _ { 1 } z _ { 1 } = y - U _ { 2 } z _ { 2 } { \mathrm { ~ f o r ~ } } z _ { 1 }$

Step 3. Set

$$
x = Q ^ {T} \left[ \begin{array}{l} z _ {1} \\ z _ {2} \end{array} \right].
$$

Setting $z _ { 2 } = 0$ is a natural choice. We have more to say about underdetermined systems in §5.6.2.

# 3.4.9 The LU Mentality

We offer three examples that illustrate how to think in terms of the LU factorization when confronted with a linear equation situation.

Example 1. Suppose A is nonsingular and n-by-n and that B is n-by-p. Consider the problem of finding $X \ ( n { \mathrm { - b y - } } p )$ so $A X = B$ . This is the multiple right hand side problem. If $X = \left[ \left. x _ { 1 } \right| \cdot \cdot \cdot \right| \left. x _ { p } \right]$ and $B = \left[ \left. b _ { 1 } \right| \cdot \cdot \cdot \right| \left. b _ { p } \right]$ are column partitions, then

Compute $P A = L U$

for k = 1:p

$\mathrm { S o l v e ~ } L y = P b _ { k } \mathrm { ~ a n d ~ t h e n ~ } U x _ { k } = y .$ (3.4.12)

end

If $B = I _ { n }$ , then we emerge with an approximation to $A ^ { - 1 }$ .

Example 2. Suppose we want to overwrite b with the solution to $A ^ { k } x = b$ where $A \in \mathbb { R } ^ { n \times n } , b \in \mathbb { R } ^ { n }$ , and k is a positive integer. One approach is to compute $C = A ^ { k }$ and then solve $C x = b$ . However, the matrix multiplications can be avoided altogether:

Compute P A = LU.

for j = 1:k

Overwrite b with the solution to Ly = P b. (3.4.13)

Overwrite b with the solution to U x = b.

end

As in Example 1, the idea is to get the LU factorization “outside the loop.”

Example 3. Suppose we are given $A \in \mathbb { R } ^ { n \times n } , \ d \in \mathbb { R } ^ { n }$ , and $c \in \mathbb { R } ^ { n }$ and that we want to compute $s = c ^ { T } A ^ { - 1 } d .$ . One approach is to compute $X = A ^ { - 1 }$ as discussed in (i) and then compute $s = c ^ { T } X d .$ However, it is more economical to proceed as follows:

Compute $P A = L U$ .

Solve Ly = P d and then $U x = y$ .

$$
s = c ^ {T} x
$$

An $^ { 6 6 } A ^ { - 1 9 }$ in a formula almost always means “solve a linear system” and almost never means “compute A−1.” $A ^ { - 1 }$

# 3.4.10 A Model Problem for Numerical Analysis

We are now in possession of a very important and well-understood algorithm (Gaussian elimination) for a very important and well-understood problem (linear equations). Let us take advantage of our position and formulate more abstractly what we mean by “problem sensitivity” and “algorithm stability.” Our discussion follows Higham (ASNA, §1.5–1.6), Stewart (MA, §4.3), and Trefethen and Bau (NLA, Lectures 12, 14, 15, and 22).

A problem is a function $f { : } D \to S$ from “data/input space” D to “solution/output space” S. A problem instance is f together with a particular $d \in D$ . We assume D and S are normed vector spaces. For linear systems, D is the set of matrix-vector pairs $( A , b )$ where $A \in \mathbb { R } ^ { n \times n }$ is nonsingular and $b \in \mathbb { R } ^ { n }$ . The function f maps $( A , b )$ to $A ^ { - 1 } b .$ , an element of S. For a particular A and $b , A x = b$ is a problem instance.

A perturbation theory for the problem f sheds light on the difference between $f ( d )$ and $f ( d + \Delta d )$ where $d \in D$ and $d + \Delta d \in D$ . For linear systems, we discussed in §2.6 the difference between the solution to $A x = b$ and the solution to $( A + \Delta A ) ( x + \Delta x ) =$ $( b + \Delta b )$ . We bounded $\| \Delta x \| / \| x \|$ in terms of $\| \Delta A \| / \| A \|$ and $\parallel \Delta b \parallel / \parallel b \parallel$ .

The conditioning of a problem refers to the behavior of f under perturbation at d. A condition number of a problem quantifies the rate of change of the solution with respect to the input data. If small changes in d induce relatively large changes in $f ( d )$ , then that problem instance is ill-conditioned. If small changes in d do not induce relatively large changes in $f ( d )$ , then that problem instance is well-conditioned. Definitions for “small” and “large” are required. For linear systems we showed in $\ S 2 . 6$ that the magnitude of the condition number $\kappa ( A ) = \| A \| \| A ^ { - 1 } \|$ determines whether an $A x = b$ problem is ill-conditioned or well-conditioned. One might say that a linear equation problem is well-conditioned if $\kappa ( A ) \approx O ( 1 )$ and ill-conditioned if $\kappa ( A ) \approx O ( 1 / \mathbf { u } )$ .

An algorithm for computing $f ( d )$ produces an approximation $\tilde { f } ( d )$ . Depending on the situation, it may be necessary to identify a particular software implementation of the underlying method. The $\tilde { f }$ function for Gaussian elimination with partial pivoting, Gaussian elimination with rook pivoting, and Gaussian elimination with complete pivoting are all different.

An algorithm for computing $f ( d )$ is stable if for some small $\Delta d ,$ the computed solution $\tilde { f } ( d )$ is close to $f ( d + \Delta d )$ . A stable algorithm nearly solves a nearby problem. $\mathrm { A n }$ algorithm for computing $f ( d )$ is backward stable if for some small $\Delta d ,$ the computed solution $\tilde { f } ( d )$ satisfies $\tilde { f } ( d ) = f ( d + \Delta d )$ . A backward stable algorithm exactly solves a nearby problem. Applied to a given linear system $A x = b $ , Gaussian elimination with complete pivoting is backward stable because the computed solution ˜x satisfies

$$
(A + \Delta) \tilde {x} = b
$$

and $\| \Delta \| / \| A \| \approx { \cal O } ( { \bf u } )$ . On the other hand, if b is specified by a matrix-vector product $b = M v$ , then

$$
(A + \Delta) \tilde {x} = M v + \delta
$$

where $\| \Delta \| / \| A \| \approx O ( \mathbf { u } )$ and $\delta / ( \lVert \mathbf { \nabla } M \rVert \lVert \mathbf { \nabla } v \rVert ) \approx { \cal O } ( \mathbf { u } )$ . Here, the underlying f is defined by $f \colon ( A , M , v ) \ \to \ A ^ { - 1 } ( M v )$ . In this case the algorithm is stable but not backward stable.

# Problems

P3.4.1 Let $A = L U$ be the LU factorization of n-by-n A with $| \ell _ { i j } | \le 1$ . Let $a _ { i } ^ { T }$ and $u _ { i } ^ { T }$ denote the ith rows of A and U , respectively. Verify the equation

$$
u _ {i} ^ {T} = a _ {i} ^ {T} - \sum_ {j = 1} ^ {i - 1} \ell_ {i j} u _ {j} ^ {T}
$$

and use it to show that $\parallel U \parallel _ { \infty } \leq 2 ^ { n - 1 } \parallel A \parallel _ { \infty }$ . (Hint: Take norms and use induction.)

P3.4.2 Show that if $P A Q = \mathrm { L U }$ is obtained via Gaussian elimination with complete pivoting, then no element of $U ( i , i ; n )$ is larger in absolute value than $| u _ { i i } |$ . Is this true with rook pivoting?

P3.4.3 Suppose $A \in \mathbb { R } ^ { n \times n }$ has an $_ { L U }$ factorization and that L and U are known. Give an algorithm which can compute the $( i , j )$ entry of $A ^ { - 1 }$ in approximately $( n - j ) ^ { 2 } + ( n - i ) ^ { 2 }$ flops.

P3.4.4 Suppose Xˆ is the computed inverse obtained via (3.4.12). Give an upper bound for $\Vert \ A \hat { X } - I \Vert _ { F } .$

P3.4.5 Extend Algorithm 3.4.3 so that it can produce the factorization (3.4.11). How many flops are required?

# Notes and References for 3.4

Papers concerned with element growth and pivoting include:

C.W. Cryer (1968). “Pivot Size in Gaussian Elimination,” Numer. Math. 12, 335–345.

J.K. Reid (1971). “A Note on the Stability of Gaussian Elimination,” J.Inst. Math. Applic. 8, 374–375.

P.A. Businger (1971). “Monitoring the Numerical Stability of Gaussian Elimination,” Numer. Math. 16, 360–361.

A.M. Cohen (1974). “A Note on Pivot Size in Gaussian Elimination,” Lin. Alg. Applic. 8, 361–68.

A.M. Erisman and J.K. Reid (1974). “Monitoring the Stability of the Triangular Factorization of a Sparse Matrix,” Numer. Math. 22, 183–186.

J. Day and B. Peterson (1988). “Growth in Gaussian Elimination,” Amer. Math. Monthly 95, 489–513.

N.J. Higham and D.J. Higham (1989). “Large Growth Factors in Gaussian Elimination with Pivoting,” SIAM J. Matrix Anal. Applic. 10, 155–164.

L.N. Trefethen and R.S. Schreiber (1990). “Average-Case Stability of Gaussian Elimination,” SIAM J. Matrix Anal. Applic. 11, 335–360.

N. Gould (1991). “On Growth in Gaussian Elimination with Complete Pivoting,” SIAM J. Matrix Anal. Applic. 12, 354–361.   
A. Edelman (1992). “The Complete Pivoting Conjecture for Gaussian Elimination is False,” Mathematica J. 2, 58–61.   
S.J. Wright (1993). “A Collection of Problems for Which Gaussian Elimination with Partial Pivoting is Unstable,” SIAM J. Sci. Stat. Comput. 14, 231–238.   
L.V. Foster (1994). “Gaussian Elimination with Partial Pivoting Can Fail in Practice,” SIAM J. Matrix Anal. Applic. 15, 1354–1362.   
A. Edelman and W. Mascarenhas (1995). “On the Complete Pivoting Conjecture for a Hadamard Matrix of Order 12,” Lin. Multilin. Alg. 38, 181–185.   
J.M. Pena (1996). “Pivoting Strategies Leading to Small Bounds of the Errors for Certain Linear Systems,” IMA J. Numer. Anal. 16, 141–153.   
J.L. Barlow and H. Zha (1998). “Growth in Gaussian Elimination, Orthogonal Matrices, and the 2-Norm,” SIAM J. Matrix Anal. Applic. 19, 807–815.   
P. Favati, M. Leoncini, and A. Martinez (2000). “On the Robustness of Gaussian Elimination with Partial Pivoting,” BIT 40, 62–73.   
As we mentioned, the size of L−1 is relevant to the growth factor. Thus, it is important to have an understanding of triangular matrix condition, see:   
D. Viswanath and L.N. Trefethen (1998). “Condition Numbers of Random Triangular Matrices,” SIAM J. Matrix Anal. Applic. 19, 564–581.   
The connection between small pivots and near singularity is reviewed in:   
T.F. Chan (1985). “On the Existence and Computation of LU Factorizations with Small Pivots,” Math. Comput. 42, 535–548.   
A pivot strategy that we did not discuss is pairwise pivoting. In this approach, 2-by-2 Gauss transformations are used to zero the lower triangular portion of A. The technique is appealing in certain multiprocessor environments because only adjacent rows are combined in each step, see:   
D. Sorensen (1985). “Analysis of Pairwise Pivoting in Gaussian Elimination,” IEEE Trans. Comput. C-34, 274–278.   
A related type of pivoting called tournament pivoting that is of interest in distributed memory computing is outlined in §3.6.3. For a discussion of rook pivoting and its properties, see:   
L.V. Foster (1997). “The Growth Factor and Efficiency of Gaussian Elimination with Rook Pivoting,” J. Comput. Appl. Math., 86, 177–194.   
G. Poole and L. Neal (2000). “The Rook’s Pivoting Strategy,” J. Comput. Appl. Math. 123, 353–369. X-W Chang (2002) “Some Features of Gaussian Elimination with Rook Pivoting,” BIT 42, 66–83.

# 3.5 Improving and Estimating Accuracy

Suppose we apply Gaussian elimination with partial pivoting to the n-by-n system Ax = b and that IEEE double precision arithmetic is used. Equation (3.4.9) essentially says that if the growth factor is modest then the computed solution ˆx satisfies

$$
(A + E) \hat {x} = b, \quad \| E \| _ {\infty} \approx \mathbf {u} \| A \| _ {\infty}. \tag {3.5.1}
$$

In this section we explore the practical ramifications of this result. We begin by stressing the distinction that should be made between residual size and accuracy. This is followed by a discussion of scaling, iterative improvement, and condition estimation. See Higham (ASNA) for a more detailed treatment of these topics.

We make two notational remarks at the outset. The infinity norm is used throughout since it is very handy in roundoff error analysis and in practical error estimation. Second, whenever we refer to “Gaussian elimination” in this section we really mean Gaussian elimination with some stabilizing pivot strategy such as partial pivoting.

# 3.5.1 Residual Size versus Accuracy

The residual of a computed solution ˆx to the linear system $A x \ = \ b$ is the vector $b - A \hat { x }$ . A small residual means that Axˆ effectively “predicts” the right hand side b. From Equation 3.5.1 we have $\parallel b - A { \hat { x } } \parallel _ { \infty } \approx \mathbf { u } \parallel A \parallel _ { \infty } \parallel { \hat { x } } \parallel _ { \infty }$ and so we obtain

Heuristic I. Gaussian elimination produces a solution xˆ with a relatively small residual.

Small residuals do not imply high accuracy. Combining Theorem 2.6.2 and (3.5.1), we see that

$$
\frac {\| \hat {x} - x \| _ {\infty}}{\| x \| _ {\infty}} \approx \mathbf {u} \kappa_ {\infty} (A). \tag {3.5.2}
$$

This justifies a second guiding principle.

Heuristic II. If the unit roundoff and condition satisfy $\mathbf { u } \approx 1 0 ^ { - d }$ and $\kappa _ { \infty } ( A ) \approx 1 0 ^ { q }$ , then Gaussian elimination produces a solution xˆ that has about $d - q$ correct decimal digits.

If u $\kappa _ { \infty } ( A )$ is large, then we say that A is ill-conditioned with respect to the machine precision.

As an illustration of the Heuristics I and II, consider the system

$$
{\left[ \begin{array}{l l}. 9 8 6 & . 5 7 9 \\ . 4 0 9 & . 2 3 7 \end{array} \right]} {\left[ \begin{array}{l} x _ {1} \\ x _ {2} \end{array} \right]} = {\left[ \begin{array}{l}. 2 3 5 \\ . 1 0 7 \end{array} \right]}
$$

in which $\kappa _ { \infty } ( A ) \approx 7 0 0$ and $x = [ 2 , - 3 ] ^ { T }$ . Here is what we find for various machine precisions:

<table><tr><td>u</td><td> $\hat{x}_{1}$ </td><td> $\hat{x}_{2}$ </td><td> $\frac{\|\hat{x}-x\|_{\infty}}{\|x\|_{\infty}}$ </td><td> $\frac{\|b-A\hat{x}\|_{\infty}}{\|A\|_{\infty}\|\hat{x}\|_{\infty}}$ </td></tr><tr><td> $10^{-3}$ </td><td>2.11</td><td>-3.17</td><td> $5\cdot10^{-2}$ </td><td> $2.0\cdot10^{-3}$ </td></tr><tr><td> $10^{-4}$ </td><td>1.986</td><td>-2.975</td><td> $8\cdot10^{-3}$ </td><td> $1.5\cdot10^{-4}$ </td></tr><tr><td> $10^{-5}$ </td><td>2.0019</td><td>-3.0032</td><td> $1\cdot10^{-3}$ </td><td> $2.1\cdot10^{-6}$ </td></tr><tr><td> $10^{-6}$ </td><td>2.00025</td><td>-3.00094</td><td> $3\cdot10^{-4}$ </td><td> $4.2\cdot10^{-7}$ </td></tr></table>

Whether or not to be content with the computed solution ˆx depends on the requirements of the underlying source problem. In many applications accuracy is not important but small residuals are. In such a situation, the ˆx produced by Gaussian elimination is probably adequate. On the other hand, if the number of correct digits in ˆx is an issue, then the situation is more complicated and the discussion in the remainder of this section is relevant.

# 3.5.2 Scaling

Let $\beta$ be the machine base $\left( \mathrm { t y p i c a l l y ~ } \beta = 2 \right)$ and define the diagonal matrices $D _ { 1 } =$ dia $\boldsymbol { \mathfrak { z } } ( \beta ^ { r _ { 1 } } , \ldots , \beta ^ { r _ { n } } )$ and $D _ { 2 } \ = \ \mathrm { d i a g } ( \beta ^ { c _ { 1 } } , \ldots , \beta ^ { c _ { n } } )$ . The solution to the n-by-n linear system $A x = b$ can be found by solving the scaled system $( D _ { 1 } ^ { - 1 } A D _ { 2 } ) y = D _ { 1 } ^ { - 1 } b$ using

Gaussian elimination and then setting $x = D _ { 2 } y$ . The scalings of A, $b ,$ and $y$ require only $O ( n ^ { 2 } )$ flops and may be accomplished without roundoff. Note that $D _ { 1 }$ scales equations and $D _ { 2 }$ scales unknowns.

It follows from Heuristic II that if ˆx and $\hat { y }$ are the computed versions of x and $y .$ , then

$$
\frac {\| D _ {2} ^ {- 1} (\hat {x} - x) \| _ {\infty}}{\| D _ {2} ^ {- 1} x \| _ {\infty}} = \frac {\| \hat {y} - y \| _ {\infty}}{\| y \| _ {\infty}} \approx \mathbf {u} \kappa_ {\infty} (D _ {1} ^ {- 1} A D _ {2}). \tag {3.5.3}
$$

Thus, if $\kappa _ { \infty } ( D _ { 1 } ^ { - 1 } A D _ { 2 } )$ can be made considerably smaller than $\kappa _ { \infty } ( A )$ , then we might expect a correspondingly more accurate ${ \hat { x } } .$ , provided errors are measured in the $" D _ { 2 } "$ norm defined by $\| z \| _ { D _ { 2 } } = \| D _ { 2 } ^ { - 1 } z \| _ { \infty }$ . This is the objective of scaling. Note that it encompasses two issues: the condition of the scaled problem and the appropriateness of appraising error in the $D _ { 2 } { \mathrm { - n o r m } }$ .

An interesting but very difficult mathematical problem concerns the exact minimization of $\kappa _ { p } ( D _ { 1 } ^ { - 1 } A D _ { 2 } )$ for general diagonal $D _ { i }$ and various $p .$ Such results as there are in this direction are not very practical. This is hardly discouraging, however, when we recall that (3.5.3) is a heuristic result, it makes little sense to minimize exactly a heuristic bound. What we seek is a fast, approximate method for improving the quality of the computed solution ${ \hat { x } } .$ .

One technique of this variety is simple row scaling. In this scheme $D _ { 2 }$ is the identity and $D _ { 1 }$ is chosen so that each row in $D _ { 1 } ^ { - 1 } A$ has approximately the same ∞- norm. Row scaling reduces the likelihood of adding a very small number to a very large number during elimination—an event that can greatly diminish accuracy.

Slightly more complicated than simple row scaling is row-column equilibration. Here, the object is to choose $D _ { 1 }$ and $D _ { 2 }$ so that the ∞-norm of each row and column of $D _ { 1 } ^ { - 1 } A D _ { 2 }$ belongs to the interval $[ 1 / \beta , 1 ]$ where $\beta$ is the base of the floating point system. For work along these lines, see McKeeman (1962).

It cannot be stressed too much that simple row scaling and row-column equilibration do not “solve” the scaling problem. Indeed, either technique can render a worse xˆ than if no scaling whatever is used. The ramifications of this point are thoroughly discussed in Forsythe and Moler (SLE, Chap. 11). The basic recommendation is that the scaling of equations and unknowns must proceed on a problem-by-problem basis. General scaling strategies are unreliable. It is best to scale (if at all) on the basis of what the source problem proclaims about the significance of each $a _ { i j }$ . Measurement units and data error may have to be considered.

# 3.5.3 Iterative Improvement

Suppose $A x = b$ has been solved via the partial pivoting factorization $P A = L U$ and that we wish to improve the accuracy of the computed solution ˆx. If we execute

$$
r = b - A \hat {x}
$$

$$
\text { Solve } L y = P r. \tag {3.5.4}
$$

$$
\text { Solve } U z = y.
$$

$$
x _ {\mathrm{new}} = \hat {x} + z
$$

then in exact arithmetic $A x _ { \mathrm { n e w } } = A { \hat { x } } + A z = ( b - r ) + r = b$ . Unfortunately, the naive floating point execution of these formulae renders an $x _ { \mathrm { n e w } }$ that is no more accurate than ˆx. This is to be expected since $\hat { r } = \mathsf { f l } ( b - A \hat { x } )$ has few, if any, correct significant digits. (Recall Heuristic I.) Consequently, $\hat { z } = \mathsf { f l } ( A ^ { - 1 } r ) \approx A ^ { - 1 }$ · noise ≈ noise is a very poor correction from the standpoint of improving the accuracy of xˆ. However, Skeel (1980) has an error analysis that indicates when (3.5.4) gives an improved $x _ { \mathrm { n e w } }$ from the standpoint of backward error. In particular, if the quantity

$$
\tau = \left(\| | A | | A ^ {- 1} | \| _ {\infty}\right) \left(\max _ {i} (| A | | x |) _ {i} / \min _ {i} (| A | | x |) _ {i}\right)
$$

is not too big, then (3.5.4) produces an $x _ { \mathrm { n e w } }$ such that $( A + E ) x _ { n e w } = b$ for very small E. Of course, if Gaussian elimination with partial pivoting is used, then the computed ˆx already solves a nearby system. However, this may not be the case for certain pivot strategies used to preserve sparsity. In this situation, the fixed precision iterative improvement step (3.5.4) can be worthwhile and cheap. See Arioli, Demmel, and Duff (1988).

In general, for (3.5.4) to produce a more accurate x, it is necessary to compute the residual $b - A \hat { x }$ with extended precision floating point arithmetic. Typically, this means that if t-digit arithmetic is used to compute $P A = L U , x$ , y, and z, then 2t-digit arithmetic is used to form $b - A \hat { x }$ . The process can be iterated. In particular, once we have computed $P A = L U$ and initialize $x = 0$ , we repeat the following:

$$
r = b - A x (\text { higher   precision })
$$

$$
\text { Solve } L y = P r \text { for } y \text { and } U z = y \text { for } z. \tag {3.5.5}
$$

$$
x = x + z
$$

We refer to this process as mixed-precision iterative improvement. The original A must be used in the high-precision computation of r. The basic result concerning the performance of (3.5.5) is summarized in the following heuristic:

Heuristic III. If the machine precision u and condition satisfy $\mathbf { u } = 1 0 ^ { - d }$ and $\kappa _ { \infty } ( A ) \approx$ 10q, then after k executions of (3.5.5), x has approximately min $\{ d , k ( d - q ) \}$ correct digits if the residual computation is performed with precision $\mathbf { u } ^ { 2 }$ .

Roughly speaking, if u $\kappa _ { \infty } ( A ) \leq 1$ , then iterative improvement can ultimately produce a solution that is correct to full (single) precision. Note that the process is relatively cheap. Each improvement costs $O ( n ^ { 2 } )$ , to be compared with the original $O ( n ^ { 3 } )$ investment in the factorization $P A = L U$ . Of course, no improvement may result if A is badly conditioned with respect to the machine precision.

# 3.5.4 Condition Estimation

Suppose that we have solved $A x = b$ via $P A = L U$ and that we now wish to ascertain the number of correct digits in the computed solution ˆx. It follows from Heuristic II that in order to do this we need an estimate of the condition $\kappa _ { \infty } ( A ) = \| A \| _ { \infty } \| A ^ { - 1 } \| _ { \infty }$ . Computing $\| A \| _ { \infty }$ poses no problem as we merely use the $O ( n ^ { 2 } )$ formula (2.3.10). The challenge is with respect to the factor $\| \ A ^ { - 1 } \ \| _ { \infty }$ . Conceivably, we could estimate this quantity by $\| { \hat { X } } \| _ { \infty }$ , where $\hat { X } ~ = ~ \left[ ~ \hat { x } _ { 1 } ~ | \cdots | ~ \hat { x } _ { n } ~ \right]$ and ${ \hat { x } } _ { i }$ is the computed solution to $A x _ { i } = e _ { i }$ . (See §3.4.9.) The trouble with this approach is its expense: ${ \hat { \kappa } } _ { \infty } = \| A \| _ { \infty } \| { \hat { X } } \| _ { \infty }$ costs about three times as much as ${ \hat { x } } .$ .

The central problem of condition estimation is how to estimate reliably the condition number in $O ( n ^ { 2 } )$ flops assuming the availability of $P A = L U$ or one of the factorizations that are presented in subsequent chapters. An approach described in Forsythe and Moler (SLE, p. 51) is based on iterative improvement and the heuristic

$$
\mathbf {u} \kappa_ {\infty} (A) \approx \| z \| _ {\infty} / \| x \| _ {\infty}
$$

where z is the first correction of x in (3.5.5).

Cline, Moler, Stewart, and Wilkinson (1979) propose an approach to the condition estimation problem thatis based on the implication

$$
A y = d \quad \Longrightarrow \quad \| A ^ {- 1} \| _ {\infty} \geq \| y \| _ {\infty} / \| d \| _ {\infty}.
$$

The idea behind their estimator is to choose d so that the solution y is large in norm and then set

$$
\hat {\kappa} _ {\infty} = \parallel A \parallel_ {\infty} \parallel y \parallel_ {\infty} / \parallel d \parallel_ {\infty}.
$$

The success of this method hinges on how close the ratio $\parallel y \parallel _ { \infty } / \parallel d \parallel _ { \infty }$ is to its maximum value $\| \ A ^ { - 1 } \ \| _ { \infty }$ .

Consider the case when $A = T$ is upper triangular. The relation between d and y is completely specified by the following column version of back substitution:

$$
\begin{array}{l} p (1: n) = 0 \\ \text { for } k = n: - 1: 1 \\ \end{array}
$$

Choose d(k).

$$
y (k) = (d (k) - p (k)) / T (k, k) \tag {3.5.6}
$$

$$
p (1: k - 1) = p (1: k - 1) + y (k) T (1: k - 1, k)
$$

end

Normally, we use this algorithm to solve a given triangular system $T y = d .$ However, in the condition estimation setting we are free to pick the right-hand side d subject to the “constraint” that y is large relative to d.

One way to encourage growth in y is to choose $d ( k )$ from the set $\{ - 1 , + 1 \}$ so as to maximize $y ( k )$ . If $p ( k ) \ge 0$ , then set $d ( k ) = - 1 . \mathrm { ~ I f ~ } p ( k ) < 0$ , then set $d ( k ) = + 1$ . In other words, $( 3 . 5 . 6 )$ is invoked with $d ( k ) = - \mathrm { s i g n } ( p ( k ) )$ . Overall, the vector d has the form $d ( 1 { : } n ) \stackrel { \cdot } { = } [ \pm 1 , \ldots , \pm 1 ] ^ { T }$ . Since this is a unit vector, we obtain the estimate $\hat { \kappa } _ { \infty } = \| \mathcal { T } \| _ { \infty } \| y \| _ { \infty }$ .

A more reliable estimator results if $d ( k ) \in \{ - 1 , + 1 \}$ is chosen so as to encourage growth both in $y ( k )$ and the running sum update $p ( 1 { : } k - 1 , k ) + T ( 1 { : } k - 1 , k ) y ( k )$ . In particular, at step k we compute

$$
\begin{array}{l} y (k) ^ {+} = (1 - p (k)) / T (k, k), \\ s (k) ^ {+} = | y (k) ^ {+} | + \| p (1: k - 1) + T (1: k - 1, k) y (k) ^ {+} \| _ {1}, \\ y (k) ^ {-} = (- 1 - p (k)) / T (k, k), \\ s (k) ^ {-} = | y (k) ^ {-} | + \| p (1: k - 1) + T (1: k - 1, k) y (k) ^ {-} \| _ {1}, \\ \end{array}
$$

and set

$$
y (k) = \left\{ \begin{array}{l l} y (k) ^ {+} & \text {if} s (k) ^ {+} \geq s (k) ^ {-}, \\ y (k) ^ {-} & \text {if} s (k) ^ {+} <   s (k) ^ {-}. \end{array} \right.
$$

This gives the following procedure.

Algorithm 3.5.1 (Condition Estimator) Let $T \in \mathbb { R } ^ { n \times n }$ be a nonsingular upper triangular matrix. This algorithm computes unit ∞-norm y and a scalar κ so $\| T y \| _ { \infty } \approx$ $\bar { 1 } / \| T ^ { - 1 } \| _ { \infty }$ and $\kappa \approx \kappa _ { \infty } ( T )$

$$
p (1: n) = 0
$$

for $k = n \colon - 1 \colon 1$

$$
y (k) ^ {+} = (1 - p (k)) / T (k, k)
$$

$$
y (k) ^ {-} = (- 1 - p (k)) / T (k, k)
$$

$$
p (k) ^ {+} = p (1: k - 1) + T (1: k - 1, k) y (k) ^ {+}
$$

$$
p (k) ^ {-} = p (1: k - 1) + T (1: k - 1, k) y (k) ^ {-}
$$

${ \bf i f } \ | y ( k ) ^ { + } | + \| \ p ( k ) ^ { + } \\\| _ { 1 } \geq | y ( k ) ^ { - } | + \| \ p ( k ) ^ { - } \| _ { 1 }$

$$
y (k) = y (k) ^ {+}
$$

$$
p (1: k - 1) = p (k) ^ {+}
$$

else

$$
y (k) = y (k) ^ {-}
$$

$$
p (1: k - 1) = p (k) ^ {-}
$$

end

end

$$
\kappa = \parallel y \parallel_ {\infty} \parallel T \parallel_ {\infty}
$$

$$
y = y / \parallel y \parallel_ {\infty}
$$

The algorithm involves several times the work of ordinary back substitution.

We are now in a position to describe a procedure for estimating the condition of a square nonsingular matrix A whose $P A = L U$ factorization is available:

Step 1. Apply the lower triangular version of Algorithm 3.5.1 to $U ^ { T }$ and obtain a large-norm solution to $U ^ { T } y = d .$

Step 2. Solve the triangular systems $L ^ { T } r = y , L w = P r$ , and $U z = w$

Step 3. Set $\begin{array} { r } { \hat { \kappa } _ { \infty } = \| \ A \| _ { \infty } \| z \| _ { \infty } / \| r \| _ { \infty } . } \end{array}$ .

Note that $\parallel z \parallel _ { \infty } \leq \parallel A ^ { - 1 } \parallel _ { \infty } \parallel r \parallel _ { \infty }$ . The method is based on several heuristics. First, if A is ill-conditioned and $P A = L U$ , then it is usually the case that U is correspondingly ill-conditioned. The lower triangle L tends to be fairly well-conditioned. Thus, it is more profitable to apply the condition estimator to U than to L. The vector r, because it solves $A ^ { T } P ^ { T } r = { \overset { \cdot } { d } } .$ tends to be rich in the direction of the left singular vector associated with $\sigma _ { \mathrm { m i n } } ( A )$ . Right-hand sides with this property render large solutions to the problem $A z = r$ .

In practice, it is found that the condition estimation technique that we have outlined produces adequate order-of-magnitude estimates of the true condition number.

# Problems

P3.5.1 Show by example that there may be more than one way to equilibrate a matrix.

P3.5.2 Suppose $P ( A + E ) = \hat { L } \hat { U }$ , where P is a permutation, Lˆ is lower triangular with $| \hat { \ell } _ { i j } | \le 1$ , and Uˆ is upper triangular. Show that $\hat { \kappa } _ { \infty } ( A ) \geq \| A \| _ { \infty } / ( \| E \| _ { \infty } + \mu )$ where $\mu = \operatorname* { m i n } \left| \hat { u } _ { i i } \right|$ . Conclude that if a small pivot is encountered when Gaussian elimination with pivoting is applied to A, then A is ill-conditioned. The converse is not true. (Hint: Let A be the matrix $B _ { n }$ defined in (2.6.9)).

P3.5.3 (Kahan (1966)) The system $A x = b$ where

$$
A = \left[ \begin{array}{c c c} 2 & - 1 & 1 \\ - 1 & 1 0 ^ {- 1 0} & 1 0 ^ {- 1 0} \\ 1 & 1 0 ^ {- 1 0} & 1 0 ^ {- 1 0} \end{array} \right], \qquad b = \left[ \begin{array}{c} 2 (1 + 1 0 ^ {- 1 0}) \\ - 1 0 ^ {- 1 0} \\ 1 0 ^ {- 1 0} \end{array} \right]
$$

has solution $x = [ 1 0 ^ { - 1 0 } ~ - 1 ~ 1 ] ^ { T }$ . (a) Show that if $( A + E ) y = b$ and $| E | \leq 1 0 ^ { - 8 } | A | .$ , then $| x - y | \leq$ $1 0 ^ { - 7 } | x |$ . That is, small relative changes in $A \mathrm { ^ { \circ } s }$ entries do not induce large changes in x even though $\kappa _ { \infty } ( \ r _ { A } ) = 1 0 ^ { 1 0 }$ . (b) Define $D = \mathrm { d i a g } ( 1 0 ^ { - 5 } , 1 0 ^ { 5 } , 1 0 ^ { 5 } )$ . Show that $\kappa _ { \infty } ( D A D ) \leq 5$ . (c) Explain what is going on using Theorem 2.6.3.

P3.5.4 Consider the matrix:

$$
T = \left[ \begin{array}{c c c c} 1 & 0 & M & - M \\ 0 & 1 & - M & M \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{array} \right] M \in \mathbb {R}.
$$

What estimate of $\kappa _ { \infty } ( T )$ is produced when (3.5.6) is applied with $d ( k ) = - \mathrm { s i g n } ( p ( k ) ) ?$ What estimate does Algorithm 3.5.1 produce? What is the true $\kappa _ { \infty } ( T ) ?$

P3.5.5 What does Algorithm 3.5.1 produce when applied to the matrix $B _ { n }$ given in (2.6.9)?

# Notes and References for 3.5

The following papers are concerned with the scaling of $A x = b$ problems:

F.L. Bauer (1963). “Optimally Scaled Matrices,” Numer. Math. 5, 73–87.

P.A. Businger (1968). “Matrices Which Can Be Optimally Scaled,” Numer. Math. 12, 346–48.

A. van der Sluis (1969). “Condition Numbers and Equilibration Matrices,” Numer. Math. 14, 14–23.

A. van der Sluis (1970). “Condition, Equilibration, and Pivoting in Linear Algebraic Systems,” Numer. Math. 15, 74–86.

C. McCarthy and G. Strang (1973). “Optimal Conditioning of Matrices,” SIAM J. Numer. Anal. 10, 370–388.

T. Fenner and G. Loizou (1974). “Some New Bounds on the Condition Numbers of Optimally Scaled Matrices,”J. ACM 21, 514–524.

G.H. Golub and J.M. Varah (1974). “On a Characterization of the Best L2-Scaling of a Matrix,” SIAM J. Numer. Anal. 11, 472–479.

R. Skeel (1979). “Scaling for Numerical Stability in Gaussian Elimination,” J. ACM 26, 494–526.

R. Skeel (1981). “Effect of Equilibration on Residual Size for Partial Pivoting,” SIAM J. Numer. Anal. 18, 449–55.

V. Balakrishnan and S. Boyd (1995). “Existence and Uniqueness of Optimal Matrix Scalings,” SIAM J. Matrix Anal. Applic. 16, 29–39.

Part of the difficulty in scaling concerns the selection of a norm in which to measure errors. An interesting discussion of this frequently overlooked point appears in:

W. Kahan (1966). “Numerical Linear Algebra,” Canadian Math. Bull. 9, 757–801.

For a rigorous analysis of iterative improvement and related matters, see:

C.B. Moler (1967). “Iterative Refinement in Floating Point,” J. ACM 14, 316-371.

M. Jankowski and M. Wozniakowski (1977). “Iterative Refinement Implies Numerical Stability,” BIT 17, 303–311.

R.D. Skeel (1980). “Iterative Refinement Implies Numerical Stability for Gaussian Elimination,” Math. Comput. 35, 817–832.

N.J. Higham (1997). “Iterative Refinement for Linear Systems and LAPACK,” IMA J. Numer. Anal. 17, 495–509.

A. Dax (2003). “A Modified Iterative Refinement Scheme,” SIAM J. Sci. Comput. 25, 1199–1213.   
J. Demmel, Y. Hida, W. Kahan, X.S. Li, S. Mukherjee, and E.J. Riedy (2006). “Error Bounds from Extra-Precise Iterative Refinement,” ACM Trans. Math. Softw. 32, 325–351.

The condition estimator that we described is given in:

A.K. Cline, C.B. Moler, G.W. Stewart, and J.H. Wilkinson (1979). “An Estimate for the Condition Number of a Matrix,” SIAM J. Numer. Anal. 16, 368-75.

Other references concerned with the condition estimation problem include:

C.G. Broyden (1973). “Some Condition Number Bounds for the Gaussian Elimination Process,” J. Inst. Math. Applic. 12, 273–286.

F. Lemeire (1973). “Bounds for Condition Numbers of Triangular Value of a Matrix,” Lin. Alg. Applic. 11, 1–2.

D.P. O’Leary (1980). “Estimating Matrix Condition Numbers,”SIAM J. Sci. Stat. Comput. 1, 205–209.

A.K. Cline, A.R. Conn, and C. Van Loan (1982). “Generalizing the LINPACK Condition Estimator,” in Numerical Analysis , J.P. Hennart (ed.), Lecture Notes in Mathematics No. 909, Springer-Verlag, New York.

A.K. Cline and R.K. Rew (1983). “A Set of Counter examples to Three Condition Number Estimators,” SIAM J. Sci. Stat. Comput. 4, 602–611.

W. Hager (1984). “Condition Estimates,” SIAM J. Sci. Stat. Comput. 5, 311–316.

N.J. Higham (1987). “A Survey of Condition Number Estimation for Triangular Matrices,” SIAM Review 29, 575–596.

N.J. Higham (1988). “FORTRAN Codes for Estimating the One-Norm of a Real or Complex Matrix with Applications to Condition Estimation (Algorithm 674),” ACM Trans. Math. Softw. 14, 381–396.

C.H. Bischof (1990). “Incremental Condition Estimation,” SIAM J. Matrix Anal. Applic. 11, 312– 322.

G. Auchmuty (1991). “A Posteriori Error Estimates for Linear Equations,” Numer. Math. 61, 1–6.

N.J. Higham (1993). “Optimization by Direct Search in Matrix Computations,” SIAM J. Matrix Anal. Applic. 14, 317–333.

D.J. Higham (1995). “Condition Numbers and Their Condition Numbers,” Lin. Alg. Applic. 214, 193–213.

G.W. Stewart (1997). “The Triangular Matrices of Gaussian Elimination and Related Decompositions,” IMA J. Numer. Anal. 17, 7–16.

# 3.6 Parallel LU

In §3.2.11 we show how to organize a block version of Gaussian elimination (without pivoting) so that the overwhelming majority of flops occur in the context of matrix multiplication. It is possible to incorporate partial pivoting and maintain the same level-3 fraction. After stepping through the derivation we proceed to show how the process can be effectively parallelized using the block-cyclic distribution ideas that were presented in §1.6.

# 3.6.1 Block LU with Pivoting

Throughout this section assume $A \in \mathbb { R } ^ { n \times n }$ and for clarity that $n = r N$

$$
A = \left[ \begin{array}{c c c} A _ {1 1} & \dots & A _ {1 N} \\ \vdots & \ddots & \vdots \\ A _ {N 1} & \dots & A _ {N N} \end{array} \right] \quad A _ {i, j} \in \mathbb {R} ^ {r \times r}. \tag {3.6.1}
$$

We revisit Algorithm 3.2.4 (nonrecursive block LU) and show how to incorporate partial pivoting.

The first step starts by applying scalar Gaussian elimination with partial pivoting to the first block column. Using an obvious rectangular matrix version of Algorithm 3.4.1 we obtain the following factorization:

$$
P _ {1} \left[ \begin{array}{c} A _ {1 1} \\ A _ {2 1} \\ \vdots \\ A _ {N 1} \end{array} \right] = \left[ \begin{array}{c} L _ {1 1} \\ L _ {2 1} \\ \vdots \\ L _ {N 1} \end{array} \right] U _ {1 1}. \tag {3.6.2}
$$

In this equation, $P _ { 1 } \in \mathbb { R } ^ { n \times n }$ is a permutation, $L _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ is unit lower triangular, and $U _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ is upper triangular.

The next task is to compute the first block row of U . To do this we set

$$
P _ {1} A = \left[ \begin{array}{c c c} \tilde {A} _ {1 1} & \dots & \tilde {A} _ {1 N} \\ \vdots & \ddots & \vdots \\ \tilde {A} _ {N 1} & \dots & \tilde {A} _ {N N} \end{array} \right], \quad \tilde {A} _ {i, j} \in \mathbb {R} ^ {r \times r}, \tag {3.6.3}
$$

and solve the lower triangular multiple-right-hand-side problem

$$
L _ {1 1} \left[ U _ {1 2} \mid \dots \mid U _ {1 N} \right] = \left[ \tilde {A} _ {1 2} \mid \dots \mid \tilde {A} _ {1 N} \right] \tag {3.6.4}
$$

for $U _ { 1 2 } , \dots , U _ { 1 N } \in \mathbb { R } ^ { r \times r }$ . At this stage it is easy to show that we have the partial factorization

$$
P _ {1} A = \left[ \begin{array}{c c c c} L _ {1 1} & 0 & \dots & 0 \\ \hline L _ {2 1} & I _ {r} & \dots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ L _ {N 1} & 0 & \dots & I _ {r} \end{array} \right] \left[ \begin{array}{c c} I _ {r} & 0 \\ \hline 0 & A ^ {(\mathrm{new})} \end{array} \right] \left[ \begin{array}{c c c c} U _ {1 1} & U _ {1 2} & \dots & U _ {1 N} \\ \hline 0 & I _ {r} & \dots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & I _ {r} \end{array} \right]
$$

where

$$
A ^ {(\text { new })} = \left[ \begin{array}{c c c} \tilde {A} _ {2 2} & \dots & \tilde {A} _ {2 N} \\ \vdots & \ddots & \vdots \\ \tilde {A} _ {N 2} & \dots & \tilde {A} _ {N N} \end{array} \right] - \left[ \begin{array}{c} L _ {2 1} \\ \vdots \\ L _ {N 1} \end{array} \right] [ U _ {1 2} | \dots | U _ {1 N} ]. \tag {3.6.5}
$$

Note that the computation of $A ^ { \mathrm { ( n e w ) } }$ is a level-3 operation as it involves one matrix multiplication per A-block.

The remaining task is to compute the pivoted LU factorization of $A ^ { ( \mathrm { n e w } ) }$ . Indeed, if

$$
P ^ {(\mathrm{new})} A ^ {(\mathrm{new})} = L ^ {(\mathrm{new})} U ^ {(\mathrm{new})}
$$

and

$$
P ^ {(\mathrm{new})} \left[ \begin{array}{c} L _ {2 1} \\ \vdots \\ L _ {N 1} \end{array} \right] = \left[ \begin{array}{c} \tilde {L} _ {2 1} \\ \vdots \\ \tilde {L} _ {N 1} \end{array} \right],
$$

then

$$
P A = \left[ \begin{array}{c c} L _ {1 1} & 0 \quad \dots \quad 0 \\ \hline \tilde {L} _ {2 1} & \\ \vdots & L ^ {(\text {new})} \\ \tilde {L} _ {N 1} & \end{array} \right] \left[ \begin{array}{c c} U _ {1 1} & U _ {1 2} \dots U _ {1 N} \\ \hline 0 & \\ \vdots & U ^ {(\text {new})} \\ 0 & \end{array} \right]
$$

is the pivoted block LU factorization of A with

$$
P = \left[ \begin{array}{c c} I _ {r} & 0 \\ 0 & P ^ {(\mathrm{new})} \end{array} \right] P _ {1}.
$$

In general, the processing of each block column in A is a four-part calculation:

Part A. Apply rectangular Gaussian Elimination with partial pivoting to a block column of A. This produces a permutation, a block column of $L _ { : }$ and a diagonal block of U . See (3.6.2).

Part B. Apply the Part A permutation to the “rest of $A . ^ { \mathfrak { n } }$ See (3.6.3).

Part C. Complete the computation of $U \mathrm { { ^ { * } s } }$ next block row by solving a lower triangular multiple right-hand-side problem. See (3.6.4).

Part D. Using the freshly computed L-blocks and U-blocks, update the “rest of $A . ^ { \mathfrak { n } }$ See (3.6.5).

The precise formulation of the method with overwriting is similar to Algorithm 3.2.4 and is left as an exercise.

# 3.6.2 Parallelizing the Pivoted Block LU Algorithm

Recall the discussion of the block-cyclic distribution in §1.6.2 where the parallel computation of the matrix multiplication update $C = C + A B$ was outlined. To provide insight into how the pivoted block LU algorithm can be parallelized, we examine a representative step in a small example that also makes use of the block-cyclic distribution.

Assume that $N = 8$ in (3.6.1) and that we have a $p _ { \mathrm { r o w } } { - } \mathrm { b y } { - } p _ { \mathrm { c o l } }$ processor network with $p _ { \mathrm { r o w } } = 2$ and $p _ { \mathrm { c o l } } ~ = ~ 2$ . At the start, the blocks of $A \ = \ ( A _ { i j } )$ are cyclically distributed as shown in Figure 3.6.1. Assume that we have carried out two steps of block $L U$ and that the computed $L _ { i j }$ and $U _ { i j }$ have overwritten the corresponding $A -$ blocks. Figure 3.6.2 displays the situation at the start of the third step. Blocks that are to participate in the Part A factorization

$$
P _ {3} \left[ \begin{array}{c} A _ {3 3} \\ \vdots \\ A _ {8 3} \end{array} \right] = \left[ \begin{array}{c} L _ {3 3} \\ \vdots \\ L _ {8 3} \end{array} \right] U _ {3 3}
$$

are highlighted. Typically, $p _ { \mathrm { r o w } }$ processors are involved and since the blocks are each $r { \mathrm { - } } \mathrm { b y } { \mathrm { - } } r$ , there are r steps as shown in (3.6.6).

<table><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>A11</td><td>A12</td><td>A13</td><td>A14</td><td>A15</td><td>A16</td><td>A17</td><td>A18</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>A21</td><td>A22</td><td>A23</td><td>A24</td><td>A25</td><td>A26</td><td>A27</td><td>A28</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>A31</td><td>A32</td><td>A33</td><td>A34</td><td>A35</td><td>A36</td><td>A37</td><td>A38</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>A41</td><td>A42</td><td>A43</td><td>A44</td><td>A45</td><td>A46</td><td>A47</td><td>A48</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>A51</td><td>A52</td><td>A53</td><td>A54</td><td>A55</td><td>A56</td><td>A57</td><td>A58</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>A61</td><td>A62</td><td>A63</td><td>A64</td><td>A65</td><td>A66</td><td>A67</td><td>A68</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>A71</td><td>A72</td><td>A73</td><td>A74</td><td>A75</td><td>A76</td><td>A77</td><td>A78</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>A81</td><td>A82</td><td>A83</td><td>A84</td><td>A85</td><td>A86</td><td>A87</td><td>A88</td></tr></table>

Figure 3.6.1.

Part A:

<table><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>U11</td><td>U12</td><td>U13</td><td>U14</td><td>U15</td><td>U16</td><td>U17</td><td>U18</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L21</td><td>U22</td><td>U23</td><td>U24</td><td>U25</td><td>U26</td><td>U27</td><td>U28</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L31</td><td>L32</td><td>A33</td><td>A34</td><td>A35</td><td>A36</td><td>A37</td><td>A38</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L41</td><td>L42</td><td>A43</td><td>A44</td><td>A45</td><td>A46</td><td>A47</td><td>A48</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L51</td><td>L52</td><td>A53</td><td>A54</td><td>A55</td><td>A56</td><td>A57</td><td>A58</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L61</td><td>L62</td><td>A63</td><td>A64</td><td>A65</td><td>A66</td><td>A67</td><td>A68</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L71</td><td>L72</td><td>A73</td><td>A74</td><td>A75</td><td>A76</td><td>A77</td><td>A78</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L81</td><td>L82</td><td>A83</td><td>A84</td><td>A85</td><td>A86</td><td>A87</td><td>A88</td></tr></table>

Figure 3.6.2.

for $j = 1 { : } r$

Columns $A _ { k k } ( : , j ) , \dotsc , A _ { N , k } ( : , j )$ are assembled in the processor housing $A _ { k k }$ , the “pivot processor”

The pivot processor determines the required row interchange and the Gauss transform vector

The swapping of the two A-rows may require the involvement of two processors in the network

The appropriate part of the Gauss vector together with (3.6.6) $A _ { k k } ( j , j { : } r )$ is sent by the pivot processor to the processors that house $A _ { k + 1 , k } , \dotsc , A _ { N , k }$

The processors that house $A _ { k k } , \ldots , A _ { N , k }$ carry out their share of the update, a local computation

# end

Upon completion, the parallel execution of Parts B and C follow. In the Part B computation, those blocks that may be involved in the row swapping have been highlighted. See Figure 3.6.3. This overhead generally engages the entire processor network, although communication is local to each processor column.

# Part B:

<table><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>U11</td><td>U12</td><td>U13</td><td>U14</td><td>U15</td><td>U16</td><td>U17</td><td>U18</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L21</td><td>U22</td><td>U23</td><td>U24</td><td>U25</td><td>U26</td><td>U27</td><td>U28</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L31</td><td>L32</td><td>U33</td><td>A34</td><td>A35</td><td>A36</td><td>A37</td><td>A38</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L41</td><td>L42</td><td>L43</td><td>A44</td><td>A45</td><td>A46</td><td>A47</td><td>A48</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L51</td><td>L52</td><td>L53</td><td>A54</td><td>A55</td><td>A56</td><td>A57</td><td>A58</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L61</td><td>L62</td><td>L63</td><td>A64</td><td>A65</td><td>A66</td><td>A67</td><td>A68</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L71</td><td>L72</td><td>L73</td><td>A74</td><td>A75</td><td>A76</td><td>A77</td><td>A78</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L81</td><td>L82</td><td>L83</td><td>A84</td><td>A85</td><td>A86</td><td>A87</td><td>A88</td></tr></table>

Figure 3.6.3.

Note that Part C involves just a single processor row while the “big” level-three update that follows typically involves the entire processor network. See Figures 3.6.4 and 3.6.5.

Part C: 

<table><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>U11</td><td>U12</td><td>U13</td><td>U14</td><td>U15</td><td>U16</td><td>U17</td><td>U18</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L21</td><td>U22</td><td>U23</td><td>U24</td><td>U25</td><td>U26</td><td>U27</td><td>U28</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L31</td><td>L32</td><td>U33</td><td>A34</td><td>A35</td><td>A36</td><td>A37</td><td>A38</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L41</td><td>L42</td><td>L43</td><td>A44</td><td>A45</td><td>A46</td><td>A47</td><td>A48</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L51</td><td>L52</td><td>L53</td><td>A54</td><td>A55</td><td>A56</td><td>A57</td><td>A58</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L61</td><td>L62</td><td>L63</td><td>A64</td><td>A65</td><td>A66</td><td>A67</td><td>A68</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L71</td><td>L72</td><td>L73</td><td>A74</td><td>A75</td><td>A76</td><td>A77</td><td>A78</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L81</td><td>L82</td><td>L83</td><td>A84</td><td>A85</td><td>A86</td><td>A87</td><td>A88</td></tr></table>

Figure 3.6.4.

Part D: 

<table><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>U11</td><td>U12</td><td>U13</td><td>U14</td><td>U15</td><td>U16</td><td>U17</td><td>U18</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L21</td><td>U22</td><td>U23</td><td>U24</td><td>U25</td><td>U26</td><td>U27</td><td>U28</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L31</td><td>L32</td><td>U33</td><td>A34</td><td>A35</td><td>A36</td><td>A37</td><td>A38</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L41</td><td>L42</td><td>L43</td><td>A44</td><td>A45</td><td>A46</td><td>A47</td><td>A48</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L51</td><td>L52</td><td>L53</td><td>A54</td><td>A55</td><td>A56</td><td>A57</td><td>A58</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L61</td><td>L62</td><td>L63</td><td>A64</td><td>A65</td><td>A66</td><td>A67</td><td>A68</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L71</td><td>L72</td><td>L73</td><td>A74</td><td>A75</td><td>A76</td><td>A77</td><td>A78</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L81</td><td>L82</td><td>L83</td><td>A84</td><td>A85</td><td>A86</td><td>A87</td><td>A88</td></tr></table>

Figure 3.6.5.

The communication overhead associated with Part D is masked by the matrix multiplications that are performed on each processor.

This completes the k = 3 step of parallel block LU with partial pivoting. The process can obviously be repeated on the trailing 5-by-5 block matrix. The virtues of the block-cyclic distribution are revealed through the schematics. In particular, the dominating level-3 step (Part D) is load balanced for all but the last few values of k. Subsets of the processor grid are used for the “smaller,” level-2 portions of the computation.

We shall not attempt to predict the fraction of time that is devoted to these computations or the propagation of the interchange permutations. Enlightenment in this direction requires benchmarking.

# 3.6.3 Tournament Pivoting

The decomposition via partial pivoting in Step A requires a lot of communication. An alternative that addresses this issue involves a strategy called tournament pivoting. Here is the main idea. Suppose we want to compute $P W = L U$ where the blocks of

$$
W = \left[ \begin{array}{l} W _ {1} \\ W _ {2} \\ W _ {3} \\ W _ {4} \end{array} \right] \in \mathbb {R} ^ {n \times r}
$$

are distributed around some network of processors. Assume that each $W _ { i }$ has many more rows than columns. The goal is to choose r rows from W that can serve as pivot rows. If we compute the “local” factorizations

$$
P _ {1} W _ {1} = L _ {1} U _ {1}, \qquad P _ {2} W _ {2} = L _ {2} U _ {2}, \qquad P _ {3} W _ {3} = L _ {3} U _ {3}, \qquad P _ {4} W _ {4} = L _ {4} U _ {4},
$$

via Gaussian elimination with partial pivoting, then the top r rows of the matrices $P _ { 1 } W _ { 1 } , P _ { 2 } W _ { 2 } , P _ { 3 } W _ { 3 }$ , are $P _ { 4 } W _ { 4 }$ are pivot row candidates. Call these square matrices $W _ { 1 } ^ { \prime } , W _ { 2 } ^ { \prime } , W _ { 3 } ^ { \prime } ,$ , and $W _ { 4 } ^ { \prime }$ and note that we have reduced the number of possible pivot rows from n to 4r.

Next we compute the factorizations

$$
P _ {1 2} W _ {1 2} ^ {\prime} = P _ {1 2} \left[ \begin{array}{c} W _ {1} ^ {\prime} \\ W _ {2} ^ {\prime} \end{array} \right] = L _ {1 2} U _ {1 2},
$$

$$
P _ {3 4} W _ {3 4} ^ {\prime} = P _ {3 4} \left[ \begin{array}{c} W _ {3} ^ {\prime} \\ W _ {4} ^ {\prime} \end{array} \right] = L _ {3 4} U _ {3 4},
$$

and recognize that the top r rows of $P _ { 1 2 } W _ { 1 2 } ^ { \prime }$ and the top r rows of $P _ { 3 4 } W _ { 3 4 } ^ { \prime }$ are even better pivot row candidates. Assemble these 2r rows into a matrix $W _ { 1 2 3 4 }$ and compute

$$
P _ {1 2 3 4} W _ {1 2 3 4} = L _ {1 2 3 4} U _ {1 2 3 4}.
$$

The top r rows of $P _ { 1 2 3 4 } W _ { 1 2 3 4 }$ are then the chosen pivot rows for the LU reduction of $W$ .

Of course, there are communication overheads associated with each round of the “tournament,” but the volume of interprocessor data transfers is much reduced. See Demmel, Grigori, and Xiang (2010).

# Problems

P3.6.1 In §3.6.1 we outlined a single step of block LU with partial pivoting. Specify a complete version of the algorithm.

P3.6.2 Regarding parallel block LU with partial pivoting, why is it better to “collect” all the permutations in Part A before applying them across the remaining block columns? In other words, why not propagate the Part A permutations as they are produced instead of having Part B, a separate permutation application step?

P3.6.3 Review the discussion about parallel shared memory computing in §1.6.5 and §1.6.6. Develop a shared memory version of Algorithm 3.2.1. Designate one processor for computation of the multipliers and a load-balanced scheme for the rank-1 update in which all the processors participate. A barrier is necessary because the rank-1 update cannot proceed until the multipliers are available. What if partial pivoting is incorporated?

# Notes and References for §3.6

See the scaLAPACK manual for a discussion of parallel Gaussian elimination as well as:

J. Ortega (1988). Introduction to Parallel and Vector Solution of Linear Systems, Plenum Press, New York.   
K. Gallivan, W. Jalby, U. Meier, and A.H. Sameh (1988). “Impact of Hierarchical Memory Systems on Linear Algebra Algorithm Design,” Int. J. Supercomput. Applic. 2, 12–48.   
J. Dongarra, I. Duff, D. Sorensen, and H. van der Vorst (1990). Solving Linear Systems on Vector and Shared Memory Computers, SIAM Publications, Philadelphia, PA.   
Y. Robert (1990). The Impact of Vector and Parallel Architectures on the Gaussian Elimination Algorithm, Halsted Press, New York.   
J. Choi, J.J. Dongarra, L.S. Osttrouchov, A.P. Petitet, D.W. Walker, and R.C. Whaley (1996). “Design and Implementation of the ScaLAPACK LU, QR, and Cholesky Factorization Routines,” Scientific Programming, 5, 173–184.   
X.S. Li (2005). “An Overview of SuperLU: Algorithms, Implementation, and User Interface,” ACM Trans. Math. Softw. 31, 302–325.   
S. Tomov, J. Dongarra, and M. Baboulin (2010). “Towards Dense Linear Algebra for Hybrid GPU Accelerated Manycore Systems,” Parallel Comput. 36, 232–240.

The tournament pivoting strategy is a central feature of the optimized LU implementation discussed in:

J. Demmel, L. Grigori, and H. Xiang (2011). “CALU: A Communication Optimal LU Factorization Algorithm,” SIAM J. Matrix Anal. Applic. 32, 1317-1350.

E. Solomonik and J. Demmel (2011). “Communication-Optimal Parallel 2.5D Matrix Multiplication and LU Factorization Algorithms,” Euro-Par 2011 Parallel Processing Lecture Notes in Computer Science, 2011, Volume 6853/2011, 90–109.

This page intentionally left blank
