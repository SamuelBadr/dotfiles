# 5.6 Square and Underdetermined Systems

The orthogonalization methods developed in this chapter can be applied to square systems and also to systems in which there are fewer equations than unknowns. In this brief section we examine the various possibilities.

# 5.6.1 Square Systems

The least squares solvers based on the QR factorization and the SVD can also be used to solve square linear systems. Figure 5.6.1 compares the associated flop counts. It is

<table><tr><td>Method</td><td>Flops</td></tr><tr><td>Gaussian elimination</td><td> $2n^{3}/3$ </td></tr><tr><td>Householder QR</td><td> $4n^{3}/3$ </td></tr><tr><td>Modified Gram-Schmidt</td><td> $2n^{3}$ </td></tr><tr><td>Singular value decomposition</td><td> $12n^{3}$ </td></tr></table>

Figure 5.6.1. Flops associated with various methods for square linear systems

assumed that the right-hand side is available at the time of factorization. Although Gaussian elimination involves the least amount of arithmetic, there are three reasons why an orthogonalization method might be considered:

• The flop counts tend to exaggerate the Gaussian elimination advantage. When memory traffic and vectorization overheads are considered, the QR approach is comparable in efficiency.

• The orthogonalization methods have guaranteed stability; there is no “growth factor” to worry about as in Gaussian elimination.

• In cases of ill-conditioning, the orthogonal methods give an added measure of reliability. QR with condition estimation is very dependable and, of course, SVD is unsurpassed when it comes to producing a meaningful solution to a nearly singular system.

We are not expressing a strong preference for orthogonalization methods but merely suggesting viable alternatives to Gaussian elimination.

We also mention that the SVD entry in the above table assumes the availability of b at the time of decomposition. Otherwise, $2 0 n ^ { 3 }$ flops are required because it then becomes necessary to accumulate the U matrix.

If the QR factorization is used to solve $A x = b$ , then we ordinarily have to carry out a back substitution: $R x = Q ^ { T } b$ . However, this can be avoided by “preprocessing” b. Suppose H is a Householder matrix such that $H b = \beta e _ { n }$ where $e _ { n }$ is the last column of $I _ { n }$ . If we compute the QR factorization of $( H A ) ^ { T }$ , then $A = H ^ { T } R ^ { T } Q ^ { T }$ and the system transforms to

$$
R ^ {T} y = \beta e _ {n}
$$

where $y = Q ^ { T } x$ . Since $R ^ { T }$ is lower triangular, $y = ( \beta / r _ { n n } ) e _ { n }$ and so

$$
x = \frac {\beta}{r _ {n n}} Q (:, n).
$$

# 5.6.2 Underdetermined Systems

In §3.4.8 we discussed how Gaussian elimination with either complete pivoting or rook pivoting can be used to solve a full-rank, underdetermined linear system

$$
A x = b, \quad A \in \mathbb {R} ^ {m \times n}, b \in \mathbb {R} ^ {m}. \tag {5.6.1}
$$

Various orthogonal factorizations can also be used to solve this problem. Notice that (5.6.1) either has no solution or has an infinity of solutions. In the second case, it is important to distinguish between algorithms that find the minimum 2-norm solution and those that do not. The first algorithm we present is in the latter category.

Assume that A has full row rank and that we apply QR with column pivoting to obtain

$$
Q ^ {T} A \Pi = \left[ R _ {1} \mid R _ {2} \right]
$$

where $R _ { 1 } \in \mathbb { R } ^ { m \times m }$ is upper triangular and $R _ { 2 } \in \mathbb { R } ^ { m \times ( n - m ) }$ . Thus, $A x = b$ transforms to

$$
(Q ^ {T} A \Pi) (\Pi^ {T} x) = [ R _ {1} \mid R _ {2} ] \left[ \begin{array}{c} z _ {1} \\ z _ {2} \end{array} \right] = Q ^ {T} b
$$

where

$$
\Pi^ {T} x = \left[ \begin{array}{l} z _ {1} \\ z _ {2} \end{array} \right]
$$

with $z _ { 1 } \in \mathbb { R } ^ { m }$ and $z _ { 2 } \in \mathbb { R } ^ { ( n - m ) }$ . By virtue of the column pivoting, $R _ { 1 }$ is nonsingular because we are assuming that A has full row rank. One solution to the problem is therefore obtained by setting $z _ { 1 } = R _ { 1 } ^ { - 1 } Q ^ { T } b$ and $z _ { 2 } = 0$ .

Algorithm 5.6.1 Given $A \in \mathbb { R } ^ { m \times n }$ with rank(A) = m and $b \in \mathbb { R } ^ { m }$ , the following algorithm finds an $\boldsymbol { x } \in \mathbb { R } ^ { n }$ such that $A x = b$ .

Compute QR-with-column-pivoting factorization: $Q ^ { T } A \Pi = R$

Solve $R ( 1 { : } m , 1 { : } m ) z _ { 1 } = Q ^ { T } b .$

$x = \Pi \left[ { \begin{array} { c } { z _ { 1 } } \\ { 0 } \end{array} } \right]$

This algorithm requires $2 m ^ { 2 } n - m ^ { 3 } / 3$ flops. The minimum norm solution is not guaranteed. (A different Π could render a smaller $z _ { 1 } . \ i )$ However, if we compute the QR factorization

$$
A ^ {T} = Q R = Q \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right]
$$

with $R _ { 1 } \in \mathbb { R } ^ { m \times m }$ , then $A x = b$ becomes

$$
(Q R) ^ {T} x   =   \left[ \begin{array}{c c} R _ {1} ^ {T} & 0 \end{array} \right] \left[ \begin{array}{c} z _ {1} \\ z _ {2} \end{array} \right]   =   b,
$$

where

$$
Q ^ {T} x = \left[ \begin{array}{l} z _ {1} \\ z _ {2} \end{array} \right], \qquad z _ {1} \in \mathbb {R} ^ {m},   z _ {2} \in \mathbb {R} ^ {n - m}.
$$

In this case the minimum norm solution does follow by setting $z _ { 2 } = 0$ .

Algorithm 5.6.2 Given $A \in \mathbb { R } ^ { m \times n }$ with rank $( A ) = m$ and $b \in \mathbb { R } ^ { m }$ , the following algorithm finds the minimum 2-norm solution to $A x = b$ .

Compute the QR factorization $A ^ { T } = Q R$ .

Solve $R ( 1 { : } m , 1 { : } m ) ^ { T } z = b$

Set $x = Q ( : , 1 { : } m ) z .$

This algorithm requires at most $2 m ^ { 2 } n - 2 m ^ { 3 } / 3$ flops.

The SVD can also be used to compute the minimum norm solution of an underdetermined $A x = b$ problem. If

$$
A = \sum_ {i = 1} ^ {r} \sigma_ {i} u _ {i} v _ {i} ^ {T}, \quad r = \operatorname{rank} (A)
$$

is the SVD of A, then

$$
x = \sum_ {i = 1} ^ {r} \frac {u _ {i} ^ {T} b}{\sigma_ {i}} v _ {i}.
$$

As in the least squares problem, the SVD approach is desirable if A is nearly rank deficient.

# 5.6.3 Perturbed Underdetermined Systems

We conclude this section with a perturbation result for full-rank underdetermined systems.

Theorem 5.6.1. Suppose rank $( A ) = m \leq n$ and that $A \in \mathbb { R } ^ { m \times n }$ , $\delta A \in \mathbb { R } ^ { m \times n }$ , $0 \neq$ $b \in \mathbb { R } ^ { m }$ , and $\delta b \in \mathbb { R } ^ { m }$ satisfy

$$
\epsilon = \max \{\epsilon_ {A}, \epsilon_ {b} \} <   \sigma_ {m} (A),
$$

where $\epsilon _ { A } = \parallel \delta A \parallel _ { 2 } / \parallel A \parallel _ { 2 }$ and $\epsilon _ { b } = \parallel \delta b \parallel _ { 2 } / \parallel b \parallel _ { 2 }$ . If x and $\hat { x }$ are minimum norm solutions that satisfy

$$
A x = b, \quad (A + \delta A) \hat {x} = b + \delta b,
$$

then

$$
\frac {\| \hat {x} - x \| _ {2}}{\| x \| _ {2}} \leq \kappa_ {2} (A) (\epsilon_ {A} \min \{2, n - m + 1 \} + \epsilon_ {b}) + O (\epsilon^ {2}).
$$

Proof. Let E and f be defined by $\delta A / \epsilon$ and $\delta b / \epsilon$ . Note that rank $( A + t E ) = m$ for all $0 < t < \epsilon$ and that

$$
x (t) = (A + t E) ^ {T} \left((A + t E) (A + t E) ^ {T}\right) ^ {- 1} (b + t f)
$$

satisfies $( A + t E ) x ( t ) = b + t f$ . By differentiating this expression with respect to t and setting t = 0 in the result we obtain

$$
\dot {x} (0) = \left(I - A ^ {T} (A A ^ {T}) ^ {- 1} A\right) E ^ {T} (A A ^ {T}) ^ {- 1} b + A ^ {T} (A A ^ {T}) ^ {- 1} (f - E x). \tag {5.6.2}
$$

Because

$$
\| x \| _ {2} = \| A ^ {T} (A A ^ {T}) ^ {- 1} b \| _ {2} \geq \sigma_ {m} (A) \| (A A ^ {T}) ^ {- 1} b \| _ {2},
$$

$$
\| I - A ^ {T} (A A ^ {T}) ^ {- 1} A \| _ {2} = \min (1, n - m),
$$

and

$$
\frac {\| f \| _ {2}}{\| x \| _ {2}} \leq \frac {\| f \| _ {2} \| A \| _ {2}}{\| b \| _ {2}},
$$

we have

$$
\begin{array}{l} \frac {\| \hat {x} - x \| _ {2}}{\| x \| _ {2}} = \frac {x (\epsilon) - x (0)}{\| x (0) \| _ {2}} = \epsilon \frac {\| \dot {x} (0) \| _ {2}}{\| x \| _ {2}} + O (\epsilon^ {2}) \\ \leq \epsilon \min (1, n - m) \left\{\frac {\| E \| _ {2}}{\| A \| _ {2}} + \frac {\| f \| _ {2}}{\| b \| _ {2}} + \frac {\| E \| _ {2}}{\| A \| _ {2}} \right\} \kappa_ {2} (A) + O (\epsilon^ {2}), \\ \end{array}
$$

from which the theorem follows.

Note that there is no $\kappa _ { 2 } { \left( A \right) } ^ { 2 }$ factor as in the case of overdetermined systems.

# Problems

P5.6.1 Derive equation (5.6.2).

P5.6.2 Find the minimal norm solution to the system Ax = b where $A = [ 1 2 3 ]$ and $b = 1$

P5.6.3 Show how triangular system solving can be avoided when using the QR factorization to solve an underdetermined system.

P5.6.4 Suppose $b , x \in \mathbb { R } ^ { n }$ are given and consider the following problems:

(a) Find an unsymmetric Toeplitz matrix T so $T x = b ,$ .   
(b) Find a symmetric Toeplitz matrix T so $T x = b ,$   
(c) Find a circulant matrix C so Cx = b.

Pose each problem in the form $A p = b$ where A is a matrix made up of entries from x and p is the vector of sought-after parameters.

# Notes and References for 5.6

For an analysis of linear equation solving via QR, see:

N.J. Higham (1991). “Iterative Refinement Enhances the Stability of QR Factorization Methods for Solving Linear Equations,” BIT 31, 447–468.

Interesting aspects concerning singular systems are discussed in:

T.F. Chan (1984). “Deflated Decomposition Solutions of Nearly Singular Systems,” SIAM J. Numer. Anal. 21, 738–754.

Papers concerned with underdetermined systems include:

R.E. Cline and R.J. Plemmons (1976). “L2-Solutions to Underdetermined Linear Systems,” SIAM Review 18, 92–106.

M.G. Cox (1981). “The Least Squares Solution of Overdetermined Linear Equations having Band or Augmented Band Structure,” IMA J. Numer. Anal. 1, 3–22.

M. Arioli and A. Laratta (1985). “Error Analysis of an Algorithm for Solving an Underdetermined System,” Numer. Math. 46, 255–268.

J.W. Demmel and N.J. Higham (1993). “Improved Error Bounds for Underdetermined System Solvers,” SIAM J. Matrix Anal. Applic. 14, 1–14.

S. Jokar and M.E. Pfetsch (2008). “Exact and Approximate Sparse Solutions of Underdetermined Linear Equations,” SIAM J. Sci. Comput. 31, 23–44.

The central matrix problem in the emerging field of compressed sensing is to solve an underdetermined system Ax = b such that the 1-norm of x is minimized, see:

E. Candes, J. Romberg, and T. Tao (2006). “Robust Uncertainty Principles: Exact Signal Reconstruction from Highly Incomplete Frequency Information,” IEEE Trans. Information Theory 52, 489–509.

D. Donoho (2006). “Compressed Sensing,” IEEE Trans. Information Theory 52, 1289–1306.

This strategy tends to produce a highly sparse solution vector x.
