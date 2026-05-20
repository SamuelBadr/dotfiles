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
