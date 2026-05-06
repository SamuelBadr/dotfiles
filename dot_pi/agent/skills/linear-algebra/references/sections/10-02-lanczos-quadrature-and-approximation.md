# 10.2 Lanczos, Quadrature, and Approximation

To deepen our understanding of the Lanczos process and to build an appreciation for its connections to other areas of applied mathematics, we consider an interesting approximation problem that has broad practical implications. Assume that A ∈ IRn×n $A \in \mathbb { R } ^ { n \times n }$ is a large, sparse, symmetric positive definite matrix whose eigenvalues reside in an interval $[ a , b ]$ . Let $f ( \lambda )$ be a given smooth function that is defined on $[ a , b ]$ . Given $u \in \mathbb { R } ^ { n }$ , our goal is to produce suitably tight lower and upper bounds b and B so that

$$
b \leq u ^ {T} \cdot f (A) \cdot u \leq B. \tag {10.2.1}
$$

In the approach we develop, the bounds are Gauss quadrature rule estimates of a certain integral and the evaluation of the rules requires the eigenvalues and eigenvectors of a Lanczos-produced tridiagonal matrix.

The $u ^ { T } f ( A ) u$ estimation problem has many applications throughout matrix computations. For example, suppose ˆx is an approximate solution to the symmetric positive definite system $A x = b$ and that we have computed the residual $r = b - A \hat { x }$ . Note that if $x _ { * } = A ^ { - 1 } b$ and $f ( \lambda ) = 1 / \lambda ^ { 2 }$ , then

$$
\| x _ {*} - \hat {x} \| _ {2} ^ {2} = (x _ {*} - \hat {x}) ^ {T} (x _ {*} - \hat {x}) = (A ^ {- 1} (b - A \hat {x})) ^ {T} (A ^ {- 1} (b - A \hat {x})) = r ^ {T} f (A) r.
$$

Thus, if we have a $u ^ { T } f ( A ) u$ estimation framework, then we can obtain $A x = b$ error bounds from residual bounds.

For an in-depth treatment of the material in this section, we refer the reader to the treatise by Golub and Meurant (2010). Our presentation is brief, informal, and stresses the linear algebra highlights.

# 10.2.1 Reformulation of the Problem

Without an integral in sight, it is mystifying as to why (10.2.1) involves quadrature at all. The key is to regard $u ^ { T } f ( A ) u$ as a Riemann-Stieltjes integral. In general, given a suitably nice integrand f(x) and weight function $w ( x )$ , the Riemann-Stieltjes integral

$$
I (f) = \int_ {a} ^ {b} f (x) d w (x)
$$

is a limit of sums of the form

$$
S _ {N} = \sum_ {\mu = 1} ^ {N} f (c _ {\mu}) (w (x _ {\mu}) - w (x _ {\mu + 1}))
$$

where $a = x _ { N } < \cdot \cdot \cdot < x _ { 1 } = b$ and $x _ { \mu + 1 } \leq c _ { \mu } \leq x _ { \mu }$ . Note that if w is piecewise constant on $[ a , b ]$ , then the only nonzero terms in $S _ { N }$ arise from subintervals that house a “w-jump.” For example, suppose $a = \lambda _ { n } < \lambda _ { 2 } < \cdots < \lambda _ { 1 } = b$ and that

$$
w (\lambda) = \left\{ \begin{array}{l l} w _ {n + 1} & \text { if } \lambda <   a, \\ w _ {\mu} & \text { if } \lambda_ {\mu} \leq \lambda <   \lambda_ {\mu - 1}, \\ w _ {1} & \text { if } b \leq \lambda , \end{array} \right. \quad \mu = 2: n, \tag {10.2.2}
$$

where $0 \leq w _ { n + 1 } \leq \cdot \cdot \cdot \leq w _ { 1 }$ . By considering the behavior of $S _ { N }$ as $N  \infty$ , we see that

$$
\int_ {a} ^ {b} f (\lambda) d w (\lambda) = \sum_ {\mu = 1} ^ {n} (w _ {\mu} - w _ {\mu + 1}) \cdot f (\lambda_ {\mu}). \tag {10.2.3}
$$

We are now set to explain why $u ^ { T } f ( A ) { \mathrm { : } }$ u is “secretly” a Riemann-Stieltjes integral. Let

$$
A = X \Lambda X ^ {T}, \quad \Lambda = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right), \tag {10.2.4}
$$

be a Schur decomposition of A with $\lambda _ { n } \leq \cdots \leq \lambda _ { 1 }$ . It follows that

$$
u ^ {T} f (A) u = (X ^ {T} u) ^ {T} \cdot f (\Lambda) \cdot (X ^ {T} u) = \sum_ {\mu = 1} ^ {n} [ X ^ {T} u ] _ {\mu} ^ {2} \cdot f (\lambda_ {\mu}).
$$

If we set

$$
w _ {\mu} = \left[ X ^ {T} u \right] _ {\mu} ^ {2} + \dots + \left[ X ^ {T} u \right] _ {n} ^ {2}, \quad \mu = 1: n + 1, \tag {10.2.5}
$$

in (10.2.2), then (10.2.3) becomes

$$
\int_ {a} ^ {b} f (\lambda) d w (\lambda) = \sum_ {\mu = 1} ^ {n} \left[ X ^ {T} u \right] _ {\mu} ^ {2} \cdot f \left(\lambda_ {\mu}\right) = u ^ {T} f (A) u. \tag {10.2.6}
$$

Our plan is to approximate this integral using Gauss quadrature.

# 10.2.2 Some Gauss-Type Quadrature Rules and Bounds

Given an accuracy-related parameter k, an interval $[ a , b ]$ , and a weight function $w ( \lambda )$ , a Gauss-type quadrature rule for the integral

$$
I (f) = \int_ {a} ^ {b} f (\lambda) d w (\lambda)
$$

involves a carefully constructed linear combination of f-evaluations across $[ a , b ]$ . The evaluation points (called nodes) and the coefficients (called weights) that define the linear combination are determined to make the rule correct for polynomials up to a certain degree that is related to k. Here are four examples:

1. Gauss. Compute weights $w _ { 1 } , \ldots , w _ { k }$ and nodes $t _ { 1 } , \ldots , t _ { k }$ so if

$$
I _ {G} (f) = \sum_ {i = 1} ^ {k} w _ {i} f (t _ {i}) \tag {10.2.7}
$$

then $I ( f ) = I _ { G } ( f )$ for all polynomials f that have degree $2 k - 1$ or less.

2. Gauss-Radau(a). Compute weights $w _ { a } , w _ { 1 } , \ldots , w _ { k }$ and nodes $t _ { 1 } , \ldots , t _ { k }$ so if

$$
I _ {G R (a)} (f) = w _ {a} f (a) + \sum_ {i = 1} ^ {k} w _ {i} f (t _ {i}) \tag {10.2.8}
$$

then $I ( f ) = I _ { G R ( a ) } ( f )$ for all polynomials f that have degree 2k or less.

3. Gauss-Radau(b). Compute weights $w _ { b } , w _ { 1 } , \dotsc , w _ { k }$ and nodes $t _ { 1 } , \ldots , t _ { k }$ so if

$$
I _ {G R (b)} (f) = w _ {b} f (b) + \sum_ {i = 1} ^ {k} w _ {i} f (t _ {i}) \tag {10.2.9}
$$

then $I ( f ) = I _ { G R ( b ) } ( f )$ for all polynomials f that have degree 2k or less.

4. Gauss-Lobatto. Compute weights $w _ { a } , w _ { b } , w _ { 1 } , \ldots , w _ { k }$ and nodes $t _ { 1 } , \ldots , t _ { k }$ so if

$$
I _ {G L} (f) = w _ {a} f (a) + w _ {b} f (b) + \sum_ {i = 1} ^ {k} w _ {i} f (t _ {i}) \tag {10.2.10}
$$

then $I ( f ) = I _ { G L } ( f )$ for all polynomials f that have degree $2 k + 1$ or less.

Each of these rules has a neatly specified error. It can be shown that

$$
\int_ {a} ^ {b} f (\lambda) d w (\lambda) = \left\{ \begin{array}{l l} I _ {G} (f) & + R _ {G} (f), \\ I _ {G R (a)} (f) & + R _ {G R (a)} (f), \\ I _ {G R (b)} (f) & + R _ {G R (b)} (f), \\ I _ {G L} (f) & + R _ {G G} (f), \end{array} \right.
$$

where

$$
R _ {G} (f) \quad = \frac {f ^ {(2 k)} (\eta)}{(2 n) !} \int_ {a} ^ {b} \left[ \prod_ {i = 1} ^ {k} (\lambda - t _ {i}) \right] ^ {2} d w (\lambda), \quad a <   \eta <   b,
$$

$$
R _ {G R (a)} (f) = \frac {f ^ {(2 k + 1)} (\eta)}{(2 k + 1) !} \int_ {a} ^ {b} (\lambda - a) \left[ \prod_ {i = 1} ^ {k} (\lambda - t _ {i}) \right] ^ {2} d w (\lambda), \quad a <   \eta <   b,
$$

$$
R _ {G R (b)} (f) = \frac {f ^ {(2 k + 1)} (\eta)}{(2 k + 1) !} \int_ {a} ^ {b} (\lambda - b) \left[ \prod_ {i = 1} ^ {k} (\lambda - t _ {i}) \right] ^ {2} d w (\lambda), \quad a <   \eta <   b,
$$

$$
R _ {G L} (f) = \frac {f ^ {(2 k + 2)} (\eta)}{(2 k + 2) !} \int_ {a} ^ {b} (\lambda - a) (\lambda - b) \left[ \prod_ {i = 1} ^ {k} (\lambda - t _ {i}) \right] ^ {2} d w (\lambda), \quad a <   \eta <   b.
$$

If the derivative in the remainder term does not change sign across [a, b], then the rule can be used to produce a bound. For example, if $f ( \lambda ) = 1 / \lambda ^ { 2 }$ and $0 < a < b$ , then $f ^ { ( 2 k ) }$ is positive, $f ^ { ( 2 k + 1 ) }$ is negative, and we have

$$
I _ {G} (f) \leq \int_ {a} ^ {b} f (\lambda) d w (\lambda) \leq I _ {G R (a)} (f).
$$

With this strategy, we can produce lower and upper bounds by selecting and evaluating the right rule. For this to be practical, the behavior of f ’s higher derivatives must be known and the required rules must be computable.

# 10.2.3 The Tridiagonal Connection

It turns out that the evaluation of a given Gauss quadrature rule involves a tridiagonal matrix and its eigenvalues and eigenvectors. To develop a strategy that is based upon this connection, we need three facts about orthogonal polynomials and Gauss quadrature.

Fact 1. Given $[ a , b ]$ and $w ( \lambda )$ , there is a sequence of polynomials $p _ { 0 } ( \lambda ) , p _ { 1 } ( \lambda ) , . . .$ . that satisfy

$$
\int_ {a} ^ {b} p _ {i} (\lambda) \cdot p _ {j} (\lambda) \cdot d w (\lambda) = \left\{ \begin{array}{l l} 1 & \mathrm{if} i = j, \\ 0 & \mathrm{if} i \neq j, \end{array} \right.
$$

with the property that the degree of $p _ { k } ( \cdot )$ is k for $k \geq 0$ . The polynomials are unique up to a factor of ±1 and they satisfy a 3-term recurrence

$$
\gamma_ {k} p _ {k} (\lambda) = (\lambda - w _ {k}) p _ {k - 1} (\lambda) - \gamma_ {k - 1} p _ {k - 2} (\lambda)
$$

where $p _ { - 1 } ( \lambda ) \equiv 0$ and $p _ { 0 } ( \lambda ) \equiv 1$ .

Fact 2. The zeros of $p _ { k } ( \lambda )$ are the eigenvalues of the tridiagonal matrix

$$
T _ {k} = \left[ \begin{array}{c c c c c} \omega_ {1} & \gamma_ {1} & 0 & \dots & 0 \\ \gamma_ {1} & \omega_ {2} & \ddots & & \vdots \\ 0 & \ddots & \ddots & \ddots & 0 \\ \vdots & & \ddots & \omega_ {k - 1} & \gamma_ {k - 1} \\ 0 & \dots & 0 & \gamma_ {k - 1} & \omega_ {k} \end{array} \right].
$$

Since the $\gamma _ { i }$ are nonzero, it follows from Theorem 8.4.1 that the eigenvalues are distinct.

Fact 3. If

$$
S ^ {T} T _ {k} S = \mathrm{diag} (\theta_ {1}, \dots , \theta_ {k}) \tag {10.2.11}
$$

is a Schur decomposition of $T _ { k } .$ , then the nodes and weights for the Gauss rule (10.2.7) are given by $t _ { i } = \theta _ { i }$ and $w _ { i } = s _ { 1 i } ^ { 2 }$ for i = 1:k. In other words,

$$
I _ {G} (f) = \sum_ {i = 1} ^ {k} s _ {1 i} ^ {2} \cdot f (\theta_ {i}). \tag {10.2.12}
$$

Thus, the only remaining issue is how to construct $T _ { k }$ so that it defines a Gauss rule for (10.2.6).

# 10.2.4 Gauss Quadrature via Lanczos

We show that if we apply the symmetric Lanczos process (Algorithm 10.1.1) with starting vector $q _ { 1 } = u / \parallel u \parallel _ { 2 }$ , then the tridiagonal matrices that the method generates are exactly what we need to compute $I _ { G } ( f )$ .

We first link the Lanczos process to a sequence of orthogonal polynomials. Recall from §10.1.1 that the kth Lanczos vector $q _ { k }$ is in the Krylov subspace $\kappa ( A , q _ { 1 } , k )$ . It follows that $q _ { k } ~ = ~ p _ { k } ( A ) q _ { 1 }$ for some degree-k polynomial. From Algorithm 10.1.1 we know that

$$
\beta_ {k} q _ {k + 1} = (A - \alpha_ {k} I) q _ {k} - \beta_ {k - 1} q _ {k - 1}
$$

where $\beta _ { 0 } q _ { 0 } \equiv 0$ and so

$$
\beta_ {k} p _ {k + 1} (A) q _ {1} = (A - \alpha_ {k} I) p _ {k} (A) q _ {1} - \beta_ {k - 1} p _ {k - 1} (A) q _ {1}.
$$

From this we conclude that the polynomials satisfy a 3-term recurrence:

$$
\beta_ {k} p _ {k + 1} (\lambda) = (\lambda - \alpha_ {k}) p _ {k} (\lambda) - \beta_ {k - 1} ^ {2} p _ {k - 1} (\lambda). \tag {10.2.13}
$$

These polynomials are orthogonal with respect to the $u ^ { T } f ( A ) u$ weight function defined in (10.2.5). To see this, note that

$$
\begin{array}{l} \int_ {a} ^ {b} p _ {i} (\lambda) p _ {j} (\lambda) d w (\lambda) = \sum_ {\mu = 1} ^ {n} [ X ^ {T} u ] _ {\mu} ^ {2} \cdot p _ {i} (\lambda_ {\mu}) \cdot p _ {j} (\lambda_ {\mu}) \\ = (X ^ {T} u) ^ {T} (p _ {i} (\Lambda) \cdot p _ {j} (\Lambda)) \cdot (X ^ {T} u) \\ = u ^ {T} \left(X \cdot p _ {i} (\Lambda) \cdot X ^ {T}\right) \left(X \cdot p _ {j} (\Lambda) \cdot X ^ {T}\right) u \\ = u ^ {T} \left(p _ {i} (A) p _ {j} (A)\right) u \\ = (p _ {i} (A) u) ^ {T} (p _ {j} (A) u) = \| u \| _ {2} ^ {2} q _ {i} ^ {T} q _ {j} = 0. \\ \end{array}
$$

Coupled with (10.2.13) and Facts 1-3, this result tells us that we can generate an approximation $\sigma = I _ { G } ( f )$ to $u ^ { T } f ( A ) u$ as follows:

Step 1: With starting vector $q _ { 1 } = u / \parallel u \parallel _ { 2 }$ , use the Lanczos process to compute the partial tridiagonalization $A Q _ { k } = Q _ { k } T _ { k } + r _ { k } e _ { k } ^ { T }$ . (See (10.1.4).)

Step 2: Compute the Schur decomposition $S ^ { T } T _ { k } S = \operatorname { d i a g } ( \theta _ { 1 } , . . . , \theta _ { k } )$

Step 3: Set $\sigma = s _ { 1 1 } ^ { 2 } f ( \theta _ { 1 } ) + \cdot \cdot \cdot + s _ { 1 k } ^ { 2 } f ( \theta _ { k } )$ .

See Golub and Welsch (1969) for a more rigorous derivation of this procedure.

# 10.2.5 Computing the Gauss-Radau Rule

Recall from (10.2.1) that we are interested in upper and lower bounds. In light of our remarks at the end of §10.2.2, we need techniques for evaluating other Gauss quadrature rules. By way of illustration, we show how to compute $I _ { G R ( a ) }$ defined in (10.2.8). Guided by Gauss quadrature theory, we run the Lanczos process for k steps as if we were setting out to compute $I _ { G } ( f )$ . We then must determine $\tilde { \alpha } _ { k + 1 }$ so that if

$$
\tilde {T} _ {k + 1} = \left[ \begin{array}{c c c c c c c} \alpha_ {1} & \beta_ {1} & 0 & \dots & 0 & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots & \vdots \\ 0 & \ddots & \ddots & \ddots & \vdots & \\ \vdots & & \ddots & \alpha_ {k - 1} & \beta_ {k - 1} & 0 \\ 0 & \dots & & \beta_ {k - 1} & \alpha_ {k} & \beta_ {k} \\ \hline 0 & \dots & \dots & 0 & \beta_ {k} & \tilde {\alpha} _ {k + 1} \end{array} \right]
$$

then $a \in \lambda ( \tilde { T } _ { k + 1 } )$ . By considering the top and bottom halves of the equation

$$
\tilde {T} _ {k + 1} \left[ \begin{array}{c} x \\ - 1 \end{array} \right] = a \left[ \begin{array}{c} x \\ - 1 \end{array} \right], \qquad x \in \mathbb {R} ^ {k},
$$

it is easy to verify that $\tilde { \alpha } _ { k + 1 } = a + \beta _ { k + 1 } ^ { 2 } e _ { k } ^ { T } ( T _ { k } - a I _ { k } ) ^ { - 1 } e _ { k }$ works.

# 10.2.6 The Overall Framework

All the necessary tools are now available to obtain sufficiently accurate upper and bounds in (10.2.1). At the bottom of the loop in Algorithm 10.1.1, we use the current tridiagonal (or an augmented version) to compute the nodes and weights for the lower bound rule. The rule is evaluated to obtain b. Likewise, we use the current tridiagonal (or an augmented version) to compute the nodes and weights for the upper bound rule. The rule is evaluated to obtain B. The while loop in Algorithm 10.1.1 can obviously be redesigned to terminate as soon as $B - b$ is sufficiently small.

# Problems

P10.2.1 The Chebyschev polynomials are generated by the recursion $p _ { k } ( x ) = 2 x p _ { k - 1 } ( x ) - p _ { k - 2 } ( x )$ and are orthonormal with respect to $w ( x ) = ( 1 - x ^ { 2 } ) ^ { - 1 / 2 }$ across [−1, 1]. What are the zeros of $p _ { k } ( x ) ?$

P10.2.2 Following the strategy used in §10.2.5, show how to compute $I _ { G R ( b ) }$ and $I _ { G L } ( f )$ .

# Notes and References for §10.2

For complete coverage of the Gauss quadrature/tridiagonal/Lanczos connection, see:

G.H. Golub and G. Meurant (2010). Matrices, Moments, and Quadrature with Applications, Princeton University Press, Princeton, NJ.

Research in this area has a long history:

G.H. Golub (1962). “Bounds for Eigenvalues of Tridiagonal Symmetric Matrices Computed by the LR Method,” Math. Comput. 16, 438–445.

G.H. Golub and J.H. Welsch (1969). “Calculation of Gauss Quadrature Rules,” Math. Comput. 23, 221–230.

G.H. Golub (1974). “Bounds for Matrix Moments,” Rocky Mountain J. Math. 4, 207–211.

C. de Boor and G.H. Golub (1978). “The Numerically Stable Reconstruction of a Jacobi Matrix from Spectral Data,” Lin. Alg. Applic. 21, 245–260.

J. Kautsky and G.H. Golub (1983). “On the Calculation of Jacobi Matrices,” Lin. Alg. Applic. 52/53, 439–455.

M. Berry and G.H. Golub (1991). “Estimating the Largest Singular Values of Large Sparse Matrices via Modified Moments,” Numer. Algs. 1, 353–374.   
D.P. Laurie (1996). “Anti-Gaussian Quadrature Rules,” Math. Comput. 65, 739–747.   
Z. Bai and G.H. Golub (1997). “Bounds for the Trace of the Inverse and the Determinant of Symmetric Positive Definite Matrices,” Annals Numer. Math. 4, 29–38.   
M. Benzi and G.H. Golub (1999). “Bounds for the Entries of Matrix Functions with Applications to Preconditioning,” BIT 39, 417–438.   
D. Calvetti, G. H. Golub, W. B. Gragg, and L. Reichel (2000). “Computation of Gauss–Kronrod Quadrature Rules,” Math. Comput. 69, 1035–1052.   
D.P. Laurie (2001). “Computation of Gauss-Type Quadrature Formulas,” J. Comput. Appl. Math. 127, 201–217.
