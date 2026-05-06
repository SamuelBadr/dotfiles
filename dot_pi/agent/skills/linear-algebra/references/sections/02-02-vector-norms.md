# 2.2 Vector Norms

A norm on a vector space plays the same role as absolute value: it furnishes a distance measure. More precisely, $\mathbb { R } ^ { n }$ together with a norm on $\mathbb { R } ^ { n }$ defines a metric space rendering the familiar notions of neighborhood, open sets, convergence, and continuity.

# 2.2.1 Definitions

A vector norm on $\mathbb { R } ^ { n }$ is a function $f { : \mathbb { R } ^ { n } } \to \mathbb { R }$ that satisfies the following properties:

$$
f (x) \geq 0, \quad x \in \mathbb {R} ^ {n}, \quad (f (x) = 0, \text {   iff   } x = 0),
$$

$$
f (x + y) \leq f (x) + f (y), \quad x, y \in \mathbb {R} ^ {n},
$$

$$
f (\alpha x) = | \alpha | f (x), \quad \alpha \in \mathbb {R}, x \in \mathbb {R} ^ {n}.
$$

We denote such a function with a double bar notation: $f ( x ) = \parallel x \parallel$ . Subscripts on the double bar are used to distinguish between various norms. A useful class of vector

norms are the p-norms defined by

$$
\| x \| _ {p} = \left(| x _ {1} | ^ {p} + \dots + | x _ {n} | ^ {p}\right) ^ {\frac {1}{p}}, \quad p \geq 1. \tag {2.2.1}
$$

The 1−, 2−, and $\infty -$ norms are the most important:

$$
\begin{array}{l} \| x \| _ {1} = | x _ {1} | + \dots + | x _ {n} |, \\ \| x \| _ {2} = \left(| x _ {1} | ^ {2} + \dots + | x _ {n} | ^ {2}\right) ^ {\frac {1}{2}} = \left(x ^ {T} x\right) ^ {\frac {1}{2}}, \\ \| x \| _ {\infty} = \max _ {1 \leq i \leq n} | x _ {i} |. \\ \end{array}
$$

A unit vector with respect to the norm $\parallel \cdot \parallel$ is a vector x that satisfies $\| { \boldsymbol { x } } \| = 1$ .

# 2.2.2 Some Vector Norm Properties

A classic result concerning p-norms is the H¨older inequality:

$$
\left| x ^ {T} y \right| \leq \| x \| _ {p} \| y \| _ {q} \quad \frac {1}{p} + \frac {1}{q} = 1. \tag {2.2.2}
$$

A very important special case of this is the Cauchy-Schwarz inequality:

$$
\left| x ^ {T} y \right| \leq \| x \| _ {2} \| y \| _ {2}. \tag {2.2.3}
$$

All norms on $\mathbb { R } ^ { n }$ are equivalent , i.e., $\operatorname { i f } \parallel \cdot \parallel _ { \alpha }$ and $\| \cdot \| _ { \beta }$ are norms on $\mathbb { R } ^ { n }$ , then there exist positive constants $c _ { 1 }$ and $c _ { 2 }$ such that

$$
c _ {1} \| x \| _ {\alpha} \leq \| x \| _ {\beta} \leq c _ {2} \| x \| _ {\alpha} \tag {2.2.4}
$$

for all $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . For example, if $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , then

$$
\| x \| _ {2} \leq \| x \| _ {1} \leq \sqrt {n} \| x \| _ {2}, \tag {2.2.5}
$$

$$
\| x \| _ {\infty} \leq \| x \| _ {2} \leq \sqrt {n} \| x \| _ {\infty}, \tag {2.2.6}
$$

$$
\| x \| _ {\infty} \leq \| x \| _ {1} \leq n \| x \| _ {\infty}. \tag {2.2.7}
$$

Finally, we mention that the 2-norm is preserved under orthogonal transformation. Indeed, if $Q \in \mathbb { R } ^ { n \times n }$ is orthogonal and $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , then

$$
\parallel Q x \parallel_ {2} ^ {2} = (Q x) ^ {T} (Q x) = (x ^ {T} Q ^ {T}) (Q x) = x ^ {T} (Q ^ {T} Q) x = x ^ {T} x = \parallel x \parallel_ {2} ^ {2}.
$$

# 2.2.3 Absolute and Relative Errors

Suppose $\hat { x } \in \mathbb { R } ^ { n }$ is an approximation to $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . For a given vector norm $\| \cdot \|$ we say that

$$
\epsilon_ {\mathrm{abs}} = \left\| \hat {x} - x \right\|
$$

is the absolute error in ˆx. If $x \neq 0$ , then

$$
\epsilon_ {\mathrm{rel}} = \frac {\parallel \hat {x} - x \parallel}{\parallel x \parallel}
$$

prescribes the relative error in ˆx. Relative error in the ∞-norm can be translated into a statement about the number of correct significant digits in ${ \hat { x } } .$ In particular, if

$$
\frac {\parallel \hat {x} - x \parallel_ {\infty}}{\parallel x \parallel_ {\infty}} \approx 1 0 ^ {- p},
$$

then the largest component of ˆx has approximately p correct significant digits. For example, if $x = [ 1 . 2 3 4 ~ . 0 5 6 7 4 ] ^ { T }$ and $\hat { x } = [ 1 . 2 3 5 . 0 5 1 2 8 ] ^ { T }$ , then $\| \hat { x } - x \| _ { \infty } / \| x \| _ { \infty } \approx$ $. 0 0 4 3 \approx 1 0 ^ { - 3 }$ . Note than $\hat { x } _ { 1 }$ ∞ ∞ has about three significant digits that are correct while only one significant digit in ${ \hat { x } } _ { 2 }$ is correct.

# 2.2.4 Convergence

We say that a sequence $\{ x ^ { ( k ) } \}$ of n-vectors converges to x if

$$
\lim _ {k \to \infty} \| x ^ {(k)} - x \| = 0.
$$

Because of (2.2.4), convergence in any particular norm implies convergence in all norms.

# Problems

P2.2.1 Show that if $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , then li $1 _ { p \to \infty } \parallel x \parallel _ { p } = \parallel x \parallel _ { \infty }$ .

P2.2.2 By considering the inequality $0 \leq ( a x + b y ) ^ { T } ( a x + b y )$ for suitable scalars a and b, prove (2.2.3).

P2.2.3 Verify that $\| { \bf \cdot } \| _ { 1 } , \| { \bf \cdot } \| _ { 2 }$ , and $\| \cdot \| _ { \infty }$ are vector norms.

P2.2.4 Verify (2.2.5)-(2.2.7). When is equality achieved in each result?

P2.2.5 Show that in $\mathbb { R } ^ { n } , x ^ { ( i ) } \to x$ if and only if $x _ { k } ^ { ( i ) } \to x _ { k }$ for $k = 1 { : } n$

P2.2.6 Show that for any vector norm on $\mathbb { R } ^ { n }$ that $\mid \parallel x \parallel - \parallel y \parallel \mid \leq \parallel x - y \parallel ,$ .

P2.2.7 Let $\| \cdot \|$ be a vector norm on $\mathbb { R } ^ { m }$ and assume $A \in \mathbb { R } ^ { m \times n }$ . Show that if rank(A) = n, then $\| x \| _ { A } = \| A x \|$ is a vector norm on $\mathbb { R } ^ { n }$ .

P2.2.8 Let x and y be in $\mathbb { R } ^ { n }$ and define $\psi : \mathbb { R } \to \mathbb { R }$ by $\psi ( \alpha ) = \| { \boldsymbol { x } } - \alpha y \| _ { 2 }$ . Show that $\psi$ is minimized if $\alpha = x ^ { T } y / y ^ { T } y$ .

P2.2.9 Prove or disprove:

$$
v \in \mathbb {R} ^ {n} \Rightarrow \| v \| _ {1} \| v \| _ {\infty} \leq \frac {1 + \sqrt {n}}{2} \| v \| _ {2} ^ {2}.
$$

P2.2.10 If $\boldsymbol { x } \in \mathbb { R } ^ { 3 }$ and $\boldsymbol { y } \in \mathbb { R } ^ { 3 }$ then it can be shown that $| x ^ { T } y | = \parallel x \parallel _ { 2 } \parallel y \parallel _ { 2 } | \cos ( \theta ) |$ where $\theta$ is the angle between x and y. An analogous result exists for the cross product defined by

$$
x \times y = \left[ \begin{array}{l} x _ {2} y _ {3} - x _ {3} y _ {2} \\ x _ {3} y _ {1} - x _ {1} y _ {3} \\ x _ {1} y _ {2} - x _ {2} y _ {1} \end{array} \right].
$$

In particular, $\parallel x \times y \parallel _ { 2 } = \parallel x \parallel _ { 2 } \parallel y \parallel _ { 2 } | \sin ( \theta ) |$ . Prove this.

P2.2.11 Suppose $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $y \in \mathbb { R } ^ { m }$ . Show that

$$
\parallel x \otimes y \parallel_ {p} = \parallel x \parallel_ {p} \parallel y \parallel_ {p}
$$

for $p = 1 , 2$ , and ∞.

# Notes and References for 2.2

Although a vector norm is “just” a generalization of the absolute value concept, there are some noteworthy subtleties:

J.D. Pryce (1984). “A New Measure of Relative Error for Vectors,” SIAM J. Numer. Anal. 21, 202–221.
