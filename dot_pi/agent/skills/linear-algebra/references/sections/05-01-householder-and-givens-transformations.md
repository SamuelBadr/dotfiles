# 5.1 Householder and Givens Transformations

Recall that $Q \in \mathbb { R } ^ { m \times m }$ is orthogonal if

$$
Q ^ {T} Q = Q Q ^ {T} = I _ {m}.
$$

Orthogonal matrices have an important role to play in least squares and eigenvalue computations. In this section we introduce Householder reflections and Givens rotations, the key players in this game.

# 5.1.1 A 2-by-2 Preview

It is instructive to examine the geometry associated with rotations and reflections at the m = 2 level. A 2-by-2 orthogonal matrix Q is a rotation if it has the form

$$
Q = \left[ \begin{array}{c c} \cos (\theta) & \sin (\theta) \\ - \sin (\theta) & \cos (\theta) \end{array} \right].
$$

If $y = Q ^ { T } x ,$ , then y is obtained by rotating x counterclockwise through an angle θ. A 2-by-2 orthogonal matrix Q is a reflection if it has the form

$$
Q = \left[ \begin{array}{c c} \cos (\theta) & \sin (\theta) \\ \sin (\theta) & - \cos (\theta) \end{array} \right].
$$

If $y = Q ^ { T } x = Q x$ , then y is obtained by reflecting the vector x across the line defined by

$$
S = \text { span } \left\{\left[ \begin{array}{c} \cos (\theta / 2) \\ \sin (\theta / 2) \end{array} \right] \right\}.
$$

Reflections and rotations are computationally attractive because they are easily constructed and because they can be used to introduce zeros in a vector by properly choosing the rotation angle or the reflection plane.

# 5.1.2 Householder Reflections

Let $v \in \mathbb { R } ^ { m }$ be nonzero. An m-by-m matrix P of the form

$$
P = I - \beta v v ^ {T}, \quad \beta = \frac {2}{v ^ {T} v} \tag {5.1.1}
$$

is a Householder reflection. (Synonyms are Householder matrix and Householder transformation.) The vector v is the Householder vector. If a vector x is multiplied by P , then it is reflected in the hyperplane span $\{ v \} ^ { \perp }$ . It is easy to verify that Householder matrices are symmetric and orthogonal.

Householder reflections are similar to Gauss transformations introduced in §3.2.1 in that they are rank-1 modifications of the identity and can be used to zero selected components of a vector. In particular, suppose we are given $0 \neq x \in \mathbb { R } ^ { m }$ and want

$$
P x = \left(I - \frac {2 v v ^ {T}}{v ^ {T} v}\right) x = x - \frac {2 v ^ {T} x}{v ^ {T} v} v
$$

to be a multiple of $e _ { 1 } = I _ { m } ( : , 1 )$ . From this we conclude that $v \in { \mathrm { s p a n } } \{ x , e _ { 1 } \}$ . Setting

$$
v = x + \alpha e _ {1}
$$

gives

$$
v ^ {T} x = x ^ {T} x + \alpha x _ {1}
$$

and

$$
v ^ {T} v = x ^ {T} x + 2 \alpha x _ {1} + \alpha^ {2}.
$$

Thus,

$$
\begin{array}{l} P x = \left(1 - 2 \frac {x ^ {T} x + \alpha x _ {1}}{x ^ {T} x + 2 \alpha x _ {1} + \alpha^ {2}}\right) x - 2 \alpha \frac {v ^ {T} x}{v ^ {T} v} e _ {1} \\ = \left(\frac {\alpha^ {2} - \| x \| _ {2} ^ {2}}{x ^ {T} x + 2 \alpha x _ {1} + \alpha^ {2}}\right) x - 2 \alpha \frac {v ^ {T} x}{v ^ {T} v} e _ {1}. \\ \end{array}
$$

In order for the coefficient of x to be zero, we set $\alpha = \pm \| \boldsymbol { x } \| _ { 2 }$ for then

$$
v = x \pm \| x \| _ {2} e _ {1} \Rightarrow P x = \left(I - 2 \frac {v v ^ {T}}{v ^ {T} v}\right) x = \mp \| x \| _ {2} e _ {1}. \tag {5.1.2}
$$

It is this simple determination of v that makes the Householder reflections so useful.

# 5.1.3 Computing the Householder Vector

There are a number of important practical details associated with the determination of a Householder matrix, i.e., the determination of a Householder vector. One concerns the choice of sign in the definition of v in (5.1.2). Setting

$$
v _ {1} = x _ {1} - \parallel x \parallel_ {2}
$$

leads to the nice property that $P x$ is a positive multiple of $e _ { 1 }$ . But this recipe is dangerous if x is close to a positive multiple of $e _ { 1 }$ because severe cancellation would occur. However, the formula

$$
v _ {1} = x _ {1} - \parallel x \parallel_ {2} = \frac {x _ {1} ^ {2} - \parallel x \parallel_ {2} ^ {2}}{x _ {1} + \parallel x \parallel_ {2}} = \frac {- (x _ {2} ^ {2} + \cdots + x _ {n} ^ {2})}{x _ {1} + \parallel x \parallel_ {2}}
$$

suggested by Parlett (1971) does not suffer from this defect in the $x _ { 1 } > 0$ case.

In practice, it is handy to normalize the Householder vector so that $v ( 1 ) = 1$ . This permits the storage of v(2:m) where the zeros have been introduced in x, i.e., x(2:m). We refer to v(2:m) as the essential part of the Householder vector. Recalling that $\beta = 2 / v ^ { T } v$ and letting length(x) specify vector dimension, we may encapsulate the overall process as follows:

Algorithm 5.1.1 (Householder Vector) $\mathrm { G i v e n } x \in \mathbb { R } ^ { m }$ , this function computes $v \in \mathbb { R } ^ { m }$ with $v ( 1 ) = 1$ and $\beta \in \mathbb { R }$ such that $\begin{array} { r } { \dot { P } = I _ { m } - \beta v v ^ { T } } \end{array}$ is orthogonal and $\ P x = \| \boldsymbol { x } \| _ { 2 } e _ { 1 }$ .

function $\left\lceil v , \beta \right\rceil = \mathsf { h o u s e } ( x )$

$$
m = \mathsf {l e n g t h} (x), \sigma = x (2: m) ^ {T} x (2: m), v = \left[ \begin{array}{c} 1 \\ x (2: m) \end{array} \right]
$$

$\mathbf { i f } \ \sigma = 0 \ \mathrm { a n d } \ x ( 1 ) > = 0$

$$
\beta = 0
$$

elseif $\sigma = 0 \mathrm { ~ } \& \mathrm { ~ } x ( 1 ) < 0$

$$
\beta = - 2
$$

else

$$
\mu = \sqrt {x (1) ^ {2} + \sigma}
$$

$\mathbf { i f } \ x ( 1 ) < = 0$

$$
v (1) = x (1) - \mu
$$

$$
v (1) = - \sigma / (x (1) + \mu)
$$

end

$$
\beta = 2 v (1) ^ {2} / (\sigma + v (1) ^ {2})
$$

$$
v = v / v (1)
$$

end

Here, length(·) returns the dimension of a vector. This algorithm involves about 3m flops. The computed Householder matrix that is orthogonal to machine precision, a concept discussed below.

# 5.1.4 Applying Householder Matrices

It is critical to exploit structure when applying $P = I - \beta v v ^ { T }$ to a matrix A. Premultiplication involves a matrix-vector product and a rank-1 update:

$$
P A = (I - \beta v v ^ {T}) A = A - (\beta v) (v ^ {T} A).
$$

The same is true for post-multiplication,

$$
A P = A (I - \beta v v ^ {T}) = A - (A v) (\beta v) ^ {T}.
$$

In either case, the update requires 4mn flops if $A \in \mathbb { R } ^ { m \times n }$ . Failure to recognize this and to treat P as a general matrix increases work by an order of magnitude. Householder updates never entail the explicit formation of the Householder matrix.

In a typical situation, house is applied to a subcolumn or subrow of a matrix and $( I - \beta v v ^ { T } )$ is applied to a submatrix. For example, if $A \in \mathbb { R } ^ { m \times n } , \ 1 \leq j < n$ , and $A ( j { : } m , 1 { : } j - 1 )$ is zero, then the sequence

$$
[ v, \beta ] = \text { house } (A (j: m, j))
$$

$$
A (j: m, j: n) = A (j: m, j: n) - (\beta v) \left(v ^ {T} A (j: m, j: n)\right)
$$

$$
A (j + 1: m, j) = v (2: m - j + 1)
$$

applies $( I _ { m - j + 1 } - \beta v v ^ { T } )$ to $A ( j { : } m , 1 { : } n )$ and stores the essential part of v where the “new” zeros are introduced.

# 5.1.5 Roundoff Properties

The roundoff properties associated with Householder matrices are very favorable. Wilkinson (AEP, pp. 152–162) shows that house produces a Householder vector $\hat { v }$ that is very close to the exact v. If $\hat { P } = I - 2 \hat { v } \hat { v } ^ { T } / \hat { v } ^ { \hat { T } } \hat { v }$ then

$$
\| \hat {P} - P \| _ {2} = O (\mathbf {u}).
$$

Moreover, the computed updates with $\hat { P }$ are close to the exact updates with $P$ :

$$
\mathfrak {f l} (\hat {P} A) = P (A + E), \quad \| E \| _ {2} = O (\mathbf {u} \| A \| _ {2}),
$$

$$
\mathsf {f l} (A \hat {P}) = (A + E) P, \quad \| E \| _ {2} = O (\mathbf {u} \| A \| _ {2}).
$$

For a more detailed analysis, see Higham(ASNA, pp. 357–361).

# 5.1.6 The Factored-Form Representation

Many Householder-based factorization algorithms that are presented in the following sections compute products of Householder matrices

$$
Q = Q _ {1} Q _ {2} \dots Q _ {n} \quad Q _ {j} = I _ {m} - \beta_ {j} v ^ {(j)} [ v ^ {(j)} ] ^ {T} \tag {5.1.3}
$$

where $n \leq m$ and each $v ^ { ( j ) }$ has the form

$$
v ^ {(j)} = [ \underbrace {0 , 0 , \ldots 0} _ {j - 1}, 1 v _ {j + 1} ^ {(j)}, \ldots , v _ {m} ^ {(j)} ] ^ {T}.
$$

It is usually not necessary to compute Q explicitly even if it is involved in subsequent calculations. For example, if $\boldsymbol { C } \in \mathbb { R } ^ { m \times p }$ and we wish to compute $Q ^ { T } C$ , then we merely execute the loop

for $j = 1 { : } n$

$$
C = Q _ {j} C
$$

end

The storage of the Householder vectors $\boldsymbol { v } ^ { ( 1 ) } \cdots \boldsymbol { v } ^ { ( n ) }$ and the corresponding $\beta _ { j }$ amounts to a factored-form representation of $Q$ .

To illustrate the economies of the factored-form representation, suppose we have an array A and that for $j = 1 { : } n , A ( j + 1 { : } m , j )$ houses $\boldsymbol { v } ^ { ( j ) } ( j + 1 { : } m )$ , the essential part of the jth Householder vector. The overwriting of $\boldsymbol { C } \in \mathbb { R } ^ { m \times p }$ with $Q ^ { T } C$ can then be implemented as follows:

for $j = 1 { : } n$

$$
v (j: m) = \left[ \begin{array}{c} 1 \\ A (j + 1: m, j) \end{array} \right]
$$

$$
\beta_ {j} = 2 / (1 + \| A (j + 1: m, j) \| _ {2} ^ {2} \tag {5.1.4}
$$

$$
C (j: m,:) = C (j: m,:) - \left(\beta_ {j} \cdot v (j: m)\right) \cdot \left(v (j: m) ^ {T} C (j: m,:)\right)
$$

end

This involves about $p n ( 2 m - n )$ flops. If Q is explicitly represented as an m-by-m matrix, then $Q ^ { T } C$ would involve $2 m ^ { 2 } p$ flops. The advantage of the factored form representation is apparant if $n < < m$ .

Of course, in some applications, it is necessary to explicitly form Q (or parts of it). There are two possible algorithms for computing the matrix $Q$ in (5.1.3):

<table><tr><td>Forward accumulation</td><td>Backward accumulation</td></tr><tr><td> $Q = I_{m}$ </td><td> $Q = I_{m}$ </td></tr><tr><td>for  $j = 1:n$ </td><td>for  $j = n: -1:1$ </td></tr><tr><td> $Q = QQ_{j}$ </td><td> $Q = Q_{j}Q$ </td></tr><tr><td>end</td><td>end</td></tr></table>

Recall that the leading $( j \mathrm { ~ - ~ } 1 ) – \mathrm { b y } – ( j \mathrm { ~ - ~ } 1 )$ portion of $Q _ { j }$ is the identity. Thus, at the beginning of backward accumulation, Q is “mostly the identity” and it gradually becomes full as the iteration progresses. This pattern can be exploited to reduce the number of required flops. In contrast, Q is full in forward accumulation after the first step. For this reason, backward accumulation is cheaper and the strategy of choice. Here are the details with the proviso that we only need Q(:, 1:k) where $1 \leq k \leq m \colon$

$$
Q = I _ {m} (:, 1: k)
$$

for $j = n \colon - 1 \colon 1$

$$
v (j: m) = \left[ \begin{array}{c} 1 \\ A (j + 1: m, j) \end{array} \right] \tag {5.1.5}
$$

$$
\beta_ {j} = 2 / (1 + \parallel A (j + 1: m, j) \parallel_ {2} ^ {2}
$$

$$
Q (j: m, j: k) = Q (j: m, j: k) - (\beta_ {j} v (j: m)) (v (j: m) ^ {T} Q (j: m, j: k))
$$

end

This involves about 4mnk $- 2 ( m + k ) n ^ { 2 } + ( 4 / 3 ) n ^ { 3 }$ flops.

# 5.1.7 The WY Representation

Suppose $Q = Q _ { 1 } \cdot \cdot \cdot Q _ { r }$ is a product of m-by-m Householder matrices. Since each $Q _ { j }$ is a rank-1 modification of the identity, it follows from the structure of the Householder vectors that $Q$ is a rank-r modification of the identity and can be written in the form

$$
Q = I _ {m} - W Y ^ {T} \tag {5.1.6}
$$

where W and Y are m-by-r matrices. The key to computing the WY representation (5.1.6) is the following lemma.

Lemma 5.1.1. Suppose $Q \ = \ I _ { m } - W Y ^ { T }$ is an m-by-m orthogonal matrix with W, $Y \in \mathbb { R } ^ { m \times j }$ . If $P = I _ { m } - \beta v v ^ { T }$ with $v \in \mathbb { R } ^ { m }$ and $z = \beta Q v$ , then

$$
Q _ {+} = Q P = I _ {m} - W _ {+} Y _ {+} ^ {T}
$$

where $W _ { + } = \left[ W \mid z \right]$ and $Y _ { + } = \left[ Y \mid v \right]$ are each $m - b y - ( j + 1 )$ .

Proof. Since

$$
Q P = \left(I _ {m} - W Y ^ {T}\right) \left(I _ {m} - \beta v v ^ {T}\right) = I _ {m} - W Y ^ {T} - \beta Q v v ^ {T}
$$

it follows from the definition of z that

$$
Q _ {+} = I _ {m} - W Y ^ {T} - z v ^ {T} = I _ {m} - \left[ W \mid z \right] \left[ Y \mid v \right] ^ {T} = I _ {m} - W _ {+} Y _ {+} ^ {T}. \quad \square
$$

By repeatedly applying the lemma, we can transition from a factored-form representation to a block representation.

Algorithm 5.1.2 Suppose Q = Q1 · · · Qr where the Qj = Im − βjv(j) v(j)T $Q = Q _ { 1 } \cdot \cdot \cdot Q _ { r }$ $Q _ { j } = I _ { m } - \beta _ { j } v ^ { ( j ) } v ^ { ( j ) ^ { T } }$ are stored in factored form. This algorithm computes matrices $W , Y \in \mathbb { R } ^ { \bar { m } \times r }$ such that $Q =$ $I _ { m } - W Y ^ { T }$ .

$$
Y = v ^ {(1)}; W = \beta_ {1} v ^ {(1)}
$$

for $j = 2 { : } r$

$$
z = \beta_ {j} (I _ {m} - W Y ^ {T}) v ^ {(j)}
$$

$$
W = [ W \mid z ]
$$

$$
Y = \left[ Y \mid v ^ {(j)} \right]
$$

end

This algorithm involves about $2 r ^ { 2 } m - 2 r ^ { 3 } / 3$ flops if the zeros in the $v ^ { ( j ) }$ are exploited. Note that Y is merely the matrix of Householder vectors and is therefore unit lower triangular. Clearly, the central task in the generation of the WY representation (5.1.6) is the computation of the matrix W .

The block representation for products of Householder matrices is attractive in situations where $Q$ must be applied to a matrix. Suppose $\boldsymbol { C } \in \mathbb { R } ^ { m \times p }$ . It follows that the operation

$$
C = Q ^ {T} C = (I _ {m} - W Y ^ {T}) ^ {T} C = C - Y (W ^ {T} C)
$$

is rich in level-3 operations. On the other hand, if Q is in factored form, then the formation of $Q ^ { T } C$ is just rich in the level-2 operations of matrix-vector multiplication and outer product updates. Of course, in this context, the distinction between level-2 and level-3 diminishes as C gets narrower.

We mention that the WY representation (5.1.6) is not a generalized Householder transformation from the geometric point of view. True block reflectors have the form

$$
Q = I - 2 V V ^ {T}
$$

where $V \in \mathbb { R } ^ { n \times r }$ satisfies $V ^ { T } V = I _ { r }$ . See Schreiber and Parlett (1987).

# 5.1.8 Givens Rotations

Householder reflections are exceedingly useful for introducing zeros on a grand scale, e.g., the annihilation of all but the first component of a vector. However, in calculations where it is necessary to zero elements more selectively, Givens rotations are the transformation of choice. These are rank-2 corrections to the identity of the form

$$
G (i, k, \theta) = \left[ \begin{array}{c c c c c c c} 1 & \dots & 0 & \dots & 0 & \dots & 0 \\ \vdots & \ddots & \vdots & & \vdots & & \vdots \\ 0 & \dots & c & \dots & s & \dots & 0 \\ \vdots & & \vdots & \ddots & \vdots & & \vdots \\ 0 & \dots & - s & \dots & c & \dots & 0 \\ \vdots & & \vdots & & \vdots & \ddots & \vdots \\ 0 & \dots & 0 & \dots & 0 & \dots & 1 \\ & & i & & k \end{array} \right] \begin{array}{l} i \\ k \end{array} \tag {5.1.7}
$$

where $c = \cos ( \theta )$ and $s = \sin ( \theta )$ for some θ. Givens rotations are clearly orthogonal.

Premultiplication by $G ( i , k , \theta ) ^ { T }$ amounts to a counterclockwise rotation of θ radians in the (i, k) coordinate plane. Indeed, if $\boldsymbol { x } \in \mathbb { R } ^ { m }$ and

$$
y = G (i, k, \theta) ^ {T} x,
$$

then

$$
y _ {j} = \left\{ \begin{array}{c c} c x _ {i} - s x _ {k}, & j = i, \\ s x _ {i} + c x _ {k}, & j = k, \\ x _ {j}, & j \neq i, k. \end{array} \right..
$$

From these formulae it is clear that we can force $y _ { k }$ to be zero by setting

$$
c = \frac {x _ {i}}{\sqrt {x _ {i} ^ {2} + x _ {k} ^ {2}}}, \quad s = \frac {- x _ {k}}{\sqrt {x _ {i} ^ {2} + x _ {k} ^ {2}}}. \tag {5.1.8}
$$

Thus, it is a simple matter to zero a specified entry in a vector by using a Givens rotation. In practice, there are better ways to compute c and s than (5.1.8), e.g.,

Algorithm 5.1.3 Given scalars a and b, this function computes $c = \cos ( \theta )$ and $s = \sin ( \theta )$ so

$$
\left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} \left[ \begin{array}{c} a \\ b \end{array} \right] = \left[ \begin{array}{c} r \\ 0 \end{array} \right].
$$

$$
\begin{array}{l} \text { function } [ c, s ] = \text { givens } (a, b) \\ \text { if } b = 0 \\ c = 1; s = 0 \\ \text { if } | b | > | a | \\ \tau = - a / b; s = 1 / \sqrt {1 + \tau^ {2}}; c = s \tau \\ \tau = - b / a; c = 1 / \sqrt {1 + \tau^ {2}}; s = c \tau \\ \end{array}
$$

This algorithm requires 5 flops and a single square root. Note that inverse trigonometric functions are not involved.

# 5.1.9 Applying Givens Rotations

It is critical that the simple structure of a Givens rotation matrix be exploited when it is involved in a matrix multiplication. Suppose $A \in \mathbb { R } ^ { m \times n } , c = \cos ( \theta )$ , and $s = \sin ( \theta )$ . If $G ( i , k , \theta ) \in \mathbb { R } ^ { m \times m }$ , then the update $A = G ( i , k , \theta ) ^ { T } A$ affects just two rows,

$$
A ([ i, k ],:) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A ([ i, k ],:),
$$

and involves 6n flops:

$$
\begin{array}{l} \text { for } j = 1: n \\ \tau_ {1} = A (i, j) \\ \tau_ {2} = A (k, j) \\ A (i, j) = c \tau_ {1} - s \tau_ {2} \\ A (k, j) = s \tau_ {1} + c \tau_ {2} \\ \end{array}
$$

end

Likewise, if $G ( i , k , \theta ) \in \mathbb { R } ^ { n \times n }$ , then the update $A = A G ( i , k , \theta )$ affects just two columns,

$$
A (:, [ i, k ]) = A (:, [ i, k ]) \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right],
$$

and involves 6m flops:

$$
\begin{array}{l} \text {   for   } j = 1: m \\ \tau_ {1} = A (j, i) \\ \tau_ {2} = A (j, k) \\ A (j, i) = c \tau_ {1} - s \tau_ {2} \\ A (j, k) = s \tau_ {1} + c \tau_ {2} \\ \end{array}
$$

end

# 5.1.10 Roundoff Properties

The numerical properties of Givens rotations are as favorable as those for Householder reflections. In particular, it can be shown that the computed ˆc and ˆs in givens satisfy

$$
\begin{array}{l} \hat {c} = c (1 + \epsilon_ {c}), \quad \epsilon_ {c} = O (\mathbf {u}), \\ \hat {s} = s (1 + \epsilon_ {s}), \quad \epsilon_ {s} = O (\mathbf {u}). \\ \end{array}
$$

If ˆc and ˆs are subsequently used in a Givens update, then the computed update is the exact update of a nearby matrix:

$$
\mathsf {f l} [ \hat {G} (i, k, \theta) ^ {T} A ] = G (i, k, \theta) ^ {T} (A + E), \quad \| E \| _ {2} \approx \mathbf {u} \| A \| _ {2},
$$

$$
\mathsf {f l} [ A \hat {G} (i, k, \theta) ] = (A + E) G (i, k, \theta), \quad \| E \| _ {2} \approx \mathbf {u} \| A \| _ {2}.
$$

Detailed error analysis of Givens rotations may be found in Wilkinson (AEP, pp. 131- 39), Higham(ASNA, pp. 366–368), and Bindel, Demmel, Kahan, and Marques (2002).

# 5.1.11 Representing Products of Givens Rotations

Suppose $Q = G _ { 1 } \cdot \cdot \cdot G _ { t }$ is a product of Givens rotations. As with Householder reflections, it is sometimes more economical to keep Q in factored form rather than to compute explicitly the product of the rotations. Stewart (1976) has shown how to do this in a very compact way. The idea is to associate a single floating point number $\rho$ with each rotation. Specifically, if

$$
Z = \left[ \begin{array}{c c} {c} & {s} \\ {- s} & {c} \end{array} \right], \qquad c ^ {2} + s ^ {2} = 1,
$$

then we define the scalar $\rho$ by

$\mathbf { i f } \ c = 0$

$$
\rho = 1
$$

elseif $| s | < | c |$

$$
\rho = \operatorname{sign} (c) \cdot s / 2 \tag {5.1.9}
$$

$$
\rho = 2 \cdot \mathrm{sign} (s) / c
$$

end

Essentially, this amounts to storing $s / 2$ if the sine is smaller and $2 / c$ if the cosine is smaller. With this encoding, it is possible to reconstruct Z (or −Z) as follows:

if $\rho = 1$

$$
c = 0; s = 1
$$

elseif $| \rho | < 1$

$$
s = 2 \rho ; c = \sqrt {1 - s ^ {2}} \tag {5.1.10}
$$

$$
c = 2 / \rho ; s = \sqrt {1 - c ^ {2}}
$$

end

Note that the reconstruction of −Z is not a problem, for if $Z$ introduces a strategic zero then so does $- Z$ . The reason for essentially storing the smaller of c and s is that the formula $\sqrt { 1 - x ^ { 2 } }$ renders poor results if x is near unity. More details may be found in Stewart (1976). Of course, to “reconstruct” $G ( i , k , \theta )$ we need i and k in addition to the associated $\rho .$ This poses no difficulty if we agree to store $\rho$ in the (i, k) entry of some array.

# 5.1.12 Error Propagation

An m-by-m floating point matrix $\hat { Q }$ is orthogonal to working precision if there exists an orthogonal $Q \in \mathbb { R } ^ { m \times m }$ such that

$$
\| \hat {Q} - Q \| = O (\mathbf {u}).
$$

A corollary of this is that

$$
\| \hat {Q} ^ {T} \hat {Q} - I _ {m} \| = O (\mathbf {u}).
$$

The matrices defined by the floating point output of house and givens are orthogonal to working precision.

In many applications, sequences of Householders and/or Given transformations are generated and applied. In these settings, the rounding errors are nicely bounded. To be precise, suppose $A = A _ { 0 } \in \mathbb { R } ^ { m \times n }$ is given and that matrices $A _ { 1 } , \dotsc , A _ { p } = B$ are generated via the formula

$$
A _ {k} = \mathsf {f l} (\hat {Q} _ {k} A _ {k - 1} \hat {Z} _ {k}), \qquad k = 1: p  .
$$

Assume that the above Householder and Givens algorithms are used for both the generation and application of the $\hat { Q } _ { k }$ and $\hat { Z } _ { k }$ . Let $Q _ { k }$ and $Z _ { k }$ be the orthogonal matrices that would be produced in the absence of roundoff. It can be shown that

$$
B = (Q _ {p} \dots Q _ {1}) (A + E) (Z _ {1} \dots Z _ {p}), \tag {5.1.11}
$$

where $\| E \| _ { 2 } \ \leq \ c \cdot \mathbf { u } \| \ A \ \| _ { 2 }$ and c is a constant that depends mildly on n, m, and p. In other words, B is an exact orthogonal update of a matrix near to A. For a comprehensive error analysis of Householder and Givens computations, see Higham (ASNA, §19.3, §19.6).

# 5.1.13 The Complex Case

Most of the algorithms that we present in this book have complex versions that are fairly straightforward to derive from their real counterparts. (This is not to say that everything is easy and obvious at the implementation level.) As an illustration we briefly discuss complex Householder and complex Givens transformations.

Recall that if $A = ( a _ { i j } ) \in \mathbb { C } ^ { m \times n }$ , then $B = A ^ { H } \in \mathbb { C } ^ { n \times m }$ is its conjugate transpose. The 2-norm of a vector $x \in \mathbb { C } ^ { n }$ is defined by

$$
\parallel x \parallel_ {2} ^ {2} = x ^ {H} x = | x _ {1} | ^ {2} + \dots + | x _ {n} | ^ {2}
$$

and $Q \in \mathbb { C } ^ { n \times n }$ is unitary if $Q ^ { H } Q = I _ { n }$ . Unitary matrices preserve the 2-norm.

A complex Householder transformation is a unitary matrix of the form

$$
P = I _ {m} - \beta v v ^ {H}, \qquad 0 \neq v \in \mathbb {C} ^ {m},
$$

where $\beta = 2 / v ^ { H } v$ . Given a nonzero vector $\boldsymbol { x } \in \mathbb { C } ^ { m }$ , it is easy to determine v so that if $y = P x$ , then $y ( 2 { : } m ) = 0$ . Indeed, if

$$
x _ {1} = r e ^ {i \theta}
$$

where $r , \theta \in \mathbb { R }$ and

$$
v = x \pm e ^ {i \theta} \| x \| _ {2} e _ {1}, \qquad e _ {1} = I _ {m} (:, 1),
$$

then $P x = \mp e ^ { i \theta } \parallel x \parallel _ { 2 } e _ { 1 }$ . The sign can be determined to maximize $\parallel v \parallel _ { 2 }$ for the sake of stability.

Regarding complex Givens rotations, it is easy to verify that a 2-by-2 matrix of the form

$$
Q = \left[ \begin{array}{c c} \cos (\theta) & \sin (\theta) e ^ {i \phi} \\ - \sin (\theta) e ^ {- i \phi} & \cos (\theta) \end{array} \right]
$$

where $\theta , \phi \in \mathbb { R }$ is unitary. We show how to compute $c = \cos ( \theta )$ and $s = \sin ( \theta ) e ^ { i \phi }$ so that

$$
\left[ \begin{array}{c c} c & s \\ - \bar {s} & c \end{array} \right] ^ {H} \left[ \begin{array}{l} u \\ v \end{array} \right] = \left[ \begin{array}{l} r \\ 0 \end{array} \right] \tag {5.1.12}
$$

where $u = u _ { 1 } + i u _ { 2 }$ and $v = v _ { 1 } + i v _ { 2 }$ are given complex numbers. First, givens is applied to compute real cosine-sine pairs $\{ c _ { \alpha } , s _ { \alpha } \} , \{ c _ { \beta } , s _ { \beta } \}$ , and $\big \{ c _ { \theta } , s _ { \theta } \big \}$ so that

$$
\left[ \begin{array}{c c} c _ {\alpha} & s _ {\alpha} \\ - s _ {\alpha} & c _ {\alpha} \end{array} \right] ^ {T} \left[ \begin{array}{c} u _ {1} \\ u _ {2} \end{array} \right] = \left[ \begin{array}{c} r _ {u} \\ 0 \end{array} \right],
$$

$$
\left[ \begin{array}{c c} c _ {\beta} & s _ {\beta} \\ - s _ {\beta} & c _ {\beta} \end{array} \right] ^ {T} \left[ \begin{array}{c} v _ {1} \\ v _ {2} \end{array} \right] = \left[ \begin{array}{c} r _ {v} \\ 0 \end{array} \right],
$$

and

$$
\left[ \begin{array}{c c} c _ {\theta} & s _ {\theta} \\ - s _ {\theta} & c _ {\theta} \end{array} \right] ^ {T} \left[ \begin{array}{c} r _ {u} \\ r _ {v} \end{array} \right] = \left[ \begin{array}{c} r \\ 0 \end{array} \right].
$$

Note that $u = r _ { u } e ^ { - i \alpha }$ and $v = r _ { v } e ^ { - i \beta }$ . If we set

$$
e ^ {i \phi} = e ^ {i (\beta - \alpha)} = (c _ {\alpha} c _ {\beta} + s _ {\alpha} s _ {\beta}) + i (c _ {\alpha} s _ {\beta} - c _ {\beta} s _ {\alpha}),
$$

$c = c _ { \theta }$ , and $s = s _ { \theta } e ^ { i \phi }$ , then

$$
\bar {s} u + c v = s _ {\theta} e ^ {- i \phi} r _ {u} e ^ {- i \alpha} + c _ {\theta} r _ {v} e ^ {- i \beta} = e ^ {- i \beta} (s _ {\theta} r _ {u} + c _ {\theta} r _ {v}) = 0
$$

which confirms (5.1.12).

# Problems

P5.1.1 Let x and y be nonzero vectors in $\mathbb { R } ^ { m }$ . Give an algorithm for determining a Householder matrix P such that $P x$ is a multiple of y.

P5.1.2 Use Householder matrices to show that det $( I + x y ^ { T } ) = 1 + x ^ { T } y$ where x and y are given m-vectors.

P5.1.3 (a) Assume that $x , y \in \mathbb { R } ^ { 2 }$ have unit 2-norm. Give an algorithm that computes a Givens rotation Q so that $y = Q ^ { T } x$ . Make effective use of givens. (b) Suppose x and y are unit vectors in $\mathbb { R } ^ { m }$ . Give an algorithm using Givens transformations which computes an orthogonal Q such that $Q ^ { T } x = y$ .

P5.1.4 By generalizing the ideas in §5.1.11, develop a compact representation scheme for complex givens rotations.

P5.1.5 Suppose that $Q = I { - } Y T Y ^ { T }$ is orthogonal where $Y \in \mathbb { R } ^ { m \times j }$ and $T \in \mathbb { R } ^ { j \times j }$ is upper triangular. Show that if $Q _ { + } = Q P$ where $P = I - 2 v v ^ { T } \big / { v ^ { T } v }$ is a Householder matrix, then $Q _ { + }$ can be expressed in the form $Q _ { + } = I - Y _ { + } T _ { + } Y _ { + } ^ { T }$ where $Y _ { + } \stackrel { \cdot } { \in } \mathbb { R } ^ { m \times ( j + 1 ) }$ and $T _ { + } \in \mathbb { R } ^ { ( j + 1 ) \times ( j + 1 ) }$ is upper triangular. This is the main idea behind the compact WY representation. See Schreiber and Van Loan (1989).

P5.1.6 Suppose $Q _ { 1 } = I _ { m } - Y _ { 1 } T _ { 1 } Y _ { 1 }$ and $Q _ { 2 } \ = \ I _ { m } - Y _ { 2 } T _ { 2 } Y _ { 2 } ^ { T }$ are orthogonal where $Y _ { 1 } \in \mathbb { R } ^ { m \times r _ { 1 } }$ , $Y _ { 2 } \in \mathbb { R } ^ { m \times { r _ { 2 } } } , T _ { 1 } \in \mathbb { R } ^ { r _ { 1 } \times { r _ { 1 } } }$ , and $T _ { 2 } \in \mathbb { R } ^ { r _ { 2 } \times r _ { 2 } }$ . Assume that T1 and $T _ { 2 }$ are upper triangular. Show how to compute $Y \in \mathbb { R } ^ { m \times r }$ and upper triangular $T \in \mathbb { R } ^ { r \times r }$ with $r = r _ { 1 } + r _ { 2 }$ so that $Q _ { 2 } { \bar { Q _ { 1 } } } = I _ { m } - Y T Y ^ { T }$ .

P5.1.7 Give a detailed implementation of Algorithm 5.1.2 with the assumption that $\boldsymbol { v } ^ { ( j ) } ( j + 1 { : } m )$ , the essential part of the jth Householder vector, is stored in $A ( j + 1 { : } m , j )$ . Since Y is effectively represented in A, your procedure need only set up the W matrix.

P5.1.8 Show that if S is skew-symmetric $( S ^ { T } = - S )$ , then $Q = ( I + S ) ( I - S ) ^ { - 1 }$ is orthogonal. (The matrix Q is called the Cayley transform of S.) Construct a rank-2 S so that if x is a vector, then Qx is zero except in the first component.

P5.1.9 Suppose $P \in \mathbb { R } ^ { m \times m }$ satisfies $\parallel P ^ { T } P - I _ { m } \parallel _ { 2 } = \epsilon < 1$ . Show that all the singular values of P are in the interval $[ 1 - \epsilon , 1 + \epsilon ]$ and that $\| \ b { P } - \ b { U V } ^ { T } \| _ { 2 } \le \epsilon$ where $P = U \Sigma V ^ { T }$ is the SVD of $P ,$ .

P5.1.10 Suppose $A \in \mathbb { R } ^ { 2 \times 2 }$ . Under what conditions is the closest rotation to A closer than the closest reflection to A? Work with the Frobenius norm.

P5.1.11 How could Algorithm 5.1.3 be modified to ensure $r \geq 0 ?$

P5.1.12 (Fast Givens Transformations) Suppose

$$
x = \left[ \begin{array}{l} x _ {1} \\ x _ {2} \end{array} \right] \qquad \text {and} \qquad D = \left[ \begin{array}{l l} d _ {1} & 0 \\ 0 & d _ {2} \end{array} \right]
$$

with $d _ { 1 }$ and $d _ { 2 }$ positive. Show how to compute

$$
M _ {1} = \left[ \begin{array}{c c} \beta_ {1} & 1 \\ 1 & \alpha_ {1} \end{array} \right]
$$

so that if $y = M _ { 1 } x$ and $\tilde { D } = M _ { 1 } ^ { T } D M _ { 1 }$ , then $y _ { 2 } = 0$ and $\tilde { D }$ is diagonal. Repeat with $M _ { 1 }$ replaced by

$$
M _ {2} = \left[ \begin{array}{c c} 1 & \alpha_ {2} \\ \beta_ {2} & 1 \end{array} \right].
$$

(b) Show that either $\Vert \ M _ { 1 } ^ { T } D M _ { 1 } \ \Vert _ { 2 } \leq 2 \Vert \ D \ \Vert _ { 2 }$ or $| { \cal M } _ { 2 } ^ { T } { \cal D } { \cal M } _ { 2 } \| _ { 2 } \leq 2 \| { \cal D } \| _ { 2 }$ . (c) Suppose $\boldsymbol { x } \in \mathbb { R } ^ { m }$ and that $D \in \mathbb { R } ^ { n \times n }$ is diagonal with positive diagonal entries. Given indices i and j with $1 \leq i < j \leq m .$ , show how to compute $M \in \mathbb { R } ^ { n \times n }$ so that if $y = M x$ and $\tilde { D } = M ^ { T } D M$ , then $y _ { j } = 0$ and $\tilde { D }$ is diagonal with $\| \tilde { D } \| _ { 2 } \leq 2 \| D \| _ { 2 }$ . (d) From part (c) conclude that $Q = D ^ { 1 / 2 } M \tilde { D } ^ { - 1 / 2 }$ is orthogonal and that the update $y = M x$ can be diagonally transformed to $( D ^ { 1 / 2 } y ) = Q ( D ^ { 1 / 2 } x )$ .

# Notes and References for 5.1

Householder matrices are named after A.S. Householder, who popularized their use in numerical analysis. However, the properties of these matrices have been known for quite some time, see:

H.W. Turnbull and A.C. Aitken (1961). An Introduction to the Theory of Canonical Matrices, Dover Publications, New York, 102–105.

Other references concerned with Householder transformations include:

A.R. Gourlay (1970). “Generalization of Elementary Hermitian Matrices,” Comput. J. 13, 411–412.   
B.N. Parlett (1971). “Analysis of Algorithms for Reflections in Bisectors,” SIAM Review 13, 197–208.   
N.K. Tsao (1975). “A Note on Implementing the Householder Transformations.” SIAM J. Numer. Anal. 12, 53–58.   
B. Danloy (1976). “On the Choice of Signs for Householder Matrices,” J. Comput. Appl. Math. 2, 67–69.   
J.J.M. Cuppen (1984). “On Updating Triangular Products of Householder Matrices,” Numer. Math. 45, 403–410.   
A.A. Dubrulle (2000). “Householder Transformations Revisited,” SIAM J. Matrix Anal. Applic. 22, 33–40.   
J.W. Demmel, M. Hoemmen, Y. Hida, and E.J. Riedy (2009). “Nonnegative Diagonals and High Performance On Low-Profile Matrices from Householder QR,” SIAM J. Sci. Comput. 31, 2832– 2841.

A detailed error analysis of Householder transformations is given in Lawson and Hanson (SLE, pp. 83–89). The basic references for block Householder representations and the associated computations include:

C.H. Bischof and C. Van Loan (1987). “The WY Representation for Products of Householder Matrices,” SIAM J. Sci. Stat. Comput. 8, s2–s13.

B.N. Parlett and R. Schreiber (1988). “Block Reflectors: Theory and Computation,” SIAM J. Numer. Anal. 25, 189–205.   
R.S. Schreiber and C. Van Loan (1989). “A Storage-Efficient WY Representation for Products of Householder Transformations,” SIAM J. Sci. Stat. Comput. 10, 52–57.   
C. Puglisi (1992). “Modification of the Householder Method Based on the Compact WY Representation,” SIAM J. Sci. Stat. Comput. 13, 723–726.   
X. Sun and C.H. Bischof (1995). “A Basis-Kernel Representation of Orthogonal Matrices,” SIAM J. Matrix Anal. Applic. 16, 1184–1196.   
T. Joffrain, T.M. Low, E.S. Quintana-Orti, R. van de Geijn, and F.G. Van Zee (2006). “Accumulating Householder Transformations, Revisited,” ACM Trans. Math. Softw. 32, 169–179.   
M. Sadkane and A. Salam (2009). “A Note on Symplectic Block Reflectors,” ETNA 33, 45–52.   
Givens rotations are named after Wallace Givens. There are some subtleties associated with their computation and representation:   
G.W. Stewart (1976). “The Economical Storage of Plane Rotations,” Numer. Math. 25, 137–138.   
D. Bindel, J. Demmel, W. Kahan, and O. Marques (2002). “On computing givens rotations reliably and efficiently,” ACM Trans. Math. Softw. 28, 206–238.   
It is possible to aggregate rotation transformations to achieve high performance, see:   
B. Lang (1998). “Using Level 3 BLAS in Rotation–Based Algorithms,” SIAM J. Sci. Comput. 19, 626–634.   
Fast Givens transformations (see P5.1.11) are also referred to as square-root-free Givens transformations. (Recall that a square root must ordinarily be computed during the formation of Givens transformation.) There are several ways fast Givens calculations can be arranged, see:   
M. Gentleman (1973). “Least Squares Computations by Givens Transformations without Square Roots,” J. Inst. Math. Appl. 12, 329–336.   
C.F. Van Loan (1973). “Generalized Singular Values With Algorithms and Applications,” PhD Thesis, University of Michigan, Ann Arbor.   
S. Hammarling (1974). “A Note on Modifications to the Givens Plane Rotation,” J. Inst. Math. Applic. 13, 215–218.   
J.H. Wilkinson (1977). “Some Recent Advances in Numerical Linear Algebra,” in The State of the Art in Numerical Analysis, D.A.H. Jacobs (ed.), Academic Press, New York, 1–53.   
A.A. Anda and H. Park (1994). “Fast Plane Rotations with Dynamic Scaling,” SIAM J. Matrix Anal. Applic. 15, 162–174.   
R.J. Hanson and T. Hopkins (2004). “Algorithm 830: Another Visit with Standard and Modified Givens Transformations and a Remark on Algorithm 539,” ACM Trans. Math. Softw. 20, 86–94.
