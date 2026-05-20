# 6.4 Subspace Computations with the SVD

It is sometimes necessary to investigate the relationship between two given subspaces. How close are they? Do they intersect? Can one be “rotated” into the other? And so on. In this section we show how questions like these can be answered using the singular value decomposition.

# 6.4.1 Rotation of Subspaces

Suppose $A \in \mathbb { R } ^ { m \times p }$ is a data matrix obtained by performing a certain set of experiments. If the same set of experiments is performed again, then a different data matrix, $B \in \mathbb { R } ^ { m \times p }$ , is obtained. In the orthogonal Procrustes problem the possibility that B can be rotated into A is explored by solving the following problem:

$$
\text { minimize } \parallel A - B Q \parallel_ {F}, \quad \text { subject   to } Q ^ {T} Q = I _ {p}. \tag {6.4.1}
$$

We show that optimizing $Q$ can be specified in terms of the SVD of $B ^ { T } A$ . The matrix trace is critical to the derivation. The trace of a matrix is the sum of its diagonal entries:

$$
\operatorname{tr} (C) = \sum_ {i = 1} ^ {n} c _ {i i}, \quad C \in \mathbb {R} ^ {n \times n}.
$$

It is easy to show that if $C _ { 1 }$ and $C _ { 2 }$ have the same row and column dimension, then

$$
\operatorname{tr} (C _ {1} ^ {T} C _ {2}) = \operatorname{tr} (C _ {2} ^ {T} C _ {1})  . \tag {6.4.2}
$$

Returning to the Procrustes problem (6.4.1), if $Q \in \mathbb { R } ^ { p \times p }$ is orthogonal, then

$$
\begin{array}{l} \left\| A - B Q \right\| _ {F} ^ {2} = \sum_ {k = 1} ^ {p} \left\| A (:, k) - B \cdot Q (:, k) \right\| _ {2} ^ {2} \\ = \sum_ {k = 1} ^ {p} \| A (:, k) \| _ {2} ^ {2} + \| B Q (:, k) \| _ {2} ^ {2} - 2 Q (:, k) ^ {T} B ^ {T} A (:, k) \\ = \| A \| _ {F} ^ {2} + \| B Q \| _ {F} ^ {2} - 2 \sum_ {k = 1} ^ {p} \left[ Q ^ {T} \left(B ^ {T} A\right) \right] _ {k k} \\ = \| A \| _ {F} ^ {2} + \| B \| _ {F} ^ {2} - 2 \operatorname{tr} \left(Q ^ {T} \left(B ^ {T} A\right)\right). \\ \end{array}
$$

Thus, (6.4.1) is equivalent to the problem

$$
\max _ {Q ^ {T} Q = I _ {p}} \operatorname{tr} (Q ^ {T} B ^ {T} A)  .
$$

If $\begin{array} { c c c c c } { { U ^ { T } ( B ^ { T } A ) V } } & { { = } } & { { \Sigma } } & { { = } } & { { \mathrm { d i a g } ( \sigma _ { 1 } , \dots , \sigma _ { p } ) } } \end{array}$ is the SVD of $B ^ { T } A$ and we define the orthogonal matrix Z by $Z ~ = ~ V ^ { T } Q ^ { T } U$ , then by using (6.4.2) we have

$$
\mathsf {t r} (Q ^ {T} B ^ {T} A) = \mathsf {t r} (Q ^ {T} U \Sigma V ^ {T}) = \mathsf {t r} (Z \Sigma) = \sum_ {i = 1} ^ {p} z _ {i i} \sigma_ {i} \leq \sum_ {i = 1} ^ {p} \sigma_ {i}.
$$

The upper bound is clearly attained by setting $Z = I _ { p } , \mathrm { i . e . , } Q = U V ^ { T }$ .

Algorithm 6.4.1 Given A and B in $\mathbb { R } ^ { m \times p }$ , the following algorithm finds an orthogonal $Q \in \mathbb { R } ^ { p \times p }$ such that $\| A - B Q \| _ { F }$ is minimum.

$$
C = B ^ {T} A
$$

Compute the SVD $U ^ { T } C V = \Sigma$ and save U and V .

$$
Q = U V ^ {T}
$$

We mention that if $B = I _ { p } ,$ , then the problem (6.4.1) is related to the polar decomposition. This decomposition states that any square matrix A has a factorization of the form $A = Q P$ where Q is orthogonal and P is symmetric and positive semidefinite. Note that if $A = U \Sigma V ^ { T }$ is the SVD of A, then $\overset { \cdot } { A } = ( U V ^ { T } ) ( V \overset { \cdot } { \Sigma } V ^ { T } )$ is its polar decomposition. For further discussion, see §9.4.3.

# 6.4.2 Intersection of Nullspaces

Let $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { p \times n }$ be given, and consider the problem of finding an orthonormal basis for null $( A ) \cap \mathsf { n u l l } ( B )$ . One approach is to compute the nullspace of the matrix

$$
C = \left[ \begin{array}{l} A \\ B \end{array} \right]
$$

since this is just what we want: $C x = 0 \Leftrightarrow x \in \mathsf { n u l l } ( A ) \cap \mathsf { n u l l } ( B )$ . However, a more economical procedure results if we exploit the following theorem.

Theorem 6.4.1. Suppose $A \in \mathbb { R } ^ { m \times n }$ and let $\{ z _ { 1 } , \ldots , z _ { t } \}$ be an orthonormal basis for null(A). Define $Z ~ = ~ \left[ ~ z _ { 1 } ~ | \cdots | ~ z _ { t } ~ \right]$ and let $\{ w _ { 1 } , \ldots , w _ { q } \}$ be an orthonormal basis for null(BZ) where $B \in \mathbb { R } ^ { p \times n }$ . $I f W = [  w _ { 1 } | \cdot \cdot \cdot | w _ { q } ]$ , then the columns of ZW form an orthonormal basis for null(A) ∩ null(B).

Proof. Since $A Z = 0$ and $( B Z ) W = 0$ , we clearly have ran $( Z W ) \subset$ null(A) ∩ null(B). Now suppose x is in both null(A) and null(B). It follows that $x \ = \ Z a$ for some $0 \neq a \in \mathbb { R } ^ { t }$ . But since $0 = B x = B Z a$ , we must have $a = W b$ for some $b \in \mathbb { R } ^ { q }$ . Thus, $x = Z W b \in \mathsf { r a n } ( Z W )$

If the SVD is used to compute the orthonormal bases in this theorem, then we obtain the following procedure:

Algorithm 6.4.2 Given $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { p \times n }$ , the following algorithm computes and integer s and a matrix $Y = [ y _ { 1 } \vert \cdot \cdot \cdot \vert y _ { s } ]$ having orthonormal columns which span $\mathsf { n u l l } ( A ) \cap \mathsf { n u l l } ( B )$ . If the intersection is trivial, then $s = 0$ .

Compute the SVD $U _ { A } ^ { T } A V _ { A } = \mathrm { d i a g } ( \sigma _ { i } )$ , save $V _ { A }$ , and set $r = \mathsf { r a n k } ( A )$ .

if $r < n$

$$
C = B V _ {A} (:, r + 1: n)
$$

Compute the SVD $U _ { c } ^ { T } C V _ { c } = \mathrm { { d i a g } } ( \gamma _ { i } )$ , save $V _ { C }$ , and set $q = \mathsf { r a n k } ( C )$ .

if q < n − r

$$
s = n - r - q
$$

$$
Y = V _ {A} (:, r + 1: n) V _ {C} (:, q + 1: n - r)
$$

else

$$
s = 0
$$

end

else

$$
s = 0
$$

end

The practical implementation of this algorithm requires an ability to reason about numerical rank. See §5.4.1.

# 6.4.3 Angles Between Subspaces

Let F and G be subspaces in $\mathbb { R } ^ { m }$ whose dimensions satisfy

$$
p = \dim (F) \geq \dim (G) = q \geq 1.
$$

The principal angles $\{ \theta _ { i } \} _ { i = 1 } ^ { q }$ between these two subspaces and the associated principal vectors $\{ f _ { 1 } , g _ { i } \} _ { i = 1 } ^ { q }$ are defined recursively by

$$
\begin{array}{l} \cos (\theta_ {k}) = f _ {k} ^ {T} g _ {k} = \max \quad \max \quad f ^ {T} g. \\ f \in F, \| f \| _ {2} = 1 \quad g \in G, \| g \| _ {2} = 1 \tag {6.4.3} \\ \end{array}
$$

$$
f ^ {T} [ f _ {1}, \dots , f _ {k - 1} ] = 0 \quad g ^ {T} [ g _ {1}, \dots , g _ {k - 1} ] = 0
$$

Note that the principal angles satisfy $0 \le \theta _ { 1 } \le \dots \le \theta _ { q } \le \pi / 2$ .. The problem of computing principal angles and vectors is oftentimes referred to as the canonical correlation problem.

Typically, the subspaces F and G are matrix ranges, e.g.,

$$
F = \operatorname{ran} (A), \qquad A \in \mathbb {R} ^ {n \times p},
$$

$$
G = \operatorname{ran} (B), \qquad B \in \mathbb {R} ^ {n \times q}.
$$

The principal vectors and angles can be computed using the QR factorization and the SVD. Let $A = Q _ { A } R _ { A }$ and $B = Q _ { B } R _ { B }$ be thin QR factorizations and assume that

$$
Q _ {A} ^ {T} Q _ {B} = Y \Sigma Z ^ {T} = \sum_ {i = 1} ^ {q} \sigma_ {i} y _ {i} z _ {i} ^ {T}
$$

is the SVD of $Q _ { A } ^ { T } Q _ { B } \in \mathbb { R } ^ { p \times q }$ . Since $\Vert Q _ { A } ^ { T } Q _ { B } \parallel _ { 2 } \leq 1$ , all the singular values are between 0 and 1 and we may write $\sigma _ { i } = \cos ( \theta _ { i } ) , i = 1 { : } q$ . Let

$$
Q _ {A} Y = \left[ f _ {1} \mid \dots \mid f _ {p} \right], \tag {6.4.4}
$$

$$
Q _ {B} Z = \left[ g _ {1} \mid \dots \mid g _ {q} \right] \tag {6.4.5}
$$

be column partitionings of the matrices $Q _ { A } Y \in \mathbb { R } ^ { n \times p }$ and $Q _ { B } Z \in \mathbb { R } ^ { n \times q }$ . These matrices have orthonormal columns. If $f \in F$ and $g \in G$ are unit vectors, then there exist unit vectors $u \in \mathbb { R } ^ { p }$ and $v \in \mathbb { R } ^ { q }$ so that $f = Q _ { A } u$ and $g = Q _ { B } v$ . Thus,

$$
\begin{array}{l} f ^ {T} g = (Q _ {A} u) ^ {T} (Q _ {B} v) = u ^ {T} (Q _ {A} ^ {T} Q _ {B}) v = u ^ {T} (Y \Sigma Z ^ {T}) v \\ = (Y ^ {T} u) ^ {T} \Sigma (Z ^ {T} v) = \sum_ {i = 1} ^ {q} \sigma_ {i} (y _ {i} ^ {T} u) (z _ {i} ^ {T} v). \tag {6.4.6} \\ \end{array}
$$

This expression attains its maximal value of $\sigma _ { 1 } = \cos ( \theta _ { 1 } )$ by setting $u = y _ { 1 }$ and $v = z _ { 1 }$ . It follows that $f = Q _ { A } y _ { 1 } = f _ { 1 }$ and $v = Q _ { B } z _ { 1 } = g _ { 1 }$ .

Now assume that k > 1 and that the first k − 1 columns of the matrices in (6.4.4) and (6.4.5) are known, $\mathrm { i . e . , ~ } f _ { 1 } , \ldots , f _ { k - 1 }$ and $g _ { 1 } , \ldots , g _ { k - 1 }$ . Consider the problem of maximizing $f ^ { T } g$ given that $f = Q _ { A } u$ and $g = Q _ { B } v$ are unit vectors that satisfy

$$
f ^ {T} \left[ f _ {1} \mid \dots \mid f _ {k - 1} \right] = 0,
$$

$$
g ^ {T} \left[ g _ {1} \mid \dots \mid g _ {k - 1} \right] = 0.
$$

It follows from (6.4.6) that

$$
f ^ {T} g = \sum_ {i = k} ^ {q} \sigma_ {i} (y _ {i} ^ {T} u) (z _ {i} ^ {T} v) \leq \sigma_ {k} \sum_ {i = k} ^ {q} | y _ {i} ^ {T} u | \cdot | z _ {i} ^ {T} v |.
$$

This expression attains its maximal value of $\sigma _ { k } = \cos ( \theta _ { k } )$ by setting $u = y _ { k }$ and $v = z _ { k }$ . It follows from (6.4.4) and (6.4.5) that $f = Q _ { A } y _ { k } = f _ { k }$ and $g = Q _ { B } z _ { k } = g _ { k }$ . Combining these observations we obtain

Algorithm 6.4.3 (Principal Angles and Vectors) Given $A \in \mathbb { R } ^ { m \times p }$ and $B \in \mathbb { R } ^ { m \times q }$ $( p \geq q )$ each with linearly independent columns, the following algorithm computes the cosines of the principal angles $\theta _ { 1 } \geq \cdots \geq \theta _ { q }$ between ran(A) and ran(B). The vectors $f _ { 1 } , \ldots , f _ { q }$ and $g _ { 1 } , \ldots , g _ { q }$ are the associated principal vectors.

Compute the thin QR factorizations $A = Q _ { A } R _ { A }$ and $B = Q _ { B } R _ { B }$ .

$$
C = Q _ {A} ^ {T} Q _ {B}
$$

Compute the SVD $Y ^ { T } C Z = \mathrm { d i a g } ( \cos ( \theta _ { k } ) )$ .

$$
Q _ {A} Y (:, 1: q) = \left[ f _ {1} \mid \dots \mid f _ {q} \right]
$$

$$
Q _ {B} Z (:, 1: q) = \left[ g _ {1} \mid \dots \mid g _ {q} \right]
$$

The idea of using the SVD to compute the principal angles and vectors is due to Bj¨orck and Golub (1973). The problem of rank deficiency in A and B is also treated in this paper. Principal angles and vectors arise in many important statistical applications. The largest principal angle is related to the notion of distance between equidimensional subspaces that we discussed in §2.5.3. If $p = q$ , then

$$
\operatorname{dist} (F, G) = \sqrt {1 - \cos (\theta_ {p}) ^ {2}} = \sin (\theta_ {p}).
$$

# 6.4.4 Intersection of Subspaces

In light of the following theorem, Algorithm 6.4.3 can also be used to compute an orthonormal basis for ran(A) ∩ ran(B) where $A \in \mathbb { R } ^ { m \times p }$ and $B \in \mathbb { R } ^ { m \times q }$

Theorem 6.4.2. Let $\{ \cos ( \theta _ { i } ) \} _ { i = 1 } ^ { q }$ and $\{ f _ { i } , g _ { i } \} _ { i = 1 } ^ { q }$ be defined by Algorithm 6.4.3. If the index s is defined by $1 = \cos ( \theta _ { 1 } ) = \cdot \cdot \cdot = \cos ( \theta _ { s } ) > \cos ( \theta _ { s + 1 } )$ , then

$$
\operatorname{ran} (A) \cap \operatorname{ran} (B) = \operatorname{span} \left\{f _ {1}, \dots , f _ {s} \right\} = \operatorname{span} \left\{g _ {1}, \dots , g _ {s} \right\}.
$$

Proof. The proof follows from the observation that if cos $( \theta _ { i } ) = 1$ , then $f _ { i } = g _ { i }$ .

The practical determination of the intersection dimension s requires a definition of what it means for a computed singular value to equal 1. For example, a computed singular value $\hat { \sigma } _ { i } = \cos ( \hat { \theta } _ { i } )$ could be regarded as a unit singular value if $\hat { \sigma } _ { i } \geq 1 - \delta$ for some intelligently chosen small parameter δ.

# Problems

P6.4.1 Show that if A and B are m-by-p matrices, with $p \leq m$ , then

$$
\min _ {Q ^ {T} Q = I _ {p}} \| A - B Q \| _ {F} ^ {2} = \sum_ {i = 1} ^ {p} (\sigma_ {i} (A) ^ {2} - 2 \sigma_ {i} (B ^ {T} A) + \sigma_ {i} (B) ^ {2}).
$$

P6.4.2 Extend Algorithm 6.4.2 so that it computes an orthonormal basis for $\mathsf { n u l l } ( A _ { 1 } ) \cap \cdots \cap \mathsf { n u l l } ( A _ { s } )$ where each matrix Ai has n columns.

P6.4.3 Extend Algorithm 6.4.3 so that it can handle the case when A and B are rank deficient.

P6.4.4 Verify Equation (6.4.2).

P6.4.5 Suppose A, $B \in \mathbb { R } ^ { m \times n }$ and that A has full column rank. Show how to compute a symmetric matrix $X \in \mathbb { R } ^ { n \times n }$ that minimizes $\| A X - B \| _ { F }$ . Hint: Compute the SVD of A.

P6.4.6 This problem is an exercise in F-norm optimization. (a) Show that if $C \in \mathbb { R } ^ { m \times n }$ and $e \in \mathbb { R } ^ { m }$ is a vector of ones, then $v = C ^ { T } e / m$ minimizes $\| \boldsymbol { C } - e \boldsymbol { v } ^ { T } \| _ { F }$ . (b) Suppose $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { m \times n }$ and that we wish to solve

$$
\min _ {Q ^ {T} Q = I _ {n}, v \in \mathbf {R} ^ {n}} \| A - (B + e v ^ {T}) Q \| _ {F}
$$

Show that $\boldsymbol { v } _ { \mathrm { o p t } } = ( A - B ) ^ { T } \boldsymbol { e } / m$ and $Q _ { \mathrm { o p t } } = U \Sigma V ^ { T }$ solve this problem where $B ^ { T } ( I - e e ^ { T } / m ) A = U V ^ { T }$ is the SVD.

P6.4.7 A 3-by-3 matrix H is ROPR matrix if $H = Q + x y ^ { T }$ where $Q \in \mathbb { R } ^ { 3 \times 3 }$ rotation and $x , y \in \mathbb { R } ^ { 3 }$ . (A rotation matrix is an orthogonal matrix with unit determinant. $\mathrm { \mathrm { \Omega ^ { * } R O P R } \mathrm { \Omega ^ { * } } }$ stands for “rank-1 perturbation of a rotation.”) ROPR matrices arise in computational photography and this problem highlights some of their properties. (a) If H is a ROPR matrix, then there exist rotations $U , V \in \mathbb { R } ^ { 3 \times 3 }$ , such that $\begin{array} { r } { U ^ { T } H V = \mathrm { d i a g } ( \sigma _ { 1 } , \sigma _ { 2 } , \sigma _ { 3 } ) } \end{array}$ satisfies $\sigma _ { 1 } \geq \sigma _ { 2 } \geq | \sigma _ { 3 } |$ . (b) Show that if $Q \in \mathbb { R } ^ { 3 \times 3 }$ is a rotation, then there exist cosine-sine pairs $( c _ { i } , s _ { i } ) = ( \cos ( \theta _ { i } ) , \sin ( \theta _ { i } ) ) , i = 1 { : } 3$ such that $Q = Q ( \theta _ { 1 } , \theta _ { 2 } , \theta _ { 3 } )$ where

$$
\begin{array}{l} Q (\theta_ {1}, \theta_ {2}, \theta_ {3}) = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & c _ {1} & s _ {1} \\ 0 & - s _ {1} & c _ {1} \end{array} \right] \left[ \begin{array}{c c c} c _ {2} & s _ {2} & 0 \\ - s _ {2} & c _ {2} & 0 \\ 0 & 0 & 1 \end{array} \right] \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & c _ {3} & s _ {3} \\ 0 & - s _ {3} & c _ {3} \end{array} \right] \\ = \left[ \begin{array}{c c c} c _ {2} & s _ {2} c _ {3} & s _ {2} s _ {3} \\ - c _ {1} s _ {2} & c _ {1} c _ {2} c _ {3} - s _ {1} s _ {3} & c _ {1} c _ {2} s _ {3} + s _ {1} c _ {3} \\ s _ {1} s _ {2} & - s _ {1} c _ {2} c _ {3} - c _ {1} s _ {3} & - s _ {1} c _ {2} s _ {3} + c _ {1} c _ {3} \end{array} \right]. \\ \end{array}
$$

Hint: The Givens QR factorization involves three rotations. (c) Show that if

$$
\left[ \begin{array}{c c c} \sigma_ {1} & 0 & 0 \\ 0 & \sigma_ {2} & 0 \\ 0 & 0 & \sigma_ {3} \end{array} \right] = Q (\theta_ {1}, \theta_ {2}, \theta_ {3}) - x y ^ {T}, \qquad x, y \in \mathbb {R} ^ {3}
$$

then $x y ^ { T }$ must have the form

$$
x y ^ {T} = \left[ \begin{array}{c} s _ {2} \\ \mu c _ {1} \\ - \mu s _ {1} \end{array} \right] \left[ \begin{array}{c} - s _ {2} / \mu \\ c _ {3} \\ s _ {3} \end{array} \right] ^ {T}
$$

for some $\mu \geq 0$ and

$$
\left[ \begin{array}{c c} c _ {2} - \mu & 1 \\ 1 & c _ {2} - \mu \end{array} \right] \left[ \begin{array}{c} c _ {1} s _ {3} \\ s _ {1} c _ {3} \end{array} \right] = \left[ \begin{array}{c} 0 \\ 0 \end{array} \right].
$$

(d) Show that the second singular value of a ROPR matrix is 1.

P6.4.8 Let $U _ { * } \in \mathbb { R } ^ { n \times d }$ be a matrix with orthonormal columns whose span is a subspace S that we wish to estimate. Assume that $U _ { c } \in \mathbb { R } ^ { n \times d }$ is a given matrix with orthonormal columns and regard $\mathsf { r a n } ( U _ { c } )$ as the “current” estimate of $S _ { \cdot }$ . This problem examines what is required to get an improved estimate of S given the availability of a vector $v \in S$ . (a) Define the vectors

$$
w = U _ {c} ^ {T} v, \qquad v _ {1} = U _ {c} U _ {c} ^ {T} v, \qquad v _ {2} = (I _ {n} - U _ {c} U _ {c} ^ {T}) v,
$$

and assume that each is nonzero. (a) Show that if

$$
z _ {\theta} = \left(\frac {\cos (\theta) - 1}{\| v _ {1} \| \| w \|}\right) v _ {1} + \left(\frac {\sin (\theta)}{\| v _ {2} \| \| w \|}\right) v _ {2}
$$

and

$$
U _ {\theta} = (I _ {n} + z _ {\theta} v ^ {T}) U _ {c},
$$

then $U _ { \theta } ^ { T } U _ { \theta } = I _ { d }$ . Thus, $U _ { \theta } U _ { \theta } ^ { T }$ is an orthogonal projection. (b) Define the distance function

$$
\operatorname{dist} _ {F} (\operatorname{ran} (V), \operatorname{ran} (W)) = \| V V ^ {T} - W W ^ {T} \| _ {F}
$$

where $V , W \in \mathbb { R } ^ { n \times d }$ have orthonormal columns and show

$$
\mathsf {d i s t} _ {F} (\mathsf {r a n} (V), \mathsf {r a n} (W)) ^ {2} = 2 (d - \| W ^ {T} V \| _ {F} ^ {2}) = 2 \sum_ {i = 1} ^ {d} (1 - \sigma_ {i} (W ^ {T} V) ^ {2}).
$$

$\operatorname { N o t e \ t h a t \ d i s t } ( \mathsf { r a n } ( V ) , \mathsf { r a n } ( W ) ) ^ { 2 } = 1 - \sigma _ { 1 } ( W ^ { T } V ) ^ { 2 } . ~ ( \mathsf { c } ) \operatorname { S h o w \ t h a t }$

$$
d _ {\theta} ^ {2} = d _ {c} ^ {2} - 2 \cdot \mathbf {t r} (U _ {*} U _ {*} ^ {T} (U _ {\theta} U _ {\theta} ^ {T} - U _ {c} U _ {c} ^ {T}))
$$

where $d _ { \theta } = \mathsf { d i s t } _ { F } ( \mathsf { r a n } ( U _ { * } ) , \mathsf { r a n } ( U _ { \theta } ) )$ and $d _ { c } = { \tt d i s t } _ { F } ( { \tt r a n } ( U _ { * } ) , { \tt r a n } ( U _ { c } ) )$ . (d) Show that if

$$
y _ {\theta} = \cos (\theta) \frac {v _ {1}}{\parallel v _ {1} \parallel} + \sin (\theta) \frac {v _ {2}}{\parallel v _ {2} \parallel},
$$

then

$$
U _ {\theta} U _ {\theta} ^ {T} - U _ {c} U _ {c} ^ {T} = y _ {\theta} y _ {\theta} ^ {T} - \frac {v _ {1} v _ {1} ^ {T}}{v _ {1} ^ {T} v _ {1}}
$$

and

$$
d _ {\theta} ^ {2} = d _ {c} ^ {2} + 2 \left(\frac {\| U _ {*} ^ {T} v _ {1} \| _ {2} ^ {2}}{\| v _ {1} \| _ {2} ^ {2}} - \| U _ {*} ^ {T} y _ {\theta} \| _ {2} ^ {2}\right).
$$

(e) Show that if θ minimizes this quantity, then

$$
\sin (2 \theta) \left(\frac {\| P _ {S} v _ {2} \| ^ {2}}{\| v _ {2} \| _ {2} ^ {2}} - \frac {\| P _ {S} v _ {1} \| ^ {2}}{\| v _ {1} \| _ {2} ^ {2}}\right) + \cos (2 \theta) \frac {v _ {1} ^ {T} P _ {S} v _ {2}}{\| v _ {1} \| _ {2} \| v _ {2} \| _ {2}} = 0, \qquad P _ {S} = U _ {*} U _ {*} ^ {T}.
$$

# Notes and References for 6.4

References for the Procrustes problem include:

Using the SVD to solve the angles-between-subspaces problem is discussed in:   
B. Green (1952). “The Orthogonal Approximation of an Oblique Structure in Factor Analysis,” Psychometrika 17, 429–40.   
P. Schonemann (1966). “A Generalized Solution of the Orthogonal Procrustes Problem,” Psychometrika 31, 1–10.   
R.J. Hanson and M.J. Norris (1981). “Analysis of Measurements Based on the Singular Value Decomposition,” SIAM J. Sci. Stat. Comput. 2, 363–374.   
N.J. Higham (1988). “The Symmetric Procrustes Problem,” BIT 28, 133–43.   
H. Park (1991). “A Parallel Algorithm for the Unbalanced Orthogonal Procrustes Problem,” Parallel Comput. 17, 913–923.   
L.E. Andersson and T. Elfving (1997). “A Constrained Procrustes Problem,” SIAM J. Matrix Anal. Applic. 18, 124–139.   
L. Eld´en and H. Park (1999). “A Procrustes Problem on the Stiefel Manifold,” Numer. Math. 82, 599–619.   
A.W. Bojanczyk and A. Lutoborski (1999). “The Procrustes Problem for Orthogonal Stiefel Matrices,” SIAM J. Sci. Comput. 21, 1291–1304.   
If $B \ = \ I ,$ , then the Procrustes problem amounts to finding the closest orthogonal matrix. This computation is related to the polar decomposition problem that we consider in §9.4.3. Here are some basic references:   
˚A. Bj¨orck and C. Bowie (1971). “An Iterative Algorithm for Computing the Best Estimate of an Orthogonal Matrix,” SIAM J. Numer. Anal. 8, 358–64.   
N.J. Higham (1986). “Computing the Polar Decomposition with Applications,” SIAM J. Sci. Stat. Comput. 7, 1160–1174.   
˚A. Bj¨orck and G.H. Golub (1973). “Numerical Methods for Computing Angles Between Linear Subspaces,” Math. Comput. 27, 579–94.   
L.M. Ewerbring and F.T. Luk (1989). “Canonical Correlations and Generalized SVD: Applications and New Algorithms,” J. Comput. Appl. Math. 27, 37–52.   
G.H. Golub and H. Zha (1994). “Perturbation Analysis of the Canonical Correlations of Matrix Pairs,” Lin. Alg. Applic. 210, 3–28.

Z. Drmac (2000). “On Principal Angles between Subspaces of Euclidean Space,” SIAM J. Matrix Anal. Applic. 22, 173–194.   
A.V. Knyazev and M.E. Argentati (2002). “Principal Angles between Subspaces in an A–Based Scalar Product: Algorithms and Perturbation Estimates,” SIAM J. Sci. Comput. 23, 2008–2040.   
P. Strobach (2008). “Updating the Principal Angle Decomposition,” Numer. Math. 110, 83–112.   
In reduced-rank regression the object is to connect a matrix of signals to a matrix of noisey observations through a matrix that has specified low rank. An svd-based computational procedure that involves principal angles is discussed in:   
L. Eld´en and B. Savas (2005). “The Maximum Likelihood Estimate in Reduced-Rank Regression,” Num. Lin. Alg. Applic. 12, 731–741,   
The SVD has many roles to play in statistical computation, see:   
S.J. Hammarling (1985). “The Singular Value Decomposition in Multivariate Statistics,” ACM SIGNUM Newsletter 20, 2-25.   
An algorithm for computing the rotation and rank-one matrix in P6.4.7 that define a given ROPR matrix is discussed in:   
R. Schreiber, Z. Li, and H. Baker (2009). “Robust Software for Computing Camera Motion Parameters,” J. Math. Imaging Vision 33, 1–9.   
For a more details about the estimation problem associated with P6.4.8, see:   
L. Balzano, R. Nowak, and B. Recht (2010). “Online Identification and Tracking of Subspaces from Highly Incomplete Information,” Proceedings of the Allerton Conference on Communication, Control, and Computing 2010.
