# 2.4 The Singular Value Decomposition

It is fitting that the first matrix decomposition that we present in the book is the singular value decomposition (SVD). The practical and theoretical importance of the SVD is hard to overestimate. It has a prominent role to play in data analysis and in the characterization of the many matrix “nearness problems.”

# 2.4.1 Derivation

The SVD is an orthogonal matrix reduction and so the 2-norm and Frobenius norm figure heavily in this section. Indeed, we can prove the existence of the decomposition using some elementary facts about the 2-norm developed in the previous two sections.

Theorem 2.4.1 (Singular Value Decomposition ). If A is a real m-by-n matrix, then there exist orthogonal matrices

$$
U = \left[ u _ {1} \mid \dots \mid u _ {m} \right] \in \mathbb {R} ^ {m \times m} \quad a n d \quad V = \left[ v _ {1} \mid \dots \mid v _ {n} \right] \in \mathbb {R} ^ {n \times n}
$$

such that

$$
U ^ {T} A V = \Sigma = \operatorname{diag} \left(\sigma_ {1}, \dots , \sigma_ {p}\right) \in \mathbb {R} ^ {m \times n}, \quad p = \min \{m, n \},
$$

where $\sigma _ { 1 } \geq \sigma _ { 2 } \geq . . . \geq \sigma _ { p } \geq 0 .$

Proof. Let $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $\boldsymbol { y } \in \mathbb { R } ^ { m }$ be unit 2-norm vectors that satisfy $A x = \sigma y$ with $\sigma = \parallel A \parallel _ { 2 }$ . From Theorem 2.1.1 there exist $V _ { 2 } \in \mathbb { R } ^ { n \times ( n - 1 ) }$ and $\bar { U _ { 2 } } \in \mathbb { R } ^ { m \times ( m - 1 ) }$ so $V = \left[ \ b { x } \mid V _ { 2 } \right] \in \mathbb { R } ^ { n \times n }$ and $U = \left[ \left. y \right| U _ { 2 } \right] \in \mathbb { R } ^ { m \times m }$ are orthogonal. It is not hard to show that

$$
U ^ {T} A V = \left[ \begin{array}{c c} \sigma & w ^ {T} \\ 0 & B \end{array} \right] \equiv A _ {1}
$$

where $w \in \mathbb { R } ^ { n - 1 }$ and $B \in \mathbb { R } ^ { ( m - 1 ) \times ( n - 1 ) }$ . Since

$$
\left\| A _ {1} \left(\left[ \begin{array}{c} \sigma \\ w \end{array} \right]\right) \right\| _ {2} ^ {2} \geq (\sigma^ {2} + w ^ {T} w) ^ {2}
$$

we have $\parallel A _ { 1 } \parallel _ { 2 } ^ { 2 } \geq ( \sigma ^ { 2 } + w ^ { T } w )$ . But $\sigma ^ { 2 } = \parallel A \parallel _ { 2 } ^ { 2 } = \parallel A _ { 1 } \parallel _ { 2 } ^ { 2 }$ , and so we must have $w = 0$ . An obvious induction argument completes the proof of the theorem.

The $\sigma _ { i }$ are the singular values of A, the $u _ { i }$ are the left singular vectors of A, and the $v _ { i }$ are right singular vectors of A. Separate visualizations of the SVD are required depending upon whether A has more rows or columns. Here are the 3-by-2 and 2-by-3 examples:

$$
\left[ \begin{array}{l l l} u _ {1 1} & u _ {1 2} & u _ {1 3} \\ u _ {2 1} & u _ {2 2} & u _ {2 3} \\ u _ {3 1} & u _ {3 2} & u _ {3 3} \end{array} \right] ^ {T} \left[ \begin{array}{l l} a _ {1 1} & a _ {1 2} \\ a _ {2 1} & a _ {2 2} \\ a _ {3 1} & a _ {3 2} \end{array} \right] \left[ \begin{array}{l l} v _ {1 1} & v _ {1 2} \\ v _ {2 1} & v _ {2 2} \end{array} \right] = \left[ \begin{array}{l l} \sigma_ {1} & 0 \\ 0 & \sigma_ {2} \\ 0 & 0 \end{array} \right],
$$

$$
\left[ \begin{array}{c c} u _ {1 1} & u _ {1 2} \\ u _ {2 1} & u _ {2 2} \end{array} \right] ^ {T} \left[ \begin{array}{c c c} a _ {1 1} & a _ {1 2} & a _ {1 3} \\ a _ {2 1} & a _ {2 2} & a _ {2 3} \end{array} \right] \left[ \begin{array}{c c c} v _ {1 1} & v _ {1 2} & v _ {1 3} \\ v _ {2 1} & v _ {2 2} & v _ {2 3} \\ v _ {3 1} & v _ {3 2} & v _ {3 3} \end{array} \right] = \left[ \begin{array}{c c c} \sigma_ {1} & 0 & 0 \\ 0 & \sigma_ {2} & 0 \end{array} \right].
$$

In later chapters, the notation $\sigma _ { i } ( A )$ is used to designate the ith largest singular value of a matrix A. The largest and smallest singular values are important and for them we also have a special notation:

$$
\sigma_ {\max} (A) = \text { the   largest   singular   value   of   matrix } A,
$$

$$
\sigma_ {\min} (A) = \text { the   smallest   singular   value   of   matrix } A.
$$

# 2.4.2 Properties

We establish a number of important corollaries to the SVD that are used throughout the book.

Corollary 2.4.2. $I f U ^ { T } A V = \Sigma$ is the SVD of $A \in \mathbb { R } ^ { m \times n }$ and $m \geq n$ , then $f o r i = 1$ :n $A v _ { i } = \sigma _ { i } u _ { i }$ and $A ^ { T } u _ { i } = \sigma _ { i } v _ { i }$ .

Proof. Compare columns in $A V = U \Sigma$ and $A ^ { T } U = V \Sigma ^ { T }$ .□

There is a nice geometry behind this result. The singular values of a matrix A are the lengths of the semiaxes of the hyperellipsoid E defined by $E = \left\{ \ A x : \| \ x \ \| _ { 2 } = 1 \right\}$ . The semiaxis directions are defined by the $u _ { i }$ and their lengths are the singular values.

It follows immediately from the corollary that

$$
A ^ {T} A v _ {i} = \sigma_ {i} ^ {2} v _ {i}, \tag {2.4.1}
$$

$$
A A ^ {T} u _ {i} = \sigma_ {i} ^ {2} u _ {i} \tag {2.4.2}
$$

for $i = 1 { : } n$ . This shows that there is an intimate connection between the SVD of A and the eigensystems of the symmetric matrices $A ^ { T } A$ and $A A ^ { T }$ . See §8.6 and §10.4.

The 2-norm and the Frobenius norm have simple SVD characterizations.

Corollary 2.4.3. If $A \in \mathbb { R } ^ { m \times n }$ , then

$$
\| A \| _ {2} = \sigma_ {1}, \qquad \| A \| _ {F} = \sqrt {\sigma_ {1} ^ {2} + \cdots + \sigma_ {p} ^ {2}},
$$

where $p = \operatorname* { m i n } \{ m , n \}$ .

Proof. These results follow immediately from the fact that $\| U ^ { T } A V \| = \| \Sigma \|$ for both the 2-norm and the Frobenius norm. □

We show in §8.6 that if A is perturbed by a matrix E, then no singular value can move by more than $\parallel E \parallel _ { 2 }$ . The following corollary identifies two useful instances of this result.

Corollary 2.4.4. If $A \in \mathbb { R } ^ { m \times n }$ and $E \in \mathbb { R } ^ { m \times n }$ , then

$$
\sigma_ {\max} (A + E) \leq \sigma_ {\max} (A) + \| E \| _ {2},
$$

$$
\sigma_ {\min} (A + E) \geq \sigma_ {\min} (A) - \parallel E \parallel_ {2}.
$$

Proof. Using Corollary 2.4.2 it is easy to show that

$$
\sigma_ {\min} (A) \cdot \| x \| _ {2} \leq \| A x \| _ {2} \leq \sigma_ {\max} (A) \cdot \| x \| _ {2}.
$$

The required inequalities follow from this result.

If a column is added to a matrix, then the largest singular value increases and the smallest singular value decreases.

Corollary 2.4.5. If $A \in \mathbb { R } ^ { m \times n }$ , $m > n$ , and $z \in \mathbb { R } ^ { m }$ , then

$$
\sigma_ {\max} \big (\left[ A \mid z \right] \big) \geq \sigma_ {\max} (A),
$$

$$
\sigma_ {\min} \left(\left[ A \mid z \right]\right) \leq \sigma_ {\min} (A).
$$

Proof. Suppose $A = U \Sigma V ^ { T }$ is the SVD of A and let $x = V ( : , 1 )$ and $\tilde { A } = [ A | z ]$ . Using Corollary 2.4.4, we have

$$
\sigma_ {\max} (A) = \| A x \| _ {2} = \left\| \tilde {A} \left[ \begin{array}{c} x \\ 0 \end{array} \right] \right\| _ {2} \leq \sigma_ {\max} (\tilde {A}).
$$

The proof that $\sigma _ { \operatorname* { m i n } } ( A ) \geq \sigma _ { \operatorname* { m i n } } ( \tilde { A } )$ is similar.

The SVD neatly characterizes the rank of a matrix and orthonormal bases for both its nullspace and its range.

Corollary 2.4.6. If A has r positive singular values, then rank(A) = r and

$$
\operatorname{null} (A) = \operatorname{span} \left\{v _ {r + 1}, \dots , v _ {n} \right\},
$$

$$
\operatorname{ran} (A) = \operatorname{span} \{u _ {1}, \dots , u _ {r} \}.
$$

Proof. The rank of a diagonal matrix equals the number of nonzero diagonal entries. Thus, rank $( A ) = { \mathsf { r a n k } } ( \Sigma ) = r$ . The assertions about the nullspace and range follow from Corollary 2.4.2.

If A has rank r, then it can be written as the sum of r rank-1 matrices. The SVD gives us a particularly nice choice for this expansion.

Corollary 2.4.7. If $A \in \mathbb { R } ^ { m \times n }$ and rank(A) = r, then

$$
A = \sum_ {i = 1} ^ {r} \sigma_ {i} u _ {i} v _ {i} ^ {T}.
$$

Proof. This is an exercise in partitioned matrix multiplication:

$$
(U \Sigma) V ^ {T} = \left(\left[ \sigma_ {1} u _ {1} \mid \sigma_ {2} u _ {2} \mid \dots \mid \sigma_ {r} u _ {r} \mid 0 \mid \dots \mid 0 \right]\right) \left[ \begin{array}{c} v _ {1} ^ {T} \\ \vdots \\ v _ {n} ^ {T} \end{array} \right] = \sum_ {i = 1} ^ {r} \sigma_ {i} u _ {i} v _ {i} ^ {T}. \quad \square
$$

The intelligent handling of rank degeneracy is an important topic that we discuss in Chapter 5. The SVD has a critical role to play because it can be used to identify nearby matrices of lesser rank.

Theorem 2.4.8 (The Eckhart-Young Theorem). If $k < r = \mathrm { r a n k } ( A )$ and

$$
A _ {k} = \sum_ {i = 1} ^ {k} \sigma_ {i} u _ {i} v _ {i} ^ {T}, \tag {2.4.3}
$$

then

$$
\min _ {\operatorname{rank} (B) = k} \| A - B \| _ {2} = \| A - A _ {k} \| _ {2} = \sigma_ {k + 1}. \tag {2.4.4}
$$

Proof. Since $U ^ { T } A _ { k } V = \mathrm { d i a g } ( \sigma _ { 1 } , \dots , \sigma _ { k } , 0 , \dots , 0 )$ it follows that $A _ { k }$ is rank k. Moreover, $U ^ { T } ( A - A _ { k } ) V = \mathrm { d i a g } ( 0 , \dots , 0 , \sigma _ { k + 1 } , \dots , \sigma _ { p } )$ and so $\parallel A - A _ { k } \parallel _ { 2 } = \sigma _ { k + 1 }$ .

Now suppose $\mathsf { r a n k } ( B ) = k$ for some $\boldsymbol { B } \in \mathbb { R } ^ { m \times n }$ . It follows that we can find orthonormal vectors $x _ { 1 } , \ldots , x _ { n - k }$ so nul $\left| \left( B \right) = { \mathsf { s p a n } } \{ x _ { 1 } , \ldots , x _ { n - k } \} \right.$ . A dimension argument shows that

$$
\operatorname{span} \left\{x _ {1}, \dots , x _ {n - k} \right\} \cap \operatorname{span} \left\{v _ {1}, \dots , v _ {k + 1} \right\} \neq \{0 \}.
$$

Let $z$ be a unit 2-norm vector in this intersection. Since $B z = 0$ and

$$
A z = \sum_ {i = 1} ^ {k + 1} \sigma_ {i} (v _ {i} ^ {T} z) u _ {i},
$$

we have

$$
\parallel A - B \parallel_ {2} ^ {2} \geq \parallel (A - B) z \parallel_ {2} ^ {2} = \parallel A z \parallel_ {2} ^ {2} = \sum_ {i = 1} ^ {k + 1} \sigma_ {i} ^ {2} (v _ {i} ^ {T} z) ^ {2} \geq \sigma_ {k + 1} ^ {2},
$$

completing the proof of the theorem.

Note that this theorem says that the smallest singular value of A is the 2-norm distance of A to the set of all rank-deficient matrices. We also mention that the matrix $A _ { k }$ defined in (2.4.3) is the closest rank-k matrix to A in the Frobenius norm.

# 2.4.3 The Thin SVD

If $A = U \Sigma V ^ { T } \in \mathbb { R } ^ { m \times n }$ is the SVD of A and $m \geq n$ , then

$$
A = U _ {1} \Sigma_ {1} V ^ {T}
$$

where

$$
U _ {1} = U (:, 1: n) = \left[ u _ {1} \mid \dots \mid u _ {n} \right] \in \mathbb {R} ^ {m \times n}
$$

and

$$
\Sigma_ {1} = \Sigma (1: n, 1: n) = \operatorname{diag} \left(\sigma_ {1}, \dots , \sigma_ {n}\right) \in \mathbb {R} ^ {n \times n}.
$$

We refer to this abbreviated version of the SVD as the thin SVD.

# 2.4.4 Unitary Matrices and the Complex SVD

Over the complex field the unitary matrices correspond to the orthogonal matrices. In particular, $Q \in \mathbb { C } ^ { n \times n }$ is unitary if $Q ^ { H } Q = Q Q ^ { \bar { H } } = I _ { n }$ . Unitary transformations preserve both the 2-norm and the Frobenius norm. The SVD of a complex matrix involves unitary matrices. If $A \in \mathbb { C } ^ { m \times n }$ , then there exist unitary matrices $U \in \mathbb { C } ^ { m \times m }$ and $V \in \mathbb { C } ^ { n \times n }$ such that

$$
U ^ {H} A V = \mathrm{diag} (\sigma_ {1}, \ldots , \sigma_ {p}) \in \mathbb {R} ^ {m \times n} \qquad p = \min \{m, n \}
$$

where $\sigma _ { 1 } \geq \sigma _ { 2 } \geq . . . \geq \sigma _ { p } \geq 0$ . All of the real SVD properties given above have obvious complex analogs.

# Problems

P2.4.1 Show that if $Q = Q _ { 1 } + i Q _ { 2 }$ is unitary with $Q _ { 1 } , Q _ { 2 } \in \mathbb { R } ^ { n \times n }$ , then the 2n-by-2n real matrix

$$
Z = \left[ \begin{array}{c c} Q _ {1} & - Q _ {2} \\ Q _ {2} & Q _ {1} \end{array} \right]
$$

is orthogonal.

P2.4.2 Prove that if $A \in \mathbb { R } ^ { m \times n }$ , then

$$
\sigma_ {\max} (A) = \max _ { \begin{array}{c} y \in \mathbb {R} ^ {m} \\ x \in \mathbb {R} ^ {n} \end{array} } \frac {y ^ {T} A x}{\| x \| _ {2} \| y \| _ {2}}.
$$

P2.4.3 For the 2-by-2 matrix $A \ = \ \left[ \begin{array} { c c } { { w } } & { { x } } \\ { { y } } & { { z } } \end{array} \right]$ , derive expressions for $\sigma _ { \operatorname* { m a x } } ( A )$ and $\sigma _ { \mathrm { m i n } } ( A )$ that are functions of $w , x , y ,$ and z.

P2.4.4 Show that any matrix in $\mathbb { R } ^ { m \times n }$ is the limit of a sequence of full rank matrices.

P2.4.5 Show that if $A \in \mathbb { R } ^ { m \times n }$ has rank n, then $\parallel A ( A ^ { T } A ) ^ { - 1 } A ^ { T } \parallel _ { 2 } = 1$ .

P2.4.6 What is the nearest rank-1 matrix to

$$
A = \left[ \begin{array}{c c} 1 & M \\ 0 & 1 \end{array} \right]
$$

in the Frobenius norm?

P2.4.7 Show that if $A \in \mathbb { R } ^ { m \times n }$ , then $\parallel A \parallel _ { F } \leq \sqrt { \mathsf { r a n k } ( A ) } \parallel A \parallel _ { 2 }$ , thereby sharpening (2.3.7).

P2.4.8 Suppose $A \in \mathbb { R } ^ { n \times n }$ . Give an SVD solution to the following problem:

$$
\min _ {\det (B) = | \det (A) |} \| A - B \| _ {F}.
$$

P2.4.9 Show that if a nonzero row is added to a matrix, then both the largest and smallest singular values increase.

P2.4.10 Show that if $\theta _ { u }$ and $\theta _ { v }$ are real numbers and

$$
A = \left[ \begin{array}{l l} \cos (\theta_ {u}) & \sin (\theta_ {u}) \\ \cos (\theta_ {v}) & \sin (\theta_ {v}) \end{array} \right],
$$

then $U ^ { T } A V = \Sigma$ where

$$
U = \left[ \begin{array}{c c} \cos (\pi / 4) & - \sin (\pi / 4) \\ \sin (\pi / 4) & \cos (\pi / 4) \end{array} \right], V = \left[ \begin{array}{c c} \cos (a) & - \sin (a) \\ \sin (a) & \cos (a) \end{array} \right],
$$

and Σ = diag(√2 cos(b), √2 sin(b)) with a = (θv + θu)/2 and $b = ( \theta _ { v } - \theta _ { u } ) / 2$ .

# Notes and References for 2.4

Forsythe and Moler (SLAS) offer a good account of the SVD’s role in the analysis of the Ax = b problem. Their proof of the decomposition is more traditional than ours in that it makes use of the eigenvalue theory for symmetric matrices. Historical SVD references include:

E. Beltrami (1873). “Sulle Funzioni Bilineari,” Gionale di Mathematiche 11, 98–106.   
C. Eckart and G. Young (1939). “A Principal Axis Transformation for Non-Hermitian Matrices,” Bull. AMS 45, 118–21.   
G.W. Stewart (1993). “On the Early History of the Singular Value Decomposition,” SIAM Review 35, 551–566.

One of the most significant developments in scientific computation has been the increased use of the SVD in application areas that require the intelligent handling of matrix rank. This work started with:

C. Eckart, and G. Young (1936). “The Approximation of One Matrix by Another of Lower Rank,” Psychometrika 1, 211–218.

For generalizations of the SVD to infinite dimensional Hilbert space, see:

I.C. Gohberg and M.G. Krein (1969). Introduction to the Theory of Linear Non-Self Adjoint Operators, Amer. Math. Soc., Providence, RI.   
F. Smithies (1970). Integral Equations, Cambridge University Press, Cambridge.

Reducing the rank of a matrix as in Corollary 2.4.6 when the perturbing matrix is constrained is discussed in:

J.W. Demmel (1987). “The Smallest Perturbation of a Submatrix which Lowers the Rank and Constrained Total Least Squares Problems, SIAM J. Numer. Anal. 24, 199–206.   
G.H. Golub, A. Hoffman, and G.W. Stewart (1988). “A Generalization of the Eckart-Young-Mirsky Approximation Theorem.” Lin. Alg. Applic. 88/89, 317–328.   
G.A. Watson (1988). “The Smallest Perturbation of a Submatrix which Lowers the Rank of the Matrix,” IMA J. Numer. Anal. 8, 295–304.
