# 2.1 Basic Ideas from Linear Algebra

This section is a quick review of linear algebra. Readers who wish a more detailed coverage should consult the references at the end of the section.

# 2.1.1 Independence, Subspace, Basis, and Dimension

A set of vectors $\{ a _ { 1 } , \ldots , a _ { n } \}$ in $\mathbb { R } ^ { m }$ is linearly independent if $\begin{array} { r } { \sum _ { j = 1 } ^ { n } \alpha _ { j } a _ { j } = 0 } \end{array}$ implies $\alpha ( 1 { : } n ) = 0$ . Otherwise, a nontrivial combination of the $a _ { i }$ is zero and $\{ a _ { 1 } , \ldots , a _ { n } \}$ is said to be linearly dependent.

A subspace of $\mathbb { R } ^ { m }$ is a subset that is also a vector space. Given a collection of vectors $\boldsymbol { a } _ { 1 } , \dots , \boldsymbol { a } _ { n } \in \mathbb { R } ^ { m }$ , the set of all linear combinations of these vectors is a subspace referred to as the span of $\{ a _ { 1 } , \ldots , a _ { n } \}$ :

$$
\operatorname{span} \left\{a _ {1}, \dots , a _ {n} \right\} = \left\{\sum_ {j = 1} ^ {n} \beta_ {j} a _ {j}: \beta_ {j} \in \mathbb {R} \right\}.
$$

If $\{ a _ { 1 } , \ldots , a _ { n } \}$ is independent and $b \in { \mathsf { s p a n } } \{ a _ { 1 } , \ldots , a _ { n } \}$ , then b is a unique linear combination of the aj. $a _ { j }$

If $S _ { 1 } , \ldots , S _ { k }$ are subspaces of $\mathbb { R } ^ { m }$ , then their sum is the subspace defined by $S = \{ \ a _ { 1 } + a _ { 2 } + \cdots + a _ { k } : a _ { i } \in S _ { i } , \ i = 1 { : } k \ \}$ . S is said to be a direct sum if each $v \in S$ has a unique representation $v = a _ { 1 } + \cdots + a _ { k }$ with $a _ { i } \in S _ { i }$ . In this case we write $S = S _ { 1 } \oplus \cdot \cdot \cdot \oplus S _ { k }$ . The intersection of the $S _ { i }$ is also a subspace, $S = S _ { 1 } \cap S _ { 2 } \cap \cdot \cdot \cdot \cap S _ { k }$ .

The subset $\{ a _ { i _ { 1 } } , \ldots , a _ { i _ { k } } \}$ is a maximal linearly independent subset of $\{ a _ { 1 } , \ldots , a _ { n } \}$ if it is linearly independent and is not properly contained in any linearly independent subset of $\{ a _ { 1 } , \ldots , a _ { n } \}$ . If $\{ a _ { i _ { 1 } } , \ldots , a _ { i _ { k } } \}$ is maximal, then span $\{ a _ { 1 } , \ldots , a _ { n } \} \ =$ span $\{ a _ { i _ { 1 } } , \ldots , a _ { i _ { k } } \}$ and $\{ a _ { i _ { 1 } } , \ldots , a _ { i _ { k } } \}$ is a basis for span $\{ a _ { 1 } , \ldots , a _ { n } \}$ . If $S \subseteq \mathbb { R } ^ { m }$ is a subspace, then it is possible to find independent basic vectors $a _ { 1 } , \dots , a _ { k } \in S$ such that $S = \mathfrak { s p a n } \{ a _ { 1 } , \dotsc , a _ { k } \}$ . All bases for a subspace $S$ have the same number of elements. This number is the dimension and is denoted by dim(S).

# 2.1.2 Range, Null Space, and Rank

There are two important subspaces associated with an m-by-n matrix A. The range of A is defined by

$$
\operatorname{ran} (A) = \{y \in \mathbb {R} ^ {m}: y = A x \text {   for   some   } x \in \mathbb {R} ^ {n} \}
$$

and the nullspace of A is defined by

$$
\operatorname{null} (A) = \{x \in \mathbb {R} ^ {n}: A x = 0 \}.
$$

If $A = \left[ { a _ { 1 } } \mid \cdots \mid { a _ { n } } \right]$ is a column partitioning, then

$$
\operatorname{ran} (A) = \operatorname{span} \left\{a _ {1}, \dots , a _ {n} \right\}.
$$

The rank of a matrix A is defined by

$$
\operatorname{rank} (A) = \dim (\operatorname{ran} (A)).
$$

If $A \in \mathbb { R } ^ { m \times n }$ , then

$$
\dim (\operatorname{null} (A)) + \operatorname{rank} (A) = n.
$$

We say that $A \in \mathbb { R } ^ { m \times n }$ is rank deficient if rank $( A ) < \operatorname* { m i n } \{ m , n \}$ . The rank of a matrix is the maximal number of linearly independent columns (or rows).

# 2.1.3 Matrix Inverse

If A and X are in $\mathbb { R } ^ { n \times n }$ and satisfy $A X = I$ , then X is the inverse of A and is denoted by $A ^ { - 1 }$ . If $A ^ { - 1 }$ exists, then A is said to be nonsingular. Otherwise, we say A is singular. The inverse of a product is the reverse product of the inverses:

$$
(A B) ^ {- 1} = B ^ {- 1} A ^ {- 1}. \tag {2.1.1}
$$

Likewise, the transpose of the inverse is the inverse of the transpose:

$$
(A ^ {- 1}) ^ {T} = (A ^ {T}) ^ {- 1} \equiv A ^ {- T}. \tag {2.1.2}
$$

# 2.1.4 The Sherman-Morrison-Woodbury Formula

The identity

$$
B ^ {- 1} = A ^ {- 1} - B ^ {- 1} (B - A) A ^ {- 1} \tag {2.1.3}
$$

shows how the inverse changes if the matrix changes. The Sherman-Morrison-Woodbury formula gives a convenient expression for the inverse of the matrix $( A + U V ^ { T } )$ where $A \in \mathbb { R } ^ { n \times n }$ and U and V are n-by-k:

$$
(A + U V ^ {T}) ^ {- 1} = A ^ {- 1} - A ^ {- 1} U (I + V ^ {T} A ^ {- 1} U) ^ {- 1} V ^ {T} A ^ {- 1}. \tag {2.1.4}
$$

A rank-k correction to a matrix results in a rank-k correction of the inverse. In (2.1.4) we assume that both A and $( I + V ^ { T } A ^ { - 1 } U )$ are nonsingular.

The $k = 1$ case is particularly useful. If $A \in \mathbb { R } ^ { n \times n }$ is nonsingular, $u , v \in \mathbb { R } ^ { n }$ , and $\alpha = 1 + v ^ { T } A ^ { - 1 } u \neq 0$ , then

$$
(A + u v ^ {T}) ^ {- 1} = A ^ {- 1} - \frac {1}{\alpha} A ^ {- 1} u v ^ {T} A ^ {- 1}. \tag {2.1.5}
$$

This is referred to as the Sherman-Morrison formula.

# 2.1.5 Orthogonality

A set of vectors $\{ x _ { 1 } , \ldots , x _ { p } \}$ in $\mathbb { R } ^ { m }$ is orthogonal if $x _ { i } ^ { T } x _ { j } = 0$ whenever $i \neq j$ and orthonormal if $x _ { i } ^ { T } x _ { j } = \delta _ { i j }$ . Intuitively, orthogonal vectors are maximally independent for they point in totally different directions.

A collection of subspaces $S _ { 1 } , \ldots , S _ { p }$ in $\mathbb { R } ^ { m }$ is mutually orthogonal if $x ^ { T } y = 0$ whenever $x \in S _ { i }$ and $y \in S _ { j }$ for $i \neq j$ . The orthogonal complement of a subspace $S \subseteq \mathbb { R } ^ { m }$ is defined by

$$
S ^ {\perp} = \{y \in \mathbb {R} ^ {m}: y ^ {T} x = 0 \text { for   all } x \in S \}.
$$

It is not hard to show that ran $( A ) ^ { \perp } = \mathsf { n u l l } ( A ^ { T } )$ . The vectors $v _ { 1 } , \ldots , v _ { k }$ form an $o r \mathrm { - }$ thonormal basis for a subspace $S \subseteq \mathbb { R } ^ { m }$ if they are orthonormal and span S.

A matrix $Q \in \mathbb { R } ^ { m \times m }$ is said to be orthogonal if $Q ^ { T } Q = I$ . If $Q = [  q _ { 1 } | \cdot \cdot \cdot | q _ { m } ]$ is orthogonal, then the $q _ { i }$ form an orthonormal basis for $\mathbb { R } ^ { m }$ . It is always possible to extend such a basis to a full orthonormal basis $\{ v _ { 1 } , \ldots , v _ { m } \}$ for $\mathbb { R } ^ { m }$ :

Theorem 2.1.1. $I f V _ { 1 } \in \mathbb { R } ^ { n \times r }$ has orthonormal columns, then there exists $V _ { 2 } \in \mathbb { R } ^ { n \times ( n - r ) }$ such that

$$
V = \left[ V _ {1} \mid V _ {2} \right]
$$

is orthogonal. Note that ran $( V _ { 1 } ) ^ { \perp } = \mathsf { r a n } ( V _ { 2 } )$

Proof. This is a standard result from introductory linear algebra. It is also a corollary of the QR factorization that we present in §5.2.

# 2.1.6 The Determinant

If $A = ( a ) \in \mathbb { R } ^ { 1 \times 1 }$ , then its determinant is given by ${ \mathsf { d e t } } ( A ) = a$ . The determinant of $A \in \mathbb { R } ^ { n \times n }$ is defined in terms of order-(n−1) determinants:

$$
\det (A) = \sum_ {j = 1} ^ {n} (- 1) ^ {j + 1} a _ {1 j} \det (A _ {1 j}).
$$

Here, $A _ { 1 j }$ is an $( n - 1 ) – \mathrm { b y } – ( n - 1 )$ matrix obtained by deleting the first row and jth column of A. Well-known properties of the determinant include det $( A B ) = \mathsf { d e t } ( A ) \mathsf { d e t } ( B )$ , de $\langle A ^ { T } \rangle = \operatorname* { d e t } ( A )$ , and det $( c A ) = c ^ { n } \mathsf { d e t } ( A )$ where A, $B \in \mathbb { R } ^ { n \times n }$ and $c \in \mathbb { R }$ . In addition, det $( A ) \neq 0$ if and only if A is nonsingular.

# 2.1.7 Eigenvalues and Eigenvectors

Until we get to the main eigenvalue part of the book (Chapters 7 and 8), we need a handful of basic properties so that we can fully appreciate the singular value decomposition (§2.4), positive definiteness (§4.2), and various fast linear equation solvers (§4.8).

The eigenvalues of $A \in \mathbb { C } ^ { n \times n }$ are the zeros of the characteristic polynomial

$$
p (x) = \det (A - x I).
$$

Thus, every n-by-n matrix has n eigenvalues. We denote the set of A’s eigenvalues by

$$
\lambda (A) = \{x: \det (A - x I) = 0 \}.
$$

If the eigenvalues of A are real, then we index them from largest to smallest as follows:

$$
\lambda_ {n} (A) \leq \dots \leq \lambda_ {2} (A) \leq \lambda_ {1} (A).
$$

In this case, we sometimes use the notation $\lambda _ { \mathrm { m a x } } ( A )$ and $\lambda _ { \mathrm { m i n } } ( A )$ to denote $\lambda _ { 1 } ( A )$ and $\lambda _ { n } ( A )$ respectively.

If $X \in \mathbb { C } ^ { n \times n }$ is nonsingular and $B = X ^ { - 1 } A X$ , then A and B are similar. If two matrices are similar, then they have exactly the same eigenvalues.

If $\lambda \in \lambda ( A )$ , then there exists a nonzero vector x so that $A x = \lambda x$ . Such a vector is said to be an eigenvector for A associated with λ. If $A \in \mathbb { C } ^ { n \times n }$ has n independent eigenvectors $x _ { 1 } , \ldots , x _ { n }$ and $A x _ { i } = \lambda _ { i } x _ { i }$ for $i = 1 { : } n$ , then A is diagonalizable. The terminology is appropriate for if

$$
X = \left[ x _ {1} \mid \dots \mid x _ {n} \right],
$$

then

$$
X ^ {- 1} A X = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right).
$$

Not all matrices are diagonalizable. However, if $A \in \mathbb { R } ^ { n \times n }$ is symmetric, then there exists an orthogonal $Q$ so that

$$
Q ^ {T} A Q = \mathrm{diag} (\lambda_ {1}, \dots , \lambda_ {n}). \tag {2.1.6}
$$

This is called the Schur decomposition. The largest and smallest eigenvalues of a symmetric matrix satisfy

$$
\lambda_ {\max} (A) = \max _ {x \neq 0} \frac {x ^ {T} A x}{x ^ {T} x} \tag {2.1.7}
$$

and

$$
\lambda_ {\min} (A) = \min _ {x \neq 0} \frac {x ^ {T} A x}{x ^ {T} x}. \tag {2.1.8}
$$

# 2.1.8 Differentiation

Suppose α is a scalar and that $A ( \alpha )$ is an m-by-n matrix with entries $a _ { i j } ( \alpha )$ . If $a _ { i j } ( \alpha )$ is a differentiable function of α for all i and $j ,$ then by $\dot { A } ( \alpha )$ we mean the matrix

$$
\dot {A} (\alpha) = \frac {d}{d \alpha} A (\alpha) = \left(\frac {d}{d \alpha} a _ {i j} (\alpha)\right) = (\dot {a} _ {i j} (\alpha)).
$$

Differentiation is a useful tool that can sometimes provide insight into the sensitivity of a matrix problem.

# Problems

P2.1.1 Show that if $A \in \mathbb { R } ^ { m \times n }$ has rank p, then there exists an $\boldsymbol { X } \in \mathbb { R } ^ { m \times p }$ and ${ \mathrm { ~ a ~ } } Y \in \mathbb { R } ^ { n \times p }$ such that $A = X Y ^ { T }$ , where rank $\langle X \rangle = \operatorname { r a n k } ( Y ) = p ,$ .

P2.1.2 Suppose $A ( \alpha ) \in \mathbb { R } ^ { m \times r }$ and $B ( \alpha ) \in \mathbb { R } ^ { r \times n }$ are matrices whose entries are differentiable functions of the scalar α. (a) Show

$$
\frac {d}{d \alpha} [ A (\alpha) B (\alpha) ] = \left[ \frac {d}{d \alpha} A (\alpha) \right] B (\alpha) + A (\alpha) \left[ \frac {d}{d \alpha} B (\alpha) \right].
$$

(b) Assuming A(α) is always nonsingular, show

$$
\frac {d}{d \alpha} \left[ A (\alpha) ^ {- 1} \right] = - A (\alpha) ^ {- 1} \left[ \frac {d}{d \alpha} A (\alpha) \right] A (\alpha) ^ {- 1}.
$$

P2.1.3 Suppose $A \in \mathbb { R } ^ { n \times n } , \ b \in \mathbb { R } ^ { n }$ and that $\begin{array} { r } { \phi ( x ) = \frac { 1 } { 2 } x ^ { T } A x - x ^ { T } b } \end{array}$ . Show that the gradient of $\phi$ is given by $\begin{array} { r } { \nabla \phi ( x ) = \frac { 1 } { 2 } ( A ^ { T } + A ) x - b . } \end{array}$ .

P2.1.4 Assume that both A and $A + u v ^ { T }$ are nonsingular where $A \in \mathbb { R } ^ { n \times n }$ and $u , v \in \mathbb { R } ^ { n }$ . Show that if x solves $( A + u v ^ { T } ) x = b ,$ then it also solves a perturbed right-hand-side problem of the form $A x = b + \alpha u$ . Give an expression for α in terms of $A , u ,$ and v.

P2.1.5 Show that a triangular orthogonal matrix is diagonal.

P2.1.6 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and nonsingular and define

$$
\tilde {A} = A + \alpha (u u ^ {T} + v v ^ {T}) + \beta (u v ^ {T} + v u ^ {T})
$$

where u, $v \in \mathbb { R } ^ { n }$ and $\alpha , \beta \in \mathbb { R }$ . Assuming that $\tilde { A }$ is nonsingular, use the Sherman-Morrison-Woodbury formula to develop a formula for $\tilde { A } ^ { - 1 }$ .

P2.1.7 Develop a symmetric version of the Sherman-Morrison-Woodbury formula that characterizes the inverse of $A + U S U ^ { T }$ where $A \in \mathbb { R } ^ { n \times n }$ and $S \in \mathbb { R } ^ { k \times k }$ are symmetric and $U \in \mathbb { R } ^ { n \times k }$ .

P2.1.8 Suppose $Q \in \mathbb { R } ^ { n \times n }$ is orthogonal and $z \in \mathbb { R } ^ { n }$ . Give an efficient algorithm for setting up an m-by-m matrix $A = \left( a _ { i j } \right)$ defined by $a _ { i j } = v ^ { T } ( Q ^ { i } ) ^ { T } ( Q ^ { j } ) v .$ -

P2.1.9 Show that if S is real and $S ^ { T } = - S$ , then I − S is nonsingular and the matrix $( I - S ) ^ { - 1 } ( I + S )$ is orthogonal. This is known as the Cayley transform of S.

P2.1.10 Refer to §1.3.10. (a) Show that if $S \in \mathbb { R } ^ { 2 n \times 2 n }$ is symplectic, then $S ^ { - 1 }$ exists and is also symplectic. (b) Show that if $M \in \mathbb { R } ^ { 2 n \times 2 n }$ is Hamiltonian and $S \in \mathbb { R } ^ { 2 n \times 2 n }$ is symplectic, then the matrix $M _ { 1 } = S ^ { - 1 } M S$ is Hamiltonian.

P2.1.11 Use (2.1.6) to prove (2.1.7) and (2.1.8).

# Notes and References for 2.1

In addition to Horn and Johnson (MA) and Horn and Johnson (TMA), the following introductory applied linear algebra texts are highly recommended:

R. Bellman (1997). Introduction to Matrix Analysis, Second Edition, SIAM Publications, Philadelphia, PA.

C. Meyer (2000). Matrix Analysis and Applied Linear Algebra, SIAM Publications, Philadelphia, PA.

D. Lay (2005). Linear Algebra and Its Applications, Third Edition, Addison-Wesley, Reading, MA.

S.J. Leon (2007). Linear Algebra with Applications, Seventh Edition, Prentice-Hall, Englewood Cliffs, NJ.

G. Strang (2009). Introduction to Linear Algebra, Fourth Edition, SIAM Publications, Philadelphia, PA.
