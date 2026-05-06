# 8.7 Generalized Eigenvalue Problems with Symmetry

This section is mostly about a pair of symmetrically structured versions of the generalized eigenvalue problem that we considered in §7.7. In the symmetric-definite problem we seek nontrivial solutions to the problem

$$
A x = \lambda B x \tag {8.7.1}
$$

where $A \in \mathbb { R } ^ { n \times n }$ is symmetric and $B \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite. The generalized singular value problem has the form

$$
A ^ {T} A x = \mu^ {2} B ^ {T} B x \tag {8.7.2}
$$

where $A \in \mathbb { R } ^ { m _ { 1 } \times n }$ and $B \in \mathbb { R } ^ { m _ { 2 } \times n }$ . By setting $B = I _ { n }$ we see that these problems are (respectively) generalizations of the symmetric eigenvalue problem and the singular value problem.

# 8.7.1 The Symmetric-Definite Generalized Eigenproblem

The generalized eigenvalues of the symmetric-definite pair {A, B} are denoted by $\lambda ( A , B )$ where

$$
\lambda (A, B) = \{\lambda \mid \det (A - \lambda B) = 0 \}.
$$

If $\lambda \in \lambda ( A , B )$ and x is a nonzero vector that satisfies $A x = \lambda B x$ , then x is a generalized eigenvector.

A symmetric-definite problem can be transformed to an equivalent symmetricdefinite problem with a congruence transformation:

$$
A - \lambda B \text {   is   singular   } \quad \Leftrightarrow \quad (X ^ {T} A X) - \lambda (X ^ {T} B X) \text {   is   singular. }
$$

Thus, if X is nonsingular, then $\lambda ( A , B ) = \lambda ( X ^ { T } A X , X ^ { T } B X )$ .

For a symmetric-definite pair $\{ A , B \}$ , it is possible to choose a real nonsingular X so that $X ^ { T } A X$ and $X ^ { T } B X$ are diagonal. This follows from the next result.

Theorem 8.7.1. Suppose A and B are $n { - } b y { - } n$ symmetric matrices, and define $C ( \mu )$ $b y$

$$
C (\mu) = \mu A + (1 - \mu) B \quad \mu \in \mathbb {R}. \tag {8.7.3}
$$

If there exists a $\mu \in [ 0 , 1 ]$ such that $C ( \mu )$ is nonnegative definite and

$$
\operatorname{null} (C (\mu)) = \operatorname{null} (A) \cap \operatorname{null} (B)
$$

then there exists a nonsingular X such that both $X ^ { T } A X$ and $X ^ { T } B X$ are diagonal.

Proof. Let $\mu \in [ 0 , 1 ]$ be chosen so that $C ( \mu )$ is nonnegative definite with the property that nul $( C ( \mu ) ) = { \mathsf { n u l l } } ( A ) \cap { \mathsf { n u l l } } ( B )$ . Let

$$
Q _ {1} ^ {T} C (\mu) Q _ {1} = \left[ \begin{array}{c c} D & 0 \\ 0 & 0 \end{array} \right], \qquad D = \mathrm{diag} (d _ {1}, \ldots , d _ {k}), d _ {i} > 0,
$$

be the Schur decomposition of $C ( \mu )$ and define $X _ { 1 } = Q _ { 1 } \cdot \mathrm { d i a g } ( D ^ { - 1 / 2 } , I _ { n - k } )$ . If

$$
A _ {1} = X _ {1} ^ {T} A X _ {1}, \quad B _ {1} = X _ {1} ^ {T} B X _ {1}, \quad C _ {1} = X _ {1} ^ {T} C (\mu) X _ {1},
$$

then

$$
C _ {1} = \left[ \begin{array}{c c} I _ {k} & 0 \\ 0 & 0 \end{array} \right] = \mu A _ {1} + (1 - \mu) B _ {1}.
$$

Since

$$
\operatorname{span} \left\{e _ {k + 1}, \dots , e _ {n} \right\} = \operatorname{null} \left(C _ {1}\right) = \operatorname{null} \left(A _ {1}\right) \cap \operatorname{null} \left(B _ {1}\right)
$$

it follows that $A _ { 1 }$ and $B _ { 1 }$ have the following block structure:

$$
A _ {1} = \left[ \begin{array}{c c} A _ {1 1} & 0 \\ 0 & 0 \end{array} \right] _ {n - k} ^ {k}, \qquad B _ {1} = \left[ \begin{array}{c c} B _ {1 1} & 0 \\ 0 & 0 \end{array} \right] _ {n - k} ^ {k}.
$$

Moreover $I _ { k } = \mu A _ { 1 1 } + ( 1 - \mu ) B _ { 1 1 }$ .

Suppose $\mu \neq 0$ . It then follows that if $Z ^ { T } B _ { 1 1 } Z \ = \ \mathrm { d i a g } ( b _ { 1 } , . . . , b _ { k } )$ is the Schur decomposition of $B _ { 1 1 }$ and we set

$$
X = X _ {1} \cdot \mathrm{diag} (Z, I _ {n - k})
$$

then

$$
X ^ {T} B X = \operatorname{diag} (b _ {1}, \dots , b _ {k}, 0, \dots , 0) \equiv D _ {B}
$$

and

$$
X ^ {T} A X = \frac {1}{\mu} X ^ {T} \left(C (\mu) - (1 - \mu) B\right) X = \frac {1}{\mu} \left(\left[ \begin{array}{c c} I _ {k} & 0 \\ 0 & 0 \end{array} \right] - (1 - \mu) D _ {B}\right) \equiv D _ {A}.
$$

On the other hand, if $\mu = 0$ , then let $Z ^ { T } A _ { 1 1 } Z = \mathrm { d i a g } ( a _ { 1 } , \dots , a _ { k } )$ be the Schur decomposition of $A _ { 1 1 }$ and set $X = X _ { 1 } \mathrm { d i a g } ( Z , I _ { n - k } )$ . It is easy to verify that in this case as well, both $X ^ { T } A X$ and $X ^ { T } B X$ are diagonal.

Frequently, the conditions in Theorem 8.7.1 are satisfied because either A or B is positive definite.

Corollary 8.7.2. If $A - \lambda B \in \mathbb { R } ^ { n \times n }$ is symmetric-definite, then there exists a nonsingular

$$
X = \left[ x _ {1} \mid \dots \mid x _ {n} \right]
$$

such that

$$
X ^ {T} A X = \operatorname{diag} (a _ {1}, \dots , a _ {n})
$$

and

$$
X ^ {T} B X = \operatorname{diag} (b _ {1}, \dots , b _ {n}).
$$

Moreover, $A x _ { i } = \lambda _ { i } B x _ { i }$ for i = 1:n where $\lambda _ { i } = a _ { i } / b _ { i }$ .

Proof. By setting $\mu = 0$ in Theorem 8.7.1 we see that symmetric-definite pencils can be simultaneously diagonalized. The rest of the corollary is easily verified.

Stewart (1979) has worked out a perturbation theory for symmetric pencils $A - \lambda B$ that satisfy

$$
c (A, B) = \min _ {\| x \| _ {2} = 1} (x ^ {T} A x) ^ {2} + (x ^ {T} B x) ^ {2} > 0. \tag {8.7.4}
$$

The scalar $c ( A , B )$ is called the Crawford number of the pencil $A - \lambda B$ .

Theorem 8.7.3. Suppose $A - \lambda B$ is an $n { - } b y { - } n$ symmetric-definite pencil with eigenvalues

$$
\lambda_ {1} \geq \lambda_ {2} \geq \dots \geq \lambda_ {n}.
$$

Suppose $E _ { A }$ and $E _ { B }$ are symmetric n-by-n matrices that satisfy

$$
\epsilon^ {2} = \parallel E _ {A} \parallel_ {2} ^ {2} + \parallel E _ {B} \parallel_ {2} ^ {2} <   c (A, B).
$$

Then $( A + E _ { A } ) - \lambda ( B + E _ { B } )$ is symmetric-definite with eigenvalues

$$
\mu_ {1} \geq \dots \geq \mu_ {n}
$$

that satisfy

$$
\left| \arctan \left(\lambda_ {i}\right) - \arctan \left(\mu_ {i}\right) \right| \leq \arctan \left(\epsilon / c (A, B)\right)
$$

for $i = 1 { : } n$

Proof. See Stewart (1979).

□

# 8.7.2 Simultaneous Reduction of A and B

Turning to algorithmic matters, we first present a method for solving the symmetricdefinite problem that utilizes both the Cholesky factorization and the symmetric QR algorithm.

Algorithm 8.7.1 Given $A = A ^ { T } \in \mathbb { R } ^ { n \times n }$ and $B = B ^ { T } \in \mathbb { R } ^ { n \times n }$ with B positive definite, the following algorithm computes a nonsingular X such that $X ^ { T } A X = \operatorname { d i a g } ( a _ { 1 } , . . . , a _ { n } )$ and $X ^ { T } B X \ = \ I _ { n }$ .

Compute the Cholesky factorization $B = G G ^ { T }$ using Algorithm 4.2.2.

Compute ${ \cal C } = { \cal G } ^ { - 1 } A { \cal G } ^ { - T }$ .

Use the symmetric QR algorithm to compute the Schur decomposition

$$
Q ^ {T} C Q = \operatorname{diag} \left(a _ {1}, \dots , a _ {n}\right).
$$

$\mathrm { S e t ~ } X = G ^ { - T } Q .$

This algorithm requires about $1 4 n ^ { 3 }$ flops. In a practical implementation, A can be overwritten by the matrix C. See Martin and Wilkinson (1968) for details. Note that

$$
\lambda (A, B) = \lambda (A, G G ^ {T}) = \lambda (G ^ {- 1} A G ^ {- T}, I) = \lambda (C) = \{a _ {1}, \dots , a _ {n} \}.
$$

If $\hat { a } _ { i }$ is a computed eigenvalue obtained by Algorithm 8.7.1, then it can be shown that

$$
\hat {a} _ {i} \in \lambda (G ^ {- 1} A G ^ {- T} + E _ {i})
$$

where

$$
\parallel E _ {i} \parallel_ {2} \approx \mathbf {u} \parallel A \parallel_ {2} \parallel B ^ {- 1} \parallel_ {2}.
$$

Thus, if B is ill-conditioned, then $\hat { a } _ { i }$ may be severely contaminated with roundoff error even if $a _ { i }$ is a well-conditioned generalized eigenvalue. The problem, of course, is that in this case, the matrix ${ \cal C } = { \cal G } ^ { - 1 } A { \cal G } ^ { - T }$ can have some very large entries if B, and hence G, is ill-conditioned. This difficulty can sometimes be overcome by replacing the matrix G in Algorithm 8.7.1 with $V D ^ { - 1 / 2 }$ where $V ^ { T } B V = D$ is the Schur decomposition of B. If the diagonal entries of D are ordered from smallest to largest, then the large entries in C are concentrated in the upper left-hand corner. The small eigenvalues of C can then be computed without excessive roundoff error contamination (or so the heuristic goes). For further discussion, consult Wilkinson (AEP, pp. 337–38).

The condition of the matrix X in Algorithm 8.7.1 can sometimes be improved by replacing B with a suitable convex combination of A and B. The connection between the eigenvalues of the modified pencil and those of the original are detailed in the proof of Theorem 8.7.1.

Other difficulties concerning Algorithm 8.7.1 relate to the fact that $G ^ { - 1 } A G ^ { - T }$ is generally full even when A and B are sparse. This is a serious problem, since many of the symmetric-definite problems arising in practice are large and sparse. Crawford (1973) has shown how to implement Algorithm 8.7.1 effectively when A and B are banded. Aside from this case, however, the simultaneous diagonalization approach is impractical for the large, sparse symmetric-definite problem. Alternate strategies are discussed in Chapter 10.

# 8.7.3 Other Methods

Many of the symmetric eigenvalue methods presented in earlier sections have symmetricdefinite generalizations. For example, the Rayleigh quotient iteration (8.2.6) can be extended as follows:

x0 given with $\parallel x _ { 0 } \parallel _ { 2 } = 1$

for $k = 0 , 1 , \ldots$

$$
\mu_ {k} = x _ {k} ^ {T} A x _ {k} / x _ {k} ^ {T} B x _ {k} \tag {8.7.5}
$$

Solve $( A - \mu _ { k } B ) z _ { k + 1 } = B x _ { k } { \mathrm { ~ f o r ~ } } z _ { k + 1 }$

$$
x _ {k + 1} = z _ {k + 1} / \left\| z _ {k + 1} \right\| _ {2}
$$

end

The main idea behind this iteration is that

$$
\lambda = \frac {x ^ {T} A x}{x ^ {T} B x} \tag {8.7.6}
$$

minimizes

$$
f (\lambda) = \left\| A x - \lambda B x \right\| _ {B} \tag {8.7.7}
$$

where $\| \cdot \| _ { B }$ is defined by $\| z \| _ { B } ^ { 2 } = z ^ { T } B ^ { - 1 } z$ . The mathematical properties of (8.7.5) are similar to those of (8.2.6). Its applicability depends on whether or not systems of the form $( A - \mu B ) z = x$ can be readily solved. Likewise, the same comment pertains to the generalized orthogonal iteration:

$Q _ { 0 } \in \mathbb { R } ^ { n \times p }$ given with $Q _ { 0 } ^ { T } Q _ { 0 } = I _ { p }$

for k = 1, 2, . . . (8.7.8)

Solve $B Z _ { k } = A Q _ { k - 1 } { \mathrm { ~ f o r ~ } } Z _ { k }$

$$
Z _ {k} = Q _ {k} R _ {k} \quad \left(\text { QR   factorization }, Q _ {k} \in \mathbb {R} ^ {n \times p}, R _ {k} \in \mathbb {R} ^ {p \times p}\right)
$$

end

This is mathematically equivalent to (7.3.6) with A replaced by $B ^ { - 1 } A$ . Its practicality strongly depends on how easy it is to solve linear systems of the form $B z = y$ .

# 8.7.4 The Generalized Singular Value Problem

We now turn our attention to the generalized singular value decomposition introduced in §6.1.6. This decomposition is concerned with the simultaneous diagonalization of two rectangular matrices A and B that are assumed to have the same number of columns. We restate the decomposition here with a simplification that both A and B have at least as many rows as columns. This assumption is not necessary, but it serves to unclutter our presentation of the GSVD algorithm.

Theorem 8.7.4 (Tall Rectangular Version). If $A \in \mathbb { R } ^ { m _ { 1 } \times n }$ and $B \in \mathbb { R } ^ { m _ { 2 } \times n }$ have at least as many rows as columns, then there exists an orthogonal matrix $U _ { 1 } \in \mathbb { R } ^ { m _ { 1 } \times m _ { 1 } }$ , an orthogonal matrix $U _ { 2 } \in \mathbb { R } ^ { m _ { 2 } \times m _ { 2 } }$ , and a nonsingular matrix $X \in \mathbb { R } ^ { n \times n }$ such that

$$
U _ {1} ^ {T} A X = \operatorname{diag} \left(\alpha_ {1}, \dots , \alpha_ {n}\right),
$$

$$
U _ {2} ^ {T} B X = \mathrm{diag} (\beta_ {1}, \ldots , \beta_ {n}).
$$

# Proof. See Theorem 6.1.1.

The generalized singular values of the matrix pair $\{ A , B \}$ are defined by

$$
\sigma (A, B) = \{\alpha_ {1} / \beta_ {1}, \dots , \alpha_ {n} / \beta_ {n} \}.
$$

We give names to the columns of X, $U _ { 1 }$ , and $U _ { 2 }$ . The columns of X are the right generalized singular vectors, the columns of $U _ { 1 }$ are the left-A generalized singular vectors, and the columns of $U _ { 2 }$ are the $l e f t { - } B$ generalized singular vectors. Note that

$$
A X (:, k) = \alpha_ {k} U _ {1} (:, k),
$$

$$
B X (:, k) = \beta_ {k} U _ {2} (:, k),
$$

for $k = 1 { : } n$

There is a connection between the GSVD of the matrix pair $\{ A , B \}$ and the “symmetric-definite-definite” pencil $A ^ { T } A - \lambda B ^ { T } B$ . Since

$$
X ^ {T} (A ^ {T} A - \lambda B ^ {T} B) X = D _ {A} ^ {T} D _ {A} - \lambda D _ {B} ^ {T} D _ {B} = \mathrm{diag} (\alpha_ {k} ^ {2} - \lambda \beta_ {k} ^ {2}),
$$

it follows that the right generalized singular vectors of $\{ A , B \}$ are the generalized eigenvectors for $A ^ { T } A - \lambda B ^ { T } B$ and the eigenvalues of the pencil $A ^ { T } A - \lambda B ^ { T } B$ are squares of the generalized singular values of $\{ A , B \}$ .

All these GSVD facts revert to familiar SVD facts by setting $B = I _ { n }$ . For example, if $B = I _ { n }$ , then we can set $X = U _ { 2 }$ and $U _ { 1 } ^ { T } A X = D _ { A }$ is the SVD.

We mention that the generalized singular values of {A, B} are the stationary values of

$$
\phi_ {A, B} (x) = \frac {\| A x \| _ {2}}{\| B x \| _ {2}}
$$

and the right generalized singular vectors are the associated stationary vectors. The left-A and left-B generalized singular vectors are stationary vectors associated with the quotient $\parallel y \parallel _ { 2 } / \parallel x \parallel _ { 2 }$ subject to the constraints

$$
A ^ {T} x = B ^ {T} y, \quad x \perp \operatorname{null} (A ^ {T}), \quad y \perp \operatorname{null} (A ^ {T}).
$$

See Chu, Funderlic, and Golub (1997).

A GSVD perturbation theory has been developed by Sun (1983, 1998, 2000), Paige (1984), and Li (1990).

# 8.7.5 Computing the GSVD Using the CS Decomposition

Our proof of the GSVD in Theorem 6.1.1 is constructive and makes use of the CS decomposition. In practice, computing the GSVD via the CS decomposition is a viable strategy.

Algorithm 8.7.2 (GSVD (Tall, Full-Rank Version)) Assume that $A \in \mathbb { R } ^ { m _ { 1 } \times n }$ and $B \in \mathbb { R } ^ { m _ { 2 } \times n }$ , with $m _ { 1 } \geq n , m _ { 2 } \geq n$ , and $\mathsf { n u } | | ( A ) \cap \mathsf { n u } | | ( B ) = \emptyset$ . The following algorithm computes an orthogonal matrix $U _ { 1 } \in \mathbb { R } ^ { m _ { 1 } \times m _ { 1 } }$ , an orthogonal matrix $U _ { 2 } \in \mathbb { R } ^ { m _ { 2 } \times m _ { 2 } }$ , a nonsingular matrix $X \in \mathbb { R } ^ { n \times n }$ , and diagonal matrices $D _ { A } \in \mathbb { R } ^ { m _ { 1 } \times n }$ and $D _ { B } \in \mathbb { R } ^ { m _ { 1 } \times n }$ such that $U _ { 1 } ^ { T } A X = D _ { A }$ and $U _ { 2 } ^ { T } B X = \bar { D _ { B } }$ .

Compute the the QR factorization

$$
\left[ \begin{array}{l} A \\ B \end{array} \right] = \left[ \begin{array}{l} Q _ {1} \\ Q _ {2} \end{array} \right] R.
$$

Compute the CS decomposition

$$
U _ {1} ^ {T} Q _ {1} V = D _ {A} = \operatorname{diag} \left(\alpha_ {1}, \dots , \alpha_ {n}\right),
$$

$$
U _ {2} ^ {T} Q _ {2} V = D _ {B} = \mathrm{diag} (\beta_ {1}, \ldots , \beta_ {n}).
$$

Solve RX = V for X.

The assumption that nul $\lvert ( A ) \cap \mathsf { n u l l } ( B ) = \emptyset$ is not essential. See Van Loan (1985). Regardless, the condition of the matrix X is an issue that affects accuracy. However, we point out that it is possible to compute designated right generalized singular vector subspaces without having to compute explicitly selected columns of the matrix $X =$ $V R ^ { - 1 }$ . For example, suppose that we wish to compute an orthonormal basis for the subspace $S = { \tt s p a n } \{ x _ { 1 } , \ldots x _ { k } \}$ where $x _ { i } = X ( : , i )$ . If we compute an orthogonal Z and upper triangular T so $T Z ^ { T } = V ^ { T } R$ , then

$$
Z T ^ {- 1} = R ^ {- 1} V = X
$$

and $S = \mathfrak { s p a n } \{ z _ { 1 } , \dots { } . . . z _ { k } \}$ where $z _ { i } = Z ( : , i )$ . See P5.2.2 concerning the computation of Z and T .

# 8.7.6 Computing the CS Decomposition

At first glance, the computation of the CS decomposition looks easy. After all, it is just a collection of SVDs. However, there are some complicating numerical issues that need to be addressed. To build an appreciation for this, we step through the “thin” version of the algorithm developed by Van Loan (1985) for the case

$$
Q = \left[ \begin{array}{c} Q _ {1} \\ \hline Q _ {2} \end{array} \right] = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \hline \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \end{array} \right].
$$

In exact arithmetic, the goal is to compute 5-by-5 orthogonal matrices $U _ { 1 } , U _ { 2 }$ , and V so that

$$
U _ {1} ^ {T} Q _ {1} V = C = \operatorname{diag} \left(c _ {1}, c _ {2}, c _ {3}, c _ {4}, c _ {5}\right),
$$

$$
U _ {2} ^ {T} Q _ {2} V = S = \mathrm{diag} (s _ {1}, s _ {2}, s _ {3}, s _ {4}, s _ {5}).
$$

In floating point, we strive to compute matrices $\widehat { U } _ { 2 } , \widehat { U } _ { 2 }$ and $\widehat { V }$ that are orthogonal to working precision and which transform $Q _ { 1 }$ and $Q _ { 2 }$ into nearly diagonal form:

$$
\mathsf {f l} (\widehat {U} _ {1} ^ {T} Q _ {1} \widehat {V}) = \mathrm{diag} (\hat {c} _ {k}) + E _ {1}, \quad \| E _ {1} \| \approx \mathbf {u}, \tag {8.7.9}
$$

$$
\mathsf {f l} (\widehat {U} _ {2} ^ {T} Q _ {2} \widehat {V}) = \mathrm{diag} (\widehat {s} _ {k}) + E _ {2}, \quad \| E _ {2} \| \approx \mathbf {u}. \tag {8.7.10}
$$

In what follows, it will be obvious that the computed versions of $U _ { 1 } , \ U _ { 2 }$ and V are orthogonal to working precision, as they will be “put together” from numerically sound QR factorizations and SVDs. The challenge is to affirm (8.7.9) and (8.7.10).

We start by computing the SVD

$$
U _ {2} ^ {T} Q _ {1} V = S
$$

followed by the QR factorization

$$
U _ {1} R = Q _ {1} V.
$$

Overwriting $Q _ { 2 }$ with S and $Q _ { 1 }$ with R gives

$$
Q = \left[ \begin{array}{c c c c c} r _ {1 1} & r _ {1 2} & r _ {1 3} & r _ {1 4} & r _ {1 5} \\ \epsilon_ {2 1} & r _ {2 2} & r _ {2 3} & r _ {2 4} & r _ {2 5} \\ \epsilon_ {3 1} & \epsilon_ {3 2} & r _ {3 3} & r _ {3 4} & r _ {3 5} \\ \epsilon_ {4 1} & \epsilon_ {4 2} & \epsilon_ {4 3} & r _ {4 4} & r _ {4 5} \\ \epsilon_ {5 1} & \epsilon_ {5 2} & \epsilon_ {5 3} & \epsilon_ {5 4} & r _ {5 5} \\ \hline s _ {1} & \delta_ {1 2} & \delta_ {1 3} & \delta_ {1 4} & \delta_ {2 5} \\ \delta_ {2 1} & s _ {2} & \delta_ {2 3} & \delta_ {2 4} & \delta_ {2 5} \\ \delta_ {3 1} & \delta_ {3 2} & s _ {3} & \delta_ {3 4} & \delta_ {3 5} \\ \delta_ {4 1} & \delta_ {4 2} & \delta_ {4 3} & s _ {4} & \delta_ {4 5} \\ \delta_ {5 1} & \delta_ {5 2} & \delta_ {5 3} & \delta_ {5 4} & s _ {5} \end{array} \right] \qquad \begin{array}{l} \epsilon_ {i j} = O (\mathbf {u}), \\ \delta_ {i j} = O (\mathbf {u}), \end{array}
$$

Since the columns of this matrix are orthonormal to machine precision, it follows that

$$
\left| r _ {1 1} r _ {1 j} \right| \approx \mathbf {u}, \qquad j = 2: 5.
$$

Note that if $| r _ { 1 1 } | = O ( 1 )$ , then we may conclude that √ $| r _ { 1 j } | \approx \mathbf { u }$ for $j = 2 { : } 5$ . This will be the case if (for example) $s _ { 1 } \leq 1 / \sqrt { 2 }$ for then

$$
| r _ {1 1} | \approx \sqrt {1 - s _ {1} ^ {2}} \geq \frac {1}{\sqrt {2}}.
$$

With this in mind, let us assume that the singular values $s _ { 1 } , \ldots , s _ { 5 }$ are ordered from little to big and that

$$
0 \leq s _ {1} \leq s _ {2} \leq \frac {1}{\sqrt {2}} <   s _ {3} \leq s _ {4} \leq s _ {5}. \tag {8.7.11}
$$

Working with the near-orthonormality of the columns of $Q ,$ we conclude that

$$
Q = \left[ \begin{array}{c c c c c} c _ {1} & \epsilon_ {1 2} & \epsilon_ {1 3} & \epsilon_ {1 4} & \epsilon_ {1 5} \\ \epsilon_ {2 1} & c _ {2} & \epsilon_ {2 3} & \epsilon_ {2 4} & \epsilon_ {2 5} \\ \hline \epsilon_ {3 1} & \epsilon_ {3 2} & r _ {3 3} & r _ {3 4} & r _ {3 5} \\ \epsilon_ {4 1} & \epsilon_ {4 2} & \epsilon_ {4 3} & r _ {4 4} & r _ {4 5} \\ \epsilon_ {5 1} & \epsilon_ {5 2} & \epsilon_ {5 3} & \epsilon_ {5 4} & r _ {5 5} \\ \hline s _ {1} & \delta_ {1 2} & \delta_ {1 3} & \delta_ {1 4} & \delta_ {2 5} \\ \delta_ {2 1} & s _ {2} & \delta_ {2 3} & \delta_ {2 4} & \delta_ {2 5} \\ \hline \delta_ {3 1} & \delta_ {3 2} & s _ {3} & \delta_ {3 4} & \delta_ {3 5} \\ \delta_ {4 1} & \delta_ {4 2} & \delta_ {4 3} & s _ {4} & \delta_ {4 5} \\ \delta_ {5 1} & \delta_ {5 2} & \delta_ {5 3} & \delta_ {5 4} & s _ {5} \end{array} \right] \qquad \begin{array}{l} \epsilon_ {i j} = O (\mathbf {u}), \\ \delta_ {i j} = O (\mathbf {u}). \end{array}
$$

Note that

$$
| r _ {3 4} | \approx \frac {\mathbf {u}}{| r _ {3 3} |} \approx \frac {\mathbf {u}}{\sqrt {1 - s _ {3} ^ {2}}}.
$$

Since $s _ { 3 }$ can be close to 1, we cannot guarantee that $r _ { 3 4 }$ is sufficiently small. Similar comments apply to $r _ { 3 5 }$ and $r _ { 4 5 }$ .

To rectify this we compute the SVD of Q(3:5, 3:5), taking care to apply the $U _ { - }$ matrix across rows 3 to 5 and the V matrix across columns 3:5. This gives

$$
Q = \left[ \begin{array}{c c c c c} c _ {1} & \epsilon_ {1 2} & \epsilon_ {1 3} & \epsilon_ {1 4} & \epsilon_ {1 5} \\ \epsilon_ {2 1} & c _ {2} & \epsilon_ {2 3} & \epsilon_ {2 4} & \epsilon_ {2 5} \\ \hline \epsilon_ {3 1} & \epsilon_ {3 2} & c _ {3} & \epsilon_ {3 4} & \epsilon_ {3 5} \\ \epsilon_ {4 1} & \epsilon_ {4 2} & \epsilon_ {4 3} & c _ {4} & \epsilon_ {4 5} \\ \epsilon_ {5 1} & \epsilon_ {5 2} & \epsilon_ {5 3} & \epsilon_ {5 4} & c _ {5} \\ \hline s _ {1} & \delta_ {1 2} & \delta_ {1 3} & \delta_ {1 4} & \delta_ {2 5} \\ \delta_ {2 1} & s _ {2} & \delta_ {2 3} & \delta_ {2 4} & \delta_ {2 5} \\ \hline \delta_ {3 1} & \delta_ {3 2} & t _ {3 3} & t _ {3 4} & t _ {3 5} \\ \delta_ {4 1} & \delta_ {4 2} & t _ {4 3} & t _ {4 4} & t _ {4 5} \\ \delta_ {5 1} & \delta_ {5 2} & t _ {5 3} & t _ {5 4} & t _ {5 5} \end{array} \right] \qquad \begin{array}{l} \epsilon_ {i j} = O (\mathbf {u}), \\ \delta_ {i j} = O (\mathbf {u}). \end{array}
$$

Thus, by diagonalizing the (2,2) block of $Q _ { 1 }$ we fill the (2,2) block of $Q _ { 2 }$ . However, if we compute the QR factorization of Q(8:10, 3:5) and apply the orthogonal factor across

rows 8:10, then we obtain

$$
Q = \left[ \begin{array}{c c c c c} c _ {1} & \epsilon_ {1 2} & \epsilon_ {1 3} & \epsilon_ {1 4} & \epsilon_ {1 5} \\ \epsilon_ {2 1} & c _ {2} & \epsilon_ {2 3} & \epsilon_ {2 4} & \epsilon_ {2 5} \\ \hline \epsilon_ {3 1} & \epsilon_ {3 2} & c _ {3} & \epsilon_ {3 4} & \epsilon_ {3 5} \\ \epsilon_ {4 1} & \epsilon_ {4 2} & \epsilon_ {4 3} & c _ {4} & \epsilon_ {4 5} \\ \epsilon_ {5 1} & \epsilon_ {5 2} & \epsilon_ {5 3} & \epsilon_ {5 4} & c _ {5} \\ \hline s _ {1} & \delta_ {1 2} & \delta_ {1 3} & \delta_ {1 4} & \delta_ {2 5} \\ \delta_ {2 1} & s _ {2} & \delta_ {2 3} & \delta_ {2 4} & \delta_ {2 5} \\ \hline \delta_ {3 1} & \delta_ {3 2} & t _ {3 3} & t _ {3 4} & t _ {3 5} \\ \delta_ {4 1} & \delta_ {4 2} & \delta_ {4 3} & t _ {4 4} & t _ {4 5} \\ \delta_ {5 1} & \delta_ {5 2} & \delta_ {5 3} & \delta_ {5 4} & t _ {5 5} \end{array} \right] \qquad \begin{array}{l} \epsilon_ {i j} = O (\mathbf {u}), \\ \delta_ {i j} = O (\mathbf {u}). \end{array}
$$

Using the near-orthonormality of the columns of Q and the fact that √ $c _ { 3 } , c _ { 4 }$ , and $c _ { 5 }$ are all less than $1 / \sqrt { 2 }$ , we can conclude (for example) that

$$
| t _ {3 4} | \approx O \left(\frac {\mathbf {u}}{| t _ {3 3} |}\right) \approx O \left(\frac {\mathbf {u}}{\sqrt {1 - c _ {3} ^ {2}}}\right) = O (\mathbf {u}).
$$

Using similar arguments we may conclude that both $t _ { 3 5 }$ and $t _ { 4 5 }$ are $O ( \mathbf { u } )$ . It follows that the updated $Q _ { 1 }$ and $Q _ { 2 }$ are diagonal to within the required tolerance and that (8.7.9) and (8.7.10) are achieved as a result.

# 8.7.7 The Kogbetliantz Approach

Paige (1986) developed a method for computing the GSVD based on the Kogbetliantz Jacobi SVD procedure. At each step a 2-by-2 GSVD problem is solved, a calculation that we briefly examine. Suppose $F$ and G are 2-by-2 and that G is nonsingular. If

$$
U _ {1} ^ {T} (F G ^ {- 1}) U _ {2} = \Sigma = \left[ \begin{array}{c c} \sigma_ {1} & 0 \\ 0 & \sigma_ {2} \end{array} \right]
$$

is the SVD of $F G ^ { - 1 }$ , then $\sigma ( F , G ) = \{ \sigma _ { 1 } , \sigma _ { 2 } \}$ and

$$
U _ {1} ^ {T} F = (U _ {2} ^ {T} G) \Sigma .
$$

This says that the rows of $U _ { 1 } ^ { T } F$ are parallel to the corresponding rows of $U _ { 2 } ^ { T } G$ . Thus, if Z is orthogonal so that $U _ { 2 } ^ { T } G Z = G _ { 1 }$ is upper triangular, then $U _ { 1 } ^ { \bar { T } } F Z = F _ { 1 }$ is also upper triangular. In the Paige algorithm, these 2-by-2 calculations resonate with the preservation of the triangular form that is key to the Kogbetliantz procedure. Moreover, the A and B input matrices are separately updated and the updates only involve orthogonal transformations. Although some of the calculations are very delicate, the overall procedure is tantamount to applying Kogbetliantz implicitly to the matrix $A B ^ { - 1 }$ .

# 8.7.8 Other Generalizations of the SVD

What we have been calling the “generalized singular value decomposition” is sometimes referred to as the quotient singular value decomposition or QSVD. A key feature of the decomposition is that it separately transforms the input matrices A and B in such a way that the generalized singular values and vectors are exposed, sometimes implicitly.

It turns out that there are other ways to generalize the SVD. In the product singular value decomposition problem we are given $A \in \mathbb { R } ^ { m \times n _ { 1 } }$ and $B \in \mathbb { R } ^ { m \times n _ { 2 } }$ and require the SVD of ${ \bar { A } } ^ { T } B$ . The challenge is to compute $U ^ { T } ( A ^ { T } B ) V = \Sigma$ without actually forming $A ^ { T } B$ as that operation can result in a significant loss of information. See Drma˘c (1998, 2000).

The restricted singular value decomposition involves three matrices and is best motivated from a a variational point of view. If $A \in \mathbb { R } ^ { m \times n } , B \in \mathbb { R } ^ { m \times q }$ , and $\ b { C } \in \mathbb { R } ^ { n \times p }$ , then the restricted singular values of the triplet $\{ A , B , C \}$ are the stationary values of

$$
\psi_ {A, B, C} (x, y) = \frac {y ^ {T} A x}{\| B y \| _ {2} \| C x \| _ {2}}.
$$

See Zha (1991), De Moor and Golub (1991), and Chu, De Lathauwer, and De Moor (2000). As with the product SVD, the challenge is to compute the required quantities without forming inverses and products.

All these ideas can be extended to chains of matrices, e.g., the computation of the SVD of a matrix product $A = A _ { 1 } A _ { 2 } \cdot \cdot \cdot A _ { k }$ without explicitly forming A. See De Moor and Zha (1991) and De Moor and Van Dooren (1992).

# 8.7.9 A Note on the Quadratic Eigenvalue Problem

We build on our §7.7.9 discussion of the polynomial eigenvalue problem and briefly consider some structured versions of the quadratic case,

$$
\left(\lambda^ {2} M + \lambda C + K\right) x = 0, \quad M, C, K \in \mathbb {R} ^ {n \times n}. \tag {8.7.12}
$$

We recommend the excellent survey by Tisseur and Meerbergen (2001) for more detail. Note that the eigenvalue in (8.7.12) solves the quadratic equation

$$
(x ^ {H} M x) \lambda^ {2} + (x ^ {H} C x) \lambda + (x ^ {H} K x) = 0. \tag {8.7.13}
$$

and thus

$$
\lambda = \frac {- (x ^ {H} C x) \pm \sqrt {(x ^ {H} C x) ^ {2} - 4 (x ^ {H} M x) (x ^ {H} K x)}}{2 (x ^ {H} M x)}, \tag {8.7.14}
$$

assuming that $x ^ { H } M x \neq 0$ . Linearized versions of (8.7.12) include

$$
\left[ \begin{array}{l l} 0 & N \\ K & C \end{array} \right] \left[ \begin{array}{l} x \\ u \end{array} \right] = \lambda \left[ \begin{array}{c c} N & 0 \\ 0 & - M \end{array} \right] \left[ \begin{array}{l} x \\ u \end{array} \right] \tag {8.7.15}
$$

and

$$
\left[ \begin{array}{c c} - K & 0 \\ 0 & N \end{array} \right] \left[ \begin{array}{l} x \\ u \end{array} \right] = \lambda \left[ \begin{array}{c c} C & M \\ N & 0 \end{array} \right] \left[ \begin{array}{l} x \\ u \end{array} \right] \tag {8.7.16}
$$

where $N \in \mathbb { R } ^ { n \times n }$ is nonsingular.

In many applications, the matrices M and C are symmetric and positive definite and K is symmetric and positive semidefinite. It follows from (8.7.14) that in this case the eigenvalues have nonpositive real part. If we set $N = K$ in (8.7.15), then we obtain the following generalized eigenvalue problem:

$$
\left[ \begin{array}{c c} 0 & K \\ K & C \end{array} \right] \left[ \begin{array}{c} x \\ u \end{array} \right] = \lambda \left[ \begin{array}{c c} K & 0 \\ 0 & - M \end{array} \right] \left[ \begin{array}{c} x \\ u \end{array} \right].
$$

This is not a symmetric-definite problem. However, if the overdamping condition

$$
\min _ {x ^ {T} x = 1} (x ^ {T} C x) ^ {2} - 4 (x ^ {T} M x) (x ^ {T} K x) = \gamma^ {2} > 0
$$

holds, then it can be shown that there is a scalar $\mu > 0$ so that

$$
A (\mu) = \left[ \begin{array}{c c} \mu K & K \\ K & C - \mu M \end{array} \right]
$$

is positive definite. It follows from Theorem 8.7.1 that (8.7.16) can be diagonalized by congruence. See Vesceli´c (1993).

A quadratic eigenvalue problem that arises in the analysis of gyroscopic systems has the property that $M = M ^ { T }$ (positive definite), $K = K ^ { T }$ , and $C = - C ^ { T }$ . It is easy to see from (8.7.14) that the eigenvalues are all purely imaginary. For this problem we have the structured linearization

$$
\left[ \begin{array}{c c} 0 & - K \\ M & 0 \end{array} \right] \left[ \begin{array}{c} u \\ x \end{array} \right] = \lambda \left[ \begin{array}{c c} M & C \\ 0 & M \end{array} \right] \left[ \begin{array}{c} u \\ x \end{array} \right].
$$

Notice that this is a Hamiltonian/skew-Hamiltonian generalized eigenvalue problem.

In the quadratic palindomic problem, $K = M ^ { T }$ and $C = C ^ { T }$ and the eigenvalues come in reciprocal pairs, i.e., if $Q ( \lambda )$ is singular then so is $Q ( 1 / \lambda )$ . In addition, we have the linearization

$$
\left[ \begin{array}{c c} M ^ {T} & M ^ {T} \\ C - M & M ^ {T} \end{array} \right] \left[ \begin{array}{l} y \\ z \end{array} \right] = \lambda \left[ \begin{array}{c c} - M & M ^ {T} - C \\ - M & - M \end{array} \right] \left[ \begin{array}{l} y \\ z \end{array} \right]. \tag {8.7.17}
$$

Note that if this equation holds, then

$$
(\lambda^ {2} M + \lambda C + M ^ {T}) (y + z) = 0. \tag {8.7.18}
$$

For a systematic treatment of linearizations for structured polynomial eigenvalue problems, see Mackey, Mackey, Mehl, and Mehrmann (2006).

# Problems

P8.7.1 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and $G \in \mathbb { R } ^ { n \times n }$ is lower triangular and nonsingular. Give an efficient algorithm for computing ${ \cal C } = { \cal G } ^ { - 1 } A { \cal G } ^ { - T }$ .

P8.7.2 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and $B \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite. Give an algorithm for computing the eigenvalues of AB that uses the Cholesky factorization and the symmetric

QR algorithm.

P8.7.3 Relate the principal angles and vectors between ran(A) and ran(B) to the eigenvalues and eigenvectors of the generalized eigenvalue problem

$$
\left[ \begin{array}{c c} 0 & A ^ {T} B \\ B ^ {T} A & 0 \end{array} \right] \left[ \begin{array}{c} y \\ z \end{array} \right] = \sigma \left[ \begin{array}{c c} A ^ {T} A & 0 \\ 0 & B ^ {T} B \end{array} \right] \left[ \begin{array}{c} y \\ z \end{array} \right].
$$

P8.7.4 Show that if C is real and diagonalizable, then there exist symmetric matrices A and B, B nonsingular, such that $C = A B ^ { - 1 }$ . This shows that symmetric pencils $A - \lambda B$ are essentially general.

P8.7.5 Show how to convert an $A x = \lambda B x$ problem into a generalized singular value problem if A and B are both symmetric and nonnegative definite.

P8.7.6 Given $Y \in \mathbb { R } ^ { n \times n }$ show how to compute Householder matrices $H _ { 2 } , \ldots , H _ { n }$ so that $Y H _ { n } \cdots H _ { 2 }$ $= T$ is upper triangular. Hint: $H _ { k }$ zeros out the kth row.

P8.7.7 Suppose

$$
\left[ \begin{array}{c c} 0 & A \\ A ^ {T} & 0 \end{array} \right] \left[ \begin{array}{c} y \\ z \end{array} \right] = \lambda \left[ \begin{array}{c c} B _ {1} & 0 \\ 0 & B _ {2} \end{array} \right] \left[ \begin{array}{c} y \\ z \end{array} \right]
$$

where $A \in \mathbb { R } ^ { m \times n } , \ B _ { 1 } \in \mathbb { R } ^ { m \times m }$ , and $B _ { 2 } \in \mathbb { R } ^ { n \times n }$ . Assume that $B _ { 1 }$ and $B _ { 2 }$ are positive definite with Cholesky triangles $G _ { 1 }$ and $G _ { 2 }$ respectively. Relate the generalized eigenvalues of this problem to the singular values of $G _ { 1 } ^ { - 1 } A G _ { 2 } ^ { - }$ T

P8.7.8 Suppose A and B are both symmetric positive definite. Show how to compute $\lambda ( A , B )$ and the corresponding eigenvectors using the Cholesky factorization and CS decomposition.

P8.7.9 Consider the problem

$$
\begin{array}{l} \min _ {2} \quad \| A x - b \| _ {2}, \quad A \in \mathbb {R} ^ {m \times n}, b \in \mathbb {R} ^ {m}, B, C \in \mathbb {R} ^ {n \times n}. \\ x ^ {T} B x = \beta^ {2} \\ x ^ {T} C x = \gamma^ {2} \\ \end{array}
$$

Assume that B and C are positive definite and that $Z \in \mathbb { R } ^ { n \times n }$ is a nonsingular matrix with the property that $Z ^ { T } B Z = \mathrm { d i a g } ( \lambda _ { 1 } , \ldots , \lambda _ { n } )$ and $Z ^ { T } C Z = I _ { n }$ . Assume that $\lambda _ { 1 } \geq \cdots \geq \lambda _ { n }$ . (a) Show that the the set of feasible x is empty unless $\lambda _ { n } \le \beta ^ { 2 } / \gamma ^ { 2 } \le \lambda _ { 1 }$ . (b) Using $Z ,$ , show how the two-constraint problem can be converted to a single-constraint problem of the form

$$
\min _ {y ^ {T} W y = \beta^ {2} - \lambda_ {n} \gamma^ {2}} \| \tilde {A} x - b \| _ {2}
$$

where $W = \mathrm { d i a g } ( \lambda _ { 1 } , \ldots , \lambda _ { n } ) - \lambda _ { n } I .$

P8.7.10 Show that (8.7.17) implies (8.7.18).

# Notes and References for §8.7

Just how far one can simplify a symmetric pencil A − λB via congruence is thoroughly discussed in:

P. Lancaster and L. Rodman (2005). “Canonical Forms for Hermitian Matrix Pairs under Strict Equivalence and Congruence,” SIAM Review 47, 407–443.

The sensitivity of the symmetric-definite eigenvalue problem is covered in Stewart and Sun (MPT, Chap. 6). See also:

C.R. Crawford (1976). “A Stable Generalized Eigenvalue Problem,” SIAM J. Numer. Anal. 13, 854–860.

C.-K. Li and R. Mathias (1998). “Generalized Eigenvalues of a Definite Hermitian Matrix Pair,” Lin. Alg. Applic. 271, 309–321.

S.H. Cheng and N.J. Higham (1999). “The Nearest Definite Pair for the Hermitian Generalized Eigenvalue Problem,” Lin. Alg. Applic. 302–3, 63–76.

C.-K. Li and R. Mathias (2006). “Distances from a Hermitian Pair to Diagonalizable and Nondiagonalizable Hermitian Pairs,” SIAM J. Matrix Anal. Applic. 28, 301–305.

Y. Nakatsukasa (2010). “Perturbed Behavior of a Multiple Eigenvalue in Generalized Hermitian Eigenvalue Problems,” BIT 50, 109–121.

R.-C. Li, Y. Nakatsukasa, N. Truhar, and S. Xu (2011). “Perturbation of Partitioned Hermitian Definite Generalized Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 32, 642–663.   
Although it is possible to diagonalize a symmetric-definite pencil, serious numerical issues arise if the congruence transformation is ill-conditioned. Various methods for “controlling the damage” have been proposed including:   
R.S. Martin and J.H. Wilkinson (1968). “Reduction of a Symmetric Eigenproblem Ax = λBx and Related Problems to Standard Form,” Numer. Math. 11, 99–110.   
G. Fix and R. Heiberger (1972). “An Algorithm for the Ill-Conditioned Generalized Eigenvalue Problem,” SIAM J. Numer. Anal. 9, 78–88.   
A. Bunse-Gerstner (1984). “An Algorithm for the Symmetric Generalized Eigenvalue Problem,” Lin. Alg. Applic. 58, 43–68.   
S. Chandrasekaran (2000). “An Efficient and Stable Algorithm for the Symmetric-Definite Generalized Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 21, 1202–1228.   
P.I. Davies, N.J. Higham, and F. Tisseur (2001). “Analysis of the Cholesky Method with Iterative Refinement for Solving the Symmetric Definite Generalized Eigenproblem,” SIAM J. Matrix Anal. Applic. 23, 472–493.   
F. Tisseur (2004). “Tridiagonal-Diagonal Reduction of Symmetric Indefinite Pairs,” SIAM J. Matrix Anal. Applic. 26, 215–232.   
Exploiting bandedness in A and B can be important, see:   
G. Peters and J.H. Wilkinson (1969). “Eigenvalues of Ax = λBx with Band Symmetric A and B,” Comput. J. 12, 398-404.   
C.R. Crawford (1973). “Reduction of a Band Symmetric Generalized Eigenvalue Problem,” Commun. ACM 16, 41–44.   
L. Kaufman (1993). “An Algorithm for the Banded Symmetric Generalized Matrix Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 14, 372–389.   
K. Li, T-Y. Li, and Z. Zeng (1994). “An Algorithm for the Generalized Symmetric Tridiagonal Eigenvalue Problem,” Numer. Algorithms 8, 269–291.   
The existence of a positive semidefinite linear combination of A and B was central to Theorem 8.7.1. Interestingly, the practical computation of such a combination has been addressed, see:   
C.R. Crawford (1986). “Algorithm 646 PDFIND: A Routine to Find a Positive Definite Linear Combination of Two Real Symmetric Matrices,” ACM Trans. Math. Softw. 12, 278–282.   
C.-H. Guo, N.J. Higham, and F. Tisseur (2009). “An Improved Arc Algorithm for Detecting Definite Hermitian Pairs,” SIAM J. Matrix Anal. Applic. 31, 1131–1151.   
As we mentioned, many techniques for the symmetric eigenvalue problem have natural extensions to the symmetric-definite problem. These include methods based on the Rayleigh quotient idea:   
E. Jiang(1990). “An Algorithm for Finding Generalized Eigenpairs of a Symmetric Definite Matrix Pencil,” Lin. Alg. Applic. 132, 65–91.   
R-C. Li (1994). “On Eigenvalue Variations of Rayleigh Quotient Matrix Pencils of a Definite Pencil,” Lin. Alg. Applic. 208/209, 471–483.   
There are also generalizations of the Jacobi method:   
K. Veseli˘c (1993). “A Jacobi Eigenreduction Algorithm for Definite Matrix Pairs,” Numer. Math. 64, 241–268.   
C. Mehl (2004). “Jacobi-like Algorithms for the Indefinite Generalized Hermitian Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 25, 964–985.   
Homotopy methods have also found application:   
K. Li and T-Y. Li (1993). “A Homotopy Algorithm for a Symmetric Generalized Eigenproblem,” Numer. Algorithms 4, 167–195.   
T. Zhang and K.H. Law, and G.H. Golub (1998). “On the Homotopy Method for Perturbed Symmetric Generalized Eigenvalue Problems,” SIAM J. Sci. Comput. 19, 1625–1645.   
We shall have more to say about symmetric-definite problems with general sparsity in Chapter 10. If the matrices are banded, then it is possible to implement an effective a generalization of simultaneous iteration, see:

H. Zhang and W.F. Moss (1994). “Using Parallel Banded Linear System Solvers in Generalized Eigenvalue Problems,” Parallel Comput. 20, 1089–1106.   
Turning our attention to the GSVD literature, the original references include:   
C.F. Van Loan (1976). “Generalizing the Singular Value Decomposition,” SIAM J. Numer. Anal. 13, 76–83.   
C.C. Paige and M. Saunders (1981). “Towards A Generalized Singular Value Decomposition,” SIAM J. Numer. Anal. 18, 398–405.   
The sensitivity of the GSVD is detailed in Stewart and Sun (MPT) as as well in the following papers:   
J.-G. Sun (1983). “Perturbation Analysis for the Generalized Singular Value Problem,” SIAM J. Numer. Anal. 20, 611–625.   
C.C. Paige (1984). “A Note on a Result of Sun J.-Guang: Sensitivity of the CS and GSV Decompositions,” SIAM J. Numer. Anal. 21, 186–191.   
R-C. Li (1993). “Bounds on Perturbations of Generalized Singular Values and of Associated Subspaces,” SIAM J. Matrix Anal. Applic. 14, 195–234.   
J.-G. Sun (1998). “Perturbation Analysis of Generalized Singular Subspaces,” Numer. Math. 79, 615–641.   
J.-G. Sun (2000). “Condition Number and Backward Error for the Generalized Singular Value Decomposition,” SIAM J. Matrix Anal. Applic. 22, 323–341.   
X.S. Chen and W. Li (2008). “A Note on Backward Error Analysis of the Generalized Singular Value Decomposition,” SIAM J. Matrix Anal. Applic. 30, 1358–1370.   
The variational characterization of the GSVD is analyzed in:   
M.T. Chu, R.F Funderlic, and G.H. Golub (1997). “On a Variational Formulation of the Generalized Singular Value Decomposition,” SIAM J. Matrix Anal. Applic. 18, 1082–1092.   
Connections between GSVD and the pencil A  λB are discussed in:   
B. K˚agstr¨om (1985). “The Generalized Singular Value Decomposition and the General A − λB Problem,” BIT 24, 568–583.   
Stable methods for computing the CS and generalized singular value decompositions are described in:   
G.W. Stewart (1982). “Computing the C-S Decomposition of a Partitioned Orthonormal Matrix,” Numer. Math. 40, 297–306.   
G.W. Stewart (1983). “A Method for Computing the Generalized Singular Value Decomposition,” in Matrix Pencils , B. K˚agstr¨om and A. Ruhe (eds.), Springer-Verlag, New York, 207–220.   
C.F. Van Loan (1985). “Computing the CS and Generalized Singular Value Decomposition,” Numer. Math. 46, 479–492.   
B.D. Sutton (2012). “Stable Computation of the CS Decomposition: Simultaneous Bidiagonalization,” SIAM. J. Matrix Anal. Applic. 33, 1–21.   
The idea of using the Kogbetliantz procedure for the GSVD problem is developed in:   
C.C. Paige (1986). “Computing the Generalized Singular Value Decomposition,” SIAM J. Sci. Stat. Comput. 7, 1126–1146.   
Z. Bai and H. Zha (1993). “A New Preprocessing Algorithm for the Computation of the Generalized Singular Value Decomposition,” SIAM J. Sci. Comp. 14, 1007–1012.   
Z. Bai and J.W. Demmel (1993). “Computing the Generalized Singular Value Decomposition,” SIAM J. Sci. Comput. 14, 1464–1486.   
Other methods for computing the GSVD include:   
Z. Drma˘c (1998). “A Tangent Algorithm for Computing the Generalized Singular Value Decomposition,” SIAM J. Numer. Anal. 35, 1804–1832.   
Z. Drma˘c and E.R. Jessup (2001). “On Accurate Quotient Singular Value Computation in Floating-Point Arithmetic,” SIAM J. Matrix Anal. Applic. 22, 853–873.   
S. Friedland (2005). “A New Approach to Generalized Singular Value Decomposition,” SIAM J. Matrix Anal. Applic. 27, 434–444.   
Stable methods for computing the product and restricted SVDs are discussed in the following papers:

M.T. Heath, A.J. Laub, C.C. Paige, and R.C. Ward (1986). “Computing the Singular Value Decomposition of a Product of Two Matrices,” SIAM J. Sci. Stat. Comput. 7, 1147–1159.   
K.V. Fernando and S. Hammarling (1988). “A Product-Induced Singular Value Decomposition for Two Matrices and Balanced Realization,” in Linear Algebra in Systems and Control, B.N. Datta et al (eds), SIAM Publications, Philadelphia, PA.   
B. De Moor and H. Zha (1991). “A Tree of Generalizations of the Ordinary Singular Value Decomposition,” Lin. Alg. Applic. 147, 469–500.   
H. Zha (1991). “The Restricted Singular Value Decomposition of Matrix Triplets,” SIAM J. Matrix Anal. Applic. 12, 172–194.   
B. De Moor and G.H. Golub (1991). “The Restricted Singular Value Decomposition: Properties and Applications,” SIAM J. Matrix Anal. Applic. 12, 401–425.   
B. De Moor and P. Van Dooren (1992). “Generalizing the Singular Value and QR Decompositions,” SIAM J. Matrix Anal. Applic. 13, 993–1014.   
H. Zha (1992). “A Numerical Algorithm for Computing the Restricted Singular Value Decomposition of Matrix Triplets,” Lin. Alg. Applic. 168, 1–25.   
G.E. Adams, A.W. Bojanczyk, and F.T. Luk (1994). “Computing the PSVD of Two 2×2 Triangular Matrices,” SIAM J. Matrix Anal. Applic. 15, 366–382.   
Z. Drma´c (1998). “Accurate Computation of the Product-Induced Singular Value Decomposition with Applications,” SIAM J. Numer. Anal. 35, 1969–1994.   
D. Chu, L. De Lathauwer, and B. De Moor (2000). “On the Computation of the Restricted Singular Value Decomposition via the Cosine-Sine Decomposition,” SIAM J. Matrix Anal. Applic. 22, 580–601.   
D. Chu and B.De Moor (2000). “On a variational formulation of the QSVD and the RSVD,” Lin. Alg. Applic. 311, 61–78.

For coverage of structured quadratic eigenvalue problems, see:

P. Lancaster (1991). “Quadratic Eigenvalue Problems,” Lin. Alg. Applic. 150, 499–506.   
F. Tisseur and N.J. Higham (2001). “Structured Pseudospectra for Polynomial Eigenvalue Problems, with Applications,” SIAM J. Matrix Anal. Applic. 23, 187–208.   
F. Tisseur and K. Meerbergen (2001). “The Quadratic Eigenvalue Problem,” SIAM Review 43, 235– 286.   
V. Mehrmann and D. Watkins (2002). “Polynomial Eigenvalue Problems with Hamiltonian Structure,” Electr. Trans. Numer. Anal. 13, 106–118.   
U.B. Holz, G.H. Golub, and K.H. Law (2004). “A Subspace Approximation Method for the Quadratic Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 26, 498–521.   
D.S. Mackey, N. Mackey, C. Mehl, and V. Mehrmann (2006). “Structured Polynomial Eigenvalue Problems: Good Vibrations from Good Linearizations,” SIAM. J. Matrix Anal. Applic. 28, 1029– 1051.   
B. Plestenjak (2006). “Numerical Methods for the Tridiagonal Hyperbolic Quadratic Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 28, 1157–1172.   
E.K.-W. Chu, T.-M. Hwang, W.-W. Lin, and C.-T. Wu (2008). “Vibration of Fast Trains, Palindromic Eigenvalue Problems, and Structure-Preserving Doubling Algorithms,” J. Comp. Appl. Math. 219, 237–252.
