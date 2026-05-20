# 9.4 The Sign, Square Root, and Log of a Matrix

The matrix logarithm problem is the inverse of the matrix exponential problem. Not surprisingly, there is an inverse of the scaling and squaring procedure given in §9.3.1 that involves repeated matrix square roots. Thus, before we can discuss log(A) we need to understand the $\sqrt { A }$ problem. This in turn has connections to the matrix sign function and the polar decomposition.

# 9.4.1 The Matrix Sign Function

For all $z \in \mathbb { C }$ that are not on the imaginary axis, we define the sign(·) function by

$$
\operatorname{sign} (z) = \left\{ \begin{array}{l l} - 1 & \text {if} \operatorname{Re} (z) <   0, \\ + 1 & \text {if} \operatorname{Re} (z) > 0. \end{array} \right.
$$

The sign of a matrix has a particularly simple form Suppose $A \in \mathbb { C } ^ { n \times n }$ has no pure imaginary eigenvalues and that the blocks in its JCF $A = X J X ^ { - 1 }$ are ordered so that

$$
J = \left[ \begin{array}{c c} J _ {1} & 0 \\ 0 & J _ {2} \end{array} \right] \begin{array}{c} m _ {1} \\ m _ {2} \end{array}
$$

where the eigenvalues of $J _ { 1 } \in \mathbb { C } ^ { m _ { 1 } \times m _ { 1 } }$ lie in the open left half plane and the eigenvalues of $J _ { 2 } \in \mathbb { C } ^ { m _ { 2 } \times m _ { 2 } }$ lie in the open right half plane. Noting that all the derivatives of the sign function are zero, it follows from Theorem 9.1.1 that

$$
\mathrm{sign} (A) = X \left[ \begin{array}{c c} \mathrm{sign} (J _ {1}) & 0 \\ 0 & \mathrm{sign} (J _ {2}) \end{array} \right] X ^ {- 1} = X \left[ \begin{array}{c c} - I _ {m _ {1}} & 0 \\ 0 & I _ {m _ {2}} \end{array} \right] X ^ {- 1}.
$$

With the partitionings

$$
X = \left[ \begin{array}{c c} X _ {1} & X _ {2} \end{array} \right] \qquad \qquad X ^ {- H} = \left[ \begin{array}{c c} Y _ {1} & Y _ {2} \end{array} \right] \quad ,
$$

we have

$$
\operatorname{sign} (A) = X _ {2} Y _ {2} ^ {H} - X _ {1} Y _ {1} ^ {H}
$$

$$
I _ {n} = X _ {1} Y _ {1} ^ {H} + X _ {2} Y _ {2} ^ {H}
$$

and so

$$
X _ {2} Y _ {2} ^ {H} = \frac {1}{2} \left(I _ {n} + \operatorname{sign} (A)\right).
$$

Suppose apply QR-with-column pivoting to this rank- $\mathbf { \nabla } \cdot m _ { 2 }$ matrix:

$$
\frac {1}{2} \left(I _ {n} + \mathrm{sign} (A)\right) \Pi = Q R.
$$

It follows that ran $\left( Q ( : , 1 { : } m _ { 2 } ) \right) = \operatorname { r a n } ( X _ { 2 } )$ , the invariant subspace associated with A’s right half-plane eigenvalues. Thus, an approximation of sign(A) yields approximate invariant subspace information.

A number of iterative methods for computing sign(A) have been proposed. The fact that sign(z) is a zero of $g ( z ) = z ^ { 2 } - 1$ suggests a matrix analogue of the Newton iteration

$$
z _ {k + 1} = z _ {k} - \frac {g (z _ {k})}{g ^ {\prime} (z _ {k})} = \frac {1}{2} \left(z _ {k} + \frac {1}{z _ {k}}\right),
$$

i.e.,

$$
S _ {0} = A
$$

for $k = 0 , 1 , \ldots$ . (9.4.1)

$$
S _ {k + 1} = \left(S _ {k} + S _ {k} ^ {- 1}\right) / 2
$$

end

We proceed to show that this iteration is well-defined and converges to sign(A), assuming that A has no eigenvalues on the imaginary axis.

Note that if $a + b i$ is an eigenvalue of $S _ { k }$ , then

$$
\frac {1}{2} \left(a + b i + \frac {1}{a + b i}\right) = \frac {a}{2} \left(1 + \frac {1}{a ^ {2} + b ^ {2}}\right) + \frac {b}{2} \left(1 - \frac {1}{a ^ {2} + b ^ {2}}\right) i
$$

is an eigenvalue of $S _ { k + 1 }$ . Thus, if $S _ { k }$ is nonsingular, then $S _ { k + 1 }$ is nonsingular. It follows by induction that (9.4.1) is defined. Moreover, sign $( S _ { k } ) = \mathrm { s i g n } ( A )$ because an eigenvalue cannot “jump” across the imaginary axis during the iteration.

To prove that $S _ { k }$ converges to $S = \mathrm { s i g n } ( A )$ , we first observe that $S S _ { k } = S _ { k } S$ since both matrices are rational functions of A. Using this commutivity result and the identity $S ^ { 2 } = S$ , it is easy to show that

$$
S _ {k + 1} - S = \frac {1}{2} S _ {k} ^ {- 1} (S _ {k} - S) ^ {2} \tag {9.4.2}
$$

and

$$
S _ {k + 1} + S = \frac {1}{2} S _ {k} ^ {- 1} (S _ {k} + S) ^ {2}. \tag {9.4.3}
$$

If M is a matrix and sign(M ) is defined, then $M + \mathrm { s i g n } ( M )$ is nonsingular because its eigenvalues have the form $\lambda + \mathrm { s i g n } ( \lambda )$ which are clearly nonzero. Thus, the matrix

$$
S _ {k} + S = S _ {k} + \operatorname{sign} (A) = S _ {k} + \operatorname{sign} (S _ {k})
$$

is nonsingular. By manipulating equations (9.4.2) and (9.4.3) we conclude that if

$$
G _ {k} = (S _ {k} - S) (S _ {k} + S) ^ {- 1}, \tag {9.4.4}
$$

then $G _ { k + 1 } = G _ { k } ^ { 2 }$ . It follows by induction that $G _ { k } = G _ { 0 } ^ { 2 ^ { k } } . { \mathrm { ~ I f ~ } } \lambda \in \lambda ( A )$ , then

$$
\mu = \frac {\lambda - \mathrm{sign} (\lambda)}{\lambda + \mathrm{sign} (\lambda)}
$$

is an eigenvalue of $G _ { 0 } = ( A - S ) ( A + S ) ^ { - 1 }$ . Since $| \mu | < 1$ it follows from Lemma 7.3.2 that $G _ { k } \to 0$ and so

$$
S _ {k} = S (I + G _ {k}) (I - G _ {k}) ^ {- 1} \rightarrow S.
$$

Taking norms in (9.4.2) we conclude that the rate of convergence is quadratic:

$$
\| S _ {k + 1} - S \| \leq \frac {1}{2} \| S _ {k} ^ {- 1} \| \cdot \| S _ {k} - S \| ^ {2}.
$$

The overall efficiency of the method in practice is a concern since $O ( n ^ { 3 } )$ flops per iteration are required. To address this issue several enhancements of the basic iteration (9.4.1) have been proposed. One idea is to incorporate the Newton approximation

$$
S _ {k} ^ {- 1} \approx S _ {k} (2 I - S _ {k} ^ {2}).
$$

(See P9.4.1.) Using this estimate instead of the actual inverse in (9.4.1) gives update step

$$
S _ {k + 1} = \frac {1}{2} (S _ {k} + S _ {k} (2 I - S _ {k} ^ {2}) = \frac {1}{2} S _ {k} (3 I - S _ {k} ^ {2}). \tag {9.4.5}
$$

This is referred to as the Newton-Schultz iteration. Another idea is to introduce a scale factor:

$$
S _ {k + 1} = \frac {1}{2} \left((\mu_ {k} S _ {k}) + (\mu_ {k} S _ {k}) ^ {- 1}\right). \tag {9.4.6}
$$

Interesting choices for $\mu _ { k }$ include $\lvert \operatorname* { d e t } ( S _ { k } ) \rvert ^ { 1 / n } , \sqrt { \rho ( S _ { k } ^ { - 1 } ) / \rho ( S _ { k } ) }$ , and $\sqrt { \parallel S _ { k } ^ { - 1 } \parallel \parallel S _ { k } \parallel }$ where $\rho ( \cdot )$ is the spectral radius. For insights into the effective computation of the matrix sign function and related stability issues, see Kenney and Laub (1991, 1992), Higham (2007), and Higham (FOM, Chap. 5).

# 9.4.2 The Matrix Square Root

Ambiguity arises in the $f ( A )$ problem if the underlying function has branches. For example, if $f ( x ) = { \sqrt { x } }$ and

$$
A = \left[ \begin{array}{c c} 4 & 1 0 \\ 0 & 9 \end{array} \right],
$$

then

$$
A = \left[ \begin{array}{c c} 2 & 2 \\ 0 & 3 \end{array} \right] ^ {2} = \left[ \begin{array}{c c} - 2 & 1 0 \\ 0 & 3 \end{array} \right] ^ {2} = \left[ \begin{array}{c c} - 2 & - 2 \\ 0 & - 3 \end{array} \right] ^ {2} = \left[ \begin{array}{c c} 2 & - 1 0 \\ 0 & - 3 \end{array} \right] ^ {2},
$$

which shows that there are at least four legitimate choices for ${ \sqrt { A } } .$ . To clarify the situation we say F is the principal square root of A if (a) $F ^ { 2 } = A$ and (b) the eigenvalues of $F$ have positive real part. We designate this matrix by $A ^ { 1 / 2 }$ .

Analogous to the Newton iteration for scalar square roots, $x _ { k + 1 } = ( x _ { k } + a / x _ { k } ) / 2$ , we have

$$
X _ {0} = A
$$

for $k = 0 , 1 , \ldots$ . (9.4.7)

$$
X _ {k + 1} = \left(X _ {k} + X _ {k} ^ {- 1} A\right) / 2
$$

end

Notice the similarity between this iteration and the Newton sign iteration (9.4.1). Indeed, by making the substitution $X _ { k } = A ^ { 1 / 2 } S _ { k }$ in (9.4.7) we obtain the Newton sign iteration for $A ^ { 1 / 2 }$ . Global convergence and local quadratic convergence follow from what we know about (9.4.1).

Another connection between the matrix sign problem and the matrix square root problem is revealed by applying the Newton sign iteration to the matrix

$$
\tilde {A} = \left[ \begin{array}{l l} 0 & A \\ I & 0 \end{array} \right].
$$

Designate the iterates by $\tilde { S } _ { k }$ . We show by induction that $\tilde { S } _ { k }$ has the form

$$
\tilde {S} _ {k} = \left[ \begin{array}{c c} 0 & X _ {k} \\ Y _ {k} & 0 \end{array} \right].
$$

This is true for $k = 0$ by setting $X _ { 0 } = A$ and $Y _ { 0 } = I$ . To see that the result holds for $k > 0$ , observe that

$$
\tilde {S} _ {k + 1} = \frac {1}{2} \left(\tilde {S} _ {k} + \tilde {S} _ {k} ^ {- 1}\right) = \frac {1}{2} \left(\left[ \begin{array}{c c} 0 & X _ {k} \\ Y _ {k} & 0 \end{array} \right] + \left[ \begin{array}{c c} 0 & Y _ {k} ^ {- 1} \\ X _ {k} ^ {- 1} & 0 \end{array} \right]\right)
$$

and thus

$$
X _ {k + 1} = \left(X _ {k} + Y _ {k} ^ {- 1}\right) / 2, \quad Y _ {k + 1} = \left(Y _ {k} + X _ {k} ^ {- 1}\right) / 2. \tag {9.4.8}
$$

Another induction argument shows that

$$
X _ {k} = A Y _ {k}, \quad k = 0, 1, \dots , \tag {9.4.9}
$$

and so

$$
X _ {k + 1} = \left(X _ {k} + A X _ {k} ^ {- 1}\right) / 2, \quad Y _ {k + 1} = \left(Y _ {k} + A ^ {- 1} Y _ {k} ^ {- 1}\right) / 2. \tag {9.4.10}
$$

It follows that $X _ { k }  A ^ { 1 / 2 }$ and $Y _ { k }  A ^ { - 1 / 2 }$ and we have established the following identity:

$$
\mathrm{sign} \left(\left[ \begin{array}{c c} 0 & A \\ I & 0 \end{array} \right]\right) = \left[ \begin{array}{c c} 0 & A ^ {1 / 2} \\ A ^ {- 1 / 2} & 0 \end{array} \right].
$$

Equation (9.4.8) defines the Denman-Beavers iteration which turns out to have better numerical properties than (9.4.7). See Meini (2004), Higham (FOM, Chap. 6), and Higham (2008) for an analysis of these and other matrix square root algorithms.

# 9.4.3 The Polar Decomposition

If $z = a + b i \in \mathbb { C }$ is a nonzero complex number, then its polar representation is a factorization of the form $z = e ^ { i \theta } r$ where $r = \sqrt { a ^ { 2 } + b ^ { 2 } }$ and $e ^ { i \theta } = \cos ( \theta ) + i \sin ( \theta )$ i s defined by $( \cos ( \theta ) , \sin ( \theta ) ) = ( a / r , b / r )$ . The polar decomposition of a matrix is similar.

Theorem 9.4.1 (Polar Decomposition). If $A \in \mathbb { R } ^ { m \times n }$ and m $\geq n$ , then there exists a matrix $U \in \mathbb { R } ^ { m \times n }$ with orthonormal columns and a symmetric positive semidefinite $P \in \mathbb { R } ^ { n \times n }$ so that $A = U P$ .

Proof. Suppose $\ U _ { A } ^ { T } A V _ { A } \ = \ \Sigma _ { A }$ is the thin SVD of A. It is easy to show that if $\begin{array} { r c l } { U } & { = } & { U _ { A } V _ { A } ^ { \bar { T } } } \end{array}$ and $\begin{array} { l l l } { P } & { = } & { { V _ { A } } { \Sigma _ { A } } { V _ { A } ^ { T } } } \end{array}$ , then $A = U P$ and U and P have the required properties.

We refer to U as the orthogonal polar factor and P as the symmetric polar factor. Note that $P = ( A ^ { T } A ) ^ { 1 / 2 }$ and if $\mathsf { r a n k } ( A ) = n$ , then $U = A ( A ^ { T } A ) ^ { - 1 / 2 }$ . An important application of the polar decomposition is the orthogonal Procrustes problem (see §6.4.1).

Various iterative methods for computing the orthogonal polar factor have been proposed. A quadratically convergent Newton iteration for the square nonsingular case proceeds by repeatedly averaging the current iterate with the inverse of its transpose:

$$
X _ {0} = A \quad \left(\text { Assume } A \in \mathbb {R} ^ {n \times n} \text { is   nonsingular }\right)
$$

$$
\text { for } k = 0, 1, \dots \tag {9.4.11}
$$

$$
X _ {k + 1} = \left(X _ {k} + X _ {k} ^ {- T}\right) / 2
$$

end

To show that this iteration is well defined we assume that for some k the matrix $X _ { k }$ is nonsingular and that $X _ { k } = U _ { k } P _ { k }$ is its polar decomposition. It follows that

$$
X _ {k + 1} = \frac {1}{2} \left(X _ {k} + X _ {k} ^ {- T}\right) = \frac {1}{2} \left(U _ {k} P _ {k} + U _ {k} P _ {k} ^ {- 1}\right) = U _ {k} \left(\frac {P _ {k} + P _ {k} ^ {- 1}}{2}\right). \tag {9.4.12}
$$

Since the average of a positive definite matrix and its inverse is also positive definite it follows that $X _ { k + 1 }$ is nonsingular. This shows by induction that (9.4.11) is well-defined and that the $P _ { k }$ satisfy

$$
P _ {k + 1} = (P _ {k} + P _ {k} ^ {- 1}) / 2, \qquad P _ {0} = P.
$$

This is precisely the Newton sign iteration (9.4.1) with starting matrix $P _ { 0 } = P$ . Since

$$
\left\| X _ {k} - U \right\| _ {2} = \left\| U (P _ {k} - I) \right\| _ {2} = \left\| P _ {k} - I \right\| _ {2}
$$

and $P _ { k }  \mathrm { s i g n } ( P ) = I$ quadratically, we conclude that $X _ { k }$ matrices in (9.4.11) converge to U quadratically.

Extensions to the rectangular case and various ways to accelerate (9.4.11) are discussed in Higham (1986), Higham and Schreiber (1990), Gander (1990), and Kenney and Laub (1992). In this regard the matrix sign function is (once again) a handy tool for deriving algorithms. Note that if $A = U _ { A } \Sigma _ { A } V _ { A } ^ { T }$ is the SVD of $A \in \bar { \mathbb { R } } ^ { n \times n }$ and

$$
Q = \frac {1}{\sqrt {2}} \left[ \begin{array}{c c} U _ {A} & 0 \\ 0 & V _ {A} \end{array} \right] \left[ \begin{array}{c c} I _ {n} & I _ {n} \\ I _ {n} & - I _ {n} \end{array} \right]
$$

then Q is orthogonal and

$$
Q ^ {T} \left[ \begin{array}{c c} 0 & A \\ A ^ {T} & 0 \end{array} \right] Q = \left[ \begin{array}{c c} \Sigma_ {A} & 0 \\ 0 & - \Sigma_ {A} \end{array} \right].
$$

It follows that

$$
\mathrm{sign} \left(\left[ \begin{array}{c c} 0 & A \\ A ^ {T} & 0 \end{array} \right]\right) = Q \left[ \begin{array}{c c} I _ {n} & 0 \\ 0 & - I _ {n} \end{array} \right] Q ^ {T} = \left[ \begin{array}{c c} 0 & U \\ U ^ {T} & 0 \end{array} \right]
$$

where $U = U _ { A } V _ { A } ^ { T }$ is the orthogonal polar factor of A.

There is a well-developed perturbation theory for the polar decomposition. A sample result for square nonsingular matrices due to Li and Sun (2003) says that the orthogonal polar factors U and $\tilde { U }$ for nonsingular $A , \tilde { A } \in \mathbb { R } ^ { n \times n }$ satisfy the bound

$$
\| U - \tilde {U} \| _ {F} \leq \frac {4 \| A - \tilde {A} \| _ {F}}{\sigma_ {n - 1} (A) + \sigma_ {n} (A) + \sigma_ {n - 1} (\tilde {A}) + \sigma_ {n} (\tilde {A})}.
$$

# 9.4.4 The Matrix Logarithm

Given $A \in \mathbb { R } ^ { n \times n }$ , a solution to the matrix equation $e ^ { X } = A$ is a logarithm of A. Note that if $X = \log ( A )$ , then $X + 2 k \pi i$ is also a logarithm. To remove this ambiguity we define the principal logarithm as follows. If the real eigenvalues of $A \in \mathbb { R } ^ { n \times n }$ are all positive then there is a unique real matrix X that satisfies $e ^ { X } = A$ with the property that its eigenvalues satisfy $\lambda ( X ) \subset \{ z \in \mathbb { C } : - \pi < \mathsf { I m } ( z ) < \pi \}$ .

Of course, the eigenvalue-based methods of §9.2 are applicable for the log(A) problem. We discuss an approximation method that is analogous to Algorithm 9.3.1, the scaling and squaring method for the matrix exponential

As with the exponential, there are a number of different series expansions for the log function that are of computational interest. The simplest is the Maclaurin expansion:

$$
\log (A) \approx M _ {q} (A) = \sum_ {k = 1} ^ {q} (- 1) ^ {k + 1} \frac {(A - I) ^ {k}}{k}.
$$

To apply this formula we must have $\rho ( A - I ) < 1$ where $\rho ( \cdot )$ is the spectral radius.

The Gregory series expansion for log(x) yields a rational approximation:

$$
\log (A) \approx G _ {q} (A) = - 2 \sum_ {k = 0} ^ {q} \frac {1}{2 k + 1} \left((I - A) (I + A) ^ {- 1}\right) ^ {2 k + 1}.
$$

For this to converge, the real parts of A’s eigenvalues must be positive.

Diagonal Pad´e approximants are also of interest. For example, the (3,3) Pad´e approximant is given by

$$
\log (A) \approx r _ {3 3} (A) = D (A) ^ {- 1} N (A)
$$

where

$$
D (A) = 6 0 I + 9 0 (A - I) + 3 6 (A - I) ^ {2} + 3 (A - I) ^ {3},
$$

$$
N (A) = 6 0 (A - I) + 6 0 (A - I) ^ {2} + 1 1 (A - I) ^ {3}.
$$

For an approximation of this type to be effective, the matrix A must be sufficiently close to the identity matrix. Repeated square roots are one way to achieve this:

$$
k = 0
$$

$$
A _ {0} = A
$$

while $\| A - I \| >$ tol

$$
k = k + 1
$$

$$
A _ {k} = A _ {k - 1} ^ {1 / 2}
$$

end

The Denman-Beavers iteration (9.4.8) can be invoked to compute the matrix square roots. If we next compute $F \approx \log ( A _ { k } )$ by using (say) an appropriately chosen Pade approximant, then log $( A ) \ = \ 2 ^ { k } \log ( A _ { k } ) \ \approx \ 2 ^ { k } F$ . This solution framework is referred to as inverse scaling and squaring. There are many details associated with the proper implementation of this procedure and we refer the reader to Cheng, Higham. Kenney, and Laub (2001), Higham (2001), and Higham (FOM, Chap. 11).

# Problems

P9.4.1 What does the Newton iteration look like when it is applied to find a root of the function $f ( x ) = 1 / x - a ?$ Develop an inverse-free Newton iteration for solving the matrix equation $X ^ { - 1 } - A$ .

P9.4.2 Show that if $\mu _ { k } > 0$ in (9.4.6), then $\mathrm { s i g n } ( S _ { k + 1 } ) = \mathrm { s i g n } ( S _ { k } )$ .

P9.4.3 Show that s $\operatorname { i g n } ( A ) = A ( A ^ { 2 } ) ^ { - 1 / 2 }$ .

P9.4.4 Verify Equation (9.4.9).

P9.4.5 In the Denman-Beavers iteration (9.4.8), define $M _ { k } = X _ { k } Y _ { k }$ and develop a recipe for $M _ { k + 1 }$

P9.4.6 Show that if we apply the Newton square root iteration (9.4.9) to a symmetric positive definite matrix A, then $A _ { k } - A _ { k + 1 }$ is positive definite for all k.

P9.4.7 Suppose A is normal. Relate the polar factors of $e ^ { A }$ to $S = ( A - A ^ { T } ) / 2$ and $T = ( A + A ^ { T } ) / 2$

P9.4.8 Show that the polar decomposition of a nonsingular matrix is unique. Hint: If $A = U _ { 1 } P _ { 1 }$ and $A = U _ { 2 } P _ { 2 }$ are two polar decompositions, then $U _ { 2 } ^ { T } U _ { 1 } = P _ { 2 } P _ { 1 } ^ { - 1 }$ and $U _ { 1 } ^ { T } \bar { U _ { 2 } } = P _ { 1 } P _ { 2 } ^ { - 1 }$ have the same eigenvalues.

P9.4.9 Give a closed-form expression for the polar decomposition $A = U P$ of a real 2-by-2 matrix. Under what conditions is U a rotation?

P9.4.10 Give a closed-form expression for log(Q) where Q is a 2-by-2 rotation matrix.

P9.4.11 Formulate an $m < n$ version of the polar decomposition for $A \in \mathbb { R } ^ { m \times n }$ .

P9.4.12 Let A by an n-by-n symmetric positive definite matrix. (a) Show that there exists a unique symmetric positive definite X such that $\dot { \boldsymbol { A } } = \boldsymbol { X } ^ { 2 }$ . (b) Show that if $X _ { 0 } = I$ and

$$
X _ {k + 1} = (X _ {k} + A X _ {k} ^ {- 1}) / 2
$$

then $X _ { k } \to { \sqrt { A } }$ quadratically where $\sqrt { A }$ denotes the matrix X in part (a).

P9.4.13 Show that

$$
X (t) = C _ {1} \cos (t \sqrt {A}) + C _ {2} \sqrt {A ^ {- 1}} \sin (t \sqrt {A})
$$

solves the initial value problem $\ddot { X } ( t ) = - A X ( t ) , X ( 0 ) = C _ { 1 } , \dot { X } ( 0 ) = C _ { 2 }$ . Assume that A is symmetric positive definite.

# Notes and References for §9.4

Everything in this section is covered in greater depth in Higham (FOM). See also:

N.J. Higham (2005). “Functions of Matrices,” in Handbook of Linear Algebra, L. Hogben (ed.), Chapman and Hall, Boca Raton, ${ \mathrm { F L } } , \{ 1 1 - 1 - \ S 1 1 - 1 3 . $ .

Papers that discuss the ubiquitous matrix sign function and its applications include:

R. Byers (1987). “Solving the Algebraic Riccati Equation with the Matrix Sign Function,” Linear Alg. Applic. 85, 267–279.   
C.S. Kenney and A.J. Laub (1991). “Rational Iterative Methods for the Matrix Sign Function,” SIAM J. Matrix Anal. Appl. 12, 273–291.   
C.S. Kenney, A.J. Laub, and P.M. Papadopouos (1992). “Matrix Sign Algorithms for Riccati Equations,” IMA J. Math. Control Info. 9, 331–344.   
C.S. Kenney and A.J. Laub (1992). “On Scaling Newton’s Method for Polar Decomposition and the Matrix Sign Function,” SIAM J. Matrix Anal. Applic. 13, 688–706.   
R. Byers, C. He, and V. Mehrmann (1997). “The Matrix Sign Function Method and the Computation of Invariant Subspaces,” SIAM J. Matrix Anal. Applic. 18, 615–632.   
Z. Bai and J.W. Demmel (1998). “Using the Matrix Sign Function to Compute Invariant Subspaces,” SIAM J. Matrix Anal. Applic. 19, 2205–2225.   
N.J. Higham (1994). “The Matrix Sign Decomposition and Its Relation to the Polar Decomposition,” Lin. Alg. Applic. 212/213, 3–20.   
N.J. Higham, D.S. Mackey, N. Mackey, and F. Tisseur (2004). “Computing the Polar Decomposition and the Matrix Sign Decomposition in Matrix Groups,” SIAM J. Matrix Anal. Applic. 25, 1178– 1192.

Various aspects of the matrix square root problem are discussed in:

E.D. Denman and A.N. Beavers (1976). “The Matrix Sign Function and Computations in Systems,” Appl. Math. Comput., 2, 63–94.

˚A. Bj¨orck and S. Hammarling (1983). “A Schur Method for the Square Root of a Matrix,” Lin. Alg. Applic. 52/53, 127–140.

N.J. Higham (1986). “Newton’s Method for the Matrix Square Root,” Math. Comput. 46, 537–550.

N.J. Higham (1987). “Computing Real Square Roots of a Real Matrix,” Lin. Alg. Applic. 88/89, 405–430.

N.J. Higham (1997). “Stable Iterations for the Matrix Square Root,” Numer. Algorithms 15, 227–242.

Y.Y. Lu (1998). “A Pad´e Approximation Method for Square Roots of Symmetric Positive Definite Matrices,” SIAM J. Matrix Anal. Applic. 19, 833–845.

N.J. Higham, D.S. Mackey, N. Mackey, and F. Tisseur (2005). “Functions Preserving Matrix Groups and Iterations for the Matrix Square Root,” SIAM J. Matrix Anal. Applic. 26, 849–877.

C.-H. Guo and N. J. Higham (2006). “A Schur–Newton Method for the Matrix pth Root and its Inverse,” SIAM J. Matrix Anal. Applic. 28, 788–804.

B. Meini (2004). “The Matrix Square Root from a New Functional Perspective: Theoretical Results and Computational Issues,” SIAM J. Matrix Anal. Applic. 26, 362–376.

A. Frommer and B. Hashemi (2009). “Verified Computation of Square Roots of a Matrix,” SIAM J. Matrix Anal. Applic. 31, 1279–1302.   
Computational aspects of the polar decomposition and its generalizations are covered in:   
N.J. Higham (1986). “Computing the Polar Decomposition with Applications,” SIAM J. Sci. Statist. Comp. 7, 1160–1174.   
R.S. Schreiber and B.N. Parlett (1988). “Block Reflectors: Theory and Computation,” SIAM J. Numer. Anal. 25, 189–205.   
N.J. Higham and R.S. Schreiber (1990). “Fast Polar Decomposition of an Arbitrary Matrix,” SIAM J. Sci. Statist. Comput. 11, 648–655.   
N.J. Higham and P. Papadimitriou (1994). “A Parallel Algorithm for Computing the Polar Decomposition,” Parallel Comput. 20, 1161–1173.   
A.A. Dubrulle (1999). “An Optimum Iteration for the Matrix Polar Decomposition,” ETNA 8, 21–25.   
A. Zanna and H. Z. Munthe-Kaas (2002). “Generalized Polar Decompositions for the Approximation of the Matrix Exponential,” SIAM J. Matrix Anal. Applic. 23, 840–862.   
B. Laszkiewicz and K. Zietak (2006). “Approximation of Matrices and a Family of Gander Methods for Polar Decomposition,” BIT 46, 345–366.   
R. Byers and H. Xu (2008). “A New Scaling for Newton’s Iteration for the Polar Decomposition and Its Backward Stability,” SIAM J. Matrix Anal. Applic. 30, 822–843.   
N.J. Higham, C. Mehl, and F. Tisseur (2010). “The Canonical Generalized Polar Decomposition,” SIAM J. Matrix Anal. Applic. 31, 2163–2180.   
For an analysis as to whether or not the polar decomposition can be computed in a finite number of steps, see:   
A. George and Kh. Ikramov (1996). “Is The Polar Decomposition Finitely Computable?,” SIAM J. Matrix Anal. Applic. 17, 348–354.   
A. George and Kh. Ikramov (1997). “Addendum: Is The Polar Decomposition Finitely Computable?,” SIAM J. Matrix Anal. Appl. 18, 264–264.   
There is a considerable literature concerned with how the polar factors change under perturbation:   
R. Mathias (1993). “Perturbation Bounds for the Polar Decomposition,” SIAM J. Matrix Anal. Applic. 14, 588–597.   
R.-C. Li (1997). “Relative Perturbation Bounds for the Unitary Polar Factor,” BIT 37, 67–75.   
F. Chaitin-Chatelin, S. Gratton (2000). “On the Condition Numbers Associated with the Polar Factorization of a Matrix,” Numer. Lin. Alg. 7, 337–354.   
W. Li and W. Sun (2003). “New Perturbation Bounds for Unitary Polar Factors,” SIAM J. Matrix Anal. Applic. 25, 362–372.   
Finally, details concerning the matrix logarithm and its computation may be found in:   
B.W. Helton (1968). “Logarithms of Matrices,” Proc. AMS 19, 733–736.   
L. Dieci (1996). “Considerations on Computing Real Logarithms of Matrices, Hamiltonian Logarithms, and Skew-Symmetric Logarithms,” Lin. Alg. Applic. 244, 35–54.   
L. Dieci, B. Morini, and A. Papini (1996). “Computational Techniques for Real Logarithms of Matrices,” SIAM J. Matrix Anal. Applic. 17, 570–593.   
C. S. Kenney and A. J. Laub (1998). “A Schur-Fr´echet Algorithm for Computing the Logarithm and Exponential of a Matrix,” SIAM J. Matrix Anal. Applic. 19, 640–663.   
L. Dieci (1998). “Real Hamiltonian Logarithm of a Symplectic Matrix,” Lin. Alg. Applic. 281, 227–246.   
L. Dieci and A. Papini (2000). “Conditioning and Pad´e Approximation of the Logarithm of a Matrix,” SIAM J. Matrix Anal. Applic. 21, 913–930.   
N.J. Higham (2001). “Evaluating Pad´e Approximants of the Matrix Logarithm,” SIAM J. Matrix Anal. Applic. 22, 1126–1135.   
S.H. Cheng, N.J. Higham, C.S. Kenney, and A.J. Laub (2001). “Approximating the Logarithm of a Matrix to Specified Accuracy,” SIAM J. Matrix Anal. Applic. 22, 1112–1125.
