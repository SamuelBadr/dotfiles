# 7.2 Perturbation Theory

The act of computing eigenvalues is the act of computing zeros of the characteristic polynomial. Galois theory tells us that such a process has to be iterative if $n > 4$ and so errors arise because of finite termination. In order to develop intelligent stopping criteria we need an informative perturbation theory that tells us how to think about approximate eigenvalues and invariant subspaces.

# 7.2.1 Eigenvalue Sensitivity

An important framework for eigenvalue computation is to produce a sequence of similarity transformations $\{ X _ { k } \}$ with the property that the matrices $X _ { k } ^ { - 1 } A X _ { k }$ are progressively “more diagonal.” The question naturally arises, how well do the diagonal elements of a matrix approximate its eigenvalues?

Theorem 7.2.1 (Gershgorin Circle Theorem). $I f X ^ { - 1 } A X = D + F$ where $D =$ dia $\mathfrak { g } ( d _ { 1 } , \ldots , d _ { n } )$ and F has zero diagonal entries, then

$$
\lambda (A) \subseteq \bigcup_ {i = 1} ^ {n} D _ {i}
$$

$w h e r e ~ D _ { i } ~ = ~ \{ z \in \mathbb { C } : | z - d _ { i } | ~ \leq ~ \sum _ { j = 1 } ^ { n } | f _ { i j } | \} .$

Proof. Suppose $\lambda \in \lambda ( A )$ and assume without loss of generality that $\lambda \neq d _ { i }$ for $i = 1 { : } n$ . Since $( D - \lambda I ) + F$ is singular, it follows from Lemma 2.3.3 that

$$
1 \leq \| (D - \lambda I) ^ {- 1} F \| _ {\infty} = \sum_ {j = 1} ^ {n} \frac {| f _ {k j} |}{| d _ {k} - \lambda |}
$$

for some k, $1 \leq k \leq n$ . But this implies that $\lambda \in D _ { k }$ .

It can also be shown that if the Gershgorin disk $D _ { i }$ is isolated from the other disks, then it contains precisely one eigenvalue of A. See Wilkinson (AEP, pp. 71ff.).

For some methods it is possible to show that the computed eigenvalues are the exact eigenvalues of a matrix $A + E$ where E is small in norm. Consequently, we should understand how the eigenvalues of a matrix can be affected by small perturbations.

Theorem 7.2.2 (Bauer-Fike). If µ is an eigenvalue $\begin{array} { r } { o f A + E \in \mathbb { C } ^ { n \times n } \ a n d X ^ { - 1 } A X = } \end{array}$ $D = \mathrm { d i a g } ( \lambda _ { 1 } , \ldots , \lambda _ { n } )$ , then

$$
\min _ {\lambda \in \lambda (A)} | \lambda - \mu | \leq \kappa_ {p} (X) \| E \| _ {p}
$$

where $\| \cdot \| _ { p }$ denotes any of the p-norms.

Proof. If $\mu \in \lambda ( A )$ , then the theorem is obviously true. Otherwise if the matrix $X ^ { - 1 } ( A + E - \mu I ) X$ is singular, then so is $I + ( D \bar { \bf \Phi } - \mu I ) ^ { - 1 } ( X ^ { - 1 } E X )$ . Thus, from

Lemma 2.3.3 we obtain

$$
1 \leq \| (D - \mu I) ^ {- 1} (X ^ {- 1} E X) \| _ {p} \leq \| (D - \mu I) ^ {- 1} \| _ {p} \| X \| _ {p} \| E \| _ {p} \| X ^ {- 1} \| _ {p}.
$$

Since $( D - \mu I ) ^ { - 1 }$ is diagonal and the p-norm of a diagonal matrix is the absolute value of the largest diagonal entry, it follows that

$$
\| (D - \mu I) ^ {- 1} \| _ {p} = \max _ {\lambda \in \lambda (A)} \frac {1}{| \lambda - \mu |},
$$

completing the proof.

An analogous result can be obtained via the Schur decomposition:

Theorem 7.2.3. Let $Q ^ { H } A Q = D + N$ be a Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ as in $( 7 . 1 . 7 ) . \ I f \mu \in \lambda ( A + E )$ and p is the smallest positive integer such that $| N | ^ { p } = 0$ , then

$$
\min _ {\lambda \in \lambda (A)} | \lambda - \mu | \leq \max \{\theta , \theta^ {1 / p} \}
$$

where

$$
\theta = \| E \| _ {2} \sum_ {k = 0} ^ {p - 1} \| N \| _ {2} ^ {k}.
$$

Proof. Define

$$
\delta = \min _ {\lambda \in \lambda (A)} | \lambda - \mu | = \frac {1}{\| (\mu I - D) ^ {- 1} \| _ {2}}.
$$

The theorem is clearly true if $\delta = 0$ . If $\delta > 0$ , then $I - ( \mu I - A ) ^ { - 1 } E$ is singular and by Lemma 2.3.3 we have

$$
1 \leq \| (\mu I - A) ^ {- 1} E \| _ {2} \leq \| (\mu I - A) ^ {- 1} \| _ {2} \| E \| _ {2} \tag {7.2.1}
$$

$$
= \left\| \left(\left(\mu I - D\right) - N\right) ^ {- 1} \right\| _ {2} \left\| E \right\| _ {2}.
$$

Since $( \mu I - D ) ^ { - 1 }$ is diagonal and $| N | ^ { p } = 0$ , it follows that $( ( \mu I - D ) ^ { - 1 } N ) ^ { p } = 0$ . Thus,

$$
\left((\mu I - D) - N\right) ^ {- 1} = \sum_ {k = 0} ^ {p - 1} \left((\mu I - D) ^ {- 1} N\right) ^ {k} (\mu I - D) ^ {- 1}
$$

and so

$$
\| \left(\left(\mu I - D\right) - N\right) ^ {- 1} \| _ {2} \leq \frac {1}{\delta} \sum_ {k = 0} ^ {p - 1} \left(\frac {\| N \| _ {2}}{\delta}\right) ^ {k}.
$$

If $\delta > 1$ , then

$$
\parallel (\mu I - D) - N) ^ {- 1} \parallel_ {2} \leq \frac {1}{\delta} \sum_ {k = 0} ^ {p - 1} \parallel N \parallel_ {2} ^ {k}
$$

and so from (7.2.1), $\delta \leq \theta$ . If $\delta \leq 1$ , then

$$
\| (\mu I - D) - N) ^ {- 1} \| _ {2} \leq \frac {1}{\delta^ {p}} \sum_ {k = 0} ^ {p - 1} \| N \| _ {2} ^ {k}.
$$

By using (7.2.1) again we have $\delta ^ { p } \leq \theta$ and so $\delta \le \operatorname* { m a x } \{ \theta , \theta ^ { 1 / p } \}$ .

Theorems 7.2.2 and 7.2.3 suggest that the eigenvalues of a nonnormal matrix may be sensitive to perturbations. In particular, if $\kappa _ { 2 } ( X )$ or $\parallel N \parallel _ { 2 } ^ { p - 1 }$ is large, then small changes in A can induce large changes in the eigenvalues.

# 7.2.2 The Condition of a Simple Eigenvalue

Extreme eigenvalue sensitivity for a matrix A cannot occur if A is normal. On the other hand, nonnormality does not necessarily imply eigenvalue sensitivity. Indeed, a nonnormal matrix can have a mixture of well-conditioned and ill-conditioned eigenvalues. For this reason, it is beneficial to refine our perturbation theory so that it is applicable to individual eigenvalues and not the spectrum as a whole.

To this end, suppose that λ is a simple eigenvalue of $A \in \mathbb { C } ^ { n \times n }$ and that x and y satisfy $A x = \lambda x$ and $y ^ { H } A = \lambda y ^ { H }$ with $\parallel x \parallel _ { 2 } = \parallel y \parallel _ { 2 } = 1$ . If $Y ^ { H } A X = J$ is the Jordan decomposition with $Y ^ { H } = X ^ { - 1 }$ , then y and x are nonzero multiples of $X ( : , i )$ and $Y ( : , i )$ for some i. It follows from $1 = Y ( : , i ) ^ { H } X ( : , i )$ that $y ^ { H } x \neq 0$ , a fact that we shall use shortly.

Using classical results from function theory, it can be shown that in a neighborhood of the origin there exist differentiable $x ( \epsilon )$ and $\lambda ( \epsilon )$ such that

$$
(A + \epsilon F) x (\epsilon) = \lambda (\epsilon) x (\epsilon), \quad \| F \| _ {2} = 1,
$$

where $\lambda ( 0 ) = \lambda$ and $x ( 0 ) = x$ . By differentiating this equation with respect to 
 and setting $\epsilon = 0$ in the result, we obtain

$$
A \dot {x} (0) + F x = \dot {\lambda} (0) x + \lambda \dot {x} (0).
$$

Applying $y ^ { H }$ to both sides of this equation, dividing by $y ^ { H } x$ , and taking absolute values gives

$$
| \dot {\lambda} (0) | = \left| \frac {y ^ {H} F x}{y ^ {H} x} \right| \leq \frac {1}{| y ^ {H} x |}.
$$

The upper bound is attained if $F = y x ^ { H }$ . For this reason we refer to the reciprocal of

$$
s (\lambda) = | y ^ {H} x | \tag {7.2.2}
$$

as the condition of the eigenvalue λ.

Roughly speaking, the above analysis shows that $O ( \epsilon )$ perturbations in A can induce $\epsilon / s ( \lambda )$ changes in an eigenvalue. Thus, if $s ( \lambda )$ is small, then λ is appropriately regarded as ill-conditioned. Note that $s ( \lambda )$ is the cosine of the angle between the left and right eigenvectors associated with λ and is unique only if λ is simple.

A small $s ( \lambda )$ implies that A is near a matrix having a multiple eigenvalue. In particular, if λ is distinct and $s ( \lambda ) < 1$ , then there exists an E such that λ is a repeated eigenvalue of $A + E$ and

$$
\frac {\parallel E \parallel_ {2}}{\parallel A \parallel_ {2}} \leq \frac {s (\lambda)}{\sqrt {1 - s (\lambda) ^ {2}}}.
$$

This result is proved by Wilkinson (1972).

# 7.2.3 Sensitivity of Repeated Eigenvalues

If λ is a repeated eigenvalue, then the eigenvalue sensitivity question is more complicated. For example, if

$$
A = \left[ \begin{array}{c c} 1 & a \\ 0 & 1 \end{array} \right] \qquad \text { and } \qquad F = \left[ \begin{array}{c c} 0 & 0 \\ 1 & 0 \end{array} \right],
$$

then $\lambda ( A + \epsilon F ) = \{ 1 \pm \sqrt { \epsilon a } \}$ . Note that if $a \neq 0$ , then it follows that the eigenvalues of $A + \epsilon F$ are not differentiable at zero; their rate of change at the origin is infinite. In general, if λ is a defective eigenvalue of A, then $O ( \epsilon )$ perturbations in A can result in $O ( \epsilon ^ { 1 / p } )$ perturbations in λ if λ is associated with a p-dimensional Jordan block. See Wilkinson (AEP, pp. 77ff.) for a more detailed discussion.

# 7.2.4 Invariant Subspace Sensitivity

A collection of sensitive eigenvectors can define an insensitive invariant subspace provided the corresponding cluster of eigenvalues is isolated. To be precise, suppose

$$
Q ^ {H} A Q = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {n - r} ^ {r} \tag {7.2.3}
$$

is a Schur decomposition of A with

$$
Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ r & n - r \end{array} \right]. \tag {7.2.4}
$$

It is clear from our discussion of eigenvector perturbation that the sensitivity of the invariant subspace $\mathsf { r a n } ( Q _ { 1 } )$ depends on the distance between $\lambda ( T _ { 1 1 } )$ and $\lambda ( T _ { 2 2 } )$ . The proper measure of this distance turns out to be the smallest singular value of the linear transformation $X  T _ { 1 1 } X - X T _ { 2 2 }$ . (Recall that this transformation figures in Lemma 7.1.5.) In particular, if we define the separation between the matrices $T _ { 1 1 }$ and $T _ { 2 2 }$ by

$$
\mathsf {s e p} (T _ {1 1}, T _ {2 2}) = \min _ {X \neq 0} \frac {\| T _ {1 1} X - X T _ {2 2} \| _ {F}}{\| X \| _ {F}}, \tag {7.2.5}
$$

then we have the following general result:

Theorem 7.2.4. Suppose that (7.2.3) and (7.2.4) hold and that for any matrix $E \in \mathbb { C } ^ { n \times n }$ we partition $\bar { Q } ^ { H } E Q$ as follows:

$$
Q ^ {H} E Q = \left[ \begin{array}{c c} E _ {1 1} & E _ {1 2} \\ E _ {2 1} & E _ {2 2} \end{array} \right] _ {n - r} ^ {r}.
$$

$I f \mathsf { s e p } ( T _ { 1 1 } , T _ { 2 2 } ) > 0$ and

$$
\| E \| _ {F} \left(1 + \frac {5 \| T _ {1 2} \| _ {F}}{\operatorname{sep} (T _ {1 1} , T _ {2 2})}\right) \leq \frac {\operatorname{sep} (T _ {1 1} , T _ {2 2})}{5},
$$

then there exists a $P \in \mathbb { C } ^ { ( n - r ) \times r } \ w i t h$

$$
\| P \| _ {F} \leq 4 \frac {\| E _ {2 1} \| _ {F}}{\operatorname{sep} \left(T _ {1 1} , T _ {2 2}\right)}
$$

such that the columns of $\widetilde { Q } _ { 1 } = ( Q _ { 1 } + Q _ { 2 } P ) ( I + P ^ { H } P ) ^ { - 1 / 2 }$ are an orthonormal basis for a subspace invariant for $A + E$ .

Proof. This result is a slight recasting of Theorem 4.11 in Stewart (1973) which should be consulted for proof details. See also Stewart and Sun (MPA, p. 230). The matrix $( I + P ^ { H } P ) ^ { - 1 / 2 }$ is the inverse of the square root of the symmetric positive definite matrix ${ \cal I } + { \cal P } ^ { H } { \cal P }$ . See §4.2.4.

Corollary 7.2.5. If the assumptions in Theorem 7.2.4 hold, then

$$
\operatorname{dist} \left(\operatorname{ran} \left(Q _ {1}\right), \operatorname{ran} \left(\widetilde {Q} _ {1}\right)\right) \leq 4 \frac {\left\| E _ {2 1} \right\| _ {F}}{\operatorname{sep} \left(T _ {1 1} , T _ {2 2}\right)}.
$$

Proof. Using the SVD of P , it can be shown that

$$
\| P (I + P ^ {H} P) ^ {- 1 / 2} \| _ {2} \leq \| P \| _ {2} \leq \| P \| _ {F}. \tag {7.2.6}
$$

Since the required distance is the 2-norm of $Q _ { 2 } ^ { H } \widetilde { Q } _ { 1 } = P ( I + P ^ { H } P ) ^ { - 1 / 2 }$ , the proof is complete.

Thus, the reciprocal of ${ \mathsf { s e p } } ( T _ { 1 1 } , T _ { 2 2 } )$ can be thought of as a condition number that measures the sensitivity of $\mathsf { r a n } ( Q _ { 1 } )$ as an invariant subspace.

# 7.2.5 Eigenvector Sensitivity

If we set $r = 1$ in the preceding subsection, then the analysis addresses the issue of eigenvector sensitivity.

Corollary 7.2.6. Suppose A, $E \in \mathbb { C } ^ { n \times n }$ and that $Q = \left[ \left. q _ { 1 } \right| Q _ { 2 } \right] \in \mathbb { C } ^ { n \times n }$ is unitary with $q _ { 1 } \in \mathbb { C } ^ { n }$ . Assume

$$
Q ^ {H} A Q = \left[ \begin{array}{c c} \lambda & v ^ {H} \\ 0 & T _ {2 2} \end{array} \right] _ {n - 1} ^ {1}, \qquad Q ^ {H} E Q = \left[ \begin{array}{c c} \epsilon & \gamma^ {H} \\ \delta & E _ {2 2} \end{array} \right] _ {n - 1} ^ {1}.
$$

(Thus, $q _ { 1 }$ is an eigenvector.) If $\sigma = \sigma _ { \operatorname* { m i n } } ( T _ { 2 2 } - \lambda I ) > 0$ and

$$
\left\| E \right\| _ {F} \left(1 + \frac {5 \| v \| _ {2}}{\sigma}\right) \leq \frac {\sigma}{5},
$$

then there exists $p \in \mathbb { C } ^ { n - 1 }$ with

$$
\| p \| _ {2} \leq 4 \frac {\| \delta \| _ {2}}{\sigma}
$$

such that $\tilde { q } _ { 1 } = ( q _ { 1 } + Q _ { 2 } p ) / \sqrt { 1 + p ^ { H } p }$ is a unit 2-norm eigenvector for $A + E$ . Moreover,

$$
\operatorname{dist} \left(\operatorname{span} \left\{q _ {1} \right\}, \operatorname{span} \left\{\tilde {q} _ {1} \right\}\right) \leq 4 \frac {\| \delta \| _ {2}}{\sigma}.
$$

Proof. The result follows from Theorem 7.2.4, Corollary 7.2.5, and the observation that if $T _ { 1 1 } = \lambda$ , then sep $( T _ { 1 1 } , T _ { 2 2 } ) = \sigma _ { \operatorname* { m i n } } ( T _ { 2 2 } - \lambda I )$ .

Note that $\sigma _ { \operatorname* { m i n } } ( T _ { 2 2 } - \lambda I )$ roughly measures the separation of λ from the eigenvalues of $T _ { 2 2 }$ . We have to say “roughly” because

$$
\operatorname{sep} (\lambda , T _ {2 2}) = \sigma_ {\min} (T _ {2 2} - \lambda I) \leq \min _ {\mu \in \lambda (T _ {2 2})} | \mu - \lambda |
$$

and the upper bound can be a gross overestimate.

That the separation of the eigenvalues should have a bearing upon eigenvector sensitivity should come as no surprise. Indeed, if λ is a nondefective, repeated eigenvalue, then there are an infinite number of possible eigenvector bases for the associated invariant subspace. The preceding analysis merely indicates that this indeterminancy begins to be felt as the eigenvalues coalesce. In other words, the eigenvectors associated with nearby eigenvalues are “wobbly.”

# Problems

P7.2.1 Suppose $Q ^ { H } A Q = \mathrm { d i a g } ( \lambda _ { 1 } ) +$ N is a Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ and define $\nu ( A ) =$ $\parallel \boldsymbol { A } ^ { H } \boldsymbol { A } - \boldsymbol { \bar { A } } \boldsymbol { \bar { A } } ^ { H } \parallel _ { F } .$ . The upper and lower bounds in

$$
\frac {\nu (A) ^ {2}}{6 \| A \| _ {F} ^ {2}} \leq \| N \| _ {F} ^ {2} \leq \sqrt {\frac {n ^ {3} - n}{1 2}} \nu (A)
$$

are established by Henrici (1962) and Eberlein (1965), respectively. Verify these results for the case $n = 2$ .

P7.2.2 Suppose $A \in \mathbb { C } ^ { n \times n }$ and $X ^ { - 1 } A X = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ with distinct $\lambda _ { i }$ . Show that if the columns of X have unit 2-norm, then $\kappa _ { F } ( X ) ^ { 2 } = n ( 1 / s ( \bar { \lambda _ { 1 } } ) ^ { 2 } + \cdot \cdot \cdot + 1 / s ( \lambda _ { n } ) ^ { 2 } )$ .

P7.2.3 Suppose $Q ^ { H } A Q = \mathrm { d i a g } ( \lambda _ { i } ) + N$ is a Schur decomposition of A and that $X ^ { - 1 } A X = \mathrm { d i a g } \left( \lambda _ { i } \right)$ . Show κ2 $( \bar { X } ) ^ { 2 } \geq 1 + ( \| \ N \| _ { F } / \| \bar { \ A } \| _ { F } ) ^ { 2 }$ . See Loizou (1969).

P7.2.4 If $X ^ { - 1 } A X = \mathrm { { d i a g } } \left( \lambda _ { i } \right)$ and $| \lambda _ { 1 } | \geq \cdots \geq | \lambda _ { n } |$ , then

$$
\frac {\sigma_ {i} (A)}{\kappa_ {2} (X)} \leq | \lambda_ {i} | \leq \kappa_ {2} (X) \sigma_ {i} (A).
$$

Prove this result for the $n = 2$ case. See Ruhe (1975).

P7.2.5 Show that if $A = { \left[ \begin{array} { l l } { a } & { c } \\ { 0 } & { b } \end{array} \right] }$ and $a \neq b ,$ then $s ( a ) = s ( b ) = ( 1 + | c / ( a - b ) | ^ { 2 } ) ^ { - 1 / 2 }$ .

# P7.2.6 Suppose

$$
A = \left[ \begin{array}{c c} \lambda & v ^ {T} \\ 0 & T _ {2 2} \end{array} \right]
$$

and that $\lambda \not \in \lambda ( T _ { 2 2 } )$ . Show that if $\sigma = { \tt s e p } ( \lambda , T _ { 2 2 } )$ , then

$$
s (\lambda) = \frac {1}{\sqrt {1 + \| (T _ {2 2} - \lambda I) ^ {- 1} v \| _ {2} ^ {2}}} \leq \frac {\sigma}{\sqrt {\sigma^ {2} + \| v \| _ {2} ^ {2}}}.
$$

where $s ( \lambda )$ is defined in (7.2.2).

P7.2.7 Show that the condition of a simple eigenvalue is preserved under unitary similarity transformations.

P7.2.8 With the same hypothesis as in the Bauer-Fike theorem (Theorem 7.2.2), show that

$$
\min _ {\lambda \in \lambda (A)} | \lambda - \mu | \leq \| | X ^ {- 1} | | E | | X | \| _ {p}.
$$

P7.2.9 Verify (7.2.6).

P7.2.10 Show that if B ∈ Cm×m $\boldsymbol { B } \in \mathbb { C } ^ { m \times m }$ and $C \in \mathbb { C } ^ { n \times n }$ , then sep(B, C) is less than or equal to $| \lambda - \mu |$ for all $\lambda \in \lambda ( B )$ and $\mu \in \lambda ( C )$ .

# Notes and References for §7.2

Many of the results presented in this section may be found in Wilkinson (AEP), Stewart and Sun (MPA) as well as:

F.L. Bauer and C.T. Fike (1960). “Norms and Exclusion Theorems,” Numer. Math. 2, 123–44.

A.S. Householder (1964). The Theory of Matrices in Numerical Analysis. Blaisdell, New York.

R. Bhatia (2007). Perturbation Bounds for Matrix Eigenvalues, SIAM Publications, Philadelphia, PA.

Early papers concerned with the effect of perturbations on the eigenvalues of a general matrix include:

A. Ruhe (1970). “Perturbation Bounds for Means of Eigenvalues and Invariant Subspaces,” BIT 10, 343–54.

A. Ruhe (1970). “Properties of a Matrix with a Very Ill-Conditioned Eigenproblem,” Numer. Math. 15, 57–60.

J.H. Wilkinson (1972). “Note on Matrices with a Very Ill-Conditioned Eigenproblem,” Numer. Math. 19, 176–78.

W. Kahan, B.N. Parlett, and E. Jiang (1982). “Residual Bounds on Approximate Eigensystems of Nonnormal Matrices,” SIAM J. Numer. Anal. 19, 470–484.

J.H. Wilkinson (1984). “On Neighboring Matrices with Quadratic Elementary Divisors,” Numer. Math. 44, 1-21.

Wilkinson’s work on nearest defective matrices is typical of a growing body of literature that is concerned with “nearness” problems, see:

A. Ruhe (1987). “Closest Normal Matrix Found!,” BIT 27, 585-598.

J.W. Demmel (1987). “On the Distance to the Nearest Ill-Posed Problem,” Numer. Math. 51, 251–289.

J.W. Demmel (1988). “The Probability that a Numerical Analysis Problem is Difficult,” Math. Comput. 50, 449–480.

N.J. Higham (1989). “Matrix Nearness Problems and Applications,” in Applications of Matrix Theory, M.J.C. Gover and S. Barnett (eds.), Oxford University Press, Oxford, 1–27.

A.N. Malyshev (1999). “A Formula for the 2-norm Distance from a Matrix to the Set of Matrices with Multiple Eigenvalues,” Numer. Math. 83, 443–454.

J.-M. Gracia (2005). “Nearest Matrix with Two Prescribed Eigenvalues,” Lin. Alg. Applic. 401, 277–294.

An important subset of this literature is concerned with nearness to the set of unstable matrices. A matrix is unstable if it has an eigenvalue with nonnegative real part. Controllability is a related notion, see:

C. Van Loan (1985). “How Near is a Stable Matrix to an Unstable Matrix?,” Contemp. Math. 47, 465–477.   
J.W. Demmel (1987). “A Counterexample for two Conjectures About Stability,” IEEE Trans. Autom. Contr. AC-32, 340–342.   
R. Byers (1988). “A Bisection Method for Measuring the distance of a Stable Matrix to the Unstable Matrices,” J. Sci. Stat. Comput. 9, 875–881.   
J.V. Burke and M.L. Overton (1992). “Stable Perturbations of Nonsymmetric Matrices,” Lin. Alg. Applic. 171, 249–273.   
C. He and G.A. Watson (1998). “An Algorithm for Computing the Distance to Instability,” SIAM J. Matrix Anal. Applic. 20, 101–116.   
M. Gu, E. Mengi, M.L. Overton, J. Xia, and J. Zhu (2006). “Fast Methods for Estimating the Distance to Uncontrollability,” SIAM J. Matrix Anal. Applic. 28, 477–502.   
Aspects of eigenvalue condition are discussed in:   
C. Van Loan (1987). “On Estimating the Condition of Eigenvalues and Eigenvectors,” Lin. Alg. Applic. 88/89, 715–732.   
C.D. Meyer and G.W. Stewart (1988). “Derivatives and Perturbations of Eigenvectors,” SIAM J. Numer. Anal. 25, 679–691.   
G.W. Stewart and G. Zhang (1991). “Eigenvalues of Graded Matrices and the Condition Numbers of Multiple Eigenvalues,” Numer. Math. 58, 703–712.   
J.-G. Sun (1992). “On Condition Numbers of a Nondefective Multiple Eigenvalue,” Numer. Math. 61, 265–276.   
S.M. Rump (2001). “Computational Error Bounds for Multiple or Nearly Multiple Eigenvalues,” Lin. Alg. Applic. 324, 209–226.   
The relationship between the eigenvalue condition number, the departure from normality, and the condition of the eigenvector matrix is discussed in:   
P. Henrici (1962). “Bounds for Iterates, Inverses, Spectral Variation and Fields of Values of Nonnormal Matrices,” Numer. Math. 4, 24–40.   
P. Eberlein (1965). “On Measures of Non-Normality for Matrices,” AMS Monthly 72, 995–996.   
R.A. Smith (1967). “The Condition Numbers of the Matrix Eigenvalue Problem,” Numer. Math. 10 232–240.   
G. Loizou (1969). “Nonnormality and Jordan Condition Numbers of Matrices,” J. ACM 16, 580–640.   
A. van der Sluis (1975). “Perturbations of Eigenvalues of Non-normal Matrices,” Commun. ACM 18, 30–36.   
S.L. Lee (1995). “A Practical Upper Bound for Departure from Normality,” SIAM J. Matrix Anal. Applic. 16, 462–468.   
Gershgorin’s theorem can be used to derive a comprehensive perturbation theory. The theorem itself can be generalized and extended in various ways, see:   
R.S. Varga (1970). “Minimal Gershgorin Sets for Partitioned Matrices,” SIAM J. Numer. Anal. 7, 493–507.   
R.J. Johnston (1971). “Gershgorin Theorems for Partitioned Matrices,” Lin. Alg. Applic. 4, 205–20.   
R.S. Varga and A. Krautstengl (1999). “On Gergorin-type Problems and Ovals of Cassini,” ETNA 8, 15–20.   
R.S. Varga (2001). “Gergorin-type Eigenvalue Inclusion Theorems and Their Sharpness,” ETNA 12, 113–133.   
C. Beattie and I.C.F. Ipsen (2003). “Inclusion Regions for Matrix Eigenvalues,” Lin. Alg. Applic. 358, 281–291.   
In our discussion, the perturbations to the A-matrix are general. More can be said when the perturbations are structured, see:   
G.W. Stewart (2001). “On the Eigensystems of Graded Matrices,” Numer. Math. 90, 349–370.   
J. Moro and F.M. Dopico (2003). “Low Rank Perturbation of Jordan Structure,” SIAM J. Matrix Anal. Applic. 25, 495–506.   
R. Byers and D. Kressner (2004). “On the Condition of a Complex Eigenvalue under Real Perturbations,” BIT 44, 209–214.   
R. Byers and D. Kressner (2006). “Structured Condition Numbers for Invariant Subspaces,” SIAM J. Matrix Anal. Applic. 28, 326–347.

An absolute perturbation bound comments on the difference between an eigenvalue λ and its perturbation λ˜. A relative perturbation bound examines the quotient $| \lambda - { \tilde { \lambda } } | / | { \bar { \lambda } } | ,$ , something that can be very important when there is a concern about a small eigenvalue. For results in this direction consult:

R.-C. Li (1997). “Relative Perturbation Theory. III. More Bounds on Eigenvalue Variation,” Lin. Alg. Applic. 266, 337–345.

S.C. Eisenstat and I.C.F. Ipsen (1998). “Three Absolute Perturbation Bounds for Matrix Eigenvalues Imply Relative Bounds,” SIAM J. Matrix Anal. Applic. 20, 149–158.

S.C. Eisenstat and I.C.F. Ipsen (1998). “Relative Perturbation Results for Eigenvalues and Eigenvectors of Diagonalisable Matrices,” BIT 38, 502–509.

I.C.F. Ipsen (1998). “Relative Perturbation Results for Matrix Eigenvalues and Singular Values,” Acta Numerica, 7, 151–201.

I.C.F. Ipsen (2000). “Absolute and Relative Perturbation Bounds for Invariant Subspaces of Matrices,” Lin. Alg. Applic. 309, 45–56.

I.C.F. Ipsen (2003). “A Note on Unifying Absolute and Relative Perturbation Bounds,” Lin. Alg. Applic. 358, 239–253.

Y. Wei, X. Li, F. Bu, and F. Zhang (2006). “Relative Perturbation Bounds for the Eigenvalues of Diagonalizable and Singular Matrices–Application to Perturbation Theory for Simple Invariant Subspaces,” Lin. Alg. Applic. 419, 765-771.

The eigenvectors and invariant subspaces of a matrix also “move” when there are perturbations. Tracking these changes is typically more challenging than tracking changes in the eigenvalues, see:

T. Kato (1966). Perturbation Theory for Linear Operators, Springer-Verlag, New York.

C. Davis and W.M. Kahan (1970). “The Rotation of Eigenvectors by a Perturbation, III,” SIAM J. Numer. Anal. 7, 1–46.

G.W. Stewart (1971). “Error Bounds for Approximate Invariant Subspaces of Closed Linear Operators,” SIAM. J. Numer. Anal. 8, 796–808.

G.W. Stewart (1973). “Error and Perturbation Bounds for Subspaces Associated with Certain Eigenvalue Problems,” SIAM Review 15, 727–764.

J. Xie (1997). “A Note on the Davis-Kahan sin(2θ) Theorem,” Lin. Alg. Applic. 258, 129–135.

S.M. Rump and J.-P.M. Zemke (2003). “On Eigenvector Bounds,” BIT 43, 823–837.

Detailed analyses of the function sep(.,.) and the map $X  A X + X A ^ { T }$ are given in:

J. Varah (1979). “On the Separation of Two Matrices,” SIAM J. Numer. Anal. 16, 216–22.

R. Byers and S.G. Nash (1987). “On the Singular Vectors of the Lyapunov Operator,” SIAM J. Alg. Disc. Methods 8, 59–66.
