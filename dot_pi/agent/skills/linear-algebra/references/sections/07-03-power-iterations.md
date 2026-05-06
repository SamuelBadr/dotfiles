# 7.3 Power Iterations

Suppose that we are given $A \in \mathbb { C } ^ { n \times n }$ and a unitary $U _ { 0 } \in \mathbb { C } ^ { n \times n }$ . Recall from §5.2.10 that the Householder QR factorization can be extended to complex matrices and consider the following iteration:

$$
T _ {0} = U _ {0} ^ {H} A U _ {0}
$$

for k = 1, 2, . . .

$$
T _ {k - 1} = U _ {k} R _ {k} \quad (\text { QR   factorization }) \tag {7.3.1}
$$

$$
T _ {k} = R _ {k} U _ {k}
$$

end

Since $T _ { k } = R _ { k } U _ { k } = U _ { k } ^ { H } ( U _ { k } R _ { k } ) U _ { k } = U _ { k } ^ { H } T _ { k - 1 } U _ { k }$ it follows by induction that

$$
T _ {k} = (U _ {0} U _ {1} \dots U _ {k}) ^ {H} A (U _ {0} U _ {1} \dots U _ {k}). \tag {7.3.2}
$$

Thus, each $T _ { k }$ is unitarily similar to A. Not so obvious, and what is a central theme of this section, is that the $T _ { k }$ almost always converge to upper triangular form, i.e., (7.3.2) almost always “converges” to a Schur decomposition of A.

Iteration (7.3.1) is called the QR iteration, and it forms the backbone of the most effective algorithm for computing a complete Schur decomposition of a dense general matrix. In order to motivate the method and to derive its convergence properties, two other eigenvalue iterations that are important in their own right are presented first: the power method and the method of orthogonal iteration.

# 7.3.1 The Power Method

Suppose $A \in \mathbb { C } ^ { n \times n }$ and $X ^ { - 1 } A X = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ with $X = \left[ \left. x _ { 1 } \right| \cdot \cdot \cdot \right| \left. x _ { n } \right]$ . Assume that

$$
\left| \lambda_ {1} \right| > \left| \lambda_ {2} \right| \geq \dots \geq \left| \lambda_ {n} \right|.
$$

Given a unit 2-norm $q ^ { ( 0 ) } \in \mathbb { C } ^ { n }$ , the power method produces a sequence of vectors $q ^ { ( k ) }$ as follows:

for k = 1, 2, . . .

$$
z ^ {(k)} = A q ^ {(k - 1)}
$$

$$
q ^ {(k)} = z ^ {(k)} / \| z ^ {(k)} \| _ {2} \tag {7.3.3}
$$

$$
\lambda^ {(k)} = [ q ^ {(k)} ] ^ {H} A q ^ {(k)}
$$

end

There is nothing special about using the 2-norm for normalization except that it imparts a greater unity on the overall discussion in this section.

Let us examine the convergence properties of the power iteration. If

$$
q ^ {(0)} = a _ {1} x _ {1} + a _ {2} x _ {2} + \dots + a _ {n} x _ {n} \tag {7.3.4}
$$

and $a _ { 1 } \neq 0$ , then

$$
A ^ {k} q ^ {(0)} = a _ {1} \lambda_ {1} ^ {k} \left(x _ {1} + \sum_ {j = 2} ^ {n} \frac {a _ {j}}{a _ {1}} \left(\frac {\lambda_ {j}}{\lambda_ {1}}\right) ^ {k} x _ {j}\right).
$$

Since $q ^ { ( k ) } \in \mathsf { s p a n } \{ A ^ { k } q ^ { ( 0 ) } \}$ we conclude that

$$
\operatorname{dist} \left(\operatorname{span} \{q ^ {(k)} \}, \operatorname{span} \{x _ {1} \}\right) = O \left(\left| \frac {\lambda_ {2}}{\lambda_ {1}} \right| ^ {k}\right).
$$

It is also easy to verify that

$$
\left| \lambda_ {1} - \lambda^ {(k)} \right| = O \left(\left| \frac {\lambda_ {2}}{\lambda_ {1}} \right| ^ {k}\right). \tag {7.3.5}
$$

Since $\lambda _ { 1 }$ is larger than all the other eigenvalues in modulus, it is referred to as a dominant eigenvalue. Thus, the power method converges if $\lambda _ { 1 }$ is dominant and if $q ^ { ( 0 ) }$ has a component in the direction of the corresponding dominant eigenvector $x _ { 1 }$ . The behavior of the iteration without these assumptions is discussed in Wilkinson (AEP, p. 570) and Parlett and Poole (1973).

In practice, the usefulness of the power method depends upon the ratio $| \lambda _ { 2 } | / | \lambda _ { 1 } |$ , since it dictates the rate of convergence. The danger that $q ^ { ( 0 ) }$ is deficient in $x _ { 1 }$ is less worrisome because rounding errors sustained during the iteration typically ensure that subsequent iterates have a component in this direction. Moreover, it is typically the case in applications that one has a reasonably good guess as to the direction of $x _ { 1 }$ . This guards against having a pathologically small coefficient $a _ { 1 }$ in (7.3.4).

Note that the only thing required to implement the power method is a procedure for matrix-vector products. It is not necessary to store A in an n-by-n array. For this reason, the algorithm is of interest when the dominant eigenpair for a large sparse matrix is required. We have much more to say about large sparse eigenvalue problems in Chapter 10.

Estimates for the error $| \lambda ^ { ( k ) } - \lambda _ { 1 } |$ can be obtained by applying the perturbation theory developed in §7.2.2. Define the vector

$$
r ^ {(k)} = A q ^ {(k)} - \lambda^ {(k)} q ^ {(k)}
$$

and observe that $( A + E ^ { ( k ) } ) q ^ { ( k ) } \ = \ \lambda ^ { ( k ) } q ^ { ( k ) }$ where $E ^ { ( k ) } = - r ^ { ( k ) } [ q ^ { ( k ) } ] ^ { H }$ . Thus $\lambda ^ { ( k ) }$ is an eigenvalue of $A + E ^ { ( k ) }$ and

$$
| \lambda^ {(k)} - \lambda_ {1} | \approx \frac {\parallel E ^ {(k)} \parallel_ {2}}{s (\lambda_ {1})} = \frac {\parallel r ^ {(k)} \parallel_ {2}}{s (\lambda_ {1})}.
$$

If we use the power method to generate approximate right and left dominant eigenvectors, then it is possible to obtain an estimate of $s ( \lambda _ { 1 } )$ . In particular, if $w ^ { ( k ) }$ is a unit 2-norm vector in the direction of $( A ^ { H } ) ^ { k } w ^ { ( 0 ) }$ , then we can use the approximation s(λ1) ≈ | w(k)H q $s ( \lambda _ { 1 } ) \approx | \boldsymbol { w } ^ { ( k ) } \boldsymbol { q } ^ { \prime } ( k ) |$ .

# 7.3.2 Orthogonal Iteration

A straightforward generalization of the power method can be used to compute higherdimensional invariant subspaces. Let r be a chosen integer satisfying $1 \leq r \leq n$ . Given $A \in \mathbb { C } ^ { n \times n }$ and an n-by-r matrix $Q _ { 0 }$ with orthonormal columns, the method of orthogonal iteration generates a sequence of matrices $\{ Q _ { k } \} \subseteq \mathbb { C } ^ { n \times r }$ and a sequence of eigenvalue estimates $\big \{ \lambda _ { 1 } ^ { ( k ) } , \dots , \lambda _ { r } ^ { ( k ) } \big \}$ as follows:

for $k = 1 , 2 , \dots$

$$
Z _ {k} = A Q _ {k - 1}
$$

$$
Q _ {k} R _ {k} = Z _ {k} \quad \text {(QR factorization)} \tag {7.3.6}
$$

$$
\lambda (Q _ {k} ^ {H} A Q _ {k}) = \{\lambda_ {1} ^ {(k)}, \dots , \lambda_ {r} ^ {(k)} \}
$$

end

Note that if $r \ = \ 1$ , then this is just the power method (7.3.3). Moreover, the sequence $\{ Q _ { k } e _ { 1 } \}$ is precisely the sequence of vectors produced by the power iteration with starting vector $q ^ { ( 0 ) } = Q _ { 0 } e _ { 1 }$ .

In order to analyze the behavior of this iteration, suppose that

$$
Q ^ {H} A Q = T = \operatorname{diag} \left(\lambda_ {i}\right) + N, \quad \left| \lambda_ {1} \right| \geq \left| \lambda_ {2} \right| \geq \dots \geq \left| \lambda_ {n} \right| \tag {7.3.7}
$$

is a Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ . Assume that $1 \leq r < n$ and partition Q and $T$ as follows:

$$
Q = \left[ \begin{array}{c c} Q _ {\alpha} & Q _ {\beta} \\ r & n - r \end{array} \right], \quad T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \\ r & n - r \end{array} \right] _ {n - r} ^ {r}. \tag {7.3.8}
$$

If $\left| \lambda _ { r } \right| > \left| \lambda _ { r + 1 } \right|$ , then the subspace $D _ { r } ( A ) \ = \ r { \mathsf { a n } } ( Q _ { \alpha } )$ is referred to as a dominant invariant subspace. It is the unique invariant subspace associated with the eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { r }$ . The following theorem shows that with reasonable assumptions, the subspaces $\mathsf { r a n } ( Q _ { k } )$ generated by (7.3.6) converge to $D _ { r } ( A )$ at a rate proportional to $| \lambda _ { r + 1 } / \lambda _ { r } | ^ { k }$ .

Theorem 7.3.1. Let the Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ be given by (7.3.7) and (7.3.8) with $n \geq 2$ . Assume that $| \lambda _ { r } | > | \lambda _ { r + 1 } |$ and that $\mu \geq 0$ satisfies

$$
(1 + \mu) | \lambda_ {r} | > \| N \| _ {F}.
$$

Suppose $Q _ { 0 } \in \mathbb { C } ^ { n \times r }$ has orthonormal columns and that $d _ { k }$ is defined by

$$
d _ {k} = \operatorname{dist} \left(D _ {r} (A), \operatorname{ran} \left(Q _ {k}\right)\right), \quad k \geq 0.
$$

If

$$
d _ {0} <   1, \tag {7.3.9}
$$

then the matrices $Q _ { k }$ generated by (7.3.6) satisfy

$$
d _ {k} \leq (1 + \mu) ^ {n - 2} \cdot \left(1 + \frac {\| T _ {1 2} \| _ {F}}{\operatorname{sep} \left(T _ {1 1} , T _ {2 2}\right)}\right) \cdot \left[ \frac {\left| \lambda_ {r + 1} \right| + \frac {\| N \| _ {F}}{1 + \mu}}{\left| \lambda_ {r} \right| - \frac {\| N \| _ {F}}{1 + \mu}} \right] ^ {k} \cdot \frac {d _ {0}}{\sqrt {1 - d _ {0} ^ {2}}}. \tag {7.3.10}
$$

Proof. The proof is given in an appendix at the end of this section.

The condition (7.3.9) ensures that the initial matrix $Q _ { 0 }$ is not deficient in certain eigendirections. In particular, no vector in the span of $Q _ { 0 } \mathrm { { ^ { \circ } s } }$ columns is orthogonal to $\bar { D _ { r } } ( A ^ { H } )$ . The theorem essentially says that if this condition holds and if $\mu$ is chosen large enough, then

$$
\operatorname{dist} \left(D _ {r} (A), \operatorname{ran} \left(Q _ {k}\right)\right) \approx c \left| \frac {\lambda_ {r + 1}}{\lambda_ {r}} \right| ^ {k}
$$

where c depends on sep $( T _ { 1 1 } , T _ { 2 2 } )$ and A’s departure from normality.

It is possible to accelerate the convergence in orthogonal iteration using a technique described in Stewart (1976). In the accelerated scheme, the approximate eigenvalue ${ \lambda } _ { i } ^ { ( k ) }$ satisfies

$$
\left| \lambda_ {i} ^ {(k)} - \lambda_ {i} \right| \approx \left| \frac {\lambda_ {r + 1}}{\lambda_ {i}} \right| ^ {k}, \quad i = 1: r.
$$

(Without the acceleration, the right-hand side is $| \lambda _ { i + 1 } / \lambda _ { i } | ^ { k } . \big )$ Stewart’s algorithm involves computing the Schur decomposition of the matrices $Q _ { k } ^ { T } A Q _ { k }$ every so often. The method can be very useful in situations where A is large and sparse and a few of its largest eigenvalues are required.

# 7.3.3 The QR Iteration

We now derive the QR iteration (7.3.1) and examine its convergence. Suppose $r = n$ in (7.3.6) and the eigenvalues of A satisfy

$$
\left| \lambda_ {1} \right| > \left| \lambda_ {2} \right| > \dots > \left| \lambda_ {n} \right|.
$$

Partition the matrix Q in (7.3.7) and $Q _ { k }$ in (7.3.6) as follows:

$$
Q = \left[ q _ {1} \mid \dots \mid q _ {n} \right], \quad Q _ {k} = \left[ q _ {1} ^ {(k)} \mid \dots \mid q _ {n} ^ {(k)} \right].
$$

If

$$
\mathsf {d i s t} (D _ {i} (A ^ {H}), \mathsf {s p a n} \{q _ {1} ^ {(0)}, \dots , q _ {i} ^ {(0)} \}) <   1, \quad i = 1: n, \tag {7.3.11}
$$

then it follows from Theorem 7.3.1 that

$$
\operatorname{dist} \left(\operatorname{span} \left\{q _ {1} ^ {(k)}, \dots , q _ {i} ^ {(k)} \right\}, \operatorname{span} \left\{q _ {1}, \dots , q _ {i} \right\}\right)\rightarrow 0
$$

for i = 1:n. This implies that the matrices $T _ { k }$ defined by

$$
T _ {k} = Q _ {k} ^ {H} A Q _ {k}
$$

are converging to upper triangular form. Thus, it can be said that the method of orthogonal iteration computes a Schur decomposition provided the original iterate $Q _ { 0 } \in \mathbb { C } ^ { n \times n }$ is not deficient in the sense of (7.3.11).

The QR iteration arises naturally by considering how to compute the matrix $T _ { k }$ directly from its predecessor $T _ { k - 1 }$ . On the one hand, we have from (7.3.6) and the definition of $T _ { k - 1 }$ that

$$
T _ {k - 1} = Q _ {k - 1} ^ {H} A Q _ {k - 1} = Q _ {k - 1} ^ {H} (A Q _ {k - 1}) = (Q _ {k - 1} ^ {H} Q _ {k}) R _ {k}.
$$

On the other hand,

$$
T _ {k} = Q _ {k} ^ {H} A Q _ {k} = (Q _ {k} ^ {H} A Q _ {k - 1}) (Q _ {k - 1} ^ {H} Q _ {k}) = R _ {k} (Q _ {k - 1} ^ {H} Q _ {k}).
$$

Thus, $T _ { k }$ is determined by computing the QR factorization of $T _ { k - 1 }$ and then multiplying the factors together in reverse order, precisely what is done in (7.3.1).

Note that a single QR iteration is an $O ( n ^ { 3 } )$ calculation. Moreover, since convergence is only linear (when it exists), it is clear that the method is a prohibitively expensive way to compute Schur decompositions. Fortunately these practical difficulties can be overcome as we show in §7.4 and §7.5.

# 7.3.4 LR Iterations

We conclude with some remarks about power iterations that rely on the LU factorization rather than the QR factorizaton. Let $G _ { 0 } \in \mathbb { C } ^ { n \times r }$ have rank r. Corresponding to (7.3.1) we have the following iteration:

for $k = 1 , 2 , \dots$

$$
Z _ {k} = A G _ {k - 1} \tag {7.3.12}
$$

$$
Z _ {k} = G _ {k} R _ {k} \quad \text {(LU factorization)}
$$

end

Suppose $r = n$ and that we define the matrices $T _ { k }$ by

$$
T _ {k} = G _ {k} ^ {- 1} A G _ {k}. \tag {7.3.13}
$$

It can be shown that if we set $L _ { 0 } = G _ { 0 }$ , then the $T _ { k }$ can be generated as follows:

$$
T _ {0} = L _ {0} ^ {- 1} A L _ {0}
$$

for $k = 1 , 2 ,$ . . .

$$
T _ {k - 1} = L _ {k} R _ {k} \quad (\text {LU factorization}) \tag {7.3.14}
$$

$$
T _ {k} = R _ {k} L _ {k}
$$

end

Iterations (7.3.12) and (7.3.14) are known as treppeniteration and the LR iteration, respectively. Under reasonable assumptions, the $T _ { k }$ converge to upper triangular form. To successfully implement either method, it is necessary to pivot. See Wilkinson (AEP, p. 602).

# Appendix

In order to establish Theorem 7.3.1 we need the following lemma that bounds powers of a matrix and powers of its inverse.

Lemma 7.3.2. Let $Q ^ { H } A Q = T = D + N $ be a Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ where D is diagonal and N strictly upper triangular. Let $\lambda _ { \mathrm { m a x } }$ and $\lambda _ { \mathrm { m i n } }$ denote the largest and smallest eigenvalues of A in absolute value. If $\mu \geq 0$ , then for all $k \geq 0$ we have

$$
\| A ^ {k} \| _ {2} \leq (1 + \mu) ^ {n - 1} \left(\left| \lambda_ {\max} \right| + \frac {\left\| N \right\| _ {F}}{1 + \mu}\right) ^ {k}. \tag {7.3.15}
$$

If A is nonsingular and $\mu \geq 0$ satisfies $( 1 + \mu ) | \lambda _ { \operatorname* { m i n } } | > \parallel N \parallel _ { _ F }$ , then for all $k \geq 0$ we also have

$$
\| A ^ {- k} \| _ {2} \leq (1 + \mu) ^ {n - 1} \left(\frac {1}{\left| \lambda_ {\min} \right| - \| N \| _ {F} / (1 + \mu)}\right) ^ {k}. \tag {7.3.16}
$$

Proof. For $\mu \geq 0$ , define the diagonal matrix $\Delta$ by

$$
\Delta = \operatorname{diag} (1, (1 + \mu), (1 + \mu) ^ {2}, \dots , (1 + \mu) ^ {n - 1})
$$

and note that $\kappa _ { 2 } ( \Delta ) = ( 1 + \mu ) ^ { n - 1 }$ . Since N is strictly upper triangular, it is easy to verify that

$$
\| \Delta N \Delta^ {- 1} \| _ {F} \leq \frac {\| N \| _ {F}}{1 + \mu}
$$

and thus

$$
\begin{array}{l} \| A ^ {k} \| _ {2} = \| T ^ {k} \| _ {2} = \| \Delta^ {- 1} (D + \Delta N \Delta^ {- 1}) ^ {k} \Delta \| _ {2} \\ \leq \kappa_ {2} (\Delta) \left(\| D \| _ {2} + \| \Delta N \Delta^ {- 1} \| _ {2}\right) ^ {k} \leq (1 + \mu) ^ {n - 1} \left(| \lambda_ {\max} | + \frac {\| N \| _ {F}}{1 + \mu}\right) ^ {k}. \\ \end{array}
$$

On the other hand, if A is nonsingular and $( 1 + \mu ) | \lambda _ { \operatorname* { m i n } } | > \parallel N \parallel _ { F }$ , then

$$
\| \Delta D ^ {- 1} N \Delta^ {- 1} \| _ {2} = \| D ^ {- 1} (\Delta N \Delta^ {- 1}) \| _ {2} \leq \frac {1}{| \lambda_ {\min} |} \| \Delta N \Delta^ {- 1} \| _ {F} <   1.
$$

Using Lemma 2.3.3 we obtain

$$
\begin{array}{l} \| A ^ {- k} \| _ {2} = \| T ^ {- k} \| _ {2} = \left\| \Delta^ {- 1} [ (I + \Delta D ^ {- 1} N \Delta^ {- 1}) ^ {- 1} D ^ {- 1} ] ^ {k} \Delta \right\| _ {2} \\ \leq \kappa_ {2} (\Delta) \left(\frac {\| D ^ {- 1} \| _ {2}}{1 - \| \Delta D ^ {- 1} N \Delta^ {- 1} \| _ {2}}\right) ^ {k} \leq (1 + \mu) ^ {n - 1} \left(\frac {1}{| \mu | - \| N \| _ {F} / (1 + \mu)}\right) ^ {k} \\ \end{array}
$$

completing the proof of the lemma.

Proof of Theorem 7.3.1. By induction it is easy to show that the matrix $Q _ { k }$ in (7.3.6) satisfies

$$
A ^ {k} Q _ {0} = Q _ {k} (R _ {k} \dots R _ {1}),
$$

a QR factorization of $A ^ { k } Q _ { 0 }$ . By substituting the Schur decomposition (7.3.7)-(7.3.8) into this equation we obtain

$$
T ^ {k} \left[ \begin{array}{c} V _ {0} \\ W _ {0} \end{array} \right] = \left[ \begin{array}{c} V _ {k} \\ W _ {k} \end{array} \right] (R _ {k} \dots R _ {1}) \tag {7.3.17}
$$

where

$$
V _ {k} = Q _ {\alpha} ^ {H} Q _ {k}, \qquad W _ {k} = Q _ {\beta} ^ {H} Q _ {k}.
$$

Our goal is to bound $\parallel W _ { k } \parallel _ { 2 }$ since by the definition of subspace distance given in §2.5.3 we have

$$
\| W _ {k} \| _ {2} = \operatorname{dist} (D _ {r} (A), \operatorname{ran} (Q _ {k})). \tag {7.3.18}
$$

Note from the thin CS decomposition (Theorem 2.5.2) that

$$
1 = d _ {k} ^ {2} + \sigma_ {\min} (V _ {k}) ^ {2}. \tag {7.3.19}
$$

Since $T _ { 1 1 }$ and $T _ { 2 2 }$ have no eigenvalues in common, Lemma 7.1.5 tells us that the Sylvester equation $T _ { 1 1 } X \ : - \ : X T _ { 2 2 } \ : = \ : - T _ { 1 2 }$ has a solution $X \in \mathbb { C } ^ { r \times ( n - r ) }$ and that

$$
\| X \| _ {F} \leq \frac {\| T _ {1 2} \| _ {F}}{\operatorname{sep} \left(T _ {1 1} , T _ {2 2}\right)}. \tag {7.3.20}
$$

It follows that

$$
\left[ \begin{array}{c c} I _ {r} & X \\ 0 & I _ {n - r} \end{array} \right] ^ {- 1} \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] \left[ \begin{array}{c c} I _ {r} & X \\ 0 & I _ {n - r} \end{array} \right] = \left[ \begin{array}{c c} T _ {1 1} & 0 \\ 0 & T _ {2 2} \end{array} \right].
$$

By substituting this into (7.3.17) we obtain

$$
\left[ \begin{array}{c c} T _ {1 1} ^ {k} & 0 \\ 0 & T _ {2 2} ^ {k} \end{array} \right] \left[ \begin{array}{c} V _ {0} - X W _ {0} \\ W _ {0} \end{array} \right] = \left[ \begin{array}{c} V _ {k} - X W _ {k} \\ W _ {k} \end{array} \right] (R _ {k} \dots R _ {1}),
$$

$$
T _ {1 1} ^ {k} (V _ {0} - X W _ {0}) = (V _ {k} - X W _ {k}) (R _ {k} \dots R _ {1}), \tag {7.3.21}
$$

$$
T _ {2 2} ^ {k} W _ {0} = W _ {k} (R _ {k} \dots R _ {1}). \tag {7.3.22}
$$

The matrix $I + X X ^ { H }$ is Hermitian positive definite and so it has a Cholesky factorization

$$
I + X X ^ {H} = G G ^ {H}. \tag {7.3.23}
$$

It is clear that

$$
\sigma_ {\min} (G) \geq 1. \tag {7.3.24}
$$

If the matrix $Z \in \mathbb { C } ^ { n \times ( n - r ) }$ is defined by

$$
Z = Q \left[ \begin{array}{c} I _ {r} \\ - X ^ {H} \end{array} \right] G ^ {- H} = \left[ Q _ {\alpha} Q _ {\beta} \right] \left[ \begin{array}{c} I _ {r} \\ - X ^ {H} \end{array} \right] G ^ {- H} = (Q _ {\alpha} - Q _ {\beta} X ^ {H}) G ^ {- H},
$$

then it follows from the equation $A ^ { H } Q = Q T ^ { H }$ that

$$
A ^ {H} (Q _ {\alpha} - Q _ {\beta} X ^ {H}) = (Q _ {\alpha} - Q _ {\beta} X ^ {H}) T _ {1 1} ^ {H}. \tag {7.3.25}
$$

Since $Z ^ { H } Z = I _ { r }$ and ran $( Z ) \ : = \ : \mathsf { r a n } ( Q _ { \alpha } - Q _ { \beta } X ^ { H } )$ , it follows that the columns of Z are an orthonormal basis for $D _ { r } ( A ^ { H } )$ . Using the CS decomposition, (7.3.19), and the fact that ran $( Q _ { \beta } ) = D _ { r } ( A ^ { H } ) ^ { \perp }$ , we have

$$
\begin{array}{l} \sigma_ {\min} (Z ^ {T} Q _ {0}) ^ {2} = 1 - \operatorname{dist} (D _ {r} (A ^ {H}), Q _ {0}) ^ {2} = 1 - \| Q _ {\beta} ^ {H} Q _ {0} \| \\ = \sigma_ {\mathrm{min}} (Q _ {\alpha} ^ {T} Q _ {0}) ^ {2} = \sigma_ {\mathrm{min}} (V _ {0}) ^ {2} = 1 - d _ {0} ^ {2} > 0. \\ \end{array}
$$

This shows that

$$
V _ {0} - X W _ {0} = \left[ \begin{array}{l} I _ {r} \end{array} \right| - X \left. \right] \left[ \begin{array}{l} Q _ {\alpha} ^ {H} Q _ {0} \\ Q _ {\beta} ^ {H} Q _ {0} \end{array} \right] = (Z G ^ {H}) ^ {H} Q _ {0} = G (Z ^ {H} Q _ {0})
$$

is nonsingular and together with (7.3.24) we obtain

$$
\| (V _ {0} - X W _ {0}) ^ {- 1} \| _ {2} \leq \| G ^ {- 1} \| _ {2} \| (Z ^ {H} Q _ {0}) ^ {- 1} \| _ {2} \leq \frac {1}{\sqrt {1 - d _ {0} ^ {2}}}. \tag {7.3.26}
$$

Manipulation of (7.3.19) and (7.3.20) yields

$$
W _ {k} = T _ {2 2} ^ {k} W _ {0} (R _ {k} \dots R _ {1}) ^ {- 1} = T _ {2 2} ^ {k} W _ {0} (V _ {0} - X W _ {0}) ^ {- 1} T _ {1 1} ^ {- k} (V _ {k} - X W _ {k}).
$$

The verification of (7.3.10) is completed by taking norms in this equation and using (7.3.18), (7.3.19), (7.3.20), (7.3.26), and the following facts:

$$
\| T _ {2 2} ^ {k} \| _ {2} \leq (1 + \mu) ^ {n - r - 1} \left(| \lambda_ {r + 1} | + \| N \| _ {F} / (1 + \mu)\right) ^ {k},
$$

$$
\left\| T _ {1 1} ^ {- k} \right\| _ {2} \leq (1 + \mu) ^ {r - 1} / \left(\left| \lambda_ {r} \right| - \left\| N \right\| _ {F} / (1 + \mu)\right) ^ {k},
$$

$$
\| V _ {k} - X W _ {k} \| _ {2} \leq \| V _ {k} \| _ {2} + \| X \| _ {2} \| W _ {k} \| _ {2} \leq 1 + \| T _ {1 2} \| _ {F} / \mathsf {s e p} (T _ {1 1}, T _ {2 2}).
$$

The bounds for $\| \ : T _ { 2 2 } ^ { k } \ : \| _ { 2 }$ and $\Vert \ T _ { 1 1 } ^ { - k } \ \Vert _ { 2 }$ follow from Lemma 7.3.2.

# Problems

# P7.3.1 Verify Equation (7.3.5).

P7.3.2 Suppose the eigenvalues of $A \in \mathbb { R } ^ { n \times n }$ satisfy $| \lambda _ { 1 } | = | \lambda _ { 2 } | > | \lambda _ { 3 } | \geq \cdot \cdot \cdot \geq | \lambda _ { n } |$ and that $\lambda _ { 1 }$ and λ2 are complex conjugates of one another. Let $S = \mathsf { s p a n } \{ y , z \}$ where $y , z \in \mathbb { R } ^ { n }$ satisfy $A ( y + i z ) =$ $\lambda _ { 1 } ( y + i z )$ . Show how the power method with a real starting vector can be used to compute an approximate basis for S.

P7.3.3 Assume $A \in \mathbb { R } ^ { n \times n }$ has eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { n }$ that satisfy

$$
\lambda = \lambda_ {1} = \lambda_ {2} = \lambda_ {3} = \lambda_ {4} > | \lambda_ {5} | \geq \dots \geq | \lambda_ {n} |
$$

where λ is positive. Assume that A has two Jordan blocks of the form.

$$
\left[ \begin{array}{c c} \lambda & 1 \\ 0 & \lambda \end{array} \right].
$$

Discuss the convergence properties of the power method when applied to this matrix and how the convergence might be accelerated.

P7.3.4 A matrix A is a positive matrix if $a _ { i j } ~ > ~ 0$ for all i and j. A vector $v \in \mathbb { R } ^ { n }$ is a positive vector if $v _ { i } > 0$ for all i. Perron’s theorem states that if A is a positive square matrix, then it has a unique dominant eigenvalue equal to its spectral radius $\rho ( A )$ and there is a positive vector x so that $A x = \rho ( A ) \cdot x$ . In this context, x is called the Perron vector and $\rho ( A )$ is called the Perron root. Assume that $A \in \mathbb { R } ^ { n \times n }$ is positive and $q \in \mathbb { R } ^ { n }$ is positive with unit 2-norm. Consider the following implementation of the power method (7.3.3):

$$
\begin{array}{l} z = A q, \lambda = q ^ {T} z \\ q = z, q = q / \| q \| _ {2}, z = A q, \lambda = q ^ {T} z \\ \end{array}
$$

(a) Adjust the termination criteria to guarantee (in principle) that the final λ and q satisfy $\tilde { A } q = \lambda q$ , where $\begin{array} { r } { \parallel \tilde { A } - A \parallel _ { 2 } \leq \delta } \end{array}$ and A˜ is positive. (b) Applied to a positive matrix $A \in \mathbb { R } ^ { n \times n }$ , the Collatz-Wielandt formula states that $\rho ( A )$ is the maximum value of the function f defined by

$$
f (x) = \min _ {1 \leq i \leq n} \frac {y _ {i}}{x _ {i}}
$$

where $\boldsymbol { x } \in \mathbb { R } ^ { n }$ is positive and $y = A x$ . Does it follow that $f ( A q ) \geq f ( q ) ?$ In other words, do the iterates $\{ q ^ { ( k ) } \}$ in the power method have the property that $f ( q ^ { ( k ) } )$ increases monotonically to the Perron root, assuming that $q ^ { ( 0 ) }$ is positive?

P7.3.5 (Read the previous problem for background.) A matrix A is a nonnegative matrix if $a _ { i j } \geq 0$ for all i and j. A matrix $A \in \mathbb { R } ^ { n \times n }$ is reducible if there is a permutation P so that $P ^ { T } A P$ is block triangular with two or more square diagonal blocks. A matrix that is not reducible is irreducible. The Perron-Frobenius theorem states that if A is a square, nonnegative, and irreducible, then $\rho ( A )$ , the Perron root, is an eigenvalue for A and there is a positive vector x, the Perron vector, so that $A x = \rho ( A ) { \cdot } x$ . Assume that $A _ { 1 } , A _ { 2 } , A _ { 3 } \in \mathbb { R } ^ { n \times n }$ are each positive and let the nonnegative matrix A be defined by

$$
A = \left[ \begin{array}{c c c} 0 & A _ {1} & 0 \\ 0 & 0 & A _ {2} \\ A _ {3} & 0 & 0 \end{array} \right].
$$

(a) Show that A is irreducible. (b) Let $B = A _ { 1 } A _ { 2 } A _ { 3 }$ . Show how to compute the Perron root and vector for A from the Perron root and vector for B. (c) Show that A has other eigenvalues with absolute value equal to the Perron root. How could those eigenvalues and the associated eigenvectors be computed?

P7.3.6 (Read the previous two problems for background.) A nonnegative matrix $P \in \mathbb { R } ^ { n \times n }$ is stochastic if the entries in each column sum to 1. A vector $v \in \mathbb { R } ^ { n }$ is a probability vector if its entries are nonnegative and sum to 1. (a) Show that if $P \in \mathbb { R } ^ { n \times n }$ is stochastic and $v \in \mathbb { R } ^ { n }$ is a probability vector, then $w = P v$ is also a probability vector. (b) The entries in a stochastic matrix $P \in \mathbb { R } ^ { n \times n }$ can be regarded as the transition probabilities associated with an n-state Markov Chain. Let $v _ { j }$ be the probability of being in state $j$ at time $t = t _ { \mathrm { c u r r e n t } }$ . In the Markov model, the probability of being in state i at time $t = t _ { \mathrm { n e x t } }$ is given by

$$
w _ {i} = \sum_ {j = 1} ^ {n} p _ {i j} v _ {j} \quad i = 1: n,
$$

$\mathrm { i } . \mathrm { e } . , w = P v$ . With the help of a biased coin, a surfer on the World Wide Web randomly jumps from page to page. Assume that the surfer is currently viewing web page $j$ and that the coin comes up heads with probability α. Here is how the surfer determines the next page to visit:

Step 1. A coin is tossed.

Step 2. If it comes up heads and web page j has at least one outlink, then the next page to visit is randomly selected from the list of outlink pages.

Step 3. Otherwise, the next page to visit is randomly selected from the list of all possible pages.

Let $P \in \mathbb { R } ^ { n \times n }$ be the matrix of transition probabilities that define this random process. Specify P in terms of $\alpha ,$ , the vector of ones $e ,$ and the link matrix $H \in \mathbb { R } ^ { n \times n }$ defined by

$$
h _ {i j} = \left\{ \begin{array}{l l} 1 & \text { if   there   is   a   link   on   web   page   } j \text {   to   web   page   } i \\ 0 & \text { otherwise } \end{array} \right.
$$

Hints: The number of nonzero components in $H ( : , j )$ is the number of outlinks on web page $j , P$ is a convex combination of a very sparse sparse matrix and a very dense rank-1 matrix. (c) Detail how the power method can be used to determine a probability vector x so that $P x = x$ . Strive to get as much computation “outside the loop” as possible. Note that in the limit we can expect to find the random surfer viewing web page i with probability $x _ { i }$ . Thus, a case can be made that more important pages are associated with the larger components of x. This is the basis of Google PageRank. If

$$
x _ {i _ {1}} \geq x _ {i _ {2}} \geq \dots \geq x _ {i _ {n}}
$$

then web page $i _ { k }$ has page rank k.

P7.3.7 (a) Show that if $\ b X \in \mathbb { C } ^ { n \times n }$ is nonsingular, then

$$
\parallel A \parallel_ {X} = \parallel X ^ {- 1} A X \parallel_ {2}
$$

defines a matrix norm with the property that

$$
\| A B \| _ {X} \leq \| A \| _ {X} \| B \| _ {X}.
$$

(b) Show that for any $\epsilon > 0$ there exists a nonsingular $\ b X \in \mathbb { C } ^ { n \times n }$ such that

$$
\parallel A \parallel_ {X} = \parallel X ^ {- 1} A X \parallel_ {2} \leq \rho (A) + \epsilon
$$

where $\rho ( A )$ is A’s spectral radius. Conclude that there is a constant M such that

$$
\| A ^ {k} \| _ {2} \leq M (\rho (A) + \epsilon) ^ {k}
$$

for all non-negative integers k. (Hint: Set $X = Q \ \mathrm { d i a g } ( 1 , a , . . . , a ^ { n - 1 } )$ where $Q ^ { H } A Q = D + N$ is A’s Schur decomposition.)

P7.3.8 Verify that (7.3.14) calculates the matrices $T _ { k }$ defined by (7.3.13).

P7.3.9 Suppose $A \in \mathbb { C } ^ { n \times n }$ is nonsingular and that $Q _ { 0 } \in \mathbb { C } ^ { n \times p }$ has orthonormal columns. The following iteration is referred to as inverse orthogonal iteration.

$$
\begin{array}{l} \text { for } k = 1, 2, \dots \\ \text { Solve } A Z _ {k} = Q _ {k - 1} \text { for } Z _ {k} \in \mathbb {C} ^ {n \times p} \\ Z _ {k} = Q _ {k} R _ {k} \quad \text {(QR factorization)} \\ \end{array}
$$

end

Explain why this iteration can usually be used to compute the p smallest eigenvalues of A in absolute value. Note that to implement this iteration it is necessary to be able to solve linear systems that involve A. If $p = 1$ , the method is referred to as the inverse power method.

# Notes and References for §7.3

For an excellent overview of the QR iteration and related procedures, see Watkins (MEP), Stewart (MAE), and Kressner (NMSE). A detailed, practical discussion of the power method is given in Wilkinson (AEP, Chap. 10). Methods are discussed for accelerating the basic iteration, for calculating nondominant eigenvalues, and for handling complex conjugate eigenvalue pairs. The connections among the various power iterations are discussed in:

B.N. Parlett and W.G. Poole (1973). “A Geometric Theory for the QR, LU, and Power Iterations,” SIAM J. Numer. Anal. 10, 389–412.

The QR iteration was concurrently developed in:

J.G.F. Francis (1961). “The QR Transformation: A Unitary Analogue to the LR Transformation,” Comput. J. 4, 265–71, 332–334.

V.N. Kublanovskaya (1961). “On Some Algorithms for the Solution of the Complete Eigenvalue Problem,” USSR Comput. Math. Phys. 3, 637–657.

As can be deduced from the title of the first paper by Francis, the LR iteration predates the QR iteration. The former very fundamental algorithm was proposed by:

H. Rutishauser (1958). “Solution of Eigenvalue Problems with the LR Transformation,” Nat. Bur. Stand. Appl. Math. Ser. 49, 47–81.

More recent, related work includes:

B.N. Parlett (1995). “The New qd Algorithms,” Acta Numerica 5, 459–491.

C. Ferreira and B.N. Parlett (2009). “Convergence of the LR Algorithm for a One-Point Spectrum Tridiagonal Matrix,” Numer. Math. 113, 417–431.

Numerous papers on the convergence and behavior of the QR iteration have appeared, see:

J.H. Wilkinson (1965). “Convergence of the LR, QR, and Related Algorithms,” Comput. J. 8, 77–84.

B.N. Parlett (1965). “Convergence of the Q-R Algorithm,” Numer. Math. 7, 187–93. (Correction in Numer. Math. 10, 163–164.)

B.N. Parlett (1966). “Singular and Invariant Matrices Under the QR Algorithm,” Math. Comput. 20, 611–615.

B.N. Parlett (1968). “Global Convergence of the Basic QR Algorithm on Hessenberg Matrices,” Math. Comput. 22, 803–817.

D.S. Watkins (1982). “Understanding the QR Algorithm,” SIAM Review 24, 427–440.

T. Nanda (1985). “Differential Equations and the QR Algorithm,” SIAM J. Numer. Anal. 22, 310–321.

D.S. Watkins (1993). “Some Perspectives on the Eigenvalue Problem,” SIAM Review 35, 430–471.

D.S. Watkins (2008). “The QR Algorithm Revisited,” SIAM Review 50, 133–145.

D.S. Watkins (2011). “Francis’s Algorithm,” AMS Monthly 118, 387–403.

A block analog of the QR iteration is discussed in:

M. Robb\`e and M. Sadkane (2005). “Convergence Analysis of the Block Householder Block Diagonalization Algorithm,” BIT 45, 181–195.

The following references are concerned with various practical and theoretical aspects of simultaneous iteration:

H. Rutishauser (1970). “Simultaneous Iteration Method for Symmetric Matrices,” Numer. Math. 16, 205–223.

M. Clint and A. Jennings (1971). “A Simultaneous Iteration Method for the Unsymmetric Eigenvalue Problem,” J. Inst. Math. Applic. 8, 111-121.

G.W. Stewart (1976). “Simultaneous Iteration for Computing Invariant Subspaces of Non-Hermitian Matrices,” Numer. Math. 25, 123–136.

A. Jennings (1977). Matrix Computation for Engineers and Scientists, John Wiley and Sons, New York.

Z. Bai and G.W. Stewart (1997). “Algorithm 776: SRRIT: a Fortran Subroutine to Calculate the Dominant Invariant Subspace of a Nonsymmetric Matrix,” ACM Trans. Math. Softw. 23, 494– 513.


---

<!-- golub_400_449 -->

Problems P7.3.4–P7.3.6 explore the relevance of the power method to the problem of computing the Perron root and vector of a nonnegative matrix. For further background and insight, see:

A. Berman and R.J. Plemmons (1994). Nonnegative Matrices in the Mathematical Sciences, SIAM Publications,Philadelphia, PA.

A.N. Langville and C.D. Meyer (2006). Google’s PageRank and Beyond, Princeton University Press, Princeton and Oxford. .

The latter volume is outstanding in how it connects the tools of numerical linear algebra to the design and analysis of Web browsers. See also:

W.J. Stewart (1994). Introduction to the Numerical Solution of Markov Chains, Princeton University Press, Princeton, NJ.

M.W. Berry, Z. Drmaˇc, and E.R. Jessup (1999). “Matrices, Vector Spaces, and Information Retrieval,” SIAM Review 41, 335–362.

A.N. Langville and C.D. Meyer (2005). “A Survey of Eigenvector Methods for Web Information Retrieval,” SIAM Review 47, 135–161.

A.N. Langville and C.D. Meyer (2006). “A Reordering for the PageRank Problem”, SIAM J. Sci. Comput. 27, 2112–2120.

A.N. Langville and C.D. Meyer (2006). “Updating Markov Chains with an Eye on Google’s PageRank,” SIAM J. Matrix Anal. Applic. 27, 968–987.
