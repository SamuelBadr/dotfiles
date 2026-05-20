# 8.1 Properties and Decompositions

In this section we summarize the mathematics required to develop and analyze algorithms for the symmetric eigenvalue problem.

# 8.1.1 Eigenvalues and Eigenvectors

Symmetry guarantees that all of A’s eigenvalues are real and that there is an orthonormal basis of eigenvectors.

Theorem 8.1.1 (Symmetric Schur Decomposition). If $A \in \mathbb { R } ^ { n \times n }$ is symmetric, then there exists a real orthogonal Q such that

$$
Q ^ {T} A Q = \Lambda = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right).
$$

Moreover, for $k = 1 { : } n , A Q ( : , k ) = \lambda _ { k } Q ( : , k )$ . Compare with Theorem 7.1.3.

Proof. Suppose $\lambda _ { 1 } \in \lambda ( A )$ and that $x \in \mathbb { C } ^ { n }$ is a unit 2-norm eigenvector with $A x =$ $\lambda _ { 1 } x$ . Since $\lambda _ { 1 } = x ^ { H } A x = x ^ { H } A ^ { H } x = { \overline { { x ^ { H } A x } } } = { \overline { { \lambda _ { 1 } } } }$ it follows that $\lambda _ { 1 } \in \mathbb { R }$ . Thus, we may assume that $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . Let $P _ { 1 } \in \mathbb { R } ^ { n \times n }$ be a Householder matrix such that $P _ { 1 } ^ { T } x = e _ { 1 } = I _ { n } ( : , 1 )$ . It follows from $A x = \lambda _ { 1 } x$ that $( P _ { 1 } ^ { T } A P _ { 1 } ) e _ { 1 } = \lambda e _ { 1 }$ . This says that the first column of $P _ { 1 } ^ { T } A P _ { 1 }$ is a multiple of $e _ { 1 }$ . But since $P _ { 1 } ^ { T } A P _ { 1 }$ is symmetric, it must have the form

$$
P _ {1} ^ {T} A P _ {1} = \left[ \begin{array}{c c} \lambda_ {1} & 0 \\ 0 & A _ {1} \end{array} \right]
$$

where $A _ { 1 } \in \mathbb { R } ^ { ( n - 1 ) \times ( n - 1 ) }$ is symmetric. By induction we may assume that there is an orthogonal $Q _ { 1 } \in \mathbb { R } ^ { ( n - 1 ) \times ( \overset {  } { n - 1 } ) }$ such that $Q _ { 1 } ^ { T } A _ { 1 } Q _ { 1 } = \Lambda _ { 1 }$ is diagonal. The theorem follows by setting

$$
Q = P _ {1} \left[ \begin{array}{c c} 1 & 0 \\ 0 & Q _ {1} \end{array} \right] \qquad \text {and} \qquad \Lambda = \left[ \begin{array}{c c} \lambda_ {1} & 0 \\ 0 & \Lambda_ {1} \end{array} \right]
$$

and comparing columns in the matrix equation $A Q = Q \Lambda$ .

For a symmetric matrix A we shall use the notation $\lambda _ { k } ( A )$ to designate the kth largest eigenvalue, i.e.,

$$
\lambda_ {n} (A) \leq \dots \leq \lambda_ {2} (A) \leq \lambda_ {1} (A).
$$

It follows from the orthogonal invariance of the 2-norm that A has singular values $\{ | \lambda _ { 1 } ( A ) | , \ldots , | \lambda _ { n } ( A ) | \}$ and

$$
\| A \| _ {2} = \max \{\left| \lambda_ {1} (A) \right|, \left| \lambda_ {n} (A) \right| \}.
$$

The eigenvalues of a symmetric matrix have a minimax characterization that revolves around the quadratic form $x ^ { T } A x / x ^ { T } x$ .

Theorem 8.1.2 (Courant-Fischer Minimax Theorem). If $A \in \mathbb { R } ^ { n \times n }$ is symmetric, then

$$
\lambda_ {k} (A) = \max _ {\dim (S) = k} \min _ {0 \neq y \in S} \frac {y ^ {T} A y}{y ^ {T} y}
$$

for $k = 1 { : } n$ .

Proof. Let $Q ^ { T } A Q \ = \ \mathrm { d i a g } ( \lambda _ { i } )$ be the Schur decomposition with $\lambda _ { k } = \lambda _ { k } ( A )$ and $Q = [  q _ { 1 } | \cdots | q _ { n } ]$ . Define

$$
S _ {k} = \operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\},
$$

the invariant subspace associated with $\lambda _ { 1 } , \ldots , \lambda _ { k }$ . It is easy to show that

$$
\max _ {\dim (S) = k} \min _ {0 \neq y \in S} \frac {y ^ {T} A y}{y ^ {T} y} \geq \min _ {0 \neq y \in S _ {k}} \frac {y ^ {T} A y}{y ^ {T} y} = q _ {k} ^ {T} A q _ {k} = \lambda_ {k} (A).
$$

To establish the reverse inequality, let S be any k-dimensional subspace and note that it must intersect span $\{ q _ { k } , \ldots , q _ { n } \}$ , a subspace that has dimension $n - k + 1$ . If $y _ { * } = \alpha _ { k } q _ { k } + \cdot \cdot \cdot + \alpha _ { n } q _ { n }$ is in this intersection, then

$$
\min _ {0 \neq y \in S} \frac {y ^ {T} A y}{y ^ {T} y} \leq \frac {y _ {*} ^ {T} A y _ {*}}{y _ {*} ^ {T} y _ {*}} \leq \lambda_ {k} (A).
$$

Since this inequality holds for all k-dimensional subspaces,

$$
\max _ {\dim (S) = k} \quad \min _ {0 \neq y \in S} \frac {y ^ {T} A y}{y ^ {T} y} \leq \lambda_ {k} (A)
$$

thereby completing the proof of the theorem.

Note that if $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite, then $\lambda _ { n } ( A ) > 0$ .

# 8.1.2 Eigenvalue Sensitivity

An important solution framework for the symmetric eigenproblem involves the production of a sequence of orthogonal transformations $\{ Q _ { k } \}$ with the property that the matrices $Q _ { k } ^ { T } A Q _ { k }$ are progressively “more diagonal.” The question naturally arises, how well do the diagonal elements of a matrix approximate its eigenvalues?

Theorem 8.1.3 (Gershgorin). Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and that $Q \in \mathbb { R } ^ { n \times n }$ is orthogonal. $I f \dot { Q } ^ { T } A Q = D + \dot { F }$ where $D = \mathrm { d i a g } ( d _ { 1 } , \ldots , d _ { n } )$ and F has zero diagonal entries, then

$$
\lambda (A) \subseteq \bigcup_ {i = 1} ^ {n} \left[ d _ {i} - r _ {i}, d _ {i} + r _ {i} \right]
$$

where $r _ { i } ~ = ~ \sum _ { j = 1 } ^ { n } | f _ { i j } | ~ f o r ~ i = 1 { : } n$ . Compare with Theorem 7.2.1.

Proof. Suppose $\lambda \in \lambda ( A )$ and assume without loss of generality that $\lambda \neq d _ { i }$ for $i = 1 { : } n$ . Since $( D - \lambda I ) + F$ is singular, it follows from Lemma 2.3.3 that

$$
1 \leq \| (D - \lambda I) ^ {- 1} F \| _ {\infty} = \sum_ {j = 1} ^ {n} \frac {| f _ {k j} |}{| d _ {k} - \lambda |} = \frac {r _ {k}}{| d _ {k} - \lambda |}
$$

for some k, $1 \leq k \leq n$ . But this implies that $\lambda \in [ d _ { k } - r _ { k } , d _ { k } + r _ { k } ]$ .

The next results show that if A is perturbed by a symmetric matrix $E _ { \mathrm { { i } } }$ then its eigenvalues do not move by more than $\| E \| _ { F }$ .

Theorem 8.1.4 (Wielandt-Hoffman). If A and $A + E$ are $n { - } b y { - } n$ symmetric matrices, then

$$
\sum_ {i = 1} ^ {n} \left(\lambda_ {i} (A + E) - \lambda_ {i} (A)\right) ^ {2} \leq \| E \| _ {F} ^ {2}.
$$

Proof. See Wilkinson (AEP, pp. 104–108), Stewart and Sun (MPT, pp. 189–191), or Lax (1997, pp. 134–136).

Theorem 8.1.5. If A and $A + E$ are $n { - } b y { - } n$ symmetric matrices, then

$$
\lambda_ {k} (A) + \lambda_ {n} (E) \leq \lambda_ {k} (A + E) \leq \lambda_ {k} (A) + \lambda_ {1} (E), \quad k = 1: n.
$$

Proof. This follows from the minimax characterization. For details see Wilkinson (AEP, pp. 101–102) or Stewart and Sun (MPT, p. 203).

Corollary 8.1.6. If A and $A + E$ are $n { - } b y { - } n$ symmetric matrices, then

$$
\left| \lambda_ {k} (A + E) - \lambda_ {k} (A) \right| \leq \| E \| _ {2}
$$

for k = 1:n.

Proof. Observe that

$$
\left| \lambda_ {k} (A + E) - \lambda_ {k} (A) \right| \leq \max \left\{\left| \lambda_ {n} (E) \right|, \left| \lambda_ {1} (E) \right\| \right\} = \| E \| _ {2}
$$

for $k = 1 { : } n$ .

A pair of additional perturbation results that are important follow from the minimax property.

Theorem 8.1.7 (Interlacing Property). If $A \in \mathbb { R } ^ { n \times n }$ is symmetric and $A _ { r } \ =$ $A ( 1 { : } r , 1 { : } r )$ , then

$$
\lambda_ {r + 1} (A _ {r + 1}) \leq \lambda_ {r} (A _ {r}) \leq \lambda_ {r} (A _ {r + 1}) \leq \dots \leq \lambda_ {2} (A _ {r + 1}) \leq \lambda_ {1} (A _ {r}) \leq \lambda_ {1} (A _ {r + 1})
$$

for $r = 1 { : } n - 1$ .

Proof. Wilkinson (AEP, pp. 103–104).

Theorem 8.1.8. Suppose $\boldsymbol { B } = \boldsymbol { A } + \tau c c ^ { T }$ where $A \in \mathbb { R } ^ { n \times n }$ is symmetric, $c \in \mathbb { R } ^ { n }$ has unit 2-norm, and $\tau \in \mathbb { R }$ . $I f \tau \geq 0$ , then

$$
\lambda_ {i} (B) \in [ \lambda_ {i} (A), \lambda_ {i - 1} (A) ], \quad i = 2: n,
$$

while if $\tau \leq 0$ then

$$
\lambda_ {i} (B) \in [ \lambda_ {i + 1} (A), \lambda_ {i} (A) ], \quad i = 1: n - 1.
$$

In either case, there exist nonnegative $m _ { 1 } , \ldots , m _ { n }$ such that

$$
\lambda_ {i} (B) = \lambda_ {i} (A) + m _ {i} \tau , \qquad i = 1: n
$$

with $m _ { 1 } + \cdots + m _ { n } = 1$ .

Proof. Wilkinson (AEP, pp. 94–97). See also P8.1.8.

# 8.1.3 Invariant Subspaces

If $S \subseteq \mathbb { R } ^ { n }$ and $x \in S \Rightarrow A x \in S$ , then S is an invariant subspace for $A \in \mathbb { R } ^ { n \times n }$ . Note that if $\boldsymbol { x } \in \mathbb { R } ^ { i } ;$ s an eigenvector for A, then $S = \mathsf { s p a n } \{ x \}$ is 1-dimensional invariant subspace. Invariant subspaces serve to “take apart” the eigenvalue problem and figure heavily in many solution frameworks. The following theorem explains why.

Theorem 8.1.9. Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and that

$$
Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ r & n - r \end{array} \right]
$$

is orthogonal. If $\mathsf { r a n } ( Q _ { 1 } )$ is an invariant subspace, then

$$
Q ^ {T} A Q = D = \left[ \begin{array}{c c} D _ {1} & 0 \\ 0 & D _ {2} \end{array} \right] _ {n - r} ^ {r} \tag {8.1.1}
$$

and $\lambda ( A ) = \lambda ( D _ { 1 } ) \cup \lambda ( D _ { 2 } )$ . Compare with Lemma 7.1.2.

Proof. If

$$
Q ^ {T} A Q = \left[ \begin{array}{c c} D _ {1} & E _ {2 1} ^ {T} \\ E _ {2 1} & D _ {2} \end{array} \right],
$$

then from $A Q = Q D$ we have $A Q _ { 1 } - Q _ { 1 } D _ { 1 } \ = \ Q _ { 2 } E _ { 2 1 }$ . Since $\mathsf { r a n } ( Q _ { 1 } )$ is invariant, the columns of $Q _ { 2 } E _ { 2 1 }$ are also in $\mathsf { r a n } ( Q _ { 1 } )$ and therefore perpendicular to the columns of $Q _ { 2 }$ . Thus,

$$
0 = Q _ {2} ^ {T} \left(A Q _ {1} - Q _ {1} D _ {1}\right) = Q _ {2} ^ {T} Q _ {2} E _ {2 1} = E _ {2 1}.
$$

and so (8.1.1) holds. It is easy to show

$$
\det (A - \lambda I _ {n}) = \det (Q ^ {T} A Q - \lambda I _ {n}) = \det (D _ {1} - \lambda I _ {r}) \cdot \det (D _ {2} - \lambda I _ {n - r})
$$

confirming that $\lambda ( A ) = \lambda ( D _ { 1 } ) \cup \lambda ( D _ { 2 } )$ .

The sensitivity to perturbation of an invariant subspace depends upon the separation of the associated eigenvalues from the rest of the spectrum. The appropriate measure of separation between the eigenvalues of two symmetric matrices B and C is given by

$$
\operatorname{sep}(B,C) = \min_{\substack{\lambda \in \lambda (B)\\ \mu \in \lambda (C)}}|\lambda -\mu |. \tag{8.1.2}
$$

With this definition we have the following result.

Theorem 8.1.10. Suppose A and $A + E$ are $n { - } b y { - } n$ symmetric matrices and that

$$
Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ r & n - r \end{array} \right]
$$

is an orthogonal matrix such that ran $\left( Q _ { 1 } \right)$ is an invariant subspace for A. Partition the matrices $Q ^ { T } A Q$ and $Q ^ { T } E Q$ as follows:

$$
Q ^ {T} A Q = \left[ \begin{array}{c c} D _ {1} & 0 \\ 0 & D _ {2} \end{array} \right] _ {n - r} ^ {r}, \quad Q ^ {T} E Q = \left[ \begin{array}{c c} E _ {1 1} & E _ {2 1} ^ {T} \\ E _ {2 1} & E _ {2 2} \end{array} \right] _ {n - r} ^ {r}.
$$

If $\mathsf { s e p } ( D _ { 1 } , D _ { 2 } ) > 0$ and

$$
\| E \| _ {F} \leq \frac {\operatorname{sep} \left(D _ {1} , D _ {2}\right)}{5},
$$

then there exists a matrix $P \in \mathbb { R } ^ { ( n - r ) \times r }$ with

$$
\| P \| _ {F} \leq \frac {4}{\operatorname{sep} \left(D _ {1} , D _ {2}\right)} \| E _ {2 1} \| _ {F}
$$

such that the columns of $\hat { Q } _ { 1 } = ( Q _ { 1 } + Q _ { 2 } P ) ( I + P ^ { T } P ) ^ { - 1 / 2 }$ define an orthonormal basis for a subspace that is invariant $f o r \ A + E$ . Compare with Theorem $\it 7 . 2 . 4$ .

Proof. This result is a slight adaptation of Theorem 4.11 in Stewart (1973). The matrix $( I + P ^ { T } P ) ^ { - 1 / 2 }$ is the inverse of the square root of $I + P ^ { T } P$ . See §4.2.4.

Corollary 8.1.11. If the conditions of the theorem hold, then

$$
\operatorname{dist} \left(\operatorname{ran} \left(Q _ {1}\right), \operatorname{ran} \left(\hat {Q} _ {1}\right)\right) \leq \frac {4}{\operatorname{sep} \left(D _ {1} , D _ {2}\right)} \| E _ {2 1} \| _ {F}.
$$

Compare with Corollary 7.2.5.

Proof. It can be shown using the SVD that

$$
\| P (I + P ^ {T} P) ^ {- 1 / 2} \| _ {2} \leq \| P \| _ {2} \leq \| P \| _ {F}. \tag {8.1.3}
$$

Since $Q _ { 2 } ^ { T } \hat { Q } _ { 1 } = P ( I + P ^ { T } P ) ^ { - 1 / 2 }$ it follows that

$$
\begin{array}{l} \operatorname{dist} \left(\operatorname{ran} \left(Q _ {1}\right), \operatorname{ran} \left(\hat {Q} _ {1}\right)\right) = \| Q _ {2} ^ {T} \hat {Q} _ {1} \| _ {2} = \| P (I + P ^ {H} P) ^ {- 1 / 2} \| _ {2} \\ \leq \| P \| _ {2} \leq 4 \| E _ {2 1} \| _ {F} / \mathsf {s e p} (D _ {1}, D _ {2}) \\ \end{array}
$$

completing the proof.

Thus, the reciprocal of $\mathsf { s e p } ( D _ { 1 } , D _ { 2 } )$ can be thought of as a condition number that measures the sensitivity of $\mathsf { r a n } ( Q _ { 1 } )$ as an invariant subspace.

The effect of perturbations on a single eigenvector is sufficiently important that we specialize the above results to this case.

Theorem 8.1.12. Suppose A and $A + E$ are $n { - } b y { - } n$ symmetric matrices and that

$$
Q = \left[ \begin{array}{c c} q _ {1} & Q _ {2} \\ 1 & n - 1 \end{array} \right]
$$

is an orthogonal matrix such that $q _ { 1 }$ is an eigenvector for A. Partition the matrices $Q ^ { T } A Q$ and $Q ^ { T } E Q$ as follows:

$$
Q ^ {T} A Q   =   \left[ \begin{array}{c c} \lambda & 0 \\ 0 & D _ {2} \end{array} \right] _ {n - 1} ^ {1}  , \qquad Q ^ {T} E Q   =   \left[ \begin{array}{c c} \epsilon & e ^ {T} \\ e & E _ {2 2} \end{array} \right] _ {n - 1} ^ {1}  .
$$

If

$$
d = \min _ {\mu \in \lambda (D _ {2})} | \lambda - \mu | > 0
$$

and

$$
\parallel E \parallel_ {F} \leq \frac {d}{5},
$$

then there exists $p \in \mathbb { R } ^ { n - 1 }$ satisfying

$$
\| p \| _ {2} \leq \frac {4}{d} \| e \| _ {2}
$$

such that $\hat { q } _ { 1 } = ( q _ { 1 } + Q _ { 2 } p ) / \sqrt { 1 + p ^ { T } p }$ is a unit 2-norm eigenvector for $A + E$ . Moreover,

$$
\operatorname{dist} \left(\operatorname{span} \left\{q _ {1} \right\}, \operatorname{span} \left\{\hat {q} _ {1} \right\}\right) = \sqrt {1 - \left(q _ {1} ^ {T} \hat {q} _ {1}\right) ^ {2}} \leq \frac {4}{d} \| e \| _ {2}.
$$

Compare with Corollary 7.2.6.

Proof. Apply Theorem 8.1.10 and Corollary 8.1.11 with $r = 1$ and observe that if $D _ { 1 } = \left( \lambda \right)$ , then $d = { \tt s e p } ( D _ { 1 } , D _ { 2 } )$ .

# 8.1.4 Approximate Invariant Subspaces

If the columns of $Q _ { 1 } \in \mathbb { R } ^ { n \times r }$ are independent and the residual matrix $R = A Q _ { 1 } - Q _ { 1 } S$ is small for some $S \in \mathbb { R } ^ { r \times r }$ , then the columns of $Q _ { 1 }$ define an approximate invariant subspace. Let us discover what we can say about the eigensystem of A when in the possession of such a matrix.

Theorem 8.1.13. Suppose $A \in \mathbb { R } ^ { n \times n }$ and $S \in \mathbb { R } ^ { r \times r }$ are symmetric and that

$$
A Q _ {1} - Q _ {1} S = E _ {1}
$$

where $Q _ { 1 } \in \mathbb { R } ^ { n \times r }$ satisfies $Q _ { 1 } ^ { T } Q _ { 1 } = I _ { r }$ . Then there exist $\mu _ { 1 } , \ldots , \mu _ { r } \in \lambda ( A )$ such that

$$
| \mu_ {k} - \lambda_ {k} (S) | \leq \sqrt {2} \| E _ {1} \| _ {2}
$$

for $k = 1 { : } r$ .

Proof. Let $Q _ { 2 } \in \mathbb { R } ^ { n \times ( n - r ) }$ be any matrix such that $Q = \left[ Q _ { 1 } \mid Q _ { 2 } \right]$ is orthogonal. It follows that

$$
Q ^ {T} A Q = \left[ \begin{array}{c c} S & 0 \\ 0 & Q _ {2} ^ {T} A Q _ {2} \end{array} \right] + \left[ \begin{array}{c c} Q _ {1} ^ {T} E _ {1} & E _ {1} ^ {T} Q _ {2} \\ Q _ {2} ^ {T} E _ {1} & 0 \end{array} \right] \equiv B + E
$$

and so by using Corollary 8.1.6 we have $| \lambda _ { k } ( A ) - \lambda _ { k } ( B ) | \ \leq \ \| \ E \| _ { 2 }$ for $k = 1 { : } n$ . Since $\lambda ( S ) \subseteq \lambda ( B )$ , there exist $\mu _ { 1 } , \ldots , \mu _ { r } \in \lambda ( A )$ such that $\begin{array} { l l l } { | \mu _ { k } - \lambda _ { k } ( S ) | } & { \leq } & { \parallel E \parallel _ { 2 } } \end{array}$ for $k = 1 { : } r$ . The theorem follows by noting that for any $\boldsymbol { x } \in \mathbb { R } ^ { r }$ and $y \in \mathbb { R } ^ { n - r }$ we have

$$
\left\| E \left[ \begin{array}{c} x \\ y \end{array} \right] \right\| _ {2} \leq \| E _ {1} x \| _ {2} + \| E _ {1} ^ {T} Q _ {2} y \| _ {2} \leq \| E _ {1} \| _ {2} \| x \| _ {2} + \| E _ {1} \| _ {2} \| y \| _ {2}
$$

from which we readily conclude that $\| E \| _ { 2 } \leq \sqrt { 2 } \| E _ { 1 } \| _ { 2 }$ .

The eigenvalue bounds in Theorem 8.1.13 depend on $\parallel A Q _ { 1 } - Q _ { 1 } S \parallel _ { 2 }$ . Given A and $Q _ { 1 }$ , the following theorem indicates how to choose S so that this quantity is minimized in the Frobenius norm.

Theorem 8.1.14. If $A \in \mathbb { R } ^ { n \times n }$ is symmetric and $Q _ { 1 } \in \mathbb { R } ^ { n \times r }$ has orthonormal columns, then

$$
\min _ {S \in \mathbb {R} ^ {r \times r}} \parallel A Q _ {1} - Q _ {1} S \parallel_ {F} = \parallel (I - Q _ {1} Q _ {1} ^ {T}) A Q _ {1} \parallel_ {F}
$$

and $S = Q _ { 1 } ^ { T } A Q _ { 1 }$ is the minimizer.

Proof. Let $Q _ { 2 } \in \mathbb { R } ^ { n \times ( n - r ) }$ be such that $Q \ = \ [ \ Q _ { 1 } , \ Q _ { 2 } \ ]$ is orthogonal. For any $S \in \mathbb { R } ^ { r \times r }$ we have

$$
\left\| A Q _ {1} - Q _ {1} S \right\| _ {F} ^ {2} = \left\| Q ^ {T} A Q _ {1} - Q ^ {T} Q _ {1} S \right\| _ {F} ^ {2} = \left\| Q _ {1} ^ {T} A Q _ {1} - S \right\| _ {F} ^ {2} + \left\| Q _ {2} ^ {T} A Q _ {1} \right\| _ {F} ^ {2}.
$$

Clearly, the minimizing S is given by $S = Q _ { 1 } ^ { T } A Q _ { 1 }$ .

This result enables us to associate any r-dimensional subspace ran $\left( Q _ { 1 } \right)$ , with a set of $r$ “optimal” eigenvalue-eigenvector approximates.

Theorem 8.1.15. Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and that $Q _ { 1 } \in \mathbb { R } ^ { n \times r }$ satisfies $Q _ { 1 } ^ { T } Q _ { 1 } = I _ { r }$ . If

$$
Z ^ {T} (Q _ {1} ^ {T} A Q _ {1}) Z = \operatorname{diag} (\theta_ {1}, \dots , \theta_ {r}) = D
$$

is the Schur decomposition of $Q _ { 1 } ^ { T } A Q _ { 1 }$ and $Q _ { 1 } Z = \left[ y _ { 1 } \ : | \cdots | \ : y _ { r } \ : \right]$ , then

$$
\left\| A y _ {k} - \theta_ {k} y _ {k} \right\| _ {2} = \left\| \left(I - Q _ {1} Q _ {1} ^ {T}\right) A Q _ {1} Z e _ {k} \right\| _ {2} \leq \left\| \left(I - Q _ {1} Q _ {1} ^ {T}\right) A Q _ {1} \right\| _ {2}
$$

for $k = 1 { : } r$ .

Proof. It is easy to show that

$$
A y _ {k} - \theta_ {k} y _ {k} = A Q _ {1} Z e _ {k} - Q _ {1} Z D e _ {k} = (A Q _ {1} - Q _ {1} (Q _ {1} ^ {T} A Q _ {1})) Z e _ {k}.
$$

The theorem follows by taking norms.

In Theorem 8.1.15, the $\theta _ { k }$ are called Ritz values, the $y _ { k }$ are called Ritz vectors, and the $( \theta _ { k } , y _ { k } )$ are called Ritz pairs.

The usefulness of Theorem 8.1.13 is enhanced if we weaken the assumption that the columns of $Q _ { 1 }$ are orthonormal. As can be expected, the bounds deteriorate with the loss of orthogonality.

Theorem 8.1.16. Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and that

$$
A X _ {1} - X _ {1} S = F _ {1},
$$

where $X _ { 1 } \in \mathbb { R } ^ { n \times r }$ and $S = X _ { 1 } ^ { T } A X _ { 1 }$ . If

$$
\| X _ {1} ^ {T} X _ {1} - I _ {r} \| _ {2} = \tau <   1, \tag {8.1.4}
$$

then there exist $\mu _ { 1 } , \ldots , \mu _ { r } \in \lambda ( A )$ such that

$$
\left| \mu_ {k} - \lambda_ {k} (S) \right| \leq \sqrt {2} \left(\left\| F _ {1} \right\| _ {2} + \tau (2 + \tau) \right\| A \| _ {2})
$$

for $k = 1 { : } r$

Proof. For any $Q \in \mathbb { R } ^ { n \times r }$ with orthonormal columns, define $E _ { 1 } \in \mathbb { R } ^ { n \times r }$ by

$$
E _ {1} = A Q - Q S.
$$

It follows that

$$
E _ {1} = A (Q - X _ {1}) - (Q - X _ {1}) S + F _ {1}
$$

and so

$$
\left\| E _ {1} \right\| _ {2} \leq \left\| F _ {1} \right\| _ {2} + \left\| Q - X \right\| _ {2} \left\| A \right\| _ {2} \left(1 + \left\| X _ {1} \right\| _ {2} ^ {2}\right). \tag {8.1.5}
$$

Note that

$$
\left\| X _ {1} \right\| _ {2} ^ {2} = \left\| X _ {1} ^ {T} X _ {1} \right\| _ {2} \leq \left\| X ^ {T} X _ {1} - I _ {r} \right\| _ {2} + \left\| I _ {r} \right\| _ {2} = 1 + \tau . \tag {8.1.6}
$$

Let $U ^ { T } X _ { 1 } V ~ = ~ \Sigma ~ = ~ \mathrm { d i a g } ( \sigma _ { 1 } , . . . , \sigma _ { r } )$ be the thin SVD of $X _ { 1 }$ . It follows from (8.1.4) that

$$
\| \Sigma^ {2} - I _ {r} \| _ {2} = \tau
$$

and thus $1 - \sigma _ { r } ^ { 2 } = \tau$ . This implies

$$
\parallel Q - X _ {1} \parallel_ {2} = \parallel U (I _ {r} - \Sigma) V ^ {T} \parallel_ {2} = \parallel I _ {r} - \Sigma \parallel_ {2} = 1 - \sigma_ {r} \leq 1 - \sigma_ {r} ^ {2} = \tau . \tag {8.1.7}
$$

The theorem is established by substituting (8.1.6) and (8.1.7) into (8.1.5) and using Theorem 8.1.13.

# 8.1.5 The Law of Inertia

The inertia of a symmetric matrix A is a triplet of nonnegative integers $( m , z , p )$ where m, z, and p are respectively the numbers of negative, zero, and positive eigenvalues.

Theorem 8.1.17 (Sylvester Law of Inertia). If $A \in \mathbb { R } ^ { n \times n }$ is symmetric and $X \in \mathbb { R } ^ { n \times n }$ is nonsingular, then A and $X ^ { T } A X$ have the same inertia.

Proof. Suppose for some r that $\lambda _ { r } ( A ) > 0$ and define the subspace $S _ { 0 } \subseteq \mathbb { R } ^ { n }$ by

$$
S _ {0} = \operatorname{span} \left\{X ^ {- 1} q _ {1}, \dots , X ^ {- 1} q _ {r} \right\}, \quad q _ {i} \neq 0,
$$

where $A q _ { i } = \lambda _ { i } ( A ) q _ { i }$ and $i = 1 { : } r$ . From the minimax characterization of $\lambda _ { r } ( X ^ { T } A X )$ we have

$$
\lambda_ {r} (X ^ {T} A X) = \max _ {\dim (S) = r} \min _ {y \in S} \frac {y ^ {T} (X ^ {T} A X) y}{y ^ {T} y} \geq \min _ {y \in S _ {0}} \frac {y ^ {T} (X ^ {T} A X) y}{y ^ {T} y}.
$$

$$
y \in \mathbb {R} ^ {n} \Rightarrow \frac {y ^ {T} (X ^ {T} X) y}{y ^ {T} y} \geq \sigma_ {n} (X) ^ {2} \qquad y \in S _ {0} \Rightarrow \frac {y ^ {T} (X ^ {T} A X) y}{y ^ {T} (X ^ {T} X) y} \geq \lambda_ {r} (A),
$$

it follows that

$$
\lambda_ {r} (X ^ {T} A X) \geq \min _ {y \in S _ {0}} \left\{\frac {y ^ {T} (X ^ {T} A X) y}{y ^ {T} (X ^ {T} X) y} \frac {y ^ {T} (X ^ {T} X) y}{y ^ {T} y} \right\} \geq \lambda_ {r} (A) \sigma_ {n} (X) ^ {2}.
$$

An analogous argument with the roles of A and $X ^ { T } A X$ reversed shows that

$$
\lambda_ {r} (A) \geq \lambda_ {r} (X ^ {T} A X) \sigma_ {n} (X ^ {- 1}) ^ {2} = \frac {\lambda_ {r} (X ^ {T} A X)}{\sigma_ {1} (X) ^ {2}}.
$$

Thus, $\lambda _ { r } ( A )$ and $\lambda _ { r } ( X ^ { T } A X )$ have the same sign and so we have shown that A and $X ^ { T } A X$ have the same number of positive eigenvalues. If we apply this result to −A, we conclude that A and $X ^ { T } A X$ have the same number of negative eigenvalues. Obviously, the number of zero eigenvalues possessed by each matrix is also the same.

A transformation of the form $A \to X ^ { T } A X$ where X is nonsingular is called a conguence transformation. Thus, a congruence transformation of a symmetric matrix preserves inertia.

# Problems

P8.1.1 Without using any of the results in this section, show that the eigenvalues of a 2-by-2 symmetric matrix must be real.

P8.1.2 Compute the Schur decomposition of $A = { \left[ \begin{array} { l l } { 1 } & { 2 } \\ { 2 } & { 3 } \end{array} \right] }$

P8.1.3 Show that the eigenvalues of a Hermitian matrix $( A ^ { H } = A )$ are real. For each theorem and corollary in this section, state and prove the corresponding result for Hermitian matrices. Which results have analogs when A is skew-symmetric? Hint: If $A ^ { \overset { \triangledown } { T } } = - A$ , then iA is Hermitian.

P8.1.4 Show that if $X \in \mathbb { R } ^ { n \times r } , r \leq n ,$ and $\| X ^ { T } X - I \| _ { 2 } = \tau < 1$ , then $\sigma _ { \operatorname* { m i n } } ( X ) \geq 1 - \tau .$ .

P8.1.5 Suppose $\mathbf { l } , E \in \mathbb { R } ^ { n \times n }$ are symmetric and consider the Schur decomposition $A + t E = Q D Q ^ { T }$ where we assume that $Q = Q ( t )$ and $D = D ( t )$ are continuously differentiable functions of $t \in \mathbb { R }$ . Show that $\dot { D } ( t ) ~ = ~ \mathrm { d i a g } ( Q ( t ) ^ { T } E Q ( t ) )$ where the matrix on the right is the diagonal part of $Q ( t ) ^ { T } E Q ( t )$ . Establish the Wielandt-Hoffman theorem by integrating both sides of this equation from 0 to 1 and taking Frobenius norms to show that

$$
\| D (1) - D (0) \| _ {F} \leq \int_ {0} ^ {1} \| \operatorname{diag} (Q (t) ^ {T} E Q (t) \| _ {F} d t \leq \| E \| _ {F}.
$$

P8.1.6 Prove Theorem 8.1.5.

P8.1.7 Prove Theorem 8.1.7.

P8.1.8 Prove Theorem 8.1.8 using the fact that the trace of a square matrix is the sum of its eigenvalues.

P8.1.9 Show that if $B \in \mathbb { R } ^ { m \times m }$ and $C \in \mathbb { R } ^ { n \times n }$ are symmetric, then $\operatorname { s e p } ( B , C ) = \operatorname* { m i n } \parallel B X - X C \parallel _ { F }$ where the min is taken over all matrices $X \in \mathbb { R } ^ { m \times n }$ .

P8.1.10 Prove the inequality (8.1.3).

P8.1.11 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and $C \in \mathbb { R } ^ { n \times r }$ has full column rank and assume that $r \ll n$ . By using Theorem 8.1.8 relate the eigenvalues of $A + C C ^ { T }$ to the eigenvalues of A.

P8.1.12 Give an algorithm for computing the solution to

$$
\min \quad \| A - S \| _ {F}.
$$

$$
\operatorname{rank} (S) = 1
$$

$$
S = S ^ {T}
$$

Note that if $S \in \mathbb { R } ^ { n \times n }$ is a symmetric rank-1 matrix then either ${ \boldsymbol { S } } = v { \boldsymbol { v } } ^ { T } ~ { \mathrm { o r } } ~ { \boldsymbol { S } } = - v { \boldsymbol { v } } ^ { T }$ for some $v \in \mathbb { R } ^ { n }$ .

P8.1.13 Give an algorithm for computing the solution to

$$
\min \quad \| A - S \| _ {F}.
$$

$$
\operatorname{rank} (S) = 2
$$

$$
S = - S ^ {T}
$$

P8.1.14 Give an example of a real 3-by-3 normal matrix with integer entries that is neither orthogonal, symmetric, nor skew-symmetric.

# Notes and References for §8.1

The perturbation theory for the symmetric eigenproblem is surveyed in Wilkinson (AEP, Chap. 2), Parlett (SEP, Chaps. 10 and 11), and Stewart and Sun (MPT, Chaps. 4 and 5). Some representative papers in this well-researched area include:

G.W. Stewart (1973). “Error and Perturbation Bounds for Subspaces Associated with Certain Eigenvalue Problems,” SIAM Review 15, 727–764.

C.C. Paige (1974). “Eigenvalues of Perturbed Hermitian Matrices,” Lin. Alg. Applic. 8, 1–10.

W. Kahan (1975). “Spectra of Nearly Hermitian Matrices,” Proc. AMS 48, 11–17.

A. Schonhage (1979). “Arbitrary Perturbations of Hermitian Matrices,” Lin. Alg. Applic. 24, 143–49.

D.S. Scott (1985). “On the Accuracy of the Gershgorin Circle Theorem for Bounding the Spread of a Real Symmetric Matrix,” Lin. Alg. Applic. 65, 147–155

J.-G. Sun (1995). “A Note on Backward Error Perturbations for the Hermitian Eigenvalue Problem,” BIT 35, 385–393.

Z. Drmaˇc (1996). On Relative Residual Bounds for the Eigenvalues of a Hermitian Matrix,” Lin. Alg. Applic. 244, 155-163.

Z. Drmaˇc and V. Hari (1997). “Relative Residual Bounds For The Eigenvalues of a Hermitian Semidefinite Matrix,” SIAM J. Matrix Anal. Applic. 18, 21–29.

R.-C. Li (1998). “Relative Perturbation Theory: I. Eigenvalue and Singular Value Variations,” SIAM J. Matrix Anal. Applic. 19, 956–982.

R.-C. Li (1998). “Relative Perturbation Theory: II. Eigenspace and Singular Subspace Variations,” SIAM J. Matrix Anal. Applic. 20, 471–492.

F.M. Dopico, J. Moro and J.M. Molera (2000). “Weyl-Type Relative Perturbation Bounds for Eigensystems of Hermitian Matrices,” Lin. Alg. Applic. 309, 3–18.

J.L. Barlow and I. Slapniˇcar (2000). “Optimal Perturbation Bounds for the Hermitian Eigenvalue Problem,” Lin. Alg. Applic. 309, 19–43.

N. Truhar and R.-C. Li (2003). “A sin(2θ) Theorem for Graded Indefinite Hermitian Matrices,” Lin. Alg. Applic. 359, 263–276.

W. Li and W. Sun (2004). “The Perturbation Bounds for Eigenvalues of Normal Matrices,” Num. Lin. Alg. 12, 89–94.

C.-K. Li and R.-C. Li (2005). “A Note on Eigenvalues of Perturbed Hermitian Matrices,” Lin. Alg. Applic. 395, 183–190.

N. Truhar (2006). “Relative Residual Bounds for Eigenvalues of Hermitian Matrices,” SIAM J. Matrix Anal. Applic. 28, 949–960.

An elementary proof of the Wielandt-Hoffman theorem is given in:

P. Lax (1997). Linear Algebra, Wiley-Interscience, New York.

For connections to optimization and differential equations, see:

P. Deift, T. Nanda, and C. Tomei (1983). “Ordinary Differential Equations and the Symmetric Eigenvalue Problem,” SIAM J. Numer. Anal. 20, 1–22.

M.L. Overton (1988). “Minimizing the Maximum Eigenvalue of a Symmetric Matrix,” SIAM J. Matrix Anal. Applic. 9, 256-268.

T. Kollo and H. Neudecker (1997). “The Derivative of an Orthogonal Matrix of Eigenvectors of a Symmetric Matrix,” Lin. Alg. Applic. 264, 489–493.
