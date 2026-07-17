# 8.2 Power Iterations

Assume that $A \in \mathbb { R } ^ { n \times n }$ is symmetric and that $U _ { 0 } \in \mathbb { R } ^ { n \times n }$ is orthogonal. Consider the following QR iteration:

$$
T _ {0} = U _ {0} ^ {T} A U _ {0}
$$

for k = 1, 2, . . .

$$
T _ {k - 1} = U _ {k} R _ {k} \quad (\text { QR   factorization }) \tag {8.2.1}
$$

$$
T _ {k} = R _ {k} U _ {k}
$$

end

Since $T _ { k } = R _ { k } U _ { k } = U _ { k } ^ { T } ( U _ { k } R _ { k } ) U _ { k } = U _ { k } ^ { T } T _ { k - 1 } U _ { k }$ it follows by induction that

$$
T _ {k} = \left(U _ {0} U _ {1} \dots U _ {k}\right) ^ {T} A \left(U _ {0} U _ {1} \dots U _ {k}\right). \tag {8.2.2}
$$

Thus, each $T _ { k }$ is orthogonally similar to A. Moreover, the $T _ { k }$ almost always converge to diagonal form and so it can be said that (8.2.1) almost always converges to a Schur decomposition of A. In order to establish this remarkable result we first consider the power method and the method of orthogonal iteration.

# 8.2.1 The Power Method

Given a unit 2-norm $q ^ { ( 0 ) } \in \mathbb { R } ^ { n }$ , the power method produces a sequence of vectors $q ^ { ( k ) }$ as follows:

for k = 1, 2, . . .

$$
z ^ {(k)} = A q ^ {(k - 1)}
$$

$$
q ^ {(k)} = z ^ {(k)} / \| z ^ {(k)} \| _ {2} \tag {8.2.3}
$$

$$
\lambda^ {(k)} = \left[ q ^ {(k)} \right] ^ {T} A q ^ {(k)}
$$

end

If $q ^ { ( 0 ) }$ is not “deficient” and A’s eigenvalue of maximum modulus is unique, then the q(k) $q ^ { ( k ) }$ converge to an eigenvector.

Theorem 8.2.1. Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and that

$$
Q ^ {T} A Q = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right)
$$

where $Q = [ q _ { 1 } | \cdots | q _ { n } ]$ is orthogonal and $\left| { \lambda } _ { 1 } \right| > \left| { \lambda } _ { 2 } \right| \geq \cdot \cdot \cdot \geq \left| { \lambda } _ { n } \right|$ . Let the vectors $q ^ { ( k ) }$ be specified by (8.2.3) and define $\theta _ { k } \in [ 0 , \pi / 2 ]$ by

$$
\cos (\theta_ {k}) = \left| q _ {1} ^ {T} q ^ {(k)} \right|.
$$

If $\cos ( \theta _ { 0 } ) \neq 0$ , then for $k = 0 , 1 , \ldots$ . we have

$$
\left| \sin (\theta_ {k}) \right| \leq \tan (\theta_ {0}) \left| \frac {\lambda_ {2}}{\lambda_ {1}} \right| ^ {k}, \tag {8.2.4}
$$

$$
\left| \lambda^ {(k)} - \lambda_ {1} \right| \leq \max _ {2 \leq i \leq n} \left| \lambda_ {1} - \lambda_ {i} \right| \tan (\theta_ {0}) ^ {2} \left| \frac {\lambda_ {2}}{\lambda_ {1}} \right| ^ {2 k}. \tag {8.2.5}
$$

Proof. From the definition of the iteration, it follows that $q ^ { ( k ) }$ is a multiple of $A ^ { k } q ^ { ( 0 ) }$ and so

$$
| \sin (\theta_ {k}) | ^ {2} = 1 - \left(q _ {1} ^ {T} q ^ {(k)}\right) ^ {2} = 1 - \left(\frac {q _ {1} ^ {T} A ^ {k} q ^ {(0)}}{\| A ^ {k} q ^ {(0)} \| _ {2}}\right) ^ {2}.
$$

If $q ^ { ( 0 ) }$ has the eigenvector expansion $q ^ { ( 0 ) } = a _ { 1 } q _ { 1 } + \cdot \cdot \cdot + a _ { n } q _ { n }$ , then

$$
\left| a _ {1} \right| = \left| q _ {1} ^ {T} q ^ {(0)} \right| = \cos \left(\theta_ {0}\right) \neq 0,
$$

$$
a _ {1} ^ {2} + \dots + a _ {n} ^ {2} = 1,
$$

and

$$
A ^ {k} q ^ {(0)} = a _ {1} \lambda_ {1} ^ {k} q _ {1} + a _ {2} \lambda_ {2} ^ {k} q _ {2} + \dots + a _ {n} \lambda_ {n} ^ {k} q _ {n}.
$$

Thus,

$$
\begin{array}{l} | \sin (\theta_ {k}) | ^ {2} = 1 - \frac {a _ {1} ^ {2} \lambda_ {1} ^ {2 k}}{\sum_ {i = 1} ^ {n} a _ {i} ^ {2} \lambda_ {i} ^ {2 k}} = \frac {\sum_ {i = 2} ^ {n} a _ {i} ^ {2} \lambda_ {i} ^ {2 k}}{\sum_ {i = 1} ^ {n} a _ {i} ^ {2} \lambda_ {i} ^ {2 k}} \leq \frac {\sum_ {i = 2} ^ {n} a _ {i} ^ {2} \lambda_ {i} ^ {2 k}}{a _ {1} ^ {2} \lambda_ {1} ^ {2 k}} \\ = \frac {1}{a _ {1} ^ {2}} \sum_ {i = 2} ^ {n} a _ {i} ^ {2} \left(\frac {\lambda_ {i}}{\lambda_ {1}}\right) ^ {2 k} \leq \frac {1}{a _ {1} ^ {2}} \left(\sum_ {i = 2} ^ {n} a _ {i} ^ {2}\right) \left(\frac {\lambda_ {2}}{\lambda_ {1}}\right) ^ {2 k} \\ = \frac {1 - a _ {1} ^ {2}}{a _ {1} ^ {2}} \left(\frac {\lambda_ {2}}{\lambda_ {1}}\right) ^ {2 k} = \tan (\theta_ {0}) ^ {2} \left(\frac {\lambda_ {2}}{\lambda_ {1}}\right) ^ {2 k}. \\ \end{array}
$$

This proves (8.2.4). Likewise,

$$
\lambda^ {(k)} = \left[ q ^ {(k)} \right] ^ {T} A q ^ {(k)} = \frac {\left[ q ^ {(0)} \right] ^ {T} A ^ {2 k + 1} q ^ {(0)}}{\left[ q ^ {(0)} \right] ^ {T} A ^ {2 k} q ^ {(0)}} = \frac {\sum_ {i = 1} ^ {n} a _ {i} ^ {2} \lambda_ {i} ^ {2 k + 1}}{\sum_ {i = 1} ^ {n} a _ {i} ^ {2} \lambda_ {i} ^ {2 k}}
$$

and so

$$
\begin{array}{l} \left| \lambda^ {(k)} - \lambda_ {1} \right| = \left| \frac {\sum_ {i = 2} ^ {n} a _ {i} ^ {2} \lambda_ {i} ^ {2 k} \left(\lambda_ {i} - \lambda_ {1}\right)}{\sum_ {i = 1} ^ {n} a _ {i} ^ {2} \lambda_ {i} ^ {2 k}} \right| \leq \max _ {2 \leq i \leq n} | \lambda_ {1} - \lambda_ {i} | \cdot \frac {1}{a _ {1} ^ {2}} \cdot \sum_ {i = 2} ^ {n} a _ {i} ^ {2} \left(\frac {\lambda_ {i}}{\lambda_ {1}}\right) ^ {2 k} \\ \leq \max _ {2 \leq i \leq n} | \lambda_ {1} - \lambda_ {n} | \cdot \tan (\theta_ {0}) ^ {2} \cdot \left(\frac {\lambda_ {2}}{\lambda_ {1}}\right) ^ {2 k}, \\ \end{array}
$$

completing the proof of the theorem.

Computable error bounds for the power method can be obtained by using Theorem 8.1.13. If

$$
\| A q ^ {(k)} - \lambda^ {(k)} q ^ {(k)} \| _ {2} = \delta ,
$$

then there exists $\lambda \in \lambda ( A )$ such that ${ | \lambda ^ { ( k ) } - \lambda | } \leq \sqrt { 2 } \delta$ .

# 8.2.2 Inverse Iteration

If the power method (8.2.3) is applied with A replaced by $( A - \lambda I ) ^ { - 1 }$ , then we obtain the method of inverse iteration. If λ is very close to a distinct eigenvalue of $A .$ , then $q ^ { ( k ) }$ will be much richer in the corresponding eigenvector direction than its predecessor $\bar { q } ^ { ( k - 1 ) }$ :

$$
\left. \begin{array}{l} x = \sum_ {i = 1} ^ {n} a _ {i} q _ {i} \\ A q _ {i} = \lambda_ {i} q _ {i}, i = 1: n \end{array} \right\} \Rightarrow (A - \lambda I) ^ {- 1} x = \sum_ {i = 1} ^ {n} \frac {a _ {i}}{\lambda_ {i} - \lambda} q _ {i}.
$$

Thus, if λ is reasonably close to a well-separated eigenvalue $\lambda _ { j }$ , then inverse iteration will produce iterates that are increasingly in the direction of $q _ { j }$ . Note that inverse iteration requires at each step the solution of a linear system with matrix of coefficients $A - \lambda I$ .

# 8.2.3 Rayleigh Quotient Iteration

Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and that x is a given nonzero n-vector. A simple differentiation reveals that

$$
\lambda = r (x) \equiv \frac {x ^ {T} A x}{x ^ {T} x}
$$

minimizes $\parallel ( A - \lambda I ) x \parallel _ { 2 }$ . (See also Theorem 8.1.14.) The scalar $r ( x )$ is called the Rayleigh quotient of x. Clearly, if x is an approximate eigenvector, then $r ( x )$ i s a reasonable choice for the corresponding eigenvalue. Combining this idea with inverse iteration gives rise to the Rayleigh quotient iteration where $x _ { 0 } \neq 0$ is given.

for $k = 0 , 1 , \ldots$

$$
\mu_ {k} = r (x _ {k}) \tag {8.2.6}
$$

${ \mathrm { S o l v e ~ } } ( A - \mu _ { k } I ) z _ { k + 1 } = x _ { k } { \mathrm { ~ f o r ~ } } z _ { k + 1 }$

$$
x _ {k + 1} = z _ {k + 1} / \left\| z _ {k + 1} \right\| _ {2}
$$

end

The Rayleigh quotient iteration almost always converges and when it does, the rate of convergence is cubic. We demonstrate this for the case $n = 2$ . Without loss of generality, we may assume that $A = \mathrm { d i a g } ( \lambda _ { 1 } , \lambda _ { 2 } )$ , with $\lambda _ { 1 } > \lambda _ { 2 }$ . Denoting $x _ { k }$ by

$$
x _ {k} = \left[ \begin{array}{l} c _ {k} \\ s _ {k} \end{array} \right], \qquad c _ {k} ^ {2} + s _ {k} ^ {2} = 1,
$$

it follows that $\mu _ { k } \ = \ \lambda _ { 1 } c _ { k } ^ { 2 } + \lambda _ { 2 } s _ { k } ^ { 2 }$ in (8.2.6) and

$$
z _ {k + 1} = \frac {1}{\lambda_ {1} - \lambda_ {2}} \left[ \begin{array}{c} c _ {k} / s _ {k} ^ {2} \\ - s _ {k} / c _ {k} ^ {2} \end{array} \right].
$$

A calculation shows that

$$
c _ {k + 1} = \frac {\left| c _ {k} \right| ^ {3}}{\sqrt {c _ {k} ^ {6} + s _ {k} ^ {6}}}, \quad s _ {k + 1} = \frac {\left| s _ {k} \right| ^ {3}}{\sqrt {c _ {k} ^ {6} + s _ {k} ^ {6}}}. \tag {8.2.7}
$$

From these equations it is clear that the $x _ { k }$ converge cubically to either span $\{ e _ { 1 } \}$ or span $\left\{ e _ { 2 } \right\}$ provided $\vert c _ { k } \vert \neq \vert s _ { k } \vert$ . Details associated with the practical implementation of the Rayleigh quotient iteration may be found in Parlett (1974).

# 8.2.4 Orthogonal Iteration

A straightforward generalization of the power method can be used to compute higherdimensional invariant subspaces. Let r be a chosen integer that satisfies $1 \leq r \leq$ n. Given an n-by-r matrix $Q _ { 0 }$ with orthonormal columns, the method of orthogonal iteration generates a sequence of matrices $\{ Q _ { k } \} \subseteq \mathbb { R } ^ { n \times r }$ as follows:

$\mathbf { f o r } \ k = 1 , 2 , \ldots$

$$
Z _ {k} = A Q _ {k - 1} \tag {8.2.8}
$$

$$
Q _ {k} R _ {k} = Z _ {k} \quad \text {(QR factorization)}
$$

end

Note that, if $r = 1$ , then this is just the power method. Moreover, the sequence $\{ Q _ { k } e _ { 1 } \}$ is precisely the sequence of vectors produced by the power iteration with starting vector $q ^ { ( \bar { 0 } ) } = Q _ { 0 } \dot { e } _ { 1 }$ .

In order to analyze the behavior of (8.2.8), assume that

$$
Q ^ {T} A Q = D = \operatorname{diag} (\lambda_ {i}), \quad | \lambda_ {1} | \geq | \lambda_ {2} | \geq \dots \geq | \lambda_ {n} | \tag {8.2.9}
$$

is a Schur decomposition of $A \in \mathbb { R } ^ { n \times n }$ . Partition $Q$ and D as follows:

$$
Q = \left[ \begin{array}{c c} Q _ {\alpha} & Q _ {\beta} \\ r & n - r \end{array} \right], \quad D = \left[ \begin{array}{c c} D _ {1} & 0 \\ 0 & D _ {2} \\ r & n - r \end{array} \right] _ {n - r} ^ {r}. \tag {8.2.10}
$$

If $| \lambda _ { r } | > | \lambda _ { r + 1 } |$ , then

$$
D _ {r} (A) = \operatorname{ran} (Q _ {\alpha})
$$

is the dominant invariant subspace of dimension r. It is the unique invariant subspace associated with the eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { r }$ .

The following theorem shows that with reasonable assumptions, the subspaces ran $\left( Q _ { k } \right)$ generated by (8.2.8) converge to $D _ { r } ( A )$ at a rate proportional to $| \lambda _ { r + 1 } / \lambda _ { r } | ^ { k }$ .

Theorem 8.2.2. Let the Schur decomposition of $A \in \mathbb { R } ^ { n \times n }$ be given by (8.2.9) and (8.2.10) with $n \geq 2$ . Assume $\left| \lambda _ { r } \right| > \left| \lambda _ { r + 1 } \right|$ and that $d _ { k }$ is defined by

$$
d _ {k} = \operatorname{dist} \left(D _ {r} (A), \operatorname{ran} \left(Q _ {k}\right)\right), \quad k \geq 0.
$$

If

$$
d _ {0} <   1, \tag {8.2.11}
$$

then the matrices $Q _ { k }$ generated by (8.2.8) satisfy

$$
d _ {k} \leq \left| \frac {\lambda_ {r + 1}}{\lambda_ {r}} \right| ^ {k} \frac {d _ {0}}{\sqrt {1 - d _ {0} ^ {2}}}. \tag {8.2.12}
$$

Compare with Theorem 7.3.1.

Proof. We mention at the start that the condition (8.2.11) means that no vector in the span of $Q _ { 0 } \mathrm { { ^ { 2 } s } }$ columns is perpendicular to $D _ { r } ( A )$ .

Using induction it can be shown that the matrix $Q _ { k }$ in (8.2.8) satisfies

$$
A ^ {k} Q _ {0} = Q _ {k} \left(R _ {k} \dots R _ {1}\right).
$$

This is a QR factorization of $A ^ { k } Q _ { 0 }$ and upon substitution of the Schur decomposition (8.2.9)-(8.2.10) we obtain

$$
\left[ \begin{array}{c c} D _ {1} ^ {k} & 0 \\ 0 & D _ {2} ^ {k} \end{array} \right] \left[ \begin{array}{c} Q _ {\alpha} ^ {T} Q _ {0} \\ Q _ {\beta} ^ {T} Q _ {0} \end{array} \right] = \left[ \begin{array}{c} Q _ {\alpha} ^ {T} Q _ {k} \\ Q _ {\beta} ^ {T} Q _ {k} \end{array} \right] (R _ {k} \dots R _ {1})  .
$$

If the matrices $V _ { k }$ and $W _ { k }$ are defined by

$$
V _ {k} = Q _ {\alpha} ^ {T} Q _ {0},
$$

$$
W _ {k} = Q _ {\beta} ^ {T} Q _ {0},
$$

then

$$
D _ {1} ^ {k} V _ {0} = V _ {k} \left(R _ {k} \dots R _ {1}\right), \tag {8.2.13}
$$

$$
D _ {2} ^ {k} W _ {0} = W _ {k} \left(R _ {k} \dots R _ {1}\right). \tag {8.2.14}
$$

Since

$$
\left[ \begin{array}{c} V _ {k} \\ W _ {k} \end{array} \right] = \left[ \begin{array}{c} Q _ {\alpha} ^ {T} Q _ {k} \\ Q _ {\beta} ^ {T} Q _ {k} \end{array} \right] = [ Q _ {\alpha} \mid Q _ {\beta} ] ^ {T} Q _ {k} = Q ^ {T} Q _ {k},
$$

it follows from the thin CS decomposition (Theorem 2.5.2) that

$$
1 = \sigma_ {\mathrm{min}} (V _ {k}) ^ {2} + \sigma_ {\mathrm{max}} (W _ {k}) ^ {2} = \sigma_ {\mathrm{min}} (V _ {k}) ^ {2} + d _ {k} ^ {2}.
$$

A consequence of this is that

$$
\sigma_ {\mathrm{min}} (V _ {0}) ^ {2} = 1 - \sigma_ {\mathrm{max}} (W _ {0}) ^ {2} = 1 - d _ {0} ^ {2} > 0.
$$

It follows from (8.2.13) that the matrices $V _ { k }$ and $( R _ { k } \cdot \cdot \cdot R _ { 1 } )$ are nonsingular. Using both that equation and (8.2.14) we obtain

$$
W _ {k} = D _ {2} ^ {k} W _ {0} (R _ {k} \dots R _ {1}) ^ {- 1} = D _ {2} ^ {k} W _ {0} (D _ {1} ^ {k} V _ {0}) ^ {- 1} V _ {k} = D _ {2} ^ {k} (W _ {0} V _ {0} ^ {- 1}) D _ {1} ^ {- k} V _ {k}
$$

and so

$$
\begin{array}{l} d _ {k} = \left\| W _ {k} \right\| _ {2} \leq \left\| D _ {2} ^ {k} \right\| _ {2} \cdot \left\| W _ {0} \right\| _ {2} \cdot \left\| V _ {0} ^ {- 1} \right\| _ {2} \cdot \left\| D _ {1} ^ {- k} \right\| _ {2} \cdot \left\| V _ {k} \right\| _ {2} \\ \leq | \lambda_ {r + 1} | ^ {k} \cdot d _ {0} \cdot \frac {1}{1 - d _ {0} ^ {2}} \cdot \frac {1}{| \lambda_ {r} | ^ {k}}, \\ \end{array}
$$

from which the theorem follows.

# 8.2.5 The QR Iteration

Consider what happens if we apply the method of orthogonal iteration (8.2.8) with $r = n$ . Let $Q ^ { T } A Q = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ be the Schur decomposition and assume

$$
\left| \lambda_ {1} \right| > \left| \lambda_ {2} \right| > \dots > \left| \lambda_ {n} \right|.
$$

${ \mathrm { I f ~ } } Q = \left[ q _ { 1 } \left| \cdots \right| q _ { n } \right] , Q _ { k } = \left[ q _ { 1 } ^ { ( k ) } \left| \cdots \right| q _ { n } ^ { ( k ) } \right] , { \mathrm { a n d } } $

$$
\operatorname{dist} (D _ {i} (A), \operatorname{span} \{q _ {1} ^ {(0)}, \dots , q _ {i} ^ {(0)} \}) <   1 \tag {8.2.15}
$$

for $i = 1 { : } n - 1$ , then it follows from Theorem 8.2.2 that

$$
\operatorname{dist} \left(\operatorname{span} \left\{q _ {1} ^ {(k)}, \dots , q _ {i} ^ {(k)} \right\}, \operatorname{span} \left\{q _ {1}, \dots , q _ {i} \right\}\right) = O \left(\left| \frac {\lambda_ {i + 1}}{\lambda_ {i}} \right| ^ {k}\right)
$$

for $i = 1 { : } n - 1$ . This implies that the matrices $T _ { k }$ defined by

$$
T _ {k} = Q _ {k} ^ {T} A Q _ {k}
$$

are converging to diagonal form. Thus, it can be said that the method of orthogonal iteration computes a Schur decomposition if $r = n$ and the original iterate $Q _ { 0 } \in \mathbb { R } ^ { n \times n }$ is not deficient in the sense of (8.2.11).

The QR iteration arises by considering how to compute the matrix $T _ { k }$ directly from its predecessor $T _ { k - 1 }$ . On the one hand, we have from (8.2.8) and the definition of $T _ { k - 1 }$ that

$$
T _ {k - 1} = Q _ {k - 1} ^ {T} A Q _ {k - 1} = Q _ {k - 1} ^ {T} (A Q _ {k - 1}) = (Q _ {k - 1} ^ {T} Q _ {k}) R _ {k}.
$$

On the other hand,

$$
T _ {k} = Q _ {k} ^ {T} A Q _ {k} = (Q _ {k} ^ {T} A Q _ {k - 1}) (Q _ {k - 1} ^ {T} Q _ {k}) = R _ {k} (Q _ {k - 1} ^ {T} Q _ {k}).
$$

Thus, $T _ { k }$ is determined by computing the QR factorization of $T _ { k - 1 }$ and then multiplying the factors together in reverse order. This is precisely what is done in (8.2.1).

Note that a single QR iteration involves $O ( n ^ { 3 } )$ flops. Moreover, since convergence is only linear (when it exists), it is clear that the method is a prohibitively expensive way to compute Schur decompositions. Fortunately, these practical difficulties can be overcome, as we show in the next section.

# Problems

P8.2.1 Suppose $A _ { 0 } \in \mathbb { R } ^ { n \times n }$ is symmetric and positive definite and consider the following iteration:

$$
\begin{array}{l} \text { for } k = 1, 2, \dots \\ A _ {k - 1} = G _ {k} G _ {k} ^ {T} \quad \text {(Cholesky factorization)} \\ A _ {k} = G _ {k} ^ {T} G _ {k} \\ \end{array}
$$

end

(a) Show that this iteration is defined. (b) Show that if

$$
A _ {0} = \left[ \begin{array}{c c} a & b \\ b & c \end{array} \right]
$$

with $a \geq c$ has eigenvalues $\lambda _ { 1 } \geq \lambda _ { 2 } > 0$ , then the $A _ { k }$ converge to $\mathrm { d i a g } ( \lambda _ { 1 } , \lambda _ { 2 } )$ .

P8.2.2 Prove (8.2.7).

P8.2.3 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and define the function $f { : \mathbb { R } ^ { n + 1 } } \to \mathbb { R } ^ { n + 1 }$ by

$$
f \left(\left[ \begin{array}{l} x \\ \lambda \end{array} \right]\right) = \left[ \begin{array}{c} A x - \lambda x \\ (x ^ {T} x - 1) / 2 \end{array} \right]
$$

where $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $\lambda \in \mathbb { R }$ . Suppose $x _ { + }$ and $\lambda _ { + }$ are produced by applying Newton’s method to f at the “current point” defined by $x _ { c }$ and $\lambda _ { c }$ . Give expressions for $x _ { + }$ and $\lambda _ { + }$ assuming that $\| \ b { x } _ { c } \| _ { 2 } = 1$ and $\lambda _ { c } = x _ { c } ^ { T } A x _ { c }$ .

# Notes and References for 8.2

The following references are concerned with the method of orthogonal iteration, which is also known as the method of simultaneous iteration:

G.W. Stewart (1969). “Accelerating The Orthogonal Iteration for the Eigenvalues of a Hermitian Matrix,” Numer. Math. 13, 362–376.   
M. Clint and A. Jennings (1970). “The Evaluation of Eigenvalues and Eigenvectors of Real Symmetric Matrices by Simultaneous Iteration,” Comput. J. 13, 76–80.   
H. Rutishauser (1970). “Simultaneous Iteration Method for Symmetric Matrices,” Numer. Math. 16, 205–223.

References for the Rayleigh quotient method include:

J. Vandergraft (1971). “Generalized Rayleigh Methods with Applications to Finding Eigenvalues of Large Matrices,” Lin. Alg. Applic. 4, 353–368.   
B.N. Parlett (1974). “The Rayleigh Quotient Iteration and Some Generalizations for Nonnormal Matrices,” Math. Comput. 28, 679-693.   
S. Batterson and J. Smillie (1989). “The Dynamics of Rayleigh Quotient Iteration,” SIAM J. Numer. Anal. 26, 624–636.   
C. Beattie and D.W. Fox (1989). “Localization Criteria and Containment for Rayleigh Quotient Iteration,” SIAM J. Matrix Anal. Applic. 10, 80–93.   
P.T.P. Tang (1994). “Dynamic Condition Estimation and Rayleigh-Ritz Approximation,” SIAM J. Matrix Anal. Applic. 15, 331–346.   
D. P. O’Leary and G. W. Stewart (1998). “On the Convergence of a New Rayleigh Quotient Method with Applications to Large Eigenproblems,” ETNA 7, 182–189.   
J.-L. Fattebert (1998). “A Block Rayleigh Quotient Iteration with Local Quadratic Convergence,” ETNA 7, 56–74.   
Z. Jia and G.W. Stewart (2001). “An Analysis of the Rayleigh-Ritz Method for Approximating Eigenspaces,” Math. Comput. 70, 637–647.   
V. Simoncini and L. Eld´en (2002). “Inexact Rayleigh Quotient-Type Methods for Eigenvalue Computations,” BIT 42, 159–182.   
P.A. Absil, R. Mahony, R. Sepulchre, and P. Van Dooren (2002). “A Grassmann-Rayleigh Quotient Iteration for Computing Invariant Subspaces,” SIAM Review 44, 57–73.   
Y. Notay (2003). “Convergence Analysis of Inexact Rayleigh Quotient Iteration,” SIAM J. Matrix Anal. Applic. 24, 627–644.   
A. Dax (2003). “The Orthogonal Rayleigh Quotient Iteration (ORQI) method,” Lin. Alg. Applic. 358, 23–43.   
R.-C. Li (2004). “Accuracy of Computed Eigenvectors Via Optimizing a Rayleigh Quotient,” BIT 44, 585–593.

Various Newton-type methods have also been derived for the symmetric eigenvalue problem, see:

R.A. Tapia and D.L. Whitley (1988). “The Projected Newton Method Has Order $1 + { \sqrt { 2 } }$ for the Symmetric Eigenvalue Problem,” SIAM J. Numer. Anal. 25, 1376–1382.

P.A. Absil, R. Sepulchre, P. Van Dooren, and R. Mahony (2004). “Cubically Convergent Iterations for Invariant Subspace Computation,” SIAM J. Matrix Anal. Applic. 26, 70–96.
