# 10.1 The Symmetric Lanczos Process

Suppose $A \in \mathbb { R } ^ { n \times n }$ is large, sparse, and symmetric and assume that a few of its largest and/or smallest eigenvalues are desired. Eigenvalues at either end of the spectrum are referred to as extremal eigenvalues. This problem can be addressed by a method attributed to Lanczos (1950). The method generates a sequence of tridiagonal matrices $\left\{ T _ { k } \right\}$ with the property that the extremal eigenvalues of $T _ { k } \in \mathbb { R } ^ { k \times k }$ are progressively better estimates of A’s extremal eigenvalues. In this section, we derive the technique and investigate some of its exact arithmetic properties.

One way to motivate the Lanczos idea is to be reminded about the shortcomings of the power method that we discussed in §8.2.1. Recall that the power method can be used to find the dominant eigenvalue $\lambda _ { 1 }$ and an associated eigenvector $x _ { 1 }$ . However, the rate of convergence is dictated by $| \lambda _ { 2 } / \lambda _ { 1 } | ^ { k }$ where $\lambda _ { 2 }$ is the second largest eigenvalue in absolute value. Unless there is a sufficient magnitude gap between these two eigenvalues, the power method is very slow. Moreover, it does not take advantage of “prior experience.” After k steps with initial vector $v ^ { ( 0 ) }$ , it has visited the directions defined by the vectors $A v ^ { ( 0 ) } , \ldots , \hat { A ^ { k } v ^ { ( 0 ) } }$ . However, instead of searching the span of these vectors for an optimal estimate of $x _ { 1 }$ , it settles for $A ^ { k } v ^ { ( 0 ) }$ . The method of orthogonal iteration with Ritz acceleration (§8.3.7) addresses some of these concerns, but it too has a certain disregard for prior iterates. What we need is a method that “learns from experience” and takes advantage of all previously computed matrix-vector products. The Lanczos method fits the bill.

# 10.1.1 Krylov Subspaces

The derivation of the Lanczos process can proceed in several ways. So that its remarkable convergence properties do not come as a complete surprise, we motivate the method by considering the optimization of the Rayleigh quotient

$$
r (x) = \frac {x ^ {T} A x}{x ^ {T} x}, \qquad x \neq 0.
$$

Recall from Theorem 8.1.2 that the maximum and minimum values of $r ( x )$ are $\lambda _ { 1 } ( A )$ and $\lambda _ { n } ( A )$ , respectively. Suppose $\{ q _ { i } \} \subseteq \mathbb { R } ^ { n }$ is a sequence of orthonormal vectors and define the scalars $M _ { k }$ and $m _ { k }$ by

$$
M _ {k} = \lambda_ {1} (Q _ {k} ^ {T} A Q _ {k}) = \max _ {y \neq 0} \frac {y ^ {T} (Q _ {k} ^ {T} A Q _ {k}) y}{y ^ {T} y} = \max _ {\| y \| _ {2} = 1} r (Q _ {k} y) \leq \lambda_ {1} (A),
$$

$$
m _ {k} = \lambda_ {k} (Q _ {k} ^ {T} A Q _ {k}) = \min _ {y \neq 0} \frac {y ^ {T} (Q _ {k} ^ {T} A Q _ {k}) y}{y ^ {T} y} = \min _ {\| y \| _ {2} = 1} r (Q _ {k} y) \geq \lambda_ {n} (A),
$$

where $Q _ { k } \ = \ \left[ \ q _ { 1 } \ | \cdot \cdot \cdot | \ q _ { k } \ \right]$ . Since

$$
\operatorname{ran} \left(Q _ {1}\right) \subset \operatorname{ran} \left(Q _ {2}\right) \subset \dots \subset \operatorname{ran} \left(Q _ {n}\right) = \mathbb {R} ^ {n}
$$

it follows that

$$
\begin{array}{l} M _ {1} \leq M _ {2} \leq \dots \leq M _ {n} = \lambda_ {1} (A), \\ m _ {1} \geq m _ {2} \geq \dots \geq m _ {n} = \lambda_ {n} (A). \\ \end{array}
$$

Thus, the proposed optimization framework will ultimately converge. However, the challenge is to choose the $q \mathrm { - }$ vectors in such a way that $M _ { k }$ and $m _ { k }$ are high-quality estimates well before $k$ equals $n .$ .

Searching for a good $q _ { k }$ prompts consideration of the gradient:

$$
\nabla r (x) = \frac {2}{x ^ {T} x} (A x - r (x) x). \tag {10.1.1}
$$

Suppose $u _ { k } \in \mathsf { s p a n } \{ q _ { 1 } , \dots , q _ { k } \}$ satisfies $M _ { k } = r ( u _ { k } )$ . If $\nabla r ( u _ { k } ) = 0$ , then $( r ( u _ { k } ) , u _ { k } )$ i s an eigenpair of A. If not, then from the standpoint of making $M _ { k + 1 }$ as large as possible it makes sense to choose the next trial vector $q _ { k + 1 }$ so that

$$
\nabla r (u _ {k}) \in \operatorname{span} \{q _ {1}, \dots , q _ {k + 1} \}. \tag {10.1.2}
$$

This is because $r ( x )$ increases most rapidly in the direction of the gradient $\nabla r ( x )$ . The strategy will guarantee that $M _ { k + 1 }$ is greater than $M _ { k }$ , hopefully by a significant amount. Likewise, if $v _ { k } \in { \mathsf { s p a n } } \{ q _ { 1 } , \dots , q _ { k } \}$ satisfies $r ( v _ { k } ) = m _ { k }$ , then it makes sense to require

$$
\nabla r (v _ {k}) \in \operatorname{span} \{q _ {1}, \dots , q _ {k + 1} \} \tag {10.1.3}
$$

since $r ( x )$ decreases most rapidly in the direction of $- \nabla r ( x )$ .

Note that for any $\boldsymbol { x } \in \mathbb { R } ^ { n }$ we have

$$
\nabla r (x) \in \operatorname{span} \{x, A x \}.
$$

Since the vectors $u _ { k }$ and $v _ { k }$ each belong to span $\{ q _ { 1 } , \ldots , q _ { k } \}$ , it follows that the inclusions (10.1.2) and (10.1.3) are satisfied if

$$
\operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\} = \operatorname{span} \left\{q _ {1}, A q _ {1}, \dots , A ^ {k - 1} q _ {1} \right\}.
$$

This suggests we choose $q _ { k + 1 }$ so that

$$
\operatorname{span} \left\{q _ {1}, \dots , q _ {k + 1} \right\} = \operatorname{span} \left\{q _ {1}, A q _ {1}, \dots , A ^ {k - 1} q _ {1}, A ^ {k} q _ {1} \right\}
$$

and thus we are led to the problem of computing orthonormal bases for the Krylov subspaces

$$
\mathcal {K} (A, q _ {1}, k) = \operatorname{span} \left\{q _ {1}, A q _ {1}, \dots , A ^ {k - 1} q _ {1} \right\}.
$$

These are just the range spaces of the Krylov matrices

$$
K (A, q _ {1}, k) = \left[ q _ {1} \mid A q _ {1} \mid A ^ {2} q _ {1} \mid \dots \mid A ^ {k - 1} q _ {1} \right]
$$

that we introduced in §8.3.2. Note that $\kappa ( A , q _ { 1 } , k )$ is precisely the subspace that the power method “overlooks” since it merely searches in the direction of $A ^ { k - 1 } q _ { 1 }$ .

# 10.1.2 Tridiagonalization

In order to generate an orthonormal basis for a Krylov subspace we exploit the connection between the tridiagonalization of A and the QR factorization of $K ( A , q _ { 1 } , n )$ . Recall from §8.3.2 that if $\bar { Q } ^ { T } A Q = T$ is tridiagonal and $Q Q ^ { T } = I _ { n }$ , then

$$
K (A, q _ {1}, n) = Q Q ^ {T} K (A, q _ {1}, n) = Q \left[ e _ {1} \mid T e _ {1} \mid T ^ {2} e _ {1} \mid \dots \mid T ^ {n - 1} e _ {1} \right]
$$

is the QR factorization of $K ( A , q _ { 1 } , n )$ where $e _ { 1 }$ and $q _ { 1 }$ are respectively the first columns of $I _ { n }$ and $Q .$ . Thus, the columns of Q can effectively be generated by tridiagonalizing A with an orthogonal matrix whose first column is $q _ { 1 }$ .

Householder tridiagonalization, discussed in §8.3.1, can be adapted for this purpose. However, this approach is impractical if A is large and sparse because Householder similarity updates almost always destroy sparsity. As a result, unacceptably large, dense matrices arise during the reduction. This suggests that we try to compute the elements of the tridiagonal matrix $T = Q ^ { T } A Q$ directly. Toward that end, designate the columns of $Q$ by

$$
Q = \left[ q _ {1} \mid \dots \mid q _ {n} \right]
$$

and the components of T by

$$
T = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \beta_ {n - 1} \\ 0 & \dots & & \beta_ {n - 1} & \alpha_ {n} \end{array} \right].
$$

Equating columns in $A Q = Q T$ , we conclude that

$$
A q _ {k} = \beta_ {k - 1} q _ {k - 1} + \alpha_ {k} q _ {k} + \beta_ {k} q _ {k + 1}, \quad (\beta_ {0} q _ {0} \equiv 0),
$$

for $k = 1 { : } n - 1$ . The orthonormality of the q-vectors implies

$$
\alpha_ {k} = q _ {k} ^ {T} A q _ {k}.
$$

(Another way to see this is that $T _ { i j } = q _ { i } ^ { T } A q _ { j } . )$ Moreover, if we define the vector $r _ { k }$ by

$$
r _ {k} = (A - \alpha_ {k} I) q _ {k} - \beta_ {k - 1} q _ {k - 1}
$$

and if it is nonzero, then

$$
q _ {k + 1} = r _ {k} / \beta_ {k}
$$

where

$$
\beta_ {k} = \pm \| r _ {k} \| _ {2}.
$$

If $r _ { k } = 0$ , then the iteration breaks down but (as we shall see) not without the acquisition of valuable invariant subspace information.

By properly sequencing the above formulae and assuming that $q _ { 1 } \in \mathbb { R } ^ { n }$ is a given unit vector, we obtain what may be regarded as “Version $0 ^ { \dag }$ of the Lanczos iteration.

Algorithm 10.1.1 (Lanczos Tridiagonalization) Given a symmetric matrix $A \in \mathbb { R } ^ { n \times n }$ and a unit 2-norm vector $q _ { 1 } \in \mathbb { R } ^ { n }$ , the following algorithm computes a matrix $Q _ { k } =$ $[ q _ { 1 } | \dots | q _ { k } ]$ with orthonormal columns and a tridiagonal matrix $T _ { k } \in \mathbb { R } ^ { k \times k }$ so that $A Q _ { k } \ = \ Q _ { k } T _ { k }$ . The diagonal and superdiagonal entries of $T _ { k }$ are $\alpha _ { 1 } , \ldots , \alpha _ { k }$ and $\beta _ { 1 } , \ldots , \beta _ { k - 1 }$ respectively. The integer k satisfies $1 \leq k \leq n$ .

$$
\begin{array}{l} k = 0, \beta_ {0} = 1, q _ {0} = 0, r _ {0} = q _ {1} \\ q _ {k + 1} = r _ {k} / \beta_ {k} \\ k = k + 1 \\ \alpha_ {k} = q _ {k} ^ {T} A q _ {k} \\ r _ {k} = (A - \alpha_ {k} I) q _ {k} - \beta_ {k - 1} q _ {k - 1} \\ \beta_ {k} = \parallel r _ {k} \parallel_ {2} \\ \end{array}
$$

end

There is no loss of generality in choosing $\beta _ { k }$ to be positive. The $q _ { k }$ vectors are called Lanczos vectors. It is important to mention that there are better ways numerically to organize the computation of the Lanczos vectors than Algorithm 10.1.1. See §10.3.1.

# 10.1.3 Termination and Error Bounds

The Lanczos iteration halts before complete tridiagonalization if $q _ { 1 }$ is contained in a proper invariant subspace. This is one of several mathematical properties of the method that we summarize in the following theorem.

Theorem 10.1.1. The Lanczos iteration (Algorithm 10.1.1) runs until $k = m$ , where

$$
m = \operatorname{rank} (K (A, q _ {1}, n)).
$$

Moreover, for k = 1:m we have

$$
A Q _ {k} = Q _ {k} T _ {k} + r _ {k} e _ {k} ^ {T} \tag {10.1.4}
$$

where $Q _ { k } = \left[ q _ { 1 } | \cdots | q _ { k } \right]$ has orthonormal columns that span $\displaystyle { \mathcal { K } } ( A , q _ { 1 } , k ) , e _ { k } = I _ { n } ( : , k )$ ,

and

$$
T _ {k} = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \beta_ {k - 1} \\ 0 & \dots & & \beta_ {k - 1} & \alpha_ {k} \end{array} \right]. \tag {10.1.5}
$$

Proof. The proof is by induction on k. It clearly holds if $k = 1$ . Suppose for some $k > 1$ that the iteration has produced $Q _ { k } \ = \ \left[ \ q _ { 1 } \ | \cdot \cdot \cdot | \ q _ { k } \ \right]$ with orthonormal columns such that

$$
\operatorname{ran} (Q _ {k}) = \mathcal {K} (A, q _ {1}, k).
$$

It is easy to see from Algorithm 10.1.1 that equation (10.1.4) holds and so

$$
Q _ {k} ^ {T} A Q _ {k} = T _ {k} + Q _ {k} ^ {T} r _ {k} e _ {k} ^ {T}. \tag {10.1.6}
$$

Suppose i and j are integers that satisfy $1 \leq i \leq j \leq k$ . From the equation

$$
q _ {j} ^ {T} A q _ {i} = q _ {j} ^ {T} \left(\beta_ {i - 1} q _ {i - 1} + \alpha_ {i} q _ {i} + \beta_ {i} q _ {i + 1}\right) = \beta_ {i - 1} q _ {j} ^ {T} q _ {i - 1} + \alpha_ {i} q _ {j} ^ {T} q _ {i} + \beta_ {i} q _ {j} ^ {T} q _ {i + 1}
$$

and the induction assumption $Q _ { k } ^ { T } Q _ { k } = I _ { k }$ , we see that

$$
q _ {i} ^ {T} A q _ {j} = q _ {j} ^ {T} A q _ {i} = \left\{ \begin{array}{l l} 0, & \text {if} i <   j - 1, \\ \beta_ {j - 1}, & \text {if} i = j - 1, \\ \alpha_ {j}, & \text {if} i = j. \end{array} \right.
$$

It follows that $Q _ { k } ^ { T } A Q _ { k } = T _ { k }$ and so from (10.1.6) we have $Q _ { k } ^ { T } r _ { k } = 0$

If $r _ { k } \neq 0$ , then $q _ { k + 1 } = r _ { k } / \parallel r _ { k } \parallel _ { 2 }$ is orthogonal to $q _ { 1 } , \ldots , q _ { k }$ . It follows that $q _ { k + 1 } \notin \mathcal { K } ( A , q _ { 1 } , k )$ and

$$
q _ {k + 1} \in \operatorname{span} \left\{A q _ {k}, q _ {k}, q _ {k - 1} \right\} \subseteq \mathcal {K} (A, q _ {1}, k + 1).
$$

Thus, $Q _ { k + 1 } ^ { T } Q _ { k + 1 } \ = \ I _ { k + 1 }$ and

$$
\operatorname{ran} (Q _ {k + 1}) = \mathcal {K} (A, q _ {1}, k + 1).
$$

On the other hand, if $r _ { k } = 0$ , then $A Q _ { k } = Q _ { k } T _ { k }$ . This says that ra $\mathsf { \Omega } _ { 1 } ( Q _ { k } ) = \mathcal { K } ( A , q _ { 1 } , k )$ is invariant for A and so $k = m = { \mathsf { r a n k } } ( K ( A , q _ { 1 } , n ) )$ .

To encounter a zero $\beta _ { k }$ in the Lanczos iteration is a welcome event in that it signals the computation of an exact invariant subspace. However, valuable approximate invariant subspace information tends to emerge long before the occurrence of a small $\beta .$ Apparently, more information can be extracted from the tridiagonal matrix $T _ { k }$ and the Krylov subspace spanned by the columns of $Q _ { k }$ .

# 10.1.4 Ritz Approximations

Recall from §8.1.4 that if S is a subspace of $\mathbb { R } ^ { n }$ , then with respect to $S$ we say that $( \theta , y )$ is a Ritz pair for $A \in \mathbb { R } ^ { n \times n }$ if $w ^ { T } ( A y - \theta y ) = 0$ for all $w \in S$ . If $S = \mathcal { K } ( A , q _ { 1 } , k )$ , then the Lanczos process can be used to compute the associated Ritz values and vectors. Suppose

$$
S _ {k} ^ {T} T _ {k} S _ {k} = \Theta_ {k} = \operatorname{diag} \left(\theta_ {1}, \dots , \theta_ {k}\right) \tag {10.1.7}
$$

is a Schur decomposition of the tridiagonal matrix $T _ { k }$ . If

$$
Y _ {k} = \left[ y _ {1} \mid \dots \mid y _ {k} \right] = Q _ {k} S _ {k} \in \mathbb {R} ^ {n \times k},
$$

then for $i = 1 { : } k$ it follows that $( \theta _ { i } , y _ { i } )$ is a Ritz pair because

$$
Q _ {k} ^ {T} (A Y _ {k} - Y _ {k} \Theta_ {k}) = (Q _ {k} ^ {T} A Q _ {k}) S _ {k} - Q _ {k} ^ {T} (Q _ {k} S _ {k}) \Theta_ {k} = T _ {k} S _ {k} - S _ {k} \Theta_ {k} = 0.
$$

Two theorems in §8.1 concern Ritz approximation and are of interest to us in the Lanczos setting. Theorem 8.1.14 tells us that the problem of minimizing $\parallel A Q _ { k } - Q _ { k } B \parallel _ { 2 }$ over all k-by-k matrices B is solved by setting $B = T _ { k } = Q _ { k } ^ { T } A Q _ { k }$ . Thus, the $\theta _ { i }$ are the eigenvalues of a “best possible matrix” that happens to be tridiagonal. Theorem 8.1.15 can be used to provide a bound for $\parallel A y _ { i } - \theta _ { i } y _ { i } \parallel _ { 2 }$ . However, we can actually do better. Using (10.1.6) we have

$$
A y _ {i} - \theta_ {i} y _ {i} = (A Q _ {k} - Q _ {k} T _ {k}) S _ {k} e _ {i} = r _ {k} (e _ {k} ^ {T} S _ {k} e _ {i})
$$

from which it follows that

$$
\left\| A y _ {i} - \theta_ {i} y _ {i} \right\| _ {2} = \left| \beta_ {k} \right| \left| s _ {k i} \right|. \tag {10.1.8}
$$

Note that since $S _ { k }$ is orthogonal, $| s _ { k i } | \le 1$ .

We can use (10.1.8) to obtain a computable error bound. If E is the rank-1 matrix

$$
E = - s _ {k i} \cdot r _ {k} y _ {i} ^ {T},
$$

then

$$
(A + E) y _ {i} = \theta_ {i} y _ {i}.
$$

It follows from Corollary 8.1.6 that

$$
\min _ {\mu \in \lambda (A)} | \theta_ {i} - \mu | \leq | \beta_ {k} | | s _ {k i} |
$$

for $i = 1 k$

Golub (1974) describes the construction of a more informative rank-1 perturbation E. Use Lanczos tridiagonalization to compute $A Q _ { k } = Q _ { k } T _ { k } + r _ { k } e _ { k } ^ { T }$ and then set $E = \tau w w ^ { T }$ , where $\tau = \pm 1$ and $w = a q _ { k } + b r _ { k }$ . It follows that

$$
(A + E) Q _ {k} = Q _ {k} (T _ {k} + \tau a ^ {2} e _ {k} e _ {k} ^ {T}) + (1 + \tau a b) r _ {k} e _ {k} ^ {T}.
$$

${ \mathrm { I f ~ } } 0 = 1 + \tau a b .$ , then

$$
\bar {T} _ {k} = T _ {k} + \tau a ^ {2} e _ {k} e _ {k} ^ {T}
$$

is a tridiagonal matrix whose eigenvalues are also eigenvalues for $A + E$ . Using Theorem 8.1.8, it can be shown that the interval $[ \lambda _ { i } ( \tilde { T } _ { k } ) , \lambda _ { i - 1 } ( \tilde { T } _ { k } ) ]$ contains an eigenvalue of A for $i = 2 { : } k$ . These bracketing intervals depend on the choice of $\tau a ^ { 2 }$ . Suppose we have an approximate eigenvalue λ of A. One possibility is to choose $\tau a ^ { 2 }$ so that

$$
\det (\tilde {T} _ {k} - \lambda I _ {k}) = (\alpha_ {k} + \tau a ^ {2} - \lambda) p _ {k - 1} (\lambda) - \beta_ {k - 1} ^ {2} p _ {k - 2} (\lambda) = 0
$$

where the polynomials $p _ { i } ( x ) = \mathsf { d e t } ( T _ { i } - x I _ { i } )$ can be evaluated at λ using the three-term recurrence (8.4.2). (This assumes that $p _ { k - 1 } ( \lambda ) \neq 0 . )$ The idea of characterizing an approximate eigenvalue λ as an exact eigenvalue of a nearby matrix $A + E$ is discussed in Lehmann (1963) and Householder (1968).

# 10.1.5 Convergence Theory

The preceding discussion indicates how eigenvalue estimates can be obtained via the Lanczos process, but it reveals nothing about the approximation quality of $T _ { k } \mathrm { { ^ { * } s } }$ eigenvalues as a function of k. Results of this variety have been developed by Kaniel, Paige, Saad, and others and the following theorem is a sample from this body of research.

Theorem 10.1.2. Let A be an n-by-n symmetric matrix with Schur decomposition

$$
Z ^ {T} A Z = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right), \quad \lambda_ {1} \geq \dots \geq \lambda_ {n}, \quad Z = \left[ z _ {1} \mid \dots \mid z _ {n} \right]. \tag {10.1.9}
$$

Suppose k steps of the Lanczos iteration (Algorithm $\it { 1 0 . 1 . 1 ) }$ are performed and that $T _ { k }$ is the tridiagonal matrix (10.1.5). If $\ d \theta _ { 1 } = \lambda _ { 1 } ( T _ { k } ) \ d t$ , then

$$
\lambda_ {1} \geq \theta_ {1} \geq \lambda_ {1} - (\lambda_ {1} - \lambda_ {n}) \left(\frac {\tan (\phi_ {1})}{c _ {k - 1} (1 + 2 \rho_ {1})}\right) ^ {2}
$$

where cos $\left( \phi _ { 1 } \right) = \left| q _ { 1 } ^ { T } z _ { 1 } \right|$ ,

$$
\rho_ {1} = \frac {\lambda_ {1} - \lambda_ {2}}{\lambda_ {2} - \lambda_ {n}}, \tag {10.1.10}
$$

and $c _ { k - 1 } ( x )$ is the Chebyshev polynomial of degree $k - 1$ .

Proof. From Theorem 8.1.2, we have

$$
\theta_ {1} = \max _ {y \neq 0} \frac {y ^ {T} T _ {k} y}{y ^ {T} y} = \max _ {y \neq 0} \frac {(Q _ {k} y) ^ {T} A (Q _ {k} y)}{(Q _ {k} y) ^ {T} (Q _ {k} y)} = \max _ {0 \neq w \in \mathcal {K} (A, q _ {1}, k)} \frac {w ^ {T} A w}{w ^ {T} w}.
$$

Since $\lambda _ { 1 }$ is the maximum of $w ^ { T } A w / w ^ { T } \boldsymbol { \imath }$ w over all nonzero $w ,$ it follows that $\theta _ { 1 } \leq \lambda _ { 1 }$ . To obtain the lower bound for $\theta _ { 1 }$ , note that

$$
\theta_ {1} = \max _ {p \in \mathbb {P} _ {k - 1}} \frac {q _ {1} ^ {T} p (A) A p (A) q _ {1}}{q _ {1} ^ {T} p (A) ^ {2} q _ {1}},
$$

where $\mathbb { P } _ { k - 1 }$ is the set of degree-(k−1) polynomials and $p ( x )$ is the amplifying polynomial. Given the eigenvector expansion $q _ { 1 } = d _ { 1 } z _ { 1 } + \cdot \cdot \cdot + d _ { n } z _ { n }$ where $d _ { i } = q _ { 1 } ^ { T } z _ { i }$ , it follows that

$$
\frac {q _ {1} ^ {T} p (A) A p (A) q _ {1}}{q _ {1} ^ {T} p (A) ^ {2} q _ {1}} = \frac {\sum_ {i = 1} ^ {n} d _ {i} ^ {2} p (\lambda_ {i}) ^ {2} \lambda_ {i}}{\sum_ {i = 1} ^ {n} d _ {i} ^ {2} p (\lambda_ {i}) ^ {2}} \geq \frac {\lambda_ {1} d _ {1} ^ {2} p (\lambda_ {1}) ^ {2} + \lambda_ {n} \delta^ {2}}{d _ {1} ^ {2} p (\lambda_ {1}) ^ {2} + \delta^ {2}} = \lambda_ {1} - \frac {(\lambda_ {1} - \lambda_ {n}) \delta^ {2}}{d _ {1} ^ {2} p (\lambda_ {1}) ^ {2} + \delta^ {2}}
$$

where

$$
\delta^ {2} = \sum_ {i = 2} ^ {n} d _ {i} ^ {2} p (\lambda_ {i}) ^ {2}.
$$

If the polynomial $p$ has the property that it is large at $x = \lambda _ { 1 }$ compared to its value at $\lambda _ { 2 } , \ldots , \lambda _ { n }$ , then we get a better lower bound for the Ritz value $\theta _ { 1 }$ . This is the act of finding an amplifying polynomial and a good choice is to set

$$
p (x) = c _ {k - 1} \left(- 1 + 2 \frac {x - \lambda_ {n}}{\lambda_ {2} - \lambda_ {n}}\right)
$$

where $c _ { k - 1 } ( z )$ is the (k−1)st Chebyshev polynomial generated via the recursion

$$
c _ {k} (z) = 2 z c _ {k - 1} (z) - c _ {k - 2} (z), \quad c _ {0} = 1, c _ {1} = z.
$$

These polynomials are bounded by unity on [−1, 1], but grow very rapidly outside this interval. By defining $p ( x )$ this way, it follows that $| p ( \lambda _ { i } ) | \leq 1$ for $i = 2 { : } n$ and $p ( \lambda _ { 1 } ) =$ $c _ { k - 1 } ( 1 + 2 \rho _ { 1 } )$ where $\rho _ { 1 }$ is defined by (10.1.10). Thus,

$$
\delta^ {2} \leq \sum_ {i = 2} ^ {n} d _ {i} ^ {2} = 1 - d _ {1} ^ {2}
$$

and so

$$
\theta_ {1} \geq \lambda_ {1} - \left(\lambda_ {1} - \lambda_ {n}\right) \frac {1 - d _ {1} ^ {2}}{d _ {1} ^ {2}} \frac {1}{\left(c _ {k - 1} \left(1 + 2 \rho_ {1}\right)\right) ^ {2}}.
$$

The desired lower bound is obtained by noting that tan $( \phi _ { 1 } ) ^ { 2 } = ( 1 - d _ { 1 } ^ { 2 } ) / d _ { 1 } ^ { 2 }$ .

An analogous result pertaining to $T _ { k }$ ’s smallest eigenvalue is an easy corollary.

Corollary 10.1.3. Using the same notation as in the theorem, if $\theta _ { k } = \lambda _ { k } ( T _ { k } )$ , then

$$
\lambda_ {n} \leq \theta_ {k} \leq \lambda_ {n} + (\lambda_ {1} - \lambda_ {n}) \left(\frac {\tan (\phi_ {n})}{c _ {k - 1} (1 + 2 \rho_ {n})}\right) ^ {2}
$$

where

$$
\rho_ {n} = \frac {\lambda_ {n - 1} - \lambda_ {n}}{\lambda_ {1} - \lambda_ {n - 1}}
$$

and cos $\left( \phi _ { n } \right) = q _ { 1 } ^ { T } z _ { n }$ .

Proof. Apply Theorem 10.1.2 with A replaced by −A.

The key idea in the proof of Theorem 10.1.2 is to take the amplifying polynomial $p ( x )$ 号 to be the translated Chebyshev polynomial, for then $p ( A ) q _ { 1 }$ amplifies the component of $q _ { 1 }$ in the direction of the eigenvector $z _ { 1 }$ . A similar idea can be used to obtain bounds for an interior Ritz value $\theta _ { i }$ . However, the results are not as satisfactory because the new amplifying polynomial involves the product of the Chebyshev polynomial $c _ { k - i }$ and the polynomial $\left( x - \lambda _ { 1 } \right) \cdot \cdot \cdot \left( x - \lambda _ { i - 1 } \right)$ . For details, see Kaniel (1966) and Paige (1971) and also Saad (1980), who improved the bounds. The main theorem is as follows.

Theorem 10.1.4. Using the same notation as Theorem 10.1.2, if $1 \leq i \leq k$ and $\theta _ { i } = \lambda _ { i } ( T _ { k } )$ , then

$$
\lambda_ {i} \geq \theta_ {i} \geq \lambda_ {i} - (\lambda_ {1} - \lambda_ {n}) \left(\frac {\kappa_ {i} \tan (\phi_ {i})}{c _ {k - i} (1 + 2 \rho_ {i})}\right) ^ {2}
$$

where

$$
\rho_ {i} = \frac {\lambda_ {i} - \lambda_ {i + 1}}{\lambda_ {i + 1} - \lambda_ {n}}, \qquad \kappa_ {i} = \prod_ {j = 1} ^ {i - 1} \frac {\theta_ {j} - \lambda_ {n}}{\theta_ {j} - \lambda_ {i}}, \qquad \cos (\phi_ {i}) = | q _ {1} ^ {T} z _ {i} |.
$$

Proof. See Saad (NMLE, p. 201).

Because of the $\kappa _ { i }$ factor and the reduced degree of the amplifying Chebyshev polynomial, it is clear that the bounds deteriorate as i increases.

# 10.1.6 The Power Method versus the Lanczos Method

It is instructive to compare $\theta _ { 1 }$ with the corresponding power method estimate of $\lambda _ { 1 }$ . (See §8.2.1.) For clarity, assume $\lambda _ { 1 } \geq \cdot \cdot \cdot \geq \lambda _ { n } \geq 0$ in the Schur decomposition (10.1.7). After $k - 1$ power method steps applied to $q _ { 1 }$ , a vector is obtained in the direction of

$$
v = A ^ {k - 1} q _ {1} = \sum_ {i = 1} ^ {n} d _ {i} \lambda_ {i} ^ {k - 1} z _ {i}
$$

along with an eigenvalue estimate

$$
\gamma_ {1} = \frac {v ^ {T} A v}{v ^ {T} v}.
$$

By setting $p ( x ) = x ^ { k - 1 }$ in the proof of Theorem 10.1.2, it is easy to show that

$$
\lambda_ {1} \geq \gamma_ {1} \geq \lambda_ {1} - (\lambda_ {1} - \lambda_ {n}) \tan (\phi_ {1}) ^ {2} \left(\frac {\lambda_ {2}}{\lambda_ {1}}\right) ^ {2 (k - 1)}. \tag {10.1.11}
$$

Thus, we can compare the quality of the lower bounds for $\theta _ { 1 }$ and $\gamma _ { 1 }$ by comparing

$$
L _ {k - 1} \equiv \frac {1}{\left[ c _ {k - 1} \left(2 \frac {\lambda_ {1}}{\lambda_ {2}} - 1\right) \right] ^ {2}} \geq \frac {1}{\left[ c _ {k - 1} (1 + 2 \rho_ {1}) \right] ^ {2}}
$$

and

$$
R _ {k - 1} = \left(\frac {\lambda_ {2}}{\lambda_ {1}}\right) ^ {2 (k - 1)}.
$$

Figure 10.1.1 compares these quantities for various values of k and $\lambda _ { 2 } / \lambda _ { 1 }$ . The superiority of the Lanczos bound is self-evident. This is not a surprise since $\theta _ { 1 }$ is the maximum of $r ( x ) = x ^ { T } A x / x ^ { T } x$ over all of $\kappa ( A , q _ { 1 } , k )$ , while $\gamma _ { 1 } = r ( \boldsymbol { v } )$ for a particular v in $\kappa ( A , q _ { 1 } , k )$ , namely, $v = A ^ { k - 1 } q _ { 1 }$ .

<table><tr><td> $\lambda_1/\lambda_2$ </td><td>k=5</td><td>k=10</td><td>k=15</td><td>k=20</td><td>k=25</td></tr><tr><td>1.50</td><td> $\frac{1.1\times 10^{-4}}{3.9\times 10^{-2}}$ </td><td> $\frac{2.0\times 10^{-10}}{6.8\times 10^{-4}}$ </td><td> $\frac{3.9\times 10^{-16}}{1.2\times 10^{-5}}$ </td><td> $\frac{7.4\times 10^{-22}}{2.0\times 10^{-7}}$ </td><td> $\frac{1.4\times 10^{-27}}{3.5\times 10^{-9}}$ </td></tr><tr><td>1.10</td><td> $\frac{2.7\times 10^{-2}}{4.7\times 10^{-1}}$ </td><td> $\frac{5.5\times 10^{-5}}{1.8\times 10^{-1}}$ </td><td> $\frac{1.1\times 10^{-7}}{6.9\times 10^{-2}}$ </td><td> $\frac{2.1\times 10^{-10}}{2.7\times 10^{-2}}$ </td><td> $\frac{4.2\times 10^{-13}}{1.0\times 10^{-2}}$ </td></tr><tr><td>1.01</td><td> $\frac{5.6\times 10^{-1}}{9.2\times 10^{-1}}$ </td><td> $\frac{1.0\times 10^{-1}}{8.4\times 10^{-1}}$ </td><td> $\frac{1.5\times 10^{-2}}{7.6\times 10^{-1}}$ </td><td> $\frac{2.0\times 10^{-3}}{6.9\times 10^{-1}}$ </td><td> $\frac{2.8\times 10^{-4}}{6.2\times 10^{-1}}$ </td></tr></table>

Figure 10.1.1. $L _ { k - 1 } / R _ { k - 1 }$

# Problems

P10.1.1 Suppose $A \in \mathbb { R } ^ { n \times n }$ is skew-symmetric. Derive a Lanczos-like algorithm for computing a skew-symmetric tridiagonal matrix $T _ { m }$ such that $A Q _ { m } = Q _ { m } T _ { m }$ , where $Q _ { m } ^ { T } Q _ { m } = I _ { m }$ .

P10.1.2 Let $A \in \mathbb { R } ^ { n \times n }$ be symmetric and define $r ( x ) = x ^ { T } A x / x ^ { T } x$ . Suppose $S \subseteq \mathbb { R } ^ { n }$ is a subspace with the property that $x \in S$ implies $\nabla r ( x ) \in S$ . Show that S is invariant for A.

P10.1.3 Show that if a symmetric matrix $A \in \mathbb { R } ^ { n \times n }$ has a multiple eigenvalue, then the Lanczos process terminates prematurely.

P10.1.4 Show that the index m in Theorem 10.1.1 is the dimension of the smallest invariant subspace for A that contains $q _ { 1 }$ .

P10.1.5 Let $A \in \mathbb { R } ^ { n \times n }$ be symmetric and consider the problem of determining an orthonormal sequence $q _ { 1 } , q _ { 2 } , . . .$ . with the property that once $Q _ { k } = [ \ : q _ { 1 } \ : | \cdots | \ : q _ { k } \ : ]$ is known, $q _ { k + 1 }$ is chosen so as to minimize $\begin{array} { r c l } { \mu _ { k } } & { = } & { \parallel ( I - Q _ { k + 1 } Q _ { k + 1 } ^ { T } ) A Q _ { k } \parallel _ { _ { F } } } \end{array}$ . Show that if span $\{ q _ { 1 } , . . . , q _ { k } \} = { \mathcal { K } } ( A , q _ { 1 } , k )$ , then it is possible to choose $q _ { k + 1 }$ so $\mu _ { k } = 0$ . Explain how this optimization problem leads to the Lanczos iteration.

P10.1.6 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and that we wish to compute its largest eigenvalue. Let η be an approximate eigenvector and set $\ d \alpha ~ = ~ \eta ^ { T } \ d A \eta / \eta ^ { T } \eta$ and $z ~ = ~ A \eta - \alpha \eta$ . (a) Show that the interval $[ \alpha - \delta , \alpha + \delta ]$ must contain an eigenvalue of A where $\delta = \parallel z \parallel _ { 2 } / \parallel \eta \parallel _ { 2 }$ . (b) Consider the new approximation $\bar { \eta } = a \eta + b z$ and determine the scalars a and b so that $\bar { \alpha } = { \bar { \eta } ^ { T } A \bar { \eta } } / { \bar { \eta } ^ { T } \bar { \eta } }$ is maximized. (c) Relate the above computations to the first two steps of the Lanczos process.

P10.1.7 Suppose $T \in \mathbb { R } ^ { n \times n }$ is tridiagonal and symmetric and that $v \in \mathbb { R } ^ { n }$ . Show how the Lanczos process can be used (in principle) to compute an orthogonal $Q \in \mathbb { R } ^ { n \times n }$ in $O ( n ^ { 2 } )$ flops such that $Q ^ { T } ( T + v v ^ { T } ) Q = \tilde { T }$ is also tridiagonal.

# Notes and References for 10.1

Detailed treatments of the symmetric Lanczos algorithm may be found in Parlett (SEP) and Meurant (LCG). The classic reference for the Lanczos method is:

C. Lanczos (1950). “An Iteration Method for the Solution of the Eigenvalue Problem of Linear Differential and Integral Operators,” J. Res. Nat. Bur. Stand. 45, 255–282.

For details about the convergence of the Ritz values, see:

S. Kaniel (1966). “Estimates for Some Computational Techniques in Linear Algebra,” Math. Comput. 20, 369–378.

C.C. Paige (1971). “The Computation of Eigenvalues and Eigenvectors of Very Large Sparse Matrices,” PhD thesis, University of London.

Y. Saad (1980). “On the Rates of Convergence of the Lanczos and the Block Lanczos Methods,” SIAM J. Numer. Anal. 17, 687–706.

The connections between Lanczos tridiagonalization, orthogonal polynomials, and the theory of moments are discussed in:

N.J. Lehmann (1963). “Optimale Eigenwerteinschliessungen,” Numer. Math. 5, 246–272.

A.S. Householder (1968). “Moments and Characteristic Roots II,” Numer. Math. 11, 126–128.

G.H. Golub (1974). “Some Uses of the Lanczos Algorithm in Numerical Linear Algebra,” in Topics in Numerical Analysis, J.J.H. Miller (ed.), Academic Press, New York.

C.C. Paige, B.N. Parlett, and H.A. van der Vorst (1995). “Approximate Solutions and Eigenvalue Bounds from Krylov Subspaces,” Numer. Lin. Alg. Applic. 2, 115–133.
