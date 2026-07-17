# Chapter 10

# Large Sparse Eigenvalue Problems

10.1 The Symmetric Lanczos Process   
10.2 Lanczos, Quadrature, and Approximation   
10.3 Practical Lanczos Procedures   
10.4 Large Sparse SVD Frameworks   
10.5 Krylov Methods for Unsymmetric Problems   
10.6 Jacobi-Davidson and Related Methods

The Lanczos process computes a sequence of partial tridiagonalizations that are orthogonally related to a given symmetric matrix A. It is of particular interest if A is large and sparse because, instead of updating A along the way as in the Householder method of §8.2, it simply relies on matrix-vector products. Equally important, information about A’s extremal eigenvalues tends to emerge fairly early during the iteration, making the method very useful in situations where just a few of A’s largest or smallest eigenvalues are desired, together with the corresponding eigenvectors.

The derivation and exact arithmetic attributes of the method are presented in §10.1, including its extraordinary convergence properties. Central to the discussion is the connection to an underlying Krylov subspace that is defined by the starting vector. In §10.2 we point out connections between Gauss quadrature and the Lanczos process that can be used to estimate expressions of the form $u ^ { T } f ( A ) u$ where f (A) is a function of a large, sparse symmetric positive definite matrix A. Unfortunately, a “math book” implementation of the Lanczos method is practically useless because of roundoff error. This makes it necessary to enlist the help of various “workarounds,” which we describe in §10.3. A sparse SVD framework based on Golub-Kahan bidiagonalization is detailed in §10.4. We also introduce the idea of a randomized SVD. The last two sections deal with the more difficult unsymmetric problem. The Arnoldi iteration is a Krylov subspace iteration like Lanczos. To make it effective, it is necessary to extract valuable “restart information” from the Hessenberg matrix sequence that it produces. This is discussed in §10.5 together with a brief presentation of the unsymmetric Lanczos framework. In the last section we derive the Jacobi-Davidson method, which combines Newton ideas with Rayleigh-Ritz refinement.

# Reading Notes

Familiarity with Chapters 5, 7, and 8 is recommended. Within this chapter there are the following dependencies:

$$
\begin{array}{c c c c c c c c} \S 1 0. 1 & \to & \S 1 0. 3 & \to & \S 1 0. 5 & \to & \S 1 0. 6 \\ \downarrow & & \downarrow & & \\ \S 1 0. 2 & & \S 1 0. 4 & & \end{array}
$$

General references for this chapter include Parlett (SEP), Stewart (MAE), Watkins (MEP), Chatelin (EOM), Cullum and Willoughby (LALSE), Meurant (LCG), Saad (NMLE), Kressner (NMSE), and EIG TEMPLATES.

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

# 10.2 Lanczos, Quadrature, and Approximation

To deepen our understanding of the Lanczos process and to build an appreciation for its connections to other areas of applied mathematics, we consider an interesting approximation problem that has broad practical implications. Assume that A ∈ IRn×n $A \in \mathbb { R } ^ { n \times n }$ is a large, sparse, symmetric positive definite matrix whose eigenvalues reside in an interval $[ a , b ]$ . Let $f ( \lambda )$ be a given smooth function that is defined on $[ a , b ]$ . Given $u \in \mathbb { R } ^ { n }$ , our goal is to produce suitably tight lower and upper bounds b and B so that

$$
b \leq u ^ {T} \cdot f (A) \cdot u \leq B. \tag {10.2.1}
$$

In the approach we develop, the bounds are Gauss quadrature rule estimates of a certain integral and the evaluation of the rules requires the eigenvalues and eigenvectors of a Lanczos-produced tridiagonal matrix.

The $u ^ { T } f ( A ) u$ estimation problem has many applications throughout matrix computations. For example, suppose ˆx is an approximate solution to the symmetric positive definite system $A x = b$ and that we have computed the residual $r = b - A \hat { x }$ . Note that if $x _ { * } = A ^ { - 1 } b$ and $f ( \lambda ) = 1 / \lambda ^ { 2 }$ , then

$$
\| x _ {*} - \hat {x} \| _ {2} ^ {2} = (x _ {*} - \hat {x}) ^ {T} (x _ {*} - \hat {x}) = (A ^ {- 1} (b - A \hat {x})) ^ {T} (A ^ {- 1} (b - A \hat {x})) = r ^ {T} f (A) r.
$$

Thus, if we have a $u ^ { T } f ( A ) u$ estimation framework, then we can obtain $A x = b$ error bounds from residual bounds.

For an in-depth treatment of the material in this section, we refer the reader to the treatise by Golub and Meurant (2010). Our presentation is brief, informal, and stresses the linear algebra highlights.

# 10.2.1 Reformulation of the Problem

Without an integral in sight, it is mystifying as to why (10.2.1) involves quadrature at all. The key is to regard $u ^ { T } f ( A ) u$ as a Riemann-Stieltjes integral. In general, given a suitably nice integrand f(x) and weight function $w ( x )$ , the Riemann-Stieltjes integral

$$
I (f) = \int_ {a} ^ {b} f (x) d w (x)
$$

is a limit of sums of the form

$$
S _ {N} = \sum_ {\mu = 1} ^ {N} f (c _ {\mu}) (w (x _ {\mu}) - w (x _ {\mu + 1}))
$$

where $a = x _ { N } < \cdot \cdot \cdot < x _ { 1 } = b$ and $x _ { \mu + 1 } \leq c _ { \mu } \leq x _ { \mu }$ . Note that if w is piecewise constant on $[ a , b ]$ , then the only nonzero terms in $S _ { N }$ arise from subintervals that house a “w-jump.” For example, suppose $a = \lambda _ { n } < \lambda _ { 2 } < \cdots < \lambda _ { 1 } = b$ and that

$$
w (\lambda) = \left\{ \begin{array}{l l} w _ {n + 1} & \text { if } \lambda <   a, \\ w _ {\mu} & \text { if } \lambda_ {\mu} \leq \lambda <   \lambda_ {\mu - 1}, \\ w _ {1} & \text { if } b \leq \lambda , \end{array} \right. \quad \mu = 2: n, \tag {10.2.2}
$$

where $0 \leq w _ { n + 1 } \leq \cdot \cdot \cdot \leq w _ { 1 }$ . By considering the behavior of $S _ { N }$ as $N  \infty$ , we see that

$$
\int_ {a} ^ {b} f (\lambda) d w (\lambda) = \sum_ {\mu = 1} ^ {n} (w _ {\mu} - w _ {\mu + 1}) \cdot f (\lambda_ {\mu}). \tag {10.2.3}
$$

We are now set to explain why $u ^ { T } f ( A ) { \mathrm { : } }$ u is “secretly” a Riemann-Stieltjes integral. Let

$$
A = X \Lambda X ^ {T}, \quad \Lambda = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right), \tag {10.2.4}
$$

be a Schur decomposition of A with $\lambda _ { n } \leq \cdots \leq \lambda _ { 1 }$ . It follows that

$$
u ^ {T} f (A) u = (X ^ {T} u) ^ {T} \cdot f (\Lambda) \cdot (X ^ {T} u) = \sum_ {\mu = 1} ^ {n} [ X ^ {T} u ] _ {\mu} ^ {2} \cdot f (\lambda_ {\mu}).
$$

If we set

$$
w _ {\mu} = \left[ X ^ {T} u \right] _ {\mu} ^ {2} + \dots + \left[ X ^ {T} u \right] _ {n} ^ {2}, \quad \mu = 1: n + 1, \tag {10.2.5}
$$

in (10.2.2), then (10.2.3) becomes

$$
\int_ {a} ^ {b} f (\lambda) d w (\lambda) = \sum_ {\mu = 1} ^ {n} \left[ X ^ {T} u \right] _ {\mu} ^ {2} \cdot f \left(\lambda_ {\mu}\right) = u ^ {T} f (A) u. \tag {10.2.6}
$$

Our plan is to approximate this integral using Gauss quadrature.

# 10.2.2 Some Gauss-Type Quadrature Rules and Bounds

Given an accuracy-related parameter k, an interval $[ a , b ]$ , and a weight function $w ( \lambda )$ , a Gauss-type quadrature rule for the integral

$$
I (f) = \int_ {a} ^ {b} f (\lambda) d w (\lambda)
$$

involves a carefully constructed linear combination of f-evaluations across $[ a , b ]$ . The evaluation points (called nodes) and the coefficients (called weights) that define the linear combination are determined to make the rule correct for polynomials up to a certain degree that is related to k. Here are four examples:

1. Gauss. Compute weights $w _ { 1 } , \ldots , w _ { k }$ and nodes $t _ { 1 } , \ldots , t _ { k }$ so if

$$
I _ {G} (f) = \sum_ {i = 1} ^ {k} w _ {i} f (t _ {i}) \tag {10.2.7}
$$

then $I ( f ) = I _ { G } ( f )$ for all polynomials f that have degree $2 k - 1$ or less.

2. Gauss-Radau(a). Compute weights $w _ { a } , w _ { 1 } , \ldots , w _ { k }$ and nodes $t _ { 1 } , \ldots , t _ { k }$ so if

$$
I _ {G R (a)} (f) = w _ {a} f (a) + \sum_ {i = 1} ^ {k} w _ {i} f (t _ {i}) \tag {10.2.8}
$$

then $I ( f ) = I _ { G R ( a ) } ( f )$ for all polynomials f that have degree 2k or less.

3. Gauss-Radau(b). Compute weights $w _ { b } , w _ { 1 } , \dotsc , w _ { k }$ and nodes $t _ { 1 } , \ldots , t _ { k }$ so if

$$
I _ {G R (b)} (f) = w _ {b} f (b) + \sum_ {i = 1} ^ {k} w _ {i} f (t _ {i}) \tag {10.2.9}
$$

then $I ( f ) = I _ { G R ( b ) } ( f )$ for all polynomials f that have degree 2k or less.

4. Gauss-Lobatto. Compute weights $w _ { a } , w _ { b } , w _ { 1 } , \ldots , w _ { k }$ and nodes $t _ { 1 } , \ldots , t _ { k }$ so if

$$
I _ {G L} (f) = w _ {a} f (a) + w _ {b} f (b) + \sum_ {i = 1} ^ {k} w _ {i} f (t _ {i}) \tag {10.2.10}
$$

then $I ( f ) = I _ { G L } ( f )$ for all polynomials f that have degree $2 k + 1$ or less.

Each of these rules has a neatly specified error. It can be shown that

$$
\int_ {a} ^ {b} f (\lambda) d w (\lambda) = \left\{ \begin{array}{l l} I _ {G} (f) & + R _ {G} (f), \\ I _ {G R (a)} (f) & + R _ {G R (a)} (f), \\ I _ {G R (b)} (f) & + R _ {G R (b)} (f), \\ I _ {G L} (f) & + R _ {G G} (f), \end{array} \right.
$$

where

$$
R _ {G} (f) \quad = \frac {f ^ {(2 k)} (\eta)}{(2 n) !} \int_ {a} ^ {b} \left[ \prod_ {i = 1} ^ {k} (\lambda - t _ {i}) \right] ^ {2} d w (\lambda), \quad a <   \eta <   b,
$$

$$
R _ {G R (a)} (f) = \frac {f ^ {(2 k + 1)} (\eta)}{(2 k + 1) !} \int_ {a} ^ {b} (\lambda - a) \left[ \prod_ {i = 1} ^ {k} (\lambda - t _ {i}) \right] ^ {2} d w (\lambda), \quad a <   \eta <   b,
$$

$$
R _ {G R (b)} (f) = \frac {f ^ {(2 k + 1)} (\eta)}{(2 k + 1) !} \int_ {a} ^ {b} (\lambda - b) \left[ \prod_ {i = 1} ^ {k} (\lambda - t _ {i}) \right] ^ {2} d w (\lambda), \quad a <   \eta <   b,
$$

$$
R _ {G L} (f) = \frac {f ^ {(2 k + 2)} (\eta)}{(2 k + 2) !} \int_ {a} ^ {b} (\lambda - a) (\lambda - b) \left[ \prod_ {i = 1} ^ {k} (\lambda - t _ {i}) \right] ^ {2} d w (\lambda), \quad a <   \eta <   b.
$$

If the derivative in the remainder term does not change sign across [a, b], then the rule can be used to produce a bound. For example, if $f ( \lambda ) = 1 / \lambda ^ { 2 }$ and $0 < a < b$ , then $f ^ { ( 2 k ) }$ is positive, $f ^ { ( 2 k + 1 ) }$ is negative, and we have

$$
I _ {G} (f) \leq \int_ {a} ^ {b} f (\lambda) d w (\lambda) \leq I _ {G R (a)} (f).
$$

With this strategy, we can produce lower and upper bounds by selecting and evaluating the right rule. For this to be practical, the behavior of f ’s higher derivatives must be known and the required rules must be computable.

# 10.2.3 The Tridiagonal Connection

It turns out that the evaluation of a given Gauss quadrature rule involves a tridiagonal matrix and its eigenvalues and eigenvectors. To develop a strategy that is based upon this connection, we need three facts about orthogonal polynomials and Gauss quadrature.

Fact 1. Given $[ a , b ]$ and $w ( \lambda )$ , there is a sequence of polynomials $p _ { 0 } ( \lambda ) , p _ { 1 } ( \lambda ) , . . .$ . that satisfy

$$
\int_ {a} ^ {b} p _ {i} (\lambda) \cdot p _ {j} (\lambda) \cdot d w (\lambda) = \left\{ \begin{array}{l l} 1 & \mathrm{if} i = j, \\ 0 & \mathrm{if} i \neq j, \end{array} \right.
$$

with the property that the degree of $p _ { k } ( \cdot )$ is k for $k \geq 0$ . The polynomials are unique up to a factor of ±1 and they satisfy a 3-term recurrence

$$
\gamma_ {k} p _ {k} (\lambda) = (\lambda - w _ {k}) p _ {k - 1} (\lambda) - \gamma_ {k - 1} p _ {k - 2} (\lambda)
$$

where $p _ { - 1 } ( \lambda ) \equiv 0$ and $p _ { 0 } ( \lambda ) \equiv 1$ .

Fact 2. The zeros of $p _ { k } ( \lambda )$ are the eigenvalues of the tridiagonal matrix

$$
T _ {k} = \left[ \begin{array}{c c c c c} \omega_ {1} & \gamma_ {1} & 0 & \dots & 0 \\ \gamma_ {1} & \omega_ {2} & \ddots & & \vdots \\ 0 & \ddots & \ddots & \ddots & 0 \\ \vdots & & \ddots & \omega_ {k - 1} & \gamma_ {k - 1} \\ 0 & \dots & 0 & \gamma_ {k - 1} & \omega_ {k} \end{array} \right].
$$

Since the $\gamma _ { i }$ are nonzero, it follows from Theorem 8.4.1 that the eigenvalues are distinct.

Fact 3. If

$$
S ^ {T} T _ {k} S = \mathrm{diag} (\theta_ {1}, \dots , \theta_ {k}) \tag {10.2.11}
$$

is a Schur decomposition of $T _ { k } .$ , then the nodes and weights for the Gauss rule (10.2.7) are given by $t _ { i } = \theta _ { i }$ and $w _ { i } = s _ { 1 i } ^ { 2 }$ for i = 1:k. In other words,

$$
I _ {G} (f) = \sum_ {i = 1} ^ {k} s _ {1 i} ^ {2} \cdot f (\theta_ {i}). \tag {10.2.12}
$$

Thus, the only remaining issue is how to construct $T _ { k }$ so that it defines a Gauss rule for (10.2.6).

# 10.2.4 Gauss Quadrature via Lanczos

We show that if we apply the symmetric Lanczos process (Algorithm 10.1.1) with starting vector $q _ { 1 } = u / \parallel u \parallel _ { 2 }$ , then the tridiagonal matrices that the method generates are exactly what we need to compute $I _ { G } ( f )$ .

We first link the Lanczos process to a sequence of orthogonal polynomials. Recall from §10.1.1 that the kth Lanczos vector $q _ { k }$ is in the Krylov subspace $\kappa ( A , q _ { 1 } , k )$ . It follows that $q _ { k } ~ = ~ p _ { k } ( A ) q _ { 1 }$ for some degree-k polynomial. From Algorithm 10.1.1 we know that

$$
\beta_ {k} q _ {k + 1} = (A - \alpha_ {k} I) q _ {k} - \beta_ {k - 1} q _ {k - 1}
$$

where $\beta _ { 0 } q _ { 0 } \equiv 0$ and so

$$
\beta_ {k} p _ {k + 1} (A) q _ {1} = (A - \alpha_ {k} I) p _ {k} (A) q _ {1} - \beta_ {k - 1} p _ {k - 1} (A) q _ {1}.
$$

From this we conclude that the polynomials satisfy a 3-term recurrence:

$$
\beta_ {k} p _ {k + 1} (\lambda) = (\lambda - \alpha_ {k}) p _ {k} (\lambda) - \beta_ {k - 1} ^ {2} p _ {k - 1} (\lambda). \tag {10.2.13}
$$

These polynomials are orthogonal with respect to the $u ^ { T } f ( A ) u$ weight function defined in (10.2.5). To see this, note that

$$
\begin{array}{l} \int_ {a} ^ {b} p _ {i} (\lambda) p _ {j} (\lambda) d w (\lambda) = \sum_ {\mu = 1} ^ {n} [ X ^ {T} u ] _ {\mu} ^ {2} \cdot p _ {i} (\lambda_ {\mu}) \cdot p _ {j} (\lambda_ {\mu}) \\ = (X ^ {T} u) ^ {T} (p _ {i} (\Lambda) \cdot p _ {j} (\Lambda)) \cdot (X ^ {T} u) \\ = u ^ {T} \left(X \cdot p _ {i} (\Lambda) \cdot X ^ {T}\right) \left(X \cdot p _ {j} (\Lambda) \cdot X ^ {T}\right) u \\ = u ^ {T} \left(p _ {i} (A) p _ {j} (A)\right) u \\ = (p _ {i} (A) u) ^ {T} (p _ {j} (A) u) = \| u \| _ {2} ^ {2} q _ {i} ^ {T} q _ {j} = 0. \\ \end{array}
$$

Coupled with (10.2.13) and Facts 1-3, this result tells us that we can generate an approximation $\sigma = I _ { G } ( f )$ to $u ^ { T } f ( A ) u$ as follows:

Step 1: With starting vector $q _ { 1 } = u / \parallel u \parallel _ { 2 }$ , use the Lanczos process to compute the partial tridiagonalization $A Q _ { k } = Q _ { k } T _ { k } + r _ { k } e _ { k } ^ { T }$ . (See (10.1.4).)

Step 2: Compute the Schur decomposition $S ^ { T } T _ { k } S = \operatorname { d i a g } ( \theta _ { 1 } , . . . , \theta _ { k } )$

Step 3: Set $\sigma = s _ { 1 1 } ^ { 2 } f ( \theta _ { 1 } ) + \cdot \cdot \cdot + s _ { 1 k } ^ { 2 } f ( \theta _ { k } )$ .

See Golub and Welsch (1969) for a more rigorous derivation of this procedure.

# 10.2.5 Computing the Gauss-Radau Rule

Recall from (10.2.1) that we are interested in upper and lower bounds. In light of our remarks at the end of §10.2.2, we need techniques for evaluating other Gauss quadrature rules. By way of illustration, we show how to compute $I _ { G R ( a ) }$ defined in (10.2.8). Guided by Gauss quadrature theory, we run the Lanczos process for k steps as if we were setting out to compute $I _ { G } ( f )$ . We then must determine $\tilde { \alpha } _ { k + 1 }$ so that if

$$
\tilde {T} _ {k + 1} = \left[ \begin{array}{c c c c c c c} \alpha_ {1} & \beta_ {1} & 0 & \dots & 0 & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots & \vdots \\ 0 & \ddots & \ddots & \ddots & \vdots & \\ \vdots & & \ddots & \alpha_ {k - 1} & \beta_ {k - 1} & 0 \\ 0 & \dots & & \beta_ {k - 1} & \alpha_ {k} & \beta_ {k} \\ \hline 0 & \dots & \dots & 0 & \beta_ {k} & \tilde {\alpha} _ {k + 1} \end{array} \right]
$$

then $a \in \lambda ( \tilde { T } _ { k + 1 } )$ . By considering the top and bottom halves of the equation

$$
\tilde {T} _ {k + 1} \left[ \begin{array}{c} x \\ - 1 \end{array} \right] = a \left[ \begin{array}{c} x \\ - 1 \end{array} \right], \qquad x \in \mathbb {R} ^ {k},
$$

it is easy to verify that $\tilde { \alpha } _ { k + 1 } = a + \beta _ { k + 1 } ^ { 2 } e _ { k } ^ { T } ( T _ { k } - a I _ { k } ) ^ { - 1 } e _ { k }$ works.

# 10.2.6 The Overall Framework

All the necessary tools are now available to obtain sufficiently accurate upper and bounds in (10.2.1). At the bottom of the loop in Algorithm 10.1.1, we use the current tridiagonal (or an augmented version) to compute the nodes and weights for the lower bound rule. The rule is evaluated to obtain b. Likewise, we use the current tridiagonal (or an augmented version) to compute the nodes and weights for the upper bound rule. The rule is evaluated to obtain B. The while loop in Algorithm 10.1.1 can obviously be redesigned to terminate as soon as $B - b$ is sufficiently small.

# Problems

P10.2.1 The Chebyschev polynomials are generated by the recursion $p _ { k } ( x ) = 2 x p _ { k - 1 } ( x ) - p _ { k - 2 } ( x )$ and are orthonormal with respect to $w ( x ) = ( 1 - x ^ { 2 } ) ^ { - 1 / 2 }$ across [−1, 1]. What are the zeros of $p _ { k } ( x ) ?$

P10.2.2 Following the strategy used in §10.2.5, show how to compute $I _ { G R ( b ) }$ and $I _ { G L } ( f )$ .

# Notes and References for §10.2

For complete coverage of the Gauss quadrature/tridiagonal/Lanczos connection, see:

G.H. Golub and G. Meurant (2010). Matrices, Moments, and Quadrature with Applications, Princeton University Press, Princeton, NJ.

Research in this area has a long history:

G.H. Golub (1962). “Bounds for Eigenvalues of Tridiagonal Symmetric Matrices Computed by the LR Method,” Math. Comput. 16, 438–445.

G.H. Golub and J.H. Welsch (1969). “Calculation of Gauss Quadrature Rules,” Math. Comput. 23, 221–230.

G.H. Golub (1974). “Bounds for Matrix Moments,” Rocky Mountain J. Math. 4, 207–211.

C. de Boor and G.H. Golub (1978). “The Numerically Stable Reconstruction of a Jacobi Matrix from Spectral Data,” Lin. Alg. Applic. 21, 245–260.

J. Kautsky and G.H. Golub (1983). “On the Calculation of Jacobi Matrices,” Lin. Alg. Applic. 52/53, 439–455.

M. Berry and G.H. Golub (1991). “Estimating the Largest Singular Values of Large Sparse Matrices via Modified Moments,” Numer. Algs. 1, 353–374.   
D.P. Laurie (1996). “Anti-Gaussian Quadrature Rules,” Math. Comput. 65, 739–747.   
Z. Bai and G.H. Golub (1997). “Bounds for the Trace of the Inverse and the Determinant of Symmetric Positive Definite Matrices,” Annals Numer. Math. 4, 29–38.   
M. Benzi and G.H. Golub (1999). “Bounds for the Entries of Matrix Functions with Applications to Preconditioning,” BIT 39, 417–438.   
D. Calvetti, G. H. Golub, W. B. Gragg, and L. Reichel (2000). “Computation of Gauss–Kronrod Quadrature Rules,” Math. Comput. 69, 1035–1052.   
D.P. Laurie (2001). “Computation of Gauss-Type Quadrature Formulas,” J. Comput. Appl. Math. 127, 201–217.

# 10.3 Practical Lanczos Procedures

Rounding errors greatly affect the behavior of the Lanczos iteration. The basic difficulty is caused by loss of orthogonality among the Lanczos vectors, a phenomenon that muddies the issue of termination and complicates the relationship between A’s eigenvalues and those of the tridiagonal matrices $T _ { k }$ . This troublesome feature, coupled with the advent of Householder’s perfectly stable method of tridiagonalization, explains why the Lanczos algorithm was disregarded by numerical analysts during the 1950’s and 1960’s. However, the pressure to solve large, sparse eigenproblems coupled with the computational insights set forth by Paige (1971) changed all that. With many fewer than n iterations typically required to get good approximate extremal eigenvalues, the Lanczos method became attractive as a sparse matrix technique rather than as a competitor of the Householder approach.

Successful implementation of the Lanczos iteration involves much more than a simple encoding of Algorithm 10.1.1. In this section we present some of the ideas that have been proposed to make the Lanczos procedure viable in practice.

# 10.3.1 Required Storage and Work

With careful overwriting in Algorithm 10.1.1 and exploitation of the formula

$$
\alpha_ {k} = q _ {k} ^ {T} (A q _ {k} - \beta_ {k - 1} q _ {k - 1}),
$$

the whole Lanczos process can be implemented with just a pair of n-vectors:

$$
w = q _ {1}, v = A w, \alpha_ {1} = w ^ {T} v, v = v - \alpha_ {1} w, \beta_ {1} = \| v \| _ {2}, k = 1
$$

while $\beta _ { k } \neq 0$

for i = 1:n

$$
t = w _ {i}, w _ {i} = v _ {i} / \beta_ {k}, v _ {i} = - \beta_ {k} t
$$

end (10.3.1)

$$
v = v + A w
$$

$$
k = k + 1, \alpha_ {k} = w ^ {T} v, v = v - \alpha_ {k} w, \beta_ {k} = \parallel v \parallel_ {2}
$$

end

At the end of the loop body, the array w houses $q _ { k }$ and v houses the residual vector $r _ { k } = A q _ { k } - \alpha _ { k } q _ { k } - \beta _ { k - 1 } q _ { k - 1 }$ . See Paige (1972) for a discussion of various Lanczos implementations and their numerical properties. Note that A is not modified during the entire process and that is what makes the procedure so useful for large sparse matrices.

If A has an average of ν nonzeros per row, then approximately $( 2 \nu + 8 ) n$ flops are involved in a single Lanczos step. Upon termination the eigenvalues of $T _ { k }$ can be found using the symmetric tridiagonal QR algorithm or any of the special methods of §8.5 such as bisection. The Lanczos vectors are generated in the n-vector w. If eigenvectors are required, then the Lanczos vectors must be saved. Typically, they are stored in secondary memory units.

# 10.3.2 Roundoff Properties

The development of a practical, easy-to-use Lanczos tridiagonalization process requires an appreciation of the fundamental error analyses of Paige (1971, 1976, 1980). An examination of his results is the best way to motivate the several modified Lanczos procedures of this section.

After j steps of the iteration we obtain the matrix of computed Lanczos vectors $\hat { Q } _ { k } = \left[ \begin{array} { l } { \hat { q } _ { 1 } } \end{array} \right| \cdot \cdot \cdot \left| \begin{array} { l } { \hat { q } _ { k } } \end{array} \right]$ and the associated tridiagonal matrix

$$
\hat {T} _ {k} = \left[ \begin{array}{c c c c c} \hat {\alpha} _ {1} & \hat {\beta} _ {1} & & \dots & 0 \\ \hat {\beta} _ {1} & \hat {\alpha} _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \hat {\beta} _ {k - 1} \\ 0 & \dots & & \hat {\beta} _ {k - 1} & \hat {\alpha} _ {k} \end{array} \right].
$$

Paige (1971, 1976) shows that if $\hat { r } _ { k }$ is the computed analog of $r _ { k }$ , then

$$
A \hat {Q} _ {k} = \hat {Q} _ {k} \hat {T} _ {k} + \hat {r} _ {k} e _ {k} ^ {T} + E _ {k} \tag {10.3.2}
$$

where

$$
\| E _ {k} \| _ {2} \approx \mathbf {u} \| A \| _ {2}. \tag {10.3.3}
$$

This shows that the equation $A Q _ { k } = Q _ { k } T _ { k } + r _ { k } e _ { k } ^ { T }$ is satisfied to working precision.

Unfortunately, the picture is much less rosy with respect to the orthogonality among the $\hat { q } _ { i }$ . (Normality is not an issue. The computed Lanczos vectors essentially have unit length.) If ${ \hat { \beta } } _ { k } = \mathsf { f l } ( \parallel { \hat { r } } _ { k } \parallel _ { 2 } )$ and we compute $\hat { q } _ { k + 1 } = \mathsf { f l } \left( \hat { r } _ { k } / \hat { \beta } _ { k } \right)$ , then a simple analysis shows that

$$
\hat {\beta} _ {k} \hat {q} _ {k + 1} \approx \hat {r} _ {k} + w _ {k}
$$

where

$$
\| w _ {k} \| _ {2} \approx \mathbf {u} \| \hat {r} _ {k} \| _ {2} \approx \mathbf {u} \| A \| _ {2}.
$$

Thus, we may conclude that

$$
| \hat {q} _ {k + 1} ^ {T} \hat {q} _ {i} | \approx \frac {| \hat {r} _ {k} ^ {T} \hat {q} _ {i} | + \mathbf {u} | | A | | _ {2}}{| \hat {\beta} _ {k} |}
$$

for $i = 1 { : } k$ . In other words, significant departures from orthogonality can be expected when $\hat { \beta } _ { k }$ is small, even in the ideal situation where $\hat { r } _ { k } ^ { T } \hat { Q } _ { k }$ is zero. A small $\hat { \beta } _ { k }$ implies cancellation in the computation of $\hat { r } _ { k }$ . We stress that loss of orthogonality is due to one or several such cancellations and is not the result of the gradual accumulation of roundoff error.

Further details of the Paige analysis are given shortly. Suffice it to say now that loss of orthogonality always occurs in practice and with it, an apparent deterioration in the quality of $\hat { T } _ { k } ^ { \phantom { } } \mathrm { ' s }$ eigenvalues. This can be quantified by combining (10.3.2) with Theorem 8.1.16. In particular, if we set

$$
F _ {1} = \hat {r} _ {k} e _ {k} ^ {T} + E _ {k}, \qquad X _ {1} = \hat {Q} _ {k}, \qquad S = \hat {T} _ {k},
$$

in that theorem and assume that

$$
\tau = \parallel \hat {Q} _ {k} ^ {T} \hat {Q} _ {k} - I _ {k} \parallel_ {2}
$$

satisfies $\tau < 1$ , then there exist eigenvalues $\mu _ { 1 } , \dots , \mu _ { k } \in \lambda ( A )$ such that

$$
| \mu_ {i} - \lambda_ {i} (T _ {k}) | \leq \sqrt {2} \left(\| \hat {r} _ {k} \| _ {2} + \| E _ {k} \| _ {2} + \tau (2 + \tau) \| A \| _ {2}\right)
$$

for $i = 1 { : } k$ . An obvious way to control the τ factor is to orthogonalize each newly computed Lanczos vector against its predecessors. This leads directly to our first “practical” Lanczos procedure.

# 10.3.3 Lanczos with Complete Reorthogonalization

Let $r _ { 0 } , \ldots , r _ { k - 1 } \in \mathbb { R } ^ { n }$ be given and suppose that Householder matrices $H _ { 0 } , \ldots , H _ { k - 1 }$ have been computed such that $( H _ { 0 } \cdot \cdot \cdot H _ { k - 1 } ) ^ { T } \left[ \begin{array} { l } { r _ { 0 } } \end{array} \right| \cdot \cdot \cdot \left| \begin{array} { l } { r _ { k - 1 } } \end{array} \right]$ is upper triangular. Let ${ \left[ \begin{array} { l } { q _ { 1 } } \end{array} \right| } \cdots { \left| \begin{array} { l } { q _ { k } } \end{array} \right] }$ denote the first k columns of the Householder product $\left( H _ { 0 } \cdot \cdot \cdot H _ { k - 1 } \right)$ . Now suppose that we are given a vector $r _ { k } \in \mathbb { R } ^ { n }$ and wish to compute a unit vector $q _ { k + 1 }$ in the direction of

$$
w = r _ {k} - \sum_ {i = 1} ^ {k} \left(q _ {i} ^ {T} r _ {k}\right) q _ {i} \in \operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\} ^ {\perp}.
$$

If a Householder matrix $H _ { k }$ is determined so $( H _ { 0 } \cdot \cdot \cdot H _ { k } ) ^ { T } \left[ \begin{array} { l } { r _ { 0 } } \end{array} | \cdot \cdot \cdot | \ r _ { k } \right]$ is upper triangular, then it follows that column $( k + 1 )$ of $H _ { 0 } \cdots H _ { k }$ is the desired unit vector.

If we incorporate these Householder computations into the Lanczos process, then we can produce Lanczos vectors that are orthogonal to machine precision:

$r _ { 0 } = q _ { 1 }$ (given unit vector)

Determine Householder $H _ { 0 }$ so $H _ { 0 } r _ { 0 } = e _ { 1 }$ .

for $k = 1 { : } n - 1$

$$
\begin{array}{l} \alpha_ {k} = q _ {k} ^ {T} A q _ {k} \\ r _ {k} = (A - \alpha_ {k} I) q _ {k} - \beta_ {k - 1} q _ {k - 1}, \quad (\beta_ {0} q _ {0} \equiv 0) \tag {10.3.4} \\ \end{array}
$$

$$
w = (H _ {k - 1} \dots H _ {0}) r _ {k}
$$

$\mathrm { D e t e r m i n e ~ H o u s e h o l d e r ~ } H _ { k } \mathrm { ~ s o ~ } H _ { k } w = [ w _ { 1 } , \dots , w _ { k } , \beta _ { k } , 0 , \dots , 0 ] ^ { T } .$

$$
q _ {k + 1} = H _ {0} \dots H _ {k} e _ {k + 1}
$$

end

This is an example of a complete reorthorgonalization Lanczos scheme. The idea of using Householder matrices to enforce orthogonality appears in Golub, Underwood, and Wilkinson (1972). That the computed $\hat { q } _ { i }$ in (10.3.4) are orthogonal to working precision follows from the roundoff properties of Householder matrices. Note that by virtue of the definition of $q _ { k + 1 }$ , it makes no difference if $\beta _ { k } = 0$ . For this reason, the algorithm may safely run until $k = n - 1$ . (However, in practice one would terminate for a much smaller value of k.)

Of course, in any implementation of (10.3.4), one stores the Householder vectors $v _ { k }$ and never explicitly forms the corresponding matrix product. Since we have $H _ { k } ( 1 { : } k , 1 { : } k ) = I _ { k }$ there is no need to compute the first k components of the vector w in (10.3.4) since we do not use them. (Ideally they are zero.)

Unfortunately, these economies make but a small dent in the computational overhead associated with complete reorthogonalization. The Householder calculations increase the work in the kth Lanczos step by $O ( k n )$ flops. Moreover, to compute $q _ { k + 1 }$ , the Householder vectors associated with $H _ { 0 } , \ldots , H _ { k }$ must be accessed. For large n and k, this usually implies a prohibitive level of memory traffic.

Thus, there is a high price associated with complete reorthogonalization. Fortunately, there are more effective courses of action to take, but these require a greater understanding of just how orthogonality is lost.

# 10.3.4 Selective Reorthogonalization

A remarkable, ironic consequence of the Paige (1971) error analysis is that loss of orthogonality goes hand in hand with convergence of a Ritz pair. To be precise, suppose the symmetric QR algorithm is applied to $\hat { T } _ { k }$ and renders computed Ritz values $\hat { \theta } _ { 1 } , \ldots , \hat { \theta } _ { k }$ and a nearly orthogonal matrix of eigenvectors $\hat { S } _ { k } = ( \hat { s } _ { p q } )$ . If

$$
\hat {Y} _ {k} = \left[ \begin{array}{c c c c} \hat {y} _ {1} & \dots & \hat {y} _ {k} \end{array} \right] = \mathsf {f l} (\hat {Q} _ {k} \hat {S} _ {k}),
$$

then it can be shown that for $i = 1 { : } k$ we have

$$
\left| \hat {q} _ {k + 1} ^ {T} \hat {y} _ {i} \right| \approx \frac {\mathbf {u} \| A \| _ {2}}{\left| \hat {\beta} _ {k} \right| \left| \hat {s} _ {k i} \right|} \tag {10.3.5}
$$

and

$$
\| A \hat {y} _ {i} - \hat {\theta} _ {i} \hat {y} _ {i} \| _ {2} \approx | \hat {\beta} _ {k} | | \hat {s} _ {k i} |. \tag {10.3.6}
$$

That is, the most recently computed Lanczos vector $\hat { q } _ { k + 1 }$ tends to have a nontrivial and unwanted component in the direction of any converged Ritz vector. Consequently, instead of orthogonalizing $\hat { q } _ { k + 1 }$ against all of the previously computed Lanczos vectors, we can achieve the same effect by orthogonalizing it against the much smaller set of converged Ritz vectors.

The practical aspects of enforcing orthogonality in this way are discussed in Parlett and Scott (1979). In their scheme, known as selective reorthogonalization, a computed Ritz pair $\{ \hat { \theta } , \hat { y } \}$ is called “good” if it satisfies

$$
\| A \hat {y} - \hat {\theta} \hat {y} \| _ {2} \leq \sqrt {\mathbf {u}} \| A \| _ {2}.
$$

As soon as $\hat { q } _ { k + 1 }$ is computed, it is orthogonalized against each good Ritz vector. This is much less costly than complete reorthogonalization, since, at least at first, there are many fewer good Ritz vectors than Lanczos vectors.

One way to implement selective reorthogonalization is to diagonalize $\hat { T } _ { k }$ at each step and then examine the $\hat { s } _ { k i }$ in light of (10.3.5) and (10.3.6). A more efficient approach for large k is to estimate the loss-of-orthogonality measure  $I _ { k } - \hat { Q } _ { k } ^ { T } \hat { Q } _ { k } \parallel _ { 2 }$ using the following result.

Lemma 10.3.1. Suppose $S _ { + } = [ S d ]$ where $S \in \mathbb { R } ^ { n \times k }$ and $d \in \mathbb { R } ^ { n }$ . If

$$
\left\| I _ {k} - S ^ {T} S \right\| _ {2} \leq \mu \quad | 1 - d ^ {T} d | \leq \delta ,
$$

then

$$
\| I _ {k + 1} - S _ {+} ^ {T} S _ {+} \| _ {2} \leq \mu_ {+}
$$

where

$$
\mu_ {+} = \frac {1}{2} \left(\mu + \delta + \sqrt {(\mu - \delta) ^ {2} + 4 \| S ^ {T} d \| _ {2} ^ {2}}\right).
$$

Proof. See Kahan and Parlett (1974) or Parlett and Scott (1979).

Thus, if we have a bound for $\parallel I _ { k } - \hat { Q } _ { k } ^ { T } \hat { Q } _ { k } \parallel _ { 2 }$ , then by applying the lemma with $S = \hat { Q } _ { k }$ and $d = \hat { q } _ { k + 1 }$ we can generate a bound for $\parallel I _ { k + 1 } - \hat { Q } _ { k + 1 } ^ { T } \hat { Q } _ { k + 1 } \parallel _ { 2 }$ . (In this case $\delta \approx \mathbf { u }$ and we assume that $\hat { q } _ { k + 1 }$ has been orthogonalized against the set of currently good Ritz vectors.) It is possible to estimate the norm of $\hat { Q } _ { k } ^ { T } \hat { q } _ { k + 1 }$ from a simple recurrence that spares one the need to access $\hat { q } _ { 1 } , \dots , \hat { q } _ { k }$ . The overhead is minimal, and when the bounds signal loss of orthogonality, it is time to contemplate the enlargement of the set of good Ritz vectors. Then and only then is $\hat { T } _ { k }$ diagonalized.

# 10.3.5 The Ghost Eigenvalue Problem

Considerable effort has been spent in trying to develop a workable Lanczos procedure that does not involve any kind of orthogonality enforcement. Research in this direction focuses on the problem of “ghost” eigenvalues. These are multiple eigenvalues of $\hat { T } _ { k }$ that correspond to simple eigenvalues of A. They arise because the iteration essentially restarts itself when orthogonality to a converged Ritz vector is lost. (By way of analogy, consider what would happen during orthogonal iteration (8.2.8) if we “forgot” to orthogonalize.)

The problem of identifying ghost eigenvalues and coping with their presence is discussed by Cullum and Willoughby (1979) and Parlett and Reid (1981). It is a particularly pressing problem in those applications where all of $A \mathrm { { } i \mathrm { { s } } }$ eigenvalues are desired, for then the above orthogonalization procedures are expensive to implement.

Difficulties with the Lanczos iteration can be expected even if A has a genuinely multiple eigenvalue. This follows because the $\hat { T } _ { k }$ are unreduced, and unreduced tridiagonal matrices cannot have multiple eigenvalues. The next practical Lanczos procedure that we discuss attempts to circumvent this difficulty.

# 10.3.6 Block Lanczos Algorithm

Just as the simple power method has a block analogue in simultaneous iteration, so does the Lanczos algorithm have a block version. Suppose $n = r p$ and consider the

decomposition

$$
Q ^ {T} A Q = \bar {T} = \left[ \begin{array}{c c c c c} M _ {1} & B _ {1} ^ {T} & & \dots & 0 \\ B _ {1} & M _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & B _ {r - 1} ^ {T} \\ 0 & \dots & & B _ {r - 1} & M _ {r} \end{array} \right] \tag {10.3.7}
$$

where

$$
Q = \left[ X _ {1} \mid \dots \mid X _ {r} \right], \quad X _ {i} \in \mathbb {R} ^ {n \times p},
$$

is orthogonal, each $M _ { i } \in \mathbb { R } ^ { p \times p }$ , and each $B _ { i } \in \mathbb { R } ^ { p \times p }$ is upper triangular. Comparison of blocks in $A Q = Q { \bar { T } }$ shows that

$$
A X _ {k} = X _ {k - 1} B _ {k - 1} ^ {T} + X _ {k} M _ {k} + X _ {k + 1} B _ {k}
$$

for $k = 1 { : } r$ assuming $X _ { 0 } B _ { 0 } ^ { T } \equiv 0$ and $X _ { r + 1 } B _ { r } \equiv 0$ . From the orthogonality of $Q$ we have

$$
M _ {k} = X _ {k} ^ {T} A X _ {k}
$$

for $k = 1 { : } r$ . Moreover, if we define

$$
R _ {k} = A X _ {k} - X _ {k} M _ {k} - X _ {k - 1} B _ {k - 1} ^ {T} \in \mathbb {R} ^ {n \times p},
$$

then

$$
X _ {k + 1} B _ {k} = R _ {k}
$$

is a QR factorization of $R _ { k }$ . These observations suggest that the block tridiagonal matrix $\bar { T }$ in (10.3.7) can be generated as follows:

$$
X _ {1} \in \mathbb {R} ^ {n \times p} \text {   given   with   } X _ {1} ^ {T} X _ {1} = I _ {p}
$$

$$
M _ {1} = X _ {1} ^ {T} A X _ {1}
$$

$$
\text { for } k = 1: r - 1 \tag {10.3.8}
$$

$$
R _ {k} = A X _ {k} - X _ {k} M _ {k} - X _ {k - 1} B _ {k - 1} ^ {T} \quad (X _ {0} B _ {0} ^ {T} \equiv 0)
$$

$$
X _ {k + 1} B _ {k} = R _ {k} \quad \text {(QR factorization of R_ {k})}
$$

$$
M _ {k + 1} = X _ {k + 1} ^ {T} A X _ {k + 1}
$$

end

At the beginning of the kth pass through the loop we have

$$
A \left[ X _ {1} \mid \dots \mid X _ {k} \right] = \left[ X _ {1} \mid \dots \mid X _ {k} \right] \bar {T} _ {k} + R _ {k} \left[ 0 \mid \dots \mid 0 \mid I _ {p} \right], \tag {10.3.9}
$$

where

$$
\bar {T} _ {k} = \left[ \begin{array}{c c c c c} M _ {1} & B _ {1} ^ {T} & & \dots & 0 \\ B _ {1} & M _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & B _ {k - 1} ^ {T} \\ 0 & \dots & & B _ {k - 1} & M _ {k} \end{array} \right].
$$

Using an argument similar to the one used in the proof of Theorem 10.1.1, we can show that the $X _ { k }$ are mutually orthogonal provided none of the $R _ { k }$ is rank-deficient. However if rank $( R _ { k } ) < p$ for some k, then it is possible to choose the columns of $X _ { k + 1 }$ such that $X _ { k + 1 } ^ { T } X _ { i } = 0$ , for $i = 1 { : } k$ . See Golub and Underwood (1977).

Because $\hat { T } _ { k }$ has bandwidth p, it can be efficiently reduced to tridiagonal form using an algorithm of Schwartz (1968). Once tridiagonal form is achieved, the Ritz values can be obtained via the symmetric QR algorithm or any of the special methods of §8.4. In order to decide intelligently when to use block Lanczos, it is necessary to understand how the block dimension affects convergence of the Ritz values. The following generalization of Theorem 10.1.2 sheds light on this issue.

Theorem 10.3.2. Let A be an n-by-n symmetric matrix with Schur decomposition

$$
Z ^ {T} A Z = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right), \quad \lambda_ {1} \geq \dots \geq \lambda_ {n}, \quad Z = \left[ z _ {1} \mid \dots \mid z _ {n} \right].
$$

Let $\mu _ { 1 } \geq \cdots \geq \mu _ { p }$ be the p largest eigenvalues of the matrix $\hat { T } _ { k }$ obtained after k steps of $( 1 0 . 3 . 8 )$ . Suppose $Z _ { 1 } = { \left[ \begin{array} { l } { z _ { 1 } } \end{array} | \cdot \cdot \cdot \ \right] } z _ { p } \ ]$ and

$$
0 <   \cos (\phi_ {p}) = \sigma_ {p} (Z _ {1} ^ {T} X _ {1}),
$$

the smallest singular value of $Z _ { 1 } ^ { T } X _ { 1 }$ . Then for $i = 1 { : } p _ { : }$ ,

$$
\lambda_ {i} \geq \mu_ {i} \geq \lambda_ {i} - (\lambda_ {1} - \lambda_ {n}) \left(\frac {\tan (\theta_ {p})}{c _ {k - 1} (1 + 2 \rho_ {i})}\right) ^ {2}
$$

where

$$
\rho_ {i} = \frac {\lambda_ {i} - \lambda_ {p + 1}}{\lambda_ {p + 1} - \lambda_ {n}}
$$

and $c _ { k - 1 } ( z )$ is the Chebyshev polynomial of degree $k - 1$ .

Proof. See Underwood (1975). Compare with Theorem 10.1.2.

Analogous inequalities can be obtained for $\bar { T } _ { k } ^ { \ , } \mathrm { s }$ smallest eigenvalues by applying the theorem with A replaced by −A. Based on the theorem and scrutiny of (10.3.8), we conclude that

• the error bounds for the Ritz values improve with increased $p$   
• the amount of work required to compute $\hat { T } _ { k }$ ’s eigenvalues is proportional to $k p ^ { 2 }$   
• the block dimension should be at least as large as the largest multiplicity of any sought-after eigenvalue.

Determination of the block dimension in the face of these trade-offs is discussed in detail by Scott (1979). We mention that loss of orthogonality also plagues the block Lanczos algorithm. However, all of the orthogonality enforcement schemes described above can be extended to the block setting.

# 10.3.7 Block Lanczos Algorithm with Restarting

The block Lanczos algorithm (10.3.8) can be used in an iterative fashion to calculate selected eigenvalues of A. To fix ideas, suppose we wish to calculate the p largest eigenvalues. If $X _ { 1 } \in \mathbb { R } ^ { n \times p }$ is a given matrix having orthonormal columns, then it can be refined as follows:

Step 1. Generate $\boldsymbol { X } _ { 2 } , \ldots , \boldsymbol { X } _ { s } \in \mathbb { R } ^ { n \times p }$ via the block Lanczos algorithm.

Step 2. Form ${ \bar { T } } _ { s } = [ X _ { 1 } | \cdots | X _ { s } ] ^ { T } A [ X _ { 1 } | \cdots | X _ { s } ]$ , an sp-by-sp matrix that has bandwidth p.

Step 3. Compute an orthogonal matrix $U = \ [ \boldsymbol { u } _ { 1 } \vert \cdot \cdot \cdot \vert \ u _ { s p } ]$ such that $U ^ { T } \bar { T } _ { s } U =$ di $\arg ( \theta _ { 1 } , \ldots , \theta _ { s p } )$ with $\theta _ { 1 } \ge \cdots \ge \theta _ { s p }$ .

Step 4. Set $X _ { 1 } ^ { ( \mathrm { n e w } ) } = [ X _ { 1 } | \cdots | X _ { s } ] [ u _ { 1 } | \cdots | u _ { p } ]$

This is the block analog of the s-step Lanczos algorithm, which has been extensively analyzed by Cullum and Donath (1974) and Underwood (1975). The same idea can be used to compute several of A’s smallest eigenvalues or a mixture of both large and small eigenvalues. See Cullum (1978). The choice of the parameters s and p depends upon storage constraints as well as upon the block-size implications that we discussed above. The value of p can be diminished as the good Ritz vectors emerge. However, this demands that orthogonality to the converged vectors be enforced.

# Problems

P10.3.1 Rearrange (10.3.4) and (10.3.8) so that they require one matrix-vector product per iteration.

P10.3.2 If rank $( R _ { k } ) < p$ in (10.3.8), does it follow that ran $( \left[ X _ { 1 } \mid \cdots \mid X _ { k } \right] )$ contains an eigenvector of A?

# Notes and References for §10.3

The behavior of the Lanczos method in the presence of roundoff error was originally reported in:

C.C. Paige (1971). “The Computation of Eigenvalues and Eigenvectors of Very Large Sparse Matrices,” PhD thesis, University of London.

Important follow-up papers include:

C.C. Paige (1972). “Computational Variants of the Lanczos Method for the Eigenproblem,” J. Inst. Math. Applic. 10, 373–381.

C.C. Paige (1976). “Error Analysis of the Lanczos Algorithm for Tridiagonalizing a Symmetric Matrix,” J. Inst. Math. Applic. 18, 341–349.

C.C. Paige (1980). “Accuracy and Effectiveness of the Lanczos Algorithm for the Symmetric Eigenproblem,” Lin. Alg. Applic. 34, 235–258.

For additional analysis of the method, see Parlett (SEP), Meurant (LCG) as well as:

D.S. Scott (1979). “How to Make the Lanczos Algorithm Converge Slowly,” Math. Comput. 33, 239–247.

B.N. Parlett, H.D. Simon, and L.M. Stringer (1982). “On Estimating the Largest Eigenvalue with the Lanczos Algorithm,” Math. Comput. 38, 153–166.

B.N. Parlett and B. Nour-Omid (1985). “The Use of a Refined Error Bound When Updating Eigenvalues of Tridiagonals,” Lin. Alg. Applic. 68, 179–220.

J. Kuczy´nski and H. Wo´zniakowski (1992). “Estimating the Largest Eigenvalue by the Power and Lanczos Algorithms with a Random Start,” SIAM J. Matrix Anal. Applic. 13, 1094–1122.   
G. Meurant and Z. Strakos (2006). “The Lanczos and Conjugate Gradient Algorithms in Finite Precision Arithmetic,” Acta Numerica 15, 471–542.   
A wealth of practical, Lanczos-related information may be found in:   
J.K. Cullum and R.A. Willoughby (2002). Lanczos Algorithms for Large Symmetric Eigenvalue Computations: Vol. I: Theory, SIAM Publications, Philadelphia, PA.   
J. Brown, M. Chu, D. Ellison, and R. Plemmons (1994). Proceedings of the Cornelius Lanczos International Centenary Conference, SIAM Publications, Philadelphia, PA.   
For a discussion about various reorthogonalization schemes, see:   
C.C. Paige (1970). “Practical Use of the Symmetric Lanczos Process with Reorthogonalization,” BIT 10, 183–195.   
G.H. Golub, R. Underwood, and J.H. Wilkinson (1972). “The Lanczos Algorithm for the Symmetric Ax = λBx Problem,” Report STAN-CS-72-270, Department of Computer Science, Stanford University, Stanford, CA.   
B.N. Parlett and D.S. Scott (1979). “The Lanczos Algorithm with Selective Orthogonalization,” Math. Comput. 33, 217–238.   
H.D. Simon (1984). “Analysis of the Symmetric Lanczos Algorithm with Reorthogonalization Methods,” Lin. Alg. Applic. 61, 101–132.   
Without any reorthogonalization it is necessary either to monitor the loss of orthogonality and quit at the appropriate instant or else to devise a scheme that will identify unconverged eigenvalues and false multiplicities, see:   
W. Kahan and B.N. Parlett (1976). “How Far Should You Go with the Lanczos Process?” in Sparse Matrix Computations, J.R. Bunch and D.J. Rose (eds.), Academic Press, New York, 131–144.   
J. Cullum and R.A. Willoughby (1979). “Lanczos and the Computation in Specified Intervals of the Spectrum of Large, Sparse Real Symmetric Matrices, in Sparse Matrix Proc., I.S. Duff and G.W. Stewart (eds.), SIAM Publications, Philadelphia, PA.   
B.N. Parlett and J.K. Reid (1981). “Tracking the Progress of the Lanczos Algorithm for Large Symmetric Eigenproblems,” IMA J. Num. Anal. 1, 135–155.   
For a restarting framework to be successful, it must exploit the approximate invariant subspace information that has been acquired by the iteration that is about to be shut down, see:   
D. Calvetti, L. Reichel, and D.C. Sorensen (1994). “An Implicitly Restarted Lanczos Method for Large Symmetric Eigenvalue Problems,” ETNA 2, 1–21.   
K. Wu and H. Simon (2000). “Thick-Restart Lanczos Method for Large Symmetric Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 22, 602–616.   
The block Lanczos algorithm is discussed in:   
J. Cullum and W.E. Donath (1974). “A Block Lanczos Algorithm for Computing the q Algebraically Largest Eigenvalues and a Corresponding Eigenspace of Large Sparse Real Symmetric Matrices,” Proceedings of the 1974 IEEE Conference on Decision and Control, Phoenix, AZ, 505–509.   
R. Underwood (1975). “An Iterative Block Lanczos Method for the Solution of Large Sparse Symmetric Eigenvalue Problems,” Report STAN-CS-75-495, Department of Computer Science, Stanford University, Stanford, CA.   
G.H. Golub and R. Underwood (1977). “The Block Lanczos Method for Computing Eigenvalues,” in Mathematical Software III , J. Rice (ed.), Academic Press, New York, pp. 364–377.   
J. Cullum (1978). “The Simultaneous Computation of a Few of the Algebraically Largest and Smallest Eigenvalues of a Large Sparse Symmetric Matrix,” BIT 18, 265–275.   
A. Ruhe (1979). “Implementation Aspects of Band Lanczos Algorithms for Computation of Eigenvalues of Large Sparse Symmetric Matrices,” Math. Comput. 33, 680–687.   
The block Lanczos algorithm generates a symmetric band matrix whose eigenvalues can be computed in any of several ways. One approach is described in:   
H.R. Schwartz (1968). “Tridiagonalization of a Symmetric Band Matrix,” Numer. Math. 12, 231–241.

In some applications it is necessary to obtain estimates of interior eigenvalues. One strategy is to apply Lanczos to the matrix $( A - \mu I ) ^ { - 1 }$ because the extremal eigenvalues of this matrix are eigenvalues close to µ. However, “shift-and-invert” strategies replace the matrix-vector product in the Lanczos iteration with a large sparse linear equation solve, see:

A.K. Cline, G.H. Golub, and G.W. Platzman (1976). “Calculation of Normal Modes of Oceans Using a Lanczos Method,” in Sparse Matrix Computations, J.R. Bunch and D.J. Rose (eds), Academic Press, New York, pp. 409–426.   
T. Ericsson and A. Ruhe (1980). “The Spectral Transformation Lanczos Method for the Numerical Solution of Large Sparse Generalized Symmetric Eigenvalue Problems,” Math. Comput. 35, 1251– 1268.   
R.B. Morgan (1991). “Computing Interior Eigenvalues of Large Matrices,” Lin. Alg. Applic. 154-156, 289–309.   
R.G. Grimes, J.G. Lewis, and H.D. Simon (1994). “A Shifted Block Lanczos Algorithm for Solving Sparse Symmetric Generalized Eigenproblems,” SIAM J. Matrix Anal. Applic. 15, 228–272.

# 10.4 Large Sparse SVD Frameworks

The connections between the SVD problem and the symmetric eigenvalue problem are discussed in §8.6.1. In light of that discussion, it is not surprising that there is a Lanczos process for computing selected singular values and vectors of a large, sparse, rectangular matrix A. The basic idea is to generate a bidiagonal matrix B that is orthogonally equivalent to A. We show how to do this in §5.4 using Householder transformations. However, to avoid large dense submatrices along the way, the Lanczos approach generates the bidiagonal entries entries directly.

# 10.4.1 Golub-Kahan Upper Bidiagonalization

Suppose $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ and recall from §5.4.8 that there exist orthogonal $U \in \mathbb { R } ^ { m \times m }$ and $V \in \mathbb { R } ^ { n \times n }$ so that

$$
U ^ {T} A V = B = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & \dots & \dots & 0 \\ 0 & \alpha_ {2} & \beta_ {2} & \dots & \vdots \\ \vdots & \ddots & \ddots & \ddots & \vdots \\ \vdots & & 0 & \alpha_ {n - 1} & \beta_ {n - 1} \\ 0 & \dots & \dots & 0 & \alpha_ {n} \\ \hline & & 0 \end{array} \right]. \tag {10.4.1}
$$

Since A and B are orthogonally related, they have the same singular values.

Analogously to our derivation of the symmetric Lanczos procedure in §10.1.1, we proceed to outline a sparse-matrix-friendly method for determining the diagonal and superdiagonal of B. The challenge is to bypass the generally full intermediate matrices associated with the Householder bidiagonalization process (Algorithm 5.4.2). We expect to extract good singular value/vector information long before the full bidiagonalization is complete.

The key is to develop useful recipes for the α’s and $\beta \mathrm { { ^ { * } s } }$ from the matrix equations $A V = U B$ and $A ^ { T } U = V { \bar { B } } ^ { T }$ . Given the column partitionings

$$
U = \left[ u _ {1} \mid \dots \mid u _ {m} \right], \quad V = \left[ v _ {1} \mid \dots \mid v _ {n} \right],
$$

we have

$$
A v _ {k} = \alpha_ {k} u _ {k} + \beta_ {k - 1} u _ {k - 1}, \tag {10.4.2}
$$

$$
A ^ {T} u _ {k} = \alpha_ {k} v _ {k} + \beta_ {k} v _ {k + 1} \tag {10.4.3}
$$

for $k = 1 { : } n$ with the convention that $\beta _ { 0 } u _ { 0 } \equiv 0$ and $\beta _ { n } v _ { n + 1 } \equiv 0$ . Define the vectors

$$
r _ {k} = A v _ {k} - \beta_ {k - 1} u _ {k - 1}, \tag {10.4.4}
$$

$$
p _ {k} = A ^ {T} u _ {k} - \alpha_ {k} v _ {k}. \tag {10.4.5}
$$

Using (10.4.2), (10.4.4), and the orthonormality of the u-vectors, we have

$$
\alpha_ {k} = \pm \| r _ {k} \| _ {2},
$$

$$
u _ {k} = r _ {k} / \alpha_ {k}, \quad (\alpha_ {k} \neq 0).
$$

Note that if $\alpha _ { k } \ = \ 0$ , then from (10.4.1) it follows that $A ( : , 1 : k )$ is rank deficient. Similarly we may conclude from (10.4.3) and (10.4.5) that

$$
\beta_ {k} = \pm \| p _ {k} \| _ {2},
$$

$$
v _ {k + 1} = p _ {k} / \beta_ {k}, \quad (\beta_ {k} \neq 0).
$$

If $\beta _ { k } = 0$ , then it follows from the equations $A V = U B$ and $A ^ { T } U = V B ^ { T }$ that

$$
A U (:, 1: k) = V (:, 1: k) B (1: k, 1: k), \tag {10.4.6}
$$

$$
A ^ {T} V (:, 1: k) = U (:, 1: k) B (1: k, 1: k) ^ {T}, \tag {10.4.7}
$$

and thus

$$
A ^ {T} A V (:, 1: k) = V (:, 1: k) B (1: k, 1: k) ^ {T} B (1: k, 1: k).
$$

It follows that $\sigma ( B ( 1 { : } k , 1 { : } k ) ) \subseteq \sigma ( A )$ .

Properly sequenced, the above equations mathematically define the Golub-Kahan process for bidiagonalizing a rectangular matrix.

Algorithm 10.4.1 (Golub-Kahan Bidiagonalization) Given a matrix $A \in \mathbb { R } ^ { m \times n }$ with full column rank and a unit 2-norm vector $v _ { c } \in \mathbb { R } ^ { n }$ , the following algorithm computes the factorizations (10.4.6) and (10.4.7) for some k with $1 \leq k \leq n$ . The first column of V is $v _ { c }$ .

$$
k = 0, p _ {0} = v _ {c}, \beta_ {0} = 1, u _ {0} = 0
$$

while $\beta _ { k } \neq 0$

$$
v _ {k + 1} = p _ {k} / \beta_ {k}
$$

$$
k = k + 1
$$

$$
r _ {k} = A v _ {k} - \beta_ {k - 1} u _ {k - 1}
$$

$$
\alpha_ {k} = \left\| r _ {k} \right\| _ {2}
$$

$$
u _ {k} = r _ {k} / \alpha_ {k}
$$

$$
p _ {k} = A ^ {T} u _ {k} - \alpha_ {k} v _ {k}
$$

$$
\beta_ {k} = \left\| p _ {k} \right\| _ {2}
$$

end

This computation was first described by Golub and Kahan (1965). If $V _ { k } = \left[ v _ { 1 } | \cdot \cdot \cdot | v _ { k } \right]$ , $U _ { k } = [  u _ { 1 } | \cdot \cdot \cdot | u _ { k } ]$ , and

$$
B _ {k} = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & \dots & \dots & 0 \\ 0 & \alpha_ {2} & \beta_ {2} & \dots & \vdots \\ \vdots & \ddots & \ddots & \ddots & 0 \\ \vdots & & 0 & \alpha_ {k - 1} & \beta_ {k - 1} \\ 0 & \dots & 0 & 0 & \alpha_ {k} \end{array} \right], \tag {10.4.8}
$$

then after the kth pass through the loop we have

$$
A V _ {k} = U _ {k} B _ {k}, \tag {10.4.9}
$$

$$
A ^ {T} U _ {k} = V _ {k} B _ {k} ^ {T} + p _ {k} e _ {k} ^ {T}, \tag {10.4.10}
$$

assuming that $\alpha _ { k } > 0$ . It can be shown that

$$
\operatorname{span} \{v _ {1}, \dots , v _ {k} \} = \mathcal {K} (A ^ {T} A, v _ {c}, k), \tag {10.4.11}
$$

$$
\operatorname{span} \{u _ {1}, \dots , u _ {k} \} = \mathcal {K} (A A ^ {T}, A v _ {c}, k). \tag {10.4.12}
$$

Thus, the symmetric Lanczos convergence theory presented in §10.1.5 can be applied. Good approximations to A’s large singular values emerge early, while the small singular values are typically more problematic, especially if there is a cluster near the origin. For further insight, see Luk (1978), Golub, Luk, and Overton (1981), and Bj¨orck (NMLS, §7.6).

# 10.4.2 Ritz Approximations

The Ritz idea can be applied to extract approximate singular values and vectors from the matrices $U _ { k } , V _ { k }$ , and $B _ { k }$ . We simply compute the SVD

$$
F _ {k} ^ {T} B _ {k} G _ {k} = \Gamma = \operatorname{diag} \left(\gamma_ {1}, \dots , \gamma_ {k}\right) \tag {10.4.13}
$$

and form the matrices

$$
Y _ {k} = V _ {k} G _ {k} = \left[ y _ {1} \mid \dots \mid y _ {k} \right],
$$

$$
Z _ {k} = U _ {k} F _ {k} = \left[ z _ {1} \mid \dots \mid z _ {k} \right].
$$

It follows from (10.4.9), (10.4.10), and (10.4.13) that

$$
A Y _ {k} = Z _ {k} \Gamma ,
$$

$$
A ^ {T} Z _ {k} = Y _ {k} \Gamma + p _ {k} e _ {k} ^ {T} F _ {k},
$$

and so for i = 1:k we have

$$
A y _ {i} = \gamma_ {i} z _ {i}, \tag {10.4.14}
$$

$$
A ^ {T} z _ {i} = \gamma_ {i} y _ {i} + [ F _ {k} ] _ {k i} \cdot p _ {k}. \tag {10.4.15}
$$

It follows that $A ^ { T } A y _ { i } = \gamma _ { i } ^ { 2 } z _ { i } + [ F _ { k } ] _ { k i } \cdot p _ { k }$ and thus, $\{ \gamma _ { i } , y _ { i } \}$ is a Ritz pair for $A ^ { T } A$ with respect to ran(Vk).

# 10.4.3 The Tridiagonal-Bidiagonal Connection

In §8.6.1 we showed that there is a connection between the SVD of a matrix $A \in \mathbb { R } ^ { m \times n }$ and the Schur decomposition of the symmetric matrix

$$
C = \left[ \begin{array}{c c} 0 & A \\ A ^ {T} & 0 \end{array} \right]. \tag {10.4.16}
$$

In particular, if $\sigma$ is a singular value of A, then both $\sigma$ and $- \sigma$ are eigenvalues of C and the corresponding singular vectors “makeup” the corresponding eigenvectors.

Likewise, a given bidiagonalization of A can be related to a tridiagonalization of C. Assume that $m \geq n$ and that

$$
[ U _ {1} \mid U _ {2} ] ^ {T} A V = \left[ \begin{array}{c} \tilde {B} \\ 0 \end{array} \right], \qquad \tilde {B} \in \mathbb {R} ^ {n \times n},
$$

is a bidiagonalization of A with $U _ { 1 } \in \mathbb { R } ^ { m \times n } , \ U _ { 2 } \in \mathbb { R } ^ { m \times ( m - n ) }$ , and $V \in \mathbb { R } ^ { n \times n }$ . Note that

$$
Q = \left[ \begin{array}{l l} U & 0 \\ 0 & V \end{array} \right]
$$

is orthogonal and

$$
\tilde {T} = Q ^ {T} C Q = \left[ \begin{array}{c c} 0 & \tilde {B} \\ \tilde {B} ^ {T} & 0 \end{array} \right].
$$

This matrix can be symmetrically permuted into tridiagonal form. For example, in the 4-by-3 case, if $P = I _ { 7 } ( : , [ 5 1 6 2 7 3 4 ] )$ , then the reordering $\tilde { T }  P \tilde { T } P ^ { T }$ has the form

$$
\left[\begin{array}{c c c c c c c}0&0&0&0&\alpha_ {1}&\beta_ {1}&0\\0&0&0&0&0&\alpha_ {2}&\beta_ {2}\\0&0&0&0&0&0&\alpha_ {3}\\0&0&0&0&0&0&0\\\hline \alpha_ {1}&0&0&0&0&0&0\\\beta_ {1}&\alpha_ {2}&0&0&0&0&0\\0&\beta_ {2}&\alpha_ {3}&0&0&0&0\end{array}\right] \quad \rightarrow \quad \left[\begin{array}{c c c c c c c}0&\alpha_ {1}&0&0&0&0&0\\\alpha_ {1}&0&\beta_ {1}&0&0&0&0\\0&\beta_ {1}&0&\alpha_ {2}&0&0&0\\0&0&\alpha_ {2}&0&\beta_ {2}&0&0\\0&0&0&\beta_ {2}&0&\alpha_ {3}&0\\0&0&0&0&\alpha_ {3}&0&0\\0&0&0&0&0&0&0\end{array}\right].
$$

This points to an interesting connection between Golub-Kahan bidiagonalization (Algorithm 10.4.1) and Lanczos tridiagonalization (Algorithm 10.1.1). Suppose we apply Algorithm 10.4.1 to $A \in \mathbb { R } ^ { m \times n }$ with starting vector $v _ { c } .$ . Assume that the procedure runs for k steps and produces the bidiagonal matrix $B _ { k }$ displayed in (10.4.8). If we apply Algorithm 10.1.1 to the matrix C defined by (10.4.16) with a starting vector

$$
q _ {1} = \left[ \begin{array}{c} 0 \\ v _ {c} \end{array} \right] \in \mathbb {R} ^ {m + n} \tag {10.4.17}
$$

then after 2k steps the resulting tridiagonal matrix $T _ { 2 k }$ has a zero diagonal and a subdiagonal specified by $\big [ \alpha _ { 1 } , \beta _ { 1 } , \alpha _ { 2 } , \beta _ { 2 } , \cdot \cdot \cdot \alpha _ { k - 1 } , \beta _ { k - 1 } , \alpha _ { k } \big ]$ .

# 10.4.4 Paige-Saunders Lower Bidiagonalization

In §11.4.2 we show how the Golub-Kahan bidiagonalization can be used to solve sparse linear systems and least squares problems. It turns out that in this context, lower bidiagonalization is more useful:

$$
U ^ {T} A V = B = \left[ \begin{array}{c c c c c} \alpha_ {1} & 0 & \dots & \dots & 0 \\ \beta_ {2} & \alpha_ {2} & 0 & \dots & \vdots \\ \vdots & \beta_ {3} & \ddots & \ddots & \vdots \\ \vdots & & \ddots & \alpha_ {n - 1} & 0 \\ 0 & \dots & \dots & \beta_ {n - 1} & \alpha_ {n} \\ 0 & \dots & \dots & 0 & \beta_ {n} \\ \hline & & 0 \end{array} \right]. \tag {10.4.18}
$$

Proceeding as in the derivation of the Golub-Kahan bidiagonalization, we compare columns in the equations $A ^ { T } U = V B ^ { T }$ and $A V = U B$ . If $U = \ [ \ u _ { 1 } \ | \cdot \cdot \cdot | \ u _ { m } \ ]$ and $V = [  v _ { 1 } | \cdot \cdot \cdot | v _ { n } ]$ are column partitionings and we define $\beta _ { 1 } v _ { 0 } \equiv 0$ and $\alpha _ { n + 1 } v _ { n + 1 } \equiv 0$ , then for $k = 1 { : } n$ we have $A ^ { T } u _ { k } = \beta _ { k } v _ { k - 1 } + \alpha _ { k } v _ { k }$ and $A v _ { k } = \alpha _ { k } u _ { k } + \beta _ { k + 1 } u _ { k + 1 }$ . Leaving the rest of the derivation to the exercises, we obtain the following.

Algorithm 10.4.2 (Paige-Saunders Bidiagonalization) Given a matrix $A \in \mathbb { R } ^ { m \times n }$ with the property that $A ( 1 { : } n , 1 { : } n )$ is nonsingular and a unit 2-norm vector $u _ { c } \in \mathbb { R } ^ { n }$ , the following algorithm computes the factorization $A V ( : , 1 : k ) = U ( : , 1 : k + 1 ) B ( 1 : k + 1 , 1 : k )$ where $U , V$ , and B are given by (10.4.18). The first column of U is $u _ { c }$ and the integer k satisfies $1 \leq k \leq n$ .

$$
k = 1, p _ {0} = u _ {c}, \beta_ {1} = 1, v _ {0} = 0
$$

while $\beta _ { k } > 0$

$$
u _ {k} = p _ {k - 1} / \beta_ {k}
$$

$$
r _ {k} = A ^ {T} u _ {k} - \beta_ {k} v _ {k - 1}
$$

$$
\alpha_ {k} = \left\| r _ {k} \right\| _ {2}
$$

$$
v _ {k} = r _ {k} / \alpha_ {k}
$$

$$
p _ {k} = A v _ {k} - \alpha_ {k} u _ {k}
$$

$$
\beta_ {k + 1} = \left\| p _ {k} \right\| _ {2}
$$

$$
k = k + 1
$$

end

It can be shown that after k passes through the loop we have

$$
A V (:, 1: k) = U (:, 1: k) B (1: k, 1: k) + p _ {k} e _ {k} ^ {T} \tag {10.4.19}
$$

where $e _ { k } = I _ { k } ( : , k )$ . See Paige and Saunders (1982) for more details. Their bidiagonalization is equivalent to Golub-Kahan bidiagonalization applied to $[ b | A ]$ .


---

<!-- golub_600_649 -->

The need to extract information from unimaginably large datasets has prompted the development of matrix methods that involve randomization. The idea is to develop matrix approximations that are very fast to compute because they rely on limited, random samplings of the given matrix. To give a snapshot of this increasingly important paradigm for large-scale matrix computations, we consider the problem of computing a rank-k approximation to a given matrix $A \in \mathbb { R } ^ { m \times n }$ . For clarity we assume that $k \leq \mathsf { r a n k } ( A )$ . Recall that if $A = \tilde { Z } \tilde { \Sigma } \tilde { Y } ^ { T }$ is the SVD of A, then

$$
\tilde {A} _ {k} = \tilde {Z} _ {1} \tilde {\Sigma} _ {1} \tilde {Y} _ {1} ^ {T} = \tilde {Z} _ {1} \tilde {Z} _ {1} ^ {T} A \tag {10.4.20}
$$

where $\tilde { Z } _ { 1 } = \tilde { Z } ( : , 1 : k ) , \tilde { \Sigma } _ { 1 } = \tilde { \Sigma } ( 1 : k , 1 : k )$ , and $\tilde { Y } _ { 1 } = \tilde { Y } ( : , 1 { : } k )$ , is the closest rank-k matrix to A as measured in either the 2-norm or Frobenius norm. We assume that A is so large that the Krylov methods just discussed are impractical.

Drineas, Kannan, and Mahoney (2006c) propose a method that approximates the intractable $\tilde { A } _ { k }$ with a rank-k matrix of the form

$$
A _ {k} = C U R, \quad C \in \mathbb {R} ^ {m \times c}, U \in \mathbb {R} ^ {c \times r}, R \in \mathbb {R} ^ {r \times n}, k \leq c, k \leq r \tag {10.4.21}
$$

where the matrices C and R are comprised of randomly chosen values taken from A. The integers c and r are parameters of the method. Discussion of the CUR decomposition (10.4.21) nicely illustrates the notion of random sampling in the matrix context and the idea of a probabilistic error bound.

The first step in the CUR framework is to determine C. Each column of this matrix is a scaled, randomly-selected column of A:

Determine column probabilities $q _ { j } = \parallel A ( : , j ) \parallel _ { 2 } / \parallel A \parallel _ { _ { F } } ^ { 2 } , j = 1 : n .$

for $t = 1 { : } c$

Randomly pick col $( t ) \in \{ 1 , 2 , \ldots , n \}$ with $q _ { \alpha }$ the probability that $c o l ( t ) = \alpha$

$$
C (:, t) = A (:, c o l (t)) / \sqrt {c q _ {c o l (t)}}
$$

end

It follows that $C = A ( : , c o l ) D _ { c }$ where $D _ { c } \in \mathbb { R } ^ { c \times c }$ is a diagonal scaling matrix.

The matrix R is similarly constructed. Each row of this matrix is a scaled, randomly-selected row of A:

Determine row probabilities $p _ { i } = \parallel A ( i , : ) \parallel _ { 2 } / \parallel A \parallel _ { F } ^ { 2 } , i = 1 { : } m .$

for $t = 1 { : } r$

Randomly pick row $( t ) \in \{ 1 , 2 , \ldots , m \}$ with $p _ { \alpha }$ the probability that $r o w ( t ) = \alpha$

$$
R (t,:) = A (r o w (t),:) / \sqrt {r p _ {r o w (t)}}
$$

end

The matrix R has the form ${ \cal R } ~ = ~ D _ { \scriptscriptstyle R } A ( r o w , \colon )$ where $D _ { R } \in \mathbb { R } ^ { r \times r }$ is a diagonal scaling matrix.

The next step is to choose a rank-k matrix U so that $A _ { k } = C U R$ is close to the best rank-k approximation $\tilde { A } _ { k }$ . In the CUR framework, this requires the SVD

$$
C = Z \Sigma Y ^ {T} = Z _ {1} \Sigma_ {1} Y _ {1} ^ {T} + Z _ {2} \Sigma_ {2} Y _ {2}
$$

where $Z _ { 1 } = Z ( : , 1 : k ) , \ : \Sigma _ { 1 } = \Sigma ( 1 : k , 1 : k )$ , and $Y _ { 1 } = Y ( : , 1 { : } k )$ . The matrix U is then given by

$$
U = \Phi \Psi^ {T}, \qquad \Phi = Y _ {1} \Sigma_ {1} ^ {- 2} Y _ {1} ^ {T}, \Psi = D _ {R} C (r o w,:).
$$

With these definitions, simple manipulations confirm that

$$
C \Phi = Z _ {1} \Sigma_ {1} ^ {- 1} Y _ {1} ^ {T}, \tag {10.4.22}
$$

$$
\Psi^ {T} R = \left(D _ {R} \left(Z _ {1} (r o w,:) \Sigma_ {1} Y _ {1} ^ {T} + Z _ {2} (r o w,:) \Sigma_ {2} Y _ {2} ^ {T})\right) ^ {T} D _ {R} A (r o w,:), \right. \tag {10.4.23}
$$

and

$$
C U R = (C \Phi) (\Psi R) = Z _ {1} \left(D _ {R} Z _ {1} (r o w,:)\right) \left(D _ {R} A (r o w,:)\right). \tag {10.4.24}
$$

An analysis that critically depends on the selection probabilities $\left\{ q _ { i } \right\}$ and $\{ p _ { i } \}$ shows that ran $( Z _ { 1 } ) \approx \mathsf { r a n } ( \tilde { Z } _ { 1 } )$ and $\left( D _ { R } Z _ { 1 } ( r o w , : ) \right) ^ { T } \left( D _ { R } A ( r o w , : ) \right) \approx Z _ { 1 } ^ { T } A$ . Upon comparison with (10.4.20) we see that $C U R \approx Z _ { 1 } Z _ { 1 } ^ { T } A \approx \tilde { Z } _ { 1 } \tilde { Z } _ { 1 } ^ { T } A = \tilde { A } _ { k }$ . Moreover, given $\epsilon > 0$ , $\delta > 0$ , and k, it is possible to choose the parameters r and c so that the inequality

$$
\left\| A - C U R \right\| _ {F} \leq \left\| A - \tilde {A} _ {k} \right\| _ {F} + \epsilon \left\| A \right\| _ {F}
$$

holds with probability 1 − δ. Lower bounds for r and c that depend inversely on 
 and $\delta$ are given by Drineas, Kannan, and Mahoney (2006c).

# Problems

P10.4.1 Verify Equations (10.4.6), (10.4.7), (10.4.9), and (10.4.10).

P10.4.2 Corresponding to (10.3.1), develop an implementation of Algorithm 10.4.1 that involves a minimum number of vector workspaces.

P10.4.3 Show that if rank(A) = n, then the bidiagonal matrix B in (10.4.18) cannot have a zero on its diagonal.

P10.4.4 Prove (10.4.19). What can you say about $U ( : , 1 : k )$ and $V ( : , 1 : k )$ if $\beta _ { k + 1 } = 0$ in Algorithm 10.4.2?

P10.4.5 Analogous to (10.4.11)-(10.4.12), show that for Algorithm 10.4.2 we have

$$
\operatorname{span} \left\{v _ {1}, \dots , v _ {k} \right\} = \mathcal {K} \left(A ^ {T} A, A ^ {T} u _ {c}, k\right), \quad \operatorname{span} \left\{u _ {1}, \dots , u _ {k} \right\} = \mathcal {K} \left(A A ^ {T}, u _ {c}, k\right).
$$

P10.4.6 Suppose C and $q _ { 1 }$ are defined by (10.4.16) and (10.4.17) respectively. (a) Show that

$$
\mathcal {K} (C, q _ {1}, 2 k) = \text { span } \left\{\left[ \begin{array}{c} 0 \\ v _ {c} \end{array} \right], \left[ \begin{array}{c} A v _ {c} \\ 0 \end{array} \right], \left[ \begin{array}{c} 0 \\ A ^ {T} A v _ {c} \end{array} \right], \ldots , \left[ \begin{array}{c} 0 \\ (A ^ {T} A) ^ {k - 1} v _ {c} \end{array} \right], \left[ \begin{array}{c} A (A ^ {T} A) ^ {k - 1} v _ {c} \\ 0 \end{array} \right] \right\}.
$$

(b) Rigorously prove the claim made in §10.4.3 about the subdiagonal of $T _ { 2 k }$ . (c) State and prove analogous results when the Paige-Saunders bidiagonalization is used.

P10.4.7 Verify Equations 10.4.22–10.4.24.

# Notes and References for 10.4

For a more comprehensive treatment of Golub-Kahan bidiagonalization, see Bj¨orck (NMLS, §7.6). The relevance of the Lanczos process to the bidiagonalization of a rectangular matrix was first presented in:

G.H. Golub and W. Kahan (1965). “Calculating the Singular Values and Pseudo-Inverse of a Matrix,” SIAM J. Numer. Anal. Ser. B, 2, 205–224.

The idea of using Golub-Kahan bidiagonalization to solve large sparse linear systems and least squares problems started with the paper:

C.C. Paige (1974). “Bidiagonalization of Matrices and Solution of Linear Equations,” SIAM J. Numer. Anal. 11, 197–209.   
We shall have more to say about this in the next chapter. It is in anticipation of that discussion that we presented the lower bidiagonal scheme, see:   
C.C. Paige and M.A. Saunders (1982). “LSQR, An Algorithm for Sparse Linear Equations and Sparse Least Squares,” ACM Trans. Math. Softw. 8, 43–71.   
For practical implementation issues, see:   
G.H. Golub, F.T. Luk, and M.L. Overton (1981). “A Block Lanczos Method for Computing the Singular Values and Corresponding Singular Vectors of a Matrix,” ACM Trans. Math. Softw. 7, 149–169.   
J. Cullum, R.A. Willoughby, and M. Lake (1983). “A Lanczos Algorithm for Computing Singular Values and Vectors of Large Matrices,” SIAM J. Sci. Stat. Comput. 4, 197–215.   
M. Berry (1992). “Large-Scale Sparse Singular Value Computations,” International J. Supercomputing Appl. 6, 13–49.   
M. Berry and R.L. Auerbach (1993). “A Block Lanczos SVD Method with Adaptive Reorthogonalization,” in Proceedings of the Cornelius Lanczos International Centenary Conference, Raleigh, NC, SIAM Publications, Philadelphia, PA.   
Z. Jia and D. Niu (2003). “An Implicitly Restarted Refined Bidiagonalization Lanczos Method for Computing a Partial Singular Value Decomposition,” SIAM J. Matrix Anal. Applic. 25, 246–265.   
Interesting applications of the Lanczos bidiagonalization include:   
D.P. OLeary and J.A. Simmons (1981). “A Bidiagonalization-Regularization Procedure for Large Scale Discretizations of Ill-Posed Problems,” SIAM J. Sci. Stat. Comput. 2, 474–489.   
D. Calvetti, G.H. Golub, and L. Reichel (1999). “Estimation of the L-curve via Lanczos Bidiagonalization,” BIT 39, 603–619.   
H.D. Simon and H. Zha (2000). “Low-Rank Matrix Approximation Using the Lanczos Bidiagonalization Process with Applications,” SIAM J. Sci. Comput. 21, 2257–2274.   
Our sketch of the CUR decomposition framework is based on:   
P. Drineas, R. Kannan, and M.W. Mahoney (2006). “Fast Monte Carlo Algorithms for Matrices III: Computing an Efficient Approximate Decomposition of a Matrix,” SIAM J. Comput. 36, 184–206.   
Additional references concerned with randomization in matrix computations include:   
P. Drineas, R. Kannan, and M. W. Mahoney (2006). “Fast Monte Carlo Algorithms for Matrices I: Approximating Matrix Multiplication,” SIAM J. Comput. 36, 132–157.   
P. Drineas, R. Kannan, and M.W. Mahoney (2006). “Fast Monte Carlo Algorithms for Matrices II: Computing Low-Rank Approximations to a Matrix,” SIAM J. Comput. 36, 158–183.   
M.W. Mahoney, M. Maggioni, and P. Drineas (2008). “Tensor-CUR Decompositions For Tensor-Based Data,” SIAM J. Mat. Anal. Applic. 30, 957–987.   
P. Drineas, M.W. Mahoney, and S. Muthukrishnan (2008). “Relative-Error CUR Matrix Decompositions,” SIAM J. Mat. Anal. Applic. 30, 844–881.   
E. Liberty, F. Woolfe, P.-G. Martinsson, V. Rokhlin, and M.Tygert (2008). “Randomized Algorithms for the Low-Rank Approximation of Matrices,” Proc. Natl. Acad. Sci. 104, 20167–20172.   
V. Rokhlin and Mark Tygert (2008). “A Fast Randomized Algorithm for Overdetermined Linear Least-Squares Regression,” Proc. Natl. Acad. Sci. 105, 13212–13217.   
M.W. Mahoney and P. Drineas (2009). “CUR Matrix Decompositions for Improved Data Analysis,” Proc. Natl. Acad. Sci. 106, 697–702.   
D. Achlioptas and F. McSherry (2007). “Fast Computation of Low-Rank Matrix Approximations,” JACM 54(2), Article No. 9.   
V. Rokhlin, A. Szlam, and M. Tygert (2010). “A Randomized Algorithm for Principal Component Analysis,” SIAM J. Mat. Anal. Applic. 31, 1100–1124   
M.W. Mahoney (2011). “Randomized Algorithms for Matrices and Data,” Foundations and Trends in Machine Learning 3, 123–224.   
N. Halko, P.G. Martinsson, and J.A. Tropp (2011). “Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions,” SIAM Review 53, 217–288   
For another perspective on the increasing important role of randomness in matrix computations, see:

A. Edelman and N. Raj Rao (2005). “Random Matrix Theory,” Acta Numerica 14, 233–297

# 10.5 Krylov Methods for Unsymmetric Problems

If A is not symmetric, then the orthogonal tridiagonalization $Q ^ { T } A Q = T$ does not exist in general. There are two ways to proceed. The Arnoldi approach involves the column-by-column generation of an orthogonal Q such that $Q ^ { T } A \bar { Q } = H$ is the Hessenberg reduction of $\ S 7 . 4$ . The unsymmetric Lanczos approach computes the columns of matrices Q and P so that $P ^ { T } A Q = T$ is tridiagonal and $P ^ { T } Q = I$ . Methods based on these ideas that are suitable for large, sparse, unsymmetric eigenvalue problems are discussed in this section.

# 10.5.1 The Basic Arnoldi Process

One way to extend the Lanczos process to unsymmetric matrices is due to Arnoldi (1951) and revolves around the Hessenberg reduction $Q ^ { T } A Q = H$ . In particular, if $Q = [  q _ { 1 } | \cdots | q _ { n } ]$ and we compare columns in $A Q = Q H$ , then

$$
A q _ {k} = \sum_ {i = 1} ^ {k + 1} h _ {i k} q _ {i}, \quad 1 \leq k \leq n - 1.
$$

Isolating the last term in the summation gives

$$
h _ {k + 1, k} q _ {k + 1} = A q _ {k} - \sum_ {i = 1} ^ {k} h _ {i k} q _ {i} \equiv r _ {k}
$$

where $h _ { i k } = q _ { i } ^ { T } A q _ { k }$ for $i = 1 { : } k$ . It follows that if $r _ { k } \neq 0$ , then $q _ { k + 1 }$ is specified by

$$
q _ {k + 1} = r _ {k} / h _ {k + 1, k}
$$

where $h _ { k + 1 , k } = \parallel r _ { k } \parallel _ { 2 }$ . These equations define the Arnoldi process and in strict analogy to the symmetric Lanczos process (Algorithm 10.1.1) we obtain the following.

Algorithm 10.5.1 (Arnoldi Process) If $A \in \mathbb { R } ^ { n \times n }$ and $q _ { 1 } \in \mathbb { R } ^ { n }$ has unit 2-norm, then the following algorithm computes a matrix $Q _ { t } = [ q _ { 1 } , \dots , q _ { t } ] \in \mathbb { R } ^ { n \times t }$ with orthonormal columns and an upper Hessenberg matrix $H _ { t } = ( h _ { i j } ) \in \mathbb { R } ^ { t \times t }$ with the property that $A Q _ { t } = Q _ { t } H _ { t }$ . The integer t satisfies $1 \leq t \leq n$ .

$$
k = 0, r _ {0} = q _ {1}, h _ {1 0} = 1
$$

while $( h _ { k + 1 , k } \neq 0 )$

$$
q _ {k + 1} = r _ {k} / h _ {k + 1, k}
$$

$$
k = k + 1
$$

$$
r _ {k} = A q _ {k}
$$

for i = 1:k

$$
h _ {i k} = q _ {i} ^ {T} r _ {k}
$$

$$
r _ {k} = r _ {k} - h _ {i k} q _ {i}
$$

end

$$
h _ {k + 1, k} = \left\| r _ {k} \right\| _ {2}
$$

end

$$
t = k
$$

The $q _ { k }$ are called Arnoldi vectors and they define an orthonormal basis for the Krylov subspace $\kappa ( A , q _ { 1 } , k )$ :

$$
\operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\} = \operatorname{span} \left\{q _ {1}, A q _ {1}, \dots , A ^ {k - 1} q _ {1} \right\}. \tag {10.5.1}
$$

The situation after k steps is summarized by the equation

$$
A Q _ {k} = Q _ {k} H _ {k} + r _ {k} e _ {k} ^ {T} \tag {10.5.2}
$$

where $Q _ { k } = \left[ q _ { 1 } | \cdots | q _ { k } \right] , e _ { k } = I _ { k } ( : , k )$ , and

$$
H _ {k} = \left[ \begin{array}{c c c c c} h _ {1 1} & h _ {1 2} & \dots & \dots & h _ {1 k} \\ h _ {2 1} & h _ {2 2} & \dots & \dots & h _ {2 k} \\ 0 & h _ {3 2} & \ddots & & \vdots \\ \vdots & & \ddots & \ddots & \vdots \\ 0 & \dots & \dots & h _ {k, k - 1} & h _ {k k} \end{array} \right].
$$

Any decomposition of the form (10.5.2) is a k-step Arnoldi decomposition if $Q _ { k } \in \mathbb { R } ^ { n \times k }$ has orthonormal columns, $H _ { k } \in \mathbb { R } ^ { k \times k }$ is upper Hessenberg, and $Q _ { k } ^ { T } r _ { k } = 0 .$ .

If $\boldsymbol { y } \in \mathbb { R } ^ { k }$ is a unit 2-norm eigenvector for $H _ { k }$ and $H _ { k } y = \lambda y$ , then from (10.5.2)

$$
(A - \lambda I) x = (e _ {k} ^ {T} y) r _ {k}
$$

where $x = Q _ { k } y$ . Since $r _ { k } \in \mathcal { K } ( A , q _ { 1 } , k ) ^ { \perp }$ , it follows that $( \lambda , x )$ is a Ritz pair for A with respect to $\kappa ( A , q _ { 1 } , k )$ . Note that if $v = ( e _ { k } ^ { T } y ) r _ { k }$ , then

$$
(A + E) x = \lambda x
$$

where $E = - v x ^ { T }$ with $\left\| \mathbf { \nabla } E \right\| _ { 2 } = \left| y _ { k } \right| \left\| r _ { k } \right\| _ { 2 }$ . Recall that in the unsymmetric case, computing an eigenvalue of a nearby matrix does not mean that it is close to an exact eigenvalue.

Some numerical properties of the Arnoldi iteration are discussed by Wilkinson (AEP, p. 382). The history of practical Arnoldi-based eigensolvers begins with Saad (1980). Two features of the method distinguish it from the symmetric Lanczos process:

• Arnoldi vectors $q _ { 1 } , \ldots , q _ { k }$ must all be referenced in step k and the computation of $q _ { k + 1 }$ involves $O ( k n )$ flops excluding the matrix-vector product $A q _ { k }$ Thus, there is a steep penalty associated with the generation of long Arnoldi sequences.   
• Extremal eigenvalue information is not as forthcoming as in the symmetric case. There is no unsymmetric Kaniel-Paige-Saad convergence theory.

These realities suggest a framework in which we use the Arnoldi iteration idea with repeated, carefully chosen restarts and a controlled iteration maximum. We described such a framework in conjunction with the block Lanczos procedure in §10.3.7.

# 10.5.2 Arnoldi with Restarting

Consider running Arnoldi for m steps and then restarting the iteration with a new initial vector $q _ { + }$ chosen from the span of the Arnoldi vectors $q _ { 1 } , \ldots , q _ { m }$ . Because of the Krylov connection (10.5.1), $q _ { + }$ has the form

$$
q _ {+} = p (A) q _ {1}
$$

for some polynomial of degree $m - 1$ . It is instructive to examine the action of $p ( A )$ in terms of A’s eigenvalues and eigenvectors. Assume for clarity that $A \in \mathbb { R } ^ { n \times n }$ is diagonalizable and that $A z _ { i } = \lambda _ { i } z _ { i }$ for $i = 1 { : } n$ . If $q _ { 1 }$ has the eigenvector expansion

$$
q _ {1} = a _ {1} z _ {1} + \dots + a _ {n} z _ {n},
$$

then $q _ { + }$ is a scalar multiple of

$$
z = a _ {1} p (\lambda_ {1}) z _ {1} + \dots + a _ {n} p (\lambda_ {n}) z _ {n}.
$$

Note that if $p ( \lambda _ { \alpha } ) \gg p ( \lambda _ { \beta } )$ , then relatively speaking, $q _ { + }$ is much richer in the direction of $z _ { \alpha }$ than in the direction of $z _ { \beta }$ . More generally, by carefully choosing $p ( \lambda )$ we can design $q _ { + }$ so that its component in certain eigenvector directions is emphasized while its component in other eigenvector directions is deemphasized. For example, if

$$
p (\lambda) = c \cdot (\lambda - \mu_ {1}) (\lambda - \mu_ {2}) \dots (\lambda - \mu_ {p}) \tag {10.5.3}
$$

where $c$ is a constant, then $q _ { + }$ is a unit vector in the direction of

$$
z = c \cdot \sum_ {k = 1} ^ {n} a _ {k} \left(\prod_ {i = 1} ^ {p} \left(\lambda_ {k} - \mu_ {i}\right)\right) z _ {k}.
$$

It follows that $z _ { \beta }$ is deemphasized relative to $z _ { \alpha }$ if $\lambda _ { \beta }$ is near to one of the “filter values” $\mu _ { 1 } , \ldots , \mu _ { p }$ and $\lambda _ { \alpha }$ is not. Thus, the act of picking a good restart vector $q _ { + }$ from $\mathcal { K } ( A , q _ { 1 } , m )$ is the act of picking a filter polynomial that tunes out unwanted portions of the spectrum. Various heuristics for doing this have been developed based on computed Ritz vectors. See Saad (1980, 1984, 1992).

# 10.5.3 Implicit Restarting

We describe an Arnoldi restarting procedure due to Sorensen (1992) that implicitly determines the filter polynomial (10.5.3) using the QR iteration with shifts. (See §7.5.2.) Suppose $H _ { c } \in \mathbb { R } ^ { m \times m }$ is upper Hessenberg, $\mu _ { 1 } , \ldots , \mu _ { p }$ are scalars, and the matrix $H _ { + }$ is obtained via the shifted QR iteration:

$$
H ^ {(0)} = H _ {c}
$$

for $i = 0 { : } p$

$$
H ^ {(i - 1)} - \mu_ {i} I = V _ {i} R _ {i} \quad \text {(Givens QR)} \tag {10.5.4}
$$

$$
H ^ {(i)} = R _ {i} V _ {i} + \mu_ {i} I
$$

end

$$
H _ {+} = H ^ {(p)}
$$

Recall from §7.4.2 that each $H ^ { ( i ) }$ is upper Hessenberg. Moreover, if

$$
V = V _ {1} \dots V _ {p}, \tag {10.5.5}
$$

then

$$
H _ {+} = V ^ {T} H _ {c} V. \tag {10.5.6}
$$

The following result shows that the filter polynomial (10.5.3) has a relationship to (10.5.4).

Theorem 10.5.1. $I f V = V _ { 1 } \cdots V _ { p }$ and $R = R _ { p } \cdot \cdot \cdot R _ { 1 }$ are defined by $( 1 0 . 5 . 4 )$ , then

$$
V R = \left(H _ {c} - \mu_ {1} I\right) \dots \left(H _ {c} - \mu_ {p} I\right). \tag {10.5.7}
$$

Proof. We use induction, noting that the theorem is obviously true if $p = 1$ . If $\tilde { V } = V _ { 1 } \cdot \cdot \cdot V _ { p - 1 }$ and $\tilde { R } = R _ { p - 1 } \cdot \cdot \cdot R _ { 1 }$ , then

$$
\begin{array}{l} V R = \tilde {V} (V _ {p} R _ {p}) \tilde {R} = \tilde {V} (H ^ {(p - 1)} - \mu_ {p} I) \tilde {R} = \tilde {V} (\tilde {V} ^ {T} H _ {c} \tilde {V} - \mu_ {p} I) \tilde {R} \\ = \left(H _ {c} - \mu_ {p} I\right) \tilde {V} \tilde {R} = \left(H _ {c} - \mu_ {p} I\right) \left(H _ {c} - \mu_ {1} I\right) \dots \left(H _ {c} - \mu_ {p - 1} I\right), \\ \end{array}
$$

where we used the fact that $H ^ { ( p - 1 ) } = \tilde { V } ^ { T } H _ { c } \tilde { V }$ .

Note that the matrix R in (10.5.7) is upper triangular and so it follows that

$$
V (:, 1) = p (H _ {c}) e _ {1}
$$

where $p ( \lambda )$ is the filter polynomial (10.5.3) with $c = 1 / R ( 1 , 1 )$ .

Now suppose that we have performed m steps of the Arnoldi iteration with starting vector $q _ { 1 }$ . The Arnoldi factorization (10.5.2) says that we have an upper Hessenberg matrix $H _ { c } \in \mathbb { R } ^ { m \times m }$ and a matrix $Q _ { c } \in \mathbb { R } ^ { n \times m }$ with orthonormal columns such that

$$
A Q _ {c} = Q _ {c} H _ {c} + r _ {c} e _ {m} ^ {T}. \tag {10.5.8}
$$

Note that $Q _ { c } ( : , 1 ) = q _ { 1 }$ and $r _ { c } \in \mathbb { R } ^ { n }$ has the property that $Q _ { c } ^ { T } r _ { c } = 0$ . If we apply (10.5.4) to $H _ { c }$ , then by using (10.5.5) and (10.5.6) the preceding Arnoldi factorization transforms to

$$
A Q _ {+} = Q _ {+} H _ {+} + r _ {c} e _ {m} ^ {T} V \tag {10.5.9}
$$

where

$$
Q _ {+} = Q _ {c} V.
$$

If $q _ { + }$ is the first column of this matrix, then

$$
q _ {+} = Q _ {+} (:, 1) = Q _ {c} V (:, 1) = c \cdot Q _ {c} \left(H _ {c} - \mu_ {1} I\right) \dots \left(H _ {c} - \mu_ {p} I\right) e _ {1}.
$$

Equation (10.5.8) implies that

$$
(A - \mu I) Q _ {c} e _ {1} = Q _ {c} (H _ {c} - \mu I) e _ {1}
$$

for any $\mu \in \mathbb { R }$ and so

$$
q _ {+} = c (A - \mu_ {1} I) \dots (A - \mu_ {p} I) Q _ {c} e _ {1} = p (A) q _ {1}.
$$

This suggests the following framework for repeated restarting:

# Repeat:

With starting vector $q _ { 1 }$ , perform m steps of the Arnoldi iteration obtaining $Q _ { c } \in \mathbb { R } ^ { n \times m }$ and $H _ { c } \in \mathbb { R } ^ { m \times m }$ .

Determine filter values µ1, . . . , µp . $\mu _ { 1 } , \ldots , \mu _ { p }$ (10.5.10)

Perform $p$ steps of the shifted QR iteration (10.5.4) obtaining the Hessenberg matrix $H _ { + }$ and the orthogonal matrix $V$ .

Replace $q _ { 1 }$ with the first column of $Q _ { c } V$ .

However, we can do better than this. The orthogonal matrices $V _ { 1 } , \ldots , V _ { p }$ that arise in (10.5.4) are each upper Hessenberg. (This is easily deduced from the structure of the Givens rotations in Algorithm 5.2.5.) Thus, V has lower bandwidth $p$ and so $V ( m , 1 : m - p - 1 ) = 0$ . It follows from (10.5.9) that if $j = m - p ,$ then

$$
A Q _ {+} (:, 1: j) = Q _ {+} (:, 1: j) H _ {+} (1: j, 1: j) + v _ {m j} r _ {c} e _ {j}
$$

is a j-step Arnoldi decomposition. In other words, we are all set to perform step $j + 1$ of the Arnoldi iteration with starting vector $q _ { + }$ . There is no need to launch the restart from step 1. This leads us to the following modification of (10.5.10):

With starting vector $q _ { 1 }$ , perform m steps of the Arnoldi iteration obtaining

$$
Q _ {c} \in \mathbb {R} ^ {n \times m}, H _ {c} \in \mathbb {R} ^ {m \times m}, \text {   and   } r _ {c} \in \mathbb {R} ^ {n} \text {   so   } A Q _ {c} = Q _ {c} H _ {c} + r _ {c} e _ {m} ^ {T}.
$$

# Repeat:

Determine filter values $\mu _ { 1 } , \ldots , \mu _ { p }$ .

Perform $p$ steps of the shifted QR iteration (10.5.4) applied to $H _ { c }$ obtaining $H _ { + } \in \mathbb { R } ^ { m \times m }$ and $V = ( v _ { i j } ) \in \mathbb { R } ^ { m \times m }$ .

Replace $Q _ { c }$ with the first j columns of $Q _ { c } V$ .

Replace $H _ { c }$ with $H _ { + } ( 1 { : } j , 1 { : } j )$ . .

Replace $r _ { c }$ with $v _ { m j } r _ { c }$ .

Starting with $A Q _ { c } = Q _ { c } H _ { c } + r _ { c } e _ { j } ^ { T }$ , perform steps $j + 1 , \ldots , j + p = m$ of the Arnoldi iteration obtaining $A Q _ { m } = Q _ { m } H _ { m } + r _ { m } e _ { m } ^ { T }$ .

Set $Q _ { c } = Q _ { m } , H _ { c } = H _ { m }$ , and $r _ { c } = r _ { m }$ .

In light of our remarks in §10.5.2, the filter values $\mu _ { 1 } , \ldots , \mu _ { p }$ should be chosen in the vicinity of $A \mathrm { { } s }$ “unwanted” eigenvalues. In this regard it is possible to formulate useful heuristics that are based on the eigenvalues of the m-by-m Hessenberg matrix $H _ { + }$ . For example, suppose the goal is to find the three smallest eigenvalues of A in absolute value. If $p = m - 3$ and $\lambda ( H _ { + } ) = \{ \tilde { \lambda } _ { 1 } , . . . , \tilde { \lambda } _ { m } \}$ with $| \tilde { \lambda } _ { 1 } \bar { | } \geq \cdots \geq | \tilde { \lambda } _ { m } |$ , then it is reasonable to set $\mu _ { i } = \tilde { \lambda } _ { i }$ for $i = 1 { : } p$ .

The Arnoldi iteration with implicit restarts has many attractive attributes. For implementation details and further analysis, see Lehoucq and Sorensen (1996), Morgan (1996), and the ARPACK manual by Lehoucq, Sorensen, and Yang (1998).

# 10.5.4 The Krylov-Schur Algorithm

An alternative restart procedure due to Stewart (2001) relies upon a carefully ordered Schur decomposition of the Hessenberg matrix $H _ { m }$ that is produced after m steps of the Arnoldi iteration. Suppose we have computed

$$
A Q _ {m} = Q _ {m} H _ {m} + r _ {m} e _ {m} ^ {T}
$$

and that $m = j + p ;$ , where $j$ is the number of A’s eigenvalues that we wish to compute. Let

$$
U ^ {T} H _ {m} U = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right]
$$

be the Schur decomposition of A and assume that the eigenvalues have been ordered so that the eigenvalues of $T _ { 1 1 } \in \mathbb { R } ^ { j \times j }$ are of interest and the eigenvalues of $T _ { 2 2 } \in \mathbb { R } ^ { p \times p }$ are not. (For clarity we ignore the possibility of complex eigenvalues.) The Arnoldi decomposition above transforms to

$$
A Q _ {+} = Q _ {+} T + r _ {c} e _ {m} ^ {T} U
$$

where $Q _ { + } = Q _ { m } U$ . It follows that

$$
A Q _ {+} (:, 1: j) = Q _ {+} (:, 1: j) T _ {1 1} + r _ {m} u ^ {T}
$$

where $u ^ { T } = U ( m , 1 { : } j )$ . It is possible to determine an orthogonal $\boldsymbol { Z } \in \mathbb { R } ^ { j \times j }$ so that $Z ^ { T } T _ { 1 1 } Z$ is upper Hessenberg and $Z ^ { T } u = \tau e _ { j }$ . (See P10.5.2.) It follows that

$$
A \left(Q _ {+} Z\right) = \left(Q _ {+} Z\right) \left(Z ^ {T} T _ {1 1} Z\right) + r _ {c} \left(Z ^ {T} u\right) ^ {T}
$$

is a j-step Arnoldi factorization. We then set $Q _ { j } , H _ { j }$ and $r _ { j }$ to be $Q _ { + } Z , Z ^ { T } T _ { 1 1 } Z$ , and $\tau r _ { m }$ respectively and perform Arnoldi steps $j + 1$ through $j + p = m$ . For more detailed discussion, see Stewart (MAE, Chap. 5) and Watkins (FMC, Chap. 9).

# 10.5.5 Unsymmetric Lanczos Tridiagonalization

Another way to extend the symmetric Lanczos process is to reduce A to tridiagonal form using a general similarity transformation. Suppose $A \in \mathbb { R } ^ { n \times n }$ and that a nonsingular matrix $Q$ exists such that

$$
Q ^ {- 1} A Q   =   T   =   \left[ \begin{array}{c c c c c} \alpha_ {1} & \gamma_ {1} & & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \gamma_ {n - 1} \\ 0 & \dots & & \beta_ {n - 1} & \alpha_ {n} \end{array} \right]  .
$$

With the column partitionings

$$
Q = \left[ q _ {1} \mid \dots \mid q _ {n} \right],
$$

$$
Q ^ {- T} = \tilde {Q} = \left[ \tilde {q} _ {1} | \dots | \tilde {q} _ {n} \right],
$$

we find upon comparing columns in $A Q { = } Q T$ and $A ^ { T } \tilde { Q } = \tilde { Q } T ^ { T }$ that

$$
A q _ {k} = \gamma_ {k - 1} q _ {k - 1} + \alpha_ {k} q _ {k} + \beta_ {k} q _ {k + 1}, \quad \gamma_ {0} q _ {0} \equiv 0,
$$

$$
A ^ {T} \tilde {q} _ {k} = \beta_ {k - 1} \tilde {q} _ {k - 1} + \alpha_ {k} \tilde {q} _ {k} + \gamma_ {k} \tilde {q} _ {k + 1}, \quad \beta_ {0} \tilde {q} _ {0} \equiv 0,
$$

for $k = 1 { : } n - 1$ . These equations together with the biorthogonality condition

$$
\tilde {Q} ^ {T} Q = I _ {n}
$$

imply

$$
\alpha_ {k} = \tilde {q} _ {k} ^ {T} A q _ {k}
$$

and

$$
\begin{array}{l} \beta_ {k} q _ {k + 1} \equiv r _ {k} = (A - \alpha_ {k} I) q _ {k} - \gamma_ {k - 1} q _ {k - 1}, \\ \gamma_ {k} \tilde {q} _ {k + 1} \equiv \tilde {r} _ {k} = (A - \alpha_ {k} I) ^ {T} \tilde {q} _ {k} - \beta_ {k - 1} \tilde {q} _ {k - 1}. \\ \end{array}
$$

There is some flexibility in choosing the scale factors $\beta _ { k }$ and $\gamma _ { k }$ . Note that

$$
1 = \tilde {q} _ {k + 1} ^ {T} q _ {k + 1} = \left(\tilde {r} _ {k} / \gamma_ {k}\right) ^ {T} \left(r _ {k} / \beta_ {k}\right).
$$

It follows that once $\beta _ { k }$ is specified, then $\gamma _ { k }$ is given by

$$
\gamma_ {k} = \tilde {r} _ {k} ^ {T} r _ {k} / \beta_ {k}.
$$

With the “canonical” choice $\beta _ { k } = \parallel r _ { k } \parallel _ { 2 }$ we obtain

$q _ { 1 } , \tilde { q } _ { 1 }$ given unit 2-norm vectors with $\tilde { q } _ { 1 } ^ { T } q _ { 1 } \neq 0$

$$
k = 0, q _ {0} = 0, r _ {0} = q _ {1}, \tilde {q} _ {0} = 0, s _ {0} = \tilde {q} _ {1}
$$

while $( r _ { k } \neq 0 )$ and $( \tilde { r } _ { k } \neq 0 )$ and $( \tilde { r } _ { k } ^ { T } r _ { k } \neq 0 )$

$$
\begin{array}{l} \beta_ {k} = \left\| r _ {k} \right\| _ {2} \\ \gamma_ {k} = \tilde {r} _ {k} ^ {T} r _ {k} / \beta_ {k} \\ q _ {k + 1} = r _ {k} / \beta_ {k} \\ \tilde {q} _ {k + 1} = \tilde {r} _ {k} / \gamma_ {k} \\ k = k + 1 \tag {10.5.11} \\ \end{array}
$$

$$
\alpha_ {k} = \tilde {q} _ {k} ^ {T} A q _ {k}
$$

$$
r _ {k} = (A - \alpha_ {k} I) q _ {k} - \gamma_ {k - 1} q _ {k - 1}
$$

$$
\tilde {r} _ {k} = (A - \alpha_ {k} I) ^ {T} \tilde {q} _ {k} - \beta_ {k - 1} \tilde {q} _ {k - 1}
$$

end

If

$$
T _ {k} = \left[ \begin{array}{c c c c c} \alpha_ {1} & \gamma_ {1} & & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \gamma_ {k - 1} \\ 0 & \dots & & \beta_ {k - 1} & \alpha_ {k} \end{array} \right],
$$

then the situation at the bottom of the loop is summarized by the equations

$$
A \left[ q _ {1} \mid \dots \mid q _ {k} \right] = \left[ q _ {1} \mid \dots \mid q _ {k} \right] T _ {k} + r _ {k} e _ {k} ^ {T}, \tag {10.5.12}
$$

$$
A ^ {T} \left[ \tilde {q} _ {1} \mid \dots \mid \tilde {q} _ {k} \right] = \left[ \tilde {q} _ {1} \mid \dots \mid \tilde {q} _ {k} \right] T _ {k} ^ {T} + \tilde {r} _ {k} e _ {k} ^ {T}. \tag {10.5.13}
$$

If $r _ { k } = 0 .$ , then the iteration terminates and span $\{ q _ { 1 } , \ldots , q _ { k } \}$ is an invariant subspace for A. If $\tilde { r } _ { k } = 0$ , then the iteration also terminates and span $\{ \tilde { q } _ { 1 } , \dots , \tilde { q } _ { k } \}$ is an invariant subspace for $A ^ { T }$ . However, if neither of these conditions is true and $\tilde { r } _ { k } ^ { T } r _ { k } = 0$ , then the tridiagonalization process ends without any invariant subspace information. This is called serious breakdown. See Wilkinson (AEP, p. 389) for an early discussion of the matter.

# 10.5.6 The Look-Ahead Idea

It is interesting to examine the serious breakdown issue in the block version of (10.5.11). For clarity assume that $A \in \mathbb { R } ^ { n \times n }$ with $n = r p$ . Consider the factorization in which we want $\tilde { Q } ^ { T } \dot { Q } = I _ { n }$ :

$$
\tilde {Q} ^ {T} A Q = \left[ \begin{array}{c c c c c} M _ {1} & C _ {1} ^ {T} & & \dots & 0 \\ B _ {1} & M _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & C _ {r - 1} ^ {T} \\ 0 & \dots & & B _ {r - 1} & M _ {r} \end{array} \right] \tag {10.5.14}
$$

where all the blocks are p-by-p. Let $Q = [  Q _ { 1 } | \cdot \cdot \cdot | Q _ { r } ]$ and $\tilde { Q } = [ \tilde { Q } _ { 1 } | \cdots | \tilde { Q } _ { r } ]$ be conformable partitionings of Q and Q˜. Comparing block columns in the equations $A Q = Q T$ and $\mathring { A } ^ { T } \tilde { Q } = \tilde { Q } T ^ { T }$ , we obtain

$$
Q _ {k + 1} B _ {k} = A Q _ {k} - Q _ {k} M _ {k} - Q _ {k - 1} C _ {k - 1} ^ {T} \equiv R _ {k},
$$

$$
\tilde {Q} _ {k + 1} C _ {k} = A ^ {T} \tilde {Q} _ {k} - \tilde {Q} _ {k} M _ {k} ^ {T} - \tilde {Q} _ {k - 1} B _ {k - 1} ^ {T} \equiv S _ {k}.
$$

Note that

$$
M _ {k} = \tilde {Q} _ {k} ^ {T} A Q _ {k}.
$$

If $S _ { k } ^ { T } R _ { k } = C _ { k } ^ { T } \tilde { Q } _ { k + 1 } ^ { T } Q _ { k + 1 } B _ { k } \in \mathbb { R } ^ { p \times p }$ is nonsingular and we compute $B _ { k } , C _ { k } \in \mathbb { R } ^ { p \times p }$ so that

$$
C _ {k} ^ {T} B _ {k} = S _ {k} ^ {T} R _ {k},
$$

then

$$
Q _ {k + 1} = R _ {k} B _ {k} ^ {- 1}, \tag {10.5.15}
$$

$$
\tilde {Q} _ {k + 1} = S _ {k} C _ {k} ^ {- 1} \tag {10.5.16}
$$

satisfy $\tilde { Q } _ { k + 1 } ^ { T } Q _ { k + 1 } = I _ { p }$ . Serious breakdown in this setting is associated with having a singular $S _ { k } ^ { T } R _ { k }$ .

One way of solving the serious breakdown problem in (10.5.11) is to go after a factorization of the form (10.5.14) in which the block sizes are dynamically determined. Roughly speaking, in this approach matrices $Q _ { k + 1 }$ and $\tilde { Q } _ { k + 1 }$ are built up column by column with special recursions that culminate in the production of a nonsingular $\tilde { Q } _ { k + 1 } ^ { T } Q _ { k + 1 }$ . The computations are arranged so that the biorthogonality conditions $\tilde { Q } _ { i } ^ { T } Q _ { k + 1 } = 0$ and $Q _ { i } ^ { T } \tilde { Q } _ { k + 1 } = 0$ hold for $i = 1 { : } k$ .

A method of this form belongs to the family of look-ahead Lanczos methods. The length of a look-ahead step is the width of the $Q _ { k + 1 }$ and $\tilde { Q } _ { k + 1 }$ that it produces. If that width is one, a conventional block Lanczos step may be taken. Length-2 lookahead steps are discussed in Parlett, Taylor, and Liu (1985). The notion of incurable breakdown is also presented by these authors. Freund, Gutknecht, and Nachtigal (1993) cover the general case along with a host of implementation details. Floating point considerations require the handling of “near” serious breakdown. In practice, each $M _ { k }$ that is 2-by-2 or larger corresponds to an instance of near serious breakdown.

# Problems

P10.5.1 Recalling how Theorem 10.1.1 establishes the orthogonality of the Lanczos vectors in Algorithm 10.1.1, state and prove an analogous theorem that does the same thing for the Arnoldi vectors in Algorithm 10.5.1.

P10.5.2 Show that if $C \in \mathbb { R } ^ { j \times j }$ and $u \in \mathbb { R } ^ { j }$ , then there exists an orthogonal $Z \in \mathbb { R } ^ { n \times n }$ so that $Z ^ { T } A Z \ = \ H$ is upper Hessenberg and the last column of Z is a multiple of u. Hint: Compute a Householder matrix P so that $P u$ is a multiple of $e _ { j }$ . Then reduce $C = \mathsf { \bar { P } } ^ { T } C P$ to upper Hessenberg form by producing a sequence of Householder updates $C = P _ { i } ^ { T } C P$ where $C ( n - i + 1 , 1 { : } n - i - 1 )$ is zeroed, $i = 1 { : } n - 2$ .

P10.5.3 Give an example of a starting vector for which the unsymmetric Lanczos iteration (10.5.11) breaks down without rendering any invariant subspace information. Use

$$
A = \left[ \begin{array}{c c c} 1 & 6 & 2 \\ 3 & 0 & 2 \\ 1 & 3 & 5 \end{array} \right].
$$

P10.5.4 Suppose $H \in \mathbb { R } ^ { n \times n }$ is upper Hessenberg. Discuss the computation of a unit upper triangular matrix U such that $H U = U T$ where T is tridiagonal.

P10.5.5 Show that the QR algorithm for eigenvalues does not preserve tridiagonal structure in the unsymmetric case.

# Notes and References for §10.5

For both analysis and implementation insight, Saad (NMLE) offers the most comprehensive treatment of unsymmetric Krylov methods. Stewart (MAE) and Watkins (MEP) devote entire chapters to the subject and are highly recommended as is the following review article:

D.C. Sorensen (2002). “Numerical Methods for Large Eigenvalue Problems,” Acta Numerica 11, 519–584.

The original Arnoldi idea first appeared in:

W.E. Arnoldi (1951). “The Principle of Minimized Iterations in the Solution of the Matrix Eigenvalue Problem,” Quarterly of Applied Mathematics 9, 17–29.

Saad set the stage for the development of practical implementations, see:

Y. Saad (1980). “Variations of Arnoldi’s Method for Computing Eigenelements of Large Unsymmetric Matrices.,” Lin. Alg. Applic. 34, 269–295.

Y. Saad (1984). “Chebyshev Acceleration Techniques for Solving Nonsymmetric Eigenvalue Problems,” Math. Comput. 42, 567–588.

Y. Saad (1989). “Krylov Subspace Methods on Supercomputers,” SIAM J. Sci. Stat. Comput., 10, 1200–1232.

References for implicit restarting in the Arnoldi context include:

D.C. Sorensen (1992). “Implicit Application of Polynomial Filters in a k-Step Arnoldi Method,” SIAM J. Matrix Anal. Applic. 13, 357–385.   
R.B. Lehoucq and D.C. Sorensen (1996). “Deflation Techniques for an Implicitly Restarted Iteration,” SIAM J. Matrix Anal. Applic. 17, 789–821.   
R.B. Morgan (1996). “On Restarting the Arnoldi Method for Large Nonsymmetric Eigenvalue Problems,” Math Comput. 65, 1213–1230.   
K. Meerbergen and A. Spence (1997). “Implicitly Restarted Arnoldi with Purification for the Shift-Invert Transformation,” Math. Comput. 218, 667–689.   
R.B. Lehoucq, D. C. Sorensen, and C. Yang (1998). ARPACK Users’ Guide: Solution of Large-Scale Eigenvalue Problems with Implicitly Restarted Arnoldi Methods, SIAM Publications, Philadelphia, PA.   
A. Stathopoulos, Y. Saad, and K. Wu (1998). “Dynamic Thick Restarting of the Davidson and the Implicitly Restarted Arnoldi Methods,” SIAM J. Sci. Comput. 19, 227–245.   
R.B. Lehoucq (2001). “Implicitly Restarted Arnoldi Methods and Subspace Iteration,” SIAM J. Matrix Anal. Applic. 23, 551–562.   
The Krylov-Schur approach to Arnoldi restarting is proposed in:   
G.W. Stewart (2001). “A Krylov-Schur Algorithm for Large Eigenproblems,” SIAM J. Matrix Anal. Applic., 23, 601–614.   
The rational Arnoldi process involves the shift-and-invert idea. In this framework Arnoldi is applied to the matrix $( A - \mu { \dot { I } } ) ^ { - 1 }$ , see:   
A. Ruhe (1984). “Rational Krylov Algorithms for Eigenvalue Computation,” Lin. Alg. Applic. 58, 391–405.   
A. Ruhe (1994). “Rational Krylov Algorithms for Nonsymmetric Eigenvalue Problems II. Matrix Pairs,” Lin. Alg. Applic. 197, 283–295.   
A. Ruhe (1994). “The Rational Krylov Algorithm for Nonsymmetric Eigenvalue Problems III: Complex Shifts for Real Matrices,” BIT 34, 165–176.   
Matrix function problems that involve large sparse matrices can be addressed using Krylov/Lanczos ideas, see:   
Y. Saad (1992). “Analysis of Some Krylov Subspace Approximations to the Matrix Exponential,” SIAM J. Numer. Anal. 29, 209–228.   
M. Hochbruck and C. Lubich (1997). “On Krylov Subspace Approximations to the Matrix Exponential Operator,” SIAM J. Numer. Anal. 34, 1911–1925.   
V. Druskin, A. Greenbaum and L. Knizhnerman (1998). “Using Nonorthogonal Lanczos Vectors in the Computation of Matrix Functions,” SIAM J. Sci. Comput. 19, 38–54.   
N. Del Buono, L. Lopez, and R. Peluso (2005). “Computation of the Exponential of Large Sparse Skew–Symmetric Matrices,” SIAM J. Sci. Comp. 27, 278–293.   
M. Eiermann and O.G. Ernst (2006). “A Restarted Krylov Subspace Method for the Evaluation of Matrix Functions,” SIAM J. Numer. Anal. 44, 2481–2504.   
J. van den Eshof and M. Hochbruck (2006). “Preconditioning Lanczos Approximations to the Matrix Exponential,” SIAM J. Sci. Comput. 27, 1438–1457.   
Other Arnoldi-related papers include:   
T. Huckle (1994). “The Arnoldi Method for Normal Matrices,” SIAM J. Matrix Anal. Applic. 15, 479–489.   
K.C. Toh and L.N. Trefethen (1996). “Calculation of Pseudospectra by the Arnoldi Iteration,” SIAM J. Sci. Comput. 17, 1–15.   
T.G. Wright and L.N. Trefethen (2001). “Large–Scale Computation of Pseudospectra Using ARPACK and Eigs,” SIAM J. Sci. Comput. 23, 591–605.   
V. Hernandez, J.E. Roman, and A. Tomas (2007). “Parallel Arnoldi Eigensolvers with Enhanced Scalability via Global Communications Rearrangement,” Parallel Computing 33, 521–540.   
The unsymmetric Lanczos process and related look ahead ideas are nicely presented in:   
B.N. Parlett, D. Taylor, and Z. Liu (1985). “A Look-Ahead Lanczos Algorithm for Unsymmetric Matrices,” Math. Comput. 44, 105–124.   
R.W. Freund, M. Gutknecht, and N. Nachtigal (1993). “An Implementation of the Look-Ahead Lanczos Algorithm for Non-Hermitian Matrices,” SIAM J. Sci. Stat. Comput. 14, 137–158.

# See also:

Y. Saad (1982). “The Lanczos Biorthogonalization Algorithm and Other Oblique Projection Methods for Solving Large Unsymmetric Eigenproblems,” SIAM J. Numer. Anal. 19, 485–506.   
D.L. Boley, S. Elhay, G.H. Golub and M.H. Gutknecht (1991) “Nonsymmetric Lanczos and Finding Orthogonal Polynomials Associated with Indefinite Weights,” Numer. Algorithms 1, 21–43.   
G.A. Geist (1991). “Reduction of a General Matrix to Tridiagonal Form,” SIAM J. Matrix Anal. Applic. 12, 362–373.   
C. Brezinski, M. Zaglia, and H. Sadok (1991). “Avoiding Breakdown and Near Breakdown in Lanczos Type Algorithms,” Numer. Algorithms 1, 261–284.   
S.K. Kim and A.T. Chronopoulos (1991). “A Class of Lanczos-Like Algorithms Implemented on Parallel Computers,” Parallel Comput. 17, 763–778.   
B.N. Parlett (1992). “Reduction to Tridiagonal Form and Minimal Realizations,” SIAM J. Matrix Anal. Applic. 13, 567–593.   
M. Gutknecht (1992). “A Completed Theory of the Unsymmetric Lanczos Process and Related Algorithms, Part I,” SIAM J. Matrix Anal. Applic. 13, 594–639.   
M. Gutknecht (1994). “A Completed Theory of the Unsymmetric Lanczos Process and Related Algorithms, Part II,” SIAM J. Matrix Anal. Applic. 15, 15–58.   
Z. Bai (1994). “Error Analysis of the Lanczos Algorithm for Nonsymmetric Eigenvalue Problem,” Math. Comput. 62, 209–226.   
T. Huckle (1995). “Low-Rank Modification of the Unsymmetric Lanczos Algorithm,” Math.Comput. 64, 1577–1588.   
Z. Jia (1995). “The Convergence of Generalized Lanczos Methods for Large Unsymmetric Eigenproblems,” SIAM J. Matrix Anal. Applic. 16, 543–562.   
M.T. Chu, R.E. Funderlic, and G.H. Golub (1995). “A Rank-One Reduction Formula and Its Applications to Matrix Factorizations,” SIAM Review 37, 512–530.

Computing eigenvalues of unsymmetric tridiagonal matrices is discussed in:

D.A. Bini, L. Gemignani, and F. Tisseur (2005). “The Ehrlich-Aberth Method for the Nonsymmetric Tridiagonal Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 27, 153–175.

# 10.6 Jacobi-Davidson and Related Methods

We close the chapter with a brief discussion of the Jacobi-Davidson method, a solution framework that involves a mix of several important ideas. The starting point is a reformulation of the eigenvalue problem as a nonlinear systems problem, a maneuver that enables us to apply Newton-like methods. This leads in a natural way to a method of Jacobi that can be used to compute eigenvalue-eigenvector pairs of symmetric matrices that have a strong diagonal dominance. Eigenproblems of this variety arise in quantum chemistry and it is in that venue where Davidson (1975) developed a very successful generalization of the Jacobi procedure. It builds a (non-Krylov) nested sequence of subspaces and incorporates Ritz approximation. By restricting the Davidson corrections to the orthogonal complement of the current subspace, we arrive at the Jacobi-Davidson method developed by Sleijpen and van der Vorst (1996). Their technique does not require symmetry or diagonal dominance. Thus, in terms of abstraction, exposition in this section starts from the general, descends to the specific, and then climbs back out to the general. All along the way we are driven by practical, algorithmic concerns. Our presentation draws upon the insightful treatments of the Jacobi-Davidson method in Sorensen (2002) and Stewart (MAE, pp. 404–420).

We mention that full appreciation of the Jacobi-Davidson method and its versatility requires an understanding of the next chapter. This is because a critical step in the method requires the approximate solution of a large sparse linear system and preconditioned iterative solvers are typically brought into play. See §11.5.

# 10.6.1 The Approximate Newton Framework

Consider the n-by-n eigenvalue problem $A x = \lambda x$ and how we might improve an approximate eigenpair $\{ x _ { c } , \lambda _ { c } \}$ . Note that if

$$
A (x _ {c} + \delta x _ {c}) = (\lambda_ {c} + \delta \lambda_ {c}) (x _ {c} + \delta x _ {c}),
$$

then

$$
\left(A - \lambda_ {c} I\right) \delta x _ {c} - \delta \lambda_ {c} x _ {c} = - r _ {c} + \delta \lambda_ {c} \cdot \delta x _ {c}, \tag {10.6.1}
$$

where

$$
r _ {c} = A x _ {c} - \lambda_ {c} x _ {c}
$$

is the current residual. By ignoring the second-order term $\delta \lambda _ { c } \cdot \delta x _ { c }$ we arrive at the following specification for the corrections $\delta x _ { c }$ and $\delta \lambda _ { c }$ :

$$
\left(A - \lambda_ {c} I\right) \delta x _ {c} - \delta \lambda_ {c} x _ {c} = - r _ {c}. \tag {10.6.2}
$$

This is an underdetermined system of nonlinear equations that has a very uninteresting solution obtained by setting $\delta x _ { c } = - x _ { c }$ and $\delta \lambda _ { c } = 0$ . To keep away from this situation we add a constraint so that if

$$
\left[ \begin{array}{l} x _ {+} \\ \lambda_ {+} \end{array} \right] = \left[ \begin{array}{l} x _ {c} \\ \lambda_ {c} \end{array} \right] + \left[ \begin{array}{l} \delta x _ {c} \\ \delta \lambda_ {c} \end{array} \right], \tag {10.6.3}
$$

then the new eigenvector approximation $x _ { + }$ is nonzero. One way to do this is to require

$$
w ^ {T} x _ {+} = 1,
$$

where $w \in \mathbb { R } ^ { n }$ is an appropriately chosen nonzero vector. Possibilities include $w = x$ , which forces $x _ { + }$ to have unit 2-norm, and $w = e _ { 1 }$ , which forces its first component to be one. Regardless, if $x _ { c }$ is also normalized with respect to $w .$ , then

$$
w ^ {T} \delta x _ {c} = w ^ {T} (x _ {+} - x _ {c}) = 0. \tag {10.6.4}
$$

By assembling (10.6.2) and (10.6.4) into a single matrix-vector equation we obtain

$$
\left[ \begin{array}{c c} A - \lambda_ {c} I & - x _ {c} \\ w ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} \delta x _ {c} \\ \delta \lambda_ {c} \end{array} \right] = - \left[ \begin{array}{l} r _ {c} \\ 0 \end{array} \right]. \tag {10.6.5}
$$

This is precisely the Jacobian system that arises if Newton’s method is used to find a zero of the function

$$
F \left(\left[ \begin{array}{l} x \\ \lambda \end{array} \right]\right) = \left[ \begin{array}{l} A x - \lambda x \\ w ^ {T} x - 1 \end{array} \right].
$$

Its solution is easy to specify:

$$
\delta \lambda_ {c} = \frac {w ^ {T} (A - \lambda_ {c} I) ^ {- 1} r _ {c}}{w ^ {T} (A - \lambda_ {c} I) ^ {- 1} x _ {c}}, \tag {10.6.6}
$$

$$
\delta x _ {c} = - \left(A - \lambda_ {c} I\right) ^ {- 1} \left(r _ {c} - \delta \lambda_ {c} x _ {c}\right). \tag {10.6.7}
$$

Unfortunately, the required linear equation solving is problematic if A is large and sparse and this prompts us to consider the approximate Newton framework.

The idea behind approximate Newton methods is to replace the Jacobian system with a nearby, look-alike system that is easier to solve. One way to do this in our problem is to approximate A with a matrix M with the proviso that systems of the form $( M - \lambda _ { c } I ) z = r$ are “easy” to solve. If $N = M - A$ , then (10.6.5) transforms to

$$
\left[ \begin{array}{c c} M - \lambda_ {c} I & - x _ {c} \\ w ^ {T} & 0 \end{array} \right] \left[ \begin{array}{c} \delta x _ {c} \\ \delta \lambda_ {c} \end{array} \right] = - \left[ \begin{array}{c} r _ {c} - N \cdot \delta x _ {c} \\ 0 \end{array} \right].
$$

Continuing with the approximate-Newton mentality, let us throw away the inconvenient $N { \cdot } \delta x _ { c }$ term that is part of the right-hand side. This leaves us with the system

$$
\left[ \begin{array}{c c} M - \lambda_ {c} I & - x _ {c} \\ w ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} \delta x _ {c} \\ \delta \lambda_ {c} \end{array} \right] = - \left[ \begin{array}{l} r _ {c} \\ 0 \end{array} \right], \tag {10.6.8}
$$

and the following compute-friendly recipes for the corrections:

$$
\delta \lambda_ {c} = \frac {w ^ {T} (M - \lambda_ {c} I) ^ {- 1} r _ {c}}{w ^ {T} (M - \lambda_ {c} I) ^ {- 1} x _ {c}}, \tag {10.6.9}
$$

$$
\delta x _ {c} = - (M - \lambda_ {c} I) ^ {- 1} \left(r _ {c} - \delta \lambda_ {c} x _ {c}\right). \tag {10.6.10}
$$

Of course, by cutting corners in Newton’s method we risk losing quadratic convergence. Thus, the design of an approximate Newton strategy must balance the efficiency of the approximate Jacobian solution procedure with a possibly degraded rate of convergence. For an excellent discussion of this tension in the context of the eigenvalue problem, see Stewart (MAE, pp. 396–404).

# 10.6.2 The Jacobi Orthogonal Component Correction Method

Now suppose

$$
A = \left[ \begin{array}{c c} \alpha & c ^ {T} \\ c & A _ {1} \end{array} \right], \quad \alpha \in \mathbb {R}, c \in \mathbb {R} ^ {n - 1}, A _ {1} \in \mathbb {R} ^ {(n - 1) \times (n - 1)} \tag {10.6.11}
$$

is symmetric and strongly diagonally dominant. Assume that α is the largest element on the diagonal in absolute value. Our ambition is to compute λ (close to $\alpha )$ and $z \in \mathbb { R } ^ { n - 1 }$ so that

$$
\left[ \begin{array}{c c} \alpha & c ^ {T} \\ c & A _ {1} \end{array} \right] \left[ \begin{array}{l} 1 \\ z \end{array} \right] = \lambda \left[ \begin{array}{l} 1 \\ z \end{array} \right]. \tag {10.6.12}
$$

Because of the dominance assumption, there is no danger in assuming that the soughtafter eigenvector is nicely normalized by setting its first component to 1. Partition $\delta \boldsymbol { x } _ { c } .$ , $x _ { c } ,$ and $x _ { + }$ as follows:

$$
\delta x _ {c} = \left[ \begin{array}{l} \delta \mu_ {c} \\ \delta z _ {c} \end{array} \right], \qquad x _ {c} = \left[ \begin{array}{l} 1 \\ z _ {c} \end{array} \right], \qquad x _ {+} = \left[ \begin{array}{l} 1 \\ z _ {+} \end{array} \right].
$$

By substituting (10.6.11) and $w = e _ { 1 }$ into the Jacobian system (10.6.5), we get

$$
\left[ \begin{array}{c|cc} \alpha - \lambda_ {c} & c ^ {T} & - 1 \\ c & A _ {1} - \lambda_ {c} I & - z _ {c} \\ \hline 1 & 0 & 0 \\ \end{array} \right]\left[ \begin{array}{c} \delta \mu_ {c} \\ \hline \delta z _ {c} \\ \delta \lambda_ {c} \\ \end{array} \right] = - \left[ \begin{array}{c} \alpha + + c ^ {T} z _ {c} - \lambda_ {c} \\ (A _ {1} - \lambda_ {c} I) z _ {c} + c \\ \hline 0 \\ \end{array} \right],
$$

i.e.,

$$
\left[ \begin{array}{c c} A _ {1} - \lambda_ {c} I & - z _ {c} \\ c ^ {T} & - 1 \end{array} \right] \left[ \begin{array}{l} \delta z _ {c} \\ \delta \lambda_ {c} \end{array} \right] = - \left[ \begin{array}{c} (A _ {1} - \lambda_ {c} I) z _ {c} + c \\ \alpha + c ^ {T} z _ {c} - \lambda_ {c} \end{array} \right]. \tag {10.6.13}
$$

It is easy to verify that this is the Jacobian system that arises if Newton’s method is used to compute a zero of

$$
f \left(\left[ \begin{array}{l} z \\ \lambda \end{array} \right]\right) = \left[ \begin{array}{l l} \alpha & c ^ {T} \\ c & A _ {1} \end{array} \right] \left[ \begin{array}{l} 1 \\ z \end{array} \right] - \lambda \left[ \begin{array}{l} 1 \\ z \end{array} \right].
$$

If $A _ { 1 } = M _ { 1 } - N _ { 1 }$ , then (10.6.13) can be rearranged as follows:

$$
(M _ {1} - \lambda_ {c} I) z _ {+} = - c + N _ {1} z _ {c} + \{\delta \lambda_ {c} \cdot z _ {c} + N _ {1} \cdot \delta z _ {c} \},
$$

$$
\lambda_ {+} = \alpha + c ^ {T} z _ {+}.
$$

The Jacobi orthogonal component correction $( J O C C )$ method is defined by ignoring the terms enclosed by the curly brackets and taking $M _ { 1 }$ to be the diagonal part of $A _ { 1 } \mathrm { { : } }$ :

$$
\lambda_ {1} = \alpha , z _ {1} = 0 _ {n - 1}, \rho_ {1} = \left\| c \right\| _ {2}, k = 1
$$

while $\rho _ { k } >$ tol

$$
(M _ {1} - \lambda_ {k} I) z _ {k + 1} = - c + N _ {1} z _ {k}
$$

$$
\lambda_ {k + 1} = \alpha + c ^ {T} z _ {k + 1} \tag {10.6.14}
$$

$$
k = k + 1
$$

$$
\rho_ {k} = \left\| A _ {1} z _ {k} - \lambda_ {k} z _ {k} + c \right\| _ {2}
$$

end

The name of the method stems from the fact that the corrections to the approximate eigenvectors

$$
x _ {k} = \left[ \begin{array}{c} 1 \\ z _ {k} \end{array} \right],
$$

are all orthogonal to $e _ { 1 }$ . Indeed, it is clear from (10.6.14) that each residual

$$
r _ {k} = (A - \lambda_ {k} I) x _ {k}
$$

has a zero first component:

$$
r _ {k} = \left[ \begin{array}{c c} \alpha & c ^ {T} \\ c & A _ {1} \end{array} \right] \left[ \begin{array}{c} 1 \\ z _ {k} \end{array} \right] - \lambda_ {k} \left[ \begin{array}{c} 1 \\ z _ {k} \end{array} \right] = \left[ \begin{array}{c} 0 \\ (A _ {1} - \lambda_ {k} I) z _ {k} + c \end{array} \right]. \tag {10.6.15}
$$

Hence, the termination criterion in (10.6.14) is based on the size of the residual.

Jacobi intended this method to be use in conjunction with his diagonalization procedure for the symmetric eigenvalue problem. As discussed in §8.5, after a sufficient number of sweeps the matrix A is very close to being diagonal. At that point, the JOCC iteration (10.6.14) can be invoked after a possible $\breve { P } A P ^ { \breve { T } }$ update to maximize the (1,1) entry.

# 10.6.3 The Davidson Method

As with the JOCC iteration, Davidson’s method is applicable to the symmetric diagonally dominant eigenvalue problem (10.6.12). However, it involves a more sophisticated placement of the residual vectors. To motivate the main idea, let M be the diagonal part of A and use (10.6.15) to rewrite the JOCC iteration as follows:

$$
x _ {1} = e _ {1}, \lambda_ {1} = x _ {1} ^ {T} A x _ {1}, r _ {1} = A x _ {1} - \lambda_ {1} x _ {1}, V _ {1} = [ e _ {1} ], k = 1
$$

while $\parallel r _ { k } \parallel > \mathsf { t o l }$

Solve the residual correction equation:

$$
(M - \lambda_ {k} I) \delta v _ {k} = - r _ {k}.
$$

Compute an improved eigenpair $\{ \lambda _ { k + 1 } , x _ { k + 1 } \}$ so $r _ { k + 1 } \in \mathsf { r a n } ( V _ { 1 } ) ^ { \perp } \colon$

$$
\delta x _ {k} = \delta v _ {k}, x _ {k + 1} = x _ {k} + \delta x _ {k}, \lambda_ {k + 1} = \lambda_ {k} + c ^ {T} \delta x _ {k}
$$

$$
k = k + 1
$$

$$
r _ {k} = A x _ {k} - \lambda_ {k} x _ {k}
$$

end

Davidson’s method uses Ritz approximation to ensure that $r _ { k }$ is orthogonal to $e _ { 1 }$ and $\delta v _ { 1 } , \ldots , \delta v _ { k - 1 }$ . To acomplish this, the boxed fragment is replaced with the following:

Expand the current subspace ran $( V _ { k } )$ :

$$
s _ {k + 1} = (I - V _ {k} V _ {k} ^ {T}) \delta v _ {k}
$$

$$
v _ {k + 1} = s _ {k + 1} / \left\| s _ {k + 1} \right\| _ {2}, V _ {k + 1} = \left[ V _ {k} \mid v _ {k + 1} \right] \tag {10.6.16}
$$

Compute an improved eigenpair {λk+1, xk+1} so rk+1 ∈ ran(Vk+1)⊥:

$$
(V _ {k + 1} ^ {T} A V _ {k + 1}) t _ {k + 1} = \theta_ {k + 1} t _ {k + 1} \quad \text {(a suitably chosen Ritz pair)}
$$

$$
\lambda_ {k + 1} = \theta_ {k + 1}, x _ {k + 1} = V _ {k + 1} t _ {k + 1}
$$

There are a number of important issues associated with this method. To begin with, $V _ { k }$ is an n-by-k matrix with orthonormal columns. The transition from $V _ { k }$ to $V _ { k + 1 }$ can be effectively carried out by a modified Gram-Schmidt process. Of course, if k gets too big, then it may be necessary to restart the process using $v _ { k }$ as the initial vector.

Because $r _ { k } = A x _ { k } - \lambda _ { k } x _ { k } = A ( V _ { k } t _ { k } ) - \theta _ { k } ( V _ { k } t _ { k } )$ , it follows that

$$
V _ {k} ^ {T} r _ {k} = (V _ {k} ^ {T} A V _ {k}) t _ {k} - \theta_ {k} t _ {k} = 0,
$$

i.e., $r _ { k }$ is orthogonal to the range of $V _ { k }$ as required.

We mention that the Davidson algorithm can be generalized by allowing M to be a more involved approximation to A than just its diagonal part. See Crouzeix, Philippe, and Sadkane (1994) for details.

# 10.6.4 The Jacobi-Davidson Framework

Instead of forcing the correction $\delta x _ { c }$ to be orthogonal to $e _ { 1 }$ as in the Davidson setting, the Jacobi-Davidson method insists that $\delta x _ { c }$ be orthogonal to the current eigenvector approximation $x _ { c } .$ . The idea is to expand the current search space in a profitable, unexplored direction.

To see what is involved computationally and to connect with Newton’s method, we consider the following modification of (10.6.5):

$$
\left[ \begin{array}{c c} A - \lambda_ {c} I & - x _ {c} \\ x _ {c} ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} \delta x _ {c} \\ \delta \lambda_ {c} \end{array} \right] = - \left[ \begin{array}{l} r _ {c} \\ 0 \end{array} \right]. \tag {10.6.17}
$$

Note that this is the Jacobian system associated with the function

$$
F \left(\left[ \begin{array}{l} x \\ \lambda \end{array} \right]\right) = \left[ \begin{array}{c} A x - \lambda x \\ (x ^ {T} x - 1) / 2 \end{array} \right]
$$

given that $x _ { c } ^ { T } x _ { c } = 1$ . If $x _ { c }$ is so normalized and $\lambda _ { c } = x _ { c } ^ { T } A x _ { c } ,$ then from (10.6.17) we have

$$
\begin{array}{l} (I - x _ {c} x _ {c} ^ {T}) (A - \lambda_ {c} I) (I - x _ {c} x _ {c} ^ {T}) \delta x _ {c} = - (I - x _ {c} x _ {c} ^ {T}) (r _ {c} - \delta \lambda_ {c} x _ {c}) \\ = - \left(I - x _ {c} x _ {c} ^ {T}\right) r _ {c} \\ = - \left(I - x _ {c} x _ {c} ^ {T}\right) \left(A x _ {c} - \lambda_ {c} x _ {c}\right) \\ = - \left(I - x _ {c} x _ {c} ^ {T}\right) A x _ {c} \\ = - \left(A x _ {c} - \lambda_ {c} x _ {c}\right) = - r _ {c}. \\ \end{array}
$$

Thus, the correction $\delta \boldsymbol { x } _ { c }$ is obtained by solving the projected system

$$
(I - x _ {c} x _ {c} ^ {T}) (A - \lambda_ {c} I) (I - x _ {c} x _ {c} ^ {T}) \delta x _ {c} = - r _ {c} \tag {10.6.18}
$$

subject to the constraint that $x _ { c } ^ { T } \delta x _ { c } = 0$ .

In Jacobi-Davidson, approximate projected systems are used to expand the current subspace. Compared to the Davidson algorithm, everything remains the same in (10.6.16) except that instead of solving $( M - \lambda _ { c } I ) \delta v _ { k } = - r _ { k }$ to determine $\delta v _ { k }$ , we solve

$$
(I - x _ {k} x _ {k} ^ {T}) (M - \lambda_ {k} I) (I - x _ {k} x _ {k} ^ {T}) \delta v _ {k} = - r _ {k}, \tag {10.6.19}
$$

subject to the constraint that $x _ { k } ^ { T } \delta v _ { k } = 0$ . The resulting framework permits greater flexibility. The initial unit vector x1 can be arbitrary and various Chapter 11 iterative solvers can be applied to (10.6.19). See Sleijpen and van der Vorst (1996) and Sorensen (2002) for details.

The Jacobi-Davidson framework can be used to solve both symmetric and nonsymmetric eigenvalue problems and is important for the way it channels sparse $A x = b$ technology to the sparse $A x = \lambda x$ problem. It can be regarded as an approximate Newton iteration that is “steered” to the eigenpair of interest by Ritz calculations. Because an ever-expanding orthonormal basis is maintained, restarting has a key role to play as in the Arnoldi setting (§10.5).

# 10.6.5 The Trace-Min Algorithm

We briefly discuss the trace-min algorithm that can be used to compute the k smallest eigenvalues and associated eigenvectors for the n-by-n symmetric-definite problem $A x = \lambda B x$ . It has similarities to the Jacobi-Davidson procedure. The starting point is to realize that if $V _ { \mathrm { o p t } } \in \mathbb { R } ^ { n \times k }$ solves

$$
\min _ {V ^ {T} B V = I _ {k}} \operatorname{tr} (V ^ {T} A V),
$$

then the required eigenvalues/eigenvectors are exposed by $V _ { \mathrm { o p t } } ^ { T } A V _ { \mathrm { o p t } } = \mathrm { d i a g } ( \mu _ { 1 } , \dots , \mu _ { k } )$ （2号 and $A V _ { \mathrm { o p t } } ( : , j ) = \mu _ { j } B V _ { \mathrm { o p t } } ( : , j )$ , for $j = 1 { : } k$ . The method produces a sequence of $V -$ matrices, each of which satisfies $V ^ { T } B V = I _ { k }$ . The transition from $V _ { c }$ to $V _ { + }$ requires the solution of a projected system

$$
(I - Q _ {c} Q _ {c} ^ {T}) A (I - Q _ {c} Q _ {c} ^ {T}) Z _ {c} = A V _ {c}
$$

where $Z _ { c } \in \mathbb { R } ^ { n \times k }$ and $Q R = B V _ { c }$ is the thin QR factorization. This system, analogous to the central Jacobi-Davidson update system (10.6.19), can be solved using a suitably preconditioned conjugate gradient iteration. For details, see Sameh and Wisniewski (1982) and Sameh and Tong (2000).

# Problems

P10.6.1 How would you solve (10.6.1) assuming that A is upper Hessenberg?

P10.6.2 Assume that

$$
A = \left[ \begin{array}{c c} \alpha & b \\ b & D + E \end{array} \right]
$$

is an n-by-n symmetric matrix. Assume that D is the diagonal of A(2:n, 2:n) and that the eigenvalue gap $\delta = \lambda _ { 1 } ( A ) - \lambda _ { 2 } ( A )$ is positive. How small must b and E be in order to ensure that (D + E) − αI is diagonally dominant? Use Theorem 8.1.4.

# Notes and References for §10.6

For deeper perspectives on the methods of this section, we recommend Stewart (MAE, 404–420) and:

D.C. Sorensen (2002). “Numerical Methods for Large Eigenvalue Problems,” Acta Numerica 11, 519–584.

Davidson method papers include:

E.R. Davidson (1975). “The Iterative Calculation of a Few of the Lowest Eigenvalues and Corresponding Eigenvectors of Large Real Symmetric Matrices,” J. Comput. Phys. 17, 87–94.

R.B. Morgan and D.S. Scott (1986). “Generalizations of Davidson’s Method for Computing Eigenvalues of Sparse Symmetric Matrices,” SIAM J. Sci. Stat. Comput. 7, 817–825.

J. Olsen, P. Jorgensen, and J. Simons (1990). “Passing the One-Billion Limit in Full-Configuration (FCI) Interactions,” Chem. Phys. Letters 169, 463–472.

R.B. Morgan (1992). “Generalizations of Davidson’s Method for Computing Eigenvalues of Large Nonsymmetric Matrices,” J. Comput. Phys. 101, 287–291.

M. Sadkane (1993) “Block-Arnoldi and Davidson Methods for Unsymmetric Large Eigenvalue Problems,” Numer. Math. 64, 195–211.   
M. Crouzeix, B. Philippe, and M. Sadkane (1994). “The Davidson Method,” SIAM J. Sci. Comput. 15, 62–76.   
A. Strathopoulos, Y. Saad, and C.F. Fischer (1995). “Robust Preconditioning for Large, Sparse, Symmetric Eigenvalue Problems,” J. Comput. Appl. Math. 64, 197–215.   
The original Jacobi-Davidson idea appears in:   
G.L.G. Sleijpen and H.A. van der Vorst (1996). “A Jacobi-Davidson Iteration Method for Linear Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 17, 401–425.   
For applications and extensions to other problems, see:   
G.L.G. Sleijpen, A.G.L. Booten, D.R. Fokkema, and H.A. van der Vorst (1996). “Jacobi-Davidson Type Methods for Generalized Eigenproblems and Polynomial Eigenproblems,” BIT 36, 595–633.   
G.L.G. Sleijpen, H.A. van der Vorst, and E. Meijerink (1998). “Efficient Expansion of Subspaces in the Jacobi-Davidson Method for Standard and Generalized Eigenproblems,” ETNA 7, 75–89.   
D.R. Fokkema, G.L.G. Sleijpen, and H.A. van der Vorst (1998). “Jacobi-Davidson Style QR and QZ Algorithms for the Reduction of Matrix Pencils,” SIAM J. Sci. Computut. 20, 94–125.   
P. Arbenz and M.E. Hochstenbach (2004). “A Jacobi-Davidson Method for Solving Complex Symmetric Eigenvalue Problems,” SIAM J. Sci. Comput. 25, 1655–1673.

The trace-min method is detailed in:

A. Sameh and J. Wisniewski (1982). “A Trace Minimization Algorithm for the Generalized Eigenproblem,” SIAM J. Numer. Anal. 19, 1243–1259.   
A. Sameh and Z. Tong (2000). “A Trace Minimization Algorithm for the Symmetric Generalized Eigenproblem,” J. Comput. Appl. Math. 123, 155–175.
