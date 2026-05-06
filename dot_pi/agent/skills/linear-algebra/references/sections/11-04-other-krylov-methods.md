# 11.4 Other Krylov Methods

The conjugate gradient method can be regarded as a clever pairing of the symmetric Lanczos process and the $L D L ^ { T }$ factorization. The “cleverness” is associated with the recursions that support an economical transition from $x _ { k - 1 }$ to $x _ { k }$ . In this section we move beyond symmetric positive definite systems and present instances of the same paradigm for more general problems:

$$
\binom{\text {Krylov}}{\text {process}} + \binom{\text {Matrix}}{\text {factorization}} + \binom{\text {Clever}}{\text {recursions}} = \left( \begin{array}{c} \text {Sparse} \\ \text {matrix} \\ \text {method} \end{array} \right).
$$

Methods for the symmetric indefinite problem (MINRES, SYMMLQ), the least squares problem (LSQR, LSMR), and the square Ax = b problem (GMRES, QMR, BiCG, CGS, BiCGStab) are briefly discussed. The Lanczos, Arnoldi, and unsymmetric Lanczos iterations are in the mix. Our goal is to communicate the main idea behind these methods. For deeper insight, practical intuition, and analysis, see Saad (ISPLA), Greenbaum (IMSL), van der Vorst (IMK), Freund, Golub, and Nachtigal (1992), and LIN TEMPLATES.

# 11.4.1 MINRES and SYMMLQ for Symmetric Systems

Assume that $A \in \mathbb { R } ^ { n \times n }$ is symmetric indefinite, i.e., $\lambda _ { \operatorname* { m i n } } ( A ) < 0 < \lambda _ { \operatorname* { m a x } } ( A )$ . A consequence of this is that we cannot recast the Ax = b problem as a minimization problem associated with $\phi ( x ) = x ^ { T } A x / 2 - x ^ { T } b .$ Indeed, this function has no lower bound. If $A x = \lambda _ { \operatorname* { m i n } } x$ , then $\phi ( \alpha x ) = \alpha ^ { 2 } \lambda _ { \operatorname* { m i n } } - \alpha x ^ { T } b$ approaches −∞ as α gets big.

This suggests that we switch to a more workable objective function. Instead of adopting the CG strategy of minimizing $\phi$ over the affine space $x _ { 0 } + \mathcal { K } ( A , r _ { 0 } , k )$ , we propose to solve

$$
\min _ {x \in x _ {0} + \mathcal {K} (A, r _ {0}, k)} \| b - A x \| _ {2}. \tag {11.4.1}
$$

at each step. As in CG, we use the Lanczos process to generate the Krylov subspaces, setting $q _ { 1 } = r _ { 0 } / \beta _ { 0 }$ where $r _ { 0 } = b - A x _ { 0 }$ and $\beta _ { 0 } = \parallel g _ { 0 } \parallel _ { 2 }$ . After k steps we have

$$
A Q _ {k} = Q _ {k} T _ {k} + \beta_ {k} q _ {k + 1} e _ {k} ^ {T}.
$$

That is,

$$
A Q _ {k} = Q _ {k + 1} H _ {k}, \tag {11.4.2}
$$

where $H _ { k } \in \mathbb { R } ^ { k + 1 \times k }$ is the Hessenberg matrix

$$
H _ {k} = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {2} & \dots & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & 0 \\ \vdots & \ddots & \ddots & & \vdots \\ \vdots & & & \ddots & \beta_ {k - 1} \\ 0 & \dots & \dots & \beta_ {k - 1} & \alpha_ {k} \\ \hline 0 & \dots & \dots & 0 & \beta_ {k} \end{array} \right]. \tag {11.4.3}
$$

Writing $x \ : = \ : x _ { 0 } + Q _ { k } y$ and recalling that ${ \mathsf { r a n } } ( Q _ { k } ) \ = \ K ( A , r _ { 0 } , k )$ , we see that the optimization (11.4.1) involves minimizing

$$
\parallel A (x _ {0} + Q _ {k} y) - b \parallel_ {2} = \parallel Q _ {k + 1} H _ {k} y - (b - A x _ {0}) \parallel_ {2} = \parallel H _ {k} y - \beta_ {0} e _ {1} \parallel_ {2}
$$

over all $\boldsymbol { y } \in \mathbb { R } ^ { k }$ . To solve this problem we take a hint from §5.2.6 and use the Givens QR factorization procedure. Suppose $G _ { 1 } , \ldots , G _ { k }$ are Givens rotations such that

$$
G _ {k} ^ {T} \cdot \cdot \cdot G _ {1} ^ {T} H _ {k} = \left[ \frac {R _ {k}}{0} \right], \qquad R _ {k} \in \mathbb {R} ^ {k \times k},
$$

is upper triangular. If

$$
G _ {k} ^ {T} \dots G _ {1} ^ {T} (\beta_ {0} e _ {1}) = \left[ \frac {p _ {k}}{\rho_ {k}} \right], \qquad p _ {k} \in \mathbb {R} ^ {k},
$$

and $y _ { k } \in \mathbb { R } ^ { k }$ solves $R _ { k } y _ { k } = p _ { k }$ , then $x _ { k } = x _ { 0 } + Q _ { k } y _ { k }$ solves (11.4.1) and the norm of the residual is given by  $\vert b - A x _ { k } \vert \vert _ { 2 } = \vert \rho _ { k } \vert$ . The transition

$$
\left\{H _ {k - 1}, R _ {k - 1}, p _ {k - 1}, \rho_ {k - 1} \right\} \quad \rightarrow \quad \left\{H _ {k}, R _ {k}, p _ {k}, \rho_ {k} \right\}
$$

can be realized with O(1) flops after the kth Lanczos step is performed. The Givens rotation $G _ { k }$ can be determined from $\beta _ { k }$ and $[ R _ { k - 1 } ] _ { k - 1 , k - 1 }$ . Note that after step k−1 we already have the first k − 2 rows of $R _ { k }$ and the first k−2 components of $p _ { k }$ . The matrix $R _ { k }$ has upper bandwidth 2 and so the triangular system that determines $y _ { k }$ can be solved with O(k) flops. Thus, in computing $x _ { k } = x _ { 0 } + Q _ { k } y _ { k }$ each step is not essential. On the other hand, it is possible to work out an $O ( n )$ transition from $x _ { k - 1 }$ to $x _ { k }$ through recursions that involve $Q _ { k }$ and the QR factorization of $H _ { k }$ . (This corresponds to the $L D L ^ { T } \mathrm { \ - p l u s - } Q _ { k }$ recursions associated with CG developed in §11.3.5.) Either way, there is no need to access all the Lanczos vectors each step. Properly implemented, we have the MINRES method of Paige and Saunders (1975).

An alternative approach developed by the same authors works with the LQ factorization of the tridiagonal matrix $T _ { k }$ . We mimic the §11.3.4 in the CG derivation leading to (11.3.14). However, the solution of the tridiagonal system

$$
T _ {k} y _ {k} = \beta_ {0} e _ {1} \tag {11.4.4}
$$

is problematic because $T _ { k }$ is no longer positive definite. This means that the $\mathrm { L D L } ^ { T }$ factorization, together with the associated recursions, is no longer safe to use.

A way around this difficulty is to work with the transpose of the matrix equation $A Q _ { k - 1 } = Q _ { k } H _ { k - 1 }$ . Suppose $x _ { k } = x _ { 0 } + Q _ { k } y _ { k }$ where $y _ { k }$ is the minimum-norm solution to the $( k - 1 ) – \mathrm { b y } – k$ underdetermined system

$$
H _ {k - 1} ^ {T} y _ {k} = \beta_ {0} e _ {1}. \tag {11.4.5}
$$

It follows from $r _ { 0 } = \beta _ { 0 } Q _ { k - 1 } e _ { 1 } , r _ { k } = r _ { 0 } - A Q _ { k - 1 } y _ { k }$ , and $Q _ { k - 1 } ^ { T } A = H _ { k - 1 } ^ { T } Q _ { k } ^ { T }$ that

$$
Q _ {k - 1} ^ {T} r _ {k} = \beta_ {0} e _ {1} - H _ {k - 1} ^ {T} y _ {k} = 0.
$$

Thus, the residual $r _ { k } = b - A x _ { k }$ is orthogonal to $q _ { 1 } , \ldots , q _ { k - 1 }$ . Note that the underdetermined system (11.4.5) has full row rank and that $y _ { k }$ can be determined via a Givens rotation lower triangularization, e.g.,

$$
\left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & 0 & 0 & 0 \\ \beta_ {1} & \alpha_ {2} & \beta_ {2} & 0 & 0 \\ 0 & \beta_ {2} & \alpha_ {3} & \beta_ {3} & 0 \\ 0 & 0 & \beta_ {3} & \alpha_ {4} & \beta_ {4} \end{array} \right] G _ {1} G _ {2} G _ {3} G _ {4} = \left[ \begin{array}{c c c c c} \times & 0 & 0 & 0 & 0 \\ \times & \times & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 \\ 0 & \times & \times & \times & 0 \end{array} \right] = \left[ \begin{array}{c c c c c} L _ {4} & 0 \end{array} \right].
$$

This is an LQ factorization and in general we have

$$
H _ {k - 1} ^ {T} G _ {1} \dots G _ {k - 1} = \left[ L _ {k - 1} \mid 0 \right]
$$

where $L _ { k - 1 }$ is lower triangular. (This is just the transpose of the Givens QR factorization of $H _ { k - 1 } . )$ If $w _ { k - 1 } \in \mathbf { \bar { R } } ^ { k - 1 }$ solves the necessarily nonsingular system $L _ { k - 1 } w _ { k - 1 } =$ $\beta _ { 0 } e _ { 1 }$ , then

$$
y _ {k} = G _ {1} \dots G _ {k - 1} \left[ \begin{array}{c} w _ {k - 1} \\ 0 \end{array} \right].
$$

The special structure of $L _ { k - 1 }$ (it has lower bandwidth equal to 2) and the Givens rotation sequence make it possible to realize the transition from $x _ { k }$ to $x _ { k + 1 }$ with $O ( n )$ work in a way that does not require access to all the Lanczos vectors. Collectively, these ideas define the SYMMLQ method of Paige and Saunders (1975).

# 11.4.2 LSQR and LSMR for Least Squares Problems

We show how the sparse least squares problem min $A x - b \parallel _ { 2 }$ can be solved using the Paige-Saunders lower bidiagonalization process described in §10.4.4. Indeed, if we apply Algorithm 10.4.2 with $u _ { 1 } = r _ { 0 } / \beta _ { 0 }$ where $r _ { 0 } = b - A x _ { 0 }$ and $\beta _ { 0 } = \parallel r _ { 0 } \parallel _ { 2 }$ , then after k steps we have a partial factorization of the form

$$
A V _ {k} = U _ {k} B _ {k} + p _ {k} e _ {k} ^ {T}
$$

where $V = [  v _ { 1 } | \cdot \cdot \cdot | v _ { k } ] \in \mathbb { R } ^ { n \times k }$ has orthonormal columns, $\boldsymbol { U } = \left[ \boldsymbol { u } _ { 1 } \vert \cdot \cdot \cdot \vert \boldsymbol { u } _ { k } \right] \in \mathbb { R } ^ { m \times k }$ has orthonormal columns, and $B _ { k } \in \mathbb { R } ^ { k \times k }$ is lower bidiagonal. If $p _ { k } \in \mathbb { R } ^ { m }$ is nonzero, then we can write

$$
A V _ {k} = U _ {k + 1} \tilde {B} _ {k}
$$

where $\tilde { B } _ { k } \in \mathbb { R } ^ { k + 1 \times k }$ is given by

$$
\tilde {B} _ {k} = \left[ \begin{array}{c c c c c} \alpha_ {1} & 0 & \dots & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & 0 \\ \vdots & \ddots & & & \vdots \\ \vdots & & & \ddots & 0 \\ 0 & \dots & \dots & \beta_ {k - 1} & \alpha_ {k} \\ \hline 0 & \dots & \dots & 0 & \beta_ {k} \end{array} \right]. \tag {11.4.6}
$$

It can be shown that span $\{ v _ { 1 } , . . . , v _ { k } \} \ = \ K ( A ^ { T } A , A ^ { T } r _ { 0 } , k )$ . In the LSQR method of Paige and Saunders (1982), the kth approximate minimizer $x _ { k }$ solves the problem

$$
\min _ {x \in x _ {0} + \mathcal {K} (A ^ {T} A, A ^ {T} r _ {0}, k)} \| A x - b \| _ {2}. \tag {11.4.7}
$$

Thus, $x _ { k } = x _ { 0 } + V _ { k } y _ { k }$ where $y _ { k } \in \mathbb { R } ^ { k }$ is the minimizer of

$$
\| A (x _ {0} + V _ {k} y) - b \| _ {2} = \| U _ {k + 1} \tilde {B} _ {k} y - (b - A x _ {0}) \| _ {2} = \| \tilde {B} _ {k} y - \beta_ {0} e _ {1} \| _ {2}.
$$

Givens QR can be used to solve this problem just as it is used in the MINRES context above. Suppose

$$
G _ {k} ^ {T} \dots G _ {1} ^ {T} \tilde {B} _ {k} = \left[ \frac {R _ {k}}{0} \right], \qquad G _ {k} ^ {T} \dots G _ {1} ^ {T} (\beta_ {1} e _ {1}) = \left[ \frac {p _ {k}}{\rho_ {k}} \right],
$$

where $G _ { 1 } , \ldots , G _ { k }$ are Givens rotations, $R _ { k } \in \mathbb { R } ^ { k \times k }$ is upper triangular, $p _ { k } \in \mathbb { R } ^ { k }$ , and $\rho _ { k } \in \mathbb { R }$ . Then, $y _ { k }$ solves $R _ { k } y = p _ { k }$ and

$$
x _ {k} = x _ {0} + V _ {k} y _ {k} = x _ {0} + W _ {k} p _ {k}
$$

where $W _ { k } = V _ { k } R _ { k } ^ { - 1 }$ . It is possible to compute $x _ { k }$ from $x _ { k - 1 }$ via a simple recursion that involves the last column of $W _ { k }$ . Overall, we obtain the LSQR method of Paige and Saunders (1982). It requires only a few vectors of storage to implement.

The LSMR method provides an alternative to the LSQR method and is mathematically equivalent to MINRES applied to the normal equations $A ^ { T } A x = A ^ { T } b$ . Like LSQR, the technique can be used to solve least squares problems, regularized least squares problems, undetermined systems, and square unsymmetric systems. The 2- norms of the vectors $r _ { k } = b - A x _ { k }$ and $A ^ { T } r _ { k }$ decrease monotonically, which allows for tractable early-termination. See Fong and Saunders (2011) for more details.

# 11.4.3 GMRES for General $A x = b$

The Paige-Saunders MINRES method (§11.4.1) is a Lanczos-based technique that can be used to solve symmetric Ax = b problems. The kth iterate $x _ { k }$ minimizes $\parallel A x - b \parallel _ { 2 }$ over $x _ { 0 } + \mathcal { K } ( A , b , k )$ . We now present an Arnoldi-based iteration that does the same thing and is applicable to general linear systems. The method is referred to as the generalized minimum residual (GMRES) method and is due to Saad and Shultz (1986).

After k steps of the Arnoldi iteration (Algorithm 10.5.1) it is easy to confirm using (10.5.2) that

$$
A Q _ {k} = Q _ {k + 1} \tilde {H} _ {k} \tag {11.4.8}
$$

where the columns of

$$
Q _ {k + 1} = \left[ Q _ {k} \mid q _ {k + 1} \right]
$$

are the orthonormal Arnoldi vectors and the upper Hessenberg matrix $\tilde { H } _ { k }$ is given by

$$
\tilde {H} _ {k} = \left[ \begin{array}{c c c c c} h _ {1 1} & h _ {1 2} & \dots & \dots & h _ {1 k} \\ h _ {2 1} & h _ {2 2} & \dots & \dots & h _ {2 k} \\ 0 & \ddots & \ddots & & \vdots \\ \vdots & & \ddots & \ddots & \vdots \\ 0 & \dots & \dots & h _ {k, k - 1} & h _ {k k} \\ 0 & \dots & \dots & 0 & h _ {k + 1, k} \end{array} \right] \in \mathbb {R} ^ {k + 1 \times k}.
$$

Moreover, if $q _ { 1 } ~ = ~ r _ { 0 } / \beta _ { 0 }$ where $r _ { 0 } = b - A x _ { 0 }$ and $\beta _ { 0 } ~ = ~ \parallel r _ { 0 } \parallel _ { 2 }$ , then

$$
\operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\} = \mathcal {K} (A, r _ {0}, k).
$$

In step k, the GMRES method requires minimization of $\parallel A x - b \parallel _ { 2 }$ over the affine space $x _ { 0 } + { \mathcal { K } } ( A , r _ { 0 } , k )$ . As with MINRES, we must find a vector $\boldsymbol { y } \in \mathbb { R } ^ { k }$ so that

$$
\parallel A (x _ {0} + Q _ {k} y) - b \parallel_ {2} = \parallel Q _ {k + 1} \tilde {H} _ {k} y - (b - A x _ {0}) \parallel_ {2} = \parallel \tilde {H} _ {k} y - \beta_ {0} e _ {1} \parallel_ {2}
$$

is minimized. If $y _ { k }$ is the solution to this $( k + 1 ) – \mathrm { b y } – k$ least squares problem, then the k-th GMRES iterate is given by $x _ { k } = x _ { 0 } + Q _ { k } y _ { k }$ . Note that if Givens rotations $G _ { 1 } , \ldots , G _ { k }$ have been determined so that

$$
G _ {k} ^ {T} \dots G _ {1} ^ {T} \tilde {H} _ {k} = \left[ \frac {R _ {k}}{0} \right], \quad R _ {k} \in \mathbb {R} ^ {k \times k}, \tag {11.4.9}
$$

is upper triangular and we set

$$
G _ {k} ^ {T} \dots G _ {1} ^ {T} \left(\beta_ {0} e _ {1}\right) = \left[ \frac {p _ {k}}{\rho_ {k}} \right], \tag {11.4.10}
$$

where $p _ { k } \in \mathbb { R } ^ { k }$ and $\rho _ { k } \in \mathbb { R }$ , then $R _ { k } y _ { k } = p _ { k }$ and

$$
| \rho_ {k} | = \parallel A x _ {k} - b \parallel_ {2}.
$$

The transition

$$
\{R _ {k - 1}, p _ {k - 1}, \rho_ {k - 1} \} \rightarrow \{R _ {k}, p _ {k}, \rho_ {k} \}
$$

is a particularly simple update that involves the generation of a single rotation $G _ { k }$ and exploitation of the identities $R _ { k - 1 } = R _ { k } ( 1 { : } k - 1 , 1 { : } k - 1 )$ and $p _ { k } ( 1 ; k - 1 ) = p _ { k - 1 }$ .

As a procedure for large sparse problems, the GMRES method inherits the usual Arnoldi concern: the computation of $H ( 1 ; k + 1 , k )$ requires $O ( k n )$ flops and access to all previously computed Arnoldi vectors. For this reason it is necesssary to build a restart strategy around the following, m-step GMRES building block:

Algorithm 11.4.2 (m-step GMRES) If $A \in \mathbb { R } ^ { n \times n }$ is nonsingular, $b \in \mathbb { R } ^ { n }$ , $A x _ { 0 } \approx b ,$ , and m is a positive iteration limit, then this algorithm computes $\widetilde { \boldsymbol { x } } \in \mathbb { R } ^ { n }$ where either $\tilde { x }$ solves $A x = b$ or minimizes $\parallel A x - b \parallel _ { 2 }$ over the affine space $x _ { 0 } + \mathcal { K } ( A , r _ { 0 } , m )$ where $r _ { 0 } = b - A x _ { 0 }$ .

$$
k = 0, r _ {0} = b - A x _ {0}, \beta_ {0} = \parallel r _ {0} \parallel_ {2}
$$

while $( \beta _ { k } > 0 )$ and $k < m$

$$
q _ {k + 1} = r _ {k} / \beta_ {k}
$$

$$
k = k + 1
$$

$$
r _ {k} = A q _ {k}
$$

for i = 1:k

$$
h _ {i k} = q _ {i} ^ {T} r _ {k} \tag {11.4.11}
$$

$$
r _ {k} = r _ {k} - h _ {i k} q _ {i}
$$

end

$$
\beta_ {k} = \parallel r _ {k} \parallel_ {2}
$$

$$
h _ {k + 1, k} = \beta_ {k}
$$

Apply $G _ { 1 } , \ldots , G _ { k - 1 }$ to $H ( 1 { : } k , k )$ and determine $G _ { k } , R _ { k } , p _ { k }$ , and $\rho _ { k }$ end

Solve $R _ { k } y _ { k } = p _ { k }$ and set $\tilde { x } = x _ { 0 } + Q _ { k } y _ { k }$

If ˜x is not good enough, then the process can be repeated with the new $x _ { 0 }$ set to ${ \tilde { x } } .$ . There are many important implementation details associated with this framework, see Saad (IMSLA, pp. 164–184) and van der Vorst (IMK, pp. 65–84).

# 11.4.4 Optimizing from the Polynomial Point of View

Before we present the next group of methods, it is instructive to connect the Krylov framework with polynomial approximation. Suppose the columns of $Q _ { k } \in \mathbb { R } ^ { n \times k }$ span $\kappa ( A , q _ { 1 } , k )$ . It follows that if $\bar { y } \in \mathbb { R } ^ { k }$ , then $Q _ { k } y = \varphi ( A ) q _ { 1 }$ for some polynomial $\varphi$ that has degree $k - 1$ or less. This is because

$$
Q _ {k} = \left[ q _ {1} \mid A q _ {1} \mid \dots \mid A ^ {k - 1} q _ {1} \right] B
$$

for some nonsingular $B \in \mathbb { R } ^ { k \times k }$ and so if $\alpha = B y$ , then

$$
Q _ {k} y = \left[ q _ {1} \mid A q _ {1} \mid \dots \mid A ^ {k - 1} q _ {1} \right] \alpha = \left(\alpha_ {1} I + \alpha_ {2} A + \dots + \alpha_ {k} A ^ {k - 1}\right) q _ {1}.
$$

Thus, the GMRES (and MINRES) optimization can be rephrased as a polynomial optimization problem. If $\mathbb { P } _ { k }$ denotes the set of all degree-k polynomials, then we have

$$
\begin{array}{l} \min \quad \| b - A x \| _ {2} = \min \quad \| b - A (x _ {0} + \varphi (A)) r _ {0} \| _ {2} \\ x \in x _ {0} + \mathcal {K} (A, r _ {0}, k) \quad \varphi \in \mathbb {P} _ {k - 1} \\ = \min _ {\varphi \in \mathbb {P} _ {k - 1}} \| (I - A \cdot \varphi (A)) r _ {0} \| _ {2} \\ = \min _ {\psi \in \mathbf {P} _ {k}, \psi (0) = 1} \| \psi (A) r _ {0} \| _ {2}. \\ \end{array}
$$

This point of view figures heavily in the analysis of various Krylov subspace methods and can also be used to suggest alternative strategies.

# 11.4.5 BiCG, CGS, BiCGstab, and QMR for General $A x = b$

Just as the Arnoldi iteration underwrites GMRES, the unsymmetric Lanczos process (10.5.11) underwrites the next cohort of methods that we present. Suppose we complete k steps of (10.5.11) with $q _ { 1 } = r _ { 0 } / \beta _ { 0 } , r _ { 0 } = b - A x _ { 0 } , \beta _ { 0 } = \| r _ { 0 } \| _ { 2 } ,$ , and $r _ { 0 } ^ { T } \tilde { r } _ { 0 } \neq 0$ . This means we have the partial factorizations

$$
A Q _ {k} = Q _ {k} T _ {k} + r _ {k} e _ {k} ^ {T}, \quad \tilde {Q} _ {k} ^ {T} r _ {k} = 0, \tag {11.4.12}
$$

$$
A ^ {T} \tilde {Q} _ {k} = \tilde {Q} _ {k} T _ {k} ^ {T} + \tilde {r} _ {k} e _ {k} ^ {T}, \quad Q _ {k} ^ {T} \tilde {r} _ {k} = 0, \tag {11.4.13}
$$

where

$$
Q _ {k} = \left[ q _ {1} \mid \dots \mid q _ {k} \right], \quad \operatorname{ran} \left(Q _ {k}\right) = \mathcal {K} \left(A, r _ {0}, k\right),
$$

$$
\tilde {Q} _ {k} = \left[ \tilde {q} _ {1} \mid \dots \mid \tilde {q} _ {k} \right], \quad \operatorname{ran} \left(\tilde {Q} _ {k}\right) = \mathcal {K} \left(A ^ {T}, \tilde {r} _ {0}, k\right).
$$

In addition, $\tilde { Q } _ { k } ^ { T } Q _ { k } \ = \ I _ { k }$ and $\tilde { Q } _ { k } ^ { T } A Q _ { k } = T _ { k } \in \mathbb { R } ^ { k \times k }$ is tridiagonal. Vectors $q _ { k + 1 }$ and $\tilde { q } _ { k + 1 }$ and scalars $\beta _ { k }$ and $\tau _ { k }$ satisfy

$$
\beta_ {k} q _ {k + 1} = r _ {k}, \quad \tau_ {k} \tilde {q} _ {k + 1} = \tilde {r} _ {k}
$$

and can be generated with access to just the last two columns of $Q _ { k }$ and $\tilde { Q } _ { k }$ .

In step k of the biconjugate gradient (BiCG) method, an iterate $x _ { k } = x _ { 0 } + Q _ { k } y _ { k }$ is produced where $y _ { k } \in \mathbb { R } ^ { k }$ solves the k-by-k tridiagonal system

$$
T _ {k} y _ {k} = \tilde {Q} _ {k} ^ {T} r _ {0}.
$$

It follows that

$$
\tilde {Q} _ {k} ^ {T} (b - A x _ {k}) = \tilde {Q} _ {k} ^ {T} (b - A (x _ {0} + Q _ {k} y _ {k})) = \tilde {Q} _ {k} ^ {T} r _ {0} - T _ {k} y _ {k} = 0.
$$

Thus, the residual associated with $x _ { k }$ is orthogonal to the range of $\tilde { Q } _ { k }$ .

Assume that $T _ { k }$ has an $L U$ factorization $T _ { k } = L _ { k } U _ { k }$ and note that $L _ { k }$ is unit lower bidiagonal and $U _ { k }$ is upper bidiagonal. It follows that

$$
x _ {k} = x _ {0} + Q _ {k} T _ {k} ^ {- 1} \tilde {Q} _ {k} ^ {T} r _ {0} = (Q _ {k} U _ {k} ^ {- 1}) (L _ {k} ^ {- 1} (\tilde {Q} _ {k} ^ {T} r _ {0})).
$$

Analogously to how we derived the CG algorithm, it is possible to develop simple connections between the matrix $( Q _ { k } U _ { k } ^ { - 1 } )$ and its predecessor and between the vector $( L _ { k } ^ { - 1 } ( \tilde { Q } _ { k } ^ { T } r _ { 0 } ) )$ and its predecessor. The end result is a procedure that can generate $x _ { k }$ through simple recursions, which we report in Figure 11.4.1. We mention that the BiCG method is subject to serious breakdown because of its dependence on the unsymmetric Lanczos process. However, with the look-ahead idea discussed in §10.5.6, it is possible to overcome some of these difficulties. Notice that BiCG collapses to CG if A is symmetric positive definite and $\tilde { r } _ { 0 } = r _ { 0 }$ . Also observe the similarity between the r and $\tilde { r }$ updates and the $p$ and $\tilde { p }$ updates.

A negative aspect of the BiCG method is that it requires procedures for both A-times-vector and $\overrightharpoon { A ^ { T } - } \mathrm { t i m e s - v e c t o r } .$ . (In some applications the latter is a challenge.)

<table><tr><td>BiCG</td><td>CGS</td><td>BiCGstab</td></tr><tr><td> $r_0 = b - Ax_0$ </td><td> $r_0 = b - Ax_0$ </td><td> $r_0 = b - Ax_0$ </td></tr><tr><td> $\tilde{r}_0^T r_0 \neq 0$ </td><td> $\tilde{r}_0^T \tilde{r}_0 \neq 0$ </td><td> $\tilde{r}_0^T \tilde{r}_0 \neq 0$ </td></tr><tr><td> $x_c = x_0$ </td><td> $x_c = x_0$ </td><td> $x_c = x_0$ </td></tr><tr><td> $p_c = r_c = r_0$ </td><td> $p_c = r_c = r_0$ </td><td> $p_c = r_c = r_0$ </td></tr><tr><td> $\tilde{p}_c = \tilde{r}_c = \tilde{r}_0$ </td><td> $u_c = r_c$ </td><td></td></tr><tr><td> $\mu = \frac{\tilde{r}_c^T r_c}{\tilde{p}_c^T A p_c}$ </td><td> $\mu = \frac{\tilde{r}_0^T r_c}{\tilde{r}_0^T A p_c}$ </td><td> $\mu = \frac{\tilde{r}_0^T r_c}{\tilde{r}_0^T A p_c}$ </td></tr><tr><td> $x_+ = x_c + \mu p_c$ </td><td> $q_c = u_c - \mu A p_c$ </td><td> $s_c = r_c - \mu A p_c$ </td></tr><tr><td> $r_+ = r_c - \mu A p_c$ </td><td> $x_+ = x_c + \mu (u_c + q_c)$ </td><td> $\omega = \frac{s_c^T A s_c}{(A s_c)^T (A s_c)}$ </td></tr><tr><td> $\tilde{r}_+ = \tilde{r}_c - \mu A^T \tilde{p}_c$ </td><td> $r_+ = r_c - \mu A (u_c + q_c)$ </td><td> $x_+ = x_c + \mu p_c + \omega s_c$ </td></tr><tr><td> $\tau = \frac{\tilde{r}_+^T r_+}{\tilde{r}_c^T r_c}$ </td><td> $\tau = \frac{\tilde{r}_0^T r_+}{\tilde{r}_0^T r_c}$ </td><td> $r_+ = s_c - \omega A s_c$ </td></tr><tr><td> $p_+ = r_+ + \tau p_c$ </td><td> $u_+ = r_+ + \tau q_c$ </td><td> $\tau = \frac{(\tilde{r}_0^T r_+) \mu}{(\tilde{r}_0^T r_c) \omega}$ </td></tr><tr><td> $\tilde{p}_+ = \tilde{r}_+ + \tau \tilde{p}_c$ </td><td> $p_+ = u_+ + \tau (q_c + \tau p_c)$ </td><td> $p_+ = r_+ + \tau (p_c - \omega A p_c)$ </td></tr></table>

Figure 11.4.1. The initializations and update formulae for the biconjugate gradient (BiCG) method, the conjugate gradient squared (CGS) method, and the biconjugate gradient stablilized (BiCGstab) method. The subscript $^ { 6 6 } c ^ { 9 9 }$ designates “current” while the subscript $^ { 6 6 } + { } ^ { , 9 }$ designates “next”.

The conjugate gradient squared (CGS) method circumvents this problem and has some interesting convergence properties as well. The derivation of the method uses the polynomial point of view that we outlined in the previous section. It is easy to conclude from Figure 11.4.1 that after k steps of the procedure we have degree-k polynomials $\psi _ { k }$ and $\varphi _ { k }$ so that

$$
\begin{array}{l} r _ {k} = \psi_ {k} (A) r _ {0}, \quad p _ {k} = \varphi_ {k} (A) r _ {0}, \\ \tilde {\tau} = \psi_ {k} (A T) \tilde {\tau}, \quad \tilde {\tau} = \varphi_ {k} (A T) \tilde {\tau}. \end{array} \tag {11.4.14}
$$

$$
\tilde {r} _ {k} = \psi_ {k} (A ^ {T}) \tilde {r} _ {0}, \quad \tilde {p} _ {k} = \varphi_ {k} (A ^ {T}) \tilde {r} _ {0},
$$

and $\psi _ { k } ( 0 ) = \varphi _ { k } ( 0 ) = 1$ . This enables us to characterize expressions like $\tilde { r } _ { k } ^ { T } r _ { k }$ and $\tilde { p } _ { k } ^ { T } A p _ { k }$ in a way that involves only A-times-vector:

$$
\tilde {r} _ {k} ^ {T} r _ {k} = \left(\psi_ {k} (A ^ {T}) \tilde {r} _ {0}\right) ^ {T} (\psi_ {k} (A) r _ {0}) = \tilde {r} _ {0} ^ {T} \left(\psi_ {k} ^ {2} (A) r _ {0}\right),
$$

$$
\tilde {p} _ {k} ^ {T} A p _ {k} = \left(\varphi_ {k} (A ^ {T}) \tilde {r} _ {0}\right) ^ {T} A \left(\varphi_ {k} (A) r _ {0}\right) = \tilde {r} _ {0} ^ {T} \left(A \varphi_ {k} ^ {2} (A) r _ {0}\right).
$$

It is possible to develop simple recursions among the polynomials $\{ \psi _ { k } \}$ and $\{ \varphi _ { k } \}$ that facilitate the transitions

$$
r _ {k - 1} = \psi_ {k - 1} ^ {2} (A) r _ {0} \rightarrow \psi_ {k} ^ {2} (A) r _ {0} = r _ {k},
$$

$$
p _ {k - 1} = \varphi_ {k - 1} ^ {2} (A) r _ {0} \rightarrow \varphi_ {k} ^ {2} (A) r _ {0} = p _ {k}.
$$

This leads to the conjugate gradient squared (CGS) method of Sonneveld (1989). It produces iterates $x _ { k }$ whose residuals $r _ { k }$ satisfy $r _ { k } = \dot { \psi } _ { k } ( A ) ^ { 2 } r _ { 0 }$ . Note from Figure 11.4.1 that the updates rely on only matrix-vector products that involve only A. Because of the squaring of the BiCG residual polynomial $\psi _ { k }$ , the method typically outperforms BiCG when it works, i.e., $( \| \psi _ { k } ( A ) ^ { 2 } r _ { 0 } \| _ { 2 } \ll \| \psi _ { k } ( A ) r _ { 0 } \| _ { 2 } )$ . By the same token, it typically underperforms when BiCG struggles.

A third member in this family of $A x = b$ solvers is the BiCGstab method of van der Vorst (1992). It addresses the sometimes erratic behavior of BiCG by producing iterates $x _ { k }$ whose residuals satisfy

$$
r _ {k} = (1 - \omega_ {k} A) \dots (1 - \omega_ {1} A) \psi_ {k} (A) r _ {0}
$$

where $\psi _ { k }$ is the BiCG residual polynomial defined in (11.4.14). The parameter $\omega _ { k }$ is chosen in step $k$ to minimize $\| r _ { k } \| _ { 2 }$ given $\omega _ { 1 } , \ldots , \omega _ { k - 1 }$ and the vector $\psi _ { k } ( A ) r _ { 0 }$ . The computations associated with this transpose-free method are given in Figure 11.4.1.

Yet another iteration that is built upon the unsymmetric Lanczos process is the quasi-minimum residual (QMR) method of Freund and Nachtigal (1991). As in BiCG, the kth iterate has the form $x _ { k } = x _ { 0 } + Q _ { k } y _ { k }$ where $Q _ { k }$ is specified by (11.4.12). This equation can be rewritten as $A Q _ { k } = Q _ { k + 1 } \tilde { T } _ { k }$ where $\tilde { T } _ { k } \in \mathbf { \bar { \mathbb { R } } } ^ { k + 1 \times k }$ is tridiagonal. It follows that if $q _ { 1 } = r _ { 0 } / \beta _ { 0 }$ where $r _ { 0 } = b - A x _ { 0 }$ and $\beta _ { 0 } = \parallel r _ { 0 } \parallel _ { 2 }$ , then

$$
b - A \left(x _ {0} + Q _ {k} y\right) = r _ {0} - A Q _ {k} y = r _ {0} - Q _ {k + 1} \tilde {T} _ {k} y = Q _ {k + 1} \left(\beta_ {0} e _ {1} - \tilde {T} _ {k} y\right).
$$

In QMR, $y$ is chosen to minimize $\parallel \beta _ { 0 } e _ { 1 } - \tilde { T } _ { k } y \parallel _ { 2 }$ . Note that GMRES minimizes the same quantity because $Q _ { k + 1 }$ has orthonormal columns in Arnoldi.

# Problems

P11.4.1 Assume that the cost of a length-n inner product or saxpy is one unit. Assume that $A \in \mathbb { R } ^ { n \times n }$ and that the matrix-vector products involving A and $A ^ { T }$ cost α and $\beta$ units, respectively. Compare the per iteration cost associated with the BiCG, CGS, and BiCGstab methods.

P11.4.2 Suppose $A \in \mathbb { R } ^ { n \times n }$ and $v \in \mathbb { R } ^ { n }$ are given. How can we choose ω to minimize $\parallel ( I - \omega A ) v \parallel _ { 2 } ?$

P11.4.3 Give an algorithm that computes $\psi _ { k } ( a )$ where $a \in \mathbb { R }$ and $\psi _ { k }$ is defined by (11.4.14).

# Notes and References for §11.4

For general systems, we have avoided the when-to-use-what-method question because there are no clear-cut answers. For guidance we recommend LIN TEMPLATES, Greenbaum (IMSL), Saad (ISPLA), and van der Vorst (IKM), each of which provides a great deal of insight. See also:

R.W. Freund, G.H. Golub, and N.M. Nachtigal (1992). “Iterative Solution of Linear Systems,” Acta Numerica 1, 57–100.

The MINRES, SYMMLQ, and LSQR frameworks due to Paige and Saunders initiated one of the most important threads of Krylov method research:

C.C. Paige and M.A. Saunders (1975). “Solution of Sparse Indefinite Systems of Linear Equations,” SIAM J. Numer. Anal. 12, 617–629.   
C.C. Paige and M.A. Saunders (1982). “LSQR: An Algorithm for Sparse Linear Equations and Sparse Least Squares,” ACM Trans. Math. Softw. 8, 43–71.   
M.A. Saunders, H.D. Simon, and E.L. Yip (1988). “Two Conjugate-Gradient Type Methods for Unsymmetric Linear Systems,” SIAM J. Numer. Anal. 25, 927–940.   
C.C. Paige, B.N. Parlett, and H.A. van der Vorst (1995). “Approximate Solutions and Eigenvalue Bounds for Krylov Subspaces,” Numer. Lin. Alg. Applic. 3, 115–133.   
M.A. Saunders (1997). “Computing Projections with LSQR,” BIT 37, 96–104.   
F.A. Dul (1998). “MINRES and MINERR Are Better than SYMMLQ in Eigenpair Computations,” SIAM J. Sci. Comput. 19, 1767–1782.   
S.J. Benbow (1999). “Solving Generalized Least-Squares Problems with LSQR,” SIAM J. Matrix Anal. Applic. 21, 166–177.   
M. Kilmer and G.W. Stewart (2000). “Iterative Regularization and MINRES,” SIAM J. Matrix Anal. Appl. 21, 613–628.   
L. Reichel and Q. Ye (2008). “A Generalized LSQR Algorithm,” Numer. Lin. Alg. Applic. 15, 643–660.   
X.-W. Chang, C.C. Paige, and D. Titley-Peloquin (2009). “Stopping Criteria for the Iterative Solution of Linear Least Squares Problems,” SIAM J. Matrix Anal. Applic. 31, 831–852.   
S.-C. Choi, C.C. Paige, and M.A. Saunders (2011). “MINRES-QLP: A Krylov Subspace Method for Indefinite or Singular Symmetric Systems,” SIAM J. Sci. Comput. 33, 1810–1836.   
D.C.-L. Fong and M.A. Saunders (2011). “LSMR: An Iterative Algorithm for Sparse Least-Squares Problems,” SIAM J. Sci. Comput. 33, 2950–2971.

The original GMRES paper is set forth in:

Y. Saad and M. Schultz (1986). “GMRES: A Generalized Minimum Residual Algorithm for Solving Unsymmetric Linear Systems,” SIAM J. Sci. Stat. Comput. 7, 856–869.

and there is a great deal of follow-up analysis:

S.L. Campbell, I.C.F. Ipsen, C.T. Kelley, and C.D. Meyer (1996). “GMRES and the Minimal Polynomial,” BIT 36, 664–675.

A. Greenbaum, V. Ptak, and Z. Strakoˇs (1996). “Any Nonincreasing Convergence Curve is Possible for GMRES,” SIAM J. Matrix Anal. Applic. 17, 465–469.

K.-C. Toh (1997). “GMRES vs. Ideal GMRES,” SIAM J. Matrix Anal. Applic. 18, 30–36.

M. Arioli, V. Ptak, and Z. Strakoˇs (1998). “Krylov Sequences of Maximal Length and Convergence of GMRES,” BIT 38, 636–643.

Y. Saad (2000). “Further Analysis of Minimum Residual Iterations,” Numer. Lin. Alg. 7, 67–93.

I.C.F. Ipsen (2000). “Expressions and Bounds for the GMRES Residual,” BIT 40, 524–535.

D. Calvetti, B. Lewis, and L. Reichel (2002). “On the Regularizing Properties of the GMRES Method,” Numer. Math. 91, 605–625.

J. Liesen, M. Rozloznik, and Z. Strakoˇs (2002). “Least Squares Residuals and Minimal Residual Methods,” SIAM J. Sci. Comput. 23, 1503–1525.

J. Liesen and P. Tich´y (2004). “The Worst-Case GMRES for Normal Matrices,” BIT 44, 79–98.

C.C. Paige, M. Rozloznik, and Z. Strakoˇs (2006). “Modified Gram-Schmidt (MGS), Least Squares, and Backward Stability of MGS-GMRES,” SIAM J. Matrix Anal. Applic 28, 264–284.

For pseudosprectral analysis of the method, see Trefethen and Embree (SAP, Chap. 26) as well as

M. Embree (1999). “Convergence of Krylov Subspace Methods for Non-Normal Matrices,” PhD Thesis, Oxford University.

References concerned with the critical issue of restarting include:

R.B. Morgan (1995). “A Restarted GMRES Method Augmented with Eigenvectors,” SIAM J. Matrix Anal. Applic. 16, 1154–1171.

A. Frommer and U. Glassner (1998). “Restarted GMRES for Shifted Linear Systems,” SIAM J. Sci. Comput. 19, 15–26.

V. Simoncini (1999). “A New Variant of Restarted GMRES,” Numer. Lin. Alg. 6, 61–77.

R.B. Morgan (2000). “Implicitly Restarted GMRES and Arnoldi Methods for Nonsymmetric Systems of Equations,” SIAM J. Matrix Anal. Applic. 21, 1112–1135.

K. Moriya and T. Nodera (2000). “The DEFLATED-GMRES(m,k) Method with Switching the Restart Frequency Dynamically,” Numer. Lin. Alg. 7, 569–584.   
J. Zitko (2000). “Generalization of Convergence Conditions for a Restarted GMRES,” Numer. Lin. Alg. 7, 117–131.   
R.B. Morgan (2002). “GMRES with Deflated Restarting,” SIAM J. Sci. Comput. 24, 20–37.   
M. Embree (2003). “The Tortoise and the Hare Restart GMRES,” SIAM Review 45, 259–266.   
J. Zitko (2004). “Convergence Conditions for a Restarted GMRES Method Augmented with Eigenspaces,” Numer. Lin. Alg. 12, 373–390.

Various practical issues concerned with GMRES implementation are covered in:

H.F. Walker (1988). “Implementation of the GMRES Method Using Householder Transformations,” SIAM J. Sci. Stat. Comput. 9, 152–163.

A. Greenbaum, M. Rozloznik, and Z. Strako (1997). “Numerical Behaviour of the Modified Gram-Schmidt GMRES Implementation,” BIT 37, 706–719.

P.N. Brown and H.F. Walker (1997). “GMRES On (Nearly) Singular Systems,” SIAM J. Matrix Anal. Applic. 18, 37–51.

K. Burrage and J. Erhel (1998). “On the Performance of Various Adaptive Preconditioned GMRES Strategies,” Numer. Lin. Alg. 5, 101–121.

Y. Saad and K. Wu (1998). “DQGMRES: a Direct Quasi-minimal Residual Algorithm Based on Incomplete Orthogonalization,” Numer. Lin. Alg. 3, 329–343.

M. Sosonkina, L.T. Watson, R.K. Kapania, and H.F. Walker (1999). “A New Adaptive GMRES Algorithm for Achieving High Accuracy,” Numer. Lin. Alg. 5, 275–297.

J. Liesen (2000). “Computable Convergence Bounds for GMRES,” SIAM J. Matrix Anal. Applic. 21, 882–903.

V. Frayss, L. Giraud, S. Gratton, and J. Langou (2005). “Algorithm 842: A Set of GMRES Routines for Real and Complex Arithmetics on High Performance Computers,” ACM Trans. Math. Softw. 31, 228–238.

A.H. Baker, E.R. Jessup and T. Manteuffel (2005). “A Technique for Accelerating the Convergence of Restarted GMRES,” SIAM J. Matrix Anal. Applic. 26, 962–984.

L. Reichel and Q. Ye (2005). “Breakdown-free GMRES for Singular Systems,” SIAM J. Matrix Anal. Applic. 26, 1001–1021.

There is a block version of the GMRES method, see:

V. Simoncini and E. Gallopoulos (1996). “Convergence Properties of Block GMRES and Matrix Polynomials,” Lin. Alg. Applic. 247, 97–119.

A.H. Baker, J.M. Dennis, and E.R. Jessup (2006). “On Improving Linear Solver Performance: A Block Variant of GMRES,” SIAM J. Sci. Comput. 27, 1608–1626.

M. Robb and M. Sadkane (2006). “Exact and Inexact Breakdowns in the Block GMRES Method,” Lin. Alg. Applic. 419, 265–285.

Original references associated with the BiCG, CGS, QMR, and BiCGstab methods include:

C. Lanczos (1952). “Solution of Systems of Linear Equations by Minimized Iterations,” J. Res. Nat. Bur. Stand. 49, 33-53.

R. Fletcher (1975). “Conjugate Gradient Methods for Indefinite Systems,” in Proceedings of the Dundee Biennial Conference on Numerical Analysis, 1974, G.A. Watson (ed), Springer-Verlag, New York.

P. Sonneveld (1989). “CGS: A Fast Lanczos-Type Solver for Nonsymmetric Linear Systems,” SIAM J. Sci. Stat. Comput. 10, 36–52.

R. Freund and N. Nachtigal (1991). “QMR: A Quasi-Minimal Residual Method for Non-Hermitian Linear Systems,” Numer. Math. 60, 315–339.

H.A. van der Vorst (1992). “Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems,” SIAM J. Sci. Stat. Comput. 13, 631–644.

Subsequent papers that pertain to these methods include:

G.L.G. Sleijpen and D.R. Fokkema (1993). “BiCGstab(l) for Linear Equations Involving Unsymmetric Matrices with Complex Spectrum,” ETNA 1, 11–32.

R. Freund (1993). “A Transpose Free Quasi-Minimum Residual Algoroithm for Non-Hermitian Linear Systems,” SIAM J. Sci. Comput. 14, 470–482.

R.W. Freund and N.M. Nachtigal (1996). “QMRPACK: a Package of QMR Algorithms,” ACM Trans. Math. Softw. 22, 46–77.   
M.-C. Yeung and T.F. Chan (1999). “ML(k)BiCGSTAB: A BiCGSTAB Variant Based on Multiple Lanczos Starting Vectors,” SIAM J. Sci. Comput. 21, 1263–1290.   
M. Kilmer, E. Miller, and C. Rappaport (2001). “QMR-Based Projection Techniques for the Solution of Non–Hermitian Systems with Multiple Right–Hand Sides,” SIAM J. Sci. Comput. 23, 761–780.   
A. El Guennouni, K. Jbilou, and H. Sadok (2003). “A Block Version of BiCGSTAB for Linear Systems with Multiple Right-Hand Sides,” ETNA 16, 129–142.   
G.L.G. Sleijpen, P. Sonneveld, and M.B. van Gijzen (2009). “BiCGSTAB as an Induced Dimension Reduction Method,” Appl. Numer. Math. 60, 1100–1114.   
M.H. Gutknecht (2010). “IDR Explained,” ETNA 36, 126–148.
