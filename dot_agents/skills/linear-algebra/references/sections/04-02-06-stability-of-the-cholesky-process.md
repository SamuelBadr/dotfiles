# 4.2.6 Stability of the Cholesky Process

In exact arithmetic, we know that a symmetric positive definite matrix has a Cholesky factorization. Conversely, if the Cholesky process runs to completion with strictly positive square roots, then A is positive definite. Thus, to find out if a matrix A is positive definite, we merely try to compute its Cholesky factorization using any of the methods given above.

The situation in the context of roundoff error is more interesting. The numerical stability of the Cholesky algorithm roughly follows from the inequality

$$
g _ {i j} ^ {2} \leq \sum_ {k = 1} ^ {i} g _ {i k} ^ {2} = a _ {i i}.
$$

This shows that the entries in the Cholesky triangle are nicely bounded. The same conclusion can be reached from the equation $\parallel G \parallel _ { 2 } ^ { 2 } = \parallel A \parallel _ { 2 }$ .

The roundoff errors associated with the Cholesky factorization have been extensively studied in a classical paper by Wilkinson (1968). Using the results in this paper, it can be shown that if ˆx is the computed solution to $A x = b$ , obtained via the Cholesky process, then ˆx solves the perturbed system

$$
(A + E) \hat {x} = b \quad \| E \| _ {2} \leq c _ {n} \mathbf {u} \| A \| _ {2},
$$

where $c _ { n }$ is a small constant that depends upon n. Moreover, Wilkinson shows that if $q _ { n } \mathbf { u } \kappa _ { 2 } ( A ) \leq 1$ where $q _ { n }$ is another small constant, then the Cholesky process runs to completion, i.e., no square roots of negative numbers arise.

It is important to remember that symmetric positive definite linear systems can be ill-conditioned. Indeed, the eigenvalues and singular values of a symmetric positive definite matrix are the same. This follows from (2.4.1) and Theorem 4.2.3. Thus,

$$
\kappa_ {2} (A) = \frac {\lambda_ {\max} (A)}{\lambda_ {\min} (A)}.
$$

The eigenvalue $\lambda _ { \mathrm { m i n } } ( A )$ is the “distance to trouble” in the Cholesky setting. This prompts us to consider a permutation strategy that steers us away from using small diagonal elements that jeopardize the factorization process.
