# 4.2.3 Symmetric Positive Definite Systems

If we apply the above results to a symmetric positive definite matrix we know that the factorization $A = L U$ exists and is stable to compute. The computation of the factorization $A = L D L ^ { T }$ via Algorithm 4.1.2 is also stable and exploits symmetry. However, for symmetric positive definite systems it is often handier to work with a variation of $L \dot { D } L ^ { T }$ .

Theorem 4.2.7 (Cholesky Factorization). If $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite, then there exists a unique lower triangular $G \in \mathbb { R } ^ { n \times n }$ with positive diagonal entries such that $A = G G ^ { T }$ .

Proof. From Theorem 4.1.3, there exists a unit lower triangular L and a diagonal

$$
D = \mathrm{diag} (d _ {1}, \ldots , d _ {n})
$$

such that $A = L D L ^ { T }$ . Theorem 4.2.1 tells us that $L ^ { - 1 } A L ^ { - T } = D$ is positive definite. Thus, the $d _ { k }$ are positive and the matrix $G = L \operatorname { d i a g } ( { \sqrt { d _ { 1 } } } , \ldots , { \sqrt { d _ { n } } } )$ is real and lower triangular with positive diagonal entries. It also satisfies $A = G G ^ { T }$ . Uniqueness follows from the uniqueness of the $\mathrm { \Delta L D L ^ { T } }$ factorization.

The factorization $A = G G ^ { T }$ is known as the Cholesky factorization and G is the Cholesky factor. Note that if we compute the Cholesky factorization and solve the triangular systems $G y = b$ and $G ^ { T } x = y .$ , then $b = G y = G ( G ^ { T } x ) = ( G G ^ { T } ) x = A x$ .
