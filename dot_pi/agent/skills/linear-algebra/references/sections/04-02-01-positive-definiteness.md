# 4.2.1 Positive Definiteness

Suppose $A \in \mathbb { R } ^ { n \times n }$ is positive definite. It is obvious that a positive definite matrix is nonsingular for otherwise we could find a nonzero x so $x ^ { T } A x = 0$ . However, much more is implied by the positivity of the quadratic form $x ^ { T } A x$ as the following results show.

Theorem 4.2.1. If $A \in \mathbb { R } ^ { n \times n }$ is positive definite and $\boldsymbol { X } \in \mathbb { R } ^ { n \times k }$ has rank k, then $B = X ^ { T } A X \in \mathbb { R } ^ { k \times k }$ is also positive definite.

Proof. If $z \in \mathbb { R } ^ { k }$ satisfies $0 \geq z ^ { T } B z = ( X z ) ^ { T } A ( X z )$ , then $X z = 0$ . But since X has full column rank, this implies that $z = 0 . \qquad \bigtriangledown$

Corollary 4.2.2. If A is positive definite, then all its principal submatrices are positive definite. In particular, all the diagonal entries are positive.

Proof. If v is an integer length-k vector with $1 \leq v _ { 1 } < \cdot \cdot \cdot < v _ { k } \leq n$ , then $X = I _ { n } ( : , v )$ （ is a rank-k matrix made up of columns $v _ { 1 } , \ldots , v _ { k }$ of the identity. It follows from Theorem 4.2.1 that $A ( v , v ) = X ^ { T } A X$ is positive definite.

Theorem 4.2.3. The matrix $A \in \mathbb { R } ^ { n \times n }$ is positive definite if and only if the symmetric matrix

$$
T = \frac {A + A ^ {T}}{2}
$$

has positive eigenvalues.

Proof. Note that $x ^ { T } A x = x ^ { T } T x$ . If $T x = \lambda x$ then $x ^ { T } A x = \lambda \cdot x ^ { T } x$ . Thus, if A is positive definite then λ is positive. Conversely, suppose T has positive eigenvalues and $Q ^ { T } T Q = \mathrm { d i a g } ( \lambda _ { i } )$ is its Schur decomposition. (See §2.1.7.) It follows that if $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $y = Q ^ { T } x$ , then

$$
x ^ {T} A x = x ^ {T} T x = y ^ {T} (Q ^ {T} T Q) y = \sum_ {k = 1} ^ {n} \lambda_ {k} y _ {k} ^ {2} > 0,
$$

completing the proof of the theorem.

□

Corollary 4.2.4. If A is positive definite, then it has an LU factorization and the diagonal entries of U are positive.

Proof. From Corollary 4.2.2, it follows that the submatrices $A ( 1 { : } k , 1 { : } k )$ are nonsingular for $k = 1 { : } n$ and so from Theorem 3.2.1 the factorization $A = L U$ exists. If we apply Theorem 4.2.1 with $X = ( L ^ { - 1 } ) ^ { T } = L ^ { - T }$ , then $B = X ^ { T } A X = L ^ { - 1 } ( L U ) L ^ { - 1 } = U \bar { L } ^ { - \bar { T } }$ is positive definite and therefore has positive diagonal entries. The corollary follows because $L ^ { - T }$ is unit upper triangular and this implies $b _ { i i } = u _ { i i } , ~ i = 1 { : } n$ . □

The mere existence of an LU factorization does not mean that its computation is advisable because the resulting factors may have unacceptably large elements. For example, if $\epsilon > 0$ , then the matrix

$$
A = \left[ \begin{array}{c c} \epsilon & m \\ - m & \epsilon \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ - m / \epsilon & 1 \end{array} \right] \left[ \begin{array}{c c} \epsilon & m \\ 0 & 1 + m ^ {2} / \epsilon \end{array} \right]
$$

is positive definite. However, if $m / \epsilon \gg 1$ , then it appears that some kind of pivoting is in order. This prompts us to pose an interesting question. Are there conditions that guarantee when it is safe to compute the LU-without-pivoting factorization of a positive definite matrix?
