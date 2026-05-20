# 4.2.7 The LDLT Factorization with Symmetric Pivoting

With an eye towards handling ill-conditioned symmetric positive definite systems, we return to the $\mathrm { L D L } ^ { T }$ factorization and develop an outer product implementation with pivoting. We first observe that if A is symmetric and $P _ { 1 }$ is a permutation, then $P _ { 1 } A$ is not symmetric. On the other hand, $P _ { 1 } A P _ { 1 } ^ { T }$ is symmetric suggesting that we consider the following factorization:

$$
P _ {1} A P _ {1} ^ {T} = \left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{c c} \alpha & 0 \\ 0 & \tilde {A} \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & I _ {n - 1} \end{array} \right] ^ {T}
$$

where

$$
\tilde {A} = B - \frac {1}{\alpha} v v ^ {T}.
$$

Note that with this kind of symmetric pivoting, the new (1,1) entry α is some diagonal entry $a _ { i i }$ . Our plan is to choose $P _ { 1 }$ so that α is the largest of A’s diagonal entries. If we apply the same strategy recursively to A˜ and compute

$$
\tilde {P} \tilde {A} \tilde {P} ^ {T} = \tilde {L} \tilde {D} \tilde {L} ^ {T},
$$

then we emerge with the factorization

$$
P A P ^ {T} = L D L ^ {T} \tag {4.2.10}
$$

where

$$
P = \left[ \begin{array}{c c} 1 & 0 \\ 0 & \tilde {P} \end{array} \right] P _ {1}, \qquad L = \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & \tilde {L} \end{array} \right], \qquad D = \left[ \begin{array}{c c} \alpha & 0 \\ 0 & \tilde {D} \end{array} \right].
$$

By virtue of this pivot strategy,

$$
d _ {1} \geq d _ {2} \geq \dots \geq d _ {n} > 0.
$$

Here is a nonrecursive implementation of the overall algorithm:

Algorithm 4.2.2 (Outer Product $\mathbf { \mathbf { \mathbf { \mathbf { L } } } } \mathbf { \mathbf { \mathbf { D } } } \mathbf { \mathbf \mathbf { \mathbf { \mathbf { L } } } } ^ { T }$ with Pivoting) Given a symmetric positive semidefinite $A \in \mathbb { R } ^ { n \times n }$ , the following algorithm computes a permutation P , a unit lower triangular L, and a diagonal matrix $D = \mathrm { d i a g } ( d _ { 1 } , \ldots , d _ { n } )$ so $P A P ^ { T } = L D L ^ { T }$ with $d _ { 1 } \geq d _ { 2 } \geq \dots \geq d _ { n } > 0 .$ The matrix element $a _ { i j }$ is overwritten by $d _ { i }$ if $i = j$ and by $\ell _ { i j } { \mathrm { ~ i f ~ } } i > j . P = P _ { 1 } \cdot \cdot \cdot P _ { n } $ where $P _ { k }$ is the identity with rows k and $p i v ( k )$ interchanged.

$$
\begin{array}{l} \text { for } k = 1: n \\ p i v (k) = j \text { where } a _ {j j} = \max \{a _ {k k}, \dots , a _ {n n} \} \\ A (k,:) \leftrightarrow A (j,:) \\ A (:, k) \leftrightarrow A (:, j) \\ \alpha = A (k, k) \\ v = A (k + 1: n, k) \\ A (k + 1: n, k) = v / \alpha \\ A (k + 1: n, k + 1: n) = A (k + 1: n, k + 1: n) - v v ^ {T} / \alpha \\ \end{array}
$$

end

If symmetry is exploited in the outer product update, then $n ^ { 3 } / 3$ flops are required. To solve $A x = b$ given $P A P ^ { T } = L D L ^ { T }$ , we proceed as follows:

$$
L w = P b, \qquad D y = w, \qquad L ^ {T} z = y, \qquad x = P ^ {T} z.
$$

We mention that Algorithm 4.2.2 can be implemented in a way that only references the lower trianglar part of A.

It is reasonable to ask why we even bother with the $\mathrm { L D L } ^ { T }$ factorization given that it appears to offer no real advantage over the Cholesky factorization. There are two reasons. First, it is more efficient in narrow band situations because it avoids square roots; see §4.3.6. Second, it is a graceful way to introduce factorizations of the form

$$
P A P ^ {T} = \left( \begin{array}{c} \mathrm{lower} \\ \mathrm{triangular} \end{array} \right) \times \left( \begin{array}{c} \mathrm{simple} \\ \mathrm{matrix} \end{array} \right) \times \left( \begin{array}{c} \mathrm{lower} \\ \mathrm{triangular} \end{array} \right) ^ {T},
$$

where P is a permutation arising from a symmetry-exploiting pivot strategy. The symmetric indefinite factorizations that we develop in §4.4 fall under this heading as does the “rank revealing” factorization that we are about to discuss for semidefinite problems.
