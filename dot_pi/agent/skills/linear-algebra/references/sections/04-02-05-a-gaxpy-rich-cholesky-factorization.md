# 4.2.5 A Gaxpy-Rich Cholesky Factorization

Our proof of the Cholesky factorization in Theorem 4.2.7 is constructive. However, we can develop a more effective procedure by comparing columns in $A = G G ^ { T }$ . If $A \in \mathbb { R } ^ { n \times n }$ and $1 \leq j \leq n$ , then

$$
A (:, j) = \sum_ {k = 1} ^ {j} G (j, k) \cdot G (:, k).
$$

This says that

$$
G (j, j) G (:, j) = A (:, j) - \sum_ {k = 1} ^ {j - 1} G (j, k) \cdot G (:, k) \equiv v. \tag {4.2.9}
$$

If the first $j - 1$ columns of G are known, then v is computable. It follows by equating components in (4.2.9) that

$$
G (j: n, j) = v (j: n) / \sqrt {v (j)}
$$

and so we obtain

for $j = 1:n$ $v(j:n) = A(j:n, j)$ for $k = 1:j - 1$ $v(j:n) = v(j:n) - G(j, k) \cdot G(j:n, k)$ end $G(j:n, j) = v(j:n) / \sqrt{v(j)}$ end

It is possible to arrange the computations so that G overwrites the lower triangle of A.

Algorithm 4.2.1 (Gaxpy Cholesky) Given a symmetric positive definite $A \in \mathbb { R } ^ { n \times n }$ , the following algorithm computes a lower triangular G such that $A = G G ^ { T }$ . For all $i \geq j , G ( i , j )$ overwrites $A ( i , j )$ .

for $j = 1:n$ if $j > 1$ $A(j:n,j) = A(j:n,j) - A(j:n,1:j - 1) \cdot A(j,1:j - 1)^T$ end $A(j:n,j) = A(j:n,j) / \sqrt{A(j,j)}$ end

This algorithm requires $n ^ { 3 } / 3$ flops.
