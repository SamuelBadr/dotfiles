# 4.2.4 The Cholesky Factor is not a Square Root

A matrix $X \in \mathbb { R } ^ { n \times n }$ that satisfies $A = X ^ { 2 }$ is a square root of A. Note that if A symmetric, positive definite, and not diagonal, then its Cholesky factor is not a square root. However, if $A = G G ^ { T }$ and $\mathbf { \mathit { X } } = { \bar { U } } { \bar { \Sigma } } U ^ { T }$ where $G = U \Sigma V ^ { T }$ is the SVD, then

$$
X ^ {2} = (U \Sigma U ^ {T}) (U \Sigma U ^ {T}) = U \Sigma^ {2} U ^ {T} = (U \Sigma V ^ {T}) (U \Sigma V ^ {T}) ^ {T} = G G ^ {T} = A.
$$

Thus, a symmetric positive definite matrix A has a symmetric positive definite square root denoted by $A ^ { 1 \bar { / } 2 }$ . We have more to say about matrix square roots in §9.4.2.
