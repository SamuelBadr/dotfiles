# 4.2.8 The Symmetric Semidefinite Case

A symmetric matrix $A \in \mathbb { R } ^ { n \times n }$ is positive semidefinite if

$$
x ^ {T} A x \geq 0
$$

for every $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . It is easy to show that if $A \in \mathbb { R } ^ { n \times n }$ is symmetric and positive semidefinite, then its eigenvalues satisfy

$$
0 = \lambda_ {n} (A) = \dots = \lambda_ {r + 1} (A) <   \lambda_ {r} (A) \leq \dots \leq \lambda_ {1} (A) \tag {4.2.11}
$$

where r is the rank of A. Our goal is to show that Algorithm 4.2.2 can be used to estimate r and produce a streamlined version of (4.2.10). But first we establish some useful properties.

Theorem 4.2.8. If $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive semidefinite, then

$$
\left| a _ {i j} \right| \leq \left(a _ {i i} + a _ {j j}\right) / 2, \tag {4.2.12}
$$

$$
\left| a _ {i j} \right| \leq \sqrt {a _ {i i} a _ {j j}}, \quad (i \neq j), \tag {4.2.13}
$$

$$
\max \left| a _ {i j} \right| = \max a _ {i i}, \tag {4.2.14}
$$

$$
a _ {i i} = 0 \Rightarrow A (i,:) = 0, A (:, i) = 0. \tag {4.2.15}
$$

Proof. Let $e _ { i }$ denote the ith column of $I _ { n }$ . Since

$$
x = e _ {i} + e _ {j} \Rightarrow 0 \leq x ^ {T} A x = a _ {i i} + 2 a _ {i j} + a _ {j j},
$$

$$
x = e _ {i} - e _ {j} \Rightarrow 0 \leq x ^ {T} A x = a _ {i i} - 2 a _ {i j} + a _ {j j},
$$

it follows that

$$
- 2 a _ {i j} \leq a _ {i i} + a _ {j j},
$$

$$
2 a _ {i j} \leq a _ {i i} + a _ {j j}.
$$

These two equations confirm (4.2.12), which in turn implies (4.2.14).

To prove (4.2.13), set $x = \tau e _ { i } + e _ { j }$ where $\tau \in \mathbb { R }$ . It follows that

$$
0 <   x ^ {T} A x = a _ {i i} \tau^ {2} + 2 a _ {i j} \tau + a _ {j j}
$$

must hold for all τ . This is a quadratic equation in τ and for the inequality to hold, the discriminant $4 a _ { i j } ^ { 2 } - 4 a _ { i i } a _ { j j }$ must be negative, i.e., $| a _ { i j } | \leq \sqrt { a _ { i i } a _ { j j } }$ . The implication in (4.2.15) follows immediately from (4.2.13).

Let us examine what happens when Algorithm 4.2.2 is applied to a rank-r positive semidefinite matrix. If $k \leq r$ , then after k steps we have the factorization

$$
\tilde {P} A \tilde {P} ^ {T} = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & I _ {n - k} \end{array} \right] \left[ \begin{array}{c c} D _ {k} & 0 \\ 0 & A _ {k} \end{array} \right] \left[ \begin{array}{c c} L _ {1 1} ^ {T} & L _ {2 1} ^ {T} \\ 0 & I _ {n - k} \end{array} \right] \tag {4.2.16}
$$

where $D _ { k } = \mathrm { d i a g } ( d _ { 1 } , \ldots , d _ { k } ) \in \mathbb { R } ^ { k \times k }$ and $d _ { 1 } \geq \dots \geq d _ { k } \geq 0$ . By virtue of the pivot strategy, if $d _ { k } = 0$ , then $A _ { k }$ has a zero diagonal. Since $A _ { k }$ is positive semidefinite, it follows from (4.2.15) that $A _ { k } = 0$ . This contradicts the assumption that A has rank r unless $k = r$ . Thus, if $k \leq r .$ , then $d _ { k } > 0$ . Moreover, we must have $A _ { r } = 0$ since A has the same rank as diag $( D _ { r } , A _ { r } )$ . It follows from (4.2.16) that

$$
P A P ^ {T} = \left[ \begin{array}{l} L _ {1 1} \\ L _ {2 1} \end{array} \right] D _ {r} \left[ \begin{array}{l} L _ {1 1} ^ {T} \mid L _ {2 1} ^ {T} \end{array} \right] \tag {4.2.17}
$$

where $D _ { r } = \operatorname { d i a g } ( d _ { 1 } , \ldots , d _ { r } )$ has positive diagonal entries, $L _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ is unit lower triangular, and $L _ { 2 1 } \in \mathbb { R } ^ { ( n - r ) \times r }$ . If $\ell _ { j }$ is the jth column of the L-matrix, then we can rewrite (4.2.17) as a sum of rank-1 matrices:

$$
P A P ^ {T} = \sum_ {j = 1} ^ {r} d _ {j} \ell_ {j} \ell_ {j} ^ {T}.
$$

This can be regarded as a relatively cheap alternative to the SVD rank-1 expansion.

It is important to note that our entire semidefinite discussion has been an exact arithmetic discussion. In practice, a threshold tolerance for small diagonal entries has to be built into Algorithm 4.2.2. If the diagonal of the computed $A _ { k }$ in (4.2.16) is sufficiently small, then the loop can be terminated and $\tilde { r }$ can be regarded as the numerical rank of A. For more details, see Higham (1989).
