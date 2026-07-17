# 11.1.2 Operations and Allocations

Consider the gaxpy operation $y = y + A x$ with A in compressed-column format. If A ∈ IRm×n $A \in \mathbb { R } ^ { m \times n }$ and the dense vectors $\boldsymbol { y } \in \mathbb { R } ^ { m }$ and $\boldsymbol { x } \in \mathbb { R } ^ { n }$ are conventionally stored, then

for $j = 1 { : } n$

$$
k = A. c (j): A. c (j + 1) - 1 \tag {11.1.1}
$$

$$
y (A. r (k)) = y (A. r (k)) + A. v a l (k) \cdot x (j)
$$

end

overwrites y with $y + A x$ . It is easy to show that $2 \cdot \mathsf { n n z } ( A )$ flops are required. Regarding memory access, x is referenced sequentially, y is referenced randomly, and A is referenced through A.r and A.c.

A second example highlights the issue of memory allocation. Consider the outerproduct update $\overset { \vartriangle } { \boldsymbol { A } } = \overset { \vartriangle } { \boldsymbol { A } } + \overset { \vartriangle } { \boldsymbol { u } } \boldsymbol { v } ^ { T }$ where $A \in \mathbb { R } ^ { m \times n } , u \in \mathbb { R } ^ { m }$ , and $v \in \mathbb { R } ^ { n }$ are each stored in compressed-column format. In general, the updated A will have more nonzeros than the original A, e.g.,

![](images/golub_600_649__dc27d5072a710cc949d629cde5f1ad6bab49845024d2bffb5404a0e191b4d796.jpg)

<details>
<summary>text_image</summary>

Mathematical diagram showing a grid with black dots and an equals sign, followed by a plus sign and a separate grid with black dots.
</details>

Thus, unlike dense matrix computations where we simply overwrite A with $A + u v ^ { T }$ without concern for additional storage, now we must increase the memory allocation for A in order to house the result. Moreover, the expansion of the vectors A.val and A.r to accommodate the new nonzero entries is a nontrivial overhead. On the other hand, if we can predict the sparsity structure of $A + u v ^ { T }$ in advance and allocate space accordingly, then the update can be carried out more efficiently. This amounts to storing zeros in locations that are destined to become nonzero, e.g.,

$$
A. v a l = \begin{array}{c c c c c c c c c c c c c c c c} \hline a _ {1 1} & a _ {4 1} & 0 & a _ {5 2} & a _ {2 3} & a _ {3 3} & a _ {6 3} & a _ {1 4} & 0 & a _ {4 4} & 0 & a _ {2 5} & a _ {5 5} & a _ {6 5} \\ \hline \end{array} ,
$$

$$
A. c = \quad \boxed {1} \quad 3 \quad 5 \quad 8 \quad 1 2 \quad 1 5 \quad ,
$$

$$
A. r = \quad \boxed {1} \quad 4 \quad 3 \quad 5 \quad 2 \quad 3 \quad 6 \quad 1 \quad 3 \quad 4 \quad 5 \quad 2 \quad 5 \quad 6 \quad .
$$

With this assumption, the outer product update can proceed as follows:

$$
\text { for } \beta = 1: \mathsf {n n z} (v)
$$

$$
j = v. r (\beta)
$$

$$
\alpha = 1
$$

$$
\text { for } \ell = A. c (j): A. c (j + 1) - 1
$$

$$
\text { if } \alpha \leq \mathsf {n n z} (u) \& \& A. r (\ell) = u. r (\alpha) \tag {11.1.2}
$$

$$
A. v a l (\ell) = A. v a l (\ell) + u. v a l (\alpha) \cdot v. v a l (\beta)
$$

$$
\alpha = \alpha + 1
$$

end end end

Note that A.val() houses $a _ { i j }$ and is updated only if $u _ { i } v _ { j }$ is nonzero. The index α is used to reference the nonzero entries of u and is incremented after every access.

The overall success of a sparse matrix procedure typically depends strongly upon how efficiently it predicts and manages the fill-in phenomenon.
