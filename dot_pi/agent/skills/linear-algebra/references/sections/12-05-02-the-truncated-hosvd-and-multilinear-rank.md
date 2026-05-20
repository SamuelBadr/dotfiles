# 12.5.2 The Truncated HOSVD and Multilinear Rank

If $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ , then its multilinear rank is a the vector of modal unfolding ranks:

$$
\operatorname{rank} _ {*} (\mathcal {A}) = \left[ \operatorname{rank} \left(\mathcal {A} _ {(1)}\right), \dots , \operatorname{rank} \left(\mathcal {A} _ {(d)}\right) \right].
$$

Note that the summation upper bounds in the HOSVD can be replaced by rank (A). For example, (12.5.7) becomes

$$
\mathcal {A} = \sum_ {\mathbf {j} = \mathbf {1}} ^ {\operatorname{rank} _ {*} (\mathcal {A})} \mathcal {S} (\mathbf {j}) U _ {1} (:, j _ {1}) \circ \dots \circ U _ {d} (:, j _ {d}).
$$

This suggests a path to low-rank approximation. If $\mathbf { r } \le$ rank (A) with inquality in at least one component, then we can regard

$$
\mathcal {A} ^ {(\mathbf {r})} = \sum_ {\mathbf {j} = \mathbf {1}} ^ {\mathbf {r}} \mathcal {S} (\mathbf {j}) U _ {1} (:, j _ {1}) \circ \dots \circ U _ {d} (:, j _ {d})
$$

as a truncated HOSVD approximation to A. It can be shown that

$$
\left\| \mathcal {A} - \mathcal {A} ^ {(\mathbf {r})} \right\| _ {F} ^ {2} \leq \min _ {1 \leq k \leq d} \sum_ {i = r _ {k} + 1} ^ {\operatorname{rank} \left(\mathcal {A} _ {(k)}\right)} \sigma_ {i} \left(A _ {(k)}\right) ^ {2}. \tag {12.5.11}
$$
