# 12.5.1 The Higher-Order SVD

Let us think about the SVD of $A \in \mathbb { R } ^ { m \times n }$ , not as

$$
A = U \Sigma V ^ {T} = \sum_ {i = 1} ^ {n} \sigma_ {i} u _ {i} v _ {i} T, \tag {12.5.1}
$$

but as $U ^ { T } A = \Sigma V ^ { T }$ . The matrix U structures the rows of $U ^ { T } A$ so that they are orthogonal to each other and monotone decreasing in norm:

$$
U ^ {T} A = \left[ \begin{array}{c} \sigma_ {1} v _ {1} ^ {T} \\ \vdots \\ \sigma_ {n} v _ {n} ^ {T} \end{array} \right]. \tag {12.5.2}
$$

The optimality of this structure can be seen by considering the following problem:

$$
\max \quad \| Q ^ {T} A \| _ {F}, \quad Q \in \mathbb {R} ^ {m \times r}. \tag {12.5.3}
$$

$$
Q ^ {T} Q = I _ {r}
$$

It is easy to verify that the maximum value is $\sigma _ { 1 } ^ { 2 } + \cdots + \sigma _ { r } ^ { 2 }$ and that it can be attained by setting $Q = U ( : , 1 : r )$ . The left singular vector matrix does the best job from the standpoint of getting as much “mass” as possible to the top of the transformed A. And that is what SVD does—it concentrates mass and supports an illuminating rank-1 expansion.

Now suppose $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times n _ { 2 } \times n _ { 3 } }$ and consider the following triplet of $\mathrm { S V D ^ { \prime } s } .$ one for each modal unfolding:

$$
U _ {1} ^ {T} \mathcal {A} _ {(1)} = \Sigma_ {1} V _ {1} ^ {T}, \quad U _ {2} ^ {T} \mathcal {A} _ {(2)} = \Sigma_ {2} V _ {2} ^ {T}, \quad U _ {3} ^ {T} \mathcal {A} _ {(3)} \Sigma_ {3} V _ {3} ^ {T}. \tag {12.5.4}
$$

These define three independent modal products:

$$
\mathcal {B} ^ {(1)} = \mathcal {A} \times_ {1} U _ {1}, \quad \mathcal {B} ^ {(2)} = \mathcal {A} \times_ {2} U _ {2}, \quad \mathcal {B} ^ {(3)} = \mathcal {A} \times_ {3} U _ {3}. \tag {12.5.5}
$$

Using Theorem 12.4.1, we have the following unfoldings:

$$
\mathcal {B} _ {(1)} ^ {(1)} = \Sigma_ {1} V _ {1} ^ {T} (U _ {3} \otimes U _ {2}) ^ {T}, \qquad \mathcal {B} _ {(2)} ^ {(2)} = \Sigma_ {2} V _ {2} ^ {T} (U _ {3} \otimes U _ {1}) ^ {T}, \qquad \mathcal {B} _ {(3)} ^ {(3)} = \Sigma_ {1} V _ {1} ^ {T} (U _ {2} \otimes U _ {1}) ^ {T}.
$$

Note that each of these matrices has the same kind singular value “grading” that is displayed in (12.5.1). Recalling from §12.4.5 that the rows of an unfolding are subtensors, it is easy to show that

$$
\left\| \mathcal {B} ^ {(1)} (i,:,:) \right\| _ {F} = \sigma_ {i} (\mathcal {A} _ {(1)}), \quad i = 1: n _ {1},
$$

$$
\| \mathcal {B} ^ {(2)} (:, i,:) \| _ {F} = \sigma_ {i} (\mathcal {A} _ {(2)}), \quad i = 1: n _ {2},
$$

$$
\left\| \mathcal {B} ^ {(3)} (:,:, i) \right\| _ {F} = \sigma_ {i} \left(\mathcal {A} _ {(3)}\right), \quad i = 1: n _ {3}.
$$

If we assemble these three modal products into a single multilinear product, then we get

$$
\mathcal {S} = \mathcal {A} \times_ {1} U _ {1} ^ {T} \times_ {2} U _ {2} ^ {T} \times_ {3} U _ {3} ^ {T}.
$$

Because the $U _ { i }$ are orthogonal, we can apply Theorem 12.4.1 and get

$$
\mathcal {A} = \mathcal {S} \times_ {1} U _ {1} \times_ {2} U _ {2} \times_ {3} U _ {3}.
$$

This is the higher-order SVD (HOSVD) developed by De Lathauwer, De Moor, and Vandewalle (2000). We summarize some of its important properties in the following theorem.

Theorem 12.5.1 (HOSVD). If $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ and

$$
\mathcal {A} _ {(k)} = U _ {k} \Sigma_ {k} V _ {k} ^ {T}, \quad k = 1: d,
$$

are the SVDs of its modal unfoldings, then its HOSVD is given by

$$
\mathcal {A} = \mathcal {S} \times_ {1} U _ {1} \times_ {2} U _ {2} \dots \times_ {d} U _ {d} \tag {12.5.6}
$$

where $\mathcal { S } = \mathcal { A } \times _ { 1 } U _ { 1 } ^ { T } \times _ { 2 } U _ { 2 } ^ { T } \cdot \cdot \cdot \times _ { d } U _ { d } ^ { T }$ . The formulation (12.5.6) is equivalent to

$$
\mathcal {A} = \sum_ {\mathbf {j} = \mathbf {1}} ^ {\mathbf {n}} \mathcal {S} (\mathbf {j}) \cdot U _ {1} (:, j _ {1}) \circ \dots \circ U _ {d} (:, j _ {d}), \tag {12.5.7}
$$

$$
\mathcal {A} (\mathbf {i}) = \sum_ {\mathbf {j} = 1} ^ {\mathbf {n}} \mathcal {S} (\mathbf {j}) \cdot U _ {1} (i _ {1}, j _ {1}) \dots U _ {d} (i _ {d}, j _ {d}), \tag {12.5.8}
$$

$$
\operatorname{vec} (\mathcal {A}) = (U _ {d} \otimes \dots \otimes U _ {1}) \cdot \operatorname{vec} (\mathcal {S}). \tag {12.5.9}
$$

Moreover,

$$
\left\| \mathcal {S} _ {(k)} (i,:) \right\| _ {F} = \sigma_ {i} \left(A _ {(k)}\right), \quad i = 1: \operatorname{rank} \left(A _ {(k)}\right) \tag {12.5.10}
$$

for $k = 1 { : } d .$

Proof. We leave the verification of (12.5.7)–(12.5.9) to the reader. To establish (12.5.10), note that

$$
\begin{array}{l} \mathcal {S} _ {(k)} = U _ {k} ^ {T} \mathcal {A} _ {(k)} \left(U _ {d} \otimes \dots \otimes U _ {k + 1} \otimes U _ {k - 1} \otimes \dots \otimes U _ {1}\right) \\ = \Sigma_ {k} V _ {k} ^ {T} \left(U _ {d} \otimes \dots \otimes U _ {k + 1} \otimes U _ {k - 1} \otimes \dots \otimes U _ {1}\right). \\ \end{array}
$$

It follows that the rows of $S _ { ( k ) }$ are mutually orthogonal and that the singular values of $\boldsymbol { \mathcal { A } } _ { ( \boldsymbol { k } ) }$ are the 2-norms of these rows.

In the HOSVD, the tensor S is called the core tensor. Note that it is not diagonal. However, the inequalities (12.5.10) tell us that, the values in $s$ tend to be smaller as “distance” from the $( 1 , 1 , \ldots , 1 )$ entry increases.
