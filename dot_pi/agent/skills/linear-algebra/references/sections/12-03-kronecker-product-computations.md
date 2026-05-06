# 12.3 Kronecker Product Computations

The Kronecker product (KP) has a rich algebra that supports a wide range of fast, practical algorithms. It also provides a bridge between matrix computations and tensor computations. This section is a compendium of its most important properties from that point of view. Recall that we introduced the KP in §1.3.6 and identified a few of its properties in §1.3.7 and §1.3.8. Our discussion of fast transforms in §1.4 and the 2-dimensional Poisson problem in §4.8.4 made heavy use of the operation.

# 12.3.1 Basic Properties

Kronecker product computations are structured block matrix computations. Basic properties are given in §1.3.6–§1.3.8, including

$$
\text { Transpose: } \quad (B \otimes C) ^ {T} \quad = B ^ {T} \otimes C ^ {T} ,
$$

$$
\text { Inverse: } \quad (B \otimes C) ^ {- 1} \quad = B ^ {- 1} \otimes C ^ {- 1},
$$

$$
\text { Product: } \quad (B \otimes C) (D \otimes F) = B D \otimes C F ,
$$

$$
\text { Associativity: } \quad B \otimes (C \otimes D) \quad = \quad (B \otimes C) \otimes D.
$$

Recall that $B \otimes C \neq C \otimes B$ , but if $B \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } }$ 1 and $C \in \mathbb { R } ^ { m _ { 2 } \times n _ { 2 } }$ , then

$$
P (B \otimes C) Q ^ {T} = C \otimes B \tag {12.3.1}
$$

where $P = \mathcal { P } _ { m _ { 1 } , m _ { 2 } }$ and $Q \ = \ \mathcal P _ { n _ { 1 } , n _ { 2 } }$ are perfect shuffle permutations, see §1.2.11.

Regarding the Kronecker product of structured matrices, if B is sparse, then $B \otimes C$ has the same sparsity pattern at the block level. If B and C are permutation matrices, then $B \otimes C$ is also a permutation matrix. Indeed, if p and q are permutations of 1:m and 1:n, then

$$
I _ {m} (p,:) \otimes I _ {n} (q,:) = I _ {m n} (w,:), \quad w = \left(1 _ {m} \otimes q\right) + n \cdot \left(p - 1 _ {m}\right) \otimes \mathbf {1} _ {n}. \tag {12.3.2}
$$

We also have

$$
(\text { orthogonal }) \otimes (\text { orthogonal }) = (\text { orthogonal }),
$$

$$
(\text { stochastic }) \otimes (\text { stochastic }) = (\text { stochastic }),
$$

$$
(\text { sym   pos   def }) \otimes (\text { sym   pos   def }) = (\text { sym   pos   def }).
$$

The inheritance of positive definiteness follows from

$$
\begin{array}{l} B = G _ {B} G _ {B} ^ {T} \\ C = G _ {C} G _ {C} ^ {T} \end{array} \Rightarrow B \otimes C = G _ {B} G _ {B} ^ {T} \otimes G _ {C} G _ {C} ^ {T} = (G _ {B} \otimes G _ {C}) (G _ {B} \otimes G _ {C}) ^ {T}.
$$

In other words, the Cholesky factor of $B \otimes C$ is the Kronecker product of the B and C Cholesky factors. Similar results apply to square LU and QR factorizations:

$$
\left. \begin{array}{l} P _ {B} B = L _ {B} U _ {B} \\ P _ {C} C = L _ {C} U _ {C} \end{array} \right\} \Rightarrow (P _ {B} \otimes P _ {C}) (B \otimes C) = (L _ {B} \otimes L _ {C}) (U _ {B} \otimes U _ {C}),
$$

$$
\left. \begin{array}{l} B = Q _ {B} R _ {B} \\ C = Q _ {C} R _ {C} \end{array} \right\} \Rightarrow B \otimes C = (Q _ {B} \otimes Q _ {C}) (R _ {B} \otimes R _ {C}).
$$

It should be noted that if B and/or C have more rows than columns, then the same can be said about the upper triangular matrices $R _ { B }$ and $R _ { C }$ . In this case, row permutations of $R _ { B } \otimes R _ { C }$ are required to achieve triangular form. On the other hand,

$$
(B \otimes C) (P _ {B} \otimes P _ {C}) = (Q _ {B} \otimes Q _ {C}) (R _ {B} \otimes R _ {C})
$$

is a thin QR factorization of $B \otimes C$ if $B P _ { B } = Q _ { B } R _ { B }$ and $C P _ { c } = Q _ { c } R _ { c }$ are thin QR factorizations.

The eigenvalues and singular values of $B \otimes C$ have a product connection to the eigenvalues and singular values of B and C:

$$
\lambda (B \otimes C) = \{\beta_ {i} \gamma_ {j}: \beta_ {i} \in \lambda (B), \gamma_ {j} \in \lambda (C) \},
$$

$$
\sigma (B \otimes C) = \{\beta_ {i} \gamma_ {j}: \beta_ {i} \in \sigma (B), \gamma_ {j} \in \sigma (C) \}.
$$

These results are a consequence of the following decompositions:

$$
\left. \begin{array}{l} Q _ {B} ^ {H} B Q _ {B} = T _ {B} \\ Q _ {C} ^ {H} C Q _ {C} = T _ {C} \end{array} \right\} \Rightarrow (Q _ {B} \otimes Q _ {C}) ^ {H} (B \otimes C) (Q _ {B} \otimes Q _ {C}) = T _ {B} \otimes T _ {C}, \tag {12.3.3}
$$

$$
\left. \begin{array}{l} U _ {B} ^ {H} B V _ {B} = \Sigma_ {B} \\ U _ {C} ^ {H} C V _ {C} = \Sigma_ {C} \end{array} \right\} \Rightarrow (U _ {B} \otimes U _ {C}) ^ {H} (B \otimes C) (V _ {B} \otimes V _ {C}) = \Sigma_ {B} \otimes \Sigma_ {C}. \tag {12.3.4}
$$

Note that if $B y = \beta y$ and $C z = \gamma z$ , then $( B \otimes C ) ( y \otimes z ) = \beta \gamma ( y \otimes z )$ . Other properties that follow from (12.3.3) and (12.3.4) include

$$
\operatorname{rank} (B \otimes C) = \operatorname{rank} (B) \cdot \operatorname{rank} (C),
$$

$$
\det (B \otimes C) = \det (B) ^ {n} \cdot \det (C) ^ {m}, \quad B \in \mathbb {R} ^ {m \times m}, C \in \mathbb {R} ^ {n \times n},
$$

$$
\operatorname{tr} (B \otimes C) = \operatorname{tr} (B) \cdot \operatorname{tr} (C),
$$

$$
\left\| B \otimes C \right\| _ {F} = \left\| B \right\| _ {F} \cdot \left\| C \right\| _ {F},
$$

$$
\| B \otimes C \| _ {2} = \| B \| _ {2} \cdot \| C \| _ {2}.
$$

See Horn and Johnson (TMA) for additional KP facts.

# 12.3.2 The Tracy-Singh Product

We can think of the Kronecker product of two matrices $B = \left( b _ { i j } \right)$ and $C = \left( c _ { i j } \right)$ as the systematic layout of all possible products $b _ { i j } c _ { k \ell } , \mathrm { e . g . }$ ,

$$
\left[ \begin{array}{c c} b _ {1 1} & b _ {1 2} \\ b _ {2 1} & b _ {2 2} \end{array} \right] \otimes \left[ \begin{array}{c c} c _ {1 1} & c _ {1 2} \\ c _ {2 1} & c _ {2 2} \end{array} \right] = \left[ \begin{array}{c c c c} b _ {1 1} c _ {1 1} & b _ {1 1} c _ {1 2} & b _ {1 2} c _ {1 1} & b _ {1 2} c _ {1 2} \\ b _ {1 1} c _ {2 1} & b _ {1 1} c _ {2 2} & b _ {1 2} c _ {2 1} & b _ {1 2} c _ {2 2} \\ \hline b _ {2 1} c _ {1 1} & b _ {2 1} c _ {1 2} & b _ {2 2} c _ {1 1} & b _ {2 2} c _ {1 2} \\ b _ {2 1} c _ {2 1} & b _ {2 1} c _ {2 2} & b _ {2 2} c _ {2 1} & b _ {2 2} c _ {2 2} \end{array} \right].
$$

However, the Kronecker product of two block matrices $B = \left( B _ { i j } \right)$ and $C _ { i j } )$ is not the corresponding layout of all possible block-level Kronecker products $B _ { i j } \otimes B _ { k \ell } \colon$

$$
\left[ \begin{array}{c c} B _ {1 1} & B _ {1 2} \\ B _ {2 1} & B _ {2 2} \end{array} \right] \otimes \left[ \begin{array}{c c} C _ {1 1} & C _ {1 2} \\ C _ {2 1} & C _ {2 2} \end{array} \right] \neq \left[ \begin{array}{c c c c} B _ {1 1} C _ {1 1} & B _ {1 1} C _ {1 2} & B _ {1 2} C _ {1 1} & B _ {1 2} C _ {1 2} \\ B _ {1 1} C _ {2 1} & B _ {1 1} C _ {2 2} & B _ {1 2} C _ {2 1} & B _ {1 2} C _ {2 2} \\ \hline B _ {2 1} C _ {1 1} & B _ {2 1} C _ {1 2} & B _ {2 2} C _ {1 1} & B _ {2 2} C _ {1 2} \\ B _ {2 1} C _ {2 1} & B _ {2 1} C _ {2 2} & B _ {2 2} C _ {2 1} & B _ {2 2} C _ {2 2} \end{array} \right].
$$

The matrix on the right is an example of the Tracy-Singh product. Formally, if we are given the blockings

$$
B = \left[ \begin{array}{c c c} B _ {1 1} & \dots & B _ {1, N _ {1}} \\ \vdots & \ddots & \vdots \\ B _ {M _ {1}, 1} & \dots & B _ {M _ {1}, N _ {1}} \end{array} \right] \quad C = \left[ \begin{array}{c c c} C _ {1 1} & \dots & C _ {1, N _ {2}} \\ \vdots & \ddots & \vdots \\ C _ {M _ {2}, 1} & \dots & C _ {M _ {2}, N _ {2}} \end{array} \right], \tag {12.3.5}
$$

with $B _ { i j } \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } }$ and $C _ { i j } \in \mathbb { R } ^ { m _ { 2 } \times n _ { 2 } }$ , then the Tracy-Singh product is an $M _ { 1 ^ { - } } \mathrm { b y } { - } N _ { 1 }$ block matrix $B \underset { \mathbf { r } \mathbf { s } } { \otimes } C$ whose $( i , j )$ block is given by

$$
[ B \underset {\mathbf {T S}} {\otimes} C ] _ {i j} = \left[ \begin{array}{c c c} B _ {i j} \otimes C _ {1 1} & \dots & B _ {i j} \otimes C _ {1, N _ {2}} \\ \vdots & \ddots & \vdots \\ B _ {i j} \otimes C _ {M _ {2}, 1} & \dots & B _ {i j} \otimes C _ {M _ {2}, N _ {2}} \end{array} \right].
$$

See Tracy and Singh (1972). Given (12.3.5), it can be shown using (12.3.1) that

$$
B \underset {\mathbf {T S}} {\otimes} C = P (B \otimes C) Q ^ {T} \tag {12.3.6}
$$

where

$$
P = \left(I _ {M _ {1} M _ {2}} \otimes \mathcal {P} _ {m _ {2}, m _ {1}}\right) \left(I _ {M _ {1}} \otimes \mathcal {P} _ {m _ {1}, M _ {2} m _ {2}}\right), \tag {12.3.7}
$$

$$
Q = \left(I _ {N _ {1} N _ {2}} \otimes \mathcal {P} _ {n _ {2}, n _ {1}}\right) \left(I _ {N _ {1}} \otimes \mathcal {P} _ {n _ {1}, N _ {2} n _ {2}}\right). \tag {12.3.8}
$$

# 12.3.3 The Hadamard and Khatri-Rao Products

There are two submatrices of $B \otimes C$ that are particularly important. The Hadamard Product is a pointwise product:

$$
B \underset {\text { HAD }} {\otimes} C = B. * C.
$$

Thus, if $B \in \mathbb { R } ^ { m \times n }$ and $C \in \mathbb { R } ^ { m \times n }$ , then

$$
\left[ \begin{array}{l l} b _ {1 1} & b _ {1 2} \\ b _ {2 1} & b _ {2 2} \\ b _ {3 1} & b _ {3 2} \end{array} \right] \underset {\mathbf {H A D}} {\otimes} \left[ \begin{array}{l l} c _ {1 1} & c _ {1 2} \\ c _ {2 1} & c _ {2 2} \\ c _ {3 1} & c _ {3 2} \end{array} \right] = \left[ \begin{array}{l l} b _ {1 1} c _ {1 1} & b _ {1 2} c _ {1 2} \\ b _ {2 1} c _ {2 1} & b _ {2 2} c _ {2 2} \\ b _ {3 1} c _ {3 1} & b _ {3 2} c _ {3 2} \end{array} \right].
$$

The block analog of this is the Khatri-Rao Product. If $B = \left( B _ { i j } \right)$ and $C = ( C _ { i j } )$ are each m-by-n block matrices, then

$$
B \underset {\mathbf {K R}} {\otimes} C = (A _ {i j}), \qquad A _ {i j} = B _ {i j} \otimes C _ {i j},
$$

e.g.,

$$
\left[ \begin{array}{c c} B _ {1 1} & B _ {1 2} \\ B _ {2 1} & B _ {2 2} \\ B _ {3 1} & B _ {3 2} \end{array} \right] \underset {\mathbf {K R}} {\otimes} \left[ \begin{array}{c c} C _ {1 1} & C _ {1 2} \\ C _ {2 1} & C _ {2 2} \\ C _ {3 1} & C _ {3 2} \end{array} \right] = \left[ \begin{array}{c c} B _ {1 1} \otimes C _ {1 1} & B _ {1 2} \otimes C _ {1 2} \\ B _ {2 1} \otimes C _ {2 1} & B _ {2 2} \otimes C _ {2 2} \\ B _ {3 1} \otimes C _ {3 1} & B _ {3 2} \otimes C _ {3 2} \end{array} \right].
$$

A particularly important instance of the Khatri-Rao product is based on column partitionings:

$$
\left[ \begin{array}{c c c c} b _ {1} & \dots & b _ {n} \end{array} \right] \underset {\mathbf {K R}} {\otimes} \left[ \begin{array}{c c c c} c _ {1} & \dots & c _ {n} \end{array} \right] = \left[ \begin{array}{c c c c} b _ {1} \otimes c _ {1} & \dots & b _ {n} \otimes c _ {n} \end{array} \right].
$$

For more details on the Khatri-Rao product, see Smilde, Bro, and Geladi (2004).

# 12.3.4 The Vec and Reshape Operations

In Kronecker product work, matrices are sometimes regarded as vectors and vectors are sometimes turned into matrices. To be precise about these reshapings, we remind the reader about the vec and reshape operations defined in §1.3.7. If $\boldsymbol { X } \in \mathbb { R } ^ { m \times n }$ , then vec(X) is an nm-by-1 vector obtained by “stacking” X’s columns:

$$
\operatorname{vec} (X) = \left[ \begin{array}{c} X (:, 1) \\ \vdots \\ X (:, n) \end{array} \right].
$$

If $B \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } }$ 1 , $C \in \mathbb { R } ^ { m _ { 2 } \times n _ { 2 } }$ , and $X \in \mathbb { R } ^ { n _ { 1 } \times m _ { 2 } }$ , then

$$
Y = C X B ^ {T} \Leftrightarrow \operatorname{vec} (Y) = (B \otimes C) \cdot \operatorname{vec} (X). \tag {12.3.9}
$$

Note that the matrix equation

$$
F _ {1} X G _ {1} ^ {T} + \dots + F _ {p} X G _ {p} ^ {T} = C \tag {12.3.10}
$$

is equivalent to

$$
\left(G _ {1} \otimes F _ {1} + \dots + G _ {p} \otimes F _ {p}\right) \operatorname{vec} (X) = \operatorname{vec} (C). \tag {12.3.11}
$$

See Lancaster (1970), Vetter (1975), and also our discussion about block diagonalization in §7.6.3.

The reshape operation takes a vector and turns it into a matrix. If $a \in \mathbb { R } ^ { m n }$ then

$$
A = \operatorname{reshape} (a, m, n) \in \mathbb {R} ^ {m \times n} \quad \Leftrightarrow \quad \operatorname{vec} (A) = a.
$$

Thus, if $u \in \mathbb { R } ^ { m }$ and $v \in \mathbb { R } ^ { n }$ , then reshape $( v \otimes u , m , n ) = u v ^ { T }$ .

# 12.3.5 Vec, Perfect Shuffles, and Transposition

There is an important connection between matrix transposition and perfect shuffle permutations. In particular, if $A \in \mathbb { R } ^ { q \times r }$ , then

$$
\mathsf {v e c} (A ^ {T}) = \mathcal {P} _ {r, q} \mathsf {v e c} (A). \tag {12.3.12}
$$

This formulation of matrix transposition provides a handy way to reason about large scale, multipass transposition algorithms that are required when $A \in \mathbb { R } ^ { q \times r }$ is too large to fit in fast memory. In this situation the transposition must proceed in stages and the overall process corresponds to a factorization of $\mathcal { P } _ { r , q }$ . For example, if

$$
\mathcal {P} _ {r, q} = \Gamma_ {t} \dots \Gamma_ {1} \tag {12.3.13}
$$

where each $\Gamma _ { k }$ is a “data-motion-friendly” permutation, then $B = A ^ { T }$ can be computed with t passes through the data:

$$
a = \operatorname{vec} (A)
$$

for k = 1:t

$$
a = \Gamma_ {k} a
$$

end

$$
B = \operatorname{reshape} (a, q, r)
$$

The idea is to choose a factorization (12.3.13) so that the data motion behind the operation kth pass, $\mathrm { i . e . , } a \gets \Gamma _ { k } a$ , is in harmony with the architecture of the underlying memory hierarchy, i.e., blocks that can fit in cache, etc.

As an illustration, suppose we want to assign $A ^ { T }$ to B where

$$
A = \left[ \begin{array}{c} A _ {1} \\ \vdots \\ A _ {r} \end{array} \right], \qquad A _ {k} \in \mathbb {R} ^ {q \times q}.
$$

We assume that A is stored by column which means that the $A _ { i }$ are not contiguous in memory. To complete the story, suppose each block comfortably fits in cache but that A cannot. Here is a 2-pass factorization of $\mathcal { P } _ { r q , q } \mathrm { : }$ :

$$
\mathcal {P} _ {q, r q} = \Gamma_ {2} \Gamma_ {1} = \left(I _ {r} \otimes \mathcal {P} _ {q, q}\right) \left(\mathcal {P} _ {r, q} \otimes I _ {q}\right).
$$

If ˜a = Γ1 · vec(A), then

$$
\operatorname{reshape} (\tilde {a}, q, r q) = \left[ \begin{array}{c c c c} A _ {1} & \dots & A _ {r} \end{array} \right].
$$

In other words, after the first pass through the data we have computed the block transpose of A. (The $A _ { i }$ are now contiguous in memory.) To complete the overall task, we must transpose each of these blocks. If $b = \Gamma _ { 2 } \tilde { a }$ , then

$$
B = \operatorname{reshape} (b, q, r q) = \left[ \begin{array}{c c c c} A _ {1} ^ {T} & \dots & A _ {r} ^ {T} \end{array} \right].
$$

See Van Loan (FFT) for more details about perfect shuffle factorizations and multipass matrix transposition algorithms.

# 12.3.6 The Kronecker Product SVD

Suppose $A \in \mathbb { R } ^ { m \times n }$ is given with $m = m _ { 1 } m _ { 2 }$ and $n = n _ { 1 } n _ { 2 }$ . For these integer factorizations the nearest Kronecker product (NKP) problem involves minimizing

$$
\phi (B, C) = \left\| A - B \otimes C \right\| _ {F} \tag {12.3.14}
$$

where $B \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } }$ and $C \in \mathbb { R } ^ { m _ { 2 } \times n _ { 2 } }$ . Van Loan and Pitsianis (1992) show how to solve the NKP problem using the singular value decomposition of a permuted version of A. A small example communicates the main idea. Suppose $m _ { 1 } = 3$ and $n _ { 1 } = m _ { 2 } = n _ { 2 } = 2$ . By carefully thinking about the sum of squares that define φ, we see that

$$
\phi (B, C) = \left\| \left[ \begin{array}{c c c c} a _ {1 1} & a _ {1 2} & a _ {1 3} & a _ {1 4} \\ a _ {2 1} & a _ {2 2} & a _ {2 3} & a _ {2 4} \\ \hline a _ {3 1} & a _ {3 2} & a _ {3 3} & a _ {3 4} \\ a _ {4 1} & a _ {4 2} & a _ {4 3} & a _ {4 4} \\ \hline a _ {5 1} & a _ {5 2} & a _ {5 3} & a _ {5 4} \\ a _ {6 1} & a _ {6 2} & a _ {6 3} & a _ {6 4} \end{array} \right] - \left[ \begin{array}{c c} b _ {1 1} & b _ {1 2} \\ b _ {2 1} & b _ {2 2} \\ b _ {3 1} & b _ {3 2} \end{array} \right] \otimes \left[ \begin{array}{c c} c _ {1 1} & c _ {1 2} \\ c _ {2 1} & c _ {2 2} \end{array} \right] \right\| _ {F}
$$

$$
= \left\| \left[ \begin{array}{l l l l} a _ {1 1} & a _ {2 1} & a _ {1 2} & a _ {2 2} \\ a _ {3 1} & a _ {4 1} & a _ {3 2} & a _ {4 2} \\ a _ {5 1} & a _ {6 1} & a _ {5 2} & a _ {6 2} \\ a _ {1 3} & a _ {2 3} & a _ {1 4} & a _ {2 4} \\ a _ {3 3} & a _ {4 3} & a _ {3 4} & a _ {4 4} \\ a _ {5 3} & a _ {6 3} & a _ {5 4} & a _ {6 4} \end{array} \right] - \left[ \begin{array}{l} b _ {1 1} \\ b _ {2 1} \\ b _ {3 1} \\ b _ {1 2} \\ b _ {2 2} \\ b _ {3 2} \end{array} \right] \left[ \begin{array}{l l l l} c _ {1 1} & c _ {2 1} & c _ {1 2} & c _ {2 2} \end{array} \right] \right\| _ {F}.
$$

Denote the preceding 6-by-4 matrix by R(A) and observe that

$$
\mathcal {R} (A) = \left[ \begin{array}{c} \mathsf {v e c} (A _ {1 1}) ^ {T} \\ \mathsf {v e c} (A _ {2 1}) ^ {T} \\ \mathsf {v e c} (A _ {3 1}) ^ {T} \\ \mathsf {v e c} (A _ {1 2}) ^ {T} \\ \mathsf {v e c} (A _ {2 2}) ^ {T} \\ \mathsf {v e c} (A _ {3 2}) ^ {T} \end{array} \right].
$$

It follows that

$$
\phi (B, C) = \left\| \mathcal {R} (A) - \operatorname{vec} (B) \operatorname{vec} (C) ^ {T} \right\| _ {F}
$$

and so the act of minimizing $\phi$ is equivalent to finding a nearest rank-1 matrix to $\mathcal { R } ( A )$ . This problem has a simple SVD solution. Referring to Theorem 2.4.8, if

$$
U ^ {T} \mathcal {R} (A) V = \Sigma \tag {12.3.15}
$$

is the SVD of $\mathcal { R } ( A )$ , then the optimizing B and C are defined by

$$
\mathsf {v e c} (B _ {\mathrm{opt}}) = \sqrt {\sigma_ {1}}   U (:, 1), \qquad \mathsf {v e c} (C _ {\mathrm{opt}}) = \sqrt {\sigma_ {1}}   V (:, 1).
$$

The scalings are arbitrary. Indeed, if $B _ { \mathrm { o p t } }$ and $C _ { \mathrm { o p t } }$ solve the NKP problem and $\alpha \neq 0$ then $\alpha \cdot B _ { \mathrm { o p t } }$ and $( 1 / \alpha ) \cdot C _ { \mathrm { o p t } }$ are also optimal.

In general, if

$$
A = \left[ \begin{array}{c c c} A _ {1 1} & \dots & A _ {1, n _ {1}} \\ \vdots & \ddots & \vdots \\ A _ {m _ {1}, 1} & \dots & A _ {m _ {1}, n _ {1}} \end{array} \right] \tag {12.3.16}
$$

where each $A _ { i j } \in \mathbb { R } ^ { m _ { 2 } \times n _ { 2 } }$ , then $\mathcal { R } ( A ) \in \mathbb { R } ^ { m _ { 1 } n _ { 1 } \times m _ { 2 } n _ { 2 } }$ is defined by

$$
\mathcal {R} (A) = \left[ \begin{array}{c} \tilde {A} _ {1} \\ \vdots \\ \tilde {A} _ {n _ {1}} \end{array} \right], \qquad \tilde {A} _ {j} = \left[ \begin{array}{c} \operatorname{vec} (A _ {1 j}) ^ {T} \\ \vdots \\ \operatorname{vec} (A _ {m _ {1}, j}) ^ {T} \end{array} \right].
$$

The SVD of $\mathcal { R } ( A )$ can be “reshaped” into a special SVD-like expansion for A.

Theorem 12.3.1 (Kronecker Product SVD). If $A \in \mathbb { R } ^ { m _ { 1 } m _ { 2 } \times n _ { 1 } n _ { 2 } }$ is blocked according to (12.3.16) and

$$
\mathcal {R} (A) = U \Sigma V ^ {T} = \sum_ {k = 1} ^ {r} \sigma_ {k} \cdot u _ {k} v _ {k} ^ {T} \tag {12.3.17}
$$

is the SVD of R(A) with $u _ { k } = U ( : , k ) , v _ { k } = V ( : , k )$ , and $\sigma _ { k } = \Sigma ( k , k )$ , then

$$
A = \sum_ {k = 1} ^ {r} \sigma_ {k} \cdot U _ {k} \otimes V _ {k} \tag {12.3.18}
$$

where $U _ { k } ~ = ~ { \mathsf { r e s h a p e } } ( u _ { k } , m _ { 1 } , n _ { 1 } )$ and $V _ { k } = { \mathsf { r e s h a p e } } ( v _ { k } , m _ { 2 } , n _ { 2 } )$ .

Proof. In light of (12.3.18), we must show that

$$
A _ {i j} = \sum_ {k = 1} ^ {r} \sigma_ {k} \cdot U _ {k} (i, j) \cdot V _ {k}.
$$

But this follows immediately from (12.3.17) which says that

$$
\mathsf {v e c} (A _ {i j}) ^ {T} = \sum_ {k = 1} ^ {r} \sigma_ {k} \cdot U _ {k} (i, j) v _ {k} ^ {T}
$$

for all i and j.

The integer r in the theorem is the Kronecker product rank of A given the blocking (12.3.16). Note that if $\tilde { r } \leq r$ , then

$$
A _ {\tilde {r}} = \sum_ {k = 1} ^ {\tilde {r}} \sigma_ {k} U _ {k} \otimes V _ {k} \tag {12.3.19}
$$

is the closest matrix to A (in the Frobenius norm) that is the sum of ˜r Kronecker products. If A is large and sparse and ˜r is small, then the Lanzcos SVD iteration can effectively be used to compute the required singular values and vectors of $\mathcal { R } ( A )$ . See §10.4.

# 12.3.7 Constrained NKP Problems

If A is structured, then it is sometimes the case that the B and C matrices that solve the NKP problem are similarly structured. For example, if A is symmetric and positive definite, then the same can be said of $B _ { \mathrm { o p t } }$ and $C _ { \mathrm { o p t } }$ (if properly normalized). Likewise, if A is nonnegative, then the optimal B and C can be chosen to be nonnegative. These and other structured NKP problems are discussed in Van Loan and Pitsianis (1992).

We mention that a problem like

$$
\min _ {B, C \text { Toeplitz }} \| A - B \otimes C \| _ {F}, \qquad B \in \mathbb {R} ^ {m \times m},   C \in \mathbb {R} ^ {n \times n},
$$

turns into a constrained nearest rank-1 problem of the form

$$
\begin{array}{l} \min \quad \left\| \mathcal {A} - b c ^ {T} \right\| _ {F} \\ F ^ {T} \operatorname{vec} (B) = 0 \\ G ^ {T} \mathsf {v e c} (C) = 0 \\ \end{array}
$$

where the nullspaces of $F ^ { T }$ and $G ^ { T }$ define the vector space of m-by-m and n-by-n Toeplitz matrices respectively. This problem can be solved by computing QR factorizations of $F$ and G followed by a reduced-dimension SVD.

# 12.3.8 Computing the Nearest $X \otimes X$

Suppose A ∈ IRm2×m2 $A \in \mathbb { R } ^ { m ^ { 2 } \times m ^ { 2 } }$ and that we want to find $\boldsymbol { X } \in \mathbb { R } ^ { m \times m }$ so that

$$
\phi_ {\mathrm{sym}} (X) = \left\| A - X \otimes X \right\| _ {F}
$$

is minimized. Proceeding as we did with the NKP problem, we can reshape this into a nearest symmetric rank-1 problem:

$$
\phi_ {\text { sym }} (X) = \| \mathcal {R} (A) - \text { vec } (X) \cdot \text { vec } (X) ^ {T} \| _ {F}. \tag {12.3.20}
$$

It turns out that the solution $X _ { \mathrm { o p t } }$ is a reshaping of an eigenvector associated with the symmetric part of $\mathcal { R } ( A )$ .

Lemma 12.3.2. Suppose $M \in \mathbb { R } ^ { n \times n }$ and that $Q ^ { T } T Q = \operatorname { d i a g } ( \alpha _ { 1 } , . . . , \alpha _ { n } )$ is a Schur decomposition of $T = ( M + M ^ { T } ) / 2$ . If

$$
| \alpha_ {k} | = \max \{| \alpha_ {1} |, \dots , | \alpha_ {n} | \}
$$

then the solution to the problem

$$
\begin{array}{l} \min \quad \| M - Z \| _ {F} \\ Z = Z ^ {T} \\ \operatorname{rank} (Z) = 1 \\ \end{array}
$$

is given by $Z _ { \mathrm { o p t } } = \alpha _ { k } q _ { k } q _ { k } ^ { T }$ where $q _ { k } = Q ( : , k )$ .

Proof. See P12.3.11.

# 12.3.9 Computing the Nearest $X \otimes Y - Y \otimes X$

Suppose $A \in \mathbb { R } ^ { n \times n } , n = m ^ { 2 }$ and that we wish to find X, $Y \in \mathbb { R } ^ { m \times m }$ so that

$$
\phi_ {\text { skew }} (X, Y) = \left\| A - \left(X \otimes Y - Y \otimes X\right) \right\| _ {F}
$$

is minimized. It can be shown that

$$
\phi_ {\mathrm{skew}} (X) = \left\| \mathcal {R} (A) - (\mathrm{vec} (X) \cdot \mathrm{vec} (Y) ^ {T} - \mathrm{vec} (Y) \cdot \mathrm{vec} (X) ^ {T} \right\| _ {F}. \tag {12.3.21}
$$

The optimizing X and Y can be determined by exploiting the following lemma.

Lemma 12.3.3. Suppose $M \in \mathbb { R } ^ { n \times n }$ with skew-symmetric part $S = ( M - M ^ { T } ) / 2$ . If

$$
S [ u \mid v ] = [ u \mid v ] \left[ \begin{array}{c c} 0 & \mu \\ - \mu & 0 \end{array} \right], \qquad u, v \in \mathbb {R} ^ {n},
$$

with $\mu = \rho ( S ) , \| u \| _ { 2 } = \| v \| _ { 2 } = 1$ , and $u ^ { T } v = 0$ , then $Z _ { \mathrm { o p t } } = \mu \left( u v ^ { T } - v u ^ { T } \right)$ minimizes $\| M - Z \| _ { F }$ over all rank-2 skew-symmetric matrices $Z \in \mathbb { R } ^ { n \times n }$ .

Proof. See P12.3.12.

# 12.3.10 Some Comments About Multiple Kronecker Products

The Kronecker product of three or more matrices results in a matrix that has a recursive block structure. For example,

$$
B \otimes C \otimes D   =   \left[ \begin{array}{l l} b _ {1 1} & b _ {1 2} \\ b _ {2 1} & b _ {2 2} \end{array} \right] \otimes \left[ \begin{array}{l l l l} c _ {1 1} & c _ {1 2} & c _ {1 3} & c _ {1 4} \\ c _ {2 1} & c _ {2 2} & c _ {2 3} & c _ {2 4} \\ c _ {3 1} & c _ {3 2} & c _ {3 3} & c _ {3 4} \\ c _ {4 1} & c _ {4 2} & c _ {4 3} & c _ {4 4} \end{array} \right] \otimes \left[ \begin{array}{l l l} d _ {1 1} & d _ {1 2} & d _ {1 3} \\ d _ {2 1} & d _ {2 2} & d _ {2 3} \\ d _ {3 1} & d _ {3 2} & d _ {3 3} \end{array} \right]
$$

is a 2-by-2 block matrix whose entries are 4-by-4 block matrices whose entries are 3-by-3 matrices.

A Kronecker product can be regarded as a data-sparse representation. If $A =$ $B _ { 1 } \otimes B _ { 2 }$ and each B-matrix is m-by-m, then $2 m ^ { 2 }$ numbers are used to encode a matrix that has $m ^ { 4 }$ entries. The data sparsity is more dramatic for multiple Kronecker products. If $A = B _ { 1 } \otimes \dots \otimes B _ { p }$ and $B _ { i } \in \mathbb { R } ^ { m \times m }$ , then $p m ^ { 2 }$ numbers fully describe A, a matrix with $m ^ { 2 p }$ entries.

Order of operation can be important when a multiple Kronecker product is involved and the participating matrices vary in dimension. Suppose $B _ { i } \in \mathbb { R } ^ { m _ { i } \times n _ { i } }$ for $i = 1 { : } p$ and that $M _ { i } = m _ { 1 } \cdot \cdot \cdot m _ { i }$ and $N _ { i } = n _ { 1 } \cdot \cdot \cdot n _ { i }$ for $i = 1 { : } p .$ The matrix-vector product

$$
y = (B _ {1} \otimes \dots B _ {p}) x \qquad x \in \mathbb {R} ^ {N _ {p}}
$$

can be evaluated in many different orders and the associated flop counts can vary tremendously. The search for an optimal ordering is a dynamic programming problem that involves the recursive analysis of calculations like

$$
\operatorname{reshape} \left(y, M _ {p} / M _ {i}, M _ {i}\right) = \left(B _ {i + 1} \otimes \dots \otimes B _ {p}\right) \cdot \operatorname{reshape} \left(x, N _ {p} / N _ {i}, N _ {i}\right) \cdot \left(B _ {1} \otimes \dots B _ {i}\right) ^ {T}.
$$

# Problems

P12.3.1 Prove (12.3.1) and (12.3.2).

P12.3.2 Assume that the matrices $A _ { 1 } , \dots , A _ { N } \in \mathbb { R } ^ { m \times n }$ . Express the summation

$$
f (x, y) = \sum_ {k = 1} ^ {N} \left(y ^ {T} A _ {k} x - b _ {k}\right) ^ {2}
$$

in matrix-vector terms given that $y \in \mathbb { R } ^ { m } , x \in \mathbb { R } ^ { m }$ , and $b \in \mathbb { R } ^ { N }$ .

P12.3.3 A total least squares solution to $( B \otimes C ) x \approx b$ requires the computation of the smallest singular value and the associated right singular vector of the augmented matrix $M = [ B \otimes C | b ]$ . Outline an efficient procedure for doing this that exploits the Kronecker structure of the data matrix.

P12.3.4 Show how to minimize $\parallel ( A _ { 1 } \otimes A _ { 2 } ) x - f \parallel$ subject to the constraint that $( B _ { 1 } \otimes B _ { 2 } ) x = g .$ . Assume that $A _ { 1 }$ and $A _ { 2 }$ have more rows than columns and that $B _ { 1 }$ and $B _ { 2 }$ have more columns than rows. Also assume that each of these four matrices has full rank. See Barrlund (1998).

P12.3.5 Suppose $B \in \mathbb { R } ^ { n \times n }$ and $C \in \mathbb { R } ^ { m \times m }$ are unsymmetric and positive definite. Does it follow that $B \otimes C$ is positive definite?

P12.3.6 Show how to construct the normalized SVD of $B \otimes C$ from the normalized SVDs of B and C. Assume that $B \in \mathbb { R } ^ { m _ { B } \times n _ { B } }$ and $C \in \mathbb { R } ^ { m _ { C } \times n _ { C } }$ with m ${ \bf \beta } _ { B } \geq n _ { B }$ and $m _ { C } \geq n _ { C }$ .

P12.3.7 Show how to solve the linear system $( A \otimes B \otimes C ) x = d$ assuming that $A , B , C \in \mathbb { R } ^ { n \times n }$ are symmetric positive definite.

P12.3.8 (a) Given $A \in \mathbb { R } ^ { m n \times m n }$ and $B \in \mathbb { R } ^ { m \times m }$ , how would you compute $X \in \mathbb { R } ^ { n \times n }$ so that

$$
\phi_ {B} (X) = \left\| A - B \otimes X \right\| _ {F}
$$

is minimized? (b) Given $A \in \mathbb { R } ^ { m n \times m n }$ and $C \in \mathbb { R } ^ { n \times n } ,$ , how would you compute $\boldsymbol { X } \in \mathbb { R } ^ { m \times m }$ so that

$$
\phi_ {C} (X) = \left\| A - X \otimes C \right\| _ {F}
$$

is minimized?

P12.3.9 What is the nearest Kronecker product to the matrix $A = I _ { n } \otimes T _ { m } ^ { D D } + T _ { n } ^ { D D } \otimes I _ { n }$ where $\mathcal { T } _ { k } ^ { D D }$ is defined in (4.8.7).

P12.3.10 If A ∈ IRmn×mn is symmetric and tridiagonal, show how to minimize $\| A - B \otimes C \| _ { F }$ subject to the constraint that $B \in \mathbb { R } ^ { m \times m }$ and $C \in \mathbb { R } ^ { n \times n }$ are symmetric and tridiagonal.

P12.3.11 Prove Lemma 12.3.2. Hint: Show

$$
\parallel M - \alpha x x ^ {T} \parallel_ {F} ^ {2} = \parallel M \parallel_ {F} ^ {2} - 2 \alpha x ^ {T} T x + \alpha^ {2}
$$

where $T = ( M + M ^ { T } ) / 2$ .

P12.3.12 Prove Lemma 12.3.3. Hint: Show

$$
\parallel M - (x y ^ {T} - y x ^ {T}) \parallel_ {F} ^ {2} = \parallel M \parallel_ {F} ^ {2} + 2 \parallel x \parallel_ {2} ^ {2} \parallel y \parallel_ {2} ^ {2} - 2 (x ^ {T} y) ^ {2} - 4 x ^ {T} S y
$$

where $S = ( M - M ^ { T } ) / 2$ and use the real Schur form of S.

P12.3.13 For a symmetric matrix $S \in \mathbb { R } ^ { n \times n }$ , the symmetric vec operation is fully defined by

$$
S = \left[ \begin{array}{l l l} s _ {1 1} & s _ {1 2} & s _ {1 3} \\ s _ {2 1} & s _ {2 2} & s _ {2 3} \\ s _ {3 1} & s _ {3 2} & s _ {3 3} \end{array} \right] \Rightarrow \mathsf {s v e c} (S) = \left[ \begin{array}{l l l l l l l} s _ {1 1} & \sqrt {2}   s _ {2 1} & \sqrt {2}   s _ {3 1} & s _ {2 2} & \sqrt {2}   s _ {3 2} & s _ {3 3} \end{array} \right] ^ {T}.
$$

For symmetric $X \in \mathbb { R } ^ { n \times n }$ and arbitrary $B , C \in \mathbb { R } ^ { n \times n }$ , the symmetric Kronecker product is defined by

$$
(B \underset {\mathbf {S Y M}} {\otimes} C) \cdot \operatorname{svec} (X) = \operatorname{svec} \left(\frac {1}{2} \left(C X B ^ {T} + B X C ^ {T}\right)\right).
$$

For the case $n \ = \ 3 ,$ show that there is a matrix $P \in \mathbb { R } ^ { 9 \times 6 }$ with orthonormal columns so that $P ^ { T } ( B \otimes C ) P = B _ { \bf \Phi _ { S Y M } } C .$ SYM . See Vandenberge and Boyd (1996).

P12.3.14 The bi-alternate product is defined by

$$
B \underset {\mathbf {B I}} {\otimes} C = \frac {1}{2} (B \otimes C + C \otimes B).
$$

If $B = I , C = A$ , then solutions to $A X + X A ^ { T } \ = \ H$ where H is symmetric or skew-symmetric shed light on A’s eigenvalue placement. See Govaerts (2000). Given a matrix M, show how to compute the nearest bi-alternate product to M .

P12.3.15 Given $f \in \mathbb { R } ^ { q }$ and $g _ { i } \in \mathbb { R } ^ { \rho _ { i } }$ for $i = 1 { : } m$ , determine a permutation P so that

$$
P \left(f \otimes \left[ \begin{array}{c} g _ {1} \\ \vdots \\ g _ {m} \end{array} \right]\right) = \left[ \begin{array}{c} f \otimes g _ {1} \\ \vdots \\ f \otimes g _ {m} \end{array} \right].
$$

Hint: What does (12.3.1) say when B and C are vectors?

# Notes and References for 12.3

The history of the Kronecker product (including why it might better be called the “Zehfuss product”) is discussed in:

H.V. Henderson, F. Pukelsheim, and S.R. Searle (1983). “On the History of the Kronecker Product,” Lin. Mult. Alg. 14, 113–120.

For general background on the operation, see:

F. Stenger (1968), “Kronecker Product Extensions of Linear Operators,” SIAM J. Numer. Anal. 5, 422–435.   
J.W. Brewer (1978). “Kronecker Products and Matrix Calculus in System Theory,” IEEE Trans. Circuits Syst. 25, 772–781.   
A. Graham (1981). Kronecker Products and Matrix Calculus with Applications, Ellis Horwood, Chichester, England.   
M. Davio (1981), “Kronecker Products and Shuffle Algebra,” IEEE Trans. Comput. c-30, 116–125.   
H.V. Henderson and S.R. Searle (1981). “The Vec-Permutation Matrix, The Vec Operator and Kronecker Products: A Review,” Lin. Multilin. Alg. 9, 271–288.   
H.V. Henderson and S.R. Searle(1998). “Vec and Vech Operators for Matrices, with Some uses in Jacobians and Multivariate Statistics,” Canadian J. of Stat. 7, 65–81.   
C. Van Loan (2000). “The Ubiquitous Kronecker Product,” J. Comput. and Appl. Math. 123, 85–100.   
References concerned with various KP-like operations include:   
C.R. Rao and S.K. Mitra (1971). Generalized Inverse of Matrices and Applications, John Wiley and Sons, New York.   
D.S. Tracy and R.P. Singh (1972). “A New Matrix Product and Its Applications in Partitioned Matrices,” Statistica Neerlandica 26, 143–157.   
P.A. Regalia and S. Mitra (1989). “Kronecker Products, Unitary Matrices, and Signal Processing Applications,” SIAM Review 31, 586–613.   
J. Seberry and X-M Zhang (1993). “Some Orthogonal Matrices Constructed by Strong Kronecker Product Multiplication,” Austral. J. Combin. 7, 213–224.   
W. De Launey and J. Seberry (1994), “The Strong Kronecker Product,” J. Combin. Theory, Ser. A 66, 192–213.   
L. Vandenberghe and S. Boyd (1996). “Semidefinite Programming,” SIAM Review 38, 27–48.   
W. Govaerts (2000). Numerical Methods for Bifurcations of Dynamical Equilibria, SIAM Publications, Philadelphia, PA.   
A. Smilde, R. Bro, and P. Geladi (2004). Multiway Analysis, John Wiley, Chichester, England.   
For background on the KP connection to Sylvester-type equations, see:   
P. Lancaster (1970). “Explicit Solution of Linear Matrix Equations,” SIAM Review 12, 544–566.   
W.J. Vetter (1975). “Vector Structures and Solutions of Linear Matrix Equations,” Lin. Alg. Applic. 10, 181–188.   
Issues associated with the efficient implementation of KP operations are discussed in:   
H.C. Andrews and J. Kane (1970). “Kronecker Matrices, Computer Implementation, and Generalized Spectra,” J. Assoc. Comput. Mach. 17, 260–268.   
V. Pereyra and G. Scherer (1973). “Efficient Computer Manipulation of Tensor Products with Applications to Multidimensional Approximation,” Math. Comput. 27, 595–604.   
C. de Boor (1979). “Efficient Computer Manipulation of Tensor Products,” ACM Trans. Math. Softw. 5, 173–182.   
P.E. Buis and W.R. Dyksen (1996). “Efficient Vector and Parallel Manipulation of Tensor Products,” ACM Trans. Math. Softw. 22, 18–23.   
P.E. Buis and W.R. Dyksen (1996). “Algorithm 753: TENPACK: An LAPACK-based Library for the Computer Manipulation of Tensor Products,” ACM Trans. Math. Softw. 22, 24–29.   
W-H. Steeb (1997). Matrix Calculus and Kronecker Product with Applications and C++ Programs, World Scientific Publishing, Singapore.   
M. Huhtanen (2006). “Real Linear Kronecker Product Operations,” Lin. Alg. Applic. 417, 347–361.   
The KP is associated with the vast majority fast linear transforms. See Van Loan (FFT) as well as:   
C-H Huang, J.R. Johnson, and R.W. Johnson (1991). “Multilinear Algebra and Parallel Programming,” J. Supercomput. 5, 189–217.   
J. Granata, M. Conner, and R. Tolimieri (1992). “‘Recursive Fast Algorithms and the Role of the Tensor Product,” IEEE Trans. Signal Process. 40, 2921–2930.   
J. Granata, M. Conner, and R. Tolimieri (1992). “The Tensor Product: A Mathematical Programming Language for FFTs and Other Fast DSP Operations,” IEEE SP Magazine, January, 40–48.   
For a discussion of the role of KP approximation in a variety of situations, see:

C. Van Loan and N.P Pitsianis (1992). “Approximation with Kronecker Products”, in Linear Algebra for Large Scale and Real Time Applications, M.S. Moonen and G.H. Golub (eds.), Kluwer Publications, Dordrecht, 293–314,   
T.F. Andre, R.D. Nowak, and B.D. Van Veen (1997). “Low Rank Estimation of Higher Order Statistics,” IEEE Trans. Signal Process. 45, 673–685.   
R.D. Nowak and B. Van Veen (1996). “Tensor Product Basis Approximations for Volterra Filters,” IEEE Trans. Signal Process. 44, 36–50.   
J. Kamm and J.G. Nagy (1998). “Kronecker Product and SVD Approximations in Image Restoration,” Lin. Alg. Applic. 284, 177–192.   
J.G. Nagy and D.P. O’Leary (1998). “Restoring Images Degraded by Spatially Variant Blur,” SIAM J. Sci. Comput. 19, 1063–1082.   
J. Kamm and J.G. Nagy (2000). “Optimal Kronecker Product Approximation of Block Toeplitz Matrices,” SIAM J. Matrix Anal. Applic. 22, 155–172.   
J.G. Nagy, M.K. Ng, and L. Perrone (2003). “Kronecker Product Approximations for Image Restoration with Reflexive Boundary Conditions,” SIAM J. Matrix Anal. Applic. 25, 829–841.   
A.N. Langville and W.J. Stewart (2004). “A Kronecker Product Approximate Preconditioner for SANs,” Num. Lin. Alg. 11, 723–752.   
E. Tyrtyshnikov (2004). “Kronecker-Product Approximations for Some Function-Related Matrices,” Lin. Alg. Applic. 379, 423–437.   
L. Perrone (2005). “Kronecker Product Approximations for Image Restoration with Anti-Reflective Boundary Conditions,” Num. Lin. Alg. 13, 1–22.   
W. Hackbusch, B.N. Khoromskij, and E.E. Tyrtyshnikov (2005). “Hierarchical Kronecker Tensor-Product Approximations,” J. Numer. Math. 13, 119–156.   
V. Olshevsky, I. Oseledets, and E. Tyrtyshnikov (2006). “Tensor Properties of Multilevel Toeplitz and Related matrices,” Lin. Alg. Applic. 412, 1–21.   
J. Leskovec and C. Faloutsos (2007). “Scalable Modeling of Real Graphs Using Kronecker Multiplication,” in Proc. of the 24th International Conference on Machine Learning, Corvallis, OR.   
J. Leskovic (2011). “Kronecker Graphs,” in Graph Algorithms in the Language of Linear Algebra, J. Kepner and J. Gilbert (eds), SIAM Publications, Philadelphia, PA, 137–204.   
For a snapshot of KP algorithms for linear systems and least squares problems, see:   
H. Sunwoo (1996). “Simple Algorithms about Kronecker Products in the Linear Model,” Lin. Alg. Applic. 237–8, 351–358.   
D.W. Fausett, C.T. Fulton, and H. Hashish (1997). “Improved Parallel QR Method for Large Least Squares Problems Involving Kronecker Products,” J. Comput. Appl. Math. 78, 63–78.   
A. Barrlund (1998). “Efficient Solution of Constrained Least Squares Problems with Kronecker Product Structure,” SIAM J. Matrix Anal. Applic. 19, 154–160.   
P. Buchholz and T.R. Dayar (2004). “Block SOR for Kronecker Structured Representations,” Lin. Alg. Applic. 386, 83–109.   
A.W. Bojanczyk and A. Lutoborski (2003). “The Procrustes Problem for Orthogonal Kronecker Products,” SIAM J. Sci. Comput. 25, 148–163.   
C.D.M. Martin and C.F. Van Loan (2006). “Shifted Kronecker Product Systems,” SIAM J. Matrix Anal. Applic. 29, 184–198.
