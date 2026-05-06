# 1.3 Block Matrices and Algorithms

A block matrix is a matrix whose entries are themselves matrices. It is a point of view. For example, an 8-by-15 matrix of scalars can be regarded as a 2-by-3 block matrix with 4-by-5 entries. Algorithms that manipulate matrices at the block level are often more efficient because they are richer in level-3 operations. The derivation of many important algorithms is often simplified by using block matrix notation.

# 1.3.1 Block Matrix Terminology

Column and row partitionings (§1.1.7) are special cases of matrix blocking. In general, we can partition both the rows and columns of an m-by-n matrix A to obtain

$$
A = \left[ \begin{array}{c c c} A _ {1 1} & \ldots & A _ {1 r} \\ \vdots & & \vdots \\ A _ {q 1} & \dots & A _ {q r} \\ n _ {1} & & n _ {r} \end{array} \right] \begin{array}{c} m _ {1} \\ m _ {q} \end{array}
$$

where $m _ { 1 } + \cdots + m _ { q } = m , n _ { 1 } + \cdots + n _ { r } = n$ , and $A _ { \alpha \beta }$ designates the $( \alpha , \beta )$ block (submatrix). With this notation, block $A _ { \alpha \beta }$ has dimension $m _ { \alpha } – \mathrm { b y } – n _ { \beta }$ and we say that $A = \left( A _ { \alpha \beta } \right)$ is a q-by-r block matrix.

Terms that we use to describe well-known band structures for matrices with scalar entries have natural block analogs. Thus,

$$
\operatorname{diag} (A _ {1 1}, A _ {2 2}, A _ {3 3}) = \left[ \begin{array}{c c c} A _ {1 1} & 0 & 0 \\ 0 & A _ {2 2} & 0 \\ 0 & 0 & A _ {3 3} \end{array} \right]
$$

is block diagonal while the matrices

$$
L = \left[ \begin{array}{c c c} L _ {1 1} & 0 & 0 \\ L _ {2 1} & L _ {2 2} & 0 \\ L _ {3 1} & L _ {3 2} & L _ {3 3} \end{array} \right], \quad U = \left[ \begin{array}{c c c} U _ {1 1} & U _ {1 2} & U _ {1 3} \\ 0 & U _ {2 2} & U _ {2 3} \\ 0 & 0 & U _ {3 3} \end{array} \right], \quad T = \left[ \begin{array}{c c c} T _ {1 1} & T _ {1 2} & 0 \\ T _ {2 1} & T _ {2 2} & T _ {2 3} \\ 0 & T _ {3 2} & T _ {3 3} \end{array} \right],
$$

are, respectively, block lower triangular, block upper triangular, and block tridiagonal. The blocks do not have to be square in order to use this block sparse terminology.

# 1.3.2 Block Matrix Operations

Block matrices can be scaled and transposed:

$$
\begin{array}{l} \mu \left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \\ A _ {3 1} & A _ {3 2} \end{array} \right] = \left[ \begin{array}{c c} \mu A _ {1 1} & \mu A _ {1 2} \\ \mu A _ {2 1} & \mu A _ {2 2} \\ \mu A _ {3 1} & \mu A _ {3 2} \end{array} \right], \\ \left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \\ A _ {3 1} & A _ {3 2} \end{array} \right] ^ {T} = \left[ \begin{array}{c c c} A _ {1 1} ^ {T} & A _ {2 1} ^ {T} & A _ {3 1} ^ {T} \\ A _ {1 2} ^ {T} & A _ {2 2} ^ {T} & A _ {3 2} ^ {T} \end{array} \right]. \\ \end{array}
$$

Note that the transpose of the original $( i , j )$ block becomes the $( j , i )$ block of the result. Identically blocked matrices can be added by summing the corresponding blocks:

$$
\left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \\ A _ {3 1} & A _ {3 2} \end{array} \right] + \left[ \begin{array}{c c} B _ {1 1} & B _ {1 2} \\ B _ {2 1} & B _ {2 2} \\ B _ {3 1} & B _ {3 2} \end{array} \right] = \left[ \begin{array}{c c} A _ {1 1} + B _ {1 1} & A _ {1 2} + B _ {1 2} \\ A _ {2 1} + B _ {2 1} & A _ {2 2} + B _ {2 2} \\ A _ {3 1} + B _ {3 1} & A _ {3 2} + B _ {3 2} \end{array} \right].
$$

Block matrix multiplication requires more stipulations about dimension. For example, if

$$
\left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \\ A _ {3 1} & A _ {3 2} \end{array} \right] \left[ \begin{array}{c c} B _ {1 1} & B _ {1 2} \\ B _ {2 1} & B _ {2 2} \end{array} \right] = \left[ \begin{array}{c c} A _ {1 1} B _ {1 1} + A _ {1 2} B _ {2 1} & A _ {1 1} B _ {1 2} + A _ {1 2} B _ {2 2} \\ A _ {2 1} B _ {1 1} + A _ {2 2} B _ {2 1} & A _ {2 1} B _ {1 2} + A _ {2 2} B _ {2 2} \\ A _ {3 1} B _ {1 1} + A _ {3 2} B _ {2 1} & A _ {3 1} B _ {1 2} + A _ {3 2} B _ {2 2} \end{array} \right]
$$

is to make sense, then the column dimensions of $A _ { 1 1 } , A _ { 2 1 }$ , and $A _ { 3 1 }$ must each be equal to the row dimension of both $B _ { 1 1 }$ and $B _ { 1 2 }$ . Likewise, the column dimensions of $A _ { 1 2 }$ , $A _ { 2 2 }$ , and $A _ { 3 2 }$ must each be equal to the row dimensions of both $B _ { 2 1 }$ and $B _ { 2 2 }$ .

Whenever a block matrix addition or multiplication is indicated, it is assumed that the row and column dimensions of the blocks satisfy all the necessary constraints. In that case we say that the operands are partitioned conformably as in the following theorem.

Theorem 1.3.1. $I f$

$$
A = \left[ \begin{array}{c c c} A _ {1 1} & \ldots & A _ {1 s} \\ \vdots & & \vdots \\ A _ {q 1} & \dots & A _ {q s} \end{array} \right] _ {p _ {1}} ^ {m _ {1}}, \qquad B = \left[ \begin{array}{c c c} B _ {1 1} & \ldots & B _ {1 r} \\ \vdots & & \vdots \\ B _ {s 1} & \dots & B _ {s r} \end{array} \right] _ {n _ {1}} ^ {p _ {1}},
$$

and we partition the product $C = A B$ as follows,

$$
C = \left[ \begin{array}{c c c} C _ {1 1} & \ldots & C _ {1 r} \\ \vdots & & \vdots \\ C _ {q 1} & \dots & C _ {q r} \\ n _ {1} & & n _ {r} \end{array} \right] _ {m _ {q}} ^ {m _ {1}},
$$

then for $\alpha = 1 { : } q$ and $\beta = 1 { : } r$ we have $C _ { \alpha \beta } = \sum _ { \gamma = 1 } ^ { s } A _ { \alpha \gamma } B _ { \gamma \beta }$ .

Proof. The proof is a tedious exercise in subscripting. Suppose $1 \leq \alpha \leq q$ and $1 \leq \beta \leq r$ . Set $M = m _ { 1 } + \cdot \cdot \cdot + m _ { \alpha - 1 }$ and $N = n _ { 1 } + \cdot \cdot \cdot n _ { \beta - 1 }$ . It follows that if $1 \leq i \leq m _ { \alpha }$ and $1 \leq j \leq n _ { \beta }$ then

$$
\begin{array}{l} [ C _ {\alpha \beta} ] _ {i j} = \sum_ {k = 1} ^ {p _ {1} + \dots p _ {s}} a _ {M + i, k} b _ {k, N + j} = \sum_ {\gamma = 1} ^ {s} \sum_ {k = p _ {1} + \dots + p _ {\gamma - 1} + 1} ^ {p _ {1} + \dots + p _ {\gamma}} a _ {M + i, k} b _ {k, N + j} \\ = \sum_ {\gamma = 1} ^ {s} \sum_ {k = 1} ^ {p _ {\gamma}} [ A _ {\alpha \gamma} ] _ {i k} [ B _ {\gamma \beta} ] _ {k j} = \sum_ {\gamma = 1} ^ {s} [ A _ {\alpha \gamma} B _ {\gamma \beta} ] _ {i j} = \left[ \sum_ {\gamma = 1} ^ {s} A _ {\alpha \gamma} B _ {\gamma \beta} \right] _ {i j}. \\ \end{array}
$$

Thus, $C _ { \alpha \beta } = A _ { \alpha , 1 } B _ { 1 , \beta } + \cdot \cdot \cdot + A _ { \alpha , s } B _ { s , \beta } .$

If you pay attention to dimension and remember that matrices do not commute, i.e., $A _ { 1 1 } B _ { 1 1 } + A _ { 1 2 } B _ { 2 1 } \ne B _ { 1 1 } A _ { 1 1 } + B _ { 2 1 } A _ { 1 2 }$ , then block matrix manipulation is just ordinary matrix manipulation with the $\arcsin { \mathrm { \Omega } } a _ { i j } { \mathrm { \Omega } } ^ { \prime } \mathrm { { s } }$ and $b _ { i j } \mathrm { ^ { \prime } s }$ written as $A _ { i j } \mathrm { ^ { , } s }$ and $B _ { i j } \mathrm { ^ { , } s ! }$

# 1.3.3 Submatrices

Suppose $A \in \mathbb { R } ^ { m \times n }$ . If $\alpha = [ \alpha _ { 1 } , \dots , \alpha _ { s } ]$ and $\beta = [ \beta _ { 1 } , \dots , \beta _ { t } ]$ are integer vectors with distinct components that satisfy $1 \leq \alpha _ { i } \leq m$ , and $1 \leq \beta _ { i } \leq n$ , then

$$
A (\alpha , \beta) = \left[ \begin{array}{c c c} a _ {\alpha_ {1}, \beta_ {1}} & \dots & a _ {\alpha_ {1}, \beta_ {t}} \\ \vdots & \ddots & \vdots \\ a _ {\alpha_ {s}, \beta_ {1}} & \dots & a _ {\alpha_ {s}, \beta_ {t}} \end{array} \right]
$$

is an s-by-t submatrix of A. For example, if $A \in \mathbb { R } ^ { 8 \times 6 } , \alpha = [ 2 4 6 8 ]$ , and $\beta = [ 4 5 6 ]$ then

$$
A (\alpha , \beta) = \left[ \begin{array}{c c c} a _ {2 4} & a _ {2 5} & a _ {2 6} \\ a _ {4 4} & a _ {4 5} & a _ {4 6} \\ a _ {6 4} & a _ {6 5} & a _ {6 6} \\ a _ {8 4} & a _ {8 5} & a _ {8 6} \end{array} \right].
$$

If $\alpha = \beta$ , then $A ( \alpha , \beta )$ is a principal submatrix. If $\alpha = \beta = 1 { : } k$ and $1 \leq k \leq \operatorname* { m i n } \{ m , n \}$ , then $A ( \alpha , \beta )$ is a leading principal submatrix.

If $A \in \mathbb { R } ^ { m \times n }$ and

$$
A = \left[ \begin{array}{c c c} A _ {1 1} & \ldots & A _ {1 s} \\ \vdots & & \vdots \\ A _ {q 1} & \dots & A _ {q s} \\ n _ {1} & & n _ {r} \end{array} \right] \begin{array}{c} m _ {1} \\ m _ {q} \end{array} ,
$$

then the colon notation can be used to specify the individual blocks. In particular,

$$
A _ {i j} = A (\tau + 1: \tau + m _ {i}, \mu + 1: \mu + n _ {j})
$$

where $\tau = m _ { 1 } + \cdot \cdot \cdot + m _ { i - 1 }$ and $\mu = n _ { 1 } + \cdot \cdot \cdot + n _ { j - 1 }$ . Block matrix notation is valuable for the way in which it hides subscript range expressions.

# 1.3.4 The Blocked Gaxpy

As an exercise in block matrix manipulation and submatrix designation, we consider two block versions of the gaxpy operation $y ~ = ~ y + A x$ where $A \in \mathbb { R } ^ { m \times n } , x \in \mathbb { R } ^ { n }$ , and $\boldsymbol { y } \in \mathbb { R } ^ { m }$ . If

$$
A = \left[ \begin{array}{c} A _ {1} \\ \vdots \\ A _ {q} \end{array} \right] _ {m _ {q}} ^ {m _ {1}} \qquad \text { and } \qquad y = \left[ \begin{array}{c} y _ {1} \\ \vdots \\ y _ {q} \end{array} \right] _ {m _ {q}} ^ {m _ {1}},
$$

then

$$
\left[ \begin{array}{c} y _ {1} \\ \vdots \\ y _ {q} \end{array} \right] = \left[ \begin{array}{c} y _ {1} \\ \vdots \\ y _ {q} \end{array} \right] + \left[ \begin{array}{c} A _ {1} \\ \vdots \\ A _ {q} \end{array} \right] x,
$$

and we obtain

$$
\alpha = 0
$$

for i = 1:q

$$
i d x = \alpha + 1: \alpha + m _ {i}
$$

$$
y (i d x) = y (i d x) + A (i d x,:) \cdot x
$$

$$
\alpha = \alpha + m _ {i}
$$

end

The assignment to $y ( i d x )$ corresponds to $y _ { i } = y _ { i } + A _ { i } x$ . This row-blocked version of the gaxpy computation breaks the given gaxpy into q “shorter” gaxpys. We refer to $A _ { i }$ as the ith block row of A.

Likewise, with the partitionings

$$
A = \left[ \begin{array}{c c} A _ {1} & \dots & A _ {r} \\ n _ {1} & n _ {r} \end{array} \right] \qquad \text {and} \qquad x = \left[ \begin{array}{c} x _ {1} \\ \vdots \\ x _ {r} \end{array} \right] \begin{array}{c} n _ {1} \\ n _ {r} \end{array} ,
$$

we see that

$$
y = y + \left[ A _ {1} \mid \dots \mid A _ {r} \right] \left[ \begin{array}{c} x _ {1} \\ \vdots \\ x _ {r} \end{array} \right] = y + \sum_ {j = 1} ^ {r} A _ {j} x _ {j}
$$

and we obtain

$$
\beta = 0
$$

for j = 1:r

$$
j d x = \beta + 1: \beta + n _ {j}
$$

$$
y = y + A (:, j d x) \cdot x (j d x)
$$

$$
\beta = \beta + n _ {j}
$$

end

The assignment to y corresponds to $y = y + A _ { j } x _ { j }$ . This column-blocked version of the gaxpy computation breaks the given gaxpy into r “thinner” gaxpys. We refer to $A _ { j }$ as the jth block column of A.


---

<!-- golub_050_099 -->

Just as ordinary, scalar-level matrix multiplication can be arranged in several possible ways, so can the multiplication of block matrices. To illustrate this with a minimum of subscript clutter, we consider the update

$$
C = C + A B
$$

where we regard $A = ( A _ { \alpha \beta } ) , B = ( B _ { \alpha \beta } )$ , and $C = ( C _ { \alpha \beta } )$ as N -by-N block matrices with -by- blocks. From Theorem 1.3.1 we have

$$
C _ {\alpha \beta} = C _ {\alpha \beta} + \sum_ {\gamma = 1} ^ {N} A _ {\alpha \gamma} B _ {\gamma \beta}, \quad \alpha = 1: N, \quad \beta = 1: N.
$$

If we organize a matrix multiplication procedure around this summation, then we obtain a block analog of Algorithm 1.1.5:

for $\alpha = 1:N$ $i = (\alpha - 1)\ell + 1:\alpha\ell$ for $\beta = 1:N$ $j = (\beta - 1)\ell + 1:\beta\ell$ for $\gamma = 1:N$ $k = (\gamma - 1)\ell + 1:\gamma\ell$ $C(i, j) = C(i, j) + A(i, k) \cdot B(k, j)$ end

end

end

Note that, if  = 1, then $\alpha \equiv i , \beta \equiv j$ , and $\gamma \equiv k$ and we revert to Algorithm 1.1.5.

Analogously to what we did in §1.1, we can obtain different variants of this procedure by playing with loop orders and blocking strategies. For example, corresponding to

$$
\left[ \begin{array}{c c c} C _ {1 1} & \dots & C _ {1 N} \\ \vdots & \ddots & \vdots \\ C _ {N 1} & \dots & C _ {N N} \end{array} \right] + \left[ \begin{array}{c} A _ {1} \\ \vdots \\ A _ {N} \end{array} \right] \left[ \begin{array}{c c c} B _ {1} & \dots & B _ {N} \end{array} \right]
$$

where $A _ { i } \in \mathbb { R } ^ { \ell \times n }$ and $B _ { j } \in \mathbb { R } ^ { n \times \ell }$ , we obtain the following block outer product computation:

for $i = 1:N$ for $j = 1:N$ $C_{ij} = C_{ij} + A_iB_j$ end end

# 1.3.6 The Kronecker Product

It is sometimes the case that the entries in a block matrix A are all scalar multiples of the same matrix. This means that A is a Kronecker product. Formally, if $B \in \mathbb { R } ^ { m _ { 1 } }$ ×n1 and $C \in \mathbb { R } ^ { m _ { 2 } \times n _ { 2 } }$ , then their Kronecker product $B \otimes C$ is an $m _ { 1 } { \mathrm { - b y } } { \mathrm { - } } n _ { 1 }$ block matrix whose (i, j) block is the $m _ { 2 } { \mathrm { - b y } } { \mathrm { - } } n _ { 2 }$ matrix $b _ { i j } C$ . Thus, if

$$
A = \left[ \begin{array}{l l} b _ {1 1} & b _ {1 2} \\ b _ {2 1} & b _ {2 2} \\ b _ {3 1} & b _ {3 2} \end{array} \right] \otimes \left[ \begin{array}{l l l} c _ {1 1} & c _ {1 2} & c _ {1 3} \\ c _ {2 1} & c _ {2 2} & c _ {2 3} \\ c _ {3 1} & c _ {3 2} & c _ {3 3} \end{array} \right]
$$

then

$$
A   =   \left[ \begin{array}{c c c c c c} b _ {1 1} c _ {1 1} & b _ {1 1} c _ {1 2} & b _ {1 1} c _ {1 3} & b _ {1 2} c _ {1 1} & b _ {1 2} c _ {1 2} & b _ {1 2} c _ {1 3} \\ b _ {1 1} c _ {2 1} & b _ {1 1} c _ {2 2} & b _ {1 1} c _ {2 3} & b _ {1 2} c _ {2 1} & b _ {1 2} c _ {2 2} & b _ {1 2} c _ {2 3} \\ b _ {1 1} c _ {3 1} & b _ {1 1} c _ {3 2} & b _ {1 1} c _ {3 3} & b _ {1 2} c _ {3 1} & b _ {1 2} c _ {3 2} & b _ {1 2} c _ {3 3} \\ \hline b _ {2 1} c _ {1 1} & b _ {2 1} c _ {1 2} & b _ {2 1} c _ {1 3} & b _ {2 2} c _ {1 1} & b _ {2 2} c _ {1 2} & b _ {2 2} c _ {1 3} \\ b _ {2 1} c _ {2 1} & b _ {2 1} c _ {2 2} & b _ {2 1} c _ {2 3} & b _ {2 2} c _ {2 1} & b _ {2 2} c _ {2 2} & b _ {2 2} c _ {2 3} \\ b _ {2 1} c _ {3 1} & b _ {2 1} c _ {3 2} & b _ {2 1} c _ {3 3} & b _ {2 2} c _ {3 1} & b _ {2 2} c _ {3 2} & b _ {2 2} c _ {3 3} \\ \hline b _ {3 1} c _ {1 1} & b _ {3 1} c _ {1 2} & b _ {3 1} c _ {1 3} & b _ {3 2} c _ {1 1} & b _ {3 2} c _ {1 2} & b _ {3 2} c _ {1 3} \\ b _ {3 1} c _ {2 1} & b _ {3 1} c _ {2 2} & b _ {3 1} c _ {2 3} & b _ {3 2} c _ {2 1} & b _ {3 2} c _ {2 2} & b _ {3 2} c _ {2 3} \\ b _ {3 1} c _ {3 1} & b _ {3 1} c _ {3 2} & b _ {3 1} c _ {3 3} & b _ {3 2} c _ {3 1} & b _ {3 2} c _ {3 2} & b _ {3 2} c _ {3 3} \end{array} \right].
$$

This type of highly structured blocking occurs in many applications and results in dramatic economies when fully exploited.

Note that if B has a band structure, then $B \otimes C$ “inherits” that structure at the block level. For example, if

$$
B \text {is} \left\{ \begin{array}{l} \text {diagonal} \\ \text {tridiagonal} \\ \text {lower triangular} \\ \text {upper triangular} \end{array} \right\} \text {then} B \otimes C \text {is} \left\{ \begin{array}{l} \text {block diagonal} \\ \text {block tridiagonal} \\ \text {block lower triangular} \\ \text {block upper triangular} \end{array} \right\}.
$$

Important Kronecker product properties include:

$$
(B \otimes C) ^ {T} = B ^ {T} \otimes C ^ {T}, \tag {1.3.1}
$$

$$
(B \otimes C) (D \otimes F) = B D \otimes C F, \tag {1.3.2}
$$

$$
(B \otimes C) ^ {- 1} = B ^ {- 1} \otimes C ^ {- 1}, \tag {1.3.3}
$$

$$
B \otimes (C \otimes D) = (B \otimes C) \otimes D. \tag {1.3.4}
$$

Of course, the products BD and CF must be defined for (1.3.2) to make sense. Likewise, the matrices B and C must be nonsingular in (1.3.3).

In general, $B \otimes C \neq C \otimes B$ . However, there is a connection between these two matrices via the perfect shuffle permutation that is defined in §1.2.11. If $B \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } }$ and $C \in \mathbb { R } ^ { m _ { 2 } \times n _ { 2 } }$ , then

$$
P (B \otimes C) Q ^ {T} = C \otimes B \tag {1.3.5}
$$

where $P = \mathcal { P } _ { m _ { 1 } , m _ { 2 } }$ and $Q \ = \ \mathcal P _ { n _ { 1 } , n _ { 2 } }$ .

# 1.3.7 Reshaping Kronecker Product Expressions

A matrix-vector product in which the matrix is a Kronecker product is “secretly” a matrix-matrix-matrix product. For example, if $B \in \mathbb { R } ^ { 3 \times 2 } , C \in \bar { \mathbb { R } } ^ { m \times n }$ , and $x _ { 1 } , x _ { 2 } \in \mathbb { R } ^ { n }$ , then

$$
\begin{array}{l} \left[ \begin{array}{c} y _ {1} \\ y _ {2} \\ y _ {3} \end{array} \right] = (B \otimes C) \left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right] = \left[ \begin{array}{c c} b _ {1 1} C & b _ {1 2} C \\ b _ {2 1} C & b _ {2 2} C \\ b _ {3 1} C & b _ {3 2} C \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right] \\ = \left[ \begin{array}{c} b _ {1 1} C x _ {1} + b _ {1 2} C x _ {2} \\ b _ {2 1} C x _ {1} + b _ {2 2} C x _ {2} \\ b _ {3 1} C x _ {1} + b _ {3 2} C x _ {2} \end{array} \right] \\ \end{array}
$$

where $y _ { 1 } , y _ { 2 } , y _ { 3 } \in \mathbb { R } ^ { m }$ . On the other hand, if we define the matrices

$$
X = \left[ \begin{array}{c c} x _ {1} & x _ {2} \end{array} \right] \quad \text { and } \quad Y = \left[ \begin{array}{c c c} y _ {1} & y _ {2} & y _ {3} \end{array} \right],
$$

then $Y = C X B ^ { T }$ .

To be precise about this reshaping, we introduce the vec operation. If $\boldsymbol { X } \in \mathbb { R } ^ { m \times n }$ , then $\mathsf { v e c } ( X )$ is an nm-by-1 vector obtained by “stacking” X’s columns:

$$
\operatorname{vec} (X) = \left[ \begin{array}{c} X (:, 1) \\ \vdots \\ X (:, n) \end{array} \right].
$$

If B ∈ IRm1×n1, C ∈ IRm2×n2, and $X \in \mathbb { R } ^ { n _ { 1 } \times m _ { 2 } }$ , then

$$
Y = C X B ^ {T} \Leftrightarrow \operatorname{vec} (Y) = (B \otimes C) \operatorname{vec} (X). \tag {1.3.6}
$$

Note that if $B , C , X \in \mathbb { R } ^ { n \times n }$ , then $Y = C X B ^ { T }$ costs $O ( n ^ { 3 } )$ to evaluate while the disregard of Kronecker structure in $y = ( B \otimes C ) x$ leads to an $O ( n ^ { 4 } )$ calculation. This is why reshaping is central for effective Kronecker product computation. The reshape operator is handy in this regard. If $A \in \mathbb { R } ^ { m \times n }$ and $m _ { 1 } n _ { 1 } = m n$ , then

$$
B = \operatorname{reshape} (A, m _ {1}, n _ {1})
$$

is the $m _ { 1 } { \mathrm { - b y } } { \mathrm { - } } n _ { 1 }$ matrix defined by $\mathsf { v e c } ( B ) = \mathsf { v e c } ( A )$ . Thus, if $\ b { A } \in \mathbb { R } ^ { 3 \times 4 }$ , then

$$
\operatorname{reshape} (A, 2, 6) = \left[ \begin{array}{c c c c c c} a _ {1 1} & a _ {3 1} & a _ {2 2} & a _ {1 3} & a _ {3 3} & a _ {2 4} \\ a _ {2 1} & a _ {1 2} & a _ {3 2} & a _ {2 3} & a _ {1 4} & a _ {3 4} \end{array} \right].
$$

# 1.3.8 Multiple Kronecker Products

Note that $A = B \otimes C \otimes D$ can be regarded as a block matrix whose entries are block matrices. In particular, $b _ { i j } c _ { k \ell } D$ is the $( k , \ell )$ block of $A ^ { \prime } \mathrm { s } \ ( i , j )$ block.

As an example of a multiple Kronecker product computation, let us consider the calculation of $y = ( B \otimes C \otimes D ) x$ where B, $C , D \in \mathbb { R } ^ { n \times n }$ and $\boldsymbol { x } \in \mathbb { R } ^ { N }$ with $N = n ^ { 3 }$ . Using (1.3.6) it follows that

$$
\operatorname{reshape} (y, n ^ {2}, n) = (C \otimes D) \cdot \operatorname{reshape} (x, n ^ {2}, n) \cdot B ^ {T}.
$$

Thus, if

$$
F = \operatorname{reshape} (x, n ^ {2}, n) \cdot B ^ {T},
$$

then $G = ( C \otimes D ) F \in \mathbb { R } ^ { n ^ { 2 } \times n }$ can computed column-by-column using (1.3.6):

$$
G (:, k) = \operatorname{reshape} (D \cdot \operatorname{reshape} (F (:, k), n, n) \cdot C ^ {T}, n ^ {2}, 1) \quad k = 1: n.
$$

It follows that $y = { \mathsf { r e s h a p e } } ( G , N , 1 )$ . A careful accounting reveals that $6 n ^ { 4 }$ flops are required. Ordinarily, a matrix-vector product of this dimension would require $2 n ^ { 6 }$ flops.

The Kronecker product has a prominent role to play in tensor computations and in §13.1 we detail more of its properties.

# 1.3.9 A Note on Complex Matrix Multiplication

Consider the complex matrix multiplication update

$$
C _ {1} + i C _ {2} = (C _ {1} + i C _ {2}) + (A _ {1} + i A _ {2}) (B _ {1} + i B _ {2})
$$

where all the matrices are real and $i ^ { 2 } = - 1$ . Comparing the real and imaginary parts we conclude that

$$
\left[ \begin{array}{l} C _ {1} \\ C _ {2} \end{array} \right] = \left[ \begin{array}{l} C _ {1} \\ C _ {2} \end{array} \right] + \left[ \begin{array}{c c} A _ {1} & - A _ {2} \\ A _ {2} & A _ {1} \end{array} \right] \left[ \begin{array}{l} B _ {1} \\ B _ {2} \end{array} \right].
$$

Thus, complex matrix multiplication corresponds to a structured real matrix multiplication that has expanded dimension.

# 1.3.10 Hamiltonian and Symplectic Matrices

While on the topic of 2-by-2 block matrices, we identify two classes of structured matrices that arise at various points later on in the text. A matrix M ∈ IR2n×2n $\boldsymbol { M } \in \mathbb { R } ^ { 2 n \times 2 n }$ is a Hamiltonian matrix if it has the form

$$
M = \left[ \begin{array}{c c} A & G \\ F & - A ^ {T} \end{array} \right]
$$

where A, $F , G \in \mathbb { R } ^ { n \times n }$ and F and G are symmetric. Hamiltonian matrices arise in optimal control and other application areas. An equivalent definition can be given in terms of the permutation matrix

$$
J = \left[ \begin{array}{c c} 0 & I _ {n} \\ - I _ {n} & 0 \end{array} \right].
$$

In particular, if

$$
J M J ^ {T} = - M ^ {T},
$$

then M is Hamiltonian. A related class of matrices are the symplectic matrices. A matrix $S \in \mathbb { R } ^ { 2 n \times 2 n }$ is symplectic if

$$
S ^ {T} J S = J.
$$

If

$$
S = \left[ \begin{array}{l l} S _ {1 1} & S _ {1 2} \\ S _ {2 1} & S _ {2 2} \end{array} \right]
$$

where the blocks are $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ , then it follows that both $S _ { 1 1 } ^ { T } S _ { 2 1 }$ and $S _ { 2 2 } ^ { T } S _ { 1 2 }$ are symmetric and $S _ { 1 1 } ^ { T } S _ { 2 2 } = I _ { n } + S _ { 2 1 } ^ { T } S _ { 1 2 }$ .

# 1.3.11 Strassen Matrix Multiplication

We conclude this section with a completely different approach to the matrix-matrix multiplication problem. The starting point in the discussion is the 2-by-2 block matrix product

$$
{\left[ \begin{array}{l l} C _ {1 1} & C _ {1 2} \\ C _ {2 1} & C _ {2 2} \end{array} \right]} = {\left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right]} {\left[ \begin{array}{l l} B _ {1 1} & B _ {1 2} \\ B _ {2 1} & B _ {2 2} \end{array} \right]}
$$

where each block is square. In the ordinary algorithm, $C _ { i j } = A _ { i 1 } B _ { 1 j } + A _ { i 2 } B _ { 2 j }$ . There are 8 multiplies and 4 adds. Strassen (1969) has shown how to compute $C$ with just 7 multiplies and 18 adds:

$$
\begin{array}{l} P _ {1} = (A _ {1 1} + A _ {2 2}) (B _ {1 1} + B _ {2 2}), \\ P _ {2} = (A _ {2 1} + A _ {2 2}) B _ {1 1}, \\ P _ {3} = A _ {1 1} (B _ {1 2} - B _ {2 2}), \\ P _ {4} = A _ {2 2} (B _ {2 1} - B _ {1 1}), \\ P _ {5} = (A _ {1 1} + A _ {1 2}) B _ {2 2}, \\ P _ {6} = (A _ {2 1} - A _ {1 1}) (B _ {1 1} + B _ {1 2}), \\ P _ {7} = \left(A _ {1 2} - A _ {2 2}\right) \left(B _ {2 1} + B _ {2 2}\right), \\ C _ {1 1} = P _ {1} + P _ {4} - P _ {5} + P _ {7}, \\ C _ {1 2} = P _ {3} + P _ {5}, \\ C _ {2 1} = P _ {2} + P _ {4}, \\ C _ {2 2} = P _ {1} + P _ {3} - P _ {2} + P _ {6}. \\ \end{array}
$$

These equations are easily confirmed by substitution. Suppose $n = 2 m$ so that the blocks are m-by-m. Counting adds and multiplies in the computation $C = A B$ , we find that conventional matrix multiplication involves $( 2 m ) ^ { 3 }$ multiplies and $( 2 m ) ^ { 3 } - ( 2 m ) ^ { 2 }$ adds. In contrast, if Strassen’s algorithm is applied with conventional multiplication at the block level, then $7 m ^ { 3 }$ multiplies and $7 m ^ { 3 } + 1 1 m ^ { 2 }$ adds are required. If $m \gg 1$ , then the Strassen method involves about $7 / 8$ the arithmetic of the fully conventional algorithm.

Now recognize that we can recur on the Strassen idea. In particular, we can apply the Strassen algorithm to each of the half-sized block multiplications associated with the $P _ { i }$ . Thus, if the original A and B are n-by-n and $n = 2 ^ { q }$ , then we can repeatedly apply the Strassen multiplication algorithm. At the bottom “level,” the blocks are 1-by-1.

Of course, there is no need to recur down to the $n = 1$ level. When the block size gets sufficiently small, $( n \leq n _ { \mathrm { m i n } } )$ , it may be sensible to use conventional matrix multiplication when finding the $P _ { i }$ . Here is the overall procedure:

Algorithm 1.3.1 (Strassen Matrix Multiplication) Suppose $n = 2 ^ { q }$ and that $A \in \mathbb { R } ^ { n \times n }$ and $B \in \mathbb { R } ^ { n \times n }$ . If $n _ { \mathrm { m i n } } = 2 ^ { d }$ with d $\leq q$ , then this algorithm computes $C = A B$ by applying Strassen procedure recursively.

function $C = { \sf s t r a s s } ( A , B , n , n _ { \sf m i n } )$

if n ≤ nmin

$$
C = A B \quad (\text { conventionally   computed })
$$

else

$$
m = n / 2; u = 1: m; v = m + 1: n
$$

$$
P _ {1} = \operatorname{strass} (A (u, u) + A (v, v), B (u, u) + B (v, v), m, n _ {\min})
$$

$$
P _ {2} = \operatorname{strass} (A (v, u) + A (v, v), B (u, u), m, n _ {\min})
$$

$$
P _ {3} = \operatorname{strass} (A (u, u), B (u, v) - B (v, v), m, n _ {\min})
$$

$$
P _ {4} = \operatorname{strass} (A (v, v), B (v, u) - B (u, u), m, n _ {\min})
$$

$$
P _ {5} = \operatorname{strass} (A (u, u) + A (u, v), B (v, v), m, n _ {\min})
$$

$$
P _ {6} = \operatorname{strass} (A (v, u) - A (u, u), B (u, u) + B (u, v), m, n _ {\min})
$$

$$
P _ {7} = \operatorname{strass} (A (u, v) - A (v, v), B (v, u) + B (v, v), m, n _ {\min})
$$

$$
C (u, u) = P _ {1} + P _ {4} - P _ {5} + P _ {7}
$$

$$
C (u, v) = P _ {3} + P _ {5}
$$

$$
C (v, u) = P _ {2} + P _ {4}
$$

$$
C (v, v) = P _ {1} + P _ {3} - P _ {2} + P _ {6}
$$

end

Unlike any of our previous algorithms, strass is recursive. Divide and conquer algorithms are often best described in this fashion. We have presented strass in the style of a Matlab function so that the recursive calls can be stated with precision.

The amount of arithmetic associated with strass is a complicated function of n and $n _ { \mathrm { m i n } } . \mathrm { ~ I f ~ } n _ { \mathrm { m i n } } \gg 1$ , then it suffices to count multiplications as the number of additions is roughly the same. If we just count the multiplications, then it suffices to examine the deepest level of the recursion as that is where all the multiplications occur. In strass there are $q - d$ subdivisions and thus $7 ^ { q - d }$ conventional matrix-matrix multiplications to perform. These multiplications have size $n _ { \mathrm { m i n } }$ and thus strass involves about $s =$ $( 2 ^ { d } ) ^ { 3 } 7 ^ { q - d }$ multiplications compared to $c = ~ ( 2 ^ { q } ) ^ { 3 }$ , the number of multiplications in the conventional approach. Notice that

$$
{\frac {s}{c}} = \left(\frac {2 ^ {d}}{2 ^ {q}}\right) ^ {3} 7 ^ {q - d} = \left(\frac {7}{8}\right) ^ {q - d}.
$$

If $d = 0$ , i.e., we recur on down to the 1-by-1 level, then

$$
s = (7 / 8) ^ {q} c = 7 ^ {q} = n ^ {\log_ {2} 7} \approx n ^ {2. 8 0 7}.
$$

Thus, asymptotically, the number of multiplications in Strassen’s method is $O ( n ^ { 2 . 8 0 7 } )$ . However, the number of additions (relative to the number of multiplications) becomes significant as $n _ { \mathrm { m i n } }$ gets small.

# Problems

P1.3.1 Rigorously prove the following block matrix equation:

$$
\left[ \begin{array}{c c c} A _ {1 1} & \dots & A _ {1 r} \\ \vdots & \ddots & \vdots \\ A _ {q 1} & \dots & A _ {q r} \end{array} \right] ^ {T} = \left[ \begin{array}{c c c} A _ {1 1} ^ {T} & \dots & A _ {q 1} ^ {T} \\ \vdots & \ddots & \vdots \\ A _ {1 r} ^ {T} & \dots & A _ {q r} ^ {T} \end{array} \right].
$$

P1.3.2 Suppose $M \in \mathbb { R } ^ { n \times n }$ is Hamiltonian. How many flops are required to compute $N = M ^ { 2 } ?$

P1.3.3 What can you say about the 2-by-2 block structure of a matrix $A \in \mathbb { R } ^ { 2 n \times 2 n }$ that satisfies $\mathcal { E } _ { 2 n } A \mathcal { E } _ { 2 n } = A ^ { T }$ where $\mathcal { E } _ { 2 n }$ is the exchange permutation defined in §1.2.11. Explain why A is symmetric about the “antidiagonal” that extends from the (2n, 1) entry to the (1, 2n) entry.

P1.3.4 Suppose

$$
A = \left[ \begin{array}{c c} 0 & B \\ B ^ {T} & 0 \end{array} \right]
$$

where $B \in \mathbb { R } ^ { n \times n }$ is upper bidiagonal. Describe the structure of $T = P A P ^ { T }$ where $P = \mathcal { P } _ { 2 , n }$ is the perfect shuffle permutation defined in §1.2.11.

P1.3.5 Show that if B and C are each permutation matrices, then $B \otimes C$ is also a permutation matrix.

P1.3.6 Verify Equation (1.3.5).

P1.3.7 Verify that if $x \in \mathbb { R } ^ { m }$ and $y \in \mathbb { R } ^ { n }$ , then $y \otimes x \ = \ \mathbf { v e c } ( x y ^ { T } )$ .

P1.3.8 Show that if $B \in \mathbb { R } ^ { p \times p } , C \in \mathbb { R } ^ { q \times q } ,$ and

$$
x = \left[ \begin{array}{c} x _ {1} \\ \vdots \\ x _ {p} \end{array} \right] \qquad x _ {i} \in \mathbb {R} ^ {q},
$$

then

$$
x ^ {T} (B \otimes C) x = \sum_ {i = 1} ^ {p} \sum_ {j = 1} ^ {p} b _ {i j} \left(x _ {i} ^ {T} C x _ {j}\right).
$$

P1.3.9 Suppose $A ^ { ( k ) } \in \mathbb { R } ^ { n _ { k } \times n _ { k } }$ for $k = 1 { : } r$ and that $\boldsymbol { x } \in \mathbb { R } ^ { n }$ where $n = n _ { 1 } \cdot \cdot \cdot n _ { r }$ . Give an efficient algorithm for computing $y = \left( A ^ { ( r ) } \otimes \cdots \otimes A ^ { ( 2 ) } \otimes A ^ { ( 1 ) } \right)$ x.

P1.3.10 Suppose n is even and define the following function from $\mathbb { R } ^ { n }$ to IR:

$$
f (x) = x (1: 2: n) ^ {T} x (2: 2: n) = \sum_ {i = 1} ^ {n / 2} x _ {2 i - 1} x _ {2 i}.
$$

(a) Show that if x, $\boldsymbol { y } \in \mathbb { R } ^ { n }$ then

$$
x ^ {T} y = \sum_ {i = 1} ^ {n / 2} (x _ {2 i - 1} + y _ {2 i}) (x _ {2 i} + y _ {2 i - 1}) - f (x) - f (y).
$$

(b) Now consider the n-by-n matrix multiplication $C = A B$ . Give an algorithm for computing this product that requires $n ^ { 3 } / 2$ multiplies once f is applied to the rows of A and the columns of B. See Winograd (1968) for details.

P1.3.12 Adapt strass so that it can handle square matrix multiplication of any order. Hint: If the “current” A has odd dimension, append a zero row and column.

P1.3.13 Adapt strass so that it can handle nonsquare products, $\mathbf { e . g . } , C = A B$ where $A \in \mathbb { R } ^ { m \times r }$ and $B \in \mathbb { R } ^ { r \times n }$ . Is it better to augment A and B with zeros so that they become square and equal in size or to “tile” A and B with square submatrices?

P1.3.14 Let $W _ { n }$ be the number of flops that strass requires to compute an $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ product where n is a power of 2. Note that $W _ { 2 } = 2 5$ and that for $n \geq 4$

$$
W _ {n} = 7 W _ {n / 2} + 1 8 (n / 2) ^ {2}
$$

Show that for every $\epsilon > 0$ there is a constant $c _ { \epsilon }$ so $W _ { n } \leq c _ { \epsilon } n ^ { \omega + \epsilon }$ where $\omega = \log _ { 2 } 7$ and n is any power of two.

P1.3.15 Suppose $B \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } } , ~ C \in \mathbb { R } ^ { m _ { 2 } \times n _ { 2 } }$ , and $D \in \mathbb { R } ^ { m _ { 3 } \times n _ { 3 } }$ . Show how to compute the vector $y = ( B \otimes C \otimes D ) x$ where $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $n = n _ { 1 } n _ { 2 } n _ { 3 }$ is given. Is the order of operations important from the flop point of view?

# Notes and References for §1.3

Useful references for the Kronecker product include Horn and Johnson (TMA, Chap. 4), Van Loan (FFT), and:

C.F. Van Loan (2000). “The Ubiquitous Kronecker Product,” J. Comput. Appl. Math., 123, 85–100.

For quite some time fast methods for matrix multiplication have attracted a lot of attention within computer science, see:

S. Winograd (1968). “A New Algorithm for Inner Product,” IEEE Trans. Comput. C-17, 693–694.

V. Strassen (1969). “Gaussian Elimination is not Optimal,” Numer. Math. 13, 354–356.

V. Pan (1984). “How Can We Speed Up Matrix Multiplication?,” SIAM Review 26, 393–416.

I. Kaporin (1999). “A Practical Algorithm for Faster Matrix Multiplication,” Num. Lin. Alg. 6, 687–700.

H. Cohn, R. Kleinberg, B. Szegedy, and C. Umans (2005). “Group-theoretic Algorithms for Matrix Multiplication,” Proceeedings of the 2005 Conference on the Foundations of Computer Science (FOCS), 379–388.

J. Demmel, I. Dumitriu, O. Holtz, and R. Kleinberg (2007). “Fast Matrix Multiplication is Stable,” Numer. Math. 106, 199–224.

P. D’Alberto and A. Nicolau (2009). “Adaptive Winograd’s Matrix Multiplication,” ACM Trans. Math. Softw. 36, Article 3.

At first glance, many of these methods do not appear to have practical value. However, this has proven not to be the case, see:

D. Bailey (1988). “Extra High Speed Matrix Multiplication on the Cray-2,” SIAM J. Sci. Stat. Comput. 9, 603–607.

N.J. Higham (1990). “Exploiting Fast Matrix Multiplication within the Level 3 BLAS,” ACM Trans. Math. Softw. 16, 352–368.

C.C. Douglas, M. Heroux, G. Slishman, and R.M. Smith (1994). “GEMMW: A Portable Level 3 BLAS Winograd Variant of Strassen’s Matrix-Matrix Multiply Algorithm,” J. Comput. Phys. 110, 1–10.

Strassen’s algorithm marked the beginning of a search for the fastest possible matrix multiplication algorithm from the complexity point of view. The exponent of matrix multiplication is the smallest number ω such that, for all $\epsilon > 0 , O ( n ^ { \omega + \epsilon } )$ work suffices. The best known value of ω has decreased over the years and is currently around 2.4. It is interesting to speculate on the existence of an $O ( n ^ { 2 + \epsilon } )$ procedure.
