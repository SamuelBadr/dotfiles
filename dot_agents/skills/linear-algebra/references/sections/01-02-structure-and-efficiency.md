# 1.2 Structure and Efficiency

The efficiency of a given matrix algorithm depends upon several factors. Most obvious and what we treat in this section is the amount of required arithmetic and storage. How to reason about these important attributes is nicely illustrated by considering examples that involve triangular matrices, diagonal matrices, banded matrices, symmetric matrices, and permutation matrices. These are among the most important types of structured matrices that arise in practice, and various economies can be realized if they are involved in a calculation.

# 1.2.1 Band Matrices

A matrix is sparse if a large fraction of its entries are zero. An important special case is the band matrix. We say that $A \in \mathbb { R } ^ { m \times n }$ has lower bandwidth p if $a _ { i j } = 0$ whenever $i > j + p$ and upper bandwidth q if $j > i + q$ implies $a _ { i j } = 0$ . Here is an example of an 8-by-5 matrix that has lower bandwidth 1 and upper bandwidth 2:

$$
A = \left[ \begin{array}{c c c c c} \times & \times & \times & 0 & 0 \\ \times & \times & \times & \times & 0 \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 \end{array} \right].
$$

The ×’s designate arbitrary nonzero entries. This notation is handy to indicate the structure of a matrix and we use it extensively. Band structures that occur frequently are tabulated in Table 1.2.1.

<table><tr><td>Type of Matrix</td><td>Lower Bandwidth</td><td>Upper Bandwidth</td></tr><tr><td>Diagonal</td><td>0</td><td>0</td></tr><tr><td>Upper triangular</td><td>0</td><td>n-1</td></tr><tr><td>Lower triangular</td><td>m-1</td><td>0</td></tr><tr><td>Tridiagonal</td><td>1</td><td>1</td></tr><tr><td>Upper bidiagonal</td><td>0</td><td>1</td></tr><tr><td>Lower bidiagonal</td><td>1</td><td>0</td></tr><tr><td>Upper Hessenberg</td><td>1</td><td>n-1</td></tr><tr><td>Lower Hessenberg</td><td>m-1</td><td>1</td></tr></table>

Table 1.2.1. Band terminology for m-by-n matrices

# 1.2.2 Triangular Matrix Multiplication

To introduce band matrix “thinking” we look at the matrix multiplication update problem $C = C + A B$ where A, B, and C are each n-by-n and upper triangular. The 3-by-3 case is illuminating:

$$
A B = \left[ \begin{array}{c c c} a _ {1 1} b _ {1 1} & a _ {1 1} b _ {1 2} + a _ {1 2} b _ {2 2} & a _ {1 1} b _ {1 3} + a _ {1 2} b _ {2 3} + a _ {1 3} b _ {3 3} \\ 0 & a _ {2 2} b _ {2 2} & a _ {2 2} b _ {2 3} + a _ {2 3} b _ {3 3} \\ 0 & 0 & a _ {3 3} b _ {3 3} \end{array} \right].
$$

It suggests that the product is upper triangular and that its upper triangular entries are the result of abbreviated inner products. Indeed, since $a _ { i k } b _ { k j } = 0$ whenever $k < i$ or $j < k$ , we see that the update has the form

$$
c _ {i j} = c _ {i j} + \sum_ {k = i} ^ {j} a _ {i k} b _ {k j}
$$

for all i and j that satisfy $i \leq j$ . This yields the following algorithm:

Algorithm 1.2.1 (Triangular Matrix Multiplication) Given upper triangular matrices $A , B , C \in \mathbb { R } ^ { n \times n }$ , this algorithm overwrites C with $C + A B$ .

for i = 1:n
    for j = i:n
    for k = i:j
    C(i,j) = C(i,j) + A(i,k)·B(k,j)
    end
    end
end

# 1.2.3 The Colon Notation—Again

The dot product that the k-loop performs in Algorithm 1.2.1 can be succinctly stated if we extend the colon notation introduced in §1.1.8. If $A \in \mathbb { R } ^ { m \times n }$ and the integers $p ,$ q, and r satisfy $1 \leq p \leq q \leq n$ and $1 \leq r \leq m$ , then

$$
A (r, p: q) = \left[ a _ {r p} \mid \dots \mid a _ {r q} \right] \in \mathbb {R} ^ {1 \times (q - p + 1)}.
$$

Likewise, if $1 \leq p \leq q \leq { \mathrm { ~ } }$ m and $1 \leq c \leq n$ , then

$$
A (p: q, c) = \left[ \begin{array}{c} a _ {p c} \\ \vdots \\ a _ {q c} \end{array} \right] \in \mathbb {R} ^ {q - p + 1}.
$$

With this notation we can rewrite Algorithm 1.2.1 as

for i = 1:n
    for j = i:n $C(i,j) = C(i,j) + A(i,i:j) \cdot B(i:j,j)$ end
end

This highlights the abbreviated inner products that are computed by the innermost loop.

# 1.2.4 Assessing Work

Obviously, upper triangular matrix multiplication involves less arithmetic than full matrix multiplication. Looking at Algorithm 1.2.1, we see that $c _ { i j }$ requires $2 ( j - i + 1 )$ flops if $( i \leq j )$ . Using the approximations

$$
\sum_ {p = 1} ^ {q} p = \frac {q (q + 1)}{2} \approx \frac {q ^ {2}}{2}
$$

and

$$
\sum_ {p = 1} ^ {q} p ^ {2} = \frac {q ^ {3}}{3} + \frac {q ^ {2}}{2} + \frac {q}{6} \approx \frac {q ^ {3}}{3},
$$

we find that triangular matrix multiplication requires one-sixth the number of flops as full matrix multiplication:

$$
\sum_ {i = 1} ^ {n} \sum_ {j = i} ^ {n} 2 (j - i + 1) = \sum_ {i = 1} ^ {n} \sum_ {j = 1} ^ {n - i + 1} 2 j \approx \sum_ {i = 1} ^ {n} \frac {2 (n - i + 1) ^ {2}}{2} = \sum_ {i = 1} ^ {n} i ^ {2} \approx \frac {n ^ {3}}{3}.
$$

We throw away the low-order terms since their inclusion does not contribute to what the flop count “says.” For example, an exact flop count of Algorithm 1.2.1 reveals that precisely $n ^ { 3 } / 3 + n ^ { 2 } + 2 n / 3$ flops are involved. For large n (the typical situation of interest) we see that the exact flop count offers no insight beyond the simple $n ^ { 3 } / 3$ accounting.

Flop counting is a necessarily crude approach to the measurement of program efficiency since it ignores subscripting, memory traffic, and other overheads associated with program execution. We must not infer too much from a comparison of flop counts. We cannot conclude, for example, that triangular matrix multiplication is six times faster than full matrix multiplication. Flop counting captures just one dimension of what makes an algorithm efficient in practice. The equally relevant issues of vectorization and data locality are taken up in §1.5.

# 1.2.5 Band Storage

Suppose $A \in \mathbb { R } ^ { n \times n }$ has lower bandwidth p and upper bandwidth q and assume that p and q are much smaller than n. Such a matrix can be stored in a $( p + q + 1 ) \ – \mathrm { b y } \ – n$ array A.band with the convention that

$$
a _ {i j} = A. \text {band} (i - j + q + 1, j) \tag {1.2.1}
$$

for all $( i , j )$ that fall inside the band, $\mathrm { e . g . }$

$$
\left[ \begin{array}{c c c c c c} a _ {1 1} & a _ {1 2} & a _ {1 3} & 0 & 0 & 0 \\ a _ {2 1} & a _ {2 2} & a _ {2 3} & a _ {2 4} & 0 & 0 \\ 0 & a _ {3 2} & a _ {3 3} & a _ {3 4} & a _ {3 5} & 0 \\ 0 & 0 & a _ {4 3} & a _ {4 4} & a _ {4 5} & a _ {4 6} \\ 0 & 0 & 0 & a _ {5 4} & a _ {5 5} & a _ {5 6} \\ 0 & 0 & 0 & 0 & a _ {6 5} & a _ {6 6} \end{array} \right] \quad \Rightarrow \quad \left[ \begin{array}{c c c c c c} * & * & a _ {1 3} & a _ {2 4} & a _ {3 5} & a _ {4 6} \\ * & a _ {1 2} & a _ {2 3} & a _ {3 4} & a _ {4 5} & a _ {5 6} \\ a _ {1 1} & a _ {2 2} & a _ {3 3} & a _ {4 4} & a _ {5 5} & a _ {6 6} \\ a _ {2 1} & a _ {3 2} & a _ {4 3} & a _ {5 4} & a _ {6 5} & * \end{array} \right].
$$

Here, the $^ { 6 6 } * ^ { 9 9 }$ entries are unused. With this data structure, our column-oriented gaxpy algorithm (Algorithm 1.1.4) transforms to the following:

Algorithm 1.2.2 (Band Storage Gaxpy) Suppose $A \in \mathbb { R } ^ { n \times n }$ has lower bandwidth p and upper bandwidth $q$ and is stored in the A.band format (1.2.1). If x, $\boldsymbol { y } \in \mathbb { R } ^ { n }$ , then this algorithm overwrites y with $y + A x$ .

for $j = 1 { : } n$

$$
\alpha_ {1} = \max (1, j - q), \alpha_ {2} = \min (n, j + p)
$$

$$
\beta_ {1} = \max (1, q + 2 - j), \beta_ {2} = \beta_ {1} + \alpha_ {2} - \alpha_ {1}
$$

$$
y \left(\alpha_ {1}: \alpha_ {2}\right) = y \left(\alpha_ {1}: \alpha_ {2}\right) + A. b a n d \left(\beta_ {1}: \beta_ {2}, j\right) x (j)
$$

end

Notice that by storing A column by column in A.band, we obtain a column-oriented saxpy procedure. Indeed, Algorithm 1.2.2 is derived from Algorithm 1.1.4 by recognizing that each saxpy involves a vector with a small number of nonzeros. Integer arithmetic is used to identify the location of these nonzeros. As a result of this careful zero/nonzero analysis, the algorithm involves just 2n $( p + q + 1 )$ flops with the assumption that p and q are much smaller than n.

# 1.2.6 Working with Diagonal Matrices

Matrices with upper and lower bandwidth zero are diagonal. If $D \in \mathbb { R } ^ { m \times n }$ is diagonal, then we use the notation

$$
D = \operatorname{diag} \left(d _ {1}, \dots , d _ {q}\right), \quad q = \min \{m, n \} \quad \Longleftrightarrow \quad d _ {i} = d _ {i i}.
$$

Shortcut notations when the dimension is clear include $\mathrm { d i a g } ( d )$ and $\mathrm { d i a g } ( d _ { i } )$ . Note that if $D = \mathrm { d i a g } ( d ) \in \mathbb { R } ^ { n \times n }$ and $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , then $D x = d . * x .$ . If $A \in \mathbb { R } ^ { m \times n }$ , then premultiplication by $D = \operatorname { d i a g } ( d _ { 1 } , \ldots , d _ { m } ) \in \mathbb { R } ^ { m \times m }$ scales rows,

$$
B = D A \quad \Longleftrightarrow \quad B (i,:) = d _ {i} \cdot A (i,:), i = 1: m
$$

while post-multiplication by $D = \mathrm { d i a g } ( d _ { 1 } , \ldots , d _ { n } ) \in \mathbb { R } ^ { n \times n }$ scales columns,

$$
B = A D \quad \Longleftrightarrow \quad B (:, j) = d _ {j} \cdot A (:, j), j = 1: n.
$$

Both of these special matrix-matrix multiplications require mn flops.

# 1.2.7 Symmetry

A matrix $A \in \mathbb { R } ^ { n \times n }$ is symmetric if $A ^ { T } = A$ and skew-symmetric if $A ^ { T } = - A$ . Likewise, a matrix $A \in \mathbb { C } ^ { n \times n }$ is Hermitian if $A ^ { H } = A$ and skew-Hermitian if $A ^ { H } = - A$ . Here are some examples:

$$
\text {Symmetric:} \qquad \left[ \begin{array}{c c c} 1 & 2 & 3 \\ 2 & 4 & 5 \\ 3 & 5 & 6 \end{array} \right], \qquad \text {Hermitian:} \qquad \left[ \begin{array}{c c c} 1 & 2 - 3 i & 4 - 5 i \\ 2 + 3 i & 6 & 7 - 8 i \\ 4 + 5 i & 7 + 8 i & 9 \end{array} \right],
$$

$$
\mathrm{Skew-Symmetric:} \left[ \begin{array}{r r r} 0 & - 2 & 3 \\ 2 & 0 & - 5 \\ - 3 & 5 & 0 \end{array} \right], \mathrm{Skew-Hermitian:} \left[ \begin{array}{r r r} i & - 2 + 3 i & - 4 + 5 i \\ 2 + 3 i & 6 i & - 7 + 8 i \\ 4 + 5 i & 7 + 8 i & 9 i \end{array} \right].
$$

For such matrices, storage requirements can be halved by simply storing the lower triangle of elements, e.g.,

$$
A = \left[ \begin{array}{l l l} 1 & 2 & 3 \\ 2 & 4 & 5 \\ 3 & 5 & 6 \end{array} \right] \quad \Leftrightarrow \quad A. v e c = \left[ \begin{array}{l l l l l l} 1 & 2 & 3 & 4 & 5 & 6 \end{array} \right].
$$

For general n, we set

$$
A. v e c ((n - j / 2) (j - 1) + i) = a _ {i j} \quad 1 \leq j \leq i \leq n. \tag {1.2.2}
$$

Here is a column-oriented gaxpy with the matrix A represented in A.vec.

Algorithm 1.2.3 (Symmetric Storage Gaxpy) Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and stored in the A.vec style (1.2.2). If $x , y \in \mathbb { R } ^ { n }$ , then this algorithm overwrites y with $y + A x$ .

for $j = 1:n$ for $i = 1:j - 1$ $y(i) = y(i) + A.vec((i - 1)n - i(i - 1)/2 + j)x(j)$ end  
    for $i = j:n$ $y(i) = y(i) + A.vec((j - 1)n - j(j - 1)/2 + i)x(j)$ end  
end

This algorithm requires the same $2 n ^ { 2 }$ flops that an ordinary gaxpy requires.

# 1.2.8 Permutation Matrices and the Identity

We denote the n-by-n identity matrix by $I _ { n } , \mathrm { e . g . }$

$$
I _ {4} = \left[ \begin{array}{c c c c} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{array} \right].
$$

We use the notation $e _ { i }$ to designate the ith column of $I _ { n }$ . If the rows of $I _ { n }$ are reordered, then the resulting matrix is said to be a permutation matrix, e.g.,

$$
P = \left[ \begin{array}{c c c c} 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 1 & 0 & 0 & 0 \end{array} \right]. \tag {1.2.3}
$$

The representation of an n-by-n permutation matrix requires just an n-vector of integers whose components specify where the 1’s occur. For example, if $v \in \mathbb { R } ^ { n }$ has the property that $v _ { i }$ specifies the column where the “1” occurs in row $i ,$ then $y = P x$ implies that $y _ { i } = x _ { v _ { i } } , \ i = 1 { : } n$ . In the example above, the underlying v-vector is $v = \left[ 2 4 3 1 \right]$ .

# 1.2.9 Specifying Integer Vectors and Submatrices

For permutation matrix work and block matrix manipulation (§1.3) it is convenient to have a method for specifying structured integer vectors of subscripts. The Matlab colon notation is again the proper vehicle and a few examples suffice to show how it works. If $n = 8$ , then

$$
\begin{array}{l} v = 1: 2: n \quad \Longrightarrow \quad v = [ 1 3 5 7 ], \\ v = n: - 1: 1 \quad \Longrightarrow \quad v = \left[ 8 7 6 5 4 3 2 1 \right], \\ v = \left[ (1: 2: n) (2: 2: n) \right] \Longrightarrow v = \left[ 1 3 5 7 2 4 6 8 \right]. \\ \end{array}
$$

Suppose $A \in \mathbb { R } ^ { m \times n }$ and that $v \in \mathbb { R } ^ { r }$ and $w \in \mathbb { R } ^ { s }$ are integer vectors with the property that $1 \leq v _ { i } \leq$ m and $1 \leq w _ { i } \leq n$ . If $B = A ( v , w )$ , then $B \in \mathbb { R } ^ { r \times s }$ is the matrix defined by $b _ { i j } = a _ { v _ { i } , w _ { j } }$ for i = 1:r and $j = 1 { : } s$ . Thus, if $A \in \mathbb { R } ^ { 8 \times 8 }$ , then

$$
A (1: 2: 8, 2: 2: 8) = \left[ \begin{array}{l l l l} a _ {1 2} & a _ {1 4} & a _ {1 6} & a _ {1 8} \\ a _ {3 2} & a _ {3 4} & a _ {3 6} & a _ {3 8} \\ a _ {5 2} & a _ {5 4} & a _ {5 6} & a _ {5 8} \\ a _ {7 2} & a _ {7 4} & a _ {7 6} & a _ {7 8} \end{array} \right].
$$

# 1.2.10 Working with Permutation Matrices

Using the colon notation, the 4-by-4 permutation matrix in (1.2.3) is defined by $P =$ $I _ { 4 } ( v , : )$ where $v ~ = ~ [ ~ 2 ~ 4 ~ 3 ~ 1 ~ ]$ . In general, if $v \in \mathbb { R } ^ { n }$ is a permutation of the vector $1 : n = [ 1 , 2 , \ldots , n ]$ and $P = I _ { n } ( v , : )$ , then

$$
\begin{array}{l} y = P x \quad \Longrightarrow \quad y = x (v) \quad \Longrightarrow \quad y _ {i} = x _ {v _ {i}}, i = 1: n \\ y = P ^ {T} x \quad \Longrightarrow \quad y (v) = x \quad \Longrightarrow \quad y _ {v _ {i}} = x _ {i}, i = 1: n \\ \end{array}
$$

The second result follows from the fact that $v _ { i }$ is the row index of the $^ { 6 6 } 1 ^ { \mathfrak { n } }$ in column i of $P ^ { T }$ . Note that $P ^ { T } ( P x ) = x$ . The inverse of a permutation matrix is its transpose.

The action of a permutation matrix on a given matrix $A \in \mathbb { R } ^ { m \times n }$ is easily described. If $P = I _ { m } ( \boldsymbol { v } , : )$ and $Q = I _ { n } ( w , : )$ , then $P A Q ^ { T } = A ( v , w )$ . It also follows that $I _ { n } ( v , : ) \cdot I _ { n } ( w , : ) = I _ { n } ( w ( v ) , : )$ . Although permutation operations involve no flops, they move data and contribute to execution time, an issue that is discussed in §1.5.

# 1.2.11 Three Famous Permutation Matrices

The exchange permutation ${ \mathcal { E } } _ { n }$ turns vectors upside down, $\mathrm { e . g . }$ ,

$$
y = \mathcal {E} _ {4} x = \left[ \begin{array}{l l l l} 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \end{array} \right] \left[ \begin{array}{l} x _ {1} \\ x _ {2} \\ x _ {3} \\ x _ {4} \end{array} \right] = \left[ \begin{array}{l} x _ {4} \\ x _ {3} \\ x _ {2} \\ x _ {1} \end{array} \right].
$$

In general, if $v = n \colon - 1 \colon 1$ , then the n-by-n exchange permutation is given by $\mathcal { E } _ { n } =$ $I _ { n } ( v , : )$ . No change results if a vector is turned upside down twice and thus, $\mathcal { E } _ { n } ^ { T } \mathcal { E } _ { n } ~ =$ $\mathcal { E } _ { n } ^ { 2 } \ = \ I _ { n }$ .

The downshift permutation $\mathcal { D } _ { n }$ pushes the components of a vector down one notch with wraparound, e.g.,

$$
y = \mathcal {D} _ {4} x = \left[ \begin{array}{l l l l} 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \end{array} \right] \left[ \begin{array}{l} x _ {1} \\ x _ {2} \\ x _ {3} \\ x _ {4} \end{array} \right] = \left[ \begin{array}{l} x _ {4} \\ x _ {1} \\ x _ {2} \\ x _ {3} \end{array} \right].
$$

In general, if $v = \left[ \left( 2 { : } n \right) 1 \right]$ , then the n-by-n downshift permutation is given by $\mathcal { D } _ { n } =$ $I _ { n } ( v , : )$ . Note that $\mathcal { D } _ { n } ^ { T }$ can be regarded as an upshift permutation.

The mod-p perfect shuffle permutation $\mathcal { P } _ { p , r }$ treats the components of the input vector $x \in \mathbb { R } ^ { n } , n = p r$ , as cards in a deck. The deck is cut into p equal “piles” and reassembled by taking one card from each pile in turn. Thus, if p = 3 and r = 4, then the piles are x(1:4), x(5:8), and x(9:12) and

$$
y = \mathcal {P} _ {3, 4} x = I _ {p r} ([ 1 5 9 2 6 1 0 3 7 1 1 4 8 1 2 ],:) x = \left[ \begin{array}{l} x (1: 4: 1 2) \\ x (2: 4: 1 2) \\ x (3: 4: 1 2) \\ x (4: 4: 1 2) \end{array} \right].
$$

In general, if $n = p r$ , then

$$
\mathcal {P} _ {p, r} = I _ {n} ([ (1: r: n) (2: r: n) \dots (r: r: n) ],:)
$$

and it can be shown that

$$
\mathcal {P} _ {p, r} ^ {T} = I _ {n} ([ (1: p: n) (2: p: n) \dots (p: p: n) ],:). \tag {1.2.4}
$$

Continuing with the card deck metaphor, $\mathcal { P } _ { p , r } ^ { T }$ reassembles the card deck by placing all the $x _ { i }$ having i mod p = 1 first, followed by all the $x _ { i }$ having i mod $p = 2$ second, and so on.

# Problems

P1.2.1 Give an algorithm that overwrites A with $A ^ { 2 }$ where $A \in \mathbb { R } ^ { n \times n }$ . How much extra storage is required? Repeat for the case when A is upper triangular.

P1.2.2 Specify an algorithm that computes the first column of the matrix $M = \left( A - \lambda _ { 1 } I \right) \cdot \cdot \cdot \left( A - \lambda _ { r } I \right)$ where $A \in \mathbb { R } ^ { n \times n }$ is upper Hessenberg and $\lambda _ { 1 } , \ldots , \lambda _ { r }$ are given scalars. How many flops are required assuming that $r \ll n ?$

P1.2.3 Give a column saxpy algorithm for the n-by-n matrix multiplication problem $C = C + A B$ where A is upper triangular and B is lower triangular.

P1.2.4 Extend Algorithm 1.2.2 so that it can handle rectangular band matrices. Be sure to describe the underlying data structure.

P1.2.5 If $A = B + i C$ is Hermitian with $B \in \mathbb { R } ^ { n \times n }$ , then it is easy to show that $B ^ { T } \ = \ B$ and $C ^ { T } = - C$ . Suppose we represent A in an array A.herm with the property that $A . h e r m ( i , j )$ houses $b _ { i j } { \mathrm { ~ i f ~ } } i \geq j$ and $c _ { i j } \ i \mathrm { f } \ j > i .$ . Using this data structure, write a matrix-vector multiply function that computes $\mathsf { R e } ( z )$ and Im(z) from Re(x) and Im(x) so that $z = A x$ .

P1.2.6 Suppose $X \in \mathbb { R } ^ { n \times p }$ and $A \in \mathbb { R } ^ { n \times n }$ are given and that A is symmetric. Give an algorithm for computing $B = X ^ { T } A X$ assuming that both A and B are to be stored using the symmetric storage scheme presented in 1.2.7.

P1.2.7 Suppose $a \in \mathbb { R } ^ { n }$ is given and that $A \in \mathbb { R } ^ { n \times n }$ has the property that $a _ { i j } = a _ { | i - j | + 1 }$ . Give an algorithm that overwrites y with $y + A x$ where x, $\boldsymbol { \mathbf { \rho } } , \boldsymbol { y } \in \mathbb { R } ^ { n }$ are given.

P1.2.8 Suppose a ∈ IRn is given and that A ∈ IRn×n has the property that aij = a((i+j−1) mod n)+1. $a \in \mathbb { R } ^ { n }$ $A \in \mathbb { R } ^ { n \times n }$ $a _ { i j } = a _ { ( ( i + j - 1 ) }$ $_ { n ) + 1 } \cdot$ Give an algorithm that overwrites y with $y + A x$ where $x , y \in \mathbb { R } ^ { n }$ are given.

P1.2.9 Develop a compact storage scheme for symmetric band matrices and write the corresponding gaxpy algorithm.

P1.2.10 Suppose $A \in \mathbb { R } ^ { n \times n } , u \in \mathbb { R } ^ { n }$ , and $v \in \mathbb { R } ^ { n }$ are given and that $k \leq n$ is an integer. Show how to compute $X \in \mathbb { R } ^ { n \times k }$ and $Y \in \mathbb { R } ^ { n \times k }$ so that $( A + u v ^ { T } ) ^ { k } = A ^ { k } + X Y ^ { T }$ . How many flops are required?

P1.2.11 Suppose $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . Write a single-loop algorithm that computes $y = \mathcal { D } _ { n } ^ { k } x$ where k is a positive integer and $\mathcal { D } _ { n }$ is defined in 1.2.11.

P1.2.12 (a) Verify (1.2.4). (b) Show that $\mathcal { P } _ { p , r } ^ { T } = \mathcal { P } _ { r , p }$

P1.2.13 The number of n-by-n permutation matrices is n!. How many of these are symmetric?

# Notes and References for §1.2

See LAPACK for a discussion about appropriate data structures when symmetry and/or bandedness is present in addition to

F.G. Gustavson (2008). “The Relevance of New Data Structure Approaches for Dense Linear Algebra in the New Multi-Core/Many-Core Environments,” in Proceedings of the 7th international Conference on Parallel Processing and Applied Mathematics, Springer-Verlag, Berlin, 618–621.

The exchange, downshift, and perfect shuffle permutations are discussed in Van Loan (FFT).
