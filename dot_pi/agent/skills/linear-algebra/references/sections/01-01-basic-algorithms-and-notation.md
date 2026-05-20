# 1.1 Basic Algorithms and Notation

Matrix computations are built upon a hierarchy of linear algebraic operations. Dot products involve the scalar operations of addition and multiplication. Matrix-vector multiplication is made up of dot products. Matrix-matrix multiplication amounts to a collection of matrix-vector products. All of these operations can be described in algorithmic form or in the language of linear algebra. One of our goals is to show how these two styles of expression complement each other. Along the way we pick up notation and acquaint the reader with the kind of thinking that underpins the matrix computation area. The discussion revolves around the matrix multiplication problem, a computation that can be organized in several ways.

# 1.1.1 Matrix Notation

Let IR designate the set of real numbers. We denote the vector space of all m-by-n real matrices by IRm×n: $\mathbb { R } ^ { m \times n }$

$$
A \in \mathbb {R} ^ {m \times n} \quad \Longleftrightarrow \quad A = (a _ {i j}) = \left[ \begin{array}{c c c} a _ {1 1} & \dots & a _ {1 n} \\ \vdots & & \vdots \\ a _ {m 1} & \dots & a _ {m n} \end{array} \right], \quad a _ {i j} \in \mathbb {R}
$$

If a capital letter is used to denote a matrix $( \mathrm { e . g . } , A , B , \Delta )$ , then the corresponding lower case letter with subscript ij refers to the (i, j) entry $( \mathrm { e . g . } , a _ { i j } , b _ { i j } , \delta _ { i j } )$ . Sometimes we designate the elements of a matrix with the notation $[ A ] _ { i j }$ or $A ( i , j )$ .

# 1.1.2 Matrix Operations

Basic matrix operations include transposition $( \mathbb { R } ^ { m \times n } \to \mathbb { R } ^ { n \times m } )$ ,

$$
C = A ^ {T} \quad \Longrightarrow \quad c _ {i j} = a _ {j i},
$$

addition $( \mathbb { R } ^ { m \times n } \times \mathbb { R } ^ { m \times n }  \mathbb { R } ^ { m \times n } )$ ,

$$
C = A + B \quad \Longrightarrow \quad c _ {i j} = a _ {i j} + b _ {i j},
$$

scalar-matrix multiplication $( \mathbb { R } \times \mathbb { R } ^ { m \times n }  \mathbb { R } ^ { m \times n } )$ ,

$$
C = \alpha A \quad \Longrightarrow \quad c _ {i j} = \alpha a _ {i j},
$$

and matrix-matrix multiplication $( \mathbb { R } ^ { m \times p } \times \mathbb { R } ^ { p \times n } \to \mathbb { R } ^ { m \times n } )$ ,

$$
C = A B \quad \Longrightarrow \quad c _ {i j} = \sum_ {k = 1} ^ {p} a _ {i k} b _ {k j}.
$$

Pointwise matrix operations are occasionally useful, especially pointwise multiplication $( \mathbb { R } ^ { m \times n } \times \mathbb { R } ^ { m \times n }  \mathbb { R } ^ { m \times n } )$ ,

$$
C = A. * B \quad \Longrightarrow \quad c _ {i j} = a _ {i j} b _ {i j}
$$

and pointwise division $( \mathbb { R } ^ { m \times n } \times \mathbb { R } ^ { m \times n }  \mathbb { R } ^ { m \times n } )$ ,

$$
C = A. / B \quad \Longrightarrow \quad c _ {i j} = a _ {i j} / b _ {i j}.
$$

Of course, for pointwise division to make sense, the “denominator matrix” must have nonzero entries.

# 1.1.3 Vector Notation

Let $\mathbb { R } ^ { n }$ denote the vector space of real n-vectors:

$$
x \in \mathbb {R} ^ {n} \qquad \Longleftrightarrow \qquad x = \left[ \begin{array}{c} x _ {1} \\ \vdots \\ x _ {n} \end{array} \right] \quad x _ {i} \in \mathbb {R}  .
$$

We refer to $x _ { i }$ as the ith component of x. Depending upon context, the alternative notations $[ x ] _ { i }$ and $x ( i )$ are sometimes used.

Notice that we are identifying $\mathbb { R } ^ { n }$ with $\mathbb { R } ^ { n \times 1 }$ and so the members of $\mathbb { R } ^ { n }$ are column vectors. On the other hand, the elements of $\mathbb { R } ^ { 1 \times n }$ are row vectors:

$$
x \in \mathbb {R} ^ {1 \times n} \quad \Longleftrightarrow \quad x = [ x _ {1}, \dots , x _ {n} ].
$$

If x is a column vector, then $y = x ^ { T }$ is a row vector.

# 1.1.4 Vector Operations

Assume that $a \in \mathbb { R } , x \in \mathbb { R } ^ { n }$ , and $\boldsymbol { y } \in \mathbb { R } ^ { n }$ . Basic vector operations include scalar-vector multiplication,

$$
z = a x \quad \Longrightarrow \quad z _ {i} = a x _ {i},
$$

vector addition,

$$
z = x + y \quad \Longrightarrow \quad z _ {i} = x _ {i} + y _ {i},
$$

and the inner product (or dot product),

$$
c = x ^ {T} y \quad \Longrightarrow \quad c = \sum_ {i = 1} ^ {n} x _ {i} y _ {i}.
$$

A particularly important operation, which we write in update form, is the saxpy:

$$
y = a x + y \quad \Longrightarrow \quad y _ {i} = a x _ {i} + y _ {i}
$$

Here, the symbol “=” is used to denote assignment, not mathematical equality. The vector y is being updated. The name “saxpy” is used in LAPACK, a software package that implements many of the algorithms in this book. “Saxpy” is a mnemonic for “scalar a x plus y.” See LAPACK.

Pointwise vector operations are also useful, including vector multiplication,

$$
z = x. * y \quad \Longrightarrow \quad z _ {i} = x _ {i} y _ {i},
$$

and vector division,

$$
z = x. / y \quad \Longrightarrow \quad z _ {i} = x _ {i} / y _ {i}.
$$

# 1.1.5 The Computation of Dot Products and Saxpys

Algorithms in the text are expressed using a stylized version of the Matlab language. Here is our first example:

Algorithm 1.1.1 (Dot Product) If x, $\boldsymbol { y } \in \mathbb { R } ^ { n }$ , then this algorithm computes their dot product $c = x ^ { T } y$ .

$$
\begin{array}{l} c = 0 \\ \text { for } i = 1: n \\ c = c + x (i) y (i) \\ \end{array}
$$

It is clear from the summation that the dot product of two n-vectors involves n multiplications and n additions. The dot product operation is an $^ { \mathfrak { o } ( n ) \mathfrak { n } }$ operation, meaning that the amount of work scales linearly with the dimension. The saxpy computation is also $O ( n )$ :

Algorithm 1.1.2 (Saxpy) If $x , y \in \mathbb { R } ^ { n }$ and $a \in \mathbb { R }$ , then this algorithm overwrites y with $y + a x$ .

$$
\text { for } i = 1: n
$$

$$
y (i) = y (i) + a x (i)
$$

end

We stress that the algorithms in this book are encapsulations of important computational ideas and are not to be regarded as “production codes.”

# 1.1.6 Matrix-Vector Multiplication and the Gaxpy

Suppose $A \in \mathbb { R } ^ { m \times n }$ and that we wish to compute the update

$$
y = y + A x
$$

where $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $\boldsymbol { y } \in \mathbb { R } ^ { m }$ are given. This generalized saxpy operation is referred to as a gaxpy. A standard way that this computation proceeds is to update the components one-at-a-time:

$$
y _ {i} = y _ {i} + \sum_ {j = 1} ^ {n} a _ {i j} x _ {j}, \quad i = 1: m.
$$

This gives the following algorithm:

Algorithm 1.1.3 (Row-Oriented Gaxpy) If $A \in \mathbb { R } ^ { m \times n } , x \in \mathbb { R } ^ { n }$ , and $\boldsymbol { y } \in \mathbb { R } ^ { m }$ , then this algorithm overwrites y with Ax + y.

for i = 1:m
    for j = 1:n
    y(i) = y(i) + A(i, j)x(j)
    end
end

Note that this involves O(mn) work. If each dimension of A is doubled, then the amount of arithmetic increases by a factor of 4.

An alternative algorithm results if we regard Ax as a linear combination of A’s columns, e.g.,

$$
\left[ \begin{array}{l l} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{array} \right] \left[ \begin{array}{l} 7 \\ 8 \end{array} \right] = \left[ \begin{array}{l} 1 \cdot 7 + 2 \cdot 8 \\ 3 \cdot 7 + 4 \cdot 8 \\ 5 \cdot 7 + 6 \cdot 8 \end{array} \right] = 7 \left[ \begin{array}{l} 1 \\ 3 \\ 5 \end{array} \right] + 8 \left[ \begin{array}{l} 2 \\ 4 \\ 6 \end{array} \right] = \left[ \begin{array}{l} 2 3 \\ 5 3 \\ 8 3 \end{array} \right].
$$

Algorithm 1.1.4 (Column-Oriented Gaxpy) If $A \in \mathbb { R } ^ { m \times n } , x \in \mathbb { R } ^ { n }$ , and $\boldsymbol { y } \in \mathbb { R } ^ { m }$ , then this algorithm overwrites y with Ax + y.

for $j = 1:n$ for $i = 1:m$ $y(i) = y(i) + A(i,j)\cdot x(j)$ end   
end

Note that the inner loop in either gaxpy algorithm carries out a saxpy operation. The column version is derived by rethinking what matrix-vector multiplication “means” at the vector level, but it could also have been obtained simply by interchanging the order of the loops in the row version.

# 1.1.7 Partitioning a Matrix into Rows and Columns

Algorithms 1.1.3 and 1.1.4 access the data in A by row and by column, respectively. To highlight these orientations more clearly, we introduce the idea of a partitioned matrix.

From one point of view, a matrix is a stack of row vectors:

$$
A \in \mathbb {R} ^ {m \times n} \quad \Longleftrightarrow \quad A = \left[ \begin{array}{c} r _ {1} ^ {T} \\ \vdots \\ r _ {m} ^ {T} \end{array} \right], \quad r _ {k} \in \mathbb {R} ^ {n}. \tag {1.1.1}
$$

This is called a row partition of A. Thus, if we row partition

$$
\left[ \begin{array}{l l} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{array} \right],
$$

then we are choosing to think of A as a collection of rows with

$$
r _ {1} ^ {T} = \left[ \begin{array}{c c} 1 & 2 \end{array} \right], \qquad r _ {2} ^ {T} = \left[ \begin{array}{c c} 3 & 4 \end{array} \right], \qquad r _ {3} ^ {T} = \left[ \begin{array}{c c} 5 & 6 \end{array} \right].
$$

With the row partitioning (1.1.1), Algorithm 1.1.3 can be expressed as follows:

$$
\begin{array}{l} \text { for } i = 1: m \\ y _ {i} = y _ {i} + r _ {i} ^ {T} x \\ \end{array}
$$

Alternatively, a matrix is a collection of column vectors:

$$
A \in \mathbb {R} ^ {m \times n} \quad \Longleftrightarrow \quad A = \left[ c _ {1} \mid \dots \mid c _ {n} \right], \quad c _ {k} \in \mathbb {R} ^ {m}. \tag {1.1.2}
$$

We refer to this as a column partition of A. In the 3-by-2 example above, we thus would set c1 and c2 to be the first and second columns of A, respectively:

$$
c _ {1} = \left[ \begin{array}{l} 1 \\ 3 \\ 5 \end{array} \right], \qquad c _ {2} = \left[ \begin{array}{l} 2 \\ 4 \\ 6 \end{array} \right].
$$

With (1.1.2) we see that Algorithm 1.1.4 is a saxpy procedure that accesses A by columns:

$$
\text { for } j = 1: n
$$

$$
y = y + x _ {j} c _ {j}
$$

end

In this formulation, we appreciate y as a running vector sum that undergoes repeated saxpy updates.

# 1.1.8 The Colon Notation

A handy way to specify a column or row of a matrix is with the “colon” notation. If $A \in \mathbb { R } ^ { m \times n }$ , then $A ( k , : )$ designates the kth row, i.e.,

$$
A (k,:) = [ a _ {k 1}, \dots , a _ {k n} ].
$$

The kth column is specified by

$$
A (:, k) = \left[ \begin{array}{c} a _ {1 k} \\ \vdots \\ a _ {m k} \end{array} \right].
$$

With these conventions we can rewrite Algorithms 1.1.3 and 1.1.4 as

for i = 1:m

$$
y (i) = y (i) + A (i,:) \cdot x
$$

end

and

for j = 1:n

$$
y = y + x (j) \cdot A (:, j)
$$

end

respectively. By using the colon notation, we are able to suppress inner loop details and encourage vector-level thinking.

# 1.1.9 The Outer Product Update

As a preliminary application of the colon notation, we use it to understand the outer product update

$$
A = A + x y ^ {T}, \qquad A \in \mathbb {R} ^ {m \times n}, x \in \mathbb {R} ^ {m}, y \in \mathbb {R} ^ {n}.
$$

The outer product operation xyT “looks funny” but is perfectly legal, e.g.,

$$
\left[ \begin{array}{l} 1 \\ 2 \\ 3 \end{array} \right] \left[ \begin{array}{l l} 4 & 5 \end{array} \right] = \left[ \begin{array}{l l} 4 & 5 \\ 8 & 1 0 \\ 1 2 & 1 5 \end{array} \right].
$$

This is because xyT is the product of two “skinny” matrices and the number of columns in the left matrix x equals the number of rows in the right matrix yT . The entries in the outer product update are prescribed by

for i = 1:m

for j = 1:n

$$
a _ {i j} = a _ {i j} + x _ {i} y _ {j}
$$

end

end

This involves O(mn) arithmetic operations. The mission of the j loop is to add a multiple of yT to the ith row of A, i.e.,

for i = 1:m

$$
A (i,:) = A (i,:) + x (i) \cdot y ^ {T}
$$

end

On the other hand, if we make the i-loop the inner loop, then its task is to add a multiple of x to the jth column of A:

$$
\begin{array}{l} \text { for } j = 1: n \\ A (:, j) = A (:, j) + y (j) \cdot x \\ \end{array}
$$

Note that both implementations amount to a set of saxpy computations.

# 1.1.10 Matrix-Matrix Multiplication

Consider the 2-by-2 matrix-matrix multiplication problem. In the dot product formulation, each entry is computed as a dot product:

$$
\left[ \begin{array}{c c} 1 & 2 \\ 3 & 4 \end{array} \right] \left[ \begin{array}{c c} 5 & 6 \\ 7 & 8 \end{array} \right] = \left[ \begin{array}{c c} 1 \cdot 5 + 2 \cdot 7 & 1 \cdot 6 + 2 \cdot 8 \\ 3 \cdot 5 + 4 \cdot 7 & 3 \cdot 6 + 4 \cdot 8 \end{array} \right].
$$

In the saxpy version, each column in the product is regarded as a linear combination of left-matrix columns:

$$
\left[ \begin{array}{c c} 1 & 2 \\ 3 & 4 \end{array} \right] \left[ \begin{array}{c c} 5 & 6 \\ 7 & 8 \end{array} \right] = \left[ \begin{array}{c} 5 \left[ \begin{array}{c} 1 \\ 3 \end{array} \right] + 7 \left[ \begin{array}{c} 2 \\ 4 \end{array} \right], \quad 6 \left[ \begin{array}{c} 1 \\ 3 \end{array} \right] + 8 \left[ \begin{array}{c} 2 \\ 4 \end{array} \right] \end{array} \right].
$$

Finally, in the outer product version, the result is regarded as the sum of outer products:

$$
\left[ \begin{array}{l l} 1 & 2 \\ 3 & 4 \end{array} \right] \left[ \begin{array}{l l} 5 & 6 \\ 7 & 8 \end{array} \right] = \left[ \begin{array}{l} 1 \\ 3 \end{array} \right] \left[ \begin{array}{l l} 5 & 6 \end{array} \right] + \left[ \begin{array}{l} 2 \\ 4 \end{array} \right] \left[ \begin{array}{l l} 7 & 8 \end{array} \right].
$$

Although equivalent mathematically, it turns out that these versions of matrix multiplication can have very different levels of performance because of their memory traffic properties. This matter is pursued in §1.5. For now, it is worth detailing the various approaches to matrix multiplication because it gives us a chance to review notation and to practice thinking at different linear algebraic levels. To fix the discussion, we focus on the matrix-matrix update computation:

$$
C = C + A B, \qquad C \in \mathbb {R} ^ {m \times n}, A \in \mathbb {R} ^ {m \times r}, B \in \mathbb {R} ^ {r \times n}.
$$

The update $C = C + A B$ is considered instead of just C = AB because it is the more typical situation in practice.

# 1.1.11 Scalar-Level Specifications

The starting point is the familiar triply nested loop algorithm:

Algorithm 1.1.5 (ijk Matrix Multiplication) If $A \in \mathbb { R } ^ { m \times r } , B \in \mathbb { R } ^ { r \times n }$ , and $C \in \mathbb { R } ^ { m \times n }$ are given, then this algorithm overwrites C with $C + A B$ .

$$
\begin{array}{l} \text { for } i = 1: m \\ \text { for } j = 1: n \\ \text { for } k = 1: r \\ C (i, j) = C (i, j) + A (i, k) \cdot B (k, j) \\ \begin{array}{c} \text {end} \\ \text {end} \\ \text {end} \end{array} \\ \end{array}
$$

This computation involves O(mnr) arithmetic. If the dimensions are doubled, then work increases by a factor of 8.

Each loop index in Algorithm 1.1.5 has a particular role. (The subscript i names the row, j names the column, and k handles the dot product.) Nevertheless, the ordering of the loops is arbitrary. Here is the (mathematically equivalent) jki variant:

for $j = 1:n$ for $k = 1:r$ for $i = 1:m$ $C(i,j) = C(i,j) + A(i,k)B(k,j)$ end  
    end  
end

Altogether, there are six (= 3!) possibilities:

$$
i j k, \quad j i k, \quad i k j, \quad j k i, \quad k i j, \quad k j i.
$$

Each features an inner loop operation (dot product or saxpy) and each has its own pattern of data flow. For example, in the ijk variant, the inner loop oversees a dot product that requires access to a row of A and a column of B. The jki variant involves a saxpy that requires access to a column of C and a column of A. These attributes are summarized in Table 1.1.1 together with an interpretation of what is going on when

<table><tr><td>Loop Order</td><td>Inner Loop</td><td>Inner Two Loops</td><td>Inner Loop Data Access</td></tr><tr><td>ijk</td><td>dot</td><td>vector × matrix</td><td>A by row, B by column</td></tr><tr><td>jik</td><td>dot</td><td>matrix × vector</td><td>A by row, B by column</td></tr><tr><td>ikj</td><td>saxpy</td><td>row gaxpy</td><td>B by row, C by row</td></tr><tr><td>jki</td><td>saxpy</td><td>column gaxpy</td><td>A by column, C by column</td></tr><tr><td>kij</td><td>saxpy</td><td>row outer product</td><td>B by row, C by row</td></tr><tr><td>kji</td><td>saxpy</td><td>column outer product</td><td>A by column, C by column</td></tr></table>

Table 1.1.1. Matrix multiplication: loop orderings and properties

the middle and inner loops are considered together. Each variant involves the same amount of arithmetic, but accesses the A, B, and C data differently. The ramifications of this are discussed in §1.5.

# 1.1.12 A Dot Product Formulation

The usual matrix multiplication procedure regards A·B as an array of dot products to be computed one at a time in left-to-right, top-to-bottom order. This is the idea behind Algorithm 1.1.5 which we rewrite using the colon notation to highlight the mission of the innermost loop:

Algorithm 1.1.6 (Dot Product Matrix Multiplication) If $A \in \mathbb { R } ^ { m \times r } , B \in \mathbb { R } ^ { r \times n }$ , and $C \in \mathbb { R } ^ { m \times n }$ are given, then this algorithm overwrites C with $C + A B$ .

for i = 1:m
    for j = 1:n $C(i,j) = C(i,j) + A(i,:) \cdot B(:,j)$ end
end

In the language of partitioned matrices, if

$$
A = \left[ \begin{array}{c} a _ {1} ^ {T} \\ \vdots \\ a _ {m} ^ {T} \end{array} \right] \qquad \text { and } \qquad B = \left[ b _ {1} \mid \dots \mid b _ {n} \right],
$$

then Algorithm 1.1.6 has this interpretation:

for $i = 1:m$ for $j = 1:n$ $c_{ij} = c_{ij} + a_i^T b_j$ end   
end

Note that the purpose of the j-loop is to compute the ith row of the update. To emphasize this we could write

for $i = 1:m$ $c_{i}^{T} = c_{i}^{T} + a_{i}^{T}B$ end

where

$$
C = \left[ \begin{array}{c} c _ {1} ^ {T} \\ \vdots \\ c _ {m} ^ {T} \end{array} \right]
$$

is a row partitioning of C. To say the same thing with the colon notation we write

for i = 1:m $C(i,:) = C(i,:) + A(i,:) \cdot B$ end

Either way we see that the inner two loops of the ijk variant define a transposed gaxpy operation.

# 1.1.13 A Saxpy Formulation

Suppose A and C are column-partitioned as follows:

$$
A = \left[ a _ {1} \mid \dots \mid a _ {r} \right], \quad C = \left[ c _ {1} \mid \dots \mid c _ {n} \right].
$$

By comparing jth columns in $C = C + A B$ we see that

$$
c _ {j} = c _ {j} + \sum_ {k = 1} ^ {r} a _ {k} b _ {k j}, \quad j = 1: n.
$$

These vector sums can be put together with a sequence of saxpy updates.

Algorithm 1.1.7 (Saxpy Matrix Multiplication) If the matrices $A \in \mathbb { R } ^ { m \times r } , B \in \mathbb { R } ^ { r \times n }$ , and $\boldsymbol { C } \in \mathbb { R } ^ { m \times n }$ are given, then this algorithm overwrites C with $C + A B$ .

for $j = 1:n$ for $k = 1:r$ $C(:,j) = C(:,j) + A(:,k) \cdot B(k,j)$ end  
end

Note that the k-loop oversees a gaxpy operation:

for $j = 1:n$ $C(:,j) = C(:,j) + AB(:,j)$ end

# 1.1.14 An Outer Product Formulation

Consider the kij variant of Algorithm 1.1.5:

```matlab
for k = 1:r
    for j = 1:n
    for i = 1:m
    C(i, j) = C(i, j) + A(i, k)·B(k, j)
    end
    end
end 
```

The inner two loops oversee the outer product update

$$
C = C + a _ {k} b _ {k} ^ {T}
$$

where

$$
A = \left[ a _ {1} \mid \dots \mid a _ {r} \right] \quad \text { and } \quad B = \left[ \begin{array}{c} b _ {1} ^ {T} \\ \vdots \\ b _ {r} ^ {T} \end{array} \right] \tag {1.1.3}
$$

with $a _ { k } \in \mathbb { R } ^ { m }$ and $b _ { k } \in \mathbb { R } ^ { n }$ . This renders the following implementation:

Algorithm 1.1.8 (Outer Product Matrix Multiplication) If the matrices $A \in \mathbb { R } ^ { m \times r }$ , $B \in \mathbb { R } ^ { r \times n }$ , and $\dot { C } \in \mathbb { R } ^ { m \times n }$ are given, then this algorithm overwrites C with $C + A B$ .

for $k = 1:r$ $C = C + A(:,k) \cdot B(k,:)$ end

Matrix-matrix multiplication is a sum of outer products.

# 1.1.15 Flops

One way to quantify the volume of work associated with a computation is to count flops. A flop is a floating point add, subtract, multiply, or divide. The number of flops in a given matrix computation is usually obtained by summing the amount of arithmetic associated with the most deeply nested statements. For matrix-matrix multiplication, e.g., Algorithm 1.1.5, this is the 2-flop statement

$$
C (i, j) = C (i, j) + A (i, k) \cdot B (k, j).
$$

If $A \in \mathbb { R } ^ { m \times r } , B \in \mathbb { R } ^ { r \times n }$ , and $C \in \mathbb { R } ^ { m \times n }$ , then this statement is executed mnr times. Table 1.1.2 summarizes the number of flops that are required for the common operations detailed above.

<table><tr><td>Operation</td><td>Dimension</td><td>Flops</td></tr><tr><td> $\alpha = x^{T}y$ </td><td> $x, y \in \mathbb{R}^{n}$ </td><td>2n</td></tr><tr><td> $y = y + ax$ </td><td> $a \in \mathbb{R}, x, y \in \mathbb{R}^{n}$ </td><td>2n</td></tr><tr><td> $y = y + Ax$ </td><td> $A \in \mathbb{R}^{m \times n}, x \in \mathbb{R}^{n}, y \in \mathbb{R}^{m}$ </td><td>2mn</td></tr><tr><td> $A = A + yx^{T}$ </td><td> $A \in \mathbb{R}^{m \times n}, x \in \mathbb{R}^{n}, y \in \mathbb{R}^{m}$ </td><td>2mn</td></tr><tr><td> $C = C + AB$ </td><td> $A \in \mathbb{R}^{m \times r}, B \in \mathbb{R}^{r \times n}, C \in \mathbb{R}^{m \times n}$ </td><td>2mnr</td></tr></table>

Table 1.1.2. Important flop counts

# 1.1.16 Big-Oh Notation/Perspective

In certain settings it is handy to use the “Big-Oh” notation when an order-of-magnitude assessment of work suffices. (We did this in §1.1.1.) Dot products are $O ( n )$ , matrixvector products are $O ( n ^ { 2 } )$ , and matrix-matrix products are $O ( n ^ { 3 } )$ . Thus, to make efficient an algorithm that involves a mix of these operations, the focus should typically be on the highest order operations that are involved as they tend to dominate the overall computation.

# 1.1.17 The Notion of “Level” and the BLAS

The dot product and saxpy operations are examples of level-1 operations. Level-1 operations involve an amount of data and an amount of arithmetic that are linear in the dimension of the operation. An m-by-n outer product update or a gaxpy operation involves a quadratic amount of data (O(mn)) and a quadratic amount of work (O(mn)). These are level-2 operations. The matrix multiplication update $C = C + A B$ is a level-3 operation. Level-3 operations are quadratic in data and cubic in work.

Important level-1, level-2, and level-3 operations are encapsulated in the “BLAS,” an acronym that stands for Basic Linear Algebra Subprograms. See LAPACK. The design of matrix algorithms that are rich in level-3 BLAS operations is a major preoccupation of the field for reasons that have to do with data reuse (§1.5).

# 1.1.18 Verifying a Matrix Equation

In striving to understand matrix multiplication via outer products, we essentially established the matrix equation

$$
A B = \sum_ {k = 1} ^ {r} a _ {k} b _ {k} ^ {T}, \tag {1.1.4}
$$

where the $a _ { k }$ and $b _ { k }$ are defined by the partitionings in (1.1.3).

Numerous matrix equations are developed in subsequent chapters. Sometimes they are established algorithmically as above and other times they are proved at the ij-component level, e.g.,

$$
\left[ \sum_ {k = 1} ^ {r} a _ {k} b _ {k} ^ {T} \right] _ {i j} = \sum_ {k = 1} ^ {r} \left[ a _ {k} b _ {k} ^ {T} \right] _ {i j} = \sum_ {k = 1} ^ {r} a _ {i k} b _ {k j} = [ A B ] _ {i j}.
$$

Scalar-level verifications such as this usually provide little insight. However, they are sometimes the only way to proceed.

# 1.1.19 Complex Matrices

On occasion we shall be concerned with computations that involve complex matrices. The vector space of m-by-n complex matrices is designated by $\mathbb { C } ^ { m \times n }$ . The scaling, addition, and multiplication of complex matrices correspond exactly to the real case. However, transposition becomes conjugate transposition:

$$
C = A ^ {H} \quad \Longrightarrow \quad c _ {i j} = \overline {{a}} _ {j i}.
$$

The vector space of complex n-vectors is designated by $\mathbb { C } ^ { n }$ . The dot product of complex n-vectors x and y is prescribed by

$$
s = x ^ {H} y = \sum_ {i = 1} ^ {n} \overline {{x}} _ {i} y _ {i}.
$$

If $A = B + i C \in \mathbb { C } ^ { m \times n }$ , then we designate the real and imaginary parts of A by $\mathsf { R e } ( A ) =$ B and Im $( A ) = C$ , respectively. The conjugate of A is the matrix $\bar { A } = ( \bar { a } _ { i j } )$ .

# Problems

P1.1.1 Suppose $A \in \mathbb { R } ^ { n \times n }$ and $\boldsymbol { x } \in \mathbb { R } ^ { r }$ are given. Give an algorithm for computing the first column of $M = ( A - x _ { 1 } I ) \cdot \cdot \cdot ( A - x _ { r } I )$ .

P1.1.2 In a conventional 2-by-2 matrix multiplication $C = A B$ , there are eight multiplications: $a _ { 1 1 } b _ { 1 1 }$ , $a _ { 1 1 } b _ { 1 2 } , a _ { 2 1 } b _ { 1 1 }$ , a21b12, a12b21, a12b22, a22b21, and $a _ { 2 2 } b _ { 2 2 }$ . Make a table that indicates the order that these multiplications are performed for the ijk, jik, kij, ikj, jki, and $k j i$ matrix multiplication algorithms.

P1.1.3 Give an $O ( n ^ { 2 } )$ algorithm for computing $C = ( x y ^ { T } ) ^ { k }$ where x and y are n-vectors.

P1.1.4 Suppose $D = A B C$ where $A \in \mathbb { R } ^ { m \times n }$ , $B \in \mathbb { R } ^ { n \times p } ,$ , and $C \in \mathbb { R } ^ { p \times q } .$ . Compare the flop count of an algorithm that computes D via the formula $D = ( A B ) C$ versus the flop count for an algorithm that computes D using $D = A ( B C )$ . Under what conditions is the former procedure more flop-efficient than the latter?

P1.1.5 Suppose we have real n-by-n matrices $C , D , E ,$ and F . Show how to compute real n-by-n matrices A and B with just three real n-by-n matrix multiplications so that

$$
A + i B = (C + i D) (E + i F).
$$

Hint: Compute $W = ( C + D ) ( E - F )$ .

P1.1.6 Suppose $W \in \mathbb { R } ^ { n \times n }$ is defined by

$$
w _ {i j} = \sum_ {p = 1} ^ {n} \sum_ {q = 1} ^ {n} x _ {i p} y _ {p q} z _ {q j}
$$

where $X , Y , Z \in \mathbb { R } ^ { n \times n }$ . If we use this formula for each $w _ { i j }$ then it would require $O ( n ^ { 4 } )$ operations to set up W . On the other hand,

$$
w _ {i j} = \sum_ {p = 1} ^ {n} x _ {i p} \left(\sum_ {q = 1} ^ {n} y _ {p q} z _ {q j}\right) = \sum_ {p = 1} ^ {n} x _ {i p} u _ {p j}
$$

where $U = Y Z$ . Thus, $W = X U = X Y Z$ and only $O ( n ^ { 3 } )$ operations are required.

Use this methodology to develop an $O ( n ^ { 3 } )$ procedure for computing the n-by-n matrix A defined by

$$
a _ {i j} = \sum_ {k _ {1} = 1} ^ {n} \sum_ {k _ {2} = 1} ^ {n} \sum_ {k _ {3} = 1} ^ {n} E (k _ {1}, i) F (k _ {1}, i) G (k _ {2}, k _ {1}) H (k _ {2}, k _ {3}) F (k _ {2}, k _ {3}) G (k _ {3}, j)
$$

where $E , F , G , H \in \mathbb { R } ^ { n \times n }$ . Hint. Transposes and pointwise products are involved.

# Notes and References for §1.1

For an appreciation of the BLAS and their foundational role, see:

C.L. Lawson, R.J. Hanson, D.R. Kincaid, and F.T. Krogh (1979). “Basic Linear Algebra Subprograms for FORTRAN Usage,” ACM Trans. Math. Softw. 5, 308–323.   
J.J. Dongarra, J. Du Croz, S. Hammarling, and R.J. Hanson (1988). “An Extended Set of Fortran Basic Linear Algebra Subprograms,” ACM Trans. Math. Softw. 14, 1–17.   
J.J. Dongarra, J. Du Croz, I.S. Duff, and S.J. Hammarling (1990). “A Set of Level 3 Basic Linear Algebra Subprograms,” ACM Trans. Math. Softw. 16, 1–17.   
B. K˚agstr¨om, P. Ling, and C. Van Loan (1991). “High-Performance Level-3 BLAS: Sample Routines for Double Precision Real Data,” in High Performance Computing II, M. Durand and F. El Dabaghi (eds.), North-Holland, Amsterdam, 269–281.   
L.S. Blackford, J. Demmel, J. Dongarra, I. Duff, S. Hammarling, G. Henry, M. Heroux, L. Kaufman, A. Lumsdaine, A. Petitet, R. Pozo, K. Remington, and R.C. Whaley (2002). “An Updated Set of Basic Linear Algebra Subprograms (BLAS)”, ACM Trans. Math. Softw. 28, 135–151.

The order in which the operations in the matrix product $A _ { 1 } \cdots A _ { r }$ are carried out affects the flop count if the matrices vary in dimension. (See P1.1.4.) Optimization in this regard requires dynamic programming, see:

T.H. Corman, C.E. Leiserson, R.L. Rivest, and C. Stein (2001). Introduction to Algorithms, MIT Press and McGraw-Hill, 331–339.
