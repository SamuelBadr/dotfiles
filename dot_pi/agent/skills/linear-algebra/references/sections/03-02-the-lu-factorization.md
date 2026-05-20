# 3.2 The LU Factorization

Triangular system solving is an easy $O ( n ^ { 2 } )$ computation. The idea behind Gaussian elimination is to convert a given system Ax = b to an equivalent triangular system. The conversion is achieved by taking appropriate linear combinations of the equations. For example, in the system

$$
3 x _ {1} + 5 x _ {2} = 9,
$$

$$
6 x _ {1} + 7 x _ {2} = 4,
$$

if we multiply the first equation by 2 and subtract it from the second we obtain

$$
3 x _ {1} + 5 x _ {2} = \quad 9,
$$

$$
- 3 x _ {2} = - 1 4.
$$

This is n = 2 Gaussian elimination. Our objective in this section is to describe the procedure in the language of matrix factorizations. This means showing that the algorithm computes a unit lower triangular matrix L and an upper triangular matrix U so that A = LU , e.g.,

$$
{\left[ \begin{array}{l l} 3 & 5 \\ 6 & 7 \end{array} \right]} = {\left[ \begin{array}{l l} 1 & 0 \\ 2 & 1 \end{array} \right]} {\left[ \begin{array}{l l} 3 & 5 \\ 0 & - 3 \end{array} \right]}.
$$

The solution to the original Ax = b problem is then found by a two-step triangular solve process:

$$
L y = b, \quad U x = y \quad \Longrightarrow \quad A x = L U x = L y = b. \tag {3.2.1}
$$

The LU factorization is a “high-level” algebraic description of Gaussian elimination. Linear equation solving is not about the matrix vector product $A ^ { - 1 } b$ but about computing LU and using it effectively; see §3.4.9. Expressing the outcome of a matrix algorithm in the “language” of matrix factorizations is a productive exercise, one that is repeated many times throughout this book. It facilitates generalization and highlights connections between algorithms that can appear very different at the scalar level.

# 3.2.1 Gauss Transformations

To obtain a factorization description of Gaussian elimination as it is traditionally presented, we need a matrix description of the zeroing process. At the $n = 2$ level, if $v _ { 1 } \neq 0$ and $\tau = v _ { 2 } / v _ { 1 }$ , then

$$
\left[ \begin{array}{c c} 1 & 0 \\ - \tau & 1 \end{array} \right] \left[ \begin{array}{c} v _ {1} \\ v _ {2} \end{array} \right] = \left[ \begin{array}{c} v _ {1} \\ 0 \end{array} \right].
$$

More generally, suppose $v \in \mathbb { R } ^ { n }$ with $v _ { k } \neq 0$ . If

$$
\tau^ {T} = [ \underbrace {0 , \ldots , 0} _ {k}, \tau_ {k + 1}, \ldots , \tau_ {n} ], \qquad \tau_ {i} = \frac {v _ {i}}{v _ {k}}, \quad i = k + 1: n,
$$

and we define

$$
M _ {k} = I _ {n} - \tau e _ {k} ^ {T}, \tag {3.2.2}
$$

then

$$
M _ {k} v = \left[ \begin{array}{c c c c c c} 1 & \dots & 0 & 0 & \dots & 0 \\ \vdots & \ddots & \vdots & \vdots & & \vdots \\ 0 & & 1 & 0 & & 0 \\ 0 & & - \tau_ {k + 1} & 1 & & 0 \\ \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & \dots & - \tau_ {n} & 0 & \dots & 1 \end{array} \right] \left[ \begin{array}{c} v _ {1} \\ \vdots \\ v _ {k} \\ v _ {k + 1} \\ \vdots \\ v _ {n} \end{array} \right] = \left[ \begin{array}{c} v _ {1} \\ \vdots \\ v _ {k} \\ 0 \\ \vdots \\ 0 \end{array} \right].
$$

A matrix of the form $M _ { k } = I _ { n } - \tau e _ { k } ^ { T } \in \mathbb { R } ^ { n \times n }$ is a Gauss transformation if the first k components of $\tau \in \mathbb { R } ^ { n }$ are zero. Such a matrix is unit lower triangular. The components of $\tau ( k + 1 { : } n )$ are called multipliers. The vector τ is called the Gauss vector.

# 3.2.2 Applying Gauss Transformations

Multiplication by a Gauss transformation is particularly simple. If $C \in \mathbb { R } ^ { n \times r }$ and $M _ { k } = I _ { n } - \tau e _ { k } ^ { T }$ is a Gauss transformation, then

$$
M _ {k} C = (I _ {n} - \tau e _ {k} ^ {T}) C = C - \tau (e _ {k} ^ {T} C) = C - \tau C (k,:) \nonumber
$$

is an outer product update. Since $\tau ( 1 { : } k ) = 0$ only $C ( k + 1 { : } n , : )$ is affected and the update $C = M _ { k } C$ can be computed row by row as follows:

for $i = k + 1 { : } n$

$$
C (i,:) = C (i,:) - \tau_ {i} \cdot C (k,:) \tag {1}
$$

end

This computation requires $2 ( n - k ) r$ flops. Here is an example:

$$
C   =   \left[ \begin{array}{c c c} 1 & 4 & 7 \\ 2 & 5 & 8 \\ 3 & 6 & 1 0 \end{array} \right],    \tau   =   \left[ \begin{array}{c} 0 \\ 1 \\ - 1 \end{array} \right] \qquad \Longrightarrow \qquad (I - \tau e _ {1} ^ {T}) C = \left[ \begin{array}{c c c} 1 & 4 & 7 \\ 1 & 1 & 1 \\ 4 & 1 0 & 1 7 \end{array} \right].
$$

# 3.2.3 Roundoff Properties of Gauss Transformations

If ˆτ is the computed version of an exact Gauss vector τ , then it is easy to verify that

$$
\hat {\tau} = \tau + e, \qquad | e | \leq \mathbf {u} | \tau |.
$$

If ˆτ is used in a Gauss transform update and $\mathrm { H } ( ( I _ { n } - \hat { \tau } e _ { k } ^ { T } ) C )$ denotes the computed result, then

$$
\mathrm{fl} \left((I _ {n} - \hat {\tau} e _ {k} ^ {T}) C\right) = (I - \tau e _ {k} ^ {T}) C + E  ,
$$

where

$$
| E | \leq 3 \mathbf {u} (| C | + | \tau | | C (k,:) |) + O (\mathbf {u} ^ {2}).
$$

Clearly, if τ has large components, then the errors in the update may be large in comparison to |C|. For this reason, care must be exercised when Gauss transformations are employed, a matter that is pursued in §3.4.

# 3.2.4 Upper Triangularizing

Assume that $A \in \mathbb { R } ^ { n \times n }$ . Gauss transformations $M _ { 1 } , \dots , M _ { n - 1 }$ can usually be found such that $M _ { n - 1 } \cdot \cdot \cdot M _ { 2 } M _ { 1 } A = U$ is upper triangular. To see this we first look at the $n = 3$ case. Suppose

$$
A = \left[ \begin{array}{l l l} 1 & 4 & 7 \\ 2 & 5 & 8 \\ 3 & 6 & 1 0 \end{array} \right]
$$

and note that

$$
M _ {1} = \left[ \begin{array}{r r r} 1 & 0 & 0 \\ - 2 & 1 & 0 \\ - 3 & 0 & 1 \end{array} \right] \quad \Rightarrow \quad M _ {1} A = \left[ \begin{array}{r r r} 1 & 4 & 7 \\ 0 & - 3 & - 6 \\ 0 & - 6 & - 1 1 \end{array} \right].
$$

Likewise, in the second step we have

$$
M _ {2} = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & - 2 & 1 \end{array} \right] \quad \Rightarrow \quad M _ {2} (M _ {1} A) = \left[ \begin{array}{c c c} 1 & 4 & 7 \\ 0 & - 3 & - 6 \\ 0 & 0 & 1 \end{array} \right].
$$

Extrapolating from this example to the general n case we conclude two things.

• At the start of the kth step we have a matrix $A ^ { ( k - 1 ) } = M _ { k - 1 } \cdot \cdot \cdot M _ { 1 } A$ that is upper triangular in columns 1 through k − 1.   
• The multipliers in the kth Gauss transform $M _ { k }$ are based on $A ^ { ( k - 1 ) } ( k + 1 { : } n , k )$ and a(k−kk $a _ { k k } ^ { ( k - 1 ) }$ must be nonzero in order to proceed.

Noting that complete upper triangularization is achieved after n − 1 steps, we obtain the following rough draft of the overall process:

$$
A ^ {(1)} = A
$$

for $k = 1 { : } n - 1$

For $i = k + 1 { : } n$ , determine the multipliers $\tau _ { i } ^ { ( k ) } = a _ { i k } ^ { ( k ) } / a _ { k k } ^ { ( k ) }$ = a ik /akk . (3.2.3)

Apply $M _ { k } = I - \tau ^ { ( k ) } e _ { k } ^ { T }$ to obtain $A ^ { ( k + 1 ) } = M _ { k } A ^ { ( k ) }$ .

end

For this process to be well-defined, the matrix entries $a _ { 1 1 } ^ { ( 1 ) } , a _ { 2 2 } ^ { ( 2 ) } , \dots , a _ { n - 1 , n - 1 } ^ { ( n - 1 ) }$ must be nonzero. These quantities are called pivots.

# 3.2.5 Existence

If no zero pivots are encountered in (3.2.3), then Gauss transformations $M _ { 1 } , \dots , M _ { n - 1 }$ are generated such that $M _ { n - 1 } \cdot \cdot \cdot M _ { 1 } A \ = \ U$ is upper triangular. It is easy to check that if $M _ { k } = I _ { n } - \tau ^ { ( k ) } e _ { k } ^ { T }$ − , then its inverse is prescribed by $\Breve { M _ { k } ^ { - 1 } } = I _ { n } + \tau ^ { ( k ) } \Breve { e _ { k } ^ { T } }$ and so

$$
A = L U \tag {3.2.4}
$$

where

$$
L = M _ {1} ^ {- 1} \dots M _ {n - 1} ^ {- 1}. \tag {3.2.5}
$$

It is clear that L is a unit lower triangular matrix because each $M _ { k } ^ { - 1 }$ is unit lower triangular. The factorization (3.2.4) is called the LU factorization.

The LU factorization may not exist. For example, it is impossible to find $l _ { i j }$ and $u _ { i j }$ so

$$
{\left[ \begin{array}{l l l} 1 & 2 & 3 \\ 2 & 4 & 7 \\ 3 & 5 & 3 \end{array} \right]} = {\left[ \begin{array}{l l l} 1 & 0 & 0 \\ \ell_ {2 1} & 1 & 0 \\ \ell_ {3 1} & \ell_ {3 2} & 1 \end{array} \right]} {\left[ \begin{array}{l l l} u _ {1 1} & u _ {1 2} & u _ {1 3} \\ 0 & u _ {2 2} & u _ {2 3} \\ 0 & 0 & u _ {3 3} \end{array} \right]}.
$$

To see this, equate entries and observe that we must have $u _ { 1 1 } = 1 , u _ { 1 2 } = 2 , \ell _ { 2 1 } = 2$ , $u _ { 2 2 } = 0$ , and $\ell _ { 3 1 } = 3$ . But then the (3,2) entry gives us the contradictory equation $5 = \ell _ { 3 1 } u _ { 1 2 } + \ell _ { 3 2 } u _ { 2 2 } = 6$ . For this example, the pivot $a _ { 2 2 } ^ { ( 1 ) } = a _ { 2 2 } - ( a _ { 2 1 } / a _ { 1 1 } ) a _ { 1 2 }$ is zero.

It turns out that the kth pivot in (3.2.3) is zero if $A ( 1 { : } k , 1 { : } k )$ is singular. A submatrix of the form A(1:k, 1:k) is called a leading principal submatrix.

Theorem 3.2.1. (LU Factorization). If $A \in \mathbb { R } ^ { n \times n }$ and det $( A ( 1 { : } k , 1 { : } k ) ) \neq 0$ for $k = 1 { : } n - 1$ , then there exists a unit lower triangular $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ and an upper triangular $U \in \mathbb { R } ^ { n \times n }$ such that $A \ = \ L U$ . If this is the case and A is nonsingular, then the factorization is unique and det $( A ) = u _ { 1 1 } \cdot \cdot \cdot u _ { n n }$ .

Proof. Suppose k − 1 steps in (3.2.3) have been executed. At the beginning of step k the matrix A has been overwritten by $M _ { k - 1 } \cdot \cdot \cdot M _ { 1 } A = A ^ { ( k - 1 ) }$ . Since Gauss transformations are unit lower triangular, it follows by looking at the leading k-by-k portion of this equation that

$$
\det (A (1: k, 1: k)) = a _ {1 1} ^ {(k - 1)} \dots a _ {k k} ^ {(k - 1)}. \tag {3.2.6}
$$

Thus, if $A ( 1 { : } k , 1 { : } k )$ is nonsingular, then the kth pivot akk $a _ { k k } ^ { ( k - 1 ) }$ is nonzero.

As for uniqueness, if $A = L _ { 1 } U _ { 1 }$ and $A = L _ { 2 } U _ { 2 }$ are two LU factorizations of a nonsingular A, then $L _ { 2 } ^ { - 1 } L _ { 1 } = U _ { 2 } U _ { 1 } ^ { - 1 }$ . Since $L _ { 2 } ^ { - 1 } L _ { 1 }$ is unit lower triangular and $U _ { 2 } U _ { 1 } ^ { - 1 }$ is upper triangular, it follows that both of these matrices must equal the identity. Hence, $L _ { 1 } = L _ { 2 }$ and $U _ { 1 } = U _ { 2 }$ . Finally, if $A = L U$ , then

$$
\det (A) = \det (L U) = \det (L) \det (U) = \det (U).
$$

It follows that det $( A ) \ = \ u _ { 1 1 } \cdot \cdot \cdot u _ { n n }$

□

# 3.2.6 L Is the Matrix of Multipliers

It turns out that the construction of L is not nearly so complicated as Equation (3.2.5) suggests. Indeed,

$$
\begin{array}{l} L = M _ {1} ^ {- 1} \dots M _ {n - 1} ^ {- 1} \\ = \left(I _ {n} - \tau^ {(1)} e _ {1} ^ {T}\right) ^ {- 1} \dots \left(I _ {n} - \tau^ {(n - 1)} e _ {n - 1} ^ {T}\right) ^ {- 1} \\ = \left(I _ {n} + \tau^ {(1)} e _ {1} ^ {T}\right) \dots \left(I _ {n} + \tau^ {(n - 1)} e _ {n - 1} ^ {T}\right) \\ = I _ {n} + \sum_ {k = 1} ^ {n - 1} \tau^ {(k)} e _ {k} ^ {T} \\ \end{array}
$$

showing that

$$
L (k + 1: n, k) = \tau^ {(k)} (k + 1: n) \quad k = 1: n - 1. \tag {3.2.7}
$$

In other words, the kth column of L is defined by the multipliers that arise in the k-th step of (3.2.3). Consider the example in §3.2.4:

$$
\tau^ {(1)} = \left[ \begin{array}{l} 0 \\ 2 \\ 3 \end{array} \right], \tau^ {(2)} = \left[ \begin{array}{l} 0 \\ 0 \\ 2 \end{array} \right] \quad \Rightarrow \quad \left[ \begin{array}{l l l} 1 & 4 & 7 \\ 2 & 5 & 8 \\ 3 & 6 & 1 0 \end{array} \right] = \left[ \begin{array}{l l l} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 3 & 2 & 1 \end{array} \right] \left[ \begin{array}{l l l} 1 & 4 & 7 \\ 0 & - 3 & - 6 \\ 0 & 0 & 1 \end{array} \right].
$$

# 3.2.7 The Outer Product Point of View

Since the application of a Gauss transformation to a matrix involves an outer product, we can regard (3.2.3) as a sequence of outer product updates. Indeed, if

$$
A = \left[ \begin{array}{c c} \alpha & w ^ {T} \\ v & B \\ 1 & n - 1 \end{array} \right] _ {n - 1} ^ {1}
$$

then the first step in Gaussian elimination results in the decomposition

$$
\left[ \begin{array}{c c} \alpha & w ^ {T} \\ z & B \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ z / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ 0 & B - z w ^ {T} / \alpha \end{array} \right] \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & I _ {n - 1} \end{array} \right].
$$

Steps 2 through n − 1 compute the LU factorization

$$
B - z w ^ {T} / \alpha = L _ {1} U _ {1}
$$

for then

$$
A = \left[ \begin{array}{l l} 1 & 0 \\ z / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{l l} 1 & 0 \\ 0 & L _ {1} U _ {1} \end{array} \right] \left[ \begin{array}{l l} \alpha & w ^ {T} \\ 0 & I _ {n - 1} \end{array} \right] = \left[ \begin{array}{l l} 1 & 0 \\ z / \alpha & L _ {1} \end{array} \right] \left[ \begin{array}{l l} \alpha & w ^ {T} \\ 0 & U _ {1} \end{array} \right] \equiv L U.
$$

# 3.2.8 Practical Implementation

Let us consider the efficient implementation of (3.2.3). First, because zeros have already been introduced in columns 1 through $k - 1$ , the Gauss transformation update need only be applied to columns k through n. Of course, we need not even apply the kth Gauss transform to $A ( : , k )$ since we know the result. So the efficient thing to do is simply to update $A ( k + 1 { : } n , k + 1 { : } n )$ . Also, the observation (3.2.7) suggests that we can overwrite $A ( k + 1 { : } n , k )$ with $L ( k + 1 { : } n , k )$ since the latter houses the multipliers that are used to zero the former. Overall we obtain:

Algorithm 3.2.1 (Outer Product LU) Suppose $A \in \mathbb { R } ^ { n \times n }$ has the property that $A ( 1 { : } k , 1 { : } k )$ is nonsingular for $k = 1 { : } n - 1$ . This algorithm computes the factorization $A = L U$ where L is unit lower triangular and U is upper triangular. For $i = 1 { : } n - 1$ , $A ( i , i { : } n )$ is overwritten by $U ( i , i ; n )$ while $A ( i + 1 { : } n , i )$ is overwritten by $L ( i + 1 { : } n , i )$ .

for $k = 1 { : } n - 1$

$$
\rho = k + 1: n
$$

$$
A (\rho , k) = A (\rho , k) / A (k, k)
$$

$$
A (\rho , \rho) = A (\rho , \rho) - A (\rho , k) \cdot A (k, \rho)
$$

end

This algorithm involves $2 n ^ { 3 } / 3$ flops and it is one of several formulations of Gaussian elimination. Note that the k-th step involves an $( n - k ) – \mathrm { b y } – ( n - k )$ outer product.

# 3.2.9 Other Versions

Similar to matrix-matrix multiplication, Gaussian elimination is a triple-loop procedure that can be arranged in several ways. Algorithm 3.2.1 corresponds to the $\ " k i j \ "$ version of Gaussian elimination if we compute the outer product update row by row:

for k = 1:n - 1
    A(k + 1:n, k) = A(k + 1:n, k)/A(k, k)
    for i = k + 1:n
    for j = k + 1:n
    A(i, j) = A(i, j) - A(i, k)·A(k, j)
    end
    end
end

There are five other versions: $k j i , i k j , i j k , j i k$ , and $j k i$ . The last of these results in an implementation that features a sequence of gaxpys and forward eliminations which we now derive at the vector level.

The plan is to compute the jth columns of L and $U$ in step j. If $j = 1$ , then by comparing the first columns in $A = L U$ we conclude that

$$
L (2: n, j) = A (2: n, 1) / A (1, 1)
$$

and $U ( 1 , 1 ) = A ( 1 , 1 )$ . Now assume that $L ( : , 1 : j - 1 )$ and $U ( 1 { : } j - 1 , 1 { : } j - 1 )$ are known. To get the jth columns of $L$ and U we equate the jth columns in the equation $A = L U$ and infer from the vector equation $A ( : , j ) = L U ( : , j )$ that

$$
A (1: j - 1, j) = L (1: j - 1, 1: j - 1) \cdot U (1: j - 1, j)
$$

and

$$
A (j: n, j) = \sum_ {k = 1} ^ {j} L (j: n, k) \cdot U (k, j).
$$

The first equation is a lower triangular linear system that can be solved for the vector $U ( 1 ; j - 1 , j )$ . Once this is accomplished, the second equation can be rearranged to produce recipes for $U ( j , j )$ and $L ( j + 1 { : } n , j )$ . Indeed, if we set

$$
\begin{array}{l} v (j: n) = A (j: n, j) - \sum_ {k = 1} ^ {j - 1} L (j: n, k) U (k, j) \\ = A (j: n, j) - L (j: n, 1: j - 1) \cdot U (1: j - 1, j), \\ \end{array}
$$

then $L ( j + 1 ; n , j ) = v ( j + 1 ; n ) / v ( j )$ and $U ( j , j ) = v ( j )$ . Thus, $L ( j + 1 { : } n , j )$ is a scaled gaxpy and we obtain the following alternative to Algorithm 3.2.1:

Algorithm 3.2.2 (Gaxpy LU) Suppose $A \in \mathbb { R } ^ { n \times n }$ has the property that $A ( 1 { : } k , 1 { : } k )$ i s nonsingular for $k = 1 { : } n - 1$ . This algorithm computes the factorization $A = L U$ where L is unit lower triangular and U is upper triangular.

Initialize L to the identity and U to the zero matrix.

for $j = 1:n$ if $j = 1$ $v = A(:,1)$ else $\tilde{a} = A(:,j)$ Solve $L(1:j-1,1:j-1) \cdot z = \tilde{a}(1:j-1)$ for $z \in \mathbb{R}^{j-1}$ . $U(1:j-1,j) = z$ $v(j:n) = \tilde{a}(j:n) - L(j:n,1:j-1) \cdot z$ end $U(j,j) = v(j)$ $L(j+1:n,j) = v(j+1:n)/v(j)$

end

(We chose to have separate arrays for L and U for clarity; it is not necessary in practice.) Algorithm 3.2.2 requires $2 n ^ { 3 } / 3$ flops, the same volume of floating point work required by Algorithm 3.2.1. However, from §1.5.2 there is less memory traffic associated with a gaxpy than with an outer product, so the two implementations could perform differently in practice. Note that in Algorithm 3.2.2, the original $A ( : , j )$ is untouched until step j.

The terms right-looking and left-looking are sometimes applied to Algorithms 3.2.1 and 3.2.2. In the outer-product implementation, after $L ( k { : } n , k )$ is determined, the columns to the right of $A ( : , k )$ are updated so it is a right-looking procedure. In contrast, subcolumns to the left of $A ( : , k )$ are accessed in gaxpy LU before $L ( k + 1 { : } n , k )$ is produced so that implementation left-looking.

# 3.2.10 The LU Factorization of a Rectangular Matrix

The LU factorization of a rectangular matrix $A \in \mathbb { R } ^ { n \times r }$ can also be performed. The $n > r$ case is illustrated by

$$
{\left[ \begin{array}{l l} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{array} \right]} = {\left[ \begin{array}{l l} 1 & 0 \\ 3 & 1 \\ 5 & 2 \end{array} \right]} {\left[ \begin{array}{l l} 1 & 2 \\ 0 & - 2 \end{array} \right]}
$$

while

$$
{\left[ \begin{array}{l l l} 1 & 2 & 3 \\ 4 & 5 & 6 \end{array} \right]} = {\left[ \begin{array}{l l} 1 & 0 \\ 4 & 1 \end{array} \right]} {\left[ \begin{array}{l l l} 1 & 2 & 3 \\ 0 & - 3 & - 6 \end{array} \right]}
$$

depicts the $n < r$ situation. The LU factorization of $A \in \mathbb { R } ^ { n \times r }$ is guaranteed to exist if $A ( 1 { : } k , 1 { : } k )$ is nonsingular for $k = 1 { : } \mathrm { m i n } \{ n , r \}$ .

The square LU factorization algorithms above needs only minor alterations to handle the rectangular case. For example, if $n > r$ , then Algorithm 3.2.1 modifies to the following:

for $k = 1 { : } r$

$$
\rho = k + 1: n
$$

$$
A (\rho , k) = A (\rho , k) / A (k, k)
$$

$\mathbf { i f } \ k < r$

$$
\mu = k + 1: r
$$

$$
A (\rho , \mu) = A (\rho , \mu) - A (\rho , k) \cdot A (k, \mu)
$$

end

end

This calculation requires $n r ^ { 2 } - r ^ { 3 } / 3$ flops. Upon completion, A is overwritten by the strictly lower triangular portion of $\ b { L } \in \mathbb { R } ^ { n \times r }$ and the upper triangular portion of $U \in \mathbb { R } ^ { r \times r }$ .

# 3.2.11 Block LU

It is possible to organize Gaussian elimination so that matrix multiplication becomes the dominant operation. Partition $A \in \mathbb { R } ^ { n \times n }$ as follows:

$$
\begin{array}{r l} {A =} & {\left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right] _ {n - r} ^ {r}} \\ & {r \quad n - r} \end{array}
$$

where $r$ is a blocking parameter. Suppose we compute the LU factorization

$$
\left[ \begin{array}{l} A _ {1 1} \\ A _ {2 1} \end{array} \right] = \left[ \begin{array}{l} L _ {1 1} \\ L _ {2 1} \end{array} \right] U _ {1 1}.
$$

Here, $L _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ is unit lower triangular and $U _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ is upper triangular and assumed to be nonsingular. If we solve $L _ { 1 1 } U _ { 1 2 } = A _ { 1 2 }$ for $U _ { 1 2 } \in \mathbb { R } ^ { r \times n - r }$ , then

$$
\left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right] = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & I _ {n - r} \end{array} \right] \left[ \begin{array}{c c} I _ {r} & 0 \\ 0 & \tilde {A} \end{array} \right] \left[ \begin{array}{c c} U _ {1 1} & U _ {1 2} \\ 0 & I _ {n - r} \end{array} \right],
$$

where

$$
\tilde {A} = A _ {2 2} - L _ {2 1} U _ {1 2} = A _ {2 2} - A _ {2 1} A _ {1 1} ^ {- 1} A _ {1 2} \tag {3.2.9}
$$

is the Schur complement of $A _ { 1 1 }$ in A. Note that if

$$
\tilde {A} = L _ {2 2} U _ {2 2}
$$

is the LU factorization of ${ \tilde { A } } .$ then

$$
A = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \end{array} \right] \left[ \begin{array}{c c} U _ {1 1} & U _ {1 2} \\ 0 & U _ {2 2} \end{array} \right]
$$

is the LU factorization of A. This lays the groundwork for a recursive implementation.

Algorithm 3.2.3 (Recursive Block LU) Suppose $A \in \mathbb { R } ^ { n \times n }$ has an LU factorization and r is a positive integer. The following algorithm computes unit lower triangular $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ and upper triangular $U \in \mathbb { R } ^ { n \times n }$ so $A = L U$ .

function $[ L , U ] = { \mathsf { B l o c k L U } } ( A , n , r )$

$n \leq r$

Compute the LU factorization A = LU using (say) Algorithm 3.2.1.

else

Use (3.2.8) to compute the LU factorization $A ( : , 1 : r ) = \left[ \begin{array} { l } { L _ { 1 1 } } \\ { L _ { 2 1 } } \end{array} \right] U _ { 1 1 }$ .

Solve $L _ { 1 1 } U _ { 1 2 } = A ( 1 { : } r , r + 1 { : } n )$ for $U _ { 1 2 }$ .

$$
\tilde {A} = A (r + 1: n, r + 1: n) - L _ {2 1} U _ {1 2}
$$

$$
\left[ L _ {2 2}, U _ {2 2} \right] = \operatorname{BlockLU} (\tilde {A}, n - r, r)
$$

$$
L = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \end{array} \right], U = \left[ \begin{array}{c c} U _ {1 1} & U _ {1 2} \\ 0 & U _ {2 2} \end{array} \right]
$$

end

The following table explains where the flops come from:

<table><tr><td>Activity</td><td>Flops</td></tr><tr><td> $L_{11}, L_{21}, U_{11}$ </td><td> $nr^{2} - r^{3}/3$ </td></tr><tr><td> $U_{12}$ </td><td> $(n - r)r^{2}$ </td></tr><tr><td> $\tilde{A}$ </td><td> $2(n - r)^{2}$ </td></tr></table>

If $n \gg r ,$ , then there are a total of about $2 n ^ { 3 } / 3$ flops, the same volume of atithmetic as Algorithms 3.2.1 and 3.2.2. The vast majority of these flops are the level-3 flops associated with the production of ${ \tilde { A } } .$ .

The actual level-3 fraction, a concept developed in §3.1.5, is more easily derived from a nonrecursive implementation. Assume for clarity that $n = N r$ where N is a positive integer and that we want to compute

$$
\left[ \begin{array}{c c c} A _ {1 1} & \dots & A _ {1 N} \\ \vdots & \ddots & \vdots \\ A _ {N 1} & \dots & A _ {N N} \end{array} \right] = \left[ \begin{array}{c c c} L _ {1 1} & \dots & 0 \\ \vdots & \ddots & \vdots \\ L _ {N 1} & \dots & L _ {N N} \end{array} \right] \left[ \begin{array}{c c c} U _ {1 1} & \dots & U _ {1 N} \\ \vdots & \ddots & \vdots \\ 0 & \dots & U _ {N N} \end{array} \right] \tag {3.2.10}
$$

where all blocks are $r { \mathrm { - } } \mathrm { b y } { \mathrm { - } } r$ . Analogously to Algorithm 3.2.3 we have the following.

Algorithm 3.2.4 (Nonrecursive Block LU) Suppose $A \in \mathbb { R } ^ { n \times n }$ has an $L U$ factorization and r is a positive integer. The following algorithm computes unit lower triangular $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ and upper triangular $U \in \mathbb { R } ^ { n \times n }$ so $A = L U$ .

for $k = 1 { : } N$

Rectangular Gaussian elimination:

$$
\left[ \begin{array}{c} A _ {k k} \\ \vdots \\ A _ {N k} \end{array} \right] = \left[ \begin{array}{c} L _ {k k} \\ \vdots \\ L _ {N k} \end{array} \right] U _ {k k}
$$

Multiple right hand side solve:

$$
L _ {k k} \left[ \begin{array}{c c c c} U _ {k, k + 1} & \ldots & U _ {k N} \end{array} \right] = \left[ \begin{array}{c c c c} A _ {k, k + 1} & \ldots & A _ {k N} \end{array} \right]
$$

Level-3 updates:

$$
A _ {i j} = A _ {i j} - L _ {i k} U _ {k j}, \quad i = k + 1: N, j = k + 1: N
$$

end

Here is the flop situation during the kth pass through the loop:

<table><tr><td>Activity</td><td>Flops</td></tr><tr><td>Gaussian elimination</td><td> $(N - k + 1)r^{3} - r^{3}/3$ </td></tr><tr><td>Multiple RHS solve</td><td> $(N - k)r^{3}$ </td></tr><tr><td>Level-3 updates</td><td> $2(N - k)^{2}r^{2}$ </td></tr></table>

Summing these quantities for $k = 1 { : } N$ we find that the level-3 fraction is approximately

$$
{\frac {2 n ^ {3} / 3}{2 n ^ {3} / 3 + n ^ {2} r}} = 1 - {\frac {3}{2 N}}.
$$

Thus, for large N almost all arithmetic takes place in the context of matrix multiplication. This ensures a favorable amount of data reuse as discussed in §1.5.4.

# Problems

P3.2.1 Verify Equation (3.2.6).

P3.2.2 Suppose the entries of $A ( \epsilon ) \in \mathbb { R } ^ { n \times n }$ are continuously differentiable functions of the scalar -. Assume that $A \equiv A ( 0 )$ and all its principal submatrices are nonsingular. Show that for sufficiently small -, the matrix $A ( \epsilon )$ has an LU factorization $A ( \epsilon ) = L ( \epsilon ) U ( \epsilon )$ and that $L ( \epsilon )$ and $U ( \epsilon )$ are both continuously differentiable.

P3.2.3 Suppose we partition $A \in \mathbb { R } ^ { n \times n }$

$$
A = \left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right]
$$

where $A _ { 1 1 }$ is $r { \mathrm { - b y } } - r$ and nonsingular. Let S be the Schur complement of $A _ { 1 1 }$ in A as defined in (3.2.9). Show that after r steps of Algorithm 3.2.1, $A ( r + 1 { : } n , r + 1 { : } n )$ houses S. How could S be obtained after r steps of Algorithm 3.2.2?

P3.2.4 Suppose $A \in \mathbb { R } ^ { n \times n }$ has an LU factorization. Show how Ax = b can be solved without storing the multipliers by computing the LU factorization of the n-by-(n + 1) matrix [A b].

P3.2.5 Describe a variant of Gaussian elimination that introduces zeros into the columns of A in the order, $n \colon - 1 { : } 2$ and which produces the factorization $A = U L$ where U is unit upper triangular and L is lower triangular.

P3.2.6 Matrices in $\mathbb { R } ^ { n \times n }$ of the form $N ( y , k ) = I - y e _ { k } ^ { T }$ where $\boldsymbol { y } \in \mathbb { R } ^ { n }$ are called Gauss-Jordan transformations. (a) Give a formula for $N ( y , k ) ^ { - 1 }$ assuming it exists. (b) Given $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , under what conditions can y be found so $N ( y , k ) x = e _ { k } ? \ ( \mathrm { c } )$ Give an algorithm using Gauss-Jordan transformations that overwrites A with $A ^ { - 1 }$ . What conditions on A ensure the success of your algorithm?

P3.2.7 Extend Algorithm 3.2.2 so that it can also handle the case when A has more rows than columns.

P3.2.8 Show how A can be overwritten with L and U in Algorithm 3.2.2. Give a 3-loop specification so that unit stride access prevails.

P3.2.9 Develop a version of Gaussian elimination in which the innermost of the three loops oversees a dot product.

# Notes and References for §3.2

The method of Gaussian elimination has a long and interesting history, see:

J.F. Grcar (2011). “How Ordinary Elimination Became Gaussian Elimination,” Historica Mathematica, 38, 163–218.   
J.F. Grcar (2011). “Mathematicians of Gaussian Elimination,” Notices of the AMS 58, 782–792.   
Schur complements (3.2.9) arise in many applications. For a survey of both practical and theoretical interest, see:   
R.W. Cottle (1974). “Manifestations of the Schur Complement,” Lin. Alg. Applic. 8, 189–211.   
Schur complements are known as “Gauss transforms” in some application areas. The use of Gauss-Jordan transformations (P3.2.6) is detailed in Fox (1964). See also:   
T. Dekker and W. Hoffman (1989). “Rehabilitation of the Gauss-Jordan Algorithm,” Numer. Math. 54, 591–599.   
As we mentioned, inner product versions of Gaussian elimination have been known and used for some time. The names of Crout and Doolittle are associated with these techniques, see:   
G.E. Forsythe (1960). “Crout with Pivoting,” Commun. ACM 3, 507–508.   
W.M. McKeeman (1962). “Crout with Equilibration and Iteration,” Commun. ACM. 5, 553–555.   
Loop orderings and block issues in LU computations are discussed in:   
J.J. Dongarra, F.G. Gustavson, and A. Karp (1984). “Implementing Linear Algebra Algorithms for Dense Matrices on a Vector Pipeline Machine,” SIAM Review 26, 91–112.   
J.M. Ortega (1988). “The ijk Forms of Factorization Methods I: Vector Computers,” Parallel Comput. 7, 135–147.   
D.H. Bailey, K.Lee, and H.D. Simon (1991). “Using Strassen’s Algorithm to Accelerate the Solution of Linear Systems,” J. Supercomput. 4, 357–371.   
J.W. Demmel, N.J. Higham, and R.S. Schreiber (1995). “Stability of Block LU Factorization,” Numer. Lin. Alg. Applic. 2, 173–190.   
Suppose $A = L U$ and $A + \Delta A = ( L + \Delta L ) ( U + \Delta U )$ are LU factorizations. Bounds on the perturbations ∆L and ∆U in terms of ∆A are given in:   
G.W. Stewart (1997). “On the Perturbation of LU and Cholesky Factors,” IMA J. Numer. Anal. 17, 1–6.   
X.-W. Chang and C.C. Paige (1998). “On the Sensitivity of the LU factorization,” BIT 38, 486–501.

In certain limited domains, it is possible to solve linear systems exactly using rational arithmetic. For a snapshot of the challenges, see:

P. Alfeld and D.J. Eyre (1991). “The Exact Analysis of Sparse Rectangular Linear Systems,” ACM Trans. Math. Softw. 17, 502–518.   
P. Alfeld (2000). “Bivariate Spline Spaces and Minimal Determining Sets,” J. Comput. Appl. Math. 119, 13–27.
