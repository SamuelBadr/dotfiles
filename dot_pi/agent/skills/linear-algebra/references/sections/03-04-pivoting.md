# 3.4 Pivoting

The analysis in the previous section shows that we must take steps to ensure that no large entries appear in the computed triangular factors $\hat { L }$ and $\hat { U }$ . The example

$$
A = \left[ \begin{array}{c c}. 0 0 0 1 & 1 \\ 1 & 1 \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ 1 0 0 0 0 & 1 \end{array} \right] \left[ \begin{array}{c c}. 0 0 0 1 & 1 \\ 0 & - 9 9 9 9 \end{array} \right] = L U
$$


---

<!-- golub_150_199 -->

correctly identifies the source of the difficulty: relatively small pivots. A way out of this difficulty is to interchange rows. For example, if P is the permutation

$$
P = \left[ \begin{array}{l l} 0 & 1 \\ 1 & 0 \end{array} \right]
$$

then

$$
P A = \left[ \begin{array}{c c} 1 & 1 \\ . 0 0 0 1 & 1 \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ . 0 0 0 1 & 1 \end{array} \right] \left[ \begin{array}{c c} 1 & 1 \\ 0 & . 9 9 9 9 \end{array} \right] = L U.
$$

Observe that the triangular factors have modestly sized entries.

In this section we show how to determine a permuted version of A that has a reasonably stable LU factorization. There are several ways to do this and they each corresponds to a different pivoting strategy. Partial pivoting, complete pivoting, and rook pivoting are considered. The efficient implementation of these strategies and their properties are discussed. We begin with a few comments about permutation matrices that can be used to swap rows or columns.

# 3.4.1 Interchange Permutations

The stabilizations of Gaussian elimination that are developed in this section involve data movements such as the interchange of two matrix rows. In keeping with our desire to describe all computations in “matrix terms,” we use permutation matrices to describe this process. (Now is a good time to review §1.2.8–§1.2.11.) Interchange permutations are particularly important. These are permutations obtained by swapping two rows in the identity, e.g.,

$$
\Pi = \left[ \begin{array}{c c c c} 0 & 0 & 0 & 1 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 1 & 0 & 0 & 0 \end{array} \right].
$$

Interchange permutations can be used to describe row and column swapping. If $\ b { A } \in \mathbb { R } ^ { 4 \times 4 }$ , then Π·A is A with rows 1 and 4 interchanged while A·Π is A with columns 1 and 4 swapped.

If $P = \Pi _ { m } \cdot \cdot \cdot \Pi _ { 1 }$ and each $\Pi _ { k }$ is the identity with rows k and $p i v ( k )$ interchanged, then $p i v ( 1 { : } m )$ encodes P . Indeed, $\boldsymbol { x } \in \mathbb { R } ^ { n }$ can be overwritten by $P x$ as follows:

$$
\begin{array}{l} \text { for } k = 1: m \\ x (k) \leftrightarrow x (p i v (k)) \\ \end{array}
$$

Here, the $^ { 6 6 }  ^ { \mathfrak { N } }$ notation means “swap contents.” Since each $\Pi _ { k }$ is symmetric, we have $P ^ { T } = \Pi _ { 1 } \cdot \cdot \cdot \Pi _ { m }$ . Thus, the piv representation can also be used to overwrite x with $P ^ { T } x \mathrm { : }$ :

$$
\begin{array}{l} \text { for } k = m: - 1: 1 \\ x (k) \leftrightarrow x (p i v (k)) \\ \end{array}
$$

We remind the reader that although no floating point arithmetic is involved in a permutation operation, permutations move data and have a nontrivial effect upon performance.

# 3.4.2 Partial Pivoting

Interchange permutations can be used in LU computations to guarantee that no multiplier is greater than 1 in absolute value. Suppose

$$
A = \left[ \begin{array}{c c c} 3 & 1 7 & 1 0 \\ 2 & 4 & - 2 \\ 6 & 1 8 & - 1 2 \end{array} \right].
$$

To get the smallest possible multipliers in the first Gauss transformation, we need $a _ { 1 1 }$ to be the largest entry in the first column. Thus, if $\Pi _ { 1 }$ is the interchange permutation

$$
\Pi_ {1} = \left[ \begin{array}{l l l} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{array} \right]
$$

then

$$
\Pi_ {1} A = \left[ \begin{array}{c c c} 6 & 1 8 & - 1 2 \\ 2 & 4 & - 2 \\ 3 & 1 7 & 1 0 \end{array} \right].
$$

It follows that

$$
M _ {1} = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ - 1 / 3 & 1 & 0 \\ - 1 / 2 & 0 & 1 \end{array} \right] \qquad \Longrightarrow \qquad M _ {1} \Pi_ {1} A = \left[ \begin{array}{c c c} 6 & 1 8 & - 1 2 \\ 0 & - 2 & 2 \\ 0 & 8 & 1 6 \end{array} \right].
$$

To obtain the smallest possible multiplier in $M _ { 2 }$ , we need to swap rows 2 and 3. Thus, if

$$
\Pi_ {2} = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{array} \right] \qquad \text {and} \qquad M _ {2} = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 1 / 4 & 1 \end{array} \right],
$$

then

$$
M _ {2} \Pi_ {2} M _ {1} \Pi_ {1} A = \left[ \begin{array}{c c c} 6 & 1 8 & - 1 2 \\ 0 & 8 & 1 6 \\ 0 & 0 & 6 \end{array} \right].
$$

For general n we have

for $k = 1 { : } n - 1$

Find an interchange permutation $\Pi _ { k } \in \mathbb { R } ^ { n \times n }$ that swaps

A(k, k) with the largest element in $\left| A ( k { : } n , k ) \right|$ .

$$
A = \Pi_ {k} A \tag {3.4.1}
$$

Determine the Gauss transformation $M _ { k } = I _ { n } - \tau ^ { ( k ) } e _ { k } ^ { T }$ such that if

v is the kth column of $M _ { k } A$ , then $v ( k + 1 { : } n ) = 0$ .

$$
A = M _ {k} A
$$

end

This particular row interchange strategy is called partial pivoting and upon completion, we have

$$
M _ {n - 1} \Pi_ {n - 1} \dots M _ {1} \Pi_ {1} A = U \tag {3.4.2}
$$

where U is upper triangular. As a consequence of the partial pivoting, no multiplier is larger than one in absolute value.

# 3.4.3 Where is L?

It turns out that (3.4.1) computes the factorization

$$
P A = L U \tag {3.4.3}
$$

where $P = \Pi _ { n - 1 } \cdot \cdot \cdot \Pi _ { 1 }$ , U is upper triangular, and L is unit lower triangular with $| \ell _ { i j } | \le 1$ . We show that $L ( k + 1 { : } n , k )$ is a permuted version of $M _ { k }$ ’s multipliers. From (3.4.2) it can be shown that

$$
\tilde {M} _ {n - 1} \dots \tilde {M} _ {1} P A = U \tag {3.4.4}
$$

where

$$
\tilde {M} _ {k} = (\Pi_ {n - 1} \dots \Pi_ {k + 1}) M _ {k} (\Pi_ {k + 1} \dots \Pi_ {n - 1}) \tag {3.4.5}
$$

for $k = 1 { : } n - 1$ . For example, in the $n = 4$ case we have

$$
\tilde {M} _ {3} \tilde {M} _ {2} \tilde {M} _ {1} P A = M _ {3} \cdot (\Pi_ {3} M _ {2} \Pi_ {3}) \cdot (\Pi_ {3} \Pi_ {2} M _ {1} \Pi_ {2} \Pi_ {3}) \cdot (\Pi_ {3} \Pi_ {2} \Pi_ {1}) A
$$

since the $\Pi _ { i }$ are symmetric. Moreover,

$$
\tilde {M} _ {k} = (\Pi_ {n - 1} \dots \Pi_ {k + 1}) \cdot (I _ {n} - \tau^ {(k)} e _ {k} ^ {T}) \cdot (\Pi_ {k + 1} \dots \Pi_ {n - 1}) = I _ {n} - \tilde {\tau} ^ {(k)} e _ {k} ^ {T}
$$

with $\tilde { \tau } ^ { ( k ) } = \Pi _ { n - 1 } \cdot \cdot \cdot \Pi _ { k + 1 } \tau ^ { ( k ) }$ . This shows that $\tilde { M _ { k } }$ is a Gauss transformation. The transformation from $\tau ^ { ( k ) }$ to $\tilde { \tau } ^ { ( k ) }$ is easy to implement in practice.

Algorithm 3.4.1 (Outer Product LU with Partial Pivoting) This algorithm computes the factorization $P A = L U$ where P is a permutation matrix encoded by $p i v ( 1 { : } n - 1 )$ , $L$ is unit lower triangular with $| \ell _ { i j } | \le 1$ , and U is upper triangular. For $i = 1 { : } n$ , $A ( i , i { : } n )$ is overwritten by $U ( i , i ; n )$ and $A ( i + 1 { : } n , i )$ is overwritten by $L ( i + 1 { : } n , i )$ . The permutation $P$ is given by $P = \Pi _ { n - 1 } \cdot \cdot \cdot \Pi _ { 1 }$ where $\Pi _ { k }$ is an interchange permutation obtained by swapping rows k and $p i v ( k )$ of $I _ { n }$ .

for $k = 1 { : } n - 1$

Determine $\mu$ with $k \leq \mu \leq n$ so |A(µ, k)| =  A(k:n, k) ∞

$$
p i v (k) = \mu
$$

$$
A (k,:) \leftrightarrow A (\mu ,:)
$$

$ { \mathbf { i } } \mathbf { f } \ A ( k , k ) \neq 0$

$$
\rho = k + 1: n
$$

$$
A (\rho , k) = A (\rho , k) / A (k, k)
$$

$$
A (\rho , \rho) = A (\rho , \rho) - A (\rho , k) A (k, \rho)
$$

end

end

The floating point overhead associated with partial pivoting is minimal from the standpoint of arithmetic as there are only $O ( n ^ { 2 } )$ comparisons associated with the search for the pivots. The overall algorithm involves $2 n ^ { 3 } / 3$ flops.

If Algorithm 3.4.1 is applied to

$$
A = \left[ \begin{array}{c c c} 3 & 1 7 & 1 0 \\ 2 & 4 & - 2 \\ 6 & 1 8 & - 1 2 \end{array} \right],
$$

then upon completion

$$
A = \left[ \begin{array}{c c c} 6 & 1 8 & - 1 2 \\ 1 / 2 & 8 & 1 6 \\ 1 / 3 & - 1 / 4 & 6 \end{array} \right]
$$

and $p i v = [ 3 , 3 ]$ . These two quantities encode all the information associated with the reduction:

$$
\left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{array} \right] \left[ \begin{array}{c c c} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{array} \right] A = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 1 / 2 & 1 & 0 \\ 1 / 3 & - 1 / 4 & 1 \end{array} \right] \left[ \begin{array}{c c c} 6 & 1 8 & - 1 2 \\ 0 & 8 & 1 6 \\ 0 & 0 & 6 \end{array} \right].
$$

To compute the solution to $A x \ = \ b$ after invoking Algorithm 3.4.1, we solve $L y = P b$ for y and $U x = y$ for x. Note that b can be overwritten by P b as follows

for $k = 1 { : } n - 1$

$$
b (k) \leftrightarrow b (p i v (k))
$$

end

We mention that if Algorithm 3.4.1 is applied to the problem,

$$
\left[ \begin{array}{c c}. 0 0 1 & 1. 0 0 \\ 1. 0 0 & 2. 0 0 \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right] = \left[ \begin{array}{c} 1. 0 0 \\ 3. 0 0 \end{array} \right],
$$

using 3-digit floating point arithmetic, then

$$
P   =   \left[ \begin{array}{c c} 0 & 1 \\ 1 & 0 \end{array} \right], \quad \hat {L}   =   \left[ \begin{array}{c c} 1. 0 0 & 0 \\ . 0 0 1 & 1. 0 0 \end{array} \right], \quad \hat {U}   =   \left[ \begin{array}{c c} 1. 0 0 & 2. 0 0 \\ 0 & 1. 0 0 \end{array} \right],
$$

and $\boldsymbol { \hat { x } } = [ 1 . 0 0 , . 9 9 6 ] ^ { T }$ . Recall from §3.3.2 that if Gaussian elimination without pivoting is applied to this problem, then the computed solution has O(1) error.

We mention that Algorithm 3.4.1 always runs to completion. If $A ( k { : } n , k ) = 0$ i n step k, then $M _ { k } = I _ { n }$ .

# 3.4.4 The Gaxpy Version

In §3.2 we developed outer product and gaxpy schemes for computing the LU factorization. Having just incorporated pivoting in the outer product version, it is equally straight forward to do the same with the gaxpy approach. Referring to Algorithm 3.2.2, we simply search the vector $\left| v ( j { : } n ) \right|$ in that algorithm for its maximal element and proceed accordingly.

Algorithm 3.4.2 (Gaxpy LU with Partial Pivoting) This algorithm computes the factorization $P A = L U$ where $P$ is a permutation matrix encoded by $p i v ( 1 : n - 1 )$ , $L$ is unit lower triangular with $| \ell _ { i j } | \le 1$ , and $U$ is upper triangular. For $i = 1 { : } n$ , $A ( i , i { : } n )$ is overwritten by $U ( i , i ; n )$ and $A ( i + 1 { : } n , i )$ is overwritten by $L ( i + 1 { : } n , i )$ . The permutation $P$ is given by $P = \Pi _ { n - 1 } \cdot \cdot \cdot \Pi _ { 1 }$ where $\Pi _ { k }$ is an interchange permutation obtained by swapping rows k and $p i v ( k )$ of $I _ { n }$ .

Initialize $L$ to the identity and $U$ to the zero matrix.

for j = 1:n
    if j = 1
    v = A(:, 1)
    else $\tilde{a} = \Pi_{j-1} \cdots \Pi_1 A(:, j)$ Solve $L(1:j-1, 1:j-1)z = \tilde{a}(1:j-1)$ for $z \in \mathbb{R}^{j-1}$ $U(1:j-1, j) = z$ , $v(j:n) = \tilde{a}(j:n) - L(j:n, 1:j-1) \cdot z$ end
    Determine $\mu$ with $j \leq \mu \leq n$ so $|v(\mu)| = \| v(j:n) \|_{\infty}$ and set $piv(j) = \mu$ $v(j) \leftrightarrow v(\mu)$ , $L(j, 1:j-1) \leftrightarrow L(\mu, 1:j-1)$ , $U(j, j) = v(j)$ if $v(j) \neq 0$ $L(j+1:n, j) = v(j+1:n)/v(j)$ end
end

As with Algorithm 3.4.1, this procedure requires $2 n ^ { 3 } / 3$ flops and $O ( n ^ { 2 } )$ comparisons.

# 3.4.5 Error Analysis and the Growth Factor

We now examine the stability that is obtained with partial pivoting. This requires an accounting of the rounding errors that are sustained during elimination and during the triangular system solving. Bearing in mind that there are no rounding errors associated with permutation, it is not hard to show using Theorem 3.3.2 that the computed solution $\hat { x }$ satisfies $( { \boldsymbol { A } } + { \boldsymbol { E } } ) { \hat { \boldsymbol { x } } } \ = \ b$ where

$$
| E | \leq n \mathbf {u} \left(2 | A | + 4 \hat {P} ^ {T} | \hat {L} | | \hat {U} |\right) + O (\mathbf {u} ^ {2}). \tag {3.4.6}
$$

Here we are assuming that $\hat { P } , \hat { L } .$ and $\hat { U }$ are the computed analogs of $P , L ,$ , and $U$ as produced by the above algorithms. Pivoting implies that the elements of $\hat { L }$ are bounded by one. Thus $\| { \hat { L } } \| _ { \infty } \leq n$ and we obtain the bound

$$
\| E \| _ {\infty} \leq n \mathbf {u} \left(2 \| A \| _ {\infty} + 4 n \| \hat {U} \| _ {\infty}\right) + O \left(\mathbf {u} ^ {2}\right). \tag {3.4.7}
$$

The problem now is to bound $\| \hat { U } \| _ { \infty }$ . Define the growth factor $\rho$ by

$$
\rho = \max _ {i, j, k} \frac {\left| \hat {a} _ {i j} ^ {(k)} \right|}{\| A \| _ {\infty}} \tag {3.4.8}
$$

where $\hat { A } ^ { ( k ) }$ is the computed version of the matrix $A ^ { ( k ) } = M _ { k } \Pi _ { k } \cdot \cdot \cdot M _ { 1 } \Pi _ { 1 } A$ . It follows that

$$
\| E \| _ {\infty} \leq 6 n ^ {3} \rho \| A \| _ {\infty} \mathbf {u} + O \left(\mathbf {u} ^ {2}\right). \tag {3.4.9}
$$

Whether or not this compares favorably with the ideal bound (2.7.20) hinges upon the size of the growth factor of $\rho .$ (The factor $n ^ { 3 }$ is not an operating factor in practice and may be ignored in this discussion.)

The growth factor measures how large the A-entries become during the process of elimination. Whether or not we regard Gaussian elimination with partial pivoting is safe to use depends upon what we can say about this quantity. From an average-case point of view, experiments by Trefethen and Schreiber (1990) suggest that $\rho$ is usually in the vicinity of $n ^ { 2 / 3 }$ . However, from the worst-case point of view, $\rho$ can be as large as $2 ^ { n - 1 }$ . In particular, if $A \in \mathbb { R } ^ { n \times n }$ is defined by

$$
a _ {i j} = \left\{ \begin{array}{r l} 1 & \text { if } i = j \text { or } j = n, \\ - 1 & \text { if } i > j, \\ 0 & \text { otherwise }, \end{array} \right.
$$

then there is no swapping of rows during Gaussian elimination with partial pivoting. We emerge with $A = L U$ and it can be shown that $u _ { n n } = 2 ^ { n - 1 }$ . For example,

$$
\left[ \begin{array}{r r r r} 1 & 0 & 0 & 1 \\ - 1 & 1 & 0 & 1 \\ - 1 & - 1 & 1 & 1 \\ - 1 & - 1 & - 1 & 1 \end{array} \right] = \left[ \begin{array}{r r r r} 1 & 0 & 0 & 0 \\ - 1 & 1 & 0 & 0 \\ - 1 & - 1 & 1 & 0 \\ - 1 & - 1 & - 1 & 1 \end{array} \right] \left[ \begin{array}{r r r r} 1 & 0 & 0 & 1 \\ 0 & 1 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 8 \end{array} \right].
$$

Understanding the behavior of $\rho$ requires an intuition about what makes the $U _ { - }$ factor large. Since $P A = L U$ implies $U = L ^ { - 1 } P A$ it would appear that the size of $L ^ { - 1 }$ is relevant. However, Stewart (1997) discusses why one can expect the L-factor to be well conditioned.

Although there is still more to understand about $\rho ,$ the consensus is that serious element growth in Gaussian elimination with partial pivoting is extremely rare. The method can be used with confidence.

# 3.4.6 Complete Pivoting

Another pivot strategy called complete pivoting has the property that the associated growth factor bound is considerably smaller than $2 ^ { n - 1 }$ . Recall that in partial pivoting, the kth pivot is determined by scanning the current subcolumn $A ( k { : } n , k )$ . In complete pivoting, the largest entry in the current submatrix $A ( k { : } n , k { : } n )$ is permuted into the $( k , k )$ position. Thus, we compute the upper triangularization

$$
M _ {n - 1} \Pi_ {n - 1} \dots M _ {1} \Pi_ {1} A \Gamma_ {1} \dots \Gamma_ {n - 1} = U.
$$

In step k we are confronted with the matrix

$$
A ^ {(k - 1)} = M _ {k - 1} \Pi_ {k - 1} \dots M _ {1} \Pi_ {1} A \Gamma_ {1} \dots \Gamma_ {k - 1}
$$

and determine interchange permutations $\Pi _ { k }$ and $\Gamma _ { k }$ such that

$$
\left| \left(\Pi_ {k} A ^ {(k - 1)} \Gamma_ {k}\right) _ {k k} \right| = \max _ {k \leq i, j \leq n} \left| \left(\Pi_ {k} A ^ {(k - 1)} \Gamma_ {k}\right) _ {i j} \right|.
$$

Algorithm 3.4.3 (Outer Product LU with Complete Pivoting) This algorithm computes the factorization $P A Q ^ { T } = L U$ where P is a permutation matrix encoded by $p i v ( 1 { : } n - 1 )$ , Q is a permutation matrix encoded by $c o l p i v ( 1 { : } n - 1 )$ , L is unit lower triangular with $| \ell _ { i j } | \le 1$ , and $U$ is upper triangular. For $i = 1 { : } n , A ( i , i { : } n )$ is overwritten by $U ( i , i ; n )$ and $A ( i + 1 { : } n , i )$ is overwritten by $L ( i + 1 { : } n , i )$ . The permutation $P$ is given by $P = \Pi _ { n - 1 } \cdot \cdot \cdot \Pi _ { 1 }$ where $\Pi _ { k }$ is an interchange permutation obtained by swapping rows k and rowpiv(k) of $I _ { n }$ . The permutation $Q$ is given by $Q = \Gamma _ { n - 1 } \cdot \cdot \cdot \Gamma _ { 1 }$ where $\Gamma _ { k }$ is an interchange permutation obtained by swapping rows k and colpiv(k) of $I _ { n }$ .

for k = 1:n - 1
    Determine $\mu$ with $k \leq \mu \leq n$ and $\lambda$ with $k \leq \lambda \leq n$ so $|A(\mu, \lambda)| = \max\{|A(i, j)| : i = k:n, j = k:n\}$ rowpiv(k) = $\mu$ $A(k, 1:n) \leftrightarrow A(\mu, 1:n)$ colpiv(k) = $\lambda$ $A(1:n, k) \leftrightarrow A(1:n, \lambda)$ if $A(k, k) \neq 0$ $\rho = k + 1:n$ $A(\rho, k) = A(\rho, k)/A(k, k)$ $A(\rho, \rho) = A(\rho, \rho) - A(\rho, k)A(k, \rho)$ end
end

This algorithm requires $2 n ^ { 3 } / 3$ flops and $O ( n ^ { 3 } )$ comparisons. Unlike partial pivoting, complete pivoting involves a significant floating point arithmetic overhead because of the two-dimensional search at each stage.

With the factorization $P A Q ^ { T } = L U$ in hand the solution to $A x = b$ proceeds as follows:

Step 1. Solve $L z = P b$ for z.

Step 2. Solve $U y = z \ { \mathrm { f o r } } \ y$ .

Step 3. Set $x = Q ^ { T } y .$ .

The rowpiv and colpiv representations can be used to form P b and $Q y$ , respectively.

Wilkinson (1961) has shown that in exact arithmetic the elements of the matrix $A ^ { ( k ) } = M _ { k } \Pi _ { k } \cdot \cdot \cdot M _ { 1 } \Pi _ { 1 } A \Gamma _ { 1 } \cdot \cdot \cdot \Gamma _ { k }$ satisfy

$$
\left| a _ {i j} ^ {(k)} \right| \leq k ^ {1 / 2} \left(2 \cdot 3 ^ {1 / 2} \dots k ^ {1 / k - 1}\right) ^ {1 / 2} \max \left| a _ {i j} \right|. \tag {3.4.10}
$$

The upper bound is a rather slow-growing function of k. This fact coupled with vast empirical evidence suggesting that $\rho$ is always modestly sized $( \mathrm { e . g } , \rho = 1 0 )$ permit us to conclude that Gaussian elimination with complete pivoting is stable. The method solves a nearby linear system $( { \boldsymbol { A } } + { \boldsymbol { E } } ) { \hat { \boldsymbol { x } } } = { \boldsymbol { b } }$ in the sense of (2.7.21). However, in general there is little reason to choose complete pivoting over partial pivoting. A possible exception is when A is rank deficient. In principal, complete pivoting can be used to reveal the rank of a matrix. Suppose rank $( A ) = r < n$ . It follows that at the beginning of step $r + 1$ , $A ( r + 1 { : } n , r + 1 { : } n ) = 0$ . This implies that $\Pi _ { k } = \Gamma _ { k } = M _ { k } = I$ for $k = r + 1$ :n and so the algorithm can be terminated after step r with the following factorization in hand:

$$
P A Q ^ {T}   =   L U   =   \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & I _ {n - r} \end{array} \right] \left[ \begin{array}{c c} U _ {1 1} & U _ {1 2} \\ 0 & 0 \end{array} \right]  .
$$

Here, $L _ { 1 1 }$ and $U _ { 1 1 }$ are $r { \mathrm { - } } \mathrm { b y } { \mathrm { - } } r$ and $L _ { 2 1 }$ and $U _ { 1 2 } ^ { T }$ are $( n - r ) – \mathrm { b y } – r$ . Thus, Gaussian elimination with complete pivoting can in principle be used to determine the rank of a matrix. Nevertheless, roundoff errors make the probability of encountering an exactly zero pivot remote. In practice one would have to “declare” A to have rank k if the pivot element in step $k + 1$ was sufficiently small. The numerical rank determination problem is discussed in detail in §5.5.

# 3.4.7 Rook Pivoting

A third type of LU stablization strategy called rook pivoting provides an interesting alternative to partial pivoting and complete pivoting. As with complete pivoting, it computes the factorization $P A Q \ = \ L U$ . However, instead of choosing as pivot the largest value in $\left| A ( k { : } n , k { : } n ) \right|$ , it searches for an element of that submatrix that is maximal in both its row and column. Thus, if

$$
A (k: n, k: n) = \left[ \begin{array}{c c c c} 2 4 & 3 6 & 1 3 & 6 1 \\ 4 2 & 6 7 & 7 2 & 5 0 \\ 3 8 & 1 1 & 3 6 & 4 3 \\ 5 2 & 3 7 & 4 8 & 1 6 \end{array} \right],
$$

then $^ { 6 6 } 7 2 ^ { 5 }$ would be identified by complete pivoting while $^ { 6 4 } { \bf 5 } 2 , ^ { 5 } \quad ^ { 6 } { \bf 7 } 2 , ^ { 5 }$ or $^ { 6 6 1 }$ would be acceptable with the rook pivoting strategy. To implement rook pivoting, the scanand-swap portion of Algorithm 3.4.3 is changed to

$$
\mu = k, \lambda = k, \tau = | a _ {\mu \lambda} |, s = 0
$$

$$
\text { while } \tau <   \| (A (k: n, \lambda) \| _ {\infty} \vee \tau <   \| (A (\mu , k: n) \| _ {\infty}
$$

$$
\text { if } \mod (s, 2) = 0
$$

$$
\text { Update } \mu \text { so   that } | a _ {\mu \lambda} | = \| (A (k: n, \lambda) \| _ {\infty} \text { with } k \leq \mu \leq n.
$$

$$
\text { Update } \lambda \text { so that } | a _ {\mu \lambda} | = \| (A (\mu , k: n) \| _ {\infty} \text { with } k \leq \lambda \leq n.
$$

end

$$
s = s + 1
$$

end

$$
\operatorname{rowpiv} (k) = \mu , A (k,:) \leftrightarrow A (\mu ,:) \operatorname{colpiv} (k) = \lambda , A (:, k) \leftrightarrow A (:, \lambda)
$$

The search for a larger $| a _ { \mu \lambda } |$ involves alternate scans of $A ( k { : } n , \lambda )$ and $A ( \mu , k { : } n )$ . The value of $\tau$ is monotone increasing and that ensures termination of the while-loop. In theory, the exit value of s could be $O ( n - k ) ^ { 2 } )$ , but in practice its value is $O ( 1 )$ . See Chang (2002). The bottom line is that rook pivoting represents the same $O ( n ^ { 2 } )$ overhead as partial pivoting, but that it induces the same level of reliability as complete pivoting.

# 3.4.8 A Note on Underdetermined Systems

If $A \in \mathbb { R } ^ { m \times n }$ with $m < n$ , rank $( A ) = m$ , and $b \in \mathbb { R } ^ { m }$ , then the linear system $A x = b$ is said to be underdetermined. Note that in this case there are an infinite number of solutions. With either complete or rook pivoting, it is possible to compute an LU factorization of the form

$$
P A Q ^ {T} = L \left[ U _ {1} \mid U _ {2} \right] \tag {3.4.11}
$$

where P and $Q$ are permutations, $\boldsymbol { L } \in \mathbb { R } ^ { m \times m }$ is unit lower triangular, and $U _ { 1 } \in \mathbb { R } ^ { m }$ ×m is nonsingular and upper triangular. Note that

$$
A x = b \Leftrightarrow (P A Q ^ {T}) (Q x) = (P b) \Leftrightarrow L \left[ U _ {1} \mid U _ {2} \right] \left[ \begin{array}{l} z _ {1} \\ z _ {2} \end{array} \right] = L (U _ {1} z _ {1} + U _ {2} z _ {2}) = c
$$

where $c = P b$ and

$$
\left[ \begin{array}{l} z _ {1} \\ z _ {2} \end{array} \right] = Q x.
$$

This suggests the following solution procedure:

Step 1. Solve $L y = P b$ for $\boldsymbol { y } \in \mathbb { R } ^ { m }$ .

Step 2. Choose $z _ { 2 } \in \mathbb { R } ^ { n - m }$ and solve $U _ { 1 } z _ { 1 } = y - U _ { 2 } z _ { 2 } { \mathrm { ~ f o r ~ } } z _ { 1 }$

Step 3. Set

$$
x = Q ^ {T} \left[ \begin{array}{l} z _ {1} \\ z _ {2} \end{array} \right].
$$

Setting $z _ { 2 } = 0$ is a natural choice. We have more to say about underdetermined systems in §5.6.2.

# 3.4.9 The LU Mentality

We offer three examples that illustrate how to think in terms of the LU factorization when confronted with a linear equation situation.

Example 1. Suppose A is nonsingular and n-by-n and that B is n-by-p. Consider the problem of finding $X \ ( n { \mathrm { - b y - } } p )$ so $A X = B$ . This is the multiple right hand side problem. If $X = \left[ \left. x _ { 1 } \right| \cdot \cdot \cdot \right| \left. x _ { p } \right]$ and $B = \left[ \left. b _ { 1 } \right| \cdot \cdot \cdot \right| \left. b _ { p } \right]$ are column partitions, then

Compute $P A = L U$

for k = 1:p

$\mathrm { S o l v e ~ } L y = P b _ { k } \mathrm { ~ a n d ~ t h e n ~ } U x _ { k } = y .$ (3.4.12)

end

If $B = I _ { n }$ , then we emerge with an approximation to $A ^ { - 1 }$ .

Example 2. Suppose we want to overwrite b with the solution to $A ^ { k } x = b$ where $A \in \mathbb { R } ^ { n \times n } , b \in \mathbb { R } ^ { n }$ , and k is a positive integer. One approach is to compute $C = A ^ { k }$ and then solve $C x = b$ . However, the matrix multiplications can be avoided altogether:

Compute P A = LU.

for j = 1:k

Overwrite b with the solution to Ly = P b. (3.4.13)

Overwrite b with the solution to U x = b.

end

As in Example 1, the idea is to get the LU factorization “outside the loop.”

Example 3. Suppose we are given $A \in \mathbb { R } ^ { n \times n } , \ d \in \mathbb { R } ^ { n }$ , and $c \in \mathbb { R } ^ { n }$ and that we want to compute $s = c ^ { T } A ^ { - 1 } d .$ . One approach is to compute $X = A ^ { - 1 }$ as discussed in (i) and then compute $s = c ^ { T } X d .$ However, it is more economical to proceed as follows:

Compute $P A = L U$ .

Solve Ly = P d and then $U x = y$ .

$$
s = c ^ {T} x
$$

An $^ { 6 6 } A ^ { - 1 9 }$ in a formula almost always means “solve a linear system” and almost never means “compute A−1.” $A ^ { - 1 }$

# 3.4.10 A Model Problem for Numerical Analysis

We are now in possession of a very important and well-understood algorithm (Gaussian elimination) for a very important and well-understood problem (linear equations). Let us take advantage of our position and formulate more abstractly what we mean by “problem sensitivity” and “algorithm stability.” Our discussion follows Higham (ASNA, §1.5–1.6), Stewart (MA, §4.3), and Trefethen and Bau (NLA, Lectures 12, 14, 15, and 22).

A problem is a function $f { : } D \to S$ from “data/input space” D to “solution/output space” S. A problem instance is f together with a particular $d \in D$ . We assume D and S are normed vector spaces. For linear systems, D is the set of matrix-vector pairs $( A , b )$ where $A \in \mathbb { R } ^ { n \times n }$ is nonsingular and $b \in \mathbb { R } ^ { n }$ . The function f maps $( A , b )$ to $A ^ { - 1 } b .$ , an element of S. For a particular A and $b , A x = b$ is a problem instance.

A perturbation theory for the problem f sheds light on the difference between $f ( d )$ and $f ( d + \Delta d )$ where $d \in D$ and $d + \Delta d \in D$ . For linear systems, we discussed in §2.6 the difference between the solution to $A x = b$ and the solution to $( A + \Delta A ) ( x + \Delta x ) =$ $( b + \Delta b )$ . We bounded $\| \Delta x \| / \| x \|$ in terms of $\| \Delta A \| / \| A \|$ and $\parallel \Delta b \parallel / \parallel b \parallel$ .

The conditioning of a problem refers to the behavior of f under perturbation at d. A condition number of a problem quantifies the rate of change of the solution with respect to the input data. If small changes in d induce relatively large changes in $f ( d )$ , then that problem instance is ill-conditioned. If small changes in d do not induce relatively large changes in $f ( d )$ , then that problem instance is well-conditioned. Definitions for “small” and “large” are required. For linear systems we showed in $\ S 2 . 6$ that the magnitude of the condition number $\kappa ( A ) = \| A \| \| A ^ { - 1 } \|$ determines whether an $A x = b$ problem is ill-conditioned or well-conditioned. One might say that a linear equation problem is well-conditioned if $\kappa ( A ) \approx O ( 1 )$ and ill-conditioned if $\kappa ( A ) \approx O ( 1 / \mathbf { u } )$ .

An algorithm for computing $f ( d )$ produces an approximation $\tilde { f } ( d )$ . Depending on the situation, it may be necessary to identify a particular software implementation of the underlying method. The $\tilde { f }$ function for Gaussian elimination with partial pivoting, Gaussian elimination with rook pivoting, and Gaussian elimination with complete pivoting are all different.

An algorithm for computing $f ( d )$ is stable if for some small $\Delta d ,$ the computed solution $\tilde { f } ( d )$ is close to $f ( d + \Delta d )$ . A stable algorithm nearly solves a nearby problem. $\mathrm { A n }$ algorithm for computing $f ( d )$ is backward stable if for some small $\Delta d ,$ the computed solution $\tilde { f } ( d )$ satisfies $\tilde { f } ( d ) = f ( d + \Delta d )$ . A backward stable algorithm exactly solves a nearby problem. Applied to a given linear system $A x = b $ , Gaussian elimination with complete pivoting is backward stable because the computed solution ˜x satisfies

$$
(A + \Delta) \tilde {x} = b
$$

and $\| \Delta \| / \| A \| \approx { \cal O } ( { \bf u } )$ . On the other hand, if b is specified by a matrix-vector product $b = M v$ , then

$$
(A + \Delta) \tilde {x} = M v + \delta
$$

where $\| \Delta \| / \| A \| \approx O ( \mathbf { u } )$ and $\delta / ( \lVert \mathbf { \nabla } M \rVert \lVert \mathbf { \nabla } v \rVert ) \approx { \cal O } ( \mathbf { u } )$ . Here, the underlying f is defined by $f \colon ( A , M , v ) \ \to \ A ^ { - 1 } ( M v )$ . In this case the algorithm is stable but not backward stable.

# Problems

P3.4.1 Let $A = L U$ be the LU factorization of n-by-n A with $| \ell _ { i j } | \le 1$ . Let $a _ { i } ^ { T }$ and $u _ { i } ^ { T }$ denote the ith rows of A and U , respectively. Verify the equation

$$
u _ {i} ^ {T} = a _ {i} ^ {T} - \sum_ {j = 1} ^ {i - 1} \ell_ {i j} u _ {j} ^ {T}
$$

and use it to show that $\parallel U \parallel _ { \infty } \leq 2 ^ { n - 1 } \parallel A \parallel _ { \infty }$ . (Hint: Take norms and use induction.)

P3.4.2 Show that if $P A Q = \mathrm { L U }$ is obtained via Gaussian elimination with complete pivoting, then no element of $U ( i , i ; n )$ is larger in absolute value than $| u _ { i i } |$ . Is this true with rook pivoting?

P3.4.3 Suppose $A \in \mathbb { R } ^ { n \times n }$ has an $_ { L U }$ factorization and that L and U are known. Give an algorithm which can compute the $( i , j )$ entry of $A ^ { - 1 }$ in approximately $( n - j ) ^ { 2 } + ( n - i ) ^ { 2 }$ flops.

P3.4.4 Suppose Xˆ is the computed inverse obtained via (3.4.12). Give an upper bound for $\Vert \ A \hat { X } - I \Vert _ { F } .$

P3.4.5 Extend Algorithm 3.4.3 so that it can produce the factorization (3.4.11). How many flops are required?

# Notes and References for 3.4

Papers concerned with element growth and pivoting include:

C.W. Cryer (1968). “Pivot Size in Gaussian Elimination,” Numer. Math. 12, 335–345.

J.K. Reid (1971). “A Note on the Stability of Gaussian Elimination,” J.Inst. Math. Applic. 8, 374–375.

P.A. Businger (1971). “Monitoring the Numerical Stability of Gaussian Elimination,” Numer. Math. 16, 360–361.

A.M. Cohen (1974). “A Note on Pivot Size in Gaussian Elimination,” Lin. Alg. Applic. 8, 361–68.

A.M. Erisman and J.K. Reid (1974). “Monitoring the Stability of the Triangular Factorization of a Sparse Matrix,” Numer. Math. 22, 183–186.

J. Day and B. Peterson (1988). “Growth in Gaussian Elimination,” Amer. Math. Monthly 95, 489–513.

N.J. Higham and D.J. Higham (1989). “Large Growth Factors in Gaussian Elimination with Pivoting,” SIAM J. Matrix Anal. Applic. 10, 155–164.

L.N. Trefethen and R.S. Schreiber (1990). “Average-Case Stability of Gaussian Elimination,” SIAM J. Matrix Anal. Applic. 11, 335–360.

N. Gould (1991). “On Growth in Gaussian Elimination with Complete Pivoting,” SIAM J. Matrix Anal. Applic. 12, 354–361.   
A. Edelman (1992). “The Complete Pivoting Conjecture for Gaussian Elimination is False,” Mathematica J. 2, 58–61.   
S.J. Wright (1993). “A Collection of Problems for Which Gaussian Elimination with Partial Pivoting is Unstable,” SIAM J. Sci. Stat. Comput. 14, 231–238.   
L.V. Foster (1994). “Gaussian Elimination with Partial Pivoting Can Fail in Practice,” SIAM J. Matrix Anal. Applic. 15, 1354–1362.   
A. Edelman and W. Mascarenhas (1995). “On the Complete Pivoting Conjecture for a Hadamard Matrix of Order 12,” Lin. Multilin. Alg. 38, 181–185.   
J.M. Pena (1996). “Pivoting Strategies Leading to Small Bounds of the Errors for Certain Linear Systems,” IMA J. Numer. Anal. 16, 141–153.   
J.L. Barlow and H. Zha (1998). “Growth in Gaussian Elimination, Orthogonal Matrices, and the 2-Norm,” SIAM J. Matrix Anal. Applic. 19, 807–815.   
P. Favati, M. Leoncini, and A. Martinez (2000). “On the Robustness of Gaussian Elimination with Partial Pivoting,” BIT 40, 62–73.   
As we mentioned, the size of L−1 is relevant to the growth factor. Thus, it is important to have an understanding of triangular matrix condition, see:   
D. Viswanath and L.N. Trefethen (1998). “Condition Numbers of Random Triangular Matrices,” SIAM J. Matrix Anal. Applic. 19, 564–581.   
The connection between small pivots and near singularity is reviewed in:   
T.F. Chan (1985). “On the Existence and Computation of LU Factorizations with Small Pivots,” Math. Comput. 42, 535–548.   
A pivot strategy that we did not discuss is pairwise pivoting. In this approach, 2-by-2 Gauss transformations are used to zero the lower triangular portion of A. The technique is appealing in certain multiprocessor environments because only adjacent rows are combined in each step, see:   
D. Sorensen (1985). “Analysis of Pairwise Pivoting in Gaussian Elimination,” IEEE Trans. Comput. C-34, 274–278.   
A related type of pivoting called tournament pivoting that is of interest in distributed memory computing is outlined in §3.6.3. For a discussion of rook pivoting and its properties, see:   
L.V. Foster (1997). “The Growth Factor and Efficiency of Gaussian Elimination with Rook Pivoting,” J. Comput. Appl. Math., 86, 177–194.   
G. Poole and L. Neal (2000). “The Rook’s Pivoting Strategy,” J. Comput. Appl. Math. 123, 353–369. X-W Chang (2002) “Some Features of Gaussian Elimination with Rook Pivoting,” BIT 42, 66–83.
