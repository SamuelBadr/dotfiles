# 4.3 Banded Systems

In many applications that involve linear systems, the matrix of coefficients is banded. This is the case whenever the equations can be ordered so that each unknown $x _ { i }$ appears in only a few equations in a “neighborhood”of the ith equation. Recall from §1.2.1 that $A = \left( a _ { i j } \right)$ has upper bandwidth q if $a _ { i j } = 0$ whenever $j > i + q$ and lower bandwidth p if $a _ { i j } = 0$ whenever $i > j + p$ . Substantial economies can be realized when solving banded systems because the triangular factors in LU, $\mathrm { G G } ^ { T }$ , and $\mathrm { L D L } ^ { T }$ are also banded.

# 4.3.1 Band LU Factorization

Our first result shows that if A is banded and $A = L U$ , then L inherits the lower bandwidth of A and U inherits the upper bandwidth of A.

Theorem 4.3.1. Suppose $A \in \mathbb { R } ^ { n \times n }$ has an LU factorization $A = L U$ . If A has upper bandwidth q and lower bandwidth p, then U has upper bandwidth q and L has lower bandwidth p.

Proof. The proof is by induction on n. Since

$$
A   =   \left[ \begin{array}{c c} \alpha & w ^ {T} \\ v & B \end{array} \right]   =   \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ 0 & B - v w ^ {T} / \alpha \end{array} \right] \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & I _ {n - 1} \end{array} \right]  .
$$

It is clear that $B - v w ^ { T } / \alpha$ has upper bandwidth q and lower bandwidth p because only the first q components of w and the first p components of v are nonzero. Let $L _ { 1 } U _ { 1 }$ be the LU factorization of this matrix. Using the induction hypothesis and the sparsity of w and v, it follows that

$$
L = \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & L _ {1} \end{array} \right], \qquad U = \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & U _ {1} \end{array} \right]
$$

have the desired bandwidth properties and satisfy $A = L U . \quad \ \bigtriangledown$

The specialization of Gaussian elimination to banded matrices having an LU factorization is straightforward.

Algorithm 4.3.1 (Band Gaussian Elimination) Given $A \in \mathbb { R } ^ { n \times n }$ with upper bandwidth q and lower bandwidth $p ,$ the following algorithm computes the factorization $A = L U$ , assuming it exists. $A ( i , j )$ is overwritten by $L ( i , j ) { \mathrm { ~ i f ~ } } i > j$ and by $U ( i , j )$ otherwise.

for k = 1:n - 1
    for i = k + 1:min{k + p, n}
    A(i, k) = A(i, k)/A(k, k)
    end
    for j = k + 1:min{k + q, n}
    for i = k + 1:min{k + p, n}
    A(i, j) = A(i, j) - A(i, k)·A(k, j)
    end
    end
end

If $n \gg p$ and $n \gg q ,$ , then this algorithm involves about $2 n p q$ flops. Effective implementations would involve band matrix data structures; see §1.2.5. A band version of Algorithm 4.1.1 (LDLT) is similar and we leave the details to the reader.

# 4.3.2 Band Triangular System Solving

Banded triangular system solving is also fast. Here are the banded analogues of Algorithms 3.1.3 and 3.1.4:

Algorithm 4.3.2 (Band Forward Substitution) Let $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ be a unit lower triangular matrix with lower bandwidth $p .$ Given $b \in \mathbb { R } ^ { n }$ , the following algorithm overwrites b with the solution to $L x = b$ .

for $j = 1:n$ for $i = j + 1:\min \{j + p,n\}$ $b(i) = b(i) - L(i,j)\cdot b(j)$ end   
end

If $n \gg p ,$ then this algorithm requires about 2np flops.

Algorithm 4.3.3 (Band Back Substitution) Let $U \in \mathbb { R } ^ { n \times n }$ be a nonsingular upper triangular matrix with upper bandwidth q. Given $b \in \mathbb { R } ^ { n }$ , the following algorithm overwrites b with the solution to $U x = b$ .

for $j = n: -1:1$ $b(j) = b(j)/U(j, j)$ for $i = \max\{1, j - q\}: j - 1$ $b(i) = b(i) - U(i, j) \cdot b(j)$ end
end

If $n \gg q$ , then this algorithm requires about 2nq flops.

# 4.3.3 Band Gaussian Elimination with Pivoting

Gaussian elimination with partial pivoting can also be specialized to exploit band structure in A. However, if $P A = L U$ , then the band properties of L and U are not quite so simple. For example, if A is tridiagonal and the first two rows are interchanged at the very first step of the algorithm, then $u _ { 1 3 }$ is nonzero. Consequently, row interchanges expand bandwidth. Precisely how the band enlarges is the subject of the following theorem.

Theorem 4.3.2. Suppose $A \in \mathbb { R } ^ { n \times n }$ is nonsingular and has upper and lower bandwidths q and $p ,$ respectively. If Gaussian elimination with partial pivoting is used to compute Gauss transformations

$$
M _ {j} = I - \alpha^ {(j)} e _ {j} ^ {T} \qquad j = 1: n - 1
$$

and permutations $P _ { 1 } , \ldots , P _ { n - 1 }$ such that $M _ { n - 1 } P _ { n - 1 } \cdot \cdot \cdot M _ { 1 } P _ { 1 } A = U$ is upper triangular, then U has upper bandwidth $p + q$ and $\alpha _ { i } ^ { ( j ) } = \theta$ whenever $i \leq j$ or $i > j + p$ .

Proof. Let $P A = L U$ be the factorization computed by Gaussian elimination with partial pivoting and recall that $P = P _ { n - 1 } \cdot \cdot \cdot P _ { 1 }$ . Write $P ^ { T } = \ [ e _ { s _ { 1 } } \ | \cdot \cdot \cdot | \ e _ { s _ { n } } \ ]$ , where $\{ s _ { 1 } , . . . , s _ { n } \}$ is a permutation of $\{ 1 , 2 , . . . , n \}$ . If $s _ { i } > i + p$ then it follows that the leading i-by-i principal submatrix of P A is singular, since $[ P A ] _ { i j } = a _ { s _ { i } , j }$ for $j = 1 \colon s _ { i } - p - 1$ and $s _ { i } - p - 1 \geq i$ . This implies that U and A are singular, a contradiction. Thus, $s _ { i } \leq i + p$ for $i = 1 { : } n$ and therefore, P A has upper bandwidth $p + q$ . It follows from Theorem 4.3.1 that U has upper bandwidth $p + q .$ . The assertion about the $\alpha ^ { ( j ) }$ can be verified by observing that $M _ { j }$ need only zero elements $( j + 1 , j ) , \dotsc , ( j + p , j )$ of the partially reduced matrix $P _ { j } M _ { j - 1 } P _ { j - 1 } \cdot \cdot \cdot _ { 1 } P _ { 1 } A$ .

Thus, pivoting destroys band structure in the sense that U becomes “wider” than A’s upper triangle, while nothing at all can be said about the bandwidth of L. However, since the jth column of L is a permutation of the jth Gauss vector $\alpha _ { j }$ , it follows that L has at most $p + 1$ nonzero elements per column.

# 4.3.4 Hessenberg LU

As an example of an unsymmetric band matrix computation, we show how Gaussian elimination with partial pivoting can be applied to factor an upper Hessenberg matrix H. (Recall that if H is upper Hessenberg then $h _ { i j } = 0 , i > j + 1 . )$ After $k - 1$ steps of Gaussian elimination with partial pivoting we are left with an upper Hessenberg matrix of the form

$$
\left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right], \qquad k = 3, n = 5.
$$

By virtue of the special structure of this matrix, we see that the next permutation, $P _ { 3 } .$ , is either the identity or the identity with rows 3 and 4 interchanged. Moreover, the next Gauss transformation $M _ { k }$ has a single nonzero multiplier in the $( k + 1 , k )$ position. This illustrates the kth step of the following algorithm.

Algorithm 4.3.4 (Hessenberg LU) Given an upper Hessenberg matrix $H \in \mathbb { R } ^ { n \times n }$ , the following algorithm computes the upper triangular matrix $M _ { n - 1 } P _ { n - 1 } \cdot \cdot \cdot M _ { 1 } P _ { 1 } H = U$ where each $P _ { k }$ is a permutation and each $M _ { k }$ is a Gauss transformation whose entries are bounded by unity. $H ( i , k )$ is overwritten with $U ( i , k )$ if $i \leq k$ and by $- [ M _ { k } ] _ { k + 1 , k }$ if $i = k + 1$ . An integer vector $p i v ( 1 { : } n - 1 )$ encodes the permutations. If $P _ { k } = I$ , then $p i v ( k ) = 0$ . If $P _ { k }$ interchanges rows k and $k + 1$ , then $p i v ( k ) = 1$ .

for k = 1:n - 1
    if $|H(k,k)| < |H(k+1,k)|$ $piv(k) = 1; H(k,k:n) \leftrightarrow H(k+1,k:n)$ else $piv(k) = 0$ end
    if $H(k,k) \neq 0$ $\tau = H(k+1,k)/H(k,k)$ $H(k+1,k+1:n) = H(k+1,k+1:n) - \tau \cdot H(k,k+1:n)$ $H(k+1,k) = \tau$ end
end

This algorithm requires $n ^ { 2 }$ flops.

# 4.3.5 Band Cholesky

The rest of this section is devoted to banded $A x = b$ problems where the matrix A is also symmetric positive definite. The fact that pivoting is unnecessary for such matrices leads to some very compact, elegant algorithms. In particular, it follows from Theorem 4.3.1 that if $A = G G ^ { T }$ is the Cholesky factorization of $A ,$ , then G has the same lower bandwidth as A. This leads to the following banded version of Algorithm 4.2.1.

Algorithm 4.3.5 (Band Cholesky) Given a symmetric positive definite $A \in \mathbb { R } ^ { n \times n }$ with bandwidth $p ,$ the following algorithm computes a lower triangular matrix $G$ with lower bandwidth $p$ such that $A = G G ^ { T }$ . For all $i \geq j , G ( i , j )$ overwrites $A ( i , j )$ .

$$
\begin{array}{l} \text { for } j = 1: n \\ \text { for } k = \max (1, j - p): j - 1 \\ \lambda = \min (k + p, n) \\ A (j: \lambda , j) = A (j: \lambda , j) - A (j, k) \cdot A (j: \lambda , k) \\ \lambda = \min (j + p, n) \\ A (j: \lambda , j) = A (j: \lambda , j) / \sqrt {A (j , j)} \\ \end{array}
$$

end

If $n \gg p .$ , then this algorithm requires about $n ( p ^ { 2 } + 3 p )$ flops and n square roots. Of course, in a serious implementation an appropriate data structure for A should be used. For example, if we just store the nonzero lower triangular part, then a $( p + 1 ) – \mathrm { b y } – n$ array would suffice.

If our band Cholesky procedure is coupled with appropriate band triangular solve routines, then approximately $n p ^ { 2 } + 7 n p + 2 n$ flops and n square roots are required to solve $A x = b$ . For small $p$ it follows that the square roots represent a significant portion of the computation and it is preferable to use the LDLT approach. Indeed, a careful flop count of the steps $A = L D L ^ { T } , L y = b , D z = y$ , and $L ^ { T } x \ = \ z$ reveals that $n p ^ { 2 } + 8 n p + n$ flops and no square roots are needed.

# 4.3.6 Tridiagonal System Solving

As a sample narrow band $\mathrm { L D L } ^ { T }$ solution procedure, we look at the case of symmetric positive definite tridiagonal systems. Setting

$$
L = \left[ \begin{array}{c c c c} 1 & 0 & \dots & 0 \\ \ell_ {1} & 1 & & \vdots \\ \vdots & \ddots & \ddots & 0 \\ 0 & \dots & \ell_ {n - 1} & 1 \end{array} \right]
$$

and $D = \operatorname { d i a g } ( d _ { 1 } , \dotsc , d _ { n } )$ , we deduce from the equation $A = L D L ^ { T }$ that

$$
\begin{array}{l} a _ {1 1} \quad = d _ {1}, \\ a _ {k, k - 1} = \ell_ {k - 1} d _ {k - 1}, \quad k = 2: n, \\ a _ {k k} \quad = d _ {k} + \ell_ {k - 1} ^ {2} d _ {k - 1} = d _ {k} + \ell_ {k - 1} a _ {k, k - 1}, \quad k = 2: n. \\ \end{array}
$$

Thus, the $d _ { i }$ and $\ell _ { i }$ can be resolved as follows:

$$
\begin{array}{l} d _ {1} = a _ {1 1} \\ \ell_ {k - 1} = a _ {k, k - 1} / d _ {k - 1} \\ d _ {k} = a _ {k k} - \ell_ {k - 1} a _ {k, k - 1} \\ \end{array}
$$

end

To obtain the solution to Ax = b we solve $L y = b , \ D z \ = \ y$ , and $L ^ { T } x = z$ . With overwriting we obtain

Algorithm 4.3.6 (Symmetric, Tridiagonal, Positive Definite System Solver) Given an n-by-n symmetric, tridiagonal, positive definite matrix A and $b \in \mathbb { R } ^ { n }$ , the following algorithm overwrites b with the solution to $A x = b $ . It is assumed that the diagonal of A is stored in α(1:n) and the superdiagonal in $\beta ( 1 : n - 1 )$ .

for k = 2:n

$$
t = \beta (k - 1), \beta (k - 1) = t / \alpha (k - 1), \alpha (k) = \alpha (k) - t \cdot \beta (k - 1)
$$

end

for k = 2:n

$$
b (k) = b (k) - \beta (k - 1) \cdot \beta (k - 1)
$$

end

$$
b (n) = b (n) / \alpha (n)
$$

for k = n − 1: − 1:1

$$
b (k) = b (k) / \alpha (k) - \beta (k) \cdot b (k + 1)
$$

end

This algorithm requires 8n flops.

# 4.3.7 Vectorization Issues

The tridiagonal example brings up a sore point: narrow band problems and vectorization do not mix. However, it is sometimes the case that large, independent sets of such problems must be solved at the same time. Let us examine how such a computation could be arranged in light of the issues raised in §1.5. For simplicity, assume that we must solve the n-by-n unit lower bidiagonal systems

$$
A ^ {(k)} x ^ {(k)} = b ^ {(k)}, \qquad k = 1: m,
$$

and that $m \gg n$ . Suppose we have arrays $E ( 1 { : } n - 1 , 1 { : } m )$ and $B ( 1 { : } n , 1 { : } m )$ with the property that $E ( 1 { : } n - 1 , k )$ houses the subdiagonal of $A ^ { ( k ) }$ and $B ( 1 { : } n , k )$ houses the kth right-hand side $b ^ { ( k ) }$ . We can overwrite $b ^ { ( k ) }$ with the solution $x ^ { ( k ) }$ as follows:

for k = 1:m

for i = 2:n

$$
B (i, k) = B (i, k) - E (i - 1, k) \cdot B (i - 1, k)
$$

end

end

This algorithm sequentially solves each bidiagonal system in turn. Note that the inner loop does not vectorize because of the dependence of $B ( i , k )$ on $B ( i - 1 , k )$ . However, if we interchange the order of the two loops, then the calculation does vectorize:

for $i = 2 { : } n$

$$
B (i,:) = B (i,:) - E (i - 1,:). * B (i - 1,:)
$$

end

A column-oriented version can be obtained simply by storing the matrix subdiagonals by row in E and the right-hand sides by row in B:

for $i = 2 { : } n$

$$
B (:, i) = B (:, i) - E (:, i - 1). * B (:, i - 1)
$$

end

Upon completion, the transpose of solution $x ^ { ( k ) }$ is housed on $B ( k , : )$ .

# 4.3.8 The Inverse of a Band Matrix

In general, the inverse of a nonsingular band matrix A is full. However, the off-diagonal blocks of $A ^ { - 1 }$ have low rank.

Theorem 4.3.3. Suppose

$$
A = \left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right]
$$

is nonsingular and has lower bandwidth p and upper bandwidth q. Assume that the diagonal blocks are square. If

$$
A ^ {- 1} = X = \left[ \begin{array}{l l} X _ {1 1} & X _ {1 2} \\ X _ {2 1} & X _ {2 2} \end{array} \right]
$$

is partitioned conformably, then

$$
\operatorname{rank} (X _ {2 1}) \leq p, \tag {4.3.1}
$$

$$
\operatorname{rank} (X _ {1 2}) \leq q. \tag {4.3.2}
$$

Proof. Assume that $A _ { 1 1 }$ and $A _ { 2 2 }$ are nonsingular. From the equation $A X = I$ we conclude that

$$
A _ {2 1} X _ {1 1} + A _ {2 2} X _ {2 1} = 0,
$$

$$
A _ {1 1} X _ {1 2} + A _ {1 2} X _ {2 2} = 0,
$$

and so

$$
\operatorname{rank} \left(X _ {2 1}\right) = \operatorname{rank} \left(A _ {2 2} ^ {- 1} A _ {2 1} X _ {1 1}\right) \leq \operatorname{rank} \left(A _ {2 1}\right)
$$

$$
\operatorname{rank} \left(X _ {1 2}\right) = \operatorname{rank} \left(A _ {1 1} ^ {- 1} A _ {1 2} X _ {2 2}\right) \leq \operatorname{rank} \left(A _ {1 2}\right).
$$

From the bandedness assumptions it follows that $A _ { 2 1 }$ has at most p nonzero rows and $A _ { 1 2 }$ has at most q nonzero rows. Thus, rank $( A _ { 2 1 } ) \le p$ and rank $( A _ { 1 2 } ) \le q$ which proves the theorem for the case when both $A _ { 1 1 }$ and $A _ { 2 2 }$ are nonsingular. A simple limit argument can be used to handle the situation when $A _ { 1 1 }$ and/or $A _ { 2 2 }$ are singular. See P4.3.11.

It can actually be shown that rank $( A _ { 2 1 } ) = \mathsf { r a n k } ( X _ { 2 1 } )$ and rank $( A _ { 1 2 } ) = \mathsf { r a n k } ( X _ { 1 2 } )$ . See Strang and Nguyen (2004). As we will see in §11.5.9 and §12.2, the low-rank, offdiagonal structure identified by the theorem has important algorithmic ramifications.

# 4.3.9 Band Matrices with Banded Inverse

If $A \in \mathbb { R } ^ { n \times n }$ is a product

$$
A = F _ {1} \dots F _ {N} \tag {4.3.3}
$$

and each $F _ { i } \in \mathbb { R } ^ { n \times n }$ is block diagonal with 1-by-1 and 2-by-2 diagonal blocks, then it follows that both A and

$$
A ^ {- 1} = F _ {N} ^ {- 1} \dots F _ {1} ^ {- 1}
$$

are banded, assuming that N is not too big. For example, if

$$
A = \left[ \begin{array}{c c c c c c c c c} \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & \times & \times & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & \times & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & \times & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & \times & \times & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times \end{array} \right] \left[ \begin{array}{c c c c c c c c c} \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & \times & \times & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & \times & \times & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & \times & \times & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & \times & \times & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times \end{array} \right]
$$

then

$$
A = \left[ \begin{array}{c c c c c c c c c} \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & \times & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & \times & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & \times & \times & \times & \times & 0 & 0 & 0 \\ 0 & 0 & \times & \times & \times & \times & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & \times & \times & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times & \times & \times & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times \end{array} \right], \quad A ^ {- 1} = \left[ \begin{array}{c c c c c c c c c} \times & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & \times & \times & \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & \times & \times & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & \times & \times & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & \times & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times \end{array} \right].
$$

Strang (2010a, 2010b) has pointed out a very important “reverse” fact. If A and $A ^ { - 1 }$ are banded, then there is a factorization of the form (4.3.3) with relatively small N . Indeed, he shows that N is very small for certain types of matrices that arise in signal processing. An important consequence of this is that both the forward transform Ax and the inverse transform $A ^ { - 1 } x$ can be computed very fast.

# Problems

P4.3.1 Develop a version of Algorithm 4.3.1 which assumes that the matrix A is stored in band format style. (See §1.2.5.)

P4.3.2 Show how the output of Algorithm 4.3.4 can be used to solve the upper Hessenberg system $H x = b $ .

P4.3.3 Show how Algorithm 4.3.4 could be used to solve a lower hessenberg system $H x = b$ .

P4.3.4 Give an algorithm for solving an unsymmetric tridiagonal system $A x = b$ that uses Gaussian elimination with partial pivoting. It should require only four n-vectors of floating point storage for the factorization.

P4.3.5 (a) For $C \in \mathbb { R } ^ { n \times n }$ define the profile indices $m ( C , i ) = \operatorname* { m i n } \{ j { : } c _ { i j } \neq 0 \}$ , where $i = 1 { : } n$ . Show that if $\mathbf { \Psi } _ { A } \stackrel { \mathbf { \tilde { = } } } { = } \mathbf { \Psi } _ { G } \mathbf { { \cal { G } } } ^ { T }$ is the Cholesky factorization of A, then $m ( A , i ) = m ( G , i )$ for $i = 1 { : } n$ (We say that G has the same profile as A.) (b) Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite with profile indices $m _ { i } = m ( A , i )$ where $i = 1 { : } n$ . Assume that A is stored in a one-dimensional array v as follows:

$$
v = (a _ {1 1}, a _ {2, m _ {2}}, \ldots , a _ {2 2}, a _ {3, m _ {3}}, \ldots , a _ {3 3}, \ldots , a _ {n, m _ {n}}, \ldots , a _ {n n}).
$$

Give an algorithm that overwrites v with the corresponding entries of the Cholesky factor G and then uses this factorization to solve $A x = b ,$ . How many flops are required? (c) For $C \in \mathbb { R } ^ { n \times n }$ define $p ( C , i )$ $= \operatorname* { m a x } \{ j { : } c _ { i j } \neq 0 \}$ . Suppose that $A \in \mathbb { R } ^ { n \times n }$ has an LU factorization $A = L U$ and that

$$
m (A, 1) \leq m (A, 2) \leq \dots \leq m (A, n),
$$

$$
p (A, 1) \leq p (A, 2) \leq \dots \leq p (A, n).
$$

Show that $m ( A , i ) = m ( L , i )$ and $p ( A , i ) = p ( U , i )$ for $i = 1 { : } n$ .

P4.3.6 Develop a gaxpy version of Algorithm 4.3.1.

P4.3.7 Develop a unit stride, vectorizable algorithm for solving the symmetric positive definite tridiagonal systems $A ^ { ( k ) } x ^ { ( k ) } = b ^ { ( \dot { k } ) }$ . Assume that the diagonals, superdiagonals, and right hand sides are stored by row in arrays $D , E ,$ and B and that $b ^ { ( k ) }$ is overwritten with $x ^ { ( k ) }$ .

P4.3.8 Give an example of a 3-by-3 symmetric positive definite matrix whose tridiagonal part is not positive definite.

P4.3.9 Suppose a symmetric positive definite matrix $A \in \mathbb { R } ^ { n \times n }$ has the “arrow structure”, e.g.,

$$
A = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & \times & 0 & 0 & 0 \\ \times & 0 & \times & 0 & 0 \\ \times & 0 & 0 & \times & 0 \\ \times & 0 & 0 & 0 & \times \end{array} \right].
$$

(a) Show how the linear system $A x = b$ can be solved with $O ( n )$ flops using the Sherman-Morrison-Woodbury formula. (b) Determine a permutation matrix P so that the Cholesky factorization

$$
P A P ^ {T} = G G ^ {T}
$$

can be computed with O(n) flops.

P4.3.10 Suppose $A \in \mathbb { R } ^ { n \times n }$ is tridiagonal, positive definite, but not symmetric. Give an efficient algorithm for computing the largest entry of $\vert S T ^ { - 1 } S \vert$ where $S = ( A - A ^ { T } ) / 2$ and $T = ( A + A ^ { T } ) / 2$ .

P4.3.11 Show that if $A \in \mathbb { R } ^ { n \times n }$ and $\epsilon > 0$ , then there is a $B \in \mathbb { R } ^ { n \times n }$ such that $\| A - B \| \leq \epsilon$ and B has the property that all its principal submatrices are nonsingular. Use this result to formally complete the proof of Theorem 4.3.3.

P4.3.12 Give an upper bound on the bandwidth of the matrix A in (4.3.3).

P4.3.13 Show that $A ^ { T }$ and $A ^ { - 1 }$ have the same upper and lower bandwidths in (4.3.3).

P4.3.14 For the $A = F _ { 1 } F _ { 2 }$ example in 4.3.9, show that A(2:3, :), A(4:5, :), A(6:7, :), . . . each consist of two singular 2-by-2 blocks.

# Notes and References for §4.3

Representative papers on the topic of banded systems include:

R.S. Martin and J.H. Wilkinson (1965). “Symmetric Decomposition of Positive Definite Band Matrices,” Numer. Math. 7, 355–61.   
R. S. Martin and J.H. Wilkinson (1967). “Solution of Symmetric and Unsymmetric Band Equations and the Calculation of Eigenvalues of Band Matrices,” Numer. Math. 9, 279–301.   
E.L. Allgower (1973). “Exact Inverses of Certain Band Matrices,” Numer. Math. 21, 279–284.   
Z. Bohte (1975). “Bounds for Rounding Errors in the Gaussian Elimination for Band Systems,” J. Inst. Math. Applic. 16, 133–142.   
L. Kaufman (2007). “The Retraction Algorithm for Factoring Banded Symmetric Matrices,” Numer. Lin. Alg. Applic. 14, 237–254.   
C. Vomel and J. Slemons (2009). “Twisted Factorization of a Banded Matrix,” BIT 49, 433–447.

Tridiagonal systems are particularly important, see:

C. Fischer and R.A. Usmani (1969). “Properties of Some Tridiagonal Matrices and Their Application to Boundary Value Problems,” SIAM J. Numer. Anal. 6, 127–142.   
D.J. Rose (1969). “An Algorithm for Solving a Special Class of Tridiagonal Systems of Linear Equations,” Commun. ACM 12, 234–236.   
M.A. Malcolm and J. Palmer (1974). “A Fast Method for Solving a Class of Tridiagonal Systems of Linear Equations,” Commun. ACM 17, 14–17.   
N.J. Higham (1986). “Efficient Algorithms for Computing the Condition Number of a Tridiagonal Matrix,” SIAM J. Sci. Stat. Comput. 7, 150–165.   
N.J. Higham (1990). “Bounding the Error in Gaussian Elimination for Tridiagonal Systems,” SIAM J. Matrix Anal. Applic. 11, 521–530.   
I.S. Dhillon (1998). “Reliable Computation of the Condition Number of a Tridiagonal Matrix in O(n) Time,” SIAM J. Matrix Anal. Applic. 19, 776–796.   
I. Bar-On and M. Leoncini (2000). “Reliable Solution of Tridiagonal Systems of Linear Equations,” SIAM J. Numer. Anal. 38, 1134–1153.   
M.I. Bueno and F.M. Dopico (2004). “Stability and Sensitivity of Tridiagonal LU Factorization without Pivoting,” BIT 44, 651–673.   
J.R. Bunch and R.F. Marcia (2006). “A Simplified Pivoting Strategy for Symmetric Tridiagonal Matrices,” Numer. Lin. Alg. 13, 865–867.

For a discussion of parallel methods for banded problems, see:

H.S. Stone (1975). “Parallel Tridiagonal Equation Solvers,” ACM Trans. Math. Softw. 1, 289–307.   
I. Bar-On, B. Codenotti and M. Leoncini (1997). “A Fast Parallel Cholesky Decomposition Algorithm for Tridiagonal Symmetric Matrices,” SIAM J. Matrix Anal. Applic. 18, 403–418.   
G.H. Golub, A.H. Sameh, and V. Sarin (2001). “A Parallel Balance Scheme for Banded Linear Systems,” Num. Lin. Alg. 8, 297–316.   
S. Rao and Sarita (2008). “Parallel Solution of Large Symmetric Tridiagonal Linear Systems,” Parallel Comput. 34, 177–197.

Papers that are concerned with the structure of the inverse of a band matrix include:

E. Asplund (1959). “Inverses of Matrices {aij } Which Satisfy aij = 0 for j > i + p,” Math. Scand. 7, 57–60.   
C.A. Micchelli (1992). “Banded Matrices with Banded Inverses,” J. Comput. Appl. Math. 41, 281-300.   
G. Strang and T. Nguyen (2004). “The Interplay of Ranks of Submatrices,” SIAM Review 46, 637–648.   
G. Strang (2010a). “Fast Transforms: Banded Matrices with Banded Inverses,” Proc. National Acad. Sciences 107, 12413-12416.   
G. Strang (2010b). “Banded Matrices with Banded Inverses and A = LP U,” Proceedings International Congress of Chinese Mathematicians, Beijing.

A pivotal result in this arena is the nullity theorem, a more general version of Theorem 4.3.3, see:   
R. Vandebril, M. Van Barel, and N. Mastronardi (2008). Matrix Computations and Semiseparable Matrices, Volume I Linear Systems, Johns Hopkins University Press, Baltimore, MD., 37–40.
