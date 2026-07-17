# 5.2 The QR Factorization

A rectangular matrix $A \in \mathbb { R } ^ { m \times n }$ can be factored into a product of an orthogonal matrix $Q \in \mathbb { R } ^ { m \times m }$ and an upper triangular matrix $R \in \mathbb { R } ^ { m \times n }$ :

$$
A = Q R.
$$

This factorization is referred to as the QR factorization and it has a central role to play in the linear least squares problem. In this section we give methods for computing QR based on Householder, block Householder, and Givens transformations. The QR factorization is related to the well-known Gram-Schmidt process.

# 5.2.1 Existence and Properties

We start with a constructive proof of the QR factorization.

Theorem 5.2.1 (QR Factorization). If $A \in \mathbb { R } ^ { m \times n }$ , then there exists an orthogonal $Q \in \mathbb { R } ^ { m \times m }$ and an upper triangular $R \in \mathbb { R } ^ { m \times n }$ so that $A = Q R$ .

Proof. We use induction. Suppose n = 1 and that Q is a Householder matrix so that if $R = Q ^ { T } A$ , then $R ( 2 { : } m ) = 0$ . It follows that $A = Q R$ is a QR factorization of A. For general n we partition A,

$$
A = \left[ \begin{array}{c c} A _ {1} & v \end{array} \right],
$$

where $v = A ( : , n )$ . By induction, there exists an orthogonal $Q _ { 1 } \in \mathbb { R } ^ { m \times m }$ so that $R _ { 1 } = Q _ { 1 } ^ { T } A _ { 1 }$ is upper triangular. Set $w = Q ^ { T } v$ and let $w ( n { : } m ) = Q _ { 2 } R _ { 2 }$ be the QR factorization of $w ( n { : } m )$ . If

$$
Q = Q _ {1} \left[ \begin{array}{c c} I _ {n - 1} & 0 \\ 0 & Q _ {2} \end{array} \right],
$$

then

$$
A = Q \left[ \begin{array}{c c c} R _ {1} & w (1: n - 1) \\ & R _ {2} \end{array} \right]
$$

is a QR factorization of A.

The columns of Q have an important connection to the range of A and its orthogonal complement.

Theorem 5.2.2. If A = QR is a QR factorization of a full column rank $A \in \mathbb { R } ^ { m \times n }$ and

$$
A = \left[ a _ {1} \mid \dots \mid a _ {n} \right],
$$

$$
Q = \left[ q _ {1} \mid \dots \mid q _ {m} \right]
$$

are column partitionings, then for k = 1:n

$$
\operatorname{span} \{a _ {1}, \dots , a _ {k} \} = \operatorname{span} \{q _ {1}, \dots , q _ {k} \} \tag {5.2.1}
$$

and $r _ { k k } \neq 0$ . Moreover, $i f Q _ { 1 } = Q ( 1 { : } m , 1 { : } n ) , Q _ { 2 } = Q ( 1 { : } m , n + 1 { : } m )$ , and $R _ { 1 } ~ =$ $R ( 1 { : } n , 1 { : } n )$ , then

$$
\operatorname{ran} (A) \quad = \operatorname{ran} (Q _ {1}),
$$

$$
\operatorname{ran} (A) ^ {\perp} = \operatorname{ran} (Q _ {2}),
$$

and

$$
A = Q _ {1} R _ {1}. \tag {5.2.2}
$$

Proof. Comparing the kth columns in A = QR we conclude that

$$
a _ {k} = \sum_ {i = 1} ^ {k} r _ {i k} q _ {i} \in \operatorname{span} \{q _ {1}, \dots , q _ {k} \}, \tag {5.2.3}
$$

and so

$$
\operatorname{span} \left\{a _ {1}, \dots , a _ {k} \right\} \subseteq \operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\}.
$$

If $r _ { k k } = 0$ , then $a _ { 1 } , \ldots , a _ { k }$ are dependent. Thus, R cannot have a zero on its diagonal and so span $\{ a _ { 1 } , \ldots , a _ { k } \}$ has dimension k. Coupled with (5.2.3) this establishes (5.2.1). To prove (5.2.2) we note that

$$
A = Q R = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \end{array} \right] \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right] = Q _ {1} R _ {1}. \quad \square
$$

The matrices $Q _ { 1 } = Q ( 1 { : } m , 1 { : } n )$ and $Q _ { 2 } = Q ( 1 { : } m , n + 1 { : } m )$ can be easily computed from a factored form representation of Q. We refer to (5.2.2) as the thin QR factorization. The next result addresses its uniqueness.

Theorem 5.2.3 (Thin QR Factorization). Suppose $A \in \mathbb { R } ^ { m \times n }$ has full column rank. The thin QR factorization

$$
A = Q _ {1} R _ {1}
$$

is unique where $Q _ { 1 } \in \mathbb { R } ^ { m \times n }$ has orthonormal columns and $R _ { 1 }$ is upper triangular with positive diagonal entries. Moreover, $R _ { 1 } = G ^ { T }$ where G is the lower triangular Cholesky factor of $A ^ { \bar { T } } A$ .

Proof. Since $A ^ { T } A = ( Q _ { 1 } R _ { 1 } ) ^ { T } ( Q _ { 1 } R _ { 1 } ) = R _ { 1 } ^ { T } R _ { 1 }$ we see that $G = R _ { 1 } ^ { T }$ is the Cholesky factor of $A ^ { T } A$ . This factor is unique by Theorem 4.2.7. Since $Q _ { 1 } \stackrel { - } { = } A R _ { 1 } ^ { - 1 }$ it follows that $Q _ { 1 }$ is also unique.

How are $Q _ { 1 }$ and $R _ { 1 }$ affected by perturbations in A? To answer this question we need to extend the notion of 2-norm condition to rectangular matrices. Recall from §2.6.2 that for square matrices, $\kappa _ { 2 } ( A )$ is the ratio of the largest to the smallest singular value. For rectangular matrices A with full column rank we continue with this definition:

$$
\kappa_ {2} (A) = \frac {\sigma_ {\max} (A)}{\sigma_ {\min} (A)}. \tag {5.2.4}
$$

If the columns of A are nearly dependent, then this quotient is large. Stewart (1993) has shown that $O ( \epsilon )$ relative error in A induces $O ( \epsilon { \cdot } \kappa _ { 2 } ( A ) )$ error in $Q _ { 1 }$ and $R _ { 1 }$ .

# 5.2.2 Householder QR

We begin with a QR factorization method that utilizes Householder transformations. The essence of the algorithm can be conveyed by a small example. Suppose $m = 6$ , $n = 5$ , and assume that Householder matrices $H _ { 1 }$ and $H _ { 2 }$ have been computed so that

$$
H _ {2} H _ {1} A = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \mathbf {x} & \times & \times \\ 0 & 0 & \mathbf {x} & \times & \times \\ 0 & 0 & \mathbf {x} & \times & \times \\ 0 & 0 & \mathbf {x} & \times & \times \end{array} \right].
$$

Concentrating on the highlighted entries, we determine a Householder matrix $\tilde { H } _ { 3 } \in \mathbb { R } ^ { 4 \times 4 }$ such that

$$
\tilde {H} _ {3} \left[ \begin{array}{l} \mathbf {x} \\ \mathbf {x} \\ \mathbf {x} \\ \mathbf {x} \end{array} \right] = \left[ \begin{array}{l} \times \\ 0 \\ 0 \\ 0 \end{array} \right].
$$

If $H _ { 3 } = \mathrm { d i a g } ( I _ { 2 } , \tilde { H } _ { 3 } )$ , then

$$
H _ {3} H _ {2} H _ {1} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right].
$$

After n such steps we obtain an upper triangular $H _ { n } H _ { n - 1 } \cdot \cdot \cdot H _ { 1 } A = R$ and so by setting $Q = H _ { 1 } \cdot \cdot \cdot H _ { n }$ we obtain $A = Q R$ .

Algorithm 5.2.1 (Householder QR) Given $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n ,$ , the following algorithm finds Householder matrices $H _ { 1 } , \ldots , H _ { n }$ such that if $Q \ = \ H _ { 1 } \cdot \cdot \cdot H _ { n } $ , then $Q ^ { \bar { T } } A = R$ is upper triangular. The upper triangular part of A is overwritten by the upper triangular part of R and components $j + 1 { : } m$ of the jth Householder vector are stored in $A ( j + 1 ; m , j ) , j < m$ .

for $j = 1 { : } n$

$$
[ v, \beta ] = \operatorname{house} (A (j: m, j))
$$

$$
A (j: m, j: n) = (I - \beta v v ^ {T}) A (j: m, j: n)
$$

$$
\text { if } j <   m
$$

$$
A (j + 1: m, j) = v (2: m - j + 1)
$$

$$
\mathrm{end}
$$

end

This algorithm requires $2 n ^ { 2 } ( m - n / 3 )$ flops.

To clarify how A is overwritten, if

$$
v ^ {(j)} = [ \underbrace {0 , \ldots , 0} _ {j - 1}, 1, v _ {j + 1} ^ {(j)}, \ldots , v _ {m} ^ {(j)} ] ^ {T}
$$

is the jth Householder vector, then upon completion

$$
A = \left[ \begin{array}{l l l l l} r _ {1 1} & r _ {1 2} & r _ {1 3} & r _ {1 4} & r _ {1 5} \\ v _ {2} ^ {(1)} & r _ {2 2} & r _ {2 3} & r _ {2 4} & r _ {2 5} \\ v _ {3} ^ {(1)} & v _ {3} ^ {(2)} & r _ {3 3} & r _ {3 4} & r _ {3 5} \\ v _ {4} ^ {(1)} & v _ {4} ^ {(2)} & v _ {4} ^ {(3)} & r _ {4 4} & r _ {4 5} \\ v _ {5} ^ {(1)} & v _ {5} ^ {(2)} & v _ {5} ^ {(3)} & v _ {5} ^ {(4)} & r _ {5 5} \\ v _ {6} ^ {(1)} & v _ {6} ^ {(2)} & v _ {6} ^ {(3)} & v _ {6} ^ {(4)} & v _ {6} ^ {(5)} \end{array} \right].
$$

If the matrix $Q = H _ { 1 } \cdot \cdot \cdot H _ { n }$ is required, then it can be accumulated using (5.1.5). This accumulation requires $4 ( m ^ { 2 } n - m n ^ { 2 } + n ^ { 3 } / 3 )$ flops. Note that the β-values that arise in Algorithm 5.2.1 can be retrieved from the stored Householder vectors:

$$
\beta_ {j} = \frac {2}{1 + \parallel A (j + 1 : m , j) \parallel^ {2}}.
$$

We mention that the computed upper triangular matrix $\hat { R }$ is the exact $R$ for a nearby A in the sense that $Z ^ { T } ( \bar { A } + E ) = \hat { R }$ where $Z$ is some exact orthogonal matrix and $\parallel E \parallel _ { 2 } \approx \mathbf { u } \parallel A \parallel _ { 2 }$ .

# 5.2.3 Block Householder QR Factorization

Algorithm 5.2.1 is rich in the level-2 operations of matrix-vector multiplication and outer product updates. By reorganizing the computation and using the WY representation discussed in §5.1.7 we can obtain a level-3 procedure. The idea is to apply the underlying Householder transformations in clusters of size r. Suppose $n = 1 2$ and $r = 3$ . The first step is to generate Householders $H _ { 1 } , H _ { 2 }$ , and $H _ { 3 }$ as in Algorithm 5.2.1. However, unlike Algorithm 5.2.1 where each $H _ { i }$ is applied across the entire remaining submatrix, we apply only $H _ { 1 } , H _ { 2 }$ , and $H _ { 3 }$ to $A ( : , 1 { : } 3 )$ . After this is accomplished we generate the block representation $H _ { 1 } H _ { 2 } H _ { 3 } = I - W _ { 1 } Y _ { 1 } ^ { T }$ and then perform the level-3 update

$$
A (:, 4: 1 2) = (I - W Y ^ {T}) A (:, 4: 1 2).
$$

Next, we generate $H _ { 4 } , H _ { 5 }$ , and $H _ { 6 }$ as in Algorithm 5.2.1. However, these transformations are not applied to $A ( : , 7 : 1 2 )$ until their block representation $H _ { 4 } H _ { 5 } H _ { 6 } = I - W _ { 2 } Y _ { 2 } ^ { T }$ is found. This illustrates the general pattern.

Algorithm 5.2.2 (Block Householder QR) If $A \in \mathbb { R } ^ { m \times n }$ and r is a positive integer, then the following algorithm computes an orthogonal $Q \in \mathbb { R } ^ { m \times m }$ and an upper triangular $R \in \mathbb { R } ^ { m \times n }$ so that $A = Q R$ .

$$
Q = I _ {m}; \lambda = 1; k = 0
$$

while $\lambda \leq n$

$$
\tau \leftarrow \min (\lambda + r - 1, n); k = k + 1
$$

Use Algorithm 5.2.1, to upper triangularize $A ( \lambda { : } m , \lambda { : } \tau )$ , generating Householder matrices $H _ { \lambda } , \ldots , H _ { \tau }$ .

Use Algorithm 5.1.2 to get the block representation

$$
I - W _ {k} Y _ {k} = H _ {\lambda} \dots H _ {\tau}.
$$

$$
A (\lambda : m, \tau + 1: n) = (I - W _ {k} Y _ {k} ^ {T}) ^ {T} A (\lambda : m, \tau + 1: n)
$$

$$
Q (:, \lambda : m) = Q (:, \lambda : m) (I - W _ {k} Y _ {k} ^ {T})
$$

$$
\lambda = \tau + 1
$$

end

The zero-nonzero structure of the Householder vectors that define $H _ { \lambda } , \ldots , H _ { \tau }$ implies that the first $\lambda - 1$ rows of $W _ { k }$ and $Y _ { k }$ are zero. This fact would be exploited in a practical implementation.

The proper way to regard Algorithm 5.2.2 is through the partitioning

$$
A = \left[ A _ {1} \mid \dots \mid A _ {N} \right], \qquad N = \operatorname{ceil} (n / r)
$$

where block column $A _ { k }$ is processed during the kth step. In the kth step of the reduction, a block Householder is formed that zeros the subdiagonal portion of $A _ { k }$ . The remaining block columns are then updated.

The roundoff properties of Algorithm 5.2.2 are essentially the same as those for Algorithm 5.2.1. There is a slight increase in the number of flops required because of the W -matrix computations. However, as a result of the blocking, all but a small fraction of the flops occur in the context of matrix multiplication. In particular, the level-3 fraction of Algorithm 5.2.2 is approximately $1 - O ( 1 / N )$ . See Bischof and Van Loan (1987) for further details.

# 5.2.4 Block Recursive QR

A more flexible approach to blocking involves recursion. Suppose $A \in \mathbb { R } ^ { m \times n }$ and assume for clarity that A has full column rank. Partition the thin QR factorization of A as follows:

$$
\left[ \begin{array}{c c} A _ {1} & A _ {2} \end{array} \right] = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \end{array} \right] \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & R _ {2 2} \end{array} \right].
$$

where $n _ { 1 } = \mathrm { H o o r } ( n / 2 ) , n _ { 2 } = n - n _ { 1 } , A _ { 1 } , Q _ { 1 } \in \mathbb { R } ^ { m \times n _ { 1 } }$ and $A _ { 2 } , Q _ { 2 } \in \mathbb { R } ^ { m \times n _ { 2 } }$ . From the equations $Q _ { 1 } R _ { 1 1 } = A _ { 1 } , R _ { 1 2 } = Q _ { 1 } ^ { T } A _ { 2 }$ , and $Q _ { 2 } R _ { 2 2 } = A _ { 2 } - Q _ { 1 } R _ { 1 2 }$ we obtain the following recursive procedure:

Algorithm 5.2.3 (Recursive Block QR) Suppose $A \in \mathbb { R } ^ { m \times n }$ has full column rank and $n _ { b }$ is a positive blocking parameter. The following algorithm computes $Q \in \mathbb { R } ^ { m \times n }$ with orthonormal columns and upper triangular $R \in \mathbb { R } ^ { n \times n }$ such that $A = Q R$ .

function [Q, R] = BlockQR(A, n, nb)

if $n \leq n _ { b }$

Use Algorithm 5.2.1 to compute the thin QR factorization $A = Q R$

$$
n _ {1} = \operatorname{floor} (n / 2)
$$

$$
[ Q _ {1}, R _ {1 1} ] = \text { BlockQR } (A (:, 1: n _ {1}), n _ {1}, n _ {b})
$$

$$
R _ {1 2} = Q _ {1} ^ {T} A (:, n _ {1} + 1: n)
$$

$$
A (:, n _ {1} + 1: n) = A (:, n _ {1} + 1: n) - Q _ {1} R _ {1 2}
$$

$$
[ Q _ {2}, R _ {2 2} ] = \operatorname{BlockQR} (A (:, n _ {1} + 1: n), n - n _ {1}, n _ {b})
$$

$$
Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \end{array} \right], R = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & R _ {2 2} \end{array} \right]
$$

end end

This divide-and-conquer approach is rich in matrix-matrix multiplication and provides a framework for the effective parallel computation of the QR factorization. See Elmroth and Gustavson (2001). Key implementation ideas concern the representation of the $Q -$ matrices and the incorporation of the §5.2.3 blocking strategies.

# 5.2.5 Givens QR Methods

Givens rotations can also be used to compute the QR factorization and the 4-by-3 case illustrates the general idea:

$$
\begin{array}{l} \left[ \begin{array}{c c c} \times & \times & \times \\ \times & \times & \times \\ \mathbf {X} & \times & \times \\ \mathbf {X} & \times & \times \end{array} \right] \xrightarrow {(3 , 4)} \left[ \begin{array}{c c c} \times & \times & \times \\ \mathbf {X} & \times & \times \\ \mathbf {X} & \times & \times \\ 0 & \times & \times \end{array} \right] \xrightarrow {(2 , 3)} \left[ \begin{array}{c c c} \mathbf {X} & \times & \times \\ \mathbf {X} & \times & \times \\ 0 & \times & \times \\ 0 & \times & \times \end{array} \right] \xrightarrow {(1 , 2)} \\ \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & \mathbf {x} & \times \\ 0 & \mathbf {x} & \times \end{array} \right] \stackrel {{(3, 4)}} {{\longrightarrow}} \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \mathbf {x} & \times \\ 0 & \mathbf {x} & \times \\ 0 & 0 & \times \end{array} \right] \stackrel {{(2, 3)}} {{\longrightarrow}} \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & 0 & \mathbf {x} \\ 0 & 0 & \mathbf {x} \end{array} \right] \stackrel {{(3, 4)}} {{\longrightarrow}} R. \\ \end{array}
$$

We highlighted the 2-vectors that define the underlying Givens rotations. If $G _ { j }$ denotes the jth Givens rotation in the reduction, then $Q ^ { \dot { T } } A \stackrel { - } { = } R$ is upper triangular, where $Q = G _ { 1 } \cdot \cdot \cdot G _ { t }$ and t is the total number of rotations. For general m and n we have:

Algorithm 5.2.4 (Givens QR) Given $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ , the following algorithm overwrites A with $Q ^ { T } A = R$ , where R is upper triangular and $Q$ is orthogonal.

for j = 1:n

for $i = m \colon - 1 { : } j + 1$

$$
[ c, s ] = \operatorname{givens} (A (i - 1, j), A (i, j))
$$

$$
A (i - 1: i, j: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A (i - 1: i, j: n)
$$

end end

This algorithm requires $3 n ^ { 2 } ( m - n / 3 )$ flops. Note that we could use the representation ideas from §5.1.11 to encode the Givens transformations that arise during the calculation. Entry $A ( i , j )$ can be overwritten with the associated representation.

With the Givens approach to the QR factorization, there is flexibility in terms of the rows that are involved in each update and also the order in which the zeros are introduced. For example, we can replace the inner loop body in Algorithm 5.2.4 with

$$
[ c, s ] = \operatorname{givens} (A (j, j), A (i, j))
$$

$$
A ([ j i ], j: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A ([ j i ], j: n)
$$

and still emerge with the QR factorization. It is also possible to introduce zeros by row. Whereas Algorithm 5.2.4 introduces zeros by column,

$$
\left[ \begin{array}{c c c} \times & \times & \times \\ 3 & \times & \times \\ 2 & 5 & \times \\ 1 & 4 & 6 \end{array} \right],
$$

the implementation

$$
\begin{array}{l} \text { for } i = 2: m \\ \text { for } j = 1: i - 1 \\ [ c, s ] = \operatorname{givens} (A (j, j), A (i, j)) \\ A ([ j i ], j: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A ([ j i ], j: n) \\ \end{array}
$$

introduces zeros by row, e.g.,

$$
\left[ \begin{array}{c c c} \times & \times & \times \\ 1 & \times & \times \\ 2 & 3 & \times \\ 4 & 5 & 6 \end{array} \right].
$$

# 5.2.6 Hessenberg QR via Givens

As an example of how Givens rotations can be used in a structured problem, we show how they can be employed to compute the QR factorization of an upper Hessenberg matrix. (Other structured QR factorizations are discussed in Chapter 6 and §11.1.8.) A small example illustrates the general idea. Suppose $n = 6$ and that after two steps we have computed

$$
G (2, 3, \theta_ {2}) ^ {T} G (1, 2, \theta_ {1}) ^ {T} A = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & \mathbf {x} & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right].
$$

Next, we compute $G ( 3 , 4 , \theta _ { 3 } )$ to zero the current (4,3) entry, thereby obtaining

$$
G (3, 4, \theta_ {3}) ^ {T} G (2, 3, \theta_ {2}) ^ {T} G (1, 2, \theta_ {1}) ^ {T} A = \left[ \begin{array}{l l l l l l} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right].
$$

Continuing in this way we obtain the following algorithm.

Algorithm 5.2.5 (Hessenberg QR) If $A \in \mathbb { R } ^ { n \times n }$ is upper Hessenberg, then the following algorithm overwrites A with $Q ^ { T } A = R$ where Q is orthogonal and R is upper triangular. $Q = G _ { 1 } \cdot \cdot \cdot G _ { n - 1 }$ is a product of Givens rotations where $G _ { j }$ has the form $G _ { j } = G ( j , j + 1 , \theta _ { j } )$ .

for $j = 1 \colon n - 1$

$$
[ c, s ] = \operatorname{givens} (A (j, j), A (j + 1, j))
$$

$$
A (j: j + 1, j: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A (j: j + 1, j: n)
$$

end

This algorithm requires about $3 n ^ { 2 }$ flops.

# 5.2.7 Classical Gram-Schmidt Algorithm

We now discuss two alternative methods that can be used to compute the thin QR factorization $A = Q _ { 1 } R _ { 1 }$ directly. If $\mathsf { r a n k } ( A ) = n$ , then equation (5.2.3) can be solved for $q _ { k }$ :

$$
q _ {k} = \left. \left(a _ {k} - \sum_ {i = 1} ^ {k - 1} r _ {i k} q _ {i}\right) \right/ r _ {k k}.
$$

Thus, we can think of $q _ { k }$ as a unit 2-norm vector in the direction of

$$
z _ {k} = a _ {k} - \sum_ {i = 1} ^ {k - 1} r _ {i k} q _ {i}
$$

where to ensure $z _ { k } \in \mathsf { s p a n } \{ q _ { 1 } , \dots , q _ { k - 1 } \} ^ { \perp }$ we choose

$$
r _ {i k} = q _ {i} ^ {T} a _ {k}, \qquad i = 1: k - 1.
$$

This leads to the classical Gram-Schmidt (CGS) algorithm for computing $A = Q _ { 1 } R _ { 1 }$ .

$$
R (1, 1) = \left\| A (:, 1) \right\| _ {2}
$$

$$
Q (:, 1) = A (:, 1) / R (1, 1)
$$

for $k = 2 { : } n$

$$
\begin{array}{l} R (1: k - 1, k) = Q (1: m, 1: k - 1) ^ {T} A (1: m, k) \\ z = A (1: m, k) - Q (1: m, 1: k - 1) \cdot R (1: k - 1, k) \\ R (k, k) = \parallel z \parallel_ {2} \\ Q (1: m, k) = z / R (k, k) \\ \end{array}
$$

end

In the kth step of CGS, the kth columns of both Q and R are generated.

# 5.2.8 Modified Gram-Schmidt Algorithm

Unfortunately, the CGS method has very poor numerical properties in that there is typically a severe loss of orthogonality among the computed $q _ { i }$ . Interestingly, a rearrangement of the calculation, known as modified Gram-Schmidt (MGS), leads to a more reliable procedure. In the kth step of MGS, the kth column of $Q$ (denoted by $q _ { k } )$ and the kth row of R (denoted by $r _ { k } ^ { T } )$ are determined. To derive the MGS method, define the matrix $A ^ { ( k ) } \in \mathbb { R } ^ { m \times ( n - k + 1 ) }$ by

$$
[ 0 \mid A ^ {(k)} ] = A - \sum_ {i = 1} ^ {k - 1} q _ {i} r _ {i} ^ {T} = \sum_ {i = k} ^ {n} q _ {i} r _ {i} ^ {T}.
$$

It follows that if

$$
A ^ {(k)} = \left[ \begin{array}{c c} z & B \\ 1 & n - k \end{array} \right]
$$

then $r _ { k k } ~ = ~ \parallel { z } \parallel _ { 2 } , q _ { k } ~ = ~ z / r _ { k k }$ , and $[ r _ { k , k + 1 } , . . . , r _ { k n } ] \ = \ q _ { k } ^ { T } B$ . We then compute the outer product $A ^ { ( k + 1 ) } = B ~ - ~ q _ { k } \left[ r _ { k , k + 1 } \cdot \cdot \cdot r _ { k n } \right]$ and proceed to the next step. This completely describes the kth step of MGS.

Algorithm 5.2.6 (Modified Gram-Schmidt) Given $A \in \mathbb { R } ^ { m \times n }$ with rank $( A ) = n$ , the following algorithm computes the thin QR factorization $A = Q _ { 1 } R _ { 1 }$ where $Q _ { 1 } \in \mathbb { R } ^ { m \times n }$ has orthonormal columns and $R _ { 1 } \in \mathbb { R } ^ { n \times n }$ is upper triangular.

$$
\begin{array}{l} R (k, k) = \left\| A (1: m, k) \right\| _ {2} \\ Q (1: m, k) = A (1: m, k) / R (k, k) \\ R (k, j) = Q (1: m, k) ^ {T} A (1: m, j) \\ A (1: m, j) = A (1: m, j) - Q (1: m, k) R (k, j) \\ \end{array}
$$

This algorithm requires $2 m n ^ { 2 }$ flops. It is not possible to overwrite A with both $Q _ { 1 }$ and $R _ { 1 }$ . Typically, the MGS computation is arranged so that A is overwritten by $Q _ { 1 }$ and the matrix $R _ { 1 }$ is stored in a separate array.

# 5.2.9 Work and Accuracy

If one is interested in computing an orthonormal basis for ran(A), then the Householder approach requires $2 m n ^ { 2 } - 2 n ^ { 3 } / 3$ flops to get $Q$ in factored form and another $2 m n ^ { 2 } -$ $2 n ^ { 3 } / 3$ flops to get the first n columns of Q. (This requires “paying attention” to just the first n columns of $Q$ in (5.1.5).) Therefore, for the problem of finding an orthonormal basis for ran(A), MGS is about twice as efficient as Householder orthogonalization. However, Bj¨orck (1967) has shown that MGS produces a computed ${ \hat { Q } } _ { 1 } = \left[ { \hat { q } } _ { 1 } | \cdots | { \hat { q } } _ { n } \right]$ that satisfies

$$
\hat {Q} _ {1} ^ {T} \hat {Q} _ {1} = I + E _ {M G S}, \quad \| E _ {M G S} \| _ {2} \approx \mathbf {u} \kappa_ {2} (A),
$$

whereas the corresponding result for the Householder approach is of the form

$$
\hat {Q} _ {1} ^ {T} \hat {Q} _ {1} = I + E _ {H}, \quad \| E _ {H} \| _ {2} \approx \mathbf {u}.
$$

Thus, if orthonormality is critical, then MGS should be used to compute orthonormal bases only when the vectors to be orthogonalized are fairly independent.

We also mention that the computed triangular factor $\hat { R }$ produced by MGS satisfies $\Vert A - \hat { Q } \hat { R } \Vert \approx \mathbf { u } \Vert A \Vert$ and that there exists a Q with perfectly orthonormal columns such that $\parallel A - Q { \hat { R } } \parallel \approx \mathbf { u } \parallel A \parallel$ . See Higham (ASNA, p. 379) and additional references given at the end of this section.

# 5.2.10 A Note on Complex Householder QR

Complex Householder transformations (§5.1.13) can be used to compute the QR factorization of a complex matrix $A \in \mathbb { C } ^ { m \times n }$ . Analogous to Algorithm 5.2.1 we have

for $j = 1 { : } n$

Compute a Householder matrix $Q _ { j }$ so that $Q _ { j } A$ is upper triangular through its first j columns.

$$
A = Q _ {j} A
$$

end

Upon termination, A has been reduced to an upper triangular matrix $R \in \mathbb { C } ^ { m \times n }$ and we have $A = Q R$ where $Q = Q _ { 1 } \cdot \cdot \cdot Q _ { n }$ is unitary. The reduction requires about four times the number of flops as the real case.

# Problems

P5.2.1 Adapt the Householder QR algorithm so that it can efficiently handle the case when $A \in \mathbb { R } ^ { m \times n }$ has lower bandwidth p and upper bandwidth q.

P5.2.2 Suppose $A \in \mathbb { R } ^ { n \times n }$ and let E be the exchange permutation $\mathcal { E } _ { n }$ obtained by reversing the order of the rows in ${ { I } _ { n } } . \mathrm { { \Gamma } ( a ) }$ Show that if $R \in \mathbb { R } ^ { n \times n }$ is upper triangular, then $L = \mathcal { E } R \mathcal { \bar { E } }$ is lower triangular. (b) Show how to compute an orthogonal $Q \in \mathbb { R } ^ { n \times n }$ and a lower triangular $\ b { L } \in \mathbb { R } ^ { n \times n }$ so that $A = Q L$ assuming the availability of a procedure for computing the QR factorization.

P5.2.3 Adapt the Givens QR factorization algorithm so that the zeros are introduced by diagonal. That is, the entries are zeroed in the order $( m , 1 ) , ( m - 1 , 1 ) , ( m , 2 ) , ( m - 2 , 1 ) , ( m - 1 , 2 ) , ( m , 3 )$ , etc.

P5.2.4 Adapt the Givens QR factorization algorithm so that it efficiently handles the case when A is n-by-n and tridiagonal. Assume that the subdiagonal, diagonal, and superdiagonal of A are stored in $e ( 1 { : } n - 1 ) , a ( 1 { : } n ) , f ( 1 { : } n - 1 )$ , respectively. Design your algorithm so that these vectors are overwritten by the nonzero portion of T .

P5.2.5 Suppose $\boldsymbol { L } \in \mathbb { R } ^ { m \times n }$ with $m \ \geq \ n$ is lower triangular. Show how Householder matrices $H _ { 1 } , \ldots , H _ { n }$ can be used to determine a lower triangular $L _ { 1 } \in \mathbb { R } ^ { n \times n }$ so that

$$
H _ {n} \dots H _ {1} L = \left[ \begin{array}{c} L _ {1} \\ 0 \end{array} \right].
$$

Hint: The second step in the 6-by-3 case involves finding $H _ { 2 }$ so that

$$
H _ {2} \left[ \begin{array}{c c c} \times & 0 & 0 \\ \times & \times & 0 \\ \times & \times & \times \\ \times & \times & 0 \\ \times & \times & 0 \\ \times & \times & 0 \end{array} \right] = \left[ \begin{array}{c c c} \times & 0 & 0 \\ \times & \times & 0 \\ \times & \times & \times \\ \times & 0 & 0 \\ \times & 0 & 0 \\ \times & 0 & 0 \end{array} \right]
$$

with the property that rows 1 and 3 are left alone.

P5.2.6 Suppose $A \in \mathbb { R } ^ { n \times n }$ and $D = \operatorname { d i a g } ( d _ { 1 } , \ldots , d _ { n } ) \in \mathbb { R } ^ { n \times n }$ . Show how to construct an orthogonal Q such that

$$
Q ^ {T} A - D Q ^ {T} = R
$$

is upper triangular. Do not worry about efficiency—this is just an exercise in QR manipulation.

P5.2.7 Show how to compute the QR factorization of the product

$$
A = A _ {p} \dots A _ {2} A _ {1}
$$

without explicitly multiplying the matrices $A _ { 1 } , \dotsc , A _ { p }$ together. Assume that each $A _ { i }$ is square. Hint: In the $p = 3 \ \mathrm { c a s e }$ , write

$$
Q _ {3} ^ {T} A = Q _ {3} ^ {T} A _ {3} Q _ {2} Q _ {2} ^ {T} A _ {2} Q _ {1} Q _ {1} ^ {T} A _ {1}
$$

and determine orthogonal $Q _ { i }$ so that $Q _ { i } ^ { T } ( A _ { i } Q _ { i - 1 } )$ is upper triangular. $( Q _ { 0 } = I . )$

P5.2.8 MGS applied to $A \in \mathbb { R } ^ { m \times n }$ is numerically equivalent to the first step in Householder QR applied to

$$
\tilde {A} = \left[ \begin{array}{c} O _ {n} \\ A \end{array} \right]
$$

where $O _ { n }$ is the $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ zero matrix. Verify that this statement is true after the first step of each method is completed.

P5.2.9 Reverse the loop orders in Algorithm 5.2.6 (MGS) so that R is computed column by column.

P5.2.10 How many flops are required by the complex QR factorization procedure outlined in §5.10?

P5.2.11 Develop a complex version of the Givens QR factorization in which the diagonal of R is nonnegative. See §5.1.13.

P5.2.12 Show that if $A \in \mathbb { R } ^ { n \times n }$ and $a _ { i } = A ( : , i )$ , then

$$
| \det (A) | \leq \| a _ {1} \| _ {2} \dots \| a _ {n} \| _ {2}.
$$

Hint: Use the QR factorization.

P5.2.13 Suppose $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ Construct an orthogonal $Q \in \mathbb { R } ^ { ( m + n ) \times ( m + n ) }$ with the property that Q(1:m, 1:n) is a scalar multiple of A. Hint. If $\alpha \in \mathbb { R }$ is chosen properly, then $I - \alpha ^ { 2 } A ^ { T } A$ has a Cholesky factorization.

P5.2.14 Suppose $A \in \mathbb { R } ^ { m \times n }$ . Analogous to Algorithm 5.2.4, show how fast Givens transformations (P5.1.12) can be used to compute $\bar { M } \in \mathbb { R } ^ { m \times m }$ and a diagonal $D \in \mathbb { R } ^ { m \times m }$ with positive diagonal entries so that $M ^ { T } A = S$ is upper triangular and $M M ^ { T } = \bar { D }$ . Relate M and S to A’s QR factors.

P5.2.15 (Parallel Givens QR) Suppose $A \in \mathbb { R } ^ { 9 \times 3 }$ and that we organize a Givens QR so that the subdiagonal entries are zeroed over the course of ten “time steps” as follows:

<table><tr><td>Step</td><td colspan="3">Entries Zeroed</td></tr><tr><td>T=1</td><td>(9,1)</td><td></td><td></td></tr><tr><td>T=2</td><td>(8,1)</td><td></td><td></td></tr><tr><td>T=3</td><td>(7,1)</td><td>(9,2)</td><td></td></tr><tr><td>T=4</td><td>(6,1)</td><td>(8,2)</td><td></td></tr><tr><td>T=5</td><td>(5,1)</td><td>(7,2)</td><td>(9,3)</td></tr><tr><td>T=6</td><td>(4,1)</td><td>(6,2)</td><td>(8,3)</td></tr><tr><td>T=7</td><td>(3,1)</td><td>(5,2)</td><td>(7,3)</td></tr><tr><td>T=8</td><td>(2,1)</td><td>(4,2)</td><td>(6,3)</td></tr><tr><td>T=9</td><td></td><td>(3,2)</td><td>(5,3)</td></tr><tr><td>T=10</td><td></td><td></td><td>(4,3)</td></tr></table>

Assume that a rotation in plane $( i - 1 , i )$ is used to zero a matrix entry $( i , j )$ . It follows that the rotations associated with any given time step involve disjoint pairs of rows and may therefore be computed in parallel. For example, during time step $T = 6 ,$ , there is a $( 3 , 4 ) , ( 5 , 6 )$ , and $^ { ( 7 , 8 ) }$ rotation. Three separate processors could oversee the three updates. Extrapolate from this example to the m-by-n case and show how the QR factorization could be computed in $O ( m + n )$ time steps. How many of those time steps would involve n “nonoverlapping” rotations?

# Notes and References for §5.2

The idea of using Householder transformations to solve the least squares problem was proposed in:

A.S. Householder (1958). “Unitary Triangularization of a Nonsymmetric Matrix,” J. ACM 5, 339–342.   
The practical details were worked out in:   
P. Businger and G.H. Golub (1965). “Linear Least Squares Solutions by Householder Transformations,” Numer. Math. 7, 269–276.   
G.H. Golub (1965). “Numerical Methods for Solving Linear Least Squares Problems,” Numer. Math. 7, 206–216.   
The basic references for Givens QR include:   
W. Givens (1958). “Computation of Plane Unitary Rotations Transforming a General Matrix to Triangular Form,” SIAM J. Appl. Math. 6, 26–50.   
M. Gentleman (1973). “Error Analysis of QR Decompositions by Givens Transformations,” Lin. Alg. Applic. 10, 189–197.   
There are modifications for the QR factorization that make it more attractive when dealing with rank deficiency. See §5.4. Nevertheless, when combined with the condition estimation ideas in §3.5.4, the traditional QR factorization can be used to address rank deficiency issues:   
L.V. Foster (1986). “Rank and Null Space Calculations Using Matrix Decomposition without Column Interchanges,” Lin. Alg. Applic. 74, 47–71.   
The behavior of the Q and R factors when A is perturbed is of interest. A main result is that the resulting changes in Q and R are bounded by the condition of A times the relative change in A, see:   
G.W. Stewart (1977). “Perturbation Bounds for the QR Factorization of a Matrix,” SIAM J. Numer. Anal. 14, 509–518.   
H. Zha (1993). “A Componentwise Perturbation Analysis of the QR Decomposition,” SIAM J. Matrix Anal. Applic. 4, 1124–1131.   
G.W. Stewart (1993). “On the Perturbation of LU Cholesky, and QR Factorizations,” SIAM J. Matrix Anal. Applic. 14, 1141–1145.   
A. Barrlund (1994). “Perturbation Bounds for the Generalized QR Factorization,” Lin. Alg. Applic. 207, 251–271.   
J.-G. Sun (1995). “On Perturbation Bounds for the QR Factorization,” Lin. Alg. Applic. 215, 95–112.   
X.-W. Chang and C.C. Paige (2001). “Componentwise Perturbation Analyses for the QR factorization,” Numer. Math. 88, 319–345.   
Organization of the computation so that the entries in Q depend continuously on the entries in A is discussed in:   
T.F. Coleman and D.C. Sorensen (1984). “A Note on the Computation of an Orthonormal Basis for the Null Space of a Matrix,” Mathematical Programming 29, 234–242.   
References for the Gram-Schmidt process and various ways to overcome its shortfalls include:   
J.R. Rice (1966). “Experiments on Gram-Schmidt Orthogonalization,” Math. Comput. 20, 325–328.   
A. Bj¨orck (1967). “Solving Linear Least Squares Problems by Gram-Schmidt Orthogonalization,” BIT 7, 1–21.   
N.N. Abdelmalek (1971). “Roundoff Error Analysis for Gram-Schmidt Method and Solution of Linear Least Squares Problems,” BIT 11, 345–368.   
A. Ruhe (1983). “Numerical Aspects of Gram-Schmidt Orthogonalization of Vectors,” Lin. Alg. Applic. 52/53, 591–601.   
W. Jalby and B. Philippe (1991). “Stability Analysis and Improvement of the Block Gram-Schmidt Algorithm,” SIAM J. Sci. Stat. Comput. 12, 1058–1073.   
A. Bj¨ ˚ orck and C.C. Paige (1992). “Loss and Recapture of Orthogonality in the Modified Gram-Schmidt Algorithm,” SIAM J. Matrix Anal. Applic. 13, 176–190.   
A. Bj¨orck (1994). “Numerics of Gram-Schmidt Orthogonalization,” Lin. Alg. Applic. 197/198, 297–316.   
L. Giraud and J. Langou (2003). “A Robust Criterion for the Modified Gram-Schmidt Algorithm with Selective Reorthogonalization,” SIAM J. Sci. Comput. 25, 417–441.   
G.W. Stewart (2005). “Error Analysis of the Quasi-Gram–Schmidt Algorithm,” SIAM J. Matrix Anal. Applic. 27, 493–506.

L. Giraud, J. Langou, M. Rozlonk, and J. van den Eshof (2005). “Rounding Error Analysis of the Classical Gram-Schmidt Orthogonalization Process,” Numer. Math. 101, 87–100.   
A. Smoktunowicz, J.L. Barlow and J. Langou (2006). “A Note on the Error Analysis of Classical Gram-Schmidt,” Numer. Math. 105, 299–313.   
Various high-performance issues pertaining to the QR factorization are discussed in:   
B. Mattingly, C. Meyer, and J. Ortega (1989). “Orthogonal Reduction on Vector Computers,” SIAM J. Sci. Stat. Comput. 10, 372–381.   
P.A. Knight (1995). “Fast Rectangular Matrix Multiplication and the QR Decomposition,” Lin. Alg. Applic. 221, 69–81.   
J.J. Carrig, Jr. and G.L. Meyer (1997). “Efficient Householder QR Factorization for Superscalar Processors,” ACM Trans. Math. Softw. 23, 362–378.   
D. Vanderstraeten (2000). “An Accurate Parallel Block Gram-Schmidt Algorithm without Reorthogonalization,” Numer. Lin. Alg. 7, 219–236.   
E. Elmroth and F.G. Gustavson (2000). “Applying Recursion to Serial and Parallel QR Factorization Leads to Better Performance,” IBM J. Res. Dev. 44, 605–624.   
Many important high-performance implementation ideas apply equally to LU, Cholesky, and QR, see:   
A. Buttari, J. Langou, J. Kurzak, and J. Dongarra (2009). “A Class of Parallel Tiled Linear Algebra Algorithms for Multicore Architectures,” Parallel Comput. 35, 38–53.   
J. Kurzak, H. Ltaief, and J. Dongarra (2010). “Scheduling Dense Linear Algebra Operations on Multicore Processors,” Concurrency Comput. Pract. Exper. 22, 15–44.   
J. Demmel, L. Grigori, M, Hoemmen, and J. Langou (2012). “Methods and Algorithms for Scientific Computing Communication-optimal Parallel and Sequential QR and LU Factorizations,” SIAM J. Sci. Comput. 34, A206-A239.   
Historical references concerned with parallel Givens QR include:   
W.M. Gentleman and H.T. Kung (1981). “Matrix Triangularization by Systolic Arrays,” SPIE Proc. 298, 19–26.   
D.E. Heller and I.C.F. Ipsen (1983). “Systolic Networks for Orthogonal Decompositions,” SIAM J. Sci. Stat. Comput. 4, 261–269.   
M. Costnard, J.M. Muller, and Y. Robert (1986). “Parallel QR Decomposition of a Rectangular Matrix,” Numer. Math. 48, 239–250.   
L. Eldin and R. Schreiber (1986). “An Application of Systolic Arrays to Linear Discrete Ill-Posed Problems,” SIAM J. Sci. Stat. Comput. 7, 892–903.   
F.T. Luk (1986). “A Rotation Method for Computing the QR Factorization,” SIAM J. Sci. Stat. Comput. 7, 452–459.   
J.J. Modi and M.R.B. Clarke (1986). “An Alternative Givens Ordering,” Numer. Math. 43, 83–90.   
The QR factorization of a structured matrix is usually structured itself, see:   
A.W. Bojanczyk, R.P. Brent, and F.R. de Hoog (1986). “QR Factorization of Toeplitz Matrices,” Numer. Math. 49, 81–94.   
S. Qiao (1986). “Hybrid Algorithm for Fast Toeplitz Orthogonalization,” Numer. Math. 53, 351–366.   
C.J. Demeure (1989). “Fast QR Factorization of Vandermonde Matrices,” Lin. Alg. Applic. 122/123/124, 165–194.   
L. Reichel (1991). “Fast QR Decomposition of Vandermonde-Like Matrices and Polynomial Least Squares Approximation,” SIAM J. Matrix Anal. Applic. 12, 552–564.   
D.R. Sweet (1991). “Fast Block Toeplitz Orthogonalization,” Numer. Math. 58, 613–629.   
Quantum computation has an interesting connection to complex Givens rotations and their application to vectors, see:   
G. Cybenko (2001). “Reducing Quantum Computations to Elementary Unitary Transformations,” Comput. Sci. Eng. 3, 27–32.   
D.P. O’Leary and S.S. Bullock (2005). “QR Factorizations Using a Restricted Set of Rotations,” ETNA 21, 20–27.   
N.D. Mermin (2007). Quantum Computer Science, Cambridge University Press, New York.
