# 7.7.4 Hessenberg-Triangular Form

The first step in computing the generalized real Schur decomposition of the pair $( A , B )$ is to reduce A to upper Hessenberg form and B to upper triangular form via orthogonal transformations. We first determine an orthogonal U such that $U ^ { T } B$ is upper triangular. Of course, to preserve eigenvalues, we must also update A in exactly the same way. Let us trace what happens in the $n = 5 ~ \mathrm { c a s e }$ .

$$
A \gets U ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \end{array} \right], \quad B \gets U ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

Next, we reduce A to upper Hessenberg form while preserving $B ^ { \prime } \mathrm { s }$ upper triangular form. First, a Givens rotation $Q _ { 4 5 }$ is determined to zero $a _ { 5 1 }$ :

$$
A \gets Q _ {4 5} ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], \quad B \gets Q _ {4 5} ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right].
$$

The nonzero entry arising in the (5,4) position in B can be zeroed by postmultiplying with an appropriate Givens rotation $Z _ { 4 5 }$ :

$$
A \gets A Z _ {4 5} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], B \gets B Z _ {4 5} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

Zeros are similarly introduced into the (4, 1) and (3, 1) positions in A:

$$
A \gets Q _ {3 4} ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], B \gets Q _ {3 4} ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
A \gets A Z _ {3 4} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], \quad B \gets B Z _ {3 4} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
A \gets Q _ {2 3} ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], \quad B \gets Q _ {2 3} ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
A \gets A Z _ {2 3} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], B \gets B Z _ {2 3} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

A is now upper Hessenberg through its first column. The reduction is completed by zeroing $a _ { 5 2 } , ~ a _ { 4 2 }$ , and $a _ { 5 3 }$ . Note that two orthogonal transformations are required for each $a _ { i j }$ that is zeroed—one to do the zeroing and the other to restore $B \mathrm { { ^ { * } s } }$ triangularity. Either Givens rotations or 2-by-2 modified Householder transformations can be used. Overall we have:

Algorithm 7.7.1 (Hessenberg-Triangular Reduction) Given A and B in $\mathbb { R } ^ { n \times n }$ , the following algorithm overwrites A with an upper Hessenberg matrix $Q ^ { T } A Z$ and $B$ with an upper triangular matrix $Q ^ { T } B Z$ where both $Q$ and $Z$ are orthogonal.

Compute the factorization $B = Q R$ using Algorithm 5.2.1 and overwrite

A with $Q ^ { T } A$ and B with $Q ^ { T } B$

for $j = 1 { : } n - 2$

for $i = n \colon - 1 { : } j + 2$

$$
[ c, s ] = \text {   givens   } (A (i - 1, j), A (i, j))
$$

$$
A (i - 1: i, j: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A (i - 1: i, j: n)
$$

$$
B (i - 1: i, i - 1: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} B (i - 1: i, i - 1: n)
$$

$$
[ c, s ] = \text {   givens } (- B (i, i), B (i, i - 1))
$$

$$
B (1: i, i - 1: i) = B (1: i, i - 1: i) \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right]
$$

$$
A (1: n, i - 1: i) = A (1: n, i - 1: i) \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right]
$$

end

end

This algorithm requires about $8 n ^ { 3 }$ flops. The accumulation of $Q$ and $Z$ requires about $4 n ^ { 3 }$ and $3 n ^ { 3 }$ flops, respectively.

The reduction of $A - \lambda B$ to Hessenberg-triangular form serves as a “front end” decomposition for a generalized QR iteration known as the $\mathrm { Q Z }$ iteration which we describe next.
