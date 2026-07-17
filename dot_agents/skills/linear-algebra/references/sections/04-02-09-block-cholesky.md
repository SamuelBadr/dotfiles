# 4.2.9 Block Cholesky

Just as there are block methods for computing the LU factorization, so are there are block methods for computing the Cholesky factorization. Paralleling the derivation of the block LU algorithm in §3.2.11, we start by blocking $A = G G ^ { T }$ as follows

$$
\left[ \begin{array}{l l} A _ {1 1} & A _ {2 1} ^ {T} \\ A _ {2 1} & A _ {2 2} \end{array} \right] = \left[ \begin{array}{c c} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right] \left[ \begin{array}{c c} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right] ^ {T}. \tag {4.2.18}
$$

Here, $A _ { 1 1 } \in \mathbb { R } ^ { r \times r } , A _ { 2 2 } \in \mathbb { R } ^ { ( n - r ) \times ( n - r ) }$ , r is a blocking parameter, and G is partitioned conformably. Comparing blocks in (4.2.18) we conclude that

$$
A _ {1 1} = G _ {1 1} G _ {1 1} ^ {T},
$$

$$
A _ {2 1} = G _ {2 1} G _ {1 1} ^ {T},
$$

$$
A _ {2 2} = G _ {2 1} G _ {2 1} ^ {T} + G _ {2 2} G _ {2 2} ^ {T},
$$

which suggests the following 3-step procedure:

Step 1: Compute the Cholesky factorization of $A _ { 1 1 }$ to get $G _ { 1 1 }$ .

Step 2: Solve a lower triangular multiple-right-hand-side system for $G _ { 2 1 }$ .

Step 3: Compute the Cholesky factor $G _ { 2 2 }$ of $A _ { 2 2 } - G _ { 2 1 } G _ { 2 1 } ^ { T } = A _ { 2 2 } - A _ { 2 1 } A _ { 1 1 } ^ { - 1 } A _ { 2 1 } ^ { T }$ . In recursive form we obtain the following algorithm.

Algorithm 4.2.3 (Recursive Block Cholesky) Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite and r is a positive integer. The following algorithm computes a lower triangular $G \in \mathbb { R } ^ { n \times n }$ so $\bar { A } = G G ^ { T }$ .

function $G = { \mathsf { B l o c k C h o l e s k y } } ( A , n , r )$

if $n \leq r$

Compute the Cholesky factorization $A = G G ^ { T }$ .

else

Compute the Cholesky factorization $A ( 1 { : } r , 1 { : } r ) = G _ { 1 1 } G _ { 1 1 } ^ { T }$

Solve $G _ { 2 1 } G _ { 1 1 } ^ { T } = A ( r + 1 { : } n , 1 { : } r )$ for $G _ { 2 1 }$

$$
\tilde {A} = A (r + 1: n, r + 1: n) - G _ {2 1} G _ {2 1} ^ {T}
$$

$$
G _ {2 2} = \text { BlockCholesky } (\tilde {A}, n - r, r)
$$

$$
G = \left[ \begin{array}{c c} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right]
$$

end

If symmetry is exploited in the computation of ${ \tilde { A } } ,$ then this algorithm requires $n ^ { 3 } / 3$ flops. A careful accounting of flops reveals that the level-3 fraction is about $1 - 1 / N ^ { 2 }$ where $N \approx n / r$ . The “small” Cholesky computation for $G _ { 1 1 }$ and the “thin” solution process for $G _ { 2 1 }$ are dominated by the “large” level-3 update for $\tilde { A }$ .

To develop a nonrecursive implementation, we assume for clarity that $n = N r$ where N is a positive integer and consider the partitioning

$$
\left[ \begin{array}{c c c} A _ {1 1} & \dots & A _ {1 N} \\ \vdots & \ddots & \vdots \\ A _ {N 1} & \dots & A _ {N N} \end{array} \right] = \left[ \begin{array}{c c c} G _ {1 1} & \dots & 0 \\ \vdots & \ddots & \vdots \\ G _ {N 1} & \dots & G _ {N N} \end{array} \right] \left[ \begin{array}{c c c} G _ {1 1} & \dots & 0 \\ \vdots & \ddots & \vdots \\ G _ {N 1} & \dots & G _ {N N} \end{array} \right] ^ {T} \tag {4.2.19}
$$

where all blocks are r-by-r. By equating $( i , j )$ blocks with $i \geq j$ it follows that

$$
A _ {i j} = \sum_ {k = 1} ^ {j} G _ {i k} G _ {j k} ^ {T}.
$$

Define

$$
S = A _ {i j} - \sum_ {k = 1} ^ {j - 1} G _ {i k} G _ {j k} ^ {T} = A _ {i j} - \left[ G _ {i 1} \left| \dots \right| G _ {i, j - 1} \right] \left[ \begin{array}{c} G _ {j 1} ^ {T} \\ \vdots \\ G _ {j, j - 1} ^ {T} \end{array} \right].
$$

If $i = j$ , then $G _ { j j }$ is the Cholesky factor of S. If $i > j$ , then $G _ { i j } G _ { j j } ^ { T } = S$ and $G _ { i j }$ is the solution to a triangular multiple right hand side problem. Properly sequenced, these equations can be arranged to compute all the G-blocks.

Algorithm 4.2.4 (Nonrecursive Block Cholesky) Given a symmetric positive definite $A \in \mathbb { R } ^ { n \times n }$ with $n = N r$ with blocking (4.2.19), the following algorithm computes a lower triangular $G \in \mathbb { R } ^ { n \times n }$ such that $A { \stackrel { \cdot } { = } } G G ^ { T }$ . The lower triangular part of A is overwritten by the lower triangular part of G.

for $j = 1:N$ for $i = j:N$ Compute $S = A_{ij} - \sum_{k=1}^{j-1} G_{ik} G_{jk}^T$ .
    if $i = j$ Compute Cholesky factorization $S = G_{jj} G_{jj}^T$ .
    else
    Solve $G_{ij} G_{jj}^T = S$ for $G_{ij}$ .
    end $A_{ij} = G_{ij}$ .
    end
end

The overall process involves $n ^ { 3 } / 3$ flops like the other Cholesky procedures that we have developed. The algorithm is rich in matrix multiplication with a level-3 fraction given by $1 - \left( 1 / N ^ { 2 } \right)$ . The algorithm can be easily modified to handle the case when r does not divide n.
