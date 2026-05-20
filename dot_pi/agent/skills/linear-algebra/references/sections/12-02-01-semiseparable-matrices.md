# 12.2.1 Semiseparable Matrices

A matrix $A \in \mathbb { R } ^ { n \times n }$ is semiseparable if every block that does not “cross” the diagonal has unit rank or less. This means

$$
j _ {2} \leq i _ {1} \text {   or   } i _ {2} \leq j _ {1} \Rightarrow \operatorname{rank} (A (i _ {1}: i _ {2}, j _ {1}: j _ {2})) \leq 1. \tag {12.2.1}
$$

The rank-1 blocks of interest in a semiseparable matrix are wholly contained in either its upper triangular part or its lower triangular part, e.g.,

$$
\left[ \begin{array}{c c c c c c} \times & \times & a _ {1 3} & a _ {1 4} & \times & \times \\ \times & \times & a _ {2 3} & a _ {2 4} & \times & \times \\ \times & \times & a _ {3 3} & a _ {3 4} & \times & \times \\ \times & \times & \times & \times & \times & \times \\ a _ {5 1} & a _ {5 2} & \times & \times & \times & \times \\ a _ {6 1} & a _ {6 2} & \times & \times & \times & \times \end{array} \right], \qquad \begin{array}{l} \operatorname{rank} (A (1: 3, 3: 4)) \leq 1, \\ \operatorname{rank} (A (5: 6, 1: 2)) \leq 1. \end{array}
$$

Semiseparable matrices are data-sparse and enormous savings can be realized when their structure is exploited. For example, we will show that the factorizations $A = L U$ and $A = Q R$ for semiseparable A require just $O ( n )$ flops to compute and $O ( n )$ flops to represent.

An important example of a semiseparable matrix is the inverse of a unit bidiagonal matrix. Given $r \in \mathbb { R } ^ { n - 1 }$ we define $B ( r ) \in \mathbb { R } ^ { n \times n }$ by

$$
B (r) = \left[ \begin{array}{c c c c c} 1 & - r _ {1} & 0 & 0 & 0 \\ 0 & 1 & - r _ {2} & 0 & 0 \\ 0 & 0 & 1 & - r _ {3} & 0 \\ 0 & 0 & 0 & 1 & - r _ {4} \\ 0 & 0 & 0 & 0 & 1 \end{array} \right]. \tag {12.2.2}
$$

Observe that any submatrix extracted from the upper triangular portion of

$$
B (r) ^ {- 1} = \left[ \begin{array}{c c c c c} 1 & r _ {1} & r _ {1} r _ {2} & r _ {1} r _ {2} r _ {3} & r _ {1} r _ {2} r _ {3} r _ {4} \\ 0 & 1 & r _ {2} & r _ {2} r _ {3} & r _ {2} r _ {3} r _ {4} \\ 0 & 0 & 1 & r _ {3} & r _ {3} r _ {4} \\ 0 & 0 & 0 & 1 & r _ {4} \\ 0 & 0 & 0 & 0 & 1 \end{array} \right] \tag {12.2.3}
$$

has unit rank. If $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $r = x ( 2 { : } n ) \cdot / x ( 1 { : } n - 1 )$ is defined, then

$$
B (r) ^ {T} x = x _ {1} e _ {1}.
$$

Thus, the matrix $B ( r )$ can (in principle) be used to introduce zeros into a vector.
