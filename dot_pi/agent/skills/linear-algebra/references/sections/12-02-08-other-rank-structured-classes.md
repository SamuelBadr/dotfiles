# 12.2.8 Other Rank-Structured Classes

We briefly mention several other rank structures that arise in applications. Fast LU and QR procedures exist in each case.

If p and q are nonnegative integers, then a matrix A is $\{ p , q \}$ -semiseparable if

$$
j _ {2} <   i _ {1} + p \Rightarrow \operatorname{rank} \left(A \left(i _ {1}: i _ {2}, j _ {1}: j _ {2}\right)\right) \leq p,
$$

$$
i _ {2} > j _ {1} + q \Rightarrow \operatorname{rank} \left(A \left(i _ {1}: i _ {2}, j _ {1}: j _ {2}\right)\right) \leq q.
$$

For example, if A is {2, 3}-semiseparable, then

$$
A = \left[ \begin{array}{c c c c c c c} \times & \times & \times & \times & \times & \times & \times \\ a _ {2 1} & a _ {2 2} & a _ {2 3} & \times & \times & \times & \times \\ a _ {3 1} & a _ {3 2} & a _ {3 3} & a _ {3 4} & a _ {3 5} & a _ {3 6} & a _ {3 7} \\ a _ {4 1} & a _ {4 2} & a _ {4 3} & a _ {4 4} & a _ {4 5} & a _ {4 6} & a _ {4 7} \\ \times & \times & \times & a _ {5 4} & a _ {5 5} & a _ {5 6} & a _ {5 7} \\ \times & \times & \times & a _ {6 4} & a _ {6 5} & a _ {6 6} & a _ {6 7} \\ \times & \times & \times & a _ {7 4} & a _ {7 5} & a _ {7 6} & a _ {7 7} \end{array} \right] \Rightarrow \quad \begin{array}{l} \operatorname{rank} (A (2: 4, 1: 3)) \leq 2, \\ \operatorname{rank} (A (3: 7, 4: 7)) \leq 3. \end{array}
$$

In general, A is $\{ p , q \}$ -generator representable if we have U, $V \in \mathbb { R } ^ { n \times p }$ and $P , Q \in \mathbb { R } ^ { n \times q }$ such that

$$
\operatorname{tril} (A, p - 1) = \operatorname{tril} (U V ^ {T}, p - 1),
$$

$$
\operatorname{triu} (A, - q + 1) = \operatorname{triu} \left(P Q ^ {T}, - q + 1\right).
$$

If such a matrix is nonsingular, then $A ^ { - 1 }$ has lower bandwidth p and upper bandwidth q. If the $\{ p , q \}$ -semiseparable definition is modified so that the rank-p blocks come from tril(A) and the rank-q blocks come from triu(A), then A belongs to the class of extended $\{ p , q \} { - } s e p a r a b l e$ matrices. If the $\{ p , q \}$ -semiseparable definition is modified so that the rank-p blocks come from tri $( A , - 1 )$ and the rank-q come from triu(A, 1), then A belongs to the class of extended $\{ p , q \}$ -quasiseparable matrices. A sequentially semiseparable matrix is a block matrix that has the following form:

$$
A = \left[ \begin{array}{c c c c} D _ {1} & P _ {1} Q _ {2} ^ {T} & P _ {1} R _ {2} Q _ {3} ^ {T} & P _ {1} R _ {2} R _ {3} Q _ {4} ^ {T} \\ U _ {2} V _ {1} ^ {T} & D _ {2} & P _ {2} Q _ {3} ^ {T} & P _ {2} R _ {3} Q _ {4} ^ {T} \\ U _ {3} T _ {2} V _ {1} ^ {T} & U _ {3} V _ {2} ^ {T} & D _ {3} & P _ {3} Q _ {4} ^ {T} \\ U _ {4} T _ {3} T _ {2} V _ {1} ^ {T} & U _ {4} T _ {3} V _ {2} ^ {T} & U _ {4} V _ {3} ^ {T} & D _ {4} \end{array} \right]. \tag {12.2.17}
$$

See Dewilde and van der Veen (1997) and Chandrasekaran et al. (2005). The blocks can be rectangular so least squares problems with this structure can be handled.

Matrices with hierarchical rank structure are based on low-rank patterns that emerge through recursive 2-by-2 blockings. (With one level of recursion we would have 2-by-2 block matrix whose diagonal blocks are 2-by-2 block matrices.) Various connections may exist between the low-rank representations of the off-diagonal blocks. The important class of hierarchically semiseparable matrices has a particularly rich and exploitable structure; see Xia (2012).
