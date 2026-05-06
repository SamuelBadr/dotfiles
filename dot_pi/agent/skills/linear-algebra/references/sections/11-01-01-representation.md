# 11.1.1 Representation

Data structures play an important role in sparse matrix computations. Typically, a real vector is used to house the nonzero entries of the matrix and one or two integer vectors are used to specify their “location.” The compressed-column representation serves as a good illustration. Using a dot-on-grid notation to display sparsity patterns, suppose

![](images/golub_600_649__4da72aea0776511f070bee87e5daaae4cf7420f31c93cc2d02ca14f46b8360ae.jpg)

<details>
<summary>text_image</summary>

A = 
.
</details>

The compressed-column representation stores the nonzero entries column by column in a real vector. If A is the matrix, then we denote this vector by A.val, e.g.,

$$
A. v a l = \boxed {a _ {1 1} \mid a _ {4 1} \mid a _ {5 2} \mid a _ {2 3} \mid a _ {3 3} \mid a _ {6 3} \mid a _ {1 4} \mid a _ {4 4} \mid a _ {2 5} \mid a _ {5 5} \mid a _ {6 5}}.
$$

An integer vector A.c is used to indicate where each column “begins” in A.val:

$$
A. c = \boxed { \begin{array}{c c c c c c} 1 & 3 & 4 & 7 & 9 & 1 2 \end{array} }.
$$

Thus, if $k = A . c ( j ) { : } A . c ( j + 1 ) - 1$ , then $v = A . v a l ( k )$ is the vector of nonzero components of $A ( : , j )$ . By convention, the last component of A.c houses $\boldsymbol { \mathsf { n n z } } ( A ) + 1$ where

$$
\mathfrak {n n z} (A) = \text {   the   number   of   nonzeros   in   } A.
$$

The row indices for the nonzero components in $A ( : , 1 ) , \ldots , A ( : , n )$ are encoded in an integer vector A.r, e.g.,

$$
A. r = \begin{array}{c c c c c c c c c c c c c c c c} \hline 1 & 4 & 5 & 2 & 3 & 6 & 1 & 4 & 2 & 5 & 6 \\ \hline \end{array} .
$$

In general, if $k = A . c ( j ) { : } A . c ( j + 1 ) - 1$ , then $A . v a l ( k ) = A ( A . r ( k ) , j )$ .

Note that the amount of storage required for A.r is comparable to the amount of storage required for the floating-point vector A.val. Index vectors represent one of the overheads that distinguish sparse from conventional dense matrix computations.
