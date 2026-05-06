# 11.1.8 Sparse QR and the Sparse Least Squares Problem

Suppose we want to minimize $\parallel A x - b \parallel _ { 2 }$ where $A \in \mathbb { R } ^ { m \times n }$ has full column rank and is sparse. If we are willing and able to form $A ^ { T } A$ , then we can apply sparse Cholesky technology to the normal equations $A ^ { T } A x = A ^ { T } b$ . In particular, we would compute a permutation $P$ so that $P ( A ^ { \bar { T } } A ) P ^ { T }$ has a sufficiently sparse Cholesky factor. However, aside from the pitfalls of normal equations, the matrix $A ^ { T } A$ can be dense even though

A is sparse. (Consider the case when A has a dense row.)

If we prefer to take the QR approach, then it still makes sense to reorder the columns of A, for if $A P ^ { T } = Q R$ is the thin QR factorization of $A P ^ { T }$ , then

$$
P (A ^ {T} A) P ^ {T} = R ^ {T} R,
$$

i.e., $R ^ { T }$ is the Cholesky factor of $P ( A ^ { T } A ) P ^ { T }$ . However, this poses serious issues that revolve around fill-in and the Q matrix. Suppose Q is determined via Householder QR. Even though P is chosen so that the final matrix R is reasonably sparse, the intermediate Householder updates $A = H _ { k } A$ tend to have high levels of fill-in. A corollary of this is that Q is almost always dense. This can be a show-stopper especially if $m \gg n$ and motivates the Sparse QR challenge:

# The Sparse QR Challenge

Given a sparse matrix $A \in \mathbb { R } ^ { m \times n }$ , efficiently determine a permutation p of 1:n so that if $P = I _ { n } ( : , p )$ , then the R-factor in the thin QR factorization $A ( : , p ) = A P ^ { T } = Q R$ is close to being optimally sparse. Use orthogonal transformations to determine R from $A ( : , p )$ .

Before we show how to address the challenge we establish its relevance to the sparse least squares problem. If $A P ^ { T } = Q R$ is the thin QR factorization of $A ( : , p )$ , then the normal equation system $A ^ { T } b = A ^ { T } A x _ { L S }$ transforms to

$$
P (A ^ {T} b) = (P (A ^ {T} A) P ^ {T}) P x _ {L S} = R ^ {T} R P x _ {L S}.
$$

Solving the normal equations with a QR-produced Cholesky factor constitutes the seminormal equations approach to least squares. Observe that it is not necessary to compute Q. If followed by a single step of iterative improvement, then it is possible to show that the computed $x _ { L S }$ is just as good as the least squares solution obtained via the QR factorization. Here is the overall solution framework:

Step 1. Determine P so that the Cholesky factor for $P ( A ^ { T } A ) P ^ { T }$ is sparse.

Step 2. Carefully compute the matrix R in the thin QR factorization $A P ^ { T } = Q R .$

Step 3. Solve: $R ^ { T } y _ { 0 } = P ( A ^ { T } b ) , R z _ { 0 } = y _ { 0 } , x _ { 0 } = P ^ { T } z _ { 0 }$ .

Step 4. Improve: $\boldsymbol { r } = \boldsymbol { b } - A x _ { 0 } , R ^ { T } y _ { 1 } = P ( A ^ { T } r ) , R z _ { 1 } = y _ { 1 } , \boldsymbol { e } = P ^ { T } z _ { 1 } , x _ { \iota s } = x _ { 0 } + e .$

To appreciate Steps 3 and 4, think of $x _ { 0 }$ as being contaminated by unacceptable levels of error due to the pitfalls of normal equations. Noting that $A ^ { T } \dot { A x _ { 0 } } = A ^ { T } \dot { b } - A ^ { T } r$ and $A ^ { T } A e = A ^ { T } r$ , we have

$$
A ^ {T} A (x _ {0} + e) = A ^ {T} b - A ^ {T} r + A ^ {T} r = A ^ {T} b.
$$

For a detailed analysis of the seminormal equation approach, see Bj¨orck (1987).

Let us return to the Sparse QR challenge and the efficient computaton of R using orthogonal transformations. Recall from §5.2.5 that with the Givens rotation approach there is considerable flexibility with respect to the zeroing order. A strategy for introducing zeros into $A \in \mathbb { R } ^ { m \times n }$ one row at a time can be organized as follows:

ajn = ain

The index i names the row that is being “rotated into” the current R matrix. Here is an example that shows how the j-loop oversees that process if $i > n \colon$

![](images/golub_600_649__0c0b1383a28f3131e7e011de1a9a7b0f863a51be911404284e50ca8545e06bdc.jpg)

Notice that the rotations can induce fill-in both in R and in the row that is currently being zeroed. Various row-ordering strategies have been proposed to minimize fill-in “along the way” to the final matrix R. See George and Heath (1980) and Bj¨orck (NMLS, p. 244). For example, before (11.1.9) is executed, the rows can be arranged so that the first nonzero in each row is never to the left of the first nonzero in the previous row. Rows where the first nonzero element occurs in the same column can be sorted according to the location of the last nonzero element.
