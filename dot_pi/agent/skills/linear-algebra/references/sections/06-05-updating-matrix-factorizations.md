# 6.5 Updating Matrix Factorizations

In many applications it is necessary to refactor a given matrix $A \in \mathbb { R } ^ { m \times n }$ after it has undergone a small modification. For example, given that we have the QR factorization of a matrix A, we may require the QR factorization of the matrix A obtained from A by appending a row or column or deleting a row or column. In this section we show that in situations like these, it is much more efficient to “update” A’s QR factorization than to generate the required QR factorization of A from scratch. Givens rotations have a prominent role to play. In addition to discussing various update-QR strategies, we show how to downdate a Cholesky factorization using hyperbolic rotations and how to update a rank-revealing ULV decomposition.

# 6.5.1 Rank-1 Changes

Suppose we have the QR factorization $Q R = A \in \mathbb { R } ^ { n \times n }$ and that we need to compute the QR factorization $\widetilde { A } = A + u v ^ { T } = Q _ { 1 } R _ { 1 }$ where u, $v \in \mathbb { R } ^ { n }$ are given. Observe that

$$
\widetilde {A} = A + u v ^ {T} = Q (R + w v ^ {T}) \tag {6.5.1}
$$

where $w = Q ^ { T } u$ . Suppose rotations $J _ { n - 1 } , \ldots , J _ { 2 } , J _ { 1 }$ are computed such that

$$
J _ {1} ^ {T} \dots J _ {n - 1} ^ {T} w = \pm \| w \| _ {2} e _ {1}.
$$

where each $J _ { k }$ is a Givens rotation in planes k and k + 1. If these same rotations are applied to R, then

$$
H = J _ {1} ^ {T} \dots J _ {n - 1} ^ {T} R \tag {6.5.2}
$$

is upper Hessenberg. For example, in the $n = 4$ case we start with

$$
w \leftarrow \left[ \begin{array}{c} \times \\ \times \\ \times \\ \times \end{array} \right], \qquad R \leftarrow \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right],
$$

and then update as follows:

$$
\begin{array}{l} w \leftarrow J _ {3} ^ {T} w = \left[ \begin{array}{l} \times \\ \times \\ \times \\ 0 \end{array} \right], \qquad R \leftarrow J _ {3} ^ {T} R = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \end{array} \right], \\ w \leftarrow J _ {2} ^ {T} w = \left[ \begin{array}{c} \times \\ \times \\ 0 \\ 0 \end{array} \right], \qquad R \leftarrow J _ {2} ^ {T} R = \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \end{array} \right], \\ w \leftarrow J _ {1} ^ {T} w = \left[ \begin{array}{c} \times \\ 0 \\ 0 \\ 0 \end{array} \right], \qquad H \leftarrow J _ {1} ^ {T} R = \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \end{array} \right]. \\ \end{array}
$$

Consequently,

$$
(J _ {1} ^ {T} \dots J _ {n - 1} ^ {T}) (R + w v ^ {T}) = H \pm \| w \| _ {2} e _ {1} v ^ {T} = H _ {1} \tag {6.5.3}
$$

is also upper Hessenberg. Following Algorithm 5.2.4, we compute Givens rotations $G _ { k }$ , $k = 1 { : } n - 1$ such that $\check { G _ { n - 1 } } \cdot \cdot \cdot G _ { 1 } ^ { T } \check { H } _ { 1 } = R _ { 1 }$ is upper triangular. Combining everything we obtain the QR factorization $\widetilde A = A + u v ^ { T } \ = \ Q _ { 1 } R _ { 1 }$ where

$$
Q _ {1} = Q J _ {n - 1} \dots J _ {1} G _ {1} \dots G _ {n - 1}.
$$

A careful assessment of the work reveals that about $2 6 n ^ { 2 }$ flops are required.

The technique readily extends to the case when A is rectangular. It can also be generalized to compute the QR factorization of $A + U V ^ { T }$ where $U \in \mathbb { R } ^ { m \times p }$ and $V \in \mathbb { R } ^ { n \times p }$ .

# 6.5.2 Appending or Deleting a Column

Assume that we have the QR factorization

$$
Q R = A = \left[ a _ {1} \mid \dots \mid a _ {n} \right], \quad a _ {i} \in \mathbb {R} ^ {m}, \tag {6.5.4}
$$

and for some k, $1 \leq k \leq n$ , partition the upper triangular matrix $R \in \mathbb { R } ^ { m \times n }$ as follows:

$$
R = \left[ \begin{array}{c c c} R _ {1 1} & v & R _ {1 3} \\ 0 & r _ {k k} & w ^ {T} \\ 0 & 0 & R _ {3 3} \end{array} \right] \begin{array}{c} k - 1 \\ 1 \\ m - k \end{array} .
$$

Now suppose that we want to compute the QR factorization of

$$
\widetilde {A} = \left[ a _ {1} \mid \dots \mid a _ {k - 1} \mid a _ {k + 1} \mid \dots \mid a _ {n} \right] \in \mathbb {R} ^ {m \times (n - 1)}.
$$

Note that $\widetilde { A }$ is just A with its kth column deleted and that

$$
Q ^ {T} \widetilde {A} = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 3} \\ 0 & w ^ {T} \\ 0 & R _ {3 3} \end{array} \right] = H
$$

is upper Hessenberg, e.g.,

$$
H = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 \end{array} \right], \qquad m = 7, n = 6, k = 3.
$$

Clearly, the unwanted subdiagonal elements $h _ { k + 1 , k } , \ldots , h _ { n , n - 1 }$ can be zeroed by a sequence of Givens rotations: $G _ { n - 1 } ^ { T } \cdot \cdot \cdot G _ { k } ^ { T } H ~ = ~ R _ { 1 }$ . Here, $G _ { i }$ is a rotation in planes i and $i + 1$ for $i = k { : } n - 1$ . Thus, if $Q _ { 1 } = Q G _ { k } \cdot \cdot \cdot G _ { n - 1 }$ then $\widetilde { A } = Q _ { 1 } R _ { 1 }$ is the QR factorization of ${ \widetilde { A } } .$ .

The above update procedure can be executed in $O ( n ^ { 2 } )$ flops and is very useful in certain least squares problems. For example, one may wish to examine the significance of the kth factor in the underlying model by deleting the kth column of the corresponding data matrix and solving the resulting LS problem.

Analogously, it is possible to update efficiently the QR factorization of a matrix after a column has been added. Assume that we have (6.5.4) but now want the QR factorization of

$$
\widetilde {A} = \left[ a _ {1} \mid \ldots \mid a _ {k} \mid z \mid a _ {k + 1} \mid \ldots \mid a _ {n} \right]
$$

where $z \in \mathbb { R } ^ { m }$ is given. Note that if $w = Q ^ { T } z$ then

$$
Q ^ {T} \widetilde {A} = \left[ Q ^ {T} a _ {1} \mid \dots \mid Q ^ {T} a _ {k} \mid w \mid Q ^ {T} a _ {k + 1} \mid \dots \mid Q ^ {T} a _ {n} \right]
$$

is upper triangular except for the presence of a “spike” in its $( k + 1 )$ st column, e.g.,

$$
\widetilde {A} \leftarrow Q ^ {T} \widetilde {A} = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & 0 & \times \\ 0 & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & \times & 0 & 0 \end{array} \right], \qquad m = 7, n = 5, k = 3.
$$

It is possible to determine a sequence of Givens rotations that restores the triangular form:

$$
\widetilde {A} \leftarrow J _ {6} ^ {T} \widetilde {A} = \left[ \begin{array}{l l l l l l} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & 0 & \times \\ 0 & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 \end{array} \right], \qquad \widetilde {A} \leftarrow J _ {5} ^ {T} \widetilde {A} = \left[ \begin{array}{l l l l l l} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & 0 & \times \\ 0 & 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 & 0 \end{array} \right],
$$

$$
\widetilde {A} \leftarrow J _ {4} ^ {T} \widetilde {A} = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 & 0 \end{array} \right].
$$

This update requires O(mn) flops.

# 6.5.3 Appending or Deleting a Row

Suppose we have the QR factorization $Q R = A \in \mathbb { R } ^ { m \times n }$ and now wish to obtain the QR factorization of

$$
\widetilde {A} = \left[ \begin{array}{c} w ^ {T} \\ A \end{array} \right]
$$

where $w \in \mathbb { R } ^ { n }$ . Note that

$$
\mathrm{diag} (1, Q ^ {T}) \widetilde {A} = \left[ \begin{array}{c} w ^ {T} \\ R \end{array} \right] = H
$$

is upper Hessenberg. Thus, rotations $J _ { 1 } , \ldots , J _ { n }$ can be determined so $J _ { n } ^ { T } \cdot \cdot \cdot J _ { 1 } ^ { T } H =$ $R _ { 1 }$ is upper triangular. It follows that $\begin{array} { r c l } { \bar { \cal A } } & { = } & { { \cal Q } _ { 1 } { \cal R } _ { 1 } } \end{array}$ is the desired QR factorization, where $Q _ { 1 } = \operatorname { d i a g } ( 1 , Q ) J _ { 1 } \cdot \cdot \cdot J _ { n }$ . See Algorithm 5.2.5.

No essential complications result if the new row is added between rows k and $k + 1$ of A. Indeed, if

$$
\left[ \begin{array}{l} A _ {1} \\ A _ {2} \end{array} \right] = Q R, \qquad A _ {1} \in \mathbb {R} ^ {k \times n},   A _ {2} \in \mathbb {R} ^ {(m - k) \times n},
$$

and

$$
P = \left[ \begin{array}{c c c} 0 & 1 & 0 \\ I _ {k} & 0 & 0 \\ 0 & 0 & I _ {m - k} \end{array} \right],
$$

then

$$
\operatorname{diag} (1, Q ^ {T}) P \left[ \begin{array}{c} A _ {1} \\ w ^ {T} \\ A _ {2} \end{array} \right] = \left[ \begin{array}{c} w ^ {T} \\ R \end{array} \right] = H
$$

is upper Hessenberg and we proceed as before.

Lastly, we consider how to update the QR factorization $Q R = A \in \mathbb { R } ^ { m \times n }$ when the first row of A is deleted. In particular, we wish to compute the QR factorization of the submatrix $A _ { 1 }$ in

$$
A = \left[ \begin{array}{l} z ^ {T} \\ A _ {1} \end{array} \right] _ {m - 1} ^ {1}.
$$

(The procedure is similar when an arbitrary row is deleted.) Let $q ^ { T }$ be the first row of Q and compute Givens rotations $G _ { 1 } , \ldots , G _ { m - 1 }$ such that $G _ { 1 } ^ { T } \cdot \cdot \cdot G _ { m - 1 } ^ { T } q = \alpha e _ { 1 }$ where $\alpha = \pm 1$ . Note that

$$
H = G _ {1} ^ {T} \dots G _ {m - 1} ^ {T} R = \left[ \begin{array}{l} v ^ {T} \\ R _ {1} \end{array} \right] _ {m - 1} ^ {1}
$$

is upper Hessenberg and that

$$
Q G _ {m - 1} \dots G _ {1} = \left[ \begin{array}{l l} \alpha & 0 \\ 0 & Q _ {1} \end{array} \right]
$$

where $Q _ { 1 } \in \mathbb { R } ^ { ( m - 1 ) \times ( m - 1 ) }$ is orthogonal. Thus,

$$
A = \left[ \begin{array}{l} z ^ {T} \\ A _ {1} \end{array} \right] = (Q G _ {m - 1} \dots G _ {1}) (G _ {1} ^ {T} \dots G _ {m - 1} ^ {T} R) = \left[ \begin{array}{l l} \alpha & 0 \\ 0 & Q _ {1} \end{array} \right] \left[ \begin{array}{l} v ^ {T} \\ R _ {1} \end{array} \right]
$$

from which we conclude that $A _ { 1 } = Q _ { 1 } R _ { 1 }$ is the desired QR factorization.

# 6.5.4 Cholesky Updating and Downdating

Suppose we are given a symmtetric positive definite matrix $A \in \mathbb { R } ^ { n \times n }$ and its Cholesky factor G. In the Cholesky updating problem, the challenge is to compute the Cholesky factorization $\widetilde { A } = \widetilde { G } \widetilde { G } ^ { T }$ where

$$
\widetilde {A} = A + z z ^ {T}, \quad z \in \mathbb {R} ^ {n}. \tag {6.5.5}
$$

Noting that

$$
\widetilde {A} = \left[ \begin{array}{c} G ^ {T} \\ z ^ {T} \end{array} \right] ^ {T} \left[ \begin{array}{c} G ^ {T} \\ z ^ {T} \end{array} \right], \tag {6.5.6}
$$

we can solve this problem by computing a product of Givens rotations $Q = Q _ { 1 } \cdot \cdot \cdot Q _ { n }$ so that

$$
Q ^ {T} \left[ \begin{array}{l} G ^ {T} \\ z ^ {T} \end{array} \right] = \left[ \begin{array}{l} R \\ 0 \end{array} \right], \quad R \in \mathbb {R} ^ {n \times n} \tag {6.5.7}
$$

is upper triangular. It follows that $\widetilde { A } = R R ^ { T }$ and so the updated Cholesky factor is given by $\widetilde { G } = \overline { { R } } ^ { T }$ . The zeroing sequence that produces R is straight forward, e.g.,

$$
\left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & 0 & \times \\ \times & \times & \times \end{array} \right] \xrightarrow {Q _ {1}} \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & 0 & \times \\ 0 & \times & \times \end{array} \right] \xrightarrow {Q _ {2}} \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & 0 & \times \\ 0 & 0 & \times \end{array} \right] \xrightarrow {Q _ {3}} \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & 0 & \times \\ 0 & 0 & 0 \end{array} \right].
$$

The $Q _ { k }$ update involves only rows k and $n + 1$ . The overall process is essentially the same as the strategy we outlined in the previous subsection for updating the QR factorization of a matrix when a row is appended.

The Cholesky downdating problem involves a different set of tools and a new set of numerical concerns. We are again given a Cholesky factorization $A = G G ^ { T }$ and a vector $z \in \mathbb { R } ^ { n }$ . However, now the challenge is to compute the Cholesky factorization $\widetilde { A } = \widetilde { G } \widetilde { G } ^ { T }$ where

$$
\widetilde {A} = A - z z ^ {T} \tag {6.5.8}
$$

is presumed to be positive definite. By introducing the notion of a hyperbolic rotation we can develop a downdating framework that corresponds to the Givens-based updating framework. Define the matrix S as follows

$$
S = \left[ \begin{array}{c c} I _ {n} & 0 \\ 0 & - 1 \end{array} \right] \tag {6.5.9}
$$

and note that

$$
\widetilde {A} = G G ^ {T} - z z ^ {T} = \left[ \begin{array}{c} G ^ {T} \\ z ^ {T} \end{array} \right] ^ {T} S \left[ \begin{array}{c} G ^ {T} \\ z ^ {T} \end{array} \right]. \tag {6.5.10}
$$

This corresponds to (6.5.6), but instead of computing the QR factorization (6.5.7), we seek a matrix $H \in \mathbb { R } ^ { ( n + 1 ) \times ( n + 1 ) }$ that satisfies two properties:

$$
H S H ^ {T} = S, \tag {6.5.11}
$$

$$
H ^ {T} \left[ \begin{array}{c} G ^ {T} \\ z ^ {T} \end{array} \right] = \left[ \begin{array}{c} R \\ 0 \end{array} \right], \quad R \in \mathbb {R} ^ {n \times n} (\text {upper triangular}). \tag {6.5.12}
$$

If this can be accomplished, then it follows from

$$
\widetilde {A} = \left(H ^ {T} \left[ \begin{array}{l} G ^ {T} \\ z ^ {T} \end{array} \right]\right) ^ {T} \left[ \begin{array}{l l} I _ {n} & 0 \\ 0 & - 1 \end{array} \right] \left(H ^ {T} \left[ \begin{array}{l} G ^ {T} \\ z ^ {T} \end{array} \right]\right) = R ^ {T} R
$$

that the Cholesky factor of $\widetilde { A } = A - z z ^ { T }$ is given by $\widetilde { G } = R ^ { T }$ . A matrix H that satisfies (6.5.11) is said to be S-orthogonal. Note that the product of S-orthogonal matrices is also S-orthogonal.

An important subset of the S-orthogonal matrices are the hyperbolic rotations and here is a 4-by-4 example:

$$
H _ {2} (\theta) = \left[ \begin{array}{c c c c} 1 & 0 & 0 & 0 \\ 0 & c & 0 & - s \\ 0 & 0 & 1 & 0 \\ 0 & - s & 0 & c \end{array} \right], \qquad c = \cosh (\theta), s = \sinh (\theta).
$$

The S-orthogonality of this matrix follows from cosh $( \theta ) ^ { 2 } - \sinh ( \theta ) ^ { 2 } = 1$ . In general, $H _ { k } \in \mathbb { R } ^ { ( n + 1 ) \times ( n + 1 ) }$ is a hyperbolic rotation if it agrees with $I _ { n + 1 }$ except in four locations:

$$
\left[ \begin{array}{c c} [ H _ {k} ] _ {k, k} & [ H _ {k} ] _ {k, n + 1} \\ [ H _ {k} ] _ {n + 1, k} & [ H _ {k} ] _ {n + 1, n + 1} \end{array} \right] = \left[ \begin{array}{c c} \cosh (\theta) & - \sinh (\theta) \\ - \sinh (\theta) & \cosh (\theta) \end{array} \right].
$$

Hyperbolic rotations look like Givens rotations and, not surprisingly, can be used to introduce zeros into a vector or matrix. However, upon consideration of the equation

$$
\left[ \begin{array}{c c} {c} & {- s} \\ {- s} & {c} \end{array} \right] \left[ \begin{array}{l} {x _ {1}} \\ {x _ {2}} \end{array} \right] = \left[ \begin{array}{l} {r} \\ {0} \end{array} \right], \qquad c ^ {2} - s ^ {2} = 1
$$

we see that the required cosh-sinh pair may not exist. Since we always have $| \cosh ( \theta ) | >$ $| \sinh ( \theta ) |$ , there is no real solution t $\mathrm { ~ o ~ } - s x _ { 1 } + c x _ { 2 } = 0 \mathrm { ~ i f ~ } | x _ { 2 } | > | x _ { 1 } |$ . On the other hand, $\mathrm { i f } \ | x _ { 1 } | > | x _ { 2 } |$ , then $\{ c , s \} = \{ \cosh ( \theta ) , \sinh ( \theta ) \}$ can be computed as follows:

$$
\tau = \frac {x _ {2}}{x _ {1}}, \quad c = \frac {1}{\sqrt {1 - \tau^ {2}}}, \quad s = c \cdot \tau . \tag {6.5.13}
$$

There are clearly numerical issues $\mathrm { i f } \ | x _ { 1 } |$ is just slightly greater than $| x _ { 2 } |$ . However, it is possible to organize hyperbolic rotation computations successfully, see Alexander, Pan, and Plemmons (1988).

Putting these concerns aside, we show how the matrix H in (6.5.12) can be computed as a product of hyperbolic rotations $H = H _ { 1 } \cdot \cdot \cdot H _ { n }$ just as the transforming Q in the updating problem is a product of Givens rotations. Consider the role of $H _ { 1 }$ in the $n = 3$ case:

$$
\left[ \begin{array}{c c c c} c & 0 & 0 & - s \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ - s & 0 & 0 & c \end{array} \right] ^ {T} \left[ \begin{array}{c c c} g _ {1 1} & g _ {2 1} & g _ {3 1} \\ 0 & g _ {2 2} & g _ {3 2} \\ 0 & 0 & g _ {3 3} \\ z _ {1} & z _ {2} & z _ {3} \end{array} \right] = \left[ \begin{array}{c c c} \tilde {g} _ {1 1} & \tilde {g} _ {2 1} & \tilde {g} _ {3 1} \\ 0 & g _ {2 2} & g _ {3 2} \\ 0 & 0 & g _ {3 3} \\ 0 & z _ {2} ^ {\prime} & z _ {3} ^ {\prime} \end{array} \right].
$$

Since $\widetilde { A } = G G ^ { T } - z z ^ { T }$ is positive definite, $[ \widetilde { A } ] _ { 1 1 } = g _ { 1 1 } ^ { 2 } - z _ { 1 } ^ { 2 } > 0$ . It follows that $\vert g _ { 1 1 } \vert > \vert z _ { 1 } \vert$ which guarantees that the cosh-sinh computations (6.5.13) go through. For the overall process to be defined, we have to guarantee that hyperbolic rotations $H _ { 2 } , \ldots , H _ { n }$ can be found to zero out the bottom row in the matrix $\big [ G ^ { T } ~ z \big ] ^ { T }$ . The following theorem ensures that this is the case.

Theorem 6.5.1. If

$$
A = \left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right] = \left[ \begin{array}{c c} g _ {1 1} & 0 \\ g _ {1} & G _ {1} \end{array} \right] \left[ \begin{array}{c c} g _ {1 1} & g _ {1} ^ {T} \\ 0 & G _ {1} ^ {T} \end{array} \right]
$$

and

$$
\widetilde {A} = A - z z ^ {T} = A - \left[ \begin{array}{l} \mu \\ w \end{array} \right] \left[ \begin{array}{l} \mu \\ w \end{array} \right] ^ {T}
$$

are positive definite, then it is possible to determine $c = \cosh ( \theta )$ and $s = \sinh ( \theta )$ so

$$
\left[ \begin{array}{c c c} c & 0 & - s \\ 0 & I _ {n - 1} & 0 \\ - s & 0 & c \end{array} \right] \left[ \begin{array}{c c} g _ {1 1} & g _ {1} ^ {T} \\ 0 & G _ {1} ^ {T} \\ \mu & w ^ {T} \end{array} \right] = \left[ \begin{array}{c c} \tilde {g} _ {1 1} & \tilde {g} _ {1} ^ {T} \\ 0 & G _ {1} ^ {T} \\ 0 & w _ {1} ^ {T} \end{array} \right].
$$

Moreover, the matrix $\widetilde { A } _ { 1 } = G _ { 1 } G _ { 1 } ^ { T } - w _ { 1 } w _ { 1 } ^ { T }$ is positive definite.

Proof. The blocks in A’s Cholesky factor are given by

$$
g _ {1 1} = \sqrt {\alpha}, \quad g _ {1} = v / g _ {1 1}, \quad G _ {1} G _ {1} ^ {T} = B - \frac {1}{\alpha} v v ^ {T}. \tag {6.5.14}
$$

Since $A - z z ^ { T }$ is positive definite, $a _ { 1 1 } - z _ { 1 } ^ { 2 } = g _ { 1 1 } ^ { 2 } - \mu ^ { 2 } > 0$ and so from (6.5.13) with $\tau = \mu / g _ { 1 1 }$ we see that

$$
c = \frac {\sqrt {\alpha}}{\sqrt {\alpha - \mu^ {2}}}, \quad s = \frac {\mu}{\sqrt {\alpha - \mu^ {2}}}. \tag {6.5.15}
$$

Since $w _ { 1 } = - s g _ { 1 } + c w$ it follows from (6.5.14) and (6.5.15) that

$$
\begin{array}{l} \widetilde {A} _ {1} = G _ {1} G _ {1} ^ {T} - w _ {1} w _ {1} ^ {T} = B - \frac {1}{\alpha} v v ^ {T} - (- s g _ {1} + c w) (- s g _ {1} + c w) ^ {T} \\ = B - \frac {c ^ {2}}{\alpha} v v ^ {T} - c ^ {2} w w ^ {T} + \frac {s c}{\sqrt {\alpha}} (v w ^ {T} + w v ^ {T}) \\ = B - \frac {1}{\alpha - \mu^ {2}} v v ^ {T} - \frac {\alpha}{\alpha - \mu^ {2}} w w ^ {T} + \frac {\mu}{\alpha - \mu^ {2}} (v w ^ {T} + w v ^ {T}). \\ \end{array}
$$

It is easy to verify that this matrix is precisely the Schur complement of α in

$$
\widetilde {A} = A - z z ^ {T} = \left[ \begin{array}{l l} \alpha - \mu^ {2} & v ^ {T} - \mu w ^ {T} \\ v - \mu w & B - w w ^ {T} \end{array} \right]
$$

and is therefore positive definite.

The theorem provides the key step in an induction proof that the factorization (6.5.12) exists.

# 6.5.5 Updating a Rank-Revealing ULV Decomposition

We close with a discussion about updating a nullspace basis after one or more rows have been appended to the underlying matrix. We work with the ULV decomposition which is much more tractable than the SVD from the updating point of view. We pattern our remarks after Stewart(1993).

A rank -revealing ULV decomposition of a matrix $A \in \mathbb { R } ^ { m \times n }$ has the form

$$
U ^ {T} A V = \left[ \begin{array}{l} L \\ 0 \end{array} \right] = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \\ 0 & 0 \end{array} \right], \quad U ^ {T} U = I _ {m}, V ^ {T} V = I _ {n} \tag {6.5.16}
$$

where $L _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ and $L _ { 2 2 } \in \mathbb { R } ^ { ( n - r ) \times ( n - r ) }$ are lower triangular and $\parallel L _ { 2 1 } \parallel _ { 2 }$ and $\parallel L _ { 2 2 } \parallel _ { 2 }$ are small compared to $\sigma _ { \mathrm { m i n } } ( L _ { 1 1 } )$ . Such a decomposition can be obtained by applying QR with column pivoting

$$
U ^ {T} A \Pi = \left[ \begin{array}{c} R \\ 0 \end{array} \right], \qquad R \in \mathbb {R} ^ {n \times n}
$$

followed by a QR factorization $V _ { 1 } ^ { T } R ^ { T } = L ^ { T }$ . In this case the matrix V in (6.5.16) is given by $V = \Pi V _ { 1 }$ . The parameter r is the estimated rank. Note that if

$$
V = \left[ \begin{array}{c c} V _ {1} & V _ {2} \\ r & n - r \end{array} \right], \qquad U = \left[ \begin{array}{c c} U _ {1} & U _ {2} \\ r & m - r \end{array} \right],
$$

then the columns of $V _ { 2 }$ define an approximate nullspace:

$$
\left\| A V _ {2} \right\| _ {2} = \left\| U _ {2} L _ {2 2} \right\| _ {2} = \left\| L _ {2 2} \right\| _ {2}.
$$

Our goal is to produce cheaply a rank-revealing ULV decomposition for the rowappended matrix

$$
\tilde {A} = \left[ \begin{array}{c} A \\ z ^ {T} \end{array} \right],
$$

In particular, we show how to revise L, V , and possibly r in $O ( n ^ { 2 } )$ flops. Note that

$$
\left[ \begin{array}{l l} U & 0 \\ 0 & 1 \end{array} \right] ^ {T} \left[ \begin{array}{l} A \\ z ^ {T} \end{array} \right] V = \left[ \begin{array}{l l} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \\ 0 & 0 \\ w ^ {T} & y ^ {T} \end{array} \right].
$$

We illustrate the key ideas through an example. Suppose $n = 7$ and $r \ = \ 4$ . By permuting the rows so that the bottom row is just underneath L, we obtain

$$
\left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \\ w ^ {T} & y ^ {T} \end{array} \right] = \left[ \begin{array}{c c c c c c c} \ell & 0 & 0 & 0 & 0 & 0 & 0 \\ \ell & \ell & 0 & 0 & 0 & 0 & 0 \\ \ell & \ell & \ell & 0 & 0 & 0 & 0 \\ \ell & \ell & \ell & \ell & 0 & 0 & 0 \\ \hline \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & 0 & 0 \\ \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & 0 \\ \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon \\ \hline w & w & w & w & y & y & y \end{array} \right].
$$

The 
 entries are small while the , w, and y entries are not. Next, a sequence of Givens rotations $G _ { 7 } , \ldots , G _ { 1 }$ are applied from the left to zero out the bottom row:

$$
\left[ \frac {\tilde {L}}{0} \right] = \left[ \begin{array}{l l l l l l l} \times & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 & 0 \\ \times & \times & \times & \times & 0 & 0 & 0 \\ \times & \times & \times & \times & \times & 0 & 0 \\ \times & \times & \times & \times & \times & \times & 0 \\ \times & \times & \times & \times & \times & \times & \times \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{array} \right] = G _ {1 7} \dots G _ {5 7} G _ {6 7} \left[ \begin{array}{l l} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \\ w ^ {T} & y ^ {T} \end{array} \right].
$$

Because this zeroing process intermingles the (presumably large) entries of the bottom row with the entries from each of the other rows, the lower triangular form is typically not rank revealing. However, and this is key, we can restore the rank-revealing structure with a combination of condition estimation and Givens zero chasing.

Let us assume that with the added row, the new nullspace has dimension 2. With a reliable condition estimator we produce a unit 2-norm vector p such that

$$
\parallel p ^ {T} \widetilde {L} \parallel_ {2} \approx \sigma_ {\mathrm{min}} (\widetilde {L}).
$$

(See §3.5.4). Rotations $\{ U _ { i , i + 1 } \} _ { i = 1 } ^ { 6 }$ can be found such that

$$
U _ {6 7} ^ {T} U _ {5 6} ^ {T} U _ {4 5} ^ {T} U _ {3 4} ^ {T} U _ {2 3} ^ {T} U _ {1 2} ^ {T} p = e _ {7} = I _ {7} (:, 7).
$$

Applying these rotations to $\widetilde { L }$ produces a lower Hessenberg matrix

$$
H = U _ {6 7} ^ {T} U _ {5 6} ^ {T} U _ {4 5} ^ {T} U _ {3 4} ^ {T} U _ {2 3} ^ {T} U _ {1 2} ^ {T} \tilde {L}.
$$

Applying more rotations from the right restores H to a lower triangular form:

$$
L _ {+} = H V _ {1 2} V _ {2 3} V _ {3 4} V _ {4 5} V _ {5 6} V _ {6 7}.
$$

It follows that

$$
e _ {7} ^ {T} L _ {+} = \left(e _ {8} ^ {T} H\right) V _ {1 2} V _ {2 3} V _ {3 4} V _ {4 5} V _ {5 6} V _ {6 7} = \left(p ^ {T} \tilde {L}\right) V _ {1 2} V _ {2 3} V _ {3 4} V _ {4 5} V _ {5 6} V _ {6 7}
$$

has approximate norm $\sigma _ { \mathrm { m i n } } ( \widetilde { L } )$ . Thus, we obtain a lower triangular matrix of the form

$$
L _ {+} = \left[ \begin{array}{c c c c c c c} \times & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 & 0 \\ \times & \times & \times & \times & 0 & 0 & 0 \\ \times & \times & \times & \times & \times & 0 & 0 \\ \times & \times & \times & \times & \times & \times & 0 \\ \hline \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon \end{array} \right]
$$

We can repeat the condition estimation and zero chasing on the leading 6-by-6 portion. Assuming that the nullspace of the augmented matrix has dimension two, this produces another row of small numbers:

$$
\left[ \begin{array}{c c c c c c c} \times & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 & 0 \\ \times & \times & \times & \times & 0 & 0 & 0 \\ \times & \times & \times & \times & \times & 0 & 0 \\ \hline \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & 0 \\ \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon \end{array} \right].
$$

This illustrates how we can restore any lower triangular matrix to rank-revealing form.

# Problems

P6.5.1 Suppose we have the QR factorization for $A \in \mathbb { R } ^ { m \times n }$ and now wish to solve

$$
\min _ {x \in \mathbf {R} ^ {n}} \| (A + u v ^ {T}) x - b \| _ {2}
$$

where $u , b \in \mathbb { R } ^ { m }$ and $v \in \mathbb { R } ^ { n }$ are given. Give an algorithm for solving this problem that requires $O ( m n )$ flops. Assume that $Q$ must be updated.

P6.5.2 Suppose

$$
A   =   \left[ \begin{array}{c} c ^ {T} \\ B \end{array} \right], \qquad c \in \mathbb {R} ^ {n},   B \in \mathbb {R} ^ {(m - 1) \times n}
$$

has full column rank and $m > n$ . Using the Sherman-Morrison-Woodbury formula show that

$$
\frac {1}{\sigma_ {\min} (B)} \leq \frac {1}{\sigma_ {\min} (A)} + \frac {\| (A ^ {T} A) ^ {- 1} c \| _ {2} ^ {2}}{1 - c ^ {T} (A ^ {T} A) ^ {- 1} c}.
$$

P6.5.3 As a function of $x _ { 1 }$ and $x _ { 2 }$ , what is the 2-norm of the hyperbolic rotation produced by (6.5.13)?

P6.5.4 Assume that

$$
A   =   \left[ \begin{array}{c c} R & H \\ 0 & E \end{array} \right], \qquad \rho   =   \frac {\|   E   \| _ {2}}{\sigma_ {\min} (R)}   <     1,
$$

where R and E are square. Show that if

$$
Q = \left[ \begin{array}{l l} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right]
$$

is orthogonal and

$$
\left[ \begin{array}{c c} R & H \\ 0 & E \end{array} \right] \left[ \begin{array}{c c} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right] = \left[ \begin{array}{c c} R _ {1} & 0 \\ H _ {1} & E _ {1} \end{array} \right],
$$

then  $H _ { 1 } \parallel _ { 2 } \leq \rho \parallel H \parallel _ { 2 }$ .

P6.5.5 Suppose $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ with $m \geq n$ . In the indefinite least squares (ILS) problem, the goal is to minimize

$$
\phi (x) = (b - A x) ^ {T} J (b - A x),
$$

where

$$
S = \left[ \begin{array}{c c} I _ {p} & 0 \\ 0 & - I _ {q} \end{array} \right], \qquad p + q = m.
$$

It is assumed that $p \geq 1$ and $q \geq 1$ . (a) By taking the gradient of $\phi ,$ , show that the ILS problem has a unique solution if and only if $A ^ { T } S A$ is positive definite. (b) Assume that the ILS problem has a unique solution. Show how it can be found by computing the Cholesky factorization of $\bar { Q } _ { 1 } ^ { T } Q _ { 1 } - Q _ { 2 } ^ { T } Q _ { 2 }$ where

$$
A   =   \left[ \begin{array}{l} Q _ {1} \\ Q _ {2} \end{array} \right], \qquad Q _ {1} \in \mathbb {R} ^ {p \times n},   Q _ {2} \in \mathbb {R} ^ {q \times n}
$$

is the thin QR factorization. (c) A matrix $Q \in \mathbb { R } ^ { m \times m }$ is S-orthogonal if $Q S Q ^ { T } = S { \mathrm { ~ I f } } $

$$
Q = \left[ \begin{array}{c c} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right] _ {q} ^ {p}
$$

is S-orthogonal, then by comparing blocks in the equation $Q ^ { T } S Q = S$ we have

$$
Q _ {1 1} ^ {T} Q _ {1 1} = I _ {p} + Q _ {2 1} ^ {T} Q _ {2 1}, \qquad Q _ {1 1} ^ {T} Q _ {1 2} = Q _ {2 1} ^ {T} Q _ {2 2}, \qquad Q _ {2 2} ^ {T} Q _ {2 2} = I _ {q} + Q _ {1 2} ^ {T} Q _ {1 2}.
$$

Thus, the singular values of $Q _ { 1 1 }$ and $Q _ { 2 2 }$ are never smaller than 1. Assume that $p \geq q$ . By analogy with how the CS decomposition is established in §2.5.4, show that there exist orthogonal matrices $U _ { 1 }$ , $U _ { 2 } , V _ { 1 }$ and $V _ { 2 }$ such that

$$
\left[\begin{array}{c c}U _ {1}&0\\0&U _ {2}\end{array}\right] ^ {T} Q \left[\begin{array}{c c}V _ {1}&0\\0&V _ {2}\end{array}\right] = \left[ \right.\begin{array}{c c}D&0\\0&I _ {p - q}\\\hline (D ^ {2} - I _ {p}) ^ {1 / 2}&0\end{array}\left. \right|\begin{array}{c}(D ^ {2} - I) ^ {1 / 2}\\0\\D\end{array}\left. \right]
$$

where $D = \operatorname { d i a g } ( d _ { 1 } , \dotsc , d _ { p } )$ with $d _ { i } \geq 1 , i = 1 { : } p$ . This is the hyperbolic CS decomposition and details can be found in Stewart and Van Dooren (2006).

# Notes and References for 6.5

The seminal matrix factorization update paper is:

P.E. Gill, G.H. Golub, W. Murray, and M.A. Saunders (1974). “Methods for Modifying Matrix Factorizations,” Math. Comput. 28, 505–535.

Initial research into the factorization update problem was prompted by the development of quasi-Newton methods and the simplex method for linear programming. In these venues, a linear system must be solved in step k that is a low-rank perturbation of the linear system solved in step k − 1, see:

R.H. Bartels (1971). “A Stabilization of the Simplex Method,” Numer. Math. 16, 414–434.

P.E. Gill, W. Murray, and M.A. Saunders (1975). “Methods for Computing and Modifying the LDV Factors of a Matrix,” Math. Comput. 29, 1051–1077.

D. Goldfarb (1976). “Factored Variable Metric Methods for Unconstrained Optimization,” Math. Comput. 30, 796–811.

J.E. Dennis and R.B. Schnabel (1983). Numerical Methods for Unconstrained Optimization and Nonlinear Equations, Prentice-Hall, Englewood Cliffs, NJ.

W.W. Hager (1989). “Updating the Inverse of a Matrix,” SIAM Review 31, 221–239.

S.K. Eldersveld and M.A. Saunders (1992). “A Block-LU Update for Large-Scale Linear Programming,” SIAM J. Matrix Anal. Applic. 13, 191–201.

Updating issues in the least squares setting are discussed in:

J. Daniel, W.B. Gragg, L. Kaufman, and G.W. Stewart (1976). “Reorthogonaization and Stable Algorithms for Updating the Gram-Schmidt QR Factorization,” Math. Comput. 30, 772–795.

S. Qiao (1988). “Recursive Least Squares Algorithm for Linear Prediction Problems,” SIAM J. Matrix Anal. Applic. 9, 323–328.

˚A. Bj¨orck, H. Park, and L. Eld´en (1994). “Accurate Downdating of Least Squares Solutions,” SIAM J. Matrix Anal. Applic. 15, 549–568.

S.J. Olszanskyj, J.M. Lebak, and A.W. Bojanczyk (1994). “Rank-k Modification Methods for Recursive Least Squares Problems,” Numer. Alg. 7, 325–354.

L. Eld´en and H. Park (1994). “Block Downdating of Least Squares Solutions,” SIAM J. Matrix Anal. Applic. 15, 1018–1034.

Kalman filtering is a very important tool for estimating the state of a linear dynamic system in the presence of noise. An illuminating, stable implementation that involves updating the QR factorization of an evolving block banded matrix is given in:

C.C. Paige and M.A. Saunders (1977). “Least Squares Estimation of Discrete Linear Dynamic Systems Using Orthogonal Transformations,” SIAM J. Numer. Anal. 14, 180–193.

The Cholesky downdating literature includes:

G.W. Stewart (1979). “The Effects of Rounding Error on an Algorithm for Downdating a Cholesky Factorization,” J. Inst. Math. Applic. 23, 203–213.

A.W. Bojanczyk, R.P. Brent, P. Van Dooren, and F.R. de Hoog (1987). “A Note on Downdating the Cholesky Factorization,” SIAM J. Sci. Stat. Comput. 8, 210–221.

C.-T. Pan (1993). “A Perturbation Analysis of the Problem of Downdating a Cholesky Factorization,” Lin. Alg. Applic. 183, 103–115.

L. Eld´en and H. Park (1994). “Perturbation Analysis for Block Downdating of a Cholesky Decomposition,” Numer. Math. 68, 457–468.

M.R. Osborne and L. Sun (1999). “A New Approach to Symmetric Rank-One Updating,” IMA J. Numer. Anal. 19, 497–507.

E.S. Quintana-Orti and R.A. Van Geijn (2008). “Updating an LU Factorization with Pivoting,” ACM Trans. Math. Softw. 35(2), Article 11.

Hyperbolic tranformations have been successfully used in a number of settings:

G.H. Golub (1969). “Matrix Decompositions and Statistical Computation,” in Statistical Computation, ed., R.C. Milton and J.A. Nelder, Academic Press, New York, pp. 365–397.

C.M. Rader and A.O. Steinhardt (1988). “Hyperbolic Householder Transforms,” SIAM J. Matrix Anal. Applic. 9, 269–290.

S.T. Alexander, C.T. Pan, and R.J. Plemmons (1988). “Analysis of a Recursive Least Squares Hyperbolic Rotation Algorithm for Signal Processing,” Lin. Alg. and Its Applic. 98, 3–40.   
G. Cybenko and M. Berry (1990). “Hyperbolic Householder Algorithms for Factoring Structured Matrices,” SIAM J. Matrix Anal. Applic. 11, 499–520.   
A.W. Bojanczyk, R. Onn, and A.O. Steinhardt (1993). “Existence of the Hyperbolic Singular Value Decomposition,” Lin. Alg. Applic. 185, 21–30.   
S. Chandrasekaran, M. Gu, and A.H. Sayad (1998). “A Stable and Efficient Algorithm for the Indefinite Linear Least Squares Problem,” SIAM J. Matrix Anal. Applic. 20, 354–362.   
A.J. Bojanczyk, N.J. Higham, and H. Patel (2003a). “Solving the Indefinite Least Squares Problem by Hyperbolic QR Factorization,” SIAM J. Matrix Anal. Applic. 24, 914–931.   
A. Bojanczyk, N.J. Higham, and H. Patel (2003b). “The Equality Constrained Indefinite Least Squares Problem: Theory and Algorithms,” BIT 43, 505–517.   
M. Stewart and P. Van Dooren (2006). “On the Factorization of Hyperbolic and Unitary Transformations into Rotations,” SIAM J. Matrix Anal. Applic. 27, 876–890.   
N.J. Higham (2003). “J-Orthogonal Matrices: Properties and Generation,” SIAM Review 45, 504–519.   
High-performance issues associated with QR updating are discussed in:   
B.C. Gunter and R.A. Van De Geijn (2005). “Parallel Out-of-Core Computation and Updating of the QR Factorization,” ACM Trans. Math. Softw. 31, 60–78.   
Updating and downdating the ULV and URV decompositions and related topics are covered in:   
C.H. Bischof and G.M. Shroff (1992). “On Updating Signal Subspaces,” IEEE Trans. Signal Proc. 40, 96–105.   
G.W. Stewart (1992). “An Updating Algorithm for Subspace Tracking,” IEEE Trans. Signal Proc. 40, 1535–1541.   
G.W. Stewart (1993). “Updating a Rank-Revealing ULV Decomposition,” SIAM J. Matrix Anal. Applic. 14, 494–499.   
G.W. Stewart (1994). “Updating URV Decompositions in Parallel,” Parallel Comp. 20, 151–172.   
H. Park and L. Eld´en (1995). “Downdating the Rank-Revealing URV Decomposition,” SIAM J. Matrix Anal. Applic. 16, 138–155.   
J.L. Barlow and H. Erbay (2009). “Modifiable Low-Rank Approximation of a Matrix,” Num. Lin. Alg. Applic. 16, 833–860.

Other interesting update-related topics include the updating of condition estimates, see:

W.R. Ferng, G.H. Golub, and R.J. Plemmons (1991). “Adaptive Lanczos Methods for Recursive Condition Estimation,” Numerical Algorithms 1, 1-20.

G. Shroff and C.H. Bischof (1992). “Adaptive Condition Estimation for Rank-One Updates of QR Factorizations,” SIAM J. Matrix Anal. Applic. 13, 1264–1278.

D.J. Pierce and R.J. Plemmons (1992). “Fast Adaptive Condition Estimation,” SIAM J. Matrix Anal. Applic. 13, 274–291.

and the updating of solutions to constrained least squares problems:

K. Schittkowski and J. Stoer (1979). “A Factorization Method for the Solution of Constrained Linear Least Squares Problems Allowing for Subsequent Data changes,” Numer. Math. 31, 431–463.

˚A. Bj¨orck (1984). “A General Updating Algorithm for Constrained Linear Least Squares Problems,” SIAM J. Sci. Stat. Comput. 5, 394–402.

Finally, we mention the following paper concerned with SVD updating:

M. Moonen, P. Van Dooren, and J. Vandewalle (1992). “A Singular Value Decomposition Updating Algorithm,” SIAM J. Matrix Anal. Applic. 13, 1015–1038.
