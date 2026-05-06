# 11.5.8 Incomplete Cholesky Preconditioners

Assume that $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite and that we are driven to consider the PCG method because A’s Cholesky factor G has many more nonzero entries than the lower triangular portion of A. A natural idea for a preconditioner is to set $M =$ $H H ^ { T }$ where H is a sufficiently sparse lower triangular matrix so that if

$$
R = H H ^ {T} - A \tag {11.5.6}
$$

then

$$
a _ {i j} \neq 0 \Rightarrow r _ {i j} = 0. \tag {11.5.7}
$$

This means that $[ H H ^ { T } ] _ { i j } = a _ { i j }$ for all nonzero $a _ { i j }$ . In this sense, $M = H H ^ { T }$ captures the essence of A. To articulate what we mean by a “sufficiently sparse” H matrix, we specify a set P of subdiagonal index pairs and insist that

$$
(i, j) \in P \Rightarrow h _ {i j} = 0. \tag {11.5.8}
$$

Given $P _ { : }$ , any matrix H that satisfies (11.5.6)–(11.5.8) is an incomplete Cholesky factor of A.

It turns out that it is not always possible to compute H given P . To see what the issues are consider the outer-product implementation of the Cholesky factorization. Recall from §4.2 that it involves repeated application of the factorization

$$
\left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right] = \left[ \begin{array}{c c} \sqrt {\alpha} & 0 \\ w & I _ {n - 1} \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ 0 & A _ {1} \end{array} \right] \left[ \begin{array}{c c} \sqrt {\alpha} & w ^ {T} \\ 0 & I _ {n - 1} \end{array} \right] \tag {11.5.9}
$$

where $w = v / \sqrt { \alpha }$ and $A _ { 1 } = B - w w ^ { T }$ . Indeed, if $G _ { 1 }$ is the Cholesky factor of $A _ { 1 }$ , then

$$
G = \left[ \begin{array}{c c} \sqrt {\alpha} & 0 \\ w & G _ {1} \end{array} \right]
$$

is the Cholesky factor of A. Now suppose $Z \in \mathbb { R } ^ { n \times n }$ is a matrix of zeros and ones with $z _ { i j } = z _ { j i } = 0$ if and only if $( i , j ) \in P$ . To ensure the existence of an incomplete Cholesky factor with respect to P , we need to guarantee that the following recursive function works:

$\mathbf { f u n c t i o n } \ H = \mathsf { i n c C h o l } ( A , Z , n )$

$\mathbf { i f } \ n = 1$

$$
H = \sqrt {A}
$$

$$
\alpha = A (1, 1), v = A (2: n, 1), B = A (2: n, 2: n)
$$

$$
w = (v / \sqrt {\alpha}) \cdot * Z (2: n, 1)
$$

$$
A _ {1} = \left(B - w w ^ {T}\right). * Z (2: n, 2: n), H _ {1} = \operatorname{incChol} \left(A _ {1}, Z (2: n, 2: n), n - 1\right)
$$

$$
H = \left[ \begin{array}{c c} \sqrt {\alpha} & 0 \\ w & H _ {1} \end{array} \right]
$$

end

If Z is the matrix of all 1’s, then this is just a recursive form of Cholesky factorization. (Set $r = 1$ in Algorithm 4.2.4). As it stands, it is Cholesky with forced zeros in both the w and $A _ { 1 }$ calculations. It is easy to show that if the algorithm runs to completion, then Equations (11.5.6), (11.5.7), and (11.5.8) are satisfied. One way to guarantee that this happens is to show that $A _ { 1 }$ is positive definite. This turns out to be the case if A is a Stieltjes matrix. A matrix $A \in \mathbb { R } ^ { n \times n }$ is a Stieltjes matrix if it is symmetric positive definite and has nonpositive off-diagonal entries. This property holds in many applications. For example, the model problem matrices in §4.8.3 are Stieltjes matrices. Using the notation $C \geq 0$ to mean that matrix C has nonnegative entries, we show that if A is a Stieltjes, then $A ^ { - 1 } \geq 0$ .

Lemma 11.5.1. If $A \in \mathbb { R } ^ { n \times n }$ is a Stieltjes matrix, then $A ^ { - 1 } \geq 0$

Proof. Write $A = D - E$ where D and −E are the diagonal and off-diagonal parts. Since $A = D ^ { 1 / 2 } ( I - F ) D ^ { 1 / 2 }$ it follows that the spectral radius of $F = D ^ { - 1 / 2 } E \hat { D } ^ { - 1 / 2 }$ satisfies $\rho ( F ) < 1$ . Thus, the entries of

$$
A ^ {- 1} = D ^ {- 1 / 2} \left(\sum_ {k = 0} ^ {\infty} F ^ {k}\right) D ^ {- 1 / 2}
$$

are clearly nonnegative.

The following result is what we need to guarantee that the function incChol does not break down.

Theorem 11.5.2. If

$$
A = \left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right], \qquad \alpha \in \mathbb {R},   v \in \mathbb {R} ^ {n - 1},   B \in \mathbb {R} ^ {(n - 1) \times (n - 1)},
$$

is a Stieltjes matrix and $\tilde { v } \in \mathbb { R } ^ { n - 1 }$ is obtained from v by setting any subset of its components to zero, then

$$
\tilde {B} = B - \frac {\tilde {v} \tilde {v} ^ {T}}{\alpha}
$$

is a Stieltjes matrix.

Proof. It is clear that $\tilde { B } = \left( \tilde { b } _ { i j } \right)$ has nonpositive off-diagonal entries since $\tilde { v } \leq 0$ and

$$
\tilde {b} _ {i j} = b _ {i j} - \frac {\tilde {v} _ {i} \tilde {v} _ {j}}{\alpha}.
$$

Our task is to show that $\tilde { B }$ is positive definite.

Since A is positive definite it follows that if

$$
x = \frac {1}{\sqrt {\alpha}} \left[ \begin{array}{c} 1 \\ - B ^ {- 1} v \end{array} \right]
$$

then

$$
0 <   x ^ {T} A x = 1 - \frac {v ^ {T} B ^ {- 1} v}{\alpha}.
$$

Since $B ^ { - 1 } \geq 0$ and $v \leq 0$ , we have $\tilde { v } ^ { T } B ^ { - 1 } \tilde { v } \leq v ^ { T } B ^ { - 1 } v$ and so

$$
\gamma \equiv 1 - \frac {\tilde {v} ^ {T} B ^ {- 1} \tilde {v}}{\alpha} \geq 1 - \frac {v ^ {T} B ^ {- 1} v}{\alpha} > 0.
$$

Using the Sherman-Morrison formula

$$
\tilde {B} ^ {- 1} = \left(B - \frac {\tilde {v} \tilde {v} ^ {T}}{\alpha}\right) ^ {- 1} = B ^ {- 1} + \frac {1}{\gamma} B ^ {- 1} \frac {\tilde {v} \tilde {v} ^ {T}}{\alpha} B ^ {- 1}
$$

we see that $\tilde { B }$ is positive definite.

A theorem of this variety can be found in the landmark paper by Meijerink and van der Vorst $( 1 9 7 7 )$ .

So far we have just discussed incomplete Cholesky by position. The sparsity pattern for the incomplete factor is determined in advance through the set $P$ and does not depend on the values in A. An alternative approach makes use of a drop tolerance $\tau > 0$ , which is used to determine whether or not a “potential” $h _ { i j }$ is set to zero. As an example of this strategy, suppose we compute the matrix $A _ { 1 }$ in incChol as follows:

$$
[ A _ {1} ] _ {i j}   =   \left\{ \begin{array}{l l} 0 & \text {if} | b _ {i j} - w _ {i} w _ {j} | <   \tau \sqrt {b _ {i i} b _ {j j}}  , \\ b _ {i j} - w _ {i} w _ {j} & \text {if} | b _ {i j} - w _ {i} w _ {j} | \geq \tau \sqrt {b _ {i i} b _ {j j}}  . \end{array} \right.
$$

The idea is to drop unimportant entries in the update if they are small in a relative sense. Care has to be exercised in the selection of τ so as not to induce an unacceptable level of fill-in. (Larger values of τ reduce fill-in.) The drop tolerance approach is an example of incomplete Cholesky by value.

Lin and Mor´e (1999) describe a strategy that combines the best features of incomplete Cholesky by position and incomplete Cholesky by value. Recall in gaxpy Cholesky (§4.2.5) that the triangular factor G is computed column by column. The idea is to adapt that procedure so that $H ( j { : } n , j )$ has at most $N _ { j } + p$ nonzeros, where $N _ { j }$ is the number of nonzeros in $A ( j { : } n , j )$ and p is a nonnegative integer:

for $j = 1 { : } n$

$$
v (j: n) = A (j: n, j) - H (j: n, 1: j - 1) H (j, 1: j - 1) ^ {T}
$$

$$
H (j, j) = \sqrt {v (j)}
$$

Nj = number of nonzeros in $A ( j { : } n , j )$

Set to zero each component of $v ( j + 1 ! n )$ that is not one of the $N _ { j } + p$

largest entries in $\left| v ( j { : } n ) \right|$ .

$$
H (j + 1: n, j) = v (j + 1: n) / H (j, j)
$$

end

It follows that the number of nonzeros in H is bounded by $p n + N _ { 1 } + \cdot \cdot \cdot + N _ { n }$ . Thus, the value of $p$ can be set in accordance with available memory. Note that $H ( j { : } n , j )$ is defined by the “most important” entries in $v ( j { : } n )$ . The gaxpy computation of this vector is a sparse gaxpy, and it is critical that this structure be exploited.

The incomplete factorization idea has been highly studied. Research themes include extension to LU, stability, and ways to increase the “mass” of the diagonal to guarantee existence. Particularly important has been the development of ILU() preconditioners, which control fill-in by bounding the number of times that an $a _ { i j }$ is allowed to be updated. See Benzi (2002).
