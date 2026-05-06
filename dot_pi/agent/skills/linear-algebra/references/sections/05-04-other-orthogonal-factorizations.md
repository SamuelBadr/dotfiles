# 5.4 Other Orthogonal Factorizations

Suppose $A \in \mathbb { R } ^ { m \times 4 }$ has a thin QR factorization of the following form:

$$
A = \left[ a _ {1}, a _ {2}, a _ {3}, a _ {4} \right] = \left[ q _ {1}, q _ {2}, q _ {3}, q _ {4} \right] \left[ \begin{array}{c c c c} 1 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 1 \end{array} \right].
$$

Note that ran(A) has dimension 3 but does not equal span $\{ q _ { 1 } , q _ { 2 } , q _ { 3 } \}$ , span $\{ q _ { 1 } , q _ { 2 } , q _ { 4 } \}$ , span $\{ q _ { 1 } , q _ { 3 } , q _ { 4 } \}$ , or span $\{ q _ { 2 } , q _ { 3 } , q _ { 4 } \}$ because $a _ { 4 }$ does not belong to any of these subspaces. In this case, the QR factorization reveals neither the range nor the nullspace of A and the number of nonzeros on R’s diagonal does not equal its rank. Moreover, the LS solution process based on the QR factorization (Algorithm 5.3.2) breaks down because the upper triangular portion of R is singular.

We start this section by introducing several decompositions that overcome these shortcomings. They all have the form $\mathbf { \bar { \phi } } _ { Q ^ { T } A Z } = T$ where T is a structured block triangular matrix that sheds light on A’s rank, range, and nullspace. We informally refer to matrix reductions of this form as rank revealing. See Chandrasekaren and Ipsen (1994) for a more precise formulation of the concept.

Our focus is on a modification of the QR factorization that involves column pivoting. The resulting R-matrix has a structure that supports rank estimation. To set the stage for updating methods, we briefly discus the $U L V$ and UT V frameworks Updating is discussed in §6.5 and refers to the efficient recomputation of a factorization after the matrix undergoes a low-rank change.

All these methods can be regarded as inexpensive alternatives to the SVD, which represents the “gold standard” in the area of rank determination. Nothing “takes apart” a matrix so conclusively as the SVD and so we include an explanation of its airtight reliability. The computation of the full SVD, which we discuss in §8.6, begins with the reduction to bidiagonal form using Householder matrices. Because this decomposition is important in its own right, we provide some details at the end of this section.

# 5.4.1 Numerical Rank and the SVD

Suppose $A \in \mathbb { R } ^ { m \times n }$ has SVD $U ^ { T } A V = \Sigma = \mathrm { d i a g } ( \sigma _ { i } )$ . If $\mathsf { r a n k } ( A ) = r < n$ , then according to the exact arithmetic discussion of §2.4 the singular values $\sigma _ { r + 1 } , \ldots , \sigma _ { n }$ are zero and

$$
A = \sum_ {i = 1} ^ {r} \sigma_ {k} u _ {k} v _ {k} ^ {T}. \tag {5.4.1}
$$

The exposure of rank degeneracy could not be more clear.

In Chapter 8 we describe the Golub-Kahan-Reinsch algorithm for computing the SVD. Properly implemented, it produces nearly orthogonal matrices $\widehat { U }$ and $\widehat { V }$ so that

$$
\widehat {U} ^ {T} A \widehat {V} \approx \widehat {\Sigma} = \mathrm{diag} (\widehat {\sigma} _ {1}, \ldots , \widehat {\sigma} _ {n}), \qquad \widehat {\sigma} _ {1} \geq \dots \geq \widehat {\sigma} _ {n} \geq 0.
$$

(Other SVD procedures have this property as well.) Unfortunately, unless remarkable cancellation occurs, none of the computed singular values will be zero because of roundoff error. This forces an issue. On the one hand, we can adhere to the strict mathematical definition of rank, count the number of nonzero computed singular values, and conclude from

$$
A \approx \sum_ {i = 1} ^ {n} \widehat {\sigma} _ {k} \widehat {u} _ {k} \widehat {v} _ {k} ^ {T} \tag {5.4.2}
$$

that A has full rank. However, working with every matrix as if it possessed full column rank is not particularly useful. It is more productive to liberalize the notion of rank by setting small computed singular values to zero in (5.4.2). This results in an approximation of the form

$$
A \approx \sum_ {i = 1} ^ {\widehat {r}} \widehat {\sigma} _ {k} \widehat {u} _ {k} \widehat {v} _ {k} ^ {T}, \quad \widehat {r} \leq \widehat {n} \tag {5.4.3}
$$

where we regard $\widehat { r }$ as the numerical rank. For this approach to make sense we need to guarantee that $| \widehat { \sigma } _ { i } - \sigma _ { i } |$ is small.

<For a properly implemented Golub-Kahan-Reinsch SVD algorithm, it can be shown that

$$
\widehat {U} = W + \Delta U, W ^ {T} W = I _ {m}, \quad \| \Delta U \| _ {2} \leq \epsilon ,
$$

$$
\widehat {V} = Z + \Delta V, Z ^ {T} Z = I _ {n}, \quad \| \Delta V \| _ {2} \leq \epsilon , \tag {5.4.4}
$$

$$
\widehat {\Sigma} = W ^ {T} (A + \Delta A) Z, \quad \| \Delta A \| _ {2} \leq \epsilon \| A \| _ {2},
$$

where $\epsilon$ is a small multiple of u, the machine precision. In other words, the SVD algorithm computes the singular values of a nearby matrix $A + \Delta A$ .


---

<!-- golub_300_349 -->

Note that $\widehat { U }$ and $\widehat { V }$ are not necessarily close to their exact counterparts. However, we can show that $\widehat { \sigma } _ { k }$ is close to $\sigma _ { k }$ as follows. Using Corollary 2.4.6 we have

$$
\sigma_ {k} = \min _ {\operatorname{rank} (B) = k - 1} \| A - B \| _ {2} = \min _ {\operatorname{rank} (B) = k - 1} \| (\widehat {\Sigma} - B) - E \| _ {2}
$$

where

$$
E = W ^ {T} (\Delta A) Z
$$

and

$$
\| E \| _ {2} \leq \epsilon \| A \| _ {2} = \epsilon \sigma_ {1}.
$$

Since

$$
\| \widehat {\Sigma} - B \| - \| E \| \leq \| \widehat {\Sigma} - B \| \leq \| \widehat {\Sigma} - B \| + \| E \|
$$

and

$$
\min _ {\operatorname{rank} (B) = k - 1} \| \widehat {\Sigma} _ {k} - B \| _ {2} = \widehat {\sigma} _ {k},
$$

it follows that

$$
\left| \sigma_ {k} - \widehat {\sigma} _ {k} \right| \leq \epsilon \sigma_ {1}
$$

for $k = 1 { : } n$ . Thus, if A has rank r, then we can expect $n - r$ of the computed singular values to be small. Near rank deficiency in A cannot escape detection if the SVD of A is computed.

Of course, all this hinges on having a definition of “small.” This amounts to choosing a tolerance $\delta > 0$ and declaring A to have numerical rank r if the computed singular values satisfy

$$
\widehat {\sigma} _ {1} \geq \dots \geq \widehat {\sigma} _ {\hat {r}} > \delta \geq \widehat {\sigma} _ {\hat {r} + 1} \geq \dots \geq \widehat {\sigma} _ {n}. \tag {5.4.5}
$$

We refer to the integer $\hat { r }$ as the δ-rank of A. The tolerance should be consistent with the machine precision, $\mathrm { e . g . , } \delta = \mathbf { u } \parallel A \parallel _ { \infty }$ . However, if the general level of relative error in the data is larger than u, then δ should be correspondingly bigger, e.g., $\delta = 1 0 ^ { - 2 } \parallel A \parallel _ { \infty }$ if the entries in A are correct to two digits.

For a given δ it is important to stress that, although the SVD provides a great deal of rank-related insight, it does not change the fact that the determination of numerical rank is a sensitive computation. If the gap between $\widehat { \sigma } _ { \widehat { r } }$ and $\widehat \sigma _ { r + 1 }$ is small, then A is also close (in the δ sense) to a matrix with rank $\widehat { r } - 1$ < <. Thus, the amount of confidence we have in the correctness of $\widehat { r }$ <and in how we proceed to use the approximation (5.4.2) depends on the gap between $\widehat { \sigma } _ { \widehat { r } }$ and $\widehat \sigma _ { r + 1 }$ .

# 5.4.2 QR with Column Pivoting

We now examine alternative rank-revealing strategies to the SVD starting with a modification of the Householder QR factorization procedure (Algorithm 5.2.1). In exact arithmetic, the modified algorithm computes the factorization

$$
Q ^ {T} A \Pi = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & 0 \\ r & n - r \end{array} \right] _ {m - r} ^ {r} \tag {5.4.6}
$$

where $r = \mathsf { r a n k } ( A )$ , Q is orthogonal, $R _ { 1 1 }$ is upper triangular and nonsingular, and Π is a permutation. If we have the column partitionings $A \Pi = \left[ \left. a _ { c _ { 1 } } \right| \cdot \cdot \cdot \right| \left. a _ { c _ { n } } \right]$ and $Q = [  q _ { 1 } | \cdots |  q _ { m }  ]$ , then for $k = 1 { : } n$ we have

$$
a _ {c _ {k}} = \sum_ {i = 1} ^ {\min \{r, k \}} r _ {i k} q _ {i} \in \operatorname{span} \left\{q _ {1}, \dots , q _ {r} \right\}
$$

implying

$$
\operatorname{ran} (A) = \operatorname{span} \left\{q _ {1}, \dots , q _ {r} \right\}.
$$

To see how to compute such a factorization, assume for some k that we have computed Householder matrices $H _ { 1 } , \ldots , H _ { k - 1 }$ and permutations $\Pi _ { 1 } , \ldots , \Pi _ { k - 1 }$ such that

$$
(H _ {k - 1} \dots H _ {1}) A (\Pi_ {1} \dots \Pi_ {k - 1}) = R ^ {(k - 1)} = \left[ \begin{array}{c c} R _ {1 1} ^ {(k - 1)} & R _ {1 2} ^ {(k - 1)} \\ 0 & R _ {2 2} ^ {(k - 1)} \end{array} \right] _ {m - k + 1} ^ {k - 1} \tag {5.4.7}
$$

$R _ { 1 1 } ^ { ( k - 1 ) }$ is a nonsingular and upper triangular matrix. Now suppose that

$$
R _ {2 2} ^ {(k - 1)} = \left[ z _ {k} ^ {(k - 1)} \mid \dots \mid z _ {n} ^ {(k - 1)} \right]
$$

is a column partitioning and let $p \geq k$ be the smallest index such that

$$
\left\| z _ {p} ^ {(k - 1)} \right\| _ {2} = \max \left\{\left\| z _ {k} ^ {(k - 1)} \right\| _ {2}, \dots , \left\| z _ {n} ^ {(k - 1)} \right\| _ {2} \right\}. \tag {5.4.8}
$$

Note that if $\mathsf { r a n k } ( A ) = k - 1$ , then this maximum is zero and we are finished. Otherwise, let $\Pi _ { k }$ be the $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ identity with columns p and k interchanged and determine a Householder matrix $H _ { k }$ such that if

$$
R ^ {(k)} = H _ {k} R ^ {(k - 1)} \Pi_ {k},
$$

then $R ^ { ( k ) } ( k + 1 { : } m , k ) = 0$ . In other words, $\Pi _ { k }$ moves the largest column in $R _ { 2 2 } ^ { ( k - 1 ) }$ the lead position and $H _ { k }$ zeroes all of its subdiagonal components.

The column norms do not have to be recomputed at each stage if we exploit the property

$$
Q ^ {T} z = \left[ \begin{array}{c} \alpha \\ w \end{array} \right] _ {s - 1} ^ {1} \qquad \Longrightarrow \qquad \| w \| _ {2} ^ {2} = \| z \| _ {2} ^ {2} - \alpha^ {2},
$$

which holds for any orthogonal matrix $Q \in \mathbb { R } ^ { s \times s }$ . This reduces the overhead associated with column pivoting from $O ( m n ^ { 2 } )$ flops to $O ( m n )$ flops because we can get the new column norms by updating the old column norms, e.g.,

$$
\parallel z _ {j} ^ {(k)} \parallel_ {2} ^ {2} = \parallel z _ {j} ^ {(k - 1)} \parallel_ {2} ^ {2} - r _ {k j} ^ {2} \qquad j = k + 1: n.
$$

Combining all of the above we obtain the following algorithm first presented by Businger and Golub (1965):

Algorithm 5.4.1 (Householder QR With Column Pivoting) Given $A \in \mathbb { R } ^ { m \times n }$ with m $\geq n$ , the following algorithm computes $r = \mathsf { r a n k } ( A )$ and the factorization (5.4.6) with $Q = H _ { 1 } \cdot \cdot \cdot H _ { r }$ and $\Pi = \Pi _ { 1 } \cdots \Pi _ { r }$ . The upper triangular part of A is overwritten by the upper triangular part of R and components $j + 1 { : } m$ of the jth Householder vector are stored in $A ( j + 1 { : } m , j )$ . The permutation Π is encoded in an integer vector piv. In particular, $\Pi _ { j }$ is the identity with rows $j$ and $p i v ( j )$ interchanged.

for $j = 1:n$ $c(j) = A(1:m,j)^T A(1:m,j)$ end $r = 0$ $\tau = \max\{c(1),\ldots,c(n)\}$ while $\tau > 0$ and $r < n$ $r = r + 1$ Find smallest $k$ with $r \leq k \leq n$ so $c(k) = \tau$ . $piv(r) = k$ $A(1:m,r) \leftrightarrow A(1:m,k)$ $c(r) \leftrightarrow c(k)$ $[v,\beta] = \text{house}(A(r:m,r))$ $A(r:m,r:n) = (I_{m-r+1} - \beta vv^T)A(:r:m,r:n)$ $A(r+1:m,r) = v(2:m-r+1)$ for $i = r+1:n$ $c(i) = c(i) - A(r,i)^2$ end $\tau = \max\{c(r+1),\ldots,c(n)\}$ end

This algorithm requires 4mnr $- 2 r ^ { 2 } ( m + n ) + 4 r ^ { 3 } / 3$ flops where r = rank(A).

# 5.4.3 Numerical Rank and A  QR

In principle, QR with column pivoting reveals rank. But how informative is the method in the context of floating point arithmetic? After k steps we have

$$
\mathsf {f l} (H _ {k} \dots H _ {1} A \Pi_ {1} \dots \Pi_ {k}) = \widehat {R} ^ {(k)} = \left[ \begin{array}{c c} \widehat {R} _ {1 1} ^ {(k)} & \widehat {R} _ {1 2} ^ {(k)} \\ 0 & \widehat {R} _ {2 2} ^ {(k)} \end{array} \right] _ {m - k} ^ {k}. \tag {5.4.9}
$$

If ${ \widehat { R } } _ { 2 2 } ^ { ( k ) }$ is suitably small in norm, then it is reasonable to terminate the reduction and declare A to have rank k. A typical termination criteria might be

$$
\| \widehat {R} _ {2 2} ^ {(k)} \| _ {2} \leq \epsilon_ {1} \| A \| _ {2}
$$

for some small machine-dependent parameter $\epsilon _ { 1 }$ . In view of the roundoff properties associated with Householder matrix computation (cf. §5.1.12), we know that $\widehat { R } ^ { ( k ) }$ is the exact R-factor of a matrix $A + E _ { k }$ , where

$$
\| E _ {k} \| _ {2} \leq \epsilon_ {2} \| A \| _ {2}, \quad \epsilon_ {2} = O (\mathbf {u}).
$$

Using Corollary 2.4.4 we have

$$
\sigma_ {k + 1} (A + E _ {k}) = \sigma_ {k + 1} (\widehat {R} ^ {(k)}) \leq \| \widehat {R} _ {2 2} ^ {(k)} \| _ {2}.
$$

Since $\sigma _ { k + 1 } ( A ) \leq \sigma _ { k + 1 } ( A + E _ { k } ) + \| E _ { k } \| _ { 2 }$ , it follows that

$$
\sigma_ {k + 1} (A) \leq (\epsilon_ {1} + \epsilon_ {2}) \| A \| _ {2}.
$$

In other words, a relative perturbation of $O ( \epsilon _ { 1 } + \epsilon _ { 2 } )$ in A can yield a rank-k matrix. With this termination criterion, we conclude that QR with column pivoting discovers rank deficiency if ${ \widehat { R } } _ { 2 2 } ^ { ( k ) }$ is small for some $k < n$ . However, it does not follow that the matrix ${ \widehat { R } } _ { 2 2 } ^ { ( k ) }$ <22 in (5.4.9) is small if rank $( A ) = k$ . There are examples of nearly rank deficient matrices whose R-factor look perfectly “normal.” A famous example is the Kahan matrix

$$
\mathsf {K a h} _ {n} (s) = \operatorname{diag} (1, s, \ldots , s ^ {n - 1}) \left[ \begin{array}{c c c c c} 1 & - c & - c & \dots & - c \\ 0 & 1 & - c & \dots & - c \\ & \ddots & & \vdots & \vdots \\ \vdots & & & 1 & - c \\ 0 & & \dots & & 1 \end{array} \right].
$$

Here, $c ^ { 2 } + s ^ { 2 } = 1$ with $c , s > 0$ . (See Lawson and Hanson (SLS, p. 31).) These matrices are unaltered by Algorithm 5.4.1 and thus $\parallel R _ { 2 2 } ^ { ( k ) } \parallel _ { 2 } \geq s ^ { n - 1 }$ for $k = 1 { : } n - 1$ . This inequality implies (for example) that the matrix $\mathsf { K a h } _ { 3 0 0 } ( . 9 9 )$ has no particularly small trailing principal submatrix since $s ^ { 2 9 9 } \approx . 0 5$ . However, a calculation shows that σ300 $= O ( 1 0 ^ { - 1 9 } )$ .

Nevertheless, in practice, small trailing R-submatrices almost always emerge that correlate well with the underlying rank. In other words, it is almost always the case that ${ \widehat { R } } _ { 2 2 } ^ { ( k ) }$ is small if A has rank k.

# 5.4.4 Finding a Good Column Ordering

It is important to appreciate that Algorithm 5.4.1 is just one way to determine the column permutation Π. The following result sets the stage for a better way.

Theorem 5.4.1. If $A \in \mathbb { R } ^ { m \times n }$ and $v \in \mathbb { R } ^ { n }$ is a unit 2-norm vector, then there exists a permutation Π so that the QR factorization

$$
A \Pi = Q R
$$

satisfies $| r _ { n n } | \leq { \sqrt { n } } \sigma$ where $\sigma = \parallel A v \parallel _ { 2 }$ .

Proof. Suppose $\Pi \in \mathbb { R } ^ { n \times n }$ is a permutation such that if $w = \Pi ^ { T } v$ , then

$$
\left| w _ {n} \right| = \max \left| v _ {i} \right|.
$$

Since $w _ { n }$ is the largest component of a unit 2-norm vector, $| w _ { n } | \geq 1 / { \sqrt { n } }$ . If $A \Pi = Q R$ is a QR factorization, then

$$
\sigma = \parallel A v \parallel_ {2} = \parallel (Q ^ {T} A \Pi) (\Pi^ {T} v) \parallel_ {2} = \parallel R (1: n, 1: n) w \parallel_ {2} \geq | r _ {n n} w _ {n} | \geq | r _ {n n} | / \sqrt {n}. \quad \square
$$

Note that if $v = v _ { n }$ is the right singular vector corresponding to $\sigma _ { \mathrm { m i n } } ( A )$ , then $| r _ { n n } | \leq$ $\sqrt { n } \sigma _ { n }$ . This suggests a framework whereby the column permutation matrix Π is based on an estimate of $v _ { n }$ :

Step 1. Compute the QR factorization $A = Q _ { 0 } R _ { 0 }$ and note that $R _ { 0 }$ has the same right singular vectors as A.

Step 2. Use condition estimation techniques to obtain a unit vector v with $\parallel R _ { 0 } v \parallel _ { 2 } \approx \sigma _ { n }$ .

Step 3. Determine Π and the QR factorization $A \Pi = Q R .$ .

See Chan (1987) for details about this approach to rank determination. The permutation Π can be generated as a sequence of swap permutations. This supports a very economical Givens rotation method for generating of Q and R from $Q _ { 0 }$ and $R _ { 0 }$ .

# 5.4.5 More General Rank-Revealing Decompositions

Additional rank-revealing strategies emerge if we allow general orthogonal recombinations of the A’s columns instead of just permutations. That is, we look for an orthogonal Z so that the QR factorization

$$
A Z = Q R
$$

produces a rank-revealing R. To impart the spirit of this type of matrix reduction, we show how the rank-revealing properties of a given $A Z = Q R$ factorization can be improved by replacing $Z , Q$ , and R with

$$
Z _ {\mathrm{new}} = Z Z _ {G}, \qquad Q _ {\mathrm{new}} = Q Q _ {G}, \qquad R _ {\mathrm{new}} = Q _ {G} ^ {T} R Z _ {G},
$$

respectively, where $Q _ { G }$ and $Z _ { G }$ are products of Givens rotations and $R _ { \mathrm { n e w } }$ is upper triangular. The rotations are generated by introducing zeros into a unit 2-norm nvector v which we assume approximates the n-th right singular vector of AZ. In particular, if $Z _ { \cal G } ^ { T } v = e _ { n } = I _ { n } ( : , n )$ and $\parallel R v \parallel _ { 2 } \approx \sigma _ { n }$ , then

$$
\parallel R _ {\mathrm{new}} e _ {n} \parallel_ {2} = \parallel Q _ {G} ^ {T} R Z _ {G} e _ {n} \parallel_ {2} = \parallel Q _ {G} ^ {T} R v \parallel_ {2} = \parallel R v \parallel_ {2} \approx \sigma_ {n}
$$

This says that the norm of the last column of $R _ { \mathrm { n e w } }$ is approximately the smallest singular value of A, which is certainly one way to reveal the underlying matrix rank.

We use the case $n = 4$ to illustrate how the Givens rotations arise and why the overall process is economical. Because we are transforming v to $e _ { n }$ and not $e _ { 1 }$ , we need to “flip” the mission of the 2-by-2 rotations in the $Z _ { G }$ computations so that top components are zeroed, i.e.,

$$
\left[ \begin{array}{c} 0 \\ \times \end{array} \right] = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] \left[ \begin{array}{c} \times \\ \times \end{array} \right].
$$

This requires only a slight modification of Algorithm 5.1.3.

In the n = 4 case we start with

$$
R = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right] \qquad v = \left[ \begin{array}{l} \times \\ \times \\ \times \\ \times \end{array} \right]
$$

and proceed to compute

$$
Z _ {G} = G _ {1 2} G _ {2 3} G _ {3 4}
$$

and

$$
Q _ {G} = H _ {1 2} H _ {2 3} H _ {3 4}
$$

as products of Givens rotations. The first step is to zero the top component of v with a “flipped” (1,2) rotation and update R accordingly:

$$
R \leftarrow R G _ {1 2} = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right], \qquad v \leftarrow G _ {1 2} ^ {T} v = \left[ \begin{array}{l} 0 \\ \times \\ \times \\ \times \end{array} \right].
$$

To remove the unwanted subdiagonal in R, we apply a conventional (nonflipped) Givens rotation from the left to R (but not v):

$$
R \leftarrow H _ {1 2} ^ {T} R = \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right], \qquad \qquad v = \left[ \begin{array}{c} 0 \\ \times \\ \times \\ \times \end{array} \right].
$$

The next step is analogous:

$$
R \leftarrow R G _ {2 3} = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right], \qquad v \leftarrow G _ {2 3} ^ {T} v = \left[ \begin{array}{l} 0 \\ 0 \\ \times \\ \times \end{array} \right].
$$

$$
R \gets H _ {2 3} ^ {T} R = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right], \qquad \qquad v = \left[ \begin{array}{l} 0 \\ 0 \\ \times \\ \times \end{array} \right].
$$

And finally,

$$
R \gets R G _ {3 4} = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \end{array} \right], \qquad v = G _ {3 4} ^ {T} v = \left[ \begin{array}{l} 0 \\ 0 \\ 0 \\ \times \end{array} \right],
$$

$$
R \gets H _ {3 4} ^ {T} R = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right], \qquad \qquad v = \left[ \begin{array}{l} 0 \\ 0 \\ 0 \\ \times \end{array} \right].
$$

The pattern is clear, for $i = 1 { : } n - 1$ , a $G _ { i , i + 1 }$ is used to zero the current $v _ { i }$ and an $H _ { i , i + 1 }$ is used to zero the current $r _ { i + 1 , i } .$ . The overall transition from $\{ Q , ~ Z , ~ R \}$ to $\{ Q _ { \mathrm { n e w } } , Z _ { \mathrm { n e w } } , R _ { \mathrm { n e w } } \}$ involves $O ( m n )$ flops. If the Givens rotations are kept in factored form, this flop count is reduced to $O ( n ^ { 2 } )$ . We mention that the ideas in this subsection can be iterated to develop matrix reductions that expose the structure of matrices whose rank is less than $n - 1$ . “Zero-chasing” with Givens rotations is at the heart of many important matrix algorithms; see §6.3, §7.5, and §8.3.

# 5.4.6 The UT V Framework

As mentioned at the start of this section, we are interested in factorizations that are cheaper than the SVD but which provide the same high quality information about rank, range, and nullspace. Factorizations of this type are referred to as UTV factorizations where the $^ { 6 6 } \mathrm { T } ^ { 9 }$ stands for triangular and the $^ { 6 } U ^ { , s }$ and ${ } ^ { 6 6 } V ^ { 5 9 }$ remind us of the SVD and orthogonal U and V matrices of singular vectors.

The matrix T can be upper triangular (these are the URV factorizations) or lower triangular (these are the ULV factorizations). It turns out that in a particular application one may favor a URV approach over a ULV approach, see §6.3. Moreover, the two reductions have different approximation properties. For example, suppose $\sigma _ { k } ( A ) > \sigma _ { k + 1 } ( A )$ and S is the subspace spanned by A’s right singular vectors $v _ { k + 1 } , \ldots , v _ { n }$ . Think of S as an approximate nullspace of A. Following Stewart (1993), if

$$
U ^ {T} A V = R = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & R _ {2 2} \end{array} \right] _ {m - k} ^ {k}
$$

and $V = \left[ \left. V _ { 1 } \right| V _ { 2 } \right]$ is partitioned conformably, then

$$
\text { dist } (\text { ran } (V _ {2}), S) \leq \frac {\| R _ {1 2} \| _ {2}}{(1 - \rho_ {R} ^ {2}) \sigma_ {\min} (R _ {1 1})} \tag {5.4.10}
$$

where

$$
\rho_ {R} = \frac {\parallel R _ {2 2} \parallel_ {2}}{\sigma_ {\mathrm{min}} (R _ {1 1})}
$$

is assumed to be less than 1. On the other hand, in the ULV setting we have

$$
U ^ {T} A V = L = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \end{array} \right] _ {m - k} ^ {k}.
$$

If $V = \left[ \left. V _ { 1 } \right| V _ { 2 } \right]$ is partitioned conformably, then

$$
\operatorname{dist} \left(\operatorname{ran} \left(V _ {2}\right), S\right) \leq \rho_ {L} \frac {\left\| L _ {1 2} \right\| _ {2}}{\left(1 - \rho_ {L} ^ {2}\right) \sigma_ {\min} \left(L _ {1 1}\right)} \tag {5.4.11}
$$

where

$$
\rho_ {L} = \frac {\parallel L _ {2 2} \parallel_ {2}}{\sigma_ {\mathrm{min}} (L _ {1 1})}
$$

is also assumed to be less than 1. However, in practice the ρ-factors in both (5.4.10) and (5.4.11) are often much less than 1. Observe that when this is the case, the upper bound in (5.4.11) is much smaller than the upper bound in (5.4.10).

# 5.4.7 Complete Orthogonal Decompositions

Related to the UTV framework is the idea of a complete orthogonal factorization. Here we compute orthogonal U and V such that

$$
U ^ {T} A V = \left[ \begin{array}{c c} T _ {1 1} & 0 \\ 0 & 0 \end{array} \right] _ {m - r} ^ {r} \tag {5.4.12}
$$

where $r = { \mathrm { r a n k } } ( A )$ . The SVD is obviously an example of a decomposition that has this structure. However, a cheaper, two-step QR process is also possible. We first use Algorithm 5.4.1 to compute

$$
U ^ {T} A \Pi = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & 0 \\ r & n - r \end{array} \right] _ {m - r} ^ {r}
$$

and then follow up with a second QR factorization

$$
Q ^ {T} \left[ \begin{array}{l} R _ {1 1} ^ {T} \\ R _ {1 2} ^ {T} \end{array} \right] = \left[ \begin{array}{l} S _ {1} \\ 0 \end{array} \right]
$$

via Algorithm 5.2.1. If we set $V = \Pi Q$ , then (5.4.12) is realized with $T _ { 1 1 } = S _ { 1 } ^ { T }$ . Note that two important subspaces are defined by selected columns of $U = [  u _ { 1 } | \cdot \cdot \cdot | u _ { m } ]$ and $V = [  v _ { 1 } | \cdot \cdot \cdot | v _ { n } ]$ :

$$
\operatorname{ran} (A) = \operatorname{span} \left\{u _ {1}, \dots , u _ {r} \right\},
$$

$$
\operatorname{null} (A) = \operatorname{span} \left\{v _ {r + 1}, \dots , v _ {n} \right\}.
$$

Of course, the computation of a complete orthogonal decomposition in practice would require the careful handling of numerical rank.

# 5.4.8 Bidiagonalization

There is one other two-sided orthogonal factorization that is important to discuss and that is the bidiagonal factorization. It is not a rank-revealing factorization per se, but it has a useful role to play because it rivals the SVD in terms of data compression.

Suppose $A \in \mathbb { R } ^ { m \times n }$ and $m \geq n$ . The idea is to compute orthogonal $U _ { B } \ ( m { \mathrm { - b y - } } m )$ and $V _ { B } \ ( n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n )$ such that

$$
U _ {B} ^ {T} A V _ {B} = \left[ \begin{array}{c c c c c} d _ {1} & f _ {1} & 0 & \dots & 0 \\ 0 & d _ {2} & f _ {2} & & 0 \\ \vdots & \ddots & \ddots & \ddots & \vdots \\ 0 & \dots & & d _ {n - 1} & f _ {n - 1} \\ 0 & \dots & & 0 & d _ {n} \\ \hline & & 0 \end{array} \right]. \tag {5.4.13}
$$

$U _ { B } = U _ { 1 } \cdot \cdot \cdot U _ { n }$ and $V _ { B } = V _ { 1 } \cdot \cdot \cdot V _ { n - 2 }$ can each be determined as a product of Householder matrices, e.g.,

$$
\left[ \begin{array}{l l l l} \times & \times & \times & \times \\ \times & \times & \times & \times \\ \times & \times & \times & \times \\ \times & \times & \times & \times \\ \times & \times & \times & \times \end{array} \right] \xrightarrow {U _ {1}} \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \end{array} \right] \xrightarrow {V _ {1}}
$$

$$
\left[ \begin{array}{c c c c} \times & \times & 0 & 0 \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \end{array} \right] \xrightarrow {U _ {2}} \left[ \begin{array}{c c c c} \times & \times & 0 & 0 \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \end{array} \right] \xrightarrow {V _ {2}}
$$

$$
\left[ \begin{array}{c c c c} \times & \times & 0 & 0 \\ 0 & \times & \times & 0 \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \end{array} \right] \xrightarrow {U _ {3}} \left[ \begin{array}{c c c c} \times & \times & 0 & 0 \\ 0 & \times & \times & 0 \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \\ 0 & 0 & 0 & \times \end{array} \right] \xrightarrow {U _ {4}} \left[ \begin{array}{c c c c} \times & \times & 0 & 0 \\ 0 & \times & \times & 0 \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 \end{array} \right].
$$

In general, $U _ { k }$ introduces zeros into the kth column, while $V _ { k }$ zeros the appropriate entries in row k. Overall we have:

Algorithm 5.4.2 (Householder Bidiagonalization) Given $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ , the following algorithm overwrites A with $U _ { s } ^ { T } A V _ { B } = \mathbf { \dot { B } }$ where B is upper bidiagonal and $U _ { B } = U _ { 1 } \cdot \cdot \cdot U _ { n }$ and $V _ { B } = V _ { 1 } \cdots V _ { n - 2 }$ . The essential part of $U _ { j } \mathrm { ^ { 5 } s }$ Householder vector is stored in $A ( j + 1 { : } m , j )$ and the essential part of $V _ { j } \mathrm { { ^ { s } } }$ Householder vector is stored in $A ( j , j + 2 { : } n )$ .

for $j = 1:n$ [ [v, \beta] = \text{house}(A(j:m, j)) ] $A(j:m, j:n) = (I_{m-j+1} - \beta vv^T)A(j:m, j:n)$ $A(j+1:m, j) = v(2:m - j+1)$ if $j \leq n-2$ [ [v, \beta] = \text{house}(A(j, j+1:n)^T) ] $A(j:m, j+1:n) = A(j:m, j+1:n)(I_{n-j} - \beta vv^T)$ $A(j, j+2:n) = v(2:n-j)^T$ end

end

This algorithm requires $4 m n ^ { 2 } - 4 n ^ { 3 } / 3$ flops. Such a technique is used by Golub and Kahan (1965), where bidiagonalization is first described. If the matrices $U _ { B }$ and $V _ { B }$ are explicitly desired, then they can be accumulated in $4 m ^ { 2 } n - 4 n ^ { 3 } / 3$ and $4 n ^ { 3 } / 3$ flops, respectively. The bidiagonalization of A is related to the tridiagonalization of $A ^ { T } A$ . See §8.3.1.

# 5.4.9 R-Bidiagonalization

If $m \gg n$ , then a faster method of bidiagonalization method results if we upper triangularize A first before applying Algorithm 5.4.2. In particular, suppose we compute an orthogonal $Q \in \mathbb { R } ^ { m \times m }$ such that

$$
Q ^ {T} A = \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right]
$$

is upper triangular. We then bidiagonalize the square matrix $R _ { 1 }$ ,

$$
U _ {R} ^ {T} R _ {1} V _ {B} = B _ {1},
$$

where $U _ { R }$ and $V _ { B }$ are orthogonal. If $U _ { B } = Q \mathrm { d i a g } \left( U _ { R } , I _ { m - n } \right)$ , then

$$
U ^ {T} A V = \left[ \begin{array}{c} B _ {1} \\ 0 \end{array} \right] \equiv B
$$

is a bidiagonalization of A.

The idea of computing the bidiagonalization in this manner is mentioned by Lawson and Hanson (SLS, p. 119) and more fully analyzed by Chan (1982). We refer to this method as R-bidiagonalization and it requires $( 2 m n ^ { 2 } + 2 n ^ { 3 } )$ flops. This is less than the flop count for Algorithm 5.4.2 whenever $m \geq 5 n / 3$ .

# Problems

P5.4.1 Let x, $y \in \mathbb { R } ^ { m }$ and $Q \in \mathbb { R } ^ { m \times m }$ be given with Q orthogonal. Show that if

$$
Q ^ {T} x = \left[ \begin{array}{c} \alpha \\ u \end{array} \right] _ {m - 1} ^ {1}, \qquad Q ^ {T} y = \left[ \begin{array}{c} \beta \\ v \end{array} \right] _ {m - 1} ^ {1}
$$

then $u ^ { T } v = x ^ { T } y - \alpha \beta$ .

P5.4.2 Let $A = [ a _ { 1 } \mid \cdots \mid a _ { n } ] \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ be given. For any column subset $\{ a _ { c _ { 1 } } , \ldots , a _ { c _ { k } } \}$ define

$$
\operatorname{res} \left(\left[ a _ {c _ {1}} \mid \dots \mid a _ {c _ {k}} \right]\right) = \min _ {x \in \mathbf {R} ^ {k}} \| \left[ a _ {c _ {1}} \mid \dots \mid a _ {c _ {k}} \right] x - b \| _ {2}
$$

Describe an alternative pivot selection procedure for Algorithm 5.4.1 such that if $Q R \ = \ A \Pi \ =$ $\left[ \boldsymbol { a } _ { c 1 } \left| \cdots \right| \boldsymbol { a } _ { c n } \right]$ in the final factorization, then for $k = 1 { : } n { : }$

$$
\operatorname{res} \left(\left[ a _ {c _ {1}} \mid \dots \mid a _ {c _ {k}} \right]\right) = \min _ {i \geq k} \operatorname{res} \left(\left[ a _ {c _ {1}}, \ldots , a _ {c _ {k - 1}}, a _ {c _ {i}} \right]\right).
$$

P5.4.3 Suppose $T \in \mathbb { R } ^ { n \times n }$ is upper triangular and $t _ { k k } = \sigma _ { m i n } ( T )$ . Show that $T ( 1 { : } k - 1 , k ) = 0$ and $T ( k , k + 1 { : } n ) = 0$ .

P5.4.4 Suppose $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ . Give an algorithm that uses Householder matrices to compute an orthogonal $Q \in \mathbb { R } ^ { m \times m }$ so that if $Q ^ { T } A = L ,$ , then $L ( n + 1 { : } m , : ) = 0$ and $L ( 1 { : } n , 1 { : } n )$ is lower triangular.

P5.4.5 Suppose $R \in \mathbb { R } ^ { n \times n }$ is upper triangular and $Y \in \mathbb { R } ^ { n \times j }$ has orthonormal columns and satisfies $\| \ R Y \| _ { 2 } = \sigma ,$ . Give an algorithm that computes orthogonal U and V , each products of Givens rotations, so that $\mathbf { \dot { \zeta } } U ^ { T } R V = R _ { \mathrm { n e w } }$ is upper triangular and $V ^ { T } \check { Y } = Y _ { \mathrm { n e w } }$ has the property that

$$
Y _ {\text { new }} (n - j + 1: n,:) = \text { diag } (\pm 1).
$$

What can you say about $R _ { \mathrm { n e w } } ( n - j + 1 { : } n , n - j + 1 { : } n ) ?$

P5.4.6 Give an algorithm for reducing a complex matrix A to real bidiagonal form using complex Householder transformations.

P5.4.7 Suppose $B \in \mathbb { R } ^ { n \times n }$ is upper bidiagonal with $b _ { n n } = 0 $ . Show how to construct orthogonal U and V (product of Givens rotations) so that $U ^ { T } B V$ is upper bidiagonal with a zero nth column.

P5.4.8 Suppose $A \in \mathbb { R } ^ { m \times n }$ with $m < n ,$ . Give an algorithm for computing the factorization

$$
U ^ {T} A V = \left[ \begin{array}{c c} B & O \end{array} \right]
$$

where B is an m-by-m upper bidiagonal matrix. (Hint: Obtain the form

$$
\left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 \\ 0 & 0 & \times & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & 0 \end{array} \right].
$$

using Householder matrices and then “chase” the $( m , m + 1 )$ entry up the (m+1)st column by applying Givens rotations from the right.)

P5.4.9 Show how to efficiently bidiagonalize an n-by-n upper triangular matrix using Givens rotations.

P5.4.10 Show how to upper bidiagonalize a tridiagonal matrix $T \in \mathbb { R } ^ { n \times n }$ using Givens rotations.

P5.4.11 Show that if $B \in \mathbb { R } ^ { n \times n }$ is an upper bidiagonal matrix having a repeated singular value, then B must have a zero on its diagonal or superdiagonal.

# Notes and References for §5.4

QR with column pivoting was first discussed in:

P.A. Businger and G.H. Golub (1965). “Linear Least Squares Solutions by Householder Transformations,” Numer. Math. 7, 269–276.

In matters that concern rank deficiency, it is helpful to obtain information about the smallest singular value of the upper triangular matrix R. This can be done using the techniques of 3.5.4 or those that are discussed in:

I. Karasalo (1974). “A Criterion for Truncation of the QR Decomposition Algorithm for the Singular Linear Least Squares Problem,” BIT 14, 156–166.

N. Anderson and I. Karasalo (1975). “On Computing Bounds for the Least Singular Value of a Triangular Matrix,” BIT 15, 1–4.

C.-T. Pan and P.T.P. Tang (1999). “Bounds on Singular Values Revealed by QR Factorizations,” BIT 39, 740–756.   
C.H. Bischof (1990). “Incremental Condition Estimation,” SIAM J. Matrix Anal. Applic., 11, 312– 322.

Revealing the rank of a matrix through a carefully implementated factorization has prompted a great deal of research, see:   
T.F. Chan (1987). “Rank Revealing QR Factorizations,” Lin. Alg. Applic. 88/89, 67–82.   
T.F. Chan and P. Hansen (1992). “Some Applications of the Rank Revealing QR Factorization,” SIAM J. Sci. Stat. Comp. 13, 727–741.   
S. Chandrasekaren and I.C.F. Ipsen (1994). “On Rank-Revealing Factorizations,” SIAM J. Matrix Anal. Applic. 15, 592–622.   
M. Gu and S.C. Eisenstat (1996). “Efficient Algorithms for Computing a Strong Rank-Revealing QR Factorization,” SIAM J. Sci. Comput. 17, 848–869.   
G.W. Stewart (1999). “The QLP Approximation to the Singular Value Decomposition,” SIAM J. Sci. Comput. 20, 1336–1348.   
D.A. Huckaby and T.F. Chan (2005). “Stewart’s Pivoted QLP Decomposition for Low-Rank Matrices,” Num. Lin. Alg. Applic. 12, 153–159.   
A. Dax (2008). “Orthogonalization via Deflation: A Minimum Norm Approach to Low-Rank Approximation of a Matrix,” SIAM J. Matrix Anal. Applic. 30, 236–260.   
Z. Drma˘c and Z. Bujanovi˘c (2008). “On the Failure of Rank-Revealing QR Factorization Software—A Case Study,” ACM Trans. Math. Softw. 35, Article 12.   
We have more to say about the UTV framework in §6.5 where updating is discussed. Basic references for what we cover in this section include:   
G.W. Stewart (1993). “UTV Decompositions,” in Numerical Analysis 1993, Proceedings of the 15th Dundee Conference, June–July 1993, Longman Scientic & Technical, Harlow, Essex, UK, 225–236.   
P.A. Yoon and J.L. Barlow (1998) “An Efficient Rank Detection Procedure for Modifying the ULV Decomposition,” BIT 38, 781–801.   
J.L. Barlow, H. Erbay, and I. Slapnicar (2005). “An Alternative Algorithm for the Refinement of ULV Decompositions,” SIAM J. Matrix Anal. Applic. 27, 198–211.   
Column-pivoting makes it more difficult to achieve high performance when computing the QR factorization. However, it can be done:   
C.H. Bischof and P.C. Hansen (1992). “A Block Algorithm for Computing Rank-Revealing QR Factorizations,” Numer. Algorithms 2, 371-392.   
C.H. Bischof and G. Quintana-Orti (1998). “Computing Rank-revealing QR factorizations of Dense Matrices,” ACM Trans. Math. Softw. 24, 226–253.   
C.H. Bischof and G. Quintana-Orti (1998). “Algorithm 782: Codes for Rank-Revealing QR factorizations of Dense Matrices,” ACM Trans. Math. Softw. 24, 254–257.   
G. Quintana-Orti, X. Sun, and C.H. Bischof (1998). “A BLAS–3 Version of the QR Factorization with Column Pivoting,” SIAM J. Sci. Comput. 19, 1486–1494.   
A carefully designed LU factorization can also be used to shed light on matrix rank:   
T-M. Hwang, W-W. Lin, and E.K. Yang (1992). “Rank-Revealing LU Factorizations,” Lin. Alg. Applic. 175, 115–141.   
T.-M. Hwang, W.-W. Lin and D. Pierce (1997). “Improved Bound for Rank Revealing LU Factorizations,” Lin. Alg. Applic. 261, 173–186.   
L. Miranian and M. Gu (2003). “Strong Rank Revealing LU Factorizations,” Lin. Alg. Applic. 367, 1–16.   
Column pivoting can be incorporated into the modified Gram-Schmidt process, see:   
A. Dax (2000). “A Modified Gram-Schmidt Algorithm with Iterative Orthogonalization and Column Pivoting,” Lin. Alg. Applic. 310, 25–42.   
M. Wei and Q. Liu (2003). “Roundoff Error Estimates of the Modified GramSchmidt Algorithm with Column Pivoting,” BIT 43, 627–645.   
Aspects of the complete orthogonal decomposition are discussed in:

R.J. Hanson and C.L. Lawson (1969). “Extensions and Applications of the Householder Algorithm for Solving Linear Least Square Problems,” Math. Comput. 23, 787–812.   
P.A. Wedin (1973). “On the Almost Rank-Deficient Case of the Least Squares Problem,” BIT 13, 344–354.   
G.H. Golub and V. Pereyra (1976). “Differentiation of Pseudo-Inverses, Separable Nonlinear Least Squares Problems and Other Tales,” in Generalized Inverses and Applications, M.Z. Nashed (ed.), Academic Press, New York, 303–324.   
The quality of the subspaces that are exposed through a complete orthogonal decomposition are analyzed in:   
R.D. Fierro and J.R. Bunch (1995). “Bounding the Subspaces from Rank Revealing Two-Sided Orthogonal Decompositions,” SIAM J. Matrix Anal. Applic. 16, 743–759.   
R.D. Fierro (1996). “Perturbation Analysis for Two-Sided (or Complete) Orthogonal Decompositions,” SIAM J. Matrix Anal. Applic. 17, 383–400.   
The bidiagonalization is a particularly important decomposition because it typically precedes the computation of the SVD as we discuss in §8.6. Thus, there has been a strong research interest in its efficient and accurate computation:   
B. Lang (1996). “Parallel Reduction of Banded Matrices to Bidiagonal Form,” Parallel Comput. 22, 1–18.   
J.L. Barlow (2002). “More Accurate Bidiagonal Reduction for Computing the Singular Value Decomposition,” SIAM J. Matrix Anal. Applic. 23, 761–798.   
J.L. Barlow, N. Bosner and Z. Drma˘c (2005). “A New Stable Bidiagonal Reduction Algorithm,” Lin. Alg. Applic. 397, 35–84.   
B.N. Parlett (2005). “A Bidiagonal Matrix Determines Its Hyperbolic SVD to Varied Relative Accuracy,” SIAM J. Matrix Anal. Applic. 26, 1022–1057.   
N. Bosner and J.L. Barlow (2007). “Block and Parallel Versions of One-Sided Bidiagonalization,” SIAM J. Matrix Anal. Applic. 29, 927–953.   
G.W. Howell, J.W. Demmel, C.T. Fulton, S. Hammarling, and K. Marmol (2008). “Cache Efficient Bidiagonalization Using BLAS 2.5 Operators,” ACM Trans. Math. Softw. 34, Article 14.   
H. Ltaief, J. Kurzak, and J. Dongarra (2010). “Parallel Two-Sided Matrix Reduction to Band Bidiagonal Form on Multicore Architectures,” IEEE Trans. Parallel Distrib. Syst. 21, 417–423.
