# 5.5 The Rank-Deficient Least Squares Problem

If A is rank deficient, then there are an infinite number of solutions to the LS problem. We must resort to techniques that incorporate numerical rank determination and identify a particular solution as “special.” In this section we focus on using the SVD to compute the minimum norm solution and QR-with-column-pivoting to compute what is called the basic solution. Both of these approaches have their merits and we conclude with a subset selection procedure that combines their positive attributes.

# 5.5.1 The Minimum Norm Solution

Suppose $A \in \mathbb { R } ^ { m \times n }$ and rank $( A ) = r < n$ . The rank-deficient LS problem has an infinite number of solutions, for if x is a minimizer and $z \in \mathsf { n u l l } ( A )$ , then $x + z$ is also a minimizer. The set of all minimizers

$$
\mathcal {X} = \left\{x \in \mathbb {R} ^ {n}: \| A x - b \| _ {2} = \min \right\}
$$

is convex and so if $x _ { 1 } , x _ { 2 } \in \mathcal { X }$ and $\lambda \in [ 0 , 1 ]$ , then

$$
\| A (\lambda x _ {1} + (1 - \lambda) x _ {2}) - b \| _ {2} \leq \lambda \| A x _ {1} - b \| _ {2} + (1 - \lambda) \| A x _ {2} - b \| _ {2} = \min _ {x \in \mathbb {R} ^ {n}} \| A x - b \| _ {2}.
$$

Thus, $\lambda x _ { 1 } + ( 1 - \lambda ) x _ { 2 } \in \mathcal { X }$ . It follows that X has a unique element having minimum 2-norm and we denote this solution by $x _ { L S }$ . (Note that in the full-rank case, there is only one LS solution and so it must have minimal 2-norm. Thus, we are consistent with the notation in §5.3.)

Any complete orthogonal factorization (§5.4.7) can be used to compute $x _ { L S }$ . In particular, if Q and $Z$ are orthogonal matrices such that

$$
Q ^ {T} A Z = T = \left[ \begin{array}{c c} T _ {1 1} & 0 \\ 0 & 0 \end{array} \right] _ {m - r} ^ {r}, r = \operatorname{rank} (A)
$$

then

$$
\| A x - b \| _ {2} ^ {2} = \| (Q ^ {T} A Z) Z ^ {T} x - Q ^ {T} b \| _ {2} ^ {2} = \| T _ {1 1} w - c \| _ {2} ^ {2} + \| d \| _ {2} ^ {2}
$$

where

$$
Z ^ {T} x = \left[ \begin{array}{l} w \\ y \end{array} \right] _ {n - r} ^ {r}, \qquad Q ^ {T} b = \left[ \begin{array}{l} c \\ d \end{array} \right] _ {m - r} ^ {r}.
$$

Clearly, if x is to minimize the sum of squares, then we must have $w = T _ { 1 1 } ^ { - 1 } c$ . For x to have minimal 2-norm, y must be zero, and thus

$$
x _ {L S} = Z \left[ \begin{array}{c} T _ {1 1} ^ {- 1} c \\ 0 \end{array} \right].
$$

Of course, the SVD is a particularly revealing complete orthogonal decomposition. It provides a neat expression for $x _ { L S }$ and the norm of the minimum residual $\rho _ { L S } =$ $\parallel A x _ { L S } - b \parallel _ { 2 } .$ .

Theorem 5.5.1. Suppose $U ^ { T } A V = \Sigma$ is the SVD of $A \in \mathbb { R } ^ { m \times n }$ with $r = r a n k ( A )$ . If $U = [  u _ { 1 } | \cdot \cdot \cdot | u _ { m } ]$ and $V = [  v _ { 1 } | \cdot \cdot \cdot | v _ { n } ]$ are column partitionings and $b \in \mathbb { R } ^ { m }$ , then

$$
x _ {L S} = \sum_ {i = 1} ^ {r} \frac {u _ {i} ^ {T} b}{\sigma_ {i}} v _ {i} \tag {5.5.1}
$$

minimizes $\parallel A x - b \parallel _ { 2 }$ and has the smallest 2-norm of all minimizers. Moreover

$$
\rho_ {L S} ^ {2} = \| A x _ {L S} - b \| _ {2} ^ {2} = \sum_ {i = r + 1} ^ {m} (u _ {i} ^ {T} b) ^ {2}. \tag {5.5.2}
$$

Proof. For any $\boldsymbol { x } \in \mathbb { R } ^ { n }$ we have

$$
\begin{array}{l} \| A x - b \| _ {2} ^ {2} = \| (U ^ {T} A V) (V ^ {T} x) - U ^ {T} b \| _ {2} ^ {2} = \| \Sigma \alpha - U ^ {T} b \| _ {2} ^ {2} \\ = \sum_ {i = 1} ^ {r} (\sigma_ {i} \alpha_ {i} - u _ {i} ^ {T} b) ^ {2} + \sum_ {i = r + 1} ^ {m} (u _ {i} ^ {T} b) ^ {2}, \\ \end{array}
$$

where $\alpha = V ^ { T } x$ . Clearly, if x solves the LS problem, then $\alpha _ { i } = ( u _ { i } ^ { T } b / \sigma _ { i } )$ for $i = 1 { : } r$ . If we set $\alpha ( r + 1 { : } n ) = 0$ , then the resulting x has minimal 2-norm.

# 5.5.2 A Note on the Pseudoinverse

If we define the matrix $A ^ { + } \in \mathbb { R } ^ { n \times m }$ by $A ^ { + } = V \Sigma ^ { + } U ^ { T }$ where

$$
\Sigma^ {+} = \mathrm{diag} \left(\frac {1}{\sigma_ {1}}, \ldots , \frac {1}{\sigma_ {r}}, 0, \ldots , 0\right) \in \mathbb {R} ^ {n \times m}, \qquad r = \mathsf {r a n k} (A),
$$

then $x _ { L S } = A ^ { + } b$ and $\rho _ { L S } = \parallel ( I - A A ^ { + } ) b \parallel _ { 2 }$ . A+ is referred to as the pseudo-inverse of A. It is the unique minimal Frobenius norm solution to the problem

$$
\min _ {X \in \mathbb {R} ^ {m \times n}} \| A X - I _ {m} \| _ {F}. \tag {5.5.3}
$$

If rank(A) = n, then $A ^ { + } = ( A ^ { T } A ) ^ { - 1 } A ^ { T }$ , while if $m = n = \operatorname { r a n k } ( A )$ , then $A ^ { + } = A ^ { - 1 }$ . Typically, A+ is defined to be the unique matrix $\boldsymbol { X } \in \mathbb { R } ^ { n \times m }$ that satisfies the four Moore-Penrose conditions:

These conditions amount to the requirement that $A A ^ { + }$ and $A ^ { + } A$ be orthogonal projections onto ran(A) and $\mathsf { r a n } ( A ^ { T } )$ , respectively. Indeed,

$$
A A ^ {+} = U _ {1} U _ {1} ^ {T}
$$

where $U _ { 1 } = U ( 1 { : } m , 1 { : } r )$ and

$$
A ^ {+} A = V _ {1} V _ {1} ^ {T}
$$

where $V _ { 1 } = V ( 1 { : } n , 1 { : } r )$ .

# 5.5.3 Some Sensitivity Issues

In §5.3 we examined the sensitivity of the full-rank LS problem. The behavior of $x _ { L S }$ in this situation is summarized in Theorem 5.3.1. If we drop the full-rank assumption, then $x _ { L S }$ is not even a continuous function of the data and small changes in A and b can induce arbitrarily large changes in $x _ { L S } = A ^ { + } b$ . The easiest way to see this is to consider the behavior of the pseudoinverse. If A and δA are in $\mathbb { R } ^ { m \times n }$ , then Wedin (1973) and Stewart (1975) show that

$$
\left\| (A + \delta A) ^ {+} - A ^ {+} \right\| _ {F} \leq 2 \| \delta A \| _ {F} \max \left\{\left\| A ^ {+} \right\| _ {2} ^ {2}, \left\| (A + \delta A) ^ {+} \right\| _ {2} ^ {2} \right\}.
$$

This inequality is a generalization of Theorem 2.3.4 in which perturbations in the matrix inverse are bounded. However, unlike the square nonsingular case, the upper bound does not necessarily tend to zero as δA tends to zero. If

$$
A = \left[ \begin{array}{l l} 1 & 0 \\ 0 & 0 \\ 0 & 0 \end{array} \right] \qquad \text {and} \qquad \delta A = \left[ \begin{array}{l l} 0 & 0 \\ 0 & \epsilon \\ 0 & 0 \end{array} \right]
$$

then

$$
A ^ {+} = \left[ \begin{array}{l l l} 1 & 0 & 0 \\ 0 & 0 & 0 \end{array} \right] \qquad \text {and} \qquad (A + \delta A) ^ {+} = \left[ \begin{array}{l l l} 1 & 0 & 0 \\ 1 & 1 / \epsilon & 0 \end{array} \right],
$$

and

$$
\| A ^ {+} - (A + \delta A) ^ {+} \| _ {2} = 1 / \epsilon .
$$

The numerical determination of an LS minimizer in the presence of such discontinuities is a major challenge.

# 5.5.4 The Truncated SVD Solution

Suppose ${ \widehat { U } } , { \widehat { \Sigma } } ,$ and $\widehat { V }$ are the computed SVD factors of a matrix A and $\hat { r }$ is accepted as its δ-rank, i.e.,

$$
\hat {\sigma} _ {n} \leq \dots \leq \hat {\sigma} _ {\hat {r}} \leq \delta <   \hat {\sigma} _ {\hat {r}} \leq \dots \leq \hat {\sigma} _ {1}.
$$

It follows that we can regard

$$
x _ {\hat {r}} = \sum_ {i = 1} ^ {\hat {r}} \frac {\hat {u} _ {i} ^ {T} b}{\hat {\sigma} _ {i}} \hat {v} _ {i}
$$

as an approximation to $x _ { L S }$ . Since $\Vert \ x _ { \hat { r } } \ \Vert _ { 2 } \approx 1 / \sigma _ { \hat { r } } \leq 1 / \delta$ , then $\delta$ may also be chosen with the intention of producing an approximate LS solution with suitably small norm. In $\ S 6 . 2 . 1$ , we discuss more sophisticated methods for doing this.

If $\hat { \sigma } _ { \hat { r } } \gg \delta$ , then we have reason to be comfortable with $x _ { \hat { r } }$ because A can then be unambiguously regarded as a $\mathrm { r a n k } ( A _ { \hat { r } } )$ matrix (modulo δ).

On the other hand, $\{ \hat { \sigma } _ { 1 } , \hdots , \hat { \sigma } _ { n } \}$ might not clearly split into subsets of small and large singular values, making the determination of $\hat { r }$ by this means somewhat arbitrary. This leads to more complicated methods for estimating rank, which we now discuss in the context of the LS problem. The issues are readily communicated by making two simplifying assumptions. Assume that $r ~ = ~ n$ , and that $\Delta A = 0$ in (5.4.4), which implies that $\boldsymbol { W } ^ { T } \boldsymbol { A } \boldsymbol { Z } = \widehat { \boldsymbol { \Sigma } } = \boldsymbol { \Sigma }$ is the SVD. Denote the ith columns of the matrices $\widehat { U }$ , $W , \widehat { V }$ , and $Z$ by $\hat { u } _ { i } , w _ { i } , \hat { v } _ { i }$ , and $z _ { i }$ , respectively. Because

$$
\begin{array}{l} x _ {L S} - x _ {\hat {r}} = \sum_ {i = 1} ^ {n} \frac {w _ {i} ^ {T} b}{\sigma_ {i}} z _ {i} - \sum_ {i = 1} ^ {\hat {r}} \frac {\hat {u} _ {i} ^ {T} b}{\sigma_ {i}} \hat {v} _ {i} \\ = \sum_ {i = 1} ^ {\hat {r}} \frac {((w _ {i} - \hat {u} _ {i}) ^ {T} b) z _ {i} + (\hat {u} _ {i} ^ {T} b) (z _ {i} - \hat {v} _ {i})}{\sigma_ {i}} + \sum_ {i = \hat {r} + 1} ^ {n} \frac {w _ {i} ^ {T} b}{\sigma_ {i}} z _ {i} \\ \end{array}
$$

it follows from $\begin{array} { r } { \| w _ { i } - \hat { u } _ { i } \| _ { 2 } \leq \epsilon , \| \hat { u } _ { i } \| _ { 2 } \leq 1 + \epsilon , } \end{array}$ and $\| z _ { i } - \hat { v } _ { i } \| _ { 2 } \leq \epsilon$ that

$$
\parallel x _ {\hat {r}} - x _ {L S} \parallel_ {2} \leq \frac {\hat {r}}{\sigma_ {\hat {r}}} 2 (1 + \epsilon) \epsilon \parallel b \parallel_ {2} + \sqrt {\sum_ {i = \hat {r} + 1} ^ {n} \left(\frac {w _ {i} ^ {T} b}{\sigma_ {i}}\right) ^ {2}}.
$$

The parameter $\hat { r }$ can be determined as that integer which minimizes the upper bound. Notice that the first term in the bound increases with ${ \hat { r } } ,$ , while the second decreases.

On occasions when minimizing the residual is more important than accuracy in the solution, we can determine $\hat { r }$ on the basis of how close we surmise $\parallel b - A x _ { \hat { r } } \parallel _ { 2 }$ is to the true minimum. Paralleling the above analysis, it can be shown that

$$
\| b - A x _ {\hat {r}} \| _ {2} \leq \| b - A x _ {L S} \| _ {2} + (n - \hat {r}) \| b \| _ {2} + \epsilon \hat {r} \| b \| _ {2} \left(1 + (1 + \epsilon) \frac {\hat {\sigma} _ {1}}{\hat {\sigma} _ {\hat {r}}}\right).
$$

Again $\hat { r }$ could be chosen to minimize the upper bound. See Varah (1973) for practical details and also LAPACK.

# 5.5.5 Basic Solutions via QR with Column Pivoting

Suppose $A \in \mathbb { R } ^ { m \times n }$ has rank r. QR with column pivoting (Algorithm 5.4.1) produces the factorization $A \Pi = Q R$ where

$$
R = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & 0 \\ r & n - r \end{array} \right] _ {m - r} ^ {r}.
$$

Given this reduction, the LS problem can be readily solved. Indeed, for any $\boldsymbol { x } \in \mathbb { R } ^ { n }$ we have

$$
\parallel A x - b \parallel_ {2} ^ {2} = \parallel (Q ^ {T} A \Pi) (\Pi^ {T} x) - (Q ^ {T} b) \parallel_ {2} ^ {2} = \parallel R _ {1 1} y - (c - R _ {1 2} z) \parallel_ {2} ^ {2} + \parallel d \parallel_ {2} ^ {2},
$$

where

$$
\Pi^ {T} x = \left[ \begin{array}{l} y \\ z \end{array} \right] _ {n - r} ^ {r} \quad \text {and} \quad Q ^ {T} b = \left[ \begin{array}{l} c \\ d \end{array} \right] _ {m - r} ^ {r}.
$$

Thus, if x is an LS minimizer, then we must have

$$
x = \Pi \left[ \begin{array}{c} R _ {1 1} ^ {- 1} (c - R _ {1 2} z) \\ z \end{array} \right].
$$

$\operatorname { I f } z$ is set to zero in this expression, then we obtain the basic solution

$$
x _ {B} = \Pi \left[ \begin{array}{c} R _ {1 1} ^ {- 1} c \\ 0 \end{array} \right].
$$

Notice that $x _ { B }$ has at most r nonzero components and so $A x _ { B }$ involves a subset of A’s columns.

The basic solution is not the minimal 2-norm solution unless the submatrix $R _ { 1 2 }$ is zero since

$$
\left\| x _ {L S} \right\| _ {2} = \min _ {z \in \mathbb {R} ^ {n - 2}} \left\| x _ {B} - \Pi \left[ \begin{array}{c} R _ {1 1} ^ {- 1} R _ {1 2} \\ - I _ {n - r} \end{array} \right] z \right\| _ {2}. \tag {5.5.4}
$$

Indeed, this characterization of $\Vert { x } _ { L S } \ \Vert _ { 2 }$ can be used to show that

$$
1 \leq \frac {\| x _ {B} \| _ {2}}{\| x _ {L S} \| _ {2}} \leq \sqrt {1 + \| R _ {1 1} ^ {- 1} R _ {1 2} \| _ {2} ^ {2}}. \tag {5.5.5}
$$

See Golub and Pereyra (1976) for details.

# 5.5.6 Some Comparisons

As we mentioned, when solving the LS problem via the SVD, only Σ and V have to be computed assuming that the right hand side b is available. The table in Figure 5.5.1 compares the flop efficiency of this approach with the other algorithms that we have presented.

<table><tr><td>LS Algorithm</td><td>Flop Count</td></tr><tr><td>Normal equations</td><td> $mn^{2} + n^{3}/3$ </td></tr><tr><td>Householder QR</td><td> $n^{3}/3$ </td></tr><tr><td>Modified Gram-Schmidt</td><td> $2mn^{2}$ </td></tr><tr><td>Givens QR</td><td> $3mn^{2} - n^{3}$ </td></tr><tr><td>Householder Bidiagonalization</td><td> $4mn^{2} - 2n^{3}$ </td></tr><tr><td>R-Bidiagonalization</td><td> $2mn^{2} + 2n^{3}$ </td></tr><tr><td>SVD</td><td> $4mn^{2} + 8n^{3}$ </td></tr><tr><td>R-SVD</td><td> $2mn^{2} + 11n^{3}$ </td></tr></table>

Figure 5.5.1. Flops associated with various least squares methods

# 5.5.7 SVD-Based Subset Selection

Replacing A by $A _ { \tilde { r } }$ in the LS problem amounts to filtering the small singular values and can make a great deal of sense in those situations where A is derived from noisy data. In other applications, however, rank deficiency implies redundancy among the factors that comprise the underlying model. In this case, the model-builder may not be interested in a predictor such as $A _ { \tilde { r } } x _ { \tilde { r } }$ that involves all n redundant factors. Instead, a predictor $A y$ may be sought where y has at most ˜r nonzero components. The position of the nonzero entries determines which columns of A, i.e., which factors in the model, are to be used in approximating the observation vector b. How to pick these columns is the problem of subset selection.

QR with column pivoting is one way to proceed. However, Golub, Klema, and Stewart (1976) have suggested a technique that heuristically identifies a more independent set of columns than are involved in the predictor $A x _ { B }$ . The method involves both the SVD and QR with column pivoting:

Step 1. Compute the SVD $A = U \Sigma V ^ { T }$ and use it to determine a rank estimate ˜r.

Step 2. Calculate a permutation matrix P such that the columns of the matrix $B _ { 1 } \in \mathbb { R } ^ { m \times \tilde { r } }$ in $A P = \left[ B _ { 1 } \mid B _ { 2 } \right]$ are “sufficiently independent.”

Step 3. Predict b with Ay where $\boldsymbol { y } = \boldsymbol { P } \left[ \begin{array} { l } { z } \\ { 0 } \end{array} \right]$ z and $z \in \mathbb { R } ^ { \tilde { r } }$ minimizes $\parallel B _ { 1 } z - b \parallel _ { 2 }$ .

The second step is key. Because

$$
\min _ {z \in \mathbf {R} ^ {\tilde {r}}} \| B _ {1} z - b \| _ {2} = \| A y - b \| _ {2} \geq \min _ {x \in \mathbf {R} ^ {n}} \| A x - b \| _ {2}
$$

it can be argued that the permutation P should be chosen to make the residual $r =$ $( I - B _ { 1 } B _ { 1 } ^ { + } ) b$ as small as possible. Unfortunately, such a solution procedure can be

unstable. For example, if

$$
A = \left[ \begin{array}{c c c} 1 & 1 & 0 \\ 1 & 1 + \epsilon & 1 \\ 0 & 0 & 1 \end{array} \right], \qquad b = \left[ \begin{array}{c} 1 \\ - 1 \\ 0 \end{array} \right],
$$

$\tilde { r } = 2$ , and $P = I$ , then min $\| \ B _ { 1 } z - b \| _ { 2 } = 0$ , but $\parallel B _ { 1 } ^ { + } b \parallel _ { 2 } = { \cal O } ( 1 / \epsilon )$ . On the other hand, any proper subset involving the third column of A is strongly independent but renders a much larger residual.

This example shows that there can be a trade-off between the independence of the chosen columns and the norm of the residual that they render. How to proceed in the face of this trade-off requires useful bounds on $\sigma _ { \tilde { r } } ( B _ { 1 } )$ , the smallest singular value of $B _ { 1 }$ .

Theorem 5.5.2. Let the SVD of $A \in \mathbb { R } ^ { m \times n }$ be given by $U ^ { T } A V = \Sigma = \mathrm { d i a g } ( \sigma _ { i } )$ and define the matrix $B _ { 1 } \in \mathbb { R } ^ { m \times \tilde { r } } , \tilde { r } \leq \mathsf { r a n k } ( A )$ , by

$$
A P = \left[ \begin{array}{c c} B _ {1} & B _ {2} \\ \tilde {r} & n - \tilde {r} \end{array} \right]
$$

where $P \in \mathbb { R } ^ { n \times n }$ is a permutation. If

$$
P ^ {T} V = \left[ \begin{array}{c c} \widetilde {V} _ {1 1} & \widetilde {V} _ {1 2} \\ \widetilde {V} _ {2 1} & \widetilde {V} _ {2 2} \end{array} \right] _ {n - \tilde {r}} ^ {\tilde {r}} \tag {5.5.6}
$$

and $\tilde { V } _ { 1 1 }$ is nonsingular, then

$$
\frac {\sigma_ {\tilde {r}} (A)}{\| \tilde {V} _ {1 1} ^ {- 1} \| _ {2}} \leq \sigma_ {\tilde {r}} (B _ {1}) \leq \sigma_ {\tilde {r}} (A).
$$

Proof. The upper bound follows from Corollary 2.4.4. To establish the lower bound, partition the diagonal matrix of singular values as follows:

$$
\Sigma = \left[ \begin{array}{c c} \Sigma_ {1} & 0 \\ 0 & \Sigma_ {2} \end{array} \right] _ {m - \tilde {r}} ^ {\tilde {r}}.
$$

If $\boldsymbol { w } \in \mathbb { R } ^ { \tilde { r } }$ is a unit vector with the property that $\parallel B _ { 1 } w \parallel _ { 2 } = \sigma _ { \tilde { r } } ( B _ { 1 } )$ , then

$$
\sigma_ {\tilde {r}} (B _ {1}) ^ {2} = \parallel B _ {1} w \parallel_ {2} ^ {2} = \left\| U \Sigma V ^ {T} P \left[ \begin{array}{c} w \\ 0 \end{array} \right] \right\| _ {2} ^ {2} = \parallel \Sigma_ {1} \widetilde {V} _ {1 1} ^ {T} w \parallel_ {2} ^ {2} + \parallel \Sigma_ {2} \widetilde {V} _ {1 2} ^ {T} w \parallel_ {2} ^ {2}.
$$

The theorem now follows because $\begin{array} { r } { \| \Sigma _ { 1 } \widetilde { V } _ { 1 1 } ^ { T } w \| _ { 2 } \geq \sigma _ { \widetilde { r } } ( A ) / \| \widetilde { V } _ { 1 1 } ^ { - 1 } \| _ { 2 } . } \end{array}$

This result suggests that in the interest of obtaining a sufficiently independent subset of columns, we choose the permutation P such that the resulting $\widetilde { V } _ { 1 1 }$ submatrix is as well-conditioned as possible. A heuristic solution to this problem can be obtained by computing the QR with column-pivoting factorization of the matrix $\left[ V _ { 1 1 } ^ { T } ~ V _ { 2 1 } ^ { T } \right]$ , where

$$
V = \left[ \begin{array}{c c} V _ {1 1} & V _ {1 2} \\ V _ {2 1} & V _ {2 2} \end{array} \right] _ {n - \tilde {r}} ^ {\tilde {r}}
$$

is a partitioning of the matrix V , A’s matrix of right singular vectors. In particular, if we apply QR with column pivoting (Algorithm 5.4.1) to compute

$$
Q ^ {T} [ V _ {1 1} ^ {T} V _ {2 1} ^ {T} ] P = [ R _ {1 1} | R _ {1 2} ] _ {\tilde {r} n - \tilde {r}}
$$

where $Q$ is orthogonal, P is a permutation matrix, and $R _ { 1 1 }$ is upper triangular, then (5.5.6) implies

$$
\left[ \begin{array}{l} \widetilde {V} _ {1 1} \\ \widetilde {V} _ {2 1} \end{array} \right] = P ^ {T} \left[ \begin{array}{l} V _ {1 1} \\ V _ {2 1} \end{array} \right] = \left[ \begin{array}{l} R _ {1 1} ^ {T} Q ^ {T} \\ R _ {1 2} ^ {T} Q ^ {T} \end{array} \right].
$$

Note that $R _ { 1 1 }$ is nonsingular and that $\parallel \widetilde { V } _ { 1 1 } ^ { - 1 } \parallel _ { 2 } = \parallel R _ { 1 1 } ^ { - 1 } \parallel _ { 2 }$ . Heuristically, column pivoting tends to produce a well-conditioned $R _ { 1 1 }$ , and so the overall process tends to produce a well-conditioned $\widetilde { V } _ { 1 1 }$ .

Algorithm 5.5.1 Given $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ the following algorithm computes a permutation $P ,$ , a rank estimate ${ \tilde { r } } _ { : }$ and a vector $z \in \mathbb { R } ^ { \tilde { r } }$ such that the first $\tilde { r }$ columns of $B = A P$ are independent and $\parallel B ( : , 1 : \tilde { r } ) z - b \parallel _ { 2 }$ is minimized.

Compute the SVD $U ^ { T } A V = \mathrm { d i a g } ( \sigma _ { 1 } , . . . , \sigma _ { n } )$ and save $V .$ .

Determine $\tilde { r } \le \mathsf { r a n k } ( A )$ .

Apply QR with column pivoting: $Q ^ { T } V ( : , 1 : \tilde { r } ) ^ { T } P = [ R _ { 1 1 } | R _ { 1 2 } ]$ and set

$$
A P = \left[ B _ {1} \mid B _ {2} \right] \text {with} B _ {1} \in \mathbb {R} ^ {m \times \tilde {r}} \text {and} B _ {2} \in \mathbb {R} ^ {m \times (n - \tilde {r})}.
$$

Determine $z \in \mathbb { R } ^ { \tilde { r } }$ such that $\parallel b - B _ { 1 } z \parallel _ { 2 } = \operatorname* { m i n }$ .

# 5.5.8 Column Independence Versus Residual Size

We return to the discussion of the trade-off between column independence and norm of the residual. In particular, to assess the above method of subset selection we need to examine the residual of the vector y that it produces

$$
r _ {y} = b - A y = b - B _ {1} z = (I - B _ {1} B _ {1} ^ {+}) b.
$$

Here, $B _ { 1 } = B ( : , 1 { : } \tilde { r } )$ with $B = A P$ . To this end, it is appropriate to compare $r _ { y }$ with

$$
r _ {x _ {\tilde {r}}} = b - A x _ {\tilde {r}}
$$

since we are regarding A as a rank-˜r matrix and since $x _ { \tilde { r } }$ solves the nearest rank-˜r LS problem min $\parallel A _ { \tilde { r } } x - b \parallel _ { 2 }$ .

Theorem 5.5.3. Assume that $U ^ { T } A V = \Sigma$ is the SVD of $A \in \mathbb { R } ^ { m \times n }$ and that $r _ { y }$ and $r _ { x _ { \tilde { r } } }$ are defined as above. $I f \tilde { V } _ { 1 1 }$ is the leading r-by-r principal submatrix of $P ^ { T } V$ , then

$$
\| r _ {x _ {\tilde {r}}} - r _ {y} \| _ {2} \leq \frac {\sigma_ {\tilde {r} + 1} (A)}{\sigma_ {\tilde {r}} (A)} \| \tilde {V} _ {1 1} ^ {- 1} \| _ {2} \| b \| _ {2}.
$$

Proof. Note that $r _ { x _ { \tilde { r } } } = ( I - U _ { 1 } U _ { 1 } ^ { T } ) l$ and $r _ { y } = ( I - Q _ { 1 } Q _ { 1 } ^ { T } ) b$ where

$$
U = \left[ \begin{array}{c c} U _ {1} & U _ {2} \\ \tilde {r} & m - \tilde {r} \end{array} \right]
$$

is a partitioning of the matrix U and $Q _ { 1 } = B _ { 1 } ( B _ { 1 } ^ { T } B _ { 1 } ) ^ { - 1 / 2 }$ . Using Theorem 2.6.1 we obtain

$$
\left\| r _ {x _ {\tilde {r}}} - r _ {y} \right\| _ {2} \leq \left\| U _ {1} U _ {1} ^ {T} - Q _ {1} Q _ {1} ^ {T} \right\| _ {2} \left\| b \right\| _ {2} = \left\| U _ {2} ^ {T} Q _ {1} \right\| _ {2} \left\| b \right\| _ {2}
$$

while Theorem 5.5.2 permits us to conclude

$$
\begin{array}{l} \| U _ {2} ^ {T} Q _ {1} \| _ {2} \leq \| U _ {2} ^ {T} B _ {1} \| _ {2} \| (B _ {1} ^ {T} B _ {1}) ^ {- 1 / 2} \| _ {2} \\ \leq \sigma_ {\tilde {r} + 1} (A) \frac {1}{\sigma_ {\tilde {r}} (B _ {1})} \leq \frac {\sigma_ {\tilde {r} + 1} (A)}{\sigma_ {\tilde {r}} (A)} \| \tilde {V} _ {1 1} ^ {- 1} \| _ {2}, \\ \end{array}
$$

and this establishes the theorem.

Noting that

$$
\left\| r _ {x _ {\bar {r}}} - r _ {y} \right\| _ {2} = \left\| B _ {1} y - \sum_ {i = 1} ^ {r} (u _ {i} ^ {T} b) u _ {i} \right\| _ {2}
$$

we see that Theorem 5.5.3 sheds light on how well $B _ { 1 } y$ can predict the “stable” component of b, i.e., $U _ { 1 } ^ { T } b$ . Any attempt to approximate $U _ { 2 } ^ { T } b$ can lead to a large norm solution. Moreover, the theorem says that if $\sigma _ { \tilde { r } + 1 } ( A ) \ll \sigma _ { \tilde { r } } ( A )$ , then any reasonably independent subset of columns produces essentially the same-sized residual. On the other hand, if there is no well-defined gap in the singular values, then the determination of $\tilde { r }$ becomes difficult and the entire subset selection problem becomes more complicated.

# Problems

P5.5.1 Show that if

$$
A = \left[ \begin{array}{c c} T & S \\ 0 & 0 \end{array} \right] _ {m - r} ^ {r}
$$

where $r = { \mathsf { r a n k } } ( A )$ and T is nonsingular, then

$$
X = \left[ \begin{array}{c c} T ^ {- 1} & 0 \\ 0 & 0 \\ r & m - r \end{array} \right] _ {n - r} ^ {r}
$$

satisfies $A X A = A$ and $( A X ) ^ { T } = ( A X )$ . In this case, we say that X is a (1,3) pseudoinverse of A. Show that for general A, $x _ { B } = X b$ where X is a (1,3) pseudoinverse of A.

P5.5.2 Define $B ( \lambda ) \in \mathbb { R } ^ { n \times m }$ by

$$
B (\lambda) = (A ^ {T} A + \lambda I) ^ {- 1} A ^ {T}
$$

where $\lambda > 0$ . Show that

$$
\parallel B (\lambda) - A ^ {+} \parallel_ {2} = \frac {\lambda}{\sigma_ {r} (A) [ \sigma_ {r} (A) ^ {2} + \lambda ]}, \qquad r = \operatorname{rank} (A),
$$

and therefore that $B ( \lambda ) \to A ^ { + }$ as $\lambda  0$ .

P5.5.3 Consider the rank-deficient LS problem

$$
\min _ {y \in \mathbf {R} ^ {r},   z \in \mathbf {R} ^ {n - r}} \left\| \left[ \begin{array}{c c} R & S \\ 0 & 0 \end{array} \right] \left[ \begin{array}{c} y \\ z \end{array} \right] - \left[ \begin{array}{c} c \\ d \end{array} \right] \right\| _ {2}
$$

where $R \in \mathbb R ^ { r \times r } , S \in \mathbb R ^ { r \times n - r } , y \in \mathbb R ^ { r }$ , and $z \in \mathbb { R } ^ { n - r }$ . Assume that R is upper triangular and nonsingular. Show how to obtain the minimum norm solution to this problem by computing an appropriate QR factorization without pivoting and then solving for the appropriate y and $z .$

P5.5.4 Show that if $A _ { k } \to A$ and $A _ { k } ^ { + } \to A ^ { + }$ , then there exists an integer k0 such that rank $\left( A _ { k } \right)$ is constant for all $k \geq k _ { 0 }$ .

P5.5.5 Show that if $A \in \mathbb { R } ^ { m \times n }$ has rank n, then so does $A + E { \mathrm { ~ i f ~ } } \| E \| _ { 2 } \| A ^ { + } \| _ { 2 } < 1$ .

P5.5.6 Suppose $A \in \mathbb { R } ^ { m \times n }$ is rank deficient and $b \in \mathbb { R } ^ { m }$ . Assume for $k = 0 , 1 , \ldots$ . that $x ^ { ( k + 1 ) }$ minimizes

$$
\phi_ {k} (x) = \parallel A x - b \parallel_ {2} ^ {2} + \lambda \parallel x - x ^ {(k)} \parallel_ {2} ^ {2}
$$

where $\lambda > 0$ and $x ^ { ( 0 ) } = 0$ . Show that $x ^ { ( k ) } \to x _ { L S }$

P5.5.8 Suppose $A \in \mathbb { R } ^ { m \times n }$ and that $\Vert u ^ { T } A \Vert _ { 2 } = \sigma$ with $u ^ { T } u = 1$ . Show that if $u ^ { T } ( A x - b ) = 0$ for $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $b \in \mathbb { R } ^ { m }$ , then $\parallel x \parallel _ { 2 } \geq | u ^ { T } b | / \sigma$ .

P5.5.9 In Equation (5.5.6) we know that the matrix $P ^ { T } V$ is orthogonal. Thus, $\parallel \tilde { V } _ { 1 1 } ^ { - 1 } \parallel _ { 2 } = \parallel \tilde { V } _ { 2 2 } ^ { - 1 } \parallel _ { 2 }$ from the CS decomposition (Theorem 2.5.3). Show how to compute P by applying the QR with column-pivoting algorithm to $[ \tilde { V } _ { 2 2 } ^ { T } | \tilde { V } _ { 1 2 } ^ { T } ]$ . (For $\tilde { r } > n / 2$ , this procedure would be more economical than the technique discussed in the text.) Incorporate this observation in Algorithm 5.5.1.

P5.5.10 Suppose $F \in \mathbb { R } ^ { m \times r }$ and $G \in \mathbb { R } ^ { n \times r }$ each have rank r. (a) Give an efficient algorithm for computing the minimum 2-norm minimizer of $\parallel F G ^ { T } x - b \parallel _ { 2 }$ where $b \in \mathbb { R } ^ { m }$ . (b) Show how to compute the vector $x _ { B }$ .

# Notes and References for §5.5

For a comprehensive treatment of the pseudoinverse and its manipulation, see:

M.Z. Nashed (1976). Generalized Inverses and Applications, Academic Press, New York.

S.L. Campbell and C.D. Meyer (2009). Generalized Inverses of Linear Transformations, SIAM Publications, Philadelphia, PA.

For an analysis of how the pseudo-inverse is affected by perturbation, see:

P.A. Wedin (1973). “Perturbation Theory for Pseudo-Inverses,” BIT 13, 217–232.

G.W. Stewart (1977). “On the Perturbation of Pseudo-Inverses, Projections, and Linear Least Squares,” SIAM Review 19, 634–662.

Even for full rank problems, column pivoting seems to produce more accurate solutions. The error analysis in the following paper attempts to explain why:

L.S. Jennings and M.R. Osborne (1974). “A Direct Error Analysis for Least Squares,” Numer. Math. 22, 322–332.

Various other aspects of the rank-deficient least squares problem are discussed in:

J.M. Varah (1973). “On the Numerical Solution of Ill-Conditioned Linear Systems with Applications to Ill-Posed Problems,” SIAM J. Numer. Anal. 10, 257–67.

G.W. Stewart (1984). “Rank Degeneracy,” SIAM J. Sci. Stat. Comput. 5, 403–413.

P.C. Hansen (1987). “The Truncated SVD as a Method for Regularization,” BIT 27, 534–553.

G.W. Stewart (1987). “Collinearity and Least Squares Regression,” Stat. Sci. 2, 68–100.

R.D. Fierro and P.C. Hansen (1995). “Accuracy of TSVD Solutions Computed from Rank-Revealing Decompositions,” Numer. Math. 70, 453–472.   
P.C. Hansen (1997). Rank-Deficient and Discrete Ill-Posed Problems: Numerical Aspects of Linear Inversion, SIAM Publications, Philadelphia, PA.   
A. Dax and L. Elden (1998). “Approximating Minimum Norm Solutions of Rank-Deficient Least Squares Problems,” Numer. Lin. Alg. 5, 79–99.   
G. Quintana-Orti, E.S. Quintana-Orti, and A. Petitet (1998). “Efficient Solution of the Rank-Deficient Linear Least Squares Problem,” SIAM J. Sci. Comput. 20, 1155–1163.   
L.V. Foster (2003). “Solving Rank-Deficient and Ill-posed Problems Using UTV and QR Factorizations,” SIAM J. Matrix Anal. Applic. 25, 582–600.   
D.A. Huckaby and T.F. Chan (2004). “Stewart’s Pivoted QLP Decomposition for Low-Rank Matrices,” Numer. Lin. Alg. 12, 153–159.   
L. Foster and R. Kommu (2006). “Algorithm 853: An Efficient Algorithm for Solving Rank-Deficient Least Squares Problems,” ACM Trans. Math. Softw. 32, 157–165.

For a sampling of the subset selection literature, we refer the reader to:

H. Hotelling (1957). “The Relations of the Newer Multivariate Statistical Methods to Factor Analysis,” Brit. J. Stat. Psych. 10, 69–79.

G.H. Golub, V. Klema and G.W. Stewart (1976). “Rank Degeneracy and Least Squares Problems,” Technical Report TR-456, Department of Computer Science, University of Maryland, College Park, MD.

S. Van Huffel and J. Vandewalle (1987). “Subset Selection Using the Total Least Squares Approach in Collinearity Problems with Errors in the Variables,” Lin. Alg. Applic. 88/89, 695–714.

M.R. Osborne, B. Presnell, and B.A. Turlach (2000). “A New Approach to Variable Selection in Least Squares Problems,” IMA J. Numer. Anal. 20, 389–403.
