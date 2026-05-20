# 8.6 Computing the SVD

If $U ^ { T } A V = B$ is the bidiagonal decomposition of $A \in \mathbb { R } ^ { m \times n }$ , then $V ^ { T } ( A ^ { T } A ) V = B ^ { T } B$ is the tridiagonal decomposition of the symmetric matrix $A ^ { T } A \in \mathbb { R } ^ { n \times { \dot { n } } }$ . Thus, there is an intimate connection between Algorithm 5.4.2 (Householder bidiagonalization) and Algorithm 8.3.1 (Householder tridiagonalization). In this section we carry this a step further and show that there is a bidiagonal SVD procedure that corresponds to the symmetric tridiagonal QR iteration. Before we get into the details, we catalog some important SVD properties that have algorithmic ramifications.

# 8.6.1 Connections to the Symmetric Eigenvalue Problem

There are important relationships between the singular value decomposition of a matrix A and the Schur decompositions of the symmetric matrices

$$
S _ {1} = A ^ {T} A, \qquad S _ {2} = A A ^ {T} \qquad S _ {3} = \left[ \begin{array}{c c} 0 & A ^ {T} \\ A & 0 \end{array} \right].
$$

Indeed, if

$$
U ^ {T} A V = \operatorname{diag} \left(\sigma_ {1}, \dots , \sigma_ {n}\right)
$$

is the SVD of $A \in \mathbb { R } ^ { m \times n } ~ ( m \geq n )$ , then

$$
V ^ {T} (A ^ {T} A) V = \operatorname{diag} (\sigma_ {1} ^ {2}, \dots , \sigma_ {n} ^ {2}) \in \mathbb {R} ^ {n \times n} \tag {8.6.1}
$$

and

$$
U ^ {T} (A A ^ {T}) U = \operatorname{diag} (\sigma_ {1} ^ {2}, \dots , \sigma_ {n} ^ {2}, \underbrace {0 , \dots , 0} _ {m - n}) \in \mathbb {R} ^ {m \times m} \tag {8.6.2}
$$

Moreover, if

$$
U = \left[ \begin{array}{c c} U _ {1} & U _ {2} \\ n & m - n \end{array} \right]
$$

and we define the orthogonal matrix $Q \in \mathbb { R } ^ { ( m + n ) \times ( m + n ) }$ by

$$
Q = \frac {1}{\sqrt {2}} \left[ \begin{array}{c c c} V & V & 0 \\ U _ {1} & - U _ {1} & \sqrt {2}   U _ {2} \end{array} \right],
$$

then

$$
Q ^ {T} \left[ \begin{array}{c c} 0 & A ^ {T} \\ A & 0 \end{array} \right] Q = \operatorname{diag} \left(\sigma_ {1}, \dots , \sigma_ {n}, - \sigma_ {1}, \dots , - \sigma_ {n}, \underbrace {0 , \dots , 0} _ {m - n}\right). \tag {8.6.3}
$$

These connections to the symmetric eigenproblem allow us to adapt the mathematical and algorithmic developments of the previous sections to the singular value problem. Good references for this section include Lawson and Hanson (SLS) and Stewart and Sun (MPT).

# 8.6.2 Perturbation Theory and Properties

We first establish perturbation results for the SVD based on the theorems of §8.1. Recall that $\sigma _ { i } ( A )$ denotes the ith largest singular value of A.

Theorem 8.6.1. If $A \in \mathbb { R } ^ { m \times n }$ , then for $k = 1 { : } \operatorname* { m i n } \{ m , n \}$

$$
\sigma_{k}(A) = \min_{\mathbf{dim}(S) = n - k + 1}\max_{\substack{x\in S\\ y\in \mathbf{R}^{m}}}\frac{y^{T}Ax}{\|x\|_{2}\|y\|_{2}} = \max_{\mathbf{dim}(S) = k}\min_{x\in S}\frac{\|Ax\|_{2}}{\|x\|_{2}}.
$$

In this expression, S is a subspace of $\mathbb { R } ^ { n }$ .

Proof. The rightmost characterization follows by applying Theorem 8.1.2 to $A ^ { T } A$ . For the remainder of the proof see Xiang (2006).

Corollary 8.6.2. If A and $A + E$ are in $\mathbb { R } ^ { m \times n }$ with $m \geq n$ , then for $k = 1 { : } n$

$$
\left| \sigma_ {k} (A + E) - \sigma_ {k} (A) \right| \leq \sigma_ {1} (E) = \| E \| _ {2}.
$$

Proof. Define $\widetilde { A }$ and $\widetilde { E }$ by

$$
\widetilde {A} = \left[ \begin{array}{c c} 0 & A ^ {T} \\ A & 0 \end{array} \right], \quad \widetilde {A} + \widetilde {E} = \left[ \begin{array}{c c} 0 & (A + E) ^ {T} \\ A + E & 0 \end{array} \right]. \tag {8.6.4}
$$

The corollary follows by applying Corollary 8.1.6 with A replaced by $\widetilde { A }$ and $A + E$ replaced by $\widetilde A + \widetilde E$ .

Corollary 8.6.3. Let $A = [ a _ { 1 } \mid \cdots \mid a _ { n } ] \in \mathbb { R } ^ { m \times n }$ be a column partitioning with $m \geq$ n. $I f A _ { r } = \left[ a _ { 1 } \ : | \cdots | a _ { r } \ : \right]$ , then for $r = 1 { : } n - 1$

$$
\sigma_ {1} (A _ {r + 1}) \geq \sigma_ {1} (A _ {r}) \geq \sigma_ {2} (A _ {r + 1}) \geq \dots \geq \sigma_ {r} (A _ {r + 1}) \geq \sigma_ {r} (A _ {r}) \geq \sigma_ {r + 1} (A _ {r + 1}).
$$

Proof. Apply Corollary 8.1.7 to $A ^ { T } A$ .

The next result is a Wielandt-Hoffman theorem for singular values:

Theorem 8.6.4. If A and $A + E$ are in $\mathbb { R } ^ { m \times n }$ with $m \geq n$ , then

$$
\sum_ {k = 1} ^ {n} \left(\sigma_ {k} (A + E) - \sigma_ {k} (A)\right) ^ {2} \leq \| E \| _ {F} ^ {2}.
$$

Proof. Apply Theorem 8.1.4 with A and E replaced by the matrices $\tilde { A }$ and $\tilde { E }$ defined by (8.6.4).

For $A \in \mathbb { R } ^ { m \times n }$ we say that the k-dimensional subspaces $S \subseteq \mathbb { R } ^ { n }$ and $T \subseteq \mathbb { R } ^ { m }$ form a singular subspace pair if $x \in S$ and $y \in T$ imply Ax $\in T$ and $A ^ { T } y \in S$ . The following result is concerned with the perturbation of singular subspace pairs.

Theorem 8.6.5. Let A, $E \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ be given and suppose that $V \in \mathbb { R } ^ { n \times n }$ and $U \in \mathbb { R } ^ { m \times m }$ are orthogonal. Assume that

$$
V = \left[ \begin{array}{c c} V _ {1} & V _ {2} \\ r & n - r \end{array} \right]    , \qquad U = \left[ \begin{array}{c c} U _ {1} & U _ {2} \\ r & m - r \end{array} \right]    ,
$$

and that ran(V1) and ran $( U _ { 1 } )$ form a singular subspace pair for A. Let

$$
U ^ {T} A V = \left[ \begin{array}{c c} A _ {1 1} & 0 \\ 0 & A _ {2 2} \end{array} \right] _ {m - r} ^ {r}, \qquad U ^ {T} E V = \left[ \begin{array}{c c} E _ {1 1} & E _ {1 2} \\ E _ {2 1} & E _ {2 2} \end{array} \right] _ {m - r} ^ {r},
$$

and assume that

$$
\delta = \min_{\substack{\sigma \in \sigma (A_{11})\\ \gamma \in \sigma (A_{22})}}|\sigma -\gamma | > 0.
$$

If

$$
\left\| E \right\| _ {F} \leq \frac {\delta}{5},
$$

then there exist matrices $P \in \mathbb { R } ^ { ( n - r ) \times r }$ and $Q \in \mathbb { R } ^ { ( m - r ) \times r }$ satisfying

$$
\left\| \left[ \begin{array}{c} Q \\ P \end{array} \right] \right\| _ {F} \leq 4 \frac {\left\| E \right\| _ {F}}{\delta}
$$

such that ran $( V _ { 1 } + V _ { 2 } Q )$ and ran $( U _ { 1 } + U _ { 2 } P )$ is a singular subspace pair for $A + E$ .

Proof. See Stewart (1973, Theorem 6.4).

Roughly speaking, the theorem says that $O ( \epsilon )$ changes in A can alter a singular subspace by an amount $\epsilon / \delta$ where δ measures the separation of the associated singular values.

# 8.6.3 The SVD Algorithm

We now show how a variant of the QR algorithm can be used to compute the SVD of an $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ . At first glance, this appears straightforward. Equation (8.6.1) suggests that we proceed as follows:

Step 1. Form $C = A ^ { T } A$ ,

Step 2. Use the symmetric QR algorithm to compute $V _ { 1 } ^ { T } C V _ { 1 } = \mathrm { { d i a g } } ( \sigma _ { i } ^ { 2 } )$

Step 3. Apply QR with column pivoting to $A V _ { 1 }$ obtaining $U ^ { T } ( A V _ { 1 } ) \Pi = R .$

Since R has orthogonal columns, it follows that $U ^ { T } A ( V _ { 1 } \Pi )$ is diagonal. However, as we saw in §5.3.2, the formation of $A ^ { T } A$ can lead to a loss of information. The situation is not quite so bad here, since the original A is used to compute U .

A preferable method for computing the SVD is described by Golub and Kahan (1965). Their technique finds U and V simultaneously by implicitly applying the symmetric QR algorithm to $A ^ { T } A$ . The first step is to reduce A to upper bidiagonal form using Algorithm 5.4.2:

$$
U _ {B} ^ {T} A V _ {B} = \left[ \begin{array}{l} B \\ 0 \end{array} \right], \qquad B = \left[ \begin{array}{l l l l l} d _ {1} & f _ {1} & & \dots & 0 \\ 0 & d _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & f _ {n - 1} \\ 0 & \dots & & 0 & d _ {n} \end{array} \right] \in \mathbb {R} ^ {n \times n}.
$$

The remaining problem is thus to compute the SVD of B. To this end, consider applying an implicit-shift QR step (Algorithm 8.3.2) to the tridiagonal matrix $T = B ^ { T } B$ :

Step 1. Compute the eigenvalue λ of

$$
T (m: n, m: n) = \left[ \begin{array}{c c} d _ {m} ^ {2} + f _ {m - 1} ^ {2} & d _ {m} f _ {m} \\ d _ {m} f _ {m} & d _ {n} ^ {2} + f _ {m} ^ {2} \end{array} \right], \qquad m = n - 1,
$$

that is closer to $d _ { n } ^ { 2 } + f _ { m } ^ { 2 }$

Step 2. Compute $c _ { 1 } = \cos ( \theta _ { 1 } )$ and $s _ { 1 } = \sin ( \theta _ { 1 } )$ such that

$$
\left[ \begin{array}{c c} {c _ {1}} & {s _ {1}} \\ {- s _ {1}} & {c _ {1}} \end{array} \right] ^ {T} \left[ \begin{array}{c} {d _ {1} ^ {2} - \lambda} \\ {d _ {1} f _ {1}} \end{array} \right] = \left[ \begin{array}{c} {\times} \\ {0} \end{array} \right]
$$

and set $G _ { 1 } = G ( 1 , 2 , \theta _ { 1 } )$ .

Step 3. Compute Givens rotations $G _ { 2 } , \ldots , G _ { n - 1 }$ so that if $Q = G _ { 1 } \cdot \cdot \cdot G _ { n - 1 }$ then $Q ^ { T } T Q$ is tridiagonal and $Q e _ { 1 } = G _ { 1 } e _ { 1 }$ .

Note that these calculations require the explicit formation of $B ^ { T } B$ , which, as we have seen, is unwise from the numerical standpoint.

Suppose instead that we apply the Givens rotation $G _ { 1 }$ above to B directly. Illustrating with the $n = 6$ case we have

$$
B \leftarrow B G _ {1} = \left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ + & \times & \times & 0 & 0 & 0 \\ 0 & 0 & \times & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

We then can determine Givens rotations $U _ { 1 } , V _ { 2 } , U _ { 2 } , . . . , V _ { n - 1 }$ , and $U _ { n - 1 }$ to chase the unwanted nonzero element down the bidiagonal:

$$
B \gets U _ {1} ^ {T} B = \left[ \begin{array}{l l l l l l} \times & \times & + & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 \\ 0 & 0 & \times & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right], \quad B \gets B V _ {2} = \left[ \begin{array}{l l l l l l} \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 \\ 0 & + & \times & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
B \gets U _ {2} ^ {T} B = \left[ \begin{array}{l l l l l l} \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & + & 0 & 0 \\ 0 & 0 & \times & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right], \quad B \gets B V _ {3} = \left[ \begin{array}{l l l l l l} \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 \\ 0 & 0 & \times & \times & 0 & 0 \\ 0 & 0 & + & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

and so on. The process terminates with a new bidiagonal $\tilde { B }$ that is related to $B$ as follows:

$$
\tilde {B} = (U _ {n - 1} ^ {T} \dots U _ {1} ^ {T}) B (G _ {1} V _ {2} \dots V _ {n - 1}) = \tilde {U} ^ {T} B \tilde {V}.
$$

Since each $V _ { i }$ has the form $V _ { i } = G ( i , i + 1 , \theta _ { i } )$ where $i = 2 { : } n - 1$ , it follows that $\bar { V } e _ { 1 } = Q e _ { 1 }$ . By the Implicit $Q$ theorem we can assert that V¯ and Q are essentially the same. Thus, we can implicitly effect the transition from $T$ to $\bar { T } = \bar { B } ^ { T } \bar { B }$ by working directly on the bidiagonal matrix $B .$ .

Of course, for these claims to hold it is necessary that the underlying tridiagonal matrices be unreduced. Since the subdiagonal entries of $B ^ { T } B$ are of the form $d _ { i } f _ { i }$ , it is clear that we must search the bidiagonal band for zeros. If $f _ { k } = 0$ for some $k$ , then

$$
B = \left[ \begin{array}{c c} B _ {1} & 0 \\ 0 & B _ {2} \end{array} \right] _ {n - k} ^ {k}
$$

and the original SVD problem decouples into two smaller problems involving the matrices $B _ { 1 }$ and $B _ { 2 }$ . If $d _ { k } = 0$ for some $k < n$ , then premultiplication by a sequence of Givens transformations can zero $f _ { k }$ . For example, if $n = 6$ and $k = 3$ , then by rotating in row planes (3,4), (3,5), and (3,6) we can zero the entire third row:

$$
B = \left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right] \stackrel {{(3, 4)}} {{\longrightarrow}} \left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & + & 0 \\ 0 & 0 & 0 & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right]
$$

$$
\stackrel {(3, 5)} {\longrightarrow} \left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & + \\ 0 & 0 & 0 & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right] \stackrel {(3, 6)} {\longrightarrow} \left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & \times & 0 \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

If $d _ { n } = 0$ , then the last column can be zeroed with a series of column rotations in planes $( n - 1 , n ) , ( n - 2 , n ) , \ldots , ( 1 , n )$ . Thus, we can decouple if $f _ { 1 } \cdot \cdot \cdot f _ { n - 1 } = 0$ or $d _ { 1 } \cdots d _ { n } =$ 0. Putting it all together we obtain the following SVD analogue of Algorithm 8.3.2.

Algorithm 8.6.1 (Golub-Kahan SVD Step) Given a bidiagonal matrix $B \in \mathbb { R } ^ { m \times n }$ having no zeros on its diagonal or superdiagonal, the following algorithm overwrites B with the bidiagonal matrix ${ \bar { B } } = { \bar { U } } ^ { T } B { \bar { V } }$ where U¯ and V¯ are orthogonal and V¯ is essentially the orthogonal matrix that would be obtained by applying Algorithm 8.3.2 to $T = B ^ { T } B$ .

Let $\mu$ be the eigenvalue of the trailing 2-by-2 submatrix of $T = B ^ { T } B$ that is closer to $t _ { n n }$ .

$$
y = t _ {1 1} - \mu
$$

$$
z = t _ {1 2}
$$

for $k = 1 { : } n - 1$

Determine $c = \cos ( \theta )$ and $s = \sin ( \theta )$ such that

$$
\left[ \begin{array}{c c} y & z \end{array} \right] \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] = \left[ \begin{array}{c c} * & 0 \end{array} \right].
$$

$$
B = B \cdot G (k, k + 1, \theta)
$$

$$
y = b _ {k k}
$$

$$
z = b _ {k + 1, k}
$$

Determine $c = \cos ( \theta )$ and $s = \sin ( \theta )$ such that

$$
\left[ \begin{array}{r r} c & s \\ - s & c \end{array} \right] ^ {T} \left[ \begin{array}{l} y \\ z \end{array} \right] = \left[ \begin{array}{l} * \\ 0 \end{array} \right].
$$

$$
B = G (k, k + 1, \theta) ^ {T} B
$$

$$
\text { if } k <   n - 1
$$

$$
y = b _ {k, k + 1}
$$

$$
z = b _ {k, k + 2}
$$

end

end

An efficient implementation of this algorithm would store $B \mathrm { { ^ { * } s } }$ diagonal and superdiagonal in vectors $d ( 1 { : } n )$ and $f ( 1 : n - 1 )$ , respectively, and would require 30n flops and 2n square roots. Accumulating U requires 6mn flops. Accumulating V requires $6 n ^ { 2 }$ flops.

Typically, after a few of the above SVD iterations, the superdiagonal entry $f _ { n - 1 }$ becomes negligible. Criteria for smallness within $B ^ { \prime } \mathrm { s }$ band are usually of the form

$$
\left| f _ {i} \right| \leq \operatorname{tol} \cdot \left(\left| d _ {i} \right| + \left| d _ {i + 1} \right|\right),
$$

$$
\left| d _ {i} \right| \leq \operatorname{tol} \cdot \| B \|,
$$

where tol is a small multiple of the unit roundoff and $\parallel \cdot \parallel$ is some computationally convenient norm. Combining Algorithm 5.4.2 (bidiagonalization), Algorithm 8.6.1, and the decoupling calculations mentioned earlier gives the following procedure.

Algorithm 8.6.2 (The SVD Algorithm) Given $A \in \mathbb { R } ^ { m \times n } \ ( m \geq n )$ and $\epsilon , \mathrm { ~ a ~ }$ small multiple of the unit roundoff, the following algorithm overwrites A with $U ^ { T } A V = D { + } E$ , where $U \in \mathbb { R } ^ { m \times m }$ is orthogonal, $V \in \mathbb { R } ^ { n \times n }$ is orthogonal, $D \in \mathbb { R } ^ { m \times n }$ is diagonal, and $E$ satisfies $\parallel E \parallel _ { 2 } \approx \mathbf { u } \parallel A \parallel _ { 2 }$ .

Use Algorithm 5.4.2 to compute the bidiagonalization.

$$
\left[ \begin{array}{c} B \\ 0 \end{array} \right] \leftarrow (U _ {1} \dots U _ {n}) ^ {T} A (V _ {1} \dots V _ {n - 2}).
$$

until $q = n$

For $i = 1 { : } n - 1$ , set $b _ { i , i + 1 }$ to zero if $| b _ { i , i + 1 } | \leq \epsilon ( | b _ { i i } | + | b _ { i + 1 , i + 1 } | )$

Find the largest q and the smallest p such that if

$$
B = \left[ \begin{array}{c c c} B _ {1 1} & 0 & 0 \\ 0 & B _ {2 2} & 0 \\ 0 & 0 & B _ {3 3} \end{array} \right] \begin{array}{c} p \\ n - p - q \\ q \end{array}
$$

then $B_{33}$ is diagonal and $B_{22}$ has a nonzero superdiagonal.
if q < n
    if any diagonal entry in $B_{22}$ is zero, then zero the superdiagonal entry in the same row.
    else
    Apply Algorithm 8.6.1 to $B_{22}$ . $B = \text{diag}(I_p, U, I_{q+m-n})^T B \text{diag}(I_p, V, I_q)$ end
    end
end

The amount of work required by this algorithm depends on how much of the SVD is required. For example, when solving the LS problem, $U ^ { T }$ need never be explicitly formed but merely applied to $b$ as it is developed. In other applications, only the matrix $U _ { 1 } = U ( : , 1 : n )$ is required. Another variable that affects the volume of work in Algorithm 8.6.2 concerns the R-bidiagonalization idea that we discussed in §5.4.9. Recall that unless A is “almost square,” it pays to reduce $A$ to triangular form via QR and before bidiagonalizing. If R-bidiagonalization is used in the SVD context, then we refer to the overall process as the R-SVD. Figure 8.6.1 summarizes the work associated with the various possibilities By comparing the entries in this table (which are meant only as approximate estimates of work), we conclude that the R-SVD approach is more efficient unless m ≈ n.

# 8.6.4 Jacobi SVD Procedures

It is straightforward to adapt the Jacobi procedures of §8.5 to the SVD problem. Instead of solving a sequence of 2-by-2 symmetric eigenproblems, we solve a sequence of 2-by-2 SVD problems. Thus, for a given index pair $( p , q )$ we compute a pair of rotations such that

<table><tr><td>Required</td><td>Golub-Reinsch SVD</td><td>R-SVD</td></tr><tr><td> $\Sigma$ </td><td> $4mn^{2} - 4n^{3}/3$ </td><td> $2mn^{2} + 2n^{3}$ </td></tr><tr><td> $\Sigma, V$ </td><td> $4mn^{2} + 8n^{3}$ </td><td> $2mn^{2} + 11n^{3}$ </td></tr><tr><td> $\Sigma, U$ </td><td> $4m^{2}n - 8mn^{2}$ </td><td> $4m^{2}n + 13n^{3}$ </td></tr><tr><td> $\Sigma, U_{1}$ </td><td> $14mn^{2} - 2n^{3}$ </td><td> $6mn^{2} + 11n^{3}$ </td></tr><tr><td> $\Sigma, U, V$ </td><td> $4m^{2}n + 8mn^{2} + 9n^{3}$ </td><td> $4m^{2}n + 22n^{3}$ </td></tr><tr><td> $\Sigma, U_{1}, V$ </td><td> $14mn^{2} + 8n^{3}$ </td><td> $6mn^{2} + 20n^{3}$ </td></tr></table>

Figure 8.6.1. Work associated with various SVD-related calculations

$$
\left[ \begin{array}{c c} {c _ {1}} & {s _ {1}} \\ {- s _ {1}} & {c _ {1}} \end{array} \right] ^ {T} \left[ \begin{array}{c c} {a _ {p p}} & {a _ {p q}} \\ {a _ {q p}} & {a _ {q q}} \end{array} \right] \left[ \begin{array}{c c} {c _ {2}} & {s _ {2}} \\ {- s _ {2}} & {c _ {2}} \end{array} \right] = \left[ \begin{array}{c c} {d _ {p}} & 0 \\ {0} & {d _ {q}} \end{array} \right].
$$

See P8.6.5. The resulting algorithm is referred to as two-sided because each update involves a pre- and a post-multiplication.

A one-sided Jacobi algorithm involves a sequence of pairwise column orthogonalizations. For a given index pair $( p , q )$ a Jacobi rotation $J ( p , q , \theta )$ is determined so that columns p and q of $A J ( p , q , \theta )$ are orthogonal to each other. See P8.6.8. Note that this corresponds to zeroing the $( p , q )$ and $( q , p )$ entries in $A ^ { T } A$ . Once $A V$ has sufficiently orthogonal columns, the rest of the SVD (U and Σ) follows from column scaling: $A V = U \Sigma$ .

# Problems

P8.6.1 Give formulae for the eigenvectors of

$$
S = \left[ \begin{array}{c c} 0 & A ^ {T} \\ A & 0 \end{array} \right]
$$

in terms of the singular vectors of $A \in \mathbb { R } ^ { m \times n }$ where m $\geq n$

P8.6.2 Relate the singular values and vectors of $A = B + i C \ ( B , C \in \mathbb { R } ^ { m \times n } )$ to those of

$$
\tilde {A} = \left[ \begin{array}{c c} B & - C \\ C & B \end{array} \right].
$$

P8.6.3 Suppose $B \in \mathbb { R } ^ { n \times n }$ is upper bidiagonal with diagonal entries $d ( 1 { : } n )$ and superdiagonal entries $f ( 1 : n - 1 )$ . State and prove a singular value version of Theorem 8.3.1.

P8.6.4 Assume that $n = 2 m$ and that $S \in \mathbb { R } ^ { n \times n }$ is skew-symmetric and tridiagonal. Show that there exists a permutation $P \in \mathbb { R } ^ { n \times n }$ such that

$$
P ^ {T} S P = \left[ \begin{array}{c c} 0 & - B ^ {T} \\ B & 0 \end{array} \right]
$$

where $B \in \mathbb { R } ^ { m \times m }$ . Describe the structure of B and show how to compute the eigenvalues and eigenvectors of S via the SVD of B. Repeat for the case $n = 2 m + 1$ .

P8.6.5 (a) Let

$$
C = \left[ \begin{array}{c c} w & x \\ y & z \end{array} \right]
$$

be real. Give a stable algorithm for computing c and s with $c ^ { 2 } + s ^ { 2 } = 1$ such that

$$
B = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] C
$$

is symmetric. (b) Combine (a) with Algorithm 8.5.1 to obtain a stable algorithm for computing the SVD of C. (c) Part (b) can be used to develop a Jacobi-like algorithm for computing the SVD of $A \in \mathbb { R } ^ { n \times n }$ . For a given (p, q) with $p \ < \ q .$ , Jacobi transformations $J ( p , q , \theta _ { 1 } )$ and $J ( p , q , \theta _ { 2 } )$ are determined such that if

$$
B = J (p, q, \theta_ {1}) ^ {T} A J (p, q, \theta_ {2}),
$$

then $b _ { p q } = b _ { q p } = 0$ . Show

$$
\operatorname{off} (B) ^ {2} = \operatorname{off} (A) ^ {2} - a _ {p q} ^ {2} - a _ {q p} ^ {2}.
$$

(d) Consider one sweep of a cyclic-by-row Jacobi SVD procedure applied to $A \in \mathbb { R } ^ { n \times n }$ :

$$
\begin{array}{l} \text { for } p = 1: n - 1 \\ \text { for } q = p + 1: n \\ A = J (p, q, \theta_ {1}) ^ {T} A J (p, q, \theta_ {2}) \\ \end{array}
$$

Assume that the Jacobi rotation matrices are chosen so that $\begin{array} { r } { a _ { p q } = a _ { q p } = 0 } \end{array}$ after the $( p , q )$ update. Show that if A is upper (lower) triangular at the beginning of the sweep, then it is lower (upper) triangular after the sweep is completed. See Kogbetliantz (1955). (e) How could these Jacobi ideas be used to compute the SVD of a rectangular matrix?

P8.6.6 Let x and y be in $\mathbb { R } ^ { m }$ and define the orthogonal matrix Q by

$$
Q   =   \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right].
$$

Give a stable algorithm for computing c and s such that the columns of $[ x \mid y ] Q$ are orthogonal to each other.

# Notes and References for 8.6

For a general perspective and overview of the SVD we recommend:

G.W. Stewart (1993). “On the Early History of the Singular Value Decomposition,” SIAM Review 35, 551–566.

A.K. Cline and I.S. Dhillon (2006). “Computation of the Singular Value Decomposition,” in Handbook of Linear Algebra, L. Hogben (ed.), Chapman and Hall, London, 45-1.

A perturbation theory for the SVD is developed in Stewart and Sun (MPT). See also:

P.A. Wedin (1972). “Perturbation Bounds in Connection with the Singular Value Decomposition,” BIT 12, 99–111.

G.W. Stewart (1973). “Error and Perturbation Bounds for Subspaces Associated with Certain Eigenvalue Problems,” SIAM Review 15, 727–764.

A. Ruhe (1975). “On the Closeness of Eigenvalues and Singular Values for Almost Normal Matrices,” Lin. Alg. Applic. 11, 87–94.

G.W. Stewart (1979). “A Note on the Perturbation of Singular Values,” Lin. Alg. Applic. 28, 213–216.

G.W. Stewart (1984). “A Second Order Perturbation Expansion for Small Singular Values,” Lin. Alg. Applic. 56, 231–236.

S. Chandrasekaren and I.C.F. Ipsen (1994). “Backward Errors for Eigenvalue and Singular Value Decompositions,” Numer. Math. 68, 215–223.

R.J. Vaccaro (1994). “A Second-Order Perturbation Expansion for the SVD,” SIAM J. Matrix Anal. Applic. 15, 661–671.

J. Sun (1996). “Perturbation Analysis of Singular Subspaces and Deflating Subspaces,” Numer. Math. 73, 235–263.   
F.M. Dopico (2000). “A Note on Sin T Theorems for Singular Subspace Variations BIT 40, 395–403.   
R.-C. Li and G. W. Stewart (2000). “A New Relative Perturbation Theorem for Singular Subspaces,” Lin. Alg. Applic. 313, 41–51.   
C.-K. Li and R. Mathias (2002). “Inequalities on Singular Values of Block Triangular Matrices,” SIAM J. Matrix Anal. Applic. 24, 126–131.   
F.M. Dopico and J. Moro (2002). “Perturbation Theory for Simultaneous Bases of Singular Subspaces,” BIT 42, 84–109.   
K.A. O’Neil (2005). “Critical Points of the Singular Value Decomposition,” SIAM J. Matrix Anal. Applic. 27, 459–473.   
M. Stewart (2006). “Perturbation of the SVD in the Presence of Small Singular Values,” Lin. Alg. Applic. 419, 53–77.   
H. Xiang (2006). “A Note on the Minimax Representation for the Subspace Distance and Singular Values,” Lin. Alg. Applic. 414, 470–473.   
W. Li and W. Sun (2007). “Combined Perturbation Bounds: I. Eigensystems and Singular Value Decompositions,” SIAM J. Matrix Anal. Applic. 29, 643–655.   
J. Matejaˇs and V. Hari (2008). “Relative Eigenvalues and Singular Value Perturbations of Scaled Diagonally Dominant Matrices,” BIT 48, 769–781.   
Classical papers that lay out the ideas behind the SVD algorithm include:   
G.H. Golub and W. Kahan (1965). “Calculating the Singular Values and Pseudo-Inverse of a Matrix,” SIAM J. Numer. Anal. 2, 205–224.   
P.A. Businger and G.H. Golub (1969). “Algorithm 358: Singular Value Decomposition of the Complex Matrix,” Commun. ACM 12, 564–565.   
G.H. Golub and C. Reinsch (1970). “Singular Value Decomposition and Least Squares Solutions,” Numer. Math. 14, 403–420.   
For related algorithmic developments and analysis, see:   
T.F. Chan (1982). “An Improved Algorithm for Computing the Singular Value Decomposition,” ACM Trans. Math. Softw. 8, 72–83.   
J.J.M. Cuppen (1983). “The Singular Value Decomposition in Product Form,” SIAM J. Sci. Stat. Comput. 4, 216–222.   
J.J. Dongarra (1983). “Improving the Accuracy of Computed Singular Values,” SIAM J. Sci. Stat. Comput. 4, 712–719.   
S. Van Huffel, J. Vandewalle, and A. Haegemans (1987). “An Efficient and Reliable Algorithm for Computing the Singular Subspace of a Matrix Associated with its Smallest Singular Values,” J. Comp. Appl. Math. 19, 313–330.   
P. Deift, J. Demmel, L.-C. Li, and C. Tomei (1991). “The Bidiagonal Singular Value Decomposition and Hamiltonian Mechanics,” SIAM J. Numer. Anal. 28, 1463–1516.   
R. Mathias and G.W. Stewart (1993). “A Block QR Algorithm and the Singular Value Decomposition,” Lin. Alg. Applic. 182, 91–100.   
V. Mehrmann and W. Rath (1993). “Numerical Methods for the Computation of Analytic Singular Value Decompositions,” ETNA 1, 72–88.   
˚A. Bj¨orck, E. Grimme, and P. Van Dooren (1994). “An Implicit Shift Bidiagonalization Algorithm for Ill-Posed Problems,” BIT 34, 510–534.   
K.V. Fernando and B.N. Parlett (1994). “Accurate Singular Values and Differential qd Algorithms,” Numer. Math. 67, 191–230.   
S. Chandrasekaran and I.C.F. Ipsen (1995). “Analysis of a QR Algorithm for Computing Singular Values,” SIAM J. Matrix Anal. Applic. 16, 520–535.   
U. von Matt (1997). “The Orthogonal qd–Algorithm,” SIAM J. Sci. Comput. 18, 1163–1186.   
K.V. Fernando (1998). “Accurately Counting Singular Values of Bidiagonal Matrices and Eigenvalues of Skew-Symmetric Tridiagonal Matrices,” SIAM J. Matrix Anal. Applic. 20, 373–399.   
N.J. Higham (2000). “QR factorization with Complete Pivoting and Accurate Computation of the SVD,” Lin. Alg. Applic. 309, 153–174.   
Divide-and-conquer methods for the bidiagonal SVD problem have been developed that are analogous to the tridiagonal eigenvalue strategies outlined in §8.4.4:   
J.W. Demmel and W. Kahan (1990). “Accurate Singular Values of Bidiagonal Matrices,” SIAM J. Sci. Stat. Comput. 11, 873–912.

E.R. Jessup and D.C. Sorensen (1994). “A Parallel Algorithm for Computing the Singular Value Decomposition of a Matrix,” SIAM J. Matrix Anal. Applic. 15, 530–548.   
M. Gu and S.C. Eisenstat (1995). “A Divide-and-Conquer Algorithm for the Bidiagonal SVD,” SIAM J. Matrix Anal. Applic. 16, 79–92.   
P.R. Willems, B. Lang, and C. V¨omel (2006). “Computing the Bidiagonal SVD Using Multiple Relatively Robust Representations,” SIAM J. Matrix Anal. Applic. 28, 907–926.   
T. Konda and Y. Nakamura (2009). “A New Algorithm for Singular Value Decomposition and Its Parallelization,” Parallel Comput. 35, 331–344.

For structured SVD problems, there are interesting, specialized results, see:

S. Van Huffel and H. Park (1994). “Parallel Tri- and Bidiagonalization of Bordered Bidiagonal Matrices,” Parallel Comput. 20, 1107–1128.   
J. Demmel and P. Koev (2004). “Accurate SVDs of Weakly Diagonally Dominant M-matrices,” Num. Math. 98, 99–104.   
N. Mastronardi, M. Van Barel, and R. Vandebril (2008). “A Fast Algorithm for the Recursive Calculation of Dominant Singular Subspaces,” J. Comp. Appl. Math. 218, 238–246.

Jacobi methods for the SVD fall into two categories. The two-sided Jacobi algorithms repeatedly perform the update $A  U ^ { T } A V$ producing a sequence of iterates that are increasingly diagonal.

E.G. Kogbetliantz (1955). “Solution of Linear Equations by Diagonalization of Coefficient Matrix,” Quart. Appl. Math. 13, 123–132. G.E. Forsythe and P. Henrici (1960). “The Cyclic Jacobi Method for Computing the Principal Values of a Complex Matrix,” Trans. AMS 94, 1–23.

C.C. Paige and P. Van Dooren (1986). “On the Quadratic Convergence of Kogbetliantz’s Algorithm for Computing the Singular Value Decomposition,” Lin. Alg. Applic. 77, 301–313.

J.P. Charlier and P. Van Dooren (1987). “On Kogbetliantz’s SVD Algorithm in the Presence of Clusters,” Lin. Alg. Applic. 95, 135–160.   
Z. Bai (1988). “Note on the Quadratic Convergence of Kogbetliantz’s Algorithm for Computing the Singular Value Decomposition,” Lin. Alg. Applic. 104, 131–140.   
J.P. Charlier, M. Vanbegin, and P. Van Dooren (1988). “On Efficient Implementation of Kogbetliantz’s Algorithm for Computing the Singular Value Decomposition,” Numer. Math. 52, 279–300.   
K.V. Fernando (1989). “Linear Convergence of the Row-Cyclic Jacobi and Kogbetliantz Methods,” Numer. Math. 56, 73–92.   
Z. Drmaˇc and K. Veseliˇc (2008). “New Fast and Accurate Jacobi SVD Algorithm I,” SIAM J. Matrix Anal. Applic. 29, 1322–1342.

The one-sided Jacobi SVD procedures repeatedly perform the update A  AV producing a sequence of iterates with columns that are increasingly orthogonal, see:

J.C. Nash (1975). “A One-Sided Tranformation Method for the Singular Value Decomposition and Algebraic Eigenproblem,” Comput. J. 18, 74–76.   
P.C. Hansen (1988). “Reducing the Number of Sweeps in Hestenes Method,” in Singular Value Decomposition and Signal Processing, E.F. Deprettere (ed.) North Holland, Amsterdam.   
K. Veseli˘c and V. Hari (1989). “A Note on a One-Sided Jacobi Algorithm,” Numer. Math. 56, 627–633.

Careful implementation and analysis has shown that Jacobi SVD has remarkably accuracy:

J. Demmel, M. Gu, S. Eisenstat, I. Slapnicar, K. Veseli´c, and Z. Drmaˇc (1999). “Computing the Singular Value Decomposition with High Relative Accuracy,” Lin. Alg. Applic. 299, 21–80.   
Z Drmaˇc (1999). “A Posteriori Computation of the Singular Vectors in a Preconditioned Jacobi SVD Algorithm,” IMA J. Numer. Anal. 19, 191–213.   
Z. Drmaˇc (1997). “Implementation of Jacobi Rotations for Accurate Singular Value Computation in Floating Point Arithmetic,” SIAM J. Sci. Comput. 18, 1200–1222.   
F.M. Dopico and J. Moro (2004). “A Note on Multiplicative Backward Errors of Accurate SVD Algorithms,” SIAM J. Matrix Anal. Applic. 25, 1021–1031.

The parallel implementation of the Jacobi SVD has a long and interesting history:

F.T. Luk (1980). “Computing the Singular Value Decomposition on the ILLIAC IV,” ACM Trans. Math. Softw. 6, 524–539.

R.P. Brent and F.T. Luk (1985). “The Solution of Singular Value and Symmetric Eigenvalue Problems on Multiprocessor Arrays,” SIAM J. Sci. Stat. Comput. 6, 69–84.   
R.P. Brent, F.T. Luk, and C. Van Loan (1985). “Computation of the Singular Value Decomposition Using Mesh Connected Processors,” J. VLSI Computer Systems 1, 242–270.   
F.T. Luk (1986). “A Triangular Processor Array for Computing Singular Values,” Lin. Alg. Applic. 77, 259–274.   
M. Berry and A. Sameh (1986). “Multiprocessor Jacobi Algorithms for Dense Symmetric Eigenvalue and Singular Value Decompositions,” in Proceedings International Conference on Parallel Processing, 433–440.   
R. Schreiber (1986). “Solving Eigenvalue and Singular Value Problems on an Undersized Systolic Array,” SIAM J. Sci. Stat. Comput. 7, 441–451.   
C.H. Bischof and C. Van Loan (1986). “Computing the SVD on a Ring of Array Processors,” in Large Scale Eigenvalue Problems, J. Cullum and R. Willoughby (eds.), North Holland, Amsterdam, 51– 66.   
C.H. Bischof (1987). “The Two-Sided Block Jacobi Method on Hypercube Architectures,” in Hypercube Multiprocessors, M.T. Heath (ed.), SIAM Publications, Philadelphia, PA.   
C.H. Bischof (1989). “Computing the Singular Value Decomposition on a Distributed System of Vector Processors,” Parallel Comput. 11, 171–186.   
M. Beˇca, G. Okˇsa, M. Vajterˇsic, and L. Grigori (2010). “On Iterative QR Pre-Processing in the Parallel Block-Jacobi SVD Algorithm,” Parallel Comput. 36, 297–307.
