# 7.4 The Hessenberg and Real Schur Forms

In this and the next section we show how to make the QR iteration (7.3.1) a fast, effective method for computing Schur decompositions. Because the majority of eigenvalue/invariant subspace problems involve real data, we concentrate on developing the real analogue of (7.3.1) which we write as follows:

$$
H _ {0} = U _ {0} ^ {T} A U _ {0}
$$

for k = 1, 2, . . .

$$
H _ {k - 1} = U _ {k} R _ {k} \quad \text {(QR factorization)} \tag {7.4.1}
$$

$$
H _ {k} = R _ {k} U _ {k}
$$

end

Here, $A \in \mathbb { R } ^ { n \times n }$ , each $U _ { k } \in \mathbb { R } ^ { n \times n }$ is orthogonal, and each $R _ { k } \in \mathbb { R } ^ { n \times n }$ is upper triangular. A difficulty associated with this real iteration is that the $H _ { k }$ can never converge to triangular form in the event that A has complex eigenvalues. For this reason, we must lower our expectations and be content with the calculation of an alternative decomposition known as the real Schur decomposition.

In order to compute the real Schur decomposition efficiently we must carefully choose the initial orthogonal similarity transformation $U _ { 0 }$ in (7.4.1). In particular, if we choose $U _ { 0 }$ so that $H _ { 0 }$ is upper Hessenberg, then the amount of work per iteration is reduced from $O ( n ^ { 3 } )$ to $O ( n ^ { 2 } )$ . The initial reduction to Hessenberg form (the $U _ { 0 }$ computation) is a very important computation in its own right and can be realized by a sequence of Householder matrix operations.

# 7.4.1 The Real Schur Decomposition

A block upper triangular matrix with either 1-by-1 or 2-by-2 diagonal blocks is upper quasi-triangular. The real Schur decomposition amounts to a real reduction to upper quasi-triangular form.

Theorem 7.4.1 (Real Schur Decomposition). If $A \in \mathbb { R } ^ { n \times n }$ , then there exists an orthogonal $Q \in \mathbb { R } ^ { n \times n }$ such that

$$
Q ^ {T} A Q = \left[ \begin{array}{c c c c} R _ {1 1} & R _ {1 2} & \dots & R _ {1 m} \\ 0 & R _ {2 2} & \dots & R _ {2 m} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & R _ {m m} \end{array} \right] \tag {7.4.2}
$$

where each $R _ { i i }$ is either a 1-by-1 matrix or a ${ \it 2 - b y - 2 }$ matrix having complex conjugate eigenvalues.

Proof. The complex eigenvalues of A occur in conjugate pairs since the characteristic polynomial det $( z I - A )$ has real coefficients. Let k be the number of complex conjugate pairs in $\lambda ( A )$ . We prove the theorem by induction on k. Observe first that Lemma 7.1.2 and Theorem 7.1.3 have obvious real analogs. Thus, the theorem holds if $k = 0$ . Now suppose that $k \geq 1$ . If $\lambda = \gamma + i \mu \in \lambda ( A )$ and $\mu \neq 0$ , then there exist vectors y and z in $\mathbb { R } ^ { n } ( z \neq 0 )$ such that $A ( y + i z ) ~ = ~ ( \gamma + i \mu ) ( y + i z )$ , i.e.,

$$
A \left[ \begin{array}{c c} y & z \end{array} \right] = \left[ \begin{array}{c c} y & z \end{array} \right] \left[ \begin{array}{c c} \gamma & \mu \\ - \mu & \gamma \end{array} \right].
$$

The assumption that $\mu \neq 0$ implies that y and z span a 2-dimensional, real invariant subspace for A. It then follows from Lemma 7.1.2 that an orthogonal $U \in \mathbb { R } ^ { n \times n }$ exists such that

$$
U ^ {T} A U = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {n - 2} ^ {2}
$$

where $\lambda ( T _ { 1 1 } ) = \{ \lambda , \bar { \lambda } \}$ . By induction, there exists an orthogonal $\tilde { U }$ so $\tilde { U } ^ { T } T _ { 2 2 } \tilde { U }$ has the required structure. The theorem follows by setting $Q = U \cdot \mathrm { d i a g } ( I _ { 2 } , \tilde { U } )$ .

The theorem shows that any real matrix is orthogonally similar to an upper quasitriangular matrix. It is clear that the real and imaginary parts of the complex eigenvalues can be easily obtained from the 2-by-2 diagonal blocks. Thus, it can be said that the real Schur decomposition is an eigenvalue-revealing decomposition.

# 7.4.2 A Hessenberg QR Step

We now turn our attention to the efficient execution of a single QR step in (7.4.1). In this regard, the most glaring shortcoming associated with (7.4.1) is that each step requires a full QR factorization costing $O ( n ^ { \bar { 3 } } )$ flops. Fortunately, the amount of work per iteration can be reduced by an order of magnitude if the orthogonal matrix $U _ { 0 }$ is judiciously chosen. In particular, if $U _ { 0 } ^ { T } A U _ { 0 } = H _ { 0 } = ( h _ { i j } )$ is upper Hessenberg $( h _ { i j } = 0 ,$ , $i > j + 1 \}$ ), then each subsequent $H _ { k }$ requires only $O ( n ^ { 2 } )$ flops to calculate. To see this we look at the computations $H = Q R$ and $H _ { + } = R Q$ when H is upper Hessenberg. As described in §5.2.5, we can upper triangularize H with a sequence of $n - 1$ Givens rotations: $Q ^ { T } H \equiv G _ { n - 1 } ^ { T } \cdot \cdot \cdot G _ { 1 } ^ { T } \bar { H = } R$ . Here, $G _ { i } = G ( i , i + 1 , \theta _ { i } )$ . For the $n = 4$ case there are three Givens premultiplications:

$$
\left[ \begin{array}{c c c c} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right].
$$

See Algorithm 5.2.5. The computation $R Q = R ( G _ { 1 } \cdot \cdot \cdot G _ { n - 1 } )$ is equally easy to implement. In the $n = 4$ case there are three Givens post-multiplications:

$$
\left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \end{array} \right].
$$

Overall we obtain the following algorithm:

Algorithm 7.4.1 If H is an n-by-n upper Hessenberg matrix, then this algorithm overwrites H with $H _ { + } = R Q$ where $H = Q R$ is the QR factorization of H.

for $k = 1 { : } n - 1$

$$
[ c _ {k}, s _ {k} ] = \operatorname{givens} (H (k, k), H (k + 1, k))
$$

$$
H (k: k + 1, k: n) = \left[ \begin{array}{c c} c _ {k} & s _ {k} \\ - s _ {k} & c _ {k} \end{array} \right] ^ {T} H (k: k + 1, k: n)
$$

end

for $k = 1 { : } n - 1$

$$
H (1: k + 1, k: k + 1) = H (1: k + 1, k: k + 1) \left[ \begin{array}{c c} c _ {k} & s _ {k} \\ - s _ {k} & c _ {k} \end{array} \right]
$$

end

Let $G _ { k } = G ( k , k { + } 1 , \theta _ { k } )$ be the kth Givens rotation. It is easy to confirm that the matrix $Q = G _ { 1 } \cdot \cdot \cdot G _ { n - 1 }$ is upper Hessenberg. Thus, $R Q = H _ { + }$ is also upper Hessenberg. The algorithm requires about $6 n ^ { 2 }$ flops, an order of magnitude more efficient than a full matrix QR step (7.3.1).

# 7.4.3 The Hessenberg Reduction

It remains for us to show how the Hessenberg decomposition

$$
U _ {0} ^ {T} A U _ {0} = H, \quad U _ {0} ^ {T} U _ {0} = I \tag {7.4.3}
$$

can be computed. The transformation $U _ { 0 }$ can be computed as a product of Householder matrices $P _ { 1 } , \ldots , P _ { n - 2 }$ . The role of $P _ { k }$ is to zero the kth column below the subdiagonal. In the $n = 6$ case, we have

$$
\left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \end{array} \right] \xrightarrow {P _ {1}} \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \end{array} \right] \xrightarrow {P _ {2}}
$$

$$
\left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\end{array}\right] \stackrel {{P _ {3}}} {{\rightarrow}} \left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&0&\times&\times&\times\\0&0&0&\times&\times&\times\end{array}\right] \stackrel {{P _ {4}}} {{\rightarrow}} \left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&0&\times&\times&\times\\0&0&0&0&\times&\times\end{array}\right].
$$

In general, after $k - 1$ steps we have computed $k - 1$ Householder matrices $P _ { 1 } , \ldots , P _ { k - 1 }$ such that

$$
(P _ {1} \dots P _ {k - 1}) ^ {T} A (P _ {1} \dots P _ {k - 1}) = \left[ \begin{array}{c c c} B _ {1 1} & B _ {1 2} & B _ {1 3} \\ B _ {2 1} & B _ {2 2} & B _ {2 3} \\ 0 & B _ {3 2} & B _ {3 3} \end{array} \right] \begin{array}{c} k - 1 \\ 1 \\ n - k \end{array}
$$

is upper Hessenberg through its first $k - 1$ columns. Suppose $\tilde { P } _ { k }$ is an order- $( n - k )$ Householder matrix such that $\tilde { P } _ { k } B _ { 3 2 }$ is a multiple of $e _ { 1 } ^ { ( n - k ) }$ . If $P _ { k } = \mathrm { d i a g } ( I _ { k } , \tilde { P } _ { k } )$ , then

$$
(P _ {1} \dots P _ {k}) ^ {T} A (P _ {1} \dots P _ {k}) = \left[ \begin{array}{c c c} B _ {1 1} & B _ {1 2} & B _ {1 3} \tilde {P} _ {k} \\ B _ {2 1} & B _ {2 2} & B _ {2 3} \tilde {P} _ {k} \\ 0 & \tilde {P} _ {k} B _ {3 2} & \tilde {P} _ {k} B _ {3 3} \tilde {P} _ {k} \end{array} \right]
$$

is upper Hessenberg through its first k columns. Repeating this for $k = 1 { : } n - 2$ we obtain

Algorithm 7.4.2 (Householder Reduction to Hessenberg Form) Given $A \in \mathbb { R } ^ { n \times n }$ , the following algorithm overwrites A with $H = U _ { 0 } ^ { T } A U _ { 0 }$ where H is upper Hessenberg and $U _ { 0 }$ is a product of Householder matrices.

$$
\begin{array}{l} \text { for } k = 1: n - 2 \\ [ v, \beta ] = \text { house } (A (k + 1: n, k)) \\ A (k + 1: n, k: n) = (I - \beta v v ^ {T}) A (k + 1: n, k: n) \\ A (1: n, k + 1: n) = A (1: n, k + 1: n) (I - \beta v v ^ {T}) \\ \end{array}
$$

end

This algorithm requires $1 0 n ^ { 3 } / 3$ flops. If $U _ { 0 }$ is explicitly formed, an additional $4 n ^ { 3 } / 3$ flops are required. The kth Householder matrix can be represented in $A ( k + 2 { : } n , k )$ . See Martin and Wilkinson (1968) for a detailed description.

The roundoff properties of this method for reducing A to Hessenberg form are very desirable. Wilkinson (AEP, p. 351) states that the computed Hessenberg matrix $\hat { H }$ satisfies

$$
\hat {H} = Q ^ {T} (A + E) Q,
$$

where $Q$ is orthogonal and $\| E \| _ { F } \leq c n ^ { 2 } \mathbf { u } \| A \| _ { F }$ with c a small constant.

# 7.4.4 Level-3 Aspects

The Hessenberg reduction (Algorithm 7.4.2) is rich in level-2 operations: half gaxpys and half outer product updates. We briefly mention two ideas for introducing level-3 computations into the process.

The first involves a block reduction to block Hessenberg form and is quite straightforward. Suppose (for clarity) that $n = r N$ and write

$$
A   =   \left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right] _ {n - r} ^ {r}    .
$$

Suppose that we have computed the QR factorization $A _ { 2 1 } = \tilde { Q } _ { 1 } R _ { 1 }$ and that $\tilde { Q } _ { 1 }$ is in WY form. That is, we have $W _ { 1 } , Y _ { 1 } \in \dot { \mathbb { R } } ^ { ( n - r ) \times r }$ such that $\tilde { Q } _ { 1 } = I + W _ { 1 } Y _ { 1 } ^ { T }$ . (See §5.2.2 for details.) If $Q _ { 1 } = \mathrm { d i a g } ( I _ { r } , \tilde { Q } _ { 1 } )$ then

$$
Q _ {1} ^ {T} A Q _ {1} = \left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \tilde {Q} _ {1} \\ R _ {1} & \tilde {Q} _ {1} ^ {T} A _ {2 2} \tilde {Q} _ {1} \end{array} \right].
$$

Notice that the updates of the (1,2) and (2,2) blocks are rich in level-3 operations given that $\tilde { Q } _ { 1 }$ is in WY form. This fully illustrates the overall process as $Q _ { 1 } ^ { T } A Q _ { 1 }$ is block upper Hessenberg through its first block column. We next repeat the computations on the first r columns of $\tilde { Q } _ { 1 } ^ { T } A _ { 2 2 } \tilde { Q } _ { 1 }$ . After $N - 1$ such steps we obtain

$$
H = U _ {0} ^ {T} A U _ {0} = \left[ \begin{array}{c c c c c} H _ {1 1} & H _ {1 2} & \dots & \dots & H _ {1 N} \\ H _ {2 1} & H _ {2 2} & \dots & \dots & H _ {2 N} \\ 0 & \ddots & \ddots & \dots & \vdots \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & \dots & H _ {N, N - 1} & H _ {N N} \end{array} \right]
$$

where each $H _ { i j }$ is $r { \mathrm { - } } \mathrm { b y } { \mathrm { - } } r$ and $U _ { 0 } = Q _ { 1 } \cdot \cdot \cdot Q _ { N - 2 }$ with each $Q _ { i }$ in WY form. The overall algorithm has a level-3 fraction of the form $1 - O ( 1 / N )$ . Note that the subdiagonal blocks in H are upper triangular and so the matrix has lower bandwidth $r .$ . It is possible to reduce H to actual Hessenberg form by using Givens rotations to zero all but the first subdiagonal.

Dongarra, Hammarling, and Sorensen (1987) have shown how to proceed directly to Hessenberg form using a mixture of gaxpys and level-3 updates. Their idea involves minimal updating after each Householder transformation is generated. For example, suppose the first Householder $P _ { 1 }$ has been computed. To generate $P _ { 2 }$ we need just the second column of $P _ { 1 } A P _ { 1 }$ , not the full outer product update. To generate $P _ { 3 }$ we need just the thirrd column of $P _ { 2 } P _ { 1 } A P _ { 1 } P _ { 2 }$ , etc. In this way, the Householder matrices can be determined using only gaxpy operations. No outer product updates are involved. Once a suitable number of Householder matrices are known they can be aggregated and applied in level-3 fashion.

For more about the challenges of organizing a high-performance Hessenberg reduction, see Karlsson (2011).

# 7.4.5 Important Hessenberg Matrix Properties

The Hessenberg decomposition is not unique. If Z is any n-by-n orthogonal matrix and we apply Algorithm 7.4.2 to $Z ^ { T } A Z$ , then $Q ^ { T } A Q = H$ is upper Hessenberg where $Q \mathrm { ~ = ~ } Z U _ { 0 }$ . However, $Q e _ { 1 } = Z ( U _ { 0 } e _ { 1 } ) = Z e _ { 1 }$ suggesting that H is unique once the first column of Q is specified. This is essentially the case provided H has no zero subdiagonal entries. Hessenberg matrices with this property are said to be unreduced. Here is important theorem that clarifies these issues.

Theorem 7.4.2 ( Implicit Q Theorem ). Suppose $Q = { \left[ \begin{array} { l } { q _ { 1 } \left| \cdots \right| q _ { n } } \end{array} \right] }$ and $V =$ $[  v _ { 1 } | \cdots | v _ { n } ]$ are orthogonal matrices with the property that the matrices $Q ^ { T } A Q \ = H$ and $V ^ { T } A V = G$ are each upper Hessenberg where $A \in \mathbb { R } ^ { n \times n }$ . Let k denote the smallest positive integer for which $h _ { k + 1 , k } = 0 ,$ , with the convention that $k = n$ if H is unreduced. $H f q _ { 1 } = v _ { 1 }$ , then $q _ { i } = \pm v _ { i }$ and $| h _ { i , i - 1 } | = | g _ { i , i - 1 } | f o r i = 2 { : } k$ . Moreover, if $k < n$ , then $g _ { k + 1 , k } = 0$ .

Proof. Define the orthogonal matrix $W = \left[ \left. w _ { 1 } \right| \cdot \cdot \cdot \mid w _ { n } \right] = V ^ { T } Q$ and observe that $G W = W H$ . By comparing column i − 1 in this equation for $i = 2 { : } k$ we see that

$$
h _ {i, i - 1} w _ {i} = G w _ {i - 1} - \sum_ {j = 1} ^ {i - 1} h _ {j, i - 1} w _ {j}.
$$

Since $w _ { 1 } = e _ { 1 }$ , it follows that $\left[ \boldsymbol { w } _ { 1 } \left| \cdots \right| \boldsymbol { w } _ { k } \right]$ is upper triangular and so for $i = 2 { : } k$ we have $w _ { i } = \pm I _ { n } ( : , i ) = \pm e _ { i }$ . Since $w _ { i } = V ^ { T } q _ { i }$ and $h _ { i , i - 1 } = w _ { i } ^ { T } G w _ { i - 1 }$ it follows that $v _ { i } = \pm q _ { i }$ and

$$
| h _ {i, i - 1} | = | q _ {i} ^ {T} A q _ {i - 1} | = | v _ {i} ^ {T} A v _ {i - 1} | = | g _ {i, i - 1} |
$$

for i = 2:k. If $k < n$ , then

$$
g _ {k + 1, k} = e _ {k + 1} ^ {T} G e _ {k} = \pm e _ {k + 1} ^ {T} G W e _ {k} = \pm e _ {k + 1} ^ {T} W H e _ {k}
$$

$$
= \pm e _ {k + 1} ^ {T} \sum_ {i = 1} ^ {k} h _ {i k} W e _ {i} = \pm \sum_ {i = 1} ^ {k} h _ {i k} e _ {k + 1} ^ {T} e _ {i} = 0,
$$

completing the proof of the theorem.

The gist of the implicit Q theorem is that if $Q ^ { T } A Q = H$ and $Z ^ { T } A Z = G$ are each unreduced upper Hessenberg matrices and Q and Z have the same first column, then G and H are “essentially equal” in the sense that $G = D ^ { - 1 } H D$ where $D = \mathrm { d i a g } ( \pm 1 , \ldots , \pm 1 )$ .

Our next theorem involves a new type of matrix called a Krylov matrix. If $A \in \mathbb { R } ^ { n \times n }$ and $v \in \mathbb { R } ^ { n }$ , then the Krylov matrix $K ( A , v , j ) \in \mathbb { R } ^ { n \times j }$ is defined by

$$
K (A, v, j) = \left[ v \mid A v \mid \dots \mid A ^ {j - 1} v \right].
$$

It turns out that there is a connection between the Hessenberg reduction $Q ^ { T } A Q = H$ and the QR factorization of the Krylov matrix $K ( A , Q ( : , 1 ) , n )$ .

Theorem 7.4.3. Suppose $Q \in \mathbb { R } ^ { n \times n }$ is an orthogonal matrix and $A \in \mathbb { R } ^ { n \times n }$ . Then $Q ^ { T } A Q = H$ is an unreduced upper Hessenberg matrix if and only i ${ } ^ { \sharp } Q ^ { T } K ( A , Q ( : , 1 ) , n ) =$ R is nonsingular and upper triangular.

Proof. Suppose $Q \in \mathbb { R } ^ { n \times n }$ is orthogonal and set $H = Q ^ { T } A Q$ . Consider the identity

$$
Q ^ {T} K (A, Q (:, 1), n) = \left[ e _ {1} \mid H e _ {1} \mid \dots \mid H ^ {n - 1} e _ {1} \right] \equiv R.
$$

If H is an unreduced upper Hessenberg matrix, then it is clear that R is upper triangular with $r _ { i i } = h _ { 2 1 } h _ { 3 2 } \cdot \cdot \cdot h _ { i , i - 1 }$ for $i = 2 { : } n$ . Since $r _ { 1 1 } = 1$ it follows that R is nonsingular.

To prove the converse, suppose R is upper triangular and nonsingular. Since $R ( : , k + 1 ) = H R ( : , k )$ it follows that $H ( : , k ) \in { \mathsf { s p a n } } { \left\{ \begin{array} { l l } { e _ { 1 } , \ldots , e _ { k + 1 } } \end{array} \right\} }$ . This implies that H is upper Hessenberg. Since $r _ { n n } = h _ { 2 1 } h _ { 3 2 } \cdot \cdot \cdot h _ { n , n - 1 } \neq 0$ it follows that H is also unreduced.

Thus, there is more or less a correspondence between nonsingular Krylov matrices and orthogonal similarity reductions to unreduced Hessenberg form.

Our last result is about the geometric multiplicity of an eigenvalue of an unreduced upper Hessenberg matrix.

Theorem 7.4.4. If λ is an eigenvalue of an unreduced upper Hessenberg matrix $H \in \mathbb { R } ^ { n \times n }$ , then its geometric multiplicity is 1.

Proof. For any $\lambda \in \mathbb { C }$ we have rank $( A - \lambda I ) \geq n - 1$ because the first $n - 1$ columns of $H - \lambda I$ are independent.

# 7.4.6 Companion Matrix Form

Just as the Schur decomposition has a nonunitary analogue in the Jordan decomposition, so does the Hessenberg decomposition have a nonunitary analog in the companion matrix decomposition. Let $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and suppose that the Krylov matrix $K = K ( A , x , n )$ is nonsingular. If $c = c ( 0 { : } n - 1 )$ solves the linear system $K c = - A ^ { n } x$ , then it follows that $A K = K C$ where C has the form

$$
C = \left[ \begin{array}{c c c c c} 0 & 0 & \dots & 0 & - c _ {0} \\ 1 & 0 & \dots & 0 & - c _ {1} \\ 0 & 1 & \dots & 0 & - c _ {2} \\ \vdots & \vdots & \vdots & \vdots & \vdots \\ 0 & 0 & \dots & 1 & - c _ {n - 1} \end{array} \right]. \tag {7.4.4}
$$

The matrix C is said to be a companion matrix. Since

$$
\det (z I - C) = c _ {0} + c _ {1} z + \dots + c _ {n - 1} z ^ {n - 1} + z ^ {n},
$$

it follows that if K is nonsingular, then the decomposition $K ^ { - 1 } A K = C$ displays A’s characteristic polynomial. This, coupled with the sparseness of C, leads to “companion matrix methods” in various application areas. These techniques typically involve:

Step 1. Compute the Hessenberg decomposition $U _ { 0 } ^ { T } A U _ { 0 } = H$ .

Step 2. Hope H is unreduced and set $Y = \left[ e _ { 1 } \mid H e _ { 1 } \mid . . . \mid H ^ { n - 1 } e _ { 1 } \right]$ .

Step 3. Solve $Y C = H Y$ for C.

Unfortunately, this calculation can be highly unstable. A is similar to an unreduced Hessenberg matrix only if each eigenvalue has unit geometric multiplicity. Matrices that have this property are called nonderogatory. It follows that the matrix Y above can be very poorly conditioned if A is close to a derogatory matrix.

A full discussion of the dangers associated with companion matrix computation can be found in Wilkinson (AEP, pp. 405ff.).

# Problems

P7.4.1 Suppose $A \in \mathbb { R } ^ { n \times n }$ and $z \in \mathbb { R } ^ { n }$ . Give a detailed algorithm for computing an orthogonal Q such that $\bar { Q } ^ { \bar { T } } A Q$ is upper Hessenberg and $Q ^ { T } .$ z is a multiple of e1. Hint: Reduce z first and then apply Algorithm 7.4.2.

P7.4.2 Develop a similarity reduction to Hessenberg form using Gauss transforms with pivoting. How many flops are required. See Businger (1969).

P7.4.3 In some situations, it is necessary to solve the linear system $( A + z I ) x = b$ for many different values of $z \in \mathbb { R }$ and $b \in \mathbb { R } ^ { n }$ . Show how this problem can be efficiently and stably solved using the Hessenberg decomposition.

P7.4.4 Suppose $H \in \mathbb { R } ^ { n \times n }$ is an unreduced upper Hessenberg matrix. Show that there exists a diagonal matrix D such that each subdiagonal element of $D ^ { - 1 } H D$ is equal to 1. What is $\kappa _ { 2 } ( D ) ?$

P7.4.5 Suppose $W , Y \in \mathbb { R } ^ { n \times n }$ and define the matrices C and B by

$$
C = W + i Y, \qquad B = \left[ \begin{array}{c c} W & - Y \\ Y & W \end{array} \right].
$$

Show that if $\lambda \in \lambda ( C )$ is real, then $\lambda \in \lambda ( B )$ . Relate the corresponding eigenvectors.

P7.4.6 Suppose

$$
A = \left[ \begin{array}{c c} w & x \\ y & z \end{array} \right]
$$

is a real matrix having eigenvalues $\lambda \pm i \mu ,$ , where µ is nonzero. Give an algorithm that stably determines $c = \cos ( \theta )$ and $s = \sin ( \theta )$ such that

$$
\left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} \left[ \begin{array}{c c} w & x \\ y & z \end{array} \right] \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] = \left[ \begin{array}{c c} \lambda & \beta \\ \alpha & \lambda \end{array} \right]
$$

where $\alpha \beta = - \mu ^ { 2 }$ .

P7.4.7 Suppose $( \lambda , x )$ is a known eigenvalue-eigenvector pair for the upper Hessenberg matrix $H \in \mathbb { R } ^ { n \times n }$ . Give an algorithm for computing an orthogonal matrix P such that

$$
P ^ {T} H P = \left[ \begin{array}{l l} \lambda & w ^ {T} \\ 0 & H _ {1} \end{array} \right]
$$

where $H _ { 1 } \in \mathbb { R } ^ { ( n - 1 ) \times ( n - 1 ) }$ is upper Hessenberg. Compute P as a product of Givens rotations.

P7.4.8 Suppose $H \in \mathbb { R } ^ { n \times n }$ has lower bandwidth p. Show how to compute $Q \in \mathbb { R } ^ { n \times n } ;$ , a product of Givens rotations, such that $Q ^ { T } H Q$ is upper Hessenberg. How many flops are required?

P7.4.9 Show that if C is a companion matrix with distinct eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { n }$ , then $V C V ^ { - 1 } =$ $\mathrm { d i a g } ( \lambda _ { 1 } , \ldots , \lambda _ { n } )$ where

$$
V = \left[ \begin{array}{c c c c} 1 & \lambda_ {1} & \dots & \lambda_ {1} ^ {n - 1} \\ 1 & \lambda_ {2} & \dots & \lambda_ {2} ^ {n - 1} \\ \vdots & \vdots & \ddots & \vdots \\ 1 & \lambda_ {n} & \dots & \lambda_ {n} ^ {n - 1} \end{array} \right].
$$

# Notes and References for 7.4

The real Schur decomposition was originally presented in:

F.D. Murnaghan and A. Wintner (1931). “A Canonical Form for Real Matrices Under Orthogonal Transformations,” Proc. Nat. Acad. Sci. 17, 417–420.

A thorough treatment of the reduction to Hessenberg form is given in Wilkinson (AEP, Chap. 6), and Algol procedures appear in:   
R.S. Martin and J.H. Wilkinson (1968). “Similarity Reduction of a General Matrix to Hessenberg Form,” Numer. Math. 12, 349–368.   
Givens rotations can also be used to compute the Hessenberg decomposition, see:   
W. Rath (1982). “Fast Givens Rotations for Orthogonal Similarity,” Numer. Math. 40, 47–56.   
The high-performance computation of the Hessenberg reduction is a major challenge because it is a two-sided factorization, see:   
J.J. Dongarra, L. Kaufman, and S. Hammarling (1986). “Squeezing the Most Out of Eigenvalue Solvers on High Performance Computers,” Lin. Alg. Applic. 77, 113–136.   
J.J. Dongarra, S. Hammarling, and D.C. Sorensen (1989). “Block Reduction of Matrices to Condensed Forms for Eigenvalue Computations,” J. ACM 27, 215–227.   
M.W. Berry, J.J. Dongarra, and Y. Kim (1995). “A Parallel Algorithm for the Reduction of a Nonsymmetric Matrix to Block Upper Hessenberg Form,” Parallel Comput. 21, 1189–1211.   
G. Quintana-Orti and R. Van De Geijn (2006). “Improving the Performance of Reduction to Hessenberg Form,” ACM Trans. Math. Softw. 32, 180–194.   
S. Tomov, R. Nath, and J. Dongarra (2010). “Accelerating the Reduction to Upper Hessenberg, Tridiagonal, and Bidiagonal Forms Through Hybrid GPU-Based Computing,” Parallel Comput. 36, 645–654.   
L. Karlsson (2011). “Scheduling of Parallel Matrix Computations and Data Layout Conversion for HPC and Multicore Architectures,” PhD Thesis, University of Ume˚a.

Reaching the Hessenberg form via Gauss transforms is discussed in:

P. Businger (1969). “Reducing a Matrix to Hessenberg Form,” Math. Comput. 23, 819–821. G.W. Howell and N. Diaa (2005). “Algorithm 841: BHESS: Gaussian Reduction to a Similar Banded Hessenberg Form,” ACM Trans. Math. Softw. 31, 166–185.

Some interesting mathematical properties of the Hessenberg form may be found in:

B.N. Parlett (1967). “Canonical Decomposition of Hessenberg Matrices,” Math. Comput. 21, 223– 227.

Although the Hessenberg decomposition is largely appreciated as a “front end” decomposition for the QR iteration, it is increasingly popular as a cheap alternative to the more expensive Schur decomposition in certain problems. For a sampling of applications where it has proven to be very useful, consult:

W. Enright (1979). “On the Efficient and Reliable Numerical Solution of Large Linear Systems of O.D.E.’s,” IEEE Trans. Autom. Contr. AC-24, 905–908.   
G.H. Golub, S. Nash and C. Van Loan (1979). “A Hessenberg-Schur Method for the Problem AX + XB = C,” IEEE Trans. Autom. Contr. AC-24, 909–913.   
A. Laub (1981). “Efficient Multivariable Frequency Response Computations,” IEEE Trans. Autom. Contr. AC-26, 407–408.   
C.C. Paige (1981). “Properties of Numerical Algorithms Related to Computing Controllability,” IEEE Trans. Auto. Contr. AC-26, 130–138.   
G. Miminis and C.C. Paige (1982). “An Algorithm for Pole Assignment of Time Invariant Linear Systems,” Int. J. Contr. 35, 341–354.   
C. Van Loan (1982). “Using the Hessenberg Decomposition in Control Theory,” in Algorithms and Theory in Filtering and Control , D.C. Sorensen and R.J. Wets (eds.), Mathematical Programming Study No. 18, North Holland, Amsterdam, 102–111.   
C.D. Martin and C.F. Van Loan (2006). “Solving Real Linear Systems with the Complex Schur Decomposition,” SIAM J. Matrix Anal. Applic. 29, 177–183.

The advisability of posing polynomial root problems as companion matrix eigenvalue problem is discussed in:

A. Edelman and H. Murakami (1995). “Polynomial Roots from Companion Matrix Eigenvalues,” Math. Comput. 64, 763–776.
