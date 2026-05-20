# 8.3 The Symmetric QR Algorithm

The symmetric QR iteration (8.2.1) can be made more efficient in two ways. First, we show how to compute an orthogonal $U _ { 0 }$ such that $U _ { 0 } ^ { T } A U _ { 0 } = T$ is tridiagonal. With this reduction, the iterates produced by (8.2.1) are all tridiagonal and this reduces the work per step to $O ( n ^ { 2 } )$ . Second, the idea of shifts are introduced and with this change the convergence to diagonal form proceeds at a cubic rate. This is far better than having the off-diagonal entries going to to zero as $| \lambda _ { i + 1 } / \lambda _ { i } | ^ { k }$ as discussed in §8.2.5.

# 8.3.1 Reduction to Tridiagonal Form

If A is symmetric, then it is possible to find an orthogonal $Q$ such that

$$
Q ^ {T} A Q = T \tag {8.3.1}
$$

is tridiagonal. We call this the tridiagonal decomposition and as a compression of data, it represents a very big step toward diagonalization.

We show how to compute (8.3.1) with Householder matrices. Suppose that Householder matrices $P _ { 1 } , \ldots , P _ { k - 1 }$ have been determined such that if

$$
A _ {k - 1} = \left(P _ {1} \dots P _ {k - 1}\right) ^ {T} A \left(P _ {1} \dots P _ {k - 1}\right),
$$

then

$$
A _ {k - 1} = \left[ \begin{array}{c c c c} B _ {1 1} & B _ {1 2} & 0 \\ B _ {2 1} & B _ {2 2} & B _ {2 3} \\ 0 & B _ {3 2} & B _ {3 3} \\ k - 1 & 1 & n - k \end{array} \right] \begin{array}{c} k - 1 \\ 1 \\ n - k \end{array}
$$

is tridiagonal through its first k − 1 columns. If $\tilde { P } _ { k }$ is an order-(n −k) Householder matrix such that $\tilde { P } _ { k } B _ { 3 2 }$ is a multiple of $I _ { n - k } ( : , 1 )$ and if $P _ { k } = \mathrm { d i a g } ( I _ { k } , \tilde { P } _ { k } )$ , then the leading k-by-k principal submatrix of

$$
A _ {k} = P _ {k} A _ {k - 1} P _ {k} = \left[ \begin{array}{c c c} B _ {1 1} & B _ {1 2} & 0 \\ B _ {2 1} & B _ {2 2} & B _ {2 3} \tilde {P} _ {k} \\ 0 & \tilde {P} _ {k} B _ {3 2} & \tilde {P} _ {k} B _ {3 3} \tilde {P} _ {k} \end{array} \right] \begin{array}{c} k - 1 \\ 1 \\ n - k \end{array}
$$

is tridiagonal. Clearly, if $U _ { 0 } = P _ { 1 } \cdots P _ { n - 2 }$ , then $U _ { 0 } ^ { T } A U _ { 0 } = T$ is tridiagonal.

In the calculation of $A _ { k }$ it is important to exploit symmetry during the formation of the matrix $\tilde { P } _ { k } B _ { 3 3 } \tilde { P } _ { k }$ . To be specific, suppose that $\tilde { P } _ { k }$ has the form

$$
\tilde {P} _ {k} = I - \beta v v ^ {T}, \quad \beta = 2 / v ^ {T} v, \quad 0 \neq v \in \mathbb {R} ^ {n - k}.
$$

Note that if $p = \beta B _ { 3 3 } v$ and $w = p - ( \beta p ^ { T } v / 2 ) v$ , then

$$
\tilde {P} _ {k} B _ {3 3} \tilde {P} _ {k} = B _ {3 3} - v w ^ {T} - w v ^ {T}.
$$

Since only the upper triangular portion of this matrix needs to be calculated, we see that the transition from $A _ { k - 1 }$ to $A _ { k }$ can be accomplished in only $4 ( n - k ) ^ { 2 }$ flops.

Algorithm 8.3.1 (Householder Tridiagonalization) Given a symmetric $A \in \mathbb { R } ^ { n \times n }$ , the following algorithm overwrites A with $T = Q ^ { T } A { \dot { Q } }$ , where $T$ is tridiagonal and $Q =$ ${ \cal H } _ { 1 } \cdots { \cal H } _ { n - 2 }$ is the product of Householder transformations.

$$
\begin{array}{l} \text { for } k = 1: n - 2 \\ [ v, \beta ] = \operatorname{house} (A (k + 1: n, k)) \\ p = \beta A (k + 1: n, k + 1: n) v \\ w = p - (\beta p ^ {T} v / 2) v \\ A (k + 1, k) = \left\| A (k + 1: n, k) \right\| _ {2}; A (k, k + 1) = A (k + 1, k) \\ A (k + 1: n, k + 1: n) = A (k + 1: n, k + 1: n) - v w ^ {T} - w v ^ {T} \\ \end{array}
$$

end

This algorithm requires $4 n ^ { 3 } / 3$ flops when symmetry is exploited in calculating the rank-2 update. The matrix $Q$ can be stored in factored form in the subdiagonal portion of A. If $Q$ is explicitly required, then it can be formed with an additional $4 n ^ { 3 } / 3$ flops.

Note that if $T$ has a zero subdiagonal, then the eigenproblem splits into a pair of smaller eigenproblems. In particular, if $t _ { k + 1 , k } = 0$ , then

$$
\lambda (T) = \lambda (T (1: k, 1: k)) \cup \lambda (T (k + 1: n, k + 1: n)).
$$

If $T$ has no zero subdiagonal entries, then it is said to be unreduced.

Let $\hat { T }$ denote the computed version of $T$ obtained by Algorithm 8.3.1. It can be shown that $\hat { T } { = } \tilde { Q } ^ { T } ( A + \tilde { E } ) \tilde { Q }$ where $\tilde { Q }$ is exactly orthogonal and $E$ is a symmetric matrix satisfying $\left\| \ E \right\| _ { F } \leq c \mathbf { u } \left\| \ A \right\| _ { F }$ where c is a small constant. See Wilkinson (AEP, p. 297).

# 8.3.2 Properties of the Tridiagonal Decomposition

We prove two theorems about the tridiagonal decomposition both of which have key roles to play in the following. The first connects (8.3.1) to the QR factorization of a certain Krylov matrix. These matrices have the form

$$
K (A, v, k) = \left[ v \mid A v \mid \dots \mid A ^ {k - 1} v \right], \quad A \in \mathbb {R} ^ {n \times n}, v \in \mathbb {R} ^ {n}.
$$

Theorem 8.3.1. $I f Q ^ { T } A Q = T$ is the tridiagonal decomposition of the symmetric matrix $A \in \mathbb { R } ^ { n \times n }$ , then $Q ^ { T } K ( A , Q ( : , 1 ) , n ) = R$ is upper triangular. If R is nonsingular, then $T$ is unreduced. If R is singular and k is the smallest index so $r _ { k k } = 0$ , then $k$ is also the smallest index so $t _ { k , k - 1 }$ is zero. Compare with Theorem $\it { 7 . 4 . 3 } .$

Proof. It is clear that if $q _ { 1 } = Q ( : , 1 )$ , then

$$
\begin{array}{l} Q ^ {T} K (A, Q (:, 1), n) = \left[ Q ^ {T} q _ {1} \mid \left(Q ^ {T} A Q\right) \left(Q ^ {T} q _ {1}\right) \mid \dots \mid \left(Q ^ {T} A Q\right) ^ {n - 1} \left(Q ^ {T} q _ {1}\right) \right] \\ = \left[ e _ {1} \mid T e _ {1} \mid \dots \mid T ^ {n - 1} e _ {1} \right] = R \\ \end{array}
$$

is upper triangular with the property that $r _ { 1 1 } = 1$ and $r _ { i i } = t _ { 2 1 } t _ { 3 2 } \cdot \cdot \cdot t _ { i , i - 1 }$ for $i = 2 { : } n$ . Clearly, if $R$ is nonsingular, then $T$ is unreduced. If R is singular and $r _ { k k }$ is its first zero diagonal entry, then $k \geq 2$ and $t _ { k , k - 1 }$ is the first zero subdiagonal entry.

The next result shows that Q is essentially unique once $Q ( : , 1 )$ is specified.

Theorem 8.3.2 (Implicit Q Theorem). Suppose $Q \ = \ [ \ q _ { 1 } \ | \cdot \cdot \cdot | \ q _ { n } \ ]$ and $V =$ $\left[ \begin{array} { l } { v _ { 1 } } \end{array} \right| \cdots \left| \begin{array} { l } { v _ { n } } \end{array} \right]$ are orthogonal matrices with the property that both $Q ^ { T } A Q \ = \ T$ and $V ^ { T } A V = S$ are tridiagonal where $A \in \mathbb { R } ^ { n \times n }$ is symmetric. Let k denote the smallest positive integer for which $t _ { k + 1 , k } = 0$ , with the convention that $k = n$ if T is unreduced. $I f v _ { 1 } = q _ { 1 }$ , then $v _ { i } = \pm q _ { i }$ and $| t _ { i , i - 1 } | = | s _ { i , i - 1 } | f o r i = 2 { : } k$ . Moreover, if $k < n$ , then $s _ { k + 1 , k } = 0$ . Compare with Theorem $\it { 7 . 4 . 2 . }$

Proof. Define the orthogonal matrix $W = Q ^ { T } V$ and observe that $W ( : , 1 ) = I _ { n } ( : , 1 ) =$ $e _ { 1 }$ and $W ^ { T } T W = S$ . By Theorem 8.3.1, $W ^ { T } { \cdot } K ( T , e _ { 1 } , k )$ is upper triangular with full column rank. But $K ( T , e _ { 1 } , k )$ is upper triangular and so by the essential uniqueness of the thin QR factorization, $W ( : , 1 { : } k ) = I _ { n } ( : , 1 { : } k ) \cdot \mathrm { d i a g } ( \pm 1 , . . . , \pm 1 )$ . This says that $Q ( : , i ) = \pm V ( : , i )$ for $i = 1 { : } k$ . The comments about the subdiagonal entries follow since $t _ { i + 1 , i } = Q ( : , i + 1 ) ^ { T } A Q ( : , i )$ and $s _ { i + 1 , i } = V ( : , i + 1 ) ^ { T } A V ( : , i )$ for $i = 1 { : } n - 1$ .

# 8.3.3 The QR Iteration and Tridiagonal Matrices

We quickly state four facts that pertain to the QR iteration and tridiagonal matrices. Complete verifications are straightforward.

• Preservation of Form. If $T = Q R$ is the QR factorization of a symmetric tridiagonal matrix $T \in \mathbb { R } ^ { n \times n }$ , then Q has lower bandwidth 1 and R has upper bandwidth 2 and it follows that $T _ { + } = R Q = Q ^ { T } ( Q R ) Q = Q ^ { T } T Q$ is also symmetric and tridiagonal.

• $S h i f t s .$ . If $s \in \mathbb { R }$ and $T - s I = Q R$ is the QR factorization, then $T _ { + } = R Q + s I =$ $Q ^ { T } T Q$ is also tridiagonal. This is called a shifted QR step.

• Perfect $S h i f t s$ . If T is unreduced, then the first $n - 1$ columns of $T - s I$ are independent regardless of s. Thus, if $s \in \lambda ( T )$ and $Q R \ : = \ : T \ : - \ : s I$ is a QR factorization, then $r _ { n n } = 0$ and the last column of $T _ { + } = R Q + s I$ equals $s I _ { n } ( : , n ) =$ $s e _ { n }$ .

• Cost. If $T \in \mathbb { R } ^ { n \times n }$ is tridiagonal, then its QR factorization can be computed by applying a sequence of $n - 1$ Givens rotations:

$$
\begin{array}{l} \text { for } k = 1: n - 1 \\ [ c, s ] = \operatorname{givens} \left(t _ {k k}, t _ {k + 1, k}\right) \\ m = \min \{k + 2, n \} \\ T (k: k + 1, k: m) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} T (k: k + 1, k: m) \\ \end{array}
$$

end

This requires $O ( n )$ flops. If the rotations are accumulated, then $O ( n ^ { 2 } )$ flops are needed.

# 8.3.4 Explicit Single-Shift QR Iteration

If s is a good approximate eigenvalue, then we suspect that the $( n , n - 1 )$ will be small after a QR step with shift s. This is the philosophy behind the following iteration:

$$
T = U _ {0} ^ {T} A U _ {0} \quad \text {(tridiagonal)}
$$

for k = 0, 1, . . .

Determine real shift $\mu .$ (8.3.2)

$$
T - \mu I = U R \quad (\text { QR   factorization })
$$

$$
T = R U + \mu I
$$

end

If

$$
T = \left[ \begin{array}{c c c c c} a _ {1} & b _ {1} & & \dots & 0 \\ b _ {1} & a _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & b _ {n - 1} \\ 0 & \dots & & b _ {n - 1} & a _ {n} \end{array} \right],
$$

then one reasonable choice for the shift is $\mu = a _ { n }$ . However, a more effective choice is to shift by the eigenvalue of

$$
T (n - 1: n, n - 1: n) = \left[ \begin{array}{c c} a _ {n - 1} & b _ {n - 1} \\ b _ {n - 1} & a _ {n} \end{array} \right]
$$

that is closer to $a _ { n }$ . This is known as the Wilkinson shift and it is given by

$$
\mu = a _ {n} + d - \operatorname{sign} (d) \sqrt {d ^ {2} + b _ {n - 1} ^ {2}} \tag {8.3.3}
$$

where $\begin{array} { l l l } { d } & { = } & { ( a _ { n - 1 } - a _ { n } ) / 2 } \end{array}$ . Wilkinson (1968) has shown that (8.3.2) is cubically convergent with either shift strategy, but gives heuristic reasons why (8.3.3) is preferred.

# 8.3.5 Implicit Shift Version

It is possible to execute the transition from $T$ to $T _ { + } ~ = ~ R U + \mu I ~ = ~ U ^ { T } T U$ without explicitly forming the matrix $T - \mu I$ . This has advantages when the shift is much larger than some of the $a _ { i }$ . Let $c = \cos ( \theta )$ and $s = \sin ( \theta )$ be computed such that

$$
\left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} \left[ \begin{array}{c} a _ {1} - \mu \\ b _ {1} \end{array} \right] = \left[ \begin{array}{c} \times \\ 0 \end{array} \right].
$$

If we set $G _ { 1 } = G ( 1 , 2 , \theta )$ , then $G _ { 1 } e _ { 1 } = U e _ { 1 }$ and

$$
T \gets G _ {1} ^ {T} T G _ {1} = \left[ \begin{array}{c c c c c c} \times & \times & + & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 \\ + & \times & \times & \times & 0 & 0 \\ 0 & 0 & \times & \times & \times & 0 \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right].
$$

We are thus in a position to apply the implicit Q theorem provided we can compute rotations $G _ { 2 } , \ldots , G _ { n - 1 }$ with the property that if $Z ~ = ~ G _ { 1 } G _ { 2 } \cdot \cdot \cdot G _ { n - 1 }$ , then $Z e _ { 1 } =$ $\boldsymbol { G } _ { 1 } \boldsymbol { e } _ { 1 } = U \boldsymbol { e } _ { 1 }$ and $Z ^ { T } T Z$ is tridiagonal. Note that the first column of Z and U are identical provided we take each $G _ { i }$ to be of the form $G _ { i } = G ( i , i + 1 , \theta _ { i } ) , i = 2 { \cdot } n - 1$ . But $G _ { i }$ of this form can be used to chase the unwanted nonzero element $" + "$ out of the matrix $G _ { 1 } ^ { T } T G _ { 1 }$ as follows:

$$
\begin{array}{l} \xrightarrow {G _ {2}} \left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ \times & \times & \times & + & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 \\ 0 & + & \times & \times & \times & 0 \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right] \xrightarrow {G _ {3}} \left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 \\ 0 & \times & \times & \times & + & 0 \\ 0 & 0 & \times & \times & \times & 0 \\ 0 & 0 & + & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right] \\ \xrightarrow {G _ {4}} \left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 \\ 0 & 0 & \times & \times & \times & + \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & + & \times & \times \end{array} \right] \xrightarrow {G _ {5}} \left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 \\ 0 & 0 & \times & \times & \times & 0 \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right]. \\ \end{array}
$$

Thus, it follows from the implicit Q theorem that the tridiagonal matrix $Z ^ { T } T Z$ produced by this zero-chasing technique is essentially the same as the tridiagonal matrix $T$ obtained by the explicit method. (We may assume that all tridiagonal matrices in question are unreduced for otherwise the problem decouples.)

Note that at any stage of the zero-chasing, there is only one nonzero entry outside the tridiagonal band. How this nonzero entry moves down the matrix during the update $T  G _ { k } ^ { T } T G _ { k }$ is illustrated in the following:

$$
\left[ \begin{array}{c c c c} 1 & 0 & 0 & 0 \\ 0 & c & s & 0 \\ 0 & - s & c & 0 \\ 0 & 0 & 0 & 1 \end{array} \right] ^ {T} \left[ \begin{array}{c c c c} a _ {k} & b _ {k} & z _ {k} & 0 \\ b _ {k} & a _ {p} & b _ {p} & 0 \\ z _ {k} & b _ {p} & a _ {q} & b _ {q} \\ 0 & 0 & b _ {q} & a _ {r} \end{array} \right] \left[ \begin{array}{c c c c} 1 & 0 & 0 & 0 \\ 0 & c & s & 0 \\ 0 & - s & c & 0 \\ 0 & 0 & 0 & 1 \end{array} \right] = \left[ \begin{array}{c c c c} a _ {k} & b _ {k} & 0 & 0 \\ b _ {k} & a _ {p} & b _ {p} & z _ {p} \\ 0 & b _ {p} & a _ {q} & b _ {q} \\ 0 & z _ {p} & b _ {q} & a _ {r} \end{array} \right].
$$

Here $( p , q , r ) = ( k + 1 , k + 2 , k + 3 )$ . This update can be performed in about 26 flops once c and s have been determined from the equation $b _ { k } s + z _ { k } c = 0$ . Overall, we obtain

Algorithm 8.3.2 (Implicit Symmetric QR Step with Wilkinson Shift) Given an unreduced symmetric tridiagonal matrix $T \in \mathbb { R } ^ { n \times n }$ , the following algorithm overwrites $T$ with $Z ^ { T } T Z$ , where $Z = G _ { 1 } \cdot \cdot \cdot G _ { n - 1 }$ is a product of Givens rotations with the property that $Z ^ { T } ( T - \mu I )$ is upper triangular and $\mu$ is that eigenvalue of $T \mathrm { { s } }$ trailing 2-by-2 principal submatrix closer to $t _ { n n }$ .

$$
\begin{array}{l} d = \left(t _ {n - 1, n - 1} - t _ {n n}\right) / 2 \\ \mu = t _ {n n} - t _ {n, n - 1} ^ {2} / \left(d + \operatorname{sign} (d) \sqrt {d ^ {2} + t _ {n , n - 1} ^ {2}}\right) \\ x = t _ {1 1} - \mu \\ z = t _ {2 1} \\ \end{array}
$$

for $k = 1:n - 1$ [ [c,s] = \text{givens}(x,z) ] $T = G_k^T TG_k$ , where $G_k = G(k,k + 1,\theta)$ if $k < n - 1$ $x = t_{k+1,k}$ $z = t_{k+2,k}$ end  
end

This algorithm requires about 30n flops and n square roots. If a given orthogonal matrix Q is overwritten with $Q G _ { 1 } \cdots G _ { n - 1 }$ , then an additional $6 n ^ { 2 }$ flops are needed. Of course, in any practical implementation the tridiagonal matrix T would be stored in a pair of n-vectors and not in an n-by-n array.

Algorithm 8.3.2 is the basis of the symmetric QR algorithm—the standard means for computing the Schur decomposition of a dense symmetric matrix.

Algorithm 8.3.3 (Symmetric QR Algorithm) Given $A \in \mathbb { R } ^ { n \times n }$ (symmetric) and a tolerance tol greater than the unit roundoff, this algorithm computes an approximate symmetric Schur decomposition $Q ^ { T } A Q = D$ . A is overwritten with the tridiagonal decomposition.

Use Algorithm 8.3.1, compute the tridiagonalization

$$
T = (P _ {1} \dots P _ {n - 2}) ^ {T} A (P _ {1} \dots P _ {n - 2})
$$

Set $D = T$ and if Q is desired, form $Q = P _ { 1 } \cdot \cdot \cdot P _ { n - 2 }$ . (See §5.1.6.)

$$
\text { until } q = n
$$

For $i = 1 { : } n - 1$ , set $d _ { i + 1 , i }$ and $d _ { i , i + 1 }$ to zero if

$$
\left| d _ {i + 1, i} \right| = \left| d _ {i, i + 1} \right| \leq \operatorname{tol} \left(\left| d _ {i i} \right| + \left| d _ {i + 1, i + 1} \right|\right)
$$

Find the largest q and the smallest p such that if

$$
D = \left[ \begin{array}{c c c} D _ {1 1} & 0 & 0 \\ 0 & D _ {2 2} & 0 \\ 0 & 0 & D _ {3 3} \end{array} \right] _ { \begin{array}{c} p \\ n - p - q \\ q \end{array} }
$$

then $D _ { 3 3 }$ is diagonal and $D _ { 2 2 }$ is unreduced.

$$
\text { if } q <   n
$$

Apply Algorithm 8.3.2 to $D _ { 2 2 }$

$$
D = \mathrm{diag} (I _ {p}, Z, I _ {q}) ^ {T} \cdot D \cdot \mathrm{diag} (I _ {p}, Z, I _ {q})
$$

If Q is desired, then $Q = Q \cdot \mathrm { d i a g } ( I _ { p } , Z , I _ { q } )$ .

end end

This algorithm requires about $4 n ^ { 3 } / 3$ flops if Q is not accumulated and about $9 n ^ { 3 }$ flops if Q is accumulated.

The computed eigenvalues $\hat { \lambda } _ { i }$ obtained via Algorithm 8.3.3 are the exact eigenvalues of a matrix that is near to A:

$$
Q _ {0} ^ {T} (A + E) Q _ {0} = \operatorname{diag} (\hat {\lambda} _ {i}), \quad Q _ {0} ^ {T} Q _ {0} = I, \quad \| E \| _ {2} \approx \mathbf {u} \| A \| _ {2}.
$$

Using Corollary 8.1.6 we know that the absolute error in each $\hat { \lambda } _ { i }$ is small in the sense that

$$
\left| \hat {\lambda} _ {i} - \lambda_ {i} \right| \approx \mathbf {u} \| A \| _ {2}.
$$

If $\hat { Q } = \left[ \begin{array} { l } { \hat { q } _ { 1 } \vert \cdots \vert \hat { q } _ { n } } \end{array} \right]$ is the computed matrix of orthonormal eigenvectors, then the accuracy of $\hat { q } _ { i }$ depends on the separation of $\lambda _ { i }$ from the remainder of the spectrum. See Theorem 8.1.12.

If all of the eigenvalues and a few of the eigenvectors are desired, then it is cheaper not to accumulate $Q$ in Algorithm 8.3.3. Instead, the desired eigenvectors can be found via inverse iteration with T . See §8.2.2. Usually just one step is sufficient to get a good eigenvector, even with a random initial vector.

If just a few eigenvalues and eigenvectors are required, then the special techniques in §8.4 are appropriate.

# 8.3.6 The Rayleigh Quotient Connection

It is interesting to identify a relationship between the Rayleigh quotient iteration and the symmetric QR algorithm. Suppose we apply the latter to the tridiagonal matrix $T \in \mathbb { R } ^ { n \times n }$ with shift $\sigma = e _ { n } ^ { T } T e _ { n } = t _ { n n }$ . If $T - \sigma I = Q R$ , then we obtain $T _ { + } = R Q + \sigma I$ . From the equation $( T - \sigma I ) Q = R ^ { T }$ it follows that

$$
(T - \sigma I) q _ {n} = r _ {n n} e _ {n},
$$

where $q _ { n }$ is the last column of the orthogonal matrix $Q .$ . Thus, if we apply (8.2.6) with $x _ { 0 } = e _ { n }$ , then $x _ { 1 } = q _ { n }$ .

# 8.3.7 Orthogonal Iteration with Ritz Acceleration

Recall from §8.2.4 that an orthogonal iteration step involves a matrix-matrix product and a QR factorization:

$$
Z _ {k} = A \tilde {Q} _ {k - 1},
$$

$$
\tilde {Q} _ {k} R _ {k} = Z _ {k} \quad (\text { QR   factorization })
$$

Theorem 8.1.14 says that we can minimize $\parallel A \tilde { Q } _ { k } - \tilde { Q } _ { k } S \parallel _ { _ F }$ by setting S equal to

$$
S _ {k} = \tilde {Q} _ {k} ^ {T} A \tilde {Q} _ {k}.
$$

If $U _ { k } ^ { T } S _ { k } U _ { k } = D _ { k }$ is the Schur decomposition of $S _ { k } \in \mathbb { R } ^ { r \times r }$ and $Q _ { k } = \tilde { Q } _ { k } U _ { k }$ , then

$$
\left\| A Q _ {k} - Q _ {k} D _ {k} \right\| _ {F} = \left\| A \tilde {Q} _ {k} - \tilde {Q} _ {k} S _ {k} \right\| _ {F}
$$

showing that the columns of $Q _ { k }$ are the best possible basis to take after k steps from the standpoint of minimizing the residual. This defines the Ritz acceleration idea:

$$
Q _ {0} \in \mathbb {R} ^ {n \times r} \text {   given   with   } Q _ {0} ^ {T} Q _ {0} = I _ {r}
$$

for $k = 1 , 2 , \dots$

$$
Z _ {k} = A Q _ {k - 1}
$$

$$
\tilde {Q} _ {k} R _ {k} = Z _ {k} \quad \text {(QR factorization)}
$$

$$
S _ {k} = \tilde {Q} _ {k} ^ {T} A \tilde {Q} _ {k} \tag {8.3.6}
$$

$$
U _ {k} ^ {T} S _ {k} U _ {k} = D _ {k} \quad \text {(Schur decomposition)}
$$

$$
Q _ {k} = \tilde {Q} _ {k} U _ {k}
$$

end

It can be shown that if

$$
D _ {k} = \operatorname{diag} \left(\theta_ {1} ^ {(k)}, \dots , \theta_ {r} ^ {(k)}\right) ], \quad | \theta_ {1} ^ {(k)} | \geq \dots \geq | \theta_ {r} ^ {(k)} |,
$$

then

$$
\left| \theta_ {i} ^ {(k)} - \lambda_ {i} (A) \right| = O \left(\left| \frac {\lambda_ {r + 1}}{\lambda_ {i}} \right| ^ {k}\right), \quad i = 1: r.
$$

Recall that Theorem 8.2.2 says the eigenvalues of $\tilde { Q } _ { k } ^ { T } A \tilde { Q } _ { k }$ converge with rate $| \lambda _ { r + 1 } / \lambda _ { r } | ^ { k }$ . Thus, the Ritz values converge at a more favorable rate. For details, see Stewart (1969).

# Problems

P8.3.1 Suppose λ is an eigenvalue of a symmetric tridiagonal matrix T . Show that if λ has algebraic multiplicity k, then at least k − 1 of T ’s subdiagonal elements are zero.

P8.3.2 Suppose A is symmetric and has bandwidth p. Show that if we perform the shifted QR step $A - \mu I = Q R , A = R Q + \mu I$ , then A has bandwidth p.

P8.3.3 Let

$$
A = \left[ \begin{array}{c c} w & x \\ x & z \end{array} \right]
$$

be real and suppose we perform the following shifted QR step: $A - z I = U R , \tilde { A } = R U + z I$ . Show that

$$
\tilde {A} = \left[ \begin{array}{c c} \tilde {w} & \tilde {x} \\ \tilde {x} & \tilde {z} \end{array} \right]
$$

where

$$
\tilde {w} = w + x ^ {2} (w - z) / [ (w - z) ^ {2} + x ^ {2} ],
$$

$$
\tilde {z} = z - x ^ {2} (w - z) / [ (w - z) ^ {2} + x ^ {2} ],
$$

$$
\tilde {x} = - x ^ {3} / [ (w - z) ^ {2} + x ^ {2} ].
$$

P8.3.4 Suppose $A \in \mathbb { C } ^ { n \times n }$ is Hermitian. Show how to construct unitary Q such that $Q ^ { H } A Q = T$ is real, symmetric, and tridiagonal.

P8.3.5 Show that if $A = B + i C$ is Hermitian, then

$$
M = \left[ \begin{array}{c c} B & - C \\ C & B \end{array} \right]
$$

is symmetric. Relate the eigenvalues and eigenvectors of A and M.

P8.3.6 Rewrite Algorithm 8.3.2 for the case when A is stored in two n-vectors. Justify the given flop count.

P8.3.7 Suppose $\boldsymbol { A } = \boldsymbol { S } + \sigma u \boldsymbol { u } ^ { T }$ where $S \in \mathbb { R } ^ { n \times n }$ is skew-symmetric $( S ^ { T } = - S ) , u \in \mathbb { R } ^ { n }$ has unit

2-norm, and $\sigma \in \mathbb { R }$ . Show how to compute an orthogonal Q such that $Q ^ { T } A Q$ is tridiagonal and $Q ^ { T } u = e _ { 1 }$ .

P8.3.8 Suppose

$$
C = \left[ \begin{array}{c c} 0 & B ^ {T} \\ B & 0 \end{array} \right]
$$

where $B \in \mathbb { R } ^ { n \times n }$ is upper bidiagonal. Determine a perfect shuffle permutation $P \in \mathbb { R } ^ { 2 n \times 2 n }$ so that $T = P C P ^ { T }$ is tridiagonal with a zero diagonal.

# Notes and References for §8.3

Historically important Algol specifications related to the algorithms in this section include:

R.S. Martin and J.H. Wilkinson (1967). “Solution of Symmetric and Unsymmetric Band Equations and the Calculation of Eigenvectors of Band Matrices,” Numer. Math. 9, 279–301.   
H. Bowdler, R.S. Martin, C. Reinsch, and J.H. Wilkinson (1968). “The QR and QL Algorithms for Symmetric Matrices,” Numer. Math. 11, 293–306.   
A. Dubrulle, R.S. Martin, and J.H. Wilkinson (1968). “The Implicit QL Algorithm,” Numer. Math. 12, 377–383.   
R.S. Martin and J.H. Wilkinson (1968). “Householder’s Tridiagonalization of a Symmetric Matrix,” Numer. Math. 11, 181–195.   
C. Reinsch and F.L. Bauer (1968). “Rational QR Transformation with Newton’s Shift for Symmetric Tridiagonal Matrices,” Numer. Math. 11, 264–272.   
R.S. Martin, C. Reinsch, and J.H. Wilkinson (1970). “The QR Algorithm for Band Symmetric Matrices,” Numer. Math. 16, 85–92.

The convergence properties of Algorithm 8.3.3 are detailed in Lawson and Hanson (SLE), see:

J.H. Wilkinson (1968). “Global Convergence of Tridiagonal QR Algorithm With Origin Shifts,” Lin. Alg. Applic. 1, 409–420.   
T.J. Dekker and J.F. Traub (1971). “The Shifted QR Algorithm for Hermitian Matrices,” Lin. Alg. Applic. 4, 137–154.   
W. Hoffman and B.N. Parlett (1978). “A New Proof of Global Convergence for the Tridiagonal QL Algorithm,” SIAM J. Numer. Anal. 15, 929–937.   
S. Batterson (1994). “Convergence of the Francis Shifted QR Algorithm on Normal Matrices,” Lin. Alg. Applic. 207, 181–195.   
T.-L. Wang (2001). “Convergence of the Tridiagonal QR Algorithm,” Lin. Alg. Applic. 322, 1–17.

Shifting and deflation are critical to the effective implementation of the symmetric QR iteration, see:

F.L. Bauer and C. Reinsch (1968). “Rational QR Transformations with Newton Shift for Symmetric Tridiagonal Matrices,” Numer. Math. 11, 264–272.

G.W. Stewart (1970). “Incorporating Origin Shifts into the QR Algorithm for Symmetric Tridiagonal Matrices,” Commun. ACM 13, 365–367.

I.S. Dhillon and A.N. Malyshev (2003). “Inner Deflation for Symmetric Tridiagonal Matrices,” Lin. Alg. Applic. 358, 139–144.

The efficient reduction of a general band symmetric matrix to tridiagonal form is a challenging computation from several standpoints:

H.R. Schwartz (1968). “Tridiagonalization of a Symmetric Band Matrix,” Numer. Math. 12, 231–241. C.H. Bischof and X. Sun (1996). “On Tridiagonalizing and Diagonalizing Symmetric Matrices with Repeated Eigenvalues,” SIAM J. Matrix Anal. Applic. 17, 869–885.

L. Kaufman (2000). “Band Reduction Algorithms Revisited,” ACM Trans. Math. Softw. 26, 551–567.

C.H. Bischof, B. Lang, and X. Sun (2000). “A Framework for Symmetric Band Reduction,” ACM Trans. Math. Softw. 26, 581–601.

Finally we mention that comparable techniques exist for skew-symmetric and general normal matrices, see:

R.C. Ward and L.J. Gray (1978). “Eigensystem Computation for Skew-Symmetric and A Class of Symmetric Matrices,” ACM Trans. Math. Softw. 4, 278–285.

C.P. Huang (1981). “On the Convergence of the QR Algorithm with Origin Shifts for Normal Matrices,” IMA J. Numer. Anal. 1, 127–133.

S. Iwata (1998). “Block Triangularization of Skew-Symmetric Matrices,” Lin. Alg. Applic. 273, 215–226.
