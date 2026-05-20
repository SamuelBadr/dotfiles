# 7.5 The Practical QR Algorithm

We return to the Hessenberg QR iteration, which we write as follows:

$$
H = U _ {0} ^ {T} A U _ {0} \quad \text {(Hessenberg reduction)}
$$

for k = 1, 2, . . .

$$
H = U R \quad (\text { QR   factorization }) \tag {7.5.1}
$$

$$
H = R U
$$

end

Our aim in this section is to describe how the H’s converge to upper quasi-triangular form and to show how the convergence rate can be accelerated by incorporating shifts.

# 7.5.1 Deflation

Without loss of generality we may assume that each Hessenberg matrix H in (7.5.1) is unreduced. If not, then at some stage we have

$$
H = \left[ \begin{array}{c c} H _ {1 1} & H _ {1 2} \\ 0 & H _ {2 2} \end{array} \right] _ {n - p} ^ {p}
$$

where $1 \leq p < n$ and the problem decouples into two smaller problems involving $H _ { 1 1 }$ and $H _ { 2 2 }$ . The term deflation is also used in this context, usually when $p = n - 1$ o r $n - 2$ .

In practice, decoupling occurs whenever a subdiagonal entry in H is suitably small. For example, if

$$
\left| h _ {p + 1, p} \right| \leq c \mathbf {u} \left(\left| h _ {p p} \right| + \left| h _ {p + 1, p + 1} \right|\right) \tag {7.5.2}
$$

for a small constant c, then $h _ { p + 1 , p }$ can justifiably be set to zero because rounding errors of order u H  are typically present throughout the matrix anyway.

# 7.5.2 The Shifted QR Iteration

Let $\mu \in \mathbb { R }$ and consider the iteration:

$$
H = U _ {0} ^ {T} A U _ {0} \quad \text {(Hessenberg reduction)}
$$

for k = 1, 2, . . .

Determine a scalar $\mu .$

$$
H - \mu I = U R \quad (\text { QR   factorization }) \tag {7.5.3}
$$

$$
H = R U + \mu I
$$

end

The scalar $\mu$ is referred to as a shift . Each matrix H generated in (7.5.3) is similar to A, since

$$
R U + \mu I = U ^ {T} (U R + \mu I) U = U ^ {T} H U.
$$

If we order the eigenvalues $\lambda _ { i }$ of A so that

$$
\left| \lambda_ {1} - \mu \right| \geq \dots \geq \left| \lambda_ {n} - \mu \right|,
$$

and $\mu$ is fixed from iteration to iteration, then the theory of §7.3 says that the pth subdiagonal entry in H converges to zero with rate

$$
\left| \frac {\lambda_ {p + 1} - \mu}{\lambda_ {p} - \mu} \right| ^ {k}.
$$

Of course, if $\lambda _ { p } = \lambda _ { p + 1 }$ , then there is no convergence at all. But if, for example, µ is much closer to $\lambda _ { n }$ than to the other eigenvalues, then the zeroing of the $( n , n - 1 )$ entry is rapid. In the extreme case we have the following:

Theorem 7.5.1. Let $\mu$ be an eigenvalue of an $n { - } b y { - } n$ unreduced Hessenberg matrix H . If

$$
\tilde {H} = R U + \mu I,
$$

where $H - \mu I = U R$ is the QR factorization of $H - \mu I$ , then $\tilde { h } _ { n , n - 1 } = 0$ and $\tilde { h } _ { n n } = \mu$ .

Proof. Since H is an unreduced Hessenberg matrix the first n − 1 columns of $H - \mu I$ are independent, regardless of µ. Thus, if $U R = \left( H - \mu I \right)$ is the QR factorization then $r _ { i i } \neq 0$ for $i = 1 { : } n - 1$ . But if $H - \mu I$ is singular, then $r _ { 1 1 } \cdot \cdot \cdot r _ { n n } = 0$ . Thus, $r _ { n n } = 0$ and $\tilde { H } ( n , : ) = [ 0 , . . . , 0 , \mu ]$ .

The theorem says that if we shift by an exact eigenvalue, then in exact arithmetic deflation occurs in one step.

# 7.5.3 The Single-Shift Strategy

Now let us consider varying µ from iteration to iteration incorporating new information about $\lambda ( A )$ as the subdiagonal entries converge to zero. A good heuristic is to regard $h _ { n n }$ as the best approximate eigenvalue along the diagonal. If we shift by this quantity during each iteration, we obtain the single-shift QR iteration:

for $k = 1 , 2 , \dots$

$$
\mu = H (n, n)
$$

$$
H - \mu I = U R \quad (\text { QR   factorization }) \tag {7.5.4}
$$

$$
H = R U + \mu I
$$

end

If the $( n , n - 1 )$ entry converges to zero, it is likely to do so at a quadratic rate. To see this, we borrow an example from Stewart (IMC, p. 366). Suppose H is an unreduced upper Hessenberg matrix of the form

$$
H = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \epsilon & h _ {n n} \end{array} \right]
$$

and that we perform one step of the single-shift QR algorithm, i.e.,

$$
U R = H - h _ {n n}
$$

$$
\tilde {H} = R U + h _ {n n} I.
$$

After $n - 2$ steps in the orthogonal reduction of $H - h _ { n n } I$ to upper triangular form we obtain a matrix with the following structure:

$$
H = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & a & b \\ 0 & 0 & 0 & \epsilon & 0 \end{array} \right].
$$

It is not hard to show that

$$
\tilde {h} _ {n, n - 1} = - \frac {\epsilon^ {2} b}{a ^ {2} + \epsilon^ {2}}.
$$

If we assume that $\epsilon \ll a$ , then it is clear that the new $( n , n - 1 )$ entry has order $\epsilon ^ { 2 }$ , precisely what we would expect of a quadratically converging algorithm.

# 7.5.4 The Double-Shift Strategy

Unfortunately, difficulties with (7.5.4) can be expected if at some stage the eigenvalues $a _ { 1 }$ and $a _ { 2 }$ of

$$
G = \left[ \begin{array}{c c} h _ {m m} & h _ {m n} \\ h _ {n m} & h _ {n n} \end{array} \right], \qquad m = n - 1, \tag {7.5.5}
$$

are complex for then $h _ { n n }$ would tend to be a poor approximate eigenvalue.

A way around this difficulty is to perform two single-shift QR steps in succession using $a _ { 1 }$ and $a _ { 2 }$ as shifts:

$$
H - a _ {1} I = U _ {1} R _ {1}
$$

$$
H _ {1} = R _ {1} U _ {1} + a _ {1} I \tag {7.5.6}
$$

$$
H _ {1} - a _ {2} I = U _ {2} R _ {2}
$$

$$
H _ {2} = R _ {2} U _ {2} + a _ {2} I
$$

These equations can be manipulated to show that

$$
(U _ {1} U _ {2}) (R _ {2} R _ {1}) = M \tag {7.5.7}
$$

where M is defined by

$$
M = (H - a _ {1} I) (H - a _ {2} I). \tag {7.5.8}
$$

Note that M is a real matrix even if G’s eigenvalues are complex since

$$
M = H ^ {2} - s H + t I
$$

where

$$
s = a _ {1} + a _ {2} = h _ {m m} + h _ {n n} = \operatorname{tr} (G) \in \mathbb {R}
$$

and

$$
t = a _ {1} a _ {2} = h _ {m m} h _ {n n} - h _ {m n} h _ {n m} = \det (G) \in \mathbb {R}.
$$

Thus, (7.5.7) is the QR factorization of a real matrix and we may choose $U _ { 1 }$ and $U _ { 2 }$ so that $Z = U _ { 1 } U _ { 2 }$ is real orthogonal. It then follows that

$$
H _ {2} = U _ {2} ^ {H} H _ {1} U _ {2} = U _ {2} ^ {H} (U _ {1} ^ {H} H U _ {1}) U _ {2} = (U _ {1} U _ {2}) ^ {H} H (U _ {1} U _ {2}) = Z ^ {T} H Z
$$

is real.

Unfortunately, roundoff error almost always prevents an exact return to the real field. A real $H _ { 2 }$ could be guaranteed if we

• explicitly form the real matrix $M = H ^ { 2 } - s H + t I$ ,   
• compute the real QR factorization $M = Z R .$ , and   
• set $H _ { 2 } = Z ^ { T } H Z .$ .

But since the first of these steps requires $O ( n ^ { 3 } )$ flops, this is not a practical course of action.

# 7.5.5 The Double-Implicit-Shift Strategy

Fortunately, it turns out that we can implement the double-shift step with $O ( n ^ { 2 } )$ flops by appealing to the implicit Q theorem of §7.4.5. In particular we can effect the transition from H to $H _ { 2 }$ in $O ( n ^ { 2 } )$ flops if we

• compute $M e _ { 1 }$ , the first column of $M$   
• determine a Householder matrix $P _ { 0 }$ such that $P _ { 0 } ( M e _ { 1 } )$ is a multiple of $e _ { 1 }$   
• compute Householder matrices $P _ { 1 } , \ldots , P _ { n - 2 }$ such that if

$$
Z _ {1} = P _ {0} P _ {1} \dots P _ {n - 2},
$$

then $Z _ { 1 } ^ { T } H Z _ { 1 }$ is upper Hessenberg and the first columns of Z and $Z _ { 1 }$ are the same.

Under these circumstances, the implicit Q theorem permits us to conclude that, if $Z ^ { T } H Z$ and $Z _ { 1 } ^ { T } H Z _ { 1 }$ are both unreduced upper Hessenberg matrices, then they are essentially equal. Note that if these Hessenberg matrices are not unreduced, then we can effect a decoupling and proceed with smaller unreduced subproblems.

Let us work out the details. Observe first that $P _ { 0 }$ can be determined in $O ( 1 )$ flops since $M e _ { 1 } = [ x , y , z , 0 , \ldots , 0 ] ^ { T }$ where

$$
\begin{array}{l} x = h _ {1 1} ^ {2} + h _ {1 2} h _ {2 1} - s h _ {1 1} + t, \\ y = h _ {2 1} (h _ {1 1} + h _ {2 2} - s), \\ z = h _ {2 1} h _ {3 2}. \\ \end{array}
$$

Since a similarity transformation with $P _ { 0 }$ only changes rows and columns 1, 2, and 3, we see that

$$
P _ {0} H P _ {0} = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right].
$$

Now the mission of the Householder matrices $P _ { 1 } , \ldots , P _ { n - 2 }$ is to restore this matrix to upper Hessenberg form. The calculation proceeds as follows:

$$
\left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right] \xrightarrow {P _ {1}} \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right] \xrightarrow {P _ {2}}
$$

$$
\left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\end{array}\right] \overset {P _ {3}} {\rightarrow} \left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&0&\times&\times&\times\\0&0&0&\times&\times&\times\end{array}\right] \overset {P _ {4}} {\rightarrow} \left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&0&\times&\times&\times\\0&0&0&0&\times&\times\end{array}\right].
$$

Each $P _ { k }$ is the identity with a 3-by-3 or 2-by-2 Householder somewhere along its diagonal, e.g.,

$$
\begin{array}{l} P _ {1} = \left[ \begin{array}{c c c c c c} 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 \end{array} \right], \quad P _ {2} = \left[ \begin{array}{c c c c c c} 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & \times & \times & \times & 0 \\ 0 & 0 & \times & \times & \times & 0 \\ 0 & 0 & \times & \times & \times & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 \end{array} \right], \\ P _ {3} = \left[ \begin{array}{c c c c c c} 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \end{array} \right], \quad P _ {4} = \left[ \begin{array}{c c c c c c} 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right]. \\ \end{array}
$$

The applicability of Theorem 7.4.3 (the implicit Q theorem) follows from the observation that $P _ { k } e _ { 1 } = e _ { 1 }$ for $k = 1 { : } n - 2$ and that $P _ { 0 }$ and $Z$ have the same first column. Hence, $Z _ { 1 } e _ { 1 } = Z e _ { 1 }$ , and we can assert that $Z _ { 1 }$ essentially equals $Z$ provided that the upper Hessenberg matrices $Z ^ { T } H Z$ and $Z _ { 1 } ^ { T } H Z _ { 1 }$ are each unreduced.

The implicit determination of $H _ { 2 }$ from H outlined above was first described by Francis (1961) and we refer to it as a Francis QR step. The complete Francis step is summarized as follows:

Algorithm 7.5.1 (Francis QR step) Given the unreduced upper Hessenberg matrix $H \in \mathbb { R } ^ { n \times n }$ whose trailing 2-by-2 principal submatrix has eigenvalues $a _ { 1 }$ and $a _ { 2 }$ , this algorithm overwrites H with $Z ^ { T } H Z$ , where $Z$ is a product of Householder matrices and $Z ^ { T } ( H - a _ { 1 } I ) ( H - a _ { 2 } I )$ is upper triangular.

$m = n - 1$ {Compute first column of $(H - a_{1}I)(H - a_{2}I)\}$ $s = H(m,m) + H(n,n)$ $t = H(m,m)\cdot H(n,n) - H(m,n)\cdot H(n,m)$ $x = H(1,1)\cdot H(1,1) + H(1,2)\cdot H(2,1) - s\cdot H(1,1) + t$ $y = H(2,1)\cdot (H(1,1) + H(2,2) - s)$ $z = H(2,1)\cdot H(3,2)$ for $k = 0:n - 3$ $[v,\beta ] = \mathsf{house}([xyz]^T)$ $q = \max \{1,k\}$ $H(k + 1:k + 3,q:n) = (I - \beta vv^T)\cdot H(k + 1:k + 3,q:n)$ $r = \min \{k + 4,n\}$ $H(1:r,k + 1:k + 3) = H(1:r,k + 1:k + 3)\cdot (I - \beta vv^T)$ $x = H(k + 2,k + 1)$ $y = H(k + 3,k + 1)$ if $k <   n - 3$ $z = H(k + 4,k + 1)$ end

end

$[v, \beta] = \text{house}([x y]^T)$ $H(n - 1:n, n - 2:n) = (I - \beta vv^T) \cdot H(n - 1:n, n - 2:n)$ $H(1:n, n - 1:n) = H(1:n, n - 1:n) \cdot (I - \beta vv^T)$

This algorithm requires $1 0 n ^ { 2 }$ flops. If Z is accumulated into a given orthogonal matrix, an additional $1 0 n ^ { 2 }$ flops are necessary.

# 7.5.6 The Overall Process

Reduction of A to Hessenberg form using Algorithm 7.4.2 and then iteration with Algorithm 7.5.1 to produce the real Schur form is the standard means by which the dense unsymmetric eigenproblem is solved. During the iteration it is necessary to monitor the subdiagonal elements in H in order to spot any possible decoupling. How this is done is illustrated in the following algorithm:

Algorithm 7.5.2 (QR Algorithm) Given $A \in \mathbb { R } ^ { n \times n }$ and a tolerance tol greater than the unit roundoff, this algorithm computes the real Schur canonical form ${ \overset { \cdot } { Q } } ^ { T } A Q = T$ . If $Q$ and $T$ are desired, then $T$ is stored in H. If only the eigenvalues are desired, then diagonal blocks in $T$ are stored in the corresponding positions in H.

Use Algorithm 7.4.2 to compute the Hessenberg reduction

$$
H = U _ {0} ^ {T} A U _ {0} \text { where } U _ {0} = P _ {1} \dots P _ {n - 2}.
$$

If $Q$ is desired form $Q = P _ { 1 } \cdots P _ { n - 2 }$ . (See §5.1.6.)

$$
\text { until } q = n
$$

Set to zero all subdiagonal elements that satisfy:

$$
\left| h _ {i, i - 1} \right| \leq \operatorname{tol} \cdot \left(\left| h _ {i i} \right| + \left| h _ {i - 1, i - 1} \right|\right).
$$

Find the largest nonnegative q and the smallest non-negative p such that

$$
H = \left[ \begin{array}{c c c} H _ {1 1} & H _ {1 2} & H _ {1 3} \\ 0 & H _ {2 2} & H _ {2 3} \\ 0 & 0 & H _ {3 3} \end{array} \right] \begin{array}{c} p \\ n - p - q \\ q \end{array}
$$

where $H _ { 3 3 }$ is upper quasi-triangular and $H _ { 2 2 }$ is unreduced.

$$
\text { if } q <   n
$$

Perform a Francis QR step on $H _ { 2 2 } \colon H _ { 2 2 } = Z ^ { T } H _ { 2 2 } Z .$

if Q is required

$$
Q = Q \cdot \operatorname{diag} \left(I _ {p}, Z, I _ {q}\right)
$$

$$
H _ {1 2} = H _ {1 2} Z
$$

$$
H _ {2 3} = Z ^ {T} H _ {2 3}
$$

end

end

end

Upper triangularize all 2-by-2 diagonal blocks in H that have real eigenvalues and accumulate the transformations (if necessary).

This algorithm requires $2 5 n ^ { 3 }$ flops if $Q$ and $T$ are computed. If only the eigenvalues are desired, then $1 0 n ^ { 3 }$ flops are necessary. These flops counts are very approximate and are based on the empirical observation that on average only two Francis iterations are required before the lower 1-by-1 or 2-by-2 decouples.

The roundoff properties of the QR algorithm are what one would expect of any orthogonal matrix technique. The computed real Schur form $\hat { T }$ is orthogonally similar to a matrix near to A, i.e.,

$$
Q ^ {T} (A + E) Q = \hat {T}
$$

where $Q ^ { T } Q = I$ and $\parallel E \parallel _ { 2 } \approx \mathbf { u } \parallel A \parallel _ { 2 }$ . The computed $\hat { Q }$ is almost orthogonal in the sense that ${ \hat { Q } } ^ { T } { \hat { Q } } = I + { \ddot { F } }$ where $\parallel F \parallel _ { 2 } \approx \mathbf { u }$ .

The order of the eigenvalues along $\hat { T }$ is somewhat arbitrary. But as we discuss in $\ S 7 . 6$ , any ordering can be achieved by using a simple procedure for swapping two adjacent diagonal entries.

# 7.5.7 Balancing

Finally, we mention that if the elements of A have widely varying magnitudes, then A should be balanced before applying the QR algorithm. This is an $O ( n ^ { 2 } )$ calculation in which a diagonal matrix D is computed so that if

$$
D ^ {- 1} A D = \left[ c _ {1} \mid \dots \mid c _ {n} \right] = \left[ \begin{array}{c} r _ {1} ^ {T} \\ \vdots \\ r _ {n} ^ {T} \end{array} \right]
$$

then $\parallel r _ { i } \parallel _ { \infty } \approx \parallel c _ { i } \parallel _ { \infty }$ for $i = 1 { : } n$ . The diagonal matrix D is chosen to have the form

$$
D = \mathrm{diag} (\beta^ {i _ {1}}, \ldots , \beta^ {i _ {n}})
$$

where $\beta$ is the floating point base. Note that $D ^ { - 1 } A D$ can be calculated without roundoff. When A is balanced, the computed eigenvalues are usually more accurate although there are exceptions. See Parlett and Reinsch (1969) and Watkins(2006).

# Problems

P7.5.1 Show that if $\bar { H } = Q ^ { T } H Q$ is obtained by performing a single-shift QR step with

$$
H = \left[ \begin{array}{c c} w & x \\ y & z \end{array} \right],
$$

then $| \bar { h } _ { 2 1 } | \le | y ^ { 2 } x | / [ ( w - z ) ^ { 2 } + y ^ { 2 } ] .$ .

P7.5.2 Given $A \in \mathbb { R } ^ { 2 \times 2 } ,$ , show how to compute a diagonal $D \in \mathbb { R } ^ { 2 \times 2 }$ so that $| \ D ^ { - 1 } A D \ | | _ { F }$ is minimized.

P7.5.3 Explain how the single-shift QR step H − $\mathbf { \nabla } \cdot \mu I = U R , \tilde { H } = R U + \mu I$ can be carried out implicitly. That is, show how the transition from H to H˜ can be carried out without subtracting the shift $\mu$ from the diagonal of H.

P7.5.4 Suppose H is upper Hessenberg and that we compute the factorization $P H = L U$ via Gaussian elimination with partial pivoting. (See Algorithm 4.3.4.) Show that $H _ { 1 } = U ( P ^ { T } L )$ is upper Hessenberg and similar to H. (This is the basis of the modified LR algorithm.)

P7.5.5 Show that if $H = H _ { 0 }$ is given and we generate the matrices $H _ { k }$ via $H _ { k } - \mu _ { k } I = U _ { k } R _ { k } , H _ { k + 1 }$ $= R _ { k } U _ { k } + \mu _ { k } I ,$ then $( U _ { 1 } \cdot \cdot \cdot U _ { j } ) ( R _ { j } \cdot \cdot \cdot R _ { 1 } ) = ( H - \mu _ { 1 } I ) \cdot \cdot \cdot ( H - \mu _ { j } I )$ .

# Notes and References for 7.5

Historically important papers associated with the QR iteration include:

H. Rutishauser (1958). “Solution of Eigenvalue Problems with the LR Transformation,” Nat. Bur. Stand. App. Math. Ser. 49, 47–81.

J.G.F. Francis (1961). “The QR Transformation: A Unitary Analogue to the LR Transformation, Parts I and II” Comput. J. 4, 265–72, 332–345.

V.N. Kublanovskaya (1961). “On Some Algorithms for the Solution of the Complete Eigenvalue Problem,” Vychisl. Mat. Mat. Fiz 1(4), 555–570.

R.S. Martin and J.H. Wilkinson (1968). “The Modified LR Algorithm for Complex Hessenberg Matrices,” Numer. Math. 12, 369–376.

R.S. Martin, G. Peters, and J.H. Wilkinson (1970). “The QR Algorithm for Real Hessenberg Matrices,” Numer. Math. 14, 219–231.

For a general insight, we recommend:

D.S. Watkins (1982). “Understanding the QR Algorithm,” SIAM Review 24, 427–440.

D.S. Watkins (1993). “Some Perspectives on the Eigenvalue Problem,” SIAM Review 35, 430–471.

D.S. Watkins (2008). “The QR Algorithm Revisited,” SIAM Review 50, 133–145.   
D.S. Watkins (2011). “Francis’s Algorithm,” Amer. Math. Monthly 118, 387–403.

Papers concerned with the convergence of the method, shifting, deflation, and related matters include:   
P.A. Businger (1971). “Numerically Stable Deflation of Hessenberg and Symmetric Tridiagonal Matrices, BIT 11, 262–270.   
D.S. Watkins and L. Elsner (1991). “Chasing Algorithms for the Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 12, 374–384.   
D.S. Watkins and L. Elsner (1991). “Convergence of Algorithms of Decomposition Type for the Eigenvalue Problem,” Lin. Alg. Applic. 143, 19–47.   
J. Erxiong (1992). “A Note on the Double-Shift QL Algorithm,” Lin. Alg. Applic. 171, 121–132.   
A.A. Dubrulle and G.H. Golub (1994). “A Multishift QR Iteration Without Computation of the Shifts,” Numer. Algorithms 7, 173–181.   
D.S. Watkins (1996). “Forward Stability and Transmission of Shifts in the QR Algorithm,” SIAM J. Matrix Anal. Applic. 16, 469–487.   
D.S. Watkins (1996). “The Transmission of Shifts and Shift Blurring in the QR algorithm,” Lin. Alg. Applic. 241–3, 877–896.   
D.S. Watkins (1998). “Bulge Exchanges in Algorithms of QR Type,” SIAM J. Matrix Anal. Applic. 19, 1074–1096.   
R. Vandebril (2011). “Chasing Bulges or Rotations? A Metamorphosis of the QR-Algorithm” SIAM. J. Matrix Anal. Applic. 32, 217–247.   
Aspects of the balancing problem are discussed in:   
E.E. Osborne (1960). “On Preconditioning of Matrices,” J. ACM 7, 338–345.   
B.N. Parlett and C. Reinsch (1969). “Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors,” Numer. Math. 13, 292–304.   
D.S. Watkins (2006). “A Case Where Balancing is Harmful,” ETNA 23, 1–4.   
Versions of the algorithm that are suitable for companion matrices are discussed in:   
D.A. Bini, F. Daddi, and L. Gemignani (2004). “On the Shifted QR iteration Applied to Companion Matrices,” ETNA 18, 137–152.   
M. Van Barel, R. Vandebril, P. Van Dooren, and K. Frederix (2010). “Implicit Double Shift QR-Algorithm for Companion Matrices,” Numer. Math. 116, 177–212.   
Papers that are concerned with the high-performance implementation of the QR iteration include:   
Z. Bai and J.W. Demmel (1989). “On a Block Implementation of Hessenberg Multishift QR Iteration,” Int. J. High Speed Comput. 1, 97–112.   
R.A. Van De Geijn (1993). “Deferred Shifting Schemes for Parallel QR Methods,” SIAM J. Matrix Anal. Applic. 14, 180–194.   
D.S. Watkins (1994). “Shifting Strategies for the Parallel QR Algorithm,” SIAM J. Sci. Comput. 15, 953–958.   
G. Henry and R. van de Geijn (1996). “Parallelizing the QR Algorithm for the Unsymmetric Algebraic Eigenvalue Problem: Myths and Reality,” SIAM J. Sci. Comput. 17, 870–883.   
Z. Bai, J. Demmel, J. Dongarra, A. Petitet, H. Robinson, and K. Stanley (1997). “The Spectral Decomposition of Nonsymmetric Matrices on Distributed Memory Parallel Computers,” SIAM J. Sci. Comput. 18, 1446–1461.   
G. Henry, D.S. Watkins, and J. Dongarra (2002). “A Parallel Implementation of the Nonsymmetric QR Algorithm for Distributed Memory Architectures,” SIAM J. Sci. Comput. 24, 284–311.   
K. Braman, R. Byers, and R. Mathias (2002). “The Multishift QR Algorithm. Part I: Maintaining Well-Focused Shifts and Level 3 Performance,” SIAM J. Matrix Anal. Applic. 23, 929–947.   
K. Braman, R. Byers, and R. Mathias (2002). “The Multishift QR Algorithm. Part II: Aggressive Early Deflation,” SIAM J. Matrix Anal. Applic. 23, 948–973.   
M.R. Fahey (2003). “Algorithm 826: A Parallel Eigenvalue Routine for Complex Hessenberg Matrices,” ACM Trans. Math. Softw. 29, 326–336.   
D. Kressner (2005). “On the Use of Larger Bulges in the QR Algorithm,” ETNA 20, 50–63.   
D. Kressner (2008). “The Effect of Aggressive Early Deflation on the Convergence of the QR Algorithm,” SIAM J. Matrix Anal. Applic. 30, 805–821.
