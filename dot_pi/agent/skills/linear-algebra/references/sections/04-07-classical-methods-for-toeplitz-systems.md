# 4.7 Classical Methods for Toeplitz Systems

Matrices whose entries are constant along each diagonal arise in many applications and are called Toeplitz matrices. Formally, $T \in \mathbb { R } ^ { n \times n }$ is Toeplitz if there exist scalars $r _ { - n + 1 } , \ldots , r _ { 0 } , \ldots , r _ { n - 1 }$ such that $a _ { i j } = r _ { j - i }$ for all i and j. Thus,

$$
T = \left[ \begin{array}{l l l l} r _ {0} & r _ {1} & r _ {2} & r _ {3} \\ r _ {- 1} & r _ {0} & r _ {1} & r _ {2} \\ r _ {- 2} & r _ {- 1} & r _ {0} & r _ {1} \\ r _ {- 3} & r _ {- 2} & r _ {- 1} & r _ {0} \end{array} \right] = \left[ \begin{array}{l l l l} 3 & 1 & 7 & 6 \\ 4 & 3 & 1 & 7 \\ 0 & 4 & 3 & 1 \\ 9 & 0 & 4 & 3 \end{array} \right]
$$

is Toeplitz. In this section we show that Toeplitz systems can be solved in $O ( n ^ { 2 } )$ flops The discussion focuses on the important case when T is also symmetric and positive definite, but we also include a few comments about general Toeplitz systems. An alternative approach to Toeplitz system solving based on displacement rank is given in §12.1.

# 4.7.1 Persymmetry

The key fact that makes it possible to solve a Toeplitz system $T x = b$ so fast has to do with the structure of $T ^ { - 1 }$ . Toeplitz matrices belong to the larger class of persymmetric matrices. We say that $B \in \mathbb { R } ^ { n \times n }$ is persymmetric if

$$
\mathcal {E} _ {n} B \mathcal {E} _ {n} = B ^ {T}
$$

where ${ \mathcal { E } } _ { n }$ is the n-by-n exchange matrix defined in §1.2.11, e.g.,

$$
\mathcal {E} _ {4} = \left[ \begin{array}{c c c c} 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \end{array} \right].
$$

If B is persymmetric, then ${ \mathcal { E } } _ { n } B$ is symmetric. This means that B is symmetric about its antidiagonal. Note that the inverse of a persymmetric matrix is also persymmetric:

$$
\mathcal {E} _ {n} B ^ {- 1} \mathcal {E} _ {n} = (\mathcal {E} _ {n} B \mathcal {E} _ {n}) ^ {- 1} = (B ^ {T}) ^ {- 1} = (B ^ {- 1}) ^ {T}.
$$

Thus, the inverse of a nonsingular Toeplitz matrix is persymmetric.

# 4.7.2 Three Problems

Assume that we have scalars $r _ { 1 } , \ldots , r _ { n }$ such that for $k = 1 { : } n$ the matrices

$$
T _ {k} = \left[ \begin{array}{c c c c c} 1 & r _ {1} & \dots & r _ {k - 2} & r _ {k - 1} \\ r _ {1} & 1 & \ddots & & r _ {k - 2} \\ \vdots & \ddots & \ddots & \ddots & \vdots \\ r _ {k - 2} & & \ddots & \ddots & r _ {1} \\ r _ {k - 1} & r _ {k - 2} & \dots & r _ {1} & 1 \end{array} \right]
$$

are positive definite. (There is no loss of generality in normalizing the diagonal.) We set out to describe three important algorithms:

• Durbin’s algorithm for the Yule-Walker problem $T _ { n } y = - [ r _ { 1 } , \ldots , r _ { n } ] ^ { T }$   
• Levinson’s algorithm for the general right-hand-side problem $T _ { n } x = b$   
• Trench’s algorithm for computing $B = T _ { n } ^ { - 1 }$

# 4.7.3 Solving the Yule-Walker Equations

We begin by presenting Durbin’s algorithm for the Yule-Walker equations which arise in conjunction with certain linear prediction problems. Suppose for some k that satisfies $1 \leq k \leq n - 1$ we have solved the kth order Yule-Walker system $T _ { k } y = - r =$ $- [ r _ { 1 } , \ldots , r _ { k } ] ^ { T }$ . We now show how the $( k + 1 ) \mathrm { s t }$ order Yule-Walker system

$$
\left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ r ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] \left[ \begin{array}{c} z \\ \alpha \end{array} \right] = - \left[ \begin{array}{c} r \\ r _ {k + 1} \end{array} \right]
$$

can be solved in O(k) flops. First observe that

$$
z = T _ {k} ^ {- 1} (- r - \alpha \mathcal {E} _ {k} r) = y - \alpha T _ {k} ^ {- 1} \mathcal {E} _ {k} r
$$

and

$$
\alpha = - r _ {k + 1} - r ^ {T} \mathcal {E} _ {k} z.
$$

Since $T _ { k } ^ { - 1 }$ is persymmetric, $T _ { k } ^ { - 1 } { \mathcal { E } } _ { k } = { \mathcal { E } } _ { k } T _ { k } ^ { - 1 }$ and thus

$$
z = y - \alpha \mathcal {E} _ {k} T _ {k} ^ {- 1} r = y + \alpha \mathcal {E} _ {k} y.
$$

By substituting this into the above expression for α we find

$$
\alpha = - r _ {k + 1} - r ^ {T} \mathcal {E} _ {k} (y + \alpha \mathcal {E} _ {k} y) = - (r _ {k + 1} + r ^ {T} \mathcal {E} _ {k} y) / (1 + r ^ {T} y).
$$

The denominator is positive because $T _ { k + 1 }$ is positive definite and because

$$
\left[ \begin{array}{c c} I & \mathcal {E} _ {k} y \\ 0 & 1 \end{array} \right] ^ {T} \left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ r ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] \left[ \begin{array}{c c} I & \mathcal {E} _ {k} y \\ 0 & 1 \end{array} \right] = \left[ \begin{array}{c c} T _ {k} & 0 \\ 0 & 1 + r ^ {T} y \end{array} \right].
$$

We have illustrated the kth step of an algorithm proposed by Durbin (1960). It proceeds by solving the Yule-Walker systems

$$
T _ {k} y ^ {(k)} = - r ^ {(k)} = - \left[ r _ {1}, \dots , r _ {k} \right] ^ {T}
$$

for k = 1:n as follows:

$$
y ^ {(1)} = - r _ {1}
$$

for $k = 1 { : } n - 1$

$$
\beta_ {k} = 1 + [ r ^ {(k)} ] ^ {T} y ^ {(k)}
$$

$$
\alpha_ {k} = - (r _ {k + 1} + r ^ {(k) ^ {T}} \mathcal {E} _ {k} y ^ {(k)}) / \beta_ {k} \tag {4.7.1}
$$

$$
z ^ {(k)} = y ^ {(k)} + \alpha_ {k} \mathcal {E} _ {k} y ^ {(k)}
$$

$$
y ^ {(k + 1)} = \left[ \begin{array}{c} z ^ {(k)} \\ \alpha_ {k} \end{array} \right]
$$

end

As it stands, this algorithm would require $3 n ^ { 2 }$ flops to generate $y = y ^ { ( n ) }$ . It is possible, however, to reduce the amount of work even further by exploiting some of the above expressions:

$$
\begin{array}{l} \beta_ {k} = 1 + \left[ r ^ {(k)} \right] ^ {T} y ^ {(k)} \\ = 1 + \left[ \begin{array}{c} r ^ {(k - 1)} \\ r _ {k} \end{array} \right] ^ {T} \left[ \begin{array}{c} y ^ {(k - 1)} + \alpha_ {k - 1} \mathcal {E} _ {k - 1} y ^ {(k - 1)} \\ \alpha_ {k - 1} \end{array} \right] \\ = \left(1 + \left[ r ^ {(k - 1)} \right] ^ {T} y ^ {(k - 1)}\right) + \alpha_ {k - 1} \left(\left[ r ^ {(k - 1)} \right] ^ {T} \mathcal {E} _ {k - 1} y ^ {(k - 1)} + r _ {k}\right) \\ = \beta_ {k - 1} + \alpha_ {k - 1} (- \beta_ {k - 1} \alpha_ {k - 1}) \\ = (1 - \alpha_ {k - 1} ^ {2}) \beta_ {k - 1}. \\ \end{array}
$$

Using this recursion we obtain the following algorithm:

Algorithm 4.7.1 (Durbin) Given real numbers $r _ { 0 } , r _ { 1 } , \ldots , r _ { n }$ with $r _ { 0 } = 1$ such that $T = ( r _ { | i - j | } ) \in \mathbb { R } ^ { n \times n }$ is positive definite, the following algorithm computes $\boldsymbol { y } \in \mathbb { R } ^ { n }$ such that $T \dot { y } = - [ r _ { 1 } , \ldots , r _ { n } ] ^ { T }$ .

$$
y (1) = - r (1); \beta = 1; \alpha = - r (1)
$$

for $k = 1 { : } n - 1$

$$
\beta = (1 - \alpha^ {2}) \beta
$$

$$
\alpha = - \left(r (k + 1) + r (k: - 1: 1) ^ {T} y (1: k)\right) / \beta
$$

$$
z (1: k) = y (1: k) + \alpha y (k: - 1: 1)
$$

$$
y (1: k + 1) = \left[ \begin{array}{c} z (1: k) \\ \alpha \end{array} \right]
$$

end

This algorithm requires $2 n ^ { 2 }$ flops. We have included an auxiliary vector z for clarity, but it can be avoided.

# 4.7.4 The General Right-Hand-Side Problem

With a little extra work, it is possible to solve a symmetric positive definite Toeplitz system that has an arbitrary right-hand side. Suppose that we have solved the system

$$
T _ {k} x = b = \left[ b _ {1}, \dots , b _ {k} \right] ^ {T} \tag {4.7.2}
$$

for some k satisfying $1 \leq k < n$ and that we now wish to solve

$$
\left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ r ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] \left[ \begin{array}{l} v \\ \mu \end{array} \right] = \left[ \begin{array}{c} b \\ b _ {k + 1} \end{array} \right]. \tag {4.7.3}
$$

Here, $\boldsymbol { r } = [ r _ { 1 } , \ldots , r _ { k } ] ^ { T }$ as above. Assume also that the solution to the order-k Yule-Walker system $T _ { k } y = - r$ is also available. From $T _ { k } v + \mu \mathcal { E } _ { k } r = b$ it follows that

$$
v = T _ {k} ^ {- 1} (b - \mu \mathcal {E} _ {k} r) = x - \mu T _ {k} ^ {- 1} \mathcal {E} _ {k} r = x + \mu \mathcal {E} _ {k} y
$$

and so

$$
\begin{array}{l} \mu = b _ {k + 1} - r ^ {T} \mathcal {E} _ {k} v \\ = b _ {k + 1} - r ^ {T} \mathcal {E} _ {k} x - \mu r ^ {T} y \\ = \left(b _ {k + 1} - r ^ {T} \mathcal {E} _ {k} x\right) / \left(1 + r ^ {T} y\right). \\ \end{array}
$$

Consequently, we can effect the transition from (4.7.2) to (4.7.3) in $O ( k )$ flops.

Overall, we can efficiently solve the system $T _ { n } x = b$ by solving the systems

$$
T _ {k} x ^ {(k)} = b ^ {(k)} = [ b _ {1}, \dots , b _ {k} ] ^ {T}
$$

and

$$
T _ {k} y ^ {(k)} = - r ^ {(k)} = - [ r _ {1}, \dots , r _ {k} ] ^ {T}
$$

“in parallel” for $k = 1 { : } n$ . This is the gist of the Levinson algorithm.

Algorithm 4.7.2 (Levinson) Given $b \in \mathbb { R } ^ { n }$ and real numbers $1 = r _ { 0 } , r _ { 1 } , . . . , r _ { n }$ such that $T = ( r _ { | i - j | } ) \in \mathbb { R } ^ { n \times n }$ is positive definite, the following algorithm computes $\boldsymbol { x } \in \mathbb { R } ^ { n }$ such that $T x = b$ .

$$
\begin{array}{l} y (1) = - r (1); x (1) = b (1); \beta = 1; \alpha = - r (1) \\ \beta = (1 - \alpha^ {2}) \beta \\ \mu = \left(b (k + 1) - r (1: k) ^ {T} x (k: - 1: 1)\right) / \beta \\ v (1: k) = x (1: k) + \mu \cdot y (k: - 1: 1) \\ x (1: k + 1) = \left[ \begin{array}{c} v (1: k) \\ \mu \end{array} \right] \\ z (1: k) = y (1: k) + \alpha \cdot y (k: - 1: 1) \\ y (1: k + 1) = \left[ \begin{array}{c} z (1: k) \\ \alpha \end{array} \right] \\ \end{array}
$$

$$
\begin{array}{l} \beta = (1 - \alpha^ {2}) \beta \\ \mu = \left(b (k + 1) - r (1: k) ^ {T} x (k: - 1: 1)\right) / \beta \\ v (1: k) = x (1: k) + \mu \cdot y (k: - 1: 1) \\ x (1: k + 1) = \left[ \begin{array}{c} v (1: k) \\ \mu \end{array} \right] \\ \end{array}
$$

$$
\begin{array}{l} \alpha = - \left(r (k + 1) + r (1: k) ^ {T} y (k: - 1: 1)\right) / \beta \\ z (1: k) = y (1: k) + \alpha \cdot y (k: - 1: 1) \\ y (1: k + 1) = \left[ \begin{array}{c} z (1: k) \\ \alpha \end{array} \right] \\ \end{array}
$$

end

This algorithm requires $4 n ^ { 2 }$ flops. The vectors z and v are for clarity and can be avoided in a detailed implementation.

# 4.7.5 Computing the Inverse

One of the most surprising properties of a symmetric positive definite Toeplitz matrix $T _ { n }$ is that its complete inverse can be calculated in $O ( n ^ { 2 } )$ flops. To derive the algorithm for doing this, partition $T _ { n } ^ { - 1 }$ as follows:

$$
T _ {n} ^ {- 1} = \left[ \begin{array}{c c} A & E r \\ r ^ {T} E & 1 \end{array} \right] ^ {- 1} = \left[ \begin{array}{c c} B & v \\ v ^ {T} & \gamma \end{array} \right] \tag {4.7.4}
$$

where $A = T _ { n - 1 } , E = \mathcal { E } _ { n - 1 }$ , and $r = [ r _ { 1 } , \ldots , r _ { n - 1 } ] ^ { T }$ . From the equation

$$
\left[ \begin{array}{c c} A & E r \\ r ^ {T} E & 1 \end{array} \right] \left[ \begin{array}{c} v \\ \gamma \end{array} \right] = \left[ \begin{array}{c} 0 \\ 1 \end{array} \right]
$$

it follows that $A v = - \gamma E r = - \gamma E ( r _ { 1 } , \ldots , r _ { n - 1 } ) ^ { T }$ and $\gamma = 1 - r ^ { T } E v$ . If y solves the order-(n−1) Yule-Walker system $A y = - r$ , then these expressions imply that

$$
\gamma = 1 / (1 + r ^ {T} y),
$$

$$
v = \gamma E y.
$$

Thus, the last row and column of $T _ { n } ^ { - 1 }$ are readily obtained.

It remains for us to develop working formulae for the entries of the submatrix B in (4.7.4). Since $A B + \mathcal { E } r v ^ { T } = I _ { n - 1 }$ , it follows that

$$
B = A ^ {- 1} - (A ^ {- 1} E r) v ^ {T} = A ^ {- 1} + \frac {v v ^ {T}}{\gamma}.
$$

Now since $A = T _ { n - 1 }$ is nonsingular and Toeplitz, its inverse is persymmetric. Thus,

$$
\begin{array}{l} b _ {i j} = (A ^ {- 1}) _ {i j} + \frac {v _ {i} v _ {j}}{\gamma} \\ = (A ^ {- 1}) _ {n - j, n - i} + \frac {v _ {i} v _ {j}}{\gamma} \tag {4.7.5} \\ = b _ {n - j, n - i} - \frac {v _ {n - j} v _ {n - i}}{\gamma} + \frac {v _ {i} v _ {j}}{\gamma} \\ = b _ {n - j, n - i} + \frac {1}{\gamma} \left(v _ {i} v _ {j} - v _ {n - j} v _ {n - i}\right). \\ \end{array}
$$

This indicates that although B is not persymmetric, we can readily compute an element $b _ { i j }$ from its reflection across the northeast-southwest axis. Coupling this with the fact that $A ^ { - 1 }$ is persymmetric enables us to determine B from its “edges” to its “interior.”

Because the order of operations is rather cumbersome to describe, we preview the formal specification of the algorithm pictorially. To this end, assume that we know the last column and row of $T _ { n } ^ { - 1 }$ :

$$
T _ {n} ^ {- 1} = \left[ \begin{array}{l l l l l l} u & u & u & u & u & k \\ u & u & u & u & u & k \\ u & u & u & u & u & k \\ u & u & u & u & u & k \\ u & u & u & u & u & k \\ k & k & k & k & k & k \end{array} \right].
$$

Here $^ { 6 } u ^ { \dag }$ and $^ { 6 } k ^ { 7 }$ denote the unknown and the known entries, respectively, and $n =$ 6. Alternately exploiting the persymmetry of $T _ { n } ^ { - 1 }$ and the recursion (4.7.5), we can compute B, the leading $( n - 1 ) – \mathrm { b y } – ( n - 1 )$ block of $T _ { n } ^ { - 1 }$ , as follows:

$$
\stackrel {\mathrm{persym}} {\longrightarrow} \left[ \begin{array}{l l l l l l} k & k & k & k & k & k \\ k & u & u & u & u & k \\ k & u & u & u & u & k \\ k & u & u & u & u & k \\ k & u & u & u & u & k \\ k & k & k & k & k & k \end{array} \right] \stackrel {(4. 7. 5)} {\longrightarrow} \left[ \begin{array}{l l l l l l} k & k & k & k & k & k \\ k & u & u & u & k & k \\ k & u & u & u & k & k \\ k & u & u & u & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \end{array} \right] \stackrel {\mathrm{persym}} {\longrightarrow} \left[ \begin{array}{l l l l l l} k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & u & u & k & k \\ k & k & u & u & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \end{array} \right]
$$

$$
\stackrel {(4. 7. 5)} {\longrightarrow} \left[ \begin{array}{c c c c c c} k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & u & k & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \end{array} \right] \stackrel {\text {persym}} {\longrightarrow} \left[ \begin{array}{c c c c c c} k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \\ k & k & k & k & k & k \end{array} \right].
$$

Of course, when computing a matrix that is both symmetric and persymmetric, such as $T _ { n } ^ { - 1 }$ , it is only necessary to compute the “upper wedge” of the matrix—e.g.,

$$
\begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ & \times & \times & \times & \times \\ & & \times & \times \end{array} \qquad (n = 6).
$$

With this last observation, we are ready to present the overall algorithm.

Algorithm 4.7.3 (Trench) Given real numbers $1 = r _ { 0 } , r _ { 1 } , . . . , r _ { n }$ such that $T =$ $( r _ { | i - j | } ) \in \mathbb { R } ^ { n \times n }$ is positive definite, the following algorithm computes $B = T _ { n } ^ { - 1 }$ . Only those $b _ { i j }$ for which $i \leq j$ and $i + j \le n + 1$ are computed.

Use Algorithm 4.7.1 to solve $T _ { n - 1 } y = - ( r _ { 1 } , \ldots , r _ { n - 1 } ) ^ { T }$

$$
\gamma = 1 / (1 + r (1: n - 1) ^ {T} y (1: n - 1))
$$

$$
v (1: n - 1) = \gamma y (n - 1: - 1: 1)
$$

$$
B (1, 1) = \gamma
$$

$$
B (1, 2: n) = v (n - 1: - 1: 1) ^ {T}
$$

for i = 2 : floor((n − 1)/2) + 1

for j = i:n − i + 1

$$
B (i, j) = B (i - 1, j - 1) + (v (n + 1 - j) v (n + 1 - i) - v (i - 1) v (j - 1)) / \gamma
$$

end

end

This algorithm requires $1 3 n ^ { 2 } / 4$ flops.

# 4.7.6 Stability Issues

Error analyses for the above algorithms have been performed by Cybenko (1978), and we briefly report on some of his findings.

The key quantities turn out to be the $\alpha _ { k }$ in (4.7.1). In exact arithmetic these scalars satisfy

$$
\left| \alpha_ {k} \right| <   1
$$

and can be used to bound $\| T ^ { - 1 } \| _ { 1 }$ :

$$
\max \left\{\frac {1}{\prod_ {j = 1} ^ {n - 1} (1 - \alpha_ {j} ^ {2})}, \frac {1}{\prod_ {j = 1} ^ {n - 1} (1 - \alpha_ {j})} \right\} \leq \| T _ {n} ^ {- 1} \| \leq \prod_ {j = 1} ^ {n - 1} \frac {1 + | \alpha_ {j} |}{1 - | \alpha_ {j} |}. \tag {4.7.6}
$$

Moreover, the solution to the Yule-Walker system $T _ { n } y = - r ( 1 { : } n )$ satisfies

$$
\| y \| _ {1} = \left(\prod_ {k = 1} ^ {n} (1 + \alpha_ {k})\right) - 1 \tag {4.7.7}
$$

provided all the $\alpha _ { k }$ are nonnegative.

Now if ˆx is the computed Durbin solution to the Yule-Walker equations, then the vector $\boldsymbol { r } _ { D } = T _ { n } \boldsymbol { \hat { x } } + \boldsymbol { r }$ can be bounded as follows

$$
\parallel r _ {D} \parallel \approx \mathbf {u} \prod_ {k = 1} ^ {n} (1 + | \hat {\alpha} _ {k} |),
$$

where $\hat { \alpha } _ { k }$ is the computed version of $\alpha _ { k }$ . By way of comparison, since each $| r _ { i } |$ is bounded by unity, it follows that $\parallel r _ { c } \parallel \approx \mathbf { u } \parallel y \parallel _ { 1 }$ where $r _ { C }$ is the residual associated with the computed solution obtained via the Cholesky factorization. Note that the two residuals are of comparable magnitude provided (4.7.7) holds. Experimental evidence suggests that this is the case even if some of the $\alpha _ { k }$ are negative. Similar comments apply to the numerical behavior of the Levinson algorithm.

For the Trench method, the computed inverse $\hat { B }$ of $T _ { n } ^ { - 1 }$ can be shown to satisfy

$$
\frac {\| T _ {n} ^ {- 1} - \hat {B} \| _ {1}}{\| T _ {n} ^ {- 1} \| _ {1}} \approx \mathbf {u} \prod_ {k = 1} ^ {n} \frac {1 + | \hat {\alpha} _ {k} |}{1 - | \hat {\alpha} _ {k} |}.
$$

In light of (4.7.7) we see that the right-hand side is an approximate upper bound for $\mathbf { u } \parallel T _ { n } ^ { - 1 } \parallel$ which is approximately the size of the relative error when $T _ { n } ^ { - 1 }$ is calculated using the Cholesky factorization.

# 4.7.7 A Toeplitz Eigenvalue Problem

Our discussion of the symmetric eigenvalue problem begins in Chapter 8. However, we are able to describe a solution procedure for an important Toeplitz eigenvalue problem that does not require the heavy machinery from that later chapter. Suppose

$$
T = \left[ \begin{array}{l l} 1 & r ^ {T} \\ r & B \end{array} \right]
$$

is symmetric, positive definite, and Toeplitz with $r \in \mathbb { R } ^ { n - 1 }$ . Cybenko and Van Loan (1986) show how to pair the Durbin algorithm with Newton’s method to compute $\lambda _ { \operatorname* { m i n } } ( T )$ assuming that

$$
\lambda_ {\min} (T) <   \lambda_ {\min} (B). \tag {4.7.8}
$$

This assumption is typically the case in practice. If

$$
\left[ \begin{array}{c c} 1 & r ^ {T} \\ r & B \end{array} \right] \left[ \begin{array}{c} \alpha \\ y \end{array} \right] = \lambda_ {\min} \left[ \begin{array}{c} \alpha \\ y \end{array} \right],
$$

then $y = - \alpha ( B - \lambda _ { \operatorname* { m i n } } I ) ^ { - 1 } r , \alpha \neq 0$ , and

$$
\alpha + r ^ {T} \left[ - \alpha (B - \lambda_ {\mathrm{min}} I) ^ {- 1} r \right] = \lambda_ {\mathrm{min}} \alpha .
$$

Thus, $\lambda _ { \mathrm { m i n } }$ is a zero of the rational function

$$
f (\lambda) = 1 - \lambda - r ^ {T} (B - \lambda I) ^ {- 1} r.
$$

Note that if $\lambda < \lambda _ { \operatorname* { m i n } } ( B )$ , then

$$
f ^ {\prime} (\lambda) = - 1 - \left\| (B - \lambda I) ^ {- 1} r \right\| _ {2} ^ {2} \leq - 1,
$$

$$
f ^ {\prime \prime} (\lambda) = - 2 r ^ {T} (B - \lambda I) ^ {- 3} r \leq 0.
$$

Using these facts it can be shown that if

$$
\lambda_ {\min} (T) \leq \lambda^ {(0)} <   \lambda_ {\min} (B), \tag {4.7.9}
$$

then the Newton iteration

$$
\lambda^ {(k + 1)} = \lambda^ {(k)} - \frac {f (\lambda^ {(k)})}{f ^ {\prime} (\lambda^ {(k)})} \tag {4.7.10}
$$

converges to $\lambda _ { \operatorname* { m i n } } ( T )$ monotonically from the right. The iteration has the form

$$
\lambda^ {(k + 1)} = \lambda^ {(k)} + \frac {1 + r ^ {T} w - \lambda^ {(k)}}{1 + w ^ {T} w},
$$

where w solves the “shifted” Yule-Walker system

$$
(B - \lambda^ {(k)} I) w = - r.
$$

Since $\lambda ^ { ( k ) } < \lambda _ { \operatorname* { m i n } } ( B )$ , this system is positive definite and the Durbin algorithm (Algorithm 4.7.1) can be applied to the normalized Toeplitz matrix $( B - \lambda ^ { ( k ) } I ) / ( 1 - \dot { \lambda } ^ { ( k ) } )$ .

The Durbin algorithm can also be used to determine a starting value $\lambda ^ { ( 0 ) }$ that satisfies (4.7.9). If that algorithm is applied to

$$
T _ {\lambda} = (T - \lambda I) / (1 - \lambda)
$$

then it runs to completion if $T _ { \lambda }$ is positive definite. In this case, the $\beta _ { k }$ defined in (4.7.1) are all positive. On the other hand, if $k \leq n - 1 , \beta _ { k } \leq 0$ and $\beta _ { 1 } , \ldots , \beta _ { k - 1 }$ are all positive, then it follows that $T _ { \lambda } ( 1 { : } k , 1 { : } k )$ is positive definite but that $T _ { \lambda } ( 1 { : } k + 1 , k + 1 )$ i s not. Let $m ( \lambda )$ be the index of the first nonpositive $\beta$ and observe that if $m ( \lambda ^ { ( 0 ) } ) = n { - } 1$ , then $B - \lambda ^ { ( 0 ) } I$ is positive definite and $T - \lambda ^ { ( 0 ) } I$ is not, thereby establishing (4.7.9). A bisection scheme can be formulated to compute $\lambda ^ { ( 0 ) }$ with this property:

$$
L = 0
$$

$$
R = 1 - \left| r _ {1} \right|
$$

$$
\mu = (L + R) / 2
$$

while $m ( \mu ) \neq n - 1$

$$
\text { if } m (\mu) <   n - 1
$$

$$
R = \mu
$$

else (4.7.11)

$$
L = \mu
$$

end

$$
\mu = (L + R) / 2
$$

end

$$
\lambda^ {(0)} = \mu
$$

At all times during the iteration we have $m ( L ) \leq n - 1 \leq m ( R )$ . The initial value for R follows from the inequality

$$
0 <   \lambda_ {\min} (T) <   \lambda_ {\min} (B) \leq \lambda_ {\min} \left(\left[ \begin{array}{c c} 1 & r _ {1} \\ r _ {1} & 1 \end{array} \right]\right) = 1 - | r _ {1} |.
$$

Note that the iterations in (4.7.10) and (4.7.11) involve at most $O ( n ^ { 2 } )$ flops per pass. A heuristic argument that O(log n) iterations are required is given by Cybenko and Van Loan (1986).

# 4.7.8 Unsymmetric Toeplitz System Solving

We close with some remarks about unsymmetric Toeplitz system-solving. Suppose we are given scalars $r _ { 1 } , \ldots , r _ { n - 1 } , p _ { 1 } , \ldots , p _ { n - 1 }$ , and $b _ { 1 } , \ldots , b _ { n }$ and that we want to solve a linear system $T x = b$ of the form

$$
\left[ \begin{array}{c c c c c} 1 & r _ {1} & r _ {2} & r _ {3} & r _ {4} \\ p _ {1} & 1 & r _ {1} & r _ {2} & r _ {3} \\ p _ {2} & p _ {1} & 1 & r _ {1} & r _ {2} \\ p _ {3} & p _ {2} & p _ {1} & 1 & r _ {1} \\ p _ {4} & p _ {3} & p _ {2} & p _ {1} & 1 \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \\ x _ {3} \\ x _ {4} \\ x _ {5} \end{array} \right] = \left[ \begin{array}{c} b _ {1} \\ b _ {2} \\ b _ {3} \\ b _ {4} \\ b _ {5} \end{array} \right] \qquad (n = 5).
$$

Assume that $T _ { k } = T ( 1 { : } k , 1 { : } k )$ is nonsingular for $k = 1 { : } n$ . It can shown that if we have the solutions to the k-by-k systems

$$
T _ {k} ^ {T} y = - r = - \left[ r _ {1} r _ {2} \dots r _ {k} \right] ^ {T},
$$

$$
T _ {k} w = - p = - \left[ p _ {1} p _ {2} \dots p _ {k} \right] ^ {T}, \tag {4.7.12}
$$

$$
T _ {k} x = b = \left[ b _ {1} b _ {2} \dots b _ {k} \right] ^ {T},
$$

then we can obtain solutions to

$$
\begin{array}{l} \left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ p ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] ^ {T} \left[ \begin{array}{c} z \\ \alpha \end{array} \right] = - \left[ \begin{array}{c} r \\ r _ {k + 1} \end{array} \right], \\ \left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ p ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] \left[ \begin{array}{l} u \\ \nu \end{array} \right] = - \left[ \begin{array}{c} p \\ p _ {k + 1} \end{array} \right], \tag {4.7.13} \\ \left[ \begin{array}{c c} T _ {k} & \mathcal {E} _ {k} r \\ p ^ {T} \mathcal {E} _ {k} & 1 \end{array} \right] \left[ \begin{array}{c} v \\ \mu \end{array} \right] = \left[ \begin{array}{c} b \\ b _ {k + 1} \end{array} \right] \\ \end{array}
$$

in $O ( k )$ flops. The update formula derivations are very similar to the Levinson algorithm derivations in §4.7.3. Thus, if the process is repeated for $k = 1 { : } n - 1$ , then we emerge with the solution to $T x = T _ { n } x = b$ . Care must be exercised if a $T _ { k }$ matrix is singular or ill-conditioned. One strategy involves a lookahead idea. In this framework, one might transition from the $T _ { k }$ problem directly to the $T _ { k + 2 }$ problem if it is deemed that the $T _ { k + 1 }$ problem is dangerously ill-conditioned. See Chan and Hansen (1992). An alternative approach based on displacement rank is given in §12.1.

# Problems

P4.7.1 For any $v \in \mathbb { R } ^ { n }$ define the vectors $v _ { + } = ( v + \mathcal { E } _ { n } v ) / 2$ and $v _ { - } = ( v - \mathcal { E } _ { n } v ) / 2$ . Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and persymmetric. Show that if $A x = b$ then $A x _ { + } = b _ { + }$ and $A x _ { - } = b _ { - }$ .

P4.7.2 Let $U \in \mathbb { R } ^ { n \times n }$ be the unit upper triangular matrix with the property that $U ( 1 { : } k - 1 , k ) =$ $\mathcal { E } _ { k - 1 } y ^ { ( k - 1 ) }$ where $y ^ { ( k ) }$ is defined by (4.7.1). Show that $U ^ { T } T _ { n } U = \operatorname { d i a g } ( 1 , \beta _ { 1 } , . . . , \beta _ { n - 1 } )$ .

P4.7.3 Suppose that $z \in \mathbb { R } ^ { n }$ and that $S \in \mathbb { R } ^ { n \times n }$ is orthogonal. Show that if $X = \left[ z , \ S z , \ . . . , S ^ { n - 1 } z \right]$ , then $X ^ { T } X$ is Toeplitz.

P4.7.4 Consider the $\mathrm { L D L ^ { T } }$ factorization of an n-by-n symmetric, tridiagonal, positive definite Toeplitz matrix. Show that $d _ { n }$ and $\ell _ { n , n - 1 }$ converge as $n \to \infty$ .

P4.7.5 Show that the product of two lower triangular Toeplitz matrices is Toeplitz.

P4.7.6 Give an algorithm for determining $\mu \in \mathbb { R }$ such that $T _ { n } + \mu \left( e _ { n } e _ { 1 } ^ { T } + e _ { 1 } e _ { n } ^ { T } \right)$ is singular. Assume $T _ { n } = ( r _ { | i - j | } )$ is positive definite, with $r _ { 0 } = 1$ .

P4.7.7 Suppose $T \in \mathbb { R } ^ { n \times n }$ is symmetric, positive definite, and Toeplitz with unit diagonal. What is the smallest perturbation of the the ith diagonal that makes T semidefinite?

P4.7.8 Rewrite Algorithm 4.7.2 so that it does not require the vectors z and v.

P4.7.9 Give an algorithm for computing $\kappa _ { \infty } ( T _ { k } )$ for $k = 1 { : } n$ .

P4.7.10 A p-by-p block matrix $A ~ = ~ ( A _ { i j } )$ with m-by-m blocks is block Toeplitz if there exist $A _ { - p + 1 } , \hdots , A _ { - 1 } , A _ { 0 } , A _ { 1 } , \hdots , A _ { p - 1 } \in \mathbb { R } ^ { m \times m }$ so that $\begin{array} { r } { A _ { i j } = A _ { i - j } , \mathrm { e . g . } } \end{array}$ ,

$$
A = \left[ \begin{array}{c c c c} A _ {0} & A _ {1} & A _ {2} & A _ {3} \\ A _ {- 1} & A _ {0} & A _ {1} & A _ {2} \\ A _ {- 2} & A _ {- 1} & A _ {0} & A _ {1} \\ A _ {- 3} & A _ {- 2} & A _ {- 1} & A _ {0} \end{array} \right].
$$

(a) Show that there is a permutation Π such that

$$
\Pi^ {T} A \Pi =: \left[ \begin{array}{c c c c} T _ {1 1} & T _ {1 2} & \dots & T _ {1 m} \\ T _ {2 1} & T _ {2 2} & & \vdots \\ \vdots & & \ddots & \vdots \\ T _ {m 1} & \dots & & T _ {m m} \end{array} \right]
$$

where each $T _ { i j }$ is p-by-p and Toeplitz. Each $T _ { i j }$ should be “made up” of $( i , j )$ entries selected from the $A _ { k }$ matrices. (b) What can you say about the $T _ { i j } \mathrm { ~ i f ~ } A _ { k } = A _ { - k } , \bar { k } = 1 : p - 1 ?$

P4.7.11 Show how to compute the solutions to the systems in (4.7.13) given that the solutions to the systems in (4.7.12) are available. Assume that all the matrices involved are nonsingular. Proceed to develop a fast unsymmetric Toeplitz solver for $T x = b$ assuming that $T \ ' _ { \mathrm { s } }$ leading principal submatrices are all nonsingular.

P4.7.12 Consider the order-k Yule-Walker system $T _ { k } y ^ { ( k ) } = - r ^ { ( k ) }$ that arises in (4.7.1). Show that if $y ^ { ( k ) } = [ y _ { k 1 } , \dots , y _ { k k } ] ^ { T }$ for k = 1:n − 1 and

$$
L = \left[ \begin{array}{c c c c c c} 1 & 0 & 0 & 0 & \dots & 0 \\ y _ {1 1} & 1 & 0 & 0 & \dots & 0 \\ y _ {2 2} & y _ {2 1} & 1 & 0 & \dots & 0 \\ \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\ y _ {n - 1, n - 1} & y _ {n - 1, n - 2} & y _ {n - 1, n - 3} & \dots & y _ {n - 1, 1} & 1 \end{array} \right],
$$

then $L ^ { T } T _ { n } L \ = \ \mathrm { d i a g } ( 1 , \beta _ { 1 } , . . . , \beta _ { n - 1 } )$ where $\beta _ { k } = 1 + r ^ { ( k ) ^ { T } } y ^ { ( k ) }$ . Thus, the Durbin algorithm can be thought of as a fast method for computing and $\mathrm { L D L } ^ { T }$ factorization of $T _ { n } ^ { - 1 }$ .

P4.7.13 Show how the Trench algorithm can be used to obtain an initial bracketing interval for the bisection scheme (4.7.11).

# Notes and References for §4.7

The original references for the three algorithms described in this section are as follows:

J. Durbin (1960). “The Fitting of Time Series Models,”Rev. Inst. Int. Stat. 28, 233–243.   
N. Levinson (1947). “The Weiner RMS Error Criterion in Filter Design and Prediction,” J. Math. Phys. 25, 261–278.   
W.F. Trench (1964). “An Algorithm for the Inversion of Finite Toeplitz Matrices,” J. SIAM 12, 515–522.   
As is true with the “fast algorithms” area in general, unstable Toeplitz techniques abound and caution must be exercised, see:   
G. Cybenko (1978). “Error Analysis of Some Signal Processing Algorithms,” PhD Thesis, Princeton University.   
G. Cybenko (1980). “The Numerical Stability of the Levinson-Durbin Algorithm for Toeplitz Systems of Equations,” SIAM J. Sci. Stat. Comput. 1, 303-319.   
J.R. Bunch (1985). “Stability of Methods for Solving Toeplitz Systems of Equations,” SIAM J. Sci. Stat. Comput. 6, 349–364.   
E. Linzer (1992). “On the Stability of Solution Methods for Band Toeplitz Systems,” Lin. Alg. Applic. 170, 1–32.   
J.M. Varah (1994). “Backward Error Estimates for Toeplitz Systems,” SIAM J. Matrix Anal. Applic. 15, 408–417.   
A.W. Bojanczyk, R.P. Brent, F.R. de Hoog, and D.R. Sweet (1995). “On the Stability of the Bareiss and Related Toeplitz Factorization Algorithms,” SIAM J. Matrix Anal. Applic. 16, 40–57.   
M.T. Chu, R.E. Funderlic, and R.J. Plemmons (2003). “Structured Low Rank Approximation,” Lin. Alg. Applic. 366, 157–172.   
A. Bottcher and S. M. Grudsky (2004). “Structured Condition Numbers of Large Toeplitz Matrices are Rarely Better than Usual Condition Numbers,” Num. Lin. Alg. 12, 95–102.   
J.-G. Sun (2005). “A Note on Backwards Errors for Structured Linear Systems,” Numer. Lin. Alg. Applic. 12, 585–603.   
P. Favati, G. Lotti, and O. Menchi (2010). “Stability of the Levinson Algorithm for Toeplitz-Like Systems,” SIAM J. Matrix Anal. Applic. 31, 2531–2552.   
Papers concerned with the lookahead idea include:   
T.F. Chan and P. Hansen (1992). “A Look-Ahead Levinson Algorithm for Indefinite Toeplitz Systems,” SIAM J. Matrix Anal. Applic. 13, 490–506.   
M. Gutknecht and M. Hochbruck (1995). “Lookahead Levinson and Schur Algorithms for Nonhermitian Toeplitz Systems,” Numer. Math. 70, 181–227.

M. Van Barel and A. Bultheel (1997). “A Lookahead Algorithm for the Solution of Block Toeplitz Systems,” Lin. Alg. Applic. 266, 291–335.

Various Toeplitz eigenvalue computations are presented in:

G. Cybenko and C. Van Loan (1986). “Computing the Minimum Eigenvalue of a Symmetric Positive Definite Toeplitz Matrix,” SIAM J. Sci. Stat. Comput. 7, 123–131.

W.F. Trench (1989). “Numerical Solution of the Eigenvalue Problem for Hermitian Toeplitz Matrices,” SIAM J. Matrix Anal. Appl. 10, 135–146.

H. Voss (1999). “Symmetric Schemes for Computing the Minimum Eigenvalue of a Symmetric Toeplitz Matrix,” Lin. Alg. Applic. 287, 359–371.

A. Melman (2004). “Computation of the Smallest Even and Odd Eigenvalues of a Symmetric Positive-Definite Toeplitz Matrix,” SIAM J. Matrix Anal. Applic. 25, 947–963.
