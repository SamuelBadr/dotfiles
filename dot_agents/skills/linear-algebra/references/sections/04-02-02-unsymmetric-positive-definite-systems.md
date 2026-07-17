# 4.2.2 Unsymmetric Positive Definite Systems

The positive definiteness of a general matrix A is inherited from its symmetric part:

$$
T = \frac {A + A ^ {T}}{2}.
$$

Note that for any square matrix we have $A = T + S$ where

$$
S = \frac {A - A ^ {T}}{2}
$$

is the skew-symmetric part of A. Recall that a matrix S is skew symmetric if $S ^ { T } = - S$ . If S is skew-symmetric, then $x ^ { T } S x = 0$ for all $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $s _ { i i } = 0 , i = 1 { : } n$ . It follows that A is positive definite if and only if its symmetric part is positive definite.

The derivation and analysis of methods for positive definite systems require an understanding about how the symmetric and skew-symmetric parts interact during the LU process.

Theorem 4.2.5. Suppose

$$
A = \left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right] + \left[ \begin{array}{c c} 0 & - w ^ {T} \\ w & C \end{array} \right]
$$

is positive definite and that $B \in \mathbb { R } ^ { ( n - 1 ) \times ( n - 1 ) }$ is symmetric and $C \in \mathbb { R } ^ { ( n - 1 ) \times ( n - 1 ) }$ is skew-symmetric. Then it follows that

$$
A = \left[ \begin{array}{c c} 1 & 0 \\ (v + w) / \alpha & I \end{array} \right] \left[ \begin{array}{c c} \alpha & (v - w) ^ {T} \\ 0 & B _ {1} + C _ {1} \end{array} \right] \tag {4.2.1}
$$

where

$$
B _ {1} = B - \frac {1}{\alpha} \left(v v ^ {T} - w w ^ {T}\right) \tag {4.2.2}
$$

is symmetric positive definite and

$$
C _ {1} = C - \frac {1}{\alpha} \left(w v ^ {T} - v w ^ {T}\right) \tag {4.2.3}
$$

is skew-symmetric.

Proof. Since $\alpha \neq 0$ it follows that (4.2.1) holds. It is obvious from their definitions that $B _ { 1 }$ is symmetric and that $C _ { 1 }$ is skew-symmetric. Thus, all we have to show is that $B _ { 1 }$ is positive definite i.e.,

$$
0 <   z ^ {T} B _ {1} z = z ^ {T} B z - \frac {1}{\alpha} (v ^ {T} z) ^ {2} + \frac {1}{\alpha} (w ^ {T} z) ^ {2} \tag {4.2.4}
$$

for all nonzero $z \in \mathbb { R } ^ { n - 1 }$ . For any $\mu \in \mathbb { R }$ and $0 \neq z \in \mathbb { R } ^ { n - 1 }$ we have

$$
\begin{array}{l} 0 <   \left[ \begin{array}{c} \mu \\ z \end{array} \right] ^ {T} A \left[ \begin{array}{c} \mu \\ z \end{array} \right] = \left[ \begin{array}{c} \mu \\ z \end{array} \right] ^ {T} \left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right] \left[ \begin{array}{c} \mu \\ z \end{array} \right] \\ = \mu^ {2} \alpha + 2 \mu v ^ {T} z + z ^ {T} B z. \\ \end{array}
$$

If $\mu = - ( v ^ { T } z ) / \alpha$ , then

$$
0 <   z ^ {T} B z - \frac {1}{\alpha} (v ^ {T} z) ^ {2},
$$

which establishes the inequality (4.2.4).

From (4.2.1) we see that if $B _ { 1 } + C _ { 1 } = L _ { 1 } U _ { 1 }$ is the LU factorization, then $A = L U$ where

$$
L   =   \left[ \begin{array}{c c} 1 & 0 \\ (v + w) / \alpha & L _ {1} \end{array} \right] \left[ \begin{array}{c c} \alpha & (v - w) ^ {T} \\ 0 & U _ {1} \end{array} \right].
$$

Thus, the theorem shows that triangular factors in $A = L U$ are nicely bounded if S is not too big compared to $T ^ { - 1 }$ . Here is a result that makes this precise:

Theorem 4.2.6. Let $A \in \mathbb { R } ^ { n \times n }$ be positive definite and set $T = ( A + A ^ { T } ) / 2$ and $S = ( A - A ^ { T } ) / 2$ . If A = LU is the LU factorization, then

$$
\left\| | L | | U | \right\| _ {F} \leq n \left(\left\| T \right\| _ {2} + \left\| S T ^ {- 1} S \right\| _ {2}\right). \tag {4.2.5}
$$

Proof. See Golub and Van Loan (1979).

The theorem suggests when it is safe not to pivot. Assume that the computed factors $\hat { L }$ and $\hat { U }$ satisfy

$$
\| | \hat {L} | | \hat {U} | \| _ {F} \leq c \| | L | | U | \| _ {F}, \tag {4.2.6}
$$

where c is a constant of modest size. It follows from (4.2.1) and the analysis in §3.3 that if these factors are used to compute a solution to $A x = b ,$ then the computed solution ˆx satisfies $( { \boldsymbol { A } } + { \boldsymbol { E } } ) { \hat { \boldsymbol { x } } } = { \boldsymbol { b } }$ with

$$
\left\| E \right\| _ {F} \leq \mathbf {u} \left(2 n \| A \| _ {F} + 4 c n ^ {2} \left(\| T \| _ {2} + \| S T ^ {- 1} S \| _ {2}\right)\right) + O \left(\mathbf {u} ^ {2}\right). \tag {4.2.7}
$$

It is easy to show that $\parallel T \parallel _ { 2 } \leq \parallel A \parallel _ { 2 }$ , and so it follows that if

$$
\Omega = \frac {\| S T ^ {- 1} S \| _ {2}}{\| A \| _ {2}} \tag {4.2.8}
$$

is not too large, then it is safe not to pivot. In other words, the norm of the skewsymmetric part S has to be modest relative to the condition of the symmetric part T . Sometimes it is possible to estimate Ω in an application. This is trivially the case when A is symmetric for then $\Omega = 0$ .
