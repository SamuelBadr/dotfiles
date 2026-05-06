# 7.7.7 The Overall QZ Process

By applying a sequence of QZ steps to the Hessenberg-triangular pencil $A - \lambda B$ , it is possible to reduce A to quasi-triangular form. In doing this it is necessary to monitor A’s subdiagonal and B’s diagonal in order to bring about decoupling whenever possible. The complete process, due to Moler and Stewart (1973), is as follows:

Algorithm 7.7.3 Given $A \in \mathbb { R } ^ { n \times n }$ and $B \in \mathbb { R } ^ { n \times n }$ , the following algorithm computes orthogonal Q and Z such that $Q ^ { T } A Z = T$ is upper quasi-triangular and $Q ^ { T } B \bar { Z } = S$ is upper triangular. A is overwritten by T and B by S.

Using Algorithm 7.7.1, overwrite A with $Q ^ { T } A Z$ (upper Hessenberg) and B with $Q ^ { T } B Z$ (upper triangular).

$$
\text { until } q = n
$$

Set to zero subdiagonal entries that satisfy $| a _ { i , i - 1 } | \leq \epsilon ( | a _ { i - 1 , i - 1 } | + | a _ { i i } | )$ .

Find the largest nonnegative q and the smallest nonnegative p such that if

$$
A = \left[ \begin{array}{c c c} A _ {1 1} & A _ {1 2} & A _ {1 3} \\ 0 & A _ {2 2} & A _ {2 3} \\ 0 & 0 & A _ {3 3} \end{array} \right] \begin{array}{c} p \\ n - p - q \\ q \end{array}
$$

then $A _ { 3 3 }$ is upper quasi-triangular and $A _ { 2 2 }$ is upper Hessenberg and unreduced.

Partition B conformably:

$$
B = \left[ \begin{array}{c c c} B _ {1 1} & B _ {1 2} & B _ {1 3} \\ 0 & B _ {2 2} & B _ {2 3} \\ 0 & 0 & B _ {3 3} \end{array} \right] \begin{array}{c} p \\ n - p - q \\ q \end{array}
$$

if $q < n$

if $B _ { 2 2 }$ is singular

Zero $a _ { n - q , n - q - 1 }$

Apply Algorithm 7.7.2 to $A _ { 2 2 }$ and $B _ { 2 2 }$ and update:

$$
A = \mathrm{diag} (I _ {p}, Q, I _ {q}) ^ {T} A \cdot \mathrm{diag} (I _ {p}, Z, I _ {q})
$$

$$
B = \mathrm{diag} (I _ {p}, Q, I _ {q}) ^ {T} B \cdot \mathrm{diag} (I _ {p}, Z, I _ {q})
$$

end end end

This algorithm requires $3 0 n ^ { 3 }$ flops. If $Q$ is desired, an additional $1 6 n ^ { 3 }$ are necessary. If Z is required, an additional $2 0 \bar { n } ^ { 3 }$ are needed. These estimates of work are based on the experience that about two $\mathrm { Q Z }$ iterations per eigenvalue are necessary. Thus, the convergence properties of $\mathrm { Q Z }$ are the same as for QR. The speed of the QZ algorithm is not affected by rank deficiency in B.

The computed S and T can be shown to satisfy

$$
Q _ {0} ^ {T} (A + E) Z _ {0} = T, \quad Q _ {0} ^ {T} (B + F) Z _ {0} = S,
$$

where $Q _ { 0 }$ and $Z _ { 0 }$ are exactly orthogonal and $\parallel E \parallel _ { 2 } \approx \mathbf { u } \parallel A \parallel _ { 2 }$ and $\parallel F \parallel _ { 2 } \approx \mathbf { u } \parallel B \parallel _ { 2 }$ .
