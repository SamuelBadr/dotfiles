# 7.7.5 Deflation

In describing the QZ iteration we may assume without loss of generality that A is an unreduced upper Hessenberg matrix and that B is a nonsingular upper triangular matrix. The first of these assertions is obvious, for if $a _ { k + 1 , k } = 0$ then

$$
A - \lambda B = \left[ \begin{array}{c c} A _ {1 1} - \lambda B _ {1 1} & A _ {1 2} - \lambda B _ {1 2} \\ 0 & A _ {2 2} - \lambda B _ {2 2} \end{array} \right] _ {n - k} ^ {k},
$$

and we may proceed to solve the two smaller problems $A _ { 1 1 } - \lambda B _ { 1 1 }$ and $A _ { 2 2 } - \lambda B _ { 2 2 }$ . On the other hand, if $b _ { k k } = 0$ for some k, then it is possible to introduce a zero in A’s $( n , n - 1 )$ position and thereby deflate. Illustrating by example, suppose $n = 5$ and $k = 3 \colon$

$$
A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right], \qquad B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

The zero on $B ^ { \prime } \mathrm { s }$ diagonal can be “pushed down” to the (5,5) position as follows using Givens rotations:

$$
A \gets Q _ {3 4} ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \gets Q _ {3 4} ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
A \leftarrow A Z _ {2 3} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \leftarrow B Z _ {2 3} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
A \gets Q _ {4 5} ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & \times & \times & \times \end{array} \right], \quad B \gets Q _ {4 5} ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 \end{array} \right],
$$

$$
A \leftarrow A Z _ {3 4} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \leftarrow B Z _ {3 4} ^ {T} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 \end{array} \right],
$$

$$
A \leftarrow A Z _ {4 5} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right], \quad B \leftarrow B Z _ {4 5} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 \end{array} \right].
$$

This zero-chasing technique is perfectly general and can be used to zero $^ { a _ { n , n - 1 } }$ regardless of where the zero appears along $B ^ { \prime } \mathrm { s }$ diagonal.
