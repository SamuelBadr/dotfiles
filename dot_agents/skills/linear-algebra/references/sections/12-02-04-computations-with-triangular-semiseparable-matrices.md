# 12.2.4 Computations with Triangular Semiseparable Matrices

Lower and upper triangular matrices that are also semiseparable can be written as follows:

$$
L \text {   lower   semiseparable   } \Rightarrow L = \mathbf {S} (u, v, t, u. * v, 0, 0, 0) = \operatorname{tril} (u v ^ {T}) \cdot * B (t) ^ {- T},
$$

$$
U \text {   upper   semiseparable   } \Rightarrow U = \mathbf {S} (0, 0, 0, p. * q, p, q, r) = \operatorname{triu} (p q ^ {T}) * B (r) ^ {- 1}.
$$

Operations with matrices that have this structure can be organized very efficiently. Consider the matrix-vector product

$$
y = \left(\operatorname{triu} \left(p q ^ {T}\right). * B (r) ^ {- 1}\right) x \tag {12.2.9}
$$

where $x , y , p , q \in \mathbb { R } ^ { n }$ and $r \in \mathbb { R } ^ { n - 1 }$ . This calculation has the form

$$
\left[ \begin{array}{c c c c} p _ {1} q _ {1} & p _ {1} r _ {1} q _ {2} & p _ {1} r _ {1} r _ {2} q _ {3} & p _ {1} r _ {1} r _ {2} r _ {3} q _ {4} \\ 0 & p _ {2} q _ {2} & p _ {2} r _ {2} q _ {3} & p _ {2} r _ {2} r _ {3} q _ {4} \\ 0 & 0 & p _ {3} q _ {3} & p _ {3} r _ {3} q _ {4} \\ 0 & 0 & 0 & p _ {4} q _ {4} \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \\ x _ {3} \\ x _ {4} \end{array} \right] = \left[ \begin{array}{c} y _ {1} \\ y _ {2} \\ y _ {3} \\ y _ {4} \end{array} \right].
$$

By grouping the $q \mathrm { ^ s }$ with the $x _ { \mathrm { ~ S ~ } } ^ { \prime }$ and extracting the $p \mathrm { { ^ { \circ } s } , }$ , we see that

$$
\mathrm{diag} (p _ {1}, p _ {2}, p _ {3}, p _ {4}) \left[ \begin{array}{c c c c} 1 & r _ {1} & r _ {1} r _ {2} & r _ {1} r _ {2} r _ {3} \\ 0 & 1 & r _ {2} & r _ {2} r _ {3} \\ 0 & 0 & 1 & r _ {3} \\ 0 & 0 & 0 & 1 \end{array} \right] \left[ \begin{array}{c} q _ {1} x _ {1} \\ q _ {2} x _ {2} \\ q _ {3} x _ {3} \\ q _ {4} x _ {4} \end{array} \right] = \left[ \begin{array}{c} y _ {1} \\ y _ {2} \\ y _ {3} \\ y _ {4} \end{array} \right].
$$

In other words, (12.2.9) is equivalent to

$$
y = p. * \left(B (r) ^ {- 1} (q. * x)\right).
$$

Given $x ,$ this is clearly an $O ( n )$ computation since bidiagonal system solving is $O ( n )$ . Indeed, y can be computed with just 4n flops.

Note that if y is given in (12.2.9) and p and q have nonzero components, then we can solve for x equally fast: $x = ( B ( r ) ( y . / p ) ) . / q$ .
