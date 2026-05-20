# 12.2.5 The LU Factorization of a Semiseparable Matrix

Suppose $A = \mathbf { S } ( u , v , t , u . * v , p , q , r )$ is an $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ semiseparable matrix that has an LU factorization. It turns out that both L and $U$ are semiseparable and their respective representations can be computed with $O ( n )$ work:

for $k = n { - } 1 { : } - 1 { : } 1$

Using A’s representation, determine $\tau _ { k }$ so that if $\tilde { A } = M _ { k } A$ , where

$$
M _ {k} = \operatorname{diag} (I _ {k - 1}, \tilde {M} _ {k}, I _ {n - k - 1}), \qquad \tilde {M} _ {k} = \left[ \begin{array}{c c} 1 & 0 \\ - \tau_ {k} & 1 \end{array} \right],
$$

then $\tilde { A } ( k + 1 , 1 { : } k ) \ \mathrm { i s ~ z e r o }$ (12.2.10)

Compute the update $A = M _ { k } A$ by updating A’s representation

end

$$
U = A
$$

Note that if $M = M _ { 1 } \cdot \cdot \cdot M _ { n - 1 }$ , then $M A = U$ and $M = B ( \tau )$ with $\boldsymbol { \tau } = [ \tau _ { 1 } , \dots , \tau _ { n - 1 } ] ^ { T }$ . It follows that if $L = M ^ { - 1 }$ , then L is semiseparable from (12.2.4) and $A = L U$ . The challenge is to show that the updates $A = M _ { k } A$ preserve semiseparability.

To see what is involved, suppose $n = 6$ and that we have computed $M _ { 5 }$ and $M _ { 4 }$ so that

$$
M _ {4} M _ {5} A = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \lambda & \lambda & \lambda & \mu & \mu & \mu \\ \lambda & \lambda & \lambda & \mu & \mu & \mu \\ \hline 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right] = \mathbf {S} (u, v, t, u. * v, p, q, r)
$$

is semiseparable. Note that the λ-block and the µ-block are given by

$$
\left[ \begin{array}{c c c} \lambda & \lambda & \lambda \\ \lambda & \lambda & \lambda \end{array} \right] = \left[ \begin{array}{c c c} u _ {3} t _ {2} t _ {1} v _ {1} & u _ {3} t _ {2} v _ {2} & u _ {3} v _ {3} \\ u _ {4} t _ {3} t _ {2} t _ {1} v _ {1} & u _ {4} t _ {3} t _ {2} v _ {2} & u _ {4} t _ {3} v _ {3} \end{array} \right],
$$

$$
\left[ \begin{array}{c c c} \mu & \mu & \mu \\ \mu & \mu & \mu \end{array} \right] = \left[ \begin{array}{c c c} p _ {3} r _ {3} q _ {4} & p _ {3} r _ {3} r _ {4} q _ {5} & p _ {3} r _ {3} r _ {4} r _ {5} q _ {6} \\ p _ {4} q _ {4} & p _ {4} r _ {4} q _ {5} & p _ {4} r _ {4} r _ {5} q _ {6} \end{array} \right].
$$

Thus, if

$$
\tilde {M} _ {3} = \left[ \begin{array}{c c} 1 & 0 \\ - \tau_ {3} & 1 \end{array} \right],
$$

then

$$
\tilde {M} _ {3} \left[ \begin{array}{c c c} \lambda & \lambda & \lambda \\ \lambda & \lambda & \lambda \end{array} \right] = \left[ \begin{array}{c c c} u _ {3} t _ {2} t _ {1} v _ {1} & u _ {3} t _ {2} v _ {2} & u _ {3} v _ {3} \\ (u _ {4} t _ {3} - \tau_ {3} u _ {3}) t _ {2} t _ {1} v _ {1} & (u _ {4} t _ {3} - \tau_ {3} u _ {3}) t _ {2} v _ {2} & (u _ {4} t _ {3} - \tau_ {3} u _ {3}) v _ {3} \end{array} \right],
$$

$$
\tilde {M} _ {3} \left[ \begin{array}{c c c} \mu & \mu & \mu \\ \mu & \mu & \mu \end{array} \right] = \left[ \begin{array}{c c c} p _ {3} r _ {3} q _ {4} & p _ {3} r _ {3} r _ {4} q _ {5} & p _ {3} r _ {3} r _ {4} r _ {5} q _ {6} \\ (p _ {4} - \tau_ {3} p _ {3} r _ {3}) q _ {4} & (p _ {4} - \tau_ {3} p _ {3} r _ {3}) r _ {4} q _ {5} & (p _ {4} - \tau_ {3} p _ {3} r _ {3}) r _ {4} r _ {5} q _ {6} \end{array} \right].
$$

If $u _ { 3 } \neq 0 , \tau _ { 3 } = u _ { 4 } t _ { 3 } / u _ { 3 }$ , and we perform the updates

$$
u _ {4} = 0, \qquad p _ {4} = p _ {4} - \tau_ {3} p _ {3} r _ {3},
$$

then

$$
M _ {3} M _ {4} M _ {5} A = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \lambda & \lambda & \lambda & \mu & \mu & \mu \\ 0 & 0 & 0 & \tilde {\mu} & \tilde {\mu} & \tilde {\mu} \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right] = \mathbf {S} (u, v, t, u. * v, p, q, r)
$$

is still semiseparable. (The tildes designate updated entries.) Picking up the pattern from this example, we obtain the following $O ( n )$ method for computing the LU factorization of a semiseparable matrix.

Algorithm 12.2.1 Assume that $u , v , p , q \in \mathbb { R } ^ { n }$ with u . $\ast v = p . \ast q$ and that $t , r \in \mathbb { R } ^ { n - 1 }$ . If $\boldsymbol { A } = \mathbf { S } ( u , t , v , u$ . ∗ $v , p , r , q )$ has an LU factorization, then the following algorithm computes $\tilde { p } \in  { \mathbb { R } } ^ { n }$ and $\tau \in \mathbb { R } ^ { n - 1 }$ so that if $L = B ( \tau ) ^ { - T }$ and $U = \mathrm { t r i u } ( \tilde { p } q ^ { T } )$ . ∗ $\boldsymbol { B } ( \boldsymbol { r } ) ^ { - 1 }$ , then $A = L U$ .

for $k = n { - } 1 { : } - 1 { : } 1$

$$
\tau_ {k} = t _ {k} u _ {k + 1} / u _ {k}
$$

$$
\tilde {p} _ {k + 1} = p _ {k + 1} - p _ {k} \tau_ {k} r _ {k}
$$

end

$$
\tilde {p} _ {1} = p _ {1}
$$

This algorithm requires about 5n flops. Given our remarks in the previous section about triangular semiseparable matrices, we see that a semiseparable system $A x = b$ can be solved with $O ( n )$ work: A = LU , Ly = b, U x = y. Note that the vectors $\tau$ and $\tilde { p }$ in algorithm 12.2.1 are given by

$$
\tau = (u (2: n). * t). / u (1: n - 1)
$$

and

$$
\tilde {p} = \left[ \begin{array}{c} p _ {1} \\ p (2: n) - p (1: n - 1)  . *   \tau  . *   r \end{array} \right].
$$

Pivoting can be incorporated in Algorithm 12.2.1 to ensure that $| \tau _ { k } | \le 1$ for $k = n - 1 \colon - 1 \colon 1$ . At the beginning of step k, if $\vert u _ { k } \vert < \vert u _ { k + 1 } \vert$ , then rows k and $k +$ 1 are interchanged. The swapping is orchestrated by updating the quasiseparable respresentation of the current A. The end result is an $O ( n )$ reduction of the form $M _ { 1 } \cdot \cdot \cdot M _ { n - 1 } A = U$ where $U$ is upper triangular and quasiseparable and $\begin{array} { r l } { M _ { k } } & { { } = } \end{array}$ diag $( I _ { k - 1 } , \tilde { M } _ { k } \tilde { P } _ { k } , I _ { n - k - 1 } )$ with

$$
\tilde {P} _ {k} = \left[ \begin{array}{l l} 1 & 0 \\ 0 & 1 \end{array} \right] \text {or} \left[ \begin{array}{l l} 0 & 1 \\ 1 & 0 \end{array} \right].
$$

See Vandebril, Van Barel, and Mastronardi (2008, pp. 165–170) for further details and also how to perform the same tasks when A is quasiseparable.
