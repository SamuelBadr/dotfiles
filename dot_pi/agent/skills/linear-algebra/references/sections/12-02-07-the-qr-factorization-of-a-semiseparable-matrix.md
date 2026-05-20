# 12.2.7 The QR Factorization of a Semiseparable Matrix

The matrix $Q$ in the QR factorization of a semiseparable matrix $A \in \mathbb { R } ^ { n \times n }$ has a very simple form. Indeed, it is a product of Givens rotations $Q ^ { T } = G _ { 1 } \cdot \cdot \cdot G _ { n - 1 }$ where the underlying cosine-sine pairs are precisely those that define Givens representation of $A _ { L }$ . To see this, consider how easy it is to compute the QR factorization of $A _ { L }$ :

$$
\left[ \begin{array}{c c c c} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & c _ {3} & s _ {3} \\ 0 & 0 & - s _ {3} & c _ {3} \end{array} \right] \left[ \begin{array}{c c c c} c _ {1} v _ {1} & 0 & 0 & 0 \\ c _ {2} s _ {1} v _ {1} & c _ {2} v _ {2} & 0 & 0 \\ c _ {3} s _ {2} s _ {1} v _ {1} & c _ {3} s _ {2} v _ {2} & c _ {3} v _ {3} & 0 \\ s _ {3} s _ {2} s _ {1} v _ {1} & s _ {3} s _ {2} v _ {2} & s _ {3} v _ {3} & v _ {4} \end{array} \right] = \left[ \begin{array}{c c c c} c _ {1} v _ {1} & 0 & 0 & 0 \\ c _ {2} s _ {1} v _ {1} & c _ {2} v _ {2} & 0 & 0 \\ s _ {2} s _ {1} v _ {1} & s _ {2} v _ {2} & v _ {3} & s _ {3} v _ {4} \\ 0 & 0 & 0 & c _ {3} v _ {4} \end{array} \right],
$$

$$
\left[ \begin{array}{c c c c} 1 & 0 & 0 & 0 \\ 0 & c _ {2} & s _ {2} & 0 \\ 0 & - s _ {2} & c _ {2} & 0 \\ 0 & 0 & 0 & 1 \end{array} \right] \left[ \begin{array}{c c c c} c _ {1} v _ {1} & 0 & 0 & 0 \\ c _ {2} s _ {1} v _ {1} & c _ {2} v _ {2} & 0 & 0 \\ s _ {2} s _ {1} v _ {1} & s _ {2} v _ {2} & v _ {3} & s _ {3} v _ {4} \\ 0 & 0 & 0 & c _ {3} v _ {4} \end{array} \right] = \left[ \begin{array}{c c c c} c _ {1} v _ {1} & 0 & 0 & 0 \\ s _ {1} v _ {1} & v _ {2} & s _ {2} v _ {3} & s _ {2} s _ {3} v _ {4} \\ 0 & 0 & c _ {2} v _ {3} & c _ {2} s _ {3} v _ {4} \\ 0 & 0 & 0 & c _ {3} v _ {4} \end{array} \right],
$$

$$
\left[ \begin{array}{c c c c} c _ {1} & s _ {1} & 0 & 0 \\ - s _ {1} & c _ {1} & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{array} \right] \left[ \begin{array}{c c c c} c _ {1} v _ {1} & 0 & 0 & 0 \\ s _ {1} v _ {1} & v _ {2} & s _ {2} v _ {3} & s _ {2} s _ {3} v _ {4} \\ 0 & 0 & c _ {2} v _ {3} & c _ {2} s _ {3} v _ {4} \\ 0 & 0 & 0 & c _ {3} v _ {4} \end{array} \right] = \left[ \begin{array}{c c c c} v _ {1} & s _ {1} v _ {2} & s _ {1} s _ {2} v _ {3} & s _ {1} s _ {2} s _ {3} v _ {4} \\ 0 & c _ {1} v _ {2} & c _ {1} s _ {2} v _ {3} & c _ {1} s _ {2} s _ {3} v _ {4} \\ 0 & 0 & c _ {2} v _ {3} & c _ {2} s _ {3} v _ {4} \\ 0 & 0 & 0 & c _ {3} v _ {4} \end{array} \right].
$$

In general, if $\mathsf { t r i l } ( A ) \ : = \ : B ( s ) ^ { - T }$ .∗ $\mathrm { t r i l } ( c v ^ { T } )$ is a Givens vector representation and

$$
Q ^ {T} = G _ {1} \dots G _ {n - 1} \tag {12.2.12}
$$

where

$$
G _ {k} = \operatorname{diag} (I _ {k - 1}, \tilde {G} _ {k}, I _ {n - k - 1}), \quad \tilde {G} _ {k} = \left[ \begin{array}{c c} c _ {k} & s _ {k} \\ - s _ {k} & c _ {k} \end{array} \right], \tag {12.2.13}
$$

for k = 1:n − 1, then

$$
Q ^ {T} \operatorname{tril} (A) = R _ {L} = \operatorname{triu} \left(\left(\mathcal {D} _ {n} c\right) v ^ {T}\right). * B (s) ^ {- 1}. \tag {12.2.14}
$$

(Recall that $\mathcal { D } _ { n }$ is the downshift permutation, see §1.3.x.) Since $Q ^ { T }$ is upper Hessenberg, it follows that

$$
Q ^ {T} \operatorname{triu} (\mathrm{A}, 1) = R _ {U}
$$

is also upper triangular. Thus,

$$
Q ^ {T} A = Q ^ {T} (A _ {L} + A _ {U}) = R _ {L} + R _ {U} = R
$$

is the QR factorization of A. Unfortunately, this is not a useful O(n) representation of R from the standpoint of solving $A x = b$ because the summation gets in the way when we try to solve $( R _ { L } + R _ { U } ) x = Q ^ { T } b$ .

Fortunately, there is a handier way to encode R. Assume for clarity that A has a generator representation

$$
A = \operatorname{tril} \left(u v ^ {T}\right) + \operatorname{triu} \left(p q ^ {T}\right), \tag {12.2.15}
$$

where $u , v , p , q \in \mathbb { R } ^ { n }$ and $u . * v = p . * q$ . We show that R is the upper triangular portion of a rank-2 matrix, i.e.,

$$
R = \operatorname{triu} \left(f g ^ {T} + h q ^ {T}\right), \quad f, g, h \in \mathbb {R} ^ {n}. \tag {12.2.16}
$$

This means that any submatrix extracted from the upper triangular part of R has rank two or less.

From (12.2.15) we see that the first column of A is a multiple of u. It follows that the Givens rotations that define Q in (12.2.12) can be determined from this vector:

$$
G _ {1} \dots G _ {n - 1} u = \left[ \begin{array}{c} \tilde {u} _ {1} \\ 0 \\ \vdots \\ 0 \end{array} \right].
$$

Suppose $n = 6$ and that we have computed $G _ { 5 } , G _ { 4 }$ and $G _ { 3 }$ so that $A ^ { ( 3 ) } = G _ { 3 } G _ { 4 } G _ { 5 } A$ has the form

$$
A ^ {(3)} = \left[ \begin{array}{c c c c c c} u _ {1} v _ {1} & p _ {1} q _ {2} & p _ {1} q _ {3} & p _ {1} q _ {4} & p _ {1} q _ {5} & p _ {1} q _ {6} \\ u _ {2} v _ {1} & u _ {2} v _ {2} & p _ {2} q _ {3} & p _ {2} q _ {4} & p _ {2} q _ {5} & p _ {2} q _ {6} \\ \tilde {u} _ {3} v _ {1} & \tilde {u} _ {3} v _ {2} & \tilde {f} _ {3} g _ {3} + \tilde {h} _ {3} q _ {3} & \tilde {f} _ {3} g _ {4} + \tilde {h} _ {3} q _ {4} & \tilde {f} _ {3} g _ {5} + \tilde {h} _ {3} q _ {5} & \tilde {f} _ {3} g _ {6} + \tilde {h} _ {3} q _ {6} \\ 0 & 0 & 0 & f _ {4} g _ {4} + h _ {4} q _ {4} & f _ {4} g _ {5} + h _ {4} q _ {5} & f _ {4} g _ {6} + h _ {4} q _ {6} \\ 0 & 0 & 0 & 0 & f _ {5} g _ {5} + h _ {5} q _ {5} & f _ {5} g _ {6} + h _ {5} q _ {6} \\ 0 & 0 & 0 & 0 & 0 & f _ {6} g _ {6} + h _ {6} q _ {6} \end{array} \right].
$$

Next, we compute the cosine-sine pair $\left\{ c _ { 2 } , s _ { 2 } \right\}$ so that

$$
\tilde {G} _ {2} \left[ \begin{array}{c} u _ {2} \\ \tilde {u} _ {3} \end{array} \right] = \left[ \begin{array}{c c} c _ {2} & s _ {2} \\ - s _ {2} & c _ {2} \end{array} \right] \left[ \begin{array}{c} u _ {2} \\ \tilde {u} _ {3} \end{array} \right] = \left[ \begin{array}{c} \tilde {u} _ {2} \\ 0 \end{array} \right].
$$

Since

$$
\left[ \begin{array}{c c} c _ {2} & s _ {2} \\ - s _ {2} & c _ {2} \end{array} \right] \left[ \begin{array}{c} p _ {2} q _ {j} \\ \tilde {f} _ {3} g _ {j} + \tilde {h} _ {3} q _ {j} \end{array} \right] = \left[ \begin{array}{c} c _ {2} p _ {2} + s _ {2} \tilde {h} _ {3} \\ - s _ {2} p _ {2} + c _ {2} \tilde {h} _ {3} \end{array} \right] q _ {j} + \left[ \begin{array}{c} s _ {2} \tilde {f} _ {3} \\ c _ {2} \tilde {f} _ {3} \end{array} \right] g _ {j},
$$

for $j = 3 { : } 6$ , it follows that $A ^ { ( 2 ) } = G _ { 2 } A ^ { ( 3 ) } = \mathrm { d i a g } ( 1 , \tilde { G } _ { 2 } , I _ { 3 } ) A ^ { ( 3 ) }$ has the form

$$
A ^ {(2)} = \left[ \begin{array}{c c c c c c} u _ {1} v _ {1} & p _ {1} q _ {2} & p _ {1} q _ {3} & p _ {1} q _ {4} & p _ {1} q _ {5} & p _ {1} q _ {6} \\ \tilde {u} _ {2} v _ {1} & \tilde {f} _ {2} g _ {2} + \tilde {h} _ {2} q _ {2} & \tilde {f} _ {2} g _ {3} + \tilde {h} _ {2} q _ {3} & \tilde {f} _ {2} g _ {4} + \tilde {h} _ {2} q _ {4} & \tilde {f} _ {2} g _ {5} + \tilde {h} _ {2} q _ {5} & \tilde {f} _ {2} g _ {6} + \tilde {h} _ {2} q _ {6} \\ 0 & 0 & f _ {3} g _ {3} + h _ {3} q _ {3} & f _ {3} g _ {4} + h _ {3} q _ {4} & f _ {3} g _ {5} + h _ {3} q _ {5} & f _ {3} g _ {6} + h _ {3} q _ {6} \\ 0 & 0 & 0 & f _ {4} g _ {4} + h _ {4} q _ {4} & f _ {4} g _ {5} + h _ {4} q _ {5} & f _ {4} g _ {6} + h _ {4} q _ {6} \\ 0 & 0 & 0 & 0 & f _ {5} g _ {5} + h _ {5} q _ {5} & f _ {5} g _ {6} + h _ {5} q _ {6} \\ 0 & 0 & 0 & 0 & 0 & f _ {6} g _ {6} + h _ {6} q _ {6} \end{array} \right]
$$

where

$$
\tilde {f} _ {2} = s _ {2} \tilde {f} _ {3}, \qquad f _ {3} = c _ {2} \tilde {f} _ {3}, \qquad \tilde {h} _ {2} = c _ {2} p _ {2} + s _ {2} \tilde {h} _ {3}, \quad h _ {3} = - s _ {2} p _ {2} + c _ {2} \tilde {h} _ {3}.
$$

By considering the transition from $A ^ { ( 3 ) }$ to $A ^ { ( 2 ) }$ via the Givens rotation $G _ { 2 }$ , we conclude that $\left[ A ^ { ( 2 ) } \right] _ { 2 2 } = \tilde { u } _ { 2 } v _ { 2 }$ . Since this must equal $\tilde { f } _ { 2 } g _ { 2 } + \tilde { h } _ { 2 } q _ { 2 }$ we have

$$
g _ {2} = \frac {\tilde {u} _ {2} v _ {2} - \tilde {h} _ {2} q _ {2}}{\tilde {f} _ {2}}.
$$

By extrapolating from this example and making certain assumptions to guard against divison by zero, we obtain the following QR factorization procedure.

Algorithm 12.2.2 Suppose u, v, p, and q are n-vectors that satisfy $u . * v = p . * q$ and $u _ { n } \neq 0$ . If $A = \mathsf { t r i l } ( u \bar { v } ^ { T } ) + \mathsf { t r i u } ( p q ^ { T } , 1 )$ , then this algorithm computes cosine-sine pairs $\left\{ c _ { 1 } , s _ { 1 } \right\} , \ldots , \left\{ c _ { n - 1 } , s _ { n - 1 } \right\}$ and vectors $f , g , h \in \mathbb { R } ^ { n }$ so that if $Q$ is defined by (12.2.12) and (12.2.13), then $Q ^ { T } A = R = \mathsf { t r i u } ( f g ^ { T } + h q ^ { T } )$ .

$$
\tilde {u} _ {n} = u _ {n}, \tilde {f} _ {n} = u _ {n}, g _ {n} = v _ {n}, h _ {n} = 0
$$

for $k = n { - } 1 { : } { - } 1 { : } 1$

$\mathrm { D e t e r m i n e } \ c _ { k } \ \mathrm { a n d } \ s _ { k } \ \mathrm { s o \ t h a t } \left[ \begin{array} { c c } { { c _ { k } } } & { { s _ { k } } } \\ { { - s _ { k } } } & { { c _ { k } } } \end{array} \right] \left[ \begin{array} { c } { { u _ { k } } } \\ { { \tilde { u } _ { k + 1 } } } \end{array} \right] = \left[ \begin{array} { c } { { \tilde { u } _ { k } } } \\ { { 0 } } \end{array} \right] .$

$$
\tilde {f} _ {k} = s _ {k} \tilde {f} _ {k + 1}, f _ {k + 1} = c _ {k} \tilde {f} _ {k + 1}
$$

$$
\left[ \begin{array}{c} h _ {k} \\ h _ {k + 1} \end{array} \right] = \left[ \begin{array}{c c} c _ {k} & s _ {k} \\ - s _ {k} & c _ {k} \end{array} \right] \left[ \begin{array}{c} p _ {k} \\ h _ {k + 1} \end{array} \right]
$$

$$
g _ {k} = (u _ {k} v _ {k} - h _ {k} q _ {k}) / \tilde {f} _ {k}
$$

end

$$
f _ {1} = \tilde {f} _ {1}
$$

Regarding the condition that $u _ { n } \neq 0$ , it is easy to show by induction that

$$
\tilde {f} _ {k} = s _ {k} \dots s _ {n - 1} u _ {n}.
$$

The $s _ { k }$ are nonzero because $\lvert \tilde { u } _ { k } \rvert = \lVert \boldsymbol { u } ( \boldsymbol { k } ; \boldsymbol { n } ) \rVert _ { 2 } \neq 0$ . This algorithm requires $O ( n )$ flops and $O ( n )$ storage. We stress that there are better ways to implement the QR factorization of a semiseparable matrix than Algorithm 12.2.2. See Van Camp, Mastronardi, and Van Barel (2004). Our goal, as stated above, is to suggest how a structured rank matrix factorization can be organized around Givens rotations. Equally efficient QR factorizations for quasiseparable and semiseparable-plus-diagonal matrices are also possible.

We mention that an n-by-n system of the form triu $( f g ^ { T } + h q ^ { T } ) x = y$ can be solved in $O ( n )$ flops. An induction argument based on the partitioning

$$
\left[ \begin{array}{c c} f _ {k} g _ {k} + h _ {k} q _ {k} & f _ {k} \tilde {g} ^ {T} + h _ {1} \tilde {q} ^ {T} \\ 0 & \tilde {f} \tilde {g} ^ {T} + \tilde {h} \tilde {q} ^ {T} \end{array} \right] \left[ \begin{array}{c} x _ {k} \\ \tilde {x} \end{array} \right] = \left[ \begin{array}{c} y _ {k} \\ \tilde {y} \end{array} \right]
$$

where all the “tilde” vectors belong to $\mathbb { R } ^ { n - k }$ shows why. If ˜x, $\alpha = \tilde { g } ^ { T } \tilde { x }$ , and $\tilde { q } ^ { T } \tilde { x }$ are available, then $x _ { k }$ and the updates $\alpha = \alpha + g _ { k } x _ { k }$ and $\beta = \beta + q _ { k } x _ { k }$ require O(1) flops.
