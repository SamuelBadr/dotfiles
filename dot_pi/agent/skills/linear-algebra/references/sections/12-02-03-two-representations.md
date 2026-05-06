# 12.2.3 Two Representations

The Matlab tril and triu notation is very handy when formulating a quasiseparable matrix computation. If $A \in \mathbb { R } ^ { m \times n }$ , then $a _ { i j }$ is on its kth diagonal if $j = i + k$ . The matrix $B = { \sf t r i l } ( A , k )$ is obtained from A by setting to zero all its entries above the kth diagonal while $B = \mathsf { t r i u } ( A , k )$ is obtained from A by setting to zero all its entries below the kth diagonal. If k = 0, then we simply write tril(A) and triu(A). We also use the notation diag(d) to designate the diagonal matrix diag $( d _ { 1 } , \ldots , d _ { n } )$ where $d \in \mathbb { R } ^ { n }$ . Note that if $u , v , d , p , q \in \mathbb { R } ^ { n }$ , then the matrix

$$
A = \operatorname{tril} \left(u v ^ {T}, - 1\right) + \operatorname{diag} (d) + \operatorname{triu} \left(p q ^ {T}, 1\right) \tag {12.2.6}
$$

is quasiseparable, e.g.,

$$
A = \left[ \begin{array}{c c c c c} d _ {1} & p _ {1} q _ {2} & p _ {1} q _ {3} & p _ {1} q _ {4} & p _ {1} q _ {5} \\ u _ {2} v _ {1} & d _ {2} & p _ {2} q _ {3} & p _ {2} q _ {4} & p _ {2} q _ {5} \\ u _ {3} v _ {1} & u _ {3} v _ {2} & d _ {3} & p _ {3} q _ {4} & p _ {3} q _ {5} \\ u _ {4} v _ {1} & u _ {4} v _ {2} & u _ {4} v _ {3} & d _ {4} & p _ {4} q _ {5} \\ u _ {5} v _ {1} & u _ {5} v _ {2} & u _ {5} v _ {3} & u _ {5} v _ {4} & d _ {5} \end{array} \right].
$$

Should it be the case that $d = u . * v = p . * q$ , then this matrix is semiseparable. The representation (12.2.6) is referred to as the generator representation.

Not every quasiseparable matrix has a generator representation. For example, if $A = B ( r )$ and r has nonzero entries, then it is impossible to find u, $v , d , p , q \in \mathbb { R } ^ { n }$ so that (12.2.6) holds. To address this shortcoming, we use the fact that

$$
\binom{\text { Quasiseparable }}{\text { Matrix }}. * \binom{\text { Quasiseparable }}{\text { Matrix }} = \binom{\text { Quasiseparable }}{\text { Matrix }}, \tag {12.2.7}
$$

and embellish (12.2.6) with a pair of inverse bidiagonal factors. It can be shown that if $A \in \mathbb { R } ^ { n \times n }$ is quasiseparable, then there exist $u , v , d , p , q \in \mathbb { R } ^ { n }$ and $t , r \in \mathbb { R } ^ { n - 1 }$ such that

$$
A = \operatorname{tril} \left(u v ^ {T}, - 1\right). * B (t) ^ {- T} + \operatorname{diag} (d) + \operatorname{triu} \left(p q ^ {T}, 1\right). * B (r) ^ {- 1} \tag {12.2.8}
$$

$$
\equiv \mathbf {S} (u, v, t, d, p, q, r),
$$

e.g.,

$$
A = \left[ \begin{array}{c c c c c} d _ {1} & p _ {1} r _ {1} q _ {2} & p _ {1} r _ {1} r _ {2} q _ {3} & p _ {1} r _ {1} r _ {2} r _ {3} q _ {4} & p _ {1} r _ {1} r _ {2} r _ {3} r _ {4} q _ {5} \\ u _ {2} t _ {1} v _ {1} & d _ {2} & p _ {2} r _ {2} q _ {3} & p _ {2} r _ {2} r _ {3} q _ {4} & p _ {2} r _ {2} r _ {3} r _ {4} q _ {5} \\ u _ {3} t _ {2} t _ {1} v _ {1} & u _ {3} t _ {2} v _ {2} & d _ {3} & p _ {3} r _ {3} q _ {4} & p _ {3} r _ {3} r _ {4} q _ {5} \\ u _ {4} t _ {3} t _ {2} t _ {1} v _ {1} & u _ {4} t _ {3} t _ {2} v _ {2} & u _ {4} t _ {3} v _ {3} & d _ {4} & p _ {4} r _ {4} q _ {5} \\ u _ {5} t _ {4} t _ {3} t _ {2} t _ {1} v _ {1} & u _ {5} t _ {4} t _ {3} t _ {2} v _ {2} & u _ {5} t _ {4} t _ {3} v _ {3} & u _ {5} t _ {4} v _ {4} & d _ {5} \end{array} \right].
$$

We refer to (12.2.8) as a quasiseparable representation and it has a number of important specializations. If $d = u . * v = p . * q .$ , then A is semiseparable. If $t = r = \mathbf { 1 } _ { n - 1 }$ , then A is generator representable. If $u \ = \ q , \ v \ = \ p .$ , and $t \ = \ r$ , then A is symmetric. The representation also supports the semiseparable-plus-diagonal structure. A matrix $\mathbf { S } ( u , v , t , d , p , q , r )$ has this form if d is arbitrary and $u . * v = p . * q$ . Here are some inverse-related facts that pertain to semiseparable, quasiseparable, and diagonal-plussemiseparable matrices:

Fact 1. If A is nonsingular and tridiagonal, then $A ^ { - 1 }$ is semiseparable. In addition, if the subdiagonal and superdiagonal entries are nonzero, then $A ^ { - 1 }$ is generator-representable.

Fact 2. If A is nonsingular and quasiseparable, then so is $A ^ { - 1 }$ .

Fact 3. If $A = D + S$ is nonsingular where D is diagonal and nonsingular and S is semiseparable, then $A ^ { - 1 } = D ^ { - 1 } + S _ { 1 }$ where $S _ { 1 }$ is semiseparable.

Aspects of the first fact were encountered in §4.3.8.
