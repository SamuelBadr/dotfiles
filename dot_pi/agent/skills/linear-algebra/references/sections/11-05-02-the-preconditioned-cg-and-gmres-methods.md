# 11.5.2 The Preconditioned CG and GMRES Methods

Before we step through the various ways that a linear system can be preconditioned, we show how the CG and GMRES iterations transform in the presence of a preconditioner. For details related to other preconditioned Krylov methods, see LIN TEMPLATES.

Suppose $M \in \mathbb { R } ^ { n \times n }$ is a symmetric positive definite matrix that we choose to regard as a preconditioner for the symmetric positive definite linear systems $A x = b$ . Recall that there is a unique symmetric positive definite matrix C such that $M = C ^ { 2 }$ . See §4.2.4. If

$$
\tilde {A} = C ^ {- 1} A C ^ {- 1}, \qquad \tilde {b} = C ^ {- 1} b,
$$

then we can solve $A x = b$ by applying CG to the symmetric positive definite system $\tilde { A } \tilde { x } = \tilde { b }$ and then solving $C x = \tilde { x }$ . For this to be a practical strategy, we must be able execute CG efficiently when it is applied to this “tilde” problem. Referring to Figure 11.3.1, here are the CG update formulae in this case:

$$
\begin{array}{l} \mu = (\tilde {r} _ {c} ^ {T} \tilde {r} _ {c}) / (\tilde {p} _ {c} ^ {T} \tilde {A} \tilde {p} _ {c}), \\ \tilde {x} _ {+} = \tilde {x} _ {c} - \mu \tilde {p} _ {c}, \\ \tilde {r} _ {+} = \tilde {r} _ {c} + \mu \tilde {A} \tilde {p} _ {c}, \tag {11.5.1} \\ \tau = (\tilde {r} _ {+} ^ {T} \tilde {r} _ {+}) / (\tilde {r} _ {c} ^ {T} \tilde {r} _ {c}), \\ \tilde {p} _ {+} = \tilde {r} _ {c} + \tau \tilde {p} _ {c}. \\ \end{array}
$$

Typically $\tilde { A }$ is dense and so we must clearly reformulate these five steps if a suitable level of efficiency is to be reached. Note that if $x _ { c } = C ^ { - 1 } \tilde { x } _ { c }$ and $r _ { c } = b - A x _ { c }$ , then

$$
\tilde {r} _ {c} = \tilde {b} - \tilde {A} \tilde {x} _ {c} = C ^ {- 1} (b - A x _ {c}) = C ^ {- 1} r _ {c}.
$$

By substituting this formula together with $\tilde { r } _ { + } = C ^ { - 1 } r _ { + }$ and the definition of $\tilde { A }$ into (11.5.1) we obtain

$$
\begin{array}{l} \mu = (r _ {c} ^ {T} M ^ {- 1} r _ {c}) / (C ^ {- 1} \tilde {p} _ {c}) ^ {T} A (C ^ {- 1} \tilde {p} _ {c}), \\ C x _ {+} = C x _ {c} - \mu \tilde {p} _ {c}, \\ C ^ {- 1} r _ {+} = C ^ {- 1} r _ {c} + \mu C ^ {- 1} A C ^ {- 1} \tilde {p} _ {c}, \\ \tau = (r _ {+} ^ {T} M ^ {- 1} r _ {+}) / (r _ {c} ^ {T} M ^ {- 1} r _ {c}), \\ \tilde {p} _ {+} = C ^ {- 1} r _ {c} + \tau \tilde {p} _ {c}. \\ \end{array}
$$

If we define $p _ { c } = C ^ { - 1 } \tilde { p } _ { c }$ and set $z _ { c } = M ^ { - 1 } r _ { c }$ , then this transforms to

Solve $\begin{array} { r } { M z _ { c } = r _ { c } , } \end{array}$

$$
\mu = (r _ {c} ^ {T} z _ {c}) / (p _ {c} ^ {T} A p _ {c}),
$$

$$
x _ {+} = x _ {c} - \mu p _ {c},
$$

$$
r _ {+} = r _ {c} + \mu A p _ {c},
$$

$$
\tau = (r _ {+} ^ {T} z _ {+}) / (r _ {c} ^ {T} z _ {c}),
$$

$$
p _ {+} = z _ {c} + \tau p _ {c},
$$

and we arrive at the method of preconditioned conjugate gradients (PCG). Note that although the square root matrix $C = M ^ { 1 / 2 }$ figured heavily in the derivation of PCG, in the end its action is felt only through the preconditioner $M = C ^ { 2 }$ .

Algorithm 11.5.1 (Preconditioned Conjugate Gradients) If $A \in \mathbb { R } ^ { n \times n }$ and $M \in \mathbb { R } ^ { n \times n }$ are symmetric positive definite, $b \in \mathbb { R } ^ { n }$ , and $A x _ { 0 } \approx b .$ , then this algorithm computes $\boldsymbol { x } _ { * } \in \mathbb { R } ^ { n }$ so that $A x _ { * } = b$ .

$$
k = 0, r _ {0} = b - A x _ {0}, \boxed {\text { Solve } M z _ {0} = r _ {0}}
$$

while $\parallel r _ { k } \parallel _ { 2 } > 0$

$$
k = k + 1
$$

$\mathbf { i f } \ k = 1$

$$
p _ {k} = z _ {0}
$$

else

$$
\tau = \left(r _ {k - 1} ^ {T} z _ {k - 1}\right) / \left(r _ {k - 2} ^ {T} z _ {k - 2}\right)
$$

$$
p _ {k} = z _ {k - 1} + \tau p _ {k - 1}
$$

end

$$
\mu = (r _ {k - 1} ^ {T} z _ {k - 1}) / (p _ {k} ^ {T} A p _ {k})
$$

$$
x _ {k} = x _ {k - 1} - \mu p _ {k}
$$

$$
r _ {k} = r _ {k - 1} - \mu A p _ {k}
$$

$$
\boxed {\text {Solve} M z _ {k} = r _ {k}}
$$

end

$$
x _ {*} = x _ {k}
$$

To highlight the difference between PCG and CG (Algorithm 11.3.2) we have boxed the preconditioner system $M z = r$ . It follows that the volume of work associated with a PCG iteration is essentially the volume of work associated with an ordinary CG iteration plus the cost of solving the preconditioner system. It can be shown that the residuals and search directions satisfy

$$
r _ {j} ^ {T} M ^ {- 1} r _ {i} = 0, \quad p _ {j} ^ {T} (C ^ {- 1} A C ^ {- 1}) p _ {i} = 0, \tag {11.5.2}
$$

for all $i \neq j$

We now turn our attention to the preconditioned GMRES method. The idea is to apply the method to the system $( M ^ { - 1 } A ) x = ( M ^ { - 1 } b )$ . Modifying Algorithm 11.4.2 in this way yields the following procedure:

Algorithm 11.5.2 (Preconditioned m-step GMRES) If $A \in \mathbb { R } ^ { n \times n }$ and $M \in \mathbb { R } ^ { n \times n }$ are nonsingular, $b \in \mathbb { R } ^ { n } , A x _ { 0 } \approx b .$ and m is a positive iteration limit, then this algorithm computes $\widetilde { \boldsymbol { x } } \in \mathbb { R } ^ { n }$ where either ˜x solves $A x = b$ or minimizes $\parallel M ^ { - 1 } ( A x - b ) \parallel _ { 2 }$ over the affine space $x _ { 0 } + \mathcal { K } ( M ^ { - 1 } A , M ^ { - 1 } r _ { 0 } , m )$ where $r _ { 0 } = b - A x _ { 0 }$ .

$$
k = 0, r _ {0} = b - A x _ {0}, \boxed {\text { Solve } M z _ {0} = r _ {0}}, \beta_ {0} = \| z _ {0} \| _ {2}
$$

while $( \beta _ { k } > 0 )$ and $k < m )$

$$
q _ {k + 1} = z _ {k} / \beta_ {k}
$$

$$
k = k + 1
$$

$$
\boxed {\text {Solve} M z _ {k} = A q _ {k}}
$$

for i = 1:k

$$
h _ {i k} = q _ {i} ^ {T} z _ {k}
$$

$$
z _ {k} = z _ {k} - h _ {i k} q _ {i}
$$

end

$$
\beta_ {k} = \left\| z _ {k} \right\| _ {2}, h _ {k + 1, k} = \beta_ {k}
$$

Apply $G _ { 1 } , \ldots , G _ { k - 1 }$ to $H ( 1 { : } k , k )$ and determine $G _ { k } , R _ { k } , p _ { k }$ , and $\rho _ { k }$ .

end

Solve $R _ { k } y _ { k } = p _ { k }$ and set $\tilde { x } = x _ { 0 } + Q _ { k } y _ { k }$

Note that $\rho _ { k } = \parallel M ^ { - 1 } ( b - A x _ { k } ) \parallel _ { 2 }$ in this formulation.
