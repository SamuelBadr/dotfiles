# 11.3 The Conjugate Gradient Method

A difficulty associated with the SOR, Chebyshev semi-iterative, and related methods is that they depend upon parameters that are sometimes hard to choose properly. For example, the Chebyshev acceleration scheme requires good estimates of the largest and smallest eigenvalues of the underlying iteration matrix $M ^ { - 1 } N$ . This can be a very challenging problem unless this matrix is sufficiently structured. In this section and the next we present various Krylov subspace methods that avoid this difficulty.

We start with the well-known conjugate gradient (CG) method due to Hestenes and Stieffel (1952) and which is applicable to symmetric positive definite systems.


---

<!-- golub_650_699 -->

There are several ways to motivate and derive the technique. Our approach involves the method of steepest descent, Krylov subspaces, the Lanczos process, and tridiagonal system solving. After developing the Lanczos implementation of the CG process, we proceed to establish its equivalence with the Hestenes-Stieffel formulation.

A brief comment about notation is in order. Most of the methods in the previous section are developed at the $( i , j )$ level and this necessitated the use of superscripts to designate vector iterates. From now on, the derivations in this chapter can proceed at the vector level. Subscripts will be used to designate vector iterates, so instead of $\{ x ^ { ( k ) } \}$ we now have $\{ x _ { k } \}$ .

# 11.3.1 An Optimization Problem

Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite, $b \in \mathbb { R } ^ { n }$ , and that we want to compute the solution $x _ { * }$ to

$$
A x = b. \tag {11.3.1}
$$

Note that this problem is equivalent to solving the optimization problem

$$
\min \phi (x) \tag {11.3.2}
$$

$$
x \in \mathbb {R} ^ {n}
$$

where

$$
\phi (x) = \frac {1}{2} x ^ {T} A x - x ^ {T} b. \tag {11.3.3}
$$

This is because $\phi$ is convex and its gradient is given by

$$
\nabla \phi (x) = A x - b.
$$

Thus, if $x _ { c }$ is an approximate minimizer of $\phi ,$ then $x _ { c }$ can be regarded as an approximate solution to $A x = b$ . To make this precise, we define the A-norm by

$$
\left\| v \right\| _ {A} = \sqrt {v ^ {T} A v}. \tag {11.3.4}
$$

Since

$$
\phi (x _ {c}) = \frac {1}{2} x _ {c} ^ {T} A x _ {c} - x _ {c} ^ {T} b = \frac {1}{2} (x _ {c} - x _ {*}) A (x _ {c} - x _ {*}) - \frac {1}{2} b ^ {T} A ^ {- 1} b
$$

and $\phi ( x _ { * } ) = - b ^ { T } A ^ { - 1 } b / 2$ , it follows that

$$
\phi (x _ {c}) = \frac {1}{2} \| x _ {c} - x _ {*} \| _ {A} ^ {2} + \phi (x _ {*}). \tag {11.3.5}
$$

Thus, an iteration that produces a sequence of ever-better approximate minimizers for $\phi$ is an iteration that produces ever-better approximate solutions to $A x = b$ as measured in the A-norm.

# 11.3.2 The Method of Steepest Descent

Let us consider the minimization of $\phi$ using the method of steepest descent with exact line searches. In this method the current approximate minimizer $x _ { c }$ is improved by searching in the direction of the negative gradient, i.e., the direction of most rapid decrease. In particular, the improved approximate minimizer $x _ { + }$ is given by

$$
x _ {+} = x _ {c} - \mu_ {c} g _ {c},
$$

where $g _ { c } = A x _ { c } - b$ is the current gradient and $\mu _ { c }$ solves

$$
\min \phi (x _ {c} - \mu g _ {c}). \tag {11.3.6}
$$

$$
\mu \in \mathbb {R}
$$

This is an exact line search framework. It is easy to show that

$$
\mu_ {c} = \frac {g _ {c} ^ {T} g _ {c}}{g _ {c} ^ {T} A g _ {c}}
$$

and

$$
\phi (x _ {+}) = \phi (x _ {c}) - \frac {1}{2} \cdot \frac {(g _ {c} ^ {T} g _ {c}) ^ {2}}{r _ {c} ^ {T} A r _ {c}}. \tag {11.3.7}
$$

Thus, the objective function is decreased if $r _ { c } \neq 0$ . To establish global convergence of the method, define

$$
\kappa_ {c} = \frac {g _ {c} ^ {T} A g _ {c}}{g _ {c} ^ {T} g _ {c}} \cdot \frac {g _ {c} ^ {T} A ^ {- 1} g _ {c}}{g _ {c} ^ {T} g _ {c}}
$$

and observe that $g _ { c } ^ { T } A ^ { - 1 } g _ { c } = 2 \phi ( x _ { c } ) + b ^ { T } A ^ { - 1 } b$ and

$$
\phi (x _ {+}) = \phi (x _ {c}) - \frac {1}{2} \frac {1}{\kappa_ {c}} g _ {c} ^ {T} A ^ {- 1} g _ {c} = \phi (x _ {c}) - \frac {1}{\kappa_ {c}} \left(\phi (x _ {c}) + \frac {1}{2} b ^ {T} A ^ {- 1} b\right). \tag {11.3.8}
$$

If $\lambda _ { \mathrm { m a x } } ( A )$ and $\lambda _ { \mathrm { m i n } } ( A )$ are the largest and smallest eigenvalues of A, then we have

$$
\kappa_ {c} = \frac {g _ {c} ^ {T} A g _ {c}}{g _ {c} ^ {T} g _ {c}} \cdot \frac {g _ {c} ^ {T} A ^ {- 1} g _ {c}}{g _ {c} ^ {T} g _ {c}} \leq \frac {\lambda_ {\max} (A)}{\lambda_ {\min} (A)} = \kappa_ {2} (A).
$$

If we subtract $\phi ( x _ { * } ) = - ( b ^ { T } A ^ { - 1 } b ) / 2$ from both sides of (11.3.8) and use (11.3.5), then we obtain

$$
\left\| x _ {+} - x _ {*} \right\| _ {A} ^ {2} \leq \left(1 - \frac {1}{\kappa_ {2} (A)}\right) \left\| x _ {c} - x _ {*} \right\| _ {A} ^ {2}. \tag {11.3.9}
$$

It follows by induction that the method of steepest descent with exact line search is globally convergent.

Algorithm 11.3.1 (Steepest Descent with Exact Line Search) Given a symmetric positive definite $A \in \mathbb { R } ^ { n \times n } , \ b \in \mathbb { R } ^ { n }$ , $A x _ { 0 } \approx b .$ and a termination tolerance $\tau _ { : }$ , the following algorithm produces $\boldsymbol { x } \in \mathbb { R } ^ { n }$ so that $\| A x - b \| _ { 2 } \leq \tau$ .

$$
x = x _ {0}, g = A x - b
$$

while $\parallel g \parallel _ { 2 } > \tau$

$$
\mu = (g ^ {T} g) / (g ^ {T} A g), x = x - \mu g, g = A x - b
$$

end

Unfortunately, a convergence rate characterized by $( 1 - 1 / \kappa _ { 2 } ( A ) ) ^ { k / 2 }$ is typically not good enough unless A is extremely well-conditioned.

# 11.3.3 A Subspace Strategy

We can improve upon the steepest descent idea by expanding the dimension of the search space each step. To pursue this idea we introduce the notion of an affine space. Formally, if $v \in \mathbb { R } ^ { n }$ and $S \subseteq \mathbb { R } ^ { n }$ is a subspace, then

$$
v + S = \{x \mid x = v + s, s \in S \}.
$$

is an $a f f i n e$ space. Note that in Algorithm 11.3.1, the step-k optimization is over the affine space $x _ { k } + { \mathsf { s p a n } } \{ \nabla \phi ( x _ { k } ) \}$ .

Given $A x _ { 0 } \approx b ,$ our plan is to produce a nested sequence of subspaces

$$
S _ {1} \subset S _ {2} \subset S _ {3} \subset \dots
$$

that satisfy dim $( S _ { k } ) = k$ and to solve the problem

$$
\min _ {x \in x _ {0} + S _ {k}} \phi (x) \tag {11.3.10}
$$

each step along the way. If $x _ { k }$ is the step-k minimizer, then because of the nesting we have $\phi ( x _ { 1 } ) \geq \phi ( x _ { 2 } ) \geq \cdot \cdot \cdot \geq \phi ( x _ { n } ) = \phi ( x _ { * } )$ . Since $S _ { n } = \mathbb { R } ^ { n }$ , we ultimately obtain $x _ { * } ~ = ~ A ^ { - 1 } b .$ Even though this is a finite-step solution framework, it may not be attractive if n is extremely large. The challenge is to find a subspace sequence that promotes rapid decrease in the value of $\phi ,$ for then we may be able to terminate the iteration long before k equals n.

With this goal in mind we note that at $x _ { k }$ the function $\phi$ decreases most rapidly in the direction of the negative gradient. Thus, it makes sense to choose $S _ { k + 1 }$ so that it includes $x _ { k }$ and the gradient $g _ { k } \ = \ \nabla \phi ( x _ { k } ) \ = \ A x _ { k } - b .$ . This strategy guarantees that $x _ { k + 1 }$ is at least as good as a steepest descent update:

$$
\min _ {x \in x _ {0} + S _ {k + 1}} \phi (x) = \phi (x _ {k + 1}) \leq \min _ {\mu \in \mathbb {R}} \phi (x _ {k} - \mu g _ {k}) \tag {11.3.11}
$$

If $x _ { 0 }$ is an initial guess and we define $g _ { 0 } = A x _ { 0 } - b$ , then since $\nabla \phi ( x _ { k } ) \in \mathsf { s p a n } \{ x _ { k } , A x _ { k } \}$ it follows that the only way to satisfy this requirement is to set

$$
S _ {k} = \mathcal {K} (A, g _ {0}, k) = \operatorname{span} \left\{g _ {0}, A g _ {0}, A ^ {2} g _ {0}, \dots , A ^ {k - 1} g _ {0} \right\}.
$$

We can use the Lanczos process (§10.1) to generate these Krylov subspaces.

# 11.3.4 The Method of Conjugate Gradients: First Version

Recall that after k steps of the Lanczos iteration (Algorithm 10.1.1) we have generated a matrix

$$
Q _ {k} = \left[ q _ {1} \mid \dots \mid q _ {k} \right] \in \mathbb {R} ^ {n \times k}
$$

with orthonormal columns, a tridiagonal matrix

$$
T _ {k} = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \beta_ {k - 1} \\ 0 & \dots & & \beta_ {k - 1} & \alpha_ {k} \end{array} \right], \tag {11.3.12}
$$

and a vector $r _ { k } \in { \mathsf { r a n } } ( Q _ { k } ) ^ { \perp }$ so that

$$
A Q _ {k} = Q _ {k} T _ {k} + r _ {k} e _ {k} ^ {T}. \tag {11.3.13}
$$

Note that the tridiagonal matrix

$$
Q _ {k} ^ {T} A Q _ {k} = T _ {k}
$$

is positive definite. The solution to the optimization problem (11.3.10) via Lanczos is particularly simple if we set $q _ { 1 } = r _ { 0 } / \beta _ { 0 }$ where $r _ { 0 } = b - A x _ { 0 } = - g _ { 0 }$ , and $\beta _ { 0 } = \parallel r _ { 0 } \parallel _ { 2 } .$ . Since the columns of $Q _ { k }$ span $S _ { k } = \mathcal { K } ( A , g _ { 0 } , k )$ , it follows that the act of minimizing $\phi$ over $x _ { 0 } + S _ { k }$ is equivalent to minimizing $\phi ( x _ { 0 } + Q _ { k } y )$ over all vectors $\boldsymbol { y } \in \mathbb { R } ^ { k }$ . Since

$$
\begin{array}{l} \phi (x _ {0} + Q _ {k} y) = \frac {1}{2} (x _ {0} + Q _ {k} y) ^ {T} A (x _ {0} + Q _ {k} y) - (x _ {0} + Q _ {k} y) ^ {T} b \\ = \frac {1}{2} y ^ {T} \left(Q _ {k} ^ {T} A Q _ {k}\right) y - y ^ {T} \left(Q _ {k} ^ {T} r _ {0}\right) + \phi \left(x _ {0}\right) \\ \end{array}
$$

and $\beta _ { 0 } Q _ { k } ( : , 1 ) = r _ { 0 }$ , it follows that the minimizer $y _ { k }$ satisfies

$$
T _ {k} y _ {k} = Q _ {k} ^ {T} r _ {0} = \beta_ {0} e _ {1}
$$

and so $x _ { k } = x _ { 0 } + Q _ { k } y _ { k }$ . Building on Algorithm 10.1.1, this leads to a preliminary version of the conjugate gradient (CG) method:

$$
k = 0, r _ {0} = b - A x _ {0}, \beta_ {0} = \parallel r _ {0} \parallel_ {2}, q _ {0} = 0
$$

while $\beta _ { k } \neq 0$

$$
q _ {k + 1} = r _ {k} / \beta_ {k}
$$

$$
k = k + 1
$$

$$
\alpha_ {k} = q _ {k} ^ {T} A q _ {k} \tag {11.3.14}
$$

$$
T _ {k} y _ {k} = \beta_ {0} e _ {1}
$$

$$
x _ {k} = Q _ {k} y _ {k}
$$

$$
r _ {k} = (A - \alpha_ {k} I) q _ {k} - \beta_ {k - 1} q _ {k - 1}
$$

$$
\beta_ {k} = \left\| r _ {k} \right\| _ {2}
$$

end

$$
x _ {*} = x _ {k}
$$

As it stands, this formulation is not suitable for large problems because $x _ { k }$ is computed as an explicit n-by-k matrix-vector product and this requires access to all previously computed Lanczos vectors. However, before we develop a slick recursion for $x _ { k }$ that circumvents this problem, we establish some important properties that are associated with the iteration.

Theorem 11.3.1. If $k _ { * }$ is the dimension of the smallest invariant subspace that contains $r _ { 0 }$ , then the conjugate gradient iteration $( 1 1 . 3 . 1 \acute { 4 } )$ terminates with $x _ { k _ { * } } = x _ { * }$ .

Proof. From Theorem 10.1.1 we know that the Lanczos iteration terminates after generating $q _ { k }$ if $\kappa ( A , q _ { 1 } , k )$ is an invariant subspace. If $q _ { 1 } ~ = ~ r _ { 0 } / \parallel r _ { 0 } \parallel _ { 2 }$ , then $q _ { k , \ast }$ must be generated for otherwise $r _ { 0 }$ would be contained in an invariant subspace with dimension less than $k _ { * }$ . Since we can write $r _ { 0 }$ as a linear combination of $k _ { * }$ eigenvectors, it follows that the Krylov matrix $\left[ r _ { 0 } | A r _ { 0 } | A ^ { 2 } r _ { 0 } | \cdot \cdot \cdot | A ^ { k _ { * } } r _ { 0 } \right]$ has rank $k _ { * }$ . This implies $\beta _ { k _ { * } } = 0$ in (11.3.14) and so the iteration terminates with $x _ { * } = x _ { k * }$ .

An important ramification is that early termination can be expected if the matrix A is a low-rank perturbation of the identity matrix.

Corollary 11.3.2. Assume that $U \in \mathbb { R } ^ { n \times r } , ~ D \in \mathbb { R } ^ { r \times r }$ is symmetric, and $r < n$ . If $A = I _ { n } + U D U ^ { T }$ is positive definite and the conjugate gradient iteration $( 1 1 . 3 . 1 \acute { 4 } )$ is applied to the problem $A x = b$ , then at most $r + 1$ iterations are required to compute $x _ { * }$ .

Proof. If $v \in \mathbb { R } ^ { n }$ is in the nullspace of $U ^ { T }$ , then $A v = v$ and $\lambda = 1$ is an eigenvalue of A with multiplicity at least $n - r$ . It follows that A cannot have more than $r + 1$ distinct eigenvalues. Thus, $r _ { 0 }$ is contained in an invariant subspace with dimension $r + 1 . \quad \ \perp$

Recall that our derivation of (11.3.14) begins with a plan to improve upon the method of steepest descent. Instead of determining $x _ { k }$ from a 1-dimensional search in the direction of the $\nabla \phi ( x _ { k - 1 } )$ , the CG method determines $x _ { k }$ by searching over a Krylov subspace that includes $\nabla \phi ( x _ { k - 1 } )$ . It follows that a CG step is at least as good as a steepest descent step, as the following theorem shows.

Theorem 11.3.3. If x is the solution to the symmetric positive definite system $A x = b$ and $x _ { k }$ and $x _ { k + 1 }$ are produced by the CG method $( 1 1 . 3 . 1 \acute { 4 } )$ , then

$$
\left\| x _ {k + 1} - x _ {*} \right\| _ {A} \leq \left(1 - \frac {1}{\kappa_ {2} (A)}\right) ^ {1 / 2} \cdot \left\| x _ {k} - x _ {*} \right\| _ {A}.
$$

Proof. Setting $x _ { c } = x _ { k }$ in (11.3.9) gives

$$
\left\| x _ {+} - x _ {*} \right\| _ {A} \leq \left(1 - \frac {1}{\kappa_ {2} (A)}\right) ^ {1 / 2} \left\| x _ {k} - x _ {*} \right\| _ {A},
$$

where $x _ { + }$ is the steepest descent successor to $x _ { c }$ . By using inequality (11.3.11) we have $\left. \mathbf { \Phi } x _ { k + 1 } - x _ { * } \right. _ { A } \leq \left. \mathbf { \Phi } x _ { + } - x _ { * } \right. _ { A }$ .

Just how these mathematical results color practical matters is detailed in §11.5. For now, we continue with our exact arithmetic derivation of the method.

# 11.3.5 The Method of Conjugate Gradients: Second Version

Returning to the initial version of the CG method in (11.3.14), we work out the details associated with the tridiagonal solve $T _ { k } y _ { k } = \beta _ { 0 } e _ { 1 }$ and the matrix-vector product $x _ { k } =$ $Q _ { k } y _ { k }$ . For the overall implementation to be attractive for large sparse A, we need a way to compute $x _ { k }$ without having to access Lanczos vectors $q _ { 1 } , \ldots , q _ { k }$ . Since the tridiagonal matrix $T _ { k } = Q _ { k } ^ { T } A Q _ { k }$ is positive definite, it has an $L D L ^ { T }$ factorization. By comparing coefficients in $\ddot { T _ { k } } = L _ { k } D _ { k } L _ { k } ^ { T }$ where

$$
L _ {k} = \left[ \begin{array}{l l l l} 1 & 0 & 0 & 0 \\ \ell_ {1} & 1 & 0 & 0 \\ \vdots & \ddots & \ddots & \vdots \\ 0 & \dots & \ell_ {k - 1} & 1 \end{array} \right], \qquad D _ {k} = \left[ \begin{array}{l l l l} d _ {1} & 0 & \dots & 0 \\ 0 & d _ {2} & & \vdots \\ \vdots & & \ddots & 0 \\ 0 & \dots & 0 & d _ {k} \end{array} \right],
$$

we find

$$
d _ {1} = \alpha_ {1}
$$

for $i = 2 { : } k$

$$
\ell_ {i - 1} = \beta_ {i - 1} / d _ {i - 1} \tag {11.3.15}
$$

$$
d _ {i} = \alpha_ {i} - \ell_ {i - 1} \beta_ {i - 1}
$$

end

Given this factorization, we see that if $v _ { k } \in \mathbb { R } ^ { k }$ solves

$$
L _ {k} D _ {k} v _ {k} = \beta_ {0} e _ {1} \tag {11.3.16}
$$

then $L _ { k } ^ { T } y _ { k } = v _ { k }$ . If $C _ { k } \in \mathbb { R } ^ { n \times k }$ satisfies

$$
C _ {k} L _ {k} ^ {T} = Q _ {k}, \tag {11.3.17}
$$

then

$$
x _ {k} = x _ {0} + Q _ {k} y _ {k} = x _ {0} + C _ {k} L _ {k} ^ {T} y _ {k} = x _ {0} + C _ {k} v _ {k}. \tag {11.3.18}
$$

This is an impractical recipe because the matrix $C _ { k }$ is full and involves all the Lanczos vectors. However, there are simple connections between $C _ { k - 1 }$ and $C _ { k }$ and between $v _ { k - 1 }$ and $v _ { k }$ that can be used to transform (11.3.18) into a very handy update recipe for $x _ { k }$ . Consider the lower bidiagonal system (11.3.16), e.g.,

$$
\left[ \begin{array}{c c c c} d _ {1} & 0 & 0 & 0 \\ d _ {1} \ell_ {1} & d _ {2} & 0 & 0 \\ 0 & d _ {2} \ell_ {2} & d _ {3} & 0 \\ \hline 0 & 0 & d _ {3} \ell_ {3} & d _ {4} \end{array} \right] \left[ \begin{array}{c} \nu_ {1} \\ \nu_ {2} \\ \nu_ {3} \\ \hline \nu_ {4} \end{array} \right] = \left[ \begin{array}{c} \beta_ {0} \\ 0 \\ 0 \\ \hline 0 \end{array} \right].
$$

We conclude that

$$
v _ {k} = \left[ \begin{array}{c} \nu_ {1} \\ \vdots \\ \frac {\nu_ {k - 1}}{\nu_ {k}} \end{array} \right] = \left[ \begin{array}{c} v _ {k - 1} \\ \hline \nu_ {k} \end{array} \right] \tag {11.3.19}
$$

where

$$
\nu_ {k} = \left\{ \begin{array}{l l} \beta_ {0} / d _ {1} & \text { if   } k = 1 \\ - d _ {k - 1} \ell_ {k - 1} \nu_ {k - 1} / d _ {k} & \text { if   } k > 1 \end{array} \right.. \tag {11.3.20}
$$

Next, we consider a column partitioning of equation (11.3.17), e.g.,

$$
\left[ \begin{array}{c c c c} c _ {1} & c _ {2} & c _ {3} & c _ {4} \end{array} \right] \left[ \begin{array}{c c c c} 1 & \ell_ {1} & 0 & 0 \\ 0 & 1 & \ell_ {2} & 0 \\ 0 & 0 & 1 & \ell_ {3} \\ 0 & 0 & 0 & 1 \end{array} \right] = \left[ \begin{array}{c c c c} q _ {1} & q _ {2} & q _ {3} & q _ {4} \end{array} \right].
$$

From this we conclude that

$$
C _ {k} = \left[ \begin{array}{c c} C _ {k - 1} & c _ {k} \end{array} \right] \tag {11.3.21}
$$

where

$$
c _ {k} = \left\{ \begin{array}{l l} q _ {1} & \text { if   } k = 1 \\ q _ {k} - \ell_ {k - 1} c _ {k - 1} & \text { if   } k > 1 \end{array} \right.. \tag {11.3.22}
$$

It follows from (11.3.19) and (11.3.21) that

$$
x _ {k} = x _ {0} + C _ {k} v _ {k} = x _ {0} + C _ {k - 1} v _ {k - 1} + \nu_ {k} c _ {k} = x _ {k - 1} + \nu_ {k} c _ {k}.
$$

This is precisely the kind of recursive formula for $x _ { k }$ that we need to make the recipe (11.3.18) attractive for large sparse problems. Combining this expression with (11.3.20) and (11.3.22), we obtain the following implementation of (11.3.14).

Algorithm 11.3.2 (Conjugate Gradients: Lanczos Version) If $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite, $b \in \mathbb { R } ^ { n }$ , and $A x _ { 0 } \approx b$ , then this algorithm computes $\boldsymbol { x } _ { * } \in \mathbb { R } ^ { n }$ so that $A x _ { * } = b$ .

$$
k = 0, r _ {0} = b - A x _ {0}, \beta_ {0} = \left\| r _ {0} \right\| _ {2}, q _ {0} = 0, c _ {0} = 0
$$

while $\beta _ { k } \neq 0$

$$
q _ {k + 1} = r _ {k} / \beta_ {k}
$$

$$
k = k + 1
$$

$$
\alpha_ {k} = q _ {k} ^ {T} A q _ {k}
$$

if k = 1

$$
d _ {1} = \alpha_ {1}, \nu_ {1} = \beta_ {0} / d _ {1}
$$

$$
c _ {k} = q _ {1}
$$

$$
\ell_ {k - 1} = \beta_ {k - 1} / d _ {k - 1}, d _ {k} = \alpha_ {k} - \beta_ {k - 1} \ell_ {k - 1}, \nu_ {k} = - \beta_ {k - 1} \nu_ {k - 1} / d _ {k}
$$

$$
c _ {k} = q _ {k} - \ell_ {k - 1} c _ {k - 1}
$$

end

$$
x _ {k} = x _ {k - 1} + \nu_ {k} c _ {k}
$$

$$
r _ {k} = A q _ {k} - \alpha_ {k} q _ {k} - \beta_ {k - 1} q _ {k - 1}
$$

$$
\beta_ {k} = \left\| r _ {k} \right\| _ {2}
$$

end

$$
x _ {*} = x _ {k}
$$

Each iteration involves a single matrix-vector product and about 13n flops. It can be implemented with just a handful of length-n storage arrays as we discuss in §11.3.8.

# 11.3.6 The Gradients Are Conjugate

We make some observations about the gradients and search directions that arise during the CG iteration. First, we show that the gradients

$$
g _ {k} = A x _ {k} - b = \nabla \phi (x _ {k})
$$

are mutually orthogonal, a fact that explains the name of the algorithm.

Theorem 11.3.4. If $x _ { 1 } , \ldots , x _ { k }$ are generated by Algorithm 11.3.2, then $g _ { i } ^ { T } g _ { j } = 0$ for all i and j that satisfy $1 \leq i < j \leq k$ . Moreover, $g _ { k } = \nu _ { k } r _ { k }$ where $\nu _ { k }$ and $r _ { k }$ are defined by the algorithm.

Proof. The partial tridiagonalization (11.3.13) permits us to write

$$
g _ {k} = A x _ {k} - b = A (x _ {0} + Q _ {k} y _ {k}) - b = - r _ {0} + (Q _ {k} T _ {k} + r _ {k} e _ {k} ^ {T}) y _ {k}.
$$

Since $Q _ { k } T _ { k } y _ { k } \ : = \ : \beta _ { 0 } Q _ { k } e _ { 1 } = r _ { 0 }$ , it follows that

$$
g _ {k} = (e _ {k} ^ {T} y _ {k}) r _ {k}.
$$

Since each $r _ { i }$ is a multiple of $q _ { i + 1 }$ , it follows that the $g _ { i }$ are mutually orthogonal. To show that $g _ { k } = \nu _ { k } r _ { k }$ , we must verify that $e _ { k } ^ { T } y _ { k } = \nu _ { k }$ . From the equation

$$
T _ {k} y _ {k} = (L _ {k} D _ {k}) (L _ {k} ^ {T} y _ {k}) = \beta_ {0} e _ {1}
$$

we know that $L _ { k } ^ { T } y _ { k } = v _ { k }$ where $( L _ { k } D _ { k } ) v _ { k } = \beta _ { 0 } e _ { 1 }$ . To complete the proof, recall from (11.3.19) that $\nu _ { k }$ is the bottom component of $v _ { k }$ and exploit the fact that $L _ { k } ^ { T }$ is unit upper bidiagonal.

The search directions $c _ { 1 } , \ldots , c _ { k }$ satisfy a different kind of orthogonality property.

Theorem 11.3.5. $I f c _ { 1 } , \ldots , c _ { k }$ are generated by Algorithm 11.3.2, then

$$
c _ {i} ^ {T} A c _ {j} = \left\{ \begin{array}{l l} 0 & \text {if} i \neq j, \\ d _ {j} & \text {if} i = j, \end{array} \right.
$$

for all i and j that satisfy $1 \leq i < j \leq k$ .

Proof. Since $Q _ { k } = C _ { k } L _ { k } ^ { T }$ and $T _ { k } = Q _ { k } ^ { T } A Q _ { k }$ , we have

$$
T _ {k} = L _ {k} (C _ {k} ^ {T} A C _ {k}) L _ {k} ^ {T}.
$$

But $T _ { k } = L _ { k } D _ { k } L _ { k } ^ { T }$ and so from the uniqueness of the $L D L ^ { T }$ factorization, we have

$$
D _ {k} = C _ {k} ^ {T} A C _ {k}.
$$

The column partitioning $C _ { k } = [ c _ { 1 } \vert \dots \vert c _ { k } ]$ implies that $c _ { i } ^ { T } A c _ { j } \ = \ [ D _ { k } ] _ { i j }$ .

The theorem tells us that the search directions $c _ { 1 } , \ldots , c _ { k }$ are A-conjugate.

# 11.3.7 The Hestenes-Stiefel Formulation

The preceding results permit us to rewrite Algorithm 11.3.2 in a way that avoids explicit reference to the Lanczos vectors and the entries in the ongoing $\dot { L } D L ^ { T }$ factorization. In addition, we will be able to formulate the termination criterion in terms of the linear system residual $b - A x _ { k }$ instead of the more obscure “Lanczos residual vector” $( A - \alpha _ { k } I ) q _ { k } - \beta _ { k - 1 } q _ { k - 1 }$ . The key idea is to think of $c _ { k }$ as a search direction and $\rho _ { k }$ as a step length and to recognize that these quantities can be scaled. Consider the search direction update recipe

$$
c _ {k} = q _ {k} - \ell_ {k - 1} c _ {k - 1}
$$

from Algorithm 11.3.2. Since $q _ { k }$ is a multiple of $g _ { k - 1 }$ we see that

$$
(\text { search   direction   } k) = g _ {k - 1} + \text { scalar } \times (\text { search   direction   } k - 1)
$$

If we write this as

$$
p _ {k} = g _ {k - 1} + \tau_ {k - 1} p _ {k - 1}, \tag {11.3.23}
$$

then it follows from

$$
A p _ {k} = A g _ {k - 1} + \tau_ {k - 1} A p _ {k - 1}
$$

and Theorem 11.3.5 that

$$
\tau_ {k - 1} = - \frac {p _ {k - 1} A g _ {k - 1}}{p _ {k - 1} ^ {T} A p _ {k - 1}} \tag {11.3.24}
$$

and

$$
p _ {k} ^ {T} A g _ {k - 1} = p _ {k} ^ {T} A p _ {k}. \tag {11.3.25}
$$

Since $p _ { k }$ is a multiple of $c _ { k }$ , the update formula $x _ { k } = x _ { k - 1 } + \rho _ { k } c _ { k }$ in Algorithm 11.3.2 has the form

$$
x _ {k} = x _ {k - 1} - \mu_ {k} p _ {k}
$$

for some scalar $\mu _ { k }$ . By applying A to both sides of this equation and subtracting b we get

$$
g _ {k} = g _ {k - 1} - \mu_ {k} A p _ {k}.
$$

Using Theorem 11.3.4 and equation (11.3.25) we see that

$$
\mu_ {k} = \frac {g _ {k - 1} ^ {T} g _ {k - 1}}{g _ {k - 1} ^ {T} A p _ {k}} = \frac {g _ {k - 1} ^ {T} g _ {k - 1}}{p _ {k} ^ {T} A p _ {k}}.
$$

From the equations $g _ { k - 1 } = g _ { k - 2 } - \mu _ { k - 1 } A p _ { k - 1 }$ and $g _ { k - 1 } ^ { T } g _ { k - 2 } = 0$ , it follows that

$$
\begin{array}{l} g _ {k - 1} ^ {T} g _ {k - 1} = - \mu_ {k - 1} g _ {k - 1} ^ {T} A p _ {k - 1}, \\ g _ {k - 2} ^ {T} g _ {k - 2} = \mu_ {k - 1} g _ {k - 2} ^ {T} A p _ {k - 1} = \mu_ {k - 1} p _ {k - 1} ^ {T} A p _ {k - 1}. \\ \end{array}
$$

Substituting these equations into (11.3.24) gives

$$
\tau_ {k - 1} = \frac {g _ {k - 1} ^ {T} g _ {k - 1}}{g _ {k - 2} ^ {T} g _ {k - 2}}.
$$

By exploiting these recipes for $p _ { k } , \ x _ { k } , \ g _ { k } , \ \mu _ { k }$ , and $\tau _ { k - 1 }$ , and redefining $r _ { k }$ to be the residual $b - A x _ { k } = - g _ { k }$ , we can rewrite Algorithm 11.3.2 as follows.

Algorithm 11.3.3 (Conjugate Gradients: Hestenes-Stiefel Version) If $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite, $b \in \mathbb { R } ^ { n }$ , and $A x _ { 0 } \approx b ;$ , then this algorithm computes $\boldsymbol { x } _ { * } \in \mathbb { R } ^ { n }$ so that $A x _ { * } = b$ .

$k = 0, r_0 = b - Ax_0$ while $\| r_k\| _2 > 0$ $k = k + 1$ if $k = 1$ $p_k = r_0$ else $\tau_{k - 1} = (r_{k - 1}^T r_{k - 1}) / (r_{k - 2}^T r_{k - 2})$ $p_k = r_{k - 1} + \tau_{k - 1}p_{k - 1}$ end $\mu_{k} = (r_{k - 1}^{T}r_{k - 1}) / (p_{k}^{T}Ap_{k})$ $x_{k} = x_{k - 1} + \mu_{k}p_{k}$ $r_k = r_{k - 1} - \mu_kAp_k$

end

$$
x _ {*} = x _ {k}
$$

This procedure is essentially the form delineated in Hestenes and Stieffel (1952).

# 11.3.8 A Few Practical Details

Rounding errors lead to a loss of orthogonality among the residuals and finite termination is not guaranteed in floating point. For an extensive analysis of this fact, see Meurant (LCG). Thus, it makes sense to have a termination criterion based on (say) the size of $\parallel r _ { k } \parallel = \parallel b - A x _ { k } \parallel$ . With that in mind and being careful about required vector workspaces, we obtain the following more practical version of Algorithm 11.3.3.

$k = 0, x = x_{0}, r = b - Ax, \rho_{c} = r^{T}r, \delta = \text{tol} \cdot \parallel b \parallel_{2}$ while $\sqrt{\rho_{c}} > \delta$ $k = k + 1$ if $k = 1$ $p = r$ else $\tau = \rho_{c}/\rho_{-}, p = r + \tau p$ end $w = Ap$ $\mu = \rho_{c}/p^{T}w, x = x + \mu p, r = r - \mu w, \rho_{-} = \rho_{c}, \rho_{c} = r^{T}r$ end

Thus, a CG step requires one matrix-vector product, three saxpys, and two inner products. Four length-n arrays are required. Note that if $x _ { c }$ is the final iterate and $x _ { * }$ is the exact solution, then

$$
\| x _ {c} - x _ {*} \| = \| A ^ {- 1} (b - A x _ {c}) \| _ {2} \leq \operatorname{tol} \cdot \| A ^ {- 1} \| _ {2} \| b \| _ {2} \leq \operatorname{tol} \cdot \kappa_ {2} (A) \| x _ {*} \|.
$$

Thus, a stopping criterion ensures a relative error that is bounded by the product of tol and the condition number.

In practice, it is desirable to terminate the iteration long before k approaches n. Trefethen and Bau (NLA, p. 299) show that

$$
\| x - x _ {k} \| _ {A} \leq 2 \| x - x _ {0} \| _ {A} \left(\frac {\sqrt {\kappa_ {2} (A)} - 1}{\sqrt {\kappa_ {2} (A)} + 1}\right) ^ {k}. \tag {11.3.27}
$$

Of course, it does not take much of a condition number for the upper bound to be hopelessly close to 1, so, by itself, this result does not provide hope for an early exit. However, as we will see in §11.5, there is a way to induce speedy convergence by applying the method to an equivalent “preconditioned” system that is designed in such a way that (11.3.27) and/or Corollary 11.3.2 predict good things.

# 11.3.9 Conjugate Gradients Applied to $A ^ { T } A$ and $A A ^ { T }$

There are two obvious ways to convert an unsymmetric $A x = b$ problem into an equivalent symmetric positive definite problem:

$$
A x = b \equiv \left\{ \begin{array}{l} A ^ {T} A x = A ^ {T} b, \\ A A ^ {T} y = b, x = A ^ {T} y. \end{array} \right.
$$

Each of these conversions creates an opportunity to apply the method of conjugate gradients.

If we apply CG to the $A ^ { T } A x = A ^ { T } b$ problem, then at the kth step a vector $x _ { k }$ is produced that minimizes

$$
\phi_ {A ^ {T} A} (x) = \frac {1}{2} x ^ {T} (A ^ {T} A) x - x ^ {T} (A ^ {T} b) = \frac {1}{2} \| A x - b \| _ {2} ^ {2} - \frac {1}{2} b ^ {T} b
$$

over the affine space

$$
S _ {k} = x _ {0} + \mathcal {K} (A ^ {T} A, A ^ {T} r _ {0}, k) \tag {11.3.28}
$$

where $r _ { 0 } = b - A x _ { 0 }$ . The resulting algorithm is the conjugate gradient normal equation residual (CGNR) method.

If we apply the CG method to the “y-problem” $A A ^ { T } y = b$ , then at the kth step a vector $y _ { k }$ is produced that minimizes

$$
\phi_ {A A ^ {T}} (y) = \frac {1}{2} y ^ {T} A A ^ {T} y - y ^ {T} b = \frac {1}{2} \| A ^ {T} y - A ^ {- 1} b \| _ {2} ^ {2} - \frac {1}{2} b ^ {T} (A A ^ {T}) ^ {- 1} b
$$

over the affine space $y _ { 0 } + \mathcal { K } ( A A ^ { T } , r _ { 0 } , k )$ where $r _ { 0 } = b - A x _ { 0 }$ . Setting $x _ { k } = A ^ { T } y _ { k }$ , this says that $x = x _ { k }$ minimizes $\Vert \ b { x } - \ b { x } _ { * } \Vert _ { 2 }$ over the affine space defined in (11.3.28).

<table><tr><td>CG</td><td>CGNR</td><td>CGNE</td></tr><tr><td> $r_c = b - Ax_0$  $p_c = r_c$ </td><td> $r_c = b - Ax_0, z_c = A^T r_c$  $p_c = z_c$ </td><td> $r_c = b - Ax_c$  $p_c = A^T r_c$ </td></tr><tr><td> $\mu = \frac{r_c^T r_c}{p_c^T Ap_c}$  $x_+ = x_c + \mu p_c$  $r_+ = r_c - \mu Ap_c$  $\tau = \frac{r_+^T r_+}{r_c^T r_c}$  $p_+ = r_+ + \tau p_c$ </td><td> $\mu = \frac{z_c^T z_c}{(Ap_c)^T(Ap_c)}$  $x_+ = x_c + \mu p_c$  $r_+ = r_c - \mu Ap_c, z_+ = A^T r_+$  $\tau = \frac{z_+^T z_+}{z_c^T z_c}$  $p_+ = z_+ + \tau p_c$ </td><td> $\mu = \frac{r_c^T r_c}{p_c^T p_c}$  $x_+ = x_c + \mu p_c$  $r_+ = r_c - \mu Ap_c$  $\tau = \frac{r_+^T r_+}{r_c^T r_c}$  $p_+ = A^T r_+ + \tau p_c$ </td></tr></table>

Figure 11.3.1. The initializations and update formulae for the conjugate gradient (CG) method, the conjugate gradient normal equation residual (CGNR) method, and the conjugate gradient normal equation error (CGNE) method. The subscript $^ { 6 6 } c ^ { 9 9 }$ designates “current” while the subscript “+” designates “next”.

The resulting method is called the conjugate gradient normal equation error (CGNE) method. It is also known as Craig’s method.

Simple modifications of the CG update formulae in Algorithm 11.3.3 are required to implement CGNR and CGNE. We tabulate the initializations and updates of the three methods in Figure 11.3.1. Notice that CGNR and CGNE require procedures for A-times-vector and $A ^ { T }$ -times-vector. See Saad (IMSLS, pp. 251–254) and Greenbaum (IMSL, Chap. 7) for details and perspective on the squaring of the condition number that is associated with these methods. The CGNR method can be applied if A is rectangular. Thus, it provides a normal equation framework for solving sparse, full rank, least squares problems. See Bj¨orck (SLE, pp. 288–293) for discussion and analysis. The CGNE method can also be applied to rectangular problems, but the underlying system must be consistent.

# Problems

P11.3.1 How many n-vectors are required to implement each of the algorithms in this section?   
P11.3.2 Let $\alpha _ { i }$ and $\beta _ { i }$ be defined by Algorithm 11.3.2. How could those tridiagonal entries be generated as the iteration in Algorithm 11.3.3 proceeds?   
P11.3.3 Derive the update formulae for the CGNR and CGNE methods displayed in Figure 11.3.1.   
P11.3.4 Show that if the while-loop condition in Algorithm 11.3.3 is changed to

$$
\left\| r _ {k} \right\| > \operatorname{tol} \left(\left\| A \right\| \left\| x _ {k} \right\| + \left\| b \right\|\right),
$$

then the algorithm produces the exact solution to a nearby Ax = b problem relative to tol.

# Notes and References for §11.3

Background texts for the material in this section include Greenbaum (IMSL), Meurant (LCG), and Saad (ISPLA). The original reference for the conjugate gradient method is:

M.R. Hestenes and E. Stiefel (1952). “Methods of Conjugate Gradients for Solving Linear Systems,” J. Res. Nat. Bur. Stand. 49, 409–436.   
The idea of regarding conjugate gradients as an iterative method began with the following paper:   
J.K. Reid (1971). “ On the Method of Conjugate Gradients for the Solution of Large Sparse Systems of Linear Equations,” in Large Sparse Sets of Linear Equations, J.K. Reid (ed.), Academic Press, New York, 231–254.   
Some historical and unifying perspectives are offered in:   
G.H. Golub and D.P. O’Leary (1989). “Some History of the Conjugate Gradient and Lanczos Methods,” SIAM Review 31, 50–102.   
M.R. Hestenes (1990). “Conjugacy and Gradients,” in A History of Scientific Computing, Addison-Wesley, Reading, MA.   
S. Ashby, T.A. Manteuffel, and P.E. Saylor (1992). “A Taxonomy for Conjugate Gradient Methods,” SIAM J. Numer. Anal. 27, 1542–1568.   
Over the years, many authors have analyzed the method:   
G.W. Stewart (1975). “The Convergence of the Method of Conjugate Gradients at Isolated Extreme Points in the Spectrum,” Numer. Math. 24, 85–93.   
A. Jennings (1977). “Influence of the Eigenvalue Spectrum on the Convergence Rate of the Conjugate Gradient Method,” J. Inst. Math. Applic. 20, 61–72.   
O. Axelsson (1977). “Solution of Linear Systems of Equations: Iterative Methods,” in Sparse Matrix Techniques: Copenhagen, 1976, V.A. Barker (ed.), Springer-Verlag, Berlin.   
M.R. Hestenes (1980). Conjugate Direction Methods in Optimization, Springer-Verlag, Berlin.   
J. Cullum and R. Willoughby (1980). “The Lanczos Phenomena: An Interpretation Based on Conjugate Gradient Optimization,” Lin. Alg. Applic. 29, 63–90.   
A. van der Sluis and H.A. van der Vorst (1986). “The Rate of Convergence of Conjugate Gradients,” Numer. Math. 48, 543–560.   
A.E. Naiman, I.M. Babuka, and H.C. Elman (1997). “A Note on Conjugate Gradient Convergence,” Numer. Math. 76, 209–230.   
A.E. Naiman and S. Engelberg (2000). “A Note on Conjugate Gradient Convergence - Part II,” Numer. Math. 85, 665–683.   
S. Engelberg and A.E. Naiman (2000). “A Note on Conjugate Gradient Convergence - Part III,” Numer. Math. 85, 685–696.   
For a floating-point discussion of CG, see Meurant (LCG) as well as:   
H. Wozniakowski (1980). “Roundoff Error Analysis of a New Class of Conjugate Gradient Algorithms,” Lin. Alg. Applic. 29, 509–529.   
A. Greenbaum and Z. Strakos (1992). “Predicting the Behavior of Finite Precision Lanczos and Conjugate Gradient Computations,” SIAM J. Matrix Anal. Applic. 13, 121–137.   
Z. Strakoˇs and P. Tich´y (2002). “On Error Estimation in the Conjugate Gradient Method and Why it Works in Finite Precision Computations,” ETNA 13, 56–80.   
G. Meurant and Z. Strakoˇs (2006). “The Lanczos and Conjugate Gradient Algorithms in Finite Precision Arithmetic,” Acta Numerica 15, 471–542.   
The family of CG-related methods is very large and the following is a small subset of the literature:   
G.W. Stewart (1973). “Conjugate Direction Methods for Solving Systems of Linear Equations,” Numer. Math. 21, 284–297.   
D.P. O’Leary (1980). “The Block Conjugate Gradient Algorithm and Related Methods,” Lin. Alg. Applic. 29, 293–322.   
J.E. Dennis Jr. and K. Turner (1987). “Generalized Conjugate Directions,” Lin. Alg. Applic. 88/89, 187–209.   
A. Bunse-Gerstner and R. Stover (1999). “On a Conjugate Gradient-Type Method for Solving Complex Symmetric Linear Systems,” Lin. Alg. Applic. 287, 105–123.   
T. Barth and T. Manteuffel (2000). “Multiple Recursion Conjugate Gradient Algorithms Part I: Sufficient Conditions,” SIAM J. Matrix Anal. Applic. 21, 768–796.   
C. Li (2001). “CGNR Is an Error Reducing Algorithm,” SIAM J. Sci. Comput. 22, 2109–2112.   
A.A. Dubrulle (2001). “Retooling the Method of Block Conjugate Gradients,” ETNA 12, 216–233.

W.W. Hager and H. Zhang (2006). “Algorithm 851: CG DESCENT, a Conjugate Gradient Method with Guaranteed Descent,” ACM Trans. Math. Softw. 32, 113–137.

Y. Saad (2006). “Filtered Conjugate Residual-type Algorithms with Applications,” SIAM J. Matrix Anal. Applic. 28, 845–870.

The use of the method to solve certain eigenvalue problems is detailed in:

A. Ruhe and T. Wiberg (1972). “The Method of Conjugate Gradients Used in Inverse Iteration,” BIT 12, 543–554.

A. Edelman and S.T. Smith (1996). “On Conjugate Gradient-Like Methods for Eigen-Like Problems,” BIT 36, 494–508.

The design of sensible stopping criteria has many subtleties, see:

S.F. Ashby, M.J. Holst, A. Manteuffel, and P.E. Saylor (2001). “The Role of the Inner Product in Stopping Criteria for Conjugate Gradient Iterations,” BIT 41, 26–52.

M. Arioli (2004). “A Stopping Criterion for the Conjugate Gradient Algorithm in a Finite Element Method Framework,” Numer. Math. 97, 1–24.
