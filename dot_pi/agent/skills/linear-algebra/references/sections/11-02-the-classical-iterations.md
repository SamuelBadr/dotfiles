# 11.2 The Classical Iterations

An iterative method for the Ax = b problem generates a sequence of approximate solutions $\{ x ^ { ( k ) } \}$ that converges to $x = A ^ { - 1 } b$ . Typically, the matrix A is involved only in the context of matrix-vector multiplication and that is what makes this framework attractive when A is large and sparse. The critical attributes of an iterative method include the rate of convergence, the amount of computation per step, the volume of required storage, and the pattern of memory access. In this section, we present a collection of classical iterative methods, discuss their practical implementation, and prove a few representative theorems that illuminate their behavior.

# 11.2.1 The Jacobi and Gauss-Seidel Iterations

The simplest iterative method for the $A x = b$ problem is the Jacobi iteration. The 3-by-3 instance of the method can be motivated by rewriting the equations as follows:

$$
\begin{array}{l} x _ {1} = (b _ {1} - a _ {1 2} x _ {2} - a _ {1 3} x _ {3}) / a _ {1 1}, \\ x _ {2} = (b _ {2} - a _ {2 1} x _ {1} - a _ {2 3} x _ {3}) / a _ {2 2}, \\ x _ {3} = \left(b _ {3} - a _ {3 1} x _ {1} - a _ {3 2} x _ {2}\right) / a _ {3 3}. \\ \end{array}
$$

Suppose $x ^ { ( k - 1 ) }$ is a “current” approximation to $x = A ^ { - 1 } b$ . A natural way to generate a new approximation $x ^ { ( k ) }$ is to compute

$$
\begin{array}{l} x _ {1} ^ {(k)} = (b _ {1} - a _ {1 2} x _ {2} ^ {(k - 1)} - a _ {1 3} x _ {3} ^ {(k - 1)}) / a _ {1 1}, \\ x _ {2} ^ {(k)} = (b _ {2} - a _ {2 1} x _ {1} ^ {(k - 1)} - a _ {2 3} x _ {3} ^ {(k - 1)}) / a _ {2 2}, \tag {11.2.1} \\ x _ {3} ^ {(k)} = (b _ {3} - a _ {3 1} x _ {1} ^ {(k - 1)} - a _ {3 2} x _ {2} ^ {(k - 1)}) / a _ {3 3}. \\ \end{array}
$$

Clearly, A must have nonzeros along its diagonal for the method to be defined. For general n we have

for i = 1:n

$$
x _ {i} ^ {(k)} = \left(b _ {i} - \sum_ {j = 1} ^ {i - 1} a _ {i j} x _ {j} ^ {(k - 1)} - \sum_ {j = i + 1} ^ {n} a _ {i j} x _ {j} ^ {(k - 1)}\right) / a _ {i i} \tag {11.2.2}
$$

end

Note that the most recent solution estia particular component. For example, $x _ { 1 } ^ { ( k - 1 ) }$ s not fully exploited in the uis used in the calculation of $x _ { 2 } ^ { ( k ) }$ ing ofeven though $x _ { 1 } ^ { ( k ) }$ is available. If we revise the process so that the most current estimates of the solution components are always used, then we obtain the Gauss-Seidel iteration:

for $i = 1 { : } n$

$$
x _ {i} ^ {(k)} = \left(b _ {i} - \sum_ {j = 1} ^ {i - 1} a _ {i j} x _ {j} ^ {(k)} - \sum_ {j = i + 1} ^ {n} a _ {i j} x _ {j} ^ {(k - 1)}\right) / a _ {i i} \tag {11.2.3}
$$

end

As with Jacobi, $a _ { 1 1 } , \ldots , a _ { n n }$ must be nonzero for the iteration to be defined.

For both of these methods, the transition from $x ^ { ( k - 1 ) } { \mathrm { ~ t o ~ } } x ^ { ( k ) }$ can be succinctly described in terms of the strictly lower triangular, diagonal, and strictly upper triangular parts of the matrix A. Denote these three matrices by $L _ { A } , D _ { A }$ , and $U _ { A }$ respectively, e.g.,

$$
L _ {A} = \left[ \begin{array}{c c c} 0 & 0 & 0 \\ a _ {2 1} & 0 & 0 \\ a _ {3 1} & a _ {3 2} & 0 \end{array} \right], D _ {A} = \left[ \begin{array}{c c c} a _ {1 1} & 0 & 0 \\ 0 & a _ {2 2} & 0 \\ 0 & 0 & a _ {3 3} \end{array} \right], U _ {A} = \left[ \begin{array}{c c c} 0 & a _ {1 2} & a _ {1 3} \\ 0 & 0 & a _ {2 3} \\ 0 & 0 & 0 \end{array} \right].
$$

It is easy to show that the Jacobi step (11.2.2) has the form

$$
M _ {\mathrm{J}} x ^ {(k)} = N _ {\mathrm{J}} x ^ {(k - 1)} + b \tag {11.2.4}
$$

where $M _ { \mathrm { J } } = D _ { \mathrm { \it A } }$ and $N _ { J } ~ = ~ - ( L _ { A } + U _ { A } )$ . On the other hand, the Gauss-Seidel step (11.2.3) is defined by

$$
M _ {\mathrm{GS}} x ^ {(k)} = N _ {\mathrm{GS}} x ^ {(k - 1)} + b \tag {11.2.5}
$$

with $M _ { \mathrm { G S } } = ( D _ { A } + L _ { A } )$ and $N _ { \mathrm { G S } } = - U _ { A }$ .

# 11.2.2 Block Versions

The Jacobi and Gauss-Seidel methods have obvious block analogs. For example, if A is a 3-by-3 block matrix with square, nonsingular diagonal blocks, then the system

$$
{\left[ \begin{array}{l l l} A _ {1 1} & A _ {1 2} & A _ {1 3} \\ A _ {2 1} & A _ {2 2} & A _ {2 3} \\ A _ {3 1} & A _ {3 2} & A _ {3 3} \end{array} \right]} {\left[ \begin{array}{l} x _ {1} \\ x _ {2} \\ x _ {3} \end{array} \right]} = {\left[ \begin{array}{l} b _ {1} \\ b _ {2} \\ b _ {3} \end{array} \right]}
$$

can be rewritten as follows:

$$
A _ {1 1} x _ {1} = b _ {1} - A _ {1 2} x _ {2} - A _ {1 3} x _ {3},
$$

$$
A _ {2 2} x _ {2} = b _ {2} - A _ {2 1} x _ {1} - A _ {2 3} x _ {3},
$$

$$
A _ {3 3} x _ {3} = b _ {3} - A _ {3 1} x _ {1} - A _ {3 2} x _ {2}.
$$

From this we obtain the block Jacobi iteration

$$
A _ {1 1} x _ {1} ^ {(k)} = b _ {1} - A _ {1 2} x _ {2} ^ {(k - 1)} - A _ {1 3} x _ {3} ^ {(k - 1)},
$$

$$
A _ {2 2} x _ {2} ^ {(k)} = b _ {2} - A _ {2 1} x _ {1} ^ {(k - 1)} - A _ {2 3} x _ {3} ^ {(k - 1)},
$$

$$
A _ {3 3} x _ {3} ^ {(k)} = b _ {3} - A _ {3 1} x _ {1} ^ {(k - 1)} - A _ {3 2} x _ {2} ^ {(k - 1)},
$$

and the block Gauss-Seidel iteration

$$
A _ {1 1} x _ {1} ^ {(k)} = b _ {1} - A _ {1 2} x _ {2} ^ {(k - 1)} - A _ {1 3} x _ {3} ^ {(k - 1)},
$$

$$
A _ {2 2} x _ {2} ^ {(k)} = b _ {2} - A _ {2 1} x _ {1} ^ {(k)} - A _ {2 3} x _ {3} ^ {(k - 1)},
$$

$$
A _ {3 3} x _ {3} ^ {(k)} = b _ {3} - A _ {3 1} x _ {1} ^ {(k)} - A _ {3 2} x _ {2} ^ {(k)}.
$$

In contrast to the point versions of these iterations, a genuine linear system must be solved for $x _ { i } ^ { ( k ) }$ . These can be solved directly using LU or Cholesky factorizations or approximately solved via some iterative method. Of course, for this framework to make sense, the diagonal blocks must be nonsingular.

# 11.2.3 Splittings and Convergence

Many iterative methods for the Ax = b problem can be written in the form

$$
M x ^ {(k)} = N x ^ {(k - 1)} + b \tag {11.2.6}
$$

where $A = M - N$ is a splitting and $x ^ { ( 0 ) }$ is a starting vector. For the iteration to be practical, it must be easy to solve linear systems that involve M . This is certainly the case for the Jacobi method where M is diagonal and the Gauss-Seidel method where M is lower triangular.

It turns out that the rate of convergence associated with (11.2.6) depends on the eigenvalues of the iteration matrix

$$
G = M ^ {- 1} N.
$$

By subtracting the equation $M x = N x + b$ from (11.2.6) we obtain

$$
M (x ^ {(k)} - x) = N (x ^ {(k - 1)} - x).
$$

Thus, there is a simple connection between the error at a given step and the error at the previous step. Indeed, if

$$
e ^ {(k)} = x ^ {(k)} - x,
$$

then

$$
e ^ {(k)} = M ^ {- 1} N e ^ {(k - 1)} = G e ^ {(k - 1)} = G ^ {k} e ^ {(0)}. \tag {11.2.7}
$$

Everything hinges on the behavior of $G ^ { k }$ as $k  \infty$ . If $\| G \| < 1$ for some choice of norm, then convergence is assured because

$$
\| e ^ {(k)} \| = \| G ^ {k} e ^ {(0)} \| \leq \| G ^ {k} \| \| e ^ {(0)} \| \leq \| G \| ^ {k} \| e ^ {(0)} \|.
$$

However, it is the largest eigenvalue of G that determines the asymptotic behavior of $G ^ { k }$ . For example, if

$$
G = \left[ \begin{array}{c c} \lambda & \alpha \\ 0 & \lambda \end{array} \right],
$$

then

$$
G ^ {k} = \left[ \begin{array}{c c} \lambda^ {k} & \alpha \lambda^ {k - 1} \\ 0 & \lambda^ {k} \end{array} \right]. \tag {11.2.8}
$$

We conclude that for this problem $G ^ { k }  0$ if and only if the eigenvalue λ satisfies $| \lambda | < 1$ . Recall from (7.1.1) the definition of spectral radius:

$$
\rho (C) = \max \{| \lambda |: \lambda \in \lambda (C) \}.
$$

The following theorem links the size of $\rho ( M ^ { - 1 } N )$ to the convergence of (11.2.6).

Theorem 11.2.1. Suppose $A \ = \ M - \ N$ is a splitting of a nonsingular matrix $A \in \mathbb { R } ^ { n \times n }$ . Assuming that M is nonsingular, the iteration $( 1 1 . 2 . 6 )$ converges to $x =$ $A ^ { - 1 } b$ for all starting n-vectors $x ^ { ( 0 ) }$ if and only if $\rho ( G ) < 1$ where $\dot { G } = M ^ { - 1 } N$ .

Proof. In light of (11.2.7), it suffices to show that $G ^ { k }  0$ if and only if $\rho ( G ) < 1$ . If $G x = \lambda x$ , then $G ^ { \dot { k } } x = \dot { \lambda } ^ { k } x$ . Thus, if $G ^ { k }  0$ , then we must have $| \lambda | < 1 , \mathrm { i . e . }$ , the spectral radius of G must be less than 1.

Now assume $\rho ( G ) ~ < ~ 1$ and let $G = Q T Q ^ { H }$ be its Schur decomposition. If $D = \operatorname { d i a g } ( t _ { 1 1 } , \dots , t _ { n n } )$ and $E = A - D$ , then it follows from (7.3.15) that

$$
\| G ^ {k} \| _ {2} \leq (1 + \mu) ^ {n - 1} \left(\rho (G) + \frac {\| E \| _ {F}}{1 + \mu}\right) ^ {k}
$$

where $\mu$ is any nonnegative real number. It is clear that we can choose this parameter so that the upper bound converges to zero. For example, if G is normal, then $E = 0$ and we can set $\mu = 0$ . Otherwise, if

$$
\mu = \frac {2 \| E \| _ {2}}{1 - \rho (G)},
$$

then it is easy to verify that

$$
\| G ^ {k} \| _ {2} \leq \left(1 + \frac {2 \| E \| _ {F}}{1 - \rho (G)}\right) ^ {n - 1} \left(\frac {1 + \rho (G)}{2}\right) ^ {k} \tag {11.2.9}
$$

and this guarantees convergence because $1 + \rho ( G ) < 2$ .

The 2-by-2 example (11.2.8) and the inequality (11.2.9) serve as a reminder that the spectral radius does not tell us everything about the powers of a nonnormal matrix. Indeed, if G is nonnormal, then is possible for $G ^ { k }$ (and the error $\parallel x ^ { ( k ) } - x \parallel )$ to grow considerably before decay sets in. The 
-pseudospectral radius introduced in §7.9.6 provides greater insight into this situation.

To summarize what we have learned so far, two attributes are critical if a method of the form (11.2.6) is to be of interest:

• The underlying splitting $A = M - N$ must have the property that linear systems of the form $M z = d$ are relatively easy to solve.

• A way must be found to guarantee that $\rho ( M ^ { - 1 } N ) < 1$ .

To give a flavor for the kind of analysis that attends the second requirement, we state and prove a pair of convergence results that apply to the Jacobi and Gauss-Seidel iterations.

# 11.2.4 Diagonal Dominance and Jacobi Iteration

One way to establish that the spectral radius of the iteration matrix G is less than one is to show that $\| G \| < 1$ for some choice of norm. This inequality ensures that all of $G ^ { \prime }$ s eigenvalues are inside the unit circle. As an example of this type of analysis, consider the situation where the Jacobi iteration is applied to a strictly diagonally dominant linear system. Recall from §4.1.1 that $A \in \mathbb { R } ^ { n \times n }$ has this property if

$$
\sum_{\substack{j = 1\\ j\neq i}}^{n}|a_{ij}| <   |a_{ii}|,\qquad i = 1:n.
$$

Theorem 11.2.2. If $A \in \mathbb { R } ^ { n \times n }$ is strictly diagonally dominant, then the Jacobi itreation $( 1 1 . 2 . \not \to )$ converges to $x = A ^ { - 1 } b$ .

Proof. Since $G _ { \mathrm { J } } = - D _ { A } ^ { - 1 } ( L _ { A } + U _ { A } )$ it follows that

$$
\| G_{\mathrm{J}}\|_{\infty} = \| D_{A}^{-1}(L_{A} + U_{A})\|_{\infty} = \max_{1\leq i\leq n}\sum_{\substack{j = 1\\ j\neq i}}^{n}\left|\frac{a_{ij}}{a_{ii}}\right| <   1.
$$

The theorem follows because no eigenvalue of A can be bigger that $\| A \| _ { \infty }$

Usually, the “more dominant” the diagonal the more rapid the convergence, but there are counterexamples. See P11.2.3.

# 11.2.5 Positive Definiteness and Gauss-Seidel Iteration

A more complicated spectral radius argument is needed to show that Gauss-Seidel converges for matrices that are symmetric positive definite.

Theorem 11.2.3. If $A \in \mathbb { R } ^ { n \times n }$ is symmetric and positive definite, then the Gauss-Seidel iteration (11.2.5) converges for any $x ^ { ( 0 ) }$ .

Proof. We must verify that the eigenvalues of $G _ { \mathrm { G S } } = - ( D _ { A } + L _ { A } ) ^ { - 1 } L _ { A } ^ { T }$ are inside the unit circle. This matrix has the same eigenvalues as the matrix

$$
G = D _ {A} ^ {1 / 2} G _ {\mathrm{GS}} D _ {A} ^ {- 1 / 2} = - (I + L) ^ {- 1} L ^ {T}
$$

where $L = D _ { A } ^ { - 1 / 2 } L _ { A } D _ { A } ^ { - 1 / 2 }$ . If

$$
- (I + L) ^ {- 1} L ^ {T} v = \lambda v \quad v ^ {H} v = 1
$$

$\begin{array} { r } { \mathrm { t h e n } \ - v ^ { H } L ^ { H } v = \lambda ( 1 + v ^ { H } L v ) . \ \mathrm { I f } \ v ^ { H } L v = a + b i , \mathrm { t h e n } } \end{array}$

$$
| \lambda | ^ {2} = \left| \frac {- a + b i}{1 + a + b i} \right| ^ {2} = \frac {a ^ {2} + b ^ {2}}{1 + 2 a + a ^ {2} + b ^ {2}}.
$$

However, since $D _ { A } ^ { - 1 / 2 } A D _ { A } ^ { - 1 / 2 } = I + L + L ^ { T }$ is positive definite, it is not hard to show that $\begin{array} { r } { 0 < 1 + v ^ { H } L v + v ^ { H } L ^ { T } v = 1 + 2 a } \end{array}$ and hence that $| \lambda | < 1$ .

We mention that to bound $\rho ( M _ { \mathrm { G S } } ^ { - 1 } N _ { \mathrm { G S } } )$ away from 1 requires additional information about A. The required analysis can be quite involved.

# 11.2.6 Discussion of a Model Problem

It is instructive to consider application of the Jacobi and Gauss-Seidel methods to the symmetric positive definite linear system

$$
\left(I _ {n _ {1}} \otimes T _ {n _ {2}} + T _ {n _ {1}} \otimes I _ {n _ {2}}\right) u = b \tag {11.2.10}
$$

where

$$
T _ {m} = \left[ \begin{array}{c c c c} 2 & - 1 & \dots & 0 \\ - 1 & 2 & \ddots & \vdots \\ \vdots & \ddots & \ddots & - 1 \\ 0 & \dots & - 1 & 2 \end{array} \right] \in \mathbb {R} ^ {m \times m}. \tag {11.2.11}
$$

Systems with this structure arise from discretization of the Poisson equation on a rectangular grid; see §4.8.3. Recall that it is convenient to think of the solution vector as doubly subscripted. Associated with grid point $( i , j )$ is the unknown $U ( i , j )$ . When the system is solved, the value of $U ( i , j )$ is the average of the values associated with its north, east, south, and west “grid neighbors.” Boundary values are known and fixed and this permits us to reformulate (11.2.10) as a 2-dimensional array averaging problem:

Given $U ( 0 ; n _ { 1 } + 1 , 0 ; n _ { 2 } + 1 )$ with fixed values in its top and bottom row and fixed values in its leftmost and rightmost columns, determine $U ( 1 { : } n _ { 1 } , 1 { : } n _ { 2 } )$ such that

$$
U (i, j) = \frac {U (i , j - 1) + U (i , j + 1) + U (i - 1 , j) + U (i + 1 , j)}{4}
$$

for i = 1:n1 and j = 1:n2.

It is much easier to reason about Jacobi and Gauss-Seidel from this point of view. For example, the update

$V = U$ for $i = 1:n_{1}$ for $j = 1:n_{2}$ $U(i,j) = (V(i - 1,j) + V(i,j + 1) + V(i + 1,j) + V(i,j - 1)) / 4$ end   
end

corresponds to one step of Jacobi while

for $i = 1:n_{1}$ for $j = 1:n_{2}$ $U(i,j) = (U(i-1,j) + U(i,j+1) + U(i+1,j) + U(i,j-1))/4$ end
end

is the corresponding update associated with Gauss-Seidel. The organization of both methods reflects the ultimate exploitation of matrix structure: The matrix A is nowhere in sight! We simply take advantage of the Kronecker structure at the block level and the 1-2-1 structure of the underlying tridiagonal matrices.

The array-update point of view for the model problem that we are considering makes it easy to appreciate why the Jacobi process is typically easier to vectorize and/or parallelize than Gauss-Seidel. The Jacobi update of $U ( 1 { : } n _ { 1 } , 1 { : } n _ { 2 } )$ is a matrix averaging:

$$
\frac {U (1 : n _ {1} , 0 : n _ {2} - 1) + U (2 : n _ {1} + 1 , 1 : n _ {2}) + U (1 : n _ {1} , 2 : n _ {2} + 1) + U (0 : n _ {1} - 1 , 1 : n _ {2})}{4}.
$$

The use-the-most-recent-estimate attribute of the Gauss-Seidel method makes it harder to describe the update at such a high level.

Now let us analyze the spectral radius $\rho ( M _ { \mathrm { J } } ^ { - 1 } N _ { \mathrm { J } } )$ . Closed-form expressions for $T _ { m } \mathrm { { ^ { ~ \circ ~ } s } }$ eigenvalues permit us to determine this important quantity. Note that

$$
\mathcal {T} _ {m} = 2 I - E _ {m}
$$

where

$$
E _ {m} = \left[ \begin{array}{c c c c} 0 & 1 & \dots & 0 \\ 1 & 0 & \ddots & \vdots \\ \vdots & \ddots & \ddots & 1 \\ 0 & \dots & 1 & 0 \end{array} \right].
$$

Since

$$
A = I _ {n _ {1}} \otimes \mathcal {T} _ {n _ {2}} + T _ {n _ {1}} \otimes I _ {n _ {2}} = 4 I _ {n _ {1} n _ {2}} - (I _ {n _ {1}} \otimes E _ {n _ {2}}) - (E _ {n _ {1}} \otimes I _ {n _ {2}}), \tag {11.2.12}
$$

the Jacobi splitting $A = M _ { \mathrm { J } } - N _ { \mathrm { J } }$ is given by

$$
\begin{array}{l} M _ {\mathrm{J}} = 4 I _ {n _ {1} n _ {2}}, \\ N _ {\mathrm{J}} = \left(I _ {n _ {1}} \otimes E _ {n _ {2}}\right) + \left(E _ {n _ {1}} \otimes I _ {n _ {2}}\right). \\ \end{array}
$$

Using results from our fast eigensystem discussion in §4.8.6, it can be shown that

$$
S _ {m} ^ {- 1} E _ {m} S _ {m} = D _ {m} = \operatorname{diag} \left(\mu_ {1} ^ {(m)}, \dots , \mu_ {m} ^ {(m)}\right) \tag {11.2.13}
$$

where $S _ { m }$ is the sine transform matrix $[ S _ { m } ] _ { k j } = \sin ( k j \pi / ( m + 1 ) )$ and

$$
\mu_ {k} ^ {(m)} = 2 \cos \left(\frac {k \pi}{m + 1}\right), \quad k = 1: m. \tag {11.2.14}
$$

It follows that

$$
\left(S _ {n _ {1}} \otimes S _ {n _ {2}}\right) ^ {- 1} \left(M _ {\mathrm{J}} ^ {- 1} N _ {\mathrm{J}}\right) \left(S _ {n _ {1}} \otimes S _ {n _ {2}}\right) = \left(I _ {n _ {1}} \otimes D _ {n _ {2}} + D _ {n _ {1}} \otimes I _ {n _ {2}}\right) / 4.
$$

By using the Kronecker structure of this diagonal matrix and (11.2.14), it is easy to verify that

$$
\rho (M _ {\mathrm{J}} ^ {- 1} N _ {\mathrm{J}}) = \frac {2 \cos (\pi / (n _ {1} + 1)) + 2 \cos (\pi / (n _ {2} + 1))}{4}. \tag {11.2.15}
$$

Note that this quantity approaches unity as $n _ { 1 }$ and $n _ { 2 }$ increase.

As a final exercise concerning the model problem, we use its special structure to develop an interesting alternative iteration. From (11.2.12) we can write $A = M _ { x } - N _ { x }$ where

$$
M _ {x} = 4 I _ {n _ {1} n _ {2}} - (I _ {n _ {1}} \otimes E _ {n _ {2}}), \quad N _ {x} = (E _ {n _ {1}} \otimes I _ {n _ {2}}).
$$

Likewise, $A = M _ { y } - N _ { y }$ where

$$
M _ {y} = 4 I _ {n _ {1} n _ {2}} - (E _ {n _ {1}} \otimes I _ {n _ {2}}), \qquad N _ {y} = (I _ {n _ {1}} \otimes E _ {n _ {2}}).
$$

These two splittings can be paired to produce the following transition from $u ^ { ( k - 1 ) }$ to $u ^ { ( k ) }$ :

$$
\begin{array}{r c l} M _ {x} v ^ {(k)} & = & N _ {x} u ^ {(k - 1)} + b, \\ M _ {x} (k) & = & N _ {x} (k) + b \end{array} \tag {11.2.16}
$$

$$
M _ {y} u ^ {(k)} = N _ {y} v ^ {(k)} + b.
$$

Each step has a natural interpretation based on the underlying partial differential equation; see §4.8.4. The first step corresponds to treating the north and south values at each grid point as fixed, while the second step corresponds to treating the east and west values at each grid point as fixed. The resulting iteration is an example of an alternating direction iteration. See Varga (1962, Chap. 7). Since

$$
u ^ {(k)} - x = (M _ {y} ^ {- 1} N _ {y}) (v ^ {(k)} - x) = (M _ {y} ^ {- 1} N _ {y}) (M _ {x} ^ {- 1} N _ {x}) (u ^ {(k - 1)} - x)
$$

it follows that $e ^ { ( k ) } = G ^ { k } e ^ { ( 0 ) }$ where

$$
\begin{array}{l} G = (M _ {y} ^ {- 1} N _ {y}) (M _ {x} ^ {- 1} N _ {x}) \\ = \left(4 I _ {n _ {1} n _ {2}} - E _ {n _ {1}} \otimes I _ {n _ {2}}\right) ^ {- 1} \left(I _ {n _ {1}} \otimes E _ {n _ {2}}\right) \left(4 I _ {n _ {1} n _ {2}} - I _ {n _ {1}} \otimes E _ {n _ {2}}\right) ^ {- 1} \left(E _ {n _ {1}} \otimes I _ {n _ {2}}\right). \\ \end{array}
$$

Using (11.2.13) and (11.2.14) it is easy to show that

$$
\left(S _ {n _ {1}} \otimes S _ {n _ {2}}\right) ^ {- 1} G \left(S _ {n _ {1}} \otimes S _ {n _ {2}}\right) =
$$

$$
(4 I _ {n _ {1} n _ {2}} - D _ {n _ {1}} \otimes I _ {n _ {2}}) ^ {- 1} (I _ {n _ {1}} \otimes D _ {n _ {2}}) (4 I _ {n _ {1} n _ {2}} - I _ {n _ {1}} \otimes D _ {n _ {2}}) ^ {- 1} (D _ {n _ {1}} \otimes I _ {n _ {2}})
$$

is diagonal and that

$$
\rho (G) = = \frac {\cos (\pi / (n _ {1} + 1)) \cos (\pi / (n _ {2} + 1))}{(2 - \cos (\pi / (n _ {1} + 1)) (2 - \cos (\pi / (n _ {2} + 1)))} <   1. \tag {11.2.17}
$$

# 11.2.7 SOR and Symmetric SOR

The Gauss-Seidel iteration is very attractive because of its simplicity. Unfortunately, if the spectral radius of $M _ { \mathrm { G S } } ^ { - 1 } { N } _ { \mathrm { G S } }$ is close to unity, then it may be prohibitively slow. To address this concern, we consider the parameterized splitting $A = M _ { \omega } - N _ { \omega }$ where

$$
M _ {\omega} = \frac {1}{\omega} D _ {A} + L _ {A} \quad N _ {\omega} = \left(\frac {1}{\omega} - 1\right) D _ {A} + U _ {A}. \tag {11.2.18}
$$

This defines the method of successive over-relaxation (SOR):

$$
\left(\frac {1}{\omega} D _ {A} + L _ {A}\right) x ^ {(k)} = \left(\left(\frac {1}{\omega} - 1\right) D _ {A} + U _ {A}\right) x ^ {(k - 1)} + b. \tag {11.2.19}
$$

At the component level we have

for $i = 1 { : } n$

$$
x _ {i} ^ {(k)} = \omega \left(b _ {i} - \sum_ {j = 1} ^ {i - 1} a _ {i j} x _ {j} ^ {(k)} - \sum_ {j = i + 1} ^ {n} a _ {i j} x _ {j} ^ {(k - 1)}\right) / a _ {i i} + (1 - \omega) x _ {i} ^ {(k - 1)}
$$

end

Note that if $\omega = 1$ , then this is just the Gauss-Seidel method. The idea is to choose ω so that $\rho ( M _ { \omega } ^ { - 1 } N _ { \omega } )$ is minimized. A detailed theory on how to do this is developed by Young (1971). For an excellent synopsis of that theory, see Greenbaum (IMSL, p. 149).

Observe that x is updated top to bottom in the SOR step. We can just as easily update from bottom to top:

for $i = n \colon - 1 \colon 1$

$$
x _ {i} ^ {(k)} = \omega \left(b _ {i} - \sum_ {j = 1} ^ {i - 1} a _ {i j} x _ {j} ^ {(k - 1)} - \sum_ {j = i + 1} ^ {n} a _ {i j} x _ {j} ^ {(k)}\right) / a _ {i i} + (1 - \omega) \cdot x _ {i} ^ {(k - 1)}
$$

end

This defines the backward SOR iteration:

$$
\left(\frac {1}{\omega} D _ {A} + U _ {A}\right) x ^ {(k)} = \left(\left(\frac {1}{\omega} - 1\right) D _ {A} + L _ {A}\right) x ^ {(k - 1)} + b. \tag {11.2.21}
$$

Note that this update can be obtained from (11.2.19) simply by interchanging the roles of L and U .

If A is symmetric $( U _ { A } = L _ { A } ^ { T } )$ , then the symmetric SOR (SSOR) method is obtained by combining the forward and backward implementations of the update as follows:

$$
\left(\frac {1}{\omega} D _ {A} + L _ {A}\right) y ^ {(k)} = \left(\left(\frac {1}{\omega} - 1\right) D _ {A} - L _ {A} ^ {T}\right) x ^ {(k - 1)} + b, \tag {11.2.22}
$$

$$
\left(\frac {1}{\omega} D _ {A} + L _ {A} ^ {T}\right) x ^ {(k)} = \left(\left(\frac {1}{\omega} - 1\right) D _ {A} - L _ {A}\right) y ^ {(k)} + b. \tag {11.2.23}
$$

It can be shown that if

$$
M _ {\text {SSOR}} = \frac {\omega}{2 - \omega} \left(\frac {1}{\omega} D _ {A} + L _ {A}\right) D _ {A} ^ {- 1} \left(\frac {1}{\omega} D _ {A} + L _ {A} ^ {T}\right) \tag {11.2.24}
$$

then the transition from $x ^ { ( k - 1 ) } \ \mathrm { t o } \ x ^ { ( k ) }$ is given by

$$
x ^ {(k)} = x ^ {(k - 1)} + M _ {\mathrm{SSOR}} ^ {- 1} (b - A x ^ {(k - 1)}). \tag {11.2.25}
$$

Note that $M _ { S S O R }$ is defined if $0 < \omega < 2$ and that it is symmetric. It is also positive definite if A has positive diagonal entries. Here is a result that shows SSOR converges if A is symmetric and positive definite.

Theorem 11.2.4. Suppose the SSOR method (11.2.22) and (11.2.23) is applied to a symmetric positive definite Ax = b problem and that $0 < \omega < 2$ . If

$$
M _ {\omega} = \frac {1}{\omega} D _ {A} + L _ {A}, \qquad N _ {\omega} = \left(\frac {1}{\omega} - 1\right) D _ {A} - L _ {A} ^ {T}, \qquad G = M _ {\omega} ^ {- T} N _ {\omega} ^ {T} M _ {\omega} ^ {- 1} N _ {\omega},
$$

then G has real eigenvalues, $\rho ( G ) < 1$ , and

$$
(x ^ {(k)} - x) = G ^ {k} (x ^ {(0)} - x). \tag {11.2.26}
$$

Proof. From (11.2.22) and (11.2.23) it follows that

$$
y ^ {(k)} - x = M _ {\omega} ^ {- 1} N _ {\omega} (x ^ {(k - 1)} - x),
$$

$$
x ^ {(k)} - x = M _ {\omega} ^ {- T} N _ {\omega} ^ {T} (y ^ {(k)} - x),
$$

from which it is easy to verify (11.2.26). Since D is a diagonal matrix with positive diagonal entries, there is a diagonal matrix $D _ { 1 }$ so $D = D _ { 1 } ^ { \bar { 2 } }$ . If ${ \cal L } _ { 1 } = D _ { 1 } ^ { - 1 } L \bar { D _ { 1 } ^ { - 1 } }$ and $G _ { 1 } = D _ { 1 } G D _ { 1 } ^ { - 1 }$ , then with a little manipulation we have

$$
G _ {1} = (I + \omega L _ {1} ^ {T}) ^ {- 1} (I + \omega L _ {1}) ^ {- 1} ((1 - \omega) I - \omega L _ {1}) ((1 - \omega) I - \omega L _ {1} ^ {T}).
$$

We show that if $\lambda \in \lambda ( G _ { 1 ) }$ , then $0 \leq \lambda < 1$ . If $G _ { 1 } v = \lambda v$ , then

$$
((1 - \omega) I - \omega L _ {1}) ((1 - \omega) I - \omega L _ {1} ^ {T}) v = \lambda (I + \omega L _ {1}) (I + \omega L _ {1} ^ {T}) v.
$$

This is a generalized singular value problem; see §8.7.4. It follows that λ is real and nonnegative. Assuming that $v \in \mathbb { R } ^ { n }$ has unit 2-norm, it is easy to show that

$$
\lambda = \frac {\| (1 - \omega) v - \omega L _ {1} ^ {T} v \| _ {2} ^ {2}}{\| v + \omega L _ {1} ^ {T} v \| _ {2} ^ {2}} = 1 - \omega (2 - \omega) \frac {1 + 2 v ^ {T} L _ {1} ^ {T} v}{\| v + \omega L _ {1} ^ {T} v \| _ {2} ^ {2}}. \tag {11.2.27}
$$

To complete the proof, note that $1 + 2 v ^ { T } L _ { 1 } ^ { T } v = ( D _ { 1 } ^ { - 1 } v ) ^ { T } A ( D _ { 1 } ^ { - 1 } v )$ and that this quantity is positive. By hypothesis, $\omega ( 2 - \omega ) > 0$ and so we have $\lambda < 1$

The original analysis of the symmetric SOR method is in Young (1970).

# 11.2.8 The Chebyshev Semi-Iterative Method

Another way to accelerate the convergence of certain iterative methods makes use of Chebyshev polynomials. Suppose the iteration $M x ^ { ( j + 1 ) } = N x ^ { ( j ) } + b$ has been used to generate $x ^ { ( 1 ) } , \ldots , x ^ { ( k ) }$ and that we wish to determine coefficients $\nu _ { j } ( k ) , j = 0 { : } k$ such that

$$
y ^ {(k)} = \sum_ {j = 0} ^ {k} \nu_ {j} (k) x ^ {(j)} \tag {11.2.28}
$$

represents an improvement over $x ^ { ( k ) }$ . If $x ^ { ( 0 ) } = \cdot \cdot \cdot = x ^ { ( k ) } = x$ , then it is reasonable to insist that $y ^ { ( k ) } = x$ . If the polynomial

$$
p _ {k} (z) = \sum_ {j = 0} ^ {k} \nu_ {j} (k) z ^ {j}
$$

satisfies $p _ { k } ( 1 ) = 1$ , then this criterion is satisfied and

$$
y ^ {(k)} - x = \sum_ {j = 0} ^ {k} \nu_ {j} (k) (x ^ {(j)} - x) = \sum_ {j = 0} ^ {k} \nu_ {j} (k) (M ^ {- 1} N) ^ {j} e ^ {(0)} = p _ {k} (G) e ^ {(0)}
$$

where $G = M ^ { - 1 } N$ . By taking norms in this equation we obtain

$$
\| y ^ {(k)} - x \| _ {2} \leq \| p _ {k} (G) \| _ {2} \| e ^ {(0)} \| _ {2}. \tag {11.2.29}
$$

This suggests that we can produce an improved approximate solution if we can find a polynomial $p _ { k } ( \cdot )$ that (a) has degree k, (b) satisfies $p _ { k } ( 1 ) = 1$ , and (c) does a good job of minimizing the upper bound.

To implement this idea, we assume for simplicity that G is symmetric. (There are ways to proceed if this is not the case; see Manteuffel (1977). Let

$$
S ^ {T} G S = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right) = \Lambda
$$

be a Schur decomposition of G and assume that

$$
- 1 <   \alpha \leq \lambda_ {n} \leq \dots \leq \lambda_ {1} \leq \beta <   1 \tag {11.2.30}
$$

where α and $\beta$ are known estimates. It follows that

$$
\| p _ {k} (G) \| _ {2} = \| p _ {k} (\Lambda) \| _ {2} = \max _ {\lambda_ {i} \in \lambda (A)} | p _ {k} (\lambda_ {i}) | \leq \max _ {\alpha \leq \lambda \leq \beta} | p _ {k} (\lambda) |.
$$

The degree-k Chebyshev polynomial $c _ { k } ( \cdot )$ can be used to design a good choice for $p _ { k } ( \cdot )$ . We want a polynomial whose value on $[ \alpha , \beta ]$ is small subject to the constraint that $p _ { k } ( 1 ) = 1$ . Recall from the discussion in §10.1.5 that the Chebyshev polynomials are bounded by unity on $[ - 1 , + 1 ]$ , but that their value is very large outside this range. As a consequence, if

$$
\mu = - 1 + 2 \frac {1 - \alpha}{\beta - \alpha} = 1 + 2 \frac {1 - \beta}{\beta - \alpha},
$$

then the polynomial

$$
p _ {k} (z) = c _ {k} \left(- 1 + 2 \frac {z - \alpha}{\beta - \alpha}\right) / c _ {k} (\mu)
$$

satisfies $p _ { k } ( 1 ) = 1$ and is bounded by $1 / | c _ { k } ( \mu ) |$ on $[ \alpha , \beta ]$ . From the definition of $p _ { k } ( z )$ and inequality (11.2.29) we see

$$
\| y ^ {(k)} - x \| _ {2} \leq \frac {\| x - x ^ {(0)} \| _ {2}}{| c _ {k} (\mu) |}.
$$

The larger the value of $\mu$ the greater the acceleration of convergence.

In order for the whole process to be effective, we need a more efficient method for calculating $y ^ { ( k ) }$ than (11.2.28). The retrieval of the vectors $x ^ { ( 0 ) } , \ldots , x ^ { ( k ) }$ becomes an unacceptable overhead as k increases. Fortunately, it is possible to derive a three-term recurrence among the $y ^ { ( k ) }$ by exploiting the three-term recurrence that exists among the Chebyshev polynomials. Assume (for simplicity) that $\alpha = - \beta$ in (11.2.30) and that we are given $\boldsymbol { x } ^ { ( 0 ) } \in \mathbb { R } ^ { n }$ . Here is how the process plays out when it is used to accelerate the iteration $M x ^ { ( j + 1 ) } = N x ^ { ( j ) } + b \colon$

$$
c _ {0} = 1; c _ {1} = 1 / \beta
$$

$$
y ^ {(0)} = x ^ {(0)}, M y ^ {(1)} = N y ^ {(0)} + b, r ^ {(1)} = b - A y ^ {(1)}, k = 1
$$

while $\parallel r ^ { ( k ) } \parallel >$ tol

$$
c _ {k + 1} = (2 / \beta) c _ {k} - c _ {k - 1}
$$

$$
\omega_ {k + 1} = 1 + c _ {k - 1} / c _ {k + 1}
$$

$$
M z ^ {(k)} = r ^ {(k)}
$$

$$
y ^ {(k + 1)} = y ^ {(k - 1)} + \omega_ {k + 1} \left(y ^ {(k)} + z ^ {(k)} - y ^ {(k - 1)}\right)
$$

$$
k = k + 1
$$

$$
r ^ {(k)} = b - A y ^ {(k)}
$$

end

Note that $y ^ { ( 0 ) } = x ^ { ( 0 ) }$ and $y ^ { ( 1 ) } = x ^ { ( 1 ) }$ , but that thereafter the $x ^ { ( k ) }$ are not involved. For the acceleration to be effective we need good lower and upper bounds in (11.2.30) and that is sometimes difficult to accomplish. The method is extensively analyzed in Golub and Varga (1961) and Varga (1962, Chap. 5).

# Problems

P11.2.1 Show that the Jacobi iteration converges for 2-by-2 symmetric positive definite systems.

P11.2.2 Show that if $A = M - N$ is singular, then we can never have $\rho ( M ^ { - 1 } N ) < 1$ even if M is nonsingular.

P11.2.3 (Supplied by R.S. Varga) Suppose that

$$
A _ {1} = \left[ \begin{array}{c c} 1 & - 1 / 2 \\ - 1 / 2 & 1 \end{array} \right], \qquad A _ {2} = \left[ \begin{array}{c c} 1 & - 3 / 4 \\ - 1 / 1 2 & 1 \end{array} \right].
$$

Let $J _ { 1 }$ and $J _ { 2 }$ be the associated Jacobi iteration matrices. Show that $\rho ( J _ { 1 } ) > \rho ( J _ { 2 } )$ , thereby refuting the claim that greater diagonal dominance implies more rapid Jacobi convergence.

P11.2.4 Suppose $A = T _ { n _ { 1 } } \otimes I _ { n _ { 2 } } \otimes I _ { n _ { 3 } } + I _ { n _ { 1 } } \otimes T _ { n _ { 2 } } \otimes I _ { n _ { 3 } } + I _ { n _ { 1 } } \otimes I _ { n _ { 2 } } \otimes T _ { n _ { 3 } }$ . If Jacobi’s method is

applied to the problem $A u = b ,$ then what is the spectral radius of the associated iteration matrix?

P11.2.5 A 5-point “stencil” is associated with the matrix $A = I _ { n _ { 1 } } \otimes T _ { n _ { 2 } } + T _ { n _ { 1 } } \otimes I _ { n _ { 2 } }$ and leads to the requirement that $U ( i , j )$ be the average of $U ( i - 1 , j ) , U ( i , j + 1 ) , \bar { U } ( i + 1 , \bar { j } )$ , and $U ( i , j - 1 )$ . Formulate a 9-point stencil procedure in which $U ( i , j )$ is a suitable average of its eight neighbors. (a) Describe the resulting matrix using Kronecker products. (b) If Jacobi’s method is used to solve $A u = b$ , then what is the spectral radius of the associated iteration matrix?

P11.2.6 Consider the linear system $( I _ { n _ { 1 } } \otimes { \mathcal { T } } _ { n _ { 2 } } + { \mathcal { T } } _ { n _ { 1 } } \otimes I _ { n _ { 2 } } ) x = b ,$ . What is the spectral radius of the iteration matrix for the block Jacobi iteration if the diagonal blocks are $n _ { 2 } { \mathrm { - b y } } { \mathrm { - } } n _ { 2 } { \mathrm { ? } }$

P11.2.7 Prove (11.2.13) and (11.2.14).

P11.2.8 Prove (11.2.15).

P11.2.9 Prove (11.2.17).

P11.2.10 Prove (11.2.24) and (11.2.25).

P11.2.11 Consider the 2-by-2 matrix

$$
A = \left[ \begin{array}{c c} 1 & \rho \\ - \rho & 1 \end{array} \right].
$$

(a) Under what conditions do we have $\rho ( M _ { \mathrm { G S } } ^ { - 1 } N _ { \mathrm { G S } } ) ~ < ~ 1 ?$ (b) For what range of $\omega$ do we have $\rho ( M _ { \omega } ^ { - 1 } N _ { \omega } ) < 1 ?$ What value of ω minimizes $\rho ( M _ { \omega } ^ { - 1 } N _ { \omega } ) ?$ (c) Repeat (a) and (b) for the matrix

$$
A = \left[ \begin{array}{c c} I _ {n} & S \\ - S ^ {T} & I _ {n} \end{array} \right]
$$

where $S \in \mathbb { R } ^ { n \times n }$ . Hint: Use the SVD of S.

P11.2.12 We want to investigate the solution of $A u \ : = \ : f$ where $A \neq A ^ { T }$ . For a model problem, consider the finite difference approximation to

$$
- u ^ {\prime \prime} + \sigma u ^ {\prime} = 0, \qquad 0 <   x <   1,
$$

where $u ( 0 ) = 1 0 \ \mathrm { a n d } \ u ( 1 ) = 1 0 \exp ^ { \sigma }$ . This leads to the difference equation

$$
- u _ {i - 1} + 2 u _ {i} - u _ {i + 1} + R (u _ {i + 1} - u _ {i - 1}) = 0, \quad i = 1: n,
$$

where $R = \sigma h / 2 , u _ { 0 } = 1 0$ , and $u _ { n + 1 } = 1 0 e ^ { \sigma }$ . The number R should be less than 1. What is the spectral radius of $M ^ { - 1 } N$ where $M \stackrel { . } { = } ( A + A ^ { T } ) / 2$ and $N = ( A ^ { T } - A ) / 2 ?$

P11.2.13 Consider the iteration

$$
y ^ {(k + 1)} = \omega (B y ^ {(k)} + d - y ^ {(k - 1)}) + y ^ {(k - 1)}
$$

where B has Schur decomposition $Q ^ { T } B Q \ = \ \mathrm { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ with $\lambda _ { 1 } \geq \cdots \geq \lambda _ { n }$ . Assume that $x = B x + d . \mathrm { ( a ) }$ Derive an equation for $\begin{array} { r } { \dot { e ^ { ( k ) } } = y ^ { ( \bar { k } ) } - x . { \mathrm { ( b ) } } } \end{array}$ Assume $\begin{array} { r } { y ^ { ( 1 ) } = B y ^ { ( 0 ) } + d . } \end{array}$ Show that $e ^ { ( k ) } = p _ { k } ( B ) e ^ { ( 0 ) }$ where $p _ { k }$ is an even polynomial if k is even and an odd polynomial if k is odd. (c) Write $f ^ { ( k ) } = Q ^ { T } e ^ { ( k ) }$ . Derive a difference equation for $f _ { j } ^ { ( k ) }$ for $j = 1 { : } n$ . Try to specify the exact solution for general $f _ { j } ^ { ( 0 ) }$ and $f _ { j } ^ { ( 1 ) }$ . (d) Show how to determine an optimal ω.

P11.2.14 Suppose we want to solve the linear least squares problem min $A x - b \parallel _ { 2 }$ where $A \in \mathbb { R } ^ { m \times n }$ , rank $( A ) = r \leq n$ , and $b \in \mathbb { R } ^ { m }$ . Consider the iterative scheme

$$
M x _ {i + 1} = N x _ {i} + A ^ {T} b
$$

where $M = ( A ^ { T } A + \lambda W ) , N = \lambda W , \lambda > 0$ and $W \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite. (a) Show that $M ^ { - 1 } N$ is diagonalizable and that $\rho ( M ^ { - 1 } N ) < 1 { \mathrm { ~ i f ~ r a n k } } ( A ) = n . { \mathrm { ~ ( b ) } }$ Suppose $x _ { 0 } = 0$ and that $\Vert \ v \ \Vert _ { W } = \left( v ^ { T } W v \right) ^ { - 1 / 2 }$ v
−1/2 , the “W -norm.” Show that regardless of A’s rank, the iterates xi converge $x _ { i }$ to the minimum W -norm solution to the least squares problem. (c) Show that if rank $( A ) = n$ then $\Vert \ b { x } _ { L S } - \ b { x } _ { i + 1 } \Vert _ { W } \leq \Vert \ b { x } _ { L S } - \ b { x } _ { i } \Vert _ { W }$ . (d) Show how to implement the iteration give the QR factorization of

$$
M = \left[ \begin{array}{c} A \\ \sqrt {\lambda} F \end{array} \right]
$$

where $W = F F ^ { T }$ is the Cholesky factorization of W .

P11.2.15 (a) Suppose $T \in \mathbb { R } ^ { n \times n }$ is tridiagonal with the property that $t _ { i , i + 1 } t _ { i + 1 , i } > 0$ for $i = 1 { : } n - 1$ .

Show that there is a diagonal matrix $D \in \mathbb { R } ^ { n \times n }$ so that $S = D T D ^ { - 1 }$ is symmetric. (b) Consider the following linear system for unknowns $u _ { 1 } , \ldots , u _ { n }$ :

$$
- u _ {i - 1} + 2 u _ {i} - u _ {i + 1} + \frac {\sigma h}{2} (u _ {i + 1} - u _ {i}) = f _ {i}, \quad i = 1: n.
$$

Assume $u _ { 0 } \equiv \alpha , u _ { n + 1 } \equiv \beta , \sigma > 0$ , and $h > 0$ . Under what conditions can this tridiagonal system be symmetrized using (a)? (c) Give formulae for the eigenvalues of the Jacobi iteration matrix.

# Notes and References for §11.2

For detailed treatment of the material in this section, see Greenbaum (IMSL, Chap. 10) or any of the following volumes:

R.S. Varga (1962). Matrix Iterative Analysis, Prentice-Hall, Englewood Cliffs, NJ.

D.M. Young (1971). Iterative Solution of Large Linear Systems, Academic Press, New York.

L.A. Hageman and D.M. Young (1981). Applied Iterative Methods, Academic Press, New York.

W. Hackbusch (1994). Iterative Solution of Large Sparse Systems of Equations, Springer-Verlag, New York.

As we mentioned, Young (1971) has the most comprehensive treatment of the SOR method. The object of SOR theory is to guide the user in choosing the relaxation parameter ω. In this setting, the ordering of equations and unknowns is critical, see:

M.J.M. Bernal and J.H. Verner (1968). “On Generalizing of the Theory of Consistent Orderings for Successive Over-Relaxation Methods,” Numer. Math. 12, 215–222.

D.M. Young (1970). “Convergence Properties of the Symmetric and Unsymmetric Over-Relaxation Methods,” Math. Comput. 24, 793–807.

D.M. Young (1972). “Generalization of Property A and Consistent Ordering,” SIAM J. Numer. Anal. 9, 454–463.

R.A. Nicolaides (1974). “On a Geometrical Aspect of SOR and the Theory of Consistent Ordering for Positive Definite Matrices,” Numer. Math. 12, 99–104.

A. Ruhe (1974). “SOR Methods for the Eigenvalue Problem with Large Sparse Matrices,” Math. Comput. 28, 695–710.

L. Adams and H. Jordan (1986). “Is SOR Color-Blind?” SIAM J. Sci. Stat. Comput. 7, 490–506.

M. Eiermann and R.S. Varga (1993). “Is the Optimal ω Best for the SOR Iteration Method,” Lin. Alg. Applic. 182, 257–277.

H. Lu (1999). “Stair Matrices and Their Generalizations with Applications to Iterative Methods I: A Generalization of the Successive Overrelaxation Method,” SIAM J. Numer. Anal. 37, 1–17.

An analysis of the Chebyshev semi-iterative method appears in:

G.H. Golub and R.S. Varga (1961). “Chebyshev Semi-Iterative Methods, Successive Over-Relaxation Iterative Methods, and Second-Order Richardson Iterative Methods, Parts I and II,” Numer. Math. 3, 147–156, 157–168.

That work is premised on the assumption that the underlying iteration matrix has real eigenvalues. How to proceed when this is not the case is discussed in:

T.A. Manteuffel (1977). “The Tchebychev Iteration for Nonsymmetric Linear Systems,” Numer. Math. 28, 307–327.

M. Eiermann and W. Niethammer (1983). “On the Construction of Semi-iterative Methods,” SIAM J. Numer. Anal. 20, 1153–1160.

W. Niethammer and R.S. Varga (1983). “The Analysis of k-step Iterative Methods for Linear Systems from Summability Theory,” Numer. Math. 41, 177–206.

G.H. Golub and M. Overton (1988). “The Convergence of Inexact Chebyshev and Richardson Iterative Methods for Solving Linear Systems,” Numer. Math. 53, 571–594.

D. Calvetti, G.H. Golub, and L. Reichel (1994). “An Adaptive Chebyshev Iterative Method for Nonsymmetric Linear Systems Based on Modified Moments,” Numer. Math. 67, 21–40.

E. Giladi, G.H. Golub, and J.B. Keller (1998). “Inner and Outer Iterations for the Chebyshev Algorithm,” SIAM J. Numer. Anal. 35, 300–319.

Other methods for unsymmetric problems are discussed in:

M. Eiermann, W. Niethammer, and R.S. Varga (1992). “Acceleration of Relaxation Methods for Non-Hermitian Linear Systems,” SIAM J. Matrix Anal. Applic. 13, 979–991.   
H. Elman and G.H. Golub (1990). “Iterative Methods for Cyclically Reduced Non-Self-Adjoint Linear Systems I,” Math. Comput. 54, 671–700.   
H. Elman and G.H. Golub (1990). “Iterative Methods for Cyclically Reduced Non-Self-Adjoint Linear Systems II,” Math. Comput. 56, 215–242.   
R. Bramley and A. Sameh (1992). “Row Projection Methods for Large Nonsymmetric Linear Systems,” SIAM J. Sci. Statist. Comput. 13, 168–193.

Iterative methods for complex symmetric systems are detailed in:

O. Axelsson and A. Kucherov (2000). “Real Valued Iterative Methods for Solving Complex Symmetric Linear Systems,” Numer. Lin. Alg. 7, 197–218.

V.E. Howle and S.A. Vavasis (2005). “An Iterative Method for Solving Complex-Symmetric Systems Arising in Electrical Power Modeling,” SIAM J. Matrix Anal. Applic. 26, 1150–1178.

Iterative methods for singular systems are discussed in:

A. Dax (1990). “The Convergence of Linear Stationary Iterative Processes for Solving Singular Unstructured Systems of Linear Equations,” SIAM Review 32, 611–635.   
Z.-H. Cao (2001). “A Note on Properties of Splittings of Singular Symmetric Positive Semidefinite Matrices,” Numer. Math. 88, 603–606.

Papers that are concerned with parallel implementation include:

D.J. Evans (1984). “Parallel SOR Iterative Methods,” Parallel Comput. 1, 3–18.

N. Patel and H. Jordan (1984). “A Parallelized Point Rowwise Successive Over-Relaxation Method on a Multiprocessor,” Parallel Comput. 1, 207–222.

R.J. Plemmons (1986). “A Parallel Block Iterative Scheme Applied to Computations in Structural Analysis,” SIAM J. Alg. Disc. Meth. 7, 337–347.

C. Kamath and A. Sameh (1989). “A Projection Method for Solving Nonsymmetric Linear Systems on Multiprocessors,” Parallel Computing 9, 291–312.

P. Amodio and F. Mazzia (1995). “A Parallel Gauss-Seidel Method for Block Tridiagonal Linear Systems,” SIAM J. Sci. Comput. 16, 1451–1461.

We have seen that the condition κ(A) is an important issue when direct methods are applied to Ax = b. However, the condition of the system also has a bearing on iterative method performance, see:

M. Arioli and F. Romani (1985). “Relations Between Condition Numbers and the Convergence of the Jacobi Method for Real Positive Definite Matrices,” Numer. Math. 46, 31–42.

M. Arioli, I.S. Duff, and D. Ruiz (1992). “Stopping Criteria for Iterative Solvers,” SIAM J. Matrix Anal. Applic. 13, 138–144.

Finally, the effect of rounding errors on the methods of this section is treated in:

H. Wozniakowski (1978). “Roundoff-Error Analysis of Iterations for Large Linear Systems,” Numer. Math. 30, 301–314.

P.A. Knight (1993). “Error Analysis of Stationary Iteration and Associated Problems,” Ph.D. thesis, Department of Mathematics, University of Manchester, England.
