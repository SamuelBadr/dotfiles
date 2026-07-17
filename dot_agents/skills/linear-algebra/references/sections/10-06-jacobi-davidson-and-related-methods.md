# 10.6 Jacobi-Davidson and Related Methods

We close the chapter with a brief discussion of the Jacobi-Davidson method, a solution framework that involves a mix of several important ideas. The starting point is a reformulation of the eigenvalue problem as a nonlinear systems problem, a maneuver that enables us to apply Newton-like methods. This leads in a natural way to a method of Jacobi that can be used to compute eigenvalue-eigenvector pairs of symmetric matrices that have a strong diagonal dominance. Eigenproblems of this variety arise in quantum chemistry and it is in that venue where Davidson (1975) developed a very successful generalization of the Jacobi procedure. It builds a (non-Krylov) nested sequence of subspaces and incorporates Ritz approximation. By restricting the Davidson corrections to the orthogonal complement of the current subspace, we arrive at the Jacobi-Davidson method developed by Sleijpen and van der Vorst (1996). Their technique does not require symmetry or diagonal dominance. Thus, in terms of abstraction, exposition in this section starts from the general, descends to the specific, and then climbs back out to the general. All along the way we are driven by practical, algorithmic concerns. Our presentation draws upon the insightful treatments of the Jacobi-Davidson method in Sorensen (2002) and Stewart (MAE, pp. 404–420).

We mention that full appreciation of the Jacobi-Davidson method and its versatility requires an understanding of the next chapter. This is because a critical step in the method requires the approximate solution of a large sparse linear system and preconditioned iterative solvers are typically brought into play. See §11.5.

# 10.6.1 The Approximate Newton Framework

Consider the n-by-n eigenvalue problem $A x = \lambda x$ and how we might improve an approximate eigenpair $\{ x _ { c } , \lambda _ { c } \}$ . Note that if

$$
A (x _ {c} + \delta x _ {c}) = (\lambda_ {c} + \delta \lambda_ {c}) (x _ {c} + \delta x _ {c}),
$$

then

$$
\left(A - \lambda_ {c} I\right) \delta x _ {c} - \delta \lambda_ {c} x _ {c} = - r _ {c} + \delta \lambda_ {c} \cdot \delta x _ {c}, \tag {10.6.1}
$$

where

$$
r _ {c} = A x _ {c} - \lambda_ {c} x _ {c}
$$

is the current residual. By ignoring the second-order term $\delta \lambda _ { c } \cdot \delta x _ { c }$ we arrive at the following specification for the corrections $\delta x _ { c }$ and $\delta \lambda _ { c }$ :

$$
\left(A - \lambda_ {c} I\right) \delta x _ {c} - \delta \lambda_ {c} x _ {c} = - r _ {c}. \tag {10.6.2}
$$

This is an underdetermined system of nonlinear equations that has a very uninteresting solution obtained by setting $\delta x _ { c } = - x _ { c }$ and $\delta \lambda _ { c } = 0$ . To keep away from this situation we add a constraint so that if

$$
\left[ \begin{array}{l} x _ {+} \\ \lambda_ {+} \end{array} \right] = \left[ \begin{array}{l} x _ {c} \\ \lambda_ {c} \end{array} \right] + \left[ \begin{array}{l} \delta x _ {c} \\ \delta \lambda_ {c} \end{array} \right], \tag {10.6.3}
$$

then the new eigenvector approximation $x _ { + }$ is nonzero. One way to do this is to require

$$
w ^ {T} x _ {+} = 1,
$$

where $w \in \mathbb { R } ^ { n }$ is an appropriately chosen nonzero vector. Possibilities include $w = x$ , which forces $x _ { + }$ to have unit 2-norm, and $w = e _ { 1 }$ , which forces its first component to be one. Regardless, if $x _ { c }$ is also normalized with respect to $w .$ , then

$$
w ^ {T} \delta x _ {c} = w ^ {T} (x _ {+} - x _ {c}) = 0. \tag {10.6.4}
$$

By assembling (10.6.2) and (10.6.4) into a single matrix-vector equation we obtain

$$
\left[ \begin{array}{c c} A - \lambda_ {c} I & - x _ {c} \\ w ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} \delta x _ {c} \\ \delta \lambda_ {c} \end{array} \right] = - \left[ \begin{array}{l} r _ {c} \\ 0 \end{array} \right]. \tag {10.6.5}
$$

This is precisely the Jacobian system that arises if Newton’s method is used to find a zero of the function

$$
F \left(\left[ \begin{array}{l} x \\ \lambda \end{array} \right]\right) = \left[ \begin{array}{l} A x - \lambda x \\ w ^ {T} x - 1 \end{array} \right].
$$

Its solution is easy to specify:

$$
\delta \lambda_ {c} = \frac {w ^ {T} (A - \lambda_ {c} I) ^ {- 1} r _ {c}}{w ^ {T} (A - \lambda_ {c} I) ^ {- 1} x _ {c}}, \tag {10.6.6}
$$

$$
\delta x _ {c} = - \left(A - \lambda_ {c} I\right) ^ {- 1} \left(r _ {c} - \delta \lambda_ {c} x _ {c}\right). \tag {10.6.7}
$$

Unfortunately, the required linear equation solving is problematic if A is large and sparse and this prompts us to consider the approximate Newton framework.

The idea behind approximate Newton methods is to replace the Jacobian system with a nearby, look-alike system that is easier to solve. One way to do this in our problem is to approximate A with a matrix M with the proviso that systems of the form $( M - \lambda _ { c } I ) z = r$ are “easy” to solve. If $N = M - A$ , then (10.6.5) transforms to

$$
\left[ \begin{array}{c c} M - \lambda_ {c} I & - x _ {c} \\ w ^ {T} & 0 \end{array} \right] \left[ \begin{array}{c} \delta x _ {c} \\ \delta \lambda_ {c} \end{array} \right] = - \left[ \begin{array}{c} r _ {c} - N \cdot \delta x _ {c} \\ 0 \end{array} \right].
$$

Continuing with the approximate-Newton mentality, let us throw away the inconvenient $N { \cdot } \delta x _ { c }$ term that is part of the right-hand side. This leaves us with the system

$$
\left[ \begin{array}{c c} M - \lambda_ {c} I & - x _ {c} \\ w ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} \delta x _ {c} \\ \delta \lambda_ {c} \end{array} \right] = - \left[ \begin{array}{l} r _ {c} \\ 0 \end{array} \right], \tag {10.6.8}
$$

and the following compute-friendly recipes for the corrections:

$$
\delta \lambda_ {c} = \frac {w ^ {T} (M - \lambda_ {c} I) ^ {- 1} r _ {c}}{w ^ {T} (M - \lambda_ {c} I) ^ {- 1} x _ {c}}, \tag {10.6.9}
$$

$$
\delta x _ {c} = - (M - \lambda_ {c} I) ^ {- 1} \left(r _ {c} - \delta \lambda_ {c} x _ {c}\right). \tag {10.6.10}
$$

Of course, by cutting corners in Newton’s method we risk losing quadratic convergence. Thus, the design of an approximate Newton strategy must balance the efficiency of the approximate Jacobian solution procedure with a possibly degraded rate of convergence. For an excellent discussion of this tension in the context of the eigenvalue problem, see Stewart (MAE, pp. 396–404).

# 10.6.2 The Jacobi Orthogonal Component Correction Method

Now suppose

$$
A = \left[ \begin{array}{c c} \alpha & c ^ {T} \\ c & A _ {1} \end{array} \right], \quad \alpha \in \mathbb {R}, c \in \mathbb {R} ^ {n - 1}, A _ {1} \in \mathbb {R} ^ {(n - 1) \times (n - 1)} \tag {10.6.11}
$$

is symmetric and strongly diagonally dominant. Assume that α is the largest element on the diagonal in absolute value. Our ambition is to compute λ (close to $\alpha )$ and $z \in \mathbb { R } ^ { n - 1 }$ so that

$$
\left[ \begin{array}{c c} \alpha & c ^ {T} \\ c & A _ {1} \end{array} \right] \left[ \begin{array}{l} 1 \\ z \end{array} \right] = \lambda \left[ \begin{array}{l} 1 \\ z \end{array} \right]. \tag {10.6.12}
$$

Because of the dominance assumption, there is no danger in assuming that the soughtafter eigenvector is nicely normalized by setting its first component to 1. Partition $\delta \boldsymbol { x } _ { c } .$ , $x _ { c } ,$ and $x _ { + }$ as follows:

$$
\delta x _ {c} = \left[ \begin{array}{l} \delta \mu_ {c} \\ \delta z _ {c} \end{array} \right], \qquad x _ {c} = \left[ \begin{array}{l} 1 \\ z _ {c} \end{array} \right], \qquad x _ {+} = \left[ \begin{array}{l} 1 \\ z _ {+} \end{array} \right].
$$

By substituting (10.6.11) and $w = e _ { 1 }$ into the Jacobian system (10.6.5), we get

$$
\left[ \begin{array}{c|cc} \alpha - \lambda_ {c} & c ^ {T} & - 1 \\ c & A _ {1} - \lambda_ {c} I & - z _ {c} \\ \hline 1 & 0 & 0 \\ \end{array} \right]\left[ \begin{array}{c} \delta \mu_ {c} \\ \hline \delta z _ {c} \\ \delta \lambda_ {c} \\ \end{array} \right] = - \left[ \begin{array}{c} \alpha + + c ^ {T} z _ {c} - \lambda_ {c} \\ (A _ {1} - \lambda_ {c} I) z _ {c} + c \\ \hline 0 \\ \end{array} \right],
$$

i.e.,

$$
\left[ \begin{array}{c c} A _ {1} - \lambda_ {c} I & - z _ {c} \\ c ^ {T} & - 1 \end{array} \right] \left[ \begin{array}{l} \delta z _ {c} \\ \delta \lambda_ {c} \end{array} \right] = - \left[ \begin{array}{c} (A _ {1} - \lambda_ {c} I) z _ {c} + c \\ \alpha + c ^ {T} z _ {c} - \lambda_ {c} \end{array} \right]. \tag {10.6.13}
$$

It is easy to verify that this is the Jacobian system that arises if Newton’s method is used to compute a zero of

$$
f \left(\left[ \begin{array}{l} z \\ \lambda \end{array} \right]\right) = \left[ \begin{array}{l l} \alpha & c ^ {T} \\ c & A _ {1} \end{array} \right] \left[ \begin{array}{l} 1 \\ z \end{array} \right] - \lambda \left[ \begin{array}{l} 1 \\ z \end{array} \right].
$$

If $A _ { 1 } = M _ { 1 } - N _ { 1 }$ , then (10.6.13) can be rearranged as follows:

$$
(M _ {1} - \lambda_ {c} I) z _ {+} = - c + N _ {1} z _ {c} + \{\delta \lambda_ {c} \cdot z _ {c} + N _ {1} \cdot \delta z _ {c} \},
$$

$$
\lambda_ {+} = \alpha + c ^ {T} z _ {+}.
$$

The Jacobi orthogonal component correction $( J O C C )$ method is defined by ignoring the terms enclosed by the curly brackets and taking $M _ { 1 }$ to be the diagonal part of $A _ { 1 } \mathrm { { : } }$ :

$$
\lambda_ {1} = \alpha , z _ {1} = 0 _ {n - 1}, \rho_ {1} = \left\| c \right\| _ {2}, k = 1
$$

while $\rho _ { k } >$ tol

$$
(M _ {1} - \lambda_ {k} I) z _ {k + 1} = - c + N _ {1} z _ {k}
$$

$$
\lambda_ {k + 1} = \alpha + c ^ {T} z _ {k + 1} \tag {10.6.14}
$$

$$
k = k + 1
$$

$$
\rho_ {k} = \left\| A _ {1} z _ {k} - \lambda_ {k} z _ {k} + c \right\| _ {2}
$$

end

The name of the method stems from the fact that the corrections to the approximate eigenvectors

$$
x _ {k} = \left[ \begin{array}{c} 1 \\ z _ {k} \end{array} \right],
$$

are all orthogonal to $e _ { 1 }$ . Indeed, it is clear from (10.6.14) that each residual

$$
r _ {k} = (A - \lambda_ {k} I) x _ {k}
$$

has a zero first component:

$$
r _ {k} = \left[ \begin{array}{c c} \alpha & c ^ {T} \\ c & A _ {1} \end{array} \right] \left[ \begin{array}{c} 1 \\ z _ {k} \end{array} \right] - \lambda_ {k} \left[ \begin{array}{c} 1 \\ z _ {k} \end{array} \right] = \left[ \begin{array}{c} 0 \\ (A _ {1} - \lambda_ {k} I) z _ {k} + c \end{array} \right]. \tag {10.6.15}
$$

Hence, the termination criterion in (10.6.14) is based on the size of the residual.

Jacobi intended this method to be use in conjunction with his diagonalization procedure for the symmetric eigenvalue problem. As discussed in §8.5, after a sufficient number of sweeps the matrix A is very close to being diagonal. At that point, the JOCC iteration (10.6.14) can be invoked after a possible $\breve { P } A P ^ { \breve { T } }$ update to maximize the (1,1) entry.

# 10.6.3 The Davidson Method

As with the JOCC iteration, Davidson’s method is applicable to the symmetric diagonally dominant eigenvalue problem (10.6.12). However, it involves a more sophisticated placement of the residual vectors. To motivate the main idea, let M be the diagonal part of A and use (10.6.15) to rewrite the JOCC iteration as follows:

$$
x _ {1} = e _ {1}, \lambda_ {1} = x _ {1} ^ {T} A x _ {1}, r _ {1} = A x _ {1} - \lambda_ {1} x _ {1}, V _ {1} = [ e _ {1} ], k = 1
$$

while $\parallel r _ { k } \parallel > \mathsf { t o l }$

Solve the residual correction equation:

$$
(M - \lambda_ {k} I) \delta v _ {k} = - r _ {k}.
$$

Compute an improved eigenpair $\{ \lambda _ { k + 1 } , x _ { k + 1 } \}$ so $r _ { k + 1 } \in \mathsf { r a n } ( V _ { 1 } ) ^ { \perp } \colon$

$$
\delta x _ {k} = \delta v _ {k}, x _ {k + 1} = x _ {k} + \delta x _ {k}, \lambda_ {k + 1} = \lambda_ {k} + c ^ {T} \delta x _ {k}
$$

$$
k = k + 1
$$

$$
r _ {k} = A x _ {k} - \lambda_ {k} x _ {k}
$$

end

Davidson’s method uses Ritz approximation to ensure that $r _ { k }$ is orthogonal to $e _ { 1 }$ and $\delta v _ { 1 } , \ldots , \delta v _ { k - 1 }$ . To acomplish this, the boxed fragment is replaced with the following:

Expand the current subspace ran $( V _ { k } )$ :

$$
s _ {k + 1} = (I - V _ {k} V _ {k} ^ {T}) \delta v _ {k}
$$

$$
v _ {k + 1} = s _ {k + 1} / \left\| s _ {k + 1} \right\| _ {2}, V _ {k + 1} = \left[ V _ {k} \mid v _ {k + 1} \right] \tag {10.6.16}
$$

Compute an improved eigenpair {λk+1, xk+1} so rk+1 ∈ ran(Vk+1)⊥:

$$
(V _ {k + 1} ^ {T} A V _ {k + 1}) t _ {k + 1} = \theta_ {k + 1} t _ {k + 1} \quad \text {(a suitably chosen Ritz pair)}
$$

$$
\lambda_ {k + 1} = \theta_ {k + 1}, x _ {k + 1} = V _ {k + 1} t _ {k + 1}
$$

There are a number of important issues associated with this method. To begin with, $V _ { k }$ is an n-by-k matrix with orthonormal columns. The transition from $V _ { k }$ to $V _ { k + 1 }$ can be effectively carried out by a modified Gram-Schmidt process. Of course, if k gets too big, then it may be necessary to restart the process using $v _ { k }$ as the initial vector.

Because $r _ { k } = A x _ { k } - \lambda _ { k } x _ { k } = A ( V _ { k } t _ { k } ) - \theta _ { k } ( V _ { k } t _ { k } )$ , it follows that

$$
V _ {k} ^ {T} r _ {k} = (V _ {k} ^ {T} A V _ {k}) t _ {k} - \theta_ {k} t _ {k} = 0,
$$

i.e., $r _ { k }$ is orthogonal to the range of $V _ { k }$ as required.

We mention that the Davidson algorithm can be generalized by allowing M to be a more involved approximation to A than just its diagonal part. See Crouzeix, Philippe, and Sadkane (1994) for details.

# 10.6.4 The Jacobi-Davidson Framework

Instead of forcing the correction $\delta x _ { c }$ to be orthogonal to $e _ { 1 }$ as in the Davidson setting, the Jacobi-Davidson method insists that $\delta x _ { c }$ be orthogonal to the current eigenvector approximation $x _ { c } .$ . The idea is to expand the current search space in a profitable, unexplored direction.

To see what is involved computationally and to connect with Newton’s method, we consider the following modification of (10.6.5):

$$
\left[ \begin{array}{c c} A - \lambda_ {c} I & - x _ {c} \\ x _ {c} ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} \delta x _ {c} \\ \delta \lambda_ {c} \end{array} \right] = - \left[ \begin{array}{l} r _ {c} \\ 0 \end{array} \right]. \tag {10.6.17}
$$

Note that this is the Jacobian system associated with the function

$$
F \left(\left[ \begin{array}{l} x \\ \lambda \end{array} \right]\right) = \left[ \begin{array}{c} A x - \lambda x \\ (x ^ {T} x - 1) / 2 \end{array} \right]
$$

given that $x _ { c } ^ { T } x _ { c } = 1$ . If $x _ { c }$ is so normalized and $\lambda _ { c } = x _ { c } ^ { T } A x _ { c } ,$ then from (10.6.17) we have

$$
\begin{array}{l} (I - x _ {c} x _ {c} ^ {T}) (A - \lambda_ {c} I) (I - x _ {c} x _ {c} ^ {T}) \delta x _ {c} = - (I - x _ {c} x _ {c} ^ {T}) (r _ {c} - \delta \lambda_ {c} x _ {c}) \\ = - \left(I - x _ {c} x _ {c} ^ {T}\right) r _ {c} \\ = - \left(I - x _ {c} x _ {c} ^ {T}\right) \left(A x _ {c} - \lambda_ {c} x _ {c}\right) \\ = - \left(I - x _ {c} x _ {c} ^ {T}\right) A x _ {c} \\ = - \left(A x _ {c} - \lambda_ {c} x _ {c}\right) = - r _ {c}. \\ \end{array}
$$

Thus, the correction $\delta \boldsymbol { x } _ { c }$ is obtained by solving the projected system

$$
(I - x _ {c} x _ {c} ^ {T}) (A - \lambda_ {c} I) (I - x _ {c} x _ {c} ^ {T}) \delta x _ {c} = - r _ {c} \tag {10.6.18}
$$

subject to the constraint that $x _ { c } ^ { T } \delta x _ { c } = 0$ .

In Jacobi-Davidson, approximate projected systems are used to expand the current subspace. Compared to the Davidson algorithm, everything remains the same in (10.6.16) except that instead of solving $( M - \lambda _ { c } I ) \delta v _ { k } = - r _ { k }$ to determine $\delta v _ { k }$ , we solve

$$
(I - x _ {k} x _ {k} ^ {T}) (M - \lambda_ {k} I) (I - x _ {k} x _ {k} ^ {T}) \delta v _ {k} = - r _ {k}, \tag {10.6.19}
$$

subject to the constraint that $x _ { k } ^ { T } \delta v _ { k } = 0$ . The resulting framework permits greater flexibility. The initial unit vector x1 can be arbitrary and various Chapter 11 iterative solvers can be applied to (10.6.19). See Sleijpen and van der Vorst (1996) and Sorensen (2002) for details.

The Jacobi-Davidson framework can be used to solve both symmetric and nonsymmetric eigenvalue problems and is important for the way it channels sparse $A x = b$ technology to the sparse $A x = \lambda x$ problem. It can be regarded as an approximate Newton iteration that is “steered” to the eigenpair of interest by Ritz calculations. Because an ever-expanding orthonormal basis is maintained, restarting has a key role to play as in the Arnoldi setting (§10.5).

# 10.6.5 The Trace-Min Algorithm

We briefly discuss the trace-min algorithm that can be used to compute the k smallest eigenvalues and associated eigenvectors for the n-by-n symmetric-definite problem $A x = \lambda B x$ . It has similarities to the Jacobi-Davidson procedure. The starting point is to realize that if $V _ { \mathrm { o p t } } \in \mathbb { R } ^ { n \times k }$ solves

$$
\min _ {V ^ {T} B V = I _ {k}} \operatorname{tr} (V ^ {T} A V),
$$

then the required eigenvalues/eigenvectors are exposed by $V _ { \mathrm { o p t } } ^ { T } A V _ { \mathrm { o p t } } = \mathrm { d i a g } ( \mu _ { 1 } , \dots , \mu _ { k } )$ （2号 and $A V _ { \mathrm { o p t } } ( : , j ) = \mu _ { j } B V _ { \mathrm { o p t } } ( : , j )$ , for $j = 1 { : } k$ . The method produces a sequence of $V -$ matrices, each of which satisfies $V ^ { T } B V = I _ { k }$ . The transition from $V _ { c }$ to $V _ { + }$ requires the solution of a projected system

$$
(I - Q _ {c} Q _ {c} ^ {T}) A (I - Q _ {c} Q _ {c} ^ {T}) Z _ {c} = A V _ {c}
$$

where $Z _ { c } \in \mathbb { R } ^ { n \times k }$ and $Q R = B V _ { c }$ is the thin QR factorization. This system, analogous to the central Jacobi-Davidson update system (10.6.19), can be solved using a suitably preconditioned conjugate gradient iteration. For details, see Sameh and Wisniewski (1982) and Sameh and Tong (2000).

# Problems

P10.6.1 How would you solve (10.6.1) assuming that A is upper Hessenberg?

P10.6.2 Assume that

$$
A = \left[ \begin{array}{c c} \alpha & b \\ b & D + E \end{array} \right]
$$

is an n-by-n symmetric matrix. Assume that D is the diagonal of A(2:n, 2:n) and that the eigenvalue gap $\delta = \lambda _ { 1 } ( A ) - \lambda _ { 2 } ( A )$ is positive. How small must b and E be in order to ensure that (D + E) − αI is diagonally dominant? Use Theorem 8.1.4.

# Notes and References for §10.6

For deeper perspectives on the methods of this section, we recommend Stewart (MAE, 404–420) and:

D.C. Sorensen (2002). “Numerical Methods for Large Eigenvalue Problems,” Acta Numerica 11, 519–584.

Davidson method papers include:

E.R. Davidson (1975). “The Iterative Calculation of a Few of the Lowest Eigenvalues and Corresponding Eigenvectors of Large Real Symmetric Matrices,” J. Comput. Phys. 17, 87–94.

R.B. Morgan and D.S. Scott (1986). “Generalizations of Davidson’s Method for Computing Eigenvalues of Sparse Symmetric Matrices,” SIAM J. Sci. Stat. Comput. 7, 817–825.

J. Olsen, P. Jorgensen, and J. Simons (1990). “Passing the One-Billion Limit in Full-Configuration (FCI) Interactions,” Chem. Phys. Letters 169, 463–472.

R.B. Morgan (1992). “Generalizations of Davidson’s Method for Computing Eigenvalues of Large Nonsymmetric Matrices,” J. Comput. Phys. 101, 287–291.

M. Sadkane (1993) “Block-Arnoldi and Davidson Methods for Unsymmetric Large Eigenvalue Problems,” Numer. Math. 64, 195–211.   
M. Crouzeix, B. Philippe, and M. Sadkane (1994). “The Davidson Method,” SIAM J. Sci. Comput. 15, 62–76.   
A. Strathopoulos, Y. Saad, and C.F. Fischer (1995). “Robust Preconditioning for Large, Sparse, Symmetric Eigenvalue Problems,” J. Comput. Appl. Math. 64, 197–215.   
The original Jacobi-Davidson idea appears in:   
G.L.G. Sleijpen and H.A. van der Vorst (1996). “A Jacobi-Davidson Iteration Method for Linear Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 17, 401–425.   
For applications and extensions to other problems, see:   
G.L.G. Sleijpen, A.G.L. Booten, D.R. Fokkema, and H.A. van der Vorst (1996). “Jacobi-Davidson Type Methods for Generalized Eigenproblems and Polynomial Eigenproblems,” BIT 36, 595–633.   
G.L.G. Sleijpen, H.A. van der Vorst, and E. Meijerink (1998). “Efficient Expansion of Subspaces in the Jacobi-Davidson Method for Standard and Generalized Eigenproblems,” ETNA 7, 75–89.   
D.R. Fokkema, G.L.G. Sleijpen, and H.A. van der Vorst (1998). “Jacobi-Davidson Style QR and QZ Algorithms for the Reduction of Matrix Pencils,” SIAM J. Sci. Computut. 20, 94–125.   
P. Arbenz and M.E. Hochstenbach (2004). “A Jacobi-Davidson Method for Solving Complex Symmetric Eigenvalue Problems,” SIAM J. Sci. Comput. 25, 1655–1673.

The trace-min method is detailed in:

A. Sameh and J. Wisniewski (1982). “A Trace Minimization Algorithm for the Generalized Eigenproblem,” SIAM J. Numer. Anal. 19, 1243–1259.   
A. Sameh and Z. Tong (2000). “A Trace Minimization Algorithm for the Symmetric Generalized Eigenproblem,” J. Comput. Appl. Math. 123, 155–175.
