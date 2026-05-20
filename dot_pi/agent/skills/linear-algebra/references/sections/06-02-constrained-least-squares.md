# 6.2 Constrained Least Squares

In the least squares setting it is sometimes natural to minimize $\parallel A x - b \parallel _ { 2 }$ over a proper subset of $\mathbb { R } ^ { n }$ . For example, we may wish to predict b as best we can with Ax subject to the constraint that x is a unit vector. Or perhaps the solution defines a fitting function f (t) which is to have prescribed values at certain points. This can lead to an equality-constrained least squares problem. In this section we show how these problems can be solved using the QR factorization, the SVD, and the GSVD.

# 6.2.1 Least Squares Minimization Over a Sphere

Given $\boldsymbol { A } \in \mathbb { R } ^ { m \times n } , \boldsymbol { b } \in \mathbb { R } ^ { m }$ , and a positive $\alpha \in \mathbb { R }$ , we consider the problem

$$
\min _ {\| x \| _ {2} \leq \alpha} \| A x - b \| _ {2}. \tag {6.2.1}
$$

This is an example of the LSQI (least squares with quadratic inequality constraint) problem. This problem arises in nonlinear optimization and other application areas. As we are soon to observe, the LSQI problem is related to the ridge regression problem discussed in §6.1.4.

Suppose

$$
A = U \Sigma V ^ {T} = \sum_ {i = 1} ^ {r} \sigma_ {i} u _ {i} v _ {i} ^ {T} \tag {6.2.2}
$$

is the SVD of A which we assume to have rank r. If the unconstrained minimum norm solution

$$
x _ {L S} = \sum_ {i = 1} ^ {r} \frac {u _ {i} ^ {T} b}{\sigma_ {i}} v _ {i}
$$

satisfies $\Vert \textbf { } x _ { L S } \textbf { } \Vert _ { 2 } \leq \alpha$ , then it obviously solves (6.2.1). Otherwise,

$$
\left\| x _ {L S} \right\| _ {2} ^ {2} = \sum_ {i = 1} ^ {r} \left(\frac {u _ {i} ^ {T} b}{\sigma_ {i}}\right) ^ {2} > \alpha^ {2}, \tag {6.2.3}
$$

and it follows that the solution to (6.2.1) is on the boundary of the constraint sphere. Thus, we can approach this constrained optimization problem using the method of Lagrange multipliers. Define the parameterized objective function $\phi$ by

$$
\phi (x, \lambda) = \frac {1}{2} \| A x - b \| _ {2} ^ {2} + \frac {\lambda}{2} \left(\| x \| _ {2} ^ {2} - \alpha^ {2}\right)
$$

and equate its gradient to zero. This gives a shifted normal equation system:

$$
(A ^ {T} A + \lambda I) \cdot x (\lambda) = A ^ {T} b.
$$

The goal is to choose λ so that $\parallel x ( \lambda ) \parallel _ { 2 } = \alpha$ . Using the SVD (6.2.2), this leads to the problem of finding a zero of the function

$$
f (\lambda) = \| x (\lambda) \| _ {2} ^ {2} - \alpha^ {2} = \sum_ {k = 1} ^ {n} \left(\frac {\sigma_ {k} u _ {k} ^ {T} b}{\sigma_ {k} ^ {2} + \lambda}\right) ^ {2} - \alpha^ {2}.
$$

This is an example of a secular equation problem. From (6.2.3), $f ( 0 ) > 0$ . Since $f ^ { \prime } ( \lambda ) < 0$ for $\lambda \geq 0$ , it follows that f has a unique positive root $\lambda _ { + }$ . It can be shown that

$$
\rho (\lambda) = \| A x (\lambda) - b \| _ {2} ^ {2} = \| A x _ {L S} - b \| _ {2} ^ {2} + \sum_ {i = 1} ^ {r} \left(\frac {\lambda u _ {i} ^ {T} b}{\sigma_ {i} ^ {2} + \lambda}\right) ^ {2}. \tag {6.2.4}
$$

It follows that $x ( \lambda _ { + } )$ solves (6.2.1).

Algorithm 6.2.1 Given $A \in \mathbb { R } ^ { m \times n }$ with $m \ge n , b \in \mathbb { R } ^ { m }$ , and $\alpha > 0$ , the following algorithm computes a vector $\boldsymbol { x } \in \mathbb { R } ^ { n }$ such that $\parallel A x - b \parallel _ { 2 }$ is minimum subject to the constraint that $\Vert \ b { x } \Vert _ { 2 } \leq \alpha$ .

Compute the SVD $A = U \Sigma V ^ { T }$ , save $V = [  v _ { 1 } | \cdot \cdot \cdot | v _ { n } ]$ , form $\tilde { b } = U ^ { T } b$ , and determine $r = \mathsf { r a n k } ( A )$ .

$\mathrm { i f } \sum _ { i = 1 } ^ { r } \left( \frac { { \tilde { b } } _ { i } } { \sigma _ { i } } \right) ^ { 2 } > \alpha ^ { 2 }$ r > α2

Find λ+ > 0 such that $\sum _ { i = 1 } ^ { r } \left( \frac { \sigma _ { i } \tilde { b } _ { i } } { \sigma _ { i } ^ { 2 } + \lambda _ { + } } \right) ^ { 2 } = \alpha ^ { 2 }$ .

$$
x = \sum_ {i = 1} ^ {r} \left(\frac {\sigma_ {i} \tilde {b} _ {i}}{\sigma_ {i} ^ {2} + \lambda_ {+}}\right) v _ {i}
$$

else

$$
x = \sum_ {i = 1} ^ {r} \left(\frac {\tilde {b} _ {i}}{\sigma_ {i}}\right) v _ {i}
$$

end

The SVD is the dominant computation in this algorithm.

# 6.2.2 More General Quadratic Constraints

A more general version of (6.2.1) results if we minimize $\parallel A x - b \parallel _ { 2 }$ over an arbitrary hyperellipsoid:

$$
\text { minimize } \parallel A x - b \parallel_ {2} \quad \text { subject   to } \parallel B x - d \parallel_ {2} \leq \alpha . \tag {6.2.5}
$$

Here we are assuming that $A \in \mathbb R ^ { m _ { 1 } \times n _ { 1 } } , b \in \mathbb R ^ { m _ { 1 } } , B \in \mathbb R ^ { m _ { 2 } \times n _ { 1 } } , d \in \mathbb R ^ { m _ { 2 } }$ , and $\alpha \geq 0$ . Just as the SVD turns (6.2.1) into an equivalent diagonal problem, we can use the GSVD to transform (6.2.5) into a diagonal problem. In particular, if the GSVD of A and B is given by (6.1.22) and (6.2.23), then (6.2.5) is equivalent to

$$
\text { minimize } \parallel D _ {A} y - \tilde {b} \parallel_ {2} \quad \text { subject   to } \parallel D _ {B} y - \tilde {d} \parallel_ {2} \leq \alpha \tag {6.2.6}
$$

where

$$
\tilde {b} = U _ {1} ^ {T} b, \qquad \tilde {d} = U _ {2} ^ {T} d, \qquad y = X ^ {- 1} x.
$$

The simple form of the objective function and the constraint equation facilitate the analysis. For example, if rank $( B ) = m _ { 2 } < n _ { 1 }$ , then

$$
\| D _ {A} y - \tilde {b} \| _ {2} ^ {2} = \sum_ {i = 1} ^ {n _ {1}} \left(\alpha_ {i} y _ {i} - \tilde {b} _ {i}\right) ^ {2} + \sum_ {i = n _ {1} + 1} ^ {m _ {1}} \tilde {b} _ {i} ^ {2} \tag {6.2.7}
$$

and

$$
\| D _ {B} y - \tilde {d} \| _ {2} ^ {2} = \sum_ {i = 1} ^ {m _ {2}} \left(\beta_ {i} y _ {i} - \tilde {d} _ {i}\right) ^ {2} + \sum_ {i = m _ {2} + 1} ^ {n _ {1}} \tilde {d} _ {i} ^ {2} \leq \alpha^ {2}. \tag {6.2.8}
$$

A Lagrange multiplier argument can be used to determine the solution to this transformed problem (if it exists).

# 6.2.3 Least Squares With Equality Constraints

We consider next the constrained least squares problem

$$
\min _ {B x = d} \| A x - b \| _ {2} \tag {6.2.9}
$$

where $A \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } }$ with $m _ { 1 } \geq n _ { 1 } , B \in \mathbb { R } ^ { m _ { 2 } \times n _ { 1 } }$ with $m _ { 2 } < n _ { 1 } , b \in \mathbb { R } ^ { m _ { 1 } }$ , and $d \in \mathbb { R } ^ { m _ { 2 } }$ . We refer to this as the LSE problem (least squares with equality constraints). By setting $\alpha = 0$ in (6.2.5) we see that the LSE problem is a special case of the LSQI problem. However, it is simpler to approach the LSE problem directly rather than through Lagrange multipliers.

For clarity, we assume that both A and B have full rank. Let

$$
Q ^ {T} B ^ {T} = \left[ \begin{array}{c} R \\ 0 \end{array} \right] _ {n _ {1} - m _ {2}} ^ {n _ {1}}
$$

be the QR factorization of $B ^ { T }$ and set

$$
A Q = \left[ \begin{array}{c c} A _ {1} & A _ {2} \\ m _ {2} & n _ {1} - m _ {2} \end{array} \right], \qquad Q ^ {T} x = \left[ \begin{array}{c} y \\ z \end{array} \right] \begin{array}{c} m _ {2} \\ n _ {1} - m _ {2} \end{array} .
$$

It is clear that with these transformations (6.2.9) becomes

$$
\min _ {R ^ {T} y = d} \| A _ {1} y + A _ {2} z - b \| _ {2}.
$$

Thus, $y$ is determined from the constraint equation $R ^ { T } y ~ = ~ d$ and the vector z is obtained by solving the unconstrained LS problem

$$
\min _ {z \in \mathbb {R} ^ {n _ {1} - m _ {2}}} \| A _ {2} z - (b - A _ {1} y) \| _ {2}.
$$

Combining the above, we see that the following vector solves the LSE problem:

$$
x = Q \left[ \begin{array}{c} y \\ z \end{array} \right].
$$

Algorithm 6.2.2 Suppose $A \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } } , \ B \in \mathbb { R } ^ { m _ { 2 } \times n _ { 1 } } , \ b \in \mathbb { R } ^ { m _ { 1 } }$ , and $d \in \mathbb { R } ^ { m _ { 2 } }$ . If rank $( A ) \ = \ n _ { 1 }$ and rank $( B ) \ = \ m _ { 2 } \ < \ n _ { 1 }$ , then the following algorithm minimizes $\parallel A x - b \parallel _ { 2 }$ subject to the constraint $B x = d$ .

Compute the QR factorization $B ^ { T } = Q R$

Solve $R ( 1 { : } m _ { 2 } , 1 { : } m _ { 2 } ) ^ { T } { \cdot } y = d$ for y.

$$
A = A Q
$$

${ \mathrm { F i n d ~ } } z { \mathrm { ~ s o ~ } } \parallel A ( : , m _ { 2 } + 1 : n _ { 1 } ) z - ( b - A ( : , 1 { : } m _ { 2 } ) \cdot y ) \parallel _ { 2 } { \mathrm { ~ i s ~ m i n i m i z e d } } .$

$$
x = Q (:, 1: m _ {2}) \cdot y + Q (:, m _ {2} + 1: n _ {1}) \cdot z.
$$

Note that this approach to the LSE problem involves two QR factorizations and a matrix multiplication. If A and/or B are rank deficient, then it is possible to devise a similar solution procedure using the SVD instead of QR. Note that there may not be a solution if rank $( B ) < m _ { 2 }$ . Also, if nul $( A ) \cap { \mathsf { n u l l } } ( B ) \neq \{ 0 \}$ and $d \in \mathsf { r a n } ( B )$ , then the LSE solution is not unique.

# 6.2.4 LSE Solution Using the Augmented System

The LSE problem can also be approached through the method of Lagrange multipliers. Define the augmented objective function

$$
f (x, \lambda) = \frac {1}{2} \| A x - b \| _ {2} ^ {2} + \lambda^ {T} (d - B x), \quad \lambda \in \mathbb {R} ^ {m _ {2}},
$$

and set to zero its gradient with respect to x:

$$
A ^ {T} A x - A ^ {T} b - B ^ {T} \lambda = 0.
$$

Combining this with the equations $r = b - A x$ and $B x = d$ we obtain the symmetric indefinite linear system

$$
\left[ \begin{array}{c c c} 0 & A ^ {T} & B ^ {T} \\ A & I & 0 \\ B & 0 & 0 \end{array} \right] \left[ \begin{array}{l} x \\ r \\ \lambda \end{array} \right] = \left[ \begin{array}{l} 0 \\ b \\ d \end{array} \right]. \tag {6.2.10}
$$

This system is nonsingular if both A and B have full rank. The augmented system presents a solution framework for the sparse LSE problem.

# 6.2.5 LSE Solution Using the GSVD

Using the GSVD given by (6.1.22) and (6.1.23), we see that the LSE problem transforms to

$$
\min _ {D _ {B} y = \tilde {d}} \| D _ {A} y - \tilde {b} \| _ {2} \tag {6.2.11}
$$

where $\tilde { b } = U _ { 1 } ^ { T } b , \tilde { d } = U _ { 2 } ^ { T } d .$ and $y = X ^ { - 1 } x$ . It follows that if $\mathsf { n u l l } ( A ) \cap \mathsf { n u l l } ( B ) = \{ 0 \}$ and $X = [ x _ { 1 } \mid \cdots \mid x _ { n } ]$ , then

$$
x = \sum_ {i = 1} ^ {m _ {2}} \left(\frac {\tilde {d} _ {i}}{\beta_ {i}}\right) x _ {i} + \sum_ {i = m _ {2} + 1} ^ {n _ {1}} \left(\frac {\tilde {b} _ {i}}{\alpha_ {i}}\right) x _ {i} \tag {6.2.12}
$$

solves the LSE problem.

# 6.2.6 LSE Solution Using Weights

An interesting way to obtain an approximate LSE solution is to solve the unconstrained LS problem

$$
\min _ {x} \left\| \left[ \begin{array}{c} A \\ \sqrt {\lambda} B \end{array} \right] x - \left[ \begin{array}{c} b \\ \sqrt {\lambda} d \end{array} \right] \right\| _ {2} \tag {6.2.13}
$$

for large λ. (Compare with the Tychanov regularization problem (6.1.21).) Since

$$
\left\| \left[ \begin{array}{c} A \\ \sqrt {\lambda} B \end{array} \right] x \right. - \left[ \begin{array}{c} b \\ \sqrt {\lambda} d \end{array} \right] \left\| \right. _ {2} ^ {2} = \left\| A x - b \right\| _ {2} ^ {2} + \lambda \left\| B x - d \right\| ^ {2},
$$

we see that there is a penalty for discrepancies among the constraint equations. To quantify this, assume that both A and B have full rank and substitute the GSVD defined by (6.1.22) and (6.1.23) into the normal equation system

$$
(A ^ {T} A + \lambda B ^ {T} B) x = A ^ {T} b + \lambda B ^ {T} d.
$$

This shows that the solution $x ( \lambda )$ is given by $x ( \lambda ) = X y ( \lambda )$ where $y ( \lambda )$ solves

$$
(D _ {A} ^ {T} D _ {A} + \lambda D _ {B} ^ {T} D _ {B}) y = D _ {A} ^ {T} \tilde {b} + \lambda D _ {B} ^ {T} \tilde {d}
$$

with $\tilde { b } = U _ { 1 } ^ { T } b$ and $\tilde { d } = U _ { 2 } ^ { T } d .$ . It follows that

$$
x (\lambda) = \sum_ {i = 1} ^ {m _ {2}} \left(\frac {\alpha_ {i} \tilde {b} _ {i} + \lambda \beta_ {i} \tilde {d} _ {i}}{\alpha_ {i} ^ {2} + \lambda \beta_ {i} ^ {2}}\right) x _ {i} + \sum_ {i = m _ {2} + 1} ^ {n _ {1}} \left(\frac {\tilde {b} _ {i}}{\alpha_ {i}}\right) x _ {i}
$$

and so from (6.2.13) we have

$$
x (\lambda) - x = \sum_ {i = 1} ^ {p} \frac {\alpha_ {i}}{\beta_ {i}} \left(\frac {\beta_ {i} u _ {i} ^ {T} b - \alpha_ {i} v _ {i} ^ {T} d}{\alpha_ {i} ^ {2} + \lambda^ {2} \beta_ {i} ^ {2}}\right) x _ {i}. \tag {6.2.14}
$$

This shows that $x ( \lambda ) \to x { \mathrm { ~ a s ~ } } \lambda \to \infty$ . The appeal of this approach to the LSE problem is that it can be implemented with unconstrained LS problem software. However, for large values of λ numerical problems can arise and it is necessary to take precautions. See Powell and Reid (1968) and Van Loan (1982).

# Problems

P6.2.1 Is the solution to (6.2.1) always unique?

P6.2.2 Let $p _ { 0 } ( x ) , \ldots , p _ { n } ( x )$ be given polynomials and $( x _ { 0 } , y _ { 0 } ) , \dots , ( x _ { m } , y _ { m } )$ be a given set of coordinate pairs with $x _ { i } \in [ a , b ]$ . It is desired to find a polynomial $\begin{array} { r } { p ( x ) \ = \ \sum _ { k = 0 } ^ { n } \alpha _ { k } p _ { k } ( x ) } \end{array}$ such that

$$
\phi (\alpha) = \sum_ {i = 0} ^ {m} (p (x _ {i}) - y _ {i}) ^ {2}
$$

is minimized subject to the constraint that

$$
\int_ {a} ^ {b} [ p ^ {\prime \prime} (x) ] ^ {2} d x \approx h \sum_ {i = 0} ^ {N} \left(\frac {p (z _ {i - 1}) - 2 p (z _ {i}) + p (z _ {i + 1})}{h ^ {2}}\right) ^ {2} \leq \alpha^ {2}
$$

where $z _ { i } = a + i h$ and $b = a + N h$ . Show that this leads to an LSQI problem of the form (6.2.5) with $d = 0$ .

P6.2.3 Suppose $Y = [ y _ { 1 } \mid \cdot \cdot \cdot \mid y _ { k } ] \in \mathbb { R } ^ { m \times k }$ has the property that

$$
Y ^ {T} Y = \mathrm{diag} (d _ {1} ^ {2}, \ldots , d _ {k} ^ {2}), \qquad d _ {1} \geq d _ {2} \geq \dots \geq d _ {k} > 0.
$$

Show that if $Y = Q R$ is the QR factorization of Y , then R is diagonal with $| r _ { i i } | = d _ { i }$ .

P6.2.4 (a) Show that if $( A ^ { T } A + \lambda I ) x = A ^ { T } b , \lambda > 0 .$ , and $\parallel x \parallel _ { 2 } = \alpha ,$ then $z = ( A x - b ) / \lambda$ solves the dual equations $( A A ^ { T } + \lambda I ) z = - b$ with $\parallel A ^ { T } z \parallel _ { 2 } = \alpha . ( \mathrm { b } )$ Show that if $( A A ^ { T } + \lambda I ) z = - b ,$ $\parallel A ^ { T } z \parallel _ { 2 } = \alpha$ , then $x = - A ^ { T } z$ satisfies $( A ^ { T } A + \lambda I ) x = { \overset { - } { A } } A ^ { T } b , \| x \| _ { 2 } = \alpha .$ .

P6.2.5 Show how to compute y (if it exists) so that both (6.2.7) and (6.2.8) are satisfied.

P6.2.6 Develop an SVD version of Algorithm 6.2.2 that can handle the situation when A and/or B are rank deficient.

P6.2.7 Suppose

$$
A = \left[ \begin{array}{l} A _ {1} \\ A _ {2} \end{array} \right]
$$

where $A _ { 1 } \in \mathbb { R } ^ { n \times n }$ is nonsingular and $A _ { 2 } \in \mathbb { R } ^ { ( m - n ) \times n }$ . Show that

$$
\sigma_ {\min} (A) \geq \sqrt {1 + \sigma_ {\min} (A _ {2} A _ {1} ^ {- 1}) ^ {2}} \sigma_ {\min} (A _ {1}).
$$

P6.2.8 Suppose $p \geq m \geq n$ and that $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { m \times p }$ Show how to compute orthogonal $Q \in \mathbb { R } ^ { m \times m }$ and orthogonal $V \in \mathbb { R } ^ { n \times n }$ so that

$$
Q ^ {T} A = \left[ \begin{array}{c} R \\ 0 \end{array} \right], \qquad Q ^ {T} B V = \left[ 0 \mid S \right]
$$

where $R \in \mathbb { R } ^ { n \times n }$ and $S \in \mathbb { R } ^ { m \times m }$ are upper triangular.

P6.2.9 Suppose $r \in \mathbb { R } ^ { m } , y \in \mathbb { R } ^ { n }$ , and $\delta > 0$ . Show how to solve the problem

$$
\min _ {E \in \mathbf {R} ^ {m \times n}, \| E \| _ {F} \leq \delta} \| E y - r \| _ {2}
$$

Repeat with “min” replaced by “max.”

P6.2.10 Show how the constrained least squares problem

$$
\min _ {B x = d} \| A x - b \| _ {2} \quad A \in \mathbb {R} ^ {m \times n}, B \in \mathbb {R} ^ {p \times n}, \operatorname{rank} (B) = p
$$

can be reduced to an unconstrained least square problem by performing p steps of Gaussian elimination on the matrix

$$
\left[ \begin{array}{l} B \\ A \end{array} \right] = \left[ \begin{array}{l l} B _ {1} & B _ {2} \\ A _ {1} & A _ {2} \end{array} \right], \qquad B _ {1} \in \mathbb {R} ^ {p \times p},   \mathsf {r a n k} (B _ {1}) = p.
$$

Explain. Hint: The Schur complement is of interest.

# Notes and References for §6.2

The LSQI problem is discussed in:

G.E. Forsythe and G.H. Golub (1965). “On the Stationary Values of a Second-Degree Polynomial on the Unit Sphere,” SIAM J. App. Math. 14, 1050–1068.

L. Eld´en (1980). “Perturbation Theory for the Least Squares Problem with Linear Equality Constraints,” SIAM J. Numer. Anal. 17, 338–350.

W. Gander (1981). “Least Squares with a Quadratic Constraint,” Numer. Math. 36, 291–307.

L. Eld´en (1983). “A Weighted Pseudoinverse, Generalized Singular Values, and Constrained Least Squares Problems,” BIT 22 , 487–502.

G.W. Stewart (1984). “On the Asymptotic Behavior of Scaled Singular Value and QR Decompositions,” Math. Comput. 43, 483–490.   
G.H. Golub and U. von Matt (1991). “Quadratically Constrained Least Squares and Quadratic Problems,” Numer. Math. 59, 561–580.   
T.F. Chan, J.A. Olkin, and D. Cooley (1992). “Solving Quadratically Constrained Least Squares Using Black Box Solvers,” BIT 32, 481–495.   
Secular equation root-finding comes up in many numerical linear algebra settings. For an algorithmic overview, see:   
O.E. Livne and A. Brandt (2002). “N Roots of the Secular Equation in O(N) Operations,” SIAM J. Matrix Anal. Applic. 24, 439–453.

For a discussion of the augmented systems approach to least squares problems, see:

˚A. Bj¨orck (1992). “Pivoting and Stability in the Augmented System Method,” Proceedings of the 14th Dundee Conference, D.F. Griffiths and G.A. Watson (eds.), Longman Scientific and Technical, Essex, U.K.   
˚A. Bj¨orck and C.C. Paige (1994). “Solution of Augmented Linear Systems Using Orthogonal Factorizations,” BIT 34, 1–24.

References that are concerned with the method of weighting for the LSE problem include:

M.J.D. Powell and J.K. Reid (1968). “On Applying Householder’s Method to Linear Least Squares Problems,” Proc. IFIP Congress, pp. 122–26.

C. Van Loan (1985). “On the Method of Weighting for Equality Constrained Least Squares Problems,” SIAM J. Numer. Anal. 22, 851–864.

J.L. Barlow and S.L. Handy (1988). “The Direct Solution of Weighted and Equality Constrained Least-Squares Problems,” SIAM J. Sci. Stat. Comput. 9, 704–716.

J.L. Barlow, N.K. Nichols, and R.J. Plemmons (1988). “Iterative Methods for Equality Constrained Least Squares Problems,” SIAM J. Sci. Stat. Comput. 9, 892–906.

J.L. Barlow (1988). “Error Analysis and Implementation Aspects of Deferred Correction for Equality Constrained Least-Squares Problems,” SIAM J. Numer. Anal. 25, 1340–1358.

J.L. Barlow and U.B. Vemulapati (1992). “A Note on Deferred Correction for Equality Constrained Least Squares Problems,” SIAM J. Numer. Anal. 29, 249–256.

M. Gulliksson and P.-˚A. Wedin (1992). “Modifying the QR-Decomposition to Constrained and Weighted Linear Least Squares,” SIAM J. Matrix Anal. Applic. 13, 1298–1313.

M. Gulliksson (1994). “Iterative Refinement for Constrained and Weighted Linear Least Squares,” BIT 34, 239–253.

G. W. Stewart (1997). “On the Weighting Method for Least Squares Problems with Linear Equality Constraints,” BIT 37, 961–967.

For the analysis of the LSE problem and related methods, see:

M. Wei (1992). “Perturbation Theory for the Rank-Deficient Equality Constrained Least Squares Problem,” SIAM J. Numer. Anal. 29, 1462–1481.

M. Wei (1992). “Algebraic Properties of the Rank-Deficient Equality-Constrained and Weighted Least Squares Problems,” Lin. Alg. Applic. 161, 27–44.

M. Gulliksson (1995). “Backward Error Analysis for the Constrained and Weighted Linear Least Squares Problem When Using the Weighted QR Factorization,” SIAM J. Matrix. Anal. Applic. 13, 675–687.

M. Gulliksson (1995). “Backward Error Analysis for the Constrained and Weighted Linear Least Squares Problem When Using the Weighted QR Factorization,” SIAM J. Matrix Anal. Applic. 16, 675–687.

J. Ding and W. Hang (1998). “New Perturbation Results for Equality-Constrained Least Squares Problems,” Lin. Alg. Applic. 272, 181–192.

A.J. Cox and N.J. Higham (1999). “Accuracy and Stability of the Null Space Method for Solving the Equality Constrained Least Squares Problem,” BIT 39, 34–50.

A.J. Cox and N.J. Higham (1999). “Row-Wise Backward Stable Elimination Methods for the Equality Constrained Least Squares Problem,” SIAM J. Matrix Anal. Applic. 21, 313–326.

A.J. Cox and Nicholas J. Higham (1999). “Backward Error Bounds for Constrained Least Squares Problems,” BIT 39, 210–227.

M. Gulliksson and P-A. Wedin (2000). “Perturbation Theory for Generalized and Constrained Linear Least Squares,” Num. Lin. Alg. 7, 181–195.   
M. Wei and A.R. De Pierro (2000). “Upper Perturbation Bounds of Weighted Projections, Weighted and Constrained Least Squares Problems,” SIAM J. Matrix Anal. Applic. 21, 931–951.   
E.Y. Bobrovnikova and S.A. Vavasis (2001). “Accurate Solution of Weighted Least Squares by Iterative Methods SIAM. J. Matrix Anal. Applic. 22, 1153–1174.   
M. Gulliksson, X-Q.Jin, and Y-M. Wei (2002). “Perturbation Bounds for Constrained and Weighted Least Squares Problems,” Lin. Alg. Applic. 349, 221–232.
