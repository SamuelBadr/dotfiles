# 8.5 Jacobi Methods

Jacobi methods for the symmetric eigenvalue problem attract current attention because they are inherently parallel. They work by performing a sequence of orthogonal similarity updates $A  Q ^ { T } A Q$ with the property that each new A, although full, is “more diagonal” than its predecessor. Eventually, the off-diagonal entries are small enough to be declared zero.

After surveying the basic ideas behind the Jacobi approach we develop a parallel Jacobi procedure.

# 8.5.1 The Jacobi Idea

The idea behind Jacobi’s method is to systematically reduce the quantity

$$
\operatorname{off}(A) = \sqrt{\sum_{i = 1}^{n}\sum_{\substack{j = 1\\ j\neq i}}^{n}a_{ij}^{2}}  ,
$$

i.e., the Frobenius norm of the off-diagonal elements. The tools for doing this are rotations of the form

$$
J (p, q, \theta) = \left[ \begin{array}{c c c c c c c} 1 & \dots & 0 & \dots & 0 & \dots & 0 \\ \vdots & \ddots & \vdots & & \vdots & & \vdots \\ 0 & \dots & c & \dots & s & \dots & 0 \\ \vdots & & \vdots & \ddots & \vdots & & \vdots \\ 0 & \dots & - s & \dots & c & \dots & 0 \\ \vdots & & \vdots & & \vdots & \ddots & \vdots \\ 0 & \dots & 0 & \dots & 0 & \dots & 1 \\ & & p & & q \end{array} \right] _ {q} ^ {p}
$$

which we call Jacobi rotations. Jacobi rotations are no different from Givens rotations; see §5.1.8. We submit to the name change in this section to honor the inventor.

The basic step in a Jacobi eigenvalue procedure involves (i) choosing an index pair $( p , q )$ that satisfies $1 \leq p < q \leq n$ , (ii) computing a cosine-sine pair $( c , s )$ such that

$$
\left[ \begin{array}{c c} b _ {p p} & b _ {p q} \\ b _ {q p} & b _ {q q} \end{array} \right] = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} \left[ \begin{array}{c c} a _ {p p} & a _ {p q} \\ a _ {q p} & a _ {q q} \end{array} \right] \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] \tag {8.5.1}
$$

is diagonal, and (iii) overwriting A with $B = J ^ { T } A J$ where $J = J ( p , q , \theta )$ . Observe that the matrix B agrees with A except in rows and columns p and q. Moreover, since the Frobenius norm is preserved by orthogonal transformations, we find that

$$
a _ {p p} ^ {2} + a _ {q q} ^ {2} + 2 a _ {p q} ^ {2} = b _ {p p} ^ {2} + b _ {q q} ^ {2} + 2 b _ {p q} ^ {2} = b _ {p p} ^ {2} + b _ {q q} ^ {2}.
$$

It follows that

$$
\begin{array}{l} \operatorname{off} (B) ^ {2} = \| B \| _ {F} ^ {2} - \sum_ {i = 1} ^ {n} b _ {i i} ^ {2} = \| A \| _ {F} ^ {2} - \sum_ {i = 1} ^ {n} a _ {i i} ^ {2} + \left(a _ {p p} ^ {2} + a _ {q q} ^ {2} - b _ {p p} ^ {2} - b _ {q q} ^ {2}\right) \tag {8.5.2} \\ = \mathrm{off} (A) ^ {2} - 2 a _ {p q} ^ {2}. \\ \end{array}
$$

It is in this sense that A moves closer to diagonal form with each Jacobi step.

Before we discuss how the index pair $( p , q )$ can be chosen, let us look at the actual computations associated with the $( p , q )$ subproblem.

# 8.5.2 The 2-by-2 Symmetric Schur Decomposition

To say that we diagonalize in (8.5.1) is to say that

$$
0 = b _ {p q} = a _ {p q} (c ^ {2} - s ^ {2}) + (a _ {p p} - a _ {q q}) c s. \tag {8.5.3}
$$

If $a _ { p q } = 0$ , then we just set c = 1 and s = 0. Otherwise, define

$$
\tau = \frac {a _ {q q} - a _ {p p}}{2 a _ {p q}} \mathrm{and} t = s / c
$$

and conclude from (8.5.3) that t = tan(θ) solves the quadratic

$$
t ^ {2} + 2 \tau t - 1 = 0.
$$

It turns out to be important to select the smaller of the two roots:

$$
t _ {\min} = \left\{ \begin{array}{l l} 1 / (\tau + \sqrt {1 + \tau^ {2}}) & \text {if} \tau \geq 0, \\ 1 / (\tau - \sqrt {1 + \tau^ {2}}) & \text {if} \tau <   0. \end{array} \right.
$$

This is implies that the rotation angle satisfies $| \theta | \leq \pi / 4$ and has the effect of maximizing c:

$$
c = 1 / \sqrt {1 + t _ {\mathrm{min}} ^ {2}}, \qquad s = t _ {\mathrm{min}} c.
$$

This in turn minimizes the difference between A and the update B:

$$
\| B - A\|_{F}^{2} = 4(1 - c)\sum_{\substack{i = 1\\ i\neq p,q}}^{n}(a_{ip}^{2} + a_{iq}^{2}) + 2a_{pq}^{2} / c^{2}.
$$

We summarize the 2-by-2 computations as follows:

Algorithm 8.5.1 Given an n-by-n symmetric A and integers p and q that satisfy $1 \leq p < q \leq n .$ , this algorithm computes a cosine-sine pair $\{ c , s \}$ such that if $B =$ $J ( p , q , \theta ) ^ { T } A J ( p , q , \theta )$ , then $b _ { p q } = b _ { q p } = 0$ .

function [c , s] = symSchur2(A, p, q)

if $A ( p , q ) \neq 0$

$$
\tau = (A (q, q) - A (p, p)) / (2 A (p, q))
$$

$\mathbf { i f } \ \tau \geq 0$

$$
t = 1 / (\tau + \sqrt {1 + \tau^ {2}})
$$

else

$$
t = 1 / (\tau - \sqrt {1 + \tau^ {2}})
$$

end

$$
c = 1 / \sqrt {1 + t ^ {2}}, s = t c
$$

else

$$
c = 1, s = 0
$$

end

# 8.5.3 The Classical Jacobi Algorithm

As we mentioned above, only rows and columns p and q are altered when the $( p , q )$ subproblem is solved. Once symSchur2 determines the 2-by-2 rotation, then the update $A  J ( p , q , \theta ) ^ { T } A J ( p , q , \theta )$ can be implemented in 6n flops if symmetry is exploited.

How do we choose the indices p and $q ?$ From the standpoint of maximizing the reduction of off(A) in (8.5.2), it makes sense to choose $( p , q )$ so that $a _ { p q } ^ { 2 }$ is maximal. This is the basis of the classical Jacobi algorithm.

Algorithm 8.5.2 (Classical Jacobi) Given a symmetric $A \in \mathbb { R } ^ { n \times n }$ and a positive tolerance tol, this algorithm overwrites A with $\overset { \cdot } { V } ^ { T } A V$ where V is orthogonal and off $( V ^ { T } A V ) \leq t o l \cdot \parallel A \parallel _ { F } ,$ .

$$
V = I _ {n}, \delta = \operatorname{tol} \cdot \| A \| _ {F}
$$

while off $( A ) > \delta$

$$
\text { Choose } (p, q) \text { so } | a _ {p q} | = \max _ {i \neq j} | a _ {i j} |
$$

$$
[ c, s ] = \operatorname{symSchur2} (A, p, q)
$$

$$
A = J (p, q, \theta) ^ {T} A J (p, q, \theta)
$$

$$
V = V J (p, q, \theta)
$$

end

Since $| a _ { p q } |$ is the largest off-diagonal entry,

$$
\mathsf {o f f} (A) ^ {2} \leq N (a _ {p q} ^ {2} + a _ {q p} ^ {2})
$$

where

$$
N = \frac {n (n - 1)}{2}.
$$

From (8.5.2) it follows that

$$
\operatorname{off} (B) ^ {2} \leq \left(1 - \frac {1}{N}\right) \operatorname{off} (A) ^ {2}.
$$

By induction, if $A ^ { ( k ) }$ denotes the matrix A after k Jacobi updates, then

$$
\operatorname{off} (A ^ {(k)}) ^ {2} \leq \left(1 - \frac {1}{N}\right) ^ {k} \operatorname{off} (A ^ {(0)}) ^ {2}.
$$

This implies that the classical Jacobi procedure converges at a linear rate.

However, the asymptotic convergence rate of the method is considerably better than linear. Schonhage (1964) and van Kempen (1966) show that for k large enough, there is a constant c such that

$$
\operatorname{off} (A ^ {(k + N)}) \leq c \cdot \operatorname{off} (A ^ {(k)}) ^ {2},
$$

i.e., quadratic convergence. An earlier paper by Henrici (1958) established the same result for the special case when A has distinct eigenvalues. In the convergence theory for the Jacobi iteration, it is critical that $| \theta | \leq \pi / 4$ . Among other things this precludes the possibility of interchanging nearly converged diagonal entries. This follows from the formulae $b _ { p p } = a _ { p p } - t a _ { p q }$ and $b _ { q q } = a _ { q q } + t a _ { p q }$ , which can be derived from Equation (8.5.1) and the definition t = sin(θ)/ cos(θ).

It is customary to refer to N Jacobi updates as a sweep. Thus, after a sufficient number of iterations, quadratic convergence is observed when examining off(A) after every sweep.

There is no rigorous theory that enables one to predict the number of sweeps that are required to achieve a specified reduction in off(A). However, Brent and Luk (1985) have argued heuristically that the number of sweeps is proportional to log(n) and this seems to be the case in practice.

# 8.5.4 The Cyclic-by-Row Algorithm

The trouble with the classical Jacobi method is that the updates involve $O ( n )$ flops while the search for the optimal $( p , q )$ is $O ( n ^ { 2 } )$ . One way to address this imbalance is to fix the sequence of subproblems to be solved in advance. A reasonable possibility is to step through all the subproblems in row-by-row fashion. For example, if $n = 4$ we cycle as follows:

$$
(p, q) = (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (1, 2), \dots .
$$

This ordering scheme is referred to as cyclic by row and it results in the following procedure:

Algorithm 8.5.3 (Cyclic Jacobi) Given a symmetric matrix $A \in \mathbb { R } ^ { n \times n }$ and a positive tolerance tol, this algorithm overwrites A with $V ^ { T } A V$ where V is orthogonal and off $( V ^ { T } A V ) \leq \mathsf { t o l } \cdot \parallel A \parallel _ { _ { F } }$ .

$V = I_{n},\quad \delta = \mathsf{tol}\cdot \| A\|_{F}$ while off $(A) > \delta$ for $p = 1:n - 1$ for $q = p + 1:n$ $[c,s] = \mathsf{symSchur2}(A,p,q)$ $A = J(p,q,\theta)^{T}AJ(p,q,\theta)$ $V = VJ(p,q,\theta)$ end   
end   
end

The cyclic Jacobi algorithm also converges quadratically. (See Wilkinson (1962) and van Kempen (1966).) However, since it does not require off-diagonal search, it is considerably faster than Jacobi’s original algorithm.

# 8.5.5 Error Analysis

Using Wilkinson’s error analysis it is possible to show that if r sweeps are required by Algorithm 8.5.3 and $d _ { 1 } , \ldots , d _ { n }$ specify the diagonal entries of the final, computed A

matrix, then

$$
\sum_ {i = 1} ^ {n} (d _ {i} - \lambda_ {i}) ^ {2} \leq (\operatorname{tol} + k _ {r} \mathbf {u}) \| A \| _ {F}
$$

for some ordering of A’s eigenvalues $\lambda _ { i }$ . The parameter $k _ { r }$ depends mildly on r.

Although the cyclic Jacobi method converges quadratically, it is not generally competitive with the symmetric QR algorithm. For example, if we just count flops, then two sweeps of Jacobi are roughly equivalent to a complete QR reduction to diagonal form with accumulation of transformations. However, for small n this liability is not very dramatic. Moreover, if an approximate eigenvector matrix V is known, then $V ^ { T } A V$ is almost diagonal, a situation that Jacobi can exploit but not QR.

Another interesting feature of the Jacobi method is that it can compute the eigenvalues with small relative error if A is positive definite. To appreciate this point, note that the Wilkinson analysis cited above coupled with the §8.1 perturbation theory ensures that the computed eigenvalues $\hat { \lambda } _ { 1 } \geq \cdots \hat { \geq } \hat { \lambda } _ { n }$ satisfy

$$
\frac {| \hat {\lambda} _ {i} - \lambda_ {i} (A) |}{\lambda_ {i} (A)} \approx \mathbf {u} \frac {\| A \| _ {2}}{\lambda_ {i} (A)} \leq \mathbf {u} \kappa_ {2} (A).
$$

However, a refined, componentwise error analysis by Demmel and Veseli´c (1992) shows that in the positive definite case

$$
\frac {| \hat {\lambda} _ {i} - \lambda_ {i} (A) |}{\lambda_ {i} (A)} \approx \mathbf {u} \kappa_ {2} (D ^ {- 1} A D ^ {- 1}) \tag {8.5.4}
$$

where $D = \operatorname { d i a g } ( { \sqrt { a _ { 1 1 } } } , \dots , { \sqrt { a _ { n n } } } )$ and this is generally a much smaller approximating bound. The key to establishing this result is some new perturbation theory and a demonstration that if $A _ { + }$ is a computed Jacobi update obtained from the current matrix $A _ { c } .$ , then the eigenvalues of $A _ { + }$ are relatively close to the eigenvalues of $A _ { c }$ in the sense of (8.5.4). To make the whole thing work in practice, the termination criterion is not based upon the comparison of off(A) with $\mathbf { u } \Vert { \cal A } \Vert _ { { \cal F } }$ but rather on the size of each $| a _ { i j } |$ compared to ${ \bf { u } } _ { \surd \overline { { a _ { i i } a _ { j j } } } }$ .

# 8.5.6 Block Jacobi Procedures

It is usually the case when solving the symmetric eigenvalue problem on a p-processor machine that $n \gg p .$ . In this case a block version of the Jacobi algorithm may be appropriate. Block versions of the above procedures are straightforward. Suppose that $n = r N$ and that we partition the n-by-n matrix A as follows:

$$
A = \left[ \begin{array}{c c c} A _ {1 1} & \dots & A _ {1 N} \\ \vdots & & \vdots \\ A _ {N 1} & \dots & A _ {N N} \end{array} \right].
$$

Here, each $A _ { i j }$ is $r { \mathrm { - } } \mathrm { b y } { \mathrm { - } } r$ . In a block Jacobi procedure the $( p , q )$ subproblem involves computing the 2r-by-2r Schur decomposition

$$
\left[ \begin{array}{c c} V _ {p p} & V _ {p q} \\ V _ {q p} & V _ {q q} \end{array} \right] ^ {T} \left[ \begin{array}{c c} A _ {p p} & A _ {p q} \\ A _ {q p} & A _ {q q} \end{array} \right] \left[ \begin{array}{c c} V _ {p p} & V _ {p q} \\ V _ {q p} & V _ {q q} \end{array} \right] = \left[ \begin{array}{c c} D _ {p p} & 0 \\ 0 & D _ {q q} \end{array} \right]
$$

and then applying to A the block Jacobi rotation made up of the $V _ { i j }$ . If we call this block rotation V , then it is easy to show that

$$
\mathsf {o f f} (V ^ {T} A V) ^ {2} = \mathsf {o f f} (A) ^ {2} - \left(2 \| A _ {p q} \| _ {F} ^ {2} + \mathsf {o f f} (A _ {p p}) ^ {2} + \mathsf {o f f} (A _ {q q}) ^ {2}\right).
$$

Block Jacobi procedures have many interesting computational aspects. For example, there are several ways to solve the subproblems, and the choice appears to be critical. See Bischof (1987).

# 8.5.7 A Note on the Parallel Ordering

The Block Jacobi approach to the symmetric eigenvalue problem has an inherent parallelism that has attracted significant attention. The key observation is that the $( i _ { 1 } , j _ { 1 } )$ subproblem is independent of the $( i _ { 2 } , j _ { 2 } )$ subproblem if the four indices $i _ { 1 } , j _ { 1 } , i _ { 2 }$ , and $j _ { 2 }$ are distinct. Moreover, if we regard the A as a 2m-by-2m block matrix, then it is possible to partition the set of off-diagonal index pairs into a collection of 2m − 1 rotation sets, each of which identifies m, nonconflicting subproblems.

A good way to visualize this is to imagine a chess tournament with 2m players in which everybody must play everybody else exactly once. Suppose m = 4. In “round 1” we have Player 1 versus Player 2, Player 3 versus Player 4, Player 5 versus Player 6, and Player 7 versus Player 8. Thus, there are four tables of action:

<table><tr><td>1</td><td>3</td><td>5</td><td>7</td></tr><tr><td>2</td><td>4</td><td>6</td><td>8</td></tr></table>

This corresponds to the first rotation set:

$$
\operatorname{rot.set} (1) = \{(1, 2), (3, 4), (5, 6), (7, 8) \}.
$$

To set up rounds 2 through 7, Player 1 stays put and Players 2 through 8 move from table to table in merry-go-round fashion:

<table><tr><td>1</td><td>2</td><td>3</td><td>5</td></tr><tr><td>4</td><td>6</td><td>8</td><td>7</td></tr></table>

$$
\operatorname{rot.set} (2) = \{(1, 4), (2, 6), (3, 8), (5, 7) \},
$$

<table><tr><td>1</td><td>4</td><td>2</td><td>3</td></tr><tr><td>6</td><td>8</td><td>7</td><td>5</td></tr></table>

$$
\operatorname{rot.set} (3) = \{(1, 6), (4, 8), (2, 7), (3, 5) \},
$$

<table><tr><td>1</td><td>6</td><td>4</td><td>2</td></tr><tr><td>8</td><td>7</td><td>5</td><td>3</td></tr></table>

$$
\operatorname{rot.set} (4) = \{(1, 8), (6, 7), (4, 5), (2, 3) \},
$$

<table><tr><td>1</td><td>8</td><td>6</td><td>4</td></tr><tr><td>7</td><td>5</td><td>3</td><td>2</td></tr></table>

$$
\operatorname{rot.set} (5) = \{(1, 7), (5, 8), (3, 6), (2, 4) \},
$$

<table><tr><td>1</td><td>7</td><td>8</td><td>6</td></tr><tr><td>5</td><td>3</td><td>2</td><td>4</td></tr></table>

$$
\operatorname{rot.set} (6) = \{(1, 5), (3, 7), (2, 8), (4, 6) \},
$$

<table><tr><td>1</td><td>5</td><td>7</td><td>8</td></tr><tr><td>3</td><td>2</td><td>4</td><td>6</td></tr></table>

$$
\operatorname{rot.set} (7) = \{(1, 3), (2, 5), (4, 7), (6, 8) \}.
$$

Taken in order, the seven rotation sets define the parallel ordering of the 28 possible off-diagonal index pairs.

For general m, a multiprocessor implementation would involve solving the subproblems within each rotation set in parallel. Although the generation of the subproblem rotations is independent, some synchronization is required to carry out the block similarity transform updates.

# Problems

P8.5.1 Let the scalar $\gamma$ be given along with the matrix

$$
A = \left[ \begin{array}{c c} w & x \\ x & z \end{array} \right].
$$

It is desired to compute an orthogonal matrix

$$
J = \left[ \begin{array}{c c} {c} & {s} \\ {- s} & {c} \end{array} \right]
$$

such that the (1, 1) entry of $J ^ { T } A J$ equals γ. Show that this requirement leads to the equation

$$
(w - \gamma) \tau^ {2} - 2 x \tau + (z - \gamma) = 0,
$$

where $\tau = c / s$ . Verify that this quadratic has real roots if $\gamma$ satisfies $\lambda _ { 2 } \leq \gamma \leq \lambda _ { 1 }$ , where $\lambda _ { 1 }$ and $\lambda _ { 2 }$ are the eigenvalues of A.

P8.5.2 Let $A \in \mathbb { R } ^ { n \times n }$ be symmetric. Give an algorithm that computes the factorization

$$
Q ^ {T} A Q = \gamma I + F
$$

where $Q$ is a product of Jacobi rotations, $\gamma = \operatorname { t r } ( A ) / n$ , and F has zero diagonal entries. Discuss the uniqueness of $Q$ .

P8.5.3 Formulate Jacobi procedures for (a) skew-symmetric matrices and (b) complex Hermitian matrices.

P8.5.4 Partition the n-by-n real symmetric matrix A as follows:

$$
A = \left[ \begin{array}{c c} a & v ^ {T} \\ v & A _ {1} \\ 1 & n - 1 \end{array} \right] _ {n - 1} ^ {1}.
$$

Let $Q$ be a Householder matrix such that if $B = Q ^ { T } A Q$ , then $B ( 3 { : } n , 1 ) = 0$ . Let $J = J ( 1 , 2 , \theta )$ b e determined such that if $C = J ^ { T } B J ,$ then $c _ { 1 2 } = 0$ and $c _ { 1 1 } \geq c _ { 2 2 }$ . Show $c _ { 1 1 } \geq a + \| v \| _ { 2 }$ . La Budde (1964) formulated an algorithm for the symmetric eigenvalue probem based upon repetition of this Householder-Jacobi computation.

P8.5.5 When implementing the cyclic Jacobi algorithm, it is sensible to skip the annihilation of $a _ { p q }$ if its modulus is less than some small, sweep-dependent parameter because the net reduction in off(A) is not worth the cost. This leads to what is called the threshold Jacobi method. Details concerning this variant of Jacobi’s algorithm may be found in Wilkinson (AEP, p. 277). Show that appropriate thresholding can guarantee convergence.

P8.5.6 Given a positive integer m, let $M = ( 2 m - 1 ) m$ . Develop an algorithm for computing integer vectors $i , j \in \mathbb { R } ^ { M }$ so that $( i _ { 1 } , j _ { 1 } ) , \dotsc , ( i _ { M } , j _ { M } )$ defines the parallel ordering.

# Notes and References for §8.5

Jacobi’s original paper is one of the earliest references found in the numerical analysis literature:

C.G.J. Jacobi (1846). “Uber ein Leichtes Verfahren Die in der Theorie der Sacularstroungen Vorkommendern Gleichungen Numerisch Aufzulosen,” Crelle’s J. 30, 51–94.

Prior to the QR algorithm, the Jacobi technique was the standard method for solving dense symmetric eigenvalue problems. Early references include:

M. Lotkin (1956). “Characteristic Values of Arbitrary Matrices,” Quart. Appl. Math. 14, 267–275.   
D.A. Pope and C. Tompkins (1957). “Maximizing Functions of Rotations: Experiments Concerning Speed of Diagonalization of Symmetric Matrices Using Jacobi’s Method,” J. ACM 4, 459–466.   
C.D. La Budde (1964). “Two Classes of Algorithms for Finding the Eigenvalues and Eigenvectors of Real Symmetric Matrices,” J. ACM 11, 53–58.   
H. Rutishauser (1966). “The Jacobi Method for Real Symmetric Matrices,” Numer. Math. 9, 1–10.

See also Wilkinson (AEP, p. 265) and:

J.H. Wilkinson (1968). “Almost Diagonal Matrices with Multiple or Close Eigenvalues,” Lin. Alg. Applic. 1, 1–12.

Papers that are concerned with quadratic convergence include:

P. Henrici (1958). “On the Speed of Convergence of Cyclic and Quasicyclic Jacobi Methods for Computing the Eigenvalues of Hermitian Matrices,” SIAM J. Appl. Math. 6, 144–162.

E.R. Hansen (1962). “On Quasicyclic Jacobi Methods,” J. ACM 9, 118–135.

J.H. Wilkinson (1962). “Note on the Quadratic Convergence of the Cyclic Jacobi Process,” Numer. Math. 6, 296–300.

E.R. Hansen (1963). “On Cyclic Jacobi Methods,” SIAM J. Appl. Math. 11, 448–459.

A. Schonhage (1964). “On the Quadratic Convergence of the Jacobi Process,” Numer. Math. 6, 410–412.

H.P.M. van Kempen (1966). “On Quadratic Convergence of the Special Cyclic Jacobi Method,” Numer. Math. 9, 19–22.

P. Henrici and K. Zimmermann (1968). “An Estimate for the Norms of Certain Cyclic Jacobi Operators,” Lin. Alg. Applic. 1, 489–501.

K.W. Brodlie and M.J.D. Powell (1975). “On the Convergence of Cyclic Jacobi Methods,” J. Inst. Math. Applic. 15, 279–287.

The ordering of the subproblems within a sweep is important:

W.F. Mascarenhas (1995). “On the Convergence of the Jacobi Method for Arbitrary Orderings,” SIAM J. Matrix Anal. Applic. 16, 1197–1209.

Z. Dramaˇc (1996). “On the Condition Behaviour in the Jacobi Method,” SIAM J. Matrix Anal. Applic. 17, 509–514.

V. Hari (2007). “Convergence of a Block-Oriented Quasi-Cyclic Jacobi Method,” SIAM J. Matrix Anal. Applic. 29, 349–369.

Z. Drmaˇc (2010). “A Global Convergence Proof for Cyclic Jacobi Methods with Block Rotations,” SIAM J. Matrix Anal. Applic. 31, 1329–1350.

Detailed error analyses that establish the high accuracy of Jacobi’s method include:

J. Barlow and J. Demmel (1990). “Computing Accurate Eigensystems of Scaled Diagonally Dominant Matrices,” SIAM J. Numer. Anal. 27, 762–791.

J.W. Demmel and K. Veseli´c (1992). “Jacobi’s Method is More Accurate than QR,” SIAM J. Matrix Anal. Applic. 13, 1204–1245.

W.F. Mascarenhas (1994). “A Note on Jacobi Being More Accurate than QR,” SIAM J. Matrix Anal. Applic. 15, 215–218.

R. Mathias (1995). “Accurate Eigensystem Computations by Jacobi Methods,” SIAM J. Matrix Anal. Applic. 16, 977–1003.

K. Veseli´c (1996). “A Note on the Accuracy of Symmetric Eigenreduction Algorithms,” ETNA 4, 37–45.

F.M. Dopico, J.M. Molera, and J. Moro (2003). “An Orthogonal High Relative Accuracy Algorithm for the Symmetric Eigenproblem,” SIAM J. Matrix Anal. Applic. 25, 301–351.

F.M. Dopico, P. Koev, and J.M. Molera (2008). “Implicit Standard Jacobi Gives High Relative Accuracy,” Numer. Math. 113, 519–553.

Attempts have been made to extend the Jacobi iteration to other classes of matrices and to push through corresponding convergence results. The case of normal matrices is discussed in:

H.H. Goldstine and L.P. Horowitz (1959). “A Procedure for the Diagonalization of Normal Matrices,” J. ACM 6, 176–195.

G. Loizou (1972). “On the Quadratic Convergence of the Jacobi Method for Normal Matrices,” Comput. J. 15, 274–276.

M.H.C. Paardekooper (1971). “An Eigenvalue Algorithm for Skew Symmetric Matrices,” Numer. Math. 17, 189–202.   
A. Ruhe (1972). “On the Quadratic Convergence of the Jacobi Method for Normal Matrices,” BIT 7, 305–313.   
D. Hacon (1993). “Jacobi’s Method for Skew-Symmetric Matrices,” SIAM J. Matrix Anal. Applic. 14, 619–628.   
Essentially, the analysis and algorithmic developments presented in the text carry over to the normal case with minor modification. For non-normal matrices, the situation is considerably more difficult:   
J. Greenstadt (1955). “A Method for Finding Roots of Arbitrary Matrices,” Math. Tables and Other Aids to Comp. 9, 47–52.   
C.E. Froberg (1965). “On Triangularization of Complex Matrices by Two Dimensional Unitary Tranformations,” BIT 5, 230–234.   
J. Boothroyd and P.J. Eberlein (1968). “Solution to the Eigenproblem by a Norm-Reducing Jacobi-Type Method (Handbook),” Numer. Math. 11, 1–12.   
A. Ruhe (1968). “On the Quadratic Convergence of a Generalization of the Jacobi Method to Arbitrary Matrices,” BIT 8, 210–231.   
A. Ruhe (1969). “The Norm of a Matrix After a Similarity Transformation,” BIT 9, 53–58.   
P.J. Eberlein (1970). “Solution to the Complex Eigenproblem by a Norm-Reducing Jacobi-type Method,” Numer. Math. 14, 232–245.   
C.P. Huang (1975). “A Jacobi-Type Method for Triangularizing an Arbitrary Matrix,” SIAM J. Numer. Anal. 12, 566–570.   
V. Hari (1982). “On the Global Convergence of the Eberlein Method for Real Matrices,” Numer. Math. 39, 361–370.   
G.W. Stewart (1985). “A Jacobi-Like Algorithm for Computing the Schur Decomposition of a Nonhermitian Matrix,” SIAM J. Sci. Stat. Comput. 6, 853–862.   
C. Mehl (2008). “On Asymptotic Convergence of Nonsymmetric Jacobi Algorithms,” SIAM J. Matrix Anal. Applic. 30, 291–311.

Jacobi methods for complex symmetric matrices have also been developed, see:

J.J. Seaton (1969). “Diagonalization of Complex Symmetric Matrices Using a Modified Jacobi Method,” Comput. J. 12, 156–157.

P.J. Eberlein (1971). “On the Diagonalization of Complex Symmetric Matrices,” J. Inst. Math. Applic. 7, 377–383.

P. Anderson and G. Loizou (1973). “On the Quadratic Convergence of an Algorithm Which Diagonalizes a Complex Symmetric Matrix,” J. Inst. Math. Applic. 12, 261–271.

P. Anderson and G. Loizou (1976). “A Jacobi-Type Method for Complex Symmetric Matrices (Handbook),” Numer. Math. 25, 347–363.

Other extensions include:

N. Mackey (1995). “Hamilton and Jacobi Meet Again: Quaternions and the Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 16, 421–435.

A.W. Bojanczyk (2003). “An Implicit Jacobi-like Method for Computing Generalized Hyperbolic SVD,” Lin. Alg. Applic. 358, 293–307.

For a sampling of papers concerned with various aspects of parallel Jacobi, see:

A. Sameh (1971). “On Jacobi and Jacobi-like Algorithms for a Parallel Computer,” Math. Comput. 25, 579–590.   
D.S. Scott, M.T. Heath, and R.C. Ward (1986). “Parallel Block Jacobi Eigenvalue Algorithms Using Systolic Arrays,” Lin. Alg. Applic. 77, 345–356.   
P.J. Eberlein (1987). “On Using the Jacobi Method on a Hypercube,” in Hypercube Multiprocessors, M.T. Heath (ed.), SIAM Publications, Philadelphia.   
G. Shroff and R. Schreiber (1989). “On the Convergence of the Cyclic Jacobi Method for Parallel Block Orderings,” SIAM J. Matrix Anal. Applic. 10, 326–346.   
M.H.C. Paardekooper (1991). “A Quadratically Convergent Parallel Jacobi Process for Diagonally Dominant Matrices with Nondistinct Eigenvalues,” Lin. Alg. Applic. 145, 71–88.   
T. Londre and N.H. Rhee (2005). “Numerical Stability of the Parallel Jacobi Method,” SIAM J. Matrix Anal. Applic. 26, 985–1000.
