# 2.6 The Sensitivity of Square Systems

We use tools developed in previous sections to analyze the linear system problem $A x = b$ where $A \in \mathbb { R } ^ { n \times n }$ is nonsingular and $b \in \mathbb { R } ^ { n }$ . Our aim is to examine how perturbations in A and b affect the solution x. Higham (ASNA) offers a more detailed treatment.

# 2.6.1 An SVD Analysis

If

$$
A = \sum_ {i = 1} ^ {n} \sigma_ {i} u _ {i} v _ {i} ^ {T} = U \Sigma V ^ {T}
$$

is the SVD of A, then

$$
x = A ^ {- 1} b = (U \Sigma V ^ {T}) ^ {- 1} b = \sum_ {i = 1} ^ {n} \frac {u _ {i} ^ {T} b}{\sigma_ {i}} v _ {i}. \tag {2.6.1}
$$

This expansion shows that small changes in A or b can induce relatively large changes in x if $\sigma _ { n }$ is small.

It should come as no surprise that the magnitude of $\sigma _ { n }$ should have a bearing on the sensitivity of the $A x = b$ problem. Recall from Theorem 2.4.8 that $\sigma _ { n }$ is the 2-norm distance from A to the set of singular matrices. As the matrix of coefficients approaches this set, it is intuitively clear that the solution x should be increasingly sensitive to perturbations.

# 2.6.2 Condition

A precise measure of linear system sensitivity can be obtained by considering the parameterized system

$$
(A + \epsilon F) x (\epsilon) = b + \epsilon f, \qquad x (0) = x,
$$

where $F \in \mathbb { R } ^ { n \times n }$ and $f \in \mathbb { R } ^ { n }$ . If A is nonsingular, then it is clear that $x ( \epsilon )$ is differentiable in a neighborhood of zero. Moreover, $\dot { x } ( 0 ) \ = \ A ^ { - 1 } ( f - F x )$ and so the Taylor series expansion for $x ( \epsilon )$ has the form

$$
x (\epsilon) = x + \epsilon \dot {x} (0) + O (\epsilon^ {2}).
$$

Using any vector norm and consistent matrix norm we obtain

$$
\frac {\| x (\epsilon) - x \|}{\| x \|} \leq | \epsilon | \| A ^ {- 1} \| \left\{\frac {\| f \|}{\| x \|} + \| F \| \right\} + O \left(\epsilon^ {2}\right). \tag {2.6.2}
$$

For square matrices A define the condition number $\kappa ( A )$ by

$$
\kappa (A) = \left\| A \right\| \left\| A ^ {- 1} \right\| \tag {2.6.3}
$$

with the convention that $\kappa ( A ) = \infty$ for singular A. From $\left\| ~ b ~ \right\| ~ \leq ~ \left\| ~ A ~ \right\| ~ \left\| ~ x ~ \right\|$ and (2.6.2) it follows that

$$
\frac {\| x (\epsilon) - x \|}{\| x \|} \leq \kappa (A) (\rho_ {A} + \rho_ {b}) + O (\epsilon^ {2}) \tag {2.6.4}
$$

where

$$
\rho_ {A} = | \epsilon | \frac {\| F \|}{\| A \|} \text { and } \rho_ {b} = | \epsilon | \frac {\| f \|}{\| b \|}
$$

represent the relative errors in A and b, respectively. Thus, the relative error in x can be $\kappa ( A )$ times the relative error in A and b. In this sense, the condition number $\kappa ( A )$ quantifies the sensitivity of the $A x = b$ problem.

Note that $\kappa ( \cdot )$ depends on the underlying norm and subscripts are used accordingly, e.g.,

$$
\kappa_ {2} (A) = \| A \| _ {2} \| A ^ {- 1} \| _ {2} = \frac {\sigma_ {\max} (A)}{\sigma_ {\min} (A)}. \tag {2.6.5}
$$

Thus, the 2-norm condition of a matrix A measures the elongation of the hyperellipsoid $\{ A x : \| x \| _ { 2 } = 1 \}$ .

We mention two other characterizations of the condition number. For p-norm condition numbers, we have

$$
\frac {1}{\kappa_ {p} (A)} = \min _ {A + \Delta A \text {   singular }} \frac {\| \Delta A \| _ {p}}{\| A \| _ {p}}. \tag {2.6.6}
$$

This result may be found in Kahan (1966) and shows that $\kappa _ { p } ( A )$ measures the relative p-norm distance from A to the set of singular matrices.

For any norm, we also have

$$
\kappa (A) = \lim _ {\epsilon \rightarrow 0} \sup _ {\| \Delta A \| \leq \epsilon \| A \|} \frac {\| (A + \Delta A) ^ {- 1} - A ^ {- 1} \|}{\epsilon} \frac {1}{\| A ^ {- 1} \|}. \tag {2.6.7}
$$

This imposing result merely says that the condition number is a normalized Fr´echet derivative of the map $A \to A ^ { - 1 }$ . Further details may be found in Rice (1966). Recall that we were initially led to $\kappa ( A )$ through differentiation.

If $\kappa ( A )$ is large, then A is said to be an ill-conditioned matrix. Note that this is a norm-dependent property.1 However, any two condition numbers $\kappa _ { \alpha } ( \cdot )$ and $\kappa _ { \beta } ( \cdot )$ o n IRn×n $\mathbb { R } ^ { n \times n }$ are equivalent in that constants $c _ { 1 }$ and $c _ { 2 }$ can be found for which

$$
c _ {1} \kappa_ {\alpha} (A) \leq \kappa_ {\beta} (A) \leq c _ {2} \kappa_ {\alpha} (A), \quad A \in \mathbb {R} ^ {n \times n}.
$$

For example, on $\mathbb { R } ^ { n \times n }$ we have

$$
\frac {1}{n} \kappa_ {2} (A) \leq \kappa_ {1} (A) \leq n \kappa_ {2} (A),
$$

$$
\frac {1}{n} \kappa_ {\infty} (A) \leq \kappa_ {2} (A) \leq n \kappa_ {\infty} (A), \tag {2.6.8}
$$

$$
\frac {1}{n ^ {2}} \kappa_ {1} (A) \leq \kappa_ {\infty} (A) \leq n ^ {2} \kappa_ {1} (A).
$$

Thus, if a matrix is ill-conditioned in the α-norm, it is ill-conditioned in the $\beta \mathrm { { - n o r m } }$ modulo the constants $c _ { 1 }$ and $c _ { 2 }$ above.

For any of the p-norms, we have $\kappa _ { p } ( A ) \geq 1$ . Matrices with small condition numbers are said to be well-conditioned. In the 2-norm, orthogonal matrices are perfectly conditioned because if $Q$ is orthogonal, then $\kappa _ { 2 } ( Q ) = \| Q \| _ { 2 } \| Q ^ { T } \| _ { 2 } = 1$ .

# 2.6.3 Determinants and Nearness to Singularity

It is natural to consider how well determinant size measures ill-conditioning. If det $( A ) =$ 0 is equivalent to singularity, is det $( A ) \approx 0$ equivalent to near singularity? Unfortunately, there is little correlation between det(A) and the condition of $A x = b$ . For example, the matrix $B _ { n }$ defined by

$$
B _ {n} = \left[ \begin{array}{c c c c} 1 & - 1 & \dots & - 1 \\ 0 & 1 & \dots & - 1 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & 1 \end{array} \right] \in \mathbb {R} ^ {n \times n} \tag {2.6.9}
$$

has unit determinant, but $\kappa _ { \infty } ( B _ { n } ) ~ = ~ n \cdot 2 ^ { n - 1 }$ . On the other hand, a very wellconditioned matrix can have a very small determinant. For example,

$$
D _ {n} = \operatorname{diag} (1 0 ^ {- 1}, \dots , 1 0 ^ {- 1}) \in \mathbb {R} ^ {n \times n}
$$

satisfies $\kappa _ { p } ( D _ { n } ) = 1$ although det $( D _ { n } ) = 1 0 ^ { - n }$ .

# 2.6.4 A Rigorous Norm Bound

Recall that the derivation of (2.6.4) was valuable because it highlighted the connection between $\kappa ( A )$ and the rate of change of $x ( \epsilon )$ at $\epsilon = 0$ . However, it is a little unsatisfying because it is contingent on 
 being “small enough” and because it sheds no light on the size of the $O ( \epsilon ^ { 2 } )$ term. In this and the next subsection we develop some additional $A x = b$ perturbation theorems that are completely rigorous.

We first establish a lemma that indicates in terms of $\kappa ( A )$ when we can expect a perturbed system to be nonsingular.

Lemma 2.6.1. Suppose

$$
A x = b, \quad A \in \mathbb {R} ^ {n \times n}, 0 \neq b \in \mathbb {R} ^ {n},
$$

$$
(A + \Delta A) y = b + \Delta b, \quad \Delta A \in \mathbb {R} ^ {n \times n}, \Delta b \in \mathbb {R} ^ {n},
$$

with $\| \Delta A \| \le \epsilon \| A \|$ and $\| \Delta b \| \leq \epsilon \| b \|$ . $I f \epsilon \kappa ( A ) = r < 1$ , then $A + \Delta A$ is nonsingular and

$$
\frac {\parallel y \parallel}{\parallel x \parallel} \leq \frac {1 + r}{1 - r}.
$$

Proof. Since $\| A ^ { - 1 } \Delta A \| ~ \le ~ \epsilon \| ~ A ^ { - 1 } \| ~ \| ~ A \| = r < 1$ it follows from Theorem 2.3.4 that $( A + \Delta A )$ is nonsingular. Using Lemma 2.3.3 and the equality

$$
(I + A ^ {- 1} \Delta A) y = x + A ^ {- 1} \Delta b
$$

we find

$$
\begin{array}{l} \| y \| \leq \| (I + A ^ {- 1} \Delta A) ^ {- 1} \| \left(\| x \| + \epsilon \| A ^ {- 1} \| \| b \|\right) \\ \leq \frac {1}{1 - r} \left(\| x \| + \epsilon \| A ^ {- 1} \| \| b \|\right) = \frac {1}{1 - r} \left(\| x \| + r \frac {\| b \|}{\| A \|}\right). \\ \end{array}
$$

Since $\parallel b \parallel = \parallel A x \parallel \leq \parallel A \parallel \parallel x \parallel$ it follows that

$$
\parallel y \parallel \leq \frac {1}{1 - r} \left(\parallel x \parallel + r \parallel x \parallel\right)
$$

and this establishes the required inequality.

We are now set to establish a rigorous $A x = b$ perturbation bound.

Theorem 2.6.2. If the conditions of Lemma 2.6.1 hold, then

$$
\frac {\| y - x \|}{\| x \|} \leq \frac {2 \epsilon}{1 - r} \kappa (A). \tag {2.6.10}
$$

Proof. Since

$$
y - x = A ^ {- 1} \Delta b - A ^ {- 1} \Delta A y \tag {2.6.11}
$$

we have

$$
\| y - x \| \leq \epsilon \| A ^ {- 1} \| \| b \| + \epsilon \| A ^ {- 1} \| \| A \| \| y \|.
$$

Thus,

$$
\frac {\parallel y - x \parallel}{\parallel x \parallel} \leq \epsilon \kappa (A) \frac {\parallel b \parallel}{\parallel A \parallel \parallel x \parallel} + \epsilon \kappa (A) \frac {\parallel y \parallel}{\parallel x \parallel} \leq \epsilon \left(1 + \frac {1 + r}{1 - r}\right) \kappa (A),
$$

from which the theorem readily follows.

A small example helps put this result in perspective. The $A x = b$ problem

$$
{\left[ \begin{array}{c c} 1 & 0 \\ 0 & 1 0 ^ {- 6} \end{array} \right]} {\left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right]} = {\left[ \begin{array}{c} 1 \\ 1 0 ^ {- 6} \end{array} \right]}
$$

has solution $x = [ 1 , 1 ] ^ { T }$ and condition $\kappa _ { \infty } ( A ) = 1 0 ^ { 6 }$ . If $\Delta { b } = [ 1 0 ^ { - 6 } , 0 ] ^ { T } , \Delta { A } = 0$ , and $( A + \Delta A ) y = b + \Delta b$ , then $y = [ 1 + 1 0 ^ { - 6 } , 1 ] ^ { T }$ and the inequality (2.6.10) says

$$
1 0 ^ {- 6} = \frac {\parallel x - y \parallel_ {\infty}}{\parallel x \parallel_ {\infty}} \ll \frac {\parallel \Delta b \parallel_ {\infty}}{\parallel b \parallel_ {\infty}} \kappa_ {\infty} (A) = 1 0 ^ {- 6} 1 0 ^ {6} = 1.
$$

Thus, the upper bound in (2.6.10) can be a gross overestimate of the error induced by the perturbation.

On the other hand, if $\Delta b = ( 0 , 1 0 ^ { - 6 } ) ^ { T } , \Delta A = 0$ , and $( A + \Delta A ) y = b + \Delta b$ , then this inequality says that

$$
\frac {1 0 ^ {0}}{1 0 ^ {0}} \leq 2 \times 1 0 ^ {- 6} 1 0 ^ {6}.
$$

Thus, there are perturbations for which the bound in (2.6.10) is essentially attained.

# 2.6.5 More Refined Bounds

An interesting refinement of Theorem 2.6.2 results if we extend the notion of absolute value to matrices:

$$
F = (f _ {i j}) \in \mathbb {R} ^ {m \times n} \qquad \Rightarrow \qquad | F | = (| f _ {i j} |) \in \mathbb {R} ^ {m \times n}.
$$

This notation together with a matrix-level version of $\stackrel { 6 6 } { \leq } \underline { { \stackrel { \triangledown } { \scriptscriptstyle \mathrm { \blacktriangledown } } } }$ makes it easy to specify componentwise error bounds. If $F , G \in \mathbb { R } ^ { m \times n }$ , then

$$
| F | \leq | G | \quad \Leftrightarrow \quad | f _ {i j} | \leq | g _ {i j} |
$$

for all i and $j .$ Also note that if $F \in \mathbb { R } ^ { m \times q }$ and $G \in \mathbb { R } ^ { q \times n }$ , then $| F G | \leq | F | \cdot | G |$ . With these definitions and facts we obtain the following refinement of Theorem 2.6.2.

Theorem 2.6.3. Suppose

$$
A x = b, \quad A \in \mathbb {R} ^ {n \times n}, 0 \neq b \in \mathbb {R} ^ {n},
$$

$$
(A + \Delta A) y, = b + \Delta b \Delta A \in \mathbb {R} ^ {n \times n}, \Delta b \in \mathbb {R} ^ {n},
$$

and that $| \Delta A | \le \epsilon | A |$ and $| \Delta b | \le \epsilon | b | . \ I f \delta \kappa _ { \infty } ( A ) = r < 1$ , then $( A { + } \Delta A )$ is nonsingular and

$$
\frac {\| y - x \| _ {\infty}}{\| x \| _ {\infty}} \leq \frac {2 \epsilon}{1 - r} \cdot \| | A ^ {- 1} | | A | \| _ {\infty}. \tag {2.6.12}
$$

Proof. Since $\| \Delta A \| _ { \infty } \leq \epsilon \| A \| _ { \infty }$ and $\| \Delta b \| _ { \infty } \leq \epsilon \| b \| _ { \infty }$ the conditions of Lemma 2.6.1 are satisfied in the infinity norm. This implies that $A + \Delta A$ is nonsingular and

$$
\frac {\parallel y \parallel_ {\infty}}{\parallel x \parallel_ {\infty}} \leq \frac {1 + r}{1 - r}.
$$

Now using (2.6.11) we find

$$
\begin{array}{l} | y - x | \leq | A ^ {- 1} | | \Delta b | + | A ^ {- 1} | | \Delta A | | y | \\ \leq \epsilon | A ^ {- 1} | | b | + \epsilon | A ^ {- 1} | | A | | y | \leq \epsilon | A ^ {- 1} | | A | (| x | + | y |). \\ \end{array}
$$

If we take norms, then

$$
\| y - x \| _ {\infty} \leq \epsilon \| | A ^ {- 1} | | A | \| _ {\infty} \left(\| x \| _ {\infty} + \frac {1 + r}{1 - r} \| x \| _ {\infty}\right).
$$

The theorem follows upon division by $\| { \boldsymbol { x } } \| _ { \infty }$ .

The quantity $\parallel | A ^ { - 1 } | | A | \parallel _ { \infty }$ is known as the Skeel condition number and there are examples where it is considerably less than $\kappa _ { \infty } ( A )$ . In these situations, (2.6.12) is more informative than (2.6.10).

Norm bounds are frequently good enough when assessing error, but sometimes it is desirable to examine error at the component level. Oettli and Prager (1964) have an interesting result that indicates if an approximate solution $\hat { x } \in \mathbb { R } ^ { n }$ to the n-by-n system $A x = b$ satisfies a perturbed system with prescribed structure. Consider the problem of finding $\Delta A \in \mathbb { R } ^ { n \times n } , \Delta b \in \mathbb { R } ^ { n }$ , and $\omega \ge 0$ such that

$$
(A + \Delta A) \hat {x} = b + \Delta b \quad | \Delta A | \leq \omega | E |, | \Delta b | \leq \omega | f |. \tag {2.6.13}
$$

where $E \in \mathbb { R } ^ { n \times n }$ and $f \in \mathbb { R } ^ { n }$ are given. With proper choice of E and f, the perturbed system can take on certain qualities. For example, if $E = A$ and $f = b$ and ω is small, then ˆx satisfies a nearby system in the componentwise sense. The authors show that for a given $A , b , { \hat { x } } , E ,$ , and $f$ the smallest ω possible in (2.6.13) is given by

$$
\omega_ {\mathrm{min}} = \max _ {1 \leq i \leq n} \frac {| A \hat {x} - b | _ {i}}{(| E | \cdot | \hat {x} | + | f |) _ {i}}.
$$

If $A \hat { x } = b$ , then $\omega _ { \mathrm { m i n } } = 0$ . On the other hand, if $\omega _ { \mathrm { m i n } } = \infty$ , then ˆx does not satisfy any system of the prescribed perturbation structure.

# Problems

P2.6.1 Show that if $\parallel I \parallel \geq 1$ , then $\kappa ( A ) \geq 1$ .

P2.6.2 Show that for a given norm, $\kappa ( A B ) \leq \kappa ( A ) \kappa ( B )$ and that $\kappa ( \alpha A ) = \kappa ( A )$ for all nonzero α.

P2.6.3 Relate the 2-norm condition of $X \in \mathbb { R } ^ { m \times n } ( m \geq n )$ to the 2-norm condition of the matrices

$$
B = \left[ \begin{array}{c c} I _ {m} & X \\ 0 & I _ {n} \end{array} \right] \qquad \mathrm{and} \qquad C = \left[ \begin{array}{c} X \\ I _ {n} \end{array} \right].
$$

P2.6.4 Suppose $A \in \mathbb { R } ^ { n \times n }$ is nonsingular. Assume for a particular i and j that there is no way to make A singular by changing the value of $a _ { i j }$ . What can you conclude about $A ^ { - 1 \mathord { \ ? } }$ Hint: Use the Sherman-Morrison formula.

P2.6.5 Suppose $A \in \mathbb { R } ^ { n \times n }$ is nonsingular, $b \in \mathbb { R } ^ { n } , A x = b$ , and $C = A ^ { - 1 }$ . Use the Sherman-Morrison formula to show that

$$
\frac {\partial x _ {k}}{\partial a _ {i j}} = - x _ {j} c _ {k i}.
$$

# Notes and References for §2.6

The condition concept is thoroughly investigated in:

J. Rice (1966). “A Theory of Condition,” SIAM J. Numer. Anal. 3, 287–310.

W. Kahan (1966). “Numerical Linear Algebra,” Canadian Math. Bull. 9, 757–801.

References for componentwise perturbation theory include:

W. Oettli and W. Prager (1964). “Compatibility of Approximate Solutions of Linear Equations with Given Error Bounds for Coefficients and Right Hand Sides,” Numer. Math. 6, 405–409.

J.E. Cope and B.W. Rust (1979). “Bounds on Solutions of Systems with Accurate Data,” SIAM J. Numer. Anal. l6, 950–63.

R.D. Skeel (1979). “Scaling for Numerical Stability in Gaussian Elimination,” J. ACM 26, 494–526.

J.W. Demmel (1992). “The Componentwise Distance to the Nearest Singular Matrix,” SIAM J. Matrix Anal. Applic. 13, 10–19.

D.J. Higham and N.J. Higham (1992). “Componentwise Perturbation Theory for Linear Systems with Multiple Right-Hand Sides,” Lin. Alg. Applic. 174, 111–129.

N.J. Higham (1994). “A Survey of Componentwise Perturbation Theory in Numerical Linear Algebra,” in Mathematics of Computation 1943–1993: A Half Century of Computational Mathematics, W. Gautschi (ed.), Volume 48 of Proceedings of Symposia in Applied Mathematics, American Mathematical Society, Providence, RI.

S. Chandrasekaren and I.C.F. Ipsen (1995). “On the Sensitivity of Solution Components in Linear Systems of Equations,” SIAM J. Matrix Anal. Applic. 16, 93–112.

S.M. Rump (1999). “Ill-Conditioned Matrices Are Componentwise Near to Singularity,” SIAM Review 41, 102–112.

The reciprocal of the condition number measures how near a given Ax = b problem is to singularity. The importance of knowing how near is a given problem to a difficult or insoluble problem has come to be appreciated in many computational settings, see:

A. Laub(1985). “Numerical Linear Algebra Aspects of Control Design Computations,” IEEE Trans. Autom. Control. AC-30, 97–108.

J.W. Demmel (1987). “On the Distance to the Nearest Ill-Posed Problem,” Numer. Math. 51, 251–289.

N.J. Higham (1989). “Matrix Nearness Problems and Applications,” in Applications of Matrix Theory, M.J.C. Gover and S. Barnett (eds.), Oxford University Press, Oxford, UK, 1–27.

Much has been written about problem sensitivity from the statistical point of view, see:

J.W. Demmel (1988). “The Probability that a Numerical Analysis Problem is Difficult,” Math. Comput. 50, 449–480.

G.W. Stewart (1990). “Stochastic Perturbation Theory,” SIAM Review 32, 579–610.

C. S. Kenney, A.J. Laub, and M.S. Reese (1998). “Statistical Condition Estimation for Linear Systems,” SIAM J. Sci. Comput. 19, 566–583.

The problem of minimizing κ2(A + U V T ) where U V T is a low-rank matrix is discussed in:

C. Greif and J.M. Varah (2006). “Minimizing the Condition Number for Small Rank Modifications,” SIAM J. Matrix Anal. Applic. 29, 82–97.
