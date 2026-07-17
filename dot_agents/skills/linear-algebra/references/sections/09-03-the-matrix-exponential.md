# 9.3 The Matrix Exponential

One of the most frequently computed matrix functions is the exponential

$$
e ^ {A t} = \sum_ {k = 0} ^ {\infty} \frac {(A t) ^ {k}}{k !}.
$$

Numerous algorithms for computing $e ^ { \boldsymbol { A } t }$ have been proposed, but most of them are of dubious numerical quality, as is pointed out in the survey articles by Moler and Van Loan (1978) and its update Moler and Van Loan (2003). In order to illustrate what the computational difficulties are, we present a “scaling and squaring” method based upon Pad´e approximation. A brief analysis of the method follows that involves some $e ^ { \boldsymbol { A } t }$ perturbation theory and includes comments about the shortcomings of eigenanalysis in settings where nonnormality prevails.

# 9.3.1 A Pad´e Approximation Method

Following the discussion in §9.2, if $g ( z ) \approx e ^ { z }$ , then $g ( A ) \approx e ^ { A }$ . A very useful class of approximants for this purpose are the Pad´e functions defined by

$$
R _ {p q} (z) = D _ {p q} (z) ^ {- 1} N _ {p q} (z),
$$

where

$$
N _ {p q} (z) = \sum_ {k = 0} ^ {p} \frac {(p + q - k) ! p !}{(p + q) ! k ! (p - k) !} z ^ {k}
$$

and

$$
D _ {p q} (z) = \sum_ {k = 0} ^ {q} \frac {(p + q - k) ! q !}{(p + q) ! k ! (q - k) !} (- z) ^ {k}.
$$

Notice that

$$
R _ {p o} (z) = 1 + z + \dots + z ^ {p} / p!
$$

is the order-p Taylor polynomial.

Unfortunately, the Pad´e approximants are good only near the origin, as the following identity reveals:

$$
e ^ {A} = R _ {p q} (A) + \frac {(- 1) ^ {q}}{(p + q) !} A ^ {p + q + 1} D _ {p q} (A) ^ {- 1} \int_ {0} ^ {1} u ^ {p} (1 - u) ^ {q} e ^ {A (1 - u)} d u. \tag {9.3.1}
$$

However, this problem can be overcome by exploiting the fact that

$$
e ^ {A} = (e ^ {A / m}) ^ {m}.
$$

In particular, we can scale A by m such that $F _ { p q } = R _ { p q } ( A / m )$ is a suitably accurate approximation to $e ^ { \boldsymbol { A } / m }$ . We then compute $F _ { p q } ^ { m }$ using Algorithm 9.2.2. If m is a power of two, then this amounts to repeated squaring and so is very efficient. The success of the overall procedure depends on the accuracy of the approximant

$$
F _ {p q} = \left(R _ {p q} \left(\frac {A}{2 ^ {j}}\right)\right) ^ {2 ^ {j}}.
$$

In Moler and Van Loan (1978) it is shown that, if

$$
\frac {\| A \| _ {\infty}}{2 ^ {j}} \leq \frac {1}{2},
$$

then there exists an $E \in \mathbb { R } ^ { n \times n }$ such that $F _ { p q } = e ^ { A + E } , A E = E A$ , and

$$
\| E \| _ {\infty} \leq \varepsilon (p, q) \| A \| _ {\infty},
$$

where

$$
\varepsilon (p, q) = 2 ^ {3 - (p + q)} \frac {p ! q !}{(p + q) ! (p + q + 1) !}.
$$

Using these results it is easy to establish the inequality

$$
\frac {\| e ^ {A} - F _ {p q} \| _ {\infty}}{\| e ^ {A} \| _ {\infty}} \leq \epsilon (p, q) \| A \| _ {\infty} e ^ {\epsilon (p, q) \| A \| _ {\infty}}.
$$

The parameters $p$ and $q$ can be determined according to some relative error tolerance. Since $F _ { p q }$ requires about $j + \operatorname* { m a x } \{ p , q \}$ matrix multiplications, it makes sense to set p $= q$ as this choice minimizes $\epsilon ( p , q )$ for a given amount of work. Overall we obtain

Algorithm 9.3.1 (Scaling and Squaring) Given $\delta > 0$ and $A \in \mathbb { R } ^ { n \times n }$ , the following algorithm computes $F = \stackrel { - } { e } ^ { A + E }$ where $\parallel E \parallel _ { \infty } \leq \delta \parallel A \parallel _ { \infty }$ .

$$
j = \max \{0, 1 + \operatorname{floor} (\log_ {2} (\| A \| _ {\infty})) \}
$$

$$
A = A / 2 ^ {j}
$$

Let $q$ be the smallest nonnegative integer such that $\epsilon ( q , q ) \leq \delta$

$$
D = I, N = I, X = I, c = 1
$$

for $k = 1 { : } q$

$$
c = c \cdot (q - k + 1) / ((2 q - k + 1) k)
$$

$$
X = A X, N = N + c \cdot X, D = D + (- 1) ^ {k} c \cdot X
$$

end

Solve $D F = N$ for $F$ using Gaussian elimination

for $k = 1 { : } j$

$$
F = F ^ {2}
$$

end

This algorithm requires about $2 ( q + j + 1 / 3 ) n ^ { 3 }$ flops. Its roundoff error properties of have been analyzed by Ward (1977). For further analysis and algorithmic improvements, see Higham (2005) and Al-Mohy and Higham (2009).

The special Horner techniques of §9.2.4 can be applied to quicken the computation of $D = D _ { q q } ( A )$ and $N = N _ { q q } ( A )$ . For example, if $q = 8$ we have $N _ { q q } ( A ) = U + A V$ and $D _ { q q } ( A ) = U - A V$ where

$$
U = c _ {0} I + c _ {2} A ^ {2} + (c _ {4} I + c _ {6} A ^ {2} + c _ {8} A ^ {4}) A ^ {4}
$$

and

$$
V = c _ {1} I + c _ {3} A ^ {2} + (c _ {5} I + c _ {7} A ^ {2}) A ^ {4}.
$$

Clearly, N and D can be computed with five matrix multiplications instead of seven as required by Algorithm 9.3.1.

# 9.3.2 Perturbation Theory

Is Algorithm 9.3.1 stable in the presence of roundoff error? To answer this question we need to understand the sensitivity of the matrix exponential to perturbations in A. The rich structure of this particular matrix function enables us to say more about the condition of the $e ^ { A }$ problem than is typically the case for a general matrix function. (See §9.1.6.)

The starting point in the discussion is the initial value problem

$$
\dot {X} (t) = A X (t), \qquad X (0) = I,
$$

where A, $\ b X ( t ) \in \mathbb R ^ { n \times n }$ . This has the unique solution $X ( t ) = e ^ { A t }$ , a characterization of the matrix exponential that can be used to establish the identity

$$
e ^ {(A + E) t} - e ^ {A t} = \int_ {0} ^ {t} e ^ {A (t - s)} E e ^ {(A + E) s} d s.
$$

From this it follows that

$$
\frac {\parallel e ^ {(A + E) t} - e ^ {A t} \parallel_ {2}}{\parallel e ^ {A t} \parallel_ {2}} \leq \frac {\parallel E \parallel_ {2}}{\parallel e ^ {A t} \parallel_ {2}} \int_ {0} ^ {t} \parallel e ^ {A (t - s)} \parallel_ {2} \parallel e ^ {(A + E) s} \parallel_ {2} d s.
$$

Further simplifications result if we bound the norms of the exponentials that appear in the integrand. One way of doing this is through the Schur decomposition. If $Q ^ { H } A Q \ =$ $\mathrm { d i a g } ( \lambda _ { i } ) + N$ is the Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ , then it can be shown that

$$
\| e ^ {A t} \| _ {2} \leq e ^ {\alpha (A) t M _ {S} (t)}, \tag {9.3.2}
$$

where

$$
\alpha (A) = \max \left\{\operatorname{Re} (\lambda): \lambda \in \lambda (A) \right\} \tag {9.3.3}
$$

is the spectral abscissa and

$$
M _ {S} (t) = \sum_ {k = 0} ^ {n - 1} \frac {\parallel N t \parallel_ {2} ^ {k}}{k !}.
$$

With a little manipulation it can be shown that

$$
\frac {\parallel e ^ {(A + E) t} - e ^ {A t} \parallel_ {2}}{\parallel e ^ {A t} \parallel_ {2}} \leq t \parallel E \parallel_ {2} M _ {S} (t) ^ {2} \exp (t M _ {S} (t) \parallel E \parallel_ {2}).
$$

Notice that $M _ { s } ( t ) \equiv 1$ if and only if A is normal, suggesting that the matrix exponential problem is “well-behaved” if A is normal. This observation is confirmed by the behavior of the matrix exponential condition number $\nu ( A , t )$ , defined by

$$
\nu (A, t) = \max _ {\| E \| \leq 1} \left\| \int_ {0} ^ {t} e ^ {A (t - s)} E e ^ {A s} d s \right\| _ {2} \frac {\| A \| _ {2}}{\| e ^ {A t} \| _ {2}}.
$$

This quantity, discussed by Van Loan (1977), measures the sensitivity of the map $A  e ^ { A t }$ in that for a given t, there is a matrix E for which

$$
\frac {\parallel e ^ {(A + E) t} - e ^ {A t} \parallel_ {2}}{\parallel e ^ {A t} \parallel_ {2}} \approx \nu (A, t) \frac {\parallel E \parallel_ {2}}{\parallel A \parallel_ {2}}.
$$

![](images/golub_550_599__13701d27df41c2cee17d499bfd734a0bca95b23cec76e74d9e672a10efe423c2.jpg)  
Figure 9.3.1. $\parallel e ^ { \boldsymbol { A } t } \parallel _ { 2 }$ can grow even if $\alpha ( A ) < 0$

Thus, if $\nu ( A , t )$ is large, small changes in A can induce relatively large changes in $e ^ { \boldsymbol { A } t }$ . Unfortunately, it is difficult to characterize precisely those A for which $\nu ( A , t )$ is large. (This is in contrast to the linear equation problem Ax = b, where the ill-conditioned A are neatly described in terms of SVD.) One thing we can say, however, is that $\nu ( A , t ) \geq t \| A \| _ { 2 }$ , with equality holding for all nonnegative t if and only if the matrix A is normal.

# 9.3.3 Pseudospectra

Dwelling a little more on the effect of nonnormality, we know from the analysis of §9.2 that approximating $e ^ { \boldsymbol { A } t }$ involves more than just approximating $e ^ { z t }$ on $\lambda ( A )$ . Another clue that eigenvalues do not “tell the whole story” in the $e ^ { \ b { A } t }$ problem has to do with the inability of the spectral abscissa (9.3.3) to predict the size of $\parallel e ^ { \boldsymbol { A } t } \parallel _ { 2 }$ as a function of time. If A is normal, then

$$
\| e ^ {A t} \| _ {2} = e ^ {\alpha (A) t}. \tag {9.3.4}
$$

Thus, there is uniform decay if the eigenvalues of A are in the open left half plane. But if A is non-normal, then $e ^ { \boldsymbol { A } \dot { \boldsymbol { t } } }$ can grow before decay sets in. The 2-by-2 example

$$
A = \left[ \begin{array}{c c} - 1 & 1 0 0 0 \\ 0 & - 1 \end{array} \right] \quad \Leftrightarrow \quad e ^ {A t} = e ^ {- t} \left[ \begin{array}{c c} 1 & 1 0 0 0 \cdot t \\ 0 & 1 \end{array} \right] \tag {9.3.5}
$$

plainly illustrates this point in Figure 9.3.1.

Pseudospectra can be used to shed light on the transient growth of $\parallel e ^ { A t } \parallel$ . For example, it can be shown that for every $\epsilon > 0$ ,

$$
\sup _ {t > 0} \| e ^ {A t} \| _ {2} \geq \frac {\alpha_ {\epsilon} (A)}{\epsilon} \tag {9.3.6}
$$

where $\alpha _ { \epsilon } ( A )$ is the 
-pseudospectral abscissa introduced in (7.8.8):

$$
\alpha_ {\epsilon} (A) = \sup _ {z \in \Lambda_ {\epsilon} (A)} \operatorname{Re} (z).
$$

For the 2-by-2 matrix in (9.3.5), it can be shown that $\alpha _ { . 0 1 } ( A ) / . 0 1 \approx 2 1 6$ , a value that is consistent with the growth curve in Figure 9.3.1. See Trefethen and Embree (SAP, Chap. 15) for more pseudospectral insights into the behavior of $\parallel e ^ { \boldsymbol { A } t } \parallel _ { 2 }$

# 9.3.4 Some Stability Issues

With this discussion we are ready to begin thinking about the stability of Algorithm 9.3.1. A potential difficulty arises during the squaring process if A is a matrix whose exponential grows before it decays. If

$$
G = R _ {q q} \left(\frac {A}{2 ^ {j}}\right) \approx e ^ {A / 2 ^ {j}},
$$

then it can be shown that rounding errors of order

$$
\gamma = \mathbf {u} \| G ^ {2} \| _ {2} \cdot \| G ^ {4} \| _ {2} \cdot \| G ^ {8} \| _ {2} \dots \| G ^ {2 ^ {j - 1}} \| _ {2}
$$

can be expected to contaminate the computed $G ^ { 2 ^ { j } } . \mathrm { ~ H ~ } \parallel \boldsymbol { e } ^ { A t } \parallel _ { 2 }$ has a substantial initial growth, then it may be the case that

$$
\gamma \gg \mathbf {u} \| G ^ {2 ^ {j}} \| _ {2} \approx \mathbf {u} \| e ^ {A} \| _ {2},
$$

thus ruling out the possibility of small relative errors.

If A is normal, then so is the matrix G and therefore $\parallel G ^ { m } \parallel _ { 2 } = \parallel G \parallel _ { 2 } ^ { m }$ for all positive integers m. Thus, $\gamma \approx \mathbf { u } \parallel G ^ { 2 ^ { j } } \parallel _ { 2 } \approx \mathbf { u } \parallel e ^ { A } \parallel _ { 2 }$ and so the initial growth problems disappear. The algorithm can essentially be guaranteed to produce small relative error when A is normal. On the other hand, it is more difficult to draw conclusions about the method when A is nonnormal because the connection between $\nu ( A , t )$ and the initial growth phenomena is unclear. However, numerical experiments suggest that Algorithm 9.3.1 fails to produce a relatively accurate $e ^ { A }$ only when v(A, 1) is correspondingly large.

# Problems

P9.3.1 Show that $e ^ { ( A + B ) t } = e ^ { A t } e ^ { B t }$ for all t if and only if AB = BA. Hint: Express both sides as a power series in t and compare the coefficient of t.

P9.3.2 Suppose that A is skew-symmetric. Show that both $e ^ { A }$ and the (1,1) Pad´e approximatant $R _ { 1 1 } ( A )$ are orthogonal. Are there any other values of p and q for which $R _ { p q } ( A )$ is orthogonal?

P9.3.3 Show that if A is nonsingular, then there exists a matrix X such that $A = e ^ { X }$ . Is X unique?

P9.3.4 Show that if

$$
\exp \left(\left[ \begin{array}{c c} - A ^ {T} & P \\ 0 & A \end{array} \right] z\right) = \left[ \begin{array}{c c} F _ {1 1} & F _ {1 2} \\ 0 & F _ {2 2} \end{array} \right] _ {n} ^ {n}
$$

then

$$
F _ {1 1} ^ {T} F _ {1 2} = \int_ {0} ^ {z} e ^ {A ^ {T} t} P e ^ {A t} d t.
$$

P9.3.5 Give an algorithm for computing $e ^ { A }$ when $A = u v ^ { T } , u , v \in \mathbb { R } ^ { n }$ .

P9.3.6 Suppose $A \in \mathbb { R } ^ { n \times n }$ and that $v \in \mathbb { R } ^ { n }$ has unit 2-norm. Define the function $\phi ( t ) = \parallel e ^ { A t } v \parallel _ { 2 } ^ { 2 } / 2$ and show that

$$
\dot {\phi} (t) \leq \mu (A) \phi (t)
$$

where $\mu ( A ) = \lambda _ { 1 } ( ( A + A ^ { T } ) / 2 )$ . Conclude that

$$
\parallel e ^ {A t} \parallel_ {2} \leq e ^ {\mu (A) t}
$$

where $t \geq 0 .$

P9.3.7 Suppose $A \in \mathbb { R } ^ { n \times n }$ has the property that its off-diagonal entries are negative and its column sums are zero. Show that for all t, $F = \exp ( A t )$ has nonnegative entries and unit column sums.

# Notes and References for §9.3

Much of what appears in this section and an extensive bibliography may be found in the following survey articles:

C.B. Moler and C.F. Van Loan (1978). “Nineteen Dubious Ways to Compute the Exponential of a Matrix,” SIAM Review 20, 801–836.

C.B. Moler and C.F.Van Loan (2003). “Nineteen Dubious Ways to Compute the Exponential of a Matrix, Twenty-Five Years Later,” SIAM Review 45, 3–49.

Scaling and squaring with Pad´e approximants (Algorithm 9.3.1) and a careful implementation of the Schur decomposition method (Algorithm 9.1.1) were found to be among the less dubious of the nineteen methods scrutinized. Various aspects of Pad´e approximation of the matrix exponential are discussed in:

W. Fair and Y. Luke (1970). “Pad´e Approximations to the Operator Exponential,” Numer. Math. 14, 379–382.

C.F. Van Loan (1977). “On the Limitation and Application of Pad´e Approximation to the Matrix Exponential,” in Pad´e and Rational Approximation, E.B. Saff and R.S. Varga (eds.), Academic Press, New York.

R.C. Ward (1977). “Numerical Computation of the Matrix Exponential with Accuracy Estimate,” SIAM J. Numer. Anal. 14, 600–614.

A. Wragg (1973). “Computation of the Exponential of a Matrix I: Theoretical Considerations,” J. Inst. Math. Applic. 11, 369–375.

A. Wragg (1975). “Computation of the Exponential of a Matrix II: Practical Considerations,” J. Inst. Math. Applic. 15, 273–278.

L. Dieci and A. Papini (2000). “Pad´e Approximation for the Exponential of a Block Triangular Matrix,” Lin. Alg. Applic. 308, 183–202.

M. Arioli, B. Codenotti and C. Fassino (1996). “The Pad´e Method for Computing the Matrix Exponential,” Lin. Alg. Applic. 240, 111–130.

N.J. Higham (2005). “The Scaling and Squaring Method for the Matrix Exponential Revisited,” SIAM J. Matrix Anal. Applic. 26, 1179–1193.

A.H. Al-Mohy and N.J. Higham (2009). “A New Scaling and Squaring Algorithm for the Matrix Exponential,” SIAM J. Matrix Anal. Applic. 31, 970–989.

A proof of Equation (9.3.1) for the scalar case appears in:

R.S. Varga (1961). “On Higher-Order Stable Implicit Methods for Solving Parabolic Partial Differential Equations,” J. Math. Phys. 40, 220–231.

There are many applications in control theory calling for the computation of the matrix exponential. In the linear optimal regular problem, for example, various integrals involving the matrix exponential are required, see:

J. Johnson and C.L. Phillips (1971). “An Algorithm for the Computation of the Integral of the State Transition Matrix,” IEEE Trans. Autom. Control AC-16, 204–205.

C.F. Van Loan (1978). “Computing Integrals Involving the Matrix Exponential,” IEEE Trans. Autom. Control AC-23, 395–404.

An understanding of the map $A \to \exp ( A t )$ and its sensitivity is helpful when assessing the performance of algorithms for computing the matrix exponential. Work in this direction includes:

B. K˚agstr¨om (1977). “Bounds and Perturbation Bounds for the Matrix Exponential,” BIT 17, 39–57.

C.F. Van Loan (1977). “The Sensitivity of the Matrix Exponential,” SIAM J. Numer. Anal. 14, 971–981.

R. Mathias (1992). “Evaluating the Fr´echet Derivative of the Matrix Exponential,” Numer. Math. 63, 213–226.

I. Najfeld and T.F. Havel (1995). “Derivatives of the Matrix Exponential and Their Computation,” Adv. Appl. Math. 16, 321–375.

A.H. Al-Mohy and N.J. Higham (2009). “Computing the Fr\`echet Derivative of the Matrix Exponential, with an Application to Condition Number Estimation,” SIAM J. Matrix Anal. Applic. 30, 1639– 1657.

A software package for computing small dense and large sparse matrix exponentials in Fortran and Matlab is presented in the following reference:

R.B. Sidje (1998) “Expokit: a Software Package for Computing Matrix Exponentials,” ACM Trans. Math. Softw. 24, 130–156.

Consideration of P9.3.2 and P9.3.5 shows that the exponential of a structured matrix can have important properties, see:

J. Xue and Q. Ye (2008). “Entrywise Relative Perturbation Bounds for Exponentials of Essentially Non-negative Matrices,” Numer. Math. 110, 393–403.

J. Cardoso and F.S. Leite (2010). “Exponentials of Skew-Symmetric Matrices and Logarithms of Orthogonal Matrices,” J. Comput. Appl. Math. 233, 2867–2875.
