# 9.2 Approximation Methods

We now consider a class of methods for computing matrix functions which at first glance do not appear to involve eigenvalues. These techniques are based on the idea that, if $g ( z )$ approximates f(z) on $\lambda ( A )$ , then f (A) approximates g(A), e.g.,

$$
e ^ {A} \approx I + A + \frac {A ^ {2}}{2 !} + \dots + \frac {A ^ {q}}{q !}.
$$

We begin by bounding $\parallel f ( A ) - g ( A ) \parallel$ using the Jordan and Schur matrix function representations. We follow this discussion with some comments on the evaluation of matrix polynomials.

# 9.2.1 A Jordan Analysis

The Jordan representation of matrix functions (Theorem 9.1.2) can be used to bound the error in an approximant $g ( A )$ of $f ( A )$ .

Theorem 9.2.1. Assume that

$$
A = X \cdot \mathrm{diag} (J _ {1}, \ldots , J _ {q}) \cdot X ^ {- 1}
$$

is the JCF of $A \in \mathbb { C } ^ { n \times n }$ with

$$
J _ {i} = \left[ \begin{array}{c c c c c} \lambda_ {i} & 1 & \dots & \dots & 0 \\ 0 & \lambda_ {i} & 1 & \vdots & \vdots \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ \vdots & \vdots & \vdots & \ddots & 1 \\ 0 & \dots & \dots & \dots & \lambda_ {i} \end{array} \right], \qquad n _ {i} \text {-by-} n _ {i},
$$

for i = 1:q. If f (z) and $g ( z )$ are analytic on an open set containing $\lambda ( A )$ , then

$$
\| f(A) - g(A)\|_{2}\leq \kappa_{2}(X)\max_{\substack{1\leq i\leq p\\ 0\leq r\leq n_{i} - 1}}n_{i}\frac{\left|f^{(r)}(\lambda_{i}) - g^{(r)}(\lambda_{i})\right|}{r!}.
$$

Proof. Defining $h ( z ) = f ( z ) - g ( z )$ we have

$$
\| f (A) - g (A) \| _ {2} = \| X \operatorname{diag} \left(h \left(J _ {1}\right), \dots , h \left(J _ {q}\right)\right) X ^ {- 1} \| _ {2} \leq \kappa_ {2} (X) \max _ {1 \leq i \leq q} \| h \left(J _ {i}\right) \| _ {2}.
$$

Using Theorem 9.1.2 and equation (2.3.8) we conclude that

$$
\| h (J _ {i}) \| _ {2} \leq n _ {i} \max _ {0 \leq r \leq n _ {i} - 1} \frac {| h ^ {(r)} (\lambda_ {i}) |}{r !}
$$

thereby proving the theorem.

# 9.2.2 A Schur Analysis

If we use the Schur decomposition $A = Q T Q ^ { H }$ instead of the Jordan decomposition, then the norm of $T \mathrm { s }$ strictly upper triangular portion is involved in the discrepancy between $f ( A )$ and $g ( A )$ .

Theorem 9.2.2. Let $Q ^ { H } A Q = T = \mathrm { d i a g } ( \lambda _ { i } ) + N$ be the Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ , with N being the strictly upper triangular portion of T . If $f ( z )$ and $g ( z )$ are analytic on a closed convex set Ω whose interior contains $\lambda ( A )$ , then

$$
\| f (A) - g (A) \| _ {F} \leq \sum_ {r = 0} ^ {n - 1} \delta_ {r} \frac {\| | N | ^ {r} \| _ {F}}{r !}
$$

where

$$
\delta_ {r} = \sup _ {z \in \Omega} \left| f ^ {(r)} (z) - g ^ {(r)} (z) \right|.
$$

Proof. Let $h ( z ) = f ( z ) - g ( z )$ and set $H = ( h _ { i j } ) = h ( A )$ . Let $S _ { i j } ^ { ( r ) }$ denote the set of strictly increasing integer sequences $\big ( s _ { 0 } , \ldots , s _ { r } \big )$ with the property that $s _ { 0 } = i$ and $s _ { r } = j$ . Notice that

$$
S _ {i j} = \bigcup_ {r = 1} ^ {j - i} S _ {i j} ^ {(r)}
$$

and so from Theorem 9.1.3, we obtain the following for all $i < j \colon$

$$
h _ {i j} = \sum_ {r = 1} ^ {j - 1} \sum_ {s \in S _ {i j} ^ {(r)}} n _ {s _ {0}, s _ {1}} n _ {s _ {1}, s _ {2}} \dots n _ {s _ {r - 1}, s _ {r}} h \left[ \lambda_ {s _ {0}}, \ldots , \lambda_ {s _ {r}} \right].
$$

Now since Ω is convex and h analytic, we have

$$
\left| h \left[ \lambda_ {s _ {0}}, \dots , \lambda_ {s _ {r}} \right] \right| \leq \sup _ {z \in \Omega} \frac {\left| h ^ {(r)} (z) \right|}{r !} = \frac {\delta_ {r}}{r !}. \tag {9.2.1}
$$

Furthermore if $| N | ^ { r } = ( n _ { i j } ^ { ( r ) } )$ for $r \geq 1$ , then it can be shown that

$$
n _ {i j} ^ {(r)} = \left\{ \begin{array}{l l} 0, & j <   i + r, \\ \sum_ {s \in S _ {i j} ^ {(r)}} \left| n _ {s _ {0}, s _ {1}} n _ {s _ {1}, s _ {2}} \dots n _ {s _ {r - 1}, s _ {r}} \right|, & j \geq i + r. \end{array} \right. \tag {9.2.2}
$$

The theorem now follows by taking absolute values in the expression for $h _ { i j }$ and then using (9.2.1) and (9.2.2).

There can be a pronounced discrepancy between the Jordan and Schur error bounds. For example, if

$$
A = \left[ \begin{array}{c c c} -. 0 1 & 1 & 1 \\ 0 & 0 & 1 \\ 0 & 0 & . 0 1 \end{array} \right].
$$

If $f ( z ) = e ^ { z }$ and $g ( z ) = 1 + z + z ^ { 2 } / 2$ , then  $f ( A ) - g ( A ) \parallel \approx \ 1 0 ^ { - 5 }$ in either the Frobenius norm or the 2-norm. Since $\kappa _ { 2 } ( X ) \approx 1 0 ^ { 7 }$ , the error predicted by Theorem 9.2.1 is $O ( 1 )$ , rather pessimistic. On the other hand, the error predicted by the Schur decomposition approach is $O ( 1 0 ^ { - 2 } )$ .

Theorems 9.2.1 and 9.2.2 remind us that approximating a function of a nonnormal matrix is more complicated than approximating a function of a scalar. In particular, we see that if the eigensystem of A is ill-conditioned and/or A’s departure from normality is large, then the discrepancy between $f ( A )$ and $g ( A )$ may be considerably larger than the maximum of $| f ( z ) - g ( z ) |$ on $\lambda ( A )$ . Thus, even though approximation methods avoid eigenvalue computations, they evidently appear to be influenced by the structure of A’s eigensystem. It is a perfect venue for pseudospectral analysis.

# 9.2.3 Taylor Approximants

A common way to approximate a matrix function such as $e ^ { A }$ is by truncating its Taylor series. The following theorem bounds the errors that arise when matrix functions such as these are approximated via truncated Taylor series.

Theorem 9.2.3. If $f ( z )$ has the Taylor series

$$
f (z) = \sum_ {k = 0} ^ {\infty} \alpha_ {k} z ^ {k}
$$

on an open disk containing the eigenvalues of $A \in \mathbb { C } ^ { n \times n }$ , then

$$
\left\| f (A) - \sum_ {k = 0} ^ {q} \alpha_ {k} A ^ {k} \right\| _ {2} \leq \frac {n}{(q + 1) !} \max _ {0 \leq s \leq 1} \| A ^ {q + 1} f ^ {(q + 1)} (A s) \| _ {2}.
$$

Proof. Define the matrix $E ( s )$ by

$$
f (A s) = \sum_ {k = 0} ^ {q} \alpha_ {k} (A s) ^ {k} + E (s), \quad 0 \leq s \leq 1. \tag {9.2.3}
$$

If $f _ { i j } ( s )$ is the $( i , j )$ entry of $f ( A s )$ , then it is necessarily analytic and so

$$
f _ {i j} (s) = \left(\sum_ {k = 0} ^ {q} \frac {f _ {i j} ^ {(k)} (0)}{k !} s ^ {k}\right) + \frac {f _ {i j} ^ {(q + 1)} (\varepsilon_ {i j})}{(q + 1) !} s ^ {q + 1} \tag {9.2.4}
$$

where $\varepsilon _ { i j }$ satisfies $0 \leq \varepsilon _ { i j } \leq s \leq 1$

By comparing powers of s in (9.2.3) and (9.2.4) we conclude that $e _ { i j } ( s )$ , the $( i , j )$ entry of $E ( s )$ , has the form

$$
e _ {i j} (s) = \frac {f _ {i j} ^ {(q + 1)} (\varepsilon_ {i j})}{(q + 1) !} s ^ {q + 1}.
$$

Now $f _ { i j } ^ { ( q - 1 ) } ( s )$ is the $( i , j )$ entry of $A ^ { q + 1 } f ^ { ( q + 1 ) } ( A s )$ and therefore

$$
| e _ {i j} (s) | \leq \max _ {0 \leq s \leq 1} \frac {f _ {i j} ^ {(q + 1)} (s)}{(q + 1) !} \leq \max _ {0 \leq s \leq 1} \frac {\| A ^ {q + 1} f ^ {(q + 1)} (A s) \| _ {2}}{(q + 1) !}.
$$

The theorem now follows by applying (2.3.8).

We mention that the factor of n in the upper bound can be removed with more careful analysis. See Mathias (1993).

In practice, it does not follow that greater accuracy results by taking a longer Taylor approximation. For example, if

$$
A = \left[ \begin{array}{l l} - 4 9 & 2 4 \\ - 6 4 & 3 1 \end{array} \right],
$$

then it can be shown that

$$
e ^ {A} = \left[ \begin{array}{l l} - 0. 7 3 5 7 5 9 & . 0 5 5 1 8 1 9 \\ - 1. 4 7 1 5 1 8 & 1. 1 0 3 6 3 8 \end{array} \right].
$$

For q = 59, Theorem 9.2.3 predicts that

$$
\left\| e ^ {A} - \sum_ {k = 0} ^ {q} \frac {A ^ {k}}{k !} \right\| _ {2} \leq \frac {n}{(q + 1) !} \max _ {0 \leq s \leq 1} \left\| A ^ {q + 1} e ^ {A s} \right\| _ {2} \leq 1 0 ^ {- 6 0}.
$$

However, if $\mathbf { u } \approx 1 0 ^ { - 7 }$ , then we find

$$
\mathsf {f l} \left(\sum_ {k = 0} ^ {5 9} \frac {A ^ {k}}{k !}\right) = \left[ \begin{array}{c c} - 2 2. 2 5 8 8 0 & - 1. 4 3 2 2 7 6 6 \\ - 6 1. 4 9 9 3 1 & - 3. 4 7 4 2 8 0 \end{array} \right].
$$

The problem is that some of the partial sums have large elements. For example, the matrix $I + A + \cdot \cdot \cdot + A ^ { 1 7 } / 1 7 !$ has entries of order $1 0 ^ { 7 }$ . Since the machine precision is approximately $1 0 ^ { - 7 }$ , rounding errors larger than the norm of the solution are sustained.


---

<!-- golub_550_599 -->

The example highlights the a well known shortcoming of truncated Taylor series approximation–it tends to be effcetive only near the origin. The problem can sometimes be circumvented through a change of scale. For example, by repeatedly using the double angle formulae

$$
\cos (2 A) = 2 \cos (A) ^ {2} - I, \quad \sin (2 A) = 2 \sin (A) \cos (A),
$$

the cosine and sine of a matrix can be built up from Taylor approximations to $\cos ( A / 2 ^ { k } )$ and sin $( A / 2 ^ { k } )$ :

$$
S _ {0} = \text { Taylor   approximate   to } \sin (A / 2 ^ {k})
$$

$$
C _ {0} = \text { Taylor   approximate   to } \cos (A / 2 ^ {k})
$$

$$
\text { for } j = 1: k
$$

$$
S _ {j} = 2 S _ {j - 1} C _ {j - 1}
$$

$$
C _ {j} = 2 C _ {j - 1} ^ {2} - I
$$

end

Here k is a positive integer chosen so that, say, $\| A \| _ { \infty } \approx 2 ^ { k }$ . See Serbin and Blalock (1979), Higham and Smith (2003), and Hargreaves and Higham (2005).

# 9.2.4 Evaluating Matrix Polynomials

Since the approximation of transcendental matrix functions usually involves the evaluation of polynomials, it is worthwhile to look at the details of computing

$$
p (A) = b _ {0} I + b _ {1} A + \dots + b _ {q} A ^ {q}
$$

where the scalars $b _ { 0 } , \dots , b _ { q } \in \mathbb { R }$ are given. The most obvious approach is to invoke Horner’s scheme:

Algorithm 9.2.1 Given a matrix A and $b ( 0 { : } q )$ , the following algorithm computes the polynomial $F = b _ { q } A ^ { q } + \cdot \cdot \cdot + b _ { 1 } A + b _ { 0 } I$ .

$$
F = b _ {q} A + b _ {q - 1} I
$$

$$
\text { for   } k = q - 2: - 1: 0
$$

$$
F = A F + b _ {k} I
$$

end

This requires $q - 1$ matrix multiplications. However, unlike the scalar case, this summation process is not optimal. To see why, suppose $q = 9$ and observe that

$$
p (A) = A ^ {3} (A ^ {3} (b _ {9} A ^ {3} + (b _ {8} A ^ {2} + b _ {7} A + b _ {6} I)) + (b _ {5} A ^ {2} + b _ {4} A + b _ {3} I)) + b _ {2} A ^ {2} + b _ {1} A + b _ {0} I.
$$

Thus, $F = p ( A )$ can be evaluated with only four matrix multiplications:

$$
\begin{array}{l} A _ {2} = A ^ {2}, \\ A _ {3} = A A _ {2}, \\ F _ {1} = b _ {9} A _ {3} + b _ {8} A _ {2} + b _ {7} A + b _ {6} I, \\ F _ {2} = A _ {3} F _ {1} + b _ {5} A _ {2} + b _ {4} A + b _ {3} I, \\ F = A _ {3} F _ {2} + b _ {2} A _ {2} + b _ {1} A + b _ {0} I. \\ \end{array}
$$

In general, if s is any integer that satisfies $1 \leq s \leq { \sqrt { q } }$ , then

$$
p (A) = \sum_ {k = 0} ^ {r} B _ {k} \cdot (A ^ {s}) ^ {k}, \quad r = \text { floor } (q / s), \tag {9.2.5}
$$

where

$$
B _ {k} = \left\{ \begin{array}{l l} b _ {s k + s - 1} A ^ {s - 1} + \dots + b _ {s k + 1} A + b _ {s k} I, & k = 0: r - 1, \\ b _ {q} A ^ {q - s r} + \dots + b _ {s r + 1} A + b _ {s r} I, & k = r. \end{array} \right.
$$

After $A ^ { 2 } , \ldots , A ^ { s }$ are computed, then Horner’s rule can be applied to (9.2.5) and the net result is that $p ( A )$ can be computed with $s + r - 1$ matrix multiplications. By choosing $s = \mathsf { f l o o r } ( { \sqrt { q } } )$ , the number of matrix multiplications is approximately minimized. This technique is discussed by Paterson and Stockmeyer (1973). Van Loan (1978) shows how the procedure can be implemented without storage arrays for $A ^ { 2 } , \ldots , A ^ { s }$ .

# 9.2.5 Computing Powers of a Matrix

The problem of raising a matrix to a given power deserves special mention. Suppose it is required to compute $A ^ { 1 3 }$ . Noting that $A ^ { \bar { 4 } } = ( A ^ { 2 } ) ^ { 2 } , A ^ { 8 } = ( A ^ { 4 } ) ^ { 2 }$ , and $A ^ { 1 3 } = \bar { A } ^ { 8 } A ^ { 4 } A$ , we see that this can be accomplished with just five matrix multiplications. In general we have

Algorithm 9.2.2 (Binary Powering) The following algorithm computes $F = A ^ { s }$ where s is a positive integer and A ∈ IRn×n. $A \in \mathbb { R } ^ { n \times n }$

Let $s = \sum _ { k = 0 } ^ { t } \beta _ { k } 2 ^ { k }$ be the binary expansion of s with $\beta _ { t } \neq 0$

$$
Z = A; q = 0
$$

while $\beta _ { q } = 0$

$$
Z = Z ^ {2}; q = q + 1
$$

end

$$
F = Z
$$

for $k = q + 1 { : } t$

$$
Z = Z ^ {2}
$$

$\mathbf { i f } \ \beta _ { k } \neq 0$

$$
F = F Z
$$

end

end

This algorithm requires at most $2 \mathsf { f l o o r } [ \log _ { 2 } ( s ) ]$ matrix multiplications. If s is a power of 2, then only $\log _ { 2 } ( s )$ matrix multiplications are needed.

# 9.2.6 Integrating Matrix Functions

We conclude this section with some remarks about the integration of a parameterized matrix function. Suppose $A \in \mathbb { R } ^ { n \times n }$ and that $f ( A t )$ is defined for all $t \in [ a , b ]$ . We can

approximate

$$
F = \int_ {a} ^ {b} f (A t) d t \quad \Leftrightarrow \quad [ F ] _ {i j} = \int_ {a} ^ {b} [ f (A t) ] _ {i j} d t
$$

by applying any suitable quadrature rule. For example, with Simpson’s rule, we have

$$
F \approx \tilde {F} = \frac {h}{3} \sum_ {k = 0} ^ {m} w _ {k} f (A (a + k h)) \tag {9.2.6}
$$

where m is even, $h = ( b - a ) / m$ , and

$$
w _ {k} = \left\{ \begin{array}{l l} 1 & k = 0, m, \\ 4 & k \text {odd}, \\ 2 & k \text {even}, k \neq 0, m. \end{array} \right.
$$

If $( d ^ { 4 } / d z ^ { 4 } ) f ( z t ) = f ^ { ( 4 ) } ( z t )$ is continuous for $t \in [ a , b ]$ and if $f ^ { ( 4 ) } ( A t )$ is defined on this same interval, then it can be shown that $\tilde { \boldsymbol { F } } = \boldsymbol { F } + \boldsymbol { E }$ where

$$
\| E \| _ {2} \leq \frac {n h ^ {4} (b - a)}{1 8 0} \max _ {a \leq t \leq b} \| f ^ {(4)} (A t) \| _ {2}. \tag {9.2.7}
$$

Let $f _ { i j }$ and $e _ { i j }$ denote the $( i , j )$ entries of F and E, respectively. Under the above assumptions we can apply the standard error bounds for Simpson’s rule and obtain

$$
| e _ {i j} | \leq \frac {h ^ {4} (b - a)}{1 8 0} \max _ {a \leq t \leq b} | e _ {i} ^ {T} f ^ {(4)} (A t) e _ {j} |.
$$

The inequality (9.2.7) now follows since $\parallel E \parallel _ { 2 } \leq n \operatorname* { m a x } | e _ { i j } |$ and

$$
\max _ {a \leq t \leq b} | e _ {i} ^ {T} f ^ {(4)} (A t) e _ {j} | \leq \max _ {a \leq t \leq b} \| f ^ {(4)} (A t) \| _ {2}.
$$

Of course, in a practical application of (9.2.6), the function evaluations $f ( A ( a + k h ) )$ normally have to be approximated. Thus, the overall error involves the error in approximating $f ( A ( a + k h )$ as well as the Simpson rule error.

# 9.2.7 A Note on the Cauchy Integral Formulation

Yet another way to define a function of a matrix $C \in \mathbb { C } ^ { n \times n }$ is through the Cauchy integral theorem. Suppose $f ( z )$ is analytic inside and on a closed contour Γ which encloses $\lambda ( A )$ . We can define $f ( A )$ to be the matrix

$$
f (A) = \frac {1}{2 \pi i} \oint_ {\Gamma} f (z) (z I - A) ^ {- 1} d z. \tag {9.2.8}
$$

The integral is defined on an element-by-element basis:

$$
f (A) = \left(f _ {k j}\right) \quad \Longrightarrow \quad f _ {k j} = \frac {1}{2 \pi i} \oint_ {\Gamma} f (z) e _ {k} ^ {T} (z I - A) ^ {- 1} e _ {j} d z.
$$

Notice that the entries of $( z I { - } A ) ^ { - 1 }$ are analytic on Γ and that $f ( A )$ is defined whenever $f ( z )$ is analytic in a neighborhood of $\lambda ( A )$ . Using quadrature and other tools, Hale, Higham, and Trefethen (2007) have shown how this characterization can be used in practice to compute certain types of matrix functions.

# Problems

P9.2.1 Verify (9.2.2).

P9.2.2 Show that if $\| \ A \| _ { 2 } < 1$ , then $\log ( I + A )$ exists and satisfies the bound

$$
\| \log (I + A) \| _ {2} \leq \| A \| _ {2} / (1 - \| A \| _ {2}).
$$

P9.2.3 Using Theorem 9.2.3, bound the error in the following approximations:

$$
\sin (A) \approx \sum_ {k = 0} ^ {q} (- 1) ^ {k} \frac {A ^ {2 k + 1}}{(2 k + 1) !}, \qquad \cos (A) \approx \sum_ {k = 0} ^ {q} (- 1) ^ {k} \frac {A ^ {2 k}}{(2 k) !}.
$$

P9.2.4 Suppose $A \in \mathbb { R } ^ { n \times n }$ is nonsingular and $X _ { 0 } \in \mathbb { R } ^ { n \times n }$ is given. The iteration defined by

$$
X _ {k + 1} = X _ {k} (2 I - A X _ {k})
$$

is the matrix analogue of Newton’s method applied to the function $f ( x ) = a - ( 1 / x )$ . Use the SVD to analyze this iteration. Do the iterates converge to $A ^ { - 1 \mathord { \ ? } }$ Discuss the choice of $X _ { 0 }$ .

P9.2.5 Assume $A \in \mathbb { R } ^ { 2 \times 2 }$ . (a) Specify real scalars α and $\beta$ so that $A ^ { 4 } = \alpha I + \beta A$ . (b) Develop recursive recipes for $\alpha _ { k }$ and $\beta _ { k }$ so that $A ^ { k } = \alpha _ { k } I + \beta _ { k } A$ for $k \geq 2$ .

# Notes and References for 9.2

The optimality of Horner’s rule for polynomial evaluation is discussed in:

M.S. Paterson and L.J. Stockmeyer (1973). “On the Number of Nonscalar Multiplications Necessary to Evaluate Polynomials,” SIAM J. Comput. 2, 60–66.   
D.E. Knuth (1981). The Art of Computer Programming, Vol. 2. Seminumerical Algorithms, second edition, Addison-Wesley, Reading, MA.   
The Horner evaluation of matrix polynomials is analyzed in:   
C.F. Van Loan (1978). “A Note on the Evaluation of Matrix Polynomials,” IEEE Trans. Autom. Control AC-24, 320–321.   
Other aspects of matrix function approximation and evaluation are discussed in:   
H. Bolz and W. Niethammer (1988). “On the Evaluation of Matrix Functions Given by Power Series,” SIAM J. Matrix Anal. Applic. 9, 202–209.   
R. Mathias (1993). “Approximation of Matrix-Valued Functions,” SIAM J. Matrix Anal. Applic. 14, 1061–1063.   
N.J. Higham and P.A. Knight (1995). “Matrix Powers in Finite Precision Arithmetic,” SIAM J. Matrix Anal. Applic. 16, 343–358.   
P. Sebastiani (1996). “On the Derivatives of Matrix Powers,” SIAM J. Matrix Anal. Applic. 17, 640–648.   
D.S. Bernstein and C.F. Van Loan (2000). “Rational Matrix Functions and Rank-One Updates,” SIAM J. Matrix Anal. Applic. 22, 145–154.   
For a discussion of methods for computing the sine and cosine of a matrix, see:   
S. Serbin and S. Blalock (1979). “An Algorithm for Computing the Matrix Cosine,” SIAM J. Sci. Stat. Comput. 1, 198–204.   
N.J. Higham and M.I. Smit (2003). “Computing the Matrix Cosine,” Numer. Algorithms 34, 13–26.   
G. Hargreaves and N.J. Higham (2005). “Efficient Algorithms for the Matrix Cosine and Sine,” Numer. Algorithms 40, 383–400.   
The computation of f(A) using contour integrals is analyzed in:   
N. Hale, N.J. Higham, and L.N. Trefethen (2007). “Computing Aα, log(A), and Related Matrix Functions by Contour Integrals,” SIAM J. Numer. Anal. 46, 2505–2523.
