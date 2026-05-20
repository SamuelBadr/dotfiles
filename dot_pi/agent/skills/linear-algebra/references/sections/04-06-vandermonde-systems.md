# 4.6 Vandermonde Systems

Suppose x $( 0 : n ) \in \mathbb { R } ^ { n + 1 }$ . A matrix $V \in \mathbb { R } ^ { ( n + 1 ) \times ( n + 1 ) }$ of the form

$$
V = V (x _ {0}, \ldots , x _ {n}) = \left[ \begin{array}{c c c c} 1 & 1 & \dots & 1 \\ x _ {0} & x _ {1} & \dots & x _ {n} \\ \vdots & \vdots & & \vdots \\ x _ {0} ^ {n} & x _ {1} ^ {n} & \dots & x _ {n} ^ {n} \end{array} \right]
$$

is said to be a Vandermonde matrix. Note that the discrete Fourier transform matrix (§1.4.1) is a very special complex Vandermonde matrix.

In this section, we show how the systems $V ^ { T } a = f = f ( 0 { : } n )$ and $V z = b = b ( 0 { : } n )$ can be solved in $O ( n ^ { 2 } )$ flops. For convenience, vectors and matrices are subscripted from 0 in this section.

# 4.6.1 Polynomial Interpolation: $V ^ { T } a = f$

Vandermonde systems arise in many approximation and interpolation problems. Indeed, the key to obtaining a fast Vandermonde solver is to recognize that solving $V ^ { T } a = f$ is equivalent to polynomial interpolation. This follows because if $V ^ { T } a = f$ and

$$
p (x) = \sum_ {j = 0} ^ {n} a _ {j} x ^ {j}, \tag {4.6.1}
$$

then $p ( x _ { i } ) = f _ { i }$ for $i = 0 { : } n$ .

Recall that if the $x _ { i }$ are distinct then there is a unique polynomial of degree n that interpolates $( x _ { 0 } , f _ { 0 } ) , \ldots , ( x _ { n } , f _ { n } )$ . Consequently, V is nonsingular as long as the $x _ { i }$ are distinct. We assume this throughout the section.

The first step in computing the $a _ { j }$ of (4.6.1) is to calculate the Newton representation of the interpolating polynomial $p \mathrm { : }$

$$
p (x) = \sum_ {k = 0} ^ {n} c _ {k} \left(\prod_ {i = 0} ^ {k - 1} (x - x _ {i})\right). \tag {4.6.2}
$$

The constants $c _ { k }$ are divided differences and may be determined as follows:

$$
c (0: n) = f (0: n)
$$

$$
\text { for } k = 0: n - 1
$$

$$
\text { for   } i = n: - 1: k + 1 \tag {4.6.3}
$$

$$
c _ {i} = \left(c _ {i} - c _ {i - 1}\right) / \left(x _ {i} - x _ {i - k - 1}\right)
$$

end

end

See Conte and deBoor (1980).

The next task is to generate the coefficients $a _ { 0 } , \ldots , a _ { n }$ in (4.6.1) from the Newton representation coefficients $c _ { 0 } , \ldots , c _ { n }$ . Define the polynomials $p _ { n } ( x ) , \ldots , p _ { 0 } ( x )$ by the iteration

$$
p _ {n} (x) = c _ {n}
$$

for k = n − 1 : −1 : 0

$$
p _ {k} (x) = c _ {k} + (x - x _ {k}) \cdot p _ {k + 1} (x)
$$

end

and observe that $p _ { 0 } ( x ) = p ( x )$ . Writing

$$
p _ {k} (x) = a _ {k} ^ {(k)} + a _ {k + 1} ^ {(k)} x + \dots + a _ {n} ^ {(k)} x ^ {n - k}
$$

and equating like powers of x in the equation $p _ { k } = c _ { k } + ( x - x _ { k } ) p _ { k + 1 }$ gives the following recursion for the coefficients $a _ { i } ^ { ( k ) }$ (k)

$$
a _ {n} ^ {(n)} = c _ {n}
$$

for k = n−1: −1 : 0

$$
a _ {k} ^ {(k)} = c _ {k} - x _ {k} a _ {k + 1} ^ {(k + 1)}
$$

for i = k + 1 : n − 1

$$
a _ {i} ^ {(k)} = a _ {i} ^ {(k + 1)} - x _ {k} a _ {i + 1} ^ {(k + 1)}
$$

end

$$
a _ {n} ^ {(k)} = a _ {n} ^ {(k + 1)}
$$

end

Consequently, the coefficients $a _ { i } = a _ { i } ^ { ( 0 ) }$ can be calculated as follows:

$$
a (0: n) = c (0: n)
$$

for k = n−1: −1 : 0

for i = k:n − 1 (4.6.4)

$$
a _ {i} = a _ {i} - x _ {k} a _ {i + 1}
$$

end

end

Combining this iteration with (4.6.3) gives the following algorithm.

Algorithm 4.6.1 Given $x ( 0 : n ) \in \mathbb { R } ^ { n + 1 }$ with distinct entries and $f = f ( 0 : n ) \in \mathbb { R } ^ { n + 1 }$ , the following algorithm overwrites f with the solution $a = a ( 0 : n )$ to the Vandermonde system $V ( x _ { 0 } , \ldots , x _ { n } ) ^ { T } a = f .$ .

for k = 0 : n − 1

for i = n: −1 :k + 1

$$
f (i) = (f (i) - f (i - 1)) / (x (i) - x (i - k - 1))
$$

end

end

for k = n − 1: −1 : 0

for i = k : n − 1

$$
f (i) = f (i) - f (i + 1) \cdot x (k)
$$

end

end

This algorithm requires $5 n ^ { 2 } / 2$ flops.

# 4.6.2 The System $V z = b$

Now consider the system $V z = b$ . To derive an efficient algorithm for this problem, we describe what Algorithm 4.6.1 does in matrix-vector language. Define the lower bidiagonal matrix $L _ { k } ( \alpha ) \in \mathbb { R } ^ { ( n + 1 ) \times ( n + 1 ) }$ by

$$
L _ {k} (\alpha) = \left[ \begin{array}{c c c c c c c} I _ {k} & & & 0 & & & \\ \hline & 1 & 0 & \dots & & & 0 \\ & - \alpha & 1 & & & & \\ & 0 & \ddots & \ddots & & & \\ 0 & \vdots & & \ddots & \ddots & & \vdots \\ & & & & \ddots & 1 & \\ & 0 & & \dots & & - \alpha & 1 \end{array} \right]
$$

and the diagonal matrix $D _ { k }$ by

$$
D _ {k} = \operatorname{diag} (\underbrace {1 , \ldots , 1} _ {k + 1}, x _ {k + 1} - x _ {0}, \ldots , x _ {n} - x _ {n - k - 1}).
$$

With these definitions it is easy to verify from (4.6.3) that, if $f = f ( 0 : n ) \operatorname { a n d } c = c ( 0 : n )$ is the vector of divided differences, then

$$
c = U ^ {T} f
$$

where U is the upper triangular matrix defined by

$$
U ^ {T} = D _ {n - 1} ^ {- 1} L _ {n - 1} (1) \dots D _ {0} ^ {- 1} L _ {0} (1).
$$

Similarly, from (4.6.4) we have

$$
a = L ^ {T} c,
$$

where L is the unit lower triangular matrix defined by

$$
L ^ {T} = L _ {0} (x _ {0}) ^ {T} \dots L _ {n - 1} (x _ {n - 1}) ^ {T}.
$$

It follows that $a = V ^ { - T } f$ is given by

$$
a = L ^ {T} U ^ {T} f.
$$

Thus,

$$
V ^ {- T} = L ^ {T} U ^ {T}
$$

which shows that Algorithm 4.6.1 solves $V ^ { T } a ~ = ~ f$ by tacitly computing the $^ { 6 6 } \mathrm { U L }$ factorization” of $V ^ { - 1 }$ . Consequently, the solution to the system $V z = b$ is given by

$$
\begin{array}{l} z = V ^ {- 1} b = U (L b) \\ = \left(L _ {0} (1) ^ {T} D _ {0} ^ {- 1} \dots L _ {n - 1} (1) ^ {T} D _ {n - 1} ^ {- 1}\right) \left(L _ {n - 1} (x _ {n - 1}) \dots L _ {0} (x _ {0}) b\right). \\ \end{array}
$$

This observation gives rise to the following algorithm:

Algorithm 4.6.2 Given $x ( 0 : n ) \in \mathbb { R } ^ { n + 1 }$ with distinct entries and $b = b ( 0 : n ) \in \mathbb { R } ^ { n + 1 }$ , the following algorithm overwrites b with the solution $z = z ( 0 : n )$ to the Vandermonde system $V ( x _ { 0 } , \ldots , x _ { n } ) z = b .$

for k = 0:n - 1
    for i = n: -1:k + 1
    b(i) = b(i) - x(k)b(i - 1)
    end
end
for k = n - 1: -1: 0
    for i = k + 1:n
    b(i) = b(i)/(x(i) - x(i - k - 1))
    end
    for i = k:n - 1
    b(i) = b(i) - b(i + 1)
    end
end

This algorithm requires $5 n ^ { 2 } / 2$ flops.

Algorithms 4.6.1 and 4.6.2 are discussed and analyzed by Bj¨orck and Pereyra (1970). Their experience is that these algorithms frequently produce surprisingly accurate solutions, even if V is ill-conditioned.

We mention that related techniques have been developed and analyzed for confluent Vandermonde systems, e.g., systems of the form

$$
\left[ \begin{array}{c c c c} 1 & 1 & 0 & 1 \\ x _ {0} & x _ {1} & 1 & x _ {3} \\ x _ {0} ^ {2} & x _ {1} ^ {2} & 2 x _ {1} & x _ {3} ^ {2} \\ x _ {0} ^ {3} & x _ {1} ^ {3} & 3 x _ {1} ^ {2} & x _ {3} ^ {3} \end{array} \right] ^ {T} \left[ \begin{array}{c} a _ {0} \\ a _ {1} \\ a _ {2} \\ a _ {3} \end{array} \right] = \left[ \begin{array}{c} f _ {0} \\ f _ {1} \\ f _ {2} \\ f _ {3} \end{array} \right].
$$

See Higham (1990).

# Problems

P4.6.1 Show that if $V = V ( x _ { 0 } , \ldots , x _ { n } )$ , then

$$
\det (V) = \prod_ {n \geq i > j \geq 0} (x _ {i} - x _ {j}).
$$

P4.6.2 (Gautschi 1975) Verify the following inequality for the n = 1 case above:

$$
\| V^{-1}\|_{\infty}\leq \max_{0\leq k\leq n}\prod_{\substack{i = 0\\ i\neq k}}^{n}\frac{1 + |x_{i}|}{|x_{k} - x_{i}|}  .
$$

Equality results if the $x _ { i }$ are all on the same ray in the complex plane.

# Notes and References for §4.6

Our discussion of Vandermonde linear systems is drawn from the following papers:

The divided difference computations we discussed are detailed in:   
Error analyses of Vandermonde system solvers include:   
Interesting theoretical results concerning the condition of Vandermonde systems may be found in:   
A. Bj¨orck and V. Pereyra (1970). “Solution of Vandermonde Systems of Equations,”Math. Comput. 24, 893–903.   
A. Bj¨orck and T. Elfving (1973). “Algorithms for Confluent Vandermonde Systems,” Numer. Math. 21, 130–37.   
S.D. Conte and C. de Boor (1980). Elementary Numerical Analysis: An Algorithmic Approach, Third Edition, McGraw-Hill, New York, Chapter 2.   
N.J. Higham (1987). “Error Analysis of the Bj¨orck-Pereyra Algorithms for Solving Vandermonde Systems,” Numer. Math. 50, 613–632.   
N.J. Higham (1988). “Fast Solution of Vandermonde-like Systems Involving Orthogonal Polynomials,” IMA J. Numer. Anal. 8, 473–486.   
N.J. Higham (1990). “Stability Analysis of Algorithms for Solving Confluent Vandermonde-like Systems,” SIAM J. Matrix Anal. Applic. 11, 23–41.   
S.G. Bartels and D.J. Higham (1992). “The Structured Sensitivity of Vandermonde-Like Systems,” Numer. Math. 62, 17–34.   
J.M. Varah (1993). “Errors and Perturbations in Vandermonde Systems,” IMA J. Numer. Anal. 13, 1–12.   
W. Gautschi (1975). “Norm Estimates for Inverses of Vandermonde Matrices,”Numer. Math. 23, 337–347.   
W. Gautschi (1975). “Optimally Conditioned Vandermonde Matrices,” Numer. Math. 24, 1–12.   
J-G. Sun (1998). “Bounds for the Structured Backward Errors of Vandermonde Systems,” SIAM J. Matrix Anal. Applic. 20, 45–59.   
B.K. Alpert (1996). “Condition Number of a Vandermonde Matrix,” SIAM Review 38, 314–314.   
B. Beckermann (2000). “The condition number of real Vandermonde, Krylov and positive definite Hankel matrices,” Numer. Math. 85, 553–577.   
The basic algorithms presented can be extended to cover confluent Vandermonde systems, block Vandermonde systems, and Vandermonde systems with other polynomial bases:   
G. Galimberti and V. Pereyra (1970). “Numerical Differentiation and the Solution of Multidimensional Vandermonde Systems,” Math. Comput. 24, 357–364.   
G. Galimberti and V. Pereyra (1971). “Solving Confluent Vandermonde Systems of Hermitian Type,” Numer. Math. 18, 44–60.   
H. Van de Vel (1977). “Numerical Treatment of Generalized Vandermonde Systems of Equations,” Lin. Alg. Applic. 17, 149–174.   
G.H. Golub and W.P Tang (1981). “The Block Decomposition of a Vandermonde Matrix and Its Applications,”BIT 21, 505–517.   
D. Calvetti and L. Reichel (1992). “A Chebychev-Vandermonde Solver,” Lin. Alg. Applic. 172, 219–229.   
D. Calvetti and L. Reichel (1993). “Fast Inversion of Vandermonde-Like Matrices Involving Orthogonal Polynomials,” BIT 33, 473–484.   
H. Lu (1994). “Fast Solution of Confluent Vandermonde Linear Systems,” SIAM J. Matrix Anal. Applic. 15, 1277–1289.   
H. Lu (1996). “Solution of Vandermonde-like Systems and Confluent Vandermonde-like Systems,” SIAM J. Matrix Anal. Applic. 17, 127–138.   
M.-R. Skrzipek (2004). “Inversion of Vandermonde-Like Matrices,” BIT 44, 291–306.   
J.W. Demmel and P. Koev (2005). “The Accurate and Efficient Solution of a Totally Positive Generalized Vandermonde Linear System,” SIAM J. Matrix Anal. Applic. 27, 142–152.   
The displacement rank idea that we discuss in §12.1 can also be used to develop fast methods for Vandermonde systems.
