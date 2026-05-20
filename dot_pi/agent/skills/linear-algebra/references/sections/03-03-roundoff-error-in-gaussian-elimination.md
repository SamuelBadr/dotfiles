# 3.3 Roundoff Error in Gaussian Elimination

We now assess the effect of rounding errors when the algorithms in the previous two sections are used to solve the linear system $A x = b$ . A much more detailed treatment of roundoff error in Gaussian elimination is given in Higham (ASNA).

# 3.3.1 Errors in the LU Factorization

Let us see how the error bounds for Gaussian elimination compare with the ideal bounds derived in §2.7.11. We work with the infinity norm for convenience and focus our attention on Algorithm 3.2.1, the outer product version. The error bounds that we derive also apply to the gaxpy formulation (Algorithm 3.2.2). Our first task is to quantify the roundoff errors associated with the computed triangular factors.

Theorem 3.3.1. Assume that A is an n-by-n matrix of floating point numbers. If no zero pivots are encountered during the execution of Algorithm 3.2.1, then the computed triangular matrices $\hat { L }$ and Uˆ satisfy

$$
\hat {L} \hat {U} = A + H, \tag {3.3.1}
$$

$$
| H | \leq 2 (n - 1) \mathbf {u} \left(| A | + | \hat {L} | | \hat {U} |\right) + O \left(\mathbf {u} ^ {2}\right). \tag {3.3.2}
$$

Proof. The proof is by induction on n. The theorem obviously holds for $n = 1$ . Assume that $n \geq 2$ and that the theorem holds for all $( n - 1 ) – \mathrm { b y } – ( n - 1 )$ floating point matrices. If A is partitioned as follows

$$
A = \left[ \begin{array}{c c} \alpha & w ^ {T} \\ v & B \end{array} \right] _ {n - 1} ^ {1}
$$

then the first step in Algorithm 3.2.1 is to compute

$$
\hat {z} = \mathsf {f l} (v / \alpha), \qquad \hat {C} = \mathsf {f l} (\hat {z} w ^ {T}), \qquad \hat {A} _ {1} = \mathsf {f l} (B - \hat {C}),
$$

from which we conclude that

$$
\hat {z} = v / \alpha + f, \tag {3.3.3}
$$

$$
| f | \leq \mathbf {u} | v / \alpha |, \tag {3.3.4}
$$

$$
\hat {C} = \hat {z} w ^ {T} + F _ {1}, \tag {3.3.5}
$$

$$
\left| F _ {1} \right| \leq \mathbf {u} \left| \hat {z} \right| \left| w ^ {T} \right|, \tag {3.3.6}
$$

$$
\hat {A} _ {1} = B - \left(\hat {z} w ^ {T} + F _ {1}\right) + F _ {2}, \tag {3.3.7}
$$

$$
\left| F _ {2} \right| \leq \mathbf {u} \left(\left| B \right| + \left| \hat {z} \right| \left| w ^ {T} \right|\right) + O \left(\mathbf {u} ^ {2}\right), \tag {3.3.8}
$$

$$
| \hat {A} _ {1} | \leq | B | + | \hat {z} | | w ^ {T} | + O (\mathbf {u}). \tag {3.3.9}
$$

The algorithm proceeds to compute the LU factorization of $\hat { A } _ { 1 }$ . By induction, the computed factors $\hat { L } _ { 1 }$ and $\hat { U } _ { 1 }$ satisfy

$$
\hat {L} _ {1} \hat {U} _ {1} = \hat {A} _ {1} + H _ {1} \tag {3.3.10}
$$

where

$$
\left| H _ {1} \right| \leq 2 (n - 2) \mathbf {u} \left(\left| \hat {A} _ {1} \right| + \left| \hat {L} _ {1} \right| \left| \hat {U} _ {1} \right|\right) + O \left(\mathbf {u} ^ {2}\right). \tag {3.3.11}
$$

If

$$
\hat {L} = \left[ \begin{array}{c c} 1 & 0 \\ \hat {z} & \hat {L} _ {1} \end{array} \right], \qquad \hat {U} = \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & \hat {U} _ {1} \end{array} \right],
$$

then it is easy to verify that

$$
\hat {L} \hat {U} = A + H
$$

where

$$
H = \left[ \begin{array}{c c} 0 & 0 \\ \alpha f & H _ {1} - F _ {1} + F _ {2} \end{array} \right]. \tag {3.3.12}
$$

To prove the theorem we must verify (3.3.2), i.e.,

$$
| H | \leq 2 (n - 1) \mathbf {u} \left[ \begin{array}{c c} 2 | \alpha | & 2 | w ^ {T} | \\ | v | + | \alpha | | f | & | B | + | \hat {L} _ {1} | | \hat {U} _ {1} | + | \hat {z} | | w ^ {T} | \end{array} \right] + O (\mathbf {u} ^ {2}).
$$

Considering (3.3.12), this is obviously the case if

$$
\left| H _ {1} \right| + \left| F _ {1} \right| + \left| F _ {2} \right| \leq 2 (n - 1) \mathbf {u} \left(\left| B \right| + \left| \hat {z} \right| \left| w ^ {T} \right| + \left| \hat {L} _ {1} \right| \left| \hat {U} _ {1} \right|\right) + O \left(\mathbf {u} ^ {2}\right). \tag {3.3.13}
$$

Using (3.3.9) and (3.3.11) we have

$$
| H _ {1} | \leq 2 (n - 2) \mathbf {u} \left(| B | + | \hat {z} | | w ^ {T} | + | \hat {L} _ {1} | | \hat {U} _ {1} |\right) + O (\mathbf {u} ^ {2}),
$$

while (3.3.6) and (3.3.8) imply

$$
| F _ {1} | + | F _ {2} | \leq \mathbf {u} (| B | + 2 | \hat {z} | | w |) + O (\mathbf {u} ^ {2}).
$$

These last two results establish (3.3.13) and therefore the theorem.

We mention that if A is m-by-n, then the theorem applies with n replaced by the smaller of n and m in Equation 3.3.2.

# 3.3.2 Triangular Solving with Inexact Triangles

We next examine the effect of roundoff error when $\hat { L }$ and $\hat { U }$ are used by the triangular system solvers of §3.1.

Theorem 3.3.2. Let $\hat { L }$ and $\hat { U }$ be the computed LU factors obtained by Algorithm 3.2.1 when it is applied to an n-by-n floating point matrix A. If the methods of §3.1 are used to produce the computed solution yˆ to Lyˆ = b and the computed solution xˆ to $\hat { U } x = \hat { y }$ , then $( { \boldsymbol { A } } + { \boldsymbol { E } } ) { \hat { \boldsymbol { x } } } = { \boldsymbol { b } }$ with

$$
| E | \leq n \mathbf {u} (2 | A | + 4 | \hat {L} | | \hat {U} |) + O \left(\mathbf {u} ^ {2}\right). \tag {3.3.14}
$$

Proof. From (3.1.1) and (3.1.2) we have

$$
(\hat {L} + F) \hat {y} = b, \quad | F | \leq n \mathbf {u} | \hat {L} | + O (\mathbf {u} ^ {2}),
$$

$$
(\hat {U} + G) \hat {x} = \hat {y}, \quad | G | \leq n \mathbf {u} | \hat {U} | + O (\mathbf {u} ^ {2}),
$$

and thus

$$
(\hat {L} + F) (\hat {U} + G) \hat {x} = (\hat {L} \hat {U} + F \hat {U} + \hat {L} G + F G) \hat {x} = b.
$$

If follows from Theorem 3.3.1 that $\hat { L } \hat { U } = A + H$ with

$$
| H | \leq 2 (n - 1) \mathbf {u} (| A | + | \hat {L} | | \hat {U} |) + O (\mathbf {u} ^ {2}),
$$

and so by defining

$$
E = H + F \hat {U} + \hat {L} G + F G
$$

we find (A + E)ˆx = b. Moreover,

$$
| E | \leq | H | + | F | | \hat {U} | + | \hat {L} | | G | + O (\mathbf {u} ^ {2})
$$

$$
\leq 2 n \mathbf {u} \left(| A | + | \hat {L} | | \hat {U} |\right) + 2 n \mathbf {u} \left(| \hat {L} | | \hat {U} |\right) + O (\mathbf {u} ^ {2}),
$$

completing the proof of the theorem.

If it were not for the possibility of a large $| \hat { L } | | \hat { U } |$ term, (3.3.14) would compare favorably with the ideal bound (2.7.21). (The factor n is of no consequence, cf. the Wilkinson quotation in §2.7.7.) Such a possibility exists, for there is nothing in Gaussian elimination to rule out the appearance of small pivots. If a small pivot is encountered, then we can expect large numbers to be present in $\hat { L }$ and $\hat { U }$ .

We stress that small pivots are not necessarily due to ill-conditioning as the example

$$
A = \left[ \begin{array}{c c} \epsilon & 1 \\ 1 & 0 \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ 1 / \epsilon & 1 \end{array} \right] \left[ \begin{array}{c c} \epsilon & 1 \\ 0 & - 1 / \epsilon \end{array} \right]
$$

shows. Thus, Gaussian elimination can give arbitrarily poor results, even for wellconditioned problems. The method is unstable. For example, suppose 3-digit floating point arithmetic is used to solve

$$
\left[ \begin{array}{l l}. 0 0 1 & 1. 0 0 \\ 1. 0 0 & 2. 0 0 \end{array} \right] \left[ \begin{array}{l} x _ {1} \\ x _ {2} \end{array} \right] = \left[ \begin{array}{l} 1. 0 0 \\ 3. 0 0 \end{array} \right].
$$

(See §2.7.1.) Applying Gaussian elimination we get

$$
\hat {L} = \left[ \begin{array}{c c} 1 & 0 \\ 1 0 0 0 & 1 \end{array} \right], \qquad \hat {U} = \left[ \begin{array}{c c}. 0 0 1 & 1 \\ 0 & - 1 0 0 0 \end{array} \right],
$$

and a calculation shows that

$$
\hat {L} \hat {U} = \left[ \begin{array}{c c}. 0 0 1 & 1 \\ 1 & 2 \end{array} \right] + \left[ \begin{array}{c c} 0 & 0 \\ 0 & - 2 \end{array} \right] \equiv A + H.
$$

If we go on to solve the problem using the triangular system solvers of §3.1, then using the same precision arithmetic we obtain a computed solution $\hat { x } = [ 0 \ , \ 1 ] ^ { T }$ . This is in contrast to the exact solution $\boldsymbol { x } = [ 1 . 0 0 2 \ldots , . 9 9 8 \ldots ] ^ { T }$ .

# Problems

P3.3.1 Show that if we drop the assumption that A is a floating point matrix in Theorem 3.3.1, then Equation 3.3.2 holds with the coefficient “2”replaced by “3.”

P3.3.2 Suppose A is an n-by-n matrix and that $\hat { L }$ and $\hat { U }$ are produced by Algorithm 3.2.1. (a) How many flops are required to compute $\| | \hat { L } | | \hat { U } | \| _ { \infty } ?$ (b) Show $\mathsf { f l } \big ( \big | \hat { L } \big | \big | \hat { U } \big | \big ) \leq ( 1 + 2 n \mathbf { u } ) \big | \hat { L } \big | \big | \hat { U } \big | + O \big ( \mathbf { u } ^ { 2 } \big )$ .

# Notes and References for 3.3

The original roundoff analysis of Gaussian elimination appears in:

J.H. Wilkinson (1961). “Error Analysis of Direct Methods of Matrix Inversion,” J. ACM 8, 281–330.

Various improvements and insights regarding the bounds and have been made over the years, see:

B.A. Chartres and J.C. Geuder (1967). “Computable Error Bounds for Direct Solution of Linear Equations,” J. ACM 14, 63–71.

J.K. Reid (1971). “A Note on the Stability of Gaussian Elimination,” J. Inst. Math. Applic. 8, 374–75.

C.C. Paige (1973). “An Error Analysis of a Method for Solving Matrix Equations,” Math. Comput. 27, 355–59.

H.H. Robertson (1977). “The Accuracy of Error Estimates for Systems of Linear Algebraic Equations,” J. Inst. Math. Applic. 20, 409–14.

J.J. Du Croz and N.J. Higham (1992). “Stability of Methods for Matrix Inversion,” IMA J. Numer. Anal. 12, 1–19.

J.M. Banoczi, N.C. Chiu, G.E. Cho, and I.C.F. Ipsen (1998). “The Lack of Influence of the Right–Hand Side on the Accuracy of Linear System Solution,” SIAM J. Sci. Comput. 20, 203–227.

P. Amodio and F. Mazzia (1999). “A New Approach to Backward Error Analysis of LU Factorization BIT 39, 385–402.

An interesting account of von Neuman’s contributions to the numerical analysis of Gaussian elimination is detailed in:

J.F. Grcar (2011). “John von Neuman’s Analysis of Gaussian Elimination and the Origins of Modern Numerical Analysis,” SIAM Review 53, 607–682.
