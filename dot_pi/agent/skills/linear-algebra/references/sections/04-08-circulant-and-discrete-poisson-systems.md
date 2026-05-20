# 4.8 Circulant and Discrete Poisson Systems

If $A \in \mathbb { C } ^ { n \times n }$ has a factorization of the form

$$
V ^ {- 1} A V = \Lambda = \mathrm{diag} (\lambda_ {1}, \dots , \lambda_ {n}), \tag {4.8.1}
$$

then the columns of V are eigenvectors and the $\lambda _ { i }$ are the corresponding eigenvalues2. In principle, such a decomposition can be used to solve a nonsingular Au = b problem:

$$
u = A ^ {- 1} b = (V \Lambda V ^ {- 1}) ^ {- 1} b = V (\Lambda^ {- 1} (V ^ {- 1} b)). \tag {4.8.2}
$$

However, if this solution framework is to rival the efficiency of Gaussian elimination or the Cholesky factorization, then V and Λ need to be very special. We say that A has a fast eigenvalue decomposition (4.8.1) if

(1) Matrix-vector products of the form $y = V x$ require O(n log n) flops to evaluate.   
(2) The eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { n }$ require O(n log n) flops to evaluate.   
(3) Matrix-vector products of the form $\tilde { b } = V ^ { - 1 } b$ require O(n log n) flops to evaluate.

If these three properties hold, then it follows from (4.8.2) that $O ( n \log n )$ flops are required to solve $A u = b$ .

Circulant systems and related discrete Poisson systems lend themselves to this strategy and are the main concern of this section. In these applications, the V -matrices are associated with the discrete Fourier transform and various sine and cosine transforms. (Now is the time to review §1.4.1 and §1.4.2 and to recall that we have n log n methods for the DFT, DST, DST2, and DCT.) It turns out that fast methods exist for the inverse of these transforms and that is important because of (3). We will not be concerned with precise flop counts because in the fast transform “business”, some n are friendlier than others from the efficiency point of view. While this issue may be important in practice, it is not something that we have to worry about in our brief, proof-of-concept introduction. Our discussion is modeled after §4.3–§4.5 in Van Loan (FFT) where the reader can find complete derivations and greater algorithmic detail. The interconnection between boundary conditions and fast transforms is a central theme and in that regard we also recommend Strang (1999).

# 4.8.1 The Inverse of the DFT Matrix

Recall from §1.4.1 that the DFT matrix $F _ { n } \in \mathbb { C } ^ { n \times n }$ is defined by

$$
[ F _ {n} ] _ {k j} = \omega_ {n} ^ {(k - 1) (j - 1)}, \qquad \omega_ {n} = \cos \left(\frac {2 \pi}{n}\right) - i \sin \left(\frac {2 \pi}{n}\right).
$$

It is easy to verify that

$$
F _ {n} ^ {H} = \bar {F} _ {n}
$$

and so for all p and q that satisfy $0 \leq p < n$ and $0 \leq q < n$ we have

$$
F _ {n} (:, p + 1) ^ {H} F _ {n} (:, q + 1) = \sum_ {k = 0} ^ {n - 1} \bar {\omega} _ {n} ^ {k p} \omega_ {n} ^ {k q} = \sum_ {k = 0} ^ {n - 1} \omega_ {n} ^ {k (q - p)}.
$$

If $q = p ,$ , then this sum equals n. Otherwise,

$$
\sum_ {k = 0} ^ {n - 1} \omega_ {n} ^ {k (q - p)} = \frac {1 - \omega_ {n} ^ {n (q - p)}}{1 - \omega_ {n} ^ {q - p}} = \frac {1 - 1}{1 - \omega_ {n} ^ {q - p}} = 0.
$$

It follows that

$$
n I _ {n} = F _ {n} ^ {H} F _ {n} = \bar {F} _ {n} F _ {n}.
$$

Thus, the DFT matrix is a scaled unitary matrix and

$$
F _ {n} ^ {- 1} = \frac {1}{n} \bar {F} _ {n}.
$$

A fast Fourier transform procedure for $F _ { n } x$ can be turned into a fast inverse Fourier transform procedure for $F _ { n } ^ { - 1 } x$ . Since

$$
y = F _ {n} ^ {- 1} x = \frac {1}{n} \bar {F} _ {n} x,
$$

simply replace each reference to $\omega _ { n }$ with a reference to $\bar { \omega } _ { n }$ and scale. See Algorithm 1.4.1.

# 4.8.2 Circulant Systems

A circulant matrix is a Toeplitz matrix with “wraparound”, e.g.,

$$
C (z) = \left[ \begin{array}{l l l l l} z _ {0} & z _ {4} & z _ {3} & z _ {2} & z _ {1} \\ z _ {1} & z _ {0} & z _ {4} & z _ {3} & z _ {2} \\ z _ {2} & z _ {1} & z _ {0} & z _ {4} & z _ {3} \\ z _ {3} & z _ {2} & z _ {1} & z _ {0} & z _ {4} \\ z _ {4} & z _ {3} & z _ {2} & z _ {1} & z _ {0} \end{array} \right].
$$

We assume that the vector z is complex. Any circulant $C ( z ) \in \mathbb { C } ^ { n \times n }$ is a linear combination of $I _ { n } , D _ { n } , \ldots , D _ { n } ^ { n - 1 }$ where $\mathcal { D } _ { n }$ is the downshift permutation defined in §1.2.11. For example, if $n = 5$ , then

$$
\mathcal {D} _ {5} = \left[ \begin{array}{l l l l l} 0 & 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \end{array} \right]
$$

and

$$
\mathcal {D} _ {5} ^ {2} = \left[ \begin{array}{l l l l l} 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \end{array} \right], \quad \mathcal {D} _ {5} ^ {3} = \left[ \begin{array}{l l l l l} 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 \end{array} \right], \quad \mathcal {D} _ {5} ^ {4} = \left[ \begin{array}{l l l l l} 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 & 0 \end{array} \right].
$$

Thus, the 5-by-5 circulant matrix displayed above is given by

$$
C (z) = z _ {0} I + z _ {1} \mathcal {D} _ {n} + z _ {2} \mathcal {D} _ {n} ^ {2} + z _ {3} \mathcal {D} _ {n} ^ {3} + z _ {4} \mathcal {D} _ {n} ^ {4}.
$$

Note that $\mathcal { D } _ { 5 } ^ { 5 } = I _ { 5 }$ . More generally,

$$
z = \left[ \begin{array}{c} z _ {0} \\ z _ {1} \\ \vdots \\ z _ {n - 1} \end{array} \right] \quad \Rightarrow \quad C (z) = \sum_ {k = 0} ^ {n - 1} z _ {k} \mathcal {D} _ {n} ^ {k}. \tag {4.8.3}
$$

Note that if $V ^ { - 1 } { \mathcal { D } } _ { n } V = \Lambda$ is diagonal, then

$$
V ^ {- 1} C (z) V = V ^ {- 1} \left(\sum_ {k = 0} ^ {n - 1} z _ {k} \mathcal {D} _ {n} ^ {k}\right) V = \sum_ {k = 0} ^ {n - 1} z _ {k} \left(V ^ {- 1} \mathcal {D} _ {n} V ^ {- 1}\right) ^ {k} = \sum_ {k = 0} ^ {n - 1} z _ {k} \Lambda^ {k} \tag {4.8.4}
$$

is diagonal. It turns out that the DFT matrix diagonalizes the downshift permutation.

Lemma 4.8.1. $I f V = F _ { n }$ , then $V ^ { - 1 } { \mathcal D } _ { n } V = \Lambda = \mathrm { d i a g } ( \lambda _ { 1 } , \dots , \lambda _ { n } )$ where

$$
\lambda_ {j + 1} = \bar {\omega} _ {n} ^ {j} = \cos \left(\frac {2 j \pi}{n}\right) + i \sin \left(\frac {2 j \pi}{n}\right)
$$

for $j = 0 { : } n - 1$ .

Proof. For $j = 0 { : } n - 1$ we have

$$
\mathcal {D} _ {n} F _ {n} (:, j + 1) = \mathcal {D} _ {n} \left[ \begin{array}{c} 1 \\ \omega_ {n} ^ {j} \\ \omega_ {n} ^ {2 j} \\ \vdots \\ \omega_ {n} ^ {(n - 1) j} \end{array} \right] = \left[ \begin{array}{c} \omega_ {n} ^ {(n - 1) j} \\ 1 \\ \omega_ {n} ^ {j} \\ \vdots \\ \omega_ {n} ^ {(n - 2) j} \end{array} \right] = \bar {\omega} _ {n} ^ {j} \left[ \begin{array}{c} 1 \\ \omega_ {n} ^ {j} \\ \omega_ {n} ^ {2 j} \\ \vdots \\ \omega_ {n} ^ {(n - 1) j} \end{array} \right].
$$

This vector is precisely $F _ { n } \Lambda ( : , j + 1 )$ . Thus, $\mathcal { D } _ { n } V = V \Lambda$ , i.e., $V ^ { - 1 } { \mathcal { D } } _ { n } V = \Lambda$ .

It follows from (4.8.4) that any circulant $C ( z )$ is diagonalized by $F _ { n }$ and the eigenvalues of $C ( z )$ can be computed fast.

Theorem 4.8.2. Suppose $z \in \mathbb { C } ^ { n }$ and C(z) are defined by (4.8.3). $I f V = F _ { n }$ and $\lambda \ : = \ : \bar { F } _ { n } z$ , then $V ^ { - 1 } C ( z ) V = \mathrm { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ .

Proof. Define

$$
f = \left[ \begin{array}{c} 1 \\ \bar {\omega} _ {n} \\ \vdots \\ \bar {\omega} _ {n} ^ {n - 1} \end{array} \right]
$$

and note that the columns of ${ \bar { F } } _ { n }$ are componentwise powers of this vector. In particular, $\bar { F } _ { n } ( : , k + 1 ) = f . \hat { \cdot } k$ where $[ f . \hat { \cdot } k ] _ { j } = f _ { j } ^ { k }$ . Since $\Lambda = \operatorname { d i a g } ( f )$ , it follows from Lemma 4.8.1 that

$$
\begin{array}{l} V ^ {- 1} C (z) V = \sum_ {k = 0} ^ {n - 1} z _ {k} \Lambda^ {k} = \sum_ {k = 0} ^ {n - 1} z _ {k} \operatorname{diag} (f) ^ {k} = \sum_ {k = 0} ^ {n - 1} z _ {k} \operatorname{diag} (f. \hat {} k) \\ = \operatorname{diag} \left(\sum_ {k = 0} ^ {n - 1} z _ {k} f. \hat {k}\right) = \operatorname{diag} \left(\bar {F} _ {n} z\right) \\ \end{array}
$$

completing the proof of the theorem

Thus, the eigenvalues of the circulant matrix $C ( z )$ are the components of the vector $\bar { F } _ { n } z$ . Using this result we obtain the following algorithm.

Algorithm 4.8.1 If $z \in \mathbb { C } ^ { n } , \quad y \in \mathbb { C } ^ { n }$ , and $C ( z )$ is nonsingular, then the following algorithm solves the linear system $C ( z ) x = y$ .

Use an FFT to compute $c = { \bar { F } } _ { n } y$ and $d = { \bar { F } } _ { n } z$ .

$$
w = c. / d
$$

Use an FFT to compute $u \ = \ F _ { n } w$

$$
x = u / n
$$

This algorithm requires $O ( n \log n )$ flops.

# 4.8.3 The Discretized Poisson Equation in One Dimension

We now turn our attention to a family of real matrices that have real, fast eigenvalue decompositions. The starting point in the discussion is the differential equation

$$
\frac {d ^ {2} u}{d x ^ {2}} = - f (x) \quad \alpha \leq u (x) \leq \beta , \tag {4.8.5}
$$

together with one of four possible specifications of $u ( x )$ on the boundary.

Dirichlet-Dirichlet (DD): $u ( \alpha ) = u _ { \alpha } , \qquad u ( \beta ) = u _ { \beta } ,$

Dirichlet-Neumann (DN): $u ( \alpha ) = u _ { \alpha } , \qquad u ^ { \prime } ( \beta ) = u _ { \beta } ^ { \prime } ,$

Neumann-Neumann (NN): $u ^ { \prime } ( \alpha ) = u _ { \alpha } ^ { \prime } , \qquad u ^ { \prime } ( \beta ) = u _ { \beta } ^ { \prime } .$

Periodic (P): $u ( \alpha ) = u ( \beta ) .$ .

By replacing the derivatives in (4.8.5) with divided differences, we obtain a system of linear equations. Indeed, if m is a positive integer and

$$
h = \frac {\beta - \alpha}{m},
$$

then for $i = 1 { : } m - 1$ we have

$$
\frac {\frac {u _ {i + 1} - u _ {i}}{h} - \frac {u _ {i} - u _ {i - 1}}{h}}{h} = \frac {u _ {i - 1} - 2 u _ {i} + u _ {i + 1}}{h ^ {2}} = - f _ {i} \tag {4.8.6}
$$

where $f _ { i } = f ( \alpha { + } i h )$ and $u _ { i } \approx u ( \alpha { + } i h )$ . To appreciate this discretization we display the linear equations that result when $m = 5$ for the various possible boundary conditions. The matrices $\mathcal { T } _ { n } ^ { ( { D D } ) } , \mathcal { T } _ { n } ^ { ( { D N } ) } , \mathcal { T } _ { n } ^ { ( { N N } ) }$ , and ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ are formally defined afterwards.

For the Dirichlet-Dirichlet problem, the system is 4-by-4 and tridiagonal:

$$
\mathcal {T} _ {4} ^ {(D D)} \cdot u (1: 4) \equiv \left[ \begin{array}{r r r r} 2 & - 1 & 0 & 0 \\ - 1 & 2 & - 1 & 0 \\ 0 & - 1 & 2 & - 1 \\ 0 & 0 & - 1 & 2 \end{array} \right] \left[ \begin{array}{l} u _ {1} \\ u _ {2} \\ u _ {3} \\ u _ {4} \end{array} \right] = \left[ \begin{array}{l} h ^ {2} f _ {1} + u _ {\alpha} \\ h ^ {2} f _ {2} \\ h ^ {2} f _ {3} \\ h ^ {2} f _ {4} + u _ {\beta} \end{array} \right].
$$

For the Dirichlet-Neumann problem the system is still tridiagonal, but $u _ { 5 }$ joins $u _ { 1 } , \ldots , u _ { 4 }$ as an unknown:

$$
\mathcal {T} _ {5} ^ {(D N)} \cdot u (1: 5) \equiv \left[ \begin{array}{r r r r r} 2 & - 1 & 0 & 0 & 0 \\ - 1 & 2 & - 1 & 0 & 0 \\ 0 & - 1 & 2 & - 1 & 0 \\ 0 & 0 & - 1 & 2 & - 1 \\ 0 & 0 & 0 & - 2 & 2 \end{array} \right] \left[ \begin{array}{l} u _ {1} \\ u _ {2} \\ u _ {3} \\ u _ {4} \\ u _ {5} \end{array} \right] = \left[ \begin{array}{c} h ^ {2} f _ {1} + u _ {\alpha} \\ h ^ {2} f _ {2} \\ h ^ {2} f _ {3} \\ h ^ {2} f _ {4} \\ 2 h u _ {\beta} ^ {\prime} \end{array} \right].
$$

The new equation on the bottom is derived from the approximation $u ^ { \prime } ( \beta ) \approx ( u _ { 5 } - u _ { 4 } ) / h$ . (The scaling of this equation by 2 simplifies some of the derivations below.) For the Neumann-Neumann problem, $u _ { 5 }$ and $u _ { 0 }$ need to be determined:

$$
\mathcal {T} _ {6} ^ {(N N)} \cdot u (0 {:} 5) \equiv \left[ \begin{array}{r r r r r r} 2 & - 2 & 0 & 0 & 0 & 0 \\ - 1 & 2 & - 1 & 0 & 0 & 0 \\ 0 & - 1 & 2 & - 1 & 0 & 0 \\ 0 & 0 & - 1 & 2 & - 1 & 0 \\ 0 & 0 & 0 & - 1 & 2 & - 1 \\ 0 & 0 & 0 & 0 & - 2 & 2 \end{array} \right] \left[ \begin{array}{l} u _ {0} \\ u _ {1} \\ u _ {2} \\ u _ {3} \\ u _ {4} \\ u _ {5} \end{array} \right] = \left[ \begin{array}{l} - 2 h u _ {\alpha} ^ {\prime} \\ h ^ {2} f _ {1} \\ h ^ {2} f _ {2} \\ h ^ {2} f _ {3} \\ h ^ {2} f _ {3} \\ 2 h u _ {\beta} ^ {\prime} \end{array} \right].
$$

Finally, for the periodic problem we have

$$
\mathcal {T} _ {5} ^ {(P)} \cdot u (1: 5) \equiv \left[ \begin{array}{r r r r r} 2 & - 1 & 0 & 0 & - 1 \\ - 1 & 2 & - 1 & 0 & 0 \\ 0 & - 1 & 2 & - 1 & 0 \\ 0 & 0 & - 1 & 2 & - 1 \\ - 1 & 0 & 0 & - 1 & 2 \end{array} \right] \left[ \begin{array}{l} u _ {1} \\ u _ {2} \\ u _ {3} \\ u _ {4} \\ u _ {5} \end{array} \right] = \left[ \begin{array}{l} h ^ {2} f _ {1} \\ h ^ {2} f _ {2} \\ h ^ {2} f _ {3} \\ h ^ {2} f _ {4} \\ h ^ {2} f _ {5} \end{array} \right].
$$

The first and last equations use the conditions $u _ { 0 } = u _ { 5 }$ and $u _ { 1 } = u _ { 6 }$ . These constraints follow from the assumption that u has period $\beta - \alpha$ .

As we show below, the n-by-n matrix

$$
\mathcal {T} _ {n} ^ {(D D)} = \left[ \begin{array}{c c c c} 2 & - 1 & \dots & 0 \\ - 1 & 2 & \ddots & \vdots \\ \vdots & \ddots & \ddots & - 1 \\ 0 & \dots & - 1 & 2 \end{array} \right] \tag {4.8.7}
$$

and its low-rank adjustments

$$
\mathcal {T} _ {n} ^ {(D N)} = \mathcal {T} _ {n} ^ {(D D)} - e _ {n} e _ {n - 1} ^ {T}, \tag {4.8.8}
$$

$$
\mathcal {T} _ {n} ^ {(N N)} = \mathcal {T} _ {n} ^ {(D D)} - e _ {n} e _ {n - 1} ^ {T} - e _ {1} e _ {2} ^ {T}, \tag {4.8.9}
$$

$$
\mathcal {T} _ {n} ^ {(P)} = \mathcal {T} _ {n} ^ {(D D)} - e _ {1} e _ {n} ^ {T} - e _ {n} e _ {1} ^ {T}. \tag {4.8.10}
$$

have fast eigenvalue decompositions. However, the existence of O(n log n) methods for these systems is not very interesting because algorithms based on Gaussian elimination are faster: $O ( n )$ versus O(n log n). Things get much more interesting when we discretize the 2-dimensional analogue of (4.8.5).

# 4.8.4 The Discretized Poisson Equation in Two Dimensions

To launch the 2D discussion, suppose $F ( x , y )$ is defined on the rectangle

$$
R = \{(x, y): \alpha_ {x} \leq x \leq \beta_ {x}, \alpha_ {y} \leq y \leq \beta_ {y} \}
$$

and that we wish to find a function u that satisfies

$$
\frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}} = - F (x, y) \tag {4.8.11}
$$

on R and has its value prescribed on the boundary of R. This is Poisson’s equation with Dirichlet boundary conditions. Our plan is to approximate u at the grid points $( \alpha _ { x } + i h _ { x } , \alpha _ { y } + j h _ { y } )$ where $i = 1 { : } m _ { 1 } - 1 , j = 1 { : } m _ { 2 } - 1$ , and

$$
h _ {x} = \frac {\beta_ {x} - \alpha_ {x}}{m _ {1}} \qquad h _ {y} = \frac {\beta_ {y} - \alpha_ {y}}{m _ {2}}.
$$

Refer to Figure 4.8.1, which displays the case when $m _ { 1 } = 6$ and $m _ { 2 } = 5$ . Notice that there are two kinds of grid points. The function u is known at the $\mathbf { \widetilde { \Gamma } } ^ { ( 6 } \bullet ^ { 9 }$ grid points on the boundary. The function u is to be determined at the “◦” grid points in the interior. The interior grid points have been indexed in a top-to-bottom, left-to-right order. The idea is to have $u _ { k }$ approximate the value of $u ( x , y )$ at grid point k.

As in the one-dimensional problem considered §4.8.3, we use divided differences to obtain a set of linear equations that define the unknowns. An interior grid point P has a north (N ), east (E), south (S), and west (W ) neighbor. Using this “compass point” notation we obtain the following approximation to (4.8.11) at $P 3$

$$
\begin{array}{r l r} \frac {u (E) - u (P)}{h _ {x}} - \frac {u (P) - u (W)}{h _ {x}} & + & \frac {u (N) - u (P)}{h _ {y}} - \frac {u (P) - u (S)}{h _ {y}} \\ \hline h _ {x} & + & \frac {h _ {y}}{h _ {y}} = - F (P) \end{array}
$$

![](images/golub_200_249__541d3ad665f3a4144be5b650e770d934fb9ac0e547bc33cba358850ba3e872af.jpg)

<details>
<summary>text_image</summary>

1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
</details>

Figure 4.8.1. A grid with $m _ { 1 } = 6$ and $m _ { 2 } = 5$ .

The x-partial and y-partial have been replaced by second-order divided differences. Assume for clarity that the horizontal and vertical grid spacings are equal, i.e., $h _ { x } =$ $h _ { y } = h$ . With this assumption, the linear equation at point $P$ has the form

$$
4 u (P) - u (N) - u (E) - u (S) - u (W) = h ^ {2} F (P).
$$

In our example, there are 20 such equations. It should be noted that some of $P ^ { * } { \mathrm { s } }$ neighbors may be on the boundary, in which case the corresponding linear equation involves fewer than 5 unknowns. For example, if P is the third grid point then we see from Figure 4.8.1 that the north neighbor N is on the boundary. It follows that the associated linear equation has the form

$$
4 u (P) - u (E) - u (S) - u (W) = h ^ {2} F (P) + u (N).
$$

Reasoning like this, we conclude that the matrix of coefficients has the following block tridiagonal form

$$
A = \left[ \begin{array}{c c c c} \mathcal {T} _ {5} ^ {(D D)} & 0 & 0 & 0 \\ 0 & \mathcal {T} _ {5} ^ {(D D)} & 0 & 0 \\ 0 & 0 & \mathcal {T} _ {5} ^ {(D D)} & 0 \\ 0 & 0 & 0 & \mathcal {T} _ {5} ^ {(D D)} \end{array} \right] + \left[ \begin{array}{c c c c} 2 I _ {5} & - I _ {5} & 0 & 0 \\ - I _ {5} & 2 I _ {5} & - I _ {5} & 0 \\ 0 & - I _ {5} & 2 I _ {5} & - I _ {5} \\ 0 & 0 & - I _ {5} & 2 I _ {5} \end{array} \right]
$$

i.e.,

$$
A = I _ {4} \otimes \mathcal {T} _ {5} ^ {(D D)} + \mathcal {T} _ {4} ^ {(D D)} \otimes I _ {5}.
$$

Notice that the first matrix is associated with the x-partials while the second matrix is associated with the y-partials. The right-hand side in $A u = b$ is made up of $F _ { - }$ evaluations and specified values of $u ( x , y )$ on the boundary.


---

<!-- golub_250_299 -->

Extrapolating from our example, we conclude that the matrix of coefficients is an $( m _ { 2 } - 1 ) – \mathrm { b y } – ( m _ { 2 } - 1 )$ block tridiagonal matrix with $( m _ { 1 } - 1 ) – \mathrm { b y } – ( m _ { 1 } - 1 )$ blocks:

$$
A = I _ {m _ {2} - 1} \otimes \mathcal {T} _ {m _ {1} - 1} ^ {(D D)} + \mathcal {T} _ {m _ {2} - 1} ^ {(D D)} \otimes I _ {m _ {1} - 1}.
$$

Alternative specifications along the boundary lead to systems with similar structure, $\mathrm { e . g . }$ ,

$$
A u \equiv \left(I _ {n _ {2}} \otimes A _ {1} + A _ {2} \otimes I _ {n _ {1}}\right) u = b. \tag {4.8.12}
$$

For example, if we impose Dirichet-Neumann, Neumann-Neumann, or periodic boundary coequal the, or and right edges of the rectangular domain accordingly. Likewise, if we impose Dir $R ,$ then et-Ne $A _ { 1 }$ willann, $\mathscr { T } _ { m _ { 1 } } ^ { ( D N ) } , \mathscr { T } _ { m _ { 1 } + 1 } ^ { ( N N ) }$ $\mathcal { T } _ { m _ { 1 } } ^ { \left( P \right) }$ Neumann-Neumann, or periodic boundary conditions along the bottom and top edges of R, then A2 will equal T (DNm2 $R ,$ $A _ { 2 }$ $\mathcal { T } _ { m _ { 2 } } ^ { ( D N ) } , \mathcal { T } _ { m _ { 2 } + 1 } ^ { ( N N ) }$ , or $\mathcal { T } _ { m _ { 2 } } ^ { \left( P \right) }$ . If the system (4.8.12) is nonsingular and $A _ { 1 }$ and $A _ { 2 }$ have fast eigenvalue decompositions, then it can be solved with just O(N log N) flops where $N = n _ { 1 } n _ { 2 }$ . To see why this is possible, assume that

$$
V ^ {- 1} A _ {1} V = D _ {1} = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n _ {1}}\right), \tag {4.8.13}
$$

$$
W ^ {- 1} A _ {2} W = D _ {2} = \operatorname{diag} \left(\mu_ {1}, \dots , \mu_ {n _ {2}}\right) \tag {4.8.14}
$$

are fast eigenvalue decompositions. Using facts about the Kronecker product that are set forth in §1.3.6–§1.3.8, we can reformulate (4.8.12) as a matrix equation

$$
A _ {1} U + U A _ {2} ^ {T} = B
$$

where $U = { \mathsf { r e s h a p e } } ( u , n _ { 1 } , n _ { 2 } )$ and $\boldsymbol { B } = \mathsf { r e s h a p e } ( b , n _ { 1 } , n _ { 2 } )$ . Substituting the above eigenvalue decompositions into this equation we obtain

$$
D _ {1} \tilde {U} + \tilde {U} D _ {2} = \tilde {B},
$$

where $\tilde { U } = ( \tilde { u } _ { i j } ) = V ^ { - 1 } U W ^ { - T }$ and $\begin{array} { r } { \tilde { B } = ( \tilde { b } _ { i j } ) = V ^ { - 1 } B W ^ { - T } } \end{array}$ . Note how easy it is to solve this transformed system because $D _ { 1 }$ and $D _ { 2 }$ are diagonal:

$$
\tilde {u} _ {i j} = \frac {\tilde {b} _ {i j}}{\lambda_ {i} + \mu_ {j}} \qquad i = 1: n _ {1}, j = 1: n _ {2}.
$$

For this to be well-defined, no eigenvalue of $A _ { 1 }$ can be the negative of an eigenvalue of $A _ { 2 }$ . In our example, all the $\lambda _ { i }$ and $\mu _ { i }$ are positive. Overall we obtain

Algorithm 4.8.2 (Fast Poisson Solver Framework) Assume that $A _ { 1 } \in \mathbb { R } ^ { n _ { 1 } \times n _ { 1 } }$ and $A _ { 2 } \in \mathbb { R } ^ { n _ { 2 } \times n _ { 2 } }$ have fast eigenvalue decompositions (4.8.13) and (4.8.14) and that the matrix $A = I _ { n _ { 2 } } \otimes A _ { 1 } + A _ { 2 } \otimes I _ { n _ { 1 } }$ is nonsingular. The following algorithm solves the linear system Au = b where b ∈ IRn1n2 . $A u = b$ $b \in \mathbb { R } ^ { n _ { 1 } n _ { 2 } }$

$$
\tilde {B} = (W ^ {- 1} (V ^ {- 1} B) ^ {T}) ^ {T} \text {   where   } B = \operatorname{reshape} (b, n _ {1}, n _ {2})
$$

for $i = 1 { : } n _ { 1 }$

$$
\text { for } j = 1: n _ {2}
$$

$$
\tilde {u} _ {i j} = \tilde {b} _ {i j} / (\lambda_ {i} + \mu_ {j})
$$

end

end

$$
u = \operatorname{reshape} (U, n _ {1} n _ {2}, 1) \text {   where   } U = (W (V \tilde {U}) ^ {T}) ^ {T}
$$

The following table accounts for the work involved:

<table><tr><td>Operation</td><td>How Many?</td><td>Work</td></tr><tr><td> $V^{-1}$  times  $n_1$ -vector</td><td> $n_2$ </td><td> $O(n_2 \cdot n_1 \cdot \log n_1)$ </td></tr><tr><td> $W^{-1}$  times  $n_2$ -vector</td><td> $n_1$ </td><td> $O(n_1 \cdot n_2 \cdot \log n_2)$ </td></tr><tr><td> $V$  times  $n_1$ -vector</td><td> $n_2$ </td><td> $O(n_2 \cdot n_1 \cdot \log n_1)$ </td></tr><tr><td> $W$  times  $n_2$ -vector</td><td> $n_1$ </td><td> $O(n_1 \cdot n_2 \cdot \log n_2)$ </td></tr></table>

Adding up the operation counts, we see that $O ( n _ { 1 } n _ { 2 } \log ( n _ { 1 } n _ { 2 } ) ) = O ( N \log N )$ flops are required where $N = n _ { 1 } n _ { 2 }$ is the size of the matrix A.

Below we show that the matrices $\mathcal { T } _ { n } ^ { ( { D D } ) } , \mathcal { T } _ { n } ^ { ( { D N } ) } , \mathcal { T } _ { n } ^ { ( { N N } ) }$ , and ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ have fast eigenvalue decompositions and this means that Algorithm 4.8.2 can be used to solve discrete Poisson systems. To appreciate the speedup over conventional methods, suppose $A _ { 1 } = \mathcal { T } _ { n _ { 1 } } ^ { \left( D D \right) }$ and $A _ { 2 } = \mathcal { T } _ { n _ { 2 } } ^ { \left( D D \right) }$ . It can be shown that A is symmetric positive definite with bandwidth $n _ { 1 } + 1$ . Solving $A u = b$ using Algorithm 4.3.5 (band Cholesky) would require $O ( n _ { 1 } ^ { 3 } n _ { 2 } ) = O ( N n _ { 1 } ^ { 2 } )$ flops.

# 4.8.5 The Inverse of the DST and DCT Matrices

The eigenvector matrices for $\mathcal { T } _ { n } ^ { ( { D D } ) } , \mathcal { T } _ { n } ^ { ( { D N } ) } , \mathcal { T } _ { n } ^ { ( { N N } ) }$ , and T (P ) n ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ are associated with the fast trigonometric transforms presented in §1.4.2. It is incumbent upon us to show that the inverse of these transforms can also be computed fast. We do this for the discrete sine transform (DST) and the discrete cosine transform (DCT) and leave similar fast inverse verifications to the exercises at the end of the section.

By considering the blocks of the DFT matrix $F _ { 2 m }$ , we can determine the inverses of the transform matrices $\mathrm { D S T } ( m - 1 )$ and $\mathrm { D C T } ( m + 1 )$ . Recall from §1.4.2 that if $C _ { r } \in \mathbb { R } ^ { r \times r }$ and $S _ { r } \in \mathbb { R } ^ { r \times r }$ are defined by

$$
[ C _ {r} ] _ {k j} = \cos \left(\frac {k j \pi}{r + 1}\right), [ S _ {r} ] _ {k j} = \sin \left(\frac {k j \pi}{r + 1}\right)
$$

then

$$
F _ {2 m} = \left[ \begin{array}{c c c c} 1 & e ^ {T} & 1 & e ^ {T} \\ e & C - i S & v & (C + i S) E \\ 1 & v ^ {T} & (- 1) ^ {m} & v ^ {T} E \\ e & E (C + i S) & E v & E (C - i S) E \end{array} \right]
$$

where $C = C _ { m - 1 } , S = S _ { m - 1 } , E = \mathcal { E } _ { m - 1 }$ , and

$$
e ^ {T} = (\underbrace {1 , 1 , \ldots , 1} _ {m - 1}) \qquad v ^ {T} = (\underbrace {- 1 , 1 , \ldots , (- 1) ^ {m - 1}} _ {m - 1}).
$$

By comparing the (2,1), (2,2), (2,3), and (2,4) blocks in the equation $2 m I = \bar { F } _ { 2 m } F _ { 2 m }$ we conclude that

$$
0 = 2 C e + e + v,
$$

$$
\begin{array}{l} 2 m I _ {m - 1} = 2 C ^ {2} + 2 S ^ {2} + e e ^ {T} + v v ^ {T}, \\ 0 = 2 C v + e + (- 1) ^ {m} v, \\ 0 = 2 C ^ {2} - 2 S ^ {2} + e e ^ {T} + v v ^ {T}. \\ \end{array}
$$

It follows that $2 S ^ { 2 } = m I _ { m - 1 }$ and $2 C ^ { 2 } = m I _ { m - 1 } - e e ^ { T } - v v ^ { T }$ . Using these equations it is easy to verify that

$$
S _ {m - 1} ^ {- 1} = \frac {2}{m} S _ {m - 1}
$$

and

$$
\left[ \begin{array}{c c c} 1 / 2 & e ^ {T} & 1 / 2 \\ e / 2 & C _ {m - 1} & v / 2 \\ 1 / 2 & v ^ {T} & (- 1) ^ {m} / 2 \end{array} \right] ^ {- 1} = \frac {2}{m} \left[ \begin{array}{c c c} 1 / 2 & e ^ {T} & 1 / 2 \\ e / 2 & C _ {m - 1} & v / 2 \\ 1 / 2 & v ^ {T} & (- 1) ^ {m} / 2 \end{array} \right].
$$

Thus, it follows from the definitions (1.4.8) and (1.4.10) that

$$
V = \mathrm{DST} (m - 1) \Rightarrow V ^ {- 1} = \frac {2}{m} \mathrm{DST} (m - 1),
$$

$$
V = \mathrm{DCT} (m + 1) \Rightarrow V ^ {- 1} = \frac {2}{m} \mathrm{DCT} (m + 1).
$$

In both cases, the inverse transform is a multiple of the “forward” transform and can be computed fast. See Algorithms 1.4.2 and 1.4.3.

# 4.8.6 Four Fast Eigenvalue Decompositions

The matrices $\mathcal { T } _ { n } ^ { ( { D D } ) } , \mathcal { T } _ { n } ^ { ( { D N } ) } , \mathcal { T } _ { n } ^ { ( { N N } ) }$ , and ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ do special things to vectors of sines and cosines.

Lemma 4.8.3. Define the real n-vectors s(θ) and $c ( \theta )$ by

$$
s (\theta) = \left[ \begin{array}{c} s _ {1} \\ \vdots \\ s _ {n} \end{array} \right], \quad c (\theta) = \left[ \begin{array}{c} c _ {0} \\ \vdots \\ c _ {n - 1} \end{array} \right], \tag {4.8.15}
$$

where $s _ { k } = \sin ( k \theta )$ and $c _ { k } = \cos ( k \theta )$ . $\textit { I f e } _ { k } = \textit { I } _ { n } ( : , k )$ and $\lambda = 4 \sin ^ { 2 } ( \theta / 2 )$ , then

$$
\mathcal {T} _ {n} ^ {(D D)} \cdot s (\theta) = \lambda \cdot s (\theta) + s _ {n + 1} e _ {n}, \tag {4.8.16}
$$

$$
\mathcal {T} _ {n} ^ {(D D)} \cdot c (\theta) = \lambda \cdot c (\theta) + c _ {1} e _ {1} + c _ {n} e _ {n}, \tag {4.8.17}
$$

$$
\mathcal {T} _ {n} ^ {(D N)} \cdot s (\theta) = \lambda \cdot s (\theta) + (s _ {n + 1} - s _ {n - 1}) e _ {n}, \tag {4.8.18}
$$

$$
\mathcal {T} _ {n} ^ {(N N)} \cdot c (\theta) = \lambda \cdot c (\theta) + (c _ {n} - c _ {n - 2}) e _ {n}, \tag {4.8.19}
$$

$$
\mathcal {T} _ {n} ^ {(P)} \cdot s (\theta) = \lambda \cdot s (\theta) - s _ {n} e _ {1} + (s _ {n + 1} - s _ {1}) e _ {n}, \tag {4.8.20}
$$

$$
\mathcal {T} _ {n} ^ {(P)} \cdot c (\theta) = \lambda \cdot c (\theta) + (c _ {1} - c _ {n - 1}) e _ {1} + (c _ {n} - 1) e _ {n}. \tag {4.8.21}
$$

Proof. The proof is mainly an exercise in using the trigonometric identities

$$
s _ {k - 1} = c _ {1} s _ {k} - s _ {1} c _ {k}, \quad c _ {k - 1} = c _ {1} c _ {k} + s _ {1} s _ {k},
$$

$$
s _ {k + 1} = c _ {1} s _ {k} + s _ {1} c _ {k}, \quad c _ {k + 1} = c _ {1} c _ {k} - s _ {1} s _ {k}.
$$

For example, if $y = { \mathcal { T } _ { n } ^ { \left( D D \right) } } s ( \theta )$ , then

$$
y _ {k} = \left\{ \begin{array}{l l} 2 s _ {1} - s _ {2} = 2 s _ {1} (1 - c _ {1}), & \text { if } k = 1, \\ - s _ {k - 1} + 2 s _ {k} - s _ {k + 1} = 2 s _ {k} (1 - c _ {1}), & \text { if } 2 \leq k \leq n - 1, \\ - s _ {n - 1} + 2 s _ {n} = 2 s _ {n} (1 - c _ {1}) + s _ {n + 1}, & \text { if } k = n. \end{array} \right.
$$

Equation (4.8.16) follows since $( 1 - c _ { 1 } ) = 1 - \cos ( \theta ) = 2 \sin ^ { 2 } ( \theta / 2 )$ . The proof of (4.8.17) is similar while the remaining equations follow from Equations (4.8.8)–(4.8.10).

Notice that (4.8.16)-(4.8.21) are eigenvector equations except for the $" e _ { 1 } \ "$ and $" e _ { n } "$ terms. By choosing the right value for θ, we can make these residuals disappear, thereby obtaining recipes for the eigensystems of $\mathcal { T } _ { n } ^ { ( { D D } ) } , \mathcal { T } _ { n } ^ { ( { D N } ) } , \mathcal { T } _ { n } ^ { ( { N N } ) }$ T n N) , and T (P ) n . ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$

The Dirichlet-Dirichlet Matrix

If j is an integer and $\theta = j \pi / ( n + 1 )$ , then $s _ { n + 1 } = \sin ( ( n + 1 ) \theta ) = 0$ . It follows from (4.8.16) that

$$
\mathcal {T} _ {n} ^ {(D D)} s (\theta_ {j}) = 4 \sin^ {2} (\theta_ {j} / 2) s (\theta_ {j}), \qquad \theta_ {j} = \frac {j \pi}{n + 1},
$$

for $j = 1 { : } n$ . Thus, the columns of the matrix $V _ { n } ^ { \left( D D \right) } \in \mathbb { R } ^ { n \times n }$ defined by

$$
[ V _ {n} ^ {(D D)} ] _ {k j} = \sin \left(\frac {k j \pi}{n + 1}\right)
$$

are eigenvectors for T (DD) n ${ \mathcal { T } } _ { n } ^ { ( D D ) }$ and the corresponding eigenvalues are given by

$$
\lambda_ {j} = 4 \sin^ {2} \left(\frac {j \pi}{2 (n + 1)}\right),
$$

for $j = 1 { : } n$ . Note that $V _ { n } ^ { \left( { D D } \right) } = \mathrm { D S T } ( n )$ . It follows that ${ \mathcal { T } } _ { n } ^ { ( D D ) }$ has a fast eigenvalue decomposition.

The Dirichlet-Neumann Matrix

If j is an integer and $\theta = ( 2 j - 1 ) \pi / ( 2 n )$ , then $s _ { n + 1 } - s _ { n - 1 } \ = \ 2 s _ { 1 } c _ { n } = 0$ . It follows from (4.8.18) that

$$
\mathcal {T} _ {n} ^ {(D N)} \cdot s (\theta_ {j}) = 4 \sin^ {2} (\theta_ {j} / 2) \cdot s (\theta_ {j}), \qquad \theta_ {j} = \frac {(2 j - 1) \pi}{2 n},
$$

for j = 1:n. Thus, the columns of the matrix V (DN) n $j = 1 { : } n$ $V _ { n } ^ { \left( D N \right) } \in \mathbb { R } ^ { n \times n }$ defined by

$$
[ V _ {n} ^ {(D N)} ] _ {k j} = \sin \left(\frac {k (2 j - 1) \pi}{2 n}\right)
$$

are eigenvectors of the matrix ${ \mathcal { T } } _ { n } ^ { ( D N ) }$ and the corresponding eigenvalues are given by

$$
\lambda_ {j} = 4 \sin^ {2} \left(\frac {(2 j - 1) \pi}{4 n}\right)
$$

for $j = 1 { : } n$ . Comparing with (1.4.13) we see that that $V _ { n } ^ { ( D N ) } = { \mathrm { D S T 2 } } ( n )$ . The inverse DST2 can be evaluated fast. See Van Loan (FFT, p. 242) for details, but also P4.8.11. It follows that $\mathcal { T } ^ { ( D N ) }$ has a fast eigenvalue decomposition.

The Neumann-Neumann Matrix

If j is an integer and $\theta = ( j - 1 ) \pi / ( n - 1 )$ , then $c _ { n } - c _ { n - 2 } = - 2 s _ { 1 } s _ { n - 1 } = 0$ . It follows from (4.8.19) that

$$
\mathcal {T} _ {n} ^ {(N N)} \cdot c (\theta_ {j}) = 4 \sin^ {2} \left(\frac {\theta_ {j}}{2}\right) \cdot c (\theta_ {j}), \qquad \theta_ {j} = \frac {(j - 1) \pi}{n - 1}.
$$

Thus, the columns of the matrix $V _ { n } ^ { \left( D N \right) } \in \mathbb { R } ^ { n \times n }$ defined by

$$
[ V _ {n} ^ {(N N)} ] _ {k j} = \cos \left(\frac {(k - 1) (j - 1) \pi}{n - 1}\right)
$$

are eigenvectors of the matrix ${ \mathcal { T } } _ { n } ^ { ( D N ) }$ and the corresponding eigenvalues are given by

$$
\lambda_ {j} = 4 \sin^ {2} \left(\frac {(j - 1) \pi}{2 (n - 1)}\right)
$$

for $j = 1 { : } n$ . Comparing with (1.4.10) we see that

$$
V _ {n} ^ {(N N)} = \mathrm{DCT} (n) \cdot \mathrm{diag} (2, I _ {n - 2}, 2)
$$

and therefore $\mathcal { T } ^ { ( N N ) }$ has a fast eigenvalue decomposition.

The Periodic Matrix

We can proceed to work out the eigenvalue decomposition for ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ as we did in the previous three cases, i.e., by zeroing the residuals in (4.8.20) and (4.8.21). However, $\mathcal { T } _ { n } ^ { \left( P \right) }$ is a circulant matrix and so we know from Theorem 4.8.2 that

$$
F _ {n} ^ {- 1} \mathcal {T} _ {n} ^ {(P)} F _ {n} = \mathrm{diag} (\lambda_ {1}, \ldots , \lambda_ {n})
$$

where

$$
\lambda = \bar {F} _ {n} \left[ \begin{array}{c} 2 \\ - 1 \\ 0 \\ \vdots \\ - 1 \end{array} \right] = 2 \bar {F} _ {n} (:, 1) - \bar {F} _ {n} (:, 2) - \bar {F} _ {n} (:, n).
$$

It can be shown that

$$
\lambda_ {j} = 4 \sin^ {2} \left(\frac {(j - 1) \pi}{n}\right)
$$

for $j = 1 { : } n$ . It follows that ${ \mathcal { T } } _ { n } ^ { \left( P \right) }$ has a fast eigenvalue decomposition. However, since this matrix is real it is preferable to have a real V -matrix. Using the facts that

$$
\lambda_ {j} = \lambda_ {n + 2 - j} \tag {4.8.22}
$$

and

$$
\bar {F} _ {n} (:, j) = F _ {n} (:, (n + 2 - j)) \tag {4.8.23}
$$

for $j = 2 { : } n$ , it can be shown that if $m = { \mathsf { c e i l } } ( ( n + 1 ) / 2 )$ and

$$
V _ {n} ^ {(P)} = \left[ \operatorname{Re} \left(F _ {n} (:, 1: m) \mid \operatorname{Im} \left(F _ {n} (:, m + 1: n)\right) \right] \right. \tag {4.8.24}
$$

then

$$
\mathcal {T} _ {n} ^ {(P)} V _ {n} ^ {(P)} (:, j) = \lambda_ {j} V _ {n} ^ {(P)} (:, j) \tag {4.8.25}
$$

for $j = 1 { : } n$ . Manipulations with this real matrix and its inverse can be carried out rapidly as discussed in Van Loan (FFT, Chap. 4).

# 4.8.7 A Note on Symmetry and Boundary Conditions

In our presentation, the matrices ${ \mathcal { T } } _ { n } ^ { ( D N ) }$ ) and T (NN) n a ${ \mathcal { T } } _ { n } ^ { ( N N ) }$ re not symmetric. However, a simple diagonal similarity transformation changes this. For example, if $D = \mathrm { d i a g } ( I _ { n - 1 } , \sqrt { 2 } )$ , then $D ^ { - 1 } \mathcal { T } _ { n } ^ { ( D N ) } D$ is symmetric. Working with symmetric second difference matrices has certain attractions, i.e., the automatic orthogonality of the eigenvector matrix. See Strang (1999).

# Problems

P4.8.1 Suppose $z \in \mathbb { R } ^ { n }$ has the property that $z ( 2 : n ) = \mathcal { E } _ { n - 1 } z ( 2 : n )$ . Show that $C ( z )$ is symmetric and $\bar { F } _ { n } z$ is real.

P4.8.2 As measured in the Frobenius norm, what is the nearest real circulant matrix to a given real Toeplitz matrix?

P4.8.3 Given $x , z \in \mathbb { C } ^ { n }$ , show how to compute $y = C ( z ) \cdot x$ in O(n log n) flops. In this case, y is the cyclic convolution of x and z.

P4.8.4 Suppose $\begin{array} { l } { { a } } \end{array} = \left[ \begin{array} { l } { { a } _ { - n + 1 } , \ldots , { a } _ { - 1 } , { a } _ { 0 } , { a } _ { 1 } , \ldots , { a } _ { n - 1 } } \end{array} \right]$ and let $T = ( t _ { k j } )$ be the n-by-n Toeplitz matrix defined by $t _ { k j } = a _ { k - j }$ . Thus, if $\boldsymbol { a } = \left[ a _ { - 2 } , a _ { - 1 } , a _ { 0 } , a _ { 1 } , a _ { 2 } \right]$ , then

$$
T = T (a) = \left[ \begin{array}{c c c} a _ {0} & a _ {- 1} & a _ {- 2} \\ a _ {1} & a _ {0} & a _ {- 1} \\ a _ {2} & a _ {1} & a _ {0} \end{array} \right].
$$

It is possible to “embed” T into a circulant, e.g.,

$$
C = \left[ \begin{array}{c c c c c c c c} a _ {0} & a _ {- 1} & a _ {- 2} & 0 & 0 & 0 & a _ {2} & a _ {1} \\ a _ {1} & a _ {0} & a _ {- 1} & a _ {- 2} & 0 & 0 & 0 & a _ {2} \\ a _ {2} & a _ {1} & a _ {0} & a _ {- 1} & a _ {- 2} & 0 & 0 & 0 \\ 0 & a _ {2} & a _ {1} & a _ {0} & a _ {- 1} & a _ {- 2} & 0 & 0 \\ 0 & 0 & a _ {2} & a _ {1} & a _ {0} & a _ {- 1} & a _ {- 2} & 0 \\ 0 & 0 & 0 & a _ {2} & a _ {1} & a _ {0} & a _ {- 1} & a _ {- 2} \\ a _ {- 2} & 0 & 0 & 0 & a _ {2} & a _ {1} & a _ {0} & a _ {- 1} \\ a _ {- 1} & a _ {- 2} & 0 & 0 & 0 & a _ {2} & a _ {1} & a _ {0} \end{array} \right].
$$

Given $a _ { - n + 1 } , \ldots , a _ { - 1 } , 1 _ { 0 } , a _ { 1 } , \ldots , a _ { n - 1 }$ and $m \geq 2 n - 1$ , show how to construct a vector $v \in \mathbb { C } ^ { m }$ so that if $C = C ( v )$ , then $C ( 1 { : } n , 1 { : } n ) = T$ . Note that v is not unique if $m > 2 n - 1$ .

P4.8.5 Complete the proof of Lemma 4.8.3.

P4.8.6 Show how to compute a Toeplitz-vector product $y = T u$ in n log n time using the embedding idea outlined in the previous problem and the fact that circulant matrices have a fast eigenvalue decomposition.

P4.8.7 Give a complete specification of the vector b in (4.8.12) if $A _ { 1 } = \mathcal { T } _ { n _ { 1 } } ^ { \left( D D \right) } , A _ { 2 } = \mathcal { T } _ { n _ { 2 } } ^ { \left( D D \right) }$ D) , A2 = T (DD) n2 , and $u ( x , y ) = 0$ on the boundary of the rectangular domain R. In terms of the underlying grid, $n _ { 1 } = m _ { 1 } - 1$ and $n _ { 2 } = m _ { 2 } - 1$ .

P4.8.8 Give a complete specification of the vector b in (4.8.12) if $A _ { 1 } \ = \ T _ { n _ { 1 } } ^ { ( D N ) } , \ A _ { 2 } \ = \ T _ { n _ { 2 } } ^ { ( D N ) }$ , $u ( x , y ) = 0$ on the bottom and left edge of $R , u _ { x } ( x , y ) = 0$ along the right edge of R, and $u _ { y } ( x , y ) = 0$ along the top edge of R. In terms of the underlying grid, $n _ { 1 } = m _ { 1 }$ and $n _ { 2 } = m _ { 2 }$ .

P4.8.9 Define a Neumann-Dirichlet matrix ${ \mathcal { T } } _ { n } ^ { ( N D ) }$ that would arise in conjunction with (4.8.5) if $u ^ { \prime } ( \alpha )$ and $u ( \beta )$ were specified. Show that ${ \mathcal { T } } _ { n } ^ { ( N D ) }$ has a fast eigenvalue decomposition.

${ \mathcal { T } } _ { n } ^ { ( N N ) }$ and  woul ${ \mathcal { T } } _ { n } ^ { ( P ) }$ are singular. (a) Assum solve the linear system b is in the range of subject to the const $A =$ In2 ⊗ T (P ) n1 + $I _ { n _ { 2 } } \otimes \mathcal { T } _ { n _ { 1 } } ^ { ( P ) } + \mathcal { T } _ { n _ { 2 } } ^ { ( P ) } \otimes I _ { n _ { 1 } }$ $A u = b$ Repeat part (a) replacing $\mathcal { T } _ { n _ { 1 } } ^ { \left( P \right) }$ with $\mathcal { T } _ { n _ { 1 } } ^ { ( N N ) }$ and $\mathcal { T } _ { n _ { 2 } } ^ { ( P ) }$ with $\mathcal { T } _ { n _ { 2 } } ^ { ( N N ) }$ .

P4.8.11 Let V be the matrix that defines the $\mathrm { D S T 2 } ( n )$ transformation in (1.4.12). (a) Show that

$$
V ^ {T} V = \frac {n}{2} I _ {n} + \frac {1}{2} v v ^ {T}
$$

where $v = [ 1 , - 1 , 1 , . . . , ( - 1 ) ^ { n } ] ^ { T }$ . (b) Verify that

$$
V ^ {- 1} = \frac {2}{n} \left(I - \frac {1}{2 n} v v ^ {T}\right) V ^ {T}.
$$

(c) Show how to compute $V ^ { - 1 } x$ rapidly.

P4.8.12 Verify (4.8.22), (4.8.23), and (4.8.25).

P4.8.13 Show that if $V = V _ { 2 m } ^ { \left( P \right) }$ V (P ) defined in (4.8.24), then

$$
V ^ {T} V = m \left(I _ {n} + e _ {1} e _ {1} ^ {T} + e _ {m + 1} e _ {m + 1} ^ {T}\right).
$$

What can you say about $V ^ { T } V { \mathrm { ~ i f ~ } } V = V _ { 2 m - 1 } ^ { ( P ) } ?$

# Notes and References for §4.8

As we mentioned, this section is based on Van Loan (FFT). For more details about fast Poisson solvers, see:

R.W. Hockney (1965). “A Fast Direct Solution of Poisson’s Equation Using Fourier Analysis,” J. Assoc. Comput. Mach. 12, 95–113.   
B. Buzbee, G. Golub, and C. Nielson (1970). “On Direct Methods for Solving Poisson’s Equation,” SIAM J. Numer. Anal. 7, 627–656.   
F. Dorr (1970). “The Direct Solution of the Discrete Poisson Equation on a Rectangle,” SIAM Review 12, 248–263.   
R. Sweet (1973). “Direct Methods for the Solution of Poisson’s Equation on a Staggered Grid,” J. Comput. Phys. 12, 422–428.   
P.N. Swarztrauber (1974). “A Direct Method for the Discrete Solution of Separable Elliptic Equations,” SIAM J. Numer. Anal. 11, 1136–1150.   
P.N. Swarztrauber (1977). “The Methods of Cyclic Reduction, Fourier Analysis and Cyclic Reduction-Fourier Analysis for the Discrete Solution of Poisson’s Equation on a Rectangle,” SIAM Review 19, 490–501.

There are actually eight variants of the discrete cosine transform each of which corresponds to the location of the Neumann conditions and how the divided difference approximations are set up. For a unified, matrix-based treatment, see:

G. Strang (1999). “The Discrete Cosine Transform,” SIAM Review 41, 135–147.
