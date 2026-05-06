# 1.4 Fast Matrix-Vector Products

In this section we refine our ability to think at the block level by examining some matrix-vector products $y = A x$ in which the n-by-n matrix A is so highly structured that the computation can be carried out with many fewer than the usual $O ( n ^ { 2 } )$ flops. These results are used in §4.8.

# 1.4.1 The Fast Fourier Transform

The discrete Fourier transform (DFT) of a vector $x \in \mathbb { C } ^ { n }$ is a matrix-vector product

$$
y = F _ {n} x
$$

where the DFT matrix $F _ { n } = ( f _ { k j } ) \in \mathbb { C } ^ { n \times n }$ is defined by

$$
f _ {k j} = \omega_ {n} ^ {(k - 1) (j - 1)} \tag {1.4.1}
$$

with

$$
\omega_ {n} = \exp (- 2 \pi i / n) = \cos (2 \pi / n) - i \cdot \sin (2 \pi / n). \tag {1.4.2}
$$

Here is an example:

$$
F _ {4} = \left[ \begin{array}{c c c c} 1 & 1 & 1 & 1 \\ 1 & \omega_ {4} & \omega_ {4} ^ {2} & \omega_ {4} ^ {3} \\ 1 & \omega_ {4} ^ {2} & \omega_ {4} ^ {4} & \omega_ {4} ^ {6} \\ 1 & \omega_ {4} ^ {3} & \omega_ {4} ^ {6} & \omega_ {4} ^ {9} \end{array} \right] = \left[ \begin{array}{c c c c} 1 & 1 & 1 & 1 \\ 1 & - i & - 1 & i \\ 1 & - 1 & 1 & - 1 \\ 1 & i & - 1 & - i \end{array} \right].
$$

The DFT is ubiquitous throughout computational science and engineering and one reason has to do with the following property:

If n is highly composite, then it is possible to carry out the DFT in many fewer than the $O ( n ^ { 2 } )$ flops required by conventional matrix-vector multiplication.

To illustrate this we set $n = 2 ^ { t }$ and proceed to develop the radix-2 fast Fourier transform.

The starting point is to examine the block structure of an even-order DFT matrix after its columns are reordered so that the odd-indexed columns come first. Consider the case

$$
F _ {8} = \left[ \begin{array}{c c c c c c c c} 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\ 1 & \omega & \omega^ {2} & \omega^ {3} & \omega^ {4} & \omega^ {5} & \omega^ {6} & \omega^ {7} \\ 1 & \omega^ {2} & \omega^ {4} & \omega^ {6} & 1 & \omega^ {2} & \omega^ {4} & \omega^ {6} \\ 1 & \omega^ {3} & \omega^ {6} & \omega & \omega^ {4} & \omega^ {7} & \omega^ {2} & \omega^ {5} \\ 1 & \omega^ {4} & 1 & \omega^ {4} & 1 & \omega^ {4} & 1 & \omega^ {4} \\ 1 & \omega^ {5} & \omega^ {2} & \omega^ {7} & \omega^ {4} & \omega & \omega^ {6} & \omega^ {3} \\ 1 & \omega^ {6} & \omega^ {4} & \omega^ {2} & 1 & \omega^ {6} & \omega^ {4} & \omega^ {2} \\ 1 & \omega^ {7} & \omega^ {6} & \omega^ {5} & \omega^ {4} & \omega^ {3} & \omega^ {2} & \omega \end{array} \right] \qquad (\omega = \omega_ {8}).
$$

(Note that $\omega _ { 8 }$ is a root of unity so that high powers simplify, e.g., $[ F _ { 8 } ] _ { 4 , 7 } = \omega ^ { 3 \cdot 6 } =$ $\dot { \omega } ^ { 1 8 } = \omega ^ { 2 } . )$ If co $s = [ 1 3 5 7 2 4 6 8 ]$ , then

$$
F _ {8} (:, c o l s) = \left[ \begin{array}{c c c c c c c c} 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\ 1 & \omega^ {2} & \omega^ {4} & \omega^ {6} & \omega & \omega^ {3} & \omega^ {5} & \omega^ {7} \\ 1 & \omega^ {4} & 1 & \omega^ {4} & \omega^ {2} & \omega^ {6} & \omega^ {2} & \omega^ {6} \\ 1 & \omega^ {6} & \omega^ {4} & \omega^ {2} & \omega^ {3} & \omega & \omega^ {7} & \omega^ {5} \\ \hline 1 & 1 & 1 & 1 & - 1 & - 1 & - 1 & - 1 \\ 1 & \omega^ {2} & \omega^ {4} & \omega^ {6} & - \omega & - \omega^ {3} & - \omega^ {5} & - \omega^ {7} \\ 1 & \omega^ {4} & 1 & \omega^ {4} & - \omega^ {2} & - \omega^ {6} & - \omega^ {2} & - \omega^ {6} \\ 1 & \omega^ {6} & \omega^ {4} & \omega^ {2} & - \omega^ {3} & - \omega & - \omega^ {7} & - \omega^ {5} \end{array} \right].
$$

The lines through the matrix are there to help us think of $F _ { 8 } ( : , c o l s )$ as a 2-by-2 matrix with 4-by-4 blocks. Noting that $\omega ^ { 2 } = \omega _ { 8 } ^ { 2 } = \omega _ { 4 }$ , we see that

$$
F _ {8} (:, c o l s) = \left[ \begin{array}{c c} F _ {4} & \Omega_ {4} F _ {4} \\ \hline F _ {4} & - \Omega_ {4} F _ {4} \end{array} \right]
$$

where $\Omega _ { 4 } = \mathrm { d i a g } ( 1 , \omega _ { 8 } , \omega _ { 8 } ^ { 2 } , \omega _ { 8 } ^ { 3 } )$ . It follows that if $\boldsymbol { x } \in \mathbb { R } ^ { 8 }$ , then

$$
F _ {8} x = F _ {8} (:, c o l s) \cdot x (c o l s) = \left[ \begin{array}{c c} F _ {4} & \Omega_ {4} F _ {4} \\ \hline F _ {4} & - \Omega_ {4} F _ {4} \end{array} \right] \left[ \begin{array}{c} x (1: 2: 8) \\ \hline x (2: 2: 8) \end{array} \right] = \left[ \begin{array}{c c} I _ {4} & \Omega_ {4} \\ \hline I _ {4} & - \Omega_ {4} \end{array} \right] \left[ \begin{array}{c} F _ {4} x (1: 2: 8) \\ \hline F _ {4} x (2: 2: 8) \end{array} \right].
$$

Thus, by simple scalings we can obtain the 8-point DFT $y = F _ { 8 } x$ from the 4-point DFTs $y _ { T } = F _ { 4 } { \cdot } x ( 1 { : } 2 { : } 8 )$ and $y _ { B } = F _ { 4 } { \cdot } x ( 2 { : } 2 { : } 8 )$ . In particular,

$$
y (1: 4) = y _ {T} + d. * y _ {B},
$$

$$
y (5: 8) = y _ {T} - d. * y _ {B}
$$

where

$$
d = \left[ \begin{array}{c} 1 \\ \omega \\ \omega^ {2} \\ \omega^ {3} \end{array} \right].
$$

More generally, if $n = 2 m$ , then $y = F _ { n } x$ is given by

$$
y (1: m) = y _ {T} + d. * y _ {B},
$$

$$
y (m + 1: n) = y _ {T} - d. * y _ {B}
$$

where $d = \left[ 1 , \omega _ { n } , \ldots , \omega _ { n } ^ { m - 1 } \right] ^ { T }$ and

$$
y _ {T} = F _ {m} x (1: 2: n),
$$

$$
y _ {B} = F _ {m} x (2: 2: n).
$$

For $n = 2 ^ { t }$ , we can recur on this process until $n = 1$ , noting that $F _ { 1 } x = x$

Algorithm 1.4.1 If $x \in \mathbb { C } ^ { n }$ and $n = 2 ^ { t }$ , then this algorithm computes the discrete Fourier transform $y = F _ { n } x$ .

function $y = \mathsf { f f t } ( x , n )$

if n = 1

$$
y = x
$$

else

$$
m = n / 2
$$

$$
y _ {T} = \operatorname{fft} (x (1: 2: n), m)
$$

$$
y _ {B} = \operatorname{fft} (x (2: 2: n), m)
$$

$$
\omega = \exp (- 2 \pi i / n)
$$

$$
d = \left[ 1, \omega , \dots , \omega^ {m - 1} \right] ^ {T}
$$

$$
z = d. * y _ {B}
$$

$$
y = \left[ \begin{array}{l} y _ {T} + z \\ y _ {T} - z \end{array} \right]
$$

end

The flop analysis of fft requires an assessment of complex arithmetic and the solution of an interesting recursion. We first observe that the multiplication of two complex numbers involves six (real) flops while the addition of two complex numbers involves two flops. Let $f _ { n }$ be the number of flops that fft needs to produce the DFT of $x \in \mathbb { C } ^ { n }$ . Scrutiny of the method reveals that

$$
\left\{ \begin{array}{l} y _ {T} \\ y _ {B} \\ d \\ z \\ y \end{array} \right\} \text {requires} \left\{ \begin{array}{l} f _ {m} \text {flops} \\ f _ {m} \text {flops} \\ 6 m \text {flops} \\ 6 m \text {flops} \\ 2 n \text {flops} \end{array} \right\}
$$

where $n = 2 m$ . Thus,

$$
f _ {n} = 2 f _ {m} + 8 n \quad (f _ {1} = 0).
$$

Conjecturing that $f _ { n } = c { \cdot } n \log _ { 2 } ( n )$ for some constant $c ,$ it follows that

$$
f _ {n} = c \cdot n \log_ {2} (n) = 2 c \cdot m \log_ {2} (m) + 8 n = c \cdot n (\log_ {2} (n) - 1) + 8 n,
$$

from which we conclude that $c = 8$ . Thus, fft requires 8n $\log _ { 2 } ( n )$ flops. Appreciate the speedup over conventional matrix-vector multiplication. If $n = 2 ^ { 2 0 }$ , it is a factor of about 10,000. We mention that the fft flop count can be reduced to $5 n \log _ { 2 } ( n )$ by precomputing $\omega _ { n } , \ldots , \omega _ { n } ^ { n / 2 - 1 }$ . See P1.4.1.

# 1.4.2 Fast Sine and Cosine Transformations

In the discrete sine transform (DST) problem, we are given real values $x _ { 1 } , \ldots , x _ { m - 1 }$ and compute

$$
y _ {k} = \sum_ {j = 1} ^ {m - 1} \sin \left(\frac {k j \pi}{m}\right) x _ {j} \tag {1.4.3}
$$

for $k = 1 { : } m - 1$ . In the discrete cosine transform (DCT) problem, we are given real values $x _ { 0 } , x _ { 1 } , \ldots , x _ { m }$ and compute

$$
y _ {k} = \frac {x _ {0}}{2} + \sum_ {j = 1} ^ {m - 1} \cos \left(\frac {k j \pi}{m}\right) x _ {j} + \frac {(- 1) ^ {k} x _ {m}}{2} \tag {1.4.4}
$$

for $k = 0 { : } m$ . Note that the sine and cosine evaluations “show $\mathrm { u p } ^ { \mathrm { , } }$ in the DFT matrix. Indeed, for $k = 0 { : } 2 m - 1$ and $j = 0 { : } 2 m - 1$ we have

$$
[ F _ {2 m} ] _ {k + 1, j + 1} = \omega_ {2 m} ^ {k j} = \cos \left(\frac {k j \pi}{m}\right) - i \sin \left(\frac {k j \pi}{m}\right). \tag {1.4.5}
$$

This suggests (correctly) that there is an exploitable connection between each of these trigonometric transforms and the DFT. The key observation is to block properly the real and imaginary parts of $F _ { 2 m }$ . To that end, define the matrices $S _ { r } \in \mathbb { R } ^ { r \times r }$ and $C _ { r } \in \mathbb { R } ^ { r \times r }$ by

$$
[ S _ {r} ] _ {k j} = \sin \left(\frac {k j \pi}{r + 1}\right), \quad k = 1: r, j = 1: r. \tag {1.4.6}
$$

$$
[ C _ {r} ] _ {k j} = \cos \left(\frac {k j \pi}{r + 1}\right),
$$

Recalling from §1.2.11 the definition of the exchange permutation ${ \mathcal { E } } _ { n }$ , we have

Theorem 1.4.1. Let m be a positive integer and define the vectors e and v by

$$
e ^ {T} = (\underbrace {1 , 1 , \ldots , 1} _ {m - 1}), \qquad v ^ {T} = (\underbrace {- 1 , 1 , \ldots , (- 1) ^ {m - 1}} _ {m - 1}).
$$

$I f E = \mathcal { E } _ { m - 1 } , C = C _ { m - 1 }$ , and $S \ = \ S _ { m - 1 }$ , then

$$
F _ {2 m} = \left[ \begin{array}{c c c c} 1 & e ^ {T} & 1 & e ^ {T} \\ e & C - i S & v & (C + i S) E \\ 1 & v ^ {T} & (- 1) ^ {m} & v ^ {T} E \\ e & E (C + i S) & E v & E (C - i S) E \end{array} \right]. \tag {1.4.7}
$$

Proof. It is clear from (1.4.5) that $F _ { 2 m } ( : , 1 ) , F _ { 2 m } ( 1 , : 1 ) , F _ { 2 m } ( : , m + 1 )$ , and $F _ { 2 m } ( m { + } 1 , : )$ （ are correctly specified. It remains for us to show that equation (1.4.7) holds in blocks positions (2,2), (2,4), (4,2), and (4,4). The (2,2) verification is straightforward:

$$
\begin{array}{l} [ F _ {2 m} (2: m, 2: m) ] _ {k j} = \cos \left(\frac {k j \pi}{m}\right) - i \sin \left(\frac {k j \pi}{m}\right) \\ = [ C - i S ] _ {k j}. \\ \end{array}
$$

A little trigonometry is required to verify correctness in the (2,4) position:

$$
\begin{array}{l} [ F _ {2 m} (2: m, m + 2: 2 m) ] _ {k j} = \cos \left(\frac {k (m + j) \pi}{m}\right) - i \sin \left(\frac {k (m + j) \pi}{m}\right) \\ = \cos \left(\frac {k j \pi}{m} + k \pi\right) - i \sin \left(\frac {k j \pi}{m} + k \pi\right) \\ = \cos \left(- \frac {k j \pi}{m} + k \pi\right) + i \sin \left(- \frac {k j \pi}{m} + k \pi\right) \\ = \cos \left(\frac {(k (m - j) \pi)}{m}\right) + i \sin \left(\frac {k (m - j) \pi}{m}\right) \\ = \left[ (C + i S) E \right] _ {k j}. \\ \end{array}
$$

We used the fact that post-multiplying a matrix by the permutation $E = \mathcal { E } _ { m - 1 }$ has the effect of reversing the order of its columns. The recipes for $F _ { 2 m } ( m + 2 { : } 2 m , 2 { : } m )$ and $F _ { 2 m } ( m + 2 ; 2 m , m + 2 ; 2 m )$ are derived similarly.

Using the notation of the theorem, we see that the sine transform (1.4.3) is a matrix-vector product

$$
y (1: m - 1) = \mathrm{DST} (m - 1) \cdot x (1: m - 1)
$$

where

$$
\mathrm{DST} (m - 1) = S _ {m - 1}. \tag {1.4.8}
$$

If $x = x ( 1 { : } m - 1 )$ and

$$
x _ {\sin} = \left[ \begin{array}{c} 0 \\ x \\ 0 \\ - E x \end{array} \right] \in \mathbb {R} ^ {2 m}, \tag {1.4.9}
$$

then since $e ^ { T } E = e$ and $E ^ { 2 } = E$ we have

$$
\frac {i}{2} F _ {2 m} x _ {\sin} = \frac {i}{2} \left[ \begin{array}{c c c c} 1 & e ^ {T} & 1 & e ^ {T} \\ e & C - i S & v & (C + i S) E \\ 1 & v ^ {T} & (- 1) ^ {m} & v ^ {T} E \\ e & E (C + i S) & E v & E (C - i S) E \end{array} \right] \left[ \begin{array}{c} 0 \\ x \\ 0 \\ - E x \end{array} \right]
$$

$$
= \frac {i}{2} \left[ \begin{array}{c} e ^ {T} x - e ^ {T} E x \\ - 2 i S x \\ v ^ {T} x - v ^ {T} E ^ {2} x \\ i (E S x + E S E ^ {2} x) \end{array} \right] = \left[ \begin{array}{c} 0 \\ S x \\ 0 \\ - E S x \end{array} \right].
$$

Thus, the DST of $x ( 1 { : } m - 1 )$ is a scaled subvector of $F _ { 2 m } x _ { \mathrm { s i n } }$ .

Algorithm 1.4.2 The following algorithm assigns the DST of $x _ { 1 } , \ldots , x _ { m - 1 }$ to $y$

Set up the vector $x _ { \mathrm { s i n } }$ defined by (1.4.9).

Use fft (e.g., Algorithm 1.4.1) to compute ${ \tilde { y } } = F _ { 2 m } x _ { \sin }$

$$
y = i \cdot \tilde {y} (2: m) / 2
$$

This computation involves $O ( m \log _ { 2 } ( m ) )$ flops. We mention that the vector $x _ { \mathrm { s i n } }$ is real and highly structured, something that would be exploited in a truly efficient implementation.

Now let us consider the discrete cosine transform defined by (1.4.4). Using the notation from Theorem 1.4.1, the DCT is a matrix-vector product

$$
y (0: m) = \mathrm{DCT} (m + 1) \cdot x (0: m)
$$

where

$$
\mathrm{DCT} (m + 1) = \left[ \begin{array}{c c c} 1 / 2 & e ^ {T} & 1 / 2 \\ e / 2 & C _ {m - 1} & v / 2 \\ 1 / 2 & v ^ {T} & (- 1) ^ {m} / 2 \end{array} \right] \tag {1.4.10}
$$

$\mathrm { I f } \ \tilde { x } = x ( 1 { : } m - 1 )$ and

$$
x _ {\cos} = \left[ \begin{array}{c} x _ {0} \\ \tilde {x} \\ x _ {m} \\ E \tilde {x} \end{array} \right] \in \mathbb {R} ^ {2 m}, \tag {1.4.11}
$$

then

$$
\frac {1}{2} F _ {2 m} x _ {\mathrm{cos}} = \frac {1}{2} \left[ \begin{array}{c c c c} 1 & e ^ {T} & 1 & e ^ {T} \\ e & C - i S & v & (C + i S) E \\ 1 & v ^ {T} & (- 1) ^ {m} & v ^ {T} E \\ e & E (C + i S) & E v & E (C - i S) E \end{array} \right] \left[ \begin{array}{c} x _ {0} \\ \tilde {x} \\ x _ {m} \\ E \tilde {x} \end{array} \right]
$$

$$
= \left[ \begin{array}{c c c c c} (x _ {0} / 2) & + & e ^ {T} \tilde {x} & + & (x _ {m} / 2) \\ (x _ {0} / 2) e & + & C \tilde {x} & + & (x _ {m} / 2) v \\ (x _ {0} / 2) & + & v ^ {T} \tilde {x} & + & (- 1) ^ {m} (x _ {m} / 2) \\ (x _ {0} / 2) e & + & E C \tilde {x} & + & (x _ {m} / 2) E v \end{array} \right].
$$

Notice that the top three components of this block vector define the DCT of $x ( 0 { : } m )$ . Thus, the DCT is a scaled subvector of $F _ { \mathrm { 2 } m } x _ { \mathrm { c o s } }$ .

Algorithm 1.4.3 The following algorithm assigns to $y \in \mathbb { R } ^ { m + 1 }$ the DCT of $x _ { 0 } , \ldots , x _ { m }$

Set up the vector $x _ { \mathrm { c o s } } \in \mathbb { R } ^ { 2 m }$ defined by (1.4.11).

Use fft (e.g., Algorithm 1.4.1) to compute ${ \tilde { y } } = F _ { 2 m } x _ { \mathrm { { c o s } } }$

$$
y = \tilde {y} (1: m + 1) / 2
$$

This algorithm requires $O ( m \log m )$ flops, but as with Algorithm 1.4.2, it can be more efficiently implemented by exploiting symmetries in the vector $x _ { c o s }$ .

We mention that there are important variants of the DST and the DCT that can be computed fast:

$$
\text { DST - II: } \quad y _ {k} = \sum_ {j = 1} ^ {m} \sin \left(\frac {k (2 j - 1) \pi}{2 m}\right) x _ {j}, \quad k = 1: m,
$$

$$
\text { DST - III: } \quad y _ {k} = \sum_ {j = 1} ^ {m} \sin \left(\frac {(2 k - 1) j \pi}{2 m}\right) x _ {j}, \quad k = 1: m,
$$

$$
\text {DST - IV:} \quad y _ {k} = \sum_ {j = 1} ^ {m} \sin \left(\frac {(2 k - 1) (2 j - 1) \pi}{2 m}\right) x _ {j}, \quad k = 1: m, \tag {1.4.12}
$$

$$
\text { DCT - II: } \quad y _ {k} = \sum_ {j = 0} ^ {m - 1} \cos \left(\frac {k (2 j - 1) \pi}{2 m}\right) x _ {j}, \quad k = 0: m - 1,
$$

$$
\text { DCT - III: } \quad y _ {k} = \frac {x _ {0}}{2} = \sum_ {j = 1} ^ {m - 1} \cos \left(\frac {(2 k - 1) j \pi}{2 m}\right) x _ {j}, \quad k = 0: m - 1,
$$

$$
\text { DCT - IV: } \quad y _ {k} = \sum_ {j = 0} ^ {m - 1} \cos \left(\frac {(2 k - 1) (2 j - 1) \pi}{2 m}\right) x _ {j}, \qquad k = 0: m - 1.
$$

For example, if $\tilde { y } \in \mathbb { R } ^ { 2 m - 1 }$ is the DST of $\tilde { \boldsymbol { x } } = \left[ x _ { 1 } , 0 , x _ { 2 } , 0 , \dots , 0 , x _ { m - 1 } , x _ { m } \right] ^ { T }$ , then ${ \tilde { y } } ( 1 { : } m )$ is the DST-II of $\boldsymbol { x } \in \mathbb { R } ^ { m }$ . See Van Loan (FFT) for further details.

# 1.4.3 The Haar Wavelet Transform

If $n = 2 ^ { t }$ , then the Haar wavelet transform $y = W _ { n } x$ is a matrix-vector product in which the transform matrix $W _ { n } \in \mathbb { R } ^ { n \times n }$ is defined recursively:

$$
W _ {n} = \left\{ \begin{array}{l l} \left[ \begin{array}{c} W _ {m} \otimes \binom{1}{1} \end{array} \right| I _ {m} \otimes \binom{1}{- 1} \Bigg ] & \text { if } n = 2 m, \\ [ 1 ] & \text { if } n = 1. \end{array} \right.
$$

Here are some examples:

$$
W _ {2} = \left[ \begin{array}{c c} 1 & 1 \\ \hline 1 & - 1 \end{array} \right],
$$

$$
W _ {4} = \left[ \begin{array}{c c c c} 1 & 1 & 1 & 0 \\ 1 & 1 & - 1 & 0 \\ \hline 1 & - 1 & 0 & 1 \\ 1 & - 1 & 0 & - 1 \end{array} \right],
$$

$$
W _ {8} = \left[ \begin{array}{c c c c c c c c} 1 & 1 & 1 & 0 & 1 & 0 & 0 & 0 \\ 1 & 1 & 1 & 0 & - 1 & 0 & 0 & 0 \\ 1 & 1 & - 1 & 0 & 0 & 1 & 0 & 0 \\ 1 & 1 & - 1 & 0 & 0 & - 1 & 0 & 0 \\ \hline 1 & - 1 & 0 & 1 & 0 & 0 & 1 & 0 \\ 1 & - 1 & 0 & 1 & 0 & 0 & - 1 & 0 \\ 1 & - 1 & 0 & - 1 & 0 & 0 & 0 & 1 \\ 1 & - 1 & 0 & - 1 & 0 & 0 & 0 & - 1 \end{array} \right].
$$

An interesting block pattern emerges if we reorder the rows of $W _ { n }$ so that the oddindexed rows come first:

$$
\mathcal {P} _ {2, m} ^ {T} W _ {n} = \left[ \begin{array}{c c} W _ {m} & I _ {m} \\ W _ {m} & - I _ {m} \end{array} \right] = (W _ {2} \otimes I _ {m}) \left[ \begin{array}{c c} W _ {m} & 0 \\ 0 & I _ {m} \end{array} \right]. \tag {1.4.13}
$$

Thus, if $x \in \mathbb { R } ^ { n } , x _ { T } = x ( 1 { : } m )$ , and $x _ { B } = x ( m + 1 { : } n )$ , then

$$
\begin{array}{l} y = W _ {n} x = \mathcal {P} _ {2, m} \left[ \begin{array}{c c} I _ {m} & I _ {m} \\ I _ {m} & - I _ {m} \end{array} \right] \left[ \begin{array}{c c} W _ {m} & 0 \\ 0 & I _ {m} \end{array} \right] \left[ \begin{array}{c} x _ {T} \\ x _ {B} \end{array} \right] \\ = \mathcal {P} _ {2, m} \left[ \begin{array}{l} W _ {m} x _ {T} + x _ {B} \\ W _ {m} x _ {T} - x _ {B} \end{array} \right]. \\ \end{array}
$$

In other words,

$$
y (1: 2: n) = W _ {m} x _ {T} + x _ {B}, \qquad y (2: 2: n) = W _ {m} x _ {T} - x _ {B}.
$$

This points the way to a fast recursive procedure for computing $y = W _ { n } x$ .

Algorithm 1.4.4 (Haar Wavelet Transform) If $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $n = 2 ^ { t }$ , then this algorithm computes the Haar transform $y = W _ { n } x$ .

function $y = \mathsf { f h t } ( x , n )$

if $n = 1$

$$
y = x
$$

else

$$
m = n / 2
$$

$$
z = \operatorname{fht} (x (1: m), m)
$$

$$
y (1: 2: m) = z + x (m + 1: n)
$$

$$
y (2: 2: m) = z - x (m + 1: n)
$$

end

It can be shown that this algorithm requires 2n flops.

# Problems

P1.4.1 Suppose $w = \left[ 1 , \omega _ { n } , \omega _ { n } ^ { 2 } , . . . , \omega _ { n } ^ { n / 2 - 1 } \right]$ where $n = 2 ^ { t }$ . Using the colon notation, express

$$
\left[ 1, \omega_ {r}, \omega_ {r} ^ {2}, \ldots , \omega_ {r} ^ {r / 2 - 1} \right]
$$

as a subvector of w where $r = 2 ^ { q } , q = 1 { : } t$ . Rewrite Algorithm 1.4.1 with the assumption that w is precomputed. Show that this maneuver reduces the flop count to 5n $\log _ { 2 } n$ .

P1.4.2 Suppose $n = 3 m$ and examine

$$
G = \left[ F _ {n} (:, 1: 3: n - 1) \mid F _ {n} (:, 2: 3: n - 1) \mid F _ {n} (:, 3: 3: n - 1) \right]
$$

as a 3-by-3 block matrix, looking for scaled copies of $F _ { m }$ . Based on what you find, develop a recursive radix-3 FFT analogous to the radix-2 implementation in the text.

P1.4.3 If $n = 2 ^ { t }$ , then it can be shown that $F _ { n } = ( A _ { t } \Gamma _ { t } ) \cdot \cdot \cdot ( A _ { 1 } \Gamma _ { 1 } )$ where for $q = 1 { : } t$

$$
L _ {q} = 2 ^ {q}, \quad r _ {q} = n / L _ {q},
$$

$$
A _ {q} = I _ {r _ {q}} \otimes \left[ \begin{array}{c c} I _ {L _ {q - 1}} & \Omega_ {q} \\ I _ {L _ {q - 1}} & - \Omega_ {q} \end{array} \right],
$$

$$
\Gamma_ {q} = \mathcal {P} _ {2, r _ {q}} \otimes I _ {L _ {q - 1}},
$$

$$
\Omega_ {q} = \mathrm{diag} (1, \omega_ {L _ {q}}, \ldots , \omega_ {L _ {q}} ^ {L _ {q - 1} - 1}).
$$

Note that with this factorization, the DFT $y = F _ { n } x$ can be computed as follows:

$$
y = x
$$

for $q = 1 { : } t$

$$
y = A _ {q} (\Gamma_ {q} y)
$$

end

Fill in the details associated with the y updates and show that a careful implementation requires 5n $\log _ { 2 } ( n )$ flops.

P1.4.4 What fraction of the components of $W _ { n }$ are zero?

P1.4.5 Using (1.4.13), verify by induction that if $n = 2 ^ { t }$ , then the Haar tranform matrix $W _ { n }$ has the factorization $W _ { n } = H _ { t } \cdot \cdot \cdot H _ { 1 }$ where

$$
H _ {q} = \left[ \begin{array}{c c} \mathcal {P} _ {2, L _ {*}} & 0 \\ 0 & I _ {n - L} \end{array} \right] \left[ \begin{array}{c c} W _ {2} \otimes I _ {L _ {*}} & 0 \\ 0 & I _ {n - L} \end{array} \right] \qquad L = 2 ^ {q}, L _ {*} = L / 2.
$$

Thus, the computation of $y = W _ { n } x$ may proceed as follows:

$$
y = x
$$

for $q = 1 { : } t$

$$
y = H _ {q} y
$$

end

Fill in the details associated with the update $y = H _ { q } y$ and confirm that $W _ { n } x$ costs 2n flops.

P1.4.6 Using (1.4.13), develop an O(n) procedure for solving $W _ { n } y = x$ where $\boldsymbol { x } \in \mathbb { R } ^ { n }$ is given and $n = 2 ^ { t }$ .

# Notes and References for §1.4

In Van Loan (FFT) the FFT family of algorithms is described in the language of matrix-factorizations. A discussion of various fast trigonometric transforms is also included. See also:

W.L. Briggs and V.E. Henson (1995). The DFT: An Owners’ Manual for the Discrete Fourier Transform, SIAM Publications, Philadelphia, PA.

The design of a high-performance FFT is a nontrivial task. An important development in this regard is a software tool known as “the fastest Fourier transform in the west”:

M. Frigo and S.G. Johnson (2005). “The Design and Implementation of FFTW3”, Proceedings of the IEEE, 93, 216–231.

It automates the search for the “right” FFT given the underlying computer architecture. FFT references that feature interesting factorization and approximation ideas include:

A. Edelman, P. McCorquodale, and S. Toledo (1998). “The Future Fast Fourier Transform?,” SIAM J. Sci. Comput. 20, 1094–1114.   
A. Dutt and and V. Rokhlin (1993). “Fast Fourier Transforms for Nonequally Spaced Data,” SIAM J. Sci. Comput. 14, 1368–1393.   
A. F. Ware (1998). “Fast Approximate Fourier Transforms for Irregularly Spaced Data,” SIAM Review 40, 838 –856.   
N. Nguyen and Q.H. Liu (1999). “The Regular Fourier Matrices and Nonuniform Fast Fourier Transforms,” SIAM J. Sci. Comput. 21, 283–293.   
A. Nieslony and G. Steidl (2003). “Approximate Factorizations of Fourier Matrices with Nonequispaced Knots,” Lin. Alg. Applic. 366, 337–351.   
L. Greengard and J.–Y. Lee (2004). “Accelerating the Nonuniform Fast Fourier Transform,” SIAM Review 46, 443–454.   
K. Ahlander and H. Munthe-Kaas (2005). “Applications of the Generalized Fourier Transform in Numerical Linear Algebra,” BIT 45, 819–850.

The fast multipole method and the fast Gauss transform represent another type of fast transform that is based on a combination of clever blocking and approximation.

L. Greengard and V. Rokhlin (1987). “A Fast Algorithm for Particle Simulation,” J. Comput. Phys. 73, 325–348.   
X. Sun and N.P. Pitsianis (2001). “A Matrix Version of the Fast Multipole Method,” SIAM Review 43, 289–300.   
L. Greengard and J. Strain (1991). “The Fast Gauss Transform,” SIAM J. Sci. Stat. Comput. 12, 79–94.   
M. Spivak, S.K. Veerapaneni, and L. Greengard (2010). “The Fast Generalized Gauss Transform,” SIAM J. Sci. Comput. 32, 3092–3107.   
X. Sun and Y. Bao (2003). “A Kronecker Product Representation of the Fast Gauss Transform,” SIAM J. Matrix Anal. Applic. 24, 768–786.

The Haar transform is a simple example of a wavelet transform. The wavelet idea has had a profound impact throughout computational science and engineering. In many applications, wavelet basis functions work better than the sines and cosines that underly the DFT. Excellent monographs on this subject include

I Daubechies (1992). Ten Lectures on Wavelets, SIAM Publications, Philadelphia, PA.   
G. Strang (1993). “Wavelet Transforms Versus Fourier Transforms,” Bull. AMS 28, 288–305.   
G. Strang and T. Nguyan (1996). Wavelets and Filter Banks, Wellesley-Cambridge Press.
