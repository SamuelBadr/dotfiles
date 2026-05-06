# 2.7 Finite Precision Matrix Computations

Rounding errors are part of what makes the field of matrix computations so challenging. In this section we describe a model of floating point arithmetic and then use it to develop error bounds for floating point dot products, saxpys, matrix-vector products, and matrix-matrix products.

# 2.7.1 A 3-digit Calculator

Suppose we have a base-10 calculator that represents nonzero numbers in the following style:

$$
x = \pm d _ {0}. d _ {1} d _ {2} \times 1 0 ^ {e} \qquad \text {where} \quad \left\{ \begin{array}{r l} 1 & \leq d _ {0} \leq 9, \\ 0 & \leq d _ {1} \leq 9, \\ 0 & \leq d _ {2} \leq 9, \\ - 9 & \leq e \leq 9. \end{array} \right.
$$

Let us call these numbers floating point numbers. After playing around a bit we make a number of observations:

• The precision of the calculator has to do with the “length” of the significand $d _ { 0 } . d _ { 1 } d _ { 2 }$ . For example, the number π would be represented as $3 . 1 4 \times 1 0 ^ { 0 }$ , which has a relative error approximately equal to $1 0 ^ { - 3 }$ .

• There is not enough “room” to store exactly the results from most arithmetic operations between floating point numbers. Sums and products like

$$
(1. 2 3 \times 1 0 ^ {6}) + (4. 5 6 \times 1 0 ^ {4}) = 1 2 7 5 6 0 0,
$$

$$
(1. 2 3 \times 1 0 ^ {1}) * (4. 5 6 \times 1 0 ^ {2}) = 5 6 0 8. 8
$$

involve more than three significant digits. Results must be rounded in order to $\mathrm { \Omega ^ { \circ } f i t { \Omega ^ { \circ } } }$ the 3-digit format, e.g., round $( 1 2 7 5 6 0 0 ) = 1 . 2 8 \times 1 0 ^ { 6 }$ , round(5608.8) = $5 . 6 1 \times 1 0 ^ { 3 }$ .

• If zero is to be a floating point number (and it must be), then we need a special convention for its representation, $\mathrm { e . g . , 0 . 0 0 \times 1 0 ^ { 0 } }$ .

• In contrast to the real numbers, there is a smallest positive floating point number $( N _ { \mathrm { m i n } } = 1 . 0 0 { \times } 1 0 ^ { - 9 } )$ and there is a largest positive floating point number $( N _ { \mathrm { m a x } } =$ $9 . 9 9 \times 1 0 ^ { 9 } )$ .

• Some operations yield answers whose exponents exceed the 1-digit allocation, $\mathrm { e . g . , ( 1 . 2 3 \times 1 0 ^ { 4 } ) * ( 4 . 5 6 \times 1 0 ^ { 7 } ) }$ and $( 1 . 2 3 \times 1 0 ^ { - 2 } ) / ( 4 . 5 6 \times 1 0 ^ { 8 } )$ .

• The set of floating point numbers is finite. For the toy calculator there are $2 \times 9 \times 1 0 \times 1 0 \times 1 9 + 1 = 3 4 2 0 1$ floating point numbers.

• The spacing between the floating point numbers varies. Between $1 . 0 0 \times 1 0 ^ { e }$ and $1 . 0 0 \times \bar { 1 0 ^ { e + 1 } }$ the spacing is $1 0 ^ { e - \hat { 2 } }$ .

The careful design and analysis of a floating point computation requires an understanding of these inexactitudes and limitations. How are results rounded? How accurate is floating point arithmetic? What can we say about a sequence of floating point operations?

# 2.7.2 IEEE Floating Point Arithmetic

To build a solid, practical understanding of finite precision computation, we set aside our toy, motivational base-10 calculator and consider the key ideas behind the widely accepted IEEE floating point standard. The IEEE standard includes a 32-bit single format and a 64-bit double format. We will illustrate concepts using the latter as an example because typical accuracy requirements make it the format of choice.

The importance of having a standard for floating point arithmetic that is upheld by hardware manufacturers cannot be overstated. After all, floating point arithmetic is the foundation upon which all of scientific computing rests. The IEEE standard promotes software reliability and enables numerical analysts to make rigorous statements about computed results. Our discussion is based on the excellent book by Overton (2001).

The 64-bit double format allocates a single bit for the sign of the floating point number, 52 bits for the mantissa , and eleven bits for the exponent:

$$
x: \boxed {\pm} \left| a _ {1} a _ {2} \dots a _ {1 1} \right| b _ {1} b _ {2} \dots b _ {5 2}. \tag {2.7.1}
$$

The “formula” for the value of this representation depends upon the exponent bits:

If $a _ { 1 } \ldots a _ { 1 1 }$ is neither all 0’s nor all 1’s, then x is a normalized floating point number with value

$$
x = \pm (1. b _ {1} b _ {2} \dots b _ {5 2}) _ {2} \times 2 ^ {(a _ {1} a _ {2} \dots a _ {1 1}) _ {2} - 1 0 2 3}. \tag {2.7.2}
$$

The “1023 bias” in the exponent supports the graceful inclusion of various “unnormalized” floating numbers which we describe shortly. Several important quantities capture the finiteness of the representation. The machine epsilon is the gap between 1 and the next largest floating point number. Its value is $2 ^ { - 5 2 } \approx 1 0 ^ { - 1 6 }$ for the double format. Among the positive normalized floating point numbers, $N _ { \mathrm { m i n } } = 2 ^ { - 1 0 2 2 } \approx 1 0 ^ { - 3 0 8 }$ is the smallest and $N _ { \mathrm { { m a x } } } = ( 2 - 2 ^ { - 5 2 } ) 2 ^ { 1 0 2 3 } \approx 1 0 ^ { 3 0 8 }$ is the largest. A real number x is within the normalized range if $N _ { \mathrm { m i n } } \le | x | \le N _ { \mathrm { m a x } }$ .

If $a _ { 1 } \ldots a _ { 1 1 }$ is all $0 \mathrm { { s } }$ , then the value of the representation (2.7.1) is

$$
x = \pm \left(0. b _ {1} b _ {2} \dots b _ {5 2}\right) _ {2} \times 2 ^ {\left(a _ {1} a _ {2} \dots a _ {1 1}\right) _ {2} - 1 0 2 2} \tag {2.7.3}
$$

This includes 0 and the subnormal floating point numbers. This feature creates a uniform spacing of the floating point numbers between $- N _ { \mathrm { m i n } }$ and $+ N _ { \mathrm { m i n } }$ .

If $a _ { 1 } \ldots a _ { 1 1 }$ is all 1’s, then the encoding (2.7.1) represents inf for +∞, -inf for $- \infty$ , or NaN for “not-a-number.” The determining factor is the value of the $b _ { i }$ . (If the $b _ { i }$ are not all zero, then the value of x is NaN.) Quotients like $1 / 0 , - 1 / 0$ , and $0 / 0$ produce these special floating point numbers instead of prompting program termination.

There are four rounding modes: round down (toward $- \infty )$ , round up (toward $+ \infty )$ , round-toward-zero, and round-toward-nearest. We focus on round-toward-nearest since it is the mode almost always used in practice.

If a real number x is outside the range of the normalized floating point numbers then

$$
\operatorname{round} (x) = \left\{ \begin{array}{l l} - \infty & \text {if} x <   - N _ {\max} , \\ + \infty & \text {if} x > N _ {\max}. \end{array} \right.
$$

Otherwise, the rounding process depends upon its floating point “neighbors”:

$x _ { - }$ is the nearest floating point number to x that is $\leq x .$

$x _ { + }$ is the nearest floating point number to x that is $\geq x$

Define $d _ { - } = x - x _ { - }$ and $d _ { + } = x _ { + } - x$ and let “lsb” stand for “least significant bit.” If $N _ { \mathrm { m i n } } \le | x | \le N _ { \mathrm { m a x } }$ , then

$$
\operatorname{round} (x) = \left\{ \begin{array}{l l} x _ {-} & \text { if } d _ {-} <   d _ {+} \text { or } d _ {-} = d _ {+} \text { and } \operatorname{lsb} (x _ {-}) = 0, \\ x _ {+} & \text { if } d _ {+} <   d _ {-} \text { or } d _ {+} = d _ {-} \text { and } \operatorname{lsb} (x _ {+}) = 0. \end{array} \right.
$$

The tie-breaking criteria is well-defined because $x _ { - }$ and $x _ { + }$ are adjacent floating point numbers and so must differ in their least significant bit.

Regarding the accuracy of the round-to-nearest strategy, suppose x is a real number that satisfies $N _ { \mathrm { m i n } } \le | x | \le N _ { \mathrm { m a x } }$ . Thus,

$$
| \operatorname{round} (x) - x | \leq \frac {2 ^ {- 5 2}}{2} 2 ^ {e} \leq \frac {2 ^ {- 5 2}}{2} | x |
$$

which says that relative error is bounded by half of the machine epsilon:

$$
\frac {| \operatorname{round} (x) - x |}{| x |} \leq 2 ^ {- 5 3}.
$$

The IEEE standard stipulates that each arithmetic operation be correctly rounded, meaning that the computed result is the rounded version of the exact result. The implementation of correct rounding is far from trivial and requires registers that are equipped with several extra bits of precision.

We mention that the IEEE standard also requires correct rounding in the square root operation, the remainder operation, and various format conversion operations.

# 2.7.3 The “fl” Notation

With intuition gleaned from the toy calculator example and an understanding of IEEE arithmetic, we are ready to move on to the roundoff analysis of some basic algebraic calculations. The challenge when presenting the effects of finite precision arithmetic in this section and throughout the book is to communicate essential behavior without excessive detail. To that end we use the notation fl(·) to identify a floating point storage and/or computation. Unless exceptions are a critical part of the picture, we freely invoke the fl notation without mentioning “−∞,” “∞,” “NaN,” etc.

If $x \in \mathbb { R }$ , then fl(x) is its floating point representation and we assume that

$$
\mathsf {f l} (x) = x (1 + \delta), \quad | \delta | \leq \mathbf {u}, \tag {2.7.4}
$$

where u is the unit roundoff defined by

$$
\mathbf {u} = \frac {1}{2} \times (\text { gap   between   1   and   next   largest   floating   point   number }). \tag {2.7.5}
$$

The unit roundoff for IEEE single format is about $1 0 ^ { - 7 }$ and for double format it is about 10−16.

If x and y are floating point numbers and “op” is any of the four arithmetic operations, then fl(x op y) is the floating point result from the floating point op. Following Trefethen and Bau (NLA), the fundamental axiom of floating point arithmetic is that

$$
\operatorname{fl} (x \text {   op   } y) = (x \text {   op   } y) (1 + \delta), \quad | \delta | \leq \mathbf {u}, \tag {2.7.6}
$$

where x and y are floating point numbers and the “op” inside the fl operation means “floating point operation.” This shows that there is small relative error associated with individual arithmetic operations:

$$
\frac {| \mathsf {f l} (x \textsf {o p} y) - (x \textsf {o p} y) |}{| x \textsf {o p} y |} \leq \mathbf {u}, \qquad x \textsf {o p} y \neq 0.
$$

Again, unless it is particularly relevant to the discussion, it will be our habit not to bring up the possibilities of an exception arising during the floating point operation.

# 2.7.4 Become a Floating Point Thinker

It is a good idea to have a healthy respect for the subleties of floating point calculation. So before we proceed with our first serious roundoff error analysis we offer three maxims to keep in mind when designing a practical matrix computation. Each reinforces the distinction between computer arithmetic and exact arithmetic.

# Maxim 1. Order is Important.

Floating point arithmetic is not associative. For example, suppose

$$
x = 1. 2 4 \times 1 0 ^ {0}, \qquad y = - 1. 2 3 \times 1 0 ^ {0}, \qquad z = 1. 0 0 \times 1 0 ^ {- 3}.
$$

Using toy calculator arithmetic we have

$$
\mathsf {f l} (\mathsf {f l} (x + y) + z)) = 1. 1 0 \times 1 0 ^ {- 2}
$$

while

$$
\operatorname{fl} (x + \operatorname{fl} (y + z)) = 1. 0 0 \times 1 0 ^ {- 2}.
$$

A consequence of this is that mathematically equivalent algorithms may produce different results in floating point.

# Maxim 2. Larger May Mean Smaller.

Suppose we want to compute the derivative of $f ( x ) = \sin ( x )$ using a divided difference. Calculus tells us that $d = ( \sin ( x + h ) - \sin ( x ) ) / h$ satisfies $| d - \cos ( x ) | = O ( h )$ which argues for making h as small as possible. On the other hand, any roundoff error sustained in the sine evaluations is magnified by $1 / h$ . By setting $h = { \sqrt { \mathbf { u } } } .$ , the sum of the calculus error and roundoff error is approximately minimized. In other words, a value of h much greater than u renders a much smaller overall error. See Overton(2001, pp. 70–72).

# Maxim 3. A Math Book Is Not Enough.

The explicit coding of a textbook formula is not always the best way to design an effective computation. As an example, we consider the quadratic equation $x ^ { 2 } - 2 p x - q =$ 0 where both p and q are positive. Here are two methods for computing the smaller (necessarily real) root:

$\mathrm { M e t h o d ~ 1 : } \quad r _ { \mathrm { m i n } } = p - \sqrt { p ^ { 2 } + q } ,$

$\mathrm { M e t h o d ~ 2 } { : } \quad r _ { \mathrm { m i n } } = \frac { q } { p + \sqrt { p ^ { 2 } + q } } .$

The first method is based on the familiar quadratic formula while the second uses the fact that −q is the product of $r _ { \mathrm { m i n } }$ and the larger root. Using IEEE double format arithmetic with input $p = 1 2 3 4 5 6 7 8$ and $q = 1$ we obtain these results:

$\mathrm { M e t h o d ~ 1 } { : } \quad r _ { \mathrm { m i n } } = - 4 . 0 9 7 8 1 9 3 2 8 3 0 8 1 0 6 \times 1 0 ^ { - 8 } ,$

$\mathrm { M e t h o d ~ 2 } { : } \quad r _ { \mathrm { m i n } } = - 4 . 0 5 0 0 0 0 0 3 3 2 1 0 0 0 2 1 \times 1 0 ^ { - 8 } \quad ( \mathrm { c o r r e c t } ) .$

Method 1 produces an answer that has almost no correct significant digits. It attempts to compute a small number by subtracting a pair of nearly equal large numbers. Almost all correct significant digits in the input data are lost during the subtraction, a phenomenon known as catastrophic cancellation. In contrast, Method 2 produces an answer that is correct to full machine precision. It computes a small number as a division of one number by a much larger number. See Forsythe (1970).

Keeping these maxims in mind does not guarantee the production of accurate, reliable software, but it helps.

# 2.7.5 Application: Storing a Real Matrix

Suppose $A \in \mathbb { R } ^ { m \times n }$ and that we wish to quantify the errors associated with its floating point representation. Denoting the stored version of A by ${ \mathsf { f l } } ( A )$ , we see that

$$
[ \mathsf {f l} (A) ] _ {i j} = \mathsf {f l} (a _ {i j}) = a _ {i j} (1 + \epsilon_ {i j}), \quad | \epsilon_ {i j} | \leq \mathbf {u}, \tag {2.7.7}
$$

for all i and j, i.e.,

$$
| \mathrm{fl} (A) - A | \leq \mathbf {u} | A |.
$$

A relation such as this can be easily turned into a norm inequality, e.g.,

$$
\| \mathsf {f l} (A) - A \| _ {1} \leq \mathbf {u} \| A \| _ {1}.
$$

However, when quantifying the rounding errors in a matrix manipulation, the absolute value notation is sometimes more informative because it provides a comment on each entry.

# 2.7.6 Roundoff in Dot Products

We begin our study of finite precision matrix computations by considering the rounding errors that result in the standard dot product algorithm:

$$
s = 0
$$

for $k = 1 { : } n$

$$
s = s + x _ {k} y _ {k} \tag {2.7.8}
$$

end

Here, x and y are n-by-1 floating point vectors.

In trying to quantify the rounding errors in this algorithm, we are immediately confronted with a notational problem: the distinction between computed and exact quantities. If the underlying computations are clear, we shall use the fl(·) operator to signify computed quantities. Thus, $\mathsf { I } ( x ^ { T } y )$ denotes the computed output of (2.7.8). Let us bound $| \mathsf { f l } ( x ^ { \hat { T } } y ) - x ^ { T } y |$ . If

$$
s _ {p} = \mathsf {f l} \left(\sum_ {k = 1} ^ {p} x _ {k} y _ {k}\right),
$$

then $s _ { 1 } = x _ { 1 } y _ { 1 } ( 1 + \delta _ { 1 } )$ with $| \delta _ { 1 } | \leq \mathbf { u }$ and for $p = 2 { : } n$

$$
\begin{array}{l} s _ {p} = \mathsf {f l} (s _ {p - 1} + \mathsf {f l} (x _ {p} y _ {p})) \\ = (s _ {p - 1} + x _ {p} y _ {p} (1 + \delta_ {p})) (1 + \epsilon_ {p}) \quad | \delta_ {p} |, | \epsilon_ {p} | \leq \mathbf {u}. \tag {2.7.9} \\ \end{array}
$$

A little algebra shows that

$$
\mathsf {f l} (x ^ {T} y) = s _ {n} = \sum_ {k = 1} ^ {n} x _ {k} y _ {k} (1 + \gamma_ {k})
$$

where

$$
(1 + \gamma_ {k}) = (1 + \delta_ {k}) \prod_ {j = k} ^ {n} (1 + \epsilon_ {j})
$$

with the convention that $\epsilon _ { 1 } = 0$ . Thus,

$$
| \mathfrak {f l} (x ^ {T} y) - x ^ {T} y | \leq \sum_ {k = 1} ^ {n} | x _ {k} y _ {k} | | \gamma_ {k} |. \tag {2.7.10}
$$

To proceed further, we must bound the quantities $| \gamma _ { k } |$ in terms of u. The following result is useful for this purpose.

Lemma 2.7.1. $I f \left( 1 + \alpha \right) = \prod _ { k = 1 } ^ { n } ( 1 + \alpha _ { k } )$ where $| \alpha _ { k } | \le \mathbf { u }$ and $n \mathbf { u } \leq . 0 1$ , then $| \alpha | \leq$ 1.01nu. k=1

Proof. See Higham (ASNA, p. 75).

Application of this result to (2.7.10) under the “reasonable” assumption $n \mathbf { u } \leq . 0 1$ gives

$$
\left| \mathrm{fl} \left(x ^ {T} y\right) - x ^ {T} y \right| \leq 1. 0 1 n \mathbf {u} | x | ^ {T} | y |. \tag {2.7.11}
$$

Notice that if $| x ^ { T } y | \ll | x | ^ { T } | y |$ , then the relative error in $\mathsf { f l } ( x ^ { T } y )$ may not be small.

# 2.7.7 Alternative Ways to Quantify Roundoff Error

An easier but less rigorous way of bounding α in Lemma 2.7.1 is to say $| \alpha | \le n \mathbf { u } \mathbf { + } O ( \mathbf { u } ^ { 2 } )$ . With this convention we have

$$
\left| \mathrm{fl} \left(x ^ {T} y\right) - x ^ {T} y \right| \leq n \mathbf {u} \left| x \right| ^ {T} \left| y \right| + O \left(\mathbf {u} ^ {2}\right). \tag {2.7.12}
$$

Other ways of expressing the same result include

$$
\left| \mathbf {f l} \left(x ^ {T} y\right) - x ^ {T} y \right| \leq \phi (n) \mathbf {u} \left| x \right| ^ {T} | y | \tag {2.7.13}
$$

and

$$
\left| \mathbf {f l} \left(x ^ {T} y\right) - x ^ {T} y \right| \leq c n \mathbf {u} | x | ^ {T} | y |, \tag {2.7.14}
$$

where $\phi ( n )$ is a “modest” function of n and c is a constant of order unity.

We shall not express a preference for any of the error bounding styles shown in (2.7.11)–(2.7.14). This spares us the necessity of translating the roundoff results that appear in the literature into a fixed format. Moreover, paying overly close attention to the details of an error bound is inconsistent with the “philosophy” of roundoff analysis. As Wilkinson (1971, p. 567) says,

There is still a tendency to attach too much importance to the precise error bounds obtained by an a priori error analysis. In my opinion, the bound itself is usually the least important part of it. The main object of such an analysis is to expose the potential instabilities, if any, of an algorithm so that hopefully from the insight thus obtained one might be led to improved algorithms. Usually the bound itself is weaker than it might have been because of the necessity of restricting the mass of detail to a reasonable level and because of the limitations imposed by expressing the errors in terms of matrix norms. A priori bounds are not, in general, quantities that should be used in practice. Practical error bounds should usually be determined by some form of a posteriori error analysis, since this takes full advantage of the statistical distribution of rounding errors and of any special features, such as sparseness, in the matrix.

It is important to keep these perspectives in mind.

# 2.7.8 Roundoff in Other Basic Matrix Computations

It is easy to show that if A and B are floating point matrices and α is a floating point number, then

$$
\mathsf {f l} (\alpha A) = \alpha A + E, \quad | E | \leq \mathbf {u} | \alpha A |, \tag {2.7.15}
$$

and

$$
\mathsf {f l} (A + B) = (A + B) + E, \quad | E | \leq \mathbf {u} | A + B |. \tag {2.7.16}
$$

As a consequence of these two results, it is easy to verify that computed saxpy’s and outer product updates satisfy

$$
\operatorname{fl} (y + \alpha x) = y + \alpha x + z, \quad | z | \leq \mathbf {u} (| y | + 2 | \alpha x |) + O \left(\mathbf {u} ^ {2}\right), \tag {2.7.17}
$$

$$
\mathsf {f l} (C + u v ^ {T}) = C + u v ^ {T} + E, \quad | E | \leq \mathbf {u} \left(| C | + 2 | u v ^ {T} |\right) + O \left(\mathbf {u} ^ {2}\right). \tag {2.7.18}
$$

Using (2.7.11) it is easy to show that a dot-product-based multiplication of two floating point matrices A and B satisfies

$$
\mathsf {f l} (A B) = A B + E, \quad | E | \leq n \mathbf {u} | A | | B | + O \left(\mathbf {u} ^ {2}\right). \tag {2.7.19}
$$

The same result applies if a gaxpy or outer product based procedure is used. Notice that matrix multiplication does not necessarily give small relative error since $| A B |$ may be much smaller than $| A | | B |$ , e.g.,

$$
\left[ \begin{array}{c c} 1 & 1 \\ 0 & 0 \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ -. 9 9 & 0 \end{array} \right] = \left[ \begin{array}{c c}. 0 1 & 0 \\ 0 & 0 \end{array} \right].
$$

It is easy to obtain norm bounds from the roundoff results developed thus far. If we look at the 1-norm error in floating point matrix multiplication, then it is easy to show from (2.7.19) that

$$
\left\| \mathfrak {f l} (A B) - A B \right\| _ {1} \leq n \mathbf {u} \| A \| _ {1} \| B \| _ {1} + O \left(\mathbf {u} ^ {2}\right). \tag {2.7.20}
$$

# 2.7.9 Forward and Backward Error Analyses

Each roundoff bound given above is the consequence of a forward error analysis. An alternative style of characterizing the roundoff errors in an algorithm is accomplished through a technique known as backward error analysis. Here, the rounding errors are related to the input data rather than the answer. By way of illustration, consider the n = 2 version of triangular matrix multiplication. It can be shown that:

$$
\mathsf {f l} (A B) = \left[ \begin{array}{c c} a _ {1 1} b _ {1 1} (1 + \epsilon_ {1}) & (a _ {1 1} b _ {1 2} (1 + \epsilon_ {2}) + a _ {1 2} b _ {2 2} (1 + \epsilon_ {3})) (1 + \epsilon_ {4}) \\ 0 & a _ {2 2} b _ {2 2} (1 + \epsilon_ {5}) \end{array} \right]
$$

where $| \epsilon _ { i } | \le \mathbf { u } .$ , for $i = 1 { : } 5$ . However, if we define

$$
\hat {A} = \left[ \begin{array}{c c} a _ {1 1} & a _ {1 2} (1 + \epsilon_ {3}) (1 + \epsilon_ {4}) \\ 0 & a _ {2 2} (1 + \epsilon_ {5}) \end{array} \right]
$$

and

$$
\hat {B} = \left[ \begin{array}{c c} b _ {1 1} (1 + \epsilon_ {1}) & b _ {1 2} (1 + \epsilon_ {2}) (1 + \epsilon_ {4}) \\ 0 & b _ {2 2} \end{array} \right],
$$

then it is easily verified that $\mathsf { f l } ( A B ) = \hat { A } \hat { B }$ . Moreover,

$$
\hat {A} = A + E, \quad | E | \leq 2 \mathbf {u} | A | + O (\mathbf {u} ^ {2}),
$$

$$
\hat {B} = B + F, \quad | F | \leq 2 \mathbf {u} | B | + O (\mathbf {u} ^ {2}).
$$

which shows that the computed product is the exact product of slightly perturbed A and B.

# 2.7.10 Error in Strassen Multiplication

In §1.3.11 we outlined a recursive matrix multiplication procedure due to Strassen. It is instructive to compare the effect of roundoff in this method with the effect of roundoff in any of the conventional matrix multiplication methods of §1.1.

It can be shown that the Strassen approach (Algorithm 1.3.1) produces a $\hat { C } =$ fl(AB) that satisfies an inequality of the form (2.7.20). This is perfectly satisfactory in many applications. However, the $\hat { C }$ that Strassen’s method produces does not always satisfy an inequality of the form (2.7.19). To see this, suppose that

$$
A = B = \left[ \begin{array}{c c}. 9 9 & . 0 0 1 0 \\ . 0 0 1 0 & . 9 9 \end{array} \right]
$$

and that we execute Algorithm 1.3.1 using 2-digit floating point arithmetic. Among other things, the following quantities are computed:

$$
\hat {P} _ {3} = \mathrm{fl} (. 9 9 (. 0 0 1 -. 9 9)) = -. 9 8,
$$

$$
\hat {P} _ {5} = \mathsf {f l} ((. 9 9 +. 0 0 1). 9 9) = . 9 8,
$$

$$
\hat {c} _ {1 2} = \mathsf {f l} (\hat {P} _ {3} + \hat {P} _ {5}) = 0. 0.
$$

In exact arithmetic $c _ { 1 2 } = 2 ( . 0 0 1 ) ( . 9 9 ) = . 0 0 1 9 8$ and thus Algorithm 1.3.1 produces a $\hat { c } _ { 1 2 }$ with no correct significant digits. The Strassen approach gets into trouble in this example because small off-diagonal entries are combined with large diagonal entries. Note that in conventional matrix multiplication the sums $b _ { 1 2 } + b _ { 2 2 }$ and $a _ { 1 1 } + a _ { 1 2 }$ do not arise. For that reason, the contribution of the small off-diagonal elements is not lost in this example. Indeed, for the above A and B a conventional matrix multiplication gives $\hat { c } _ { 1 2 } = . 0 0 2 0$ .

Failure to produce a componentwise accurate $\hat { C }$ can be a serious shortcoming in some applications. For example, in Markov processes the $a _ { i j } , b _ { i j }$ , and $c _ { i j }$ are transition probabilities and are therefore nonnegative. It may be critical to compute $c _ { i j }$ accurately if it reflects a particularly important probability in the modeled phenomenon. Note that if $A \geq 0$ and $B \geq 0$ , then conventional matrix multiplication produces a product $\hat { C }$ that has small componentwise relative error:

$$
| \hat {C} - C | \leq n \mathbf {u} | A | | B | + O \left(\mathbf {u} ^ {2}\right) = n \mathbf {u} | C | + O \left(\mathbf {u} ^ {2}\right).
$$

This follows from (2.7.19). Because we cannot say the same for the Strassen approach, we conclude that Algorithm 1.3.1 is not attractive for certain nonnegative matrix multiplication problems if relatively accurate $\hat { c } _ { i j }$ are required.

Extrapolating from this discussion we reach two fairly obvious but important conclusions:

• Different methods for computing the same quantity can produce substantially different results.   
• Whether or not an algorithm produces satisfactory results depends upon the type of problem solved and the goals of the user.

These observations are clarified in subsequent chapters and are intimately related to the concepts of algorithm stability and problem condition. See §3.4.10.

# 2.7.11 Analysis of an Ideal Equation Solver

A nice way to conclude this chapter and to anticipate the next is to analyze the quality of a “make-believe” $A x = b$ solution process in which all floating point operations are performed exactly except the storage of the matrix A and the right-hand-side b. It follows that the computed solution ˆx satisfies

$$
(A + E) \hat {x} = (b + e), \quad \| E \| _ {\infty} \leq \mathbf {u} \| A \| _ {\infty}, \quad \| e \| _ {\infty} \leq \mathbf {u} \| b \| _ {\infty}. \tag {2.7.21}
$$

where

$$
\mathsf {f l} (b) = b + e, \quad \mathsf {f l} (A) = A + E.
$$

If u $\kappa _ { \infty } ( A ) \leq { \frac { 1 } { 2 } }$ (say), then by Theorem 2.6.2 it can be shown that

$$
\frac {\| x - \hat {x} \| _ {\infty}}{\| x \| _ {\infty}} \leq 4 \mathbf {u} \kappa_ {\infty} (A). \tag {2.7.22}
$$

The bounds (2.7.21) and (2.7.22) are “best possible” norm bounds. No general ∞- norm error analysis of a linear equation solver that requires the storage of A and b can render sharper bounds. As a consequence, we cannot justifiably criticize an algorithm for returning an inaccurate ˆx if A is ill-conditioned relative to the unit roundoff, e.g., u $\kappa _ { \infty } ( A ) \approx 1$ . On the other hand, we have every “right” to pursue the development of a linear equation solver that renders the exact solution to a nearby problem in the style of (2.7.21).

# Problems

P2.7.1 Show that if (2.7.8) is applied with $y = x$ , then $\mathsf { f l } ( x ^ { T } x ) = x ^ { T } x ( 1 + \alpha )$ where $| \alpha | \le n \mathbf { u } + O ( \mathbf { u } ^ { 2 } )$

P2.7.2 Prove (2.7.4) assuming that $\operatorname { f l } ( x )$ is the nearest floating point number to $x \in \mathbb { R }$

P2.7.3 Show that if $E \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ , then $\| \mathbf { \epsilon } | \| _ { 2 } \leq { \sqrt { n } } \| E \| _ { 2 }$ . This result is useful when deriving norm bounds from absolute value bounds.

P2.7.4 Assume the existence of a square root function satisfying ${ \mathfrak { f l } } ( { \sqrt { x } } ) = { \sqrt { x } } ( 1 + \epsilon )$ with $| \epsilon | \le \mathbf { u }$ Give an algorithm for computing $\parallel x \parallel _ { 2 }$ and bound the rounding errors.

P2.7.5 Suppose A and B are n-by-n upper triangular floating point matrices. If $\hat { C } = \mathfrak { f l } ( A B )$ is computed using one of the conventional §1.1 algorithms, does it follow that $\hat { C } = \hat { A } \hat { B }$ where Aˆ and Bˆ are close to A and B?

P2.7.6 Suppose A and B are n-by-n floating point matrices and that $\parallel | A ^ { - 1 } | | A | \parallel _ { \infty } = \tau .$ . Show that if $\hat { C } = \mathsf { f l } ( A B )$ is obtained using any of the §1.1 algorithms, then there exists a Bˆ so that $\hat { C } = A \hat { B }$ and $\| \hat { B } - B \| _ { \infty } \leq n \mathbf { u } \tau \| B \| _ { \infty } + O ( \mathbf { u } ^ { 2 } )$ .

P2.7.7 Prove (2.7.19).

P2.7.8 For the IEEE double format, what is the largest power of 10 that can be represented exactly? What is the largest integer that can be represented exactly?

P2.7.9 For $k = 1 { : } 6 2$ , what is the largest power of 10 that can be stored exactly if k bits are are allocated for the mantissa and $6 3 - k$ are allocated for the exponent?

P2.7.10 Consider the quadratic equation

$$
q (\lambda) = \det \left(\left[ \begin{array}{c c} w - \lambda & x \\ x & z - \lambda \end{array} \right]\right).
$$

This quadratic has two real roots $r _ { 1 }$ and $r _ { 2 }$ . Assume that $| r _ { 1 } - z | \leq | r _ { 2 } - z |$ . Give an algorithm that computes $r _ { 1 }$ to full machine precision.

# Notes and References for §2.7

For an excellent, comprehensive treatment of IEEE arithmetic and its implications, see:

M.L. Overton (2001). Numerical Computing with IEEE Arithmetic, SIAM Publications, Philadelphia, PA.

The following basic references are notable for the floating point insights that they offer: Wilkinson (AEP), Stewart (IMC), Higham (ASNA), and Demmel (ANLA). For high-level perspectives we recommend:

J.H. Wilkinson (1963). Rounding Errors in Algebraic Processes, Prentice-Hall, Englewood Cliffs, NJ.

G.E. Forsythe (1970). “Pitfalls in Computation or Why a Math Book is Not Enough,” Amer. Math. Monthly 77, 931–956.

J.H. Wilkinson (1971). “Modern Error Analysis,” SIAM Review 13, 548–68.

U.W. Kulisch and W.L. Miranker (1986). “The Arithmetic of the Digital Computer,” SIAM Review 28, 1–40.

F. Chaitin-Chatelin and V. Frayse´e (1996). Lectures on Finite Precision Computations, SIAM Publications, Philadelphia, PA.

The design of production software for matrix computations requires a detailed understanding of finite precision arithmetic, see:

J.W. Demmel (1984). “Underflow and the Reliability of Numerical Software,” SIAM J. Sci. Stat. Comput. 5, 887–919.

W.J. Cody (1988). “ALGORITHM 665 MACHAR: A Subroutine to Dynamically Determine Machine Parameters,” ACM Trans. Math. Softw. 14, 303–311.

D. Goldberg (1991). “What Every Computer Scientist Should Know About Floating Point Arithmetic,” ACM Surveys 23, 5–48.

Other developments in error analysis involve interval analysis, the building of statistical models of roundoff error, and the automating of the analysis itself:

J. Larson and A. Sameh (1978). “Efficient Calculation of the Effects of Roundoff Errors,” ACM Trans. Math. Softw. 4, 228–36.

W. Miller and D. Spooner (1978). “Software for Roundoff Analysis, II,” ACM Trans. Math. Softw. 4, 369–90.

R.E. Moore (1979). Methods and Applications of Interval Analysis, SIAM Publications, Philadelphia, PA.

J.M. Yohe (1979). “Software for Interval Arithmetic: A Reasonable Portable Package,” ACM Trans. Math. Softw. 5, 50–63.

The accuracy of floating point summation is detailed in:

S.M. Rump, T. Ogita, and S. Oishi (2008). “Accurate Floating-Point Summation Part I: Faithful Rounding,” SIAM J. Sci. Comput. 31, 189–224.

S.M. Rump, T. Ogita, and S. Oishi (2008). “Accurate Floating-Point Summation Part II: Sign, K-fold Faithful and Rounding to Nearest,” SIAM J. Sci. Comput. 31, 1269–1302.   
For an analysis of the Strassen algorithm and other “fast” linear algebra procedures, see:   
R.P. Brent (1970). “Error Analysis of Algorithms for Matrix Multiplication and Triangular Decomposition Using Winograd’s Identity,” Numer. Math. 16, 145–156.   
W. Miller (1975). “Computational Complexity and Numerical Stability,” SIAM J. Comput. 4, 97–107.   
N.J. Higham (1992). “Stability of a Method for Multiplying Complex Matrices with Three Real Matrix Multiplications,” SIAM J. Matrix Anal. Applic. 13, 681–687.   
J.W. Demmel and N.J. Higham (1992). “Stability of Block Algorithms with Fast Level-3 BLAS,” ACM Trans. Math. Softw. 18, 274–291.   
B. Dumitrescu (1998). “Improving and Estimating the Accuracy of Strassen’s Algorithm,” Numer. Math. 79, 485–499.   
The issue of extended precision has received considerable attention. For example, a superaccurate dot product results if the summation can be accumulated in a register that is “twice as wide” as the floating representation of vector components. The overhead may be tolerable in a given algorithm if extended precision is needed in only a few critical steps. For insights into this topic, see:   
R.P. Brent (1978). “A Fortran Multiple Precision Arithmetic Package,” ACM Trans. Math. Softw. 4, 57–70.   
R.P. Brent (1978). “Algorithm 524 MP, a Fortran Multiple Precision Arithmetic Package,” ACM Trans. Math. Softw. 4, 71–81.   
D.H. Bailey (1993). “Algorithm 719: Multiprecision Translation and Execution of FORTRAN Programs,” ACM Trans. Math. Softw. 19, 288–319.   
X.S. Li, J.W. Demmel, D.H. Bailey, G. Henry, Y. Hida, J. Iskandar, W. Kahan, S.Y. Kang, A. Kapur, M.C. Martin, B.J. Thompson, T. Tung, and D.J. Yoo (2002). “Design, Implementation and Testing of Extended and Mixed Precision BLAS,” ACM Trans. Math. Softw. 28, 152–205.   
J.W. Demmel and Y. Hida (2004). “Accurate and Efficient Floating Point Summation,” SIAM J. Sci. Comput. 25, 1214–1248.   
M. Baboulin, A. Buttari, J. Dongarra, J. Kurzak, J. Langou, J. Langou, P. Luszczek, and S. Tomov (2009). “Accelerating Scientific Computations with Mixed Precision Algorithms,” Comput. Phys. Commun. 180, 2526–2533.
