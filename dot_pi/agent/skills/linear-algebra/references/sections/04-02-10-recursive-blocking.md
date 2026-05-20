# 4.2.10 Recursive Blocking

It is instructive to look a little more deeply into the implementation of a block Cholesky factorization as it is an occasion to stress the importance of designing data structures that are tailored to the problem at hand. High-performance matrix computations are filled with tensions and tradeoffs. For example, a successful pivot strategy might balance concerns about stability and memory traffic. Another tension is between performance and memory constraints. As an example of this, we consider how to achieve level-3 performance in a Cholesky implementation given that the matrix is represented in packed format. This data structure houses the lower (or upper) triangular portion of a matrix $A \in \mathbb { R } ^ { n \times n }$ in a vector of length $N = n ( n + 1 ) / 2$ . The symvec arrangement stacks the lower triangular subcolumns, e.g.,

$$
\operatorname{symvec} (A) = \left[ a _ {1 1} a _ {2 1} a _ {3 1} a _ {4 1} a _ {2 2} a _ {3 2} a _ {4 2} a _ {3 3} a _ {4 3} a _ {4 4} \right] ^ {T}. \tag {4.2.20}
$$

This layout is not very friendly when it comes to block Cholesky calculations because the assembly of an A-block (say $A ( i _ { 1 } { : } i _ { 2 } , j _ { 1 } { : } j _ { 2 } ) )$ involves irregular memory access patterns. To realize a high-performance matrix multiplication it is usually necessary to have the matrices laid out conventionally as full rectangular arrays that are contiguous in memory, e.g.,

$$
\operatorname{vec} (A) = \left[ a _ {1 1} a _ {2 1} a _ {3 1} a _ {4 1} a _ {1 2} a _ {2 2} a _ {3 2} a _ {4 2} a _ {1 3} a _ {2 3} a _ {3 3} a _ {4 3} a _ {1 4} a _ {2 4} a _ {3 4} a _ {4 4} \right] ^ {T}. \tag {4.2.21}
$$

(Recall that we introduced the vec operation in §1.3.7.) Thus, the challenge is to develop a high performance block algorithm that overwrites a symmetric positive definite A in packed format with its Cholesky factor G in packed format. Toward that end, we present the main ideas behind a recursive data structure that supports level-3 computation and is storage efficient. As memory hierarchies get deeper and more complex, recursive data structures are an interesting way to address the problem of blocking for performance.

The starting point is once again a 2-by-2 blocking of the equation $A = G G ^ { T }$ :

$$
{\left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right]} = {\left[ \begin{array}{l l} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right]} {\left[ \begin{array}{l l} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right]} ^ {T}.
$$

However, unlike in (4.2.18) where $A _ { 1 1 }$ has a chosen block size, we now assume that $A _ { 1 1 } \in \mathbb { R } ^ { m \times m }$ where $m = \operatorname { c e i l } ( n / 2 )$ . In other words, the four blocks are roughly the same size. As before, we equate entries and identify the key subcomputations:

<table><tr><td> $G_{11}G_{11}^{T} = A_{11}$ </td><td>half-sized Cholesky.</td></tr><tr><td> $G_{21}G_{11}^{T} = A_{21}$ </td><td>multiple-right-hand-side triangular solve.</td></tr><tr><td> $\tilde{A}_{22} = A_{22} - G_{21}G_{21}^{T}$ </td><td>symmetric matrix multiplication update.</td></tr><tr><td> $G_{22}G_{22}^{T} = \tilde{A}_{22}$ </td><td>half-sized Cholesky.</td></tr></table>

Our goal is to develop a symmetry-exploiting, level-3-rich procedure that overwrites A with its Cholesky factor G. To do this we introduce the mixed packed format. An n = 9 example with A11 ∈ IR5×5 $n = 9$ $A _ { 1 1 } \in \mathbb { R } ^ { 5 \times 5 }$ serves to distinguish this layout from the conventional packed format layout:

<table><tr><td>1</td><td></td><td></td><td></td><td></td></tr><tr><td>2 10</td><td></td><td></td><td></td><td></td></tr><tr><td>3 11 18</td><td></td><td></td><td></td><td></td></tr><tr><td>4 12 19 25</td><td></td><td></td><td></td><td></td></tr><tr><td>5 13 20 26 31</td><td></td><td></td><td></td><td></td></tr><tr><td>6 14 21 27 32</td><td></td><td></td><td></td><td>36</td></tr><tr><td>7 15 22 28 33</td><td></td><td></td><td></td><td>37 40</td></tr><tr><td>8 16 23 29 34</td><td></td><td></td><td></td><td>38 41 43</td></tr><tr><td>9 17 24 30 35</td><td></td><td></td><td></td><td>39 42 44 45</td></tr></table>

Packed format

<table><tr><td>1</td><td></td><td></td><td></td><td></td><td></td></tr><tr><td>2</td><td>6</td><td></td><td></td><td></td><td></td></tr><tr><td>3</td><td> $\overline{t}$ </td><td>10</td><td></td><td></td><td></td></tr><tr><td>4</td><td>8</td><td>11</td><td>13</td><td></td><td></td></tr><tr><td>5</td><td>9</td><td>12</td><td>14</td><td>15</td><td></td></tr><tr><td>16</td><td>20</td><td>24</td><td>28</td><td>32</td><td>36</td></tr><tr><td>17</td><td>21</td><td>25</td><td>29</td><td>33</td><td>37 40</td></tr><tr><td>18</td><td>22</td><td>26</td><td>30</td><td>34</td><td>38 41 43</td></tr><tr><td>19</td><td>23</td><td>27</td><td>31</td><td>35</td><td>39 42 44 45</td></tr></table>

Mixed packed format

Notice how the entries from $A _ { 1 1 }$ and $A _ { 2 1 }$ are shuffled with the conventional packed format layout. On the other hand, with the mixed packed format layout, the 15 entries that define $A _ { 1 1 }$ are followed by the 20 numbers that define $A _ { 2 1 }$ which in turn are followed by the 10 numbers that define $A _ { 2 2 }$ . The process can be repeated on $A _ { 1 1 }$ and

<table><tr><td colspan="3">12 43 5 6</td><td colspan="2"></td><td rowspan="2" colspan="2"></td></tr><tr><td colspan="3">7 9 118 10 12</td><td colspan="2">1314 15</td></tr><tr><td rowspan="3" colspan="3">16 20 2417 21 2518 22 2619 23 27</td><td rowspan="3" colspan="2">28 3229 3330 3431 35</td><td colspan="2">3637 38</td></tr><tr><td rowspan="2" colspan="2">39 4140 42</td></tr><tr></tr></table>

Thus, the key to this recursively defined data layout is the idea of representing square diagonal blocks in a mixed packed format. To be precise, recall the definition of vec and symvec in (4.2.20) and (4.2.21). If $C \in \mathbb { R } ^ { q \times q }$ is such a block, then

$$
\operatorname{mixvec} (C) = \left[ \begin{array}{c} \operatorname{symvec} \left(C _ {1 1}\right) \\ \operatorname{vec} \left(C _ {2 1}\right) \\ \operatorname{symvec} \left(C _ {2 2}\right) \end{array} \right] \tag {4.2.22}
$$

where $m = \mathrm { c e i l } ( q / 2 ) , C _ { 1 1 } = C ( 1 { : } m , 1 { : } m ) , C _ { 2 2 } = C ( m + 1 { : } n , m + 1 { : } n )$ , and $C _ { 2 1 } =$ $C ( m + 1 { : } n , 1 { : } m )$ . Notice that since $C _ { 2 1 }$ is conventionally stored, it is ready to be engaged in a high-performance matrix multiplication.

We now outline a recursive, divide-and-conquer block Cholesky procedure that works with A in packed format. To achieve high performance the incoming A is converted to mixed format at each level of the recursion. Assuming the existence of a triangular system solve procedure TriSol (for the system $G _ { 2 1 } G _ { 1 1 } ^ { T } = A _ { 2 1 } )$ and a symmetric update procedure SymUpdate (for $A _ { 2 2 }  A _ { 2 2 } - G _ { 2 1 } G _ { 2 1 } ^ { T } )$ we have the following framework:

function G = PackedBlockCholesky(A)

{A and G in packed format}

$$
n = \operatorname{size} (A)
$$

if $n \leq n$ min

G is obtained via any level-2, packed-format Cholesky method .

else

Set $m = { \mathsf { c e i l } } ( n / 2 )$ and overwrite A’s packed-format representation with its mixed-format representation.

$$
G _ {1 1} = \text { PackedBlockCholesky } (A _ {1 1})
$$

$$
G _ {2 1} = \operatorname{TriSol} \left(G _ {1 1}, A _ {2 1}\right)
$$

$$
A _ {2 2} = \text { SymUpdate } (A _ {2 2}, G _ {2 1})
$$

$$
G _ {2 2} = \text { PackedBlockCholesky } (A _ {2 2})
$$

end

Here, $n _ { \mathrm { m i n } }$ is a threshold dimension below which it is not possible to achieve level-3 performance. To take full advantage of the mixed format, the procedures TriSol and SymUpdate require a recursive design based on blockings that halve problem size. For example, TriSol should take the incoming packed format $A _ { 1 1 }$ , convert it to mixed format, and solve a 2-by-2 blocked system of the form

$$
\left[ \begin{array}{c c} X _ {1} & X _ {2} \end{array} \right] \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \end{array} \right] ^ {T} = \left[ \begin{array}{c c} B _ {1} & B _ {2} \end{array} \right].
$$

This sets up a recursive solution based on the half-sized problems

$$
\begin{array}{l} X _ {1} L _ {1 1} ^ {T} = B _ {1}, \\ X _ {2} L _ {2 2} ^ {T} = B _ {2} - X _ {1} L _ {2 1} ^ {T}. \\ \end{array}
$$

Likewise, SymUpdate should take the incoming packed format $A _ { 2 2 }$ , convert it to mixed format, and block the required update as follows:

$$
\left[ \begin{array}{c c} C _ {1 1} & C _ {2 1} ^ {T} \\ C _ {2 1} & C _ {2 2} \end{array} \right] = \left[ \begin{array}{c c} C _ {1 1} & C _ {2 1} ^ {T} \\ C _ {2 1} & C _ {2 2} \end{array} \right] - \left[ \begin{array}{c} Y _ {1} \\ Y _ {2} \end{array} \right] \left[ \begin{array}{c} Y _ {1} \\ Y _ {2} \end{array} \right] ^ {T}.
$$

The evaluation is recursive and based on the half-sized updates

$$
\begin{array}{l} C _ {1 1} = C _ {1 1} - Y _ {1} Y _ {1} ^ {T}, \\ C _ {2 1} = C _ {2 1} - Y _ {2} Y _ {1} ^ {T}, \\ C _ {2 2} = C _ {2 2} - Y _ {2} Y _ {2} ^ {T}. \\ \end{array}
$$

Of course, if the incoming matrices are small enough relative to $n _ { \mathrm { m i n } }$ , then TriSol and SymUpdate carry out their tasks conventionally without any further subdivisions.

Overall, it can be shown that PackedBlockCholesky has a level-3 fraction approximately equal to $1 - O ( n _ { \operatorname* { m i n } } / n )$ .

# Problems

P4.2.1 Suppose that $H = A + i B$ is Hermitian and positive definite with A, $B \in \mathbb { R } ^ { n \times n }$ . This means that $x ^ { H } H x > 0$ whenever x 
= 0. (a) Show that

$$
C = \left[ \begin{array}{c c} A & - B \\ B & A \end{array} \right]
$$

is symmetric and positive definite. (b) Formulate an algorithm for solving $( A + i B ) ( x + i y ) = ( b + i c )$ , where b, c, x, and y are in $\mathbb { R } ^ { n }$ . It should involve $8 n ^ { 3 } / 3$ flops. How much storage is required?

P4.2.2 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and positive definite. Give an algorithm for computing an upper triangular matrix $R \in \mathbb { R } ^ { n \times n }$ such that $\dot { A } = R R ^ { T }$ .

P4.2.3 Let $A \in \mathbb { R } ^ { n \times n }$ be positive definite and set $T = ( A + A ^ { T } ) / 2$ and $S = ( A - A ^ { T } ) / 2$ . (a) Show that $\parallel A ^ { - 1 } \parallel _ { 2 } \leq \parallel T ^ { - 1 } \parallel _ { 2 }$ and $x ^ { T } A ^ { - 1 } x \leq x ^ { T } T ^ { - 1 } x$ for all $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . (b) Show that if $A = L D M ^ { T }$ , then $d _ { k } \geq 1 / \Vert \ T ^ { - 1 } \Vert _ { 2 }$ for $k = 1 { : } n$ .

P4.2.4 Find a 2-by-2 real matrix A with the property that $x ^ { T } A x > 0$ for all real nonzero 2-vectors but which is not positive definite when regarded as a member of C2×2 . $\mathbb { C } ^ { 2 \times 2 }$

P4.2.5 Suppose $A \in \mathbb { R } ^ { n \times n }$ has a positive diagonal. Show that if both A and $A ^ { T }$ are strictly diagonally

dominant, then A is positive definite.

P4.2.6 Show that the function $f ( x ) = \sqrt { x ^ { T } A x } / 2$ is a vector norm on $\mathbb { R } ^ { n }$ if and only if A is positive definite.

P4.2.7 Modify Algorithm 4.2.1 so that if the square root of a negative number is encountered, then the algorithm finds a unit vector x so that $x ^ { T } A x < 0$ and terminates.

P4.2.8 Develop an outer product implementation of Algorithm 4.2.1 and a gaxpy implementation of Algorithm 4.2.2.

P4.2.9 Assume that $A \in \mathbb { C } ^ { n \times n }$ is Hermitian and positive definite. Show that if $a _ { 1 1 } = \cdots = a _ { n n } = 1$ and $| a _ { i j } | < 1$ for all $i \neq j$ , then $\mathrm { d i a g } ( A ^ { - 1 } ) \geq \mathrm { d i a g } ( ( \mathsf { R e } ( A ) ) ^ { - 1 } )$ .

P4.2.10 Suppose $A = I + u u ^ { T }$ where $A \in \mathbb { R } ^ { n \times n }$ and $\parallel u \parallel _ { 2 } = 1$ . Give explicit formulae for the diagonal and subdiagonal of A’s Cholesky factor.

P4.2.11 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric positive definite and that its Cholesky factor is available. Let $e _ { k } = I _ { n } ( : , k )$ . For $1 \leq i < j \leq n$ , let $\alpha _ { i j }$ be the smallest real that makes $\mathbf { \bar { \Phi } } A + \alpha ( e _ { i } e _ { j } ^ { T } + e _ { j } e _ { i } ^ { T } )$ singular. Likewise, let $\alpha _ { i i }$ be the smallest real that makes $( A + \alpha e _ { i } e _ { i } ^ { T } )$ singular. Show how to compute these quantities using the Sherman-Morrison-Woodbury formula. How many flops are required to find all the $\alpha _ { i j } ?$

P4.2.12 Show that if

$$
M = \left[ \begin{array}{c c} A & B \\ B ^ {T} & C \end{array} \right]
$$

is symmetric positive definite and A and C are square, then

$$
M ^ {- 1} = \left[ \begin{array}{c c} A ^ {- 1} + A ^ {- 1} B S ^ {- 1} B ^ {T} A ^ {- 1} & - A ^ {- 1} B S ^ {- 1} \\ S ^ {- 1} B ^ {T} A ^ {- 1} & S ^ {- 1} \end{array} \right], \qquad S = C - B ^ {T} A ^ {- 1} B.
$$

P4.2.13 Suppose $\sigma \in \mathbb { R }$ and $u \in \mathbb { R } ^ { n }$ . Under what conditions can we find a matrix $X \in \mathbb { R } ^ { n \times n }$ so that $X ( I + \sigma u u ^ { T } ) X = I _ { n } ?$ Give an efficient algorithm for computing X if it exists.

P4.2.14 Suppose $D = \operatorname { d i a g } ( d _ { 1 } , \dotsc , d _ { n } )$ with $d _ { i } > 0$ for all i. Give an efficient algorithm for computing the largest entry in the matrix $( D + C C ^ { T } ) ^ { - 1 }$ where $C \in \mathbb { R } ^ { n \times r }$ . Hint: Use the Sherman-Morrison-Woodbury formula.

P4.2.15 Suppose $A ( \lambda )$ has continuously differentiable entries and is always symmetric and positive definite. If $f ( \lambda ) = \log ( { \mathsf { d e t } } ( A ( \lambda ) ) )$ , then how would you compute $f ^ { \prime } ( 0 ) \smash { \ ? }$

P4.2.16 Suppose $A \in \mathbb { R } ^ { n \times n }$ is a rank-r symmetric positive semidefinite matrix. Assume that it costs one dollar to evaluate each $\boldsymbol { a } _ { i j }$ . Show how to compute the factorization (4.2.17) spending only $O ( n r )$ dollars on $a _ { i j }$ evaluation.

P4.2.17 The point of this problem is to show that from the complexity point of view, if you have a fast matrix multiplication algorithm, then you have an equally fast matrix inversion algorithm, and vice versa. (a) Suppose $F _ { n }$ is the number of flops required by some method to form the inverse of an n-by-n matrix. Assume that there exists a constant $c _ { 1 }$ and a real number α such that $F _ { n } \leq c _ { 1 } n ^ { \alpha }$ for all n. Show that there is a method that can compute the $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ matrix product $A B$ with fewer than $c _ { 2 } n ^ { \alpha }$ flops where $c _ { 2 }$ is a constant independent of n. Hint: Consider the inverse of

$$
C = \left[ \begin{array}{c c c} I _ {n} & A & 0 \\ 0 & I _ {n} & B \\ 0 & 0 & I _ {n} \end{array} \right].
$$

(b) Let $G _ { n }$ be the number of flops required by some method to form the $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ matrix product $A B$ . Assume that there exists a constant $c _ { 1 }$ and a real number α such that $G _ { n } \leq c _ { 1 } n ^ { \alpha }$ for all n. Show that there is a method that can invert a nonsingular n-by-n matrix A with fewer than $c _ { 2 } n ^ { \alpha }$ flops where $c _ { 2 }$ is a constant. Hint: First show that the result applies for triangular matrices by applying recursion to

$$
\left[ \begin{array}{c c} G _ {1 1} & 0 \\ G _ {2 1} & G _ {2 2} \end{array} \right] ^ {- 1} = \left[ \begin{array}{c c} G _ {1 1} ^ {- 1} & 0 \\ - G _ {2 2} ^ {- 1} G _ {2 1} G _ {1 1} ^ {- 1} & G _ {2 2} ^ {- 1} \end{array} \right].
$$

Then observe that for general A, $A ^ { - 1 } = A ^ { T } ( A A ^ { T } ) ^ { - 1 } = A ^ { T } G ^ { - T } G ^ { - 1 }$ where $A A ^ { T } = G G ^ { T }$ is the Cholesky factorization.

# Notes and References for §4.2

For an in-depth theoretical treatment of positive definiteness, see:

R. Bhatia (2007). Positive Definite Matrices, Princeton University Press, Princeton, NJ.

The definiteness of the quadratic form $x ^ { T }$ Ax can frequently be established by considering the mathematics of the underlying problem. For example, the discretization of certain partial differential operators gives rise to provably positive definite matrices. Aspects of the unsymmetric positive definite problem are discussed in:

A. Buckley (1974). “A Note on Matrices A = I + H, H Skew-Symmetric,” Z. Angew. Math. Mech. 54, 125–126.   
A. Buckley (1977). “On the Solution of Certain Skew-Symmetric Linear Systems,” SIAM J. Numer. Anal. 14, 566–570.   
G.H. Golub and C. Van Loan (1979). “Unsymmetric Positive Definite Linear Systems,” Lin. Alg. Applic. 28, 85–98.   
R. Mathias (1992). “Matrices with Positive Definite Hermitian Part: Inequalities and Linear Systems,” SIAM J. Matrix Anal. Applic. 13, 640–654.   
K.D. Ikramov and A.B. Kucherov (2000). “Bounding the growth factor in Gaussian elimination for Buckley’s class of complex symmetric matrices,” Numer. Lin. Alg. 7, 269–274.

Complex symmetric matrices have the property that their real and imaginary parts are each symmetric. The following paper shows that if they are also positive definite, then the $\mathrm { \dot { L } D L ^ { T } }$ factorization is safe to compute without pivoting:

S. Serbin (1980). “On Factoring a Class of Complex Symmetric Matrices Without Pivoting,” Math. Comput. 35, 1231–1234.

Historically important Algol implementations of the Cholesky factorization include:

R.S. Martin, G. Peters, and J.H. Wilkinson (1965). “Symmetric Decomposition of a Positive Definite Matrix,” Numer. Math. 7, 362–83.   
R.S. Martin, G. Peters, and J.H. Wilkinson (1966). “Iterative Refinement of the Solution of a Positive Definite System of Equations,” Numer. Math. 8, 203–16.   
F.L. Bauer and C. Reinsch (1971). “Inversion of Positive Definite Matrices by the Gauss-Jordan Method,” in Handbook for Automatic Computation Vol. 2, Linear Algebra, J.H. Wilkinson and C. Reinsch (eds.), Springer-Verlag, New York, 45–49.

For roundoff error analysis of Cholesky, see:

J.H. Wilkinson (1968). “A Priori Error Analysis of Algebraic Processes,” Proceedings of the International Congress on Mathematics, Izdat. Mir, 1968, Moscow, 629–39.   
J. Meinguet (1983). “Refined Error Analyses of Cholesky Factorization,” SIAM J. Numer. Anal. 20, 1243–1250.   
A. Kielbasinski (1987). “A Note on Rounding Error Analysis of Cholesky Factorization,” Lin. Alg. Applic. 88/89, 487–494.   
N.J. Higham (1990). “Analysis of the Cholesky Decomposition of a Semidefinite Matrix,” in Reliable Numerical Computation, M.G. Cox and S.J. Hammarling (eds.), Oxford University Press, Oxford, U.K., 161–185.   
J-Guang Sun (1992). “Rounding Error and Perturbation Bounds for the Cholesky and $L D L ^ { T }$ Factorizations,” Lin. Alg. Applic. 173, 77–97.

The floating point determination of positive definiteness is an interesting problem, see:

S.M. Rump (2006). “Verification of Positive Definiteness,” BIT 46, 433–452.

The question of how the Cholesky triangle G changes when $A = G G ^ { T }$ is perturbed is analyzed in:

G.W. Stewart (1977). “Perturbation Bounds for the QR Factorization of a Matrix,” SIAM J. Num. Anal. 14, 509–18.

Z. Dram˘ac, M. Omladi˘c, and K. Veseli˘c (1994). “On the Perturbation of the Cholesky Factorization,” SIAM J. Matrix Anal. Applic. 15, 1319–1332.

X-W. Chang, C.C. Paige, and G.W. Stewart (1996). “New Perturbation Analyses for the Cholesky Factorization,” IMA J. Numer. Anal. 16, 457–484.


---

<!-- golub_200_249 -->

G.W. Stewart (1997) “On the Perturbation of LU and Cholesky Factors,” IMA J. Numer. Anal. 17, 1–6.   
Nearness/sensitivity issues associated with positive semidefiniteness are presented in:   
N.J. Higham (1988). “Computing a Nearest Symmetric Positive Semidefinite Matrix,” Lin. Alg. Applic. 103, 103–118.   
The numerical issues associated with semi-definite rank determination are covered in:   
P.C. Hansen and P.Y. Yalamov (2001). “Computing Symmetric Rank-Revealing Decompositions via Triangular Factorization,” SIAM J. Matrix Anal. Applic. 23, 443–458.   
M. Gu and L. Miranian (2004). “Strong Rank-Revealing Cholesky Factorization,” ETNA 17, 76–92.   
The issues that surround level-3 performance of packed-format Cholesky are discussed in:   
F.G. Gustavson (1997). “Recursion Leads to Automatic Variable Blocking for Dense Linear-Algebra Algorithms,” IBM J. Res. Dev. 41, 737–756.   
F.G. Gustavson, A. Henriksson, I. Jonsson, B. K˚agstr¨om, , and P. Ling (1998). “Recursive Blocked Data Formats and BLAS’s for Dense Linear Algebra Algorithms,” Applied Parallel Computing Large Scale Scientific and Industrial Problems, Lecture Notes in Computer Science, Springer-Verlag, 1541/1998, 195–206.   
F.G. Gustavson and I. Jonsson (2000). “Minimal Storage High-Performance Cholesky Factorization via Blocking and Recursion,” IBM J. Res. Dev. 44, 823–849.   
B.S. Andersen, J. Wasniewski, and F.G. Gustavson (2001). “A Recursive Formulation of Cholesky Factorization of a Matrix in Packed Storage,” ACM Trans. Math. Softw. 27, 214–244.   
E. Elmroth, F. Gustavson, I. Jonsson, and B. K˚agstr¨om, (2004). “Recursive Blocked Algorithms and Hybrid Data Structures for Dense Matrix Library Software,” SIAM Review 46, 3–45.   
F.G. Gustavson, J. Wasniewski, J.J. Dongarra, and J. Langou (2010). “Rectangular Full Packed Format for Cholesky’s Algorithm: Factorization, Solution, and Inversion,” ACM Trans. Math. Softw. 37, Article 19.   
Other high-performance Cholesky implementations include:   
F.G. Gustavson, L. Karlsson, and B. K˚agstr¨om, (2009). “Distributed SBP Cholesky Factorization Algorithms with Near-Optimal Scheduling,” ACM Trans. Math. Softw. 36, Article 11.   
G. Ballard, J. Demmel, O. Holtz, and O. Schwartz (2010). “Communication-Optimal Parallel and Sequential Cholesky,” SIAM J. Sci. Comput. 32, 3495–3523.   
P. Bientinesi, B. Gunter, and R.A. van de Geijn (2008). “Families of Algorithms Related to the Inversion of a Symmetric Positive Definite Matrix,” ACM Trans. Math. Softw. 35, Article 3.   
M.D. Petkovi´c and P.S. Stanimirovi´c (2009). “Generalized Matrix Inversion is not Harder than Matrix Multiplication,” J. Comput. Appl. Math. 230, 270–282.
