# 5.3 The Full-Rank Least Squares Problem

Consider the problem of finding a vector $\boldsymbol { x } \in \mathbb { R } ^ { n }$ such that $A x = b$ where the data matrix $A \in \mathbb { R } ^ { m \times n }$ and the observation vector $b \in \mathbb { R } ^ { m }$ are given and $m \geq n$ . When there are more equations than unknowns, we say that the system $A x = b$ is overdetermined. Usually an overdetermined system has no exact solution since b must be an element of ran(A), a proper subspace of $\mathbb { R } ^ { m }$ .

This suggests that we strive to minimize $\| A x - b \| _ { p }$ for some suitable choice of $p .$ Different norms render different optimum solutions. For example, if $A = [ 1 , 1 , 1 ] ^ { T }$ and $\boldsymbol { b } = [ b _ { 1 } , b _ { 2 } , b _ { 3 } ] ^ { T }$ with $b _ { 1 } \geq b _ { 2 } \geq b _ { 3 } \geq 0$ , then it can be verified that

$$
p = 1 \Rightarrow x _ {\text {opt}} = b _ {2},
$$

$$
p = 2 \Rightarrow x _ {\text {opt}} = (b _ {1} + b _ {2} + b _ {3}) / 3,
$$

$$
p = \infty \Rightarrow x _ {\text {opt}} = (b _ {1} + b _ {3}) / 2.
$$

Minimization in the 1-norm and infinity-norm is complicated by the fact that the function $f ( x ) = \parallel A x - b \parallel _ { p }$ is not differentiable for these values of $p .$ However, there are several good techniques available for 1-norm and ∞-norm minimization. See Coleman and Li (1992), Li (1993), and Zhang (1993).

In contrast to general p-norm minimization, the least squares (LS) problem

$$
\min _ {x \in \mathbb {R} ^ {n}} \| A x - b \| _ {2} \tag {5.3.1}
$$

is more tractable for two reasons:

• $\begin{array} { r } { \phi ( x ) = \frac { 1 } { 2 } \| A x - b \| _ { 2 } ^ { 2 } } \end{array}$ is a differentiable function of x and so the minimizers of $\phi$ satisfy the gradient equation $\nabla \phi ( x ) = 0$ . This turns out to be an easily constructed symmetric linear system which is positive definite if A has full column rank.   
• The 2-norm is preserved under orthogonal transformation. This means that we can seek an orthogonal Q such that the equivalent problem of minimizing $\parallel ( Q ^ { T } A ) x - ( Q ^ { T } b ) \parallel _ { 2 }$ is “easy” to solve.

In this section we pursue these two solution approaches for the case when A has full column rank. Methods based on normal equations and the QR factorization are detailed and compared.

# 5.3.1 Implications of Full Rank

Suppose $x \in \mathbb { R } ^ { n } , z \in \mathbb { R } ^ { n } , \alpha \in \mathbb { R }$ , and consider the equality

$$
\| A (x + \alpha z) - b \| _ {2} ^ {2} = \| A x - b \| _ {2} ^ {2} + 2 \alpha z ^ {T} A ^ {T} (A x - b) + \alpha^ {2} \| A z \| _ {2} ^ {2}
$$

where $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ . If x solves the LS problem (5.3.1), then we must have $A ^ { T } ( A x - b ) = 0 .$ . Otherwise, i $\mathrm { ~ f ~ } z = - A ^ { T } ( A x - b )$ and we make α small enough, then we obtain the contradictory inequality $\parallel A ( x + \alpha z ) - b \parallel _ { 2 } < \parallel A x - b \parallel _ { 2 }$ . We may also conclude that if x and $x + \alpha z$ are LS minimizers, then $z \in \mathsf { n u l l } ( A )$ .

Thus, if A has full column rank, then there is a unique LS solution $x _ { \mathrm { L S } }$ and it solves the symmetric positive definite linear system

$$
A ^ {T} A x _ {\mathrm{LS}} = A ^ {T} b.
$$

These are called the normal equations. Note that if

$$
\phi (x) = \frac {1}{2} \| A x - b \| _ {2} ^ {2},
$$

then

$$
\nabla \phi (x) = A ^ {T} (A x - b),
$$

so solving the normal equations is tantamount to solving the gradient equation $\nabla \phi = 0$ . We call

$$
r _ {\mathrm{LS}} = b - A x _ {\mathrm{LS}}
$$

the minimum residual and we use the notation

$$
\rho_ {\mathrm{LS}} = \left\| A x _ {\mathrm{LS}} - b \right\| _ {2}
$$

to denote its size. Note that if $\rho _ { \mathrm { L S } }$ is small, then we can do a good job “predicting” b by using the columns of A.

Thus far we have been assuming that $A \in \mathbb { R } ^ { m \times n }$ has full column rank, an assumption that is dropped in §5.5. However, even if rank(A) = n, trouble can be expected if A is nearly rank deficient. The SVD can be used to substantiate this remark. If

$$
A = U \Sigma V ^ {T} = \sum_ {i = 1} ^ {n} \sigma_ {i} u _ {i} v _ {i} ^ {T}
$$

is the SVD of a full rank matrix $A \in \mathbb { R } ^ { m \times n }$ , then

$$
\parallel A x - b \parallel_ {2} ^ {2} = \parallel (U ^ {T} A V) (V ^ {T} x) - U ^ {T} b \parallel_ {2} ^ {2} = \sum_ {i = 1} ^ {n} (\sigma_ {i} y _ {i} - (u _ {i} ^ {T} b)) ^ {2} + \sum_ {i = n + 1} ^ {m} (u _ {i} ^ {T} b) ^ {2}
$$

where $y = V ^ { T } x$ . It follows that this summation is minimized by setting $y _ { i } = u _ { i } ^ { T } b / \sigma _ { i }$ , $i = 1 { : } n$ . Thus,

$$
x _ {\mathrm{LS}} = \sum_ {i = 1} ^ {n} \frac {u _ {i} ^ {T} b}{\sigma_ {i}} v _ {i} \tag {5.3.2}
$$

and

$$
\rho_ {\mathrm{LS}} ^ {2} = \sum_ {i = n + 1} ^ {2} (u _ {i} ^ {T} b) ^ {2}. \tag {5.3.3}
$$

It is clear that the presence of small singular values means LS solution sensitivity. The effect of perturbations on the minimum sum of squares is less clear and requires further analysis which we offer below.

When assessing the quality of a computed LS solution $\hat { x } _ { \mathrm { L S } }$ , there are two important issues to bear in mind:

• How close is $\hat { x } _ { \mathrm { L S } }$ to $x _ { \mathrm { L S } } \mathrm { ^ { 2 } }$   
• How small is $\hat { r } _ { \mathrm { L S } } = b - A \hat { x } _ { \mathrm { L S } }$ compared to $r _ { \mathrm { L S } } = b - A x _ { \mathrm { L S } } ?$

The relative importance of these two criteria varies from application to application. In any case it is important to understand how $x _ { \mathrm { L S } }$ and $r _ { \mathrm { L S } }$ are affected by perturbations in A and b. Our intuition tells us that if the columns of A are nearly dependent, then these quantities may be quite sensitive. For example, suppose

$$
A = \left[ \begin{array}{l l} 1 & 0 \\ 0 & 1 0 ^ {- 6} \\ 0 & 0 \end{array} \right], \delta A = \left[ \begin{array}{l l} 0 & 0 \\ 0 & 0 \\ 0 & 1 0 ^ {- 8} \end{array} \right], b = \left[ \begin{array}{l} 1 \\ 0 \\ 1 \end{array} \right], \delta b = \left[ \begin{array}{l} 0 \\ 0 \\ 0 \end{array} \right],
$$

and that $x _ { \mathrm { L S } }$ and $\hat { x } _ { \mathrm { L S } }$ minimize $\parallel A x - b \parallel _ { 2 }$ and $\parallel ( A + \delta A ) x - ( b + \delta b ) \parallel _ { 2 } :$ , respectively. If $r _ { \mathrm { L S } }$ and $\hat { r } _ { \mathrm { L S } }$ are the corresponding minimum residuals, then it can be shown that

$$
x _ {\mathrm{LS}} = \left[ \begin{array}{c} 1 \\ 0 \end{array} \right], \hat {x} _ {\mathrm{LS}} = \left[ \begin{array}{c} 1 \\ . 9 9 9 9 \cdot 1 0 ^ {4} \end{array} \right], r _ {\mathrm{LS}} = \left[ \begin{array}{c} 0 \\ 0 \\ 1 \end{array} \right], \hat {r} _ {\mathrm{LS}} = \left[ \begin{array}{c} 0 \\ -. 9 9 9 9 \cdot 1 0 ^ {- 2} \\ . 9 9 9 9 \cdot 1 0 ^ {0} \end{array} \right].
$$

Recall that the 2-norm condition of a rectangular matrix is the ratio of its largest to smallest singular values. Since $\kappa _ { 2 } ( A ) = 1 0 ^ { 6 }$ we have

$$
\frac {\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \| _ {2}}{\| x _ {\mathrm{LS}} \| _ {2}} \approx . 9 9 9 9 \cdot 1 0 ^ {4} \leq \kappa_ {2} (A) ^ {2} \frac {\| \delta A \| _ {2}}{\| A \| _ {2}} = 1 0 ^ {1 2} \cdot 1 0 ^ {- 8}
$$

and

$$
\frac {\| \hat {r} _ {\mathrm{LS}} - r _ {\mathrm{LS}} \| _ {2}}{\| b \| _ {2}} \approx . 7 0 7 0 \cdot 1 0 ^ {- 2} \leq \kappa_ {2} (A) \frac {\| \delta A \| _ {2}}{\| A \| _ {2}} = 1 0 ^ {6} \cdot 1 0 ^ {- 8}.
$$

The example suggests that the sensitivity of $x _ { \mathrm { L S } }$ can depend upon $\kappa _ { 2 } ( A ) ^ { 2 }$ . Below we offer an LS perturbation theory that confirms the possibility.

# 5.3.2 The Method of Normal Equations

A widely-used method for solving the full-rank LS problem is the method of normal equations.

Algorithm 5.3.1 (Normal Equations) Given $A \in \mathbb { R } ^ { m \times n }$ with the property that rank(A) = n and $b \in \mathbb { R } ^ { m }$ , this algorithm computes a vector $x _ { \mathrm { L S } }$ that minimizes $\parallel A x - b \parallel _ { 2 }$ .

Compute the lower triangular portion of $C = A ^ { T } A$ .

Form the matrix-vector product $d = A ^ { T } b$ .

Compute the Cholesky factorization $C = G G ^ { T }$ .

Solve $G y = d$ and $G ^ { T } x _ { \mathrm { L S } } = y .$

This algorithm requires $( m + n / 3 ) n ^ { 2 }$ flops. The normal equation approach is convenient because it relies on standard algorithms: Cholesky factorization, matrix-matrix multiplication, and matrix-vector multiplication. The compression of the m-by-n data matrix A into the (typically) much smaller n-by-n cross-product matrix C is attractive.

Let us consider the accuracy of the computed normal equations solution $\hat { x } _ { \mathrm { L S } }$ . For clarity, assume that no roundoff errors occur during the formation of $C = A ^ { T } A$ and $d = A ^ { T } b$ . It follows from what we know about the roundoff properties of the Cholesky factorization (§4.2.6) that

$$
(A ^ {T} A + E) \hat {x} _ {\mathrm{LS}} = A ^ {T} b
$$

where

$$
\| E \| _ {2} \approx \mathbf {u} \| A ^ {T} \| _ {2} \| A \| _ {2} = \mathbf {u} \| A ^ {T} A \| _ {2}.
$$

Thus, we can expect

$$
\frac {\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \| _ {2}}{\| x _ {\mathrm{LS}} \| _ {2}} \approx \mathbf {u} \kappa_ {2} (A ^ {T} A) = \mathbf {u} \kappa_ {2} (A) ^ {2}. \tag {5.3.4}
$$

In other words, the accuracy of the computed normal equations solution depends on the square of the condition. See Higham (ASNA, §20.4) for a detailed roundoff analysis of the normal equations approach.

It should be noted that the formation of $A ^ { T } A$ can result in a significant loss of information. If

$$
A = \left[ \begin{array}{c c} 1 & 1 \\ \sqrt {\mathbf {u}} & 0 \\ 0 & \sqrt {\mathbf {u}} \end{array} \right],
$$

then $\kappa _ { 2 } ( A ) \approx \sqrt { \bf { u } }$ . However,

$$
\mathsf {f l} (A ^ {T} A) = \left[ \begin{array}{l l} 1 & 1 \\ 1 & 1 \end{array} \right]
$$

is exactly singular. Thus, the method of normal equations can break down on matrices that are not particularly close to being numerically rank deficient.

# 5.3.3 LS Solution Via QR Factorization

Let $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ and $b \in \mathbb { R } ^ { m }$ be given and suppose that an orthogonal matrix $Q \in \mathbb { R } ^ { m \times m }$ has been computed such that

$$
Q ^ {T} A = R = \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right] _ {m - n} ^ {n} \tag {5.3.5}
$$

is upper triangular. If

$$
Q ^ {T} b = \left[ \begin{array}{l} c \\ d \end{array} \right] _ {m - n} ^ {n}
$$

then

$$
\parallel A x - b \parallel_ {2} ^ {2} = \parallel Q ^ {T} A x - Q ^ {T} b \parallel_ {2} ^ {2} = \parallel R _ {1} x - c \parallel_ {2} ^ {2} + \parallel d \parallel_ {2} ^ {2}
$$

for any $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . Since $\mathsf { r a n k } ( A ) = \mathsf { r a n k } ( R _ { 1 } ) = n$ , it follows that $x _ { \mathrm { L S } }$ is defined by the upper triangular system

$$
R _ {1} x _ {\mathrm{LS}} = c.
$$

Note that

$$
\rho_ {\mathrm{LS}} = \parallel d \parallel_ {2}.
$$

We conclude that the full-rank LS problem can be readily solved once we have computed the QR factorization of A. Details depend on the exact QR procedure. If Householder matrices are used and $Q ^ { T }$ is applied in factored form to b, then we obtain

Algorithm 5.3.2 (Householder LS Solution) If $A \in \mathbb { R } ^ { m \times n }$ has full column rank and $b \in \mathbb { R } ^ { m }$ , then the following algorithm computes a vector $\boldsymbol { x } _ { \mathrm { L S } } \in \mathbb { R } ^ { n }$ such that $\parallel A x _ { \mathrm { L S } } - b \parallel _ { 2 }$ is minimum.

Use Algorithm 5.2.1 to overwrite A with its QR factorization.

for $j = 1 { : } n$

$$
v = \left[ \begin{array}{c} 1 \\ A (j + 1: m, j) \end{array} \right]
$$

$$
\beta = 2 / v ^ {T} v
$$

$$
b (j: m) = b (j: m) - \beta (v ^ {T} b (j: m)) v
$$

end

Solve $R ( 1 : n , 1 : n ) \cdot x _ { \mathrm { { L S } } } = b ( 1 { : } n )$ .

This method for solving the full-rank LS problem requires $2 n ^ { 2 } ( m - n / 3 )$ flops. The O(mn) flops associated with the updating of b and the $O ( n ^ { 2 } )$ flops associated with the back substitution are not significant compared to the work required to factor A.

It can be shown that the computed $\hat { x } _ { \mathrm { L S } }$ solves

$$
\min \| (A + \delta A) x - (b + \delta b) \| _ {2} \tag {5.3.6}
$$

where

$$
\left\| \delta A \right\| _ {F} \leq (6 m - 3 n + 4 1) n \mathbf {u} \left\| A \right\| _ {F} + O \left(\mathbf {u} ^ {2}\right) \tag {5.3.7}
$$

and

$$
\| \delta b \| _ {2} \leq (6 m - 3 n + 4 0) n \mathbf {u} \| b \| _ {2} + O (\mathbf {u} ^ {2}). \tag {5.3.8}
$$

These inequalities are established in Lawson and Hanson (SLS, p. 90ff) and show that $\hat { x } _ { \mathrm { L S } }$ satisfies a “nearby” LS problem. (We cannot address the relative error in $\hat { x } _ { \mathrm { L S } }$ without an LS perturbation theory, to be discussed shortly.) We mention that similar results hold if Givens QR is used.

# 5.3.4 Breakdown in Near-Rank-Deficient Case

As with the method of normal equations, the Householder method for solving the LS problem breaks down in the back-substitution phase if rank $( A ) ~ < ~ n$ . Numerically, trouble can be expected if $\kappa _ { 2 } ( A ) = \kappa _ { 2 } ( R ) \approx 1 / \mathbf { u }$ . This is in contrast to the normal equations approach, where completion of the Cholesky factorization becomes problematical once $\kappa _ { 2 } ( A )$ is in the neighborhood of $1 / \sqrt { \mathbf { u } }$ as we showed above. Hence the claim in Lawson and Hanson (SLS, pp. 126–127) that for a fixed machine precision, a wider class of LS problems can be solved using Householder orthogonalization.

# 5.3.5 A Note on the MGS Approach

In principle, MGS computes the thin QR factorization $A = Q _ { 1 } R _ { 1 }$ . This is enough to solve the full-rank LS problem because it transforms the normal equation system $( A ^ { T } A ) x = A ^ { T } b$ to the upper triangular system $R _ { 1 } x = Q _ { 1 } ^ { T } b$ . But an analysis of this approach when $Q _ { 1 } ^ { T } b$ is explicitly formed introduces a $\kappa _ { 2 } ( A ) ^ { 2 }$ term. This is because the computed factor $\hat { Q } _ { 1 }$ satisfies $\Vert \hat { Q } _ { 1 } ^ { T } \hat { Q } _ { 1 } - I _ { n } \Vert _ { 2 } \approx \mathbf { u } \kappa _ { 2 } ( A )$ as we mentioned in §5.2.9.

However, if MGS is applied to the augmented matrix

$$
A _ {+} = \left[ \begin{array}{c c} A & b \end{array} \right] = \left[ \begin{array}{c c} Q _ {1} & q _ {n + 1} \end{array} \right] \left[ \begin{array}{c c} R _ {1} & z \\ 0 & \rho \end{array} \right],
$$

then $z = Q _ { 1 } ^ { T } b$ . Computing $Q _ { 1 } ^ { T } b$ in this fashion and solving $R _ { 1 } x _ { \mathrm { L S } } = z$ produces an LS solution $\hat { x } _ { \mathrm { L S } }$ that is “just as good” as the Householder QR method. That is to say, a result of the form (5.3.6)–(5.3.8) applies. See Bj¨orck and Paige (1992).

It should be noted that the MGS method is slightly more expensive than Householder QR because it always manipulates m-vectors whereas the latter procedure deals with vectors that become shorter in length as the algorithm progresses.

# 5.3.6 The Sensitivity of the LS Problem

We now develop a perturbation theory for the full-rank LS problem that assists in the comparison of the normal equations and QR approaches. LS sensitivity analysis has a long and fascinating history. Grcar (2009, 2010) compares about a dozen different results that have appeared in the literature over the decades and the theorem below follows his analysis. It examines how the LS solution and its residual are affected by changes in A and b and thereby sheds light on the condition of the LS problem. Four facts about $A \in \mathbb { R } ^ { m \times n }$ are used in the proof, where it is assumed that $m > n \colon$

$$
1 = \| A (A ^ {T} A) ^ {- 1} A ^ {T} \| _ {2}, \quad \frac {1}{\sigma_ {n} (A)} = \| (A ^ {T} A) ^ {- 1} A ^ {T} \| _ {2}, \tag {5.3.9}
$$

$$
1 = \| I - A (A ^ {T} A) ^ {- 1} A ^ {T} \| _ {2}, \quad \frac {1}{\sigma_ {n} (A) ^ {2}} = \| (A ^ {T} A) ^ {- 1} \| _ {2}.
$$

These equations are easily verified using the SVD.

Theorem 5.3.1. Suppose that $x _ { \mathrm { L S } } , \ r _ { \mathrm { L S } } , \ \hat { x } _ { \mathrm { L S } }$ , and $\hat { r } _ { \mathrm { L S } }$ satisfy

$$
\left\| A x _ {\mathrm{LS}} - b \right\| _ {2} = \min, \quad r _ {\mathrm{LS}} = b - A x _ {\mathrm{LS}},
$$

$$
\| (A + \delta A) \hat {x} _ {\mathrm{LS}} - (b + \delta b) \| _ {2} = \min, \quad \hat {r} _ {\mathrm{LS}} = (b + \delta b) - (A + \delta A) \hat {x} _ {\mathrm{LS}},
$$

where A has rank n and $\parallel \delta A \parallel _ { 2 } < \sigma _ { n } ( A )$ . Assume that b, $r _ { \mathrm { L S } }$ , and $x _ { \mathrm { L S } }$ are not zero. Let $\theta _ { \mathrm { L S } } \in ( 0 , \pi / 2 )$ be defined by

$$
\sin (\theta_ {\mathrm{LS}}) = \frac {\parallel r _ {\mathrm{LS}} \parallel_ {2}}{\parallel b \parallel_ {2}}.
$$

If

$$
\epsilon = \max \left\{\frac {\parallel \delta A \parallel_ {2}}{\parallel A \parallel_ {2}}, \frac {\parallel \delta b \parallel_ {2}}{\parallel b \parallel_ {2}} \right\}
$$

and

$$
\nu_ {\mathrm{LS}} = \frac {\| A x _ {\mathrm{LS}} \| _ {2}}{\sigma_ {n} (A) \| x _ {\mathrm{LS}} \| _ {2}}, \tag {5.3.10}
$$

then

$$
\frac {\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \| _ {2}}{\| x \| _ {2}} \leq \epsilon \left\{\frac {v _ {\mathrm{LS}}}{\cos (\theta_ {\mathrm{LS}})} + [ 1 + \nu_ {\mathrm{LS}} \tan (\theta_ {\mathrm{LS}}) ] \kappa_ {2} (A) \right\} + O (\epsilon^ {2}) \tag {5.3.11}
$$

and

$$
\frac {\| \hat {r} _ {\mathrm{LS}} - r _ {\mathrm{LS}} \| _ {2}}{\| r _ {\mathrm{LS}} \| _ {2}} \leq \epsilon \left\{\frac {1}{\sin (\theta_ {\mathrm{LS}})} + \left[ \frac {1}{\nu_ {\mathrm{LS}} \tan (\theta_ {\mathrm{LS}})} + 1 \right] \kappa_ {2} (A) \right\} + O (\epsilon^ {2}). \tag {5.3.12}
$$

Proof. Let E and f be defined by $E = \delta A / \epsilon$ and $f = \delta b / \epsilon$ . By Theorem 2.5.2 we have rank $( A + t E ) = n$ for all $t \in [ 0 , \epsilon ]$ . It follows that the solution $x ( t )$ t o

$$
(A + t E) ^ {T} (A + t E) x (t) = (A + t E) ^ {T} (b + t f) \tag {5.3.13}
$$

is continuously differentiable for all $t \in [ 0 , \epsilon ]$ . Since $x _ { \mathrm { L S } } = x ( 0 )$ and $\hat { x } _ { \mathrm { L S } } = x ( \epsilon )$ , we have

$$
\hat {x} _ {\mathrm{LS}} = x _ {\mathrm{LS}} + \epsilon \dot {x} (0) + O (\epsilon^ {2}).
$$

By taking norms and dividing by $\Vert \boldsymbol { x } _ { \mathrm { L S } } \Vert _ { 2 }$ we obtain

$$
\frac {\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \| _ {2}}{\| x _ {\mathrm{LS}} \| _ {2}} = \epsilon \frac {\| \dot {x} (0) \| _ {2}}{\| x _ {\mathrm{LS}} \| _ {2}} + O (\epsilon^ {2}). \tag {5.3.14}
$$

In order to bound $\Vert \dot { x } ( 0 ) \Vert _ { 2 }$ , we differentiate (5.3.13) and set t = 0 in the result. This gives

$$
E ^ {T} A x _ {\mathrm{LS}} + A ^ {T} E x _ {\mathrm{LS}} + A ^ {T} A \dot {x} (0) = A ^ {T} f + E ^ {T} b,
$$

i.e.,

$$
\dot {x} (0) = (A ^ {T} A) ^ {- 1} A ^ {T} (f - E x _ {\mathrm{LS}}) + (A ^ {T} A) ^ {- 1} E ^ {T} r _ {\mathrm{LS}}. \tag {5.3.15}
$$

Using (5.3.9) and the inequalities $\left\| \ f \right\| _ { 2 } \leq \left\| \ b \right\| _ { 2 }$ and $\left. \ E \right. _ { 2 } \leq \left. \ A \right. _ { 2 }$ , it follows that

$$
\| \dot {x} (0) \| \leq \| (A ^ {T} A) ^ {- 1} A ^ {T} f \| _ {2} + \| (A ^ {T} A) ^ {- 1} A ^ {T} E x _ {\mathrm{LS}} \| _ {2} + \| (A ^ {T} A) ^ {- 1} E ^ {T} r _ {\mathrm{LS}} \| _ {2}
$$

$$
\leq \frac {\| b \| _ {2}}{\sigma_ {n} (A)} + \frac {\| A \| _ {2} \| x _ {\mathrm{LS}} \| _ {2}}{\sigma_ {n} (A)} + \frac {\| A \| _ {2} \| r _ {\mathrm{LS}} \| _ {2}}{\sigma_ {n} (A) ^ {2}}.
$$

By substituting this into (5.3.14) we obtain

$$
\frac {\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \| _ {2}}{\| x _ {\mathrm{LS}} \| _ {2}} \leq \epsilon \left(\frac {\| b \| _ {2}}{\sigma_ {n} (A) \| x _ {\mathrm{LS}} \| _ {2}} + \frac {\| A \| _ {2}}{\sigma_ {n} (A)} + \frac {\| A \| _ {2} \| r _ {\mathrm{LS}} \| _ {2}}{\sigma_ {n} (A) ^ {2} \| x _ {\mathrm{LS}} \| _ {2}}\right) + O (\epsilon^ {2}).
$$

Inequality (5.3.11) follows from the definitions of $\kappa _ { 2 } ( A )$ and $\nu _ { \mathrm { L S } }$ and the identities

$$
\cos (\theta_ {\mathrm{LS}}) = \frac {\parallel A x _ {\mathrm{LS}} \parallel_ {2}}{\parallel b \parallel_ {2}}, \quad \tan (\theta_ {\mathrm{LS}}) = \frac {\parallel r _ {\mathrm{LS}} \parallel_ {2}}{\parallel A x _ {\mathrm{LS}} \parallel_ {2}}. \tag {5.3.16}
$$

The proof of the residual bound (5.3.12) is similar. Define the differentiable vector function r(t) by

$$
r (t) = (b + t f) - (A + t E) x (t)
$$

and observe that $r _ { \mathrm { L S } } = r ( 0 )$ and $\hat { r } _ { \mathrm { L S } } = r ( \epsilon )$ . Thus,

$$
\frac {\left\| \hat {r} _ {\mathrm{LS}} - r _ {\mathrm{LS}} \right\| _ {2}}{\left\| r _ {\mathrm{LS}} \right\| _ {2}} = \epsilon \frac {\left\| \dot {r} (0) \right\| _ {2}}{\left\| r _ {\mathrm{LS}} \right\| _ {2}} + O \left(\epsilon^ {2}\right). \tag {5.3.17}
$$

From (5.3.15) we have

$$
\dot {r} (0) = \left(I - A (A ^ {T} A) ^ {- 1} A ^ {T}\right) (f - E x _ {\mathrm{LS}}) - A (A ^ {T} A) ^ {- 1} E ^ {T} r _ {\mathrm{LS}}.
$$

By taking norms, using (5.3.9) and the inequalities  $f \parallel _ { 2 } \leq \parallel b \parallel _ { 2 }$ and $\| E \| _ { 2 } \leq \| A \| _ { 2 }$ we obtain

$$
\| \dot {r} (0) \| _ {2} \leq \| b \| _ {2} + \| A \| _ {2} \| x _ {\mathrm{LS}} \| _ {2} + \frac {\| A \| _ {2} \| r _ {\mathrm{LS}} \| _ {2}}{\sigma_ {n} (A)}
$$

and thus from (5.3.17) we have

$$
\frac {\| \hat {r} _ {\mathrm{LS}} - r _ {\mathrm{LS}} \| _ {2}}{\| r _ {\mathrm{LS}} \| _ {2}} \leq \frac {\| b \| _ {2}}{\| r _ {\mathrm{LS}} \| _ {2}} + \frac {\| A \| _ {2} \| x _ {\mathrm{LS}} \| _ {2}}{\| r _ {\mathrm{LS}} \| _ {2}} + \frac {\| A \| _ {2}}{\sigma_ {n} (A)}.
$$

The inequality (5.3.12) follows from the definitions of $\kappa _ { 2 } ( A )$ and $\nu _ { \mathrm { L S } }$ and the identities (5.3.16).

It is instructive to identify conditions that turn the upper bound in (5.3.11) into a bound that involves $\kappa _ { 2 } ( A ) ^ { 2 }$ . The example in §5.3.1 suggests that this factor might figure in the definition of an LS condition number. However, the theorem shows that the situation is more subtle. Note that

$$
\nu_ {\mathrm{LS}} = \frac {\parallel A x _ {\mathrm{LS}} \parallel_ {2}}{\sigma_ {n} (A) \parallel x _ {\mathrm{LS}} \parallel_ {2}} \leq \frac {\parallel A \parallel_ {2}}{\sigma_ {n} (A)} = \kappa_ {2} (A).
$$

The SVD expansion (5.3.2) suggests that if b has a modest component in the direction of the left singular vector $u _ { n }$ , then

$$
\nu_ {\mathrm{LS}} \approx \kappa_ {2} (A).
$$

If this is the case and $\theta _ { \mathrm { L S } }$ is sufficiently bounded away from $\pi / 2$ , then the inequality (5.3.11) essentially says that

$$
\frac {\left\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \right\| _ {2}}{\left\| x _ {\mathrm{LS}} \right\| _ {2}} \approx \epsilon \left(\kappa_ {2} (A) + \frac {\rho_ {\mathrm{LS}}}{\left\| b \right\| _ {2}} \kappa_ {2} (A) ^ {2}\right). \tag {5.3.18}
$$

Although this simple heuristic assessment of LS sensitivity is almost always applicable, it important to remember that the true condition of a particular LS problem depends on $\nu _ { \mathrm { L S } } , \theta _ { \mathrm { L S } }$ , and $\kappa _ { 2 } ( A )$ .

Regarding the perturbation of the residual, observe that the upper bound in the residual result (5.3.12) is less than the upper bound in the solution result (5.3.11) by a factor of $\nu _ { \mathrm { L S } } \tan ( \theta _ { \mathrm { L S } } )$ . We also observe that if $\theta _ { \mathrm { L S } }$ is sufficiently bounded away from both 0 and $\pi / 2$ , then (5.3.12) essentially says that

$$
\frac {\left\| \hat {r} _ {\mathrm{LS}} - r _ {\mathrm{LS}} \right\| _ {2}}{\left\| r _ {\mathrm{LS}} \right\| _ {2}} \approx \epsilon \cdot \kappa_ {2} (A). \tag {5.3.19}
$$

For more insights into the subtleties behind Theorem 5.3.1., see Wedin (1973), Vandersluis (1975), Bj¨orck (NMLS, p. 30), Higham (ASNA, p. 382), and Grcar(2010).

# 5.3.7 Normal Equations Versus QR

It is instructive to compare the normal equation and QR approaches to the full-rank LS problem in light of Theorem 5.3.1.

• The method of normal equations produces an $\hat { x } _ { \mathrm { L S } }$ whose relative error depends on $\kappa _ { 2 } ( A ) ^ { 2 }$ , a factor that can be considerably larger than the condition number associated with a “small residual” LS problem.   
• The QR approach (Householder, Givens, careful MGS) solves a nearby LS problem. Therefore, these methods produce a computed solution with relative error that is “predicted” by the condition of the underlying LS problem.

Thus, the QR approach is more appealing in situations where b is close to the span of A’s columns.

Finally, we mention two other factors that figure in the debate about QR versus normal equations. First, the normal equations approach involves about half of the arithmetic when m 
 n and does not require as much storage, assuming that $Q ( : , 1 : n )$ is required. Second, QR approaches are applicable to a wider class of LS problems. This is because the Cholesky solve in the method of normal equations is “in trouble” if $\kappa _ { 2 } ( A ) \approx 1 / \sqrt { \mathbf { u } }$ while the R-solve step in a QR approach is in trouble only if $\kappa _ { 2 } ( A ) \approx$ 1/u. Choosing the “right” algorithm requires having an appreciation for these tradeoffs.

# 5.3.8 Iterative Improvement

A technique for refining an approximate LS solution has been analyzed by Bj¨orck (1967, 1968). It is based on the idea that if

$$
\left[ \begin{array}{l l} I _ {m} & A \\ A ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} r \\ x \end{array} \right] = \left[ \begin{array}{l} b \\ 0 \end{array} \right], \quad A \in \mathbb {R} ^ {m \times n}, b \in \mathbb {R} ^ {m}, \tag {5.3.20}
$$

then $\parallel b - A x \parallel _ { 2 } = \operatorname* { m i n }$ . This follows because $r + A x = b$ and $A ^ { T } r = 0$ imply $A ^ { T } A x =$ $A ^ { T } b$ . The above augmented system is nonsingular if $\mathsf { r a n k } ( A ) = n$ , which we hereafter assume. By casting the LS problem in the form of a square linear system, the iterative improvement scheme §3.5.3 can be applied:

$$
r ^ {(0)} = 0, x ^ {(0)} = 0
$$

$\mathbf { f o r } \ k = 0 , 1 , \ldots$

$$
\begin{array}{l} \left[ \begin{array}{c} f ^ {(k)} \\ g ^ {(k)} \end{array} \right] = \left[ \begin{array}{c} b \\ 0 \end{array} \right] - \left[ \begin{array}{c c} I & A \\ A ^ {T} & 0 \end{array} \right] \left[ \begin{array}{c} r ^ {(k)} \\ x ^ {(k)} \end{array} \right] \\ \left[ \begin{array}{c c} I & A \\ A ^ {T} & 0 \end{array} \right] \left[ \begin{array}{c} p ^ {(k)} \\ z ^ {(k)} \end{array} \right] = \left[ \begin{array}{c} f ^ {(k)} \\ g ^ {(k)} \end{array} \right] \\ \left[ \begin{array}{l} r ^ {(k + 1)} \\ x ^ {(k + 1)} \end{array} \right] = \left[ \begin{array}{l} r ^ {(k)} \\ x ^ {(k)} \end{array} \right] + \left[ \begin{array}{l} p ^ {(k)} \\ z ^ {(k)} \end{array} \right] \\ \end{array}
$$

end

The residuals $f ^ { ( k ) }$ and $g ^ { ( k ) }$ must be computed in higher precision, and an original copy of A must be around for this purpose.

If the QR factorization of A is available, then the solution of the augmented system is readily obtained. In particular, if $A = Q R$ and $R _ { 1 } = R ( 1 { : } n , 1 { : } n )$ , then a system of the form

$$
{\left[ \begin{array}{c c} I & A \\ A ^ {T} & 0 \end{array} \right]} {\left[ \begin{array}{c} p \\ z \end{array} \right]} = {\left[ \begin{array}{c} f \\ g \end{array} \right]}
$$

transforms to

$$
\left[ \begin{array}{c c c} I _ {n} & 0 & R _ {1} \\ 0 & I _ {m - n} & 0 \\ R _ {1} ^ {T} & 0 & 0 \end{array} \right] \left[ \begin{array}{c} h \\ f _ {2} \\ z \end{array} \right] = \left[ \begin{array}{c} f _ {1} \\ f _ {2} \\ g \end{array} \right]
$$

where

$$
Q ^ {T} f = \left[ \begin{array}{l} f _ {1} \\ f _ {2} \end{array} \right] _ {m - n} ^ {n}, \qquad Q ^ {T} p = \left[ \begin{array}{l} h \\ f _ {2} \end{array} \right] _ {m - n} ^ {n}.
$$

Thus, p and z can be determined by solving the triangular systems $R _ { 1 } ^ { T } h = g$ and $R _ { 1 } z = f _ { 1 } - h$ and setting

$$
p = Q \left[ \begin{array}{c} h \\ f _ {2} \end{array} \right].
$$

Assuming that Q is stored in factored form, each iteration requires 8mn $- 2 n ^ { 2 }$ flops.

The key to the iteration’s success is that both the LS residual and solution are updated—not just the solution. Bj¨orck (1968) shows that if $\kappa _ { 2 } ( A ) \approx \beta ^ { q }$ and t-digit, β-base arithmetic is used, then $x ^ { ( k ) }$ has approximately $k ( t - q )$ correct base-β digits, provided the residuals are computed in double precision. Notice that it is $\kappa _ { 2 } ( A )$ , not $\kappa _ { 2 } { \left( A \right) } ^ { 2 }$ , that appears in this heuristic.

# 5.3.9 Some Point/Line/Plane Nearness Problems in 3-Space

The fields of computer graphics and computer vision are replete with many interesting matrix problems. Below we pose three geometric “nearness” problems that involve points, lines, and planes in 3-space. Each is a highly structured least squares problem with a simple, closed-form solution. The underlying trigonometry leads rather naturally to the vector cross product, so we start with a quick review of this important operation.

The cross product of a vector $\boldsymbol { p } \in \mathbb { R } ^ { 3 }$ with a vector $q \in \mathbb { R } ^ { 3 }$ is defined by

$$
p \times q = \left[ \begin{array}{l} p _ {2} q _ {3} - p _ {3} q _ {2} \\ p _ {3} q _ {1} - p _ {1} q _ {3} \\ p _ {1} q _ {2} - p _ {2} q _ {1} \end{array} \right].
$$

This operation can be framed as a matrix-vector product. For any $v \in \mathbb { R } ^ { 3 }$ , define the skew-symmetric matrix $v ^ { c }$ by

$$
v ^ {c} = \left[ \begin{array}{c c c} 0 & - v _ {3} & v _ {2} \\ v _ {3} & 0 & - v _ {1} \\ - v _ {2} & v _ {1} & 0 \end{array} \right].
$$

It follows that

$$
p \times q = p ^ {c} \cdot q = - q ^ {c} \cdot p = - (q \times p).
$$

Using the skew-symmetry of $p ^ { c }$ and $q ^ { c }$ , it is easy to show that

$$
p \times q \in \operatorname{span} \{p, q \} ^ {\perp}. \tag {5.3.21}
$$

Other properties include

$$
(p \times q) \times r = (p ^ {c} \cdot q) ^ {c} r = (q p ^ {T} - p q ^ {T}) r = (p ^ {T} r) \cdot q - (q ^ {T} r) \cdot p, \tag {5.3.22}
$$

$$
(p \times q) ^ {T} (r \times s) = (p ^ {c} q) ^ {T} \cdot (r ^ {c} s) = \det ([ p q ] ^ {T} [ r s ]), \tag {5.3.23}
$$

$$
p ^ {c} p ^ {c} = p p ^ {T} - \parallel p \parallel_ {2} ^ {2} \cdot I _ {3}, \tag {5.3.24}
$$

$$
\parallel p ^ {c} q \parallel_ {2} ^ {2} = \parallel p \parallel_ {2} ^ {2} \cdot \parallel q \parallel_ {2} ^ {2} \cdot \left(1 - \left(\frac {p ^ {T} q}{\parallel p \parallel_ {2} \cdot \parallel q \parallel_ {2}}\right) ^ {2}\right). \tag {5.3.25}
$$

We are now set to state the three problems and specify their theoretical solutions. For hints at how to establish the correctness of the solutions, see P5.3.13–P5.3.15.

Problem 1. Given a line L and a point y, find the point $z ^ { \mathrm { o p t } }$ on $L$ that is closest to $y .$ i.e., solve

$$
\min _ {z \in L} \| z - y \| _ {2}.
$$

If L passes through distinct points $p _ { 1 }$ and $p _ { 2 }$ , then it can be shown that

$$
z ^ {\mathrm{opt}} = y + \frac {1}{v ^ {T} v} v ^ {c} v ^ {c} (y - p _ {1}), \quad v = p _ {2} - p _ {1}. \tag {5.3.26}
$$

Problem 2. Given lines $L _ { 1 }$ and $L _ { 2 }$ , find the point $z _ { 1 } ^ { \mathrm { o p t } }$ on $L _ { 1 }$ that is closest to $L _ { 2 }$ and the point $z _ { 2 } ^ { \mathrm { o p t } }$ on $L _ { 2 }$ that is closest to $L _ { \mathrm { 1 } } , \mathrm { i . e . }$ , solve

$$
\min _ {z _ {1} \in L _ {1}, z _ {2} \in L _ {2}} \| z _ {1} - z _ {2} \| _ {2}.
$$

If $L _ { 1 }$ passes through distinct points $p _ { 1 }$ and $p _ { 2 }$ and $L _ { 2 }$ passes through distinct points $q _ { 1 }$ and $q _ { 2 }$ , then it can be shown that

$$
z _ {1} ^ {\mathrm{opt}} = p _ {1} + \frac {1}{r ^ {T} r} \cdot v w ^ {T} \cdot r ^ {c} (q _ {1} - p _ {1}), \tag {5.3.27}
$$

$$
z _ {2} ^ {\mathrm{opt}} = q _ {1} + \frac {1}{r ^ {T} r} \cdot w v ^ {T} \cdot r ^ {c} (q _ {1} - p _ {1}), \tag {5.3.28}
$$

where v = p2 − p1, w = q2 − q1, and $r = v ^ { c } w$ .

Problem 3. Given a plane $P$ and a point y, find the point $z ^ { \mathrm { o p t } }$ on $P$ that is closest to $y , { \mathrm { i . e . } }$ , solve

$$
\min _ {z \in P} \| z - y \| _ {2}.
$$

If P passes through three distinct points $p _ { 1 } , p _ { 2 }$ , and $p _ { 3 }$ , then it can be shown that

$$
z ^ {\mathrm{opt}} = p _ {1} - \frac {1}{v ^ {T} v} \cdot v ^ {c} v ^ {c} (y - p _ {1}) \tag {5.3.29}
$$

where $v = ( p _ { 2 } - p _ { 1 } ) ^ { c } ( p _ { 3 } - p _ { 1 } )$ .

The nice closed-form solutions (5.3.26)–(5.3.29) are deceptively simple and great care must be exercised when computing with these formulae or their mathematical equivalents. See Kahan (2011).

# Problems

P5.3.1 Assume $A ^ { T } A x = A ^ { T } b , ( A ^ { T } A + F ) \hat { x } = A ^ { T } b ,$ , and $2 \| F \| _ { 2 } \leq \sigma _ { n } ( A ) ^ { 2 }$ . Show that if $r = b - A x$ and $\hat { r } = b - A \hat { x }$ , then $\hat { r } - r = A ( A ^ { T } A + F ) ^ { - 1 } F x$ and

$$
\| \hat {r} - r \| _ {2} \leq 2 \kappa_ {2} (A) \frac {\| F \| _ {2}}{\| A \| _ {2}} \| x \| _ {2}.
$$

P5.3.2 Assume that $A ^ { T } A x = A ^ { T } b$ and that $A ^ { T } A \hat { x } = A ^ { T } b + f$ where $\parallel f \parallel _ { 2 } \leq c \mathbf { u } \parallel A ^ { T } \parallel _ { 2 } \parallel b \parallel _ { 2 }$ and A has full column rank. Show that

$$
\frac {\| x - \hat {x} \| _ {2}}{\| x \| _ {2}} \leq c \mathbf {u} \kappa_ {2} (A) ^ {2} \frac {\| A ^ {T} \| _ {2} \| b \| _ {2}}{\| A ^ {T} b \|}.
$$

P5.3.3 Let $A \in \mathbb { R } ^ { m \times n } ( m \geq n ) , w \in \mathbb { R } ^ { n }$ , and define

$$
B = \left[ \begin{array}{c} A \\ w ^ {T} \end{array} \right].
$$

Show that $\sigma _ { n } ( B ) \geq \sigma _ { n } ( A )$ and $\sigma _ { 1 } ( B ) \leq \sqrt { \| A \| _ { 2 } ^ { 2 } + \| w \| _ { 2 } ^ { 2 } }$ . Thus, the condition of a matrix may increase or decrease if a row is added.

P5.3.4 (Cline 1973) Suppose that $A \in \mathbb { R } ^ { m \times n }$ has rank n and that Gaussian elimination with partial pivoting is used to compute the factorization $P A = L U$ , where $L \in \mathbb { R } ^ { m \times n }$ is unit lower triangular, $\bar { U } \in \mathbb { R } ^ { n \times n }$ is upper triangular, and $P \in \mathbb { R } ^ { m \times m }$ is a permutation. Explain how the decomposition in P5.2.5 can be used to find a vector $z \in \mathbb { R } ^ { n }$ such that  $L z - P b \| _ { 2 }$ is minimized. Show that if $U x = z ,$ then $\parallel A x - b \parallel _ { 2 }$ is minimum. Show that this method of solving the LS problem is more efficient than Householder QR from the flop point of view whenever m $\leq 5 n / 3$ .

P5.3.5 The matrix $C = ( A ^ { T } A ) ^ { - 1 }$ , where rank $( A ) = n$ , arises in many statistical applications. Assume that the factorization $A = Q R$ is available. (a) Show $C = ( R ^ { T } \dot { R } ) ^ { - 1 }$ . (b) Give an algorithm for computing the diagonal of C that requires $n ^ { 3 } / 3$ flops. (c) Show that

$$
R = \left[ \begin{array}{c c} \alpha & v ^ {T} \\ 0 & S \end{array} \right] \qquad \Rightarrow \qquad C = (R ^ {T} R) ^ {- 1} = \left[ \begin{array}{c c} (1 + v ^ {T} C _ {1} v) / \alpha^ {2} & - v ^ {T} C _ {1} / \alpha \\ - C _ {1} v / \alpha & C _ {1} \end{array} \right]
$$

where $C _ { 1 } = ( S ^ { T } S ) ^ { - 1 }$ . (d) Using (c), give an algorithm that overwrites the upper triangular portion of R with the upper triangular portion of C. Your algorithm should require $2 { n } ^ { 3 } / 3$ flops.

P5.3.6 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and that $r = b - A \acute { \iota }$ x where $r , b , x \in \mathbb { R } ^ { n }$ and x is nonzero. Show how to compute a symmetric $E \in \mathbb { R } ^ { n \times n }$ with minimal Frobenius norm so that $( A + E ) x = b$ . Hint: Use the QR factorization of $\left[ \boldsymbol { x } \mid \boldsymbol { r } \right]$ and note that $E x = r \Rightarrow ( Q ^ { T } E Q ) ( Q ^ { T } x ) = \dot { Q } ^ { T } r$ .

P5.3.7 Points $P _ { 1 } , \ldots , P _ { n }$ on the x-axis have x-coordinates $x _ { 1 } , \ldots , x _ { n }$ . We know that $x _ { 1 } = 0$ and wish to compute $x _ { 2 } , \ldots , x _ { n }$ given that we have estimates $d _ { i j }$ of the separations:

$$
x _ {i} - x _ {j} \approx d _ {i j}, \qquad 1 \leq i <   j \leq n.
$$

Using the method of normal equations, show how to minimize

$$
\phi (x _ {1}, \ldots , x _ {n}) = \sum_ {i = 1} ^ {n - 1} \sum_ {j = i + 1} ^ {n} (x _ {i} - x _ {j} - d _ {i j}) ^ {2}
$$

subject to the constraint $x _ { 1 } = 0$ .

P5.3.8 Suppose $A \in \mathbb { R } ^ { m \times n }$ has full rank and that $b \in \mathbb { R } ^ { m }$ and $c \in \mathbb { R } ^ { n }$ are given. Show how to compute $\alpha = c ^ { T } x _ { \mathrm { L S } }$ without computing $x _ { \mathrm { L S } }$ explicitly. Hint: Suppose $Z$ is a Householder matrix such that $Z ^ { T } c$ is a multiple of $I _ { n } ( : , n )$ . It follows that $\alpha = ( Z ^ { T } c ) ^ { \bar { T } } y _ { \mathrm { L S } }$ where $y _ { \mathrm { L S } }$ minimizes $\parallel \tilde { A } y - b \parallel _ { 2 }$ with $y = Z ^ { T } x$ and $\tilde { A } = A Z$ .

P5.3.9 Suppose $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ with $m \geq n$ . How would you solve the full rank least squares problem given the availability of a matrix $\boldsymbol { M } \in \mathbb { R } ^ { m \times m }$ such that $M ^ { T } A = S$ is upper triangular and $\mathbf { \overset { \cdot } { M } } ^ { T } \mathbf { \overset { \cdot } { M } } = \mathbf { \overset { \cdot } { D } }$ is diagonal?

P5.3.10 Let $A \in \mathbb { R } ^ { m \times n }$ have rank n and for $\alpha \geq 0$ define

$$
M (\alpha) = \left[ \begin{array}{c c} \alpha I _ {m} & A \\ A ^ {T} & 0 \end{array} \right].
$$

Show that

$$
\sigma_ {m + n} (M (\alpha)) = \min \left\{\alpha , - \frac {\alpha}{2} + \sqrt {\sigma_ {n} (A) ^ {2} + \left(\frac {\alpha}{2}\right) ^ {2}} \right\}
$$

and determine the value of α that minimizes $\kappa _ { 2 } ( M ( \alpha ) )$ .

P5.3.11 Another iterative improvement method for LS problems is the following:

$$
x ^ {(0)} = 0
$$

$\mathbf { f o r } \ k = 0 , 1 , \ldots$

$$
r ^ {(k)} = b - A x ^ {(k)} \quad (\text { double   precision })
$$

$$
\| A z ^ {(k)} - r ^ {(k)} \| _ {2} = \min
$$

$$
x ^ {(k + 1)} = x ^ {(k)} + z ^ {(k)}
$$

end

(a) Assuming that the QR factorization of A is available, how many flops per iteration are required? (b) Show that the above iteration results by setting $g ^ { ( k ) } = 0$ in the iterative improvement scheme given in 5.3.8.

P5.3.12 Verify (5.3.21)–(5.3.25).

P5.3.13 Verify (5.3.26) noting that $L = \left\{ p _ { 1 } + \tau ( p _ { 2 } - p _ { 1 } ) : \tau \in \mathbb { R } \right\}$

P5.3.14 Verify (5.3.27) noting that the minimizer $\tau ^ { \mathrm { o p t } } \in \mathbb { R } ^ { 2 } \ \mathrm { o f } \ \| \ ( p _ { 1 } - q _ { 1 } ) \ - \ \big [ \ : p _ { 2 } - p _ { 1 } \ : | \ : q _ { 2 } - q _ { 1 } \ : \big ] \tau \| _ { 2 }$ is relevant.

P5.3.15 Verify (5.3.29) noting that $P = \{ x : x ^ { T } ( ( p _ { 2 } - p _ { 1 } ) \times ( p _ { 3 } - p _ { 1 } ) ) = 0 .$

# Notes and References for 5.3

Some classical references for the least squares problem include:

F.L. Bauer (1965). “Elimination with Weighted Row Combinations for Solving Linear Equations and Least Squares Problems,” Numer. Math. 7, 338–352.

G.H. Golub and J.H. Wilkinson (1966). “Note on the Iterative Refinement of Least Squares Solution,” Numer. Math. 9, 139–148.

A. van der Sluis (1975). “Stability of the Solutions of Linear Least Squares Problem,” Numer. Math. 23, 241–254.

The use of Gauss transformations to solve the LS problem has attracted some attention because they are cheaper to use than Householder or Givens matrices, see:

G. Peters and J.H. Wilkinson (1970). “The Least Squares Problem and Pseudo-Inverses,” Comput. J. 13, 309–16.

A.K. Cline (1973). “An Elimination Method for the Solution of Linear Least Squares Problems,” SIAM J. Numer. Anal. 10, 283–289.

R.J. Plemmons (1974). “Linear Least Squares by Elimination and MGS,” J. ACM 21, 581–585.

The seminormal equations are given by $R ^ { T } R x = A ^ { T } b$ where $A = Q R .$ . It can be shown that by solving the seminormal equations an acceptable $\mathrm { L S }$ solution is obtained if one step of fixed precision iterative improvement is performed, see:

˚A. Bj¨orck (1987). “Stability Analysis of the Method of Seminormal Equations,” Lin. Alg. Applic. 88/89, 31–48.

Survey treatments of LS perturbation theory include Lawson and Hanson (SLS), Stewart and Sun (MPT), and Bj¨orck (NMLS). See also:

P.-A. Wedin (1973). “Perturbation Theory for Pseudoinverses,” BIT 13, 217–232.

˚A. Bj¨orck (1991). “Component-wise Perturbation Analysis and Error Bounds for Linear Least Squares Solutions,” BIT 31, 238–244.

B. Wald´en, R. Karlson, J. Sun (1995). “Optimal Backward Perturbation Bounds for the Linear Least Squares Problem,” Numerical Lin. Alg. Applic. 2, 271–286.

J.-G. Sun (1996). “Optimal Backward Perturbation Bounds for the Linear Least-Squares Problem with Multiple Right-Hand Sides,” IMA J. Numer. Anal. 16, 1–11.

J.-G. Sun (1997). “On Optimal Backward Perturbation Bounds for the Linear Least Squares Problem,” BIT 37, 179–188.

R. Karlson and B. Wald´en (1997). “Estimation of Optimal Backward Perturbation Bounds for the Linear Least Squares Problem,” BIT 37, 862–869.

J.-G. Sun (1997). “On Optimal Backward Perturbation Bounds for the Linear Least Squares Problem,” BIT 37, 179–188.

M. Gu (1998). “Backward Perturbation Bounds for Linear Least Squares Problems,” SIAM J. Matrix Anal. Applic. 20, 363–372.

M. Arioli, M. Baboulin and S. Gratton (2007). “A Partial Condition Number for Linear Least Squares Problems,” SIAM J. Matrix Anal. Applic. 29, 413–433.

M. Baboulin, J. Dongarra, S. Gratton, and J. Langou (2009). “Computing the Conditioning of the Components of a Linear Least-Squares Solution,” Num. Lin. Alg. Applic. 16, 517–533.

M. Baboulin and S. Gratton (2009). “Using Dual Techniques to Derive Componentwise and Mixed Condition Numbers for a Linear Function of a Least Squares Solution,” BIT 49, 3–19.

J. Grcar (2009). “Nuclear Norms of Rank-2 Matrices for Spectral Condition Numbers of Rank Least Squares Solutions,” ArXiv:1003.2733v4.

J. Grcar (2010). “Spectral Condition Numbers of Orthogonal Projections and Full Rank Linear Least Squares Residuals,” SIAM J. Matrix Anal. Applic. 31, 2934–2949.

Practical insights into the accuracy of a computed least squares solution can be obtained by applying the condition estimation ideas of 3.5. to the R matrix in $A = Q R$ or the Cholesky factor of $\bar { A } ^ { T } \bar { A }$ should a normal equation approach be used. For a discussion of LS-specific condition estimation, see:

G.W. Stewart (1980). “The Efficient Generation of Random Orthogonal Matrices with an Application to Condition Estimators,” SIAM J. Numer. Anal. 17, 403–9.

S. Gratton (1996). “On the Condition Number of Linear Least Squares Problems in a Weighted Frobenius Norm,” BIT 36, 523–530.

C.S. Kenney, A.J. Laub, and M.S. Reese (1998). “Statistical Condition Estimation for Linear Least Squares,” SIAM J. Matrix Anal. Applic. 19, 906–923.

Our restriction to least squares approximation is not a vote against minimization in other norms. There are occasions when it is advisable to minimize $\| A x - b \| _ { p }$ for $p = 1$ and ∞. Some algorithms for doing this are described in:

A.K. Cline (1976). “A Descent Method for the Uniform Solution to Overdetermined Systems of Equations,” SIAM J. Numer. Anal. 13, 293–309.

R.H. Bartels, A.R. Conn, and C. Charalambous (1978). “On Cline’s Direct Method for Solving Overdetermined Linear Systems in the $L _ { \infty }$ Sense,” SIAM J. Numer. Anal. 15, 255–270.

T.F. Coleman and Y. Li (1992). “A Globally and Quadratically Convergent Affine Scaling Method for Linear $L _ { 1 }$ Problems,” Mathematical Programming 56, Series A, 189–222.

Y. Li (1993). $^ { 6 6 } \mathrm { A }$ Globally Convergent Method for $L _ { p }$ Problems,” SIAM J. Optim. 3, 609–629.

Y. Zhang (1993). “A Primal-Dual Interior Point Approach for Computing the $L _ { 1 }$ and $L _ { \infty }$ Solutions of Overdetermined Linear Systems,” J. Optim. Theory Applic. $^ { 7 7 , }$ 323–341.

Iterative improvement in the least squares context is discussed in:

G.H. Golub and J.H. Wilkinson (1966). “Note on Iterative Refinement of Least Squares Solutions,” Numer. Math. 9, 139–148.

˚A. Bj¨orck and G.H. Golub (1967). “Iterative Refinement of Linear Least Squares Solutions by Householder Transformation,” BIT 7, 322–337.

˚A. Bj¨orck (1967). “Iterative Refinement of Linear Least Squares Solutions I,” BIT 7, 257–278.

˚A. Bj¨orck (1968). “Iterative Refinement of Linear Least Squares Solutions II,”BIT 8, 8–30.

J. Gluchowska and A. Smoktunowicz (1999). “Solving the Linear Least Squares Problem with Very High Relative Acuracy,” Computing 45, 345–354.

J. Demmel, Y. Hida, and E.J. Riedy (2009). “Extra-Precise Iterative Refinement for Overdetermined Least Squares Problems,” ACM Trans. Math. Softw. 35, Article 28.

The following texts treat various geometric matrix problems that arise in computer graphics and vision:

A.S. Glassner (1989). An Introduction to Ray Tracing, Morgan Kaufmann, Burlington, MA.

R. Hartley and A. Zisserman (2004). Multiple View Geometry in Computer Vision, Second Edition, Cambridge University Press, New York.

M. Pharr and M. Humphreys (2010). Physically Based Rendering, from Theory to Implementation, Second Edition, Morgan Kaufmann, Burlington, MA.

For a numerical perspective, see:

W. Kahan (2008). “Computing Cross-Products and Rotations in 2- and 3-dimensional Euclidean Spaces,” http://www.cs.berkeley.edu/ wkahan/MathH110/Cross.pdf.
