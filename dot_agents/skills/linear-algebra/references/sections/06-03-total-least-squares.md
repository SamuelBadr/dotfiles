# 6.3 Total Least Squares

The problem of minimizing $\mid A x - b \mid \mid _ { 2 }$ where $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ can be recast as follows:

$$
\min _ {b + r \in \operatorname{ran} (A)} \| r \| _ {2}. \tag {6.3.1}
$$

In this problem, there is a tacit assumption that the errors are confined to the vector of observations b. If error is also present in the data matrix A, then it may be more natural to consider the problem

$$
\min _ {b + r \in \operatorname{ran} (A + E)} \| [ E \mid r ] \| _ {F}. \tag {6.3.2}
$$

This problem, discussed by Golub and Van Loan (1980), is referred to as the total least squares (TLS) problem. If a minimizing $\left[ E _ { 0 } \right| r _ { 0 } ]$ can be found for (6.3.2), then any x satisfying $( A + E _ { 0 } ) x = b + r _ { 0 }$ is called a TLS solution. However, it should be realized that (6.3.2) may fail to have a solution altogether. For example, if

$$
A = \left[ \begin{array}{l l} 1 & 0 \\ 0 & 0 \\ 0 & 0 \end{array} \right], b = \left[ \begin{array}{l} 1 \\ 1 \\ 1 \end{array} \right], E _ {\epsilon} = \left[ \begin{array}{l l} 0 & 0 \\ 0 & \epsilon \\ 0 & \epsilon \end{array} \right],
$$

then for all $\epsilon > 0 , b \in \mathsf { r a n } ( A + E \epsilon )$ . However, there is no smallest value of $\| \left[ E , r \right] \| _ { F }$ for which $b + r \in \mathsf { r a n } ( A + E )$ .

A generalization of (6.3.2) results if we allow multiple right-hand sides and use a weighted Frobenius norm. In particular, if $B \in \mathbb { R } ^ { m \times k }$ and the matrices

$$
D = \operatorname{diag} (d _ {1}, \dots , d _ {m}),
$$

$$
T = \mathrm{diag} (t _ {1}, \ldots , t _ {n + k})
$$

are nonsingular, then we are led to an optimization problem of the form

$$
\min _ {B + R \in \operatorname{ran} (A + E)} \| D [ E \mid R ] T \| _ {F} \tag {6.3.3}
$$

where $E \in \mathbb { R } ^ { m \times n }$ and $\boldsymbol { R } \in \mathbb { R } ^ { m \times k }$ . $\mathrm { I f } \ \big [ \ E _ { 0 } \ | \ R _ { 0 } \ \big ]$ solves (6.3.3), then any $\boldsymbol { X } \in \mathbb { R } ^ { n \times k }$ that satisfies

$$
(A + E _ {0}) X = (B + R _ {0})
$$

is said to be a TLS solution to (6.3.3).

In this section we discuss some of the mathematical properties of the total least squares problem and show how it can be solved using the SVD. For a more detailed introduction, see Van Huffel and Vanderwalle (1991).

# 6.3.1 Mathematical Background

The following theorem gives conditions for the uniqueness and existence of a TLS solution to the multiple-right-hand-side problem.

Theorem 6.3.1. Suppose $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { m \times k }$ and that $D = \operatorname { d i a g } ( d _ { 1 } , \dots , d _ { m } )$ and $T = \mathrm { d i a g } ( t _ { 1 } , \dots , t _ { n + k } )$ are nonsingular. Assume $m \geq n + k$ and let the $S V D o f$

$$
C = D [ A \mid B ] T = \left[ \begin{array}{c c} C _ {1} & C _ {2} \\ n & k \end{array} \right]
$$

be specified by $U ^ { T } C V = \operatorname { d i a g } ( \sigma _ { 1 } , \dots , \sigma _ { n + k } ) = \Sigma$ where U , V , and Σ are partitioned as follows:

$$
U = \left[ \begin{array}{c c} U _ {1} & U _ {2} \\ n & k \end{array} \right] \quad , \qquad V = \left[ \begin{array}{c c} V _ {1 1} & V _ {1 2} \\ V _ {2 1} & V _ {2 2} \\ n & k \end{array} \right] _ {k} ^ {n} \quad , \qquad \Sigma = \left[ \begin{array}{c c} \Sigma_ {1} & 0 \\ 0 & \Sigma_ {2} \\ n & k \end{array} \right] _ {k} ^ {n} \quad .
$$

$I f \sigma _ { n } ( C _ { 1 } ) > \sigma _ { n + 1 } ( C )$ , then the matrix $\left[ E _ { 0 } \mid R _ { 0 } \right]$ defined by

$$
D [ E _ {0} \mid R _ {0} ] T = - U _ {2} \Sigma_ {2} [ V _ {1 2} ^ {T} \mid V _ {2 2} ^ {T} ] \tag {6.3.4}
$$

solves (6.3.3). If $T _ { 1 } = \mathrm { d i a g } ( t _ { 1 } , \dots , t _ { n } )$ and $T _ { 2 } = \mathrm { d i a g } ( t _ { n + 1 } , \ldots , t _ { n + k } )$ , then the matrix

$$
X _ {T L S} = - T _ {1} V _ {1 2} V _ {2 2} ^ {- 1} T _ {2} ^ {- 1}
$$

exists and is the unique TLS solution to $( A + E _ { 0 } ) X = B + R _ { 0 }$ .

Proof. We first establish two results that follow from the assumption $\sigma _ { n } ( C _ { 1 } ) > \sigma _ { n + 1 } ( C )$ . From the equation $C V = U \Sigma$ we have

$$
C _ {1} V _ {1 2} + C _ {2} V _ {2 2} = U _ {2} \Sigma_ {2}.
$$

We wish to show that $V _ { 2 2 }$ is nonsingular. Suppose $V _ { 2 2 } x = 0$ for some unit 2-norm x. It follows from

$$
V _ {1 2} ^ {T} V _ {1 2} + V _ {2 2} ^ {T} V _ {2 2} = I
$$

that $\| \ V _ { 1 2 } x \| _ { 2 } = 1$ . But then

$$
\sigma_ {n + 1} (C) \geq \| U _ {2} \Sigma_ {2} x \| _ {2} = \| C _ {1} V _ {1 2} x \| _ {2} \geq \sigma_ {n} (C _ {1}),
$$

a contradiction. Thus, the submatrix $V _ { 2 2 }$ is nonsingular. The second fact concerns the strict separation of $\sigma _ { n } ( C )$ and $\sigma _ { n + 1 } ( C )$ . From Corollary 2.4.5, we have $\sigma _ { n } ( C ) \geq \sigma _ { n } ( C _ { 1 } )$ and so

$$
\sigma_ {n} (C) \geq \sigma_ {n} (C _ {1}) > \sigma_ {n + 1} (C).
$$

We are now set to prove the theorem. If ran $( B + R ) \subset \mathsf { r a n } ( A + E )$ , then there is an $X \ ( n { \mathrm { - b y - } } k ) \ { \mathrm { s o } } \ ( A + E ) X = B + R , { \mathrm { i . e } }$ .,

$$
\{D [ A \mid B ] T + D [ E \mid R ] T \} T ^ {- 1} \left[ \begin{array}{c} X \\ - I _ {k} \end{array} \right] = 0. \tag {6.3.5}
$$

Thus, the rank of the matrix in curly brackets is at most equal to n. By following the argument in the proof of Theorem 2.4.8, it can be shown that

$$
\| D [ E | R ] T \| _ {F} ^ {2} \geq \sum_ {i = n + 1} ^ {n + k} \sigma_ {i} (C) ^ {2}.
$$

Moreover, the lower bound is realized by setting $[ \ E \ | \ R \ ] \ = \ \big [ \ E _ { 0 } \ | \ R _ { 0 } \ \big ]$ . Using the inequality $\sigma _ { n } ( C ) > \sigma _ { n + 1 } ( C )$ , we may infer that $\left[ E _ { 0 } \mid R _ { 0 } \right]$ is the unique minimizer.

To identify the TLS solution $X _ { T L S }$ , we observe that the nullspace of

$$
\left\{D \left[ A \mid B \right] T + D \left[ E _ {0} \mid R _ {0} \right] T \right\} = U _ {1} \Sigma_ {1} \left[ V _ {1 1} ^ {T} \mid V _ {2 1} ^ {T} \right]
$$

is the range of $\left[ \begin{array} { l } { V _ { 1 2 } } \\ { V _ { 2 2 } } \end{array} \right]$ . Thus, from (6.3.5)

$$
T ^ {- 1} \left[ \begin{array}{c} X \\ - I _ {k} \end{array} \right] = \left[ \begin{array}{c} V _ {1 2} \\ V _ {2 2} \end{array} \right] S
$$

for some k-by-k matrix $S .$ From the equations $T _ { 1 } ^ { - 1 } X = V _ { 1 2 } S$ and $- T _ { 2 } ^ { - 1 } = V _ { 2 2 } S$ we see that $S = - V _ { 2 2 } ^ { - 1 } T _ { 2 } ^ { - 1 }$ and so

$$
X = T _ {1} V _ {1 2} S = - T _ {1} V _ {1 2} V _ {2 2} ^ {- 1} T _ {2} ^ {- 1} = X _ {\mathrm{TLS}}. \quad \square
$$

Note from the thin CS decomposition (Theorem 2.5.2) that

$$
\parallel X \parallel_ {\tau} ^ {2} = \parallel V _ {1 2} V _ {2 2} ^ {- 1} \parallel_ {2} ^ {2} = \frac {1 - \sigma_ {k} (V _ {2 2}) ^ {2}}{\sigma_ {k} (V _ {2 2}) ^ {2}}
$$

where we define the $\ " \tau \mathrm { - n o r m } ^ { \prime \mathrm { - } }$ on $\mathbb { R } ^ { n \times k }$ by $\lVert Z \rVert _ { \tau } = \lVert \boldsymbol { T } _ { 1 } ^ { - 1 } Z \boldsymbol { T } _ { 2 } \rVert _ { 2 } .$

If $\sigma _ { n } ( C _ { 1 } ) = \sigma _ { n + 1 } ( C )$ , then the solution procedure implicit in the above proof is problematic. The TLS problem may have no solution or an infinite number of solutions. See §6.3.4 for suggestions as to how one might proceed.

# 6.3.2 Solving the Single Right Hand Side Case

We show how to maximize $\sigma _ { k } ( V _ { 2 2 } )$ in the important $k = 1$ case. Suppose the singular values of C satisfy $\sigma _ { n - p } > \sigma _ { n - p + 1 } = \cdot \cdot \cdot = \sigma _ { n + 1 }$ and let $V = { \left[ \begin{array} { l } { v _ { 1 } \mid \cdots \mid v _ { n + 1 } } \end{array} \right] }$ b e a column partitioning of V . If $\widetilde { Q }$ is a Householder matrix such that

$$
V (:, n + 1 - p: n + 1) \widetilde {Q} = \left[ \begin{array}{c c} W & z \\ 0 & \alpha \end{array} \right] _ {1} ^ {n},
$$

then the last column of this matrix has the largest $( n + 1 ) \mathrm { s t }$ component of all the vectors in span $\{ v _ { n + 1 - p } , \ldots , v _ { n + 1 } \}$ . If $\alpha = 0$ , then the TLS problem has no solution. Otherwise

$$
x _ {\mathrm{TLS}} = - T _ {1} z / (t _ {n + 1} \alpha).
$$

Moreover,

$$
\left[ \begin{array}{c c} I _ {n - 1} & 0 \\ 0 & \widetilde {Q} \end{array} \right] U ^ {T} (D [ A | b ] T) V \left[ \begin{array}{c c} I _ {n - p} & 0 \\ 0 & \widetilde {Q} \end{array} \right] = \Sigma
$$

and so

$$
D \left[ E _ {0} \mid r _ {0} \right] T = - D \left[ A \mid b \right] T \left[ \begin{array}{c} z \\ \alpha \end{array} \right] \left[ z ^ {T} \mid \alpha \right].
$$

Overall, we have the following algorithm:

Algorithm 6.3.1 Given $A \in \mathbb { R } ^ { m \times n } ( m > n ) , b \in \mathbb { R } ^ { m }$ , nonsingular $D = \operatorname* { d i a g } ( d _ { 1 } , \dots , d _ { m } )$ and nonsingular $T = \mathrm { d i a g } ( t _ { 1 } , \dots , t _ { n + 1 } )$ , the following algorithm computes (if possible) a vector $\boldsymbol { x } _ { \mathrm { T L S } } \in \mathbb { R } ^ { n }$ such that $( A + E _ { 0 } ) x _ { \mathrm { T L S } } = ( b + r _ { 0 } )$ and $\parallel D [ E _ { 0 } \mid r _ { 0 } ] T \parallel _ { F }$ is minimal.

Compute the SVD $U ^ { T } ( D [ \ : A \ : | \ : b \ : | T ) V = \ : \mathrm { d i a g } ( \sigma _ { 1 } , \ldots , \sigma _ { n + 1 } )$ and save V .

Determine p such that $\sigma _ { 1 } \geq \cdot \cdot \cdot \geq \sigma _ { n - p } > \sigma _ { n - p + 1 } = \cdot \cdot \cdot = \sigma _ { n + 1 } $ .

Compute a Householder P such that if $\tilde { V } = V P$ , then $\tilde { V } ( n + 1 , n - p + 1 { : } n ) = 0$ .

if $\tilde { v } _ { n + 1 , n + 1 } \neq 0$

for $i = 1 { : } n$

$$
x _ {i} = - t _ {i} \tilde {v} _ {i, n + 1} / (t _ {n + 1} \tilde {v} _ {n + 1, n + 1})
$$

end

$$
x _ {\mathrm{TLS}} = x
$$

end

This algorithm requires about $2 m n ^ { 2 } + 1 2 n ^ { 3 }$ flops and most of these are associated with the SVD computation.

# 6.3.3 A Geometric Interpretation

It can be shown that the TLS solution $x _ { T L S }$ minimizes

$$
\psi (x) = \sum_ {i = 1} ^ {m} d _ {i} ^ {2} \left(\frac {\left| a _ {i} ^ {T} x - b _ {i} \right| ^ {2}}{x ^ {T} T _ {1} ^ {- 2} x + t _ {n + 1} ^ {- 2}}\right) \tag {6.3.6}
$$

where $a _ { i } ^ { T }$ is the ith row of A and $b _ { i }$ is the ith component of b. A geometrical interpretation of the TLS problem is made possible by this observation. Indeed,

$$
\delta_ {i} = \frac {| a _ {i} ^ {T} x - b _ {i} | ^ {2}}{x ^ {T} T _ {1} ^ {- 2} x + t _ {n + 1} ^ {- 2}}
$$

is the square of the distance from

$$
\left[ \begin{array}{l} a _ {i} \\ b _ {i} \end{array} \right] \in \mathbb {R} ^ {n + 1}
$$

to the nearest point in the subspace

$$
P _ {x} = \left\{\left[ \begin{array}{c} a \\ b \end{array} \right]: a \in \mathbb {R} ^ {n}, b \in \mathbb {R}, b = x ^ {T} a \right\}
$$

where the distance in $\mathbb { R } ^ { n + 1 }$ is measured by the norm $\| z \| = \| T z \| _ { 2 }$ . The TLS problem is essentially the problem of orthogonal regression, a topic with a long history. See Pearson (1901) and Madansky (1959).

# 6.3.4 Variations of the Basic TLS Problem

We briefly mention some modified TLS problems that address situations when additional constraints are imposed on the optimizing E and R and the associated TLS solution.

In the restricted TLS problem, we are given $A \in \mathbb { R } ^ { m \times n } , B \in \mathbb { R } ^ { m \times k } , P _ { 1 } \in \mathbb { R } ^ { m \times q }$ , and $P _ { 2 } \in \mathbb { R } ^ { n + k \times r }$ , and solve

$$
\min _ {B + R \subset \operatorname{ran} (A + E)} \| P _ {1} ^ {T} [ E \mid R ] P _ {2} \| _ {F}. \tag {6.3.7}
$$

We assume that $q \leq m$ and $r \leq n + k$ . An important application arises if some of the columns of A are error-free. For example, if the first s columns of A are error-free, then it makes sense to force the optimizing E to satisfy $E ( : , 1 : s ) = 0$ . This goal is achieved by setting $P _ { 1 } = I _ { m }$ and $P _ { 2 } = I _ { m + k } ( : , s + 1 { : } n + k )$ in the restricted TLS problem.

If a particular TLS problem has no solution, then it is referred to as a nongeneric TLS problem. By adding a constraint it is possible to produce a meaningful solution. For example, let $U ^ { T } [ \stackrel { \bar { A ( } ) } { A ( } b ] V = \Sigma$ be the SVD and let p be the largest index so $V ( n + 1 , p ) \neq 0$ . It can be shown that the problem

$$
\begin{array}{l} \min \quad \left\| \left[ E \mid r \right] \right\| _ {F} \\ (A + E) x = b + r \tag {6.3.8} \\ [ E \mid r ] V (:, p + 1: n + 1) = 0 \\ \end{array}
$$

has a solution $\left[ E _ { 0 } \right| r _ { 0 } ]$ and the nongeneric TLS solution satisfies $( A + E _ { 0 } ) x + b + r _ { 0 }$ . See Van Huffel (1992).

In the regularized TLS problem additional constraints are imposed to ensure that the solution x is properly constrained/smoothed:

$$
\begin{array}{l} \min \quad \| [ E | r ] \| _ {F}. \\ (A + E) x = b + r \tag {6.3.9} \\ \| L x \| _ {2} \leq \delta \\ \end{array}
$$

The matrix $\boldsymbol { L } \in \mathbb { R } ^ { n \times n }$ could be the identity or a discretized second-derivative operator. The regularized TLS problem leads to a Lagrange multiplier system of the form

$$
(A ^ {T} A + \lambda_ {1} I + \lambda_ {2} L ^ {T} L) x = A ^ {T} b.
$$

See Golub, Hansen, and O’Leary (1999) for more details. Another regularization approach involves setting the small singular values of $\left[ A \mid b \right]$ to zero. This is the truncated TLS problem discussed in Fierro, Golub, Hansen, and O’Leary (1997).

# Problems

P6.3.1 Consider the TLS problem (6.3.2) with nonsingular D and T . (a) Show that if rank $( A ) < n$ , then (6.3.2) has a solution if and only if $b \in \mathsf { r a n } ( A )$ . (b) Show that if rank $( A ) = n$ , then (6.3.2) has no solution if $A ^ { T } D ^ { 2 } b = 0$ and $\left| t _ { n + 1 } \right| \left\| \ D b \right\| _ { 2 } \geq \sigma _ { n } ( D A T _ { 1 } )$ where $T _ { 1 } = \mathrm { d i a g } ( t _ { 1 } , \dots , t _ { n } )$ .

P6.3.2 Show that if $C = D [ A \mid b ] T = [ A _ { 1 } \mid d ]$ and $\sigma _ { n } ( C ) > \sigma _ { n + 1 } ( C )$ , then $x _ { T L S }$ satisfies

$$
(A _ {1} ^ {T} A _ {1} - \sigma_ {n + 1} (C) ^ {2} I) x _ {\mathrm{TLS}} = A _ {1} ^ {T} d.
$$

Appreciate this as a “negatively shifted” system of normal equations.

P6.3.3 Show how to solve (6.3.2) with the added constraint that the first p columns of the minimizing E are zero. Hint: Compute the QR factorization of $A ( : , 1 { : } p )$ .

P6.3.4 Show how to solve (6.3.3) given that D and $T$ are general nonsingular matrices.

P6.3.5 Verify Equation (6.3.6).

P6.3.6 If $A \in \mathbb { R } ^ { m \times n }$ has full column rank and $B \in \mathbb { R } ^ { p \times n }$ has full row rank, show how to minimize

$$
f (x) = \frac {\parallel A x - b \parallel_ {2} ^ {2}}{1 + x ^ {T} x}
$$

subject to the constraint that $B x = 0 .$

P6.3.7 In the data least squares problem, we are given $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ and minimize $\Vert E \Vert _ { F }$ subject to the constraint that $b \in \mathsf { r a n } ( A + E )$ . Show how to solve this problem. See Paige and Strakoˇs (2002b).

# Notes and References for 6.3

Much of this section is based on:

G.H. Golub and C.F. Van Loan (1980). “An Analysis of the Total Least Squares Problem,” SIAM J. Numer. Anal. 17, 883–93.

The idea of using the SVD to solve the TLS problem is set forth in:

G.H. Golub and C. Reinsch (1970). “Singular Value Decomposition and Least Squares Solutions,” Numer. Math. 14, 403–420.

G.H. Golub (1973). “Some Modified Matrix Eigenvalue Problems,” SIAM Review 15, 318–334.

The most comprehensive treatment of the TLS problem is:

S. Van Huffel and J. Vandewalle (1991). The Total Least Squares Problem: Computational Aspects and Analysis, SIAM Publications, Philadelphia, PA.

There are two excellent conference proceedings that cover just about everything you would like to know about TLS algorithms, generalizations, applications, and the associated statistical foundations:

S. Van Huffel (ed.) (1996). Recent Advances in Total Least Squares Techniques and Errors in Variables Modeling, SIAM Publications, Philadelphia, PA.

S. Van Huffel and P. Lemmerling (eds.) (2002) Total Least Squares and Errors-in-Variables Modeling: Analysis, Algorithms, and Applications, Kluwer Academic, Dordrecht, The Netherlands.

TLS is but one approach to the errors-in-variables problem, a subject that has a long and important history in statistics:

K. Pearson (1901). “On Lines and Planes of Closest Fit to Points in Space,” Phil. Mag. 2, 559–72.

A. Wald (1940). “The Fitting of Straight Lines if Both Variables are Subject to Error,” Annals of Mathematical Statistics 11, 284–300.

G.W. Stewart (2002). “Errors in Variables for Numerical Analysts,” in Recent Advances in Total Least Squares Techniques and Errors-in-Variables Modelling, S. Van Huffel (ed.), SIAM Publications, Philadelphia PA, pp. 3–10,

In certain settings there are more economical ways to solve the TLS problem than the Golub-Kahan-Reinsch SVD algorithm:

S. Van Huffel and H. Zha (1993). “An Efficient Total Least Squares Algorithm Based On a Rank-Revealing Two-Sided Orthogonal Decomposition,” Numer. Alg. 4, 101–133.

˚A. Bj¨orck, P. Heggernes, and P. Matstoms (2000). “Methods for Large Scale Total Least Squares Problems,” SIAM J. Matrix Anal. Applic. 22, 413–429.


---

<!-- golub_350_399 -->

H. Guo and R.A. Renaut (2005). “Parallel Variable Distribution for Total Least Squares,” Num. Lin. Alg. 12, 859–876.   
The condition of the TLS problem is analyzed in:   
M. Baboulin and S. Gratton (2011). “A Contribution to the Conditioning of the Total Least-Squares Problem,” SIAM J. Matrix Anal. Applic. 32, 685–699.   
Efforts to connect the LS and TLS paradigms have lead to nice treatments that unify the presentation of both approaches:   
B.D. Rao (1997). “Unified Treatment of LS, TLS, and Truncated SVD Methods Using a Weighted TLS Framework,” in Recent Advances in Total Least Squares Techniques and Errors-in-Variables Modelling, S. Van Huffel (ed.), SIAM Publications, Philadelphia, PA., pp. 11–20.   
C.C. Paige and Z. Strakoˇs (2002a). “Bounds for the Least Squares Distance Using Scaled Total Least Squares,” Numer. Math. 91, 93–115.   
C.C. Paige and Z. Strakoˇs (2002b). “Scaled Total Least Squares Fundamentals,” Numer. Math. 91, 117–146.   
X.-W. Chang, G.H. Golub, and C.C. Paige (2008). “Towards a Backward Perturbation Analysis for Data Least Squares Problems,” SIAM J. Matrix Anal. Applic. 30, 1281–1301.   
X.-W. Chang and D. Titley-Peloquin (2009). “Backward Perturbation Analysis for Scaled Total Least-Squares,” Num. Lin. Alg. Applic. 16, 627–648.   
For a discussion of the situation when there is no TLS solution or when there are multiple solutions, see:   
S. Van Huffel and J. Vandewalle (1988). “Analysis and Solution of the Nongeneric Total Least Squares Problem,” SIAM J. Matrix Anal. Appl. 9, 360–372.   
S. Van Huffel (1992). “On the Significance of Nongeneric Total Least Squares Problems,” SIAM J. Matrix Anal. Appl. 13, 20–35.   
M. Wei (1992). “The Analysis for the Total Least Squares Problem with More than One Solution,” SIAM J. Matrix Anal. Appl. 13, 746–763.   
For a treatment of the multiple right hand side TLS problem, see:   
I. Hn˘etynkov˜a, M. Ple˘singer, D.M. Sima, Z. Strako˘s, and S. Van Huffel (2011). “The Total Least Squares Problem in AX  B: A New Classification with the Relationship to the Classical Works,” SIAM J. Matrix Anal. Applic. 32, 748–770.   
If some of the columns of A are known exactly then it is sensible to force the TLS perturbation matrix E to be zero in the same columns. Aspects of this constrained TLS problem are discussed in:   
J.W. Demmel (1987). “The Smallest Perturbation of a Submatrix which Lowers the Rank and Constrained Total Least Squares Problems,” SIAM J. Numer. Anal. 24, 199–206.   
S. Van Huffel and J. Vandewalle (1988). “The Partial Total Least Squares Algorithm,” J. Comput. App. Math. 21, 333–342.   
S. Van Huffel and J. Vandewalle (1989). “Analysis and Properties of the Generalized Total Least Squares Problem AX ≈ B When Some or All Columns in A are Subject to Error,” SIAM J. Matrix Anal. Applic. 10, 294–315.   
S. Van Huffel and H. Zha (1991). “The Restricted Total Least Squares Problem: Formulation, Algorithm, and Properties,” SIAM J. Matrix Anal. Applic. 12, 292–309.   
C.C. Paige and M. Wei (1993). “Analysis of the Generalized Total Least Squares Problem AX = B when Some of the Columns are Free of Error,” Numer. Math. 65, 177–202.   
Another type of constraint that can be imposed in the TLS setting is to insist that the optimum perturbation of A have the same structure as A. For examples and related strategies, see:   
J. Kamm and J.G. Nagy (1998). “A Total Least Squares Method for Toeplitz Systems of Equations,” BIT 38, 560–582.   
P. Lemmerling, S. Van Huffel, and B. De Moor (2002). “The Structured Total Least Squares Approach for Nonlinearly Structured Matrices,” Num. Lin. Alg. 9, 321–332.   
P. Lemmerling, N. Mastronardi, and S. Van Huffel (2003). “Efficient Implementation of a Structured Total Least Squares Based Speech Compression Method,” Lin. Alg. Applic. 366, 295–315.   
N. Mastronardi, P. Lemmerling, and S. Van Huffel (2004). “Fast Regularized Structured Total Least Squares Algorithm for Solving the Basic Deconvolution Problem,” Num. Lin. Alg. 12, 201–209.

I. Markovsky, S. Van Huffel, and R. Pintelon (2005). “Block-Toeplitz/Hankel Structured Total Least Squares,” SIAM J. Matrix Anal. Applic. 26, 1083–1099.   
A. Beck and A. Ben-Tal (2005). “A Global Solution for the Structured Total Least Squares Problem with Block Circulant Matrices,” SIAM J. Matrix Anal. Applic. 27, 238–255.   
H. Fu, M.K. Ng, and J.L. Barlow (2006). “Structured Total Least Squares for Color Image Restoration,” SIAM J. Sci. Comput. 28, 1100–1119.   
As in the least squares problem, there are techniques that can be used to regularlize an otherwise “wild” TLS solution:   
R.D. Fierro and J.R. Bunch (1994). “Collinearity and Total Least Squares,” SIAM J. Matrix Anal. Applic. 15, 1167–1181.   
R.D. Fierro, G.H. Golub, P.C. Hansen and D.P. O’Leary (1997). “Regularization by Truncated Total Least Squares,” SIAM J. Sci. Comput. 18, 1223–1241.   
G.H. Golub, P.C. Hansen, and D.P. O’Leary (1999). “Tikhonov Regularization and Total Least Squares,” SIAM J. Matrix Anal. Applic. 21, 185–194.   
R.A. Renaut and H. Guo (2004). “Efficient Algorithms for Solution of Regularized Total Least Squares,” SIAM J. Matrix Anal. Applic. 26, 457–476.   
D.M. Sima, S. Van Huffel, and G.H. Golub (2004). “Regularized Total Least Squares Based on Quadratic Eigenvalue Problem Solvers,” BIT 44, 793–812.   
N. Mastronardi, P. Lemmerling, and S. Van Huffel (2005). “Fast Regularized Structured Total Least Squares Algorithm for Solving the Basic Deconvolution Problem,” Num. Lin. Alg. Applic. 12, 201–209.   
S. Lu, S.V. Pereverzev, and U. Tautenhahn (2009). “Regularized Total Least Squares: Computational Aspects and Error Bounds,” SIAM J. Matrix Anal. Applic. 31, 918–941.   
Finally, we mention an interesting TLS problem where the solution is subject to a unitary constraint:   
K.S. Arun (1992). “A Unitarily Constrained Total Least Squares Problem in Signal Processing,” SIAM J. Matrix Anal. Applic. 13, 729–745.
