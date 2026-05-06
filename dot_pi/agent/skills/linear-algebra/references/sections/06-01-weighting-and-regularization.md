# 6.1 Weighting and Regularization

We consider two basic modifications to the linear least squares problem. The first concerns how much each equation “counts” in the $\parallel A x - b \parallel _ { 2 }$ minimization. Some equations may be more important than others and there are ways to produce approximate minimzers that reflect this. Another situation arises when A is ill-conditioned. Instead of minimizing $\parallel A x - b \parallel _ { 2 }$ with a possibly wild, large norm x-vector, we settle for a predictor Ax in which x is “nice” according to some regularizing metric.

# 6.1.1 Row Weighting

In ordinary least squares, the minimization of $\parallel A x - b \parallel _ { 2 }$ amounts to minimizing the sum of the squared discrepancies in each equation:

$$
\parallel A x - b \parallel^ {2} = \sum_ {i = 1} ^ {m} \left(a _ {i} ^ {T} x - b _ {i}\right) ^ {2}.
$$

We assume that $A \in \mathbb { R } ^ { m \times n } , \ b \in \mathbb { R } ^ { m }$ , and $a _ { i } ^ { T } = A ( i , : )$ . In the weighted least squares problem the discrepancies are scaled and we solve

$$
\min _ {x \in \mathbb {R} ^ {n}} \| D (A x - b) \| ^ {2} = \min _ {x \in \mathbb {R} ^ {n}} \sum_ {i = 1} ^ {m} d _ {i} ^ {2} \left(a _ {i} ^ {T} x - b _ {i}\right) ^ {2} \tag {6.1.1}
$$

where $D = \operatorname* { d i a g } ( d _ { 1 } , \dots , d _ { m } )$ is nonsingular. Note that if $x _ { D }$ minimizes this summation, then it minimizes $\lVert \tilde { A } \boldsymbol { x } - \tilde { \boldsymbol { b } } \rVert _ { 2 }$ where $\widetilde { A } = D A$ and $\tilde { b } = D b$ . Although there can be numerical issues associated with disparate weight values, it is generally possible to solve the weighted least squares problem by applying any Chapter 5 method to the “tilde problem.” For example, if A has full column rank and we apply the method of normal equations, then we are led to the following positive definite system:

$$
(A ^ {T} D ^ {2} A) x _ {D} = A ^ {T} D ^ {2} b. \tag {6.1.2}
$$

Subtracting the unweighted system $A ^ { T } A x _ { L S } = A ^ { T } b$ we see that

$$
x _ {D} - x _ {L S} = (A ^ {T} D ^ {2} A) ^ {- 1} A ^ {T} (D ^ {2} - I) (b - A x _ {L S}). \tag {6.1.3}
$$

Note that weighting has less effect if b is almost in the range of A.

At the component level, increasing $d _ { k }$ relative to the other weights stresses the importance of the kth equation and the resulting residual $r = b - A x _ { D }$ tends to be smaller in that component. To make this precise, define

$$
D (\delta) = \mathrm{diag} (d _ {1}, \ldots , d _ {k - 1}, d _ {k} \sqrt {1 + \delta}, d _ {k + 1}, \ldots , d _ {m})
$$

where $\delta > - 1$ . Assume that $x ( \delta )$ minimizes $\parallel D ( \delta ) ( A x - b ) \parallel _ { 2 }$ and set

$$
r _ {k} (\delta) = e _ {k} ^ {T} (b - A x (\delta)) = b _ {k} - a _ {k} ^ {T} (A ^ {T} D (\delta) ^ {2} A) ^ {- 1} A ^ {T} D (\delta) ^ {2} b
$$

where $e _ { k } = I _ { m } ( : , k )$ . We show that the penalty for disagreement between $a _ { k } ^ { T } x$ and $b _ { k }$ increases with δ. Since

$$
\frac {d}{d \delta} \left[ D (\delta) ^ {2} \right] = d _ {k} ^ {2} e _ {k} e _ {k} ^ {T}
$$

and

$$
\frac {d}{d \delta} \left[ (A ^ {T} D (\delta) ^ {2} A) ^ {- 1} \right] = - (A ^ {T} D (\delta) ^ {2} A) ^ {- 1} (A ^ {T} (d _ {k} ^ {2} e _ {k} e _ {k} ^ {T}) A) (A ^ {T} D (\delta) ^ {2} A) ^ {- 1},
$$

it can be shown that

$$
\frac {d}{d \delta} r _ {k} (\delta) = - d _ {k} ^ {2} \left(a _ {k} ^ {T} (A ^ {T} D (\delta) ^ {2} A) ^ {- 1} a _ {k}\right) r _ {k} (\delta). \tag {6.1.4}
$$

Assuming that A has full rank, the matrix $( A ^ { T } D ( \delta ) A ) ^ { - 1 }$ is positive definite and so

$$
\frac {d}{d \delta} [ r _ {k} (\delta) ^ {2} ] = 2 r _ {k} (\delta) \cdot \frac {d}{d \delta} r _ {k} (\delta) = - 2 d _ {k} ^ {2} \left(a _ {k} ^ {T} (A ^ {T} D (\delta) ^ {2} A) ^ {- 1} a _ {k}\right) r _ {k} (\delta) ^ {2} <   0.
$$

It follows that $| r _ { k } ( \delta ) |$ is a monotone decreasing function of δ. Of course, the change in $r _ { k }$ when all the weights are varied at the same time is much more complicated.

Before we move on to a more general type of row weighting, we mention that (6.1.1) can be framed as a symmetric indefinite linear system. In particular, if

$$
\left[ \begin{array}{c c} D ^ {- 2} & A \\ A ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} r \\ x \end{array} \right] = \left[ \begin{array}{l} b \\ 0 \end{array} \right], \tag {6.1.5}
$$

then x minimizes (6.1.1). Compare with (5.3.20).

# 6.1.2 Generalized Least Squares

In statistical data-fitting applications, the weights in (6.1.1) are often chosen to increase the relative importance of accurate measurements. For example, suppose the vector of observations b has the form $b _ { \mathrm { t r u e } } + \Delta$ where $\Delta _ { i }$ is normally distributed with mean zero and standard deviation $\sigma _ { i }$ . If the errors are uncorrelated, then it makes statistical sense to minimize (6.1.1) with $d _ { i } = 1 / \sigma _ { i }$ .

In more general estimation problems, the vector b is related to x through the equation

$$
b = A x + w \tag {6.1.6}
$$

where the noise vector w has zero mean and a symmetric positive definite covariance matrix $\sigma ^ { 2 } W$ . Assume that W is known and that $W = \bar { B } B ^ { T }$ for some $B \in \mathbb { R } ^ { m \times m }$ . The matrix B might be given or it might be $W \mathrm { { s } }$ Cholesky triangle. In order that all the equations in (6.1.6) contribute equally to the determination of x, statisticians frequently solve the LS problem

$$
\min _ {x \in \mathbb {R} ^ {n}} \| B ^ {- 1} (A x - b) \| _ {2}. \tag {6.1.7}
$$

An obvious computational approach to this problem is to form $\widetilde { A } = B ^ { - 1 } A$ and $\tilde { b } = B ^ { - 1 } b$ and then apply any of our previous techniques to minimize $\parallel \tilde { A } x - \tilde { b } \parallel _ { 2 }$ . Unfortunately, if B is ill-conditioned, then x will be poorly determined by such a procedure.

A more stable way of solving (6.1.7) using orthogonal transformations has been suggested by Paige (1979a, 1979b). It is based on the idea that (6.1.7) is equivalent to the generalized least squares problem,

$$
\min _ {b = A x + B v} v ^ {T} v. \tag {6.1.8}
$$

Notice that this problem is defined even if A and B are rank deficient. Although in the Paige technique can be applied when this is the case, we shall describe it under the assumption that both these matrices have full rank.

The first step is to compute the QR factorization of A:

$$
Q ^ {T} A = \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right], \qquad Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ n & m - n \end{array} \right].
$$

Next, an orthogonal matrix $\boldsymbol { Z } \in \mathbb { R } ^ { m \times m }$ is determined such that

$$
(Q _ {2} ^ {T} B) Z = \left[ \begin{array}{c c} 0 & S \\ n & m - n \end{array} \right], \qquad Z = \left[ \begin{array}{c c} Z _ {1} & Z _ {2} \\ n & m - n \end{array} \right]
$$

where S is upper triangular. With the use of these orthogonal matrices, the constraint in (6.1.8) transforms to

$$
\left[ \begin{array}{c} Q _ {1} ^ {T} b \\ Q _ {2} ^ {T} b \end{array} \right] = \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right] x + \left[ \begin{array}{c c} Q _ {1} ^ {T} B Z _ {1} & Q _ {1} ^ {T} B Z _ {2} \\ 0 & S \end{array} \right] \left[ \begin{array}{c} Z _ {1} ^ {T} v \\ Z _ {2} ^ {T} v \end{array} \right].
$$

The bottom half of this equation determines v while the top half prescribes x:

$$
S u = Q _ {2} ^ {T} b, \quad v = Z _ {2} u, \tag {6.1.9}
$$

$$
R _ {1} x = Q _ {1} ^ {T} b - \left(Q _ {1} ^ {T} B Z _ {1} Z _ {1} ^ {T} + Q _ {1} ^ {T} B Z _ {2} Z _ {2} ^ {T}\right) v = Q _ {1} ^ {T} b - Q _ {1} ^ {T} B Z _ {2} u. \tag {6.1.10}
$$

The attractiveness of this method is that all potential ill-conditioning is concentrated in the triangular systems (6.1.9) and (6.1.10). Moreover, Paige (1979b) shows that the above procedure is numerically stable, something that is not true of any method that explicitly forms $B ^ { - 1 } A$ .

# 6.1.3 A Note on Column Weighting

Suppose $G \in \mathbb { R } ^ { n \times n }$ is nonsingular and define the G-norm $\| \cdot \| _ { G }$ on $\mathbb { R } ^ { n }$ by

$$
\left\| z \right\| _ {G} = \left\| G z \right\| _ {2}.
$$

If $\boldsymbol { A } \in \mathbb { R } ^ { m \times n } , \boldsymbol { b } \in \mathbb { R } ^ { m }$ , and we compute the minimum 2-norm solution $y _ { L S }$ to

$$
\min _ {x \in \mathbb {R} ^ {n}} \| (A G ^ {- 1}) y - b \| _ {2}  ,
$$

then $x _ { G } = G ^ { - 1 } y _ { L S }$ is a minimizer of $\parallel A x - b \parallel _ { 2 }$ . If rank $( A ) < n$ , then within the set of minimizers, $x _ { G }$ has the smallest G-norm.

The choice of G is important. Sometimes its selection can be based upon a priori knowledge of the uncertainties in A. On other occasions, it may be desirable to normalize the columns of A by setting

$$
G = G _ {0} \equiv \operatorname{diag} (\| A (:, 1) \| _ {2}, \dots , \| A (:, n) \| _ {2}).
$$

Van der Sluis (1969) has shown that with this choice, $\kappa _ { 2 } ( A G ^ { - 1 } )$ is approximately minimized. Since the computed accuracy of $y _ { L S }$ depends on $\kappa _ { 2 } ( A G ^ { - 1 } )$ , a case can be made for setting $G = G _ { 0 }$ .

We remark that column weighting affects singular values. Consequently, a scheme for determining numerical rank may not return the same estimate when applied to A and $A G ^ { - 1 }$ . See Stewart (1984).

# 6.1.4 Ridge Regression

In the ridge regression problem we are given $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ and proceed to solve

$$
\min _ {x} \left\| \left[ \begin{array}{c} A \\ \sqrt {\lambda} I \end{array} \right] x - \left[ \begin{array}{l} b \\ 0 \end{array} \right] \right\| _ {2} ^ {2} = \min _ {x} \| A x - b \| _ {2} ^ {2} + \lambda \| x \| _ {2} ^ {2}. \tag {6.1.11}
$$

where the value of the ridge parameter λ is chosen to “shape” the solution $x = x ( \lambda )$ （2号 in some meaningful way. Notice that the normal equation system for this problem is given by

$$
(A ^ {T} A + \lambda I) x = A ^ {T} b. \tag {6.1.12}
$$

It follows that if

$$
A = U \Sigma V ^ {T} = \sum_ {i = 1} ^ {r} \sigma_ {i} u _ {i} v _ {i} ^ {T} \tag {6.1.13}
$$

is the SVD of A, then (6.1.12) converts to

$$
(\Sigma^ {T} \Sigma + \lambda I _ {n}) (V ^ {T} x) = \Sigma^ {T} U ^ {T} b
$$

and so

$$
x (\lambda) = \sum_ {i = 1} ^ {r} \frac {\sigma_ {i} u _ {i} ^ {T} b}{\sigma_ {i} ^ {2} + \lambda} v _ {i}. \tag {6.1.14}
$$

By inspection, it is clear that

$$
\lim _ {\lambda \to 0} x (\lambda) = x _ {L S}
$$

and $\parallel x ( \lambda ) \parallel _ { 2 }$ is a monotone decreasing function of λ. These two facts show how an ill-conditioned least squares solution can be regularized by judiciously choosing λ. The idea is to get sufficiently close to $x _ { L S }$ subject to the constraint that the norm of the ridge regression minimzer $x ( \lambda )$ is sufficiently modest. Regularization in this context is all about the intelligent balancing of these two tensions.

The ridge parameter can also be chosen with an eye toward balancing the “impact” of each equation in the overdetermined system $A x = b$ . We describe a λ-selection procedure due to Golub, Heath, and Wahba (1979). Set

$$
D _ {k} = I - e _ {k} e _ {k} ^ {T} = \operatorname{diag} (1, \dots , 1, 0, 1, \dots , 1) \in \mathbb {R} ^ {m \times m}
$$

and let $x _ { k } ( \lambda )$ solve

$$
\min _ {x \in \mathbb {R} ^ {n}} \| D _ {k} (A x - b) \| _ {2} ^ {2} + \lambda \| x \| _ {2} ^ {2}. \tag {6.1.15}
$$

Thus, $x _ { k } ( \lambda )$ is the solution to the ridge regression problem with the kth row of A and kth component of b deleted, i.e., the kth equation in the overdetermined system Ax = b is deleted. Now consider choosing λ so as to minimize the cross-validation weighted square error $C ( \lambda )$ defined by

$$
C (\lambda) = \frac {1}{m} \sum_ {k = 1} ^ {m} w _ {k} (a _ {k} ^ {T} x _ {k} (\lambda) - b _ {k}) ^ {2}.
$$

Here, $w _ { 1 } , \ldots , w _ { m }$ are nonnegative weights and $a _ { k } ^ { T }$ is the kth row of A. Noting that

$$
\parallel A x _ {k} (\lambda) - b \parallel_ {2} ^ {2} = \parallel D _ {k} (A x _ {k} (\lambda) - b) \parallel_ {2} ^ {2} + (a _ {k} ^ {T} x _ {k} (\lambda) - b _ {k}) ^ {2},
$$

we see that $\left( a _ { k } ^ { T } x _ { k } ( \lambda ) - b _ { k } \right) ^ { 2 }$ is the increase in the sum of squares that results when the kth row is “reinstated.” Minimizing $C ( \lambda )$ is tantamount to choosing λ such that the final model is not overly dependent on any one experiment.

A more rigorous analysis can make this statement precise and also suggest a method for minimizing $C ( \lambda )$ . Assuming that $\lambda > 0$ , an algebraic manipulation shows that

$$
x _ {k} (\lambda) = x (\lambda) + \frac {a _ {k} ^ {T} x (\lambda) - b _ {k}}{1 - z _ {k} ^ {T} a _ {k}} z _ {k} \tag {6.1.16}
$$

where $\begin{array} { l l l } { z _ { k } } & { = } & { ( A ^ { T } A + \lambda I ) ^ { - 1 } { a } _ { k } } \end{array}$ and $x ( \lambda ) ~ = ~ ( A ^ { T } A + \lambda I ) ^ { - 1 } A ^ { T } b$ . Applying $- a _ { k } ^ { T }$ to (6.1.16) and then adding $b _ { k }$ to each side of the resulting equation gives

$$
r _ {k} = b _ {k} - a _ {k} ^ {T} x _ {k} (\lambda) = \frac {e _ {k} ^ {T} (I - A (A ^ {T} A + \lambda I) ^ {- 1} A ^ {T}) b}{e _ {k} ^ {T} (I - A (A ^ {T} A + \lambda I) ^ {- 1} A ^ {T}) e _ {k}}. \tag {6.1.17}
$$

Noting that the residual $\boldsymbol { r } = [ r _ { 1 } , \ldots , r _ { m } ] ^ { T } = \boldsymbol { b } - \boldsymbol { A } \boldsymbol { x } ( \boldsymbol { \lambda } )$ is given by the formula

$$
r = [ I - A (A ^ {T} A + \lambda I) ^ {- 1} A ^ {T} ] b,
$$

we see that

$$
C (\lambda) = \frac {1}{m} \sum_ {k = 1} ^ {m} w _ {k} \left(\frac {r _ {k}}{\partial r _ {k} / \partial b _ {k}}\right) ^ {2}. \tag {6.1.18}
$$

The quotient $r _ { k } / ( \partial r _ { k } / \partial b _ { k } )$ may be regarded as an inverse measure of the “impact” of the kth observation $b _ { k }$ on the model. If $\partial { r _ { k } } / \partial { b _ { k } }$ is small, then this says that the error in the model’s prediction of $b _ { k }$ is somewhat independent of $b _ { k }$ . The tendency for this to be true is lessened by basing the model on the $\lambda ^ { * }$ that minimizes $C ( \lambda )$ .

The actual determination of $\lambda ^ { * }$ is simplified by computing the SVD of A. Using the SVD (6.1.13) and Equations (6.1.17) and (6.1.18), it can be shown that

$$
C (\lambda) = \frac {1}{m} \sum_ {k = 1} ^ {m} w _ {k} \left[ \frac {\tilde {b} _ {k} - \sum_ {j = 1} ^ {r} u _ {k j} \tilde {b} _ {j} \left(\frac {\sigma_ {j} ^ {2}}{\sigma_ {j} ^ {2} + \lambda}\right)}{1 - \sum_ {j = 1} ^ {r} u _ {k j} ^ {2} \left(\frac {\sigma_ {j} ^ {2}}{\sigma_ {j} ^ {2} + \lambda}\right)} \right] ^ {2} \tag {6.1.19}
$$

where $\tilde { b } = U ^ { T } b$ . The minimization of this expression is discussed in Golub, Heath, and Wahba (1979).

# 6.1.5 Tikhonov Regularization

In the Tikhonov regularization problem, we are given $A \in \mathbb { R } ^ { m \times n } , B \in \mathbb { R } ^ { n \times n }$ , and $b \in \mathbb { R } ^ { m }$ and solve

$$
\min _ {x} \left\| \left[ \begin{array}{c} A \\ \sqrt {\lambda} B \end{array} \right] x - \left[ \begin{array}{c} b \\ 0 \end{array} \right] \right\| _ {2} ^ {2} = \min _ {x} \| A x - b \| _ {2} ^ {2} + \lambda \| B x \| _ {2} ^ {2}. \tag {6.1.20}
$$

The normal equations for this problem have the form

$$
(A ^ {T} A + \lambda B ^ {T} B) x = A ^ {T} b. \tag {6.1.21}
$$

This system is nonsingular if $\mathsf { n u l l } ( A ) \cap \mathsf { n u l l } ( B ) = \{ 0 \}$ . The matrix B can be chosen in several ways. For example, in certain data-fitting applications second derivative smoothness can be promoted by setting $B = \mathcal { T } _ { D D }$ , the second difference matrix defined in Equation 4.8.7.

To analyze how A and B interact in the Tikhonov problem, it would be handy to transform (6.1.21) into an equivalent diagonal problem. For the ridge regression problem $\left( B = I _ { n } \right)$ the SVD accomplishes this task. For the Tikhonov problem, we need a generalization of the SVD that simultaneously diagonalizes both A and B.

# 6.1.6 The Generalized Singular Value Decomposition

The generalized singular value decomposition (GSVD) set forth in Van Loan (1974) provides a useful way to simplify certain two-matrix problems such as the Tychanov regularization problem.

Theorem 6.1.1 (Generalized Singular Value Decomposition). Assume that $A \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } }$ and $B \in \mathbb { R } ^ { m _ { 2 } \times n _ { 1 } }$ with $m _ { 1 } \geq n _ { 1 }$ and

$$
r = \operatorname{rank} \left(\left[ \begin{array}{c} A \\ B \end{array} \right]\right).
$$

There exist orthogonal $U _ { 1 } \in \mathbb { R } ^ { m _ { 1 } \times m _ { 1 } }$ and $U _ { 2 } \in \mathbb { R } ^ { m _ { 2 } \times m _ { 2 } }$ and invertible $X \in \mathbb { R } ^ { n _ { 1 } \times n _ { 1 } }$ such that

$$
U _ {1} ^ {T} A X = D _ {A} = \left[ \begin{array}{c c c} I & 0 & 0 \\ 0 & \operatorname{diag} \left(\alpha_ {p + 1}, \dots , \alpha_ {r}\right) & 0 \\ 0 & 0 & 0 \\ p & r - p & n _ {1} - r \end{array} \right] \begin{array}{c} p \\ r - p \\ m _ {1} - r \end{array} , \tag {6.1.22}
$$

$$
U _ {2} ^ {T} B X = D _ {B} = \left[ \begin{array}{c c c} 0 & 0 & 0 \\ 0 & \operatorname{diag} \left(\beta_ {p + 1}, \dots , \beta_ {r}\right) & 0 \\ 0 & 0 & 0 \\ p & r - p & n _ {1} - r \end{array} \right] \begin{array}{c} p \\ r - p \\ m _ {2} - r \end{array} , \tag {6.1.23}
$$

where $p = \operatorname* { m a x } \{ r - m _ { 2 } , 0 \}$ .

Proof. The proof makes use of the SVD and the CS decomposition (Theorem 2.5.3). Let

$$
\left[ \begin{array}{l} A \\ B \end{array} \right] = \left[ \begin{array}{l l} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right] \left[ \begin{array}{c c} \Sigma_ {r} & 0 \\ 0 & 0 \end{array} \right] Z ^ {T} \tag {6.1.24}
$$

be the SVD where $\Sigma _ { r } \in \mathbb { R } ^ { r \times r }$ is nonsingular, $Q _ { 1 1 } \in \mathbb { R } ^ { m _ { 1 } \times r }$ , and $Q _ { 2 1 } \in \mathbb { R } ^ { m _ { 2 } \times r }$ . Using the CS decomposition, there exist orthogonal matrices $U _ { 1 } ( m _ { 1 } – \mathrm { b y } – m _ { 1 } ) , U _ { 2 } ( m _ { 2 } – \mathrm { b y } – m _ { 2 } )$ , and $V _ { 1 } \ ( r – \mathrm { b y } – r )$ such that

$$
\left[ \begin{array}{c c} U _ {1} & 0 \\ 0 & U _ {2} \end{array} \right] ^ {T} \left[ \begin{array}{l} Q _ {1 1} \\ Q _ {2 1} \end{array} \right] V _ {1} = \left[ \begin{array}{l} D _ {A} (:, 1: r) \\ D _ {B} (:, 1: r) \end{array} \right] \tag {6.1.25}
$$

where $D _ { A }$ and $D _ { B }$ have the forms specified by (6.1.21) and (6.1.22). It follows from (6.1.24) and (6.1.25) that

$$
\begin{array}{l} \left[ \begin{array}{c c} U _ {1} & 0 \\ 0 & U _ {2} \end{array} \right] ^ {T} \left[ \begin{array}{c} A \\ B \end{array} \right] Z = \left[ \begin{array}{c c} D _ {A} (:, 1: r) & U _ {1} Q _ {1 2} \\ D _ {B} (:, 1: r) & U _ {2} Q _ {2 2} \end{array} \right] \left[ \begin{array}{c c} V _ {1} ^ {T} \Sigma_ {r} & 0 \\ 0 & 0 \end{array} \right] \\ = \left[ \begin{array}{c c} D _ {A} (:, 1: r) & 0 \\ D _ {B} (:, 1: r) & 0 \end{array} \right] \left[ \begin{array}{c c} V _ {1} ^ {T} \Sigma_ {r} & 0 \\ 0 & I _ {n _ {1} - r} \end{array} \right] \\ = \left[ \begin{array}{c} D _ {A} \\ D _ {B} \end{array} \right] \left[ \begin{array}{c c} V _ {1} ^ {T} \Sigma_ {r} & 0 \\ 0 & I _ {n _ {1} - r} \end{array} \right]. \\ \end{array}
$$

By setting

$$
X = Z \left[ \begin{array}{c c} V _ {1} ^ {T} \Sigma_ {r} & 0 \\ 0 & I _ {n _ {1} - r} \end{array} \right] ^ {- 1}
$$

the proof is complete.

Note that if $B = I _ { n _ { 1 } }$ and we set $X = U _ { 2 }$ , then we obtain the SVD of A. The GSVD is related to the generalized eigenvalue problem

$$
A ^ {T} A x = \mu^ {2} B ^ {T} B x
$$

which is considered in §8.7.4. As with the SVD, algorithmic issues cannot be addressed until we develop procedures for the symmetric eigenvalue problem in Chapter 8.

To illustrate the insight that can be provided by the GSVD, we return to the Tikhonov regularization problem (6.1.20). If B is square and nonsingular, then the GSVD defined by (6.1.22) and (6.1.23) transforms the system (6.1.21) to

$$
(D _ {A} ^ {T} D _ {A} + \lambda D _ {B} ^ {T} D _ {B}) y = D _ {A} ^ {T} \tilde {b}
$$

where $x = X y , \tilde { b } = U _ { 1 } ^ { T } b .$ , and

$$
\left(D _ {A} ^ {T} D _ {A} + \lambda D _ {B} ^ {T} D _ {B}\right) = \operatorname{diag} \left(\alpha_ {1} ^ {2} + \lambda \beta_ {1} ^ {2}, \dots , \alpha_ {n} ^ {2} + \lambda \beta_ {n} ^ {2}\right).
$$

Thus, if

$$
X = \left[ x _ {1} \mid \dots \mid x _ {n} \right]
$$

is a column partitioning, then

$$
x (\lambda) = \sum_ {k = 1} ^ {n} \left(\frac {\alpha_ {k} \tilde {b} _ {k}}{\alpha_ {k} ^ {2} + \lambda \beta_ {k} ^ {2}}\right) x _ {k} \tag {6.1.26}
$$

solves (6.1.20). The “calming influence” of the regularization is revealed through this representation. Use of λ to manage “trouble” in the direction of $x _ { k }$ depends on the values of $\alpha _ { k }$ and $\beta _ { k }$ .

# Problems

P6.1.1 Verify (6.1.4).

P6.1.2 What is the inverse of the matrix in (6.1.5)?

P6.1.3 Show how the SVD can be used to solve the generalized LS problem (6.1.8) if the matrices A and B are rank deficient.

P6.1.4 Suppose A is the m-by-1 matrix of 1’s and let $b \in \mathbb { R } ^ { m }$ . Show that the cross-validation technique with unit weights prescribes an optimal λ given by

$$
\lambda = \left(\left(\frac {\tilde {b}}{s}\right) ^ {2} - \frac {1}{m}\right) ^ {- 1}
$$

where $\tilde { b } = ( b _ { 1 } + \cdot \cdot \cdot + b _ { m } ) / m$ and

$$
s = \sum_ {i = 1} ^ {m} (b _ {i} - \tilde {b}) ^ {2} / (m - 1).
$$

P6.1.5 Using the GSVD, give bounds for $\parallel x ( \lambda ) - x ( 0 ) \parallel$  and $\parallel A x ( \lambda ) - b \parallel _ { 2 } ^ { 2 } - \parallel A x ( 0 ) - b \parallel _ { 2 } ^ { 2 }$ where x(λ) is defined by (6.1.26).

# Notes and References for §6.1

Row and column weighting in the LS problem is discussed in Lawson and Hanson (SLS, pp. 180-88). Other analyses include:

A. van der Sluis (1969). “Condition Numbers and Equilibration of Matrices,” Numer. Math. 14, 14–23.   
G.W. Stewart (1984). “On the Asymptotic Behavior of Scaled Singular Value and QR Decompositions,” Math. Comput. 43, 483–490.   
A. Forsgren (1996). “On Linear Least-Squares Problems with Diagonally Dominant Weight Matrices,” SIAM J. Matrix Anal. Applic. 17, 763–788.   
P.D. Hough and S.A. Vavasis (1997). “Complete Orthogonal Decomposition for Weighted Least Squares,” SIAM J. Matrix Anal. Applic. 18, 551–555.   
J.K. Reid (2000). “Implicit Scaling of Linear Least Squares Problems,” BIT 40, 146–157.

For a discussion of cross-validation issues, see:

G.H. Golub, M. Heath, and G. Wahba (1979). “Generalized Cross-Validation as a Method for Choosing a Good Ridge Parameter,” Technometrics 21, 215–23.

L. Eld´en (1985). “A Note on the Computation of the Generalized Cross-Validation Function for Ill-Conditioned Least Squares Problems,” BIT 24, 467–472.

Early references concerned with the generalized singular value decomposition include:

C.F. Van Loan (1976). “Generalizing the Singular Value Decomposition,” SIAM J. Numer. Anal. 13, 76–83.

C.C. Paige and M.A. Saunders (1981). “Towards A Generalized Singular Value Decomposition,” SIAM J. Numer. Anal. 18, 398–405.   
The theoretical and computational aspects of the generalized least squares problem appear in:   
C.C. Paige (1979). “Fast Numerically Stable Computations for Generalized Linear Least Squares Problems,” SIAM J. Numer. Anal. 16, 165–171.   
C.C. Paige (1979b). “Computer Solution and Perturbation Analysis of Generalized Least Squares Problems,” Math. Comput. 33, 171–84.   
S. Kourouklis and C.C. Paige (1981). “A Constrained Least Squares Approach to the General Gauss-Markov Linear Model,” J. Amer. Stat. Assoc. 76, 620–625.   
C.C. Paige (1985). “The General Limit Model and the Generalized Singular Value Decomposition,” Lin. Alg. Applic. 70, 269–284.   
Generalized factorizations have an important bearing on generalized least squares problems, see:   
C.C. Paige (1990). “Some Aspects of Generalized QR Factorization,” in Reliable Numerical Computations, M. Cox and S. Hammarling (eds.), Clarendon Press, Oxford.   
E. Anderson, Z. Bai, and J. Dongarra (1992). “Generalized QR Factorization and Its Applications,” Lin. Alg. Applic. 162/163/164, 243–271.   
The development of regularization techniques has a long history, see:   
L. Eld´en (1977). “Algorithms for the Regularization of Ill-Conditioned Least Squares Problems,” BIT 17, 134–45.   
D.P. O’Leary and J.A. Simmons (1981). “A Bidiagonalization-Regularization Procedure for Large Scale Discretizations of Ill-Posed Problems,” SIAM J. Sci. Stat. Comput. 2, 474–489.   
L. Eld´en (1984). “An Algorithm for the Regularization of Ill-Conditioned, Banded Least Squares Problems,” SIAM J. Sci. Stat. Comput. 5, 237–254.   
P.C. Hansen (1990). “Relations Between SVD and GSVD of Discrete Regularization Problems in Standard and General Form,” Lin.Alg. Applic. 141, 165–176.   
P.C. Hansen (1995). “Test Matrices for Regularization Methods,” SIAM J. Sci. Comput. 16, 506–512.   
A. Neumaier (1998). “Solving Ill–Conditioned and Singular Linear Systems: A Tutorial on Regularization,” SIAM Review 40, 636–666.   
P.C. Hansen (1998). Rank-Deficient and Discrete Ill-Posed Problems: Numerical Aspects of Linear Inversion, SIAM Publications, Philadelphia, PA.   
M.E. Gulliksson and P.-A. Wedin (2000). “The Use and Properties of Tikhonov Filter Matrices,” SIAM J. Matrix Anal. Applic. 22, 276–281.   
M.E. Gulliksson, P.-A. Wedin, and Y. Wei (2000). “Perturbation Identities for Regularized Tikhonov Inverses and Weighted Pseudoinverses,” BIT 40, 513–523.   
T. Kitagawa, S. Nakata, and Y. Hosoda (2001). “Regularization Using QR Factorization and the Estimation of the Optimal Parameter,” BIT 41, 1049–1058.   
M.E. Kilmer and D.P. O’Leary. (2001). “Choosing Regularization Parameters in Iterative Methods for Ill-Posed Problems,” SIAM J. Matrix Anal. Applic. 22, 1204–1221.   
A. N. Malyshev (2003). “A Unified Theory of Conditioning for Linear Least Squares and Tikhonov Regularization Solutions,” SIAM J. Matrix Anal. Applic. 24, 1186–1196.   
M. Hanke (2006). “A Note on Tikhonov Regularization of Large Linear Problems,” BIT 43, 449–451.   
P.C. Hansen, J.G. Nagy, and D.P. OLeary (2006). Deblurring Images: Matrices, Spectra, and Filtering, SIAM Publications, Philadelphia, PA.   
M.E. Kilmer, P.C. Hansen, and M.I. Espa˜nol (2007). “A Projection-Based Approach to General-Form Tikhonov Regularization,” SIAM J. Sci. Comput. 29, 315–330.   
T. Elfving and I. Skoglund (2009). “A Direct Method for a Regularized Least-Squares Problem,” Num. Lin. Alg. Applic. 16, 649–675.   
I. Hn˘etynkov´a and M. Ple˘singer (2009). “The Regularizing Effect of the Golub-Kahan Iterative Bidiagonalization and revealing the Noise level in Data,” BIT 49, 669–696.   
P.C. Hansen (2010). Discrete Inverse Problems: Insight and Algorithms, SIAM Publications, Philadelphia, PA.
