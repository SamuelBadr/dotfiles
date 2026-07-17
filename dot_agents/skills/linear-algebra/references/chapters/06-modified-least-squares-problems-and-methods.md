# Chapter 6

# Modified Least Squares Problems and Methods

6.1 Weighting and Regularization   
6.2 Constrained Least Squares   
6.3 Total Least Squares   
6.4 Subspace Computations with the SVD   
6.5 Updating Matrix Factorizations

In this chapter we discuss an assortment of least square problems that can be solved using QR and SVD. We also introduce a generalization of the SVD that can be used to simultaneously diagonalize a pair of matrices, a maneuver that is useful in certain applications.

The first three sections deal with variations of the ordinary least squares problem that we treated in Chapter 5. The unconstrained minimization of $\parallel A x - b \parallel _ { 2 }$ does not always make a great deal of sense. How do we balance the importance of each equation in $\boldsymbol { A } \boldsymbol { x } = \boldsymbol { b } \boldsymbol { ? }$ How might we control the size of x if A is ill-conditioned? How might we minimize $\parallel A x - b \parallel _ { 2 }$ over a proper subspace of $\mathbb { R } ^ { n } ?$ What if there are errors in the “data matrix” A in addition to the usual errors in the “vector of observations” b?

In §6.4 we consider a number of multidimensional subspace computations including the problem of determining the principal angles between a pair of given subspaces. The SVD plays a prominent role.

The final section is concerned with the updating of matrix factorizations. In many applications, one is confronted with a succession of least squares (or linear equation) problems where the matrix associated with the current step is highly related to the matrix associated with the previous step. This opens the door to updating strategies that can reduce factorization overheads by an order of magnitude.

# Reading Notes

Knowledge of Chapter 5 is assumed. The sections in this chapter are independent of each other except that §6.1 should be read before §6.2. Excellent global references include Bj¨orck (NMLS) and Lawson and Hansen (SLS).

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

# 6.2 Constrained Least Squares

In the least squares setting it is sometimes natural to minimize $\parallel A x - b \parallel _ { 2 }$ over a proper subset of $\mathbb { R } ^ { n }$ . For example, we may wish to predict b as best we can with Ax subject to the constraint that x is a unit vector. Or perhaps the solution defines a fitting function f (t) which is to have prescribed values at certain points. This can lead to an equality-constrained least squares problem. In this section we show how these problems can be solved using the QR factorization, the SVD, and the GSVD.

# 6.2.1 Least Squares Minimization Over a Sphere

Given $\boldsymbol { A } \in \mathbb { R } ^ { m \times n } , \boldsymbol { b } \in \mathbb { R } ^ { m }$ , and a positive $\alpha \in \mathbb { R }$ , we consider the problem

$$
\min _ {\| x \| _ {2} \leq \alpha} \| A x - b \| _ {2}. \tag {6.2.1}
$$

This is an example of the LSQI (least squares with quadratic inequality constraint) problem. This problem arises in nonlinear optimization and other application areas. As we are soon to observe, the LSQI problem is related to the ridge regression problem discussed in §6.1.4.

Suppose

$$
A = U \Sigma V ^ {T} = \sum_ {i = 1} ^ {r} \sigma_ {i} u _ {i} v _ {i} ^ {T} \tag {6.2.2}
$$

is the SVD of A which we assume to have rank r. If the unconstrained minimum norm solution

$$
x _ {L S} = \sum_ {i = 1} ^ {r} \frac {u _ {i} ^ {T} b}{\sigma_ {i}} v _ {i}
$$

satisfies $\Vert \textbf { } x _ { L S } \textbf { } \Vert _ { 2 } \leq \alpha$ , then it obviously solves (6.2.1). Otherwise,

$$
\left\| x _ {L S} \right\| _ {2} ^ {2} = \sum_ {i = 1} ^ {r} \left(\frac {u _ {i} ^ {T} b}{\sigma_ {i}}\right) ^ {2} > \alpha^ {2}, \tag {6.2.3}
$$

and it follows that the solution to (6.2.1) is on the boundary of the constraint sphere. Thus, we can approach this constrained optimization problem using the method of Lagrange multipliers. Define the parameterized objective function $\phi$ by

$$
\phi (x, \lambda) = \frac {1}{2} \| A x - b \| _ {2} ^ {2} + \frac {\lambda}{2} \left(\| x \| _ {2} ^ {2} - \alpha^ {2}\right)
$$

and equate its gradient to zero. This gives a shifted normal equation system:

$$
(A ^ {T} A + \lambda I) \cdot x (\lambda) = A ^ {T} b.
$$

The goal is to choose λ so that $\parallel x ( \lambda ) \parallel _ { 2 } = \alpha$ . Using the SVD (6.2.2), this leads to the problem of finding a zero of the function

$$
f (\lambda) = \| x (\lambda) \| _ {2} ^ {2} - \alpha^ {2} = \sum_ {k = 1} ^ {n} \left(\frac {\sigma_ {k} u _ {k} ^ {T} b}{\sigma_ {k} ^ {2} + \lambda}\right) ^ {2} - \alpha^ {2}.
$$

This is an example of a secular equation problem. From (6.2.3), $f ( 0 ) > 0$ . Since $f ^ { \prime } ( \lambda ) < 0$ for $\lambda \geq 0$ , it follows that f has a unique positive root $\lambda _ { + }$ . It can be shown that

$$
\rho (\lambda) = \| A x (\lambda) - b \| _ {2} ^ {2} = \| A x _ {L S} - b \| _ {2} ^ {2} + \sum_ {i = 1} ^ {r} \left(\frac {\lambda u _ {i} ^ {T} b}{\sigma_ {i} ^ {2} + \lambda}\right) ^ {2}. \tag {6.2.4}
$$

It follows that $x ( \lambda _ { + } )$ solves (6.2.1).

Algorithm 6.2.1 Given $A \in \mathbb { R } ^ { m \times n }$ with $m \ge n , b \in \mathbb { R } ^ { m }$ , and $\alpha > 0$ , the following algorithm computes a vector $\boldsymbol { x } \in \mathbb { R } ^ { n }$ such that $\parallel A x - b \parallel _ { 2 }$ is minimum subject to the constraint that $\Vert \ b { x } \Vert _ { 2 } \leq \alpha$ .

Compute the SVD $A = U \Sigma V ^ { T }$ , save $V = [  v _ { 1 } | \cdot \cdot \cdot | v _ { n } ]$ , form $\tilde { b } = U ^ { T } b$ , and determine $r = \mathsf { r a n k } ( A )$ .

$\mathrm { i f } \sum _ { i = 1 } ^ { r } \left( \frac { { \tilde { b } } _ { i } } { \sigma _ { i } } \right) ^ { 2 } > \alpha ^ { 2 }$ r > α2

Find λ+ > 0 such that $\sum _ { i = 1 } ^ { r } \left( \frac { \sigma _ { i } \tilde { b } _ { i } } { \sigma _ { i } ^ { 2 } + \lambda _ { + } } \right) ^ { 2 } = \alpha ^ { 2 }$ .

$$
x = \sum_ {i = 1} ^ {r} \left(\frac {\sigma_ {i} \tilde {b} _ {i}}{\sigma_ {i} ^ {2} + \lambda_ {+}}\right) v _ {i}
$$

else

$$
x = \sum_ {i = 1} ^ {r} \left(\frac {\tilde {b} _ {i}}{\sigma_ {i}}\right) v _ {i}
$$

end

The SVD is the dominant computation in this algorithm.

# 6.2.2 More General Quadratic Constraints

A more general version of (6.2.1) results if we minimize $\parallel A x - b \parallel _ { 2 }$ over an arbitrary hyperellipsoid:

$$
\text { minimize } \parallel A x - b \parallel_ {2} \quad \text { subject   to } \parallel B x - d \parallel_ {2} \leq \alpha . \tag {6.2.5}
$$

Here we are assuming that $A \in \mathbb R ^ { m _ { 1 } \times n _ { 1 } } , b \in \mathbb R ^ { m _ { 1 } } , B \in \mathbb R ^ { m _ { 2 } \times n _ { 1 } } , d \in \mathbb R ^ { m _ { 2 } }$ , and $\alpha \geq 0$ . Just as the SVD turns (6.2.1) into an equivalent diagonal problem, we can use the GSVD to transform (6.2.5) into a diagonal problem. In particular, if the GSVD of A and B is given by (6.1.22) and (6.2.23), then (6.2.5) is equivalent to

$$
\text { minimize } \parallel D _ {A} y - \tilde {b} \parallel_ {2} \quad \text { subject   to } \parallel D _ {B} y - \tilde {d} \parallel_ {2} \leq \alpha \tag {6.2.6}
$$

where

$$
\tilde {b} = U _ {1} ^ {T} b, \qquad \tilde {d} = U _ {2} ^ {T} d, \qquad y = X ^ {- 1} x.
$$

The simple form of the objective function and the constraint equation facilitate the analysis. For example, if rank $( B ) = m _ { 2 } < n _ { 1 }$ , then

$$
\| D _ {A} y - \tilde {b} \| _ {2} ^ {2} = \sum_ {i = 1} ^ {n _ {1}} \left(\alpha_ {i} y _ {i} - \tilde {b} _ {i}\right) ^ {2} + \sum_ {i = n _ {1} + 1} ^ {m _ {1}} \tilde {b} _ {i} ^ {2} \tag {6.2.7}
$$

and

$$
\| D _ {B} y - \tilde {d} \| _ {2} ^ {2} = \sum_ {i = 1} ^ {m _ {2}} \left(\beta_ {i} y _ {i} - \tilde {d} _ {i}\right) ^ {2} + \sum_ {i = m _ {2} + 1} ^ {n _ {1}} \tilde {d} _ {i} ^ {2} \leq \alpha^ {2}. \tag {6.2.8}
$$

A Lagrange multiplier argument can be used to determine the solution to this transformed problem (if it exists).

# 6.2.3 Least Squares With Equality Constraints

We consider next the constrained least squares problem

$$
\min _ {B x = d} \| A x - b \| _ {2} \tag {6.2.9}
$$

where $A \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } }$ with $m _ { 1 } \geq n _ { 1 } , B \in \mathbb { R } ^ { m _ { 2 } \times n _ { 1 } }$ with $m _ { 2 } < n _ { 1 } , b \in \mathbb { R } ^ { m _ { 1 } }$ , and $d \in \mathbb { R } ^ { m _ { 2 } }$ . We refer to this as the LSE problem (least squares with equality constraints). By setting $\alpha = 0$ in (6.2.5) we see that the LSE problem is a special case of the LSQI problem. However, it is simpler to approach the LSE problem directly rather than through Lagrange multipliers.

For clarity, we assume that both A and B have full rank. Let

$$
Q ^ {T} B ^ {T} = \left[ \begin{array}{c} R \\ 0 \end{array} \right] _ {n _ {1} - m _ {2}} ^ {n _ {1}}
$$

be the QR factorization of $B ^ { T }$ and set

$$
A Q = \left[ \begin{array}{c c} A _ {1} & A _ {2} \\ m _ {2} & n _ {1} - m _ {2} \end{array} \right], \qquad Q ^ {T} x = \left[ \begin{array}{c} y \\ z \end{array} \right] \begin{array}{c} m _ {2} \\ n _ {1} - m _ {2} \end{array} .
$$

It is clear that with these transformations (6.2.9) becomes

$$
\min _ {R ^ {T} y = d} \| A _ {1} y + A _ {2} z - b \| _ {2}.
$$

Thus, $y$ is determined from the constraint equation $R ^ { T } y ~ = ~ d$ and the vector z is obtained by solving the unconstrained LS problem

$$
\min _ {z \in \mathbb {R} ^ {n _ {1} - m _ {2}}} \| A _ {2} z - (b - A _ {1} y) \| _ {2}.
$$

Combining the above, we see that the following vector solves the LSE problem:

$$
x = Q \left[ \begin{array}{c} y \\ z \end{array} \right].
$$

Algorithm 6.2.2 Suppose $A \in \mathbb { R } ^ { m _ { 1 } \times n _ { 1 } } , \ B \in \mathbb { R } ^ { m _ { 2 } \times n _ { 1 } } , \ b \in \mathbb { R } ^ { m _ { 1 } }$ , and $d \in \mathbb { R } ^ { m _ { 2 } }$ . If rank $( A ) \ = \ n _ { 1 }$ and rank $( B ) \ = \ m _ { 2 } \ < \ n _ { 1 }$ , then the following algorithm minimizes $\parallel A x - b \parallel _ { 2 }$ subject to the constraint $B x = d$ .

Compute the QR factorization $B ^ { T } = Q R$

Solve $R ( 1 { : } m _ { 2 } , 1 { : } m _ { 2 } ) ^ { T } { \cdot } y = d$ for y.

$$
A = A Q
$$

${ \mathrm { F i n d ~ } } z { \mathrm { ~ s o ~ } } \parallel A ( : , m _ { 2 } + 1 : n _ { 1 } ) z - ( b - A ( : , 1 { : } m _ { 2 } ) \cdot y ) \parallel _ { 2 } { \mathrm { ~ i s ~ m i n i m i z e d } } .$

$$
x = Q (:, 1: m _ {2}) \cdot y + Q (:, m _ {2} + 1: n _ {1}) \cdot z.
$$

Note that this approach to the LSE problem involves two QR factorizations and a matrix multiplication. If A and/or B are rank deficient, then it is possible to devise a similar solution procedure using the SVD instead of QR. Note that there may not be a solution if rank $( B ) < m _ { 2 }$ . Also, if nul $( A ) \cap { \mathsf { n u l l } } ( B ) \neq \{ 0 \}$ and $d \in \mathsf { r a n } ( B )$ , then the LSE solution is not unique.

# 6.2.4 LSE Solution Using the Augmented System

The LSE problem can also be approached through the method of Lagrange multipliers. Define the augmented objective function

$$
f (x, \lambda) = \frac {1}{2} \| A x - b \| _ {2} ^ {2} + \lambda^ {T} (d - B x), \quad \lambda \in \mathbb {R} ^ {m _ {2}},
$$

and set to zero its gradient with respect to x:

$$
A ^ {T} A x - A ^ {T} b - B ^ {T} \lambda = 0.
$$

Combining this with the equations $r = b - A x$ and $B x = d$ we obtain the symmetric indefinite linear system

$$
\left[ \begin{array}{c c c} 0 & A ^ {T} & B ^ {T} \\ A & I & 0 \\ B & 0 & 0 \end{array} \right] \left[ \begin{array}{l} x \\ r \\ \lambda \end{array} \right] = \left[ \begin{array}{l} 0 \\ b \\ d \end{array} \right]. \tag {6.2.10}
$$

This system is nonsingular if both A and B have full rank. The augmented system presents a solution framework for the sparse LSE problem.

# 6.2.5 LSE Solution Using the GSVD

Using the GSVD given by (6.1.22) and (6.1.23), we see that the LSE problem transforms to

$$
\min _ {D _ {B} y = \tilde {d}} \| D _ {A} y - \tilde {b} \| _ {2} \tag {6.2.11}
$$

where $\tilde { b } = U _ { 1 } ^ { T } b , \tilde { d } = U _ { 2 } ^ { T } d .$ and $y = X ^ { - 1 } x$ . It follows that if $\mathsf { n u l l } ( A ) \cap \mathsf { n u l l } ( B ) = \{ 0 \}$ and $X = [ x _ { 1 } \mid \cdots \mid x _ { n } ]$ , then

$$
x = \sum_ {i = 1} ^ {m _ {2}} \left(\frac {\tilde {d} _ {i}}{\beta_ {i}}\right) x _ {i} + \sum_ {i = m _ {2} + 1} ^ {n _ {1}} \left(\frac {\tilde {b} _ {i}}{\alpha_ {i}}\right) x _ {i} \tag {6.2.12}
$$

solves the LSE problem.

# 6.2.6 LSE Solution Using Weights

An interesting way to obtain an approximate LSE solution is to solve the unconstrained LS problem

$$
\min _ {x} \left\| \left[ \begin{array}{c} A \\ \sqrt {\lambda} B \end{array} \right] x - \left[ \begin{array}{c} b \\ \sqrt {\lambda} d \end{array} \right] \right\| _ {2} \tag {6.2.13}
$$

for large λ. (Compare with the Tychanov regularization problem (6.1.21).) Since

$$
\left\| \left[ \begin{array}{c} A \\ \sqrt {\lambda} B \end{array} \right] x \right. - \left[ \begin{array}{c} b \\ \sqrt {\lambda} d \end{array} \right] \left\| \right. _ {2} ^ {2} = \left\| A x - b \right\| _ {2} ^ {2} + \lambda \left\| B x - d \right\| ^ {2},
$$

we see that there is a penalty for discrepancies among the constraint equations. To quantify this, assume that both A and B have full rank and substitute the GSVD defined by (6.1.22) and (6.1.23) into the normal equation system

$$
(A ^ {T} A + \lambda B ^ {T} B) x = A ^ {T} b + \lambda B ^ {T} d.
$$

This shows that the solution $x ( \lambda )$ is given by $x ( \lambda ) = X y ( \lambda )$ where $y ( \lambda )$ solves

$$
(D _ {A} ^ {T} D _ {A} + \lambda D _ {B} ^ {T} D _ {B}) y = D _ {A} ^ {T} \tilde {b} + \lambda D _ {B} ^ {T} \tilde {d}
$$

with $\tilde { b } = U _ { 1 } ^ { T } b$ and $\tilde { d } = U _ { 2 } ^ { T } d .$ . It follows that

$$
x (\lambda) = \sum_ {i = 1} ^ {m _ {2}} \left(\frac {\alpha_ {i} \tilde {b} _ {i} + \lambda \beta_ {i} \tilde {d} _ {i}}{\alpha_ {i} ^ {2} + \lambda \beta_ {i} ^ {2}}\right) x _ {i} + \sum_ {i = m _ {2} + 1} ^ {n _ {1}} \left(\frac {\tilde {b} _ {i}}{\alpha_ {i}}\right) x _ {i}
$$

and so from (6.2.13) we have

$$
x (\lambda) - x = \sum_ {i = 1} ^ {p} \frac {\alpha_ {i}}{\beta_ {i}} \left(\frac {\beta_ {i} u _ {i} ^ {T} b - \alpha_ {i} v _ {i} ^ {T} d}{\alpha_ {i} ^ {2} + \lambda^ {2} \beta_ {i} ^ {2}}\right) x _ {i}. \tag {6.2.14}
$$

This shows that $x ( \lambda ) \to x { \mathrm { ~ a s ~ } } \lambda \to \infty$ . The appeal of this approach to the LSE problem is that it can be implemented with unconstrained LS problem software. However, for large values of λ numerical problems can arise and it is necessary to take precautions. See Powell and Reid (1968) and Van Loan (1982).

# Problems

P6.2.1 Is the solution to (6.2.1) always unique?

P6.2.2 Let $p _ { 0 } ( x ) , \ldots , p _ { n } ( x )$ be given polynomials and $( x _ { 0 } , y _ { 0 } ) , \dots , ( x _ { m } , y _ { m } )$ be a given set of coordinate pairs with $x _ { i } \in [ a , b ]$ . It is desired to find a polynomial $\begin{array} { r } { p ( x ) \ = \ \sum _ { k = 0 } ^ { n } \alpha _ { k } p _ { k } ( x ) } \end{array}$ such that

$$
\phi (\alpha) = \sum_ {i = 0} ^ {m} (p (x _ {i}) - y _ {i}) ^ {2}
$$

is minimized subject to the constraint that

$$
\int_ {a} ^ {b} [ p ^ {\prime \prime} (x) ] ^ {2} d x \approx h \sum_ {i = 0} ^ {N} \left(\frac {p (z _ {i - 1}) - 2 p (z _ {i}) + p (z _ {i + 1})}{h ^ {2}}\right) ^ {2} \leq \alpha^ {2}
$$

where $z _ { i } = a + i h$ and $b = a + N h$ . Show that this leads to an LSQI problem of the form (6.2.5) with $d = 0$ .

P6.2.3 Suppose $Y = [ y _ { 1 } \mid \cdot \cdot \cdot \mid y _ { k } ] \in \mathbb { R } ^ { m \times k }$ has the property that

$$
Y ^ {T} Y = \mathrm{diag} (d _ {1} ^ {2}, \ldots , d _ {k} ^ {2}), \qquad d _ {1} \geq d _ {2} \geq \dots \geq d _ {k} > 0.
$$

Show that if $Y = Q R$ is the QR factorization of Y , then R is diagonal with $| r _ { i i } | = d _ { i }$ .

P6.2.4 (a) Show that if $( A ^ { T } A + \lambda I ) x = A ^ { T } b , \lambda > 0 .$ , and $\parallel x \parallel _ { 2 } = \alpha ,$ then $z = ( A x - b ) / \lambda$ solves the dual equations $( A A ^ { T } + \lambda I ) z = - b$ with $\parallel A ^ { T } z \parallel _ { 2 } = \alpha . ( \mathrm { b } )$ Show that if $( A A ^ { T } + \lambda I ) z = - b ,$ $\parallel A ^ { T } z \parallel _ { 2 } = \alpha$ , then $x = - A ^ { T } z$ satisfies $( A ^ { T } A + \lambda I ) x = { \overset { - } { A } } A ^ { T } b , \| x \| _ { 2 } = \alpha .$ .

P6.2.5 Show how to compute y (if it exists) so that both (6.2.7) and (6.2.8) are satisfied.

P6.2.6 Develop an SVD version of Algorithm 6.2.2 that can handle the situation when A and/or B are rank deficient.

P6.2.7 Suppose

$$
A = \left[ \begin{array}{l} A _ {1} \\ A _ {2} \end{array} \right]
$$

where $A _ { 1 } \in \mathbb { R } ^ { n \times n }$ is nonsingular and $A _ { 2 } \in \mathbb { R } ^ { ( m - n ) \times n }$ . Show that

$$
\sigma_ {\min} (A) \geq \sqrt {1 + \sigma_ {\min} (A _ {2} A _ {1} ^ {- 1}) ^ {2}} \sigma_ {\min} (A _ {1}).
$$

P6.2.8 Suppose $p \geq m \geq n$ and that $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { m \times p }$ Show how to compute orthogonal $Q \in \mathbb { R } ^ { m \times m }$ and orthogonal $V \in \mathbb { R } ^ { n \times n }$ so that

$$
Q ^ {T} A = \left[ \begin{array}{c} R \\ 0 \end{array} \right], \qquad Q ^ {T} B V = \left[ 0 \mid S \right]
$$

where $R \in \mathbb { R } ^ { n \times n }$ and $S \in \mathbb { R } ^ { m \times m }$ are upper triangular.

P6.2.9 Suppose $r \in \mathbb { R } ^ { m } , y \in \mathbb { R } ^ { n }$ , and $\delta > 0$ . Show how to solve the problem

$$
\min _ {E \in \mathbf {R} ^ {m \times n}, \| E \| _ {F} \leq \delta} \| E y - r \| _ {2}
$$

Repeat with “min” replaced by “max.”

P6.2.10 Show how the constrained least squares problem

$$
\min _ {B x = d} \| A x - b \| _ {2} \quad A \in \mathbb {R} ^ {m \times n}, B \in \mathbb {R} ^ {p \times n}, \operatorname{rank} (B) = p
$$

can be reduced to an unconstrained least square problem by performing p steps of Gaussian elimination on the matrix

$$
\left[ \begin{array}{l} B \\ A \end{array} \right] = \left[ \begin{array}{l l} B _ {1} & B _ {2} \\ A _ {1} & A _ {2} \end{array} \right], \qquad B _ {1} \in \mathbb {R} ^ {p \times p},   \mathsf {r a n k} (B _ {1}) = p.
$$

Explain. Hint: The Schur complement is of interest.

# Notes and References for §6.2

The LSQI problem is discussed in:

G.E. Forsythe and G.H. Golub (1965). “On the Stationary Values of a Second-Degree Polynomial on the Unit Sphere,” SIAM J. App. Math. 14, 1050–1068.

L. Eld´en (1980). “Perturbation Theory for the Least Squares Problem with Linear Equality Constraints,” SIAM J. Numer. Anal. 17, 338–350.

W. Gander (1981). “Least Squares with a Quadratic Constraint,” Numer. Math. 36, 291–307.

L. Eld´en (1983). “A Weighted Pseudoinverse, Generalized Singular Values, and Constrained Least Squares Problems,” BIT 22 , 487–502.

G.W. Stewart (1984). “On the Asymptotic Behavior of Scaled Singular Value and QR Decompositions,” Math. Comput. 43, 483–490.   
G.H. Golub and U. von Matt (1991). “Quadratically Constrained Least Squares and Quadratic Problems,” Numer. Math. 59, 561–580.   
T.F. Chan, J.A. Olkin, and D. Cooley (1992). “Solving Quadratically Constrained Least Squares Using Black Box Solvers,” BIT 32, 481–495.   
Secular equation root-finding comes up in many numerical linear algebra settings. For an algorithmic overview, see:   
O.E. Livne and A. Brandt (2002). “N Roots of the Secular Equation in O(N) Operations,” SIAM J. Matrix Anal. Applic. 24, 439–453.

For a discussion of the augmented systems approach to least squares problems, see:

˚A. Bj¨orck (1992). “Pivoting and Stability in the Augmented System Method,” Proceedings of the 14th Dundee Conference, D.F. Griffiths and G.A. Watson (eds.), Longman Scientific and Technical, Essex, U.K.   
˚A. Bj¨orck and C.C. Paige (1994). “Solution of Augmented Linear Systems Using Orthogonal Factorizations,” BIT 34, 1–24.

References that are concerned with the method of weighting for the LSE problem include:

M.J.D. Powell and J.K. Reid (1968). “On Applying Householder’s Method to Linear Least Squares Problems,” Proc. IFIP Congress, pp. 122–26.

C. Van Loan (1985). “On the Method of Weighting for Equality Constrained Least Squares Problems,” SIAM J. Numer. Anal. 22, 851–864.

J.L. Barlow and S.L. Handy (1988). “The Direct Solution of Weighted and Equality Constrained Least-Squares Problems,” SIAM J. Sci. Stat. Comput. 9, 704–716.

J.L. Barlow, N.K. Nichols, and R.J. Plemmons (1988). “Iterative Methods for Equality Constrained Least Squares Problems,” SIAM J. Sci. Stat. Comput. 9, 892–906.

J.L. Barlow (1988). “Error Analysis and Implementation Aspects of Deferred Correction for Equality Constrained Least-Squares Problems,” SIAM J. Numer. Anal. 25, 1340–1358.

J.L. Barlow and U.B. Vemulapati (1992). “A Note on Deferred Correction for Equality Constrained Least Squares Problems,” SIAM J. Numer. Anal. 29, 249–256.

M. Gulliksson and P.-˚A. Wedin (1992). “Modifying the QR-Decomposition to Constrained and Weighted Linear Least Squares,” SIAM J. Matrix Anal. Applic. 13, 1298–1313.

M. Gulliksson (1994). “Iterative Refinement for Constrained and Weighted Linear Least Squares,” BIT 34, 239–253.

G. W. Stewart (1997). “On the Weighting Method for Least Squares Problems with Linear Equality Constraints,” BIT 37, 961–967.

For the analysis of the LSE problem and related methods, see:

M. Wei (1992). “Perturbation Theory for the Rank-Deficient Equality Constrained Least Squares Problem,” SIAM J. Numer. Anal. 29, 1462–1481.

M. Wei (1992). “Algebraic Properties of the Rank-Deficient Equality-Constrained and Weighted Least Squares Problems,” Lin. Alg. Applic. 161, 27–44.

M. Gulliksson (1995). “Backward Error Analysis for the Constrained and Weighted Linear Least Squares Problem When Using the Weighted QR Factorization,” SIAM J. Matrix. Anal. Applic. 13, 675–687.

M. Gulliksson (1995). “Backward Error Analysis for the Constrained and Weighted Linear Least Squares Problem When Using the Weighted QR Factorization,” SIAM J. Matrix Anal. Applic. 16, 675–687.

J. Ding and W. Hang (1998). “New Perturbation Results for Equality-Constrained Least Squares Problems,” Lin. Alg. Applic. 272, 181–192.

A.J. Cox and N.J. Higham (1999). “Accuracy and Stability of the Null Space Method for Solving the Equality Constrained Least Squares Problem,” BIT 39, 34–50.

A.J. Cox and N.J. Higham (1999). “Row-Wise Backward Stable Elimination Methods for the Equality Constrained Least Squares Problem,” SIAM J. Matrix Anal. Applic. 21, 313–326.

A.J. Cox and Nicholas J. Higham (1999). “Backward Error Bounds for Constrained Least Squares Problems,” BIT 39, 210–227.

M. Gulliksson and P-A. Wedin (2000). “Perturbation Theory for Generalized and Constrained Linear Least Squares,” Num. Lin. Alg. 7, 181–195.   
M. Wei and A.R. De Pierro (2000). “Upper Perturbation Bounds of Weighted Projections, Weighted and Constrained Least Squares Problems,” SIAM J. Matrix Anal. Applic. 21, 931–951.   
E.Y. Bobrovnikova and S.A. Vavasis (2001). “Accurate Solution of Weighted Least Squares by Iterative Methods SIAM. J. Matrix Anal. Applic. 22, 1153–1174.   
M. Gulliksson, X-Q.Jin, and Y-M. Wei (2002). “Perturbation Bounds for Constrained and Weighted Least Squares Problems,” Lin. Alg. Applic. 349, 221–232.

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

# 6.4 Subspace Computations with the SVD

It is sometimes necessary to investigate the relationship between two given subspaces. How close are they? Do they intersect? Can one be “rotated” into the other? And so on. In this section we show how questions like these can be answered using the singular value decomposition.

# 6.4.1 Rotation of Subspaces

Suppose $A \in \mathbb { R } ^ { m \times p }$ is a data matrix obtained by performing a certain set of experiments. If the same set of experiments is performed again, then a different data matrix, $B \in \mathbb { R } ^ { m \times p }$ , is obtained. In the orthogonal Procrustes problem the possibility that B can be rotated into A is explored by solving the following problem:

$$
\text { minimize } \parallel A - B Q \parallel_ {F}, \quad \text { subject   to } Q ^ {T} Q = I _ {p}. \tag {6.4.1}
$$

We show that optimizing $Q$ can be specified in terms of the SVD of $B ^ { T } A$ . The matrix trace is critical to the derivation. The trace of a matrix is the sum of its diagonal entries:

$$
\operatorname{tr} (C) = \sum_ {i = 1} ^ {n} c _ {i i}, \quad C \in \mathbb {R} ^ {n \times n}.
$$

It is easy to show that if $C _ { 1 }$ and $C _ { 2 }$ have the same row and column dimension, then

$$
\operatorname{tr} (C _ {1} ^ {T} C _ {2}) = \operatorname{tr} (C _ {2} ^ {T} C _ {1})  . \tag {6.4.2}
$$

Returning to the Procrustes problem (6.4.1), if $Q \in \mathbb { R } ^ { p \times p }$ is orthogonal, then

$$
\begin{array}{l} \left\| A - B Q \right\| _ {F} ^ {2} = \sum_ {k = 1} ^ {p} \left\| A (:, k) - B \cdot Q (:, k) \right\| _ {2} ^ {2} \\ = \sum_ {k = 1} ^ {p} \| A (:, k) \| _ {2} ^ {2} + \| B Q (:, k) \| _ {2} ^ {2} - 2 Q (:, k) ^ {T} B ^ {T} A (:, k) \\ = \| A \| _ {F} ^ {2} + \| B Q \| _ {F} ^ {2} - 2 \sum_ {k = 1} ^ {p} \left[ Q ^ {T} \left(B ^ {T} A\right) \right] _ {k k} \\ = \| A \| _ {F} ^ {2} + \| B \| _ {F} ^ {2} - 2 \operatorname{tr} \left(Q ^ {T} \left(B ^ {T} A\right)\right). \\ \end{array}
$$

Thus, (6.4.1) is equivalent to the problem

$$
\max _ {Q ^ {T} Q = I _ {p}} \operatorname{tr} (Q ^ {T} B ^ {T} A)  .
$$

If $\begin{array} { c c c c c } { { U ^ { T } ( B ^ { T } A ) V } } & { { = } } & { { \Sigma } } & { { = } } & { { \mathrm { d i a g } ( \sigma _ { 1 } , \dots , \sigma _ { p } ) } } \end{array}$ is the SVD of $B ^ { T } A$ and we define the orthogonal matrix Z by $Z ~ = ~ V ^ { T } Q ^ { T } U$ , then by using (6.4.2) we have

$$
\mathsf {t r} (Q ^ {T} B ^ {T} A) = \mathsf {t r} (Q ^ {T} U \Sigma V ^ {T}) = \mathsf {t r} (Z \Sigma) = \sum_ {i = 1} ^ {p} z _ {i i} \sigma_ {i} \leq \sum_ {i = 1} ^ {p} \sigma_ {i}.
$$

The upper bound is clearly attained by setting $Z = I _ { p } , \mathrm { i . e . , } Q = U V ^ { T }$ .

Algorithm 6.4.1 Given A and B in $\mathbb { R } ^ { m \times p }$ , the following algorithm finds an orthogonal $Q \in \mathbb { R } ^ { p \times p }$ such that $\| A - B Q \| _ { F }$ is minimum.

$$
C = B ^ {T} A
$$

Compute the SVD $U ^ { T } C V = \Sigma$ and save U and V .

$$
Q = U V ^ {T}
$$

We mention that if $B = I _ { p } ,$ , then the problem (6.4.1) is related to the polar decomposition. This decomposition states that any square matrix A has a factorization of the form $A = Q P$ where Q is orthogonal and P is symmetric and positive semidefinite. Note that if $A = U \Sigma V ^ { T }$ is the SVD of A, then $\overset { \cdot } { A } = ( U V ^ { T } ) ( V \overset { \cdot } { \Sigma } V ^ { T } )$ is its polar decomposition. For further discussion, see §9.4.3.

# 6.4.2 Intersection of Nullspaces

Let $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { p \times n }$ be given, and consider the problem of finding an orthonormal basis for null $( A ) \cap \mathsf { n u l l } ( B )$ . One approach is to compute the nullspace of the matrix

$$
C = \left[ \begin{array}{l} A \\ B \end{array} \right]
$$

since this is just what we want: $C x = 0 \Leftrightarrow x \in \mathsf { n u l l } ( A ) \cap \mathsf { n u l l } ( B )$ . However, a more economical procedure results if we exploit the following theorem.

Theorem 6.4.1. Suppose $A \in \mathbb { R } ^ { m \times n }$ and let $\{ z _ { 1 } , \ldots , z _ { t } \}$ be an orthonormal basis for null(A). Define $Z ~ = ~ \left[ ~ z _ { 1 } ~ | \cdots | ~ z _ { t } ~ \right]$ and let $\{ w _ { 1 } , \ldots , w _ { q } \}$ be an orthonormal basis for null(BZ) where $B \in \mathbb { R } ^ { p \times n }$ . $I f W = [  w _ { 1 } | \cdot \cdot \cdot | w _ { q } ]$ , then the columns of ZW form an orthonormal basis for null(A) ∩ null(B).

Proof. Since $A Z = 0$ and $( B Z ) W = 0$ , we clearly have ran $( Z W ) \subset$ null(A) ∩ null(B). Now suppose x is in both null(A) and null(B). It follows that $x \ = \ Z a$ for some $0 \neq a \in \mathbb { R } ^ { t }$ . But since $0 = B x = B Z a$ , we must have $a = W b$ for some $b \in \mathbb { R } ^ { q }$ . Thus, $x = Z W b \in \mathsf { r a n } ( Z W )$

If the SVD is used to compute the orthonormal bases in this theorem, then we obtain the following procedure:

Algorithm 6.4.2 Given $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { p \times n }$ , the following algorithm computes and integer s and a matrix $Y = [ y _ { 1 } \vert \cdot \cdot \cdot \vert y _ { s } ]$ having orthonormal columns which span $\mathsf { n u l l } ( A ) \cap \mathsf { n u l l } ( B )$ . If the intersection is trivial, then $s = 0$ .

Compute the SVD $U _ { A } ^ { T } A V _ { A } = \mathrm { d i a g } ( \sigma _ { i } )$ , save $V _ { A }$ , and set $r = \mathsf { r a n k } ( A )$ .

if $r < n$

$$
C = B V _ {A} (:, r + 1: n)
$$

Compute the SVD $U _ { c } ^ { T } C V _ { c } = \mathrm { { d i a g } } ( \gamma _ { i } )$ , save $V _ { C }$ , and set $q = \mathsf { r a n k } ( C )$ .

if q < n − r

$$
s = n - r - q
$$

$$
Y = V _ {A} (:, r + 1: n) V _ {C} (:, q + 1: n - r)
$$

else

$$
s = 0
$$

end

else

$$
s = 0
$$

end

The practical implementation of this algorithm requires an ability to reason about numerical rank. See §5.4.1.

# 6.4.3 Angles Between Subspaces

Let F and G be subspaces in $\mathbb { R } ^ { m }$ whose dimensions satisfy

$$
p = \dim (F) \geq \dim (G) = q \geq 1.
$$

The principal angles $\{ \theta _ { i } \} _ { i = 1 } ^ { q }$ between these two subspaces and the associated principal vectors $\{ f _ { 1 } , g _ { i } \} _ { i = 1 } ^ { q }$ are defined recursively by

$$
\begin{array}{l} \cos (\theta_ {k}) = f _ {k} ^ {T} g _ {k} = \max \quad \max \quad f ^ {T} g. \\ f \in F, \| f \| _ {2} = 1 \quad g \in G, \| g \| _ {2} = 1 \tag {6.4.3} \\ \end{array}
$$

$$
f ^ {T} [ f _ {1}, \dots , f _ {k - 1} ] = 0 \quad g ^ {T} [ g _ {1}, \dots , g _ {k - 1} ] = 0
$$

Note that the principal angles satisfy $0 \le \theta _ { 1 } \le \dots \le \theta _ { q } \le \pi / 2$ .. The problem of computing principal angles and vectors is oftentimes referred to as the canonical correlation problem.

Typically, the subspaces F and G are matrix ranges, e.g.,

$$
F = \operatorname{ran} (A), \qquad A \in \mathbb {R} ^ {n \times p},
$$

$$
G = \operatorname{ran} (B), \qquad B \in \mathbb {R} ^ {n \times q}.
$$

The principal vectors and angles can be computed using the QR factorization and the SVD. Let $A = Q _ { A } R _ { A }$ and $B = Q _ { B } R _ { B }$ be thin QR factorizations and assume that

$$
Q _ {A} ^ {T} Q _ {B} = Y \Sigma Z ^ {T} = \sum_ {i = 1} ^ {q} \sigma_ {i} y _ {i} z _ {i} ^ {T}
$$

is the SVD of $Q _ { A } ^ { T } Q _ { B } \in \mathbb { R } ^ { p \times q }$ . Since $\Vert Q _ { A } ^ { T } Q _ { B } \parallel _ { 2 } \leq 1$ , all the singular values are between 0 and 1 and we may write $\sigma _ { i } = \cos ( \theta _ { i } ) , i = 1 { : } q$ . Let

$$
Q _ {A} Y = \left[ f _ {1} \mid \dots \mid f _ {p} \right], \tag {6.4.4}
$$

$$
Q _ {B} Z = \left[ g _ {1} \mid \dots \mid g _ {q} \right] \tag {6.4.5}
$$

be column partitionings of the matrices $Q _ { A } Y \in \mathbb { R } ^ { n \times p }$ and $Q _ { B } Z \in \mathbb { R } ^ { n \times q }$ . These matrices have orthonormal columns. If $f \in F$ and $g \in G$ are unit vectors, then there exist unit vectors $u \in \mathbb { R } ^ { p }$ and $v \in \mathbb { R } ^ { q }$ so that $f = Q _ { A } u$ and $g = Q _ { B } v$ . Thus,

$$
\begin{array}{l} f ^ {T} g = (Q _ {A} u) ^ {T} (Q _ {B} v) = u ^ {T} (Q _ {A} ^ {T} Q _ {B}) v = u ^ {T} (Y \Sigma Z ^ {T}) v \\ = (Y ^ {T} u) ^ {T} \Sigma (Z ^ {T} v) = \sum_ {i = 1} ^ {q} \sigma_ {i} (y _ {i} ^ {T} u) (z _ {i} ^ {T} v). \tag {6.4.6} \\ \end{array}
$$

This expression attains its maximal value of $\sigma _ { 1 } = \cos ( \theta _ { 1 } )$ by setting $u = y _ { 1 }$ and $v = z _ { 1 }$ . It follows that $f = Q _ { A } y _ { 1 } = f _ { 1 }$ and $v = Q _ { B } z _ { 1 } = g _ { 1 }$ .

Now assume that k > 1 and that the first k − 1 columns of the matrices in (6.4.4) and (6.4.5) are known, $\mathrm { i . e . , ~ } f _ { 1 } , \ldots , f _ { k - 1 }$ and $g _ { 1 } , \ldots , g _ { k - 1 }$ . Consider the problem of maximizing $f ^ { T } g$ given that $f = Q _ { A } u$ and $g = Q _ { B } v$ are unit vectors that satisfy

$$
f ^ {T} \left[ f _ {1} \mid \dots \mid f _ {k - 1} \right] = 0,
$$

$$
g ^ {T} \left[ g _ {1} \mid \dots \mid g _ {k - 1} \right] = 0.
$$

It follows from (6.4.6) that

$$
f ^ {T} g = \sum_ {i = k} ^ {q} \sigma_ {i} (y _ {i} ^ {T} u) (z _ {i} ^ {T} v) \leq \sigma_ {k} \sum_ {i = k} ^ {q} | y _ {i} ^ {T} u | \cdot | z _ {i} ^ {T} v |.
$$

This expression attains its maximal value of $\sigma _ { k } = \cos ( \theta _ { k } )$ by setting $u = y _ { k }$ and $v = z _ { k }$ . It follows from (6.4.4) and (6.4.5) that $f = Q _ { A } y _ { k } = f _ { k }$ and $g = Q _ { B } z _ { k } = g _ { k }$ . Combining these observations we obtain

Algorithm 6.4.3 (Principal Angles and Vectors) Given $A \in \mathbb { R } ^ { m \times p }$ and $B \in \mathbb { R } ^ { m \times q }$ $( p \geq q )$ each with linearly independent columns, the following algorithm computes the cosines of the principal angles $\theta _ { 1 } \geq \cdots \geq \theta _ { q }$ between ran(A) and ran(B). The vectors $f _ { 1 } , \ldots , f _ { q }$ and $g _ { 1 } , \ldots , g _ { q }$ are the associated principal vectors.

Compute the thin QR factorizations $A = Q _ { A } R _ { A }$ and $B = Q _ { B } R _ { B }$ .

$$
C = Q _ {A} ^ {T} Q _ {B}
$$

Compute the SVD $Y ^ { T } C Z = \mathrm { d i a g } ( \cos ( \theta _ { k } ) )$ .

$$
Q _ {A} Y (:, 1: q) = \left[ f _ {1} \mid \dots \mid f _ {q} \right]
$$

$$
Q _ {B} Z (:, 1: q) = \left[ g _ {1} \mid \dots \mid g _ {q} \right]
$$

The idea of using the SVD to compute the principal angles and vectors is due to Bj¨orck and Golub (1973). The problem of rank deficiency in A and B is also treated in this paper. Principal angles and vectors arise in many important statistical applications. The largest principal angle is related to the notion of distance between equidimensional subspaces that we discussed in §2.5.3. If $p = q$ , then

$$
\operatorname{dist} (F, G) = \sqrt {1 - \cos (\theta_ {p}) ^ {2}} = \sin (\theta_ {p}).
$$

# 6.4.4 Intersection of Subspaces

In light of the following theorem, Algorithm 6.4.3 can also be used to compute an orthonormal basis for ran(A) ∩ ran(B) where $A \in \mathbb { R } ^ { m \times p }$ and $B \in \mathbb { R } ^ { m \times q }$

Theorem 6.4.2. Let $\{ \cos ( \theta _ { i } ) \} _ { i = 1 } ^ { q }$ and $\{ f _ { i } , g _ { i } \} _ { i = 1 } ^ { q }$ be defined by Algorithm 6.4.3. If the index s is defined by $1 = \cos ( \theta _ { 1 } ) = \cdot \cdot \cdot = \cos ( \theta _ { s } ) > \cos ( \theta _ { s + 1 } )$ , then

$$
\operatorname{ran} (A) \cap \operatorname{ran} (B) = \operatorname{span} \left\{f _ {1}, \dots , f _ {s} \right\} = \operatorname{span} \left\{g _ {1}, \dots , g _ {s} \right\}.
$$

Proof. The proof follows from the observation that if cos $( \theta _ { i } ) = 1$ , then $f _ { i } = g _ { i }$ .

The practical determination of the intersection dimension s requires a definition of what it means for a computed singular value to equal 1. For example, a computed singular value $\hat { \sigma } _ { i } = \cos ( \hat { \theta } _ { i } )$ could be regarded as a unit singular value if $\hat { \sigma } _ { i } \geq 1 - \delta$ for some intelligently chosen small parameter δ.

# Problems

P6.4.1 Show that if A and B are m-by-p matrices, with $p \leq m$ , then

$$
\min _ {Q ^ {T} Q = I _ {p}} \| A - B Q \| _ {F} ^ {2} = \sum_ {i = 1} ^ {p} (\sigma_ {i} (A) ^ {2} - 2 \sigma_ {i} (B ^ {T} A) + \sigma_ {i} (B) ^ {2}).
$$

P6.4.2 Extend Algorithm 6.4.2 so that it computes an orthonormal basis for $\mathsf { n u l l } ( A _ { 1 } ) \cap \cdots \cap \mathsf { n u l l } ( A _ { s } )$ where each matrix Ai has n columns.

P6.4.3 Extend Algorithm 6.4.3 so that it can handle the case when A and B are rank deficient.

P6.4.4 Verify Equation (6.4.2).

P6.4.5 Suppose A, $B \in \mathbb { R } ^ { m \times n }$ and that A has full column rank. Show how to compute a symmetric matrix $X \in \mathbb { R } ^ { n \times n }$ that minimizes $\| A X - B \| _ { F }$ . Hint: Compute the SVD of A.

P6.4.6 This problem is an exercise in F-norm optimization. (a) Show that if $C \in \mathbb { R } ^ { m \times n }$ and $e \in \mathbb { R } ^ { m }$ is a vector of ones, then $v = C ^ { T } e / m$ minimizes $\| \boldsymbol { C } - e \boldsymbol { v } ^ { T } \| _ { F }$ . (b) Suppose $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { m \times n }$ and that we wish to solve

$$
\min _ {Q ^ {T} Q = I _ {n}, v \in \mathbf {R} ^ {n}} \| A - (B + e v ^ {T}) Q \| _ {F}
$$

Show that $\boldsymbol { v } _ { \mathrm { o p t } } = ( A - B ) ^ { T } \boldsymbol { e } / m$ and $Q _ { \mathrm { o p t } } = U \Sigma V ^ { T }$ solve this problem where $B ^ { T } ( I - e e ^ { T } / m ) A = U V ^ { T }$ is the SVD.

P6.4.7 A 3-by-3 matrix H is ROPR matrix if $H = Q + x y ^ { T }$ where $Q \in \mathbb { R } ^ { 3 \times 3 }$ rotation and $x , y \in \mathbb { R } ^ { 3 }$ . (A rotation matrix is an orthogonal matrix with unit determinant. $\mathrm { \mathrm { \Omega ^ { * } R O P R } \mathrm { \Omega ^ { * } } }$ stands for “rank-1 perturbation of a rotation.”) ROPR matrices arise in computational photography and this problem highlights some of their properties. (a) If H is a ROPR matrix, then there exist rotations $U , V \in \mathbb { R } ^ { 3 \times 3 }$ , such that $\begin{array} { r } { U ^ { T } H V = \mathrm { d i a g } ( \sigma _ { 1 } , \sigma _ { 2 } , \sigma _ { 3 } ) } \end{array}$ satisfies $\sigma _ { 1 } \geq \sigma _ { 2 } \geq | \sigma _ { 3 } |$ . (b) Show that if $Q \in \mathbb { R } ^ { 3 \times 3 }$ is a rotation, then there exist cosine-sine pairs $( c _ { i } , s _ { i } ) = ( \cos ( \theta _ { i } ) , \sin ( \theta _ { i } ) ) , i = 1 { : } 3$ such that $Q = Q ( \theta _ { 1 } , \theta _ { 2 } , \theta _ { 3 } )$ where

$$
\begin{array}{l} Q (\theta_ {1}, \theta_ {2}, \theta_ {3}) = \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & c _ {1} & s _ {1} \\ 0 & - s _ {1} & c _ {1} \end{array} \right] \left[ \begin{array}{c c c} c _ {2} & s _ {2} & 0 \\ - s _ {2} & c _ {2} & 0 \\ 0 & 0 & 1 \end{array} \right] \left[ \begin{array}{c c c} 1 & 0 & 0 \\ 0 & c _ {3} & s _ {3} \\ 0 & - s _ {3} & c _ {3} \end{array} \right] \\ = \left[ \begin{array}{c c c} c _ {2} & s _ {2} c _ {3} & s _ {2} s _ {3} \\ - c _ {1} s _ {2} & c _ {1} c _ {2} c _ {3} - s _ {1} s _ {3} & c _ {1} c _ {2} s _ {3} + s _ {1} c _ {3} \\ s _ {1} s _ {2} & - s _ {1} c _ {2} c _ {3} - c _ {1} s _ {3} & - s _ {1} c _ {2} s _ {3} + c _ {1} c _ {3} \end{array} \right]. \\ \end{array}
$$

Hint: The Givens QR factorization involves three rotations. (c) Show that if

$$
\left[ \begin{array}{c c c} \sigma_ {1} & 0 & 0 \\ 0 & \sigma_ {2} & 0 \\ 0 & 0 & \sigma_ {3} \end{array} \right] = Q (\theta_ {1}, \theta_ {2}, \theta_ {3}) - x y ^ {T}, \qquad x, y \in \mathbb {R} ^ {3}
$$

then $x y ^ { T }$ must have the form

$$
x y ^ {T} = \left[ \begin{array}{c} s _ {2} \\ \mu c _ {1} \\ - \mu s _ {1} \end{array} \right] \left[ \begin{array}{c} - s _ {2} / \mu \\ c _ {3} \\ s _ {3} \end{array} \right] ^ {T}
$$

for some $\mu \geq 0$ and

$$
\left[ \begin{array}{c c} c _ {2} - \mu & 1 \\ 1 & c _ {2} - \mu \end{array} \right] \left[ \begin{array}{c} c _ {1} s _ {3} \\ s _ {1} c _ {3} \end{array} \right] = \left[ \begin{array}{c} 0 \\ 0 \end{array} \right].
$$

(d) Show that the second singular value of a ROPR matrix is 1.

P6.4.8 Let $U _ { * } \in \mathbb { R } ^ { n \times d }$ be a matrix with orthonormal columns whose span is a subspace S that we wish to estimate. Assume that $U _ { c } \in \mathbb { R } ^ { n \times d }$ is a given matrix with orthonormal columns and regard $\mathsf { r a n } ( U _ { c } )$ as the “current” estimate of $S _ { \cdot }$ . This problem examines what is required to get an improved estimate of S given the availability of a vector $v \in S$ . (a) Define the vectors

$$
w = U _ {c} ^ {T} v, \qquad v _ {1} = U _ {c} U _ {c} ^ {T} v, \qquad v _ {2} = (I _ {n} - U _ {c} U _ {c} ^ {T}) v,
$$

and assume that each is nonzero. (a) Show that if

$$
z _ {\theta} = \left(\frac {\cos (\theta) - 1}{\| v _ {1} \| \| w \|}\right) v _ {1} + \left(\frac {\sin (\theta)}{\| v _ {2} \| \| w \|}\right) v _ {2}
$$

and

$$
U _ {\theta} = (I _ {n} + z _ {\theta} v ^ {T}) U _ {c},
$$

then $U _ { \theta } ^ { T } U _ { \theta } = I _ { d }$ . Thus, $U _ { \theta } U _ { \theta } ^ { T }$ is an orthogonal projection. (b) Define the distance function

$$
\operatorname{dist} _ {F} (\operatorname{ran} (V), \operatorname{ran} (W)) = \| V V ^ {T} - W W ^ {T} \| _ {F}
$$

where $V , W \in \mathbb { R } ^ { n \times d }$ have orthonormal columns and show

$$
\mathsf {d i s t} _ {F} (\mathsf {r a n} (V), \mathsf {r a n} (W)) ^ {2} = 2 (d - \| W ^ {T} V \| _ {F} ^ {2}) = 2 \sum_ {i = 1} ^ {d} (1 - \sigma_ {i} (W ^ {T} V) ^ {2}).
$$

$\operatorname { N o t e \ t h a t \ d i s t } ( \mathsf { r a n } ( V ) , \mathsf { r a n } ( W ) ) ^ { 2 } = 1 - \sigma _ { 1 } ( W ^ { T } V ) ^ { 2 } . ~ ( \mathsf { c } ) \operatorname { S h o w \ t h a t }$

$$
d _ {\theta} ^ {2} = d _ {c} ^ {2} - 2 \cdot \mathbf {t r} (U _ {*} U _ {*} ^ {T} (U _ {\theta} U _ {\theta} ^ {T} - U _ {c} U _ {c} ^ {T}))
$$

where $d _ { \theta } = \mathsf { d i s t } _ { F } ( \mathsf { r a n } ( U _ { * } ) , \mathsf { r a n } ( U _ { \theta } ) )$ and $d _ { c } = { \tt d i s t } _ { F } ( { \tt r a n } ( U _ { * } ) , { \tt r a n } ( U _ { c } ) )$ . (d) Show that if

$$
y _ {\theta} = \cos (\theta) \frac {v _ {1}}{\parallel v _ {1} \parallel} + \sin (\theta) \frac {v _ {2}}{\parallel v _ {2} \parallel},
$$

then

$$
U _ {\theta} U _ {\theta} ^ {T} - U _ {c} U _ {c} ^ {T} = y _ {\theta} y _ {\theta} ^ {T} - \frac {v _ {1} v _ {1} ^ {T}}{v _ {1} ^ {T} v _ {1}}
$$

and

$$
d _ {\theta} ^ {2} = d _ {c} ^ {2} + 2 \left(\frac {\| U _ {*} ^ {T} v _ {1} \| _ {2} ^ {2}}{\| v _ {1} \| _ {2} ^ {2}} - \| U _ {*} ^ {T} y _ {\theta} \| _ {2} ^ {2}\right).
$$

(e) Show that if θ minimizes this quantity, then

$$
\sin (2 \theta) \left(\frac {\| P _ {S} v _ {2} \| ^ {2}}{\| v _ {2} \| _ {2} ^ {2}} - \frac {\| P _ {S} v _ {1} \| ^ {2}}{\| v _ {1} \| _ {2} ^ {2}}\right) + \cos (2 \theta) \frac {v _ {1} ^ {T} P _ {S} v _ {2}}{\| v _ {1} \| _ {2} \| v _ {2} \| _ {2}} = 0, \qquad P _ {S} = U _ {*} U _ {*} ^ {T}.
$$

# Notes and References for 6.4

References for the Procrustes problem include:

Using the SVD to solve the angles-between-subspaces problem is discussed in:   
B. Green (1952). “The Orthogonal Approximation of an Oblique Structure in Factor Analysis,” Psychometrika 17, 429–40.   
P. Schonemann (1966). “A Generalized Solution of the Orthogonal Procrustes Problem,” Psychometrika 31, 1–10.   
R.J. Hanson and M.J. Norris (1981). “Analysis of Measurements Based on the Singular Value Decomposition,” SIAM J. Sci. Stat. Comput. 2, 363–374.   
N.J. Higham (1988). “The Symmetric Procrustes Problem,” BIT 28, 133–43.   
H. Park (1991). “A Parallel Algorithm for the Unbalanced Orthogonal Procrustes Problem,” Parallel Comput. 17, 913–923.   
L.E. Andersson and T. Elfving (1997). “A Constrained Procrustes Problem,” SIAM J. Matrix Anal. Applic. 18, 124–139.   
L. Eld´en and H. Park (1999). “A Procrustes Problem on the Stiefel Manifold,” Numer. Math. 82, 599–619.   
A.W. Bojanczyk and A. Lutoborski (1999). “The Procrustes Problem for Orthogonal Stiefel Matrices,” SIAM J. Sci. Comput. 21, 1291–1304.   
If $B \ = \ I ,$ , then the Procrustes problem amounts to finding the closest orthogonal matrix. This computation is related to the polar decomposition problem that we consider in §9.4.3. Here are some basic references:   
˚A. Bj¨orck and C. Bowie (1971). “An Iterative Algorithm for Computing the Best Estimate of an Orthogonal Matrix,” SIAM J. Numer. Anal. 8, 358–64.   
N.J. Higham (1986). “Computing the Polar Decomposition with Applications,” SIAM J. Sci. Stat. Comput. 7, 1160–1174.   
˚A. Bj¨orck and G.H. Golub (1973). “Numerical Methods for Computing Angles Between Linear Subspaces,” Math. Comput. 27, 579–94.   
L.M. Ewerbring and F.T. Luk (1989). “Canonical Correlations and Generalized SVD: Applications and New Algorithms,” J. Comput. Appl. Math. 27, 37–52.   
G.H. Golub and H. Zha (1994). “Perturbation Analysis of the Canonical Correlations of Matrix Pairs,” Lin. Alg. Applic. 210, 3–28.

Z. Drmac (2000). “On Principal Angles between Subspaces of Euclidean Space,” SIAM J. Matrix Anal. Applic. 22, 173–194.   
A.V. Knyazev and M.E. Argentati (2002). “Principal Angles between Subspaces in an A–Based Scalar Product: Algorithms and Perturbation Estimates,” SIAM J. Sci. Comput. 23, 2008–2040.   
P. Strobach (2008). “Updating the Principal Angle Decomposition,” Numer. Math. 110, 83–112.   
In reduced-rank regression the object is to connect a matrix of signals to a matrix of noisey observations through a matrix that has specified low rank. An svd-based computational procedure that involves principal angles is discussed in:   
L. Eld´en and B. Savas (2005). “The Maximum Likelihood Estimate in Reduced-Rank Regression,” Num. Lin. Alg. Applic. 12, 731–741,   
The SVD has many roles to play in statistical computation, see:   
S.J. Hammarling (1985). “The Singular Value Decomposition in Multivariate Statistics,” ACM SIGNUM Newsletter 20, 2-25.   
An algorithm for computing the rotation and rank-one matrix in P6.4.7 that define a given ROPR matrix is discussed in:   
R. Schreiber, Z. Li, and H. Baker (2009). “Robust Software for Computing Camera Motion Parameters,” J. Math. Imaging Vision 33, 1–9.   
For a more details about the estimation problem associated with P6.4.8, see:   
L. Balzano, R. Nowak, and B. Recht (2010). “Online Identification and Tracking of Subspaces from Highly Incomplete Information,” Proceedings of the Allerton Conference on Communication, Control, and Computing 2010.

# 6.5 Updating Matrix Factorizations

In many applications it is necessary to refactor a given matrix $A \in \mathbb { R } ^ { m \times n }$ after it has undergone a small modification. For example, given that we have the QR factorization of a matrix A, we may require the QR factorization of the matrix A obtained from A by appending a row or column or deleting a row or column. In this section we show that in situations like these, it is much more efficient to “update” A’s QR factorization than to generate the required QR factorization of A from scratch. Givens rotations have a prominent role to play. In addition to discussing various update-QR strategies, we show how to downdate a Cholesky factorization using hyperbolic rotations and how to update a rank-revealing ULV decomposition.

# 6.5.1 Rank-1 Changes

Suppose we have the QR factorization $Q R = A \in \mathbb { R } ^ { n \times n }$ and that we need to compute the QR factorization $\widetilde { A } = A + u v ^ { T } = Q _ { 1 } R _ { 1 }$ where u, $v \in \mathbb { R } ^ { n }$ are given. Observe that

$$
\widetilde {A} = A + u v ^ {T} = Q (R + w v ^ {T}) \tag {6.5.1}
$$

where $w = Q ^ { T } u$ . Suppose rotations $J _ { n - 1 } , \ldots , J _ { 2 } , J _ { 1 }$ are computed such that

$$
J _ {1} ^ {T} \dots J _ {n - 1} ^ {T} w = \pm \| w \| _ {2} e _ {1}.
$$

where each $J _ { k }$ is a Givens rotation in planes k and k + 1. If these same rotations are applied to R, then

$$
H = J _ {1} ^ {T} \dots J _ {n - 1} ^ {T} R \tag {6.5.2}
$$

is upper Hessenberg. For example, in the $n = 4$ case we start with

$$
w \leftarrow \left[ \begin{array}{c} \times \\ \times \\ \times \\ \times \end{array} \right], \qquad R \leftarrow \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right],
$$

and then update as follows:

$$
\begin{array}{l} w \leftarrow J _ {3} ^ {T} w = \left[ \begin{array}{l} \times \\ \times \\ \times \\ 0 \end{array} \right], \qquad R \leftarrow J _ {3} ^ {T} R = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \end{array} \right], \\ w \leftarrow J _ {2} ^ {T} w = \left[ \begin{array}{c} \times \\ \times \\ 0 \\ 0 \end{array} \right], \qquad R \leftarrow J _ {2} ^ {T} R = \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \end{array} \right], \\ w \leftarrow J _ {1} ^ {T} w = \left[ \begin{array}{c} \times \\ 0 \\ 0 \\ 0 \end{array} \right], \qquad H \leftarrow J _ {1} ^ {T} R = \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \end{array} \right]. \\ \end{array}
$$

Consequently,

$$
(J _ {1} ^ {T} \dots J _ {n - 1} ^ {T}) (R + w v ^ {T}) = H \pm \| w \| _ {2} e _ {1} v ^ {T} = H _ {1} \tag {6.5.3}
$$

is also upper Hessenberg. Following Algorithm 5.2.4, we compute Givens rotations $G _ { k }$ , $k = 1 { : } n - 1$ such that $\check { G _ { n - 1 } } \cdot \cdot \cdot G _ { 1 } ^ { T } \check { H } _ { 1 } = R _ { 1 }$ is upper triangular. Combining everything we obtain the QR factorization $\widetilde A = A + u v ^ { T } \ = \ Q _ { 1 } R _ { 1 }$ where

$$
Q _ {1} = Q J _ {n - 1} \dots J _ {1} G _ {1} \dots G _ {n - 1}.
$$

A careful assessment of the work reveals that about $2 6 n ^ { 2 }$ flops are required.

The technique readily extends to the case when A is rectangular. It can also be generalized to compute the QR factorization of $A + U V ^ { T }$ where $U \in \mathbb { R } ^ { m \times p }$ and $V \in \mathbb { R } ^ { n \times p }$ .

# 6.5.2 Appending or Deleting a Column

Assume that we have the QR factorization

$$
Q R = A = \left[ a _ {1} \mid \dots \mid a _ {n} \right], \quad a _ {i} \in \mathbb {R} ^ {m}, \tag {6.5.4}
$$

and for some k, $1 \leq k \leq n$ , partition the upper triangular matrix $R \in \mathbb { R } ^ { m \times n }$ as follows:

$$
R = \left[ \begin{array}{c c c} R _ {1 1} & v & R _ {1 3} \\ 0 & r _ {k k} & w ^ {T} \\ 0 & 0 & R _ {3 3} \end{array} \right] \begin{array}{c} k - 1 \\ 1 \\ m - k \end{array} .
$$

Now suppose that we want to compute the QR factorization of

$$
\widetilde {A} = \left[ a _ {1} \mid \dots \mid a _ {k - 1} \mid a _ {k + 1} \mid \dots \mid a _ {n} \right] \in \mathbb {R} ^ {m \times (n - 1)}.
$$

Note that $\widetilde { A }$ is just A with its kth column deleted and that

$$
Q ^ {T} \widetilde {A} = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 3} \\ 0 & w ^ {T} \\ 0 & R _ {3 3} \end{array} \right] = H
$$

is upper Hessenberg, e.g.,

$$
H = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 \end{array} \right], \qquad m = 7, n = 6, k = 3.
$$

Clearly, the unwanted subdiagonal elements $h _ { k + 1 , k } , \ldots , h _ { n , n - 1 }$ can be zeroed by a sequence of Givens rotations: $G _ { n - 1 } ^ { T } \cdot \cdot \cdot G _ { k } ^ { T } H ~ = ~ R _ { 1 }$ . Here, $G _ { i }$ is a rotation in planes i and $i + 1$ for $i = k { : } n - 1$ . Thus, if $Q _ { 1 } = Q G _ { k } \cdot \cdot \cdot G _ { n - 1 }$ then $\widetilde { A } = Q _ { 1 } R _ { 1 }$ is the QR factorization of ${ \widetilde { A } } .$ .

The above update procedure can be executed in $O ( n ^ { 2 } )$ flops and is very useful in certain least squares problems. For example, one may wish to examine the significance of the kth factor in the underlying model by deleting the kth column of the corresponding data matrix and solving the resulting LS problem.

Analogously, it is possible to update efficiently the QR factorization of a matrix after a column has been added. Assume that we have (6.5.4) but now want the QR factorization of

$$
\widetilde {A} = \left[ a _ {1} \mid \ldots \mid a _ {k} \mid z \mid a _ {k + 1} \mid \ldots \mid a _ {n} \right]
$$

where $z \in \mathbb { R } ^ { m }$ is given. Note that if $w = Q ^ { T } z$ then

$$
Q ^ {T} \widetilde {A} = \left[ Q ^ {T} a _ {1} \mid \dots \mid Q ^ {T} a _ {k} \mid w \mid Q ^ {T} a _ {k + 1} \mid \dots \mid Q ^ {T} a _ {n} \right]
$$

is upper triangular except for the presence of a “spike” in its $( k + 1 )$ st column, e.g.,

$$
\widetilde {A} \leftarrow Q ^ {T} \widetilde {A} = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & 0 & \times \\ 0 & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & \times & 0 & 0 \end{array} \right], \qquad m = 7, n = 5, k = 3.
$$

It is possible to determine a sequence of Givens rotations that restores the triangular form:

$$
\widetilde {A} \leftarrow J _ {6} ^ {T} \widetilde {A} = \left[ \begin{array}{l l l l l l} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & 0 & \times \\ 0 & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 \end{array} \right], \qquad \widetilde {A} \leftarrow J _ {5} ^ {T} \widetilde {A} = \left[ \begin{array}{l l l l l l} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & 0 & \times \\ 0 & 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 & 0 \end{array} \right],
$$

$$
\widetilde {A} \leftarrow J _ {4} ^ {T} \widetilde {A} = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 & 0 \end{array} \right].
$$

This update requires O(mn) flops.

# 6.5.3 Appending or Deleting a Row

Suppose we have the QR factorization $Q R = A \in \mathbb { R } ^ { m \times n }$ and now wish to obtain the QR factorization of

$$
\widetilde {A} = \left[ \begin{array}{c} w ^ {T} \\ A \end{array} \right]
$$

where $w \in \mathbb { R } ^ { n }$ . Note that

$$
\mathrm{diag} (1, Q ^ {T}) \widetilde {A} = \left[ \begin{array}{c} w ^ {T} \\ R \end{array} \right] = H
$$

is upper Hessenberg. Thus, rotations $J _ { 1 } , \ldots , J _ { n }$ can be determined so $J _ { n } ^ { T } \cdot \cdot \cdot J _ { 1 } ^ { T } H =$ $R _ { 1 }$ is upper triangular. It follows that $\begin{array} { r c l } { \bar { \cal A } } & { = } & { { \cal Q } _ { 1 } { \cal R } _ { 1 } } \end{array}$ is the desired QR factorization, where $Q _ { 1 } = \operatorname { d i a g } ( 1 , Q ) J _ { 1 } \cdot \cdot \cdot J _ { n }$ . See Algorithm 5.2.5.

No essential complications result if the new row is added between rows k and $k + 1$ of A. Indeed, if

$$
\left[ \begin{array}{l} A _ {1} \\ A _ {2} \end{array} \right] = Q R, \qquad A _ {1} \in \mathbb {R} ^ {k \times n},   A _ {2} \in \mathbb {R} ^ {(m - k) \times n},
$$

and

$$
P = \left[ \begin{array}{c c c} 0 & 1 & 0 \\ I _ {k} & 0 & 0 \\ 0 & 0 & I _ {m - k} \end{array} \right],
$$

then

$$
\operatorname{diag} (1, Q ^ {T}) P \left[ \begin{array}{c} A _ {1} \\ w ^ {T} \\ A _ {2} \end{array} \right] = \left[ \begin{array}{c} w ^ {T} \\ R \end{array} \right] = H
$$

is upper Hessenberg and we proceed as before.

Lastly, we consider how to update the QR factorization $Q R = A \in \mathbb { R } ^ { m \times n }$ when the first row of A is deleted. In particular, we wish to compute the QR factorization of the submatrix $A _ { 1 }$ in

$$
A = \left[ \begin{array}{l} z ^ {T} \\ A _ {1} \end{array} \right] _ {m - 1} ^ {1}.
$$

(The procedure is similar when an arbitrary row is deleted.) Let $q ^ { T }$ be the first row of Q and compute Givens rotations $G _ { 1 } , \ldots , G _ { m - 1 }$ such that $G _ { 1 } ^ { T } \cdot \cdot \cdot G _ { m - 1 } ^ { T } q = \alpha e _ { 1 }$ where $\alpha = \pm 1$ . Note that

$$
H = G _ {1} ^ {T} \dots G _ {m - 1} ^ {T} R = \left[ \begin{array}{l} v ^ {T} \\ R _ {1} \end{array} \right] _ {m - 1} ^ {1}
$$

is upper Hessenberg and that

$$
Q G _ {m - 1} \dots G _ {1} = \left[ \begin{array}{l l} \alpha & 0 \\ 0 & Q _ {1} \end{array} \right]
$$

where $Q _ { 1 } \in \mathbb { R } ^ { ( m - 1 ) \times ( m - 1 ) }$ is orthogonal. Thus,

$$
A = \left[ \begin{array}{l} z ^ {T} \\ A _ {1} \end{array} \right] = (Q G _ {m - 1} \dots G _ {1}) (G _ {1} ^ {T} \dots G _ {m - 1} ^ {T} R) = \left[ \begin{array}{l l} \alpha & 0 \\ 0 & Q _ {1} \end{array} \right] \left[ \begin{array}{l} v ^ {T} \\ R _ {1} \end{array} \right]
$$

from which we conclude that $A _ { 1 } = Q _ { 1 } R _ { 1 }$ is the desired QR factorization.

# 6.5.4 Cholesky Updating and Downdating

Suppose we are given a symmtetric positive definite matrix $A \in \mathbb { R } ^ { n \times n }$ and its Cholesky factor G. In the Cholesky updating problem, the challenge is to compute the Cholesky factorization $\widetilde { A } = \widetilde { G } \widetilde { G } ^ { T }$ where

$$
\widetilde {A} = A + z z ^ {T}, \quad z \in \mathbb {R} ^ {n}. \tag {6.5.5}
$$

Noting that

$$
\widetilde {A} = \left[ \begin{array}{c} G ^ {T} \\ z ^ {T} \end{array} \right] ^ {T} \left[ \begin{array}{c} G ^ {T} \\ z ^ {T} \end{array} \right], \tag {6.5.6}
$$

we can solve this problem by computing a product of Givens rotations $Q = Q _ { 1 } \cdot \cdot \cdot Q _ { n }$ so that

$$
Q ^ {T} \left[ \begin{array}{l} G ^ {T} \\ z ^ {T} \end{array} \right] = \left[ \begin{array}{l} R \\ 0 \end{array} \right], \quad R \in \mathbb {R} ^ {n \times n} \tag {6.5.7}
$$

is upper triangular. It follows that $\widetilde { A } = R R ^ { T }$ and so the updated Cholesky factor is given by $\widetilde { G } = \overline { { R } } ^ { T }$ . The zeroing sequence that produces R is straight forward, e.g.,

$$
\left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & 0 & \times \\ \times & \times & \times \end{array} \right] \xrightarrow {Q _ {1}} \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & 0 & \times \\ 0 & \times & \times \end{array} \right] \xrightarrow {Q _ {2}} \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & 0 & \times \\ 0 & 0 & \times \end{array} \right] \xrightarrow {Q _ {3}} \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & 0 & \times \\ 0 & 0 & 0 \end{array} \right].
$$

The $Q _ { k }$ update involves only rows k and $n + 1$ . The overall process is essentially the same as the strategy we outlined in the previous subsection for updating the QR factorization of a matrix when a row is appended.

The Cholesky downdating problem involves a different set of tools and a new set of numerical concerns. We are again given a Cholesky factorization $A = G G ^ { T }$ and a vector $z \in \mathbb { R } ^ { n }$ . However, now the challenge is to compute the Cholesky factorization $\widetilde { A } = \widetilde { G } \widetilde { G } ^ { T }$ where

$$
\widetilde {A} = A - z z ^ {T} \tag {6.5.8}
$$

is presumed to be positive definite. By introducing the notion of a hyperbolic rotation we can develop a downdating framework that corresponds to the Givens-based updating framework. Define the matrix S as follows

$$
S = \left[ \begin{array}{c c} I _ {n} & 0 \\ 0 & - 1 \end{array} \right] \tag {6.5.9}
$$

and note that

$$
\widetilde {A} = G G ^ {T} - z z ^ {T} = \left[ \begin{array}{c} G ^ {T} \\ z ^ {T} \end{array} \right] ^ {T} S \left[ \begin{array}{c} G ^ {T} \\ z ^ {T} \end{array} \right]. \tag {6.5.10}
$$

This corresponds to (6.5.6), but instead of computing the QR factorization (6.5.7), we seek a matrix $H \in \mathbb { R } ^ { ( n + 1 ) \times ( n + 1 ) }$ that satisfies two properties:

$$
H S H ^ {T} = S, \tag {6.5.11}
$$

$$
H ^ {T} \left[ \begin{array}{c} G ^ {T} \\ z ^ {T} \end{array} \right] = \left[ \begin{array}{c} R \\ 0 \end{array} \right], \quad R \in \mathbb {R} ^ {n \times n} (\text {upper triangular}). \tag {6.5.12}
$$

If this can be accomplished, then it follows from

$$
\widetilde {A} = \left(H ^ {T} \left[ \begin{array}{l} G ^ {T} \\ z ^ {T} \end{array} \right]\right) ^ {T} \left[ \begin{array}{l l} I _ {n} & 0 \\ 0 & - 1 \end{array} \right] \left(H ^ {T} \left[ \begin{array}{l} G ^ {T} \\ z ^ {T} \end{array} \right]\right) = R ^ {T} R
$$

that the Cholesky factor of $\widetilde { A } = A - z z ^ { T }$ is given by $\widetilde { G } = R ^ { T }$ . A matrix H that satisfies (6.5.11) is said to be S-orthogonal. Note that the product of S-orthogonal matrices is also S-orthogonal.

An important subset of the S-orthogonal matrices are the hyperbolic rotations and here is a 4-by-4 example:

$$
H _ {2} (\theta) = \left[ \begin{array}{c c c c} 1 & 0 & 0 & 0 \\ 0 & c & 0 & - s \\ 0 & 0 & 1 & 0 \\ 0 & - s & 0 & c \end{array} \right], \qquad c = \cosh (\theta), s = \sinh (\theta).
$$

The S-orthogonality of this matrix follows from cosh $( \theta ) ^ { 2 } - \sinh ( \theta ) ^ { 2 } = 1$ . In general, $H _ { k } \in \mathbb { R } ^ { ( n + 1 ) \times ( n + 1 ) }$ is a hyperbolic rotation if it agrees with $I _ { n + 1 }$ except in four locations:

$$
\left[ \begin{array}{c c} [ H _ {k} ] _ {k, k} & [ H _ {k} ] _ {k, n + 1} \\ [ H _ {k} ] _ {n + 1, k} & [ H _ {k} ] _ {n + 1, n + 1} \end{array} \right] = \left[ \begin{array}{c c} \cosh (\theta) & - \sinh (\theta) \\ - \sinh (\theta) & \cosh (\theta) \end{array} \right].
$$

Hyperbolic rotations look like Givens rotations and, not surprisingly, can be used to introduce zeros into a vector or matrix. However, upon consideration of the equation

$$
\left[ \begin{array}{c c} {c} & {- s} \\ {- s} & {c} \end{array} \right] \left[ \begin{array}{l} {x _ {1}} \\ {x _ {2}} \end{array} \right] = \left[ \begin{array}{l} {r} \\ {0} \end{array} \right], \qquad c ^ {2} - s ^ {2} = 1
$$

we see that the required cosh-sinh pair may not exist. Since we always have $| \cosh ( \theta ) | >$ $| \sinh ( \theta ) |$ , there is no real solution t $\mathrm { ~ o ~ } - s x _ { 1 } + c x _ { 2 } = 0 \mathrm { ~ i f ~ } | x _ { 2 } | > | x _ { 1 } |$ . On the other hand, $\mathrm { i f } \ | x _ { 1 } | > | x _ { 2 } |$ , then $\{ c , s \} = \{ \cosh ( \theta ) , \sinh ( \theta ) \}$ can be computed as follows:

$$
\tau = \frac {x _ {2}}{x _ {1}}, \quad c = \frac {1}{\sqrt {1 - \tau^ {2}}}, \quad s = c \cdot \tau . \tag {6.5.13}
$$

There are clearly numerical issues $\mathrm { i f } \ | x _ { 1 } |$ is just slightly greater than $| x _ { 2 } |$ . However, it is possible to organize hyperbolic rotation computations successfully, see Alexander, Pan, and Plemmons (1988).

Putting these concerns aside, we show how the matrix H in (6.5.12) can be computed as a product of hyperbolic rotations $H = H _ { 1 } \cdot \cdot \cdot H _ { n }$ just as the transforming Q in the updating problem is a product of Givens rotations. Consider the role of $H _ { 1 }$ in the $n = 3$ case:

$$
\left[ \begin{array}{c c c c} c & 0 & 0 & - s \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ - s & 0 & 0 & c \end{array} \right] ^ {T} \left[ \begin{array}{c c c} g _ {1 1} & g _ {2 1} & g _ {3 1} \\ 0 & g _ {2 2} & g _ {3 2} \\ 0 & 0 & g _ {3 3} \\ z _ {1} & z _ {2} & z _ {3} \end{array} \right] = \left[ \begin{array}{c c c} \tilde {g} _ {1 1} & \tilde {g} _ {2 1} & \tilde {g} _ {3 1} \\ 0 & g _ {2 2} & g _ {3 2} \\ 0 & 0 & g _ {3 3} \\ 0 & z _ {2} ^ {\prime} & z _ {3} ^ {\prime} \end{array} \right].
$$

Since $\widetilde { A } = G G ^ { T } - z z ^ { T }$ is positive definite, $[ \widetilde { A } ] _ { 1 1 } = g _ { 1 1 } ^ { 2 } - z _ { 1 } ^ { 2 } > 0$ . It follows that $\vert g _ { 1 1 } \vert > \vert z _ { 1 } \vert$ which guarantees that the cosh-sinh computations (6.5.13) go through. For the overall process to be defined, we have to guarantee that hyperbolic rotations $H _ { 2 } , \ldots , H _ { n }$ can be found to zero out the bottom row in the matrix $\big [ G ^ { T } ~ z \big ] ^ { T }$ . The following theorem ensures that this is the case.

Theorem 6.5.1. If

$$
A = \left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right] = \left[ \begin{array}{c c} g _ {1 1} & 0 \\ g _ {1} & G _ {1} \end{array} \right] \left[ \begin{array}{c c} g _ {1 1} & g _ {1} ^ {T} \\ 0 & G _ {1} ^ {T} \end{array} \right]
$$

and

$$
\widetilde {A} = A - z z ^ {T} = A - \left[ \begin{array}{l} \mu \\ w \end{array} \right] \left[ \begin{array}{l} \mu \\ w \end{array} \right] ^ {T}
$$

are positive definite, then it is possible to determine $c = \cosh ( \theta )$ and $s = \sinh ( \theta )$ so

$$
\left[ \begin{array}{c c c} c & 0 & - s \\ 0 & I _ {n - 1} & 0 \\ - s & 0 & c \end{array} \right] \left[ \begin{array}{c c} g _ {1 1} & g _ {1} ^ {T} \\ 0 & G _ {1} ^ {T} \\ \mu & w ^ {T} \end{array} \right] = \left[ \begin{array}{c c} \tilde {g} _ {1 1} & \tilde {g} _ {1} ^ {T} \\ 0 & G _ {1} ^ {T} \\ 0 & w _ {1} ^ {T} \end{array} \right].
$$

Moreover, the matrix $\widetilde { A } _ { 1 } = G _ { 1 } G _ { 1 } ^ { T } - w _ { 1 } w _ { 1 } ^ { T }$ is positive definite.

Proof. The blocks in A’s Cholesky factor are given by

$$
g _ {1 1} = \sqrt {\alpha}, \quad g _ {1} = v / g _ {1 1}, \quad G _ {1} G _ {1} ^ {T} = B - \frac {1}{\alpha} v v ^ {T}. \tag {6.5.14}
$$

Since $A - z z ^ { T }$ is positive definite, $a _ { 1 1 } - z _ { 1 } ^ { 2 } = g _ { 1 1 } ^ { 2 } - \mu ^ { 2 } > 0$ and so from (6.5.13) with $\tau = \mu / g _ { 1 1 }$ we see that

$$
c = \frac {\sqrt {\alpha}}{\sqrt {\alpha - \mu^ {2}}}, \quad s = \frac {\mu}{\sqrt {\alpha - \mu^ {2}}}. \tag {6.5.15}
$$

Since $w _ { 1 } = - s g _ { 1 } + c w$ it follows from (6.5.14) and (6.5.15) that

$$
\begin{array}{l} \widetilde {A} _ {1} = G _ {1} G _ {1} ^ {T} - w _ {1} w _ {1} ^ {T} = B - \frac {1}{\alpha} v v ^ {T} - (- s g _ {1} + c w) (- s g _ {1} + c w) ^ {T} \\ = B - \frac {c ^ {2}}{\alpha} v v ^ {T} - c ^ {2} w w ^ {T} + \frac {s c}{\sqrt {\alpha}} (v w ^ {T} + w v ^ {T}) \\ = B - \frac {1}{\alpha - \mu^ {2}} v v ^ {T} - \frac {\alpha}{\alpha - \mu^ {2}} w w ^ {T} + \frac {\mu}{\alpha - \mu^ {2}} (v w ^ {T} + w v ^ {T}). \\ \end{array}
$$

It is easy to verify that this matrix is precisely the Schur complement of α in

$$
\widetilde {A} = A - z z ^ {T} = \left[ \begin{array}{l l} \alpha - \mu^ {2} & v ^ {T} - \mu w ^ {T} \\ v - \mu w & B - w w ^ {T} \end{array} \right]
$$

and is therefore positive definite.

The theorem provides the key step in an induction proof that the factorization (6.5.12) exists.

# 6.5.5 Updating a Rank-Revealing ULV Decomposition

We close with a discussion about updating a nullspace basis after one or more rows have been appended to the underlying matrix. We work with the ULV decomposition which is much more tractable than the SVD from the updating point of view. We pattern our remarks after Stewart(1993).

A rank -revealing ULV decomposition of a matrix $A \in \mathbb { R } ^ { m \times n }$ has the form

$$
U ^ {T} A V = \left[ \begin{array}{l} L \\ 0 \end{array} \right] = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \\ 0 & 0 \end{array} \right], \quad U ^ {T} U = I _ {m}, V ^ {T} V = I _ {n} \tag {6.5.16}
$$

where $L _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ and $L _ { 2 2 } \in \mathbb { R } ^ { ( n - r ) \times ( n - r ) }$ are lower triangular and $\parallel L _ { 2 1 } \parallel _ { 2 }$ and $\parallel L _ { 2 2 } \parallel _ { 2 }$ are small compared to $\sigma _ { \mathrm { m i n } } ( L _ { 1 1 } )$ . Such a decomposition can be obtained by applying QR with column pivoting

$$
U ^ {T} A \Pi = \left[ \begin{array}{c} R \\ 0 \end{array} \right], \qquad R \in \mathbb {R} ^ {n \times n}
$$

followed by a QR factorization $V _ { 1 } ^ { T } R ^ { T } = L ^ { T }$ . In this case the matrix V in (6.5.16) is given by $V = \Pi V _ { 1 }$ . The parameter r is the estimated rank. Note that if

$$
V = \left[ \begin{array}{c c} V _ {1} & V _ {2} \\ r & n - r \end{array} \right], \qquad U = \left[ \begin{array}{c c} U _ {1} & U _ {2} \\ r & m - r \end{array} \right],
$$

then the columns of $V _ { 2 }$ define an approximate nullspace:

$$
\left\| A V _ {2} \right\| _ {2} = \left\| U _ {2} L _ {2 2} \right\| _ {2} = \left\| L _ {2 2} \right\| _ {2}.
$$

Our goal is to produce cheaply a rank-revealing ULV decomposition for the rowappended matrix

$$
\tilde {A} = \left[ \begin{array}{c} A \\ z ^ {T} \end{array} \right],
$$

In particular, we show how to revise L, V , and possibly r in $O ( n ^ { 2 } )$ flops. Note that

$$
\left[ \begin{array}{l l} U & 0 \\ 0 & 1 \end{array} \right] ^ {T} \left[ \begin{array}{l} A \\ z ^ {T} \end{array} \right] V = \left[ \begin{array}{l l} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \\ 0 & 0 \\ w ^ {T} & y ^ {T} \end{array} \right].
$$

We illustrate the key ideas through an example. Suppose $n = 7$ and $r \ = \ 4$ . By permuting the rows so that the bottom row is just underneath L, we obtain

$$
\left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \\ w ^ {T} & y ^ {T} \end{array} \right] = \left[ \begin{array}{c c c c c c c} \ell & 0 & 0 & 0 & 0 & 0 & 0 \\ \ell & \ell & 0 & 0 & 0 & 0 & 0 \\ \ell & \ell & \ell & 0 & 0 & 0 & 0 \\ \ell & \ell & \ell & \ell & 0 & 0 & 0 \\ \hline \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & 0 & 0 \\ \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & 0 \\ \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon \\ \hline w & w & w & w & y & y & y \end{array} \right].
$$

The 
 entries are small while the , w, and y entries are not. Next, a sequence of Givens rotations $G _ { 7 } , \ldots , G _ { 1 }$ are applied from the left to zero out the bottom row:

$$
\left[ \frac {\tilde {L}}{0} \right] = \left[ \begin{array}{l l l l l l l} \times & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 & 0 \\ \times & \times & \times & \times & 0 & 0 & 0 \\ \times & \times & \times & \times & \times & 0 & 0 \\ \times & \times & \times & \times & \times & \times & 0 \\ \times & \times & \times & \times & \times & \times & \times \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{array} \right] = G _ {1 7} \dots G _ {5 7} G _ {6 7} \left[ \begin{array}{l l} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \\ w ^ {T} & y ^ {T} \end{array} \right].
$$

Because this zeroing process intermingles the (presumably large) entries of the bottom row with the entries from each of the other rows, the lower triangular form is typically not rank revealing. However, and this is key, we can restore the rank-revealing structure with a combination of condition estimation and Givens zero chasing.

Let us assume that with the added row, the new nullspace has dimension 2. With a reliable condition estimator we produce a unit 2-norm vector p such that

$$
\parallel p ^ {T} \widetilde {L} \parallel_ {2} \approx \sigma_ {\mathrm{min}} (\widetilde {L}).
$$

(See §3.5.4). Rotations $\{ U _ { i , i + 1 } \} _ { i = 1 } ^ { 6 }$ can be found such that

$$
U _ {6 7} ^ {T} U _ {5 6} ^ {T} U _ {4 5} ^ {T} U _ {3 4} ^ {T} U _ {2 3} ^ {T} U _ {1 2} ^ {T} p = e _ {7} = I _ {7} (:, 7).
$$

Applying these rotations to $\widetilde { L }$ produces a lower Hessenberg matrix

$$
H = U _ {6 7} ^ {T} U _ {5 6} ^ {T} U _ {4 5} ^ {T} U _ {3 4} ^ {T} U _ {2 3} ^ {T} U _ {1 2} ^ {T} \tilde {L}.
$$

Applying more rotations from the right restores H to a lower triangular form:

$$
L _ {+} = H V _ {1 2} V _ {2 3} V _ {3 4} V _ {4 5} V _ {5 6} V _ {6 7}.
$$

It follows that

$$
e _ {7} ^ {T} L _ {+} = \left(e _ {8} ^ {T} H\right) V _ {1 2} V _ {2 3} V _ {3 4} V _ {4 5} V _ {5 6} V _ {6 7} = \left(p ^ {T} \tilde {L}\right) V _ {1 2} V _ {2 3} V _ {3 4} V _ {4 5} V _ {5 6} V _ {6 7}
$$

has approximate norm $\sigma _ { \mathrm { m i n } } ( \widetilde { L } )$ . Thus, we obtain a lower triangular matrix of the form

$$
L _ {+} = \left[ \begin{array}{c c c c c c c} \times & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 & 0 \\ \times & \times & \times & \times & 0 & 0 & 0 \\ \times & \times & \times & \times & \times & 0 & 0 \\ \times & \times & \times & \times & \times & \times & 0 \\ \hline \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon \end{array} \right]
$$

We can repeat the condition estimation and zero chasing on the leading 6-by-6 portion. Assuming that the nullspace of the augmented matrix has dimension two, this produces another row of small numbers:

$$
\left[ \begin{array}{c c c c c c c} \times & 0 & 0 & 0 & 0 & 0 & 0 \\ \times & \times & 0 & 0 & 0 & 0 & 0 \\ \times & \times & \times & 0 & 0 & 0 & 0 \\ \times & \times & \times & \times & 0 & 0 & 0 \\ \times & \times & \times & \times & \times & 0 & 0 \\ \hline \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & 0 \\ \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon & \epsilon \end{array} \right].
$$

This illustrates how we can restore any lower triangular matrix to rank-revealing form.

# Problems

P6.5.1 Suppose we have the QR factorization for $A \in \mathbb { R } ^ { m \times n }$ and now wish to solve

$$
\min _ {x \in \mathbf {R} ^ {n}} \| (A + u v ^ {T}) x - b \| _ {2}
$$

where $u , b \in \mathbb { R } ^ { m }$ and $v \in \mathbb { R } ^ { n }$ are given. Give an algorithm for solving this problem that requires $O ( m n )$ flops. Assume that $Q$ must be updated.

P6.5.2 Suppose

$$
A   =   \left[ \begin{array}{c} c ^ {T} \\ B \end{array} \right], \qquad c \in \mathbb {R} ^ {n},   B \in \mathbb {R} ^ {(m - 1) \times n}
$$

has full column rank and $m > n$ . Using the Sherman-Morrison-Woodbury formula show that

$$
\frac {1}{\sigma_ {\min} (B)} \leq \frac {1}{\sigma_ {\min} (A)} + \frac {\| (A ^ {T} A) ^ {- 1} c \| _ {2} ^ {2}}{1 - c ^ {T} (A ^ {T} A) ^ {- 1} c}.
$$

P6.5.3 As a function of $x _ { 1 }$ and $x _ { 2 }$ , what is the 2-norm of the hyperbolic rotation produced by (6.5.13)?

P6.5.4 Assume that

$$
A   =   \left[ \begin{array}{c c} R & H \\ 0 & E \end{array} \right], \qquad \rho   =   \frac {\|   E   \| _ {2}}{\sigma_ {\min} (R)}   <     1,
$$

where R and E are square. Show that if

$$
Q = \left[ \begin{array}{l l} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right]
$$

is orthogonal and

$$
\left[ \begin{array}{c c} R & H \\ 0 & E \end{array} \right] \left[ \begin{array}{c c} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right] = \left[ \begin{array}{c c} R _ {1} & 0 \\ H _ {1} & E _ {1} \end{array} \right],
$$

then  $H _ { 1 } \parallel _ { 2 } \leq \rho \parallel H \parallel _ { 2 }$ .

P6.5.5 Suppose $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ with $m \geq n$ . In the indefinite least squares (ILS) problem, the goal is to minimize

$$
\phi (x) = (b - A x) ^ {T} J (b - A x),
$$

where

$$
S = \left[ \begin{array}{c c} I _ {p} & 0 \\ 0 & - I _ {q} \end{array} \right], \qquad p + q = m.
$$

It is assumed that $p \geq 1$ and $q \geq 1$ . (a) By taking the gradient of $\phi ,$ , show that the ILS problem has a unique solution if and only if $A ^ { T } S A$ is positive definite. (b) Assume that the ILS problem has a unique solution. Show how it can be found by computing the Cholesky factorization of $\bar { Q } _ { 1 } ^ { T } Q _ { 1 } - Q _ { 2 } ^ { T } Q _ { 2 }$ where

$$
A   =   \left[ \begin{array}{l} Q _ {1} \\ Q _ {2} \end{array} \right], \qquad Q _ {1} \in \mathbb {R} ^ {p \times n},   Q _ {2} \in \mathbb {R} ^ {q \times n}
$$

is the thin QR factorization. (c) A matrix $Q \in \mathbb { R } ^ { m \times m }$ is S-orthogonal if $Q S Q ^ { T } = S { \mathrm { ~ I f } } $

$$
Q = \left[ \begin{array}{c c} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right] _ {q} ^ {p}
$$

is S-orthogonal, then by comparing blocks in the equation $Q ^ { T } S Q = S$ we have

$$
Q _ {1 1} ^ {T} Q _ {1 1} = I _ {p} + Q _ {2 1} ^ {T} Q _ {2 1}, \qquad Q _ {1 1} ^ {T} Q _ {1 2} = Q _ {2 1} ^ {T} Q _ {2 2}, \qquad Q _ {2 2} ^ {T} Q _ {2 2} = I _ {q} + Q _ {1 2} ^ {T} Q _ {1 2}.
$$

Thus, the singular values of $Q _ { 1 1 }$ and $Q _ { 2 2 }$ are never smaller than 1. Assume that $p \geq q$ . By analogy with how the CS decomposition is established in §2.5.4, show that there exist orthogonal matrices $U _ { 1 }$ , $U _ { 2 } , V _ { 1 }$ and $V _ { 2 }$ such that

$$
\left[\begin{array}{c c}U _ {1}&0\\0&U _ {2}\end{array}\right] ^ {T} Q \left[\begin{array}{c c}V _ {1}&0\\0&V _ {2}\end{array}\right] = \left[ \right.\begin{array}{c c}D&0\\0&I _ {p - q}\\\hline (D ^ {2} - I _ {p}) ^ {1 / 2}&0\end{array}\left. \right|\begin{array}{c}(D ^ {2} - I) ^ {1 / 2}\\0\\D\end{array}\left. \right]
$$

where $D = \operatorname { d i a g } ( d _ { 1 } , \dotsc , d _ { p } )$ with $d _ { i } \geq 1 , i = 1 { : } p$ . This is the hyperbolic CS decomposition and details can be found in Stewart and Van Dooren (2006).

# Notes and References for 6.5

The seminal matrix factorization update paper is:

P.E. Gill, G.H. Golub, W. Murray, and M.A. Saunders (1974). “Methods for Modifying Matrix Factorizations,” Math. Comput. 28, 505–535.

Initial research into the factorization update problem was prompted by the development of quasi-Newton methods and the simplex method for linear programming. In these venues, a linear system must be solved in step k that is a low-rank perturbation of the linear system solved in step k − 1, see:

R.H. Bartels (1971). “A Stabilization of the Simplex Method,” Numer. Math. 16, 414–434.

P.E. Gill, W. Murray, and M.A. Saunders (1975). “Methods for Computing and Modifying the LDV Factors of a Matrix,” Math. Comput. 29, 1051–1077.

D. Goldfarb (1976). “Factored Variable Metric Methods for Unconstrained Optimization,” Math. Comput. 30, 796–811.

J.E. Dennis and R.B. Schnabel (1983). Numerical Methods for Unconstrained Optimization and Nonlinear Equations, Prentice-Hall, Englewood Cliffs, NJ.

W.W. Hager (1989). “Updating the Inverse of a Matrix,” SIAM Review 31, 221–239.

S.K. Eldersveld and M.A. Saunders (1992). “A Block-LU Update for Large-Scale Linear Programming,” SIAM J. Matrix Anal. Applic. 13, 191–201.

Updating issues in the least squares setting are discussed in:

J. Daniel, W.B. Gragg, L. Kaufman, and G.W. Stewart (1976). “Reorthogonaization and Stable Algorithms for Updating the Gram-Schmidt QR Factorization,” Math. Comput. 30, 772–795.

S. Qiao (1988). “Recursive Least Squares Algorithm for Linear Prediction Problems,” SIAM J. Matrix Anal. Applic. 9, 323–328.

˚A. Bj¨orck, H. Park, and L. Eld´en (1994). “Accurate Downdating of Least Squares Solutions,” SIAM J. Matrix Anal. Applic. 15, 549–568.

S.J. Olszanskyj, J.M. Lebak, and A.W. Bojanczyk (1994). “Rank-k Modification Methods for Recursive Least Squares Problems,” Numer. Alg. 7, 325–354.

L. Eld´en and H. Park (1994). “Block Downdating of Least Squares Solutions,” SIAM J. Matrix Anal. Applic. 15, 1018–1034.

Kalman filtering is a very important tool for estimating the state of a linear dynamic system in the presence of noise. An illuminating, stable implementation that involves updating the QR factorization of an evolving block banded matrix is given in:

C.C. Paige and M.A. Saunders (1977). “Least Squares Estimation of Discrete Linear Dynamic Systems Using Orthogonal Transformations,” SIAM J. Numer. Anal. 14, 180–193.

The Cholesky downdating literature includes:

G.W. Stewart (1979). “The Effects of Rounding Error on an Algorithm for Downdating a Cholesky Factorization,” J. Inst. Math. Applic. 23, 203–213.

A.W. Bojanczyk, R.P. Brent, P. Van Dooren, and F.R. de Hoog (1987). “A Note on Downdating the Cholesky Factorization,” SIAM J. Sci. Stat. Comput. 8, 210–221.

C.-T. Pan (1993). “A Perturbation Analysis of the Problem of Downdating a Cholesky Factorization,” Lin. Alg. Applic. 183, 103–115.

L. Eld´en and H. Park (1994). “Perturbation Analysis for Block Downdating of a Cholesky Decomposition,” Numer. Math. 68, 457–468.

M.R. Osborne and L. Sun (1999). “A New Approach to Symmetric Rank-One Updating,” IMA J. Numer. Anal. 19, 497–507.

E.S. Quintana-Orti and R.A. Van Geijn (2008). “Updating an LU Factorization with Pivoting,” ACM Trans. Math. Softw. 35(2), Article 11.

Hyperbolic tranformations have been successfully used in a number of settings:

G.H. Golub (1969). “Matrix Decompositions and Statistical Computation,” in Statistical Computation, ed., R.C. Milton and J.A. Nelder, Academic Press, New York, pp. 365–397.

C.M. Rader and A.O. Steinhardt (1988). “Hyperbolic Householder Transforms,” SIAM J. Matrix Anal. Applic. 9, 269–290.

S.T. Alexander, C.T. Pan, and R.J. Plemmons (1988). “Analysis of a Recursive Least Squares Hyperbolic Rotation Algorithm for Signal Processing,” Lin. Alg. and Its Applic. 98, 3–40.   
G. Cybenko and M. Berry (1990). “Hyperbolic Householder Algorithms for Factoring Structured Matrices,” SIAM J. Matrix Anal. Applic. 11, 499–520.   
A.W. Bojanczyk, R. Onn, and A.O. Steinhardt (1993). “Existence of the Hyperbolic Singular Value Decomposition,” Lin. Alg. Applic. 185, 21–30.   
S. Chandrasekaran, M. Gu, and A.H. Sayad (1998). “A Stable and Efficient Algorithm for the Indefinite Linear Least Squares Problem,” SIAM J. Matrix Anal. Applic. 20, 354–362.   
A.J. Bojanczyk, N.J. Higham, and H. Patel (2003a). “Solving the Indefinite Least Squares Problem by Hyperbolic QR Factorization,” SIAM J. Matrix Anal. Applic. 24, 914–931.   
A. Bojanczyk, N.J. Higham, and H. Patel (2003b). “The Equality Constrained Indefinite Least Squares Problem: Theory and Algorithms,” BIT 43, 505–517.   
M. Stewart and P. Van Dooren (2006). “On the Factorization of Hyperbolic and Unitary Transformations into Rotations,” SIAM J. Matrix Anal. Applic. 27, 876–890.   
N.J. Higham (2003). “J-Orthogonal Matrices: Properties and Generation,” SIAM Review 45, 504–519.   
High-performance issues associated with QR updating are discussed in:   
B.C. Gunter and R.A. Van De Geijn (2005). “Parallel Out-of-Core Computation and Updating of the QR Factorization,” ACM Trans. Math. Softw. 31, 60–78.   
Updating and downdating the ULV and URV decompositions and related topics are covered in:   
C.H. Bischof and G.M. Shroff (1992). “On Updating Signal Subspaces,” IEEE Trans. Signal Proc. 40, 96–105.   
G.W. Stewart (1992). “An Updating Algorithm for Subspace Tracking,” IEEE Trans. Signal Proc. 40, 1535–1541.   
G.W. Stewart (1993). “Updating a Rank-Revealing ULV Decomposition,” SIAM J. Matrix Anal. Applic. 14, 494–499.   
G.W. Stewart (1994). “Updating URV Decompositions in Parallel,” Parallel Comp. 20, 151–172.   
H. Park and L. Eld´en (1995). “Downdating the Rank-Revealing URV Decomposition,” SIAM J. Matrix Anal. Applic. 16, 138–155.   
J.L. Barlow and H. Erbay (2009). “Modifiable Low-Rank Approximation of a Matrix,” Num. Lin. Alg. Applic. 16, 833–860.

Other interesting update-related topics include the updating of condition estimates, see:

W.R. Ferng, G.H. Golub, and R.J. Plemmons (1991). “Adaptive Lanczos Methods for Recursive Condition Estimation,” Numerical Algorithms 1, 1-20.

G. Shroff and C.H. Bischof (1992). “Adaptive Condition Estimation for Rank-One Updates of QR Factorizations,” SIAM J. Matrix Anal. Applic. 13, 1264–1278.

D.J. Pierce and R.J. Plemmons (1992). “Fast Adaptive Condition Estimation,” SIAM J. Matrix Anal. Applic. 13, 274–291.

and the updating of solutions to constrained least squares problems:

K. Schittkowski and J. Stoer (1979). “A Factorization Method for the Solution of Constrained Linear Least Squares Problems Allowing for Subsequent Data changes,” Numer. Math. 31, 431–463.

˚A. Bj¨orck (1984). “A General Updating Algorithm for Constrained Linear Least Squares Problems,” SIAM J. Sci. Stat. Comput. 5, 394–402.

Finally, we mention the following paper concerned with SVD updating:

M. Moonen, P. Van Dooren, and J. Vandewalle (1992). “A Singular Value Decomposition Updating Algorithm,” SIAM J. Matrix Anal. Applic. 13, 1015–1038.
