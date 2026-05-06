# 2.5 Subspace Metrics

If the object of a computation is to compute a matrix or a vector, then norms are useful for assessing the accuracy of the answer or for measuring progress during an iteration. If the object of a computation is to compute a subspace, then to make similar comments we need to be able to quantify the distance between two subspaces. Orthogonal projections are critical in this regard. After the elementary concepts are established we discuss the CS decomposition. This is an SVD-like decomposition that is handy when we have to compare a pair of subspaces.

# 2.5.1 Orthogonal Projections

Let $S \subseteq \mathbb { R } ^ { n }$ be a subspace. $P \in \mathbb { R } ^ { n \times n }$ is the orthogonal projection onto S if ran $( P ) = S$ , $P ^ { 2 } = P ,$ , and $P ^ { T } = P$ . From this definition it is easy to show that if $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , then $P x \in S$ and $( I - P ) x \in S ^ { \perp }$ .

If $P _ { 1 }$ and $P _ { 2 }$ are each orthogonal projections, then for any $z \in \mathbb { R } ^ { n }$ we have

$$
\| (P _ {1} - P _ {2}) z \| _ {2} ^ {2} = (P _ {1} z) ^ {T} (I - P _ {2}) z + (P _ {2} z) ^ {T} (I - P _ {1}) z.
$$

If $\mathsf { r a n } ( P _ { 1 } ) = \mathsf { r a n } ( P _ { 2 } ) = S$ , then the right-hand side of this expression is zero, showing that the orthogonal projection for a subspace is unique. If the columns of $V =$ $\left[ \left. v _ { 1 } \right| \cdots \right| \left. v _ { k } \right]$ are an orthonormal basis for a subspace S, then it is easy to show that ${ \cal P } = { \cal V } { \cal V } ^ { T }$ is the unique orthogonal projection onto S. Note that if $v \in \mathbb { R } ^ { n }$ , then $P = v v ^ { T } / v ^ { T } v$ is the orthogonal projection onto $S = \mathsf { s p a n } \{ v \}$ .

# 2.5.2 SVD-Related Projections

There are several important orthogonal projections associated with the singular value decomposition. Suppose $A = U \Sigma V ^ { T } \in \mathbb { R } ^ { m \times n }$ is the SVD of A and that $r = \mathsf { r a n k } ( A )$ . If we have the U and V partitionings

$$
U   =   \left[ \begin{array}{c c} U _ {r} & \tilde {U} _ {r} \\ r & m - r \end{array} \right]    , \qquad V   =   \left[ \begin{array}{c c} V _ {r} & \tilde {V} _ {r} \\ r & n - r \end{array} \right]    ,
$$

then

$$
V _ {r} V _ {r} ^ {T} = \text { projection   on   to   } \text { null } (A) ^ {\perp} = \text { ran } (A ^ {T}),
$$

$$
\tilde {V} _ {r} \tilde {V} _ {r} ^ {T} = \text { projection   on   to   } \mathsf {n u l l} (A),
$$

$$
U _ {r} U _ {r} ^ {T} = \text { projection   on   to   } \operatorname{ran} (A),
$$

$$
\tilde {U} _ {r} \tilde {U} _ {r} ^ {T} = \text { projection   on   to } \operatorname{ran} (A) ^ {\perp} = \operatorname{null} (A ^ {T}).
$$

# 2.5.3 Distance Between Subspaces

The one-to-one correspondence between subspaces and orthogonal projections enables us to devise a notion of distance between subspaces. Suppose $S _ { 1 }$ and $S _ { 2 }$ are subspaces of $\mathbb { R } ^ { n }$ and that dim $( S _ { 1 } ) = \mathsf { d i m } ( S _ { 2 } )$ . We define the distance between these two spaces by

$$
\mathrm{dist} (S _ {1}, S _ {2}) = \| P _ {1} - P _ {2} \| _ {2} \tag {2.5.1}
$$

where $P _ { i }$ is the orthogonal projection onto $S _ { i }$ . The distance between a pair of subspaces can be characterized in terms of the blocks of a certain orthogonal matrix.

Theorem 2.5.1. Suppose

$$
W = \left[ \begin{array}{c c} W _ {1} & W _ {2} \\ k & n - k \end{array} \right], Z = \left[ \begin{array}{c c} Z _ {1} & Z _ {2} \\ k & n - k \end{array} \right],
$$

are $n { - } b y { - } n$ orthogonal matrices. $I f S _ { 1 } = \mathsf { r a n } ( W _ { 1 } )$ and $S _ { 2 } = \mathsf { r a n } ( Z _ { 1 } )$ , then

$$
\operatorname{dist} \left(S _ {1}, S _ {2}\right) = \left\| W _ {1} ^ {T} Z _ {2} \right\| _ {2} = \left\| Z _ {1} ^ {T} W _ {2} \right\| _ {2}.
$$

Proof. We first observe that

$$
\begin{array}{l} \operatorname{dist} \left(S _ {1}, S _ {2}\right) = \left\| W _ {1} W _ {1} ^ {T} - Z _ {1} Z _ {1} ^ {T} \right\| _ {2} \\ = \left\| W ^ {T} \left(W _ {1} W _ {1} ^ {T} - Z _ {1} Z _ {1} ^ {T}\right) Z \right\| _ {2} \\ = \left\| \left[ \begin{array}{c c} 0 & W _ {1} ^ {T} Z _ {2} \\ - W _ {2} ^ {T} Z _ {1} & 0 \end{array} \right] \right\| _ {2}. \\ \end{array}
$$

Note that the matrices $W _ { 2 } ^ { T } Z _ { 1 }$ and $W _ { 1 } ^ { T } Z _ { 2 }$ are submatrices of the orthogonal matrix

$$
Q = \left[ \begin{array}{l l} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right] \equiv \left[ \begin{array}{l l} W _ {1} ^ {T} Z _ {1} & W _ {1} ^ {T} Z _ {2} \\ W _ {2} ^ {T} Z _ {1} & W _ {2} ^ {T} Z _ {2} \end{array} \right] = W ^ {T} Z. \tag {2.5.2}
$$

Our goal is to show that $\parallel Q _ { 2 1 } \parallel _ { 2 } = \parallel Q _ { 1 2 } \parallel _ { 2 }$ . Since Q is orthogonal it follows from

$$
Q \left[ \begin{array}{l} x \\ 0 \end{array} \right] = \left[ \begin{array}{l} Q _ {1 1} x \\ Q _ {2 1} x \end{array} \right]
$$

that $1 = \parallel Q _ { 1 1 } x \parallel _ { 2 } ^ { 2 } + \parallel Q _ { 2 1 } x \parallel _ { 2 } ^ { 2 }$ for all unit 2-norm $\boldsymbol { x } \in \mathbb { R } ^ { k }$ . Thus,

$$
\| Q _ {2 1} \| _ {2} ^ {2} = \max _ {\| x \| _ {2} = 1} \| Q _ {2 1} x \| _ {2} ^ {2} = 1 - \min _ {\| x \| _ {2} = 1} \| Q _ {1 1} x \| _ {2} ^ {2} = 1 - \sigma_ {\min} (Q _ {1 1}) ^ {2}.
$$

Analogously, by working with $Q ^ { T }$ (which is also orthogonal) it is possible to show that

$$
\parallel Q _ {1 2} ^ {T} \parallel_ {2} ^ {2} = 1 - \sigma_ {\min} (Q _ {1 1} ^ {T}) ^ {2},
$$

and therefore

$$
\parallel Q _ {1 2} \parallel_ {2} ^ {2} = 1 - \sigma_ {\min} (Q _ {1 1}) ^ {2}.
$$

Thus,  Q21 2 =  Q12 2.

Note that if $S _ { 1 }$ and $S _ { 2 }$ are subspaces in $\mathbb { R } ^ { n }$ with the same dimension, then

$$
0 \leq \operatorname{dist} \left(S _ {1}, S _ {2}\right) \leq 1.
$$

It is easy to show that

$$
\begin{array}{l} \operatorname{dist} \left(S _ {1}, S _ {2}\right) = 0 \Rightarrow S _ {1} = S _ {2}, \\ \operatorname{dist} \left(S _ {1}, S _ {2}\right) = 1 \Rightarrow S _ {1} \cap S _ {2} ^ {\perp} \neq \{0 \}. \\ \end{array}
$$

A more refined analysis of the blocks of the matrix Q in (2.5.2) sheds light on the difference between a pair of subspaces. A special, SVD-like decomposition for orthogonal matrices is required.

# 2.5.4 The CS Decomposition

The blocks of an orthogonal matrix partitioned into 2-by-2 form have highly related SVDs. This is the gist of the CS decomposition. We prove a very useful special case first.

Theorem 2.5.2 (The CS Decomposition (Thin Version)). Consider the matrix

$$
Q = \left[ \begin{array}{l} Q _ {1} \\ Q _ {2} \end{array} \right], \qquad \qquad Q _ {1} \in \mathbb {R} ^ {m _ {1} \times n _ {1}},   Q _ {2} \in \mathbb {R} ^ {m _ {2} \times n _ {1}},
$$

where $m _ { 1 } \geq n _ { 1 }$ and $m _ { 2 } \geq n _ { 1 }$ . If the columns of Q are orthonormal, then there exist orthogonal matrices $U _ { 1 } \in \mathbb { R } ^ { m _ { 1 } \times m _ { 1 } } , U _ { 2 } \in \mathbb { R } ^ { m _ { 2 } \times m _ { 2 } }$ , and $V _ { 1 } \in \mathbb { R } ^ { n _ { 1 } \times n _ { 1 } }$ such that

$$
\left[ \begin{array}{c c} U _ {1} & 0 \\ 0 & U _ {2} \end{array} \right] ^ {T} \left[ \begin{array}{c} Q _ {1} \\ Q _ {2} \end{array} \right] V _ {1} = \left[ \begin{array}{c} C \\ S \end{array} \right]
$$

where

$$
C _ {0} = \operatorname{diag} \left(\cos \left(\theta_ {1}\right), \dots , \cos \left(\theta_ {n _ {1}}\right)\right) \in \mathbb {R} ^ {m _ {1} \times n _ {1}},
$$

$$
S _ {0} = \operatorname{diag} \left(\sin \left(\theta_ {1}\right), \dots , \sin \left(\theta_ {n _ {1}}\right)\right) \in \mathbb {R} ^ {m _ {2} \times n _ {1}},
$$

and

$$
0 \leq \theta_ {1} \leq \theta_ {2} \leq \dots \leq \theta_ {n _ {1}} \leq \frac {\pi}{2}.
$$

Proof. Since $\| Q _ { 1 } \| _ { 2 } \leq \| Q \| _ { 2 } = 1$ , the singular values of $Q _ { 1 }$ are all in the interval [0, 1]. Let

$$
U _ {1} ^ {T} Q _ {1} V _ {1} = C _ {0} = \mathrm{diag} (c _ {1}, \ldots , c _ {n _ {1}}) = \left[ \begin{array}{c c} I _ {t} & 0 \\ 0 & \Sigma \end{array} \right] _ {m _ {1} - t} ^ {t}
$$

be the SVD of $Q _ { 1 }$ where we assume

$$
1 = c _ {1} = \dots = c _ {t} > c _ {t + 1} \geq \dots \geq c _ {n _ {1}} \geq 0.
$$

To complete the proof of the theorem we must construct the orthogonal matrix $U _ { 2 }$ . If

$$
Q _ {2} V _ {1} = \left[ \begin{array}{c c} W _ {1} & W _ {2} \\ t & n _ {1} - t \end{array} \right] \quad ,
$$

then

$$
\left[ \begin{array}{c c} U _ {1} & 0 \\ 0 & I _ {m _ {2}} \end{array} \right] ^ {T} \left[ \begin{array}{c} Q _ {1} \\ Q _ {2} \end{array} \right] V _ {1} = \left[ \begin{array}{c c} I _ {t} & 0 \\ 0 & \Sigma \\ W _ {1} & W _ {2} \end{array} \right].
$$

Since the columns of this matrix have unit 2-norm, $W _ { 1 } = 0$ . The columns of $W _ { 2 }$ are nonzero and mutually orthogonal because

$$
W _ {2} ^ {T} W _ {2} = I _ {n _ {1} - t} - \Sigma^ {T} \Sigma \equiv \mathrm{diag} (1 - c _ {t + 1} ^ {2}, \ldots , 1 - c _ {n _ {1}} ^ {2})
$$

is nonsingular. If $s _ { k } = \sqrt { 1 - c _ { k } ^ { 2 } }$ for $k = 1 { : } n _ { 1 }$ , then the columns of

$$
Z = W _ {2} \operatorname{diag} (1 / s _ {t + 1}, \dots , 1 / s _ {n})
$$

are orthonormal. By Theorem 2.1.1 there exists an orthogonal matrix $U _ { 2 } \in \mathbb { R } ^ { m _ { 2 } \times m _ { 2 } }$ with $U _ { 2 } ( : , t + 1 { : } n _ { 1 } ) = Z$ . It is easy to verify that

$$
U _ {2} ^ {T} Q _ {2} V _ {1} = \mathrm{diag} (s _ {1}, \ldots , s _ {n _ {1}}) \equiv S _ {0}.
$$

Since $c _ { k } ^ { 2 } + s _ { k } ^ { 2 } = 1$ for $k = 1 { : } n _ { 1 }$ , it follows that these quantities are the required cosines and sines. □

By using the same techniques it is possible to prove the following, more general version of the decomposition:

Theorem 2.5.3 (CS Decomposition). Suppose

$$
Q = = \left[ \begin{array}{c c} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right] _ {m _ {2}} ^ {m _ {1}}
$$

is a square orthogonal matrix and that $m _ { 1 } \geq n _ { 1 }$ and $m _ { 1 } \geq m _ { 2 }$ . Define the nonnegative integers p and q by $\textit { p } = \operatorname* { m a x } \{ 0 , n _ { 1 } - m _ { 2 } \}$ and q = max $\{ 0 , m _ { 2 } - n _ { 1 } \}$ . There exist orthogonal $U _ { 1 } \in \mathbb { R } ^ { m _ { 1 } \times m _ { 1 } }$ , $U _ { 2 } \in \mathbb { R } ^ { m _ { 2 } \times m _ { 2 } }$ , $V _ { 1 } \in \mathbb { R } ^ { n _ { 1 } \times n _ { 1 } }$ , and $V _ { 2 } \in \mathbb { R } ^ { n _ { 2 } \times n _ { 2 } }$ such that if

$$
U = \left[ \begin{array}{c c} U _ {1} & 0 \\ \hline 0 & U _ {2} \end{array} \right] \qquad \text { and } \qquad V = \left[ \begin{array}{c c} V _ {1} & 0 \\ \hline 0 & V _ {2} \end{array} \right],
$$

then

$$
U ^ {T} Q V = \left[ \begin{array}{c c c c c} I & 0 & 0 & 0 & 0 \\ 0 & C & S & 0 & 0 \\ 0 & 0 & 0 & 0 & I \\ \hline 0 & S & - C & 0 & 0 \\ 0 & 0 & 0 & I & 0 \\ \hline p & n _ {1} - p & n _ {1} - p & q & m _ {1} - n _ {1} \end{array} \right] \begin{array}{l} p \\ n _ {1} - p \\ m _ {1} - n _ {1} \\ n _ {1} - p \\ q \end{array}
$$

where

$$
C = \operatorname{diag} \left(\cos \left(\theta_ {p + 1}\right), \dots , \cos \left(\theta_ {n _ {1}}\right)\right) = \operatorname{diag} \left(c _ {p + 1}, \dots , c _ {n _ {1}}\right),
$$

$$
S = \operatorname{diag} \left(\sin \left(\theta_ {p + 1}\right), \dots , \sin \left(\theta_ {n _ {1}}\right)\right) = \operatorname{diag} \left(s _ {p + 1}, \dots , s _ {n _ {1}}\right),
$$

and $0 \le \theta _ { p + 1 } \le \cdots \le \theta _ { n _ { 1 } } \le \pi / 2$ .

Proof. See Paige and Saunders (1981) for details.

We made the assumptions $m _ { 1 } \geq n _ { 1 }$ and $m _ { 1 } \geq m _ { 2 }$ for clarity. Through permutation and transposition, any 2-by-2 block orthogonal matrix can be put into the form required by the theorem. Note that the blocks in the transformed Q, i.e., the $U _ { i } ^ { T } Q _ { i j } V _ { j }$ , are diagonal-like but not necessarily diagonal. Indeed, as we have presented it, the CS decomposition gives us four unnormalized SVDs. If $Q _ { 2 1 }$ has more rows than columns, then $p = 0$ and the reduction looks like this (for example):

$$
U ^ {T} Q V = \left[ \begin{array}{c c c c c c c} c _ {1} & 0 & s _ {1} & 0 & 0 & 0 & 0 \\ 0 & c _ {2} & 0 & s _ {2} & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 1 \\ \hline s _ {1} & 0 & - c _ {1} & 0 & 0 & 0 & 0 \\ 0 & s _ {2} & 0 & - c _ {2} & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 1 & 0 & 0 \end{array} \right].
$$

On the other hand, if $Q _ { 2 1 }$ has more columns than rows, then $q = 0$ and the decomposition has the form

$$
U ^ {T} Q V = \left[ \begin{array}{c c c c c} 1 & 0 & 0 & 0 & 0 \\ 0 & c _ {2} & 0 & s _ {2} & 0 \\ 0 & 0 & c _ {3} & 0 & s _ {3} \\ \hline 0 & s _ {2} & 0 & - c _ {2} & 0 \\ 0 & 0 & s _ {3} & 0 & - c _ {3} \end{array} \right].
$$

Regardless of the partitioning, the essential message of the CS decomposition is that the SVDs of the Q-blocks are highly related.

# Problems

P2.5.1 Show that if P is an orthogonal projection, then $Q = I - 2 P$ is orthogonal.

P2.5.2 What are the singular values of an orthogonal projection?

P2.5.3 Suppose $S _ { 1 } = \mathsf { s p a n } \{ x \}$ and $S _ { 2 } = \mathsf { s p a n } \{ y \}$ , where x and y are unit 2-norm vectors in $\mathbb { R } ^ { 2 }$ . Working only with the definition of dist $( \cdot , \cdot )$ , show that dis $\langle S _ { 1 } , S _ { 2 } \rangle = \sqrt { 1 - ( x ^ { T } y ) ^ { 2 } }$ , verifying that the distance between $S _ { 1 }$ and $S _ { 2 }$ equals the sine of the angle between x and y.

P2.5.4 Refer to §1.3.10. Show that if $Q \in \mathbb { R } ^ { 2 n \times 2 n }$ is orthogonal and symplectic, then Q has the form

$$
Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ - Q _ {2} & Q _ {1} \end{array} \right], \qquad Q _ {1}, Q _ {2} \in \mathbb {R} ^ {n \times n}.
$$

P2.5.5 Suppose $P \in \mathbb { R } ^ { n \times n }$ and $P ^ { 2 } = P$ . Show that $\| \ P \| _ { 2 } > 1$ if null(P ) is not a subspace of $\mathsf { r a n } ( A ) ^ { \perp }$ . Such a matrix is called an oblique projector. See Stewart (2011).

# Notes and References for 2.5

The computation of the CS decomposition is discussed in 8.7.6. For a discussion of its analytical properties, see:   
C. Davis and W. Kahan (1970). “The Rotation of Eigenvectors by a Perturbation III,” SIAM J. Numer. Anal. 7, 1–46.   
G.W. Stewart (1977). “On the Perturbation of Pseudo-Inverses, Projections and Linear Least Squares Problems,” SIAM Review 19, 634–662.   
C.C. Paige and M. Saunders (1981). “Toward a Generalized Singular Value Decomposition,” SIAM J. Numer. Anal. 18, 398–405.   
C.C. Paige and M. Wei (1994). “History and Generality of the CS Decomposition,” Lin. Alg. Applic. 208/209, 303–326.   
A detailed numerical discussion of oblique projectors (P2.5.5) is given in:   
G.W. Stewart (2011). “On the Numerical Analysis of Oblique Projectors,” SIAM J. Matrix Anal. Applic. 32, 309–348.
