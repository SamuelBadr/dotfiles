# 7.6 Invariant Subspace Computations

Several important invariant subspace problems can be solved once the real Schur decomposition $Q ^ { T } A Q = T$ has been computed. In this section we discuss how to

• compute the eigenvectors associated with some subset of λ(A),   
• compute an orthonormal basis for a given invariant subspace,   
• block-diagonalize A using well-conditioned similarity transformations,   
• compute a basis of eigenvectors regardless of their condition, and   
• compute an approximate Jordan canonical form of A.

Eigenvector/invariant subspace computation for sparse matrices is discussed in §7.3.1 and §7.3.2 as well as portions of Chapters 8 and 10.

# 7.6.1 Selected Eigenvectors via Inverse Iteration

Let $q ^ { ( 0 ) } \in \mathbb { R } ^ { n }$ be a given unit 2-norm vector and assume that $A - \mu I \in \mathbb { R } ^ { n \times n }$ is nonsingular. The following is referred to as inverse iteration:

for k = 1, 2, . . .

$$
\text { Solve } (A - \mu I) z ^ {(k)} = q ^ {(k - 1)}.
$$

$$
q ^ {(k)} = z ^ {(k)} / \| z ^ {(k)} \| _ {2} \tag {7.6.1}
$$

$$
\lambda^ {(k)} = q ^ {(k) ^ {T}} A q ^ {(k)}
$$

end

Inverse iteration is just the power method applied to $( A - \mu I ) ^ { - 1 }$ .

To analyze the behavior of (7.6.1), assume that A has a basis of eigenvectors $\{ x _ { 1 } , \ldots , x _ { n } \}$ and that $A x _ { i } = \lambda _ { i } x _ { i }$ for i = 1:n. If

$$
q ^ {(0)} = \sum_ {i = 1} ^ {n} \beta_ {i} x _ {i}
$$

then $q ^ { ( k ) }$ is a unit vector in the direction of

$$
(A - \mu I) ^ {- k} q ^ {(0)} = \sum_ {i = 1} ^ {n} \frac {\beta_ {i}}{(\lambda_ {i} - \mu) ^ {k}} x _ {i}.
$$

Clearly, if $\mu$ is much closer to an eigenvalue $\lambda _ { j }$ than to the other eigenvalues, then $q ^ { ( k ) }$ is rich in the direction of $x _ { j }$ provided $\beta _ { j } \neq 0$ .

A sample stopping criterion for (7.6.1) might be to quit as soon as the residual

$$
r ^ {(k)} = (A - \mu I) q ^ {(k)}
$$

satisfies

$$
\| r ^ {(k)} \| _ {\infty} \leq c \mathbf {u} \| A \| _ {\infty} \tag {7.6.2}
$$

where c is a constant of order unity. Since

$$
(A + E _ {k}) q ^ {(k)} = \mu q ^ {(k)}
$$

with $E _ { k } = - r ^ { ( k ) } q ^ { ( k ) ^ { T } }$ , it follows that (7.6.2) forces $\mu$ and $q ^ { ( k ) }$ to be an exact eigenpair for a nearby matrix.

Inverse iteration can be used in conjunction with Hessenberg reduction and the QR algorithm as follows:

Step 1. Compute the Hessenberg decomposition $U _ { 0 } ^ { T } A U _ { 0 } = H .$ .

Step 2. Apply the double-implicit-shift Francis iteration to H without accumulating transformations.

Step 3. For each computed eigenvalue λ whose corresponding eigenvector x is sought, apply (7.6.1) with A = H and $\mu = \lambda$ to produce a vector z such that $H z \approx \mu z$ .

Step 4. Set $x = U _ { 0 } z$

Inverse iteration with H is very economical because we do not have to accumulate transformations during the double Francis iteration. Moreover, we can factor matrices of the form $H - \lambda I$ in $O ( n ^ { 2 } )$ flops, and (3) only one iteration is typically required to produce an adequate approximate eigenvector.

This last point is perhaps the most interesting aspect of inverse iteration and requires some justification since λ can be comparatively inaccurate if it is ill-conditioned. Assume for simplicity that λ is real and let

$$
H - \lambda I = \sum_ {i = 1} ^ {n} \sigma_ {i} u _ {i} v _ {i} ^ {T} = U \Sigma V ^ {T}
$$

be the SVD of $H - \lambda I$ . From what we said about the roundoff properties of the QR algorithm in $\ S 7 . 5 . 6 \AA$ , there exists a matrix $E \in \mathbb { R } ^ { n \times n }$ such that $H + E - \lambda I$ is singular and $\parallel E \parallel _ { 2 } \approx \mathbf { u } \parallel H \parallel _ { 2 }$ . It follows that $\sigma _ { n } \approx \mathbf { u } \sigma _ { 1 }$ and

$$
\| (H - \hat {\lambda} I) v _ {n} \| _ {2} \approx \mathbf {u} \sigma_ {1},
$$

i.e., $v _ { n }$ is a good approximate eigenvector. Clearly if the starting vector $q ^ { ( 0 ) }$ has the expansion

$$
q ^ {(0)} = \sum_ {i = 1} ^ {n} \gamma_ {i} u _ {i}
$$

then

$$
z ^ {(1)} = \sum_ {i = 1} ^ {n} \frac {\gamma_ {i}}{\sigma_ {i}} v _ {i}
$$

is “rich” in the direction $v _ { n }$ . Note that if $s ( \lambda ) \approx | u _ { n } ^ { T } v _ { n } |$ is small, then $z ^ { ( 1 ) }$ is rather deficient in the direction $u _ { n }$ . This explains (heuristically) why another step of inverse iteration is not likely to produce an improved eigenvector approximate, especially if λ is ill-conditioned. For more details, see Peters and Wilkinson (1979).

# 7.6.2 Ordering Eigenvalues in the Real Schur Form

Recall that the real Schur decomposition provides information about invariant subspaces. If

$$
Q ^ {T} A Q = T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {q} ^ {p}
$$

and

$$
\lambda (T _ {1 1}) \cap \lambda (T _ {2 2}) = \emptyset ,
$$

then the first $p$ columns of $Q$ span the unique invariant subspace associated with $\lambda ( T _ { 1 1 } )$ . (See $\ S 7 . 1 . 4 . )$ Unfortunately, the Francis iteration supplies us with a real Schur decomposition $Q _ { F } ^ { T } A Q _ { F } = T _ { F }$ in which the eigenvalues appear somewhat randomly along the diagonal of $T _ { F }$ . This poses a problem if we want an orthonormal basis for an invariant subspace whose associated eigenvalues are not at the top of $T _ { F } \mathrm { \ ' } _ { \mathrm { s } }$ diagonal. Clearly, we need a method for computing an orthogonal matrix $Q _ { D }$ such that $Q _ { D } ^ { T } T _ { F } Q _ { D }$ is upper quasi-triangular with appropriate eigenvalue ordering.

A look at the 2-by-2 case suggests how this can be accomplished. Suppose

$$
Q _ {F} ^ {T} A Q _ {F} = T _ {F} = \left[ \begin{array}{c c} \lambda_ {1} & t _ {1 2} \\ 0 & \lambda_ {2} \end{array} \right], \qquad \lambda_ {1} \neq \lambda_ {2}
$$

and that we wish to reverse the order of the eigenvalues. Note that

$$
T _ {F} x = \lambda_ {2} x
$$

where

$$
x = \left[ \begin{array}{c} t _ {1 2} \\ \lambda_ {2} - \lambda_ {1} \end{array} \right].
$$

Let $Q _ { D }$ be a Givens rotation such that the second component of $Q _ { D } ^ { T } x$ is zero. If

$$
Q = Q _ {F} Q _ {D},
$$

then

$$
(Q ^ {T} A Q) e _ {1} = Q _ {D} ^ {T} T _ {F} (Q _ {D} e _ {1}) = \lambda_ {2} Q _ {D} ^ {T} (Q _ {D} e _ {1}) = \lambda_ {2} e _ {1}.
$$

The matrices A and $Q ^ { T } A Q$ have the same Frobenius norm and so it follows that the latter must have the following form:

$$
Q ^ {T} A Q = \left[ \begin{array}{c c} \lambda_ {2} & \pm t _ {1 2} \\ 0 & \lambda_ {1} \end{array} \right].
$$

The swapping gets a little more complicated if $T$ has 2-by-2 blocks along its diagonal. See Ruhe (1970) and Stewart (1976) for details.

By systematically interchanging adjacent pairs of eigenvalues (or 2-by-2 blocks), we can move any subset of $\lambda ( A )$ to the top of $T \mathrm { { s } }$ diagonal. Here is the overall procedure for the case when there are no 2-by-2 bumps:

Algorithm 7.6.1 Given an orthogonal matrix $Q \in \mathbb { R } ^ { n \times n }$ , an upper triangular matrix ${ \bar { T } } { \bar { = } } Q ^ { T } A Q$ , and a subset $\Delta = \{ \lambda _ { 1 } , \ldots , \lambda _ { p } \}$ of $\lambda ( A )$ , the following algorithm computes an orthogonal matrix $Q _ { D }$ such that $Q _ { D } ^ { T } T \dot { Q } _ { D } = S$ is upper triangular and $\{ s _ { 1 1 } , \dotsc , s _ { p p } \}$ $= \Delta$ . The matrices $Q$ and $T$ are overwritten by $Q Q _ { D }$ and $S ,$ respectively.

while $\{ t _ { 1 1 } , \hdots , t _ { p p } \} \neq \Delta$

for $k = 1 { : } n - 1$

$\mathbf { i f } \ t _ { k k } \notin \Delta \ \mathrm { a n d } \ t _ { k + 1 , k + 1 } \in \Delta$

$$
[ c, s ] = \text { g   i   v   e   n   s } (T (k, k + 1), T (k + 1, k + 1) - T (k, k))
$$

$$
T (k: k + 1, k: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} T (k: k + 1, k: n)
$$

$$
T (1: k + 1, k: k + 1) = T (1: k + 1, k: k + 1) \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right]
$$

$$
Q (1: n, k: k + 1) = Q (1: n, k: k + 1) \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right]
$$

end end end

This algorithm requires $k ( 1 2 n )$ flops, where k is the total number of required swaps. The integer k is never greater than $( n - p ) p$ .

Computation of invariant subspaces by manipulating the real Schur decomposition is extremely stable. If $\hat { Q } = [ \hat { q } _ { 1 } | \cdot \cdot \cdot | \hat { q } _ { n } ]$ denotes the computed orthogonal matrix $Q .$ then $\Vert \hat { Q } ^ { T } \hat { Q } - I \Vert _ { 2 } \approx$ u and there exists a matrix E satisfying  $E \parallel _ { 2 } \approx \mathbf { u } \parallel A$ 2 such that $( A + E ) \hat { q } _ { i } \ \in \mathsf { s p a n } \{ \hat { q } _ { 1 } , \dots , \hat { q } _ { p } \}$ for $i = 1 { : } p$ .

# 7.6.3 Block Diagonalization

Let

$$
T = \left[ \begin{array}{c c c c} T _ {1 1} & T _ {1 2} & \dots & T _ {1 q} \\ 0 & T _ {2 2} & \dots & T _ {2 q} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & T _ {q q} \end{array} \right] \begin{array}{c} n _ {1} \\ n _ {2} \\ n _ {q} \end{array} \tag {7.6.3}
$$

be a partitioning of some real Schur canonical form $Q ^ { T } A Q \ = \ T \in \mathbb { R } ^ { n \times n }$ such that $\lambda ( T _ { 1 1 } ) , \dots , \lambda ( T _ { q q } )$ are disjoint. By Theorem 7.1.6 there exists a matrix Y such that

$$
Y ^ {- 1} T Y = \mathrm{diag} (T _ {1 1}, \ldots , T _ {q q}).
$$

A practical procedure for determining $Y$ is now given together with an analysis of $Y \mathrm { { } s }$ sensitivity as a function of the above partitioning.

Partition $I _ { n } = [ E _ { 1 } | \cdot \cdot \cdot | E _ { q } ]$ conformably with T and define the matrix $Y _ { i j } \in \mathbb { R } ^ { n \times n }$ as follows:

$$
Y _ {i j} = I _ {n} + E _ {i} Z _ {i j} E _ {j} ^ {T}, \quad i <   j, Z _ {i j} \in \mathbb {R} ^ {n _ {i} \times n _ {j}}.
$$

In other words, $Y _ { i j }$ looks just like the identity except that $Z _ { i j }$ occupies the $( i , j )$ block position. It follows that if $Y _ { i j } ^ { - 1 } T Y _ { i j } \ = \ \bar { T } = ( \bar { T } _ { i j } )$ , then $T$ and $\bar { T }$ are identical except that

$$
\begin{array}{l} \bar {T} _ {i j} = T _ {i i} Z _ {i j} - Z _ {i j} T _ {j j} + T _ {i j}, \\ \bar {T} _ {i k} = T _ {i k} - Z _ {i j} T _ {j k}, \quad (k = j + 1: q), \\ \bar {T} _ {k j} = T _ {k i} Z _ {i j} + T _ {k j}, \quad (k = 1: i - 1). \\ \end{array}
$$

Thus, $T _ { i j }$ can be zeroed provided we have an algorithm for solving the Sylvester equation

$$
F Z - Z G = C \tag {7.6.4}
$$

where $F \in \mathbb { R } ^ { p \times p }$ and $G \in \mathbb { R } ^ { r \times r }$ are given upper quasi-triangular matrices and $C \in \mathbb { R } ^ { p \times r }$

Bartels and Stewart (1972) have devised a method for doing this. Let $C =$ $\left[ c _ { 1 } \mid \cdots \mid c _ { r } \right]$ and $Z = \left[ \ : z _ { 1 } \ : | \cdots | \ : z _ { r } \ : \right]$ be column partitionings. If $g _ { k + 1 , k } = 0$ , then by comparing columns in (7.6.4) we find

$$
F z _ {k} - \sum_ {i = 1} ^ {k} g _ {i k} z _ {i} = c _ {k}.
$$

Thus, once we know $z _ { 1 } , \ldots , z _ { k - 1 }$ , then we can solve the quasi-triangular system

$$
(F - g _ {k k} I) z _ {k} = c _ {k} + \sum_ {i = 1} ^ {k - 1} g _ {i k} z _ {i}
$$

for $z _ { k }$ . If $g _ { k + 1 , k } \neq 0$ , then $z _ { k }$ and $z _ { k + 1 }$ can be simultaneously found by solving the $2 p { \mathrm { - } } \mathrm { b y - } 2 p$ system

$$
\left[ \begin{array}{c c} F - g _ {k k} I & - g _ {m k} I \\ - g _ {k m} I & F - g _ {m m} I \end{array} \right] \left[ \begin{array}{l} z _ {k} \\ z _ {m} \end{array} \right] = \left[ \begin{array}{l} c _ {k} \\ c _ {m} \end{array} \right] + \sum_ {i = 1} ^ {k - 1} \left[ \begin{array}{l} g _ {i k} z _ {i} \\ g _ {i m} z _ {i} \end{array} \right] \tag {7.6.5}
$$

where $m = k + 1$ . By reordering the equations according to the perfect shuffle permutation $( 1 , p + 1 , 2 , p + 2 , \ldots , p , 2 p )$ , a banded system is obtained that can be solved in $O ( p ^ { 2 } )$ flops. The details may be found in Bartels and Stewart (1972). Here is the overall process for the case when $F$ and $G$ are each triangular.

Algorithm 7.6.2 (Bartels-Stewart Algorithm) Given $C \in \mathbb { R } ^ { p \times r }$ and upper triangular matrices $\boldsymbol { F } \in \mathbb { R } ^ { p \times p }$ and $G \in \mathbb { R } ^ { r \times r }$ that satisfy $\lambda ( F ) \cap \lambda ( G ) = \emptyset$ , the following algorithm overwrites $C$ with the solution to the equation $F Z - Z G = C$ .

$$
\begin{array}{l} \text { for } k = 1: r \\ C (1: p, k) = C (1: p, k) + C (1: p, 1: k - 1) \cdot G (1: k - 1, k) \\ \text { Solve } (F - G (k, k) I) z = C (1: p, k) \text { for } z. \\ C (1: p, k) = z \\ \end{array}
$$

end

This algorithm requires $p r ( p + r )$ flops. By zeroing the superdiagonal blocks in $T$ in the appropriate order, the entire matrix can be reduced to block diagonal form.

Algorithm 7.6.3 Given an orthogonal matrix $Q \in \mathbb { R } ^ { n \times n }$ , an upper quasi-triangular matrix $T = Q ^ { T } A Q$ , and the partitioning (7.6.3), the following algorithm overwrites Q with QY where $Y ^ { - 1 } T Y = \mathrm { d i a g } ( T _ { 1 1 } , \dots , T _ { q q } )$ .

for j = 2:q

for i = 1:j − 1

Solve $T _ { i i } Z - Z T _ { j j } = - T _ { i j }$ for Z using the Bartels-Stewart algorithm.

for k = j + 1:q

$$
T _ {i k} = T _ {i k} - Z T _ {j k}
$$

end

for k = 1:q

$$
Q _ {k j} = Q _ {k i} Z + Q _ {k j}
$$

end

end end

The number of flops required by this algorithm is a complicated function of the block sizes in (7.6.3).

The choice of the real Schur form T and its partitioning in (7.6.3) determines the sensitivity of the Sylvester equations that must be solved in Algorithm 7.6.3. This in turn affects the condition of the matrix Y and the overall usefulness of the block diagonalization. The reason for these dependencies is that the relative error of the computed solution $\hat { Z }$ to

$$
T _ {i i} Z - Z T _ {j j} = - T _ {i j} \tag {7.6.6}
$$

satisfies

$$
\frac {\parallel \hat {Z} - Z \parallel_ {F}}{\parallel Z \parallel_ {F}} \approx \mathbf {u} \frac {\parallel T \parallel_ {F}}{\mathsf {s e p} (T _ {i i} , T _ {j j})}.
$$

For details, see Golub, Nash, and Van Loan (1979). Since

$$
\mathsf{sep}(T_{ii},T_{jj}) = \min_{X\neq 0}\frac{\parallel T_{ii}X - XT_{jj}\parallel_{F}}{\parallel X\parallel_{F}}\leq \min_{\substack{\lambda \in \lambda (T_{ii})\\ \mu \in \lambda (T_{jj})}}|\lambda -\mu |
$$

there can be a substantial loss of accuracy whenever the subsets $\lambda ( T _ { i i } )$ are insufficiently separated. Moreover, if Z satisfies (7.6.6) then

$$
\| Z \| _ {F} \leq \frac {\| T _ {i j} \| _ {F}}{\mathsf {s e p} (T _ {i i} , T _ {j j})}.
$$

Thus, large norm solutions can be expected if $\mathsf { s e p } ( T _ { i i } , T _ { j j } )$ is small. This tends to make the matrix Y in Algorithm 7.6.3 ill-conditioned since it is the product of the matrices

$$
Y _ {i j} = \left[ \begin{array}{c c} I _ {n _ {i}} & Z \\ 0 & I _ {n _ {j}} \end{array} \right].
$$

Note that $\kappa _ { F } ( Y _ { i j } ) = n _ { i } ^ { 2 } + n _ { j } ^ { 2 } + \left\| Z \right\| _ { F } ^ { 2 } .$ .

Confronted with these difficulties, Bavely and Stewart (1979) develop an algorithm for block diagonalizing that dynamically determines the eigenvalue ordering and partitioning in (7.6.3) so that all the Z matrices in Algorithm 7.6.3 are bounded in norm by some user-supplied tolerance. Their research suggests that the condition of Y can be controlled by controlling the condition of the $Y _ { i j }$ .

# 7.6.4 Eigenvector Bases

If the blocks in the partitioning (7.6.3) are all 1-by-1, then Algorithm 7.6.3 produces a basis of eigenvectors. As with the method of inverse iteration, the computed eigenvalueeigenvector pairs are exact for some “nearby” matrix. A widely followed rule of thumb for deciding upon a suitable eigenvector method is to use inverse iteration whenever fewer than 25% of the eigenvectors are desired.

We point out, however, that the real Schur form can be used to determine selected eigenvectors. Suppose

$$
Q ^ {T} A Q = \left[ \begin{array}{c c c} T _ {1 1} & u & T _ {1 3} \\ 0 & \lambda & v ^ {T} \\ 0 & 0 & T _ {3 3} \end{array} \right] \begin{array}{c} k - 1 \\ 1 \\ n - k \end{array}
$$

is upper quasi-triangular and that $\lambda \not \in \lambda ( T _ { 1 1 } ) \cup \lambda ( T _ { 3 3 } )$ . It follows that if we solve the linear systems $( T _ { 1 1 } - \lambda I ) w = - u$ and $( T _ { 3 3 } - \lambda I ) ^ { T } z = - v$ then

$$
x = Q {\left[ \begin{array}{l} w \\ 1 \\ 0 \end{array} \right]} \quad {\mathrm{and}} \quad y = Q {\left[ \begin{array}{l} 0 \\ 1 \\ z \end{array} \right]}
$$

are the associated right and left eigenvectors, respectively. Note that the condition of λ is prescribed by

$$
1 / s (\lambda) = \sqrt {(1 + w ^ {T} w) (1 + z ^ {T} z)}.
$$

# 7.6.5 Ascertaining Jordan Block Structures

Suppose that we have computed the real Schur decomposition $A = Q T Q ^ { T }$ , identified clusters of “equal” eigenvalues, and calculated the corresponding block diagonalization $T = Y { \cdot } \mathrm { d i a g } ( T _ { 1 1 } , \dots , T _ { q q } ) Y ^ { - 1 }$ . As we have seen, this can be a formidable task. However, even greater numerical problems confront us if we attempt to ascertain the Jordan block structure of each $T _ { i i }$ . A brief examination of these difficulties will serve to highlight the limitations of the Jordan decomposition.

Assume for clarity that $\lambda ( T _ { i i } )$ is real. The reduction of $T _ { i i }$ to Jordan form begins by replacing it with a matrix of the form $C = \lambda I + N$ , where N is the strictly upper triangular portion of $T _ { i i }$ and where λ, say, is the mean of its eigenvalues.

Recall that the dimension of a Jordan block $J ( \lambda )$ is the smallest nonnegative integer k for which $[ J ( \lambda ) - \lambda I ] ^ { k } = 0$ . Thus, if $p _ { i } = \mathsf { d i m } [ \mathsf { n u l l } ( N ^ { i } ) ]$ , for $i = 0 { : } n$ , then $p _ { i } - p _ { i - 1 }$ equals the number of blocks in $C ^ { \mathrm { { * } } } \mathrm { { s } }$ Jordan form that have dimension i or greater. A concrete example helps to make this assertion clear and to illustrate the role of the SVD in Jordan form computations.

Assume that C is 7-by-7. Suppose we compute the SVD $U _ { 1 } ^ { T } N V _ { 1 } = \Sigma _ { 1 }$ and “discover” that N has rank 3. If we order the singular values from small to large then it follows that the matrix $N _ { 1 } = V _ { 1 } ^ { T } N V _ { 1 }$ has the form

$$
N _ {1} = \left[ \begin{array}{c c} 0 & K \\ 0 & L \end{array} \right] _ {3} ^ {4}.
$$

At this point, we know that the geometric multiplicity of λ is $4 \mathrm { - } \mathrm { i . e } , C \mathrm { s }$ Jordan form has four blocks $( p _ { 1 } - p _ { 0 } = 4 - 0 = 4 )$ .

Now suppose $\tilde { U } _ { 2 } ^ { T } L \tilde { V } _ { 2 } = \Sigma _ { 2 }$ is the SVD of L and that we find that L has unit rank. If we again order the singular values from small to large, then $L _ { 2 } = \tilde { V } _ { 2 } ^ { T } L \tilde { V } _ { 2 }$ clearly has the following structure:

$$
L _ {2} = \left[ \begin{array}{c c c} 0 & 0 & a \\ 0 & 0 & b \\ 0 & 0 & c \end{array} \right].
$$

However, $\lambda ( L _ { 2 } ) = \lambda ( L ) = \{ 0 , 0 , 0 \}$ and so $c = 0$ . Thus, if

$$
V _ {2} = \mathrm{diag} (I _ {4}, \tilde {V} _ {2})
$$

then $N _ { 2 } = V _ { 2 } ^ { T } N _ { 1 } V _ { 2 }$ has the following form:

$$
N _ {2} = \left[ \begin{array}{c c c c c c c} 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & a \\ 0 & 0 & 0 & 0 & 0 & 0 & b \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{array} \right].
$$

Besides allowing us to introduce more zeros into the upper triangle, the SVD of L also enables us to deduce the dimension of the nullspace of $N ^ { 2 }$ . Since

$$
N _ {1} ^ {2} = \left[ \begin{array}{c c} 0 & K L \\ 0 & L ^ {2} \end{array} \right] = \left[ \begin{array}{c c} 0 & K \\ 0 & L \end{array} \right] \left[ \begin{array}{c c} 0 & K \\ 0 & L \end{array} \right]
$$

$\left[ \begin{array} { l } { K } \\ { L } \end{array} \right]$ has full column rank,

$$
p _ {2} = \dim (\operatorname{null} (N ^ {2})) = \dim (\operatorname{null} (N _ {1} ^ {2})) = 4 + \dim (\operatorname{null} (L)) = p _ {1} + 2.
$$

Hence, we can conclude at this stage that the Jordan form of C has at least two blocks of dimension 2 or greater.

Finally, it is easy to see that $N _ { 1 } ^ { 3 } = 0$ , from which we conclude that there is $p _ { 3 } - p _ { 2 }$ $= 7 - 6 = 1$ block of dimension 3 or larger. If we define $V = V _ { 1 } V _ { 2 }$ then it follows that

the decomposition

$$
V ^ {T} C V = \left[ \begin{array}{c c c c c c c} \lambda & 0 & 0 & 0 & \times & \times & \times \\ 0 & \lambda & 0 & 0 & \times & \times & \times \\ 0 & 0 & \lambda & 0 & \times & \times & \times \\ 0 & 0 & 0 & \lambda & \times & \times & \times \\ 0 & 0 & 0 & 0 & \lambda & \times & a \\ 0 & 0 & 0 & 0 & 0 & \lambda & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & \lambda \end{array} \right] \left\{ \begin{array}{l} \text {four blocks of order 1 or larger} \\ \text {two blocks of order 2 or larger} \\ \text {one block of order 3 or larger} \end{array} \right.
$$

displays $C ^ { \mathrm { { * } } } \mathrm { { s } }$ Jordan block structure: two blocks of order 1, one block of order 2, and one block of order 3.

To compute the Jordan decomposition it is necessary to resort to nonorthogonal transformations. We refer the reader to Golub and Wilkinson (1976), K˚agstr¨om and Ruhe (1980a, 1980b), and Demmel (1983) for more details. The above calculations with the SVD amply illustrate that difficult rank decisions must be made at each stage and that the final computed block structure depends critically on those decisions.

# Problems

P7.6.1 Give a complete algorithm for solving a real, n-by-n, upper quasi-triangular system $T x = b$

P7.6.2 Suppose $U ^ { - 1 } A U = \mathrm { d i a g } ( \alpha _ { 1 } , . . . , \alpha _ { m } )$ and $V ^ { - 1 } B V = \mathrm { d i a g } ( \beta _ { 1 } , . ~ . ~ . , \beta _ { n } )$ . Show that if

$$
\phi (X) = A X - X B,
$$

then

$$
\lambda (\phi) = \{\alpha_ {i} - \beta_ {j}: i = 1: m, j = 1: n \}.
$$

What are the corresponding eigenvectors? How can these facts be used to solve $A X - X B = C ?$

P7.6.3 Show that if $Z \in \mathbb { C } ^ { p \times q }$ and

$$
Y = \left[ \begin{array}{c c} I _ {p} & Z \\ 0 & I _ {q} \end{array} \right],
$$

then $\kappa _ { 2 } ( Y ) = [ 2 + \sigma ^ { 2 } + \sqrt { 4 \sigma ^ { 2 } + \sigma ^ { 4 } } ] / 2 \mathrm { ~ w h e r e ~ } \sigma = \| Z \| _ { 2 } .$

P7.6.4 Derive the system $_ { ( 7 . 6 . 5 ) }$ .

P7.6.5 Assume that $T \in \mathbb { R } ^ { n \times n }$ is block upper triangular and partitioned as follows:

$$
T = \left[ \begin{array}{c c c} T _ {1 1} & T _ {1 2} & T _ {1 3} \\ 0 & T _ {2 2} & T _ {2 3} \\ 0 & 0 & T _ {3 3} \end{array} \right], \qquad T \in \mathbb {R} ^ {n \times n}.
$$

Suppose that the diagonal block $T _ { 2 2 }$ is 2-by-2 with complex eigenvalues that are disjoint from $\lambda ( T _ { 1 1 } )$ and $\lambda ( T _ { 3 3 } )$ . Give an algorithm for computing the 2-dimensional real invariant subspace associated with $T _ { 2 2 } \mathrm { { ' s } }$ eigenvalues.

P7.6.6 Suppose $H \in \mathbb { R } ^ { n \times n }$ is upper Hessenberg with a complex eigenvalue $\lambda + i \cdot \mu$ . How could inverse iteration be used to compute x, $\boldsymbol { y } \in \mathbb { R } ^ { n }$ so that $H ( x + i y ) = ( \lambda + i \mu ) ( x + i y ) ?$ Hint: Compare real and imaginary parts in this equation and obtain a 2n-by-2n real system.

# Notes and References for §7.6

Much of the material discussed in this section may be found in the following survey paper:

G.H. Golub and J.H. Wilkinson (1976). “Ill-Conditioned Eigensystems and the Computation of the Jordan Canonical Form,” SIAM Review 18, 578–619.

The problem of ordering the eigenvalues in the real Schur form is the subject of:

A. Ruhe (1970). “An Algorithm for Numerical Determination of the Structure of a General Matrix,” BIT 10, 196–216.   
G.W. Stewart (1976). “Algorithm 406: HQR3 and EXCHNG: Fortran Subroutines for Calculating and Ordering the Eigenvalues of a Real Upper Hessenberg Matrix,” ACM Trans. Math. Softw. 2, 275–280.   
J.J. Dongarra, S. Hammarling, and J.H. Wilkinson (1992). “Numerical Considerations in Computing Invariant Subspaces,” SIAM J. Matrix Anal. Applic. 13, 145–161.   
Z. Bai and J.W. Demmel (1993). “On Swapping Diagonal Blocks in Real Schur Form,” Lin. Alg. Applic. 186, 73–95

Procedures for block diagonalization including the Jordan form are described in:

C. Bavely and G.W. Stewart (1979). “An Algorithm for Computing Reducing Subspaces by Block Diagonalization,” SIAM J. Numer. Anal. 16, 359–367.   
B. K˚agstr¨om and A. Ruhe (1980a). “An Algorithm for Numerical Computation of the Jordan Normal Form of a Complex Matrix,” ACM Trans. Math. Softw. 6, 398–419.   
B. K˚agstr¨om and A. Ruhe (1980b). “Algorithm 560 JNF: An Algorithm for Numerical Computation of the Jordan Normal Form of a Complex Matrix,” ACM Trans. Math. Softw. 6, 437–443.   
J.W. Demmel (1983). “A Numerical Analyst’s Jordan Canonical Form,” PhD Thesis, Berkeley.   
N. Ghosh, W.W. Hager, and P. Sarmah (1997). “The Application of Eigenpair Stability to Block Diagonalization,” SIAM J. Numer. Anal. 34, 1255–1268.   
S. Serra-Capizzano, D. Bertaccini, and G.H. Golub (2005). “How to Deduce a Proper Eigenvalue Cluster from a Proper Singular Value Cluster in the Nonnormal Case,” SIAM J. Matrix Anal. Applic. 27, 82–86.

Before we offer pointers to the literature associated with invariant subspace computation, we remind the reader that in 7.3 we discussed the power method for computing the dominant eigenpair and the method of orthogonal iteration that can be used to compute dominant invariant subspaces. Inverse iteration is a related idea and is the concern of the following papers:

J. Varah (1968). “The Calculation of the Eigenvectors of a General Complex Matrix by Inverse Iteration,” Math. Comput. 22, 785–791.   
J. Varah (1970). “Computing Invariant Subspaces of a General Matrix When the Eigensystem is Poorly Determined,” Math. Comput. 24, 137–149.   
G. Peters and J.H. Wilkinson (1979). “Inverse Iteration, Ill-Conditioned Equations, and Newton’s Method,” SIAM Review 21, 339–360.

I.C.F. Ipsen (1997). “Computing an Eigenvector with Inverse Iteration,” SIAM Review 39, 254–291. In certain applications it is necessary to track an invariant subspace as the matrix changes, see:

L. Dieci and M.J. Friedman (2001). “Continuation of Invariant Subspaces,” Num. Lin. Alg. 8, 317–327.   
D. Bindel, J.W. Demmel, and M. Friedman (2008). “Continuation of Invariant Subsapces in Large Bifurcation Problems,” SIAM J. Sci. Comput. 30, 637–656.

Papers concerned with estimating the error in a computed eigenvalue and/or eigenvector include:

S.P. Chan and B.N. Parlett (1977). “Algorithm 517: A Program for Computing the Condition Numbers of Matrix Eigenvalues Without Computing Eigenvectors,” ACM Trans. Math. Softw. 3, 186–203.   
H.J. Symm and J.H. Wilkinson (1980). “Realistic Error Bounds for a Simple Eigenvalue and Its Associated Eigenvector,” Numer. Math. 35, 113–126.   
C. Van Loan (1987). “On Estimating the Condition of Eigenvalues and Eigenvectors,” Lin. Alg. Applic. 88/89, 715–732.   
Z. Bai, J. Demmel, and A. McKenney (1993). “On Computing Condition Numbers for the Nonsymmetric Eigenproblem,” ACM Trans. Math. Softw. 19, 202–223.

Some ideas about improving computed eigenvalues, eigenvectors, and invariant subspaces may be found in:

J. Varah (1968). “Rigorous Machine Bounds for the Eigensystem of a General Complex Matrix,” Math. Comp. 22, 793–801.   
J.J. Dongarra, C.B. Moler, and J.H. Wilkinson (1983). “Improving the Accuracy of Computed Eigenvalues and Eigenvectors,” SIAM J. Numer. Anal. 20, 23–46.

J.W. Demmel (1987). “Three Methods for Refining Estimates of Invariant Subspaces,” Comput. 38, 43–57.   
As we have seen, the sep(.,.) function is of great importance in the assessment of a computed invariant subspace. Aspects of this quantity and the associated Sylvester equation are discussed in:   
J. Varah (1979). “On the Separation of Two Matrices,” SIAM J. Numer. Anal. 16, 212–222.   
R. Byers (1984). “A Linpack-Style Condition Estimator for the Equation $A X - X B ^ { T } = C , ^ { , }$ IEEE Trans. Autom. Contr. AC-29, 926–928.   
M. Gu and M.L. Overton (2006). “An Algorithm to Compute Sepλ,” SIAM J. Matrix Anal. Applic. 28, 348–359.   
N.J. Higham (1993). “Perturbation Theory and Backward Error for $\mathrm { A X - X B } = \mathrm { C } , \mathrm { " } \ B I T \ 3 3 , 1 2 4 \mathrm { - } 1 3 6 .$ .   
Sylvester equations arise in many settings, and there are many solution frameworks, see:   
R.H. Bartels and G.W. Stewart (1972). “Solution of the Equation $A X + X B = C , "$ Commun. ACM 15, 820–826.   
G.H. Golub, S. Nash, and C. Van Loan (1979). “A Hessenberg-Schur Method for the Matrix Problem AX + XB = C,” IEEE Trans. Autom. Contr. AC-24, 909–913.   
K. Datta (1988). “The Matrix Equation XA − BX = R and Its Applications,” Lin. Alg. Applic. 109, 91–105.   
B. K˚agstr¨om and P. Poromaa (1992). “Distributed and Shared Memory Block Algorithms for the Triangular Sylvester Equation with $\mathrm { s e p } ^ { - 1 }$ Estimators,” SIAM J. Matrix Anal. Applic. 13, 90– 101.   
J. Gardiner, M.R. Wette, A.J. Laub, J.J. Amato, and C.B. Moler (1992). “Algorithm 705: A FORTRAN-77 Software Package for Solving the Sylvester Matrix Equation $A X B ^ { T } + C X D ^ { T } = E ,$ ACM Trans. Math. Softw. 18, 232–238.   
V. Simoncini (1996). “On the Numerical Solution of AX -XB =C,” BIT 36, 814–830.   
C.H. Bischof, B.N Datta, and A. Purkayastha (1996). “A Parallel Algorithm for the Sylvester Observer Equation,” SIAM J. Sci. Comput. 17, 686–698.   
D. Calvetti, B. Lewis, L. Reichel (2001). “On the Solution of Large Sylvester-Observer Equations,” Num. Lin. Alg. 8, 435–451.

The constrained Sylvester equation problem is considered in:

J.B. Barlow, M.M. Monahemi, and D.P. O’Leary (1992). “Constrained Matrix Sylvester Equations,” SIAM J. Matrix Anal. Applic. 13, 1–9.

A.R. Ghavimi and A.J. Laub (1996). “Numerical Methods for Nearly Singular Constrained Matrix Sylvester Equations.” SIAM J. Matrix Anal. Applic. 17, 212–221.

The Lyapunov problem F X $+ X F ^ { T } = - C$ where C is non-negative definite has a very important role to play in control theory, see:

G. Hewer and C. Kenney (1988). “The Sensitivity of the Stable Lyapunov Equation,” SIAM J. Control Optim 26, 321–344.

A.R. Ghavimi and A.J. Laub (1995). “Residual Bounds for Discrete-Time Lyapunov Equations,” IEEE Trans. Autom. Contr. 40, 1244–1249.

J.-R. Li and J. White (2004). “Low-Rank Solution of Lyapunov Equations,” SIAM Review 46, 693– 713.

Several authors have considered generalizations of the Sylvester equation, i.e., $\Sigma F _ { i } X G _ { i } = C .$ These include:

P. Lancaster (1970). “Explicit Solution of Linear Matrix Equations,” SIAM Review 12, 544–566.

H. Wimmer and A.D. Ziebur (1972). “Solving the Matrix Equations $\Sigma f _ { p } ( A ) g _ { p } ( A ) = C , "$ SIAM Review 14, 318–323.

W.J. Vetter (1975). “Vector Structures and Solutions of Linear Matrix Equations,” Lin. Alg. Applic. 10, 181–188.
