# 12.1 Linear Systems with Displacement Structure

If $A \in \mathbb { R } ^ { n \times n }$ has rank r, then it has a (non-unique) product representation of the form $U V ^ { T }$ where $U , V \in \mathbb { R } ^ { n \times r }$ . Note that if $r \ll n$ , then the product representation is much more compact than the explicit representation that encodes each $a _ { i j }$ . In addition to the obvious storage economies, the product representation supports fast computation. If the product representation is fully utilized, then the n-by-n matrix-matrix product $A B \ = \ U ( V ^ { T } B )$ is $O ( n ^ { 2 } r )$ instead of $O ( n ^ { 3 } )$ . Likewise, by applying the Sherman-Morrison-Woodbury formula, the solution to a linear system of the form $( I + U V ^ { T } ) x = b$ is $O ( n r + r ^ { 3 } )$ instead of $O ( n ^ { 3 } )$ . The message is simple in both cases: work with U and V and not their explicit product $U V ^ { T }$ .

In this section we continue in this direction by discussing “low-rank” way to represent Cauchy, Toeplitz, and Hankel matrices together with some of their generalizations. The data-sparse representation supports fast stable linear equation solving. The key idea is to turn explicit rank-1 updates that are at the heart of Gaussian elimination into equivalent, inexpensive updates of their representation. Our presentation is based on Gohberg, Kailath, and Olshevsky (1995) and Gu (1998).

# 12.1.1 Displacement Rank

If $F , G \in \mathbb { R } ^ { n \times n }$ and the Sylvester map

$$
X \rightarrow F X - X G \tag {12.1.1}
$$

is nonsingular, then the {F, G}-displacement rank of $A \in \mathbb { R } ^ { n \times n }$ is defined by

$$
\operatorname{rank} _ {\{F, G \}} (A) = \operatorname{rank} (F A - A G). \tag {12.1.2}
$$

Recall from §7.6.3 that the Sylvester map is nonsingular provided $\lambda ( F ) \cup \lambda ( G ) = \varnothing$ Note that if rank $\{ F , G \}  ( A ) = r$ , then we can write

$$
F A - A G = R S ^ {T}, \quad R, S \in \mathbb {R} ^ {n \times r}. \tag {12.1.3}
$$

The matrices R and S are generators for A with respect to F and G, a term that makes sense since we can generate A (or part of A) by working with this equation. If $r \ll n$ , then R and S define a data-sparse representation for A. Of course, for this representation to be of interest F and G must be sufficiently simple so that the reconstruction of A via (12.1.3) is cheap.

# 12.1.2 Cauchy-Like Matrices

If $\boldsymbol \omega \in \mathbb { R } ^ { n }$ and $\lambda \in \mathbb { R } ^ { n }$ and $\omega _ { k } \neq \lambda _ { j }$ for all k and j, then the n-by-n matrix $A = \left( a _ { k j } \right)$ defined by

$$
a _ {k j} = \frac {1}{\omega_ {k} - \lambda_ {j}}
$$

is a Cauchy matrix. Note that if

$$
\Omega = \operatorname{diag} \left(\omega_ {1}, \dots , \omega_ {n}\right), \quad \Lambda = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right),
$$

then

$$
[ \Omega A - A \Lambda ] _ {k j} = \frac {\omega_ {k}}{\omega_ {k} - \lambda_ {j}} - \frac {\lambda_ {j}}{\omega_ {k} - \lambda_ {j}} = 1.
$$

If $e \in \mathbb { R } ^ { n }$ is the vector of all $\mathrm { 1 ^ { \circ } s } .$ then

$$
\Omega A - A \Lambda = e e ^ {T}
$$

and thus ran $\mathfrak { c } _ { \{ \Omega , \Lambda \} } ( A ) = 1$ .

More generally, if $R \in \mathbb { R } ^ { n \times r }$ and $S \in \mathbb { R } ^ { n \times r }$ have rank r, then any matrix A that satisfies

$$
\Omega A - A \Lambda = R S ^ {T} \tag {12.1.4}
$$

is a Cauchy-like matrix. This just means that

$$
a _ {k j} = \frac {r _ {k} ^ {T} s _ {j}}{\omega_ {k} - \lambda_ {j}}
$$

where

$$
R ^ {T} = \left[ \begin{array}{c c c} r _ {1} & \dots & r _ {n} \end{array} \right], \qquad S ^ {T} = \left[ \begin{array}{c c c} s _ {1} & \dots & s _ {n} \end{array} \right]
$$

are column partitionings. Note that R and S are generators with respect to Ω and Λ and that $O ( r )$ flops are required to reconstruct a matrix entry $a _ { k j }$ from (12.1.4).

# 12.1.3 The Apparent Loss of Structure

Suppose

$$
A = \left[ \begin{array}{l l} \alpha & g ^ {T} \\ f & B \end{array} \right], \qquad \alpha \in \mathbb {R}, f, g \in \mathbb {R} ^ {n - 1}, B \in \mathbb {R} ^ {(n - 1) \times (n - 1)},
$$

and assume $\alpha \neq 0$ . The first step in Gaussian elimination produces

$$
A _ {1} = B - \frac {1}{\alpha} f g ^ {T}
$$

and the factorization

$$
A = \left[ \begin{array}{c c} 1 & 0 \\ f / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{c c} \alpha & g ^ {T} \\ 0 & A _ {1} \end{array} \right].
$$

Let us examine the structure of $A _ { 1 }$ given that A is a Cauchy matrix. If $n = 4$ and $a _ { k j } = 1 / ( \omega _ { k } - \lambda _ { j } )$ , then

$$
A _ {1} = \left[ \begin{array}{c c c} \frac {1}{\omega_ {2} - \lambda_ {2}} & \frac {1}{\omega_ {2} - \lambda_ {3}} & \frac {1}{\omega_ {2} - \lambda_ {4}} \\ \frac {1}{\omega_ {3} - \lambda_ {2}} & \frac {1}{\omega_ {3} - \lambda_ {3}} & \frac {1}{\omega_ {3} - \lambda_ {4}} \\ \frac {1}{\omega_ {4} - \lambda_ {2}} & \frac {1}{\omega_ {4} - \lambda_ {3}} & \frac {1}{\omega_ {4} - \lambda_ {4}} \end{array} \right] - \left[ \begin{array}{c} \frac {\omega_ {1} - \lambda_ {1}}{\omega_ {2} - \lambda_ {1}} \\ \frac {\omega_ {1} - \lambda_ {1}}{\omega_ {3} - \lambda_ {1}} \\ \frac {\omega_ {1} - \lambda_ {1}}{\omega_ {4} - \lambda_ {1}} \end{array} \right] \left[ \begin{array}{c} \frac {1}{\omega_ {1} - \lambda_ {2}} \\ \frac {1}{\omega_ {1} - \lambda_ {3}} \\ \frac {1}{\omega_ {1} - \lambda_ {4}} \end{array} \right] ^ {T}.
$$

If we choose to work with the explicit representation of A, then for general n this update requires $O ( n ^ { 2 } )$ work even though it is highly structured and involves $O ( n )$ data. And worse, all subsequent steps in the factorization process essentially deal with general matrices rendering an LU computation that is $O ( n ^ { 3 } )$ .

# 12.1.4 Displacement Rank and Rank-1 Updates

The situation is much happier if we replace the explicit transition from A to $A _ { 1 }$ with a transition that involves updating data sparse representations. The key to developing a fast LU factorization for a Cauchy-like matrix is to recognize that rank-1 updates preserve displacement rank. Here is the result that makes it all possible.

Theorem 12.1.1. Suppose $A \in \mathbb { R } ^ { n \times n }$ satisfies

$$
\Omega A - A \Lambda = R S ^ {T} \tag {12.1.5}
$$

where R, $S \in \mathbb { R } ^ { n \times r }$ and

$$
\Omega = \operatorname{diag} \left(\omega_ {1}, \dots , \omega_ {n}\right), \quad \Lambda = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right)
$$

have no common diagonal entries. If

$$
A = \left[ \begin{array}{c c} \alpha & g ^ {T} \\ f & B \end{array} \right], \qquad R = \left[ \begin{array}{c} r _ {1} ^ {T} \\ R _ {1} \end{array} \right], \qquad S = \left[ \begin{array}{c} s _ {1} ^ {T} \\ S _ {1} \end{array} \right]
$$

are conformably partitioned, $\alpha \neq 0$ , and

$$
\Omega_ {1} = \operatorname{diag} \left(\omega_ {2}, \dots , \omega_ {n}\right), \quad \Lambda_ {1} = \operatorname{diag} \left(\lambda_ {2}, \dots , \lambda_ {n}\right),
$$

then

$$
\Omega_ {1} A _ {1} - A _ {1} \Lambda_ {1} = \tilde {R} _ {1} \tilde {S} _ {1} ^ {T} \tag {12.1.6}
$$

where

$$
A _ {1} = B - \frac {f g ^ {T}}{\alpha}, \qquad \tilde {R} _ {1} = R _ {1} - \frac {1}{\alpha} f r _ {1} ^ {T}, \qquad \tilde {S} _ {1} = S _ {1} - \frac {1}{\alpha} g s _ {1} ^ {T}.
$$

Proof. By comparing blocks in (12.1.5) we see that

$$
\begin{array}{l} (1, 1): (\omega_ {1} - \lambda_ {1}) \alpha = r _ {1} ^ {T} s _ {1}, \\ (1, 2): g ^ {T} \Lambda_ {1} = \omega_ {1} g ^ {T} - r _ {1} ^ {T} S _ {1} ^ {T}, \\ (2, 1): \Omega_ {1} f = R _ {1} s _ {1} + \lambda_ {1} f, \\ (2, 2): \Omega_ {1} B - B \Lambda_ {1} = R _ {1} S _ {1} ^ {T}, \\ \end{array}
$$

and so

$$
\begin{array}{l} \Omega_ {1} A _ {1} - A _ {1} \Lambda_ {1} = \Omega_ {1} \left(B - \frac {1}{\alpha} f g ^ {T}\right) - \left(B - \frac {1}{\alpha} f g ^ {T}\right) \Lambda_ {1} \\ = \left(\Omega_ {1} B - B \Lambda_ {1}\right) - \frac {1}{\alpha} \left(\left(\Omega_ {1} f\right) g ^ {T} - f \left(g ^ {T} \Lambda_ {1}\right)\right) \\ = R _ {1} S _ {1} ^ {T} - \frac {1}{\alpha} \left((R _ {1} s _ {1} + \lambda_ {1} f) g ^ {T} - f (\omega_ {1} g ^ {T} - r _ {1} ^ {T} S _ {1} ^ {T})\right) \\ = R _ {1} S _ {1} ^ {T} - \frac {1}{\alpha} \left((R _ {1} s _ {1}) g ^ {T} + f (r _ {1} ^ {T} S _ {1} ^ {T}) - \frac {r _ {1} ^ {T} s _ {1}}{\alpha} f g ^ {T}\right) \\ = \left(R _ {1} - \frac {1}{\alpha} f r _ {1} ^ {T}\right) \left(S _ {1} - \frac {1}{\alpha} g s _ {1} ^ {T}\right) ^ {T} = \tilde {R} _ {1} \tilde {S} _ {1} ^ {T}. \\ \end{array}
$$

This confirms (12.1.6) and completes the proof of the theorem.

The theorem says that

$$
\operatorname{rank} _ {\{\Omega , \Lambda \}} (A) \leq r \quad \Rightarrow \quad \operatorname{rank} _ {\{\Omega_ {1}, \Lambda_ {1} \}} (A _ {1}) \leq r.
$$

This suggests that instead of updating A explicitly to get $A _ { 1 }$ at a cost of $O ( n ^ { 2 } )$ flops, we should update A’s representation {Ω, Λ, R, S} at a cost of $O ( n r )$ flops to get $A _ { 1 } \mathrm { { ' } s }$ representation $\{ \Omega _ { 1 } , \Lambda _ { 1 } , \tilde { R } _ { 1 } , \tilde { S } _ { 1 } \}$ .

# 12.1.5 Fast LU for Cauchy-Like Matrices

Based on Theorem 12.1.1 we can specify a fast LU procedure for Cauchy-like matrices. If A satisfies (12.1.5) and has an LU factorization, then it can be computed using the function LUdisp defined as follows:

Algorithm 12.1.1 If $\omega \in \mathbb { R } ^ { n }$ and $\lambda \in \mathbb { R } ^ { n }$ have no common components, R, $S \in \mathbb { R } ^ { n \times r }$ , and $\Omega A - A \Lambda = R S ^ { T }$ where $\boldsymbol { \Omega } = \mathrm { d i a g } ( \omega _ { 1 } , \dots , \omega _ { n } )$ and $\boldsymbol { \Lambda } = \operatorname { d i a g } ( \lambda _ { 1 } , \ldots , \lambda _ { n } )$ , then the following function computes the LU factorization $A = L U$ .

function $[ L , U ] = { \mathsf { L U d i s p } } ( \omega , \lambda , R , S , n )$

$$
r _ {1} ^ {T} = R (1,:), R _ {1} = R (2: n,:) \tag {1.1}
$$

$$
s _ {1} ^ {T} = S (1,:), S _ {1} = S (2: n,:) \tag {1.1}
$$

$\mathbf { i f } \ n = 1$

$$
L = 1
$$

$$
U = r _ {1} ^ {T} s _ {1} / (\omega_ {1} - \lambda_ {1})
$$

else

$$
a = \left(R s _ {1}\right). / \left(\omega - \lambda_ {1}\right)
$$

$$
\alpha = a _ {1 1}
$$

$$
f = a (2: n)
$$

$$
g = \left(S _ {1} r _ {1}\right). / \left(\omega_ {1} - \lambda (2: n)\right)
$$

$$
\tilde {R} _ {1} = R _ {1} - f r _ {1} ^ {T} / \alpha
$$

$$
\tilde {S} _ {1} = S _ {1} - g s _ {1} ^ {T} / \alpha
$$

$$
\left[ L _ {1}, U _ {1} \right] = \mathsf {L U d i s p} (\omega (2: n), \lambda (2: n), \tilde {R} _ {1}, \tilde {S} _ {1}, n - 1)
$$

$$
L = \left[ \begin{array}{c c} 1 & 0 \\ f / \alpha & L _ {1} \end{array} \right]
$$

$$
U = \left[ \begin{array}{l l} \alpha & g ^ {T} \\ 0 & U _ {1} \end{array} \right]
$$

end

The nonrecursive version would have the following structure:

Let $R ^ { ( 1 ) }$ and $S ^ { ( 1 ) }$ be the generators of $A = A ^ { ( 1 ) }$ with respect to diag(ω) and diag(λ).

for k = 1:n − 1

Use ω(k:n), λ(k:n), R(k) and S(k) to compute the first row and column of

$$
A ^ {(k)} = \left[ \begin{array}{c c} \alpha & g ^ {T} \\ f & B \end{array} \right].
$$

$$
L (k + 1: n, k) = f / \alpha , U (k, k) = \alpha , U (k, k + 1: n) = g ^ {T}
$$

$\begin{array} { r l } & { \mathrm { D e t e r m i n e ~ t h e ~ g e n e r a t o r s ~ } R ^ { ( k + 1 ) } \mathrm { ~ a n d ~ } S ^ { ( k + 1 ) } \mathrm { ~ o f ~ } A ^ { ( k + 1 ) } = B - f g ^ { T } / \alpha } \\ & { \qquad \mathrm { w i t h ~ r e s p e c t ~ t o ~ d i a g } ( \omega ( k \cdot n ) ) \mathrm { ~ a n d ~ d i a g } ( \lambda ( k \cdot n ) ) . } \\ & { \mathrm { e n d ~ } } \end{array}$

$$
U (n, n) = R ^ {(n)} \cdot S ^ {(n)} / (\omega_ {n} - \lambda_ {n})
$$

A careful accounting reveals that $2 n ^ { 2 } r$ flops are required.

# 12.1.6 Pivoting

The procedure just developed has numerical difficulties if a small α shows up during the recursion. To guard against this we show how to incorporate a pivoting strategy. Suppose $A \in \mathbb { R } ^ { n \times n }$ is a Cauchy-like matrix that satisfies the displacement equation

$$
\Omega A - A \Lambda = R S ^ {T}
$$

for diagonal matrices Ω and Λ and n-by-r matrices R and S. If P and Q are n-by-n permutations, then

$$
(P \Omega P ^ {T} (P A Q ^ {T}) - (P A Q ^ {T}) (Q \Lambda Q ^ {T}) = (P R) (Q S) ^ {T}.
$$

This shows that

$$
\tilde {A} = P A Q ^ {T}
$$

is a Cauchy-like matrix having generators

$$
\tilde {R} = P R, \quad \tilde {S} = Q S
$$

with respect to the diagonal matrices

$$
\tilde {\Omega} = P \Omega P ^ {T}, \qquad \tilde {\Lambda} = Q \Lambda Q ^ {T}.
$$

Thus, it is easy to track row and column permutations in the the displacement representation:

$$
A \rightarrow P A Q ^ {T}, \equiv \{\Omega , \Lambda , R, S \} \rightarrow \{P \Omega P ^ {T}, Q \Lambda Q ^ {T}, P R, Q S \}.
$$

By taking advantage of this, it is a simple matter to incorporate partial pivoting in LUdisp and to emerge with the factorization $P A = L U$ :

Algorithm 12.1.2 If $\omega \in \mathbb { R } ^ { n }$ and $\lambda \in \mathbb { R } ^ { n }$ have no common components, R, $S \in \mathbb { R } ^ { n \times r }$ , and $\Omega A - A \Lambda = R S ^ { T }$ , then the following function computes the LU-with-pivoting factorization $P A = L U$ , where $\boldsymbol \Omega = \mathrm { d i a g } ( \omega _ { 1 } , \dots , \omega _ { n } )$ and $\boldsymbol { \Lambda } = \operatorname { d i a g } ( \lambda _ { 1 } , \ldots , \lambda _ { n } )$ .

function $[ L , U , P ] = \mathsf { L U d i s p P i v } ( \omega , \lambda , R , S , n )$

Define $r _ { 1 } , R _ { 1 } , s _ { 1 }$ and $S _ { 1 }$ by $R = { \left[ \begin{array} { l } { r _ { 1 } ^ { T } } \\ { R _ { 1 } } \end{array} \right] } { \mathrm { ~ a n d ~ } } S = { \left[ \begin{array} { l } { s _ { 1 } ^ { T } } \\ { S _ { 1 } } \end{array} \right] }$

if $n = 1$

$$
L = 1
$$

$$
U = r _ {1} ^ {T} s _ {1} / (\omega_ {1} - \lambda_ {1})
$$

else

$$
a = \left(R s _ {1}\right). / \left(\omega - \lambda_ {1}\right)
$$

Determine permutation $P \in \mathbb { R } ^ { n \times n }$ so that $[ P a ] _ { 1 }$ is maximal and

update: $a = P a , R = P R , \omega = P \omega .$ .

$$
\alpha = a _ {1}
$$

$$
f = a (2: n)
$$

$$
g = \left(S _ {1} r _ {1}\right). / \left(\omega_ {1} - \lambda (2: n)\right)
$$

$$
\tilde {R} _ {1} = R _ {1} - f r _ {1} ^ {T} / \alpha
$$

$$
\tilde {S} _ {1} = S _ {1} - g s _ {1} ^ {T} / \alpha
$$

$$
\left[ L _ {1}, U _ {1}, P _ {1} \right] = \text { LUdispPiv } (\omega (2: n), \lambda (2: n), \tilde {R} _ {1}, \tilde {S} _ {1}, n - 1)
$$

$$
L = \left[ \begin{array}{c c} 1 & 0 \\ P _ {1} f / \alpha & L _ {1} \end{array} \right]
$$

$$
U = \left[ \begin{array}{l l} \alpha & g ^ {T} \\ 0 & U _ {1} \end{array} \right]
$$

$$
P = \left[ \begin{array}{c c} 1 & 0 \\ 0 & P _ {1} \end{array} \right] P
$$

end

The processing of the recursive call is based on the fact that if

$$
P A = \left[ \begin{array}{c c} \alpha & g ^ {T} \\ f & B \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ f / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{c c} \alpha & g ^ {T} \\ 0 & A _ {1} \end{array} \right], \qquad A _ {1} = B - \frac {1}{\alpha} f g ^ {T},
$$

and $P _ { 1 } A _ { 1 } = L _ { 1 } U _ { 1 }$ , then

$$
\left[ \begin{array}{c c} 1 & 0 \\ 0 & P _ {1} \end{array} \right] P A = \left[ \begin{array}{c c} 1 & 0 \\ P _ {1} f / \alpha & L _ {1} \end{array} \right] \left[ \begin{array}{c c} \alpha & g ^ {T} \\ 0 & U _ {1} \end{array} \right].
$$

For LUdispPiv implementation details and a proof of its stability, see Gu (1998).

# 12.1.7 Toeplitz-Like Matrices and Hankel-Like Matrices

Recall from §4.7 that a Toeplitz matrix is constant along each of its diagonals. For example, if $c \in \mathbb { R } ^ { n - 1 } , \ \tau \in \mathbb { R }$ , and $r \in \mathbb { R } ^ { n - 1 }$ are given, then the matrix $T \in \mathbb { R } ^ { n \times n }$ defined by

$$
t _ {i j} = \left\{ \begin{array}{l l} c _ {i - j} & \text {if} i > j, \\ \tau & \text {if} i = j, \\ r _ {j - i} & \text {if} j > i, \end{array} \right.
$$

is Toeplitz, e.g.,

$$
T = \left[ \begin{array}{c c c c c} \tau & r _ {1} & r _ {2} & r _ {3} & r _ {4} \\ c _ {1} & \tau & r _ {1} & r _ {2} & r _ {3} \\ c _ {2} & c _ {1} & \tau & r _ {1} & r _ {2} \\ c _ {3} & c _ {2} & c _ {1} & \tau & r _ {1} \\ c _ {4} & c _ {3} & c _ {2} & c _ {1} & \tau \end{array} \right].
$$

To expose the low-displacement-rank structure of a Toeplitz matrix, we define matrices $Z _ { \phi }$ and $Y _ { \gamma , \delta }$ analogously to their $n = 5$ instances:

$$
Z _ {\phi} = \left[ \begin{array}{c c c c c} 0 & 0 & 0 & 0 & \phi \\ 1 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \end{array} \right], \quad Y _ {\gamma , \delta} = \left[ \begin{array}{c c c c c} \gamma & 1 & 0 & 0 & 0 \\ 1 & 0 & 1 & 0 & 0 \\ 0 & 1 & 0 & 1 & 0 \\ 0 & 0 & 1 & 0 & 1 \\ 0 & 0 & 0 & 1 & \delta \end{array} \right]. \tag {12.1.7}
$$

It can be shown that

$$
Z _ {1} T - T Z _ {- 1} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right], \quad \operatorname{rank} _ {\{Z _ {1}, Z _ {- 1} \}} (T) \leq 2, \tag {12.1.8}
$$

$$
Y _ {0 0} T - T Y _ {1 1} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & 0 & 0 & 0 & \times \\ \times & 0 & 0 & 0 & \times \\ \times & 0 & 0 & 0 & \times \\ \times & \times & \times & \times & \times \end{array} \right], \quad \operatorname{rank} _ {\{Y _ {0 0}, Y _ {1 1} \}} (T) \leq 4. \tag {12.1.9}
$$

Furthermore, $\lambda ( Z _ { - 1 } ) \cup \lambda ( Z _ { 1 } ) = \emptyset \ \mathrm { a n d } \ \lambda ( Y _ { 0 0 } ) \cup \lambda ( Y _ { 1 1 } ) = \emptyset .$

A Hankel matrix is constant along its antidiagonals, e.g.,

$$
H = \left[ \begin{array}{c c c c c} c _ {4} & c _ {3} & c _ {2} & c _ {1} & \tau \\ c _ {3} & c _ {2} & c _ {1} & \tau & r _ {1} \\ c _ {2} & c _ {1} & \tau & r _ {1} & r _ {2} \\ c _ {1} & \tau & r _ {1} & r _ {2} & r _ {3} \\ \tau & r _ {1} & r _ {2} & r _ {3} & r _ {4} \end{array} \right].
$$

Note that if $H \in \mathbb { R } ^ { n \times n }$ is Hankel, then $\mathcal { E } _ { n } H$ is Toeplitz, and so it is not surprising that

Hankel and Toeplitz matrices have similar displacement rank properties:

$$
Z _ {1} ^ {T} H - H Z _ {- 1} = \left[ \begin{array}{c c c c c} 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & \times \\ \times & \times & \times & \times & \times \end{array} \right], \quad \operatorname{rank} _ {\{Z _ {1} ^ {T}, Z _ {- 1} \}} (H) \leq 2, \tag {12.1.10}
$$

$$
Y _ {0 0} H - H Y _ {1 1} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & 0 & 0 & 0 & \times \\ \times & 0 & 0 & 0 & \times \\ \times & 0 & 0 & 0 & \times \\ \times & \times & \times & \times & \times \end{array} \right], \quad \operatorname{rank} _ {\left\{Y _ {0 0}, Y _ {1 1} \right\}} (H) \leq 4. \tag {12.1.11}
$$

It follows from (12.1.9) and (12.1.11) that if $A = T + H$ is the sum of a Toeplitz matrix and a Hankel matrix, then rank $\{ Y _ { 0 0 } , Y _ { 1 1 } \} { } ^ { ( A ) } \leq 4$ .

The classes of Toeplitz, Hankel, and Toeplitz-plus-Hankel matrices can be expanded through the notion of low displacement rank. Analogous to how we defined Cauchy-like matrices in (12.1.4) we have the following, assuming that $R \in \mathbb { R } ^ { n \times r }$ , $S \in \mathbb { R } ^ { n \times r }$ , and $r \ll n \colon$

$$
\left\{ \begin{array}{l} Z _ {1} A - A Z _ {- 1} = R S ^ {T} \\ Z _ {1} ^ {T} A - A Z _ {- 1} = R S ^ {T} \\ Y _ {0 0} A - A Y _ {1 1} = R S ^ {T} \end{array} \right\} \text {means that} A \text {is} \left\{ \begin{array}{l} \text {Toeplitz - like} \\ \text {Hankel - like} \\ \text {Toeplitz - plus - Hankel - like} \end{array} \right\}.
$$

Our next task is to show that a linear system with any of these properties can be efficiently converted to a Cauchy-like system and solved with $O ( n ^ { 2 } r )$ work.

# 12.1.8 Fast Solvers via Conversion to Cauchy-Like Form

Suppose

$$
F A - A G = R S ^ {T}, \quad A, F, G \in \mathbb {R} ^ {n \times n}, R, S \in \mathbb {R} ^ {n \times r}, r \ll n,
$$

and that F and G are diagonalizable:

$$
X _ {F} ^ {- 1} F X _ {F} = \operatorname{diag} \left(\omega_ {1}, \dots , \omega_ {n}\right) = \Omega ,
$$

$$
X _ {G} ^ {- 1} G X _ {G} = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right) = \Lambda .
$$

For clarity we assume that F and G have real eigenvalues. It follows from

$$
(X _ {F} ^ {- 1} F X _ {F}) (X _ {F} ^ {- 1} A X _ {G}) - (X _ {F} ^ {- 1} A X _ {G}) (X _ {G} ^ {- 1} G X _ {F}) = (X _ {F} ^ {- 1} R) (X _ {G} ^ {T} S) ^ {T}
$$

that

$$
\Omega \tilde {A} - \tilde {A} \Lambda = \tilde {R} \tilde {S} ^ {T}
$$

where $\tilde { A } = X _ { r } ^ { - 1 } A X _ { G } , \tilde { R } = X _ { r } ^ { - 1 } R$ , and $\tilde { S } = X _ { G } ^ { T } S$ Thus, $\tilde { A }$ is Cauchy-like and we can go about solving the given linear system $A x = b$ as follows:

Step 1. Compute $\tilde { R } = X _ { \scriptscriptstyle F } ^ { - 1 } R , \tilde { S } = X _ { \scriptscriptstyle G } ^ { T } S , \tilde { b } = X _ { \scriptscriptstyle F } ^ { - 1 } b ,$ and $\tilde { A } = X _ { \scriptscriptstyle F } ^ { - 1 } A X _ { \scriptscriptstyle G }$

Step 2. Use Algorithm 12.1.2 to compute $P \tilde { A } = L U$

Step 3. Use $P \tilde { A } = L U$ to solve ${ \tilde { A } } { \tilde { x } } = { \tilde { b } } .$

Step 4. Compute $x = X _ { G } \tilde { x }$ .

This will not be an attractive framework unless the matrices F and G have fast eigensystems, a concept introduced in §4.8. Fortunately, this is the case for the matrices $Z _ { 1 }$ , $Z _ { - 1 } , Y _ { 0 0 }$ and $Y _ { 1 1 }$ . For example,

$$
\mathcal {S} _ {n} ^ {T} Y _ {0 0} \mathcal {S} _ {n} = 2 \cdot \mathrm{diag} \left(\cos \left(\frac {\pi}{n + 1}\right), \ldots , \cos \left(\frac {n \pi}{n + 1}\right)\right), \tag {12.1.12}
$$

$$
\mathcal {C} _ {n} ^ {T} Y _ {1 1} \mathcal {C} _ {n} = 2 \cdot \mathrm{diag} \left(1, \cos \left(\frac {\pi}{n}\right), \ldots , \cos \left(\frac {(n - 1) \pi}{n}\right)\right), \tag {12.1.13}
$$

where $S _ { n }$ is the sine transform (DST-I) matrix

$$
[ \mathcal {S} _ {n} ] _ {k j} = \sqrt {\frac {2}{n + 1}} \cdot \sin \left(\frac {k j \pi}{n + 1}\right),
$$

and $\mathcal { C } _ { n }$ is the cosine transform (DCT-II) matrix

$$
\left[ \mathcal {C} _ {n} \right] _ {k j} = \sqrt {\frac {2}{n}} \cdot \cos \left(\frac {(2 k - 1) (j - 1) \pi}{2 n}\right) q _ {j}, \qquad q _ {j} = \left\{ \begin{array}{l l} 1 / \sqrt {2} & \text {if} j = 1, \\ 1 & \text {if} j > 1. \end{array} \right.
$$

This allows products like $\scriptstyle { S _ { n } R }$ and $\mathcal { C } _ { n } ^ { T } S$ to be computed with $O ( r n \log n )$ flops. In short, Step 3 in the above framework is the most expensive step in the process and it involves $O ( n ^ { 2 } r )$ work. See Gohberg, Kailath, and Olshevsky (1995) and Gu (1998) for details and related references.

# Problems

P12.1.1 Refer to (12.1.8) and (12.1.9). (a) Show that if $Z _ { 1 } X - X Z _ { - 1 } = 0$ , then X = 0. (b) Show that if $Y _ { 0 0 } X - X Y _ { 1 1 } = 0$ , then $\dot { X } = 0 .$

P12.1.2 Develop a nonrecursive version of Algorithm 12.1.2.

P12.1.3 (a) If $T \in \mathbb { R } ^ { n \times n }$ is Toeplitz, show how to compute $R , S \in \mathbb { R } ^ { n \times 2 }$ so that $Z _ { 1 } T - T Z _ { - 1 } = R S ^ { T }$ . (b) Suppose $R , S \in \mathbb { R } ^ { n \times r }$ and $\bar { T } \in \mathbb { R } ^ { n \times n }$ satisfy $Z _ { 1 } T - \bar { T } Z _ { - 1 } = R S ^ { T }$ − −  . Give an algorithm that computes $u = T ( : , 1 )$ and $v = T ( 1 , : ) ^ { T }$ .

P12.1.4 (a) If $T \in \mathbb { R } ^ { n \times n }$ is Toeplitz, show how to compute R, $S \in \mathbb { R } ^ { n \times 4 }$ so that $Y _ { 0 0 } T - T Y _ { 1 1 } = R S ^ { T }$ . (b) Suppose $R , S \in \mathbb { R } ^ { n \times r }$ and $\bar { T } \in \mathbb { R } ^ { n \times n }$ satisfy $Y _ { 0 0 } T _ { - } T Y _ { 1 1 } = R S ^ { T }$ − . Give an algorithm that computes $u = T ( : , 1 )$ and $v = T ( 1 , : ) ^ { T }$ .

P12.1.5 Verify(12.1.13).

P12.1.6 Show that if $A \in \mathbb { R } ^ { n \times n }$ is defined by

$$
a _ {i j} = \int_ {a} ^ {b} \cos (k \theta) \cos (j \theta) d \theta
$$

then A is the sum of a Hankel matrix and Toeplitz matrix. Hint: Make use of the identity $\cos ( u + v ) =$ $\cos ( u ) \cos ( v ) - \sin ( u ) \sin ( v )$ .

# Notes and References for §12.1

For a general introduction to the area of fast algorithms for structured matrices we recommend:

T. Kailath and A. H. Sayed (eds) (1999). Fast Reliable Algorithms for Matrices with Structure, SIAM Publications, Philadelphia, PA.   
V. Olshevsky (ed.) (2000). Structured Matrices in Mathematics, Computer Science, and Engineering I and II, AMS Contemporary Mathematics Vol. 280/281, AMS, Providence, RI.   
D.A. Bini, V. Mehrmann, V. Olshevsky, E.E. Tyrtyshnikov, and M. Van Barel (eds.) (2010). Structured Matrices and Applications–The Georg Heinig Memorial Volume, Birkhauser-Springer, Basel, Switzerland.

Papers concerned with the development of fast stable solvers for structured matrices include:

T. Kailath, S. Kung, and M. Morf (1979). “Displacement Ranks of Matrices and Linear Equations,” J. Math. Anal. Applic. 68, 395–407.   
J. Chun and T. Kailath (1991). “Displacement Structure for Hankel, Vandermonde, and Related Matrices,” Lin. Alg. Applic. 151, 199–227.   
T. Kailath and A.H. Sayed (1995). “Displacement Structure: Theory and Applications,” SIAM Review 37, 297–386.   
I. Gohberg, T. Kailath, and V. Olshevsky (1995). “Fast Gaussian Elimination with Partial Pivoting for Matrices with Displacement Structure,” Math. Comput. 212, 1557–1576.   
T. Kailath and V. Olshevsky (1997). “Displacement-Structure Approach to Polynomial Vandermonde and Related Matrices,” Lin. Alg. Applic. 261, 49–90.   
G. Heinig (1997). “Matrices with Higher-Order Displacement Structure,” Lin. Alg. Applic. 278, 295–301.   
M. Gu (1998). “Stable and Efficient Algorithms for Structured Systems of Linear Systems,” SIAM J. Matrix Anal. Applic. 19, 279–306.   
S. Chandrasekaran, M. Gu, X. Sun, J. Xia, and J. Zhu (2007). “A Superfast Algorithm for Toeplitz Systems of Linear Equations,” SIAM J. Matrix Anal. Applic. 29, 1247–1266.

Displacement rank ideas can be extended to least squares problems:

R.H. Chan, J.G. Nagy, and R.J. Plemmons (1994). “Displacement Preconditioner for Toeplitz Least Squares Iterations,” ETNA 2, 44–56.   
M. Gu (1998). “New Fast Algorithms for Structured Linear Least Squares Problems,” SIAM J. Matrix Anal. Applic. 20, 244–269.   
G. Rodriguez (2006). “Fast Solution of Toeplitz- and Cauchy-Like Least-Squares Problems,” SIAM J. Matrix Anal. Applic. 28, 724–748.

For insight into the application low-displacement-rank preconditioners, see:

I. Gohberg and V. Olshevsky (1994). “Complexity of Multiplication with Vectors for Structured Matrices,” Linear Alg. Applic. 202, 163–192.   
M.E. Kilmer and D.P. O’Leary (1999). “Pivoted Cauchy-like Preconditioners for Regularized Solution of Ill-Posed Problems,” SIAM J. Sci. Comput. 21, 88–110.   
T. Kailath and V. Olshevsky (2005). “Displacement Structure Approach to Discrete-Trigonometric-Transform Based Preconditioners of G. Strang Type and of T. Chan Type,” SIAM J. Matrix Anal. Applic. 26, 706–734.
