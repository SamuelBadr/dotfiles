# 10.5 Krylov Methods for Unsymmetric Problems

If A is not symmetric, then the orthogonal tridiagonalization $Q ^ { T } A Q = T$ does not exist in general. There are two ways to proceed. The Arnoldi approach involves the column-by-column generation of an orthogonal Q such that $Q ^ { T } A \bar { Q } = H$ is the Hessenberg reduction of $\ S 7 . 4$ . The unsymmetric Lanczos approach computes the columns of matrices Q and P so that $P ^ { T } A Q = T$ is tridiagonal and $P ^ { T } Q = I$ . Methods based on these ideas that are suitable for large, sparse, unsymmetric eigenvalue problems are discussed in this section.

# 10.5.1 The Basic Arnoldi Process

One way to extend the Lanczos process to unsymmetric matrices is due to Arnoldi (1951) and revolves around the Hessenberg reduction $Q ^ { T } A Q = H$ . In particular, if $Q = [  q _ { 1 } | \cdots | q _ { n } ]$ and we compare columns in $A Q = Q H$ , then

$$
A q _ {k} = \sum_ {i = 1} ^ {k + 1} h _ {i k} q _ {i}, \quad 1 \leq k \leq n - 1.
$$

Isolating the last term in the summation gives

$$
h _ {k + 1, k} q _ {k + 1} = A q _ {k} - \sum_ {i = 1} ^ {k} h _ {i k} q _ {i} \equiv r _ {k}
$$

where $h _ { i k } = q _ { i } ^ { T } A q _ { k }$ for $i = 1 { : } k$ . It follows that if $r _ { k } \neq 0$ , then $q _ { k + 1 }$ is specified by

$$
q _ {k + 1} = r _ {k} / h _ {k + 1, k}
$$

where $h _ { k + 1 , k } = \parallel r _ { k } \parallel _ { 2 }$ . These equations define the Arnoldi process and in strict analogy to the symmetric Lanczos process (Algorithm 10.1.1) we obtain the following.

Algorithm 10.5.1 (Arnoldi Process) If $A \in \mathbb { R } ^ { n \times n }$ and $q _ { 1 } \in \mathbb { R } ^ { n }$ has unit 2-norm, then the following algorithm computes a matrix $Q _ { t } = [ q _ { 1 } , \dots , q _ { t } ] \in \mathbb { R } ^ { n \times t }$ with orthonormal columns and an upper Hessenberg matrix $H _ { t } = ( h _ { i j } ) \in \mathbb { R } ^ { t \times t }$ with the property that $A Q _ { t } = Q _ { t } H _ { t }$ . The integer t satisfies $1 \leq t \leq n$ .

$$
k = 0, r _ {0} = q _ {1}, h _ {1 0} = 1
$$

while $( h _ { k + 1 , k } \neq 0 )$

$$
q _ {k + 1} = r _ {k} / h _ {k + 1, k}
$$

$$
k = k + 1
$$

$$
r _ {k} = A q _ {k}
$$

for i = 1:k

$$
h _ {i k} = q _ {i} ^ {T} r _ {k}
$$

$$
r _ {k} = r _ {k} - h _ {i k} q _ {i}
$$

end

$$
h _ {k + 1, k} = \left\| r _ {k} \right\| _ {2}
$$

end

$$
t = k
$$

The $q _ { k }$ are called Arnoldi vectors and they define an orthonormal basis for the Krylov subspace $\kappa ( A , q _ { 1 } , k )$ :

$$
\operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\} = \operatorname{span} \left\{q _ {1}, A q _ {1}, \dots , A ^ {k - 1} q _ {1} \right\}. \tag {10.5.1}
$$

The situation after k steps is summarized by the equation

$$
A Q _ {k} = Q _ {k} H _ {k} + r _ {k} e _ {k} ^ {T} \tag {10.5.2}
$$

where $Q _ { k } = \left[ q _ { 1 } | \cdots | q _ { k } \right] , e _ { k } = I _ { k } ( : , k )$ , and

$$
H _ {k} = \left[ \begin{array}{c c c c c} h _ {1 1} & h _ {1 2} & \dots & \dots & h _ {1 k} \\ h _ {2 1} & h _ {2 2} & \dots & \dots & h _ {2 k} \\ 0 & h _ {3 2} & \ddots & & \vdots \\ \vdots & & \ddots & \ddots & \vdots \\ 0 & \dots & \dots & h _ {k, k - 1} & h _ {k k} \end{array} \right].
$$

Any decomposition of the form (10.5.2) is a k-step Arnoldi decomposition if $Q _ { k } \in \mathbb { R } ^ { n \times k }$ has orthonormal columns, $H _ { k } \in \mathbb { R } ^ { k \times k }$ is upper Hessenberg, and $Q _ { k } ^ { T } r _ { k } = 0 .$ .

If $\boldsymbol { y } \in \mathbb { R } ^ { k }$ is a unit 2-norm eigenvector for $H _ { k }$ and $H _ { k } y = \lambda y$ , then from (10.5.2)

$$
(A - \lambda I) x = (e _ {k} ^ {T} y) r _ {k}
$$

where $x = Q _ { k } y$ . Since $r _ { k } \in \mathcal { K } ( A , q _ { 1 } , k ) ^ { \perp }$ , it follows that $( \lambda , x )$ is a Ritz pair for A with respect to $\kappa ( A , q _ { 1 } , k )$ . Note that if $v = ( e _ { k } ^ { T } y ) r _ { k }$ , then

$$
(A + E) x = \lambda x
$$

where $E = - v x ^ { T }$ with $\left\| \mathbf { \nabla } E \right\| _ { 2 } = \left| y _ { k } \right| \left\| r _ { k } \right\| _ { 2 }$ . Recall that in the unsymmetric case, computing an eigenvalue of a nearby matrix does not mean that it is close to an exact eigenvalue.

Some numerical properties of the Arnoldi iteration are discussed by Wilkinson (AEP, p. 382). The history of practical Arnoldi-based eigensolvers begins with Saad (1980). Two features of the method distinguish it from the symmetric Lanczos process:

• Arnoldi vectors $q _ { 1 } , \ldots , q _ { k }$ must all be referenced in step k and the computation of $q _ { k + 1 }$ involves $O ( k n )$ flops excluding the matrix-vector product $A q _ { k }$ Thus, there is a steep penalty associated with the generation of long Arnoldi sequences.   
• Extremal eigenvalue information is not as forthcoming as in the symmetric case. There is no unsymmetric Kaniel-Paige-Saad convergence theory.

These realities suggest a framework in which we use the Arnoldi iteration idea with repeated, carefully chosen restarts and a controlled iteration maximum. We described such a framework in conjunction with the block Lanczos procedure in §10.3.7.

# 10.5.2 Arnoldi with Restarting

Consider running Arnoldi for m steps and then restarting the iteration with a new initial vector $q _ { + }$ chosen from the span of the Arnoldi vectors $q _ { 1 } , \ldots , q _ { m }$ . Because of the Krylov connection (10.5.1), $q _ { + }$ has the form

$$
q _ {+} = p (A) q _ {1}
$$

for some polynomial of degree $m - 1$ . It is instructive to examine the action of $p ( A )$ in terms of A’s eigenvalues and eigenvectors. Assume for clarity that $A \in \mathbb { R } ^ { n \times n }$ is diagonalizable and that $A z _ { i } = \lambda _ { i } z _ { i }$ for $i = 1 { : } n$ . If $q _ { 1 }$ has the eigenvector expansion

$$
q _ {1} = a _ {1} z _ {1} + \dots + a _ {n} z _ {n},
$$

then $q _ { + }$ is a scalar multiple of

$$
z = a _ {1} p (\lambda_ {1}) z _ {1} + \dots + a _ {n} p (\lambda_ {n}) z _ {n}.
$$

Note that if $p ( \lambda _ { \alpha } ) \gg p ( \lambda _ { \beta } )$ , then relatively speaking, $q _ { + }$ is much richer in the direction of $z _ { \alpha }$ than in the direction of $z _ { \beta }$ . More generally, by carefully choosing $p ( \lambda )$ we can design $q _ { + }$ so that its component in certain eigenvector directions is emphasized while its component in other eigenvector directions is deemphasized. For example, if

$$
p (\lambda) = c \cdot (\lambda - \mu_ {1}) (\lambda - \mu_ {2}) \dots (\lambda - \mu_ {p}) \tag {10.5.3}
$$

where $c$ is a constant, then $q _ { + }$ is a unit vector in the direction of

$$
z = c \cdot \sum_ {k = 1} ^ {n} a _ {k} \left(\prod_ {i = 1} ^ {p} \left(\lambda_ {k} - \mu_ {i}\right)\right) z _ {k}.
$$

It follows that $z _ { \beta }$ is deemphasized relative to $z _ { \alpha }$ if $\lambda _ { \beta }$ is near to one of the “filter values” $\mu _ { 1 } , \ldots , \mu _ { p }$ and $\lambda _ { \alpha }$ is not. Thus, the act of picking a good restart vector $q _ { + }$ from $\mathcal { K } ( A , q _ { 1 } , m )$ is the act of picking a filter polynomial that tunes out unwanted portions of the spectrum. Various heuristics for doing this have been developed based on computed Ritz vectors. See Saad (1980, 1984, 1992).

# 10.5.3 Implicit Restarting

We describe an Arnoldi restarting procedure due to Sorensen (1992) that implicitly determines the filter polynomial (10.5.3) using the QR iteration with shifts. (See §7.5.2.) Suppose $H _ { c } \in \mathbb { R } ^ { m \times m }$ is upper Hessenberg, $\mu _ { 1 } , \ldots , \mu _ { p }$ are scalars, and the matrix $H _ { + }$ is obtained via the shifted QR iteration:

$$
H ^ {(0)} = H _ {c}
$$

for $i = 0 { : } p$

$$
H ^ {(i - 1)} - \mu_ {i} I = V _ {i} R _ {i} \quad \text {(Givens QR)} \tag {10.5.4}
$$

$$
H ^ {(i)} = R _ {i} V _ {i} + \mu_ {i} I
$$

end

$$
H _ {+} = H ^ {(p)}
$$

Recall from §7.4.2 that each $H ^ { ( i ) }$ is upper Hessenberg. Moreover, if

$$
V = V _ {1} \dots V _ {p}, \tag {10.5.5}
$$

then

$$
H _ {+} = V ^ {T} H _ {c} V. \tag {10.5.6}
$$

The following result shows that the filter polynomial (10.5.3) has a relationship to (10.5.4).

Theorem 10.5.1. $I f V = V _ { 1 } \cdots V _ { p }$ and $R = R _ { p } \cdot \cdot \cdot R _ { 1 }$ are defined by $( 1 0 . 5 . 4 )$ , then

$$
V R = \left(H _ {c} - \mu_ {1} I\right) \dots \left(H _ {c} - \mu_ {p} I\right). \tag {10.5.7}
$$

Proof. We use induction, noting that the theorem is obviously true if $p = 1$ . If $\tilde { V } = V _ { 1 } \cdot \cdot \cdot V _ { p - 1 }$ and $\tilde { R } = R _ { p - 1 } \cdot \cdot \cdot R _ { 1 }$ , then

$$
\begin{array}{l} V R = \tilde {V} (V _ {p} R _ {p}) \tilde {R} = \tilde {V} (H ^ {(p - 1)} - \mu_ {p} I) \tilde {R} = \tilde {V} (\tilde {V} ^ {T} H _ {c} \tilde {V} - \mu_ {p} I) \tilde {R} \\ = \left(H _ {c} - \mu_ {p} I\right) \tilde {V} \tilde {R} = \left(H _ {c} - \mu_ {p} I\right) \left(H _ {c} - \mu_ {1} I\right) \dots \left(H _ {c} - \mu_ {p - 1} I\right), \\ \end{array}
$$

where we used the fact that $H ^ { ( p - 1 ) } = \tilde { V } ^ { T } H _ { c } \tilde { V }$ .

Note that the matrix R in (10.5.7) is upper triangular and so it follows that

$$
V (:, 1) = p (H _ {c}) e _ {1}
$$

where $p ( \lambda )$ is the filter polynomial (10.5.3) with $c = 1 / R ( 1 , 1 )$ .

Now suppose that we have performed m steps of the Arnoldi iteration with starting vector $q _ { 1 }$ . The Arnoldi factorization (10.5.2) says that we have an upper Hessenberg matrix $H _ { c } \in \mathbb { R } ^ { m \times m }$ and a matrix $Q _ { c } \in \mathbb { R } ^ { n \times m }$ with orthonormal columns such that

$$
A Q _ {c} = Q _ {c} H _ {c} + r _ {c} e _ {m} ^ {T}. \tag {10.5.8}
$$

Note that $Q _ { c } ( : , 1 ) = q _ { 1 }$ and $r _ { c } \in \mathbb { R } ^ { n }$ has the property that $Q _ { c } ^ { T } r _ { c } = 0$ . If we apply (10.5.4) to $H _ { c }$ , then by using (10.5.5) and (10.5.6) the preceding Arnoldi factorization transforms to

$$
A Q _ {+} = Q _ {+} H _ {+} + r _ {c} e _ {m} ^ {T} V \tag {10.5.9}
$$

where

$$
Q _ {+} = Q _ {c} V.
$$

If $q _ { + }$ is the first column of this matrix, then

$$
q _ {+} = Q _ {+} (:, 1) = Q _ {c} V (:, 1) = c \cdot Q _ {c} \left(H _ {c} - \mu_ {1} I\right) \dots \left(H _ {c} - \mu_ {p} I\right) e _ {1}.
$$

Equation (10.5.8) implies that

$$
(A - \mu I) Q _ {c} e _ {1} = Q _ {c} (H _ {c} - \mu I) e _ {1}
$$

for any $\mu \in \mathbb { R }$ and so

$$
q _ {+} = c (A - \mu_ {1} I) \dots (A - \mu_ {p} I) Q _ {c} e _ {1} = p (A) q _ {1}.
$$

This suggests the following framework for repeated restarting:

# Repeat:

With starting vector $q _ { 1 }$ , perform m steps of the Arnoldi iteration obtaining $Q _ { c } \in \mathbb { R } ^ { n \times m }$ and $H _ { c } \in \mathbb { R } ^ { m \times m }$ .

Determine filter values µ1, . . . , µp . $\mu _ { 1 } , \ldots , \mu _ { p }$ (10.5.10)

Perform $p$ steps of the shifted QR iteration (10.5.4) obtaining the Hessenberg matrix $H _ { + }$ and the orthogonal matrix $V$ .

Replace $q _ { 1 }$ with the first column of $Q _ { c } V$ .

However, we can do better than this. The orthogonal matrices $V _ { 1 } , \ldots , V _ { p }$ that arise in (10.5.4) are each upper Hessenberg. (This is easily deduced from the structure of the Givens rotations in Algorithm 5.2.5.) Thus, V has lower bandwidth $p$ and so $V ( m , 1 : m - p - 1 ) = 0$ . It follows from (10.5.9) that if $j = m - p ,$ then

$$
A Q _ {+} (:, 1: j) = Q _ {+} (:, 1: j) H _ {+} (1: j, 1: j) + v _ {m j} r _ {c} e _ {j}
$$

is a j-step Arnoldi decomposition. In other words, we are all set to perform step $j + 1$ of the Arnoldi iteration with starting vector $q _ { + }$ . There is no need to launch the restart from step 1. This leads us to the following modification of (10.5.10):

With starting vector $q _ { 1 }$ , perform m steps of the Arnoldi iteration obtaining

$$
Q _ {c} \in \mathbb {R} ^ {n \times m}, H _ {c} \in \mathbb {R} ^ {m \times m}, \text {   and   } r _ {c} \in \mathbb {R} ^ {n} \text {   so   } A Q _ {c} = Q _ {c} H _ {c} + r _ {c} e _ {m} ^ {T}.
$$

# Repeat:

Determine filter values $\mu _ { 1 } , \ldots , \mu _ { p }$ .

Perform $p$ steps of the shifted QR iteration (10.5.4) applied to $H _ { c }$ obtaining $H _ { + } \in \mathbb { R } ^ { m \times m }$ and $V = ( v _ { i j } ) \in \mathbb { R } ^ { m \times m }$ .

Replace $Q _ { c }$ with the first j columns of $Q _ { c } V$ .

Replace $H _ { c }$ with $H _ { + } ( 1 { : } j , 1 { : } j )$ . .

Replace $r _ { c }$ with $v _ { m j } r _ { c }$ .

Starting with $A Q _ { c } = Q _ { c } H _ { c } + r _ { c } e _ { j } ^ { T }$ , perform steps $j + 1 , \ldots , j + p = m$ of the Arnoldi iteration obtaining $A Q _ { m } = Q _ { m } H _ { m } + r _ { m } e _ { m } ^ { T }$ .

Set $Q _ { c } = Q _ { m } , H _ { c } = H _ { m }$ , and $r _ { c } = r _ { m }$ .

In light of our remarks in §10.5.2, the filter values $\mu _ { 1 } , \ldots , \mu _ { p }$ should be chosen in the vicinity of $A \mathrm { { } s }$ “unwanted” eigenvalues. In this regard it is possible to formulate useful heuristics that are based on the eigenvalues of the m-by-m Hessenberg matrix $H _ { + }$ . For example, suppose the goal is to find the three smallest eigenvalues of A in absolute value. If $p = m - 3$ and $\lambda ( H _ { + } ) = \{ \tilde { \lambda } _ { 1 } , . . . , \tilde { \lambda } _ { m } \}$ with $| \tilde { \lambda } _ { 1 } \bar { | } \geq \cdots \geq | \tilde { \lambda } _ { m } |$ , then it is reasonable to set $\mu _ { i } = \tilde { \lambda } _ { i }$ for $i = 1 { : } p$ .

The Arnoldi iteration with implicit restarts has many attractive attributes. For implementation details and further analysis, see Lehoucq and Sorensen (1996), Morgan (1996), and the ARPACK manual by Lehoucq, Sorensen, and Yang (1998).

# 10.5.4 The Krylov-Schur Algorithm

An alternative restart procedure due to Stewart (2001) relies upon a carefully ordered Schur decomposition of the Hessenberg matrix $H _ { m }$ that is produced after m steps of the Arnoldi iteration. Suppose we have computed

$$
A Q _ {m} = Q _ {m} H _ {m} + r _ {m} e _ {m} ^ {T}
$$

and that $m = j + p ;$ , where $j$ is the number of A’s eigenvalues that we wish to compute. Let

$$
U ^ {T} H _ {m} U = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right]
$$

be the Schur decomposition of A and assume that the eigenvalues have been ordered so that the eigenvalues of $T _ { 1 1 } \in \mathbb { R } ^ { j \times j }$ are of interest and the eigenvalues of $T _ { 2 2 } \in \mathbb { R } ^ { p \times p }$ are not. (For clarity we ignore the possibility of complex eigenvalues.) The Arnoldi decomposition above transforms to

$$
A Q _ {+} = Q _ {+} T + r _ {c} e _ {m} ^ {T} U
$$

where $Q _ { + } = Q _ { m } U$ . It follows that

$$
A Q _ {+} (:, 1: j) = Q _ {+} (:, 1: j) T _ {1 1} + r _ {m} u ^ {T}
$$

where $u ^ { T } = U ( m , 1 { : } j )$ . It is possible to determine an orthogonal $\boldsymbol { Z } \in \mathbb { R } ^ { j \times j }$ so that $Z ^ { T } T _ { 1 1 } Z$ is upper Hessenberg and $Z ^ { T } u = \tau e _ { j }$ . (See P10.5.2.) It follows that

$$
A \left(Q _ {+} Z\right) = \left(Q _ {+} Z\right) \left(Z ^ {T} T _ {1 1} Z\right) + r _ {c} \left(Z ^ {T} u\right) ^ {T}
$$

is a j-step Arnoldi factorization. We then set $Q _ { j } , H _ { j }$ and $r _ { j }$ to be $Q _ { + } Z , Z ^ { T } T _ { 1 1 } Z$ , and $\tau r _ { m }$ respectively and perform Arnoldi steps $j + 1$ through $j + p = m$ . For more detailed discussion, see Stewart (MAE, Chap. 5) and Watkins (FMC, Chap. 9).

# 10.5.5 Unsymmetric Lanczos Tridiagonalization

Another way to extend the symmetric Lanczos process is to reduce A to tridiagonal form using a general similarity transformation. Suppose $A \in \mathbb { R } ^ { n \times n }$ and that a nonsingular matrix $Q$ exists such that

$$
Q ^ {- 1} A Q   =   T   =   \left[ \begin{array}{c c c c c} \alpha_ {1} & \gamma_ {1} & & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \gamma_ {n - 1} \\ 0 & \dots & & \beta_ {n - 1} & \alpha_ {n} \end{array} \right]  .
$$

With the column partitionings

$$
Q = \left[ q _ {1} \mid \dots \mid q _ {n} \right],
$$

$$
Q ^ {- T} = \tilde {Q} = \left[ \tilde {q} _ {1} | \dots | \tilde {q} _ {n} \right],
$$

we find upon comparing columns in $A Q { = } Q T$ and $A ^ { T } \tilde { Q } = \tilde { Q } T ^ { T }$ that

$$
A q _ {k} = \gamma_ {k - 1} q _ {k - 1} + \alpha_ {k} q _ {k} + \beta_ {k} q _ {k + 1}, \quad \gamma_ {0} q _ {0} \equiv 0,
$$

$$
A ^ {T} \tilde {q} _ {k} = \beta_ {k - 1} \tilde {q} _ {k - 1} + \alpha_ {k} \tilde {q} _ {k} + \gamma_ {k} \tilde {q} _ {k + 1}, \quad \beta_ {0} \tilde {q} _ {0} \equiv 0,
$$

for $k = 1 { : } n - 1$ . These equations together with the biorthogonality condition

$$
\tilde {Q} ^ {T} Q = I _ {n}
$$

imply

$$
\alpha_ {k} = \tilde {q} _ {k} ^ {T} A q _ {k}
$$

and

$$
\begin{array}{l} \beta_ {k} q _ {k + 1} \equiv r _ {k} = (A - \alpha_ {k} I) q _ {k} - \gamma_ {k - 1} q _ {k - 1}, \\ \gamma_ {k} \tilde {q} _ {k + 1} \equiv \tilde {r} _ {k} = (A - \alpha_ {k} I) ^ {T} \tilde {q} _ {k} - \beta_ {k - 1} \tilde {q} _ {k - 1}. \\ \end{array}
$$

There is some flexibility in choosing the scale factors $\beta _ { k }$ and $\gamma _ { k }$ . Note that

$$
1 = \tilde {q} _ {k + 1} ^ {T} q _ {k + 1} = \left(\tilde {r} _ {k} / \gamma_ {k}\right) ^ {T} \left(r _ {k} / \beta_ {k}\right).
$$

It follows that once $\beta _ { k }$ is specified, then $\gamma _ { k }$ is given by

$$
\gamma_ {k} = \tilde {r} _ {k} ^ {T} r _ {k} / \beta_ {k}.
$$

With the “canonical” choice $\beta _ { k } = \parallel r _ { k } \parallel _ { 2 }$ we obtain

$q _ { 1 } , \tilde { q } _ { 1 }$ given unit 2-norm vectors with $\tilde { q } _ { 1 } ^ { T } q _ { 1 } \neq 0$

$$
k = 0, q _ {0} = 0, r _ {0} = q _ {1}, \tilde {q} _ {0} = 0, s _ {0} = \tilde {q} _ {1}
$$

while $( r _ { k } \neq 0 )$ and $( \tilde { r } _ { k } \neq 0 )$ and $( \tilde { r } _ { k } ^ { T } r _ { k } \neq 0 )$

$$
\begin{array}{l} \beta_ {k} = \left\| r _ {k} \right\| _ {2} \\ \gamma_ {k} = \tilde {r} _ {k} ^ {T} r _ {k} / \beta_ {k} \\ q _ {k + 1} = r _ {k} / \beta_ {k} \\ \tilde {q} _ {k + 1} = \tilde {r} _ {k} / \gamma_ {k} \\ k = k + 1 \tag {10.5.11} \\ \end{array}
$$

$$
\alpha_ {k} = \tilde {q} _ {k} ^ {T} A q _ {k}
$$

$$
r _ {k} = (A - \alpha_ {k} I) q _ {k} - \gamma_ {k - 1} q _ {k - 1}
$$

$$
\tilde {r} _ {k} = (A - \alpha_ {k} I) ^ {T} \tilde {q} _ {k} - \beta_ {k - 1} \tilde {q} _ {k - 1}
$$

end

If

$$
T _ {k} = \left[ \begin{array}{c c c c c} \alpha_ {1} & \gamma_ {1} & & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \gamma_ {k - 1} \\ 0 & \dots & & \beta_ {k - 1} & \alpha_ {k} \end{array} \right],
$$

then the situation at the bottom of the loop is summarized by the equations

$$
A \left[ q _ {1} \mid \dots \mid q _ {k} \right] = \left[ q _ {1} \mid \dots \mid q _ {k} \right] T _ {k} + r _ {k} e _ {k} ^ {T}, \tag {10.5.12}
$$

$$
A ^ {T} \left[ \tilde {q} _ {1} \mid \dots \mid \tilde {q} _ {k} \right] = \left[ \tilde {q} _ {1} \mid \dots \mid \tilde {q} _ {k} \right] T _ {k} ^ {T} + \tilde {r} _ {k} e _ {k} ^ {T}. \tag {10.5.13}
$$

If $r _ { k } = 0 .$ , then the iteration terminates and span $\{ q _ { 1 } , \ldots , q _ { k } \}$ is an invariant subspace for A. If $\tilde { r } _ { k } = 0$ , then the iteration also terminates and span $\{ \tilde { q } _ { 1 } , \dots , \tilde { q } _ { k } \}$ is an invariant subspace for $A ^ { T }$ . However, if neither of these conditions is true and $\tilde { r } _ { k } ^ { T } r _ { k } = 0$ , then the tridiagonalization process ends without any invariant subspace information. This is called serious breakdown. See Wilkinson (AEP, p. 389) for an early discussion of the matter.

# 10.5.6 The Look-Ahead Idea

It is interesting to examine the serious breakdown issue in the block version of (10.5.11). For clarity assume that $A \in \mathbb { R } ^ { n \times n }$ with $n = r p$ . Consider the factorization in which we want $\tilde { Q } ^ { T } \dot { Q } = I _ { n }$ :

$$
\tilde {Q} ^ {T} A Q = \left[ \begin{array}{c c c c c} M _ {1} & C _ {1} ^ {T} & & \dots & 0 \\ B _ {1} & M _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & C _ {r - 1} ^ {T} \\ 0 & \dots & & B _ {r - 1} & M _ {r} \end{array} \right] \tag {10.5.14}
$$

where all the blocks are p-by-p. Let $Q = [  Q _ { 1 } | \cdot \cdot \cdot | Q _ { r } ]$ and $\tilde { Q } = [ \tilde { Q } _ { 1 } | \cdots | \tilde { Q } _ { r } ]$ be conformable partitionings of Q and Q˜. Comparing block columns in the equations $A Q = Q T$ and $\mathring { A } ^ { T } \tilde { Q } = \tilde { Q } T ^ { T }$ , we obtain

$$
Q _ {k + 1} B _ {k} = A Q _ {k} - Q _ {k} M _ {k} - Q _ {k - 1} C _ {k - 1} ^ {T} \equiv R _ {k},
$$

$$
\tilde {Q} _ {k + 1} C _ {k} = A ^ {T} \tilde {Q} _ {k} - \tilde {Q} _ {k} M _ {k} ^ {T} - \tilde {Q} _ {k - 1} B _ {k - 1} ^ {T} \equiv S _ {k}.
$$

Note that

$$
M _ {k} = \tilde {Q} _ {k} ^ {T} A Q _ {k}.
$$

If $S _ { k } ^ { T } R _ { k } = C _ { k } ^ { T } \tilde { Q } _ { k + 1 } ^ { T } Q _ { k + 1 } B _ { k } \in \mathbb { R } ^ { p \times p }$ is nonsingular and we compute $B _ { k } , C _ { k } \in \mathbb { R } ^ { p \times p }$ so that

$$
C _ {k} ^ {T} B _ {k} = S _ {k} ^ {T} R _ {k},
$$

then

$$
Q _ {k + 1} = R _ {k} B _ {k} ^ {- 1}, \tag {10.5.15}
$$

$$
\tilde {Q} _ {k + 1} = S _ {k} C _ {k} ^ {- 1} \tag {10.5.16}
$$

satisfy $\tilde { Q } _ { k + 1 } ^ { T } Q _ { k + 1 } = I _ { p }$ . Serious breakdown in this setting is associated with having a singular $S _ { k } ^ { T } R _ { k }$ .

One way of solving the serious breakdown problem in (10.5.11) is to go after a factorization of the form (10.5.14) in which the block sizes are dynamically determined. Roughly speaking, in this approach matrices $Q _ { k + 1 }$ and $\tilde { Q } _ { k + 1 }$ are built up column by column with special recursions that culminate in the production of a nonsingular $\tilde { Q } _ { k + 1 } ^ { T } Q _ { k + 1 }$ . The computations are arranged so that the biorthogonality conditions $\tilde { Q } _ { i } ^ { T } Q _ { k + 1 } = 0$ and $Q _ { i } ^ { T } \tilde { Q } _ { k + 1 } = 0$ hold for $i = 1 { : } k$ .

A method of this form belongs to the family of look-ahead Lanczos methods. The length of a look-ahead step is the width of the $Q _ { k + 1 }$ and $\tilde { Q } _ { k + 1 }$ that it produces. If that width is one, a conventional block Lanczos step may be taken. Length-2 lookahead steps are discussed in Parlett, Taylor, and Liu (1985). The notion of incurable breakdown is also presented by these authors. Freund, Gutknecht, and Nachtigal (1993) cover the general case along with a host of implementation details. Floating point considerations require the handling of “near” serious breakdown. In practice, each $M _ { k }$ that is 2-by-2 or larger corresponds to an instance of near serious breakdown.

# Problems

P10.5.1 Recalling how Theorem 10.1.1 establishes the orthogonality of the Lanczos vectors in Algorithm 10.1.1, state and prove an analogous theorem that does the same thing for the Arnoldi vectors in Algorithm 10.5.1.

P10.5.2 Show that if $C \in \mathbb { R } ^ { j \times j }$ and $u \in \mathbb { R } ^ { j }$ , then there exists an orthogonal $Z \in \mathbb { R } ^ { n \times n }$ so that $Z ^ { T } A Z \ = \ H$ is upper Hessenberg and the last column of Z is a multiple of u. Hint: Compute a Householder matrix P so that $P u$ is a multiple of $e _ { j }$ . Then reduce $C = \mathsf { \bar { P } } ^ { T } C P$ to upper Hessenberg form by producing a sequence of Householder updates $C = P _ { i } ^ { T } C P$ where $C ( n - i + 1 , 1 { : } n - i - 1 )$ is zeroed, $i = 1 { : } n - 2$ .

P10.5.3 Give an example of a starting vector for which the unsymmetric Lanczos iteration (10.5.11) breaks down without rendering any invariant subspace information. Use

$$
A = \left[ \begin{array}{c c c} 1 & 6 & 2 \\ 3 & 0 & 2 \\ 1 & 3 & 5 \end{array} \right].
$$

P10.5.4 Suppose $H \in \mathbb { R } ^ { n \times n }$ is upper Hessenberg. Discuss the computation of a unit upper triangular matrix U such that $H U = U T$ where T is tridiagonal.

P10.5.5 Show that the QR algorithm for eigenvalues does not preserve tridiagonal structure in the unsymmetric case.

# Notes and References for §10.5

For both analysis and implementation insight, Saad (NMLE) offers the most comprehensive treatment of unsymmetric Krylov methods. Stewart (MAE) and Watkins (MEP) devote entire chapters to the subject and are highly recommended as is the following review article:

D.C. Sorensen (2002). “Numerical Methods for Large Eigenvalue Problems,” Acta Numerica 11, 519–584.

The original Arnoldi idea first appeared in:

W.E. Arnoldi (1951). “The Principle of Minimized Iterations in the Solution of the Matrix Eigenvalue Problem,” Quarterly of Applied Mathematics 9, 17–29.

Saad set the stage for the development of practical implementations, see:

Y. Saad (1980). “Variations of Arnoldi’s Method for Computing Eigenelements of Large Unsymmetric Matrices.,” Lin. Alg. Applic. 34, 269–295.

Y. Saad (1984). “Chebyshev Acceleration Techniques for Solving Nonsymmetric Eigenvalue Problems,” Math. Comput. 42, 567–588.

Y. Saad (1989). “Krylov Subspace Methods on Supercomputers,” SIAM J. Sci. Stat. Comput., 10, 1200–1232.

References for implicit restarting in the Arnoldi context include:

D.C. Sorensen (1992). “Implicit Application of Polynomial Filters in a k-Step Arnoldi Method,” SIAM J. Matrix Anal. Applic. 13, 357–385.   
R.B. Lehoucq and D.C. Sorensen (1996). “Deflation Techniques for an Implicitly Restarted Iteration,” SIAM J. Matrix Anal. Applic. 17, 789–821.   
R.B. Morgan (1996). “On Restarting the Arnoldi Method for Large Nonsymmetric Eigenvalue Problems,” Math Comput. 65, 1213–1230.   
K. Meerbergen and A. Spence (1997). “Implicitly Restarted Arnoldi with Purification for the Shift-Invert Transformation,” Math. Comput. 218, 667–689.   
R.B. Lehoucq, D. C. Sorensen, and C. Yang (1998). ARPACK Users’ Guide: Solution of Large-Scale Eigenvalue Problems with Implicitly Restarted Arnoldi Methods, SIAM Publications, Philadelphia, PA.   
A. Stathopoulos, Y. Saad, and K. Wu (1998). “Dynamic Thick Restarting of the Davidson and the Implicitly Restarted Arnoldi Methods,” SIAM J. Sci. Comput. 19, 227–245.   
R.B. Lehoucq (2001). “Implicitly Restarted Arnoldi Methods and Subspace Iteration,” SIAM J. Matrix Anal. Applic. 23, 551–562.   
The Krylov-Schur approach to Arnoldi restarting is proposed in:   
G.W. Stewart (2001). “A Krylov-Schur Algorithm for Large Eigenproblems,” SIAM J. Matrix Anal. Applic., 23, 601–614.   
The rational Arnoldi process involves the shift-and-invert idea. In this framework Arnoldi is applied to the matrix $( A - \mu { \dot { I } } ) ^ { - 1 }$ , see:   
A. Ruhe (1984). “Rational Krylov Algorithms for Eigenvalue Computation,” Lin. Alg. Applic. 58, 391–405.   
A. Ruhe (1994). “Rational Krylov Algorithms for Nonsymmetric Eigenvalue Problems II. Matrix Pairs,” Lin. Alg. Applic. 197, 283–295.   
A. Ruhe (1994). “The Rational Krylov Algorithm for Nonsymmetric Eigenvalue Problems III: Complex Shifts for Real Matrices,” BIT 34, 165–176.   
Matrix function problems that involve large sparse matrices can be addressed using Krylov/Lanczos ideas, see:   
Y. Saad (1992). “Analysis of Some Krylov Subspace Approximations to the Matrix Exponential,” SIAM J. Numer. Anal. 29, 209–228.   
M. Hochbruck and C. Lubich (1997). “On Krylov Subspace Approximations to the Matrix Exponential Operator,” SIAM J. Numer. Anal. 34, 1911–1925.   
V. Druskin, A. Greenbaum and L. Knizhnerman (1998). “Using Nonorthogonal Lanczos Vectors in the Computation of Matrix Functions,” SIAM J. Sci. Comput. 19, 38–54.   
N. Del Buono, L. Lopez, and R. Peluso (2005). “Computation of the Exponential of Large Sparse Skew–Symmetric Matrices,” SIAM J. Sci. Comp. 27, 278–293.   
M. Eiermann and O.G. Ernst (2006). “A Restarted Krylov Subspace Method for the Evaluation of Matrix Functions,” SIAM J. Numer. Anal. 44, 2481–2504.   
J. van den Eshof and M. Hochbruck (2006). “Preconditioning Lanczos Approximations to the Matrix Exponential,” SIAM J. Sci. Comput. 27, 1438–1457.   
Other Arnoldi-related papers include:   
T. Huckle (1994). “The Arnoldi Method for Normal Matrices,” SIAM J. Matrix Anal. Applic. 15, 479–489.   
K.C. Toh and L.N. Trefethen (1996). “Calculation of Pseudospectra by the Arnoldi Iteration,” SIAM J. Sci. Comput. 17, 1–15.   
T.G. Wright and L.N. Trefethen (2001). “Large–Scale Computation of Pseudospectra Using ARPACK and Eigs,” SIAM J. Sci. Comput. 23, 591–605.   
V. Hernandez, J.E. Roman, and A. Tomas (2007). “Parallel Arnoldi Eigensolvers with Enhanced Scalability via Global Communications Rearrangement,” Parallel Computing 33, 521–540.   
The unsymmetric Lanczos process and related look ahead ideas are nicely presented in:   
B.N. Parlett, D. Taylor, and Z. Liu (1985). “A Look-Ahead Lanczos Algorithm for Unsymmetric Matrices,” Math. Comput. 44, 105–124.   
R.W. Freund, M. Gutknecht, and N. Nachtigal (1993). “An Implementation of the Look-Ahead Lanczos Algorithm for Non-Hermitian Matrices,” SIAM J. Sci. Stat. Comput. 14, 137–158.

# See also:

Y. Saad (1982). “The Lanczos Biorthogonalization Algorithm and Other Oblique Projection Methods for Solving Large Unsymmetric Eigenproblems,” SIAM J. Numer. Anal. 19, 485–506.   
D.L. Boley, S. Elhay, G.H. Golub and M.H. Gutknecht (1991) “Nonsymmetric Lanczos and Finding Orthogonal Polynomials Associated with Indefinite Weights,” Numer. Algorithms 1, 21–43.   
G.A. Geist (1991). “Reduction of a General Matrix to Tridiagonal Form,” SIAM J. Matrix Anal. Applic. 12, 362–373.   
C. Brezinski, M. Zaglia, and H. Sadok (1991). “Avoiding Breakdown and Near Breakdown in Lanczos Type Algorithms,” Numer. Algorithms 1, 261–284.   
S.K. Kim and A.T. Chronopoulos (1991). “A Class of Lanczos-Like Algorithms Implemented on Parallel Computers,” Parallel Comput. 17, 763–778.   
B.N. Parlett (1992). “Reduction to Tridiagonal Form and Minimal Realizations,” SIAM J. Matrix Anal. Applic. 13, 567–593.   
M. Gutknecht (1992). “A Completed Theory of the Unsymmetric Lanczos Process and Related Algorithms, Part I,” SIAM J. Matrix Anal. Applic. 13, 594–639.   
M. Gutknecht (1994). “A Completed Theory of the Unsymmetric Lanczos Process and Related Algorithms, Part II,” SIAM J. Matrix Anal. Applic. 15, 15–58.   
Z. Bai (1994). “Error Analysis of the Lanczos Algorithm for Nonsymmetric Eigenvalue Problem,” Math. Comput. 62, 209–226.   
T. Huckle (1995). “Low-Rank Modification of the Unsymmetric Lanczos Algorithm,” Math.Comput. 64, 1577–1588.   
Z. Jia (1995). “The Convergence of Generalized Lanczos Methods for Large Unsymmetric Eigenproblems,” SIAM J. Matrix Anal. Applic. 16, 543–562.   
M.T. Chu, R.E. Funderlic, and G.H. Golub (1995). “A Rank-One Reduction Formula and Its Applications to Matrix Factorizations,” SIAM Review 37, 512–530.

Computing eigenvalues of unsymmetric tridiagonal matrices is discussed in:

D.A. Bini, L. Gemignani, and F. Tisseur (2005). “The Ehrlich-Aberth Method for the Nonsymmetric Tridiagonal Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 27, 153–175.
