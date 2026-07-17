# 10.4 Large Sparse SVD Frameworks

The connections between the SVD problem and the symmetric eigenvalue problem are discussed in §8.6.1. In light of that discussion, it is not surprising that there is a Lanczos process for computing selected singular values and vectors of a large, sparse, rectangular matrix A. The basic idea is to generate a bidiagonal matrix B that is orthogonally equivalent to A. We show how to do this in §5.4 using Householder transformations. However, to avoid large dense submatrices along the way, the Lanczos approach generates the bidiagonal entries entries directly.

# 10.4.1 Golub-Kahan Upper Bidiagonalization

Suppose $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ and recall from §5.4.8 that there exist orthogonal $U \in \mathbb { R } ^ { m \times m }$ and $V \in \mathbb { R } ^ { n \times n }$ so that

$$
U ^ {T} A V = B = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & \dots & \dots & 0 \\ 0 & \alpha_ {2} & \beta_ {2} & \dots & \vdots \\ \vdots & \ddots & \ddots & \ddots & \vdots \\ \vdots & & 0 & \alpha_ {n - 1} & \beta_ {n - 1} \\ 0 & \dots & \dots & 0 & \alpha_ {n} \\ \hline & & 0 \end{array} \right]. \tag {10.4.1}
$$

Since A and B are orthogonally related, they have the same singular values.

Analogously to our derivation of the symmetric Lanczos procedure in §10.1.1, we proceed to outline a sparse-matrix-friendly method for determining the diagonal and superdiagonal of B. The challenge is to bypass the generally full intermediate matrices associated with the Householder bidiagonalization process (Algorithm 5.4.2). We expect to extract good singular value/vector information long before the full bidiagonalization is complete.

The key is to develop useful recipes for the α’s and $\beta \mathrm { { ^ { * } s } }$ from the matrix equations $A V = U B$ and $A ^ { T } U = V { \bar { B } } ^ { T }$ . Given the column partitionings

$$
U = \left[ u _ {1} \mid \dots \mid u _ {m} \right], \quad V = \left[ v _ {1} \mid \dots \mid v _ {n} \right],
$$

we have

$$
A v _ {k} = \alpha_ {k} u _ {k} + \beta_ {k - 1} u _ {k - 1}, \tag {10.4.2}
$$

$$
A ^ {T} u _ {k} = \alpha_ {k} v _ {k} + \beta_ {k} v _ {k + 1} \tag {10.4.3}
$$

for $k = 1 { : } n$ with the convention that $\beta _ { 0 } u _ { 0 } \equiv 0$ and $\beta _ { n } v _ { n + 1 } \equiv 0$ . Define the vectors

$$
r _ {k} = A v _ {k} - \beta_ {k - 1} u _ {k - 1}, \tag {10.4.4}
$$

$$
p _ {k} = A ^ {T} u _ {k} - \alpha_ {k} v _ {k}. \tag {10.4.5}
$$

Using (10.4.2), (10.4.4), and the orthonormality of the u-vectors, we have

$$
\alpha_ {k} = \pm \| r _ {k} \| _ {2},
$$

$$
u _ {k} = r _ {k} / \alpha_ {k}, \quad (\alpha_ {k} \neq 0).
$$

Note that if $\alpha _ { k } \ = \ 0$ , then from (10.4.1) it follows that $A ( : , 1 : k )$ is rank deficient. Similarly we may conclude from (10.4.3) and (10.4.5) that

$$
\beta_ {k} = \pm \| p _ {k} \| _ {2},
$$

$$
v _ {k + 1} = p _ {k} / \beta_ {k}, \quad (\beta_ {k} \neq 0).
$$

If $\beta _ { k } = 0$ , then it follows from the equations $A V = U B$ and $A ^ { T } U = V B ^ { T }$ that

$$
A U (:, 1: k) = V (:, 1: k) B (1: k, 1: k), \tag {10.4.6}
$$

$$
A ^ {T} V (:, 1: k) = U (:, 1: k) B (1: k, 1: k) ^ {T}, \tag {10.4.7}
$$

and thus

$$
A ^ {T} A V (:, 1: k) = V (:, 1: k) B (1: k, 1: k) ^ {T} B (1: k, 1: k).
$$

It follows that $\sigma ( B ( 1 { : } k , 1 { : } k ) ) \subseteq \sigma ( A )$ .

Properly sequenced, the above equations mathematically define the Golub-Kahan process for bidiagonalizing a rectangular matrix.

Algorithm 10.4.1 (Golub-Kahan Bidiagonalization) Given a matrix $A \in \mathbb { R } ^ { m \times n }$ with full column rank and a unit 2-norm vector $v _ { c } \in \mathbb { R } ^ { n }$ , the following algorithm computes the factorizations (10.4.6) and (10.4.7) for some k with $1 \leq k \leq n$ . The first column of V is $v _ { c }$ .

$$
k = 0, p _ {0} = v _ {c}, \beta_ {0} = 1, u _ {0} = 0
$$

while $\beta _ { k } \neq 0$

$$
v _ {k + 1} = p _ {k} / \beta_ {k}
$$

$$
k = k + 1
$$

$$
r _ {k} = A v _ {k} - \beta_ {k - 1} u _ {k - 1}
$$

$$
\alpha_ {k} = \left\| r _ {k} \right\| _ {2}
$$

$$
u _ {k} = r _ {k} / \alpha_ {k}
$$

$$
p _ {k} = A ^ {T} u _ {k} - \alpha_ {k} v _ {k}
$$

$$
\beta_ {k} = \left\| p _ {k} \right\| _ {2}
$$

end

This computation was first described by Golub and Kahan (1965). If $V _ { k } = \left[ v _ { 1 } | \cdot \cdot \cdot | v _ { k } \right]$ , $U _ { k } = [  u _ { 1 } | \cdot \cdot \cdot | u _ { k } ]$ , and

$$
B _ {k} = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & \dots & \dots & 0 \\ 0 & \alpha_ {2} & \beta_ {2} & \dots & \vdots \\ \vdots & \ddots & \ddots & \ddots & 0 \\ \vdots & & 0 & \alpha_ {k - 1} & \beta_ {k - 1} \\ 0 & \dots & 0 & 0 & \alpha_ {k} \end{array} \right], \tag {10.4.8}
$$

then after the kth pass through the loop we have

$$
A V _ {k} = U _ {k} B _ {k}, \tag {10.4.9}
$$

$$
A ^ {T} U _ {k} = V _ {k} B _ {k} ^ {T} + p _ {k} e _ {k} ^ {T}, \tag {10.4.10}
$$

assuming that $\alpha _ { k } > 0$ . It can be shown that

$$
\operatorname{span} \{v _ {1}, \dots , v _ {k} \} = \mathcal {K} (A ^ {T} A, v _ {c}, k), \tag {10.4.11}
$$

$$
\operatorname{span} \{u _ {1}, \dots , u _ {k} \} = \mathcal {K} (A A ^ {T}, A v _ {c}, k). \tag {10.4.12}
$$

Thus, the symmetric Lanczos convergence theory presented in §10.1.5 can be applied. Good approximations to A’s large singular values emerge early, while the small singular values are typically more problematic, especially if there is a cluster near the origin. For further insight, see Luk (1978), Golub, Luk, and Overton (1981), and Bj¨orck (NMLS, §7.6).

# 10.4.2 Ritz Approximations

The Ritz idea can be applied to extract approximate singular values and vectors from the matrices $U _ { k } , V _ { k }$ , and $B _ { k }$ . We simply compute the SVD

$$
F _ {k} ^ {T} B _ {k} G _ {k} = \Gamma = \operatorname{diag} \left(\gamma_ {1}, \dots , \gamma_ {k}\right) \tag {10.4.13}
$$

and form the matrices

$$
Y _ {k} = V _ {k} G _ {k} = \left[ y _ {1} \mid \dots \mid y _ {k} \right],
$$

$$
Z _ {k} = U _ {k} F _ {k} = \left[ z _ {1} \mid \dots \mid z _ {k} \right].
$$

It follows from (10.4.9), (10.4.10), and (10.4.13) that

$$
A Y _ {k} = Z _ {k} \Gamma ,
$$

$$
A ^ {T} Z _ {k} = Y _ {k} \Gamma + p _ {k} e _ {k} ^ {T} F _ {k},
$$

and so for i = 1:k we have

$$
A y _ {i} = \gamma_ {i} z _ {i}, \tag {10.4.14}
$$

$$
A ^ {T} z _ {i} = \gamma_ {i} y _ {i} + [ F _ {k} ] _ {k i} \cdot p _ {k}. \tag {10.4.15}
$$

It follows that $A ^ { T } A y _ { i } = \gamma _ { i } ^ { 2 } z _ { i } + [ F _ { k } ] _ { k i } \cdot p _ { k }$ and thus, $\{ \gamma _ { i } , y _ { i } \}$ is a Ritz pair for $A ^ { T } A$ with respect to ran(Vk).

# 10.4.3 The Tridiagonal-Bidiagonal Connection

In §8.6.1 we showed that there is a connection between the SVD of a matrix $A \in \mathbb { R } ^ { m \times n }$ and the Schur decomposition of the symmetric matrix

$$
C = \left[ \begin{array}{c c} 0 & A \\ A ^ {T} & 0 \end{array} \right]. \tag {10.4.16}
$$

In particular, if $\sigma$ is a singular value of A, then both $\sigma$ and $- \sigma$ are eigenvalues of C and the corresponding singular vectors “makeup” the corresponding eigenvectors.

Likewise, a given bidiagonalization of A can be related to a tridiagonalization of C. Assume that $m \geq n$ and that

$$
[ U _ {1} \mid U _ {2} ] ^ {T} A V = \left[ \begin{array}{c} \tilde {B} \\ 0 \end{array} \right], \qquad \tilde {B} \in \mathbb {R} ^ {n \times n},
$$

is a bidiagonalization of A with $U _ { 1 } \in \mathbb { R } ^ { m \times n } , \ U _ { 2 } \in \mathbb { R } ^ { m \times ( m - n ) }$ , and $V \in \mathbb { R } ^ { n \times n }$ . Note that

$$
Q = \left[ \begin{array}{l l} U & 0 \\ 0 & V \end{array} \right]
$$

is orthogonal and

$$
\tilde {T} = Q ^ {T} C Q = \left[ \begin{array}{c c} 0 & \tilde {B} \\ \tilde {B} ^ {T} & 0 \end{array} \right].
$$

This matrix can be symmetrically permuted into tridiagonal form. For example, in the 4-by-3 case, if $P = I _ { 7 } ( : , [ 5 1 6 2 7 3 4 ] )$ , then the reordering $\tilde { T }  P \tilde { T } P ^ { T }$ has the form

$$
\left[\begin{array}{c c c c c c c}0&0&0&0&\alpha_ {1}&\beta_ {1}&0\\0&0&0&0&0&\alpha_ {2}&\beta_ {2}\\0&0&0&0&0&0&\alpha_ {3}\\0&0&0&0&0&0&0\\\hline \alpha_ {1}&0&0&0&0&0&0\\\beta_ {1}&\alpha_ {2}&0&0&0&0&0\\0&\beta_ {2}&\alpha_ {3}&0&0&0&0\end{array}\right] \quad \rightarrow \quad \left[\begin{array}{c c c c c c c}0&\alpha_ {1}&0&0&0&0&0\\\alpha_ {1}&0&\beta_ {1}&0&0&0&0\\0&\beta_ {1}&0&\alpha_ {2}&0&0&0\\0&0&\alpha_ {2}&0&\beta_ {2}&0&0\\0&0&0&\beta_ {2}&0&\alpha_ {3}&0\\0&0&0&0&\alpha_ {3}&0&0\\0&0&0&0&0&0&0\end{array}\right].
$$

This points to an interesting connection between Golub-Kahan bidiagonalization (Algorithm 10.4.1) and Lanczos tridiagonalization (Algorithm 10.1.1). Suppose we apply Algorithm 10.4.1 to $A \in \mathbb { R } ^ { m \times n }$ with starting vector $v _ { c } .$ . Assume that the procedure runs for k steps and produces the bidiagonal matrix $B _ { k }$ displayed in (10.4.8). If we apply Algorithm 10.1.1 to the matrix C defined by (10.4.16) with a starting vector

$$
q _ {1} = \left[ \begin{array}{c} 0 \\ v _ {c} \end{array} \right] \in \mathbb {R} ^ {m + n} \tag {10.4.17}
$$

then after 2k steps the resulting tridiagonal matrix $T _ { 2 k }$ has a zero diagonal and a subdiagonal specified by $\big [ \alpha _ { 1 } , \beta _ { 1 } , \alpha _ { 2 } , \beta _ { 2 } , \cdot \cdot \cdot \alpha _ { k - 1 } , \beta _ { k - 1 } , \alpha _ { k } \big ]$ .

# 10.4.4 Paige-Saunders Lower Bidiagonalization

In §11.4.2 we show how the Golub-Kahan bidiagonalization can be used to solve sparse linear systems and least squares problems. It turns out that in this context, lower bidiagonalization is more useful:

$$
U ^ {T} A V = B = \left[ \begin{array}{c c c c c} \alpha_ {1} & 0 & \dots & \dots & 0 \\ \beta_ {2} & \alpha_ {2} & 0 & \dots & \vdots \\ \vdots & \beta_ {3} & \ddots & \ddots & \vdots \\ \vdots & & \ddots & \alpha_ {n - 1} & 0 \\ 0 & \dots & \dots & \beta_ {n - 1} & \alpha_ {n} \\ 0 & \dots & \dots & 0 & \beta_ {n} \\ \hline & & 0 \end{array} \right]. \tag {10.4.18}
$$

Proceeding as in the derivation of the Golub-Kahan bidiagonalization, we compare columns in the equations $A ^ { T } U = V B ^ { T }$ and $A V = U B$ . If $U = \ [ \ u _ { 1 } \ | \cdot \cdot \cdot | \ u _ { m } \ ]$ and $V = [  v _ { 1 } | \cdot \cdot \cdot | v _ { n } ]$ are column partitionings and we define $\beta _ { 1 } v _ { 0 } \equiv 0$ and $\alpha _ { n + 1 } v _ { n + 1 } \equiv 0$ , then for $k = 1 { : } n$ we have $A ^ { T } u _ { k } = \beta _ { k } v _ { k - 1 } + \alpha _ { k } v _ { k }$ and $A v _ { k } = \alpha _ { k } u _ { k } + \beta _ { k + 1 } u _ { k + 1 }$ . Leaving the rest of the derivation to the exercises, we obtain the following.

Algorithm 10.4.2 (Paige-Saunders Bidiagonalization) Given a matrix $A \in \mathbb { R } ^ { m \times n }$ with the property that $A ( 1 { : } n , 1 { : } n )$ is nonsingular and a unit 2-norm vector $u _ { c } \in \mathbb { R } ^ { n }$ , the following algorithm computes the factorization $A V ( : , 1 : k ) = U ( : , 1 : k + 1 ) B ( 1 : k + 1 , 1 : k )$ where $U , V$ , and B are given by (10.4.18). The first column of U is $u _ { c }$ and the integer k satisfies $1 \leq k \leq n$ .

$$
k = 1, p _ {0} = u _ {c}, \beta_ {1} = 1, v _ {0} = 0
$$

while $\beta _ { k } > 0$

$$
u _ {k} = p _ {k - 1} / \beta_ {k}
$$

$$
r _ {k} = A ^ {T} u _ {k} - \beta_ {k} v _ {k - 1}
$$

$$
\alpha_ {k} = \left\| r _ {k} \right\| _ {2}
$$

$$
v _ {k} = r _ {k} / \alpha_ {k}
$$

$$
p _ {k} = A v _ {k} - \alpha_ {k} u _ {k}
$$

$$
\beta_ {k + 1} = \left\| p _ {k} \right\| _ {2}
$$

$$
k = k + 1
$$

end

It can be shown that after k passes through the loop we have

$$
A V (:, 1: k) = U (:, 1: k) B (1: k, 1: k) + p _ {k} e _ {k} ^ {T} \tag {10.4.19}
$$

where $e _ { k } = I _ { k } ( : , k )$ . See Paige and Saunders (1982) for more details. Their bidiagonalization is equivalent to Golub-Kahan bidiagonalization applied to $[ b | A ]$ .


---

<!-- golub_600_649 -->

The need to extract information from unimaginably large datasets has prompted the development of matrix methods that involve randomization. The idea is to develop matrix approximations that are very fast to compute because they rely on limited, random samplings of the given matrix. To give a snapshot of this increasingly important paradigm for large-scale matrix computations, we consider the problem of computing a rank-k approximation to a given matrix $A \in \mathbb { R } ^ { m \times n }$ . For clarity we assume that $k \leq \mathsf { r a n k } ( A )$ . Recall that if $A = \tilde { Z } \tilde { \Sigma } \tilde { Y } ^ { T }$ is the SVD of A, then

$$
\tilde {A} _ {k} = \tilde {Z} _ {1} \tilde {\Sigma} _ {1} \tilde {Y} _ {1} ^ {T} = \tilde {Z} _ {1} \tilde {Z} _ {1} ^ {T} A \tag {10.4.20}
$$

where $\tilde { Z } _ { 1 } = \tilde { Z } ( : , 1 : k ) , \tilde { \Sigma } _ { 1 } = \tilde { \Sigma } ( 1 : k , 1 : k )$ , and $\tilde { Y } _ { 1 } = \tilde { Y } ( : , 1 { : } k )$ , is the closest rank-k matrix to A as measured in either the 2-norm or Frobenius norm. We assume that A is so large that the Krylov methods just discussed are impractical.

Drineas, Kannan, and Mahoney (2006c) propose a method that approximates the intractable $\tilde { A } _ { k }$ with a rank-k matrix of the form

$$
A _ {k} = C U R, \quad C \in \mathbb {R} ^ {m \times c}, U \in \mathbb {R} ^ {c \times r}, R \in \mathbb {R} ^ {r \times n}, k \leq c, k \leq r \tag {10.4.21}
$$

where the matrices C and R are comprised of randomly chosen values taken from A. The integers c and r are parameters of the method. Discussion of the CUR decomposition (10.4.21) nicely illustrates the notion of random sampling in the matrix context and the idea of a probabilistic error bound.

The first step in the CUR framework is to determine C. Each column of this matrix is a scaled, randomly-selected column of A:

Determine column probabilities $q _ { j } = \parallel A ( : , j ) \parallel _ { 2 } / \parallel A \parallel _ { _ { F } } ^ { 2 } , j = 1 : n .$

for $t = 1 { : } c$

Randomly pick col $( t ) \in \{ 1 , 2 , \ldots , n \}$ with $q _ { \alpha }$ the probability that $c o l ( t ) = \alpha$

$$
C (:, t) = A (:, c o l (t)) / \sqrt {c q _ {c o l (t)}}
$$

end

It follows that $C = A ( : , c o l ) D _ { c }$ where $D _ { c } \in \mathbb { R } ^ { c \times c }$ is a diagonal scaling matrix.

The matrix R is similarly constructed. Each row of this matrix is a scaled, randomly-selected row of A:

Determine row probabilities $p _ { i } = \parallel A ( i , : ) \parallel _ { 2 } / \parallel A \parallel _ { F } ^ { 2 } , i = 1 { : } m .$

for $t = 1 { : } r$

Randomly pick row $( t ) \in \{ 1 , 2 , \ldots , m \}$ with $p _ { \alpha }$ the probability that $r o w ( t ) = \alpha$

$$
R (t,:) = A (r o w (t),:) / \sqrt {r p _ {r o w (t)}}
$$

end

The matrix R has the form ${ \cal R } ~ = ~ D _ { \scriptscriptstyle R } A ( r o w , \colon )$ where $D _ { R } \in \mathbb { R } ^ { r \times r }$ is a diagonal scaling matrix.

The next step is to choose a rank-k matrix U so that $A _ { k } = C U R$ is close to the best rank-k approximation $\tilde { A } _ { k }$ . In the CUR framework, this requires the SVD

$$
C = Z \Sigma Y ^ {T} = Z _ {1} \Sigma_ {1} Y _ {1} ^ {T} + Z _ {2} \Sigma_ {2} Y _ {2}
$$

where $Z _ { 1 } = Z ( : , 1 : k ) , \ : \Sigma _ { 1 } = \Sigma ( 1 : k , 1 : k )$ , and $Y _ { 1 } = Y ( : , 1 { : } k )$ . The matrix U is then given by

$$
U = \Phi \Psi^ {T}, \qquad \Phi = Y _ {1} \Sigma_ {1} ^ {- 2} Y _ {1} ^ {T}, \Psi = D _ {R} C (r o w,:).
$$

With these definitions, simple manipulations confirm that

$$
C \Phi = Z _ {1} \Sigma_ {1} ^ {- 1} Y _ {1} ^ {T}, \tag {10.4.22}
$$

$$
\Psi^ {T} R = \left(D _ {R} \left(Z _ {1} (r o w,:) \Sigma_ {1} Y _ {1} ^ {T} + Z _ {2} (r o w,:) \Sigma_ {2} Y _ {2} ^ {T})\right) ^ {T} D _ {R} A (r o w,:), \right. \tag {10.4.23}
$$

and

$$
C U R = (C \Phi) (\Psi R) = Z _ {1} \left(D _ {R} Z _ {1} (r o w,:)\right) \left(D _ {R} A (r o w,:)\right). \tag {10.4.24}
$$

An analysis that critically depends on the selection probabilities $\left\{ q _ { i } \right\}$ and $\{ p _ { i } \}$ shows that ran $( Z _ { 1 } ) \approx \mathsf { r a n } ( \tilde { Z } _ { 1 } )$ and $\left( D _ { R } Z _ { 1 } ( r o w , : ) \right) ^ { T } \left( D _ { R } A ( r o w , : ) \right) \approx Z _ { 1 } ^ { T } A$ . Upon comparison with (10.4.20) we see that $C U R \approx Z _ { 1 } Z _ { 1 } ^ { T } A \approx \tilde { Z } _ { 1 } \tilde { Z } _ { 1 } ^ { T } A = \tilde { A } _ { k }$ . Moreover, given $\epsilon > 0$ , $\delta > 0$ , and k, it is possible to choose the parameters r and c so that the inequality

$$
\left\| A - C U R \right\| _ {F} \leq \left\| A - \tilde {A} _ {k} \right\| _ {F} + \epsilon \left\| A \right\| _ {F}
$$

holds with probability 1 − δ. Lower bounds for r and c that depend inversely on 
 and $\delta$ are given by Drineas, Kannan, and Mahoney (2006c).

# Problems

P10.4.1 Verify Equations (10.4.6), (10.4.7), (10.4.9), and (10.4.10).

P10.4.2 Corresponding to (10.3.1), develop an implementation of Algorithm 10.4.1 that involves a minimum number of vector workspaces.

P10.4.3 Show that if rank(A) = n, then the bidiagonal matrix B in (10.4.18) cannot have a zero on its diagonal.

P10.4.4 Prove (10.4.19). What can you say about $U ( : , 1 : k )$ and $V ( : , 1 : k )$ if $\beta _ { k + 1 } = 0$ in Algorithm 10.4.2?

P10.4.5 Analogous to (10.4.11)-(10.4.12), show that for Algorithm 10.4.2 we have

$$
\operatorname{span} \left\{v _ {1}, \dots , v _ {k} \right\} = \mathcal {K} \left(A ^ {T} A, A ^ {T} u _ {c}, k\right), \quad \operatorname{span} \left\{u _ {1}, \dots , u _ {k} \right\} = \mathcal {K} \left(A A ^ {T}, u _ {c}, k\right).
$$

P10.4.6 Suppose C and $q _ { 1 }$ are defined by (10.4.16) and (10.4.17) respectively. (a) Show that

$$
\mathcal {K} (C, q _ {1}, 2 k) = \text { span } \left\{\left[ \begin{array}{c} 0 \\ v _ {c} \end{array} \right], \left[ \begin{array}{c} A v _ {c} \\ 0 \end{array} \right], \left[ \begin{array}{c} 0 \\ A ^ {T} A v _ {c} \end{array} \right], \ldots , \left[ \begin{array}{c} 0 \\ (A ^ {T} A) ^ {k - 1} v _ {c} \end{array} \right], \left[ \begin{array}{c} A (A ^ {T} A) ^ {k - 1} v _ {c} \\ 0 \end{array} \right] \right\}.
$$

(b) Rigorously prove the claim made in §10.4.3 about the subdiagonal of $T _ { 2 k }$ . (c) State and prove analogous results when the Paige-Saunders bidiagonalization is used.

P10.4.7 Verify Equations 10.4.22–10.4.24.

# Notes and References for 10.4

For a more comprehensive treatment of Golub-Kahan bidiagonalization, see Bj¨orck (NMLS, §7.6). The relevance of the Lanczos process to the bidiagonalization of a rectangular matrix was first presented in:

G.H. Golub and W. Kahan (1965). “Calculating the Singular Values and Pseudo-Inverse of a Matrix,” SIAM J. Numer. Anal. Ser. B, 2, 205–224.

The idea of using Golub-Kahan bidiagonalization to solve large sparse linear systems and least squares problems started with the paper:

C.C. Paige (1974). “Bidiagonalization of Matrices and Solution of Linear Equations,” SIAM J. Numer. Anal. 11, 197–209.   
We shall have more to say about this in the next chapter. It is in anticipation of that discussion that we presented the lower bidiagonal scheme, see:   
C.C. Paige and M.A. Saunders (1982). “LSQR, An Algorithm for Sparse Linear Equations and Sparse Least Squares,” ACM Trans. Math. Softw. 8, 43–71.   
For practical implementation issues, see:   
G.H. Golub, F.T. Luk, and M.L. Overton (1981). “A Block Lanczos Method for Computing the Singular Values and Corresponding Singular Vectors of a Matrix,” ACM Trans. Math. Softw. 7, 149–169.   
J. Cullum, R.A. Willoughby, and M. Lake (1983). “A Lanczos Algorithm for Computing Singular Values and Vectors of Large Matrices,” SIAM J. Sci. Stat. Comput. 4, 197–215.   
M. Berry (1992). “Large-Scale Sparse Singular Value Computations,” International J. Supercomputing Appl. 6, 13–49.   
M. Berry and R.L. Auerbach (1993). “A Block Lanczos SVD Method with Adaptive Reorthogonalization,” in Proceedings of the Cornelius Lanczos International Centenary Conference, Raleigh, NC, SIAM Publications, Philadelphia, PA.   
Z. Jia and D. Niu (2003). “An Implicitly Restarted Refined Bidiagonalization Lanczos Method for Computing a Partial Singular Value Decomposition,” SIAM J. Matrix Anal. Applic. 25, 246–265.   
Interesting applications of the Lanczos bidiagonalization include:   
D.P. OLeary and J.A. Simmons (1981). “A Bidiagonalization-Regularization Procedure for Large Scale Discretizations of Ill-Posed Problems,” SIAM J. Sci. Stat. Comput. 2, 474–489.   
D. Calvetti, G.H. Golub, and L. Reichel (1999). “Estimation of the L-curve via Lanczos Bidiagonalization,” BIT 39, 603–619.   
H.D. Simon and H. Zha (2000). “Low-Rank Matrix Approximation Using the Lanczos Bidiagonalization Process with Applications,” SIAM J. Sci. Comput. 21, 2257–2274.   
Our sketch of the CUR decomposition framework is based on:   
P. Drineas, R. Kannan, and M.W. Mahoney (2006). “Fast Monte Carlo Algorithms for Matrices III: Computing an Efficient Approximate Decomposition of a Matrix,” SIAM J. Comput. 36, 184–206.   
Additional references concerned with randomization in matrix computations include:   
P. Drineas, R. Kannan, and M. W. Mahoney (2006). “Fast Monte Carlo Algorithms for Matrices I: Approximating Matrix Multiplication,” SIAM J. Comput. 36, 132–157.   
P. Drineas, R. Kannan, and M.W. Mahoney (2006). “Fast Monte Carlo Algorithms for Matrices II: Computing Low-Rank Approximations to a Matrix,” SIAM J. Comput. 36, 158–183.   
M.W. Mahoney, M. Maggioni, and P. Drineas (2008). “Tensor-CUR Decompositions For Tensor-Based Data,” SIAM J. Mat. Anal. Applic. 30, 957–987.   
P. Drineas, M.W. Mahoney, and S. Muthukrishnan (2008). “Relative-Error CUR Matrix Decompositions,” SIAM J. Mat. Anal. Applic. 30, 844–881.   
E. Liberty, F. Woolfe, P.-G. Martinsson, V. Rokhlin, and M.Tygert (2008). “Randomized Algorithms for the Low-Rank Approximation of Matrices,” Proc. Natl. Acad. Sci. 104, 20167–20172.   
V. Rokhlin and Mark Tygert (2008). “A Fast Randomized Algorithm for Overdetermined Linear Least-Squares Regression,” Proc. Natl. Acad. Sci. 105, 13212–13217.   
M.W. Mahoney and P. Drineas (2009). “CUR Matrix Decompositions for Improved Data Analysis,” Proc. Natl. Acad. Sci. 106, 697–702.   
D. Achlioptas and F. McSherry (2007). “Fast Computation of Low-Rank Matrix Approximations,” JACM 54(2), Article No. 9.   
V. Rokhlin, A. Szlam, and M. Tygert (2010). “A Randomized Algorithm for Principal Component Analysis,” SIAM J. Mat. Anal. Applic. 31, 1100–1124   
M.W. Mahoney (2011). “Randomized Algorithms for Matrices and Data,” Foundations and Trends in Machine Learning 3, 123–224.   
N. Halko, P.G. Martinsson, and J.A. Tropp (2011). “Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions,” SIAM Review 53, 217–288   
For another perspective on the increasing important role of randomness in matrix computations, see:

A. Edelman and N. Raj Rao (2005). “Random Matrix Theory,” Acta Numerica 14, 233–297
