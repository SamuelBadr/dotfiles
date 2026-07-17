# 7.1 Properties and Decompositions

In this section the background necessary to develop and analyze the eigenvalue algorithms that follow are surveyed. For further details, see Horn and Johnson (MA).

# 7.1.1 Eigenvalues and Invariant Subspaces

The eigenvalues of a matrix $A \in \mathbb { C } ^ { n \times n }$ are the n roots of its characteristic polynomial $p ( z ) = \mathsf { d e t } ( z I - A )$ . The set of these roots is called the spectrum of A and is denoted by

$$
\lambda (A) = \{z: \det (z I - A) = 0 \}.
$$

If $\lambda ( A ) = \left\{ \lambda _ { 1 } , \ldots , \lambda _ { n } \right\}$ , then

$$
\det (A) = \lambda_ {1} \lambda_ {2} \dots \lambda_ {n}
$$

and

$$
\operatorname{tr} (A) = \lambda_ {1} + \dots + \lambda_ {n}
$$

where the trace function, introduced in §6.4.1, is the sum of the diagonal entries, i.e.,

$$
\operatorname{tr} (A) = \sum_ {i = 1} ^ {n} a _ {i i}.
$$

These characterizations of the determinant and the trace follow by looking at the constant term and the coefficient of $z ^ { n - 1 }$ in the characteristic polynomial.

Four other attributes associated with the spectrum of $A \in \mathbb { C } ^ { n \times n }$ include the

$$
\text { Spectral   Radius }: \quad \rho (A) = \max _ {\lambda \in \lambda (A)} | \lambda |, \tag {7.1.1}
$$

$$
\text { Spectral   Abscissa }: \quad \alpha (A) = \max _ {\lambda \in \lambda (A)} \operatorname{Re} (\lambda), \tag {7.1.2}
$$

$$
\text { Numerical   Radius }: \quad r (A) = \max _ {\lambda \in \lambda (A)} \left\{\left| x ^ {H} A x \right|: \| x \| _ {2} = 1 \right\}, \tag {7.1.3}
$$

$$
\text { Numerical   Range: } \quad W (A) = \{x ^ {H} A x: \| x \| _ {2} = 1 \}. \tag {7.1.4}
$$

The numerical range, which is sometimes referred to as the field of values, obviously includes λ(A). It can be shown that $W ( A )$ is convex.

If $\lambda \in \lambda ( A )$ , then the nonzero vectors $x \in \mathbb { C } ^ { n }$ that satisfy $A x = \lambda x$ are eigenvectors. More precisely, x is a right eigenvector for λ if $A x = \lambda x$ and a left eigenvector if $x ^ { H } A = \lambda x ^ { H }$ . Unless otherwise stated, “eigenvector” means “right eigenvector.”

An eigenvector defines a 1-dimensional subspace that is invariant with respect to premultiplication by A. A subspace $S \subseteq \mathbb { C } ^ { n }$ with the property that

$$
x \in S \Longrightarrow A x \in S
$$

is said to be invariant (for A). Note that if

$$
A X = X B, \qquad B \in \mathbb {C} ^ {k \times k}, X \in \mathbb {C} ^ {n \times k},
$$

then $\mathsf { r a n } ( X )$ is invariant and $B y = \lambda y \Rightarrow A ( X y ) = \lambda ( X y )$ . Thus, if X has full column rank, then $A X = X B$ implies that $\lambda ( B ) \subseteq \lambda ( A )$ . If X is square and nonsingular, then A and $B = X ^ { - 1 } A X$ are similar, X is a similarity transformation, and $\lambda ( A ) = \lambda ( B )$ .

# 7.1.2 Decoupling

Many eigenvalue computations involve breaking the given problem down into a collection of smaller eigenproblems. The following result is the basis for these reductions.

Lemma 7.1.1. If $T \in \mathbb { C } ^ { n \times n }$ is partitioned as follows,

$$
T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {q} ^ {p}
$$

then $\lambda ( T ) = \lambda ( T _ { 1 1 } ) \cup \lambda ( T _ { 2 2 } )$ .

Proof. Suppose

$$
T x = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right] = \lambda \left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right]
$$

where $x _ { 1 } \in \mathbb { C } ^ { p }$ and $x _ { 2 } \in \mathbb { C } ^ { q }$ . If $x _ { 2 } \neq 0$ , then $T _ { 2 2 } x _ { 2 } = \lambda x _ { 2 }$ and so $\lambda \in \lambda ( T _ { 2 2 } )$ . If $x _ { 2 } = 0$ , then $T _ { 1 1 } x _ { 1 } = \lambda x _ { 1 }$ and so $\lambda \in \lambda ( T _ { 1 1 } )$ . It follows that $\lambda ( T ) \subset \lambda ( T _ { 1 1 } ) \cup \lambda ( T _ { 2 2 } )$ . But since both λ(T ) and $\lambda ( T _ { 1 1 } ) \cup \lambda ( T _ { 2 2 } )$ have the same cardinality, the two sets are equal.

# 7.1.3 Basic Unitary Decompositions

By using similarity transformations, it is possible to reduce a given matrix to any one of several canonical forms. The canonical forms differ in how they display the eigenvalues and in the kind of invariant subspace information that they provide. Because of their numerical stability we begin by discussing the reductions that can be achieved with unitary similarity.

Lemma 7.1.2. If $A \in \mathbb { C } ^ { n \times n } , B \in \mathbb { C } ^ { p \times p }$ , and $\ b X \in \mathbb { C } ^ { n \times p }$ satisfy

$$
A X = X B, \quad \operatorname{rank} (X) = p, \tag {7.1.5}
$$

then there exists a unitary $Q \in \mathbb { C } ^ { n \times n }$ such that

$$
Q ^ {H} A Q = T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {n - p} ^ {p} \tag {7.1.6}
$$

and $\lambda ( T _ { 1 1 } ) = \lambda ( A ) \cap \lambda ( B )$ .

Proof. Let

$$
X = Q \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right], \qquad Q \in \mathbb {C} ^ {n \times n}, R _ {1} \in \mathbb {C} ^ {p \times p}
$$

be a QR factorization of X. By substituting this into (7.1.5) and rearranging we have

$$
\left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ T _ {2 1} & T _ {2 2} \end{array} \right] \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right] = \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right] B
$$

where

$$
Q ^ {H} A Q = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ T _ {2 1} & T _ {2 2} \end{array} \right] _ {n - p} ^ {p}.
$$

By using the nonsingularity of $R _ { 1 }$ and the equations $T _ { 2 1 } R _ { 1 } = 0$ and $T _ { 1 1 } R _ { 1 } = R _ { 1 } B$ , we can conclude that $T _ { 2 1 } = 0$ and $\lambda ( T _ { 1 1 } ) = \lambda ( B )$ . The lemma follows because from Lemma 7.1.1 we have $\lambda ( A ) = \lambda ( T ) = \lambda ( T _ { 1 1 } ) \cup \lambda ( T _ { 2 2 } )$ .

Lemma 7.1.2 says that a matrix can be reduced to block triangular form using unitary similarity transformations if we know one of its invariant subspaces. By induction we can readily establish the decomposition of Schur (1909).

Theorem 7.1.3 (Schur Decomposition). If $A \in \mathbb { C } ^ { n \times n }$ , then there exists a unitary $Q \in \mathbb { C } ^ { n \times n }$ such that

$$
Q ^ {H} A Q = T = D + N \tag {7.1.7}
$$

where $D = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ and $N \in \mathbb { C } ^ { n \times n }$ is strictly upper triangular. Furthermore, $Q$ can be chosen so that the eigenvalues $\lambda _ { i }$ appear in any order along the diagonal.

Proof. The theorem obviously holds if $n = 1$ . Suppose it holds for all matrices of order $n - 1$ or less. If $A x = \lambda x$ and $x \neq 0$ , then by Lemma 7.1.2 (with $B = ( \lambda ) )$ there exists a unitary $U$ such that

$$
U ^ {H} A U = \left[ \begin{array}{c c} \lambda & w ^ {H} \\ 0 & C \end{array} \right] _ {n - 1} ^ {1}.
$$

By induction there is a unitary $\tilde { U }$ such that $\tilde { U } ^ { H } C \tilde { U }$ is upper triangular. Thus, if $Q = U { \cdot } \mathrm { d i a g } ( 1 , \tilde { U } )$ , then $Q ^ { H } A Q$ is upper triangular.

If $Q = [  q _ { 1 } | \cdots |  q _ { n }  ]$ is a column partitioning of the unitary matrix $Q$ in (7.1.7), then the $q _ { i }$ are referred to as Schur vectors. By equating columns in the equations $A Q = Q T$ , we see that the Schur vectors satisfy

$$
A q _ {k} = \lambda_ {k} q _ {k} + \sum_ {i = 1} ^ {k - 1} n _ {i k} q _ {i}, \quad k = 1: n. \tag {7.1.8}
$$

From this we conclude that the subspaces

$$
S _ {k} = \operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\}, \quad k = 1: n,
$$

are invariant. Moreover, it is not hard to show that if $Q _ { k } \ = \ \left[ \ q _ { 1 } \ | \cdot \cdot \cdot | \ q _ { k } \ \right]$ , then $\lambda ( Q _ { k } ^ { H } A Q _ { k } ) = \{ \lambda _ { 1 } , . . . , \lambda _ { k } \}$ . Since the eigenvalues in (7.1.7) can be arbitrarily ordered, it follows that there is at least one k-dimensional invariant subspace associated with each subset of k eigenvalues. Another conclusion to be drawn from (7.1.8) is that the Schur vector $q _ { k }$ is an eigenvector if and only if the kth column of N is zero. This turns out to be the case for k = 1:n whenever $A ^ { H } A = A A ^ { H }$ . Matrices that satisfy this property are called normal.

Corollary 7.1.4. $A \in \mathbb { C } ^ { n \times n }$ is normal if and only if there exists a unitary $Q \in \mathbb { C } ^ { n \times n }$ such that $Q ^ { H } A Q = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ .

Proof. See P7.1.1.

Note that if $Q ^ { H } A Q = T = \mathrm { d i a g } ( \lambda _ { i } ) + N$ is a Schur decomposition of a general n-by-n matrix A, then $\Vert N \Vert _ { F }$ is independent of the choice of $Q \mathrm { : }$

$$
\| N \| _ {F} ^ {2} = \| A \| _ {F} ^ {2} - \sum_ {i = 1} ^ {n} | \lambda_ {i} | ^ {2} \equiv \Delta^ {2} (A).
$$

This quantity is referred to as A’s departure from normality. Thus, to make T “more diagonal,” it is necessary to rely on nonunitary similarity transformations.

# 7.1.4 Nonunitary Reductions

To see what is involved in nonunitary similarity reduction, we consider the block diagonalization of a 2-by-2 block triangular matrix.

Lemma 7.1.5. Let $T \in \mathbb { C } ^ { n \times n }$ be partitioned as follows:

$$
T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {q} ^ {p}.
$$

Define the linear transformation $\phi { : } \mathbb { C } ^ { p \times q } \to \mathbb { C } ^ { p \times q }$ by

$$
\phi (X) = T _ {1 1} X - X T _ {2 2}
$$

where $X \in \mathbb { C } ^ { p \times q }$ . Then $\phi$ is nonsingular if and only if $\lambda ( T _ { 1 1 } ) \cap \lambda ( T _ { 2 2 } ) = \emptyset$ . If φ is nonsingular and Y is defined by

$$
Y = \left[ \begin{array}{c c} I _ {p} & Z \\ 0 & I _ {q} \end{array} \right]
$$

where $\phi ( Z ) = - T _ { 1 2 }$ , then $Y ^ { - 1 } T Y = \mathrm { d i a g } ( T _ { 1 1 } , T _ { 2 2 } )$

Proof. Suppose φ(X) = 0 for $X \neq 0$ and that

$$
U ^ {H} X V = \left[ \begin{array}{c c} \Sigma_ {r} & 0 \\ 0 & 0 \end{array} \right] _ {p - r} ^ {r}
$$

is the SVD of X with $\Sigma _ { r } = \mathrm { d i a g } ( \sigma _ { i } ) , r = \mathsf { r a n k } ( X )$ . Substituting this into the equation $T _ { 1 1 } X = X T _ { 2 2 }$ gives

$$
\left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right] \left[ \begin{array}{c c} \Sigma_ {r} & 0 \\ 0 & 0 \end{array} \right] = \left[ \begin{array}{c c} \Sigma_ {r} & 0 \\ 0 & 0 \end{array} \right] \left[ \begin{array}{c c} B _ {1 1} & B _ {1 2} \\ B _ {2 1} & B _ {2 2} \end{array} \right]
$$

where ${ U ^ { H } T _ { 1 1 } U = ( A _ { i j } ) }$ and $V ^ { H } T _ { 2 2 } V = \left( B _ { i j } \right)$ . By comparing blocks in this equation it is clear that $A _ { 2 1 } = 0 , B _ { 1 2 } = 0$ , and $\lambda ( A _ { 1 1 } ) = \lambda ( B _ { 1 1 } )$ . Consequently, $A _ { 1 1 }$ and $B _ { 1 1 }$ have an eigenvalue in common and that eigenvalue is in $\lambda ( T _ { 1 1 } ) \cap \lambda ( T _ { 2 2 } )$ . Thus, if φ is singular, then $T _ { 1 1 }$ and $T _ { 2 2 }$ have an eigenvalue in common. On the other hand, if $\lambda \in \lambda ( T _ { 1 1 } ) \cap \lambda ( T _ { 2 2 } )$ , then we have eigenvector equations $T _ { 1 1 } x = \lambda x$ and $y ^ { H } T _ { 2 2 } = \lambda y ^ { H }$ . A calculation shows that $\phi ( x y ^ { H } ) = 0$ confirming that φ is singular.

Finally, if $\phi$ is nonsingular, then $\phi ( Z ) = - T _ { 1 2 }$ has a solution and

$$
Y ^ {- 1} T Y = \left[ \begin{array}{c c} I _ {p} & - Z \\ 0 & I _ {q} \end{array} \right] \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] \left[ \begin{array}{c c} I _ {p} & Z \\ 0 & I _ {q} \end{array} \right] = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 1} Z - Z T _ {2 2} + T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right]
$$

has the required block diagonal form.

By repeatedly applying this lemma, we can establish the following more general result.

Theorem 7.1.6 (Block Diagonal Decomposition). Suppose

$$
Q ^ {H} A Q = T = \left[ \begin{array}{c c c c} T _ {1 1} & T _ {1 2} & \dots & T _ {1 q} \\ 0 & T _ {2 2} & \dots & T _ {2 q} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & T _ {q q} \end{array} \right] \tag {7.1.9}
$$

is a Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ and that the $T _ { i i }$ are square. $I f \lambda ( T _ { i i } ) \cap \lambda ( T _ { j j } ) = \emptyset$ whenever $i \neq j$ , then there exists a nonsingular matrix $Y \in \mathbb { C } ^ { n \times n }$ such that

$$
(Q Y) ^ {- 1} A (Q Y) = \operatorname{diag} \left(T _ {1 1}, \dots , T _ {q q}\right). \tag {7.1.10}
$$

Proof. See P7.1.2.

If each diagonal block $T _ { i i }$ is associated with a distinct eigenvalue, then we obtain

Corollary 7.1.7. If $A \in \mathbb { C } ^ { n \times n }$ , then there exists a nonsingular X such that

$$
X ^ {- 1} A X = \operatorname{diag} \left(\lambda_ {1} I + N _ {1}, \dots , \lambda_ {q} I + N _ {q}\right) \quad N _ {i} \in \mathbb {C} ^ {n _ {i} \times n _ {i}} \tag {7.1.11}
$$

where $\lambda _ { 1 } , \ldots , \lambda _ { q }$ are distinct, the integers $n _ { 1 } , \ldots , n _ { q }$ satisfy $n _ { 1 } + \cdot \cdot \cdot + n _ { q } = n$ , and each $N _ { i }$ is strictly upper triangular.

A number of important terms are connected with decomposition (7.1.11). The integer $n _ { i }$ is referred to as the algebraic multiplicity of $\lambda _ { i }$ . If $n _ { i } = 1$ , then $\lambda _ { i }$ is said to be simple. The geometric multiplicity of $\lambda _ { i }$ equals the dimensions of null $( N _ { i } )$ , i.e., the number of linearly independent eigenvectors associated with $\lambda _ { i }$ . If the algebraic multiplicity of $\lambda _ { i }$ exceeds its geometric multiplicity, then $\lambda _ { i }$ is said to be a defective eigenvalue. A matrix with a defective eigenvalue is referred to as a defective matrix. Nondefective matrices are also said to be diagonalizable.

Corollary 7.1.8 (Diagonal Form). $A \in \mathbb { C } ^ { n \times n }$ is nondefective if and only if there exists a nonsingular $\ b X \in \mathbb { C } ^ { n \times n }$ such that

$$
X ^ {- 1} A X = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right). \tag {7.1.12}
$$

Proof. A is nondefective if and only if there exist independent vectors $\boldsymbol { x } _ { 1 } \ldots \boldsymbol { x } _ { n } \in \mathbb { C } ^ { n }$ and scalars $\lambda _ { 1 } , \ldots , \lambda _ { n }$ such that $A x _ { i } = \lambda _ { i } x _ { i }$ for $i = 1 { : } n$ . This is equivalent to the existence of a nonsingular $\ b X = [ \ b x _ { 1 } \ | \cdot \ b \cdot \cdot | \ \ b x _ { n } ] \in \mathbb { C } ^ { n \times n }$ such that $A X \ = \ X D$ where $D = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ .□

Note that if $y _ { i } ^ { H }$ is the ith row of $X ^ { - 1 }$ , then $y _ { i } ^ { H } A = \lambda _ { i } y _ { i } ^ { H }$ . Thus, the columns of $X ^ { - H }$ are left eigenvectors and the columns of X are right eigenvectors.

If we partition the matrix X in (7.1.11),

$$
X = \left[ \begin{array}{c} X _ {1} \mid \dots \mid X _ {q} \\ n _ {1} \end{array} \right]
$$

then $\mathbb { C } ^ { n } = \mathsf { r a n } ( X _ { 1 } ) \oplus \ldots \oplus \mathsf { r a n } ( X _ { q } )$ , a direct sum of invariant subspaces. If the bases for these subspaces are chosen in a special way, then it is possible to introduce even more zeroes into the upper triangular portion of $X ^ { - 1 } A X$ .

Theorem 7.1.9 (Jordan Decomposition). If $A \in \mathbb { C } ^ { n \times n }$ , then there exists a nonsingular $\ b X \in \mathbb { C } ^ { n \times n }$ such that $X ^ { - 1 } A X = \operatorname { d i a g } ( J _ { 1 } , \dots , J _ { q } )$ where

$$
J _ {i} = \left[ \begin{array}{c c c c c} \lambda_ {i} & 1 & & \dots & 0 \\ 0 & \lambda_ {i} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & 1 \\ 0 & \dots & & 0 & \lambda_ {i} \end{array} \right] \in \mathbb {C} ^ {n _ {i} \times n _ {i}}
$$

and $n _ { 1 } + \cdots + n _ { q } = n$ .

Proof. See Horn and Johnson (MA, p. 330)

The $J _ { i }$ are referred to as Jordan blocks. The number and dimensions of the Jordan blocks associated with each distinct eigenvalue are unique, although their ordering along the diagonal is not.

# 7.1.5 Some Comments on Nonunitary Similarity

The Jordan block structure of a defective matrix is difficult to determine numerically. The set of n-by-n diagonalizable matrices is dense in Cn×n, and thus, small changes in a defective matrix can radically alter its Jordan form. We have more to say about this in §7.6.5.

A related difficulty that arises in the eigenvalue problem is that a nearly defective matrix can have a poorly conditioned matrix of eigenvectors. For example, any matrix X that diagonalizes

$$
A = \left[ \begin{array}{c c} 1 + \epsilon & 1 \\ 0 & 1 - \epsilon \end{array} \right], \qquad 0 <   \epsilon \ll 1, \tag {7.1.13}
$$

has a 2-norm condition of order $1 / \epsilon .$

These observations serve to highlight the difficulties associated with ill-conditioned similarity transformations. Since

$$
\mathsf {f l} (X ^ {- 1} A X) = X ^ {- 1} A X + E, \tag {7.1.14}
$$

where

$$
\| E \| _ {2} \approx \mathbf {u} \cdot \kappa_ {2} (X) \| A \| _ {2}, \tag {7.1.15}
$$

it is clear that large errors can be introduced into an eigenvalue calculation when we depart from unitary similarity.

# 7.1.6 Singular Values and Eigenvalues

Since the singular values of A and its Schur decomposition $Q ^ { H } A Q = \mathrm { d i a g } ( \lambda _ { i } ) + N$ are the same, it follows that

$$
\sigma_ {\min} (A) \leq \min _ {1 \leq i \leq n} | \lambda_ {i} | \leq \max _ {1 \leq i \leq n} | \lambda_ {i} | \leq \sigma_ {\max} (A).
$$

From what we know about the condition of triangular matrices, it may be the case that

$$
\max _ {1 \leq i, j \leq n} \frac {| \lambda_ {i} |}{| \lambda_ {j} |} \ll \kappa_ {2} (A).
$$

See §5.4.3. This is a reminder that for nonnormal matrices, eigenvalues do not have the “predictive power” of singular values when it comes to $A x = b$ sensitivity matters. Eigenvalues of nonnormal matrices have other shortcomings, a topic that is the focus of §7.9.

# Problems

P7.1.1 (a) Show that if $T \in \mathbb { C } ^ { n \times n }$ is upper triangular and normal, then T is diagonal. (b) Show that if A is normal and $Q ^ { H } A Q = T$ is a Schur decomposition, then T is diagonal. (c) Use (a) and (b) to complete the proof of Corollary 7.1.4.

P7.1.2 Prove Theorem 7.1.6 by using induction and Lemma 7.1.5.

P7.1.3 Suppose $A \in \mathbb { C } ^ { n \times n }$ has distinct eigenvalues. Show that if $Q ^ { H } A Q = T$ is its Schur decomposition and $A B = B A$ , then QH BQ is upper triangular.

P7.1.4 Show that if A and $B ^ { H }$ are in $\mathbb { C } ^ { m \times n }$ with $m \geq n$ , then

$$
\lambda (A B) = \lambda (B A) \cup \{\underbrace {0 , \ldots , 0} _ {m - n} \}.
$$

P7.1.5 Given $A \in \mathbb { C } ^ { n \times n }$ , use the Schur decomposition to show that for every $\epsilon > 0$ , there exists a diagonalizable matrix B such that $\| A - B \| _ { 2 } \leq \epsilon$ . This shows that the set of diagonalizable matrices is dense in Cn×n $\mathbb { C } ^ { n \times n }$ and that the Jordan decomposition is not a continuous matrix decomposition.

P7.1.6 Suppose $A _ { k } \to A$ and that $Q _ { k } ^ { H } A _ { k } Q _ { k } = T _ { k }$ is a Schur decomposition of $A _ { k }$ . Show that $\left\{ Q _ { k } \right\}$ has a converging subsequence $\{ Q _ { k _ { i } } \}$ with the property that

$$
\lim _ {i \to \infty} Q _ {k _ {i}} = Q
$$

where $Q ^ { H } A Q = T$ is upper triangular. This shows that the eigenvalues of a matrix are continuous functions of its entries.

P7.1.7 Justify (7.1.14) and (7.1.15).

P7.1.8 Show how to compute the eigenvalues of

$$
M = \left[ \begin{array}{c c} A & C \\ B & D \\ k & j \end{array} \right] _ {j} ^ {k}
$$

where A, B, C, and D are given real diagonal matrices.

P7.1.9 Use the Jordan decomposition to show that if all the eigenvalues of a matrix A are strictly less than unity, then $\scriptstyle \operatorname* { l i m } _ { k \to \infty } A ^ { k } = 0$ .

P7.1.10 The initial value problem

$$
\begin{array}{r c l} \dot {x} (t) & = & y (t), \qquad x (0) = 1, \\ \dot {y} (t) & = & - x (t), \qquad y (0) = 0, \end{array}
$$

has solution $x ( t ) = \cos ( t )$ and $y ( t ) = \sin ( t )$ . Let $h > 0$ . Here are three reasonable iterations that can be used to compute approximations $x _ { k } \approx x ( k h )$ and $y _ { k } \approx y ( k h )$ assuming that $x _ { 0 } = 1$ and $y _ { k } = 0 ;$ :

$$
\text { Method   1: } \quad \begin{array}{r c l} x _ {k + 1} & = & x _ {k} + h y _ {k}, \\ y _ {k + 1} & = & y _ {k} - h x _ {k}, \end{array}
$$

$$
\text { Method   2: } \quad \begin{array}{r c l} x _ {k + 1} & = & x _ {k} + h y _ {k}, \\ y _ {k + 1} & = & y _ {k} - h x _ {k + 1}, \end{array}
$$

$$
\text { Method   3: } \quad \begin{array}{r c l} x _ {k + 1} & = & x _ {k} + h y _ {k + 1}, \\ y _ {k + 1} & = & y _ {k} - h x _ {k + 1}. \end{array}
$$

Express each method in the form

$$
\left[ \begin{array}{l} x _ {k + 1} \\ y _ {k + 1} \end{array} \right] = A _ {h} \left[ \begin{array}{l} x _ {k} \\ y _ {k} \end{array} \right]
$$

where $A _ { h }$ is a 2-by-2 matrix. For each case, compute $\lambda ( A _ { h } )$ and use the previous problem to discuss lim $x _ { k }$ and lim yk as $k \to \infty$ .

P7.1.11 If $J \in \mathbb { R } ^ { d \times d }$ is a Jordan block, what is $\kappa _ { \infty } ( J ) ?$

P7.1.12 Suppose A, $B \in \mathbb { C } ^ { n \times n }$ . Show that the 2n-by-2n matrices

$$
M _ {1} = \left[ \begin{array}{c c} A B & 0 \\ B & 0 \end{array} \right] \quad \text {and} \quad M _ {2} = \left[ \begin{array}{c c} 0 & 0 \\ B & B A \end{array} \right]
$$

are similar thereby showing that $\lambda ( A B ) = \lambda ( B A )$

P7.1.13 Suppose $A \in \mathbb { R } ^ { n \times n }$ . We say that $B \in \mathbb { R } ^ { n \times n }$ is the Drazin inverse of A if (i) $A B = B A$ , (ii) $B A B = B$ , and (iii) the spectral radius of $A - A B A$ is zero. Give a formula for B in terms of the Jordan decomposition of A paying particular attention to the blocks associated with A’s zero eigenvalues.

P7.1.14 Show that if $A \in \mathbb { R } ^ { n \times n }$ , then $\rho ( A ) \geq ( \sigma _ { 1 } \cdot \cdot \cdot \sigma _ { n } ) ^ { 1 / n }$ where $\sigma _ { 1 } , \ldots , \sigma _ { n }$ are the singular values of A.

P7.1.15 Consider the polynomial $q ( x ) = \operatorname* { d e t } ( I _ { n } + x A )$ where $A \in \mathbb { R } ^ { n \times n }$ . We wish to compute the coefficient of $x ^ { 2 }$ . (a) Specify the coefficient in terms of the eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { n }$ of A. (b) Give a simple formula for the coefficient in terms of $\operatorname { t r } ( A )$ and $\operatorname { t r } ( A ^ { 2 } )$ .

P7.1.16 Given $A \in \mathbb { R } ^ { 2 \times 2 } ,$ , show that there exists a nonsingular $\begin{array} { r } { X \in \mathbb { R } ^ { 2 \times 2 } \mathrm { ~ s o ~ } X ^ { - 1 } A X = A ^ { T } } \end{array}$ . See Dubrulle and Parlett (2007).

# Notes and References for 7.1

For additional discussion about the linear algebra behind the eigenvalue problem, see Horn and Johnson (MA) and:

L. Mirsky (1963). An Introduction to Linear Algebra, Oxford University Press, Oxford, U.K.

M. Marcus and H. Minc (1964). A Survey of Matrix Theory and Matrix Inequalities, Allyn and Bacon, Boston.

R. Bellman (1970). Introduction to Matrix Analysis, second edition, McGraw-Hill, New York.

I. Gohberg, P. Lancaster, and L. Rodman (2006). Invariant Subspaces of Matrices with Applications, SIAM Publications, Philadelphia, PA.

For a general discussion about the similarity connection between a matrix and its transpose, see:

A.A. Dubrulle and B.N. Parlett (2010). “Revelations of a Transposition Matrix,” J. Comp. and Appl. Math. 233, 1217–1219.

The Schur decomposition originally appeared in:

I. Schur (1909). “On the Characteristic Roots of a Linear Substitution with an Application to the Theory of Integral Equations.” Math. Ann. 66, 488-510 (German).

A proof very similar to ours is given in:

H.W. Turnbull and A.C. Aitken (1961). An Introduction to the Theory of Canonical Forms, Dover, New York, 105.
