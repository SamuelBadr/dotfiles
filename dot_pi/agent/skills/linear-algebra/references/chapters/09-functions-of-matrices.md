# Chapter 9

# Functions of Matrices

9.1 Eigenvalue Methods   
9.2 Approximation Methods   
9.3 The Matrix Exponential   
9.4 The Sign, Square Root, and Log of a Matrix

Computing a function $f ( A )$ of an n-by-n matrix A is a common problem in many application areas. Roughly speaking, if the scalar function $f ( z )$ is defined on $\lambda ( A )$ , then $f ( A )$ is defined by substituting $^ { 6 6 } A ^ { \dag }$ for $^ { 6 6 } z ^ { 5 9 }$ in the “formula” for $f ( z )$ . For example, if $f ( z ) = ( 1 + z ) / ( 1 - z )$ and $1 \not \in \lambda ( A )$ , then $f ( A ) = ( I + A ) ( I - A ) ^ { - 1 }$ .

The computations get particularly interesting when the function f is transcendental. One approach in this more complicated situation is to compute an eigenvalue decomposition $A = Y B Y ^ { - 1 }$ and use the formula $f ( A ) = Y f ( B ) Y ^ { - 1 }$ . If B is sufficiently simple, then it is often possible to calculate f (B) directly. This is illustrated in §9.1 for the Jordan and Schur decompositions.

Another class of methods involves the approximation of the desired function $f ( A )$ with an easy-to-calculate function $g ( A )$ . For example, g might be a truncated Taylor series approximation to f. Error bounds associated with the approximation of matrix functions are given in §9.2.

In §9.3 we discuss the special and very important problem of computing the matrix exponential $e ^ { A }$ . The matrix sign, square root, and logarithm functions and connections to the polar decomposition are treated in §9.4.

# Reading Notes

Knowledge of Chapters 3 and 7 is assumed. Within this chapter there are the following dependencies:

$$
\begin{array}{c c c c c} \S 9. 1 & \to & \S 9. 2 & \to & \S 9. 3 \\ & & \downarrow & & \\ & & \S 9. 4 & & \end{array}
$$

Complementary references include Horn and Johnson (TMA) and the definitive text by Higham (FOM). We mention that aspects of the f(A)-times-a-vector problem are treated in §10.2.

# 9.1 Eigenvalue Methods

Here are some examples of matrix functions:

$$
\begin{array}{l} p (A) = I + A, \\ r (A) = \left(I - \frac {A}{2}\right) ^ {- 1} \left(I + \frac {A}{2}\right), \quad 2 \not \in \lambda (A), \\ e ^ {A} = \sum_ {k = 0} ^ {\infty} \frac {A ^ {k}}{k !}. \\ \end{array}
$$

Obviously, these are the matrix versions of the scalar-valued functions

$$
\begin{array}{l} p (z) = 1 + z, \\ r (z) = (1 - (z / 2)) ^ {- 1} (1 + (z / 2)), \quad 2 \neq z, \\ e ^ {z} = \sum_ {k = 0} ^ {\infty} \frac {z ^ {k}}{k !}. \\ \end{array}
$$

Given an n-by-n matrix A, it appears that all we have to do to define f (A) is to substitute A into the formula for f. However, to make subsequent algorithmic developments precise, we need to be a little more formal. It turns out that there are several equivalent ways to define a function of a matrix. See Higham (FOM, §1.2). Because of its prominence in the literature and its simplicity, we take as our “base” definition one that involves the Jordan canonical form (JCF).

# 9.1.1 A Jordan-Based Definition

Suppose $A \in \mathbb { C } ^ { n \times n }$ and let

$$
A = X \cdot \operatorname{diag} \left(J _ {1}, \dots , J _ {q}\right) \cdot X ^ {- 1} \tag {9.1.1}
$$

be its JCF with

$$
J _ {i} = \left[ \begin{array}{c c c c c} \lambda_ {i} & 1 & \dots & \dots & 0 \\ 0 & \lambda_ {i} & 1 & \dots & \vdots \\ \vdots & \ddots & \ddots & \ddots & \vdots \\ \vdots & \vdots & \ddots & \ddots & 1 \\ 0 & \dots & \dots & 0 & \lambda_ {i} \end{array} \right] \in \mathbb {C} ^ {n _ {i} \times n _ {i}}, \quad i = 1: q. \tag {9.1.2}
$$

The matrix function f(A) is defined by

$$
f (A) = X \cdot \operatorname{diag} \left(F _ {1}, \dots , F _ {q}\right) \cdot X ^ {- 1} \tag {9.1.3}
$$

where

$$
F _ {i} = \left[ \begin{array}{c c c c c} f \left(\lambda_ {i}\right) & f ^ {(1)} \left(\lambda_ {i}\right) & \dots & \dots & \frac {f ^ {\left(n _ {i} - 1\right)} \left(\lambda_ {i}\right)}{\left(n _ {i} - 1\right) !} \\ 0 & f \left(\lambda_ {i}\right) & \ddots & \dots & \vdots \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ \vdots & \vdots & \vdots & \ddots & f ^ {(1)} \left(\lambda_ {i}\right) \\ 0 & \dots & \dots & \dots & f \left(\lambda_ {i}\right) \end{array} \right], \quad i = 1: q, \tag {9.1.4}
$$

assuming that all the required derivative evaluations exist.

# 9.1.2 The Taylor Series Representation

If f can be represented by a Taylor series on A’s spectrum, then $f ( A )$ can be represented by the same Taylor series in A. To fix ideas, assume that f is analytic in a neighborhood of $z _ { 0 } \in \mathbb { C }$ and that for some $r > 0$ we have

$$
f (z) = \sum_ {k = 0} ^ {\infty} \frac {f ^ {(k)} (z _ {0})}{k !} (z - z _ {0}) ^ {k}, \quad | z - z _ {0} | <   r. \tag {9.1.5}
$$

Our first result applies to a single Jordan block.

Lemma 9.1.1. Suppose $B \in \mathbb { C } ^ { m \times m }$ is a Jordan block and write $B = \lambda I _ { m } + E$ where E is its strictly upper bidiagonal part. Given (9.1.5), $i f \left| \lambda - z _ { 0 } \right| < r$ , then

$$
f (B) = \sum_ {k = 0} ^ {\infty} \frac {f ^ {(k)} (z _ {0})}{k !} (B - z _ {0} I _ {m}) ^ {k}.
$$

Proof. Note that powers of E are highly structured, e.g.,

$$
E = \left[ \begin{array}{l l l l} 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 \end{array} \right], \quad E ^ {2} = \left[ \begin{array}{l l l l} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{array} \right], \quad E ^ {3} = \left[ \begin{array}{l l l l} 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{array} \right].
$$

In terms of the Kronecker delta, if $0 \leq p \leq m - 1$ , then $[ { \cal E } ^ { p } ] _ { i j } = ( \delta _ { i , j - p } )$ . It follows from (9.1.4) that

$$
f (B) = \sum_ {p = 0} ^ {m - 1} f ^ {(p)} (\lambda) \frac {E ^ {p}}{p !}. \tag {9.1.6}
$$

On the other hand, if $p > m$ , then $E ^ { p } = 0$ . Thus, for any $k \geq 0$ we have

$$
(B - z _ {0} I) ^ {k} = ((\lambda - z _ {0}) I + E) ^ {k} = \sum_ {p = 0} ^ {k} \frac {k (k - 1) \cdots (k - p + 1)}{p !} \cdot (\lambda - z _ {0}) ^ {k - p} \cdot E ^ {p}
$$

$$
= \sum_ {p = 0} ^ {\min \{k, m - 1 \}} \left[ \frac {d ^ {p}}{d \lambda^ {p}} (\lambda - z _ {0}) ^ {k} \right] \frac {E ^ {p}}{p !}.
$$

If N is a nonnegative integer, then

$$
\sum_ {k = 0} ^ {N} \frac {f ^ {(k)} (z _ {0})}{k !} (B - z _ {0} I) ^ {k} = \sum_ {p = 0} ^ {\min \{k, m - 1 \}} \frac {d ^ {p}}{d \lambda^ {p}} \left(\sum_ {k = 0} ^ {N} \frac {f ^ {(k)} (z _ {0})}{k !} (\lambda - z _ {0}) ^ {k}\right) \frac {E ^ {p}}{p !}.
$$

The lemma follows by taking limits with respect to N and using both (9.1.6) and the Taylor series representation of $f ( z )$ .□

A similar result holds for general matrices.

Theorem 9.1.2. If f has the Taylor series representation (9.1.5) and $| \lambda - z _ { 0 } | < r$ for all $\lambda \in \lambda ( A )$ where $A \in \mathbb { C } ^ { n \times n }$ , then

$$
f (A) = \sum_ {k = 0} ^ {\infty} \frac {f ^ {(k)} (z _ {0})}{k !} (A - z _ {0} I) ^ {k}.
$$

Proof. Let the JCF of A be given by (9.1.1) and (9.1.2). From Lemma 9.1.1 we have

$$
f (J _ {i}) = \sum_ {k = 0} ^ {\infty} \alpha_ {k} (J _ {i} - z _ {0} I) ^ {k}, \quad \alpha_ {k} = \frac {f ^ {(k)} (z _ {0})}{k !},
$$

for $i = 1 { : } q$ . Using the definition (9.1.3) and (9.1.4) we see that

$$
f (A) = X \cdot \operatorname{diag} \left(\sum_ {k = 0} ^ {\infty} \alpha_ {k} \left(J _ {1} - z _ {0} I _ {n _ {1}}\right) ^ {k}, \dots , \sum_ {k = 0} ^ {\infty} \alpha_ {k} \left(J _ {q} - z _ {0} I _ {n _ {q}}\right) ^ {k}\right) \cdot X ^ {- 1}
$$

$$
= X \cdot \left(\sum_ {k = 0} ^ {\infty} \alpha_ {k} \left(J - z _ {0} I _ {n}\right) ^ {k}\right) \cdot X ^ {- 1}
$$

$$
= \sum_ {k = 0} ^ {\infty} \alpha_ {k} \left(X \left(J - z _ {0} I _ {n}\right) X ^ {- 1}\right) ^ {k} = \sum_ {k = 0} ^ {\infty} \alpha_ {k} \left(A - z _ {0} I _ {n}\right) ^ {k},
$$

completing the proof of the theorem.

Important matrix functions that have simple Taylor series definitions include

$$
\begin{array}{l} \exp (A) = \sum_ {k = 0} ^ {\infty} \frac {A ^ {k}}{k !}, \\ \log (I - A) = \sum_ {k = 1} ^ {\infty} \frac {A ^ {k}}{k}, \quad | \lambda | <   1, \lambda \in \lambda (A), \\ \sin (A) = \sum_ {k = 0} ^ {\infty} (- 1) ^ {k} \frac {A ^ {2 k + 1}}{(2 k + 1) !}, \\ \cos (A) = \sum_ {k = 0} ^ {\infty} (- 1) ^ {k} \frac {A ^ {2 k}}{(2 k) !}. \\ \end{array}
$$

For clarity in this section and the next, we consider only matrix functions that have a Taylor series representation. In that case it is easy to verify that

$$
A \cdot f (A) = f (A) \cdot A \tag {9.1.7}
$$

and

$$
f (X ^ {- 1} A X) = X \cdot f (A) \cdot X ^ {- 1}. \tag {9.1.8}
$$

# 9.1.3 An Eigenvector Approach

If $A \in \mathbb { C } ^ { n \times n }$ is diagonalizable, then it is particularly easy to specify f (A) in terms of A’s eigenvalues and eigenvectors.

Corollary 9.1.3. If $A \in \mathbb { C } ^ { n \times n } , A = X \cdot \mathrm { d i a g } ( \lambda _ { 1 } , \ldots , \lambda _ { n } ) \cdot X ^ { - 1 }$ , and f(A) is defined, then

$$
f (A) = X \cdot \mathrm{diag} (f (\lambda_ {1}), \dots , f (\lambda_ {n})) \cdot X ^ {- 1}. \tag {9.1.9}
$$

Proof. This result is an easy consequence of Theorem 9.1.2 since all the Jordan blocks are 1-by-1.

Unfortunately, if the matrix of eigenvectors is ill-conditioned, then computing $f ( A )$ via (9.1.8) is likely introduce errors of order u $\kappa _ { 2 } ( X )$ because of the required solution of a linear system that involves the eigenvector matrix X. For example, if

$$
A = \left[ \begin{array}{c c} 1 + 1 0 ^ {- 5} & 1 \\ 0 & 1 - 1 0 ^ {- 5} \end{array} \right],
$$

then any matrix of eigenvectors is a column-scaled version of

$$
X = \left[ \begin{array}{c c} 1 & - 1 \\ 0 & 2 (1 - 1 0 ^ {- 5}) \end{array} \right]
$$

and has a 2-norm condition number of order $1 0 ^ { 5 }$ . Using a computer with machine precision $\mathbf { u } \approx 1 0 ^ { - 7 }$ , we find

$$
\mathsf {f l} \left(X ^ {- 1} \mathrm{diag} (\exp (1 + 1 0 ^ {- 5}), \exp (1 - 1 0 ^ {- 5})) X\right) = \left[ \begin{array}{l l} 2. 7 1 8 3 0 7 & 2. 7 5 0 0 0 0 \\ 0. 0 0 0 0 0 0 & 2. 7 1 8 2 5 4 \end{array} \right]
$$

while

$$
e ^ {A} = \left[ \begin{array}{l l} 2. 7 1 8 3 0 9 & 2. 7 1 8 2 8 2 \\ 0. 0 0 0 0 0 0 & 2. 7 1 8 2 5 5 \end{array} \right].
$$

The example suggests that ill-conditioned similarity transformations should be avoided when computing a function of a matrix. On the other hand, if A is a normal matrix, then it has a perfectly conditioned matrix of eigenvectors. In this situation, computation of f(A) via diagonalization is a recommended strategy.

# 9.1.4 A Schur Decomposition Approach

Some of the difficulties associated with the Jordan approach to the matrix function problem can be circumvented by relying upon the Schur decomposition. If $A = Q T Q ^ { H }$ is the Schur decomposition of A, then by (9.1.8),

$$
f (A) = Q f (T) Q ^ {H}.
$$

For this to be effective, we need an algorithm for computing functions of upper triangular matrices. Unfortunately, an explicit expression for $f ( T )$ is very complicated.

Theorem 9.1.4. Let $T = \left( t _ { i j } \right)$ be an n-by-n upper triangular matrix with $\lambda _ { i } = t _ { i i }$ and assume $f ( T )$ is defined. I $\textit { f f } ( T ) = ( f _ { i j } )$ , then $f _ { i j } = 0 \ i f \ i > j , \ f _ { i j } \ = f ( \lambda _ { i } ) \ f o r \ i = j ,$ , and for all $i < j$ we have

$$
f _ {i j} = \sum_ {(s _ {0}, \dots , s _ {k}) \in S _ {i j}} t _ {s _ {0}, s _ {1}} t _ {s _ {1}, s _ {2}} \dots t _ {s _ {k - 1}, s _ {k}} f \left[ \lambda_ {s _ {0}}, \dots , \lambda_ {s _ {k}} \right], \tag {9.1.10}
$$

where $S _ { i j }$ is the set of all strictly increasing sequences of integers that start at i and end at j, and $f \left[ \lambda _ { s _ { 0 } } , \ldots , \lambda _ { s _ { k } } \right]$ is the kth order divided difference of f at $\{ \lambda _ { s _ { 0 } } , \ldots , \lambda _ { s _ { k } } \}$ .

Proof. See Descloux (1963), Davis (1973), or Van Loan (1975).

To illustrate the theorem, if

$$
T = \left[ \begin{array}{c c c} \lambda_ {1} & t _ {1 2} & t _ {1 3} \\ 0 & \lambda_ {2} & t _ {2 3} \\ 0 & 0 & \lambda_ {3} \end{array} \right]
$$

then

$$
f (T) = \left[ \begin{array}{c c c} f (\lambda_ {1}) & t _ {1 2} \cdot \frac {f (\lambda_ {2}) - f (\lambda_ {1})}{\lambda_ {2} - \lambda_ {1}} & F _ {1 3} \\ 0 & f (\lambda_ {2}) & t _ {2 3} \cdot \frac {f (\lambda_ {3}) - f (\lambda_ {2})}{\lambda_ {3} - \lambda_ {2}} \\ 0 & 0 & f (\lambda_ {3}) \end{array} \right],
$$

where

$$
F _ {1 3} = t _ {1 3} \cdot \frac {f (\lambda_ {3}) - f (\lambda_ {1})}{\lambda_ {3} - \lambda_ {1}} + t _ {1 2} t _ {2 3} \cdot \frac {\frac {f (\lambda_ {3}) - f (\lambda_ {2})}{\lambda_ {3} - \lambda_ {2}} - \frac {f (\lambda_ {2}) - f (\lambda_ {1})}{\lambda_ {2} - \lambda_ {1}}}{\lambda_ {3} - \lambda_ {1}}.
$$

The recipes for the upper triangular entries get increasing complicated as we move away from the diagonal. Indeed, if we explicitly use (9.1.10) to evaluate $f ( T )$ , then $O ( 2 ^ { n } )$ flops are required. However, Parlett (1974) has derived an elegant recursive method for determining the strictly upper triangular portion of the matrix $F = f ( T )$ . It requires only $2 n ^ { 3 } / 3$ flops and can be derived from the commutivity equation $F T = T F$ . Indeed, by comparing $( i , j )$ entries in this equation, we find

$$
\sum_ {k = i} ^ {j} f _ {i k} t _ {k j} = \sum_ {k = i} ^ {j} t _ {i k} f _ {k j}, \qquad j > i,
$$

and thus, if $t _ { i i }$ and $t _ { j j }$ are distinct,

$$
f _ {i j} = t _ {i j} \frac {f _ {j j} - f _ {i i}}{t _ {j j} - t _ {i i}} + \sum_ {k = i + 1} ^ {j - 1} \frac {t _ {i k} f _ {k j} - f _ {i k} t _ {k j}}{t _ {j j} - t _ {i i}}. \tag {9.1.11}
$$

From this we conclude that $f _ { i j }$ is a linear combination of its neighbors in the matrix $F$ that are to its left and below. For example, the entry $f _ { 2 5 }$ depends upon $f _ { 2 2 }$ , f23, $f _ { 2 4 } , f _ { 5 5 } , f _ { 4 5 }$ , and $f _ { 3 5 }$ . Because of this, the entire upper triangular portion of F can be computed superdiagonal by superdiagonal beginning with diag $( f ( t _ { 1 1 } ) , \ldots , f ( t _ { n n } ) )$ . The complete procedure is as follows:

Algorithm 9.1.1 (Schur-Parlett) This algorithm computes the matrix function $F =$ $f ( T )$ where T is upper triangular with distinct eigenvalues and f is defined on $\lambda ( T )$ .

for i = 1:n

$$
f _ {i i} = f (t _ {i i})
$$

end

for p = 1:n − 1

for i = 1:n − p

$$
j = i + p
$$

$$
s = t _ {i j} (f _ {j j} - f _ {i i})
$$

for k = i + 1:j − 1

$$
s = s + t _ {i k} f _ {k j} - f _ {i k} t _ {k j}
$$

end

$$
f _ {i j} = s / (t _ {j j} - t _ {i i})
$$

end

end

This algorithm requires $2 n ^ { 3 } / 3$ flops. Assuming that $A = Q T Q ^ { H }$ is the Schur decomposition of A, $f ( A ) = Q F Q ^ { H }$ where $F = f ( T )$ . Clearly, most of the work in computing $f ( A )$ by this approach is in the computation of the Schur decomposition, unless f is extremely expensive to evaluate.

# 9.1.5 A Block Schur-Parlett Approach

If A has multiple or nearly multiple eigenvalues, then the divided differences associated with Algorithm 9.1.1 become problematic and it is advisable to use a block version of the method. We outline such a procedure due to Parlett (1974). The first step is to choose Q in the Schur decomposition so that we have a partitioning

$$
T = \left[ \begin{array}{c c c c} T _ {1 1} & T _ {1 2} & \dots & T _ {1 p} \\ 0 & T _ {2 2} & \dots & T _ {2 p} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & T _ {p p} \end{array} \right]
$$

where $\lambda ( T _ { i i } ) \cap \lambda ( T _ { j j } ) ~ = ~ \emptyset$ and each diagonal block is associated with an eigenvalue cluster. The methods of §7.6 are applicable for this stage of the calculation.

Partition $F = f ( T )$ conformably

$$
F = \left[ \begin{array}{c c c c} F _ {1 1} & F _ {1 2} & \dots & F _ {1 p} \\ 0 & F _ {2 2} & \dots & F _ {2 p} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & F _ {p p} \end{array} \right],
$$

and notice that

$$
F _ {i i} = f (T _ {i i}), \qquad i = 1: p.
$$

Since the eigenvalues of $T _ { i i }$ are clustered, these calculations require special methods. Some possibilities are discussed in the next section.

Once the diagonal blocks of F are known, the blocks in the strict upper triangle of $F$ can be found recursively, as in the scalar case. To derive the governing equations, we equate $( i , j )$ blocks in $F T = T F$ for $i < j$ and obtain the following generalization of (9.1.11):

$$
F _ {i j} T _ {j j} - T _ {i i} F _ {i j} = T _ {i j} F _ {j j} - F _ {i i} T _ {i j} + \sum_ {k = i + 1} ^ {j - 1} (T _ {i k} F _ {k j} - F _ {i k} T _ {k j}). \tag {9.1.12}
$$

This is a Sylvester system whose unknowns are the elements of the block $F _ { i j }$ and whose right-hand side is “known” if we compute the $F _ { i j }$ one block superdiagonal at a time. We can solve (9.1.12) using the Bartels-Stewart algorithm (Algorithm 7.6.2). For more details see Higham (FOM, Chap. 9).

# 9.1.6 Sensitivity of Matrix Functions

Does the Schur-Parlett algorithm avoid the pitfalls associated with the diagonalization approach when the matrix of eigenvectors is ill-conditioned? The proper comparison of the two solution frameworks requires an appreciation for the notion of condition as applied to the $f ( A )$ problem. Toward that end we define the relative condition of $f$ at matrix $A \in \mathbb { C } ^ { n \times n }$ is given as

$$
\operatorname{cond} _ {\operatorname{rel}} (f, A) = \lim _ {\epsilon \to 0} \sup _ {\| E \| \leq \epsilon \| A \|} \frac {\| f (A + E) - f (A) \|}{\epsilon \| f (A) \|}.
$$

This quantity is essentially a normalized Frechet derivative of the mapping $A \to f ( A )$ and various heuristic methods have been developed for estimating its value.

It turns out that the careful implementation of the block Schur-Parlett algorithm is usually forward stable in the sense that

$$
\frac {\parallel \hat {F} - f (A) \parallel}{\parallel f (A) \parallel} \approx \mathbf {u} \cdot \mathrm{cond} _ {\mathrm{rel}} (f, A)
$$

where $\hat { F }$ is the computed version of $f ( A )$ . The same cannot be said of the diagonalization framework when the matrix of eigenvectors is ill-conditioned. For more details, see Higham (FOM, Chap. 3).

# Problems

P9.1.1 Suppose

$$
A = \left[ \begin{array}{c c} \lambda & \mu_ {1} \\ \mu_ {2} & \lambda \end{array} \right], \qquad \mu_ {1} \mu_ {2} <   0.
$$

Use the power series definitions to develop closed form expressions for exp(A), sin(A), and cos(A).

P9.1.2 Rewrite Algorithm 9.1.1 so that f(T ) is computed column by column.

P9.1.3 Suppose $A = X \mathrm { d i a g } ( \lambda _ { i } ) X ^ { - 1 }$ where $X = \ [ \ x _ { 1 } \ | \cdot \cdot \cdot | \ x _ { n } \ ]$ and $X ^ { - 1 } = [  y _ { 1 } | \cdot \cdot \cdot | y _ { n } ] ^ { H }$ . Show that if f (A) is defined, then

$$
f (A) = \sum_ {k = 1} ^ {n} f (\lambda_ {i}) x _ {i} y _ {i} ^ {H}.
$$

P9.1.4 Show that

$$
T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {q} ^ {p} \quad \Rightarrow \quad f (T) = \left[ \begin{array}{c c} F _ {1 1} & F _ {1 2} \\ 0 & F _ {2 2} \end{array} \right] _ {q} ^ {p}
$$

where $F _ { 1 1 } = f ( T _ { 1 1 } )$ and $F _ { 2 2 } = f ( T _ { 2 2 } )$ . Assume f (T ) is defined.

# Notes and References for §9.1

As we discussed, other definitions of f (A) are possible. However, for the matrix functions typically encountered in practice, all these definitions are equivalent, see:

R.F. Rinehart (1955). “The Equivalence of Definitions of a Matric Function,” Amer. Math. Monthly 62, 395–414.

The following papers are concerned with the Schur decomposition and its relationship to the f(A) problem:

C. Davis (1973). “Explicit Functional Calculus,” Lin. Alg. Applic. 6, 193–199.

J. Descloux (1963). “Bounds for the Spectral Norm of Functions of Matrices,” Numer. Math. 5, 185–190.

C.F. Van Loan (1975). “A Study of the Matrix Exponential,” Numerical Analysis Report No. 10, Department of Mathematics, University of Manchester, England. Available as Report 2006.397 from http://eprints.ma.man.ac.uk/.

Algorithm 9.1.1 and the various computational difficulties that arise when it is applied to a matrix having close or repeated eigenvalues are discuss

B.N. Parlett (1976). “A Recurrence among the Elements of Functions of Triangular Matrices,” Lin. Alg. Applic. 14, 117–121.

P.I. Davies and N.J. Higham (2003). “A Schur-Parlett Algorithm for Computing Matrix Functions,” SIAM J. Matrix Anal. Applic. 25, 464–485.

A compromise between the Jordan and Schur approaches to the f(A) problem results if A is reduced to block diagonal form as described in 7.6.3, see:

B. K˚agstr¨om (1977). “Numerical Computation of Matrix Functions,” Department of Information Processing Report UMINF-58.77, University of Ume˚ ¨ a, Sweden.

E.B. Davies (2007). “Approximate Diagonalization,” SIAM J. Matrix Anal. Applic. 29, 1051–1064.

The sensitivity of matrix functions to perturbation is discussed in:

C.S. Kenney and A.J. Laub (1989). “Condition Estimates for Matrix Functions,” SIAM J. Matrix Anal. Applic. 10, 191–209.

C.S. Kenney and A.J. Laub (1994). “Small-Sample Statistical Condition Estimates for General Matrix Functions,” SIAM J. Sci. Comput. 15, 36–61.

R. Mathias (1995). “Condition Estimation for Matrix Functions via the Schur Decomposition,” SIAM J. Matrix Anal. Applic. 16, 565–578.

# 9.2 Approximation Methods

We now consider a class of methods for computing matrix functions which at first glance do not appear to involve eigenvalues. These techniques are based on the idea that, if $g ( z )$ approximates f(z) on $\lambda ( A )$ , then f (A) approximates g(A), e.g.,

$$
e ^ {A} \approx I + A + \frac {A ^ {2}}{2 !} + \dots + \frac {A ^ {q}}{q !}.
$$

We begin by bounding $\parallel f ( A ) - g ( A ) \parallel$ using the Jordan and Schur matrix function representations. We follow this discussion with some comments on the evaluation of matrix polynomials.

# 9.2.1 A Jordan Analysis

The Jordan representation of matrix functions (Theorem 9.1.2) can be used to bound the error in an approximant $g ( A )$ of $f ( A )$ .

Theorem 9.2.1. Assume that

$$
A = X \cdot \mathrm{diag} (J _ {1}, \ldots , J _ {q}) \cdot X ^ {- 1}
$$

is the JCF of $A \in \mathbb { C } ^ { n \times n }$ with

$$
J _ {i} = \left[ \begin{array}{c c c c c} \lambda_ {i} & 1 & \dots & \dots & 0 \\ 0 & \lambda_ {i} & 1 & \vdots & \vdots \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ \vdots & \vdots & \vdots & \ddots & 1 \\ 0 & \dots & \dots & \dots & \lambda_ {i} \end{array} \right], \qquad n _ {i} \text {-by-} n _ {i},
$$

for i = 1:q. If f (z) and $g ( z )$ are analytic on an open set containing $\lambda ( A )$ , then

$$
\| f(A) - g(A)\|_{2}\leq \kappa_{2}(X)\max_{\substack{1\leq i\leq p\\ 0\leq r\leq n_{i} - 1}}n_{i}\frac{\left|f^{(r)}(\lambda_{i}) - g^{(r)}(\lambda_{i})\right|}{r!}.
$$

Proof. Defining $h ( z ) = f ( z ) - g ( z )$ we have

$$
\| f (A) - g (A) \| _ {2} = \| X \operatorname{diag} \left(h \left(J _ {1}\right), \dots , h \left(J _ {q}\right)\right) X ^ {- 1} \| _ {2} \leq \kappa_ {2} (X) \max _ {1 \leq i \leq q} \| h \left(J _ {i}\right) \| _ {2}.
$$

Using Theorem 9.1.2 and equation (2.3.8) we conclude that

$$
\| h (J _ {i}) \| _ {2} \leq n _ {i} \max _ {0 \leq r \leq n _ {i} - 1} \frac {| h ^ {(r)} (\lambda_ {i}) |}{r !}
$$

thereby proving the theorem.

# 9.2.2 A Schur Analysis

If we use the Schur decomposition $A = Q T Q ^ { H }$ instead of the Jordan decomposition, then the norm of $T \mathrm { s }$ strictly upper triangular portion is involved in the discrepancy between $f ( A )$ and $g ( A )$ .

Theorem 9.2.2. Let $Q ^ { H } A Q = T = \mathrm { d i a g } ( \lambda _ { i } ) + N$ be the Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ , with N being the strictly upper triangular portion of T . If $f ( z )$ and $g ( z )$ are analytic on a closed convex set Ω whose interior contains $\lambda ( A )$ , then

$$
\| f (A) - g (A) \| _ {F} \leq \sum_ {r = 0} ^ {n - 1} \delta_ {r} \frac {\| | N | ^ {r} \| _ {F}}{r !}
$$

where

$$
\delta_ {r} = \sup _ {z \in \Omega} \left| f ^ {(r)} (z) - g ^ {(r)} (z) \right|.
$$

Proof. Let $h ( z ) = f ( z ) - g ( z )$ and set $H = ( h _ { i j } ) = h ( A )$ . Let $S _ { i j } ^ { ( r ) }$ denote the set of strictly increasing integer sequences $\big ( s _ { 0 } , \ldots , s _ { r } \big )$ with the property that $s _ { 0 } = i$ and $s _ { r } = j$ . Notice that

$$
S _ {i j} = \bigcup_ {r = 1} ^ {j - i} S _ {i j} ^ {(r)}
$$

and so from Theorem 9.1.3, we obtain the following for all $i < j \colon$

$$
h _ {i j} = \sum_ {r = 1} ^ {j - 1} \sum_ {s \in S _ {i j} ^ {(r)}} n _ {s _ {0}, s _ {1}} n _ {s _ {1}, s _ {2}} \dots n _ {s _ {r - 1}, s _ {r}} h \left[ \lambda_ {s _ {0}}, \ldots , \lambda_ {s _ {r}} \right].
$$

Now since Ω is convex and h analytic, we have

$$
\left| h \left[ \lambda_ {s _ {0}}, \dots , \lambda_ {s _ {r}} \right] \right| \leq \sup _ {z \in \Omega} \frac {\left| h ^ {(r)} (z) \right|}{r !} = \frac {\delta_ {r}}{r !}. \tag {9.2.1}
$$

Furthermore if $| N | ^ { r } = ( n _ { i j } ^ { ( r ) } )$ for $r \geq 1$ , then it can be shown that

$$
n _ {i j} ^ {(r)} = \left\{ \begin{array}{l l} 0, & j <   i + r, \\ \sum_ {s \in S _ {i j} ^ {(r)}} \left| n _ {s _ {0}, s _ {1}} n _ {s _ {1}, s _ {2}} \dots n _ {s _ {r - 1}, s _ {r}} \right|, & j \geq i + r. \end{array} \right. \tag {9.2.2}
$$

The theorem now follows by taking absolute values in the expression for $h _ { i j }$ and then using (9.2.1) and (9.2.2).

There can be a pronounced discrepancy between the Jordan and Schur error bounds. For example, if

$$
A = \left[ \begin{array}{c c c} -. 0 1 & 1 & 1 \\ 0 & 0 & 1 \\ 0 & 0 & . 0 1 \end{array} \right].
$$

If $f ( z ) = e ^ { z }$ and $g ( z ) = 1 + z + z ^ { 2 } / 2$ , then  $f ( A ) - g ( A ) \parallel \approx \ 1 0 ^ { - 5 }$ in either the Frobenius norm or the 2-norm. Since $\kappa _ { 2 } ( X ) \approx 1 0 ^ { 7 }$ , the error predicted by Theorem 9.2.1 is $O ( 1 )$ , rather pessimistic. On the other hand, the error predicted by the Schur decomposition approach is $O ( 1 0 ^ { - 2 } )$ .

Theorems 9.2.1 and 9.2.2 remind us that approximating a function of a nonnormal matrix is more complicated than approximating a function of a scalar. In particular, we see that if the eigensystem of A is ill-conditioned and/or A’s departure from normality is large, then the discrepancy between $f ( A )$ and $g ( A )$ may be considerably larger than the maximum of $| f ( z ) - g ( z ) |$ on $\lambda ( A )$ . Thus, even though approximation methods avoid eigenvalue computations, they evidently appear to be influenced by the structure of A’s eigensystem. It is a perfect venue for pseudospectral analysis.

# 9.2.3 Taylor Approximants

A common way to approximate a matrix function such as $e ^ { A }$ is by truncating its Taylor series. The following theorem bounds the errors that arise when matrix functions such as these are approximated via truncated Taylor series.

Theorem 9.2.3. If $f ( z )$ has the Taylor series

$$
f (z) = \sum_ {k = 0} ^ {\infty} \alpha_ {k} z ^ {k}
$$

on an open disk containing the eigenvalues of $A \in \mathbb { C } ^ { n \times n }$ , then

$$
\left\| f (A) - \sum_ {k = 0} ^ {q} \alpha_ {k} A ^ {k} \right\| _ {2} \leq \frac {n}{(q + 1) !} \max _ {0 \leq s \leq 1} \| A ^ {q + 1} f ^ {(q + 1)} (A s) \| _ {2}.
$$

Proof. Define the matrix $E ( s )$ by

$$
f (A s) = \sum_ {k = 0} ^ {q} \alpha_ {k} (A s) ^ {k} + E (s), \quad 0 \leq s \leq 1. \tag {9.2.3}
$$

If $f _ { i j } ( s )$ is the $( i , j )$ entry of $f ( A s )$ , then it is necessarily analytic and so

$$
f _ {i j} (s) = \left(\sum_ {k = 0} ^ {q} \frac {f _ {i j} ^ {(k)} (0)}{k !} s ^ {k}\right) + \frac {f _ {i j} ^ {(q + 1)} (\varepsilon_ {i j})}{(q + 1) !} s ^ {q + 1} \tag {9.2.4}
$$

where $\varepsilon _ { i j }$ satisfies $0 \leq \varepsilon _ { i j } \leq s \leq 1$

By comparing powers of s in (9.2.3) and (9.2.4) we conclude that $e _ { i j } ( s )$ , the $( i , j )$ entry of $E ( s )$ , has the form

$$
e _ {i j} (s) = \frac {f _ {i j} ^ {(q + 1)} (\varepsilon_ {i j})}{(q + 1) !} s ^ {q + 1}.
$$

Now $f _ { i j } ^ { ( q - 1 ) } ( s )$ is the $( i , j )$ entry of $A ^ { q + 1 } f ^ { ( q + 1 ) } ( A s )$ and therefore

$$
| e _ {i j} (s) | \leq \max _ {0 \leq s \leq 1} \frac {f _ {i j} ^ {(q + 1)} (s)}{(q + 1) !} \leq \max _ {0 \leq s \leq 1} \frac {\| A ^ {q + 1} f ^ {(q + 1)} (A s) \| _ {2}}{(q + 1) !}.
$$

The theorem now follows by applying (2.3.8).

We mention that the factor of n in the upper bound can be removed with more careful analysis. See Mathias (1993).

In practice, it does not follow that greater accuracy results by taking a longer Taylor approximation. For example, if

$$
A = \left[ \begin{array}{l l} - 4 9 & 2 4 \\ - 6 4 & 3 1 \end{array} \right],
$$

then it can be shown that

$$
e ^ {A} = \left[ \begin{array}{l l} - 0. 7 3 5 7 5 9 & . 0 5 5 1 8 1 9 \\ - 1. 4 7 1 5 1 8 & 1. 1 0 3 6 3 8 \end{array} \right].
$$

For q = 59, Theorem 9.2.3 predicts that

$$
\left\| e ^ {A} - \sum_ {k = 0} ^ {q} \frac {A ^ {k}}{k !} \right\| _ {2} \leq \frac {n}{(q + 1) !} \max _ {0 \leq s \leq 1} \left\| A ^ {q + 1} e ^ {A s} \right\| _ {2} \leq 1 0 ^ {- 6 0}.
$$

However, if $\mathbf { u } \approx 1 0 ^ { - 7 }$ , then we find

$$
\mathsf {f l} \left(\sum_ {k = 0} ^ {5 9} \frac {A ^ {k}}{k !}\right) = \left[ \begin{array}{c c} - 2 2. 2 5 8 8 0 & - 1. 4 3 2 2 7 6 6 \\ - 6 1. 4 9 9 3 1 & - 3. 4 7 4 2 8 0 \end{array} \right].
$$

The problem is that some of the partial sums have large elements. For example, the matrix $I + A + \cdot \cdot \cdot + A ^ { 1 7 } / 1 7 !$ has entries of order $1 0 ^ { 7 }$ . Since the machine precision is approximately $1 0 ^ { - 7 }$ , rounding errors larger than the norm of the solution are sustained.


---

<!-- golub_550_599 -->

The example highlights the a well known shortcoming of truncated Taylor series approximation–it tends to be effcetive only near the origin. The problem can sometimes be circumvented through a change of scale. For example, by repeatedly using the double angle formulae

$$
\cos (2 A) = 2 \cos (A) ^ {2} - I, \quad \sin (2 A) = 2 \sin (A) \cos (A),
$$

the cosine and sine of a matrix can be built up from Taylor approximations to $\cos ( A / 2 ^ { k } )$ and sin $( A / 2 ^ { k } )$ :

$$
S _ {0} = \text { Taylor   approximate   to } \sin (A / 2 ^ {k})
$$

$$
C _ {0} = \text { Taylor   approximate   to } \cos (A / 2 ^ {k})
$$

$$
\text { for } j = 1: k
$$

$$
S _ {j} = 2 S _ {j - 1} C _ {j - 1}
$$

$$
C _ {j} = 2 C _ {j - 1} ^ {2} - I
$$

end

Here k is a positive integer chosen so that, say, $\| A \| _ { \infty } \approx 2 ^ { k }$ . See Serbin and Blalock (1979), Higham and Smith (2003), and Hargreaves and Higham (2005).

# 9.2.4 Evaluating Matrix Polynomials

Since the approximation of transcendental matrix functions usually involves the evaluation of polynomials, it is worthwhile to look at the details of computing

$$
p (A) = b _ {0} I + b _ {1} A + \dots + b _ {q} A ^ {q}
$$

where the scalars $b _ { 0 } , \dots , b _ { q } \in \mathbb { R }$ are given. The most obvious approach is to invoke Horner’s scheme:

Algorithm 9.2.1 Given a matrix A and $b ( 0 { : } q )$ , the following algorithm computes the polynomial $F = b _ { q } A ^ { q } + \cdot \cdot \cdot + b _ { 1 } A + b _ { 0 } I$ .

$$
F = b _ {q} A + b _ {q - 1} I
$$

$$
\text { for   } k = q - 2: - 1: 0
$$

$$
F = A F + b _ {k} I
$$

end

This requires $q - 1$ matrix multiplications. However, unlike the scalar case, this summation process is not optimal. To see why, suppose $q = 9$ and observe that

$$
p (A) = A ^ {3} (A ^ {3} (b _ {9} A ^ {3} + (b _ {8} A ^ {2} + b _ {7} A + b _ {6} I)) + (b _ {5} A ^ {2} + b _ {4} A + b _ {3} I)) + b _ {2} A ^ {2} + b _ {1} A + b _ {0} I.
$$

Thus, $F = p ( A )$ can be evaluated with only four matrix multiplications:

$$
\begin{array}{l} A _ {2} = A ^ {2}, \\ A _ {3} = A A _ {2}, \\ F _ {1} = b _ {9} A _ {3} + b _ {8} A _ {2} + b _ {7} A + b _ {6} I, \\ F _ {2} = A _ {3} F _ {1} + b _ {5} A _ {2} + b _ {4} A + b _ {3} I, \\ F = A _ {3} F _ {2} + b _ {2} A _ {2} + b _ {1} A + b _ {0} I. \\ \end{array}
$$

In general, if s is any integer that satisfies $1 \leq s \leq { \sqrt { q } }$ , then

$$
p (A) = \sum_ {k = 0} ^ {r} B _ {k} \cdot (A ^ {s}) ^ {k}, \quad r = \text { floor } (q / s), \tag {9.2.5}
$$

where

$$
B _ {k} = \left\{ \begin{array}{l l} b _ {s k + s - 1} A ^ {s - 1} + \dots + b _ {s k + 1} A + b _ {s k} I, & k = 0: r - 1, \\ b _ {q} A ^ {q - s r} + \dots + b _ {s r + 1} A + b _ {s r} I, & k = r. \end{array} \right.
$$

After $A ^ { 2 } , \ldots , A ^ { s }$ are computed, then Horner’s rule can be applied to (9.2.5) and the net result is that $p ( A )$ can be computed with $s + r - 1$ matrix multiplications. By choosing $s = \mathsf { f l o o r } ( { \sqrt { q } } )$ , the number of matrix multiplications is approximately minimized. This technique is discussed by Paterson and Stockmeyer (1973). Van Loan (1978) shows how the procedure can be implemented without storage arrays for $A ^ { 2 } , \ldots , A ^ { s }$ .

# 9.2.5 Computing Powers of a Matrix

The problem of raising a matrix to a given power deserves special mention. Suppose it is required to compute $A ^ { 1 3 }$ . Noting that $A ^ { \bar { 4 } } = ( A ^ { 2 } ) ^ { 2 } , A ^ { 8 } = ( A ^ { 4 } ) ^ { 2 }$ , and $A ^ { 1 3 } = \bar { A } ^ { 8 } A ^ { 4 } A$ , we see that this can be accomplished with just five matrix multiplications. In general we have

Algorithm 9.2.2 (Binary Powering) The following algorithm computes $F = A ^ { s }$ where s is a positive integer and A ∈ IRn×n. $A \in \mathbb { R } ^ { n \times n }$

Let $s = \sum _ { k = 0 } ^ { t } \beta _ { k } 2 ^ { k }$ be the binary expansion of s with $\beta _ { t } \neq 0$

$$
Z = A; q = 0
$$

while $\beta _ { q } = 0$

$$
Z = Z ^ {2}; q = q + 1
$$

end

$$
F = Z
$$

for $k = q + 1 { : } t$

$$
Z = Z ^ {2}
$$

$\mathbf { i f } \ \beta _ { k } \neq 0$

$$
F = F Z
$$

end

end

This algorithm requires at most $2 \mathsf { f l o o r } [ \log _ { 2 } ( s ) ]$ matrix multiplications. If s is a power of 2, then only $\log _ { 2 } ( s )$ matrix multiplications are needed.

# 9.2.6 Integrating Matrix Functions

We conclude this section with some remarks about the integration of a parameterized matrix function. Suppose $A \in \mathbb { R } ^ { n \times n }$ and that $f ( A t )$ is defined for all $t \in [ a , b ]$ . We can

approximate

$$
F = \int_ {a} ^ {b} f (A t) d t \quad \Leftrightarrow \quad [ F ] _ {i j} = \int_ {a} ^ {b} [ f (A t) ] _ {i j} d t
$$

by applying any suitable quadrature rule. For example, with Simpson’s rule, we have

$$
F \approx \tilde {F} = \frac {h}{3} \sum_ {k = 0} ^ {m} w _ {k} f (A (a + k h)) \tag {9.2.6}
$$

where m is even, $h = ( b - a ) / m$ , and

$$
w _ {k} = \left\{ \begin{array}{l l} 1 & k = 0, m, \\ 4 & k \text {odd}, \\ 2 & k \text {even}, k \neq 0, m. \end{array} \right.
$$

If $( d ^ { 4 } / d z ^ { 4 } ) f ( z t ) = f ^ { ( 4 ) } ( z t )$ is continuous for $t \in [ a , b ]$ and if $f ^ { ( 4 ) } ( A t )$ is defined on this same interval, then it can be shown that $\tilde { \boldsymbol { F } } = \boldsymbol { F } + \boldsymbol { E }$ where

$$
\| E \| _ {2} \leq \frac {n h ^ {4} (b - a)}{1 8 0} \max _ {a \leq t \leq b} \| f ^ {(4)} (A t) \| _ {2}. \tag {9.2.7}
$$

Let $f _ { i j }$ and $e _ { i j }$ denote the $( i , j )$ entries of F and E, respectively. Under the above assumptions we can apply the standard error bounds for Simpson’s rule and obtain

$$
| e _ {i j} | \leq \frac {h ^ {4} (b - a)}{1 8 0} \max _ {a \leq t \leq b} | e _ {i} ^ {T} f ^ {(4)} (A t) e _ {j} |.
$$

The inequality (9.2.7) now follows since $\parallel E \parallel _ { 2 } \leq n \operatorname* { m a x } | e _ { i j } |$ and

$$
\max _ {a \leq t \leq b} | e _ {i} ^ {T} f ^ {(4)} (A t) e _ {j} | \leq \max _ {a \leq t \leq b} \| f ^ {(4)} (A t) \| _ {2}.
$$

Of course, in a practical application of (9.2.6), the function evaluations $f ( A ( a + k h ) )$ normally have to be approximated. Thus, the overall error involves the error in approximating $f ( A ( a + k h )$ as well as the Simpson rule error.

# 9.2.7 A Note on the Cauchy Integral Formulation

Yet another way to define a function of a matrix $C \in \mathbb { C } ^ { n \times n }$ is through the Cauchy integral theorem. Suppose $f ( z )$ is analytic inside and on a closed contour Γ which encloses $\lambda ( A )$ . We can define $f ( A )$ to be the matrix

$$
f (A) = \frac {1}{2 \pi i} \oint_ {\Gamma} f (z) (z I - A) ^ {- 1} d z. \tag {9.2.8}
$$

The integral is defined on an element-by-element basis:

$$
f (A) = \left(f _ {k j}\right) \quad \Longrightarrow \quad f _ {k j} = \frac {1}{2 \pi i} \oint_ {\Gamma} f (z) e _ {k} ^ {T} (z I - A) ^ {- 1} e _ {j} d z.
$$

Notice that the entries of $( z I { - } A ) ^ { - 1 }$ are analytic on Γ and that $f ( A )$ is defined whenever $f ( z )$ is analytic in a neighborhood of $\lambda ( A )$ . Using quadrature and other tools, Hale, Higham, and Trefethen (2007) have shown how this characterization can be used in practice to compute certain types of matrix functions.

# Problems

P9.2.1 Verify (9.2.2).

P9.2.2 Show that if $\| \ A \| _ { 2 } < 1$ , then $\log ( I + A )$ exists and satisfies the bound

$$
\| \log (I + A) \| _ {2} \leq \| A \| _ {2} / (1 - \| A \| _ {2}).
$$

P9.2.3 Using Theorem 9.2.3, bound the error in the following approximations:

$$
\sin (A) \approx \sum_ {k = 0} ^ {q} (- 1) ^ {k} \frac {A ^ {2 k + 1}}{(2 k + 1) !}, \qquad \cos (A) \approx \sum_ {k = 0} ^ {q} (- 1) ^ {k} \frac {A ^ {2 k}}{(2 k) !}.
$$

P9.2.4 Suppose $A \in \mathbb { R } ^ { n \times n }$ is nonsingular and $X _ { 0 } \in \mathbb { R } ^ { n \times n }$ is given. The iteration defined by

$$
X _ {k + 1} = X _ {k} (2 I - A X _ {k})
$$

is the matrix analogue of Newton’s method applied to the function $f ( x ) = a - ( 1 / x )$ . Use the SVD to analyze this iteration. Do the iterates converge to $A ^ { - 1 \mathord { \ ? } }$ Discuss the choice of $X _ { 0 }$ .

P9.2.5 Assume $A \in \mathbb { R } ^ { 2 \times 2 }$ . (a) Specify real scalars α and $\beta$ so that $A ^ { 4 } = \alpha I + \beta A$ . (b) Develop recursive recipes for $\alpha _ { k }$ and $\beta _ { k }$ so that $A ^ { k } = \alpha _ { k } I + \beta _ { k } A$ for $k \geq 2$ .

# Notes and References for 9.2

The optimality of Horner’s rule for polynomial evaluation is discussed in:

M.S. Paterson and L.J. Stockmeyer (1973). “On the Number of Nonscalar Multiplications Necessary to Evaluate Polynomials,” SIAM J. Comput. 2, 60–66.   
D.E. Knuth (1981). The Art of Computer Programming, Vol. 2. Seminumerical Algorithms, second edition, Addison-Wesley, Reading, MA.   
The Horner evaluation of matrix polynomials is analyzed in:   
C.F. Van Loan (1978). “A Note on the Evaluation of Matrix Polynomials,” IEEE Trans. Autom. Control AC-24, 320–321.   
Other aspects of matrix function approximation and evaluation are discussed in:   
H. Bolz and W. Niethammer (1988). “On the Evaluation of Matrix Functions Given by Power Series,” SIAM J. Matrix Anal. Applic. 9, 202–209.   
R. Mathias (1993). “Approximation of Matrix-Valued Functions,” SIAM J. Matrix Anal. Applic. 14, 1061–1063.   
N.J. Higham and P.A. Knight (1995). “Matrix Powers in Finite Precision Arithmetic,” SIAM J. Matrix Anal. Applic. 16, 343–358.   
P. Sebastiani (1996). “On the Derivatives of Matrix Powers,” SIAM J. Matrix Anal. Applic. 17, 640–648.   
D.S. Bernstein and C.F. Van Loan (2000). “Rational Matrix Functions and Rank-One Updates,” SIAM J. Matrix Anal. Applic. 22, 145–154.   
For a discussion of methods for computing the sine and cosine of a matrix, see:   
S. Serbin and S. Blalock (1979). “An Algorithm for Computing the Matrix Cosine,” SIAM J. Sci. Stat. Comput. 1, 198–204.   
N.J. Higham and M.I. Smit (2003). “Computing the Matrix Cosine,” Numer. Algorithms 34, 13–26.   
G. Hargreaves and N.J. Higham (2005). “Efficient Algorithms for the Matrix Cosine and Sine,” Numer. Algorithms 40, 383–400.   
The computation of f(A) using contour integrals is analyzed in:   
N. Hale, N.J. Higham, and L.N. Trefethen (2007). “Computing Aα, log(A), and Related Matrix Functions by Contour Integrals,” SIAM J. Numer. Anal. 46, 2505–2523.

# 9.3 The Matrix Exponential

One of the most frequently computed matrix functions is the exponential

$$
e ^ {A t} = \sum_ {k = 0} ^ {\infty} \frac {(A t) ^ {k}}{k !}.
$$

Numerous algorithms for computing $e ^ { \boldsymbol { A } t }$ have been proposed, but most of them are of dubious numerical quality, as is pointed out in the survey articles by Moler and Van Loan (1978) and its update Moler and Van Loan (2003). In order to illustrate what the computational difficulties are, we present a “scaling and squaring” method based upon Pad´e approximation. A brief analysis of the method follows that involves some $e ^ { \boldsymbol { A } t }$ perturbation theory and includes comments about the shortcomings of eigenanalysis in settings where nonnormality prevails.

# 9.3.1 A Pad´e Approximation Method

Following the discussion in §9.2, if $g ( z ) \approx e ^ { z }$ , then $g ( A ) \approx e ^ { A }$ . A very useful class of approximants for this purpose are the Pad´e functions defined by

$$
R _ {p q} (z) = D _ {p q} (z) ^ {- 1} N _ {p q} (z),
$$

where

$$
N _ {p q} (z) = \sum_ {k = 0} ^ {p} \frac {(p + q - k) ! p !}{(p + q) ! k ! (p - k) !} z ^ {k}
$$

and

$$
D _ {p q} (z) = \sum_ {k = 0} ^ {q} \frac {(p + q - k) ! q !}{(p + q) ! k ! (q - k) !} (- z) ^ {k}.
$$

Notice that

$$
R _ {p o} (z) = 1 + z + \dots + z ^ {p} / p!
$$

is the order-p Taylor polynomial.

Unfortunately, the Pad´e approximants are good only near the origin, as the following identity reveals:

$$
e ^ {A} = R _ {p q} (A) + \frac {(- 1) ^ {q}}{(p + q) !} A ^ {p + q + 1} D _ {p q} (A) ^ {- 1} \int_ {0} ^ {1} u ^ {p} (1 - u) ^ {q} e ^ {A (1 - u)} d u. \tag {9.3.1}
$$

However, this problem can be overcome by exploiting the fact that

$$
e ^ {A} = (e ^ {A / m}) ^ {m}.
$$

In particular, we can scale A by m such that $F _ { p q } = R _ { p q } ( A / m )$ is a suitably accurate approximation to $e ^ { \boldsymbol { A } / m }$ . We then compute $F _ { p q } ^ { m }$ using Algorithm 9.2.2. If m is a power of two, then this amounts to repeated squaring and so is very efficient. The success of the overall procedure depends on the accuracy of the approximant

$$
F _ {p q} = \left(R _ {p q} \left(\frac {A}{2 ^ {j}}\right)\right) ^ {2 ^ {j}}.
$$

In Moler and Van Loan (1978) it is shown that, if

$$
\frac {\| A \| _ {\infty}}{2 ^ {j}} \leq \frac {1}{2},
$$

then there exists an $E \in \mathbb { R } ^ { n \times n }$ such that $F _ { p q } = e ^ { A + E } , A E = E A$ , and

$$
\| E \| _ {\infty} \leq \varepsilon (p, q) \| A \| _ {\infty},
$$

where

$$
\varepsilon (p, q) = 2 ^ {3 - (p + q)} \frac {p ! q !}{(p + q) ! (p + q + 1) !}.
$$

Using these results it is easy to establish the inequality

$$
\frac {\| e ^ {A} - F _ {p q} \| _ {\infty}}{\| e ^ {A} \| _ {\infty}} \leq \epsilon (p, q) \| A \| _ {\infty} e ^ {\epsilon (p, q) \| A \| _ {\infty}}.
$$

The parameters $p$ and $q$ can be determined according to some relative error tolerance. Since $F _ { p q }$ requires about $j + \operatorname* { m a x } \{ p , q \}$ matrix multiplications, it makes sense to set p $= q$ as this choice minimizes $\epsilon ( p , q )$ for a given amount of work. Overall we obtain

Algorithm 9.3.1 (Scaling and Squaring) Given $\delta > 0$ and $A \in \mathbb { R } ^ { n \times n }$ , the following algorithm computes $F = \stackrel { - } { e } ^ { A + E }$ where $\parallel E \parallel _ { \infty } \leq \delta \parallel A \parallel _ { \infty }$ .

$$
j = \max \{0, 1 + \operatorname{floor} (\log_ {2} (\| A \| _ {\infty})) \}
$$

$$
A = A / 2 ^ {j}
$$

Let $q$ be the smallest nonnegative integer such that $\epsilon ( q , q ) \leq \delta$

$$
D = I, N = I, X = I, c = 1
$$

for $k = 1 { : } q$

$$
c = c \cdot (q - k + 1) / ((2 q - k + 1) k)
$$

$$
X = A X, N = N + c \cdot X, D = D + (- 1) ^ {k} c \cdot X
$$

end

Solve $D F = N$ for $F$ using Gaussian elimination

for $k = 1 { : } j$

$$
F = F ^ {2}
$$

end

This algorithm requires about $2 ( q + j + 1 / 3 ) n ^ { 3 }$ flops. Its roundoff error properties of have been analyzed by Ward (1977). For further analysis and algorithmic improvements, see Higham (2005) and Al-Mohy and Higham (2009).

The special Horner techniques of §9.2.4 can be applied to quicken the computation of $D = D _ { q q } ( A )$ and $N = N _ { q q } ( A )$ . For example, if $q = 8$ we have $N _ { q q } ( A ) = U + A V$ and $D _ { q q } ( A ) = U - A V$ where

$$
U = c _ {0} I + c _ {2} A ^ {2} + (c _ {4} I + c _ {6} A ^ {2} + c _ {8} A ^ {4}) A ^ {4}
$$

and

$$
V = c _ {1} I + c _ {3} A ^ {2} + (c _ {5} I + c _ {7} A ^ {2}) A ^ {4}.
$$

Clearly, N and D can be computed with five matrix multiplications instead of seven as required by Algorithm 9.3.1.

# 9.3.2 Perturbation Theory

Is Algorithm 9.3.1 stable in the presence of roundoff error? To answer this question we need to understand the sensitivity of the matrix exponential to perturbations in A. The rich structure of this particular matrix function enables us to say more about the condition of the $e ^ { A }$ problem than is typically the case for a general matrix function. (See §9.1.6.)

The starting point in the discussion is the initial value problem

$$
\dot {X} (t) = A X (t), \qquad X (0) = I,
$$

where A, $\ b X ( t ) \in \mathbb R ^ { n \times n }$ . This has the unique solution $X ( t ) = e ^ { A t }$ , a characterization of the matrix exponential that can be used to establish the identity

$$
e ^ {(A + E) t} - e ^ {A t} = \int_ {0} ^ {t} e ^ {A (t - s)} E e ^ {(A + E) s} d s.
$$

From this it follows that

$$
\frac {\parallel e ^ {(A + E) t} - e ^ {A t} \parallel_ {2}}{\parallel e ^ {A t} \parallel_ {2}} \leq \frac {\parallel E \parallel_ {2}}{\parallel e ^ {A t} \parallel_ {2}} \int_ {0} ^ {t} \parallel e ^ {A (t - s)} \parallel_ {2} \parallel e ^ {(A + E) s} \parallel_ {2} d s.
$$

Further simplifications result if we bound the norms of the exponentials that appear in the integrand. One way of doing this is through the Schur decomposition. If $Q ^ { H } A Q \ =$ $\mathrm { d i a g } ( \lambda _ { i } ) + N$ is the Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ , then it can be shown that

$$
\| e ^ {A t} \| _ {2} \leq e ^ {\alpha (A) t M _ {S} (t)}, \tag {9.3.2}
$$

where

$$
\alpha (A) = \max \left\{\operatorname{Re} (\lambda): \lambda \in \lambda (A) \right\} \tag {9.3.3}
$$

is the spectral abscissa and

$$
M _ {S} (t) = \sum_ {k = 0} ^ {n - 1} \frac {\parallel N t \parallel_ {2} ^ {k}}{k !}.
$$

With a little manipulation it can be shown that

$$
\frac {\parallel e ^ {(A + E) t} - e ^ {A t} \parallel_ {2}}{\parallel e ^ {A t} \parallel_ {2}} \leq t \parallel E \parallel_ {2} M _ {S} (t) ^ {2} \exp (t M _ {S} (t) \parallel E \parallel_ {2}).
$$

Notice that $M _ { s } ( t ) \equiv 1$ if and only if A is normal, suggesting that the matrix exponential problem is “well-behaved” if A is normal. This observation is confirmed by the behavior of the matrix exponential condition number $\nu ( A , t )$ , defined by

$$
\nu (A, t) = \max _ {\| E \| \leq 1} \left\| \int_ {0} ^ {t} e ^ {A (t - s)} E e ^ {A s} d s \right\| _ {2} \frac {\| A \| _ {2}}{\| e ^ {A t} \| _ {2}}.
$$

This quantity, discussed by Van Loan (1977), measures the sensitivity of the map $A  e ^ { A t }$ in that for a given t, there is a matrix E for which

$$
\frac {\parallel e ^ {(A + E) t} - e ^ {A t} \parallel_ {2}}{\parallel e ^ {A t} \parallel_ {2}} \approx \nu (A, t) \frac {\parallel E \parallel_ {2}}{\parallel A \parallel_ {2}}.
$$

![](images/golub_550_599__13701d27df41c2cee17d499bfd734a0bca95b23cec76e74d9e672a10efe423c2.jpg)  
Figure 9.3.1. $\parallel e ^ { \boldsymbol { A } t } \parallel _ { 2 }$ can grow even if $\alpha ( A ) < 0$

Thus, if $\nu ( A , t )$ is large, small changes in A can induce relatively large changes in $e ^ { \boldsymbol { A } t }$ . Unfortunately, it is difficult to characterize precisely those A for which $\nu ( A , t )$ is large. (This is in contrast to the linear equation problem Ax = b, where the ill-conditioned A are neatly described in terms of SVD.) One thing we can say, however, is that $\nu ( A , t ) \geq t \| A \| _ { 2 }$ , with equality holding for all nonnegative t if and only if the matrix A is normal.

# 9.3.3 Pseudospectra

Dwelling a little more on the effect of nonnormality, we know from the analysis of §9.2 that approximating $e ^ { \boldsymbol { A } t }$ involves more than just approximating $e ^ { z t }$ on $\lambda ( A )$ . Another clue that eigenvalues do not “tell the whole story” in the $e ^ { \ b { A } t }$ problem has to do with the inability of the spectral abscissa (9.3.3) to predict the size of $\parallel e ^ { \boldsymbol { A } t } \parallel _ { 2 }$ as a function of time. If A is normal, then

$$
\| e ^ {A t} \| _ {2} = e ^ {\alpha (A) t}. \tag {9.3.4}
$$

Thus, there is uniform decay if the eigenvalues of A are in the open left half plane. But if A is non-normal, then $e ^ { \boldsymbol { A } \dot { \boldsymbol { t } } }$ can grow before decay sets in. The 2-by-2 example

$$
A = \left[ \begin{array}{c c} - 1 & 1 0 0 0 \\ 0 & - 1 \end{array} \right] \quad \Leftrightarrow \quad e ^ {A t} = e ^ {- t} \left[ \begin{array}{c c} 1 & 1 0 0 0 \cdot t \\ 0 & 1 \end{array} \right] \tag {9.3.5}
$$

plainly illustrates this point in Figure 9.3.1.

Pseudospectra can be used to shed light on the transient growth of $\parallel e ^ { A t } \parallel$ . For example, it can be shown that for every $\epsilon > 0$ ,

$$
\sup _ {t > 0} \| e ^ {A t} \| _ {2} \geq \frac {\alpha_ {\epsilon} (A)}{\epsilon} \tag {9.3.6}
$$

where $\alpha _ { \epsilon } ( A )$ is the 
-pseudospectral abscissa introduced in (7.8.8):

$$
\alpha_ {\epsilon} (A) = \sup _ {z \in \Lambda_ {\epsilon} (A)} \operatorname{Re} (z).
$$

For the 2-by-2 matrix in (9.3.5), it can be shown that $\alpha _ { . 0 1 } ( A ) / . 0 1 \approx 2 1 6$ , a value that is consistent with the growth curve in Figure 9.3.1. See Trefethen and Embree (SAP, Chap. 15) for more pseudospectral insights into the behavior of $\parallel e ^ { \boldsymbol { A } t } \parallel _ { 2 }$

# 9.3.4 Some Stability Issues

With this discussion we are ready to begin thinking about the stability of Algorithm 9.3.1. A potential difficulty arises during the squaring process if A is a matrix whose exponential grows before it decays. If

$$
G = R _ {q q} \left(\frac {A}{2 ^ {j}}\right) \approx e ^ {A / 2 ^ {j}},
$$

then it can be shown that rounding errors of order

$$
\gamma = \mathbf {u} \| G ^ {2} \| _ {2} \cdot \| G ^ {4} \| _ {2} \cdot \| G ^ {8} \| _ {2} \dots \| G ^ {2 ^ {j - 1}} \| _ {2}
$$

can be expected to contaminate the computed $G ^ { 2 ^ { j } } . \mathrm { ~ H ~ } \parallel \boldsymbol { e } ^ { A t } \parallel _ { 2 }$ has a substantial initial growth, then it may be the case that

$$
\gamma \gg \mathbf {u} \| G ^ {2 ^ {j}} \| _ {2} \approx \mathbf {u} \| e ^ {A} \| _ {2},
$$

thus ruling out the possibility of small relative errors.

If A is normal, then so is the matrix G and therefore $\parallel G ^ { m } \parallel _ { 2 } = \parallel G \parallel _ { 2 } ^ { m }$ for all positive integers m. Thus, $\gamma \approx \mathbf { u } \parallel G ^ { 2 ^ { j } } \parallel _ { 2 } \approx \mathbf { u } \parallel e ^ { A } \parallel _ { 2 }$ and so the initial growth problems disappear. The algorithm can essentially be guaranteed to produce small relative error when A is normal. On the other hand, it is more difficult to draw conclusions about the method when A is nonnormal because the connection between $\nu ( A , t )$ and the initial growth phenomena is unclear. However, numerical experiments suggest that Algorithm 9.3.1 fails to produce a relatively accurate $e ^ { A }$ only when v(A, 1) is correspondingly large.

# Problems

P9.3.1 Show that $e ^ { ( A + B ) t } = e ^ { A t } e ^ { B t }$ for all t if and only if AB = BA. Hint: Express both sides as a power series in t and compare the coefficient of t.

P9.3.2 Suppose that A is skew-symmetric. Show that both $e ^ { A }$ and the (1,1) Pad´e approximatant $R _ { 1 1 } ( A )$ are orthogonal. Are there any other values of p and q for which $R _ { p q } ( A )$ is orthogonal?

P9.3.3 Show that if A is nonsingular, then there exists a matrix X such that $A = e ^ { X }$ . Is X unique?

P9.3.4 Show that if

$$
\exp \left(\left[ \begin{array}{c c} - A ^ {T} & P \\ 0 & A \end{array} \right] z\right) = \left[ \begin{array}{c c} F _ {1 1} & F _ {1 2} \\ 0 & F _ {2 2} \end{array} \right] _ {n} ^ {n}
$$

then

$$
F _ {1 1} ^ {T} F _ {1 2} = \int_ {0} ^ {z} e ^ {A ^ {T} t} P e ^ {A t} d t.
$$

P9.3.5 Give an algorithm for computing $e ^ { A }$ when $A = u v ^ { T } , u , v \in \mathbb { R } ^ { n }$ .

P9.3.6 Suppose $A \in \mathbb { R } ^ { n \times n }$ and that $v \in \mathbb { R } ^ { n }$ has unit 2-norm. Define the function $\phi ( t ) = \parallel e ^ { A t } v \parallel _ { 2 } ^ { 2 } / 2$ and show that

$$
\dot {\phi} (t) \leq \mu (A) \phi (t)
$$

where $\mu ( A ) = \lambda _ { 1 } ( ( A + A ^ { T } ) / 2 )$ . Conclude that

$$
\parallel e ^ {A t} \parallel_ {2} \leq e ^ {\mu (A) t}
$$

where $t \geq 0 .$

P9.3.7 Suppose $A \in \mathbb { R } ^ { n \times n }$ has the property that its off-diagonal entries are negative and its column sums are zero. Show that for all t, $F = \exp ( A t )$ has nonnegative entries and unit column sums.

# Notes and References for §9.3

Much of what appears in this section and an extensive bibliography may be found in the following survey articles:

C.B. Moler and C.F. Van Loan (1978). “Nineteen Dubious Ways to Compute the Exponential of a Matrix,” SIAM Review 20, 801–836.

C.B. Moler and C.F.Van Loan (2003). “Nineteen Dubious Ways to Compute the Exponential of a Matrix, Twenty-Five Years Later,” SIAM Review 45, 3–49.

Scaling and squaring with Pad´e approximants (Algorithm 9.3.1) and a careful implementation of the Schur decomposition method (Algorithm 9.1.1) were found to be among the less dubious of the nineteen methods scrutinized. Various aspects of Pad´e approximation of the matrix exponential are discussed in:

W. Fair and Y. Luke (1970). “Pad´e Approximations to the Operator Exponential,” Numer. Math. 14, 379–382.

C.F. Van Loan (1977). “On the Limitation and Application of Pad´e Approximation to the Matrix Exponential,” in Pad´e and Rational Approximation, E.B. Saff and R.S. Varga (eds.), Academic Press, New York.

R.C. Ward (1977). “Numerical Computation of the Matrix Exponential with Accuracy Estimate,” SIAM J. Numer. Anal. 14, 600–614.

A. Wragg (1973). “Computation of the Exponential of a Matrix I: Theoretical Considerations,” J. Inst. Math. Applic. 11, 369–375.

A. Wragg (1975). “Computation of the Exponential of a Matrix II: Practical Considerations,” J. Inst. Math. Applic. 15, 273–278.

L. Dieci and A. Papini (2000). “Pad´e Approximation for the Exponential of a Block Triangular Matrix,” Lin. Alg. Applic. 308, 183–202.

M. Arioli, B. Codenotti and C. Fassino (1996). “The Pad´e Method for Computing the Matrix Exponential,” Lin. Alg. Applic. 240, 111–130.

N.J. Higham (2005). “The Scaling and Squaring Method for the Matrix Exponential Revisited,” SIAM J. Matrix Anal. Applic. 26, 1179–1193.

A.H. Al-Mohy and N.J. Higham (2009). “A New Scaling and Squaring Algorithm for the Matrix Exponential,” SIAM J. Matrix Anal. Applic. 31, 970–989.

A proof of Equation (9.3.1) for the scalar case appears in:

R.S. Varga (1961). “On Higher-Order Stable Implicit Methods for Solving Parabolic Partial Differential Equations,” J. Math. Phys. 40, 220–231.

There are many applications in control theory calling for the computation of the matrix exponential. In the linear optimal regular problem, for example, various integrals involving the matrix exponential are required, see:

J. Johnson and C.L. Phillips (1971). “An Algorithm for the Computation of the Integral of the State Transition Matrix,” IEEE Trans. Autom. Control AC-16, 204–205.

C.F. Van Loan (1978). “Computing Integrals Involving the Matrix Exponential,” IEEE Trans. Autom. Control AC-23, 395–404.

An understanding of the map $A \to \exp ( A t )$ and its sensitivity is helpful when assessing the performance of algorithms for computing the matrix exponential. Work in this direction includes:

B. K˚agstr¨om (1977). “Bounds and Perturbation Bounds for the Matrix Exponential,” BIT 17, 39–57.

C.F. Van Loan (1977). “The Sensitivity of the Matrix Exponential,” SIAM J. Numer. Anal. 14, 971–981.

R. Mathias (1992). “Evaluating the Fr´echet Derivative of the Matrix Exponential,” Numer. Math. 63, 213–226.

I. Najfeld and T.F. Havel (1995). “Derivatives of the Matrix Exponential and Their Computation,” Adv. Appl. Math. 16, 321–375.

A.H. Al-Mohy and N.J. Higham (2009). “Computing the Fr\`echet Derivative of the Matrix Exponential, with an Application to Condition Number Estimation,” SIAM J. Matrix Anal. Applic. 30, 1639– 1657.

A software package for computing small dense and large sparse matrix exponentials in Fortran and Matlab is presented in the following reference:

R.B. Sidje (1998) “Expokit: a Software Package for Computing Matrix Exponentials,” ACM Trans. Math. Softw. 24, 130–156.

Consideration of P9.3.2 and P9.3.5 shows that the exponential of a structured matrix can have important properties, see:

J. Xue and Q. Ye (2008). “Entrywise Relative Perturbation Bounds for Exponentials of Essentially Non-negative Matrices,” Numer. Math. 110, 393–403.

J. Cardoso and F.S. Leite (2010). “Exponentials of Skew-Symmetric Matrices and Logarithms of Orthogonal Matrices,” J. Comput. Appl. Math. 233, 2867–2875.

# 9.4 The Sign, Square Root, and Log of a Matrix

The matrix logarithm problem is the inverse of the matrix exponential problem. Not surprisingly, there is an inverse of the scaling and squaring procedure given in §9.3.1 that involves repeated matrix square roots. Thus, before we can discuss log(A) we need to understand the $\sqrt { A }$ problem. This in turn has connections to the matrix sign function and the polar decomposition.

# 9.4.1 The Matrix Sign Function

For all $z \in \mathbb { C }$ that are not on the imaginary axis, we define the sign(·) function by

$$
\operatorname{sign} (z) = \left\{ \begin{array}{l l} - 1 & \text {if} \operatorname{Re} (z) <   0, \\ + 1 & \text {if} \operatorname{Re} (z) > 0. \end{array} \right.
$$

The sign of a matrix has a particularly simple form Suppose $A \in \mathbb { C } ^ { n \times n }$ has no pure imaginary eigenvalues and that the blocks in its JCF $A = X J X ^ { - 1 }$ are ordered so that

$$
J = \left[ \begin{array}{c c} J _ {1} & 0 \\ 0 & J _ {2} \end{array} \right] \begin{array}{c} m _ {1} \\ m _ {2} \end{array}
$$

where the eigenvalues of $J _ { 1 } \in \mathbb { C } ^ { m _ { 1 } \times m _ { 1 } }$ lie in the open left half plane and the eigenvalues of $J _ { 2 } \in \mathbb { C } ^ { m _ { 2 } \times m _ { 2 } }$ lie in the open right half plane. Noting that all the derivatives of the sign function are zero, it follows from Theorem 9.1.1 that

$$
\mathrm{sign} (A) = X \left[ \begin{array}{c c} \mathrm{sign} (J _ {1}) & 0 \\ 0 & \mathrm{sign} (J _ {2}) \end{array} \right] X ^ {- 1} = X \left[ \begin{array}{c c} - I _ {m _ {1}} & 0 \\ 0 & I _ {m _ {2}} \end{array} \right] X ^ {- 1}.
$$

With the partitionings

$$
X = \left[ \begin{array}{c c} X _ {1} & X _ {2} \end{array} \right] \qquad \qquad X ^ {- H} = \left[ \begin{array}{c c} Y _ {1} & Y _ {2} \end{array} \right] \quad ,
$$

we have

$$
\operatorname{sign} (A) = X _ {2} Y _ {2} ^ {H} - X _ {1} Y _ {1} ^ {H}
$$

$$
I _ {n} = X _ {1} Y _ {1} ^ {H} + X _ {2} Y _ {2} ^ {H}
$$

and so

$$
X _ {2} Y _ {2} ^ {H} = \frac {1}{2} \left(I _ {n} + \operatorname{sign} (A)\right).
$$

Suppose apply QR-with-column pivoting to this rank- $\mathbf { \nabla } \cdot m _ { 2 }$ matrix:

$$
\frac {1}{2} \left(I _ {n} + \mathrm{sign} (A)\right) \Pi = Q R.
$$

It follows that ran $\left( Q ( : , 1 { : } m _ { 2 } ) \right) = \operatorname { r a n } ( X _ { 2 } )$ , the invariant subspace associated with A’s right half-plane eigenvalues. Thus, an approximation of sign(A) yields approximate invariant subspace information.

A number of iterative methods for computing sign(A) have been proposed. The fact that sign(z) is a zero of $g ( z ) = z ^ { 2 } - 1$ suggests a matrix analogue of the Newton iteration

$$
z _ {k + 1} = z _ {k} - \frac {g (z _ {k})}{g ^ {\prime} (z _ {k})} = \frac {1}{2} \left(z _ {k} + \frac {1}{z _ {k}}\right),
$$

i.e.,

$$
S _ {0} = A
$$

for $k = 0 , 1 , \ldots$ . (9.4.1)

$$
S _ {k + 1} = \left(S _ {k} + S _ {k} ^ {- 1}\right) / 2
$$

end

We proceed to show that this iteration is well-defined and converges to sign(A), assuming that A has no eigenvalues on the imaginary axis.

Note that if $a + b i$ is an eigenvalue of $S _ { k }$ , then

$$
\frac {1}{2} \left(a + b i + \frac {1}{a + b i}\right) = \frac {a}{2} \left(1 + \frac {1}{a ^ {2} + b ^ {2}}\right) + \frac {b}{2} \left(1 - \frac {1}{a ^ {2} + b ^ {2}}\right) i
$$

is an eigenvalue of $S _ { k + 1 }$ . Thus, if $S _ { k }$ is nonsingular, then $S _ { k + 1 }$ is nonsingular. It follows by induction that (9.4.1) is defined. Moreover, sign $( S _ { k } ) = \mathrm { s i g n } ( A )$ because an eigenvalue cannot “jump” across the imaginary axis during the iteration.

To prove that $S _ { k }$ converges to $S = \mathrm { s i g n } ( A )$ , we first observe that $S S _ { k } = S _ { k } S$ since both matrices are rational functions of A. Using this commutivity result and the identity $S ^ { 2 } = S$ , it is easy to show that

$$
S _ {k + 1} - S = \frac {1}{2} S _ {k} ^ {- 1} (S _ {k} - S) ^ {2} \tag {9.4.2}
$$

and

$$
S _ {k + 1} + S = \frac {1}{2} S _ {k} ^ {- 1} (S _ {k} + S) ^ {2}. \tag {9.4.3}
$$

If M is a matrix and sign(M ) is defined, then $M + \mathrm { s i g n } ( M )$ is nonsingular because its eigenvalues have the form $\lambda + \mathrm { s i g n } ( \lambda )$ which are clearly nonzero. Thus, the matrix

$$
S _ {k} + S = S _ {k} + \operatorname{sign} (A) = S _ {k} + \operatorname{sign} (S _ {k})
$$

is nonsingular. By manipulating equations (9.4.2) and (9.4.3) we conclude that if

$$
G _ {k} = (S _ {k} - S) (S _ {k} + S) ^ {- 1}, \tag {9.4.4}
$$

then $G _ { k + 1 } = G _ { k } ^ { 2 }$ . It follows by induction that $G _ { k } = G _ { 0 } ^ { 2 ^ { k } } . { \mathrm { ~ I f ~ } } \lambda \in \lambda ( A )$ , then

$$
\mu = \frac {\lambda - \mathrm{sign} (\lambda)}{\lambda + \mathrm{sign} (\lambda)}
$$

is an eigenvalue of $G _ { 0 } = ( A - S ) ( A + S ) ^ { - 1 }$ . Since $| \mu | < 1$ it follows from Lemma 7.3.2 that $G _ { k } \to 0$ and so

$$
S _ {k} = S (I + G _ {k}) (I - G _ {k}) ^ {- 1} \rightarrow S.
$$

Taking norms in (9.4.2) we conclude that the rate of convergence is quadratic:

$$
\| S _ {k + 1} - S \| \leq \frac {1}{2} \| S _ {k} ^ {- 1} \| \cdot \| S _ {k} - S \| ^ {2}.
$$

The overall efficiency of the method in practice is a concern since $O ( n ^ { 3 } )$ flops per iteration are required. To address this issue several enhancements of the basic iteration (9.4.1) have been proposed. One idea is to incorporate the Newton approximation

$$
S _ {k} ^ {- 1} \approx S _ {k} (2 I - S _ {k} ^ {2}).
$$

(See P9.4.1.) Using this estimate instead of the actual inverse in (9.4.1) gives update step

$$
S _ {k + 1} = \frac {1}{2} (S _ {k} + S _ {k} (2 I - S _ {k} ^ {2}) = \frac {1}{2} S _ {k} (3 I - S _ {k} ^ {2}). \tag {9.4.5}
$$

This is referred to as the Newton-Schultz iteration. Another idea is to introduce a scale factor:

$$
S _ {k + 1} = \frac {1}{2} \left((\mu_ {k} S _ {k}) + (\mu_ {k} S _ {k}) ^ {- 1}\right). \tag {9.4.6}
$$

Interesting choices for $\mu _ { k }$ include $\lvert \operatorname* { d e t } ( S _ { k } ) \rvert ^ { 1 / n } , \sqrt { \rho ( S _ { k } ^ { - 1 } ) / \rho ( S _ { k } ) }$ , and $\sqrt { \parallel S _ { k } ^ { - 1 } \parallel \parallel S _ { k } \parallel }$ where $\rho ( \cdot )$ is the spectral radius. For insights into the effective computation of the matrix sign function and related stability issues, see Kenney and Laub (1991, 1992), Higham (2007), and Higham (FOM, Chap. 5).

# 9.4.2 The Matrix Square Root

Ambiguity arises in the $f ( A )$ problem if the underlying function has branches. For example, if $f ( x ) = { \sqrt { x } }$ and

$$
A = \left[ \begin{array}{c c} 4 & 1 0 \\ 0 & 9 \end{array} \right],
$$

then

$$
A = \left[ \begin{array}{c c} 2 & 2 \\ 0 & 3 \end{array} \right] ^ {2} = \left[ \begin{array}{c c} - 2 & 1 0 \\ 0 & 3 \end{array} \right] ^ {2} = \left[ \begin{array}{c c} - 2 & - 2 \\ 0 & - 3 \end{array} \right] ^ {2} = \left[ \begin{array}{c c} 2 & - 1 0 \\ 0 & - 3 \end{array} \right] ^ {2},
$$

which shows that there are at least four legitimate choices for ${ \sqrt { A } } .$ . To clarify the situation we say F is the principal square root of A if (a) $F ^ { 2 } = A$ and (b) the eigenvalues of $F$ have positive real part. We designate this matrix by $A ^ { 1 / 2 }$ .

Analogous to the Newton iteration for scalar square roots, $x _ { k + 1 } = ( x _ { k } + a / x _ { k } ) / 2$ , we have

$$
X _ {0} = A
$$

for $k = 0 , 1 , \ldots$ . (9.4.7)

$$
X _ {k + 1} = \left(X _ {k} + X _ {k} ^ {- 1} A\right) / 2
$$

end

Notice the similarity between this iteration and the Newton sign iteration (9.4.1). Indeed, by making the substitution $X _ { k } = A ^ { 1 / 2 } S _ { k }$ in (9.4.7) we obtain the Newton sign iteration for $A ^ { 1 / 2 }$ . Global convergence and local quadratic convergence follow from what we know about (9.4.1).

Another connection between the matrix sign problem and the matrix square root problem is revealed by applying the Newton sign iteration to the matrix

$$
\tilde {A} = \left[ \begin{array}{l l} 0 & A \\ I & 0 \end{array} \right].
$$

Designate the iterates by $\tilde { S } _ { k }$ . We show by induction that $\tilde { S } _ { k }$ has the form

$$
\tilde {S} _ {k} = \left[ \begin{array}{c c} 0 & X _ {k} \\ Y _ {k} & 0 \end{array} \right].
$$

This is true for $k = 0$ by setting $X _ { 0 } = A$ and $Y _ { 0 } = I$ . To see that the result holds for $k > 0$ , observe that

$$
\tilde {S} _ {k + 1} = \frac {1}{2} \left(\tilde {S} _ {k} + \tilde {S} _ {k} ^ {- 1}\right) = \frac {1}{2} \left(\left[ \begin{array}{c c} 0 & X _ {k} \\ Y _ {k} & 0 \end{array} \right] + \left[ \begin{array}{c c} 0 & Y _ {k} ^ {- 1} \\ X _ {k} ^ {- 1} & 0 \end{array} \right]\right)
$$

and thus

$$
X _ {k + 1} = \left(X _ {k} + Y _ {k} ^ {- 1}\right) / 2, \quad Y _ {k + 1} = \left(Y _ {k} + X _ {k} ^ {- 1}\right) / 2. \tag {9.4.8}
$$

Another induction argument shows that

$$
X _ {k} = A Y _ {k}, \quad k = 0, 1, \dots , \tag {9.4.9}
$$

and so

$$
X _ {k + 1} = \left(X _ {k} + A X _ {k} ^ {- 1}\right) / 2, \quad Y _ {k + 1} = \left(Y _ {k} + A ^ {- 1} Y _ {k} ^ {- 1}\right) / 2. \tag {9.4.10}
$$

It follows that $X _ { k }  A ^ { 1 / 2 }$ and $Y _ { k }  A ^ { - 1 / 2 }$ and we have established the following identity:

$$
\mathrm{sign} \left(\left[ \begin{array}{c c} 0 & A \\ I & 0 \end{array} \right]\right) = \left[ \begin{array}{c c} 0 & A ^ {1 / 2} \\ A ^ {- 1 / 2} & 0 \end{array} \right].
$$

Equation (9.4.8) defines the Denman-Beavers iteration which turns out to have better numerical properties than (9.4.7). See Meini (2004), Higham (FOM, Chap. 6), and Higham (2008) for an analysis of these and other matrix square root algorithms.

# 9.4.3 The Polar Decomposition

If $z = a + b i \in \mathbb { C }$ is a nonzero complex number, then its polar representation is a factorization of the form $z = e ^ { i \theta } r$ where $r = \sqrt { a ^ { 2 } + b ^ { 2 } }$ and $e ^ { i \theta } = \cos ( \theta ) + i \sin ( \theta )$ i s defined by $( \cos ( \theta ) , \sin ( \theta ) ) = ( a / r , b / r )$ . The polar decomposition of a matrix is similar.

Theorem 9.4.1 (Polar Decomposition). If $A \in \mathbb { R } ^ { m \times n }$ and m $\geq n$ , then there exists a matrix $U \in \mathbb { R } ^ { m \times n }$ with orthonormal columns and a symmetric positive semidefinite $P \in \mathbb { R } ^ { n \times n }$ so that $A = U P$ .

Proof. Suppose $\ U _ { A } ^ { T } A V _ { A } \ = \ \Sigma _ { A }$ is the thin SVD of A. It is easy to show that if $\begin{array} { r c l } { U } & { = } & { U _ { A } V _ { A } ^ { \bar { T } } } \end{array}$ and $\begin{array} { l l l } { P } & { = } & { { V _ { A } } { \Sigma _ { A } } { V _ { A } ^ { T } } } \end{array}$ , then $A = U P$ and U and P have the required properties.

We refer to U as the orthogonal polar factor and P as the symmetric polar factor. Note that $P = ( A ^ { T } A ) ^ { 1 / 2 }$ and if $\mathsf { r a n k } ( A ) = n$ , then $U = A ( A ^ { T } A ) ^ { - 1 / 2 }$ . An important application of the polar decomposition is the orthogonal Procrustes problem (see §6.4.1).

Various iterative methods for computing the orthogonal polar factor have been proposed. A quadratically convergent Newton iteration for the square nonsingular case proceeds by repeatedly averaging the current iterate with the inverse of its transpose:

$$
X _ {0} = A \quad \left(\text { Assume } A \in \mathbb {R} ^ {n \times n} \text { is   nonsingular }\right)
$$

$$
\text { for } k = 0, 1, \dots \tag {9.4.11}
$$

$$
X _ {k + 1} = \left(X _ {k} + X _ {k} ^ {- T}\right) / 2
$$

end

To show that this iteration is well defined we assume that for some k the matrix $X _ { k }$ is nonsingular and that $X _ { k } = U _ { k } P _ { k }$ is its polar decomposition. It follows that

$$
X _ {k + 1} = \frac {1}{2} \left(X _ {k} + X _ {k} ^ {- T}\right) = \frac {1}{2} \left(U _ {k} P _ {k} + U _ {k} P _ {k} ^ {- 1}\right) = U _ {k} \left(\frac {P _ {k} + P _ {k} ^ {- 1}}{2}\right). \tag {9.4.12}
$$

Since the average of a positive definite matrix and its inverse is also positive definite it follows that $X _ { k + 1 }$ is nonsingular. This shows by induction that (9.4.11) is well-defined and that the $P _ { k }$ satisfy

$$
P _ {k + 1} = (P _ {k} + P _ {k} ^ {- 1}) / 2, \qquad P _ {0} = P.
$$

This is precisely the Newton sign iteration (9.4.1) with starting matrix $P _ { 0 } = P$ . Since

$$
\left\| X _ {k} - U \right\| _ {2} = \left\| U (P _ {k} - I) \right\| _ {2} = \left\| P _ {k} - I \right\| _ {2}
$$

and $P _ { k }  \mathrm { s i g n } ( P ) = I$ quadratically, we conclude that $X _ { k }$ matrices in (9.4.11) converge to U quadratically.

Extensions to the rectangular case and various ways to accelerate (9.4.11) are discussed in Higham (1986), Higham and Schreiber (1990), Gander (1990), and Kenney and Laub (1992). In this regard the matrix sign function is (once again) a handy tool for deriving algorithms. Note that if $A = U _ { A } \Sigma _ { A } V _ { A } ^ { T }$ is the SVD of $A \in \bar { \mathbb { R } } ^ { n \times n }$ and

$$
Q = \frac {1}{\sqrt {2}} \left[ \begin{array}{c c} U _ {A} & 0 \\ 0 & V _ {A} \end{array} \right] \left[ \begin{array}{c c} I _ {n} & I _ {n} \\ I _ {n} & - I _ {n} \end{array} \right]
$$

then Q is orthogonal and

$$
Q ^ {T} \left[ \begin{array}{c c} 0 & A \\ A ^ {T} & 0 \end{array} \right] Q = \left[ \begin{array}{c c} \Sigma_ {A} & 0 \\ 0 & - \Sigma_ {A} \end{array} \right].
$$

It follows that

$$
\mathrm{sign} \left(\left[ \begin{array}{c c} 0 & A \\ A ^ {T} & 0 \end{array} \right]\right) = Q \left[ \begin{array}{c c} I _ {n} & 0 \\ 0 & - I _ {n} \end{array} \right] Q ^ {T} = \left[ \begin{array}{c c} 0 & U \\ U ^ {T} & 0 \end{array} \right]
$$

where $U = U _ { A } V _ { A } ^ { T }$ is the orthogonal polar factor of A.

There is a well-developed perturbation theory for the polar decomposition. A sample result for square nonsingular matrices due to Li and Sun (2003) says that the orthogonal polar factors U and $\tilde { U }$ for nonsingular $A , \tilde { A } \in \mathbb { R } ^ { n \times n }$ satisfy the bound

$$
\| U - \tilde {U} \| _ {F} \leq \frac {4 \| A - \tilde {A} \| _ {F}}{\sigma_ {n - 1} (A) + \sigma_ {n} (A) + \sigma_ {n - 1} (\tilde {A}) + \sigma_ {n} (\tilde {A})}.
$$

# 9.4.4 The Matrix Logarithm

Given $A \in \mathbb { R } ^ { n \times n }$ , a solution to the matrix equation $e ^ { X } = A$ is a logarithm of A. Note that if $X = \log ( A )$ , then $X + 2 k \pi i$ is also a logarithm. To remove this ambiguity we define the principal logarithm as follows. If the real eigenvalues of $A \in \mathbb { R } ^ { n \times n }$ are all positive then there is a unique real matrix X that satisfies $e ^ { X } = A$ with the property that its eigenvalues satisfy $\lambda ( X ) \subset \{ z \in \mathbb { C } : - \pi < \mathsf { I m } ( z ) < \pi \}$ .

Of course, the eigenvalue-based methods of §9.2 are applicable for the log(A) problem. We discuss an approximation method that is analogous to Algorithm 9.3.1, the scaling and squaring method for the matrix exponential

As with the exponential, there are a number of different series expansions for the log function that are of computational interest. The simplest is the Maclaurin expansion:

$$
\log (A) \approx M _ {q} (A) = \sum_ {k = 1} ^ {q} (- 1) ^ {k + 1} \frac {(A - I) ^ {k}}{k}.
$$

To apply this formula we must have $\rho ( A - I ) < 1$ where $\rho ( \cdot )$ is the spectral radius.

The Gregory series expansion for log(x) yields a rational approximation:

$$
\log (A) \approx G _ {q} (A) = - 2 \sum_ {k = 0} ^ {q} \frac {1}{2 k + 1} \left((I - A) (I + A) ^ {- 1}\right) ^ {2 k + 1}.
$$

For this to converge, the real parts of A’s eigenvalues must be positive.

Diagonal Pad´e approximants are also of interest. For example, the (3,3) Pad´e approximant is given by

$$
\log (A) \approx r _ {3 3} (A) = D (A) ^ {- 1} N (A)
$$

where

$$
D (A) = 6 0 I + 9 0 (A - I) + 3 6 (A - I) ^ {2} + 3 (A - I) ^ {3},
$$

$$
N (A) = 6 0 (A - I) + 6 0 (A - I) ^ {2} + 1 1 (A - I) ^ {3}.
$$

For an approximation of this type to be effective, the matrix A must be sufficiently close to the identity matrix. Repeated square roots are one way to achieve this:

$$
k = 0
$$

$$
A _ {0} = A
$$

while $\| A - I \| >$ tol

$$
k = k + 1
$$

$$
A _ {k} = A _ {k - 1} ^ {1 / 2}
$$

end

The Denman-Beavers iteration (9.4.8) can be invoked to compute the matrix square roots. If we next compute $F \approx \log ( A _ { k } )$ by using (say) an appropriately chosen Pade approximant, then log $( A ) \ = \ 2 ^ { k } \log ( A _ { k } ) \ \approx \ 2 ^ { k } F$ . This solution framework is referred to as inverse scaling and squaring. There are many details associated with the proper implementation of this procedure and we refer the reader to Cheng, Higham. Kenney, and Laub (2001), Higham (2001), and Higham (FOM, Chap. 11).

# Problems

P9.4.1 What does the Newton iteration look like when it is applied to find a root of the function $f ( x ) = 1 / x - a ?$ Develop an inverse-free Newton iteration for solving the matrix equation $X ^ { - 1 } - A$ .

P9.4.2 Show that if $\mu _ { k } > 0$ in (9.4.6), then $\mathrm { s i g n } ( S _ { k + 1 } ) = \mathrm { s i g n } ( S _ { k } )$ .

P9.4.3 Show that s $\operatorname { i g n } ( A ) = A ( A ^ { 2 } ) ^ { - 1 / 2 }$ .

P9.4.4 Verify Equation (9.4.9).

P9.4.5 In the Denman-Beavers iteration (9.4.8), define $M _ { k } = X _ { k } Y _ { k }$ and develop a recipe for $M _ { k + 1 }$

P9.4.6 Show that if we apply the Newton square root iteration (9.4.9) to a symmetric positive definite matrix A, then $A _ { k } - A _ { k + 1 }$ is positive definite for all k.

P9.4.7 Suppose A is normal. Relate the polar factors of $e ^ { A }$ to $S = ( A - A ^ { T } ) / 2$ and $T = ( A + A ^ { T } ) / 2$

P9.4.8 Show that the polar decomposition of a nonsingular matrix is unique. Hint: If $A = U _ { 1 } P _ { 1 }$ and $A = U _ { 2 } P _ { 2 }$ are two polar decompositions, then $U _ { 2 } ^ { T } U _ { 1 } = P _ { 2 } P _ { 1 } ^ { - 1 }$ and $U _ { 1 } ^ { T } \bar { U _ { 2 } } = P _ { 1 } P _ { 2 } ^ { - 1 }$ have the same eigenvalues.

P9.4.9 Give a closed-form expression for the polar decomposition $A = U P$ of a real 2-by-2 matrix. Under what conditions is U a rotation?

P9.4.10 Give a closed-form expression for log(Q) where Q is a 2-by-2 rotation matrix.

P9.4.11 Formulate an $m < n$ version of the polar decomposition for $A \in \mathbb { R } ^ { m \times n }$ .

P9.4.12 Let A by an n-by-n symmetric positive definite matrix. (a) Show that there exists a unique symmetric positive definite X such that $\dot { \boldsymbol { A } } = \boldsymbol { X } ^ { 2 }$ . (b) Show that if $X _ { 0 } = I$ and

$$
X _ {k + 1} = (X _ {k} + A X _ {k} ^ {- 1}) / 2
$$

then $X _ { k } \to { \sqrt { A } }$ quadratically where $\sqrt { A }$ denotes the matrix X in part (a).

P9.4.13 Show that

$$
X (t) = C _ {1} \cos (t \sqrt {A}) + C _ {2} \sqrt {A ^ {- 1}} \sin (t \sqrt {A})
$$

solves the initial value problem $\ddot { X } ( t ) = - A X ( t ) , X ( 0 ) = C _ { 1 } , \dot { X } ( 0 ) = C _ { 2 }$ . Assume that A is symmetric positive definite.

# Notes and References for §9.4

Everything in this section is covered in greater depth in Higham (FOM). See also:

N.J. Higham (2005). “Functions of Matrices,” in Handbook of Linear Algebra, L. Hogben (ed.), Chapman and Hall, Boca Raton, ${ \mathrm { F L } } , \{ 1 1 - 1 - \ S 1 1 - 1 3 . $ .

Papers that discuss the ubiquitous matrix sign function and its applications include:

R. Byers (1987). “Solving the Algebraic Riccati Equation with the Matrix Sign Function,” Linear Alg. Applic. 85, 267–279.   
C.S. Kenney and A.J. Laub (1991). “Rational Iterative Methods for the Matrix Sign Function,” SIAM J. Matrix Anal. Appl. 12, 273–291.   
C.S. Kenney, A.J. Laub, and P.M. Papadopouos (1992). “Matrix Sign Algorithms for Riccati Equations,” IMA J. Math. Control Info. 9, 331–344.   
C.S. Kenney and A.J. Laub (1992). “On Scaling Newton’s Method for Polar Decomposition and the Matrix Sign Function,” SIAM J. Matrix Anal. Applic. 13, 688–706.   
R. Byers, C. He, and V. Mehrmann (1997). “The Matrix Sign Function Method and the Computation of Invariant Subspaces,” SIAM J. Matrix Anal. Applic. 18, 615–632.   
Z. Bai and J.W. Demmel (1998). “Using the Matrix Sign Function to Compute Invariant Subspaces,” SIAM J. Matrix Anal. Applic. 19, 2205–2225.   
N.J. Higham (1994). “The Matrix Sign Decomposition and Its Relation to the Polar Decomposition,” Lin. Alg. Applic. 212/213, 3–20.   
N.J. Higham, D.S. Mackey, N. Mackey, and F. Tisseur (2004). “Computing the Polar Decomposition and the Matrix Sign Decomposition in Matrix Groups,” SIAM J. Matrix Anal. Applic. 25, 1178– 1192.

Various aspects of the matrix square root problem are discussed in:

E.D. Denman and A.N. Beavers (1976). “The Matrix Sign Function and Computations in Systems,” Appl. Math. Comput., 2, 63–94.

˚A. Bj¨orck and S. Hammarling (1983). “A Schur Method for the Square Root of a Matrix,” Lin. Alg. Applic. 52/53, 127–140.

N.J. Higham (1986). “Newton’s Method for the Matrix Square Root,” Math. Comput. 46, 537–550.

N.J. Higham (1987). “Computing Real Square Roots of a Real Matrix,” Lin. Alg. Applic. 88/89, 405–430.

N.J. Higham (1997). “Stable Iterations for the Matrix Square Root,” Numer. Algorithms 15, 227–242.

Y.Y. Lu (1998). “A Pad´e Approximation Method for Square Roots of Symmetric Positive Definite Matrices,” SIAM J. Matrix Anal. Applic. 19, 833–845.

N.J. Higham, D.S. Mackey, N. Mackey, and F. Tisseur (2005). “Functions Preserving Matrix Groups and Iterations for the Matrix Square Root,” SIAM J. Matrix Anal. Applic. 26, 849–877.

C.-H. Guo and N. J. Higham (2006). “A Schur–Newton Method for the Matrix pth Root and its Inverse,” SIAM J. Matrix Anal. Applic. 28, 788–804.

B. Meini (2004). “The Matrix Square Root from a New Functional Perspective: Theoretical Results and Computational Issues,” SIAM J. Matrix Anal. Applic. 26, 362–376.

A. Frommer and B. Hashemi (2009). “Verified Computation of Square Roots of a Matrix,” SIAM J. Matrix Anal. Applic. 31, 1279–1302.   
Computational aspects of the polar decomposition and its generalizations are covered in:   
N.J. Higham (1986). “Computing the Polar Decomposition with Applications,” SIAM J. Sci. Statist. Comp. 7, 1160–1174.   
R.S. Schreiber and B.N. Parlett (1988). “Block Reflectors: Theory and Computation,” SIAM J. Numer. Anal. 25, 189–205.   
N.J. Higham and R.S. Schreiber (1990). “Fast Polar Decomposition of an Arbitrary Matrix,” SIAM J. Sci. Statist. Comput. 11, 648–655.   
N.J. Higham and P. Papadimitriou (1994). “A Parallel Algorithm for Computing the Polar Decomposition,” Parallel Comput. 20, 1161–1173.   
A.A. Dubrulle (1999). “An Optimum Iteration for the Matrix Polar Decomposition,” ETNA 8, 21–25.   
A. Zanna and H. Z. Munthe-Kaas (2002). “Generalized Polar Decompositions for the Approximation of the Matrix Exponential,” SIAM J. Matrix Anal. Applic. 23, 840–862.   
B. Laszkiewicz and K. Zietak (2006). “Approximation of Matrices and a Family of Gander Methods for Polar Decomposition,” BIT 46, 345–366.   
R. Byers and H. Xu (2008). “A New Scaling for Newton’s Iteration for the Polar Decomposition and Its Backward Stability,” SIAM J. Matrix Anal. Applic. 30, 822–843.   
N.J. Higham, C. Mehl, and F. Tisseur (2010). “The Canonical Generalized Polar Decomposition,” SIAM J. Matrix Anal. Applic. 31, 2163–2180.   
For an analysis as to whether or not the polar decomposition can be computed in a finite number of steps, see:   
A. George and Kh. Ikramov (1996). “Is The Polar Decomposition Finitely Computable?,” SIAM J. Matrix Anal. Applic. 17, 348–354.   
A. George and Kh. Ikramov (1997). “Addendum: Is The Polar Decomposition Finitely Computable?,” SIAM J. Matrix Anal. Appl. 18, 264–264.   
There is a considerable literature concerned with how the polar factors change under perturbation:   
R. Mathias (1993). “Perturbation Bounds for the Polar Decomposition,” SIAM J. Matrix Anal. Applic. 14, 588–597.   
R.-C. Li (1997). “Relative Perturbation Bounds for the Unitary Polar Factor,” BIT 37, 67–75.   
F. Chaitin-Chatelin, S. Gratton (2000). “On the Condition Numbers Associated with the Polar Factorization of a Matrix,” Numer. Lin. Alg. 7, 337–354.   
W. Li and W. Sun (2003). “New Perturbation Bounds for Unitary Polar Factors,” SIAM J. Matrix Anal. Applic. 25, 362–372.   
Finally, details concerning the matrix logarithm and its computation may be found in:   
B.W. Helton (1968). “Logarithms of Matrices,” Proc. AMS 19, 733–736.   
L. Dieci (1996). “Considerations on Computing Real Logarithms of Matrices, Hamiltonian Logarithms, and Skew-Symmetric Logarithms,” Lin. Alg. Applic. 244, 35–54.   
L. Dieci, B. Morini, and A. Papini (1996). “Computational Techniques for Real Logarithms of Matrices,” SIAM J. Matrix Anal. Applic. 17, 570–593.   
C. S. Kenney and A. J. Laub (1998). “A Schur-Fr´echet Algorithm for Computing the Logarithm and Exponential of a Matrix,” SIAM J. Matrix Anal. Applic. 19, 640–663.   
L. Dieci (1998). “Real Hamiltonian Logarithm of a Symplectic Matrix,” Lin. Alg. Applic. 281, 227–246.   
L. Dieci and A. Papini (2000). “Conditioning and Pad´e Approximation of the Logarithm of a Matrix,” SIAM J. Matrix Anal. Applic. 21, 913–930.   
N.J. Higham (2001). “Evaluating Pad´e Approximants of the Matrix Logarithm,” SIAM J. Matrix Anal. Applic. 22, 1126–1135.   
S.H. Cheng, N.J. Higham, C.S. Kenney, and A.J. Laub (2001). “Approximating the Logarithm of a Matrix to Specified Accuracy,” SIAM J. Matrix Anal. Applic. 22, 1112–1125.
