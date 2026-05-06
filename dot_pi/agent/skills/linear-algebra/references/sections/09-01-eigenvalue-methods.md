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
