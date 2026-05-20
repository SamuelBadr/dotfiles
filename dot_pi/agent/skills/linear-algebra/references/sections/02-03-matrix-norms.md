# 2.3 Matrix Norms

The analysis of matrix algorithms requires use of matrix norms. For example, the quality of a linear system solution may be poor if the matrix of coefficients is “nearly singular.” To quantify the notion of near-singularity, we need a measure of distance on the space of matrices. Matrix norms can be used to provide that measure.

# 2.3.1 Definitions

Since $\mathbb { R } ^ { m \times n }$ is isomorphic to $\mathbb { R } ^ { m n }$ , the definition of a matrix norm should be equivalent to the definition of a vector norm. In particular, $f { : \mathbb { R } ^ { m \times n } } \to \mathbb { R }$ is a matrix norm if the following three properties hold:

$$
f (A) \geq 0, \quad A \in \mathbb {R} ^ {m \times n}, \quad (f (A) = 0 \text {   iff   } A = 0)
$$

$$
f (A + B) \leq f (A) + f (B), \quad A, B \in \mathbb {R} ^ {m \times n},
$$

$$
f (\alpha A) = | \alpha | f (A), \quad \alpha \in \mathbb {R}, A \in \mathbb {R} ^ {m \times n}.
$$

As with vector norms, we use a double bar notation with subscripts to designate matrix norms, ${ \mathrm { i . e . , ~ } } \| \ A \| = f ( A )$ .

The most frequently used matrix norms in numerical linear algebra are the Frobenius norm

$$
\| A \| _ {F} = \sqrt {\sum_ {i = 1} ^ {m} \sum_ {j = 1} ^ {n} \left| a _ {i j} \right| ^ {2}} \tag {2.3.1}
$$

and the p-norms

$$
\| A \| _ {p} = \sup _ {x \neq 0} \frac {\| A x \| _ {p}}{\| x \| _ {p}}. \tag {2.3.2}
$$

Note that the matrix p-norms are defined in terms of the vector p-norms discussed in the previous section. The verification that (2.3.1) and (2.3.2) are matrix norms is left as an exercise. It is clear that $\| A \| _ { p }$ is the p-norm of the largest vector obtained by applying A to a unit p-norm vector:

$$
\| A \| _ {p} = \sup _ {x \neq 0} \left\| A \left(\frac {x}{\| x \| _ {p}}\right) \right\| _ {p} = \max _ {\| x \| _ {p} = 1} \| A x \| _ {p}.
$$

It is important to understand that (2.3.2) defines a family of norms—the 2-norm on $\mathbb { R } ^ { 3 \times 2 }$ is a different function from the 2-norm on $\mathbb { R } ^ { 5 \times 6 }$ . Thus, the easily verified inequality

$$
\| A B \| _ {p} \leq \| A \| _ {p} \| B \| _ {p}, \quad A \in \mathbb {R} ^ {m \times n}, B \in \mathbb {R} ^ {n \times q} \tag {2.3.3}
$$

is really an observation about the relationship between three different norms. Formally, we say that norms $f _ { 1 } , \ f _ { 2 }$ , and $f _ { 3 }$ on $\mathbb { R } ^ { m \times q }$ , $\mathbb { R } ^ { m \times n }$ , and $\mathbb { R } ^ { n \times q }$ are mutually consistent if for all matrices $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { n \times q }$ we have $f _ { 1 } ( A B ) \le f _ { 2 } ( A ) f _ { 3 } ( B )$ , or, in subscript-free norm notation:

$$
\| A B \| \leq \| A \| \| B \|. \tag {2.3.4}
$$

Not all matrix norms satisfy this property. For example, if $\left\| \ A \right\| _ { \Delta } = \operatorname* { m a x } \left| a _ { i j } \right|$ and

$$
A = B = \left[ \begin{array}{l l} 1 & 1 \\ 1 & 1 \end{array} \right],
$$

then $\left\| \left. A B \right\| _ { \Delta } > \right\| \left. A \right\| _ { \Delta } \right\| B \left\| _ { \Delta }$ . For the most part, we work with norms that satisfy (2.3.4).

The p-norms have the important property that for every $A \in \mathbb { R } ^ { m \times n }$ and $\boldsymbol { x } \in \mathbb { R } ^ { n }$ we have

$$
\| A x \| _ {p} \leq \| A \| _ {p} \| x \| _ {p}.
$$

More generally, for any vector norm $\| \cdot \| _ { \alpha }$ on $\mathbb { R } ^ { n }$ and $\| \cdot \| _ { \beta }$ on $\mathbb { R } ^ { m }$ we have $\Vert { \mathbf { \Omega } } A x \Vert _ { \beta } \ \leq$ $\| A \| _ { \alpha , \beta } \| \boldsymbol { x } \| _ { \alpha }$ where $\| A \| _ { \alpha , \beta }$ is a matrix norm defined by

$$
\| A \| _ {\alpha , \beta} = \sup _ {x \neq 0} \frac {\| A x \| _ {\beta}}{\| x \| _ {\alpha}}. \tag {2.3.5}
$$

We say that $\| \cdot \| _ { \alpha , \beta }$ is subordinate to the vector norms $\| \cdot \| _ { \alpha }$ and $\| \cdot \| _ { \beta }$ . Since the set $\{ x \in \mathbb { R } ^ { n } : \| x \| _ { \alpha } ^ { - \gamma } = 1 \}$ is compact and $\| \cdot \| _ { \beta }$ is continuous, it follows that

$$
\| A \| _ {\alpha , \beta} = \max _ {\| x \| _ {\alpha} = 1} \| A x \| _ {\beta} = \| A x _ {*} \| _ {\beta} \tag {2.3.6}
$$

for some $\boldsymbol { x } _ { * } \in \mathbb { R } ^ { n }$ having unit α-norm.

# 2.3.2 Some Matrix Norm Properties

The Frobenius and p-norms (especially $p = 1 , 2 , \infty )$ satisfy certain inequalities that are frequently used in the analysis of a matrix computation. If $A \in \mathbb { R } ^ { m \times n }$ we have

$$
\| A \| _ {2} \leq \| A \| _ {F} \leq \sqrt {\min \{m , n \}} \| A \| _ {2}, \tag {2.3.7}
$$

$$
\max _ {i, j} \left| a _ {i j} \right| \leq \| A \| _ {2} \leq \sqrt {m n} \max _ {i, j} \left| a _ {i j} \right|, \tag {2.3.8}
$$

$$
\| A \| _ {1} = \max _ {1 \leq j \leq n} \sum_ {i = 1} ^ {m} | a _ {i j} |, \tag {2.3.9}
$$

$$
\| A \| _ {\infty} = \max _ {1 \leq i \leq m} \sum_ {j = 1} ^ {n} | a _ {i j} |, \tag {2.3.10}
$$

$$
\frac {1}{\sqrt {n}} \parallel A \parallel_ {\infty} \leq \parallel A \parallel_ {2} \leq \sqrt {m} \parallel A \parallel_ {\infty}, \tag {2.3.11}
$$

$$
\frac {1}{\sqrt {m}} \parallel A \parallel_ {1} \leq \parallel A \parallel_ {2} \leq \sqrt {n} \parallel A \parallel_ {1}. \tag {2.3.12}
$$

If $A \in \mathbb { R } ^ { m \times n } , 1 \leq i _ { 1 } \leq i _ { 2 } \leq m$ , and $1 \leq j _ { 1 } \leq j _ { 2 } \leq n$ , then

$$
\left\| A \left(i _ {1}: i _ {2}, j _ {1}: j _ {2}\right) \right\| _ {p} \leq \left\| A \right\| _ {p}. \tag {2.3.13}
$$

The proofs of these relationships are left as exercises. We mention that a sequence $\{ A ^ { ( k ) } \} \in \mathbb { R } ^ { m \times n }$ converges if there exists a matrix $A \in \mathbb { R } ^ { m \times n }$ such that

$$
\lim _ {k \to \infty} \| A ^ {(k)} - A \| = 0.
$$

The choice of norm is immaterial since all norms on $\mathbb { R } ^ { m \times n }$ are equivalent.

# 2.3.3 The Matrix 2-Norm

A nice feature of the matrix 1-norm and the matrix ∞-norm is that they are easy, $O ( n ^ { 2 } )$ computations. (See (2.3.9) and (2.3.10).) The calculation of the 2-norm is considerably more complicated.

Theorem 2.3.1. If $A \in \mathbb { R } ^ { m \times n }$ , then there exists a unit 2-norm n-vector z such that $A ^ { T } A z = \mu ^ { 2 } z$ where $\mu = \parallel A \parallel _ { 2 }$ .

Proof. Suppose $z \in \mathbb { R } ^ { n }$ is a unit vector such that $\parallel A z \parallel _ { 2 } = \parallel A \parallel _ { 2 }$ . Since z maximizes the function

$$
g (x) = \frac {1}{2} \frac {\parallel A x \parallel_ {2} ^ {2}}{\parallel x \parallel_ {2} ^ {2}} = \frac {1}{2} \frac {x ^ {T} A ^ {T} A x}{x ^ {T} x}
$$

it follows that it satisfies $\nabla g ( z ) = 0$ where $\nabla g$ is the gradient of g. A tedious differentiation shows that for $i = 1 { : } n$

$$
\frac {\partial g (z)}{\partial z _ {i}} = \left[ (z ^ {T} z) \sum_ {j = 1} ^ {n} (A ^ {T} A) _ {i j} z _ {j} - (z ^ {T} A ^ {T} A z) z _ {i} \right] \bigg / (z ^ {T} z) ^ {2}.
$$

In vector notation this says that $A ^ { T } A z = ( z ^ { T } A ^ { T } A z ) z$ . The theorem follows by setting $\mu = \parallel A z \parallel _ { 2 }$ .

The theorem implies that $\parallel A \parallel _ { 2 } ^ { 2 }$ is a zero of $p ( \lambda ) = \mathsf { d e t } ( A ^ { T } A - \lambda I )$ . In particular,

$$
\parallel A \parallel_ {2} = \sqrt {\lambda_ {\max} (A ^ {T} A)}
$$

We have much more to say about eigenvalues in Chapters 7 and 8. For now, we merely observe that 2-norm computation is iterative and a more involved calculation than those of the matrix 1-norm or ∞-norm. Fortunately, if the object is to obtain an order-of-magnitude estimate of $\parallel A \parallel _ { 2 }$ , then (2.3.7), (2.3.8), (2.3.11), or (2.3.12) can be used.

As another example of norm analysis, here is a handy result for 2-norm estimation.

Corollary 2.3.2. If $A \in \mathbb { R } ^ { m \times n }$ , then $\| A \| _ { 2 } \leq \sqrt { \| A \| _ { 1 } \| A \| _ { \infty } }$

Proof. $\mathrm { ~ I f ~ } z ~ \ne ~ 0$ is such that $A ^ { T } A z \ = \ \mu ^ { 2 } z$ with $\mu ~ = ~ \parallel A \parallel _ { 2 }$ , then $\textstyle \mu ^ { 2 } \parallel { \boldsymbol { z } } \parallel _ { 1 } =$ $\begin{array} { r } { \| A ^ { T } A z \| _ { 1 } \leq \| A ^ { T } \| _ { 1 } \| A \| _ { 1 } \| z \| _ { 1 } = \| A \| _ { \infty } \| A \| _ { 1 } \| z \| _ { 1 } . \quad \mathbb { D } } \end{array}$

# 2.3.4 Perturbations and the Inverse

We frequently use norms to quantify the effect of perturbations or to prove that a sequence of matrices converges to a specified limit. As an illustration of these norm applications, let us quantify the change in $A ^ { - 1 }$ as a function of change in A.

Lemma 2.3.3. If $F \in \mathbb { R } ^ { n \times n }$ and $\| \boldsymbol { F } \| _ { p } < 1$ , then $I - F$ is nonsingular and

$$
(I - F) ^ {- 1} = \sum_ {k = 0} ^ {\infty} F ^ {k}
$$

with

$$
\| (I - F) ^ {- 1} \| _ {p} \leq \frac {1}{1 - \| F \| _ {p}}.
$$

Proof. Suppose $I - F$ is singular. It follows that $( I - F ) x = 0$ for some nonzero x. But then $\| x \| _ { p } = \| F x \| _ { p }$ implies  $F \parallel _ { p } \geq 1$ , a contradiction. Thus, $I - F$ is nonsingular. To obtain an expression for its inverse consider the identity

$$
\left(\sum_ {k = 0} ^ {N} F ^ {k}\right) (I - F) = I - F ^ {N + 1}.
$$

Since $\| \ F \| _ { p } < 1$ it follows that $\operatorname* { l i m } _ { k \to \infty } F ^ { k } = 0$ because $\begin{array} { r } { \| \boldsymbol { F } ^ { k } \| _ { p } \leq \| \boldsymbol { F } \| _ { p } ^ { k } } \end{array}$ . Thus,

$$
\left(\lim _ {N \rightarrow \infty} \sum_ {k = 0} ^ {N} F ^ {k}\right) (I - F) = I.
$$

It follows that $( I - F ) ^ { - 1 } = \operatorname* { l i m } _ { N \to \infty } \sum _ { k = 0 } ^ { N } F ^ { k }$ . From this it is easy to show that

$$
\| (I - F) ^ {- 1} \| _ {p} \leq \sum_ {k = 0} ^ {\infty} \| F \| _ {p} ^ {k} = \frac {1}{1 - \| F \| _ {p}}
$$

completing the proof of the theorem.

![](images/golub_050_099__6cfa699ebce92d9436c4211b9f5d8bbea0e71ee4f692d4d62cbe202a66daa824.jpg)

Note that $\begin{array} { r } { \| ( I - F ) ^ { - 1 } - I \| _ { p } \ \leq \ \| F \| _ { p } / ( 1 - \| F \| _ { p } ) } \end{array}$ is a consequence of the lemma. Thus, $ { \mathrm { i f } } \epsilon \ll 1$ , then $O ( \epsilon )$ perturbations to the identity matrix induce $O ( \epsilon )$ perturbations in the inverse. In general, we have

Theorem 2.3.4. If A is nonsingular and $r \equiv \parallel A ^ { - 1 } E \parallel _ { p } < 1$ , then $A { + } E$ is nonsingular and

$$
\parallel (A + E) ^ {- 1} - A ^ {- 1} \parallel_ {p} \leq \frac {\parallel E \parallel_ {p} \parallel A ^ {- 1} \parallel_ {p} ^ {2}}{1 - r}.
$$

Proof. Note that $A + E = ( I + F ) A$ where $F = - E A ^ { - 1 }$ . Since  $F \parallel _ { p } = r < 1$ , it follows from Lemma 2.3.3 that $I + F$ is nonsingular and $\| { \bf \zeta } ( I + { \cal F } ) ^ { - 1 } \| _ { p } \dot { \leq } 1 / ( 1 - r )$ .

Thus, $( A + E ) ^ { - 1 } = A ^ { - 1 } ( I + F ) ^ { - 1 }$ is nonsingular and

$$
(A + E) ^ {- 1} - A ^ {- 1} = A ^ {- 1} (A - (A + E)) (A + E) ^ {- 1} = - A ^ {- 1} E A ^ {- 1} (I + F) ^ {- 1}.
$$

The theorem follows by taking norms.

# 2.3.5 Orthogonal Invariance

If $A \in \mathbb { R } ^ { m \times n }$ and the matrices $Q \in \mathbb { R } ^ { m \times m }$ and $Z \in \mathbb { R } ^ { n \times n }$ are orthogonal, then

$$
\left\| Q A Z \right\| _ {F} = \left\| A \right\| _ {F} \tag {2.3.14}
$$

and

$$
\| Q A Z \| _ {2} = \| A \| _ {2}. \tag {2.3.15}
$$

These properties readily follow from the orthogonal invariance of the vector 2-norm. For example,

$$
\parallel Q A \parallel_ {F} ^ {2} = \sum_ {j = 1} ^ {n} \parallel Q A (:, j) \parallel_ {2} ^ {2} = \sum_ {j = 1} ^ {n} \parallel A (:, j) \parallel_ {2} ^ {2} = \parallel A \parallel_ {F} ^ {2}
$$

and so $\left\| Q ( A Z ) \right\| _ { F } ^ { 2 } = \left\| \left( A Z \right) \right\| _ { F } ^ { 2 } = \left\| Z ^ { T } A ^ { T } \right\| _ { F } ^ { 2 } = \left\| A ^ { T } \right\| _ { F } ^ { 2 } = \left\| A \right\| _ { F } ^ { 2 } .$

# Problems

P2.3.1 Show $\left. A B \right. _ { p } \leq \left. A \right. _ { p } \left. B \right. _ { p }$ where $1 \leq p \leq \infty$ .

P2.3.2 Let B be any submatrix of A. Show that $\| B \| _ { p } \leq \| A \| _ { p } .$

P2.3.3 Show that if $D = \operatorname { d i a g } ( \mu _ { 1 } , \dots , \mu _ { k } ) \in \mathbb { R } ^ { m \times n }$ with $k = \operatorname* { m i n } \{ m , n \}$ , then $\| D \| _ { p } = \operatorname* { m a x } | \mu _ { i } |$

P2.3.4 Verify (2.3.7) and (2.3.8).

P2.3.5 Verify (2.3.9) and (2.3.10).

P2.3.6 Verify (2.3.11) and (2.3.12).

P2.3.7 Show that if $0 \neq s \in \mathbb { R } ^ { n }$ and $E \in \mathbb { R } ^ { n \times n }$ , then

$$
\left\| E \left(I - \frac {s s ^ {T}}{s ^ {T} s}\right) \right\| _ {F} ^ {2} = \parallel E \parallel_ {F} ^ {2} - \frac {\parallel E s \parallel_ {2} ^ {2}}{s ^ {T} s}.
$$

P2.3.8 Suppose $u \in \mathbb { R } ^ { m }$ and $v \in \mathbb { R } ^ { n }$ . Show that if $E = u v ^ { T }$ , then $\parallel E \parallel _ { F } = \parallel E \parallel _ { 2 } = \parallel u \parallel _ { 2 } \parallel v \parallel _ { 2 }$ and $\left\| \textbf { } E \right\| _ { \infty } \leq \left\| \textbf { } u \right\| _ { \infty } \left\| \textbf { } v \right\| _ { 1 } .$

P2.3.9 Suppose $A \in \mathbb { R } ^ { m \times n } , ~ y \in \mathbb { R } ^ { m }$ , and $0 \neq s \in \mathbb { R } ^ { n }$ . Show that $E = ( y - A s ) s ^ { T } / s ^ { T } ;$ s has the smallest 2-norm of all m-by-n matrices E that satisfy $( A + E ) s = y .$ .

P2.3.10 Verify that there exists a scalar $c > 0$ such that

$$
\| A \| _ {\Delta , c} = \max _ {i, j} c | a _ {i j} |
$$

satisfies the submultiplicative property (2.3.4) for matrix norms on $\mathbb { R } ^ { n \times n }$ . What is the smallest value for such a constant? Referring to this value as $c _ { * }$ , exhibit nonzero matrices B and $C$ with the property that $\| B C \| _ { \Delta , c _ { * } } = \| B \| _ { \Delta , c _ { * } } \| C \| _ { \Delta , c _ { * } }$ .

P2.3.11 Show that if A and B are matrices, then $\parallel A \otimes B \parallel _ { F } = \parallel A \parallel _ { F } \parallel B \parallel _ { F }$ .


---

<!-- golub_100_149 -->

For further discussion of matrix norms, see Stewart (IMC) as well as:

F.L. Bauer and C.T. Fike (1960). “Norms and Exclusion Theorems,” Numer. Math. 2, 137–144.

L. Mirsky (1960). “Symmetric Gauge Functions and Unitarily Invariant Norms,” Quart. J. Math. 11, 50–59.

A.S. Householder (1964). The Theory of Matrices in Numerical Analysis, Dover Publications, New York.

N.J. Higham (1992). “Estimating the Matrix p-Norm,” Numer. Math. 62, 539–556.
