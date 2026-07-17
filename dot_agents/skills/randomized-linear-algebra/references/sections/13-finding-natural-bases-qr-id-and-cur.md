# Finding natural bases: QR, ID, and CUR {#sec:natural}

In Section [11](#sec:random-rangefinder){reference-type="ref"
reference="sec:random-rangefinder"}, we explored a number of efficient
techniques for building a tall thin matrix $\bm{\mathsf{Q}}$ whose
columns form an approximate basis for the range of an input matrix
$\bm{\mathsf{A}}$ that is numerically rank-deficient. The columns of
$\bm{\mathsf{Q}}$ are orthonormal, and they are formed as linear
combinations of many columns from the matrix $\bm{\mathsf{A}}$.

It is sometimes desirable to work with a basis for the range that
consists of a subset of the columns of $\bm{\mathsf{A}}$ itself. In this
case, one typically has to give up on the requirement that the basis
vectors be orthogonal. We gain the advantage of a basis that shares
properties with the original matrix, such as sparsity or nonnegativity.
Moreover, for purposes of data interpretation and analysis, it can be
very useful to identify a subset of the columns that distills the
information in the matrix.

In this section, we start by describing some popular matrix
decompositions that use "natural" basis vectors for the column space,
for the row space, or for both. We show how these matrices can be
computed somewhat efficiently by means of slight modifications to
classical deterministic techniques. Then we describe how to combine
deterministic and randomized methods to obtain algorithms with superior
performance.

## The CUR decomposition, and three flavors of interpolative decompositions {#sec:CURandIDdef}

To introduce the low-rank factorizations that we investigate in this
section, we describe how they can be used to represent an $m\times n$
matrix $\bm{\mathsf{A}}$ of *exact* rank $k$, where $k < \min(m,n)$.
This is an artificial setting, but it allows us to convey the key ideas
using a minimum of notational overhead.

A basic interpolative decomposition (ID) of a matrix $\bm{\mathsf{A}}$
with exact rank $k$ takes the form $$\begin{equation}
\label{eq:defID1}
\begin{array}{cccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{C}} &\bm{\mathsf{Z}},\\
m\times n && m\times k & k\times n
\end{array}
\end{equation}$$ where the matrix $\bm{\mathsf{C}}$ is given by a subset
of the columns of $\bm{\mathsf{A}}$ and where $\bm{\mathsf{Z}}$ is a
matrix that contains the $k\times k$ identity matrix as a submatrix. The
fact that the decomposition
[\[eq:defID1\]](#eq:defID1){reference-type="eqref"
reference="eq:defID1"} exists is an immediate consequence of the
definition of rank. A more significant observation is that there exists
a factorization of the form
[\[eq:defID1\]](#eq:defID1){reference-type="eqref"
reference="eq:defID1"} that is *well-conditioned*, in the sense that no
entry of $\bm{\mathsf{Z}}$ is larger than one in modulus. This claim can
be established through an application of Cramer's rule. See [@tyrt1997]
and [@2006_martinsson_skeletonization].

The factorization [\[eq:defID1\]](#eq:defID1){reference-type="eqref"
reference="eq:defID1"} uses a subset of the columns of $\bm{\mathsf{A}}$
to span its column space. Of course, there is an analog factorization
that uses a subset of the rows of $\bm{\mathsf{A}}$ to span the row
space. We write this as $$\begin{equation}
\label{eq:defID2}
\begin{array}{cccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{X}} &\bm{\mathsf{R}},\\
m\times n && m\times k & k\times n
\end{array}
\end{equation}$$ where $\bm{\mathsf{R}}$ is a matrix consisting of $k$
rows of $\bm{\mathsf{A}}$, and where $\bm{\mathsf{X}}$ is a matrix that
contains the $k\times k$ identity matrix.

For bookkeeping purposes, we introduce index vectors $J_{\mathrm{s}}$
and $I_{\mathrm{s}}$ that identify the columns and rows chosen in the
factorizations [\[eq:defID1\]](#eq:defID1){reference-type="eqref"
reference="eq:defID1"} and
[\[eq:defID2\]](#eq:defID2){reference-type="eqref"
reference="eq:defID2"}. To be precise, let
$J_{\mathrm{s}} \subset \{1,2,\dots,n\}$ denote the index vector of
length $k$ such that
$$\bm{\mathsf{C}} = \bm{\mathsf{A}}(:,J_{\mathrm{s}}).$$ Analogously, we
let $I_{\mathrm{s}} \subset \{1,2,\dots,m\}$ denote the index vector for
which $$\bm{\mathsf{R}} = \bm{\mathsf{A}}(I_{\mathrm{s}},:).$$ The index
vectors $I_{\mathrm{s}}$ and $J_{\mathrm{s}}$ are often referred to as
*skeleton* index vectors, whence the subscript "s". This terminology
arises from the original literature about these
factorizations [@tyrt1997].

A related two-sided factorization is based on extracting a row/column
submatrix. In this case, the basis vectors for the row and column space
are less interpretable. More precisely, $$\begin{equation}
\label{eq:defID3}
\begin{array}{ccccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{X}} & \bm{\mathsf{A}}_{\mathrm{s}} & \bm{\mathsf{Z}},\\
m\times n && m\times k & k\times k & k\times n
\end{array}
\end{equation}$$ where $\bm{\mathsf{X}}$ and $\bm{\mathsf{Z}}$ are the
same matrices as those that appear in
[\[eq:defID1\]](#eq:defID1){reference-type="eqref"
reference="eq:defID1"} and
[\[eq:defID2\]](#eq:defID2){reference-type="eqref"
reference="eq:defID2"}, and where $\bm{\mathsf{A}}_{\mathrm{s}}$ is the
$k\times k$ submatrix of $\bm{\mathsf{A}}$ given by
$$\bm{\mathsf{A}}_{\mathrm{s}} = \bm{\mathsf{A}}(I_{\mathrm{s}},J_{\mathrm{s}}).$$
To distinguish among these variants, we refer to
[\[eq:defID1\]](#eq:defID1){reference-type="eqref"
reference="eq:defID1"} as a *column ID*, to
[\[eq:defID2\]](#eq:defID2){reference-type="eqref"
reference="eq:defID2"} as a *row ID*, and to
[\[eq:defID3\]](#eq:defID3){reference-type="eqref"
reference="eq:defID3"} as a *double-sided ID*.

We introduce a fourth factorization, often called the *CUR
decomposition*. For a matrix of exact rank $k$, it takes the form
$$\begin{equation}
\label{eq:CUR}
\begin{array}{ccccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{C}} & \bm{\mathsf{U}} & \bm{\mathsf{R}},\\
m\times n && m\times k & k\times k & k\times n
\end{array}
\end{equation}$$ where $\bm{\mathsf{C}}$ and $\bm{\mathsf{R}}$ are the
matrices that appeared in
[\[eq:defID1\]](#eq:defID1){reference-type="eqref"
reference="eq:defID1"} and
[\[eq:defID2\]](#eq:defID2){reference-type="eqref"
reference="eq:defID2"}, which consist of $k$ columns and $k$ rows of
$\bm{\mathsf{A}}$, and where $\bm{\mathsf{U}}$ is a small matrix that
links them together. In the present case, where $\bm{\mathsf{A}}$ has
exact rank $k$, the matrix $\bm{\mathsf{U}}$ must take the form
$$\begin{equation}
\label{eq:CURformulaeasy}
\bm{\mathsf{U}} = \bigl(\bm{\mathsf{A}}(I_{\mathrm{s}},J_{\mathrm{s}})\bigr)^{-1}.
\end{equation}$$ The factorizations
[\[eq:defID3\]](#eq:defID3){reference-type="eqref"
reference="eq:defID3"} and [\[eq:CUR\]](#eq:CUR){reference-type="eqref"
reference="eq:CUR"} are related through the formula $$\bm{\mathsf{A}} =
\bm{\mathsf{X}}\bm{\mathsf{A}}_{\mathrm{s}}\bm{\mathsf{Z}} =
\underbrace{\bigl(\bm{\mathsf{X}}\bm{\mathsf{A}}_{\mathrm{s}}\bigr)}_{=\bm{\mathsf{C}}}\,
\underbrace{\bm{\mathsf{A}}_{\mathrm{s}}^{-1}}_{=\bm{\mathsf{U}}}\,
\underbrace{\bigl(\bm{\mathsf{A}}_{\mathrm{s}}\bm{\mathsf{Z}}\bigr)}_{=\bm{\mathsf{R}}}.$$

Comparing the two formats, we see that the CUR
[\[eq:CUR\]](#eq:CUR){reference-type="eqref" reference="eq:CUR"} has an
advantage in that it requires very little storage. As long as
$\bm{\mathsf{A}}$ is stored explicitly (or is easy to retrieve), the CUR
factorization [\[eq:CUR\]](#eq:CUR){reference-type="eqref"
reference="eq:CUR"} is determined by the index vectors $I_{\mathrm{s}}$
and $J_{\mathrm{s}}$ and the linking matrix $\bm{\mathsf{U}}$. If
$\bm{\mathsf{A}}$ is not readily available, then, in order to use the
CUR, we need to evaluate and store the matrices $\bm{\mathsf{C}}$ and
$\bm{\mathsf{R}}$. When $\bm{\mathsf{A}}$ is sparse, the latter approach
can still be more efficient than storing the matrices $\bm{\mathsf{X}}$
and $\bm{\mathsf{Z}}$.

A disadvantage of the CUR factorization
[\[eq:CUR\]](#eq:CUR){reference-type="eqref" reference="eq:CUR"} is
that, when the singular values of $\bm{\mathsf{A}}$ decay rapidly, the
factorization [\[eq:CUR\]](#eq:CUR){reference-type="eqref"
reference="eq:CUR"} is typically numerically ill-conditioned. The reason
is that, whenever the factorization is a good representation of
$\bm{\mathsf{A}}$, the singular values of $\bm{\mathsf{A}}_{\mathrm{s}}$
should approximate the $k$ dominant singular values of
$\bm{\mathsf{A}}$, so the singular values of $\bm{\mathsf{U}}$ end up
approximating the *inverses* of these singular values. This means that
$\bm{\mathsf{U}}$ will have elements of magnitude $1/\sigma_{k}$, which
is clearly undesirable when $\sigma_{k}$ is small. In contrast, the ID
[\[eq:defID3\]](#eq:defID3){reference-type="eqref"
reference="eq:defID3"} is numerically benign.

In the numerical analysis literature, what we refer to as an
interpolative decomposition is often called a *skeleton factorization*
of $\bm{\mathsf{A}}$. This term dates back at least as far as , where
the term *pseudo-skeleton* was used for the CUR decomposition
[\[eq:CUR\]](#eq:CUR){reference-type="eqref" reference="eq:CUR"}.

::: remark
**Remark 43** (Storage-efficient ID). *We mentioned that the matrices
$\bm{\mathsf{X}}$ and $\bm{\mathsf{Z}}$ that appear in the ID are almost
invariably dense, which appears to necessitate the storage of $(m+n)k$
floating-point numbers for the double-sided ID. Observe, however, that
these matrices satisfy the relations
$$\bm{\mathsf{X}} = \bm{\mathsf{C}}\bm{\mathsf{A}}_{\mathrm{s}}^{-1},
\qquad\mbox{and}\qquad
\bm{\mathsf{Z}} = \bm{\mathsf{A}}_{\mathrm{s}}^{-1}\bm{\mathsf{R}}.$$
This means that as long as we store the index vectors $I_{\mathrm{s}}$
and $J_{\mathrm{s}}$, the matrices $\bm{\mathsf{X}}$ and
$\bm{\mathsf{Z}}$ can be applied on the fly whenever needed, and do not
need to be explicitly formed.*
:::

## Approximate rank

In practical applications, the situation we considered in Section
[13.1](#sec:CURandIDdef){reference-type="ref"
reference="sec:CURandIDdef"} where a matrix has exact rank $k$ is rare.
Instead, we typically work with a matrix whose singular values decay
fast enough that it is advantageous to form a low-rank approximation.
Both the ID and the CUR can be used in this environment, but now the
discussion becomes slightly more involved.

To illustrate, let us consider a situation where we are given a
tolerance $\varepsilon$, and we seek to compute an approximation
$\bm{\mathsf{A}}_{k}$ of rank $k$, with $k$ as small as possible, such
that
$\Vert  \bm{\mathsf{A}} - \bm{\mathsf{A}}_{k}  \Vert \leq \varepsilon$.
If $\bm{\mathsf{A}}_k$ is a truncated singular value decomposition, then
the Eckart--Young theorem implies that the rank $k$ of the approximation
$\bm{\mathsf{A}}_k$ will be minimal. When we use an approximate
algorithm, such as the RSVD
(Section [11.2](#sec:rsvd){reference-type="ref" reference="sec:rsvd"}),
we may not find the exact optimum, but we typically get very close.

What happens if we seek an ID $\bm{\mathsf{A}}_{k}$ that approximates
$\bm{\mathsf{A}}$ to a fixed tolerance? There is no guarantee that the
rank $k$ (that is, the number of rows or columns involved) will be close
to the rank of the truncated SVD. How close can we get in practice?

When the singular values of $\bm{\mathsf{A}}$ decay rapidly, then the
minimal rank attainable by an approximate ID is close to what is
attainable with an SVD. Moreover, the algorithms we will describe for
computing an ID produce an answer that is close to the optimal one.

When the singular values decay slowly, however, the difference in rank
between the optimal ID and the optimal SVD can be quite substantial
[@gu1996]. On top of that, the algorithms used to compute the ID can
result in answers that are still further away from the optimal value
[@2005_martinsson_skel].

When the CUR decomposition is used in an environment of approximate
rank, standard algorithms start by determining index sets
$I_{\mathrm{s}}$ and $J_{\mathrm{s}}$ that identify the spanning rows
and columns, and then proceed to the problem of finding a "good" linking
matrix $\bm{\mathsf{U}}_{\mathrm{s}}$. One could still use the formula
[\[eq:CURformulaeasy\]](#eq:CURformulaeasy){reference-type="eqref"
reference="eq:CURformulaeasy"}, but this is rarely a good idea. The most
obvious reason is that the matrix
$\bm{\mathsf{A}}(I_{\mathrm{s}},J_{\mathrm{s}})$ need not be invertible
in this situation. Indeed, when randomized sampling is used to find the
index sets, it is common practice to compute index vectors that hold
substantially more elements than is theoretically necessary, which can
easily make $\bm{\mathsf{A}}(I_{\mathrm{s}},J_{\mathrm{s}})$ singular
or, at the very least, highly ill-conditioned. In this case, a better
approximation is given by $$\begin{equation}
\label{eq:CURformula}
\bm{\mathsf{U}} = \bm{\mathsf{C}}^{\dagger} \bm{\mathsf{A}} \bm{\mathsf{R}}^{\dagger},
\end{equation}$$ with $\bm{\mathsf{C}}^{\dagger}$ and
$\bm{\mathsf{R}}^{\dagger}$ the pseudoinverses of $\bm{\mathsf{C}}$ and
$\bm{\mathsf{R}}$. (As always, pseudoinverses should be applied
numerically by computing a QR or SVD factorization.)

## Deterministic methods, and the connection to column-pivoted QR {#sec:detID}

A substantial amount of research effort has been dedicated to the
question of how to find a set of good spanning columns and/or rows of a
given matrix. It is known that the task of finding the absolutely
optimal one is combinatorially hard, but efficient algorithms exist that
are guaranteed to produce a close-to-optimal answer [@gu1996]. In this
subsection, we briefly discuss some deterministic methods that work well
for dense matrices of modest size. In Section
[13.4](#sec:randID){reference-type="ref" reference="sec:randID"}, we
will show how these methods can be combined with randomized techniques
to arrive at algorithms that work well for general matrices, whether
they are small or huge, sparse or dense, available explicitly or not,
etc.

Perhaps the most obvious deterministic method for computing an ID is the
classical Gram--Schmidt process, which selects the columns or rows in a
greedy fashion. Say we are interested in the column ID
[\[eq:defID1\]](#eq:defID1){reference-type="eqref"
reference="eq:defID1"} of a given matrix $\bm{\mathsf{A}}$. The
Gram--Schmidt procedure first grabs the largest column and places it in
the first column of $\bm{\mathsf{C}}$. Then it projects the remaining
columns onto the orthogonal complement of the one that was picked. It
places the largest of the resulting columns in the second column of
$\bm{\mathsf{C}}$, and so on.

In the traditional numerical linear algebra literature, it is customary
to formulate the Gram--Schmidt process as a column-pivoted QR (CPQR)
decomposition. After $k$ steps, this factorization results in a partial
decomposition of $\bm{\mathsf{A}}$ such that $$\begin{equation}
\label{eq:partialCPQR}
\begin{array}{cccccccccccccccc}
\bm{\mathsf{A}} & \bm{\mathsf{\Pi}} &=& \bm{\mathsf{Q}} & \bm{\mathsf{S}} & + & \bm{\mathsf{E}}, \\
m\times n & n\times n&& m\times k & k\times n && m\times n
\end{array}
\end{equation}$$ where the columns of $\bm{\mathsf{Q}}$ form an
orthonormal basis for the space spanned by the $k$ selected columns of
$\bm{\mathsf{A}}$, where $\bm{\mathsf{S}}$ is upper-triangular, where
$\bm{\mathsf{E}}$ is a "remainder matrix" holding what remains of the
$n-k$ columns of $\bm{\mathsf{A}}$ that have not yet been picked, and
where $\bm{\mathsf{\Pi}}$ is a permutation matrix that reorders the
columns of $\bm{\mathsf{A}}$ in such a way that the $k$ columns picked
are the first $k$ columns of $\bm{\mathsf{A}}\bm{\mathsf{\Pi}}$. (We use
the letter for the upper-triangular factor in lieu of the more
traditional to avoid confusion with the matrix $\bm{\mathsf{R}}$ holding
spanning rows in [\[eq:defID2\]](#eq:defID2){reference-type="eqref"
reference="eq:defID2"} and [\[eq:CUR\]](#eq:CUR){reference-type="eqref"
reference="eq:CUR"}.)

In order to convert
[\[eq:partialCPQR\]](#eq:partialCPQR){reference-type="eqref"
reference="eq:partialCPQR"} into the ID
[\[eq:defID1\]](#eq:defID1){reference-type="eqref"
reference="eq:defID1"}, we split off the first $k$ columns of
$\bm{\mathsf{S}}$ into a $k\times k$ upper-triangular matrix
$\bm{\mathsf{S}}_{11}$, so that $$\bm{\mathsf{S}} =
\kbordermatrix{&
k & n-k\\
k &\bm{\mathsf{S}}_{11} & \bm{\mathsf{S}}_{12}}.$$ Upon multiplying
[\[eq:partialCPQR\]](#eq:partialCPQR){reference-type="eqref"
reference="eq:partialCPQR"} by $\bm{\mathsf{\Pi}}^{*}$ from the right,
we obtain $$\begin{equation}
\label{eq:postCPQR}
\bm{\mathsf{A}} =
\underbrace{\bm{\mathsf{Q}}\bm{\mathsf{S}}_{11}}_{=:\bm{\mathsf{C}}}\,
\underbrace{\bigl[\bm{\mathsf{I}}_{k}\qquad \bm{\mathsf{S}}_{11}^{-1}\bm{\mathsf{S}}_{12}\bigr]\bm{\mathsf{\Pi}}^{*}}_{=:\bm{\mathsf{Z}}}
+
\bm{\mathsf{E}}\bm{\mathsf{\Pi}}^{*}.
\end{equation}$$ We recognize equation
[\[eq:postCPQR\]](#eq:postCPQR){reference-type="eqref"
reference="eq:postCPQR"} as the ID
[\[eq:defID1\]](#eq:defID1){reference-type="eqref"
reference="eq:defID1"}, with the only difference that there is now a
remainder term that results from the fact that $\bm{\mathsf{A}}$ is only
approximately rank-deficient. (Observe that the remainder terms in
[\[eq:partialCPQR\]](#eq:partialCPQR){reference-type="eqref"
reference="eq:partialCPQR"} and in
[\[eq:postCPQR\]](#eq:postCPQR){reference-type="eqref"
reference="eq:postCPQR"} are identical, up to a permutation of the
columns.)

A row ID can obviously be computed by applying Gram--Schmidt to the rows
of $\bm{\mathsf{A}}$ instead of the columns. Alternatively, one may
express this as a column-pivoted QR factorization of
$\bm{\mathsf{A}}^{*}$ instead of $\bm{\mathsf{A}}$.

In order to build a double-sided ID, one starts by computing a
single-sided ID. If $m\geq n$, it is best to start with a column ID of
$\bm{\mathsf{A}}$ to determine $J_{\mathrm{s}}$ and $\bm{\mathsf{Z}}$.
Then we perform a row ID *on the rows of
$\bm{\mathsf{A}}(:,J_{\mathrm{s}})$* to determine $I_{\mathrm{s}}$ and
$\bm{\mathsf{X}}$.

Finally, in order to build a CUR factorization of $\bm{\mathsf{A}}$, we
can easily convert the double-sided ID to a CUR factorization using
[\[eq:CURformulaeasy\]](#eq:CURformulaeasy){reference-type="eqref"
reference="eq:CURformulaeasy"} or
[\[eq:CURformula\]](#eq:CURformula){reference-type="eqref"
reference="eq:CURformula"}.

Detailed descriptions of all algorithms can be found in Sections 10 & 11
of , while analysis and numerical results are given in . A related set
of deterministic techniques that are efficient and often result in
slightly higher quality spanning sets than column pivoting are described
in . Techniques based on optimized spanning volumes of submatrices are
described in [@2010_oseledets_tyrtyshnikov_find_submatrix] and
[@2012_thurau_deterministic_CUR].

::: remark
**Remark 44** (Quality of ID). *In this section, we have described
simple methods based on the column-pivoted QR factorization for
computing a CUR decomposition, as well as all three flavors of
interpolatory decompositions. In discussing the quality of the resulting
factorizations, we will address two questions: (1) How close to minimal
is the resulting approximation error? (2) How well-conditioned are the
basis matrices?*

*The NLA literature contains a detailed study of both questions. This
inquiry was instigated by Kahan's construction of matrices for which
CPQR performs very poorly [@1966_kahan_NLA Sec. 5]. provided a
comprehensive analysis of the situation and presented an algorithm whose
asymptotic complexity in typical environments is only slightly worse
than that of CPQR and that is guaranteed to produce near-optimal
results.*

*In practice, CPQR works well. In almost all cases, it yields
factorizations that are close to optimal. Moreover, it gives
well-conditioned factorizations as long as orthonormality of the basis
is scrupulously maintained
[@2006_martinsson_skeletonization; @2005_martinsson_skel].*

*A more serious problem with the ID and the CUR is that these
decompositions can exhibit much larger approximation errors than the SVD
when the input matrix has slowly decaying singular values. This issue
persists even when the optimal index sets are used.*
:::

## Randomized methods for finding natural bases {#sec:randID}

The deterministic techniques for computing an ID or a CUR decomposition
in Section [13.3](#sec:detID){reference-type="ref"
reference="sec:detID"} work very well for small, dense matrices. In this
section, we describe randomized methods that work much better for
matrices that are sparse or are just very large.

To be concrete, we consider the problem of finding a vector
$I_{\mathrm{s}}$ that identifies a set of rows that form a good basis
for the row space of a given matrix $\bm{\mathsf{A}}$. To do so, we use
the randomized rangefinder to build a matrix $\bm{\mathsf{Y}}$ whose
columns accurately span the column space of $\bm{\mathsf{A}}$ as in
Section [11](#sec:random-rangefinder){reference-type="ref"
reference="sec:random-rangefinder"}. Since $\bm{\mathsf{Y}}$ is far
smaller than $\bm{\mathsf{A}}$, we can use the deterministic methods in
Section [13.3](#sec:detID){reference-type="ref" reference="sec:detID"}
to find a set $I_{\mathrm{s}}$ of rows of $\bm{\mathsf{Y}}$ that form a
basis for the row space of $\bm{\mathsf{Y}}$. Next, we establish a
simple but perhaps non-obvious fact: the set $I_{\mathrm{s}}$ also
identifies a set of rows of $\bm{\mathsf{A}}$ that form a good basis for
the row space of $\bm{\mathsf{A}}$.

To simplify the argument, let us first suppose that we are given an
$m\times n$ matrix $\bm{\mathsf{A}}$ of *exact* rank $k$, and that we
have determined by some means (say, the randomized rangefinder) an
$m\times k$ matrix $\bm{\mathsf{Y}}$ whose columns span the column space
of $\bm{\mathsf{A}}$. Then $\bm{\mathsf{A}}$ admits by definition a
factorization $$\begin{equation}
\label{eq:YF}
\begin{array}{cccccccccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{Y}} & \bm{\mathsf{F}}, \\
m\times n && m\times k & k\times n
\end{array}
\end{equation}$$ for some matrix $\bm{\mathsf{F}}$. Now compute a row ID
of $\bm{\mathsf{Y}}$, by performing Gram--Schmidt on its rows, as
described in Section [13.3](#sec:detID){reference-type="ref"
reference="sec:detID"}. The result is a matrix $\bm{\mathsf{X}}$ and an
index vector $I_{\mathrm{s}}$ such that $$\begin{equation}
\label{eq:Yid}
\begin{array}{cccccccccc}
\bm{\mathsf{Y}} &=& \bm{\mathsf{X}} & \bm{\mathsf{Y}}(I_{\mathrm{s}},:). \\
m\times k && m\times k & k\times k
\end{array}
\end{equation}$$ The claim is now that
$\{I_{\mathrm{s}},\bm{\mathsf{X}}\}$ is automatically a row ID of
$\bm{\mathsf{A}}$ as well. To prove this, observe that $$\begin{align*}
\bm{\mathsf{X}}\bm{\mathsf{A}}(I_{\mathrm{s}},:) &= \bm{\mathsf{X}}\bm{\mathsf{Y}}(I_{\mathrm{s}},:)\bm{\mathsf{F}} &&\mbox{\{use \eqref{eq:YF} restricted to the rows in $I_{\mathrm{s}}$\}} \\
&= \bm{\mathsf{Y}}\bm{\mathsf{F}} && \mbox{\{use \eqref{eq:Yid}\}} \\
&= \bm{\mathsf{A}} &&\mbox{\{use \eqref{eq:YF}\}}
\end{align*}$$ The key insight here is simple and powerful: *In order to
compute a row ID of a matrix $\bm{\mathsf{A}}$, the only information
needed is a matrix $\bm{\mathsf{Y}}$ whose columns span the column space
of $\bm{\mathsf{A}}$.*

The task of finding a matrix $\bm{\mathsf{Y}}$ such that
[\[eq:YF\]](#eq:YF){reference-type="eqref" reference="eq:YF"} holds to
high accuracy is particularly well suited for the randomized rangefinder
described in Section [11](#sec:random-rangefinder){reference-type="ref"
reference="sec:random-rangefinder"}. Putting everything together, we
obtain Algorithm
[\[alg:randomizedID\]](#alg:randomizedID){reference-type="ref"
reference="alg:randomizedID"}. When a Gaussian random matrix is used,
::: {#alg:randomizedID .algorithm}
**Algorithm. Randomized ID.**

```
**Input:** Matrix $\mtx{A} \in \F^{m\times n}$, target rank $k$, oversampling parameter $p$.
**Output:** An $m\times k$ interpolation matrix $\mtx{X}$ and an index vector $I_{\rm s}$ such that $\mtx{A} \approx \mtx{X}\mtx{A}(I_{\rm s},\colon)$.

function RandomizedID($\mtx{A}$, k, p):
    Draw an $n\times (k+p)$ test matrix $\mtx{\Omega}$, e.g., from a Gaussian distribution
    Form the sample matrix $\mtx{Y} = \mtx{A}\mtx{\Omega}$    # Powering may be used
    Form an ID of the $n\times (k+p)$ sample matrix: $[I_{\rm s},\mtx{X}] = ID_row(\mtx{Y},k)$
```
:::


the method has complexity $O(mnk)$.

::: remark
**Remark 45** ($O(mn\log k)$ complexity methods). *An interesting thing
happens if we replace the Gaussian random matrix $\bm{\mathsf{\Omega}}$
in Algorithm
[\[alg:randomizedID\]](#alg:randomizedID){reference-type="ref"
reference="alg:randomizedID"} with a structured random matrix, as
described in Section [9](#sec:dimension-reduction){reference-type="ref"
reference="sec:dimension-reduction"}: Then $\bm{\mathsf{Y}}$ is computed
at cost $O(mn\log k)$, and every step after that has cost
$O((m+n)k^{2})$ or less.*
:::

:::: algorithm
::: algorithmic
Matrix $\bm{\mathsf{A}} \in \mathbb{F}^{m\times n}$, target rank $k$,
oversampling parameter $p$. An $m\times k$ interpolation matrix
$\bm{\mathsf{X}}$ and an index vector $I_{\mathrm{s}}$ such that
$\bm{\mathsf{A}} \approx \bm{\mathsf{X}}\bm{\mathsf{A}}(I_{\mathrm{s}},\colon)$.

Draw an $n\times (k+p)$ test matrix $\bm{\mathsf{\Omega}}$, e.g., from a
Gaussian distribution Form the sample matrix
$\bm{\mathsf{Y}} = \bm{\mathsf{A}}\bm{\mathsf{\Omega}}$ Powering may be
used Form an ID of the $n\times (k+p)$ sample matrix:
$[I_{\mathrm{s}},\bm{\mathsf{X}}] = \texttt{ID\_row}(\bm{\mathsf{Y}},k)$
:::
::::

## Techniques based on coordinate sampling

To find natural bases for a matrix, it is tempting just to sample
coordinates from some probability distribution on the full index vector.
Some advantages and disadvantages of this approach were discussed in
Section [9.6](#sec:coord-embed){reference-type="ref"
reference="sec:coord-embed"}.

In the current context, the main appeal of coordinate sampling is that
the cost is potentially lower than the techniques described in this
section---provided that we do not need to expend much effort to compute
the sampling probabilities. This advantage can be decisive in
applications where mixing random embeddings are too expensive.

Coordinate sampling has several disadvantages in comparison to using
mixing random embeddings. Coordinate sampling typically results in worse
approximations for a given budget of rows or columns. Moreover, the
quality of the approximation obtained from coordinate sampling tends to
be highly variable. These vulnerabilities are less pronounced when the
matrix has very low coherence, so that uniform sampling works well.
There are also a few specialized situations where we can compute
subspace leverage scores efficiently. (Section
[9.6.1](#sec:coherence_leverage){reference-type="ref"
reference="sec:coherence_leverage"} defines coherence and leverage
scores.)

In certain applications, a hybrid approach can work well. First, form an
initial approximation by drawing a very large subset of columns using a
cheap coordinate sampling method. Then slim it down using the techniques
described here, based on mixing random embeddings. An example of this
methodology appears in .

There is a distinct class of methods, based on coresets, that explicitly
takes advantage of coordinate structure for computing matrix
approximations. For example, see . These techniques can be useful for
processing enormous matrices that are very sparse. On the other hand,
they may require larger sets of basis vectors to achieve the same
quality of approximation.

