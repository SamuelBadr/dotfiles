# Linear algebra preliminaries {#sec:lin-alg}

This section contains an overview of the linear algebra tools that arise
in this survey. It collects the basic notation, along with some standard
and not-so-standard definitions. It also contains a discussion about the
role of the spectral norm.

Background references for linear algebra and matrix analysis include and
. For a comprehensive treatment of matrix computations, we refer to
[@GVL13:Matrix-Computations-4ed; @1997_trefethen_bau; @1998_stewart_volume1; @1998_stewart_volume2].

## Basics

We will work in the real field ($\mathbb{R}$) or the complex field
($\mathbb{C}$). The symbol $\mathbb{F}$ refers to either the real or
complex field, in cases where the precise choice is unimportant. As
usual, scalars are denoted by lowercase italic Roman ($a, b, c$) or
Greek ($\alpha, \beta$) letters.

Vectors are elements of $\mathbb{F}^n$, where $n$ is a natural number.
We always denote vectors with lowercase bold Roman
($\bm{\mathsf{a}}, \bm{\mathsf{b}}, \bm{\mathsf{u}}, \bm{\mathsf{v}}$)
or Greek ($\bm{\mathsf{\alpha}}, \bm{\mathsf{\beta}}$) letters. We write
$\bm{\mathsf{0}}$ for the zero vector and $\bm{\mathsf{1}}$ for the
vector of ones. The standard basis vectors are denoted as
$\bm{\mathsf{\delta}}_1, \dots, \bm{\mathsf{\delta}}_n$. The dimensions
of these special vectors are determined by context.

A general matrix is an element of $\mathbb{F}^{m \times n}$, where
$m, n$ are natural numbers. We always denote matrices with uppercase
bold Roman ($\bm{\mathsf{A}}, \bm{\mathsf{B}}, \bm{\mathsf{C}}$) or
Greek ($\bm{\mathsf{\Delta}}, \bm{\mathsf{\Lambda}}$) letters. We write
$\bm{\mathsf{0}}$ for the zero matrix and $\bm{\mathsf{I}}$ for the
identity matrix; their dimensions are determined by a subscript or by
context.

The parenthesis notation is used for indexing into vectors and matrices:
$(\bm{\mathsf{a}})_i$ is the $i$th coordinate of vector
$\bm{\mathsf{a}}$, while $(\bm{\mathsf{A}})_{ij}$ is the $(i, j)$th
coordinate of matrix $\bm{\mathsf{A}}$. In some cases it is more
convenient to invoke the functional form of indexing. For example,
$\bm{\mathsf{A}}(i, j)$ also refers to the $(i, j)$th coordinate of the
matrix $\bm{\mathsf{A}}$.

The colon notation is used to specify ranges of coordinates. For
example, $(\bm{\mathsf{a}})_{1:i}$ and $\bm{\mathsf{a}}(1:i)$ refer to
the vector comprising the first $i$ coordinates of $\bm{\mathsf{a}}$.
The colon by itself refers to the entire range of coordinates. For
instance, $(\bm{\mathsf{A}})_{i:}$ denotes the $i$th row of
$\bm{\mathsf{A}}$, while $(\bm{\mathsf{A}})_{:j}$ denotes the $j$th
column.

The symbol ${}^*$ is reserved for the (conjugate) transpose of a matrix
of vector. A matrix that satisfies $\bm{\mathsf{A}} = \bm{\mathsf{A}}^*$
is said to be self-adjoint. It is convenient to distinguish the space
$\mathbb{H}_n$ of self-adjoint $n \times n$ matrices over the scalar
field. We may write $\mathbb{H}_n(\mathbb{F})$ if it is necessary to
specify the field.

The operator ${}^\dagger$ extracts the Moore--Penrose pseudoinverse of a
matrix. More precisely, for
$\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$, the pseudoinverse
$\bm{\mathsf{A}}^\dagger\in \mathbb{F}^{n \times m}$ is the unique
matrix that satisfies the following:

1.  $\bm{\mathsf{AA}}^\dagger$ is self-adjoint.

2.  $\bm{\mathsf{A}}^\dagger\bm{\mathsf{A}}$ is self-adjoint.

3.  $\bm{\mathsf{AA}}^\dagger\bm{\mathsf{A}} = \bm{\mathsf{A}}$.

4.  $\bm{\mathsf{A}}^\dagger\bm{\mathsf{A}} \bm{\mathsf{A}}^\dagger= \bm{\mathsf{A}}^\dagger$.

If $\bm{\mathsf{A}}$ has full column rank, then
$\bm{\mathsf{A}}^\dagger= (\bm{\mathsf{A}}^* \bm{\mathsf{A}})^{-1} \bm{\mathsf{A}}^*$,
where $(\cdot)^{-1}$ denotes the ordinary matrix inverse.

## Eigenvalues and singular values

A positive semidefinite matrix is a self-adjoint matrix with nonnegative
eigenvalues. We will generally abbreviate positive semidefinite to PSD.
Likewise, a positive definite (PD) matrix is a self-adjoint matrix with
positive eigenvalues.

The symbol $\preccurlyeq$ denotes the semidefinite order on self-adjoint
matrices. The relation $\bm{\mathsf{A}} \preccurlyeq\bm{\mathsf{B}}$
means that $\bm{\mathsf{B}} - \bm{\mathsf{A}}$ is psd.

We write $\lambda_1 \geq \lambda_2 \geq \dots$ for the eigenvalues of a
self-adjoint matrix. We write $\sigma_1 \geq \sigma_2 \geq \dots$ for
the singular values of a general matrix. If the matrix is not clear from
context, we may include it in the notation so that
$\sigma_{j}(\bm{\mathsf{A}})$ is the $j$th singular value of
$\bm{\mathsf{A}}$.

Let $f : \mathbb{R}\to \mathbb{R}$ be a function on the real line. We
can extend $f$ to a spectral function
$f : \mathbb{H}_n \to \mathbb{H}_n$ on (conjugate) symmetric matrices.
Indeed, for a matrix $\bm{\mathsf{A}} \in \mathbb{H}_n$ with eigenvalue
decomposition
$$\bm{\mathsf{A}} = \sum_{i=1}^n \lambda_i \, \bm{\mathsf{u}}_i \bm{\mathsf{u}}_i^*,
\quad\text{we define}\quad
f(\bm{\mathsf{A}}) := \sum_{i=1}^n f(\lambda_i) \, \bm{\mathsf{u}}_i \bm{\mathsf{u}}_i^*.$$
The Pascal notation for definitions ($:=$ and $=:$) is used sparingly,
when we need to emphasize that the definition is taking place.

## Inner product geometry

We equip $\mathbb{F}^n$ with the standard inner product and the
associated $\ell_2$ norm. For all vectors
$\bm{\mathsf{a}}, \bm{\mathsf{b}} \in \mathbb{F}^{n}$,
$$\langle  \bm{\mathsf{a}} , \,  \bm{\mathsf{b}}  \rangle := \bm{\mathsf{a}} \cdot \bm{\mathsf{b}} := \sum_{i=1}^n (\bm{\mathsf{a}})_i^* (\bm{\mathsf{b}})_i
\quad\text{and}\quad
\Vert  \bm{\mathsf{a}}  \Vert^2 := \langle  \bm{\mathsf{a}} , \,  \bm{\mathsf{a}}  \rangle.$$
We write $\mathbb{S}^{n-1}$ for the set of vectors in $\mathbb{F}^n$
with unit $\ell_2$ norm. If needed, we may specify the field:
$\mathbb{S}^{n-1}(\mathbb{F})$.

The trace of a square matrix is the sum of its diagonal entries:
$$\operatorname{trace}(\bm{\mathsf{A}}) := \sum_{i=1}^n (\bm{\mathsf{A}})_{ii}
\quad\text{for $\bm{\mathsf{A}} \in \mathbb{F}^{n \times n}$.}$$
Nonlinear functions bind before the trace. We equip
$\mathbb{F}^{m \times n}$ with the standard trace inner product and the
Frobenius norm. For all matrices
$\bm{\mathsf{A}}, \bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$,
$$\langle  \bm{\mathsf{A}} , \,  \bm{\mathsf{B}}  \rangle := %
    \operatorname{trace}(\bm{\mathsf{A}}^* \bm{\mathsf{B}})
\quad\text{and}\quad
\Vert  \bm{\mathsf{A}}  \Vert_{\mathrm{F}}^2 := \langle  \bm{\mathsf{A}} , \,  \bm{\mathsf{A}}  \rangle.$$
For vectors, these definitions coincide with the ones in the last
paragraph.

We say that a matrix $\bm{\mathsf{U}}$ is *orthonormal* when its columns
are orthonormal with respect to the standard inner product. That is,
$\bm{\mathsf{U}}^* \bm{\mathsf{U}} = \bm{\mathsf{I}}$. If
$\bm{\mathsf{U}}$ is also square, we say instead that $\bm{\mathsf{U}}$
is *orthogonal* ($\mathbb{F} = \mathbb{R}$) or *unitary*
($\mathbb{F} = \mathbb{C}$).

## Norms on matrices

Several different norms on matrices arise during this survey. We use
consistent notation for these norms.

- The unadorned norm $\Vert  \cdot  \Vert$ refers to the spectral norm
  of a matrix, also known as the $\ell_2$ operator norm. It reports the
  maximum singular value of its argument. For vectors, it coincides with
  the $\ell_2$ norm.

- The norm $\Vert  \cdot  \Vert_{*}$ is the nuclear norm of a matrix,
  which is the dual of the spectral norm. It reports the sum of the
  singular values of its argument.

- The symbol $\Vert \cdot \Vert_{\mathrm{F}}$ refers to the Frobenius
  norm, defined in the last subsection. The Frobenius norm coincides
  with the $\ell_2$ norm of the singular values of its argument.

- The notation $\Vert  \cdot  \Vert_{p}$ denotes the Schatten $p$-norm
  for each $p \in [1, \infty]$. The Schatten $p$-norm is the $\ell_p$
  norm of the singular values of its argument. Special cases with their
  own notation include the nuclear norm (Schatten $1$), the Frobenius
  norm (Schatten $2$), and the spectral norm (Schatten ${\infty}$).

Occasionally, other norms may arise, and we will define them explicitly
when they do.

## Approximation in the spectral norm

Throughout this survey, we will almost exclusively use the spectral norm
to measure the error in matrix computations. Let us recall some of the
implications that follow from spectral norm bounds.

Suppose that $\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$ is a matrix,
and $\widehat{\bm{\mathsf{A}}} \in \mathbb{F}^{m \times n}$ is an
approximation. If the approximation satisfies the spectral norm error
bound
$$\Vert  \bm{\mathsf{A}} - \widehat{\bm{\mathsf{A}}}  \Vert \leq \varepsilon,$$
then we can transfer the following information:

- **Linear functionals:**
  $\vert  \langle \bm{\mathsf{F}}, \, \bm{\mathsf{A}} \rangle - \langle \bm{\mathsf{F}}, \, \widehat{\bm{\mathsf{A}}} \rangle  \vert \leq \varepsilon\Vert \bm{\mathsf{F}} \Vert_{*}$
  for every matrix $\bm{\mathsf{F}} \in \mathbb{F}^{m \times n}$

- **Singular values:**
  $\vert  \sigma_j(\bm{\mathsf{A}}) - \sigma_j(\widehat{\bm{\mathsf{A}}})  \vert \leq \varepsilon$
  for each index $j$.

- **Singular vectors:** if the $j$th singular value
  $\sigma_j(\bm{\mathsf{A}})$ is well separated from the other singular
  values, then the $j$th right singular vector of $\bm{\mathsf{A}}$ is
  well approximated by the $j$th right singular vector of
  $\widehat{\bm{\mathsf{A}}}$; a similar statement holds for the left
  singular vectors.

Detailed statements about the singular vectors are complicated, so we
refer the reader to  for his treatment of perturbation of spectral
subspaces.

For one-pass and streaming data models, it may not be possible to obtain
good error bounds in the spectral norm. In this case, we may retrench to
Frobenius norm or nuclear norm error bounds. These estimates give weaker
information about linear functionals, singular values, and singular
vectors.

::: {#rem:frob .remark}
**Remark 1** (Frobenius norm approximation). *In the literature on
randomized NLA, some authors prefer to bound errors with respect to the
Frobenius norm because the arguments are technically simpler. In many
instances, these bounds are less valuable because the error can have the
same scale as the matrix that we wish to approximate.*

*For example, let us consider a variant of the spiked covariance model
that is common in statistics applications [@Joh01:Distribution-Largest].
Suppose we need to approximate a rank-one matrix contaminated with
additive noise:
$\bm{\mathsf{A}} = \bm{\mathsf{uu}}^* + \varepsilon\bm{\mathsf{G}} \in \mathbb{R}^{n \times n}$,
where $\Vert \bm{\mathsf{u}} \Vert = 1$ and
$\bm{\mathsf{G}} \in \mathbb{R}^{n \times n}$ has independent
$\textsc{normal}(0, n^{-1})$ entries. It is well known that
$\Vert \bm{\mathsf{G}} \Vert \approx 2$, while
$\Vert \bm{\mathsf{G}} \Vert_{\mathrm{F}} \approx \sqrt{n}$. With
respect to the Frobenius norm, the zero matrix is almost as good an
approximation of $\bm{\mathsf{A}}$ as the rank-one matrix
$\bm{\mathsf{uu}}^*$:
$$\operatorname{\mathbb{E}}\Vert  \bm{\mathsf{A}} - \bm{\mathsf{uu}}^*  \Vert_{\mathrm{F}}^2 = \varepsilon^2 n
\quad\text{and}\quad
\operatorname{\mathbb{E}}\Vert  \bm{\mathsf{A}} - \bm{\mathsf{0}}  \Vert_{\mathrm{F}}^2 = \varepsilon^2 n + 1.$$
The difference is visible only when the size of the perturbation
$\varepsilon\approx n^{-1/2}$. In contrast, the spectral norm error can
easily distinguish between the good approximation $\bm{\mathsf{uu}}^*$
and the vacuous approximation $\bm{\mathsf{0}}$, even when
$\varepsilon= O(1)$.*

*For additional discussion, see and .*
:::

## Intrinsic dimension and stable rank

Let $\bm{\mathsf{A}} \in \mathbb{H}_n$ be a psd matrix. We define its
*intrinsic dimension*: $$\begin{equation}
 \label{eqn:intdim}
\operatorname{intdim}(\bm{\mathsf{A}}) := \frac{\operatorname{trace}(\bm{\mathsf{A}})}{\Vert \bm{\mathsf{A}} \Vert}.
\end{equation}$$ The intrinsic dimension of a nonzero matrix satisfies
the inequalities
$1 \leq \operatorname{intdim}(\bm{\mathsf{A}}) \leq \operatorname{rank}(\bm{\mathsf{A}})$;
the upper bound is saturated when $\bm{\mathsf{A}}$ is an orthogonal
projector. We can interpret the intrinsic dimension as a continuous
measure of the rank, or the number of energetic dimensions in the
matrix.

Let $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$ be a rectangular
matrix. Its *stable rank* is $$\begin{equation}
 \label{eqn:stable-rank}
\operatorname{srank}(\bm{\mathsf{B}}) := \operatorname{intdim}(\bm{\mathsf{B}}^*\bm{\mathsf{B}}) = \frac{\Vert \bm{\mathsf{B}} \Vert_{\mathrm{F}}^2}{\Vert \bm{\mathsf{B}} \Vert^2}.
\end{equation}$$ Similar to the intrinsic dimension, the stable rank
provides a continuous measure of the rank of $\bm{\mathsf{B}}$.

## Schur complements

Schur complements arise from partial Gaussian elimination and partial
least-squares. They also play a key role in several parts of randomized
NLA. We give the basic definitions here, referring to  for a more
complete treatment.

Let $\bm{\mathsf{A}} \in \mathbb{F}^{n \times n}$ be a psd matrix, and
let $\bm{\mathsf{X}} \in \mathbb{F}^{n \times k}$ be a fixed matrix.
First, define the psd matrix $$\begin{equation}
 \label{eqn:nys-def}
\bm{\mathsf{A}}\langle \bm{\mathsf{X}} \rangle
    := (\bm{\mathsf{AX}})(\bm{\mathsf{X}}^* \bm{\mathsf{A}} \bm{\mathsf{X}})^\dagger(\bm{\mathsf{AX}})^*.
\end{equation}$$ The *Schur complement* of $\bm{\mathsf{A}}$ with
respect to $\bm{\mathsf{X}}$ is the psd matrix $$\begin{equation}
 \label{eqn:schur-complement}
\bm{\mathsf{A}} / \bm{\mathsf{X}} := \bm{\mathsf{A}} - \bm{\mathsf{A}}\langle \bm{\mathsf{X}} \rangle.
\end{equation}$$ The matrices
$\bm{\mathsf{A}}\langle \bm{\mathsf{X}} \rangle$ and
$\bm{\mathsf{A}}/\bm{\mathsf{X}}$ depend on $\bm{\mathsf{X}}$ only
through its range. They also enjoy geometric interpretations in terms of
orthogonal projections with respect to the $\bm{\mathsf{A}}$
semi-inner-product.

## Miscellaneous

We use big-$O$ notation following standard computer science convention.
For instance, we say that a method has (arithmetic) complexity
$O(n^{\omega})$ if there is a finite $C$ for which the number of
floating point operations (flops) expended is bounded by $Cn^{\omega}$
as the problem size $n\rightarrow \infty$.

We use MATLAB-inspired syntax in summarizing algorithms. For instance,
the task of computing an SVD
$\bm{\mathsf{A}} = \bm{\mathsf{U\Sigma V}}^{*}$ of a given matrix
$\bm{\mathsf{A}}$ is written as
$[\bm{\mathsf{U}},\bm{\mathsf{\Sigma}},\bm{\mathsf{V}}]=\texttt{svd}(\bm{\mathsf{A}})$.
We have taken the liberty to modify the syntax when we believe that this
improves clarity. For instance, we write
$[\bm{\mathsf{Q}},\bm{\mathsf{R}}] = \texttt{qr\_econ}(\bm{\mathsf{A}})$
to denote the *economy-size* QR factorization where the matrix
$\bm{\mathsf{Q}}$ has size $m\times \min(m,n)$ for an input matrix
$\bm{\mathsf{A}} \in \mathbb{F}^{m\times n}$. Arguments that are not
needed are replaced by "$\sim$", so that, for example,
$[\bm{\mathsf{Q}},\sim] = \texttt{qr\_econ}(\bm{\mathsf{A}})$ returns
only the matrix $\bm{\mathsf{Q}}$ whose columns form an ON basis for the
range of $\bm{\mathsf{A}}$.

