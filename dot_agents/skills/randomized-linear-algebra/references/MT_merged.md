:::: center
**Randomized Numerical Linear Algebra:\
Foundations & Algorithms**

Per-Gunnar Martinsson, University of Texas at Austin

Joel A. Tropp, California Institute of Technology

::: minipage
**Abstract:** This survey describes probabilistic algorithms for linear
algebra computations, such as factorizing matrices and solving linear
systems. It focuses on techniques that have a proven track record for
real-world problems. The paper treats both the theoretical foundations
of the subject and practical computational issues.

Topics include norm estimation; matrix approximation by sampling;
structured and unstructured random embeddings; linear regression
problems; low-rank approximation; subspace iteration and Krylov methods;
error estimation and adaptivity; interpolatory and CUR factorizations;
Nyström approximation of positive semidefinite matrices; single-view
("streaming") algorithms; full rank-revealing factorizations; solvers
for linear systems; and approximation of kernel matrices that arise in
machine learning and in scientific computing.
:::
::::

# Introduction

Numerical linear algebra (NLA) is one of the great achievements of
scientific computing. On most computational platforms, we can now
routinely and automatically solve small- and medium-scale linear algebra
problems to high precision. The purpose of this survey is to describe a
set of probabilistic techniques that have joined the mainstream of NLA
over the last decade. These new techniques have accelerated everyday
computations for small- and medium-size problems, and they have enabled
large-scale computations that were beyond the reach of classical
methods.

## Classical numerical linear algebra

NLA definitively treats several major classes of problems, including

- solution of dense and sparse linear systems;

- orthogonalization, least-squares, and Tikhonov regularization;

- determination of eigenvalues, eigenvectors, and invariant subspaces;

- singular value decomposition (SVD) and total least-squares.

In spite of this catalog of successes, important challenges remain. The
sheer scale of certain datasets (terabytes and beyond) makes them
impervious to classical NLA algorithms. Modern computing architectures
(GPUs, multi-core CPUs, massively distributed systems) are powerful, but
this power can only be unleashed by algorithms that minimize data
movement and that are designed *ab initio* with parallel computation in
mind. New ways to organize and present data (out-of-core, distributed,
streaming) also demand alternative techniques.

Randomization offers novel tools for addressing all of these challenges.
This paper surveys these new ideas, provides detailed descriptions of
algorithms with a proven track record, and outlines the mathematical
techniques used to analyze these methods.

## Randomized algorithms emerge

Probabilistic algorithms have held a central place in scientific
computing ever since Ulam and von Neumann's groundbreaking work on Monte
Carlo methods in the 1940s. For instance, Monte Carlo algorithms are
essential for high-dimensional integration and for solving PDEs set in
high-dimensional spaces. They also play a major role in modern machine
learning and uncertainty quantification.

For many decades, however, numerical analysts regarded randomized
algorithms as a method of last resort---to be invoked only in the
absence of an effective deterministic alternative. Indeed, probabilistic
techniques have several undesirable features. First, Monte Carlo methods
often produce output with low accuracy. This is a consequence of the
central limit theorem, and in many situations it cannot be avoided.
Second, many computational scientists have a strong attachment to the
engineering principle that two successive runs of the same algorithm
should produce identical results. This requirement aids with debugging,
and it can be critical for applications where safety is paramount, for
example simulation of infrastructure or control of aircraft. Randomized
methods do not generally offer this guarantee. (Controlling the seed of
the random number generator can provide a partial work-around, but this
necessarily involves additional complexity.)

Nevertheless, in the 1980s, randomized algorithms started to make
inroads into NLA. Some of the early work concerns spectral computations,
where it was already traditional to use random initialization.
recognized that (a variant of) the power method with a random start
*provably* approximates the largest eigenvalue of a positive
semidefinite (PSD) matrix, even without a gap between the first and
second eigenvalue. provided a sharp analysis of this phenomenon for both
the power method and the Lanczos algorithm. Around the same time, and
proposed Monte Carlo methods for estimating the trace of a large psd
matrix. Soon after, demonstrated that randomized transformations can be
used to avoid pivoting steps in Gaussian elimination.

Starting in the late 1990s, researchers in theoretical computer science
identified other ways to apply probabilistic algorithms in NLA. and
showed that randomized embeddings allow for computations on streaming
data with limited storage. [@PRTV00:Latent-Semantic] and
[@FKV04:Fast-Monte-Carlo] proposed Monte Carlo methods for low-rank
matrix approximation. [@DKM06:Fast-Monte-Carlo-I],
[@DKM06:Fast-Monte-Carlo-II], and [@DKM06:Fast-Monte-Carlo-III] wrote
the first statement of theoretical principles for randomized NLA. showed
how subspace embeddings support linear algebra computations.

In the mid-2000s, numerical analysts introduced practical randomized
algorithms for low-rank matrix approximation and least-squares problems.
This work includes the first computational evidence that randomized
algorithms outperform classical NLA algorithms for particular classes of
problems. Early contributions include
[@2006_martinsson_random1_orig; @2007_martinsson_PNAS; @RT08:Fast-Randomized; @WLRT08:Fast-Randomized].
These papers inspired later work, such as
[@2010_avron_BLENDENPIK; @HMT11:Finding-Structure; @HMST11:Algorithm-Principal],
that has made a direct impact in applications

Parallel with the advances in numerical analysis, a tide of enthusiasm
for randomized algorithms has flooded into cognate fields. In
particular, stochastic gradient descent [@Bot10:Large-Scale-Machine] has
become a standard algorithm for solving large optimization problems in
machine learning.

At the time of writing, in late 2019, randomized algorithms have joined
the mainstream of NLA. They now appear in major reference works and
textbooks [@GVL13:Matrix-Computations-4ed; @2019_strang_LA_and_data].
Key methods are being incorporated into standard software libraries
[@2019_NAG_library27; @2017_gu_langou; @2017_ghysels_robust].

## What does randomness accomplish? {#sec:accomplishments}

Over the course of this survey we will explore a number of different
ways that randomization can be used to design effective NLA algorithms.
For the moment, let us just summarize the most important benefits.

Randomized methods can handle certain NLA problems faster than any
classical algorithm. In
Section [10](#sec:overdet-ls){reference-type="ref"
reference="sec:overdet-ls"}, we describe a randomized algorithm that can
solve a dense $m \times n$ least-squares problem with $m \gg n$ using
about $O(mn + n^3)$ arithmetic operations [@RT08:Fast-Randomized].
Meanwhile, classical methods require $O(mn^2)$ operations. In
Section [18](#sec:sparse-cholesky){reference-type="ref"
reference="sec:sparse-cholesky"}, we present an algorithm called
[SparseCholesky]{.smallcaps} that can solve the Poisson problem on a
dense undirected graph in time that is roughly *quadratic* in the number
of vertices [@KS16:Approximate-Gaussian]. Standard methods have cost
that is *cubic* in the number of vertices. The improvements can be even
larger for sparse graphs.

Randomization allows us to tackle problems that otherwise seem
impossible. Section [15](#sec:singlepass){reference-type="ref"
reference="sec:singlepass"} contains an algorithm that can compute a
rank-$r$ truncated SVD of an $m \times n$ matrix in a single pass over
the data using working storage $O(r (m + n))$. The first reference for
this kind of algorithm is . We know of no classical method with this
computational profile.

From an engineering point of view, randomization has another crucial
advantage: it allows us to restructure NLA computations in a
fundamentally different way. In
Section [11](#sec:random-rangefinder){reference-type="ref"
reference="sec:random-rangefinder"}, we will introduce the randomized
SVD algorithm [@2006_martinsson_random1_orig; @HMT11:Finding-Structure].
Essentially all the arithmetic in this procedure takes place in a short
sequence of matrix--matrix multiplications. Matrix multiplication is a
highly optimized primitive on most computer systems; it parallelizes
easily; and it performs particularly well on modern hardware such as
GPUs. In contrast, classical SVD algorithms require either random access
to the data or sequential matrix--vector multiplications. As a
consequence, the randomized SVD can process matrices that are beyond the
reach of classical SVD algorithms.

## Algorithm design considerations

Before we decide what algorithm to use for a linear algebra computation,
we must ask how we are permitted to interact with the data. A recurring
theme of this survey is that randomization allows us to reorganize
algorithms so that they control whichever computational resource is the
most scarce (flops, communication, matrix entry evaluation, etc.). Let
us illustrate with some representative examples:

- *Streaming computations ("single-view"):* There is rising demand for
  algorithms that can treat matrices that are so large that they cannot
  be stored at all; other applications involve matrices that are
  presented dynamically. In the streaming setting, the input matrix
  $\bm{\mathsf{A}}$ is given by a sequence of simple linear updates that
  can viewed only once: $$\begin{equation}
   \label{eqn:stream-intro}
  \bm{\mathsf{A}} = \bm{\mathsf{H}}_1 + \bm{\mathsf{H}}_2 + \bm{\mathsf{H}}_3 + \dots.
  \end{equation}$$ We must discard each innovation $\bm{\mathsf{H}}_i$
  after it has been processed. As it happens, the *only* type of
  algorithm that can handle the
  model [\[eqn:stream-intro\]](#eqn:stream-intro){reference-type="eqref"
  reference="eqn:stream-intro"} is one based on randomized linear
  dimension reduction [@LNW14:Turnstile-Streaming]. Our survey describes
  a number of algorithms that can operate in the streaming setting; see
  Sections [4](#sec:trace-est){reference-type="ref"
  reference="sec:trace-est"}, [5](#sec:schatten-p){reference-type="ref"
  reference="sec:schatten-p"}, [14](#sec:nystrom){reference-type="ref"
  reference="sec:nystrom"},
  and [15](#sec:singlepass){reference-type="ref"
  reference="sec:singlepass"}.

- *Dense matrices stored in RAM:* One traditional computational model
  for NLA assumes that the input matrix is stored in fast memory, so
  that any entry can quickly be read and/or overwritten as needed. The
  ability of CPUs to perform arithmetic operations keeps growing
  rapidly, but memory latency has not kept up. Thus, it has become
  essential to formulate blocked algorithms that operate on submatrices.
  Section [16](#sec:full){reference-type="ref" reference="sec:full"}
  shows how randomization can help.

- *Large sparse matrices:* For sparse matrices, it is natural to search
  for techniques that interact with a matrix only through its
  application to vectors, such as Krylov methods or subspace iteration.
  Randomization expands the design space for these methods. When the
  iteration is initialized with a random matrix, we can reach provably
  correct and highly accurate results after a few iterations; see
  Sections [11.6](#sec:rrf-subspace){reference-type="ref"
  reference="sec:rrf-subspace"} and
  [11.7](#sec:rrf-krylov){reference-type="ref"
  reference="sec:rrf-krylov"}.

  Another idea is to apply randomized sampling to control sparsity
  levels. This technique arises in Section
  [18](#sec:sparse-cholesky){reference-type="ref"
  reference="sec:sparse-cholesky"}, which contains a randomized
  algorithm that accelerates incomplete Cholesky preconditioning for
  sparse graph Laplacians.

- *Matrices for which entry evaluation is expensive:* In machine
  learning and computational physics, it is often desirable to solve
  linear systems where it is too expensive to evaluate the full
  coefficient matrix. Randomization offers a systematic way to extract
  data from the matrix and to compute approximations that serve for the
  downstream applications. See
  Sections [19](#sec:kernel){reference-type="ref"
  reference="sec:kernel"} and
  [20](#sec:rankstructured){reference-type="ref"
  reference="sec:rankstructured"}.

## Overview

This paper covers fundamental mathematical ideas, as well as algorithms
that have proved to be effective in practice. The balance shifts from
theory at the beginning toward computational practice at the end. With
the practitioner in mind, we have attempted to make the algorithmic
sections self-contained, so that they can be read with a minimum of
references to other parts of the paper.

After introducing notation and covering preliminaries from linear
algebra and probability in
Sections [2](#sec:lin-alg){reference-type="ref"
reference="sec:lin-alg"}--[3](#sec:prob){reference-type="ref"
reference="sec:prob"}, the survey covers the following topics.

- Sections [4](#sec:trace-est){reference-type="ref"
  reference="sec:trace-est"}--[5](#sec:schatten-p){reference-type="ref"
  reference="sec:schatten-p"} discuss algorithms for trace estimation
  and Schatten $p$-norm estimation based on randomized sampling (i.e.,
  Monte Carlo methods). Section [6](#sec:max-eig){reference-type="ref"
  reference="sec:max-eig"} shows how iteration can improve the quality
  of estimates for maximum eigenvalues, maximum singular values, and
  trace functions.

- Section [7](#sec:matrix-mc){reference-type="ref"
  reference="sec:matrix-mc"} develops randomized sampling methods for
  approximating matrices, including applications to matrix
  multiplication and approximation of combinatorial graphs.

- Sections [8](#sec:gauss){reference-type="ref"
  reference="sec:gauss"}--[9](#sec:dimension-reduction){reference-type="ref"
  reference="sec:dimension-reduction"} introduce the notion of a
  randomized linear embedding. These maps are frequently used to reduce
  the dimension of a set of vectors, while preserving their geometry. In
  Section [10](#sec:overdet-ls){reference-type="ref"
  reference="sec:overdet-ls"}, we explore several ways to use randomized
  embeddings in the context of an overdetermined least-squares problem.

- Sections [11](#sec:random-rangefinder){reference-type="ref"
  reference="sec:random-rangefinder"}--[12](#sec:error-est){reference-type="ref"
  reference="sec:error-est"} demonstrate how randomized methods can be
  used to find a subspace that is aligned with the range of a matrix and
  to assess the quality of this subspace.
  Sections [13](#sec:natural){reference-type="ref"
  reference="sec:natural"}, [14](#sec:nystrom){reference-type="ref"
  reference="sec:nystrom"},
  and [15](#sec:singlepass){reference-type="ref"
  reference="sec:singlepass"} show how to use this subspace to compute a
  variety of low-rank matrix approximations.

- Section [16](#sec:full){reference-type="ref" reference="sec:full"}
  develops randomized algorithms for computing a factorization of a
  full-rank matrix, such as a pivoted QR decomposition or a URV
  decomposition.

- Section [17](#sec:linear-solve){reference-type="ref"
  reference="sec:linear-solve"} describes some general approaches to
  solving linear systems using randomized techniques.
  Section [18](#sec:sparse-cholesky){reference-type="ref"
  reference="sec:sparse-cholesky"} presents the
  [SparseCholesky]{.smallcaps} algorithm for solving the Poisson problem
  on an undirected graph (i.e., a linear system in a graph Laplacian).

- Last, Sections [19](#sec:kernel){reference-type="ref"
  reference="sec:kernel"}
  and [20](#sec:rankstructured){reference-type="ref"
  reference="sec:rankstructured"} show how to use randomized methods to
  approximate kernel matrices that arise in machine learning,
  computational physics, and scientific computing.

## Omissions

While randomized NLA was a niche topic 15 years ago, we have seen an
explosion of research over the last decade. This survey can only cover a
small subset of the many important and interesting ideas that have
emerged.

Among many other omissions, we do not discuss spectral computations in
detail. There have been interesting and very recent developments,
especially for the challenging problem of computing a spectral
decomposition of a nonnormal matrix [@BGKS19:Pseudospectral-Shattering].
We also had to leave out a treatment of tensors and the rapidly
developing field of randomized multilinear algebra.

There is no better way to demonstrate the value of a numerical method
than careful numerical experiments that measure its speed and accuracy
against state-of-the-art implementations of competing methods. Our
selection of topics to cover is heavily influenced by such comparisons;
for reasons of space, we have often had to settle for citations to the
literature, instead of including the numerical evidence.

The intersection between optimization and linear algebra is of crucial
importance in applications, and it remains a fertile ground for
theoretical work. Randomized algorithms are invaluable in this context,
but we realized early on that the paper would double in length if we
included just the essential facts about randomized optimization
algorithms.

There is a complementary perspective on randomized numerical analysis
algorithms, called *probabilistic numerics*. See the
website [@HO:Probabilistic-Numerics] for a comprehensive bibliography.

For the topics that we do cover, we have made every effort to include
all essential citations. Nevertheless, the literature is vast, and we
are sure to have overlooked important work; we apologize in advance for
these oversights.

## Other surveys

There are a number of other survey papers on randomized algorithms for
randomized NLA and related topics.

- develop and analyse computational methods for low-rank matrix
  approximation (including the "randomized SVD" algorithm) from the
  point of view of a numerical analyst. The main idea is that
  randomization can furnish a subspace that captures the action of a
  matrix, and this subspace can be used to build structured low-rank
  matrix approximations.

- treats randomized methods for least-squares computations and for
  low-rank matrix approximation. He emphasizes the useful principle that
  we can often decouple the linear algebra and the probability when
  analyzing randomized NLA algorithms.

- describes how to use subspace embeddings as a primitive for developing
  randomized linear algebra algorithms. A distinctive feature is the
  development of lower bounds.

- gives an introduction to matrix concentration inequalities, and
  includes some applications to randomized NLA algorithms.

- The survey of appeared in a previous volume of *Acta Numerica*. A
  unique aspect is the discussion of randomized tensor computations.

- have updated the presentation in , and include an introduction to
  linear algebra and probability that is directed toward NLA
  applications.

- focuses on computational aspects of randomized NLA. A distinctive
  feature is the discussion of efficient algorithms for factorizing
  matrices of full, or nearly full, rank.

- gives a mathematical treatment of how matrix concentration supports a
  few randomized NLA algorithms, and includes a complete proof of
  correctness for the [SparseCholesky]{.smallcaps} algorithm described
  in Section [18](#sec:sparse-cholesky){reference-type="ref"
  reference="sec:sparse-cholesky"}.

## Acknowledgments

We are grateful to Arieh Iserles for proposing that we write this
survey. Both authors have benefited greatly from our collaborations with
Vladimir Rokhlin and Mark Tygert. Most of all, we would like to thank
Richard Kueng for his critical reading of the entire manuscript, which
has improved the presentation in many places. Madeleine Udell, Riley
Murray, James Levitt, and Abinand Gopal also gave us useful feedback on
parts of the paper. Lorenzo Rosasco offered invaluable assistance with
the section on kernel methods for machine learning. Navid Azizan, Babak
Hassibi, and Peter Richtárik helped with citations to the literature on
SGD. Finally, we would like to thank our ONR programme managers, Reza
Malek-Madani and John Tague, for supporting research on randomized
numerical linear algebra.

JAT acknowledges support from the Office of Naval Research (awards
N-00014-17-1-2146 and N-00014-18-1-2363). PGM acknowledges support from
the Office of Naval Research (award N00014-18-1-2354), from the National
Science Foundation (award DMS-1620472), and from Nvidia Corp.

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

# Probability preliminaries {#sec:prob}

This section summarizes the key definitions and background from
probability and high-dimensional probability. Later in the survey, we
will present more complete statements of foundational results, as they
are needed.

provide an accessible overview of applied probability. introduces the
field of high-dimensional probability. For more mathematical
presentations, see the classic book of or the lecture notes of .

## Basics

We work in a master probability space that is rich enough to support all
of the random variables that are defined. We will not comment further
about the underlying model.

In this paper, the unqualified term *random variable* encompasses random
scalars, vectors, and matrices. Scalar-valued random variables are
usually (but not always) denoted by uppercase italic Roman letters
($X, Y, Z$). A random vector is denoted by a lowercase bold letter
($\bm{\mathsf{x}}, \bm{\mathsf{\omega}}$). A random matrix is denoted by
an uppercase bold letter
($\bm{\mathsf{X}}, \bm{\mathsf{Y}}, \bm{\mathsf{\Gamma}}, \bm{\mathsf{\Omega}}$).
This notation works in concert with the notation for deterministic
vectors and matrices.

The map $\mathbb{P}(E)$ returns the probability of an event $E$. We
usually specify the event using the compact set builder notation that is
standard in probability. For example,
$\mathbb{P}\left\{  X > t  \right\}$ is the probability that the scalar
random variable $X$ exceeds a level $t$.

The operator $\operatorname{\mathbb{E}}$ returns the expectation of a
random variable. For vectors and matrices, the expectation can be
computed coordinate by coordinate. The expectation is linear, which
justifies relations like
$$\operatorname{\mathbb{E}}[ \bm{\mathsf{AX}} ] = \bm{\mathsf{A}} \operatorname{\mathbb{E}}[ \bm{\mathsf{X}} ]
\quad\text{when $\bm{\mathsf{A}}$ is deterministic and $\bm{\mathsf{X}}$ is random.}$$
We use the convention that nonlinear functions bind before the
expectation; for instance,
$\operatorname{\mathbb{E}}X^2 = \operatorname{\mathbb{E}}[ X^2 ]$. The
operator $\operatorname{Var}[\cdot]$ returns the variance of a scalar
random variable.

We say that a random variable is *centered* when its expectation equals
zero. A random vector $\bm{\mathsf{x}}$ is *isotropic* when
$\operatorname{\mathbb{E}}[ \bm{\mathsf{xx}}^* ] = \bm{\mathsf{I}}$. A
random vector is *standardized* when it is both centered and isotropic.
In particular, a scalar random variable is standardized when it has
expectation zero and variance one.

When referring to independent random variables, we often include the
qualification "statistically independent" to make a distinction with
"linearly independent." We abbreviate the term (statistically)
"independent and identically distributed" as i.i.d.

## Distributions

To refer to a named distribution, we use small capitals. In this
context, the symbol $\sim$ means "has the same distribution as."

We write $\textsc{unif}$ for the uniform distribution over a finite set
(with counting measure). In particular, a scalar Rademacher random
variable has the distribution $\textsc{unif}\{ \pm 1 \}$. A Rademacher
random vector has iid coordinates, each distributed as a scalar
Rademacher random variable. We sometimes require the uniform
distribution over a Borel subset of $\mathbb{F}^n$, equipped with
Lebesgue measure.

We write $\textsc{normal}(\bm{\mathsf{\mu}}, \bm{\mathsf{C}})$ for the
normal distribution on $\mathbb{F}^n$ with expectation
$\bm{\mathsf{\mu}} \in \mathbb{F}^n$ and psd covariance matrix
$\bm{\mathsf{C}} \in \mathbb{H}_n(\mathbb{F})$. A *standard normal*
random variable or random vector has expectation zero and covariance
matrix equal to the identity. We often use the term *Gaussian* to refer
to normal distributions.

## Concentration inequalities

Concentration inequalities provide bounds on the probability that a
random variable is close to its expectation. A good reference for the
scalar case is the book by . In the matrix setting, closeness is
measured in the spectral norm. For an introduction to matrix
concentration, see
[@Tro15:Introduction-Matrix; @Tro19:Matrix-Concentration-LN]. These
results play an important role in randomized linear algebra.

## Gaussian random matrix theory

On several occasions, we use comparison principles to study the action
of a random matrix with iid Gaussian entries. In particular, these
methods can be used to control the largest and smallest singular values.
The main classical comparison theorems are associated with the names
Slepian, Chevet, and Gordon. For accounts of this, see or or . More
recently, it has been observed that Gordon's inequality can be reversed
in certain settings [@TOH14:Gaussian-Minmax].

In several instances, we require more detailed information about
Gaussian random matrices. General resources include and . Most of the
specific results we need are presented in .

# Trace estimation by sampling {#sec:trace-est}

We commence with a treatment of matrix trace estimation problems. These
questions stand among the simplest linear algebra problems because the
desired result is just a scalar. Even so, the algorithms have a vast
sweep of applications, ranging from computational statistics to quantum
chemistry. They also serve as building blocks for more complicated
randomized NLA algorithms. We have chosen to begin our presentation here
because many of the techniques that drive algorithms for more difficult
problems already appear in a nascent---and especially pellucid---form in
this section.

Randomized methods for trace estimation depend on a natural technical
idea: One may construct an unbiased estimator for the trace and then
average independent copies to reduce the variance of the estimate.
Algorithms of this type are often called *Monte Carlo methods*. We
describe how to use standard methods from probability and statistics to
develop *a priori* and *a posteriori* guarantees for Monte Carlo trace
estimators. We show how to use structured random distributions to
improve the computational profile of the estimators. Last, we
demonstrate that trace estimators also yield approximations for the
Frobenius norm and the Schatten 4-norm of a general matrix.

In Section [5](#sec:schatten-p){reference-type="ref"
reference="sec:schatten-p"}, we present more involved Monte Carlo
methods that are required to estimate Schatten $p$-norms for larger
values of $p$, which give better approximations for the spectral norm.
Section [6](#sec:max-eig){reference-type="ref" reference="sec:max-eig"}
describes iterative algorithms that lead to much higher accuracy than
Monte Carlo methods. In Section [6.6](#sec:slq){reference-type="ref"
reference="sec:slq"}, we touch on related probabilistic techniques for
evaluating trace functions.

## Overview

We will focus on the problem of estimating the trace of a nonzero psd
matrix $\bm{\mathsf{A}} \in \mathbb{H}_n$. Our goal is to to produce an
approximation of $\operatorname{trace}(\bm{\mathsf{A}})$, along with a
measure of quality.

Trace estimation is easy in the case where we have inexpensive access to
the entries of the matrix $\bm{\mathsf{A}}$ because we can simply read
off the $n$ diagonal entries. But there are many environments where the
primitive operation is the matrix--vector product
$\bm{\mathsf{u}} \mapsto \bm{\mathsf{A}}\bm{\mathsf{u}}$. For example,
$\bm{\mathsf{u}} \mapsto \bm{\mathsf{A}}\bm{\mathsf{u}}$ might be the
solution of a (discretized) linear differential equation with initial
condition $\bm{\mathsf{u}}$, implemented by some computer program. In
this case, we would really prefer to avoid $n$ applications of the
primitive. (We can obviously compute the trace by applying the primitive
to each standard basis vector $\bm{\mathsf{\delta}}_i$.)

The methods in this section all use linear information about the matrix
$\bm{\mathsf{A}}$. In other words, we will extract data from the input
matrix by computing the product
$\bm{\mathsf{Y}} = \bm{\mathsf{A\Omega}}$, where
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times k}$ is a (random) test
matrix. All subsequent operations involve only the sample matrix
$\bm{\mathsf{Y}}$ and the test matrix $\bm{\mathsf{\Omega}}$. Since the
data collection process is linear, we can apply randomized trace
estimators in the one-pass or streaming environments. Moreover, parts of
these algorithms are trivially parallelizable.

The original application of randomized trace estimation was to perform
*a posteriori* error estimation for large least-squares computations.
More specifically, it was used to accelerate cross-validation procedures
for estimating the optimal regularization parameter in a smoothing
spline [@Gir89:Fast-Monte-Carlo; @Hut90:Stochastic-Estimator]. See  for
a list of contemporary applications in machine learning, uncertainty
quantification, and other fields.

## Trace estimation by randomized sampling {#sec:trace-est-basic}

Randomized trace estimation is based on the insight that it is easy to
construct a random variable whose expectation equals the trace of the
input matrix.

Consider a random test vector $\bm{\mathsf{\omega}} \in \mathbb{F}^n$
that is isotropic:
$\operatorname{\mathbb{E}}[ \bm{\mathsf{\omega}} \bm{\mathsf{\omega}}^* ] = \bm{\mathsf{I}}$.
By the cyclicity of the trace and by linearity, $$\begin{equation}
 \label{eqn:trace-est}
X = \bm{\mathsf{\omega}}^* (\bm{\mathsf{A}} \bm{\mathsf{\omega}})
\quad\text{satisfies}\quad
\operatorname{\mathbb{E}}X = \operatorname{trace}(\bm{\mathsf{A}}).
\end{equation}$$ In other words, the random variable $X$ is an unbiased
estimator of the trace. Note that the distribution of $X$ depends on the
unknown matrix $\bm{\mathsf{A}}$.

A single sample of $X$ is rarely adequate because its variance,
$\operatorname{Var}[X]$, will be large. The most common mechanism for
reducing the variance is to average $k$ independent copies of $X$. For
$k \in \mathbb{N}$, define $$\begin{equation}
 \label{eqn:trace-est-avg}
\bar{X}_k = \frac{1}{k} \sum\nolimits_{i=1}^k X_i
\quad\text{where $X_i \sim X$ are iid.}
\end{equation}$$ By linearity, $\bar{X}_k$ is also an unbiased estimator
of the trace. The individual samples are statistically independent, so
the variance decreases. Indeed,
$$\operatorname{\mathbb{E}}[ \bar{X}_k ] = \operatorname{trace}(\bm{\mathsf{A}})
\quad\text{and}\quad
\operatorname{Var}[ \bar{X}_k ] = \frac{1}{k} \operatorname{Var}[X].$$
The
estimator [\[eqn:trace-est-avg\]](#eqn:trace-est-avg){reference-type="eqref"
reference="eqn:trace-est-avg"} can be regarded as the most elementary
method in randomized linear algebra. See
Algorithm [\[alg:trace-est\]](#alg:trace-est){reference-type="ref"
reference="alg:trace-est"}.
::: {#alg:trace-est .algorithm}
**Algorithm. Trace estimation by random sampling.**

```
**Input:** Psd input matrix $\mtx{A} \in \Sym_n$, number $k$ of samples
**Output:** Trace estimate $\bar{X}_k$ and sample variance $S_k$

function TraceEstimate($\mtx{A}$, $k$):
    for $i = 1, \dots, k$:
        # Compute trace samples
        Draw isotropic test vector $\vct{\omega}_i \in \F^n$
        Compute $X_i = \vct{\omega}_i^* (\mtx{A} \vct{\omega}_i)$
    Form trace estimator: $\bar{X}_k = k^{-1} \sum_{i=1}^k X_i$
    Form sample variance: $S_k = (k-1)^{-1} \sum_{i=1}^k (X_i - \bar{X}_k)^2$
    
    # Use compensated summation techniques for large $k$!
```
:::



To compute $\bar{X}_k$, we must simulate $k$ independent copies of the
random vector $\bm{\mathsf{\omega}} \in \mathbb{F}^n$ and perform $k$
matrix--vector products with $\bm{\mathsf{A}}$, plus $O(kn)$ additional
arithmetic.

::: example
**Example 2** (). *Consider a standard normal random vector
$\bm{\mathsf{\omega}} \sim \textsc{normal}(\bm{\mathsf{0}}, \bm{\mathsf{I}})$.
The variance of the resulting trace
estimator [\[eqn:trace-est\]](#eqn:trace-est){reference-type="eqref"
reference="eqn:trace-est"}--[\[eqn:trace-est-avg\]](#eqn:trace-est-avg){reference-type="eqref"
reference="eqn:trace-est-avg"} satisfies $$\begin{equation}
 \label{eqn:trace-est-girard-var}
\operatorname{Var}[ \bar{X}_k ] = \frac{2}{k} \sum_{i,j = 1}^n \vert  (\bm{\mathsf{A}})_{ij}  \vert^2 = \frac{2}{k} \Vert \bm{\mathsf{A}} \Vert_{\mathrm{F}}^2
    \leq \frac{2}{k} \Vert \bm{\mathsf{A}} \Vert \operatorname{trace}(\bm{\mathsf{A}}).
\end{equation}$$ The rotational invariance of the standard normal
distribution allows us to characterize the behavior of this estimator in
full detail.*
:::

::: example
**Example 3** (). *Consider a Rademacher random vector
$\bm{\mathsf{\omega}} \sim \textsc{unif}\{\pm 1\}^n$. The variance of
the resulting trace
estimator [\[eqn:trace-est\]](#eqn:trace-est){reference-type="eqref"
reference="eqn:trace-est"}--[\[eqn:trace-est-avg\]](#eqn:trace-est-avg){reference-type="eqref"
reference="eqn:trace-est-avg"} satisfies
$$\operatorname{Var}[ \bar{X}_k ] =
\frac{4}{k} \sum_{1 \leq i < j \leq n} \vert  (\bm{\mathsf{A}})_{ij}  \vert^2 <
\frac{2}{k} \Vert \bm{\mathsf{A}} \Vert_{\mathrm{F}}^2 \leq
\frac{2}{k}  \Vert \bm{\mathsf{A}} \Vert \operatorname{trace}(\bm{\mathsf{A}}).$$
This is the minimum variance trace estimator generated by an isotropic
random vector $\bm{\mathsf{\omega}}$ with statistically independent
coordinates. It also avoids the simulation of normal variables.*
:::

also studied the estimator obtained by drawing
$\bm{\mathsf{\omega}} \in \mathbb{F}^n$ uniformly at random from the
sphere $\sqrt{n} \, \mathbb{S}^{n-1}(\mathbb{F})$ for
$\mathbb{F}=\mathbb{R}$. When $\mathbb{F}= \mathbb{C}$, this approach
has the minimax variance among all trace estimators of the
form [\[eqn:trace-est\]](#eqn:trace-est){reference-type="eqref"
reference="eqn:trace-est"}--[\[eqn:trace-est-avg\]](#eqn:trace-est-avg){reference-type="eqref"
reference="eqn:trace-est-avg"}. We return to this example in
Section [4.7.1](#sec:near-isotropic){reference-type="ref"
reference="sec:near-isotropic"}.

:::: algorithm
::: algorithmic
Psd input matrix $\bm{\mathsf{A}} \in \mathbb{H}_n$, number $k$ of
samples Trace estimate $\bar{X}_k$ and sample variance $S_k$

Compute trace samples Draw isotropic test vector
$\bm{\mathsf{\omega}}_i \in \mathbb{F}^n$ Compute
$X_i = \bm{\mathsf{\omega}}_i^* (\bm{\mathsf{A}} \bm{\mathsf{\omega}}_i)$

Form trace estimator: $\bar{X}_k = k^{-1} \sum_{i=1}^k X_i$

Form sample variance:
$S_k = (k-1)^{-1} \sum_{i=1}^k (X_i - \bar{X}_k)^2$

Use compensated summation techniques for large $k$!
:::
::::

::: remark
**Remark 4** (General matrices). *The assumption that $\bm{\mathsf{A}}$
is psd allows us to conclude that the standard deviation of the
randomized trace estimate is smaller than the trace of the matrix. The
same methods allow us to estimate the trace of a general square matrix,
but the variance of the estimator may no longer be comparable with the
trace.*
:::

## A priori error estimates

We can use theoretical analysis to obtain prior guarantees on the
performance of the trace estimator. These results illuminate what
features of the input matrix affect the quality of the trace estimate,
and they tell us how many samples $k$ suffice to achieve a given error
tolerance. Note, however, that these bounds depend on properties of the
input matrix that are often unknown to the user of the trace estimator.

Regardless of the distribution of the isotropic test vector
$\bm{\mathsf{\omega}}$, Chebyshev's inequality delivers a simple
probability bound for the trace estimator: $$\begin{equation}
 \label{eqn:trace-est-chebyshev}
\mathbb{P}\left\{  \vert  \bar{X}_k - \operatorname{trace}(\bm{\mathsf{A}})  \vert \geq t  \right\}
    \leq \frac{ \operatorname{Var}[X] }{k t^2}
    \quad\text{for $t > 0$.}
\end{equation}$$ We can specialize this result to specific trace
estimators by inserting the variance.

::: example
**Example 5** (Girard trace estimator). *If the test vector
$\bm{\mathsf{\omega}}$ is standard normal, the trace estimator
$\bar{X}_k$ satisfies
$$\mathbb{P}\left\{  \vert  \bar{X}_k - \operatorname{trace}(\bm{\mathsf{A}})  \vert \geq t \cdot \operatorname{trace}(\bm{\mathsf{A}})  \right\}
    \leq \frac{ 2 }{k \operatorname{intdim}(\bm{\mathsf{A}}) \, t^2}.$$
The bound follows
from [\[eqn:trace-est-girard-var\]](#eqn:trace-est-girard-var){reference-type="eqref"
reference="eqn:trace-est-girard-var"}, [\[eqn:trace-est-chebyshev\]](#eqn:trace-est-chebyshev){reference-type="eqref"
reference="eqn:trace-est-chebyshev"},
and [\[eqn:intdim\]](#eqn:intdim){reference-type="eqref"
reference="eqn:intdim"}. In words, the trace estimator achieves a
relative error bound that is sharpest when the intrinsic
dimension [\[eqn:intdim\]](#eqn:intdim){reference-type="eqref"
reference="eqn:intdim"} of $\bm{\mathsf{A}}$ is large.*
:::

For specific distributions of the random test vector
$\bm{\mathsf{\omega}}$, we can obtain much stronger probability bounds
for the resulting trace estimator using exponential concentration
inequalities. Here is a recent analysis for Girard's estimator based on
fine properties of the standard normal distribution.

::: {#thm:trace-est-exp .theorem}
**Theorem 6** (). *Let $\bm{\mathsf{A}} \in \mathbb{H}_n(\mathbb{R})$ be
a nonzero psd matrix. Consider the trace
estimator [\[eqn:trace-est\]](#eqn:trace-est){reference-type="eqref"
reference="eqn:trace-est"}--[\[eqn:trace-est-avg\]](#eqn:trace-est-avg){reference-type="eqref"
reference="eqn:trace-est-avg"} obtained from a standard normal test
vector $\bm{\mathsf{\omega}} \in \mathbb{R}^n$. For $\tau > 1$ and
$k \leq n$, $$\begin{aligned}
\mathbb{P}\left\{  \bar{X}_k \geq \tau \operatorname{trace}(\bm{\mathsf{A}})  \right\} &\leq \exp\left( -\tfrac{1}{2} k \operatorname{intdim}(\bm{\mathsf{A}}) \big(\sqrt{\tau} - 1 \big)^2 \right); \\
\mathbb{P}\left\{  \bar{X}_k \leq \tau^{-1} \operatorname{trace}(\bm{\mathsf{A}})  \right\} &\leq \exp\left( -\tfrac{1}{4} k \operatorname{intdim}(\bm{\mathsf{A}}) \big(\tau^{-1} - 1 \big)^2 \right). \\
\end{aligned}$$ When $\bm{\mathsf{A}} \in \mathbb{H}_n(\mathbb{C})$ is
psd and $\bm{\mathsf{\omega}} \in \mathbb{C}^n$ is complex standard
normal, the same bounds hold with an extra factor two in the exponent.
(So the estimator works better in the complex setting.)*
:::

::: proof
*Proof.* (Sketch) Carefully estimate the moment generating function of
the random variable $X$, and use the Cramér--Chernoff method to obtain
the probability inequalities. ◻
:::

## Universality

Empirically, for a large sample, the performance of the trace estimator
$\bar{X}_k$ only depends on the distribution of the test vector
$\bm{\mathsf{\omega}}$ through the variance of the resulting sample $X$.
In a word, the estimator exhibits *universality*. As a consequence, we
can select the distribution that is most convenient for computational
purposes.

Classical probability theory furnishes justification for these claims.
The strong law of large numbers tells us that
$$\bar{X}_k \to \operatorname{trace}(\bm{\mathsf{A}})
\quad\text{almost surely as $k \to \infty$.}$$ Concentration
inequalities [@BLM13:Concentration-Inequalities] allow us to derive
rates of convergence akin to
Theorem [6](#thm:trace-est-exp){reference-type="ref"
reference="thm:trace-est-exp"}.

To understand the sampling distribution of the estimator $\bar{X}_k$, we
can invoke the central limit theorem:
$$\sqrt{k} ( \bar{X}_k - \operatorname{trace}(\bm{\mathsf{A}}) )
    \to \textsc{normal}(0, \operatorname{Var}[X])
    \quad\text{in distribution as $k \to \infty$.}$$ We can obtain
estimates for the rate of convergence to normality using the
Berry--Esséen theorem and its
variants [@Ros11:Fundamentals-Steins; @CGS11:Normal-Approximation].

Owing to the universality phenomenon, we can formally use the normal
limit to obtain heuristic error estimates and insights for trace
estimators constructed with test vectors from any distribution. This
strategy becomes even more valuable when the linear algebra problem is
more complicated.

::: warning
**Warning 7** (CLT). *Estimators based on averaging independent samples
cannot overcome the central limit theorem. Their accuracy will always be
limited by fluctuations on the scale of $\sqrt{\operatorname{Var}[X]}$.
In other words, we must extract $\varepsilon^{-2}$ samples to reduce the
error to $\varepsilon\sqrt{\operatorname{Var}[X]}$ for small
$\varepsilon> 0$. This is the curse of Monte Carlo.*
:::

## A posteriori error estimates {#sec:trace-est-post}

In practice, we rarely have access to all the information required to
activate *a priori* error bounds. It is wiser to assess the quality of
the estimate from the information that we actually collect. Since we
have full knowledge of the random process that generates the trace
estimate, we can confidently use approaches from classical statistics.

At the most basic level, the sample variance is an unbiased estimator
for the variance of the individual samples:
$$S_k = \frac{1}{k-1} \sum_{i=1}^k (X_i - \bar{X}_k)^2
\quad\text{satisfies}\quad
\operatorname{\mathbb{E}}[ S_k ] = \operatorname{Var}[X].$$ The variance
of $S_k$ depends on the fourth moment of the random variable $X$. A
standard estimate is
$$\operatorname{Var}[S_k] \leq \frac{1}{k} \operatorname{\mathbb{E}}[ (X - \operatorname{\mathbb{E}}X)^4 ].$$
Bounds and empirical estimates for the variance of $S_k$ can also be
obtained using the Efron--Stein
inequality [@BLM13:Concentration-Inequalities Sec. 3.1].

For $\alpha \in (0, 1/2)$, we can construct the (symmetric, Student's
$t$) confidence interval at level $1 - 2\alpha$:
$$\operatorname{trace}(\bm{\mathsf{A}}) \in \bar{X}_k \pm t_{\alpha, k-1} \sqrt{S_k}$$
where $t_{\alpha, k-1}$ is the $\alpha$ quantile of the Student's
$t$-distribution with $k-1$ degrees of freedom. We interpret this result
as saying that $\operatorname{trace}(\bm{\mathsf{A}})$ lies in the
specified interval with probability roughly $1 - 2\alpha$ (over the
randomness in the trace estimator). The usual rule of thumb is that the
sample size should be moderate (say, $k \geq 30$), while $\alpha$ cannot
be too small (say, $\alpha \geq 0.025$).

## Bootstrapping the sampling distribution {#sec:trace-est-boot}

Miles Lopes has proposed a sweeping program that uses the bootstrap to
construct data-driven confidence sets for randomized NLA
algorithms [@Lop19:Estimating-Algorithmic]. For trace estimation, this
approach is straightforward to describe and implement; see
Algorithm [\[alg:trace-est-boot\]](#alg:trace-est-boot){reference-type="ref"
reference="alg:trace-est-boot"}.
::: {#alg:trace-est-boot .algorithm}
**Algorithm. Bootstrap confidence interval for trace estimation.**

```
**Input:** Psd input matrix $\mtx{A} \in \Sym_n$, number $k$ of trace samples, number $B$ of bootstrap replicates, parameter $\alpha$ for level of confidence
**Output:** Confidence interval $[\bar{X}_k + q_{\alpha}, \bar{X}_k + q_{1-\alpha}]$ at level $1-2\alpha$

function BootstrapTraceEstimate($\mtx{A}$, $k$, $b$, $\alpha$):
    for $i = 1, \dots, k$:
        # Compute trace estimators
        Draw isotropic test vector $\vct{\omega}_i \in \F^n$
        Form $X_i = \vct{\omega}_i^* (\mtx{A} \vct{\omega}_i)$
    $\mathcal{X} = ( X_1, \dots, X_k )$
    # Collate sample
    $\bar{X}_k = k^{-1} \sum_{i=1}^k X_i$
    # Trace estimate
    for $b = 1, \dots, B$:
        # Bootstrap replicates
        Draw $(X_1^*, \dots, X_k^*)$ from $\mathcal{X}$ with replacement
        Compute $e_b^* = (k^{-1} \sum_{i=1}^k X_i^*) - \bar{X}_k$
    Find $q_\alpha$ and $q_{1-\alpha}$ quantiles of errors $(e_1^*, \dots, e_b^*)$
```
:::



Let $\mathcal{X} = ( X_1, \dots, X_k )$ be the empirical sample
from [\[eqn:trace-est-avg\]](#eqn:trace-est-avg){reference-type="eqref"
reference="eqn:trace-est-avg"}. The bootstrap draws further random
samples from $\mathcal{X}$ to elicit more information about the sampling
distribution of the trace estimator $X$, such as confidence sets.

1.  For each $b = 1, \dots, B$,

    1.  Draw a bootstrap replicate $( X_1^*, \dots, X_k^* )$ uniformly
        from $\mathcal{X}$ with replacement.

    2.  Compute the error estimate $e_b^* = \bar{X}_k^* - \bar{X}_k$,
        where $\bar{X}_k^*$ is the sample average of the bootstrap
        replicate.

2.  Compute quantiles $q_\alpha$ and $q_{1-\alpha}$ of the error
    distribution $( e_1^*, \dots, e_B^* )$.

3.  Report the $1-2\alpha$ confidence set
    $[\bar{X}_k + q_{\alpha}, \bar{X}_k + q_{1 - \alpha}]$.

Typical values are $k \geq 30$ samples and $B \geq 1000$ bootstrap
replicates when $\alpha \geq 0.025$. This method is effective for a wide
range of distributions on the test vector, and it extends to other
problems. See  for an introduction to resampling methods.

:::: algorithm
::: algorithmic
Psd input matrix $\bm{\mathsf{A}} \in \mathbb{H}_n$, number $k$ of trace
samples, number $B$ of bootstrap replicates, parameter $\alpha$ for
level of confidence Confidence interval
$[\bar{X}_k + q_{\alpha}, \bar{X}_k + q_{1-\alpha}]$ at level
$1-2\alpha$

Compute trace estimators Draw isotropic test vector
$\bm{\mathsf{\omega}}_i \in \mathbb{F}^n$ Form
$X_i = \bm{\mathsf{\omega}}_i^* (\bm{\mathsf{A}} \bm{\mathsf{\omega}}_i)$

$\mathcal{X} = ( X_1, \dots, X_k )$ Collate sample
$\bar{X}_k = k^{-1} \sum_{i=1}^k X_i$ Trace estimate

Bootstrap replicates Draw $(X_1^*, \dots, X_k^*)$ from $\mathcal{X}$
with replacement Compute
$e_b^* = (k^{-1} \sum_{i=1}^k X_i^*) - \bar{X}_k$

Find $q_\alpha$ and $q_{1-\alpha}$ quantiles of errors
$(e_1^*, \dots, e_b^*)$
:::
::::

## Structured distributions for test vectors

As discussed, there is a lot of flexibility in designing the
distribution of the test vector. We can exploit this freedom to achieve
additional computational goals. For example, we might:

- minimize the variance $\operatorname{Var}[X]$ of each sample;

- minimize the number of random bits required to construct
  $\bm{\mathsf{\omega}}$;

- design a test vector $\bm{\mathsf{\omega}}$ that is "compatible" with
  the input matrix $\bm{\mathsf{A}}$ to facilitate the matrix--vector
  product. For example, if $\bm{\mathsf{A}}$ has a tensor product
  structure, we might require $\bm{\mathsf{\omega}}$ to share the tensor
  structure.

Let us describe a general construction for test vectors that can help
achieve these desiderata. The ideas come from frame theory and quantum
information theory. The approach here extends the work in .

### Optimal measurement systems {#sec:near-isotropic}

In this section, we work in the complex field. Consider a discrete set
$\mathcal{U} := \{ \bm{\mathsf{u}}_1, \dots, \bm{\mathsf{u}}_m \} \subset \mathbb{C}^n$
of vectors, each with unit $\ell_2$ norm. We say that the $\mathcal{U}$
is an *optimal measurement system* when $$\begin{equation}
 \label{eqn:near-isotropic}
\frac{1}{m} \sum\nolimits_{i=1}^m (\bm{\mathsf{u}}_i^*\bm{\mathsf{M}} \bm{\mathsf{u}}_i) \, \bm{\mathsf{u}}_i \bm{\mathsf{u}}_i^*
    = \frac{1}{(n+1) n} \left[ \bm{\mathsf{M}} + \operatorname{trace}(\bm{\mathsf{M}}) \, \bm{\mathsf{I}}\right]
\end{equation}$$ for all $\bm{\mathsf{M}} \in \mathbb{H}_n(\mathbb{C})$.
The reproducing
property [\[eqn:near-isotropic\]](#eqn:near-isotropic){reference-type="eqref"
reference="eqn:near-isotropic"} shows that the system of vectors
acquires enough information to reconstruct an arbitrary self-adjoint
matrix. A similar definition is valid for an infinite system of unit
vectors, provided we replace the sum
in [\[eqn:near-isotropic\]](#eqn:near-isotropic){reference-type="eqref"
reference="eqn:near-isotropic"} with an integral.

Now, suppose that we draw a random test vector
$\bm{\mathsf{\omega}} = \sqrt{n} \, \bm{\mathsf{u}}$, where
$\bm{\mathsf{u}}$ is drawn uniformly at random from an optimal
measurement system $\mathcal{U}$. Then the resulting trace estimator is
unbiased: $$\begin{equation}
\label{eq:unbi1}
X = \bm{\mathsf{\omega}}^* (\bm{\mathsf{A}} \bm{\mathsf{\omega}})
    \quad\text{satisfies}\quad
    \operatorname{\mathbb{E}}X = \operatorname{trace}(\bm{\mathsf{A}}).
\end{equation}$$ The variance of this trace estimator satisfies
$$\begin{equation}
\label{eq:unbi2}
\operatorname{Var}[X] = \frac{n}{n+1} \left[ \Vert \bm{\mathsf{A}} \Vert_{\mathrm{F}}^2 %
    - \frac{1}{n} \operatorname{trace}(\bm{\mathsf{A}})^2 \right].
\end{equation}$$ The identities
([\[eq:unbi1\]](#eq:unbi1){reference-type="ref" reference="eq:unbi1"})
and ([\[eq:unbi2\]](#eq:unbi2){reference-type="ref"
reference="eq:unbi2"}) follow quickly
from [\[eqn:near-isotropic\]](#eqn:near-isotropic){reference-type="eqref"
reference="eqn:near-isotropic"}. As it happens, this is the minimax
variance achievable for a best isotropic distribution on test
vectors [@Kue19:2-Designs-Minimize]

### Examples

Optimal measurement systems arise in quantum information theory (as
near-isotropic measurement systems), in approximation theory (as
projective 2-designs), and in frame theory (as tight fusion frames). The
core examples are as follows:

1.  A set of $n^2$ equiangular lines in $\mathbb{C}^n$, each spanned by
    a unit vector in
    $\{\bm{\mathsf{u}}_1, \dots, \bm{\mathsf{u}}_{n^2}\}$, gives an
    optimal measurement system. In this case, equiangularity means that
    $\vert \langle \bm{\mathsf{u}}_i, \, \bm{\mathsf{u}}_j \rangle \vert^2 = (d+1)^{-1}$
    whenever $i \neq j$. It is conjectured that these systems exist for
    every natural number $n$.

2.  The columns of a family of $(n + 1)$ mutually unbiased bases in
    $\mathbb{F}^n$ compose an optimal measurement system with $n(n+1)$
    unit vectors. For reference, a pair
    $(\bm{\mathsf{U}}, \bm{\mathsf{V}})$ of $n \times n$ unitary
    matrices is called *mutually unbiased* if
    $\bm{\mathsf{\delta}}_i^* (\bm{\mathsf{U}}^* \bm{\mathsf{V}}) \bm{\mathsf{\delta}}_j = n^{-1}$
    for all $i, j$. (For instance, consider the identity matrix and the
    discrete Fourier transform matrix.) These systems exist whenever $n$
    is a power of a prime number [@WF89:Optimal-State-Determination].

3.  The $\ell_2$ unit sphere $\mathbb{S}^{n-1}(\mathbb{C})$ in
    $\mathbb{C}^n$, equipped with the uniform measure, is a continuous
    optimal measurement system. The real case was studied by , but the
    complex case is actually more natural.

See  for more discussion and an application to quantum state tomography.
provides a good survey of what is currently known about finite optimal
measurement systems.

## Extension: The Frobenius norm and the Schatten 4-norm {#sec:schatten-2-4}

The randomized trace estimators developed in this section can also be
deployed to estimate a couple of matrix norms.

Consider a rectangular matrix
$\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$, accessed via the
matrix--vector product
$\bm{\mathsf{u}} \mapsto \bm{\mathsf{B}} \bm{\mathsf{u}}$. Let us
demonstrate how to estimate the Frobenius norm (i.e., Schatten 2-norm)
and the Schatten 4-norm of $\bm{\mathsf{B}}$.

For concreteness, suppose that we extract test vectors from the standard
normal distribution. Draw a standard normal matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times k}$ with columns
$\bm{\mathsf{\omega}}_i$. Construct the random variable
$$\bar{X}_k := \frac{1}{k} \Vert  \bm{\mathsf{B}} \bm{\mathsf{\Omega}}  \Vert_{\mathrm{F}}^2
    = \frac{1}{k} \sum_{i=1}^k \bm{\mathsf{\omega}}_i^* (\bm{\mathsf{B}}^* \bm{\mathsf{B}}) \bm{\mathsf{\omega}}_i
    =: \frac{1}{k} \sum_{i=1}^k X_i.$$ We compute $\bar{X}_k$ by
simulating $nk$ standard normal variables, taking $k$ matrix--vector
products with $\bm{\mathsf{B}}$, and performing $O(km)$ additional
arithmetic.

To analyze $\bar{X}_k$, note that it is an instance of the randomized
trace estimator, where
$\bm{\mathsf{A}} = \bm{\mathsf{B}}^* \bm{\mathsf{B}}$. In particular,
its statistics are
$$\operatorname{\mathbb{E}}[ \bar{X}_k ] = \Vert \bm{\mathsf{B}} \Vert_{\mathrm{F}}^2
\quad\text{and}\quad
\operatorname{Var}[ \bar{X}_k ] = \frac{2}{k} \Vert \bm{\mathsf{B}} \Vert_{4}^4.$$
We see that $\bar{X}_k$ provides an unbiased estimate for the squared
Frobenius norm of the matrix $\bm{\mathsf{B}}$. Meanwhile, by rescaling
the sample variance $S_k^2$ of the data $( X_1, \dots, X_k)$, we obtain
an unbiased estimate for the fourth power of the Schatten 4-norm of
$\bm{\mathsf{B}}$.

Our discussion shows how we can obtain *a priori* guarantees and *a
posteriori* error estimates for these norm computations; the results can
be expressed in terms of the stable
rank [\[eqn:stable-rank\]](#eqn:stable-rank){reference-type="eqref"
reference="eqn:stable-rank"} of $\bm{\mathsf{B}}$. We can also obtain
norm estimators that are more computationally efficient using structured
distributions for the test vectors, such as elements from an optimal
measurement system.

# Schatten $p$-norm estimation by sampling {#sec:schatten-p}

As we saw in Section [4.8](#sec:schatten-2-4){reference-type="ref"
reference="sec:schatten-2-4"}, we can easily construct unbiased
estimators for the Schatten $2$-norm and the Schatten $4$-norm of a
matrix by randomized sampling. What about the other Schatten norms? This
type of computation can be used to obtain better approximations of the
spectral norm, or it can be combined with the method of moments to
approximate the spectral density of the input matrix
[@KV17:Spectrum-Estimation].

In this section, we show that it is possible to use randomized sampling
to construct unbiased estimators for the Schatten $2p$-norm for each
natural number $p \in \mathbb{N}$. In contrast with the case
$p \in \{1, 2\}$, estimators for $p \geq 3$ are combinatorial. They may
also require a large number of samples to ensure that the variance of
the estimator is controlled.

In the next section, we explain how to use iterative methods to
approximate the spectral norm (i.e., the Schatten $\infty$-norm).
Iterative algorithms also lead to much more reliable estimators for
general Schatten norms.

## Overview

Consider a general matrix $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$,
accessed via the matrix--vector product
$\bm{\mathsf{u}} \mapsto \bm{\mathsf{B}}\bm{\mathsf{u}}$. For a sample
size $k$, let $\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times k}$ be a
(random) test matrix that does not depend on $\bm{\mathsf{B}}$. For a
natural number $p \geq 3$, we consider the problem of estimating the
Schatten $2p$-norm $\Vert \bm{\mathsf{B}} \Vert_{2p}$ from the sample
matrix $\bm{\mathsf{Y}} = \bm{\mathsf{B}}\bm{\mathsf{\Omega}}$. The key
idea is to cleverly process the sample matrix $\bm{\mathsf{Y}}$ to form
an unbiased estimator for the $2p$-th power of the norm. Since the
methods in this section use linear information, they can be
parallelized, and they are applicable in the one-pass and streaming
environments.

Of course, if we are given a singular value decomposition (SVD) of
$\bm{\mathsf{B}}$, it is straightforward to extract the Schatten
$2p$-norm. We are interested in methods that are much less expensive
than the $O(\min\{mn^2, nm^2\})$ cost of computing an SVD with a
classical direct algorithm. Randomized SVD or URV algorithms
(Sections [11.2](#sec:rsvd){reference-type="ref" reference="sec:rsvd"},
[15](#sec:singlepass){reference-type="ref" reference="sec:singlepass"},
and [16](#sec:full){reference-type="ref" reference="sec:full"}) can also
be used for Schatten norm estimation, but this approach is typically
overkill.

## Interlude: Lower bounds

How large a sample size $k$ is required to estimate
$\Vert \bm{\mathsf{B}} \Vert_{2p}$ up to a fixed constant factor with
75% probability over the randomness in $\bm{\mathsf{\Omega}}$? The
answer, unfortunately, turns out to be
$k \gtrsim \min\{m,n\}^{1 - 2/p}$. In other words, for a general matrix,
we cannot approximate the Schatten $2p$-norm for any $p > 2$ unless the
sample size $k$ grows polynomially with the dimension.

A more detailed version of this statement appears in . The authors
exhibit a particularly difficult type of input matrix (a standard normal
matrix with a rank-one spike, as in Remark
[1](#rem:frob){reference-type="ref" reference="rem:frob"}) to arrive at
the negative conclusion.

If you have a more optimistic nature, you can also take the inspiration
that other types of input matrices might be easier to handle. For
example, it is possible to use a small random sample to compute the
Schatten norm of a matrix that enjoys some decay in the singular value
spectrum.

## Estimating Schatten norms the hard way

First, let us describe a technique from classical statistics that leads
to an unbiased estimator of $\Vert \bm{\mathsf{B}} \Vert_{2p}^{2p}$.
This estimator is both highly variable and computationally expensive, so
we must proceed with caution.

For the rest of this section, we assume that the random test matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times k}$ has isotropic columns
$\bm{\mathsf{\omega}}_i$ that are iid. Form the sample matrix
$\bm{\mathsf{Y}} = \bm{\mathsf{B}}\bm{\mathsf{\Omega}}$. Abbreviate
$\bm{\mathsf{A}} = \bm{\mathsf{B}}^*\bm{\mathsf{B}}$ and
$\bm{\mathsf{X}} = \bm{\mathsf{Y}}*\bm{\mathsf{Y}}$. Observe that
$$(\bm{\mathsf{X}})_{ij} = (\bm{\mathsf{Y}}^* \bm{\mathsf{Y}})_{ij} = \bm{\mathsf{\omega}}_i^* \bm{\mathsf{A}} \bm{\mathsf{\omega}}_j.$$
Therefore, for any natural numbers that satisfy
$1 \leq i_1, \dots, i_p \leq k$,
$$(\bm{\mathsf{X}})_{i_1 i_2} (\bm{\mathsf{X}})_{i_2 i_3} \dots (\bm{\mathsf{X}})_{i_p i_1}
    = \operatorname{trace}\big[ \bm{\mathsf{\omega}}_{i_1} \bm{\mathsf{\omega}}_{i_1}^* \bm{\mathsf{A}} \dots \bm{\mathsf{\omega}}_{i_p} \bm{\mathsf{\omega}}_{i_p}^* \bm{\mathsf{A}} \big].$$
If we assume that $i_1, \dots, i_p$ are distinct, we can use
independence and isotropy to compute the expectation:
$$\operatorname{\mathbb{E}}[ (\bm{\mathsf{X}})_{i_1 i_2} (\bm{\mathsf{X}})_{i_2 i_3} \dots (\bm{\mathsf{X}})_{i_p i_1} ]
    = \operatorname{trace}[ \bm{\mathsf{A}}^p ]
    = \Vert  \bm{\mathsf{B}}  \Vert_{2p}^{2p}$$ By averaging over all
sequences of distinct indices, we obtain an unbiased estimator:
$$U_p = \frac{(k-p)!}{k!} \sum_{1 \leq i_1, \dots, i_p \leq k}^{\circ} (\bm{\mathsf{X}})_{i_1 i_2} (\bm{\mathsf{X}})_{i_2 i_3} \dots (\bm{\mathsf{X}})_{i_p i_1}.$$
The circle over the sum indicates that the indices must be distinct.
Since
$\operatorname{\mathbb{E}}U_p = \Vert \bm{\mathsf{B}} \Vert_{2p}^{2p}$,
our hope is that
$$U_p^{1/(2p)} \approx \Vert \bm{\mathsf{B}} \Vert_{2p}. %$$ To ensure
that the approximation is precise, the standard deviation of $U_p$
should be somewhat smaller than the mean of $U_p$.

To understand the statistic $U_p$, we can use tools from the theory of
$U$-statistics [@KB94:Theory-U-Statistics]. For instance, when $p$ is
fixed, we have the limit
$$k \operatorname{Var}[ U_p ] \to p^2 \operatorname{Var}[ \bm{\mathsf{\omega}}^* \bm{\mathsf{A}}^p \bm{\mathsf{\omega}} ]
    \quad\text{as $k \to \infty$.}$$ In particular, if the test vector
$\bm{\mathsf{\omega}} \sim \textsc{normal}(\bm{\mathsf{0}}, \bm{\mathsf{I}})$,
then
$$k \operatorname{Var}[ U_p ] \to 2p^2 \Vert  \bm{\mathsf{B}}  \Vert_{4p}^{4p}
    \quad\text{as $k \to \infty$.}$$ We can reduce the variance further
by using test vectors from an optimal measurement system.

Unfortunately, it is quite expensive to compute the statistic $U_p$
because it involves almost $k^p$ summands. When $p$ is a small constant
(say, $p = 4$ or $p = 5$), it is not too onerous to form $U_p$. On the
other hand, in the worst case, we cannot beat the lower bound
$k \gtrsim \min\{m,n\}^{1-2/p}$, where computation of $U_p$ requires
$O(\min\{m,n\}^p)$ operations. proposes some small economies in this
computation.

## Estimating Schatten norms the easy way {#sec:valiant-kong}

Next, we describe a more recent approach that was proposed by . This
method suffers from even higher variance, but it is computationally
efficient. As a consequence, we can target larger values of $p$ and
exploit a larger sample size $k$.

Superficially, the Kong & Valiant estimator appears similar to the
statistic $U_p$. With the same notation, they restrict their attention
to *increasing* sequences $i_1 < i_2 < \dots < i_p$ of indices. In other
words, $$V_p = \binom{k}{p}^{-1} \sum_{1 \leq i_1 < \dots < i_p \leq k}
    (\bm{\mathsf{X}})_{i_1 i_2} (\bm{\mathsf{X}})_{i_2 i_3} \dots (\bm{\mathsf{X}})_{i_p i_1}.$$
Much as before, $V_p$ gives an unbiased estimator for
$\Vert \bm{\mathsf{B}} \Vert_{2p}^{2p}$.

Although $V_p$ still appears to be combinatorial, the restriction to
increasing sequences allows for a linear algebraic reformulation of the
statistic. Let $\mathcal{T} : \mathbb{H}_k \to \mathbb{F}^{k \times k}$
be the linear map that reports the strict upper triangle of a
(conjugate) symmetric matrix. Then
$$V_p = \binom{k}{p}^{-1} \operatorname{trace}[ \mathcal{T}(\bm{\mathsf{X}})^{p-1} \bm{\mathsf{X}} ].$$
The cost of computing $V_p$ is usually dominated by the $O(k^2 n)$
arithmetic required to form $\bm{\mathsf{X}}$ given $\bm{\mathsf{Y}}$.
See
Algorithm [\[alg:valiant-kong\]](#alg:valiant-kong){reference-type="ref"
reference="alg:valiant-kong"} for the procedure.
::: {#alg:valiant-kong .algorithm}
**Algorithm. Schatten $2p$-norm estimation by random sampling.**

```
**Input:** Input matrix $\mtx{B} \in \F^{m \times n}$, order $p$ of norm, number $k$ of samples
**Output:** Schatten $2p$-norm estimate $V_p$

function SchattenEstimate($\mtx{B}$, $p$, $k$):
    Draw test matrix $\mtx{\Omega} \in \F^{n \times k}$ with iid isotropic columns
    Compute the sample matrix $\mtx{Y} = \mtx{B\Omega}$
    Form the Gram matrix $\mtx{X} = \mtx{Y}^* \mtx{Y} \in \F^{k \times k}$
    Extract the strict upper triangle $\mtx{T} = \mathcal{T}(\mtx{X})$
    Compute $\mtx{T}^{p-1}$ by repeated squaring
    Return $V_p = \trace( \mtx{T}^{p-1} \mtx{X} )$
```
:::



obtain bounds for the variance of the estimator $V_p$ to justify its
employment when the number $k$ of samples satisfies
$k \gtrsim \min\{m,n\}^{1 - 2/p}$. This bound is probably substantially
pessimistic for matrices that exhibit spectral decay, but these
theoretical and computational questions remain open.

:::: algorithm
::: algorithmic
Input matrix $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$, order $p$ of
norm, number $k$ of samples Schatten $2p$-norm estimate $V_p$

Draw test matrix $\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times k}$ with
iid isotropic columns Compute the sample matrix
$\bm{\mathsf{Y}} = \bm{\mathsf{B\Omega}}$

Form the Gram matrix
$\bm{\mathsf{X}} = \bm{\mathsf{Y}}^* \bm{\mathsf{Y}} \in \mathbb{F}^{k \times k}$
Extract the strict upper triangle
$\bm{\mathsf{T}} = \mathcal{T}(\bm{\mathsf{X}})$ Compute
$\bm{\mathsf{T}}^{p-1}$ by repeated squaring Return
$V_p = \operatorname{trace}( \bm{\mathsf{T}}^{p-1} \bm{\mathsf{X}} )$
:::
::::

## Bootstrapping the sampling distribution {#bootstrapping-the-sampling-distribution}

Given the lack of prior guarantees for the estimators $U_p$ and $V_p$
described in this section, we recommend that users apply resampling
methods to obtain empirical information about the sampling distribution.
See for reliable bootstrap procedures for $U$-statistics.

## Extension: Estimating the spectral norm by random sampling

Recall that the spectral norm of $\bm{\mathsf{B}}$ is comparable with
the Schatten $2p$-norm of $\bm{\mathsf{B}}$ for an appropriate choice of
$p$. Indeed, $$%
\Vert  \bm{\mathsf{B}}  \Vert_{2p} \geq \Vert  \bm{\mathsf{B}}  \Vert
\geq \min\{m,n\}^{-1/(2p)} \Vert \bm{\mathsf{B}} \Vert_{2p}
\quad\text{for $p \geq 1/2$.}$$ Thus, the Schatten $2p$-norm is
equivalent to the spectral norm when $p \gtrsim \log \min\{m,n\}$. In
fact, when the matrix $\bm{\mathsf{B}}$ has a decaying singular value
spectrum, the Schatten $2p$-norm may already be comparable with the
spectral norm for much smaller values of $p$.

Thus, we can try to approximate the spectral norm by estimating the
Schatten $2p$-norm for a sufficiently large value of $p$. Resampling
techniques can help ensure that the estimate is reliable. Nevertheless,
this method should be used with caution.

# Maximum eigenvalues and trace functions {#sec:max-eig}

Our discussion of estimating the spectral norm from a random sample
indicates that there is no straightforward way to construct an unbiased
estimator for the maximum eigenvalue of a PSD matrix. Instead, we turn
to Krylov methods, which repeatedly apply the matrix to appropriate
vectors to extract more information.

The power method and the Lanczos method are two classic algorithms of
this species. Historically, these algorithms have been initialized with
a random vector to ensure that the starting vector has a component in
the direction of a maximum eigenvector. Later, researchers recognized
that the randomness has ancillary benefits. In particular, randomized
algorithms can produce reliable estimates of the maximum eigenvalue,
even when it is not separated from the rest of the spectrum
[@Dix83:Estimating-Extremal; @KW92:Estimating-Largest].

In this section, we summarize theoretical results on the randomized
power method and randomized Krylov methods for computing the maximum
eigenvalue of a psd matrix. These results also have implications for
estimating the minimum eigenvalue of a psd matrix and the spectral norm
of a general matrix. Last, we explain how the Lanczos method leads to
accurate estimates for trace functions.

## Overview

Consider a psd matrix $\bm{\mathsf{A}} \in \mathbb{H}_n$ with
decreasingly ordered eigenvalues
$\lambda_1 \geq \lambda_2 \geq \dots \geq \lambda_n \geq 0$. We are
interested in the problem of estimating the maximum eigenvalue
$\lambda_1$. As in the last two sections, we assume access to
$\bm{\mathsf{A}}$ via the matrix--vector product
$\bm{\mathsf{u}} \mapsto \bm{\mathsf{A}} \bm{\mathsf{u}}$.

In contrast to the methods in
Sections [4](#sec:trace-est){reference-type="ref"
reference="sec:trace-est"} and [5](#sec:schatten-p){reference-type="ref"
reference="sec:schatten-p"}, the algorithms in this section require
sequential applications of the matrix--vector product. In other words,
we now demand nonlinear information about the input matrix
$\bm{\mathsf{A}}$. As a consequence, these algorithms resist
parallelization, and they cannot be used in the one-pass or streaming
environments.

The theoretical treatment in this section also covers the case of
estimating the spectral norm $\Vert \bm{\mathsf{B}} \Vert$ of a
rectangular matrix $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$.
Indeed, we can simply pass to the psd matrix
$\bm{\mathsf{A}} = \bm{\mathsf{B}}^* \bm{\mathsf{B}}$. From an applied
point of view, however, it is important to develop separate algorithms
that avoid squaring the matrix $\bm{\mathsf{B}}$. For brevity, we omit a
discussion about estimating spectral norms; see .

## The randomized power method {#sec:rand-power}

The randomized power method is a simple iterative algorithm for
estimating the maximum eigenvalue.

### Procedure

First, draw a random test vector
$\bm{\mathsf{\omega}} \in \mathbb{F}^n$. Since the maximum eigenvalue is
unitarily invariant, it is most natural to draw the test vector
$\bm{\mathsf{\omega}}$ from a rotationally invariant distribution, such
as
$\bm{\mathsf{\omega}} \sim \textsc{normal}(\bm{\mathsf{0}}, \bm{\mathsf{I}})$.
The power method iteratively constructs the sequence
$$\bm{\mathsf{y}}_0 = \frac{\bm{\mathsf{\omega}}}{\Vert \bm{\mathsf{\omega}} \Vert}
\quad\text{and}\quad
\bm{\mathsf{y}}_q = \frac{\bm{\mathsf{A}} \bm{\mathsf{y}}_{q-1}}{\Vert \bm{\mathsf{A}} \bm{\mathsf{y}}_{q-1} \Vert}
\quad\text{for $q \geq 1$.}$$ At each step, we obtain an eigenvalue
estimate $$\xi_q = \bm{\mathsf{y}}_q^* \bm{\mathsf{A}} \bm{\mathsf{y}}_q
    = \frac{\bm{\mathsf{\omega}}^* \bm{\mathsf{A}}^{2q+1} \bm{\mathsf{\omega}}}{\bm{\mathsf{\omega}}^* \bm{\mathsf{A}}^{2q} \bm{\mathsf{\omega}}}
\quad\text{for $q \geq 0$.}$$ The randomized power method requires
simulation of a single random test vector $\bm{\mathsf{\omega}}$. To
perform $q$ iterations, it takes $q$ sequential matrix--vector products
with $\bm{\mathsf{A}}$ and lower-order arithmetic. It operates with
storage $O(n)$. See
Algorithm [\[alg:rand-power\]](#alg:rand-power){reference-type="ref"
reference="alg:rand-power"} for pseudocode.
::: {#alg:rand-power .algorithm}
**Algorithm. Randomized power method.**

```
**Input:** Psd input matrix $\mtx{A} \in \Sym_n$, maximum number $q$ of iterations,
stopping tolerance $\eps$
**Output:** Estimate $\xi$ of maximum eigenvalue of $\mtx{A}$

function RandomizedPower($\mtx{A}$, $q$, $\eps$):
    $\vct{\omega} = randn(n,1)$
    # Starting vector is random
    $\vct{y}_0 = \vct{\omega} / \norm{\vct{\omega}}$
    for $i = 1, \dots, q$:
        $\vct{y}_{i} = \mtx{A} \vct{y}_{i-1}$
        $\xi_{i-1} = \vct{y}_{i-1}^* \vct{y}_i$
        $\vct{y}_{i} = \vct{y}_{i} / \norm{\vct{y}_i}$
        if $\abs{\xi_{i-1} - \xi_{i-2}} \leq \eps \xi_{i-1}$:
            break
            # [opt] Stopping rule
    return $\xi_{i-1}$
```
:::



:::: algorithm
::: algorithmic
Psd input matrix $\bm{\mathsf{A}} \in \mathbb{H}_n$, maximum number $q$
of iterations, stopping tolerance $\varepsilon$ Estimate $\xi$ of
maximum eigenvalue of $\bm{\mathsf{A}}$

$\bm{\mathsf{\omega}} = \texttt{randn}(n,1)$ Starting vector is random
$\bm{\mathsf{y}}_0 = \bm{\mathsf{\omega}} / \Vert \bm{\mathsf{\omega}} \Vert$

$\bm{\mathsf{y}}_{i} = \bm{\mathsf{A}} \bm{\mathsf{y}}_{i-1}$
$\xi_{i-1} = \bm{\mathsf{y}}_{i-1}^* \bm{\mathsf{y}}_i$
$\bm{\mathsf{y}}_{i} = \bm{\mathsf{y}}_{i} / \Vert \bm{\mathsf{y}}_i \Vert$
**break** Stopping rule

$\xi_{i-1}$
:::
::::

### Analysis

The question is how many iterations $q$ suffice to make $\xi_q$ close to
the maximum eigenvalue $\lambda_1$. More precisely, we aim to control
the relative error $e_q$ in the eigenvalue estimate $\xi_q$:
$$e_q = \frac{\lambda_1 - \xi_q}{\lambda_1}.$$ The error $e_q$ is always
nonnegative because of the Rayleigh theorem. Note that the computed
vector $\bm{\mathsf{y}}_q$ always has a substantial component in the
invariant subspace associated with eigenvalues larger than $\xi_q$, but
it may not be close to any maximum eigenvector, even when
$\xi_q \approx \lambda_1$.

have established several remarkable results about the evolution of the
error.

::: {#thm:rand-power .theorem}
**Theorem 8** (Randomized power method). *Let
$\bm{\mathsf{A}} \in \mathbb{H}_n(\mathbb{R})$ be a real psd matrix.
Draw a real test vector
$\bm{\mathsf{\omega}} \sim \textsc{normal}(\bm{\mathsf{0}}, \bm{\mathsf{I}})$.
After $q$ iterations of the randomized power method, the error $e_q$ in
the maximum eigenvalue estimate $\xi_q$ satisfies $$\begin{equation}
 \label{eqn:rand-power-nogap}
\operatorname{\mathbb{E}}e_q \leq 0.871 \cdot \frac{\log n}{q - 1}
\quad\text{for $q \geq 2$.}
\end{equation}$$ Furthermore, if
$\gamma = (\lambda_1 - \lambda_2) / \lambda_1$ is the relative spectral
gap, then $$\begin{equation}
 \label{eqn:rand-power-gap}
\operatorname{\mathbb{E}}e_q \leq 1.254 \cdot \sqrt{n} \, \gamma \mathrm{e}^{-q \gamma}
\quad\text{for $q \geq 1$.}
\end{equation}$$*
:::

The second
result [\[eqn:rand-power-gap\]](#eqn:rand-power-gap){reference-type="eqref"
reference="eqn:rand-power-gap"} does not appear explicitly in , but it
follows from related ideas [@Tro18:Analysis-Randomized].

### Discussion

Theorem [8](#thm:rand-power){reference-type="ref"
reference="thm:rand-power"} makes two different claims. First, the power
method can exhibit a burn-in period of $q \approx \log n$ iterations
before it produces a nontrivial estimate of the maximum eigenvalue;
after this point, it always decreases the error in proportion to the
number $q$ of iterations. The second claim concerns the situation where
the matrix has a spectral gap $\gamma$ bounded away from zero. In the
latter case, after the burn-in period, the error decreases at each
iteration by a constant factor that depends on the spectral gap. The
burn-in period of $q \approx \log n$ iterations is necessary for *any*
algorithm that estimates the maximum eigenvalue of $\bm{\mathsf{A}}$
from $q$ matrix--vector products with the
matrix [@SER17:Gap-Strict-Saddles].

Whereas classical analyses of the power method depend on the spectral
gap $\gamma$, Theorem [8](#thm:rand-power){reference-type="ref"
reference="thm:rand-power"} comprehends that we can estimate the maximum
eigenvalue even when $\gamma \approx 0$. On the other hand, it is
generally not possible to obtain a reliable estimate of the maximum
eigenvector in this extreme [@LW98:Estimating-Largest].

Theorem [8](#thm:rand-power){reference-type="ref"
reference="thm:rand-power"} can be improved in several respects. First,
we can develop variants where the dimension $n$ of the matrix
$\bm{\mathsf{A}}$ is replaced by its intrinsic
dimension [\[eqn:intdim\]](#eqn:intdim){reference-type="eqref"
reference="eqn:intdim"}, or by smaller quantities that reflect spectral
decay. Second, when the maximum eigenvalue has multiplicity greater than
one, the power method estimates the maximum eigenvalue faster. Third,
the result can be extended to the complex setting. See  for further
discussion.

Although the power method is often deprecated because it converges
slowly, it is numerically stable, and it enjoys the (minimal) storage
cost of $O(n)$.

## Randomized Krylov methods

The power method uses the $q$th power of the matrix to estimate the
maximum eigenvalue. A more sophisticated approach allows *any*
polynomial with degree $q$. Algorithms based on this general technique
are often referred to as *Krylov subspace methods*. The most famous
instantiation is the *Lanczos method*, which is an efficient
implementation of a Krylov subspace method for estimating the
eigenvalues of a self-adjoint matrix.

### Abstract procedure

Draw a random test vector $\bm{\mathsf{\omega}} \in \mathbb{F}^n$. It is
natural to use a rotationally invariant distribution, such as
$\bm{\mathsf{\omega}} \sim \textsc{normal}(\bm{\mathsf{0}}, \bm{\mathsf{I}})$.
For a depth parameter $q \in \mathbb{N}$, a randomized Krylov subspace
method implicitly constructs the subspace $$K_{q+1} %
    := \operatorname{span}\{ \bm{\mathsf{\omega}}, \bm{\mathsf{A}}\bm{\mathsf{\omega}}, \dots, \bm{\mathsf{A}}^q \bm{\mathsf{\omega}} \}.$$
We can estimate the maximum eigenvalue of $\bm{\mathsf{A}}$ as
$$\xi_q = \max_{\bm{\mathsf{u}} \in K_{q+1}} \frac{\bm{\mathsf{u}}^* \bm{\mathsf{A}} \bm{\mathsf{u}}}{\bm{\mathsf{u}}^*\bm{\mathsf{u}}}
    = \max_{\deg p \leq q} \frac{\bm{\mathsf{\omega}}^* \bm{\mathsf{A}} p^2(\bm{\mathsf{A}}) \bm{\mathsf{\omega}}}{\Vert  p(\bm{\mathsf{A}}) \bm{\mathsf{\omega}}  \Vert^2}.$$
The maximum occurs over polynomials $p$ with coefficients in
$\mathbb{F}$ and with degree at most $q$. The notation
$p(\bm{\mathsf{A}})$ refers to the spectral function induced by the
polynomial $p$. We will discuss implementations of Krylov methods below
in Section [6.3.4](#sec:rand-lanczos){reference-type="ref"
reference="sec:rand-lanczos"}.

### Analysis

Since the Krylov subspace is invariant to shifts in the spectrum of
$\bm{\mathsf{A}}$, it is more natural to compute the error relative to
the spectral range of the matrix:
$$f_q = \frac{\lambda_1 - \xi_q}{\lambda_1 - \lambda_n}.$$ The error
$f_q$ is always nonnegative because it is a Rayleigh quotient.

have established striking results for the maximum eigenvalue estimate
obtained via a randomized Krylov subspace method.

::: {#thm:rand-krylov .theorem}
**Theorem 9** (Randomized Krylov method). *Let
$\bm{\mathsf{A}} \in \mathbb{H}_n(\mathbb{R})$ be a real psd matrix.
Draw a real test vector
$\bm{\mathsf{\omega}} \sim \textsc{normal}(\bm{\mathsf{0}}, \bm{\mathsf{I}})$.
After $q$ iterations of the randomized Krylov method, the error $f_q$ in
the maximum eigenvalue estimate $\xi_q$ satisfies $$\begin{equation}
 \label{eqn:rand-krylov-nogap}
\operatorname{\mathbb{E}}f_q \leq 2.575 \cdot \left(\frac{\log n}{q - 1} \right)^2
\quad\text{for $q \geq 4$.}
\end{equation}$$ Furthermore, if
$\gamma = (\lambda_1 - \lambda_2) / \lambda_1$ is the relative spectral
gap, then $$\begin{equation}
 \label{eqn:rand-krylov-gap}
\operatorname{\mathbb{E}}f_q \leq 2.589 \cdot \sqrt{n} \, \mathrm{e}^{-2(q-1) \sqrt{\gamma}}
\quad\text{for $q \geq 1$.}
\end{equation}$$*
:::

The second
result [\[eqn:rand-krylov-gap\]](#eqn:rand-krylov-gap){reference-type="eqref"
reference="eqn:rand-krylov-gap"} is a direct consequence of a more
detailed formula reported in .

### Discussion

Like the randomized power method, a randomized Krylov method can also
exhibit a burn-in period of $q \approx \log n$ steps. Afterwards, the
result [\[eqn:rand-krylov-nogap\]](#eqn:rand-krylov-nogap){reference-type="eqref"
reference="eqn:rand-krylov-nogap"} shows that the error decreases in
proportion to $1/q^2$, which is much faster than the $1/q$ rate achieved
by the power method. Furthermore, the
result [\[eqn:rand-krylov-gap\]](#eqn:rand-krylov-gap){reference-type="eqref"
reference="eqn:rand-krylov-gap"} shows that each iteration decreases the
error by a constant factor $\mathrm{e}^{-2\sqrt{\gamma}}$ where $\gamma$
is the spectral gap. In contrast, the power method only decreases the
error by a constant factor $\mathrm{e}^{-\gamma}$.

Theorem [9](#thm:rand-krylov){reference-type="ref"
reference="thm:rand-krylov"} admits the same kind of refinements as
Theorem [8](#thm:rand-power){reference-type="ref"
reference="thm:rand-power"}. In particular, we can replace the dimension
$n$ with measures that reflect the spectral decay of the input matrix.
See for details.

### Implementing Krylov methods {#sec:rand-lanczos}

What do we have to pay for the superior performance of the randomized
Krylov method? If we only need an estimate of the maximum eigenvalue,
without an associated eigenvector estimate, the cost is almost the same
as for the randomized power method! On the other hand, if we desire the
eigenvector estimate, it is common practice to store a basis for the
Krylov subspace $K_q$. This is a classic example of a time--data
tradeoff in computation.

We present pseudocode for the randomized Lanczos method, which is an
efficient formulation of the Krylov method.
Algorithm [\[alg:rand-lanczos\]](#alg:rand-lanczos){reference-type="ref"
reference="alg:rand-lanczos"} is a direct implementation of the Lanczos
::: {#alg:rand-lanczos .algorithm}
**Algorithm. Randomized Lanczos method (with full reorthogonalization).**

```
**Input:** Psd input matrix $\mtx{A} \in \Sym_n$, maximum number $q$ of iterations
**Output:** Estimate $(\xi, \vct{y})$ for a maximum eigenpair of $\mtx{A}$

function RandomizedLanczos($\mtx{A}$, $q$):
    $q = \min(q, n-1)$
    $\mtx{Q}(:,1) = randn(n, 1)$
    # Starting vector $\vct{\omega}$ is random
    $\mtx{Q}(:,1) = \mtx{Q}(:,1) / \norm{\mtx{Q}(:,1)}$
    for $i = 1, \dots, q$:
        $\mtx{Q}(:, i+1) = \mtx{A} \mtx{Q}(:, i)$
        $\alpha_i = \real(\mtx{Q}(:, i)^* \mtx{Q}(:, i+1))$
        if $i = 1$:
            $\mtx{Q}(:, i+1) = \mtx{Q}(:, i+1) - \alpha_i \mtx{Q}(:, i)$
        else:
            $\mtx{Q}(:, i+1) = \mtx{Q}(:, i+1) - \alpha_i \mtx{Q}(:, i) - \beta_{i-1} \mtx{Q}(:, i-1)$
        
        
        # [opt] Reorthogonalize via double Gram--Schmidt
        $\mtx{Q}(:, i+1) = \mtx{Q}(:, i+1) - \mtx{Q}(:, 1:i) (\mtx{Q}(:, 1:i)^*\mtx{Q}(:, i+1))$
        $\mtx{Q}(:, i+1) = \mtx{Q}(:, i+1) - \mtx{Q}(:, 1:i) (\mtx{Q}(:, 1:i)^*\mtx{Q}(:, i+1))$
        
        $\beta_i = \norm{ \mtx{Q}(:, i+1) }$
        if $\beta_i < \mu \sqrt{n}$:
            break
            # $\mu$ is machine precision
        $\mtx{Q}(:, i+1) = \mtx{Q}(:, i+1) / \beta_i$
    $\mtx{T} = tridiag(\beta(1:i-1), \alpha(1:i), \beta(1:i-1))$
    $[\mtx{V}, \mtx{D}] = eig(\mtx{T})$
    $[\xi, ind] = \min(\diag(\mtx{D}))$
    $\vct{y} = \mtx{Q}(:, 1:i) \, \mtx{V}(:, ind)$
    # [opt] Estimate max eigenvector
```
:::


recursion, but it exhibits complicated performance in floating-point
arithmetic.
Algorithm [\[alg:rand-lanczos\]](#alg:rand-lanczos){reference-type="ref"
reference="alg:rand-lanczos"} includes the option to add full
reorthogonalization; this step removes the numerical shortcomings at a
substantial price in arithmetic (and storage).

If we use the Lanczos method without orthogonalization, then $q$
iterations require $q$ matrix--vector multiplies with $\bm{\mathsf{A}}$
plus $O(qn)$ additional arithmetic. The orthogonalization step adds a
total of $O(q^2 n)$ additional arithmetic. Computing the maximum
eigenvalue (and eigenvector) of the tridiagonal matrix $\bm{\mathsf{T}}$
can be performed with $O(q)$ arithmetic [@GVL13:Matrix-Computations-4ed
Sec. 8.4].

If we do not require the maximum eigenvector, the Lanczos method without
orthogonalization operates with storage $O(n)$. If we need the maximum
eigenvector or we add orthogonalization, the storage cost grows to
$O(qn)$. It is possible to avoid the extra storage by recomputing the
Lanczos vectors, but this approach requires great care. One of the main
thrusts in the literature on Krylov methods is to reduce these storage
costs while maintaining rapid convergence.

::: warning
**Warning 10** (Lanczos method).
*Algorithm [\[alg:rand-lanczos\]](#alg:rand-lanczos){reference-type="ref"
reference="alg:rand-lanczos"} should be used with caution. For a proper
discussion about designing Krylov methods, we recommend the books
[@Par98:Symmetric-Eigenvalue; @BDD+00:Templates-Solution; @GVL13:Matrix-Computations-4ed].
There is also some recent theoretical work on the numerical stability of
Lanczos methods [@MMS18:Stability-Lanczos; @CDST19:Rank-1-Sketch].*
:::

:::: algorithm
::: algorithmic
Psd input matrix $\bm{\mathsf{A}} \in \mathbb{H}_n$, maximum number $q$
of iterations Estimate $(\xi, \bm{\mathsf{y}})$ for a maximum eigenpair
of $\bm{\mathsf{A}}$

$q = \min(q, n-1)$ $\bm{\mathsf{Q}}(:,1) = \texttt{randn}(n, 1)$
Starting vector $\bm{\mathsf{\omega}}$ is random
$\bm{\mathsf{Q}}(:,1) = \bm{\mathsf{Q}}(:,1) / \Vert \bm{\mathsf{Q}}(:,1) \Vert$

$\bm{\mathsf{Q}}(:, i+1) = \bm{\mathsf{A}} \bm{\mathsf{Q}}(:, i)$
$\alpha_i = \operatorname{real}(\bm{\mathsf{Q}}(:, i)^* \bm{\mathsf{Q}}(:, i+1))$

$\bm{\mathsf{Q}}(:, i+1) = \bm{\mathsf{Q}}(:, i+1) - \alpha_i \bm{\mathsf{Q}}(:, i)$
$\bm{\mathsf{Q}}(:, i+1) = \bm{\mathsf{Q}}(:, i+1) - \alpha_i \bm{\mathsf{Q}}(:, i) - \beta_{i-1} \bm{\mathsf{Q}}(:, i-1)$

Reorthogonalize via double Gram--Schmidt

$\bm{\mathsf{Q}}(:, i+1) = \bm{\mathsf{Q}}(:, i+1) - \bm{\mathsf{Q}}(:, 1:i) (\bm{\mathsf{Q}}(:, 1:i)^*\bm{\mathsf{Q}}(:, i+1))$
$\bm{\mathsf{Q}}(:, i+1) = \bm{\mathsf{Q}}(:, i+1) - \bm{\mathsf{Q}}(:, 1:i) (\bm{\mathsf{Q}}(:, 1:i)^*\bm{\mathsf{Q}}(:, i+1))$

$\beta_i = \Vert  \bm{\mathsf{Q}}(:, i+1)  \Vert$

**break** $\mu$ is machine precision

$\bm{\mathsf{Q}}(:, i+1) = \bm{\mathsf{Q}}(:, i+1) / \beta_i$

$\bm{\mathsf{T}} = \texttt{tridiag}(\beta(1:i-1), \alpha(1:i), \beta(1:i-1))$

$[\bm{\mathsf{V}}, \bm{\mathsf{D}}] = \texttt{eig}(\bm{\mathsf{T}})$
$[\xi, \texttt{ind}] = \min(\operatorname{diag}(\bm{\mathsf{D}}))$
$\bm{\mathsf{y}} = \bm{\mathsf{Q}}(:, 1:i) \, \bm{\mathsf{V}}(:, \texttt{ind})$
Estimate max eigenvector
:::
::::

## The minimum eigenvalue

The randomized power method and the randomized Krylov method can be used
to estimate the minimum eigenvalue $\lambda_{n}$ of the psd matrix
$\bm{\mathsf{A}} \in \mathbb{H}_n$.

The first approach is to apply the randomized power method to the
shifted matrix $\nu \bm{\mathsf{I}}- \bm{\mathsf{A}}$, where the shift
is chosen so that $\nu \geq \lambda_1$. In this case, the algorithm
produces an approximation for $\nu - \lambda_n$. Note that the error in
Theorem [8](#thm:rand-power){reference-type="ref"
reference="thm:rand-power"} is relative to $\nu - \lambda_n$, rather
than $\lambda_n$.

The second approach begins with the computation of the Krylov subspace
$K_{q+1}$. Instead of maximizing the Rayleigh quotient over the Krylov
subspace, we minimize it:
$$\zeta_q = \min_{\bm{\mathsf{u}} \in K_{q+1}} \frac{ \bm{\mathsf{u}}^*\bm{\mathsf{A}}\bm{\mathsf{u}} }{\bm{\mathsf{u}}^*\bm{\mathsf{u}}}.$$
This approach directly produces an estimate $\zeta_q$ for $\lambda_n$.
The Krylov subspace is invariant under affine transformations of the
spectrum of $\bm{\mathsf{A}}$, so we can obtain an error bound for
$\zeta_q$ by applying Theorem [9](#thm:rand-krylov){reference-type="ref"
reference="thm:rand-krylov"} formally to
$\lambda_n \bm{\mathsf{I}}- \bm{\mathsf{A}}$.

::: remark
**Remark 11** (Inverses). *If we can apply the matrix inverse
$\bm{\mathsf{A}}^{-1}$ to vectors, then we gain access to a wider class
of algorithms for computing the minimum eigenvalue, including (shifted)
inverse iteration and the Rayleigh quotient iteration. See
[@Par98:Symmetric-Eigenvalue] and [@GVL13:Matrix-Computations-4ed].*
:::

## Block methods

The basic power method and Krylov method can be extended by applying the
iteration simultaneously to a larger number of (random) test vectors.
The resulting algorithms are called *subspace iteration* and *block
Krylov methods*, respectively. Historically, the reason for developing
block methods was to resolve repeated or clustered eigenvalues.

In randomized linear algebra, we discover additional motivations for
developing block methods. When the test vectors are drawn at random,
block methods may converge slightly faster, and they succeed with much
higher probability. On modern computer architectures, the cost of a
block method may be comparable with the cost of a simple vector
iteration, which makes this modification appealing. We will treat this
class of algorithm more thoroughly in
Section [11](#sec:random-rangefinder){reference-type="ref"
reference="sec:random-rangefinder"}, so we postpone a full discussion.
See also [@Tro18:Analysis-Randomized].

## Estimating trace functions {#sec:slq}

Finally, we turn to the problem of estimating the trace of a spectral
function of a psd matrix.

### Overview

Consider a psd matrix $\bm{\mathsf{A}} \in \mathbb{H}_n$ with eigenpairs
$(\lambda_j, \bm{\mathsf{u}}_j)$ for $j = 1, \dots, n$. Let
$f : \mathbb{R}_+ \to \mathbb{R}$ be a function, and suppose that we
wish to approximate
$$\operatorname{trace}f(\bm{\mathsf{A}}) = \sum\nolimits_{j=1}^n f(\lambda_j).$$
We outline an incredible approach to this problem, called *stochastic
Lanczos quadrature* (SLQ), that marries the randomized trace estimator
(Section [4](#sec:trace-est){reference-type="ref"
reference="sec:trace-est"}) to the Lanczos iteration
(Algorithm [\[alg:rand-lanczos\]](#alg:rand-lanczos){reference-type="ref"
reference="alg:rand-lanczos"}).

This algorithm was devised by [@GM94:Matrices-Moments]. Our presentation
is based on  and . provide a more complete treatment, and give a
theoretical discussion about stability.

Related ideas can be used to estimate the trace of a spectral function
of a rectangular matrix; that is, the sum of a function of each singular
value of the matrix. For brevity, we omit all details on the rectangular
case.

### Examples

Computing the trace of a spectral function is a ubiquitous problem with
a huge number of applications. Let us mention some of the key examples.

1.  For $f(\lambda) = \lambda^{-1}$, the resulting trace function is the
    trace of the matrix inverse. This computation arises in electronic
    structure calculations.

2.  For $f(t) = \log t$, the associated trace function is the
    log-determinant. This computation arises in Gaussian process
    regression.

3.  For $f(t) = t^p$ with $p \geq 1$, the trace function is the $p$th
    power of the Schatten $p$-norm. SLQ offers a more powerful
    alternative to the methods in
    Section [5](#sec:schatten-p){reference-type="ref"
    reference="sec:schatten-p"}.

We refer to  for additional discussion.

### Procedure

Let us summarize the mathematical ideas that lead to SLQ. As usual, we
draw an isotropic random vector $\bm{\mathsf{\omega}} \in \mathbb{F}^n$.
Then the random variable
$$X = \bm{\mathsf{\omega}}^* f(\bm{\mathsf{A}}) \bm{\mathsf{\omega}}
\quad\text{satisfies}\quad
\operatorname{\mathbb{E}}X = \operatorname{trace}f(\bm{\mathsf{A}}).$$
Using the spectral resolution of $\bm{\mathsf{A}}$, we can rewrite $X$
in the form
$$X = \sum\nolimits_{j=1}^n f(\lambda_j) \, \vert \bm{\mathsf{u}}_j^* \bm{\mathsf{\omega}} \vert^2
    = \int_{\mathbb{R}_+} f(\lambda) \, {\nu}(\mathrm{d}{\lambda})$$ for
an appropriate measure $\nu$ on $\mathbb{R}_+$ that depends on
$\bm{\mathsf{A}}$ and $\bm{\mathsf{\omega}}$. Although we cannot
generally compute the integral directly, we can approximate it by using
a numerical quadrature rule:
$$X \approx \sum\nolimits_{\ell=1}^{q+1} \tau_\ell^2 f(\theta_\ell) =: Z.$$
What is truly amazing is that the weights $\tau_\ell^2$ and the nodes
$\theta_\ell$ for the quadrature rule can be extracted from the
tridiagonal matrix $\bm{\mathsf{T}} \in \mathbb{H}_{q+1}$ produced by
$q$ iterations of the Lanczos iteration with starting vector
$\bm{\mathsf{\omega}}$. This point is not obvious, but a full
explanation exceeds our scope.

SLQ approximates the trace function by averaging independent copies of
the simple approximation:
$$\operatorname{trace}f(\bm{\mathsf{A}}) \approx \frac{1}{k} \sum\nolimits_{i=1}^k Z_i
\quad\text{where}\quad
\text{$Z_i \sim Z$ are iid.}$$ The analysis of the SLQ approximation
requires heavy machinery from approximation theory. See  and  for more
details.

Algorithm [\[alg:slq\]](#alg:slq){reference-type="ref"
reference="alg:slq"} contains pseudocode for SLQ. The dominant cost is
::: {#alg:slq .algorithm}
**Algorithm. Stochastic Lanczos quadrature.**

```
**Input:** Psd input matrix $\mtx{A} \in \Sym_n$, function $f$, number $k$ of samples, number $q$ of Lanczos iterations
**Output:** Estimate $\bar{Z}_k$ for $\trace f(\mtx{A})$.

function StochasticLanczosQuadrature($\mtx{A}$, $f$, $k$, $q$):
    for $i = 1, \dots, k$:
        # Extract $k$ independent samples
        Draw a random isotropic vector $\vct{\omega}_i \in \F^n$
        Form $\mtx{T} = RandomizedLanczos(\mtx{A}, \vct{\omega}_i, q)$
        
        # Apply $q$ steps of Lanczos with starting vector $\vct{\omega}_i$
        $[\mtx{V}, \mtx{\Theta}] = eig(\mtx{T})$
        # Tridiagonal eigenproblem
        Extract nodes $\mtx{\Theta} = \diag(\theta_1, \dots, \theta_{q+1})$
        Extract weights $\vct{\delta}_1^* \mtx{V} = (\tau_1, \dots, \tau_{q+1})$
        Form the approximation $Z_i = \sum_{\ell=1}^{q+1} \tau_\ell^2 f(\theta_\ell)$
    Return $\bar{Z}_k = k^{-1} \sum_{i=1}^k Z_i$
```
:::


$O(kq)$ matrix--vector multiplies with $\bm{\mathsf{A}}$, plus
$O(k q^2)$ additional arithmetic. We recommend using structured random
test vectors to reduce the variance of the resulting approximation. The
storage cost is $O(qn)$ numbers.

:::: algorithm
::: algorithmic
Psd input matrix $\bm{\mathsf{A}} \in \mathbb{H}_n$, function $f$,
number $k$ of samples, number $q$ of Lanczos iterations Estimate
$\bar{Z}_k$ for $\operatorname{trace}f(\bm{\mathsf{A}})$.

Extract $k$ independent samples Draw a random isotropic vector
$\bm{\mathsf{\omega}}_i \in \mathbb{F}^n$

Form
$\bm{\mathsf{T}} = \textsc{RandomizedLanczos}(\bm{\mathsf{A}}, \bm{\mathsf{\omega}}_i, q)$
Apply $q$ steps of Lanczos with starting vector $\bm{\mathsf{\omega}}_i$

$[\bm{\mathsf{V}}, \bm{\mathsf{\Theta}}] = \texttt{eig}(\bm{\mathsf{T}})$
Tridiagonal eigenproblem

Extract nodes
$\bm{\mathsf{\Theta}} = \operatorname{diag}(\theta_1, \dots, \theta_{q+1})$
Extract weights
$\bm{\mathsf{\delta}}_1^* \bm{\mathsf{V}} = (\tau_1, \dots, \tau_{q+1})$

Form the approximation
$Z_i = \sum_{\ell=1}^{q+1} \tau_\ell^2 f(\theta_\ell)$

Return $\bar{Z}_k = k^{-1} \sum_{i=1}^k Z_i$
:::
::::

# Matrix approximation by sampling {#sec:matrix-mc}

In Section [4](#sec:trace-est){reference-type="ref"
reference="sec:trace-est"}, we have seen that it is easy to form an
unbiased estimator for the trace of a matrix. By averaging multiple
copies of the simple estimator, we can improve the quality of the
estimate. The same idea applies in the context of matrix approximation.
In this setting, the goal is to produce a matrix approximation that has
more "structure" than the target matrix. The basic approach is to find a
structured unbiased estimator for the matrix and to average multiple
copies of the simple estimator to improve the approximation quality.

This section outlines two instances of this methodology. First, we
develop a toy algorithm for approximate multiplication of matrices.
Second, we show how to approximate a dense graph Laplacian by a sparse
graph Laplacian; this construction plays a role in the
[SparseCholesky]{.smallcaps} algorithm presented in
Section [18](#sec:sparse-cholesky){reference-type="ref"
reference="sec:sparse-cholesky"}.
Section [19](#sec:kernel){reference-type="ref" reference="sec:kernel"}
contains another example, the method of random features in kernel
learning.

The material in this section is summarized from the treatments of matrix
concentration in
[@Tro15:Introduction-Matrix; @Tro19:Matrix-Concentration-LN].

## Empirical approximation

We begin with a high-level discussion of the method of empirical
approximation of a matrix. The examples in this section are all
instances of this general idea.

Let $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$ be a target matrix
that we wish to approximate by a more "structured" matrix. Suppose that
we can express the matrix $\bm{\mathsf{B}}$ as a sum of "simple"
matrices: $$\begin{equation}
 \label{eqn:B-simple-sum}
\bm{\mathsf{B}} = \sum\nolimits_{i=1}^I \bm{\mathsf{B}}_i.
\end{equation}$$ In the cases we will study, the summands
$\bm{\mathsf{B}}_i$ will be sparse or low-rank.

Next, consider a probability distribution $\{ p_i : i = 1, \dots, I \}$
over the indices in the
sum [\[eqn:B-simple-sum\]](#eqn:B-simple-sum){reference-type="eqref"
reference="eqn:B-simple-sum"}. For now, we treat this distribution as
given. Construct a random matrix
$\bm{\mathsf{X}} \in \mathbb{F}^{m \times n}$ that takes values
$$\bm{\mathsf{X}} = p_i^{-1} \bm{\mathsf{B}}_i
\quad\text{with probability $p_i$ for each $i = 1, \dots, I$.}$$
(Enforce the convention that $0/0 = 0$.) It is clear that
$\bm{\mathsf{X}}$ is an unbiased estimator for $\bm{\mathsf{B}}$:
$$\operatorname{\mathbb{E}}\bm{\mathsf{X}} = \sum\nolimits_{i=1}^I (p_i^{-1} \bm{\mathsf{B}}_i) p_i
    = \bm{\mathsf{B}}.$$ The random matrix $\bm{\mathsf{X}}$ enjoys the
same kind of structure as the summands $\bm{\mathsf{B}}_i$. On the other
hand, a single draw of the matrix $\bm{\mathsf{X}}$ is rarely a good
approximation of the matrix $\bm{\mathsf{B}}$.

To obtain a better estimator for $\bm{\mathsf{B}}$, we average
independent copies of the initial estimator:
$$\bar{\bm{\mathsf{X}}}_k = \frac{1}{k} \sum\nolimits_{i=1}^k \bm{\mathsf{X}}_i
\quad\text{where $\bm{\mathsf{X}}_i \sim \bm{\mathsf{X}}$ are iid.}$$ By
linearity of expectation, $\bar{\bm{\mathsf{X}}}_k$ is also unbiased:
$$\operatorname{\mathbb{E}}\bar{\bm{\mathsf{X}}}_k = \bm{\mathsf{B}}.$$
If the parameter $k$ remains small, then $\bar{\bm{\mathsf{X}}}_k$
inherits some of the structure from the $\bm{\mathsf{B}}_i$. The
question is how many samples $k$ we need to ensure that
$\bar{\bm{\mathsf{X}}}_k$ also approximates $\bm{\mathsf{B}}$ well.

::: remark
**Remark 12** (History). *Empirical approximation was first developed by
Maurey in unpublished work from the late 1970s on the geometry of Banach
spaces. The idea was first broadcast by , in a paper on approximation
theory. Applications to randomized linear algebra were proposed by and
by
[@DKM06:Fast-Monte-Carlo-I; @DKM06:Fast-Monte-Carlo-II; @DKM06:Fast-Monte-Carlo-III].
More refined analyses of empirical matrix approximation were obtained
in [@RV07:Sampling-Large; @Tro15:Introduction-Matrix]. Many other papers
consider specific applications of the same methodology.*
:::

## The matrix Bernstein inequality

The main tool for analyzing the sample average estimator
$\bar{\bm{\mathsf{X}}}_k$ from the last section is a variant of the
matrix Bernstein inequality. The following result is drawn from .

::: {#thm:mtx-sampling .theorem}
**Theorem 13** (Matrix Monte Carlo). *Let
$\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$ be a fixed matrix.
Construct a random matrix $\bm{\mathsf{X}} \in \mathbb{F}^{m \times n}$
that satisfies
$$\operatorname{\mathbb{E}}\bm{\mathsf{X}} = \bm{\mathsf{B}}
\quad\text{and}\quad
\Vert \bm{\mathsf{X}} \Vert \leq R.$$ Define the per-sample second
moment:
$$v(\bm{\mathsf{X}}) := \max\{ \Vert  \operatorname{\mathbb{E}}[\bm{\mathsf{X}} \bm{\mathsf{X}}^*]  \Vert, \Vert \operatorname{\mathbb{E}}[\bm{\mathsf{X}}^*\bm{\mathsf{X}}] \Vert \}.$$
Form the matrix sampling estimator
$$\bar{\bm{\mathsf{X}}}_k = \frac{1}{k} \sum\nolimits_{i=1}^k \bm{\mathsf{X}}_i
\quad\text{where $\bm{\mathsf{X}}_i \sim \bm{\mathsf{X}}$ are iid.}$$
Then
$$\operatorname{\mathbb{E}}\Vert  \bar{\bm{\mathsf{X}}}_k - \bm{\mathsf{B}}  \Vert \leq \sqrt{\frac{2 v(\bm{\mathsf{X}}) \log(m+n)}{k}} + \frac{2R \log(m+n)}{3k}.$$
Furthermore, for all $t \geq 0$,
$$\mathbb{P}\left\{  \Vert  \bar{\bm{\mathsf{X}}}_k - \bm{\mathsf{B}}  \Vert \geq t  \right\} \leq (m + n) \exp\left( \frac{-kt^2/2}{v(\bm{\mathsf{X}}) + 2Rt/3} \right).$$*
:::

To explain the meaning of this result, let us determine how large $k$
should be to ensure that the expected approximation error lies below a
positive threshold $\varepsilon$. The bound
$$k \geq \max\left\{ \frac{2v(\bm{\mathsf{X}})\log(m+n)}{\varepsilon^2},\ \frac{2R \log(m+n)}{3\varepsilon} \right\}$$
implies that
$\operatorname{\mathbb{E}}\Vert \bar{\bm{\mathsf{X}}}_k - \bm{\mathsf{B}}  \Vert \leq \varepsilon+\varepsilon^{2}$.
In other words, the number $k$ of samples should be proportional to the
larger of the second moment $v(\bm{\mathsf{X}})$ and the upper bound
$R$.

This fact points toward a disappointing feature of empirical matrix
approximation: To make $\varepsilon$ small, the number $k$ of samples
must increase with $\varepsilon^{-2}$ and often with $\log (m + n)$ as
well. This phenomenon is an unavoidable consequence of the central limit
theorem and the geometry induced by the spectral norm. It means that
matrix sampling estimators are not suitable for achieving high-precision
approximations. See [@Tro15:Introduction-Matrix Sec. 6.2.3] for further
discussion. The logarithmic terms are also necessary in the worst case.

::: remark
**Remark 14** (History). *The matrix Bernstein inequality has a long
history, outlined in . The earliest related results concern uniform
smoothness estimates [@TJ74:Moduli-Smoothness] and Khintchine
inequalities [@LP86:Inegalites-Khintchine] for the Schatten classes. A
first application to statistics appears in . Modern approaches are based
on the matrix Laplace transform
method [@AW02:Strong-Converse; @Oli09:Concentration-Adjacency; @Tro12:User-Friendly]
or on the method of exchangeable
pairs [@MJCFT14:Matrix-Concentration; @Tro16:Expected-Norm].*
:::

## Approximate matrix multiplication {#sec:approx-mtx-mult}

A first application of empirical matrix approximation is to approximate
the product of two matrices: $\bm{\mathsf{M}} = \bm{\mathsf{BC}}$ where
$\bm{\mathsf{B}} \in \mathbb{F}^{m \times I}$ and
$\bm{\mathsf{C}} \in \mathbb{F}^{I \times n}$. Computing the product by
direct matrix--matrix multiplication requires $O(mnI)$ arithmetic
operations. When the inner dimension $I$ is very large as compared with
$m$ and $n$, we might try to reduce the cost by sampling.

In our own work, we have not encountered situations where approximate
matrix multiplication is practical because the quality of the
approximation is very low. Nevertheless, the theory serves a dual
purpose as the foundation for subspace embedding by discrete sampling
(Section [9.6](#sec:coord-embed){reference-type="ref"
reference="sec:coord-embed"}).

### Matrix multiplication by sampling

To simplify the analysis, let us pre-scale the matrices
$\bm{\mathsf{B}}$ and $\bm{\mathsf{C}}$ so that each one has spectral
norm equal to one:
$$\Vert \bm{\mathsf{B}} \Vert = \Vert \bm{\mathsf{C}} \Vert = 1.$$ This
step can be performed efficiently using the spectral norm estimators
outlined in Section [6](#sec:max-eig){reference-type="ref"
reference="sec:max-eig"}. This normalization will remain in force for
the rest of Section [7.3](#sec:approx-mtx-mult){reference-type="ref"
reference="sec:approx-mtx-mult"}.

Observe that the matrix--matrix product satisfies
$$\bm{\mathsf{BC}} = \sum\nolimits_{i=1}^I (\bm{\mathsf{B}})_{:i} (\bm{\mathsf{C}})_{i:},$$
where $(\bm{\mathsf{B}})_{:i}$ is the $i$th column of $\bm{\mathsf{B}}$
and $(\bm{\mathsf{C}})_{i:}$ is the $i$th row of $\bm{\mathsf{C}}$. This
expression motivates us to approximate the product by sampling terms at
random from the sum.

Let $\{ p_i : i = 1, \dots, I \}$ be a sampling distribution, to be
specified later. Form a random rank-one unbiased estimator for the
product:
$$\bm{\mathsf{X}} = p_i^{-1} \cdot (\bm{\mathsf{B}})_{:i} (\bm{\mathsf{C}})_{i:}
\quad\text{with probability $p_i$ for each $i = 1, \dots, I$.}$$ By
construction,
$\operatorname{\mathbb{E}}\bm{\mathsf{X}} = \bm{\mathsf{BC}}$. We can
average $k$ independent copies of $\bm{\mathsf{X}}$ to obtain a better
approximation $\bar{\bm{\mathsf{X}}}_k$ to the matrix product. The cost
of computing the estimator $\bar{\bm{\mathsf{X}}}_k$ of the matrix
product explicitly is only $O(mnk)$ operations, so it is more efficient
than the full multiplication when $k \ll I$.

Theorem [13](#thm:mtx-sampling){reference-type="ref"
reference="thm:mtx-sampling"} yields an easy analysis of this approach.
The conclusion depends heavily on the choice of sampling distribution.
Regardless, we can never expect to attain a high-accuracy approximation
of the product by sampling because the number $k$ of samples must scale
proportionally with the inverse square $\varepsilon^{-2}$ of the
accuracy parameter $\varepsilon$.

### Uniform sampling {#sec:mtx-mult-unif}

The easiest way to approximate matrix multiplication is to choose
uniform sampling probabilities: $p_i = 1/I$ for each $i = 1, \dots, I$.
To analyze this case, we introduce the *coherence* parameter:
$$\mu(\bm{\mathsf{B}}) := I \cdot \max_{i=1, \dots, I} \Vert (\bm{\mathsf{B}})_{:i} \Vert^2.
%$$ Up to scaling, this is the maximum squared norm of a column of
$\bm{\mathsf{B}}$. Since $\bm{\mathsf{B}} \in \mathbb{F}^{m \times I}$
and $\Vert \bm{\mathsf{B}} \Vert=1$, the coherence lies in the range
$[m, I]$. The difficulty of approximating matrix multiplication by
uniform sampling increases with the coherence of $\bm{\mathsf{B}}$ and
$\bm{\mathsf{C}}^*$.

To invoke Theorem [13](#thm:mtx-sampling){reference-type="ref"
reference="thm:mtx-sampling"}, observe that the per-sample second moment
and the spectral norm of the estimator $\bm{\mathsf{X}}$ satisfy the
bound
$$\max\{ v(\bm{\mathsf{X}}), \Vert \bm{\mathsf{X}} \Vert \} \leq \max\{\mu(\bm{\mathsf{B}}), \mu(\bm{\mathsf{C}}^*)\}.$$
Let $\varepsilon\in (0, 1]$ be an accuracy parameter. If
$$k \geq 2 \varepsilon^{-2} \max\{ \mu(\bm{\mathsf{B}}), \mu(\bm{\mathsf{C}}^*) \} \log(m + n),$$
then
$$\operatorname{\mathbb{E}}\Vert  \bar{\bm{\mathsf{X}}}_k - \bm{\mathsf{BC}}  \Vert \leq 2\varepsilon.$$
In other words, the number $k$ of rank-one factors we need to obtain a
relative approximation of the matrix product is proportional to the
maximum coherence of $\bm{\mathsf{B}}$ and $\bm{\mathsf{C}}^*$. In the
best scenario, the number of samples is proportional to
$\max\{m,n\} \log(m+n)$; in the worst case, the sample complexity can be
as large as $I$.

### Importance sampling {#sec:mtx-mult-import}

If the norms of the columns of the factors $\bm{\mathsf{B}}$ and
$\bm{\mathsf{C}}^*$ vary wildly, we may need to use importance sampling
to obtain a nontrivial approximation bound. Define the sampling
probabilities
$$p_i = \frac{\Vert (\bm{\mathsf{B}})_{:i} \Vert^2 + \Vert (\bm{\mathsf{C}})_{:i} \Vert^2}{\Vert \bm{\mathsf{B}} \Vert_{\mathrm{F}}^2 + \Vert \bm{\mathsf{C}} \Vert_{\mathrm{F}}^2}
\quad\text{for $i = 1, \dots, I$.}$$ These probabilities are designed to
balance terms arising from
Theorem [13](#thm:mtx-sampling){reference-type="ref"
reference="thm:mtx-sampling"}. In most cases, we compute the sampling
distribution by directly evaluating the formula, at the cost of
$O((m + n)I)$ operations.

With the importance sampling distribution, the per-sample second moment
and the spectral norm of the estimator $\bm{\mathsf{X}}$ satisfy
$$\max\{ v(\bm{\mathsf{X}}), \Vert \bm{\mathsf{X}} \Vert \} \leq \frac{1}{2} (\operatorname{srank}(\bm{\mathsf{B}}) + \operatorname{srank}(\bm{\mathsf{C}})).
%
%$$ The stable rank is defined
in [\[eqn:stable-rank\]](#eqn:stable-rank){reference-type="eqref"
reference="eqn:stable-rank"}. Let $\varepsilon\in (0, 1]$ be an accuracy
parameter. If
$$k \geq \varepsilon^{-2} \left( \operatorname{srank}(\bm{\mathsf{B}}) + \operatorname{srank}(\bm{\mathsf{C}}) \right) \log(m + n),$$
then
$$\operatorname{\mathbb{E}}\Vert  \bar{\bm{\mathsf{X}}}_k - \bm{\mathsf{BC}}  \Vert \leq 2 \varepsilon.$$
In other words, the number of rank-one factors we need to approximate
the matrix product by importance sampling is $\log(m+n)$ times the total
stable rank of the matrices.

The sample complexity bound for importance sampling always improves over
the bound for uniform sampling. Indeed,
$\operatorname{srank}(\bm{\mathsf{B}}) \leq \mu(\bm{\mathsf{B}})$ under
the normalization $\Vert \bm{\mathsf{B}} \Vert=1$. There are also many
situations where the stable rank is smaller than either dimension of the
matrix. These are the cases where one might consider using approximate
matrix multiplication by importance sampling.

### History

Randomized matrix multiplication was proposed in . It is implicit in ,
while give a more explicit treatment. The analysis here has its origins
in , and the detailed presentation is adapted from . Another interesting
approach to randomized matrix multiplication appears in
[@Pag11:Compressed-Matrix]. See  for further references.

## Approximating a graph by a sparse graph {#sec:graph-sparsification}

As a second application of empirical matrix approximation, we will show
how to take a dense graph and find a sparse graph that serves as a
proxy. This procedure operates by replacing the Laplacian matrix of the
graph with a sparse Laplacian matrix. Beyond its intrinsic interest as a
fact about graphs, this technique plays a central role in randomized
solvers for Laplacian linear systems
(Section [18](#sec:sparse-cholesky){reference-type="ref"
reference="sec:sparse-cholesky"}).

### Graphs and Laplacians

We will consider weighted, loop-free, undirected graphs on the vertex
set $V = \{1, \dots, n\}$. We can specify a weighted graph $G$ on $V$ by
means of a nonnegative weight function
$w : V \times V \to \mathbb{R}_+$. The graph is *loop-free* when
$w_{ii} = 0$ for each vertex $i \in V$. The graph is *undirected* if and
only if $w_{ij} = w_{ji}$ for each pair $(i, j)$. The *sparsity* of an
undirected graph is the number of strictly positive weights $w_{ij}$
with $i \leq j$.

Alternatively, we can work with graph Laplacians. The *elementary
Laplacian* $\bm{\mathsf{\Delta}}_{ij}$ on a vertex pair
$(i, j) \in V \times V$ is the psd matrix
$$\bm{\mathsf{\Delta}}_{ij} = (\bm{\mathsf{\delta}}_{i} - \bm{\mathsf{\delta}}_j)(\bm{\mathsf{\delta}}_i - \bm{\mathsf{\delta}}_j)^* \in \mathbb{H}_n(\mathbb{R}).$$
Let $w$ be the weight function of a loop-free, undirected graph $G$. The
Laplacian associated with the graph $G$ is the matrix $$\begin{equation}
 \label{eqn:graph-laplacian}
\bm{\mathsf{L}}_{G} = \sum\nolimits_{1 \leq i < j \leq n} w_{ij} \bm{\mathsf{\Delta}}_{ij} \in \mathbb{H}_n(\mathbb{R}).
\end{equation}$$ The Laplacian $\bm{\mathsf{L}}_G$ is psd because it is
a nonnegative linear combination of psd matrices.

The Laplacian of a graph is analogous to the Laplacian differential
operator. You can think about $\bm{\mathsf{L}}_G$ as an analog of the
heat kernel that models the diffusion of a particle on the graph. The
Poisson problem $\bm{\mathsf{L}}_G \bm{\mathsf{x}} = \bm{\mathsf{f}}$
serves as a primitive for answering a wide range of questions involving
undirected graphs [@Ten10:Laplacian-Paradigm]. Applications include
clustering and partitioning data, studying random walks on graphs, and
solving finite-element discretizations of elliptic PDEs.

### Spectral approximation

Fix a parameter $\varepsilon\in (0, 1)$. We say that a graph $H$ is an
$\varepsilon$-spectral approximation of a graph $G$ if their Laplacians
are comparable in the psd order: $$\begin{equation}
 \label{eqn:spectral-approximation}
(1 - \varepsilon) \, \bm{\mathsf{L}}_G \preccurlyeq\bm{\mathsf{L}}_H \preccurlyeq(1+\varepsilon) \, \bm{\mathsf{L}}_G.
\end{equation}$$ The
relation [\[eqn:spectral-approximation\]](#eqn:spectral-approximation){reference-type="eqref"
reference="eqn:spectral-approximation"} ensures that the graph $H$ and
the graph $G$ are close cousins.

In particular,
under [\[eqn:spectral-approximation\]](#eqn:spectral-approximation){reference-type="eqref"
reference="eqn:spectral-approximation"}, the matrix $\bm{\mathsf{L}}_H$
serves as an excellent preconditioner for the Laplacian
$\bm{\mathsf{L}}_G$. In other words, if we can easily solve (consistent)
linear systems of the form
$\bm{\mathsf{L}}_H \bm{\mathsf{y}} = \bm{\mathsf{b}}$, then we can just
as easily solve the Poisson problem
$\bm{\mathsf{L}}_G \bm{\mathsf{x}} = \bm{\mathsf{f}}$.

For an arbitrary input graph $G$, we will demonstrate that there is a
*sparse* graph $H$ that is a good spectral approximation of $G$. For
several reasons, this construction does not immediately lead to
effective methods for designing preconditioners. Nevertheless, related
ideas have resulted in practical, fast solvers for Poisson problems on
undirected graphs. We will make this connection in
Section [18](#sec:sparse-cholesky){reference-type="ref"
reference="sec:sparse-cholesky"}.

### The normalizing map

It is convenient to present a few more concepts from spectral graph
theory. Let us introduce the *normalizing map*
$$\bm{\mathsf{K}}_G(\bm{\mathsf{A}}) := (\bm{\mathsf{L}}_G^{\dagger})^{1/2} \bm{\mathsf{A}} (\bm{\mathsf{L}}_G^{\dagger})^{1/2}
    \quad\text{for $\bm{\mathsf{A}} \in \mathbb{H}_n(\mathbb{R})$.}$$ As
usual, $(\cdot)^\dagger$ is the pseudoinverse and $(\cdot)^{1/2}$ is the
psd square root of a psd matrix. We use the normalizing map to compare
graph Laplacians. Indeed, $$\begin{equation}
 \label{eqn:graph-comparison}
\Vert  \bm{\mathsf{K}}_G(\bm{\mathsf{L}}_H - \bm{\mathsf{L}}_G)  \Vert \leq \varepsilon
\end{equation}$$ implies that the graph $H$ is an $\varepsilon$-spectral
approximation of the graph $G$. This claim follows easily
from [\[eqn:spectral-approximation\]](#eqn:spectral-approximation){reference-type="eqref"
reference="eqn:spectral-approximation"}.

### Effective resistance

Next, define the *effective resistance* $\varrho_{ij}$ of a vertex pair
$(i, j)$ in the graph $G$ as $$\begin{equation}
 \label{eqn:effective-resistance}
\varrho_{ij} %
    := \operatorname{trace}[ \bm{\mathsf{K}}_G(\bm{\mathsf{\Delta}}_{ij}) ] \geq 0
    \quad\text{for each $1 \leq i < j \leq n$.}
\end{equation}$$ We can compute the family of effective resistances in
time $O(n^3)$ by means of a Cholesky factorization of
$\bm{\mathsf{L}}_G$, although faster algorithms are now available
[@Kyn17:Approximate-Gaussian].

To understand the terminology, let us regard the graph $G$ as an
electrical network where $w_{ij}$ is the conductivity of the wire
connecting the vertex pair $(i, j)$. The effective resistance
$\varrho_{ij}$ is the resistance of the entire electrical network $G$
against passing a unit of current from vertex $i$ to vertex $j$.

### Sparsification by sampling

We are now prepared to construct a sparse approximation of a loop-free,
undirected graph $G$ specified by the weight function $w$. The
representation [\[eqn:graph-laplacian\]](#eqn:graph-laplacian){reference-type="eqref"
reference="eqn:graph-laplacian"} of the graph Laplacian immediately
suggests that we can apply the empirical approximation paradigm. To do
so, we must design an appropriate sampling distribution.

Define the sampling probabilities
$$p_{ij} = \frac{w_{ij} \, \varrho_{ij}}{\operatorname{rank}(\bm{\mathsf{L}}_G)}
\quad\text{for each $1 \leq i < j \leq n$.}$$ It is clear that
$p_{ij} \geq 0$. With some basic matrix algebra, we can also confirm
that $\sum_{i < j} p_{ij} = 1$.

Following our standard approach, we construct the random matrix
$$\bm{\mathsf{X}} = \frac{w_{ij}}{p_{ij}} \, \bm{\mathsf{\Delta}}_{ij}
\quad\text{with probability $p_{ij}$ for each $i < j$.}$$ Next, we
average $k$ independent copies of the estimator:
$$\bar{\bm{\mathsf{X}}}_k = \frac{1}{k} \sum\nolimits_{i=1}^k \bm{\mathsf{X}}_i
\quad\text{where $\bm{\mathsf{X}}_i \sim \bm{\mathsf{X}}$ are iid.}$$ By
construction, the estimator is unbiased:
$\operatorname{\mathbb{E}}\bar{\bm{\mathsf{X}}}_k = \bm{\mathsf{L}}_G$.
Moreover, $\bar{\bm{\mathsf{X}}}_k$ is itself the graph Laplacian of a
(random) graph $H$ with sparsity at most $k$. We just need to determine
the number $k$ of samples that are sufficient to make $H$ a good
spectral approximation of $G$.

Given the effective resistances, computation of the sampling
probabilities requires $O(s)$ time, where $s$ is the sparsity of the
graph. The cost of sampling $k$ copies of $\bm{\mathsf{X}}$ is
$\tilde{O}(s + k)$. It is natural to represent $\bar{\bm{\mathsf{X}}}_k$
using a sparse matrix data structure, with storage cost $O(k \log n)$.

### Analysis

The analysis involves a small twist. In view
of [\[eqn:graph-comparison\]](#eqn:graph-comparison){reference-type="eqref"
reference="eqn:graph-comparison"}, we need to demonstrate that
$\bm{\mathsf{K}}_G(\bm{\mathsf{L}}_H) \approx \bm{\mathsf{K}}_G(\bm{\mathsf{L}}_G)$.
Therefore, instead of considering the random matrix
$\bar{\bm{\mathsf{X}}}_k$, we pass to the random matrix
$\bm{\mathsf{K}}_G(\bar{\bm{\mathsf{X}}}_k)$, which is an unbiased
estimator for $\bm{\mathsf{K}}_G(\bm{\mathsf{L}}_G)$.

Note that the random matrix $\bm{\mathsf{K}}_G(\bm{\mathsf{X}})$
satisfies
$$\max\{ v( \bm{\mathsf{K}}_G(\bm{\mathsf{X}}) ), \Vert  \bm{\mathsf{K}}_G(\bm{\mathsf{X}})  \Vert_2 \} \leq \operatorname{rank}(\bm{\mathsf{L}}_G) < n.$$
Suppose that we choose $$k \geq 3 \varepsilon^{-2} n \log(2n)
\quad\text{where $\varepsilon\in (0, 1)$.}$$ Then
Theorem [13](#thm:mtx-sampling){reference-type="ref"
reference="thm:mtx-sampling"} implies
$$\operatorname{\mathbb{E}}\Vert  \bm{\mathsf{K}}_G( \bar{\bm{\mathsf{X}}}_k - \bm{\mathsf{L}}_G )  \Vert_2 \leq \varepsilon.$$
We conclude that the random graph $H$ with Laplacian
$\bm{\mathsf{L}}_H = \bar{\bm{\mathsf{X}}}_k$ has at most $k$ nonzero
weights, and it is an $\varepsilon$-spectral approximation to the graph
$G$.

In other words, every graph on $n$ vertices has a $(1/2)$-spectral
approximation with at most $12 n \log n$ nonzero weights. Modulo the
precise constant, this is the tightest result that can be obtained if we
form the graph $H$ via random sampling.

### History

The idea of approximating a matrix in the spectral norm by means of
random sampling of entries was proposed by and . This work initiated a
line of literature on matrix sparsification in the randomized NLA
community; see for more references. Let us emphasize that these general
approaches achieve much weaker approximation guarantees
than [\[eqn:spectral-approximation\]](#eqn:spectral-approximation){reference-type="eqref"
reference="eqn:spectral-approximation"}.

The idea of sparsifying a graph by randomly sampling edges in proportion
to the effective resistances was developed by ; the analysis above is
drawn from . A deterministic method for graph sparsification, with
superior guarantees, appears in .

# Randomized embeddings {#sec:gauss}

One of the core tools in randomized linear algebra is *randomized linear
embedding*, often assigned the misnomer *random projection*. The
application of randomized embeddings is often referred to as
*sketching*.

This section begins with a formal definition of a randomized embedding.
Then we introduce the Gaussian embedding, which is the simplest
construction, and we summarize its analysis. Randomized embeddings have
a wide range of applications in randomized linear algebra. Some
implications of this theory include the Johnson--Lindenstrauss lemma and
a simple construction of a subspace embedding. We also explain why
results for Gaussians transfer to a wider setting. Last, we give a short
description of random partial isometries, a close cousin of Gaussian
embeddings.

## What is a random embedding?

Let $E \subseteq \mathbb{F}^n$ be a set, and let $\varepsilon\in (0, 1)$
be a distortion parameter. We say that a linear map
$\bm{\mathsf{S}} : \mathbb{F}^{n} \to \mathbb{F}^d$ is an ($\ell_2$)
*embedding* of $E$ with distortion $\varepsilon$ when $$\begin{equation}
 \label{eqn:embed-abstract}
(1 - \varepsilon) \, \Vert \bm{\mathsf{x}} \Vert \leq \Vert  \bm{\mathsf{S}} \bm{\mathsf{x}}  \Vert \leq (1 + \varepsilon) \, \Vert  \bm{\mathsf{x}}  \Vert
\quad\text{for all $\bm{\mathsf{x}} \in E$.}
\end{equation}$$ It is sometimes convenient to abbreviate this kind of
two-sided inequality as
$\Vert \bm{\mathsf{S}}\bm{\mathsf{x}} \Vert = (1 \pm \varepsilon) \Vert \bm{\mathsf{x}} \Vert$.

We usually think about the case where $d \ll n$, so the map
$\bm{\mathsf{S}}$ enacts a dimension reduction. In other words,
$\bm{\mathsf{S}}$ transfers data from the high-dimensional space
$\mathbb{F}^n$ to the low-dimensional space $\mathbb{F}^d$. As we will
discuss, the low-dimensional representation of the data can be used to
obtain fast, approximate solutions to computational problems.

The
relation [\[eqn:embed-abstract\]](#eqn:embed-abstract){reference-type="eqref"
reference="eqn:embed-abstract"} expresses the idea that the embedding
$\bm{\mathsf{S}}$ should preserve the geometry of the set $E$.
Unfortunately, we do not always know the set $E$ in advance. Moreover,
we would like the map $\bm{\mathsf{S}}$ to be easy to construct, and it
should be computationally efficient to apply $\bm{\mathsf{S}}$ to the
data. These goals may be in tension.

We can resolve this dilemma by drawing the embedding $\bm{\mathsf{S}}$
from a probability distribution. Many types of probability distributions
serve. In particular, we can use highly structured random matrices
(Section [9](#sec:dimension-reduction){reference-type="ref"
reference="sec:dimension-reduction"}) that are easy to build, to store,
and to apply to vectors.
Section [10](#sec:overdet-ls){reference-type="ref"
reference="sec:overdet-ls"} presents a case study about how random
embeddings can be applied to solve overdetermined least-squares
problems.

## Restricted singular values

Our initial goal is to understand something about the theoretical
behavior of randomized embeddings. To that end, let us introduce
quantities that measure how much an embedding distorts a set. Let
$\bm{\mathsf{S}} \in \mathbb{F}^{d \times n}$ be a linear map, and let
$E \subseteq \mathbb{S}^{n-1}(\mathbb{F})$ be an arbitrary subset of the
unit sphere in $\mathbb{F}^n$. The *minimum* and *maximum restricted
singular value* are, respectively, defined as $$\begin{equation}
 \label{eqn:rsv}
\sigma_{\min}(\bm{\mathsf{S}}; E) := \min_{\bm{\mathsf{x}} \in E}\ \Vert  \bm{\mathsf{S}} \bm{\mathsf{x}}  \Vert
\quad\text{and}\quad
\sigma_{\max}(\bm{\mathsf{S}}; E) := \max_{\bm{\mathsf{x}} \in E}\ \Vert  \bm{\mathsf{S}} \bm{\mathsf{x}}  \Vert.
\end{equation}$$ If $E$ composes the entire unit sphere, then these
quantities coincide with the ordinary minimum and maximum singular value
of $\bm{\mathsf{S}}$. More generally, the restricted singular values
describe how much the linear map $\bm{\mathsf{S}}$ can contract or
expand a point in $E$.

::: remark
**Remark 15** (General sets). *In this treatment, we require $E$ to be a
subset of the unit sphere. Related, but more involved, results hold when
$E$ is a general set. See  and for more results and applications.*
:::

## Gaussian embeddings

Our theoretical treatment of random embeddings focuses on the most
highly structured case. A Gaussian embedding is a random matrix of the
form $$\bm{\mathsf{\Gamma}} \in \mathbb{F}^{d \times n}
\quad\text{with iid entries $(\bm{\mathsf{\Gamma}})_{ij} \sim \textsc{normal}(0, d^{-1})$.}$$
The cost of explicitly storing a Gaussian embedding is $O(dn)$, and the
cost of applying it to a vector is $O(dn)$.

The scaling of the matrix ensures that
$$\operatorname{\mathbb{E}}\Vert  \bm{\mathsf{\Gamma}} \bm{\mathsf{x}}  \Vert^2 = \Vert  \bm{\mathsf{x}}  \Vert^2
\quad\text{for each $\bm{\mathsf{x}} \in \mathbb{F}^n$.}$$ We wish to
understand how large to choose the embedding dimension $d$ so that the
map $\bm{\mathsf{\Gamma}}$ approximately preserves the norms of all
points in a given set $E$. We can do so by obtaining bounds for the
restricted singular values. For a Gaussian embedding, we will see that
$\sigma_{\min}(\bm{\mathsf{\Gamma}}; E)$ and
$\sigma_{\max}(\bm{\mathsf{\Gamma}}; E)$ are controlled by the geometry
of the set $E$.

::: remark
**Remark 16** (Why Gaussians?). *Gaussian embeddings admit a simple and
beautiful analysis. In our computational experience, many other
embeddings exhibit the same (universal) behavior as a Gaussian map. In
spite of that, the rigorous analysis of other types of embeddings tends
to be difficult, even while it yields rather imprecise results. The
confluence of these facts motivates us to argue that the Gaussian
analysis provides enough insight for many practical purposes.*
:::

## The Gaussian width

For the remainder of Section [8](#sec:gauss){reference-type="ref"
reference="sec:gauss"}, we will work in the real field
($\mathbb{F}= \mathbb{R}$). Given a set
$E \subseteq \mathbb{S}^{n-1}(\mathbb{R})$, define the *Gaussian width*
$w(E)$ via
$$w(E) := \operatorname{\mathbb{E}}\sup_{\bm{\mathsf{x}} \in E}\ \langle  \bm{\mathsf{g}} , \,  \bm{\mathsf{x}}  \rangle
\quad\text{where $\bm{\mathsf{g}} \in \mathbb{R}^n$ is standard normal.}$$
The Gaussian width is a measure of the content of the set $E$. It plays
a fundamental role in the performance of randomized embeddings.

Here are some basic properties of the Gaussian width.

- The width is invariant under rotations: $w(\bm{\mathsf{Q}}E) = w(E)$
  for each orthogonal matrix $\bm{\mathsf{Q}}$.

- The width is increasing with respect to set inclusion: $E \subseteq F$
  implies that $w(E) \leq w(F)$.

- The width lies in the range
  $0 \leq w(E) \leq \operatorname{\mathbb{E}}\Vert \bm{\mathsf{g}} \Vert < \sqrt{n}$.

The width can be calculated accurately for many sets of interest. In
particular, if $L$ is an arbitrary $k$-dimensional subspace of
$\mathbb{R}^n$, then $$\begin{equation}
 \label{eqn:width-subspace}
\sqrt{k-1} < w( L \cap \mathbb{S}^{n-1} ) < \sqrt{k}.
\end{equation}$$ Indeed, it is productive to think about the *squared*
width $w^2(E)$ as a measure of the "dimension" of the set $E$.

::: remark
**Remark 17** (Statistical dimension). *The *statistical dimension* is
another measure of content that is closely related to the squared
Gaussian width. The statistical dimension has additional geometric
properties that make it easier to work with in some contexts.
See [@ALMT14:Living-Edge; @MT13:Achievable-Performance; @MT14:Steiner-Formulas]
and [@GNP17:Gaussian-Phase] for more information.*
:::

## Restricted singular values of Gaussian matrices {#sec:rsv-gauss}

In a classic work on Banach space geometry, showed that the Gaussian
width controls both the minimum and maximum restricted singular values
of a subset of the sphere.

::: {#thm:rsv-gauss .theorem}
**Theorem 18** (Restricted singular values: Gaussian matrix). *Fix a
subset $E \subseteq \mathbb{S}^{n-1}(\mathbb{R})$ of the unit sphere.
Draw a Gaussian matrix
$\bm{\mathsf{\Gamma}} \in \mathbb{R}^{d \times n}$ whose entries are iid
$\textsc{normal}(0, d^{-1})$. For all $t > 0$, $$\begin{aligned}
&\mathbb{P}\left\{  \sigma_{\min}(\bm{\mathsf{\Gamma}}; E ) \leq 1 - \frac{w(E) + 1}{\sqrt{d}} - t  \right\}
    &\leq \mathrm{e}^{-dt^2/2}; \\
&\mathbb{P}\left\{  \sigma_{\max}(\bm{\mathsf{\Gamma}}; E) \geq 1 + \frac{w(E)}{\sqrt{d}} + t  \right\}
    &\leq \mathrm{e}^{-dt^2/2}.
\end{aligned}$$*
:::

::: proof
*Proof.* (Sketch) The first inequality is a consequence of Gordon's
minimax theorem and Gaussian concentration. The second inequality is
essentially Chevet's theorem, which follows from Slepian's lemma. See 
for an overview of these ideas. ◻
:::

Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} yields the relations
$$1 - \frac{w(E) + 1}{\sqrt{d}} \lessapprox \sigma_{\min}(\bm{\mathsf{\Gamma}}; E) \leq \sigma_{\max}(\bm{\mathsf{\Gamma}}; E) \lessapprox 1 + \frac{w(E)}{\sqrt{d}}.$$
In other words, the embedding dimension should satisfy
$d > (w(E) + 1)^2$ to ensure that the map $\bm{\mathsf{\Gamma}}$ is
unlikely to annihilate any point in $E$. For this choice of $d$, the
random embedding is unlikely to dilate any point in $E$ by more than a
factor of two.

As a consequence, we have reduced the problem of computing embedding
dimensions for Gaussian maps to the problem of computing Gaussian
widths. In the next two subsections, we work out two important examples.

::: remark
**Remark 19** (Optimality). *The statements in
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} are nearly optimal. One way to see this is to
consider the set $E = \mathbb{S}^{n-1}(\mathbb{R})$, for which the
theorem implies that
$$1 - \sqrt{n/d} \lessapprox \operatorname{\mathbb{E}}\sigma_{\min}(\bm{\mathsf{\Gamma}})
    \leq \operatorname{\mathbb{E}}\sigma_{\max}(\bm{\mathsf{\Gamma}}) \leq 1 + \sqrt{n/d}.$$
The Bai--Yin law [@BS10:Spectral-Analysis Sec. 5.2] confirms that the
first and last inequality are sharp as $n, d \to \infty$ with
$n/d \to \mathrm{const} \in [0,1]$.*

*Moreover, if $E$ is spherically convex (i.e., the intersection of a
convex cone with the unit sphere), then the minimum restricted singular
value satisfies the reverse inequality
$$\mathbb{P}\left\{  \sigma_{\min}(\bm{\mathsf{\Gamma}}; E) \geq 1 - \frac{w(E)}{\sqrt{d}} + t  \right\}
    \leq 2\mathrm{e}^{-dt^2/2}.$$ This result is adapted
from [@TOH14:Gaussian-Minmax].*

*In addition, fifteen years of computational experiments have also shown
that the predictions from
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} are frequently sharp. See  for some examples
and references.*
:::

::: remark
**Remark 20** (History). *The application of Gaussian comparison
theorems in numerical analysis can be traced to work in mathematical
signal processing. used a corollary of Gordon's minimax theorem to study
$\ell_1$ minimization problems. Significant extensions and improvements
of this argument were made by and . seem to have been the first to
recognize that Gordon's minimax theorem can be reversed in the presence
of convexity. A substantial refinement of this observation appeared in .
There is a long series of follow-up works by Babak Hassibi's group that
apply this insight to other problems in signal processing and
communications.*
:::

## Example: Johnson--Lindenstrauss

As a first application of
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"}, let us explain how it implies the classic
dimension reduction result of .

### Overview

Let
$\{ \bm{\mathsf{a}}_1, \dots, \bm{\mathsf{a}}_N \} \subset \mathbb{R}^n$
be a discrete point set. We would like to know when a Gaussian embedding
$\bm{\mathsf{\Gamma}} \in \mathbb{R}^{d \times n}$ approximately
preserves all the pairwise distances between these points:
$$\begin{equation}
 \label{eqn:jl}
1 - \varepsilon\leq \frac{\Vert  \bm{\mathsf{\Gamma}}(\bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j)  \Vert}{ \Vert  \bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j  \Vert }
    \leq 1 + \varepsilon
    \quad\text{for all $i \neq j$.}
\end{equation}$$ The question is how large we must set the embedding
dimension $d$ to achieve distortion $\varepsilon\in (0, 1)$.

### Analysis

We can solve this problem using the machinery described in
Section [8.5](#sec:rsv-gauss){reference-type="ref"
reference="sec:rsv-gauss"}. Consider the set $E$ of normalized chords:
$$E = \left\{ \frac{\bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j}{\Vert \bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j \Vert} : 1 \leq i < j \leq N \right\}.$$
By the definition [\[eqn:rsv\]](#eqn:rsv){reference-type="eqref"
reference="eqn:rsv"} of the restricted singular values,
$$\sigma_{\min}(\bm{\mathsf{\Gamma}}; E) \leq \frac{\Vert  \bm{\mathsf{\Gamma}}(\bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j)  \Vert}{ \Vert  \bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j  \Vert }
    \leq \sigma_{\max}(\bm{\mathsf{\Gamma}}; E).$$ Therefore, we can
invoke Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} to determine how the embedding dimension
controls the distortion.

Let us summarize the argument. First, observe that the Gaussian width of
the set $E$ satisfies
$$w(E) = \operatorname{\mathbb{E}}\max_{\bm{\mathsf{x}} \in E}\ \langle  \bm{\mathsf{g}} , \,  \bm{\mathsf{x}}  \rangle
    \leq \sqrt{2 \log \# E}
    < 2 \sqrt{\log (N/2)}.$$ As a consequence, $$\begin{aligned}
&\mathbb{P}\left\{  \sigma_{\min}(\bm{\mathsf{\Gamma}}; E) \leq 1 - (1 + 2 \sqrt{\log(N/2)})/\sqrt{d} - t  \right\} &\leq \mathrm{e}^{-d t^2 /2}; \\
&\mathbb{P}\left\{  \sigma_{\max}(\bm{\mathsf{\Gamma}}; E) \geq 1 + 2 \sqrt{\log(N/2)}/\sqrt{d} + t  \right\} &\leq \mathrm{e}^{-d t^2 /2}.
\end{aligned}$$ To achieve distortion $\varepsilon$ with high
probability, it is sufficient to choose
$$d \geq 8 \varepsilon^{-2} \log N.$$ In other words, the embedding
dimension only needs to be *logarithmic* in the cardinality $N$ of the
point set. With some additional calculation, we can also extract precise
failure probabilities from this analysis.

### Discussion

Let us close this example with a few comments. In spite of its
prominence, the Johnson--Lindenstrauss embedding lemma is somewhat
impractical. Indeed, since the embedding dimension $d$ is proportional
to $\varepsilon^{-2}$, it is a challenge to achieve small distortions.
Even if we consider the setting where $\varepsilon\approx 1$, the
uniform bound [\[eqn:jl\]](#eqn:jl){reference-type="eqref"
reference="eqn:jl"} may require the embedding dimension to be
prohibitively large.

As a step toward more applicable results, note that the bound on the
*minimum* restricted singular value is more crucial than the bound on
the maximum restricted singular value, because the former ensures that
no two points coalesce after the random embedding. Similarly, it is
often more valuable to preserve the distances between nearby points than
between far-flung points. This observation is the starting point for the
theory of locality sensitive hashing [@GIM99:Similarity-Search].

### History

were concerned with a problem in Banach space geometry, namely the
prospect of extending a Lipschitz function from a finite metric space
into a Hilbert space. The famous lemma from their paper took on a life
of its own when  used it to design efficient approximation algorithms
for some graph problems. used random embeddings to develop new
algorithms for the approximate nearest neighbor problem.
[@AMS99:Space-Complexity; @AGMS02:Tracking-Join] introduced the term
*sketching*, and they showed how to use sketches to track streaming
data. Soon afterwards, and proposed using random embeddings and matrix
sampling for low-rank matrix approximation, bringing these ideas into
the realm of computational linear algebra.

## Example: Subspace embedding {#sec:subspace-embedding}

Next, we consider a question at the heart of randomized linear algebra.
Can we embed an unknown subspace into a lower-dimensional space?

### Overview

Suppose that $L$ is a $k$-dimensional subspace in $\mathbb{R}^n$. We say
that a dimension reduction map $\bm{\mathsf{S}}$ is a *subspace
embedding* for $L$ with distortion $\varepsilon\in (0, 1)$ if
$$\begin{equation}
 \label{eqn:subspace-embedding}
 (1 - \varepsilon) \, \Vert  \bm{\mathsf{x}}  \Vert
    \leq \Vert  \bm{\mathsf{S}}\bm{\mathsf{x}}  \Vert
    \leq (1 + \varepsilon) \, \Vert  \bm{\mathsf{x}}  \Vert
    \quad\text{for every $\bm{\mathsf{x}} \in L$.}
\end{equation}$$ We say that $\bm{\mathsf{S}}$ is *oblivious* if it can
be constructed without knowledge of the subspace $L$, except for its
dimension.

Two questions arise. First, what types of dimension reduction maps yield
(oblivious) subspace embeddings? Second, how large must we choose the
embedding dimension to achieve this outcome?

### Analysis

Gaussian dimension reduction maps yield very good oblivious subspace
embeddings. Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} easily furnishes the justification. Consider
the unit sphere in the subspace:
$E = L \cap \mathbb{S}^{n-1}(\mathbb{R})$. Then construct the Gaussian
dimension reduction map
$\bm{\mathsf{\Gamma}} \in \mathbb{R}^{d \times n}$. In view
of [\[eqn:width-subspace\]](#eqn:width-subspace){reference-type="eqref"
reference="eqn:width-subspace"}, we have $$\begin{aligned}
&\mathbb{P}\left\{  \sigma_{\min}( \bm{\mathsf{\Gamma}}; E ) \leq 1 - (1+\sqrt{k})/\sqrt{d} - t  \right\} &\leq \mathrm{e}^{-d t^2/2}; \\
&\mathbb{P}\left\{  \sigma_{\max}( \bm{\mathsf{\Gamma}}; E ) \geq 1 + \sqrt{k}/\sqrt{d} + t  \right\} &\leq \mathrm{e}^{-d t^2/2}. \\
\end{aligned}$$ As a specific example, we can set the embedding
dimension $d = 2k$ to ensure that
$\Vert \bm{\mathsf{\Gamma}}\bm{\mathsf{x}} \Vert = (1 \pm 0.8) \Vert \bm{\mathsf{x}} \Vert$
simultaneously for all points $\bm{\mathsf{x}} \in L$, except with
probability $\mathrm{e}^{-\mathrm{c} k}$. In some applications of
subspace embeddings, we can even choose the dimension as small as
$d = k + 5$ or $d = k + 10$.

Many theoretical papers on randomized NLA use subspace embeddings as a
primitive for designing algorithms for other linear algebra problems.
For example, Section [10](#sec:overdet-ls){reference-type="ref"
reference="sec:overdet-ls"} describes several ways to use subspace
embeddings to solve overdetermined least-squares problems.

### History

Subspace embeddings were explicitly introduced by ; see also . As work
on randomized NLA accelerated, researchers became interested in more
structured types of subspace embeddings; an early reference is .
Section [9](#sec:dimension-reduction){reference-type="ref"
reference="sec:dimension-reduction"} covers these extensions. See  for a
theoretical perspective on randomized NLA where subspace embeddings take
pride of place.

## Universality of the minimum restricted singular value

We have seen how to apply Gaussian dimension reduction for embedding
discrete point sets and for embedding subspaces.
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} contains precise theoretical results on the
behavior of Gaussian maps in terms of the Gaussian width. To what extent
can we transfer this analysis to other types of random embeddings?

The following theorem [@OT18:Universality-Laws Thm. 9.1] shows that the
bound on the *minimum* restricted singular value in
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} is universal for a large class of random
embeddings. In particular, this class includes sparse random matrices,
whose nonzero entries compose a vanishing proportion of the total.

::: {#thm:universality .theorem}
**Theorem 21** (Universality). *Fix a set
$E \subseteq \mathbb{S}^{n-1}$. Let
$\bm{\mathsf{S}} \in \mathbb{R}^{d \times n}$ be a random matrix whose
entries are independent random variables that satisfy $$%
\operatorname{\mathbb{E}}[ (\bm{\mathsf{S}})_{ij} ] = 0, \quad
\operatorname{\mathbb{E}}[ (\bm{\mathsf{S}})_{ij}^2 ] = d^{-1}, \quad
\operatorname{\mathbb{E}}[ (\bm{\mathsf{S}})_{ij}^5 ] \leq R. %$$ When
$d \leq n$, with high probability, $$\sigma_{\min}(\bm{\mathsf{S}}; E)
    \geq 1 - \frac{w(E)}{\sqrt{d}} - o(\sqrt{n/d}).$$ The constant in
$o(\sqrt{n/d})$ depends only on $R$. A matching lower bound for
$\sigma_{\min}(\bm{\mathsf{S}}; E)$ holds when $E$ is spherically
convex.*
:::

In other words, if $E$ is a moderately large set, the distribution of
the entries of the random map $\bm{\mathsf{S}}$ does not have an impact
on the embedding dimension $d$ sufficient to ensure no point in $E$ is
annihilated.

Theorem [21](#thm:universality){reference-type="ref"
reference="thm:universality"} is confirmed by extensive numerical
experiments [@OT18:Universality-Laws], which demonstrate that dimension
reduction maps with independent, standardized entries have identical
performance for a wide range of examples.

It is perhaps surprising that the bound on the *maximum* restricted
singular value from Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} is not universal. For some sets $E$, the
quantity $\sigma_{\max}(\bm{\mathsf{S}}; E)$ depends heavily on the
distribution of the entries of $\bm{\mathsf{S}}$.

::: remark
**Remark 22** (Universality for least-squares). *give some asymptotic
universality results for random embeddings in the context of
least-squares problems.*
:::

## Random partial isometries

Last, we consider a variant of Gaussian embedding that is more suitable
when the embedding dimension $d$ is close to the ambient dimension $n$.
In this section, we allow the field $\mathbb{F}$ to be real or complex.

First, suppose that $d \leq n$, and let
$\bm{\mathsf{\Gamma}} \in \mathbb{F}^{d \times n}$ be a Gaussian
embedding. Almost surely, the co-range of $\bm{\mathsf{\Gamma}}$ is a
uniformly random $d$-dimensional subspace of $\mathbb{F}^n$. Construct
an embedding $\bm{\mathsf{S}} \in \mathbb{F}^{d \times n}$ with
orthonormal rows that span the co-range of $\bm{\mathsf{\Gamma}}$, for
example by QR factorization.

Similarly, we can consider a Gaussian embedding
$\bm{\mathsf{\Gamma}} \in \mathbb{F}^{d \times n}$ with $d \geq n$. In
this case, the range of $\bm{\mathsf{\Gamma}}$ is almost surely a
uniformly random $n$-dimensional subspace of $\mathbb{F}^d$. Construct
an embedding $\bm{\mathsf{S}} \in \mathbb{F}^{d \times n}$ with
orthonormal columns that span the range of $\bm{\mathsf{\Gamma}}$, for
example by QR factorization.

In each case, we call $\bm{\mathsf{S}}$ a *random partial isometry*. The
cost of storing a random partial isometry is $O(dn)$, and the cost of
applying it to a vector is $O(dn)$. (We should warn the punctilious
reader that QR factorization of $\bm{\mathsf{\Gamma}}$ may not produce a
matrix $\bm{\mathsf{S}}$ that is Haar-distributed on the Stiefel
manifold. To achieve this guarantee, use the algorithms from .)

When $d \approx n$, random partial isometries are better embeddings than
Gaussian maps (because the nonzero singular values of a partial isometry
are all equal). When $d$ and $n$ are significantly different, the two
models are quite similar to each other.

have established some theoretical results on the embedding behavior of
real partial isometries $(\mathbb{F}= \mathbb{R})$. Unfortunately, the
situation is more complicated than in
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"}. More relations between Gaussian matrices and
partial isometries follow from the Marcus--Pisier comparison
theorem [@MP81:Random-Fourier]; see also .

# Structured random embeddings {#sec:dimension-reduction}

Gaussian embeddings and random partial isometries work extremely well.
But they are not suitable for all practical applications because they
are expensive to construct, to store, and to apply to vectors. Instead,
we may prefer to implement more structured embedding matrices that
alleviate these burdens.

This section summarizes a number of constructions that have appeared in
the literature, with a focus on methods that have been useful in
applications. Except as noted, these approaches have the same practical
performance as either a Gaussian embedding or a random partial isometry.

Although many of these approaches are supported by theoretical analysis,
the results are far less precise than for Gaussian embeddings. As such,
we will not give detailed mathematical statements about structured
embeddings. See Section [9.7](#sec:how-in-theory){reference-type="ref"
reference="sec:how-in-theory"} for a short discussion about how to
manage the lack of theoretical guarantees.

## General techniques

Many types of structured random embeddings operate on the same
principle, articulated in [@AC09:Fast-Johnson-Lindenstrauss]. When we
apply the random embedding to a fixed vector, it should homogenize
("mix") the coordinates so that each one carries about the same amount
of energy. Then the embedding can sample coordinates at random to
extract a lower-dimensional vector whose norm is proportional to the
norm of the original vector and has low variance. Random embeddings
differ in how they perform the initial mixing step. Regardless of how it
is done, mixing is very important for obtaining embeddings that work
well in practice.

With this intuition at hand, let us introduce a number of pre- and
post-processing transforms that help us design effective random
embeddings. These approaches are used in many of the constructions
below.

We say that a random variable is a *random sign* if it is
$\textsc{uniform}\{ z \in \mathbb{F}: \vert z \vert = 1 \}$. A random
matrix $\bm{\mathsf{E}} \in \mathbb{F}^{n \times n}$ is called a *random
sign flip* if it is diagonal, and the diagonal entries are iid random
sign variables.

A *random permutation* $\bm{\mathsf{\Pi}} \in \mathbb{F}^{n \times n}$
is a matrix drawn uniformly at random from the set of permutation
matrices. That is, each row and column of $\bm{\mathsf{\Pi}}$ has a
single nonzero entry, which equals one, and all such matrices are
equally likely.

For $d \leq n$, a random matrix
$\bm{\mathsf{R}} \in \mathbb{F}^{d \times n}$ is called a *random
restriction* if it selects $d$ uniformly random entries from its input.
With an abuse of terminology, we extend the definition to the case
$d \geq n$ by making $\bm{\mathsf{R}} \in \mathbb{F}^{d \times n}$ the
matrix that embeds its input into the first $n$ coordinates of the
output. That is, $(\bm{\mathsf{R}})_{ij} = 1$ when $i = j$ and zero
otherwise.

Random sign flips and permutations are useful for preconditioning the
input to a random embedding. Random restrictions are useful for reducing
the dimension of a vector that has already been homogenized.

## Sparse sign matrices {#sec:sparse-map}

Among the earliest proposals for non-Gaussian embedding is to use a
sparse random matrix whose entries are random signs.

Here is an effective construction of a sparse sign matrix
$\bm{\mathsf{S}} \in \mathbb{F}^{d \times n}$. Fix a sparsity parameter
$\zeta$ in the range $2 \leq \zeta \leq d$. The random embedding takes
the form
$$\bm{\mathsf{S}} = \sqrt{\frac{n}{\zeta}} \begin{bmatrix} \bm{\mathsf{s}}_1 & \dots & \bm{\mathsf{s}}_n \end{bmatrix} \in \mathbb{F}^{d \times n}.$$
The columns $\bm{\mathsf{s}}_i \in \mathbb{F}^d$ are iid random vectors.
To construct each column, we draw $\zeta$ iid random signs, and we
situate them in $\zeta$ uniformly random coordinates. recommend choosing
$\zeta = \min\{d, 8\}$ in practice.

We can store a sparse embedding using about $O(\zeta n \log d)$ numbers.
We can apply it to a vector in $\mathbb{F}^n$ with $O(\zeta n)$
arithmetic operations. The main disadvantage is that we must use sparse
data structures and arithmetic to achieve these benefits. Sparse sign
matrices have similar performance to Gaussian embeddings.

has shown that a sparse sign matrix serves as an oblivious subspace
embedding with constant distortion for an arbitrary $k$-dimensional
subspace of $\mathbb{R}^n$ when the embedding dimension
$d = O(k \log k)$ and the per-column sparsity $\zeta = O(\log k)$. It is
conjectured that improvements are still possible.

::: remark
**Remark 23** (History). *Sparse random embeddings emerged from the work
of and . For randomized linear algebra applications, sparse embeddings
were promoted
in [@CW13:Low-Rank-Approximation; @MM13:Low-Distortion-Subspace; @NN13:OSNAP-Faster]
and [@Ura13:Fast-Randomized]. Analyses of the embedding behavior of a
sparse map appear in [@BDN15:Toward-Unified]
and [@Coh16:Nearly-Tight-Oblivious].*
:::

## Subsampled trigonometric transforms {#sec:srtt}

Another type of structured randomized embeddings is designed to mimic
the performance of a random partial isometry. One important class of
examples consists of the subsampled randomized trigonometric transforms
(SRTTs).

To construct a random embedding
$\bm{\mathsf{S}} \in \mathbb{F}^{d \times n}$ with $d \leq n$, we select
a unitary trigonometric transform
$\bm{\mathsf{F}} \in \mathbb{F}^{n \times n}$. Then we form
$$\bm{\mathsf{S}} = \sqrt{\frac{n}{d}} \bm{\mathsf{RFE\Pi}},$$ where
$\bm{\mathsf{R}} \in \mathbb{F}^{d\times n}$ is a random restriction,
$\bm{\mathsf{E}} \in \mathbb{F}^{n \times n}$ is a random sign flip, and
$\bm{\mathsf{\Pi}} \in \mathbb{F}^{n \times n}$ is a random permutation.
Note that $\bm{\mathsf{S}}$ is a partial isometry.

The trigonometric transform $\bm{\mathsf{F}}$ can be any one of the
usual suspects. In the complex case ($\mathbb{F}= \mathbb{C}$), we often
use a discrete Fourier transform. In the real case
($\mathbb{F}= \mathbb{R}$), common choices are the discrete cosine
transform (DCT2) or the discrete Hartley transform (DHT). When $n$ is a
power of two, we can consider a Walsh--Hadamard transform (WHT). The
paper [@2010_avron_BLENDENPIK] reports that the DHT is the best option
in the real case.

The cost of storing a SRTT is $O(n \log n)$, and it can be applied to a
vector in $O(n \log d)$ operations using a fast subsampled trigonometric
transform algorithm. The main disadvantage is that it requires a good
implementation of the fast transform.

has shown that an SRTT serves as an oblivious subspace embedding with
constant distortion for an arbitrary $k$-dimensional subspace of
$\mathbb{F}^n$ provided that $d = O(k \log k)$. This paper focuses on
the Walsh--Hadamard transform, but the analysis extends to other SRTTs.
In practice, it often suffices to choose $d = O(k)$, but no rigorous
justification is available.

::: remark
**Remark 24** (Rerandomization). *It is also common to repeat the
randomization and trigonometric transformations:
$$\bm{\mathsf{S}} = \sqrt{\frac{n}{d}} \bm{\mathsf{R}}\bm{\mathsf{F}}\bm{\mathsf{E}}'\bm{\mathsf{\Pi}}'\bm{\mathsf{FE\Pi}},$$
with an independent sign flip $\bm{\mathsf{E}}'$ and an independent
permutation $\bm{\mathsf{\Pi}}'$. This enhancement can make the
embedding more robust, although it is not always necessary.*
:::

::: remark
**Remark 25** (History). *proposed the use of randomized trigonometric
transforms to precondition linear systems. The idea of applying an SRTT
for dimension reduction appears in [@2006_ailon_chazelle_FJLT] and
[@AC09:Fast-Johnson-Lindenstrauss]. develop algorithms for low-rank
matrix approximation based on SRTTs. Embedding properties of an SRTT for
general sets follow from [@RV08:Sparse-Reconstruction] and
[@KW11:New-Improved]; see or .*
:::

## Tensor random projections {#sec:trp}

Next, we describe a class of random embeddings that are useful for very
large linear algebra and multilinear algebra problems. This approach
invokes tensor products to form a random embedding for a
high-dimensional space from a family of random embeddings for
lower-dimensional spaces.

Let $\bm{\mathsf{S}}_1 \in \mathbb{F}^{d \times m_1}$ and
$\bm{\mathsf{S}}_2 \in \mathbb{F}^{d \times m_2}$ be statistically
independent random embeddings. We define the *tensor random embedding*
$$\bm{\mathsf{S}} := \bm{\mathsf{S}}_1 \odot \bm{\mathsf{S}}_2 \in \mathbb{F}^{d \times n}
\quad\text{where}\quad
n = m_1 m_2$$ to be the Khatri--Rao product of $\bm{\mathsf{S}}_1$ and
$\bm{\mathsf{S}}_2$. That is, the $i$th row of $\bm{\mathsf{S}}$ is
$$(\bm{\mathsf{S}})_{i:} = \begin{bmatrix} (\bm{\mathsf{S}}_1)_{i1} (\bm{\mathsf{S}}_2)_{i:} & \dots & (\bm{\mathsf{S}}_1)_{im} (\bm{\mathsf{S}}_2)_{i:} \end{bmatrix}
    \quad\text{for $i=1,\dots,d$.}$$ Under moderate assumptions on the
component embeddings $\bm{\mathsf{S}}_1$ and $\bm{\mathsf{S}}_2$, the
tensor random embedding $\bm{\mathsf{S}}$ preserves the squared
Euclidean norm of an arbitrary vector in $\mathbb{F}^{n}$.

A natural extension of this idea is to draw many component embeddings
$\bm{\mathsf{S}}_i \in \mathbb{F}^{d \times m_i}$ for $i = 1, \dots, k$
and to form the tensor random embedding
$$\bm{\mathsf{S}} := \bm{\mathsf{S}}_1 \odot \bm{\mathsf{S}}_2 \odot \dots \odot \bm{\mathsf{S}}_k \in \mathbb{F}^{d \times n}
\quad\text{where}\quad
n = \prod\nolimits_{i=1}^k m_i.$$ This embedding also inherits nice
properties from its components.

The striking thing about this construction is that the tensor product
embedding operates on a *much* larger space than the component
embeddings. The storage cost for the component embeddings is
$O(d (\sum_{i=1}^k m_i))$, or less. We can apply the tensor random
embedding to a vector directly with $O(dn)$ arithmetic. We can
accelerate the process by using component embeddings that have fast
transforms, and we can obtain improvements for vectors that have a
compatible tensor product structure. Some theoretical analysis is
available, but results are not yet complete.

::: remark
**Remark 26** (History). *Tensor random embeddings were introduced by 
on differential privacy. They were first analyzed by . The
paper [@SGTU18:Tensor-Random] proposed the application of tensor random
embeddings for randomized linear algebra; some extensions appear
in [@JKW19:Faster-Johnson-Lindenstrauss] and
[@MB19:Guarantees-Kronecker]. See [@BV19:Polynomial-Threshold] and
[@Ver19:Concentration-Inequalities] for related theoretical results.*
:::

## Other types of structured random embeddings

We have described the random embeddings that have received the most
attention in the NLA literature. Yet there are other types of random
embeddings that may be useful in special circumstances. Some examples
include random filters
[@TWDBB06:Random-Filters; @KW11:New-Improved; @RRT12:Restricted-Isometries; @MRW18:Improved-Bounds],
the Kac random
walk [@Kac56:Foundations-Kinetic; @Ros94:Random-Rotations; @Oli09:Convergence-Equilibrium; @PS17:Kacs-Walk],
and sequences of random
reflections [@Slo83:Encrypting-Random; @Por96:Cutoff-Phenomenon].
discusses a number of other instances.

## Coordinate sampling {#sec:coord-embed}

So far, we have discussed random embeddings that mix up the coordinates
of a vector. It is sometimes possible to construct embeddings just by
sampling coordinates at random. Coordinate sampling can be appealing in
specialized situations (e.g., kernel computations), where we only have
access to individual entries of the data. On the other hand, this
approach requires strong assumptions, and it is far less reliable than
random embeddings that mix coordinates. In this section, we summarize
the basic facts about subspace embedding via random coordinate sampling.

A note on terminology: we will use the term *coordinate sampling* to
distinguish these maps from random embeddings that mix coordinates.

### Coherence and leverage {#sec:coherence_leverage}

Let $L \subset \mathbb{F}^{n}$ be a $k$-dimensional subspace. The
*coherence* $\mu(L)$ of the subspace with respect to the standard basis
is
$$\mu(L) := n \cdot \max_{i = 1, \dots, n} \Vert  \bm{\mathsf{P}}_L \bm{\mathsf{\delta}}_i  \Vert^2,
%$$ where $\bm{\mathsf{P}}_L \in \mathbb{H}_n$ is the orthogonal
projector onto $L$ and $\bm{\mathsf{\delta}}_i$ is the $i$th standard
basis vector. The coherence $\mu(L)$ lies in the range $[k, n]$. The
behavior of coordinate sampling methods degrades as the coherence
increases.

Next, define the (subspace) *leverage score* distribution with respect
to the standard coordinate basis:
$$p_i = \frac{1}{k} \Vert  \bm{\mathsf{P}}_L \bm{\mathsf{\delta}}_i  \Vert^2
\quad\text{for $i = 1, \dots, n$.}$$ It is straightforward to verify
that $(p_1, \dots, p_n)$ is a probability distribution. In most
applications, it is expensive to compute or estimate subspace leverage
scores because we typically do not have a basis for the subspace $L$ at
hand.

### Uniform sampling {#uniform-sampling}

We can construct an embedding
$\bm{\mathsf{S}} \in \mathbb{F}^{d \times n}$ by sampling each output
coordinate uniformly at random. That is, the rows of $\bm{\mathsf{S}}$
are iid, and each row takes values $\bm{\mathsf{\delta}}_i / \sqrt{d}$,
each with probability $1/n$. (We can also sample coordinates uniformly
*without* replacement; this approach performs slightly better but
requires more work to analyze.)

The embedding dimension $d$ must be chosen to ensure that
$$\begin{equation}
 \label{eqn:sample-embed}
\Vert \bm{\mathsf{S}} \bm{\mathsf{x}}  \Vert^2 = (1 \pm \varepsilon) \Vert  \bm{\mathsf{x}}  \Vert^2
\quad\text{for all $\bm{\mathsf{x}} \in L$.}
\end{equation}$$ To achieve this goal, it suffices that
$$d \geq 2 \varepsilon^{-2} \mu(L) \log(2k).$$ In other words, the
embedding dimension is proportional to the coherence of the subspace, up
to a logarithmic factor. We expect uniform sampling to work well
precisely when the coherence is small ($\mu(L) \approx k$).

To prove this result, let $\bm{\mathsf{U}} \in \mathbb{F}^{n \times k}$
be an orthonormal basis for the subspace $L$. We can approximate the
product $\bm{\mathsf{I}}_k = \bm{\mathsf{U}}^*\bm{\mathsf{U}}$ by
sampling columns of $\bm{\mathsf{U}}$ uniformly at random. The analysis
in Section [7.3.2](#sec:mtx-mult-unif){reference-type="ref"
reference="sec:mtx-mult-unif"} furnishes the conclusion.

### Leverage score sampling

When the coherence is large, it seems more natural to sample with
respect to the leverage score distribution $(p_1, \dots, p_n)$ described
above. That is, the embedding
$\bm{\mathsf{S}} \in \mathbb{F}^{d \times n}$ has iid rows, and each row
takes value $\bm{\mathsf{\delta}}_i / \sqrt{d}$ with probability $p_i$.

To achieve the embedding
guarantee [\[eqn:sample-embed\]](#eqn:sample-embed){reference-type="eqref"
reference="eqn:sample-embed"} with this sampling distribution, we should
choose the embedding dimension $$d \geq 2 \varepsilon^{-2} k \log(2k).$$
In other words, it suffices that the embedding dimension is proportional
to the dimension $k$ of the subspace, up to the logarithmic factor. This
result follows from the analysis of matrix multiplication by importance
sampling (Section [7.3.3](#sec:mtx-mult-import){reference-type="ref"
reference="sec:mtx-mult-import"}).

### Discussion

Uniform sampling leads to an oblivious subspace embedding, although it
is not obvious how to select the embedding dimension in advance because
the coherence is usually not available. Leverage score sampling is
definitely not oblivious, because we need to compute the sampling
probabilities (potentially at great cost).

In practice, uniform sampling works better than one might anticipate,
and it has an appealing computational profile. As a consequence, it has
become a workhorse for large-scale kernel computation;
see [@KMT12:Sampling-Methods; @2013_bach_sharp] and
[@RCR17:FALKON-Optimal].

Our experience suggests that leverage score sampling is rarely a
competitive method for constructing subspace embeddings, especially once
we take account of the effort required to compute the sampling
probabilities. We recommend using other types of random embeddings
(Gaussians, sparse maps, SRTTs) in lieu of coordinate sampling whenever
possible.

Although coordinate sampling may seem like a natural approach to
construct matrix approximations involving rows or columns, we can obtain
better algorithms for this problem by using (mixing) random embeddings.
See Section [13](#sec:natural){reference-type="ref"
reference="sec:natural"} for details.

Coordinate sampling can be a compelling choice in situations where other
types of random embeddings are simply unaffordable. For instance, a
variant of leverage score sampling leads to effective algorithms for
kernel ridge regression; see  for evidence. Indeed, in the context of
kernel computations, coordinate sampling and random features may be the
only tractable methods for extracting information from the kernel
matrix. We discuss these ideas in
Section [19](#sec:kernel){reference-type="ref" reference="sec:kernel"}.

An emerging research direction uses coordinate sampling for solving
certain kinds of continuous problems, such as function interpolation. In
this setting, sampling corresponds to function evaluation, while mixing
embeddings may lead to operations that are impossible to implement in
the continuous space. For some examples, see
[@RW12:Sparse-Legendre; @CDL13:Stability-Accuracy; @HD15:Coherence-Motivated; @RW16:Interpolation-Weighted; @CM17:Optimal-Weighted; @ABC19:Sequential-Sampling; @AKM+19:Universal-Sampling]
and [@CP19:Active-Regression].

### History

Most of the early theoretical computer science papers on randomized NLA
rely on coordinate sampling methods. These approaches typically
construct an importance sampling distribution using the norms of the
rows or columns of a matrix. For examples,
see [@FKV04:Fast-Monte-Carlo; @2005_drineas_nystrom]
and [@DKM06:Fast-Monte-Carlo-I; @DKM06:Fast-Monte-Carlo-II; @DKM06:Fast-Monte-Carlo-III].
These papers measure errors in the Frobenius norm. The first spectral
norm analysis of coordinate sampling appears in .

Leverage scores are a classical tool in statistical regression, used to
identify influential data points. proposed the subspace leverage scores
as a sampling distribution for constructing low-rank matrix
approximations. identified the connection with regression. made a
theoretical case for using leverage scores as the basis for randomized
NLA algorithms. introduced an alternative definition of leverage scores
for kernel ridge regression.

It is somewhat harder to trace the application of uniform sampling in
randomized NLA. Several authors have studied the behavior of uniform
sampling in the context of Nyström approximation; see
[@WS01:Using-Nystrom; @KMT12:Sampling-Methods] and
[@Git13:Topics-Randomized]. An analysis of uniform coordinate sampling
is implicit in the theory on SRTTs; see . See for more discussion about
sampling methods in NLA.

## But how does it work in theory? {#sec:how-in-theory}

Structured random embeddings and random coordinate sampling lack the
precise guarantees that we can attribute to Gaussian embedding matrices.
So how can we apply them with confidence?

First, we advocate using *a posteriori* error estimators to assess the
quality of the output of a randomized linear algebra computation. These
error estimators are often quite cheap, yet they can give (statistical)
evidence that the computation was performed correctly. We also recommend
adaptive algorithms that can detect when the accuracy is insufficient
and make refinements. With this approach, it is not pressing to produce
theory that justifies all of the internal choices (e.g., the specific
type of random embedding) in the NLA algorithm. See
Section [12](#sec:error-est){reference-type="ref"
reference="sec:error-est"} for further discussion.

Even so, we would like to have *a priori* predictions about how our
algorithms will behave. Beyond that, we need reliable methods for
selecting algorithm parameters, especially in the streaming setting
where we cannot review the data and repeat the computation.

Here is one answer to these concerns. As a practical matter, we can
simply invoke the lessons from the Gaussian theory, even when we are
using a different type of random embedding. The universality result,
Theorem [21](#thm:universality){reference-type="ref"
reference="thm:universality"}, gives a rationale for this approach in
one special case. We also recommend undertaking computational
experiments to verify that the Gaussian theory gives an adequate
description of the observed behavior of an algorithm.

::: warning
**Warning 27** (Coordinate sampling). *Mixing random embeddings perform
similarly to Gaussian embeddings, but coordinate sampling methods
typically exhibit behavior that is markedly worse.*
:::

# How to use random embeddings {#sec:overdet-ls}

Algorithm designers have employed random embeddings for many tasks in
linear algebra, optimization, and related areas. Methods based on random
embedding fall into three rough categories: (1) sketch and solve, (2)
iterative sketching, and (3) sketch and precondition. To draw
distinctions among these paradigms, we use each one to derive an
algorithm for solving an overdetermined least-squares problem.

## Overdetermined least-squares

Overdetermined least-squares problems sometimes arise in statistics and
data-analysis applications. We may imagine that some of the data in
these problems is redundant. As such, it seems plausible that we could
reduce the size of the problem to accelerate computation without too
much loss in accuracy.

Consider a matrix $\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$ with
$m \gg n$ and a vector $\bm{\mathsf{b}} \in \mathbb{F}^m$. An
overdetermined least-squares problem has the form $$\begin{equation}
 \label{eqn:overdet-ls}
\underset{\bm{\mathsf{x}} \in \mathbb{F}^n}{\text{minimize}} \
    \frac{1}{2} \Vert  \bm{\mathsf{A}} \bm{\mathsf{x}} - \bm{\mathsf{b}}  \Vert^2.
\end{equation}$$ Following , we can also rewrite the least-squares
problem to emphasize the role of the matrix: $$\begin{equation}
 \label{eqn:overdet-ls2}
\underset{\bm{\mathsf{x}} \in \mathbb{F}^n}{\text{minimize}} \
    \frac{1}{2} \Vert  \bm{\mathsf{A}} \bm{\mathsf{x}}  \Vert^2 - \langle  \bm{\mathsf{x}} , \,  \bm{\mathsf{A}}^* \bm{\mathsf{b}}  \rangle.
\end{equation}$$ Write $\bm{\mathsf{x}}_{\star}$ for an arbitrary
solution to the
problem [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"}.

To make a clear comparison among algorithm design templates, we will
assume that $\bm{\mathsf{A}}$ is dense and unstructured. In this case,
the classical approach to
solving [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"} is based on factorization of the coefficient
matrix (such as QR or SVD) at a cost of $O(mn^2)$ arithmetic operations.

When $\bm{\mathsf{A}}$ is sparse, we would typically use iterative
methods (such as CG), which have a different computational profile. For
sparse matrices, we would also make different design choices in a
sketch-based algorithm. Nevertheless, for simplicity, we will not
discuss the sparse case.

## Subspace embeddings for least-squares

To design a sketching algorithm for the overdetermined least-squares
problem [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"}, we need to construct a subspace embedding
$\bm{\mathsf{S}} \in \mathbb{F}^{d \times m}$ that preserves the
geometry of the range of the matrix
$\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$. In some cases, we may
also need the embedding to preserve the range of the bordered matrix
$\begin{bmatrix} \bm{\mathsf{A}} & \bm{\mathsf{b}} \end{bmatrix} \in \mathbb{F}^{m \times (n+1)}$.

Since the matrix $\bm{\mathsf{A}}$ is dense and unstructured, we will
work with a structured subspace embedding, such as an SRTT
(Section [9.3](#sec:srtt){reference-type="ref" reference="sec:srtt"}).
More precisely, we will assume that evaluating the product
$\bm{\mathsf{SA}}$ costs only $O(mn \log d)$ arithmetic operations. The
best theoretical results for these structured sketches require that the
embedding dimension $d \sim n \log(n) / \varepsilon^2$ to achieve
distortion $\varepsilon$, although the logarithmic factor seems to be
unnecessary in practice.

Throughout this section, we use the heuristic notation $\sim$ to
indicate quantities that are proportional. We also write $\ll$ to mean
"much smaller than."

## Sketch and solve {#sec:sketchandsolve}

The sketch-and-solve paradigm maps the overdetermined least-squares
problem [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"} into a smaller space. Then it uses the
solution to the reduced problem as a proxy for the solution to the
original problem. This approach can be very fast, and we only need one
view of the matrix $\bm{\mathsf{A}}$. On the other hand, the results
tend to be very inaccurate.

Let $\bm{\mathsf{S}} \in \mathbb{F}^{d \times m}$ be a subspace
embedding for the range of
$\begin{bmatrix} \bm{\mathsf{A}} & \bm{\mathsf{b}} \end{bmatrix}$ with
distortion $\varepsilon$. Consider the compressed least-squares problem
$$\begin{equation}
 \label{eqn:sketch-solve-ls}
\underset{\bm{\mathsf{x}} \in \mathbb{R}^n}{\text{minimize}} \
    \frac{1}{2} \Vert  \bm{\mathsf{S}} (\bm{\mathsf{A}} \bm{\mathsf{x}} - \bm{\mathsf{b}})  \Vert^2.
\end{equation}$$ Since $\bm{\mathsf{S}}$ preserves geometry, we may hope
that the solution $\widehat{\bm{\mathsf{x}}}$ to the sketched
problem [\[eqn:sketch-solve-ls\]](#eqn:sketch-solve-ls){reference-type="eqref"
reference="eqn:sketch-solve-ls"} can replace the solution
$\bm{\mathsf{x}}_{\star}$ to the original
problem [\[eqn:overdet-ls2\]](#eqn:overdet-ls2){reference-type="eqref"
reference="eqn:overdet-ls2"}. A typical theoretical bound is
$$\Vert  \bm{\mathsf{A}} \widehat{\bm{\mathsf{x}}} - \bm{\mathsf{b}}  \Vert \leq (1 + \varepsilon) \Vert  \bm{\mathsf{A}} \bm{\mathsf{x}}_{\star} - \bm{\mathsf{b}}  \Vert
\quad\text{when $d \sim n \log(n) / \varepsilon^2$.}$$ See . Although
the residuals are comparable, it need *not* be the case that
$\widehat{\bm{\mathsf{x}}} \approx \bm{\mathsf{x}}_{\star}$, even when
the solution
to [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"} is unique.

The sketch-and-solve paradigm requires us to form the matrix
$\bm{\mathsf{SA}}$ at a cost of $O(mn \log d)$ operations. We would
typically solve the (dense) reduced problem with a direct method, using
$O(dn^2)$ operations. Assuming $d \sim n \log(n)/\varepsilon^2$, the
total arithmetic cost is
$O(mn \log(n/\varepsilon^2) + n^3 \log(n)/\varepsilon^2)$.

In summary, we witness an improvement in computational cost over
classical methods if $\log n \ll n \ll m/\log n$ and $\varepsilon$ is
constant. But we must also be willing to accept large errors, because we
cannot make $\varepsilon$ small.

::: remark
**Remark 28** (History). *The sketch-and-solve paradigm is attributed to
. It plays a major role in the theoretical algorithms literature; see 
for advocacy. It has also been proposed for enormous problems that might
otherwise be entirely hopeless [@2017_weare_randomized_iteration].*
:::

## Iterative sketching

Iterative sketching attempts to remediate the poor accuracy of the
sketch-and-solve paradigm by applying it repeatedly to reduce the
residual error.

First, we construct an initial solution
$\bm{\mathsf{x}}_0 \in \mathbb{F}^n$ using the sketch-and-solve paradigm
with a constant distortion embedding. For each iteration $i$, draw a
fresh random subspace embedding
$\bm{\mathsf{S}}_i \in \mathbb{F}^{d\times m}$ for
$\operatorname{range}(\bm{\mathsf{A}})$, with constant distortion. We
can solve a sequence of least-squares problems $$\begin{equation}
 \label{eqn:iterative-sketch-ls}
\underset{\bm{\mathsf{x}} \in \mathbb{R}^n}{\text{minimize}} \
    \frac{1}{2} \Vert  \bm{\mathsf{S}}_i \bm{\mathsf{A}}(\bm{\mathsf{x}} - \bm{\mathsf{x}}_{i-1})  \Vert^2
    + \langle  \bm{\mathsf{x}} - \bm{\mathsf{x}}_{i-1} , \,  \bm{\mathsf{A}}^* (\bm{\mathsf{b}} - \bm{\mathsf{A}}\bm{\mathsf{x}}_{i-1})  \rangle.
\end{equation}$$ The solution $\bm{\mathsf{x}}_i$ to this subproblem is
fed into the next subproblem. Without the sketch, each subproblem is
equivalent to
solving [\[eqn:overdet-ls2\]](#eqn:overdet-ls2){reference-type="eqref"
reference="eqn:overdet-ls2"} with $\bm{\mathsf{b}}$ replaced by the
residual
$\bm{\mathsf{r}}_{i-1} = \bm{\mathsf{b}} - \bm{\mathsf{A}} \bm{\mathsf{x}}_{i-1}$.
The sketch $\bm{\mathsf{S}}_i$ preserves the geometry while reducing the
problem size. A typical theoretical error bound would be
$$\Vert  \bm{\mathsf{A}}\bm{\mathsf{x}}_j - \bm{\mathsf{b}}  \Vert \leq (1 + \varepsilon) \Vert \bm{\mathsf{A}}\bm{\mathsf{x}}_{\star} - \bm{\mathsf{b}} \Vert
\quad\text{when $j \sim \log(1/\varepsilon)$ and $d \sim n \log n$.}$$
See [@PW16:Iterative-Hessian] for related results.

In each iteration, the iterative sketching approach requires us to form
$\bm{\mathsf{S}}_i \bm{\mathsf{A}}$ at a cost of $O(mn \log d)$. We
compute $\bm{\mathsf{A}}^* \bm{\mathsf{r}}_{i-1}$ at a cost of $O(mn)$.
Although it is unnecessary to solve each subproblem accurately, we
cannot obtain reliable behavior without using a dense method at a cost
of $O(dn^2)$ per iteration. With the theoretical parameter choices, the
total arithmetic is $O((mn + n^3) \log(n) \log(1/\varepsilon))$ to
achieve relative error $\varepsilon$.

The interesting parameter regime is $\log n \ll n \ll m / \log n$, but
we can now allow $\varepsilon$ to be tiny. In this setting, iterative
sketching costs slightly more than the sketch-and-solve paradigm to
achieve constant relative error, while it is faster than the classical
approach. At the same time, it can produce errors as small as
traditional least-squares algorithms. A shortcoming is that this method
requires repeated sketches of the matrix $\bm{\mathsf{A}}$.

For overdetermined least-squares, we can short-circuit the iterative
sketching approach. In this setting, we can sketch the input matrix just
once and factorize it. We can use the same factorized sketch in each
iteration to solve the subproblems faster. For problems more general
than least-squares, it may be necessary to extract a fresh sketch at
each iteration, as we have done here.

::: remark
**Remark 29** (History). *Iterative sketching can be viewed as an
extension of stochastic approximation methods from optimization, for
example stochastic gradient descent [@Bot10:Large-Scale-Machine]. In the
context of randomized NLA, these algorithms first appeared in the guise
of the randomized Kaczmarz iteration [@SV09:Randomized-Kaczmarz]; see
Section [17.4](#sec:rk){reference-type="ref" reference="sec:rk"}.
reinterpreted randomized Kaczmarz as an iterative sketching method and
developed generalizations. proposed a similar method for solving
overdetermined least-squares problems with constraints; they observed
that better numerical performance is obtained by
sketching [\[eqn:overdet-ls2\]](#eqn:overdet-ls2){reference-type="eqref"
reference="eqn:overdet-ls2"} instead
of [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"}.*
:::

## Sketch and precondition {#sec:sketchandprecond}

The sketch-and-precondition paradigm uses random embedding to find a
proxy for the input matrix. We can use this proxy to precondition a
classical iterative algorithm so it converges in a minimal number of
iterations.

Let $\bm{\mathsf{S}} \in \mathbb{F}^{d \times m}$ be a subspace
embedding for $\operatorname{range}(\bm{\mathsf{A}})$ with constant
distortion. Compress the input matrix $\bm{\mathsf{A}}$, and then
compute a (pivoted) QR factorization:
$$\bm{\mathsf{Y}} = \bm{\mathsf{SA}}
\quad\text{and}\quad
\bm{\mathsf{Y}} = \bm{\mathsf{QR}}.$$ Since $\bm{\mathsf{S}}$ preserves
the range of $\bm{\mathsf{A}}$ when $d \sim n \log n$, we anticipate
that
$\bm{\mathsf{Y}}^* \bm{\mathsf{Y}} \approx \bm{\mathsf{A}}^*\bm{\mathsf{A}}$.
As a consequence, $\bm{\mathsf{A}}\bm{\mathsf{R}}^\dagger$ should be
close to an isometry. Thus, we can pass to the preconditioned problem
$$\begin{equation}
 \label{eqn:sketch-precond-ls}
\underset{\bm{\mathsf{x}} \in \mathbb{R}^n}{\text{minimize}} \
    \frac{1}{2} \Vert  (\bm{\mathsf{A}}\bm{\mathsf{R}}^\dagger)(\bm{\mathsf{R}}\bm{\mathsf{x}}) - \bm{\mathsf{b}}  \Vert^2.
\end{equation}$$ Construct an initial solution
$\bm{\mathsf{x}}_0 \in \mathbb{F}^n$ using the sketch-and-solve paradigm
with the embedding $\bm{\mathsf{S}}$. From this starting point, we
solve [\[eqn:sketch-precond-ls\]](#eqn:sketch-precond-ls){reference-type="eqref"
reference="eqn:sketch-precond-ls"} using preconditioned LQSR. The $j$th
iterate satisfies
$$\Vert  \bm{\mathsf{A}} \bm{\mathsf{x}}_j - \bm{\mathsf{b}}  \Vert \leq (1 + \varepsilon) \Vert  \bm{\mathsf{A}} \bm{\mathsf{x}}_{\star} - \bm{\mathsf{b}}  \Vert
\quad\text{when $j \sim \log(1/\varepsilon)$.} %$$ This statement is a
reinterpretation of the theory in .

The cost of sketching the input matrix and performing the QR
decomposition is $O(mn \log d + dn^2)$. Afterwards, we pay $O(mn)$ for
each iteration of PCG. With the theoretical parameter settings, the
total cost is $O(mn \log(nB/\varepsilon) + n^3 \log n)$ operations.

Once again, the interesting regime is $\log n \ll n \ll m / \log n$, and
the value of $\varepsilon$ can be very small. For overdetermined
least-squares, this approach is faster than both the sketch-and-solve
paradigm and the iterative sketching paradigm. The
sketch-and-precondition approach leads to errors that are comparable
with classical linear algebra algorithms, but it may be a factor of
$n / \log(n)$ faster. On the other hand, it requires repeated
applications of the matrix $\bm{\mathsf{A}}$.

::: remark
**Remark 30** (History). *The randomized preconditioning idea was
proposed by . demonstrate that least-squares algorithms based on
randomized preconditioning can beat the highly engineered software in
LAPACK. The same method drives the algorithms in . contains a recent
summary of existing randomized preconditioning methods.*
:::

## Comparisons

If we seek a high-precision solution to a dense, unstructured,
overdetermined least-squares problem, randomized preconditioning leads
to the most efficient existing algorithm. For the same problem, if we
can only view the input matrix once, then the sketch-and-solve paradigm
still allows us to obtain a low-accuracy solution. Although iterative
sketching is less efficient than its competitors in this setting, it
remains useful for solving constrained least-squares problems, and it
has further connections with optimization.

## Summary

From the perspective of a numerical analyst, randomized preconditioning
and iterative sketching should be the preferred methods for designing
sketching algorithms because they allow for high precision. The
sketch-and-solve approach is appropriate only when data access is
severely constrained.

In spite of this fact, a majority of the literature on randomized NLA
develops algorithms based on the sketch-and-solve paradigm. There are
far fewer works on randomized preconditioning or iterative sketching.
This discrepancy points to an opportunity for further research.

# The randomized rangefinder {#sec:random-rangefinder}

A core challenge in linear algebra is to find a subspace that captures a
lot of the action of a matrix. We call this the *rangefinder* problem.
As motivation for considering this problem, we will use the rangefinder
to derive the randomized SVD algorithm. Then we will introduce several
randomized algorithms for computing the rangefinder primitive, along
with theoretical guarantees for these methods. These algorithms all make
use of random embeddings, but their performance depends on more subtle
features than the basic subspace embedding property.

In Section [12](#sec:error-est){reference-type="ref"
reference="sec:error-est"}, we will complement the algorithmic
discussion with details about error estimation and adaptivity for the
rangefinder primitive. In
Sections [13](#sec:natural){reference-type="ref"
reference="sec:natural"}--[16](#sec:full){reference-type="ref"
reference="sec:full"}, we will see that the subspace produced by the
rangefinder can be used as a primitive for other linear algebra
computations.

Most of the material in this section is adapted from our
papers [@HMT11:Finding-Structure] and [@HMST11:Algorithm-Principal]. We
have also incorporated more recent perspectives.

## The rangefinder: Problem statement {#sec:rrf-overview}

Let $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$ be an input matrix,
and let $\ell \leq \min\{m, n\}$ be the subspace dimension. The goal of
the rangefinder problem is to produce an orthonormal matrix
$\bm{\mathsf{Q}} \in \mathbb{F}^{m \times \ell}$ whose range aligns with
the dominant left singular vectors of $\bm{\mathsf{B}}$.

To measure the quality of $\bm{\mathsf{Q}}$, we use the spectral norm
error $$\begin{equation}
 \label{eqn:rrf-error}
\Vert  \bm{\mathsf{B}} - \bm{\mathsf{QQ}}^* \bm{\mathsf{B}}  \Vert
    = \Vert  (\bm{\mathsf{I}}- \bm{\mathsf{QQ}}^*) \bm{\mathsf{B}}  \Vert.
\end{equation}$$ If the error
measure [\[eqn:rrf-error\]](#eqn:rrf-error){reference-type="eqref"
reference="eqn:rrf-error"} is small, then the rank-$\ell$ matrix
$\widehat{\bm{\mathsf{B}}} = \bm{\mathsf{QQ}}^* \bm{\mathsf{B}}$ can
serve as a proxy for $\bm{\mathsf{B}}$. See
Section [11.2](#sec:rsvd){reference-type="ref" reference="sec:rsvd"} for
an important application.

:::: algorithm
::: algorithmic
Input matrix $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$, subspace
dimension $\ell$ Orthonormal matrix
$\bm{\mathsf{Q}} \in \mathbb{F}^{m \times \ell}$

Draw a random matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times \ell}$ Form
$\bm{\mathsf{Y}} = \bm{\mathsf{B\Omega}}$ Compute
$[\bm{\mathsf{Q}}, \sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}})$
:::
::::

### The randomized rangefinder: A pseudoalgorithm {#sec:rrf-pseudo}

Using randomized methods, it is remarkably easy to find an initial
solution to the rangefinder problem. We simply multiply the target
matrix by a random embedding and then orthogonalize the resulting
matrix.

More rigorously: consider a target matrix
$\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$ and a subspace dimension
$\ell$. We draw a random test matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times \ell}$, where
$\bm{\mathsf{\Omega}}^*$ is a mixing random embedding. We form the
product
$\bm{\mathsf{Y}} = \bm{\mathsf{B}}\bm{\mathsf{\Omega}} \in \mathbb{F}^{m \times \ell}$.
Then we compute an orthonormal basis
$\bm{\mathsf{Q}} \in \mathbb{F}^{m \times \ell}$ for the range of
$\bm{\mathsf{Y}}$ using a QR factorization method. See
Algorithm [\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
reference="alg:random-rangefinder"} for pseudocode.
::: {#alg:random-rangefinder .algorithm}
**Algorithm. The randomized rangefinder.**

```
**Input:** Input matrix $\mtx{B} \in \F^{m \times n}$, subspace dimension $\ell$
**Output:** Orthonormal matrix $\mtx{Q} \in \F^{m \times \ell}$

function RandomRangefinder($\mtx{B}$, $\ell$):
    Draw a random matrix $\mtx{\Omega} \in \F^{n \times \ell}$
    Form $\mtx{Y} = \mtx{B\Omega}$
    Compute $[\mtx{Q}, \sim] = qr_econ(\mtx{Y})$
```
:::



In a general setting, the arithmetic cost of this procedure is dominated
by $O(mn\ell)$ operations for the matrix--matrix multiplication. The QR
factorization of $\bm{\mathsf{Y}}$ requires $O(m\ell^2)$ arithmetic, and
we also need to simulate the $n \times \ell$ random matrix
$\bm{\mathsf{\Omega}}$. Economies are possible when either
$\bm{\mathsf{B}}$ or $\bm{\mathsf{\Omega}}$ admits fast multiplication.

### Practicalities

To implement
Algorithm [\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
reference="alg:random-rangefinder"} effectively, several computational
aspects require attention.

- **How do we choose the subspace dimension?** If we have advance
  knowledge of the "effective rank" $r$ of the target matrix
  $\bm{\mathsf{B}}$, the theory
  (Theorem [35](#thm:rrf-gauss){reference-type="ref"
  reference="thm:rrf-gauss"} and
  Corollary [39](#cor:rrf-power-gauss){reference-type="ref"
  reference="cor:rrf-power-gauss"}) indicates that we can select the
  subspace dimension $\ell$ to be just slightly larger, say,
  $\ell = r + p$ where $p = 5$ or $p = 10$. The value $p$ is called the
  *oversampling* parameter. Alternatively, we can use an error estimator
  to decide when the computed subspace $\bm{\mathsf{Q}}$ is sufficiently
  accurate; see Section [12.1](#sec:rrf-error){reference-type="ref"
  reference="sec:rrf-error"}.

- **What kind of random matrix?** We can use most types of mixing random
  embeddings to implement
  Algorithm [\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
  reference="alg:random-rangefinder"}. We highly recommend Gaussians and
  random partial isometries
  (Section [8](#sec:gauss){reference-type="ref" reference="sec:gauss"}).
  Sparse maps, SRTTs, and tensor random embeddings
  (Section [9](#sec:dimension-reduction){reference-type="ref"
  reference="sec:dimension-reduction"}) also work very well. In
  practice, all these approaches exhibit similar behavior; see
  Section [11.5](#sec:rrf-dimred){reference-type="ref"
  reference="sec:rrf-dimred"} for more discussion. We present analysis
  only for Gaussian dimension reduction because it is both simple and
  precise.

- **Matrix multiplication.** The randomized rangefinder is powerful
  because most of the computation takes place in the matrix
  multiplication step, which is a highly optimized primitive on most
  computer systems. When the target matrix admits fast matrix--vector
  multiplications (e.g., due to sparsity), the rangefinder can exploit
  this property.

- **Powering.** As we will discuss in
  Sections [11.6](#sec:rrf-subspace){reference-type="ref"
  reference="sec:rrf-subspace"}
  and [11.7](#sec:rrf-krylov){reference-type="ref"
  reference="sec:rrf-krylov"}, it is often beneficial to enhance
  Algorithm [\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
  reference="alg:random-rangefinder"} by means of powering or Krylov
  subspace techniques.

- **Orthogonalization.** The columns of the matrix $\bm{\mathsf{Y}}$
  tend to be strongly aligned, so it is important to use a numerically
  stable orthogonalization procedure [@GVL13:Matrix-Computations-4ed
  Chap. 5], such as Householder reflectors, double Gram--Schmidt, or
  rank-revealing QR. The rangefinder algorithm is also a natural place
  to invoke a TSQR algorithm [@DGHL12:Communication-Optimal-Parallel].

See  for much more information.

## The randomized singular value decomposition (RSVD) {#sec:rsvd}

Before we continue with our discussion of the rangefinder, let us
summarize one of the key applications: the randomized SVD algorithm.

Low-rank approximation problems often arise when a user seeks an
incomplete matrix factorization that exposes structure, such as a
truncated eigenvalue decomposition or a partial QR factorization. The
randomized rangefinder, described in Section
[11.1.1](#sec:rrf-pseudo){reference-type="ref"
reference="sec:rrf-pseudo"}, can be used to perform the heavy lifting in
these computations. Afterwards, we perform some light post-processing to
reach the desired factorization.

To illustrate how this works, suppose that we want to compute an
approximate rank-$\ell$ truncated singular value decomposition of the
input matrix $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$. That is,
$$\bm{\mathsf{B}} \approx \bm{\mathsf{U}}\bm{\mathsf{\Sigma}}\bm{\mathsf{V}}^{*},$$
where $\bm{\mathsf{U}} \in \mathbb{F}^{m \times \ell}$ and
$\bm{\mathsf{V}} \in \mathbb{F}^{n \times \ell}$ are orthonormal
matrices and $\bm{\mathsf{\Sigma}} \in \mathbb{F}^{\ell \times \ell}$ is
a diagonal matrix whose diagonal entries approximate the largest
singular values of $\bm{\mathsf{B}}$.

Choose a target rank $\ell$, and suppose that
$\bm{\mathsf{Q}} \in \mathbb{F}^{m \times \ell}$ is a computed solution
to the rangefinder problem. The rangefinder furnishes an approximate
rank-$\ell$ factorization of the input matrix:
$\bm{\mathsf{B}} \approx \bm{\mathsf{Q}}\bigl(\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bigr)$.
To convert this representation into a truncated SVD, we just compute an
economy-size SVD of the matrix
$\bm{\mathsf{C}} := \bm{\mathsf{Q}}^{*}\bm{\mathsf{A}} \in \mathbb{F}^{\ell\times n}$
and consolidate the factors.

In symbols, once $\bm{\mathsf{Q}}$ is available, the computation
proceeds as follows: $$\begin{align*}
\bm{\mathsf{B}} \approx&\ \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{B}}  &&\{\mbox{matrix--matrix multiplication: }\bm{\mathsf{C}} = \bm{\mathsf{Q}}^{*}\bm{\mathsf{B}}\} \\
=&\ \bm{\mathsf{Q}}\bm{\mathsf{C}}  &&\{\mbox{economy-size SVD: }\bm{\mathsf{C}} = \widehat{\bm{\mathsf{U}}}\bm{\mathsf{\Sigma}}\bm{\mathsf{V}}^{*}\} \\
=&\ \bm{\mathsf{Q}}\widehat{\bm{\mathsf{U}}}\bm{\mathsf{\Sigma}}\bm{\mathsf{V}}^{*}  &&\{\mbox{matrix--matrix multiplication: }\bm{\mathsf{U}} = \bm{\mathsf{Q}}\widehat{\bm{\mathsf{U}}}\} \\
=&\ \bm{\mathsf{U}}\bm{\mathsf{\Sigma}}\bm{\mathsf{V}}^{*}.
\end{align*}$$ After the rangefinder step, the remaining computations
are all exact (modulo floating-point arithmetic errors). Therefore,
$$\|\bm{\mathsf{B}} - \bm{\mathsf{U}}\bm{\mathsf{\Sigma}}\bm{\mathsf{V}}^{*}\| = \|\bm{\mathsf{B}} - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{B}}\|.$$
In words, the accuracy of the approximate SVD is determined entirely by
the error in the rangefinder computation!

Empirically, the smallest singular values and singular vectors of the
approximate SVD contribute to the accuracy of the approximation, but
they are not good estimates for the true singular values and vectors of
the matrix. Therefore, it can be valuable to truncate the rank by
zeroing out the smallest computed singular values. We omit the details.
See [@HMT11:Finding-Structure; @Gu15:Subspace-Iteration] and
[@TYUC19:Streaming-Low-Rank] for further discussion.

Algorithm [\[alg:rsvd\]](#alg:rsvd){reference-type="ref"
reference="alg:rsvd"} contains pseudocode for the randomized SVD. In a
::: {#alg:rsvd .algorithm}
**Algorithm. Randomized singular value decomposition (RSVD).**

```
**Input:** Input matrix $\mtx{B} \in \F^{m \times n}$, factorization rank $\ell$
**Output:** Orthonormal matrices $\mtx{U} \in \F^{m \times \ell}$, $\mtx{V} \in \F^{n \times \ell}$ and a diagonal matrix $\mtx{\Sigma} \in \F^{\ell \times \ell}$
such that $\mtx{B} \approx \mtx{U}\mtx{\Sigma}\mtx{V}^{*}$.

function RSVD($\mtx{B}$, $\ell$):
    $\mtx{Q} = RandomRangefinder(\mtx{B}, \ell)$
    # Algorithm \ref{alg:random-rangefinder}
    $\mtx{C} = \mtx{Q}^{*}\mtx{B}$
    $[\widehat{\mtx{U}}, \mtx{\Sigma}, \mtx{V}] = svd_econ(\mtx{C})$
    $\mtx{U} = \mtx{Q}\widehat{\mtx{U}}$
    [optional] Truncate the factorization to rank $r \leq \ell$
```
:::


general setting, the dominant cost after the rangefinder step is the
matrix--matrix multiply, which requires $O(mn\ell)$ operations. The
storage requirements are $O((m+n)\ell)$ numbers.

:::: algorithm
::: algorithmic
Input matrix $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$,
factorization rank $\ell$ Orthonormal matrices
$\bm{\mathsf{U}} \in \mathbb{F}^{m \times \ell}$,
$\bm{\mathsf{V}} \in \mathbb{F}^{n \times \ell}$ and a diagonal matrix
$\bm{\mathsf{\Sigma}} \in \mathbb{F}^{\ell \times \ell}$ such that
$\bm{\mathsf{B}} \approx \bm{\mathsf{U}}\bm{\mathsf{\Sigma}}\bm{\mathsf{V}}^{*}$.

$\bm{\mathsf{Q}} = \textsc{RandomRangefinder}(\bm{\mathsf{B}}, \ell)$
Algorithm
[\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
reference="alg:random-rangefinder"}
$\bm{\mathsf{C}} = \bm{\mathsf{Q}}^{*}\bm{\mathsf{B}}$
$[\widehat{\bm{\mathsf{U}}}, \bm{\mathsf{\Sigma}}, \bm{\mathsf{V}}] = \texttt{svd\_econ}(\bm{\mathsf{C}})$
$\bm{\mathsf{U}} = \bm{\mathsf{Q}}\widehat{\bm{\mathsf{U}}}$ Truncate
the factorization to rank $r \leq \ell$
:::
::::

The rangefinder primitive allows us to perform other matrix computations
as well. For example, in Section [13](#sec:natural){reference-type="ref"
reference="sec:natural"}, we explain how to use the rangefinder to
construct matrix factorizations where a subset of the rows/columns are
picked to form a basis for the row/column spaces. This approach gives
far better results than the more obvious randomized algorithms based on
coordinate sampling.

## The rangefinder and Schur complements

Why does the randomized rangefinder work? We will demonstrate that the
procedure has its most natural expression in the language of Schur
complements. This point is implicit in the analysis in , and it
occasionally appears more overtly in the literature (e.g., in  and ).
Nevertheless, this connection has not been explored in a systematic way.

::: {#prop:rrf-schur .proposition}
**Proposition 31** (Rangefinder: Schur complements). *Let
$\bm{\mathsf{Y}} = \bm{\mathsf{BX}}$ for an arbitrary test matrix
$\bm{\mathsf{X}} \in \mathbb{F}^{n \times \ell}$, and let
$\bm{\mathsf{P}}_{\bm{\mathsf{Y}}}$ be the orthogonal projector onto the
range of $\bm{\mathsf{Y}}$. Define the approximation error as
$$\bm{\mathsf{E}} := \bm{\mathsf{E}}(\bm{\mathsf{B}}, \bm{\mathsf{X}}) := (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}}.$$
Then the squared error can be written as a Schur
complement [\[eqn:schur-complement\]](#eqn:schur-complement){reference-type="eqref"
reference="eqn:schur-complement"}:
$$\vert \bm{\mathsf{E}} \vert^2 := \bm{\mathsf{E}}^* \bm{\mathsf{E}} %
    = (\bm{\mathsf{B}}^* \bm{\mathsf{B}}) / \bm{\mathsf{X}}.$$ We
emphasize that $\vert \bm{\mathsf{E}} \vert^2$ is a psd matrix, not a
scalar.*
:::

::: proof
*Proof.* This result follows from a short calculation. We can write the
orthogonal projector $\bm{\mathsf{P}}_{\bm{\mathsf{Y}}}$ in the form
$$\bm{\mathsf{P}}_{\bm{\mathsf{Y}}} = (\bm{\mathsf{BX}})((\bm{\mathsf{B}}\bm{\mathsf{X}})^* (\bm{\mathsf{BX}}))^\dagger(\bm{\mathsf{BX}})^*.$$
Abbreviating $\bm{\mathsf{A}} = \bm{\mathsf{B}}^* \bm{\mathsf{B}}$, we
have
$$\bm{\mathsf{E}}^* \bm{\mathsf{E}} = \bm{\mathsf{B}}^* (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}}
    = \bm{\mathsf{A}} - (\bm{\mathsf{A}} \bm{\mathsf{X}})(\bm{\mathsf{X}}^* \bm{\mathsf{A}} \bm{\mathsf{X}})^\dagger(\bm{\mathsf{AX}})^*.$$
This is precisely the
definition [\[eqn:schur-complement\]](#eqn:schur-complement){reference-type="eqref"
reference="eqn:schur-complement"} of the Schur complement
$\bm{\mathsf{A}} / \bm{\mathsf{X}}$. ◻
:::

Proposition [31](#prop:rrf-schur){reference-type="ref"
reference="prop:rrf-schur"} gives us access to the deep theory of Schur
complements [@Zha05:Schur-Complement]. In particular, we have a
beautiful monotonicity property that follows instantly from .

::: {#cor:rrf-monotone .corollary}
**Corollary 32** (Monotonicity). *Suppose that
$\bm{\mathsf{B}}^* \bm{\mathsf{B}} \preccurlyeq\bm{\mathsf{C}}^* \bm{\mathsf{C}}$
with respect to the semidefinite order $\preccurlyeq$. For each fixed
test matrix $\bm{\mathsf{X}}$, $$\begin{aligned}
\vert  \bm{\mathsf{E}}(\bm{\mathsf{B}}, \bm{\mathsf{X}})  \vert^2
%
    &= (\bm{\mathsf{B}}^* \bm{\mathsf{B}}) / \bm{\mathsf{X}} \\
    &\preccurlyeq(\bm{\mathsf{C}}^* \bm{\mathsf{C}}) / \bm{\mathsf{X}}
    = \vert  \bm{\mathsf{E}}(\bm{\mathsf{C}}, \bm{\mathsf{X}})  \vert^2
%
\end{aligned}$$*
:::

In particular, the error increases if we increase any singular value of
$\bm{\mathsf{B}}$ while retaining the same right singular vectors; the
error decreases if we decrease any singular value of $\bm{\mathsf{B}}$.
The left singular vectors do not play a role here. This observation
allows us to identify which target matrices are hardest to approximate.

::: {#ex:extremal .example}
**Example 33** (Extremals). *Consider the parameterized matrix
$\bm{\mathsf{B}}(\bm{\mathsf{\sigma}}) = \bm{\mathsf{U}} \operatorname{diag}(\bm{\mathsf{\sigma}}) \bm{\mathsf{V}}^* \in \mathbb{F}^{n \times n}$,
where $\bm{\mathsf{U}}, \bm{\mathsf{V}}$ are unitary. Suppose that we
fix $\sigma_1$ and $\sigma_{k+1}$. For each test matrix
$\bm{\mathsf{X}}$, the error
$\vert \bm{\mathsf{E}}(\bm{\mathsf{B}}(\bm{\mathsf{\sigma}}), \bm{\mathsf{X}}) \vert^2$
is maximal in the semidefinite order when
$$\bm{\mathsf{\sigma}} = (\underbrace{\sigma_1, \dots, \sigma_1}_k, \underbrace{\sigma_{k+1}, \dots, \sigma_{k+1}}_{n-k}).$$*
:::

It has long been appreciated that
Example [33](#ex:extremal){reference-type="ref" reference="ex:extremal"}
is the hardest matrix to approximate; cf. . The justification of this
insight is new.

## A priori error bounds

Proposition [31](#prop:rrf-schur){reference-type="ref"
reference="prop:rrf-schur"} shows that the error in the rangefinder
procedure can be written as a Schur complement. Incredibly, the Schur
complement of a psd matrix with respect to a random subspace tends to be
quite small. In this section, we summarize a theoretical analysis, due
to , that explains why this claim is true.

### Master error bound

First, we present a deterministic upper bound on the error incurred by
the rangefinder procedure. This requires some notation.

Without loss of generality, we may assume that $m = n$ by extending
$\bm{\mathsf{B}}$ with zeros. For any $k \leq \ell$, construct a
partitioned SVD of the target matrix:
$$\bm{\mathsf{B}} = \bm{\mathsf{U}} %
    \begin{bmatrix} \bm{\mathsf{\Sigma}}_1 & \\ & \bm{\mathsf{\Sigma}}_{2} \end{bmatrix} %
    \begin{bmatrix} \bm{\mathsf{V}}_1 & \bm{\mathsf{V}}_{2} \end{bmatrix}^*
    \quad\text{with $\bm{\mathsf{\Sigma}}_1 \in \mathbb{R}^{k \times k}$ and
    $\bm{\mathsf{V}}_1 \in \mathbb{F}^{n \times k}$.}$$ The factors
$\bm{\mathsf{U}}, \bm{\mathsf{\Sigma}}, \bm{\mathsf{V}}$ are all square
matrices. As usual, the entries of
$\bm{\mathsf{\Sigma}} = \operatorname{diag}(\sigma_1, \sigma_2, \dots)$
are arranged in weakly decreasing order. So $\bm{\mathsf{\Sigma}}_1$
lists the first $k$ singular values, and $\bm{\mathsf{\Sigma}}_2$ lists
the remaining $n - k$ singular values. The matrix $\bm{\mathsf{V}}_1$
contains the first $k$ right singular vectors; the matrix
$\bm{\mathsf{V}}_2$ contains the remaining $n - k$ right singular
vectors.

For any test matrix $\bm{\mathsf{X}} \in \mathbb{F}^{n \times \ell}$,
define $$\bm{\mathsf{X}}_1 = \bm{\mathsf{V}}_1^* \bm{\mathsf{X}}
\quad\text{and}\quad
\bm{\mathsf{X}}_{2} = \bm{\mathsf{V}}_{2}^* \bm{\mathsf{X}}.$$ These
matrices reflect the alignment of the test matrix $\bm{\mathsf{X}}$ with
the matrix $\bm{\mathsf{V}}_1$ of dominant right singular vectors of
$\bm{\mathsf{B}}$. We assume $\bm{\mathsf{X}}_1$ has full row rank.

With this notation, we can present a strong deterministic bound on the
error in
Algorithm [\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
reference="alg:random-rangefinder"}.

::: {#thm:rrf-deterministic .theorem}
**Theorem 34** (Rangefinder: Deterministic bound). *Let
$\bm{\mathsf{Y}} = \bm{\mathsf{BX}}$ be the sample matrix obtained by
testing $\bm{\mathsf{B}}$ with $\bm{\mathsf{X}}$. With the notation and
assumptions above, for all $k \leq \ell$, $$\begin{equation}
 \label{eqn:rrf-determ}
\Vert  (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}}  \Vert
    \leq \sigma_{k+1} + \Vert  \bm{\mathsf{\Sigma}}_{2} \bm{\mathsf{X}}_{2} \bm{\mathsf{X}}_1^\dagger \Vert.
%
\end{equation}$$ A related inequality holds for every quadratic
unitarily invariant norm.*
:::

The result and its proof are drawn from . The same bound was obtained
independently in  by means of a different technique.

Theorem [34](#thm:rrf-deterministic){reference-type="ref"
reference="thm:rrf-deterministic"} leads to sharp bounds on the
performance of the rangefinder in most situations of practical interest.
Let us present a sketch of the argument. Our approach can be modified to
obtain matching lower and upper bounds, but they do not give any
additional insight into the performance.

::: proof
*Proof.* In view of
Proposition [31](#prop:rrf-schur){reference-type="ref"
reference="prop:rrf-schur"}, we want to bound the spectral norm of
$$\vert (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}} \vert^2
    = (\bm{\mathsf{B}}^*\bm{\mathsf{B}}) / \bm{\mathsf{X}}.$$ First,
change coordinates so that the right singular vectors of
$\bm{\mathsf{B}}$ are the identity: $\bm{\mathsf{V}} = \bm{\mathsf{I}}$.
In particular,
$\bm{\mathsf{B}}^* \bm{\mathsf{B}} = \bm{\mathsf{\Sigma}}^2$ is
diagonal. By homogeneity
of [\[eqn:rrf-determ\]](#eqn:rrf-determ){reference-type="eqref"
reference="eqn:rrf-determ"}, we may assume that $\sigma_1 = 1$. Next,
using Corollary [32](#cor:rrf-monotone){reference-type="ref"
reference="cor:rrf-monotone"}, we may also assume that
$\sigma_1 = \dots = \sigma_k = 1$, which leads to the worst-case error.
Thus, it suffices to bound the spectral norm of the psd matrix
$$\bm{\mathsf{S}} := \begin{bmatrix} \bm{\mathsf{I}}_k & \bm{\mathsf{0}} \\ \bm{\mathsf{0}} & \bm{\mathsf{\Sigma}}_{2}^2\end{bmatrix} / \bm{\mathsf{X}}.$$

To accomplish this task, we may as well take the Schur complement of the
diagonal matrix with respect to a test matrix that has a smaller range
than $\bm{\mathsf{X}}$; see . Define
$\widetilde{\bm{\mathsf{X}}} := \bm{\mathsf{X}} \bm{\mathsf{X}}_{1}^\dagger= [\bm{\mathsf{I}}_k; \bm{\mathsf{X}}_{2}\bm{\mathsf{X}}_1^\dagger]$.
Since
$\operatorname{range}(\widetilde{\bm{\mathsf{X}}}) \subseteq \operatorname{range}(\bm{\mathsf{X}})$,
$$\bm{\mathsf{S}} \preccurlyeq\begin{bmatrix} \bm{\mathsf{I}}_k & \bm{\mathsf{0}} \\ \bm{\mathsf{0}} & \bm{\mathsf{\Sigma}}_{2}^2\end{bmatrix}  / \widetilde{\bm{\mathsf{X}}}
    =: \widetilde{\bm{\mathsf{S}}}.$$ Using the definition of the Schur
complement, we can write out the matrix on the right-hand side in block
form. With the abbreviation
$\bm{\mathsf{F}} := \bm{\mathsf{\Sigma}}_{2} \bm{\mathsf{X}}_{2} \bm{\mathsf{X}}_1^\dagger$,
$$\widetilde{\bm{\mathsf{S}}}
    = \begin{bmatrix} \bm{\mathsf{I}}- (\bm{\mathsf{I}}+ \bm{\mathsf{FF}}^*)^{-1} & \star \\ \star & \bm{\mathsf{\Sigma}}_{2}^2 - \bm{\mathsf{F}}(\bm{\mathsf{I}}+ \bm{\mathsf{FF}}^*)^{-1} \bm{\mathsf{F}}^* \end{bmatrix}.
%
%$$ The $\star$ symbol denotes matrices that do not play a role in the
rest of the argument. We can bound the block matrix above in the psd
order:
$$\widetilde{\bm{\mathsf{S}}} \preccurlyeq\begin{bmatrix} \bm{\mathsf{FF}}^* & \star \\ \star & \bm{\mathsf{\Sigma}}_{2}^2 \end{bmatrix}$$
The inequality for the top-left block holds because
$1 - (1+a)^{-1} \leq a$ for all numbers $a \geq 0$. Last, take the
spectral norm: $$%
\Vert  \bm{\mathsf{S}}  \Vert
    \leq \Vert  \widetilde{\bm{\mathsf{S}}}  \Vert
    \leq \left\Vert  \begin{bmatrix} \bm{\mathsf{FF}}^* & \star \\ \star & \bm{\mathsf{\Sigma}}_{2}^2 \end{bmatrix}  \right\Vert
    \leq \Vert \bm{\mathsf{FF}}^* \Vert + \Vert \bm{\mathsf{\Sigma}}_{2}^2 \Vert.$$
This bound is stronger than the stated result. ◻
:::

### Gaussian test matrices

We can obtain precise results for the behavior of the randomized
rangefinder when the test matrix is (real) standard normal. Let us
present a variant of .

::: {#thm:rrf-gauss .theorem}
**Theorem 35** (Rangefinder: Gaussian analysis). *Fix a matrix
$\bm{\mathsf{B}} \in \mathbb{R}^{m \times n}$ with singular values
$\sigma_1 \geq \sigma_2 \geq \dots$. Draw a standard normal test matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times \ell}$, and construct the
sample matrix $\bm{\mathsf{Y}} = \bm{\mathsf{B\Omega}}$. Choose
$k < \ell - 1$, and introduce the random variable
$$Z = \Vert  \bm{\mathsf{\Gamma}}^\dagger \Vert
\quad\text{where $\bm{\mathsf{\Gamma}} \in \mathbb{R}^{k \times \ell}$ is standard normal.}$$
Then the expected error in the random rangefinder satisfies
$$\operatorname{\mathbb{E}}\Vert (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}} \Vert
    \leq \left( 1 + \sqrt{\frac{k}{\ell-k-1}} \right) \sigma_{k+1}
        + (\operatorname{\mathbb{E}}Z) \left(\sum\nolimits_{j > k} \sigma_{j}^2\right)^{1/2}.$$*
:::

In other words, the randomized rangefinder computes an
$\ell$-dimensional subspace that captures as much of the action of the
matrix $\bm{\mathsf{B}}$ as the best $k$-dimensional subspace. If we
think about $k$ as fixed and $\ell$ as the variable, we only need to
choose $\ell$ slightly larger than $k$ to enjoy this outcome.

The error is comparable with $\sigma_{k+1}$, the error in the best
rank-$k$ approximation, provided that the tail singular values
$\sigma_{j}$ for $j > k$ have small $\ell_2$ norm. This situation
occurs, for example, when $\bm{\mathsf{B}}$ has a rapidly decaying
spectrum.

::: proof
*Proof.* Here is a sketch of the argument. Since the test matrix
$\bm{\mathsf{\Omega}}$ is standard normal, the matrices
$\bm{\mathsf{\Omega}}_{1} := \bm{\mathsf{V}}_1^* \bm{\mathsf{\Omega}}$
and
$\bm{\mathsf{\Omega}}_{2} := \bm{\mathsf{V}}_2^*\bm{\mathsf{\Omega}}$
are independent standard normal matrices because $\bm{\mathsf{V}}_1$ and
$\bm{\mathsf{V}}_2$ are orthonormal and mutually orthogonal. Using
Chevet's theorem [@HMT11:Finding-Structure Prop. 10.1],
$$\begin{aligned}
\operatorname{\mathbb{E}}\Vert  \bm{\mathsf{\Sigma}}_{2} \bm{\mathsf{\Omega}}_{2} \bm{\mathsf{\Omega}}_1^\dagger \Vert
    &= \operatorname{\mathbb{E}}_{\bm{\mathsf{\Omega}}_{1}} \operatorname{\mathbb{E}}_{\bm{\mathsf{\Omega}}_{2}} \big[ \Vert  \bm{\mathsf{\Sigma}}_{2} \bm{\mathsf{\Omega}}_{2} \bm{\mathsf{\Omega}}_1^\dagger \Vert \big] \\
    &\leq \operatorname{\mathbb{E}}\big[ \Vert \bm{\mathsf{\Sigma}}_{2} \Vert \Vert \bm{\mathsf{\Omega}}_1^\dagger \Vert_{\mathrm{F}} + \Vert \bm{\mathsf{\Sigma}}_2 \Vert_{\mathrm{F}} \Vert \bm{\mathsf{\Omega}}_1^\dagger \Vert \big] \\
    &\leq \sqrt{\frac{k}{\ell-k-1}} \Vert  \bm{\mathsf{\Sigma}}_{2}  \Vert + (\operatorname{\mathbb{E}}Z) \Vert \bm{\mathsf{\Sigma}}_{2} \Vert_{\mathrm{F}}.
\end{aligned}$$ The last inequality involves a well-known estimate for
the trace of an inverted Wishart matrix [@HMT11:Finding-Structure
Prop. 10.2]. ◻
:::

To make use of the result, we simply insert estimates for the
expectation of the random variable $Z$. For instance, $$\begin{aligned}
\operatorname{\mathbb{E}}Z \leq \frac{\mathrm{e}\sqrt{\ell}}{\ell - k}
\quad\text{when $2 \leq k < \ell$}
\quad\text{and}\quad
\operatorname{\mathbb{E}}Z \approx \frac{1}{\sqrt{\ell} - \sqrt{k}}
\quad\text{for $k \ll \ell$.}
\end{aligned}$$ These estimates lead to very accurate performance bounds
across a wide selection of matrices and parameters.

The rangefinder also operates in the regime
$k \in \{ \ell - 1, \ell \}$. In this case, it attains significantly
larger errors. A heuristic is
$$\Vert (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}} \Vert
    \lessapprox (1 + k) \sigma_{k+1} + \sqrt{k} \left( \sum\nolimits_{j > k} \sigma_{j}^2 \right)^{1/2}.$$
This point follows because
$\Vert \bm{\mathsf{\Omega}}_1^\dagger \Vert_{\mathrm{F}} \approx k$ and
$\Vert \bm{\mathsf{\Omega}}_1^\dagger \Vert \approx \sqrt{k}$ when
$k \approx \ell$.

For relevant results about Gaussian random matrices, we refer the reader
to [@Ede89:Eigenvalues-Condition; @DS01:Local-Operator; @CD05:Condition-Numbers; @BS10:Spectral-Analysis]
and [@HMT11:Finding-Structure].

## Other test matrices {#sec:rrf-dimred}

In many cases, it is too expensive to use Gaussian test matrices to
implement
Algorithm [\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
reference="alg:random-rangefinder"}. Instead, we may prefer to apply
(the adjoint of) one of the structured random embeddings discussed in
Section [9](#sec:dimension-reduction){reference-type="ref"
reference="sec:dimension-reduction"}.

### Random embeddings for the rangefinder

Good alternatives to Gaussian test matrices include the following.

- **Sparse maps.** Sparse dimension reduction maps work well in the
  rangefinder procedure, even if the input matrix is sparse. The primary
  shortcoming is the need to use sparse data structures and arithmetic.
  See Section [9.2](#sec:sparse-map){reference-type="ref"
  reference="sec:sparse-map"}.

- **SRTTs.** In practice, subsampled randomized trigonometric transforms
  perform slightly better than Gaussian maps. The main difficulty is
  that the implementation requires fast trigonometric transforms. See
  Section [9.3](#sec:srtt){reference-type="ref" reference="sec:srtt"}.

- **Tensor product maps.** Emerging evidence suggests that tensor
  product random projections are also effective in practice. See
  Section [9.4](#sec:trp){reference-type="ref" reference="sec:trp"}.

Some authors have proposed using random coordinate sampling to solve the
rangefinder problem. We cannot recommend this approach unless it is
impossible to use one of the random embeddings described above. See
Section [9.6](#sec:coord-embed){reference-type="ref"
reference="sec:coord-embed"} for a discussion of random coordinate
sampling and situations where it may be appropriate.

### Universality

In practice, if the test matrix is a mixing random embedding, the error
in the rangefinder is somewhat insensitive to the precise distribution
of the test matrix. In this case, we can use the Gaussian theory to
obtain good heuristics about the performance of other types of
embeddings. Regardless, we always recommend using *a posteriori* error
estimates to validate the performance of the rangefinder method, as well
as downstream matrix approximations; see
Section [12](#sec:error-est){reference-type="ref"
reference="sec:error-est"}.

### Aside: Subspace embeddings

If we merely assume that the test matrix is a subspace embedding, then
we can still perform a theoretical analysis of the rangefinder
algorithm. Here is a typical result, adapted from .

::: {#thm:rrf-srtt .theorem}
**Theorem 36** (Rangefinder: SRTT). *Fix a matrix
$\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$ with singular values
$\sigma_1 \geq \sigma_2 \geq \dots$. Choose a natural number $k$, and
draw an SRTT embedding matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times \ell}$ where
$$\ell \geq 8 (k + 8 \log(kn)) \log k.$$ Construct the sample matrix
$\bm{\mathsf{Y}} = \bm{\mathsf{B\Omega}}$. Then
$$\Vert  (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}}  \Vert \leq (1 + 3 \sqrt{n / \ell}) \cdot \sigma_{k+1},$$
with failure probability at most $O(k^{-1})$.*
:::

::: proof
*Proof.* (Sketch) By
Theorem [34](#thm:rrf-deterministic){reference-type="ref"
reference="thm:rrf-deterministic"},
$$\Vert  (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}}  \Vert
    \leq \left[ 1 + \Vert  \bm{\mathsf{\Omega}}_{2}  \Vert \Vert  \bm{\mathsf{\Omega}}_1^\dagger \Vert \right] \sigma_{k+1}.$$
With the specified choice of $\ell$, the test matrix
$\bm{\mathsf{\Omega}}$ is likely to be an oblivious subspace embedding
of the $k$-dimensional subspace $\bm{\mathsf{V}}_1$ with distortion
$1/3$. Thus, the matrix $\bm{\mathsf{\Omega}}_1^{\dagger}$ has spectral
norm bounded by $3$. The matrix $\sqrt{\ell/n} \, \bm{\mathsf{\Omega}}$
is orthonormal, so the spectral norm of $\bm{\mathsf{\Omega}}_{2}$ is
bounded by $\sqrt{n/\ell}$. ◻
:::

The lower bound in the subspace embedding
property [\[eqn:subspace-embedding\]](#eqn:subspace-embedding){reference-type="eqref"
reference="eqn:subspace-embedding"} is the primary fact about the SRTT
used in the proof. But this is only part of the reason that the
rangefinder works. Accordingly, the outcome of this "soft" analysis is
qualitatively weaker than the "hard" analysis in
Theorem [35](#thm:rrf-gauss){reference-type="ref"
reference="thm:rrf-gauss"}. The resulting bound does not explain the
actual (excellent) performance of
Algorithm [\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
reference="alg:random-rangefinder"} when implemented with an SRTT.

## Subspace iteration {#sec:rrf-subspace}

Theorem [35](#thm:rrf-gauss){reference-type="ref"
reference="thm:rrf-gauss"} shows that the basic randomized rangefinder
procedure,
Algorithm [\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
reference="alg:random-rangefinder"}, can be effective for target
matrices $\bm{\mathsf{B}}$ with a rapidly decaying spectrum.
Nevertheless, in many applications, we encounter matrices that do not
meet this criterion. As in the case of spectral norm estimation
(Section [6](#sec:max-eig){reference-type="ref"
reference="sec:max-eig"}), we can resolve the problem by powering the
matrix.

### Rangefinder with powering {#sec:poweringq}

Let $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$ be a fixed input
matrix, and let $q$ be a natural number. Let
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{m \times \ell}$ be a random test
matrix. We form the sample matrix
$$\bm{\mathsf{Y}} = (\bm{\mathsf{BB}}^*)^q \bm{\mathsf{\Omega}}$$ by
repeated multiplication. Then we compute an orthobasis
$\bm{\mathsf{Q}} \in \mathbb{F}^{m \times \ell}$ for the range of
$\bm{\mathsf{Y}}$ using a QR factorization method.

See
Algorithm [\[alg:power-rangefinder\]](#alg:power-rangefinder){reference-type="ref"
reference="alg:power-rangefinder"} for pseudocode. In general, the
::: {#alg:power-rangefinder .algorithm}
**Algorithm. The powered randomized rangefinder.**

```
**Input:** Input matrix $\mtx{B} \in \F^{m \times n}$, target rank $\ell$, depth $q$.
**Output:** Orthonormal matrix $\mtx{Q} \in \F^{m \times \ell}$

function PowerRangefinder($\mtx{B}$, $\ell$, $q$):
    Draw a random matrix $\mtx{\Omega} \in \F^{m \times \ell}$
    $\mtx{Y}_0 = \mtx{\Omega}$
    for $i = 1, \dots, q$:
        $[\mtx{Y}_{i-1}, \sim] = qr_econ(\mtx{Y}_{i-1})$
        $\mtx{Y}_i = \mtx{B} (\mtx{B}^* \mtx{Y}_{i-1})$
    $[\mtx{Q}, \sim] = qr_econ(\mtx{Y}_q)$
```
:::


arithmetic cost is dominated by the $O(qmn\ell)$ cost of the
matrix--matrix multiplications. Economies are possible when
$\bm{\mathsf{B}}$ admits fast multiplication. Unfortunately, when
powering is used, it is not possible to fundamentally accelerate the
computation by using a structured test matrix $\bm{\mathsf{\Omega}}$.

Algorithm [\[alg:power-rangefinder\]](#alg:power-rangefinder){reference-type="ref"
reference="alg:power-rangefinder"} coincides with the classic subspace
iteration algorithm with a random start. Historically, subspace
iteration was regarded as a method for spectral computations. The block
size $\ell$ was often chosen to be quite small, say $\ell = 3$ or
$\ell = 4$, because the intention was simply to resolve singular values
with multiplicity greater than one.

The randomized NLA literature contains several new insights about the
behavior of randomized subspace iteration. It is now recognized that
*iteration is not required*. In practice, $q = 2$ or $q = 3$ is entirely
adequate to solve the rangefinder problem to fairly high accuracy. The
modern perspective also emphasizes the value of running subspace
iteration with a very large block size $\ell$ to obtain matrix
approximations.

::: remark
**Remark 37** (History). *introduced the idea of using randomized
subspace iteration to obtain matrix approximations by selecting a large
block size $\ell$ and a small power $q$. refactored and simplified the
algorithm, and they presented a complete theoretical justification for
the approach. Subsequent analysis appears in .*
:::

### Analysis

The analysis of the powered rangefinder is an easy consequence of the
following lemma [@HMT11:Finding-Structure Prop. 8.6].

::: {#lem:powering .lemma}
**Lemma 38** (Powering). *Let
$\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$ be a fixed matrix, and let
$\bm{\mathsf{P}} \in \mathbb{F}^{m \times m}$ be an orthogonal
projector. For any number $q \geq 1$,
$$\Vert  (\bm{\mathsf{I}}- \bm{\mathsf{P}}) \bm{\mathsf{B}}  \Vert^{2q}
    \leq \Vert  (\bm{\mathsf{I}}- \bm{\mathsf{P}}) (\bm{\mathsf{BB}}^*)^q  \Vert.$$*
:::

::: proof
*Proof.* This bound follows immediately from the Araki--Lieb--Thirring
inequality . ◻
:::

Theorem [35](#thm:rrf-gauss){reference-type="ref"
reference="thm:rrf-gauss"} gives us bounds for the right-hand side of
the inequality in Lemma [38](#lem:powering){reference-type="ref"
reference="lem:powering"} when $\bm{\mathsf{P}}$ is the orthogonal
projector onto the subspace generated by the powered rangefinder.

::: {#cor:rrf-power-gauss .corollary}
**Corollary 39** (Powered rangefinder: Gaussian analysis). *Under the
same conditions as Theorem [35](#thm:rrf-gauss){reference-type="ref"
reference="thm:rrf-gauss"}, let
$\bm{\mathsf{Y}} = (\bm{\mathsf{BB}}^*)^q \bm{\mathsf{\Omega}}$ be the
sample matrix computed by
Algorithm [\[alg:power-rangefinder\]](#alg:power-rangefinder){reference-type="ref"
reference="alg:power-rangefinder"}. Then $$\begin{multline*}
\operatorname{\mathbb{E}}\Vert  (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}}  \Vert
    \leq \left( \operatorname{\mathbb{E}}\Vert  (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}}  \Vert^{2q} \right)^{1/(2q)} \\
    \leq \left[ \left(1 + \sqrt{\frac{k}{\ell-k-1}} \right) \sigma_{k+1}^{2q}
        + (\operatorname{\mathbb{E}}Z) \left( \sum\nolimits_{j > k} \sigma_j^{4q} \right)^{1/2} \right]^{1/(2q)}.
\end{multline*}$$*
:::

In other words, powering the matrix $\bm{\mathsf{B}}$ drives the error
in the rangefinder to $\sigma_{k+1}$ exponentially fast as the parameter
$q$ increases. In cases where the target matrix has some spectral decay,
it suffices to take $q = 2$ or $q = 3$ to achieve satisfactory results.
For matrices with a flat spectral tail, however, we may need to set
$q \approx \log \min\{m,n\}$ to make the error a constant multiple of
$\sigma_{k+1}$.

:::: algorithm
::: algorithmic
Input matrix $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$, target rank
$\ell$, depth $q$. Orthonormal matrix
$\bm{\mathsf{Q}} \in \mathbb{F}^{m \times \ell}$

Draw a random matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{m \times \ell}$
$\bm{\mathsf{Y}}_0 = \bm{\mathsf{\Omega}}$
$[\bm{\mathsf{Y}}_{i-1}, \sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}}_{i-1})$
$\bm{\mathsf{Y}}_i = \bm{\mathsf{B}} (\bm{\mathsf{B}}^* \bm{\mathsf{Y}}_{i-1})$
$[\bm{\mathsf{Q}}, \sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}}_q)$
:::
::::

## Block Krylov methods {#sec:rrf-krylov}

As in the case of spectral norm estimation
(Section [6](#sec:max-eig){reference-type="ref"
reference="sec:max-eig"}), we can achieve more accurate results with
Krylov subspace methods. Nevertheless, this improvement comes at the
cost of additional storage and more complicated algorithms.

### Rangefinder with a Krylov subspace

Let $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$ be a fixed input
matrix, and let $q$ be a natural number. Let
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{m \times \ell}$ be a random test
matrix. We can form the extended sample matrix
$$\bm{\mathsf{Y}} = \begin{bmatrix} \bm{\mathsf{\Omega}} & (\bm{\mathsf{B}}^*\bm{\mathsf{B}}) \bm{\mathsf{\Omega}} & \dots & (\bm{\mathsf{B}}^*\bm{\mathsf{B}})^q \bm{\mathsf{\Omega}} \end{bmatrix}.$$
Then we compute an orthobasis
$\bm{\mathsf{Q}} \in \mathbb{F}^{m \times (q+1)\ell}$ using a QR
factorization method.

See
Algorithm [\[alg:krylov-rangefinder\]](#alg:krylov-rangefinder){reference-type="ref"
reference="alg:krylov-rangefinder"} for pseudocode. This dominant source
::: {#alg:krylov-rangefinder .algorithm}
**Algorithm. The Krylov randomized rangefinder.**

```
**Input:** Input matrix $\mtx{B} \in \F^{m \times n}$, target rank $\ell$, depth $q$
**Output:** Orthonormal matrix $\mtx{Q} \in \F^{m \times 2(q+1)\ell}$

function KrylovRangefinder($\mtx{B}$, $\ell$, $q$):
    Draw a random matrix $\mtx{\Omega} \in \F^{n \times \ell}$
    $[\mtx{Y}_0,\sim] = qr_econ(\mtx{\Omega})$
    for $i = 1, \dots, q$:
        $[\mtx{Y}_{i-1}, \sim] = qr_econ(\mtx{Y}_{i-1})$
        $\mtx{Y}_i = \mtx{B} (\mtx{B}^* \mtx{Y}_{i-1})$
    $[\mtx{Q}, \sim] = qr_econ([\mtx{Y}_0, \dots, \mtx{Y}_q])$
```
:::


of arithmetic is the $O(qmn\ell)$ cost of matrix--matrix multiplication.
The QR factorization now requires $O(q^2 \ell^2 m)$ operations, a factor
of $q^2$ more than the powered rangefinder. We also need to store
$O(qm\ell)$ numbers, which is roughly a factor $q$ more than the powered
rangefinder. Nevertheless, there is evidence that we can balance the
values of $\ell$ and $q$ to make the computational cost of the Krylov
method comparable with the cost of the power method---and still achieve
higher accuracy.

Historically, block Lanczos methods were used for spectral computations
and for SVD computations. The block size $\ell$ was typically chosen to
be fairly small, say $\ell = 3$ or $\ell = 4$, with the goal of
resolving singular values with multiplicity greater than one. The depth
$q$ of the iteration was usually chosen to be quite large. There is
theoretical and empirical evidence that this parameter regime is the
most efficient for resolving the largest singular values to high
accuracy [@YGL18:Superlinear-Convergence].

The randomized NLA literature has recognized that there are still
potential advantages to choosing the block size $\ell$ to be very large
and to choose the depth $q$ to be quite small
[@HMST11:Algorithm-Principal; @MM15:Randomized-Block]. This parameter
regime leads to algorithms that are more efficient on modern computer
architectures, and it still works extremely well for matrices with a
modest amount of spectral decay [@Tro18:Analysis-Randomized]. The
contemporary literature also places a greater emphasis on the role of
block Krylov methods for computing matrix approximations.

::: remark
**Remark 40** (History). *Block Lanczos methods, which are an efficient
implementation of the block Krylov method for a symmetric matrix, were
proposed by and by . The extension to the rectangular case appears in .
These algorithms have received renewed attention, beginning with
[@HMST11:Algorithm-Principal]. The theoretical analysis of randomized
block Krylov methods is much more difficult than the analysis of
randomized subspace iteration. See
[@MM15:Randomized-Block; @YGL18:Superlinear-Convergence]
and [@Tro18:Analysis-Randomized] for some results.*
:::

:::: algorithm
::: algorithmic
Input matrix $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$, target rank
$\ell$, depth $q$ Orthonormal matrix
$\bm{\mathsf{Q}} \in \mathbb{F}^{m \times 2(q+1)\ell}$

Draw a random matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times \ell}$
$[\bm{\mathsf{Y}}_0,\sim] = \texttt{qr\_econ}(\bm{\mathsf{\Omega}})$
$[\bm{\mathsf{Y}}_{i-1}, \sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}}_{i-1})$
$\bm{\mathsf{Y}}_i = \bm{\mathsf{B}} (\bm{\mathsf{B}}^* \bm{\mathsf{Y}}_{i-1})$
$[\bm{\mathsf{Q}}, \sim] = \texttt{qr\_econ}([\bm{\mathsf{Y}}_0, \dots, \bm{\mathsf{Y}}_q])$
:::
::::

### Alternative bases

Algorithm [\[alg:krylov-rangefinder\]](#alg:krylov-rangefinder){reference-type="ref"
reference="alg:krylov-rangefinder"} computes a monomial basis for the
Krylov subspace, which is a poor choice numerically. Let us mention two
more practical alternatives.

1.  **Block Lanczos.** The classical approach uses the block Lanczos
    iteration with reorthogonalization to compute a Lanczos-type basis
    for the Krylov subspace.
    Algorithm [\[alg:rrf-lanczos\]](#alg:rrf-lanczos){reference-type="ref"
    reference="alg:rrf-lanczos"} contains pseudocode for this approach
::: {#alg:rrf-lanczos .algorithm}
**Algorithm. The Lanczos randomized rangefinder.**

```
**Input:** Input matrix $\mtx{B} \in \F^{m \times n}$, target rank $\ell$, depth $q$
**Output:** Orthonormal matrix $\mtx{Q} \in \F^{m \times (q+1)\ell}$

function LanczosRangefinder($\mtx{B}$, $\ell$, $q$):
    Draw a random matrix $\mtx{\Omega} \in \F^{m \times \ell}$
    $[\mtx{Q}_0, \sim] = qr_econ(\mtx{\Omega})$
    $\mtx{W}_0 = \mtx{B}\mtx{Q}_0$
    $[\mtx{P}_0, \mtx{R}_1] = qr_econ(\mtx{W}_0)$
    for $i = 1, \dots, q$:
        $\mtx{Z}_i = \mtx{B}^* \mtx{P}_{i-1} - \mtx{Q}_{i-1} \mtx{R}_{2i-1}^*$
        # Lanczos recursion, part 1
        for $j = 0, \dots, i-1$:
            # Double Gram--Schmidt
            $\mtx{Z}_i = \mtx{Z}_i - \mtx{Q}_{j} (\mtx{Q}_{j}^* \mtx{Z}_i)$
            $\mtx{Z}_i = \mtx{Z}_i - \mtx{Q}_{j} (\mtx{Q}_{j}^* \mtx{Z}_i)$
        $[\mtx{Q}_i, \mtx{R}_{2i}] = qr_econ(\mtx{Z}_i)$
        $\mtx{W}_i = \mtx{B}\mtx{Q}_i - \mtx{P}_{i-1} \mtx{R}_{2i}^*$
        # Lanczos recursion, part 2
        for $j = 0, \dots, i-1$:
            # Double Gram--Schmidt
            $\mtx{W}_i = \mtx{W}_i - \mtx{P}_{j} (\mtx{P}_j^* \mtx{W}_i)$
            $\mtx{W}_i = \mtx{W}_i - \mtx{P}_{j} (\mtx{P}_j^* \mtx{W}_i)$
        $[\mtx{P}_i, \mtx{R}_{2i+1}] = qr_econ(\mtx{W}_i)$
    $\mtx{Q} = [ \mtx{Q}_0, \dots, \mtx{Q}_{q} ]$
```
:::


    adapted from . The basis computed by this algorithm has the property
    that the blocks $\bm{\mathsf{Q}}_i \in \mathbb{F}^{m \times \ell}$
    of the sample matrix are mutually orthogonal. The cost is similar to
    the cost of computing the monomial basis, but there are advantages
    when the subspace is used for spectral computations
    (Section [6](#sec:max-eig){reference-type="ref"
    reference="sec:max-eig"}). Indeed, we can approximate the largest
    singular values of $\bm{\mathsf{A}}$ by means of the largest
    singular values of the band matrix
    $$\bm{\mathsf{R}} = \begin{bmatrix}
    \bm{\mathsf{R}}_1 & \bm{\mathsf{R}}_2^* \\
    & \bm{\mathsf{R}}_3 & \bm{\mathsf{R}}_4^* \\
    && \ddots & \ddots \\
    &&& \bm{\mathsf{R}}_{2q-1} & \bm{\mathsf{R}}_{2q}^* \\
    &&&& \bm{\mathsf{R}}_{2q+1}
    \end{bmatrix}.$$ (The notation in this paragraph corresponds to the
    quantities computed in
    Algorithm [\[alg:rrf-lanczos\]](#alg:rrf-lanczos){reference-type="ref"
    reference="alg:rrf-lanczos"}.)

2.  **Chebyshev.** Suppose that we have a good upper bound for the
    spectral norm of the target matrix. (For example, we can obtain one
    using techniques from Section [6](#sec:max-eig){reference-type="ref"
    reference="sec:max-eig"}.) Then we can compute a Chebyshev basis for
    the Krylov subspace. The advantage of this approach is that we can
    postpone all normalization and orthogonalization steps to the end of
    the computation, which is beneficial for distributed computation.
    This approach, which is being presented for the first time, is
    inspired by . See
    Algorithm [\[alg:rrf-chebyshev\]](#alg:rrf-chebyshev){reference-type="ref"
    reference="alg:rrf-chebyshev"} for pseudocode.
::: {#alg:rrf-chebyshev .algorithm}
**Algorithm. The Chebyshev randomized rangefinder.**

```
**Input:** Input matrix $\mtx{B} \in \F^{m \times n}$, target rank $\ell$, depth $q$, and norm bound $\norm{\mtx{B}} \leq \nu$
**Output:** Orthonormal matrix $\mtx{Q} \in \F^{m \times (q+1)\ell}$

function ChebyshevRangefinder($\mtx{B}$, $\ell$, $q$):
    Draw a random matrix $\mtx{\Omega} \in \F^{m \times \ell}$
    $[\mtx{Y}_0, \sim] = qr_econ(\mtx{\Omega})$
    $\mtx{Y}_1 = (2/\nu) \mtx{B} (\mtx{B}^*\mtx{Y}_0) - \mtx{Y}_0$
    for $i = 2, \dots, q$:
        $\mtx{Y}_i = (4/\nu)\mtx{B} (\mtx{B}^*\mtx{Y}_{i-1}) - 2\mtx{Y}_{i-1} - \mtx{Y}_{i-2}$
        # Chebyshev recursion
    $[\mtx{Q}, \sim] = qr_econ([\mtx{Y}_0, \dots, \mtx{Y}_q])$
```
:::



### Analysis

We are not aware of a direct analysis of Krylov subspace methods for
solving the rangefinder problem. One may extract some bounds from
analysis of randomized SVD algorithms. The following result is adapted
from .

::: theorem
**Theorem 41** (Krylov rangefinder: Gaussian analysis). *Under the
conditions in Theorem [35](#thm:rrf-gauss){reference-type="ref"
reference="thm:rrf-gauss"}, let $\bm{\mathsf{Y}}$ be the sample matrix
computed by
Algorithm [\[alg:krylov-rangefinder\]](#alg:krylov-rangefinder){reference-type="ref"
reference="alg:krylov-rangefinder"}. For $0 \leq \varepsilon\leq 1/2$,
$$\operatorname{\mathbb{E}}\Vert  (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{Y}}}) \bm{\mathsf{B}}  \Vert^2
    \leq \left[ 1 + 2\varepsilon+ \frac{9nk(\ell - k)}{\ell - k - 2} \cdot \mathrm{e}^{-4q\sqrt{\varepsilon}} \right] \sigma_{k+1}^2.$$*
:::

In other words, the Krylov method can drive the error bound for the
rangefinder to $O(\varepsilon)$ by using a Krylov subspace with depth
$q \approx \log (n/\varepsilon) / \sqrt{\varepsilon}$. In contrast, the
power method needs about $q \approx (\log n) / \varepsilon$ iterations
to reach the same target. The difference can be very substantial when
$\varepsilon$ is small.

### Spectral computations

The Krylov rangefinder can also be used for highly accurate computation
of the singular values of a general matrix and the eigenvalues of a self
adjoint matrix.
See [@MM15:Randomized-Block; @YGL18:Superlinear-Convergence]
and [@Tro18:Analysis-Randomized] for discussion and analysis.

:::: algorithm
::: algorithmic
Input matrix $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$, target rank
$\ell$, depth $q$ Orthonormal matrix
$\bm{\mathsf{Q}} \in \mathbb{F}^{m \times (q+1)\ell}$

Draw a random matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{m \times \ell}$

$[\bm{\mathsf{Q}}_0, \sim] = \texttt{qr\_econ}(\bm{\mathsf{\Omega}})$
$\bm{\mathsf{W}}_0 = \bm{\mathsf{B}}\bm{\mathsf{Q}}_0$
$[\bm{\mathsf{P}}_0, \bm{\mathsf{R}}_1] = \texttt{qr\_econ}(\bm{\mathsf{W}}_0)$
$\bm{\mathsf{Z}}_i = \bm{\mathsf{B}}^* \bm{\mathsf{P}}_{i-1} - \bm{\mathsf{Q}}_{i-1} \bm{\mathsf{R}}_{2i-1}^*$
Lanczos recursion, part 1 Double Gram--Schmidt
$\bm{\mathsf{Z}}_i = \bm{\mathsf{Z}}_i - \bm{\mathsf{Q}}_{j} (\bm{\mathsf{Q}}_{j}^* \bm{\mathsf{Z}}_i)$
$\bm{\mathsf{Z}}_i = \bm{\mathsf{Z}}_i - \bm{\mathsf{Q}}_{j} (\bm{\mathsf{Q}}_{j}^* \bm{\mathsf{Z}}_i)$
$[\bm{\mathsf{Q}}_i, \bm{\mathsf{R}}_{2i}] = \texttt{qr\_econ}(\bm{\mathsf{Z}}_i)$
$\bm{\mathsf{W}}_i = \bm{\mathsf{B}}\bm{\mathsf{Q}}_i - \bm{\mathsf{P}}_{i-1} \bm{\mathsf{R}}_{2i}^*$
Lanczos recursion, part 2 Double Gram--Schmidt
$\bm{\mathsf{W}}_i = \bm{\mathsf{W}}_i - \bm{\mathsf{P}}_{j} (\bm{\mathsf{P}}_j^* \bm{\mathsf{W}}_i)$
$\bm{\mathsf{W}}_i = \bm{\mathsf{W}}_i - \bm{\mathsf{P}}_{j} (\bm{\mathsf{P}}_j^* \bm{\mathsf{W}}_i)$
$[\bm{\mathsf{P}}_i, \bm{\mathsf{R}}_{2i+1}] = \texttt{qr\_econ}(\bm{\mathsf{W}}_i)$

$\bm{\mathsf{Q}} = [ \bm{\mathsf{Q}}_0, \dots, \bm{\mathsf{Q}}_{q} ]$
:::
::::

:::: algorithm
::: algorithmic
Input matrix $\bm{\mathsf{B}} \in \mathbb{F}^{m \times n}$, target rank
$\ell$, depth $q$, and norm bound $\Vert \bm{\mathsf{B}} \Vert \leq \nu$
Orthonormal matrix $\bm{\mathsf{Q}} \in \mathbb{F}^{m \times (q+1)\ell}$

Draw a random matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{m \times \ell}$

$[\bm{\mathsf{Y}}_0, \sim] = \texttt{qr\_econ}(\bm{\mathsf{\Omega}})$
$\bm{\mathsf{Y}}_1 = (2/\nu) \bm{\mathsf{B}} (\bm{\mathsf{B}}^*\bm{\mathsf{Y}}_0) - \bm{\mathsf{Y}}_0$

$\bm{\mathsf{Y}}_i = (4/\nu)\bm{\mathsf{B}} (\bm{\mathsf{B}}^*\bm{\mathsf{Y}}_{i-1}) - 2\bm{\mathsf{Y}}_{i-1} - \bm{\mathsf{Y}}_{i-2}$
Chebyshev recursion

$[\bm{\mathsf{Q}}, \sim] = \texttt{qr\_econ}([\bm{\mathsf{Y}}_0, \dots, \bm{\mathsf{Y}}_q])$
:::
::::

# Error estimation and adaptivity {#sec:error-est}

The theoretical analysis of the randomized rangefinder in Section
[11](#sec:random-rangefinder){reference-type="ref"
reference="sec:random-rangefinder"} describes with great precision when
the procedure is effective, and what errors to expect. From a practical
point of view, however, the usefulness of this analysis is limited by
the fact that we rarely have advance knowledge of the singular values of
the matrix to be approximated. In this section, we consider the more
typical situation where we are given a matrix
$\bm{\mathsf{A}} \in \mathbb{F}^{m\times n}$ and a tolerance
$\varepsilon$, and it is our job to find a low-rank factorization of
$\bm{\mathsf{A}}$ that is accurate to within precision $\varepsilon$. In
other words, part of the task to be solved is to determine the
$\varepsilon$-rank of $\bm{\mathsf{A}}$.

To complete this job, we must equip the rangefinder with an *a
posteriori* error estimator: given an orthonormal matrix
$\bm{\mathsf{Q}} \in \mathbb{F}^{m\times \ell}$ whose columns form a
putative basis for the range of $\bm{\mathsf{A}}$, it estimates the
corresponding approximation error
$\|\bm{\mathsf{A}} - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\|$.

To solve the fixed-error approximation problem, the idea is to start
with a lowball "guess" at the rank, run the rangefinder, and then check
to see if we are within the requested tolerance. If the answer is
negative, then there are several strategies for how to proceed.
Typically, we would just draw more samples to enrich the basis we
already have on hand. But in some circumstances, it may be better to
start over or to try an alternative approach, such as increasing the
amount of powering that is done.

## A posteriori error estimation {#sec:rrf-error}

Let $\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$ be a target matrix,
and let $\bm{\mathsf{Q}} \in \mathbb{F}^{m \times \ell}$ be an
orthonormal matrix whose columns may or may not form a good basis for
the range of $\bm{\mathsf{A}}$. We typically think of $\bm{\mathsf{Q}}$
as being the output of one of the rangefinder algorithms described in
Section [11](#sec:random-rangefinder){reference-type="ref"
reference="sec:random-rangefinder"}. Our goal is now to produce an
inexpensive and reliable estimate of the error
$\vert\!\vert\!\vert  (\bm{\mathsf{I}}- \bm{\mathsf{QQ}}^*) \bm{\mathsf{A}}  \vert\!\vert\!\vert$
with respect to some norm
$\vert\!\vert\!\vert \cdot \vert\!\vert\!\vert$.

To do so, we draw on the techniques for norm estimation by sampling
(Sections [4](#sec:trace-est){reference-type="ref"
reference="sec:trace-est"}--[5](#sec:schatten-p){reference-type="ref"
reference="sec:schatten-p"}). The basic idea is to collect a (small)
auxiliary sample $$\begin{equation}
\label{eq:auxiliarysample}
\bm{\mathsf{Z}} = \bm{\mathsf{A}} \bm{\mathsf{\Phi}},
\end{equation}$$ where the test matrix
$\bm{\mathsf{\Phi}} \in \mathbb{F}^{n \times s}$ is drawn from a
Gaussian distribution. We assume that $\bm{\mathsf{\Phi}}$ is
statistically independent from whatever process was used to compute
$\bm{\mathsf{Q}}$. Then $$\begin{equation}
\label{eq:errorsample}
(\bm{\mathsf{I}}- \bm{\mathsf{QQ}}^*) \bm{\mathsf{Z}} = (\bm{\mathsf{I}}- \bm{\mathsf{QQ}}^*) \bm{\mathsf{A \Phi}} \in \mathbb{F}^{m \times s}.
\end{equation}$$ is a random sample of the error in the approximation.
We can now use any of the methods from
Sections [4](#sec:trace-est){reference-type="ref"
reference="sec:trace-est"}--[5](#sec:schatten-p){reference-type="ref"
reference="sec:schatten-p"} to estimate the error from this sample. For
example,
$$\Vert  (\bm{\mathsf{I}}- \bm{\mathsf{QQ}}^*) \bm{\mathsf{A}}  \Vert_{\mathrm{F}}^2
    \approx \frac{1}{s} \Vert (\bm{\mathsf{I}}- \bm{\mathsf{QQ}}^*) \bm{\mathsf{Z}} \Vert_{\mathrm{F}}^2.$$
The theory in Section [4](#sec:trace-est){reference-type="ref"
reference="sec:trace-est"} gives precise tail bounds for this estimator.
Similarly, we can approximate the Schatten $4$-norm by computing a
sample variance. Beyond that, the Valiant--Kong estimator
(Section [5.4](#sec:valiant-kong){reference-type="ref"
reference="sec:valiant-kong"}) allows us to approximate higher-order
Schatten norms.

The cost of extracting the auxiliary sample
([\[eq:auxiliarysample\]](#eq:auxiliarysample){reference-type="ref"
reference="eq:auxiliarysample"}) is almost always much smaller than the
cost of running the rangefinder itself; it involves only $s$ additional
matrix--vector multiplications, where $s$ can be thought of as a small
fixed number, say $s=10$.

## A certificate of accuracy for structured random matrices {#sec:certificate}

The idea of drawing an auxiliary sample
([\[eq:auxiliarysample\]](#eq:auxiliarysample){reference-type="ref"
reference="eq:auxiliarysample"}) for purposes of error estimation is
particularly appealing when the rangefinder implementation involves a
structured random test matrix (see Section
[11.5](#sec:rrf-dimred){reference-type="ref"
reference="sec:rrf-dimred"}). These structured random maps can be much
faster than Gaussian random matrices, while producing errors that are
just as small [@HMT11:Finding-Structure Sec. 7.4]. Their main weakness
is that they come with far weaker *a priori* error guarantees; see
Theorem [36](#thm:rrf-srtt){reference-type="ref"
reference="thm:rrf-srtt"}.

Now, consider a situation where we use a structured random matrix to
compute an approximate basis $\bm{\mathsf{Q}}$ for the range of a matrix
$\bm{\mathsf{A}}$ (via Algorithm
[\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
reference="alg:random-rangefinder"}), and a small Gaussian random matrix
$\bm{\mathsf{\Phi}}$ to draw an auxiliary sample
([\[eq:auxiliarysample\]](#eq:auxiliarysample){reference-type="ref"
reference="eq:auxiliarysample"}) from $\bm{\mathsf{A}}$. The additional
cost of extracting the "extra" sample is small, both in terms of
practical execution time and in terms of the asymptotic cost estimates
(which would generally remain unchanged). However, the computed output
can now be relied on with supreme confidence, since it is backed up by
the strong theoretical results that govern Gaussian matrices.

One may even push these ideas further and use the "certificate of
accuracy" to gain confidence in randomized sampling methods based on
heuristics or educated guesses about the matrix being estimated. For
instance, one may observe that matrices that arise in some application
typically have low coherence
(Section [9.6](#sec:coord-embed){reference-type="ref"
reference="sec:coord-embed"}), and one may implement a fast uniform
sampling strategy that only works in this setting. No matter how unsafe
the sampling strategy, we can trust the *a posteriori* error estimator
when it promises that the computed factorization is sufficiently
accurate.

The general idea of observing the action of a residual on random vectors
in order to get an estimate for its magnitude can be traced back at
least as far as . Here, we followed the discussion in and .

::: remark
**Remark 42** (Rank doubling). *Let us consider what should be done if
the *a posteriori* error estimator tells us that a requested tolerance
has not yet been met. Since many structured random matrices have the
unfortunate property that it is not an easy matter to recycle the sample
already computed in order to build a larger one, it often makes sense to
simply start from scratch, but doubling the number of columns in the
test matrix. Such a strategy of doubling the rank at each attempt
typically does not change the order of the dominant term in the
asymptotic cost. It may, however, be somewhat wasteful from a practical
point of view.*
:::

## Adaptive rank determination using Gaussian test matrices {#sec:adaptivebasic}

When a Gaussian test matrix is used, we can incorporate *a posteriori*
error estimation into the rangefinder algorithm with negligible increase
in the amount of computation.

To illustrate, let us first consider a situation where we are given a
matrix $\bm{\mathsf{A}}$, and we use the randomized rangefinder
(Algorithm
[\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
reference="alg:random-rangefinder"}) with a Gaussian random matrix to
find an orthonormal basis for its approximate range. We do not know the
numerical rank of $\bm{\mathsf{A}}$ in advance, so we use a number
$\ell$ of samples that we believe is likely to be more than enough.

In order to include *a posteriori* error estimation, we *conceptually*
split the test matrix so that $$\begin{equation}
\label{eq:courteney1}
\bm{\mathsf{\Omega}} = \bigl[\bm{\mathsf{\Omega}}_{1}\quad\bm{\mathsf{\Omega}}_{2}\bigr],
\end{equation}$$ where $\bm{\mathsf{\Omega}}_{1}$ holds the first
$\ell-s$ columns of $\bm{\mathsf{\Omega}}$. The thin sliver
$\bm{\mathsf{\Omega}}_{2}$ that holds the last $s$ columns will
temporarily play the part of the independent test matrix
$\bm{\mathsf{\Phi}}$. (Note that the matrices $\bm{\mathsf{\Omega}}_i$
defined here are different from those in the proof of
Theorem [35](#thm:rrf-gauss){reference-type="ref"
reference="thm:rrf-gauss"}.)

The sample matrix inherits a corresponding split $$\begin{equation}
\label{eq:courteney2}
\bm{\mathsf{Y}} =
\bm{\mathsf{A}}\bm{\mathsf{\Omega}} =
\bigl[\bm{\mathsf{A}}\bm{\mathsf{\Omega}}_{1}\quad\bm{\mathsf{A}}\bm{\mathsf{\Omega}}_{2}\bigr] =:
\bigl[\bm{\mathsf{Y}}_{1}\quad\bm{\mathsf{Y}}_{2}\bigr].
\end{equation}$$ In order to orthonormalize the columns of
$\bm{\mathsf{Y}}$ to form an approximate basis for the column space, we
perform an unpivoted QR factorization: $$\bm{\mathsf{Y}} =
\bigl[\bm{\mathsf{Q}}_{1}\quad\bm{\mathsf{Q}}_{2}\bigr]\,
\left[\begin{array}{cc}
\bm{\mathsf{R}}_{11} & \bm{\mathsf{R}}_{12} \\
\bm{\mathsf{0}}      & \bm{\mathsf{R}}_{22}
\end{array}\right],$$ where the partitioning conforms with
([\[eq:courteney2\]](#eq:courteney2){reference-type="ref"
reference="eq:courteney2"}). Thus,
$\bm{\mathsf{Y}}_{1} = \bm{\mathsf{Q}}_{1}\bm{\mathsf{R}}_{11}$ and
$\bm{\mathsf{Y}}_{2} = \bm{\mathsf{Q}}_{1}\bm{\mathsf{R}}_{12} + \bm{\mathsf{Q}}_{2}\bm{\mathsf{R}}_{22}$.
Now, observe that $$\bm{\mathsf{Q}}_{2}\bm{\mathsf{R}}_{22} =
\bm{\mathsf{Y}}_{2} - \bm{\mathsf{Q}}_{1}\bm{\mathsf{R}}_{12} =
\bm{\mathsf{Y}}_{2} - \bm{\mathsf{Q}}_{1}\bm{\mathsf{Q}}_{1}^{*}\bm{\mathsf{Y}}_{2} =
\bigl(\bm{\mathsf{I}} - \bm{\mathsf{Q}}_{1}\bm{\mathsf{Q}}_{1}^{*})\bm{\mathsf{A}}\bm{\mathsf{\Omega}}_{2}.$$
In other words, the matrix $\bm{\mathsf{Q}}_{2}\bm{\mathsf{R}}_{22}$ is
a sample of the residual error resulting from using
$\bm{\mathsf{\Omega}}_{1}$ as the test matrix. We can analyze the sample
$\bm{\mathsf{Q}}_{2}\bm{\mathsf{R}}_{22}$ using the techniques described
in Section [12.1](#sec:rrf-error){reference-type="ref"
reference="sec:rrf-error"} to derive an estimate on the norm of the
residual error. If the resulting estimate is small enough, we can
confidently trust the computed factorization. (When the rangefinder is
used as a preliminary step towards computing a partial SVD of the
matrix, we may as well use the full orthonormal basis in
$\bm{\mathsf{Q}}$ in any steps that follow.)

If the error estimate computed is larger than what is acceptable, then
additional work must be done, typically by drawing additional samples to
enrich the basis already computed, as described in the next section.

## An incremental algorithm based on Gaussian test matrices {#sec:randQB}

The basic error estimation procedure outlined in Section
[12.3](#sec:adaptivebasic){reference-type="ref"
reference="sec:adaptivebasic"} is appropriate when it is not onerous to
draw a large number $\ell$ of samples and when the *a posteriori* error
merely serves as an insurance policy against a rare situation where the
singular values decay more slowly than expected. In this section, we
describe a technique that is designed for situations where we have no
notion what the rank may be in advance. The idea is to build the
approximate basis incrementally by drawing and processing one batch of
samples at a time, while monitoring the errors as we go. The upshot is
that this computation can be organized in such a way that the total cost
is essentially the same as it would have been had we known the numerical
rank in advance.

We frame the rangefinder problem as usual: for a given matrix
$\bm{\mathsf{A}}$ and a given tolerance $\tau$, we seek to build an
orthonormal matrix $\bm{\mathsf{Q}}$ such that
$$\|\bm{\mathsf{A}} - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\| \leq \tau.$$
The procedure that we describe will be controlled by a tuning parameter
$b$ that specifies how many columns we process at a time. If we choose
$b$ too large, we may overshoot the numerical rank of $\bm{\mathsf{A}}$
and perform more work than necessary. If $b$ is too small, computational
efficiency may suffer; see Section
[16.2](#sec:blocking){reference-type="ref" reference="sec:blocking"}. In
many environments, picking $b$ between $10$ and $100$ would be about
right.

While the *a posteriori* error estimator signals that the error
tolerance has not been met, the incremental rangefinder successively
draws blocks of $b$ Gaussian random vectors, computes the corresponding
samples, and adds them to the basis. See Algorithm
[\[alg:IncRangeFinder\]](#alg:IncRangeFinder){reference-type="ref"
reference="alg:IncRangeFinder"}. To understand how the method works,
::: {#alg:IncRangeFinder .algorithm}
**Algorithm. Incremental rangefinder.**

```
**Input:** Target matrix $\mtx{A} \in \F^{m \times n}$, tolerance $\tau \in \mathbb{R}_{+}$, block size $b$
**Output:** Orthonormal matrix $\mtx{Q}$ such that $\|\mtx{A} - \mtx{Q}\mtx{Q}^{*}\mtx{A}\| \leq \tau$ with high probability

function IncrementalRangefinder($\mtx{A}$, $\tau$, b):
    $\mtx{Y} = \mtx{A}\mtx{\Omega}$    # Draw $\mtx{\Omega} \in \F^{n\times b}$ from a Gaussian distribution
    $[\mtx{Q}_{1},\sim] = qr_econ(\mtx{Y})$
    $i=1$
    while $norm_est(\mtx{Y}) > \tau$:
        # Norm estimator in Section \ref{sec:rrf-error}
        $i = i+1$
        $\mtx{Y} = \mtx{A}\mtx{\Omega}$    # Draw $\mtx{\Omega} \in \F^{n\times b}$ from a Gaussian distribution
        $\mtx{Y} = \mtx{Y} - \sum_{j=1}^{i-1}\mtx{Q}_{j}\bigl(\mtx{Q}_{j}^{*}\mtx{Y}\bigr)$
        $[\mtx{Q}_{i},\sim] = qr_econ(\mtx{Y})$
    $\mtx{Q} = \bigl[\mtx{Q}_{1}\quad\mtx{Q}_{2}\quad \cdots \quad \mtx{Q}_{i-1}\bigr]$
```
:::


observe that, after line 8 has been executed, the matrix
$\bm{\mathsf{Y}}$ holds the sample $$\begin{equation}
\label{eq:logitech7}
\bm{\mathsf{Y}} = \bigl(\bm{\mathsf{I}}  - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bigr)\bm{\mathsf{A}}\bm{\mathsf{\Omega}}
\end{equation}$$ from the residual
$\bigl(\bm{\mathsf{I}}  - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bigr)\bm{\mathsf{A}}$,
where $\bm{\mathsf{Q}}$ is the cumulative basis that has been built at
that point and where $\bm{\mathsf{\Omega}}$ is an $n\times b$ matrix
drawn from a Gaussian distribution. Since
([\[eq:logitech7\]](#eq:logitech7){reference-type="ref"
reference="eq:logitech7"}) holds, we can estimate
$\|(\bm{\mathsf{I}}  - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bigr)\bm{\mathsf{A}}\|$
using the techniques described in Section
[12.1](#sec:rrf-error){reference-type="ref" reference="sec:rrf-error"}.

:::: algorithm
::: algorithmic
Target matrix $\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$, tolerance
$\tau \in \mathbb{R}_{+}$, block size $b$ Orthonormal matrix
$\bm{\mathsf{Q}}$ such that
$\|\bm{\mathsf{A}} - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\| \leq \tau$
with high probability

$\bm{\mathsf{Y}} = \bm{\mathsf{A}}\bm{\mathsf{\Omega}}$Draw
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n\times b}$ from a Gaussian
distribution
$[\bm{\mathsf{Q}}_{1},\sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}})$ $i=1$
Norm estimator in Section [12.1](#sec:rrf-error){reference-type="ref"
reference="sec:rrf-error"} $i = i+1$
$\bm{\mathsf{Y}} = \bm{\mathsf{A}}\bm{\mathsf{\Omega}}$Draw
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n\times b}$ from a Gaussian
distribution
$\bm{\mathsf{Y}} = \bm{\mathsf{Y}} - \sum_{j=1}^{i-1}\bm{\mathsf{Q}}_{j}\bigl(\bm{\mathsf{Q}}_{j}^{*}\bm{\mathsf{Y}}\bigr)$
[]{#line:Yafter label="line:Yafter"}
$[\bm{\mathsf{Q}}_{i},\sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}})$
$\bm{\mathsf{Q}} = \bigl[\bm{\mathsf{Q}}_{1}\quad\bm{\mathsf{Q}}_{2}\quad \cdots \quad \bm{\mathsf{Q}}_{i-1}\bigr]$
:::
::::

In situations where the matrix $\bm{\mathsf{A}}$ is small enough to fit
in RAM, it often makes sense to explicitly update it after every step.
The benefit to doing so is that one can then determine the norm of the
remainder matrix explicitly. Algorithm
[\[alg:IncRangeFinderUpdating\]](#alg:IncRangeFinderUpdating){reference-type="ref"
reference="alg:IncRangeFinderUpdating"} summarizes the resulting
::: {#alg:IncRangeFinderUpdating .algorithm}
**Algorithm. Incremental rangefinder with updating.**

```
**Input:** Target matrix $\mtx{A} \in \F^{m \times n}$, tolerance $\tau \in \mathbb{R}_{+}$, block size $b$.
**Output:** Orthonormal matrix $\mtx{Q}$ and a matrix $\mtx{B}$ such that $\|\mtx{A} - \mtx{Q}\mtx{B}\| \leq \tau$.

function IncrementalRangefinderWithUpdating($\mtx{A}$, $\tau$, b):
    $\mtx{Y} = \mtx{A}\mtx{\Omega}$    # Draw $\mtx{\Omega} \in \F^{n\times b}$ from a Gaussian distribution.
    $[\mtx{Q}_{1},\sim] = qr_econ(\mtx{Y})$
    $\mtx{B}_{1} = \mtx{Q}_{1}^{*}\mtx{A}$
    $\mtx{A} = \mtx{A} - \mtx{Q}_{1}\mtx{B}_{1}$
    $i=1$
    while $\|\mtx{A}\| > \tau$:
        # Use an inexpensive norm such as Frobenius
        $i           = i+1$
        $\mtx{Y}     = \mtx{A}\mtx{\Omega}$    # Draw $\mtx{\Omega} \in \F^{n\times b}$ from a Gaussian distribution.
        $[\mtx{Q}_{i},\sim] = qr_econ(\mtx{Y})$
        $\mtx{B}_{i} = \mtx{Q}_{i}^{*}\mtx{A}$
        $\mtx{A}     = \mtx{A} - \mtx{Q}_{i}\mtx{B}_{i}$
    $\mtx{Q} = \bigl[\mtx{Q}_{1}\quad\mtx{Q}_{2}\quad \cdots \quad \mtx{Q}_{i}\bigr]$
    $\mtx{B} = \bigl[\mtx{B}_{1}^{*}\quad\mtx{B}_{2}^{*}\quad \cdots \quad \mtx{B}_{i}^{*}\bigr]^{*}$
```
:::


procedure. After line 9 of the algorithm has been executed, the formula
([\[eq:logitech7\]](#eq:logitech7){reference-type="ref"
reference="eq:logitech7"}) holds because, at this point in the
computation, $\bm{\mathsf{A}}$ has been overwritten by
$\bigl(\bm{\mathsf{I}}  - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bigr)\bm{\mathsf{A}}$.
Algorithm
[\[alg:IncRangeFinderUpdating\]](#alg:IncRangeFinderUpdating){reference-type="ref"
reference="alg:IncRangeFinderUpdating"} relates to Algorithm
[\[alg:IncRangeFinder\]](#alg:IncRangeFinder){reference-type="ref"
reference="alg:IncRangeFinder"} in the same way that modified
Gram--Schmidt relates to classical Gram--Schmidt.

:::: algorithm
::: algorithmic
Target matrix $\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$, tolerance
$\tau \in \mathbb{R}_{+}$, block size $b$. Orthonormal matrix
$\bm{\mathsf{Q}}$ and a matrix $\bm{\mathsf{B}}$ such that
$\|\bm{\mathsf{A}} - \bm{\mathsf{Q}}\bm{\mathsf{B}}\| \leq \tau$.

$\bm{\mathsf{Y}} = \bm{\mathsf{A}}\bm{\mathsf{\Omega}}$Draw
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n\times b}$ from a Gaussian
distribution.
$[\bm{\mathsf{Q}}_{1},\sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}})$
$\bm{\mathsf{B}}_{1} = \bm{\mathsf{Q}}_{1}^{*}\bm{\mathsf{A}}$
$\bm{\mathsf{A}} = \bm{\mathsf{A}} - \bm{\mathsf{Q}}_{1}\bm{\mathsf{B}}_{1}$
$i=1$ Use an inexpensive norm such as Frobenius $i           = i+1$
$\bm{\mathsf{Y}}     = \bm{\mathsf{A}}\bm{\mathsf{\Omega}}$Draw
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n\times b}$ from a Gaussian
distribution.
$[\bm{\mathsf{Q}}_{i},\sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}})$
$\bm{\mathsf{B}}_{i} = \bm{\mathsf{Q}}_{i}^{*}\bm{\mathsf{A}}$
$\bm{\mathsf{A}}     = \bm{\mathsf{A}} - \bm{\mathsf{Q}}_{i}\bm{\mathsf{B}}_{i}$
$\bm{\mathsf{Q}} = \bigl[\bm{\mathsf{Q}}_{1}\quad\bm{\mathsf{Q}}_{2}\quad \cdots \quad \bm{\mathsf{Q}}_{i}\bigr]$
$\bm{\mathsf{B}} = \bigl[\bm{\mathsf{B}}_{1}^{*}\quad\bm{\mathsf{B}}_{2}^{*}\quad \cdots \quad \bm{\mathsf{B}}_{i}^{*}\bigr]^{*}$
:::
::::

For matrices whose singular values decay slowly, incorporating a few
steps of power iteration (as described in Section
[11.6](#sec:rrf-subspace){reference-type="ref"
reference="sec:rrf-subspace"}) is very beneficial. In this environment,
it is often necessary to incorporate additional reorthonormalizations to
combat loss of orthogonality due to round-off errors. Full descriptions
of the resulting techniques can be found in .

A version of Algorithm
[\[alg:IncRangeFinderUpdating\]](#alg:IncRangeFinderUpdating){reference-type="ref"
reference="alg:IncRangeFinderUpdating"} suitable for sparse matrices or
matrices stored out-of-core is described in . This variant reorganizes
the computation to avoid the explicit updating step and to reduce the
communication requirements overall. Observe that it is still possible to
evaluate the Frobenius norm of the residual exactly (without randomized
estimation) by using the identity $$\begin{multline*}
\Vert  \bm{\mathsf{A}}  \Vert_{\mathrm{F}}^2 =
\Vert  \bigl(\bm{\mathsf{I}} - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bigr)\bm{\mathsf{A}}  \Vert_{\mathrm{F}}^{2} +
\Vert  \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}  \Vert_{\mathrm{F}}^{2} \\
=
\Vert  \bigl(\bm{\mathsf{I}} - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bigr)\bm{\mathsf{A}}  \Vert_{\mathrm{F}}^{2} +
\Vert  \bm{\mathsf{Q}}\bm{\mathsf{B}}  \Vert_{\mathrm{F}}^{2} =
\Vert  \bigl(\bm{\mathsf{I}} - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bigr)\bm{\mathsf{A}}  \Vert_{\mathrm{F}}^{2} +
\Vert  \bm{\mathsf{B}}  \Vert_{\mathrm{F}}^{2},
\end{multline*}$$ where the first equality holds since the column spaces
of
$\bigl(\bm{\mathsf{I}} - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bigr)\bm{\mathsf{A}}$
and $\bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}$ are orthogonal
and where the third equality holds since $\bm{\mathsf{Q}}$ is
orthonormal. For an error measure, we can use the resulting
relationship:
$$\Vert  \bigl(\bm{\mathsf{I}} - \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bigr)\bm{\mathsf{A}}  \Vert_{\mathrm{F}} =
\sqrt{ \Vert \bm{\mathsf{A}}  \Vert_{\mathrm{F}}^{2} - \Vert \bm{\mathsf{B}}  \Vert_{\mathrm{F}}^{2}},$$
observing that both $\Vert  \bm{\mathsf{A}}  \Vert_{\mathrm{F}}$ and
$\Vert  \bm{\mathsf{B}}  \Vert_{\mathrm{F}}$ can be computed explicitly.

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

# Nyström approximation {#sec:nystrom}

We continue our discussion of matrix approximation with the problem of
finding a low-rank approximation of a positive semidefinite (PSD)
matrix. There is an elegant randomized method for accomplishing this
goal that is related to our solution to the rangefinder problem.

## Low-rank PSD approximation {#sec:nystrom-overview}

Let $\bm{\mathsf{A}} \in \mathbb{H}_n$ be a PSD matrix. For a rank
parameter $k$, the goal is to produce a rank-$k$ PSD matrix
$\widehat{\bm{\mathsf{A}}}_k \in \mathbb{H}_n$ that approximates
$\bm{\mathsf{A}}$ nearly as well as the best rank-$k$ matrix:
$$\Vert  \bm{\mathsf{A}} - \widehat{\bm{\mathsf{A}}}_k  \Vert \lesssim \sigma_{k+1}.$$
To obtain the approximation, we will adapt the randomized rangefinder
method
(Algorithm [\[alg:random-rangefinder\]](#alg:random-rangefinder){reference-type="ref"
reference="alg:random-rangefinder"}).

## The Nyström approximation

The most natural way to construct a low-rank approximation of a PSD
matrix is via the Nyström method. Let
$\bm{\mathsf{X}} \in \mathbb{F}^{n \times \ell}$ be an arbitrary test
matrix. The *Nyström approximation* of $\bm{\mathsf{A}}$ with respect to
$\bm{\mathsf{X}}$ is the PSD matrix $$\begin{equation}
 \label{eqn:nystrom-def}
\bm{\mathsf{A}}\langle \bm{\mathsf{X}} \rangle
    := (\bm{\mathsf{AX}})(\bm{\mathsf{X}}^* \bm{\mathsf{AX}})^\dagger(\bm{\mathsf{AX}})^*.
\end{equation}$$ An alternative presentation of this formula is
$$\bm{\mathsf{A}}\langle \bm{\mathsf{X}} \rangle
    = \bm{\mathsf{A}}^{1/2} \bm{\mathsf{P}}_{\bm{\mathsf{A}}^{1/2} \bm{\mathsf{X}}} \bm{\mathsf{A}}^{1/2},$$
where $\bm{\mathsf{P}}_{\bm{\mathsf{Y}}}$ is the orthogonal projector
onto the range of $\bm{\mathsf{Y}}$. In particular, the Nyström
approximation only depends on the range of the matrix $\bm{\mathsf{X}}$.

The Nyström
approximation [\[eqn:nystrom-def\]](#eqn:nystrom-def){reference-type="eqref"
reference="eqn:nystrom-def"} is closely related to the Schur
complement [\[eqn:schur-complement\]](#eqn:schur-complement){reference-type="eqref"
reference="eqn:schur-complement"} of $\bm{\mathsf{A}}$ with respect to
$\bm{\mathsf{X}}$. Indeed, $$\bm{\mathsf{A}} / \bm{\mathsf{X}}
    = \bm{\mathsf{A}} - \bm{\mathsf{A}} \langle \bm{\mathsf{X}} \rangle
    = \bm{\mathsf{A}}^{1/2} (\bm{\mathsf{I}}- \bm{\mathsf{P}}_{\bm{\mathsf{A}}^{1/2} \bm{\mathsf{X}}}) \bm{\mathsf{A}}^{1/2}.$$
That is, the Schur complement of $\bm{\mathsf{A}}$ with respect to
$\bm{\mathsf{X}}$ is precisely the error in the Nyström approximation.

Proposition [31](#prop:rrf-schur){reference-type="ref"
reference="prop:rrf-schur"} indicates that the Nyström decomposition is
also connected with our approach to solving the rangefinder problem. We
immediately perceive the opportunity to use a random test matrix
$\bm{\mathsf{X}}$ to form the Nyström approximation. Let us describe how
this choice leads to algorithms for computing a near-optimal low-rank
approximation of the matrix $\bm{\mathsf{A}}$.

## Randomized Nyström approximation algorithms {#sec:nystrom-alg}

Here is a simple and effective procedure for computing a rank-$k$ PSD
approximation of the PSD matrix $\bm{\mathsf{A}} \in \mathbb{H}_n$.

First, draw a random test matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times \ell}$, where
$\ell \geq k$. Form the sample matrix
$\bm{\mathsf{Y}} = \bm{\mathsf{A\Omega}} \in \mathbb{F}^{n \times \ell}$.
Then compute the Nyström approximation $$\begin{equation}
 \label{eqn:Ahat-nys-bad}
\widehat{\bm{\mathsf{A}}} = \bm{\mathsf{A}}\langle  \bm{\mathsf{\Omega}}  \rangle = \bm{\mathsf{Y}} (\bm{\mathsf{\Omega}}^* \bm{\mathsf{Y}})^\dagger\bm{\mathsf{Y}}^*.
\end{equation}$$ The initial approximation $\widehat{\bm{\mathsf{A}}}$
has rank $\ell$. To truncate the rank to $k$, we just report a best
rank-$k$ approximation $\widehat{\bm{\mathsf{A}}}_k$ of the initial
approximation $\widehat{\bm{\mathsf{A}}}$ with respect to the Frobenius
norm; see [@TYUC17:Fixed-Rank-Approximation; @PB19:Improved-Fixed-Rank]
and [@WGM19:Scalable-k-Means].

Let us warn the reader that the
formula [\[eqn:Ahat-nys-bad\]](#eqn:Ahat-nys-bad){reference-type="eqref"
reference="eqn:Ahat-nys-bad"} is not suitable for numerical computation.
See
Algorithm [\[alg:random-nystrom\]](#alg:random-nystrom){reference-type="ref"
reference="alg:random-nystrom"} for a numerically stable implementation
::: {#alg:random-nystrom .algorithm}
**Algorithm. Randomized Nyström approximation.**

```
**Input:** Psd target matrix $\mtx{A} \in \Sym_n(\F)$, rank $k$ for approximation,
number $\ell$ of samples
**Output:** Rank-$k$ PSD approximation $\widehat{\mtx{A}}_k \in \Sym_n(\F)$ expressed in factored form
$\widehat{\mtx{A}}_k = \mtx{U \Lambda U}^*$ where $\mtx{U} \in \F^{n \times k}$ is orthonormal
and $\mtx{\Lambda} \in \Sym_k(\F)$ is nonnegative and diagonal.

function RandomNyström($\mtx{A}$, $k$, $\ell$):
    Draw random matrix $\mtx{\Omega} \in \F^{n \times \ell}$
    Form $\mtx{Y} = \mtx{A\Omega}$
    $\nu = \sqrt{n} \, eps(norm(\mtx{Y}))$
    # Compute shift
    $\mtx{Y}_{\nu} = \mtx{Y} + \nu \mtx{\Omega}$
    # Samples of shifted matrix
    $\mtx{C} = chol(\mtx{\Omega}^* \mtx{Y}_{\nu})$
    $\mtx{B} = \mtx{Y}_{\nu} \mtx{C}^{-1}$
    # Triangular solve!
    $[\mtx{U}, \mtx{\Sigma}, \sim] = svd(\mtx{B})$
    # Dense SVD
    $\mtx{\Lambda} = \max\{0, \mtx{\Sigma}^2 - \nu \Id\}$
    # Remove shift
    $\mtx{U} = \mtx{U}(:, 1:k)$ and $\mtx{\Lambda} = \mtx{\Lambda}(1:k, 1:k)$
    # Truncate rank to $k$
```
:::


adapted from [@LLS+17:Algorithm-971] and
[@TYUC17:Fixed-Rank-Approximation]. In general, the matrix--matrix
multiply with $\bm{\mathsf{A}}$ dominates the cost, with $O(n^2 \ell)$
arithmetic operations; this expense can be reduced if either
$\bm{\mathsf{A}}$ or $\bm{\mathsf{\Omega}}$ admits fast multiplication.
The approximation steps involve $O(n \ell^2)$ arithmetic. Meanwhile,
storage costs are $O(n\ell)$.

Another interesting aspect of
Algorithm [\[alg:random-nystrom\]](#alg:random-nystrom){reference-type="ref"
reference="alg:random-nystrom"} is that it only uses linear information
about the matrix $\bm{\mathsf{A}}$. Therefore, it can be implemented in
the one-pass or the streaming data model.
Remark [47](#rem:streaming){reference-type="ref"
reference="rem:streaming"} gives more details about matrix approximation
in the streaming model. See
Section [19.3.5](#sec:streaming-kpca){reference-type="ref"
reference="sec:streaming-kpca"} for an application to kernel principal
component analysis.

:::: algorithm
::: algorithmic
Psd target matrix $\bm{\mathsf{A}} \in \mathbb{H}_n(\mathbb{F})$, rank
$k$ for approximation, number $\ell$ of samples Rank-$k$ PSD
approximation $\widehat{\bm{\mathsf{A}}}_k \in \mathbb{H}_n(\mathbb{F})$
expressed in factored form
$\widehat{\bm{\mathsf{A}}}_k = \bm{\mathsf{U \Lambda U}}^*$ where
$\bm{\mathsf{U}} \in \mathbb{F}^{n \times k}$ is orthonormal and
$\bm{\mathsf{\Lambda}} \in \mathbb{H}_k(\mathbb{F})$ is nonnegative and
diagonal.

Draw random matrix $\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times \ell}$
Form $\bm{\mathsf{Y}} = \bm{\mathsf{A\Omega}}$

$\nu = \sqrt{n} \, \texttt{eps}(\texttt{norm}(\bm{\mathsf{Y}}))$ Compute
shift

$\bm{\mathsf{Y}}_{\nu} = \bm{\mathsf{Y}} + \nu \bm{\mathsf{\Omega}}$
Samples of shifted matrix

$\bm{\mathsf{C}} = \texttt{chol}(\bm{\mathsf{\Omega}}^* \bm{\mathsf{Y}}_{\nu})$

$\bm{\mathsf{B}} = \bm{\mathsf{Y}}_{\nu} \bm{\mathsf{C}}^{-1}$
Triangular solve!

$[\bm{\mathsf{U}}, \bm{\mathsf{\Sigma}}, \sim] = \texttt{svd}(\bm{\mathsf{B}})$
Dense SVD

$\bm{\mathsf{\Lambda}} = \max\{0, \bm{\mathsf{\Sigma}}^2 - \nu \bm{\mathsf{I}}\}$
Remove shift

$\bm{\mathsf{U}} = \bm{\mathsf{U}}(:, 1:k)$ and
$\bm{\mathsf{\Lambda}} = \bm{\mathsf{\Lambda}}(1:k, 1:k)$ Truncate rank
to $k$
:::
::::

## Analysis

The randomized Nyström method enjoys the same kind of guarantees as the
randomized rangefinder. The following
result [@TYUC17:Fixed-Rank-Approximation Thm. 4.1] extends earlier
contributions from [@HMT11:Finding-Structure] and
[@Git13:Topics-Randomized].

::: {#thm:rand-nys .theorem}
**Theorem 46** (Nyström: Gaussian analysis). *Fix a PSD matrix
$\bm{\mathsf{A}} \in \mathbb{H}_n(\mathbb{F})$ with eigenvalues
$\lambda_1 \geq \lambda_2 \geq \dots$. Let $1 \leq k < \ell \leq n$.
Draw a random test matrix
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times \ell}$ that is standard
normal. Then the rank-$k$ PSD approximation
$\widehat{\bm{\mathsf{A}}}_k$ computed by
Algorithm [\[alg:random-nystrom\]](#alg:random-nystrom){reference-type="ref"
reference="alg:random-nystrom"} satisfies
$$\operatorname{\mathbb{E}}\Vert  \bm{\mathsf{A}} - \widehat{\bm{\mathsf{A}}}_k  \Vert
    \leq \lambda_{k+1} + \frac{k}{\ell - k - 1} \left( \sum\nolimits_{j > k} \lambda_j \right).$$*
:::

In other words, the computed rank-$k$ approximation
$\widehat{\bm{\mathsf{A}}}_k$ achieves almost the error $\lambda_{k+1}$
in the optimal rank-$k$ approximation of $\bm{\mathsf{A}}$. The error
declines as the number $\ell$ of samples increases and as the $\ell_1$
norm of the tail eigenvalues decreases.

When the target matrix $\bm{\mathsf{A}}$ has a sharply decaying
spectrum, Theorem [46](#thm:rand-nys){reference-type="ref"
reference="thm:rand-nys"} can be pessimistic. See  for additional
theoretical results.

## Powering {#sec:nystrom-power}

We can reduce the error in the randomized Nyström approximation by
powering the input matrix, much as subspace iteration improves the
performance in the randomized rangefinder
(Section [11.6](#sec:rrf-subspace){reference-type="ref"
reference="sec:rrf-subspace"}).

Let $\bm{\mathsf{A}} \in \mathbb{H}_n$ be a PSD matrix. Draw a random
test matrix $\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times \ell}$. For a
natural number $q$, we compute
$\bm{\mathsf{Y}} = \bm{\mathsf{A}}^q \bm{\mathsf{\Omega}}$ by repeated
multiplication. Then the Nyström approximation of the input matrix takes
the form
$$\widehat{\bm{\mathsf{A}}} = \big[ (\bm{\mathsf{A}}^q) \langle  \bm{\mathsf{\Omega}}  \rangle \big]^{1/q}
    = \big[ \bm{\mathsf{Y}} (\bm{\mathsf{\Omega}}^* \bm{\mathsf{Y}})^\dagger\bm{\mathsf{Y}}^* \big]^{1/q}.$$
This approach requires very careful numerical implementation. As in the
case of subspace iteration, it drives down the error exponentially fast
as $q$ increases. We omit the details.

## History

The Nyström approximation was developed in the context of integral
equations [@Nys30:Uber-Praktische]. It has had a substantial impact in
machine learning, beginning with the work of   on randomized low-rank
approximation of kernel matrices.
Section [19](#sec:kernel){reference-type="ref" reference="sec:kernel"}
contains a discussion of this literature. Note that the Nyström
approximation of a kernel matrix is almost always computed with respect
to a random coordinate subspace, in contrast to the uniformly random
subspace induced by a Gaussian test matrix.

Algorithmic and theoretical results on Nyström approximation with
respect to general test matrices have appeared in a number of papers,
including
[@HMT11:Finding-Structure; @Git13:Topics-Randomized; @LLS+17:Algorithm-971]
and [@2017_tropp_practical_sketching].

# Single-view algorithms {#sec:singlepass}

In this section, we will describe a remarkable class of algorithms that
are capable of computing a low-rank approximation of a matrix that is so
large that it cannot be stored at all.

We will consider the specific problem of computing an approximate
singular value decomposition of a matrix
$\bm{\mathsf{A}} \in \mathbb{F}^{m\times n}$ under the assumption that
we are allowed to view each entry of $\bm{\mathsf{A}}$ only once and
that we cannot specify the order in which they are viewed. To the best
of our knowledge, no deterministic techniques can carry off such a
computation without *a priori* information about the singular vectors of
the matrix.

For the case where $\bm{\mathsf{A}}$ is psd, we have already seen a
single-view algorithm: the Nyström technique of
Algorithm [\[alg:random-nystrom\]](#alg:random-nystrom){reference-type="ref"
reference="alg:random-nystrom"}. Here, we concentrate on the more
difficult case of general matrices. This presentation is adapted from
and the papers [@TYUC17:Practical-Sketching] and
[@TYUC19:Streaming-Low-Rank].

## Algorithms {#sec:singleviewalgorithms}

In the basic RSVD algorithm (Section
[11.2](#sec:rsvd){reference-type="ref" reference="sec:rsvd"}), we view
each element of the given matrix $\bm{\mathsf{A}}$ at least twice. In
the first view, we form a sample matrix
$\bm{\mathsf{Y}} = \bm{\mathsf{A}}\bm{\mathsf{\Omega}}$ for a given test
matrix $\bm{\mathsf{\Omega}}$. We orthonormalize the columns of
$\bm{\mathsf{Y}}$ to form the matrix $\bm{\mathsf{Q}}$ and then visit
$\bm{\mathsf{A}}$ again to form a second sample
$\bm{\mathsf{C}} = \bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}$. The columns of
$\bm{\mathsf{Y}}$ form an approximate basis for the column space of
$\bm{\mathsf{A}}$, and the columns of $\bm{\mathsf{C}}$ form an
approximate basis for the row space.

In the single-view framework, we can only visit $\bm{\mathsf{A}}$ once,
which means that we must sample both the row and the column space
simultaneously. To this end, let us draw tall thin random matrices
$$\begin{equation}
\label{eq:single1}
\bm{\mathsf{\Upsilon}} \in \mathbb{F}^{m\times \ell}
\qquad\mbox{and}\qquad
\bm{\mathsf{\Omega}} \in \mathbb{F}^{n\times \ell}
\end{equation}$$ and then form the two corresponding sample matrices
$$\begin{equation}
\label{eq:single2}
\bm{\mathsf{X}} = \bm{\mathsf{A}}^{*}\bm{\mathsf{\Upsilon}} \in \mathbb{F}^{n\times \ell}
\qquad\mbox{and}\qquad
\bm{\mathsf{Y}} = \bm{\mathsf{A}}\bm{\mathsf{\Omega}} \in \mathbb{F}^{m\times \ell}.
\end{equation}$$ In ([\[eq:single1\]](#eq:single1){reference-type="ref"
reference="eq:single1"}), we draw a number $\ell$ of samples that is
slightly larger than the rank $k$ of the low-rank approximation that we
seek. (Section [15.2](#sec:singleapriori){reference-type="ref"
reference="sec:singleapriori"} gives details about how to choose
$\ell$.) Observe that both $\bm{\mathsf{X}}$ and $\bm{\mathsf{Y}}$ can
be formed in a single pass over the matrix $\bm{\mathsf{A}}$.

Once we have seen the entire matrix, the next step is to orthonormalize
the columns of $\bm{\mathsf{X}}$ and $\bm{\mathsf{Y}}$ to obtain
orthonormal matrices $$\begin{equation}
\label{eq:single3}
[\bm{\mathsf{P}},\sim] = \texttt{qr\_econ}(\bm{\mathsf{X}}),
\qquad\mbox{and}\qquad
[\bm{\mathsf{Q}},\sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}}).
\end{equation}$$ At this point, $\bm{\mathsf{P}}$ and $\bm{\mathsf{Q}}$
hold approximate bases for the row and column spaces of
$\bm{\mathsf{A}}$, so we anticipate that $$\begin{equation}
\label{eq:single4}
\bm{\mathsf{A}} \approx \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{P}}\bm{\mathsf{P}}^{*} = \bm{\mathsf{Q}}\bm{\mathsf{C}}\bm{\mathsf{P}}^{*},
\end{equation}$$ where we defined the "core" matrix $$\begin{equation}
\label{eq:single5}
\bm{\mathsf{C}} := \bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{P}} \in \mathbb{F}^{\ell\times \ell}.
\end{equation}$$ Unfortunately, since we cannot revisit
$\bm{\mathsf{A}}$, we are not allowed to form $\bm{\mathsf{C}}$ directly
by applying formula [\[eq:single5\]](#eq:single5){reference-type="eqref"
reference="eq:single5"}.

Instead, we develop a relation that $\bm{\mathsf{C}}$ must satisfy
approximately, which allows us to estimate $\bm{\mathsf{C}}$ from the
quantities we have on hand. To do so, right-multiply the definition
[\[eq:single5\]](#eq:single5){reference-type="eqref"
reference="eq:single5"} by $\bm{\mathsf{P}}^{*}\bm{\mathsf{\Omega}}$ to
obtain $$\begin{equation}
\label{eq:single6}
\bm{\mathsf{C}}\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Omega}}\bigr) = \bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{P}}\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Omega}}\bigr).
\end{equation}$$ Inserting the approximation
$\bm{\mathsf{A}}\bm{\mathsf{P}}\bm{\mathsf{P}}^{*} \approx \bm{\mathsf{A}}$
into ([\[eq:single6\]](#eq:single6){reference-type="ref"
reference="eq:single6"}), we find that $$\begin{equation}
\label{eq:single7}
\bm{\mathsf{C}}\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Omega}}\bigr) \approx \bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{\Omega}} = \bm{\mathsf{Q}}^{*}\bm{\mathsf{Y}}.
\end{equation}$$ In ([\[eq:single7\]](#eq:single7){reference-type="ref"
reference="eq:single7"}), all quantities except $\bm{\mathsf{C}}$ are
known explicitly, which means that we can solve it, in the least-squares
sense, to arrive at an estimate $$\begin{equation}
\label{eq:single8}
\bm{\mathsf{C}}_{\mathrm{approx}} = \bigl(\bm{\mathsf{Q}}^{*}\bm{\mathsf{Y}}\bigr)\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Omega}}\bigr)^{\dagger},
\end{equation}$$ where $\dagger$ denotes the Moore-Penrose
pseudoinverse. As always, the pseudoinverse is applied by means of an
orthogonal factorization. Once $\bm{\mathsf{C}}_{\mathrm{approx}}$ has
been computed via ([\[eq:single8\]](#eq:single8){reference-type="ref"
reference="eq:single8"}), we obtain the rank-$\ell$ approximation
$$\begin{equation}
\label{eq:single8b}
\bm{\mathsf{A}} \approx \bm{\mathsf{Q}}\bm{\mathsf{C}}_{\mathrm{approx}}\bm{\mathsf{P}}^{*},
\end{equation}$$ which we can convert into an approximate SVD using the
standard postprocessing steps. For additional implementation details,
see and . Extensions of this approach, with theoretical analysis, appear
in .

Recently, have demonstrated that the numerical performance of the
single-view algorithm can be improved by extracting a third sketch of
$\bm{\mathsf{A}}$ that is independent from $\bm{\mathsf{X}}$ and
$\bm{\mathsf{Y}}$. The idea is to draw tall thin random matrices
$$\bm{\mathsf{\Phi}} \in \mathbb{F}^{m\times s},
\qquad\mbox{and}\qquad
\bm{\mathsf{\Psi}} \in \mathbb{F}^{n\times s},$$ where $s$ is another
oversampling parameter. Then we form a "core sketch" $$\begin{equation}
\label{eq:single9}
\bm{\mathsf{Z}} = \bm{\mathsf{\Phi}}^{*}\bm{\mathsf{A}}\bm{\mathsf{\Psi}} \in \mathbb{F}^{s\times s}.
\end{equation}$$ This extra data allows us to derive an alternative
equation for the core matrix $\bm{\mathsf{C}}$. We left- and
right-multiply the definition
[\[eq:single5\]](#eq:single5){reference-type="eqref"
reference="eq:single5"} by $\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}}$ and
$\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}}$ to obtain the relation
$$\begin{equation}
\label{eq:single10}
\bigl(\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}}\bigr)\bm{\mathsf{C}}\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}}\bigr) =
\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{P}}\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}}.
\end{equation}$$ Inserting the approximation
$\bm{\mathsf{A}} \approx \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{P}}\bm{\mathsf{P}}^{*}$
into [\[eq:single10\]](#eq:single10){reference-type="eqref"
reference="eq:single10"}, we find that $$\begin{equation}
\label{eq:single11}
\bigl(\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}}\bigr)\bm{\mathsf{C}}\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}}\bigr) \approx
\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{A}}\bm{\mathsf{\Psi}} = \bm{\mathsf{Z}}.
\end{equation}$$ An improved approximation to the core matrix
$\bm{\mathsf{C}}$ results by solving
[\[eq:single11\]](#eq:single11){reference-type="eqref"
reference="eq:single11"} in a least-squares sense; to wit,
$\bm{\mathsf{C}} = (\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}})^{\dagger}\bm{\mathsf{Z}}(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}})^{\dagger}$.

::: {#rem:streaming .remark}
**Remark 47** (Streaming algorithms). *Single-view algorithms are
related to the streaming model of computation
[@2005_muthukrishnan_stream]. were the first to explicitly study matrix
computations in streaming data models.*

*One important streaming model poses the assumption that the input
matrix $\bm{\mathsf{A}}$ is presented as a sequence of innovations:
$$\bm{\mathsf{A}} = \bm{\mathsf{H}}_{1} + \bm{\mathsf{H}}_{2} + \bm{\mathsf{H}}_{3} + \cdots$$
Typically, each update $\bm{\mathsf{H}}_{i}$ is simple; for instance, it
may be sparse or low-rank. The challenge is that the full matrix
$\bm{\mathsf{A}}$ is too large to be stored. Once an innovation
$\bm{\mathsf{H}}_{i}$ has been processed, it cannot be retained. This is
called the "turnstile" model in the theoretical computer science
literature.*

*The algorithms described in this section handle this difficulty by
creating a random linear transform $\mathcal{S}$ that maps
$\bm{\mathsf{A}}$ down to a low-dimensional sketch that is small enough
to store. What we actually retain in memory is the evolving sketch of
the input: $$\mathcal{S}(\bm{\mathsf{A}}) =
\mathcal{S}(\bm{\mathsf{H}}_{1}) +
\mathcal{S}(\bm{\mathsf{H}}_{2}) +
\mathcal{S}(\bm{\mathsf{H}}_{3}) + \cdots.$$ In Algorithm
[\[alg:SingleViewSVD\]](#alg:SingleViewSVD){reference-type="ref"
reference="alg:SingleViewSVD"}, we instantiate $\mathcal{S}$ by drawing
::: {#alg:SingleViewSVD .algorithm}
**Algorithm. Single-view SVD.**

```
**Input:** Target matrix $\mtx{A} \in \F^{m \times n}$, rank $k$, sampling sizes $\ell$ and $s$.
**Output:** Orthonormal matrices $\mtx{U} \in \F^{m\times k}$ and $\mtx{V} = \F^{n\times k}$, and a diagonal matrix $\mtx{\Sigma} \in \F^{k\times k}$ such that
$\mtx{A} \approx \mtx{U}\mtx{\Sigma}\mtx{V}^{*}$.

function SingleViewSVD($\mtx{A}$, $k$, $\ell$, $s$):
    Draw test matrices $\mtx{\Upsilon} \in \F^{m\times \ell}$, $\mtx{\Omega} \in \F^{n\times \ell}$,
    $\mtx{\Phi} \in \F^{m\times s}$, $\mtx{\Psi} \in \F^{n\times s}$
    Form $\mtx{X} = \mtx{A}^{*}\mtx{\Upsilon}$, $\mtx{Y} = \mtx{A}\mtx{\Omega}$, $\mtx{Z} = \mtx{\Phi}^{*}\mtx{A}\mtx{\Psi}$
    # Viewing $\mtx{A}$ only once!
    $[\mtx{P},\sim] = qr_econ(\mtx{X})$, $[\mtx{Q},\sim] = qr_econ(\mtx{Y})$
    $\mtx{C} = \bigl(\mtx{\Phi}^{*}\mtx{Q}\bigr)^{\pinv}\,\mtx{Z}\,\bigl(\mtx{P}^{*}\mtx{\Psi}\bigr)^{\pinv}$
    # Execute using a least-squares solver
    $[\widehat{\mtx{U}},\widehat{\mtx{\Sigma}},\widehat{\mtx{V}}] = svd(\mtx{C})$
    # A full SVD
    $\mtx{U} = \mtx{Q}\widehat{\mtx{U}}(\colon,1:k)$,
    $\mtx{V} = \mtx{P}\widehat{\mtx{V}}(\colon,1:k)$,
    $\mtx{\Sigma} = \widehat{\mtx{\Sigma}}(1:k,1:k)$
```
:::


the random matrices $\bm{\mathsf{\Upsilon}}$, $\bm{\mathsf{\Omega}}$,
$\bm{\mathsf{\Phi}}$, and $\bm{\mathsf{\Psi}}$, and then work with the
sketch $$\mathcal{S}(\bm{\mathsf{H}}) =
\big(\bm{\mathsf{\Upsilon}}^{*}\bm{\mathsf{H}},\,\bm{\mathsf{H}}\bm{\mathsf{\Omega}},\,\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{H}}\bm{\mathsf{\Psi}} \big).$$
The fact that the sketch is a *linear* map is essential here. prove that
randomized linear embeddings are essentially the only kind of algorithm
for handling the turnstile model. In contrast, the sketch implicit in
the RSVD algorithm from Section [11.2](#sec:rsvd){reference-type="ref"
reference="sec:rsvd"} is a quadratic or higher-order polynomial in the
input matrix.*
:::

:::: algorithm
::: algorithmic
Target matrix $\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$, rank $k$,
sampling sizes $\ell$ and $s$. Orthonormal matrices
$\bm{\mathsf{U}} \in \mathbb{F}^{m\times k}$ and
$\bm{\mathsf{V}} = \mathbb{F}^{n\times k}$, and a diagonal matrix
$\bm{\mathsf{\Sigma}} \in \mathbb{F}^{k\times k}$ such that
$\bm{\mathsf{A}} \approx \bm{\mathsf{U}}\bm{\mathsf{\Sigma}}\bm{\mathsf{V}}^{*}$.

Draw test matrices
$\bm{\mathsf{\Upsilon}} \in \mathbb{F}^{m\times \ell}$,
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n\times \ell}$,
$\bm{\mathsf{\Phi}} \in \mathbb{F}^{m\times s}$,
$\bm{\mathsf{\Psi}} \in \mathbb{F}^{n\times s}$ Form
$\bm{\mathsf{X}} = \bm{\mathsf{A}}^{*}\bm{\mathsf{\Upsilon}}$,
$\bm{\mathsf{Y}} = \bm{\mathsf{A}}\bm{\mathsf{\Omega}}$,
$\bm{\mathsf{Z}} = \bm{\mathsf{\Phi}}^{*}\bm{\mathsf{A}}\bm{\mathsf{\Psi}}$
Viewing $\bm{\mathsf{A}}$ only once!
$[\bm{\mathsf{P}},\sim] = \texttt{qr\_econ}(\bm{\mathsf{X}})$,
$[\bm{\mathsf{Q}},\sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}})$
$\bm{\mathsf{C}} = \bigl(\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}}\bigr)^{\dagger}\,\bm{\mathsf{Z}}\,\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}}\bigr)^{\dagger}$
Execute using a least-squares solver
$[\widehat{\bm{\mathsf{U}}},\widehat{\bm{\mathsf{\Sigma}}},\widehat{\bm{\mathsf{V}}}] = \texttt{svd}(\bm{\mathsf{C}})$
A full SVD
$\bm{\mathsf{U}} = \bm{\mathsf{Q}}\widehat{\bm{\mathsf{U}}}(\colon,1:k)$,
$\bm{\mathsf{V}} = \bm{\mathsf{P}}\widehat{\bm{\mathsf{V}}}(\colon,1:k)$,
$\bm{\mathsf{\Sigma}} = \widehat{\bm{\mathsf{\Sigma}}}(1:k,1:k)$
:::
::::

::: remark
**Remark 48** (Single-view versus out-of-core algorithms). *In
principle, the methods discussed in this section can also be used in
situations where a matrix is stored in slow memory, such as a spinning
disk hard drive, or on a distributed memory system. However, one has to
carefully weigh whether the decrease in accuracy and increase in
uncertainty that is inherent to single-view algorithms is worth the cost
savings. As a general matter, revisiting the matrix at least once is
advisable whenever it is possible.*
:::

## Error estimation, parameter choices and truncation {#sec:singleapriori}

In the single-view computing environment, one must choose sampling
parameters before the computation starts, and there is no way to revisit
these choices after data has been gathered. This constraint makes *a
priori* error analysis particularly important, because we need guidance
on how large to make the sketches given some prior knowledge about the
spectral decay of the input matrix. To illustrate how this may work, let
us cite :

::: {#thm:singleview .theorem}
**Theorem 49** (Single-view SVD: Gaussian analysis). *Suppose that
Algorithm
[\[alg:SingleViewSVD\]](#alg:SingleViewSVD){reference-type="ref"
reference="alg:SingleViewSVD"} is executed for an input matrix
$\bm{\mathsf{A}} \in \mathbb{C}^{m\times n}$ and for sampling parameters
$s$ and $\ell$ that satisfy $s \geq 2\ell$. When the test matrices are
drawn from a standard normal distribution, the computed matrices
$\bm{\mathsf{P}}$, $\bm{\mathsf{C}}$, and $\bm{\mathsf{Q}}$ satisfy
$$\operatorname{\mathbb{E}}\Vert  \bm{\mathsf{A}} - \bm{\mathsf{Q}}\bm{\mathsf{C}}\bm{\mathsf{P}}^{*}  \Vert_{\mathrm{F}}^{2}
\leq
\frac{s}{s-\ell}\,
\min_{k < \ell}\left(\frac{\ell+k}{\ell-k}\sum_{j=k+1}^{\min(m,n)}\sigma_{j}^{2}\right).$$
As usual, $\sigma_j$ is the $j$th largest singular value of
$\bm{\mathsf{A}}$. A very similar bound holds for the real field.*
:::

This result suggests that more aggressive oversampling is called for in
the single-view setting, as compared to the basic rangefinder problem.
For instance, if we aim for an approximation error that is comparable to
the best possible approximation with rank $k$, then we might choose
$\ell = 4k$ and $s=8k$ to obtain
$$\mathbb{E}\|\bm{\mathsf{A}} - \bm{\mathsf{Q}}\bm{\mathsf{C}}\bm{\mathsf{P}}^{*}\|_{\mathrm{F}}^{2}
\leq
\frac{10}{3}\sum_{j=k+1}^{\min(m,n)}\sigma_{j}^{2} =
\frac{10}{3}\|\bm{\mathsf{A}} - \bm{\mathsf{A}}_{k}\|_{\mathrm{F}}^{2},$$
where $\bm{\mathsf{A}}_{k}$ is the best possible rank-$k$ approximation
of $\bm{\mathsf{A}}$. As usual, the likelihood of large deviations from
the expectation is negligible. (For contrast, recall that the basic
rangefinder algorithm often works well when we select $\ell = k+5$ or
$\ell=k+10$.)

Besides computing $\bm{\mathsf{P}}$, $\bm{\mathsf{C}}$, and
$\bm{\mathsf{Q}}$, Algorithm
[\[alg:SingleViewSVD\]](#alg:SingleViewSVD){reference-type="ref"
reference="alg:SingleViewSVD"} also prunes the approximation
$\bm{\mathsf{A}} \approx \bm{\mathsf{Q}}\bm{\mathsf{C}}\bm{\mathsf{P}}^{*}$
by computing an SVD of $\bm{\mathsf{C}}$ (line 6) and then throwing out
the trailing $\ell-k$ modes (line 7). The motivation for this truncation
is that the approximation
$\bm{\mathsf{A}} \approx \bm{\mathsf{Q}}\bm{\mathsf{C}}\bm{\mathsf{P}}^{*}$,
tends to capture the dominant singular modes of $\bm{\mathsf{A}}$ well,
but the trailing ones have very low accuracy. The same thing happens
with the basic RSVD (Section [11.2](#sec:rsvd){reference-type="ref"
reference="sec:rsvd"}), but the phenomenon is more pronounced in the
single-view environment, in part because $\ell$ is substantially larger
than $k$. Theorem [49](#thm:singleview){reference-type="ref"
reference="thm:singleview"} can be applied to prove that the truncated
factorization is as accurate as one can reasonably hope for; see  for
details.

::: remark
**Remark 50** (Spectral norm bounds?).
*Theorem [49](#thm:singleview){reference-type="ref"
reference="thm:singleview"} provides a Frobenius norm error bound for a
matrix approximation algorithm. For our survey, this is a *rara avis in
terra*. Unfortunately, relative error spectral norm error bounds are not
generally possible in the streaming setting [@2014_woodruff_sketching
Chap. 6].*
:::

## Structured test matrices

Algorithm [\[alg:SingleViewSVD\]](#alg:SingleViewSVD){reference-type="ref"
reference="alg:SingleViewSVD"} can -- and should -- be implemented with
structured test matrices, rather than Gaussian test matrices. This
modification is especially appealing in the single-view environment,
where storage is often the main bottleneck.

For instance, consider the parameter selections $\ell = 4k$ and $s = 8k$
that we referenced above. Then the four test matrices consist of
$12k(m+n)$ floats that must be stored, and the sketches add another
$4k(m+n) + 64k^{2}$ floats. Since $m$ and $n$ can be huge, these numbers
could severely limit the rank $k$ of the final matrix approximation.

If we swap out the Gaussian matrices for structured random matrices, we
can almost remove the cost associated with storing the test matrices. In
particular, the addition of the core sketch
[\[eq:single9\]](#eq:single9){reference-type="eqref"
reference="eq:single9"} has a very light memory footprint because the
sketch itself only uses $O(k^{2})$ floats. Empirically, when we use a
structured random matrix, such as a sparse sign matrix
(Section [9.2](#sec:sparse-map){reference-type="ref"
reference="sec:sparse-map"}) or an SRTT
(Section [9.3](#sec:srtt){reference-type="ref" reference="sec:srtt"}),
the observed errors are more or less indistinguishable from the errors
attained with Gaussian test matrices. See .

## A posteriori error estimation {#sec:singleaposteriori}

In order to reduce the uncertainty associated with the single-view
algorithms described in this section, the "certificate of accuracy"
technique described in Section
[12.2](#sec:certificate){reference-type="ref"
reference="sec:certificate"} is very useful.

Recall that the idea is to draw a separate test matrix whose only
purpose is to provide an independent estimate of the error in the
computed solution. This additional test matrix can be *very* thin (say 5
or 10 columns wide), and it still provides a dependable bound on the
computed error. These techniques can be incorporated without any
difficulty in the single-view environment, as outlined in .

Let us mention one caveat. In the single-view environment, we have no
recourse when the *a posteriori* error estimator signals that the
approximation error is unacceptable. On the other hand, it is reassuring
that the algorithm can sound a warning that it has not met the desired
accuracy.

## History {#sec:singlehistory}

To the best of our knowledge, described the first algorithm that can
compute a low-rank matrix approximation in the single-view computational
model. Their paper introduced the idea of independently sampling the
row- and the column-space of a matrix, as summarized in formulas
[\[eq:single1\]](#eq:single1){reference-type="eqref"
reference="eq:single1"}--[\[eq:single8b\]](#eq:single8b){reference-type="eqref"
reference="eq:single8b"}. This approach inspired the single-view
algorithms presented in . It is interesting that the primary objective
of was to reduce the asymptotic flop count of the computation through
the use of structured random test matrices.

gave an explicit discussion of randomized NLA in a streaming
computational model. They independently proposed a variant of the
algorithm from [@HMT11:Finding-Structure Sec. 5.5]. Later contributions
to the field appeared in
[@LNW14:Turnstile-Streaming; @2016_boutsidis_optimal; @2016_feldman_dimensionality; @GLPW16:Frequent-Directions]
and [@2017_tropp_practical_sketching]. The idea of introducing an
additional sketch such as
[\[eq:single9\]](#eq:single9){reference-type="eqref"
reference="eq:single9"} to capture the "core" matrix was proposed by .
have provided improvements of his approach, further analysis, and
computational considerations.

# Factoring matrices of full or nearly full rank {#sec:full}

So far, we have focused on techniques for computing low-rank
approximations of an input matrix. We will now upgrade to techniques for
computing *full* rank-revealing factorizations such as the
column-pivoted QR (CPQR) decomposition.

Classical deterministic techniques for computing these factorizations
proceed through a sequence of rank-one updates to the matrix, making
them communication-intensive and slow when executed on modern hardware.
Randomization allows the algorithms to be reorganized so that the vast
majority of the arithmetic takes place inside matrix--matrix
multiplications, which greatly accelerates the execution speed.

When applied to an $n \times n$ matrix, most of the algorithms described
in this section have the same $O(n^{3})$ asymptotic complexity as
traditional methods; the objective is to improve the practical speed by
reducing communication costs. However, randomization also allows us to
incorporate Strassen-type techniques to accelerate the matrix
multiplications in a numerically stable manner, attaining an overall
cost of $O(n^{\omega})$ where $\omega$ is the exponent of square
matrix--matrix multiplication.

As well as the CPQR decomposition, we will consider algorithms for
computing factorizations of the form
$\bm{\mathsf{A}} = \bm{\mathsf{U}}\bm{\mathsf{R}}\bm{\mathsf{V}}^{*}$
where $\bm{\mathsf{U}}$ and $\bm{\mathsf{V}}$ are unitary matrices and
$\bm{\mathsf{R}}$ is upper-triangular. Factorizations of this form can
be used for almost any task where either the CPQR or the SVD is
currently used. The additional flexibility allows us to improve on both
the computational speed and on the rank-revealing qualities of the
factorization.

Sections [16.1](#sec:rankrevealing){reference-type="ref"
reference="sec:rankrevealing"}--[16.4](#sec:demmelURV){reference-type="ref"
reference="sec:demmelURV"} introduce the key concepts by describing a
simple algorithm for computing a rank-revealing factorization of a
matrix. This method is both faster than traditional column-pivoted QR
and better at revealing the spectral properties of the matrix. Sections
[16.5](#sec:classicCPQR){reference-type="ref"
reference="sec:classicCPQR"}--[16.7](#sec:randUTV){reference-type="ref"
reference="sec:randUTV"} are more technical; they describe how
randomization can be used to resolve a longstanding challenge of how to
*block* a classical algorithm for computing a column-pivoted QR
decomposition by applying groups of Householder reflectors
simultaneously. They also describe how these ideas can be extended to
the task of computing a URV factorization.

## Rank-revealing factorizations {#sec:rankrevealing}

Before we discuss algorithms, let us first define what we mean when we
say that a factorization is *rank-revealing.* Given an $m\times n$
matrix $\bm{\mathsf{A}}$, we will consider factorizations of the form
$$\begin{equation}
\label{eq:fullfactintro}
\begin{array}{ccccccccccccccccccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{U}} &\bm{\mathsf{R}}&\bm{\mathsf{V}}^{*},\\
m\times n && m\times c & c\times n & n\times n
\end{array}
\end{equation}$$ where $c = \min(m,n)$, where $\bm{\mathsf{U}}$ and
$\bm{\mathsf{V}}$ are orthonormal, and where $\bm{\mathsf{R}}$ is
upper-triangular (or banded upper-triangular). We want the factorization
to reveal the numerical rank of $\bm{\mathsf{A}}$ in the sense that we
obtain a near-optimal approximation of $\bm{\mathsf{A}}$ when we
truncate [\[eq:fullfactintro\]](#eq:fullfactintro){reference-type="eqref"
reference="eq:fullfactintro"} to any level $k$. That is,
$$\begin{equation}
\label{eq:rankreveal}
\|\bm{\mathsf{A}} - \bm{\mathsf{U}}(\colon,1:k)\bm{\mathsf{R}}(1:k,:)\bm{\mathsf{V}}^{*}\|
\approx
\inf\{\|\bm{\mathsf{A}} - \bm{\mathsf{B}}\|\,\colon\,\bm{\mathsf{B}}\mbox{ has rank }k\}
\end{equation}$$ for $k \in \{1,2,\dots,c\}$. The factorization
[\[eq:fullfactintro\]](#eq:fullfactintro){reference-type="eqref"
reference="eq:fullfactintro"} can be viewed either as a generalization
of the SVD (for which $\bm{\mathsf{R}}$ is diagonal) or as a
generalization of the column-pivoted QR factorization (for which
$\bm{\mathsf{V}}$ is a permutation matrix).

A factorization such as
[\[eq:fullfactintro\]](#eq:fullfactintro){reference-type="eqref"
reference="eq:fullfactintro"} that satisfies
[\[eq:rankreveal\]](#eq:rankreveal){reference-type="eqref"
reference="eq:rankreveal"} is very handy. It can be used to solve
ill-conditioned linear systems or least-squares problems; it can be used
for estimating the singular spectrum of $\bm{\mathsf{A}}$ (and all
Schatten $p$-norms); and it provides orthonormal bases for
approximations to the four fundamental subspaces of the matrix. Finally,
it can be used to compute approximate low-rank factorizations
efficiently in situations where the numerical rank of the matrix is not
that much smaller than $m$ or $n$. In contrast, all techniques described
up to now are efficient only when the numerical rank $k$ satisfies
$k \ll \min(m,n)$.

## Blocking of matrix computations {#sec:blocking}

A well-known feature of modern computing is that we can execute
increasingly many floating-point operations because CPUs are gaining
more cores while GPUs and other accelerators are becoming more
affordable and more energy-efficient. In contrast, the cost of
communication (data transfer up and down levels of a memory hierarchy,
among servers, and across networks, etc.) is declining very slowly. As a
result, reducing communication is often the key to accelerating
numerical algorithms in the real world.

In the context of matrix computations, the main reaction to this
development has been to cast linear-algebraic operations as operating on
blocks of the matrix, rather than on individual entries or individual
columns and rows
[@1998_dongarra_blocking; @1993_stewart_blockQR_SVD; @2015_blockQR_SISC].
The objective is to reorganize an algorithm so that the majority of
flops can be executed using highly efficient algorithms for
matrix--matrix multiplication (BLAS3), rather than the slower methods
for matrix--vector multiplications (BLAS2).

Unfortunately, it turns out that classical algorithms for computing
rank-revealing factorizations of matrices are very challenging to block.
Column-pivoted QR proceeds through a sequence of rank-1 updates to the
matrix. The next pivot cannot be found until the previous update has
been applied. Techniques for computing an SVD of a matrix start by
reducing the matrix to bidiagonal form. Then they iterate on the
bidiagonal matrix to drive it towards diagonal form. Both steps are
challenging to block.

To emphasize just how much of a difference blocking makes, let us peek
ahead at the plot in Figure [1](#fig:demmel_URV){reference-type="ref"
reference="fig:demmel_URV"}. This graph shows computational times for
computing certain matrix factorizations versus matrix size on a standard
desktop PC. In particular, look at the times for column-pivoted QR (red
solid line) and for unpivoted QR (red dashed line). The asymptotic flop
counts of these two algorithms are identical in the dominant term. Yet
the unpivoted factorization can easily be blocked, which means that it
executes one order of magnitude faster than the pivoted one.

## The powerURV algorithm {#sec:powerURV}

There is a simple randomized algorithm for computing a rank-revealing
factorization of a matrix that perfectly illustrates the power of
randomization to reduce communication. Our starting point is an
algorithm proposed by . Given an $m\times n$ matrix $\bm{\mathsf{A}}$,
typically with $m \geq n$, it proceeds as follows.

1.  Draw an $n\times n$ matrix $\bm{\mathsf{\Omega}}$ from a standard
    normal distribution.

2.  Perform an unpivoted QR factorization of $\bm{\mathsf{\Omega}}$ so
    that $[\bm{\mathsf{V}},\sim] = \texttt{qr}(\bm{\mathsf{\Omega}})$.

3.  Perform an unpivoted QR factorization of
    $\bm{\mathsf{A}}\bm{\mathsf{V}}$ so that
    $[\bm{\mathsf{U}},\bm{\mathsf{R}}] = \texttt{qr}(\bm{\mathsf{A}}\bm{\mathsf{V}})$.

Observe that the purpose of steps (1) and (2) is simply to generate a
matrix $\bm{\mathsf{V}}$ whose columns serve as a "random" orthonormal
basis. It is easily verified that the matrices $\bm{\mathsf{U}}$,
$\bm{\mathsf{R}}$, and $\bm{\mathsf{V}}$ satisfy $$\begin{equation}
\label{eq:basicURV}
\bm{\mathsf{A}} = \bm{\mathsf{U}}\bm{\mathsf{R}}\bm{\mathsf{V}}^{*}.
\end{equation}$$ The factorization
([\[eq:basicURV\]](#eq:basicURV){reference-type="ref"
reference="eq:basicURV"}) is rank-revealing in theory (see and ), but it
does not reveal the rank particularly well in practice.

The cost to compute
([\[eq:basicURV\]](#eq:basicURV){reference-type="ref"
reference="eq:basicURV"}) is dominated by the cost to perform two
unpivoted QR factorizations, and one matrix--matrix multiplication.
(Simulating the random matrix $\bm{\mathsf{\Omega}}$ requires only
$O(n^{2})$ flops.)

To improve the rank-revealing ability of the factorization, one can
incorporate a small number of power iteration steps
[@2018_martinsson_powerurv], so that step (2) in the recipe gets
modified to $$\begin{equation}
\label{eq:powerURV}
[\bm{\mathsf{V}},\sim] = \texttt{qr}((\bm{\mathsf{A}}^{*}\bm{\mathsf{A}})^{q}\bm{\mathsf{\Omega}}),
\end{equation}$$ where $q$ is a small integer. In practice, $q=1$ or
$q=2$ is often enough to dramatically improve the accuracy of the
computation. Of course, incorporating power iteration increases the cost
of the procedure by adding $2q$ additional matrix--matrix
multiplications. When the singular values of $\bm{\mathsf{A}}$ decay
rapidly, reorthonormalization in between each application of
$\bm{\mathsf{A}}$ is sometimes required to avoid loss of accuracy due to
floating-point arithmetic.

Algorithm [\[alg:powerURV\]](#alg:powerURV){reference-type="ref"
reference="alg:powerURV"} summarizes the techniques introduced in this
::: {#alg:powerURV .algorithm}
**Algorithm. powerURV.**

```
**Input:** Target matrix $\mtx{A} \in \F^{m \times n}$ for $m\geq n$, power parameter $q$
**Output:** Orthonormal matrices $\mtx{U} \in \F^{m\times n}$ and $\mtx{V} = \F^{n\times n}$, and upper triangular $\mtx{R} \in \F^{n\times n}$ such that
$\mtx{A} = \mtx{U}\mtx{R}\mtx{V}^{*}$

function powerURV($\mtx{A}$, $q$):
    Draw a test matrix $\mtx{\Omega} \in \F^{n\times n}$ from a standard normal distribution
    $[\mtx{V},\sim] = qr_econ(\bigl(\mtx{A}^{*}\mtx{A}\bigr)^{q}\mtx{\Omega})$    # Unpivoted QR
    $[\mtx{U},\mtx{R}] = qr_econ(\mtx{A}\mtx{V})$    # Unpivoted QR
```
:::


section. The method is simple, and easy to code. It requires far more
flops than traditional methods for computing rank-revealing
factorizations, yet it is faster in practice. For instance, if $m=n$ and
$q=2$, then powerURV requires $\approx 5n^{3}$ flops versus $0.5\,n^{3}$
flops for CPQR, but Figure [1](#fig:demmel_URV){reference-type="ref"
reference="fig:demmel_URV"} shows that powerURV is still faster. This is
noteworthy, since powerURV with $q=2$ does a *far* better job at
revealing the numerical rank of $\bm{\mathsf{A}}$ than CPQR, as shown in
Figure [2](#fig:powerURVaccuracy){reference-type="ref"
reference="fig:powerURVaccuracy"}. See for details.

:::: algorithm
::: algorithmic
Target matrix $\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$ for
$m\geq n$, power parameter $q$ Orthonormal matrices
$\bm{\mathsf{U}} \in \mathbb{F}^{m\times n}$ and
$\bm{\mathsf{V}} = \mathbb{F}^{n\times n}$, and upper triangular
$\bm{\mathsf{R}} \in \mathbb{F}^{n\times n}$ such that
$\bm{\mathsf{A}} = \bm{\mathsf{U}}\bm{\mathsf{R}}\bm{\mathsf{V}}^{*}$

Draw a test matrix $\bm{\mathsf{\Omega}} \in \mathbb{F}^{n\times n}$
from a standard normal distribution
$[\bm{\mathsf{V}},\sim] = \texttt{qr\_econ}(\bigl(\bm{\mathsf{A}}^{*}\bm{\mathsf{A}}\bigr)^{q}\bm{\mathsf{\Omega}})$
Unpivoted QR
$[\bm{\mathsf{U}},\bm{\mathsf{R}}] = \texttt{qr\_econ}(\bm{\mathsf{A}}\bm{\mathsf{V}})$
Unpivoted QR
:::
::::

<figure id="fig:demmel_URV">
<div class="picture">
<p>(120,53) (05,02)<span><embed src="Pics/fig_powerURVspeed.pdf"
style="width:115mm" /></span> (00,04) (60,00)<span><span
class="math inline"><em>n</em></span></span></p>
</div>
<figcaption>Computational times required for column-pivoted QR (CPQR)
and unpivoted QR (QR) of an <span
class="math inline"><em>n</em> × <em>n</em></span> real matrix using
MATLAB on an Intel i7-8700k CPU. We see that the unpivoted factorization
is an order of magnitude faster, despite having the identical asymptotic
flop count. The graph also shows the times required for the randomized
rank-revealing factorization described in Section <a
href="#sec:demmelURV" data-reference-type="ref"
data-reference="sec:demmelURV">16.4</a>, executed both on a CPU (solid
lines) and an Nvidia Titan V GPU (dashed lines).</figcaption>
</figure>

<figure id="fig:powerURVaccuracy">
<div class="picture">
<p>(123,45) (-3,04)<span><embed src="Pics/fig_powerURVaccuracy_fast.pdf"
style="width:40mm" /></span> (40,04)<span><embed
src="Pics/fig_powerURVaccuracy_S.pdf" style="width:40mm" /></span>
(83,04)<span><embed src="Pics/fig_powerURVaccuracy_BIE.pdf"
style="width:40mm" /></span> (18,00)<span>(a)</span>
(58,00)<span>(b)</span> (98,00)<span>(c)</span></p>
</div>
<figcaption>The rank-revealing ability of CPQR and powerURV for
different values of the power parameter <span
class="math inline"><em>q</em></span>, as discussed in Section <a
href="#sec:powerURV" data-reference-type="ref"
data-reference="sec:powerURV">16.3</a>. The error <span
class="math inline"><em>e</em><sub><em>k</em></sub> = ∥A − U(:,1 : <em>k</em>)R(1 : <em>k</em>,  : )V<sup>*</sup>∥</span>
(see (<a href="#eq:rankreveal" data-reference-type="ref"
data-reference="eq:rankreveal">[eq:rankreveal]</a>)), is plotted versus
<span class="math inline"><em>k</em></span> for three different matrices
<span class="math inline">A</span> of size <span
class="math inline">400 × 400</span>. The black lines plot the
theoretical minimal values <span
class="math inline"><em>σ</em><sub><em>k</em> + 1</sub></span>. (a) A
matrix whose singular values decay rapidly; we see that all methods
perform well. (b) A matrix whose singular values plateau; we see that
CPQR performs poorly, and so does the randomized method unless powering
is used. (c) A discretized boundary integral operator whose singular
values decay slowly; we again see the high precision of powerURV for
<span class="math inline"><em>q</em> = 1</span> and <span
class="math inline"><em>q</em> = 2</span>.</figcaption>
</figure>

## Computing a rank-revealing factorization of an $n\times n$ matrix in less than $O(n^{3})$ operations {#sec:demmelURV}

The basic version of the randomized URV factorization algorithm
described in Section [16.3](#sec:powerURV){reference-type="ref"
reference="sec:powerURV"} was originally proposed by for purposes
loftier than practical acceleration. Indeed, randomization allows us to
exploit fast matrix--matrix multiplication primitives to design
accelerated algorithms for other NLA problems, such as constructing
rank-revealing factorizations. The main point of this research is that,
whenever the fast matrix--matrix multiplication is stable, the
computation of a rank-revealing factorization is stable too.

Let us be more precise. embark from the observation that there exist
algorithms [^1] for multiplying two $n\times n$ matrices using
$O(n^{\omega})$ flops, where $\omega < 3$. Once such an algorithm is
available, one can stably perform a whole range of other standard matrix
operations at the same asymptotic complexity. The idea is to apply a
divide-and-conquer approach that moves essentially all flops into the
matrix--matrix multiplication. This approach turns out to be relatively
straightforward for decompositions that do not reveal the numerical rank
such as the unpivoted QR factorization. It is harder to implement,
however, for (pivoted) rank-revealing factorizations.

## Classical column-pivoted QR {#sec:classicCPQR}

The powerURV algorithm described in Section
[16.3](#sec:powerURV){reference-type="ref" reference="sec:powerURV"} can
be very effective, but it operates on the whole matrix at once, and it
cannot be used to compute a partial factorization. In the remainder of
this section, we describe algorithms that build a rank-revealing
factorization incrementally. These methods enjoy the property that the
factorization can be halted once a specified tolerance has been met.

We start off this discussion by reviewing a classical (deterministic)
method for computing a column-pivoted QR factorization. This material is
elementary, but the discussion serves to set up a notational framework
that lets us describe the randomized version succinctly in Section
[16.6](#sec:randCPQR){reference-type="ref" reference="sec:randCPQR"}.
Suppose that we are given an $m\times n$ matrix $\bm{\mathsf{A}}$ with
$m\geq n$. We seek a factorization of the form $$\begin{equation}
\label{eq:CPQRbasic}
\begin{array}{cccccccccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{Q}} & \bm{\mathsf{R}} & \bm{\mathsf{\Pi}}^{*},\\
m\times n && m\times n & n\times n & n\times n
\end{array}
\end{equation}$$ where $\bm{\mathsf{Q}}$ has orthonormal columns, where
$\bm{\mathsf{\Pi}}$ is a permutation matrix, and where $\bm{\mathsf{R}}$
is upper-triangular with diagonal elements that decay in magnitude so
that
$|\bm{\mathsf{R}}(1,1)| \geq |\bm{\mathsf{R}}(2,2)| \geq |\bm{\mathsf{R}}(3,3)| \geq \cdots$.
The factors are typically built through a sequence of steps, where
$\bm{\mathsf{A}}$ is driven to upper-triangular form one column at a
time.

To be precise, we start by forming the matrix
$\bm{\mathsf{A}}_{0} = \bm{\mathsf{A}}$. Then we proceed using the
iteration formula
$$\bm{\mathsf{A}}_{j} = \bm{\mathsf{Q}}_{j}^{*}\bm{\mathsf{A}}_{j-1}\bm{\mathsf{\Pi}}_{j},$$
where $\bm{\mathsf{\Pi}}_{j}$ is a permutation matrix that swaps the
$j$th column of $\bm{\mathsf{A}}_{j-1}$ with the column in
$\bm{\mathsf{A}}_{j-1}(:,j:n)$ that has the largest magnitude, and where
$\bm{\mathsf{Q}}_{j}$ is a Householder reflector that zeros out all
elements below the diagonal in the $j$th column of
$\bm{\mathsf{A}}_{j-1}\bm{\mathsf{\Pi}}_{j}$; see Figure
[3](#fig:classicCPQR){reference-type="ref" reference="fig:classicCPQR"}.
Once the process concludes, the relation
[\[eq:CPQRbasic\]](#eq:CPQRbasic){reference-type="eqref"
reference="eq:CPQRbasic"} holds for
$$\bm{\mathsf{Q}}   = \bm{\mathsf{Q}}_{n}\bm{\mathsf{Q}}_{n-1}\bm{\mathsf{Q}}_{n-2}\dots\bm{\mathsf{Q}}_{1},\quad
\bm{\mathsf{R}}   = \bm{\mathsf{A}}_{n},\quad
\bm{\mathsf{\Pi}} = \bm{\mathsf{\Pi}}_{n}\bm{\mathsf{\Pi}}_{n-1}\bm{\mathsf{\Pi}}_{n-2}\dots\bm{\mathsf{\Pi}}_{1}.$$

This algorithm is well understood, and it is ubiquitous in numerical
computations. For exotic matrices, it can produce factorizations that
are quite far from optimal [@1966_kahan_NLA Sec. 5], but it typically
works very well for many tasks. For instance, it serves for revealing
the numerical rank of a matrix or for solving an ill-conditioned linear
system. However, a serious drawback to this algorithm is that it
fundamentally consists of a sequence of $n-1$ steps (or $n$ steps if
$m > n$), where a large part of the matrix is updated in each step.

<figure id="fig:classicCPQR">
<table>
<tbody>
<tr>
<td style="text-align: center;"><embed src="Pics/fig_qr0.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_qr1.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_qr2.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_qr3.pdf"
style="width:20mm" /></td>
</tr>
<tr>
<td style="text-align: center;"><span
class="math inline">A<sub>0</sub> = A</span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>1</sub> = Q<sub>1</sub><sup>*</sup>A<sub>0</sub>Π<sub>1</sub></span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>2</sub> = Q<sub>2</sub><sup>*</sup>A<sub>1</sub>Π<sub>2</sub></span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>3</sub> = Q<sub>3</sub><sup>*</sup>A<sub>2</sub>Π<sub>3</sub></span></td>
</tr>
</tbody>
</table>
<figcaption>The figure shows the sparsity pattern of a <span
class="math inline">4 × 3</span> matrix <span
class="math inline">A</span> as it is driven to upper-triangular form in
the column-pivoted QR factorization algorithm described in Section <a
href="#sec:classicCPQR" data-reference-type="ref"
data-reference="sec:classicCPQR">16.5</a>. The process takes three steps
in this case, and step <span class="math inline"><em>j</em></span>
involves the application of a permutation matrix <span
class="math inline">Π<sub><em>j</em></sub></span> from the right, and by
a Householder reflector <span
class="math inline">Q<sub><em>j</em></sub></span> from the
left.</figcaption>
</figure>

## A randomized algorithm for computing a CPQR decomposition {#sec:randCPQR}

Our next objective is to recast the algorithm for computing a CPQR
decomposition that was introduced in Section
[16.5](#sec:classicCPQR){reference-type="ref"
reference="sec:classicCPQR"} so that it works with "panels" of $b$
contiguous columns, as shown in Figure
[4](#fig:randCPQR){reference-type="ref" reference="fig:randCPQR"}. The
difficulty is to find a set of $b$ pivot vectors without updating the
matrix between each selection. Fortuitously, the randomized algorithm
for interpolatory decomposition
(Section [13.4](#sec:randID){reference-type="ref"
reference="sec:randID"}) is well adapted for this task. Indeed, a set of
$b$ columns that forms a good basis for the column space also forms a
good set of pivot vectors.

To be specific, let us describe how to pick the first group of $b$ pivot
columns for an $m\times n$ matrix $\bm{\mathsf{A}}$. Adapting the ideas
in Section [13.4](#sec:randID){reference-type="ref"
reference="sec:randID"}, we draw a Gaussian random matrix
$\bm{\mathsf{\Omega}}$ of size $(b+p) \times m$, where $p$ is a small
oversampling parameter. We form a sample matrix
$$\begin{array}{ccccccccccc}
\bm{\mathsf{Y}} &=& \bm{\mathsf{\Omega}} & \bm{\mathsf{A}}, \\
(b+p)\times n && (b+p) \times m & m\times n
\end{array},$$ and then we execute $b$ steps of column-pivoted QR on the
matrix $\bm{\mathsf{Y}}$ (either Householder or Gram--Schmidt is fine
for this step). The resulting $b$ pivot columns turn out to be good
pivot columns for $\bm{\mathsf{A}}$ as well. Once these $b$ columns have
been moved to the front of $\bm{\mathsf{A}}$, we perform a local CPQR
factorization of this panel. We update the remaining $n-b$ columns using
the computed Householder reflectors.

We could then proceed using exactly the same method a second time: draw
a $(b+p) \times (m-b)$ Gaussian random matrix $\bm{\mathsf{\Omega}}$,
form a $(b+p) \times (n-b)$ sample matrix $\bm{\mathsf{Y}}$, perform
classical CPQR on $\bm{\mathsf{Y}}$, and so on. However, there is a
shortcut. We can update the sample matrix that was used in the first
step, which renders the overhead cost induced by randomization almost
negligible [@2017_blockQR_ming_article Sec. 4].

<figure id="fig:randCPQR">
<table>
<tbody>
<tr>
<td style="text-align: center;"><embed src="Pics/fig_blockqr0.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_blockqr1.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_blockqr2.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_blockqr3.pdf"
style="width:20mm" /></td>
</tr>
<tr>
<td style="text-align: center;"><span
class="math inline">A<sub>0</sub> = A</span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>1</sub> = Q<sub>1</sub><sup>*</sup>A<sub>0</sub>Π<sub>1</sub></span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>2</sub> = Q<sub>2</sub><sup>*</sup>A<sub>1</sub>Π<sub>2</sub></span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>3</sub> = Q<sub>3</sub><sup>*</sup>A<sub>2</sub>Π<sub>3</sub></span></td>
</tr>
</tbody>
</table>
<figcaption>The sparsity pattern of a matrix <span
class="math inline">A</span> consisting of <span
class="math inline">4 × 3</span> blocks, each of size <span
class="math inline">3 × 3</span>, as it undergoes the blocked version of
the Householder QR algorithm described in Section <a
href="#sec:randCPQR" data-reference-type="ref"
data-reference="sec:randCPQR">16.6</a>. Each matrix <span
class="math inline">Q<sub><em>j</em></sub></span> is a product of three
Householder reflectors. The difficulty in building an algorithm of this
type is to find groups of pivot vectors <em>before</em> applying the
corresponding Householder reflectors.</figcaption>
</figure>

<figure id="fig:HQRRPspeedup" data-latex-placement="b">
<div class="picture">
<p>(100,54) (07,05)<span><embed src="Pics/fig_HQRRP_comparison_new.pdf"
style="width:87mm" /></span> (40,00)<span><em>Matrix size <span
class="math inline"><em>n</em></span>.</em></span> (00,17)</p>
</div>
<figcaption>Speedup of the randomized algorithm for computing a
column-pivoted QR decomposition described in Section <a
href="#sec:randCPQR" data-reference-type="ref"
data-reference="sec:randCPQR">16.6</a>, relative to LAPACK’s faster
routine (dgeqp3) as implemented in the Intel MKL library (version
11.2.3), running on an Intel Xeon E5-2695 v3 processor.</figcaption>
</figure>

Extensive numerical work has demonstrated dramatic acceleration over
deterministic algorithms. Figure
[5](#fig:HQRRPspeedup){reference-type="ref"
reference="fig:HQRRPspeedup"} draws on data from that illustrates the
acceleration over a state-of-the-art software implementation of the
classical CPQR method. Computer experiments also show that the
randomized scheme chooses pivot columns whose quality is almost
indistinguishable from those chosen by traditional pivoting, in the
sense that the relation
[\[eq:rankreveal\]](#eq:rankreveal){reference-type="eqref"
reference="eq:rankreveal"} holds to about the same accuracy. (However,
the diagonal entries of $\bm{\mathsf{R}}$ do not strictly decay in
magnitude across the block boundaries.)

To understand the behavior of the algorithm, it is helpful to think
about two extreme cases. In the first, suppose that the singular values
of $\bm{\mathsf{A}}$ decay very rapidly. Here, the analysis in Section
[11](#sec:random-rangefinder){reference-type="ref"
reference="sec:random-rangefinder"} can be modified to show that, for
any $j$, the first $j$ pivot columns chosen by the randomized algorithm
is likely to span the column space nearly as well as the optimal set of
$j$ columns. Therefore, they are excellent pivot vectors. At the other
extreme, suppose that the singular values of $\bm{\mathsf{A}}$ hardly
decay at all. In this case, the randomized method may pick a completely
different set of pivot vectors than the deterministic method, but this
outcome is unproblematic because we can take any group of columns as
pivot vectors. Of course, the interesting cases are intermediate between
these two extremes. It turns out that the randomized methods work well
regardless of how rapidly the singular values decay. For a detailed
analysis, see
[@2017_blockQR_ming_article; @2017_gu_langou; @2015_gu_melgaard; @2018_gu_flipflop].

::: remark
**Remark 51** (History). *Finding a blocked version of the CPQR
algorithm described in Section
[16.5](#sec:classicCPQR){reference-type="ref"
reference="sec:classicCPQR"} has remained an open challenge in NLA for
some time
[@1998_quintana_panelRRQR1; @2015_demmel_communication_avoiding]. The
randomized technique described in this section was introduced in , while
the updating technique was described in . For full details, see
[@2015_blockQR_SISC] and [@2017_blockQR_ming_article].*
:::

## A randomized algorithm for computing a URV decomposition {#sec:randUTV}

In this section, we describe an incremental randomized algorithm for
computing the URV factorization
([\[eq:fullfactintro\]](#eq:fullfactintro){reference-type="ref"
reference="eq:fullfactintro"}). Let us recall that for
$\bm{\mathsf{A}} \in \mathbb{F}^{m\times n}$, with $m \geq n$, this
factorization takes the form $$\begin{equation}
\label{eq:UTVrepeat}
\begin{array}{cccccccccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{U}} & \bm{\mathsf{R}} & \bm{\mathsf{V}}^{*},\\
m\times n && m\times n & n\times n & n\times n
\end{array}
\end{equation}$$ where $\bm{\mathsf{U}}$ and $\bm{\mathsf{V}}$ have
orthonormal columns and where $\bm{\mathsf{R}}$ is upper-triangular.

The algorithm we describe is blocked, and executes efficiently on modern
computing platforms. It is similar in speed to the randomized CPQR
described in Section [16.6](#sec:randCPQR){reference-type="ref"
reference="sec:randCPQR"}, and it shares the advantage that the
decomposition is built incrementally so that the process can be stopped
once a requested accuracy has been met. However, the URV factorization
offers compelling advantages: (1) It is almost as good at revealing the
numerical rank as the SVD (unlike the CPQR). (2) The URV factorization
provides us with orthonormal basis vectors for both the column and the
row spaces. (3) The off-diagonal entries of $\bm{\mathsf{R}}$ are very
small in magnitude. (4) The diagonal entries of $\bm{\mathsf{R}}$ form
excellent approximations to the singular values of $\bm{\mathsf{A}}$.

The randomized algorithm for computing a URV factorization we present
follows the same algorithmic template as the randomized CPQR described
in Section [16.6](#sec:randCPQR){reference-type="ref"
reference="sec:randCPQR"}. It drives $\bm{\mathsf{A}}$ to
upper-triangular form one block at a time, but it replaces the
permutation matrices $\bm{\mathsf{\Pi}}_{j}$ in the CPQR with general
unitary matrices $\bm{\mathsf{V}}_{j}$. Using this increased freedom, we
can obtain a factorization where all diagonal blocks are themselves
diagonal matrices and all off-diagonal elements have small magnitude.
Figure [6](#fig:blockUTV){reference-type="ref" reference="fig:blockUTV"}
summarizes the process.

<figure id="fig:blockUTV">
<div class="picture">
<p>(125,40) (000,04)<span><embed src="Pics/fig_blockUTV0.pdf"
style="width:24mm" /></span> (032,04)<span><embed
src="Pics/fig_blockUTV1.pdf" style="width:24mm" /></span>
(064,04)<span><embed src="Pics/fig_blockUTV2.pdf"
style="width:24mm" /></span> (096,04)<span><embed
src="Pics/fig_blockUTV3.pdf" style="width:24mm" /></span>
(007,00)<span><span class="math inline">A<sub>0</sub> = A</span></span>
(035,00)<span><span
class="math inline">A<sub>1</sub> = U<sub>1</sub><sup>*</sup>A<sub>0</sub>V<sub>1</sub></span></span>
(067,00)<span><span
class="math inline">A<sub>2</sub> = U<sub>2</sub><sup>*</sup>A<sub>1</sub>V<sub>2</sub></span></span>
(099,00)<span><span
class="math inline">A<sub>3</sub> = U<sub>3</sub><sup>*</sup>A<sub>2</sub>V<sub>3</sub></span></span></p>
</div>
<figcaption>Sparsity pattern of a matrix being driven to
upper-triangular form in the randomized URV factorization algorithm
described in Section <a href="#sec:randUTV" data-reference-type="ref"
data-reference="sec:randUTV">16.7</a>; see Figure <a
href="#fig:randCPQR" data-reference-type="ref"
data-reference="fig:randCPQR">4</a>. The matrices <span
class="math inline">U<sub><em>i</em></sub></span> and <span
class="math inline">V<sub><em>i</em></sub></span> are now more general
unitary transforms (consisting in the bulk of Householder reflectors).
The entries shown as gray are nonzero, but are typically very small in
magnitude.</figcaption>
</figure>

To provide details on how the algorithm works, suppose that we are given
an $m\times n$ matrix $\bm{\mathsf{A}}$ and a block size $b$. In the
first step of the process, our objective is then to build unitary
matrices $\bm{\mathsf{U}}_{1}$ and $\bm{\mathsf{V}}_{1}$ such that
$$\bm{\mathsf{A}} = \bm{\mathsf{U}}_{1}\bm{\mathsf{A}}_{1}\bm{\mathsf{V}}_{1}^{*},$$
where $\bm{\mathsf{A}}_{1}$ has the block structure
$$\bm{\mathsf{A}}_{1} =
\left[\begin{array}{cc}
\bm{\mathsf{A}}_{1,11} & \bm{\mathsf{A}}_{1,12} \\
\bm{\mathsf{0}}        & \bm{\mathsf{A}}_{1,22}
\end{array}\right] =
\raisebox{-9mm}{\includegraphics[height=18mm]{Pics/fig_UTV_step1.png}},$$
so that the $b\times b$ matrix $\bm{\mathsf{A}}_{1,11}$ is diagonal and
the entries of $\bm{\mathsf{A}}_{1,12}$ are small in magnitude. To build
$\bm{\mathsf{V}}_{1}$, we use the randomized power iteration described
in Section [11.6](#sec:rrf-subspace){reference-type="ref"
reference="sec:rrf-subspace"} to find a basis that approximately spans
the same space as the top $b$ right singular vectors of
$\bm{\mathsf{A}}$. To be precise, we form the sample matrix
$$\bm{\mathsf{Y}} = \bm{\mathsf{\Omega}}\bm{\mathsf{A}}\bigl(\bm{\mathsf{A}}^{*}\bm{\mathsf{A}}\bigr)^{q},$$
where $\bm{\mathsf{\Omega}}$ is a Gaussian random matrix of size
$b \times m$ and where $q$ is a parameter indicating the number of steps
of power iteration taken. We then perform an unpivoted QR factorization
of the *rows* of $\bm{\mathsf{Y}}$ to form a matrix
$\widetilde{\bm{\mathsf{V}}}$ whose first $b$ columns form an
orthonormal basis for the column space of $\bm{\mathsf{Y}}$. We then
apply $\widetilde{\bm{\mathsf{V}}}$ from the right to form the matrix
$\bm{\mathsf{A}}\widetilde{\bm{\mathsf{V}}}$, and we perform an
unpivoted QR factorization of the first $b$ columns of
$\bm{\mathsf{A}}\widetilde{\bm{\mathsf{V}}}$. This results in a new
matrix
$$\widetilde{\bm{\mathsf{A}}} = \bigl(\widetilde{\bm{\mathsf{U}}}\bigr)^{*}\,\bm{\mathsf{A}}\,\widetilde{\bm{\mathsf{V}}}$$
that has the block structure $$\widetilde{\bm{\mathsf{A}}} =
\left[\begin{array}{cc}
\widetilde{\bm{\mathsf{A}}}_{1,11} & \widetilde{\bm{\mathsf{A}}}_{1,12} \\
\bm{\mathsf{0}}                & \widetilde{\bm{\mathsf{A}}}_{1,22}
\end{array}\right] =
\raisebox{-9mm}{\includegraphics[height=18mm]{Pics/fig_UTV_step0.png}}.$$
The top left $b\times b$ block $\widetilde{\bm{\mathsf{A}}}_{1,11}$ is
upper-triangular, and the bottom left block is zero. The entries of
$\widetilde{\bm{\mathsf{A}}}_{1,12}$ are typically small in magnitude.
Next, we compute a full SVD of the block
$\widetilde{\bm{\mathsf{A}}}_{1,11}$:
$$\widetilde{\bm{\mathsf{A}}}_{1,11} = \widehat{\bm{\mathsf{U}}}\bm{\mathsf{D}}_{11}\widehat{\bm{\mathsf{V}}}^{*}.$$
This step is inexpensive because $\widetilde{\bm{\mathsf{A}}}_{1,11}$
has size $b\times b$, where $b$ is small. As a final step, we form the
transformation matrices $$\bm{\mathsf{U}}_{1} =
\widetilde{\bm{\mathsf{U}}}
\left[\begin{array}{cc}
\widehat{\bm{\mathsf{U}}} & \bm{\mathsf{0}} \\
\bm{\mathsf{0}} & \bm{\mathsf{I}}_{m-b}
\end{array}\right],
\qquad\mbox{and}\qquad
\bm{\mathsf{V}}_{1} =
\widetilde{\bm{\mathsf{V}}}
\left[\begin{array}{cc}
\widehat{\bm{\mathsf{V}}} & \bm{\mathsf{0}} \\
\bm{\mathsf{0}} & \bm{\mathsf{I}}_{n-b}
\end{array}\right]$$ and set
$$\bm{\mathsf{A}}_{1} = \bm{\mathsf{U}}_{1}^{*}\bm{\mathsf{A}}\bm{\mathsf{V}}_{1}.$$
The result of this process is that the diagonal entries of
$\bm{\mathsf{D}}_{11}$ typically give accurate approximations to the
first $b$ singular values of $\bm{\mathsf{A}}$, and the "remainder"
matrix $\bm{\mathsf{A}}_{1,22}$ has spectral norm that is similar to
$\sigma_{b+1}$. Thus,
$$\|\bm{\mathsf{A}}_{1,22}\| \approx \inf\{\|\bm{\mathsf{A}} - \bm{\mathsf{B}}\|\,\colon\,\bm{\mathsf{B}}\mbox{ has rank }b\}.$$
This process corresponds to the first step in Figure
[6](#fig:blockUTV){reference-type="ref" reference="fig:blockUTV"}. The
succeeding iterations execute the same process, at each step working on
the remaining lower-right part of the matrix that has not yet been
driven to upper-triangular form. We refer to and for details.

From a theoretical point of view, the first step of the URV
factorization described is well understood since it is mathematically
equivalent to the randomized power iteration described in Section
[11.6](#sec:rrf-subspace){reference-type="ref"
reference="sec:rrf-subspace"}. We observe that the first $b$ columns of
$\bm{\mathsf{U}}_{1}$ do a better job of spanning the column space of
$\bm{\mathsf{A}}$ than the first $b$ columns of $\bm{\mathsf{V}}_{1}$ do
for spanning the row space; the reason for this asymmetry is that by
forming the product $\bm{\mathsf{A}}\bm{\mathsf{V}}_{1}$, we in effect
perform an additional step of the power iteration
[@2018_martinsson_powerurv Sec. 6].

An important feature of the method described in this section is that it
is incremental, and it can be halted once a given computational
tolerance has been met. This feature has been a key competitive
advantage of the column-pivoted QR decomposition, and it is often cited
as the motivation for using CPQR. The method described in this section
has almost all the advantages of the randomized Householder CPQR
factorization (it is blocked, it is incremental, and it executes very
fast in practice), while resulting in a factorization that is far closer
to optimal in revealing the rank.

::: remark
**Remark 52** (Related work). *The idea of loosening the requirements on
the factors in a rank-revealing factorization and searching for a
decomposition such as
([\[eq:UTVrepeat\]](#eq:UTVrepeat){reference-type="ref"
reference="eq:UTVrepeat"}) is well explored in the literature
[@1999_hansen_UTVtools; @1994_stewart_UTV; @1995_elden_downdate_UTV] and
[@1998_stewart_volume1]. Deterministic techniques for computing the URV
decomposition are described in
[@1999_hansen_UTVtools; @1999_stewart_QLP]; these algorithms combine
some of the appealing qualities of the SVD (high accuracy in revealing
the rank) with some of the appealing qualities of CPQR (the possibility
of halting the execution once a requested tolerance has been met).
However, they were not blocked, and therefore they were subject to the
same liabilities as deterministic algorithms for computing the SVD and
the CPQR. A more recent use of randomization in this context is
described in [@2018_gu_flipflop].*
:::

# General linear solvers {#sec:linear-solve}

Researchers are currently exploring randomized algorithms for solving
linear systems, such as $$\begin{equation}
\label{eq:Ax=b}
\bm{\mathsf{A}}\bm{\mathsf{x}} = \bm{\mathsf{b}},
\end{equation}$$ where $\bm{\mathsf{A}}$ is a given coefficient matrix
and $\bm{\mathsf{b}}$ is a given vector. This section describes a few
probabilistic approaches for solving
[\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref" reference="eq:Ax=b"}. For
the most part, we restrict our attention to the case where
$\bm{\mathsf{A}}$ is square and the system is consistent, but we will
also touch on linear regression problems. Research on randomized linear
solvers has not progressed as rapidly as some other areas of randomized
NLA, so the discussion here is more preliminary than other parts of this
survey.

## Background: Iterative solvers

It is important to keep in mind that existing iterative solvers often
work exceptionally well. Whenever $\bm{\mathsf{A}}$ is well-conditioned
or, more generally, whenever its spectrum is "clustered," Krylov solvers
such as the conjugate gradient (CG) algorithm or GMRES tend to converge
very rapidly. For practical purposes, the cost of solving
[\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref" reference="eq:Ax=b"} is
no larger than the cost of a handful of matrix--vector multiplications
with $\bm{\mathsf{A}}$. In terms of speed, it is very difficult to beat
these techniques. Consequently, we focus on the cases where known
iterative methods converge slowly and where we cannot deploy standard
preconditioners to resolve the problem.

Having limited ourselves to this situation, the choice of solver for
[\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref" reference="eq:Ax=b"} will
depend on properties of the coefficient matrix: Is it dense or sparse?
Does it fit in RAM? Do we have access to individual matrix entries? Can
we apply $\bm{\mathsf{A}}$ to a vector? We will consider several of
these environments.

## Accelerating solvers based on dense matrix factorizations {#sec:parkeretc}

As it happens, one of the early examples of randomization in NLA was a
method for accelerating the solution of a dense linear system
[\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref" reference="eq:Ax=b"}.
observed that we can precondition a linear system by left and right
multiplying the coefficient matrix by random unitary matrices
$\bm{\mathsf{U}}$ and $\bm{\mathsf{V}}$. With probability $1$, we can
solve the resulting system $$\begin{equation}
\label{eq:Ax=b_precond}
\bigl(\bm{\mathsf{U}}\bm{\mathsf{A}}\bm{\mathsf{V}}^{*}\bigr)\,\bigl(\bm{\mathsf{V}}\bm{\mathsf{x}}\bigr) = \bm{\mathsf{U}}\bm{\mathsf{b}}
\end{equation}$$ by Gaussian elimination *without pivoting*. More
precisely, Parker proved that, almost surely, blocked Gaussian
elimination will not encounter a degenerate diagonal block.

Blocked Gaussian elimination without pivoting is substantially faster
than ordinary Gaussian elimination for two reasons: matrix operations
are more efficient than vector operations on modern computers, and we
avoid the substantial communication costs that arise when we search for
pivots. (Section [16.2](#sec:blocking){reference-type="ref"
reference="sec:blocking"} contains more discussion about blocking.)
Parker also observed that structured random matrices (such as the
randomized trigonometric transforms from
Section [9.3](#sec:srtt){reference-type="ref" reference="sec:srtt"})
allow us to perform the preconditioning step at lower cost than the
subsequent Gaussian elimination procedure.

inspired many subsequent papers, including
[@2014_li_random_butterfly_pivoting; @2017_trogdon_random_butterfly; @DGHL12:Communication-Optimal-Parallel; @2017_baboulin_GPU]
and [@2017_pan_randomized_gaussian_elim]. Another related direction
concerns the smoothed analysis of Gaussian elimination undertaken in
[@SST06:Smoothed-Analysis].

As we saw in Section [16](#sec:full){reference-type="ref"
reference="sec:full"}, randomization can be used to accelerate the
computation of rank-revealing factorizations of the matrix
$\bm{\mathsf{A}}$. In this context, randomness allows us to block the
factorization method, which increases its practical speed, even though
the overall arithmetic cost remains at $O(n^3)$. Randomized
rank-revealing factorizations are ideal for solving ill-conditioned
linear systems because they allow the user to stabilize the computation
by avoiding subspaces associated with small singular values.

For instance, suppose that we have computed a singular value
decomposition (SVD): $$\bm{\mathsf{A}} =
\bm{\mathsf{U}}\bm{\mathsf{D}}\bm{\mathsf{V}}^{*} =
\sum_{j=1}^{n}\sigma_{j}\,\bm{\mathsf{u}}_{j}\bm{\mathsf{v}}_{j}^{*}.$$
Let us introduce a truncation parameter $\varepsilon$ and ignore all
singular modes where $\sigma_j \leq \varepsilon$. Then the stabilized
solution to [\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref"
reference="eq:Ax=b"} is
$$\bm{\mathsf{x}}_{\varepsilon} = \sum_{j \,:\, \sigma_{j} > \varepsilon}\frac{1}{\sigma_{j}}\,\bm{\mathsf{v}}_{j}\bm{\mathsf{u}}_{j}^{*}\bm{\mathsf{b}}.$$
By allowing the residual to take a nonzero value, we can ensure that
$\bm{\mathsf{x}}_\varepsilon$ does not include large components that
contribute little toward satisfying the original equation. The
randomized URV decomposition, described in
Section [16.7](#sec:randUTV){reference-type="ref"
reference="sec:randUTV"}, can also be used for stabilization, and we can
compute it much faster than an SVD.

::: remark
**Remark 53** (Are rank-revealing factorizations needed?). *In some
applications, computing a rank-revealing factorization is overkill for
purposes of solving the linear system
([\[eq:Ax=b\]](#eq:Ax=b){reference-type="ref" reference="eq:Ax=b"}). In
particular, if we compute an unpivoted QR decomposition of
$\bm{\mathsf{A}}$, then it is easy to block both the factorization and
the solve stages so that very high speed is attained. This process is
provably backwards stable, which is sometimes all that is needed. (In
practice, partially pivoted LU can often be used in an analogous manner,
despite being theoretically unstable.)*

*In contrast, when the actual entries of the computed solution
$\bm{\mathsf{x}}_{\mathrm{approx}}$ matter (as opposed to the value of
$\bm{\mathsf{A}}\bm{\mathsf{x}}_{\mathrm{approx}}$), a stabilized solver
is generally preferred. As a consequence, column-pivoted QR is often
cited as a method of choice for ill-conditioned problems in situations
where an SVD is not affordable.*
:::

::: remark
**Remark 54** (Strassen accelerated solvers). *We saw in Section
[16.4](#sec:demmelURV){reference-type="ref" reference="sec:demmelURV"}
that randomization has enabled us to compute a rank-revealing
factorization of an $n\times n$ matrix in less than $\mathcal{O}(n^{3})$
operations. The idea was to use randomized preconditioning as in
([\[eq:Ax=b_precond\]](#eq:Ax=b_precond){reference-type="ref"
reference="eq:Ax=b_precond"}), and then accelerate an unpivoted
factorization of the resulting coefficient matrix using fast algorithms
for the matrix-matrix multiplication such as Strassen
[@2007_demmel_fast_linear_algebra_is_stable]. This methodology can of
course be immediately applied to the task of solving ill-conditioned
linear systems. For improved numerical stability, a few steps of power
iteration can be incorporated to this approach; see
[\[eq:powerURV\]](#eq:powerURV){reference-type="eqref"
reference="eq:powerURV"}.*
:::

## Sketch and precondition {#sec:sketchtoprecond}

Another approach to preconditioning is to look for a random
transformation of the linear system that makes an iterative linear
solver converge more quickly. Typically, these preconditioning
transforms need to cluster the eigenvalues of the matrix.

The most successful example of this type of randomized preconditioning
does not concern square systems, but rather highly overdetermined
least-squares problems. See
Section [10.5](#sec:sketchandprecond){reference-type="ref"
reference="sec:sketchandprecond"} *et seq.* for a discussion of this
idea. This type of randomized preconditioning can greatly enhance the
robustness and power of "asynchronous" solvers for
communication-constrained environments [@2015_avron_revisiting]. Related
techniques for kernel ridge regression are described in . For linear
systems involving high-dimensional tensors, see .

For square linear systems, the search for randomized preconditioners has
been less fruitful.
Section [18](#sec:sparse-cholesky){reference-type="ref"
reference="sec:sparse-cholesky"} outlines the main success story.
Nevertheless, techniques already at hand can be very helpful for solving
linear systems in special situations, which we illustrate with a small
example.

Consider the task of solving
[\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref" reference="eq:Ax=b"} for
a positive-definite (PD) coefficient matrix $\bm{\mathsf{A}}$. In this
environment, the iterative method of choice is the conjugate gradient
(CG) algorithm [@1952_hestenes_CG]. A detailed convergence analysis for
CG is available; for example, see [@1997_trefethen_bau Sec. 38]. In a
nutshell, CG converges rapidly when the eigenvalues of $\bm{\mathsf{A}}$
are clustered, as in
Figure [7](#fig:eval_distributions_for_CG){reference-type="ref"
reference="fig:eval_distributions_for_CG"}(a). Therefore, our task is to
find a matrix $\bm{\mathsf{M}}$ for which $\bm{\mathsf{M}}^{-1}$ can be
applied rapidly to vectors and for which
$\bm{\mathsf{M}}^{-1/2} \bm{\mathsf{A}} \bm{\mathsf{M}}^{-1/2}$ has a
tightly clustered spectrum.

In a situation where $\bm{\mathsf{A}}$ has a few eigenvalues that are
larger than the others
(Figure [7](#fig:eval_distributions_for_CG){reference-type="ref"
reference="fig:eval_distributions_for_CG"}(b)), randomized algorithms
for low-rank approximation provide excellent preconditioners. For
instance, we can use the randomized Nyström method
(Section [14](#sec:nystrom){reference-type="ref"
reference="sec:nystrom"}) to compute an approximation $$\begin{equation}
\label{eq:spd_exact_lowrank}
\bm{\mathsf{A}} \approx \bm{\mathsf{U}}\bm{\mathsf{D}}\bm{\mathsf{U}}^{*}
\end{equation}$$ where $\bm{\mathsf{D}} \in \mathbb{R}_{+}^{k\times k}$
is a diagonal matrix whose entries hold approximations to the largest
$k$ eigenvalues of $\bm{\mathsf{A}}$, and where
$\bm{\mathsf{U}} \in \mathbb{F}^{m\times k}$ is an orthonormal matrix
holding the corresponding approximate eigenvectors. We then form a
preconditioner for $\bm{\mathsf{A}}$ by setting
$$\bm{\mathsf{M}} = (1/\alpha)\,\bm{\mathsf{U}}\bm{\mathsf{D}}\bm{\mathsf{U}}^{*} + \bigl(\bm{\mathsf{I}} - \bm{\mathsf{U}}\bm{\mathsf{U}}^{*}\bigr).$$
It is trivial to invert $\bm{\mathsf{M}}$ because
$\bm{\mathsf{M}}^{-1} = \alpha\bm{\mathsf{U}}\bm{\mathsf{D}}^{-1}\bm{\mathsf{U}}^{*} + \bigl(\bm{\mathsf{I}} - \bm{\mathsf{U}}\bm{\mathsf{U}}^{*}\bigr)$.
Now, if
[\[eq:spd_exact_lowrank\]](#eq:spd_exact_lowrank){reference-type="eqref"
reference="eq:spd_exact_lowrank"} captured the top $k$ eigenmodes of
$\bm{\mathsf{A}}$ exactly, then the preconditioned coefficient matrix
$\bm{\mathsf{M}}^{-1/2}\bm{\mathsf{A}}\bm{\mathsf{M}}^{-1/2}$ would have
the same eigenvectors as $\bm{\mathsf{A}}$, but with the top $k$
eigenvalues replaced by $\alpha$ and the remaining eigenvalues
unchanged. By setting $\alpha = \lambda_{k}$, say, the spectrum of
$\bm{\mathsf{M}}^{-1/2}\bm{\mathsf{A}}\bm{\mathsf{M}}^{-1/2}$ would
become far more tightly clustered. In reality, the columns of
$\bm{\mathsf{U}}$ do not exactly align with the eigenvectors of
$\bm{\mathsf{A}}$. Even so, the accuracy will be good for the
eigenvectors associated with the top eigenvalues, which is what matters.

<figure id="fig:eval_distributions_for_CG">
<table>
<tbody>
<tr>
<td style="text-align: left;"><embed src="Pics/fig_spectra_1.pdf"
style="width:30mm" /></td>
<td style="text-align: left;"><embed src="Pics/fig_spectra_2.pdf"
style="width:30mm" /></td>
<td style="text-align: left;"><embed src="Pics/fig_spectra_3.pdf"
style="width:30mm" /></td>
<td style="text-align: left;"><embed src="Pics/fig_spectra_4.pdf"
style="width:30mm" /></td>
</tr>
<tr>
<td style="text-align: left;">(a)</td>
<td style="text-align: left;">(b)</td>
<td style="text-align: left;">(c)</td>
<td style="text-align: left;">(d)</td>
</tr>
</tbody>
</table>
<figcaption>The eigenvalues of four different PD matrices that all have
condition number <span class="math inline">10</span> (since <span
class="math inline"><em>λ</em><sub>max</sub> = 1</span> and <span
class="math inline"><em>λ</em><sub>min</sub> = 0.1</span>). As discussed
in Section <a href="#sec:sketchtoprecond" data-reference-type="ref"
data-reference="sec:sketchtoprecond">17.3</a>, the difficulty of solving
the corresponding linear systems using conjugate gradients differ
significantly between these cases. (a) For this matrix, CG converges in
two iterations, without the need for preconditioners. (b) When the
spectrum has some large outliers, the randomized preconditioner outlined
in Section <a href="#sec:sketchtoprecond" data-reference-type="ref"
data-reference="sec:sketchtoprecond">17.3</a> works well. (c,d) Finding
randomized preconditioners for matrices with spectra like these remains
an open research problem.</figcaption>
</figure>

## The randomized Kaczmarz method and its relatives {#sec:rk}

The Kaczmarz method is an iterative algorithm for solving linear systems
that is typically used for large, overdetermined problems with
inconsistent equations. Randomized variants of the Kaczmarz method have
received a lot of attention in recent years, in part because of close
connections to stochastic gradient descent (SGD) algorithms for solving
least-squares problems.

To explain the idea, consider a (possibly inconsistent) linear system
$$\begin{equation}
 \label{eq:Ax=b_kacz}
\bm{\mathsf{A}}^* \bm{\mathsf{x}} \approx \bm{\mathsf{b}}
\quad\text{where}\quad
\bm{\mathsf{A}}^* \in \mathbb{F}^{m \times n}.
\end{equation}$$ The basic Kaczmarz algorithm starts with an initial
guess $\bm{\mathsf{x}}_0 \in \mathbb{F}^n$ for the solution. At each
iteration $t$, we select a new index $j = j(t) \in \{1, \dots, m\}$, and
we make the update $$\begin{equation}
 \label{eq:kacz_trad}
\bm{\mathsf{x}}_{t+1} = \bm{\mathsf{x}}_t + \frac{\bm{\mathsf{b}}(j) - \langle \bm{\mathsf{A}}(:, j), \, \bm{\mathsf{x}}_t \rangle}{\Vert \bm{\mathsf{A}}(:,j) \Vert^2} \bm{\mathsf{A}}(:, j).
\end{equation}$$ The
rule [\[eq:kacz_trad\]](#eq:kacz_trad){reference-type="eqref"
reference="eq:kacz_trad"} has a simple interpretation: it ensures that
$\bm{\mathsf{x}}_{t+1}$ is the closest point to $\bm{\mathsf{x}}_t$ in
the hyperplane containing solutions to the linear equation determined by
the $j$th equation in the system.

In implementing this method, we must choose a control mechanism that
determines the next index. A simple and robust approach is to cycle
through the rows consecutively; that is, $j(t) = t \bmod m$. Another
effective, but expensive, option is to select the equation with the
largest violation.

The randomized Kaczmarz (RK) algorithm uses a probabilistic control
mechanism instead. This kind of approach also has a long history, and it
is useful in cases where cyclic control is ineffective. The RK method
has received renewed attention owing to work of . They proposed sampling
each $j(t)$ independently at random, with the probability of choosing
the $i$th equation proportional to the squared $\ell_2$ norm of the
$i$th column of $\bm{\mathsf{A}}$. They proved that this version of RK
converges linearly with a rate determined by the (Demmel) condition
number of the matrix $\bm{\mathsf{A}}$. Later, it was recognized that
this approach is just a particular instantiation of SGD for the
least-squares
problem [\[eq:Ax=b_kacz\]](#eq:Ax=b_kacz){reference-type="eqref"
reference="eq:Ax=b_kacz"}. See , which draws on results from .

There are many subsequent papers that have built on the RK approach for
solving inconsistent linear systems. observed that related ideas can be
used to design randomized Gauss--Seidel and randomized Jacobi
iterations. studied a blocked version of the RK algorithm, which is
practically more efficient for many of the same reasons that other
blocked algorithms work well
(Section [16.2](#sec:blocking){reference-type="ref"
reference="sec:blocking"}).

observed that the RK algorithm is a particular type of iterative
sketching method. Based on this connection, they proposed a
generalization. At each iteration, we draw an independent random
embedding $\bm{\mathsf{S}}_t \in \mathbb{F}^{\ell \times m}$. The next
iterate is chosen by solving the least-squares problem
$$\begin{equation}
\label{eq:Ax=b_kacz_proj0}
\bm{\mathsf{x}}_{t} = \arg\min\nolimits_{\bm{\mathsf{y}}} \Vert  \bm{\mathsf{x}}_{t-1} - \bm{\mathsf{y}}  \Vert^2
\quad\text{subject to}\quad
\bm{\mathsf{S}}_t \bm{\mathsf{A}}^* \bm{\mathsf{y}} = \bm{\mathsf{S}}_t \bm{\mathsf{b}}.
\end{equation}$$ The idea is to choose the dimension $\ell$ sufficiently
small that the sketched least-squares problem can be solved explicitly
using a direct method (e.g., QR factorization). This flexibility leads
to algorithms that converge more rapidly in practice because the sketch
$\bm{\mathsf{S}}_t$ can mix equations instead of just sampling. Later,
showed that this procedure can be accelerated to achieve rates that
depend on the *square root* of an appropriate condition number; see also
.

# Linear solvers for graph Laplacians {#sec:sparse-cholesky}

In this section we describe the randomized algorithm,
[SparseCholesky]{.smallcaps}, which can efficiently solve a linear
system whose coefficient matrix is a graph Laplacian matrix. Up to a
polylogarithmic factor, this method achieves the minimum possible
runtime and storage costs. This algorithm has the potential to
accelerate many types of computations involving graph Laplacian
matrices.

The [SparseCholesky]{.smallcaps} algorithm was developed by and further
refined by . These approaches are based on earlier work from Dan
Spielman's group, notably the paper of . The presentation here is
adapted from .

## Overview

We begin with a high-level approach for solving Laplacian linear
systems. The basic idea is to construct a preconditioner using a
randomized variant of the incomplete Cholesky method. Then we can solve
the original linear system by means of the preconditioned conjugate
gradient (PCG) algorithm.

### Approximate solutions to the Poisson problem

Let $\bm{\mathsf{L}} \in \mathbb{H}_n(\mathbb{R})$ be the Laplacian
matrix of a weighted, loop-free, undirected graph on the vertex set
$V = \{1, \dots, n\}$. We write $m$ for the number of edges in the
graph, i.e., the sparsity of the graph. For simplicity, we will also
assume that the graph is *connected*; equivalently,
$\ker(\bm{\mathsf{L}}) = \operatorname{span}\{ \mathbf{1} \}$ where
$\mathbf{1} \in \mathbb{R}^n$ is the vector of ones. See
Section [7.4](#sec:graph-sparsification){reference-type="ref"
reference="sec:graph-sparsification"} for definitions. See
Figure [8](#fig:big-dipper){reference-type="ref"
reference="fig:big-dipper"} for an illustration of an unweighted,
undirected graph.

The basic goal is to find the unique solution $\bm{\mathsf{x}}_{\star}$
to the Poisson problem
$$\bm{\mathsf{L}} \bm{\mathsf{x}} = \bm{\mathsf{f}}
\quad\text{where}\quad
\mathbf{1}^* \bm{\mathsf{f}} = {0}
\quad\text{and}\quad
\mathbf{1}^* \bm{\mathsf{x}} = {0}.$$ For a parameter $\varepsilon> 0$,
we can relax this requirement by asking instead for an approximate
solution $\bm{\mathsf{x}}_{\varepsilon}$ that satisfies
$$\Vert  \bm{\mathsf{x}}_{\varepsilon} - \bm{\mathsf{x}}_{\star}  \Vert_{\bm{\mathsf{L}}}
    \leq \varepsilon\, \Vert  \bm{\mathsf{x}}_{\star}  \Vert_{\bm{\mathsf{L}}}.$$
We have written
$\Vert  \bm{\mathsf{x}}  \Vert_{\bm{\mathsf{L}}} := (\bm{\mathsf{x}}^* \bm{\mathsf{L}} \bm{\mathsf{x}})^{1/2}$
for the energy seminorm induced by the Laplacian matrix.

### Approximate Cholesky decomposition

Imagine that we can efficiently construct a sparse, approximate Cholesky
decomposition of the Laplacian matrix $\bm{\mathsf{L}}$. More precisely,
we seek a morally lower-triangular matrix $\bm{\mathsf{C}}$ that
satisfies $$\begin{equation}
 \label{eqn:approx-cholesky}
0.5 \, \bm{\mathsf{L}} \preccurlyeq\bm{\mathsf{CC}}^* \preccurlyeq 1.5 \, \bm{\mathsf{L}}
    \quad\text{where}\quad
    \texttt{nnz}(\bm{\mathsf{C}}) = O(m \log^2 n).
\end{equation}$$ In other words, there is a known permutation of the
rows that brings the matrix $\bm{\mathsf{C}}$ into lower-triangular
form. As usual, $\preccurlyeq$ is the semidefinite order, and `nnz`
returns the number of nonzero entries in a matrix.

This section describes an algorithm, called
[SparseCholesky]{.smallcaps}, that can complete the task outlined in the
previous paragraph. This algorithm is motivated by the insight that we
can produce a sparse approximation of a Laplacian matrix by random
sampling (Section [7.4](#sec:graph-sparsification){reference-type="ref"
reference="sec:graph-sparsification"}). The main challenge is to obtain
sampling probabilities without extra computation. The resulting method
can be viewed as a randomized variant of the incomplete Cholesky
factorization [@GVL13:Matrix-Computations-4ed Sec. 11.5.8].

### Preconditioning

Given the sparse, approximate Cholesky factor $\bm{\mathsf{C}}$, we can
precondition the Poisson problem: $$\begin{equation}
 \label{eqn:precond-poisson}
(\bm{\mathsf{C}}^\dagger\bm{\mathsf{L}} \bm{\mathsf{C}}^{*\dagger})(\bm{\mathsf{C}}^* \bm{\mathsf{x}}) = (\bm{\mathsf{C}}^\dagger\bm{\mathsf{f}}).
\end{equation}$$
When [\[eqn:approx-cholesky\]](#eqn:approx-cholesky){reference-type="eqref"
reference="eqn:approx-cholesky"} holds, the matrix
$\bm{\mathsf{C}}^\dagger\bm{\mathsf{L}} \bm{\mathsf{C}}^{*\dagger}$ has
condition number $\kappa \leq 3$.

Therefore, we can solve the preconditioned
system [\[eqn:precond-poisson\]](#eqn:precond-poisson){reference-type="eqref"
reference="eqn:precond-poisson"} quickly using the PCG algorithm
[@GVL13:Matrix-Computations-4ed Sec. 11.5]. If the initial iterate
$\bm{\mathsf{x}}_0 = \bm{\mathsf{0}}$, then $j$ steps of PCG produce an
iterate $\bm{\mathsf{x}}_j$ that satisfies
$$\Vert  \bm{\mathsf{x}}_j - \bm{\mathsf{x}}_{\star}  \Vert_{\bm{\mathsf{L}}}
    \leq 2 \left[ \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right]^j
    \Vert  \bm{\mathsf{x}}_{\star}  \Vert_{\bm{\mathsf{L}}}
    < 3^{1-j} \Vert  \bm{\mathsf{x}}_{\star}  \Vert_{\bm{\mathsf{L}}}.$$
As a consequence, we can achieve relative error $\varepsilon$ after
$1 + \log_3 (1/\varepsilon)$ iterations. Each iteration requires a
matrix--vector product with $\bm{\mathsf{L}}$ and the solution of a
(consistent) linear system
$\bm{\mathsf{CC}}^* \bm{\mathsf{u}} = \bm{\mathsf{y}}$. We can perform
these steps in $O(m \log^2 n)$ operations per iteration. Indeed,
$\bm{\mathsf{L}}$ has only $O(m)$ nonzero entries. The matrix
$\bm{\mathsf{C}}$ is morally triangular with $O(m \log^2 n)$ nonzero
entries, so we can apply $(\bm{\mathsf{CC}}^*)^{\dagger}$ using two
triangular solves.

::: remark
**Remark 55** (Triangular solve). *To solve a consistent linear system
$(\bm{\mathsf{CC}}^*) \bm{\mathsf{u}} = \bm{\mathsf{y}}$, we can apply
$\bm{\mathsf{C}}^\dagger$ by triangular elimination and then apply
$\bm{\mathsf{C}}^{*\dagger}$ by triangular elimination. The four
subspace theorem ensures that solution to the first problem renders the
second problem consistent too. Since
$\ker(\bm{\mathsf{CC}}^*) = \operatorname{span}\{\bm{\mathsf{1}}\}$, we
can enforce consistency numerically by removing the constant component
of the input $\bm{\mathsf{y}}$. Similarly, we can remove the constant
component of the output $\bm{\mathsf{u}}$ to ensure it belongs to the
correct space.*
:::

### Main results

The following theorem describes the performance of the
[SparseCholesky]{.smallcaps} procedure. This is the main result from .

::: {#thm:sparse-cholesky .theorem}
**Theorem 56** (SparseCholesky). *Let $\bm{\mathsf{L}}$ be the Laplacian
of a connected graph on $n$ vertices, with $m$ weighted edges. With high
probability, the [SparseCholesky]{.smallcaps} algorithm produces a
morally lower-triangular matrix $\bm{\mathsf{C}}$ that satisfies
$$0.5 \, \bm{\mathsf{L}} \preccurlyeq\bm{\mathsf{CC}}^* \preccurlyeq 1.5 \, \bm{\mathsf{L}}.$$
The matrix $\bm{\mathsf{C}}$ has $O(m \log^2 n)$ nonzero entries. The
expected running time is $O(m \log^3 n)$ operations.*
:::

In view of our discussion about PCG, we arrive at the following
statement about solving the Poisson problem.

::: corollary
**Corollary 57** (Poisson problem). *Suppose the
[SparseCholesky]{.smallcaps} algorithm delivers an approximation
$\bm{\mathsf{L}} \approx \bm{\mathsf{CC}}^*$ that
satisfies [\[eqn:approx-cholesky\]](#eqn:approx-cholesky){reference-type="eqref"
reference="eqn:approx-cholesky"}. Then we can solve each consistent
linear system $\bm{\mathsf{L}}\bm{\mathsf{x}} = \bm{\mathsf{f}}$ to
relative error $\varepsilon$ in the seminorm
$\Vert \cdot \Vert_{\bm{\mathsf{L}}}$ using at most
$1 + \log_3(1/\varepsilon)$ iterations of PCG, each with a cost of
$O(m \log^2 n)$ arithmetic operations.*
:::

### Discussion

As we have mentioned, the Poisson problem serves as a primitive for
undertaking many computations on undirected
graphs [@Ten10:Laplacian-Paradigm]. Potential applications include
clustering, analysis of random walks, and finite-element discretizations
of elliptic PDEs.

The [SparseCholesky]{.smallcaps} algorithm achieves a near-optimal
runtime and storage guarantee for the Poisson problem on a graph.
Indeed, for a general graph with $m$ edges, any algorithm must use
$O(m)$ storage and arithmetic. There is a proof that the cost
$O(m \log^{1/2}(n))$ is achievable in theory [@CKM+14:Solving-SDD], but
the resulting methods are currently impractical. Meanwhile, the
simplicity of the [SparseCholesky]{.smallcaps} method makes it a
candidate for real-world computation.

For particular classes of Laplacian matrices, existing solvers can be
very efficient. Optimized sparse direct solvers
[@2016_acta_sparse_direct_survey] work very well for small- and
medium-size problems, but they typically have superlinear scaling, which
renders them unsuitable for truly large-scale problems. Iterative
methods such as multigrid or preconditioned Krylov solvers can attain
linear complexity for important classes of problems, in particular for
sparse systems arising from the discretization of elliptic PDEs.
However, we are not aware of competing methods that provably enjoy
near-optimal complexity for all problems.

We regard the [SparseCholesky]{.smallcaps} algorithm as one of the most
dramatic examples of how randomization has the potential to accelerate
basic linear algebra computations, both in theory and in practice.

## Cholesky decomposition of a graph Laplacian

To begin our explanation of the [SparseCholesky]{.smallcaps} algorithm,
let us summarize what happens when we apply the standard Cholesky
decomposition method to a graph Laplacian.

### The Laplacian of a multigraph

For technical reasons, related to the design and analysis of the
algorithm, we need to work with multigraphs instead of ordinary graphs.
In the discussion, we will point out specific places where this
generality is important.

Consider a weighted, undirected *multigraph* $G$, defined on the vertex
set $V = \{1, \dots, n\}$. Each edge $e = \{u, v\}$ is an unordered pair
of vertices; we typically use the abbreviated notation $e = uv = vu$. We
introduce a weight function $w_G$ that assigns a positive weight to each
edge $e$ in the multigraph $G$. Since $G$ is a multigraph, there may be
many multiedges, with distinct weights, connecting the same two
vertices.

Taking some notational liberties, we will identify the multigraph $G$
with its Laplacian matrix $\bm{\mathsf{L}}$, which we express in the
form $$\begin{equation}
 \label{eqn:multi-laplacian}
\bm{\mathsf{L}} = \sum\nolimits_{e \in \bm{\mathsf{L}}} w_{\bm{\mathsf{L}}}(e) \, \bm{\mathsf{\Delta}}_{e}.
\end{equation}$$ The matrix $\bm{\mathsf{\Delta}}_{e}$ is the elementary
Laplacian on the vertex pair $(u, v)$ that composes the edge $e = uv$.
That is,
$$\bm{\mathsf{\Delta}}_e := \bm{\mathsf{\Delta}}_{uv} := (\bm{\mathsf{\delta}}_{u} - \bm{\mathsf{\delta}}_v)(\bm{\mathsf{\delta}}_{u} - \bm{\mathsf{\delta}}_v)^*
\quad\text{where $e = uv$.}$$ The sum
in [\[eqn:multi-laplacian\]](#eqn:multi-laplacian){reference-type="eqref"
reference="eqn:multi-laplacian"} takes place over all the multiedges $e$
in the multigraph $\bm{\mathsf{L}}$, so the same elementary Laplacian
may appear multiple times with different weights.

### Stars and cliques {#sec:star-clique}

To describe the Cholesky algorithm on a graph, we need to introduce a
few more concepts from graph theory. Define the *degree* and the *total
weight* of a vertex $u$ in the multigraph $\bm{\mathsf{L}}$ to be
$$\operatorname{deg}_{\bm{\mathsf{L}}}(u) := \sum\nolimits_{e=uv \in \bm{\mathsf{L}}} 1
\quad\text{and}\quad
w_{\bm{\mathsf{L}}}(u) := \sum\nolimits_{e = uv \in \bm{\mathsf{L}}} w_{\bm{\mathsf{L}}}(e).$$
In other words, the degree of $u$ is the total number of multiedges $e$
that contain $u$. The total weight of $u$ is the sum of the weights of
the multiedges $e$ that contain $u$.

Let $u$ be a fixed vertex. The *star* induced by $u$ is the Laplacian
$$\textsc{star}(u, \bm{\mathsf{L}}) := \sum\nolimits_{e = uv \in \bm{\mathsf{L}}} w_{\bm{\mathsf{L}}}(e) \, \bm{\mathsf{\Delta}}_e.$$
In words, the star includes precisely those multiedges $e$ in the
multigraph $\bm{\mathsf{L}}$ that contain the vertex $u$.

The *clique* induced by $u$ is defined implicitly as the correction that
occurs when we take the Schur
complement [\[eqn:schur-complement\]](#eqn:schur-complement){reference-type="eqref"
reference="eqn:schur-complement"} of the Laplacian with respect to the
coordinate $u$: $$\bm{\mathsf{L}} / \bm{\mathsf{\delta}}_u %
    =: (\bm{\mathsf{L}} - \textsc{star}(u, \bm{\mathsf{L}})) + \textsc{clique}(u, \bm{\mathsf{L}}).$$
Recall that $\bm{\mathsf{\delta}}_u$ is the standard basis vector in
coordinate $u$. By direct calculation, one may verify that
$$\textsc{clique}(u, \bm{\mathsf{L}}) = \frac{1}{2 w_{\bm{\mathsf{L}}}(u)} \sum\limits_{e_1 = uv_1 \in \bm{\mathsf{L}}} \sum\limits_{e_2 = uv_2 \in \bm{\mathsf{L}}}
    w_{\bm{\mathsf{L}}}(e_1) w_{\bm{\mathsf{L}}}(e_2) \bm{\mathsf{\Delta}}_{v_1v_2}.$$
Each sum takes place over all multiedges $e$ in $\bm{\mathsf{L}}$ that
contain the vertex $u$. It can be verified that the clique is also the
Laplacian of a weighted multigraph.

Figures [8](#fig:big-dipper){reference-type="ref"
reference="fig:big-dipper"}
and [9](#fig:star-clique){reference-type="ref"
reference="fig:star-clique"} contain an illustration of a (simple)
graph, along with the star and clique induced by eliminating a vertex.
In our more general setting, the edges in the star and clique would have
associated weights. These diagrams are courtesy of Richard Kueng.

<figure id="fig:big-dipper">

<figcaption><em>A combinatorial graph.</em> The Ursa Major graph with a
distinguished vertex (the star Megrez) highlighted in red.</figcaption>
</figure>

<figure id="fig:star-clique">
<table>
<tbody>
<tr>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
</tr>
</tbody>
</table>
<figcaption><em>Illustration of a star and clique in the Ursa Major
graph.</em> (a) The star induced by the red vertex consists of the three
solid black edges. (b) The clique induced by the red vertex consists of
the three solid black edges. These edges are added when the red vertex
is eliminated. </figcaption>
</figure>

### Graphs and Cholesky

With the notation introduced in
Section [18.2.2](#sec:star-clique){reference-type="ref"
reference="sec:star-clique"}, we can present the graph-theoretic
interpretation of the Cholesky algorithm as it applies to the Laplacian
$\bm{\mathsf{L}}$ of a weighted multigraph.

Define the initial Laplacian $\bm{\mathsf{S}}_0 := \bm{\mathsf{L}}$. In
each step $i = 1, 2, \dots, n$, we select a new vertex $u_i$. We extract
the associated column of the current Laplacian:
$$\bm{\mathsf{c}}_i := \frac{1}{\sqrt{(\bm{\mathsf{S}}_{i-1})_{u_i u_i}}} \bm{\mathsf{S}}_{i-1} \bm{\mathsf{\delta}}_{u_i}.$$
We compute the Schur complement with respect to the vertex $u_i$:
$$\bm{\mathsf{S}}_i := \bm{\mathsf{S}}_{i-1} / \bm{\mathsf{\delta}}_{u_i} = (\bm{\mathsf{S}}_{i-1} - \textsc{star}(u_i, \bm{\mathsf{S}}_{i-1})) + \textsc{clique}(u_i, \bm{\mathsf{S}}_{i-1}).$$
In other words, we remove the star induced by $u_i$ and replace it with
the clique induced by $u_i$. Each $\bm{\mathsf{S}}_i$ is the Laplacian
of a multigraph; it has no multiedge that contains any one of the
vertices $u_1, \dots, u_i$. Therefore, we have reduced the size of the
problem.

After $n$ steps, the Cholesky factorization is determined by a vector
$\bm{\mathsf{\pi}} = (u_1, u_2, \dots, u_n)$ that holds the chosen
indices and a matrix
$\bm{\mathsf{C}} = \begin{bmatrix} \bm{\mathsf{c}}_1 & \dots & \bm{\mathsf{c}}_n \end{bmatrix}$
such that $\bm{\mathsf{C}}(\bm{\mathsf{\pi}}, :)$ is lower-triangular.
The algorithm ensures that we have the exact decomposition
$$\bm{\mathsf{L}} = \bm{\mathsf{CC}}^*.$$ This is the (pivoted) Cholesky
factorization of the Laplacian matrix.

To choose a vertex to eliminate, the classical approach is to find a
vertex with minimum degree or with minimum total weight. Or one may
simply select one of the remaining vertices at random.

### Computational costs

The cost to compute a Cholesky factorization
$\bm{\mathsf{L}} = \bm{\mathsf{C}}\bm{\mathsf{C}}^{*}$ of a Laplacian
matrix $\bm{\mathsf{L}}$ is typically superlinear in $n$, with a
worst-case cost of $O(n^{3})$ arithmetic and $O(n^{2})$ storage. This is
the reason that $\bm{\mathsf{C}}$ is less sparse than $\bm{\mathsf{L}}$:
the clique that is introduced at an elimination step has more edges than
the star that it replaces, a phenomenon referred to as *fill-in*. The
exact growth in the number of nonzero entries depends on the sparsity
pattern of $\bm{\mathsf{L}}$ and on the chosen elimination order. For
special cases, using a nested dissection ordering can provably improve
on the worst-case estimates [@2016_acta_sparse_direct_survey]. For
instance, if $\bm{\mathsf{L}}$ results from the finite-difference or
finite-element discretization of an elliptic PDE, then
$\texttt{nnz}(\bm{\mathsf{C}}) = O(n\log(n))$ in two dimensions and
$\texttt{nnz}(\bm{\mathsf{C}}) = O(n^{4/3})$ in three.

For general graphs, one path towards improving the efficiency of the
Cholesky factorization procedure is to randomly approximate the clique
by sampling, in order to curb the fill-in.
Section [7.4](#sec:graph-sparsification){reference-type="ref"
reference="sec:graph-sparsification"} already indicates that this
innovation may be possible, provided that we can find a way to obtain
sampling probabilities.

## The [SparseCholesky]{.smallcaps} algorithm

We are now prepared to present the [SparseCholesky]{.smallcaps}
procedure, which uses randomized sampling to compute a sparse,
approximate Cholesky factorization.

### Procedure

Let $\bm{\mathsf{L}}$ be the Laplacian of a weighted multigraph on
$V = \{1, \dots, n\}$. We perform the following steps.

1.  **Preprocessing.** Split each multiedge $e = uv$ in
    $\bm{\mathsf{L}}$ into $R = \lceil 64 \log^2(\mathrm{e}n) \rceil$
    multiedges, each connecting $\{u, v\}$, and each with weight
    $w_{\bm{\mathsf{L}}}(e) / R$. The purpose of this step is to control
    the effective
    resistance [\[eqn:effective-resistance\]](#eqn:effective-resistance){reference-type="eqref"
    reference="eqn:effective-resistance"} of each multiedge at the
    outset of the algorithm. Note that this splitting results in a
    weighted multigraph, even if we begin with a simple graph.

2.  **Initialization.** Form the initial Laplacian
    $\bm{\mathsf{S}}_0 = \bm{\mathsf{L}}$ and the list of remaining
    vertices $F_0 = V$.

3.  **Iteration.** For each $i = 1, 2, \dots, n$:

    1.  **Select a vertex.** Choose a vertex $u_i$ uniformly at random
        from $F_{i-1}$. Remove this vertex from the list:
        $F_i = F_{i-1} \setminus \{ u_i \}$.

    2.  **Extract the column.** Copy the normalized $u_i$th column from
        the current Laplacian:
        $$\bm{\mathsf{c}}_i = \frac{1}{\sqrt{(\bm{\mathsf{S}}_{i-1})_{u_i u_i}}} \bm{\mathsf{S}}_{i-1} \bm{\mathsf{\delta}}_{u_i}.$$
        Set $\bm{\mathsf{c}}_i = \bm{\mathsf{0}}$ if the denominator
        equals zero.

    3.  **Sampling the clique.** Construct the Laplacian
        $\bm{\mathsf{K}}_i$ of a random sparse approximation of
        $\textsc{clique}(u_i, \bm{\mathsf{S}}_{i-1})$. We will detail
        this procedure in
        Section [18.3.2](#sec:clique-sample){reference-type="ref"
        reference="sec:clique-sample"}.

    4.  **Approximate Schur complement.** Form
        $$\bm{\mathsf{S}}_i = (\bm{\mathsf{S}}_{i-1} - \textsc{star}(u_i, \bm{\mathsf{S}}_{i-1})) + \bm{\mathsf{K}}_i.$$

4.  **Decomposition.** Collate the columns $\bm{\mathsf{c}}_i$ into the
    Cholesky factor
    $$\bm{\mathsf{C}} = \begin{bmatrix} \bm{\mathsf{c}}_1 & \dots & \bm{\mathsf{c}}_n \end{bmatrix}.$$
    Define the row permutation $\pi(i) = u_i$ for each $i$.

Once these operations are complete, $\bm{\mathsf{C}}$ is a sparse,
morally lower-triangular matrix. It is also very likely that
$\bm{\mathsf{L}} \approx \bm{\mathsf{CC}}^*$.
Theorem [56](#thm:sparse-cholesky){reference-type="ref"
reference="thm:sparse-cholesky"} makes a rigorous accounting of these
claims.

### Clique sampling {#sec:clique-sample}

The remaining question is how to construct a random approximation
$\bm{\mathsf{K}}$ of a clique $\textsc{clique}(u, \bm{\mathsf{S}})$.
Here is the procedure:

1.  **Probabilities.** Construct a probability mass $\bm{\mathsf{p}}$
    such that
    $$p(e) = \frac{w_{\bm{\mathsf{S}}}(e)}{w_{\bm{\mathsf{S}}}(u)}
    \quad\text{for each $e \in \textsc{star}(u, \bm{\mathsf{S}})$.}$$

2.  **Sampling.** For each
    $i = 1, \dots, d = \operatorname{deg}_{\bm{\mathsf{S}}}(u)$,

    1.  Draw a random multiedge $e_1 = uv_1$ from the multiedges in
        $\textsc{star}(u, \bm{\mathsf{S}})$ according to the probability
        mass $\bm{\mathsf{p}}$.

    2.  Draw a second random multiedge $e_2 = uv_2$ from the multiedges
        in $\textsc{star}(u, \bm{\mathsf{S}})$ according to the uniform
        distribution.

    3.  Form the random Laplacian matrix of a new multiedge:
        $$\bm{\mathsf{X}}_i = \frac{w_{\bm{\mathsf{S}}}(e_1) \, w_{\bm{\mathsf{S}}}(e_2)}{w_{\bm{\mathsf{S}}}(e_1) + w_{\bm{\mathsf{S}}}(e_2)} \bm{\mathsf{\Delta}}_{v_1 v_2}.$$

3.  **Approximation.** Return
    $\bm{\mathsf{K}} = \sum_{i=1}^d \bm{\mathsf{X}}_i$.

The key fact about this construction is that it produces an unbiased
estimator $\bm{\mathsf{K}}$ of the clique:
$$\operatorname{\mathbb{E}}\bm{\mathsf{K}} = \textsc{clique}(u, \bm{\mathsf{S}}).$$
Furthermore, each summand $\bm{\mathsf{X}}_i$ creates a multiedge with
uniformly bounded effective
resistance [\[eqn:effective-resistance\]](#eqn:effective-resistance){reference-type="eqref"
reference="eqn:effective-resistance"}. This property persists as the
[SparseCholesky]{.smallcaps} algorithm executes, and it ensures that the
random matrix $\bm{\mathsf{K}}$ has controlled variance. Note that the
sampling procedure can result in several edges between the same pair of
vertices, which is another reason we need the multigraph formalism.

Next, observe that the number $d$ of multiedges in $\bm{\mathsf{K}}$ is
no greater than the number $d$ of multiedges in the star that we are
removing from $\bm{\mathsf{S}}$. (For comparison, note that the full
clique has $d^2$ multiedges.) As a consequence, the clique approximation
is inexpensive to construct. Moreover, the total number of multiedges in
the Laplacian can only decrease as the [SparseCholesky]{.smallcaps}
algorithm proceeds.

### Analysis

The analysis of [SparseCholesky]{.smallcaps} is well beyond the scope of
this paper. The key technical tool is a concentration inequality for
matrix-valued martingales that was derived in
[@Oli09:Concentration-Adjacency] and [@Tro11:Freedmans-Inequality]. To
activate this result, [@KS16:Approximate-Gaussian] use the fact that the
random clique approximation is unbiased and low-variance, conditional on
previous choices made by the algorithm. The proof also relies heavily on
the fact that we eliminate a random vertex at each step of the
iteration. For the technical details,
see [@KS16:Approximate-Gaussian; @Kyn17:Approximate-Gaussian]
and [@Tro19:Matrix-Concentration-LN].

### Implementation

The [SparseCholesky]{.smallcaps} algorithm is fairly simple to describe,
but it demands some care to develop an implementation that achieves the
runtime guarantees stated in
Theorem [56](#thm:sparse-cholesky){reference-type="ref"
reference="thm:sparse-cholesky"}.

The most important point is that we need to use data structures for
weighted multigraphs. One method is to maintain the vertex--multiedge
adjacency matrix, along with a list of weights. This approach requires
sparse matrix libraries, including efficient iterators over the rows and
columns.

A secondary point is that we need fast methods for constructing finite
probability distributions and sampling from them repeatedly. See .

It is unlikely that the [SparseCholesky]{.smallcaps} procedure will fail
to produce a factor $\bm{\mathsf{C}}$ that
satisfies [\[eqn:approx-cholesky\]](#eqn:approx-cholesky){reference-type="eqref"
reference="eqn:approx-cholesky"}. Even so, it is reassuring to know that
we can detect failures. Indeed, we can estimate the extreme singular
values of the preconditioned linear
system [\[eqn:precond-poisson\]](#eqn:precond-poisson){reference-type="eqref"
reference="eqn:precond-poisson"} using the methods from
Section [6](#sec:max-eig){reference-type="ref" reference="sec:max-eig"}.

The failure probability can also be reduced by modifying the random
vertex selection rule. Instead, we can draw a random vertex whose total
weight is at most twice the average total weight of the remaining
vertices [@Kyn17:Approximate-Gaussian]. In practice, it may suffice to
use the classical elimination rules based on minimum degree or minimum
total weight.

The main shortcoming of the [SparseCholesky]{.smallcaps} procedure
arises from the initialization step, where we split each multiedge into
$O(\log^2 n)$ pieces. This step increases the storage and computation
costs of the algorithm enough to make it uncompetitive (e.g. with fast
direct Poisson solvers) for some problem instances. At present, it is
unclear whether the initialization step can be omitted or relaxed, while
maintaining the reliability and correctness of the algorithm.

# Kernel matrices in machine learning {#sec:kernel}

Randomized NLA algorithms have played a major role in developing
scalable kernel methods for large-scale machine learning. This section
contains a brief introduction to kernels. Then it treats two
probabilistic techniques that have had an impact on kernel matrix
computations: Nyström approximation by random coordinate sampling
[@WS01:Using-Nystrom] and empirical approximation by random features
[@RR08:Random-Features].

The literature on kernel methods is truly vast, so we cannot hope to
achieve comprehensive coverage within our survey. There are also many
computational considerations and learning-theoretic aspects that fall
outside the realm of NLA. Our goal is simply to give a taste of the
ideas, along with a small selection of key references.

## Kernels in machine learning

We commence with a crash course on kernels and their applications in
machine learning. The reader may refer to for a more complete treatment.

### Kernel functions and kernel matrices {#sec:kernel-functions}

Let $\mathcal{X}$ be a set, called the *input space* or *data space*.
Suppose that we acquire a finite set of observations
$\{ \bm{\mathsf{x}}_1, \dots, \bm{\mathsf{x}}_n \} \subset \mathcal{X}$.
We would like to use the observed data to perform learning tasks.

One approach is to introduce a *kernel function*:
$$k : \mathcal{X} \times \mathcal{X} \to \mathbb{F}.$$ The value
$k(\bm{\mathsf{x}}, \bm{\mathsf{y}})$ of the kernel function is
interpreted as a measure of similarity between two data points
$\bm{\mathsf{x}}$ and $\bm{\mathsf{y}}$. We can tabulate the pairwise
similarities of the observed data points in a *kernel matrix*:
$$(\bm{\mathsf{K}})_{ij} := k( \bm{\mathsf{x}}_i, \bm{\mathsf{x}}_j )
\quad\text{for $i, j = 1, \dots, n$.}$$ The kernel matrix is an analog
of the Gram matrix of a set of vectors in a Euclidean space. In
Sections [19.1.5](#sec:kpca){reference-type="ref" reference="sec:kpca"}
and [19.1.6](#sec:krr){reference-type="ref" reference="sec:krr"}, we
will explain how to use the matrix $\bm{\mathsf{K}}$ to solve some core
problems in data analysis.

The kernel function is required to be *positive definite*. That is, for
each natural number $n$ and each set
$\{\bm{\mathsf{x}}_1, \dots, \bm{\mathsf{x}}_n\} \subset \mathcal{X}$ of
observations, the associated kernel matrix
$$\bm{\mathsf{K}} = \big[ k(\bm{\mathsf{x}}_i, \bm{\mathsf{x}}_j) \big]_{i, j = 1, \dots, n} \in \mathbb{H}_n
\quad\text{is psd.}$$ In particular,
$k(\bm{\mathsf{x}}, \bm{\mathsf{x}}) \geq 0$ for all
$\bm{\mathsf{x}} \in \mathcal{X}$. The kernel must also be (conjugate)
symmetric in its arguments:
$k(\bm{\mathsf{x}}, \bm{\mathsf{y}}) = k(\bm{\mathsf{y}}, \bm{\mathsf{x}})^*$
for all $\bm{\mathsf{x}}, \bm{\mathsf{y}} \in \mathcal{X}$. These
properties mirror the properties of a Gram matrix. In
Section [19.1.3](#sec:kernel-examples){reference-type="ref"
reference="sec:kernel-examples"}, we give some examples of positive
definite kernel functions.

### The feature space

It is common to present kernel functions using the theory of
*reproducing kernel Hilbert spaces*. This approach gives an alternative
interpretation of the kernel function as the inner product defined on a
feature space. We give a very brief treatment, omitting all technical
details.

Let $\mathcal{F}$ be a Hilbert space, called the *feature space*. We
introduce a *feature map* $\Phi: \mathcal{X} \to \mathcal{F}$, which
maps a point in the input space to a point in the feature space.
Heuristically, the feature map extracts information from a data point
that is relevant for learning applications.

Under mild conditions, we can construct a positive definite kernel
function $k$ from the feature map:
$$k(\bm{\mathsf{x}}, \bm{\mathsf{y}}) = \langle  \Phi(\bm{\mathsf{x}}) , \,  \Phi(\bm{\mathsf{y}})  \rangle
\quad\text{for all $\bm{\mathsf{x}}, \bm{\mathsf{y}} \in \mathcal{X}$.}$$
In other words, the kernel function reports the inner product between
the features associated with the data points $\bm{\mathsf{x}}$ and
$\bm{\mathsf{y}}$. Conversely, a positive definite kernel always induces
a feature map into an appropriate feature space.

### Examples of kernels {#sec:kernel-examples}

Kernel methods are powerful because we can select or design a kernel
that automatically extracts relevant feature information from our data.
This approach applies in all sorts of domains, including images and text
and DNA sequences. Let us present a few kernels that commonly arise in
applications. See  for many additional examples and references.

::: example
**Example 58** (Inner product kernel). *The simplest example of a kernel
is the ordinary inner product. Let $\mathcal{X} = \mathbb{F}^d$.
Evidently,
$$k(\bm{\mathsf{x}}, \bm{\mathsf{y}}) = \langle \bm{\mathsf{x}}, \, \bm{\mathsf{y}} \rangle
\quad\text{for $\bm{\mathsf{x}}, \bm{\mathsf{y}} \in \mathbb{F}^d$}$$ is
a positive definite kernel.*
:::

::: {#ex:angsim .example}
**Example 59** (Angular similarity). *Another simple example is the
angular similarity map. Let
$\mathcal{X} = \mathbb{S}^{d-1}(\mathbb{R}) \subset \mathbb{R}^d$. This
kernel is given by the formula
$$k(\bm{\mathsf{x}}, \bm{\mathsf{y}}) = \frac{2}{\pi} \arcsin\ \langle  \bm{\mathsf{x}} , \,  \bm{\mathsf{y}}  \rangle
\quad\text{for $\bm{\mathsf{x}}, \bm{\mathsf{y}} \in \mathcal{X}$.}$$
This kernel is positive definite because of Schoenberg's
theorem [@Sch42:Positive-Definite]. We will give a short direct proof in
Example [65](#ex:angsim-rf){reference-type="ref"
reference="ex:angsim-rf"}.*
:::

::: example
**Example 60** (Polynomial kernels). *Let $\mathcal{X}$ be a subset of
$\mathbb{F}^d$. For a natural number $p$, the inhomogeneous polynomial
kernel is
$$k(\bm{\mathsf{x}}, \bm{\mathsf{y}}) = (1 + \langle \bm{\mathsf{x}}, \, \bm{\mathsf{y}} \rangle)^p.$$
This kernel is also positive definite because of Schoenberg's theorem;
see . There is a short direct proof using the Schur product theorem.*
:::

::: {#ex:rbf .example}
**Example 61** (Gaussian kernel). *An important example is the Gaussian
kernel. Let $\mathcal{X} = \mathbb{F}^d$. For a bandwidth parameter
$\sigma > 0$, define
$$k( \bm{\mathsf{x}}, \bm{\mathsf{y}} ) = \exp\left(-\frac{\Vert  \bm{\mathsf{x}} - \bm{\mathsf{y}}  \Vert^2 }{ 2 \sigma^2 } \right)
    \quad\text{for $\bm{\mathsf{x}}, \bm{\mathsf{y}} \in \mathbb{F}^d$.}$$
This kernel is positive definite because of Bochner's
theorem [@Boc33:Monotone-Funktionen]. We will give a short direct proof
in Example [66](#ex:rbf-rf){reference-type="ref"
reference="ex:rbf-rf"}.*
:::

### The kernel trick

As we have mentioned, kernels can be used for a wide range of tasks in
machine learning. state the key idea succinctly:

> Given an algorithm which is formulated in terms of a positive definite
> kernel $k$, one can construct an alternative algorithm by replacing
> $k$ with another positive definite kernel $\tilde{k}$.

In particular, any algorithm that can be formulated in terms of the
inner product kernel applies to every other kernel. That is to say, an
algorithm for Euclidean data that depends only on the Gram matrix can be
implemented with a kernel matrix instead. The next two subsections give
two specific examples of this methodology; there are many other
applications.

### Kernel PCA {#sec:kpca}

Given a set of observations in a Euclidean space, principal component
analysis (PCA) searches for orthogonal directions in which the data has
the maximum variability. The nonlinear extension, kernel PCA (KPCA), was
proposed in ; see also .

Let
$\{\bm{\mathsf{x}}_1, \dots, \bm{\mathsf{x}}_n \} \subset \mathcal{X}$
be a set of observations. For a kernel $k$ associated with a feature map
$\Phi$, construct the kernel matrix $\bm{\mathsf{K}} \in \mathbb{H}_n$
associated with the observations. For a natural number $\ell$, we
compute a truncated eigenvalue decomposition of the kernel matrix:
$$\bm{\mathsf{K}} = \sum_{i=1}^{\ell} \lambda_i \, \bm{\mathsf{u}}_i \bm{\mathsf{u}}_i^*.$$
Each unit-norm eigenvector $\bm{\mathsf{u}}_i$ determines a direction
$(n\lambda_{i})^{-1/2}\sum_{j=1}^n \bm{\mathsf{u}}_i(j) \, \Phi(\bm{\mathsf{x}}_j)$
of high variability in the feature space, called the $i$th kernel
principal component.

To find the projection of a new point $\bm{\mathsf{x}} \in \mathcal{X}$
onto the $i$th kernel principal component, we embed it into the feature
space via $\Phi(\bm{\mathsf{x}})$ and compute the inner product with the
$i$th kernel principal component. In terms of the kernel function,
$$\mathrm{PC}_i(\bm{\mathsf{x}}) :=
\frac{1}{\sqrt{n\lambda_{i}}}\sum_{j=1}^n \bm{\mathsf{u}}_i(j) \, k(\bm{\mathsf{x}}, \bm{\mathsf{x}}_j).$$
We can summarize the observation $\bm{\mathsf{x}}$ with the vector
$$(\mathrm{PC}_1(\bm{\mathsf{x}}), \dots, \mathrm{PC}_{\ell}(\bm{\mathsf{x}})) \in \mathbb{F}^\ell.$$
This representation provides a data-driven feature that can be used for
downstream learning tasks.

In practice, it is valuable to center the feature space representation
of the data, which requires a simple modification of the kernel matrix.
We also need to center each observation before computing its projection
onto the kernel principal components. See  for details.

### Kernel ridge regression {#sec:krr}

Given a set of labeled observations in a Euclidean space, ridge
regression uses regularized least-squares to model the labels as a
linear functional of the observations. The nonlinear extension of this
approach is called *kernel ridge regression* (KRR). We refer to for a
more detailed treatment, including an interpretation in terms of a
nonlinear feature map.

Let
$\{ (\bm{\mathsf{x}}_i, y_i) : i = 1, \dots, n \} \subset \mathcal{X} \times \mathbb{F}$
be a set of paired observations. For a kernel $k$, construct the kernel
matrix $\bm{\mathsf{K}} \in \mathbb{H}_n$ associated with the
observations $\bm{\mathsf{x}}_i$ (but not the numerical values $y_i$).
For a regularization parameter $\tau > 0$, the kernel ridge regression
problem takes the form
$$\underset{\bm{\mathsf{\alpha}} \in \mathbb{F}^n}{\text{minimize}} \quad \frac{1}{n} \sum_{i=1}^n \left[ y_i - (\bm{\mathsf{K}} \bm{\mathsf{\alpha}})_i \right]^2
    + \frac{\tau}{2} \bm{\mathsf{\alpha}}^* \bm{\mathsf{K}} \bm{\mathsf{\alpha}}.$$
The solution to this optimization problem is obtained by solving an
ordinary linear system:
$$( \bm{\mathsf{K}} + \tau n \, \bm{\mathsf{I}}) \bm{\mathsf{\alpha}} = \bm{\mathsf{y}}
\quad\text{where}\quad
\bm{\mathsf{y}} = (y_1, \dots, y_n).$$ Let
$\widehat{\bm{\mathsf{\alpha}}}$ be the solution to this system.

Given a new observation $\bm{\mathsf{x}} \in \mathcal{X}$, we can make a
prediction $\widehat{y} \in \mathbb{F}$ for its label via the formula
$$\widehat{y}(\bm{\mathsf{x}}) %
    := \sum_{j=1}^n \widehat{\alpha}_j \, k(\bm{\mathsf{x}}, \bm{\mathsf{x}}_j).$$
In practice, the regularization parameter $\tau$ is chosen by
cross-validation with a holdout set of the paired observations.

### The issue

Kernel methods are powerful tools for data analysis. Nevertheless, in
their native form, they suffer from two weaknesses.

First, it is very expensive to compute the kernel matrix explicitly. For
example, if points in the data space $\mathcal{X}$ have a
$d$-dimensional parameterization, we may expect that it will cost $O(d)$
arithmetic operations to evaluate the kernel a single time. Therefore,
the cost of forming the kernel matrix $\bm{\mathsf{K}}$ for $n$
observations is $O(n^2 d)$.

Second, after computing the kernel matrix $\bm{\mathsf{K}}$, it remains
expensive to perform the linear algebra required by kernel methods. Both
KPCA and KRR require $O(n^3)$ operations if we use direct methods.

The poor computational profile of kernel methods limits our ability to
use them directly for large-scale data applications.

### The solution

Fortunately, there is a path forward. To implement kernel methods, we
simply need to *approximate* the kernel matrix [@SS01:Learning-Kernels
Sec. 10.2]. Surprisingly, using the approximation often results in
*better* learning outcomes than using the exact kernel matrix. Even a
poor approximation of the kernel can suffice to achieve near-optimal
performance, both in theory and in practice
[@2013_bach_sharp; @RCR15:Less-More; @RCR17:FALKON-Optimal]. Last,
working with a structured approximation of the kernel can accelerate the
linear algebra computations dramatically.

Randomized algorithms provide several effective tools for approximating
kernel matrices. Since we pay a steep price for each kernel evaluation,
we need to develop algorithms that explicitly control this cost. The
rest of this section describes two independent approaches. In
Section [19.2](#sec:nystrom-kernel){reference-type="ref"
reference="sec:nystrom-kernel"}, we present coordinate sampling
algorithms for Nyström approximations, while
Section [19.3](#sec:random-features){reference-type="ref"
reference="sec:random-features"} develops the method of random features.

::: remark
**Remark 62** (Function approximation). *Let us remark that
approximation of the kernel matrix is usually incidental to the goals of
learning theory. In many applications, such as KRR, we actually need to
approximate a function on the input space. The sampling complexity of
the latter task may be strictly lower than the complexity of
approximating the full kernel matrix. We cannot discuss this issue in
detail because it falls outside the scope of NLA.*
:::

## Coordinate Nyström approximation of kernel matrices {#sec:nystrom-kernel}

One way to approximate a kernel matrix is to form a Nyström
decomposition with respect to a judiciously chosen coordinate subspace.
A natural idea is to draw these coordinates at random. This basic
technique was proposed by .

Coordinates play a key role here because we only have access to
individual entries of the kernel matrix. There is no direct way to
compute a matrix--vector product with the kernel matrix, so we cannot
easily apply the more effective constructions of random embeddings
(e.g., Gaussians or sparse maps or SRTTs). Indeed, kernel computation is
the primary setting where coordinate sampling is a practical idea.

### Coordinate Nyström approximation

Suppose that
$\{\bm{\mathsf{x}}_1, \dots, \bm{\mathsf{x}}_n \} \subset \mathcal{X}$
is a collection of observations. Let $\bm{\mathsf{K}}$ be the psd kernel
matrix associated with some kernel function $k$.

Given a set $I \subseteq \{1, \dots, n\}$ consisting of $r$ indices, we
can form a psd Nyström approximation of the kernel matrix:
$$\bm{\mathsf{K}}\langle I \rangle := \bm{\mathsf{K}}(:, I) \, \bm{\mathsf{K}}(I, I)^\dagger\, \bm{\mathsf{K}}(I, :).$$
This matrix is equivalent to the Nyström
decomposition [\[eqn:nystrom-def\]](#eqn:nystrom-def){reference-type="eqref"
reference="eqn:nystrom-def"} with respect to a test matrix
$\bm{\mathsf{X}}$ whose range is
$\operatorname{span}\{ \bm{\mathsf{\delta}}_i : i \in I \}$.

To obtain $\bm{\mathsf{K}}\langle I \rangle$, the basic cost is $nr$
kernel evaluations, which typically require $\mathcal{O}(nrd)$
operations for a $d$-dimensional input space $\mathcal{X}$. We typically
do not form the pseudoinverse directly, but rather use the factored form
of the Nyström approximation for downstream calculations.

For kernel problems, it is common to regularize the coordinate Nyström
approximation. One approach replaces the core matrix
$\bm{\mathsf{K}}(I, I)$ with its truncated eigenvalue decomposition
before computing the pseudoinverse; for example, see . The RSVD
algorithm (Section [11.2](#sec:rsvd){reference-type="ref"
reference="sec:rsvd"}) has been proposed for this purpose
[@LBKL15:Large-Scale-Nystrom]. When it is computationally feasible, we
recommend taking a truncated eigenvalue decomposition of the full
Nyström approximation $\bm{\mathsf{K}}\langle I \rangle$, rather than
just the core; see .

### Greedy selection of coordinates

Recall from Section [14](#sec:nystrom){reference-type="ref"
reference="sec:nystrom"} that the error
$\bm{\mathsf{K}} / I := \bm{\mathsf{K}} - \bm{\mathsf{K}}\langle I \rangle$
in the Nyström decomposition is simply the Schur complement of
$\bm{\mathsf{K}}$ with respect to the coordinates in $I$.

This connection suggests that we should use a pivoted Cholesky method or
a pivoted QR algorithm to select the coordinates $I$. These techniques
lead to Nyström approximations with superior learning performance; for
example, see [@FS01:Efficient-SVM; @BJ05:Predictive-Low-Rank]
and [@2013_bach_sharp]. Unfortunately the $\mathcal{O}(n^2 r)$ cost is
prohibitive in applications. See  for some randomized strategies that
can reduce the expense.

### Ridge leverage scores

A natural approach to selecting the coordinate set $I$ is to perform
randomized sampling. To describe these approaches, we need to take a
short detour.

Fix a regularization parameter $\tau > 0$. Consider the smoothed
projector:
$$h_{\tau}(\bm{\mathsf{K}}) := \bm{\mathsf{K}} (\bm{\mathsf{K}} + \tau n \bm{\mathsf{I}})^{-1}.$$
The number of effective degrees of freedom at regularization level
$\tau$ is
$$\nu_{\mathrm{eff}} := \operatorname{trace}h_{\tau}(\bm{\mathsf{K}}). %$$
The maximum marginal number of degrees of freedom at regularization
level $\tau$ is
$$\nu_{\mathrm{mof}} := n \cdot \max_{i = 1, \dots, n} (h_{\tau}(\bm{\mathsf{K}}))_{ii}.$$
Observe that $\nu_{\mathrm{eff}} \leq \nu_{\mathrm{mof}}$. The statistic
$\nu_{\mathrm{mof}}$ is analogous with the coherence that appears in our
initial discussion of coordinate sampling
(Section [9.6](#sec:coord-embed){reference-type="ref"
reference="sec:coord-embed"}).

The *ridge leverage scores* at regularization level $\tau$ are the
(normalized) diagonal entries of the smoothed projector:
$$p_i = \frac{(h_{\tau}(\bm{\mathsf{K}}))_{ii}}{\nu_{\mathrm{eff}}}
\quad\text{for $i = 1, \dots, n$.}$$ Evidently, $(p_1, \dots, p_n)$ is a
probability distribution. The ridge leverage scores and related
quantities are expensive to compute directly, but there are efficient
algorithms for approximating them well. These approximations suffice for
applications. See
Section [19.2.5](#sec:ridge-leverage){reference-type="ref"
reference="sec:ridge-leverage"} for more discussion.

::: remark
**Remark 63** (History). *identified the core role of the smoothed
projector for KRR. proposed the definition of the ridge leverage scores
and described a simple method for approximating them. At present, the
most practical algorithm for approximating ridge leverage scores appears
in . recognize that ridge leverage scores also have relevance for KPCA.*
:::

### Uniform sampling {#uniform-sampling-1}

The simplest way to select a set $I$ of $r$ coordinates for the Nyström
approximation $\bm{\mathsf{K}}\langle I \rangle$ is to draw the set
uniformly at random. Although this approach seems naïve, it can be
surprisingly effective in practice. The main failure mode occurs when
there are a few significant observations that make outsize contributions
to the kernel matrix; uniform sampling is likely to miss these
influential data points.

proves that we can achieve optimal learning guarantees for KRR with a
uniformly sampled Nyström approximation. It suffices that the number $r$
of coordinates is proportional to $\nu_{\mathrm{mof}} \log n$. A similar
result holds for KPCA.

Nyström approximation with uniform coordinate sampling was proposed by .
The theoretical and numerical performance of this approach has been
studied in many subsequent works, including
[@KMT12:Sampling-Methods; @Git13:Topics-Randomized; @2013_bach_sharp]
and [@RCR17:FALKON-Optimal].

### Sampling with ridge leverage scores {#sec:ridge-leverage}

Suppose that we have computed an approximation of the ridge leverage
score distribution. We can construct a coordinate set $I$ for the
Nyström approximation $\bm{\mathsf{K}}\langle I \rangle$ by sampling $r$
coordinates independently at random from the ridge leverage score
distribution. Properly implemented, this method is unlikely to miss
influential observations.

prove that we can achieve optimal learning guarantees for KRR by ridge
leverage score sampling. It suffices that the number $r$ of sampled
coordinates is proportional to $\nu_{\mathrm{eff}} \log n$. This bound
improves over the uniform sampling bound. give related theoretical
results for KPCA.

Effective algorithms for estimating ridge leverage scores are based on
multilevel procedures that sequentially improve the ridge leverage score
estimates. The basic idea is to start with a small uniform sample of
coordinates, which we use to approximate the smoothed projector for a
very large regularization parameter $\tau_0$. From this smoothed
projector, we estimate the ridge leverage scores at level $\tau_0$. We
then sample a larger set of coordinates non-uniformly using the
approximate ridge leverage score distribution at level $\tau_0$. These
samples allow us to approximate the smoothed projector at level
$\tau_1 = \mathrm{const.} \tau_0$ for a constant smaller than one. We
obtain an estimate for the ridge leverage score distribution at level
$\tau_1$. This process is repeated. In this way, the sampling and the
matrix approximation are intertwined. See [@MM17:Recursive-Sampling] and
[@RCCR18:Fast-Leverage].

provide empirical evidence that ridge leverage score sampling is more
efficient than uniform sampling for KPCA, including the cost of the
ridge leverage score approximations. Likewise, report empirical evidence
that ridge leverage score sampling is more efficient than uniform
sampling for KRR.

## Random features approximation of kernels {#sec:random-features}

A second approach to kernel approximation is based on the method of
empirical approximation
(Section [7](#sec:matrix-mc){reference-type="ref"
reference="sec:matrix-mc"}). This technique constructs a random rank-one
matrix that serves as an unbiased estimator for the kernel matrix. By
averaging many copies of the estimator, we can obtain superior
approximations of the kernel matrix. The individual rank-one components
are called *random features*.

proposed the idea of using empirical approximation for kernels arising
in Gaussian process regression. Later, and developed empirical
approximations for translation-invariant kernels and Mercer kernels, and
they coined the term "random features." Our presentation is based on an
abstract formulation of the random feature method from ; see also .

In this subsection, we introduce the idea of a random feature map, along
with some basic examples. We explain how to use random feature maps to
construct empirical approximations of a kernel matrix, and we give a
short analysis. Afterwards, we summarize two randomized NLA methods for
improving the computational profile of random features.

### Random feature maps

In many cases, a kernel function $k$ on a domain $\mathcal{X}$ can be
written as an expectation, and we can exploit this representation to
obtain empirical approximations of the kernel matrix.

Let $\mathcal{W}$ be a probability space equipped with a probability
measure $\rho$. Assume that there is a bounded function
$$\psi : \mathcal{X} \times \mathcal{W} \to \{ z \in \mathbb{C}: \vert z \vert \leq b \}$$
with the reproducing property $$\begin{equation}
 \label{eqn:repro-prop}
k(\bm{\mathsf{x}}, \bm{\mathsf{y}}) = \int \psi(\bm{\mathsf{x}}; \bm{\mathsf{w}}) \, \psi(\bm{\mathsf{y}}; \bm{\mathsf{w}})^* \, \rho(\mathrm{d}{\bm{\mathsf{w}}})
\quad\text{for all $\bm{\mathsf{x}}, \bm{\mathsf{y}} \in \mathcal{X}$.}
\end{equation}$$ The star ${}^*$ denotes the conjugate of a complex
number. We call $(\psi, \rho)$ a *random feature map* for the kernel
function $k$. As we will see in
Section [19.3.2](#sec:rf-kernel){reference-type="ref"
reference="sec:rf-kernel"}, a kernel that admits a random feature map
must be positive definite.

It is not obvious that we can equip kernels of practical interest with
random feature maps, so let us offer a few concrete examples.

::: {#ex:innerprodkernel .example}
**Example 64** (The inner product kernel).

*There are many ways to construct a random feature map for the inner
product kernel on $\mathbb{F}^d$. One simple example is
$$\psi(\bm{\mathsf{x}}; \bm{\mathsf{w}}) = \langle  \bm{\mathsf{x}} , \,  \bm{\mathsf{w}}  \rangle
\quad\text{with $\bm{\mathsf{w}} \sim \textsc{normal}(\bm{\mathsf{0}}, \bm{\mathsf{I}}_d)$.}$$
To check that this map satisfies the reproducing
property [\[eqn:repro-prop\]](#eqn:repro-prop){reference-type="eqref"
reference="eqn:repro-prop"}, just note that $\bm{\mathsf{w}}$ is
isotropic:
$\operatorname{\mathbb{E}}[ \bm{\mathsf{ww}}^* ] = \bm{\mathsf{I}}$.
This formulation is closely related to the theory of random embeddings
(Sections [8](#sec:gauss){reference-type="ref" reference="sec:gauss"}
and [9](#sec:dimension-reduction){reference-type="ref"
reference="sec:dimension-reduction"}) and to approximate matrix
multiplication (Section [7.3](#sec:approx-mtx-mult){reference-type="ref"
reference="sec:approx-mtx-mult"}).*
:::

::: {#ex:angsim-rf .example}
**Example 65** (Angular similarity kernel).

*For the angular similarity map defined in
Example [59](#ex:angsim){reference-type="ref" reference="ex:angsim"}, we
can construct a random feature map using an elegant fact from geometry.
Indeed, the function
$$\psi(\bm{\mathsf{x}}; \bm{\mathsf{w}}) = \operatorname{sgn}\ \langle \bm{\mathsf{x}}, \, \bm{\mathsf{w}} \rangle
%
\quad\text{with $\bm{\mathsf{w}} \sim \textsc{uniform}(\mathbb{S}^{d-1}(\mathbb{R}))$}$$
gives a random feature map for the angular similarity kernel. As a
consequence, the angular similarity kernel is positive definite.*
:::

::: {#ex:rbf-rf .example}
**Example 66** (Translation-invariant kernels).

*A kernel function $k$ on $\mathbb{F}^d$ is called
*translation-invariant* if it has the form
$$k( \bm{\mathsf{x}}, \bm{\mathsf{y}} ) = \phi( \bm{\mathsf{x}} - \bm{\mathsf{y}} )
\quad\text{for all $\bm{\mathsf{x}}, \bm{\mathsf{y}} \in \mathbb{F}^d$.}$$
A classic result from analysis, Bochner's
theorem [@Boc33:Monotone-Funktionen], gives a characterization of these
kernels. A kernel is continuous, positive definite and
translation-invariant if and only if it is the Fourier transform of a
positive probability measure $\rho$ on $\mathbb{F}^d$:
$$\phi( \bm{\mathsf{x}} - \bm{\mathsf{y}} ) = c_{\phi} \int \mathrm{e}^{\mathrm{i}\langle \bm{\mathsf{x}}, \, \bm{\mathsf{w}} \rangle} \mathrm{e}^{-\mathrm{i}\langle \bm{\mathsf{y}}, \, \bm{\mathsf{w}} \rangle} \, \rho(\mathrm{d}{\bm{\mathsf{w}}}).$$
The constant $c_{\phi}$ is a normalizing factor that depends only on
$\phi$, and $\mathrm{i}$ is the imaginary unit.*

*Bochner's theorem immediately delivers a random feature map for the
translation-invariant kernel $k$:
$$\psi( \bm{\mathsf{x}}; \bm{\mathsf{w}} ) = \sqrt{c_{\phi}} \, \mathrm{e}^{\mathrm{i}\langle \bm{\mathsf{x}}, \, \bm{\mathsf{w}} \rangle}
\quad\text{where $\bm{\mathsf{w}} \sim \rho$.}$$ This was one of the
original examples of a random feature map [@RR08:Random-Features]. When
working with data in $\mathbb{R}^d$, the construction can also be
modified to avoid complex values.*

*The key example of a positive definite, translation-invariant kernel is
the Gaussian kernel on $\mathbb{F}^d$, defined in
Example [61](#ex:rbf){reference-type="ref" reference="ex:rbf"}. The
Gaussian kernel is derived from the function
$$\phi( \bm{\mathsf{x}} ) = \mathrm{e}^{-\Vert \bm{\mathsf{x}} \Vert^2 / (2 \sigma^2)}
\quad\text{where the bandwidth $\sigma > 0$.}$$ The associated random
feature map is
$$\psi( \bm{\mathsf{x}}; \bm{\mathsf{w}} ) = \mathrm{e}^{\mathrm{i}\langle \bm{\mathsf{x}}, \, \bm{\mathsf{w}} \rangle}
\quad\text{where}\quad
\bm{\mathsf{w}} \sim \textsc{normal}(\bm{\mathsf{0}}, \sigma^{-2} \bm{\mathsf{I}}) \in \mathbb{F}^d.$$
This fact is both beautiful and useful because of the ubiquity of the
Gaussian kernel in data analysis.*
:::

There are many other kinds of kernels that admit random feature maps.
Random feature maps for dot product kernels were obtained in
[@KK12:Random-Feature; @PP13:Fast-Scalable] and
[@HXGD14:Compact-Random]. For nonstationary kernels, see
[@SR15:Generalized-Spectral] and [@TFSB18:Spatial-Mapping]. Catalogs of
examples appear in and .

### Random features and kernel matrix approximation {#sec:rf-kernel}

We can use the random feature map to construct an empirical
approximation of the kernel matrix $\bm{\mathsf{K}} \in \mathbb{H}_n$
induced by the dataset
$\{ \bm{\mathsf{x}}_1, \dots, \bm{\mathsf{x}}_n \}$. To do so, we draw a
random variable $\bm{\mathsf{w}} \in \mathcal{W}$ with the distribution
$\rho$. Then we form a random vector
$$\bm{\mathsf{z}} = \begin{bmatrix} z_1 \\ \vdots \\ z_n \end{bmatrix}
    = \begin{bmatrix} \psi(\bm{\mathsf{x}}_1; \bm{\mathsf{w}}) \\ \vdots \\ \psi(\bm{\mathsf{x}}_n; \bm{\mathsf{w}}) \end{bmatrix} \in \mathbb{F}^n.$$
Note that we are using the *same* random variable $\bm{\mathsf{w}}$ for
each data point. A realization of the random vector $\bm{\mathsf{z}}$ is
called a *random feature*. The reproducing
property [\[eqn:repro-prop\]](#eqn:repro-prop){reference-type="eqref"
reference="eqn:repro-prop"} ensures that the random feature verifies the
identity
$$(\bm{\mathsf{K}})_{ij} = k(\bm{\mathsf{x}}_i, \bm{\mathsf{x}}_j)
    = \int \psi(\bm{\mathsf{x}}_i;\bm{\mathsf{w}}) \, \psi(\bm{\mathsf{x}}_j;\bm{\mathsf{w}})^* \, \rho(\mathrm{d}{\bm{\mathsf{w}}})
    = \operatorname{\mathbb{E}}[ z_i \cdot z_j^* ].$$ In matrix form,
$$\bm{\mathsf{K}} = \operatorname{\mathbb{E}}[ \bm{\mathsf{zz}}^* ].$$
Therefore, the random rank-one psd matrix
$\bm{\mathsf{Z}} = \bm{\mathsf{zz}}^*$ is an unbiased estimator for the
kernel matrix. The latter display proves that a kernel $k$ must be
positive definite if it admits a random feature map.

To approximate the kernel matrix, we can average $r$ copies of the
rank-one estimator:
$$\bar{\bm{\mathsf{Z}}}_r:= \frac{1}{r} \sum\nolimits_{i=1}^r \bm{\mathsf{Z}}_i
\quad\text{where $\bm{\mathsf{Z}}_i \sim \bm{\mathsf{Z}}$ are iid.}$$ If
the points in the input space $\mathcal{X}$ are parameterized by $d$
numbers, each random feature typically requires $O(nd)$ arithmetic. The
total cost of forming the approximation $\bar{\bm{\mathsf{Z}}}_r$ is
thus $O(rnd)$. When $r \ll n$, we can obtain substantial improvements
over the direct approach of computing the kernel matrix
$\bm{\mathsf{K}}$ explicitly at a cost of $O(n^2 d)$.

### Analysis of the random feature approximation

How many random features are enough to approximate the kernel matrix in
spectral norm? Theorem [13](#thm:mtx-sampling){reference-type="ref"
reference="thm:mtx-sampling"} delivers bounds.

For simplicity, assume that the kernel satisfies
$k(\bm{\mathsf{x}}, \bm{\mathsf{x}}) = 1$ for all
$\bm{\mathsf{x}} \in \mathcal{X}$; the angular similarity kernel and the
Gaussian kernel both enjoy this property. For an accuracy parameter
$\varepsilon> 0$, suppose that we select
$$r \geq 2b \varepsilon^{-2} \operatorname{intdim}(\bm{\mathsf{K}}) \log(2n).$$
The number $b$ is the uniform bound on the feature map $\psi$ defined
in [\[eqn:repro-prop\]](#eqn:repro-prop){reference-type="eqref"
reference="eqn:repro-prop"}, and the intrinsic dimension is defined
in [\[eqn:intdim\]](#eqn:intdim){reference-type="eqref"
reference="eqn:intdim"}.
Theorem [13](#thm:mtx-sampling){reference-type="ref"
reference="thm:mtx-sampling"} implies that the empirical approximation
$\bar{\bm{\mathsf{Z}}}_r$ satisfies
$$\frac{\operatorname{\mathbb{E}}\Vert  \bar{\bm{\mathsf{Z}}}_r - \bm{\mathsf{K}}  \Vert}{\Vert \bm{\mathsf{K}} \Vert}
    \leq \varepsilon+ \varepsilon^2.$$ In other words, we achieve a
relative error approximation of the kernel matrix in spectral norm when
the number $r$ of random features is proportional to the number of
energetic dimensions in the range of the matrix $\bm{\mathsf{K}}$.

This analysis is due to ; see also [@Tro15:Introduction-Matrix
Sec. 6.5]. For learning applications, such as KPCA or KRR, this result
suggests that we need about $O(n\log n)$ random features to obtain
optimal generalization guarantees (where $\varepsilon= n^{-1/2}$). In
fact, roughly $O(\sqrt{n} \, \log n)$ random features are sufficient to
achieve optimal learning rates. This claim depends on involved arguments
from learning theory that are outside the realm of linear algebra. For
example, see
[@SS15:Optimal-Rates; @RR17:Generalization-Properties; @UMMA18:Streaming-Kernel; @SS19:Kernel-Derivative]
and [@Wan19:Simple-Almost].

### Randomized embeddings and random features {#sec:fastfood}

Random feature approximations are faster than explicit computation of a
kernel matrix. Even so, it takes a significant amount of effort to
extract $r$ random features and form an empirical approximation
$\bar{\bm{\mathsf{Z}}}_r$ of the kernel matrix. Several groups have
proposed using structured random embeddings
(Section [9](#sec:dimension-reduction){reference-type="ref"
reference="sec:dimension-reduction"}) to accelerate this process; see
[@PP13:Fast-Scalable; @LSS13:Fastfood-Computing] and
[@HXGD14:Compact-Random].

As an example, let us summarize a heuristic method, called *FFT
Fastfood* [@LSS13:Fastfood-Computing], for speeding up the computation
of random features for the complex Gaussian kernel with bandwidth
$\sigma^2$. Consider the matrix $\bm{\mathsf{X}}$ formed from $n$
observations in $\mathbb{C}^d$:
$$\bm{\mathsf{X}} = \begin{bmatrix} \bm{\mathsf{x}}_1^* \\ \vdots \\ \bm{\mathsf{x}}_n^* \end{bmatrix} \in \mathbb{C}^{n \times d}.$$
Let $\bm{\mathsf{\Gamma}} \in \mathbb{C}^{d \times d}$ be a matrix with
iid complex $\textsc{normal}(0, \sigma^{-2})$ entries. Then we can
simultaneously compute $d$ random features
$\bm{\mathsf{z}}_1, \dots, \bm{\mathsf{z}}_d$ for the Gaussian kernel by
forming a matrix product and applying the exponential map:
$$\exp\cdot\,(\mathrm{i}\bm{\mathsf{X}} \bm{\mathsf{\Gamma}}) = \begin{bmatrix} \bm{\mathsf{z}}_1 & \dots & \bm{\mathsf{z}}_d \end{bmatrix} \in \mathbb{C}^{n \times d}.$$
We have written $\exp\cdot$ for the *entrywise* exponential. This
procedure typically involves $O(nd^2)$ arithmetic.

The idea behind FFT Fastfood is to accelerate this computation by
replacing the Gaussian matrix with a structured random matrix. This
exchange is motivated by the observed universality properties of random
embeddings. Consider a random matrix of the form
$$\bm{\mathsf{S}} = \frac{1}{\sigma} \, \bm{\mathsf{E}} \bm{\mathsf{\Pi}} \bm{\mathsf{F}} \in \mathbb{C}^{d \times d},$$
where $\bm{\mathsf{E}}$ is a random sign flip, $\bm{\mathsf{\Pi}}$ is a
random permutation, and $\bm{\mathsf{F}}$ is the discrete DFT. The FFT
algorithm supports efficient matrix products with $\bm{\mathsf{S}}$.
Therefore, we can simultaneously extract $d$ random features by
computing
$\exp\cdot\,(\mathrm{i}\bm{\mathsf{XS}}) \in \mathbb{C}^{n \times d}$.
This procedure uses only $O(nd \log d)$ operations. To form $r$ random
features where $r > d$, we simply repeat the same process
$\lceil r/d \rceil$ times.

Compared with using a Gaussian matrix product, FFT Fastfood gives a
substantial reduction in arithmetic. Even so, the performance for
learning is almost identical to a direct application of the random
feature approximation.

### Random features and streaming matrix approximation {#sec:streaming-kpca}

Suppose that we wish to perform KPCA. The direct random features
approach requires us to form the empirical approximation
$\bar{\bm{\mathsf{Z}}}_r$ of the kernel matrix $\bm{\mathsf{K}}$ and to
compute its rank-$\ell$ truncated eigenvalue decomposition. It is often
the case that the desired number $\ell$ of principal components is far
smaller than the number $r$ of random features we need to obtain a
suitable approximation of the kernel matrix. In this case, we can
combine random features with streaming matrix approximation to make
economies in storage and computation.

Let
$\{ \bm{\mathsf{z}}_1, \bm{\mathsf{z}}_2, \bm{\mathsf{z}}_3, \dots \} \subset \mathbb{F}^n$
be an iid sequence of random features for the kernel matrix
$\bm{\mathsf{K}} \in \mathbb{H}_n$. The empirical approximation
$\bar{\bm{\mathsf{Z}}}_r$ of the kernel matrix, obtained from the first
$r$ random features, follows the recursion
$$\bar{\bm{\mathsf{Z}}}_0 = \bm{\mathsf{0}}
\quad\text{and}\quad
\bar{\bm{\mathsf{Z}}}_t = (1 - t^{-1}) \bar{\bm{\mathsf{Z}}}_{t-1} + t^{-1} \, \bm{\mathsf{z}}_t \bm{\mathsf{z}}_t^*
\quad\text{for $t = 1, 2, 3, \dots$.}$$ This is a psd matrix, generated
by a stream of linear updates. Therefore, we can track the evolution
using a streaming Nyström approximation
(Section [14](#sec:nystrom){reference-type="ref"
reference="sec:nystrom"}).

Let $\bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times s}$ be a random test
matrix, with $s > \ell$. By performing rank-one updates, we can
efficiently maintain the sample matrices
$$\bm{\mathsf{Y}}_t = \bar{\bm{\mathsf{Z}}}_t \bm{\mathsf{\Omega}} \in \mathbb{F}^{n \times s}
\quad\text{for $t = 1, 2, 3, \dots$.}$$ After collecting a sufficient
number $r$ of samples, we can apply
Algorithm [\[alg:random-nystrom\]](#alg:random-nystrom){reference-type="ref"
reference="alg:random-nystrom"} to $\bm{\mathsf{Y}}_t$ to obtain a
near-optimal rank-$\ell$ eigenvalue decomposition of the empirical
approximation $\bar{\bm{\mathsf{Z}}}_r$.

It usually suffices to take the sketch size $s$ to be proportional to
the rank $\ell$ of the truncated eigenvalue decomposition. In this case,
the overall approach uses $O(\ell n)$ storage. We can generate and
process $r$ random features using $O((d+\ell)rn)$ arithmetic, where $d$
is the dimension of $\mathcal{X}$. The subsequent cost of the Nyström
approximation is $O(\ell^2 n)$ operations. The streaming random features
approach has storage and arithmetic costs roughly $\ell / r$ times those
of the direct random features approach. The streaming method can be
combined with dimension reduction techniques
(Section [19.3.4](#sec:fastfood){reference-type="ref"
reference="sec:fastfood"}) for further acceleration.

::: remark
**Remark 67** (History). *proposed using a stream of random features to
perform KPCA; their algorithm tracks the stream with the (deterministic)
frequent directions sketch [@GLPW16:Frequent-Directions]. We have
presented a new variant, based on randomized Nyström approximation, that
is motivated by the work in . have develop a somewhat different
streaming KPCA algorithm based on Oja's
method [@Oja82:Simplified-Neuron]. At present, we lack a full empirical
comparison of these alternatives.*
:::

# High-accuracy approximation of kernel matrices {#sec:rankstructured}

In this section, we continue the discussion of kernel matrices that we
started in Section [19](#sec:kernel){reference-type="ref"
reference="sec:kernel"}, but we now consider the high-accuracy regime.
In particular, given a kernel matrix $\bm{\mathsf{K}}$, we seek an
approximation $\bm{\mathsf{K}}_{\mathrm{approx}}$ for which
$\|\bm{\mathsf{K}} - \bm{\mathsf{K}}_{\mathrm{approx}}\|$ is small, say
of relative accuracy $10^{-3}$ or $10^{-6}$. This objective was not
realistic for the applications discussed in Section
[19](#sec:kernel){reference-type="ref" reference="sec:kernel"}, but it
can be achieved in situations where we have access to fast techniques
for evaluating the matrix-vector product
$\bm{\mathsf{x}} \mapsto \bm{\mathsf{K}}\bm{\mathsf{x}}$ (and also
$\bm{\mathsf{x}} \mapsto \bm{\mathsf{K}}^{*}\bm{\mathsf{x}}$ when
$\bm{\mathsf{K}}$ is not self-adjoint). The algorithms that we describe
will build a data sparse approximation to $\bm{\mathsf{K}}$ by using
information in samples such as $\bm{\mathsf{K}}\bm{\mathsf{x}}$ and
$\bm{\mathsf{K}}^{*}\bm{\mathsf{x}}$ for random vectors
$\bm{\mathsf{x}}$. These techniques are particularly well suited to
problems that arise in modeling physical phenomena such as
electromagnetic scattering, or the deformation of solid bodies; we
describe how the fast matrix-vector application we need can be realized
in Section [20.3](#sec:cpkernel_application){reference-type="ref"
reference="sec:cpkernel_application"}.

As in Section [19](#sec:kernel){reference-type="ref"
reference="sec:kernel"}, we say that a matrix
$\bm{\mathsf{K}} \in \mathbb{C}^{n\times n}$ is a *kernel matrix* if its
entries are given by a formula such as $$\begin{equation}
\label{eq:Aij_kernel}
\bm{\mathsf{K}}(i,j) = k(\bm{\mathsf{x}}_{i},\bm{\mathsf{x}}_{j}),
\end{equation}$$ where $\{\bm{\mathsf{x}}_{i}\}_{i=1}^{n}$ is a set of
points in $\mathbb{R}^{d}$, and where
$k\,\colon\,\mathbb{R}^{d} \times \mathbb{R}^{d} \rightarrow \mathbb{C}$
is a kernel function. The kernel matrices that we consider tend to have
singular values that decay slowly or not at all, which rules out the
possibility that $\bm{\mathsf{K}}_{\mathrm{approx}}$ could have low
rank. Instead, we build an approximation
$\bm{\mathsf{K}}_{\mathrm{approx}}$ that is tessellated into $O(n)$
blocks in such as way that each off-diagonal block has low rank. Figure
[11](#fig:tree){reference-type="ref" reference="fig:tree"}(b) shows a
representative tessellation pattern. We say that a matrix of this type
is a *rank-structured hierarchical matrix.*

The purpose of determining a rank-structured approximation to an
operator for which we already have fast matrix-vector multiplication
techniques available is that the new representation can be used to
rapidly execute a whole range of linear algebraic operations: matrix
inversion, LU factorization, and even full spectral decompositions in
certain cases.

This section is structured to provide a high-level description of the
core ideas in Sections [20.1](#sec:sepvar){reference-type="ref"
reference="sec:sepvar"} --
[20.3](#sec:cpkernel_application){reference-type="ref"
reference="sec:cpkernel_application"}. Additional details follow in
Sections [20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"}
-- [20.9](#sec:ASKIT){reference-type="ref" reference="sec:ASKIT"}.

## Separation of variables and low-rank approximation {#sec:sepvar}

The reason that many kernel matrices can be tessellated into blocks that
have low numerical rank is that the function
$(\bm{\mathsf{x}},\bm{\mathsf{y}}) \mapsto k(\bm{\mathsf{x}},\bm{\mathsf{y}})$
is typically smooth as long as $\bm{\mathsf{x}}$ and $\bm{\mathsf{y}}$
are not close. To illustrate the connection, let us consider a
computational domain $D$ that holds a set of points
$\{\bm{\mathsf{x}}_{i}\}_{i=1}^{n} \subset \mathbb{R}^{2}$, as shown in
Figure [10](#fig:sources_targets){reference-type="ref"
reference="fig:sources_targets"}. Suppose further that $D_{\mathrm{s}}$
and $D_{\mathrm{t}}$ are two subdomains of $D$ that are located a bit
apart from each other, as shown in the figure. When the kernel function
$k$ is smooth, we can typically approximate it to high accuracy through
an approximate separation of variables of the form $$\begin{equation}
\label{eq:kernel_sepvar}
k(\bm{\mathsf{x}},\bm{\mathsf{y}}) \approx \sum_{p=1}^{P}b_{p}(\bm{\mathsf{x}})\,c_{p}(\bm{\mathsf{y}}),\qquad\bm{\mathsf{x}} \in D_{\mathrm{t}},\ \bm{\mathsf{y}} \in D_{\mathrm{s}}.
\end{equation}$$ Let $I_{\mathrm{s}}$ and $I_{\mathrm{t}}$ denote two
index vectors that identify the points located in $D_{\mathrm{s}}$ and
$D_{\mathrm{t}}$, respectively (so that, for example,
$i\in I_{\mathrm{s}}$ if and only if
$\bm{\mathsf{x}}_{i} \in D_{\mathrm{s}}$). Then combining the formula
([\[eq:Aij_kernel\]](#eq:Aij_kernel){reference-type="ref"
reference="eq:Aij_kernel"}) with the separation of variables
([\[eq:kernel_sepvar\]](#eq:kernel_sepvar){reference-type="ref"
reference="eq:kernel_sepvar"}), we get $$\begin{equation}
\label{eq:A_sepvar}
\bm{\mathsf{K}}(i,j) \approx \sum_{p=1}^{P}b_{p}(\bm{\mathsf{x}}_{i})\,c_{p}(\bm{\mathsf{x}}_{j}),\qquad i \in I_{\mathrm{t}},\ j \in I_{\mathrm{s}}.
\end{equation}$$ Equation
([\[eq:A_sepvar\]](#eq:A_sepvar){reference-type="ref"
reference="eq:A_sepvar"}) is exactly the low-rank approximation to the
block $\bm{\mathsf{K}}(I_{\mathrm{t}},I_{\mathrm{s}})$ that we seek. To
be precise, ([\[eq:A_sepvar\]](#eq:A_sepvar){reference-type="ref"
reference="eq:A_sepvar"}) can be written as
$$\bm{\mathsf{K}}(I_{\mathrm{t}},I_{\mathrm{s}}) \approx \bm{\mathsf{B}}\bm{\mathsf{C}},$$
where $\bm{\mathsf{B}}$ and $\bm{\mathsf{C}}$ are defined via
$\bm{\mathsf{B}}(i,p) = b_{p}(\bm{\mathsf{x}}_{i})$ and
$\bm{\mathsf{C}}(p,j) = c_{p}(\bm{\mathsf{x}}_{j})$.

A separation of variables such as
([\[eq:kernel_sepvar\]](#eq:kernel_sepvar){reference-type="ref"
reference="eq:kernel_sepvar"}) is sometimes provided through analytic
knowledge of the kernel function, as illustrated in Example
[68](#example:multipole){reference-type="ref"
reference="example:multipole"}. Perhaps more typically, all we know is
that such a formula *should* in principle exist, for instance because we
know that the matrix approximates a singular integral operator for which
a Calderón-Zygmund decomposition must exist. It is then the task of the
randomized algorithm to explicitly build the factors $\bm{\mathsf{B}}$
and $\bm{\mathsf{C}}$, given a computational tolerance.

<figure id="fig:sources_targets">
<div class="picture">
<p>(123,67) (05,00)<span><embed src="Pics/fig_sources_targets.pdf"
style="width:65mm" /></span> (00,03)<span><span
class="math inline"><em>D</em></span></span> (26,15)<span> <span
class="math inline"><strong>D</strong><sub>s</sub></span></span>
(57,34)<span><span
class="math inline"><strong>D</strong><sub>t</sub></span></span>
(72,34)</p>
<div class="minipage">
<p>A box <span class="math inline"><em>D</em></span> holding points
<span
class="math inline">{x<sub><em>i</em></sub>}<sub><em>i</em> = 1</sub><sup><em>n</em></sup></span>
(the blue, red, and gray dots) that define a kernel matrix <span
class="math inline">K(<em>i</em>, <em>j</em>) = <em>k</em>(x<sub><em>i</em></sub>, x<sub><em>j</em></sub>)</span>.
The regions <span class="math inline"><em>D</em><sub>s</sub></span> (the
red box) and <span class="math inline"><em>D</em><sub>t</sub></span>
(the blue box) are separated enough that <span
class="math inline"><em>k</em>(x, y)</span> is smooth when <span
class="math inline">x ∈ <em>D</em><sub>t</sub></span> and <span
class="math inline">y ∈ <em>D</em><sub>s</sub></span>. In consequence,
<span
class="math inline">K(<em>I</em><sub>t</sub>, <em>I</em><sub>s</sub>)</span>
has low numerical rank, where <span
class="math inline"><em>I</em><sub>s</sub></span> identifies the red
points and <span class="math inline"><em>I</em><sub>t</sub></span> the
blue.</p>
</div>
</div>
<figcaption>The geometry discussed in Section <a
href="#eq:kernel_sepvar" data-reference-type="ref"
data-reference="eq:kernel_sepvar">[eq:kernel_sepvar]</a>.</figcaption>
</figure>

::: {#example:multipole .example}
**Example 68** (Laplace kernel in two dimensions). *A standard example
of a kernel matrix in mathematical physics is the matrix
$\bm{\mathsf{K}}$ that maps a vector of electric source strengths
$\bm{\mathsf{q}} = (q_{j})_{j=1}^{n}$ to a vector of potentials
$\bm{\mathsf{u}} = (u_{i})_{i=1}^{n}$. When the set
$\{\bm{\mathsf{x}}_{i}\}_{i=1}^{n} \subset \mathbb{R}^{2}$ identifies
the locations of both the source and the target points,
$\bm{\mathsf{K}}$ takes the form
([\[eq:Aij_kernel\]](#eq:Aij_kernel){reference-type="ref"
reference="eq:Aij_kernel"}), for $$k(\bm{\mathsf{x}},\bm{\mathsf{y}}) =
\left\{\begin{aligned}
\log|\bm{\mathsf{x}}-\bm{\mathsf{y}}|,&\quad\mbox{when}\ \bm{\mathsf{x}}\neq \bm{\mathsf{y}},\\
0,&\quad\mbox{when}\ \bm{\mathsf{x}} = \bm{\mathsf{y}}.
\end{aligned}\right.$$ There is a well-known result from potential
theory that provides the required separation of variables. Expressing
$\bm{\mathsf{x}}$ and $\bm{\mathsf{y}}$ in polar coordinates with
respect to some expansion center $\bm{\mathsf{c}}$, so that
$\bm{\mathsf{x}} - \bm{\mathsf{c}} = r(\cos\theta,\,\sin\theta)$ and
$\bm{\mathsf{y}} - \bm{\mathsf{c}} = r'(\cos\theta',\,\sin\theta')$, the
separation of variables (known as a *multipole expansion*)
$$\begin{equation}
\label{eq:mpole}
k(\bm{\mathsf{x}},\bm{\mathsf{y}}) = \log(r) - \sum_{p=1}^{\infty}\frac{1}{p}\left(\frac{r'}{r}\right)^{p}\bigl(\cos(p\theta)\cos(p\theta') + \sin(p\theta)\sin(p\theta')\bigr)
\end{equation}$$ is valid whenever $r' < r$.*
:::

## Rank-structured matrices and randomized compression {#sec:bestiary}

The observation that the off-diagonal blocks of a kernel matrix often
have low numerical rank underpins many "fast" algorithms in
computational physics. In particular, it is the foundation of the
Barnes-Hut [@1986_barnes_hut] and fast multipole methods [@rokhlin1987]
for evaluating all $O(n^{2})$ pairwise interactions between $n$
electrically charged particles in linear or close-to-linear complexity.
These methods were generalized by Hackbusch and co-workers, who
developed the $\mathcal{H}$- and $\mathcal{H}^{2}$-matrix frameworks
[@hackbusch; @2003_hackbusch; @2002_hackbusch_H2]. These explicitly
linear algebraic formulations enable fast algorithms not only for
matrix-vector multiplication but also for matrix inversion, LU
factorization, matrix-matrix multiplication and many more.

A fundamental challenge that arises when rank-structured matrix
formulations are used is how to find the data sparse representation of
the operator in the first place. The straightforward approach would be
to form the full matrix, and then loop over all the compressible off
diagonal blocks and compress them using, for example, a singular value
decomposition. The cost of such a process is necessarily at least
$O(n^{2})$, which is rarely affordable. Using randomized compression
techniques, it turns out to be possible to compress all the off-diagonal
blocks *jointly*, without any need for sampling each block individually.
In this section, we describe two such methods. Both require the user to
supply a fast algorithm for applying the full matrix (and its adjoint)
to vectors; see Section
[20.3](#sec:cpkernel_application){reference-type="ref"
reference="sec:cpkernel_application"}. The two methods have different
computational profiles.

- The technique described in Sections
  [20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"} and
  [20.5](#sec:HODLRcompression){reference-type="ref"
  reference="sec:HODLRcompression"} is a true "black-box" technique that
  interacts with $\bm{\mathsf{K}}$ only through the matrix-vector
  multiplication. It has storage complexity $O(n\log n)$ and it requires
  $O(\log n)$ applications of $\bm{\mathsf{K}}$ to a random matrix of
  size $n\times (r+p)$ where $r$ is an upper bound on the ranks of the
  off-diagonal blocks, and $p$ is a small over-sampling parameter.

- The technique described in Sections
  [20.6](#sec:HBS){reference-type="ref" reference="sec:HBS"} and
  [20.7](#sec:HBScompression){reference-type="ref"
  reference="sec:HBScompression"} attains true linear $O(n)$ complexity,
  and requires only a single application of $\bm{\mathsf{K}}$ and
  $\bm{\mathsf{K}}^{*}$ to a random matrix of size $n\times (r+p)$,
  where $r$ and $p$ are as in (a). Its drawbacks are that it requires
  evaluation of $O(rn)$ individual matrix entries, and that it works
  only for a smaller class of matrices.

::: remark
**Remark 69** (Scope of Section
[20](#sec:rankstructured){reference-type="ref"
reference="sec:rankstructured"}). *To keep the presentation as
uncluttered by burdensome notation as possible, in this survey we
restrict attention to two basic "formats" for representing a
rank-structured hierarchical matrix. In Sections
[20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"} and
[20.5](#sec:HODLRcompression){reference-type="ref"
reference="sec:HODLRcompression"}, we use the *hierarchically
off-diagonal low rank (HODLR)* format and in Sections
[20.6](#sec:HBS){reference-type="ref" reference="sec:HBS"} and
[20.7](#sec:HBScompression){reference-type="ref"
reference="sec:HBScompression"} we use the *hierarchically block
separable (HBS)* format (sometimes referred to as *hierarchically semi
separable (HSS)* matrices). The main limitation of these formats is that
they require all off-diagonal blocks of the matrix to have low numerical
rank. This is realistic only when the points
$\{\bm{\mathsf{x}}_{i}\}_{i=1}^{n}$ are restricted to a low-dimensional
manifold. In practical applications, one sometimes has to leave a larger
part of the matrix uncompressed, to avoid attempting to impose a
separation of variables such as
([\[eq:kernel_sepvar\]](#eq:kernel_sepvar){reference-type="ref"
reference="eq:kernel_sepvar"}) on a kernel
$k = k(\bm{\mathsf{x}},\bm{\mathsf{y}})$ when $\bm{\mathsf{x}}$ and
$\bm{\mathsf{y}}$ are too close. This is done through enforcing what is
called a "strong admissibility condition" (in contrast to the "weak
admissibility condition" of the HODLR and HBS formats), as was done in
the original Barnes-Hut and fast multipole methods.*
:::

::: {#remark:compressioncompetition .remark}
**Remark 70** (Alternative compression strategies). *Let us briefly
describe what alternatives to randomized compression exist. The original
papers on $\mathcal{H}$-matrices used Taylor approximations to derive a
separation of variables, but this works only when the kernel is given
analytically. It also tends to be quite expensive. The *adaptive
cross-approximation (ACA)* technique of [@2002_rjasanow_ACA] and
[@2006_bebendorf_ACA] relies on using "natural basis" vectors (see
Section [13.1](#sec:CURandIDdef){reference-type="ref"
reference="sec:CURandIDdef"}) that are found using semi-heuristic
techniques. The method can work very well in practice, but is not
guaranteed to provide an accurate factorization. When the kernel matrix
comes from mathematical physics, specialized techniques that exploit
mathematical properties of the kernel function often perform well
[@2005_martinsson_fastdirect]. For a detailed discussion of compression
of rank-structured matrices, we refer to [@2019_book Ch. 17].*
:::

## Computational environments {#sec:cpkernel_application}

The compression techniques that we describe rely on the user providing a
fast algorithm for applying the operator to be compressed to vectors.
Representative environments where such fast algorithms are available
include the following.

*Matrix-matrix multiplication:* Suppose that
$\bm{\mathsf{K}} = \bm{\mathsf{B}}\bm{\mathsf{C}}$, where
$\bm{\mathsf{B}}$ and $\bm{\mathsf{C}}$ are matrices that can rapidly be
applied to vectors. Then the randomized compression techniques allow us
to compute their product.

*Compression of boundary integral operators:* It is often possible to
reformulate a boundary value problem involving an elliptic partial
different operator, such as the Laplace or the Helmholtz operators, as
an equivalent boundary integral equation (BIE), see [@2019_book Part
III]. When such a BIE is discretized, the result is a kernel matrix
which can be applied rapidly to vectors using fast summation techniques
such as the fast multipole method [@rokhlin1987]. Randomized compression
techniques allow us to build an approximation to the matrix that can be
factorized or inverted, thus enabling direct (as opposed to iterative)
solvers.

*Dirichlet-to-Neumann (DtN) operators:* A singular integral operator of
central importance in engineering, physics, and scientific computing is
the DtN operator which maps given Dirichlet data for an elliptic
boundary value problem to the boundary fluxes of the corresponding
solution. The kernel of the DtN operator is rarely known analytically,
but the operator can be applied through a fast solver for the PDE, such
as, e.g., a finite-element discretization combined with a multigrid
solver. The randomized techniques described allow for the DtN operator
to be built and stored explicitly.

*Frontal matrices in sparse direct solvers:* Suppose that
$\bm{\mathsf{S}}$ is a large sparse matrix arising from the
discretization of an elliptic PDEs. A common technique for solving
$\bm{\mathsf{S}}\bm{\mathsf{x}} = \bm{\mathsf{b}}$ is to compute an LU
factorization of the matrix $\bm{\mathsf{S}}$. There are ways to do this
that preserve sparsity as far as possible, but in the course of the
factorization procedure, certain dense matrices of increasing size will
need to be factorized. These matrices turn out to be kernel matrices
with kernels that are not known explicitly, but that can be built using
the techniques described here
[@2013_xia_randomized; @2017_ghysels_robust].

## Hierarchically off-diagonal low-rank matrices {#sec:HODLR}

In order to illustrate how randomized methods can be used to construct
data sparse representations of rank-structured matrices, we will
describe a particularly simple 'format' in this section that is often
referred to as the hierarchically off-diagonal low-rank (HODLR) format.
This is a basic format that works well when the points are organized on
a one- or two-dimensional manifold.

The first step towards defining the HODLR format is to build a binary
tree on the index vector $I = [1,2,3,\dots,n]$ through a process that is
illustrated in Figure [11](#fig:tree){reference-type="ref"
reference="fig:tree"}(a). With any node $\tau$ in the tree, we associate
an index vector $I_{\tau} \subseteq I$. The root of the tree is given
the index $\tau=1$, and we associate it with the full index vector, so
that $I_{1} = I$. At the next finer level of the tree, we split $I_{1}$
into two parts $I_{2}$ and $I_{3}$ so that $I_{1} = I_{2} \cup I_{3}$
forms a disjoint partition. Then continue splitting the index vectors
until each remaining index vector is "short". (Exactly what "short"
means is application-dependent, but one may think of a short vector as
holding a few hundred indices or so.) We let $\ell$ denote a *level* of
the tree, with $\ell=0$ denoting the root, so that level $\ell$ holds
$2^{\ell}$ nodes. We use the terms *parent* and *child* in the obvious
way, and say that a pair of nodes $\{\alpha,\beta\}$ forms a *sibling
pair* if they have the same parent. A *leaf* node is of course a node
that has no children.

The binary tree that we defined induces a natural tessellation of the
kernel matrix $\bm{\mathsf{K}}$ into $O(n)$ blocks. Figure
[11](#fig:tree){reference-type="ref" reference="fig:tree"}(b) shows the
tessellation that follows from the tree in Figure
[11](#fig:tree){reference-type="ref" reference="fig:tree"}(a). Each
parent node $\tau$ in the tree gives rise to two off-diagonal blocks
that both have low numerical rank. Letting $\{\alpha,\beta\}$ denote the
children of $\tau$, these two blocks are $$\begin{equation}
\label{eq:siblingboxes}
\bm{\mathsf{K}}_{\alpha,\beta} = \bm{\mathsf{K}}(I_{\alpha},I_{\beta})
\qquad\mbox{and}\qquad
\bm{\mathsf{K}}_{\beta,\alpha} = \bm{\mathsf{K}}(I_{\beta},I_{\alpha}).
\end{equation}$$ For each leaf node $\tau$, we define a corresponding
diagonal block as $$\begin{equation}
\label{eq:diagbox}
\bm{\mathsf{D}}_{\tau} = \bm{\mathsf{K}}(I_{\tau},I_{\tau}).
\end{equation}$$ The disjoint partition of $\bm{\mathsf{K}}$ into blocks
is now formed by all sibling pairs, as defined in
([\[eq:siblingboxes\]](#eq:siblingboxes){reference-type="ref"
reference="eq:siblingboxes"}), together with all diagonal blocks, as
defined by ([\[eq:diagbox\]](#eq:diagbox){reference-type="ref"
reference="eq:diagbox"}). When each off-diagonal block in this
tessellation has low rank, we say that $\bm{\mathsf{K}}$ is a
*hierarchically off-diagonal low-rank (HODLR)* matrix.

<figure id="fig:tree">
<div class="picture">
<p>(120,120) (00,115)<span>(a)</span> (25,85)</p>
<div class="picture">
<p>(122,44) ( 20, 0)<span><embed src="Pics/fig_tree.pdf"
style="width:75mm" /></span> ( 0,30)<span><em>Level <span
class="math inline">0</span>:</em></span> ( 0,21)<span><em>Level <span
class="math inline">1</span>:</em></span> ( 0,12)<span><em>Level <span
class="math inline">2</span>:</em></span> ( 0,03)<span><em>Level <span
class="math inline">3</span>:</em></span></p>
</div>
<p>(00,77)<span>(b)</span> (40,02)</p>
<div class="picture">
<p>(80,80) (-2,-2)<span><img src="Pics/fig_8x8_Smatrix.png"
style="width:82mm" alt="image" /></span> (00,73)<span><span
class="math inline">D<sub>8</sub></span></span> (10,63)<span><span
class="math inline">D<sub>9</sub></span></span> (20,53)<span><span
class="math inline">D<sub>10</sub></span></span> (30,43)<span><span
class="math inline">D<sub>11</sub></span></span> (40,33)<span><span
class="math inline">D<sub>12</sub></span></span> (50,23)<span><span
class="math inline">D<sub>13</sub></span></span> (60,13)<span><span
class="math inline">D<sub>14</sub></span></span> (70,03)<span><span
class="math inline">D<sub>15</sub></span></span> (00,63)<span><span
class="math inline">K<sub>9, 8</sub></span></span> (10,73)<span><span
class="math inline">K<sub>8, 9</sub></span></span> (60,03)<span><span
class="math inline">K<sub>15, 14</sub></span></span> (70,13)<span><span
class="math inline">K<sub>14, 15</sub></span></span> (40,23)<span><span
class="math inline">K<sub>13, 12</sub></span></span> (50,33)<span><span
class="math inline">K<sub>12, 13</sub></span></span> (20,43)<span><span
class="math inline">K<sub>11, 10</sub></span></span> (30,53)<span><span
class="math inline">K<sub>10, 11</sub></span></span> (16,18)<span><span
class="math inline">K<sub>3, 2</sub></span></span> (56,58)<span><span
class="math inline">K<sub>2, 3</sub></span></span> (06,48)<span><span
class="math inline">K<sub>5, 4</sub></span></span> (26,68)<span><span
class="math inline">K<sub>4, 5</sub></span></span> (46,08)<span><span
class="math inline">K<sub>7, 6</sub></span></span> (66,28)<span><span
class="math inline">K<sub>6, 7</sub></span></span> (-10,38)<span><span
class="math inline">K=</span></span></p>
</div>
</div>
<figcaption>(a) A binary tree with three levels; see Section <a
href="#sec:HODLR" data-reference-type="ref"
data-reference="sec:HODLR">20.4</a>. Each node <span
class="math inline"><em>τ</em></span> in the tree owns an index vector
<span class="math inline"><em>I</em><sub><em>τ</em></sub></span> that is
a subset of the full index vector <span
class="math inline"><em>I</em> = 1 : <em>n</em></span>. For the root
node, we set <span
class="math inline"><em>I</em><sub>1</sub> = <em>I</em></span>. Each
split in the tree represents a disjoint partition of the corresponding
index vectors, so that, for example, <span
class="math inline"><em>I</em><sub>1</sub> = <em>I</em><sub>2</sub> ∪ <em>I</em><sub>3</sub></span>
and <span
class="math inline"><em>I</em><sub>7</sub> = <em>I</em><sub>14</sub> ∪ <em>I</em><sub>15</sub></span>.
(b) A matrix <span class="math inline">K</span> tessellated according to
the tree shown in (a). </figcaption>
</figure>

## Compressing a rank-structured hierarchical matrix through the matrix-vector multiplication only {#sec:HODLRcompression}

In this section, we describe a randomized technique for computing a data
sparse representation of a matrix that is compressible in the "HODLR"
format that we introduced in Section
[20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"}. This
technique interacts with $\bm{\mathsf{K}}$ only through the application
of $\bm{\mathsf{K}}$ and $\bm{\mathsf{K}}^{*}$ to vectors, and is in
this sense a "black-box" technique. In particular, we do not need the
ability to evaluate individual entries of $\bm{\mathsf{K}}$. The
following theorem summarizes the main result:

::: theorem
**Theorem 71**. *Let $\bm{\mathsf{K}}$ be a HODLR matrix associated with
a fully populated binary tree on the index vector, as described in
Section [20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"}.
Suppose that the tree has $L$ levels, that each off-diagonal block has
rank at most $k$, and that each leaf node in the tree holds at most $ck$
indices for some fixed number $c$. Then the diagonal blocks
$\bm{\mathsf{D}}_{\tau}$, as well as rank-$k$ factorizations of all
sibling interaction matrices $\bm{\mathsf{K}}_{\alpha,\beta}$ can be
computed by a randomized algorithm with cost at most
$$T_{\mathrm{total}} = T_{\mathrm{matvec}} \times (4L+c)k + T_{\mathrm{flop}} \times O(L^{2}k^{2}n),$$
where $T_{\mathrm{matvec}}$ is the cost of applying either
$\bm{\mathsf{K}}$ or $\bm{\mathsf{K}}^{*}$ to a vector, and
$T_{\mathrm{flop}}$ is the cost of a floating point operation.*
:::

The proof of the theorem consists of an explicit algorithm for building
all matrices that is often referred to as the "peeling algorithm". It
was originally published in [@2011_lin_lu_ying], with later
modifications proposed in [@2016_martinsson_hudson2]. The proof below is
based on [@2019_book Sec. 17.4]. We give the proof for the case
described in the theorem involving a matrix whose off-diagonal blocks
have exact rank $k$. In practical applications, it is of course more
typical to have the off-diagonal blocks be only of approximate low rank.
In this case, the exact same algorithm can be used, but the number of
samples drawn should be increased from $k$ to $k+p$ for some modest
oversampling parameter $p$.

::: proof
*Proof.* The algorithm is a "top-down" technique that compresses the
largest blocks first, and then moves on to process one level at a time
of successively smaller blocks. In describing the technique, we assume
that the blocks are numbered as shown in Figure
[11](#fig:tree){reference-type="ref" reference="fig:tree"}(a).

In the first step of the algorithm, we build approximations to the two
largest blocks $\bm{\mathsf{K}}_{2,3}$ and $\bm{\mathsf{K}}_{3,2}$,
shown in red in Figure [12](#fig:peeling){reference-type="ref"
reference="fig:peeling"}(a). To do this, we form a random matrix
$\bm{\mathsf{\Omega}}$ of size $n\times 2k$ and then use the
matrix-vector multiplication to form a sample matrix
$$\begin{array}{cccccccccccc}
\bm{\mathsf{Y}} &=& \bm{\mathsf{K}}&\bm{\mathsf{\Omega}}.\\
n \times 2k && n \times n & n\times 2k
\end{array}$$ The key idea of the peeling algorithm is to insert zero
blocks in the random matrix $\bm{\mathsf{\Omega}}$ in the pattern shown
in Figure [12](#fig:peeling){reference-type="ref"
reference="fig:peeling"}(a). The zero blocks permit us to extract "pure"
samples from the column spaces of the blocks $\bm{\mathsf{K}}_{2,3}$ and
$\bm{\mathsf{K}}_{3,2}$. For instance, the blocks labeled
$\bm{\mathsf{Y}}_{2}$ and $\bm{\mathsf{Y}}_{3}$ in the figure are given
by the formulas $$\begin{array}{cccccccccccc}
\bm{\mathsf{Y}}_{2} &=& \bm{\mathsf{K}}_{2,3}&\bm{\mathsf{\Omega}}_{3}\\
\tfrac{n}{2} \times k && \tfrac{n}{2} \times \tfrac{n}{2} & \tfrac{n}{2}\times k
\end{array}
\qquad\mbox{and}\qquad
\begin{array}{cccccccccccc}
\bm{\mathsf{Y}}_{3} &=& \bm{\mathsf{K}}_{3,2}&\bm{\mathsf{\Omega}}_{2}\\
\tfrac{n}{2} \times k && \tfrac{n}{2} \times \tfrac{n}{2} & \tfrac{n}{2}\times k
\end{array}$$ By orthonormalizing the matrices $\bm{\mathsf{Y}}_{2}$ and
$\bm{\mathsf{Y}}_{3}$, we obtain ON bases $\bm{\mathsf{U}}_{2}$ and
$\bm{\mathsf{U}}_{3}$ for the off-diagonal blocks
$\bm{\mathsf{K}}_{2,3}$ and $\bm{\mathsf{K}}_{3,2}$. In order to
complete the factorization of $\bm{\mathsf{K}}_{2,3}$ and
$\bm{\mathsf{K}}_{3,2}$, we will perform an operation that is the
equivalent of "Stage B" in Section
[11.2](#sec:rsvd){reference-type="ref" reference="sec:rsvd"}. To this
end, we form a test matrix by interlacing the matrices
$\bm{\mathsf{U}}_{2}$ and $\bm{\mathsf{U}}_{3}$ with zero blocks, to
form the $n\times 2k$ matrix shown in Figure
[13](#fig:peeling_trans){reference-type="ref"
reference="fig:peeling_trans"}. Applying $\bm{\mathsf{K}}^{*}$ to this
test matrix, we get the sample matrices
$$\bm{\mathsf{Z}}_{2} = \bm{\mathsf{K}}_{3,2}^{*}\bm{\mathsf{U}}_{3},
\qquad\mbox{and}\qquad
\bm{\mathsf{Z}}_{3} = \bm{\mathsf{K}}_{2,3}^{*}\bm{\mathsf{U}}_{2}.$$
Since $\bm{\mathsf{U}}_{2}$ holds an ON-basis for the column space of
$\bm{\mathsf{K}}_{2,3}$, it follows that
$$\bm{\mathsf{K}}_{2,3} = \bm{\mathsf{U}}_{2}\bm{\mathsf{U}}_{2}^{*}\bm{\mathsf{K}}_{2,3} = \bm{\mathsf{U}}_{2}\bm{\mathsf{Z}}_{3}^{*},$$
which establishes the low rank factorization of
$\bm{\mathsf{K}}_{2,3}$.[^2] The block $\bm{\mathsf{K}}_{3,2}$ is of
course factorized analogously.

<figure id="fig:peeling">
<table>
<tbody>
<tr>
<td style="text-align: left;">(a) In the first step of the algorithm,
the sibling interaction matrices <span
class="math inline">K<sub>2, 3</sub></span> and <span
class="math inline">K<sub>3, 2</sub></span> on level 1, shown in red,
are compressed.</td>
<td style="text-align: left;"></td>
</tr>
<tr>
<td style="text-align: left;">(b) In the second step, the four sibling
interaction matrices on level 2, shown in red, are compressed. In this
step, we exploit that we now possess factorizations of the gray
blocks.</td>
<td style="text-align: left;"></td>
</tr>
<tr>
<td style="text-align: left;">(c) In the third step, the eight sibling
interaction matrices on level 8, shown in red, are compressed. We again
exploit that we possess factorizations of the gray blocks.</td>
<td style="text-align: left;"></td>
</tr>
</tbody>
</table>
<figcaption>The “peeling algorithm” for computing a HODLR representation
of a matrix described in Section <a href="#sec:HODLRcompression"
data-reference-type="ref"
data-reference="sec:HODLRcompression">20.5</a>.</figcaption>
</figure>

At the second step of the algorithm, the objective is to build
approximations to the sibling matrices at the next finer level, shown in
red in Figure [12](#fig:peeling){reference-type="ref"
reference="fig:peeling"}(b). We use a random matrix
$\bm{\mathsf{\Omega}}$ of the same size, $n\times 2k$, as in the
previous step, but now with four zero blocks, as shown in the figure.
Looking at the sample matrix $\bm{\mathsf{Y}}_{4}$, we find that it
takes the form $$\begin{array}{cccccccccccc}
\bm{\mathsf{Y}}_{4} &=& \bm{\mathsf{K}}_{4,5}&\bm{\mathsf{\Omega}}_{5} &+& \bm{\mathsf{B}}_{4},\\
\tfrac{n}{4} \times k && \tfrac{n}{4} \times \tfrac{n}{4} & \tfrac{n}{4}\times k && \tfrac{n}{4}\times k
\end{array}$$ where $\bm{\mathsf{B}}_{4}$ represents contributions from
interactions between the block $\bm{\mathsf{K}}_{2,3}$ and the random
matrix and $\bm{\mathsf{\Omega}}_{7}$. But now observe that we at this
point are in possession of a low rank approximation to
$\bm{\mathsf{K}}_{2,3}$. This allows us to compute $\bm{\mathsf{B}}_{4}$
explicitly, and then obtain a "pure" sample from the column space of
$\bm{\mathsf{K}}_{4,5}$ by subtracting this contribution out:
$$\bm{\mathsf{Y}}_{4} - \bm{\mathsf{B}}_{4} = \bm{\mathsf{K}}_{4,5}\bm{\mathsf{\Omega}}_{5}.$$
By orthonormalizing the matrix
$\bm{\mathsf{Y}}_{4} - \bm{\mathsf{B}}_{4}$, we obtain the ON matrix
$\bm{\mathsf{U}}_{4}$ whose columns span the column space of
$\bm{\mathsf{K}}_{4,5}$. The same procedure of course also allows us to
build bases for the column spaces of all matrices at this level, and we
can then extract bases for the row spaces and the cores of the matrices
by sampling $\bm{\mathsf{K}}^{*}$.

At the third step of the algorithm, the object is to compress the blocks
marked in red in Figure [12](#fig:peeling){reference-type="ref"
reference="fig:peeling"}(c). The random matrix used remains of size
$n\times 2k$ but now contains 8 zero blocks, as shown in the figure. The
sample matrix $\bm{\mathsf{Y}}_{8}$ takes the form
$$\begin{array}{cccccccccccc}
\bm{\mathsf{Y}}_{8} &=& \bm{\mathsf{K}}_{8,9}&\bm{\mathsf{\Omega_{9}}} &+& \bm{\mathsf{B}}_{8}\\
\tfrac{n}{8} \times k && \tfrac{n}{8} \times \tfrac{n}{8} & \tfrac{n}{8}\times k && \tfrac{n}{8}\times k
\end{array}$$ where $\bm{\mathsf{B}}_{8}$ represents contributions from
off-diagonal blocks at the coarser levels, shown in gray in the figure.
Since we already possess data sparse representations of these blocks, we
can subtract their contributions out, just as at the previous level.

Once all levels have been processed, we hold low-rank factorizations of
all off-diagonal blocks, and all that remains is to extract the diagonal
matrices $\bm{\mathsf{D}}_{\tau}$ for all leaf nodes $\tau$. Let us for
simplicity assume that all such block are of the same size $m\times m$.
We then form a test matrix by stacking $n/m$ copies of an $m\times m$
identity matrix atop of each other and applying $\bm{\mathsf{K}}$ to
this test matrix. We subtract off the contributions from all the
off-diagonal blocks using the representation we have on hand, and are
then left with a sample matrix consisting of all the blocks
$\bm{\mathsf{D}}_{\tau}$ stacked on top of each other. ◻
:::

![The process used to complete the factorization of the blocks
$\bm{\mathsf{K}}_{2,3}$ and $\bm{\mathsf{K}}_{3,2}$ once the basis
matrices $\bm{\mathsf{U}}_{2}$ and $\bm{\mathsf{U}}_{3}$ have been
constructed. We put $\bm{\mathsf{U}}_{2}$ and $\bm{\mathsf{U}}_{3}$ into
the $n\times 2k$ test matrix shown, and form the sample matrix that
holds $\bm{\mathsf{Z}}_{2}$ and $\bm{\mathsf{Z}}_{3}$ by applying
$\bm{\mathsf{K}}^{*}$.](Pics/fig_peeling_level0_trans.pdf){#fig:peeling_trans
width="85mm"}

## A linear complexity data sparse format {#sec:HBS}

The "HODLR" data sparse matrix format discussed in Sections
[20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"} and
[20.5](#sec:HODLRcompression){reference-type="ref"
reference="sec:HODLRcompression"} is simple to use, and is in many
applications efficient enough. However, its storage complexity is at
least $O(n\log n)$ as the matrix size $n$ grows. We will next describe
how the logarithmic factor can be eliminated. To keep the presentation
concise, we focus on the asymptotic storage requirements, although
similar estimates hold for the asymptotic flop count.

To start, let us first investigate why there are logarithmic factors in
the storage requirement in the HODLR format. Suppose that for a sibling
pair $\{\alpha,\,\beta\}$ the corresponding off-diagonal block
$\bm{\mathsf{K}}_{\alpha,\beta}$ is of size $m\times m$, has rank $k$,
and is stored in the form of a factorization
$$\begin{array}{ccccccccccccccc}
\bm{\mathsf{K}}_{\alpha,\beta} &=& \bm{\mathsf{U}}^{\mathrm{long}}_{\alpha}&\tilde{\bm{\mathsf{K}}}_{\alpha,\beta}&\bigl(\bm{\mathsf{V}}_{\beta}^{\mathrm{long}}\bigr)^{*}\\
m\times m && m\times k & k\times k& k\times m
\end{array}$$ where $\bm{\mathsf{U}}_{\alpha}^{\mathrm{long}}$ and
$\bm{\mathsf{V}}_{\beta}^{\mathrm{long}}$ are two matrices whose columns
form bases for the column and row space of
$\bm{\mathsf{K}}_{\alpha,\beta}$, respectively.[^3] We see that
$\bm{\mathsf{K}}_{\alpha,\beta}$ can be stored using about $2mk$
floating point numbers (where we ignore $O(k^{2})$ terms, since
typically $k \ll m$). Now let us consider how many floats are required
to store all sibling interaction matrices on a given level $\ell$ in the
tree. There are $2^{\ell}$ such matrices, and they each have size
$m \approx 2^{-\ell}n$, resulting in a total of
$2^{\ell}\times 2\times 2^{-\ell}nk = 2nk$ floats. Since there are
$O(\log(n))$ levels in the tree, it follows that just the task of
storing all the basis matrices requires $O(kn\log(n))$ storage.

A standard technique for overcoming the problem of storing all the
"long" basis matrices is to assume that the basis matrices on one level
can be formed through slight modifications to the basis matrices on the
next finer level. To be precise, we assume that if $\tau$ is a node in
the tree with children $\alpha$ and $\beta$, then there exists a "short"
basis matrix $\bm{\mathsf{U}}_{\tau}$ of size $2k\times k$ such that
$$\begin{equation}
\label{eq:nestedbasis}
\begin{array}{ccccccccccccc}
\bm{\mathsf{U}}_{\tau}^{\mathrm{long}} &=&
\left[\begin{array}{cc} \bm{\mathsf{U}}_{\alpha}^{\mathrm{long}} & \bm{\mathsf{0}} \\ \bm{\mathsf{0}} & \bm{\mathsf{U}}_{\beta}^{\mathrm{long}}\end{array}\right]
&
\bm{\mathsf{U}}_{\tau}\\
m\times m && m\times 2k & 2k \times k
\end{array}
\end{equation}$$ (with of course an analogous statement holding for the
"V" matrices holding bases for the row spaces). The point is that if we
already have the "long" basis matrices for the children available, then
the long basis matrices for the parent can be formed using the
information in the small matrix $\bm{\mathsf{U}}_{\tau}$. By applying
this argument recursively, we see that for any parent node $\tau$, all
that needs to be stored explicitly are small matrices of size
$2k\times k$.

To illustrate the idea, suppose that we use this technique for storing
basis matrices to the tree with three levels described in Section
[20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"}. Then for
each of the 8 leaf boxes, we compute basis matrices of size
$n/8 \times k$ and store these. At the next coarser level, the long
basis matrix for a box such as $\tau=4$ can be expressed through the
formula $$\begin{array}{ccccccccccccc}
\bm{\mathsf{U}}_{4}^{\mathrm{long}} &=&
\left[\begin{array}{cc} \bm{\mathsf{U}}_{8}^{\mathrm{long}} & \bm{\mathsf{0}} \\ \bm{\mathsf{0}} & \bm{\mathsf{U}}_{9}^{\mathrm{long}}\end{array}\right]
&
\bm{\mathsf{U}}_{4},\\
(n/4)\times m && (n/4)\times 2k & 2k \times k
\end{array}$$ so that only the $2k\times k$ matrix $\bm{\mathsf{U}}_{4}$
needs to be stored. At the next coarser level still, the long basis
matrix $\bm{\mathsf{U}}_{2}^{\mathrm{long}}$ is expressed as
$$\begin{array}{ccccccccccccc}
\bm{\mathsf{U}}_{2}^{\mathrm{long}} &=&
\left[\begin{array}{cccc}
\bm{\mathsf{U}}_{8}^{\mathrm{long}}  & \bm{\mathsf{0}}                 & \bm{\mathsf{0}}                 & \bm{\mathsf{0}} \\
\bm{\mathsf{0}}                 & \bm{\mathsf{U}}_{9}^{\mathrm{long}}  & \bm{\mathsf{0}}                 & \bm{\mathsf{0}} \\
\bm{\mathsf{0}}                 & \bm{\mathsf{0}}                 & \bm{\mathsf{U}}_{10}^{\mathrm{long}} & \bm{\mathsf{0}} \\
\bm{\mathsf{0}}                 & \bm{\mathsf{0}}                 & \bm{\mathsf{0}}                 & \bm{\mathsf{U}}_{11}^{\mathrm{long}}
\end{array}\right]
&
\left[\begin{array}{cc} \bm{\mathsf{U}}_{4} & \bm{\mathsf{0}} \\ \bm{\mathsf{0}} & \bm{\mathsf{U}}_{5}\end{array}\right]
&
\bm{\mathsf{U}}_{2}.\\
(n/2)\times m && (n/2)\times 4k & 4k\times 2k & 2k \times k
\end{array}$$

We say that a matrix that can be represented using nested basis matrices
in the manner described is a *hierarchically block separable (HBS)*
matrix. (This format is very closely related to the *hierarchically semi
separable (HSS)* matrix format described in
[@2010_gu_xia_HSS; @2005_gu_HSS].) When the HBS format is used, it is
natural to assume that each leaf index vector holds $O(k)$ indices,
which means that there are overall $O(n/k)$ nodes in the tree. Since we
only need $O(k^2)$ storage per node, we see that the overall storage
requirements improve from $O(kn\log(n))$ to $O(kn)$.

There is a price to be paid for using "nested" basis matrices, and it is
that the rank $k$ required to reach a requested precision $\varepsilon$
typically is higher when nested basis matrices are used, in comparison
to the HODLR format. To be precise, one can show that for the relation
([\[eq:nestedbasis\]](#eq:nestedbasis){reference-type="ref"
reference="eq:nestedbasis"}) to hold, it is necessary and sufficient
that for a node $\alpha$, the columns of
$\bm{\mathsf{U}}_{\alpha}^{\mathrm{long}}$ must span the column space of
the matrix $\bm{\mathsf{K}}(I_{\alpha},I_{\alpha}^{\mathrm{c}})$ (where
$I_{\alpha}^{\mathrm{c}} = I \backslash I_{\alpha}$). In contrast, in
the HODLR format, they only need to span the columns of
$\bm{\mathsf{K}}(I_{\alpha},I_{\beta})$, where $\beta$ is the sibling of
$\alpha$. This is easier to do since $I_{\beta}$ is a subset of the
index vector $I_{\alpha}^{\mathrm{c}}$.

## A linear complexity randomized compression technique {#sec:HBScompression}

We saw in Section [20.6](#sec:HBS){reference-type="ref"
reference="sec:HBS"} that the so called "HBS" data sparse format allows
us to store a matrix using $O(kn)$ floating point numbers, without any
factor of $\log(n)$. But is it also possible to compute such a
representation in optimal linear complexity? The answer is yes, and a
number of deterministic techniques that work in specific environments
have been proposed [@2019_book Ch. 17]. A linear complexity randomized
technique was proposed in [@2011_martinsson_randomhudson]. This method
is based on analyzing samples of the matrix obtained through the
application of $\bm{\mathsf{K}}$ and $\bm{\mathsf{K}}^{*}$ to Gaussian
random vectors, but is not quite a black-box method, as it also requires
the ability to evaluate $O(kn)$ individual elements of $\bm{\mathsf{K}}$
itself. This is not always possible, of course, but when it is, the
technique is lightning fast in practice. To be precise,
[@2011_martinsson_randomhudson] establishes the following:

::: {#thm:randHBS .theorem}
**Theorem 72**. *Let $\bm{\mathsf{K}}$ be an $n\times n$ matrix that is
compressible in the HBS format described in Section
[20.6](#sec:HBS){reference-type="ref" reference="sec:HBS"}, for some
rank $k$. Let $T_{\mathrm{apply}}$ denote the time it takes to evaluate
the two products $$\begin{align}
\label{eq:randomhudsonsample1}
\bm{\mathsf{Y}} =&\ \bm{\mathsf{K}}\bm{\mathsf{G}},\\
\label{eq:randomhudsonsample2}
\bm{\mathsf{Z}} =&\ \bm{\mathsf{K}}^{*}\bm{\mathsf{G}},
\end{align}$$ where $\bm{\mathsf{G}}$ is an $n\times k$ matrix drawn
from a Gaussian distribution. Then a full HBS representation of
$\bm{\mathsf{K}}$ can be computed at cost bounded by
$$T_{\mathrm{total}} = T_{\mathrm{apply}} + T_{\mathrm{entry}} \times O(nk) + T_{\mathrm{flop}} \times O(nk^{2}),$$
where $T_{\mathrm{flop}}$ is the cost of a floating point operation, and
where $T_{\mathrm{entry}}$ is the time it takes to evaluate an
individual entry of $\bm{\mathsf{K}}$*
:::

When the off-diagonal blocks are only approximately of rank $k$, we do
some over sampling as usual, and replace the number $k$ in the theorem
by $k+p$ for some small number $p$. This of course results in an
approximate HBS representation of $\bm{\mathsf{K}}$.

The proof of Theorem [72](#thm:randHBS){reference-type="ref"
reference="thm:randHBS"} is an algorithm that explicitly builds the data
sparse representation of the matrix within the specified time budget.
The algorithm consists of a pass through all nodes in the hierarchical
tree, going from smaller boxes to larger; "bottom-up", as opposed to the
"top-down" method for HODLR matrices in Section
[20.5](#sec:HODLRcompression){reference-type="ref"
reference="sec:HODLRcompression"}. We provide an outline of the proof
below. This outline is written to convey the main ideas without
introducing cumbersome notation; for full details, see the original
article [@2011_martinsson_randomhudson] or [@2019_book Sec. 17.4].

::: proof
*Proof.* To describe the algorithm that establishes Theorem
[72](#thm:randHBS){reference-type="ref" reference="thm:randHBS"}, we
will walk through how it applies to a matrix tessellated in accordance
with the hierarchical tree in Figure
[11](#fig:tree){reference-type="ref" reference="fig:tree"}(a). We start
by showing how we can use the information in the matrix
$\bm{\mathsf{Y}}$ defined by
([\[eq:randomhudsonsample1\]](#eq:randomhudsonsample1){reference-type="ref"
reference="eq:randomhudsonsample1"}) to build the basis matrix
$\bm{\mathsf{U}}_{8}$. As we saw in Section
[20.6](#sec:HBS){reference-type="ref" reference="sec:HBS"}, we need the
columns of $\bm{\mathsf{U}}_{8}$ to span the columns in the submatrix
$\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}})$. We achieve this by
constructing the sample matrix $$\bm{\mathsf{Y}}_{8} :=
\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}})\bm{\mathsf{G}}(I_{8}^{\mathrm{c}},:).$$
Since $\bm{\mathsf{G}}(I_{8}^{\mathrm{c}},:)$ is a Gaussian random
matrix of the appropriate size, the columns of $\bm{\mathsf{Y}}_{8}$
will form an approximate basis for the column space of
$\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}})$. But how do we form
$\bm{\mathsf{Y}}_{8}$ from the sample matrix $\bm{\mathsf{Y}}$ defined
by
([\[eq:randomhudsonsample1\]](#eq:randomhudsonsample1){reference-type="ref"
reference="eq:randomhudsonsample1"})? The method we use is illustrated
in Figure [14](#fig:hbs_randomized){reference-type="ref"
reference="fig:hbs_randomized"}(a), where the block
$\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}})$ consists of the gray block
inside the red rectangle. We see that in order to isolate the product
between $\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}})$ and
$\bm{\mathsf{G}}(I_{8}^{\mathrm{c}},:)$, all we need to do is to
subtract the contribution from the diagonal block
$\bm{\mathsf{K}}(I_{8},I_{8})$ (marked in white in the figure). In other
words, $$\begin{equation}
\label{eq:Y8}
\bm{\mathsf{Y}}_{8} =
\bm{\mathsf{K}}(:,I_{8})\bm{\mathsf{G}} - \bm{\mathsf{K}}(I_{8},I_{8})\bm{\mathsf{G}}(I_{8},:) =
\bm{\mathsf{Y}}(I_{8},:) - \bm{\mathsf{K}}(I_{8},I_{8})\bm{\mathsf{G}}(I_{8},:).
\end{equation}$$ The formula ([\[eq:Y8\]](#eq:Y8){reference-type="ref"
reference="eq:Y8"}) can be evaluated inexpensively since the block
$\bm{\mathsf{K}}(I_{8},I_{8})$ is small, and can be explicitly formed
since we assume that we have access to individual entries of the matrix
$\bm{\mathsf{K}}$. Once $\bm{\mathsf{Y}}_{8}$ has been formed, we
compute its row interpolatory decomposition to form a basis matrix
$\bm{\mathsf{U}}_{8}$ and an index vector $\tilde{I}_{8}^{\mathrm{row}}$
such that
$$\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}}) = \bm{\mathsf{U}}_{8}\,\bm{\mathsf{K}}(\tilde{I}_{8}^{\mathrm{row}},I_{8}^{\mathrm{c}}).$$
In an entirely analogous way, we form a sample matrix
$\bm{\mathsf{Y}}_{9}$ whose columns span the block
$\bm{\mathsf{K}}(I_{9},I_{9}^{\mathrm{c}})$ through the formula
$$\bm{\mathsf{Y}}_{9} =
\bm{\mathsf{Y}}(I_{9},:) - \bm{\mathsf{K}}(I_{9},I_{9})\bm{\mathsf{G}}(I_{9},:),$$
cf. Figure [14](#fig:hbs_randomized){reference-type="ref"
reference="fig:hbs_randomized"}(b), and then compute the row ID
$\{\bm{\mathsf{U}}_{9},\tilde{I}_{9}^{\mathrm{row}}\}$ of
$\bm{\mathsf{Y}}_{9}$.

<figure id="fig:hbs_randomized">
<embed src="Pics/fig_hbs_randomized.pdf" />
<p>(a) (b) (c)</p>
<figcaption>The algorithm that establishes Theorem <a
href="#thm:randHBS" data-reference-type="ref"
data-reference="thm:randHBS">72</a> is illustrated using a matrix <span
class="math inline">K</span> tessellated in accordance with the tree
shown in Figure <a href="#fig:tree" data-reference-type="ref"
data-reference="fig:tree">11</a>(a). (a) The block <span
class="math inline">K(<em>I</em><sub>8</sub>,:)</span> is marked with a
red rectangle. Within the red rectangle, the shaded area represents the
block <span
class="math inline">K(<em>I</em><sub>8</sub>, <em>I</em><sub>8</sub><sup>c</sup>)</span>.
(b) The block <span
class="math inline">K(<em>I</em><sub>9</sub>,:)</span> is marked with a
red rectangle, with <span
class="math inline">K(<em>I</em><sub>9</sub>, <em>I</em><sub>9</sub><sup>c</sup>)</span>
shaded. (c) The block <span
class="math inline">K(<em>I</em><sub>4</sub>,:)</span> is marked with a
red rectangle, with <span
class="math inline">K(<em>I</em><sub>4</sub>, <em>I</em><sub>4</sub><sup>c</sup>)</span>
shaded. </figcaption>
</figure>

The basis matrices $\bm{\mathsf{V}}_{\tau}$ that span the row spaces of
the offdiagonal blocks for any leaf $\tau$ are built through the same
procedure, but starting with the sample matrix $\bm{\mathsf{Z}}$ defined
by
([\[eq:randomhudsonsample2\]](#eq:randomhudsonsample2){reference-type="ref"
reference="eq:randomhudsonsample2"}). For instance, we form the sample
matrix $$\begin{equation}
\label{eq:Z8}
\bm{\mathsf{Z}}_{8} =
\bm{\mathsf{K}}(I_{8},:)^{*}\bm{\mathsf{G}} - \bm{\mathsf{K}}(I_{8},I_{8})^{*}\bm{\mathsf{G}}(I_{8},:) =
\bm{\mathsf{Z}}(I_{8},:) - \bm{\mathsf{K}}(I_{8},I_{8})^{*}\bm{\mathsf{G}}(I_{8},:),
\end{equation}$$ and then compute the ID of $\bm{\mathsf{Z}}_{8}$ to
build a basis matrix $\bm{\mathsf{V}}_{8}$ and an index vector
$\tilde{I}_{8}^{\mathrm{col}}$ such that
$$\bm{\mathsf{K}}(I_{8}^{\mathrm{c}},I_{8}) = \bm{\mathsf{K}}(I_{8}^{\mathrm{c}},\tilde{I}_{8}^{\mathrm{col}})\bm{\mathsf{V}}_{8}^{*}.$$

The choice to use the interpolatory decomposition in factorizing the
off-diagonal blocks is essential to making the linear complexity
compression scheme work. In particular, observe that once we have formed
the row and column IDs of all the leaf boxes, we automatically obtain
rank-$k$ factorizations of all the corresponding sibling interaction
matrices. For instance, if $\{\alpha,\beta\}$ is a sibling pair
consisting of two leaf nodes, then the $m\times m$ sibling interaction
matrix $\bm{\mathsf{K}}(I_{\alpha},I_{\beta})$ admits the factorization
$$\begin{equation}
\label{eq:randHBSsibs}
\begin{array}{cccccccccccccc}
\bm{\mathsf{K}}(I_{\alpha},I_{\beta}) &=& \bm{\mathsf{U}}_{\alpha} & \bm{\mathsf{K}}(\tilde{I}_{\alpha}^{\mathrm{row}},\tilde{I}_{\beta}^{\mathrm{col}}) & \bm{\mathsf{V}}_{\beta}^{*}.\\
m\times m && m\times k & k\times k & k\times m
\end{array}
\end{equation}$$ In order to evaluate
([\[eq:randHBSsibs\]](#eq:randHBSsibs){reference-type="ref"
reference="eq:randHBSsibs"}), we merely need to form the matrix
$\bm{\mathsf{K}}(\tilde{I}_{\alpha}^{\mathrm{row}},\tilde{I}_{\beta}^{\mathrm{col}})$.

Once the leaf nodes have all been processed, we proceed to the next
coarser level. Consider for instance the node $\tau=4$ with children
$\alpha=8$ and $\beta=9$. Our task is in principle to build the long
basis matrix $\bm{\mathsf{U}}_{4}^{\mathrm{long}}$ that spans the
columns of $\bm{\mathsf{K}}(I_{4},I_{4}^{\mathrm{c}})$, which is the
shaded matrix inside the red rectangle in Figure
[14](#fig:hbs_randomized){reference-type="ref"
reference="fig:hbs_randomized"}(c). (We say "in principle" since it will
not actually be explicitly formed.) To this end, we define the sample
matrix
$$\bm{\mathsf{Y}}_{4} := \bm{\mathsf{K}}(I_{4},I_{4}^{\mathrm{c}})\bm{\mathsf{G}}(I_{4}^{\mathrm{c}},:).$$
The idea is now to repeat the same technique that we used for the leaf
boxes, and write $$\begin{equation}
\label{eq:Y4}
\bm{\mathsf{Y}}_{4} = \bm{\mathsf{Y}}(I_{4},:) - \bm{\mathsf{K}}(I_{4},I_{4})\bm{\mathsf{G}}(I_{4},:).
\end{equation}$$ It turns out that we can evaluate
([\[eq:Y4\]](#eq:Y4){reference-type="ref" reference="eq:Y4"}) very
efficiently: The block $\bm{\mathsf{K}}(I_{4},I_{4})$ is made up of the
two diagonal blocks $\bm{\mathsf{K}}(I_{8},I_{8})$ and
$\bm{\mathsf{K}}(I_{9},I_{9})$, and the two off-diagonal blocks
$\bm{\mathsf{K}}(I_{8},I_{9})$ and $\bm{\mathsf{K}}(I_{9},I_{8})$. Now
observe that at this point in the execution of the algorithm, we have
compressed representations of all these blocks available. Using this
information, we form a short sample matrix $\tilde{\bm{\mathsf{Y}}}_{4}$
of size $2k\times k$ that we can compress to build the basis matrix
$\bm{\mathsf{U}}_{4}$ and the associated index vector
$\tilde{I}_{4}^{\mathrm{row}}$ in an ID of
$\tilde{\bm{\mathsf{Y}}}_{4}$. ◻
:::

## Butterfly matrices {#sec:butterfly}

An interesting class of rank-structured matrices arises from a
generalization of the discrete Fourier transform. These so called
"butterfly matrices" provide a data sparse format for many matrices that
appear in the analysis of neural networks, wave propagation problems,
and signal processing
[@dao2019learning; @2007_oneil_thesis; @2009_candes_demanet_ying_butterfly].
This format involves an additional complication in comparison to the
HODLR and HBS formats in that it requires $O(\log n)$ different
tessellations of the matrix, as illustrated for a simple case in Figure
[15](#fig:butterfly){reference-type="ref" reference="fig:butterfly"}. It
can be demonstrated that when all of the resulting submatrices are of
numerically low rank, the matrix as a whole can be written
(approximately) as a product of $O(\log n)$ sparse matrices, in a manner
analogous to the butterfly representation of an FFT [@1995_briggs_DFT
Sec. 10.4].

![Illustration of the *butterfly* rank-structured matrix format
described in Section [20.8](#sec:butterfly){reference-type="ref"
reference="sec:butterfly"}. The figures show the blocks that must all be
of numerically low rank for a butterfly matrix arising from a binary
tree with 16 leaf nodes.](Pics/fig_butterfly.pdf){#fig:butterfly
width="\\textwidth"}

Randomization has proven to be a powerful tool for finding butterfly
representations of matrices. The techniques involved are more complex
than the methods described in Sections
[20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"} --
[20.6](#sec:HBS){reference-type="ref" reference="sec:HBS"} due to the
multiplicative nature of the representation, and typically involve
iterative refinement rather than direct approximation
[@2015_ying_butterfly; @2017_interpolative_butterfly; @2018_ying_multidimensional_butterfly].
Similar techniques played an essential role in a recent ground breaking
paper [@2017_michielssen_butterfly_IE] that exploits butterfly
representations to directly solve linear systems arising from the
modeling of scattering problems in the high frequency regime.

## Applications of rank-structured matrices in data analysis {#sec:ASKIT}

The machinery for working with rank-structured hierarchical matrices
that we have described can be used also for the kernel matrices
discussed in Section [19](#sec:kernel){reference-type="ref"
reference="sec:kernel"} that arise in machine learning and computational
statistics. For instance, demonstrate that data sparse formats of this
type can be very effective for simulating Gaussian processes in low
dimensional spaces.

When the underlying dimension grows, the techniques that we have
described so far become uncompetitive. Even $d=4$ would be considered a
stretch. Fortunately, significant progress has recently been made
towards extending the essential ideas to higher dimensions. For
instance, describe a technique that is designed to uncover intrinsic
lower dimensional structures that are often present in sets of points
$\{\bm{\mathsf{x}}_{i}\}_{i=1}^{n}$ that ostensibly live in higher
dimensional spaces. The idea is to use randomized algorithms both for
organizing the points into a hierarchical tree (that induces the
tessellation of the matrix) and for computing low rank approximations to
the resulting admissible blocks. The authors report promising numerical
results for a wide selection of kernel matrices.

In this context, it is as we saw in Section
[19](#sec:kernel){reference-type="ref" reference="sec:kernel"} rarely
possible to execute a matrix-vector multiplication, or even to evaluate
more than a tiny fraction of the entries of the matrix. This means on
the one hand that sampling must form an integral part of the compression
strategy, and on the other that firm performance guarantees are
typically not available. The saving grace is that when a kernel matrix
is used for learning and data analysis, a rough approximation is often
sufficient.

::: remark
**Remark 73** (Geometry oblivious methods). *A curious observation is
that techniques developed for kernel matrices appear to also be
applicable for certain symmetric positive definite (pd) matrices that
are not explicitly presented as kernel matrices. This is a consequence
of the well known fact that any pd matrix $\bm{\mathsf{K}}$ admits a
factorization $$\begin{equation}
\label{eq:pd_K_fact}
\bm{\mathsf{K}} = \bm{\mathsf{G}}^{*}\bm{\mathsf{G}}
\end{equation}$$ for a so called "Gramian matrix" $\bm{\mathsf{G}}$. (If
the eigenvalue decomposition of $\bm{\mathsf{K}}$ takes the form
$\bm{\mathsf{K}} = \bm{\mathsf{U}}\bm{\mathsf{\Lambda}}\bm{\mathsf{U}}^{*}$,
then a matrix of the form
$\bm{\mathsf{G}} = \bm{\mathsf{V}}\bm{\mathsf{\Lambda}}^{1/2}\bm{\mathsf{U}}^{*}$
is a Gramian if and only if $\bm{\mathsf{V}}$ is unitary.) The
factorization ([\[eq:pd_K_fact\]](#eq:pd_K_fact){reference-type="ref"
reference="eq:pd_K_fact"}) says that the entries of $\bm{\mathsf{K}}$
are formed by the inner products between the columns of
$\bm{\mathsf{G}}$, $$\begin{equation}
\label{eq:gofmm_kernel}
\bm{\mathsf{K}}(i,j) = \langle \bm{\mathsf{g}}_{i}, \, \bm{\mathsf{g}}_{j} \rangle = k(\bm{\mathsf{g}}_{i},\bm{\mathsf{g}}_{i}),
\end{equation}$$ where $\bm{\mathsf{g}}_{i}$ is the $i$'th column of
$\bm{\mathsf{G}}$ (the "Gram vector") and where the $k$ is the inner
product kernel of Example [64](#ex:innerprodkernel){reference-type="ref"
reference="ex:innerprodkernel"}. At this point, it becomes plausible
that the techniques of [@2015_biros_ASKIT] for kernel matrices
associated with points in high dimensional spaces may apply to certain
pd matrices. The key to make this work is the observation that it is not
necessary to explicitly form the Gram factors $\bm{\mathsf{G}}$. All
that is needed in order to organize the points
$\{\bm{\mathsf{g}}_{i}\}_{i=1}^{n}$ are relative distances and angles
between the points, and we can evaluate these from the matrix entries of
$\bm{\mathsf{K}}$, via the formula
$$\|\bm{\mathsf{g}}_{i} - \bm{\mathsf{g}}_{j}\|^{2} =
\|\bm{\mathsf{g}}_{i}\|^{2} - 2\mbox{Re}\,\langle \bm{\mathsf{g}}_{i}, \, \bm{\mathsf{g}}_{j} \rangle + \|\bm{\mathsf{g}}_{j}\|^{2} =
\bm{\mathsf{K}}(i,i) - 2\mbox{Re}\,\bm{\mathsf{K}}(i,j) + \bm{\mathsf{K}}(j,j).$$
The resulting technique was presented in [@2017_biros_GOFMM] as a
"geometry oblivious FMM (GOFMM)", along with numerical evidence of its
usefulness for important classes of matrices.*
:::

[^1]: The celebrated method of has exponent
    $\omega = \log_{2}(7) = 2.807\cdots$. It is a compelling algorithm
    in terms of both its numerical stability and its practical speed,
    even for modest matrix sizes. More exotic algorithms, such as the
    Coppersmith--Winograd method and variants, attain complexity of
    about $\omega \approx 2.37$, but they are not considered to be
    practically useful.

[^2]: One may if desired continue the factorization process and form an
    "economy size" SVD of $\bm{\mathsf{Z}}_{3}^{*}$ so that
    $\bm{\mathsf{Z}}_{3}^{*} = \hat{\bm{\mathsf{U}}}_{2}\bm{\mathsf{\Sigma}}_{2,3}\bm{\mathsf{V}}_{3}^{*}$.
    This results in the factorization
    $\bm{\mathsf{K}}_{2,3} = \bigl(\bm{\mathsf{U}}_{2}\hat{\bm{\mathsf{U}}}_{2}\bigr)\bm{\mathsf{\Sigma}}_{2,3}\bm{\mathsf{V}}_{3}^{*}$,
    which is an SVD of $\bm{\mathsf{K}}_{2,3}$. However, for purposes of
    establishing the theorem, we may leave $\bm{\mathsf{Z}}_{3}$ alone,
    and let $\{\bm{\mathsf{U}}_{2},\bm{\mathsf{Z}}_{3}\}$ be our
    "compressed" representation of $\bm{\mathsf{K}}_{2,3}$.

[^3]: The matrices $\bm{\mathsf{U}}_{\alpha}^{\mathrm{long}}$ and
    $\bm{\mathsf{V}}_{\beta}^{\mathrm{long}}$ were denoted
    $\bm{\mathsf{U}}_{\alpha}$ and $\bm{\mathsf{V}}_{\beta}$ in Sections
    [20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"} and
    [20.5](#sec:HODLRcompression){reference-type="ref"
    reference="sec:HODLRcompression"}.
