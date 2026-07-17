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

