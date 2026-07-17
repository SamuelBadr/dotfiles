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

