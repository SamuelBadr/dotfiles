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

