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

