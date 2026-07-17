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

