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

