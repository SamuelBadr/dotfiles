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

