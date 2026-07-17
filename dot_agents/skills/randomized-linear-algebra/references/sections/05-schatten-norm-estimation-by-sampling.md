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

