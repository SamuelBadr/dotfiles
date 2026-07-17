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

