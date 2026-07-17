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

