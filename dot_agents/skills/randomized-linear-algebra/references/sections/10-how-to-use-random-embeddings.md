# How to use random embeddings {#sec:overdet-ls}

Algorithm designers have employed random embeddings for many tasks in
linear algebra, optimization, and related areas. Methods based on random
embedding fall into three rough categories: (1) sketch and solve, (2)
iterative sketching, and (3) sketch and precondition. To draw
distinctions among these paradigms, we use each one to derive an
algorithm for solving an overdetermined least-squares problem.

## Overdetermined least-squares

Overdetermined least-squares problems sometimes arise in statistics and
data-analysis applications. We may imagine that some of the data in
these problems is redundant. As such, it seems plausible that we could
reduce the size of the problem to accelerate computation without too
much loss in accuracy.

Consider a matrix $\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$ with
$m \gg n$ and a vector $\bm{\mathsf{b}} \in \mathbb{F}^m$. An
overdetermined least-squares problem has the form $$\begin{equation}
 \label{eqn:overdet-ls}
\underset{\bm{\mathsf{x}} \in \mathbb{F}^n}{\text{minimize}} \
    \frac{1}{2} \Vert  \bm{\mathsf{A}} \bm{\mathsf{x}} - \bm{\mathsf{b}}  \Vert^2.
\end{equation}$$ Following , we can also rewrite the least-squares
problem to emphasize the role of the matrix: $$\begin{equation}
 \label{eqn:overdet-ls2}
\underset{\bm{\mathsf{x}} \in \mathbb{F}^n}{\text{minimize}} \
    \frac{1}{2} \Vert  \bm{\mathsf{A}} \bm{\mathsf{x}}  \Vert^2 - \langle  \bm{\mathsf{x}} , \,  \bm{\mathsf{A}}^* \bm{\mathsf{b}}  \rangle.
\end{equation}$$ Write $\bm{\mathsf{x}}_{\star}$ for an arbitrary
solution to the
problem [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"}.

To make a clear comparison among algorithm design templates, we will
assume that $\bm{\mathsf{A}}$ is dense and unstructured. In this case,
the classical approach to
solving [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"} is based on factorization of the coefficient
matrix (such as QR or SVD) at a cost of $O(mn^2)$ arithmetic operations.

When $\bm{\mathsf{A}}$ is sparse, we would typically use iterative
methods (such as CG), which have a different computational profile. For
sparse matrices, we would also make different design choices in a
sketch-based algorithm. Nevertheless, for simplicity, we will not
discuss the sparse case.

## Subspace embeddings for least-squares

To design a sketching algorithm for the overdetermined least-squares
problem [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"}, we need to construct a subspace embedding
$\bm{\mathsf{S}} \in \mathbb{F}^{d \times m}$ that preserves the
geometry of the range of the matrix
$\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$. In some cases, we may
also need the embedding to preserve the range of the bordered matrix
$\begin{bmatrix} \bm{\mathsf{A}} & \bm{\mathsf{b}} \end{bmatrix} \in \mathbb{F}^{m \times (n+1)}$.

Since the matrix $\bm{\mathsf{A}}$ is dense and unstructured, we will
work with a structured subspace embedding, such as an SRTT
(Section [9.3](#sec:srtt){reference-type="ref" reference="sec:srtt"}).
More precisely, we will assume that evaluating the product
$\bm{\mathsf{SA}}$ costs only $O(mn \log d)$ arithmetic operations. The
best theoretical results for these structured sketches require that the
embedding dimension $d \sim n \log(n) / \varepsilon^2$ to achieve
distortion $\varepsilon$, although the logarithmic factor seems to be
unnecessary in practice.

Throughout this section, we use the heuristic notation $\sim$ to
indicate quantities that are proportional. We also write $\ll$ to mean
"much smaller than."

## Sketch and solve {#sec:sketchandsolve}

The sketch-and-solve paradigm maps the overdetermined least-squares
problem [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"} into a smaller space. Then it uses the
solution to the reduced problem as a proxy for the solution to the
original problem. This approach can be very fast, and we only need one
view of the matrix $\bm{\mathsf{A}}$. On the other hand, the results
tend to be very inaccurate.

Let $\bm{\mathsf{S}} \in \mathbb{F}^{d \times m}$ be a subspace
embedding for the range of
$\begin{bmatrix} \bm{\mathsf{A}} & \bm{\mathsf{b}} \end{bmatrix}$ with
distortion $\varepsilon$. Consider the compressed least-squares problem
$$\begin{equation}
 \label{eqn:sketch-solve-ls}
\underset{\bm{\mathsf{x}} \in \mathbb{R}^n}{\text{minimize}} \
    \frac{1}{2} \Vert  \bm{\mathsf{S}} (\bm{\mathsf{A}} \bm{\mathsf{x}} - \bm{\mathsf{b}})  \Vert^2.
\end{equation}$$ Since $\bm{\mathsf{S}}$ preserves geometry, we may hope
that the solution $\widehat{\bm{\mathsf{x}}}$ to the sketched
problem [\[eqn:sketch-solve-ls\]](#eqn:sketch-solve-ls){reference-type="eqref"
reference="eqn:sketch-solve-ls"} can replace the solution
$\bm{\mathsf{x}}_{\star}$ to the original
problem [\[eqn:overdet-ls2\]](#eqn:overdet-ls2){reference-type="eqref"
reference="eqn:overdet-ls2"}. A typical theoretical bound is
$$\Vert  \bm{\mathsf{A}} \widehat{\bm{\mathsf{x}}} - \bm{\mathsf{b}}  \Vert \leq (1 + \varepsilon) \Vert  \bm{\mathsf{A}} \bm{\mathsf{x}}_{\star} - \bm{\mathsf{b}}  \Vert
\quad\text{when $d \sim n \log(n) / \varepsilon^2$.}$$ See . Although
the residuals are comparable, it need *not* be the case that
$\widehat{\bm{\mathsf{x}}} \approx \bm{\mathsf{x}}_{\star}$, even when
the solution
to [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"} is unique.

The sketch-and-solve paradigm requires us to form the matrix
$\bm{\mathsf{SA}}$ at a cost of $O(mn \log d)$ operations. We would
typically solve the (dense) reduced problem with a direct method, using
$O(dn^2)$ operations. Assuming $d \sim n \log(n)/\varepsilon^2$, the
total arithmetic cost is
$O(mn \log(n/\varepsilon^2) + n^3 \log(n)/\varepsilon^2)$.

In summary, we witness an improvement in computational cost over
classical methods if $\log n \ll n \ll m/\log n$ and $\varepsilon$ is
constant. But we must also be willing to accept large errors, because we
cannot make $\varepsilon$ small.

::: remark
**Remark 28** (History). *The sketch-and-solve paradigm is attributed to
. It plays a major role in the theoretical algorithms literature; see 
for advocacy. It has also been proposed for enormous problems that might
otherwise be entirely hopeless [@2017_weare_randomized_iteration].*
:::

## Iterative sketching

Iterative sketching attempts to remediate the poor accuracy of the
sketch-and-solve paradigm by applying it repeatedly to reduce the
residual error.

First, we construct an initial solution
$\bm{\mathsf{x}}_0 \in \mathbb{F}^n$ using the sketch-and-solve paradigm
with a constant distortion embedding. For each iteration $i$, draw a
fresh random subspace embedding
$\bm{\mathsf{S}}_i \in \mathbb{F}^{d\times m}$ for
$\operatorname{range}(\bm{\mathsf{A}})$, with constant distortion. We
can solve a sequence of least-squares problems $$\begin{equation}
 \label{eqn:iterative-sketch-ls}
\underset{\bm{\mathsf{x}} \in \mathbb{R}^n}{\text{minimize}} \
    \frac{1}{2} \Vert  \bm{\mathsf{S}}_i \bm{\mathsf{A}}(\bm{\mathsf{x}} - \bm{\mathsf{x}}_{i-1})  \Vert^2
    + \langle  \bm{\mathsf{x}} - \bm{\mathsf{x}}_{i-1} , \,  \bm{\mathsf{A}}^* (\bm{\mathsf{b}} - \bm{\mathsf{A}}\bm{\mathsf{x}}_{i-1})  \rangle.
\end{equation}$$ The solution $\bm{\mathsf{x}}_i$ to this subproblem is
fed into the next subproblem. Without the sketch, each subproblem is
equivalent to
solving [\[eqn:overdet-ls2\]](#eqn:overdet-ls2){reference-type="eqref"
reference="eqn:overdet-ls2"} with $\bm{\mathsf{b}}$ replaced by the
residual
$\bm{\mathsf{r}}_{i-1} = \bm{\mathsf{b}} - \bm{\mathsf{A}} \bm{\mathsf{x}}_{i-1}$.
The sketch $\bm{\mathsf{S}}_i$ preserves the geometry while reducing the
problem size. A typical theoretical error bound would be
$$\Vert  \bm{\mathsf{A}}\bm{\mathsf{x}}_j - \bm{\mathsf{b}}  \Vert \leq (1 + \varepsilon) \Vert \bm{\mathsf{A}}\bm{\mathsf{x}}_{\star} - \bm{\mathsf{b}} \Vert
\quad\text{when $j \sim \log(1/\varepsilon)$ and $d \sim n \log n$.}$$
See [@PW16:Iterative-Hessian] for related results.

In each iteration, the iterative sketching approach requires us to form
$\bm{\mathsf{S}}_i \bm{\mathsf{A}}$ at a cost of $O(mn \log d)$. We
compute $\bm{\mathsf{A}}^* \bm{\mathsf{r}}_{i-1}$ at a cost of $O(mn)$.
Although it is unnecessary to solve each subproblem accurately, we
cannot obtain reliable behavior without using a dense method at a cost
of $O(dn^2)$ per iteration. With the theoretical parameter choices, the
total arithmetic is $O((mn + n^3) \log(n) \log(1/\varepsilon))$ to
achieve relative error $\varepsilon$.

The interesting parameter regime is $\log n \ll n \ll m / \log n$, but
we can now allow $\varepsilon$ to be tiny. In this setting, iterative
sketching costs slightly more than the sketch-and-solve paradigm to
achieve constant relative error, while it is faster than the classical
approach. At the same time, it can produce errors as small as
traditional least-squares algorithms. A shortcoming is that this method
requires repeated sketches of the matrix $\bm{\mathsf{A}}$.

For overdetermined least-squares, we can short-circuit the iterative
sketching approach. In this setting, we can sketch the input matrix just
once and factorize it. We can use the same factorized sketch in each
iteration to solve the subproblems faster. For problems more general
than least-squares, it may be necessary to extract a fresh sketch at
each iteration, as we have done here.

::: remark
**Remark 29** (History). *Iterative sketching can be viewed as an
extension of stochastic approximation methods from optimization, for
example stochastic gradient descent [@Bot10:Large-Scale-Machine]. In the
context of randomized NLA, these algorithms first appeared in the guise
of the randomized Kaczmarz iteration [@SV09:Randomized-Kaczmarz]; see
Section [17.4](#sec:rk){reference-type="ref" reference="sec:rk"}.
reinterpreted randomized Kaczmarz as an iterative sketching method and
developed generalizations. proposed a similar method for solving
overdetermined least-squares problems with constraints; they observed
that better numerical performance is obtained by
sketching [\[eqn:overdet-ls2\]](#eqn:overdet-ls2){reference-type="eqref"
reference="eqn:overdet-ls2"} instead
of [\[eqn:overdet-ls\]](#eqn:overdet-ls){reference-type="eqref"
reference="eqn:overdet-ls"}.*
:::

## Sketch and precondition {#sec:sketchandprecond}

The sketch-and-precondition paradigm uses random embedding to find a
proxy for the input matrix. We can use this proxy to precondition a
classical iterative algorithm so it converges in a minimal number of
iterations.

Let $\bm{\mathsf{S}} \in \mathbb{F}^{d \times m}$ be a subspace
embedding for $\operatorname{range}(\bm{\mathsf{A}})$ with constant
distortion. Compress the input matrix $\bm{\mathsf{A}}$, and then
compute a (pivoted) QR factorization:
$$\bm{\mathsf{Y}} = \bm{\mathsf{SA}}
\quad\text{and}\quad
\bm{\mathsf{Y}} = \bm{\mathsf{QR}}.$$ Since $\bm{\mathsf{S}}$ preserves
the range of $\bm{\mathsf{A}}$ when $d \sim n \log n$, we anticipate
that
$\bm{\mathsf{Y}}^* \bm{\mathsf{Y}} \approx \bm{\mathsf{A}}^*\bm{\mathsf{A}}$.
As a consequence, $\bm{\mathsf{A}}\bm{\mathsf{R}}^\dagger$ should be
close to an isometry. Thus, we can pass to the preconditioned problem
$$\begin{equation}
 \label{eqn:sketch-precond-ls}
\underset{\bm{\mathsf{x}} \in \mathbb{R}^n}{\text{minimize}} \
    \frac{1}{2} \Vert  (\bm{\mathsf{A}}\bm{\mathsf{R}}^\dagger)(\bm{\mathsf{R}}\bm{\mathsf{x}}) - \bm{\mathsf{b}}  \Vert^2.
\end{equation}$$ Construct an initial solution
$\bm{\mathsf{x}}_0 \in \mathbb{F}^n$ using the sketch-and-solve paradigm
with the embedding $\bm{\mathsf{S}}$. From this starting point, we
solve [\[eqn:sketch-precond-ls\]](#eqn:sketch-precond-ls){reference-type="eqref"
reference="eqn:sketch-precond-ls"} using preconditioned LQSR. The $j$th
iterate satisfies
$$\Vert  \bm{\mathsf{A}} \bm{\mathsf{x}}_j - \bm{\mathsf{b}}  \Vert \leq (1 + \varepsilon) \Vert  \bm{\mathsf{A}} \bm{\mathsf{x}}_{\star} - \bm{\mathsf{b}}  \Vert
\quad\text{when $j \sim \log(1/\varepsilon)$.} %$$ This statement is a
reinterpretation of the theory in .

The cost of sketching the input matrix and performing the QR
decomposition is $O(mn \log d + dn^2)$. Afterwards, we pay $O(mn)$ for
each iteration of PCG. With the theoretical parameter settings, the
total cost is $O(mn \log(nB/\varepsilon) + n^3 \log n)$ operations.

Once again, the interesting regime is $\log n \ll n \ll m / \log n$, and
the value of $\varepsilon$ can be very small. For overdetermined
least-squares, this approach is faster than both the sketch-and-solve
paradigm and the iterative sketching paradigm. The
sketch-and-precondition approach leads to errors that are comparable
with classical linear algebra algorithms, but it may be a factor of
$n / \log(n)$ faster. On the other hand, it requires repeated
applications of the matrix $\bm{\mathsf{A}}$.

::: remark
**Remark 30** (History). *The randomized preconditioning idea was
proposed by . demonstrate that least-squares algorithms based on
randomized preconditioning can beat the highly engineered software in
LAPACK. The same method drives the algorithms in . contains a recent
summary of existing randomized preconditioning methods.*
:::

## Comparisons

If we seek a high-precision solution to a dense, unstructured,
overdetermined least-squares problem, randomized preconditioning leads
to the most efficient existing algorithm. For the same problem, if we
can only view the input matrix once, then the sketch-and-solve paradigm
still allows us to obtain a low-accuracy solution. Although iterative
sketching is less efficient than its competitors in this setting, it
remains useful for solving constrained least-squares problems, and it
has further connections with optimization.

## Summary

From the perspective of a numerical analyst, randomized preconditioning
and iterative sketching should be the preferred methods for designing
sketching algorithms because they allow for high precision. The
sketch-and-solve approach is appropriate only when data access is
severely constrained.

In spite of this fact, a majority of the literature on randomized NLA
develops algorithms based on the sketch-and-solve paradigm. There are
far fewer works on randomized preconditioning or iterative sketching.
This discrepancy points to an opportunity for further research.

