# Single-view algorithms {#sec:singlepass}

In this section, we will describe a remarkable class of algorithms that
are capable of computing a low-rank approximation of a matrix that is so
large that it cannot be stored at all.

We will consider the specific problem of computing an approximate
singular value decomposition of a matrix
$\bm{\mathsf{A}} \in \mathbb{F}^{m\times n}$ under the assumption that
we are allowed to view each entry of $\bm{\mathsf{A}}$ only once and
that we cannot specify the order in which they are viewed. To the best
of our knowledge, no deterministic techniques can carry off such a
computation without *a priori* information about the singular vectors of
the matrix.

For the case where $\bm{\mathsf{A}}$ is psd, we have already seen a
single-view algorithm: the Nyström technique of
Algorithm [\[alg:random-nystrom\]](#alg:random-nystrom){reference-type="ref"
reference="alg:random-nystrom"}. Here, we concentrate on the more
difficult case of general matrices. This presentation is adapted from
and the papers [@TYUC17:Practical-Sketching] and
[@TYUC19:Streaming-Low-Rank].

## Algorithms {#sec:singleviewalgorithms}

In the basic RSVD algorithm (Section
[11.2](#sec:rsvd){reference-type="ref" reference="sec:rsvd"}), we view
each element of the given matrix $\bm{\mathsf{A}}$ at least twice. In
the first view, we form a sample matrix
$\bm{\mathsf{Y}} = \bm{\mathsf{A}}\bm{\mathsf{\Omega}}$ for a given test
matrix $\bm{\mathsf{\Omega}}$. We orthonormalize the columns of
$\bm{\mathsf{Y}}$ to form the matrix $\bm{\mathsf{Q}}$ and then visit
$\bm{\mathsf{A}}$ again to form a second sample
$\bm{\mathsf{C}} = \bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}$. The columns of
$\bm{\mathsf{Y}}$ form an approximate basis for the column space of
$\bm{\mathsf{A}}$, and the columns of $\bm{\mathsf{C}}$ form an
approximate basis for the row space.

In the single-view framework, we can only visit $\bm{\mathsf{A}}$ once,
which means that we must sample both the row and the column space
simultaneously. To this end, let us draw tall thin random matrices
$$\begin{equation}
\label{eq:single1}
\bm{\mathsf{\Upsilon}} \in \mathbb{F}^{m\times \ell}
\qquad\mbox{and}\qquad
\bm{\mathsf{\Omega}} \in \mathbb{F}^{n\times \ell}
\end{equation}$$ and then form the two corresponding sample matrices
$$\begin{equation}
\label{eq:single2}
\bm{\mathsf{X}} = \bm{\mathsf{A}}^{*}\bm{\mathsf{\Upsilon}} \in \mathbb{F}^{n\times \ell}
\qquad\mbox{and}\qquad
\bm{\mathsf{Y}} = \bm{\mathsf{A}}\bm{\mathsf{\Omega}} \in \mathbb{F}^{m\times \ell}.
\end{equation}$$ In ([\[eq:single1\]](#eq:single1){reference-type="ref"
reference="eq:single1"}), we draw a number $\ell$ of samples that is
slightly larger than the rank $k$ of the low-rank approximation that we
seek. (Section [15.2](#sec:singleapriori){reference-type="ref"
reference="sec:singleapriori"} gives details about how to choose
$\ell$.) Observe that both $\bm{\mathsf{X}}$ and $\bm{\mathsf{Y}}$ can
be formed in a single pass over the matrix $\bm{\mathsf{A}}$.

Once we have seen the entire matrix, the next step is to orthonormalize
the columns of $\bm{\mathsf{X}}$ and $\bm{\mathsf{Y}}$ to obtain
orthonormal matrices $$\begin{equation}
\label{eq:single3}
[\bm{\mathsf{P}},\sim] = \texttt{qr\_econ}(\bm{\mathsf{X}}),
\qquad\mbox{and}\qquad
[\bm{\mathsf{Q}},\sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}}).
\end{equation}$$ At this point, $\bm{\mathsf{P}}$ and $\bm{\mathsf{Q}}$
hold approximate bases for the row and column spaces of
$\bm{\mathsf{A}}$, so we anticipate that $$\begin{equation}
\label{eq:single4}
\bm{\mathsf{A}} \approx \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{P}}\bm{\mathsf{P}}^{*} = \bm{\mathsf{Q}}\bm{\mathsf{C}}\bm{\mathsf{P}}^{*},
\end{equation}$$ where we defined the "core" matrix $$\begin{equation}
\label{eq:single5}
\bm{\mathsf{C}} := \bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{P}} \in \mathbb{F}^{\ell\times \ell}.
\end{equation}$$ Unfortunately, since we cannot revisit
$\bm{\mathsf{A}}$, we are not allowed to form $\bm{\mathsf{C}}$ directly
by applying formula [\[eq:single5\]](#eq:single5){reference-type="eqref"
reference="eq:single5"}.

Instead, we develop a relation that $\bm{\mathsf{C}}$ must satisfy
approximately, which allows us to estimate $\bm{\mathsf{C}}$ from the
quantities we have on hand. To do so, right-multiply the definition
[\[eq:single5\]](#eq:single5){reference-type="eqref"
reference="eq:single5"} by $\bm{\mathsf{P}}^{*}\bm{\mathsf{\Omega}}$ to
obtain $$\begin{equation}
\label{eq:single6}
\bm{\mathsf{C}}\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Omega}}\bigr) = \bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{P}}\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Omega}}\bigr).
\end{equation}$$ Inserting the approximation
$\bm{\mathsf{A}}\bm{\mathsf{P}}\bm{\mathsf{P}}^{*} \approx \bm{\mathsf{A}}$
into ([\[eq:single6\]](#eq:single6){reference-type="ref"
reference="eq:single6"}), we find that $$\begin{equation}
\label{eq:single7}
\bm{\mathsf{C}}\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Omega}}\bigr) \approx \bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{\Omega}} = \bm{\mathsf{Q}}^{*}\bm{\mathsf{Y}}.
\end{equation}$$ In ([\[eq:single7\]](#eq:single7){reference-type="ref"
reference="eq:single7"}), all quantities except $\bm{\mathsf{C}}$ are
known explicitly, which means that we can solve it, in the least-squares
sense, to arrive at an estimate $$\begin{equation}
\label{eq:single8}
\bm{\mathsf{C}}_{\mathrm{approx}} = \bigl(\bm{\mathsf{Q}}^{*}\bm{\mathsf{Y}}\bigr)\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Omega}}\bigr)^{\dagger},
\end{equation}$$ where $\dagger$ denotes the Moore-Penrose
pseudoinverse. As always, the pseudoinverse is applied by means of an
orthogonal factorization. Once $\bm{\mathsf{C}}_{\mathrm{approx}}$ has
been computed via ([\[eq:single8\]](#eq:single8){reference-type="ref"
reference="eq:single8"}), we obtain the rank-$\ell$ approximation
$$\begin{equation}
\label{eq:single8b}
\bm{\mathsf{A}} \approx \bm{\mathsf{Q}}\bm{\mathsf{C}}_{\mathrm{approx}}\bm{\mathsf{P}}^{*},
\end{equation}$$ which we can convert into an approximate SVD using the
standard postprocessing steps. For additional implementation details,
see and . Extensions of this approach, with theoretical analysis, appear
in .

Recently, have demonstrated that the numerical performance of the
single-view algorithm can be improved by extracting a third sketch of
$\bm{\mathsf{A}}$ that is independent from $\bm{\mathsf{X}}$ and
$\bm{\mathsf{Y}}$. The idea is to draw tall thin random matrices
$$\bm{\mathsf{\Phi}} \in \mathbb{F}^{m\times s},
\qquad\mbox{and}\qquad
\bm{\mathsf{\Psi}} \in \mathbb{F}^{n\times s},$$ where $s$ is another
oversampling parameter. Then we form a "core sketch" $$\begin{equation}
\label{eq:single9}
\bm{\mathsf{Z}} = \bm{\mathsf{\Phi}}^{*}\bm{\mathsf{A}}\bm{\mathsf{\Psi}} \in \mathbb{F}^{s\times s}.
\end{equation}$$ This extra data allows us to derive an alternative
equation for the core matrix $\bm{\mathsf{C}}$. We left- and
right-multiply the definition
[\[eq:single5\]](#eq:single5){reference-type="eqref"
reference="eq:single5"} by $\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}}$ and
$\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}}$ to obtain the relation
$$\begin{equation}
\label{eq:single10}
\bigl(\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}}\bigr)\bm{\mathsf{C}}\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}}\bigr) =
\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{P}}\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}}.
\end{equation}$$ Inserting the approximation
$\bm{\mathsf{A}} \approx \bm{\mathsf{Q}}\bm{\mathsf{Q}}^{*}\bm{\mathsf{A}}\bm{\mathsf{P}}\bm{\mathsf{P}}^{*}$
into [\[eq:single10\]](#eq:single10){reference-type="eqref"
reference="eq:single10"}, we find that $$\begin{equation}
\label{eq:single11}
\bigl(\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}}\bigr)\bm{\mathsf{C}}\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}}\bigr) \approx
\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{A}}\bm{\mathsf{\Psi}} = \bm{\mathsf{Z}}.
\end{equation}$$ An improved approximation to the core matrix
$\bm{\mathsf{C}}$ results by solving
[\[eq:single11\]](#eq:single11){reference-type="eqref"
reference="eq:single11"} in a least-squares sense; to wit,
$\bm{\mathsf{C}} = (\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}})^{\dagger}\bm{\mathsf{Z}}(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}})^{\dagger}$.

::: {#rem:streaming .remark}
**Remark 47** (Streaming algorithms). *Single-view algorithms are
related to the streaming model of computation
[@2005_muthukrishnan_stream]. were the first to explicitly study matrix
computations in streaming data models.*

*One important streaming model poses the assumption that the input
matrix $\bm{\mathsf{A}}$ is presented as a sequence of innovations:
$$\bm{\mathsf{A}} = \bm{\mathsf{H}}_{1} + \bm{\mathsf{H}}_{2} + \bm{\mathsf{H}}_{3} + \cdots$$
Typically, each update $\bm{\mathsf{H}}_{i}$ is simple; for instance, it
may be sparse or low-rank. The challenge is that the full matrix
$\bm{\mathsf{A}}$ is too large to be stored. Once an innovation
$\bm{\mathsf{H}}_{i}$ has been processed, it cannot be retained. This is
called the "turnstile" model in the theoretical computer science
literature.*

*The algorithms described in this section handle this difficulty by
creating a random linear transform $\mathcal{S}$ that maps
$\bm{\mathsf{A}}$ down to a low-dimensional sketch that is small enough
to store. What we actually retain in memory is the evolving sketch of
the input: $$\mathcal{S}(\bm{\mathsf{A}}) =
\mathcal{S}(\bm{\mathsf{H}}_{1}) +
\mathcal{S}(\bm{\mathsf{H}}_{2}) +
\mathcal{S}(\bm{\mathsf{H}}_{3}) + \cdots.$$ In Algorithm
[\[alg:SingleViewSVD\]](#alg:SingleViewSVD){reference-type="ref"
reference="alg:SingleViewSVD"}, we instantiate $\mathcal{S}$ by drawing
::: {#alg:SingleViewSVD .algorithm}
**Algorithm. Single-view SVD.**

```
**Input:** Target matrix $\mtx{A} \in \F^{m \times n}$, rank $k$, sampling sizes $\ell$ and $s$.
**Output:** Orthonormal matrices $\mtx{U} \in \F^{m\times k}$ and $\mtx{V} = \F^{n\times k}$, and a diagonal matrix $\mtx{\Sigma} \in \F^{k\times k}$ such that
$\mtx{A} \approx \mtx{U}\mtx{\Sigma}\mtx{V}^{*}$.

function SingleViewSVD($\mtx{A}$, $k$, $\ell$, $s$):
    Draw test matrices $\mtx{\Upsilon} \in \F^{m\times \ell}$, $\mtx{\Omega} \in \F^{n\times \ell}$,
    $\mtx{\Phi} \in \F^{m\times s}$, $\mtx{\Psi} \in \F^{n\times s}$
    Form $\mtx{X} = \mtx{A}^{*}\mtx{\Upsilon}$, $\mtx{Y} = \mtx{A}\mtx{\Omega}$, $\mtx{Z} = \mtx{\Phi}^{*}\mtx{A}\mtx{\Psi}$
    # Viewing $\mtx{A}$ only once!
    $[\mtx{P},\sim] = qr_econ(\mtx{X})$, $[\mtx{Q},\sim] = qr_econ(\mtx{Y})$
    $\mtx{C} = \bigl(\mtx{\Phi}^{*}\mtx{Q}\bigr)^{\pinv}\,\mtx{Z}\,\bigl(\mtx{P}^{*}\mtx{\Psi}\bigr)^{\pinv}$
    # Execute using a least-squares solver
    $[\widehat{\mtx{U}},\widehat{\mtx{\Sigma}},\widehat{\mtx{V}}] = svd(\mtx{C})$
    # A full SVD
    $\mtx{U} = \mtx{Q}\widehat{\mtx{U}}(\colon,1:k)$,
    $\mtx{V} = \mtx{P}\widehat{\mtx{V}}(\colon,1:k)$,
    $\mtx{\Sigma} = \widehat{\mtx{\Sigma}}(1:k,1:k)$
```
:::


the random matrices $\bm{\mathsf{\Upsilon}}$, $\bm{\mathsf{\Omega}}$,
$\bm{\mathsf{\Phi}}$, and $\bm{\mathsf{\Psi}}$, and then work with the
sketch $$\mathcal{S}(\bm{\mathsf{H}}) =
\big(\bm{\mathsf{\Upsilon}}^{*}\bm{\mathsf{H}},\,\bm{\mathsf{H}}\bm{\mathsf{\Omega}},\,\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{H}}\bm{\mathsf{\Psi}} \big).$$
The fact that the sketch is a *linear* map is essential here. prove that
randomized linear embeddings are essentially the only kind of algorithm
for handling the turnstile model. In contrast, the sketch implicit in
the RSVD algorithm from Section [11.2](#sec:rsvd){reference-type="ref"
reference="sec:rsvd"} is a quadratic or higher-order polynomial in the
input matrix.*
:::

:::: algorithm
::: algorithmic
Target matrix $\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$, rank $k$,
sampling sizes $\ell$ and $s$. Orthonormal matrices
$\bm{\mathsf{U}} \in \mathbb{F}^{m\times k}$ and
$\bm{\mathsf{V}} = \mathbb{F}^{n\times k}$, and a diagonal matrix
$\bm{\mathsf{\Sigma}} \in \mathbb{F}^{k\times k}$ such that
$\bm{\mathsf{A}} \approx \bm{\mathsf{U}}\bm{\mathsf{\Sigma}}\bm{\mathsf{V}}^{*}$.

Draw test matrices
$\bm{\mathsf{\Upsilon}} \in \mathbb{F}^{m\times \ell}$,
$\bm{\mathsf{\Omega}} \in \mathbb{F}^{n\times \ell}$,
$\bm{\mathsf{\Phi}} \in \mathbb{F}^{m\times s}$,
$\bm{\mathsf{\Psi}} \in \mathbb{F}^{n\times s}$ Form
$\bm{\mathsf{X}} = \bm{\mathsf{A}}^{*}\bm{\mathsf{\Upsilon}}$,
$\bm{\mathsf{Y}} = \bm{\mathsf{A}}\bm{\mathsf{\Omega}}$,
$\bm{\mathsf{Z}} = \bm{\mathsf{\Phi}}^{*}\bm{\mathsf{A}}\bm{\mathsf{\Psi}}$
Viewing $\bm{\mathsf{A}}$ only once!
$[\bm{\mathsf{P}},\sim] = \texttt{qr\_econ}(\bm{\mathsf{X}})$,
$[\bm{\mathsf{Q}},\sim] = \texttt{qr\_econ}(\bm{\mathsf{Y}})$
$\bm{\mathsf{C}} = \bigl(\bm{\mathsf{\Phi}}^{*}\bm{\mathsf{Q}}\bigr)^{\dagger}\,\bm{\mathsf{Z}}\,\bigl(\bm{\mathsf{P}}^{*}\bm{\mathsf{\Psi}}\bigr)^{\dagger}$
Execute using a least-squares solver
$[\widehat{\bm{\mathsf{U}}},\widehat{\bm{\mathsf{\Sigma}}},\widehat{\bm{\mathsf{V}}}] = \texttt{svd}(\bm{\mathsf{C}})$
A full SVD
$\bm{\mathsf{U}} = \bm{\mathsf{Q}}\widehat{\bm{\mathsf{U}}}(\colon,1:k)$,
$\bm{\mathsf{V}} = \bm{\mathsf{P}}\widehat{\bm{\mathsf{V}}}(\colon,1:k)$,
$\bm{\mathsf{\Sigma}} = \widehat{\bm{\mathsf{\Sigma}}}(1:k,1:k)$
:::
::::

::: remark
**Remark 48** (Single-view versus out-of-core algorithms). *In
principle, the methods discussed in this section can also be used in
situations where a matrix is stored in slow memory, such as a spinning
disk hard drive, or on a distributed memory system. However, one has to
carefully weigh whether the decrease in accuracy and increase in
uncertainty that is inherent to single-view algorithms is worth the cost
savings. As a general matter, revisiting the matrix at least once is
advisable whenever it is possible.*
:::

## Error estimation, parameter choices and truncation {#sec:singleapriori}

In the single-view computing environment, one must choose sampling
parameters before the computation starts, and there is no way to revisit
these choices after data has been gathered. This constraint makes *a
priori* error analysis particularly important, because we need guidance
on how large to make the sketches given some prior knowledge about the
spectral decay of the input matrix. To illustrate how this may work, let
us cite :

::: {#thm:singleview .theorem}
**Theorem 49** (Single-view SVD: Gaussian analysis). *Suppose that
Algorithm
[\[alg:SingleViewSVD\]](#alg:SingleViewSVD){reference-type="ref"
reference="alg:SingleViewSVD"} is executed for an input matrix
$\bm{\mathsf{A}} \in \mathbb{C}^{m\times n}$ and for sampling parameters
$s$ and $\ell$ that satisfy $s \geq 2\ell$. When the test matrices are
drawn from a standard normal distribution, the computed matrices
$\bm{\mathsf{P}}$, $\bm{\mathsf{C}}$, and $\bm{\mathsf{Q}}$ satisfy
$$\operatorname{\mathbb{E}}\Vert  \bm{\mathsf{A}} - \bm{\mathsf{Q}}\bm{\mathsf{C}}\bm{\mathsf{P}}^{*}  \Vert_{\mathrm{F}}^{2}
\leq
\frac{s}{s-\ell}\,
\min_{k < \ell}\left(\frac{\ell+k}{\ell-k}\sum_{j=k+1}^{\min(m,n)}\sigma_{j}^{2}\right).$$
As usual, $\sigma_j$ is the $j$th largest singular value of
$\bm{\mathsf{A}}$. A very similar bound holds for the real field.*
:::

This result suggests that more aggressive oversampling is called for in
the single-view setting, as compared to the basic rangefinder problem.
For instance, if we aim for an approximation error that is comparable to
the best possible approximation with rank $k$, then we might choose
$\ell = 4k$ and $s=8k$ to obtain
$$\mathbb{E}\|\bm{\mathsf{A}} - \bm{\mathsf{Q}}\bm{\mathsf{C}}\bm{\mathsf{P}}^{*}\|_{\mathrm{F}}^{2}
\leq
\frac{10}{3}\sum_{j=k+1}^{\min(m,n)}\sigma_{j}^{2} =
\frac{10}{3}\|\bm{\mathsf{A}} - \bm{\mathsf{A}}_{k}\|_{\mathrm{F}}^{2},$$
where $\bm{\mathsf{A}}_{k}$ is the best possible rank-$k$ approximation
of $\bm{\mathsf{A}}$. As usual, the likelihood of large deviations from
the expectation is negligible. (For contrast, recall that the basic
rangefinder algorithm often works well when we select $\ell = k+5$ or
$\ell=k+10$.)

Besides computing $\bm{\mathsf{P}}$, $\bm{\mathsf{C}}$, and
$\bm{\mathsf{Q}}$, Algorithm
[\[alg:SingleViewSVD\]](#alg:SingleViewSVD){reference-type="ref"
reference="alg:SingleViewSVD"} also prunes the approximation
$\bm{\mathsf{A}} \approx \bm{\mathsf{Q}}\bm{\mathsf{C}}\bm{\mathsf{P}}^{*}$
by computing an SVD of $\bm{\mathsf{C}}$ (line 6) and then throwing out
the trailing $\ell-k$ modes (line 7). The motivation for this truncation
is that the approximation
$\bm{\mathsf{A}} \approx \bm{\mathsf{Q}}\bm{\mathsf{C}}\bm{\mathsf{P}}^{*}$,
tends to capture the dominant singular modes of $\bm{\mathsf{A}}$ well,
but the trailing ones have very low accuracy. The same thing happens
with the basic RSVD (Section [11.2](#sec:rsvd){reference-type="ref"
reference="sec:rsvd"}), but the phenomenon is more pronounced in the
single-view environment, in part because $\ell$ is substantially larger
than $k$. Theorem [49](#thm:singleview){reference-type="ref"
reference="thm:singleview"} can be applied to prove that the truncated
factorization is as accurate as one can reasonably hope for; see  for
details.

::: remark
**Remark 50** (Spectral norm bounds?).
*Theorem [49](#thm:singleview){reference-type="ref"
reference="thm:singleview"} provides a Frobenius norm error bound for a
matrix approximation algorithm. For our survey, this is a *rara avis in
terra*. Unfortunately, relative error spectral norm error bounds are not
generally possible in the streaming setting [@2014_woodruff_sketching
Chap. 6].*
:::

## Structured test matrices

Algorithm [\[alg:SingleViewSVD\]](#alg:SingleViewSVD){reference-type="ref"
reference="alg:SingleViewSVD"} can -- and should -- be implemented with
structured test matrices, rather than Gaussian test matrices. This
modification is especially appealing in the single-view environment,
where storage is often the main bottleneck.

For instance, consider the parameter selections $\ell = 4k$ and $s = 8k$
that we referenced above. Then the four test matrices consist of
$12k(m+n)$ floats that must be stored, and the sketches add another
$4k(m+n) + 64k^{2}$ floats. Since $m$ and $n$ can be huge, these numbers
could severely limit the rank $k$ of the final matrix approximation.

If we swap out the Gaussian matrices for structured random matrices, we
can almost remove the cost associated with storing the test matrices. In
particular, the addition of the core sketch
[\[eq:single9\]](#eq:single9){reference-type="eqref"
reference="eq:single9"} has a very light memory footprint because the
sketch itself only uses $O(k^{2})$ floats. Empirically, when we use a
structured random matrix, such as a sparse sign matrix
(Section [9.2](#sec:sparse-map){reference-type="ref"
reference="sec:sparse-map"}) or an SRTT
(Section [9.3](#sec:srtt){reference-type="ref" reference="sec:srtt"}),
the observed errors are more or less indistinguishable from the errors
attained with Gaussian test matrices. See .

## A posteriori error estimation {#sec:singleaposteriori}

In order to reduce the uncertainty associated with the single-view
algorithms described in this section, the "certificate of accuracy"
technique described in Section
[12.2](#sec:certificate){reference-type="ref"
reference="sec:certificate"} is very useful.

Recall that the idea is to draw a separate test matrix whose only
purpose is to provide an independent estimate of the error in the
computed solution. This additional test matrix can be *very* thin (say 5
or 10 columns wide), and it still provides a dependable bound on the
computed error. These techniques can be incorporated without any
difficulty in the single-view environment, as outlined in .

Let us mention one caveat. In the single-view environment, we have no
recourse when the *a posteriori* error estimator signals that the
approximation error is unacceptable. On the other hand, it is reassuring
that the algorithm can sound a warning that it has not met the desired
accuracy.

## History {#sec:singlehistory}

To the best of our knowledge, described the first algorithm that can
compute a low-rank matrix approximation in the single-view computational
model. Their paper introduced the idea of independently sampling the
row- and the column-space of a matrix, as summarized in formulas
[\[eq:single1\]](#eq:single1){reference-type="eqref"
reference="eq:single1"}--[\[eq:single8b\]](#eq:single8b){reference-type="eqref"
reference="eq:single8b"}. This approach inspired the single-view
algorithms presented in . It is interesting that the primary objective
of was to reduce the asymptotic flop count of the computation through
the use of structured random test matrices.

gave an explicit discussion of randomized NLA in a streaming
computational model. They independently proposed a variant of the
algorithm from [@HMT11:Finding-Structure Sec. 5.5]. Later contributions
to the field appeared in
[@LNW14:Turnstile-Streaming; @2016_boutsidis_optimal; @2016_feldman_dimensionality; @GLPW16:Frequent-Directions]
and [@2017_tropp_practical_sketching]. The idea of introducing an
additional sketch such as
[\[eq:single9\]](#eq:single9){reference-type="eqref"
reference="eq:single9"} to capture the "core" matrix was proposed by .
have provided improvements of his approach, further analysis, and
computational considerations.

