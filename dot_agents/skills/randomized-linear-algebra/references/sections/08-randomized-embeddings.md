# Randomized embeddings {#sec:gauss}

One of the core tools in randomized linear algebra is *randomized linear
embedding*, often assigned the misnomer *random projection*. The
application of randomized embeddings is often referred to as
*sketching*.

This section begins with a formal definition of a randomized embedding.
Then we introduce the Gaussian embedding, which is the simplest
construction, and we summarize its analysis. Randomized embeddings have
a wide range of applications in randomized linear algebra. Some
implications of this theory include the Johnson--Lindenstrauss lemma and
a simple construction of a subspace embedding. We also explain why
results for Gaussians transfer to a wider setting. Last, we give a short
description of random partial isometries, a close cousin of Gaussian
embeddings.

## What is a random embedding?

Let $E \subseteq \mathbb{F}^n$ be a set, and let $\varepsilon\in (0, 1)$
be a distortion parameter. We say that a linear map
$\bm{\mathsf{S}} : \mathbb{F}^{n} \to \mathbb{F}^d$ is an ($\ell_2$)
*embedding* of $E$ with distortion $\varepsilon$ when $$\begin{equation}
 \label{eqn:embed-abstract}
(1 - \varepsilon) \, \Vert \bm{\mathsf{x}} \Vert \leq \Vert  \bm{\mathsf{S}} \bm{\mathsf{x}}  \Vert \leq (1 + \varepsilon) \, \Vert  \bm{\mathsf{x}}  \Vert
\quad\text{for all $\bm{\mathsf{x}} \in E$.}
\end{equation}$$ It is sometimes convenient to abbreviate this kind of
two-sided inequality as
$\Vert \bm{\mathsf{S}}\bm{\mathsf{x}} \Vert = (1 \pm \varepsilon) \Vert \bm{\mathsf{x}} \Vert$.

We usually think about the case where $d \ll n$, so the map
$\bm{\mathsf{S}}$ enacts a dimension reduction. In other words,
$\bm{\mathsf{S}}$ transfers data from the high-dimensional space
$\mathbb{F}^n$ to the low-dimensional space $\mathbb{F}^d$. As we will
discuss, the low-dimensional representation of the data can be used to
obtain fast, approximate solutions to computational problems.

The
relation [\[eqn:embed-abstract\]](#eqn:embed-abstract){reference-type="eqref"
reference="eqn:embed-abstract"} expresses the idea that the embedding
$\bm{\mathsf{S}}$ should preserve the geometry of the set $E$.
Unfortunately, we do not always know the set $E$ in advance. Moreover,
we would like the map $\bm{\mathsf{S}}$ to be easy to construct, and it
should be computationally efficient to apply $\bm{\mathsf{S}}$ to the
data. These goals may be in tension.

We can resolve this dilemma by drawing the embedding $\bm{\mathsf{S}}$
from a probability distribution. Many types of probability distributions
serve. In particular, we can use highly structured random matrices
(Section [9](#sec:dimension-reduction){reference-type="ref"
reference="sec:dimension-reduction"}) that are easy to build, to store,
and to apply to vectors.
Section [10](#sec:overdet-ls){reference-type="ref"
reference="sec:overdet-ls"} presents a case study about how random
embeddings can be applied to solve overdetermined least-squares
problems.

## Restricted singular values

Our initial goal is to understand something about the theoretical
behavior of randomized embeddings. To that end, let us introduce
quantities that measure how much an embedding distorts a set. Let
$\bm{\mathsf{S}} \in \mathbb{F}^{d \times n}$ be a linear map, and let
$E \subseteq \mathbb{S}^{n-1}(\mathbb{F})$ be an arbitrary subset of the
unit sphere in $\mathbb{F}^n$. The *minimum* and *maximum restricted
singular value* are, respectively, defined as $$\begin{equation}
 \label{eqn:rsv}
\sigma_{\min}(\bm{\mathsf{S}}; E) := \min_{\bm{\mathsf{x}} \in E}\ \Vert  \bm{\mathsf{S}} \bm{\mathsf{x}}  \Vert
\quad\text{and}\quad
\sigma_{\max}(\bm{\mathsf{S}}; E) := \max_{\bm{\mathsf{x}} \in E}\ \Vert  \bm{\mathsf{S}} \bm{\mathsf{x}}  \Vert.
\end{equation}$$ If $E$ composes the entire unit sphere, then these
quantities coincide with the ordinary minimum and maximum singular value
of $\bm{\mathsf{S}}$. More generally, the restricted singular values
describe how much the linear map $\bm{\mathsf{S}}$ can contract or
expand a point in $E$.

::: remark
**Remark 15** (General sets). *In this treatment, we require $E$ to be a
subset of the unit sphere. Related, but more involved, results hold when
$E$ is a general set. See  and for more results and applications.*
:::

## Gaussian embeddings

Our theoretical treatment of random embeddings focuses on the most
highly structured case. A Gaussian embedding is a random matrix of the
form $$\bm{\mathsf{\Gamma}} \in \mathbb{F}^{d \times n}
\quad\text{with iid entries $(\bm{\mathsf{\Gamma}})_{ij} \sim \textsc{normal}(0, d^{-1})$.}$$
The cost of explicitly storing a Gaussian embedding is $O(dn)$, and the
cost of applying it to a vector is $O(dn)$.

The scaling of the matrix ensures that
$$\operatorname{\mathbb{E}}\Vert  \bm{\mathsf{\Gamma}} \bm{\mathsf{x}}  \Vert^2 = \Vert  \bm{\mathsf{x}}  \Vert^2
\quad\text{for each $\bm{\mathsf{x}} \in \mathbb{F}^n$.}$$ We wish to
understand how large to choose the embedding dimension $d$ so that the
map $\bm{\mathsf{\Gamma}}$ approximately preserves the norms of all
points in a given set $E$. We can do so by obtaining bounds for the
restricted singular values. For a Gaussian embedding, we will see that
$\sigma_{\min}(\bm{\mathsf{\Gamma}}; E)$ and
$\sigma_{\max}(\bm{\mathsf{\Gamma}}; E)$ are controlled by the geometry
of the set $E$.

::: remark
**Remark 16** (Why Gaussians?). *Gaussian embeddings admit a simple and
beautiful analysis. In our computational experience, many other
embeddings exhibit the same (universal) behavior as a Gaussian map. In
spite of that, the rigorous analysis of other types of embeddings tends
to be difficult, even while it yields rather imprecise results. The
confluence of these facts motivates us to argue that the Gaussian
analysis provides enough insight for many practical purposes.*
:::

## The Gaussian width

For the remainder of Section [8](#sec:gauss){reference-type="ref"
reference="sec:gauss"}, we will work in the real field
($\mathbb{F}= \mathbb{R}$). Given a set
$E \subseteq \mathbb{S}^{n-1}(\mathbb{R})$, define the *Gaussian width*
$w(E)$ via
$$w(E) := \operatorname{\mathbb{E}}\sup_{\bm{\mathsf{x}} \in E}\ \langle  \bm{\mathsf{g}} , \,  \bm{\mathsf{x}}  \rangle
\quad\text{where $\bm{\mathsf{g}} \in \mathbb{R}^n$ is standard normal.}$$
The Gaussian width is a measure of the content of the set $E$. It plays
a fundamental role in the performance of randomized embeddings.

Here are some basic properties of the Gaussian width.

- The width is invariant under rotations: $w(\bm{\mathsf{Q}}E) = w(E)$
  for each orthogonal matrix $\bm{\mathsf{Q}}$.

- The width is increasing with respect to set inclusion: $E \subseteq F$
  implies that $w(E) \leq w(F)$.

- The width lies in the range
  $0 \leq w(E) \leq \operatorname{\mathbb{E}}\Vert \bm{\mathsf{g}} \Vert < \sqrt{n}$.

The width can be calculated accurately for many sets of interest. In
particular, if $L$ is an arbitrary $k$-dimensional subspace of
$\mathbb{R}^n$, then $$\begin{equation}
 \label{eqn:width-subspace}
\sqrt{k-1} < w( L \cap \mathbb{S}^{n-1} ) < \sqrt{k}.
\end{equation}$$ Indeed, it is productive to think about the *squared*
width $w^2(E)$ as a measure of the "dimension" of the set $E$.

::: remark
**Remark 17** (Statistical dimension). *The *statistical dimension* is
another measure of content that is closely related to the squared
Gaussian width. The statistical dimension has additional geometric
properties that make it easier to work with in some contexts.
See [@ALMT14:Living-Edge; @MT13:Achievable-Performance; @MT14:Steiner-Formulas]
and [@GNP17:Gaussian-Phase] for more information.*
:::

## Restricted singular values of Gaussian matrices {#sec:rsv-gauss}

In a classic work on Banach space geometry, showed that the Gaussian
width controls both the minimum and maximum restricted singular values
of a subset of the sphere.

::: {#thm:rsv-gauss .theorem}
**Theorem 18** (Restricted singular values: Gaussian matrix). *Fix a
subset $E \subseteq \mathbb{S}^{n-1}(\mathbb{R})$ of the unit sphere.
Draw a Gaussian matrix
$\bm{\mathsf{\Gamma}} \in \mathbb{R}^{d \times n}$ whose entries are iid
$\textsc{normal}(0, d^{-1})$. For all $t > 0$, $$\begin{aligned}
&\mathbb{P}\left\{  \sigma_{\min}(\bm{\mathsf{\Gamma}}; E ) \leq 1 - \frac{w(E) + 1}{\sqrt{d}} - t  \right\}
    &\leq \mathrm{e}^{-dt^2/2}; \\
&\mathbb{P}\left\{  \sigma_{\max}(\bm{\mathsf{\Gamma}}; E) \geq 1 + \frac{w(E)}{\sqrt{d}} + t  \right\}
    &\leq \mathrm{e}^{-dt^2/2}.
\end{aligned}$$*
:::

::: proof
*Proof.* (Sketch) The first inequality is a consequence of Gordon's
minimax theorem and Gaussian concentration. The second inequality is
essentially Chevet's theorem, which follows from Slepian's lemma. See 
for an overview of these ideas. ◻
:::

Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} yields the relations
$$1 - \frac{w(E) + 1}{\sqrt{d}} \lessapprox \sigma_{\min}(\bm{\mathsf{\Gamma}}; E) \leq \sigma_{\max}(\bm{\mathsf{\Gamma}}; E) \lessapprox 1 + \frac{w(E)}{\sqrt{d}}.$$
In other words, the embedding dimension should satisfy
$d > (w(E) + 1)^2$ to ensure that the map $\bm{\mathsf{\Gamma}}$ is
unlikely to annihilate any point in $E$. For this choice of $d$, the
random embedding is unlikely to dilate any point in $E$ by more than a
factor of two.

As a consequence, we have reduced the problem of computing embedding
dimensions for Gaussian maps to the problem of computing Gaussian
widths. In the next two subsections, we work out two important examples.

::: remark
**Remark 19** (Optimality). *The statements in
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} are nearly optimal. One way to see this is to
consider the set $E = \mathbb{S}^{n-1}(\mathbb{R})$, for which the
theorem implies that
$$1 - \sqrt{n/d} \lessapprox \operatorname{\mathbb{E}}\sigma_{\min}(\bm{\mathsf{\Gamma}})
    \leq \operatorname{\mathbb{E}}\sigma_{\max}(\bm{\mathsf{\Gamma}}) \leq 1 + \sqrt{n/d}.$$
The Bai--Yin law [@BS10:Spectral-Analysis Sec. 5.2] confirms that the
first and last inequality are sharp as $n, d \to \infty$ with
$n/d \to \mathrm{const} \in [0,1]$.*

*Moreover, if $E$ is spherically convex (i.e., the intersection of a
convex cone with the unit sphere), then the minimum restricted singular
value satisfies the reverse inequality
$$\mathbb{P}\left\{  \sigma_{\min}(\bm{\mathsf{\Gamma}}; E) \geq 1 - \frac{w(E)}{\sqrt{d}} + t  \right\}
    \leq 2\mathrm{e}^{-dt^2/2}.$$ This result is adapted
from [@TOH14:Gaussian-Minmax].*

*In addition, fifteen years of computational experiments have also shown
that the predictions from
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} are frequently sharp. See  for some examples
and references.*
:::

::: remark
**Remark 20** (History). *The application of Gaussian comparison
theorems in numerical analysis can be traced to work in mathematical
signal processing. used a corollary of Gordon's minimax theorem to study
$\ell_1$ minimization problems. Significant extensions and improvements
of this argument were made by and . seem to have been the first to
recognize that Gordon's minimax theorem can be reversed in the presence
of convexity. A substantial refinement of this observation appeared in .
There is a long series of follow-up works by Babak Hassibi's group that
apply this insight to other problems in signal processing and
communications.*
:::

## Example: Johnson--Lindenstrauss

As a first application of
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"}, let us explain how it implies the classic
dimension reduction result of .

### Overview

Let
$\{ \bm{\mathsf{a}}_1, \dots, \bm{\mathsf{a}}_N \} \subset \mathbb{R}^n$
be a discrete point set. We would like to know when a Gaussian embedding
$\bm{\mathsf{\Gamma}} \in \mathbb{R}^{d \times n}$ approximately
preserves all the pairwise distances between these points:
$$\begin{equation}
 \label{eqn:jl}
1 - \varepsilon\leq \frac{\Vert  \bm{\mathsf{\Gamma}}(\bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j)  \Vert}{ \Vert  \bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j  \Vert }
    \leq 1 + \varepsilon
    \quad\text{for all $i \neq j$.}
\end{equation}$$ The question is how large we must set the embedding
dimension $d$ to achieve distortion $\varepsilon\in (0, 1)$.

### Analysis

We can solve this problem using the machinery described in
Section [8.5](#sec:rsv-gauss){reference-type="ref"
reference="sec:rsv-gauss"}. Consider the set $E$ of normalized chords:
$$E = \left\{ \frac{\bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j}{\Vert \bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j \Vert} : 1 \leq i < j \leq N \right\}.$$
By the definition [\[eqn:rsv\]](#eqn:rsv){reference-type="eqref"
reference="eqn:rsv"} of the restricted singular values,
$$\sigma_{\min}(\bm{\mathsf{\Gamma}}; E) \leq \frac{\Vert  \bm{\mathsf{\Gamma}}(\bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j)  \Vert}{ \Vert  \bm{\mathsf{a}}_i - \bm{\mathsf{a}}_j  \Vert }
    \leq \sigma_{\max}(\bm{\mathsf{\Gamma}}; E).$$ Therefore, we can
invoke Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} to determine how the embedding dimension
controls the distortion.

Let us summarize the argument. First, observe that the Gaussian width of
the set $E$ satisfies
$$w(E) = \operatorname{\mathbb{E}}\max_{\bm{\mathsf{x}} \in E}\ \langle  \bm{\mathsf{g}} , \,  \bm{\mathsf{x}}  \rangle
    \leq \sqrt{2 \log \# E}
    < 2 \sqrt{\log (N/2)}.$$ As a consequence, $$\begin{aligned}
&\mathbb{P}\left\{  \sigma_{\min}(\bm{\mathsf{\Gamma}}; E) \leq 1 - (1 + 2 \sqrt{\log(N/2)})/\sqrt{d} - t  \right\} &\leq \mathrm{e}^{-d t^2 /2}; \\
&\mathbb{P}\left\{  \sigma_{\max}(\bm{\mathsf{\Gamma}}; E) \geq 1 + 2 \sqrt{\log(N/2)}/\sqrt{d} + t  \right\} &\leq \mathrm{e}^{-d t^2 /2}.
\end{aligned}$$ To achieve distortion $\varepsilon$ with high
probability, it is sufficient to choose
$$d \geq 8 \varepsilon^{-2} \log N.$$ In other words, the embedding
dimension only needs to be *logarithmic* in the cardinality $N$ of the
point set. With some additional calculation, we can also extract precise
failure probabilities from this analysis.

### Discussion

Let us close this example with a few comments. In spite of its
prominence, the Johnson--Lindenstrauss embedding lemma is somewhat
impractical. Indeed, since the embedding dimension $d$ is proportional
to $\varepsilon^{-2}$, it is a challenge to achieve small distortions.
Even if we consider the setting where $\varepsilon\approx 1$, the
uniform bound [\[eqn:jl\]](#eqn:jl){reference-type="eqref"
reference="eqn:jl"} may require the embedding dimension to be
prohibitively large.

As a step toward more applicable results, note that the bound on the
*minimum* restricted singular value is more crucial than the bound on
the maximum restricted singular value, because the former ensures that
no two points coalesce after the random embedding. Similarly, it is
often more valuable to preserve the distances between nearby points than
between far-flung points. This observation is the starting point for the
theory of locality sensitive hashing [@GIM99:Similarity-Search].

### History

were concerned with a problem in Banach space geometry, namely the
prospect of extending a Lipschitz function from a finite metric space
into a Hilbert space. The famous lemma from their paper took on a life
of its own when  used it to design efficient approximation algorithms
for some graph problems. used random embeddings to develop new
algorithms for the approximate nearest neighbor problem.
[@AMS99:Space-Complexity; @AGMS02:Tracking-Join] introduced the term
*sketching*, and they showed how to use sketches to track streaming
data. Soon afterwards, and proposed using random embeddings and matrix
sampling for low-rank matrix approximation, bringing these ideas into
the realm of computational linear algebra.

## Example: Subspace embedding {#sec:subspace-embedding}

Next, we consider a question at the heart of randomized linear algebra.
Can we embed an unknown subspace into a lower-dimensional space?

### Overview

Suppose that $L$ is a $k$-dimensional subspace in $\mathbb{R}^n$. We say
that a dimension reduction map $\bm{\mathsf{S}}$ is a *subspace
embedding* for $L$ with distortion $\varepsilon\in (0, 1)$ if
$$\begin{equation}
 \label{eqn:subspace-embedding}
 (1 - \varepsilon) \, \Vert  \bm{\mathsf{x}}  \Vert
    \leq \Vert  \bm{\mathsf{S}}\bm{\mathsf{x}}  \Vert
    \leq (1 + \varepsilon) \, \Vert  \bm{\mathsf{x}}  \Vert
    \quad\text{for every $\bm{\mathsf{x}} \in L$.}
\end{equation}$$ We say that $\bm{\mathsf{S}}$ is *oblivious* if it can
be constructed without knowledge of the subspace $L$, except for its
dimension.

Two questions arise. First, what types of dimension reduction maps yield
(oblivious) subspace embeddings? Second, how large must we choose the
embedding dimension to achieve this outcome?

### Analysis

Gaussian dimension reduction maps yield very good oblivious subspace
embeddings. Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} easily furnishes the justification. Consider
the unit sphere in the subspace:
$E = L \cap \mathbb{S}^{n-1}(\mathbb{R})$. Then construct the Gaussian
dimension reduction map
$\bm{\mathsf{\Gamma}} \in \mathbb{R}^{d \times n}$. In view
of [\[eqn:width-subspace\]](#eqn:width-subspace){reference-type="eqref"
reference="eqn:width-subspace"}, we have $$\begin{aligned}
&\mathbb{P}\left\{  \sigma_{\min}( \bm{\mathsf{\Gamma}}; E ) \leq 1 - (1+\sqrt{k})/\sqrt{d} - t  \right\} &\leq \mathrm{e}^{-d t^2/2}; \\
&\mathbb{P}\left\{  \sigma_{\max}( \bm{\mathsf{\Gamma}}; E ) \geq 1 + \sqrt{k}/\sqrt{d} + t  \right\} &\leq \mathrm{e}^{-d t^2/2}. \\
\end{aligned}$$ As a specific example, we can set the embedding
dimension $d = 2k$ to ensure that
$\Vert \bm{\mathsf{\Gamma}}\bm{\mathsf{x}} \Vert = (1 \pm 0.8) \Vert \bm{\mathsf{x}} \Vert$
simultaneously for all points $\bm{\mathsf{x}} \in L$, except with
probability $\mathrm{e}^{-\mathrm{c} k}$. In some applications of
subspace embeddings, we can even choose the dimension as small as
$d = k + 5$ or $d = k + 10$.

Many theoretical papers on randomized NLA use subspace embeddings as a
primitive for designing algorithms for other linear algebra problems.
For example, Section [10](#sec:overdet-ls){reference-type="ref"
reference="sec:overdet-ls"} describes several ways to use subspace
embeddings to solve overdetermined least-squares problems.

### History

Subspace embeddings were explicitly introduced by ; see also . As work
on randomized NLA accelerated, researchers became interested in more
structured types of subspace embeddings; an early reference is .
Section [9](#sec:dimension-reduction){reference-type="ref"
reference="sec:dimension-reduction"} covers these extensions. See  for a
theoretical perspective on randomized NLA where subspace embeddings take
pride of place.

## Universality of the minimum restricted singular value

We have seen how to apply Gaussian dimension reduction for embedding
discrete point sets and for embedding subspaces.
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} contains precise theoretical results on the
behavior of Gaussian maps in terms of the Gaussian width. To what extent
can we transfer this analysis to other types of random embeddings?

The following theorem [@OT18:Universality-Laws Thm. 9.1] shows that the
bound on the *minimum* restricted singular value in
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} is universal for a large class of random
embeddings. In particular, this class includes sparse random matrices,
whose nonzero entries compose a vanishing proportion of the total.

::: {#thm:universality .theorem}
**Theorem 21** (Universality). *Fix a set
$E \subseteq \mathbb{S}^{n-1}$. Let
$\bm{\mathsf{S}} \in \mathbb{R}^{d \times n}$ be a random matrix whose
entries are independent random variables that satisfy $$%
\operatorname{\mathbb{E}}[ (\bm{\mathsf{S}})_{ij} ] = 0, \quad
\operatorname{\mathbb{E}}[ (\bm{\mathsf{S}})_{ij}^2 ] = d^{-1}, \quad
\operatorname{\mathbb{E}}[ (\bm{\mathsf{S}})_{ij}^5 ] \leq R. %$$ When
$d \leq n$, with high probability, $$\sigma_{\min}(\bm{\mathsf{S}}; E)
    \geq 1 - \frac{w(E)}{\sqrt{d}} - o(\sqrt{n/d}).$$ The constant in
$o(\sqrt{n/d})$ depends only on $R$. A matching lower bound for
$\sigma_{\min}(\bm{\mathsf{S}}; E)$ holds when $E$ is spherically
convex.*
:::

In other words, if $E$ is a moderately large set, the distribution of
the entries of the random map $\bm{\mathsf{S}}$ does not have an impact
on the embedding dimension $d$ sufficient to ensure no point in $E$ is
annihilated.

Theorem [21](#thm:universality){reference-type="ref"
reference="thm:universality"} is confirmed by extensive numerical
experiments [@OT18:Universality-Laws], which demonstrate that dimension
reduction maps with independent, standardized entries have identical
performance for a wide range of examples.

It is perhaps surprising that the bound on the *maximum* restricted
singular value from Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"} is not universal. For some sets $E$, the
quantity $\sigma_{\max}(\bm{\mathsf{S}}; E)$ depends heavily on the
distribution of the entries of $\bm{\mathsf{S}}$.

::: remark
**Remark 22** (Universality for least-squares). *give some asymptotic
universality results for random embeddings in the context of
least-squares problems.*
:::

## Random partial isometries

Last, we consider a variant of Gaussian embedding that is more suitable
when the embedding dimension $d$ is close to the ambient dimension $n$.
In this section, we allow the field $\mathbb{F}$ to be real or complex.

First, suppose that $d \leq n$, and let
$\bm{\mathsf{\Gamma}} \in \mathbb{F}^{d \times n}$ be a Gaussian
embedding. Almost surely, the co-range of $\bm{\mathsf{\Gamma}}$ is a
uniformly random $d$-dimensional subspace of $\mathbb{F}^n$. Construct
an embedding $\bm{\mathsf{S}} \in \mathbb{F}^{d \times n}$ with
orthonormal rows that span the co-range of $\bm{\mathsf{\Gamma}}$, for
example by QR factorization.

Similarly, we can consider a Gaussian embedding
$\bm{\mathsf{\Gamma}} \in \mathbb{F}^{d \times n}$ with $d \geq n$. In
this case, the range of $\bm{\mathsf{\Gamma}}$ is almost surely a
uniformly random $n$-dimensional subspace of $\mathbb{F}^d$. Construct
an embedding $\bm{\mathsf{S}} \in \mathbb{F}^{d \times n}$ with
orthonormal columns that span the range of $\bm{\mathsf{\Gamma}}$, for
example by QR factorization.

In each case, we call $\bm{\mathsf{S}}$ a *random partial isometry*. The
cost of storing a random partial isometry is $O(dn)$, and the cost of
applying it to a vector is $O(dn)$. (We should warn the punctilious
reader that QR factorization of $\bm{\mathsf{\Gamma}}$ may not produce a
matrix $\bm{\mathsf{S}}$ that is Haar-distributed on the Stiefel
manifold. To achieve this guarantee, use the algorithms from .)

When $d \approx n$, random partial isometries are better embeddings than
Gaussian maps (because the nonzero singular values of a partial isometry
are all equal). When $d$ and $n$ are significantly different, the two
models are quite similar to each other.

have established some theoretical results on the embedding behavior of
real partial isometries $(\mathbb{F}= \mathbb{R})$. Unfortunately, the
situation is more complicated than in
Theorem [18](#thm:rsv-gauss){reference-type="ref"
reference="thm:rsv-gauss"}. More relations between Gaussian matrices and
partial isometries follow from the Marcus--Pisier comparison
theorem [@MP81:Random-Fourier]; see also .

