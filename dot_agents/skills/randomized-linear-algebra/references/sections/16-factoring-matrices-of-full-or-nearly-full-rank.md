# Factoring matrices of full or nearly full rank {#sec:full}

So far, we have focused on techniques for computing low-rank
approximations of an input matrix. We will now upgrade to techniques for
computing *full* rank-revealing factorizations such as the
column-pivoted QR (CPQR) decomposition.

Classical deterministic techniques for computing these factorizations
proceed through a sequence of rank-one updates to the matrix, making
them communication-intensive and slow when executed on modern hardware.
Randomization allows the algorithms to be reorganized so that the vast
majority of the arithmetic takes place inside matrix--matrix
multiplications, which greatly accelerates the execution speed.

When applied to an $n \times n$ matrix, most of the algorithms described
in this section have the same $O(n^{3})$ asymptotic complexity as
traditional methods; the objective is to improve the practical speed by
reducing communication costs. However, randomization also allows us to
incorporate Strassen-type techniques to accelerate the matrix
multiplications in a numerically stable manner, attaining an overall
cost of $O(n^{\omega})$ where $\omega$ is the exponent of square
matrix--matrix multiplication.

As well as the CPQR decomposition, we will consider algorithms for
computing factorizations of the form
$\bm{\mathsf{A}} = \bm{\mathsf{U}}\bm{\mathsf{R}}\bm{\mathsf{V}}^{*}$
where $\bm{\mathsf{U}}$ and $\bm{\mathsf{V}}$ are unitary matrices and
$\bm{\mathsf{R}}$ is upper-triangular. Factorizations of this form can
be used for almost any task where either the CPQR or the SVD is
currently used. The additional flexibility allows us to improve on both
the computational speed and on the rank-revealing qualities of the
factorization.

Sections [16.1](#sec:rankrevealing){reference-type="ref"
reference="sec:rankrevealing"}--[16.4](#sec:demmelURV){reference-type="ref"
reference="sec:demmelURV"} introduce the key concepts by describing a
simple algorithm for computing a rank-revealing factorization of a
matrix. This method is both faster than traditional column-pivoted QR
and better at revealing the spectral properties of the matrix. Sections
[16.5](#sec:classicCPQR){reference-type="ref"
reference="sec:classicCPQR"}--[16.7](#sec:randUTV){reference-type="ref"
reference="sec:randUTV"} are more technical; they describe how
randomization can be used to resolve a longstanding challenge of how to
*block* a classical algorithm for computing a column-pivoted QR
decomposition by applying groups of Householder reflectors
simultaneously. They also describe how these ideas can be extended to
the task of computing a URV factorization.

## Rank-revealing factorizations {#sec:rankrevealing}

Before we discuss algorithms, let us first define what we mean when we
say that a factorization is *rank-revealing.* Given an $m\times n$
matrix $\bm{\mathsf{A}}$, we will consider factorizations of the form
$$\begin{equation}
\label{eq:fullfactintro}
\begin{array}{ccccccccccccccccccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{U}} &\bm{\mathsf{R}}&\bm{\mathsf{V}}^{*},\\
m\times n && m\times c & c\times n & n\times n
\end{array}
\end{equation}$$ where $c = \min(m,n)$, where $\bm{\mathsf{U}}$ and
$\bm{\mathsf{V}}$ are orthonormal, and where $\bm{\mathsf{R}}$ is
upper-triangular (or banded upper-triangular). We want the factorization
to reveal the numerical rank of $\bm{\mathsf{A}}$ in the sense that we
obtain a near-optimal approximation of $\bm{\mathsf{A}}$ when we
truncate [\[eq:fullfactintro\]](#eq:fullfactintro){reference-type="eqref"
reference="eq:fullfactintro"} to any level $k$. That is,
$$\begin{equation}
\label{eq:rankreveal}
\|\bm{\mathsf{A}} - \bm{\mathsf{U}}(\colon,1:k)\bm{\mathsf{R}}(1:k,:)\bm{\mathsf{V}}^{*}\|
\approx
\inf\{\|\bm{\mathsf{A}} - \bm{\mathsf{B}}\|\,\colon\,\bm{\mathsf{B}}\mbox{ has rank }k\}
\end{equation}$$ for $k \in \{1,2,\dots,c\}$. The factorization
[\[eq:fullfactintro\]](#eq:fullfactintro){reference-type="eqref"
reference="eq:fullfactintro"} can be viewed either as a generalization
of the SVD (for which $\bm{\mathsf{R}}$ is diagonal) or as a
generalization of the column-pivoted QR factorization (for which
$\bm{\mathsf{V}}$ is a permutation matrix).

A factorization such as
[\[eq:fullfactintro\]](#eq:fullfactintro){reference-type="eqref"
reference="eq:fullfactintro"} that satisfies
[\[eq:rankreveal\]](#eq:rankreveal){reference-type="eqref"
reference="eq:rankreveal"} is very handy. It can be used to solve
ill-conditioned linear systems or least-squares problems; it can be used
for estimating the singular spectrum of $\bm{\mathsf{A}}$ (and all
Schatten $p$-norms); and it provides orthonormal bases for
approximations to the four fundamental subspaces of the matrix. Finally,
it can be used to compute approximate low-rank factorizations
efficiently in situations where the numerical rank of the matrix is not
that much smaller than $m$ or $n$. In contrast, all techniques described
up to now are efficient only when the numerical rank $k$ satisfies
$k \ll \min(m,n)$.

## Blocking of matrix computations {#sec:blocking}

A well-known feature of modern computing is that we can execute
increasingly many floating-point operations because CPUs are gaining
more cores while GPUs and other accelerators are becoming more
affordable and more energy-efficient. In contrast, the cost of
communication (data transfer up and down levels of a memory hierarchy,
among servers, and across networks, etc.) is declining very slowly. As a
result, reducing communication is often the key to accelerating
numerical algorithms in the real world.

In the context of matrix computations, the main reaction to this
development has been to cast linear-algebraic operations as operating on
blocks of the matrix, rather than on individual entries or individual
columns and rows
[@1998_dongarra_blocking; @1993_stewart_blockQR_SVD; @2015_blockQR_SISC].
The objective is to reorganize an algorithm so that the majority of
flops can be executed using highly efficient algorithms for
matrix--matrix multiplication (BLAS3), rather than the slower methods
for matrix--vector multiplications (BLAS2).

Unfortunately, it turns out that classical algorithms for computing
rank-revealing factorizations of matrices are very challenging to block.
Column-pivoted QR proceeds through a sequence of rank-1 updates to the
matrix. The next pivot cannot be found until the previous update has
been applied. Techniques for computing an SVD of a matrix start by
reducing the matrix to bidiagonal form. Then they iterate on the
bidiagonal matrix to drive it towards diagonal form. Both steps are
challenging to block.

To emphasize just how much of a difference blocking makes, let us peek
ahead at the plot in Figure [1](#fig:demmel_URV){reference-type="ref"
reference="fig:demmel_URV"}. This graph shows computational times for
computing certain matrix factorizations versus matrix size on a standard
desktop PC. In particular, look at the times for column-pivoted QR (red
solid line) and for unpivoted QR (red dashed line). The asymptotic flop
counts of these two algorithms are identical in the dominant term. Yet
the unpivoted factorization can easily be blocked, which means that it
executes one order of magnitude faster than the pivoted one.

## The powerURV algorithm {#sec:powerURV}

There is a simple randomized algorithm for computing a rank-revealing
factorization of a matrix that perfectly illustrates the power of
randomization to reduce communication. Our starting point is an
algorithm proposed by . Given an $m\times n$ matrix $\bm{\mathsf{A}}$,
typically with $m \geq n$, it proceeds as follows.

1.  Draw an $n\times n$ matrix $\bm{\mathsf{\Omega}}$ from a standard
    normal distribution.

2.  Perform an unpivoted QR factorization of $\bm{\mathsf{\Omega}}$ so
    that $[\bm{\mathsf{V}},\sim] = \texttt{qr}(\bm{\mathsf{\Omega}})$.

3.  Perform an unpivoted QR factorization of
    $\bm{\mathsf{A}}\bm{\mathsf{V}}$ so that
    $[\bm{\mathsf{U}},\bm{\mathsf{R}}] = \texttt{qr}(\bm{\mathsf{A}}\bm{\mathsf{V}})$.

Observe that the purpose of steps (1) and (2) is simply to generate a
matrix $\bm{\mathsf{V}}$ whose columns serve as a "random" orthonormal
basis. It is easily verified that the matrices $\bm{\mathsf{U}}$,
$\bm{\mathsf{R}}$, and $\bm{\mathsf{V}}$ satisfy $$\begin{equation}
\label{eq:basicURV}
\bm{\mathsf{A}} = \bm{\mathsf{U}}\bm{\mathsf{R}}\bm{\mathsf{V}}^{*}.
\end{equation}$$ The factorization
([\[eq:basicURV\]](#eq:basicURV){reference-type="ref"
reference="eq:basicURV"}) is rank-revealing in theory (see and ), but it
does not reveal the rank particularly well in practice.

The cost to compute
([\[eq:basicURV\]](#eq:basicURV){reference-type="ref"
reference="eq:basicURV"}) is dominated by the cost to perform two
unpivoted QR factorizations, and one matrix--matrix multiplication.
(Simulating the random matrix $\bm{\mathsf{\Omega}}$ requires only
$O(n^{2})$ flops.)

To improve the rank-revealing ability of the factorization, one can
incorporate a small number of power iteration steps
[@2018_martinsson_powerurv], so that step (2) in the recipe gets
modified to $$\begin{equation}
\label{eq:powerURV}
[\bm{\mathsf{V}},\sim] = \texttt{qr}((\bm{\mathsf{A}}^{*}\bm{\mathsf{A}})^{q}\bm{\mathsf{\Omega}}),
\end{equation}$$ where $q$ is a small integer. In practice, $q=1$ or
$q=2$ is often enough to dramatically improve the accuracy of the
computation. Of course, incorporating power iteration increases the cost
of the procedure by adding $2q$ additional matrix--matrix
multiplications. When the singular values of $\bm{\mathsf{A}}$ decay
rapidly, reorthonormalization in between each application of
$\bm{\mathsf{A}}$ is sometimes required to avoid loss of accuracy due to
floating-point arithmetic.

Algorithm [\[alg:powerURV\]](#alg:powerURV){reference-type="ref"
reference="alg:powerURV"} summarizes the techniques introduced in this
::: {#alg:powerURV .algorithm}
**Algorithm. powerURV.**

```
**Input:** Target matrix $\mtx{A} \in \F^{m \times n}$ for $m\geq n$, power parameter $q$
**Output:** Orthonormal matrices $\mtx{U} \in \F^{m\times n}$ and $\mtx{V} = \F^{n\times n}$, and upper triangular $\mtx{R} \in \F^{n\times n}$ such that
$\mtx{A} = \mtx{U}\mtx{R}\mtx{V}^{*}$

function powerURV($\mtx{A}$, $q$):
    Draw a test matrix $\mtx{\Omega} \in \F^{n\times n}$ from a standard normal distribution
    $[\mtx{V},\sim] = qr_econ(\bigl(\mtx{A}^{*}\mtx{A}\bigr)^{q}\mtx{\Omega})$    # Unpivoted QR
    $[\mtx{U},\mtx{R}] = qr_econ(\mtx{A}\mtx{V})$    # Unpivoted QR
```
:::


section. The method is simple, and easy to code. It requires far more
flops than traditional methods for computing rank-revealing
factorizations, yet it is faster in practice. For instance, if $m=n$ and
$q=2$, then powerURV requires $\approx 5n^{3}$ flops versus $0.5\,n^{3}$
flops for CPQR, but Figure [1](#fig:demmel_URV){reference-type="ref"
reference="fig:demmel_URV"} shows that powerURV is still faster. This is
noteworthy, since powerURV with $q=2$ does a *far* better job at
revealing the numerical rank of $\bm{\mathsf{A}}$ than CPQR, as shown in
Figure [2](#fig:powerURVaccuracy){reference-type="ref"
reference="fig:powerURVaccuracy"}. See for details.

:::: algorithm
::: algorithmic
Target matrix $\bm{\mathsf{A}} \in \mathbb{F}^{m \times n}$ for
$m\geq n$, power parameter $q$ Orthonormal matrices
$\bm{\mathsf{U}} \in \mathbb{F}^{m\times n}$ and
$\bm{\mathsf{V}} = \mathbb{F}^{n\times n}$, and upper triangular
$\bm{\mathsf{R}} \in \mathbb{F}^{n\times n}$ such that
$\bm{\mathsf{A}} = \bm{\mathsf{U}}\bm{\mathsf{R}}\bm{\mathsf{V}}^{*}$

Draw a test matrix $\bm{\mathsf{\Omega}} \in \mathbb{F}^{n\times n}$
from a standard normal distribution
$[\bm{\mathsf{V}},\sim] = \texttt{qr\_econ}(\bigl(\bm{\mathsf{A}}^{*}\bm{\mathsf{A}}\bigr)^{q}\bm{\mathsf{\Omega}})$
Unpivoted QR
$[\bm{\mathsf{U}},\bm{\mathsf{R}}] = \texttt{qr\_econ}(\bm{\mathsf{A}}\bm{\mathsf{V}})$
Unpivoted QR
:::
::::

<figure id="fig:demmel_URV">
<div class="picture">
<p>(120,53) (05,02)<span><embed src="Pics/fig_powerURVspeed.pdf"
style="width:115mm" /></span> (00,04) (60,00)<span><span
class="math inline"><em>n</em></span></span></p>
</div>
<figcaption>Computational times required for column-pivoted QR (CPQR)
and unpivoted QR (QR) of an <span
class="math inline"><em>n</em> × <em>n</em></span> real matrix using
MATLAB on an Intel i7-8700k CPU. We see that the unpivoted factorization
is an order of magnitude faster, despite having the identical asymptotic
flop count. The graph also shows the times required for the randomized
rank-revealing factorization described in Section <a
href="#sec:demmelURV" data-reference-type="ref"
data-reference="sec:demmelURV">16.4</a>, executed both on a CPU (solid
lines) and an Nvidia Titan V GPU (dashed lines).</figcaption>
</figure>

<figure id="fig:powerURVaccuracy">
<div class="picture">
<p>(123,45) (-3,04)<span><embed src="Pics/fig_powerURVaccuracy_fast.pdf"
style="width:40mm" /></span> (40,04)<span><embed
src="Pics/fig_powerURVaccuracy_S.pdf" style="width:40mm" /></span>
(83,04)<span><embed src="Pics/fig_powerURVaccuracy_BIE.pdf"
style="width:40mm" /></span> (18,00)<span>(a)</span>
(58,00)<span>(b)</span> (98,00)<span>(c)</span></p>
</div>
<figcaption>The rank-revealing ability of CPQR and powerURV for
different values of the power parameter <span
class="math inline"><em>q</em></span>, as discussed in Section <a
href="#sec:powerURV" data-reference-type="ref"
data-reference="sec:powerURV">16.3</a>. The error <span
class="math inline"><em>e</em><sub><em>k</em></sub> = ∥A − U(:,1 : <em>k</em>)R(1 : <em>k</em>,  : )V<sup>*</sup>∥</span>
(see (<a href="#eq:rankreveal" data-reference-type="ref"
data-reference="eq:rankreveal">[eq:rankreveal]</a>)), is plotted versus
<span class="math inline"><em>k</em></span> for three different matrices
<span class="math inline">A</span> of size <span
class="math inline">400 × 400</span>. The black lines plot the
theoretical minimal values <span
class="math inline"><em>σ</em><sub><em>k</em> + 1</sub></span>. (a) A
matrix whose singular values decay rapidly; we see that all methods
perform well. (b) A matrix whose singular values plateau; we see that
CPQR performs poorly, and so does the randomized method unless powering
is used. (c) A discretized boundary integral operator whose singular
values decay slowly; we again see the high precision of powerURV for
<span class="math inline"><em>q</em> = 1</span> and <span
class="math inline"><em>q</em> = 2</span>.</figcaption>
</figure>

## Computing a rank-revealing factorization of an $n\times n$ matrix in less than $O(n^{3})$ operations {#sec:demmelURV}

The basic version of the randomized URV factorization algorithm
described in Section [16.3](#sec:powerURV){reference-type="ref"
reference="sec:powerURV"} was originally proposed by for purposes
loftier than practical acceleration. Indeed, randomization allows us to
exploit fast matrix--matrix multiplication primitives to design
accelerated algorithms for other NLA problems, such as constructing
rank-revealing factorizations. The main point of this research is that,
whenever the fast matrix--matrix multiplication is stable, the
computation of a rank-revealing factorization is stable too.

Let us be more precise. embark from the observation that there exist
algorithms [^1] for multiplying two $n\times n$ matrices using
$O(n^{\omega})$ flops, where $\omega < 3$. Once such an algorithm is
available, one can stably perform a whole range of other standard matrix
operations at the same asymptotic complexity. The idea is to apply a
divide-and-conquer approach that moves essentially all flops into the
matrix--matrix multiplication. This approach turns out to be relatively
straightforward for decompositions that do not reveal the numerical rank
such as the unpivoted QR factorization. It is harder to implement,
however, for (pivoted) rank-revealing factorizations.

## Classical column-pivoted QR {#sec:classicCPQR}

The powerURV algorithm described in Section
[16.3](#sec:powerURV){reference-type="ref" reference="sec:powerURV"} can
be very effective, but it operates on the whole matrix at once, and it
cannot be used to compute a partial factorization. In the remainder of
this section, we describe algorithms that build a rank-revealing
factorization incrementally. These methods enjoy the property that the
factorization can be halted once a specified tolerance has been met.

We start off this discussion by reviewing a classical (deterministic)
method for computing a column-pivoted QR factorization. This material is
elementary, but the discussion serves to set up a notational framework
that lets us describe the randomized version succinctly in Section
[16.6](#sec:randCPQR){reference-type="ref" reference="sec:randCPQR"}.
Suppose that we are given an $m\times n$ matrix $\bm{\mathsf{A}}$ with
$m\geq n$. We seek a factorization of the form $$\begin{equation}
\label{eq:CPQRbasic}
\begin{array}{cccccccccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{Q}} & \bm{\mathsf{R}} & \bm{\mathsf{\Pi}}^{*},\\
m\times n && m\times n & n\times n & n\times n
\end{array}
\end{equation}$$ where $\bm{\mathsf{Q}}$ has orthonormal columns, where
$\bm{\mathsf{\Pi}}$ is a permutation matrix, and where $\bm{\mathsf{R}}$
is upper-triangular with diagonal elements that decay in magnitude so
that
$|\bm{\mathsf{R}}(1,1)| \geq |\bm{\mathsf{R}}(2,2)| \geq |\bm{\mathsf{R}}(3,3)| \geq \cdots$.
The factors are typically built through a sequence of steps, where
$\bm{\mathsf{A}}$ is driven to upper-triangular form one column at a
time.

To be precise, we start by forming the matrix
$\bm{\mathsf{A}}_{0} = \bm{\mathsf{A}}$. Then we proceed using the
iteration formula
$$\bm{\mathsf{A}}_{j} = \bm{\mathsf{Q}}_{j}^{*}\bm{\mathsf{A}}_{j-1}\bm{\mathsf{\Pi}}_{j},$$
where $\bm{\mathsf{\Pi}}_{j}$ is a permutation matrix that swaps the
$j$th column of $\bm{\mathsf{A}}_{j-1}$ with the column in
$\bm{\mathsf{A}}_{j-1}(:,j:n)$ that has the largest magnitude, and where
$\bm{\mathsf{Q}}_{j}$ is a Householder reflector that zeros out all
elements below the diagonal in the $j$th column of
$\bm{\mathsf{A}}_{j-1}\bm{\mathsf{\Pi}}_{j}$; see Figure
[3](#fig:classicCPQR){reference-type="ref" reference="fig:classicCPQR"}.
Once the process concludes, the relation
[\[eq:CPQRbasic\]](#eq:CPQRbasic){reference-type="eqref"
reference="eq:CPQRbasic"} holds for
$$\bm{\mathsf{Q}}   = \bm{\mathsf{Q}}_{n}\bm{\mathsf{Q}}_{n-1}\bm{\mathsf{Q}}_{n-2}\dots\bm{\mathsf{Q}}_{1},\quad
\bm{\mathsf{R}}   = \bm{\mathsf{A}}_{n},\quad
\bm{\mathsf{\Pi}} = \bm{\mathsf{\Pi}}_{n}\bm{\mathsf{\Pi}}_{n-1}\bm{\mathsf{\Pi}}_{n-2}\dots\bm{\mathsf{\Pi}}_{1}.$$

This algorithm is well understood, and it is ubiquitous in numerical
computations. For exotic matrices, it can produce factorizations that
are quite far from optimal [@1966_kahan_NLA Sec. 5], but it typically
works very well for many tasks. For instance, it serves for revealing
the numerical rank of a matrix or for solving an ill-conditioned linear
system. However, a serious drawback to this algorithm is that it
fundamentally consists of a sequence of $n-1$ steps (or $n$ steps if
$m > n$), where a large part of the matrix is updated in each step.

<figure id="fig:classicCPQR">
<table>
<tbody>
<tr>
<td style="text-align: center;"><embed src="Pics/fig_qr0.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_qr1.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_qr2.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_qr3.pdf"
style="width:20mm" /></td>
</tr>
<tr>
<td style="text-align: center;"><span
class="math inline">A<sub>0</sub> = A</span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>1</sub> = Q<sub>1</sub><sup>*</sup>A<sub>0</sub>Π<sub>1</sub></span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>2</sub> = Q<sub>2</sub><sup>*</sup>A<sub>1</sub>Π<sub>2</sub></span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>3</sub> = Q<sub>3</sub><sup>*</sup>A<sub>2</sub>Π<sub>3</sub></span></td>
</tr>
</tbody>
</table>
<figcaption>The figure shows the sparsity pattern of a <span
class="math inline">4 × 3</span> matrix <span
class="math inline">A</span> as it is driven to upper-triangular form in
the column-pivoted QR factorization algorithm described in Section <a
href="#sec:classicCPQR" data-reference-type="ref"
data-reference="sec:classicCPQR">16.5</a>. The process takes three steps
in this case, and step <span class="math inline"><em>j</em></span>
involves the application of a permutation matrix <span
class="math inline">Π<sub><em>j</em></sub></span> from the right, and by
a Householder reflector <span
class="math inline">Q<sub><em>j</em></sub></span> from the
left.</figcaption>
</figure>

## A randomized algorithm for computing a CPQR decomposition {#sec:randCPQR}

Our next objective is to recast the algorithm for computing a CPQR
decomposition that was introduced in Section
[16.5](#sec:classicCPQR){reference-type="ref"
reference="sec:classicCPQR"} so that it works with "panels" of $b$
contiguous columns, as shown in Figure
[4](#fig:randCPQR){reference-type="ref" reference="fig:randCPQR"}. The
difficulty is to find a set of $b$ pivot vectors without updating the
matrix between each selection. Fortuitously, the randomized algorithm
for interpolatory decomposition
(Section [13.4](#sec:randID){reference-type="ref"
reference="sec:randID"}) is well adapted for this task. Indeed, a set of
$b$ columns that forms a good basis for the column space also forms a
good set of pivot vectors.

To be specific, let us describe how to pick the first group of $b$ pivot
columns for an $m\times n$ matrix $\bm{\mathsf{A}}$. Adapting the ideas
in Section [13.4](#sec:randID){reference-type="ref"
reference="sec:randID"}, we draw a Gaussian random matrix
$\bm{\mathsf{\Omega}}$ of size $(b+p) \times m$, where $p$ is a small
oversampling parameter. We form a sample matrix
$$\begin{array}{ccccccccccc}
\bm{\mathsf{Y}} &=& \bm{\mathsf{\Omega}} & \bm{\mathsf{A}}, \\
(b+p)\times n && (b+p) \times m & m\times n
\end{array},$$ and then we execute $b$ steps of column-pivoted QR on the
matrix $\bm{\mathsf{Y}}$ (either Householder or Gram--Schmidt is fine
for this step). The resulting $b$ pivot columns turn out to be good
pivot columns for $\bm{\mathsf{A}}$ as well. Once these $b$ columns have
been moved to the front of $\bm{\mathsf{A}}$, we perform a local CPQR
factorization of this panel. We update the remaining $n-b$ columns using
the computed Householder reflectors.

We could then proceed using exactly the same method a second time: draw
a $(b+p) \times (m-b)$ Gaussian random matrix $\bm{\mathsf{\Omega}}$,
form a $(b+p) \times (n-b)$ sample matrix $\bm{\mathsf{Y}}$, perform
classical CPQR on $\bm{\mathsf{Y}}$, and so on. However, there is a
shortcut. We can update the sample matrix that was used in the first
step, which renders the overhead cost induced by randomization almost
negligible [@2017_blockQR_ming_article Sec. 4].

<figure id="fig:randCPQR">
<table>
<tbody>
<tr>
<td style="text-align: center;"><embed src="Pics/fig_blockqr0.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_blockqr1.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_blockqr2.pdf"
style="width:20mm" /></td>
<td style="text-align: center;"><embed src="Pics/fig_blockqr3.pdf"
style="width:20mm" /></td>
</tr>
<tr>
<td style="text-align: center;"><span
class="math inline">A<sub>0</sub> = A</span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>1</sub> = Q<sub>1</sub><sup>*</sup>A<sub>0</sub>Π<sub>1</sub></span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>2</sub> = Q<sub>2</sub><sup>*</sup>A<sub>1</sub>Π<sub>2</sub></span></td>
<td style="text-align: center;"><span
class="math inline">A<sub>3</sub> = Q<sub>3</sub><sup>*</sup>A<sub>2</sub>Π<sub>3</sub></span></td>
</tr>
</tbody>
</table>
<figcaption>The sparsity pattern of a matrix <span
class="math inline">A</span> consisting of <span
class="math inline">4 × 3</span> blocks, each of size <span
class="math inline">3 × 3</span>, as it undergoes the blocked version of
the Householder QR algorithm described in Section <a
href="#sec:randCPQR" data-reference-type="ref"
data-reference="sec:randCPQR">16.6</a>. Each matrix <span
class="math inline">Q<sub><em>j</em></sub></span> is a product of three
Householder reflectors. The difficulty in building an algorithm of this
type is to find groups of pivot vectors <em>before</em> applying the
corresponding Householder reflectors.</figcaption>
</figure>

<figure id="fig:HQRRPspeedup" data-latex-placement="b">
<div class="picture">
<p>(100,54) (07,05)<span><embed src="Pics/fig_HQRRP_comparison_new.pdf"
style="width:87mm" /></span> (40,00)<span><em>Matrix size <span
class="math inline"><em>n</em></span>.</em></span> (00,17)</p>
</div>
<figcaption>Speedup of the randomized algorithm for computing a
column-pivoted QR decomposition described in Section <a
href="#sec:randCPQR" data-reference-type="ref"
data-reference="sec:randCPQR">16.6</a>, relative to LAPACK’s faster
routine (dgeqp3) as implemented in the Intel MKL library (version
11.2.3), running on an Intel Xeon E5-2695 v3 processor.</figcaption>
</figure>

Extensive numerical work has demonstrated dramatic acceleration over
deterministic algorithms. Figure
[5](#fig:HQRRPspeedup){reference-type="ref"
reference="fig:HQRRPspeedup"} draws on data from that illustrates the
acceleration over a state-of-the-art software implementation of the
classical CPQR method. Computer experiments also show that the
randomized scheme chooses pivot columns whose quality is almost
indistinguishable from those chosen by traditional pivoting, in the
sense that the relation
[\[eq:rankreveal\]](#eq:rankreveal){reference-type="eqref"
reference="eq:rankreveal"} holds to about the same accuracy. (However,
the diagonal entries of $\bm{\mathsf{R}}$ do not strictly decay in
magnitude across the block boundaries.)

To understand the behavior of the algorithm, it is helpful to think
about two extreme cases. In the first, suppose that the singular values
of $\bm{\mathsf{A}}$ decay very rapidly. Here, the analysis in Section
[11](#sec:random-rangefinder){reference-type="ref"
reference="sec:random-rangefinder"} can be modified to show that, for
any $j$, the first $j$ pivot columns chosen by the randomized algorithm
is likely to span the column space nearly as well as the optimal set of
$j$ columns. Therefore, they are excellent pivot vectors. At the other
extreme, suppose that the singular values of $\bm{\mathsf{A}}$ hardly
decay at all. In this case, the randomized method may pick a completely
different set of pivot vectors than the deterministic method, but this
outcome is unproblematic because we can take any group of columns as
pivot vectors. Of course, the interesting cases are intermediate between
these two extremes. It turns out that the randomized methods work well
regardless of how rapidly the singular values decay. For a detailed
analysis, see
[@2017_blockQR_ming_article; @2017_gu_langou; @2015_gu_melgaard; @2018_gu_flipflop].

::: remark
**Remark 51** (History). *Finding a blocked version of the CPQR
algorithm described in Section
[16.5](#sec:classicCPQR){reference-type="ref"
reference="sec:classicCPQR"} has remained an open challenge in NLA for
some time
[@1998_quintana_panelRRQR1; @2015_demmel_communication_avoiding]. The
randomized technique described in this section was introduced in , while
the updating technique was described in . For full details, see
[@2015_blockQR_SISC] and [@2017_blockQR_ming_article].*
:::

## A randomized algorithm for computing a URV decomposition {#sec:randUTV}

In this section, we describe an incremental randomized algorithm for
computing the URV factorization
([\[eq:fullfactintro\]](#eq:fullfactintro){reference-type="ref"
reference="eq:fullfactintro"}). Let us recall that for
$\bm{\mathsf{A}} \in \mathbb{F}^{m\times n}$, with $m \geq n$, this
factorization takes the form $$\begin{equation}
\label{eq:UTVrepeat}
\begin{array}{cccccccccc}
\bm{\mathsf{A}} &=& \bm{\mathsf{U}} & \bm{\mathsf{R}} & \bm{\mathsf{V}}^{*},\\
m\times n && m\times n & n\times n & n\times n
\end{array}
\end{equation}$$ where $\bm{\mathsf{U}}$ and $\bm{\mathsf{V}}$ have
orthonormal columns and where $\bm{\mathsf{R}}$ is upper-triangular.

The algorithm we describe is blocked, and executes efficiently on modern
computing platforms. It is similar in speed to the randomized CPQR
described in Section [16.6](#sec:randCPQR){reference-type="ref"
reference="sec:randCPQR"}, and it shares the advantage that the
decomposition is built incrementally so that the process can be stopped
once a requested accuracy has been met. However, the URV factorization
offers compelling advantages: (1) It is almost as good at revealing the
numerical rank as the SVD (unlike the CPQR). (2) The URV factorization
provides us with orthonormal basis vectors for both the column and the
row spaces. (3) The off-diagonal entries of $\bm{\mathsf{R}}$ are very
small in magnitude. (4) The diagonal entries of $\bm{\mathsf{R}}$ form
excellent approximations to the singular values of $\bm{\mathsf{A}}$.

The randomized algorithm for computing a URV factorization we present
follows the same algorithmic template as the randomized CPQR described
in Section [16.6](#sec:randCPQR){reference-type="ref"
reference="sec:randCPQR"}. It drives $\bm{\mathsf{A}}$ to
upper-triangular form one block at a time, but it replaces the
permutation matrices $\bm{\mathsf{\Pi}}_{j}$ in the CPQR with general
unitary matrices $\bm{\mathsf{V}}_{j}$. Using this increased freedom, we
can obtain a factorization where all diagonal blocks are themselves
diagonal matrices and all off-diagonal elements have small magnitude.
Figure [6](#fig:blockUTV){reference-type="ref" reference="fig:blockUTV"}
summarizes the process.

<figure id="fig:blockUTV">
<div class="picture">
<p>(125,40) (000,04)<span><embed src="Pics/fig_blockUTV0.pdf"
style="width:24mm" /></span> (032,04)<span><embed
src="Pics/fig_blockUTV1.pdf" style="width:24mm" /></span>
(064,04)<span><embed src="Pics/fig_blockUTV2.pdf"
style="width:24mm" /></span> (096,04)<span><embed
src="Pics/fig_blockUTV3.pdf" style="width:24mm" /></span>
(007,00)<span><span class="math inline">A<sub>0</sub> = A</span></span>
(035,00)<span><span
class="math inline">A<sub>1</sub> = U<sub>1</sub><sup>*</sup>A<sub>0</sub>V<sub>1</sub></span></span>
(067,00)<span><span
class="math inline">A<sub>2</sub> = U<sub>2</sub><sup>*</sup>A<sub>1</sub>V<sub>2</sub></span></span>
(099,00)<span><span
class="math inline">A<sub>3</sub> = U<sub>3</sub><sup>*</sup>A<sub>2</sub>V<sub>3</sub></span></span></p>
</div>
<figcaption>Sparsity pattern of a matrix being driven to
upper-triangular form in the randomized URV factorization algorithm
described in Section <a href="#sec:randUTV" data-reference-type="ref"
data-reference="sec:randUTV">16.7</a>; see Figure <a
href="#fig:randCPQR" data-reference-type="ref"
data-reference="fig:randCPQR">4</a>. The matrices <span
class="math inline">U<sub><em>i</em></sub></span> and <span
class="math inline">V<sub><em>i</em></sub></span> are now more general
unitary transforms (consisting in the bulk of Householder reflectors).
The entries shown as gray are nonzero, but are typically very small in
magnitude.</figcaption>
</figure>

To provide details on how the algorithm works, suppose that we are given
an $m\times n$ matrix $\bm{\mathsf{A}}$ and a block size $b$. In the
first step of the process, our objective is then to build unitary
matrices $\bm{\mathsf{U}}_{1}$ and $\bm{\mathsf{V}}_{1}$ such that
$$\bm{\mathsf{A}} = \bm{\mathsf{U}}_{1}\bm{\mathsf{A}}_{1}\bm{\mathsf{V}}_{1}^{*},$$
where $\bm{\mathsf{A}}_{1}$ has the block structure
$$\bm{\mathsf{A}}_{1} =
\left[\begin{array}{cc}
\bm{\mathsf{A}}_{1,11} & \bm{\mathsf{A}}_{1,12} \\
\bm{\mathsf{0}}        & \bm{\mathsf{A}}_{1,22}
\end{array}\right] =
\raisebox{-9mm}{\includegraphics[height=18mm]{Pics/fig_UTV_step1.png}},$$
so that the $b\times b$ matrix $\bm{\mathsf{A}}_{1,11}$ is diagonal and
the entries of $\bm{\mathsf{A}}_{1,12}$ are small in magnitude. To build
$\bm{\mathsf{V}}_{1}$, we use the randomized power iteration described
in Section [11.6](#sec:rrf-subspace){reference-type="ref"
reference="sec:rrf-subspace"} to find a basis that approximately spans
the same space as the top $b$ right singular vectors of
$\bm{\mathsf{A}}$. To be precise, we form the sample matrix
$$\bm{\mathsf{Y}} = \bm{\mathsf{\Omega}}\bm{\mathsf{A}}\bigl(\bm{\mathsf{A}}^{*}\bm{\mathsf{A}}\bigr)^{q},$$
where $\bm{\mathsf{\Omega}}$ is a Gaussian random matrix of size
$b \times m$ and where $q$ is a parameter indicating the number of steps
of power iteration taken. We then perform an unpivoted QR factorization
of the *rows* of $\bm{\mathsf{Y}}$ to form a matrix
$\widetilde{\bm{\mathsf{V}}}$ whose first $b$ columns form an
orthonormal basis for the column space of $\bm{\mathsf{Y}}$. We then
apply $\widetilde{\bm{\mathsf{V}}}$ from the right to form the matrix
$\bm{\mathsf{A}}\widetilde{\bm{\mathsf{V}}}$, and we perform an
unpivoted QR factorization of the first $b$ columns of
$\bm{\mathsf{A}}\widetilde{\bm{\mathsf{V}}}$. This results in a new
matrix
$$\widetilde{\bm{\mathsf{A}}} = \bigl(\widetilde{\bm{\mathsf{U}}}\bigr)^{*}\,\bm{\mathsf{A}}\,\widetilde{\bm{\mathsf{V}}}$$
that has the block structure $$\widetilde{\bm{\mathsf{A}}} =
\left[\begin{array}{cc}
\widetilde{\bm{\mathsf{A}}}_{1,11} & \widetilde{\bm{\mathsf{A}}}_{1,12} \\
\bm{\mathsf{0}}                & \widetilde{\bm{\mathsf{A}}}_{1,22}
\end{array}\right] =
\raisebox{-9mm}{\includegraphics[height=18mm]{Pics/fig_UTV_step0.png}}.$$
The top left $b\times b$ block $\widetilde{\bm{\mathsf{A}}}_{1,11}$ is
upper-triangular, and the bottom left block is zero. The entries of
$\widetilde{\bm{\mathsf{A}}}_{1,12}$ are typically small in magnitude.
Next, we compute a full SVD of the block
$\widetilde{\bm{\mathsf{A}}}_{1,11}$:
$$\widetilde{\bm{\mathsf{A}}}_{1,11} = \widehat{\bm{\mathsf{U}}}\bm{\mathsf{D}}_{11}\widehat{\bm{\mathsf{V}}}^{*}.$$
This step is inexpensive because $\widetilde{\bm{\mathsf{A}}}_{1,11}$
has size $b\times b$, where $b$ is small. As a final step, we form the
transformation matrices $$\bm{\mathsf{U}}_{1} =
\widetilde{\bm{\mathsf{U}}}
\left[\begin{array}{cc}
\widehat{\bm{\mathsf{U}}} & \bm{\mathsf{0}} \\
\bm{\mathsf{0}} & \bm{\mathsf{I}}_{m-b}
\end{array}\right],
\qquad\mbox{and}\qquad
\bm{\mathsf{V}}_{1} =
\widetilde{\bm{\mathsf{V}}}
\left[\begin{array}{cc}
\widehat{\bm{\mathsf{V}}} & \bm{\mathsf{0}} \\
\bm{\mathsf{0}} & \bm{\mathsf{I}}_{n-b}
\end{array}\right]$$ and set
$$\bm{\mathsf{A}}_{1} = \bm{\mathsf{U}}_{1}^{*}\bm{\mathsf{A}}\bm{\mathsf{V}}_{1}.$$
The result of this process is that the diagonal entries of
$\bm{\mathsf{D}}_{11}$ typically give accurate approximations to the
first $b$ singular values of $\bm{\mathsf{A}}$, and the "remainder"
matrix $\bm{\mathsf{A}}_{1,22}$ has spectral norm that is similar to
$\sigma_{b+1}$. Thus,
$$\|\bm{\mathsf{A}}_{1,22}\| \approx \inf\{\|\bm{\mathsf{A}} - \bm{\mathsf{B}}\|\,\colon\,\bm{\mathsf{B}}\mbox{ has rank }b\}.$$
This process corresponds to the first step in Figure
[6](#fig:blockUTV){reference-type="ref" reference="fig:blockUTV"}. The
succeeding iterations execute the same process, at each step working on
the remaining lower-right part of the matrix that has not yet been
driven to upper-triangular form. We refer to and for details.

From a theoretical point of view, the first step of the URV
factorization described is well understood since it is mathematically
equivalent to the randomized power iteration described in Section
[11.6](#sec:rrf-subspace){reference-type="ref"
reference="sec:rrf-subspace"}. We observe that the first $b$ columns of
$\bm{\mathsf{U}}_{1}$ do a better job of spanning the column space of
$\bm{\mathsf{A}}$ than the first $b$ columns of $\bm{\mathsf{V}}_{1}$ do
for spanning the row space; the reason for this asymmetry is that by
forming the product $\bm{\mathsf{A}}\bm{\mathsf{V}}_{1}$, we in effect
perform an additional step of the power iteration
[@2018_martinsson_powerurv Sec. 6].

An important feature of the method described in this section is that it
is incremental, and it can be halted once a given computational
tolerance has been met. This feature has been a key competitive
advantage of the column-pivoted QR decomposition, and it is often cited
as the motivation for using CPQR. The method described in this section
has almost all the advantages of the randomized Householder CPQR
factorization (it is blocked, it is incremental, and it executes very
fast in practice), while resulting in a factorization that is far closer
to optimal in revealing the rank.

::: remark
**Remark 52** (Related work). *The idea of loosening the requirements on
the factors in a rank-revealing factorization and searching for a
decomposition such as
([\[eq:UTVrepeat\]](#eq:UTVrepeat){reference-type="ref"
reference="eq:UTVrepeat"}) is well explored in the literature
[@1999_hansen_UTVtools; @1994_stewart_UTV; @1995_elden_downdate_UTV] and
[@1998_stewart_volume1]. Deterministic techniques for computing the URV
decomposition are described in
[@1999_hansen_UTVtools; @1999_stewart_QLP]; these algorithms combine
some of the appealing qualities of the SVD (high accuracy in revealing
the rank) with some of the appealing qualities of CPQR (the possibility
of halting the execution once a requested tolerance has been met).
However, they were not blocked, and therefore they were subject to the
same liabilities as deterministic algorithms for computing the SVD and
the CPQR. A more recent use of randomization in this context is
described in [@2018_gu_flipflop].*
:::

