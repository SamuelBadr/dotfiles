# High-accuracy approximation of kernel matrices {#sec:rankstructured}

In this section, we continue the discussion of kernel matrices that we
started in Section [19](#sec:kernel){reference-type="ref"
reference="sec:kernel"}, but we now consider the high-accuracy regime.
In particular, given a kernel matrix $\bm{\mathsf{K}}$, we seek an
approximation $\bm{\mathsf{K}}_{\mathrm{approx}}$ for which
$\|\bm{\mathsf{K}} - \bm{\mathsf{K}}_{\mathrm{approx}}\|$ is small, say
of relative accuracy $10^{-3}$ or $10^{-6}$. This objective was not
realistic for the applications discussed in Section
[19](#sec:kernel){reference-type="ref" reference="sec:kernel"}, but it
can be achieved in situations where we have access to fast techniques
for evaluating the matrix-vector product
$\bm{\mathsf{x}} \mapsto \bm{\mathsf{K}}\bm{\mathsf{x}}$ (and also
$\bm{\mathsf{x}} \mapsto \bm{\mathsf{K}}^{*}\bm{\mathsf{x}}$ when
$\bm{\mathsf{K}}$ is not self-adjoint). The algorithms that we describe
will build a data sparse approximation to $\bm{\mathsf{K}}$ by using
information in samples such as $\bm{\mathsf{K}}\bm{\mathsf{x}}$ and
$\bm{\mathsf{K}}^{*}\bm{\mathsf{x}}$ for random vectors
$\bm{\mathsf{x}}$. These techniques are particularly well suited to
problems that arise in modeling physical phenomena such as
electromagnetic scattering, or the deformation of solid bodies; we
describe how the fast matrix-vector application we need can be realized
in Section [20.3](#sec:cpkernel_application){reference-type="ref"
reference="sec:cpkernel_application"}.

As in Section [19](#sec:kernel){reference-type="ref"
reference="sec:kernel"}, we say that a matrix
$\bm{\mathsf{K}} \in \mathbb{C}^{n\times n}$ is a *kernel matrix* if its
entries are given by a formula such as $$\begin{equation}
\label{eq:Aij_kernel}
\bm{\mathsf{K}}(i,j) = k(\bm{\mathsf{x}}_{i},\bm{\mathsf{x}}_{j}),
\end{equation}$$ where $\{\bm{\mathsf{x}}_{i}\}_{i=1}^{n}$ is a set of
points in $\mathbb{R}^{d}$, and where
$k\,\colon\,\mathbb{R}^{d} \times \mathbb{R}^{d} \rightarrow \mathbb{C}$
is a kernel function. The kernel matrices that we consider tend to have
singular values that decay slowly or not at all, which rules out the
possibility that $\bm{\mathsf{K}}_{\mathrm{approx}}$ could have low
rank. Instead, we build an approximation
$\bm{\mathsf{K}}_{\mathrm{approx}}$ that is tessellated into $O(n)$
blocks in such as way that each off-diagonal block has low rank. Figure
[11](#fig:tree){reference-type="ref" reference="fig:tree"}(b) shows a
representative tessellation pattern. We say that a matrix of this type
is a *rank-structured hierarchical matrix.*

The purpose of determining a rank-structured approximation to an
operator for which we already have fast matrix-vector multiplication
techniques available is that the new representation can be used to
rapidly execute a whole range of linear algebraic operations: matrix
inversion, LU factorization, and even full spectral decompositions in
certain cases.

This section is structured to provide a high-level description of the
core ideas in Sections [20.1](#sec:sepvar){reference-type="ref"
reference="sec:sepvar"} --
[20.3](#sec:cpkernel_application){reference-type="ref"
reference="sec:cpkernel_application"}. Additional details follow in
Sections [20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"}
-- [20.9](#sec:ASKIT){reference-type="ref" reference="sec:ASKIT"}.

## Separation of variables and low-rank approximation {#sec:sepvar}

The reason that many kernel matrices can be tessellated into blocks that
have low numerical rank is that the function
$(\bm{\mathsf{x}},\bm{\mathsf{y}}) \mapsto k(\bm{\mathsf{x}},\bm{\mathsf{y}})$
is typically smooth as long as $\bm{\mathsf{x}}$ and $\bm{\mathsf{y}}$
are not close. To illustrate the connection, let us consider a
computational domain $D$ that holds a set of points
$\{\bm{\mathsf{x}}_{i}\}_{i=1}^{n} \subset \mathbb{R}^{2}$, as shown in
Figure [10](#fig:sources_targets){reference-type="ref"
reference="fig:sources_targets"}. Suppose further that $D_{\mathrm{s}}$
and $D_{\mathrm{t}}$ are two subdomains of $D$ that are located a bit
apart from each other, as shown in the figure. When the kernel function
$k$ is smooth, we can typically approximate it to high accuracy through
an approximate separation of variables of the form $$\begin{equation}
\label{eq:kernel_sepvar}
k(\bm{\mathsf{x}},\bm{\mathsf{y}}) \approx \sum_{p=1}^{P}b_{p}(\bm{\mathsf{x}})\,c_{p}(\bm{\mathsf{y}}),\qquad\bm{\mathsf{x}} \in D_{\mathrm{t}},\ \bm{\mathsf{y}} \in D_{\mathrm{s}}.
\end{equation}$$ Let $I_{\mathrm{s}}$ and $I_{\mathrm{t}}$ denote two
index vectors that identify the points located in $D_{\mathrm{s}}$ and
$D_{\mathrm{t}}$, respectively (so that, for example,
$i\in I_{\mathrm{s}}$ if and only if
$\bm{\mathsf{x}}_{i} \in D_{\mathrm{s}}$). Then combining the formula
([\[eq:Aij_kernel\]](#eq:Aij_kernel){reference-type="ref"
reference="eq:Aij_kernel"}) with the separation of variables
([\[eq:kernel_sepvar\]](#eq:kernel_sepvar){reference-type="ref"
reference="eq:kernel_sepvar"}), we get $$\begin{equation}
\label{eq:A_sepvar}
\bm{\mathsf{K}}(i,j) \approx \sum_{p=1}^{P}b_{p}(\bm{\mathsf{x}}_{i})\,c_{p}(\bm{\mathsf{x}}_{j}),\qquad i \in I_{\mathrm{t}},\ j \in I_{\mathrm{s}}.
\end{equation}$$ Equation
([\[eq:A_sepvar\]](#eq:A_sepvar){reference-type="ref"
reference="eq:A_sepvar"}) is exactly the low-rank approximation to the
block $\bm{\mathsf{K}}(I_{\mathrm{t}},I_{\mathrm{s}})$ that we seek. To
be precise, ([\[eq:A_sepvar\]](#eq:A_sepvar){reference-type="ref"
reference="eq:A_sepvar"}) can be written as
$$\bm{\mathsf{K}}(I_{\mathrm{t}},I_{\mathrm{s}}) \approx \bm{\mathsf{B}}\bm{\mathsf{C}},$$
where $\bm{\mathsf{B}}$ and $\bm{\mathsf{C}}$ are defined via
$\bm{\mathsf{B}}(i,p) = b_{p}(\bm{\mathsf{x}}_{i})$ and
$\bm{\mathsf{C}}(p,j) = c_{p}(\bm{\mathsf{x}}_{j})$.

A separation of variables such as
([\[eq:kernel_sepvar\]](#eq:kernel_sepvar){reference-type="ref"
reference="eq:kernel_sepvar"}) is sometimes provided through analytic
knowledge of the kernel function, as illustrated in Example
[68](#example:multipole){reference-type="ref"
reference="example:multipole"}. Perhaps more typically, all we know is
that such a formula *should* in principle exist, for instance because we
know that the matrix approximates a singular integral operator for which
a Calderón-Zygmund decomposition must exist. It is then the task of the
randomized algorithm to explicitly build the factors $\bm{\mathsf{B}}$
and $\bm{\mathsf{C}}$, given a computational tolerance.

<figure id="fig:sources_targets">
<div class="picture">
<p>(123,67) (05,00)<span><embed src="Pics/fig_sources_targets.pdf"
style="width:65mm" /></span> (00,03)<span><span
class="math inline"><em>D</em></span></span> (26,15)<span> <span
class="math inline"><strong>D</strong><sub>s</sub></span></span>
(57,34)<span><span
class="math inline"><strong>D</strong><sub>t</sub></span></span>
(72,34)</p>
<div class="minipage">
<p>A box <span class="math inline"><em>D</em></span> holding points
<span
class="math inline">{x<sub><em>i</em></sub>}<sub><em>i</em> = 1</sub><sup><em>n</em></sup></span>
(the blue, red, and gray dots) that define a kernel matrix <span
class="math inline">K(<em>i</em>, <em>j</em>) = <em>k</em>(x<sub><em>i</em></sub>, x<sub><em>j</em></sub>)</span>.
The regions <span class="math inline"><em>D</em><sub>s</sub></span> (the
red box) and <span class="math inline"><em>D</em><sub>t</sub></span>
(the blue box) are separated enough that <span
class="math inline"><em>k</em>(x, y)</span> is smooth when <span
class="math inline">x ∈ <em>D</em><sub>t</sub></span> and <span
class="math inline">y ∈ <em>D</em><sub>s</sub></span>. In consequence,
<span
class="math inline">K(<em>I</em><sub>t</sub>, <em>I</em><sub>s</sub>)</span>
has low numerical rank, where <span
class="math inline"><em>I</em><sub>s</sub></span> identifies the red
points and <span class="math inline"><em>I</em><sub>t</sub></span> the
blue.</p>
</div>
</div>
<figcaption>The geometry discussed in Section <a
href="#eq:kernel_sepvar" data-reference-type="ref"
data-reference="eq:kernel_sepvar">[eq:kernel_sepvar]</a>.</figcaption>
</figure>

::: {#example:multipole .example}
**Example 68** (Laplace kernel in two dimensions). *A standard example
of a kernel matrix in mathematical physics is the matrix
$\bm{\mathsf{K}}$ that maps a vector of electric source strengths
$\bm{\mathsf{q}} = (q_{j})_{j=1}^{n}$ to a vector of potentials
$\bm{\mathsf{u}} = (u_{i})_{i=1}^{n}$. When the set
$\{\bm{\mathsf{x}}_{i}\}_{i=1}^{n} \subset \mathbb{R}^{2}$ identifies
the locations of both the source and the target points,
$\bm{\mathsf{K}}$ takes the form
([\[eq:Aij_kernel\]](#eq:Aij_kernel){reference-type="ref"
reference="eq:Aij_kernel"}), for $$k(\bm{\mathsf{x}},\bm{\mathsf{y}}) =
\left\{\begin{aligned}
\log|\bm{\mathsf{x}}-\bm{\mathsf{y}}|,&\quad\mbox{when}\ \bm{\mathsf{x}}\neq \bm{\mathsf{y}},\\
0,&\quad\mbox{when}\ \bm{\mathsf{x}} = \bm{\mathsf{y}}.
\end{aligned}\right.$$ There is a well-known result from potential
theory that provides the required separation of variables. Expressing
$\bm{\mathsf{x}}$ and $\bm{\mathsf{y}}$ in polar coordinates with
respect to some expansion center $\bm{\mathsf{c}}$, so that
$\bm{\mathsf{x}} - \bm{\mathsf{c}} = r(\cos\theta,\,\sin\theta)$ and
$\bm{\mathsf{y}} - \bm{\mathsf{c}} = r'(\cos\theta',\,\sin\theta')$, the
separation of variables (known as a *multipole expansion*)
$$\begin{equation}
\label{eq:mpole}
k(\bm{\mathsf{x}},\bm{\mathsf{y}}) = \log(r) - \sum_{p=1}^{\infty}\frac{1}{p}\left(\frac{r'}{r}\right)^{p}\bigl(\cos(p\theta)\cos(p\theta') + \sin(p\theta)\sin(p\theta')\bigr)
\end{equation}$$ is valid whenever $r' < r$.*
:::

## Rank-structured matrices and randomized compression {#sec:bestiary}

The observation that the off-diagonal blocks of a kernel matrix often
have low numerical rank underpins many "fast" algorithms in
computational physics. In particular, it is the foundation of the
Barnes-Hut [@1986_barnes_hut] and fast multipole methods [@rokhlin1987]
for evaluating all $O(n^{2})$ pairwise interactions between $n$
electrically charged particles in linear or close-to-linear complexity.
These methods were generalized by Hackbusch and co-workers, who
developed the $\mathcal{H}$- and $\mathcal{H}^{2}$-matrix frameworks
[@hackbusch; @2003_hackbusch; @2002_hackbusch_H2]. These explicitly
linear algebraic formulations enable fast algorithms not only for
matrix-vector multiplication but also for matrix inversion, LU
factorization, matrix-matrix multiplication and many more.

A fundamental challenge that arises when rank-structured matrix
formulations are used is how to find the data sparse representation of
the operator in the first place. The straightforward approach would be
to form the full matrix, and then loop over all the compressible off
diagonal blocks and compress them using, for example, a singular value
decomposition. The cost of such a process is necessarily at least
$O(n^{2})$, which is rarely affordable. Using randomized compression
techniques, it turns out to be possible to compress all the off-diagonal
blocks *jointly*, without any need for sampling each block individually.
In this section, we describe two such methods. Both require the user to
supply a fast algorithm for applying the full matrix (and its adjoint)
to vectors; see Section
[20.3](#sec:cpkernel_application){reference-type="ref"
reference="sec:cpkernel_application"}. The two methods have different
computational profiles.

- The technique described in Sections
  [20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"} and
  [20.5](#sec:HODLRcompression){reference-type="ref"
  reference="sec:HODLRcompression"} is a true "black-box" technique that
  interacts with $\bm{\mathsf{K}}$ only through the matrix-vector
  multiplication. It has storage complexity $O(n\log n)$ and it requires
  $O(\log n)$ applications of $\bm{\mathsf{K}}$ to a random matrix of
  size $n\times (r+p)$ where $r$ is an upper bound on the ranks of the
  off-diagonal blocks, and $p$ is a small over-sampling parameter.

- The technique described in Sections
  [20.6](#sec:HBS){reference-type="ref" reference="sec:HBS"} and
  [20.7](#sec:HBScompression){reference-type="ref"
  reference="sec:HBScompression"} attains true linear $O(n)$ complexity,
  and requires only a single application of $\bm{\mathsf{K}}$ and
  $\bm{\mathsf{K}}^{*}$ to a random matrix of size $n\times (r+p)$,
  where $r$ and $p$ are as in (a). Its drawbacks are that it requires
  evaluation of $O(rn)$ individual matrix entries, and that it works
  only for a smaller class of matrices.

::: remark
**Remark 69** (Scope of Section
[20](#sec:rankstructured){reference-type="ref"
reference="sec:rankstructured"}). *To keep the presentation as
uncluttered by burdensome notation as possible, in this survey we
restrict attention to two basic "formats" for representing a
rank-structured hierarchical matrix. In Sections
[20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"} and
[20.5](#sec:HODLRcompression){reference-type="ref"
reference="sec:HODLRcompression"}, we use the *hierarchically
off-diagonal low rank (HODLR)* format and in Sections
[20.6](#sec:HBS){reference-type="ref" reference="sec:HBS"} and
[20.7](#sec:HBScompression){reference-type="ref"
reference="sec:HBScompression"} we use the *hierarchically block
separable (HBS)* format (sometimes referred to as *hierarchically semi
separable (HSS)* matrices). The main limitation of these formats is that
they require all off-diagonal blocks of the matrix to have low numerical
rank. This is realistic only when the points
$\{\bm{\mathsf{x}}_{i}\}_{i=1}^{n}$ are restricted to a low-dimensional
manifold. In practical applications, one sometimes has to leave a larger
part of the matrix uncompressed, to avoid attempting to impose a
separation of variables such as
([\[eq:kernel_sepvar\]](#eq:kernel_sepvar){reference-type="ref"
reference="eq:kernel_sepvar"}) on a kernel
$k = k(\bm{\mathsf{x}},\bm{\mathsf{y}})$ when $\bm{\mathsf{x}}$ and
$\bm{\mathsf{y}}$ are too close. This is done through enforcing what is
called a "strong admissibility condition" (in contrast to the "weak
admissibility condition" of the HODLR and HBS formats), as was done in
the original Barnes-Hut and fast multipole methods.*
:::

::: {#remark:compressioncompetition .remark}
**Remark 70** (Alternative compression strategies). *Let us briefly
describe what alternatives to randomized compression exist. The original
papers on $\mathcal{H}$-matrices used Taylor approximations to derive a
separation of variables, but this works only when the kernel is given
analytically. It also tends to be quite expensive. The *adaptive
cross-approximation (ACA)* technique of [@2002_rjasanow_ACA] and
[@2006_bebendorf_ACA] relies on using "natural basis" vectors (see
Section [13.1](#sec:CURandIDdef){reference-type="ref"
reference="sec:CURandIDdef"}) that are found using semi-heuristic
techniques. The method can work very well in practice, but is not
guaranteed to provide an accurate factorization. When the kernel matrix
comes from mathematical physics, specialized techniques that exploit
mathematical properties of the kernel function often perform well
[@2005_martinsson_fastdirect]. For a detailed discussion of compression
of rank-structured matrices, we refer to [@2019_book Ch. 17].*
:::

## Computational environments {#sec:cpkernel_application}

The compression techniques that we describe rely on the user providing a
fast algorithm for applying the operator to be compressed to vectors.
Representative environments where such fast algorithms are available
include the following.

*Matrix-matrix multiplication:* Suppose that
$\bm{\mathsf{K}} = \bm{\mathsf{B}}\bm{\mathsf{C}}$, where
$\bm{\mathsf{B}}$ and $\bm{\mathsf{C}}$ are matrices that can rapidly be
applied to vectors. Then the randomized compression techniques allow us
to compute their product.

*Compression of boundary integral operators:* It is often possible to
reformulate a boundary value problem involving an elliptic partial
different operator, such as the Laplace or the Helmholtz operators, as
an equivalent boundary integral equation (BIE), see [@2019_book Part
III]. When such a BIE is discretized, the result is a kernel matrix
which can be applied rapidly to vectors using fast summation techniques
such as the fast multipole method [@rokhlin1987]. Randomized compression
techniques allow us to build an approximation to the matrix that can be
factorized or inverted, thus enabling direct (as opposed to iterative)
solvers.

*Dirichlet-to-Neumann (DtN) operators:* A singular integral operator of
central importance in engineering, physics, and scientific computing is
the DtN operator which maps given Dirichlet data for an elliptic
boundary value problem to the boundary fluxes of the corresponding
solution. The kernel of the DtN operator is rarely known analytically,
but the operator can be applied through a fast solver for the PDE, such
as, e.g., a finite-element discretization combined with a multigrid
solver. The randomized techniques described allow for the DtN operator
to be built and stored explicitly.

*Frontal matrices in sparse direct solvers:* Suppose that
$\bm{\mathsf{S}}$ is a large sparse matrix arising from the
discretization of an elliptic PDEs. A common technique for solving
$\bm{\mathsf{S}}\bm{\mathsf{x}} = \bm{\mathsf{b}}$ is to compute an LU
factorization of the matrix $\bm{\mathsf{S}}$. There are ways to do this
that preserve sparsity as far as possible, but in the course of the
factorization procedure, certain dense matrices of increasing size will
need to be factorized. These matrices turn out to be kernel matrices
with kernels that are not known explicitly, but that can be built using
the techniques described here
[@2013_xia_randomized; @2017_ghysels_robust].

## Hierarchically off-diagonal low-rank matrices {#sec:HODLR}

In order to illustrate how randomized methods can be used to construct
data sparse representations of rank-structured matrices, we will
describe a particularly simple 'format' in this section that is often
referred to as the hierarchically off-diagonal low-rank (HODLR) format.
This is a basic format that works well when the points are organized on
a one- or two-dimensional manifold.

The first step towards defining the HODLR format is to build a binary
tree on the index vector $I = [1,2,3,\dots,n]$ through a process that is
illustrated in Figure [11](#fig:tree){reference-type="ref"
reference="fig:tree"}(a). With any node $\tau$ in the tree, we associate
an index vector $I_{\tau} \subseteq I$. The root of the tree is given
the index $\tau=1$, and we associate it with the full index vector, so
that $I_{1} = I$. At the next finer level of the tree, we split $I_{1}$
into two parts $I_{2}$ and $I_{3}$ so that $I_{1} = I_{2} \cup I_{3}$
forms a disjoint partition. Then continue splitting the index vectors
until each remaining index vector is "short". (Exactly what "short"
means is application-dependent, but one may think of a short vector as
holding a few hundred indices or so.) We let $\ell$ denote a *level* of
the tree, with $\ell=0$ denoting the root, so that level $\ell$ holds
$2^{\ell}$ nodes. We use the terms *parent* and *child* in the obvious
way, and say that a pair of nodes $\{\alpha,\beta\}$ forms a *sibling
pair* if they have the same parent. A *leaf* node is of course a node
that has no children.

The binary tree that we defined induces a natural tessellation of the
kernel matrix $\bm{\mathsf{K}}$ into $O(n)$ blocks. Figure
[11](#fig:tree){reference-type="ref" reference="fig:tree"}(b) shows the
tessellation that follows from the tree in Figure
[11](#fig:tree){reference-type="ref" reference="fig:tree"}(a). Each
parent node $\tau$ in the tree gives rise to two off-diagonal blocks
that both have low numerical rank. Letting $\{\alpha,\beta\}$ denote the
children of $\tau$, these two blocks are $$\begin{equation}
\label{eq:siblingboxes}
\bm{\mathsf{K}}_{\alpha,\beta} = \bm{\mathsf{K}}(I_{\alpha},I_{\beta})
\qquad\mbox{and}\qquad
\bm{\mathsf{K}}_{\beta,\alpha} = \bm{\mathsf{K}}(I_{\beta},I_{\alpha}).
\end{equation}$$ For each leaf node $\tau$, we define a corresponding
diagonal block as $$\begin{equation}
\label{eq:diagbox}
\bm{\mathsf{D}}_{\tau} = \bm{\mathsf{K}}(I_{\tau},I_{\tau}).
\end{equation}$$ The disjoint partition of $\bm{\mathsf{K}}$ into blocks
is now formed by all sibling pairs, as defined in
([\[eq:siblingboxes\]](#eq:siblingboxes){reference-type="ref"
reference="eq:siblingboxes"}), together with all diagonal blocks, as
defined by ([\[eq:diagbox\]](#eq:diagbox){reference-type="ref"
reference="eq:diagbox"}). When each off-diagonal block in this
tessellation has low rank, we say that $\bm{\mathsf{K}}$ is a
*hierarchically off-diagonal low-rank (HODLR)* matrix.

<figure id="fig:tree">
<div class="picture">
<p>(120,120) (00,115)<span>(a)</span> (25,85)</p>
<div class="picture">
<p>(122,44) ( 20, 0)<span><embed src="Pics/fig_tree.pdf"
style="width:75mm" /></span> ( 0,30)<span><em>Level <span
class="math inline">0</span>:</em></span> ( 0,21)<span><em>Level <span
class="math inline">1</span>:</em></span> ( 0,12)<span><em>Level <span
class="math inline">2</span>:</em></span> ( 0,03)<span><em>Level <span
class="math inline">3</span>:</em></span></p>
</div>
<p>(00,77)<span>(b)</span> (40,02)</p>
<div class="picture">
<p>(80,80) (-2,-2)<span><img src="Pics/fig_8x8_Smatrix.png"
style="width:82mm" alt="image" /></span> (00,73)<span><span
class="math inline">D<sub>8</sub></span></span> (10,63)<span><span
class="math inline">D<sub>9</sub></span></span> (20,53)<span><span
class="math inline">D<sub>10</sub></span></span> (30,43)<span><span
class="math inline">D<sub>11</sub></span></span> (40,33)<span><span
class="math inline">D<sub>12</sub></span></span> (50,23)<span><span
class="math inline">D<sub>13</sub></span></span> (60,13)<span><span
class="math inline">D<sub>14</sub></span></span> (70,03)<span><span
class="math inline">D<sub>15</sub></span></span> (00,63)<span><span
class="math inline">K<sub>9, 8</sub></span></span> (10,73)<span><span
class="math inline">K<sub>8, 9</sub></span></span> (60,03)<span><span
class="math inline">K<sub>15, 14</sub></span></span> (70,13)<span><span
class="math inline">K<sub>14, 15</sub></span></span> (40,23)<span><span
class="math inline">K<sub>13, 12</sub></span></span> (50,33)<span><span
class="math inline">K<sub>12, 13</sub></span></span> (20,43)<span><span
class="math inline">K<sub>11, 10</sub></span></span> (30,53)<span><span
class="math inline">K<sub>10, 11</sub></span></span> (16,18)<span><span
class="math inline">K<sub>3, 2</sub></span></span> (56,58)<span><span
class="math inline">K<sub>2, 3</sub></span></span> (06,48)<span><span
class="math inline">K<sub>5, 4</sub></span></span> (26,68)<span><span
class="math inline">K<sub>4, 5</sub></span></span> (46,08)<span><span
class="math inline">K<sub>7, 6</sub></span></span> (66,28)<span><span
class="math inline">K<sub>6, 7</sub></span></span> (-10,38)<span><span
class="math inline">K=</span></span></p>
</div>
</div>
<figcaption>(a) A binary tree with three levels; see Section <a
href="#sec:HODLR" data-reference-type="ref"
data-reference="sec:HODLR">20.4</a>. Each node <span
class="math inline"><em>τ</em></span> in the tree owns an index vector
<span class="math inline"><em>I</em><sub><em>τ</em></sub></span> that is
a subset of the full index vector <span
class="math inline"><em>I</em> = 1 : <em>n</em></span>. For the root
node, we set <span
class="math inline"><em>I</em><sub>1</sub> = <em>I</em></span>. Each
split in the tree represents a disjoint partition of the corresponding
index vectors, so that, for example, <span
class="math inline"><em>I</em><sub>1</sub> = <em>I</em><sub>2</sub> ∪ <em>I</em><sub>3</sub></span>
and <span
class="math inline"><em>I</em><sub>7</sub> = <em>I</em><sub>14</sub> ∪ <em>I</em><sub>15</sub></span>.
(b) A matrix <span class="math inline">K</span> tessellated according to
the tree shown in (a). </figcaption>
</figure>

## Compressing a rank-structured hierarchical matrix through the matrix-vector multiplication only {#sec:HODLRcompression}

In this section, we describe a randomized technique for computing a data
sparse representation of a matrix that is compressible in the "HODLR"
format that we introduced in Section
[20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"}. This
technique interacts with $\bm{\mathsf{K}}$ only through the application
of $\bm{\mathsf{K}}$ and $\bm{\mathsf{K}}^{*}$ to vectors, and is in
this sense a "black-box" technique. In particular, we do not need the
ability to evaluate individual entries of $\bm{\mathsf{K}}$. The
following theorem summarizes the main result:

::: theorem
**Theorem 71**. *Let $\bm{\mathsf{K}}$ be a HODLR matrix associated with
a fully populated binary tree on the index vector, as described in
Section [20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"}.
Suppose that the tree has $L$ levels, that each off-diagonal block has
rank at most $k$, and that each leaf node in the tree holds at most $ck$
indices for some fixed number $c$. Then the diagonal blocks
$\bm{\mathsf{D}}_{\tau}$, as well as rank-$k$ factorizations of all
sibling interaction matrices $\bm{\mathsf{K}}_{\alpha,\beta}$ can be
computed by a randomized algorithm with cost at most
$$T_{\mathrm{total}} = T_{\mathrm{matvec}} \times (4L+c)k + T_{\mathrm{flop}} \times O(L^{2}k^{2}n),$$
where $T_{\mathrm{matvec}}$ is the cost of applying either
$\bm{\mathsf{K}}$ or $\bm{\mathsf{K}}^{*}$ to a vector, and
$T_{\mathrm{flop}}$ is the cost of a floating point operation.*
:::

The proof of the theorem consists of an explicit algorithm for building
all matrices that is often referred to as the "peeling algorithm". It
was originally published in [@2011_lin_lu_ying], with later
modifications proposed in [@2016_martinsson_hudson2]. The proof below is
based on [@2019_book Sec. 17.4]. We give the proof for the case
described in the theorem involving a matrix whose off-diagonal blocks
have exact rank $k$. In practical applications, it is of course more
typical to have the off-diagonal blocks be only of approximate low rank.
In this case, the exact same algorithm can be used, but the number of
samples drawn should be increased from $k$ to $k+p$ for some modest
oversampling parameter $p$.

::: proof
*Proof.* The algorithm is a "top-down" technique that compresses the
largest blocks first, and then moves on to process one level at a time
of successively smaller blocks. In describing the technique, we assume
that the blocks are numbered as shown in Figure
[11](#fig:tree){reference-type="ref" reference="fig:tree"}(a).

In the first step of the algorithm, we build approximations to the two
largest blocks $\bm{\mathsf{K}}_{2,3}$ and $\bm{\mathsf{K}}_{3,2}$,
shown in red in Figure [12](#fig:peeling){reference-type="ref"
reference="fig:peeling"}(a). To do this, we form a random matrix
$\bm{\mathsf{\Omega}}$ of size $n\times 2k$ and then use the
matrix-vector multiplication to form a sample matrix
$$\begin{array}{cccccccccccc}
\bm{\mathsf{Y}} &=& \bm{\mathsf{K}}&\bm{\mathsf{\Omega}}.\\
n \times 2k && n \times n & n\times 2k
\end{array}$$ The key idea of the peeling algorithm is to insert zero
blocks in the random matrix $\bm{\mathsf{\Omega}}$ in the pattern shown
in Figure [12](#fig:peeling){reference-type="ref"
reference="fig:peeling"}(a). The zero blocks permit us to extract "pure"
samples from the column spaces of the blocks $\bm{\mathsf{K}}_{2,3}$ and
$\bm{\mathsf{K}}_{3,2}$. For instance, the blocks labeled
$\bm{\mathsf{Y}}_{2}$ and $\bm{\mathsf{Y}}_{3}$ in the figure are given
by the formulas $$\begin{array}{cccccccccccc}
\bm{\mathsf{Y}}_{2} &=& \bm{\mathsf{K}}_{2,3}&\bm{\mathsf{\Omega}}_{3}\\
\tfrac{n}{2} \times k && \tfrac{n}{2} \times \tfrac{n}{2} & \tfrac{n}{2}\times k
\end{array}
\qquad\mbox{and}\qquad
\begin{array}{cccccccccccc}
\bm{\mathsf{Y}}_{3} &=& \bm{\mathsf{K}}_{3,2}&\bm{\mathsf{\Omega}}_{2}\\
\tfrac{n}{2} \times k && \tfrac{n}{2} \times \tfrac{n}{2} & \tfrac{n}{2}\times k
\end{array}$$ By orthonormalizing the matrices $\bm{\mathsf{Y}}_{2}$ and
$\bm{\mathsf{Y}}_{3}$, we obtain ON bases $\bm{\mathsf{U}}_{2}$ and
$\bm{\mathsf{U}}_{3}$ for the off-diagonal blocks
$\bm{\mathsf{K}}_{2,3}$ and $\bm{\mathsf{K}}_{3,2}$. In order to
complete the factorization of $\bm{\mathsf{K}}_{2,3}$ and
$\bm{\mathsf{K}}_{3,2}$, we will perform an operation that is the
equivalent of "Stage B" in Section
[11.2](#sec:rsvd){reference-type="ref" reference="sec:rsvd"}. To this
end, we form a test matrix by interlacing the matrices
$\bm{\mathsf{U}}_{2}$ and $\bm{\mathsf{U}}_{3}$ with zero blocks, to
form the $n\times 2k$ matrix shown in Figure
[13](#fig:peeling_trans){reference-type="ref"
reference="fig:peeling_trans"}. Applying $\bm{\mathsf{K}}^{*}$ to this
test matrix, we get the sample matrices
$$\bm{\mathsf{Z}}_{2} = \bm{\mathsf{K}}_{3,2}^{*}\bm{\mathsf{U}}_{3},
\qquad\mbox{and}\qquad
\bm{\mathsf{Z}}_{3} = \bm{\mathsf{K}}_{2,3}^{*}\bm{\mathsf{U}}_{2}.$$
Since $\bm{\mathsf{U}}_{2}$ holds an ON-basis for the column space of
$\bm{\mathsf{K}}_{2,3}$, it follows that
$$\bm{\mathsf{K}}_{2,3} = \bm{\mathsf{U}}_{2}\bm{\mathsf{U}}_{2}^{*}\bm{\mathsf{K}}_{2,3} = \bm{\mathsf{U}}_{2}\bm{\mathsf{Z}}_{3}^{*},$$
which establishes the low rank factorization of
$\bm{\mathsf{K}}_{2,3}$.[^2] The block $\bm{\mathsf{K}}_{3,2}$ is of
course factorized analogously.

<figure id="fig:peeling">
<table>
<tbody>
<tr>
<td style="text-align: left;">(a) In the first step of the algorithm,
the sibling interaction matrices <span
class="math inline">K<sub>2, 3</sub></span> and <span
class="math inline">K<sub>3, 2</sub></span> on level 1, shown in red,
are compressed.</td>
<td style="text-align: left;"></td>
</tr>
<tr>
<td style="text-align: left;">(b) In the second step, the four sibling
interaction matrices on level 2, shown in red, are compressed. In this
step, we exploit that we now possess factorizations of the gray
blocks.</td>
<td style="text-align: left;"></td>
</tr>
<tr>
<td style="text-align: left;">(c) In the third step, the eight sibling
interaction matrices on level 8, shown in red, are compressed. We again
exploit that we possess factorizations of the gray blocks.</td>
<td style="text-align: left;"></td>
</tr>
</tbody>
</table>
<figcaption>The “peeling algorithm” for computing a HODLR representation
of a matrix described in Section <a href="#sec:HODLRcompression"
data-reference-type="ref"
data-reference="sec:HODLRcompression">20.5</a>.</figcaption>
</figure>

At the second step of the algorithm, the objective is to build
approximations to the sibling matrices at the next finer level, shown in
red in Figure [12](#fig:peeling){reference-type="ref"
reference="fig:peeling"}(b). We use a random matrix
$\bm{\mathsf{\Omega}}$ of the same size, $n\times 2k$, as in the
previous step, but now with four zero blocks, as shown in the figure.
Looking at the sample matrix $\bm{\mathsf{Y}}_{4}$, we find that it
takes the form $$\begin{array}{cccccccccccc}
\bm{\mathsf{Y}}_{4} &=& \bm{\mathsf{K}}_{4,5}&\bm{\mathsf{\Omega}}_{5} &+& \bm{\mathsf{B}}_{4},\\
\tfrac{n}{4} \times k && \tfrac{n}{4} \times \tfrac{n}{4} & \tfrac{n}{4}\times k && \tfrac{n}{4}\times k
\end{array}$$ where $\bm{\mathsf{B}}_{4}$ represents contributions from
interactions between the block $\bm{\mathsf{K}}_{2,3}$ and the random
matrix and $\bm{\mathsf{\Omega}}_{7}$. But now observe that we at this
point are in possession of a low rank approximation to
$\bm{\mathsf{K}}_{2,3}$. This allows us to compute $\bm{\mathsf{B}}_{4}$
explicitly, and then obtain a "pure" sample from the column space of
$\bm{\mathsf{K}}_{4,5}$ by subtracting this contribution out:
$$\bm{\mathsf{Y}}_{4} - \bm{\mathsf{B}}_{4} = \bm{\mathsf{K}}_{4,5}\bm{\mathsf{\Omega}}_{5}.$$
By orthonormalizing the matrix
$\bm{\mathsf{Y}}_{4} - \bm{\mathsf{B}}_{4}$, we obtain the ON matrix
$\bm{\mathsf{U}}_{4}$ whose columns span the column space of
$\bm{\mathsf{K}}_{4,5}$. The same procedure of course also allows us to
build bases for the column spaces of all matrices at this level, and we
can then extract bases for the row spaces and the cores of the matrices
by sampling $\bm{\mathsf{K}}^{*}$.

At the third step of the algorithm, the object is to compress the blocks
marked in red in Figure [12](#fig:peeling){reference-type="ref"
reference="fig:peeling"}(c). The random matrix used remains of size
$n\times 2k$ but now contains 8 zero blocks, as shown in the figure. The
sample matrix $\bm{\mathsf{Y}}_{8}$ takes the form
$$\begin{array}{cccccccccccc}
\bm{\mathsf{Y}}_{8} &=& \bm{\mathsf{K}}_{8,9}&\bm{\mathsf{\Omega_{9}}} &+& \bm{\mathsf{B}}_{8}\\
\tfrac{n}{8} \times k && \tfrac{n}{8} \times \tfrac{n}{8} & \tfrac{n}{8}\times k && \tfrac{n}{8}\times k
\end{array}$$ where $\bm{\mathsf{B}}_{8}$ represents contributions from
off-diagonal blocks at the coarser levels, shown in gray in the figure.
Since we already possess data sparse representations of these blocks, we
can subtract their contributions out, just as at the previous level.

Once all levels have been processed, we hold low-rank factorizations of
all off-diagonal blocks, and all that remains is to extract the diagonal
matrices $\bm{\mathsf{D}}_{\tau}$ for all leaf nodes $\tau$. Let us for
simplicity assume that all such block are of the same size $m\times m$.
We then form a test matrix by stacking $n/m$ copies of an $m\times m$
identity matrix atop of each other and applying $\bm{\mathsf{K}}$ to
this test matrix. We subtract off the contributions from all the
off-diagonal blocks using the representation we have on hand, and are
then left with a sample matrix consisting of all the blocks
$\bm{\mathsf{D}}_{\tau}$ stacked on top of each other. ◻
:::

![The process used to complete the factorization of the blocks
$\bm{\mathsf{K}}_{2,3}$ and $\bm{\mathsf{K}}_{3,2}$ once the basis
matrices $\bm{\mathsf{U}}_{2}$ and $\bm{\mathsf{U}}_{3}$ have been
constructed. We put $\bm{\mathsf{U}}_{2}$ and $\bm{\mathsf{U}}_{3}$ into
the $n\times 2k$ test matrix shown, and form the sample matrix that
holds $\bm{\mathsf{Z}}_{2}$ and $\bm{\mathsf{Z}}_{3}$ by applying
$\bm{\mathsf{K}}^{*}$.](Pics/fig_peeling_level0_trans.pdf){#fig:peeling_trans
width="85mm"}

## A linear complexity data sparse format {#sec:HBS}

The "HODLR" data sparse matrix format discussed in Sections
[20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"} and
[20.5](#sec:HODLRcompression){reference-type="ref"
reference="sec:HODLRcompression"} is simple to use, and is in many
applications efficient enough. However, its storage complexity is at
least $O(n\log n)$ as the matrix size $n$ grows. We will next describe
how the logarithmic factor can be eliminated. To keep the presentation
concise, we focus on the asymptotic storage requirements, although
similar estimates hold for the asymptotic flop count.

To start, let us first investigate why there are logarithmic factors in
the storage requirement in the HODLR format. Suppose that for a sibling
pair $\{\alpha,\,\beta\}$ the corresponding off-diagonal block
$\bm{\mathsf{K}}_{\alpha,\beta}$ is of size $m\times m$, has rank $k$,
and is stored in the form of a factorization
$$\begin{array}{ccccccccccccccc}
\bm{\mathsf{K}}_{\alpha,\beta} &=& \bm{\mathsf{U}}^{\mathrm{long}}_{\alpha}&\tilde{\bm{\mathsf{K}}}_{\alpha,\beta}&\bigl(\bm{\mathsf{V}}_{\beta}^{\mathrm{long}}\bigr)^{*}\\
m\times m && m\times k & k\times k& k\times m
\end{array}$$ where $\bm{\mathsf{U}}_{\alpha}^{\mathrm{long}}$ and
$\bm{\mathsf{V}}_{\beta}^{\mathrm{long}}$ are two matrices whose columns
form bases for the column and row space of
$\bm{\mathsf{K}}_{\alpha,\beta}$, respectively.[^3] We see that
$\bm{\mathsf{K}}_{\alpha,\beta}$ can be stored using about $2mk$
floating point numbers (where we ignore $O(k^{2})$ terms, since
typically $k \ll m$). Now let us consider how many floats are required
to store all sibling interaction matrices on a given level $\ell$ in the
tree. There are $2^{\ell}$ such matrices, and they each have size
$m \approx 2^{-\ell}n$, resulting in a total of
$2^{\ell}\times 2\times 2^{-\ell}nk = 2nk$ floats. Since there are
$O(\log(n))$ levels in the tree, it follows that just the task of
storing all the basis matrices requires $O(kn\log(n))$ storage.

A standard technique for overcoming the problem of storing all the
"long" basis matrices is to assume that the basis matrices on one level
can be formed through slight modifications to the basis matrices on the
next finer level. To be precise, we assume that if $\tau$ is a node in
the tree with children $\alpha$ and $\beta$, then there exists a "short"
basis matrix $\bm{\mathsf{U}}_{\tau}$ of size $2k\times k$ such that
$$\begin{equation}
\label{eq:nestedbasis}
\begin{array}{ccccccccccccc}
\bm{\mathsf{U}}_{\tau}^{\mathrm{long}} &=&
\left[\begin{array}{cc} \bm{\mathsf{U}}_{\alpha}^{\mathrm{long}} & \bm{\mathsf{0}} \\ \bm{\mathsf{0}} & \bm{\mathsf{U}}_{\beta}^{\mathrm{long}}\end{array}\right]
&
\bm{\mathsf{U}}_{\tau}\\
m\times m && m\times 2k & 2k \times k
\end{array}
\end{equation}$$ (with of course an analogous statement holding for the
"V" matrices holding bases for the row spaces). The point is that if we
already have the "long" basis matrices for the children available, then
the long basis matrices for the parent can be formed using the
information in the small matrix $\bm{\mathsf{U}}_{\tau}$. By applying
this argument recursively, we see that for any parent node $\tau$, all
that needs to be stored explicitly are small matrices of size
$2k\times k$.

To illustrate the idea, suppose that we use this technique for storing
basis matrices to the tree with three levels described in Section
[20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"}. Then for
each of the 8 leaf boxes, we compute basis matrices of size
$n/8 \times k$ and store these. At the next coarser level, the long
basis matrix for a box such as $\tau=4$ can be expressed through the
formula $$\begin{array}{ccccccccccccc}
\bm{\mathsf{U}}_{4}^{\mathrm{long}} &=&
\left[\begin{array}{cc} \bm{\mathsf{U}}_{8}^{\mathrm{long}} & \bm{\mathsf{0}} \\ \bm{\mathsf{0}} & \bm{\mathsf{U}}_{9}^{\mathrm{long}}\end{array}\right]
&
\bm{\mathsf{U}}_{4},\\
(n/4)\times m && (n/4)\times 2k & 2k \times k
\end{array}$$ so that only the $2k\times k$ matrix $\bm{\mathsf{U}}_{4}$
needs to be stored. At the next coarser level still, the long basis
matrix $\bm{\mathsf{U}}_{2}^{\mathrm{long}}$ is expressed as
$$\begin{array}{ccccccccccccc}
\bm{\mathsf{U}}_{2}^{\mathrm{long}} &=&
\left[\begin{array}{cccc}
\bm{\mathsf{U}}_{8}^{\mathrm{long}}  & \bm{\mathsf{0}}                 & \bm{\mathsf{0}}                 & \bm{\mathsf{0}} \\
\bm{\mathsf{0}}                 & \bm{\mathsf{U}}_{9}^{\mathrm{long}}  & \bm{\mathsf{0}}                 & \bm{\mathsf{0}} \\
\bm{\mathsf{0}}                 & \bm{\mathsf{0}}                 & \bm{\mathsf{U}}_{10}^{\mathrm{long}} & \bm{\mathsf{0}} \\
\bm{\mathsf{0}}                 & \bm{\mathsf{0}}                 & \bm{\mathsf{0}}                 & \bm{\mathsf{U}}_{11}^{\mathrm{long}}
\end{array}\right]
&
\left[\begin{array}{cc} \bm{\mathsf{U}}_{4} & \bm{\mathsf{0}} \\ \bm{\mathsf{0}} & \bm{\mathsf{U}}_{5}\end{array}\right]
&
\bm{\mathsf{U}}_{2}.\\
(n/2)\times m && (n/2)\times 4k & 4k\times 2k & 2k \times k
\end{array}$$

We say that a matrix that can be represented using nested basis matrices
in the manner described is a *hierarchically block separable (HBS)*
matrix. (This format is very closely related to the *hierarchically semi
separable (HSS)* matrix format described in
[@2010_gu_xia_HSS; @2005_gu_HSS].) When the HBS format is used, it is
natural to assume that each leaf index vector holds $O(k)$ indices,
which means that there are overall $O(n/k)$ nodes in the tree. Since we
only need $O(k^2)$ storage per node, we see that the overall storage
requirements improve from $O(kn\log(n))$ to $O(kn)$.

There is a price to be paid for using "nested" basis matrices, and it is
that the rank $k$ required to reach a requested precision $\varepsilon$
typically is higher when nested basis matrices are used, in comparison
to the HODLR format. To be precise, one can show that for the relation
([\[eq:nestedbasis\]](#eq:nestedbasis){reference-type="ref"
reference="eq:nestedbasis"}) to hold, it is necessary and sufficient
that for a node $\alpha$, the columns of
$\bm{\mathsf{U}}_{\alpha}^{\mathrm{long}}$ must span the column space of
the matrix $\bm{\mathsf{K}}(I_{\alpha},I_{\alpha}^{\mathrm{c}})$ (where
$I_{\alpha}^{\mathrm{c}} = I \backslash I_{\alpha}$). In contrast, in
the HODLR format, they only need to span the columns of
$\bm{\mathsf{K}}(I_{\alpha},I_{\beta})$, where $\beta$ is the sibling of
$\alpha$. This is easier to do since $I_{\beta}$ is a subset of the
index vector $I_{\alpha}^{\mathrm{c}}$.

## A linear complexity randomized compression technique {#sec:HBScompression}

We saw in Section [20.6](#sec:HBS){reference-type="ref"
reference="sec:HBS"} that the so called "HBS" data sparse format allows
us to store a matrix using $O(kn)$ floating point numbers, without any
factor of $\log(n)$. But is it also possible to compute such a
representation in optimal linear complexity? The answer is yes, and a
number of deterministic techniques that work in specific environments
have been proposed [@2019_book Ch. 17]. A linear complexity randomized
technique was proposed in [@2011_martinsson_randomhudson]. This method
is based on analyzing samples of the matrix obtained through the
application of $\bm{\mathsf{K}}$ and $\bm{\mathsf{K}}^{*}$ to Gaussian
random vectors, but is not quite a black-box method, as it also requires
the ability to evaluate $O(kn)$ individual elements of $\bm{\mathsf{K}}$
itself. This is not always possible, of course, but when it is, the
technique is lightning fast in practice. To be precise,
[@2011_martinsson_randomhudson] establishes the following:

::: {#thm:randHBS .theorem}
**Theorem 72**. *Let $\bm{\mathsf{K}}$ be an $n\times n$ matrix that is
compressible in the HBS format described in Section
[20.6](#sec:HBS){reference-type="ref" reference="sec:HBS"}, for some
rank $k$. Let $T_{\mathrm{apply}}$ denote the time it takes to evaluate
the two products $$\begin{align}
\label{eq:randomhudsonsample1}
\bm{\mathsf{Y}} =&\ \bm{\mathsf{K}}\bm{\mathsf{G}},\\
\label{eq:randomhudsonsample2}
\bm{\mathsf{Z}} =&\ \bm{\mathsf{K}}^{*}\bm{\mathsf{G}},
\end{align}$$ where $\bm{\mathsf{G}}$ is an $n\times k$ matrix drawn
from a Gaussian distribution. Then a full HBS representation of
$\bm{\mathsf{K}}$ can be computed at cost bounded by
$$T_{\mathrm{total}} = T_{\mathrm{apply}} + T_{\mathrm{entry}} \times O(nk) + T_{\mathrm{flop}} \times O(nk^{2}),$$
where $T_{\mathrm{flop}}$ is the cost of a floating point operation, and
where $T_{\mathrm{entry}}$ is the time it takes to evaluate an
individual entry of $\bm{\mathsf{K}}$*
:::

When the off-diagonal blocks are only approximately of rank $k$, we do
some over sampling as usual, and replace the number $k$ in the theorem
by $k+p$ for some small number $p$. This of course results in an
approximate HBS representation of $\bm{\mathsf{K}}$.

The proof of Theorem [72](#thm:randHBS){reference-type="ref"
reference="thm:randHBS"} is an algorithm that explicitly builds the data
sparse representation of the matrix within the specified time budget.
The algorithm consists of a pass through all nodes in the hierarchical
tree, going from smaller boxes to larger; "bottom-up", as opposed to the
"top-down" method for HODLR matrices in Section
[20.5](#sec:HODLRcompression){reference-type="ref"
reference="sec:HODLRcompression"}. We provide an outline of the proof
below. This outline is written to convey the main ideas without
introducing cumbersome notation; for full details, see the original
article [@2011_martinsson_randomhudson] or [@2019_book Sec. 17.4].

::: proof
*Proof.* To describe the algorithm that establishes Theorem
[72](#thm:randHBS){reference-type="ref" reference="thm:randHBS"}, we
will walk through how it applies to a matrix tessellated in accordance
with the hierarchical tree in Figure
[11](#fig:tree){reference-type="ref" reference="fig:tree"}(a). We start
by showing how we can use the information in the matrix
$\bm{\mathsf{Y}}$ defined by
([\[eq:randomhudsonsample1\]](#eq:randomhudsonsample1){reference-type="ref"
reference="eq:randomhudsonsample1"}) to build the basis matrix
$\bm{\mathsf{U}}_{8}$. As we saw in Section
[20.6](#sec:HBS){reference-type="ref" reference="sec:HBS"}, we need the
columns of $\bm{\mathsf{U}}_{8}$ to span the columns in the submatrix
$\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}})$. We achieve this by
constructing the sample matrix $$\bm{\mathsf{Y}}_{8} :=
\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}})\bm{\mathsf{G}}(I_{8}^{\mathrm{c}},:).$$
Since $\bm{\mathsf{G}}(I_{8}^{\mathrm{c}},:)$ is a Gaussian random
matrix of the appropriate size, the columns of $\bm{\mathsf{Y}}_{8}$
will form an approximate basis for the column space of
$\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}})$. But how do we form
$\bm{\mathsf{Y}}_{8}$ from the sample matrix $\bm{\mathsf{Y}}$ defined
by
([\[eq:randomhudsonsample1\]](#eq:randomhudsonsample1){reference-type="ref"
reference="eq:randomhudsonsample1"})? The method we use is illustrated
in Figure [14](#fig:hbs_randomized){reference-type="ref"
reference="fig:hbs_randomized"}(a), where the block
$\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}})$ consists of the gray block
inside the red rectangle. We see that in order to isolate the product
between $\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}})$ and
$\bm{\mathsf{G}}(I_{8}^{\mathrm{c}},:)$, all we need to do is to
subtract the contribution from the diagonal block
$\bm{\mathsf{K}}(I_{8},I_{8})$ (marked in white in the figure). In other
words, $$\begin{equation}
\label{eq:Y8}
\bm{\mathsf{Y}}_{8} =
\bm{\mathsf{K}}(:,I_{8})\bm{\mathsf{G}} - \bm{\mathsf{K}}(I_{8},I_{8})\bm{\mathsf{G}}(I_{8},:) =
\bm{\mathsf{Y}}(I_{8},:) - \bm{\mathsf{K}}(I_{8},I_{8})\bm{\mathsf{G}}(I_{8},:).
\end{equation}$$ The formula ([\[eq:Y8\]](#eq:Y8){reference-type="ref"
reference="eq:Y8"}) can be evaluated inexpensively since the block
$\bm{\mathsf{K}}(I_{8},I_{8})$ is small, and can be explicitly formed
since we assume that we have access to individual entries of the matrix
$\bm{\mathsf{K}}$. Once $\bm{\mathsf{Y}}_{8}$ has been formed, we
compute its row interpolatory decomposition to form a basis matrix
$\bm{\mathsf{U}}_{8}$ and an index vector $\tilde{I}_{8}^{\mathrm{row}}$
such that
$$\bm{\mathsf{K}}(I_{8},I_{8}^{\mathrm{c}}) = \bm{\mathsf{U}}_{8}\,\bm{\mathsf{K}}(\tilde{I}_{8}^{\mathrm{row}},I_{8}^{\mathrm{c}}).$$
In an entirely analogous way, we form a sample matrix
$\bm{\mathsf{Y}}_{9}$ whose columns span the block
$\bm{\mathsf{K}}(I_{9},I_{9}^{\mathrm{c}})$ through the formula
$$\bm{\mathsf{Y}}_{9} =
\bm{\mathsf{Y}}(I_{9},:) - \bm{\mathsf{K}}(I_{9},I_{9})\bm{\mathsf{G}}(I_{9},:),$$
cf. Figure [14](#fig:hbs_randomized){reference-type="ref"
reference="fig:hbs_randomized"}(b), and then compute the row ID
$\{\bm{\mathsf{U}}_{9},\tilde{I}_{9}^{\mathrm{row}}\}$ of
$\bm{\mathsf{Y}}_{9}$.

<figure id="fig:hbs_randomized">
<embed src="Pics/fig_hbs_randomized.pdf" />
<p>(a) (b) (c)</p>
<figcaption>The algorithm that establishes Theorem <a
href="#thm:randHBS" data-reference-type="ref"
data-reference="thm:randHBS">72</a> is illustrated using a matrix <span
class="math inline">K</span> tessellated in accordance with the tree
shown in Figure <a href="#fig:tree" data-reference-type="ref"
data-reference="fig:tree">11</a>(a). (a) The block <span
class="math inline">K(<em>I</em><sub>8</sub>,:)</span> is marked with a
red rectangle. Within the red rectangle, the shaded area represents the
block <span
class="math inline">K(<em>I</em><sub>8</sub>, <em>I</em><sub>8</sub><sup>c</sup>)</span>.
(b) The block <span
class="math inline">K(<em>I</em><sub>9</sub>,:)</span> is marked with a
red rectangle, with <span
class="math inline">K(<em>I</em><sub>9</sub>, <em>I</em><sub>9</sub><sup>c</sup>)</span>
shaded. (c) The block <span
class="math inline">K(<em>I</em><sub>4</sub>,:)</span> is marked with a
red rectangle, with <span
class="math inline">K(<em>I</em><sub>4</sub>, <em>I</em><sub>4</sub><sup>c</sup>)</span>
shaded. </figcaption>
</figure>

The basis matrices $\bm{\mathsf{V}}_{\tau}$ that span the row spaces of
the offdiagonal blocks for any leaf $\tau$ are built through the same
procedure, but starting with the sample matrix $\bm{\mathsf{Z}}$ defined
by
([\[eq:randomhudsonsample2\]](#eq:randomhudsonsample2){reference-type="ref"
reference="eq:randomhudsonsample2"}). For instance, we form the sample
matrix $$\begin{equation}
\label{eq:Z8}
\bm{\mathsf{Z}}_{8} =
\bm{\mathsf{K}}(I_{8},:)^{*}\bm{\mathsf{G}} - \bm{\mathsf{K}}(I_{8},I_{8})^{*}\bm{\mathsf{G}}(I_{8},:) =
\bm{\mathsf{Z}}(I_{8},:) - \bm{\mathsf{K}}(I_{8},I_{8})^{*}\bm{\mathsf{G}}(I_{8},:),
\end{equation}$$ and then compute the ID of $\bm{\mathsf{Z}}_{8}$ to
build a basis matrix $\bm{\mathsf{V}}_{8}$ and an index vector
$\tilde{I}_{8}^{\mathrm{col}}$ such that
$$\bm{\mathsf{K}}(I_{8}^{\mathrm{c}},I_{8}) = \bm{\mathsf{K}}(I_{8}^{\mathrm{c}},\tilde{I}_{8}^{\mathrm{col}})\bm{\mathsf{V}}_{8}^{*}.$$

The choice to use the interpolatory decomposition in factorizing the
off-diagonal blocks is essential to making the linear complexity
compression scheme work. In particular, observe that once we have formed
the row and column IDs of all the leaf boxes, we automatically obtain
rank-$k$ factorizations of all the corresponding sibling interaction
matrices. For instance, if $\{\alpha,\beta\}$ is a sibling pair
consisting of two leaf nodes, then the $m\times m$ sibling interaction
matrix $\bm{\mathsf{K}}(I_{\alpha},I_{\beta})$ admits the factorization
$$\begin{equation}
\label{eq:randHBSsibs}
\begin{array}{cccccccccccccc}
\bm{\mathsf{K}}(I_{\alpha},I_{\beta}) &=& \bm{\mathsf{U}}_{\alpha} & \bm{\mathsf{K}}(\tilde{I}_{\alpha}^{\mathrm{row}},\tilde{I}_{\beta}^{\mathrm{col}}) & \bm{\mathsf{V}}_{\beta}^{*}.\\
m\times m && m\times k & k\times k & k\times m
\end{array}
\end{equation}$$ In order to evaluate
([\[eq:randHBSsibs\]](#eq:randHBSsibs){reference-type="ref"
reference="eq:randHBSsibs"}), we merely need to form the matrix
$\bm{\mathsf{K}}(\tilde{I}_{\alpha}^{\mathrm{row}},\tilde{I}_{\beta}^{\mathrm{col}})$.

Once the leaf nodes have all been processed, we proceed to the next
coarser level. Consider for instance the node $\tau=4$ with children
$\alpha=8$ and $\beta=9$. Our task is in principle to build the long
basis matrix $\bm{\mathsf{U}}_{4}^{\mathrm{long}}$ that spans the
columns of $\bm{\mathsf{K}}(I_{4},I_{4}^{\mathrm{c}})$, which is the
shaded matrix inside the red rectangle in Figure
[14](#fig:hbs_randomized){reference-type="ref"
reference="fig:hbs_randomized"}(c). (We say "in principle" since it will
not actually be explicitly formed.) To this end, we define the sample
matrix
$$\bm{\mathsf{Y}}_{4} := \bm{\mathsf{K}}(I_{4},I_{4}^{\mathrm{c}})\bm{\mathsf{G}}(I_{4}^{\mathrm{c}},:).$$
The idea is now to repeat the same technique that we used for the leaf
boxes, and write $$\begin{equation}
\label{eq:Y4}
\bm{\mathsf{Y}}_{4} = \bm{\mathsf{Y}}(I_{4},:) - \bm{\mathsf{K}}(I_{4},I_{4})\bm{\mathsf{G}}(I_{4},:).
\end{equation}$$ It turns out that we can evaluate
([\[eq:Y4\]](#eq:Y4){reference-type="ref" reference="eq:Y4"}) very
efficiently: The block $\bm{\mathsf{K}}(I_{4},I_{4})$ is made up of the
two diagonal blocks $\bm{\mathsf{K}}(I_{8},I_{8})$ and
$\bm{\mathsf{K}}(I_{9},I_{9})$, and the two off-diagonal blocks
$\bm{\mathsf{K}}(I_{8},I_{9})$ and $\bm{\mathsf{K}}(I_{9},I_{8})$. Now
observe that at this point in the execution of the algorithm, we have
compressed representations of all these blocks available. Using this
information, we form a short sample matrix $\tilde{\bm{\mathsf{Y}}}_{4}$
of size $2k\times k$ that we can compress to build the basis matrix
$\bm{\mathsf{U}}_{4}$ and the associated index vector
$\tilde{I}_{4}^{\mathrm{row}}$ in an ID of
$\tilde{\bm{\mathsf{Y}}}_{4}$. ◻
:::

## Butterfly matrices {#sec:butterfly}

An interesting class of rank-structured matrices arises from a
generalization of the discrete Fourier transform. These so called
"butterfly matrices" provide a data sparse format for many matrices that
appear in the analysis of neural networks, wave propagation problems,
and signal processing
[@dao2019learning; @2007_oneil_thesis; @2009_candes_demanet_ying_butterfly].
This format involves an additional complication in comparison to the
HODLR and HBS formats in that it requires $O(\log n)$ different
tessellations of the matrix, as illustrated for a simple case in Figure
[15](#fig:butterfly){reference-type="ref" reference="fig:butterfly"}. It
can be demonstrated that when all of the resulting submatrices are of
numerically low rank, the matrix as a whole can be written
(approximately) as a product of $O(\log n)$ sparse matrices, in a manner
analogous to the butterfly representation of an FFT [@1995_briggs_DFT
Sec. 10.4].

![Illustration of the *butterfly* rank-structured matrix format
described in Section [20.8](#sec:butterfly){reference-type="ref"
reference="sec:butterfly"}. The figures show the blocks that must all be
of numerically low rank for a butterfly matrix arising from a binary
tree with 16 leaf nodes.](Pics/fig_butterfly.pdf){#fig:butterfly
width="\\textwidth"}

Randomization has proven to be a powerful tool for finding butterfly
representations of matrices. The techniques involved are more complex
than the methods described in Sections
[20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"} --
[20.6](#sec:HBS){reference-type="ref" reference="sec:HBS"} due to the
multiplicative nature of the representation, and typically involve
iterative refinement rather than direct approximation
[@2015_ying_butterfly; @2017_interpolative_butterfly; @2018_ying_multidimensional_butterfly].
Similar techniques played an essential role in a recent ground breaking
paper [@2017_michielssen_butterfly_IE] that exploits butterfly
representations to directly solve linear systems arising from the
modeling of scattering problems in the high frequency regime.

## Applications of rank-structured matrices in data analysis {#sec:ASKIT}

The machinery for working with rank-structured hierarchical matrices
that we have described can be used also for the kernel matrices
discussed in Section [19](#sec:kernel){reference-type="ref"
reference="sec:kernel"} that arise in machine learning and computational
statistics. For instance, demonstrate that data sparse formats of this
type can be very effective for simulating Gaussian processes in low
dimensional spaces.

When the underlying dimension grows, the techniques that we have
described so far become uncompetitive. Even $d=4$ would be considered a
stretch. Fortunately, significant progress has recently been made
towards extending the essential ideas to higher dimensions. For
instance, describe a technique that is designed to uncover intrinsic
lower dimensional structures that are often present in sets of points
$\{\bm{\mathsf{x}}_{i}\}_{i=1}^{n}$ that ostensibly live in higher
dimensional spaces. The idea is to use randomized algorithms both for
organizing the points into a hierarchical tree (that induces the
tessellation of the matrix) and for computing low rank approximations to
the resulting admissible blocks. The authors report promising numerical
results for a wide selection of kernel matrices.

In this context, it is as we saw in Section
[19](#sec:kernel){reference-type="ref" reference="sec:kernel"} rarely
possible to execute a matrix-vector multiplication, or even to evaluate
more than a tiny fraction of the entries of the matrix. This means on
the one hand that sampling must form an integral part of the compression
strategy, and on the other that firm performance guarantees are
typically not available. The saving grace is that when a kernel matrix
is used for learning and data analysis, a rough approximation is often
sufficient.

::: remark
**Remark 73** (Geometry oblivious methods). *A curious observation is
that techniques developed for kernel matrices appear to also be
applicable for certain symmetric positive definite (pd) matrices that
are not explicitly presented as kernel matrices. This is a consequence
of the well known fact that any pd matrix $\bm{\mathsf{K}}$ admits a
factorization $$\begin{equation}
\label{eq:pd_K_fact}
\bm{\mathsf{K}} = \bm{\mathsf{G}}^{*}\bm{\mathsf{G}}
\end{equation}$$ for a so called "Gramian matrix" $\bm{\mathsf{G}}$. (If
the eigenvalue decomposition of $\bm{\mathsf{K}}$ takes the form
$\bm{\mathsf{K}} = \bm{\mathsf{U}}\bm{\mathsf{\Lambda}}\bm{\mathsf{U}}^{*}$,
then a matrix of the form
$\bm{\mathsf{G}} = \bm{\mathsf{V}}\bm{\mathsf{\Lambda}}^{1/2}\bm{\mathsf{U}}^{*}$
is a Gramian if and only if $\bm{\mathsf{V}}$ is unitary.) The
factorization ([\[eq:pd_K_fact\]](#eq:pd_K_fact){reference-type="ref"
reference="eq:pd_K_fact"}) says that the entries of $\bm{\mathsf{K}}$
are formed by the inner products between the columns of
$\bm{\mathsf{G}}$, $$\begin{equation}
\label{eq:gofmm_kernel}
\bm{\mathsf{K}}(i,j) = \langle \bm{\mathsf{g}}_{i}, \, \bm{\mathsf{g}}_{j} \rangle = k(\bm{\mathsf{g}}_{i},\bm{\mathsf{g}}_{i}),
\end{equation}$$ where $\bm{\mathsf{g}}_{i}$ is the $i$'th column of
$\bm{\mathsf{G}}$ (the "Gram vector") and where the $k$ is the inner
product kernel of Example [64](#ex:innerprodkernel){reference-type="ref"
reference="ex:innerprodkernel"}. At this point, it becomes plausible
that the techniques of [@2015_biros_ASKIT] for kernel matrices
associated with points in high dimensional spaces may apply to certain
pd matrices. The key to make this work is the observation that it is not
necessary to explicitly form the Gram factors $\bm{\mathsf{G}}$. All
that is needed in order to organize the points
$\{\bm{\mathsf{g}}_{i}\}_{i=1}^{n}$ are relative distances and angles
between the points, and we can evaluate these from the matrix entries of
$\bm{\mathsf{K}}$, via the formula
$$\|\bm{\mathsf{g}}_{i} - \bm{\mathsf{g}}_{j}\|^{2} =
\|\bm{\mathsf{g}}_{i}\|^{2} - 2\mbox{Re}\,\langle \bm{\mathsf{g}}_{i}, \, \bm{\mathsf{g}}_{j} \rangle + \|\bm{\mathsf{g}}_{j}\|^{2} =
\bm{\mathsf{K}}(i,i) - 2\mbox{Re}\,\bm{\mathsf{K}}(i,j) + \bm{\mathsf{K}}(j,j).$$
The resulting technique was presented in [@2017_biros_GOFMM] as a
"geometry oblivious FMM (GOFMM)", along with numerical evidence of its
usefulness for important classes of matrices.*
:::

[^1]: The celebrated method of has exponent
    $\omega = \log_{2}(7) = 2.807\cdots$. It is a compelling algorithm
    in terms of both its numerical stability and its practical speed,
    even for modest matrix sizes. More exotic algorithms, such as the
    Coppersmith--Winograd method and variants, attain complexity of
    about $\omega \approx 2.37$, but they are not considered to be
    practically useful.

[^2]: One may if desired continue the factorization process and form an
    "economy size" SVD of $\bm{\mathsf{Z}}_{3}^{*}$ so that
    $\bm{\mathsf{Z}}_{3}^{*} = \hat{\bm{\mathsf{U}}}_{2}\bm{\mathsf{\Sigma}}_{2,3}\bm{\mathsf{V}}_{3}^{*}$.
    This results in the factorization
    $\bm{\mathsf{K}}_{2,3} = \bigl(\bm{\mathsf{U}}_{2}\hat{\bm{\mathsf{U}}}_{2}\bigr)\bm{\mathsf{\Sigma}}_{2,3}\bm{\mathsf{V}}_{3}^{*}$,
    which is an SVD of $\bm{\mathsf{K}}_{2,3}$. However, for purposes of
    establishing the theorem, we may leave $\bm{\mathsf{Z}}_{3}$ alone,
    and let $\{\bm{\mathsf{U}}_{2},\bm{\mathsf{Z}}_{3}\}$ be our
    "compressed" representation of $\bm{\mathsf{K}}_{2,3}$.

[^3]: The matrices $\bm{\mathsf{U}}_{\alpha}^{\mathrm{long}}$ and
    $\bm{\mathsf{V}}_{\beta}^{\mathrm{long}}$ were denoted
    $\bm{\mathsf{U}}_{\alpha}$ and $\bm{\mathsf{V}}_{\beta}$ in Sections
    [20.4](#sec:HODLR){reference-type="ref" reference="sec:HODLR"} and
    [20.5](#sec:HODLRcompression){reference-type="ref"
    reference="sec:HODLRcompression"}.
