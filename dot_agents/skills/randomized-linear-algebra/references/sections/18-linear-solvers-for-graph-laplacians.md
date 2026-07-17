# Linear solvers for graph Laplacians {#sec:sparse-cholesky}

In this section we describe the randomized algorithm,
[SparseCholesky]{.smallcaps}, which can efficiently solve a linear
system whose coefficient matrix is a graph Laplacian matrix. Up to a
polylogarithmic factor, this method achieves the minimum possible
runtime and storage costs. This algorithm has the potential to
accelerate many types of computations involving graph Laplacian
matrices.

The [SparseCholesky]{.smallcaps} algorithm was developed by and further
refined by . These approaches are based on earlier work from Dan
Spielman's group, notably the paper of . The presentation here is
adapted from .

## Overview

We begin with a high-level approach for solving Laplacian linear
systems. The basic idea is to construct a preconditioner using a
randomized variant of the incomplete Cholesky method. Then we can solve
the original linear system by means of the preconditioned conjugate
gradient (PCG) algorithm.

### Approximate solutions to the Poisson problem

Let $\bm{\mathsf{L}} \in \mathbb{H}_n(\mathbb{R})$ be the Laplacian
matrix of a weighted, loop-free, undirected graph on the vertex set
$V = \{1, \dots, n\}$. We write $m$ for the number of edges in the
graph, i.e., the sparsity of the graph. For simplicity, we will also
assume that the graph is *connected*; equivalently,
$\ker(\bm{\mathsf{L}}) = \operatorname{span}\{ \mathbf{1} \}$ where
$\mathbf{1} \in \mathbb{R}^n$ is the vector of ones. See
Section [7.4](#sec:graph-sparsification){reference-type="ref"
reference="sec:graph-sparsification"} for definitions. See
Figure [8](#fig:big-dipper){reference-type="ref"
reference="fig:big-dipper"} for an illustration of an unweighted,
undirected graph.

The basic goal is to find the unique solution $\bm{\mathsf{x}}_{\star}$
to the Poisson problem
$$\bm{\mathsf{L}} \bm{\mathsf{x}} = \bm{\mathsf{f}}
\quad\text{where}\quad
\mathbf{1}^* \bm{\mathsf{f}} = {0}
\quad\text{and}\quad
\mathbf{1}^* \bm{\mathsf{x}} = {0}.$$ For a parameter $\varepsilon> 0$,
we can relax this requirement by asking instead for an approximate
solution $\bm{\mathsf{x}}_{\varepsilon}$ that satisfies
$$\Vert  \bm{\mathsf{x}}_{\varepsilon} - \bm{\mathsf{x}}_{\star}  \Vert_{\bm{\mathsf{L}}}
    \leq \varepsilon\, \Vert  \bm{\mathsf{x}}_{\star}  \Vert_{\bm{\mathsf{L}}}.$$
We have written
$\Vert  \bm{\mathsf{x}}  \Vert_{\bm{\mathsf{L}}} := (\bm{\mathsf{x}}^* \bm{\mathsf{L}} \bm{\mathsf{x}})^{1/2}$
for the energy seminorm induced by the Laplacian matrix.

### Approximate Cholesky decomposition

Imagine that we can efficiently construct a sparse, approximate Cholesky
decomposition of the Laplacian matrix $\bm{\mathsf{L}}$. More precisely,
we seek a morally lower-triangular matrix $\bm{\mathsf{C}}$ that
satisfies $$\begin{equation}
 \label{eqn:approx-cholesky}
0.5 \, \bm{\mathsf{L}} \preccurlyeq\bm{\mathsf{CC}}^* \preccurlyeq 1.5 \, \bm{\mathsf{L}}
    \quad\text{where}\quad
    \texttt{nnz}(\bm{\mathsf{C}}) = O(m \log^2 n).
\end{equation}$$ In other words, there is a known permutation of the
rows that brings the matrix $\bm{\mathsf{C}}$ into lower-triangular
form. As usual, $\preccurlyeq$ is the semidefinite order, and `nnz`
returns the number of nonzero entries in a matrix.

This section describes an algorithm, called
[SparseCholesky]{.smallcaps}, that can complete the task outlined in the
previous paragraph. This algorithm is motivated by the insight that we
can produce a sparse approximation of a Laplacian matrix by random
sampling (Section [7.4](#sec:graph-sparsification){reference-type="ref"
reference="sec:graph-sparsification"}). The main challenge is to obtain
sampling probabilities without extra computation. The resulting method
can be viewed as a randomized variant of the incomplete Cholesky
factorization [@GVL13:Matrix-Computations-4ed Sec. 11.5.8].

### Preconditioning

Given the sparse, approximate Cholesky factor $\bm{\mathsf{C}}$, we can
precondition the Poisson problem: $$\begin{equation}
 \label{eqn:precond-poisson}
(\bm{\mathsf{C}}^\dagger\bm{\mathsf{L}} \bm{\mathsf{C}}^{*\dagger})(\bm{\mathsf{C}}^* \bm{\mathsf{x}}) = (\bm{\mathsf{C}}^\dagger\bm{\mathsf{f}}).
\end{equation}$$
When [\[eqn:approx-cholesky\]](#eqn:approx-cholesky){reference-type="eqref"
reference="eqn:approx-cholesky"} holds, the matrix
$\bm{\mathsf{C}}^\dagger\bm{\mathsf{L}} \bm{\mathsf{C}}^{*\dagger}$ has
condition number $\kappa \leq 3$.

Therefore, we can solve the preconditioned
system [\[eqn:precond-poisson\]](#eqn:precond-poisson){reference-type="eqref"
reference="eqn:precond-poisson"} quickly using the PCG algorithm
[@GVL13:Matrix-Computations-4ed Sec. 11.5]. If the initial iterate
$\bm{\mathsf{x}}_0 = \bm{\mathsf{0}}$, then $j$ steps of PCG produce an
iterate $\bm{\mathsf{x}}_j$ that satisfies
$$\Vert  \bm{\mathsf{x}}_j - \bm{\mathsf{x}}_{\star}  \Vert_{\bm{\mathsf{L}}}
    \leq 2 \left[ \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right]^j
    \Vert  \bm{\mathsf{x}}_{\star}  \Vert_{\bm{\mathsf{L}}}
    < 3^{1-j} \Vert  \bm{\mathsf{x}}_{\star}  \Vert_{\bm{\mathsf{L}}}.$$
As a consequence, we can achieve relative error $\varepsilon$ after
$1 + \log_3 (1/\varepsilon)$ iterations. Each iteration requires a
matrix--vector product with $\bm{\mathsf{L}}$ and the solution of a
(consistent) linear system
$\bm{\mathsf{CC}}^* \bm{\mathsf{u}} = \bm{\mathsf{y}}$. We can perform
these steps in $O(m \log^2 n)$ operations per iteration. Indeed,
$\bm{\mathsf{L}}$ has only $O(m)$ nonzero entries. The matrix
$\bm{\mathsf{C}}$ is morally triangular with $O(m \log^2 n)$ nonzero
entries, so we can apply $(\bm{\mathsf{CC}}^*)^{\dagger}$ using two
triangular solves.

::: remark
**Remark 55** (Triangular solve). *To solve a consistent linear system
$(\bm{\mathsf{CC}}^*) \bm{\mathsf{u}} = \bm{\mathsf{y}}$, we can apply
$\bm{\mathsf{C}}^\dagger$ by triangular elimination and then apply
$\bm{\mathsf{C}}^{*\dagger}$ by triangular elimination. The four
subspace theorem ensures that solution to the first problem renders the
second problem consistent too. Since
$\ker(\bm{\mathsf{CC}}^*) = \operatorname{span}\{\bm{\mathsf{1}}\}$, we
can enforce consistency numerically by removing the constant component
of the input $\bm{\mathsf{y}}$. Similarly, we can remove the constant
component of the output $\bm{\mathsf{u}}$ to ensure it belongs to the
correct space.*
:::

### Main results

The following theorem describes the performance of the
[SparseCholesky]{.smallcaps} procedure. This is the main result from .

::: {#thm:sparse-cholesky .theorem}
**Theorem 56** (SparseCholesky). *Let $\bm{\mathsf{L}}$ be the Laplacian
of a connected graph on $n$ vertices, with $m$ weighted edges. With high
probability, the [SparseCholesky]{.smallcaps} algorithm produces a
morally lower-triangular matrix $\bm{\mathsf{C}}$ that satisfies
$$0.5 \, \bm{\mathsf{L}} \preccurlyeq\bm{\mathsf{CC}}^* \preccurlyeq 1.5 \, \bm{\mathsf{L}}.$$
The matrix $\bm{\mathsf{C}}$ has $O(m \log^2 n)$ nonzero entries. The
expected running time is $O(m \log^3 n)$ operations.*
:::

In view of our discussion about PCG, we arrive at the following
statement about solving the Poisson problem.

::: corollary
**Corollary 57** (Poisson problem). *Suppose the
[SparseCholesky]{.smallcaps} algorithm delivers an approximation
$\bm{\mathsf{L}} \approx \bm{\mathsf{CC}}^*$ that
satisfies [\[eqn:approx-cholesky\]](#eqn:approx-cholesky){reference-type="eqref"
reference="eqn:approx-cholesky"}. Then we can solve each consistent
linear system $\bm{\mathsf{L}}\bm{\mathsf{x}} = \bm{\mathsf{f}}$ to
relative error $\varepsilon$ in the seminorm
$\Vert \cdot \Vert_{\bm{\mathsf{L}}}$ using at most
$1 + \log_3(1/\varepsilon)$ iterations of PCG, each with a cost of
$O(m \log^2 n)$ arithmetic operations.*
:::

### Discussion

As we have mentioned, the Poisson problem serves as a primitive for
undertaking many computations on undirected
graphs [@Ten10:Laplacian-Paradigm]. Potential applications include
clustering, analysis of random walks, and finite-element discretizations
of elliptic PDEs.

The [SparseCholesky]{.smallcaps} algorithm achieves a near-optimal
runtime and storage guarantee for the Poisson problem on a graph.
Indeed, for a general graph with $m$ edges, any algorithm must use
$O(m)$ storage and arithmetic. There is a proof that the cost
$O(m \log^{1/2}(n))$ is achievable in theory [@CKM+14:Solving-SDD], but
the resulting methods are currently impractical. Meanwhile, the
simplicity of the [SparseCholesky]{.smallcaps} method makes it a
candidate for real-world computation.

For particular classes of Laplacian matrices, existing solvers can be
very efficient. Optimized sparse direct solvers
[@2016_acta_sparse_direct_survey] work very well for small- and
medium-size problems, but they typically have superlinear scaling, which
renders them unsuitable for truly large-scale problems. Iterative
methods such as multigrid or preconditioned Krylov solvers can attain
linear complexity for important classes of problems, in particular for
sparse systems arising from the discretization of elliptic PDEs.
However, we are not aware of competing methods that provably enjoy
near-optimal complexity for all problems.

We regard the [SparseCholesky]{.smallcaps} algorithm as one of the most
dramatic examples of how randomization has the potential to accelerate
basic linear algebra computations, both in theory and in practice.

## Cholesky decomposition of a graph Laplacian

To begin our explanation of the [SparseCholesky]{.smallcaps} algorithm,
let us summarize what happens when we apply the standard Cholesky
decomposition method to a graph Laplacian.

### The Laplacian of a multigraph

For technical reasons, related to the design and analysis of the
algorithm, we need to work with multigraphs instead of ordinary graphs.
In the discussion, we will point out specific places where this
generality is important.

Consider a weighted, undirected *multigraph* $G$, defined on the vertex
set $V = \{1, \dots, n\}$. Each edge $e = \{u, v\}$ is an unordered pair
of vertices; we typically use the abbreviated notation $e = uv = vu$. We
introduce a weight function $w_G$ that assigns a positive weight to each
edge $e$ in the multigraph $G$. Since $G$ is a multigraph, there may be
many multiedges, with distinct weights, connecting the same two
vertices.

Taking some notational liberties, we will identify the multigraph $G$
with its Laplacian matrix $\bm{\mathsf{L}}$, which we express in the
form $$\begin{equation}
 \label{eqn:multi-laplacian}
\bm{\mathsf{L}} = \sum\nolimits_{e \in \bm{\mathsf{L}}} w_{\bm{\mathsf{L}}}(e) \, \bm{\mathsf{\Delta}}_{e}.
\end{equation}$$ The matrix $\bm{\mathsf{\Delta}}_{e}$ is the elementary
Laplacian on the vertex pair $(u, v)$ that composes the edge $e = uv$.
That is,
$$\bm{\mathsf{\Delta}}_e := \bm{\mathsf{\Delta}}_{uv} := (\bm{\mathsf{\delta}}_{u} - \bm{\mathsf{\delta}}_v)(\bm{\mathsf{\delta}}_{u} - \bm{\mathsf{\delta}}_v)^*
\quad\text{where $e = uv$.}$$ The sum
in [\[eqn:multi-laplacian\]](#eqn:multi-laplacian){reference-type="eqref"
reference="eqn:multi-laplacian"} takes place over all the multiedges $e$
in the multigraph $\bm{\mathsf{L}}$, so the same elementary Laplacian
may appear multiple times with different weights.

### Stars and cliques {#sec:star-clique}

To describe the Cholesky algorithm on a graph, we need to introduce a
few more concepts from graph theory. Define the *degree* and the *total
weight* of a vertex $u$ in the multigraph $\bm{\mathsf{L}}$ to be
$$\operatorname{deg}_{\bm{\mathsf{L}}}(u) := \sum\nolimits_{e=uv \in \bm{\mathsf{L}}} 1
\quad\text{and}\quad
w_{\bm{\mathsf{L}}}(u) := \sum\nolimits_{e = uv \in \bm{\mathsf{L}}} w_{\bm{\mathsf{L}}}(e).$$
In other words, the degree of $u$ is the total number of multiedges $e$
that contain $u$. The total weight of $u$ is the sum of the weights of
the multiedges $e$ that contain $u$.

Let $u$ be a fixed vertex. The *star* induced by $u$ is the Laplacian
$$\textsc{star}(u, \bm{\mathsf{L}}) := \sum\nolimits_{e = uv \in \bm{\mathsf{L}}} w_{\bm{\mathsf{L}}}(e) \, \bm{\mathsf{\Delta}}_e.$$
In words, the star includes precisely those multiedges $e$ in the
multigraph $\bm{\mathsf{L}}$ that contain the vertex $u$.

The *clique* induced by $u$ is defined implicitly as the correction that
occurs when we take the Schur
complement [\[eqn:schur-complement\]](#eqn:schur-complement){reference-type="eqref"
reference="eqn:schur-complement"} of the Laplacian with respect to the
coordinate $u$: $$\bm{\mathsf{L}} / \bm{\mathsf{\delta}}_u %
    =: (\bm{\mathsf{L}} - \textsc{star}(u, \bm{\mathsf{L}})) + \textsc{clique}(u, \bm{\mathsf{L}}).$$
Recall that $\bm{\mathsf{\delta}}_u$ is the standard basis vector in
coordinate $u$. By direct calculation, one may verify that
$$\textsc{clique}(u, \bm{\mathsf{L}}) = \frac{1}{2 w_{\bm{\mathsf{L}}}(u)} \sum\limits_{e_1 = uv_1 \in \bm{\mathsf{L}}} \sum\limits_{e_2 = uv_2 \in \bm{\mathsf{L}}}
    w_{\bm{\mathsf{L}}}(e_1) w_{\bm{\mathsf{L}}}(e_2) \bm{\mathsf{\Delta}}_{v_1v_2}.$$
Each sum takes place over all multiedges $e$ in $\bm{\mathsf{L}}$ that
contain the vertex $u$. It can be verified that the clique is also the
Laplacian of a weighted multigraph.

Figures [8](#fig:big-dipper){reference-type="ref"
reference="fig:big-dipper"}
and [9](#fig:star-clique){reference-type="ref"
reference="fig:star-clique"} contain an illustration of a (simple)
graph, along with the star and clique induced by eliminating a vertex.
In our more general setting, the edges in the star and clique would have
associated weights. These diagrams are courtesy of Richard Kueng.

<figure id="fig:big-dipper">

<figcaption><em>A combinatorial graph.</em> The Ursa Major graph with a
distinguished vertex (the star Megrez) highlighted in red.</figcaption>
</figure>

<figure id="fig:star-clique">
<table>
<tbody>
<tr>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
</tr>
</tbody>
</table>
<figcaption><em>Illustration of a star and clique in the Ursa Major
graph.</em> (a) The star induced by the red vertex consists of the three
solid black edges. (b) The clique induced by the red vertex consists of
the three solid black edges. These edges are added when the red vertex
is eliminated. </figcaption>
</figure>

### Graphs and Cholesky

With the notation introduced in
Section [18.2.2](#sec:star-clique){reference-type="ref"
reference="sec:star-clique"}, we can present the graph-theoretic
interpretation of the Cholesky algorithm as it applies to the Laplacian
$\bm{\mathsf{L}}$ of a weighted multigraph.

Define the initial Laplacian $\bm{\mathsf{S}}_0 := \bm{\mathsf{L}}$. In
each step $i = 1, 2, \dots, n$, we select a new vertex $u_i$. We extract
the associated column of the current Laplacian:
$$\bm{\mathsf{c}}_i := \frac{1}{\sqrt{(\bm{\mathsf{S}}_{i-1})_{u_i u_i}}} \bm{\mathsf{S}}_{i-1} \bm{\mathsf{\delta}}_{u_i}.$$
We compute the Schur complement with respect to the vertex $u_i$:
$$\bm{\mathsf{S}}_i := \bm{\mathsf{S}}_{i-1} / \bm{\mathsf{\delta}}_{u_i} = (\bm{\mathsf{S}}_{i-1} - \textsc{star}(u_i, \bm{\mathsf{S}}_{i-1})) + \textsc{clique}(u_i, \bm{\mathsf{S}}_{i-1}).$$
In other words, we remove the star induced by $u_i$ and replace it with
the clique induced by $u_i$. Each $\bm{\mathsf{S}}_i$ is the Laplacian
of a multigraph; it has no multiedge that contains any one of the
vertices $u_1, \dots, u_i$. Therefore, we have reduced the size of the
problem.

After $n$ steps, the Cholesky factorization is determined by a vector
$\bm{\mathsf{\pi}} = (u_1, u_2, \dots, u_n)$ that holds the chosen
indices and a matrix
$\bm{\mathsf{C}} = \begin{bmatrix} \bm{\mathsf{c}}_1 & \dots & \bm{\mathsf{c}}_n \end{bmatrix}$
such that $\bm{\mathsf{C}}(\bm{\mathsf{\pi}}, :)$ is lower-triangular.
The algorithm ensures that we have the exact decomposition
$$\bm{\mathsf{L}} = \bm{\mathsf{CC}}^*.$$ This is the (pivoted) Cholesky
factorization of the Laplacian matrix.

To choose a vertex to eliminate, the classical approach is to find a
vertex with minimum degree or with minimum total weight. Or one may
simply select one of the remaining vertices at random.

### Computational costs

The cost to compute a Cholesky factorization
$\bm{\mathsf{L}} = \bm{\mathsf{C}}\bm{\mathsf{C}}^{*}$ of a Laplacian
matrix $\bm{\mathsf{L}}$ is typically superlinear in $n$, with a
worst-case cost of $O(n^{3})$ arithmetic and $O(n^{2})$ storage. This is
the reason that $\bm{\mathsf{C}}$ is less sparse than $\bm{\mathsf{L}}$:
the clique that is introduced at an elimination step has more edges than
the star that it replaces, a phenomenon referred to as *fill-in*. The
exact growth in the number of nonzero entries depends on the sparsity
pattern of $\bm{\mathsf{L}}$ and on the chosen elimination order. For
special cases, using a nested dissection ordering can provably improve
on the worst-case estimates [@2016_acta_sparse_direct_survey]. For
instance, if $\bm{\mathsf{L}}$ results from the finite-difference or
finite-element discretization of an elliptic PDE, then
$\texttt{nnz}(\bm{\mathsf{C}}) = O(n\log(n))$ in two dimensions and
$\texttt{nnz}(\bm{\mathsf{C}}) = O(n^{4/3})$ in three.

For general graphs, one path towards improving the efficiency of the
Cholesky factorization procedure is to randomly approximate the clique
by sampling, in order to curb the fill-in.
Section [7.4](#sec:graph-sparsification){reference-type="ref"
reference="sec:graph-sparsification"} already indicates that this
innovation may be possible, provided that we can find a way to obtain
sampling probabilities.

## The [SparseCholesky]{.smallcaps} algorithm

We are now prepared to present the [SparseCholesky]{.smallcaps}
procedure, which uses randomized sampling to compute a sparse,
approximate Cholesky factorization.

### Procedure

Let $\bm{\mathsf{L}}$ be the Laplacian of a weighted multigraph on
$V = \{1, \dots, n\}$. We perform the following steps.

1.  **Preprocessing.** Split each multiedge $e = uv$ in
    $\bm{\mathsf{L}}$ into $R = \lceil 64 \log^2(\mathrm{e}n) \rceil$
    multiedges, each connecting $\{u, v\}$, and each with weight
    $w_{\bm{\mathsf{L}}}(e) / R$. The purpose of this step is to control
    the effective
    resistance [\[eqn:effective-resistance\]](#eqn:effective-resistance){reference-type="eqref"
    reference="eqn:effective-resistance"} of each multiedge at the
    outset of the algorithm. Note that this splitting results in a
    weighted multigraph, even if we begin with a simple graph.

2.  **Initialization.** Form the initial Laplacian
    $\bm{\mathsf{S}}_0 = \bm{\mathsf{L}}$ and the list of remaining
    vertices $F_0 = V$.

3.  **Iteration.** For each $i = 1, 2, \dots, n$:

    1.  **Select a vertex.** Choose a vertex $u_i$ uniformly at random
        from $F_{i-1}$. Remove this vertex from the list:
        $F_i = F_{i-1} \setminus \{ u_i \}$.

    2.  **Extract the column.** Copy the normalized $u_i$th column from
        the current Laplacian:
        $$\bm{\mathsf{c}}_i = \frac{1}{\sqrt{(\bm{\mathsf{S}}_{i-1})_{u_i u_i}}} \bm{\mathsf{S}}_{i-1} \bm{\mathsf{\delta}}_{u_i}.$$
        Set $\bm{\mathsf{c}}_i = \bm{\mathsf{0}}$ if the denominator
        equals zero.

    3.  **Sampling the clique.** Construct the Laplacian
        $\bm{\mathsf{K}}_i$ of a random sparse approximation of
        $\textsc{clique}(u_i, \bm{\mathsf{S}}_{i-1})$. We will detail
        this procedure in
        Section [18.3.2](#sec:clique-sample){reference-type="ref"
        reference="sec:clique-sample"}.

    4.  **Approximate Schur complement.** Form
        $$\bm{\mathsf{S}}_i = (\bm{\mathsf{S}}_{i-1} - \textsc{star}(u_i, \bm{\mathsf{S}}_{i-1})) + \bm{\mathsf{K}}_i.$$

4.  **Decomposition.** Collate the columns $\bm{\mathsf{c}}_i$ into the
    Cholesky factor
    $$\bm{\mathsf{C}} = \begin{bmatrix} \bm{\mathsf{c}}_1 & \dots & \bm{\mathsf{c}}_n \end{bmatrix}.$$
    Define the row permutation $\pi(i) = u_i$ for each $i$.

Once these operations are complete, $\bm{\mathsf{C}}$ is a sparse,
morally lower-triangular matrix. It is also very likely that
$\bm{\mathsf{L}} \approx \bm{\mathsf{CC}}^*$.
Theorem [56](#thm:sparse-cholesky){reference-type="ref"
reference="thm:sparse-cholesky"} makes a rigorous accounting of these
claims.

### Clique sampling {#sec:clique-sample}

The remaining question is how to construct a random approximation
$\bm{\mathsf{K}}$ of a clique $\textsc{clique}(u, \bm{\mathsf{S}})$.
Here is the procedure:

1.  **Probabilities.** Construct a probability mass $\bm{\mathsf{p}}$
    such that
    $$p(e) = \frac{w_{\bm{\mathsf{S}}}(e)}{w_{\bm{\mathsf{S}}}(u)}
    \quad\text{for each $e \in \textsc{star}(u, \bm{\mathsf{S}})$.}$$

2.  **Sampling.** For each
    $i = 1, \dots, d = \operatorname{deg}_{\bm{\mathsf{S}}}(u)$,

    1.  Draw a random multiedge $e_1 = uv_1$ from the multiedges in
        $\textsc{star}(u, \bm{\mathsf{S}})$ according to the probability
        mass $\bm{\mathsf{p}}$.

    2.  Draw a second random multiedge $e_2 = uv_2$ from the multiedges
        in $\textsc{star}(u, \bm{\mathsf{S}})$ according to the uniform
        distribution.

    3.  Form the random Laplacian matrix of a new multiedge:
        $$\bm{\mathsf{X}}_i = \frac{w_{\bm{\mathsf{S}}}(e_1) \, w_{\bm{\mathsf{S}}}(e_2)}{w_{\bm{\mathsf{S}}}(e_1) + w_{\bm{\mathsf{S}}}(e_2)} \bm{\mathsf{\Delta}}_{v_1 v_2}.$$

3.  **Approximation.** Return
    $\bm{\mathsf{K}} = \sum_{i=1}^d \bm{\mathsf{X}}_i$.

The key fact about this construction is that it produces an unbiased
estimator $\bm{\mathsf{K}}$ of the clique:
$$\operatorname{\mathbb{E}}\bm{\mathsf{K}} = \textsc{clique}(u, \bm{\mathsf{S}}).$$
Furthermore, each summand $\bm{\mathsf{X}}_i$ creates a multiedge with
uniformly bounded effective
resistance [\[eqn:effective-resistance\]](#eqn:effective-resistance){reference-type="eqref"
reference="eqn:effective-resistance"}. This property persists as the
[SparseCholesky]{.smallcaps} algorithm executes, and it ensures that the
random matrix $\bm{\mathsf{K}}$ has controlled variance. Note that the
sampling procedure can result in several edges between the same pair of
vertices, which is another reason we need the multigraph formalism.

Next, observe that the number $d$ of multiedges in $\bm{\mathsf{K}}$ is
no greater than the number $d$ of multiedges in the star that we are
removing from $\bm{\mathsf{S}}$. (For comparison, note that the full
clique has $d^2$ multiedges.) As a consequence, the clique approximation
is inexpensive to construct. Moreover, the total number of multiedges in
the Laplacian can only decrease as the [SparseCholesky]{.smallcaps}
algorithm proceeds.

### Analysis

The analysis of [SparseCholesky]{.smallcaps} is well beyond the scope of
this paper. The key technical tool is a concentration inequality for
matrix-valued martingales that was derived in
[@Oli09:Concentration-Adjacency] and [@Tro11:Freedmans-Inequality]. To
activate this result, [@KS16:Approximate-Gaussian] use the fact that the
random clique approximation is unbiased and low-variance, conditional on
previous choices made by the algorithm. The proof also relies heavily on
the fact that we eliminate a random vertex at each step of the
iteration. For the technical details,
see [@KS16:Approximate-Gaussian; @Kyn17:Approximate-Gaussian]
and [@Tro19:Matrix-Concentration-LN].

### Implementation

The [SparseCholesky]{.smallcaps} algorithm is fairly simple to describe,
but it demands some care to develop an implementation that achieves the
runtime guarantees stated in
Theorem [56](#thm:sparse-cholesky){reference-type="ref"
reference="thm:sparse-cholesky"}.

The most important point is that we need to use data structures for
weighted multigraphs. One method is to maintain the vertex--multiedge
adjacency matrix, along with a list of weights. This approach requires
sparse matrix libraries, including efficient iterators over the rows and
columns.

A secondary point is that we need fast methods for constructing finite
probability distributions and sampling from them repeatedly. See .

It is unlikely that the [SparseCholesky]{.smallcaps} procedure will fail
to produce a factor $\bm{\mathsf{C}}$ that
satisfies [\[eqn:approx-cholesky\]](#eqn:approx-cholesky){reference-type="eqref"
reference="eqn:approx-cholesky"}. Even so, it is reassuring to know that
we can detect failures. Indeed, we can estimate the extreme singular
values of the preconditioned linear
system [\[eqn:precond-poisson\]](#eqn:precond-poisson){reference-type="eqref"
reference="eqn:precond-poisson"} using the methods from
Section [6](#sec:max-eig){reference-type="ref" reference="sec:max-eig"}.

The failure probability can also be reduced by modifying the random
vertex selection rule. Instead, we can draw a random vertex whose total
weight is at most twice the average total weight of the remaining
vertices [@Kyn17:Approximate-Gaussian]. In practice, it may suffice to
use the classical elimination rules based on minimum degree or minimum
total weight.

The main shortcoming of the [SparseCholesky]{.smallcaps} procedure
arises from the initialization step, where we split each multiedge into
$O(\log^2 n)$ pieces. This step increases the storage and computation
costs of the algorithm enough to make it uncompetitive (e.g. with fast
direct Poisson solvers) for some problem instances. At present, it is
unclear whether the initialization step can be omitted or relaxed, while
maintaining the reliability and correctness of the algorithm.

