# General linear solvers {#sec:linear-solve}

Researchers are currently exploring randomized algorithms for solving
linear systems, such as $$\begin{equation}
\label{eq:Ax=b}
\bm{\mathsf{A}}\bm{\mathsf{x}} = \bm{\mathsf{b}},
\end{equation}$$ where $\bm{\mathsf{A}}$ is a given coefficient matrix
and $\bm{\mathsf{b}}$ is a given vector. This section describes a few
probabilistic approaches for solving
[\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref" reference="eq:Ax=b"}. For
the most part, we restrict our attention to the case where
$\bm{\mathsf{A}}$ is square and the system is consistent, but we will
also touch on linear regression problems. Research on randomized linear
solvers has not progressed as rapidly as some other areas of randomized
NLA, so the discussion here is more preliminary than other parts of this
survey.

## Background: Iterative solvers

It is important to keep in mind that existing iterative solvers often
work exceptionally well. Whenever $\bm{\mathsf{A}}$ is well-conditioned
or, more generally, whenever its spectrum is "clustered," Krylov solvers
such as the conjugate gradient (CG) algorithm or GMRES tend to converge
very rapidly. For practical purposes, the cost of solving
[\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref" reference="eq:Ax=b"} is
no larger than the cost of a handful of matrix--vector multiplications
with $\bm{\mathsf{A}}$. In terms of speed, it is very difficult to beat
these techniques. Consequently, we focus on the cases where known
iterative methods converge slowly and where we cannot deploy standard
preconditioners to resolve the problem.

Having limited ourselves to this situation, the choice of solver for
[\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref" reference="eq:Ax=b"} will
depend on properties of the coefficient matrix: Is it dense or sparse?
Does it fit in RAM? Do we have access to individual matrix entries? Can
we apply $\bm{\mathsf{A}}$ to a vector? We will consider several of
these environments.

## Accelerating solvers based on dense matrix factorizations {#sec:parkeretc}

As it happens, one of the early examples of randomization in NLA was a
method for accelerating the solution of a dense linear system
[\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref" reference="eq:Ax=b"}.
observed that we can precondition a linear system by left and right
multiplying the coefficient matrix by random unitary matrices
$\bm{\mathsf{U}}$ and $\bm{\mathsf{V}}$. With probability $1$, we can
solve the resulting system $$\begin{equation}
\label{eq:Ax=b_precond}
\bigl(\bm{\mathsf{U}}\bm{\mathsf{A}}\bm{\mathsf{V}}^{*}\bigr)\,\bigl(\bm{\mathsf{V}}\bm{\mathsf{x}}\bigr) = \bm{\mathsf{U}}\bm{\mathsf{b}}
\end{equation}$$ by Gaussian elimination *without pivoting*. More
precisely, Parker proved that, almost surely, blocked Gaussian
elimination will not encounter a degenerate diagonal block.

Blocked Gaussian elimination without pivoting is substantially faster
than ordinary Gaussian elimination for two reasons: matrix operations
are more efficient than vector operations on modern computers, and we
avoid the substantial communication costs that arise when we search for
pivots. (Section [16.2](#sec:blocking){reference-type="ref"
reference="sec:blocking"} contains more discussion about blocking.)
Parker also observed that structured random matrices (such as the
randomized trigonometric transforms from
Section [9.3](#sec:srtt){reference-type="ref" reference="sec:srtt"})
allow us to perform the preconditioning step at lower cost than the
subsequent Gaussian elimination procedure.

inspired many subsequent papers, including
[@2014_li_random_butterfly_pivoting; @2017_trogdon_random_butterfly; @DGHL12:Communication-Optimal-Parallel; @2017_baboulin_GPU]
and [@2017_pan_randomized_gaussian_elim]. Another related direction
concerns the smoothed analysis of Gaussian elimination undertaken in
[@SST06:Smoothed-Analysis].

As we saw in Section [16](#sec:full){reference-type="ref"
reference="sec:full"}, randomization can be used to accelerate the
computation of rank-revealing factorizations of the matrix
$\bm{\mathsf{A}}$. In this context, randomness allows us to block the
factorization method, which increases its practical speed, even though
the overall arithmetic cost remains at $O(n^3)$. Randomized
rank-revealing factorizations are ideal for solving ill-conditioned
linear systems because they allow the user to stabilize the computation
by avoiding subspaces associated with small singular values.

For instance, suppose that we have computed a singular value
decomposition (SVD): $$\bm{\mathsf{A}} =
\bm{\mathsf{U}}\bm{\mathsf{D}}\bm{\mathsf{V}}^{*} =
\sum_{j=1}^{n}\sigma_{j}\,\bm{\mathsf{u}}_{j}\bm{\mathsf{v}}_{j}^{*}.$$
Let us introduce a truncation parameter $\varepsilon$ and ignore all
singular modes where $\sigma_j \leq \varepsilon$. Then the stabilized
solution to [\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref"
reference="eq:Ax=b"} is
$$\bm{\mathsf{x}}_{\varepsilon} = \sum_{j \,:\, \sigma_{j} > \varepsilon}\frac{1}{\sigma_{j}}\,\bm{\mathsf{v}}_{j}\bm{\mathsf{u}}_{j}^{*}\bm{\mathsf{b}}.$$
By allowing the residual to take a nonzero value, we can ensure that
$\bm{\mathsf{x}}_\varepsilon$ does not include large components that
contribute little toward satisfying the original equation. The
randomized URV decomposition, described in
Section [16.7](#sec:randUTV){reference-type="ref"
reference="sec:randUTV"}, can also be used for stabilization, and we can
compute it much faster than an SVD.

::: remark
**Remark 53** (Are rank-revealing factorizations needed?). *In some
applications, computing a rank-revealing factorization is overkill for
purposes of solving the linear system
([\[eq:Ax=b\]](#eq:Ax=b){reference-type="ref" reference="eq:Ax=b"}). In
particular, if we compute an unpivoted QR decomposition of
$\bm{\mathsf{A}}$, then it is easy to block both the factorization and
the solve stages so that very high speed is attained. This process is
provably backwards stable, which is sometimes all that is needed. (In
practice, partially pivoted LU can often be used in an analogous manner,
despite being theoretically unstable.)*

*In contrast, when the actual entries of the computed solution
$\bm{\mathsf{x}}_{\mathrm{approx}}$ matter (as opposed to the value of
$\bm{\mathsf{A}}\bm{\mathsf{x}}_{\mathrm{approx}}$), a stabilized solver
is generally preferred. As a consequence, column-pivoted QR is often
cited as a method of choice for ill-conditioned problems in situations
where an SVD is not affordable.*
:::

::: remark
**Remark 54** (Strassen accelerated solvers). *We saw in Section
[16.4](#sec:demmelURV){reference-type="ref" reference="sec:demmelURV"}
that randomization has enabled us to compute a rank-revealing
factorization of an $n\times n$ matrix in less than $\mathcal{O}(n^{3})$
operations. The idea was to use randomized preconditioning as in
([\[eq:Ax=b_precond\]](#eq:Ax=b_precond){reference-type="ref"
reference="eq:Ax=b_precond"}), and then accelerate an unpivoted
factorization of the resulting coefficient matrix using fast algorithms
for the matrix-matrix multiplication such as Strassen
[@2007_demmel_fast_linear_algebra_is_stable]. This methodology can of
course be immediately applied to the task of solving ill-conditioned
linear systems. For improved numerical stability, a few steps of power
iteration can be incorporated to this approach; see
[\[eq:powerURV\]](#eq:powerURV){reference-type="eqref"
reference="eq:powerURV"}.*
:::

## Sketch and precondition {#sec:sketchtoprecond}

Another approach to preconditioning is to look for a random
transformation of the linear system that makes an iterative linear
solver converge more quickly. Typically, these preconditioning
transforms need to cluster the eigenvalues of the matrix.

The most successful example of this type of randomized preconditioning
does not concern square systems, but rather highly overdetermined
least-squares problems. See
Section [10.5](#sec:sketchandprecond){reference-type="ref"
reference="sec:sketchandprecond"} *et seq.* for a discussion of this
idea. This type of randomized preconditioning can greatly enhance the
robustness and power of "asynchronous" solvers for
communication-constrained environments [@2015_avron_revisiting]. Related
techniques for kernel ridge regression are described in . For linear
systems involving high-dimensional tensors, see .

For square linear systems, the search for randomized preconditioners has
been less fruitful.
Section [18](#sec:sparse-cholesky){reference-type="ref"
reference="sec:sparse-cholesky"} outlines the main success story.
Nevertheless, techniques already at hand can be very helpful for solving
linear systems in special situations, which we illustrate with a small
example.

Consider the task of solving
[\[eq:Ax=b\]](#eq:Ax=b){reference-type="eqref" reference="eq:Ax=b"} for
a positive-definite (PD) coefficient matrix $\bm{\mathsf{A}}$. In this
environment, the iterative method of choice is the conjugate gradient
(CG) algorithm [@1952_hestenes_CG]. A detailed convergence analysis for
CG is available; for example, see [@1997_trefethen_bau Sec. 38]. In a
nutshell, CG converges rapidly when the eigenvalues of $\bm{\mathsf{A}}$
are clustered, as in
Figure [7](#fig:eval_distributions_for_CG){reference-type="ref"
reference="fig:eval_distributions_for_CG"}(a). Therefore, our task is to
find a matrix $\bm{\mathsf{M}}$ for which $\bm{\mathsf{M}}^{-1}$ can be
applied rapidly to vectors and for which
$\bm{\mathsf{M}}^{-1/2} \bm{\mathsf{A}} \bm{\mathsf{M}}^{-1/2}$ has a
tightly clustered spectrum.

In a situation where $\bm{\mathsf{A}}$ has a few eigenvalues that are
larger than the others
(Figure [7](#fig:eval_distributions_for_CG){reference-type="ref"
reference="fig:eval_distributions_for_CG"}(b)), randomized algorithms
for low-rank approximation provide excellent preconditioners. For
instance, we can use the randomized Nyström method
(Section [14](#sec:nystrom){reference-type="ref"
reference="sec:nystrom"}) to compute an approximation $$\begin{equation}
\label{eq:spd_exact_lowrank}
\bm{\mathsf{A}} \approx \bm{\mathsf{U}}\bm{\mathsf{D}}\bm{\mathsf{U}}^{*}
\end{equation}$$ where $\bm{\mathsf{D}} \in \mathbb{R}_{+}^{k\times k}$
is a diagonal matrix whose entries hold approximations to the largest
$k$ eigenvalues of $\bm{\mathsf{A}}$, and where
$\bm{\mathsf{U}} \in \mathbb{F}^{m\times k}$ is an orthonormal matrix
holding the corresponding approximate eigenvectors. We then form a
preconditioner for $\bm{\mathsf{A}}$ by setting
$$\bm{\mathsf{M}} = (1/\alpha)\,\bm{\mathsf{U}}\bm{\mathsf{D}}\bm{\mathsf{U}}^{*} + \bigl(\bm{\mathsf{I}} - \bm{\mathsf{U}}\bm{\mathsf{U}}^{*}\bigr).$$
It is trivial to invert $\bm{\mathsf{M}}$ because
$\bm{\mathsf{M}}^{-1} = \alpha\bm{\mathsf{U}}\bm{\mathsf{D}}^{-1}\bm{\mathsf{U}}^{*} + \bigl(\bm{\mathsf{I}} - \bm{\mathsf{U}}\bm{\mathsf{U}}^{*}\bigr)$.
Now, if
[\[eq:spd_exact_lowrank\]](#eq:spd_exact_lowrank){reference-type="eqref"
reference="eq:spd_exact_lowrank"} captured the top $k$ eigenmodes of
$\bm{\mathsf{A}}$ exactly, then the preconditioned coefficient matrix
$\bm{\mathsf{M}}^{-1/2}\bm{\mathsf{A}}\bm{\mathsf{M}}^{-1/2}$ would have
the same eigenvectors as $\bm{\mathsf{A}}$, but with the top $k$
eigenvalues replaced by $\alpha$ and the remaining eigenvalues
unchanged. By setting $\alpha = \lambda_{k}$, say, the spectrum of
$\bm{\mathsf{M}}^{-1/2}\bm{\mathsf{A}}\bm{\mathsf{M}}^{-1/2}$ would
become far more tightly clustered. In reality, the columns of
$\bm{\mathsf{U}}$ do not exactly align with the eigenvectors of
$\bm{\mathsf{A}}$. Even so, the accuracy will be good for the
eigenvectors associated with the top eigenvalues, which is what matters.

<figure id="fig:eval_distributions_for_CG">
<table>
<tbody>
<tr>
<td style="text-align: left;"><embed src="Pics/fig_spectra_1.pdf"
style="width:30mm" /></td>
<td style="text-align: left;"><embed src="Pics/fig_spectra_2.pdf"
style="width:30mm" /></td>
<td style="text-align: left;"><embed src="Pics/fig_spectra_3.pdf"
style="width:30mm" /></td>
<td style="text-align: left;"><embed src="Pics/fig_spectra_4.pdf"
style="width:30mm" /></td>
</tr>
<tr>
<td style="text-align: left;">(a)</td>
<td style="text-align: left;">(b)</td>
<td style="text-align: left;">(c)</td>
<td style="text-align: left;">(d)</td>
</tr>
</tbody>
</table>
<figcaption>The eigenvalues of four different PD matrices that all have
condition number <span class="math inline">10</span> (since <span
class="math inline"><em>λ</em><sub>max</sub> = 1</span> and <span
class="math inline"><em>λ</em><sub>min</sub> = 0.1</span>). As discussed
in Section <a href="#sec:sketchtoprecond" data-reference-type="ref"
data-reference="sec:sketchtoprecond">17.3</a>, the difficulty of solving
the corresponding linear systems using conjugate gradients differ
significantly between these cases. (a) For this matrix, CG converges in
two iterations, without the need for preconditioners. (b) When the
spectrum has some large outliers, the randomized preconditioner outlined
in Section <a href="#sec:sketchtoprecond" data-reference-type="ref"
data-reference="sec:sketchtoprecond">17.3</a> works well. (c,d) Finding
randomized preconditioners for matrices with spectra like these remains
an open research problem.</figcaption>
</figure>

## The randomized Kaczmarz method and its relatives {#sec:rk}

The Kaczmarz method is an iterative algorithm for solving linear systems
that is typically used for large, overdetermined problems with
inconsistent equations. Randomized variants of the Kaczmarz method have
received a lot of attention in recent years, in part because of close
connections to stochastic gradient descent (SGD) algorithms for solving
least-squares problems.

To explain the idea, consider a (possibly inconsistent) linear system
$$\begin{equation}
 \label{eq:Ax=b_kacz}
\bm{\mathsf{A}}^* \bm{\mathsf{x}} \approx \bm{\mathsf{b}}
\quad\text{where}\quad
\bm{\mathsf{A}}^* \in \mathbb{F}^{m \times n}.
\end{equation}$$ The basic Kaczmarz algorithm starts with an initial
guess $\bm{\mathsf{x}}_0 \in \mathbb{F}^n$ for the solution. At each
iteration $t$, we select a new index $j = j(t) \in \{1, \dots, m\}$, and
we make the update $$\begin{equation}
 \label{eq:kacz_trad}
\bm{\mathsf{x}}_{t+1} = \bm{\mathsf{x}}_t + \frac{\bm{\mathsf{b}}(j) - \langle \bm{\mathsf{A}}(:, j), \, \bm{\mathsf{x}}_t \rangle}{\Vert \bm{\mathsf{A}}(:,j) \Vert^2} \bm{\mathsf{A}}(:, j).
\end{equation}$$ The
rule [\[eq:kacz_trad\]](#eq:kacz_trad){reference-type="eqref"
reference="eq:kacz_trad"} has a simple interpretation: it ensures that
$\bm{\mathsf{x}}_{t+1}$ is the closest point to $\bm{\mathsf{x}}_t$ in
the hyperplane containing solutions to the linear equation determined by
the $j$th equation in the system.

In implementing this method, we must choose a control mechanism that
determines the next index. A simple and robust approach is to cycle
through the rows consecutively; that is, $j(t) = t \bmod m$. Another
effective, but expensive, option is to select the equation with the
largest violation.

The randomized Kaczmarz (RK) algorithm uses a probabilistic control
mechanism instead. This kind of approach also has a long history, and it
is useful in cases where cyclic control is ineffective. The RK method
has received renewed attention owing to work of . They proposed sampling
each $j(t)$ independently at random, with the probability of choosing
the $i$th equation proportional to the squared $\ell_2$ norm of the
$i$th column of $\bm{\mathsf{A}}$. They proved that this version of RK
converges linearly with a rate determined by the (Demmel) condition
number of the matrix $\bm{\mathsf{A}}$. Later, it was recognized that
this approach is just a particular instantiation of SGD for the
least-squares
problem [\[eq:Ax=b_kacz\]](#eq:Ax=b_kacz){reference-type="eqref"
reference="eq:Ax=b_kacz"}. See , which draws on results from .

There are many subsequent papers that have built on the RK approach for
solving inconsistent linear systems. observed that related ideas can be
used to design randomized Gauss--Seidel and randomized Jacobi
iterations. studied a blocked version of the RK algorithm, which is
practically more efficient for many of the same reasons that other
blocked algorithms work well
(Section [16.2](#sec:blocking){reference-type="ref"
reference="sec:blocking"}).

observed that the RK algorithm is a particular type of iterative
sketching method. Based on this connection, they proposed a
generalization. At each iteration, we draw an independent random
embedding $\bm{\mathsf{S}}_t \in \mathbb{F}^{\ell \times m}$. The next
iterate is chosen by solving the least-squares problem
$$\begin{equation}
\label{eq:Ax=b_kacz_proj0}
\bm{\mathsf{x}}_{t} = \arg\min\nolimits_{\bm{\mathsf{y}}} \Vert  \bm{\mathsf{x}}_{t-1} - \bm{\mathsf{y}}  \Vert^2
\quad\text{subject to}\quad
\bm{\mathsf{S}}_t \bm{\mathsf{A}}^* \bm{\mathsf{y}} = \bm{\mathsf{S}}_t \bm{\mathsf{b}}.
\end{equation}$$ The idea is to choose the dimension $\ell$ sufficiently
small that the sketched least-squares problem can be solved explicitly
using a direct method (e.g., QR factorization). This flexibility leads
to algorithms that converge more rapidly in practice because the sketch
$\bm{\mathsf{S}}_t$ can mix equations instead of just sampling. Later,
showed that this procedure can be accelerated to achieve rates that
depend on the *square root* of an appropriate condition number; see also
.

