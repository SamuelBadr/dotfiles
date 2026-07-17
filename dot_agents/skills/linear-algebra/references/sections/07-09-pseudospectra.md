# 7.9 Pseudospectra

If the purpose of computing is insight, then it is easy to see why the well-conditioned eigenvector basis is such a valued commodity, for in many matrix problems, replacement of A with its diagonalization $X ^ { - 1 } A X$ leads to powerful, analytic simplifications. However, the insight-through-eigensystem paradigm has diminished impact in problems where the matrix of eigenvectors is ill-conditioned or nonexistent. Intelligent invariant subspace computation as discussed in §7.6 is one way to address the shortfall; pseudospectra are another. In this brief section we discuss the essential ideas behind the theory and computation of pseudospectra. The central message is simple: if you are working with a nonnormal matrix, then a graphical pseudospectral analysis effectively tells you just how much to trust the eigenvalue/eigenvector “story.”

A slightly awkward feature of our presentation has to do with the positioning of this section in the text. As we will see, SVD calculations are an essential part of the pseudospectra scene and we do not detail dense matrix algorithms for that important decomposition until the next chapter. However, it makes sense to introduce the pseudospectra concept here at the end of Chapter 7 while the challenges of the unsymmetric eigenvalue problem are fresh in mind. Moreover, with this “early” foundation we can subsequently present various pseudospectra insights that concern the behavior of the matrix exponential (§9.3), the Arnoldi method for sparse unsymmetric eigenvalue problems (§10.5), and the GMRES method for sparse unsymmetric linear systems (§11.4).

For maximum generality, we investigate the pseudospectra of complex, nonnormal matrices. The definitive pseudospectra reference is Trefethen and Embree (SAP). Virtually everything we discuss is presented in greater detail in that excellent volume.

# 7.9.1 Motivation

In many settings, the eigenvalues of a matrix “say something” about an underlying phenomenon. For example, if

$$
A = \left[ \begin{array}{c c} \lambda_ {1} & M \\ 0 & \lambda_ {2} \end{array} \right], \qquad M > 0,
$$

then

$$
\lim _ {k \to \infty} \| A ^ {k} \| _ {2} = 0
$$

if and only if $| \lambda _ { 1 } | ~ < ~ 1$ and $| \lambda _ { 2 } | ~ < ~ 1$ . This follows from Lemma 7.3.1, a result that we needed to establish the convergence of the QR iteration. Applied to our 2-by-2 example, the lemma can be used to show that

$$
\left\| A ^ {k} \right\| _ {2} \leq \frac {M}{\epsilon} (\rho (A) + \epsilon) ^ {k}
$$

for any $\epsilon > 0$ where $\rho ( A ) = \operatorname* { m a x } \{ | \lambda _ { 1 } | , | \lambda _ { 2 } | \}$ is the spectral radius. By making 
 small enough in this inequality, we can draw a conclusion about the asymptotic behavior of $A ^ { k }$ :

$$
\text {   If   } \rho (A) <   1, \text {   then   asymptotically   } A ^ {k} \text {   converges   to   zero   as   } \rho (A) ^ {k}. \tag {7.9.1}
$$

However, while the eigenvalues adequately predict the limiting behavior of $\| \ b { A } ^ { k } \| _ { 2 }$ , they do not (by themselves) tell us much about what is happening if k is small. Indeed, if $\lambda _ { 1 } \neq \lambda _ { 2 }$ , then using the diagonalization

$$
A = \left[ \begin{array}{c c} 1 & M / (\lambda_ {2} - \lambda_ {1}) \\ 0 & 1 \end{array} \right] \left[ \begin{array}{c c} \lambda_ {1} & 0 \\ 0 & \lambda_ {2} \end{array} \right] \left[ \begin{array}{c c} 1 & M / (\lambda_ {2} - \lambda_ {1}) \\ 0 & 1 \end{array} \right] ^ {- 1} \tag {7.9.2}
$$

we can show that

$$
A ^ {k} = \left[ \begin{array}{c c} \lambda_ {1} ^ {k} & M \sum_ {j = 0} ^ {k - 1} \lambda_ {1} ^ {k - 1 - j} \lambda_ {2} ^ {j} \\ \hline 0 & \lambda_ {2} ^ {k} \end{array} \right]. \tag {7.9.3}
$$

Consideration of the (1,2) entry suggests that $A ^ { k }$ may grow before decay sets in. This is affirmed in Figure 7.9.1 where the size of $\parallel A ^ { k } \parallel _ { 2 }$ is tracked for the example

$$
A = \left[ \begin{array}{c c} 0. 9 9 9 & 1 0 0 0 \\ 0. 0 & 0. 9 9 8 \end{array} \right].
$$

![](images/golub_450_499__64cff00864bd65fd006fdbe2008a35bfb0adeafb211c74efbf56bf4d477cde7d.jpg)  
Figure 7.9.1. $\| \ b { A } ^ { k } \| _ { 2 }$ can grow even if $\rho ( A ) < 1$

Thus, it is perhaps better to augment (7.9.1) as follows:

$$
\begin{array}{l} \text {   If   } \rho (A) <   1, \text {   then   asymptotically   } A ^ {k} \text {   converges   to   zero   like   } \rho (A) ^ {k}. \\ \text {   However,   } A ^ {k} \text {   may   grow   substantially   before   exponential   decay   sets   in.   } \end{array} \tag {7.9.4}
$$

This example with its ill-conditioned eigenvector matrix displayed in (7.9.2), points to just why classical eigenvalue analysis is not so informative for nonnormal matrices. Ill-conditioned eigenvector bases create a discrepancy between how A behaves and how its diagonalization $X A X ^ { - 1 }$ behaves. Pseudospectra analysis and computation narrow this gap.

# 7.9.2 Definitions

The pseudospectra idea is a generalization of the eigenvalue idea. Whereas the spectrum $\Lambda ( A )$ is the set of all $z \in \mathbb { C }$ that make $\sigma _ { m i n } ( A - \lambda I )$ zero, the 
-pseudospectrum of a matrix $A \in \mathbb { C } ^ { n \times n }$ is the subset of the complex plane defined by

$$
\Lambda_ {\epsilon} (A) = \{z \in \mathbb {C}: \sigma_ {\min} (A - \lambda I) \leq \epsilon \}. \tag {7.9.5}
$$

If $\lambda \in \Lambda _ { \epsilon } ( A )$ , then λ is an 
-pseudoeigenvalue of A. A unit 2-norm vector v that satisfies $\parallel ( A - \lambda I ) v \parallel _ { 2 } = \epsilon$ is a corresponding 
-pseudoeigenvector. Note that if 
 is zero, then $\Lambda _ { \epsilon } ( A )$ is just the set of A’s eigenvalues, i.e., $\Lambda _ { 0 } ( A ) = \Lambda ( A )$ .

We mention that because of their interest in what pseudospectra say about general linear operators, Trefethen and Embree (2005) use a strict inequality in the definition (7.9.5). The distinction has no impact in the matrix case.

Equivalent definitions of $\Lambda _ { \epsilon } ( \cdot )$ include

$$
\Lambda_ {\epsilon} (A) = \left\{z \in \mathbb {C}: \| (z I - A) ^ {- 1} \| _ {2} \geq \frac {1}{\epsilon} \right\} \tag {7.9.6}
$$

which highlights the resolvent $( z I - A ) ^ { - 1 }$ and

$$
\Lambda_ {\epsilon} (A) = \{z \in \mathbb {C}: z \in \Lambda (A + E), \| E \| _ {2} \leq \epsilon \} \tag {7.9.7}
$$

which characterize pseudspectra as (traditional) eigenvalues of nearby matrices. The equivalence of these three definitions is a straightforward verification that makes use of Chapter 2 facts about singular values, 2-norms, and matrix inverses. We mention that greater generality can be achieved in (7.9.6) and (7.9.7) by replacing the 2-norm with an arbitrary matrix norm.

# 7.9.3 Display

The pseudospectrum of a matrix is a visible subset of the complex plane so graphical display has a critical role to play in pseudospectra analysis. The Matlab-based Eigtool system developed by Wright(2002) can be used to produce pseudospectra plots that are as pleasing to the eye as they are informative. Eigtool’s pseudospectra plots are contour plots where each contour displays the z-values associated with a specified value of 
. Since

$$
\epsilon_ {1} \leq \epsilon_ {2} \quad \Rightarrow \quad \Lambda_ {\epsilon_ {1}} \subseteq \Lambda_ {\epsilon_ {2}}
$$

the typical pseudospectral plot is basically a topographical map that depicts the function $f ( z ) = \sigma _ { \operatorname* { m i n } } ( z I - A )$ in the vicinity of the eigenvalues.

We present three Eigtool-produced plots that serve as illuminating examples. The first involves the n-by-n Kahan matrix $\mathsf { K a h } _ { n } ( s ) , \mathsf { e . g . }$ .,

$$
\mathsf {K a h} _ {5} (s) = \left[ \begin{array}{l l l l l} 1 & - c & - c & - c & - c \\ 0 & s & - s c & - s c & - s c \\ 0 & 0 & s ^ {2} & - s ^ {2} c & - s ^ {2} c \\ 0 & 0 & 0 & s ^ {3} & - s ^ {3} c \\ 0 & 0 & 0 & 0 & s ^ {4} \end{array} \right], \qquad c ^ {2} + s ^ {2} = 1.
$$

Recall that we used these matrices in §5.4.3 to show that QR with column pivoting can fail to detect rank deficiency. The eigenvalues $\{ 1 , s , s ^ { 2 } , \ldots , s ^ { n - 1 } \}$ of $\mathsf { K a h } _ { n } ( s )$ are extremely sensitive to perturbation. This is revealed by considering the $\epsilon = 1 0 ^ { - 6 }$ contour that is displayed in Figure 7.9.2 together with $\Lambda ( \mathsf { K a h } _ { n } ( s ) )$ .

The second example is the Demmel matrix ${ \mathsf { D e m } } _ { n } ( \beta ) , { \mathsf { e . g . } }$ ,

$$
\mathsf {D e m} _ {5} (\beta) = - \left[ \begin{array}{l l l l l} 1 & \beta & \beta^ {2} & \beta^ {3} & \beta^ {4} \\ 0 & 1 & \beta & \beta^ {2} & \beta^ {3} \\ 0 & 0 & 1 & \beta & \beta^ {2} \\ 0 & 0 & 0 & 1 & \beta \\ 0 & 0 & 0 & 0 & 1 \end{array} \right].
$$

![](images/golub_450_499__79c19c3fb42559549d64a72fd03e87ece5e3b03f173e18bdb2bb76bc4ed0f4bb.jpg)

<details>
<summary>contour</summary>

| x    | y    |
| ---- | ---- |
| 0.0  | 0.0  |
| 0.1  | 0.05 |
| 0.2  | 0.1  |
| 0.3  | 0.15 |
| 0.4  | 0.2  |
| 0.5  | 0.25 |
| 0.6  | 0.3  |
| 0.7  | 0.35 |
| 0.8  | 0.4  |
| 0.9  | 0.45 |
| 1.0  | 0.5  |
</details>

Figure 7.9.2. Λ $\left( \mathsf { K a h } _ { 3 0 } \left( s \right) \right)$ with $s ^ { 2 9 } = 0 . 1$ and contours for $\epsilon = 1 0 ^ { - 2 } , \dots , 1 0 ^ { - 6 }$

The matrix $\mathrm { D e m } _ { n } ( \beta )$ is defective and has the property that very small perturbations can move an original eigenvalue to a position that are relatively far out on the imaginary axis. See Figure 7.9.3. The example is used to illuminate the nearness-to-instability problem presented in P7.9.13.

![](images/golub_450_499__a666206744e55fe2f955ed8b375ed5efb43c91d358243526585a90b185aa74d4.jpg)

<details>
<summary>contour</summary>

| x    | y    |
| ---- | ---- |
| -10  | 6    |
| -5   | 4    |
| 0    | 0    |
| 5    | -4   |
| 10   | -8   |
</details>

Figure 7.9.3. $\Lambda _ { \epsilon } ( \mathsf { D e m } _ { 5 0 } ( \beta ) )$ with $\beta ^ { 4 9 } = 1 0 ^ { 8 }$ and contours for $\epsilon = 1 0 ^ { - 2 } , \dots , 1 0 ^ { - 6 }$

The last example concerns the pseudospectra of the Matlab “Gallery(5)” matrix:

$$
G _ {5} = \left[ \begin{array}{r r r r r} - 9 & 1 1 & - 2 1 & 6 3 & - 2 5 2 \\ 7 0 & - 6 9 & 1 4 1 & - 4 2 1 & 1 6 8 4 \\ - 5 7 5 & 5 7 5 & - 1 1 4 9 & 3 4 5 1 & - 1 3 8 0 1 \\ 3 8 9 1 & - 3 8 9 1 & 7 7 8 2 & - 2 3 3 4 5 & 9 3 3 6 5 \\ 1 0 2 4 & - 1 0 2 4 & 2 0 4 8 & - 6 1 4 4 & 2 4 5 7 2 \end{array} \right].
$$

Notice in Figure 7.9.4 that $\Lambda _ { 1 0 ^ { - 1 3 . 5 } } ( G _ { 5 } )$ has five components. In general, it can be

![](images/golub_450_499__f635dbcd6e60007c42e81916c3ad5f9b36cdad48a69dc2e5677e3a7207fbbd71.jpg)

<details>
<summary>contour</summary>

| x       | y       |
| ------- | ------- |
| -0.02   | 0.02    |
| 0.01    | 0.01    |
| 0.02    | 0.01    |
| -0.01   | 0.00    |
| -0.03   | -0.01   |
| 0.00    | -0.02   |
</details>

Figure 7.9.4. $\Lambda _ { \epsilon } ( G _ { 5 } )$ with contours for $\epsilon = 1 0 ^ { - 1 1 . 5 } , 1 0 ^ { - 1 2 } , \allowbreak . \dots , 1 0 ^ { - 1 3 . 5 } , 1 0 ^ { - 1 4 }$

shown that each connected component of $\Lambda _ { \epsilon } ( A )$ contains at least one eigenvalue of A.

# 7.9.4 Some Elementary Properties

Pseudospectra are subsets of the complex plane so we start with a quick summary of notation. If $S _ { 1 }$ and $S _ { 2 }$ are subsets of the complex plane, then their sum $S _ { 1 } + S _ { 2 }$ is defined by

$$
S _ {1} + S _ {2} = \left\{s: s = s _ {1} + s _ {2}, s _ {1} \in S _ {1}, s _ {2} \in S _ {2} \right\}.
$$

If $S _ { 1 }$ consists of a single complex number α, then we write $\alpha + S _ { 2 }$ . If S is a subset of the complex plane and $\beta$ is a complex number, then $\beta \cdot S$ is defined by

$$
\beta \cdot S = \{\beta z: z \in S \}.
$$

The disk of radius 
 centered at the origin is denoted by

$$
\Delta_ {\epsilon} = \{z: | z | \leq \epsilon \}.
$$

Finally, the distance from a complex number $z _ { \mathrm { 0 } }$ to a set of complex numbers S is defined by

$$
\operatorname{dist} \left(z _ {0}, S\right) = \min \left\{\left| z _ {0} - z \right|: z \in S \right\}.
$$

Our first result is about the effect of translation and scaling. For eigenvalues we have

$$
\Lambda (\alpha I + \beta A) = \alpha + \beta \cdot \Lambda (A).
$$

The following theorem establishes an analogous result for pseudospectra.

Theorem 7.9.1. If $\alpha , \beta \in \mathbb { C }$ and $A \in \mathbb { C } ^ { n \times n }$ , then $\Lambda _ { \epsilon | \beta | } ( \alpha I + \beta A ) = \alpha + \beta \cdot \Lambda _ { \epsilon } ( A )$ .

Proof. Note that

$$
\begin{array}{l} \Lambda_ {\epsilon} (\alpha I + A) = \left\{z: \| (z I - (\alpha I + A)) ^ {- 1} \| \geq 1 / \epsilon \right\} \\ = \left\{z: \| ((z - \alpha) I - A) ^ {- 1} \| \geq 1 / \epsilon \right\} \\ = \alpha + \left\{z - \alpha : \| ((z - \alpha) I - A) ^ {- 1} \| \geq 1 / \epsilon \right\} \\ = \alpha + \left\{z: \| (z I - A) ^ {- 1} \| \geq 1 / \epsilon \right\} = \Lambda_ {\epsilon} (A) \\ \end{array}
$$

and

$$
\begin{array}{l} \Lambda_ {\epsilon | \beta |} (\beta \cdot A) = \left\{z: \| (z I - \beta A) ^ {- 1} \| \geq 1 / | \beta | \epsilon \right\} \\ = \left\{z: \parallel (z / \beta) I - A) ^ {- 1} \parallel \geq 1 / \epsilon \right\} \\ = \beta \cdot \left\{z / \beta : \| (z / \beta) I - A) ^ {- 1} \| \geq 1 / \epsilon \right\} \\ = \beta \cdot \left\{z: \| z I - A) ^ {- 1} \| \geq 1 / \epsilon \right\} = \beta \cdot \Lambda_ {\epsilon} (A). \\ \end{array}
$$

The theorem readily follows by composing these two results.

General similarity transforms preserve eigenvalues but not 
-pseudoeigenvalues. However, a simple inclusion property holds in the pseudospectra case.

Theorem 7.9.2. If $B = X ^ { - 1 } A X$ , then $\Lambda _ { \epsilon } ( B ) \subseteq \Lambda _ { \epsilon \kappa _ { 2 } ( X ) } ( A )$ .

Proof. $\mathrm { I f } ~ z \in \Lambda _ { \epsilon } ( B )$ , then

$$
\frac {1}{\epsilon} \leq \| (z I - B) ^ {- 1} \| = \| X ^ {- 1} (z I - A) ^ {- 1} X ^ {- 1} \| \leq \kappa_ {2} (X) \| (z I - A) ^ {- 1} \|,
$$

from which the theorem follows.

Corollary 7.9.3. If $X \in \mathbb { C } ^ { n \times n }$ is unitary and $A \in \mathbb { C } ^ { n \times n }$ , then $\Lambda _ { \epsilon } ( X ^ { - 1 } A X ) = \Lambda _ { \epsilon } ( A )$ .

Proof. The proof is left as an exercise.

The 
-pseudospectrum of a diagonal matrix is the union of 
-disks.

Theorem 7.9.4. If $D = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ , then $\Lambda _ { \epsilon } ( D ) = \{ \lambda _ { 1 } , \dots , \lambda _ { n } \} + \Delta _ { \epsilon }$ .

Proof. The proof is left as an exercise.

Corollary 7.9.5. If $A \in \mathbb { C } ^ { n \times n }$ is normal, then $\Lambda _ { \epsilon } ( A ) = \Lambda ( A ) + \Delta _ { \epsilon }$ .

Proof. Since A is normal, it has a diagonal Schur form $Q ^ { H } A Q = \operatorname { d i a g } ( \lambda _ { 1 } , \ldots , \lambda _ { n } ) = D$ with unitary Q. The proof follows from Theorem 7.9.4.

If $T = \left( T _ { i j } \right)$ is a 2-by-2 block triangular matrix, then $\Lambda ( T ) = \Lambda ( T _ { 1 1 } ) \cup \Lambda ( T _ { 2 2 } )$ . Here is the pseudospectral analog:

Theorem 7.9.6. If

$$
T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right]
$$

with square diagonal blocks, then $\Lambda _ { \epsilon } ( T _ { 1 1 } ) \cup \Lambda _ { \epsilon } ( T _ { 2 2 } ) \subseteq \Lambda _ { \epsilon } ( T )$ .

Proof. The proof is left as an exercise.

Corollary 7.9.7. If

$$
T = \left[ \begin{array}{c c} T _ {1 1} & 0 \\ 0 & T _ {2 2} \end{array} \right]
$$

with square diagonal blocks, then $\Lambda _ { \epsilon } ( T ) = \Lambda _ { \epsilon } ( T _ { 1 1 } ) \cup \Lambda _ { \epsilon } ( T _ { 2 2 } )$ .

Proof. The proof is left as an exercise.

The last property in our gallery of facts connects the resolvant $( z _ { 0 } I - A ) ^ { - 1 }$ to the distance that separates $z _ { \mathrm { 0 } }$ from $\Lambda _ { \epsilon } ( A )$ .

Theorem 7.9.8. If $z _ { 0 } \in \mathbb { C }$ and $A \in \mathbb { C } ^ { n \times n }$ , then

$$
\operatorname{dist} \left(z _ {0}, \Lambda_ {\epsilon} (A)\right) \geq \frac {1}{\| \left(z _ {0} I - A\right) ^ {- 1} \| _ {2}} - \epsilon .
$$

Proof. For any $z \in \Lambda _ { \epsilon } ( A )$ we have from Corollary 2.4.4 and (7.9.6) that

$$
\epsilon \geq \sigma_ {\min} (z I - A) = \sigma_ {\min} ((z _ {0} I - A) - (z - z _ {0}) I) \geq \sigma_ {\min} (z _ {0} I - A) - | z - z _ {0} |
$$

and thus

$$
| z - z _ {0} | \geq \frac {1}{\| (z _ {0} I - A) ^ {- 1} \|} - \epsilon .
$$

The proof is completed by minimizing over all $z \in \Lambda _ { \epsilon ( A ) }$ .

# 7.9.5 Computing Pseudospectra

The production of a pseudospectral contour plot such as those displayed above requires sufficiently accurate approximations of $\sigma _ { \operatorname* { m i n } } ( z I - A )$ on a grid that consists of (perhaps)

1000’s of z-values. As we will see in §8.6, the computation of the complete SVD of an n-by-n dense matrix is an $O ( n ^ { 3 } )$ endeavor. Fortunately, steps can be taken to reduce each grid point calculation to ${ \dot { O ( n ^ { 2 } ) } }$ or less by exploiting the following ideas:

1. Avoid SVD-type computations in regions where $\sigma _ { \operatorname* { m i n } } ( z I - A )$ is slowly varying. See Gallestey (1998).

2. Exploit Theorem 7.9.6 by ordering the eigenvalues so that the invariant subspace associated with $\Lambda ( T _ { 1 1 } )$ captures the essential behavior of $( z I - A ) ^ { - 1 }$ . See Reddy, Schmid, and Henningson (1993).

3. Precompute the Schur decomposition $Q ^ { H } A Q = T$ and apply a $\sigma _ { \mathrm { m i n } }$ algorithm that is efficient for triangular matrices. See Lui (1997).

We offer a few comments on the last strategy since it has much in common with the condition estimation problem that we discussed in §3.5.4. The starting point is to recognize that since Q is unitary,

$$
\sigma_ {\min} (z I - A) = \sigma_ {\min} (z I - T).
$$

The triangular structure of the transformed problem makes it possible to obtain a satisfactory estimate of $\sigma _ { \operatorname* { m i n } } ( z I - A )$ in $O ( n ^ { 2 } )$ flops. If d is a unit 2-norm vector and $( z I - T ) y = d { \mathrm { . } }$ , then it follows from the SVD of $z I - T$ that

$$
\sigma_ {\min} (z I - T) \leq \frac {1}{\| y \| _ {2}}.
$$

Let $u _ { \mathrm { m i n } }$ be a left singular vector associated with $\sigma _ { \operatorname* { m i n } } ( z I - T )$ . If d is has a significant component in the direction of $u _ { \mathrm { m i n } }$ , then

$$
\sigma_ {\min} (z I - T) \approx \frac {1}{\| y \| _ {2}}.
$$

Recall that Algorithm 3.5.1 is a cheap heuristic procedure that dynamically determines the right hand side vector d so that the solution to a given triangular system is large in norm. This is tantamount to choosing d so that it is rich in the direction of $u _ { \mathrm { m i n } }$ . A complex arithmetic, 2-norm variant of Algorithm 3.5.1 is outlined in P7.9.13. It can be applied to $z I - T$ . The resulting d-vector can be refined using inverse iteration ideas, see Toh and Trefethen (1996) and §8.2.2. Other approaches are discussed by Wright and Trefethen (2001).

# 7.9.6 Computing the -Pseudospectral Abscissa and Radius

The 
-pseudospectral abscissa of a matrix $A \in \mathbb { C } ^ { n \times n }$ is the rightmost point on the boundary of $\Lambda _ { \epsilon }$ :

$$
\alpha_ {\epsilon} (A) = \max \operatorname{Re} (z). \tag {7.9.8}
$$

$$
z \in \Lambda_ {\epsilon} (A)
$$

Likewise, the 
-pseudospectral radius is the point of largest magnitude on the boundary of $\Lambda _ { \epsilon }$ :

$$
\rho_ {\epsilon} (A) = \max _ {z \in \Lambda_ {\epsilon} (A)} | z |. \tag {7.9.9}
$$

These quantities arise in the analysis of dynamical systems and effective iterative algorithms for their estimation have been proposed by Burke, Lewis, and Overton (2003) and Mengi and Overton (2005). A complete presentation and analysis of their very clever optimization procedures, which build on the work of Byers (1988), is beyond the scope of the text. However, at their core they involve interesting intersection problems that can be reformulated as structured eigenvalue problems. For example, if i·r is an eigenvalue of the matrix

$$
M = \left[ \begin{array}{c c} i e ^ {i \theta} A ^ {H} & - \epsilon I \\ \epsilon I & i e ^ {- i \theta} A \end{array} \right], \tag {7.9.10}
$$

then 
 is a singular value of $A - r e ^ { i \theta } I$ . To see this, observe that if

$$
\left[ \begin{array}{c c} i e ^ {i \theta} A ^ {H} & - \epsilon I \\ \epsilon I & i e ^ {- i \theta} A \end{array} \right] \left[ \begin{array}{c} f \\ g \end{array} \right] = i \cdot r \left[ \begin{array}{c} f \\ g \end{array} \right],
$$

then

$$
(A - r e ^ {i \theta} I) ^ {H} (A - r e ^ {i \theta} I) g = \epsilon^ {2} g.
$$

The complex version of the SVD (§2.4.4) says that 
 is a singular value of $A - r e ^ { 1 \theta } I$ . It can be shown that if $i r _ { \mathrm { m a x } }$ is the largest pure imaginary eigenvalue of M, then

$$
\epsilon = \sigma_ {\mathrm{min}} (A - r _ {\mathrm{max}} e ^ {1 \theta} I).
$$

This result can be used to compute the intersection of the ray $\lbrace r e ^ { i \theta } : R \geq 0 \rbrace$ and the boundary of $\Lambda _ { \epsilon } ( A )$ . This computation is at the heart of computing the 
-pseudospectral radius. See Mengi and Overton (2005).

# 7.9.7 Matrix Powers and the -Pseudospectral Radius

At the start of this section we used the example

$$
A = \left[ \begin{array}{l l} 0. 9 9 9 & 1 0 0 0 \\ 0. 0 0 0 & 0. 9 9 8 \end{array} \right]
$$

to show that $\parallel A ^ { k } \parallel _ { 2 }$ can grow even though $\rho ( A ) < 1$ . This kind of transient behavior can be anticipated by the pseudospectral radius. Indeed, it can be shown that for any $\epsilon > 0$ ,

$$
\sup _ {k \geq 0} \| A ^ {k} \| _ {2} \geq \frac {\rho_ {\epsilon} (A) - 1}{\epsilon}. \tag {7.9.11}
$$

See Trefethen and Embree (SAP, pp. 160–161). This says that transient growth will occur if there is a contour $\left\{ z [ \right| ( \mathrm {  { \left| \left| \left| \boldsymbol { z } \right.  \right\right\} } A ) ^ { - 1 } = \mathrm { 1 } / \epsilon. }$ that extends beyond the unit disk. For the above 2-by-2 example, if $\epsilon = 1 0 ^ { - 8 }$ , then $\rho _ { \epsilon } ( A ) \approx 1 . 0 0 1 7$ and the inequality (7.9.11) says that for some k,  $A ^ { k } \parallel _ { 2 } \geq 1 . 7 \times 1 0 ^ { 5 }$ . This is consistent with what is displayed in Figure 7.9.1.

# Problems

P7.9.1 Show that the definitions (7.9.5), (7.9.6), and (7.9.7) are equivalent.

P7.9.2 Prove Corollary 7.9.3.

P7.9.3 Prove Theorem 7.9.4.

P7.9.4 Prove Theorem 7.9.6.

P7.9.5 Prove Corollary 7.9.7.

P7.9.6 Show that if $A , E \in \mathbb { C } ^ { n \times n }$ , then $\Lambda _ { \epsilon } ( A + E ) \subseteq \Lambda _ { \epsilon + \| E \| _ { 2 } } ( A )$

P7.9.7 Suppose $\sigma _ { \mathrm { m i n } } ( z _ { 1 } I - A ) = \epsilon _ { 1 }$ and $\sigma _ { \mathrm { m i n } } ( z _ { 2 } I - A ) = \epsilon _ { 2 }$ . Prove that there exists a real number µ so that if $z _ { 3 } = ( 1 - \mu ) z _ { 1 } + \mu z _ { 2 }$ , then $\sigma _ { \mathrm { m i n } } ( z _ { 3 } I - A ) = ( \epsilon _ { 1 } + \epsilon _ { 2 } ) / 2 ?$

P7.9.8 Suppose $A \in \mathbb { C } ^ { n \times n }$ is normal and $E \in \mathbb { C } ^ { n \times n }$ is nonnormal. State and prove a theorem about $\Lambda _ { \epsilon } ( A + E )$ .

P7.9.9 Explain the connection between Theorem 7.9.2 and the Bauer-Fike Theorem (Theorem 7.2.2).

P7.9.10 Define the matrix $J \in \mathbb { R } ^ { 2 n \times 2 n }$ by

$$
J = \left[ \begin{array}{c c} 0 & I _ {n} \\ - I _ {n} & 0 \end{array} \right].
$$

(a) The matrix $H \in \mathbb { R } ^ { 2 n \times 2 n }$ is a Hamiltonian matrix if $J ^ { T } H J = - H ^ { T }$ . It is easy to show that if H is Hamiltonian and $\lambda \in \Lambda ( H )$ , then $\cdot \lambda \in \Lambda ( H )$ . Does it follow that if $\lambda \in \Lambda _ { \epsilon } ( H )$ , then $\begin{array} { r } { \cdot \lambda \in \Lambda _ { \epsilon } ( H ) ? } \end{array}$ (b) The matrix $S \in \mathbb { R } ^ { 2 n \times 2 n }$ is a symplectic matrix if $J ^ { T } S J = S ^ { - T }$ . It is easy to show that if S is symplectic and $\lambda \in \Lambda ( S )$ , then $1 / \lambda \in \Lambda ( S )$ . Does it follow that if $\lambda \in \Lambda _ { \epsilon } ( S )$ , then $1 / \lambda \in \Lambda _ { \epsilon } ( S ) ?$

P7.9.11 Unsymmetric Toeplitz matrices tend to have very ill-conditioned eigensystems and thus have interesting pseudospectral properties. Suppose

$$
A = \left[ \begin{array}{c c c c} 0 & 1 & \dots & 0 \\ \alpha & 0 & \ddots & \vdots \\ \vdots & \ddots & \ddots & 1 \\ 0 & \dots & \alpha & 0 \end{array} \right].
$$

(a) Construct a diagonal matrix S so that $S ^ { - 1 } A S = B$ is symmetric and tridiagonal with 1’s on its subdiagonal and superdiagonal. (b) What can you say about the condition of A’s eigenvector matrix?

P7.9.12 A matrix $A \in \mathbb { C } ^ { n \times n }$ is stable if all of its eigenvalues have negative real parts. Consider the problem of minimizing $\parallel E \parallel _ { 2 }$ subject to the constraint that $A + E$ has an eigenvalue on the imaginary axis. Explain why this optimization problem is equivalent to minimizing $\sigma _ { \operatorname* { m i n } } ( i r I - A )$ over all $r \in \mathbb { R }$ . If $E _ { * }$ is a minimizing E, then $\parallel E \parallel _ { 2 }$ can be regarded as measure of A’s nearness to instability. What is the connection between A’s nearness to instability and $\alpha _ { \epsilon } ( A ) ?$

P7.9.13 This problem is about the cheap estimation of the minimum singular value of a matrix, a critical computation that is performed over an over again during the course of displaying the pseudospectrum of a matrix. In light of the discussion in 7.9.5, the challenge is to estimate the smallest singular value of an upper triangular matrix $U = T - z I$ where T is the Schur form of $A \in \mathbb { R } ^ { n \times n }$ . The condition estimation ideas of §3.5.4 are relevant. We want to determine a unit 2-norm vector $d \in \mathbb { C } ^ { n }$ such that the solution to $U y = d$ has a large 2-norm for then $\sigma _ { \mathrm { m i n } } ( U ) \approx 1 / \Vert \ y \Vert _ { 2 } .$ (a) Suppose

$$
U = \left[ \begin{array}{c c} u _ {1 1} & u ^ {H} \\ 0 & U _ {1} \end{array} \right] \qquad y = \left[ \begin{array}{c} \tau \\ z \end{array} \right] \qquad d = \left[ \begin{array}{c} c \\ s d _ {1} \end{array} \right]
$$

where $\begin{array} { r } { u _ { 1 1 } , \tau \in \mathbb { C } , u , z , d _ { 1 } \in \mathbb { C } ^ { n - 1 } , U _ { 1 } \in \mathbb { C } ^ { ( n - 1 ) \times ( n - 1 ) } , \parallel d _ { 1 } \parallel _ { 2 } = 1 , U _ { 1 } y _ { 1 } = d _ { 1 } } \end{array}$ , and $c ^ { 2 } + s ^ { 2 } = 1$ Give an algorithm that determines c and s so that if $U y = d ,$ then $\parallel y \parallel _ { 2 }$ is as large as possible. Hint: This is a 2-by-2 SVD problem. (b) Using part (a), develop a nonrecursive method for estimating $\sigma _ { \operatorname* { m i n } } ( U ( k { : } n , k { : } n ) )$ for $k = n \colon - 1 \colon 1$ .

# Notes and References for $\mathrm { 8 7 . 7 }$

Besides Trefethen and Embree (SAP), the following papers provide a nice introduction to the pseudospectra idea:

M. Embree and L.N. Trefethen (2001). “Generalizing Eigenvalue Theorems to Pseudospectra Theorems,” SIAM J. Sci. Comput. 23, 583–590.   
L.N. Trefethen (1997). “Pseudospectra of Linear Operators,” SIAM Review 39, 383–406.   
For more details concerning the computation and display of pseudoeigenvalues, see:   
S.C. Reddy, P.J. Schmid, and D.S. Henningson (1993). “Pseudospectra of the Orr-Sommerfeld Operator,” SIAM J. Applic. Math. 53, 15–47.   
S.-H. Lui (1997). “Computation of Pseudospectra by Continuation,” SIAM J. Sci. Comput. 18, 565–573.   
E. Gallestey (1998). “Computing Spectral Value Sets Using the Subharmonicity of the Norm of Rational Matrices,” BIT, 38, 22–33.   
L.N. Trefethen (1999). “Computation of Pseudospectra,” Acta Numerica 8, 247–295.   
T.G. Wright (2002). Eigtool, http://www.comlab.ox.ac.uk/pseudospectra/eigtool/.   
Interesting extensions/generalizations/applications of the pseudospectra idea include:   
L. Reichel and L.N. Trefethen (1992). “Eigenvalues and Pseudo-Eigenvalues of Toeplitz Matrices,” Lin. Alg. Applic. 164–164, 153–185.   
K-C. Toh and L.N. Trefethen (1994). “Pseudozeros of Polynomials and Pseudospectra of Companion Matrices,” Numer. Math. 68, 403–425.   
F. Kittaneh (1995). “Singular Values of Companion Matrices and Bounds on Zeros of Polynomials,” SIAM J. Matrix Anal. Applic. 16, 333–340.   
N.J. Higham and F. Tisseur (2000). “A Block Algorithm for Matrix 1-Norm Estimation, with an Application to 1-Norm Pseudospectra,” SIAM J. Matrix Anal. Applic. 21, 1185–1201.   
T.G. Wright and L.N. Trefethen (2002). “Pseudospectra of Rectangular matrices,” IMA J. Numer. Anal. 22, 501–519.   
R. Alam and S. Bora (2005). “On Stable Eigendecompositions of Matrices,” SIAM J. Matrix Anal. Applic. 26, 830–848.   
Pseudospectra papers that relate to the notions of controllability and stability of linear systems include:   
J.V. Burke and A.S. Lewis. and M.L. Overton (2003). “Optimization and Pseudospectra, with Applications to Robust Stability,” SIAM J. Matrix Anal. Applic. 25, 80–104.   
J.V. Burke, A.S. Lewis, and M.L. Overton (2003). “Robust Stability and a Criss–Cross Algorithm for Pseudospectra,” IMA J. Numer. Anal. 23, 359–375.   
J.V. Burke, A.S. Lewis and M.L. Overton (2004). “Pseudospectral Components and the Distance to Uncontrollability,” SIAM J. Matrix Anal. Applic. 26, 350–361.   
The following papers are concerned with the computation of the numerical radius, spectral radius, and field of values:   
C. He and G.A. Watson (1997). “An Algorithm for Computing the Numerical Radius,” IMA J. Numer. Anal. 17, 329–342.   
G.A. Watson (1996). “Computing the Numerical Radius” Lin. Alg. Applic. 234, 163–172.   
T. Braconnier and N.J. Higham (1996). “Computing the Field of Values and Pseudospectra Using the Lanczos Method with Continuation,” BIT 36, 422–440.   
E. Mengi and M.L. Overton (2005). “Algorithms for the Computation of the Pseudospectral Radius and the Numerical Radius of a Matrix,” IMA J. Numer. Anal. 25, 648–669.   
N. Guglielmi and M. Overton (2011). “Fast Algorithms for the Approximation of the Pseudospectral Abscissa and Pseudospectral Radius of a Matrix,” SIAM J. Matrix Anal. Applic. 32, 1166–1192.   
For more insight into the behavior of matrix powers, see:   
P. Henrici (1962). “Bounds for Iterates, Inverses, Spectral Variation, and Fields of Values of Nonnormal Matrices,” Numer. Math.4, 24–40.   
J. Descloux (1963). “Bounds for the Spectral Norm of Functions of Matrices,” Numer. Math. 5, 185–90.   
T. Ransford (2007). “On Pseudospectra and Power Growth,” SIAM J. Matrix Anal. Applic. 29, 699–711.   
As an example of what pseudospectra can tell us about highly structured matrices, see:   
L. Reichel and L.N. Trefethen (1992). “Eigenvalues and Pseudo-eigenvalues of Toeplitz Matrices,” Lin. Alg. Applic. 162/163/164, 153–186.

This page intentionally left blank
