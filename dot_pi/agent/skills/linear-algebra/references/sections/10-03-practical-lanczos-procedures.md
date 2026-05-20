# 10.3 Practical Lanczos Procedures

Rounding errors greatly affect the behavior of the Lanczos iteration. The basic difficulty is caused by loss of orthogonality among the Lanczos vectors, a phenomenon that muddies the issue of termination and complicates the relationship between A’s eigenvalues and those of the tridiagonal matrices $T _ { k }$ . This troublesome feature, coupled with the advent of Householder’s perfectly stable method of tridiagonalization, explains why the Lanczos algorithm was disregarded by numerical analysts during the 1950’s and 1960’s. However, the pressure to solve large, sparse eigenproblems coupled with the computational insights set forth by Paige (1971) changed all that. With many fewer than n iterations typically required to get good approximate extremal eigenvalues, the Lanczos method became attractive as a sparse matrix technique rather than as a competitor of the Householder approach.

Successful implementation of the Lanczos iteration involves much more than a simple encoding of Algorithm 10.1.1. In this section we present some of the ideas that have been proposed to make the Lanczos procedure viable in practice.

# 10.3.1 Required Storage and Work

With careful overwriting in Algorithm 10.1.1 and exploitation of the formula

$$
\alpha_ {k} = q _ {k} ^ {T} (A q _ {k} - \beta_ {k - 1} q _ {k - 1}),
$$

the whole Lanczos process can be implemented with just a pair of n-vectors:

$$
w = q _ {1}, v = A w, \alpha_ {1} = w ^ {T} v, v = v - \alpha_ {1} w, \beta_ {1} = \| v \| _ {2}, k = 1
$$

while $\beta _ { k } \neq 0$

for i = 1:n

$$
t = w _ {i}, w _ {i} = v _ {i} / \beta_ {k}, v _ {i} = - \beta_ {k} t
$$

end (10.3.1)

$$
v = v + A w
$$

$$
k = k + 1, \alpha_ {k} = w ^ {T} v, v = v - \alpha_ {k} w, \beta_ {k} = \parallel v \parallel_ {2}
$$

end

At the end of the loop body, the array w houses $q _ { k }$ and v houses the residual vector $r _ { k } = A q _ { k } - \alpha _ { k } q _ { k } - \beta _ { k - 1 } q _ { k - 1 }$ . See Paige (1972) for a discussion of various Lanczos implementations and their numerical properties. Note that A is not modified during the entire process and that is what makes the procedure so useful for large sparse matrices.

If A has an average of ν nonzeros per row, then approximately $( 2 \nu + 8 ) n$ flops are involved in a single Lanczos step. Upon termination the eigenvalues of $T _ { k }$ can be found using the symmetric tridiagonal QR algorithm or any of the special methods of §8.5 such as bisection. The Lanczos vectors are generated in the n-vector w. If eigenvectors are required, then the Lanczos vectors must be saved. Typically, they are stored in secondary memory units.

# 10.3.2 Roundoff Properties

The development of a practical, easy-to-use Lanczos tridiagonalization process requires an appreciation of the fundamental error analyses of Paige (1971, 1976, 1980). An examination of his results is the best way to motivate the several modified Lanczos procedures of this section.

After j steps of the iteration we obtain the matrix of computed Lanczos vectors $\hat { Q } _ { k } = \left[ \begin{array} { l } { \hat { q } _ { 1 } } \end{array} \right| \cdot \cdot \cdot \left| \begin{array} { l } { \hat { q } _ { k } } \end{array} \right]$ and the associated tridiagonal matrix

$$
\hat {T} _ {k} = \left[ \begin{array}{c c c c c} \hat {\alpha} _ {1} & \hat {\beta} _ {1} & & \dots & 0 \\ \hat {\beta} _ {1} & \hat {\alpha} _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \hat {\beta} _ {k - 1} \\ 0 & \dots & & \hat {\beta} _ {k - 1} & \hat {\alpha} _ {k} \end{array} \right].
$$

Paige (1971, 1976) shows that if $\hat { r } _ { k }$ is the computed analog of $r _ { k }$ , then

$$
A \hat {Q} _ {k} = \hat {Q} _ {k} \hat {T} _ {k} + \hat {r} _ {k} e _ {k} ^ {T} + E _ {k} \tag {10.3.2}
$$

where

$$
\| E _ {k} \| _ {2} \approx \mathbf {u} \| A \| _ {2}. \tag {10.3.3}
$$

This shows that the equation $A Q _ { k } = Q _ { k } T _ { k } + r _ { k } e _ { k } ^ { T }$ is satisfied to working precision.

Unfortunately, the picture is much less rosy with respect to the orthogonality among the $\hat { q } _ { i }$ . (Normality is not an issue. The computed Lanczos vectors essentially have unit length.) If ${ \hat { \beta } } _ { k } = \mathsf { f l } ( \parallel { \hat { r } } _ { k } \parallel _ { 2 } )$ and we compute $\hat { q } _ { k + 1 } = \mathsf { f l } \left( \hat { r } _ { k } / \hat { \beta } _ { k } \right)$ , then a simple analysis shows that

$$
\hat {\beta} _ {k} \hat {q} _ {k + 1} \approx \hat {r} _ {k} + w _ {k}
$$

where

$$
\| w _ {k} \| _ {2} \approx \mathbf {u} \| \hat {r} _ {k} \| _ {2} \approx \mathbf {u} \| A \| _ {2}.
$$

Thus, we may conclude that

$$
| \hat {q} _ {k + 1} ^ {T} \hat {q} _ {i} | \approx \frac {| \hat {r} _ {k} ^ {T} \hat {q} _ {i} | + \mathbf {u} | | A | | _ {2}}{| \hat {\beta} _ {k} |}
$$

for $i = 1 { : } k$ . In other words, significant departures from orthogonality can be expected when $\hat { \beta } _ { k }$ is small, even in the ideal situation where $\hat { r } _ { k } ^ { T } \hat { Q } _ { k }$ is zero. A small $\hat { \beta } _ { k }$ implies cancellation in the computation of $\hat { r } _ { k }$ . We stress that loss of orthogonality is due to one or several such cancellations and is not the result of the gradual accumulation of roundoff error.

Further details of the Paige analysis are given shortly. Suffice it to say now that loss of orthogonality always occurs in practice and with it, an apparent deterioration in the quality of $\hat { T } _ { k } ^ { \phantom { } } \mathrm { ' s }$ eigenvalues. This can be quantified by combining (10.3.2) with Theorem 8.1.16. In particular, if we set

$$
F _ {1} = \hat {r} _ {k} e _ {k} ^ {T} + E _ {k}, \qquad X _ {1} = \hat {Q} _ {k}, \qquad S = \hat {T} _ {k},
$$

in that theorem and assume that

$$
\tau = \parallel \hat {Q} _ {k} ^ {T} \hat {Q} _ {k} - I _ {k} \parallel_ {2}
$$

satisfies $\tau < 1$ , then there exist eigenvalues $\mu _ { 1 } , \dots , \mu _ { k } \in \lambda ( A )$ such that

$$
| \mu_ {i} - \lambda_ {i} (T _ {k}) | \leq \sqrt {2} \left(\| \hat {r} _ {k} \| _ {2} + \| E _ {k} \| _ {2} + \tau (2 + \tau) \| A \| _ {2}\right)
$$

for $i = 1 { : } k$ . An obvious way to control the τ factor is to orthogonalize each newly computed Lanczos vector against its predecessors. This leads directly to our first “practical” Lanczos procedure.

# 10.3.3 Lanczos with Complete Reorthogonalization

Let $r _ { 0 } , \ldots , r _ { k - 1 } \in \mathbb { R } ^ { n }$ be given and suppose that Householder matrices $H _ { 0 } , \ldots , H _ { k - 1 }$ have been computed such that $( H _ { 0 } \cdot \cdot \cdot H _ { k - 1 } ) ^ { T } \left[ \begin{array} { l } { r _ { 0 } } \end{array} \right| \cdot \cdot \cdot \left| \begin{array} { l } { r _ { k - 1 } } \end{array} \right]$ is upper triangular. Let ${ \left[ \begin{array} { l } { q _ { 1 } } \end{array} \right| } \cdots { \left| \begin{array} { l } { q _ { k } } \end{array} \right] }$ denote the first k columns of the Householder product $\left( H _ { 0 } \cdot \cdot \cdot H _ { k - 1 } \right)$ . Now suppose that we are given a vector $r _ { k } \in \mathbb { R } ^ { n }$ and wish to compute a unit vector $q _ { k + 1 }$ in the direction of

$$
w = r _ {k} - \sum_ {i = 1} ^ {k} \left(q _ {i} ^ {T} r _ {k}\right) q _ {i} \in \operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\} ^ {\perp}.
$$

If a Householder matrix $H _ { k }$ is determined so $( H _ { 0 } \cdot \cdot \cdot H _ { k } ) ^ { T } \left[ \begin{array} { l } { r _ { 0 } } \end{array} | \cdot \cdot \cdot | \ r _ { k } \right]$ is upper triangular, then it follows that column $( k + 1 )$ of $H _ { 0 } \cdots H _ { k }$ is the desired unit vector.

If we incorporate these Householder computations into the Lanczos process, then we can produce Lanczos vectors that are orthogonal to machine precision:

$r _ { 0 } = q _ { 1 }$ (given unit vector)

Determine Householder $H _ { 0 }$ so $H _ { 0 } r _ { 0 } = e _ { 1 }$ .

for $k = 1 { : } n - 1$

$$
\begin{array}{l} \alpha_ {k} = q _ {k} ^ {T} A q _ {k} \\ r _ {k} = (A - \alpha_ {k} I) q _ {k} - \beta_ {k - 1} q _ {k - 1}, \quad (\beta_ {0} q _ {0} \equiv 0) \tag {10.3.4} \\ \end{array}
$$

$$
w = (H _ {k - 1} \dots H _ {0}) r _ {k}
$$

$\mathrm { D e t e r m i n e ~ H o u s e h o l d e r ~ } H _ { k } \mathrm { ~ s o ~ } H _ { k } w = [ w _ { 1 } , \dots , w _ { k } , \beta _ { k } , 0 , \dots , 0 ] ^ { T } .$

$$
q _ {k + 1} = H _ {0} \dots H _ {k} e _ {k + 1}
$$

end

This is an example of a complete reorthorgonalization Lanczos scheme. The idea of using Householder matrices to enforce orthogonality appears in Golub, Underwood, and Wilkinson (1972). That the computed $\hat { q } _ { i }$ in (10.3.4) are orthogonal to working precision follows from the roundoff properties of Householder matrices. Note that by virtue of the definition of $q _ { k + 1 }$ , it makes no difference if $\beta _ { k } = 0$ . For this reason, the algorithm may safely run until $k = n - 1$ . (However, in practice one would terminate for a much smaller value of k.)

Of course, in any implementation of (10.3.4), one stores the Householder vectors $v _ { k }$ and never explicitly forms the corresponding matrix product. Since we have $H _ { k } ( 1 { : } k , 1 { : } k ) = I _ { k }$ there is no need to compute the first k components of the vector w in (10.3.4) since we do not use them. (Ideally they are zero.)

Unfortunately, these economies make but a small dent in the computational overhead associated with complete reorthogonalization. The Householder calculations increase the work in the kth Lanczos step by $O ( k n )$ flops. Moreover, to compute $q _ { k + 1 }$ , the Householder vectors associated with $H _ { 0 } , \ldots , H _ { k }$ must be accessed. For large n and k, this usually implies a prohibitive level of memory traffic.

Thus, there is a high price associated with complete reorthogonalization. Fortunately, there are more effective courses of action to take, but these require a greater understanding of just how orthogonality is lost.

# 10.3.4 Selective Reorthogonalization

A remarkable, ironic consequence of the Paige (1971) error analysis is that loss of orthogonality goes hand in hand with convergence of a Ritz pair. To be precise, suppose the symmetric QR algorithm is applied to $\hat { T } _ { k }$ and renders computed Ritz values $\hat { \theta } _ { 1 } , \ldots , \hat { \theta } _ { k }$ and a nearly orthogonal matrix of eigenvectors $\hat { S } _ { k } = ( \hat { s } _ { p q } )$ . If

$$
\hat {Y} _ {k} = \left[ \begin{array}{c c c c} \hat {y} _ {1} & \dots & \hat {y} _ {k} \end{array} \right] = \mathsf {f l} (\hat {Q} _ {k} \hat {S} _ {k}),
$$

then it can be shown that for $i = 1 { : } k$ we have

$$
\left| \hat {q} _ {k + 1} ^ {T} \hat {y} _ {i} \right| \approx \frac {\mathbf {u} \| A \| _ {2}}{\left| \hat {\beta} _ {k} \right| \left| \hat {s} _ {k i} \right|} \tag {10.3.5}
$$

and

$$
\| A \hat {y} _ {i} - \hat {\theta} _ {i} \hat {y} _ {i} \| _ {2} \approx | \hat {\beta} _ {k} | | \hat {s} _ {k i} |. \tag {10.3.6}
$$

That is, the most recently computed Lanczos vector $\hat { q } _ { k + 1 }$ tends to have a nontrivial and unwanted component in the direction of any converged Ritz vector. Consequently, instead of orthogonalizing $\hat { q } _ { k + 1 }$ against all of the previously computed Lanczos vectors, we can achieve the same effect by orthogonalizing it against the much smaller set of converged Ritz vectors.

The practical aspects of enforcing orthogonality in this way are discussed in Parlett and Scott (1979). In their scheme, known as selective reorthogonalization, a computed Ritz pair $\{ \hat { \theta } , \hat { y } \}$ is called “good” if it satisfies

$$
\| A \hat {y} - \hat {\theta} \hat {y} \| _ {2} \leq \sqrt {\mathbf {u}} \| A \| _ {2}.
$$

As soon as $\hat { q } _ { k + 1 }$ is computed, it is orthogonalized against each good Ritz vector. This is much less costly than complete reorthogonalization, since, at least at first, there are many fewer good Ritz vectors than Lanczos vectors.

One way to implement selective reorthogonalization is to diagonalize $\hat { T } _ { k }$ at each step and then examine the $\hat { s } _ { k i }$ in light of (10.3.5) and (10.3.6). A more efficient approach for large k is to estimate the loss-of-orthogonality measure  $I _ { k } - \hat { Q } _ { k } ^ { T } \hat { Q } _ { k } \parallel _ { 2 }$ using the following result.

Lemma 10.3.1. Suppose $S _ { + } = [ S d ]$ where $S \in \mathbb { R } ^ { n \times k }$ and $d \in \mathbb { R } ^ { n }$ . If

$$
\left\| I _ {k} - S ^ {T} S \right\| _ {2} \leq \mu \quad | 1 - d ^ {T} d | \leq \delta ,
$$

then

$$
\| I _ {k + 1} - S _ {+} ^ {T} S _ {+} \| _ {2} \leq \mu_ {+}
$$

where

$$
\mu_ {+} = \frac {1}{2} \left(\mu + \delta + \sqrt {(\mu - \delta) ^ {2} + 4 \| S ^ {T} d \| _ {2} ^ {2}}\right).
$$

Proof. See Kahan and Parlett (1974) or Parlett and Scott (1979).

Thus, if we have a bound for $\parallel I _ { k } - \hat { Q } _ { k } ^ { T } \hat { Q } _ { k } \parallel _ { 2 }$ , then by applying the lemma with $S = \hat { Q } _ { k }$ and $d = \hat { q } _ { k + 1 }$ we can generate a bound for $\parallel I _ { k + 1 } - \hat { Q } _ { k + 1 } ^ { T } \hat { Q } _ { k + 1 } \parallel _ { 2 }$ . (In this case $\delta \approx \mathbf { u }$ and we assume that $\hat { q } _ { k + 1 }$ has been orthogonalized against the set of currently good Ritz vectors.) It is possible to estimate the norm of $\hat { Q } _ { k } ^ { T } \hat { q } _ { k + 1 }$ from a simple recurrence that spares one the need to access $\hat { q } _ { 1 } , \dots , \hat { q } _ { k }$ . The overhead is minimal, and when the bounds signal loss of orthogonality, it is time to contemplate the enlargement of the set of good Ritz vectors. Then and only then is $\hat { T } _ { k }$ diagonalized.

# 10.3.5 The Ghost Eigenvalue Problem

Considerable effort has been spent in trying to develop a workable Lanczos procedure that does not involve any kind of orthogonality enforcement. Research in this direction focuses on the problem of “ghost” eigenvalues. These are multiple eigenvalues of $\hat { T } _ { k }$ that correspond to simple eigenvalues of A. They arise because the iteration essentially restarts itself when orthogonality to a converged Ritz vector is lost. (By way of analogy, consider what would happen during orthogonal iteration (8.2.8) if we “forgot” to orthogonalize.)

The problem of identifying ghost eigenvalues and coping with their presence is discussed by Cullum and Willoughby (1979) and Parlett and Reid (1981). It is a particularly pressing problem in those applications where all of $A \mathrm { { } i \mathrm { { s } } }$ eigenvalues are desired, for then the above orthogonalization procedures are expensive to implement.

Difficulties with the Lanczos iteration can be expected even if A has a genuinely multiple eigenvalue. This follows because the $\hat { T } _ { k }$ are unreduced, and unreduced tridiagonal matrices cannot have multiple eigenvalues. The next practical Lanczos procedure that we discuss attempts to circumvent this difficulty.

# 10.3.6 Block Lanczos Algorithm

Just as the simple power method has a block analogue in simultaneous iteration, so does the Lanczos algorithm have a block version. Suppose $n = r p$ and consider the

decomposition

$$
Q ^ {T} A Q = \bar {T} = \left[ \begin{array}{c c c c c} M _ {1} & B _ {1} ^ {T} & & \dots & 0 \\ B _ {1} & M _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & B _ {r - 1} ^ {T} \\ 0 & \dots & & B _ {r - 1} & M _ {r} \end{array} \right] \tag {10.3.7}
$$

where

$$
Q = \left[ X _ {1} \mid \dots \mid X _ {r} \right], \quad X _ {i} \in \mathbb {R} ^ {n \times p},
$$

is orthogonal, each $M _ { i } \in \mathbb { R } ^ { p \times p }$ , and each $B _ { i } \in \mathbb { R } ^ { p \times p }$ is upper triangular. Comparison of blocks in $A Q = Q { \bar { T } }$ shows that

$$
A X _ {k} = X _ {k - 1} B _ {k - 1} ^ {T} + X _ {k} M _ {k} + X _ {k + 1} B _ {k}
$$

for $k = 1 { : } r$ assuming $X _ { 0 } B _ { 0 } ^ { T } \equiv 0$ and $X _ { r + 1 } B _ { r } \equiv 0$ . From the orthogonality of $Q$ we have

$$
M _ {k} = X _ {k} ^ {T} A X _ {k}
$$

for $k = 1 { : } r$ . Moreover, if we define

$$
R _ {k} = A X _ {k} - X _ {k} M _ {k} - X _ {k - 1} B _ {k - 1} ^ {T} \in \mathbb {R} ^ {n \times p},
$$

then

$$
X _ {k + 1} B _ {k} = R _ {k}
$$

is a QR factorization of $R _ { k }$ . These observations suggest that the block tridiagonal matrix $\bar { T }$ in (10.3.7) can be generated as follows:

$$
X _ {1} \in \mathbb {R} ^ {n \times p} \text {   given   with   } X _ {1} ^ {T} X _ {1} = I _ {p}
$$

$$
M _ {1} = X _ {1} ^ {T} A X _ {1}
$$

$$
\text { for } k = 1: r - 1 \tag {10.3.8}
$$

$$
R _ {k} = A X _ {k} - X _ {k} M _ {k} - X _ {k - 1} B _ {k - 1} ^ {T} \quad (X _ {0} B _ {0} ^ {T} \equiv 0)
$$

$$
X _ {k + 1} B _ {k} = R _ {k} \quad \text {(QR factorization of R_ {k})}
$$

$$
M _ {k + 1} = X _ {k + 1} ^ {T} A X _ {k + 1}
$$

end

At the beginning of the kth pass through the loop we have

$$
A \left[ X _ {1} \mid \dots \mid X _ {k} \right] = \left[ X _ {1} \mid \dots \mid X _ {k} \right] \bar {T} _ {k} + R _ {k} \left[ 0 \mid \dots \mid 0 \mid I _ {p} \right], \tag {10.3.9}
$$

where

$$
\bar {T} _ {k} = \left[ \begin{array}{c c c c c} M _ {1} & B _ {1} ^ {T} & & \dots & 0 \\ B _ {1} & M _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & B _ {k - 1} ^ {T} \\ 0 & \dots & & B _ {k - 1} & M _ {k} \end{array} \right].
$$

Using an argument similar to the one used in the proof of Theorem 10.1.1, we can show that the $X _ { k }$ are mutually orthogonal provided none of the $R _ { k }$ is rank-deficient. However if rank $( R _ { k } ) < p$ for some k, then it is possible to choose the columns of $X _ { k + 1 }$ such that $X _ { k + 1 } ^ { T } X _ { i } = 0$ , for $i = 1 { : } k$ . See Golub and Underwood (1977).

Because $\hat { T } _ { k }$ has bandwidth p, it can be efficiently reduced to tridiagonal form using an algorithm of Schwartz (1968). Once tridiagonal form is achieved, the Ritz values can be obtained via the symmetric QR algorithm or any of the special methods of §8.4. In order to decide intelligently when to use block Lanczos, it is necessary to understand how the block dimension affects convergence of the Ritz values. The following generalization of Theorem 10.1.2 sheds light on this issue.

Theorem 10.3.2. Let A be an n-by-n symmetric matrix with Schur decomposition

$$
Z ^ {T} A Z = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right), \quad \lambda_ {1} \geq \dots \geq \lambda_ {n}, \quad Z = \left[ z _ {1} \mid \dots \mid z _ {n} \right].
$$

Let $\mu _ { 1 } \geq \cdots \geq \mu _ { p }$ be the p largest eigenvalues of the matrix $\hat { T } _ { k }$ obtained after k steps of $( 1 0 . 3 . 8 )$ . Suppose $Z _ { 1 } = { \left[ \begin{array} { l } { z _ { 1 } } \end{array} | \cdot \cdot \cdot \ \right] } z _ { p } \ ]$ and

$$
0 <   \cos (\phi_ {p}) = \sigma_ {p} (Z _ {1} ^ {T} X _ {1}),
$$

the smallest singular value of $Z _ { 1 } ^ { T } X _ { 1 }$ . Then for $i = 1 { : } p _ { : }$ ,

$$
\lambda_ {i} \geq \mu_ {i} \geq \lambda_ {i} - (\lambda_ {1} - \lambda_ {n}) \left(\frac {\tan (\theta_ {p})}{c _ {k - 1} (1 + 2 \rho_ {i})}\right) ^ {2}
$$

where

$$
\rho_ {i} = \frac {\lambda_ {i} - \lambda_ {p + 1}}{\lambda_ {p + 1} - \lambda_ {n}}
$$

and $c _ { k - 1 } ( z )$ is the Chebyshev polynomial of degree $k - 1$ .

Proof. See Underwood (1975). Compare with Theorem 10.1.2.

Analogous inequalities can be obtained for $\bar { T } _ { k } ^ { \ , } \mathrm { s }$ smallest eigenvalues by applying the theorem with A replaced by −A. Based on the theorem and scrutiny of (10.3.8), we conclude that

• the error bounds for the Ritz values improve with increased $p$   
• the amount of work required to compute $\hat { T } _ { k }$ ’s eigenvalues is proportional to $k p ^ { 2 }$   
• the block dimension should be at least as large as the largest multiplicity of any sought-after eigenvalue.

Determination of the block dimension in the face of these trade-offs is discussed in detail by Scott (1979). We mention that loss of orthogonality also plagues the block Lanczos algorithm. However, all of the orthogonality enforcement schemes described above can be extended to the block setting.

# 10.3.7 Block Lanczos Algorithm with Restarting

The block Lanczos algorithm (10.3.8) can be used in an iterative fashion to calculate selected eigenvalues of A. To fix ideas, suppose we wish to calculate the p largest eigenvalues. If $X _ { 1 } \in \mathbb { R } ^ { n \times p }$ is a given matrix having orthonormal columns, then it can be refined as follows:

Step 1. Generate $\boldsymbol { X } _ { 2 } , \ldots , \boldsymbol { X } _ { s } \in \mathbb { R } ^ { n \times p }$ via the block Lanczos algorithm.

Step 2. Form ${ \bar { T } } _ { s } = [ X _ { 1 } | \cdots | X _ { s } ] ^ { T } A [ X _ { 1 } | \cdots | X _ { s } ]$ , an sp-by-sp matrix that has bandwidth p.

Step 3. Compute an orthogonal matrix $U = \ [ \boldsymbol { u } _ { 1 } \vert \cdot \cdot \cdot \vert \ u _ { s p } ]$ such that $U ^ { T } \bar { T } _ { s } U =$ di $\arg ( \theta _ { 1 } , \ldots , \theta _ { s p } )$ with $\theta _ { 1 } \ge \cdots \ge \theta _ { s p }$ .

Step 4. Set $X _ { 1 } ^ { ( \mathrm { n e w } ) } = [ X _ { 1 } | \cdots | X _ { s } ] [ u _ { 1 } | \cdots | u _ { p } ]$

This is the block analog of the s-step Lanczos algorithm, which has been extensively analyzed by Cullum and Donath (1974) and Underwood (1975). The same idea can be used to compute several of A’s smallest eigenvalues or a mixture of both large and small eigenvalues. See Cullum (1978). The choice of the parameters s and p depends upon storage constraints as well as upon the block-size implications that we discussed above. The value of p can be diminished as the good Ritz vectors emerge. However, this demands that orthogonality to the converged vectors be enforced.

# Problems

P10.3.1 Rearrange (10.3.4) and (10.3.8) so that they require one matrix-vector product per iteration.

P10.3.2 If rank $( R _ { k } ) < p$ in (10.3.8), does it follow that ran $( \left[ X _ { 1 } \mid \cdots \mid X _ { k } \right] )$ contains an eigenvector of A?

# Notes and References for §10.3

The behavior of the Lanczos method in the presence of roundoff error was originally reported in:

C.C. Paige (1971). “The Computation of Eigenvalues and Eigenvectors of Very Large Sparse Matrices,” PhD thesis, University of London.

Important follow-up papers include:

C.C. Paige (1972). “Computational Variants of the Lanczos Method for the Eigenproblem,” J. Inst. Math. Applic. 10, 373–381.

C.C. Paige (1976). “Error Analysis of the Lanczos Algorithm for Tridiagonalizing a Symmetric Matrix,” J. Inst. Math. Applic. 18, 341–349.

C.C. Paige (1980). “Accuracy and Effectiveness of the Lanczos Algorithm for the Symmetric Eigenproblem,” Lin. Alg. Applic. 34, 235–258.

For additional analysis of the method, see Parlett (SEP), Meurant (LCG) as well as:

D.S. Scott (1979). “How to Make the Lanczos Algorithm Converge Slowly,” Math. Comput. 33, 239–247.

B.N. Parlett, H.D. Simon, and L.M. Stringer (1982). “On Estimating the Largest Eigenvalue with the Lanczos Algorithm,” Math. Comput. 38, 153–166.

B.N. Parlett and B. Nour-Omid (1985). “The Use of a Refined Error Bound When Updating Eigenvalues of Tridiagonals,” Lin. Alg. Applic. 68, 179–220.

J. Kuczy´nski and H. Wo´zniakowski (1992). “Estimating the Largest Eigenvalue by the Power and Lanczos Algorithms with a Random Start,” SIAM J. Matrix Anal. Applic. 13, 1094–1122.   
G. Meurant and Z. Strakos (2006). “The Lanczos and Conjugate Gradient Algorithms in Finite Precision Arithmetic,” Acta Numerica 15, 471–542.   
A wealth of practical, Lanczos-related information may be found in:   
J.K. Cullum and R.A. Willoughby (2002). Lanczos Algorithms for Large Symmetric Eigenvalue Computations: Vol. I: Theory, SIAM Publications, Philadelphia, PA.   
J. Brown, M. Chu, D. Ellison, and R. Plemmons (1994). Proceedings of the Cornelius Lanczos International Centenary Conference, SIAM Publications, Philadelphia, PA.   
For a discussion about various reorthogonalization schemes, see:   
C.C. Paige (1970). “Practical Use of the Symmetric Lanczos Process with Reorthogonalization,” BIT 10, 183–195.   
G.H. Golub, R. Underwood, and J.H. Wilkinson (1972). “The Lanczos Algorithm for the Symmetric Ax = λBx Problem,” Report STAN-CS-72-270, Department of Computer Science, Stanford University, Stanford, CA.   
B.N. Parlett and D.S. Scott (1979). “The Lanczos Algorithm with Selective Orthogonalization,” Math. Comput. 33, 217–238.   
H.D. Simon (1984). “Analysis of the Symmetric Lanczos Algorithm with Reorthogonalization Methods,” Lin. Alg. Applic. 61, 101–132.   
Without any reorthogonalization it is necessary either to monitor the loss of orthogonality and quit at the appropriate instant or else to devise a scheme that will identify unconverged eigenvalues and false multiplicities, see:   
W. Kahan and B.N. Parlett (1976). “How Far Should You Go with the Lanczos Process?” in Sparse Matrix Computations, J.R. Bunch and D.J. Rose (eds.), Academic Press, New York, 131–144.   
J. Cullum and R.A. Willoughby (1979). “Lanczos and the Computation in Specified Intervals of the Spectrum of Large, Sparse Real Symmetric Matrices, in Sparse Matrix Proc., I.S. Duff and G.W. Stewart (eds.), SIAM Publications, Philadelphia, PA.   
B.N. Parlett and J.K. Reid (1981). “Tracking the Progress of the Lanczos Algorithm for Large Symmetric Eigenproblems,” IMA J. Num. Anal. 1, 135–155.   
For a restarting framework to be successful, it must exploit the approximate invariant subspace information that has been acquired by the iteration that is about to be shut down, see:   
D. Calvetti, L. Reichel, and D.C. Sorensen (1994). “An Implicitly Restarted Lanczos Method for Large Symmetric Eigenvalue Problems,” ETNA 2, 1–21.   
K. Wu and H. Simon (2000). “Thick-Restart Lanczos Method for Large Symmetric Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 22, 602–616.   
The block Lanczos algorithm is discussed in:   
J. Cullum and W.E. Donath (1974). “A Block Lanczos Algorithm for Computing the q Algebraically Largest Eigenvalues and a Corresponding Eigenspace of Large Sparse Real Symmetric Matrices,” Proceedings of the 1974 IEEE Conference on Decision and Control, Phoenix, AZ, 505–509.   
R. Underwood (1975). “An Iterative Block Lanczos Method for the Solution of Large Sparse Symmetric Eigenvalue Problems,” Report STAN-CS-75-495, Department of Computer Science, Stanford University, Stanford, CA.   
G.H. Golub and R. Underwood (1977). “The Block Lanczos Method for Computing Eigenvalues,” in Mathematical Software III , J. Rice (ed.), Academic Press, New York, pp. 364–377.   
J. Cullum (1978). “The Simultaneous Computation of a Few of the Algebraically Largest and Smallest Eigenvalues of a Large Sparse Symmetric Matrix,” BIT 18, 265–275.   
A. Ruhe (1979). “Implementation Aspects of Band Lanczos Algorithms for Computation of Eigenvalues of Large Sparse Symmetric Matrices,” Math. Comput. 33, 680–687.   
The block Lanczos algorithm generates a symmetric band matrix whose eigenvalues can be computed in any of several ways. One approach is described in:   
H.R. Schwartz (1968). “Tridiagonalization of a Symmetric Band Matrix,” Numer. Math. 12, 231–241.

In some applications it is necessary to obtain estimates of interior eigenvalues. One strategy is to apply Lanczos to the matrix $( A - \mu I ) ^ { - 1 }$ because the extremal eigenvalues of this matrix are eigenvalues close to µ. However, “shift-and-invert” strategies replace the matrix-vector product in the Lanczos iteration with a large sparse linear equation solve, see:

A.K. Cline, G.H. Golub, and G.W. Platzman (1976). “Calculation of Normal Modes of Oceans Using a Lanczos Method,” in Sparse Matrix Computations, J.R. Bunch and D.J. Rose (eds), Academic Press, New York, pp. 409–426.   
T. Ericsson and A. Ruhe (1980). “The Spectral Transformation Lanczos Method for the Numerical Solution of Large Sparse Generalized Symmetric Eigenvalue Problems,” Math. Comput. 35, 1251– 1268.   
R.B. Morgan (1991). “Computing Interior Eigenvalues of Large Matrices,” Lin. Alg. Applic. 154-156, 289–309.   
R.G. Grimes, J.G. Lewis, and H.D. Simon (1994). “A Shifted Block Lanczos Algorithm for Solving Sparse Symmetric Generalized Eigenproblems,” SIAM J. Matrix Anal. Applic. 15, 228–272.
