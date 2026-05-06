# 12.5.5 Tensor Rank

The choice of r in the CP approximation problem brings us to the complicated issue of tensor rank. If

$$
\mathcal {A} = \sum_ {j = 1} ^ {r} \lambda_ {j} f _ {j} ^ {(1)} \circ \dots \circ f _ {j} ^ {(d)}
$$

and no shorter sum-of-rank-1’s exists, then we say that A is a rank-r tensor. Thus, we see that in the CP approximation problem is a problem of finding the best rank-r approximation. Using the CP framework to discover the rank of a tensor is problematic because of the following complications.

Complication 1. The tensor rank problem is NP-hard. See and Hillar and Lim (2012).

Complication 2. The largest rank attainable for an $n _ { 1 } \substack { \mathrm { - b y - } \cdot \cdot \cdot n _ { d } }$ tensor is called the maximum rank. There is no simple formula like min $\{ n _ { 1 } , \ldots , n _ { d } \}$ . Indeed, maximum rank is known for only a handful of special cases.

Complication 3. If the set of rank-k tensors in $\mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ has positive measure, then k is a typical rank. The space of $n _ { 1 } \times \cdots \times n _ { d }$ can have more than one typical rank. For example, the probability that a random 2-by-2-by-2 tensor has rank 2 is .79, while the probability that it has rank 3 is .21, assuming that the $a _ { i j k }$ are normally distributed with mean 0 and variance 1. See de Silva and Lim (2008) and Martin (2011) for detailed analysis of the 2-by-2-by 2 case.

Complication 4. The rank of a particular tensor over the real field may be different than its rank over the complex field.

Complication 5. There exist tensors that can be approximated with arbitrary precision by a tensor of lower rank. Such a tensor is said to be degenerate.

Complication 6. If

$$
\mathcal {X} _ {r} = \sum_ {j = 1} ^ {r + 1} \lambda_ {j} U _ {1} (:, j) \circ \dots \circ U _ {d} (:, j)
$$

is the best rank-(r + 1) approximation of A, then it does not follow that

$$
\mathcal {X} _ {r + 1} = \sum_ {j = 1} ^ {r} \lambda_ {j} \hat {U} _ {1} (:, j) \circ \dots \circ \hat {U} _ {d} (:, j)
$$

is the best rank-r approximation of A. See Kolda (2003) for an example. Subtracting the best rank-1 approximation can even increase the rank! See Stegeman and Comon (2009).

See Kolda and Bader (2009) for references on tensor rank and its implications for computation. Examples that illuminate the subtleties associated with tensor rank can be found in the the paper by de Silva and Lim (2008).
