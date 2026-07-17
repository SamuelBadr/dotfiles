# 7.7.8 Generalized Invariant Subspace Computations

Many of the invariant subspace computations discussed in §7.6 carry over to the generalized eigenvalue problem. For example, approximate eigenvectors can be found via inverse iteration:

$$
q ^ {(0)} \in \mathbb {C} ^ {n \times n} \text {   given.   }
$$

for k = 1, 2, . . .

$$
\text { Solve } (A - \mu B) z ^ {(k)} = B q ^ {(k - 1)}.
$$

$$
\text { Normalize: } q ^ {(k)} = z ^ {(k)} / \| z ^ {(k)} \| _ {2}.
$$

$$
\lambda^ {(k)} = [ q ^ {(k)} ] ^ {H} A q ^ {(k)} / [ q ^ {(k)} ] ^ {H} A q ^ {(k)}
$$

end

If B is nonsingular, then this is equivalent to applying (7.6.1) with the matrix $B ^ { - 1 } A .$ . Typically, only a single iteration is required if µ is an approximate eigenvalue computed by the QZ algorithm. By inverse iterating with the Hessenberg-triangular pencil, costly accumulation of the Z-transformations during the QZ iteration can be avoided.

Corresponding to the notion of an invariant subspace for a single matrix, we have the notion of a deflating subspace for the pencil $A - \lambda B$ . In particular, we say that a k-dimensional subspace $S \subseteq \mathbb { C } ^ { n }$ is deflating for the pencil $A - \lambda B$ if the subspace $\{ A x + B y : x , y \in S \}$ has dimension k or less. Note that if

$$
Q ^ {H} A Z = T, \quad Q ^ {H} B Z = S
$$

is a generalized Schur decomposition of A−λB, then the columns of Z in the generalized Schur decomposition define a family of deflating subspaces. Indeed, if

$$
Q = \left[ q _ {1} \mid \dots \mid q _ {n} \right], \qquad Z = \left[ z _ {1} \mid \dots \mid z _ {n} \right]
$$

are column partitionings, then

$$
\operatorname{span} \left\{A z _ {1}, \dots , A z _ {k} \right\} \subseteq \operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\},
$$

$$
\operatorname{span} \left\{B z _ {1}, \dots , B z _ {k} \right\} \subseteq \operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\},
$$

for $k = 1 { : } n$ . Properties of deflating subspaces and their behavior under perturbation are described in Stewart (1972).
