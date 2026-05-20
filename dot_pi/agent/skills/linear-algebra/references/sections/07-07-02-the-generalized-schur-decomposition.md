# 7.7.2 The Generalized Schur Decomposition

From the numerical point of view, it makes to insist that the transformation matrices Q and Z be unitary. This leads to the following decomposition described in Moler and Stewart (1973).

Theorem 7.7.1 (Generalized Schur Decomposition). If A and B are in $\mathbb { C } ^ { n \times n }$ , then there exist unitary Q and Z such that $Q ^ { H } A Z = T$ and $Q ^ { H } B Z = S$ are upper triangular. If for some $k , t _ { k k }$ and $s _ { k k }$ are both zero, then $\lambda ( A , B ) = \mathbb { C }$ . Otherwise

$$
\lambda (A, B) = \{t _ {i i} / s _ {i i}: s _ {i i} \neq 0 \}.
$$

Proof. Let $\{ B _ { k } \}$ be a sequence of nonsingular matrices that converge to B. For each k, let

$$
Q _ {k} ^ {H} (A B _ {k} ^ {- 1}) Q _ {k} = R _ {k}
$$

be a Schur decomposition of $A B _ { k } ^ { - 1 }$ . Let $Z _ { k }$ be unitary such that

$$
Z _ {k} ^ {H} (B _ {k} ^ {- 1} Q _ {k}) = S _ {k} ^ {- 1}
$$

is upper triangular. It follows that $Q _ { k } ^ { H } A Z _ { k } \ : = \ : R _ { k } S _ { k }$ and $Q _ { k } ^ { H } B _ { k } Z _ { k } \ = \ S _ { k }$ are also upper triangular. Using the Bolzano-Weierstrass theorem, we know that the bounded sequence $\{ ( Q _ { k } , Z _ { k } ) \}$ has a converging subsequence,

$$
\lim _ {i \to \infty} (Q _ {k _ {i}}, Z _ {k _ {i}}) = (Q, Z).
$$

It is easy to show that Q and Z are unitary and that $Q ^ { H } A Z$ and $Q ^ { H } B Z$ are upper triangular. The assertions about $\lambda ( A , B )$ follow from the identity

$$
\det (A - \lambda B) = \det (Q Z ^ {H}) \prod_ {i = 1} ^ {n} \left(t _ {i i} - \lambda s _ {i i}\right)
$$

and that completes the proof of the theorem.

If A and B are real then the following decomposition, which corresponds to the real Schur decomposition (Theorem 7.4.1), is of interest.

Theorem 7.7.2 (Generalized Real Schur Decomposition). If A and B are in $\mathbb { R } ^ { n \times n }$ then there exist orthogonal matrices Q and Z such that $Q ^ { T } A Z$ is upper quasitriangular and $Q ^ { T } B Z$ is upper triangular.

Proof. See Stewart (1972).

In the remainder of this section we are concerned with the computation of this decomposition and the mathematical insight that it provides.
