# 7.7 The Generalized Eigenvalue Problem

If A, $B \in \mathbb { C } ^ { n \times n }$ , then the set of all matrices of the form $A - \lambda B$ with $\lambda \in \mathbb { C }$ is a pencil. The generalized eigenvalues of $A - \lambda B$ are elements of the set $\lambda ( A , B )$ defined by

$$
\lambda (A, B) = \{z \in \mathbb {C}: \det (A - z B) = 0 \}.
$$

If $\lambda \in \lambda ( A , B )$ and $0 \neq x \in \mathbb { C } ^ { n }$ satisfies

$$
A x = \lambda B x, \tag {7.7.1}
$$

then x is an eigenvector of A − λB. The problem of finding nontrivial solutions to (7.7.1) is the generalized eigenvalue problem and in this section we survey some of its mathematical properties and derive a stable method for its solution. We briefly discuss how a polynomial eigenvalue problem can be converted into an equivalent generalized eigenvalue problem through a linearization process.
