# Chapter 7

# Unsymmetric Eigenvalue Problems

7.1 Properties and Decompositions   
7.2 Perturbation Theory   
7.3 Power Iterations   
7.4 The Hessenberg and Real Schur Forms   
7.5 The Practical QR Algorithm   
7.6 Invariant Subspace Computations   
7.7 The Generalized Eigenvalue Problem   
7.8 Hamiltonian and Product Eigenvalue Problems   
7.9 Pseudospectra

Having discussed linear equations and least squares, we now direct our attention to the third major problem area in matrix computations, the algebraic eigenvalue problem. The unsymmetric problem is considered in this chapter and the more agreeable symmetric case in the next.

Our first task is to present the decompositions of Schur and Jordan along with the basic properties of eigenvalues and invariant subspaces. The contrasting behavior of these two decompositions sets the stage for §7.2 in which we investigate how the eigenvalues and invariant subspaces of a matrix are affected by perturbation. Condition numbers are developed that permit estimation of the errors induced by roundoff.

The key algorithm of the chapter is the justly famous QR algorithm. This procedure is one of the most complex algorithms presented in the book and its development is spread over three sections. We derive the basic QR iteration in §7.3 as a natural generalization of the simple power method. The next two sections are devoted to making this basic iteration computationally feasible. This involves the introduction of the Hessenberg decomposition in §7.4 and the notion of origin shifts in §7.5.

The QR algorithm computes the real Schur form of a matrix, a canonical form that displays eigenvalues but not eigenvectors. Consequently, additional computations usually must be performed if information regarding invariant subspaces is desired. In §7.6, which could be subtitled, “What to Do after the Real Schur Form is Calculated,” we discuss various invariant subspace calculations that can be performed after the QR algorithm has done its job.

The next two sections are about Schur decomposition challenges. The generalized eigenvalue problem $A x = \lambda B x$ is the subject of §7.7. The challenge is to compute the Schur decomposition of $B ^ { - 1 } A$ without actually forming the indicated inverse or the product. The product eigenvalue problem is similar, only arbitrarily long sequences of products are considered. This is treated in §7.8 along with the Hamiltonian eigenproblem where the challenge is to compute a Schur form that has a special 2-by-2 block structure.

In the last section the important notion of pseudospectra is introduced. It is sometimes the case in unsymmetric matrix problems that traditional eigenvalue analysis fails to tell the “whole story” because the eigenvector basis is ill-conditioned. The pseudospectra framework effectively deals with this issue.

We mention that it is handy to work with complex matrices and vectors in the more theoretical passages that follow. Complex versions of the QR factorization, the singular value decomposition, and the CS decomposition surface in the discussion.

# Reading Notes

Knowledge of Chapters 1–3 and §§5.1–§5.2 are assumed. Within this chapter there are the following dependencies:

$$
\begin{array}c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
$$

Excellent texts for the dense eigenproblem include Chatelin (EOM), Kressner (NMSE), Stewart (MAE), Stewart and Sun (MPA), Watkins (MEP), and Wilkinson (AEP).

# 7.1 Properties and Decompositions

In this section the background necessary to develop and analyze the eigenvalue algorithms that follow are surveyed. For further details, see Horn and Johnson (MA).

# 7.1.1 Eigenvalues and Invariant Subspaces

The eigenvalues of a matrix $A \in \mathbb { C } ^ { n \times n }$ are the n roots of its characteristic polynomial $p ( z ) = \mathsf { d e t } ( z I - A )$ . The set of these roots is called the spectrum of A and is denoted by

$$
\lambda (A) = \{z: \det (z I - A) = 0 \}.
$$

If $\lambda ( A ) = \left\{ \lambda _ { 1 } , \ldots , \lambda _ { n } \right\}$ , then

$$
\det (A) = \lambda_ {1} \lambda_ {2} \dots \lambda_ {n}
$$

and

$$
\operatorname{tr} (A) = \lambda_ {1} + \dots + \lambda_ {n}
$$

where the trace function, introduced in §6.4.1, is the sum of the diagonal entries, i.e.,

$$
\operatorname{tr} (A) = \sum_ {i = 1} ^ {n} a _ {i i}.
$$

These characterizations of the determinant and the trace follow by looking at the constant term and the coefficient of $z ^ { n - 1 }$ in the characteristic polynomial.

Four other attributes associated with the spectrum of $A \in \mathbb { C } ^ { n \times n }$ include the

$$
\text { Spectral   Radius }: \quad \rho (A) = \max _ {\lambda \in \lambda (A)} | \lambda |, \tag {7.1.1}
$$

$$
\text { Spectral   Abscissa }: \quad \alpha (A) = \max _ {\lambda \in \lambda (A)} \operatorname{Re} (\lambda), \tag {7.1.2}
$$

$$
\text { Numerical   Radius }: \quad r (A) = \max _ {\lambda \in \lambda (A)} \left\{\left| x ^ {H} A x \right|: \| x \| _ {2} = 1 \right\}, \tag {7.1.3}
$$

$$
\text { Numerical   Range: } \quad W (A) = \{x ^ {H} A x: \| x \| _ {2} = 1 \}. \tag {7.1.4}
$$

The numerical range, which is sometimes referred to as the field of values, obviously includes λ(A). It can be shown that $W ( A )$ is convex.

If $\lambda \in \lambda ( A )$ , then the nonzero vectors $x \in \mathbb { C } ^ { n }$ that satisfy $A x = \lambda x$ are eigenvectors. More precisely, x is a right eigenvector for λ if $A x = \lambda x$ and a left eigenvector if $x ^ { H } A = \lambda x ^ { H }$ . Unless otherwise stated, “eigenvector” means “right eigenvector.”

An eigenvector defines a 1-dimensional subspace that is invariant with respect to premultiplication by A. A subspace $S \subseteq \mathbb { C } ^ { n }$ with the property that

$$
x \in S \Longrightarrow A x \in S
$$

is said to be invariant (for A). Note that if

$$
A X = X B, \qquad B \in \mathbb {C} ^ {k \times k}, X \in \mathbb {C} ^ {n \times k},
$$

then $\mathsf { r a n } ( X )$ is invariant and $B y = \lambda y \Rightarrow A ( X y ) = \lambda ( X y )$ . Thus, if X has full column rank, then $A X = X B$ implies that $\lambda ( B ) \subseteq \lambda ( A )$ . If X is square and nonsingular, then A and $B = X ^ { - 1 } A X$ are similar, X is a similarity transformation, and $\lambda ( A ) = \lambda ( B )$ .

# 7.1.2 Decoupling

Many eigenvalue computations involve breaking the given problem down into a collection of smaller eigenproblems. The following result is the basis for these reductions.

Lemma 7.1.1. If $T \in \mathbb { C } ^ { n \times n }$ is partitioned as follows,

$$
T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {q} ^ {p}
$$

then $\lambda ( T ) = \lambda ( T _ { 1 1 } ) \cup \lambda ( T _ { 2 2 } )$ .

Proof. Suppose

$$
T x = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right] = \lambda \left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right]
$$

where $x _ { 1 } \in \mathbb { C } ^ { p }$ and $x _ { 2 } \in \mathbb { C } ^ { q }$ . If $x _ { 2 } \neq 0$ , then $T _ { 2 2 } x _ { 2 } = \lambda x _ { 2 }$ and so $\lambda \in \lambda ( T _ { 2 2 } )$ . If $x _ { 2 } = 0$ , then $T _ { 1 1 } x _ { 1 } = \lambda x _ { 1 }$ and so $\lambda \in \lambda ( T _ { 1 1 } )$ . It follows that $\lambda ( T ) \subset \lambda ( T _ { 1 1 } ) \cup \lambda ( T _ { 2 2 } )$ . But since both λ(T ) and $\lambda ( T _ { 1 1 } ) \cup \lambda ( T _ { 2 2 } )$ have the same cardinality, the two sets are equal.

# 7.1.3 Basic Unitary Decompositions

By using similarity transformations, it is possible to reduce a given matrix to any one of several canonical forms. The canonical forms differ in how they display the eigenvalues and in the kind of invariant subspace information that they provide. Because of their numerical stability we begin by discussing the reductions that can be achieved with unitary similarity.

Lemma 7.1.2. If $A \in \mathbb { C } ^ { n \times n } , B \in \mathbb { C } ^ { p \times p }$ , and $\ b X \in \mathbb { C } ^ { n \times p }$ satisfy

$$
A X = X B, \quad \operatorname{rank} (X) = p, \tag {7.1.5}
$$

then there exists a unitary $Q \in \mathbb { C } ^ { n \times n }$ such that

$$
Q ^ {H} A Q = T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {n - p} ^ {p} \tag {7.1.6}
$$

and $\lambda ( T _ { 1 1 } ) = \lambda ( A ) \cap \lambda ( B )$ .

Proof. Let

$$
X = Q \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right], \qquad Q \in \mathbb {C} ^ {n \times n}, R _ {1} \in \mathbb {C} ^ {p \times p}
$$

be a QR factorization of X. By substituting this into (7.1.5) and rearranging we have

$$
\left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ T _ {2 1} & T _ {2 2} \end{array} \right] \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right] = \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right] B
$$

where

$$
Q ^ {H} A Q = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ T _ {2 1} & T _ {2 2} \end{array} \right] _ {n - p} ^ {p}.
$$

By using the nonsingularity of $R _ { 1 }$ and the equations $T _ { 2 1 } R _ { 1 } = 0$ and $T _ { 1 1 } R _ { 1 } = R _ { 1 } B$ , we can conclude that $T _ { 2 1 } = 0$ and $\lambda ( T _ { 1 1 } ) = \lambda ( B )$ . The lemma follows because from Lemma 7.1.1 we have $\lambda ( A ) = \lambda ( T ) = \lambda ( T _ { 1 1 } ) \cup \lambda ( T _ { 2 2 } )$ .

Lemma 7.1.2 says that a matrix can be reduced to block triangular form using unitary similarity transformations if we know one of its invariant subspaces. By induction we can readily establish the decomposition of Schur (1909).

Theorem 7.1.3 (Schur Decomposition). If $A \in \mathbb { C } ^ { n \times n }$ , then there exists a unitary $Q \in \mathbb { C } ^ { n \times n }$ such that

$$
Q ^ {H} A Q = T = D + N \tag {7.1.7}
$$

where $D = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ and $N \in \mathbb { C } ^ { n \times n }$ is strictly upper triangular. Furthermore, $Q$ can be chosen so that the eigenvalues $\lambda _ { i }$ appear in any order along the diagonal.

Proof. The theorem obviously holds if $n = 1$ . Suppose it holds for all matrices of order $n - 1$ or less. If $A x = \lambda x$ and $x \neq 0$ , then by Lemma 7.1.2 (with $B = ( \lambda ) )$ there exists a unitary $U$ such that

$$
U ^ {H} A U = \left[ \begin{array}{c c} \lambda & w ^ {H} \\ 0 & C \end{array} \right] _ {n - 1} ^ {1}.
$$

By induction there is a unitary $\tilde { U }$ such that $\tilde { U } ^ { H } C \tilde { U }$ is upper triangular. Thus, if $Q = U { \cdot } \mathrm { d i a g } ( 1 , \tilde { U } )$ , then $Q ^ { H } A Q$ is upper triangular.

If $Q = [  q _ { 1 } | \cdots |  q _ { n }  ]$ is a column partitioning of the unitary matrix $Q$ in (7.1.7), then the $q _ { i }$ are referred to as Schur vectors. By equating columns in the equations $A Q = Q T$ , we see that the Schur vectors satisfy

$$
A q _ {k} = \lambda_ {k} q _ {k} + \sum_ {i = 1} ^ {k - 1} n _ {i k} q _ {i}, \quad k = 1: n. \tag {7.1.8}
$$

From this we conclude that the subspaces

$$
S _ {k} = \operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\}, \quad k = 1: n,
$$

are invariant. Moreover, it is not hard to show that if $Q _ { k } \ = \ \left[ \ q _ { 1 } \ | \cdot \cdot \cdot | \ q _ { k } \ \right]$ , then $\lambda ( Q _ { k } ^ { H } A Q _ { k } ) = \{ \lambda _ { 1 } , . . . , \lambda _ { k } \}$ . Since the eigenvalues in (7.1.7) can be arbitrarily ordered, it follows that there is at least one k-dimensional invariant subspace associated with each subset of k eigenvalues. Another conclusion to be drawn from (7.1.8) is that the Schur vector $q _ { k }$ is an eigenvector if and only if the kth column of N is zero. This turns out to be the case for k = 1:n whenever $A ^ { H } A = A A ^ { H }$ . Matrices that satisfy this property are called normal.

Corollary 7.1.4. $A \in \mathbb { C } ^ { n \times n }$ is normal if and only if there exists a unitary $Q \in \mathbb { C } ^ { n \times n }$ such that $Q ^ { H } A Q = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ .

Proof. See P7.1.1.

Note that if $Q ^ { H } A Q = T = \mathrm { d i a g } ( \lambda _ { i } ) + N$ is a Schur decomposition of a general n-by-n matrix A, then $\Vert N \Vert _ { F }$ is independent of the choice of $Q \mathrm { : }$

$$
\| N \| _ {F} ^ {2} = \| A \| _ {F} ^ {2} - \sum_ {i = 1} ^ {n} | \lambda_ {i} | ^ {2} \equiv \Delta^ {2} (A).
$$

This quantity is referred to as A’s departure from normality. Thus, to make T “more diagonal,” it is necessary to rely on nonunitary similarity transformations.

# 7.1.4 Nonunitary Reductions

To see what is involved in nonunitary similarity reduction, we consider the block diagonalization of a 2-by-2 block triangular matrix.

Lemma 7.1.5. Let $T \in \mathbb { C } ^ { n \times n }$ be partitioned as follows:

$$
T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {q} ^ {p}.
$$

Define the linear transformation $\phi { : } \mathbb { C } ^ { p \times q } \to \mathbb { C } ^ { p \times q }$ by

$$
\phi (X) = T _ {1 1} X - X T _ {2 2}
$$

where $X \in \mathbb { C } ^ { p \times q }$ . Then $\phi$ is nonsingular if and only if $\lambda ( T _ { 1 1 } ) \cap \lambda ( T _ { 2 2 } ) = \emptyset$ . If φ is nonsingular and Y is defined by

$$
Y = \left[ \begin{array}{c c} I _ {p} & Z \\ 0 & I _ {q} \end{array} \right]
$$

where $\phi ( Z ) = - T _ { 1 2 }$ , then $Y ^ { - 1 } T Y = \mathrm { d i a g } ( T _ { 1 1 } , T _ { 2 2 } )$

Proof. Suppose φ(X) = 0 for $X \neq 0$ and that

$$
U ^ {H} X V = \left[ \begin{array}{c c} \Sigma_ {r} & 0 \\ 0 & 0 \end{array} \right] _ {p - r} ^ {r}
$$

is the SVD of X with $\Sigma _ { r } = \mathrm { d i a g } ( \sigma _ { i } ) , r = \mathsf { r a n k } ( X )$ . Substituting this into the equation $T _ { 1 1 } X = X T _ { 2 2 }$ gives

$$
\left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right] \left[ \begin{array}{c c} \Sigma_ {r} & 0 \\ 0 & 0 \end{array} \right] = \left[ \begin{array}{c c} \Sigma_ {r} & 0 \\ 0 & 0 \end{array} \right] \left[ \begin{array}{c c} B _ {1 1} & B _ {1 2} \\ B _ {2 1} & B _ {2 2} \end{array} \right]
$$

where ${ U ^ { H } T _ { 1 1 } U = ( A _ { i j } ) }$ and $V ^ { H } T _ { 2 2 } V = \left( B _ { i j } \right)$ . By comparing blocks in this equation it is clear that $A _ { 2 1 } = 0 , B _ { 1 2 } = 0$ , and $\lambda ( A _ { 1 1 } ) = \lambda ( B _ { 1 1 } )$ . Consequently, $A _ { 1 1 }$ and $B _ { 1 1 }$ have an eigenvalue in common and that eigenvalue is in $\lambda ( T _ { 1 1 } ) \cap \lambda ( T _ { 2 2 } )$ . Thus, if φ is singular, then $T _ { 1 1 }$ and $T _ { 2 2 }$ have an eigenvalue in common. On the other hand, if $\lambda \in \lambda ( T _ { 1 1 } ) \cap \lambda ( T _ { 2 2 } )$ , then we have eigenvector equations $T _ { 1 1 } x = \lambda x$ and $y ^ { H } T _ { 2 2 } = \lambda y ^ { H }$ . A calculation shows that $\phi ( x y ^ { H } ) = 0$ confirming that φ is singular.

Finally, if $\phi$ is nonsingular, then $\phi ( Z ) = - T _ { 1 2 }$ has a solution and

$$
Y ^ {- 1} T Y = \left[ \begin{array}{c c} I _ {p} & - Z \\ 0 & I _ {q} \end{array} \right] \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] \left[ \begin{array}{c c} I _ {p} & Z \\ 0 & I _ {q} \end{array} \right] = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 1} Z - Z T _ {2 2} + T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right]
$$

has the required block diagonal form.

By repeatedly applying this lemma, we can establish the following more general result.

Theorem 7.1.6 (Block Diagonal Decomposition). Suppose

$$
Q ^ {H} A Q = T = \left[ \begin{array}{c c c c} T _ {1 1} & T _ {1 2} & \dots & T _ {1 q} \\ 0 & T _ {2 2} & \dots & T _ {2 q} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & T _ {q q} \end{array} \right] \tag {7.1.9}
$$

is a Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ and that the $T _ { i i }$ are square. $I f \lambda ( T _ { i i } ) \cap \lambda ( T _ { j j } ) = \emptyset$ whenever $i \neq j$ , then there exists a nonsingular matrix $Y \in \mathbb { C } ^ { n \times n }$ such that

$$
(Q Y) ^ {- 1} A (Q Y) = \operatorname{diag} \left(T _ {1 1}, \dots , T _ {q q}\right). \tag {7.1.10}
$$

Proof. See P7.1.2.

If each diagonal block $T _ { i i }$ is associated with a distinct eigenvalue, then we obtain

Corollary 7.1.7. If $A \in \mathbb { C } ^ { n \times n }$ , then there exists a nonsingular X such that

$$
X ^ {- 1} A X = \operatorname{diag} \left(\lambda_ {1} I + N _ {1}, \dots , \lambda_ {q} I + N _ {q}\right) \quad N _ {i} \in \mathbb {C} ^ {n _ {i} \times n _ {i}} \tag {7.1.11}
$$

where $\lambda _ { 1 } , \ldots , \lambda _ { q }$ are distinct, the integers $n _ { 1 } , \ldots , n _ { q }$ satisfy $n _ { 1 } + \cdot \cdot \cdot + n _ { q } = n$ , and each $N _ { i }$ is strictly upper triangular.

A number of important terms are connected with decomposition (7.1.11). The integer $n _ { i }$ is referred to as the algebraic multiplicity of $\lambda _ { i }$ . If $n _ { i } = 1$ , then $\lambda _ { i }$ is said to be simple. The geometric multiplicity of $\lambda _ { i }$ equals the dimensions of null $( N _ { i } )$ , i.e., the number of linearly independent eigenvectors associated with $\lambda _ { i }$ . If the algebraic multiplicity of $\lambda _ { i }$ exceeds its geometric multiplicity, then $\lambda _ { i }$ is said to be a defective eigenvalue. A matrix with a defective eigenvalue is referred to as a defective matrix. Nondefective matrices are also said to be diagonalizable.

Corollary 7.1.8 (Diagonal Form). $A \in \mathbb { C } ^ { n \times n }$ is nondefective if and only if there exists a nonsingular $\ b X \in \mathbb { C } ^ { n \times n }$ such that

$$
X ^ {- 1} A X = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right). \tag {7.1.12}
$$

Proof. A is nondefective if and only if there exist independent vectors $\boldsymbol { x } _ { 1 } \ldots \boldsymbol { x } _ { n } \in \mathbb { C } ^ { n }$ and scalars $\lambda _ { 1 } , \ldots , \lambda _ { n }$ such that $A x _ { i } = \lambda _ { i } x _ { i }$ for $i = 1 { : } n$ . This is equivalent to the existence of a nonsingular $\ b X = [ \ b x _ { 1 } \ | \cdot \ b \cdot \cdot | \ \ b x _ { n } ] \in \mathbb { C } ^ { n \times n }$ such that $A X \ = \ X D$ where $D = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ .□

Note that if $y _ { i } ^ { H }$ is the ith row of $X ^ { - 1 }$ , then $y _ { i } ^ { H } A = \lambda _ { i } y _ { i } ^ { H }$ . Thus, the columns of $X ^ { - H }$ are left eigenvectors and the columns of X are right eigenvectors.

If we partition the matrix X in (7.1.11),

$$
X = \left[ \begin{array}{c} X _ {1} \mid \dots \mid X _ {q} \\ n _ {1} \end{array} \right]
$$

then $\mathbb { C } ^ { n } = \mathsf { r a n } ( X _ { 1 } ) \oplus \ldots \oplus \mathsf { r a n } ( X _ { q } )$ , a direct sum of invariant subspaces. If the bases for these subspaces are chosen in a special way, then it is possible to introduce even more zeroes into the upper triangular portion of $X ^ { - 1 } A X$ .

Theorem 7.1.9 (Jordan Decomposition). If $A \in \mathbb { C } ^ { n \times n }$ , then there exists a nonsingular $\ b X \in \mathbb { C } ^ { n \times n }$ such that $X ^ { - 1 } A X = \operatorname { d i a g } ( J _ { 1 } , \dots , J _ { q } )$ where

$$
J _ {i} = \left[ \begin{array}{c c c c c} \lambda_ {i} & 1 & & \dots & 0 \\ 0 & \lambda_ {i} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & 1 \\ 0 & \dots & & 0 & \lambda_ {i} \end{array} \right] \in \mathbb {C} ^ {n _ {i} \times n _ {i}}
$$

and $n _ { 1 } + \cdots + n _ { q } = n$ .

Proof. See Horn and Johnson (MA, p. 330)

The $J _ { i }$ are referred to as Jordan blocks. The number and dimensions of the Jordan blocks associated with each distinct eigenvalue are unique, although their ordering along the diagonal is not.

# 7.1.5 Some Comments on Nonunitary Similarity

The Jordan block structure of a defective matrix is difficult to determine numerically. The set of n-by-n diagonalizable matrices is dense in Cn×n, and thus, small changes in a defective matrix can radically alter its Jordan form. We have more to say about this in §7.6.5.

A related difficulty that arises in the eigenvalue problem is that a nearly defective matrix can have a poorly conditioned matrix of eigenvectors. For example, any matrix X that diagonalizes

$$
A = \left[ \begin{array}{c c} 1 + \epsilon & 1 \\ 0 & 1 - \epsilon \end{array} \right], \qquad 0 <   \epsilon \ll 1, \tag {7.1.13}
$$

has a 2-norm condition of order $1 / \epsilon .$

These observations serve to highlight the difficulties associated with ill-conditioned similarity transformations. Since

$$
\mathsf {f l} (X ^ {- 1} A X) = X ^ {- 1} A X + E, \tag {7.1.14}
$$

where

$$
\| E \| _ {2} \approx \mathbf {u} \cdot \kappa_ {2} (X) \| A \| _ {2}, \tag {7.1.15}
$$

it is clear that large errors can be introduced into an eigenvalue calculation when we depart from unitary similarity.

# 7.1.6 Singular Values and Eigenvalues

Since the singular values of A and its Schur decomposition $Q ^ { H } A Q = \mathrm { d i a g } ( \lambda _ { i } ) + N$ are the same, it follows that

$$
\sigma_ {\min} (A) \leq \min _ {1 \leq i \leq n} | \lambda_ {i} | \leq \max _ {1 \leq i \leq n} | \lambda_ {i} | \leq \sigma_ {\max} (A).
$$

From what we know about the condition of triangular matrices, it may be the case that

$$
\max _ {1 \leq i, j \leq n} \frac {| \lambda_ {i} |}{| \lambda_ {j} |} \ll \kappa_ {2} (A).
$$

See §5.4.3. This is a reminder that for nonnormal matrices, eigenvalues do not have the “predictive power” of singular values when it comes to $A x = b$ sensitivity matters. Eigenvalues of nonnormal matrices have other shortcomings, a topic that is the focus of §7.9.

# Problems

P7.1.1 (a) Show that if $T \in \mathbb { C } ^ { n \times n }$ is upper triangular and normal, then T is diagonal. (b) Show that if A is normal and $Q ^ { H } A Q = T$ is a Schur decomposition, then T is diagonal. (c) Use (a) and (b) to complete the proof of Corollary 7.1.4.

P7.1.2 Prove Theorem 7.1.6 by using induction and Lemma 7.1.5.

P7.1.3 Suppose $A \in \mathbb { C } ^ { n \times n }$ has distinct eigenvalues. Show that if $Q ^ { H } A Q = T$ is its Schur decomposition and $A B = B A$ , then QH BQ is upper triangular.

P7.1.4 Show that if A and $B ^ { H }$ are in $\mathbb { C } ^ { m \times n }$ with $m \geq n$ , then

$$
\lambda (A B) = \lambda (B A) \cup \{\underbrace {0 , \ldots , 0} _ {m - n} \}.
$$

P7.1.5 Given $A \in \mathbb { C } ^ { n \times n }$ , use the Schur decomposition to show that for every $\epsilon > 0$ , there exists a diagonalizable matrix B such that $\| A - B \| _ { 2 } \leq \epsilon$ . This shows that the set of diagonalizable matrices is dense in Cn×n $\mathbb { C } ^ { n \times n }$ and that the Jordan decomposition is not a continuous matrix decomposition.

P7.1.6 Suppose $A _ { k } \to A$ and that $Q _ { k } ^ { H } A _ { k } Q _ { k } = T _ { k }$ is a Schur decomposition of $A _ { k }$ . Show that $\left\{ Q _ { k } \right\}$ has a converging subsequence $\{ Q _ { k _ { i } } \}$ with the property that

$$
\lim _ {i \to \infty} Q _ {k _ {i}} = Q
$$

where $Q ^ { H } A Q = T$ is upper triangular. This shows that the eigenvalues of a matrix are continuous functions of its entries.

P7.1.7 Justify (7.1.14) and (7.1.15).

P7.1.8 Show how to compute the eigenvalues of

$$
M = \left[ \begin{array}{c c} A & C \\ B & D \\ k & j \end{array} \right] _ {j} ^ {k}
$$

where A, B, C, and D are given real diagonal matrices.

P7.1.9 Use the Jordan decomposition to show that if all the eigenvalues of a matrix A are strictly less than unity, then $\scriptstyle \operatorname* { l i m } _ { k \to \infty } A ^ { k } = 0$ .

P7.1.10 The initial value problem

$$
\begin{array}{r c l} \dot {x} (t) & = & y (t), \qquad x (0) = 1, \\ \dot {y} (t) & = & - x (t), \qquad y (0) = 0, \end{array}
$$

has solution $x ( t ) = \cos ( t )$ and $y ( t ) = \sin ( t )$ . Let $h > 0$ . Here are three reasonable iterations that can be used to compute approximations $x _ { k } \approx x ( k h )$ and $y _ { k } \approx y ( k h )$ assuming that $x _ { 0 } = 1$ and $y _ { k } = 0 ;$ :

$$
\text { Method   1: } \quad \begin{array}{r c l} x _ {k + 1} & = & x _ {k} + h y _ {k}, \\ y _ {k + 1} & = & y _ {k} - h x _ {k}, \end{array}
$$

$$
\text { Method   2: } \quad \begin{array}{r c l} x _ {k + 1} & = & x _ {k} + h y _ {k}, \\ y _ {k + 1} & = & y _ {k} - h x _ {k + 1}, \end{array}
$$

$$
\text { Method   3: } \quad \begin{array}{r c l} x _ {k + 1} & = & x _ {k} + h y _ {k + 1}, \\ y _ {k + 1} & = & y _ {k} - h x _ {k + 1}. \end{array}
$$

Express each method in the form

$$
\left[ \begin{array}{l} x _ {k + 1} \\ y _ {k + 1} \end{array} \right] = A _ {h} \left[ \begin{array}{l} x _ {k} \\ y _ {k} \end{array} \right]
$$

where $A _ { h }$ is a 2-by-2 matrix. For each case, compute $\lambda ( A _ { h } )$ and use the previous problem to discuss lim $x _ { k }$ and lim yk as $k \to \infty$ .

P7.1.11 If $J \in \mathbb { R } ^ { d \times d }$ is a Jordan block, what is $\kappa _ { \infty } ( J ) ?$

P7.1.12 Suppose A, $B \in \mathbb { C } ^ { n \times n }$ . Show that the 2n-by-2n matrices

$$
M _ {1} = \left[ \begin{array}{c c} A B & 0 \\ B & 0 \end{array} \right] \quad \text {and} \quad M _ {2} = \left[ \begin{array}{c c} 0 & 0 \\ B & B A \end{array} \right]
$$

are similar thereby showing that $\lambda ( A B ) = \lambda ( B A )$

P7.1.13 Suppose $A \in \mathbb { R } ^ { n \times n }$ . We say that $B \in \mathbb { R } ^ { n \times n }$ is the Drazin inverse of A if (i) $A B = B A$ , (ii) $B A B = B$ , and (iii) the spectral radius of $A - A B A$ is zero. Give a formula for B in terms of the Jordan decomposition of A paying particular attention to the blocks associated with A’s zero eigenvalues.

P7.1.14 Show that if $A \in \mathbb { R } ^ { n \times n }$ , then $\rho ( A ) \geq ( \sigma _ { 1 } \cdot \cdot \cdot \sigma _ { n } ) ^ { 1 / n }$ where $\sigma _ { 1 } , \ldots , \sigma _ { n }$ are the singular values of A.

P7.1.15 Consider the polynomial $q ( x ) = \operatorname* { d e t } ( I _ { n } + x A )$ where $A \in \mathbb { R } ^ { n \times n }$ . We wish to compute the coefficient of $x ^ { 2 }$ . (a) Specify the coefficient in terms of the eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { n }$ of A. (b) Give a simple formula for the coefficient in terms of $\operatorname { t r } ( A )$ and $\operatorname { t r } ( A ^ { 2 } )$ .

P7.1.16 Given $A \in \mathbb { R } ^ { 2 \times 2 } ,$ , show that there exists a nonsingular $\begin{array} { r } { X \in \mathbb { R } ^ { 2 \times 2 } \mathrm { ~ s o ~ } X ^ { - 1 } A X = A ^ { T } } \end{array}$ . See Dubrulle and Parlett (2007).

# Notes and References for 7.1

For additional discussion about the linear algebra behind the eigenvalue problem, see Horn and Johnson (MA) and:

L. Mirsky (1963). An Introduction to Linear Algebra, Oxford University Press, Oxford, U.K.

M. Marcus and H. Minc (1964). A Survey of Matrix Theory and Matrix Inequalities, Allyn and Bacon, Boston.

R. Bellman (1970). Introduction to Matrix Analysis, second edition, McGraw-Hill, New York.

I. Gohberg, P. Lancaster, and L. Rodman (2006). Invariant Subspaces of Matrices with Applications, SIAM Publications, Philadelphia, PA.

For a general discussion about the similarity connection between a matrix and its transpose, see:

A.A. Dubrulle and B.N. Parlett (2010). “Revelations of a Transposition Matrix,” J. Comp. and Appl. Math. 233, 1217–1219.

The Schur decomposition originally appeared in:

I. Schur (1909). “On the Characteristic Roots of a Linear Substitution with an Application to the Theory of Integral Equations.” Math. Ann. 66, 488-510 (German).

A proof very similar to ours is given in:

H.W. Turnbull and A.C. Aitken (1961). An Introduction to the Theory of Canonical Forms, Dover, New York, 105.

# 7.2 Perturbation Theory

The act of computing eigenvalues is the act of computing zeros of the characteristic polynomial. Galois theory tells us that such a process has to be iterative if $n > 4$ and so errors arise because of finite termination. In order to develop intelligent stopping criteria we need an informative perturbation theory that tells us how to think about approximate eigenvalues and invariant subspaces.

# 7.2.1 Eigenvalue Sensitivity

An important framework for eigenvalue computation is to produce a sequence of similarity transformations $\{ X _ { k } \}$ with the property that the matrices $X _ { k } ^ { - 1 } A X _ { k }$ are progressively “more diagonal.” The question naturally arises, how well do the diagonal elements of a matrix approximate its eigenvalues?

Theorem 7.2.1 (Gershgorin Circle Theorem). $I f X ^ { - 1 } A X = D + F$ where $D =$ dia $\mathfrak { g } ( d _ { 1 } , \ldots , d _ { n } )$ and F has zero diagonal entries, then

$$
\lambda (A) \subseteq \bigcup_ {i = 1} ^ {n} D _ {i}
$$

$w h e r e ~ D _ { i } ~ = ~ \{ z \in \mathbb { C } : | z - d _ { i } | ~ \leq ~ \sum _ { j = 1 } ^ { n } | f _ { i j } | \} .$

Proof. Suppose $\lambda \in \lambda ( A )$ and assume without loss of generality that $\lambda \neq d _ { i }$ for $i = 1 { : } n$ . Since $( D - \lambda I ) + F$ is singular, it follows from Lemma 2.3.3 that

$$
1 \leq \| (D - \lambda I) ^ {- 1} F \| _ {\infty} = \sum_ {j = 1} ^ {n} \frac {| f _ {k j} |}{| d _ {k} - \lambda |}
$$

for some k, $1 \leq k \leq n$ . But this implies that $\lambda \in D _ { k }$ .

It can also be shown that if the Gershgorin disk $D _ { i }$ is isolated from the other disks, then it contains precisely one eigenvalue of A. See Wilkinson (AEP, pp. 71ff.).

For some methods it is possible to show that the computed eigenvalues are the exact eigenvalues of a matrix $A + E$ where E is small in norm. Consequently, we should understand how the eigenvalues of a matrix can be affected by small perturbations.

Theorem 7.2.2 (Bauer-Fike). If µ is an eigenvalue $\begin{array} { r } { o f A + E \in \mathbb { C } ^ { n \times n } \ a n d X ^ { - 1 } A X = } \end{array}$ $D = \mathrm { d i a g } ( \lambda _ { 1 } , \ldots , \lambda _ { n } )$ , then

$$
\min _ {\lambda \in \lambda (A)} | \lambda - \mu | \leq \kappa_ {p} (X) \| E \| _ {p}
$$

where $\| \cdot \| _ { p }$ denotes any of the p-norms.

Proof. If $\mu \in \lambda ( A )$ , then the theorem is obviously true. Otherwise if the matrix $X ^ { - 1 } ( A + E - \mu I ) X$ is singular, then so is $I + ( D \bar { \bf \Phi } - \mu I ) ^ { - 1 } ( X ^ { - 1 } E X )$ . Thus, from

Lemma 2.3.3 we obtain

$$
1 \leq \| (D - \mu I) ^ {- 1} (X ^ {- 1} E X) \| _ {p} \leq \| (D - \mu I) ^ {- 1} \| _ {p} \| X \| _ {p} \| E \| _ {p} \| X ^ {- 1} \| _ {p}.
$$

Since $( D - \mu I ) ^ { - 1 }$ is diagonal and the p-norm of a diagonal matrix is the absolute value of the largest diagonal entry, it follows that

$$
\| (D - \mu I) ^ {- 1} \| _ {p} = \max _ {\lambda \in \lambda (A)} \frac {1}{| \lambda - \mu |},
$$

completing the proof.

An analogous result can be obtained via the Schur decomposition:

Theorem 7.2.3. Let $Q ^ { H } A Q = D + N$ be a Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ as in $( 7 . 1 . 7 ) . \ I f \mu \in \lambda ( A + E )$ and p is the smallest positive integer such that $| N | ^ { p } = 0$ , then

$$
\min _ {\lambda \in \lambda (A)} | \lambda - \mu | \leq \max \{\theta , \theta^ {1 / p} \}
$$

where

$$
\theta = \| E \| _ {2} \sum_ {k = 0} ^ {p - 1} \| N \| _ {2} ^ {k}.
$$

Proof. Define

$$
\delta = \min _ {\lambda \in \lambda (A)} | \lambda - \mu | = \frac {1}{\| (\mu I - D) ^ {- 1} \| _ {2}}.
$$

The theorem is clearly true if $\delta = 0$ . If $\delta > 0$ , then $I - ( \mu I - A ) ^ { - 1 } E$ is singular and by Lemma 2.3.3 we have

$$
1 \leq \| (\mu I - A) ^ {- 1} E \| _ {2} \leq \| (\mu I - A) ^ {- 1} \| _ {2} \| E \| _ {2} \tag {7.2.1}
$$

$$
= \left\| \left(\left(\mu I - D\right) - N\right) ^ {- 1} \right\| _ {2} \left\| E \right\| _ {2}.
$$

Since $( \mu I - D ) ^ { - 1 }$ is diagonal and $| N | ^ { p } = 0$ , it follows that $( ( \mu I - D ) ^ { - 1 } N ) ^ { p } = 0$ . Thus,

$$
\left((\mu I - D) - N\right) ^ {- 1} = \sum_ {k = 0} ^ {p - 1} \left((\mu I - D) ^ {- 1} N\right) ^ {k} (\mu I - D) ^ {- 1}
$$

and so

$$
\| \left(\left(\mu I - D\right) - N\right) ^ {- 1} \| _ {2} \leq \frac {1}{\delta} \sum_ {k = 0} ^ {p - 1} \left(\frac {\| N \| _ {2}}{\delta}\right) ^ {k}.
$$

If $\delta > 1$ , then

$$
\parallel (\mu I - D) - N) ^ {- 1} \parallel_ {2} \leq \frac {1}{\delta} \sum_ {k = 0} ^ {p - 1} \parallel N \parallel_ {2} ^ {k}
$$

and so from (7.2.1), $\delta \leq \theta$ . If $\delta \leq 1$ , then

$$
\| (\mu I - D) - N) ^ {- 1} \| _ {2} \leq \frac {1}{\delta^ {p}} \sum_ {k = 0} ^ {p - 1} \| N \| _ {2} ^ {k}.
$$

By using (7.2.1) again we have $\delta ^ { p } \leq \theta$ and so $\delta \le \operatorname* { m a x } \{ \theta , \theta ^ { 1 / p } \}$ .

Theorems 7.2.2 and 7.2.3 suggest that the eigenvalues of a nonnormal matrix may be sensitive to perturbations. In particular, if $\kappa _ { 2 } ( X )$ or $\parallel N \parallel _ { 2 } ^ { p - 1 }$ is large, then small changes in A can induce large changes in the eigenvalues.

# 7.2.2 The Condition of a Simple Eigenvalue

Extreme eigenvalue sensitivity for a matrix A cannot occur if A is normal. On the other hand, nonnormality does not necessarily imply eigenvalue sensitivity. Indeed, a nonnormal matrix can have a mixture of well-conditioned and ill-conditioned eigenvalues. For this reason, it is beneficial to refine our perturbation theory so that it is applicable to individual eigenvalues and not the spectrum as a whole.

To this end, suppose that λ is a simple eigenvalue of $A \in \mathbb { C } ^ { n \times n }$ and that x and y satisfy $A x = \lambda x$ and $y ^ { H } A = \lambda y ^ { H }$ with $\parallel x \parallel _ { 2 } = \parallel y \parallel _ { 2 } = 1$ . If $Y ^ { H } A X = J$ is the Jordan decomposition with $Y ^ { H } = X ^ { - 1 }$ , then y and x are nonzero multiples of $X ( : , i )$ and $Y ( : , i )$ for some i. It follows from $1 = Y ( : , i ) ^ { H } X ( : , i )$ that $y ^ { H } x \neq 0$ , a fact that we shall use shortly.

Using classical results from function theory, it can be shown that in a neighborhood of the origin there exist differentiable $x ( \epsilon )$ and $\lambda ( \epsilon )$ such that

$$
(A + \epsilon F) x (\epsilon) = \lambda (\epsilon) x (\epsilon), \quad \| F \| _ {2} = 1,
$$

where $\lambda ( 0 ) = \lambda$ and $x ( 0 ) = x$ . By differentiating this equation with respect to 
 and setting $\epsilon = 0$ in the result, we obtain

$$
A \dot {x} (0) + F x = \dot {\lambda} (0) x + \lambda \dot {x} (0).
$$

Applying $y ^ { H }$ to both sides of this equation, dividing by $y ^ { H } x$ , and taking absolute values gives

$$
| \dot {\lambda} (0) | = \left| \frac {y ^ {H} F x}{y ^ {H} x} \right| \leq \frac {1}{| y ^ {H} x |}.
$$

The upper bound is attained if $F = y x ^ { H }$ . For this reason we refer to the reciprocal of

$$
s (\lambda) = | y ^ {H} x | \tag {7.2.2}
$$

as the condition of the eigenvalue λ.

Roughly speaking, the above analysis shows that $O ( \epsilon )$ perturbations in A can induce $\epsilon / s ( \lambda )$ changes in an eigenvalue. Thus, if $s ( \lambda )$ is small, then λ is appropriately regarded as ill-conditioned. Note that $s ( \lambda )$ is the cosine of the angle between the left and right eigenvectors associated with λ and is unique only if λ is simple.

A small $s ( \lambda )$ implies that A is near a matrix having a multiple eigenvalue. In particular, if λ is distinct and $s ( \lambda ) < 1$ , then there exists an E such that λ is a repeated eigenvalue of $A + E$ and

$$
\frac {\parallel E \parallel_ {2}}{\parallel A \parallel_ {2}} \leq \frac {s (\lambda)}{\sqrt {1 - s (\lambda) ^ {2}}}.
$$

This result is proved by Wilkinson (1972).

# 7.2.3 Sensitivity of Repeated Eigenvalues

If λ is a repeated eigenvalue, then the eigenvalue sensitivity question is more complicated. For example, if

$$
A = \left[ \begin{array}{c c} 1 & a \\ 0 & 1 \end{array} \right] \qquad \text { and } \qquad F = \left[ \begin{array}{c c} 0 & 0 \\ 1 & 0 \end{array} \right],
$$

then $\lambda ( A + \epsilon F ) = \{ 1 \pm \sqrt { \epsilon a } \}$ . Note that if $a \neq 0$ , then it follows that the eigenvalues of $A + \epsilon F$ are not differentiable at zero; their rate of change at the origin is infinite. In general, if λ is a defective eigenvalue of A, then $O ( \epsilon )$ perturbations in A can result in $O ( \epsilon ^ { 1 / p } )$ perturbations in λ if λ is associated with a p-dimensional Jordan block. See Wilkinson (AEP, pp. 77ff.) for a more detailed discussion.

# 7.2.4 Invariant Subspace Sensitivity

A collection of sensitive eigenvectors can define an insensitive invariant subspace provided the corresponding cluster of eigenvalues is isolated. To be precise, suppose

$$
Q ^ {H} A Q = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {n - r} ^ {r} \tag {7.2.3}
$$

is a Schur decomposition of A with

$$
Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ r & n - r \end{array} \right]. \tag {7.2.4}
$$

It is clear from our discussion of eigenvector perturbation that the sensitivity of the invariant subspace $\mathsf { r a n } ( Q _ { 1 } )$ depends on the distance between $\lambda ( T _ { 1 1 } )$ and $\lambda ( T _ { 2 2 } )$ . The proper measure of this distance turns out to be the smallest singular value of the linear transformation $X  T _ { 1 1 } X - X T _ { 2 2 }$ . (Recall that this transformation figures in Lemma 7.1.5.) In particular, if we define the separation between the matrices $T _ { 1 1 }$ and $T _ { 2 2 }$ by

$$
\mathsf {s e p} (T _ {1 1}, T _ {2 2}) = \min _ {X \neq 0} \frac {\| T _ {1 1} X - X T _ {2 2} \| _ {F}}{\| X \| _ {F}}, \tag {7.2.5}
$$

then we have the following general result:

Theorem 7.2.4. Suppose that (7.2.3) and (7.2.4) hold and that for any matrix $E \in \mathbb { C } ^ { n \times n }$ we partition $\bar { Q } ^ { H } E Q$ as follows:

$$
Q ^ {H} E Q = \left[ \begin{array}{c c} E _ {1 1} & E _ {1 2} \\ E _ {2 1} & E _ {2 2} \end{array} \right] _ {n - r} ^ {r}.
$$

$I f \mathsf { s e p } ( T _ { 1 1 } , T _ { 2 2 } ) > 0$ and

$$
\| E \| _ {F} \left(1 + \frac {5 \| T _ {1 2} \| _ {F}}{\operatorname{sep} (T _ {1 1} , T _ {2 2})}\right) \leq \frac {\operatorname{sep} (T _ {1 1} , T _ {2 2})}{5},
$$

then there exists a $P \in \mathbb { C } ^ { ( n - r ) \times r } \ w i t h$

$$
\| P \| _ {F} \leq 4 \frac {\| E _ {2 1} \| _ {F}}{\operatorname{sep} \left(T _ {1 1} , T _ {2 2}\right)}
$$

such that the columns of $\widetilde { Q } _ { 1 } = ( Q _ { 1 } + Q _ { 2 } P ) ( I + P ^ { H } P ) ^ { - 1 / 2 }$ are an orthonormal basis for a subspace invariant for $A + E$ .

Proof. This result is a slight recasting of Theorem 4.11 in Stewart (1973) which should be consulted for proof details. See also Stewart and Sun (MPA, p. 230). The matrix $( I + P ^ { H } P ) ^ { - 1 / 2 }$ is the inverse of the square root of the symmetric positive definite matrix ${ \cal I } + { \cal P } ^ { H } { \cal P }$ . See §4.2.4.

Corollary 7.2.5. If the assumptions in Theorem 7.2.4 hold, then

$$
\operatorname{dist} \left(\operatorname{ran} \left(Q _ {1}\right), \operatorname{ran} \left(\widetilde {Q} _ {1}\right)\right) \leq 4 \frac {\left\| E _ {2 1} \right\| _ {F}}{\operatorname{sep} \left(T _ {1 1} , T _ {2 2}\right)}.
$$

Proof. Using the SVD of P , it can be shown that

$$
\| P (I + P ^ {H} P) ^ {- 1 / 2} \| _ {2} \leq \| P \| _ {2} \leq \| P \| _ {F}. \tag {7.2.6}
$$

Since the required distance is the 2-norm of $Q _ { 2 } ^ { H } \widetilde { Q } _ { 1 } = P ( I + P ^ { H } P ) ^ { - 1 / 2 }$ , the proof is complete.

Thus, the reciprocal of ${ \mathsf { s e p } } ( T _ { 1 1 } , T _ { 2 2 } )$ can be thought of as a condition number that measures the sensitivity of $\mathsf { r a n } ( Q _ { 1 } )$ as an invariant subspace.

# 7.2.5 Eigenvector Sensitivity

If we set $r = 1$ in the preceding subsection, then the analysis addresses the issue of eigenvector sensitivity.

Corollary 7.2.6. Suppose A, $E \in \mathbb { C } ^ { n \times n }$ and that $Q = \left[ \left. q _ { 1 } \right| Q _ { 2 } \right] \in \mathbb { C } ^ { n \times n }$ is unitary with $q _ { 1 } \in \mathbb { C } ^ { n }$ . Assume

$$
Q ^ {H} A Q = \left[ \begin{array}{c c} \lambda & v ^ {H} \\ 0 & T _ {2 2} \end{array} \right] _ {n - 1} ^ {1}, \qquad Q ^ {H} E Q = \left[ \begin{array}{c c} \epsilon & \gamma^ {H} \\ \delta & E _ {2 2} \end{array} \right] _ {n - 1} ^ {1}.
$$

(Thus, $q _ { 1 }$ is an eigenvector.) If $\sigma = \sigma _ { \operatorname* { m i n } } ( T _ { 2 2 } - \lambda I ) > 0$ and

$$
\left\| E \right\| _ {F} \left(1 + \frac {5 \| v \| _ {2}}{\sigma}\right) \leq \frac {\sigma}{5},
$$

then there exists $p \in \mathbb { C } ^ { n - 1 }$ with

$$
\| p \| _ {2} \leq 4 \frac {\| \delta \| _ {2}}{\sigma}
$$

such that $\tilde { q } _ { 1 } = ( q _ { 1 } + Q _ { 2 } p ) / \sqrt { 1 + p ^ { H } p }$ is a unit 2-norm eigenvector for $A + E$ . Moreover,

$$
\operatorname{dist} \left(\operatorname{span} \left\{q _ {1} \right\}, \operatorname{span} \left\{\tilde {q} _ {1} \right\}\right) \leq 4 \frac {\| \delta \| _ {2}}{\sigma}.
$$

Proof. The result follows from Theorem 7.2.4, Corollary 7.2.5, and the observation that if $T _ { 1 1 } = \lambda$ , then sep $( T _ { 1 1 } , T _ { 2 2 } ) = \sigma _ { \operatorname* { m i n } } ( T _ { 2 2 } - \lambda I )$ .

Note that $\sigma _ { \operatorname* { m i n } } ( T _ { 2 2 } - \lambda I )$ roughly measures the separation of λ from the eigenvalues of $T _ { 2 2 }$ . We have to say “roughly” because

$$
\operatorname{sep} (\lambda , T _ {2 2}) = \sigma_ {\min} (T _ {2 2} - \lambda I) \leq \min _ {\mu \in \lambda (T _ {2 2})} | \mu - \lambda |
$$

and the upper bound can be a gross overestimate.

That the separation of the eigenvalues should have a bearing upon eigenvector sensitivity should come as no surprise. Indeed, if λ is a nondefective, repeated eigenvalue, then there are an infinite number of possible eigenvector bases for the associated invariant subspace. The preceding analysis merely indicates that this indeterminancy begins to be felt as the eigenvalues coalesce. In other words, the eigenvectors associated with nearby eigenvalues are “wobbly.”

# Problems

P7.2.1 Suppose $Q ^ { H } A Q = \mathrm { d i a g } ( \lambda _ { 1 } ) +$ N is a Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ and define $\nu ( A ) =$ $\parallel \boldsymbol { A } ^ { H } \boldsymbol { A } - \boldsymbol { \bar { A } } \boldsymbol { \bar { A } } ^ { H } \parallel _ { F } .$ . The upper and lower bounds in

$$
\frac {\nu (A) ^ {2}}{6 \| A \| _ {F} ^ {2}} \leq \| N \| _ {F} ^ {2} \leq \sqrt {\frac {n ^ {3} - n}{1 2}} \nu (A)
$$

are established by Henrici (1962) and Eberlein (1965), respectively. Verify these results for the case $n = 2$ .

P7.2.2 Suppose $A \in \mathbb { C } ^ { n \times n }$ and $X ^ { - 1 } A X = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ with distinct $\lambda _ { i }$ . Show that if the columns of X have unit 2-norm, then $\kappa _ { F } ( X ) ^ { 2 } = n ( 1 / s ( \bar { \lambda _ { 1 } } ) ^ { 2 } + \cdot \cdot \cdot + 1 / s ( \lambda _ { n } ) ^ { 2 } )$ .

P7.2.3 Suppose $Q ^ { H } A Q = \mathrm { d i a g } ( \lambda _ { i } ) + N$ is a Schur decomposition of A and that $X ^ { - 1 } A X = \mathrm { d i a g } \left( \lambda _ { i } \right)$ . Show κ2 $( \bar { X } ) ^ { 2 } \geq 1 + ( \| \ N \| _ { F } / \| \bar { \ A } \| _ { F } ) ^ { 2 }$ . See Loizou (1969).

P7.2.4 If $X ^ { - 1 } A X = \mathrm { { d i a g } } \left( \lambda _ { i } \right)$ and $| \lambda _ { 1 } | \geq \cdots \geq | \lambda _ { n } |$ , then

$$
\frac {\sigma_ {i} (A)}{\kappa_ {2} (X)} \leq | \lambda_ {i} | \leq \kappa_ {2} (X) \sigma_ {i} (A).
$$

Prove this result for the $n = 2$ case. See Ruhe (1975).

P7.2.5 Show that if $A = { \left[ \begin{array} { l l } { a } & { c } \\ { 0 } & { b } \end{array} \right] }$ and $a \neq b ,$ then $s ( a ) = s ( b ) = ( 1 + | c / ( a - b ) | ^ { 2 } ) ^ { - 1 / 2 }$ .

# P7.2.6 Suppose

$$
A = \left[ \begin{array}{c c} \lambda & v ^ {T} \\ 0 & T _ {2 2} \end{array} \right]
$$

and that $\lambda \not \in \lambda ( T _ { 2 2 } )$ . Show that if $\sigma = { \tt s e p } ( \lambda , T _ { 2 2 } )$ , then

$$
s (\lambda) = \frac {1}{\sqrt {1 + \| (T _ {2 2} - \lambda I) ^ {- 1} v \| _ {2} ^ {2}}} \leq \frac {\sigma}{\sqrt {\sigma^ {2} + \| v \| _ {2} ^ {2}}}.
$$

where $s ( \lambda )$ is defined in (7.2.2).

P7.2.7 Show that the condition of a simple eigenvalue is preserved under unitary similarity transformations.

P7.2.8 With the same hypothesis as in the Bauer-Fike theorem (Theorem 7.2.2), show that

$$
\min _ {\lambda \in \lambda (A)} | \lambda - \mu | \leq \| | X ^ {- 1} | | E | | X | \| _ {p}.
$$

P7.2.9 Verify (7.2.6).

P7.2.10 Show that if B ∈ Cm×m $\boldsymbol { B } \in \mathbb { C } ^ { m \times m }$ and $C \in \mathbb { C } ^ { n \times n }$ , then sep(B, C) is less than or equal to $| \lambda - \mu |$ for all $\lambda \in \lambda ( B )$ and $\mu \in \lambda ( C )$ .

# Notes and References for §7.2

Many of the results presented in this section may be found in Wilkinson (AEP), Stewart and Sun (MPA) as well as:

F.L. Bauer and C.T. Fike (1960). “Norms and Exclusion Theorems,” Numer. Math. 2, 123–44.

A.S. Householder (1964). The Theory of Matrices in Numerical Analysis. Blaisdell, New York.

R. Bhatia (2007). Perturbation Bounds for Matrix Eigenvalues, SIAM Publications, Philadelphia, PA.

Early papers concerned with the effect of perturbations on the eigenvalues of a general matrix include:

A. Ruhe (1970). “Perturbation Bounds for Means of Eigenvalues and Invariant Subspaces,” BIT 10, 343–54.

A. Ruhe (1970). “Properties of a Matrix with a Very Ill-Conditioned Eigenproblem,” Numer. Math. 15, 57–60.

J.H. Wilkinson (1972). “Note on Matrices with a Very Ill-Conditioned Eigenproblem,” Numer. Math. 19, 176–78.

W. Kahan, B.N. Parlett, and E. Jiang (1982). “Residual Bounds on Approximate Eigensystems of Nonnormal Matrices,” SIAM J. Numer. Anal. 19, 470–484.

J.H. Wilkinson (1984). “On Neighboring Matrices with Quadratic Elementary Divisors,” Numer. Math. 44, 1-21.

Wilkinson’s work on nearest defective matrices is typical of a growing body of literature that is concerned with “nearness” problems, see:

A. Ruhe (1987). “Closest Normal Matrix Found!,” BIT 27, 585-598.

J.W. Demmel (1987). “On the Distance to the Nearest Ill-Posed Problem,” Numer. Math. 51, 251–289.

J.W. Demmel (1988). “The Probability that a Numerical Analysis Problem is Difficult,” Math. Comput. 50, 449–480.

N.J. Higham (1989). “Matrix Nearness Problems and Applications,” in Applications of Matrix Theory, M.J.C. Gover and S. Barnett (eds.), Oxford University Press, Oxford, 1–27.

A.N. Malyshev (1999). “A Formula for the 2-norm Distance from a Matrix to the Set of Matrices with Multiple Eigenvalues,” Numer. Math. 83, 443–454.

J.-M. Gracia (2005). “Nearest Matrix with Two Prescribed Eigenvalues,” Lin. Alg. Applic. 401, 277–294.

An important subset of this literature is concerned with nearness to the set of unstable matrices. A matrix is unstable if it has an eigenvalue with nonnegative real part. Controllability is a related notion, see:

C. Van Loan (1985). “How Near is a Stable Matrix to an Unstable Matrix?,” Contemp. Math. 47, 465–477.   
J.W. Demmel (1987). “A Counterexample for two Conjectures About Stability,” IEEE Trans. Autom. Contr. AC-32, 340–342.   
R. Byers (1988). “A Bisection Method for Measuring the distance of a Stable Matrix to the Unstable Matrices,” J. Sci. Stat. Comput. 9, 875–881.   
J.V. Burke and M.L. Overton (1992). “Stable Perturbations of Nonsymmetric Matrices,” Lin. Alg. Applic. 171, 249–273.   
C. He and G.A. Watson (1998). “An Algorithm for Computing the Distance to Instability,” SIAM J. Matrix Anal. Applic. 20, 101–116.   
M. Gu, E. Mengi, M.L. Overton, J. Xia, and J. Zhu (2006). “Fast Methods for Estimating the Distance to Uncontrollability,” SIAM J. Matrix Anal. Applic. 28, 477–502.   
Aspects of eigenvalue condition are discussed in:   
C. Van Loan (1987). “On Estimating the Condition of Eigenvalues and Eigenvectors,” Lin. Alg. Applic. 88/89, 715–732.   
C.D. Meyer and G.W. Stewart (1988). “Derivatives and Perturbations of Eigenvectors,” SIAM J. Numer. Anal. 25, 679–691.   
G.W. Stewart and G. Zhang (1991). “Eigenvalues of Graded Matrices and the Condition Numbers of Multiple Eigenvalues,” Numer. Math. 58, 703–712.   
J.-G. Sun (1992). “On Condition Numbers of a Nondefective Multiple Eigenvalue,” Numer. Math. 61, 265–276.   
S.M. Rump (2001). “Computational Error Bounds for Multiple or Nearly Multiple Eigenvalues,” Lin. Alg. Applic. 324, 209–226.   
The relationship between the eigenvalue condition number, the departure from normality, and the condition of the eigenvector matrix is discussed in:   
P. Henrici (1962). “Bounds for Iterates, Inverses, Spectral Variation and Fields of Values of Nonnormal Matrices,” Numer. Math. 4, 24–40.   
P. Eberlein (1965). “On Measures of Non-Normality for Matrices,” AMS Monthly 72, 995–996.   
R.A. Smith (1967). “The Condition Numbers of the Matrix Eigenvalue Problem,” Numer. Math. 10 232–240.   
G. Loizou (1969). “Nonnormality and Jordan Condition Numbers of Matrices,” J. ACM 16, 580–640.   
A. van der Sluis (1975). “Perturbations of Eigenvalues of Non-normal Matrices,” Commun. ACM 18, 30–36.   
S.L. Lee (1995). “A Practical Upper Bound for Departure from Normality,” SIAM J. Matrix Anal. Applic. 16, 462–468.   
Gershgorin’s theorem can be used to derive a comprehensive perturbation theory. The theorem itself can be generalized and extended in various ways, see:   
R.S. Varga (1970). “Minimal Gershgorin Sets for Partitioned Matrices,” SIAM J. Numer. Anal. 7, 493–507.   
R.J. Johnston (1971). “Gershgorin Theorems for Partitioned Matrices,” Lin. Alg. Applic. 4, 205–20.   
R.S. Varga and A. Krautstengl (1999). “On Gergorin-type Problems and Ovals of Cassini,” ETNA 8, 15–20.   
R.S. Varga (2001). “Gergorin-type Eigenvalue Inclusion Theorems and Their Sharpness,” ETNA 12, 113–133.   
C. Beattie and I.C.F. Ipsen (2003). “Inclusion Regions for Matrix Eigenvalues,” Lin. Alg. Applic. 358, 281–291.   
In our discussion, the perturbations to the A-matrix are general. More can be said when the perturbations are structured, see:   
G.W. Stewart (2001). “On the Eigensystems of Graded Matrices,” Numer. Math. 90, 349–370.   
J. Moro and F.M. Dopico (2003). “Low Rank Perturbation of Jordan Structure,” SIAM J. Matrix Anal. Applic. 25, 495–506.   
R. Byers and D. Kressner (2004). “On the Condition of a Complex Eigenvalue under Real Perturbations,” BIT 44, 209–214.   
R. Byers and D. Kressner (2006). “Structured Condition Numbers for Invariant Subspaces,” SIAM J. Matrix Anal. Applic. 28, 326–347.

An absolute perturbation bound comments on the difference between an eigenvalue λ and its perturbation λ˜. A relative perturbation bound examines the quotient $| \lambda - { \tilde { \lambda } } | / | { \bar { \lambda } } | ,$ , something that can be very important when there is a concern about a small eigenvalue. For results in this direction consult:

R.-C. Li (1997). “Relative Perturbation Theory. III. More Bounds on Eigenvalue Variation,” Lin. Alg. Applic. 266, 337–345.

S.C. Eisenstat and I.C.F. Ipsen (1998). “Three Absolute Perturbation Bounds for Matrix Eigenvalues Imply Relative Bounds,” SIAM J. Matrix Anal. Applic. 20, 149–158.

S.C. Eisenstat and I.C.F. Ipsen (1998). “Relative Perturbation Results for Eigenvalues and Eigenvectors of Diagonalisable Matrices,” BIT 38, 502–509.

I.C.F. Ipsen (1998). “Relative Perturbation Results for Matrix Eigenvalues and Singular Values,” Acta Numerica, 7, 151–201.

I.C.F. Ipsen (2000). “Absolute and Relative Perturbation Bounds for Invariant Subspaces of Matrices,” Lin. Alg. Applic. 309, 45–56.

I.C.F. Ipsen (2003). “A Note on Unifying Absolute and Relative Perturbation Bounds,” Lin. Alg. Applic. 358, 239–253.

Y. Wei, X. Li, F. Bu, and F. Zhang (2006). “Relative Perturbation Bounds for the Eigenvalues of Diagonalizable and Singular Matrices–Application to Perturbation Theory for Simple Invariant Subspaces,” Lin. Alg. Applic. 419, 765-771.

The eigenvectors and invariant subspaces of a matrix also “move” when there are perturbations. Tracking these changes is typically more challenging than tracking changes in the eigenvalues, see:

T. Kato (1966). Perturbation Theory for Linear Operators, Springer-Verlag, New York.

C. Davis and W.M. Kahan (1970). “The Rotation of Eigenvectors by a Perturbation, III,” SIAM J. Numer. Anal. 7, 1–46.

G.W. Stewart (1971). “Error Bounds for Approximate Invariant Subspaces of Closed Linear Operators,” SIAM. J. Numer. Anal. 8, 796–808.

G.W. Stewart (1973). “Error and Perturbation Bounds for Subspaces Associated with Certain Eigenvalue Problems,” SIAM Review 15, 727–764.

J. Xie (1997). “A Note on the Davis-Kahan sin(2θ) Theorem,” Lin. Alg. Applic. 258, 129–135.

S.M. Rump and J.-P.M. Zemke (2003). “On Eigenvector Bounds,” BIT 43, 823–837.

Detailed analyses of the function sep(.,.) and the map $X  A X + X A ^ { T }$ are given in:

J. Varah (1979). “On the Separation of Two Matrices,” SIAM J. Numer. Anal. 16, 216–22.

R. Byers and S.G. Nash (1987). “On the Singular Vectors of the Lyapunov Operator,” SIAM J. Alg. Disc. Methods 8, 59–66.

# 7.3 Power Iterations

Suppose that we are given $A \in \mathbb { C } ^ { n \times n }$ and a unitary $U _ { 0 } \in \mathbb { C } ^ { n \times n }$ . Recall from §5.2.10 that the Householder QR factorization can be extended to complex matrices and consider the following iteration:

$$
T _ {0} = U _ {0} ^ {H} A U _ {0}
$$

for k = 1, 2, . . .

$$
T _ {k - 1} = U _ {k} R _ {k} \quad (\text { QR   factorization }) \tag {7.3.1}
$$

$$
T _ {k} = R _ {k} U _ {k}
$$

end

Since $T _ { k } = R _ { k } U _ { k } = U _ { k } ^ { H } ( U _ { k } R _ { k } ) U _ { k } = U _ { k } ^ { H } T _ { k - 1 } U _ { k }$ it follows by induction that

$$
T _ {k} = (U _ {0} U _ {1} \dots U _ {k}) ^ {H} A (U _ {0} U _ {1} \dots U _ {k}). \tag {7.3.2}
$$

Thus, each $T _ { k }$ is unitarily similar to A. Not so obvious, and what is a central theme of this section, is that the $T _ { k }$ almost always converge to upper triangular form, i.e., (7.3.2) almost always “converges” to a Schur decomposition of A.

Iteration (7.3.1) is called the QR iteration, and it forms the backbone of the most effective algorithm for computing a complete Schur decomposition of a dense general matrix. In order to motivate the method and to derive its convergence properties, two other eigenvalue iterations that are important in their own right are presented first: the power method and the method of orthogonal iteration.

# 7.3.1 The Power Method

Suppose $A \in \mathbb { C } ^ { n \times n }$ and $X ^ { - 1 } A X = \operatorname { d i a g } ( \lambda _ { 1 } , . . . , \lambda _ { n } )$ with $X = \left[ \left. x _ { 1 } \right| \cdot \cdot \cdot \right| \left. x _ { n } \right]$ . Assume that

$$
\left| \lambda_ {1} \right| > \left| \lambda_ {2} \right| \geq \dots \geq \left| \lambda_ {n} \right|.
$$

Given a unit 2-norm $q ^ { ( 0 ) } \in \mathbb { C } ^ { n }$ , the power method produces a sequence of vectors $q ^ { ( k ) }$ as follows:

for k = 1, 2, . . .

$$
z ^ {(k)} = A q ^ {(k - 1)}
$$

$$
q ^ {(k)} = z ^ {(k)} / \| z ^ {(k)} \| _ {2} \tag {7.3.3}
$$

$$
\lambda^ {(k)} = [ q ^ {(k)} ] ^ {H} A q ^ {(k)}
$$

end

There is nothing special about using the 2-norm for normalization except that it imparts a greater unity on the overall discussion in this section.

Let us examine the convergence properties of the power iteration. If

$$
q ^ {(0)} = a _ {1} x _ {1} + a _ {2} x _ {2} + \dots + a _ {n} x _ {n} \tag {7.3.4}
$$

and $a _ { 1 } \neq 0$ , then

$$
A ^ {k} q ^ {(0)} = a _ {1} \lambda_ {1} ^ {k} \left(x _ {1} + \sum_ {j = 2} ^ {n} \frac {a _ {j}}{a _ {1}} \left(\frac {\lambda_ {j}}{\lambda_ {1}}\right) ^ {k} x _ {j}\right).
$$

Since $q ^ { ( k ) } \in \mathsf { s p a n } \{ A ^ { k } q ^ { ( 0 ) } \}$ we conclude that

$$
\operatorname{dist} \left(\operatorname{span} \{q ^ {(k)} \}, \operatorname{span} \{x _ {1} \}\right) = O \left(\left| \frac {\lambda_ {2}}{\lambda_ {1}} \right| ^ {k}\right).
$$

It is also easy to verify that

$$
\left| \lambda_ {1} - \lambda^ {(k)} \right| = O \left(\left| \frac {\lambda_ {2}}{\lambda_ {1}} \right| ^ {k}\right). \tag {7.3.5}
$$

Since $\lambda _ { 1 }$ is larger than all the other eigenvalues in modulus, it is referred to as a dominant eigenvalue. Thus, the power method converges if $\lambda _ { 1 }$ is dominant and if $q ^ { ( 0 ) }$ has a component in the direction of the corresponding dominant eigenvector $x _ { 1 }$ . The behavior of the iteration without these assumptions is discussed in Wilkinson (AEP, p. 570) and Parlett and Poole (1973).

In practice, the usefulness of the power method depends upon the ratio $| \lambda _ { 2 } | / | \lambda _ { 1 } |$ , since it dictates the rate of convergence. The danger that $q ^ { ( 0 ) }$ is deficient in $x _ { 1 }$ is less worrisome because rounding errors sustained during the iteration typically ensure that subsequent iterates have a component in this direction. Moreover, it is typically the case in applications that one has a reasonably good guess as to the direction of $x _ { 1 }$ . This guards against having a pathologically small coefficient $a _ { 1 }$ in (7.3.4).

Note that the only thing required to implement the power method is a procedure for matrix-vector products. It is not necessary to store A in an n-by-n array. For this reason, the algorithm is of interest when the dominant eigenpair for a large sparse matrix is required. We have much more to say about large sparse eigenvalue problems in Chapter 10.

Estimates for the error $| \lambda ^ { ( k ) } - \lambda _ { 1 } |$ can be obtained by applying the perturbation theory developed in §7.2.2. Define the vector

$$
r ^ {(k)} = A q ^ {(k)} - \lambda^ {(k)} q ^ {(k)}
$$

and observe that $( A + E ^ { ( k ) } ) q ^ { ( k ) } \ = \ \lambda ^ { ( k ) } q ^ { ( k ) }$ where $E ^ { ( k ) } = - r ^ { ( k ) } [ q ^ { ( k ) } ] ^ { H }$ . Thus $\lambda ^ { ( k ) }$ is an eigenvalue of $A + E ^ { ( k ) }$ and

$$
| \lambda^ {(k)} - \lambda_ {1} | \approx \frac {\parallel E ^ {(k)} \parallel_ {2}}{s (\lambda_ {1})} = \frac {\parallel r ^ {(k)} \parallel_ {2}}{s (\lambda_ {1})}.
$$

If we use the power method to generate approximate right and left dominant eigenvectors, then it is possible to obtain an estimate of $s ( \lambda _ { 1 } )$ . In particular, if $w ^ { ( k ) }$ is a unit 2-norm vector in the direction of $( A ^ { H } ) ^ { k } w ^ { ( 0 ) }$ , then we can use the approximation s(λ1) ≈ | w(k)H q $s ( \lambda _ { 1 } ) \approx | \boldsymbol { w } ^ { ( k ) } \boldsymbol { q } ^ { \prime } ( k ) |$ .

# 7.3.2 Orthogonal Iteration

A straightforward generalization of the power method can be used to compute higherdimensional invariant subspaces. Let r be a chosen integer satisfying $1 \leq r \leq n$ . Given $A \in \mathbb { C } ^ { n \times n }$ and an n-by-r matrix $Q _ { 0 }$ with orthonormal columns, the method of orthogonal iteration generates a sequence of matrices $\{ Q _ { k } \} \subseteq \mathbb { C } ^ { n \times r }$ and a sequence of eigenvalue estimates $\big \{ \lambda _ { 1 } ^ { ( k ) } , \dots , \lambda _ { r } ^ { ( k ) } \big \}$ as follows:

for $k = 1 , 2 , \dots$

$$
Z _ {k} = A Q _ {k - 1}
$$

$$
Q _ {k} R _ {k} = Z _ {k} \quad \text {(QR factorization)} \tag {7.3.6}
$$

$$
\lambda (Q _ {k} ^ {H} A Q _ {k}) = \{\lambda_ {1} ^ {(k)}, \dots , \lambda_ {r} ^ {(k)} \}
$$

end

Note that if $r \ = \ 1$ , then this is just the power method (7.3.3). Moreover, the sequence $\{ Q _ { k } e _ { 1 } \}$ is precisely the sequence of vectors produced by the power iteration with starting vector $q ^ { ( 0 ) } = Q _ { 0 } e _ { 1 }$ .

In order to analyze the behavior of this iteration, suppose that

$$
Q ^ {H} A Q = T = \operatorname{diag} \left(\lambda_ {i}\right) + N, \quad \left| \lambda_ {1} \right| \geq \left| \lambda_ {2} \right| \geq \dots \geq \left| \lambda_ {n} \right| \tag {7.3.7}
$$

is a Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ . Assume that $1 \leq r < n$ and partition Q and $T$ as follows:

$$
Q = \left[ \begin{array}{c c} Q _ {\alpha} & Q _ {\beta} \\ r & n - r \end{array} \right], \quad T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \\ r & n - r \end{array} \right] _ {n - r} ^ {r}. \tag {7.3.8}
$$

If $\left| \lambda _ { r } \right| > \left| \lambda _ { r + 1 } \right|$ , then the subspace $D _ { r } ( A ) \ = \ r { \mathsf { a n } } ( Q _ { \alpha } )$ is referred to as a dominant invariant subspace. It is the unique invariant subspace associated with the eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { r }$ . The following theorem shows that with reasonable assumptions, the subspaces $\mathsf { r a n } ( Q _ { k } )$ generated by (7.3.6) converge to $D _ { r } ( A )$ at a rate proportional to $| \lambda _ { r + 1 } / \lambda _ { r } | ^ { k }$ .

Theorem 7.3.1. Let the Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ be given by (7.3.7) and (7.3.8) with $n \geq 2$ . Assume that $| \lambda _ { r } | > | \lambda _ { r + 1 } |$ and that $\mu \geq 0$ satisfies

$$
(1 + \mu) | \lambda_ {r} | > \| N \| _ {F}.
$$

Suppose $Q _ { 0 } \in \mathbb { C } ^ { n \times r }$ has orthonormal columns and that $d _ { k }$ is defined by

$$
d _ {k} = \operatorname{dist} \left(D _ {r} (A), \operatorname{ran} \left(Q _ {k}\right)\right), \quad k \geq 0.
$$

If

$$
d _ {0} <   1, \tag {7.3.9}
$$

then the matrices $Q _ { k }$ generated by (7.3.6) satisfy

$$
d _ {k} \leq (1 + \mu) ^ {n - 2} \cdot \left(1 + \frac {\| T _ {1 2} \| _ {F}}{\operatorname{sep} \left(T _ {1 1} , T _ {2 2}\right)}\right) \cdot \left[ \frac {\left| \lambda_ {r + 1} \right| + \frac {\| N \| _ {F}}{1 + \mu}}{\left| \lambda_ {r} \right| - \frac {\| N \| _ {F}}{1 + \mu}} \right] ^ {k} \cdot \frac {d _ {0}}{\sqrt {1 - d _ {0} ^ {2}}}. \tag {7.3.10}
$$

Proof. The proof is given in an appendix at the end of this section.

The condition (7.3.9) ensures that the initial matrix $Q _ { 0 }$ is not deficient in certain eigendirections. In particular, no vector in the span of $Q _ { 0 } \mathrm { { ^ { \circ } s } }$ columns is orthogonal to $\bar { D _ { r } } ( A ^ { H } )$ . The theorem essentially says that if this condition holds and if $\mu$ is chosen large enough, then

$$
\operatorname{dist} \left(D _ {r} (A), \operatorname{ran} \left(Q _ {k}\right)\right) \approx c \left| \frac {\lambda_ {r + 1}}{\lambda_ {r}} \right| ^ {k}
$$

where c depends on sep $( T _ { 1 1 } , T _ { 2 2 } )$ and A’s departure from normality.

It is possible to accelerate the convergence in orthogonal iteration using a technique described in Stewart (1976). In the accelerated scheme, the approximate eigenvalue ${ \lambda } _ { i } ^ { ( k ) }$ satisfies

$$
\left| \lambda_ {i} ^ {(k)} - \lambda_ {i} \right| \approx \left| \frac {\lambda_ {r + 1}}{\lambda_ {i}} \right| ^ {k}, \quad i = 1: r.
$$

(Without the acceleration, the right-hand side is $| \lambda _ { i + 1 } / \lambda _ { i } | ^ { k } . \big )$ Stewart’s algorithm involves computing the Schur decomposition of the matrices $Q _ { k } ^ { T } A Q _ { k }$ every so often. The method can be very useful in situations where A is large and sparse and a few of its largest eigenvalues are required.

# 7.3.3 The QR Iteration

We now derive the QR iteration (7.3.1) and examine its convergence. Suppose $r = n$ in (7.3.6) and the eigenvalues of A satisfy

$$
\left| \lambda_ {1} \right| > \left| \lambda_ {2} \right| > \dots > \left| \lambda_ {n} \right|.
$$

Partition the matrix Q in (7.3.7) and $Q _ { k }$ in (7.3.6) as follows:

$$
Q = \left[ q _ {1} \mid \dots \mid q _ {n} \right], \quad Q _ {k} = \left[ q _ {1} ^ {(k)} \mid \dots \mid q _ {n} ^ {(k)} \right].
$$

If

$$
\mathsf {d i s t} (D _ {i} (A ^ {H}), \mathsf {s p a n} \{q _ {1} ^ {(0)}, \dots , q _ {i} ^ {(0)} \}) <   1, \quad i = 1: n, \tag {7.3.11}
$$

then it follows from Theorem 7.3.1 that

$$
\operatorname{dist} \left(\operatorname{span} \left\{q _ {1} ^ {(k)}, \dots , q _ {i} ^ {(k)} \right\}, \operatorname{span} \left\{q _ {1}, \dots , q _ {i} \right\}\right)\rightarrow 0
$$

for i = 1:n. This implies that the matrices $T _ { k }$ defined by

$$
T _ {k} = Q _ {k} ^ {H} A Q _ {k}
$$

are converging to upper triangular form. Thus, it can be said that the method of orthogonal iteration computes a Schur decomposition provided the original iterate $Q _ { 0 } \in \mathbb { C } ^ { n \times n }$ is not deficient in the sense of (7.3.11).

The QR iteration arises naturally by considering how to compute the matrix $T _ { k }$ directly from its predecessor $T _ { k - 1 }$ . On the one hand, we have from (7.3.6) and the definition of $T _ { k - 1 }$ that

$$
T _ {k - 1} = Q _ {k - 1} ^ {H} A Q _ {k - 1} = Q _ {k - 1} ^ {H} (A Q _ {k - 1}) = (Q _ {k - 1} ^ {H} Q _ {k}) R _ {k}.
$$

On the other hand,

$$
T _ {k} = Q _ {k} ^ {H} A Q _ {k} = (Q _ {k} ^ {H} A Q _ {k - 1}) (Q _ {k - 1} ^ {H} Q _ {k}) = R _ {k} (Q _ {k - 1} ^ {H} Q _ {k}).
$$

Thus, $T _ { k }$ is determined by computing the QR factorization of $T _ { k - 1 }$ and then multiplying the factors together in reverse order, precisely what is done in (7.3.1).

Note that a single QR iteration is an $O ( n ^ { 3 } )$ calculation. Moreover, since convergence is only linear (when it exists), it is clear that the method is a prohibitively expensive way to compute Schur decompositions. Fortunately these practical difficulties can be overcome as we show in §7.4 and §7.5.

# 7.3.4 LR Iterations

We conclude with some remarks about power iterations that rely on the LU factorization rather than the QR factorizaton. Let $G _ { 0 } \in \mathbb { C } ^ { n \times r }$ have rank r. Corresponding to (7.3.1) we have the following iteration:

for $k = 1 , 2 , \dots$

$$
Z _ {k} = A G _ {k - 1} \tag {7.3.12}
$$

$$
Z _ {k} = G _ {k} R _ {k} \quad \text {(LU factorization)}
$$

end

Suppose $r = n$ and that we define the matrices $T _ { k }$ by

$$
T _ {k} = G _ {k} ^ {- 1} A G _ {k}. \tag {7.3.13}
$$

It can be shown that if we set $L _ { 0 } = G _ { 0 }$ , then the $T _ { k }$ can be generated as follows:

$$
T _ {0} = L _ {0} ^ {- 1} A L _ {0}
$$

for $k = 1 , 2 ,$ . . .

$$
T _ {k - 1} = L _ {k} R _ {k} \quad (\text {LU factorization}) \tag {7.3.14}
$$

$$
T _ {k} = R _ {k} L _ {k}
$$

end

Iterations (7.3.12) and (7.3.14) are known as treppeniteration and the LR iteration, respectively. Under reasonable assumptions, the $T _ { k }$ converge to upper triangular form. To successfully implement either method, it is necessary to pivot. See Wilkinson (AEP, p. 602).

# Appendix

In order to establish Theorem 7.3.1 we need the following lemma that bounds powers of a matrix and powers of its inverse.

Lemma 7.3.2. Let $Q ^ { H } A Q = T = D + N $ be a Schur decomposition of $A \in \mathbb { C } ^ { n \times n }$ where D is diagonal and N strictly upper triangular. Let $\lambda _ { \mathrm { m a x } }$ and $\lambda _ { \mathrm { m i n } }$ denote the largest and smallest eigenvalues of A in absolute value. If $\mu \geq 0$ , then for all $k \geq 0$ we have

$$
\| A ^ {k} \| _ {2} \leq (1 + \mu) ^ {n - 1} \left(\left| \lambda_ {\max} \right| + \frac {\left\| N \right\| _ {F}}{1 + \mu}\right) ^ {k}. \tag {7.3.15}
$$

If A is nonsingular and $\mu \geq 0$ satisfies $( 1 + \mu ) | \lambda _ { \operatorname* { m i n } } | > \parallel N \parallel _ { _ F }$ , then for all $k \geq 0$ we also have

$$
\| A ^ {- k} \| _ {2} \leq (1 + \mu) ^ {n - 1} \left(\frac {1}{\left| \lambda_ {\min} \right| - \| N \| _ {F} / (1 + \mu)}\right) ^ {k}. \tag {7.3.16}
$$

Proof. For $\mu \geq 0$ , define the diagonal matrix $\Delta$ by

$$
\Delta = \operatorname{diag} (1, (1 + \mu), (1 + \mu) ^ {2}, \dots , (1 + \mu) ^ {n - 1})
$$

and note that $\kappa _ { 2 } ( \Delta ) = ( 1 + \mu ) ^ { n - 1 }$ . Since N is strictly upper triangular, it is easy to verify that

$$
\| \Delta N \Delta^ {- 1} \| _ {F} \leq \frac {\| N \| _ {F}}{1 + \mu}
$$

and thus

$$
\begin{array}{l} \| A ^ {k} \| _ {2} = \| T ^ {k} \| _ {2} = \| \Delta^ {- 1} (D + \Delta N \Delta^ {- 1}) ^ {k} \Delta \| _ {2} \\ \leq \kappa_ {2} (\Delta) \left(\| D \| _ {2} + \| \Delta N \Delta^ {- 1} \| _ {2}\right) ^ {k} \leq (1 + \mu) ^ {n - 1} \left(| \lambda_ {\max} | + \frac {\| N \| _ {F}}{1 + \mu}\right) ^ {k}. \\ \end{array}
$$

On the other hand, if A is nonsingular and $( 1 + \mu ) | \lambda _ { \operatorname* { m i n } } | > \parallel N \parallel _ { F }$ , then

$$
\| \Delta D ^ {- 1} N \Delta^ {- 1} \| _ {2} = \| D ^ {- 1} (\Delta N \Delta^ {- 1}) \| _ {2} \leq \frac {1}{| \lambda_ {\min} |} \| \Delta N \Delta^ {- 1} \| _ {F} <   1.
$$

Using Lemma 2.3.3 we obtain

$$
\begin{array}{l} \| A ^ {- k} \| _ {2} = \| T ^ {- k} \| _ {2} = \left\| \Delta^ {- 1} [ (I + \Delta D ^ {- 1} N \Delta^ {- 1}) ^ {- 1} D ^ {- 1} ] ^ {k} \Delta \right\| _ {2} \\ \leq \kappa_ {2} (\Delta) \left(\frac {\| D ^ {- 1} \| _ {2}}{1 - \| \Delta D ^ {- 1} N \Delta^ {- 1} \| _ {2}}\right) ^ {k} \leq (1 + \mu) ^ {n - 1} \left(\frac {1}{| \mu | - \| N \| _ {F} / (1 + \mu)}\right) ^ {k} \\ \end{array}
$$

completing the proof of the lemma.

Proof of Theorem 7.3.1. By induction it is easy to show that the matrix $Q _ { k }$ in (7.3.6) satisfies

$$
A ^ {k} Q _ {0} = Q _ {k} (R _ {k} \dots R _ {1}),
$$

a QR factorization of $A ^ { k } Q _ { 0 }$ . By substituting the Schur decomposition (7.3.7)-(7.3.8) into this equation we obtain

$$
T ^ {k} \left[ \begin{array}{c} V _ {0} \\ W _ {0} \end{array} \right] = \left[ \begin{array}{c} V _ {k} \\ W _ {k} \end{array} \right] (R _ {k} \dots R _ {1}) \tag {7.3.17}
$$

where

$$
V _ {k} = Q _ {\alpha} ^ {H} Q _ {k}, \qquad W _ {k} = Q _ {\beta} ^ {H} Q _ {k}.
$$

Our goal is to bound $\parallel W _ { k } \parallel _ { 2 }$ since by the definition of subspace distance given in §2.5.3 we have

$$
\| W _ {k} \| _ {2} = \operatorname{dist} (D _ {r} (A), \operatorname{ran} (Q _ {k})). \tag {7.3.18}
$$

Note from the thin CS decomposition (Theorem 2.5.2) that

$$
1 = d _ {k} ^ {2} + \sigma_ {\min} (V _ {k}) ^ {2}. \tag {7.3.19}
$$

Since $T _ { 1 1 }$ and $T _ { 2 2 }$ have no eigenvalues in common, Lemma 7.1.5 tells us that the Sylvester equation $T _ { 1 1 } X \ : - \ : X T _ { 2 2 } \ : = \ : - T _ { 1 2 }$ has a solution $X \in \mathbb { C } ^ { r \times ( n - r ) }$ and that

$$
\| X \| _ {F} \leq \frac {\| T _ {1 2} \| _ {F}}{\operatorname{sep} \left(T _ {1 1} , T _ {2 2}\right)}. \tag {7.3.20}
$$

It follows that

$$
\left[ \begin{array}{c c} I _ {r} & X \\ 0 & I _ {n - r} \end{array} \right] ^ {- 1} \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] \left[ \begin{array}{c c} I _ {r} & X \\ 0 & I _ {n - r} \end{array} \right] = \left[ \begin{array}{c c} T _ {1 1} & 0 \\ 0 & T _ {2 2} \end{array} \right].
$$

By substituting this into (7.3.17) we obtain

$$
\left[ \begin{array}{c c} T _ {1 1} ^ {k} & 0 \\ 0 & T _ {2 2} ^ {k} \end{array} \right] \left[ \begin{array}{c} V _ {0} - X W _ {0} \\ W _ {0} \end{array} \right] = \left[ \begin{array}{c} V _ {k} - X W _ {k} \\ W _ {k} \end{array} \right] (R _ {k} \dots R _ {1}),
$$

$$
T _ {1 1} ^ {k} (V _ {0} - X W _ {0}) = (V _ {k} - X W _ {k}) (R _ {k} \dots R _ {1}), \tag {7.3.21}
$$

$$
T _ {2 2} ^ {k} W _ {0} = W _ {k} (R _ {k} \dots R _ {1}). \tag {7.3.22}
$$

The matrix $I + X X ^ { H }$ is Hermitian positive definite and so it has a Cholesky factorization

$$
I + X X ^ {H} = G G ^ {H}. \tag {7.3.23}
$$

It is clear that

$$
\sigma_ {\min} (G) \geq 1. \tag {7.3.24}
$$

If the matrix $Z \in \mathbb { C } ^ { n \times ( n - r ) }$ is defined by

$$
Z = Q \left[ \begin{array}{c} I _ {r} \\ - X ^ {H} \end{array} \right] G ^ {- H} = \left[ Q _ {\alpha} Q _ {\beta} \right] \left[ \begin{array}{c} I _ {r} \\ - X ^ {H} \end{array} \right] G ^ {- H} = (Q _ {\alpha} - Q _ {\beta} X ^ {H}) G ^ {- H},
$$

then it follows from the equation $A ^ { H } Q = Q T ^ { H }$ that

$$
A ^ {H} (Q _ {\alpha} - Q _ {\beta} X ^ {H}) = (Q _ {\alpha} - Q _ {\beta} X ^ {H}) T _ {1 1} ^ {H}. \tag {7.3.25}
$$

Since $Z ^ { H } Z = I _ { r }$ and ran $( Z ) \ : = \ : \mathsf { r a n } ( Q _ { \alpha } - Q _ { \beta } X ^ { H } )$ , it follows that the columns of Z are an orthonormal basis for $D _ { r } ( A ^ { H } )$ . Using the CS decomposition, (7.3.19), and the fact that ran $( Q _ { \beta } ) = D _ { r } ( A ^ { H } ) ^ { \perp }$ , we have

$$
\begin{array}{l} \sigma_ {\min} (Z ^ {T} Q _ {0}) ^ {2} = 1 - \operatorname{dist} (D _ {r} (A ^ {H}), Q _ {0}) ^ {2} = 1 - \| Q _ {\beta} ^ {H} Q _ {0} \| \\ = \sigma_ {\mathrm{min}} (Q _ {\alpha} ^ {T} Q _ {0}) ^ {2} = \sigma_ {\mathrm{min}} (V _ {0}) ^ {2} = 1 - d _ {0} ^ {2} > 0. \\ \end{array}
$$

This shows that

$$
V _ {0} - X W _ {0} = \left[ \begin{array}{l} I _ {r} \end{array} \right| - X \left. \right] \left[ \begin{array}{l} Q _ {\alpha} ^ {H} Q _ {0} \\ Q _ {\beta} ^ {H} Q _ {0} \end{array} \right] = (Z G ^ {H}) ^ {H} Q _ {0} = G (Z ^ {H} Q _ {0})
$$

is nonsingular and together with (7.3.24) we obtain

$$
\| (V _ {0} - X W _ {0}) ^ {- 1} \| _ {2} \leq \| G ^ {- 1} \| _ {2} \| (Z ^ {H} Q _ {0}) ^ {- 1} \| _ {2} \leq \frac {1}{\sqrt {1 - d _ {0} ^ {2}}}. \tag {7.3.26}
$$

Manipulation of (7.3.19) and (7.3.20) yields

$$
W _ {k} = T _ {2 2} ^ {k} W _ {0} (R _ {k} \dots R _ {1}) ^ {- 1} = T _ {2 2} ^ {k} W _ {0} (V _ {0} - X W _ {0}) ^ {- 1} T _ {1 1} ^ {- k} (V _ {k} - X W _ {k}).
$$

The verification of (7.3.10) is completed by taking norms in this equation and using (7.3.18), (7.3.19), (7.3.20), (7.3.26), and the following facts:

$$
\| T _ {2 2} ^ {k} \| _ {2} \leq (1 + \mu) ^ {n - r - 1} \left(| \lambda_ {r + 1} | + \| N \| _ {F} / (1 + \mu)\right) ^ {k},
$$

$$
\left\| T _ {1 1} ^ {- k} \right\| _ {2} \leq (1 + \mu) ^ {r - 1} / \left(\left| \lambda_ {r} \right| - \left\| N \right\| _ {F} / (1 + \mu)\right) ^ {k},
$$

$$
\| V _ {k} - X W _ {k} \| _ {2} \leq \| V _ {k} \| _ {2} + \| X \| _ {2} \| W _ {k} \| _ {2} \leq 1 + \| T _ {1 2} \| _ {F} / \mathsf {s e p} (T _ {1 1}, T _ {2 2}).
$$

The bounds for $\| \ : T _ { 2 2 } ^ { k } \ : \| _ { 2 }$ and $\Vert \ T _ { 1 1 } ^ { - k } \ \Vert _ { 2 }$ follow from Lemma 7.3.2.

# Problems

# P7.3.1 Verify Equation (7.3.5).

P7.3.2 Suppose the eigenvalues of $A \in \mathbb { R } ^ { n \times n }$ satisfy $| \lambda _ { 1 } | = | \lambda _ { 2 } | > | \lambda _ { 3 } | \geq \cdot \cdot \cdot \geq | \lambda _ { n } |$ and that $\lambda _ { 1 }$ and λ2 are complex conjugates of one another. Let $S = \mathsf { s p a n } \{ y , z \}$ where $y , z \in \mathbb { R } ^ { n }$ satisfy $A ( y + i z ) =$ $\lambda _ { 1 } ( y + i z )$ . Show how the power method with a real starting vector can be used to compute an approximate basis for S.

P7.3.3 Assume $A \in \mathbb { R } ^ { n \times n }$ has eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { n }$ that satisfy

$$
\lambda = \lambda_ {1} = \lambda_ {2} = \lambda_ {3} = \lambda_ {4} > | \lambda_ {5} | \geq \dots \geq | \lambda_ {n} |
$$

where λ is positive. Assume that A has two Jordan blocks of the form.

$$
\left[ \begin{array}{c c} \lambda & 1 \\ 0 & \lambda \end{array} \right].
$$

Discuss the convergence properties of the power method when applied to this matrix and how the convergence might be accelerated.

P7.3.4 A matrix A is a positive matrix if $a _ { i j } ~ > ~ 0$ for all i and j. A vector $v \in \mathbb { R } ^ { n }$ is a positive vector if $v _ { i } > 0$ for all i. Perron’s theorem states that if A is a positive square matrix, then it has a unique dominant eigenvalue equal to its spectral radius $\rho ( A )$ and there is a positive vector x so that $A x = \rho ( A ) \cdot x$ . In this context, x is called the Perron vector and $\rho ( A )$ is called the Perron root. Assume that $A \in \mathbb { R } ^ { n \times n }$ is positive and $q \in \mathbb { R } ^ { n }$ is positive with unit 2-norm. Consider the following implementation of the power method (7.3.3):

$$
\begin{array}{l} z = A q, \lambda = q ^ {T} z \\ q = z, q = q / \| q \| _ {2}, z = A q, \lambda = q ^ {T} z \\ \end{array}
$$

(a) Adjust the termination criteria to guarantee (in principle) that the final λ and q satisfy $\tilde { A } q = \lambda q$ , where $\begin{array} { r } { \parallel \tilde { A } - A \parallel _ { 2 } \leq \delta } \end{array}$ and A˜ is positive. (b) Applied to a positive matrix $A \in \mathbb { R } ^ { n \times n }$ , the Collatz-Wielandt formula states that $\rho ( A )$ is the maximum value of the function f defined by

$$
f (x) = \min _ {1 \leq i \leq n} \frac {y _ {i}}{x _ {i}}
$$

where $\boldsymbol { x } \in \mathbb { R } ^ { n }$ is positive and $y = A x$ . Does it follow that $f ( A q ) \geq f ( q ) ?$ In other words, do the iterates $\{ q ^ { ( k ) } \}$ in the power method have the property that $f ( q ^ { ( k ) } )$ increases monotonically to the Perron root, assuming that $q ^ { ( 0 ) }$ is positive?

P7.3.5 (Read the previous problem for background.) A matrix A is a nonnegative matrix if $a _ { i j } \geq 0$ for all i and j. A matrix $A \in \mathbb { R } ^ { n \times n }$ is reducible if there is a permutation P so that $P ^ { T } A P$ is block triangular with two or more square diagonal blocks. A matrix that is not reducible is irreducible. The Perron-Frobenius theorem states that if A is a square, nonnegative, and irreducible, then $\rho ( A )$ , the Perron root, is an eigenvalue for A and there is a positive vector x, the Perron vector, so that $A x = \rho ( A ) { \cdot } x$ . Assume that $A _ { 1 } , A _ { 2 } , A _ { 3 } \in \mathbb { R } ^ { n \times n }$ are each positive and let the nonnegative matrix A be defined by

$$
A = \left[ \begin{array}{c c c} 0 & A _ {1} & 0 \\ 0 & 0 & A _ {2} \\ A _ {3} & 0 & 0 \end{array} \right].
$$

(a) Show that A is irreducible. (b) Let $B = A _ { 1 } A _ { 2 } A _ { 3 }$ . Show how to compute the Perron root and vector for A from the Perron root and vector for B. (c) Show that A has other eigenvalues with absolute value equal to the Perron root. How could those eigenvalues and the associated eigenvectors be computed?

P7.3.6 (Read the previous two problems for background.) A nonnegative matrix $P \in \mathbb { R } ^ { n \times n }$ is stochastic if the entries in each column sum to 1. A vector $v \in \mathbb { R } ^ { n }$ is a probability vector if its entries are nonnegative and sum to 1. (a) Show that if $P \in \mathbb { R } ^ { n \times n }$ is stochastic and $v \in \mathbb { R } ^ { n }$ is a probability vector, then $w = P v$ is also a probability vector. (b) The entries in a stochastic matrix $P \in \mathbb { R } ^ { n \times n }$ can be regarded as the transition probabilities associated with an n-state Markov Chain. Let $v _ { j }$ be the probability of being in state $j$ at time $t = t _ { \mathrm { c u r r e n t } }$ . In the Markov model, the probability of being in state i at time $t = t _ { \mathrm { n e x t } }$ is given by

$$
w _ {i} = \sum_ {j = 1} ^ {n} p _ {i j} v _ {j} \quad i = 1: n,
$$

$\mathrm { i } . \mathrm { e } . , w = P v$ . With the help of a biased coin, a surfer on the World Wide Web randomly jumps from page to page. Assume that the surfer is currently viewing web page $j$ and that the coin comes up heads with probability α. Here is how the surfer determines the next page to visit:

Step 1. A coin is tossed.

Step 2. If it comes up heads and web page j has at least one outlink, then the next page to visit is randomly selected from the list of outlink pages.

Step 3. Otherwise, the next page to visit is randomly selected from the list of all possible pages.

Let $P \in \mathbb { R } ^ { n \times n }$ be the matrix of transition probabilities that define this random process. Specify P in terms of $\alpha ,$ , the vector of ones $e ,$ and the link matrix $H \in \mathbb { R } ^ { n \times n }$ defined by

$$
h _ {i j} = \left\{ \begin{array}{l l} 1 & \text { if   there   is   a   link   on   web   page   } j \text {   to   web   page   } i \\ 0 & \text { otherwise } \end{array} \right.
$$

Hints: The number of nonzero components in $H ( : , j )$ is the number of outlinks on web page $j , P$ is a convex combination of a very sparse sparse matrix and a very dense rank-1 matrix. (c) Detail how the power method can be used to determine a probability vector x so that $P x = x$ . Strive to get as much computation “outside the loop” as possible. Note that in the limit we can expect to find the random surfer viewing web page i with probability $x _ { i }$ . Thus, a case can be made that more important pages are associated with the larger components of x. This is the basis of Google PageRank. If

$$
x _ {i _ {1}} \geq x _ {i _ {2}} \geq \dots \geq x _ {i _ {n}}
$$

then web page $i _ { k }$ has page rank k.

P7.3.7 (a) Show that if $\ b X \in \mathbb { C } ^ { n \times n }$ is nonsingular, then

$$
\parallel A \parallel_ {X} = \parallel X ^ {- 1} A X \parallel_ {2}
$$

defines a matrix norm with the property that

$$
\| A B \| _ {X} \leq \| A \| _ {X} \| B \| _ {X}.
$$

(b) Show that for any $\epsilon > 0$ there exists a nonsingular $\ b X \in \mathbb { C } ^ { n \times n }$ such that

$$
\parallel A \parallel_ {X} = \parallel X ^ {- 1} A X \parallel_ {2} \leq \rho (A) + \epsilon
$$

where $\rho ( A )$ is A’s spectral radius. Conclude that there is a constant M such that

$$
\| A ^ {k} \| _ {2} \leq M (\rho (A) + \epsilon) ^ {k}
$$

for all non-negative integers k. (Hint: Set $X = Q \ \mathrm { d i a g } ( 1 , a , . . . , a ^ { n - 1 } )$ where $Q ^ { H } A Q = D + N$ is A’s Schur decomposition.)

P7.3.8 Verify that (7.3.14) calculates the matrices $T _ { k }$ defined by (7.3.13).

P7.3.9 Suppose $A \in \mathbb { C } ^ { n \times n }$ is nonsingular and that $Q _ { 0 } \in \mathbb { C } ^ { n \times p }$ has orthonormal columns. The following iteration is referred to as inverse orthogonal iteration.

$$
\begin{array}{l} \text { for } k = 1, 2, \dots \\ \text { Solve } A Z _ {k} = Q _ {k - 1} \text { for } Z _ {k} \in \mathbb {C} ^ {n \times p} \\ Z _ {k} = Q _ {k} R _ {k} \quad \text {(QR factorization)} \\ \end{array}
$$

end

Explain why this iteration can usually be used to compute the p smallest eigenvalues of A in absolute value. Note that to implement this iteration it is necessary to be able to solve linear systems that involve A. If $p = 1$ , the method is referred to as the inverse power method.

# Notes and References for §7.3

For an excellent overview of the QR iteration and related procedures, see Watkins (MEP), Stewart (MAE), and Kressner (NMSE). A detailed, practical discussion of the power method is given in Wilkinson (AEP, Chap. 10). Methods are discussed for accelerating the basic iteration, for calculating nondominant eigenvalues, and for handling complex conjugate eigenvalue pairs. The connections among the various power iterations are discussed in:

B.N. Parlett and W.G. Poole (1973). “A Geometric Theory for the QR, LU, and Power Iterations,” SIAM J. Numer. Anal. 10, 389–412.

The QR iteration was concurrently developed in:

J.G.F. Francis (1961). “The QR Transformation: A Unitary Analogue to the LR Transformation,” Comput. J. 4, 265–71, 332–334.

V.N. Kublanovskaya (1961). “On Some Algorithms for the Solution of the Complete Eigenvalue Problem,” USSR Comput. Math. Phys. 3, 637–657.

As can be deduced from the title of the first paper by Francis, the LR iteration predates the QR iteration. The former very fundamental algorithm was proposed by:

H. Rutishauser (1958). “Solution of Eigenvalue Problems with the LR Transformation,” Nat. Bur. Stand. Appl. Math. Ser. 49, 47–81.

More recent, related work includes:

B.N. Parlett (1995). “The New qd Algorithms,” Acta Numerica 5, 459–491.

C. Ferreira and B.N. Parlett (2009). “Convergence of the LR Algorithm for a One-Point Spectrum Tridiagonal Matrix,” Numer. Math. 113, 417–431.

Numerous papers on the convergence and behavior of the QR iteration have appeared, see:

J.H. Wilkinson (1965). “Convergence of the LR, QR, and Related Algorithms,” Comput. J. 8, 77–84.

B.N. Parlett (1965). “Convergence of the Q-R Algorithm,” Numer. Math. 7, 187–93. (Correction in Numer. Math. 10, 163–164.)

B.N. Parlett (1966). “Singular and Invariant Matrices Under the QR Algorithm,” Math. Comput. 20, 611–615.

B.N. Parlett (1968). “Global Convergence of the Basic QR Algorithm on Hessenberg Matrices,” Math. Comput. 22, 803–817.

D.S. Watkins (1982). “Understanding the QR Algorithm,” SIAM Review 24, 427–440.

T. Nanda (1985). “Differential Equations and the QR Algorithm,” SIAM J. Numer. Anal. 22, 310–321.

D.S. Watkins (1993). “Some Perspectives on the Eigenvalue Problem,” SIAM Review 35, 430–471.

D.S. Watkins (2008). “The QR Algorithm Revisited,” SIAM Review 50, 133–145.

D.S. Watkins (2011). “Francis’s Algorithm,” AMS Monthly 118, 387–403.

A block analog of the QR iteration is discussed in:

M. Robb\`e and M. Sadkane (2005). “Convergence Analysis of the Block Householder Block Diagonalization Algorithm,” BIT 45, 181–195.

The following references are concerned with various practical and theoretical aspects of simultaneous iteration:

H. Rutishauser (1970). “Simultaneous Iteration Method for Symmetric Matrices,” Numer. Math. 16, 205–223.

M. Clint and A. Jennings (1971). “A Simultaneous Iteration Method for the Unsymmetric Eigenvalue Problem,” J. Inst. Math. Applic. 8, 111-121.

G.W. Stewart (1976). “Simultaneous Iteration for Computing Invariant Subspaces of Non-Hermitian Matrices,” Numer. Math. 25, 123–136.

A. Jennings (1977). Matrix Computation for Engineers and Scientists, John Wiley and Sons, New York.

Z. Bai and G.W. Stewart (1997). “Algorithm 776: SRRIT: a Fortran Subroutine to Calculate the Dominant Invariant Subspace of a Nonsymmetric Matrix,” ACM Trans. Math. Softw. 23, 494– 513.


---

<!-- golub_400_449 -->

Problems P7.3.4–P7.3.6 explore the relevance of the power method to the problem of computing the Perron root and vector of a nonnegative matrix. For further background and insight, see:

A. Berman and R.J. Plemmons (1994). Nonnegative Matrices in the Mathematical Sciences, SIAM Publications,Philadelphia, PA.

A.N. Langville and C.D. Meyer (2006). Google’s PageRank and Beyond, Princeton University Press, Princeton and Oxford. .

The latter volume is outstanding in how it connects the tools of numerical linear algebra to the design and analysis of Web browsers. See also:

W.J. Stewart (1994). Introduction to the Numerical Solution of Markov Chains, Princeton University Press, Princeton, NJ.

M.W. Berry, Z. Drmaˇc, and E.R. Jessup (1999). “Matrices, Vector Spaces, and Information Retrieval,” SIAM Review 41, 335–362.

A.N. Langville and C.D. Meyer (2005). “A Survey of Eigenvector Methods for Web Information Retrieval,” SIAM Review 47, 135–161.

A.N. Langville and C.D. Meyer (2006). “A Reordering for the PageRank Problem”, SIAM J. Sci. Comput. 27, 2112–2120.

A.N. Langville and C.D. Meyer (2006). “Updating Markov Chains with an Eye on Google’s PageRank,” SIAM J. Matrix Anal. Applic. 27, 968–987.

# 7.4 The Hessenberg and Real Schur Forms

In this and the next section we show how to make the QR iteration (7.3.1) a fast, effective method for computing Schur decompositions. Because the majority of eigenvalue/invariant subspace problems involve real data, we concentrate on developing the real analogue of (7.3.1) which we write as follows:

$$
H _ {0} = U _ {0} ^ {T} A U _ {0}
$$

for k = 1, 2, . . .

$$
H _ {k - 1} = U _ {k} R _ {k} \quad \text {(QR factorization)} \tag {7.4.1}
$$

$$
H _ {k} = R _ {k} U _ {k}
$$

end

Here, $A \in \mathbb { R } ^ { n \times n }$ , each $U _ { k } \in \mathbb { R } ^ { n \times n }$ is orthogonal, and each $R _ { k } \in \mathbb { R } ^ { n \times n }$ is upper triangular. A difficulty associated with this real iteration is that the $H _ { k }$ can never converge to triangular form in the event that A has complex eigenvalues. For this reason, we must lower our expectations and be content with the calculation of an alternative decomposition known as the real Schur decomposition.

In order to compute the real Schur decomposition efficiently we must carefully choose the initial orthogonal similarity transformation $U _ { 0 }$ in (7.4.1). In particular, if we choose $U _ { 0 }$ so that $H _ { 0 }$ is upper Hessenberg, then the amount of work per iteration is reduced from $O ( n ^ { 3 } )$ to $O ( n ^ { 2 } )$ . The initial reduction to Hessenberg form (the $U _ { 0 }$ computation) is a very important computation in its own right and can be realized by a sequence of Householder matrix operations.

# 7.4.1 The Real Schur Decomposition

A block upper triangular matrix with either 1-by-1 or 2-by-2 diagonal blocks is upper quasi-triangular. The real Schur decomposition amounts to a real reduction to upper quasi-triangular form.

Theorem 7.4.1 (Real Schur Decomposition). If $A \in \mathbb { R } ^ { n \times n }$ , then there exists an orthogonal $Q \in \mathbb { R } ^ { n \times n }$ such that

$$
Q ^ {T} A Q = \left[ \begin{array}{c c c c} R _ {1 1} & R _ {1 2} & \dots & R _ {1 m} \\ 0 & R _ {2 2} & \dots & R _ {2 m} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & R _ {m m} \end{array} \right] \tag {7.4.2}
$$

where each $R _ { i i }$ is either a 1-by-1 matrix or a ${ \it 2 - b y - 2 }$ matrix having complex conjugate eigenvalues.

Proof. The complex eigenvalues of A occur in conjugate pairs since the characteristic polynomial det $( z I - A )$ has real coefficients. Let k be the number of complex conjugate pairs in $\lambda ( A )$ . We prove the theorem by induction on k. Observe first that Lemma 7.1.2 and Theorem 7.1.3 have obvious real analogs. Thus, the theorem holds if $k = 0$ . Now suppose that $k \geq 1$ . If $\lambda = \gamma + i \mu \in \lambda ( A )$ and $\mu \neq 0$ , then there exist vectors y and z in $\mathbb { R } ^ { n } ( z \neq 0 )$ such that $A ( y + i z ) ~ = ~ ( \gamma + i \mu ) ( y + i z )$ , i.e.,

$$
A \left[ \begin{array}{c c} y & z \end{array} \right] = \left[ \begin{array}{c c} y & z \end{array} \right] \left[ \begin{array}{c c} \gamma & \mu \\ - \mu & \gamma \end{array} \right].
$$

The assumption that $\mu \neq 0$ implies that y and z span a 2-dimensional, real invariant subspace for A. It then follows from Lemma 7.1.2 that an orthogonal $U \in \mathbb { R } ^ { n \times n }$ exists such that

$$
U ^ {T} A U = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {n - 2} ^ {2}
$$

where $\lambda ( T _ { 1 1 } ) = \{ \lambda , \bar { \lambda } \}$ . By induction, there exists an orthogonal $\tilde { U }$ so $\tilde { U } ^ { T } T _ { 2 2 } \tilde { U }$ has the required structure. The theorem follows by setting $Q = U \cdot \mathrm { d i a g } ( I _ { 2 } , \tilde { U } )$ .

The theorem shows that any real matrix is orthogonally similar to an upper quasitriangular matrix. It is clear that the real and imaginary parts of the complex eigenvalues can be easily obtained from the 2-by-2 diagonal blocks. Thus, it can be said that the real Schur decomposition is an eigenvalue-revealing decomposition.

# 7.4.2 A Hessenberg QR Step

We now turn our attention to the efficient execution of a single QR step in (7.4.1). In this regard, the most glaring shortcoming associated with (7.4.1) is that each step requires a full QR factorization costing $O ( n ^ { \bar { 3 } } )$ flops. Fortunately, the amount of work per iteration can be reduced by an order of magnitude if the orthogonal matrix $U _ { 0 }$ is judiciously chosen. In particular, if $U _ { 0 } ^ { T } A U _ { 0 } = H _ { 0 } = ( h _ { i j } )$ is upper Hessenberg $( h _ { i j } = 0 ,$ , $i > j + 1 \}$ ), then each subsequent $H _ { k }$ requires only $O ( n ^ { 2 } )$ flops to calculate. To see this we look at the computations $H = Q R$ and $H _ { + } = R Q$ when H is upper Hessenberg. As described in §5.2.5, we can upper triangularize H with a sequence of $n - 1$ Givens rotations: $Q ^ { T } H \equiv G _ { n - 1 } ^ { T } \cdot \cdot \cdot G _ { 1 } ^ { T } \bar { H = } R$ . Here, $G _ { i } = G ( i , i + 1 , \theta _ { i } )$ . For the $n = 4$ case there are three Givens premultiplications:

$$
\left[ \begin{array}{c c c c} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right].
$$

See Algorithm 5.2.5. The computation $R Q = R ( G _ { 1 } \cdot \cdot \cdot G _ { n - 1 } )$ is equally easy to implement. In the $n = 4$ case there are three Givens post-multiplications:

$$
\left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right] \to \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \end{array} \right].
$$

Overall we obtain the following algorithm:

Algorithm 7.4.1 If H is an n-by-n upper Hessenberg matrix, then this algorithm overwrites H with $H _ { + } = R Q$ where $H = Q R$ is the QR factorization of H.

for $k = 1 { : } n - 1$

$$
[ c _ {k}, s _ {k} ] = \operatorname{givens} (H (k, k), H (k + 1, k))
$$

$$
H (k: k + 1, k: n) = \left[ \begin{array}{c c} c _ {k} & s _ {k} \\ - s _ {k} & c _ {k} \end{array} \right] ^ {T} H (k: k + 1, k: n)
$$

end

for $k = 1 { : } n - 1$

$$
H (1: k + 1, k: k + 1) = H (1: k + 1, k: k + 1) \left[ \begin{array}{c c} c _ {k} & s _ {k} \\ - s _ {k} & c _ {k} \end{array} \right]
$$

end

Let $G _ { k } = G ( k , k { + } 1 , \theta _ { k } )$ be the kth Givens rotation. It is easy to confirm that the matrix $Q = G _ { 1 } \cdot \cdot \cdot G _ { n - 1 }$ is upper Hessenberg. Thus, $R Q = H _ { + }$ is also upper Hessenberg. The algorithm requires about $6 n ^ { 2 }$ flops, an order of magnitude more efficient than a full matrix QR step (7.3.1).

# 7.4.3 The Hessenberg Reduction

It remains for us to show how the Hessenberg decomposition

$$
U _ {0} ^ {T} A U _ {0} = H, \quad U _ {0} ^ {T} U _ {0} = I \tag {7.4.3}
$$

can be computed. The transformation $U _ { 0 }$ can be computed as a product of Householder matrices $P _ { 1 } , \ldots , P _ { n - 2 }$ . The role of $P _ { k }$ is to zero the kth column below the subdiagonal. In the $n = 6$ case, we have

$$
\left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \end{array} \right] \xrightarrow {P _ {1}} \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \end{array} \right] \xrightarrow {P _ {2}}
$$

$$
\left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\end{array}\right] \stackrel {{P _ {3}}} {{\rightarrow}} \left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&0&\times&\times&\times\\0&0&0&\times&\times&\times\end{array}\right] \stackrel {{P _ {4}}} {{\rightarrow}} \left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&0&\times&\times&\times\\0&0&0&0&\times&\times\end{array}\right].
$$

In general, after $k - 1$ steps we have computed $k - 1$ Householder matrices $P _ { 1 } , \ldots , P _ { k - 1 }$ such that

$$
(P _ {1} \dots P _ {k - 1}) ^ {T} A (P _ {1} \dots P _ {k - 1}) = \left[ \begin{array}{c c c} B _ {1 1} & B _ {1 2} & B _ {1 3} \\ B _ {2 1} & B _ {2 2} & B _ {2 3} \\ 0 & B _ {3 2} & B _ {3 3} \end{array} \right] \begin{array}{c} k - 1 \\ 1 \\ n - k \end{array}
$$

is upper Hessenberg through its first $k - 1$ columns. Suppose $\tilde { P } _ { k }$ is an order- $( n - k )$ Householder matrix such that $\tilde { P } _ { k } B _ { 3 2 }$ is a multiple of $e _ { 1 } ^ { ( n - k ) }$ . If $P _ { k } = \mathrm { d i a g } ( I _ { k } , \tilde { P } _ { k } )$ , then

$$
(P _ {1} \dots P _ {k}) ^ {T} A (P _ {1} \dots P _ {k}) = \left[ \begin{array}{c c c} B _ {1 1} & B _ {1 2} & B _ {1 3} \tilde {P} _ {k} \\ B _ {2 1} & B _ {2 2} & B _ {2 3} \tilde {P} _ {k} \\ 0 & \tilde {P} _ {k} B _ {3 2} & \tilde {P} _ {k} B _ {3 3} \tilde {P} _ {k} \end{array} \right]
$$

is upper Hessenberg through its first k columns. Repeating this for $k = 1 { : } n - 2$ we obtain

Algorithm 7.4.2 (Householder Reduction to Hessenberg Form) Given $A \in \mathbb { R } ^ { n \times n }$ , the following algorithm overwrites A with $H = U _ { 0 } ^ { T } A U _ { 0 }$ where H is upper Hessenberg and $U _ { 0 }$ is a product of Householder matrices.

$$
\begin{array}{l} \text { for } k = 1: n - 2 \\ [ v, \beta ] = \text { house } (A (k + 1: n, k)) \\ A (k + 1: n, k: n) = (I - \beta v v ^ {T}) A (k + 1: n, k: n) \\ A (1: n, k + 1: n) = A (1: n, k + 1: n) (I - \beta v v ^ {T}) \\ \end{array}
$$

end

This algorithm requires $1 0 n ^ { 3 } / 3$ flops. If $U _ { 0 }$ is explicitly formed, an additional $4 n ^ { 3 } / 3$ flops are required. The kth Householder matrix can be represented in $A ( k + 2 { : } n , k )$ . See Martin and Wilkinson (1968) for a detailed description.

The roundoff properties of this method for reducing A to Hessenberg form are very desirable. Wilkinson (AEP, p. 351) states that the computed Hessenberg matrix $\hat { H }$ satisfies

$$
\hat {H} = Q ^ {T} (A + E) Q,
$$

where $Q$ is orthogonal and $\| E \| _ { F } \leq c n ^ { 2 } \mathbf { u } \| A \| _ { F }$ with c a small constant.

# 7.4.4 Level-3 Aspects

The Hessenberg reduction (Algorithm 7.4.2) is rich in level-2 operations: half gaxpys and half outer product updates. We briefly mention two ideas for introducing level-3 computations into the process.

The first involves a block reduction to block Hessenberg form and is quite straightforward. Suppose (for clarity) that $n = r N$ and write

$$
A   =   \left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right] _ {n - r} ^ {r}    .
$$

Suppose that we have computed the QR factorization $A _ { 2 1 } = \tilde { Q } _ { 1 } R _ { 1 }$ and that $\tilde { Q } _ { 1 }$ is in WY form. That is, we have $W _ { 1 } , Y _ { 1 } \in \dot { \mathbb { R } } ^ { ( n - r ) \times r }$ such that $\tilde { Q } _ { 1 } = I + W _ { 1 } Y _ { 1 } ^ { T }$ . (See §5.2.2 for details.) If $Q _ { 1 } = \mathrm { d i a g } ( I _ { r } , \tilde { Q } _ { 1 } )$ then

$$
Q _ {1} ^ {T} A Q _ {1} = \left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \tilde {Q} _ {1} \\ R _ {1} & \tilde {Q} _ {1} ^ {T} A _ {2 2} \tilde {Q} _ {1} \end{array} \right].
$$

Notice that the updates of the (1,2) and (2,2) blocks are rich in level-3 operations given that $\tilde { Q } _ { 1 }$ is in WY form. This fully illustrates the overall process as $Q _ { 1 } ^ { T } A Q _ { 1 }$ is block upper Hessenberg through its first block column. We next repeat the computations on the first r columns of $\tilde { Q } _ { 1 } ^ { T } A _ { 2 2 } \tilde { Q } _ { 1 }$ . After $N - 1$ such steps we obtain

$$
H = U _ {0} ^ {T} A U _ {0} = \left[ \begin{array}{c c c c c} H _ {1 1} & H _ {1 2} & \dots & \dots & H _ {1 N} \\ H _ {2 1} & H _ {2 2} & \dots & \dots & H _ {2 N} \\ 0 & \ddots & \ddots & \dots & \vdots \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & \dots & H _ {N, N - 1} & H _ {N N} \end{array} \right]
$$

where each $H _ { i j }$ is $r { \mathrm { - } } \mathrm { b y } { \mathrm { - } } r$ and $U _ { 0 } = Q _ { 1 } \cdot \cdot \cdot Q _ { N - 2 }$ with each $Q _ { i }$ in WY form. The overall algorithm has a level-3 fraction of the form $1 - O ( 1 / N )$ . Note that the subdiagonal blocks in H are upper triangular and so the matrix has lower bandwidth $r .$ . It is possible to reduce H to actual Hessenberg form by using Givens rotations to zero all but the first subdiagonal.

Dongarra, Hammarling, and Sorensen (1987) have shown how to proceed directly to Hessenberg form using a mixture of gaxpys and level-3 updates. Their idea involves minimal updating after each Householder transformation is generated. For example, suppose the first Householder $P _ { 1 }$ has been computed. To generate $P _ { 2 }$ we need just the second column of $P _ { 1 } A P _ { 1 }$ , not the full outer product update. To generate $P _ { 3 }$ we need just the thirrd column of $P _ { 2 } P _ { 1 } A P _ { 1 } P _ { 2 }$ , etc. In this way, the Householder matrices can be determined using only gaxpy operations. No outer product updates are involved. Once a suitable number of Householder matrices are known they can be aggregated and applied in level-3 fashion.

For more about the challenges of organizing a high-performance Hessenberg reduction, see Karlsson (2011).

# 7.4.5 Important Hessenberg Matrix Properties

The Hessenberg decomposition is not unique. If Z is any n-by-n orthogonal matrix and we apply Algorithm 7.4.2 to $Z ^ { T } A Z$ , then $Q ^ { T } A Q = H$ is upper Hessenberg where $Q \mathrm { ~ = ~ } Z U _ { 0 }$ . However, $Q e _ { 1 } = Z ( U _ { 0 } e _ { 1 } ) = Z e _ { 1 }$ suggesting that H is unique once the first column of Q is specified. This is essentially the case provided H has no zero subdiagonal entries. Hessenberg matrices with this property are said to be unreduced. Here is important theorem that clarifies these issues.

Theorem 7.4.2 ( Implicit Q Theorem ). Suppose $Q = { \left[ \begin{array} { l } { q _ { 1 } \left| \cdots \right| q _ { n } } \end{array} \right] }$ and $V =$ $[  v _ { 1 } | \cdots | v _ { n } ]$ are orthogonal matrices with the property that the matrices $Q ^ { T } A Q \ = H$ and $V ^ { T } A V = G$ are each upper Hessenberg where $A \in \mathbb { R } ^ { n \times n }$ . Let k denote the smallest positive integer for which $h _ { k + 1 , k } = 0 ,$ , with the convention that $k = n$ if H is unreduced. $H f q _ { 1 } = v _ { 1 }$ , then $q _ { i } = \pm v _ { i }$ and $| h _ { i , i - 1 } | = | g _ { i , i - 1 } | f o r i = 2 { : } k$ . Moreover, if $k < n$ , then $g _ { k + 1 , k } = 0$ .

Proof. Define the orthogonal matrix $W = \left[ \left. w _ { 1 } \right| \cdot \cdot \cdot \mid w _ { n } \right] = V ^ { T } Q$ and observe that $G W = W H$ . By comparing column i − 1 in this equation for $i = 2 { : } k$ we see that

$$
h _ {i, i - 1} w _ {i} = G w _ {i - 1} - \sum_ {j = 1} ^ {i - 1} h _ {j, i - 1} w _ {j}.
$$

Since $w _ { 1 } = e _ { 1 }$ , it follows that $\left[ \boldsymbol { w } _ { 1 } \left| \cdots \right| \boldsymbol { w } _ { k } \right]$ is upper triangular and so for $i = 2 { : } k$ we have $w _ { i } = \pm I _ { n } ( : , i ) = \pm e _ { i }$ . Since $w _ { i } = V ^ { T } q _ { i }$ and $h _ { i , i - 1 } = w _ { i } ^ { T } G w _ { i - 1 }$ it follows that $v _ { i } = \pm q _ { i }$ and

$$
| h _ {i, i - 1} | = | q _ {i} ^ {T} A q _ {i - 1} | = | v _ {i} ^ {T} A v _ {i - 1} | = | g _ {i, i - 1} |
$$

for i = 2:k. If $k < n$ , then

$$
g _ {k + 1, k} = e _ {k + 1} ^ {T} G e _ {k} = \pm e _ {k + 1} ^ {T} G W e _ {k} = \pm e _ {k + 1} ^ {T} W H e _ {k}
$$

$$
= \pm e _ {k + 1} ^ {T} \sum_ {i = 1} ^ {k} h _ {i k} W e _ {i} = \pm \sum_ {i = 1} ^ {k} h _ {i k} e _ {k + 1} ^ {T} e _ {i} = 0,
$$

completing the proof of the theorem.

The gist of the implicit Q theorem is that if $Q ^ { T } A Q = H$ and $Z ^ { T } A Z = G$ are each unreduced upper Hessenberg matrices and Q and Z have the same first column, then G and H are “essentially equal” in the sense that $G = D ^ { - 1 } H D$ where $D = \mathrm { d i a g } ( \pm 1 , \ldots , \pm 1 )$ .

Our next theorem involves a new type of matrix called a Krylov matrix. If $A \in \mathbb { R } ^ { n \times n }$ and $v \in \mathbb { R } ^ { n }$ , then the Krylov matrix $K ( A , v , j ) \in \mathbb { R } ^ { n \times j }$ is defined by

$$
K (A, v, j) = \left[ v \mid A v \mid \dots \mid A ^ {j - 1} v \right].
$$

It turns out that there is a connection between the Hessenberg reduction $Q ^ { T } A Q = H$ and the QR factorization of the Krylov matrix $K ( A , Q ( : , 1 ) , n )$ .

Theorem 7.4.3. Suppose $Q \in \mathbb { R } ^ { n \times n }$ is an orthogonal matrix and $A \in \mathbb { R } ^ { n \times n }$ . Then $Q ^ { T } A Q = H$ is an unreduced upper Hessenberg matrix if and only i ${ } ^ { \sharp } Q ^ { T } K ( A , Q ( : , 1 ) , n ) =$ R is nonsingular and upper triangular.

Proof. Suppose $Q \in \mathbb { R } ^ { n \times n }$ is orthogonal and set $H = Q ^ { T } A Q$ . Consider the identity

$$
Q ^ {T} K (A, Q (:, 1), n) = \left[ e _ {1} \mid H e _ {1} \mid \dots \mid H ^ {n - 1} e _ {1} \right] \equiv R.
$$

If H is an unreduced upper Hessenberg matrix, then it is clear that R is upper triangular with $r _ { i i } = h _ { 2 1 } h _ { 3 2 } \cdot \cdot \cdot h _ { i , i - 1 }$ for $i = 2 { : } n$ . Since $r _ { 1 1 } = 1$ it follows that R is nonsingular.

To prove the converse, suppose R is upper triangular and nonsingular. Since $R ( : , k + 1 ) = H R ( : , k )$ it follows that $H ( : , k ) \in { \mathsf { s p a n } } { \left\{ \begin{array} { l l } { e _ { 1 } , \ldots , e _ { k + 1 } } \end{array} \right\} }$ . This implies that H is upper Hessenberg. Since $r _ { n n } = h _ { 2 1 } h _ { 3 2 } \cdot \cdot \cdot h _ { n , n - 1 } \neq 0$ it follows that H is also unreduced.

Thus, there is more or less a correspondence between nonsingular Krylov matrices and orthogonal similarity reductions to unreduced Hessenberg form.

Our last result is about the geometric multiplicity of an eigenvalue of an unreduced upper Hessenberg matrix.

Theorem 7.4.4. If λ is an eigenvalue of an unreduced upper Hessenberg matrix $H \in \mathbb { R } ^ { n \times n }$ , then its geometric multiplicity is 1.

Proof. For any $\lambda \in \mathbb { C }$ we have rank $( A - \lambda I ) \geq n - 1$ because the first $n - 1$ columns of $H - \lambda I$ are independent.

# 7.4.6 Companion Matrix Form

Just as the Schur decomposition has a nonunitary analogue in the Jordan decomposition, so does the Hessenberg decomposition have a nonunitary analog in the companion matrix decomposition. Let $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and suppose that the Krylov matrix $K = K ( A , x , n )$ is nonsingular. If $c = c ( 0 { : } n - 1 )$ solves the linear system $K c = - A ^ { n } x$ , then it follows that $A K = K C$ where C has the form

$$
C = \left[ \begin{array}{c c c c c} 0 & 0 & \dots & 0 & - c _ {0} \\ 1 & 0 & \dots & 0 & - c _ {1} \\ 0 & 1 & \dots & 0 & - c _ {2} \\ \vdots & \vdots & \vdots & \vdots & \vdots \\ 0 & 0 & \dots & 1 & - c _ {n - 1} \end{array} \right]. \tag {7.4.4}
$$

The matrix C is said to be a companion matrix. Since

$$
\det (z I - C) = c _ {0} + c _ {1} z + \dots + c _ {n - 1} z ^ {n - 1} + z ^ {n},
$$

it follows that if K is nonsingular, then the decomposition $K ^ { - 1 } A K = C$ displays A’s characteristic polynomial. This, coupled with the sparseness of C, leads to “companion matrix methods” in various application areas. These techniques typically involve:

Step 1. Compute the Hessenberg decomposition $U _ { 0 } ^ { T } A U _ { 0 } = H$ .

Step 2. Hope H is unreduced and set $Y = \left[ e _ { 1 } \mid H e _ { 1 } \mid . . . \mid H ^ { n - 1 } e _ { 1 } \right]$ .

Step 3. Solve $Y C = H Y$ for C.

Unfortunately, this calculation can be highly unstable. A is similar to an unreduced Hessenberg matrix only if each eigenvalue has unit geometric multiplicity. Matrices that have this property are called nonderogatory. It follows that the matrix Y above can be very poorly conditioned if A is close to a derogatory matrix.

A full discussion of the dangers associated with companion matrix computation can be found in Wilkinson (AEP, pp. 405ff.).

# Problems

P7.4.1 Suppose $A \in \mathbb { R } ^ { n \times n }$ and $z \in \mathbb { R } ^ { n }$ . Give a detailed algorithm for computing an orthogonal Q such that $\bar { Q } ^ { \bar { T } } A Q$ is upper Hessenberg and $Q ^ { T } .$ z is a multiple of e1. Hint: Reduce z first and then apply Algorithm 7.4.2.

P7.4.2 Develop a similarity reduction to Hessenberg form using Gauss transforms with pivoting. How many flops are required. See Businger (1969).

P7.4.3 In some situations, it is necessary to solve the linear system $( A + z I ) x = b$ for many different values of $z \in \mathbb { R }$ and $b \in \mathbb { R } ^ { n }$ . Show how this problem can be efficiently and stably solved using the Hessenberg decomposition.

P7.4.4 Suppose $H \in \mathbb { R } ^ { n \times n }$ is an unreduced upper Hessenberg matrix. Show that there exists a diagonal matrix D such that each subdiagonal element of $D ^ { - 1 } H D$ is equal to 1. What is $\kappa _ { 2 } ( D ) ?$

P7.4.5 Suppose $W , Y \in \mathbb { R } ^ { n \times n }$ and define the matrices C and B by

$$
C = W + i Y, \qquad B = \left[ \begin{array}{c c} W & - Y \\ Y & W \end{array} \right].
$$

Show that if $\lambda \in \lambda ( C )$ is real, then $\lambda \in \lambda ( B )$ . Relate the corresponding eigenvectors.

P7.4.6 Suppose

$$
A = \left[ \begin{array}{c c} w & x \\ y & z \end{array} \right]
$$

is a real matrix having eigenvalues $\lambda \pm i \mu ,$ , where µ is nonzero. Give an algorithm that stably determines $c = \cos ( \theta )$ and $s = \sin ( \theta )$ such that

$$
\left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} \left[ \begin{array}{c c} w & x \\ y & z \end{array} \right] \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] = \left[ \begin{array}{c c} \lambda & \beta \\ \alpha & \lambda \end{array} \right]
$$

where $\alpha \beta = - \mu ^ { 2 }$ .

P7.4.7 Suppose $( \lambda , x )$ is a known eigenvalue-eigenvector pair for the upper Hessenberg matrix $H \in \mathbb { R } ^ { n \times n }$ . Give an algorithm for computing an orthogonal matrix P such that

$$
P ^ {T} H P = \left[ \begin{array}{l l} \lambda & w ^ {T} \\ 0 & H _ {1} \end{array} \right]
$$

where $H _ { 1 } \in \mathbb { R } ^ { ( n - 1 ) \times ( n - 1 ) }$ is upper Hessenberg. Compute P as a product of Givens rotations.

P7.4.8 Suppose $H \in \mathbb { R } ^ { n \times n }$ has lower bandwidth p. Show how to compute $Q \in \mathbb { R } ^ { n \times n } ;$ , a product of Givens rotations, such that $Q ^ { T } H Q$ is upper Hessenberg. How many flops are required?

P7.4.9 Show that if C is a companion matrix with distinct eigenvalues $\lambda _ { 1 } , \ldots , \lambda _ { n }$ , then $V C V ^ { - 1 } =$ $\mathrm { d i a g } ( \lambda _ { 1 } , \ldots , \lambda _ { n } )$ where

$$
V = \left[ \begin{array}{c c c c} 1 & \lambda_ {1} & \dots & \lambda_ {1} ^ {n - 1} \\ 1 & \lambda_ {2} & \dots & \lambda_ {2} ^ {n - 1} \\ \vdots & \vdots & \ddots & \vdots \\ 1 & \lambda_ {n} & \dots & \lambda_ {n} ^ {n - 1} \end{array} \right].
$$

# Notes and References for 7.4

The real Schur decomposition was originally presented in:

F.D. Murnaghan and A. Wintner (1931). “A Canonical Form for Real Matrices Under Orthogonal Transformations,” Proc. Nat. Acad. Sci. 17, 417–420.

A thorough treatment of the reduction to Hessenberg form is given in Wilkinson (AEP, Chap. 6), and Algol procedures appear in:   
R.S. Martin and J.H. Wilkinson (1968). “Similarity Reduction of a General Matrix to Hessenberg Form,” Numer. Math. 12, 349–368.   
Givens rotations can also be used to compute the Hessenberg decomposition, see:   
W. Rath (1982). “Fast Givens Rotations for Orthogonal Similarity,” Numer. Math. 40, 47–56.   
The high-performance computation of the Hessenberg reduction is a major challenge because it is a two-sided factorization, see:   
J.J. Dongarra, L. Kaufman, and S. Hammarling (1986). “Squeezing the Most Out of Eigenvalue Solvers on High Performance Computers,” Lin. Alg. Applic. 77, 113–136.   
J.J. Dongarra, S. Hammarling, and D.C. Sorensen (1989). “Block Reduction of Matrices to Condensed Forms for Eigenvalue Computations,” J. ACM 27, 215–227.   
M.W. Berry, J.J. Dongarra, and Y. Kim (1995). “A Parallel Algorithm for the Reduction of a Nonsymmetric Matrix to Block Upper Hessenberg Form,” Parallel Comput. 21, 1189–1211.   
G. Quintana-Orti and R. Van De Geijn (2006). “Improving the Performance of Reduction to Hessenberg Form,” ACM Trans. Math. Softw. 32, 180–194.   
S. Tomov, R. Nath, and J. Dongarra (2010). “Accelerating the Reduction to Upper Hessenberg, Tridiagonal, and Bidiagonal Forms Through Hybrid GPU-Based Computing,” Parallel Comput. 36, 645–654.   
L. Karlsson (2011). “Scheduling of Parallel Matrix Computations and Data Layout Conversion for HPC and Multicore Architectures,” PhD Thesis, University of Ume˚a.

Reaching the Hessenberg form via Gauss transforms is discussed in:

P. Businger (1969). “Reducing a Matrix to Hessenberg Form,” Math. Comput. 23, 819–821. G.W. Howell and N. Diaa (2005). “Algorithm 841: BHESS: Gaussian Reduction to a Similar Banded Hessenberg Form,” ACM Trans. Math. Softw. 31, 166–185.

Some interesting mathematical properties of the Hessenberg form may be found in:

B.N. Parlett (1967). “Canonical Decomposition of Hessenberg Matrices,” Math. Comput. 21, 223– 227.

Although the Hessenberg decomposition is largely appreciated as a “front end” decomposition for the QR iteration, it is increasingly popular as a cheap alternative to the more expensive Schur decomposition in certain problems. For a sampling of applications where it has proven to be very useful, consult:

W. Enright (1979). “On the Efficient and Reliable Numerical Solution of Large Linear Systems of O.D.E.’s,” IEEE Trans. Autom. Contr. AC-24, 905–908.   
G.H. Golub, S. Nash and C. Van Loan (1979). “A Hessenberg-Schur Method for the Problem AX + XB = C,” IEEE Trans. Autom. Contr. AC-24, 909–913.   
A. Laub (1981). “Efficient Multivariable Frequency Response Computations,” IEEE Trans. Autom. Contr. AC-26, 407–408.   
C.C. Paige (1981). “Properties of Numerical Algorithms Related to Computing Controllability,” IEEE Trans. Auto. Contr. AC-26, 130–138.   
G. Miminis and C.C. Paige (1982). “An Algorithm for Pole Assignment of Time Invariant Linear Systems,” Int. J. Contr. 35, 341–354.   
C. Van Loan (1982). “Using the Hessenberg Decomposition in Control Theory,” in Algorithms and Theory in Filtering and Control , D.C. Sorensen and R.J. Wets (eds.), Mathematical Programming Study No. 18, North Holland, Amsterdam, 102–111.   
C.D. Martin and C.F. Van Loan (2006). “Solving Real Linear Systems with the Complex Schur Decomposition,” SIAM J. Matrix Anal. Applic. 29, 177–183.

The advisability of posing polynomial root problems as companion matrix eigenvalue problem is discussed in:

A. Edelman and H. Murakami (1995). “Polynomial Roots from Companion Matrix Eigenvalues,” Math. Comput. 64, 763–776.

# 7.5 The Practical QR Algorithm

We return to the Hessenberg QR iteration, which we write as follows:

$$
H = U _ {0} ^ {T} A U _ {0} \quad \text {(Hessenberg reduction)}
$$

for k = 1, 2, . . .

$$
H = U R \quad (\text { QR   factorization }) \tag {7.5.1}
$$

$$
H = R U
$$

end

Our aim in this section is to describe how the H’s converge to upper quasi-triangular form and to show how the convergence rate can be accelerated by incorporating shifts.

# 7.5.1 Deflation

Without loss of generality we may assume that each Hessenberg matrix H in (7.5.1) is unreduced. If not, then at some stage we have

$$
H = \left[ \begin{array}{c c} H _ {1 1} & H _ {1 2} \\ 0 & H _ {2 2} \end{array} \right] _ {n - p} ^ {p}
$$

where $1 \leq p < n$ and the problem decouples into two smaller problems involving $H _ { 1 1 }$ and $H _ { 2 2 }$ . The term deflation is also used in this context, usually when $p = n - 1$ o r $n - 2$ .

In practice, decoupling occurs whenever a subdiagonal entry in H is suitably small. For example, if

$$
\left| h _ {p + 1, p} \right| \leq c \mathbf {u} \left(\left| h _ {p p} \right| + \left| h _ {p + 1, p + 1} \right|\right) \tag {7.5.2}
$$

for a small constant c, then $h _ { p + 1 , p }$ can justifiably be set to zero because rounding errors of order u H  are typically present throughout the matrix anyway.

# 7.5.2 The Shifted QR Iteration

Let $\mu \in \mathbb { R }$ and consider the iteration:

$$
H = U _ {0} ^ {T} A U _ {0} \quad \text {(Hessenberg reduction)}
$$

for k = 1, 2, . . .

Determine a scalar $\mu .$

$$
H - \mu I = U R \quad (\text { QR   factorization }) \tag {7.5.3}
$$

$$
H = R U + \mu I
$$

end

The scalar $\mu$ is referred to as a shift . Each matrix H generated in (7.5.3) is similar to A, since

$$
R U + \mu I = U ^ {T} (U R + \mu I) U = U ^ {T} H U.
$$

If we order the eigenvalues $\lambda _ { i }$ of A so that

$$
\left| \lambda_ {1} - \mu \right| \geq \dots \geq \left| \lambda_ {n} - \mu \right|,
$$

and $\mu$ is fixed from iteration to iteration, then the theory of §7.3 says that the pth subdiagonal entry in H converges to zero with rate

$$
\left| \frac {\lambda_ {p + 1} - \mu}{\lambda_ {p} - \mu} \right| ^ {k}.
$$

Of course, if $\lambda _ { p } = \lambda _ { p + 1 }$ , then there is no convergence at all. But if, for example, µ is much closer to $\lambda _ { n }$ than to the other eigenvalues, then the zeroing of the $( n , n - 1 )$ entry is rapid. In the extreme case we have the following:

Theorem 7.5.1. Let $\mu$ be an eigenvalue of an $n { - } b y { - } n$ unreduced Hessenberg matrix H . If

$$
\tilde {H} = R U + \mu I,
$$

where $H - \mu I = U R$ is the QR factorization of $H - \mu I$ , then $\tilde { h } _ { n , n - 1 } = 0$ and $\tilde { h } _ { n n } = \mu$ .

Proof. Since H is an unreduced Hessenberg matrix the first n − 1 columns of $H - \mu I$ are independent, regardless of µ. Thus, if $U R = \left( H - \mu I \right)$ is the QR factorization then $r _ { i i } \neq 0$ for $i = 1 { : } n - 1$ . But if $H - \mu I$ is singular, then $r _ { 1 1 } \cdot \cdot \cdot r _ { n n } = 0$ . Thus, $r _ { n n } = 0$ and $\tilde { H } ( n , : ) = [ 0 , . . . , 0 , \mu ]$ .

The theorem says that if we shift by an exact eigenvalue, then in exact arithmetic deflation occurs in one step.

# 7.5.3 The Single-Shift Strategy

Now let us consider varying µ from iteration to iteration incorporating new information about $\lambda ( A )$ as the subdiagonal entries converge to zero. A good heuristic is to regard $h _ { n n }$ as the best approximate eigenvalue along the diagonal. If we shift by this quantity during each iteration, we obtain the single-shift QR iteration:

for $k = 1 , 2 , \dots$

$$
\mu = H (n, n)
$$

$$
H - \mu I = U R \quad (\text { QR   factorization }) \tag {7.5.4}
$$

$$
H = R U + \mu I
$$

end

If the $( n , n - 1 )$ entry converges to zero, it is likely to do so at a quadratic rate. To see this, we borrow an example from Stewart (IMC, p. 366). Suppose H is an unreduced upper Hessenberg matrix of the form

$$
H = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \epsilon & h _ {n n} \end{array} \right]
$$

and that we perform one step of the single-shift QR algorithm, i.e.,

$$
U R = H - h _ {n n}
$$

$$
\tilde {H} = R U + h _ {n n} I.
$$

After $n - 2$ steps in the orthogonal reduction of $H - h _ { n n } I$ to upper triangular form we obtain a matrix with the following structure:

$$
H = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & a & b \\ 0 & 0 & 0 & \epsilon & 0 \end{array} \right].
$$

It is not hard to show that

$$
\tilde {h} _ {n, n - 1} = - \frac {\epsilon^ {2} b}{a ^ {2} + \epsilon^ {2}}.
$$

If we assume that $\epsilon \ll a$ , then it is clear that the new $( n , n - 1 )$ entry has order $\epsilon ^ { 2 }$ , precisely what we would expect of a quadratically converging algorithm.

# 7.5.4 The Double-Shift Strategy

Unfortunately, difficulties with (7.5.4) can be expected if at some stage the eigenvalues $a _ { 1 }$ and $a _ { 2 }$ of

$$
G = \left[ \begin{array}{c c} h _ {m m} & h _ {m n} \\ h _ {n m} & h _ {n n} \end{array} \right], \qquad m = n - 1, \tag {7.5.5}
$$

are complex for then $h _ { n n }$ would tend to be a poor approximate eigenvalue.

A way around this difficulty is to perform two single-shift QR steps in succession using $a _ { 1 }$ and $a _ { 2 }$ as shifts:

$$
H - a _ {1} I = U _ {1} R _ {1}
$$

$$
H _ {1} = R _ {1} U _ {1} + a _ {1} I \tag {7.5.6}
$$

$$
H _ {1} - a _ {2} I = U _ {2} R _ {2}
$$

$$
H _ {2} = R _ {2} U _ {2} + a _ {2} I
$$

These equations can be manipulated to show that

$$
(U _ {1} U _ {2}) (R _ {2} R _ {1}) = M \tag {7.5.7}
$$

where M is defined by

$$
M = (H - a _ {1} I) (H - a _ {2} I). \tag {7.5.8}
$$

Note that M is a real matrix even if G’s eigenvalues are complex since

$$
M = H ^ {2} - s H + t I
$$

where

$$
s = a _ {1} + a _ {2} = h _ {m m} + h _ {n n} = \operatorname{tr} (G) \in \mathbb {R}
$$

and

$$
t = a _ {1} a _ {2} = h _ {m m} h _ {n n} - h _ {m n} h _ {n m} = \det (G) \in \mathbb {R}.
$$

Thus, (7.5.7) is the QR factorization of a real matrix and we may choose $U _ { 1 }$ and $U _ { 2 }$ so that $Z = U _ { 1 } U _ { 2 }$ is real orthogonal. It then follows that

$$
H _ {2} = U _ {2} ^ {H} H _ {1} U _ {2} = U _ {2} ^ {H} (U _ {1} ^ {H} H U _ {1}) U _ {2} = (U _ {1} U _ {2}) ^ {H} H (U _ {1} U _ {2}) = Z ^ {T} H Z
$$

is real.

Unfortunately, roundoff error almost always prevents an exact return to the real field. A real $H _ { 2 }$ could be guaranteed if we

• explicitly form the real matrix $M = H ^ { 2 } - s H + t I$ ,   
• compute the real QR factorization $M = Z R .$ , and   
• set $H _ { 2 } = Z ^ { T } H Z .$ .

But since the first of these steps requires $O ( n ^ { 3 } )$ flops, this is not a practical course of action.

# 7.5.5 The Double-Implicit-Shift Strategy

Fortunately, it turns out that we can implement the double-shift step with $O ( n ^ { 2 } )$ flops by appealing to the implicit Q theorem of §7.4.5. In particular we can effect the transition from H to $H _ { 2 }$ in $O ( n ^ { 2 } )$ flops if we

• compute $M e _ { 1 }$ , the first column of $M$   
• determine a Householder matrix $P _ { 0 }$ such that $P _ { 0 } ( M e _ { 1 } )$ is a multiple of $e _ { 1 }$   
• compute Householder matrices $P _ { 1 } , \ldots , P _ { n - 2 }$ such that if

$$
Z _ {1} = P _ {0} P _ {1} \dots P _ {n - 2},
$$

then $Z _ { 1 } ^ { T } H Z _ { 1 }$ is upper Hessenberg and the first columns of Z and $Z _ { 1 }$ are the same.

Under these circumstances, the implicit Q theorem permits us to conclude that, if $Z ^ { T } H Z$ and $Z _ { 1 } ^ { T } H Z _ { 1 }$ are both unreduced upper Hessenberg matrices, then they are essentially equal. Note that if these Hessenberg matrices are not unreduced, then we can effect a decoupling and proceed with smaller unreduced subproblems.

Let us work out the details. Observe first that $P _ { 0 }$ can be determined in $O ( 1 )$ flops since $M e _ { 1 } = [ x , y , z , 0 , \ldots , 0 ] ^ { T }$ where

$$
\begin{array}{l} x = h _ {1 1} ^ {2} + h _ {1 2} h _ {2 1} - s h _ {1 1} + t, \\ y = h _ {2 1} (h _ {1 1} + h _ {2 2} - s), \\ z = h _ {2 1} h _ {3 2}. \\ \end{array}
$$

Since a similarity transformation with $P _ { 0 }$ only changes rows and columns 1, 2, and 3, we see that

$$
P _ {0} H P _ {0} = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right].
$$

Now the mission of the Householder matrices $P _ { 1 } , \ldots , P _ { n - 2 }$ is to restore this matrix to upper Hessenberg form. The calculation proceeds as follows:

$$
\left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right] \xrightarrow {P _ {1}} \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right] \xrightarrow {P _ {2}}
$$

$$
\left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\end{array}\right] \overset {P _ {3}} {\rightarrow} \left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&0&\times&\times&\times\\0&0&0&\times&\times&\times\end{array}\right] \overset {P _ {4}} {\rightarrow} \left[\begin{array}{c c c c c c}\times&\times&\times&\times&\times&\times\\\times&\times&\times&\times&\times&\times\\0&\times&\times&\times&\times&\times\\0&0&\times&\times&\times&\times\\0&0&0&\times&\times&\times\\0&0&0&0&\times&\times\end{array}\right].
$$

Each $P _ { k }$ is the identity with a 3-by-3 or 2-by-2 Householder somewhere along its diagonal, e.g.,

$$
\begin{array}{l} P _ {1} = \left[ \begin{array}{c c c c c c} 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 \end{array} \right], \quad P _ {2} = \left[ \begin{array}{c c c c c c} 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & \times & \times & \times & 0 \\ 0 & 0 & \times & \times & \times & 0 \\ 0 & 0 & \times & \times & \times & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 \end{array} \right], \\ P _ {3} = \left[ \begin{array}{c c c c c c} 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \end{array} \right], \quad P _ {4} = \left[ \begin{array}{c c c c c c} 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right]. \\ \end{array}
$$

The applicability of Theorem 7.4.3 (the implicit Q theorem) follows from the observation that $P _ { k } e _ { 1 } = e _ { 1 }$ for $k = 1 { : } n - 2$ and that $P _ { 0 }$ and $Z$ have the same first column. Hence, $Z _ { 1 } e _ { 1 } = Z e _ { 1 }$ , and we can assert that $Z _ { 1 }$ essentially equals $Z$ provided that the upper Hessenberg matrices $Z ^ { T } H Z$ and $Z _ { 1 } ^ { T } H Z _ { 1 }$ are each unreduced.

The implicit determination of $H _ { 2 }$ from H outlined above was first described by Francis (1961) and we refer to it as a Francis QR step. The complete Francis step is summarized as follows:

Algorithm 7.5.1 (Francis QR step) Given the unreduced upper Hessenberg matrix $H \in \mathbb { R } ^ { n \times n }$ whose trailing 2-by-2 principal submatrix has eigenvalues $a _ { 1 }$ and $a _ { 2 }$ , this algorithm overwrites H with $Z ^ { T } H Z$ , where $Z$ is a product of Householder matrices and $Z ^ { T } ( H - a _ { 1 } I ) ( H - a _ { 2 } I )$ is upper triangular.

$m = n - 1$ {Compute first column of $(H - a_{1}I)(H - a_{2}I)\}$ $s = H(m,m) + H(n,n)$ $t = H(m,m)\cdot H(n,n) - H(m,n)\cdot H(n,m)$ $x = H(1,1)\cdot H(1,1) + H(1,2)\cdot H(2,1) - s\cdot H(1,1) + t$ $y = H(2,1)\cdot (H(1,1) + H(2,2) - s)$ $z = H(2,1)\cdot H(3,2)$ for $k = 0:n - 3$ $[v,\beta ] = \mathsf{house}([xyz]^T)$ $q = \max \{1,k\}$ $H(k + 1:k + 3,q:n) = (I - \beta vv^T)\cdot H(k + 1:k + 3,q:n)$ $r = \min \{k + 4,n\}$ $H(1:r,k + 1:k + 3) = H(1:r,k + 1:k + 3)\cdot (I - \beta vv^T)$ $x = H(k + 2,k + 1)$ $y = H(k + 3,k + 1)$ if $k <   n - 3$ $z = H(k + 4,k + 1)$ end

end

$[v, \beta] = \text{house}([x y]^T)$ $H(n - 1:n, n - 2:n) = (I - \beta vv^T) \cdot H(n - 1:n, n - 2:n)$ $H(1:n, n - 1:n) = H(1:n, n - 1:n) \cdot (I - \beta vv^T)$

This algorithm requires $1 0 n ^ { 2 }$ flops. If Z is accumulated into a given orthogonal matrix, an additional $1 0 n ^ { 2 }$ flops are necessary.

# 7.5.6 The Overall Process

Reduction of A to Hessenberg form using Algorithm 7.4.2 and then iteration with Algorithm 7.5.1 to produce the real Schur form is the standard means by which the dense unsymmetric eigenproblem is solved. During the iteration it is necessary to monitor the subdiagonal elements in H in order to spot any possible decoupling. How this is done is illustrated in the following algorithm:

Algorithm 7.5.2 (QR Algorithm) Given $A \in \mathbb { R } ^ { n \times n }$ and a tolerance tol greater than the unit roundoff, this algorithm computes the real Schur canonical form ${ \overset { \cdot } { Q } } ^ { T } A Q = T$ . If $Q$ and $T$ are desired, then $T$ is stored in H. If only the eigenvalues are desired, then diagonal blocks in $T$ are stored in the corresponding positions in H.

Use Algorithm 7.4.2 to compute the Hessenberg reduction

$$
H = U _ {0} ^ {T} A U _ {0} \text { where } U _ {0} = P _ {1} \dots P _ {n - 2}.
$$

If $Q$ is desired form $Q = P _ { 1 } \cdots P _ { n - 2 }$ . (See §5.1.6.)

$$
\text { until } q = n
$$

Set to zero all subdiagonal elements that satisfy:

$$
\left| h _ {i, i - 1} \right| \leq \operatorname{tol} \cdot \left(\left| h _ {i i} \right| + \left| h _ {i - 1, i - 1} \right|\right).
$$

Find the largest nonnegative q and the smallest non-negative p such that

$$
H = \left[ \begin{array}{c c c} H _ {1 1} & H _ {1 2} & H _ {1 3} \\ 0 & H _ {2 2} & H _ {2 3} \\ 0 & 0 & H _ {3 3} \end{array} \right] \begin{array}{c} p \\ n - p - q \\ q \end{array}
$$

where $H _ { 3 3 }$ is upper quasi-triangular and $H _ { 2 2 }$ is unreduced.

$$
\text { if } q <   n
$$

Perform a Francis QR step on $H _ { 2 2 } \colon H _ { 2 2 } = Z ^ { T } H _ { 2 2 } Z .$

if Q is required

$$
Q = Q \cdot \operatorname{diag} \left(I _ {p}, Z, I _ {q}\right)
$$

$$
H _ {1 2} = H _ {1 2} Z
$$

$$
H _ {2 3} = Z ^ {T} H _ {2 3}
$$

end

end

end

Upper triangularize all 2-by-2 diagonal blocks in H that have real eigenvalues and accumulate the transformations (if necessary).

This algorithm requires $2 5 n ^ { 3 }$ flops if $Q$ and $T$ are computed. If only the eigenvalues are desired, then $1 0 n ^ { 3 }$ flops are necessary. These flops counts are very approximate and are based on the empirical observation that on average only two Francis iterations are required before the lower 1-by-1 or 2-by-2 decouples.

The roundoff properties of the QR algorithm are what one would expect of any orthogonal matrix technique. The computed real Schur form $\hat { T }$ is orthogonally similar to a matrix near to A, i.e.,

$$
Q ^ {T} (A + E) Q = \hat {T}
$$

where $Q ^ { T } Q = I$ and $\parallel E \parallel _ { 2 } \approx \mathbf { u } \parallel A \parallel _ { 2 }$ . The computed $\hat { Q }$ is almost orthogonal in the sense that ${ \hat { Q } } ^ { T } { \hat { Q } } = I + { \ddot { F } }$ where $\parallel F \parallel _ { 2 } \approx \mathbf { u }$ .

The order of the eigenvalues along $\hat { T }$ is somewhat arbitrary. But as we discuss in $\ S 7 . 6$ , any ordering can be achieved by using a simple procedure for swapping two adjacent diagonal entries.

# 7.5.7 Balancing

Finally, we mention that if the elements of A have widely varying magnitudes, then A should be balanced before applying the QR algorithm. This is an $O ( n ^ { 2 } )$ calculation in which a diagonal matrix D is computed so that if

$$
D ^ {- 1} A D = \left[ c _ {1} \mid \dots \mid c _ {n} \right] = \left[ \begin{array}{c} r _ {1} ^ {T} \\ \vdots \\ r _ {n} ^ {T} \end{array} \right]
$$

then $\parallel r _ { i } \parallel _ { \infty } \approx \parallel c _ { i } \parallel _ { \infty }$ for $i = 1 { : } n$ . The diagonal matrix D is chosen to have the form

$$
D = \mathrm{diag} (\beta^ {i _ {1}}, \ldots , \beta^ {i _ {n}})
$$

where $\beta$ is the floating point base. Note that $D ^ { - 1 } A D$ can be calculated without roundoff. When A is balanced, the computed eigenvalues are usually more accurate although there are exceptions. See Parlett and Reinsch (1969) and Watkins(2006).

# Problems

P7.5.1 Show that if $\bar { H } = Q ^ { T } H Q$ is obtained by performing a single-shift QR step with

$$
H = \left[ \begin{array}{c c} w & x \\ y & z \end{array} \right],
$$

then $| \bar { h } _ { 2 1 } | \le | y ^ { 2 } x | / [ ( w - z ) ^ { 2 } + y ^ { 2 } ] .$ .

P7.5.2 Given $A \in \mathbb { R } ^ { 2 \times 2 } ,$ , show how to compute a diagonal $D \in \mathbb { R } ^ { 2 \times 2 }$ so that $| \ D ^ { - 1 } A D \ | | _ { F }$ is minimized.

P7.5.3 Explain how the single-shift QR step H − $\mathbf { \nabla } \cdot \mu I = U R , \tilde { H } = R U + \mu I$ can be carried out implicitly. That is, show how the transition from H to H˜ can be carried out without subtracting the shift $\mu$ from the diagonal of H.

P7.5.4 Suppose H is upper Hessenberg and that we compute the factorization $P H = L U$ via Gaussian elimination with partial pivoting. (See Algorithm 4.3.4.) Show that $H _ { 1 } = U ( P ^ { T } L )$ is upper Hessenberg and similar to H. (This is the basis of the modified LR algorithm.)

P7.5.5 Show that if $H = H _ { 0 }$ is given and we generate the matrices $H _ { k }$ via $H _ { k } - \mu _ { k } I = U _ { k } R _ { k } , H _ { k + 1 }$ $= R _ { k } U _ { k } + \mu _ { k } I ,$ then $( U _ { 1 } \cdot \cdot \cdot U _ { j } ) ( R _ { j } \cdot \cdot \cdot R _ { 1 } ) = ( H - \mu _ { 1 } I ) \cdot \cdot \cdot ( H - \mu _ { j } I )$ .

# Notes and References for 7.5

Historically important papers associated with the QR iteration include:

H. Rutishauser (1958). “Solution of Eigenvalue Problems with the LR Transformation,” Nat. Bur. Stand. App. Math. Ser. 49, 47–81.

J.G.F. Francis (1961). “The QR Transformation: A Unitary Analogue to the LR Transformation, Parts I and II” Comput. J. 4, 265–72, 332–345.

V.N. Kublanovskaya (1961). “On Some Algorithms for the Solution of the Complete Eigenvalue Problem,” Vychisl. Mat. Mat. Fiz 1(4), 555–570.

R.S. Martin and J.H. Wilkinson (1968). “The Modified LR Algorithm for Complex Hessenberg Matrices,” Numer. Math. 12, 369–376.

R.S. Martin, G. Peters, and J.H. Wilkinson (1970). “The QR Algorithm for Real Hessenberg Matrices,” Numer. Math. 14, 219–231.

For a general insight, we recommend:

D.S. Watkins (1982). “Understanding the QR Algorithm,” SIAM Review 24, 427–440.

D.S. Watkins (1993). “Some Perspectives on the Eigenvalue Problem,” SIAM Review 35, 430–471.

D.S. Watkins (2008). “The QR Algorithm Revisited,” SIAM Review 50, 133–145.   
D.S. Watkins (2011). “Francis’s Algorithm,” Amer. Math. Monthly 118, 387–403.

Papers concerned with the convergence of the method, shifting, deflation, and related matters include:   
P.A. Businger (1971). “Numerically Stable Deflation of Hessenberg and Symmetric Tridiagonal Matrices, BIT 11, 262–270.   
D.S. Watkins and L. Elsner (1991). “Chasing Algorithms for the Eigenvalue Problem,” SIAM J. Matrix Anal. Applic. 12, 374–384.   
D.S. Watkins and L. Elsner (1991). “Convergence of Algorithms of Decomposition Type for the Eigenvalue Problem,” Lin. Alg. Applic. 143, 19–47.   
J. Erxiong (1992). “A Note on the Double-Shift QL Algorithm,” Lin. Alg. Applic. 171, 121–132.   
A.A. Dubrulle and G.H. Golub (1994). “A Multishift QR Iteration Without Computation of the Shifts,” Numer. Algorithms 7, 173–181.   
D.S. Watkins (1996). “Forward Stability and Transmission of Shifts in the QR Algorithm,” SIAM J. Matrix Anal. Applic. 16, 469–487.   
D.S. Watkins (1996). “The Transmission of Shifts and Shift Blurring in the QR algorithm,” Lin. Alg. Applic. 241–3, 877–896.   
D.S. Watkins (1998). “Bulge Exchanges in Algorithms of QR Type,” SIAM J. Matrix Anal. Applic. 19, 1074–1096.   
R. Vandebril (2011). “Chasing Bulges or Rotations? A Metamorphosis of the QR-Algorithm” SIAM. J. Matrix Anal. Applic. 32, 217–247.   
Aspects of the balancing problem are discussed in:   
E.E. Osborne (1960). “On Preconditioning of Matrices,” J. ACM 7, 338–345.   
B.N. Parlett and C. Reinsch (1969). “Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors,” Numer. Math. 13, 292–304.   
D.S. Watkins (2006). “A Case Where Balancing is Harmful,” ETNA 23, 1–4.   
Versions of the algorithm that are suitable for companion matrices are discussed in:   
D.A. Bini, F. Daddi, and L. Gemignani (2004). “On the Shifted QR iteration Applied to Companion Matrices,” ETNA 18, 137–152.   
M. Van Barel, R. Vandebril, P. Van Dooren, and K. Frederix (2010). “Implicit Double Shift QR-Algorithm for Companion Matrices,” Numer. Math. 116, 177–212.   
Papers that are concerned with the high-performance implementation of the QR iteration include:   
Z. Bai and J.W. Demmel (1989). “On a Block Implementation of Hessenberg Multishift QR Iteration,” Int. J. High Speed Comput. 1, 97–112.   
R.A. Van De Geijn (1993). “Deferred Shifting Schemes for Parallel QR Methods,” SIAM J. Matrix Anal. Applic. 14, 180–194.   
D.S. Watkins (1994). “Shifting Strategies for the Parallel QR Algorithm,” SIAM J. Sci. Comput. 15, 953–958.   
G. Henry and R. van de Geijn (1996). “Parallelizing the QR Algorithm for the Unsymmetric Algebraic Eigenvalue Problem: Myths and Reality,” SIAM J. Sci. Comput. 17, 870–883.   
Z. Bai, J. Demmel, J. Dongarra, A. Petitet, H. Robinson, and K. Stanley (1997). “The Spectral Decomposition of Nonsymmetric Matrices on Distributed Memory Parallel Computers,” SIAM J. Sci. Comput. 18, 1446–1461.   
G. Henry, D.S. Watkins, and J. Dongarra (2002). “A Parallel Implementation of the Nonsymmetric QR Algorithm for Distributed Memory Architectures,” SIAM J. Sci. Comput. 24, 284–311.   
K. Braman, R. Byers, and R. Mathias (2002). “The Multishift QR Algorithm. Part I: Maintaining Well-Focused Shifts and Level 3 Performance,” SIAM J. Matrix Anal. Applic. 23, 929–947.   
K. Braman, R. Byers, and R. Mathias (2002). “The Multishift QR Algorithm. Part II: Aggressive Early Deflation,” SIAM J. Matrix Anal. Applic. 23, 948–973.   
M.R. Fahey (2003). “Algorithm 826: A Parallel Eigenvalue Routine for Complex Hessenberg Matrices,” ACM Trans. Math. Softw. 29, 326–336.   
D. Kressner (2005). “On the Use of Larger Bulges in the QR Algorithm,” ETNA 20, 50–63.   
D. Kressner (2008). “The Effect of Aggressive Early Deflation on the Convergence of the QR Algorithm,” SIAM J. Matrix Anal. Applic. 30, 805–821.

# 7.6 Invariant Subspace Computations

Several important invariant subspace problems can be solved once the real Schur decomposition $Q ^ { T } A Q = T$ has been computed. In this section we discuss how to

• compute the eigenvectors associated with some subset of λ(A),   
• compute an orthonormal basis for a given invariant subspace,   
• block-diagonalize A using well-conditioned similarity transformations,   
• compute a basis of eigenvectors regardless of their condition, and   
• compute an approximate Jordan canonical form of A.

Eigenvector/invariant subspace computation for sparse matrices is discussed in §7.3.1 and §7.3.2 as well as portions of Chapters 8 and 10.

# 7.6.1 Selected Eigenvectors via Inverse Iteration

Let $q ^ { ( 0 ) } \in \mathbb { R } ^ { n }$ be a given unit 2-norm vector and assume that $A - \mu I \in \mathbb { R } ^ { n \times n }$ is nonsingular. The following is referred to as inverse iteration:

for k = 1, 2, . . .

$$
\text { Solve } (A - \mu I) z ^ {(k)} = q ^ {(k - 1)}.
$$

$$
q ^ {(k)} = z ^ {(k)} / \| z ^ {(k)} \| _ {2} \tag {7.6.1}
$$

$$
\lambda^ {(k)} = q ^ {(k) ^ {T}} A q ^ {(k)}
$$

end

Inverse iteration is just the power method applied to $( A - \mu I ) ^ { - 1 }$ .

To analyze the behavior of (7.6.1), assume that A has a basis of eigenvectors $\{ x _ { 1 } , \ldots , x _ { n } \}$ and that $A x _ { i } = \lambda _ { i } x _ { i }$ for i = 1:n. If

$$
q ^ {(0)} = \sum_ {i = 1} ^ {n} \beta_ {i} x _ {i}
$$

then $q ^ { ( k ) }$ is a unit vector in the direction of

$$
(A - \mu I) ^ {- k} q ^ {(0)} = \sum_ {i = 1} ^ {n} \frac {\beta_ {i}}{(\lambda_ {i} - \mu) ^ {k}} x _ {i}.
$$

Clearly, if $\mu$ is much closer to an eigenvalue $\lambda _ { j }$ than to the other eigenvalues, then $q ^ { ( k ) }$ is rich in the direction of $x _ { j }$ provided $\beta _ { j } \neq 0$ .

A sample stopping criterion for (7.6.1) might be to quit as soon as the residual

$$
r ^ {(k)} = (A - \mu I) q ^ {(k)}
$$

satisfies

$$
\| r ^ {(k)} \| _ {\infty} \leq c \mathbf {u} \| A \| _ {\infty} \tag {7.6.2}
$$

where c is a constant of order unity. Since

$$
(A + E _ {k}) q ^ {(k)} = \mu q ^ {(k)}
$$

with $E _ { k } = - r ^ { ( k ) } q ^ { ( k ) ^ { T } }$ , it follows that (7.6.2) forces $\mu$ and $q ^ { ( k ) }$ to be an exact eigenpair for a nearby matrix.

Inverse iteration can be used in conjunction with Hessenberg reduction and the QR algorithm as follows:

Step 1. Compute the Hessenberg decomposition $U _ { 0 } ^ { T } A U _ { 0 } = H .$ .

Step 2. Apply the double-implicit-shift Francis iteration to H without accumulating transformations.

Step 3. For each computed eigenvalue λ whose corresponding eigenvector x is sought, apply (7.6.1) with A = H and $\mu = \lambda$ to produce a vector z such that $H z \approx \mu z$ .

Step 4. Set $x = U _ { 0 } z$

Inverse iteration with H is very economical because we do not have to accumulate transformations during the double Francis iteration. Moreover, we can factor matrices of the form $H - \lambda I$ in $O ( n ^ { 2 } )$ flops, and (3) only one iteration is typically required to produce an adequate approximate eigenvector.

This last point is perhaps the most interesting aspect of inverse iteration and requires some justification since λ can be comparatively inaccurate if it is ill-conditioned. Assume for simplicity that λ is real and let

$$
H - \lambda I = \sum_ {i = 1} ^ {n} \sigma_ {i} u _ {i} v _ {i} ^ {T} = U \Sigma V ^ {T}
$$

be the SVD of $H - \lambda I$ . From what we said about the roundoff properties of the QR algorithm in $\ S 7 . 5 . 6 \AA$ , there exists a matrix $E \in \mathbb { R } ^ { n \times n }$ such that $H + E - \lambda I$ is singular and $\parallel E \parallel _ { 2 } \approx \mathbf { u } \parallel H \parallel _ { 2 }$ . It follows that $\sigma _ { n } \approx \mathbf { u } \sigma _ { 1 }$ and

$$
\| (H - \hat {\lambda} I) v _ {n} \| _ {2} \approx \mathbf {u} \sigma_ {1},
$$

i.e., $v _ { n }$ is a good approximate eigenvector. Clearly if the starting vector $q ^ { ( 0 ) }$ has the expansion

$$
q ^ {(0)} = \sum_ {i = 1} ^ {n} \gamma_ {i} u _ {i}
$$

then

$$
z ^ {(1)} = \sum_ {i = 1} ^ {n} \frac {\gamma_ {i}}{\sigma_ {i}} v _ {i}
$$

is “rich” in the direction $v _ { n }$ . Note that if $s ( \lambda ) \approx | u _ { n } ^ { T } v _ { n } |$ is small, then $z ^ { ( 1 ) }$ is rather deficient in the direction $u _ { n }$ . This explains (heuristically) why another step of inverse iteration is not likely to produce an improved eigenvector approximate, especially if λ is ill-conditioned. For more details, see Peters and Wilkinson (1979).

# 7.6.2 Ordering Eigenvalues in the Real Schur Form

Recall that the real Schur decomposition provides information about invariant subspaces. If

$$
Q ^ {T} A Q = T = \left[ \begin{array}{c c} T _ {1 1} & T _ {1 2} \\ 0 & T _ {2 2} \end{array} \right] _ {q} ^ {p}
$$

and

$$
\lambda (T _ {1 1}) \cap \lambda (T _ {2 2}) = \emptyset ,
$$

then the first $p$ columns of $Q$ span the unique invariant subspace associated with $\lambda ( T _ { 1 1 } )$ . (See $\ S 7 . 1 . 4 . )$ Unfortunately, the Francis iteration supplies us with a real Schur decomposition $Q _ { F } ^ { T } A Q _ { F } = T _ { F }$ in which the eigenvalues appear somewhat randomly along the diagonal of $T _ { F }$ . This poses a problem if we want an orthonormal basis for an invariant subspace whose associated eigenvalues are not at the top of $T _ { F } \mathrm { \ ' } _ { \mathrm { s } }$ diagonal. Clearly, we need a method for computing an orthogonal matrix $Q _ { D }$ such that $Q _ { D } ^ { T } T _ { F } Q _ { D }$ is upper quasi-triangular with appropriate eigenvalue ordering.

A look at the 2-by-2 case suggests how this can be accomplished. Suppose

$$
Q _ {F} ^ {T} A Q _ {F} = T _ {F} = \left[ \begin{array}{c c} \lambda_ {1} & t _ {1 2} \\ 0 & \lambda_ {2} \end{array} \right], \qquad \lambda_ {1} \neq \lambda_ {2}
$$

and that we wish to reverse the order of the eigenvalues. Note that

$$
T _ {F} x = \lambda_ {2} x
$$

where

$$
x = \left[ \begin{array}{c} t _ {1 2} \\ \lambda_ {2} - \lambda_ {1} \end{array} \right].
$$

Let $Q _ { D }$ be a Givens rotation such that the second component of $Q _ { D } ^ { T } x$ is zero. If

$$
Q = Q _ {F} Q _ {D},
$$

then

$$
(Q ^ {T} A Q) e _ {1} = Q _ {D} ^ {T} T _ {F} (Q _ {D} e _ {1}) = \lambda_ {2} Q _ {D} ^ {T} (Q _ {D} e _ {1}) = \lambda_ {2} e _ {1}.
$$

The matrices A and $Q ^ { T } A Q$ have the same Frobenius norm and so it follows that the latter must have the following form:

$$
Q ^ {T} A Q = \left[ \begin{array}{c c} \lambda_ {2} & \pm t _ {1 2} \\ 0 & \lambda_ {1} \end{array} \right].
$$

The swapping gets a little more complicated if $T$ has 2-by-2 blocks along its diagonal. See Ruhe (1970) and Stewart (1976) for details.

By systematically interchanging adjacent pairs of eigenvalues (or 2-by-2 blocks), we can move any subset of $\lambda ( A )$ to the top of $T \mathrm { { s } }$ diagonal. Here is the overall procedure for the case when there are no 2-by-2 bumps:

Algorithm 7.6.1 Given an orthogonal matrix $Q \in \mathbb { R } ^ { n \times n }$ , an upper triangular matrix ${ \bar { T } } { \bar { = } } Q ^ { T } A Q$ , and a subset $\Delta = \{ \lambda _ { 1 } , \ldots , \lambda _ { p } \}$ of $\lambda ( A )$ , the following algorithm computes an orthogonal matrix $Q _ { D }$ such that $Q _ { D } ^ { T } T \dot { Q } _ { D } = S$ is upper triangular and $\{ s _ { 1 1 } , \dotsc , s _ { p p } \}$ $= \Delta$ . The matrices $Q$ and $T$ are overwritten by $Q Q _ { D }$ and $S ,$ respectively.

while $\{ t _ { 1 1 } , \hdots , t _ { p p } \} \neq \Delta$

for $k = 1 { : } n - 1$

$\mathbf { i f } \ t _ { k k } \notin \Delta \ \mathrm { a n d } \ t _ { k + 1 , k + 1 } \in \Delta$

$$
[ c, s ] = \text { g   i   v   e   n   s } (T (k, k + 1), T (k + 1, k + 1) - T (k, k))
$$

$$
T (k: k + 1, k: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} T (k: k + 1, k: n)
$$

$$
T (1: k + 1, k: k + 1) = T (1: k + 1, k: k + 1) \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right]
$$

$$
Q (1: n, k: k + 1) = Q (1: n, k: k + 1) \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right]
$$

end end end

This algorithm requires $k ( 1 2 n )$ flops, where k is the total number of required swaps. The integer k is never greater than $( n - p ) p$ .

Computation of invariant subspaces by manipulating the real Schur decomposition is extremely stable. If $\hat { Q } = [ \hat { q } _ { 1 } | \cdot \cdot \cdot | \hat { q } _ { n } ]$ denotes the computed orthogonal matrix $Q .$ then $\Vert \hat { Q } ^ { T } \hat { Q } - I \Vert _ { 2 } \approx$ u and there exists a matrix E satisfying  $E \parallel _ { 2 } \approx \mathbf { u } \parallel A$ 2 such that $( A + E ) \hat { q } _ { i } \ \in \mathsf { s p a n } \{ \hat { q } _ { 1 } , \dots , \hat { q } _ { p } \}$ for $i = 1 { : } p$ .

# 7.6.3 Block Diagonalization

Let

$$
T = \left[ \begin{array}{c c c c} T _ {1 1} & T _ {1 2} & \dots & T _ {1 q} \\ 0 & T _ {2 2} & \dots & T _ {2 q} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & T _ {q q} \end{array} \right] \begin{array}{c} n _ {1} \\ n _ {2} \\ n _ {q} \end{array} \tag {7.6.3}
$$

be a partitioning of some real Schur canonical form $Q ^ { T } A Q \ = \ T \in \mathbb { R } ^ { n \times n }$ such that $\lambda ( T _ { 1 1 } ) , \dots , \lambda ( T _ { q q } )$ are disjoint. By Theorem 7.1.6 there exists a matrix Y such that

$$
Y ^ {- 1} T Y = \mathrm{diag} (T _ {1 1}, \ldots , T _ {q q}).
$$

A practical procedure for determining $Y$ is now given together with an analysis of $Y \mathrm { { } s }$ sensitivity as a function of the above partitioning.

Partition $I _ { n } = [ E _ { 1 } | \cdot \cdot \cdot | E _ { q } ]$ conformably with T and define the matrix $Y _ { i j } \in \mathbb { R } ^ { n \times n }$ as follows:

$$
Y _ {i j} = I _ {n} + E _ {i} Z _ {i j} E _ {j} ^ {T}, \quad i <   j, Z _ {i j} \in \mathbb {R} ^ {n _ {i} \times n _ {j}}.
$$

In other words, $Y _ { i j }$ looks just like the identity except that $Z _ { i j }$ occupies the $( i , j )$ block position. It follows that if $Y _ { i j } ^ { - 1 } T Y _ { i j } \ = \ \bar { T } = ( \bar { T } _ { i j } )$ , then $T$ and $\bar { T }$ are identical except that

$$
\begin{array}{l} \bar {T} _ {i j} = T _ {i i} Z _ {i j} - Z _ {i j} T _ {j j} + T _ {i j}, \\ \bar {T} _ {i k} = T _ {i k} - Z _ {i j} T _ {j k}, \quad (k = j + 1: q), \\ \bar {T} _ {k j} = T _ {k i} Z _ {i j} + T _ {k j}, \quad (k = 1: i - 1). \\ \end{array}
$$

Thus, $T _ { i j }$ can be zeroed provided we have an algorithm for solving the Sylvester equation

$$
F Z - Z G = C \tag {7.6.4}
$$

where $F \in \mathbb { R } ^ { p \times p }$ and $G \in \mathbb { R } ^ { r \times r }$ are given upper quasi-triangular matrices and $C \in \mathbb { R } ^ { p \times r }$

Bartels and Stewart (1972) have devised a method for doing this. Let $C =$ $\left[ c _ { 1 } \mid \cdots \mid c _ { r } \right]$ and $Z = \left[ \ : z _ { 1 } \ : | \cdots | \ : z _ { r } \ : \right]$ be column partitionings. If $g _ { k + 1 , k } = 0$ , then by comparing columns in (7.6.4) we find

$$
F z _ {k} - \sum_ {i = 1} ^ {k} g _ {i k} z _ {i} = c _ {k}.
$$

Thus, once we know $z _ { 1 } , \ldots , z _ { k - 1 }$ , then we can solve the quasi-triangular system

$$
(F - g _ {k k} I) z _ {k} = c _ {k} + \sum_ {i = 1} ^ {k - 1} g _ {i k} z _ {i}
$$

for $z _ { k }$ . If $g _ { k + 1 , k } \neq 0$ , then $z _ { k }$ and $z _ { k + 1 }$ can be simultaneously found by solving the $2 p { \mathrm { - } } \mathrm { b y - } 2 p$ system

$$
\left[ \begin{array}{c c} F - g _ {k k} I & - g _ {m k} I \\ - g _ {k m} I & F - g _ {m m} I \end{array} \right] \left[ \begin{array}{l} z _ {k} \\ z _ {m} \end{array} \right] = \left[ \begin{array}{l} c _ {k} \\ c _ {m} \end{array} \right] + \sum_ {i = 1} ^ {k - 1} \left[ \begin{array}{l} g _ {i k} z _ {i} \\ g _ {i m} z _ {i} \end{array} \right] \tag {7.6.5}
$$

where $m = k + 1$ . By reordering the equations according to the perfect shuffle permutation $( 1 , p + 1 , 2 , p + 2 , \ldots , p , 2 p )$ , a banded system is obtained that can be solved in $O ( p ^ { 2 } )$ flops. The details may be found in Bartels and Stewart (1972). Here is the overall process for the case when $F$ and $G$ are each triangular.

Algorithm 7.6.2 (Bartels-Stewart Algorithm) Given $C \in \mathbb { R } ^ { p \times r }$ and upper triangular matrices $\boldsymbol { F } \in \mathbb { R } ^ { p \times p }$ and $G \in \mathbb { R } ^ { r \times r }$ that satisfy $\lambda ( F ) \cap \lambda ( G ) = \emptyset$ , the following algorithm overwrites $C$ with the solution to the equation $F Z - Z G = C$ .

$$
\begin{array}{l} \text { for } k = 1: r \\ C (1: p, k) = C (1: p, k) + C (1: p, 1: k - 1) \cdot G (1: k - 1, k) \\ \text { Solve } (F - G (k, k) I) z = C (1: p, k) \text { for } z. \\ C (1: p, k) = z \\ \end{array}
$$

end

This algorithm requires $p r ( p + r )$ flops. By zeroing the superdiagonal blocks in $T$ in the appropriate order, the entire matrix can be reduced to block diagonal form.

Algorithm 7.6.3 Given an orthogonal matrix $Q \in \mathbb { R } ^ { n \times n }$ , an upper quasi-triangular matrix $T = Q ^ { T } A Q$ , and the partitioning (7.6.3), the following algorithm overwrites Q with QY where $Y ^ { - 1 } T Y = \mathrm { d i a g } ( T _ { 1 1 } , \dots , T _ { q q } )$ .

for j = 2:q

for i = 1:j − 1

Solve $T _ { i i } Z - Z T _ { j j } = - T _ { i j }$ for Z using the Bartels-Stewart algorithm.

for k = j + 1:q

$$
T _ {i k} = T _ {i k} - Z T _ {j k}
$$

end

for k = 1:q

$$
Q _ {k j} = Q _ {k i} Z + Q _ {k j}
$$

end

end end

The number of flops required by this algorithm is a complicated function of the block sizes in (7.6.3).

The choice of the real Schur form T and its partitioning in (7.6.3) determines the sensitivity of the Sylvester equations that must be solved in Algorithm 7.6.3. This in turn affects the condition of the matrix Y and the overall usefulness of the block diagonalization. The reason for these dependencies is that the relative error of the computed solution $\hat { Z }$ to

$$
T _ {i i} Z - Z T _ {j j} = - T _ {i j} \tag {7.6.6}
$$

satisfies

$$
\frac {\parallel \hat {Z} - Z \parallel_ {F}}{\parallel Z \parallel_ {F}} \approx \mathbf {u} \frac {\parallel T \parallel_ {F}}{\mathsf {s e p} (T _ {i i} , T _ {j j})}.
$$

For details, see Golub, Nash, and Van Loan (1979). Since

$$
\mathsf{sep}(T_{ii},T_{jj}) = \min_{X\neq 0}\frac{\parallel T_{ii}X - XT_{jj}\parallel_{F}}{\parallel X\parallel_{F}}\leq \min_{\substack{\lambda \in \lambda (T_{ii})\\ \mu \in \lambda (T_{jj})}}|\lambda -\mu |
$$

there can be a substantial loss of accuracy whenever the subsets $\lambda ( T _ { i i } )$ are insufficiently separated. Moreover, if Z satisfies (7.6.6) then

$$
\| Z \| _ {F} \leq \frac {\| T _ {i j} \| _ {F}}{\mathsf {s e p} (T _ {i i} , T _ {j j})}.
$$

Thus, large norm solutions can be expected if $\mathsf { s e p } ( T _ { i i } , T _ { j j } )$ is small. This tends to make the matrix Y in Algorithm 7.6.3 ill-conditioned since it is the product of the matrices

$$
Y _ {i j} = \left[ \begin{array}{c c} I _ {n _ {i}} & Z \\ 0 & I _ {n _ {j}} \end{array} \right].
$$

Note that $\kappa _ { F } ( Y _ { i j } ) = n _ { i } ^ { 2 } + n _ { j } ^ { 2 } + \left\| Z \right\| _ { F } ^ { 2 } .$ .

Confronted with these difficulties, Bavely and Stewart (1979) develop an algorithm for block diagonalizing that dynamically determines the eigenvalue ordering and partitioning in (7.6.3) so that all the Z matrices in Algorithm 7.6.3 are bounded in norm by some user-supplied tolerance. Their research suggests that the condition of Y can be controlled by controlling the condition of the $Y _ { i j }$ .

# 7.6.4 Eigenvector Bases

If the blocks in the partitioning (7.6.3) are all 1-by-1, then Algorithm 7.6.3 produces a basis of eigenvectors. As with the method of inverse iteration, the computed eigenvalueeigenvector pairs are exact for some “nearby” matrix. A widely followed rule of thumb for deciding upon a suitable eigenvector method is to use inverse iteration whenever fewer than 25% of the eigenvectors are desired.

We point out, however, that the real Schur form can be used to determine selected eigenvectors. Suppose

$$
Q ^ {T} A Q = \left[ \begin{array}{c c c} T _ {1 1} & u & T _ {1 3} \\ 0 & \lambda & v ^ {T} \\ 0 & 0 & T _ {3 3} \end{array} \right] \begin{array}{c} k - 1 \\ 1 \\ n - k \end{array}
$$

is upper quasi-triangular and that $\lambda \not \in \lambda ( T _ { 1 1 } ) \cup \lambda ( T _ { 3 3 } )$ . It follows that if we solve the linear systems $( T _ { 1 1 } - \lambda I ) w = - u$ and $( T _ { 3 3 } - \lambda I ) ^ { T } z = - v$ then

$$
x = Q {\left[ \begin{array}{l} w \\ 1 \\ 0 \end{array} \right]} \quad {\mathrm{and}} \quad y = Q {\left[ \begin{array}{l} 0 \\ 1 \\ z \end{array} \right]}
$$

are the associated right and left eigenvectors, respectively. Note that the condition of λ is prescribed by

$$
1 / s (\lambda) = \sqrt {(1 + w ^ {T} w) (1 + z ^ {T} z)}.
$$

# 7.6.5 Ascertaining Jordan Block Structures

Suppose that we have computed the real Schur decomposition $A = Q T Q ^ { T }$ , identified clusters of “equal” eigenvalues, and calculated the corresponding block diagonalization $T = Y { \cdot } \mathrm { d i a g } ( T _ { 1 1 } , \dots , T _ { q q } ) Y ^ { - 1 }$ . As we have seen, this can be a formidable task. However, even greater numerical problems confront us if we attempt to ascertain the Jordan block structure of each $T _ { i i }$ . A brief examination of these difficulties will serve to highlight the limitations of the Jordan decomposition.

Assume for clarity that $\lambda ( T _ { i i } )$ is real. The reduction of $T _ { i i }$ to Jordan form begins by replacing it with a matrix of the form $C = \lambda I + N$ , where N is the strictly upper triangular portion of $T _ { i i }$ and where λ, say, is the mean of its eigenvalues.

Recall that the dimension of a Jordan block $J ( \lambda )$ is the smallest nonnegative integer k for which $[ J ( \lambda ) - \lambda I ] ^ { k } = 0$ . Thus, if $p _ { i } = \mathsf { d i m } [ \mathsf { n u l l } ( N ^ { i } ) ]$ , for $i = 0 { : } n$ , then $p _ { i } - p _ { i - 1 }$ equals the number of blocks in $C ^ { \mathrm { { * } } } \mathrm { { s } }$ Jordan form that have dimension i or greater. A concrete example helps to make this assertion clear and to illustrate the role of the SVD in Jordan form computations.

Assume that C is 7-by-7. Suppose we compute the SVD $U _ { 1 } ^ { T } N V _ { 1 } = \Sigma _ { 1 }$ and “discover” that N has rank 3. If we order the singular values from small to large then it follows that the matrix $N _ { 1 } = V _ { 1 } ^ { T } N V _ { 1 }$ has the form

$$
N _ {1} = \left[ \begin{array}{c c} 0 & K \\ 0 & L \end{array} \right] _ {3} ^ {4}.
$$

At this point, we know that the geometric multiplicity of λ is $4 \mathrm { - } \mathrm { i . e } , C \mathrm { s }$ Jordan form has four blocks $( p _ { 1 } - p _ { 0 } = 4 - 0 = 4 )$ .

Now suppose $\tilde { U } _ { 2 } ^ { T } L \tilde { V } _ { 2 } = \Sigma _ { 2 }$ is the SVD of L and that we find that L has unit rank. If we again order the singular values from small to large, then $L _ { 2 } = \tilde { V } _ { 2 } ^ { T } L \tilde { V } _ { 2 }$ clearly has the following structure:

$$
L _ {2} = \left[ \begin{array}{c c c} 0 & 0 & a \\ 0 & 0 & b \\ 0 & 0 & c \end{array} \right].
$$

However, $\lambda ( L _ { 2 } ) = \lambda ( L ) = \{ 0 , 0 , 0 \}$ and so $c = 0$ . Thus, if

$$
V _ {2} = \mathrm{diag} (I _ {4}, \tilde {V} _ {2})
$$

then $N _ { 2 } = V _ { 2 } ^ { T } N _ { 1 } V _ { 2 }$ has the following form:

$$
N _ {2} = \left[ \begin{array}{c c c c c c c} 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & a \\ 0 & 0 & 0 & 0 & 0 & 0 & b \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{array} \right].
$$

Besides allowing us to introduce more zeros into the upper triangle, the SVD of L also enables us to deduce the dimension of the nullspace of $N ^ { 2 }$ . Since

$$
N _ {1} ^ {2} = \left[ \begin{array}{c c} 0 & K L \\ 0 & L ^ {2} \end{array} \right] = \left[ \begin{array}{c c} 0 & K \\ 0 & L \end{array} \right] \left[ \begin{array}{c c} 0 & K \\ 0 & L \end{array} \right]
$$

$\left[ \begin{array} { l } { K } \\ { L } \end{array} \right]$ has full column rank,

$$
p _ {2} = \dim (\operatorname{null} (N ^ {2})) = \dim (\operatorname{null} (N _ {1} ^ {2})) = 4 + \dim (\operatorname{null} (L)) = p _ {1} + 2.
$$

Hence, we can conclude at this stage that the Jordan form of C has at least two blocks of dimension 2 or greater.

Finally, it is easy to see that $N _ { 1 } ^ { 3 } = 0$ , from which we conclude that there is $p _ { 3 } - p _ { 2 }$ $= 7 - 6 = 1$ block of dimension 3 or larger. If we define $V = V _ { 1 } V _ { 2 }$ then it follows that

the decomposition

$$
V ^ {T} C V = \left[ \begin{array}{c c c c c c c} \lambda & 0 & 0 & 0 & \times & \times & \times \\ 0 & \lambda & 0 & 0 & \times & \times & \times \\ 0 & 0 & \lambda & 0 & \times & \times & \times \\ 0 & 0 & 0 & \lambda & \times & \times & \times \\ 0 & 0 & 0 & 0 & \lambda & \times & a \\ 0 & 0 & 0 & 0 & 0 & \lambda & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & \lambda \end{array} \right] \left\{ \begin{array}{l} \text {four blocks of order 1 or larger} \\ \text {two blocks of order 2 or larger} \\ \text {one block of order 3 or larger} \end{array} \right.
$$

displays $C ^ { \mathrm { { * } } } \mathrm { { s } }$ Jordan block structure: two blocks of order 1, one block of order 2, and one block of order 3.

To compute the Jordan decomposition it is necessary to resort to nonorthogonal transformations. We refer the reader to Golub and Wilkinson (1976), K˚agstr¨om and Ruhe (1980a, 1980b), and Demmel (1983) for more details. The above calculations with the SVD amply illustrate that difficult rank decisions must be made at each stage and that the final computed block structure depends critically on those decisions.

# Problems

P7.6.1 Give a complete algorithm for solving a real, n-by-n, upper quasi-triangular system $T x = b$

P7.6.2 Suppose $U ^ { - 1 } A U = \mathrm { d i a g } ( \alpha _ { 1 } , . . . , \alpha _ { m } )$ and $V ^ { - 1 } B V = \mathrm { d i a g } ( \beta _ { 1 } , . ~ . ~ . , \beta _ { n } )$ . Show that if

$$
\phi (X) = A X - X B,
$$

then

$$
\lambda (\phi) = \{\alpha_ {i} - \beta_ {j}: i = 1: m, j = 1: n \}.
$$

What are the corresponding eigenvectors? How can these facts be used to solve $A X - X B = C ?$

P7.6.3 Show that if $Z \in \mathbb { C } ^ { p \times q }$ and

$$
Y = \left[ \begin{array}{c c} I _ {p} & Z \\ 0 & I _ {q} \end{array} \right],
$$

then $\kappa _ { 2 } ( Y ) = [ 2 + \sigma ^ { 2 } + \sqrt { 4 \sigma ^ { 2 } + \sigma ^ { 4 } } ] / 2 \mathrm { ~ w h e r e ~ } \sigma = \| Z \| _ { 2 } .$

P7.6.4 Derive the system $_ { ( 7 . 6 . 5 ) }$ .

P7.6.5 Assume that $T \in \mathbb { R } ^ { n \times n }$ is block upper triangular and partitioned as follows:

$$
T = \left[ \begin{array}{c c c} T _ {1 1} & T _ {1 2} & T _ {1 3} \\ 0 & T _ {2 2} & T _ {2 3} \\ 0 & 0 & T _ {3 3} \end{array} \right], \qquad T \in \mathbb {R} ^ {n \times n}.
$$

Suppose that the diagonal block $T _ { 2 2 }$ is 2-by-2 with complex eigenvalues that are disjoint from $\lambda ( T _ { 1 1 } )$ and $\lambda ( T _ { 3 3 } )$ . Give an algorithm for computing the 2-dimensional real invariant subspace associated with $T _ { 2 2 } \mathrm { { ' s } }$ eigenvalues.

P7.6.6 Suppose $H \in \mathbb { R } ^ { n \times n }$ is upper Hessenberg with a complex eigenvalue $\lambda + i \cdot \mu$ . How could inverse iteration be used to compute x, $\boldsymbol { y } \in \mathbb { R } ^ { n }$ so that $H ( x + i y ) = ( \lambda + i \mu ) ( x + i y ) ?$ Hint: Compare real and imaginary parts in this equation and obtain a 2n-by-2n real system.

# Notes and References for §7.6

Much of the material discussed in this section may be found in the following survey paper:

G.H. Golub and J.H. Wilkinson (1976). “Ill-Conditioned Eigensystems and the Computation of the Jordan Canonical Form,” SIAM Review 18, 578–619.

The problem of ordering the eigenvalues in the real Schur form is the subject of:

A. Ruhe (1970). “An Algorithm for Numerical Determination of the Structure of a General Matrix,” BIT 10, 196–216.   
G.W. Stewart (1976). “Algorithm 406: HQR3 and EXCHNG: Fortran Subroutines for Calculating and Ordering the Eigenvalues of a Real Upper Hessenberg Matrix,” ACM Trans. Math. Softw. 2, 275–280.   
J.J. Dongarra, S. Hammarling, and J.H. Wilkinson (1992). “Numerical Considerations in Computing Invariant Subspaces,” SIAM J. Matrix Anal. Applic. 13, 145–161.   
Z. Bai and J.W. Demmel (1993). “On Swapping Diagonal Blocks in Real Schur Form,” Lin. Alg. Applic. 186, 73–95

Procedures for block diagonalization including the Jordan form are described in:

C. Bavely and G.W. Stewart (1979). “An Algorithm for Computing Reducing Subspaces by Block Diagonalization,” SIAM J. Numer. Anal. 16, 359–367.   
B. K˚agstr¨om and A. Ruhe (1980a). “An Algorithm for Numerical Computation of the Jordan Normal Form of a Complex Matrix,” ACM Trans. Math. Softw. 6, 398–419.   
B. K˚agstr¨om and A. Ruhe (1980b). “Algorithm 560 JNF: An Algorithm for Numerical Computation of the Jordan Normal Form of a Complex Matrix,” ACM Trans. Math. Softw. 6, 437–443.   
J.W. Demmel (1983). “A Numerical Analyst’s Jordan Canonical Form,” PhD Thesis, Berkeley.   
N. Ghosh, W.W. Hager, and P. Sarmah (1997). “The Application of Eigenpair Stability to Block Diagonalization,” SIAM J. Numer. Anal. 34, 1255–1268.   
S. Serra-Capizzano, D. Bertaccini, and G.H. Golub (2005). “How to Deduce a Proper Eigenvalue Cluster from a Proper Singular Value Cluster in the Nonnormal Case,” SIAM J. Matrix Anal. Applic. 27, 82–86.

Before we offer pointers to the literature associated with invariant subspace computation, we remind the reader that in 7.3 we discussed the power method for computing the dominant eigenpair and the method of orthogonal iteration that can be used to compute dominant invariant subspaces. Inverse iteration is a related idea and is the concern of the following papers:

J. Varah (1968). “The Calculation of the Eigenvectors of a General Complex Matrix by Inverse Iteration,” Math. Comput. 22, 785–791.   
J. Varah (1970). “Computing Invariant Subspaces of a General Matrix When the Eigensystem is Poorly Determined,” Math. Comput. 24, 137–149.   
G. Peters and J.H. Wilkinson (1979). “Inverse Iteration, Ill-Conditioned Equations, and Newton’s Method,” SIAM Review 21, 339–360.

I.C.F. Ipsen (1997). “Computing an Eigenvector with Inverse Iteration,” SIAM Review 39, 254–291. In certain applications it is necessary to track an invariant subspace as the matrix changes, see:

L. Dieci and M.J. Friedman (2001). “Continuation of Invariant Subspaces,” Num. Lin. Alg. 8, 317–327.   
D. Bindel, J.W. Demmel, and M. Friedman (2008). “Continuation of Invariant Subsapces in Large Bifurcation Problems,” SIAM J. Sci. Comput. 30, 637–656.

Papers concerned with estimating the error in a computed eigenvalue and/or eigenvector include:

S.P. Chan and B.N. Parlett (1977). “Algorithm 517: A Program for Computing the Condition Numbers of Matrix Eigenvalues Without Computing Eigenvectors,” ACM Trans. Math. Softw. 3, 186–203.   
H.J. Symm and J.H. Wilkinson (1980). “Realistic Error Bounds for a Simple Eigenvalue and Its Associated Eigenvector,” Numer. Math. 35, 113–126.   
C. Van Loan (1987). “On Estimating the Condition of Eigenvalues and Eigenvectors,” Lin. Alg. Applic. 88/89, 715–732.   
Z. Bai, J. Demmel, and A. McKenney (1993). “On Computing Condition Numbers for the Nonsymmetric Eigenproblem,” ACM Trans. Math. Softw. 19, 202–223.

Some ideas about improving computed eigenvalues, eigenvectors, and invariant subspaces may be found in:

J. Varah (1968). “Rigorous Machine Bounds for the Eigensystem of a General Complex Matrix,” Math. Comp. 22, 793–801.   
J.J. Dongarra, C.B. Moler, and J.H. Wilkinson (1983). “Improving the Accuracy of Computed Eigenvalues and Eigenvectors,” SIAM J. Numer. Anal. 20, 23–46.

J.W. Demmel (1987). “Three Methods for Refining Estimates of Invariant Subspaces,” Comput. 38, 43–57.   
As we have seen, the sep(.,.) function is of great importance in the assessment of a computed invariant subspace. Aspects of this quantity and the associated Sylvester equation are discussed in:   
J. Varah (1979). “On the Separation of Two Matrices,” SIAM J. Numer. Anal. 16, 212–222.   
R. Byers (1984). “A Linpack-Style Condition Estimator for the Equation $A X - X B ^ { T } = C , ^ { , }$ IEEE Trans. Autom. Contr. AC-29, 926–928.   
M. Gu and M.L. Overton (2006). “An Algorithm to Compute Sepλ,” SIAM J. Matrix Anal. Applic. 28, 348–359.   
N.J. Higham (1993). “Perturbation Theory and Backward Error for $\mathrm { A X - X B } = \mathrm { C } , \mathrm { " } \ B I T \ 3 3 , 1 2 4 \mathrm { - } 1 3 6 .$ .   
Sylvester equations arise in many settings, and there are many solution frameworks, see:   
R.H. Bartels and G.W. Stewart (1972). “Solution of the Equation $A X + X B = C , "$ Commun. ACM 15, 820–826.   
G.H. Golub, S. Nash, and C. Van Loan (1979). “A Hessenberg-Schur Method for the Matrix Problem AX + XB = C,” IEEE Trans. Autom. Contr. AC-24, 909–913.   
K. Datta (1988). “The Matrix Equation XA − BX = R and Its Applications,” Lin. Alg. Applic. 109, 91–105.   
B. K˚agstr¨om and P. Poromaa (1992). “Distributed and Shared Memory Block Algorithms for the Triangular Sylvester Equation with $\mathrm { s e p } ^ { - 1 }$ Estimators,” SIAM J. Matrix Anal. Applic. 13, 90– 101.   
J. Gardiner, M.R. Wette, A.J. Laub, J.J. Amato, and C.B. Moler (1992). “Algorithm 705: A FORTRAN-77 Software Package for Solving the Sylvester Matrix Equation $A X B ^ { T } + C X D ^ { T } = E ,$ ACM Trans. Math. Softw. 18, 232–238.   
V. Simoncini (1996). “On the Numerical Solution of AX -XB =C,” BIT 36, 814–830.   
C.H. Bischof, B.N Datta, and A. Purkayastha (1996). “A Parallel Algorithm for the Sylvester Observer Equation,” SIAM J. Sci. Comput. 17, 686–698.   
D. Calvetti, B. Lewis, L. Reichel (2001). “On the Solution of Large Sylvester-Observer Equations,” Num. Lin. Alg. 8, 435–451.

The constrained Sylvester equation problem is considered in:

J.B. Barlow, M.M. Monahemi, and D.P. O’Leary (1992). “Constrained Matrix Sylvester Equations,” SIAM J. Matrix Anal. Applic. 13, 1–9.

A.R. Ghavimi and A.J. Laub (1996). “Numerical Methods for Nearly Singular Constrained Matrix Sylvester Equations.” SIAM J. Matrix Anal. Applic. 17, 212–221.

The Lyapunov problem F X $+ X F ^ { T } = - C$ where C is non-negative definite has a very important role to play in control theory, see:

G. Hewer and C. Kenney (1988). “The Sensitivity of the Stable Lyapunov Equation,” SIAM J. Control Optim 26, 321–344.

A.R. Ghavimi and A.J. Laub (1995). “Residual Bounds for Discrete-Time Lyapunov Equations,” IEEE Trans. Autom. Contr. 40, 1244–1249.

J.-R. Li and J. White (2004). “Low-Rank Solution of Lyapunov Equations,” SIAM Review 46, 693– 713.

Several authors have considered generalizations of the Sylvester equation, i.e., $\Sigma F _ { i } X G _ { i } = C .$ These include:

P. Lancaster (1970). “Explicit Solution of Linear Matrix Equations,” SIAM Review 12, 544–566.

H. Wimmer and A.D. Ziebur (1972). “Solving the Matrix Equations $\Sigma f _ { p } ( A ) g _ { p } ( A ) = C , "$ SIAM Review 14, 318–323.

W.J. Vetter (1975). “Vector Structures and Solutions of Linear Matrix Equations,” Lin. Alg. Applic. 10, 181–188.

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

# 7.7.1 Background

The first thing to observe about the generalized eigenvalue problem is that there are n eigenvalues if and only if rank $( B ) = n$ . If B is rank deficient then $\lambda ( A , B )$ may be finite, empty, or infinite:

$$
A   =   \left[ \begin{array}{c c} 1 & 2 \\ 0 & 3 \end{array} \right],    B   =   \left[ \begin{array}{c c} 1 & 0 \\ 0 & 0 \end{array} \right]    \Rightarrow    \lambda (A, B) = \{1 \},
$$

$$
A   =   \left[ \begin{array}{c c} 1 & 2 \\ 0 & 3 \end{array} \right], \quad B   =   \left[ \begin{array}{c c} 0 & 1 \\ 0 & 0 \end{array} \right] \quad \Rightarrow \quad \lambda (A, B) = \emptyset ,
$$

$$
A   =   \left[ \begin{array}{c c} 1 & 2 \\ 0 & 0 \end{array} \right],    B   =   \left[ \begin{array}{c c} 1 & 0 \\ 0 & 0 \end{array} \right]    \Rightarrow    \lambda (A, B) = \mathbb {C}.
$$

Note that if $0 \neq \lambda \in \lambda ( A , B )$ , then $( 1 / \lambda ) \in \lambda ( B , A )$ . Moreover, if B is nonsingular, then $\lambda ( A , B ) = \lambda ( B ^ { - 1 } A , I ) = \lambda ( B ^ { - 1 } A )$ . This last observation suggests one method for solving the $A - \lambda B$ problem if B is nonsingular:

Step 1. Solve $B C = A$ for C using (say) Gaussian elimination with pivoting.

Step 2. Use the QR algorithm to compute the eigenvalues of C.

In this framework, C is affected by roundoff errors of order $\mathbf { u } \parallel A \parallel _ { 2 } \parallel B ^ { - 1 } \parallel _ { 2 }$ . If B is illconditioned, then this precludes the possibility of computing any generalized eigenvalue accurately—even those eigenvalues that may be regarded as well-conditioned. For example, if

$$
A = \left[ \begin{array}{l l} 1. 7 4 6 & . 9 4 0 \\ 1. 2 4 6 & 1. 8 9 8 \end{array} \right] \qquad \text {and} \qquad B = \left[ \begin{array}{l l}. 7 8 0 & . 5 6 3 \\ . 9 1 3 & . 6 5 9 \end{array} \right],
$$

then $\lambda ( A , B ) = \{ 2 , 1 . 0 7 \times 1 0 ^ { 6 } \}$ . With 7-digit floating point arithmetic, we find $\lambda ( \mathrm { f l } ( A B ^ { - 1 } ) ) = \{ 1 . 5 6 2 5 3 9 , 1 . 0 1 \times 1 0 ^ { 6 } \}$ . The poor quality of the small eigenvalue is because $\kappa _ { 2 } ( B ) \approx 2 \times 1 0 ^ { 6 }$ . On the other hand, we find that

$$
\lambda (I, \mathsf {f l} (A ^ {- 1} B)) \approx \{2. 0 0 0 0 0 1, 1. 0 6 \times 1 0 ^ {6} \}.
$$

The accuracy of the small eigenvalue is improved because $\kappa _ { 2 } ( A ) \approx 4$

The example suggests that we seek an alternative approach to the generalized eigenvalue problem. One idea is to compute well-conditioned Q and Z such that the matrices

$$
A _ {1} = Q ^ {- 1} A Z, \quad B _ {1} = Q ^ {- 1} B Z \tag {7.7.2}
$$

are each in canonical form. Note that $\lambda ( A , B ) { = } \lambda ( A _ { 1 } , B _ { 1 } )$ since

$$
A x = \lambda B x \Leftrightarrow A _ {1} y = \lambda B _ {1} y, x = Z y.
$$

We say that the pencils A − λB and $A _ { 1 } - \lambda B _ { 1 }$ are equivalent if (7.7.2) holds with nonsingular Q and Z.

As in the standard eigenproblem A − λI there is a choice between canonical forms. Corresponding to the Jordan form is a decomposition of Kronecker in which both $A _ { 1 }$ and $B _ { 1 }$ are block diagonal with blocks that are similar in structure to Jordan blocks. The Kronecker canonical form poses the same numerical challenges as the Jordan form, but it provides insight into the mathematical properties of the pencil $A - \lambda B$ . See Wilkinson (1978) and Demmel and K˚agstr¨om (1987) for details.

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

# 7.7.3 Sensitivity Issues

The generalized Schur decomposition sheds light on the issue of eigenvalue sensitivity for the $A - \lambda B$ problem. Clearly, small changes in A and B can induce large changes in the eigenvalue $\lambda _ { i } = t _ { i i } / s _ { i i }$ if $s _ { i i }$ is small. However, as Stewart (1978) argues, it may not be appropriate to regard such an eigenvalue as “ill-conditioned.” The reason is that the reciprocal $\mu _ { i } = s _ { i i } / t _ { i i }$ might be a very well-behaved eigenvalue for the pencil $\mu A - B$ . In the Stewart analysis, A and B are treated symmetrically and the eigenvalues are regarded more as ordered pairs $( t _ { i i } , s _ { i i } )$ than as quotients. With this point of view it becomes appropriate to measure eigenvalue perturbations in the chordal metric chord(a, b) defined by

$$
\operatorname{chord} (a, b) = \frac {| a - b |}{\sqrt {1 + a ^ {2}} \sqrt {1 + b ^ {2}}}.
$$

Stewart shows that if λ is a distinct eigenvalue of A − λB and $\lambda _ { \epsilon }$ is the corresponding eigenvalue of the perturbed pencil $\tilde { A } - \lambda \tilde { B }$ with $\parallel A - \tilde { A } \parallel _ { 2 } \approx \parallel B - \tilde { B } \parallel _ { 2 } \approx \epsilon ,$ , then

$$
\operatorname{chord} (\lambda , \lambda_ {\epsilon}) \leq \frac {\epsilon}{\sqrt {(y ^ {H} A x) ^ {2} + (y ^ {H} B x) ^ {2}}} + O (\epsilon^ {2})
$$

where x and y have unit 2-norm and satisfy Ax = λBx and $y ^ { H } A { = } \lambda y ^ { H } B$ . Note that the denominator in the upper bound is symmetric in A and B. The “truly” ill-conditioned eigenvalues are those for which this denominator is small.

The extreme case when both $t _ { k k }$ and $s _ { k k }$ are zero for some k has been studied by Wilkinson (1979). In this case, the remaining quotients $t _ { i i } / s _ { i i }$ can take on arbitrary values.

# 7.7.4 Hessenberg-Triangular Form

The first step in computing the generalized real Schur decomposition of the pair $( A , B )$ is to reduce A to upper Hessenberg form and B to upper triangular form via orthogonal transformations. We first determine an orthogonal U such that $U ^ { T } B$ is upper triangular. Of course, to preserve eigenvalues, we must also update A in exactly the same way. Let us trace what happens in the $n = 5 ~ \mathrm { c a s e }$ .

$$
A \gets U ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \end{array} \right], \quad B \gets U ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

Next, we reduce A to upper Hessenberg form while preserving $B ^ { \prime } \mathrm { s }$ upper triangular form. First, a Givens rotation $Q _ { 4 5 }$ is determined to zero $a _ { 5 1 }$ :

$$
A \gets Q _ {4 5} ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], \quad B \gets Q _ {4 5} ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right].
$$

The nonzero entry arising in the (5,4) position in B can be zeroed by postmultiplying with an appropriate Givens rotation $Z _ { 4 5 }$ :

$$
A \gets A Z _ {4 5} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], B \gets B Z _ {4 5} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

Zeros are similarly introduced into the (4, 1) and (3, 1) positions in A:

$$
A \gets Q _ {3 4} ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], B \gets Q _ {3 4} ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
A \gets A Z _ {3 4} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], \quad B \gets B Z _ {3 4} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
A \gets Q _ {2 3} ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], \quad B \gets Q _ {2 3} ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
A \gets A Z _ {2 3} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \end{array} \right], B \gets B Z _ {2 3} = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

A is now upper Hessenberg through its first column. The reduction is completed by zeroing $a _ { 5 2 } , ~ a _ { 4 2 }$ , and $a _ { 5 3 }$ . Note that two orthogonal transformations are required for each $a _ { i j }$ that is zeroed—one to do the zeroing and the other to restore $B \mathrm { { ^ { * } s } }$ triangularity. Either Givens rotations or 2-by-2 modified Householder transformations can be used. Overall we have:

Algorithm 7.7.1 (Hessenberg-Triangular Reduction) Given A and B in $\mathbb { R } ^ { n \times n }$ , the following algorithm overwrites A with an upper Hessenberg matrix $Q ^ { T } A Z$ and $B$ with an upper triangular matrix $Q ^ { T } B Z$ where both $Q$ and $Z$ are orthogonal.

Compute the factorization $B = Q R$ using Algorithm 5.2.1 and overwrite

A with $Q ^ { T } A$ and B with $Q ^ { T } B$

for $j = 1 { : } n - 2$

for $i = n \colon - 1 { : } j + 2$

$$
[ c, s ] = \text {   givens   } (A (i - 1, j), A (i, j))
$$

$$
A (i - 1: i, j: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A (i - 1: i, j: n)
$$

$$
B (i - 1: i, i - 1: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} B (i - 1: i, i - 1: n)
$$

$$
[ c, s ] = \text {   givens } (- B (i, i), B (i, i - 1))
$$

$$
B (1: i, i - 1: i) = B (1: i, i - 1: i) \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right]
$$

$$
A (1: n, i - 1: i) = A (1: n, i - 1: i) \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right]
$$

end

end

This algorithm requires about $8 n ^ { 3 }$ flops. The accumulation of $Q$ and $Z$ requires about $4 n ^ { 3 }$ and $3 n ^ { 3 }$ flops, respectively.

The reduction of $A - \lambda B$ to Hessenberg-triangular form serves as a “front end” decomposition for a generalized QR iteration known as the $\mathrm { Q Z }$ iteration which we describe next.

# 7.7.5 Deflation

In describing the QZ iteration we may assume without loss of generality that A is an unreduced upper Hessenberg matrix and that B is a nonsingular upper triangular matrix. The first of these assertions is obvious, for if $a _ { k + 1 , k } = 0$ then

$$
A - \lambda B = \left[ \begin{array}{c c} A _ {1 1} - \lambda B _ {1 1} & A _ {1 2} - \lambda B _ {1 2} \\ 0 & A _ {2 2} - \lambda B _ {2 2} \end{array} \right] _ {n - k} ^ {k},
$$

and we may proceed to solve the two smaller problems $A _ { 1 1 } - \lambda B _ { 1 1 }$ and $A _ { 2 2 } - \lambda B _ { 2 2 }$ . On the other hand, if $b _ { k k } = 0$ for some k, then it is possible to introduce a zero in A’s $( n , n - 1 )$ position and thereby deflate. Illustrating by example, suppose $n = 5$ and $k = 3 \colon$

$$
A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right], \qquad B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

The zero on $B ^ { \prime } \mathrm { s }$ diagonal can be “pushed down” to the (5,5) position as follows using Givens rotations:

$$
A \gets Q _ {3 4} ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \gets Q _ {3 4} ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
A \leftarrow A Z _ {2 3} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \leftarrow B Z _ {2 3} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right],
$$

$$
A \gets Q _ {4 5} ^ {T} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & \times & \times & \times \end{array} \right], \quad B \gets Q _ {4 5} ^ {T} B = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 \end{array} \right],
$$

$$
A \leftarrow A Z _ {3 4} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \leftarrow B Z _ {3 4} ^ {T} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 \end{array} \right],
$$

$$
A \leftarrow A Z _ {4 5} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times \end{array} \right], \quad B \leftarrow B Z _ {4 5} = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 \end{array} \right].
$$

This zero-chasing technique is perfectly general and can be used to zero $^ { a _ { n , n - 1 } }$ regardless of where the zero appears along $B ^ { \prime } \mathrm { s }$ diagonal.

# 7.7.6 The QZ Step

We are now in a position to describe a QZ step. The basic idea is to update A and B as follows

$$
(\bar {A} - \lambda \bar {B}) = \bar {Q} ^ {T} (A - \lambda B) \bar {Z},
$$

where $\bar { A }$ is upper Hessenberg, $\bar { B }$ is upper triangular, $\bar { Q }$ and $\bar { Z }$ are each orthogonal, and $\bar { A } \bar { B } ^ { - 1 }$ is essentially the same matrix that would result if a Francis QR step (Algorithm 7.5.1) were explicitly applied to $A B ^ { - 1 }$ . This can be done with some clever zero-chasing and an appeal to the implicit Q theorem.

Let $M = A B ^ { - 1 }$ (upper Hessenberg) and let v be the first column of the matrix $( M - a I ) ( M - b I )$ , where a and b are the eigenvalues of $M \mathrm { { s } }$ lower 2-by-2 submatrix. Note that v can be calculated in $O ( 1 )$ flops. If $P _ { 0 }$ is a Householder matrix such that $P _ { 0 } v$ is a multiple of $e _ { 1 }$ , then

$$
A \leftarrow P _ {0} A = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \leftarrow P _ {0} B = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

The idea now is to restore these matrices to Hessenberg-triangular form by chasing the unwanted nonzero elements down the diagonal.

To this end, we first determine a pair of Householder matrices $Z _ { 1 }$ and $Z _ { 2 }$ to zero $b _ { 3 1 } , \ b _ { 3 2 }$ , and $b _ { 2 1 }$ :

$$
A \leftarrow A Z _ {1} Z _ {2} = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \leftarrow B Z _ {1} Z _ {2} = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

Then a Householder matrix $P _ { 1 }$ is used to zero $a _ { 3 1 }$ and $a _ { 4 1 }$ :

$$
A \gets P _ {1} A = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \gets P _ {1} B = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

Notice that with this step the unwanted nonzero elements have been shifted down and to the right from their original position. This illustrates a typical step in the QZ iteration. Notice that $Q = Q _ { 0 } Q _ { 1 } \cdot \cdot \cdot Q _ { n - 2 }$ has the same first column as $Q _ { 0 }$ . By the way the initial Householder matrix was determined, we can apply the implicit Q theorem and assert that $A B ^ { - 1 } = Q ^ { T } ( A B ^ { - 1 } ) Q$ is indeed essentially the same matrix that we would obtain by applying the Francis iteration to $M = A B ^ { - 1 }$ directly. Overall we have the following algorithm.

Algorithm 7.7.2 (The QZ Step) Given an unreduced upper Hessenberg matrix $A \in \mathbb { R } ^ { n \times n }$ and a nonsingular upper triangular matrix $B \in \mathbb { R } ^ { n \times n }$ , the following algorithm overwrites A with the upper Hessenberg matrix $Q ^ { T } A Z$ and $B$ with the upper triangular matrix $Q ^ { T } B Z$ where $Q$ and $Z$ are orthogonal and $Q$ has the same first column as the orthogonal similarity transformation in Algorithm 7.5.1 when it is applied to $A B ^ { - 1 }$ .

Let $M = A B ^ { - 1 }$ and compute $( M - a I ) ( M - b I ) e _ { 1 } = [ x , y , z , 0 , \ldots , 0 ] ^ { T }$

where a and b are the eigenvalues of $M \mathrm { { s } }$ lower 2-by-2.

for $k = 1 { : } n - 2$

$\mathrm { F i n d ~ H o u s e h o l d e r } Q _ { k } \mathrm { \ s o } Q _ { k } \left[ \begin{array} { l } { x } \\ { y } \\ { z } \end{array} \right] = \left[ \begin{array} { l } { * } \\ { 0 } \\ { 0 } \end{array} \right] .$

$$
A = \operatorname{diag} (I _ {k - 1}, Q _ {k}, I _ {n - k - 2}) \cdot A
$$

$$
B = \operatorname{diag} (I _ {k - 1}, Q _ {k}, I _ {n - k - 2}) \cdot B
$$

$\mathrm { F i n d ~ H o u s e h o l d e r ~ } Z _ { k 1 } \mathrm { ~ s o ~ } \left[ \begin{array} { l } { b _ { k + 2 , k } \ \big \vert \ b _ { k + 2 , k + 1 } \ \big \vert \ b _ { k + 2 , k + 2 } } \end{array} \right] Z _ { k 1 } = \left[ \begin{array} { l }  0 \ \big \vert \ 0 \ \big \vert \ \ast \ \right] . \end{array}$

$$
A = A \cdot \operatorname{diag} \left(I _ {k - 1}, Z _ {k 1}, I _ {n - k - 2}\right)
$$

$$
B = B \cdot \operatorname{diag} \left(I _ {k - 1}, Z _ {k 1}, I _ {n - k - 2}\right)
$$

$\mathrm { F i n d ~ H o u s e h o l d e r ~ } Z _ { k 2 } \mathrm { ~ s o ~ \ } \left[ \begin{array} { l } { b _ { k + 1 , k } \ | \ b _ { k + 1 , k + 1 } } \end{array} \right] Z _ { k 2 } = \left[ \begin{array} { l } { 0 \ | \ * \ } \end{array} \right] .$

$$
A = A \cdot \operatorname{diag} \left(I _ {k - 1}, Z _ {k 2}, I _ {n - k - 1}\right)
$$

$$
B = B \cdot \operatorname{diag} \left(I _ {k - 1}, Z _ {k 2}, I _ {n - k - 1}\right)
$$

$$
x = a _ {k + 1, k}; y = a _ {k + 2, k}
$$

$\mathbf { i f } \ k < n - 2$

$$
z = a _ {k + 3, k}
$$

end

end

Find Householder Qn−1 so Qn−1  xy 	 =  ∗0 	 .

$$
A = \operatorname{diag} (I _ {n - 2}, Q _ {n - 1}) \cdot A
$$

$$
B = \operatorname{diag} (I _ {n - 2}, Q _ {n - 1}) \cdot B.
$$

$\mathrm { F i n d ~ H o u s e h o l d e r ~ } Z _ { n - 1 } \mathrm { ~ s o ~ } \left[ \textit { b } _ { n , n - 1 } \mid \textit { b } _ { n n } \right] Z _ { n - 1 } = \left[ \textit { 0 } \mid * \right] .$

$$
A = A \cdot \operatorname{diag} \left(I _ {n - 2}, Z _ {n - 1}\right)
$$

$$
B = B \cdot \mathrm{diag} (I _ {n - 2}, Z _ {n - 1})
$$

This algorithm requires $2 2 n ^ { 2 }$ flops. $Q$ and $Z$ can be accumulated for an additional $8 n ^ { 2 }$ flops and $1 3 n ^ { 2 }$ flops, respectively.

# 7.7.7 The Overall QZ Process

By applying a sequence of QZ steps to the Hessenberg-triangular pencil $A - \lambda B$ , it is possible to reduce A to quasi-triangular form. In doing this it is necessary to monitor A’s subdiagonal and B’s diagonal in order to bring about decoupling whenever possible. The complete process, due to Moler and Stewart (1973), is as follows:

Algorithm 7.7.3 Given $A \in \mathbb { R } ^ { n \times n }$ and $B \in \mathbb { R } ^ { n \times n }$ , the following algorithm computes orthogonal Q and Z such that $Q ^ { T } A Z = T$ is upper quasi-triangular and $Q ^ { T } B \bar { Z } = S$ is upper triangular. A is overwritten by T and B by S.

Using Algorithm 7.7.1, overwrite A with $Q ^ { T } A Z$ (upper Hessenberg) and B with $Q ^ { T } B Z$ (upper triangular).

$$
\text { until } q = n
$$

Set to zero subdiagonal entries that satisfy $| a _ { i , i - 1 } | \leq \epsilon ( | a _ { i - 1 , i - 1 } | + | a _ { i i } | )$ .

Find the largest nonnegative q and the smallest nonnegative p such that if

$$
A = \left[ \begin{array}{c c c} A _ {1 1} & A _ {1 2} & A _ {1 3} \\ 0 & A _ {2 2} & A _ {2 3} \\ 0 & 0 & A _ {3 3} \end{array} \right] \begin{array}{c} p \\ n - p - q \\ q \end{array}
$$

then $A _ { 3 3 }$ is upper quasi-triangular and $A _ { 2 2 }$ is upper Hessenberg and unreduced.

Partition B conformably:

$$
B = \left[ \begin{array}{c c c} B _ {1 1} & B _ {1 2} & B _ {1 3} \\ 0 & B _ {2 2} & B _ {2 3} \\ 0 & 0 & B _ {3 3} \end{array} \right] \begin{array}{c} p \\ n - p - q \\ q \end{array}
$$

if $q < n$

if $B _ { 2 2 }$ is singular

Zero $a _ { n - q , n - q - 1 }$

Apply Algorithm 7.7.2 to $A _ { 2 2 }$ and $B _ { 2 2 }$ and update:

$$
A = \mathrm{diag} (I _ {p}, Q, I _ {q}) ^ {T} A \cdot \mathrm{diag} (I _ {p}, Z, I _ {q})
$$

$$
B = \mathrm{diag} (I _ {p}, Q, I _ {q}) ^ {T} B \cdot \mathrm{diag} (I _ {p}, Z, I _ {q})
$$

end end end

This algorithm requires $3 0 n ^ { 3 }$ flops. If $Q$ is desired, an additional $1 6 n ^ { 3 }$ are necessary. If Z is required, an additional $2 0 \bar { n } ^ { 3 }$ are needed. These estimates of work are based on the experience that about two $\mathrm { Q Z }$ iterations per eigenvalue are necessary. Thus, the convergence properties of $\mathrm { Q Z }$ are the same as for QR. The speed of the QZ algorithm is not affected by rank deficiency in B.

The computed S and T can be shown to satisfy

$$
Q _ {0} ^ {T} (A + E) Z _ {0} = T, \quad Q _ {0} ^ {T} (B + F) Z _ {0} = S,
$$

where $Q _ { 0 }$ and $Z _ { 0 }$ are exactly orthogonal and $\parallel E \parallel _ { 2 } \approx \mathbf { u } \parallel A \parallel _ { 2 }$ and $\parallel F \parallel _ { 2 } \approx \mathbf { u } \parallel B \parallel _ { 2 }$ .

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

# 7.7.9 A Note on the Polynomial Eigenvalue Problem

More general than the generalized eigenvalue problem is the polynomial eigenvalue problem. Here we are given matrices $A _ { 0 } , \dots , A _ { d } \in \mathbb { C } ^ { n \times n }$ and determine $\lambda \in \mathbb { C }$ and $0 \neq x \in \mathbb { C } ^ { n }$ so that

$$
P (\lambda) x = 0 \tag {7.7.3}
$$

where the λ-matrix P (λ) is defined by

$$
P (\lambda) = A _ {0} + \lambda A _ {1} + \dots + \lambda^ {d} A _ {d}. \tag {7.7.4}
$$

We assume $A _ { d } \neq 0$ and regard d as the degree of $P ( \lambda )$ . The theory behind the polynomial eigenvalue problem is nicely developed in Lancaster (1966).

It is possible to convert (7.7.3) into an equivalent linear eigenvalue problem with larger dimension. For example, suppose d = 3 and

$$
L (\lambda) = \left[ \begin{array}{c c c} 0 & 0 & A _ {0} \\ - I & 0 & A _ {1} \\ 0 & - I & A _ {2} \end{array} \right] + \lambda \left[ \begin{array}{c c c} I & 0 & 0 \\ 0 & I & 0 \\ 0 & 0 & A _ {3} \end{array} \right]. \tag {7.7.5}
$$

If

$$
L (\lambda) \left[ \begin{array}{l} u _ {1} \\ u _ {2} \\ x \end{array} \right] = \left[ \begin{array}{l} 0 \\ 0 \\ 0 \end{array} \right],
$$

then

$$
0 = A _ {0} x + \lambda u _ {1} = A _ {0} + \lambda \left(A _ {1} x + \lambda u _ {2}\right) = A _ {0} + \lambda \left(A _ {1} x + \lambda \left(A _ {2} + \lambda A _ {3}\right)\right) x = P (\lambda) x.
$$

In general, we say that L(λ) is a linearization of P (λ) if there are dn-by-dn λ-matrices S(λ) and T (λ), each with constant nonzero determinants, so that

$$
S (\lambda) \left[ \begin{array}{c c} P (\lambda) & 0 \\ 0 & I _ {(d - 1) n} \end{array} \right] T (\lambda) = L (\lambda) \tag {7.7.6}
$$

has unit degree. With this conversion, the A − λB methods just discussed can be applied to find the required eigenvalues and eigenvectors.

Recent work has focused on how to choose the λ-transformations S(λ) and $T ( \lambda )$ so that special structure in P (λ) is reflected in L(λ). See Mackey, Mackey, Mehl, and Mehrmann (2006). The idea is to think of (7.7.6) as a factorization and to identify the transformations that produce a properly structured L(λ). To appreciate this solution framework it is necessary to have a facility with λ-matrix manipulation and to that end we briefly examine the λ-matrix transformations behind the above linearization. If

$$
P _ {1} (\lambda) = A _ {1} + \lambda A _ {2} + \dots + \lambda^ {d - 1} A _ {d}
$$

then

$$
P (\lambda) = A _ {0} + \lambda P _ {1} (\lambda)
$$

and it is easy to verify that

$$
\left[ \begin{array}{c c} I _ {n} & - \lambda I _ {n} \\ 0 & I _ {n} \end{array} \right] \left[ \begin{array}{c c} A _ {0} + \lambda P _ {1} (\lambda) & 0 \\ 0 & I _ {n} \end{array} \right] \left[ \begin{array}{c c} 0 & I _ {n} \\ - I _ {n} & P _ {1} (\lambda) \end{array} \right] = \left[ \begin{array}{c c} \lambda I _ {n} & A _ {0} \\ - I _ {n} & P _ {1} (\lambda) \end{array} \right].
$$

Notice that the transformation matrices have unit determinant and that the λ-matrix on the right-hand side has degree d − 1. The process can be repeated. If

$$
P _ {2} (\lambda) = A _ {2} + \lambda A _ {3} + \dots + \lambda^ {d - 2} A _ {d}
$$

then

$$
P _ {1} (\lambda) = A _ {1} + \lambda P _ {2} (\lambda)
$$

and

$$
\left[ \begin{array}{c c c} I _ {n} & 0 & 0 \\ \hline 0 & I _ {n} & - \lambda I _ {n} \\ 0 & 0 & I _ {n} \end{array} \right] \left[ \begin{array}{c c c} \lambda I _ {n} & A _ {0} & 0 \\ \hline - I _ {n} & P _ {1} (\lambda) & 0 \\ 0 & 0 & I _ {n} \end{array} \right] \left[ \begin{array}{c c c} I _ {n} & 0 & 0 \\ \hline 0 & 0 & I _ {n} \\ 0 & - I _ {n} & P _ {2} (\lambda) \end{array} \right] =
$$

$$
\left[ \begin{array}{c c c} \lambda I _ {n} & 0 & A _ {0} \\ \hline - I _ {n} & \lambda I _ {n} & A _ {1} \\ 0 & - I _ {n} & P _ {2} (\lambda) \end{array} \right].
$$

Note that the matrix on the right has degree d − 2. A straightforward induction argument can be assembled to establish that if the dn-by-dn matrices $S ( \lambda )$ and $T ( \lambda )$ are defined by

$$
S (\lambda) = \left[ \begin{array}{c c c c c} I _ {n} & - \lambda I _ {n} & 0 & \dots & 0 \\ 0 & I _ {n} & - \lambda I _ {n} & & \vdots \\ 0 & & \ddots & \ddots & \\ \vdots & & & I _ {n} & - \lambda I _ {n} \\ 0 & 0 & \dots & 0 & I _ {n} \end{array} \right], \quad T (\lambda) = \left[ \begin{array}{c c c c c} 0 & 0 & 0 & \dots & I \\ - I _ {n} & 0 & & & P _ {1} (\lambda) \\ 0 & - I _ {n} & \ddots & & \vdots \\ \vdots & & \ddots & \ddots & P _ {d - 2} (\lambda) \\ 0 & 0 & \dots & - I _ {n} & P _ {d - 1} (\lambda) \end{array} \right]
$$

where

$$
P _ {k} (\lambda) = A _ {k} + \lambda A _ {k + 1} + \dots + \lambda^ {d - k} A _ {d},
$$

then

$$
S (\lambda) \left[ \begin{array}{c c} P (\lambda) & 0 \\ 0 & I _ {(d - 1) n} \end{array} \right] T (\lambda) = \left[ \begin{array}{c c c c c} \lambda I _ {n} & 0 & 0 & \dots & A _ {0} \\ - I _ {n} & \lambda I _ {n} & & & A _ {1} \\ 0 & - I _ {n} & \ddots & & \vdots \\ \vdots & & \ddots & \lambda I _ {n} & A _ {d - 2} \\ 0 & 0 & \dots & - I _ {n} & A _ {d - 1} + \lambda A _ {d} \end{array} \right].
$$

Note that, if we solve the linearized problem using the QZ algorithm, then $O ( ( d n ) ^ { 3 } )$ flops are required.

# Problems

P7.7.1 Suppose A and B are in $\mathbb { R } ^ { n \times n }$ and that

$$
U ^ {T} B V = \left[ \begin{array}{c c} D & 0 \\ 0 & 0 \\ r & n - r \end{array} \right] _ {n - r} ^ {r}, \quad U = \left[ \begin{array}{c c} U _ {1} & U _ {2} \\ r & n - r \end{array} \right], \quad V = \left[ \begin{array}{c c} V _ {1} & V _ {2} \\ r & n - r \end{array} \right],
$$

is the SVD of B, where D is r-by-r and $r = { \mathsf { r a n k } } ( B )$ . Show that if $\lambda ( A , B ) = \mathbb { C }$ then $U _ { 2 } ^ { T } A V _ { 2 }$ is singular.

P7.7.2 Suppose A and B are in $\mathbb { R } ^ { n \times n }$ . Give an algorithm for computing orthogonal Q and Z such that $Q ^ { T } A \dot { Z }$ is upper Hessenberg and $Z ^ { T } B Q$ is upper triangular.

# P7.7.3 Suppose

$$
A = \left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ 0 & A _ {2 2} \end{array} \right] \quad \text {and} \quad B = \left[ \begin{array}{c c} B _ {1 1} & B _ {1 2} \\ 0 & B _ {2 2} \end{array} \right]
$$

with $A _ { 1 1 } , B _ { 1 1 } \in \mathbb { R } ^ { k \times k }$ and $A _ { 2 2 } , B _ { 2 2 } \in \mathbb { R } ^ { j \times j }$ . Under what circumstances do there exist

$$
X = \left[ \begin{array}{c c} I _ {k} & X _ {1 2} \\ 0 & I _ {j} \end{array} \right] \quad \text {and} \quad Y = \left[ \begin{array}{c c} I _ {k} & Y _ {1 2} \\ 0 & I _ {j} \end{array} \right]
$$

so that $Y ^ { - 1 } A X$ and $Y ^ { - 1 } B X$ are both block diagonal? This is the generalized Sylvester equation problem. Specify an algorithm for the case when $A _ { 1 1 } , A _ { 2 2 } , B _ { 1 1 }$ , and $B _ { 2 2 }$ are upper triangular. See K˚agstr¨om (1994).

P7.7.4 Suppose $\mu \not \in \lambda ( A , B )$ . Relate the eigenvalues and eigenvectors of $A _ { 1 } = ( A - \mu B ) ^ { - 1 } A$ and $B _ { 1 } = ( A - \mu B ) ^ { - 1 } B$ to the generalized eigenvalues and eigenvectors of $A - \lambda B$ .

P7.7.5 What does the generalized Schur decomposition say about the pencil $A - \lambda A ^ { T } ?$ Hint: If $T \in \mathbb { R } ^ { n \times n }$ is upper triangular, then $\mathcal { E } _ { n } T \mathcal { E } _ { n }$ is lower triangular where $\mathcal { E } _ { n }$ is the exchange permutation defined in §1.2.11.

# P7.7.6 Prove that

$$
L _ {1} (\lambda) = \left[ \begin{array}{c c c c} A _ {3} + \lambda A _ {4} & A _ {2} & A _ {1} & A _ {0} \\ - I _ {n} & 0 & 0 & 0 \\ 0 & - I _ {n} & 0 & 0 \\ 0 & 0 & - I _ {n} & 0 \end{array} \right], L _ {2} (\lambda) = \left[ \begin{array}{c c c c} A _ {3} + \lambda A _ {4} & - I _ {n} & 0 & 0 \\ A _ {2} & 0 & - I _ {n} & 0 \\ A _ {1} & 0 & 0 & - I _ {n} \\ A _ {0} & 0 & 0 & 0 \end{array} \right]
$$

are linearizations of

$$
P (\lambda) = A _ {0} + \lambda A _ {1} + \lambda^ {2} A _ {2} + \lambda^ {3} A _ {3} + \lambda^ {4} A _ {4}.
$$

Specify the λ-matrix transformations that relate diag $( P ( \lambda ) , I _ { 3 n } )$ to both $L _ { 1 } ( \lambda )$ and $L _ { 2 } ( \lambda )$ .

# Notes and References for §7.7

For background to the generalized eigenvalue problem we recommend Stewart(IMC), Stewart and Sun (MPT), and Watkins (MEP) and:

B. K˚agstr¨om and A. Ruhe (1983). Matrix Pencils, Proceedings Pite Havsbad, 1982, Lecture Notes in Mathematics Vol. 973, Springer-Verlag, New York.

QZ-related papers include:

C.B. Moler and G.W. Stewart (1973). “An Algorithm for Generalized Matrix Eigenvalue Problems,” SIAM J. Numer. Anal. 10, 241–256.

L. Kaufman (1974). “The LZ Algorithm to Solve the Generalized Eigenvalue Problem,” SIAM J. Numer. Anal. 11, 997–1024.

R.C. Ward (1975). “The Combination Shift QZ Algorithm,” SIAM J. Numer. Anal. 12, 835–853.

C.F. Van Loan (1975). “A General Matrix Eigenvalue Algorithm,” SIAM J. Numer. Anal. 12, 819–834.

L. Kaufman (1977). “Some Thoughts on the QZ Algorithm for Solving the Generalized Eigenvalue Problem,” ACM Trans. Math. Softw. 3, 65–75.

R.C. Ward (1981). “Balancing the Generalized Eigenvalue Problem,” SIAM J. Sci. Stat. Comput. 2, 141–152.

P. Van Dooren (1982). “Algorithm 590: DSUBSP and EXCHQZ: Fortran Routines for Computing Deflating Subspaces with Specified Spectrum,” ACM Trans. Math. Softw. 8, 376–382.

K. Dackland and B. K˚agstr¨om (1999). “Blocked Algorithms and Software for Reduction of a Regular Matrix Pair to Generalized Schur Form,” ACM Trans. Math. Softw. 25, 425–454.

D.S. Watkins (2000). “Performance of the QZ Algorithm in the Presence of Infinite Eigenvalues,” SIAM J. Matrix Anal. Applic. 22, 364–375.

B. K˚agstr¨om, D. Kressner, E.S. Quintana-Orti, and G. Quintana-Orti (2008). “Blocked Algorithms for the Reduction to Hessenberg-Triangular Form Revisited,” BIT 48, 563–584.

Many algorithmic ideas associated with the A − λI problem extend to the A − λB problem:

A. Jennings and M.R. Osborne (1977). “Generalized Eigenvalue Problems for Certain Unsymmetric Band Matrices,” Lin. Alg. Applic. 29, 139–150.

V.N. Kublanovskaya (1984). “AB Algorithm and Its Modifications for the Spectral Problem of Linear Pencils of Matrices,” Numer. Math. 43, 329–342.   
Z. Bai, J. Demmel, and M. Gu (1997). “An Inverse Free Parallel Spectral Divide and Conquer Algorithm for Nonsymmetric Eigenproblems,” Numer. Math. 76, 279–308.   
G.H. Golub and Q. Ye (2000). “Inexact Inverse Iteration for Generalized Eigenvalue Problems,” BIT 40, 671–684.   
F. Tisseur (2001). “Newton’s Method in Floating Point Arithmetic and Iterative Refinement of Generalized Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 22, 1038–1057.   
D. Lemonnier and P. Van Dooren (2006). “Balancing Regular Matrix Pencils,” SIAM J. Matrix Anal. Applic. 28, 253–263.   
R. Granat, B. K˚agstr¨om, and D. Kressner (2007). “Computing Periodic Deflating Subspaces Associated with a Specified Set of Eigenvalues,” BIT 47, 763–791.   
The perturbation theory for the generalized eigenvalue problem is treated in:   
G.W. Stewart (1972). “On the Sensitivity of the Eigenvalue Problem Ax = λBx,” SIAM J. Numer. Anal. 9, 669–686.   
G.W. Stewart (1973). “Error and Perturbation Bounds for Subspaces Associated with Certain Eigenvalue Problems,” SIAM Review 15, 727–764.   
G.W. Stewart (1975). “Gershgorin Theory for the Generalized Eigenvalue Problem Ax = λBx,” Math. Comput. 29, 600–606.   
A. Pokrzywa (1986). “On Perturbations and the Equivalence Orbit of a Matrix Pencil,” Lin. Alg. Applic. 82, 99–121.   
J. Sun (1995). “Perturbation Bounds for the Generalized Schur Decomposition,” SIAM J. Matrix Anal. Applic. 16, 1328–1340.   
R. Bhatia and R.–C. Li (1996). “On Perturbations of Matrix Pencils with Real Spectra. II,” Math. Comput. 65, 637–645.   
J.-P. Dedieu (1997). “Condition Operators, Condition Numbers, and Condition Number Theorem for the Generalized Eigenvalue Problem,” Lin. Alg. Applic. 263, 1–24.   
D.J. Higham and N.J. Higham (1998). “Structured Backward Error and Condition of Generalized Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 20, 493–512.   
R. Byers, C. He, and V. Mehrmann (1998). “Where is the Nearest Non-Regular Pencil?,” Lin. Alg. Applic. 285, 81–105.   
V. Frayss and V. Toumazou (1998). “A Note on the Normwise Perturbation Theory for the Regular Generalized Eigenproblem,” Numer. Lin. Alg. 5, 1–10.   
R.–C. Li (2003). “On Perturbations of Matrix Pencils with Real Spectra, A Revisit,” Math. Comput. 72, 715–728.   
S. Bora and V. Mehrmann (2006). “Linear Perturbation Theory for Structured Matrix Pencils Arising in Control Theory,” SIAM J. Matrix Anal. Applic. 28, 148–169.   
X.S. Chen (2007). “On Perturbation Bounds of Generalized Eigenvalues for Diagonalizable Pairs,” Numer. Math. 107, 79–86.   
The Kronecker structure of the pencil A − λB is analogous to Jordan structure of A − λI and it can provide useful information about the underlying application. Papers concerned with this important decomposition include:   
J.H. Wilkinson (1978). “Linear Differential Equations and Kronecker’s Canonical Form,” in Recent Advances in Numerical Analysis, C. de Boor and G.H. Golub (eds.), Academic Press, New York, 231–265.   
J.H. Wilkinson (1979). “Kronecker’s Canonical Form and the QZ Algorithm,” Lin. Alg. Applic. 28, 285–303.   
P. Van Dooren (1979). “The Computation of Kronecker’s Canonical Form of a Singular Pencil,” Lin. Alg. Applic. 27, 103–140.   
J.W. Demmel (1983). “The Condition Number of Equivalence Transformations that Block Diagonalize Matrix Pencils,” SIAM J. Numer. Anal. 20, 599–610.   
J.W. Demmel and B. K˚agstr¨om (1987). “Computing Stable Eigendecompositions of Matrix Pencils,” Linear Alg. Applic. 88/89, 139–186.   
B. K˚agstr¨om (1985). “The Generalized Singular Value Decomposition and the General A − λB Problem,” BIT 24, 568–583.   
B. K˚agstr¨om (1986). “RGSVD: An Algorithm for Computing the Kronecker Structure and Reducing Subspaces of Singular A λB Pencils,” SIAM J. Sci. Stat. Comput. 7, 185–211.

J. Demmel and B. K˚agstr¨om (1986). “Stably Computing the Kronecker Structure and Reducing Subspaces of Singular Pencils A  λB for Uncertain Data,” in Large Scale Eigenvalue Problems, J. Cullum and R.A. Willoughby (eds.), North-Holland, Amsterdam.   
T. Beelen and P. Van Dooren (1988). “An Improved Algorithm for the Computation of Kronecker’s Canonical Form of a Singular Pencil,” Lin. Alg. Applic. 105, 9–65.   
E. Elmroth and B. K˚agstr¨om(1996). “The Set of 2-by-3 Matrix Pencils — Kronecker Structures and Their Transitions under Perturbations,” SIAM J. Matrix Anal. Applic. 17, 1–34.   
A. Edelman, E. Elmroth, and B. K˚agstr¨om (1997). “A Geometric Approach to Perturbation Theory of Matrices and Matrix Pencils Part I: Versal Defformations,” SIAM J. Matrix Anal. Applic. 18, 653–692.   
E. Elmroth, P. Johansson, and B. K˚agstr¨om (2001). “Computation and Presentation of Graphs Displaying Closure Hierarchies of Jordan and Kronecker Structures,” Num. Lin. Alg. 8, 381–399.   
Just as the Schur decomposition can be used to solve the Sylvester equation problem $A _ { 1 } X - X A _ { 2 } = B _ { 1 }$ , the generalized Schur decomposition can be used to solve the generalized Sylvester equation problem where matrices X and Y are sought so that $A _ { 1 } X - Y A _ { 2 } = B _ { 1 }$ and $A _ { 3 } X - Y A _ { 4 } = B _ { 2 }$ , see:   
W. Enright and S. Serbin (1978). “A Note on the Efficient Solution of Matrix Pencil Systems,” BIT 18, 276–81.   
B. K˚agstr¨om and L. Westin (1989). “Generalized Schur Methods with Condition Estimators for Solving the Generalized Sylvester Equation,” IEEE Trans. Autom. Contr. AC-34, 745–751.   
B. K˚agstr¨om (1994). “A Perturbation Analysis of the Generalized Sylvester Equation (AR−LB, DR− LE) = (C, F ),” SIAM J. Matrix Anal. Applic. 15, 1045–1060.   
J.-G. Sun (1996). “Perturbation Analysis of System Hessenberg and Hessenberg-Triangular Forms,” Lin. Alg. Applic. 241–3, 811–849.   
B. K˚agstr¨om and P. Poromaa (1996). “LAPACK-style Algorithms and Software for Solving the Generalized Sylvester Equation and Estimating the Separation Between Regular Matrix Pairs,” ACM Trans. Math. Softw. 22, 78–103.   
I. Jonsson and B. K˚agstr¨om (2002). “Recursive Blocked Algorithms for Solving Triangular Systems– Part II: Two-sided and Generalized Sylvester and Lyapunov Matrix Equations,” ACM Trans. Math. Softw. 28, 416–435.   
R. Granat and B. K˚agstr¨om (2010). “Parallel Solvers for Sylvester-Type Matrix Equations with Applications in Condition Estimation, Part I: Theory and Algorithms,” ACM Trans. Math. Softw. 37, Article 32.   
Rectangular generalized eigenvalue problems also arise. In this setting the goal is to reduce the rank of A − λB, see:   
G.W. Stewart (1994). “Perturbation Theory for Rectangular Matrix Pencils,” Lin. Alg. Applic. 208/209, 297–301.   
G. Boutry, M. Elad, G.H. Golub, and P. Milanfar (2005). “The Generalized Eigenvalue Problem for Nonsquare Pencils Using a Minimal Perturbation Approach,” SIAM J. Matrix Anal. Applic. 27, 582–601.   
D. Chu and G.H. Golub (2006). “On a Generalized Eigenvalue Problem for Nonsquare Pencils,” SIAM J. Matrix Anal. Applic. 28, 770–787.   
References for the polynomial eigenvalue problem include:   
P. Lancaster (1966). Lambda-Matrices and Vibrating Systems, Pergamon Press, Oxford, U.K.   
I. Gohberg, P. Lancaster, and L. Rodman (1982). Matrix Polynomials, Academic Press, New York.   
F. Tisseur (2000). “Backward Error and Condition of Polynomial Eigenvalue Problems,” Lin. Alg. Applic. 309, 339–361.   
J.-P. Dedieu and F. Tisseur (2003). “Perturbation Theory for Homogeneous Polynomial Eigenvalue Problems,” Lin. Alg. Applic. 358, 71–94.   
N.J. Higham, D.S. Mackey, and F. Tisseur (2006). “The Conditioning of Linearizations of Matrix Polynomials,” SIAM J. Matrix Anal. Applic. 28, 1005–1028.   
D.S. Mackey, N. Mackey, C. Mehl, V. Mehrmann (2006). “Vector Spaces of Linearizations for Matrix Polynomials,” SIAM J. Matrix Anal. Applic. 28, 971–1004.   
The structured quadratic eigenvalue problem is discussed briefly in §8.7.9.

# 7.8 Hamiltonian and Product Eigenvalue Problems

Two structured unsymmetric eigenvalue problems are considered. The Hamiltonian matrix eigenvalue problem comes with its own special Schur decomposition. Orthogonal symplectic similarity transformations are used to bring about the required reduction. The product eigenvalue problem involves computing the eigenvalues of a product like $A _ { 1 } A _ { 2 } ^ { - 1 } A _ { 3 }$ without actually forming the product or the designated inverses. For detailed background to these problems, see Kressner (NMGS) and Watkins (MEP).

# 7.8.1 Hamiltonian Matrix Eigenproblems

Hamiltonian and symplectic matrices are introduced in §1.3.10. Their 2-by-2 block structure provide a nice framework for practicing block matrix manipulation, see P1.3.2 and P2.5.4. We now describe some interesting eigenvalue problems that involve these matrices. For a given n, we define the matrix $\check { J } \in \mathbb { R } ^ { 2 n \times 2 n }$ by

$$
J = \left[ \begin{array}{c c} 0 & I _ {n} \\ - I _ {n} & 0 \end{array} \right]
$$

and proceed to work with the families of 2-by-2 block structured matrices that are displayed in Figure 7.8.1. We mention four important facts concerning these matrices.

<table><tr><td>Family</td><td>Definition</td><td colspan="2">What They Look Like</td></tr><tr><td>Hamiltonian</td><td> $JM = (JM)^T$ </td><td> $M = \begin{bmatrix} A & G \\ F & -A^T \end{bmatrix}$ </td><td> $G$  symmetric $F$  symmetric</td></tr><tr><td>SkewHamiltonian</td><td> $JN = -(JN)^T$ </td><td> $N = \begin{bmatrix} A & G \\ F & A^T \end{bmatrix}$ </td><td> $G$  skew-symmetric $F$  skew-symmetric</td></tr><tr><td>Symplectic</td><td> $JS = S^{-T}J$ </td><td> $S = \begin{bmatrix} S_{11} & S_{12} \\ S_{21} & S_{22} \end{bmatrix}$ </td><td> $S_{11}^T S_{21}$  symmetric $S_{22}^T S_{12}$  symmetric $S_{11}^T S_{22} = I + S_{21}^T S_{12}$ </td></tr><tr><td>OrthogonalSymplectic</td><td> $JQ = QJ$ </td><td> $Q = \begin{bmatrix} Q_1 & Q_2 \\ -Q_2 & Q_1 \end{bmatrix}$ </td><td> $Q_1^T Q_2$  symmetric $I = Q_1^T Q_1 + Q_2^T Q_2$ </td></tr></table>

Figure 7.8.1. Hamiltonian and symplectic structures

(1) Symplectic similarity transformations preserve Hamiltonian structure:

$$
J (S ^ {- 1} M S) = (J S ^ {- 1} J ^ {T}) (J M J ^ {T}) (J S) = - S ^ {T} M ^ {T} S ^ {- T} J = (J (S ^ {- 1} M S)) ^ {T}.
$$

(2) The square of a Hamiltonian matrix is skew-Hamiltonian:

$$
J M ^ {2} = (J M J ^ {T}) (J M) = - M ^ {T} (J M) ^ {T} = - M ^ {2 T} J ^ {T} = - (J M ^ {2}) ^ {T}.
$$

(3) If M is a Hamiltonian matrix and $\lambda \in \lambda ( M )$ , then $- \lambda \in \lambda ( M )$ :

$$
M \left[ \begin{array}{c} u \\ v \end{array} \right] = \lambda \left[ \begin{array}{c} u \\ v \end{array} \right] \qquad \Rightarrow \qquad M ^ {T} \left[ \begin{array}{c} v \\ - u \end{array} \right] = - \lambda \left[ \begin{array}{c} v \\ - u \end{array} \right].
$$

(4) If S is symplectic and $\lambda \in \lambda ( S )$ , then $1 / \lambda \in \lambda ( S )$ :

$$
S \left[ \begin{array}{l} u \\ v \end{array} \right] = \lambda \left[ \begin{array}{l} u \\ v \end{array} \right] \quad \Rightarrow \quad S ^ {T} \left[ \begin{array}{l} v \\ - u \end{array} \right] = \frac {1}{\lambda} \left[ \begin{array}{l} v \\ - u \end{array} \right].
$$

Symplectic versions of Householder and Givens transformations have a prominanent role to play in Hamiltonian matrix computations. If $P = I _ { n } - 2 v v ^ { T }$ is a Householder matrix, then diag(P, P ) is a symplectic orthogonal matrix. Likewise, if $G \in \mathbb { R } ^ { 2 n \times 2 n }$ is a Givens rotation that involves planes i and $i + n$ , then G is a symplectic orthogonal matrix. Combinations of these transformations can be used to introduce zeros. For example, a Householder-Givens-Householder sequence can do this:

$$
\left[ \begin{array}{c} \times \\ \times \\ \times \\ \times \\ \hline \times \\ \times \\ \times \\ \times \\ \end{array} \right] \quad \stackrel {{\mathrm{diag} (P _ {1}, P _ {1})}} {{\longrightarrow}} \quad \left[ \begin{array}{c} \times \\ \times \\ \times \\ \times \\ \hline \times \\ 0 \\ 0 \\ 0 \end{array} \right] \quad \stackrel {{G _ {1, 5}}} {{\longrightarrow}} \quad \left[ \begin{array}{c} \times \\ \times \\ \times \\ \times \\ \hline 0 \\ 0 \\ 0 \\ 0 \end{array} \right] \quad \stackrel {{\mathrm{diag} (P _ {2}, P _ {2})}} {{\longrightarrow}} \quad \left[ \begin{array}{c} \times \\ 0 \\ 0 \\ 0 \\ \hline 0 \\ 0 \\ 0 \\ 0 \end{array} \right].
$$

This kind of vector reduction can be sequenced to produce a constructive proof of a structured Schur decomposition for Hamiltonian matrices. Suppose λ is a real eigenvalue of a Hamiltonian matrix M and that $x \in \mathbb { R } ^ { 2 n }$ is a unit 2-norm vector with $M x = \lambda x$ . If $Q _ { 1 } \in \mathbb { R } ^ { 2 n \times 2 n }$ is an orthogonal symplectic matrix and $Q _ { 1 } ^ { T } x = e _ { 1 }$ , then it follows from $( Q _ { 1 } ^ { T } M Q _ { 1 } ) ( Q _ { 1 } ^ { T } x ) = \lambda ( Q _ { 1 } ^ { T } x )$ that

$$
Q _ {1} ^ {T} M Q _ {1} = \left[ \begin{array}{c c c c c c c c} \lambda & \times & \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times & \times & \times \\ \hline 0 & 0 & 0 & 0 & - \lambda & 0 & 0 & 0 \\ 0 & \times & \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times & \times & \times \end{array} \right].
$$

The “extra” zeros follow from the Hamiltonian structure of $Q _ { 1 } ^ { T } M Q _ { 1 }$ . The process can be repeated on the 6-by-6 Hamiltonian submatrix defined by rows and columns 2-3-4- 6-7-8. Together with the assumption that M has no purely imaginary eigenvalues, it is possible to show that an orthogonal symplectic matrix Q exists so that

$$
Q ^ {T} M Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ - Q _ {2} & Q _ {1} \end{array} \right] ^ {T} \left[ \begin{array}{c c} A & F \\ G & - A ^ {T} \end{array} \right] \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ - Q _ {2} & Q _ {1} \end{array} \right] = \left[ \begin{array}{c c} T & R \\ 0 & - T ^ {T} \end{array} \right] \tag {7.8.1}
$$

where $T \in \mathbb { R } ^ { n \times n }$ is upper quasi-triangular. This is the real Hamiltonian-Schur decomposition. See Paige and Van Loan (1981) and, for a more general version, Lin, Mehrmann, and Xu (1999).

One reason that the Hamiltonian eigenvalue problem is so important is its connection to the algebraic Ricatti equation

$$
G + X A + A ^ {T} X - X F X = 0. \tag {7.8.2}
$$

This quadratic matrix problem arises in optimal control and a symmetric solution is sought so that the eigenvalues of $A - F X$ are in the open left half plane. Modest assumptions typically ensure that M has no eigenvalues on the imaginary axis and that the matrix $Q _ { 1 }$ in (7.8.1) is nonsingular. If we compare (2,1) blocks in (7.8.1), then

$$
Q _ {2} ^ {T} A Q _ {1} - Q _ {2} ^ {T} F Q _ {2} + Q _ {1} ^ {T} G Q _ {1} + Q _ {1} ^ {T} A ^ {T} Q _ {2} = 0.
$$

It follows from $I _ { n } = Q _ { 1 } ^ { T } Q _ { 1 } + Q _ { 2 } ^ { T } Q _ { 2 }$ that $X = Q _ { 2 } Q _ { 1 } ^ { - 1 }$ is symmetric and that it satisfies (7.8.2). From (7.8.1) it is easy to show that $A { \bf - } ^ { \prime } F X = Q _ { 1 } T Q _ { 1 } ^ { - 1 }$ and so the eigenvalues of $A - F X$ are the eigenvalues of T . It follows that the desired solution to the algebraic Ricatti equation can be obtained by computing the real Hamiltonian-Schur decomposition and ordering the eigenvalues so that $\lambda ( T )$ is in the left half plane.

How might the real Hamiltonian-Schur form be computed? One idea is to reduce $M$ to some condensed Hamiltonian form and then devise a structure-preserving QRiteration. Regarding the former task, it is easy to compute an orthogonal symplectic $U _ { 0 }$ so that

$$
U _ {0} ^ {T} M U _ {0} = \left[ \begin{array}{c c} H & R \\ D & - H ^ {T} \end{array} \right] \tag {7.8.3}
$$

where $H \in \mathbb { R } ^ { n \times n }$ is upper Hessenberg and D is diagonal. Unfortunately, a structurepreserving QR iteration that maintains this condensed form has yet to be devised. This impasse prompts consideration of methods that involve the skew-Hamiltonian matrix $N = M ^ { 2 }$ . Because the (2,1) block of a skew-Hamiltonian matrix is skew-symmetric, it has a zero diagonal. Symplectic similarity transforms preserve skew-Hamiltonian structure, and it is straightforward to compute an orthogonal symplectic matrix $V _ { 0 }$ such that

$$
V _ {0} ^ {T} M ^ {2} V _ {0} = \left[ \begin{array}{c c} H & R \\ 0 & H ^ {T} \end{array} \right], \tag {7.8.4}
$$

where H is upper Hessenberg. If $U ^ { T } H U = T$ is the real Schur form of H and and $Q = V _ { 0 } \cdot \mathrm { d i a g } ( U , U )$ , then

$$
Q ^ {T} M ^ {2} Q = \left[ \begin{array}{c c} T & U ^ {T} R U \\ 0 & T ^ {T} \end{array} \right]
$$

is the real skew-Hamiltonian Schur form. See Van Loan (1984). It does not follow that $Q ^ { T } M Q$ is in Schur-Hamiltonian form. Moreover, the quality of the computed small eigenvalues is not good because of the explicit squaring of M. However, these shortfalls can be overcome in an efficient numerically sound way, see Chu, Lie, and Mehrmann (2007) and the references therein. Kressner (NMSE, p. 175–208) and Watkins (MEP, p. 319–341) have in-depth treatments of the Hamiltonian eigenvalue problem.

# 7.8.2 Product Eigenvalue Problems

Using SVD and QZ, we can compute the eigenvalues of $A ^ { T } A$ and $B ^ { - 1 } A$ without forming products or inverses. The intelligent computation of the Hamiltonian-Schur decomposition involves a correspondingly careful handling of the product M -times-M . In this subsection we further develop this theme by discussing various product decompositions. Here is an example that suggests how we might compute the Hessenberg decomposition of

$$
A = A _ {3} A _ {2} A _ {1}
$$

where $A _ { 1 } , A _ { 2 } , A _ { 3 } \in \mathbb { R } ^ { n \times n }$ . Instead of forming this product explicitly, we compute orthogonal $U _ { 1 } , U _ { 2 } , U _ { 3 } \in \mathbb { R } ^ { n \times n }$ such that

$$
U _ {1} ^ {T} A _ {3} U _ {3} = H _ {3} \quad \text {(upper Hessenberg)},
$$

$$
U _ {3} ^ {T} A _ {2} U _ {2} = T _ {2} \quad (\text { upper   triangular }), \tag {7.8.5}
$$

$$
U _ {2} ^ {T} A _ {1} U _ {1} = T _ {1} \quad \text {(upper triangular)}.
$$

It follows that

$$
U _ {1} ^ {T} A U _ {1} = (U _ {1} ^ {T} A _ {3} U _ {3}) (U _ {3} ^ {T} A _ {2} U _ {2}) (U _ {2} ^ {T} A _ {1} U _ {1}) = H _ {3} T _ {2} T _ {1}
$$

is upper Hessenberg. A procedure for doing this would start by computing the QR factorizations

$$
Q _ {2} ^ {T} A _ {1} = R _ {1}, \quad Q _ {3} ^ {T} (A _ {2} Q _ {2}) = R _ {2}.
$$

If $\tilde { A } _ { 3 } = A _ { 3 } Q _ { 3 }$ , then $A = \tilde { A } _ { 3 } R _ { 2 } R _ { 1 }$ . The next phase involves reducing ${ \tilde { A } } _ { 3 }$ to Hessenberg form with Givens transformations coupled with “bulge chasing” to preserve the triangular structures already obtained. The process is similar to the reduction of $A - \lambda B$ to Hessenberg-triangular form; see §7.7.4.

Now suppose we want to compute the real Schur form of A

$$
Q _ {1} ^ {T} A _ {3} Q _ {3} = T _ {3} \quad \text {(upper quasi - triangular)},
$$

$$
Q _ {3} ^ {T} A _ {2} Q _ {2} = T _ {2} \quad \text {(upper triangular)}, \tag {7.8.6}
$$

$$
Q _ {2} ^ {T} A _ {1} Q _ {1} = T _ {1} \quad \text {(upper triangular)},
$$

where $Q _ { 1 } , Q _ { 2 } , Q _ { 3 } \in \mathbb { R } ^ { n \times n }$ are orthogonal. Without loss of generality we may assume that $\{ A _ { 3 } , A _ { 2 } , A _ { 1 } \}$ is in Hessenberg-triangular-triangular form. Analogous to the QZ iteration, the next phase is to produce a sequence of converging triplets

$$
\{A _ {3} ^ {(k)}, A _ {2} ^ {(k)}, A _ {1} ^ {(k)} \} \rightarrow \{T _ {3}, T _ {2}, T _ {1} \} \tag {7.8.7}
$$

with the property that all the iterates are in Hessenberg-triangular-triangular form.

Product decompositions (7.8.5) and (7.8.6) can be framed as structured decompositions of block-cyclic 3-by-3 matrices. For example, if

$$
U = \left[ \begin{array}{c c c} U _ {1} & 0 & 0 \\ 0 & U _ {2} & 0 \\ 0 & 0 & U _ {3} \end{array} \right]
$$

then we have the following restatement of (7.8.5):

$$
U ^ {T} \left[ \begin{array}{c c c} 0 & 0 & A _ {3} \\ A _ {1} & 0 & 0 \\ 0 & A _ {2} & 0 \end{array} \right] U = \left[ \begin{array}{c c c} 0 & 0 & H _ {3} \\ T _ {1} & 0 & 0 \\ 0 & T _ {2} & 0 \end{array} \right] = \tilde {H}.
$$

Consider the zero-nonzero structure of this matrix for the case $n = 4$ :

$$
\tilde {H} = \left[ \begin{array}{c c c c c c c c c c c c} 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times \\ \hline \times & \times & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & \times & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \hline 0 & 0 & 0 & 0 & \times & \times & \times & \times & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & \times & \times & \times & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & 0 & 0 \end{array} \right].
$$

Using the perfect shuffle $\mathcal { P } _ { 3 4 }$ (see §1.2.11) we also have

$$
\mathcal {P} _ {3 4} \tilde {H} \mathcal {P} _ {3 4} = \left[ \begin{array}{c c c c c c c c c c c c} 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times \\ \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 \\ 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 \\ \hline 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times \\ 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 \\ \hline 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 \end{array} \right].
$$

Note that this is a highly structured 12-by-12 upper Hessenberg matrix. This connection makes it possible to regard the product-QR iteration as a structure-preserving

QR iteration. For a detailed discussion about this connection and its implications for both analysis and computation, see Kressner (NMSE, pp. 146–174) and Watkins(MEP, pp. 293–303). We mention that with the “technology” that has been developed, it is possible to solve product eigenvalue problems where the factor matrices that define A are rectangular. Square nonsingular factors can also participate through their inverses, $\mathrm { e . g . , } A = A _ { 3 } A _ { 2 } ^ { - 1 } A _ { 1 }$ .

# Problems

P7.8.1 What can you say about the eigenvalues and eigenvectors of a symplectic matrix?

P7.8.2 Suppose $S _ { 1 } , S _ { 2 } \in \mathbb { R } ^ { n \times n }$ are both skew-symmetric and let $A = S _ { 1 } S _ { 2 }$ . Show that the nonzero eigenvalues of A are not simple. How would you compute these eigenvalues?

P7.8.3 Relate the eigenvalues and eigenvectors of

$$
A = \left[ \begin{array}{c c c c} 0 & A _ {1} & 0 & 0 \\ 0 & 0 & A _ {2} & 0 \\ 0 & 0 & 0 & A _ {3} \\ A _ {4} & 0 & 0 & 0 \end{array} \right].
$$

to the eigenvalues and eigenvectors of $\tilde { A } = A _ { 1 } A _ { 2 } A _ { 3 } A _ { 4 }$ . Assume that the diagonal blocks are square.

# Notes and References for §7.8

The books by Kressner(NMSE) and Watkins (MEP) have chapters on product eigenvalue problems and Hamiltonian eigenvalue problems. The sometimes bewildering network of interconnections that exist among various structured classes of matrices is clarified in:

A. Bunse-Gerstner, R. Byers, and V. Mehrmann (1992). “A Chart of Numerical Methods for Structured Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 13, 419–453.

Papers concerned with the Hamiltonian Schur decomposition include:

A.J. Laub and K. Meyer (1974). “Canonical Forms for Symplectic and Hamiltonian Matrices,” J. Celestial Mechanics 9, 213–238.

C.C. Paige and C. Van Loan (1981). “A Schur Decomposition for Hamiltonian Matrices,” Lin. Alg. Applic. 41, 11–32.

V. Mehrmann (1991). Autonomous Linear Quadratic Control Problems, Theory and Numerical Solution, Lecture Notes in Control and Information Sciences No. 163, Springer-Verlag, Heidelberg.

W.-W. Lin, V. Mehrmann, and H. Xu (1999). “Canonical Forms for Hamiltonian and Symplectic Matrices and Pencils,” Lin. Alg. Applic. 302/303, 469–533.

Various methods for Hamiltonian eigenvalue problems have been devised that exploit the rich underlying structure, see:

C. Van Loan (1984). “A Symplectic Method for Approximating All the Eigenvalues of a Hamiltonian Matrix,” Lin. Alg. Applic. 61, 233–252.

R. Byers (1986) “A Hamiltonian QR Algorithm,” SIAM J. Sci. Stat. Comput. 7, 212–229.

P. Benner, R. Byers, and E. Barth (2000). “Algorithm 800: Fortran 77 Subroutines for Computing the Eigenvalues of Hamiltonian Matrices. I: the Square-Reduced Method,” ACM Trans. Math. Softw. 26, 49–77.

H. Fassbender, D.S. Mackey and N. Mackey (2001). “Hamilton and Jacobi Come Full Circle: Jacobi Algorithms for Structured Hamiltonian Eigenproblems,” Lin. Alg. Applic. 332-4, 37–80.

D.S. Watkins (2006). “On the Reduction of a Hamiltonian Matrix to Hamiltonian Schur Form,” ETNA 23, 141–157.

D.S. Watkins (2004). “On Hamiltonian and Symplectic Lanczos Processes,” Lin. Alg. Applic. 385, 23–45.

D. Chu, X. Liu, and V. Mehrmann (2007). “A Numerical Method for Computing the Hamiltonian Schur Form,” Numer. Math. 105, 375–412.

Generalized eigenvalue problems that involve Hamiltonian matrices also arise:

P. Benner, V. Mehrmann, and H. Xu (1998). “A Numerically Stable, Structure Preserving Method for Computing the Eigenvalues of Real Hamiltonian or Symplectic Pencils,” Numer. Math. 78, 329–358.


---

<!-- golub_450_499 -->

C. Mehl (2000). “Condensed Forms for Skew-Hamiltonian/Hamiltonian Pencils,” SIAM J. Matrix Anal. Applic. 21, 454–476.   
V. Mehrmann and D.S. Watkins (2001). “Structure–Preserving Methods for Computing Eigenpairs of Large Sparse Skew–Hamiltonian/Hamiltonian Pencils,” SIAM J. Sci. Comput. 22, 1905–1925.   
P. Benner and R. Byers, V. Mehrmann, and H. Xu (2002). “Numerical Computation of Deflating Subspaces of Skew-Hamiltonian/Hamiltonian Pencils,” SIAM J. Matrix Anal. Applic. 24, 165– 190.

Methods for symplectic eigenvalue problems are discussed in:

P. Benner, H. Fassbender and D.S. Watkins (1999). “SR and SZ Algorithms for the Symplectic (Butterfly) Eigenproblem,” Lin. Alg. Applic. 287, 41–76.

The Golub-Kahan SVD algorithm that we discuss in the next chapter does not form $A ^ { T } A$ or $A A ^ { T }$ despite the rich connection to the Schur decompositions of those matrices. From that point on there has been an appreciation for the numerical dangers associated with explicit products. Here is a sampling of the literature:

C. Van Loan (1975). “A General Matrix Eigenvalue Algorithm,” SIAM J. Numer. Anal. 12, 819–834.

M.T. Heath, A.J. Laub, C.C. Paige, and R.C. Ward (1986). “Computing the SVD of a Product of Two Matrices,” SIAM J. Sci. Stat. Comput. 7, 1147–1159.

R. Mathias (1998). “Analysis of Algorithms for Orthogonalizing Products of Unitary Matrices,” Num. Lin. Alg. 3, 125–145.

G. Golub, K. Solna, and P. Van Dooren (2000). “Computing the SVD of a General Matrix Product/Quotient,” SIAM J. Matrix Anal. Applic. 22, 1–19.

D.S. Watkins (2005). “Product Eigenvalue Problems,” SIAM Review 47, 3–40.

R. Granat and B. Kgstrom (2006). “Direct Eigenvalue Reordering in a Product of Matrices in Periodic Schur Form,” SIAM J. Matrix Anal. Applic. 28, 285–300.

Finally we mention that there is a substantial body of work concerned with structured error analysis and structured perturbation theory for structured matrix problems, see:

F. Tisseur (2003). “A Chart of Backward Errors for Singly and Doubly Structured Eigenvalue Problems,” SIAM J. Matrix Anal. Applic. 24, 877–897.

R. Byers and D. Kressner (2006). “Structured Condition Numbers for Invariant Subspaces,” SIAM J. Matrix Anal. Applic. 28, 326–347.

M. Karow, D. Kressner, and F. Tisseur (2006). “Structured Eigenvalue Condition Numbers,” SIAM J. Matrix Anal. Applic. 28, 1052–1068.

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
