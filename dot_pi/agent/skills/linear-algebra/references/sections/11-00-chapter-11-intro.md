# Chapter 11

# Large Sparse Linear System Problems

11.1 Direct Methods   
11.2 The Classical Iterations   
11.3 The Conjugate Gradient Method   
11.4 Other Krylov Methods   
11.5 Preconditioning   
11.6 The Multigrid Framework

This chapter is about solving linear systems and least squares problems when the matrix in question is so large and sparse that we have to rethink our powerful dense factorization strategies. The basic challenge is to live without the standard 2- dimensional array representation where there is a 1:1 correspondence between matrix entries and storage cells.

There is sometimes sufficient structure to actually compute an LU, Cholesky, or QR factorization by using a sparse matrix data structure and by carefully reordering equations and unknowns to control the fill-in of nonzero entries during the factorization process. Methods of this variety are called direct methods and they are the subject of §11.1. Our treatment is brief, touching only some of the high points of this well-developed area. A deeper presentation requires much more graph theory and implementation-based insight than we can provide in these few pages.

The rest of the chapter is concerned with the iterative method framework. These methods produce a sequence of vectors that typically converge to the solution at a reasonable rate. The matrix A “shows up” only in the context of matrix/vector multiplication. We introduce the strategy in §11.2 through discussion of the “classical” methods of Jacobi, Gauss-Seidel, successive over-relaxation, and Chebyshev. The discrete Poisson problem from §4.8.3 is used to reinforce the major ideas.

Krylov subspace methods are treated in the next two sections. In §11.3 we derive the method of conjugate gradients that is suitable for symmetric positive definite linear systems. The derivation involves the Lanczos process, the method of steepest descent, and the idea of optimizing over a nested sequence of subspaces. Related methods for symmetric indefinite systems, general systems, and least squares problems are covered in §11.4.

It is generally the case that Krylov subspace methods are successful only if there is an effective preconditioner. For a given Ax = b problem this essentially requires the design of a matrix M that has two properties. It must capture key features of A and it must be relatively easy to solve systems of the form $M z = r$ . There are several major families of preconditioners and these are surveyed in §11.5 and §11.6, the latter being dedicated to the mesh-coarsening/multigrid framework.

# Reading Path

An understanding of the basics about LU, Cholesky, and QR factorizations is essential. Eigenvalue theory and functions of matrices have a prominent role to play in the analysis of iterative Ax = b solvers. The Krylov methods make use of the Lanczos and Arnoldi iterations that we developed in Chapter 10.

Within this chapter, there are the following dependencies:

$$
\begin{array}{c c c c c c c c} \S 1 1. 2 & \to & \S 1 1. 3 & \to & \S 1 1. 4 & \to & \S 1 1. 5 \\ \downarrow & & & & & & \\ \S 1 1. 6 & & & & & & \end{array}
$$

§11.1 is independent of the others. The books by Axelsson (ISM), Greenbaum (IMSL), Saad (ISPLA), and van der Vorst (IMK) provide excellent background. The software “templates” volume LIN TEMPLATES (1993) is very useful for its concise presentation of all the major iterative strategies and for the guidance it provides in choosing a suitable method.
