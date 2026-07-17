# Chapter 3

# General Linear Systems

3.1 Triangular Systems   
3.2 The LU Factorization   
3.3 Roundoff Error in Gaussian Elimination   
3.4 Pivoting   
3.5 Improving and Estimating Accuracy   
3.6 Parallel LU

The problem of solving a linear system Ax = b is central to scientific computation. In this chapter we focus on the method of Gaussian elimination, the algorithm of choice if A is square, dense, and unstructured. Other methods are applicable if A does not fall into this category, see Chapter 4, Chapter 11, §12.1, and §12.2. Solution procedures for triangular systems are discussed first. These are followed by a derivation of Gaussian elimination that makes use of Gauss transformations. The process of eliminating unknowns from equations is described in terms of the factorization $A = L U$ where L is lower triangular and U is upper triangular. Unfortunately, the derived method behaves poorly on a nontrivial class of problems. An error analysis pinpoints the difficulty and sets the stage for a discussion of pivoting, a permutation strategy that keeps the numbers “nice” during the elimination. Practical issues associated with scaling, iterative improvement, and condition estimation are covered. A framework for computing the LU factorization in parallel is developed in the final section.

# Reading Notes

Familiarity with Chapter 1, §§2.1–2.5, and §2.7 is assumed. The sections within this chapter depend upon each other as follows:

$$
\begin{array}{c c c c c c c c} & & & & & \S 3. 5 \\ & & & & & \uparrow \\ \S 3. 1 & \to & \S 3. 2 & \to & \S 3. 3 & \to & \S 3. 4 \\ & & & & & \downarrow \\ & & & & & \S 3. 6 \end{array}
$$

Useful global references include Forsythe and Moler (SLAS), Stewart( MABD), Higham (ASNA), Watkins (FMC), Trefethen and Bau (NLA), Demmel (ANLA), and Ipsen (NMA).
