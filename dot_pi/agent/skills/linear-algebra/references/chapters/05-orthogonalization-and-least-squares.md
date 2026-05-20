# Chapter 5

# Orthogonalization and Least Squares

5.1 Householder and Givens Transformations   
5.2 The QR Factorization   
5.3 The Full-Rank Least Squares Problem   
5.4 Other Orthogonal Factorizations   
5.5 The Rank-Deficient Least Squares Problem   
5.6 Square and Underdetermined Systems

This chapter is primarily concerned with the least squares solution of overdetermined systems of equations, i.e., the minimization of $\parallel A x - b \parallel _ { 2 }$ where $A \in \mathbb { R } ^ { m \times n }$ , $b \in \mathbb { R } ^ { m }$ , and $m \geq n$ . The most reliable solution procedures for this problem involve the reduction of A to various canonical forms via orthogonal transformations. Householder reflections and Givens rotations are central to this process and we begin the chapter with a discussion of these important transformations. In §5.2 we show how to compute the factorization $A = Q R$ where Q is orthogonal and R is upper triangular. This amounts to finding an orthonormal basis for the range of A. The QR factorization can be used to solve the full-rank least squares problem as we show in §5.3. The technique is compared with the method of normal equations after a perturbation theory is developed. In §5.4 and §5.5 we consider methods for handling the difficult situation when A is (nearly) rank deficient. QR with column pivoting and other rank-revealing procedures including the SVD are featured. Some remarks about underdetermined systems are offered in §5.6.

# Reading Notes

Knowledge of chapters 1, 2, and 3 and §§4.1–§4.3 is assumed. Within this chapter there are the following dependencies:

$$
\begin{array}c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
$$

For more comprehensive treatments of the least squares problem, see Bj¨orck (NMLS) and Lawson and Hansen (SLS). Other useful global references include Stewart ( MABD), Higham (ASNA), Watkins (FMC), Trefethen and Bau (NLA), Demmel (ANLA), and Ipsen (NMA).

# 5.1 Householder and Givens Transformations

Recall that $Q \in \mathbb { R } ^ { m \times m }$ is orthogonal if

$$
Q ^ {T} Q = Q Q ^ {T} = I _ {m}.
$$

Orthogonal matrices have an important role to play in least squares and eigenvalue computations. In this section we introduce Householder reflections and Givens rotations, the key players in this game.

# 5.1.1 A 2-by-2 Preview

It is instructive to examine the geometry associated with rotations and reflections at the m = 2 level. A 2-by-2 orthogonal matrix Q is a rotation if it has the form

$$
Q = \left[ \begin{array}{c c} \cos (\theta) & \sin (\theta) \\ - \sin (\theta) & \cos (\theta) \end{array} \right].
$$

If $y = Q ^ { T } x ,$ , then y is obtained by rotating x counterclockwise through an angle θ. A 2-by-2 orthogonal matrix Q is a reflection if it has the form

$$
Q = \left[ \begin{array}{c c} \cos (\theta) & \sin (\theta) \\ \sin (\theta) & - \cos (\theta) \end{array} \right].
$$

If $y = Q ^ { T } x = Q x$ , then y is obtained by reflecting the vector x across the line defined by

$$
S = \text { span } \left\{\left[ \begin{array}{c} \cos (\theta / 2) \\ \sin (\theta / 2) \end{array} \right] \right\}.
$$

Reflections and rotations are computationally attractive because they are easily constructed and because they can be used to introduce zeros in a vector by properly choosing the rotation angle or the reflection plane.

# 5.1.2 Householder Reflections

Let $v \in \mathbb { R } ^ { m }$ be nonzero. An m-by-m matrix P of the form

$$
P = I - \beta v v ^ {T}, \quad \beta = \frac {2}{v ^ {T} v} \tag {5.1.1}
$$

is a Householder reflection. (Synonyms are Householder matrix and Householder transformation.) The vector v is the Householder vector. If a vector x is multiplied by P , then it is reflected in the hyperplane span $\{ v \} ^ { \perp }$ . It is easy to verify that Householder matrices are symmetric and orthogonal.

Householder reflections are similar to Gauss transformations introduced in §3.2.1 in that they are rank-1 modifications of the identity and can be used to zero selected components of a vector. In particular, suppose we are given $0 \neq x \in \mathbb { R } ^ { m }$ and want

$$
P x = \left(I - \frac {2 v v ^ {T}}{v ^ {T} v}\right) x = x - \frac {2 v ^ {T} x}{v ^ {T} v} v
$$

to be a multiple of $e _ { 1 } = I _ { m } ( : , 1 )$ . From this we conclude that $v \in { \mathrm { s p a n } } \{ x , e _ { 1 } \}$ . Setting

$$
v = x + \alpha e _ {1}
$$

gives

$$
v ^ {T} x = x ^ {T} x + \alpha x _ {1}
$$

and

$$
v ^ {T} v = x ^ {T} x + 2 \alpha x _ {1} + \alpha^ {2}.
$$

Thus,

$$
\begin{array}{l} P x = \left(1 - 2 \frac {x ^ {T} x + \alpha x _ {1}}{x ^ {T} x + 2 \alpha x _ {1} + \alpha^ {2}}\right) x - 2 \alpha \frac {v ^ {T} x}{v ^ {T} v} e _ {1} \\ = \left(\frac {\alpha^ {2} - \| x \| _ {2} ^ {2}}{x ^ {T} x + 2 \alpha x _ {1} + \alpha^ {2}}\right) x - 2 \alpha \frac {v ^ {T} x}{v ^ {T} v} e _ {1}. \\ \end{array}
$$

In order for the coefficient of x to be zero, we set $\alpha = \pm \| \boldsymbol { x } \| _ { 2 }$ for then

$$
v = x \pm \| x \| _ {2} e _ {1} \Rightarrow P x = \left(I - 2 \frac {v v ^ {T}}{v ^ {T} v}\right) x = \mp \| x \| _ {2} e _ {1}. \tag {5.1.2}
$$

It is this simple determination of v that makes the Householder reflections so useful.

# 5.1.3 Computing the Householder Vector

There are a number of important practical details associated with the determination of a Householder matrix, i.e., the determination of a Householder vector. One concerns the choice of sign in the definition of v in (5.1.2). Setting

$$
v _ {1} = x _ {1} - \parallel x \parallel_ {2}
$$

leads to the nice property that $P x$ is a positive multiple of $e _ { 1 }$ . But this recipe is dangerous if x is close to a positive multiple of $e _ { 1 }$ because severe cancellation would occur. However, the formula

$$
v _ {1} = x _ {1} - \parallel x \parallel_ {2} = \frac {x _ {1} ^ {2} - \parallel x \parallel_ {2} ^ {2}}{x _ {1} + \parallel x \parallel_ {2}} = \frac {- (x _ {2} ^ {2} + \cdots + x _ {n} ^ {2})}{x _ {1} + \parallel x \parallel_ {2}}
$$

suggested by Parlett (1971) does not suffer from this defect in the $x _ { 1 } > 0$ case.

In practice, it is handy to normalize the Householder vector so that $v ( 1 ) = 1$ . This permits the storage of v(2:m) where the zeros have been introduced in x, i.e., x(2:m). We refer to v(2:m) as the essential part of the Householder vector. Recalling that $\beta = 2 / v ^ { T } v$ and letting length(x) specify vector dimension, we may encapsulate the overall process as follows:

Algorithm 5.1.1 (Householder Vector) $\mathrm { G i v e n } x \in \mathbb { R } ^ { m }$ , this function computes $v \in \mathbb { R } ^ { m }$ with $v ( 1 ) = 1$ and $\beta \in \mathbb { R }$ such that $\begin{array} { r } { \dot { P } = I _ { m } - \beta v v ^ { T } } \end{array}$ is orthogonal and $\ P x = \| \boldsymbol { x } \| _ { 2 } e _ { 1 }$ .

function $\left\lceil v , \beta \right\rceil = \mathsf { h o u s e } ( x )$

$$
m = \mathsf {l e n g t h} (x), \sigma = x (2: m) ^ {T} x (2: m), v = \left[ \begin{array}{c} 1 \\ x (2: m) \end{array} \right]
$$

$\mathbf { i f } \ \sigma = 0 \ \mathrm { a n d } \ x ( 1 ) > = 0$

$$
\beta = 0
$$

elseif $\sigma = 0 \mathrm { ~ } \& \mathrm { ~ } x ( 1 ) < 0$

$$
\beta = - 2
$$

else

$$
\mu = \sqrt {x (1) ^ {2} + \sigma}
$$

$\mathbf { i f } \ x ( 1 ) < = 0$

$$
v (1) = x (1) - \mu
$$

$$
v (1) = - \sigma / (x (1) + \mu)
$$

end

$$
\beta = 2 v (1) ^ {2} / (\sigma + v (1) ^ {2})
$$

$$
v = v / v (1)
$$

end

Here, length(·) returns the dimension of a vector. This algorithm involves about 3m flops. The computed Householder matrix that is orthogonal to machine precision, a concept discussed below.

# 5.1.4 Applying Householder Matrices

It is critical to exploit structure when applying $P = I - \beta v v ^ { T }$ to a matrix A. Premultiplication involves a matrix-vector product and a rank-1 update:

$$
P A = (I - \beta v v ^ {T}) A = A - (\beta v) (v ^ {T} A).
$$

The same is true for post-multiplication,

$$
A P = A (I - \beta v v ^ {T}) = A - (A v) (\beta v) ^ {T}.
$$

In either case, the update requires 4mn flops if $A \in \mathbb { R } ^ { m \times n }$ . Failure to recognize this and to treat P as a general matrix increases work by an order of magnitude. Householder updates never entail the explicit formation of the Householder matrix.

In a typical situation, house is applied to a subcolumn or subrow of a matrix and $( I - \beta v v ^ { T } )$ is applied to a submatrix. For example, if $A \in \mathbb { R } ^ { m \times n } , \ 1 \leq j < n$ , and $A ( j { : } m , 1 { : } j - 1 )$ is zero, then the sequence

$$
[ v, \beta ] = \text { house } (A (j: m, j))
$$

$$
A (j: m, j: n) = A (j: m, j: n) - (\beta v) \left(v ^ {T} A (j: m, j: n)\right)
$$

$$
A (j + 1: m, j) = v (2: m - j + 1)
$$

applies $( I _ { m - j + 1 } - \beta v v ^ { T } )$ to $A ( j { : } m , 1 { : } n )$ and stores the essential part of v where the “new” zeros are introduced.

# 5.1.5 Roundoff Properties

The roundoff properties associated with Householder matrices are very favorable. Wilkinson (AEP, pp. 152–162) shows that house produces a Householder vector $\hat { v }$ that is very close to the exact v. If $\hat { P } = I - 2 \hat { v } \hat { v } ^ { T } / \hat { v } ^ { \hat { T } } \hat { v }$ then

$$
\| \hat {P} - P \| _ {2} = O (\mathbf {u}).
$$

Moreover, the computed updates with $\hat { P }$ are close to the exact updates with $P$ :

$$
\mathfrak {f l} (\hat {P} A) = P (A + E), \quad \| E \| _ {2} = O (\mathbf {u} \| A \| _ {2}),
$$

$$
\mathsf {f l} (A \hat {P}) = (A + E) P, \quad \| E \| _ {2} = O (\mathbf {u} \| A \| _ {2}).
$$

For a more detailed analysis, see Higham(ASNA, pp. 357–361).

# 5.1.6 The Factored-Form Representation

Many Householder-based factorization algorithms that are presented in the following sections compute products of Householder matrices

$$
Q = Q _ {1} Q _ {2} \dots Q _ {n} \quad Q _ {j} = I _ {m} - \beta_ {j} v ^ {(j)} [ v ^ {(j)} ] ^ {T} \tag {5.1.3}
$$

where $n \leq m$ and each $v ^ { ( j ) }$ has the form

$$
v ^ {(j)} = [ \underbrace {0 , 0 , \ldots 0} _ {j - 1}, 1 v _ {j + 1} ^ {(j)}, \ldots , v _ {m} ^ {(j)} ] ^ {T}.
$$

It is usually not necessary to compute Q explicitly even if it is involved in subsequent calculations. For example, if $\boldsymbol { C } \in \mathbb { R } ^ { m \times p }$ and we wish to compute $Q ^ { T } C$ , then we merely execute the loop

for $j = 1 { : } n$

$$
C = Q _ {j} C
$$

end

The storage of the Householder vectors $\boldsymbol { v } ^ { ( 1 ) } \cdots \boldsymbol { v } ^ { ( n ) }$ and the corresponding $\beta _ { j }$ amounts to a factored-form representation of $Q$ .

To illustrate the economies of the factored-form representation, suppose we have an array A and that for $j = 1 { : } n , A ( j + 1 { : } m , j )$ houses $\boldsymbol { v } ^ { ( j ) } ( j + 1 { : } m )$ , the essential part of the jth Householder vector. The overwriting of $\boldsymbol { C } \in \mathbb { R } ^ { m \times p }$ with $Q ^ { T } C$ can then be implemented as follows:

for $j = 1 { : } n$

$$
v (j: m) = \left[ \begin{array}{c} 1 \\ A (j + 1: m, j) \end{array} \right]
$$

$$
\beta_ {j} = 2 / (1 + \| A (j + 1: m, j) \| _ {2} ^ {2} \tag {5.1.4}
$$

$$
C (j: m,:) = C (j: m,:) - \left(\beta_ {j} \cdot v (j: m)\right) \cdot \left(v (j: m) ^ {T} C (j: m,:)\right)
$$

end

This involves about $p n ( 2 m - n )$ flops. If Q is explicitly represented as an m-by-m matrix, then $Q ^ { T } C$ would involve $2 m ^ { 2 } p$ flops. The advantage of the factored form representation is apparant if $n < < m$ .

Of course, in some applications, it is necessary to explicitly form Q (or parts of it). There are two possible algorithms for computing the matrix $Q$ in (5.1.3):

<table><tr><td>Forward accumulation</td><td>Backward accumulation</td></tr><tr><td> $Q = I_{m}$ </td><td> $Q = I_{m}$ </td></tr><tr><td>for  $j = 1:n$ </td><td>for  $j = n: -1:1$ </td></tr><tr><td> $Q = QQ_{j}$ </td><td> $Q = Q_{j}Q$ </td></tr><tr><td>end</td><td>end</td></tr></table>

Recall that the leading $( j \mathrm { ~ - ~ } 1 ) – \mathrm { b y } – ( j \mathrm { ~ - ~ } 1 )$ portion of $Q _ { j }$ is the identity. Thus, at the beginning of backward accumulation, Q is “mostly the identity” and it gradually becomes full as the iteration progresses. This pattern can be exploited to reduce the number of required flops. In contrast, Q is full in forward accumulation after the first step. For this reason, backward accumulation is cheaper and the strategy of choice. Here are the details with the proviso that we only need Q(:, 1:k) where $1 \leq k \leq m \colon$

$$
Q = I _ {m} (:, 1: k)
$$

for $j = n \colon - 1 \colon 1$

$$
v (j: m) = \left[ \begin{array}{c} 1 \\ A (j + 1: m, j) \end{array} \right] \tag {5.1.5}
$$

$$
\beta_ {j} = 2 / (1 + \parallel A (j + 1: m, j) \parallel_ {2} ^ {2}
$$

$$
Q (j: m, j: k) = Q (j: m, j: k) - (\beta_ {j} v (j: m)) (v (j: m) ^ {T} Q (j: m, j: k))
$$

end

This involves about 4mnk $- 2 ( m + k ) n ^ { 2 } + ( 4 / 3 ) n ^ { 3 }$ flops.

# 5.1.7 The WY Representation

Suppose $Q = Q _ { 1 } \cdot \cdot \cdot Q _ { r }$ is a product of m-by-m Householder matrices. Since each $Q _ { j }$ is a rank-1 modification of the identity, it follows from the structure of the Householder vectors that $Q$ is a rank-r modification of the identity and can be written in the form

$$
Q = I _ {m} - W Y ^ {T} \tag {5.1.6}
$$

where W and Y are m-by-r matrices. The key to computing the WY representation (5.1.6) is the following lemma.

Lemma 5.1.1. Suppose $Q \ = \ I _ { m } - W Y ^ { T }$ is an m-by-m orthogonal matrix with W, $Y \in \mathbb { R } ^ { m \times j }$ . If $P = I _ { m } - \beta v v ^ { T }$ with $v \in \mathbb { R } ^ { m }$ and $z = \beta Q v$ , then

$$
Q _ {+} = Q P = I _ {m} - W _ {+} Y _ {+} ^ {T}
$$

where $W _ { + } = \left[ W \mid z \right]$ and $Y _ { + } = \left[ Y \mid v \right]$ are each $m - b y - ( j + 1 )$ .

Proof. Since

$$
Q P = \left(I _ {m} - W Y ^ {T}\right) \left(I _ {m} - \beta v v ^ {T}\right) = I _ {m} - W Y ^ {T} - \beta Q v v ^ {T}
$$

it follows from the definition of z that

$$
Q _ {+} = I _ {m} - W Y ^ {T} - z v ^ {T} = I _ {m} - \left[ W \mid z \right] \left[ Y \mid v \right] ^ {T} = I _ {m} - W _ {+} Y _ {+} ^ {T}. \quad \square
$$

By repeatedly applying the lemma, we can transition from a factored-form representation to a block representation.

Algorithm 5.1.2 Suppose Q = Q1 · · · Qr where the Qj = Im − βjv(j) v(j)T $Q = Q _ { 1 } \cdot \cdot \cdot Q _ { r }$ $Q _ { j } = I _ { m } - \beta _ { j } v ^ { ( j ) } v ^ { ( j ) ^ { T } }$ are stored in factored form. This algorithm computes matrices $W , Y \in \mathbb { R } ^ { \bar { m } \times r }$ such that $Q =$ $I _ { m } - W Y ^ { T }$ .

$$
Y = v ^ {(1)}; W = \beta_ {1} v ^ {(1)}
$$

for $j = 2 { : } r$

$$
z = \beta_ {j} (I _ {m} - W Y ^ {T}) v ^ {(j)}
$$

$$
W = [ W \mid z ]
$$

$$
Y = \left[ Y \mid v ^ {(j)} \right]
$$

end

This algorithm involves about $2 r ^ { 2 } m - 2 r ^ { 3 } / 3$ flops if the zeros in the $v ^ { ( j ) }$ are exploited. Note that Y is merely the matrix of Householder vectors and is therefore unit lower triangular. Clearly, the central task in the generation of the WY representation (5.1.6) is the computation of the matrix W .

The block representation for products of Householder matrices is attractive in situations where $Q$ must be applied to a matrix. Suppose $\boldsymbol { C } \in \mathbb { R } ^ { m \times p }$ . It follows that the operation

$$
C = Q ^ {T} C = (I _ {m} - W Y ^ {T}) ^ {T} C = C - Y (W ^ {T} C)
$$

is rich in level-3 operations. On the other hand, if Q is in factored form, then the formation of $Q ^ { T } C$ is just rich in the level-2 operations of matrix-vector multiplication and outer product updates. Of course, in this context, the distinction between level-2 and level-3 diminishes as C gets narrower.

We mention that the WY representation (5.1.6) is not a generalized Householder transformation from the geometric point of view. True block reflectors have the form

$$
Q = I - 2 V V ^ {T}
$$

where $V \in \mathbb { R } ^ { n \times r }$ satisfies $V ^ { T } V = I _ { r }$ . See Schreiber and Parlett (1987).

# 5.1.8 Givens Rotations

Householder reflections are exceedingly useful for introducing zeros on a grand scale, e.g., the annihilation of all but the first component of a vector. However, in calculations where it is necessary to zero elements more selectively, Givens rotations are the transformation of choice. These are rank-2 corrections to the identity of the form

$$
G (i, k, \theta) = \left[ \begin{array}{c c c c c c c} 1 & \dots & 0 & \dots & 0 & \dots & 0 \\ \vdots & \ddots & \vdots & & \vdots & & \vdots \\ 0 & \dots & c & \dots & s & \dots & 0 \\ \vdots & & \vdots & \ddots & \vdots & & \vdots \\ 0 & \dots & - s & \dots & c & \dots & 0 \\ \vdots & & \vdots & & \vdots & \ddots & \vdots \\ 0 & \dots & 0 & \dots & 0 & \dots & 1 \\ & & i & & k \end{array} \right] \begin{array}{l} i \\ k \end{array} \tag {5.1.7}
$$

where $c = \cos ( \theta )$ and $s = \sin ( \theta )$ for some θ. Givens rotations are clearly orthogonal.

Premultiplication by $G ( i , k , \theta ) ^ { T }$ amounts to a counterclockwise rotation of θ radians in the (i, k) coordinate plane. Indeed, if $\boldsymbol { x } \in \mathbb { R } ^ { m }$ and

$$
y = G (i, k, \theta) ^ {T} x,
$$

then

$$
y _ {j} = \left\{ \begin{array}{c c} c x _ {i} - s x _ {k}, & j = i, \\ s x _ {i} + c x _ {k}, & j = k, \\ x _ {j}, & j \neq i, k. \end{array} \right..
$$

From these formulae it is clear that we can force $y _ { k }$ to be zero by setting

$$
c = \frac {x _ {i}}{\sqrt {x _ {i} ^ {2} + x _ {k} ^ {2}}}, \quad s = \frac {- x _ {k}}{\sqrt {x _ {i} ^ {2} + x _ {k} ^ {2}}}. \tag {5.1.8}
$$

Thus, it is a simple matter to zero a specified entry in a vector by using a Givens rotation. In practice, there are better ways to compute c and s than (5.1.8), e.g.,

Algorithm 5.1.3 Given scalars a and b, this function computes $c = \cos ( \theta )$ and $s = \sin ( \theta )$ so

$$
\left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} \left[ \begin{array}{c} a \\ b \end{array} \right] = \left[ \begin{array}{c} r \\ 0 \end{array} \right].
$$

$$
\begin{array}{l} \text { function } [ c, s ] = \text { givens } (a, b) \\ \text { if } b = 0 \\ c = 1; s = 0 \\ \text { if } | b | > | a | \\ \tau = - a / b; s = 1 / \sqrt {1 + \tau^ {2}}; c = s \tau \\ \tau = - b / a; c = 1 / \sqrt {1 + \tau^ {2}}; s = c \tau \\ \end{array}
$$

This algorithm requires 5 flops and a single square root. Note that inverse trigonometric functions are not involved.

# 5.1.9 Applying Givens Rotations

It is critical that the simple structure of a Givens rotation matrix be exploited when it is involved in a matrix multiplication. Suppose $A \in \mathbb { R } ^ { m \times n } , c = \cos ( \theta )$ , and $s = \sin ( \theta )$ . If $G ( i , k , \theta ) \in \mathbb { R } ^ { m \times m }$ , then the update $A = G ( i , k , \theta ) ^ { T } A$ affects just two rows,

$$
A ([ i, k ],:) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A ([ i, k ],:),
$$

and involves 6n flops:

$$
\begin{array}{l} \text { for } j = 1: n \\ \tau_ {1} = A (i, j) \\ \tau_ {2} = A (k, j) \\ A (i, j) = c \tau_ {1} - s \tau_ {2} \\ A (k, j) = s \tau_ {1} + c \tau_ {2} \\ \end{array}
$$

end

Likewise, if $G ( i , k , \theta ) \in \mathbb { R } ^ { n \times n }$ , then the update $A = A G ( i , k , \theta )$ affects just two columns,

$$
A (:, [ i, k ]) = A (:, [ i, k ]) \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right],
$$

and involves 6m flops:

$$
\begin{array}{l} \text {   for   } j = 1: m \\ \tau_ {1} = A (j, i) \\ \tau_ {2} = A (j, k) \\ A (j, i) = c \tau_ {1} - s \tau_ {2} \\ A (j, k) = s \tau_ {1} + c \tau_ {2} \\ \end{array}
$$

end

# 5.1.10 Roundoff Properties

The numerical properties of Givens rotations are as favorable as those for Householder reflections. In particular, it can be shown that the computed ˆc and ˆs in givens satisfy

$$
\begin{array}{l} \hat {c} = c (1 + \epsilon_ {c}), \quad \epsilon_ {c} = O (\mathbf {u}), \\ \hat {s} = s (1 + \epsilon_ {s}), \quad \epsilon_ {s} = O (\mathbf {u}). \\ \end{array}
$$

If ˆc and ˆs are subsequently used in a Givens update, then the computed update is the exact update of a nearby matrix:

$$
\mathsf {f l} [ \hat {G} (i, k, \theta) ^ {T} A ] = G (i, k, \theta) ^ {T} (A + E), \quad \| E \| _ {2} \approx \mathbf {u} \| A \| _ {2},
$$

$$
\mathsf {f l} [ A \hat {G} (i, k, \theta) ] = (A + E) G (i, k, \theta), \quad \| E \| _ {2} \approx \mathbf {u} \| A \| _ {2}.
$$

Detailed error analysis of Givens rotations may be found in Wilkinson (AEP, pp. 131- 39), Higham(ASNA, pp. 366–368), and Bindel, Demmel, Kahan, and Marques (2002).

# 5.1.11 Representing Products of Givens Rotations

Suppose $Q = G _ { 1 } \cdot \cdot \cdot G _ { t }$ is a product of Givens rotations. As with Householder reflections, it is sometimes more economical to keep Q in factored form rather than to compute explicitly the product of the rotations. Stewart (1976) has shown how to do this in a very compact way. The idea is to associate a single floating point number $\rho$ with each rotation. Specifically, if

$$
Z = \left[ \begin{array}{c c} {c} & {s} \\ {- s} & {c} \end{array} \right], \qquad c ^ {2} + s ^ {2} = 1,
$$

then we define the scalar $\rho$ by

$\mathbf { i f } \ c = 0$

$$
\rho = 1
$$

elseif $| s | < | c |$

$$
\rho = \operatorname{sign} (c) \cdot s / 2 \tag {5.1.9}
$$

$$
\rho = 2 \cdot \mathrm{sign} (s) / c
$$

end

Essentially, this amounts to storing $s / 2$ if the sine is smaller and $2 / c$ if the cosine is smaller. With this encoding, it is possible to reconstruct Z (or −Z) as follows:

if $\rho = 1$

$$
c = 0; s = 1
$$

elseif $| \rho | < 1$

$$
s = 2 \rho ; c = \sqrt {1 - s ^ {2}} \tag {5.1.10}
$$

$$
c = 2 / \rho ; s = \sqrt {1 - c ^ {2}}
$$

end

Note that the reconstruction of −Z is not a problem, for if $Z$ introduces a strategic zero then so does $- Z$ . The reason for essentially storing the smaller of c and s is that the formula $\sqrt { 1 - x ^ { 2 } }$ renders poor results if x is near unity. More details may be found in Stewart (1976). Of course, to “reconstruct” $G ( i , k , \theta )$ we need i and k in addition to the associated $\rho .$ This poses no difficulty if we agree to store $\rho$ in the (i, k) entry of some array.

# 5.1.12 Error Propagation

An m-by-m floating point matrix $\hat { Q }$ is orthogonal to working precision if there exists an orthogonal $Q \in \mathbb { R } ^ { m \times m }$ such that

$$
\| \hat {Q} - Q \| = O (\mathbf {u}).
$$

A corollary of this is that

$$
\| \hat {Q} ^ {T} \hat {Q} - I _ {m} \| = O (\mathbf {u}).
$$

The matrices defined by the floating point output of house and givens are orthogonal to working precision.

In many applications, sequences of Householders and/or Given transformations are generated and applied. In these settings, the rounding errors are nicely bounded. To be precise, suppose $A = A _ { 0 } \in \mathbb { R } ^ { m \times n }$ is given and that matrices $A _ { 1 } , \dotsc , A _ { p } = B$ are generated via the formula

$$
A _ {k} = \mathsf {f l} (\hat {Q} _ {k} A _ {k - 1} \hat {Z} _ {k}), \qquad k = 1: p  .
$$

Assume that the above Householder and Givens algorithms are used for both the generation and application of the $\hat { Q } _ { k }$ and $\hat { Z } _ { k }$ . Let $Q _ { k }$ and $Z _ { k }$ be the orthogonal matrices that would be produced in the absence of roundoff. It can be shown that

$$
B = (Q _ {p} \dots Q _ {1}) (A + E) (Z _ {1} \dots Z _ {p}), \tag {5.1.11}
$$

where $\| E \| _ { 2 } \ \leq \ c \cdot \mathbf { u } \| \ A \ \| _ { 2 }$ and c is a constant that depends mildly on n, m, and p. In other words, B is an exact orthogonal update of a matrix near to A. For a comprehensive error analysis of Householder and Givens computations, see Higham (ASNA, §19.3, §19.6).

# 5.1.13 The Complex Case

Most of the algorithms that we present in this book have complex versions that are fairly straightforward to derive from their real counterparts. (This is not to say that everything is easy and obvious at the implementation level.) As an illustration we briefly discuss complex Householder and complex Givens transformations.

Recall that if $A = ( a _ { i j } ) \in \mathbb { C } ^ { m \times n }$ , then $B = A ^ { H } \in \mathbb { C } ^ { n \times m }$ is its conjugate transpose. The 2-norm of a vector $x \in \mathbb { C } ^ { n }$ is defined by

$$
\parallel x \parallel_ {2} ^ {2} = x ^ {H} x = | x _ {1} | ^ {2} + \dots + | x _ {n} | ^ {2}
$$

and $Q \in \mathbb { C } ^ { n \times n }$ is unitary if $Q ^ { H } Q = I _ { n }$ . Unitary matrices preserve the 2-norm.

A complex Householder transformation is a unitary matrix of the form

$$
P = I _ {m} - \beta v v ^ {H}, \qquad 0 \neq v \in \mathbb {C} ^ {m},
$$

where $\beta = 2 / v ^ { H } v$ . Given a nonzero vector $\boldsymbol { x } \in \mathbb { C } ^ { m }$ , it is easy to determine v so that if $y = P x$ , then $y ( 2 { : } m ) = 0$ . Indeed, if

$$
x _ {1} = r e ^ {i \theta}
$$

where $r , \theta \in \mathbb { R }$ and

$$
v = x \pm e ^ {i \theta} \| x \| _ {2} e _ {1}, \qquad e _ {1} = I _ {m} (:, 1),
$$

then $P x = \mp e ^ { i \theta } \parallel x \parallel _ { 2 } e _ { 1 }$ . The sign can be determined to maximize $\parallel v \parallel _ { 2 }$ for the sake of stability.

Regarding complex Givens rotations, it is easy to verify that a 2-by-2 matrix of the form

$$
Q = \left[ \begin{array}{c c} \cos (\theta) & \sin (\theta) e ^ {i \phi} \\ - \sin (\theta) e ^ {- i \phi} & \cos (\theta) \end{array} \right]
$$

where $\theta , \phi \in \mathbb { R }$ is unitary. We show how to compute $c = \cos ( \theta )$ and $s = \sin ( \theta ) e ^ { i \phi }$ so that

$$
\left[ \begin{array}{c c} c & s \\ - \bar {s} & c \end{array} \right] ^ {H} \left[ \begin{array}{l} u \\ v \end{array} \right] = \left[ \begin{array}{l} r \\ 0 \end{array} \right] \tag {5.1.12}
$$

where $u = u _ { 1 } + i u _ { 2 }$ and $v = v _ { 1 } + i v _ { 2 }$ are given complex numbers. First, givens is applied to compute real cosine-sine pairs $\{ c _ { \alpha } , s _ { \alpha } \} , \{ c _ { \beta } , s _ { \beta } \}$ , and $\big \{ c _ { \theta } , s _ { \theta } \big \}$ so that

$$
\left[ \begin{array}{c c} c _ {\alpha} & s _ {\alpha} \\ - s _ {\alpha} & c _ {\alpha} \end{array} \right] ^ {T} \left[ \begin{array}{c} u _ {1} \\ u _ {2} \end{array} \right] = \left[ \begin{array}{c} r _ {u} \\ 0 \end{array} \right],
$$

$$
\left[ \begin{array}{c c} c _ {\beta} & s _ {\beta} \\ - s _ {\beta} & c _ {\beta} \end{array} \right] ^ {T} \left[ \begin{array}{c} v _ {1} \\ v _ {2} \end{array} \right] = \left[ \begin{array}{c} r _ {v} \\ 0 \end{array} \right],
$$

and

$$
\left[ \begin{array}{c c} c _ {\theta} & s _ {\theta} \\ - s _ {\theta} & c _ {\theta} \end{array} \right] ^ {T} \left[ \begin{array}{c} r _ {u} \\ r _ {v} \end{array} \right] = \left[ \begin{array}{c} r \\ 0 \end{array} \right].
$$

Note that $u = r _ { u } e ^ { - i \alpha }$ and $v = r _ { v } e ^ { - i \beta }$ . If we set

$$
e ^ {i \phi} = e ^ {i (\beta - \alpha)} = (c _ {\alpha} c _ {\beta} + s _ {\alpha} s _ {\beta}) + i (c _ {\alpha} s _ {\beta} - c _ {\beta} s _ {\alpha}),
$$

$c = c _ { \theta }$ , and $s = s _ { \theta } e ^ { i \phi }$ , then

$$
\bar {s} u + c v = s _ {\theta} e ^ {- i \phi} r _ {u} e ^ {- i \alpha} + c _ {\theta} r _ {v} e ^ {- i \beta} = e ^ {- i \beta} (s _ {\theta} r _ {u} + c _ {\theta} r _ {v}) = 0
$$

which confirms (5.1.12).

# Problems

P5.1.1 Let x and y be nonzero vectors in $\mathbb { R } ^ { m }$ . Give an algorithm for determining a Householder matrix P such that $P x$ is a multiple of y.

P5.1.2 Use Householder matrices to show that det $( I + x y ^ { T } ) = 1 + x ^ { T } y$ where x and y are given m-vectors.

P5.1.3 (a) Assume that $x , y \in \mathbb { R } ^ { 2 }$ have unit 2-norm. Give an algorithm that computes a Givens rotation Q so that $y = Q ^ { T } x$ . Make effective use of givens. (b) Suppose x and y are unit vectors in $\mathbb { R } ^ { m }$ . Give an algorithm using Givens transformations which computes an orthogonal Q such that $Q ^ { T } x = y$ .

P5.1.4 By generalizing the ideas in §5.1.11, develop a compact representation scheme for complex givens rotations.

P5.1.5 Suppose that $Q = I { - } Y T Y ^ { T }$ is orthogonal where $Y \in \mathbb { R } ^ { m \times j }$ and $T \in \mathbb { R } ^ { j \times j }$ is upper triangular. Show that if $Q _ { + } = Q P$ where $P = I - 2 v v ^ { T } \big / { v ^ { T } v }$ is a Householder matrix, then $Q _ { + }$ can be expressed in the form $Q _ { + } = I - Y _ { + } T _ { + } Y _ { + } ^ { T }$ where $Y _ { + } \stackrel { \cdot } { \in } \mathbb { R } ^ { m \times ( j + 1 ) }$ and $T _ { + } \in \mathbb { R } ^ { ( j + 1 ) \times ( j + 1 ) }$ is upper triangular. This is the main idea behind the compact WY representation. See Schreiber and Van Loan (1989).

P5.1.6 Suppose $Q _ { 1 } = I _ { m } - Y _ { 1 } T _ { 1 } Y _ { 1 }$ and $Q _ { 2 } \ = \ I _ { m } - Y _ { 2 } T _ { 2 } Y _ { 2 } ^ { T }$ are orthogonal where $Y _ { 1 } \in \mathbb { R } ^ { m \times r _ { 1 } }$ , $Y _ { 2 } \in \mathbb { R } ^ { m \times { r _ { 2 } } } , T _ { 1 } \in \mathbb { R } ^ { r _ { 1 } \times { r _ { 1 } } }$ , and $T _ { 2 } \in \mathbb { R } ^ { r _ { 2 } \times r _ { 2 } }$ . Assume that T1 and $T _ { 2 }$ are upper triangular. Show how to compute $Y \in \mathbb { R } ^ { m \times r }$ and upper triangular $T \in \mathbb { R } ^ { r \times r }$ with $r = r _ { 1 } + r _ { 2 }$ so that $Q _ { 2 } { \bar { Q _ { 1 } } } = I _ { m } - Y T Y ^ { T }$ .

P5.1.7 Give a detailed implementation of Algorithm 5.1.2 with the assumption that $\boldsymbol { v } ^ { ( j ) } ( j + 1 { : } m )$ , the essential part of the jth Householder vector, is stored in $A ( j + 1 { : } m , j )$ . Since Y is effectively represented in A, your procedure need only set up the W matrix.

P5.1.8 Show that if S is skew-symmetric $( S ^ { T } = - S )$ , then $Q = ( I + S ) ( I - S ) ^ { - 1 }$ is orthogonal. (The matrix Q is called the Cayley transform of S.) Construct a rank-2 S so that if x is a vector, then Qx is zero except in the first component.

P5.1.9 Suppose $P \in \mathbb { R } ^ { m \times m }$ satisfies $\parallel P ^ { T } P - I _ { m } \parallel _ { 2 } = \epsilon < 1$ . Show that all the singular values of P are in the interval $[ 1 - \epsilon , 1 + \epsilon ]$ and that $\| \ b { P } - \ b { U V } ^ { T } \| _ { 2 } \le \epsilon$ where $P = U \Sigma V ^ { T }$ is the SVD of $P ,$ .

P5.1.10 Suppose $A \in \mathbb { R } ^ { 2 \times 2 }$ . Under what conditions is the closest rotation to A closer than the closest reflection to A? Work with the Frobenius norm.

P5.1.11 How could Algorithm 5.1.3 be modified to ensure $r \geq 0 ?$

P5.1.12 (Fast Givens Transformations) Suppose

$$
x = \left[ \begin{array}{l} x _ {1} \\ x _ {2} \end{array} \right] \qquad \text {and} \qquad D = \left[ \begin{array}{l l} d _ {1} & 0 \\ 0 & d _ {2} \end{array} \right]
$$

with $d _ { 1 }$ and $d _ { 2 }$ positive. Show how to compute

$$
M _ {1} = \left[ \begin{array}{c c} \beta_ {1} & 1 \\ 1 & \alpha_ {1} \end{array} \right]
$$

so that if $y = M _ { 1 } x$ and $\tilde { D } = M _ { 1 } ^ { T } D M _ { 1 }$ , then $y _ { 2 } = 0$ and $\tilde { D }$ is diagonal. Repeat with $M _ { 1 }$ replaced by

$$
M _ {2} = \left[ \begin{array}{c c} 1 & \alpha_ {2} \\ \beta_ {2} & 1 \end{array} \right].
$$

(b) Show that either $\Vert \ M _ { 1 } ^ { T } D M _ { 1 } \ \Vert _ { 2 } \leq 2 \Vert \ D \ \Vert _ { 2 }$ or $| { \cal M } _ { 2 } ^ { T } { \cal D } { \cal M } _ { 2 } \| _ { 2 } \leq 2 \| { \cal D } \| _ { 2 }$ . (c) Suppose $\boldsymbol { x } \in \mathbb { R } ^ { m }$ and that $D \in \mathbb { R } ^ { n \times n }$ is diagonal with positive diagonal entries. Given indices i and j with $1 \leq i < j \leq m .$ , show how to compute $M \in \mathbb { R } ^ { n \times n }$ so that if $y = M x$ and $\tilde { D } = M ^ { T } D M$ , then $y _ { j } = 0$ and $\tilde { D }$ is diagonal with $\| \tilde { D } \| _ { 2 } \leq 2 \| D \| _ { 2 }$ . (d) From part (c) conclude that $Q = D ^ { 1 / 2 } M \tilde { D } ^ { - 1 / 2 }$ is orthogonal and that the update $y = M x$ can be diagonally transformed to $( D ^ { 1 / 2 } y ) = Q ( D ^ { 1 / 2 } x )$ .

# Notes and References for 5.1

Householder matrices are named after A.S. Householder, who popularized their use in numerical analysis. However, the properties of these matrices have been known for quite some time, see:

H.W. Turnbull and A.C. Aitken (1961). An Introduction to the Theory of Canonical Matrices, Dover Publications, New York, 102–105.

Other references concerned with Householder transformations include:

A.R. Gourlay (1970). “Generalization of Elementary Hermitian Matrices,” Comput. J. 13, 411–412.   
B.N. Parlett (1971). “Analysis of Algorithms for Reflections in Bisectors,” SIAM Review 13, 197–208.   
N.K. Tsao (1975). “A Note on Implementing the Householder Transformations.” SIAM J. Numer. Anal. 12, 53–58.   
B. Danloy (1976). “On the Choice of Signs for Householder Matrices,” J. Comput. Appl. Math. 2, 67–69.   
J.J.M. Cuppen (1984). “On Updating Triangular Products of Householder Matrices,” Numer. Math. 45, 403–410.   
A.A. Dubrulle (2000). “Householder Transformations Revisited,” SIAM J. Matrix Anal. Applic. 22, 33–40.   
J.W. Demmel, M. Hoemmen, Y. Hida, and E.J. Riedy (2009). “Nonnegative Diagonals and High Performance On Low-Profile Matrices from Householder QR,” SIAM J. Sci. Comput. 31, 2832– 2841.

A detailed error analysis of Householder transformations is given in Lawson and Hanson (SLE, pp. 83–89). The basic references for block Householder representations and the associated computations include:

C.H. Bischof and C. Van Loan (1987). “The WY Representation for Products of Householder Matrices,” SIAM J. Sci. Stat. Comput. 8, s2–s13.

B.N. Parlett and R. Schreiber (1988). “Block Reflectors: Theory and Computation,” SIAM J. Numer. Anal. 25, 189–205.   
R.S. Schreiber and C. Van Loan (1989). “A Storage-Efficient WY Representation for Products of Householder Transformations,” SIAM J. Sci. Stat. Comput. 10, 52–57.   
C. Puglisi (1992). “Modification of the Householder Method Based on the Compact WY Representation,” SIAM J. Sci. Stat. Comput. 13, 723–726.   
X. Sun and C.H. Bischof (1995). “A Basis-Kernel Representation of Orthogonal Matrices,” SIAM J. Matrix Anal. Applic. 16, 1184–1196.   
T. Joffrain, T.M. Low, E.S. Quintana-Orti, R. van de Geijn, and F.G. Van Zee (2006). “Accumulating Householder Transformations, Revisited,” ACM Trans. Math. Softw. 32, 169–179.   
M. Sadkane and A. Salam (2009). “A Note on Symplectic Block Reflectors,” ETNA 33, 45–52.   
Givens rotations are named after Wallace Givens. There are some subtleties associated with their computation and representation:   
G.W. Stewart (1976). “The Economical Storage of Plane Rotations,” Numer. Math. 25, 137–138.   
D. Bindel, J. Demmel, W. Kahan, and O. Marques (2002). “On computing givens rotations reliably and efficiently,” ACM Trans. Math. Softw. 28, 206–238.   
It is possible to aggregate rotation transformations to achieve high performance, see:   
B. Lang (1998). “Using Level 3 BLAS in Rotation–Based Algorithms,” SIAM J. Sci. Comput. 19, 626–634.   
Fast Givens transformations (see P5.1.11) are also referred to as square-root-free Givens transformations. (Recall that a square root must ordinarily be computed during the formation of Givens transformation.) There are several ways fast Givens calculations can be arranged, see:   
M. Gentleman (1973). “Least Squares Computations by Givens Transformations without Square Roots,” J. Inst. Math. Appl. 12, 329–336.   
C.F. Van Loan (1973). “Generalized Singular Values With Algorithms and Applications,” PhD Thesis, University of Michigan, Ann Arbor.   
S. Hammarling (1974). “A Note on Modifications to the Givens Plane Rotation,” J. Inst. Math. Applic. 13, 215–218.   
J.H. Wilkinson (1977). “Some Recent Advances in Numerical Linear Algebra,” in The State of the Art in Numerical Analysis, D.A.H. Jacobs (ed.), Academic Press, New York, 1–53.   
A.A. Anda and H. Park (1994). “Fast Plane Rotations with Dynamic Scaling,” SIAM J. Matrix Anal. Applic. 15, 162–174.   
R.J. Hanson and T. Hopkins (2004). “Algorithm 830: Another Visit with Standard and Modified Givens Transformations and a Remark on Algorithm 539,” ACM Trans. Math. Softw. 20, 86–94.

# 5.2 The QR Factorization

A rectangular matrix $A \in \mathbb { R } ^ { m \times n }$ can be factored into a product of an orthogonal matrix $Q \in \mathbb { R } ^ { m \times m }$ and an upper triangular matrix $R \in \mathbb { R } ^ { m \times n }$ :

$$
A = Q R.
$$

This factorization is referred to as the QR factorization and it has a central role to play in the linear least squares problem. In this section we give methods for computing QR based on Householder, block Householder, and Givens transformations. The QR factorization is related to the well-known Gram-Schmidt process.

# 5.2.1 Existence and Properties

We start with a constructive proof of the QR factorization.

Theorem 5.2.1 (QR Factorization). If $A \in \mathbb { R } ^ { m \times n }$ , then there exists an orthogonal $Q \in \mathbb { R } ^ { m \times m }$ and an upper triangular $R \in \mathbb { R } ^ { m \times n }$ so that $A = Q R$ .

Proof. We use induction. Suppose n = 1 and that Q is a Householder matrix so that if $R = Q ^ { T } A$ , then $R ( 2 { : } m ) = 0$ . It follows that $A = Q R$ is a QR factorization of A. For general n we partition A,

$$
A = \left[ \begin{array}{c c} A _ {1} & v \end{array} \right],
$$

where $v = A ( : , n )$ . By induction, there exists an orthogonal $Q _ { 1 } \in \mathbb { R } ^ { m \times m }$ so that $R _ { 1 } = Q _ { 1 } ^ { T } A _ { 1 }$ is upper triangular. Set $w = Q ^ { T } v$ and let $w ( n { : } m ) = Q _ { 2 } R _ { 2 }$ be the QR factorization of $w ( n { : } m )$ . If

$$
Q = Q _ {1} \left[ \begin{array}{c c} I _ {n - 1} & 0 \\ 0 & Q _ {2} \end{array} \right],
$$

then

$$
A = Q \left[ \begin{array}{c c c} R _ {1} & w (1: n - 1) \\ & R _ {2} \end{array} \right]
$$

is a QR factorization of A.

The columns of Q have an important connection to the range of A and its orthogonal complement.

Theorem 5.2.2. If A = QR is a QR factorization of a full column rank $A \in \mathbb { R } ^ { m \times n }$ and

$$
A = \left[ a _ {1} \mid \dots \mid a _ {n} \right],
$$

$$
Q = \left[ q _ {1} \mid \dots \mid q _ {m} \right]
$$

are column partitionings, then for k = 1:n

$$
\operatorname{span} \{a _ {1}, \dots , a _ {k} \} = \operatorname{span} \{q _ {1}, \dots , q _ {k} \} \tag {5.2.1}
$$

and $r _ { k k } \neq 0$ . Moreover, $i f Q _ { 1 } = Q ( 1 { : } m , 1 { : } n ) , Q _ { 2 } = Q ( 1 { : } m , n + 1 { : } m )$ , and $R _ { 1 } ~ =$ $R ( 1 { : } n , 1 { : } n )$ , then

$$
\operatorname{ran} (A) \quad = \operatorname{ran} (Q _ {1}),
$$

$$
\operatorname{ran} (A) ^ {\perp} = \operatorname{ran} (Q _ {2}),
$$

and

$$
A = Q _ {1} R _ {1}. \tag {5.2.2}
$$

Proof. Comparing the kth columns in A = QR we conclude that

$$
a _ {k} = \sum_ {i = 1} ^ {k} r _ {i k} q _ {i} \in \operatorname{span} \{q _ {1}, \dots , q _ {k} \}, \tag {5.2.3}
$$

and so

$$
\operatorname{span} \left\{a _ {1}, \dots , a _ {k} \right\} \subseteq \operatorname{span} \left\{q _ {1}, \dots , q _ {k} \right\}.
$$

If $r _ { k k } = 0$ , then $a _ { 1 } , \ldots , a _ { k }$ are dependent. Thus, R cannot have a zero on its diagonal and so span $\{ a _ { 1 } , \ldots , a _ { k } \}$ has dimension k. Coupled with (5.2.3) this establishes (5.2.1). To prove (5.2.2) we note that

$$
A = Q R = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \end{array} \right] \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right] = Q _ {1} R _ {1}. \quad \square
$$

The matrices $Q _ { 1 } = Q ( 1 { : } m , 1 { : } n )$ and $Q _ { 2 } = Q ( 1 { : } m , n + 1 { : } m )$ can be easily computed from a factored form representation of Q. We refer to (5.2.2) as the thin QR factorization. The next result addresses its uniqueness.

Theorem 5.2.3 (Thin QR Factorization). Suppose $A \in \mathbb { R } ^ { m \times n }$ has full column rank. The thin QR factorization

$$
A = Q _ {1} R _ {1}
$$

is unique where $Q _ { 1 } \in \mathbb { R } ^ { m \times n }$ has orthonormal columns and $R _ { 1 }$ is upper triangular with positive diagonal entries. Moreover, $R _ { 1 } = G ^ { T }$ where G is the lower triangular Cholesky factor of $A ^ { \bar { T } } A$ .

Proof. Since $A ^ { T } A = ( Q _ { 1 } R _ { 1 } ) ^ { T } ( Q _ { 1 } R _ { 1 } ) = R _ { 1 } ^ { T } R _ { 1 }$ we see that $G = R _ { 1 } ^ { T }$ is the Cholesky factor of $A ^ { T } A$ . This factor is unique by Theorem 4.2.7. Since $Q _ { 1 } \stackrel { - } { = } A R _ { 1 } ^ { - 1 }$ it follows that $Q _ { 1 }$ is also unique.

How are $Q _ { 1 }$ and $R _ { 1 }$ affected by perturbations in A? To answer this question we need to extend the notion of 2-norm condition to rectangular matrices. Recall from §2.6.2 that for square matrices, $\kappa _ { 2 } ( A )$ is the ratio of the largest to the smallest singular value. For rectangular matrices A with full column rank we continue with this definition:

$$
\kappa_ {2} (A) = \frac {\sigma_ {\max} (A)}{\sigma_ {\min} (A)}. \tag {5.2.4}
$$

If the columns of A are nearly dependent, then this quotient is large. Stewart (1993) has shown that $O ( \epsilon )$ relative error in A induces $O ( \epsilon { \cdot } \kappa _ { 2 } ( A ) )$ error in $Q _ { 1 }$ and $R _ { 1 }$ .

# 5.2.2 Householder QR

We begin with a QR factorization method that utilizes Householder transformations. The essence of the algorithm can be conveyed by a small example. Suppose $m = 6$ , $n = 5$ , and assume that Householder matrices $H _ { 1 }$ and $H _ { 2 }$ have been computed so that

$$
H _ {2} H _ {1} A = \left[ \begin{array}{c c c c c} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \mathbf {x} & \times & \times \\ 0 & 0 & \mathbf {x} & \times & \times \\ 0 & 0 & \mathbf {x} & \times & \times \\ 0 & 0 & \mathbf {x} & \times & \times \end{array} \right].
$$

Concentrating on the highlighted entries, we determine a Householder matrix $\tilde { H } _ { 3 } \in \mathbb { R } ^ { 4 \times 4 }$ such that

$$
\tilde {H} _ {3} \left[ \begin{array}{l} \mathbf {x} \\ \mathbf {x} \\ \mathbf {x} \\ \mathbf {x} \end{array} \right] = \left[ \begin{array}{l} \times \\ 0 \\ 0 \\ 0 \end{array} \right].
$$

If $H _ { 3 } = \mathrm { d i a g } ( I _ { 2 } , \tilde { H } _ { 3 } )$ , then

$$
H _ {3} H _ {2} H _ {1} A = \left[ \begin{array}{l l l l l} \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times & \times \end{array} \right].
$$

After n such steps we obtain an upper triangular $H _ { n } H _ { n - 1 } \cdot \cdot \cdot H _ { 1 } A = R$ and so by setting $Q = H _ { 1 } \cdot \cdot \cdot H _ { n }$ we obtain $A = Q R$ .

Algorithm 5.2.1 (Householder QR) Given $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n ,$ , the following algorithm finds Householder matrices $H _ { 1 } , \ldots , H _ { n }$ such that if $Q \ = \ H _ { 1 } \cdot \cdot \cdot H _ { n } $ , then $Q ^ { \bar { T } } A = R$ is upper triangular. The upper triangular part of A is overwritten by the upper triangular part of R and components $j + 1 { : } m$ of the jth Householder vector are stored in $A ( j + 1 ; m , j ) , j < m$ .

for $j = 1 { : } n$

$$
[ v, \beta ] = \operatorname{house} (A (j: m, j))
$$

$$
A (j: m, j: n) = (I - \beta v v ^ {T}) A (j: m, j: n)
$$

$$
\text { if } j <   m
$$

$$
A (j + 1: m, j) = v (2: m - j + 1)
$$

$$
\mathrm{end}
$$

end

This algorithm requires $2 n ^ { 2 } ( m - n / 3 )$ flops.

To clarify how A is overwritten, if

$$
v ^ {(j)} = [ \underbrace {0 , \ldots , 0} _ {j - 1}, 1, v _ {j + 1} ^ {(j)}, \ldots , v _ {m} ^ {(j)} ] ^ {T}
$$

is the jth Householder vector, then upon completion

$$
A = \left[ \begin{array}{l l l l l} r _ {1 1} & r _ {1 2} & r _ {1 3} & r _ {1 4} & r _ {1 5} \\ v _ {2} ^ {(1)} & r _ {2 2} & r _ {2 3} & r _ {2 4} & r _ {2 5} \\ v _ {3} ^ {(1)} & v _ {3} ^ {(2)} & r _ {3 3} & r _ {3 4} & r _ {3 5} \\ v _ {4} ^ {(1)} & v _ {4} ^ {(2)} & v _ {4} ^ {(3)} & r _ {4 4} & r _ {4 5} \\ v _ {5} ^ {(1)} & v _ {5} ^ {(2)} & v _ {5} ^ {(3)} & v _ {5} ^ {(4)} & r _ {5 5} \\ v _ {6} ^ {(1)} & v _ {6} ^ {(2)} & v _ {6} ^ {(3)} & v _ {6} ^ {(4)} & v _ {6} ^ {(5)} \end{array} \right].
$$

If the matrix $Q = H _ { 1 } \cdot \cdot \cdot H _ { n }$ is required, then it can be accumulated using (5.1.5). This accumulation requires $4 ( m ^ { 2 } n - m n ^ { 2 } + n ^ { 3 } / 3 )$ flops. Note that the β-values that arise in Algorithm 5.2.1 can be retrieved from the stored Householder vectors:

$$
\beta_ {j} = \frac {2}{1 + \parallel A (j + 1 : m , j) \parallel^ {2}}.
$$

We mention that the computed upper triangular matrix $\hat { R }$ is the exact $R$ for a nearby A in the sense that $Z ^ { T } ( \bar { A } + E ) = \hat { R }$ where $Z$ is some exact orthogonal matrix and $\parallel E \parallel _ { 2 } \approx \mathbf { u } \parallel A \parallel _ { 2 }$ .

# 5.2.3 Block Householder QR Factorization

Algorithm 5.2.1 is rich in the level-2 operations of matrix-vector multiplication and outer product updates. By reorganizing the computation and using the WY representation discussed in §5.1.7 we can obtain a level-3 procedure. The idea is to apply the underlying Householder transformations in clusters of size r. Suppose $n = 1 2$ and $r = 3$ . The first step is to generate Householders $H _ { 1 } , H _ { 2 }$ , and $H _ { 3 }$ as in Algorithm 5.2.1. However, unlike Algorithm 5.2.1 where each $H _ { i }$ is applied across the entire remaining submatrix, we apply only $H _ { 1 } , H _ { 2 }$ , and $H _ { 3 }$ to $A ( : , 1 { : } 3 )$ . After this is accomplished we generate the block representation $H _ { 1 } H _ { 2 } H _ { 3 } = I - W _ { 1 } Y _ { 1 } ^ { T }$ and then perform the level-3 update

$$
A (:, 4: 1 2) = (I - W Y ^ {T}) A (:, 4: 1 2).
$$

Next, we generate $H _ { 4 } , H _ { 5 }$ , and $H _ { 6 }$ as in Algorithm 5.2.1. However, these transformations are not applied to $A ( : , 7 : 1 2 )$ until their block representation $H _ { 4 } H _ { 5 } H _ { 6 } = I - W _ { 2 } Y _ { 2 } ^ { T }$ is found. This illustrates the general pattern.

Algorithm 5.2.2 (Block Householder QR) If $A \in \mathbb { R } ^ { m \times n }$ and r is a positive integer, then the following algorithm computes an orthogonal $Q \in \mathbb { R } ^ { m \times m }$ and an upper triangular $R \in \mathbb { R } ^ { m \times n }$ so that $A = Q R$ .

$$
Q = I _ {m}; \lambda = 1; k = 0
$$

while $\lambda \leq n$

$$
\tau \leftarrow \min (\lambda + r - 1, n); k = k + 1
$$

Use Algorithm 5.2.1, to upper triangularize $A ( \lambda { : } m , \lambda { : } \tau )$ , generating Householder matrices $H _ { \lambda } , \ldots , H _ { \tau }$ .

Use Algorithm 5.1.2 to get the block representation

$$
I - W _ {k} Y _ {k} = H _ {\lambda} \dots H _ {\tau}.
$$

$$
A (\lambda : m, \tau + 1: n) = (I - W _ {k} Y _ {k} ^ {T}) ^ {T} A (\lambda : m, \tau + 1: n)
$$

$$
Q (:, \lambda : m) = Q (:, \lambda : m) (I - W _ {k} Y _ {k} ^ {T})
$$

$$
\lambda = \tau + 1
$$

end

The zero-nonzero structure of the Householder vectors that define $H _ { \lambda } , \ldots , H _ { \tau }$ implies that the first $\lambda - 1$ rows of $W _ { k }$ and $Y _ { k }$ are zero. This fact would be exploited in a practical implementation.

The proper way to regard Algorithm 5.2.2 is through the partitioning

$$
A = \left[ A _ {1} \mid \dots \mid A _ {N} \right], \qquad N = \operatorname{ceil} (n / r)
$$

where block column $A _ { k }$ is processed during the kth step. In the kth step of the reduction, a block Householder is formed that zeros the subdiagonal portion of $A _ { k }$ . The remaining block columns are then updated.

The roundoff properties of Algorithm 5.2.2 are essentially the same as those for Algorithm 5.2.1. There is a slight increase in the number of flops required because of the W -matrix computations. However, as a result of the blocking, all but a small fraction of the flops occur in the context of matrix multiplication. In particular, the level-3 fraction of Algorithm 5.2.2 is approximately $1 - O ( 1 / N )$ . See Bischof and Van Loan (1987) for further details.

# 5.2.4 Block Recursive QR

A more flexible approach to blocking involves recursion. Suppose $A \in \mathbb { R } ^ { m \times n }$ and assume for clarity that A has full column rank. Partition the thin QR factorization of A as follows:

$$
\left[ \begin{array}{c c} A _ {1} & A _ {2} \end{array} \right] = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \end{array} \right] \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & R _ {2 2} \end{array} \right].
$$

where $n _ { 1 } = \mathrm { H o o r } ( n / 2 ) , n _ { 2 } = n - n _ { 1 } , A _ { 1 } , Q _ { 1 } \in \mathbb { R } ^ { m \times n _ { 1 } }$ and $A _ { 2 } , Q _ { 2 } \in \mathbb { R } ^ { m \times n _ { 2 } }$ . From the equations $Q _ { 1 } R _ { 1 1 } = A _ { 1 } , R _ { 1 2 } = Q _ { 1 } ^ { T } A _ { 2 }$ , and $Q _ { 2 } R _ { 2 2 } = A _ { 2 } - Q _ { 1 } R _ { 1 2 }$ we obtain the following recursive procedure:

Algorithm 5.2.3 (Recursive Block QR) Suppose $A \in \mathbb { R } ^ { m \times n }$ has full column rank and $n _ { b }$ is a positive blocking parameter. The following algorithm computes $Q \in \mathbb { R } ^ { m \times n }$ with orthonormal columns and upper triangular $R \in \mathbb { R } ^ { n \times n }$ such that $A = Q R$ .

function [Q, R] = BlockQR(A, n, nb)

if $n \leq n _ { b }$

Use Algorithm 5.2.1 to compute the thin QR factorization $A = Q R$

$$
n _ {1} = \operatorname{floor} (n / 2)
$$

$$
[ Q _ {1}, R _ {1 1} ] = \text { BlockQR } (A (:, 1: n _ {1}), n _ {1}, n _ {b})
$$

$$
R _ {1 2} = Q _ {1} ^ {T} A (:, n _ {1} + 1: n)
$$

$$
A (:, n _ {1} + 1: n) = A (:, n _ {1} + 1: n) - Q _ {1} R _ {1 2}
$$

$$
[ Q _ {2}, R _ {2 2} ] = \operatorname{BlockQR} (A (:, n _ {1} + 1: n), n - n _ {1}, n _ {b})
$$

$$
Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \end{array} \right], R = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & R _ {2 2} \end{array} \right]
$$

end end

This divide-and-conquer approach is rich in matrix-matrix multiplication and provides a framework for the effective parallel computation of the QR factorization. See Elmroth and Gustavson (2001). Key implementation ideas concern the representation of the $Q -$ matrices and the incorporation of the §5.2.3 blocking strategies.

# 5.2.5 Givens QR Methods

Givens rotations can also be used to compute the QR factorization and the 4-by-3 case illustrates the general idea:

$$
\begin{array}{l} \left[ \begin{array}{c c c} \times & \times & \times \\ \times & \times & \times \\ \mathbf {X} & \times & \times \\ \mathbf {X} & \times & \times \end{array} \right] \xrightarrow {(3 , 4)} \left[ \begin{array}{c c c} \times & \times & \times \\ \mathbf {X} & \times & \times \\ \mathbf {X} & \times & \times \\ 0 & \times & \times \end{array} \right] \xrightarrow {(2 , 3)} \left[ \begin{array}{c c c} \mathbf {X} & \times & \times \\ \mathbf {X} & \times & \times \\ 0 & \times & \times \\ 0 & \times & \times \end{array} \right] \xrightarrow {(1 , 2)} \\ \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & \mathbf {x} & \times \\ 0 & \mathbf {x} & \times \end{array} \right] \stackrel {{(3, 4)}} {{\longrightarrow}} \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \mathbf {x} & \times \\ 0 & \mathbf {x} & \times \\ 0 & 0 & \times \end{array} \right] \stackrel {{(2, 3)}} {{\longrightarrow}} \left[ \begin{array}{c c c} \times & \times & \times \\ 0 & \times & \times \\ 0 & 0 & \mathbf {x} \\ 0 & 0 & \mathbf {x} \end{array} \right] \stackrel {{(3, 4)}} {{\longrightarrow}} R. \\ \end{array}
$$

We highlighted the 2-vectors that define the underlying Givens rotations. If $G _ { j }$ denotes the jth Givens rotation in the reduction, then $Q ^ { \dot { T } } A \stackrel { - } { = } R$ is upper triangular, where $Q = G _ { 1 } \cdot \cdot \cdot G _ { t }$ and t is the total number of rotations. For general m and n we have:

Algorithm 5.2.4 (Givens QR) Given $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ , the following algorithm overwrites A with $Q ^ { T } A = R$ , where R is upper triangular and $Q$ is orthogonal.

for j = 1:n

for $i = m \colon - 1 { : } j + 1$

$$
[ c, s ] = \operatorname{givens} (A (i - 1, j), A (i, j))
$$

$$
A (i - 1: i, j: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A (i - 1: i, j: n)
$$

end end

This algorithm requires $3 n ^ { 2 } ( m - n / 3 )$ flops. Note that we could use the representation ideas from §5.1.11 to encode the Givens transformations that arise during the calculation. Entry $A ( i , j )$ can be overwritten with the associated representation.

With the Givens approach to the QR factorization, there is flexibility in terms of the rows that are involved in each update and also the order in which the zeros are introduced. For example, we can replace the inner loop body in Algorithm 5.2.4 with

$$
[ c, s ] = \operatorname{givens} (A (j, j), A (i, j))
$$

$$
A ([ j i ], j: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A ([ j i ], j: n)
$$

and still emerge with the QR factorization. It is also possible to introduce zeros by row. Whereas Algorithm 5.2.4 introduces zeros by column,

$$
\left[ \begin{array}{c c c} \times & \times & \times \\ 3 & \times & \times \\ 2 & 5 & \times \\ 1 & 4 & 6 \end{array} \right],
$$

the implementation

$$
\begin{array}{l} \text { for } i = 2: m \\ \text { for } j = 1: i - 1 \\ [ c, s ] = \operatorname{givens} (A (j, j), A (i, j)) \\ A ([ j i ], j: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A ([ j i ], j: n) \\ \end{array}
$$

introduces zeros by row, e.g.,

$$
\left[ \begin{array}{c c c} \times & \times & \times \\ 1 & \times & \times \\ 2 & 3 & \times \\ 4 & 5 & 6 \end{array} \right].
$$

# 5.2.6 Hessenberg QR via Givens

As an example of how Givens rotations can be used in a structured problem, we show how they can be employed to compute the QR factorization of an upper Hessenberg matrix. (Other structured QR factorizations are discussed in Chapter 6 and §11.1.8.) A small example illustrates the general idea. Suppose $n = 6$ and that after two steps we have computed

$$
G (2, 3, \theta_ {2}) ^ {T} G (1, 2, \theta_ {1}) ^ {T} A = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & \mathbf {x} & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right].
$$

Next, we compute $G ( 3 , 4 , \theta _ { 3 } )$ to zero the current (4,3) entry, thereby obtaining

$$
G (3, 4, \theta_ {3}) ^ {T} G (2, 3, \theta_ {2}) ^ {T} G (1, 2, \theta_ {1}) ^ {T} A = \left[ \begin{array}{l l l l l l} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right].
$$

Continuing in this way we obtain the following algorithm.

Algorithm 5.2.5 (Hessenberg QR) If $A \in \mathbb { R } ^ { n \times n }$ is upper Hessenberg, then the following algorithm overwrites A with $Q ^ { T } A = R$ where Q is orthogonal and R is upper triangular. $Q = G _ { 1 } \cdot \cdot \cdot G _ { n - 1 }$ is a product of Givens rotations where $G _ { j }$ has the form $G _ { j } = G ( j , j + 1 , \theta _ { j } )$ .

for $j = 1 \colon n - 1$

$$
[ c, s ] = \operatorname{givens} (A (j, j), A (j + 1, j))
$$

$$
A (j: j + 1, j: n) = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] ^ {T} A (j: j + 1, j: n)
$$

end

This algorithm requires about $3 n ^ { 2 }$ flops.

# 5.2.7 Classical Gram-Schmidt Algorithm

We now discuss two alternative methods that can be used to compute the thin QR factorization $A = Q _ { 1 } R _ { 1 }$ directly. If $\mathsf { r a n k } ( A ) = n$ , then equation (5.2.3) can be solved for $q _ { k }$ :

$$
q _ {k} = \left. \left(a _ {k} - \sum_ {i = 1} ^ {k - 1} r _ {i k} q _ {i}\right) \right/ r _ {k k}.
$$

Thus, we can think of $q _ { k }$ as a unit 2-norm vector in the direction of

$$
z _ {k} = a _ {k} - \sum_ {i = 1} ^ {k - 1} r _ {i k} q _ {i}
$$

where to ensure $z _ { k } \in \mathsf { s p a n } \{ q _ { 1 } , \dots , q _ { k - 1 } \} ^ { \perp }$ we choose

$$
r _ {i k} = q _ {i} ^ {T} a _ {k}, \qquad i = 1: k - 1.
$$

This leads to the classical Gram-Schmidt (CGS) algorithm for computing $A = Q _ { 1 } R _ { 1 }$ .

$$
R (1, 1) = \left\| A (:, 1) \right\| _ {2}
$$

$$
Q (:, 1) = A (:, 1) / R (1, 1)
$$

for $k = 2 { : } n$

$$
\begin{array}{l} R (1: k - 1, k) = Q (1: m, 1: k - 1) ^ {T} A (1: m, k) \\ z = A (1: m, k) - Q (1: m, 1: k - 1) \cdot R (1: k - 1, k) \\ R (k, k) = \parallel z \parallel_ {2} \\ Q (1: m, k) = z / R (k, k) \\ \end{array}
$$

end

In the kth step of CGS, the kth columns of both Q and R are generated.

# 5.2.8 Modified Gram-Schmidt Algorithm

Unfortunately, the CGS method has very poor numerical properties in that there is typically a severe loss of orthogonality among the computed $q _ { i }$ . Interestingly, a rearrangement of the calculation, known as modified Gram-Schmidt (MGS), leads to a more reliable procedure. In the kth step of MGS, the kth column of $Q$ (denoted by $q _ { k } )$ and the kth row of R (denoted by $r _ { k } ^ { T } )$ are determined. To derive the MGS method, define the matrix $A ^ { ( k ) } \in \mathbb { R } ^ { m \times ( n - k + 1 ) }$ by

$$
[ 0 \mid A ^ {(k)} ] = A - \sum_ {i = 1} ^ {k - 1} q _ {i} r _ {i} ^ {T} = \sum_ {i = k} ^ {n} q _ {i} r _ {i} ^ {T}.
$$

It follows that if

$$
A ^ {(k)} = \left[ \begin{array}{c c} z & B \\ 1 & n - k \end{array} \right]
$$

then $r _ { k k } ~ = ~ \parallel { z } \parallel _ { 2 } , q _ { k } ~ = ~ z / r _ { k k }$ , and $[ r _ { k , k + 1 } , . . . , r _ { k n } ] \ = \ q _ { k } ^ { T } B$ . We then compute the outer product $A ^ { ( k + 1 ) } = B ~ - ~ q _ { k } \left[ r _ { k , k + 1 } \cdot \cdot \cdot r _ { k n } \right]$ and proceed to the next step. This completely describes the kth step of MGS.

Algorithm 5.2.6 (Modified Gram-Schmidt) Given $A \in \mathbb { R } ^ { m \times n }$ with rank $( A ) = n$ , the following algorithm computes the thin QR factorization $A = Q _ { 1 } R _ { 1 }$ where $Q _ { 1 } \in \mathbb { R } ^ { m \times n }$ has orthonormal columns and $R _ { 1 } \in \mathbb { R } ^ { n \times n }$ is upper triangular.

$$
\begin{array}{l} R (k, k) = \left\| A (1: m, k) \right\| _ {2} \\ Q (1: m, k) = A (1: m, k) / R (k, k) \\ R (k, j) = Q (1: m, k) ^ {T} A (1: m, j) \\ A (1: m, j) = A (1: m, j) - Q (1: m, k) R (k, j) \\ \end{array}
$$

This algorithm requires $2 m n ^ { 2 }$ flops. It is not possible to overwrite A with both $Q _ { 1 }$ and $R _ { 1 }$ . Typically, the MGS computation is arranged so that A is overwritten by $Q _ { 1 }$ and the matrix $R _ { 1 }$ is stored in a separate array.

# 5.2.9 Work and Accuracy

If one is interested in computing an orthonormal basis for ran(A), then the Householder approach requires $2 m n ^ { 2 } - 2 n ^ { 3 } / 3$ flops to get $Q$ in factored form and another $2 m n ^ { 2 } -$ $2 n ^ { 3 } / 3$ flops to get the first n columns of Q. (This requires “paying attention” to just the first n columns of $Q$ in (5.1.5).) Therefore, for the problem of finding an orthonormal basis for ran(A), MGS is about twice as efficient as Householder orthogonalization. However, Bj¨orck (1967) has shown that MGS produces a computed ${ \hat { Q } } _ { 1 } = \left[ { \hat { q } } _ { 1 } | \cdots | { \hat { q } } _ { n } \right]$ that satisfies

$$
\hat {Q} _ {1} ^ {T} \hat {Q} _ {1} = I + E _ {M G S}, \quad \| E _ {M G S} \| _ {2} \approx \mathbf {u} \kappa_ {2} (A),
$$

whereas the corresponding result for the Householder approach is of the form

$$
\hat {Q} _ {1} ^ {T} \hat {Q} _ {1} = I + E _ {H}, \quad \| E _ {H} \| _ {2} \approx \mathbf {u}.
$$

Thus, if orthonormality is critical, then MGS should be used to compute orthonormal bases only when the vectors to be orthogonalized are fairly independent.

We also mention that the computed triangular factor $\hat { R }$ produced by MGS satisfies $\Vert A - \hat { Q } \hat { R } \Vert \approx \mathbf { u } \Vert A \Vert$ and that there exists a Q with perfectly orthonormal columns such that $\parallel A - Q { \hat { R } } \parallel \approx \mathbf { u } \parallel A \parallel$ . See Higham (ASNA, p. 379) and additional references given at the end of this section.

# 5.2.10 A Note on Complex Householder QR

Complex Householder transformations (§5.1.13) can be used to compute the QR factorization of a complex matrix $A \in \mathbb { C } ^ { m \times n }$ . Analogous to Algorithm 5.2.1 we have

for $j = 1 { : } n$

Compute a Householder matrix $Q _ { j }$ so that $Q _ { j } A$ is upper triangular through its first j columns.

$$
A = Q _ {j} A
$$

end

Upon termination, A has been reduced to an upper triangular matrix $R \in \mathbb { C } ^ { m \times n }$ and we have $A = Q R$ where $Q = Q _ { 1 } \cdot \cdot \cdot Q _ { n }$ is unitary. The reduction requires about four times the number of flops as the real case.

# Problems

P5.2.1 Adapt the Householder QR algorithm so that it can efficiently handle the case when $A \in \mathbb { R } ^ { m \times n }$ has lower bandwidth p and upper bandwidth q.

P5.2.2 Suppose $A \in \mathbb { R } ^ { n \times n }$ and let E be the exchange permutation $\mathcal { E } _ { n }$ obtained by reversing the order of the rows in ${ { I } _ { n } } . \mathrm { { \Gamma } ( a ) }$ Show that if $R \in \mathbb { R } ^ { n \times n }$ is upper triangular, then $L = \mathcal { E } R \mathcal { \bar { E } }$ is lower triangular. (b) Show how to compute an orthogonal $Q \in \mathbb { R } ^ { n \times n }$ and a lower triangular $\ b { L } \in \mathbb { R } ^ { n \times n }$ so that $A = Q L$ assuming the availability of a procedure for computing the QR factorization.

P5.2.3 Adapt the Givens QR factorization algorithm so that the zeros are introduced by diagonal. That is, the entries are zeroed in the order $( m , 1 ) , ( m - 1 , 1 ) , ( m , 2 ) , ( m - 2 , 1 ) , ( m - 1 , 2 ) , ( m , 3 )$ , etc.

P5.2.4 Adapt the Givens QR factorization algorithm so that it efficiently handles the case when A is n-by-n and tridiagonal. Assume that the subdiagonal, diagonal, and superdiagonal of A are stored in $e ( 1 { : } n - 1 ) , a ( 1 { : } n ) , f ( 1 { : } n - 1 )$ , respectively. Design your algorithm so that these vectors are overwritten by the nonzero portion of T .

P5.2.5 Suppose $\boldsymbol { L } \in \mathbb { R } ^ { m \times n }$ with $m \ \geq \ n$ is lower triangular. Show how Householder matrices $H _ { 1 } , \ldots , H _ { n }$ can be used to determine a lower triangular $L _ { 1 } \in \mathbb { R } ^ { n \times n }$ so that

$$
H _ {n} \dots H _ {1} L = \left[ \begin{array}{c} L _ {1} \\ 0 \end{array} \right].
$$

Hint: The second step in the 6-by-3 case involves finding $H _ { 2 }$ so that

$$
H _ {2} \left[ \begin{array}{c c c} \times & 0 & 0 \\ \times & \times & 0 \\ \times & \times & \times \\ \times & \times & 0 \\ \times & \times & 0 \\ \times & \times & 0 \end{array} \right] = \left[ \begin{array}{c c c} \times & 0 & 0 \\ \times & \times & 0 \\ \times & \times & \times \\ \times & 0 & 0 \\ \times & 0 & 0 \\ \times & 0 & 0 \end{array} \right]
$$

with the property that rows 1 and 3 are left alone.

P5.2.6 Suppose $A \in \mathbb { R } ^ { n \times n }$ and $D = \operatorname { d i a g } ( d _ { 1 } , \ldots , d _ { n } ) \in \mathbb { R } ^ { n \times n }$ . Show how to construct an orthogonal Q such that

$$
Q ^ {T} A - D Q ^ {T} = R
$$

is upper triangular. Do not worry about efficiency—this is just an exercise in QR manipulation.

P5.2.7 Show how to compute the QR factorization of the product

$$
A = A _ {p} \dots A _ {2} A _ {1}
$$

without explicitly multiplying the matrices $A _ { 1 } , \dotsc , A _ { p }$ together. Assume that each $A _ { i }$ is square. Hint: In the $p = 3 \ \mathrm { c a s e }$ , write

$$
Q _ {3} ^ {T} A = Q _ {3} ^ {T} A _ {3} Q _ {2} Q _ {2} ^ {T} A _ {2} Q _ {1} Q _ {1} ^ {T} A _ {1}
$$

and determine orthogonal $Q _ { i }$ so that $Q _ { i } ^ { T } ( A _ { i } Q _ { i - 1 } )$ is upper triangular. $( Q _ { 0 } = I . )$

P5.2.8 MGS applied to $A \in \mathbb { R } ^ { m \times n }$ is numerically equivalent to the first step in Householder QR applied to

$$
\tilde {A} = \left[ \begin{array}{c} O _ {n} \\ A \end{array} \right]
$$

where $O _ { n }$ is the $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ zero matrix. Verify that this statement is true after the first step of each method is completed.

P5.2.9 Reverse the loop orders in Algorithm 5.2.6 (MGS) so that R is computed column by column.

P5.2.10 How many flops are required by the complex QR factorization procedure outlined in §5.10?

P5.2.11 Develop a complex version of the Givens QR factorization in which the diagonal of R is nonnegative. See §5.1.13.

P5.2.12 Show that if $A \in \mathbb { R } ^ { n \times n }$ and $a _ { i } = A ( : , i )$ , then

$$
| \det (A) | \leq \| a _ {1} \| _ {2} \dots \| a _ {n} \| _ {2}.
$$

Hint: Use the QR factorization.

P5.2.13 Suppose $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ Construct an orthogonal $Q \in \mathbb { R } ^ { ( m + n ) \times ( m + n ) }$ with the property that Q(1:m, 1:n) is a scalar multiple of A. Hint. If $\alpha \in \mathbb { R }$ is chosen properly, then $I - \alpha ^ { 2 } A ^ { T } A$ has a Cholesky factorization.

P5.2.14 Suppose $A \in \mathbb { R } ^ { m \times n }$ . Analogous to Algorithm 5.2.4, show how fast Givens transformations (P5.1.12) can be used to compute $\bar { M } \in \mathbb { R } ^ { m \times m }$ and a diagonal $D \in \mathbb { R } ^ { m \times m }$ with positive diagonal entries so that $M ^ { T } A = S$ is upper triangular and $M M ^ { T } = \bar { D }$ . Relate M and S to A’s QR factors.

P5.2.15 (Parallel Givens QR) Suppose $A \in \mathbb { R } ^ { 9 \times 3 }$ and that we organize a Givens QR so that the subdiagonal entries are zeroed over the course of ten “time steps” as follows:

<table><tr><td>Step</td><td colspan="3">Entries Zeroed</td></tr><tr><td>T=1</td><td>(9,1)</td><td></td><td></td></tr><tr><td>T=2</td><td>(8,1)</td><td></td><td></td></tr><tr><td>T=3</td><td>(7,1)</td><td>(9,2)</td><td></td></tr><tr><td>T=4</td><td>(6,1)</td><td>(8,2)</td><td></td></tr><tr><td>T=5</td><td>(5,1)</td><td>(7,2)</td><td>(9,3)</td></tr><tr><td>T=6</td><td>(4,1)</td><td>(6,2)</td><td>(8,3)</td></tr><tr><td>T=7</td><td>(3,1)</td><td>(5,2)</td><td>(7,3)</td></tr><tr><td>T=8</td><td>(2,1)</td><td>(4,2)</td><td>(6,3)</td></tr><tr><td>T=9</td><td></td><td>(3,2)</td><td>(5,3)</td></tr><tr><td>T=10</td><td></td><td></td><td>(4,3)</td></tr></table>

Assume that a rotation in plane $( i - 1 , i )$ is used to zero a matrix entry $( i , j )$ . It follows that the rotations associated with any given time step involve disjoint pairs of rows and may therefore be computed in parallel. For example, during time step $T = 6 ,$ , there is a $( 3 , 4 ) , ( 5 , 6 )$ , and $^ { ( 7 , 8 ) }$ rotation. Three separate processors could oversee the three updates. Extrapolate from this example to the m-by-n case and show how the QR factorization could be computed in $O ( m + n )$ time steps. How many of those time steps would involve n “nonoverlapping” rotations?

# Notes and References for §5.2

The idea of using Householder transformations to solve the least squares problem was proposed in:

A.S. Householder (1958). “Unitary Triangularization of a Nonsymmetric Matrix,” J. ACM 5, 339–342.   
The practical details were worked out in:   
P. Businger and G.H. Golub (1965). “Linear Least Squares Solutions by Householder Transformations,” Numer. Math. 7, 269–276.   
G.H. Golub (1965). “Numerical Methods for Solving Linear Least Squares Problems,” Numer. Math. 7, 206–216.   
The basic references for Givens QR include:   
W. Givens (1958). “Computation of Plane Unitary Rotations Transforming a General Matrix to Triangular Form,” SIAM J. Appl. Math. 6, 26–50.   
M. Gentleman (1973). “Error Analysis of QR Decompositions by Givens Transformations,” Lin. Alg. Applic. 10, 189–197.   
There are modifications for the QR factorization that make it more attractive when dealing with rank deficiency. See §5.4. Nevertheless, when combined with the condition estimation ideas in §3.5.4, the traditional QR factorization can be used to address rank deficiency issues:   
L.V. Foster (1986). “Rank and Null Space Calculations Using Matrix Decomposition without Column Interchanges,” Lin. Alg. Applic. 74, 47–71.   
The behavior of the Q and R factors when A is perturbed is of interest. A main result is that the resulting changes in Q and R are bounded by the condition of A times the relative change in A, see:   
G.W. Stewart (1977). “Perturbation Bounds for the QR Factorization of a Matrix,” SIAM J. Numer. Anal. 14, 509–518.   
H. Zha (1993). “A Componentwise Perturbation Analysis of the QR Decomposition,” SIAM J. Matrix Anal. Applic. 4, 1124–1131.   
G.W. Stewart (1993). “On the Perturbation of LU Cholesky, and QR Factorizations,” SIAM J. Matrix Anal. Applic. 14, 1141–1145.   
A. Barrlund (1994). “Perturbation Bounds for the Generalized QR Factorization,” Lin. Alg. Applic. 207, 251–271.   
J.-G. Sun (1995). “On Perturbation Bounds for the QR Factorization,” Lin. Alg. Applic. 215, 95–112.   
X.-W. Chang and C.C. Paige (2001). “Componentwise Perturbation Analyses for the QR factorization,” Numer. Math. 88, 319–345.   
Organization of the computation so that the entries in Q depend continuously on the entries in A is discussed in:   
T.F. Coleman and D.C. Sorensen (1984). “A Note on the Computation of an Orthonormal Basis for the Null Space of a Matrix,” Mathematical Programming 29, 234–242.   
References for the Gram-Schmidt process and various ways to overcome its shortfalls include:   
J.R. Rice (1966). “Experiments on Gram-Schmidt Orthogonalization,” Math. Comput. 20, 325–328.   
A. Bj¨orck (1967). “Solving Linear Least Squares Problems by Gram-Schmidt Orthogonalization,” BIT 7, 1–21.   
N.N. Abdelmalek (1971). “Roundoff Error Analysis for Gram-Schmidt Method and Solution of Linear Least Squares Problems,” BIT 11, 345–368.   
A. Ruhe (1983). “Numerical Aspects of Gram-Schmidt Orthogonalization of Vectors,” Lin. Alg. Applic. 52/53, 591–601.   
W. Jalby and B. Philippe (1991). “Stability Analysis and Improvement of the Block Gram-Schmidt Algorithm,” SIAM J. Sci. Stat. Comput. 12, 1058–1073.   
A. Bj¨ ˚ orck and C.C. Paige (1992). “Loss and Recapture of Orthogonality in the Modified Gram-Schmidt Algorithm,” SIAM J. Matrix Anal. Applic. 13, 176–190.   
A. Bj¨orck (1994). “Numerics of Gram-Schmidt Orthogonalization,” Lin. Alg. Applic. 197/198, 297–316.   
L. Giraud and J. Langou (2003). “A Robust Criterion for the Modified Gram-Schmidt Algorithm with Selective Reorthogonalization,” SIAM J. Sci. Comput. 25, 417–441.   
G.W. Stewart (2005). “Error Analysis of the Quasi-Gram–Schmidt Algorithm,” SIAM J. Matrix Anal. Applic. 27, 493–506.

L. Giraud, J. Langou, M. Rozlonk, and J. van den Eshof (2005). “Rounding Error Analysis of the Classical Gram-Schmidt Orthogonalization Process,” Numer. Math. 101, 87–100.   
A. Smoktunowicz, J.L. Barlow and J. Langou (2006). “A Note on the Error Analysis of Classical Gram-Schmidt,” Numer. Math. 105, 299–313.   
Various high-performance issues pertaining to the QR factorization are discussed in:   
B. Mattingly, C. Meyer, and J. Ortega (1989). “Orthogonal Reduction on Vector Computers,” SIAM J. Sci. Stat. Comput. 10, 372–381.   
P.A. Knight (1995). “Fast Rectangular Matrix Multiplication and the QR Decomposition,” Lin. Alg. Applic. 221, 69–81.   
J.J. Carrig, Jr. and G.L. Meyer (1997). “Efficient Householder QR Factorization for Superscalar Processors,” ACM Trans. Math. Softw. 23, 362–378.   
D. Vanderstraeten (2000). “An Accurate Parallel Block Gram-Schmidt Algorithm without Reorthogonalization,” Numer. Lin. Alg. 7, 219–236.   
E. Elmroth and F.G. Gustavson (2000). “Applying Recursion to Serial and Parallel QR Factorization Leads to Better Performance,” IBM J. Res. Dev. 44, 605–624.   
Many important high-performance implementation ideas apply equally to LU, Cholesky, and QR, see:   
A. Buttari, J. Langou, J. Kurzak, and J. Dongarra (2009). “A Class of Parallel Tiled Linear Algebra Algorithms for Multicore Architectures,” Parallel Comput. 35, 38–53.   
J. Kurzak, H. Ltaief, and J. Dongarra (2010). “Scheduling Dense Linear Algebra Operations on Multicore Processors,” Concurrency Comput. Pract. Exper. 22, 15–44.   
J. Demmel, L. Grigori, M, Hoemmen, and J. Langou (2012). “Methods and Algorithms for Scientific Computing Communication-optimal Parallel and Sequential QR and LU Factorizations,” SIAM J. Sci. Comput. 34, A206-A239.   
Historical references concerned with parallel Givens QR include:   
W.M. Gentleman and H.T. Kung (1981). “Matrix Triangularization by Systolic Arrays,” SPIE Proc. 298, 19–26.   
D.E. Heller and I.C.F. Ipsen (1983). “Systolic Networks for Orthogonal Decompositions,” SIAM J. Sci. Stat. Comput. 4, 261–269.   
M. Costnard, J.M. Muller, and Y. Robert (1986). “Parallel QR Decomposition of a Rectangular Matrix,” Numer. Math. 48, 239–250.   
L. Eldin and R. Schreiber (1986). “An Application of Systolic Arrays to Linear Discrete Ill-Posed Problems,” SIAM J. Sci. Stat. Comput. 7, 892–903.   
F.T. Luk (1986). “A Rotation Method for Computing the QR Factorization,” SIAM J. Sci. Stat. Comput. 7, 452–459.   
J.J. Modi and M.R.B. Clarke (1986). “An Alternative Givens Ordering,” Numer. Math. 43, 83–90.   
The QR factorization of a structured matrix is usually structured itself, see:   
A.W. Bojanczyk, R.P. Brent, and F.R. de Hoog (1986). “QR Factorization of Toeplitz Matrices,” Numer. Math. 49, 81–94.   
S. Qiao (1986). “Hybrid Algorithm for Fast Toeplitz Orthogonalization,” Numer. Math. 53, 351–366.   
C.J. Demeure (1989). “Fast QR Factorization of Vandermonde Matrices,” Lin. Alg. Applic. 122/123/124, 165–194.   
L. Reichel (1991). “Fast QR Decomposition of Vandermonde-Like Matrices and Polynomial Least Squares Approximation,” SIAM J. Matrix Anal. Applic. 12, 552–564.   
D.R. Sweet (1991). “Fast Block Toeplitz Orthogonalization,” Numer. Math. 58, 613–629.   
Quantum computation has an interesting connection to complex Givens rotations and their application to vectors, see:   
G. Cybenko (2001). “Reducing Quantum Computations to Elementary Unitary Transformations,” Comput. Sci. Eng. 3, 27–32.   
D.P. O’Leary and S.S. Bullock (2005). “QR Factorizations Using a Restricted Set of Rotations,” ETNA 21, 20–27.   
N.D. Mermin (2007). Quantum Computer Science, Cambridge University Press, New York.

# 5.3 The Full-Rank Least Squares Problem

Consider the problem of finding a vector $\boldsymbol { x } \in \mathbb { R } ^ { n }$ such that $A x = b$ where the data matrix $A \in \mathbb { R } ^ { m \times n }$ and the observation vector $b \in \mathbb { R } ^ { m }$ are given and $m \geq n$ . When there are more equations than unknowns, we say that the system $A x = b$ is overdetermined. Usually an overdetermined system has no exact solution since b must be an element of ran(A), a proper subspace of $\mathbb { R } ^ { m }$ .

This suggests that we strive to minimize $\| A x - b \| _ { p }$ for some suitable choice of $p .$ Different norms render different optimum solutions. For example, if $A = [ 1 , 1 , 1 ] ^ { T }$ and $\boldsymbol { b } = [ b _ { 1 } , b _ { 2 } , b _ { 3 } ] ^ { T }$ with $b _ { 1 } \geq b _ { 2 } \geq b _ { 3 } \geq 0$ , then it can be verified that

$$
p = 1 \Rightarrow x _ {\text {opt}} = b _ {2},
$$

$$
p = 2 \Rightarrow x _ {\text {opt}} = (b _ {1} + b _ {2} + b _ {3}) / 3,
$$

$$
p = \infty \Rightarrow x _ {\text {opt}} = (b _ {1} + b _ {3}) / 2.
$$

Minimization in the 1-norm and infinity-norm is complicated by the fact that the function $f ( x ) = \parallel A x - b \parallel _ { p }$ is not differentiable for these values of $p .$ However, there are several good techniques available for 1-norm and ∞-norm minimization. See Coleman and Li (1992), Li (1993), and Zhang (1993).

In contrast to general p-norm minimization, the least squares (LS) problem

$$
\min _ {x \in \mathbb {R} ^ {n}} \| A x - b \| _ {2} \tag {5.3.1}
$$

is more tractable for two reasons:

• $\begin{array} { r } { \phi ( x ) = \frac { 1 } { 2 } \| A x - b \| _ { 2 } ^ { 2 } } \end{array}$ is a differentiable function of x and so the minimizers of $\phi$ satisfy the gradient equation $\nabla \phi ( x ) = 0$ . This turns out to be an easily constructed symmetric linear system which is positive definite if A has full column rank.   
• The 2-norm is preserved under orthogonal transformation. This means that we can seek an orthogonal Q such that the equivalent problem of minimizing $\parallel ( Q ^ { T } A ) x - ( Q ^ { T } b ) \parallel _ { 2 }$ is “easy” to solve.

In this section we pursue these two solution approaches for the case when A has full column rank. Methods based on normal equations and the QR factorization are detailed and compared.

# 5.3.1 Implications of Full Rank

Suppose $x \in \mathbb { R } ^ { n } , z \in \mathbb { R } ^ { n } , \alpha \in \mathbb { R }$ , and consider the equality

$$
\| A (x + \alpha z) - b \| _ {2} ^ {2} = \| A x - b \| _ {2} ^ {2} + 2 \alpha z ^ {T} A ^ {T} (A x - b) + \alpha^ {2} \| A z \| _ {2} ^ {2}
$$

where $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ . If x solves the LS problem (5.3.1), then we must have $A ^ { T } ( A x - b ) = 0 .$ . Otherwise, i $\mathrm { ~ f ~ } z = - A ^ { T } ( A x - b )$ and we make α small enough, then we obtain the contradictory inequality $\parallel A ( x + \alpha z ) - b \parallel _ { 2 } < \parallel A x - b \parallel _ { 2 }$ . We may also conclude that if x and $x + \alpha z$ are LS minimizers, then $z \in \mathsf { n u l l } ( A )$ .

Thus, if A has full column rank, then there is a unique LS solution $x _ { \mathrm { L S } }$ and it solves the symmetric positive definite linear system

$$
A ^ {T} A x _ {\mathrm{LS}} = A ^ {T} b.
$$

These are called the normal equations. Note that if

$$
\phi (x) = \frac {1}{2} \| A x - b \| _ {2} ^ {2},
$$

then

$$
\nabla \phi (x) = A ^ {T} (A x - b),
$$

so solving the normal equations is tantamount to solving the gradient equation $\nabla \phi = 0$ . We call

$$
r _ {\mathrm{LS}} = b - A x _ {\mathrm{LS}}
$$

the minimum residual and we use the notation

$$
\rho_ {\mathrm{LS}} = \left\| A x _ {\mathrm{LS}} - b \right\| _ {2}
$$

to denote its size. Note that if $\rho _ { \mathrm { L S } }$ is small, then we can do a good job “predicting” b by using the columns of A.

Thus far we have been assuming that $A \in \mathbb { R } ^ { m \times n }$ has full column rank, an assumption that is dropped in §5.5. However, even if rank(A) = n, trouble can be expected if A is nearly rank deficient. The SVD can be used to substantiate this remark. If

$$
A = U \Sigma V ^ {T} = \sum_ {i = 1} ^ {n} \sigma_ {i} u _ {i} v _ {i} ^ {T}
$$

is the SVD of a full rank matrix $A \in \mathbb { R } ^ { m \times n }$ , then

$$
\parallel A x - b \parallel_ {2} ^ {2} = \parallel (U ^ {T} A V) (V ^ {T} x) - U ^ {T} b \parallel_ {2} ^ {2} = \sum_ {i = 1} ^ {n} (\sigma_ {i} y _ {i} - (u _ {i} ^ {T} b)) ^ {2} + \sum_ {i = n + 1} ^ {m} (u _ {i} ^ {T} b) ^ {2}
$$

where $y = V ^ { T } x$ . It follows that this summation is minimized by setting $y _ { i } = u _ { i } ^ { T } b / \sigma _ { i }$ , $i = 1 { : } n$ . Thus,

$$
x _ {\mathrm{LS}} = \sum_ {i = 1} ^ {n} \frac {u _ {i} ^ {T} b}{\sigma_ {i}} v _ {i} \tag {5.3.2}
$$

and

$$
\rho_ {\mathrm{LS}} ^ {2} = \sum_ {i = n + 1} ^ {2} (u _ {i} ^ {T} b) ^ {2}. \tag {5.3.3}
$$

It is clear that the presence of small singular values means LS solution sensitivity. The effect of perturbations on the minimum sum of squares is less clear and requires further analysis which we offer below.

When assessing the quality of a computed LS solution $\hat { x } _ { \mathrm { L S } }$ , there are two important issues to bear in mind:

• How close is $\hat { x } _ { \mathrm { L S } }$ to $x _ { \mathrm { L S } } \mathrm { ^ { 2 } }$   
• How small is $\hat { r } _ { \mathrm { L S } } = b - A \hat { x } _ { \mathrm { L S } }$ compared to $r _ { \mathrm { L S } } = b - A x _ { \mathrm { L S } } ?$

The relative importance of these two criteria varies from application to application. In any case it is important to understand how $x _ { \mathrm { L S } }$ and $r _ { \mathrm { L S } }$ are affected by perturbations in A and b. Our intuition tells us that if the columns of A are nearly dependent, then these quantities may be quite sensitive. For example, suppose

$$
A = \left[ \begin{array}{l l} 1 & 0 \\ 0 & 1 0 ^ {- 6} \\ 0 & 0 \end{array} \right], \delta A = \left[ \begin{array}{l l} 0 & 0 \\ 0 & 0 \\ 0 & 1 0 ^ {- 8} \end{array} \right], b = \left[ \begin{array}{l} 1 \\ 0 \\ 1 \end{array} \right], \delta b = \left[ \begin{array}{l} 0 \\ 0 \\ 0 \end{array} \right],
$$

and that $x _ { \mathrm { L S } }$ and $\hat { x } _ { \mathrm { L S } }$ minimize $\parallel A x - b \parallel _ { 2 }$ and $\parallel ( A + \delta A ) x - ( b + \delta b ) \parallel _ { 2 } :$ , respectively. If $r _ { \mathrm { L S } }$ and $\hat { r } _ { \mathrm { L S } }$ are the corresponding minimum residuals, then it can be shown that

$$
x _ {\mathrm{LS}} = \left[ \begin{array}{c} 1 \\ 0 \end{array} \right], \hat {x} _ {\mathrm{LS}} = \left[ \begin{array}{c} 1 \\ . 9 9 9 9 \cdot 1 0 ^ {4} \end{array} \right], r _ {\mathrm{LS}} = \left[ \begin{array}{c} 0 \\ 0 \\ 1 \end{array} \right], \hat {r} _ {\mathrm{LS}} = \left[ \begin{array}{c} 0 \\ -. 9 9 9 9 \cdot 1 0 ^ {- 2} \\ . 9 9 9 9 \cdot 1 0 ^ {0} \end{array} \right].
$$

Recall that the 2-norm condition of a rectangular matrix is the ratio of its largest to smallest singular values. Since $\kappa _ { 2 } ( A ) = 1 0 ^ { 6 }$ we have

$$
\frac {\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \| _ {2}}{\| x _ {\mathrm{LS}} \| _ {2}} \approx . 9 9 9 9 \cdot 1 0 ^ {4} \leq \kappa_ {2} (A) ^ {2} \frac {\| \delta A \| _ {2}}{\| A \| _ {2}} = 1 0 ^ {1 2} \cdot 1 0 ^ {- 8}
$$

and

$$
\frac {\| \hat {r} _ {\mathrm{LS}} - r _ {\mathrm{LS}} \| _ {2}}{\| b \| _ {2}} \approx . 7 0 7 0 \cdot 1 0 ^ {- 2} \leq \kappa_ {2} (A) \frac {\| \delta A \| _ {2}}{\| A \| _ {2}} = 1 0 ^ {6} \cdot 1 0 ^ {- 8}.
$$

The example suggests that the sensitivity of $x _ { \mathrm { L S } }$ can depend upon $\kappa _ { 2 } ( A ) ^ { 2 }$ . Below we offer an LS perturbation theory that confirms the possibility.

# 5.3.2 The Method of Normal Equations

A widely-used method for solving the full-rank LS problem is the method of normal equations.

Algorithm 5.3.1 (Normal Equations) Given $A \in \mathbb { R } ^ { m \times n }$ with the property that rank(A) = n and $b \in \mathbb { R } ^ { m }$ , this algorithm computes a vector $x _ { \mathrm { L S } }$ that minimizes $\parallel A x - b \parallel _ { 2 }$ .

Compute the lower triangular portion of $C = A ^ { T } A$ .

Form the matrix-vector product $d = A ^ { T } b$ .

Compute the Cholesky factorization $C = G G ^ { T }$ .

Solve $G y = d$ and $G ^ { T } x _ { \mathrm { L S } } = y .$

This algorithm requires $( m + n / 3 ) n ^ { 2 }$ flops. The normal equation approach is convenient because it relies on standard algorithms: Cholesky factorization, matrix-matrix multiplication, and matrix-vector multiplication. The compression of the m-by-n data matrix A into the (typically) much smaller n-by-n cross-product matrix C is attractive.

Let us consider the accuracy of the computed normal equations solution $\hat { x } _ { \mathrm { L S } }$ . For clarity, assume that no roundoff errors occur during the formation of $C = A ^ { T } A$ and $d = A ^ { T } b$ . It follows from what we know about the roundoff properties of the Cholesky factorization (§4.2.6) that

$$
(A ^ {T} A + E) \hat {x} _ {\mathrm{LS}} = A ^ {T} b
$$

where

$$
\| E \| _ {2} \approx \mathbf {u} \| A ^ {T} \| _ {2} \| A \| _ {2} = \mathbf {u} \| A ^ {T} A \| _ {2}.
$$

Thus, we can expect

$$
\frac {\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \| _ {2}}{\| x _ {\mathrm{LS}} \| _ {2}} \approx \mathbf {u} \kappa_ {2} (A ^ {T} A) = \mathbf {u} \kappa_ {2} (A) ^ {2}. \tag {5.3.4}
$$

In other words, the accuracy of the computed normal equations solution depends on the square of the condition. See Higham (ASNA, §20.4) for a detailed roundoff analysis of the normal equations approach.

It should be noted that the formation of $A ^ { T } A$ can result in a significant loss of information. If

$$
A = \left[ \begin{array}{c c} 1 & 1 \\ \sqrt {\mathbf {u}} & 0 \\ 0 & \sqrt {\mathbf {u}} \end{array} \right],
$$

then $\kappa _ { 2 } ( A ) \approx \sqrt { \bf { u } }$ . However,

$$
\mathsf {f l} (A ^ {T} A) = \left[ \begin{array}{l l} 1 & 1 \\ 1 & 1 \end{array} \right]
$$

is exactly singular. Thus, the method of normal equations can break down on matrices that are not particularly close to being numerically rank deficient.

# 5.3.3 LS Solution Via QR Factorization

Let $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ and $b \in \mathbb { R } ^ { m }$ be given and suppose that an orthogonal matrix $Q \in \mathbb { R } ^ { m \times m }$ has been computed such that

$$
Q ^ {T} A = R = \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right] _ {m - n} ^ {n} \tag {5.3.5}
$$

is upper triangular. If

$$
Q ^ {T} b = \left[ \begin{array}{l} c \\ d \end{array} \right] _ {m - n} ^ {n}
$$

then

$$
\parallel A x - b \parallel_ {2} ^ {2} = \parallel Q ^ {T} A x - Q ^ {T} b \parallel_ {2} ^ {2} = \parallel R _ {1} x - c \parallel_ {2} ^ {2} + \parallel d \parallel_ {2} ^ {2}
$$

for any $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . Since $\mathsf { r a n k } ( A ) = \mathsf { r a n k } ( R _ { 1 } ) = n$ , it follows that $x _ { \mathrm { L S } }$ is defined by the upper triangular system

$$
R _ {1} x _ {\mathrm{LS}} = c.
$$

Note that

$$
\rho_ {\mathrm{LS}} = \parallel d \parallel_ {2}.
$$

We conclude that the full-rank LS problem can be readily solved once we have computed the QR factorization of A. Details depend on the exact QR procedure. If Householder matrices are used and $Q ^ { T }$ is applied in factored form to b, then we obtain

Algorithm 5.3.2 (Householder LS Solution) If $A \in \mathbb { R } ^ { m \times n }$ has full column rank and $b \in \mathbb { R } ^ { m }$ , then the following algorithm computes a vector $\boldsymbol { x } _ { \mathrm { L S } } \in \mathbb { R } ^ { n }$ such that $\parallel A x _ { \mathrm { L S } } - b \parallel _ { 2 }$ is minimum.

Use Algorithm 5.2.1 to overwrite A with its QR factorization.

for $j = 1 { : } n$

$$
v = \left[ \begin{array}{c} 1 \\ A (j + 1: m, j) \end{array} \right]
$$

$$
\beta = 2 / v ^ {T} v
$$

$$
b (j: m) = b (j: m) - \beta (v ^ {T} b (j: m)) v
$$

end

Solve $R ( 1 : n , 1 : n ) \cdot x _ { \mathrm { { L S } } } = b ( 1 { : } n )$ .

This method for solving the full-rank LS problem requires $2 n ^ { 2 } ( m - n / 3 )$ flops. The O(mn) flops associated with the updating of b and the $O ( n ^ { 2 } )$ flops associated with the back substitution are not significant compared to the work required to factor A.

It can be shown that the computed $\hat { x } _ { \mathrm { L S } }$ solves

$$
\min \| (A + \delta A) x - (b + \delta b) \| _ {2} \tag {5.3.6}
$$

where

$$
\left\| \delta A \right\| _ {F} \leq (6 m - 3 n + 4 1) n \mathbf {u} \left\| A \right\| _ {F} + O \left(\mathbf {u} ^ {2}\right) \tag {5.3.7}
$$

and

$$
\| \delta b \| _ {2} \leq (6 m - 3 n + 4 0) n \mathbf {u} \| b \| _ {2} + O (\mathbf {u} ^ {2}). \tag {5.3.8}
$$

These inequalities are established in Lawson and Hanson (SLS, p. 90ff) and show that $\hat { x } _ { \mathrm { L S } }$ satisfies a “nearby” LS problem. (We cannot address the relative error in $\hat { x } _ { \mathrm { L S } }$ without an LS perturbation theory, to be discussed shortly.) We mention that similar results hold if Givens QR is used.

# 5.3.4 Breakdown in Near-Rank-Deficient Case

As with the method of normal equations, the Householder method for solving the LS problem breaks down in the back-substitution phase if rank $( A ) ~ < ~ n$ . Numerically, trouble can be expected if $\kappa _ { 2 } ( A ) = \kappa _ { 2 } ( R ) \approx 1 / \mathbf { u }$ . This is in contrast to the normal equations approach, where completion of the Cholesky factorization becomes problematical once $\kappa _ { 2 } ( A )$ is in the neighborhood of $1 / \sqrt { \mathbf { u } }$ as we showed above. Hence the claim in Lawson and Hanson (SLS, pp. 126–127) that for a fixed machine precision, a wider class of LS problems can be solved using Householder orthogonalization.

# 5.3.5 A Note on the MGS Approach

In principle, MGS computes the thin QR factorization $A = Q _ { 1 } R _ { 1 }$ . This is enough to solve the full-rank LS problem because it transforms the normal equation system $( A ^ { T } A ) x = A ^ { T } b$ to the upper triangular system $R _ { 1 } x = Q _ { 1 } ^ { T } b$ . But an analysis of this approach when $Q _ { 1 } ^ { T } b$ is explicitly formed introduces a $\kappa _ { 2 } ( A ) ^ { 2 }$ term. This is because the computed factor $\hat { Q } _ { 1 }$ satisfies $\Vert \hat { Q } _ { 1 } ^ { T } \hat { Q } _ { 1 } - I _ { n } \Vert _ { 2 } \approx \mathbf { u } \kappa _ { 2 } ( A )$ as we mentioned in §5.2.9.

However, if MGS is applied to the augmented matrix

$$
A _ {+} = \left[ \begin{array}{c c} A & b \end{array} \right] = \left[ \begin{array}{c c} Q _ {1} & q _ {n + 1} \end{array} \right] \left[ \begin{array}{c c} R _ {1} & z \\ 0 & \rho \end{array} \right],
$$

then $z = Q _ { 1 } ^ { T } b$ . Computing $Q _ { 1 } ^ { T } b$ in this fashion and solving $R _ { 1 } x _ { \mathrm { L S } } = z$ produces an LS solution $\hat { x } _ { \mathrm { L S } }$ that is “just as good” as the Householder QR method. That is to say, a result of the form (5.3.6)–(5.3.8) applies. See Bj¨orck and Paige (1992).

It should be noted that the MGS method is slightly more expensive than Householder QR because it always manipulates m-vectors whereas the latter procedure deals with vectors that become shorter in length as the algorithm progresses.

# 5.3.6 The Sensitivity of the LS Problem

We now develop a perturbation theory for the full-rank LS problem that assists in the comparison of the normal equations and QR approaches. LS sensitivity analysis has a long and fascinating history. Grcar (2009, 2010) compares about a dozen different results that have appeared in the literature over the decades and the theorem below follows his analysis. It examines how the LS solution and its residual are affected by changes in A and b and thereby sheds light on the condition of the LS problem. Four facts about $A \in \mathbb { R } ^ { m \times n }$ are used in the proof, where it is assumed that $m > n \colon$

$$
1 = \| A (A ^ {T} A) ^ {- 1} A ^ {T} \| _ {2}, \quad \frac {1}{\sigma_ {n} (A)} = \| (A ^ {T} A) ^ {- 1} A ^ {T} \| _ {2}, \tag {5.3.9}
$$

$$
1 = \| I - A (A ^ {T} A) ^ {- 1} A ^ {T} \| _ {2}, \quad \frac {1}{\sigma_ {n} (A) ^ {2}} = \| (A ^ {T} A) ^ {- 1} \| _ {2}.
$$

These equations are easily verified using the SVD.

Theorem 5.3.1. Suppose that $x _ { \mathrm { L S } } , \ r _ { \mathrm { L S } } , \ \hat { x } _ { \mathrm { L S } }$ , and $\hat { r } _ { \mathrm { L S } }$ satisfy

$$
\left\| A x _ {\mathrm{LS}} - b \right\| _ {2} = \min, \quad r _ {\mathrm{LS}} = b - A x _ {\mathrm{LS}},
$$

$$
\| (A + \delta A) \hat {x} _ {\mathrm{LS}} - (b + \delta b) \| _ {2} = \min, \quad \hat {r} _ {\mathrm{LS}} = (b + \delta b) - (A + \delta A) \hat {x} _ {\mathrm{LS}},
$$

where A has rank n and $\parallel \delta A \parallel _ { 2 } < \sigma _ { n } ( A )$ . Assume that b, $r _ { \mathrm { L S } }$ , and $x _ { \mathrm { L S } }$ are not zero. Let $\theta _ { \mathrm { L S } } \in ( 0 , \pi / 2 )$ be defined by

$$
\sin (\theta_ {\mathrm{LS}}) = \frac {\parallel r _ {\mathrm{LS}} \parallel_ {2}}{\parallel b \parallel_ {2}}.
$$

If

$$
\epsilon = \max \left\{\frac {\parallel \delta A \parallel_ {2}}{\parallel A \parallel_ {2}}, \frac {\parallel \delta b \parallel_ {2}}{\parallel b \parallel_ {2}} \right\}
$$

and

$$
\nu_ {\mathrm{LS}} = \frac {\| A x _ {\mathrm{LS}} \| _ {2}}{\sigma_ {n} (A) \| x _ {\mathrm{LS}} \| _ {2}}, \tag {5.3.10}
$$

then

$$
\frac {\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \| _ {2}}{\| x \| _ {2}} \leq \epsilon \left\{\frac {v _ {\mathrm{LS}}}{\cos (\theta_ {\mathrm{LS}})} + [ 1 + \nu_ {\mathrm{LS}} \tan (\theta_ {\mathrm{LS}}) ] \kappa_ {2} (A) \right\} + O (\epsilon^ {2}) \tag {5.3.11}
$$

and

$$
\frac {\| \hat {r} _ {\mathrm{LS}} - r _ {\mathrm{LS}} \| _ {2}}{\| r _ {\mathrm{LS}} \| _ {2}} \leq \epsilon \left\{\frac {1}{\sin (\theta_ {\mathrm{LS}})} + \left[ \frac {1}{\nu_ {\mathrm{LS}} \tan (\theta_ {\mathrm{LS}})} + 1 \right] \kappa_ {2} (A) \right\} + O (\epsilon^ {2}). \tag {5.3.12}
$$

Proof. Let E and f be defined by $E = \delta A / \epsilon$ and $f = \delta b / \epsilon$ . By Theorem 2.5.2 we have rank $( A + t E ) = n$ for all $t \in [ 0 , \epsilon ]$ . It follows that the solution $x ( t )$ t o

$$
(A + t E) ^ {T} (A + t E) x (t) = (A + t E) ^ {T} (b + t f) \tag {5.3.13}
$$

is continuously differentiable for all $t \in [ 0 , \epsilon ]$ . Since $x _ { \mathrm { L S } } = x ( 0 )$ and $\hat { x } _ { \mathrm { L S } } = x ( \epsilon )$ , we have

$$
\hat {x} _ {\mathrm{LS}} = x _ {\mathrm{LS}} + \epsilon \dot {x} (0) + O (\epsilon^ {2}).
$$

By taking norms and dividing by $\Vert \boldsymbol { x } _ { \mathrm { L S } } \Vert _ { 2 }$ we obtain

$$
\frac {\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \| _ {2}}{\| x _ {\mathrm{LS}} \| _ {2}} = \epsilon \frac {\| \dot {x} (0) \| _ {2}}{\| x _ {\mathrm{LS}} \| _ {2}} + O (\epsilon^ {2}). \tag {5.3.14}
$$

In order to bound $\Vert \dot { x } ( 0 ) \Vert _ { 2 }$ , we differentiate (5.3.13) and set t = 0 in the result. This gives

$$
E ^ {T} A x _ {\mathrm{LS}} + A ^ {T} E x _ {\mathrm{LS}} + A ^ {T} A \dot {x} (0) = A ^ {T} f + E ^ {T} b,
$$

i.e.,

$$
\dot {x} (0) = (A ^ {T} A) ^ {- 1} A ^ {T} (f - E x _ {\mathrm{LS}}) + (A ^ {T} A) ^ {- 1} E ^ {T} r _ {\mathrm{LS}}. \tag {5.3.15}
$$

Using (5.3.9) and the inequalities $\left\| \ f \right\| _ { 2 } \leq \left\| \ b \right\| _ { 2 }$ and $\left. \ E \right. _ { 2 } \leq \left. \ A \right. _ { 2 }$ , it follows that

$$
\| \dot {x} (0) \| \leq \| (A ^ {T} A) ^ {- 1} A ^ {T} f \| _ {2} + \| (A ^ {T} A) ^ {- 1} A ^ {T} E x _ {\mathrm{LS}} \| _ {2} + \| (A ^ {T} A) ^ {- 1} E ^ {T} r _ {\mathrm{LS}} \| _ {2}
$$

$$
\leq \frac {\| b \| _ {2}}{\sigma_ {n} (A)} + \frac {\| A \| _ {2} \| x _ {\mathrm{LS}} \| _ {2}}{\sigma_ {n} (A)} + \frac {\| A \| _ {2} \| r _ {\mathrm{LS}} \| _ {2}}{\sigma_ {n} (A) ^ {2}}.
$$

By substituting this into (5.3.14) we obtain

$$
\frac {\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \| _ {2}}{\| x _ {\mathrm{LS}} \| _ {2}} \leq \epsilon \left(\frac {\| b \| _ {2}}{\sigma_ {n} (A) \| x _ {\mathrm{LS}} \| _ {2}} + \frac {\| A \| _ {2}}{\sigma_ {n} (A)} + \frac {\| A \| _ {2} \| r _ {\mathrm{LS}} \| _ {2}}{\sigma_ {n} (A) ^ {2} \| x _ {\mathrm{LS}} \| _ {2}}\right) + O (\epsilon^ {2}).
$$

Inequality (5.3.11) follows from the definitions of $\kappa _ { 2 } ( A )$ and $\nu _ { \mathrm { L S } }$ and the identities

$$
\cos (\theta_ {\mathrm{LS}}) = \frac {\parallel A x _ {\mathrm{LS}} \parallel_ {2}}{\parallel b \parallel_ {2}}, \quad \tan (\theta_ {\mathrm{LS}}) = \frac {\parallel r _ {\mathrm{LS}} \parallel_ {2}}{\parallel A x _ {\mathrm{LS}} \parallel_ {2}}. \tag {5.3.16}
$$

The proof of the residual bound (5.3.12) is similar. Define the differentiable vector function r(t) by

$$
r (t) = (b + t f) - (A + t E) x (t)
$$

and observe that $r _ { \mathrm { L S } } = r ( 0 )$ and $\hat { r } _ { \mathrm { L S } } = r ( \epsilon )$ . Thus,

$$
\frac {\left\| \hat {r} _ {\mathrm{LS}} - r _ {\mathrm{LS}} \right\| _ {2}}{\left\| r _ {\mathrm{LS}} \right\| _ {2}} = \epsilon \frac {\left\| \dot {r} (0) \right\| _ {2}}{\left\| r _ {\mathrm{LS}} \right\| _ {2}} + O \left(\epsilon^ {2}\right). \tag {5.3.17}
$$

From (5.3.15) we have

$$
\dot {r} (0) = \left(I - A (A ^ {T} A) ^ {- 1} A ^ {T}\right) (f - E x _ {\mathrm{LS}}) - A (A ^ {T} A) ^ {- 1} E ^ {T} r _ {\mathrm{LS}}.
$$

By taking norms, using (5.3.9) and the inequalities  $f \parallel _ { 2 } \leq \parallel b \parallel _ { 2 }$ and $\| E \| _ { 2 } \leq \| A \| _ { 2 }$ we obtain

$$
\| \dot {r} (0) \| _ {2} \leq \| b \| _ {2} + \| A \| _ {2} \| x _ {\mathrm{LS}} \| _ {2} + \frac {\| A \| _ {2} \| r _ {\mathrm{LS}} \| _ {2}}{\sigma_ {n} (A)}
$$

and thus from (5.3.17) we have

$$
\frac {\| \hat {r} _ {\mathrm{LS}} - r _ {\mathrm{LS}} \| _ {2}}{\| r _ {\mathrm{LS}} \| _ {2}} \leq \frac {\| b \| _ {2}}{\| r _ {\mathrm{LS}} \| _ {2}} + \frac {\| A \| _ {2} \| x _ {\mathrm{LS}} \| _ {2}}{\| r _ {\mathrm{LS}} \| _ {2}} + \frac {\| A \| _ {2}}{\sigma_ {n} (A)}.
$$

The inequality (5.3.12) follows from the definitions of $\kappa _ { 2 } ( A )$ and $\nu _ { \mathrm { L S } }$ and the identities (5.3.16).

It is instructive to identify conditions that turn the upper bound in (5.3.11) into a bound that involves $\kappa _ { 2 } ( A ) ^ { 2 }$ . The example in §5.3.1 suggests that this factor might figure in the definition of an LS condition number. However, the theorem shows that the situation is more subtle. Note that

$$
\nu_ {\mathrm{LS}} = \frac {\parallel A x _ {\mathrm{LS}} \parallel_ {2}}{\sigma_ {n} (A) \parallel x _ {\mathrm{LS}} \parallel_ {2}} \leq \frac {\parallel A \parallel_ {2}}{\sigma_ {n} (A)} = \kappa_ {2} (A).
$$

The SVD expansion (5.3.2) suggests that if b has a modest component in the direction of the left singular vector $u _ { n }$ , then

$$
\nu_ {\mathrm{LS}} \approx \kappa_ {2} (A).
$$

If this is the case and $\theta _ { \mathrm { L S } }$ is sufficiently bounded away from $\pi / 2$ , then the inequality (5.3.11) essentially says that

$$
\frac {\left\| \hat {x} _ {\mathrm{LS}} - x _ {\mathrm{LS}} \right\| _ {2}}{\left\| x _ {\mathrm{LS}} \right\| _ {2}} \approx \epsilon \left(\kappa_ {2} (A) + \frac {\rho_ {\mathrm{LS}}}{\left\| b \right\| _ {2}} \kappa_ {2} (A) ^ {2}\right). \tag {5.3.18}
$$

Although this simple heuristic assessment of LS sensitivity is almost always applicable, it important to remember that the true condition of a particular LS problem depends on $\nu _ { \mathrm { L S } } , \theta _ { \mathrm { L S } }$ , and $\kappa _ { 2 } ( A )$ .

Regarding the perturbation of the residual, observe that the upper bound in the residual result (5.3.12) is less than the upper bound in the solution result (5.3.11) by a factor of $\nu _ { \mathrm { L S } } \tan ( \theta _ { \mathrm { L S } } )$ . We also observe that if $\theta _ { \mathrm { L S } }$ is sufficiently bounded away from both 0 and $\pi / 2$ , then (5.3.12) essentially says that

$$
\frac {\left\| \hat {r} _ {\mathrm{LS}} - r _ {\mathrm{LS}} \right\| _ {2}}{\left\| r _ {\mathrm{LS}} \right\| _ {2}} \approx \epsilon \cdot \kappa_ {2} (A). \tag {5.3.19}
$$

For more insights into the subtleties behind Theorem 5.3.1., see Wedin (1973), Vandersluis (1975), Bj¨orck (NMLS, p. 30), Higham (ASNA, p. 382), and Grcar(2010).

# 5.3.7 Normal Equations Versus QR

It is instructive to compare the normal equation and QR approaches to the full-rank LS problem in light of Theorem 5.3.1.

• The method of normal equations produces an $\hat { x } _ { \mathrm { L S } }$ whose relative error depends on $\kappa _ { 2 } ( A ) ^ { 2 }$ , a factor that can be considerably larger than the condition number associated with a “small residual” LS problem.   
• The QR approach (Householder, Givens, careful MGS) solves a nearby LS problem. Therefore, these methods produce a computed solution with relative error that is “predicted” by the condition of the underlying LS problem.

Thus, the QR approach is more appealing in situations where b is close to the span of A’s columns.

Finally, we mention two other factors that figure in the debate about QR versus normal equations. First, the normal equations approach involves about half of the arithmetic when m 
 n and does not require as much storage, assuming that $Q ( : , 1 : n )$ is required. Second, QR approaches are applicable to a wider class of LS problems. This is because the Cholesky solve in the method of normal equations is “in trouble” if $\kappa _ { 2 } ( A ) \approx 1 / \sqrt { \mathbf { u } }$ while the R-solve step in a QR approach is in trouble only if $\kappa _ { 2 } ( A ) \approx$ 1/u. Choosing the “right” algorithm requires having an appreciation for these tradeoffs.

# 5.3.8 Iterative Improvement

A technique for refining an approximate LS solution has been analyzed by Bj¨orck (1967, 1968). It is based on the idea that if

$$
\left[ \begin{array}{l l} I _ {m} & A \\ A ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} r \\ x \end{array} \right] = \left[ \begin{array}{l} b \\ 0 \end{array} \right], \quad A \in \mathbb {R} ^ {m \times n}, b \in \mathbb {R} ^ {m}, \tag {5.3.20}
$$

then $\parallel b - A x \parallel _ { 2 } = \operatorname* { m i n }$ . This follows because $r + A x = b$ and $A ^ { T } r = 0$ imply $A ^ { T } A x =$ $A ^ { T } b$ . The above augmented system is nonsingular if $\mathsf { r a n k } ( A ) = n$ , which we hereafter assume. By casting the LS problem in the form of a square linear system, the iterative improvement scheme §3.5.3 can be applied:

$$
r ^ {(0)} = 0, x ^ {(0)} = 0
$$

$\mathbf { f o r } \ k = 0 , 1 , \ldots$

$$
\begin{array}{l} \left[ \begin{array}{c} f ^ {(k)} \\ g ^ {(k)} \end{array} \right] = \left[ \begin{array}{c} b \\ 0 \end{array} \right] - \left[ \begin{array}{c c} I & A \\ A ^ {T} & 0 \end{array} \right] \left[ \begin{array}{c} r ^ {(k)} \\ x ^ {(k)} \end{array} \right] \\ \left[ \begin{array}{c c} I & A \\ A ^ {T} & 0 \end{array} \right] \left[ \begin{array}{c} p ^ {(k)} \\ z ^ {(k)} \end{array} \right] = \left[ \begin{array}{c} f ^ {(k)} \\ g ^ {(k)} \end{array} \right] \\ \left[ \begin{array}{l} r ^ {(k + 1)} \\ x ^ {(k + 1)} \end{array} \right] = \left[ \begin{array}{l} r ^ {(k)} \\ x ^ {(k)} \end{array} \right] + \left[ \begin{array}{l} p ^ {(k)} \\ z ^ {(k)} \end{array} \right] \\ \end{array}
$$

end

The residuals $f ^ { ( k ) }$ and $g ^ { ( k ) }$ must be computed in higher precision, and an original copy of A must be around for this purpose.

If the QR factorization of A is available, then the solution of the augmented system is readily obtained. In particular, if $A = Q R$ and $R _ { 1 } = R ( 1 { : } n , 1 { : } n )$ , then a system of the form

$$
{\left[ \begin{array}{c c} I & A \\ A ^ {T} & 0 \end{array} \right]} {\left[ \begin{array}{c} p \\ z \end{array} \right]} = {\left[ \begin{array}{c} f \\ g \end{array} \right]}
$$

transforms to

$$
\left[ \begin{array}{c c c} I _ {n} & 0 & R _ {1} \\ 0 & I _ {m - n} & 0 \\ R _ {1} ^ {T} & 0 & 0 \end{array} \right] \left[ \begin{array}{c} h \\ f _ {2} \\ z \end{array} \right] = \left[ \begin{array}{c} f _ {1} \\ f _ {2} \\ g \end{array} \right]
$$

where

$$
Q ^ {T} f = \left[ \begin{array}{l} f _ {1} \\ f _ {2} \end{array} \right] _ {m - n} ^ {n}, \qquad Q ^ {T} p = \left[ \begin{array}{l} h \\ f _ {2} \end{array} \right] _ {m - n} ^ {n}.
$$

Thus, p and z can be determined by solving the triangular systems $R _ { 1 } ^ { T } h = g$ and $R _ { 1 } z = f _ { 1 } - h$ and setting

$$
p = Q \left[ \begin{array}{c} h \\ f _ {2} \end{array} \right].
$$

Assuming that Q is stored in factored form, each iteration requires 8mn $- 2 n ^ { 2 }$ flops.

The key to the iteration’s success is that both the LS residual and solution are updated—not just the solution. Bj¨orck (1968) shows that if $\kappa _ { 2 } ( A ) \approx \beta ^ { q }$ and t-digit, β-base arithmetic is used, then $x ^ { ( k ) }$ has approximately $k ( t - q )$ correct base-β digits, provided the residuals are computed in double precision. Notice that it is $\kappa _ { 2 } ( A )$ , not $\kappa _ { 2 } { \left( A \right) } ^ { 2 }$ , that appears in this heuristic.

# 5.3.9 Some Point/Line/Plane Nearness Problems in 3-Space

The fields of computer graphics and computer vision are replete with many interesting matrix problems. Below we pose three geometric “nearness” problems that involve points, lines, and planes in 3-space. Each is a highly structured least squares problem with a simple, closed-form solution. The underlying trigonometry leads rather naturally to the vector cross product, so we start with a quick review of this important operation.

The cross product of a vector $\boldsymbol { p } \in \mathbb { R } ^ { 3 }$ with a vector $q \in \mathbb { R } ^ { 3 }$ is defined by

$$
p \times q = \left[ \begin{array}{l} p _ {2} q _ {3} - p _ {3} q _ {2} \\ p _ {3} q _ {1} - p _ {1} q _ {3} \\ p _ {1} q _ {2} - p _ {2} q _ {1} \end{array} \right].
$$

This operation can be framed as a matrix-vector product. For any $v \in \mathbb { R } ^ { 3 }$ , define the skew-symmetric matrix $v ^ { c }$ by

$$
v ^ {c} = \left[ \begin{array}{c c c} 0 & - v _ {3} & v _ {2} \\ v _ {3} & 0 & - v _ {1} \\ - v _ {2} & v _ {1} & 0 \end{array} \right].
$$

It follows that

$$
p \times q = p ^ {c} \cdot q = - q ^ {c} \cdot p = - (q \times p).
$$

Using the skew-symmetry of $p ^ { c }$ and $q ^ { c }$ , it is easy to show that

$$
p \times q \in \operatorname{span} \{p, q \} ^ {\perp}. \tag {5.3.21}
$$

Other properties include

$$
(p \times q) \times r = (p ^ {c} \cdot q) ^ {c} r = (q p ^ {T} - p q ^ {T}) r = (p ^ {T} r) \cdot q - (q ^ {T} r) \cdot p, \tag {5.3.22}
$$

$$
(p \times q) ^ {T} (r \times s) = (p ^ {c} q) ^ {T} \cdot (r ^ {c} s) = \det ([ p q ] ^ {T} [ r s ]), \tag {5.3.23}
$$

$$
p ^ {c} p ^ {c} = p p ^ {T} - \parallel p \parallel_ {2} ^ {2} \cdot I _ {3}, \tag {5.3.24}
$$

$$
\parallel p ^ {c} q \parallel_ {2} ^ {2} = \parallel p \parallel_ {2} ^ {2} \cdot \parallel q \parallel_ {2} ^ {2} \cdot \left(1 - \left(\frac {p ^ {T} q}{\parallel p \parallel_ {2} \cdot \parallel q \parallel_ {2}}\right) ^ {2}\right). \tag {5.3.25}
$$

We are now set to state the three problems and specify their theoretical solutions. For hints at how to establish the correctness of the solutions, see P5.3.13–P5.3.15.

Problem 1. Given a line L and a point y, find the point $z ^ { \mathrm { o p t } }$ on $L$ that is closest to $y .$ i.e., solve

$$
\min _ {z \in L} \| z - y \| _ {2}.
$$

If L passes through distinct points $p _ { 1 }$ and $p _ { 2 }$ , then it can be shown that

$$
z ^ {\mathrm{opt}} = y + \frac {1}{v ^ {T} v} v ^ {c} v ^ {c} (y - p _ {1}), \quad v = p _ {2} - p _ {1}. \tag {5.3.26}
$$

Problem 2. Given lines $L _ { 1 }$ and $L _ { 2 }$ , find the point $z _ { 1 } ^ { \mathrm { o p t } }$ on $L _ { 1 }$ that is closest to $L _ { 2 }$ and the point $z _ { 2 } ^ { \mathrm { o p t } }$ on $L _ { 2 }$ that is closest to $L _ { \mathrm { 1 } } , \mathrm { i . e . }$ , solve

$$
\min _ {z _ {1} \in L _ {1}, z _ {2} \in L _ {2}} \| z _ {1} - z _ {2} \| _ {2}.
$$

If $L _ { 1 }$ passes through distinct points $p _ { 1 }$ and $p _ { 2 }$ and $L _ { 2 }$ passes through distinct points $q _ { 1 }$ and $q _ { 2 }$ , then it can be shown that

$$
z _ {1} ^ {\mathrm{opt}} = p _ {1} + \frac {1}{r ^ {T} r} \cdot v w ^ {T} \cdot r ^ {c} (q _ {1} - p _ {1}), \tag {5.3.27}
$$

$$
z _ {2} ^ {\mathrm{opt}} = q _ {1} + \frac {1}{r ^ {T} r} \cdot w v ^ {T} \cdot r ^ {c} (q _ {1} - p _ {1}), \tag {5.3.28}
$$

where v = p2 − p1, w = q2 − q1, and $r = v ^ { c } w$ .

Problem 3. Given a plane $P$ and a point y, find the point $z ^ { \mathrm { o p t } }$ on $P$ that is closest to $y , { \mathrm { i . e . } }$ , solve

$$
\min _ {z \in P} \| z - y \| _ {2}.
$$

If P passes through three distinct points $p _ { 1 } , p _ { 2 }$ , and $p _ { 3 }$ , then it can be shown that

$$
z ^ {\mathrm{opt}} = p _ {1} - \frac {1}{v ^ {T} v} \cdot v ^ {c} v ^ {c} (y - p _ {1}) \tag {5.3.29}
$$

where $v = ( p _ { 2 } - p _ { 1 } ) ^ { c } ( p _ { 3 } - p _ { 1 } )$ .

The nice closed-form solutions (5.3.26)–(5.3.29) are deceptively simple and great care must be exercised when computing with these formulae or their mathematical equivalents. See Kahan (2011).

# Problems

P5.3.1 Assume $A ^ { T } A x = A ^ { T } b , ( A ^ { T } A + F ) \hat { x } = A ^ { T } b ,$ , and $2 \| F \| _ { 2 } \leq \sigma _ { n } ( A ) ^ { 2 }$ . Show that if $r = b - A x$ and $\hat { r } = b - A \hat { x }$ , then $\hat { r } - r = A ( A ^ { T } A + F ) ^ { - 1 } F x$ and

$$
\| \hat {r} - r \| _ {2} \leq 2 \kappa_ {2} (A) \frac {\| F \| _ {2}}{\| A \| _ {2}} \| x \| _ {2}.
$$

P5.3.2 Assume that $A ^ { T } A x = A ^ { T } b$ and that $A ^ { T } A \hat { x } = A ^ { T } b + f$ where $\parallel f \parallel _ { 2 } \leq c \mathbf { u } \parallel A ^ { T } \parallel _ { 2 } \parallel b \parallel _ { 2 }$ and A has full column rank. Show that

$$
\frac {\| x - \hat {x} \| _ {2}}{\| x \| _ {2}} \leq c \mathbf {u} \kappa_ {2} (A) ^ {2} \frac {\| A ^ {T} \| _ {2} \| b \| _ {2}}{\| A ^ {T} b \|}.
$$

P5.3.3 Let $A \in \mathbb { R } ^ { m \times n } ( m \geq n ) , w \in \mathbb { R } ^ { n }$ , and define

$$
B = \left[ \begin{array}{c} A \\ w ^ {T} \end{array} \right].
$$

Show that $\sigma _ { n } ( B ) \geq \sigma _ { n } ( A )$ and $\sigma _ { 1 } ( B ) \leq \sqrt { \| A \| _ { 2 } ^ { 2 } + \| w \| _ { 2 } ^ { 2 } }$ . Thus, the condition of a matrix may increase or decrease if a row is added.

P5.3.4 (Cline 1973) Suppose that $A \in \mathbb { R } ^ { m \times n }$ has rank n and that Gaussian elimination with partial pivoting is used to compute the factorization $P A = L U$ , where $L \in \mathbb { R } ^ { m \times n }$ is unit lower triangular, $\bar { U } \in \mathbb { R } ^ { n \times n }$ is upper triangular, and $P \in \mathbb { R } ^ { m \times m }$ is a permutation. Explain how the decomposition in P5.2.5 can be used to find a vector $z \in \mathbb { R } ^ { n }$ such that  $L z - P b \| _ { 2 }$ is minimized. Show that if $U x = z ,$ then $\parallel A x - b \parallel _ { 2 }$ is minimum. Show that this method of solving the LS problem is more efficient than Householder QR from the flop point of view whenever m $\leq 5 n / 3$ .

P5.3.5 The matrix $C = ( A ^ { T } A ) ^ { - 1 }$ , where rank $( A ) = n$ , arises in many statistical applications. Assume that the factorization $A = Q R$ is available. (a) Show $C = ( R ^ { T } \dot { R } ) ^ { - 1 }$ . (b) Give an algorithm for computing the diagonal of C that requires $n ^ { 3 } / 3$ flops. (c) Show that

$$
R = \left[ \begin{array}{c c} \alpha & v ^ {T} \\ 0 & S \end{array} \right] \qquad \Rightarrow \qquad C = (R ^ {T} R) ^ {- 1} = \left[ \begin{array}{c c} (1 + v ^ {T} C _ {1} v) / \alpha^ {2} & - v ^ {T} C _ {1} / \alpha \\ - C _ {1} v / \alpha & C _ {1} \end{array} \right]
$$

where $C _ { 1 } = ( S ^ { T } S ) ^ { - 1 }$ . (d) Using (c), give an algorithm that overwrites the upper triangular portion of R with the upper triangular portion of C. Your algorithm should require $2 { n } ^ { 3 } / 3$ flops.

P5.3.6 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and that $r = b - A \acute { \iota }$ x where $r , b , x \in \mathbb { R } ^ { n }$ and x is nonzero. Show how to compute a symmetric $E \in \mathbb { R } ^ { n \times n }$ with minimal Frobenius norm so that $( A + E ) x = b$ . Hint: Use the QR factorization of $\left[ \boldsymbol { x } \mid \boldsymbol { r } \right]$ and note that $E x = r \Rightarrow ( Q ^ { T } E Q ) ( Q ^ { T } x ) = \dot { Q } ^ { T } r$ .

P5.3.7 Points $P _ { 1 } , \ldots , P _ { n }$ on the x-axis have x-coordinates $x _ { 1 } , \ldots , x _ { n }$ . We know that $x _ { 1 } = 0$ and wish to compute $x _ { 2 } , \ldots , x _ { n }$ given that we have estimates $d _ { i j }$ of the separations:

$$
x _ {i} - x _ {j} \approx d _ {i j}, \qquad 1 \leq i <   j \leq n.
$$

Using the method of normal equations, show how to minimize

$$
\phi (x _ {1}, \ldots , x _ {n}) = \sum_ {i = 1} ^ {n - 1} \sum_ {j = i + 1} ^ {n} (x _ {i} - x _ {j} - d _ {i j}) ^ {2}
$$

subject to the constraint $x _ { 1 } = 0$ .

P5.3.8 Suppose $A \in \mathbb { R } ^ { m \times n }$ has full rank and that $b \in \mathbb { R } ^ { m }$ and $c \in \mathbb { R } ^ { n }$ are given. Show how to compute $\alpha = c ^ { T } x _ { \mathrm { L S } }$ without computing $x _ { \mathrm { L S } }$ explicitly. Hint: Suppose $Z$ is a Householder matrix such that $Z ^ { T } c$ is a multiple of $I _ { n } ( : , n )$ . It follows that $\alpha = ( Z ^ { T } c ) ^ { \bar { T } } y _ { \mathrm { L S } }$ where $y _ { \mathrm { L S } }$ minimizes $\parallel \tilde { A } y - b \parallel _ { 2 }$ with $y = Z ^ { T } x$ and $\tilde { A } = A Z$ .

P5.3.9 Suppose $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ with $m \geq n$ . How would you solve the full rank least squares problem given the availability of a matrix $\boldsymbol { M } \in \mathbb { R } ^ { m \times m }$ such that $M ^ { T } A = S$ is upper triangular and $\mathbf { \overset { \cdot } { M } } ^ { T } \mathbf { \overset { \cdot } { M } } = \mathbf { \overset { \cdot } { D } }$ is diagonal?

P5.3.10 Let $A \in \mathbb { R } ^ { m \times n }$ have rank n and for $\alpha \geq 0$ define

$$
M (\alpha) = \left[ \begin{array}{c c} \alpha I _ {m} & A \\ A ^ {T} & 0 \end{array} \right].
$$

Show that

$$
\sigma_ {m + n} (M (\alpha)) = \min \left\{\alpha , - \frac {\alpha}{2} + \sqrt {\sigma_ {n} (A) ^ {2} + \left(\frac {\alpha}{2}\right) ^ {2}} \right\}
$$

and determine the value of α that minimizes $\kappa _ { 2 } ( M ( \alpha ) )$ .

P5.3.11 Another iterative improvement method for LS problems is the following:

$$
x ^ {(0)} = 0
$$

$\mathbf { f o r } \ k = 0 , 1 , \ldots$

$$
r ^ {(k)} = b - A x ^ {(k)} \quad (\text { double   precision })
$$

$$
\| A z ^ {(k)} - r ^ {(k)} \| _ {2} = \min
$$

$$
x ^ {(k + 1)} = x ^ {(k)} + z ^ {(k)}
$$

end

(a) Assuming that the QR factorization of A is available, how many flops per iteration are required? (b) Show that the above iteration results by setting $g ^ { ( k ) } = 0$ in the iterative improvement scheme given in 5.3.8.

P5.3.12 Verify (5.3.21)–(5.3.25).

P5.3.13 Verify (5.3.26) noting that $L = \left\{ p _ { 1 } + \tau ( p _ { 2 } - p _ { 1 } ) : \tau \in \mathbb { R } \right\}$

P5.3.14 Verify (5.3.27) noting that the minimizer $\tau ^ { \mathrm { o p t } } \in \mathbb { R } ^ { 2 } \ \mathrm { o f } \ \| \ ( p _ { 1 } - q _ { 1 } ) \ - \ \big [ \ : p _ { 2 } - p _ { 1 } \ : | \ : q _ { 2 } - q _ { 1 } \ : \big ] \tau \| _ { 2 }$ is relevant.

P5.3.15 Verify (5.3.29) noting that $P = \{ x : x ^ { T } ( ( p _ { 2 } - p _ { 1 } ) \times ( p _ { 3 } - p _ { 1 } ) ) = 0 .$

# Notes and References for 5.3

Some classical references for the least squares problem include:

F.L. Bauer (1965). “Elimination with Weighted Row Combinations for Solving Linear Equations and Least Squares Problems,” Numer. Math. 7, 338–352.

G.H. Golub and J.H. Wilkinson (1966). “Note on the Iterative Refinement of Least Squares Solution,” Numer. Math. 9, 139–148.

A. van der Sluis (1975). “Stability of the Solutions of Linear Least Squares Problem,” Numer. Math. 23, 241–254.

The use of Gauss transformations to solve the LS problem has attracted some attention because they are cheaper to use than Householder or Givens matrices, see:

G. Peters and J.H. Wilkinson (1970). “The Least Squares Problem and Pseudo-Inverses,” Comput. J. 13, 309–16.

A.K. Cline (1973). “An Elimination Method for the Solution of Linear Least Squares Problems,” SIAM J. Numer. Anal. 10, 283–289.

R.J. Plemmons (1974). “Linear Least Squares by Elimination and MGS,” J. ACM 21, 581–585.

The seminormal equations are given by $R ^ { T } R x = A ^ { T } b$ where $A = Q R .$ . It can be shown that by solving the seminormal equations an acceptable $\mathrm { L S }$ solution is obtained if one step of fixed precision iterative improvement is performed, see:

˚A. Bj¨orck (1987). “Stability Analysis of the Method of Seminormal Equations,” Lin. Alg. Applic. 88/89, 31–48.

Survey treatments of LS perturbation theory include Lawson and Hanson (SLS), Stewart and Sun (MPT), and Bj¨orck (NMLS). See also:

P.-A. Wedin (1973). “Perturbation Theory for Pseudoinverses,” BIT 13, 217–232.

˚A. Bj¨orck (1991). “Component-wise Perturbation Analysis and Error Bounds for Linear Least Squares Solutions,” BIT 31, 238–244.

B. Wald´en, R. Karlson, J. Sun (1995). “Optimal Backward Perturbation Bounds for the Linear Least Squares Problem,” Numerical Lin. Alg. Applic. 2, 271–286.

J.-G. Sun (1996). “Optimal Backward Perturbation Bounds for the Linear Least-Squares Problem with Multiple Right-Hand Sides,” IMA J. Numer. Anal. 16, 1–11.

J.-G. Sun (1997). “On Optimal Backward Perturbation Bounds for the Linear Least Squares Problem,” BIT 37, 179–188.

R. Karlson and B. Wald´en (1997). “Estimation of Optimal Backward Perturbation Bounds for the Linear Least Squares Problem,” BIT 37, 862–869.

J.-G. Sun (1997). “On Optimal Backward Perturbation Bounds for the Linear Least Squares Problem,” BIT 37, 179–188.

M. Gu (1998). “Backward Perturbation Bounds for Linear Least Squares Problems,” SIAM J. Matrix Anal. Applic. 20, 363–372.

M. Arioli, M. Baboulin and S. Gratton (2007). “A Partial Condition Number for Linear Least Squares Problems,” SIAM J. Matrix Anal. Applic. 29, 413–433.

M. Baboulin, J. Dongarra, S. Gratton, and J. Langou (2009). “Computing the Conditioning of the Components of a Linear Least-Squares Solution,” Num. Lin. Alg. Applic. 16, 517–533.

M. Baboulin and S. Gratton (2009). “Using Dual Techniques to Derive Componentwise and Mixed Condition Numbers for a Linear Function of a Least Squares Solution,” BIT 49, 3–19.

J. Grcar (2009). “Nuclear Norms of Rank-2 Matrices for Spectral Condition Numbers of Rank Least Squares Solutions,” ArXiv:1003.2733v4.

J. Grcar (2010). “Spectral Condition Numbers of Orthogonal Projections and Full Rank Linear Least Squares Residuals,” SIAM J. Matrix Anal. Applic. 31, 2934–2949.

Practical insights into the accuracy of a computed least squares solution can be obtained by applying the condition estimation ideas of 3.5. to the R matrix in $A = Q R$ or the Cholesky factor of $\bar { A } ^ { T } \bar { A }$ should a normal equation approach be used. For a discussion of LS-specific condition estimation, see:

G.W. Stewart (1980). “The Efficient Generation of Random Orthogonal Matrices with an Application to Condition Estimators,” SIAM J. Numer. Anal. 17, 403–9.

S. Gratton (1996). “On the Condition Number of Linear Least Squares Problems in a Weighted Frobenius Norm,” BIT 36, 523–530.

C.S. Kenney, A.J. Laub, and M.S. Reese (1998). “Statistical Condition Estimation for Linear Least Squares,” SIAM J. Matrix Anal. Applic. 19, 906–923.

Our restriction to least squares approximation is not a vote against minimization in other norms. There are occasions when it is advisable to minimize $\| A x - b \| _ { p }$ for $p = 1$ and ∞. Some algorithms for doing this are described in:

A.K. Cline (1976). “A Descent Method for the Uniform Solution to Overdetermined Systems of Equations,” SIAM J. Numer. Anal. 13, 293–309.

R.H. Bartels, A.R. Conn, and C. Charalambous (1978). “On Cline’s Direct Method for Solving Overdetermined Linear Systems in the $L _ { \infty }$ Sense,” SIAM J. Numer. Anal. 15, 255–270.

T.F. Coleman and Y. Li (1992). “A Globally and Quadratically Convergent Affine Scaling Method for Linear $L _ { 1 }$ Problems,” Mathematical Programming 56, Series A, 189–222.

Y. Li (1993). $^ { 6 6 } \mathrm { A }$ Globally Convergent Method for $L _ { p }$ Problems,” SIAM J. Optim. 3, 609–629.

Y. Zhang (1993). “A Primal-Dual Interior Point Approach for Computing the $L _ { 1 }$ and $L _ { \infty }$ Solutions of Overdetermined Linear Systems,” J. Optim. Theory Applic. $^ { 7 7 , }$ 323–341.

Iterative improvement in the least squares context is discussed in:

G.H. Golub and J.H. Wilkinson (1966). “Note on Iterative Refinement of Least Squares Solutions,” Numer. Math. 9, 139–148.

˚A. Bj¨orck and G.H. Golub (1967). “Iterative Refinement of Linear Least Squares Solutions by Householder Transformation,” BIT 7, 322–337.

˚A. Bj¨orck (1967). “Iterative Refinement of Linear Least Squares Solutions I,” BIT 7, 257–278.

˚A. Bj¨orck (1968). “Iterative Refinement of Linear Least Squares Solutions II,”BIT 8, 8–30.

J. Gluchowska and A. Smoktunowicz (1999). “Solving the Linear Least Squares Problem with Very High Relative Acuracy,” Computing 45, 345–354.

J. Demmel, Y. Hida, and E.J. Riedy (2009). “Extra-Precise Iterative Refinement for Overdetermined Least Squares Problems,” ACM Trans. Math. Softw. 35, Article 28.

The following texts treat various geometric matrix problems that arise in computer graphics and vision:

A.S. Glassner (1989). An Introduction to Ray Tracing, Morgan Kaufmann, Burlington, MA.

R. Hartley and A. Zisserman (2004). Multiple View Geometry in Computer Vision, Second Edition, Cambridge University Press, New York.

M. Pharr and M. Humphreys (2010). Physically Based Rendering, from Theory to Implementation, Second Edition, Morgan Kaufmann, Burlington, MA.

For a numerical perspective, see:

W. Kahan (2008). “Computing Cross-Products and Rotations in 2- and 3-dimensional Euclidean Spaces,” http://www.cs.berkeley.edu/ wkahan/MathH110/Cross.pdf.

# 5.4 Other Orthogonal Factorizations

Suppose $A \in \mathbb { R } ^ { m \times 4 }$ has a thin QR factorization of the following form:

$$
A = \left[ a _ {1}, a _ {2}, a _ {3}, a _ {4} \right] = \left[ q _ {1}, q _ {2}, q _ {3}, q _ {4} \right] \left[ \begin{array}{c c c c} 1 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 1 \end{array} \right].
$$

Note that ran(A) has dimension 3 but does not equal span $\{ q _ { 1 } , q _ { 2 } , q _ { 3 } \}$ , span $\{ q _ { 1 } , q _ { 2 } , q _ { 4 } \}$ , span $\{ q _ { 1 } , q _ { 3 } , q _ { 4 } \}$ , or span $\{ q _ { 2 } , q _ { 3 } , q _ { 4 } \}$ because $a _ { 4 }$ does not belong to any of these subspaces. In this case, the QR factorization reveals neither the range nor the nullspace of A and the number of nonzeros on R’s diagonal does not equal its rank. Moreover, the LS solution process based on the QR factorization (Algorithm 5.3.2) breaks down because the upper triangular portion of R is singular.

We start this section by introducing several decompositions that overcome these shortcomings. They all have the form $\mathbf { \bar { \phi } } _ { Q ^ { T } A Z } = T$ where T is a structured block triangular matrix that sheds light on A’s rank, range, and nullspace. We informally refer to matrix reductions of this form as rank revealing. See Chandrasekaren and Ipsen (1994) for a more precise formulation of the concept.

Our focus is on a modification of the QR factorization that involves column pivoting. The resulting R-matrix has a structure that supports rank estimation. To set the stage for updating methods, we briefly discus the $U L V$ and UT V frameworks Updating is discussed in §6.5 and refers to the efficient recomputation of a factorization after the matrix undergoes a low-rank change.

All these methods can be regarded as inexpensive alternatives to the SVD, which represents the “gold standard” in the area of rank determination. Nothing “takes apart” a matrix so conclusively as the SVD and so we include an explanation of its airtight reliability. The computation of the full SVD, which we discuss in §8.6, begins with the reduction to bidiagonal form using Householder matrices. Because this decomposition is important in its own right, we provide some details at the end of this section.

# 5.4.1 Numerical Rank and the SVD

Suppose $A \in \mathbb { R } ^ { m \times n }$ has SVD $U ^ { T } A V = \Sigma = \mathrm { d i a g } ( \sigma _ { i } )$ . If $\mathsf { r a n k } ( A ) = r < n$ , then according to the exact arithmetic discussion of §2.4 the singular values $\sigma _ { r + 1 } , \ldots , \sigma _ { n }$ are zero and

$$
A = \sum_ {i = 1} ^ {r} \sigma_ {k} u _ {k} v _ {k} ^ {T}. \tag {5.4.1}
$$

The exposure of rank degeneracy could not be more clear.

In Chapter 8 we describe the Golub-Kahan-Reinsch algorithm for computing the SVD. Properly implemented, it produces nearly orthogonal matrices $\widehat { U }$ and $\widehat { V }$ so that

$$
\widehat {U} ^ {T} A \widehat {V} \approx \widehat {\Sigma} = \mathrm{diag} (\widehat {\sigma} _ {1}, \ldots , \widehat {\sigma} _ {n}), \qquad \widehat {\sigma} _ {1} \geq \dots \geq \widehat {\sigma} _ {n} \geq 0.
$$

(Other SVD procedures have this property as well.) Unfortunately, unless remarkable cancellation occurs, none of the computed singular values will be zero because of roundoff error. This forces an issue. On the one hand, we can adhere to the strict mathematical definition of rank, count the number of nonzero computed singular values, and conclude from

$$
A \approx \sum_ {i = 1} ^ {n} \widehat {\sigma} _ {k} \widehat {u} _ {k} \widehat {v} _ {k} ^ {T} \tag {5.4.2}
$$

that A has full rank. However, working with every matrix as if it possessed full column rank is not particularly useful. It is more productive to liberalize the notion of rank by setting small computed singular values to zero in (5.4.2). This results in an approximation of the form

$$
A \approx \sum_ {i = 1} ^ {\widehat {r}} \widehat {\sigma} _ {k} \widehat {u} _ {k} \widehat {v} _ {k} ^ {T}, \quad \widehat {r} \leq \widehat {n} \tag {5.4.3}
$$

where we regard $\widehat { r }$ as the numerical rank. For this approach to make sense we need to guarantee that $| \widehat { \sigma } _ { i } - \sigma _ { i } |$ is small.

<For a properly implemented Golub-Kahan-Reinsch SVD algorithm, it can be shown that

$$
\widehat {U} = W + \Delta U, W ^ {T} W = I _ {m}, \quad \| \Delta U \| _ {2} \leq \epsilon ,
$$

$$
\widehat {V} = Z + \Delta V, Z ^ {T} Z = I _ {n}, \quad \| \Delta V \| _ {2} \leq \epsilon , \tag {5.4.4}
$$

$$
\widehat {\Sigma} = W ^ {T} (A + \Delta A) Z, \quad \| \Delta A \| _ {2} \leq \epsilon \| A \| _ {2},
$$

where $\epsilon$ is a small multiple of u, the machine precision. In other words, the SVD algorithm computes the singular values of a nearby matrix $A + \Delta A$ .


---

<!-- golub_300_349 -->

Note that $\widehat { U }$ and $\widehat { V }$ are not necessarily close to their exact counterparts. However, we can show that $\widehat { \sigma } _ { k }$ is close to $\sigma _ { k }$ as follows. Using Corollary 2.4.6 we have

$$
\sigma_ {k} = \min _ {\operatorname{rank} (B) = k - 1} \| A - B \| _ {2} = \min _ {\operatorname{rank} (B) = k - 1} \| (\widehat {\Sigma} - B) - E \| _ {2}
$$

where

$$
E = W ^ {T} (\Delta A) Z
$$

and

$$
\| E \| _ {2} \leq \epsilon \| A \| _ {2} = \epsilon \sigma_ {1}.
$$

Since

$$
\| \widehat {\Sigma} - B \| - \| E \| \leq \| \widehat {\Sigma} - B \| \leq \| \widehat {\Sigma} - B \| + \| E \|
$$

and

$$
\min _ {\operatorname{rank} (B) = k - 1} \| \widehat {\Sigma} _ {k} - B \| _ {2} = \widehat {\sigma} _ {k},
$$

it follows that

$$
\left| \sigma_ {k} - \widehat {\sigma} _ {k} \right| \leq \epsilon \sigma_ {1}
$$

for $k = 1 { : } n$ . Thus, if A has rank r, then we can expect $n - r$ of the computed singular values to be small. Near rank deficiency in A cannot escape detection if the SVD of A is computed.

Of course, all this hinges on having a definition of “small.” This amounts to choosing a tolerance $\delta > 0$ and declaring A to have numerical rank r if the computed singular values satisfy

$$
\widehat {\sigma} _ {1} \geq \dots \geq \widehat {\sigma} _ {\hat {r}} > \delta \geq \widehat {\sigma} _ {\hat {r} + 1} \geq \dots \geq \widehat {\sigma} _ {n}. \tag {5.4.5}
$$

We refer to the integer $\hat { r }$ as the δ-rank of A. The tolerance should be consistent with the machine precision, $\mathrm { e . g . , } \delta = \mathbf { u } \parallel A \parallel _ { \infty }$ . However, if the general level of relative error in the data is larger than u, then δ should be correspondingly bigger, e.g., $\delta = 1 0 ^ { - 2 } \parallel A \parallel _ { \infty }$ if the entries in A are correct to two digits.

For a given δ it is important to stress that, although the SVD provides a great deal of rank-related insight, it does not change the fact that the determination of numerical rank is a sensitive computation. If the gap between $\widehat { \sigma } _ { \widehat { r } }$ and $\widehat \sigma _ { r + 1 }$ is small, then A is also close (in the δ sense) to a matrix with rank $\widehat { r } - 1$ < <. Thus, the amount of confidence we have in the correctness of $\widehat { r }$ <and in how we proceed to use the approximation (5.4.2) depends on the gap between $\widehat { \sigma } _ { \widehat { r } }$ and $\widehat \sigma _ { r + 1 }$ .

# 5.4.2 QR with Column Pivoting

We now examine alternative rank-revealing strategies to the SVD starting with a modification of the Householder QR factorization procedure (Algorithm 5.2.1). In exact arithmetic, the modified algorithm computes the factorization

$$
Q ^ {T} A \Pi = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & 0 \\ r & n - r \end{array} \right] _ {m - r} ^ {r} \tag {5.4.6}
$$

where $r = \mathsf { r a n k } ( A )$ , Q is orthogonal, $R _ { 1 1 }$ is upper triangular and nonsingular, and Π is a permutation. If we have the column partitionings $A \Pi = \left[ \left. a _ { c _ { 1 } } \right| \cdot \cdot \cdot \right| \left. a _ { c _ { n } } \right]$ and $Q = [  q _ { 1 } | \cdots |  q _ { m }  ]$ , then for $k = 1 { : } n$ we have

$$
a _ {c _ {k}} = \sum_ {i = 1} ^ {\min \{r, k \}} r _ {i k} q _ {i} \in \operatorname{span} \left\{q _ {1}, \dots , q _ {r} \right\}
$$

implying

$$
\operatorname{ran} (A) = \operatorname{span} \left\{q _ {1}, \dots , q _ {r} \right\}.
$$

To see how to compute such a factorization, assume for some k that we have computed Householder matrices $H _ { 1 } , \ldots , H _ { k - 1 }$ and permutations $\Pi _ { 1 } , \ldots , \Pi _ { k - 1 }$ such that

$$
(H _ {k - 1} \dots H _ {1}) A (\Pi_ {1} \dots \Pi_ {k - 1}) = R ^ {(k - 1)} = \left[ \begin{array}{c c} R _ {1 1} ^ {(k - 1)} & R _ {1 2} ^ {(k - 1)} \\ 0 & R _ {2 2} ^ {(k - 1)} \end{array} \right] _ {m - k + 1} ^ {k - 1} \tag {5.4.7}
$$

$R _ { 1 1 } ^ { ( k - 1 ) }$ is a nonsingular and upper triangular matrix. Now suppose that

$$
R _ {2 2} ^ {(k - 1)} = \left[ z _ {k} ^ {(k - 1)} \mid \dots \mid z _ {n} ^ {(k - 1)} \right]
$$

is a column partitioning and let $p \geq k$ be the smallest index such that

$$
\left\| z _ {p} ^ {(k - 1)} \right\| _ {2} = \max \left\{\left\| z _ {k} ^ {(k - 1)} \right\| _ {2}, \dots , \left\| z _ {n} ^ {(k - 1)} \right\| _ {2} \right\}. \tag {5.4.8}
$$

Note that if $\mathsf { r a n k } ( A ) = k - 1$ , then this maximum is zero and we are finished. Otherwise, let $\Pi _ { k }$ be the $n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n$ identity with columns p and k interchanged and determine a Householder matrix $H _ { k }$ such that if

$$
R ^ {(k)} = H _ {k} R ^ {(k - 1)} \Pi_ {k},
$$

then $R ^ { ( k ) } ( k + 1 { : } m , k ) = 0$ . In other words, $\Pi _ { k }$ moves the largest column in $R _ { 2 2 } ^ { ( k - 1 ) }$ the lead position and $H _ { k }$ zeroes all of its subdiagonal components.

The column norms do not have to be recomputed at each stage if we exploit the property

$$
Q ^ {T} z = \left[ \begin{array}{c} \alpha \\ w \end{array} \right] _ {s - 1} ^ {1} \qquad \Longrightarrow \qquad \| w \| _ {2} ^ {2} = \| z \| _ {2} ^ {2} - \alpha^ {2},
$$

which holds for any orthogonal matrix $Q \in \mathbb { R } ^ { s \times s }$ . This reduces the overhead associated with column pivoting from $O ( m n ^ { 2 } )$ flops to $O ( m n )$ flops because we can get the new column norms by updating the old column norms, e.g.,

$$
\parallel z _ {j} ^ {(k)} \parallel_ {2} ^ {2} = \parallel z _ {j} ^ {(k - 1)} \parallel_ {2} ^ {2} - r _ {k j} ^ {2} \qquad j = k + 1: n.
$$

Combining all of the above we obtain the following algorithm first presented by Businger and Golub (1965):

Algorithm 5.4.1 (Householder QR With Column Pivoting) Given $A \in \mathbb { R } ^ { m \times n }$ with m $\geq n$ , the following algorithm computes $r = \mathsf { r a n k } ( A )$ and the factorization (5.4.6) with $Q = H _ { 1 } \cdot \cdot \cdot H _ { r }$ and $\Pi = \Pi _ { 1 } \cdots \Pi _ { r }$ . The upper triangular part of A is overwritten by the upper triangular part of R and components $j + 1 { : } m$ of the jth Householder vector are stored in $A ( j + 1 { : } m , j )$ . The permutation Π is encoded in an integer vector piv. In particular, $\Pi _ { j }$ is the identity with rows $j$ and $p i v ( j )$ interchanged.

for $j = 1:n$ $c(j) = A(1:m,j)^T A(1:m,j)$ end $r = 0$ $\tau = \max\{c(1),\ldots,c(n)\}$ while $\tau > 0$ and $r < n$ $r = r + 1$ Find smallest $k$ with $r \leq k \leq n$ so $c(k) = \tau$ . $piv(r) = k$ $A(1:m,r) \leftrightarrow A(1:m,k)$ $c(r) \leftrightarrow c(k)$ $[v,\beta] = \text{house}(A(r:m,r))$ $A(r:m,r:n) = (I_{m-r+1} - \beta vv^T)A(:r:m,r:n)$ $A(r+1:m,r) = v(2:m-r+1)$ for $i = r+1:n$ $c(i) = c(i) - A(r,i)^2$ end $\tau = \max\{c(r+1),\ldots,c(n)\}$ end

This algorithm requires 4mnr $- 2 r ^ { 2 } ( m + n ) + 4 r ^ { 3 } / 3$ flops where r = rank(A).

# 5.4.3 Numerical Rank and A  QR

In principle, QR with column pivoting reveals rank. But how informative is the method in the context of floating point arithmetic? After k steps we have

$$
\mathsf {f l} (H _ {k} \dots H _ {1} A \Pi_ {1} \dots \Pi_ {k}) = \widehat {R} ^ {(k)} = \left[ \begin{array}{c c} \widehat {R} _ {1 1} ^ {(k)} & \widehat {R} _ {1 2} ^ {(k)} \\ 0 & \widehat {R} _ {2 2} ^ {(k)} \end{array} \right] _ {m - k} ^ {k}. \tag {5.4.9}
$$

If ${ \widehat { R } } _ { 2 2 } ^ { ( k ) }$ is suitably small in norm, then it is reasonable to terminate the reduction and declare A to have rank k. A typical termination criteria might be

$$
\| \widehat {R} _ {2 2} ^ {(k)} \| _ {2} \leq \epsilon_ {1} \| A \| _ {2}
$$

for some small machine-dependent parameter $\epsilon _ { 1 }$ . In view of the roundoff properties associated with Householder matrix computation (cf. §5.1.12), we know that $\widehat { R } ^ { ( k ) }$ is the exact R-factor of a matrix $A + E _ { k }$ , where

$$
\| E _ {k} \| _ {2} \leq \epsilon_ {2} \| A \| _ {2}, \quad \epsilon_ {2} = O (\mathbf {u}).
$$

Using Corollary 2.4.4 we have

$$
\sigma_ {k + 1} (A + E _ {k}) = \sigma_ {k + 1} (\widehat {R} ^ {(k)}) \leq \| \widehat {R} _ {2 2} ^ {(k)} \| _ {2}.
$$

Since $\sigma _ { k + 1 } ( A ) \leq \sigma _ { k + 1 } ( A + E _ { k } ) + \| E _ { k } \| _ { 2 }$ , it follows that

$$
\sigma_ {k + 1} (A) \leq (\epsilon_ {1} + \epsilon_ {2}) \| A \| _ {2}.
$$

In other words, a relative perturbation of $O ( \epsilon _ { 1 } + \epsilon _ { 2 } )$ in A can yield a rank-k matrix. With this termination criterion, we conclude that QR with column pivoting discovers rank deficiency if ${ \widehat { R } } _ { 2 2 } ^ { ( k ) }$ is small for some $k < n$ . However, it does not follow that the matrix ${ \widehat { R } } _ { 2 2 } ^ { ( k ) }$ <22 in (5.4.9) is small if rank $( A ) = k$ . There are examples of nearly rank deficient matrices whose R-factor look perfectly “normal.” A famous example is the Kahan matrix

$$
\mathsf {K a h} _ {n} (s) = \operatorname{diag} (1, s, \ldots , s ^ {n - 1}) \left[ \begin{array}{c c c c c} 1 & - c & - c & \dots & - c \\ 0 & 1 & - c & \dots & - c \\ & \ddots & & \vdots & \vdots \\ \vdots & & & 1 & - c \\ 0 & & \dots & & 1 \end{array} \right].
$$

Here, $c ^ { 2 } + s ^ { 2 } = 1$ with $c , s > 0$ . (See Lawson and Hanson (SLS, p. 31).) These matrices are unaltered by Algorithm 5.4.1 and thus $\parallel R _ { 2 2 } ^ { ( k ) } \parallel _ { 2 } \geq s ^ { n - 1 }$ for $k = 1 { : } n - 1$ . This inequality implies (for example) that the matrix $\mathsf { K a h } _ { 3 0 0 } ( . 9 9 )$ has no particularly small trailing principal submatrix since $s ^ { 2 9 9 } \approx . 0 5$ . However, a calculation shows that σ300 $= O ( 1 0 ^ { - 1 9 } )$ .

Nevertheless, in practice, small trailing R-submatrices almost always emerge that correlate well with the underlying rank. In other words, it is almost always the case that ${ \widehat { R } } _ { 2 2 } ^ { ( k ) }$ is small if A has rank k.

# 5.4.4 Finding a Good Column Ordering

It is important to appreciate that Algorithm 5.4.1 is just one way to determine the column permutation Π. The following result sets the stage for a better way.

Theorem 5.4.1. If $A \in \mathbb { R } ^ { m \times n }$ and $v \in \mathbb { R } ^ { n }$ is a unit 2-norm vector, then there exists a permutation Π so that the QR factorization

$$
A \Pi = Q R
$$

satisfies $| r _ { n n } | \leq { \sqrt { n } } \sigma$ where $\sigma = \parallel A v \parallel _ { 2 }$ .

Proof. Suppose $\Pi \in \mathbb { R } ^ { n \times n }$ is a permutation such that if $w = \Pi ^ { T } v$ , then

$$
\left| w _ {n} \right| = \max \left| v _ {i} \right|.
$$

Since $w _ { n }$ is the largest component of a unit 2-norm vector, $| w _ { n } | \geq 1 / { \sqrt { n } }$ . If $A \Pi = Q R$ is a QR factorization, then

$$
\sigma = \parallel A v \parallel_ {2} = \parallel (Q ^ {T} A \Pi) (\Pi^ {T} v) \parallel_ {2} = \parallel R (1: n, 1: n) w \parallel_ {2} \geq | r _ {n n} w _ {n} | \geq | r _ {n n} | / \sqrt {n}. \quad \square
$$

Note that if $v = v _ { n }$ is the right singular vector corresponding to $\sigma _ { \mathrm { m i n } } ( A )$ , then $| r _ { n n } | \leq$ $\sqrt { n } \sigma _ { n }$ . This suggests a framework whereby the column permutation matrix Π is based on an estimate of $v _ { n }$ :

Step 1. Compute the QR factorization $A = Q _ { 0 } R _ { 0 }$ and note that $R _ { 0 }$ has the same right singular vectors as A.

Step 2. Use condition estimation techniques to obtain a unit vector v with $\parallel R _ { 0 } v \parallel _ { 2 } \approx \sigma _ { n }$ .

Step 3. Determine Π and the QR factorization $A \Pi = Q R .$ .

See Chan (1987) for details about this approach to rank determination. The permutation Π can be generated as a sequence of swap permutations. This supports a very economical Givens rotation method for generating of Q and R from $Q _ { 0 }$ and $R _ { 0 }$ .

# 5.4.5 More General Rank-Revealing Decompositions

Additional rank-revealing strategies emerge if we allow general orthogonal recombinations of the A’s columns instead of just permutations. That is, we look for an orthogonal Z so that the QR factorization

$$
A Z = Q R
$$

produces a rank-revealing R. To impart the spirit of this type of matrix reduction, we show how the rank-revealing properties of a given $A Z = Q R$ factorization can be improved by replacing $Z , Q$ , and R with

$$
Z _ {\mathrm{new}} = Z Z _ {G}, \qquad Q _ {\mathrm{new}} = Q Q _ {G}, \qquad R _ {\mathrm{new}} = Q _ {G} ^ {T} R Z _ {G},
$$

respectively, where $Q _ { G }$ and $Z _ { G }$ are products of Givens rotations and $R _ { \mathrm { n e w } }$ is upper triangular. The rotations are generated by introducing zeros into a unit 2-norm nvector v which we assume approximates the n-th right singular vector of AZ. In particular, if $Z _ { \cal G } ^ { T } v = e _ { n } = I _ { n } ( : , n )$ and $\parallel R v \parallel _ { 2 } \approx \sigma _ { n }$ , then

$$
\parallel R _ {\mathrm{new}} e _ {n} \parallel_ {2} = \parallel Q _ {G} ^ {T} R Z _ {G} e _ {n} \parallel_ {2} = \parallel Q _ {G} ^ {T} R v \parallel_ {2} = \parallel R v \parallel_ {2} \approx \sigma_ {n}
$$

This says that the norm of the last column of $R _ { \mathrm { n e w } }$ is approximately the smallest singular value of A, which is certainly one way to reveal the underlying matrix rank.

We use the case $n = 4$ to illustrate how the Givens rotations arise and why the overall process is economical. Because we are transforming v to $e _ { n }$ and not $e _ { 1 }$ , we need to “flip” the mission of the 2-by-2 rotations in the $Z _ { G }$ computations so that top components are zeroed, i.e.,

$$
\left[ \begin{array}{c} 0 \\ \times \end{array} \right] = \left[ \begin{array}{c c} c & s \\ - s & c \end{array} \right] \left[ \begin{array}{c} \times \\ \times \end{array} \right].
$$

This requires only a slight modification of Algorithm 5.1.3.

In the n = 4 case we start with

$$
R = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right] \qquad v = \left[ \begin{array}{l} \times \\ \times \\ \times \\ \times \end{array} \right]
$$

and proceed to compute

$$
Z _ {G} = G _ {1 2} G _ {2 3} G _ {3 4}
$$

and

$$
Q _ {G} = H _ {1 2} H _ {2 3} H _ {3 4}
$$

as products of Givens rotations. The first step is to zero the top component of v with a “flipped” (1,2) rotation and update R accordingly:

$$
R \leftarrow R G _ {1 2} = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ \times & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right], \qquad v \leftarrow G _ {1 2} ^ {T} v = \left[ \begin{array}{l} 0 \\ \times \\ \times \\ \times \end{array} \right].
$$

To remove the unwanted subdiagonal in R, we apply a conventional (nonflipped) Givens rotation from the left to R (but not v):

$$
R \leftarrow H _ {1 2} ^ {T} R = \left[ \begin{array}{c c c c} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right], \qquad \qquad v = \left[ \begin{array}{c} 0 \\ \times \\ \times \\ \times \end{array} \right].
$$

The next step is analogous:

$$
R \leftarrow R G _ {2 3} = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right], \qquad v \leftarrow G _ {2 3} ^ {T} v = \left[ \begin{array}{l} 0 \\ 0 \\ \times \\ \times \end{array} \right].
$$

$$
R \gets H _ {2 3} ^ {T} R = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right], \qquad \qquad v = \left[ \begin{array}{l} 0 \\ 0 \\ \times \\ \times \end{array} \right].
$$

And finally,

$$
R \gets R G _ {3 4} = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \end{array} \right], \qquad v = G _ {3 4} ^ {T} v = \left[ \begin{array}{l} 0 \\ 0 \\ 0 \\ \times \end{array} \right],
$$

$$
R \gets H _ {3 4} ^ {T} R = \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \end{array} \right], \qquad \qquad v = \left[ \begin{array}{l} 0 \\ 0 \\ 0 \\ \times \end{array} \right].
$$

The pattern is clear, for $i = 1 { : } n - 1$ , a $G _ { i , i + 1 }$ is used to zero the current $v _ { i }$ and an $H _ { i , i + 1 }$ is used to zero the current $r _ { i + 1 , i } .$ . The overall transition from $\{ Q , ~ Z , ~ R \}$ to $\{ Q _ { \mathrm { n e w } } , Z _ { \mathrm { n e w } } , R _ { \mathrm { n e w } } \}$ involves $O ( m n )$ flops. If the Givens rotations are kept in factored form, this flop count is reduced to $O ( n ^ { 2 } )$ . We mention that the ideas in this subsection can be iterated to develop matrix reductions that expose the structure of matrices whose rank is less than $n - 1$ . “Zero-chasing” with Givens rotations is at the heart of many important matrix algorithms; see §6.3, §7.5, and §8.3.

# 5.4.6 The UT V Framework

As mentioned at the start of this section, we are interested in factorizations that are cheaper than the SVD but which provide the same high quality information about rank, range, and nullspace. Factorizations of this type are referred to as UTV factorizations where the $^ { 6 6 } \mathrm { T } ^ { 9 }$ stands for triangular and the $^ { 6 } U ^ { , s }$ and ${ } ^ { 6 6 } V ^ { 5 9 }$ remind us of the SVD and orthogonal U and V matrices of singular vectors.

The matrix T can be upper triangular (these are the URV factorizations) or lower triangular (these are the ULV factorizations). It turns out that in a particular application one may favor a URV approach over a ULV approach, see §6.3. Moreover, the two reductions have different approximation properties. For example, suppose $\sigma _ { k } ( A ) > \sigma _ { k + 1 } ( A )$ and S is the subspace spanned by A’s right singular vectors $v _ { k + 1 } , \ldots , v _ { n }$ . Think of S as an approximate nullspace of A. Following Stewart (1993), if

$$
U ^ {T} A V = R = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & R _ {2 2} \end{array} \right] _ {m - k} ^ {k}
$$

and $V = \left[ \left. V _ { 1 } \right| V _ { 2 } \right]$ is partitioned conformably, then

$$
\text { dist } (\text { ran } (V _ {2}), S) \leq \frac {\| R _ {1 2} \| _ {2}}{(1 - \rho_ {R} ^ {2}) \sigma_ {\min} (R _ {1 1})} \tag {5.4.10}
$$

where

$$
\rho_ {R} = \frac {\parallel R _ {2 2} \parallel_ {2}}{\sigma_ {\mathrm{min}} (R _ {1 1})}
$$

is assumed to be less than 1. On the other hand, in the ULV setting we have

$$
U ^ {T} A V = L = \left[ \begin{array}{c c} L _ {1 1} & 0 \\ L _ {2 1} & L _ {2 2} \end{array} \right] _ {m - k} ^ {k}.
$$

If $V = \left[ \left. V _ { 1 } \right| V _ { 2 } \right]$ is partitioned conformably, then

$$
\operatorname{dist} \left(\operatorname{ran} \left(V _ {2}\right), S\right) \leq \rho_ {L} \frac {\left\| L _ {1 2} \right\| _ {2}}{\left(1 - \rho_ {L} ^ {2}\right) \sigma_ {\min} \left(L _ {1 1}\right)} \tag {5.4.11}
$$

where

$$
\rho_ {L} = \frac {\parallel L _ {2 2} \parallel_ {2}}{\sigma_ {\mathrm{min}} (L _ {1 1})}
$$

is also assumed to be less than 1. However, in practice the ρ-factors in both (5.4.10) and (5.4.11) are often much less than 1. Observe that when this is the case, the upper bound in (5.4.11) is much smaller than the upper bound in (5.4.10).

# 5.4.7 Complete Orthogonal Decompositions

Related to the UTV framework is the idea of a complete orthogonal factorization. Here we compute orthogonal U and V such that

$$
U ^ {T} A V = \left[ \begin{array}{c c} T _ {1 1} & 0 \\ 0 & 0 \end{array} \right] _ {m - r} ^ {r} \tag {5.4.12}
$$

where $r = { \mathrm { r a n k } } ( A )$ . The SVD is obviously an example of a decomposition that has this structure. However, a cheaper, two-step QR process is also possible. We first use Algorithm 5.4.1 to compute

$$
U ^ {T} A \Pi = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & 0 \\ r & n - r \end{array} \right] _ {m - r} ^ {r}
$$

and then follow up with a second QR factorization

$$
Q ^ {T} \left[ \begin{array}{l} R _ {1 1} ^ {T} \\ R _ {1 2} ^ {T} \end{array} \right] = \left[ \begin{array}{l} S _ {1} \\ 0 \end{array} \right]
$$

via Algorithm 5.2.1. If we set $V = \Pi Q$ , then (5.4.12) is realized with $T _ { 1 1 } = S _ { 1 } ^ { T }$ . Note that two important subspaces are defined by selected columns of $U = [  u _ { 1 } | \cdot \cdot \cdot | u _ { m } ]$ and $V = [  v _ { 1 } | \cdot \cdot \cdot | v _ { n } ]$ :

$$
\operatorname{ran} (A) = \operatorname{span} \left\{u _ {1}, \dots , u _ {r} \right\},
$$

$$
\operatorname{null} (A) = \operatorname{span} \left\{v _ {r + 1}, \dots , v _ {n} \right\}.
$$

Of course, the computation of a complete orthogonal decomposition in practice would require the careful handling of numerical rank.

# 5.4.8 Bidiagonalization

There is one other two-sided orthogonal factorization that is important to discuss and that is the bidiagonal factorization. It is not a rank-revealing factorization per se, but it has a useful role to play because it rivals the SVD in terms of data compression.

Suppose $A \in \mathbb { R } ^ { m \times n }$ and $m \geq n$ . The idea is to compute orthogonal $U _ { B } \ ( m { \mathrm { - b y - } } m )$ and $V _ { B } \ ( n { \mathrm { - } } \mathrm { b y } { \mathrm { - } } n )$ such that

$$
U _ {B} ^ {T} A V _ {B} = \left[ \begin{array}{c c c c c} d _ {1} & f _ {1} & 0 & \dots & 0 \\ 0 & d _ {2} & f _ {2} & & 0 \\ \vdots & \ddots & \ddots & \ddots & \vdots \\ 0 & \dots & & d _ {n - 1} & f _ {n - 1} \\ 0 & \dots & & 0 & d _ {n} \\ \hline & & 0 \end{array} \right]. \tag {5.4.13}
$$

$U _ { B } = U _ { 1 } \cdot \cdot \cdot U _ { n }$ and $V _ { B } = V _ { 1 } \cdot \cdot \cdot V _ { n - 2 }$ can each be determined as a product of Householder matrices, e.g.,

$$
\left[ \begin{array}{l l l l} \times & \times & \times & \times \\ \times & \times & \times & \times \\ \times & \times & \times & \times \\ \times & \times & \times & \times \\ \times & \times & \times & \times \end{array} \right] \xrightarrow {U _ {1}} \left[ \begin{array}{l l l l} \times & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \end{array} \right] \xrightarrow {V _ {1}}
$$

$$
\left[ \begin{array}{c c c c} \times & \times & 0 & 0 \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \\ 0 & \times & \times & \times \end{array} \right] \xrightarrow {U _ {2}} \left[ \begin{array}{c c c c} \times & \times & 0 & 0 \\ 0 & \times & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \end{array} \right] \xrightarrow {V _ {2}}
$$

$$
\left[ \begin{array}{c c c c} \times & \times & 0 & 0 \\ 0 & \times & \times & 0 \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \\ 0 & 0 & \times & \times \end{array} \right] \xrightarrow {U _ {3}} \left[ \begin{array}{c c c c} \times & \times & 0 & 0 \\ 0 & \times & \times & 0 \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \\ 0 & 0 & 0 & \times \end{array} \right] \xrightarrow {U _ {4}} \left[ \begin{array}{c c c c} \times & \times & 0 & 0 \\ 0 & \times & \times & 0 \\ 0 & 0 & \times & \times \\ 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 \end{array} \right].
$$

In general, $U _ { k }$ introduces zeros into the kth column, while $V _ { k }$ zeros the appropriate entries in row k. Overall we have:

Algorithm 5.4.2 (Householder Bidiagonalization) Given $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ , the following algorithm overwrites A with $U _ { s } ^ { T } A V _ { B } = \mathbf { \dot { B } }$ where B is upper bidiagonal and $U _ { B } = U _ { 1 } \cdot \cdot \cdot U _ { n }$ and $V _ { B } = V _ { 1 } \cdots V _ { n - 2 }$ . The essential part of $U _ { j } \mathrm { ^ { 5 } s }$ Householder vector is stored in $A ( j + 1 { : } m , j )$ and the essential part of $V _ { j } \mathrm { { ^ { s } } }$ Householder vector is stored in $A ( j , j + 2 { : } n )$ .

for $j = 1:n$ [ [v, \beta] = \text{house}(A(j:m, j)) ] $A(j:m, j:n) = (I_{m-j+1} - \beta vv^T)A(j:m, j:n)$ $A(j+1:m, j) = v(2:m - j+1)$ if $j \leq n-2$ [ [v, \beta] = \text{house}(A(j, j+1:n)^T) ] $A(j:m, j+1:n) = A(j:m, j+1:n)(I_{n-j} - \beta vv^T)$ $A(j, j+2:n) = v(2:n-j)^T$ end

end

This algorithm requires $4 m n ^ { 2 } - 4 n ^ { 3 } / 3$ flops. Such a technique is used by Golub and Kahan (1965), where bidiagonalization is first described. If the matrices $U _ { B }$ and $V _ { B }$ are explicitly desired, then they can be accumulated in $4 m ^ { 2 } n - 4 n ^ { 3 } / 3$ and $4 n ^ { 3 } / 3$ flops, respectively. The bidiagonalization of A is related to the tridiagonalization of $A ^ { T } A$ . See §8.3.1.

# 5.4.9 R-Bidiagonalization

If $m \gg n$ , then a faster method of bidiagonalization method results if we upper triangularize A first before applying Algorithm 5.4.2. In particular, suppose we compute an orthogonal $Q \in \mathbb { R } ^ { m \times m }$ such that

$$
Q ^ {T} A = \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right]
$$

is upper triangular. We then bidiagonalize the square matrix $R _ { 1 }$ ,

$$
U _ {R} ^ {T} R _ {1} V _ {B} = B _ {1},
$$

where $U _ { R }$ and $V _ { B }$ are orthogonal. If $U _ { B } = Q \mathrm { d i a g } \left( U _ { R } , I _ { m - n } \right)$ , then

$$
U ^ {T} A V = \left[ \begin{array}{c} B _ {1} \\ 0 \end{array} \right] \equiv B
$$

is a bidiagonalization of A.

The idea of computing the bidiagonalization in this manner is mentioned by Lawson and Hanson (SLS, p. 119) and more fully analyzed by Chan (1982). We refer to this method as R-bidiagonalization and it requires $( 2 m n ^ { 2 } + 2 n ^ { 3 } )$ flops. This is less than the flop count for Algorithm 5.4.2 whenever $m \geq 5 n / 3$ .

# Problems

P5.4.1 Let x, $y \in \mathbb { R } ^ { m }$ and $Q \in \mathbb { R } ^ { m \times m }$ be given with Q orthogonal. Show that if

$$
Q ^ {T} x = \left[ \begin{array}{c} \alpha \\ u \end{array} \right] _ {m - 1} ^ {1}, \qquad Q ^ {T} y = \left[ \begin{array}{c} \beta \\ v \end{array} \right] _ {m - 1} ^ {1}
$$

then $u ^ { T } v = x ^ { T } y - \alpha \beta$ .

P5.4.2 Let $A = [ a _ { 1 } \mid \cdots \mid a _ { n } ] \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ be given. For any column subset $\{ a _ { c _ { 1 } } , \ldots , a _ { c _ { k } } \}$ define

$$
\operatorname{res} \left(\left[ a _ {c _ {1}} \mid \dots \mid a _ {c _ {k}} \right]\right) = \min _ {x \in \mathbf {R} ^ {k}} \| \left[ a _ {c _ {1}} \mid \dots \mid a _ {c _ {k}} \right] x - b \| _ {2}
$$

Describe an alternative pivot selection procedure for Algorithm 5.4.1 such that if $Q R \ = \ A \Pi \ =$ $\left[ \boldsymbol { a } _ { c 1 } \left| \cdots \right| \boldsymbol { a } _ { c n } \right]$ in the final factorization, then for $k = 1 { : } n { : }$

$$
\operatorname{res} \left(\left[ a _ {c _ {1}} \mid \dots \mid a _ {c _ {k}} \right]\right) = \min _ {i \geq k} \operatorname{res} \left(\left[ a _ {c _ {1}}, \ldots , a _ {c _ {k - 1}}, a _ {c _ {i}} \right]\right).
$$

P5.4.3 Suppose $T \in \mathbb { R } ^ { n \times n }$ is upper triangular and $t _ { k k } = \sigma _ { m i n } ( T )$ . Show that $T ( 1 { : } k - 1 , k ) = 0$ and $T ( k , k + 1 { : } n ) = 0$ .

P5.4.4 Suppose $A \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ . Give an algorithm that uses Householder matrices to compute an orthogonal $Q \in \mathbb { R } ^ { m \times m }$ so that if $Q ^ { T } A = L ,$ , then $L ( n + 1 { : } m , : ) = 0$ and $L ( 1 { : } n , 1 { : } n )$ is lower triangular.

P5.4.5 Suppose $R \in \mathbb { R } ^ { n \times n }$ is upper triangular and $Y \in \mathbb { R } ^ { n \times j }$ has orthonormal columns and satisfies $\| \ R Y \| _ { 2 } = \sigma ,$ . Give an algorithm that computes orthogonal U and V , each products of Givens rotations, so that $\mathbf { \dot { \zeta } } U ^ { T } R V = R _ { \mathrm { n e w } }$ is upper triangular and $V ^ { T } \check { Y } = Y _ { \mathrm { n e w } }$ has the property that

$$
Y _ {\text { new }} (n - j + 1: n,:) = \text { diag } (\pm 1).
$$

What can you say about $R _ { \mathrm { n e w } } ( n - j + 1 { : } n , n - j + 1 { : } n ) ?$

P5.4.6 Give an algorithm for reducing a complex matrix A to real bidiagonal form using complex Householder transformations.

P5.4.7 Suppose $B \in \mathbb { R } ^ { n \times n }$ is upper bidiagonal with $b _ { n n } = 0 $ . Show how to construct orthogonal U and V (product of Givens rotations) so that $U ^ { T } B V$ is upper bidiagonal with a zero nth column.

P5.4.8 Suppose $A \in \mathbb { R } ^ { m \times n }$ with $m < n ,$ . Give an algorithm for computing the factorization

$$
U ^ {T} A V = \left[ \begin{array}{c c} B & O \end{array} \right]
$$

where B is an m-by-m upper bidiagonal matrix. (Hint: Obtain the form

$$
\left[ \begin{array}{c c c c c c} \times & \times & 0 & 0 & 0 & 0 \\ 0 & \times & \times & 0 & 0 & 0 \\ 0 & 0 & \times & \times & 0 & 0 \\ 0 & 0 & 0 & \times & \times & 0 \end{array} \right].
$$

using Householder matrices and then “chase” the $( m , m + 1 )$ entry up the (m+1)st column by applying Givens rotations from the right.)

P5.4.9 Show how to efficiently bidiagonalize an n-by-n upper triangular matrix using Givens rotations.

P5.4.10 Show how to upper bidiagonalize a tridiagonal matrix $T \in \mathbb { R } ^ { n \times n }$ using Givens rotations.

P5.4.11 Show that if $B \in \mathbb { R } ^ { n \times n }$ is an upper bidiagonal matrix having a repeated singular value, then B must have a zero on its diagonal or superdiagonal.

# Notes and References for §5.4

QR with column pivoting was first discussed in:

P.A. Businger and G.H. Golub (1965). “Linear Least Squares Solutions by Householder Transformations,” Numer. Math. 7, 269–276.

In matters that concern rank deficiency, it is helpful to obtain information about the smallest singular value of the upper triangular matrix R. This can be done using the techniques of 3.5.4 or those that are discussed in:

I. Karasalo (1974). “A Criterion for Truncation of the QR Decomposition Algorithm for the Singular Linear Least Squares Problem,” BIT 14, 156–166.

N. Anderson and I. Karasalo (1975). “On Computing Bounds for the Least Singular Value of a Triangular Matrix,” BIT 15, 1–4.

C.-T. Pan and P.T.P. Tang (1999). “Bounds on Singular Values Revealed by QR Factorizations,” BIT 39, 740–756.   
C.H. Bischof (1990). “Incremental Condition Estimation,” SIAM J. Matrix Anal. Applic., 11, 312– 322.

Revealing the rank of a matrix through a carefully implementated factorization has prompted a great deal of research, see:   
T.F. Chan (1987). “Rank Revealing QR Factorizations,” Lin. Alg. Applic. 88/89, 67–82.   
T.F. Chan and P. Hansen (1992). “Some Applications of the Rank Revealing QR Factorization,” SIAM J. Sci. Stat. Comp. 13, 727–741.   
S. Chandrasekaren and I.C.F. Ipsen (1994). “On Rank-Revealing Factorizations,” SIAM J. Matrix Anal. Applic. 15, 592–622.   
M. Gu and S.C. Eisenstat (1996). “Efficient Algorithms for Computing a Strong Rank-Revealing QR Factorization,” SIAM J. Sci. Comput. 17, 848–869.   
G.W. Stewart (1999). “The QLP Approximation to the Singular Value Decomposition,” SIAM J. Sci. Comput. 20, 1336–1348.   
D.A. Huckaby and T.F. Chan (2005). “Stewart’s Pivoted QLP Decomposition for Low-Rank Matrices,” Num. Lin. Alg. Applic. 12, 153–159.   
A. Dax (2008). “Orthogonalization via Deflation: A Minimum Norm Approach to Low-Rank Approximation of a Matrix,” SIAM J. Matrix Anal. Applic. 30, 236–260.   
Z. Drma˘c and Z. Bujanovi˘c (2008). “On the Failure of Rank-Revealing QR Factorization Software—A Case Study,” ACM Trans. Math. Softw. 35, Article 12.   
We have more to say about the UTV framework in §6.5 where updating is discussed. Basic references for what we cover in this section include:   
G.W. Stewart (1993). “UTV Decompositions,” in Numerical Analysis 1993, Proceedings of the 15th Dundee Conference, June–July 1993, Longman Scientic & Technical, Harlow, Essex, UK, 225–236.   
P.A. Yoon and J.L. Barlow (1998) “An Efficient Rank Detection Procedure for Modifying the ULV Decomposition,” BIT 38, 781–801.   
J.L. Barlow, H. Erbay, and I. Slapnicar (2005). “An Alternative Algorithm for the Refinement of ULV Decompositions,” SIAM J. Matrix Anal. Applic. 27, 198–211.   
Column-pivoting makes it more difficult to achieve high performance when computing the QR factorization. However, it can be done:   
C.H. Bischof and P.C. Hansen (1992). “A Block Algorithm for Computing Rank-Revealing QR Factorizations,” Numer. Algorithms 2, 371-392.   
C.H. Bischof and G. Quintana-Orti (1998). “Computing Rank-revealing QR factorizations of Dense Matrices,” ACM Trans. Math. Softw. 24, 226–253.   
C.H. Bischof and G. Quintana-Orti (1998). “Algorithm 782: Codes for Rank-Revealing QR factorizations of Dense Matrices,” ACM Trans. Math. Softw. 24, 254–257.   
G. Quintana-Orti, X. Sun, and C.H. Bischof (1998). “A BLAS–3 Version of the QR Factorization with Column Pivoting,” SIAM J. Sci. Comput. 19, 1486–1494.   
A carefully designed LU factorization can also be used to shed light on matrix rank:   
T-M. Hwang, W-W. Lin, and E.K. Yang (1992). “Rank-Revealing LU Factorizations,” Lin. Alg. Applic. 175, 115–141.   
T.-M. Hwang, W.-W. Lin and D. Pierce (1997). “Improved Bound for Rank Revealing LU Factorizations,” Lin. Alg. Applic. 261, 173–186.   
L. Miranian and M. Gu (2003). “Strong Rank Revealing LU Factorizations,” Lin. Alg. Applic. 367, 1–16.   
Column pivoting can be incorporated into the modified Gram-Schmidt process, see:   
A. Dax (2000). “A Modified Gram-Schmidt Algorithm with Iterative Orthogonalization and Column Pivoting,” Lin. Alg. Applic. 310, 25–42.   
M. Wei and Q. Liu (2003). “Roundoff Error Estimates of the Modified GramSchmidt Algorithm with Column Pivoting,” BIT 43, 627–645.   
Aspects of the complete orthogonal decomposition are discussed in:

R.J. Hanson and C.L. Lawson (1969). “Extensions and Applications of the Householder Algorithm for Solving Linear Least Square Problems,” Math. Comput. 23, 787–812.   
P.A. Wedin (1973). “On the Almost Rank-Deficient Case of the Least Squares Problem,” BIT 13, 344–354.   
G.H. Golub and V. Pereyra (1976). “Differentiation of Pseudo-Inverses, Separable Nonlinear Least Squares Problems and Other Tales,” in Generalized Inverses and Applications, M.Z. Nashed (ed.), Academic Press, New York, 303–324.   
The quality of the subspaces that are exposed through a complete orthogonal decomposition are analyzed in:   
R.D. Fierro and J.R. Bunch (1995). “Bounding the Subspaces from Rank Revealing Two-Sided Orthogonal Decompositions,” SIAM J. Matrix Anal. Applic. 16, 743–759.   
R.D. Fierro (1996). “Perturbation Analysis for Two-Sided (or Complete) Orthogonal Decompositions,” SIAM J. Matrix Anal. Applic. 17, 383–400.   
The bidiagonalization is a particularly important decomposition because it typically precedes the computation of the SVD as we discuss in §8.6. Thus, there has been a strong research interest in its efficient and accurate computation:   
B. Lang (1996). “Parallel Reduction of Banded Matrices to Bidiagonal Form,” Parallel Comput. 22, 1–18.   
J.L. Barlow (2002). “More Accurate Bidiagonal Reduction for Computing the Singular Value Decomposition,” SIAM J. Matrix Anal. Applic. 23, 761–798.   
J.L. Barlow, N. Bosner and Z. Drma˘c (2005). “A New Stable Bidiagonal Reduction Algorithm,” Lin. Alg. Applic. 397, 35–84.   
B.N. Parlett (2005). “A Bidiagonal Matrix Determines Its Hyperbolic SVD to Varied Relative Accuracy,” SIAM J. Matrix Anal. Applic. 26, 1022–1057.   
N. Bosner and J.L. Barlow (2007). “Block and Parallel Versions of One-Sided Bidiagonalization,” SIAM J. Matrix Anal. Applic. 29, 927–953.   
G.W. Howell, J.W. Demmel, C.T. Fulton, S. Hammarling, and K. Marmol (2008). “Cache Efficient Bidiagonalization Using BLAS 2.5 Operators,” ACM Trans. Math. Softw. 34, Article 14.   
H. Ltaief, J. Kurzak, and J. Dongarra (2010). “Parallel Two-Sided Matrix Reduction to Band Bidiagonal Form on Multicore Architectures,” IEEE Trans. Parallel Distrib. Syst. 21, 417–423.

# 5.5 The Rank-Deficient Least Squares Problem

If A is rank deficient, then there are an infinite number of solutions to the LS problem. We must resort to techniques that incorporate numerical rank determination and identify a particular solution as “special.” In this section we focus on using the SVD to compute the minimum norm solution and QR-with-column-pivoting to compute what is called the basic solution. Both of these approaches have their merits and we conclude with a subset selection procedure that combines their positive attributes.

# 5.5.1 The Minimum Norm Solution

Suppose $A \in \mathbb { R } ^ { m \times n }$ and rank $( A ) = r < n$ . The rank-deficient LS problem has an infinite number of solutions, for if x is a minimizer and $z \in \mathsf { n u l l } ( A )$ , then $x + z$ is also a minimizer. The set of all minimizers

$$
\mathcal {X} = \left\{x \in \mathbb {R} ^ {n}: \| A x - b \| _ {2} = \min \right\}
$$

is convex and so if $x _ { 1 } , x _ { 2 } \in \mathcal { X }$ and $\lambda \in [ 0 , 1 ]$ , then

$$
\| A (\lambda x _ {1} + (1 - \lambda) x _ {2}) - b \| _ {2} \leq \lambda \| A x _ {1} - b \| _ {2} + (1 - \lambda) \| A x _ {2} - b \| _ {2} = \min _ {x \in \mathbb {R} ^ {n}} \| A x - b \| _ {2}.
$$

Thus, $\lambda x _ { 1 } + ( 1 - \lambda ) x _ { 2 } \in \mathcal { X }$ . It follows that X has a unique element having minimum 2-norm and we denote this solution by $x _ { L S }$ . (Note that in the full-rank case, there is only one LS solution and so it must have minimal 2-norm. Thus, we are consistent with the notation in §5.3.)

Any complete orthogonal factorization (§5.4.7) can be used to compute $x _ { L S }$ . In particular, if Q and $Z$ are orthogonal matrices such that

$$
Q ^ {T} A Z = T = \left[ \begin{array}{c c} T _ {1 1} & 0 \\ 0 & 0 \end{array} \right] _ {m - r} ^ {r}, r = \operatorname{rank} (A)
$$

then

$$
\| A x - b \| _ {2} ^ {2} = \| (Q ^ {T} A Z) Z ^ {T} x - Q ^ {T} b \| _ {2} ^ {2} = \| T _ {1 1} w - c \| _ {2} ^ {2} + \| d \| _ {2} ^ {2}
$$

where

$$
Z ^ {T} x = \left[ \begin{array}{l} w \\ y \end{array} \right] _ {n - r} ^ {r}, \qquad Q ^ {T} b = \left[ \begin{array}{l} c \\ d \end{array} \right] _ {m - r} ^ {r}.
$$

Clearly, if x is to minimize the sum of squares, then we must have $w = T _ { 1 1 } ^ { - 1 } c$ . For x to have minimal 2-norm, y must be zero, and thus

$$
x _ {L S} = Z \left[ \begin{array}{c} T _ {1 1} ^ {- 1} c \\ 0 \end{array} \right].
$$

Of course, the SVD is a particularly revealing complete orthogonal decomposition. It provides a neat expression for $x _ { L S }$ and the norm of the minimum residual $\rho _ { L S } =$ $\parallel A x _ { L S } - b \parallel _ { 2 } .$ .

Theorem 5.5.1. Suppose $U ^ { T } A V = \Sigma$ is the SVD of $A \in \mathbb { R } ^ { m \times n }$ with $r = r a n k ( A )$ . If $U = [  u _ { 1 } | \cdot \cdot \cdot | u _ { m } ]$ and $V = [  v _ { 1 } | \cdot \cdot \cdot | v _ { n } ]$ are column partitionings and $b \in \mathbb { R } ^ { m }$ , then

$$
x _ {L S} = \sum_ {i = 1} ^ {r} \frac {u _ {i} ^ {T} b}{\sigma_ {i}} v _ {i} \tag {5.5.1}
$$

minimizes $\parallel A x - b \parallel _ { 2 }$ and has the smallest 2-norm of all minimizers. Moreover

$$
\rho_ {L S} ^ {2} = \| A x _ {L S} - b \| _ {2} ^ {2} = \sum_ {i = r + 1} ^ {m} (u _ {i} ^ {T} b) ^ {2}. \tag {5.5.2}
$$

Proof. For any $\boldsymbol { x } \in \mathbb { R } ^ { n }$ we have

$$
\begin{array}{l} \| A x - b \| _ {2} ^ {2} = \| (U ^ {T} A V) (V ^ {T} x) - U ^ {T} b \| _ {2} ^ {2} = \| \Sigma \alpha - U ^ {T} b \| _ {2} ^ {2} \\ = \sum_ {i = 1} ^ {r} (\sigma_ {i} \alpha_ {i} - u _ {i} ^ {T} b) ^ {2} + \sum_ {i = r + 1} ^ {m} (u _ {i} ^ {T} b) ^ {2}, \\ \end{array}
$$

where $\alpha = V ^ { T } x$ . Clearly, if x solves the LS problem, then $\alpha _ { i } = ( u _ { i } ^ { T } b / \sigma _ { i } )$ for $i = 1 { : } r$ . If we set $\alpha ( r + 1 { : } n ) = 0$ , then the resulting x has minimal 2-norm.

# 5.5.2 A Note on the Pseudoinverse

If we define the matrix $A ^ { + } \in \mathbb { R } ^ { n \times m }$ by $A ^ { + } = V \Sigma ^ { + } U ^ { T }$ where

$$
\Sigma^ {+} = \mathrm{diag} \left(\frac {1}{\sigma_ {1}}, \ldots , \frac {1}{\sigma_ {r}}, 0, \ldots , 0\right) \in \mathbb {R} ^ {n \times m}, \qquad r = \mathsf {r a n k} (A),
$$

then $x _ { L S } = A ^ { + } b$ and $\rho _ { L S } = \parallel ( I - A A ^ { + } ) b \parallel _ { 2 }$ . A+ is referred to as the pseudo-inverse of A. It is the unique minimal Frobenius norm solution to the problem

$$
\min _ {X \in \mathbb {R} ^ {m \times n}} \| A X - I _ {m} \| _ {F}. \tag {5.5.3}
$$

If rank(A) = n, then $A ^ { + } = ( A ^ { T } A ) ^ { - 1 } A ^ { T }$ , while if $m = n = \operatorname { r a n k } ( A )$ , then $A ^ { + } = A ^ { - 1 }$ . Typically, A+ is defined to be the unique matrix $\boldsymbol { X } \in \mathbb { R } ^ { n \times m }$ that satisfies the four Moore-Penrose conditions:

These conditions amount to the requirement that $A A ^ { + }$ and $A ^ { + } A$ be orthogonal projections onto ran(A) and $\mathsf { r a n } ( A ^ { T } )$ , respectively. Indeed,

$$
A A ^ {+} = U _ {1} U _ {1} ^ {T}
$$

where $U _ { 1 } = U ( 1 { : } m , 1 { : } r )$ and

$$
A ^ {+} A = V _ {1} V _ {1} ^ {T}
$$

where $V _ { 1 } = V ( 1 { : } n , 1 { : } r )$ .

# 5.5.3 Some Sensitivity Issues

In §5.3 we examined the sensitivity of the full-rank LS problem. The behavior of $x _ { L S }$ in this situation is summarized in Theorem 5.3.1. If we drop the full-rank assumption, then $x _ { L S }$ is not even a continuous function of the data and small changes in A and b can induce arbitrarily large changes in $x _ { L S } = A ^ { + } b$ . The easiest way to see this is to consider the behavior of the pseudoinverse. If A and δA are in $\mathbb { R } ^ { m \times n }$ , then Wedin (1973) and Stewart (1975) show that

$$
\left\| (A + \delta A) ^ {+} - A ^ {+} \right\| _ {F} \leq 2 \| \delta A \| _ {F} \max \left\{\left\| A ^ {+} \right\| _ {2} ^ {2}, \left\| (A + \delta A) ^ {+} \right\| _ {2} ^ {2} \right\}.
$$

This inequality is a generalization of Theorem 2.3.4 in which perturbations in the matrix inverse are bounded. However, unlike the square nonsingular case, the upper bound does not necessarily tend to zero as δA tends to zero. If

$$
A = \left[ \begin{array}{l l} 1 & 0 \\ 0 & 0 \\ 0 & 0 \end{array} \right] \qquad \text {and} \qquad \delta A = \left[ \begin{array}{l l} 0 & 0 \\ 0 & \epsilon \\ 0 & 0 \end{array} \right]
$$

then

$$
A ^ {+} = \left[ \begin{array}{l l l} 1 & 0 & 0 \\ 0 & 0 & 0 \end{array} \right] \qquad \text {and} \qquad (A + \delta A) ^ {+} = \left[ \begin{array}{l l l} 1 & 0 & 0 \\ 1 & 1 / \epsilon & 0 \end{array} \right],
$$

and

$$
\| A ^ {+} - (A + \delta A) ^ {+} \| _ {2} = 1 / \epsilon .
$$

The numerical determination of an LS minimizer in the presence of such discontinuities is a major challenge.

# 5.5.4 The Truncated SVD Solution

Suppose ${ \widehat { U } } , { \widehat { \Sigma } } ,$ and $\widehat { V }$ are the computed SVD factors of a matrix A and $\hat { r }$ is accepted as its δ-rank, i.e.,

$$
\hat {\sigma} _ {n} \leq \dots \leq \hat {\sigma} _ {\hat {r}} \leq \delta <   \hat {\sigma} _ {\hat {r}} \leq \dots \leq \hat {\sigma} _ {1}.
$$

It follows that we can regard

$$
x _ {\hat {r}} = \sum_ {i = 1} ^ {\hat {r}} \frac {\hat {u} _ {i} ^ {T} b}{\hat {\sigma} _ {i}} \hat {v} _ {i}
$$

as an approximation to $x _ { L S }$ . Since $\Vert \ x _ { \hat { r } } \ \Vert _ { 2 } \approx 1 / \sigma _ { \hat { r } } \leq 1 / \delta$ , then $\delta$ may also be chosen with the intention of producing an approximate LS solution with suitably small norm. In $\ S 6 . 2 . 1$ , we discuss more sophisticated methods for doing this.

If $\hat { \sigma } _ { \hat { r } } \gg \delta$ , then we have reason to be comfortable with $x _ { \hat { r } }$ because A can then be unambiguously regarded as a $\mathrm { r a n k } ( A _ { \hat { r } } )$ matrix (modulo δ).

On the other hand, $\{ \hat { \sigma } _ { 1 } , \hdots , \hat { \sigma } _ { n } \}$ might not clearly split into subsets of small and large singular values, making the determination of $\hat { r }$ by this means somewhat arbitrary. This leads to more complicated methods for estimating rank, which we now discuss in the context of the LS problem. The issues are readily communicated by making two simplifying assumptions. Assume that $r ~ = ~ n$ , and that $\Delta A = 0$ in (5.4.4), which implies that $\boldsymbol { W } ^ { T } \boldsymbol { A } \boldsymbol { Z } = \widehat { \boldsymbol { \Sigma } } = \boldsymbol { \Sigma }$ is the SVD. Denote the ith columns of the matrices $\widehat { U }$ , $W , \widehat { V }$ , and $Z$ by $\hat { u } _ { i } , w _ { i } , \hat { v } _ { i }$ , and $z _ { i }$ , respectively. Because

$$
\begin{array}{l} x _ {L S} - x _ {\hat {r}} = \sum_ {i = 1} ^ {n} \frac {w _ {i} ^ {T} b}{\sigma_ {i}} z _ {i} - \sum_ {i = 1} ^ {\hat {r}} \frac {\hat {u} _ {i} ^ {T} b}{\sigma_ {i}} \hat {v} _ {i} \\ = \sum_ {i = 1} ^ {\hat {r}} \frac {((w _ {i} - \hat {u} _ {i}) ^ {T} b) z _ {i} + (\hat {u} _ {i} ^ {T} b) (z _ {i} - \hat {v} _ {i})}{\sigma_ {i}} + \sum_ {i = \hat {r} + 1} ^ {n} \frac {w _ {i} ^ {T} b}{\sigma_ {i}} z _ {i} \\ \end{array}
$$

it follows from $\begin{array} { r } { \| w _ { i } - \hat { u } _ { i } \| _ { 2 } \leq \epsilon , \| \hat { u } _ { i } \| _ { 2 } \leq 1 + \epsilon , } \end{array}$ and $\| z _ { i } - \hat { v } _ { i } \| _ { 2 } \leq \epsilon$ that

$$
\parallel x _ {\hat {r}} - x _ {L S} \parallel_ {2} \leq \frac {\hat {r}}{\sigma_ {\hat {r}}} 2 (1 + \epsilon) \epsilon \parallel b \parallel_ {2} + \sqrt {\sum_ {i = \hat {r} + 1} ^ {n} \left(\frac {w _ {i} ^ {T} b}{\sigma_ {i}}\right) ^ {2}}.
$$

The parameter $\hat { r }$ can be determined as that integer which minimizes the upper bound. Notice that the first term in the bound increases with ${ \hat { r } } ,$ , while the second decreases.

On occasions when minimizing the residual is more important than accuracy in the solution, we can determine $\hat { r }$ on the basis of how close we surmise $\parallel b - A x _ { \hat { r } } \parallel _ { 2 }$ is to the true minimum. Paralleling the above analysis, it can be shown that

$$
\| b - A x _ {\hat {r}} \| _ {2} \leq \| b - A x _ {L S} \| _ {2} + (n - \hat {r}) \| b \| _ {2} + \epsilon \hat {r} \| b \| _ {2} \left(1 + (1 + \epsilon) \frac {\hat {\sigma} _ {1}}{\hat {\sigma} _ {\hat {r}}}\right).
$$

Again $\hat { r }$ could be chosen to minimize the upper bound. See Varah (1973) for practical details and also LAPACK.

# 5.5.5 Basic Solutions via QR with Column Pivoting

Suppose $A \in \mathbb { R } ^ { m \times n }$ has rank r. QR with column pivoting (Algorithm 5.4.1) produces the factorization $A \Pi = Q R$ where

$$
R = \left[ \begin{array}{c c} R _ {1 1} & R _ {1 2} \\ 0 & 0 \\ r & n - r \end{array} \right] _ {m - r} ^ {r}.
$$

Given this reduction, the LS problem can be readily solved. Indeed, for any $\boldsymbol { x } \in \mathbb { R } ^ { n }$ we have

$$
\parallel A x - b \parallel_ {2} ^ {2} = \parallel (Q ^ {T} A \Pi) (\Pi^ {T} x) - (Q ^ {T} b) \parallel_ {2} ^ {2} = \parallel R _ {1 1} y - (c - R _ {1 2} z) \parallel_ {2} ^ {2} + \parallel d \parallel_ {2} ^ {2},
$$

where

$$
\Pi^ {T} x = \left[ \begin{array}{l} y \\ z \end{array} \right] _ {n - r} ^ {r} \quad \text {and} \quad Q ^ {T} b = \left[ \begin{array}{l} c \\ d \end{array} \right] _ {m - r} ^ {r}.
$$

Thus, if x is an LS minimizer, then we must have

$$
x = \Pi \left[ \begin{array}{c} R _ {1 1} ^ {- 1} (c - R _ {1 2} z) \\ z \end{array} \right].
$$

$\operatorname { I f } z$ is set to zero in this expression, then we obtain the basic solution

$$
x _ {B} = \Pi \left[ \begin{array}{c} R _ {1 1} ^ {- 1} c \\ 0 \end{array} \right].
$$

Notice that $x _ { B }$ has at most r nonzero components and so $A x _ { B }$ involves a subset of A’s columns.

The basic solution is not the minimal 2-norm solution unless the submatrix $R _ { 1 2 }$ is zero since

$$
\left\| x _ {L S} \right\| _ {2} = \min _ {z \in \mathbb {R} ^ {n - 2}} \left\| x _ {B} - \Pi \left[ \begin{array}{c} R _ {1 1} ^ {- 1} R _ {1 2} \\ - I _ {n - r} \end{array} \right] z \right\| _ {2}. \tag {5.5.4}
$$

Indeed, this characterization of $\Vert { x } _ { L S } \ \Vert _ { 2 }$ can be used to show that

$$
1 \leq \frac {\| x _ {B} \| _ {2}}{\| x _ {L S} \| _ {2}} \leq \sqrt {1 + \| R _ {1 1} ^ {- 1} R _ {1 2} \| _ {2} ^ {2}}. \tag {5.5.5}
$$

See Golub and Pereyra (1976) for details.

# 5.5.6 Some Comparisons

As we mentioned, when solving the LS problem via the SVD, only Σ and V have to be computed assuming that the right hand side b is available. The table in Figure 5.5.1 compares the flop efficiency of this approach with the other algorithms that we have presented.

<table><tr><td>LS Algorithm</td><td>Flop Count</td></tr><tr><td>Normal equations</td><td> $mn^{2} + n^{3}/3$ </td></tr><tr><td>Householder QR</td><td> $n^{3}/3$ </td></tr><tr><td>Modified Gram-Schmidt</td><td> $2mn^{2}$ </td></tr><tr><td>Givens QR</td><td> $3mn^{2} - n^{3}$ </td></tr><tr><td>Householder Bidiagonalization</td><td> $4mn^{2} - 2n^{3}$ </td></tr><tr><td>R-Bidiagonalization</td><td> $2mn^{2} + 2n^{3}$ </td></tr><tr><td>SVD</td><td> $4mn^{2} + 8n^{3}$ </td></tr><tr><td>R-SVD</td><td> $2mn^{2} + 11n^{3}$ </td></tr></table>

Figure 5.5.1. Flops associated with various least squares methods

# 5.5.7 SVD-Based Subset Selection

Replacing A by $A _ { \tilde { r } }$ in the LS problem amounts to filtering the small singular values and can make a great deal of sense in those situations where A is derived from noisy data. In other applications, however, rank deficiency implies redundancy among the factors that comprise the underlying model. In this case, the model-builder may not be interested in a predictor such as $A _ { \tilde { r } } x _ { \tilde { r } }$ that involves all n redundant factors. Instead, a predictor $A y$ may be sought where y has at most ˜r nonzero components. The position of the nonzero entries determines which columns of A, i.e., which factors in the model, are to be used in approximating the observation vector b. How to pick these columns is the problem of subset selection.

QR with column pivoting is one way to proceed. However, Golub, Klema, and Stewart (1976) have suggested a technique that heuristically identifies a more independent set of columns than are involved in the predictor $A x _ { B }$ . The method involves both the SVD and QR with column pivoting:

Step 1. Compute the SVD $A = U \Sigma V ^ { T }$ and use it to determine a rank estimate ˜r.

Step 2. Calculate a permutation matrix P such that the columns of the matrix $B _ { 1 } \in \mathbb { R } ^ { m \times \tilde { r } }$ in $A P = \left[ B _ { 1 } \mid B _ { 2 } \right]$ are “sufficiently independent.”

Step 3. Predict b with Ay where $\boldsymbol { y } = \boldsymbol { P } \left[ \begin{array} { l } { z } \\ { 0 } \end{array} \right]$ z and $z \in \mathbb { R } ^ { \tilde { r } }$ minimizes $\parallel B _ { 1 } z - b \parallel _ { 2 }$ .

The second step is key. Because

$$
\min _ {z \in \mathbf {R} ^ {\tilde {r}}} \| B _ {1} z - b \| _ {2} = \| A y - b \| _ {2} \geq \min _ {x \in \mathbf {R} ^ {n}} \| A x - b \| _ {2}
$$

it can be argued that the permutation P should be chosen to make the residual $r =$ $( I - B _ { 1 } B _ { 1 } ^ { + } ) b$ as small as possible. Unfortunately, such a solution procedure can be

unstable. For example, if

$$
A = \left[ \begin{array}{c c c} 1 & 1 & 0 \\ 1 & 1 + \epsilon & 1 \\ 0 & 0 & 1 \end{array} \right], \qquad b = \left[ \begin{array}{c} 1 \\ - 1 \\ 0 \end{array} \right],
$$

$\tilde { r } = 2$ , and $P = I$ , then min $\| \ B _ { 1 } z - b \| _ { 2 } = 0$ , but $\parallel B _ { 1 } ^ { + } b \parallel _ { 2 } = { \cal O } ( 1 / \epsilon )$ . On the other hand, any proper subset involving the third column of A is strongly independent but renders a much larger residual.

This example shows that there can be a trade-off between the independence of the chosen columns and the norm of the residual that they render. How to proceed in the face of this trade-off requires useful bounds on $\sigma _ { \tilde { r } } ( B _ { 1 } )$ , the smallest singular value of $B _ { 1 }$ .

Theorem 5.5.2. Let the SVD of $A \in \mathbb { R } ^ { m \times n }$ be given by $U ^ { T } A V = \Sigma = \mathrm { d i a g } ( \sigma _ { i } )$ and define the matrix $B _ { 1 } \in \mathbb { R } ^ { m \times \tilde { r } } , \tilde { r } \leq \mathsf { r a n k } ( A )$ , by

$$
A P = \left[ \begin{array}{c c} B _ {1} & B _ {2} \\ \tilde {r} & n - \tilde {r} \end{array} \right]
$$

where $P \in \mathbb { R } ^ { n \times n }$ is a permutation. If

$$
P ^ {T} V = \left[ \begin{array}{c c} \widetilde {V} _ {1 1} & \widetilde {V} _ {1 2} \\ \widetilde {V} _ {2 1} & \widetilde {V} _ {2 2} \end{array} \right] _ {n - \tilde {r}} ^ {\tilde {r}} \tag {5.5.6}
$$

and $\tilde { V } _ { 1 1 }$ is nonsingular, then

$$
\frac {\sigma_ {\tilde {r}} (A)}{\| \tilde {V} _ {1 1} ^ {- 1} \| _ {2}} \leq \sigma_ {\tilde {r}} (B _ {1}) \leq \sigma_ {\tilde {r}} (A).
$$

Proof. The upper bound follows from Corollary 2.4.4. To establish the lower bound, partition the diagonal matrix of singular values as follows:

$$
\Sigma = \left[ \begin{array}{c c} \Sigma_ {1} & 0 \\ 0 & \Sigma_ {2} \end{array} \right] _ {m - \tilde {r}} ^ {\tilde {r}}.
$$

If $\boldsymbol { w } \in \mathbb { R } ^ { \tilde { r } }$ is a unit vector with the property that $\parallel B _ { 1 } w \parallel _ { 2 } = \sigma _ { \tilde { r } } ( B _ { 1 } )$ , then

$$
\sigma_ {\tilde {r}} (B _ {1}) ^ {2} = \parallel B _ {1} w \parallel_ {2} ^ {2} = \left\| U \Sigma V ^ {T} P \left[ \begin{array}{c} w \\ 0 \end{array} \right] \right\| _ {2} ^ {2} = \parallel \Sigma_ {1} \widetilde {V} _ {1 1} ^ {T} w \parallel_ {2} ^ {2} + \parallel \Sigma_ {2} \widetilde {V} _ {1 2} ^ {T} w \parallel_ {2} ^ {2}.
$$

The theorem now follows because $\begin{array} { r } { \| \Sigma _ { 1 } \widetilde { V } _ { 1 1 } ^ { T } w \| _ { 2 } \geq \sigma _ { \widetilde { r } } ( A ) / \| \widetilde { V } _ { 1 1 } ^ { - 1 } \| _ { 2 } . } \end{array}$

This result suggests that in the interest of obtaining a sufficiently independent subset of columns, we choose the permutation P such that the resulting $\widetilde { V } _ { 1 1 }$ submatrix is as well-conditioned as possible. A heuristic solution to this problem can be obtained by computing the QR with column-pivoting factorization of the matrix $\left[ V _ { 1 1 } ^ { T } ~ V _ { 2 1 } ^ { T } \right]$ , where

$$
V = \left[ \begin{array}{c c} V _ {1 1} & V _ {1 2} \\ V _ {2 1} & V _ {2 2} \end{array} \right] _ {n - \tilde {r}} ^ {\tilde {r}}
$$

is a partitioning of the matrix V , A’s matrix of right singular vectors. In particular, if we apply QR with column pivoting (Algorithm 5.4.1) to compute

$$
Q ^ {T} [ V _ {1 1} ^ {T} V _ {2 1} ^ {T} ] P = [ R _ {1 1} | R _ {1 2} ] _ {\tilde {r} n - \tilde {r}}
$$

where $Q$ is orthogonal, P is a permutation matrix, and $R _ { 1 1 }$ is upper triangular, then (5.5.6) implies

$$
\left[ \begin{array}{l} \widetilde {V} _ {1 1} \\ \widetilde {V} _ {2 1} \end{array} \right] = P ^ {T} \left[ \begin{array}{l} V _ {1 1} \\ V _ {2 1} \end{array} \right] = \left[ \begin{array}{l} R _ {1 1} ^ {T} Q ^ {T} \\ R _ {1 2} ^ {T} Q ^ {T} \end{array} \right].
$$

Note that $R _ { 1 1 }$ is nonsingular and that $\parallel \widetilde { V } _ { 1 1 } ^ { - 1 } \parallel _ { 2 } = \parallel R _ { 1 1 } ^ { - 1 } \parallel _ { 2 }$ . Heuristically, column pivoting tends to produce a well-conditioned $R _ { 1 1 }$ , and so the overall process tends to produce a well-conditioned $\widetilde { V } _ { 1 1 }$ .

Algorithm 5.5.1 Given $A \in \mathbb { R } ^ { m \times n }$ and $b \in \mathbb { R } ^ { m }$ the following algorithm computes a permutation $P ,$ , a rank estimate ${ \tilde { r } } _ { : }$ and a vector $z \in \mathbb { R } ^ { \tilde { r } }$ such that the first $\tilde { r }$ columns of $B = A P$ are independent and $\parallel B ( : , 1 : \tilde { r } ) z - b \parallel _ { 2 }$ is minimized.

Compute the SVD $U ^ { T } A V = \mathrm { d i a g } ( \sigma _ { 1 } , . . . , \sigma _ { n } )$ and save $V .$ .

Determine $\tilde { r } \le \mathsf { r a n k } ( A )$ .

Apply QR with column pivoting: $Q ^ { T } V ( : , 1 : \tilde { r } ) ^ { T } P = [ R _ { 1 1 } | R _ { 1 2 } ]$ and set

$$
A P = \left[ B _ {1} \mid B _ {2} \right] \text {with} B _ {1} \in \mathbb {R} ^ {m \times \tilde {r}} \text {and} B _ {2} \in \mathbb {R} ^ {m \times (n - \tilde {r})}.
$$

Determine $z \in \mathbb { R } ^ { \tilde { r } }$ such that $\parallel b - B _ { 1 } z \parallel _ { 2 } = \operatorname* { m i n }$ .

# 5.5.8 Column Independence Versus Residual Size

We return to the discussion of the trade-off between column independence and norm of the residual. In particular, to assess the above method of subset selection we need to examine the residual of the vector y that it produces

$$
r _ {y} = b - A y = b - B _ {1} z = (I - B _ {1} B _ {1} ^ {+}) b.
$$

Here, $B _ { 1 } = B ( : , 1 { : } \tilde { r } )$ with $B = A P$ . To this end, it is appropriate to compare $r _ { y }$ with

$$
r _ {x _ {\tilde {r}}} = b - A x _ {\tilde {r}}
$$

since we are regarding A as a rank-˜r matrix and since $x _ { \tilde { r } }$ solves the nearest rank-˜r LS problem min $\parallel A _ { \tilde { r } } x - b \parallel _ { 2 }$ .

Theorem 5.5.3. Assume that $U ^ { T } A V = \Sigma$ is the SVD of $A \in \mathbb { R } ^ { m \times n }$ and that $r _ { y }$ and $r _ { x _ { \tilde { r } } }$ are defined as above. $I f \tilde { V } _ { 1 1 }$ is the leading r-by-r principal submatrix of $P ^ { T } V$ , then

$$
\| r _ {x _ {\tilde {r}}} - r _ {y} \| _ {2} \leq \frac {\sigma_ {\tilde {r} + 1} (A)}{\sigma_ {\tilde {r}} (A)} \| \tilde {V} _ {1 1} ^ {- 1} \| _ {2} \| b \| _ {2}.
$$

Proof. Note that $r _ { x _ { \tilde { r } } } = ( I - U _ { 1 } U _ { 1 } ^ { T } ) l$ and $r _ { y } = ( I - Q _ { 1 } Q _ { 1 } ^ { T } ) b$ where

$$
U = \left[ \begin{array}{c c} U _ {1} & U _ {2} \\ \tilde {r} & m - \tilde {r} \end{array} \right]
$$

is a partitioning of the matrix U and $Q _ { 1 } = B _ { 1 } ( B _ { 1 } ^ { T } B _ { 1 } ) ^ { - 1 / 2 }$ . Using Theorem 2.6.1 we obtain

$$
\left\| r _ {x _ {\tilde {r}}} - r _ {y} \right\| _ {2} \leq \left\| U _ {1} U _ {1} ^ {T} - Q _ {1} Q _ {1} ^ {T} \right\| _ {2} \left\| b \right\| _ {2} = \left\| U _ {2} ^ {T} Q _ {1} \right\| _ {2} \left\| b \right\| _ {2}
$$

while Theorem 5.5.2 permits us to conclude

$$
\begin{array}{l} \| U _ {2} ^ {T} Q _ {1} \| _ {2} \leq \| U _ {2} ^ {T} B _ {1} \| _ {2} \| (B _ {1} ^ {T} B _ {1}) ^ {- 1 / 2} \| _ {2} \\ \leq \sigma_ {\tilde {r} + 1} (A) \frac {1}{\sigma_ {\tilde {r}} (B _ {1})} \leq \frac {\sigma_ {\tilde {r} + 1} (A)}{\sigma_ {\tilde {r}} (A)} \| \tilde {V} _ {1 1} ^ {- 1} \| _ {2}, \\ \end{array}
$$

and this establishes the theorem.

Noting that

$$
\left\| r _ {x _ {\bar {r}}} - r _ {y} \right\| _ {2} = \left\| B _ {1} y - \sum_ {i = 1} ^ {r} (u _ {i} ^ {T} b) u _ {i} \right\| _ {2}
$$

we see that Theorem 5.5.3 sheds light on how well $B _ { 1 } y$ can predict the “stable” component of b, i.e., $U _ { 1 } ^ { T } b$ . Any attempt to approximate $U _ { 2 } ^ { T } b$ can lead to a large norm solution. Moreover, the theorem says that if $\sigma _ { \tilde { r } + 1 } ( A ) \ll \sigma _ { \tilde { r } } ( A )$ , then any reasonably independent subset of columns produces essentially the same-sized residual. On the other hand, if there is no well-defined gap in the singular values, then the determination of $\tilde { r }$ becomes difficult and the entire subset selection problem becomes more complicated.

# Problems

P5.5.1 Show that if

$$
A = \left[ \begin{array}{c c} T & S \\ 0 & 0 \end{array} \right] _ {m - r} ^ {r}
$$

where $r = { \mathsf { r a n k } } ( A )$ and T is nonsingular, then

$$
X = \left[ \begin{array}{c c} T ^ {- 1} & 0 \\ 0 & 0 \\ r & m - r \end{array} \right] _ {n - r} ^ {r}
$$

satisfies $A X A = A$ and $( A X ) ^ { T } = ( A X )$ . In this case, we say that X is a (1,3) pseudoinverse of A. Show that for general A, $x _ { B } = X b$ where X is a (1,3) pseudoinverse of A.

P5.5.2 Define $B ( \lambda ) \in \mathbb { R } ^ { n \times m }$ by

$$
B (\lambda) = (A ^ {T} A + \lambda I) ^ {- 1} A ^ {T}
$$

where $\lambda > 0$ . Show that

$$
\parallel B (\lambda) - A ^ {+} \parallel_ {2} = \frac {\lambda}{\sigma_ {r} (A) [ \sigma_ {r} (A) ^ {2} + \lambda ]}, \qquad r = \operatorname{rank} (A),
$$

and therefore that $B ( \lambda ) \to A ^ { + }$ as $\lambda  0$ .

P5.5.3 Consider the rank-deficient LS problem

$$
\min _ {y \in \mathbf {R} ^ {r},   z \in \mathbf {R} ^ {n - r}} \left\| \left[ \begin{array}{c c} R & S \\ 0 & 0 \end{array} \right] \left[ \begin{array}{c} y \\ z \end{array} \right] - \left[ \begin{array}{c} c \\ d \end{array} \right] \right\| _ {2}
$$

where $R \in \mathbb R ^ { r \times r } , S \in \mathbb R ^ { r \times n - r } , y \in \mathbb R ^ { r }$ , and $z \in \mathbb { R } ^ { n - r }$ . Assume that R is upper triangular and nonsingular. Show how to obtain the minimum norm solution to this problem by computing an appropriate QR factorization without pivoting and then solving for the appropriate y and $z .$

P5.5.4 Show that if $A _ { k } \to A$ and $A _ { k } ^ { + } \to A ^ { + }$ , then there exists an integer k0 such that rank $\left( A _ { k } \right)$ is constant for all $k \geq k _ { 0 }$ .

P5.5.5 Show that if $A \in \mathbb { R } ^ { m \times n }$ has rank n, then so does $A + E { \mathrm { ~ i f ~ } } \| E \| _ { 2 } \| A ^ { + } \| _ { 2 } < 1$ .

P5.5.6 Suppose $A \in \mathbb { R } ^ { m \times n }$ is rank deficient and $b \in \mathbb { R } ^ { m }$ . Assume for $k = 0 , 1 , \ldots$ . that $x ^ { ( k + 1 ) }$ minimizes

$$
\phi_ {k} (x) = \parallel A x - b \parallel_ {2} ^ {2} + \lambda \parallel x - x ^ {(k)} \parallel_ {2} ^ {2}
$$

where $\lambda > 0$ and $x ^ { ( 0 ) } = 0$ . Show that $x ^ { ( k ) } \to x _ { L S }$

P5.5.8 Suppose $A \in \mathbb { R } ^ { m \times n }$ and that $\Vert u ^ { T } A \Vert _ { 2 } = \sigma$ with $u ^ { T } u = 1$ . Show that if $u ^ { T } ( A x - b ) = 0$ for $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $b \in \mathbb { R } ^ { m }$ , then $\parallel x \parallel _ { 2 } \geq | u ^ { T } b | / \sigma$ .

P5.5.9 In Equation (5.5.6) we know that the matrix $P ^ { T } V$ is orthogonal. Thus, $\parallel \tilde { V } _ { 1 1 } ^ { - 1 } \parallel _ { 2 } = \parallel \tilde { V } _ { 2 2 } ^ { - 1 } \parallel _ { 2 }$ from the CS decomposition (Theorem 2.5.3). Show how to compute P by applying the QR with column-pivoting algorithm to $[ \tilde { V } _ { 2 2 } ^ { T } | \tilde { V } _ { 1 2 } ^ { T } ]$ . (For $\tilde { r } > n / 2$ , this procedure would be more economical than the technique discussed in the text.) Incorporate this observation in Algorithm 5.5.1.

P5.5.10 Suppose $F \in \mathbb { R } ^ { m \times r }$ and $G \in \mathbb { R } ^ { n \times r }$ each have rank r. (a) Give an efficient algorithm for computing the minimum 2-norm minimizer of $\parallel F G ^ { T } x - b \parallel _ { 2 }$ where $b \in \mathbb { R } ^ { m }$ . (b) Show how to compute the vector $x _ { B }$ .

# Notes and References for §5.5

For a comprehensive treatment of the pseudoinverse and its manipulation, see:

M.Z. Nashed (1976). Generalized Inverses and Applications, Academic Press, New York.

S.L. Campbell and C.D. Meyer (2009). Generalized Inverses of Linear Transformations, SIAM Publications, Philadelphia, PA.

For an analysis of how the pseudo-inverse is affected by perturbation, see:

P.A. Wedin (1973). “Perturbation Theory for Pseudo-Inverses,” BIT 13, 217–232.

G.W. Stewart (1977). “On the Perturbation of Pseudo-Inverses, Projections, and Linear Least Squares,” SIAM Review 19, 634–662.

Even for full rank problems, column pivoting seems to produce more accurate solutions. The error analysis in the following paper attempts to explain why:

L.S. Jennings and M.R. Osborne (1974). “A Direct Error Analysis for Least Squares,” Numer. Math. 22, 322–332.

Various other aspects of the rank-deficient least squares problem are discussed in:

J.M. Varah (1973). “On the Numerical Solution of Ill-Conditioned Linear Systems with Applications to Ill-Posed Problems,” SIAM J. Numer. Anal. 10, 257–67.

G.W. Stewart (1984). “Rank Degeneracy,” SIAM J. Sci. Stat. Comput. 5, 403–413.

P.C. Hansen (1987). “The Truncated SVD as a Method for Regularization,” BIT 27, 534–553.

G.W. Stewart (1987). “Collinearity and Least Squares Regression,” Stat. Sci. 2, 68–100.

R.D. Fierro and P.C. Hansen (1995). “Accuracy of TSVD Solutions Computed from Rank-Revealing Decompositions,” Numer. Math. 70, 453–472.   
P.C. Hansen (1997). Rank-Deficient and Discrete Ill-Posed Problems: Numerical Aspects of Linear Inversion, SIAM Publications, Philadelphia, PA.   
A. Dax and L. Elden (1998). “Approximating Minimum Norm Solutions of Rank-Deficient Least Squares Problems,” Numer. Lin. Alg. 5, 79–99.   
G. Quintana-Orti, E.S. Quintana-Orti, and A. Petitet (1998). “Efficient Solution of the Rank-Deficient Linear Least Squares Problem,” SIAM J. Sci. Comput. 20, 1155–1163.   
L.V. Foster (2003). “Solving Rank-Deficient and Ill-posed Problems Using UTV and QR Factorizations,” SIAM J. Matrix Anal. Applic. 25, 582–600.   
D.A. Huckaby and T.F. Chan (2004). “Stewart’s Pivoted QLP Decomposition for Low-Rank Matrices,” Numer. Lin. Alg. 12, 153–159.   
L. Foster and R. Kommu (2006). “Algorithm 853: An Efficient Algorithm for Solving Rank-Deficient Least Squares Problems,” ACM Trans. Math. Softw. 32, 157–165.

For a sampling of the subset selection literature, we refer the reader to:

H. Hotelling (1957). “The Relations of the Newer Multivariate Statistical Methods to Factor Analysis,” Brit. J. Stat. Psych. 10, 69–79.

G.H. Golub, V. Klema and G.W. Stewart (1976). “Rank Degeneracy and Least Squares Problems,” Technical Report TR-456, Department of Computer Science, University of Maryland, College Park, MD.

S. Van Huffel and J. Vandewalle (1987). “Subset Selection Using the Total Least Squares Approach in Collinearity Problems with Errors in the Variables,” Lin. Alg. Applic. 88/89, 695–714.

M.R. Osborne, B. Presnell, and B.A. Turlach (2000). “A New Approach to Variable Selection in Least Squares Problems,” IMA J. Numer. Anal. 20, 389–403.

# 5.6 Square and Underdetermined Systems

The orthogonalization methods developed in this chapter can be applied to square systems and also to systems in which there are fewer equations than unknowns. In this brief section we examine the various possibilities.

# 5.6.1 Square Systems

The least squares solvers based on the QR factorization and the SVD can also be used to solve square linear systems. Figure 5.6.1 compares the associated flop counts. It is

<table><tr><td>Method</td><td>Flops</td></tr><tr><td>Gaussian elimination</td><td> $2n^{3}/3$ </td></tr><tr><td>Householder QR</td><td> $4n^{3}/3$ </td></tr><tr><td>Modified Gram-Schmidt</td><td> $2n^{3}$ </td></tr><tr><td>Singular value decomposition</td><td> $12n^{3}$ </td></tr></table>

Figure 5.6.1. Flops associated with various methods for square linear systems

assumed that the right-hand side is available at the time of factorization. Although Gaussian elimination involves the least amount of arithmetic, there are three reasons why an orthogonalization method might be considered:

• The flop counts tend to exaggerate the Gaussian elimination advantage. When memory traffic and vectorization overheads are considered, the QR approach is comparable in efficiency.

• The orthogonalization methods have guaranteed stability; there is no “growth factor” to worry about as in Gaussian elimination.

• In cases of ill-conditioning, the orthogonal methods give an added measure of reliability. QR with condition estimation is very dependable and, of course, SVD is unsurpassed when it comes to producing a meaningful solution to a nearly singular system.

We are not expressing a strong preference for orthogonalization methods but merely suggesting viable alternatives to Gaussian elimination.

We also mention that the SVD entry in the above table assumes the availability of b at the time of decomposition. Otherwise, $2 0 n ^ { 3 }$ flops are required because it then becomes necessary to accumulate the U matrix.

If the QR factorization is used to solve $A x = b$ , then we ordinarily have to carry out a back substitution: $R x = Q ^ { T } b$ . However, this can be avoided by “preprocessing” b. Suppose H is a Householder matrix such that $H b = \beta e _ { n }$ where $e _ { n }$ is the last column of $I _ { n }$ . If we compute the QR factorization of $( H A ) ^ { T }$ , then $A = H ^ { T } R ^ { T } Q ^ { T }$ and the system transforms to

$$
R ^ {T} y = \beta e _ {n}
$$

where $y = Q ^ { T } x$ . Since $R ^ { T }$ is lower triangular, $y = ( \beta / r _ { n n } ) e _ { n }$ and so

$$
x = \frac {\beta}{r _ {n n}} Q (:, n).
$$

# 5.6.2 Underdetermined Systems

In §3.4.8 we discussed how Gaussian elimination with either complete pivoting or rook pivoting can be used to solve a full-rank, underdetermined linear system

$$
A x = b, \quad A \in \mathbb {R} ^ {m \times n}, b \in \mathbb {R} ^ {m}. \tag {5.6.1}
$$

Various orthogonal factorizations can also be used to solve this problem. Notice that (5.6.1) either has no solution or has an infinity of solutions. In the second case, it is important to distinguish between algorithms that find the minimum 2-norm solution and those that do not. The first algorithm we present is in the latter category.

Assume that A has full row rank and that we apply QR with column pivoting to obtain

$$
Q ^ {T} A \Pi = \left[ R _ {1} \mid R _ {2} \right]
$$

where $R _ { 1 } \in \mathbb { R } ^ { m \times m }$ is upper triangular and $R _ { 2 } \in \mathbb { R } ^ { m \times ( n - m ) }$ . Thus, $A x = b$ transforms to

$$
(Q ^ {T} A \Pi) (\Pi^ {T} x) = [ R _ {1} \mid R _ {2} ] \left[ \begin{array}{c} z _ {1} \\ z _ {2} \end{array} \right] = Q ^ {T} b
$$

where

$$
\Pi^ {T} x = \left[ \begin{array}{l} z _ {1} \\ z _ {2} \end{array} \right]
$$

with $z _ { 1 } \in \mathbb { R } ^ { m }$ and $z _ { 2 } \in \mathbb { R } ^ { ( n - m ) }$ . By virtue of the column pivoting, $R _ { 1 }$ is nonsingular because we are assuming that A has full row rank. One solution to the problem is therefore obtained by setting $z _ { 1 } = R _ { 1 } ^ { - 1 } Q ^ { T } b$ and $z _ { 2 } = 0$ .

Algorithm 5.6.1 Given $A \in \mathbb { R } ^ { m \times n }$ with rank(A) = m and $b \in \mathbb { R } ^ { m }$ , the following algorithm finds an $\boldsymbol { x } \in \mathbb { R } ^ { n }$ such that $A x = b$ .

Compute QR-with-column-pivoting factorization: $Q ^ { T } A \Pi = R$

Solve $R ( 1 { : } m , 1 { : } m ) z _ { 1 } = Q ^ { T } b .$

$x = \Pi \left[ { \begin{array} { c } { z _ { 1 } } \\ { 0 } \end{array} } \right]$

This algorithm requires $2 m ^ { 2 } n - m ^ { 3 } / 3$ flops. The minimum norm solution is not guaranteed. (A different Π could render a smaller $z _ { 1 } . \ i )$ However, if we compute the QR factorization

$$
A ^ {T} = Q R = Q \left[ \begin{array}{c} R _ {1} \\ 0 \end{array} \right]
$$

with $R _ { 1 } \in \mathbb { R } ^ { m \times m }$ , then $A x = b$ becomes

$$
(Q R) ^ {T} x   =   \left[ \begin{array}{c c} R _ {1} ^ {T} & 0 \end{array} \right] \left[ \begin{array}{c} z _ {1} \\ z _ {2} \end{array} \right]   =   b,
$$

where

$$
Q ^ {T} x = \left[ \begin{array}{l} z _ {1} \\ z _ {2} \end{array} \right], \qquad z _ {1} \in \mathbb {R} ^ {m},   z _ {2} \in \mathbb {R} ^ {n - m}.
$$

In this case the minimum norm solution does follow by setting $z _ { 2 } = 0$ .

Algorithm 5.6.2 Given $A \in \mathbb { R } ^ { m \times n }$ with rank $( A ) = m$ and $b \in \mathbb { R } ^ { m }$ , the following algorithm finds the minimum 2-norm solution to $A x = b$ .

Compute the QR factorization $A ^ { T } = Q R$ .

Solve $R ( 1 { : } m , 1 { : } m ) ^ { T } z = b$

Set $x = Q ( : , 1 { : } m ) z .$

This algorithm requires at most $2 m ^ { 2 } n - 2 m ^ { 3 } / 3$ flops.

The SVD can also be used to compute the minimum norm solution of an underdetermined $A x = b$ problem. If

$$
A = \sum_ {i = 1} ^ {r} \sigma_ {i} u _ {i} v _ {i} ^ {T}, \quad r = \operatorname{rank} (A)
$$

is the SVD of A, then

$$
x = \sum_ {i = 1} ^ {r} \frac {u _ {i} ^ {T} b}{\sigma_ {i}} v _ {i}.
$$

As in the least squares problem, the SVD approach is desirable if A is nearly rank deficient.

# 5.6.3 Perturbed Underdetermined Systems

We conclude this section with a perturbation result for full-rank underdetermined systems.

Theorem 5.6.1. Suppose rank $( A ) = m \leq n$ and that $A \in \mathbb { R } ^ { m \times n }$ , $\delta A \in \mathbb { R } ^ { m \times n }$ , $0 \neq$ $b \in \mathbb { R } ^ { m }$ , and $\delta b \in \mathbb { R } ^ { m }$ satisfy

$$
\epsilon = \max \{\epsilon_ {A}, \epsilon_ {b} \} <   \sigma_ {m} (A),
$$

where $\epsilon _ { A } = \parallel \delta A \parallel _ { 2 } / \parallel A \parallel _ { 2 }$ and $\epsilon _ { b } = \parallel \delta b \parallel _ { 2 } / \parallel b \parallel _ { 2 }$ . If x and $\hat { x }$ are minimum norm solutions that satisfy

$$
A x = b, \quad (A + \delta A) \hat {x} = b + \delta b,
$$

then

$$
\frac {\| \hat {x} - x \| _ {2}}{\| x \| _ {2}} \leq \kappa_ {2} (A) (\epsilon_ {A} \min \{2, n - m + 1 \} + \epsilon_ {b}) + O (\epsilon^ {2}).
$$

Proof. Let E and f be defined by $\delta A / \epsilon$ and $\delta b / \epsilon$ . Note that rank $( A + t E ) = m$ for all $0 < t < \epsilon$ and that

$$
x (t) = (A + t E) ^ {T} \left((A + t E) (A + t E) ^ {T}\right) ^ {- 1} (b + t f)
$$

satisfies $( A + t E ) x ( t ) = b + t f$ . By differentiating this expression with respect to t and setting t = 0 in the result we obtain

$$
\dot {x} (0) = \left(I - A ^ {T} (A A ^ {T}) ^ {- 1} A\right) E ^ {T} (A A ^ {T}) ^ {- 1} b + A ^ {T} (A A ^ {T}) ^ {- 1} (f - E x). \tag {5.6.2}
$$

Because

$$
\| x \| _ {2} = \| A ^ {T} (A A ^ {T}) ^ {- 1} b \| _ {2} \geq \sigma_ {m} (A) \| (A A ^ {T}) ^ {- 1} b \| _ {2},
$$

$$
\| I - A ^ {T} (A A ^ {T}) ^ {- 1} A \| _ {2} = \min (1, n - m),
$$

and

$$
\frac {\| f \| _ {2}}{\| x \| _ {2}} \leq \frac {\| f \| _ {2} \| A \| _ {2}}{\| b \| _ {2}},
$$

we have

$$
\begin{array}{l} \frac {\| \hat {x} - x \| _ {2}}{\| x \| _ {2}} = \frac {x (\epsilon) - x (0)}{\| x (0) \| _ {2}} = \epsilon \frac {\| \dot {x} (0) \| _ {2}}{\| x \| _ {2}} + O (\epsilon^ {2}) \\ \leq \epsilon \min (1, n - m) \left\{\frac {\| E \| _ {2}}{\| A \| _ {2}} + \frac {\| f \| _ {2}}{\| b \| _ {2}} + \frac {\| E \| _ {2}}{\| A \| _ {2}} \right\} \kappa_ {2} (A) + O (\epsilon^ {2}), \\ \end{array}
$$

from which the theorem follows.

Note that there is no $\kappa _ { 2 } { \left( A \right) } ^ { 2 }$ factor as in the case of overdetermined systems.

# Problems

P5.6.1 Derive equation (5.6.2).

P5.6.2 Find the minimal norm solution to the system Ax = b where $A = [ 1 2 3 ]$ and $b = 1$

P5.6.3 Show how triangular system solving can be avoided when using the QR factorization to solve an underdetermined system.

P5.6.4 Suppose $b , x \in \mathbb { R } ^ { n }$ are given and consider the following problems:

(a) Find an unsymmetric Toeplitz matrix T so $T x = b ,$ .   
(b) Find a symmetric Toeplitz matrix T so $T x = b ,$   
(c) Find a circulant matrix C so Cx = b.

Pose each problem in the form $A p = b$ where A is a matrix made up of entries from x and p is the vector of sought-after parameters.

# Notes and References for 5.6

For an analysis of linear equation solving via QR, see:

N.J. Higham (1991). “Iterative Refinement Enhances the Stability of QR Factorization Methods for Solving Linear Equations,” BIT 31, 447–468.

Interesting aspects concerning singular systems are discussed in:

T.F. Chan (1984). “Deflated Decomposition Solutions of Nearly Singular Systems,” SIAM J. Numer. Anal. 21, 738–754.

Papers concerned with underdetermined systems include:

R.E. Cline and R.J. Plemmons (1976). “L2-Solutions to Underdetermined Linear Systems,” SIAM Review 18, 92–106.

M.G. Cox (1981). “The Least Squares Solution of Overdetermined Linear Equations having Band or Augmented Band Structure,” IMA J. Numer. Anal. 1, 3–22.

M. Arioli and A. Laratta (1985). “Error Analysis of an Algorithm for Solving an Underdetermined System,” Numer. Math. 46, 255–268.

J.W. Demmel and N.J. Higham (1993). “Improved Error Bounds for Underdetermined System Solvers,” SIAM J. Matrix Anal. Applic. 14, 1–14.

S. Jokar and M.E. Pfetsch (2008). “Exact and Approximate Sparse Solutions of Underdetermined Linear Equations,” SIAM J. Sci. Comput. 31, 23–44.

The central matrix problem in the emerging field of compressed sensing is to solve an underdetermined system Ax = b such that the 1-norm of x is minimized, see:

E. Candes, J. Romberg, and T. Tao (2006). “Robust Uncertainty Principles: Exact Signal Reconstruction from Highly Incomplete Frequency Information,” IEEE Trans. Information Theory 52, 489–509.

D. Donoho (2006). “Compressed Sensing,” IEEE Trans. Information Theory 52, 1289–1306.

This strategy tends to produce a highly sparse solution vector x.
