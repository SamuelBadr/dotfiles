# Chapter 9

# Functions of Matrices

9.1 Eigenvalue Methods   
9.2 Approximation Methods   
9.3 The Matrix Exponential   
9.4 The Sign, Square Root, and Log of a Matrix

Computing a function $f ( A )$ of an n-by-n matrix A is a common problem in many application areas. Roughly speaking, if the scalar function $f ( z )$ is defined on $\lambda ( A )$ , then $f ( A )$ is defined by substituting $^ { 6 6 } A ^ { \dag }$ for $^ { 6 6 } z ^ { 5 9 }$ in the “formula” for $f ( z )$ . For example, if $f ( z ) = ( 1 + z ) / ( 1 - z )$ and $1 \not \in \lambda ( A )$ , then $f ( A ) = ( I + A ) ( I - A ) ^ { - 1 }$ .

The computations get particularly interesting when the function f is transcendental. One approach in this more complicated situation is to compute an eigenvalue decomposition $A = Y B Y ^ { - 1 }$ and use the formula $f ( A ) = Y f ( B ) Y ^ { - 1 }$ . If B is sufficiently simple, then it is often possible to calculate f (B) directly. This is illustrated in §9.1 for the Jordan and Schur decompositions.

Another class of methods involves the approximation of the desired function $f ( A )$ with an easy-to-calculate function $g ( A )$ . For example, g might be a truncated Taylor series approximation to f. Error bounds associated with the approximation of matrix functions are given in §9.2.

In §9.3 we discuss the special and very important problem of computing the matrix exponential $e ^ { A }$ . The matrix sign, square root, and logarithm functions and connections to the polar decomposition are treated in §9.4.

# Reading Notes

Knowledge of Chapters 3 and 7 is assumed. Within this chapter there are the following dependencies:

$$
\begin{array}{c c c c c} \S 9. 1 & \to & \S 9. 2 & \to & \S 9. 3 \\ & & \downarrow & & \\ & & \S 9. 4 & & \end{array}
$$

Complementary references include Horn and Johnson (TMA) and the definitive text by Higham (FOM). We mention that aspects of the f(A)-times-a-vector problem are treated in §10.2.
