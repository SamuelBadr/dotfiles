# Chapter 2

# Matrix Analysis

2.1 Basic Ideas from Linear Algebra   
2.2 Vector Norms   
2.3 Matrix Norms   
2.4 The Singular Value Decomposition   
2.5 Subspace Metrics   
2.6 The Sensitivity of Square Systems   
2.7 Finite Precision Matrix Computations

The analysis and derivation of algorithms in the matrix computation area requires a facility with linear algebra. Some of the basics are reviewed in §2.1. Norms are particularly important, and we step through the vector and matrix cases in §2.2 and §2.3. The ubiquitous singular value decomposition is introduced in §2.4 and then used in the next section to define the CS decomposition and its ramifications for the measurement of subspace separation. In §2.6 we examine how the solution to a linear system Ax = b changes if A and b are perturbed. It is the ideal setting for introducing the concepts of problem sensitivity, backward error analysis, and condition number. These ideas are central throughout the text. To complete the chapter we develop a model of finite-precision floating point arithmetic based on the IEEE standard. Several canonical examples of roundoff error analysis are offered.

# Reading Notes

Familiarity with matrix manipulation consistent with §1.1–§1.3 is essential. The sections within this chapter depend upon each other as follows:

$$
\begin{array}c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
$$

Complementary references include Forsythe and Moler (SLAS), Stewart (IMC), Horn and Johnson (MA), Stewart (MABD), Ipsen (NMA), and Watkins (FMC). Fundamentals of matrix analysis that are specific to least squares problems and eigenvalue problems appear in later chapters.

# 2.1 Basic Ideas from Linear Algebra

This section is a quick review of linear algebra. Readers who wish a more detailed coverage should consult the references at the end of the section.

# 2.1.1 Independence, Subspace, Basis, and Dimension

A set of vectors $\{ a _ { 1 } , \ldots , a _ { n } \}$ in $\mathbb { R } ^ { m }$ is linearly independent if $\begin{array} { r } { \sum _ { j = 1 } ^ { n } \alpha _ { j } a _ { j } = 0 } \end{array}$ implies $\alpha ( 1 { : } n ) = 0$ . Otherwise, a nontrivial combination of the $a _ { i }$ is zero and $\{ a _ { 1 } , \ldots , a _ { n } \}$ is said to be linearly dependent.

A subspace of $\mathbb { R } ^ { m }$ is a subset that is also a vector space. Given a collection of vectors $\boldsymbol { a } _ { 1 } , \dots , \boldsymbol { a } _ { n } \in \mathbb { R } ^ { m }$ , the set of all linear combinations of these vectors is a subspace referred to as the span of $\{ a _ { 1 } , \ldots , a _ { n } \}$ :

$$
\operatorname{span} \left\{a _ {1}, \dots , a _ {n} \right\} = \left\{\sum_ {j = 1} ^ {n} \beta_ {j} a _ {j}: \beta_ {j} \in \mathbb {R} \right\}.
$$

If $\{ a _ { 1 } , \ldots , a _ { n } \}$ is independent and $b \in { \mathsf { s p a n } } \{ a _ { 1 } , \ldots , a _ { n } \}$ , then b is a unique linear combination of the aj. $a _ { j }$

If $S _ { 1 } , \ldots , S _ { k }$ are subspaces of $\mathbb { R } ^ { m }$ , then their sum is the subspace defined by $S = \{ \ a _ { 1 } + a _ { 2 } + \cdots + a _ { k } : a _ { i } \in S _ { i } , \ i = 1 { : } k \ \}$ . S is said to be a direct sum if each $v \in S$ has a unique representation $v = a _ { 1 } + \cdots + a _ { k }$ with $a _ { i } \in S _ { i }$ . In this case we write $S = S _ { 1 } \oplus \cdot \cdot \cdot \oplus S _ { k }$ . The intersection of the $S _ { i }$ is also a subspace, $S = S _ { 1 } \cap S _ { 2 } \cap \cdot \cdot \cdot \cap S _ { k }$ .

The subset $\{ a _ { i _ { 1 } } , \ldots , a _ { i _ { k } } \}$ is a maximal linearly independent subset of $\{ a _ { 1 } , \ldots , a _ { n } \}$ if it is linearly independent and is not properly contained in any linearly independent subset of $\{ a _ { 1 } , \ldots , a _ { n } \}$ . If $\{ a _ { i _ { 1 } } , \ldots , a _ { i _ { k } } \}$ is maximal, then span $\{ a _ { 1 } , \ldots , a _ { n } \} \ =$ span $\{ a _ { i _ { 1 } } , \ldots , a _ { i _ { k } } \}$ and $\{ a _ { i _ { 1 } } , \ldots , a _ { i _ { k } } \}$ is a basis for span $\{ a _ { 1 } , \ldots , a _ { n } \}$ . If $S \subseteq \mathbb { R } ^ { m }$ is a subspace, then it is possible to find independent basic vectors $a _ { 1 } , \dots , a _ { k } \in S$ such that $S = \mathfrak { s p a n } \{ a _ { 1 } , \dotsc , a _ { k } \}$ . All bases for a subspace $S$ have the same number of elements. This number is the dimension and is denoted by dim(S).

# 2.1.2 Range, Null Space, and Rank

There are two important subspaces associated with an m-by-n matrix A. The range of A is defined by

$$
\operatorname{ran} (A) = \{y \in \mathbb {R} ^ {m}: y = A x \text {   for   some   } x \in \mathbb {R} ^ {n} \}
$$

and the nullspace of A is defined by

$$
\operatorname{null} (A) = \{x \in \mathbb {R} ^ {n}: A x = 0 \}.
$$

If $A = \left[ { a _ { 1 } } \mid \cdots \mid { a _ { n } } \right]$ is a column partitioning, then

$$
\operatorname{ran} (A) = \operatorname{span} \left\{a _ {1}, \dots , a _ {n} \right\}.
$$

The rank of a matrix A is defined by

$$
\operatorname{rank} (A) = \dim (\operatorname{ran} (A)).
$$

If $A \in \mathbb { R } ^ { m \times n }$ , then

$$
\dim (\operatorname{null} (A)) + \operatorname{rank} (A) = n.
$$

We say that $A \in \mathbb { R } ^ { m \times n }$ is rank deficient if rank $( A ) < \operatorname* { m i n } \{ m , n \}$ . The rank of a matrix is the maximal number of linearly independent columns (or rows).

# 2.1.3 Matrix Inverse

If A and X are in $\mathbb { R } ^ { n \times n }$ and satisfy $A X = I$ , then X is the inverse of A and is denoted by $A ^ { - 1 }$ . If $A ^ { - 1 }$ exists, then A is said to be nonsingular. Otherwise, we say A is singular. The inverse of a product is the reverse product of the inverses:

$$
(A B) ^ {- 1} = B ^ {- 1} A ^ {- 1}. \tag {2.1.1}
$$

Likewise, the transpose of the inverse is the inverse of the transpose:

$$
(A ^ {- 1}) ^ {T} = (A ^ {T}) ^ {- 1} \equiv A ^ {- T}. \tag {2.1.2}
$$

# 2.1.4 The Sherman-Morrison-Woodbury Formula

The identity

$$
B ^ {- 1} = A ^ {- 1} - B ^ {- 1} (B - A) A ^ {- 1} \tag {2.1.3}
$$

shows how the inverse changes if the matrix changes. The Sherman-Morrison-Woodbury formula gives a convenient expression for the inverse of the matrix $( A + U V ^ { T } )$ where $A \in \mathbb { R } ^ { n \times n }$ and U and V are n-by-k:

$$
(A + U V ^ {T}) ^ {- 1} = A ^ {- 1} - A ^ {- 1} U (I + V ^ {T} A ^ {- 1} U) ^ {- 1} V ^ {T} A ^ {- 1}. \tag {2.1.4}
$$

A rank-k correction to a matrix results in a rank-k correction of the inverse. In (2.1.4) we assume that both A and $( I + V ^ { T } A ^ { - 1 } U )$ are nonsingular.

The $k = 1$ case is particularly useful. If $A \in \mathbb { R } ^ { n \times n }$ is nonsingular, $u , v \in \mathbb { R } ^ { n }$ , and $\alpha = 1 + v ^ { T } A ^ { - 1 } u \neq 0$ , then

$$
(A + u v ^ {T}) ^ {- 1} = A ^ {- 1} - \frac {1}{\alpha} A ^ {- 1} u v ^ {T} A ^ {- 1}. \tag {2.1.5}
$$

This is referred to as the Sherman-Morrison formula.

# 2.1.5 Orthogonality

A set of vectors $\{ x _ { 1 } , \ldots , x _ { p } \}$ in $\mathbb { R } ^ { m }$ is orthogonal if $x _ { i } ^ { T } x _ { j } = 0$ whenever $i \neq j$ and orthonormal if $x _ { i } ^ { T } x _ { j } = \delta _ { i j }$ . Intuitively, orthogonal vectors are maximally independent for they point in totally different directions.

A collection of subspaces $S _ { 1 } , \ldots , S _ { p }$ in $\mathbb { R } ^ { m }$ is mutually orthogonal if $x ^ { T } y = 0$ whenever $x \in S _ { i }$ and $y \in S _ { j }$ for $i \neq j$ . The orthogonal complement of a subspace $S \subseteq \mathbb { R } ^ { m }$ is defined by

$$
S ^ {\perp} = \{y \in \mathbb {R} ^ {m}: y ^ {T} x = 0 \text { for   all } x \in S \}.
$$

It is not hard to show that ran $( A ) ^ { \perp } = \mathsf { n u l l } ( A ^ { T } )$ . The vectors $v _ { 1 } , \ldots , v _ { k }$ form an $o r \mathrm { - }$ thonormal basis for a subspace $S \subseteq \mathbb { R } ^ { m }$ if they are orthonormal and span S.

A matrix $Q \in \mathbb { R } ^ { m \times m }$ is said to be orthogonal if $Q ^ { T } Q = I$ . If $Q = [  q _ { 1 } | \cdot \cdot \cdot | q _ { m } ]$ is orthogonal, then the $q _ { i }$ form an orthonormal basis for $\mathbb { R } ^ { m }$ . It is always possible to extend such a basis to a full orthonormal basis $\{ v _ { 1 } , \ldots , v _ { m } \}$ for $\mathbb { R } ^ { m }$ :

Theorem 2.1.1. $I f V _ { 1 } \in \mathbb { R } ^ { n \times r }$ has orthonormal columns, then there exists $V _ { 2 } \in \mathbb { R } ^ { n \times ( n - r ) }$ such that

$$
V = \left[ V _ {1} \mid V _ {2} \right]
$$

is orthogonal. Note that ran $( V _ { 1 } ) ^ { \perp } = \mathsf { r a n } ( V _ { 2 } )$

Proof. This is a standard result from introductory linear algebra. It is also a corollary of the QR factorization that we present in §5.2.

# 2.1.6 The Determinant

If $A = ( a ) \in \mathbb { R } ^ { 1 \times 1 }$ , then its determinant is given by ${ \mathsf { d e t } } ( A ) = a$ . The determinant of $A \in \mathbb { R } ^ { n \times n }$ is defined in terms of order-(n−1) determinants:

$$
\det (A) = \sum_ {j = 1} ^ {n} (- 1) ^ {j + 1} a _ {1 j} \det (A _ {1 j}).
$$

Here, $A _ { 1 j }$ is an $( n - 1 ) – \mathrm { b y } – ( n - 1 )$ matrix obtained by deleting the first row and jth column of A. Well-known properties of the determinant include det $( A B ) = \mathsf { d e t } ( A ) \mathsf { d e t } ( B )$ , de $\langle A ^ { T } \rangle = \operatorname* { d e t } ( A )$ , and det $( c A ) = c ^ { n } \mathsf { d e t } ( A )$ where A, $B \in \mathbb { R } ^ { n \times n }$ and $c \in \mathbb { R }$ . In addition, det $( A ) \neq 0$ if and only if A is nonsingular.

# 2.1.7 Eigenvalues and Eigenvectors

Until we get to the main eigenvalue part of the book (Chapters 7 and 8), we need a handful of basic properties so that we can fully appreciate the singular value decomposition (§2.4), positive definiteness (§4.2), and various fast linear equation solvers (§4.8).

The eigenvalues of $A \in \mathbb { C } ^ { n \times n }$ are the zeros of the characteristic polynomial

$$
p (x) = \det (A - x I).
$$

Thus, every n-by-n matrix has n eigenvalues. We denote the set of A’s eigenvalues by

$$
\lambda (A) = \{x: \det (A - x I) = 0 \}.
$$

If the eigenvalues of A are real, then we index them from largest to smallest as follows:

$$
\lambda_ {n} (A) \leq \dots \leq \lambda_ {2} (A) \leq \lambda_ {1} (A).
$$

In this case, we sometimes use the notation $\lambda _ { \mathrm { m a x } } ( A )$ and $\lambda _ { \mathrm { m i n } } ( A )$ to denote $\lambda _ { 1 } ( A )$ and $\lambda _ { n } ( A )$ respectively.

If $X \in \mathbb { C } ^ { n \times n }$ is nonsingular and $B = X ^ { - 1 } A X$ , then A and B are similar. If two matrices are similar, then they have exactly the same eigenvalues.

If $\lambda \in \lambda ( A )$ , then there exists a nonzero vector x so that $A x = \lambda x$ . Such a vector is said to be an eigenvector for A associated with λ. If $A \in \mathbb { C } ^ { n \times n }$ has n independent eigenvectors $x _ { 1 } , \ldots , x _ { n }$ and $A x _ { i } = \lambda _ { i } x _ { i }$ for $i = 1 { : } n$ , then A is diagonalizable. The terminology is appropriate for if

$$
X = \left[ x _ {1} \mid \dots \mid x _ {n} \right],
$$

then

$$
X ^ {- 1} A X = \operatorname{diag} \left(\lambda_ {1}, \dots , \lambda_ {n}\right).
$$

Not all matrices are diagonalizable. However, if $A \in \mathbb { R } ^ { n \times n }$ is symmetric, then there exists an orthogonal $Q$ so that

$$
Q ^ {T} A Q = \mathrm{diag} (\lambda_ {1}, \dots , \lambda_ {n}). \tag {2.1.6}
$$

This is called the Schur decomposition. The largest and smallest eigenvalues of a symmetric matrix satisfy

$$
\lambda_ {\max} (A) = \max _ {x \neq 0} \frac {x ^ {T} A x}{x ^ {T} x} \tag {2.1.7}
$$

and

$$
\lambda_ {\min} (A) = \min _ {x \neq 0} \frac {x ^ {T} A x}{x ^ {T} x}. \tag {2.1.8}
$$

# 2.1.8 Differentiation

Suppose α is a scalar and that $A ( \alpha )$ is an m-by-n matrix with entries $a _ { i j } ( \alpha )$ . If $a _ { i j } ( \alpha )$ is a differentiable function of α for all i and $j ,$ then by $\dot { A } ( \alpha )$ we mean the matrix

$$
\dot {A} (\alpha) = \frac {d}{d \alpha} A (\alpha) = \left(\frac {d}{d \alpha} a _ {i j} (\alpha)\right) = (\dot {a} _ {i j} (\alpha)).
$$

Differentiation is a useful tool that can sometimes provide insight into the sensitivity of a matrix problem.

# Problems

P2.1.1 Show that if $A \in \mathbb { R } ^ { m \times n }$ has rank p, then there exists an $\boldsymbol { X } \in \mathbb { R } ^ { m \times p }$ and ${ \mathrm { ~ a ~ } } Y \in \mathbb { R } ^ { n \times p }$ such that $A = X Y ^ { T }$ , where rank $\langle X \rangle = \operatorname { r a n k } ( Y ) = p ,$ .

P2.1.2 Suppose $A ( \alpha ) \in \mathbb { R } ^ { m \times r }$ and $B ( \alpha ) \in \mathbb { R } ^ { r \times n }$ are matrices whose entries are differentiable functions of the scalar α. (a) Show

$$
\frac {d}{d \alpha} [ A (\alpha) B (\alpha) ] = \left[ \frac {d}{d \alpha} A (\alpha) \right] B (\alpha) + A (\alpha) \left[ \frac {d}{d \alpha} B (\alpha) \right].
$$

(b) Assuming A(α) is always nonsingular, show

$$
\frac {d}{d \alpha} \left[ A (\alpha) ^ {- 1} \right] = - A (\alpha) ^ {- 1} \left[ \frac {d}{d \alpha} A (\alpha) \right] A (\alpha) ^ {- 1}.
$$

P2.1.3 Suppose $A \in \mathbb { R } ^ { n \times n } , \ b \in \mathbb { R } ^ { n }$ and that $\begin{array} { r } { \phi ( x ) = \frac { 1 } { 2 } x ^ { T } A x - x ^ { T } b } \end{array}$ . Show that the gradient of $\phi$ is given by $\begin{array} { r } { \nabla \phi ( x ) = \frac { 1 } { 2 } ( A ^ { T } + A ) x - b . } \end{array}$ .

P2.1.4 Assume that both A and $A + u v ^ { T }$ are nonsingular where $A \in \mathbb { R } ^ { n \times n }$ and $u , v \in \mathbb { R } ^ { n }$ . Show that if x solves $( A + u v ^ { T } ) x = b ,$ then it also solves a perturbed right-hand-side problem of the form $A x = b + \alpha u$ . Give an expression for α in terms of $A , u ,$ and v.

P2.1.5 Show that a triangular orthogonal matrix is diagonal.

P2.1.6 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and nonsingular and define

$$
\tilde {A} = A + \alpha (u u ^ {T} + v v ^ {T}) + \beta (u v ^ {T} + v u ^ {T})
$$

where u, $v \in \mathbb { R } ^ { n }$ and $\alpha , \beta \in \mathbb { R }$ . Assuming that $\tilde { A }$ is nonsingular, use the Sherman-Morrison-Woodbury formula to develop a formula for $\tilde { A } ^ { - 1 }$ .

P2.1.7 Develop a symmetric version of the Sherman-Morrison-Woodbury formula that characterizes the inverse of $A + U S U ^ { T }$ where $A \in \mathbb { R } ^ { n \times n }$ and $S \in \mathbb { R } ^ { k \times k }$ are symmetric and $U \in \mathbb { R } ^ { n \times k }$ .

P2.1.8 Suppose $Q \in \mathbb { R } ^ { n \times n }$ is orthogonal and $z \in \mathbb { R } ^ { n }$ . Give an efficient algorithm for setting up an m-by-m matrix $A = \left( a _ { i j } \right)$ defined by $a _ { i j } = v ^ { T } ( Q ^ { i } ) ^ { T } ( Q ^ { j } ) v .$ -

P2.1.9 Show that if S is real and $S ^ { T } = - S$ , then I − S is nonsingular and the matrix $( I - S ) ^ { - 1 } ( I + S )$ is orthogonal. This is known as the Cayley transform of S.

P2.1.10 Refer to §1.3.10. (a) Show that if $S \in \mathbb { R } ^ { 2 n \times 2 n }$ is symplectic, then $S ^ { - 1 }$ exists and is also symplectic. (b) Show that if $M \in \mathbb { R } ^ { 2 n \times 2 n }$ is Hamiltonian and $S \in \mathbb { R } ^ { 2 n \times 2 n }$ is symplectic, then the matrix $M _ { 1 } = S ^ { - 1 } M S$ is Hamiltonian.

P2.1.11 Use (2.1.6) to prove (2.1.7) and (2.1.8).

# Notes and References for 2.1

In addition to Horn and Johnson (MA) and Horn and Johnson (TMA), the following introductory applied linear algebra texts are highly recommended:

R. Bellman (1997). Introduction to Matrix Analysis, Second Edition, SIAM Publications, Philadelphia, PA.

C. Meyer (2000). Matrix Analysis and Applied Linear Algebra, SIAM Publications, Philadelphia, PA.

D. Lay (2005). Linear Algebra and Its Applications, Third Edition, Addison-Wesley, Reading, MA.

S.J. Leon (2007). Linear Algebra with Applications, Seventh Edition, Prentice-Hall, Englewood Cliffs, NJ.

G. Strang (2009). Introduction to Linear Algebra, Fourth Edition, SIAM Publications, Philadelphia, PA.

# 2.2 Vector Norms

A norm on a vector space plays the same role as absolute value: it furnishes a distance measure. More precisely, $\mathbb { R } ^ { n }$ together with a norm on $\mathbb { R } ^ { n }$ defines a metric space rendering the familiar notions of neighborhood, open sets, convergence, and continuity.

# 2.2.1 Definitions

A vector norm on $\mathbb { R } ^ { n }$ is a function $f { : \mathbb { R } ^ { n } } \to \mathbb { R }$ that satisfies the following properties:

$$
f (x) \geq 0, \quad x \in \mathbb {R} ^ {n}, \quad (f (x) = 0, \text {   iff   } x = 0),
$$

$$
f (x + y) \leq f (x) + f (y), \quad x, y \in \mathbb {R} ^ {n},
$$

$$
f (\alpha x) = | \alpha | f (x), \quad \alpha \in \mathbb {R}, x \in \mathbb {R} ^ {n}.
$$

We denote such a function with a double bar notation: $f ( x ) = \parallel x \parallel$ . Subscripts on the double bar are used to distinguish between various norms. A useful class of vector

norms are the p-norms defined by

$$
\| x \| _ {p} = \left(| x _ {1} | ^ {p} + \dots + | x _ {n} | ^ {p}\right) ^ {\frac {1}{p}}, \quad p \geq 1. \tag {2.2.1}
$$

The 1−, 2−, and $\infty -$ norms are the most important:

$$
\begin{array}{l} \| x \| _ {1} = | x _ {1} | + \dots + | x _ {n} |, \\ \| x \| _ {2} = \left(| x _ {1} | ^ {2} + \dots + | x _ {n} | ^ {2}\right) ^ {\frac {1}{2}} = \left(x ^ {T} x\right) ^ {\frac {1}{2}}, \\ \| x \| _ {\infty} = \max _ {1 \leq i \leq n} | x _ {i} |. \\ \end{array}
$$

A unit vector with respect to the norm $\parallel \cdot \parallel$ is a vector x that satisfies $\| { \boldsymbol { x } } \| = 1$ .

# 2.2.2 Some Vector Norm Properties

A classic result concerning p-norms is the H¨older inequality:

$$
\left| x ^ {T} y \right| \leq \| x \| _ {p} \| y \| _ {q} \quad \frac {1}{p} + \frac {1}{q} = 1. \tag {2.2.2}
$$

A very important special case of this is the Cauchy-Schwarz inequality:

$$
\left| x ^ {T} y \right| \leq \| x \| _ {2} \| y \| _ {2}. \tag {2.2.3}
$$

All norms on $\mathbb { R } ^ { n }$ are equivalent , i.e., $\operatorname { i f } \parallel \cdot \parallel _ { \alpha }$ and $\| \cdot \| _ { \beta }$ are norms on $\mathbb { R } ^ { n }$ , then there exist positive constants $c _ { 1 }$ and $c _ { 2 }$ such that

$$
c _ {1} \| x \| _ {\alpha} \leq \| x \| _ {\beta} \leq c _ {2} \| x \| _ {\alpha} \tag {2.2.4}
$$

for all $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . For example, if $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , then

$$
\| x \| _ {2} \leq \| x \| _ {1} \leq \sqrt {n} \| x \| _ {2}, \tag {2.2.5}
$$

$$
\| x \| _ {\infty} \leq \| x \| _ {2} \leq \sqrt {n} \| x \| _ {\infty}, \tag {2.2.6}
$$

$$
\| x \| _ {\infty} \leq \| x \| _ {1} \leq n \| x \| _ {\infty}. \tag {2.2.7}
$$

Finally, we mention that the 2-norm is preserved under orthogonal transformation. Indeed, if $Q \in \mathbb { R } ^ { n \times n }$ is orthogonal and $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , then

$$
\parallel Q x \parallel_ {2} ^ {2} = (Q x) ^ {T} (Q x) = (x ^ {T} Q ^ {T}) (Q x) = x ^ {T} (Q ^ {T} Q) x = x ^ {T} x = \parallel x \parallel_ {2} ^ {2}.
$$

# 2.2.3 Absolute and Relative Errors

Suppose $\hat { x } \in \mathbb { R } ^ { n }$ is an approximation to $\boldsymbol { x } \in \mathbb { R } ^ { n }$ . For a given vector norm $\| \cdot \|$ we say that

$$
\epsilon_ {\mathrm{abs}} = \left\| \hat {x} - x \right\|
$$

is the absolute error in ˆx. If $x \neq 0$ , then

$$
\epsilon_ {\mathrm{rel}} = \frac {\parallel \hat {x} - x \parallel}{\parallel x \parallel}
$$

prescribes the relative error in ˆx. Relative error in the ∞-norm can be translated into a statement about the number of correct significant digits in ${ \hat { x } } .$ In particular, if

$$
\frac {\parallel \hat {x} - x \parallel_ {\infty}}{\parallel x \parallel_ {\infty}} \approx 1 0 ^ {- p},
$$

then the largest component of ˆx has approximately p correct significant digits. For example, if $x = [ 1 . 2 3 4 ~ . 0 5 6 7 4 ] ^ { T }$ and $\hat { x } = [ 1 . 2 3 5 . 0 5 1 2 8 ] ^ { T }$ , then $\| \hat { x } - x \| _ { \infty } / \| x \| _ { \infty } \approx$ $. 0 0 4 3 \approx 1 0 ^ { - 3 }$ . Note than $\hat { x } _ { 1 }$ ∞ ∞ has about three significant digits that are correct while only one significant digit in ${ \hat { x } } _ { 2 }$ is correct.

# 2.2.4 Convergence

We say that a sequence $\{ x ^ { ( k ) } \}$ of n-vectors converges to x if

$$
\lim _ {k \to \infty} \| x ^ {(k)} - x \| = 0.
$$

Because of (2.2.4), convergence in any particular norm implies convergence in all norms.

# Problems

P2.2.1 Show that if $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , then li $1 _ { p \to \infty } \parallel x \parallel _ { p } = \parallel x \parallel _ { \infty }$ .

P2.2.2 By considering the inequality $0 \leq ( a x + b y ) ^ { T } ( a x + b y )$ for suitable scalars a and b, prove (2.2.3).

P2.2.3 Verify that $\| { \bf \cdot } \| _ { 1 } , \| { \bf \cdot } \| _ { 2 }$ , and $\| \cdot \| _ { \infty }$ are vector norms.

P2.2.4 Verify (2.2.5)-(2.2.7). When is equality achieved in each result?

P2.2.5 Show that in $\mathbb { R } ^ { n } , x ^ { ( i ) } \to x$ if and only if $x _ { k } ^ { ( i ) } \to x _ { k }$ for $k = 1 { : } n$

P2.2.6 Show that for any vector norm on $\mathbb { R } ^ { n }$ that $\mid \parallel x \parallel - \parallel y \parallel \mid \leq \parallel x - y \parallel ,$ .

P2.2.7 Let $\| \cdot \|$ be a vector norm on $\mathbb { R } ^ { m }$ and assume $A \in \mathbb { R } ^ { m \times n }$ . Show that if rank(A) = n, then $\| x \| _ { A } = \| A x \|$ is a vector norm on $\mathbb { R } ^ { n }$ .

P2.2.8 Let x and y be in $\mathbb { R } ^ { n }$ and define $\psi : \mathbb { R } \to \mathbb { R }$ by $\psi ( \alpha ) = \| { \boldsymbol { x } } - \alpha y \| _ { 2 }$ . Show that $\psi$ is minimized if $\alpha = x ^ { T } y / y ^ { T } y$ .

P2.2.9 Prove or disprove:

$$
v \in \mathbb {R} ^ {n} \Rightarrow \| v \| _ {1} \| v \| _ {\infty} \leq \frac {1 + \sqrt {n}}{2} \| v \| _ {2} ^ {2}.
$$

P2.2.10 If $\boldsymbol { x } \in \mathbb { R } ^ { 3 }$ and $\boldsymbol { y } \in \mathbb { R } ^ { 3 }$ then it can be shown that $| x ^ { T } y | = \parallel x \parallel _ { 2 } \parallel y \parallel _ { 2 } | \cos ( \theta ) |$ where $\theta$ is the angle between x and y. An analogous result exists for the cross product defined by

$$
x \times y = \left[ \begin{array}{l} x _ {2} y _ {3} - x _ {3} y _ {2} \\ x _ {3} y _ {1} - x _ {1} y _ {3} \\ x _ {1} y _ {2} - x _ {2} y _ {1} \end{array} \right].
$$

In particular, $\parallel x \times y \parallel _ { 2 } = \parallel x \parallel _ { 2 } \parallel y \parallel _ { 2 } | \sin ( \theta ) |$ . Prove this.

P2.2.11 Suppose $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $y \in \mathbb { R } ^ { m }$ . Show that

$$
\parallel x \otimes y \parallel_ {p} = \parallel x \parallel_ {p} \parallel y \parallel_ {p}
$$

for $p = 1 , 2$ , and ∞.

# Notes and References for 2.2

Although a vector norm is “just” a generalization of the absolute value concept, there are some noteworthy subtleties:

J.D. Pryce (1984). “A New Measure of Relative Error for Vectors,” SIAM J. Numer. Anal. 21, 202–221.

# 2.3 Matrix Norms

The analysis of matrix algorithms requires use of matrix norms. For example, the quality of a linear system solution may be poor if the matrix of coefficients is “nearly singular.” To quantify the notion of near-singularity, we need a measure of distance on the space of matrices. Matrix norms can be used to provide that measure.

# 2.3.1 Definitions

Since $\mathbb { R } ^ { m \times n }$ is isomorphic to $\mathbb { R } ^ { m n }$ , the definition of a matrix norm should be equivalent to the definition of a vector norm. In particular, $f { : \mathbb { R } ^ { m \times n } } \to \mathbb { R }$ is a matrix norm if the following three properties hold:

$$
f (A) \geq 0, \quad A \in \mathbb {R} ^ {m \times n}, \quad (f (A) = 0 \text {   iff   } A = 0)
$$

$$
f (A + B) \leq f (A) + f (B), \quad A, B \in \mathbb {R} ^ {m \times n},
$$

$$
f (\alpha A) = | \alpha | f (A), \quad \alpha \in \mathbb {R}, A \in \mathbb {R} ^ {m \times n}.
$$

As with vector norms, we use a double bar notation with subscripts to designate matrix norms, ${ \mathrm { i . e . , ~ } } \| \ A \| = f ( A )$ .

The most frequently used matrix norms in numerical linear algebra are the Frobenius norm

$$
\| A \| _ {F} = \sqrt {\sum_ {i = 1} ^ {m} \sum_ {j = 1} ^ {n} \left| a _ {i j} \right| ^ {2}} \tag {2.3.1}
$$

and the p-norms

$$
\| A \| _ {p} = \sup _ {x \neq 0} \frac {\| A x \| _ {p}}{\| x \| _ {p}}. \tag {2.3.2}
$$

Note that the matrix p-norms are defined in terms of the vector p-norms discussed in the previous section. The verification that (2.3.1) and (2.3.2) are matrix norms is left as an exercise. It is clear that $\| A \| _ { p }$ is the p-norm of the largest vector obtained by applying A to a unit p-norm vector:

$$
\| A \| _ {p} = \sup _ {x \neq 0} \left\| A \left(\frac {x}{\| x \| _ {p}}\right) \right\| _ {p} = \max _ {\| x \| _ {p} = 1} \| A x \| _ {p}.
$$

It is important to understand that (2.3.2) defines a family of norms—the 2-norm on $\mathbb { R } ^ { 3 \times 2 }$ is a different function from the 2-norm on $\mathbb { R } ^ { 5 \times 6 }$ . Thus, the easily verified inequality

$$
\| A B \| _ {p} \leq \| A \| _ {p} \| B \| _ {p}, \quad A \in \mathbb {R} ^ {m \times n}, B \in \mathbb {R} ^ {n \times q} \tag {2.3.3}
$$

is really an observation about the relationship between three different norms. Formally, we say that norms $f _ { 1 } , \ f _ { 2 }$ , and $f _ { 3 }$ on $\mathbb { R } ^ { m \times q }$ , $\mathbb { R } ^ { m \times n }$ , and $\mathbb { R } ^ { n \times q }$ are mutually consistent if for all matrices $A \in \mathbb { R } ^ { m \times n }$ and $B \in \mathbb { R } ^ { n \times q }$ we have $f _ { 1 } ( A B ) \le f _ { 2 } ( A ) f _ { 3 } ( B )$ , or, in subscript-free norm notation:

$$
\| A B \| \leq \| A \| \| B \|. \tag {2.3.4}
$$

Not all matrix norms satisfy this property. For example, if $\left\| \ A \right\| _ { \Delta } = \operatorname* { m a x } \left| a _ { i j } \right|$ and

$$
A = B = \left[ \begin{array}{l l} 1 & 1 \\ 1 & 1 \end{array} \right],
$$

then $\left\| \left. A B \right\| _ { \Delta } > \right\| \left. A \right\| _ { \Delta } \right\| B \left\| _ { \Delta }$ . For the most part, we work with norms that satisfy (2.3.4).

The p-norms have the important property that for every $A \in \mathbb { R } ^ { m \times n }$ and $\boldsymbol { x } \in \mathbb { R } ^ { n }$ we have

$$
\| A x \| _ {p} \leq \| A \| _ {p} \| x \| _ {p}.
$$

More generally, for any vector norm $\| \cdot \| _ { \alpha }$ on $\mathbb { R } ^ { n }$ and $\| \cdot \| _ { \beta }$ on $\mathbb { R } ^ { m }$ we have $\Vert { \mathbf { \Omega } } A x \Vert _ { \beta } \ \leq$ $\| A \| _ { \alpha , \beta } \| \boldsymbol { x } \| _ { \alpha }$ where $\| A \| _ { \alpha , \beta }$ is a matrix norm defined by

$$
\| A \| _ {\alpha , \beta} = \sup _ {x \neq 0} \frac {\| A x \| _ {\beta}}{\| x \| _ {\alpha}}. \tag {2.3.5}
$$

We say that $\| \cdot \| _ { \alpha , \beta }$ is subordinate to the vector norms $\| \cdot \| _ { \alpha }$ and $\| \cdot \| _ { \beta }$ . Since the set $\{ x \in \mathbb { R } ^ { n } : \| x \| _ { \alpha } ^ { - \gamma } = 1 \}$ is compact and $\| \cdot \| _ { \beta }$ is continuous, it follows that

$$
\| A \| _ {\alpha , \beta} = \max _ {\| x \| _ {\alpha} = 1} \| A x \| _ {\beta} = \| A x _ {*} \| _ {\beta} \tag {2.3.6}
$$

for some $\boldsymbol { x } _ { * } \in \mathbb { R } ^ { n }$ having unit α-norm.

# 2.3.2 Some Matrix Norm Properties

The Frobenius and p-norms (especially $p = 1 , 2 , \infty )$ satisfy certain inequalities that are frequently used in the analysis of a matrix computation. If $A \in \mathbb { R } ^ { m \times n }$ we have

$$
\| A \| _ {2} \leq \| A \| _ {F} \leq \sqrt {\min \{m , n \}} \| A \| _ {2}, \tag {2.3.7}
$$

$$
\max _ {i, j} \left| a _ {i j} \right| \leq \| A \| _ {2} \leq \sqrt {m n} \max _ {i, j} \left| a _ {i j} \right|, \tag {2.3.8}
$$

$$
\| A \| _ {1} = \max _ {1 \leq j \leq n} \sum_ {i = 1} ^ {m} | a _ {i j} |, \tag {2.3.9}
$$

$$
\| A \| _ {\infty} = \max _ {1 \leq i \leq m} \sum_ {j = 1} ^ {n} | a _ {i j} |, \tag {2.3.10}
$$

$$
\frac {1}{\sqrt {n}} \parallel A \parallel_ {\infty} \leq \parallel A \parallel_ {2} \leq \sqrt {m} \parallel A \parallel_ {\infty}, \tag {2.3.11}
$$

$$
\frac {1}{\sqrt {m}} \parallel A \parallel_ {1} \leq \parallel A \parallel_ {2} \leq \sqrt {n} \parallel A \parallel_ {1}. \tag {2.3.12}
$$

If $A \in \mathbb { R } ^ { m \times n } , 1 \leq i _ { 1 } \leq i _ { 2 } \leq m$ , and $1 \leq j _ { 1 } \leq j _ { 2 } \leq n$ , then

$$
\left\| A \left(i _ {1}: i _ {2}, j _ {1}: j _ {2}\right) \right\| _ {p} \leq \left\| A \right\| _ {p}. \tag {2.3.13}
$$

The proofs of these relationships are left as exercises. We mention that a sequence $\{ A ^ { ( k ) } \} \in \mathbb { R } ^ { m \times n }$ converges if there exists a matrix $A \in \mathbb { R } ^ { m \times n }$ such that

$$
\lim _ {k \to \infty} \| A ^ {(k)} - A \| = 0.
$$

The choice of norm is immaterial since all norms on $\mathbb { R } ^ { m \times n }$ are equivalent.

# 2.3.3 The Matrix 2-Norm

A nice feature of the matrix 1-norm and the matrix ∞-norm is that they are easy, $O ( n ^ { 2 } )$ computations. (See (2.3.9) and (2.3.10).) The calculation of the 2-norm is considerably more complicated.

Theorem 2.3.1. If $A \in \mathbb { R } ^ { m \times n }$ , then there exists a unit 2-norm n-vector z such that $A ^ { T } A z = \mu ^ { 2 } z$ where $\mu = \parallel A \parallel _ { 2 }$ .

Proof. Suppose $z \in \mathbb { R } ^ { n }$ is a unit vector such that $\parallel A z \parallel _ { 2 } = \parallel A \parallel _ { 2 }$ . Since z maximizes the function

$$
g (x) = \frac {1}{2} \frac {\parallel A x \parallel_ {2} ^ {2}}{\parallel x \parallel_ {2} ^ {2}} = \frac {1}{2} \frac {x ^ {T} A ^ {T} A x}{x ^ {T} x}
$$

it follows that it satisfies $\nabla g ( z ) = 0$ where $\nabla g$ is the gradient of g. A tedious differentiation shows that for $i = 1 { : } n$

$$
\frac {\partial g (z)}{\partial z _ {i}} = \left[ (z ^ {T} z) \sum_ {j = 1} ^ {n} (A ^ {T} A) _ {i j} z _ {j} - (z ^ {T} A ^ {T} A z) z _ {i} \right] \bigg / (z ^ {T} z) ^ {2}.
$$

In vector notation this says that $A ^ { T } A z = ( z ^ { T } A ^ { T } A z ) z$ . The theorem follows by setting $\mu = \parallel A z \parallel _ { 2 }$ .

The theorem implies that $\parallel A \parallel _ { 2 } ^ { 2 }$ is a zero of $p ( \lambda ) = \mathsf { d e t } ( A ^ { T } A - \lambda I )$ . In particular,

$$
\parallel A \parallel_ {2} = \sqrt {\lambda_ {\max} (A ^ {T} A)}
$$

We have much more to say about eigenvalues in Chapters 7 and 8. For now, we merely observe that 2-norm computation is iterative and a more involved calculation than those of the matrix 1-norm or ∞-norm. Fortunately, if the object is to obtain an order-of-magnitude estimate of $\parallel A \parallel _ { 2 }$ , then (2.3.7), (2.3.8), (2.3.11), or (2.3.12) can be used.

As another example of norm analysis, here is a handy result for 2-norm estimation.

Corollary 2.3.2. If $A \in \mathbb { R } ^ { m \times n }$ , then $\| A \| _ { 2 } \leq \sqrt { \| A \| _ { 1 } \| A \| _ { \infty } }$

Proof. $\mathrm { ~ I f ~ } z ~ \ne ~ 0$ is such that $A ^ { T } A z \ = \ \mu ^ { 2 } z$ with $\mu ~ = ~ \parallel A \parallel _ { 2 }$ , then $\textstyle \mu ^ { 2 } \parallel { \boldsymbol { z } } \parallel _ { 1 } =$ $\begin{array} { r } { \| A ^ { T } A z \| _ { 1 } \leq \| A ^ { T } \| _ { 1 } \| A \| _ { 1 } \| z \| _ { 1 } = \| A \| _ { \infty } \| A \| _ { 1 } \| z \| _ { 1 } . \quad \mathbb { D } } \end{array}$

# 2.3.4 Perturbations and the Inverse

We frequently use norms to quantify the effect of perturbations or to prove that a sequence of matrices converges to a specified limit. As an illustration of these norm applications, let us quantify the change in $A ^ { - 1 }$ as a function of change in A.

Lemma 2.3.3. If $F \in \mathbb { R } ^ { n \times n }$ and $\| \boldsymbol { F } \| _ { p } < 1$ , then $I - F$ is nonsingular and

$$
(I - F) ^ {- 1} = \sum_ {k = 0} ^ {\infty} F ^ {k}
$$

with

$$
\| (I - F) ^ {- 1} \| _ {p} \leq \frac {1}{1 - \| F \| _ {p}}.
$$

Proof. Suppose $I - F$ is singular. It follows that $( I - F ) x = 0$ for some nonzero x. But then $\| x \| _ { p } = \| F x \| _ { p }$ implies  $F \parallel _ { p } \geq 1$ , a contradiction. Thus, $I - F$ is nonsingular. To obtain an expression for its inverse consider the identity

$$
\left(\sum_ {k = 0} ^ {N} F ^ {k}\right) (I - F) = I - F ^ {N + 1}.
$$

Since $\| \ F \| _ { p } < 1$ it follows that $\operatorname* { l i m } _ { k \to \infty } F ^ { k } = 0$ because $\begin{array} { r } { \| \boldsymbol { F } ^ { k } \| _ { p } \leq \| \boldsymbol { F } \| _ { p } ^ { k } } \end{array}$ . Thus,

$$
\left(\lim _ {N \rightarrow \infty} \sum_ {k = 0} ^ {N} F ^ {k}\right) (I - F) = I.
$$

It follows that $( I - F ) ^ { - 1 } = \operatorname* { l i m } _ { N \to \infty } \sum _ { k = 0 } ^ { N } F ^ { k }$ . From this it is easy to show that

$$
\| (I - F) ^ {- 1} \| _ {p} \leq \sum_ {k = 0} ^ {\infty} \| F \| _ {p} ^ {k} = \frac {1}{1 - \| F \| _ {p}}
$$

completing the proof of the theorem.

![](images/golub_050_099__6cfa699ebce92d9436c4211b9f5d8bbea0e71ee4f692d4d62cbe202a66daa824.jpg)

Note that $\begin{array} { r } { \| ( I - F ) ^ { - 1 } - I \| _ { p } \ \leq \ \| F \| _ { p } / ( 1 - \| F \| _ { p } ) } \end{array}$ is a consequence of the lemma. Thus, $ { \mathrm { i f } } \epsilon \ll 1$ , then $O ( \epsilon )$ perturbations to the identity matrix induce $O ( \epsilon )$ perturbations in the inverse. In general, we have

Theorem 2.3.4. If A is nonsingular and $r \equiv \parallel A ^ { - 1 } E \parallel _ { p } < 1$ , then $A { + } E$ is nonsingular and

$$
\parallel (A + E) ^ {- 1} - A ^ {- 1} \parallel_ {p} \leq \frac {\parallel E \parallel_ {p} \parallel A ^ {- 1} \parallel_ {p} ^ {2}}{1 - r}.
$$

Proof. Note that $A + E = ( I + F ) A$ where $F = - E A ^ { - 1 }$ . Since  $F \parallel _ { p } = r < 1$ , it follows from Lemma 2.3.3 that $I + F$ is nonsingular and $\| { \bf \zeta } ( I + { \cal F } ) ^ { - 1 } \| _ { p } \dot { \leq } 1 / ( 1 - r )$ .

Thus, $( A + E ) ^ { - 1 } = A ^ { - 1 } ( I + F ) ^ { - 1 }$ is nonsingular and

$$
(A + E) ^ {- 1} - A ^ {- 1} = A ^ {- 1} (A - (A + E)) (A + E) ^ {- 1} = - A ^ {- 1} E A ^ {- 1} (I + F) ^ {- 1}.
$$

The theorem follows by taking norms.

# 2.3.5 Orthogonal Invariance

If $A \in \mathbb { R } ^ { m \times n }$ and the matrices $Q \in \mathbb { R } ^ { m \times m }$ and $Z \in \mathbb { R } ^ { n \times n }$ are orthogonal, then

$$
\left\| Q A Z \right\| _ {F} = \left\| A \right\| _ {F} \tag {2.3.14}
$$

and

$$
\| Q A Z \| _ {2} = \| A \| _ {2}. \tag {2.3.15}
$$

These properties readily follow from the orthogonal invariance of the vector 2-norm. For example,

$$
\parallel Q A \parallel_ {F} ^ {2} = \sum_ {j = 1} ^ {n} \parallel Q A (:, j) \parallel_ {2} ^ {2} = \sum_ {j = 1} ^ {n} \parallel A (:, j) \parallel_ {2} ^ {2} = \parallel A \parallel_ {F} ^ {2}
$$

and so $\left\| Q ( A Z ) \right\| _ { F } ^ { 2 } = \left\| \left( A Z \right) \right\| _ { F } ^ { 2 } = \left\| Z ^ { T } A ^ { T } \right\| _ { F } ^ { 2 } = \left\| A ^ { T } \right\| _ { F } ^ { 2 } = \left\| A \right\| _ { F } ^ { 2 } .$

# Problems

P2.3.1 Show $\left. A B \right. _ { p } \leq \left. A \right. _ { p } \left. B \right. _ { p }$ where $1 \leq p \leq \infty$ .

P2.3.2 Let B be any submatrix of A. Show that $\| B \| _ { p } \leq \| A \| _ { p } .$

P2.3.3 Show that if $D = \operatorname { d i a g } ( \mu _ { 1 } , \dots , \mu _ { k } ) \in \mathbb { R } ^ { m \times n }$ with $k = \operatorname* { m i n } \{ m , n \}$ , then $\| D \| _ { p } = \operatorname* { m a x } | \mu _ { i } |$

P2.3.4 Verify (2.3.7) and (2.3.8).

P2.3.5 Verify (2.3.9) and (2.3.10).

P2.3.6 Verify (2.3.11) and (2.3.12).

P2.3.7 Show that if $0 \neq s \in \mathbb { R } ^ { n }$ and $E \in \mathbb { R } ^ { n \times n }$ , then

$$
\left\| E \left(I - \frac {s s ^ {T}}{s ^ {T} s}\right) \right\| _ {F} ^ {2} = \parallel E \parallel_ {F} ^ {2} - \frac {\parallel E s \parallel_ {2} ^ {2}}{s ^ {T} s}.
$$

P2.3.8 Suppose $u \in \mathbb { R } ^ { m }$ and $v \in \mathbb { R } ^ { n }$ . Show that if $E = u v ^ { T }$ , then $\parallel E \parallel _ { F } = \parallel E \parallel _ { 2 } = \parallel u \parallel _ { 2 } \parallel v \parallel _ { 2 }$ and $\left\| \textbf { } E \right\| _ { \infty } \leq \left\| \textbf { } u \right\| _ { \infty } \left\| \textbf { } v \right\| _ { 1 } .$

P2.3.9 Suppose $A \in \mathbb { R } ^ { m \times n } , ~ y \in \mathbb { R } ^ { m }$ , and $0 \neq s \in \mathbb { R } ^ { n }$ . Show that $E = ( y - A s ) s ^ { T } / s ^ { T } ;$ s has the smallest 2-norm of all m-by-n matrices E that satisfy $( A + E ) s = y .$ .

P2.3.10 Verify that there exists a scalar $c > 0$ such that

$$
\| A \| _ {\Delta , c} = \max _ {i, j} c | a _ {i j} |
$$

satisfies the submultiplicative property (2.3.4) for matrix norms on $\mathbb { R } ^ { n \times n }$ . What is the smallest value for such a constant? Referring to this value as $c _ { * }$ , exhibit nonzero matrices B and $C$ with the property that $\| B C \| _ { \Delta , c _ { * } } = \| B \| _ { \Delta , c _ { * } } \| C \| _ { \Delta , c _ { * } }$ .

P2.3.11 Show that if A and B are matrices, then $\parallel A \otimes B \parallel _ { F } = \parallel A \parallel _ { F } \parallel B \parallel _ { F }$ .


---

<!-- golub_100_149 -->

For further discussion of matrix norms, see Stewart (IMC) as well as:

F.L. Bauer and C.T. Fike (1960). “Norms and Exclusion Theorems,” Numer. Math. 2, 137–144.

L. Mirsky (1960). “Symmetric Gauge Functions and Unitarily Invariant Norms,” Quart. J. Math. 11, 50–59.

A.S. Householder (1964). The Theory of Matrices in Numerical Analysis, Dover Publications, New York.

N.J. Higham (1992). “Estimating the Matrix p-Norm,” Numer. Math. 62, 539–556.

# 2.4 The Singular Value Decomposition

It is fitting that the first matrix decomposition that we present in the book is the singular value decomposition (SVD). The practical and theoretical importance of the SVD is hard to overestimate. It has a prominent role to play in data analysis and in the characterization of the many matrix “nearness problems.”

# 2.4.1 Derivation

The SVD is an orthogonal matrix reduction and so the 2-norm and Frobenius norm figure heavily in this section. Indeed, we can prove the existence of the decomposition using some elementary facts about the 2-norm developed in the previous two sections.

Theorem 2.4.1 (Singular Value Decomposition ). If A is a real m-by-n matrix, then there exist orthogonal matrices

$$
U = \left[ u _ {1} \mid \dots \mid u _ {m} \right] \in \mathbb {R} ^ {m \times m} \quad a n d \quad V = \left[ v _ {1} \mid \dots \mid v _ {n} \right] \in \mathbb {R} ^ {n \times n}
$$

such that

$$
U ^ {T} A V = \Sigma = \operatorname{diag} \left(\sigma_ {1}, \dots , \sigma_ {p}\right) \in \mathbb {R} ^ {m \times n}, \quad p = \min \{m, n \},
$$

where $\sigma _ { 1 } \geq \sigma _ { 2 } \geq . . . \geq \sigma _ { p } \geq 0 .$

Proof. Let $\boldsymbol { x } \in \mathbb { R } ^ { n }$ and $\boldsymbol { y } \in \mathbb { R } ^ { m }$ be unit 2-norm vectors that satisfy $A x = \sigma y$ with $\sigma = \parallel A \parallel _ { 2 }$ . From Theorem 2.1.1 there exist $V _ { 2 } \in \mathbb { R } ^ { n \times ( n - 1 ) }$ and $\bar { U _ { 2 } } \in \mathbb { R } ^ { m \times ( m - 1 ) }$ so $V = \left[ \ b { x } \mid V _ { 2 } \right] \in \mathbb { R } ^ { n \times n }$ and $U = \left[ \left. y \right| U _ { 2 } \right] \in \mathbb { R } ^ { m \times m }$ are orthogonal. It is not hard to show that

$$
U ^ {T} A V = \left[ \begin{array}{c c} \sigma & w ^ {T} \\ 0 & B \end{array} \right] \equiv A _ {1}
$$

where $w \in \mathbb { R } ^ { n - 1 }$ and $B \in \mathbb { R } ^ { ( m - 1 ) \times ( n - 1 ) }$ . Since

$$
\left\| A _ {1} \left(\left[ \begin{array}{c} \sigma \\ w \end{array} \right]\right) \right\| _ {2} ^ {2} \geq (\sigma^ {2} + w ^ {T} w) ^ {2}
$$

we have $\parallel A _ { 1 } \parallel _ { 2 } ^ { 2 } \geq ( \sigma ^ { 2 } + w ^ { T } w )$ . But $\sigma ^ { 2 } = \parallel A \parallel _ { 2 } ^ { 2 } = \parallel A _ { 1 } \parallel _ { 2 } ^ { 2 }$ , and so we must have $w = 0$ . An obvious induction argument completes the proof of the theorem.

The $\sigma _ { i }$ are the singular values of A, the $u _ { i }$ are the left singular vectors of A, and the $v _ { i }$ are right singular vectors of A. Separate visualizations of the SVD are required depending upon whether A has more rows or columns. Here are the 3-by-2 and 2-by-3 examples:

$$
\left[ \begin{array}{l l l} u _ {1 1} & u _ {1 2} & u _ {1 3} \\ u _ {2 1} & u _ {2 2} & u _ {2 3} \\ u _ {3 1} & u _ {3 2} & u _ {3 3} \end{array} \right] ^ {T} \left[ \begin{array}{l l} a _ {1 1} & a _ {1 2} \\ a _ {2 1} & a _ {2 2} \\ a _ {3 1} & a _ {3 2} \end{array} \right] \left[ \begin{array}{l l} v _ {1 1} & v _ {1 2} \\ v _ {2 1} & v _ {2 2} \end{array} \right] = \left[ \begin{array}{l l} \sigma_ {1} & 0 \\ 0 & \sigma_ {2} \\ 0 & 0 \end{array} \right],
$$

$$
\left[ \begin{array}{c c} u _ {1 1} & u _ {1 2} \\ u _ {2 1} & u _ {2 2} \end{array} \right] ^ {T} \left[ \begin{array}{c c c} a _ {1 1} & a _ {1 2} & a _ {1 3} \\ a _ {2 1} & a _ {2 2} & a _ {2 3} \end{array} \right] \left[ \begin{array}{c c c} v _ {1 1} & v _ {1 2} & v _ {1 3} \\ v _ {2 1} & v _ {2 2} & v _ {2 3} \\ v _ {3 1} & v _ {3 2} & v _ {3 3} \end{array} \right] = \left[ \begin{array}{c c c} \sigma_ {1} & 0 & 0 \\ 0 & \sigma_ {2} & 0 \end{array} \right].
$$

In later chapters, the notation $\sigma _ { i } ( A )$ is used to designate the ith largest singular value of a matrix A. The largest and smallest singular values are important and for them we also have a special notation:

$$
\sigma_ {\max} (A) = \text { the   largest   singular   value   of   matrix } A,
$$

$$
\sigma_ {\min} (A) = \text { the   smallest   singular   value   of   matrix } A.
$$

# 2.4.2 Properties

We establish a number of important corollaries to the SVD that are used throughout the book.

Corollary 2.4.2. $I f U ^ { T } A V = \Sigma$ is the SVD of $A \in \mathbb { R } ^ { m \times n }$ and $m \geq n$ , then $f o r i = 1$ :n $A v _ { i } = \sigma _ { i } u _ { i }$ and $A ^ { T } u _ { i } = \sigma _ { i } v _ { i }$ .

Proof. Compare columns in $A V = U \Sigma$ and $A ^ { T } U = V \Sigma ^ { T }$ .□

There is a nice geometry behind this result. The singular values of a matrix A are the lengths of the semiaxes of the hyperellipsoid E defined by $E = \left\{ \ A x : \| \ x \ \| _ { 2 } = 1 \right\}$ . The semiaxis directions are defined by the $u _ { i }$ and their lengths are the singular values.

It follows immediately from the corollary that

$$
A ^ {T} A v _ {i} = \sigma_ {i} ^ {2} v _ {i}, \tag {2.4.1}
$$

$$
A A ^ {T} u _ {i} = \sigma_ {i} ^ {2} u _ {i} \tag {2.4.2}
$$

for $i = 1 { : } n$ . This shows that there is an intimate connection between the SVD of A and the eigensystems of the symmetric matrices $A ^ { T } A$ and $A A ^ { T }$ . See §8.6 and §10.4.

The 2-norm and the Frobenius norm have simple SVD characterizations.

Corollary 2.4.3. If $A \in \mathbb { R } ^ { m \times n }$ , then

$$
\| A \| _ {2} = \sigma_ {1}, \qquad \| A \| _ {F} = \sqrt {\sigma_ {1} ^ {2} + \cdots + \sigma_ {p} ^ {2}},
$$

where $p = \operatorname* { m i n } \{ m , n \}$ .

Proof. These results follow immediately from the fact that $\| U ^ { T } A V \| = \| \Sigma \|$ for both the 2-norm and the Frobenius norm. □

We show in §8.6 that if A is perturbed by a matrix E, then no singular value can move by more than $\parallel E \parallel _ { 2 }$ . The following corollary identifies two useful instances of this result.

Corollary 2.4.4. If $A \in \mathbb { R } ^ { m \times n }$ and $E \in \mathbb { R } ^ { m \times n }$ , then

$$
\sigma_ {\max} (A + E) \leq \sigma_ {\max} (A) + \| E \| _ {2},
$$

$$
\sigma_ {\min} (A + E) \geq \sigma_ {\min} (A) - \parallel E \parallel_ {2}.
$$

Proof. Using Corollary 2.4.2 it is easy to show that

$$
\sigma_ {\min} (A) \cdot \| x \| _ {2} \leq \| A x \| _ {2} \leq \sigma_ {\max} (A) \cdot \| x \| _ {2}.
$$

The required inequalities follow from this result.

If a column is added to a matrix, then the largest singular value increases and the smallest singular value decreases.

Corollary 2.4.5. If $A \in \mathbb { R } ^ { m \times n }$ , $m > n$ , and $z \in \mathbb { R } ^ { m }$ , then

$$
\sigma_ {\max} \big (\left[ A \mid z \right] \big) \geq \sigma_ {\max} (A),
$$

$$
\sigma_ {\min} \left(\left[ A \mid z \right]\right) \leq \sigma_ {\min} (A).
$$

Proof. Suppose $A = U \Sigma V ^ { T }$ is the SVD of A and let $x = V ( : , 1 )$ and $\tilde { A } = [ A | z ]$ . Using Corollary 2.4.4, we have

$$
\sigma_ {\max} (A) = \| A x \| _ {2} = \left\| \tilde {A} \left[ \begin{array}{c} x \\ 0 \end{array} \right] \right\| _ {2} \leq \sigma_ {\max} (\tilde {A}).
$$

The proof that $\sigma _ { \operatorname* { m i n } } ( A ) \geq \sigma _ { \operatorname* { m i n } } ( \tilde { A } )$ is similar.

The SVD neatly characterizes the rank of a matrix and orthonormal bases for both its nullspace and its range.

Corollary 2.4.6. If A has r positive singular values, then rank(A) = r and

$$
\operatorname{null} (A) = \operatorname{span} \left\{v _ {r + 1}, \dots , v _ {n} \right\},
$$

$$
\operatorname{ran} (A) = \operatorname{span} \{u _ {1}, \dots , u _ {r} \}.
$$

Proof. The rank of a diagonal matrix equals the number of nonzero diagonal entries. Thus, rank $( A ) = { \mathsf { r a n k } } ( \Sigma ) = r$ . The assertions about the nullspace and range follow from Corollary 2.4.2.

If A has rank r, then it can be written as the sum of r rank-1 matrices. The SVD gives us a particularly nice choice for this expansion.

Corollary 2.4.7. If $A \in \mathbb { R } ^ { m \times n }$ and rank(A) = r, then

$$
A = \sum_ {i = 1} ^ {r} \sigma_ {i} u _ {i} v _ {i} ^ {T}.
$$

Proof. This is an exercise in partitioned matrix multiplication:

$$
(U \Sigma) V ^ {T} = \left(\left[ \sigma_ {1} u _ {1} \mid \sigma_ {2} u _ {2} \mid \dots \mid \sigma_ {r} u _ {r} \mid 0 \mid \dots \mid 0 \right]\right) \left[ \begin{array}{c} v _ {1} ^ {T} \\ \vdots \\ v _ {n} ^ {T} \end{array} \right] = \sum_ {i = 1} ^ {r} \sigma_ {i} u _ {i} v _ {i} ^ {T}. \quad \square
$$

The intelligent handling of rank degeneracy is an important topic that we discuss in Chapter 5. The SVD has a critical role to play because it can be used to identify nearby matrices of lesser rank.

Theorem 2.4.8 (The Eckhart-Young Theorem). If $k < r = \mathrm { r a n k } ( A )$ and

$$
A _ {k} = \sum_ {i = 1} ^ {k} \sigma_ {i} u _ {i} v _ {i} ^ {T}, \tag {2.4.3}
$$

then

$$
\min _ {\operatorname{rank} (B) = k} \| A - B \| _ {2} = \| A - A _ {k} \| _ {2} = \sigma_ {k + 1}. \tag {2.4.4}
$$

Proof. Since $U ^ { T } A _ { k } V = \mathrm { d i a g } ( \sigma _ { 1 } , \dots , \sigma _ { k } , 0 , \dots , 0 )$ it follows that $A _ { k }$ is rank k. Moreover, $U ^ { T } ( A - A _ { k } ) V = \mathrm { d i a g } ( 0 , \dots , 0 , \sigma _ { k + 1 } , \dots , \sigma _ { p } )$ and so $\parallel A - A _ { k } \parallel _ { 2 } = \sigma _ { k + 1 }$ .

Now suppose $\mathsf { r a n k } ( B ) = k$ for some $\boldsymbol { B } \in \mathbb { R } ^ { m \times n }$ . It follows that we can find orthonormal vectors $x _ { 1 } , \ldots , x _ { n - k }$ so nul $\left| \left( B \right) = { \mathsf { s p a n } } \{ x _ { 1 } , \ldots , x _ { n - k } \} \right.$ . A dimension argument shows that

$$
\operatorname{span} \left\{x _ {1}, \dots , x _ {n - k} \right\} \cap \operatorname{span} \left\{v _ {1}, \dots , v _ {k + 1} \right\} \neq \{0 \}.
$$

Let $z$ be a unit 2-norm vector in this intersection. Since $B z = 0$ and

$$
A z = \sum_ {i = 1} ^ {k + 1} \sigma_ {i} (v _ {i} ^ {T} z) u _ {i},
$$

we have

$$
\parallel A - B \parallel_ {2} ^ {2} \geq \parallel (A - B) z \parallel_ {2} ^ {2} = \parallel A z \parallel_ {2} ^ {2} = \sum_ {i = 1} ^ {k + 1} \sigma_ {i} ^ {2} (v _ {i} ^ {T} z) ^ {2} \geq \sigma_ {k + 1} ^ {2},
$$

completing the proof of the theorem.

Note that this theorem says that the smallest singular value of A is the 2-norm distance of A to the set of all rank-deficient matrices. We also mention that the matrix $A _ { k }$ defined in (2.4.3) is the closest rank-k matrix to A in the Frobenius norm.

# 2.4.3 The Thin SVD

If $A = U \Sigma V ^ { T } \in \mathbb { R } ^ { m \times n }$ is the SVD of A and $m \geq n$ , then

$$
A = U _ {1} \Sigma_ {1} V ^ {T}
$$

where

$$
U _ {1} = U (:, 1: n) = \left[ u _ {1} \mid \dots \mid u _ {n} \right] \in \mathbb {R} ^ {m \times n}
$$

and

$$
\Sigma_ {1} = \Sigma (1: n, 1: n) = \operatorname{diag} \left(\sigma_ {1}, \dots , \sigma_ {n}\right) \in \mathbb {R} ^ {n \times n}.
$$

We refer to this abbreviated version of the SVD as the thin SVD.

# 2.4.4 Unitary Matrices and the Complex SVD

Over the complex field the unitary matrices correspond to the orthogonal matrices. In particular, $Q \in \mathbb { C } ^ { n \times n }$ is unitary if $Q ^ { H } Q = Q Q ^ { \bar { H } } = I _ { n }$ . Unitary transformations preserve both the 2-norm and the Frobenius norm. The SVD of a complex matrix involves unitary matrices. If $A \in \mathbb { C } ^ { m \times n }$ , then there exist unitary matrices $U \in \mathbb { C } ^ { m \times m }$ and $V \in \mathbb { C } ^ { n \times n }$ such that

$$
U ^ {H} A V = \mathrm{diag} (\sigma_ {1}, \ldots , \sigma_ {p}) \in \mathbb {R} ^ {m \times n} \qquad p = \min \{m, n \}
$$

where $\sigma _ { 1 } \geq \sigma _ { 2 } \geq . . . \geq \sigma _ { p } \geq 0$ . All of the real SVD properties given above have obvious complex analogs.

# Problems

P2.4.1 Show that if $Q = Q _ { 1 } + i Q _ { 2 }$ is unitary with $Q _ { 1 } , Q _ { 2 } \in \mathbb { R } ^ { n \times n }$ , then the 2n-by-2n real matrix

$$
Z = \left[ \begin{array}{c c} Q _ {1} & - Q _ {2} \\ Q _ {2} & Q _ {1} \end{array} \right]
$$

is orthogonal.

P2.4.2 Prove that if $A \in \mathbb { R } ^ { m \times n }$ , then

$$
\sigma_ {\max} (A) = \max _ { \begin{array}{c} y \in \mathbb {R} ^ {m} \\ x \in \mathbb {R} ^ {n} \end{array} } \frac {y ^ {T} A x}{\| x \| _ {2} \| y \| _ {2}}.
$$

P2.4.3 For the 2-by-2 matrix $A \ = \ \left[ \begin{array} { c c } { { w } } & { { x } } \\ { { y } } & { { z } } \end{array} \right]$ , derive expressions for $\sigma _ { \operatorname* { m a x } } ( A )$ and $\sigma _ { \mathrm { m i n } } ( A )$ that are functions of $w , x , y ,$ and z.

P2.4.4 Show that any matrix in $\mathbb { R } ^ { m \times n }$ is the limit of a sequence of full rank matrices.

P2.4.5 Show that if $A \in \mathbb { R } ^ { m \times n }$ has rank n, then $\parallel A ( A ^ { T } A ) ^ { - 1 } A ^ { T } \parallel _ { 2 } = 1$ .

P2.4.6 What is the nearest rank-1 matrix to

$$
A = \left[ \begin{array}{c c} 1 & M \\ 0 & 1 \end{array} \right]
$$

in the Frobenius norm?

P2.4.7 Show that if $A \in \mathbb { R } ^ { m \times n }$ , then $\parallel A \parallel _ { F } \leq \sqrt { \mathsf { r a n k } ( A ) } \parallel A \parallel _ { 2 }$ , thereby sharpening (2.3.7).

P2.4.8 Suppose $A \in \mathbb { R } ^ { n \times n }$ . Give an SVD solution to the following problem:

$$
\min _ {\det (B) = | \det (A) |} \| A - B \| _ {F}.
$$

P2.4.9 Show that if a nonzero row is added to a matrix, then both the largest and smallest singular values increase.

P2.4.10 Show that if $\theta _ { u }$ and $\theta _ { v }$ are real numbers and

$$
A = \left[ \begin{array}{l l} \cos (\theta_ {u}) & \sin (\theta_ {u}) \\ \cos (\theta_ {v}) & \sin (\theta_ {v}) \end{array} \right],
$$

then $U ^ { T } A V = \Sigma$ where

$$
U = \left[ \begin{array}{c c} \cos (\pi / 4) & - \sin (\pi / 4) \\ \sin (\pi / 4) & \cos (\pi / 4) \end{array} \right], V = \left[ \begin{array}{c c} \cos (a) & - \sin (a) \\ \sin (a) & \cos (a) \end{array} \right],
$$

and Σ = diag(√2 cos(b), √2 sin(b)) with a = (θv + θu)/2 and $b = ( \theta _ { v } - \theta _ { u } ) / 2$ .

# Notes and References for 2.4

Forsythe and Moler (SLAS) offer a good account of the SVD’s role in the analysis of the Ax = b problem. Their proof of the decomposition is more traditional than ours in that it makes use of the eigenvalue theory for symmetric matrices. Historical SVD references include:

E. Beltrami (1873). “Sulle Funzioni Bilineari,” Gionale di Mathematiche 11, 98–106.   
C. Eckart and G. Young (1939). “A Principal Axis Transformation for Non-Hermitian Matrices,” Bull. AMS 45, 118–21.   
G.W. Stewart (1993). “On the Early History of the Singular Value Decomposition,” SIAM Review 35, 551–566.

One of the most significant developments in scientific computation has been the increased use of the SVD in application areas that require the intelligent handling of matrix rank. This work started with:

C. Eckart, and G. Young (1936). “The Approximation of One Matrix by Another of Lower Rank,” Psychometrika 1, 211–218.

For generalizations of the SVD to infinite dimensional Hilbert space, see:

I.C. Gohberg and M.G. Krein (1969). Introduction to the Theory of Linear Non-Self Adjoint Operators, Amer. Math. Soc., Providence, RI.   
F. Smithies (1970). Integral Equations, Cambridge University Press, Cambridge.

Reducing the rank of a matrix as in Corollary 2.4.6 when the perturbing matrix is constrained is discussed in:

J.W. Demmel (1987). “The Smallest Perturbation of a Submatrix which Lowers the Rank and Constrained Total Least Squares Problems, SIAM J. Numer. Anal. 24, 199–206.   
G.H. Golub, A. Hoffman, and G.W. Stewart (1988). “A Generalization of the Eckart-Young-Mirsky Approximation Theorem.” Lin. Alg. Applic. 88/89, 317–328.   
G.A. Watson (1988). “The Smallest Perturbation of a Submatrix which Lowers the Rank of the Matrix,” IMA J. Numer. Anal. 8, 295–304.

# 2.5 Subspace Metrics

If the object of a computation is to compute a matrix or a vector, then norms are useful for assessing the accuracy of the answer or for measuring progress during an iteration. If the object of a computation is to compute a subspace, then to make similar comments we need to be able to quantify the distance between two subspaces. Orthogonal projections are critical in this regard. After the elementary concepts are established we discuss the CS decomposition. This is an SVD-like decomposition that is handy when we have to compare a pair of subspaces.

# 2.5.1 Orthogonal Projections

Let $S \subseteq \mathbb { R } ^ { n }$ be a subspace. $P \in \mathbb { R } ^ { n \times n }$ is the orthogonal projection onto S if ran $( P ) = S$ , $P ^ { 2 } = P ,$ , and $P ^ { T } = P$ . From this definition it is easy to show that if $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , then $P x \in S$ and $( I - P ) x \in S ^ { \perp }$ .

If $P _ { 1 }$ and $P _ { 2 }$ are each orthogonal projections, then for any $z \in \mathbb { R } ^ { n }$ we have

$$
\| (P _ {1} - P _ {2}) z \| _ {2} ^ {2} = (P _ {1} z) ^ {T} (I - P _ {2}) z + (P _ {2} z) ^ {T} (I - P _ {1}) z.
$$

If $\mathsf { r a n } ( P _ { 1 } ) = \mathsf { r a n } ( P _ { 2 } ) = S$ , then the right-hand side of this expression is zero, showing that the orthogonal projection for a subspace is unique. If the columns of $V =$ $\left[ \left. v _ { 1 } \right| \cdots \right| \left. v _ { k } \right]$ are an orthonormal basis for a subspace S, then it is easy to show that ${ \cal P } = { \cal V } { \cal V } ^ { T }$ is the unique orthogonal projection onto S. Note that if $v \in \mathbb { R } ^ { n }$ , then $P = v v ^ { T } / v ^ { T } v$ is the orthogonal projection onto $S = \mathsf { s p a n } \{ v \}$ .

# 2.5.2 SVD-Related Projections

There are several important orthogonal projections associated with the singular value decomposition. Suppose $A = U \Sigma V ^ { T } \in \mathbb { R } ^ { m \times n }$ is the SVD of A and that $r = \mathsf { r a n k } ( A )$ . If we have the U and V partitionings

$$
U   =   \left[ \begin{array}{c c} U _ {r} & \tilde {U} _ {r} \\ r & m - r \end{array} \right]    , \qquad V   =   \left[ \begin{array}{c c} V _ {r} & \tilde {V} _ {r} \\ r & n - r \end{array} \right]    ,
$$

then

$$
V _ {r} V _ {r} ^ {T} = \text { projection   on   to   } \text { null } (A) ^ {\perp} = \text { ran } (A ^ {T}),
$$

$$
\tilde {V} _ {r} \tilde {V} _ {r} ^ {T} = \text { projection   on   to   } \mathsf {n u l l} (A),
$$

$$
U _ {r} U _ {r} ^ {T} = \text { projection   on   to   } \operatorname{ran} (A),
$$

$$
\tilde {U} _ {r} \tilde {U} _ {r} ^ {T} = \text { projection   on   to } \operatorname{ran} (A) ^ {\perp} = \operatorname{null} (A ^ {T}).
$$

# 2.5.3 Distance Between Subspaces

The one-to-one correspondence between subspaces and orthogonal projections enables us to devise a notion of distance between subspaces. Suppose $S _ { 1 }$ and $S _ { 2 }$ are subspaces of $\mathbb { R } ^ { n }$ and that dim $( S _ { 1 } ) = \mathsf { d i m } ( S _ { 2 } )$ . We define the distance between these two spaces by

$$
\mathrm{dist} (S _ {1}, S _ {2}) = \| P _ {1} - P _ {2} \| _ {2} \tag {2.5.1}
$$

where $P _ { i }$ is the orthogonal projection onto $S _ { i }$ . The distance between a pair of subspaces can be characterized in terms of the blocks of a certain orthogonal matrix.

Theorem 2.5.1. Suppose

$$
W = \left[ \begin{array}{c c} W _ {1} & W _ {2} \\ k & n - k \end{array} \right], Z = \left[ \begin{array}{c c} Z _ {1} & Z _ {2} \\ k & n - k \end{array} \right],
$$

are $n { - } b y { - } n$ orthogonal matrices. $I f S _ { 1 } = \mathsf { r a n } ( W _ { 1 } )$ and $S _ { 2 } = \mathsf { r a n } ( Z _ { 1 } )$ , then

$$
\operatorname{dist} \left(S _ {1}, S _ {2}\right) = \left\| W _ {1} ^ {T} Z _ {2} \right\| _ {2} = \left\| Z _ {1} ^ {T} W _ {2} \right\| _ {2}.
$$

Proof. We first observe that

$$
\begin{array}{l} \operatorname{dist} \left(S _ {1}, S _ {2}\right) = \left\| W _ {1} W _ {1} ^ {T} - Z _ {1} Z _ {1} ^ {T} \right\| _ {2} \\ = \left\| W ^ {T} \left(W _ {1} W _ {1} ^ {T} - Z _ {1} Z _ {1} ^ {T}\right) Z \right\| _ {2} \\ = \left\| \left[ \begin{array}{c c} 0 & W _ {1} ^ {T} Z _ {2} \\ - W _ {2} ^ {T} Z _ {1} & 0 \end{array} \right] \right\| _ {2}. \\ \end{array}
$$

Note that the matrices $W _ { 2 } ^ { T } Z _ { 1 }$ and $W _ { 1 } ^ { T } Z _ { 2 }$ are submatrices of the orthogonal matrix

$$
Q = \left[ \begin{array}{l l} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right] \equiv \left[ \begin{array}{l l} W _ {1} ^ {T} Z _ {1} & W _ {1} ^ {T} Z _ {2} \\ W _ {2} ^ {T} Z _ {1} & W _ {2} ^ {T} Z _ {2} \end{array} \right] = W ^ {T} Z. \tag {2.5.2}
$$

Our goal is to show that $\parallel Q _ { 2 1 } \parallel _ { 2 } = \parallel Q _ { 1 2 } \parallel _ { 2 }$ . Since Q is orthogonal it follows from

$$
Q \left[ \begin{array}{l} x \\ 0 \end{array} \right] = \left[ \begin{array}{l} Q _ {1 1} x \\ Q _ {2 1} x \end{array} \right]
$$

that $1 = \parallel Q _ { 1 1 } x \parallel _ { 2 } ^ { 2 } + \parallel Q _ { 2 1 } x \parallel _ { 2 } ^ { 2 }$ for all unit 2-norm $\boldsymbol { x } \in \mathbb { R } ^ { k }$ . Thus,

$$
\| Q _ {2 1} \| _ {2} ^ {2} = \max _ {\| x \| _ {2} = 1} \| Q _ {2 1} x \| _ {2} ^ {2} = 1 - \min _ {\| x \| _ {2} = 1} \| Q _ {1 1} x \| _ {2} ^ {2} = 1 - \sigma_ {\min} (Q _ {1 1}) ^ {2}.
$$

Analogously, by working with $Q ^ { T }$ (which is also orthogonal) it is possible to show that

$$
\parallel Q _ {1 2} ^ {T} \parallel_ {2} ^ {2} = 1 - \sigma_ {\min} (Q _ {1 1} ^ {T}) ^ {2},
$$

and therefore

$$
\parallel Q _ {1 2} \parallel_ {2} ^ {2} = 1 - \sigma_ {\min} (Q _ {1 1}) ^ {2}.
$$

Thus,  Q21 2 =  Q12 2.

Note that if $S _ { 1 }$ and $S _ { 2 }$ are subspaces in $\mathbb { R } ^ { n }$ with the same dimension, then

$$
0 \leq \operatorname{dist} \left(S _ {1}, S _ {2}\right) \leq 1.
$$

It is easy to show that

$$
\begin{array}{l} \operatorname{dist} \left(S _ {1}, S _ {2}\right) = 0 \Rightarrow S _ {1} = S _ {2}, \\ \operatorname{dist} \left(S _ {1}, S _ {2}\right) = 1 \Rightarrow S _ {1} \cap S _ {2} ^ {\perp} \neq \{0 \}. \\ \end{array}
$$

A more refined analysis of the blocks of the matrix Q in (2.5.2) sheds light on the difference between a pair of subspaces. A special, SVD-like decomposition for orthogonal matrices is required.

# 2.5.4 The CS Decomposition

The blocks of an orthogonal matrix partitioned into 2-by-2 form have highly related SVDs. This is the gist of the CS decomposition. We prove a very useful special case first.

Theorem 2.5.2 (The CS Decomposition (Thin Version)). Consider the matrix

$$
Q = \left[ \begin{array}{l} Q _ {1} \\ Q _ {2} \end{array} \right], \qquad \qquad Q _ {1} \in \mathbb {R} ^ {m _ {1} \times n _ {1}},   Q _ {2} \in \mathbb {R} ^ {m _ {2} \times n _ {1}},
$$

where $m _ { 1 } \geq n _ { 1 }$ and $m _ { 2 } \geq n _ { 1 }$ . If the columns of Q are orthonormal, then there exist orthogonal matrices $U _ { 1 } \in \mathbb { R } ^ { m _ { 1 } \times m _ { 1 } } , U _ { 2 } \in \mathbb { R } ^ { m _ { 2 } \times m _ { 2 } }$ , and $V _ { 1 } \in \mathbb { R } ^ { n _ { 1 } \times n _ { 1 } }$ such that

$$
\left[ \begin{array}{c c} U _ {1} & 0 \\ 0 & U _ {2} \end{array} \right] ^ {T} \left[ \begin{array}{c} Q _ {1} \\ Q _ {2} \end{array} \right] V _ {1} = \left[ \begin{array}{c} C \\ S \end{array} \right]
$$

where

$$
C _ {0} = \operatorname{diag} \left(\cos \left(\theta_ {1}\right), \dots , \cos \left(\theta_ {n _ {1}}\right)\right) \in \mathbb {R} ^ {m _ {1} \times n _ {1}},
$$

$$
S _ {0} = \operatorname{diag} \left(\sin \left(\theta_ {1}\right), \dots , \sin \left(\theta_ {n _ {1}}\right)\right) \in \mathbb {R} ^ {m _ {2} \times n _ {1}},
$$

and

$$
0 \leq \theta_ {1} \leq \theta_ {2} \leq \dots \leq \theta_ {n _ {1}} \leq \frac {\pi}{2}.
$$

Proof. Since $\| Q _ { 1 } \| _ { 2 } \leq \| Q \| _ { 2 } = 1$ , the singular values of $Q _ { 1 }$ are all in the interval [0, 1]. Let

$$
U _ {1} ^ {T} Q _ {1} V _ {1} = C _ {0} = \mathrm{diag} (c _ {1}, \ldots , c _ {n _ {1}}) = \left[ \begin{array}{c c} I _ {t} & 0 \\ 0 & \Sigma \end{array} \right] _ {m _ {1} - t} ^ {t}
$$

be the SVD of $Q _ { 1 }$ where we assume

$$
1 = c _ {1} = \dots = c _ {t} > c _ {t + 1} \geq \dots \geq c _ {n _ {1}} \geq 0.
$$

To complete the proof of the theorem we must construct the orthogonal matrix $U _ { 2 }$ . If

$$
Q _ {2} V _ {1} = \left[ \begin{array}{c c} W _ {1} & W _ {2} \\ t & n _ {1} - t \end{array} \right] \quad ,
$$

then

$$
\left[ \begin{array}{c c} U _ {1} & 0 \\ 0 & I _ {m _ {2}} \end{array} \right] ^ {T} \left[ \begin{array}{c} Q _ {1} \\ Q _ {2} \end{array} \right] V _ {1} = \left[ \begin{array}{c c} I _ {t} & 0 \\ 0 & \Sigma \\ W _ {1} & W _ {2} \end{array} \right].
$$

Since the columns of this matrix have unit 2-norm, $W _ { 1 } = 0$ . The columns of $W _ { 2 }$ are nonzero and mutually orthogonal because

$$
W _ {2} ^ {T} W _ {2} = I _ {n _ {1} - t} - \Sigma^ {T} \Sigma \equiv \mathrm{diag} (1 - c _ {t + 1} ^ {2}, \ldots , 1 - c _ {n _ {1}} ^ {2})
$$

is nonsingular. If $s _ { k } = \sqrt { 1 - c _ { k } ^ { 2 } }$ for $k = 1 { : } n _ { 1 }$ , then the columns of

$$
Z = W _ {2} \operatorname{diag} (1 / s _ {t + 1}, \dots , 1 / s _ {n})
$$

are orthonormal. By Theorem 2.1.1 there exists an orthogonal matrix $U _ { 2 } \in \mathbb { R } ^ { m _ { 2 } \times m _ { 2 } }$ with $U _ { 2 } ( : , t + 1 { : } n _ { 1 } ) = Z$ . It is easy to verify that

$$
U _ {2} ^ {T} Q _ {2} V _ {1} = \mathrm{diag} (s _ {1}, \ldots , s _ {n _ {1}}) \equiv S _ {0}.
$$

Since $c _ { k } ^ { 2 } + s _ { k } ^ { 2 } = 1$ for $k = 1 { : } n _ { 1 }$ , it follows that these quantities are the required cosines and sines. □

By using the same techniques it is possible to prove the following, more general version of the decomposition:

Theorem 2.5.3 (CS Decomposition). Suppose

$$
Q = = \left[ \begin{array}{c c} Q _ {1 1} & Q _ {1 2} \\ Q _ {2 1} & Q _ {2 2} \end{array} \right] _ {m _ {2}} ^ {m _ {1}}
$$

is a square orthogonal matrix and that $m _ { 1 } \geq n _ { 1 }$ and $m _ { 1 } \geq m _ { 2 }$ . Define the nonnegative integers p and q by $\textit { p } = \operatorname* { m a x } \{ 0 , n _ { 1 } - m _ { 2 } \}$ and q = max $\{ 0 , m _ { 2 } - n _ { 1 } \}$ . There exist orthogonal $U _ { 1 } \in \mathbb { R } ^ { m _ { 1 } \times m _ { 1 } }$ , $U _ { 2 } \in \mathbb { R } ^ { m _ { 2 } \times m _ { 2 } }$ , $V _ { 1 } \in \mathbb { R } ^ { n _ { 1 } \times n _ { 1 } }$ , and $V _ { 2 } \in \mathbb { R } ^ { n _ { 2 } \times n _ { 2 } }$ such that if

$$
U = \left[ \begin{array}{c c} U _ {1} & 0 \\ \hline 0 & U _ {2} \end{array} \right] \qquad \text { and } \qquad V = \left[ \begin{array}{c c} V _ {1} & 0 \\ \hline 0 & V _ {2} \end{array} \right],
$$

then

$$
U ^ {T} Q V = \left[ \begin{array}{c c c c c} I & 0 & 0 & 0 & 0 \\ 0 & C & S & 0 & 0 \\ 0 & 0 & 0 & 0 & I \\ \hline 0 & S & - C & 0 & 0 \\ 0 & 0 & 0 & I & 0 \\ \hline p & n _ {1} - p & n _ {1} - p & q & m _ {1} - n _ {1} \end{array} \right] \begin{array}{l} p \\ n _ {1} - p \\ m _ {1} - n _ {1} \\ n _ {1} - p \\ q \end{array}
$$

where

$$
C = \operatorname{diag} \left(\cos \left(\theta_ {p + 1}\right), \dots , \cos \left(\theta_ {n _ {1}}\right)\right) = \operatorname{diag} \left(c _ {p + 1}, \dots , c _ {n _ {1}}\right),
$$

$$
S = \operatorname{diag} \left(\sin \left(\theta_ {p + 1}\right), \dots , \sin \left(\theta_ {n _ {1}}\right)\right) = \operatorname{diag} \left(s _ {p + 1}, \dots , s _ {n _ {1}}\right),
$$

and $0 \le \theta _ { p + 1 } \le \cdots \le \theta _ { n _ { 1 } } \le \pi / 2$ .

Proof. See Paige and Saunders (1981) for details.

We made the assumptions $m _ { 1 } \geq n _ { 1 }$ and $m _ { 1 } \geq m _ { 2 }$ for clarity. Through permutation and transposition, any 2-by-2 block orthogonal matrix can be put into the form required by the theorem. Note that the blocks in the transformed Q, i.e., the $U _ { i } ^ { T } Q _ { i j } V _ { j }$ , are diagonal-like but not necessarily diagonal. Indeed, as we have presented it, the CS decomposition gives us four unnormalized SVDs. If $Q _ { 2 1 }$ has more rows than columns, then $p = 0$ and the reduction looks like this (for example):

$$
U ^ {T} Q V = \left[ \begin{array}{c c c c c c c} c _ {1} & 0 & s _ {1} & 0 & 0 & 0 & 0 \\ 0 & c _ {2} & 0 & s _ {2} & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 1 \\ \hline s _ {1} & 0 & - c _ {1} & 0 & 0 & 0 & 0 \\ 0 & s _ {2} & 0 & - c _ {2} & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 1 & 0 & 0 \end{array} \right].
$$

On the other hand, if $Q _ { 2 1 }$ has more columns than rows, then $q = 0$ and the decomposition has the form

$$
U ^ {T} Q V = \left[ \begin{array}{c c c c c} 1 & 0 & 0 & 0 & 0 \\ 0 & c _ {2} & 0 & s _ {2} & 0 \\ 0 & 0 & c _ {3} & 0 & s _ {3} \\ \hline 0 & s _ {2} & 0 & - c _ {2} & 0 \\ 0 & 0 & s _ {3} & 0 & - c _ {3} \end{array} \right].
$$

Regardless of the partitioning, the essential message of the CS decomposition is that the SVDs of the Q-blocks are highly related.

# Problems

P2.5.1 Show that if P is an orthogonal projection, then $Q = I - 2 P$ is orthogonal.

P2.5.2 What are the singular values of an orthogonal projection?

P2.5.3 Suppose $S _ { 1 } = \mathsf { s p a n } \{ x \}$ and $S _ { 2 } = \mathsf { s p a n } \{ y \}$ , where x and y are unit 2-norm vectors in $\mathbb { R } ^ { 2 }$ . Working only with the definition of dist $( \cdot , \cdot )$ , show that dis $\langle S _ { 1 } , S _ { 2 } \rangle = \sqrt { 1 - ( x ^ { T } y ) ^ { 2 } }$ , verifying that the distance between $S _ { 1 }$ and $S _ { 2 }$ equals the sine of the angle between x and y.

P2.5.4 Refer to §1.3.10. Show that if $Q \in \mathbb { R } ^ { 2 n \times 2 n }$ is orthogonal and symplectic, then Q has the form

$$
Q = \left[ \begin{array}{c c} Q _ {1} & Q _ {2} \\ - Q _ {2} & Q _ {1} \end{array} \right], \qquad Q _ {1}, Q _ {2} \in \mathbb {R} ^ {n \times n}.
$$

P2.5.5 Suppose $P \in \mathbb { R } ^ { n \times n }$ and $P ^ { 2 } = P$ . Show that $\| \ P \| _ { 2 } > 1$ if null(P ) is not a subspace of $\mathsf { r a n } ( A ) ^ { \perp }$ . Such a matrix is called an oblique projector. See Stewart (2011).

# Notes and References for 2.5

The computation of the CS decomposition is discussed in 8.7.6. For a discussion of its analytical properties, see:   
C. Davis and W. Kahan (1970). “The Rotation of Eigenvectors by a Perturbation III,” SIAM J. Numer. Anal. 7, 1–46.   
G.W. Stewart (1977). “On the Perturbation of Pseudo-Inverses, Projections and Linear Least Squares Problems,” SIAM Review 19, 634–662.   
C.C. Paige and M. Saunders (1981). “Toward a Generalized Singular Value Decomposition,” SIAM J. Numer. Anal. 18, 398–405.   
C.C. Paige and M. Wei (1994). “History and Generality of the CS Decomposition,” Lin. Alg. Applic. 208/209, 303–326.   
A detailed numerical discussion of oblique projectors (P2.5.5) is given in:   
G.W. Stewart (2011). “On the Numerical Analysis of Oblique Projectors,” SIAM J. Matrix Anal. Applic. 32, 309–348.

# 2.6 The Sensitivity of Square Systems

We use tools developed in previous sections to analyze the linear system problem $A x = b$ where $A \in \mathbb { R } ^ { n \times n }$ is nonsingular and $b \in \mathbb { R } ^ { n }$ . Our aim is to examine how perturbations in A and b affect the solution x. Higham (ASNA) offers a more detailed treatment.

# 2.6.1 An SVD Analysis

If

$$
A = \sum_ {i = 1} ^ {n} \sigma_ {i} u _ {i} v _ {i} ^ {T} = U \Sigma V ^ {T}
$$

is the SVD of A, then

$$
x = A ^ {- 1} b = (U \Sigma V ^ {T}) ^ {- 1} b = \sum_ {i = 1} ^ {n} \frac {u _ {i} ^ {T} b}{\sigma_ {i}} v _ {i}. \tag {2.6.1}
$$

This expansion shows that small changes in A or b can induce relatively large changes in x if $\sigma _ { n }$ is small.

It should come as no surprise that the magnitude of $\sigma _ { n }$ should have a bearing on the sensitivity of the $A x = b$ problem. Recall from Theorem 2.4.8 that $\sigma _ { n }$ is the 2-norm distance from A to the set of singular matrices. As the matrix of coefficients approaches this set, it is intuitively clear that the solution x should be increasingly sensitive to perturbations.

# 2.6.2 Condition

A precise measure of linear system sensitivity can be obtained by considering the parameterized system

$$
(A + \epsilon F) x (\epsilon) = b + \epsilon f, \qquad x (0) = x,
$$

where $F \in \mathbb { R } ^ { n \times n }$ and $f \in \mathbb { R } ^ { n }$ . If A is nonsingular, then it is clear that $x ( \epsilon )$ is differentiable in a neighborhood of zero. Moreover, $\dot { x } ( 0 ) \ = \ A ^ { - 1 } ( f - F x )$ and so the Taylor series expansion for $x ( \epsilon )$ has the form

$$
x (\epsilon) = x + \epsilon \dot {x} (0) + O (\epsilon^ {2}).
$$

Using any vector norm and consistent matrix norm we obtain

$$
\frac {\| x (\epsilon) - x \|}{\| x \|} \leq | \epsilon | \| A ^ {- 1} \| \left\{\frac {\| f \|}{\| x \|} + \| F \| \right\} + O \left(\epsilon^ {2}\right). \tag {2.6.2}
$$

For square matrices A define the condition number $\kappa ( A )$ by

$$
\kappa (A) = \left\| A \right\| \left\| A ^ {- 1} \right\| \tag {2.6.3}
$$

with the convention that $\kappa ( A ) = \infty$ for singular A. From $\left\| ~ b ~ \right\| ~ \leq ~ \left\| ~ A ~ \right\| ~ \left\| ~ x ~ \right\|$ and (2.6.2) it follows that

$$
\frac {\| x (\epsilon) - x \|}{\| x \|} \leq \kappa (A) (\rho_ {A} + \rho_ {b}) + O (\epsilon^ {2}) \tag {2.6.4}
$$

where

$$
\rho_ {A} = | \epsilon | \frac {\| F \|}{\| A \|} \text { and } \rho_ {b} = | \epsilon | \frac {\| f \|}{\| b \|}
$$

represent the relative errors in A and b, respectively. Thus, the relative error in x can be $\kappa ( A )$ times the relative error in A and b. In this sense, the condition number $\kappa ( A )$ quantifies the sensitivity of the $A x = b$ problem.

Note that $\kappa ( \cdot )$ depends on the underlying norm and subscripts are used accordingly, e.g.,

$$
\kappa_ {2} (A) = \| A \| _ {2} \| A ^ {- 1} \| _ {2} = \frac {\sigma_ {\max} (A)}{\sigma_ {\min} (A)}. \tag {2.6.5}
$$

Thus, the 2-norm condition of a matrix A measures the elongation of the hyperellipsoid $\{ A x : \| x \| _ { 2 } = 1 \}$ .

We mention two other characterizations of the condition number. For p-norm condition numbers, we have

$$
\frac {1}{\kappa_ {p} (A)} = \min _ {A + \Delta A \text {   singular }} \frac {\| \Delta A \| _ {p}}{\| A \| _ {p}}. \tag {2.6.6}
$$

This result may be found in Kahan (1966) and shows that $\kappa _ { p } ( A )$ measures the relative p-norm distance from A to the set of singular matrices.

For any norm, we also have

$$
\kappa (A) = \lim _ {\epsilon \rightarrow 0} \sup _ {\| \Delta A \| \leq \epsilon \| A \|} \frac {\| (A + \Delta A) ^ {- 1} - A ^ {- 1} \|}{\epsilon} \frac {1}{\| A ^ {- 1} \|}. \tag {2.6.7}
$$

This imposing result merely says that the condition number is a normalized Fr´echet derivative of the map $A \to A ^ { - 1 }$ . Further details may be found in Rice (1966). Recall that we were initially led to $\kappa ( A )$ through differentiation.

If $\kappa ( A )$ is large, then A is said to be an ill-conditioned matrix. Note that this is a norm-dependent property.1 However, any two condition numbers $\kappa _ { \alpha } ( \cdot )$ and $\kappa _ { \beta } ( \cdot )$ o n IRn×n $\mathbb { R } ^ { n \times n }$ are equivalent in that constants $c _ { 1 }$ and $c _ { 2 }$ can be found for which

$$
c _ {1} \kappa_ {\alpha} (A) \leq \kappa_ {\beta} (A) \leq c _ {2} \kappa_ {\alpha} (A), \quad A \in \mathbb {R} ^ {n \times n}.
$$

For example, on $\mathbb { R } ^ { n \times n }$ we have

$$
\frac {1}{n} \kappa_ {2} (A) \leq \kappa_ {1} (A) \leq n \kappa_ {2} (A),
$$

$$
\frac {1}{n} \kappa_ {\infty} (A) \leq \kappa_ {2} (A) \leq n \kappa_ {\infty} (A), \tag {2.6.8}
$$

$$
\frac {1}{n ^ {2}} \kappa_ {1} (A) \leq \kappa_ {\infty} (A) \leq n ^ {2} \kappa_ {1} (A).
$$

Thus, if a matrix is ill-conditioned in the α-norm, it is ill-conditioned in the $\beta \mathrm { { - n o r m } }$ modulo the constants $c _ { 1 }$ and $c _ { 2 }$ above.

For any of the p-norms, we have $\kappa _ { p } ( A ) \geq 1$ . Matrices with small condition numbers are said to be well-conditioned. In the 2-norm, orthogonal matrices are perfectly conditioned because if $Q$ is orthogonal, then $\kappa _ { 2 } ( Q ) = \| Q \| _ { 2 } \| Q ^ { T } \| _ { 2 } = 1$ .

# 2.6.3 Determinants and Nearness to Singularity

It is natural to consider how well determinant size measures ill-conditioning. If det $( A ) =$ 0 is equivalent to singularity, is det $( A ) \approx 0$ equivalent to near singularity? Unfortunately, there is little correlation between det(A) and the condition of $A x = b$ . For example, the matrix $B _ { n }$ defined by

$$
B _ {n} = \left[ \begin{array}{c c c c} 1 & - 1 & \dots & - 1 \\ 0 & 1 & \dots & - 1 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & 1 \end{array} \right] \in \mathbb {R} ^ {n \times n} \tag {2.6.9}
$$

has unit determinant, but $\kappa _ { \infty } ( B _ { n } ) ~ = ~ n \cdot 2 ^ { n - 1 }$ . On the other hand, a very wellconditioned matrix can have a very small determinant. For example,

$$
D _ {n} = \operatorname{diag} (1 0 ^ {- 1}, \dots , 1 0 ^ {- 1}) \in \mathbb {R} ^ {n \times n}
$$

satisfies $\kappa _ { p } ( D _ { n } ) = 1$ although det $( D _ { n } ) = 1 0 ^ { - n }$ .

# 2.6.4 A Rigorous Norm Bound

Recall that the derivation of (2.6.4) was valuable because it highlighted the connection between $\kappa ( A )$ and the rate of change of $x ( \epsilon )$ at $\epsilon = 0$ . However, it is a little unsatisfying because it is contingent on 
 being “small enough” and because it sheds no light on the size of the $O ( \epsilon ^ { 2 } )$ term. In this and the next subsection we develop some additional $A x = b$ perturbation theorems that are completely rigorous.

We first establish a lemma that indicates in terms of $\kappa ( A )$ when we can expect a perturbed system to be nonsingular.

Lemma 2.6.1. Suppose

$$
A x = b, \quad A \in \mathbb {R} ^ {n \times n}, 0 \neq b \in \mathbb {R} ^ {n},
$$

$$
(A + \Delta A) y = b + \Delta b, \quad \Delta A \in \mathbb {R} ^ {n \times n}, \Delta b \in \mathbb {R} ^ {n},
$$

with $\| \Delta A \| \le \epsilon \| A \|$ and $\| \Delta b \| \leq \epsilon \| b \|$ . $I f \epsilon \kappa ( A ) = r < 1$ , then $A + \Delta A$ is nonsingular and

$$
\frac {\parallel y \parallel}{\parallel x \parallel} \leq \frac {1 + r}{1 - r}.
$$

Proof. Since $\| A ^ { - 1 } \Delta A \| ~ \le ~ \epsilon \| ~ A ^ { - 1 } \| ~ \| ~ A \| = r < 1$ it follows from Theorem 2.3.4 that $( A + \Delta A )$ is nonsingular. Using Lemma 2.3.3 and the equality

$$
(I + A ^ {- 1} \Delta A) y = x + A ^ {- 1} \Delta b
$$

we find

$$
\begin{array}{l} \| y \| \leq \| (I + A ^ {- 1} \Delta A) ^ {- 1} \| \left(\| x \| + \epsilon \| A ^ {- 1} \| \| b \|\right) \\ \leq \frac {1}{1 - r} \left(\| x \| + \epsilon \| A ^ {- 1} \| \| b \|\right) = \frac {1}{1 - r} \left(\| x \| + r \frac {\| b \|}{\| A \|}\right). \\ \end{array}
$$

Since $\parallel b \parallel = \parallel A x \parallel \leq \parallel A \parallel \parallel x \parallel$ it follows that

$$
\parallel y \parallel \leq \frac {1}{1 - r} \left(\parallel x \parallel + r \parallel x \parallel\right)
$$

and this establishes the required inequality.

We are now set to establish a rigorous $A x = b$ perturbation bound.

Theorem 2.6.2. If the conditions of Lemma 2.6.1 hold, then

$$
\frac {\| y - x \|}{\| x \|} \leq \frac {2 \epsilon}{1 - r} \kappa (A). \tag {2.6.10}
$$

Proof. Since

$$
y - x = A ^ {- 1} \Delta b - A ^ {- 1} \Delta A y \tag {2.6.11}
$$

we have

$$
\| y - x \| \leq \epsilon \| A ^ {- 1} \| \| b \| + \epsilon \| A ^ {- 1} \| \| A \| \| y \|.
$$

Thus,

$$
\frac {\parallel y - x \parallel}{\parallel x \parallel} \leq \epsilon \kappa (A) \frac {\parallel b \parallel}{\parallel A \parallel \parallel x \parallel} + \epsilon \kappa (A) \frac {\parallel y \parallel}{\parallel x \parallel} \leq \epsilon \left(1 + \frac {1 + r}{1 - r}\right) \kappa (A),
$$

from which the theorem readily follows.

A small example helps put this result in perspective. The $A x = b$ problem

$$
{\left[ \begin{array}{c c} 1 & 0 \\ 0 & 1 0 ^ {- 6} \end{array} \right]} {\left[ \begin{array}{c} x _ {1} \\ x _ {2} \end{array} \right]} = {\left[ \begin{array}{c} 1 \\ 1 0 ^ {- 6} \end{array} \right]}
$$

has solution $x = [ 1 , 1 ] ^ { T }$ and condition $\kappa _ { \infty } ( A ) = 1 0 ^ { 6 }$ . If $\Delta { b } = [ 1 0 ^ { - 6 } , 0 ] ^ { T } , \Delta { A } = 0$ , and $( A + \Delta A ) y = b + \Delta b$ , then $y = [ 1 + 1 0 ^ { - 6 } , 1 ] ^ { T }$ and the inequality (2.6.10) says

$$
1 0 ^ {- 6} = \frac {\parallel x - y \parallel_ {\infty}}{\parallel x \parallel_ {\infty}} \ll \frac {\parallel \Delta b \parallel_ {\infty}}{\parallel b \parallel_ {\infty}} \kappa_ {\infty} (A) = 1 0 ^ {- 6} 1 0 ^ {6} = 1.
$$

Thus, the upper bound in (2.6.10) can be a gross overestimate of the error induced by the perturbation.

On the other hand, if $\Delta b = ( 0 , 1 0 ^ { - 6 } ) ^ { T } , \Delta A = 0$ , and $( A + \Delta A ) y = b + \Delta b$ , then this inequality says that

$$
\frac {1 0 ^ {0}}{1 0 ^ {0}} \leq 2 \times 1 0 ^ {- 6} 1 0 ^ {6}.
$$

Thus, there are perturbations for which the bound in (2.6.10) is essentially attained.

# 2.6.5 More Refined Bounds

An interesting refinement of Theorem 2.6.2 results if we extend the notion of absolute value to matrices:

$$
F = (f _ {i j}) \in \mathbb {R} ^ {m \times n} \qquad \Rightarrow \qquad | F | = (| f _ {i j} |) \in \mathbb {R} ^ {m \times n}.
$$

This notation together with a matrix-level version of $\stackrel { 6 6 } { \leq } \underline { { \stackrel { \triangledown } { \scriptscriptstyle \mathrm { \blacktriangledown } } } }$ makes it easy to specify componentwise error bounds. If $F , G \in \mathbb { R } ^ { m \times n }$ , then

$$
| F | \leq | G | \quad \Leftrightarrow \quad | f _ {i j} | \leq | g _ {i j} |
$$

for all i and $j .$ Also note that if $F \in \mathbb { R } ^ { m \times q }$ and $G \in \mathbb { R } ^ { q \times n }$ , then $| F G | \leq | F | \cdot | G |$ . With these definitions and facts we obtain the following refinement of Theorem 2.6.2.

Theorem 2.6.3. Suppose

$$
A x = b, \quad A \in \mathbb {R} ^ {n \times n}, 0 \neq b \in \mathbb {R} ^ {n},
$$

$$
(A + \Delta A) y, = b + \Delta b \Delta A \in \mathbb {R} ^ {n \times n}, \Delta b \in \mathbb {R} ^ {n},
$$

and that $| \Delta A | \le \epsilon | A |$ and $| \Delta b | \le \epsilon | b | . \ I f \delta \kappa _ { \infty } ( A ) = r < 1$ , then $( A { + } \Delta A )$ is nonsingular and

$$
\frac {\| y - x \| _ {\infty}}{\| x \| _ {\infty}} \leq \frac {2 \epsilon}{1 - r} \cdot \| | A ^ {- 1} | | A | \| _ {\infty}. \tag {2.6.12}
$$

Proof. Since $\| \Delta A \| _ { \infty } \leq \epsilon \| A \| _ { \infty }$ and $\| \Delta b \| _ { \infty } \leq \epsilon \| b \| _ { \infty }$ the conditions of Lemma 2.6.1 are satisfied in the infinity norm. This implies that $A + \Delta A$ is nonsingular and

$$
\frac {\parallel y \parallel_ {\infty}}{\parallel x \parallel_ {\infty}} \leq \frac {1 + r}{1 - r}.
$$

Now using (2.6.11) we find

$$
\begin{array}{l} | y - x | \leq | A ^ {- 1} | | \Delta b | + | A ^ {- 1} | | \Delta A | | y | \\ \leq \epsilon | A ^ {- 1} | | b | + \epsilon | A ^ {- 1} | | A | | y | \leq \epsilon | A ^ {- 1} | | A | (| x | + | y |). \\ \end{array}
$$

If we take norms, then

$$
\| y - x \| _ {\infty} \leq \epsilon \| | A ^ {- 1} | | A | \| _ {\infty} \left(\| x \| _ {\infty} + \frac {1 + r}{1 - r} \| x \| _ {\infty}\right).
$$

The theorem follows upon division by $\| { \boldsymbol { x } } \| _ { \infty }$ .

The quantity $\parallel | A ^ { - 1 } | | A | \parallel _ { \infty }$ is known as the Skeel condition number and there are examples where it is considerably less than $\kappa _ { \infty } ( A )$ . In these situations, (2.6.12) is more informative than (2.6.10).

Norm bounds are frequently good enough when assessing error, but sometimes it is desirable to examine error at the component level. Oettli and Prager (1964) have an interesting result that indicates if an approximate solution $\hat { x } \in \mathbb { R } ^ { n }$ to the n-by-n system $A x = b$ satisfies a perturbed system with prescribed structure. Consider the problem of finding $\Delta A \in \mathbb { R } ^ { n \times n } , \Delta b \in \mathbb { R } ^ { n }$ , and $\omega \ge 0$ such that

$$
(A + \Delta A) \hat {x} = b + \Delta b \quad | \Delta A | \leq \omega | E |, | \Delta b | \leq \omega | f |. \tag {2.6.13}
$$

where $E \in \mathbb { R } ^ { n \times n }$ and $f \in \mathbb { R } ^ { n }$ are given. With proper choice of E and f, the perturbed system can take on certain qualities. For example, if $E = A$ and $f = b$ and ω is small, then ˆx satisfies a nearby system in the componentwise sense. The authors show that for a given $A , b , { \hat { x } } , E ,$ , and $f$ the smallest ω possible in (2.6.13) is given by

$$
\omega_ {\mathrm{min}} = \max _ {1 \leq i \leq n} \frac {| A \hat {x} - b | _ {i}}{(| E | \cdot | \hat {x} | + | f |) _ {i}}.
$$

If $A \hat { x } = b$ , then $\omega _ { \mathrm { m i n } } = 0$ . On the other hand, if $\omega _ { \mathrm { m i n } } = \infty$ , then ˆx does not satisfy any system of the prescribed perturbation structure.

# Problems

P2.6.1 Show that if $\parallel I \parallel \geq 1$ , then $\kappa ( A ) \geq 1$ .

P2.6.2 Show that for a given norm, $\kappa ( A B ) \leq \kappa ( A ) \kappa ( B )$ and that $\kappa ( \alpha A ) = \kappa ( A )$ for all nonzero α.

P2.6.3 Relate the 2-norm condition of $X \in \mathbb { R } ^ { m \times n } ( m \geq n )$ to the 2-norm condition of the matrices

$$
B = \left[ \begin{array}{c c} I _ {m} & X \\ 0 & I _ {n} \end{array} \right] \qquad \mathrm{and} \qquad C = \left[ \begin{array}{c} X \\ I _ {n} \end{array} \right].
$$

P2.6.4 Suppose $A \in \mathbb { R } ^ { n \times n }$ is nonsingular. Assume for a particular i and j that there is no way to make A singular by changing the value of $a _ { i j }$ . What can you conclude about $A ^ { - 1 \mathord { \ ? } }$ Hint: Use the Sherman-Morrison formula.

P2.6.5 Suppose $A \in \mathbb { R } ^ { n \times n }$ is nonsingular, $b \in \mathbb { R } ^ { n } , A x = b$ , and $C = A ^ { - 1 }$ . Use the Sherman-Morrison formula to show that

$$
\frac {\partial x _ {k}}{\partial a _ {i j}} = - x _ {j} c _ {k i}.
$$

# Notes and References for §2.6

The condition concept is thoroughly investigated in:

J. Rice (1966). “A Theory of Condition,” SIAM J. Numer. Anal. 3, 287–310.

W. Kahan (1966). “Numerical Linear Algebra,” Canadian Math. Bull. 9, 757–801.

References for componentwise perturbation theory include:

W. Oettli and W. Prager (1964). “Compatibility of Approximate Solutions of Linear Equations with Given Error Bounds for Coefficients and Right Hand Sides,” Numer. Math. 6, 405–409.

J.E. Cope and B.W. Rust (1979). “Bounds on Solutions of Systems with Accurate Data,” SIAM J. Numer. Anal. l6, 950–63.

R.D. Skeel (1979). “Scaling for Numerical Stability in Gaussian Elimination,” J. ACM 26, 494–526.

J.W. Demmel (1992). “The Componentwise Distance to the Nearest Singular Matrix,” SIAM J. Matrix Anal. Applic. 13, 10–19.

D.J. Higham and N.J. Higham (1992). “Componentwise Perturbation Theory for Linear Systems with Multiple Right-Hand Sides,” Lin. Alg. Applic. 174, 111–129.

N.J. Higham (1994). “A Survey of Componentwise Perturbation Theory in Numerical Linear Algebra,” in Mathematics of Computation 1943–1993: A Half Century of Computational Mathematics, W. Gautschi (ed.), Volume 48 of Proceedings of Symposia in Applied Mathematics, American Mathematical Society, Providence, RI.

S. Chandrasekaren and I.C.F. Ipsen (1995). “On the Sensitivity of Solution Components in Linear Systems of Equations,” SIAM J. Matrix Anal. Applic. 16, 93–112.

S.M. Rump (1999). “Ill-Conditioned Matrices Are Componentwise Near to Singularity,” SIAM Review 41, 102–112.

The reciprocal of the condition number measures how near a given Ax = b problem is to singularity. The importance of knowing how near is a given problem to a difficult or insoluble problem has come to be appreciated in many computational settings, see:

A. Laub(1985). “Numerical Linear Algebra Aspects of Control Design Computations,” IEEE Trans. Autom. Control. AC-30, 97–108.

J.W. Demmel (1987). “On the Distance to the Nearest Ill-Posed Problem,” Numer. Math. 51, 251–289.

N.J. Higham (1989). “Matrix Nearness Problems and Applications,” in Applications of Matrix Theory, M.J.C. Gover and S. Barnett (eds.), Oxford University Press, Oxford, UK, 1–27.

Much has been written about problem sensitivity from the statistical point of view, see:

J.W. Demmel (1988). “The Probability that a Numerical Analysis Problem is Difficult,” Math. Comput. 50, 449–480.

G.W. Stewart (1990). “Stochastic Perturbation Theory,” SIAM Review 32, 579–610.

C. S. Kenney, A.J. Laub, and M.S. Reese (1998). “Statistical Condition Estimation for Linear Systems,” SIAM J. Sci. Comput. 19, 566–583.

The problem of minimizing κ2(A + U V T ) where U V T is a low-rank matrix is discussed in:

C. Greif and J.M. Varah (2006). “Minimizing the Condition Number for Small Rank Modifications,” SIAM J. Matrix Anal. Applic. 29, 82–97.

# 2.7 Finite Precision Matrix Computations

Rounding errors are part of what makes the field of matrix computations so challenging. In this section we describe a model of floating point arithmetic and then use it to develop error bounds for floating point dot products, saxpys, matrix-vector products, and matrix-matrix products.

# 2.7.1 A 3-digit Calculator

Suppose we have a base-10 calculator that represents nonzero numbers in the following style:

$$
x = \pm d _ {0}. d _ {1} d _ {2} \times 1 0 ^ {e} \qquad \text {where} \quad \left\{ \begin{array}{r l} 1 & \leq d _ {0} \leq 9, \\ 0 & \leq d _ {1} \leq 9, \\ 0 & \leq d _ {2} \leq 9, \\ - 9 & \leq e \leq 9. \end{array} \right.
$$

Let us call these numbers floating point numbers. After playing around a bit we make a number of observations:

• The precision of the calculator has to do with the “length” of the significand $d _ { 0 } . d _ { 1 } d _ { 2 }$ . For example, the number π would be represented as $3 . 1 4 \times 1 0 ^ { 0 }$ , which has a relative error approximately equal to $1 0 ^ { - 3 }$ .

• There is not enough “room” to store exactly the results from most arithmetic operations between floating point numbers. Sums and products like

$$
(1. 2 3 \times 1 0 ^ {6}) + (4. 5 6 \times 1 0 ^ {4}) = 1 2 7 5 6 0 0,
$$

$$
(1. 2 3 \times 1 0 ^ {1}) * (4. 5 6 \times 1 0 ^ {2}) = 5 6 0 8. 8
$$

involve more than three significant digits. Results must be rounded in order to $\mathrm { \Omega ^ { \circ } f i t { \Omega ^ { \circ } } }$ the 3-digit format, e.g., round $( 1 2 7 5 6 0 0 ) = 1 . 2 8 \times 1 0 ^ { 6 }$ , round(5608.8) = $5 . 6 1 \times 1 0 ^ { 3 }$ .

• If zero is to be a floating point number (and it must be), then we need a special convention for its representation, $\mathrm { e . g . , 0 . 0 0 \times 1 0 ^ { 0 } }$ .

• In contrast to the real numbers, there is a smallest positive floating point number $( N _ { \mathrm { m i n } } = 1 . 0 0 { \times } 1 0 ^ { - 9 } )$ and there is a largest positive floating point number $( N _ { \mathrm { m a x } } =$ $9 . 9 9 \times 1 0 ^ { 9 } )$ .

• Some operations yield answers whose exponents exceed the 1-digit allocation, $\mathrm { e . g . , ( 1 . 2 3 \times 1 0 ^ { 4 } ) * ( 4 . 5 6 \times 1 0 ^ { 7 } ) }$ and $( 1 . 2 3 \times 1 0 ^ { - 2 } ) / ( 4 . 5 6 \times 1 0 ^ { 8 } )$ .

• The set of floating point numbers is finite. For the toy calculator there are $2 \times 9 \times 1 0 \times 1 0 \times 1 9 + 1 = 3 4 2 0 1$ floating point numbers.

• The spacing between the floating point numbers varies. Between $1 . 0 0 \times 1 0 ^ { e }$ and $1 . 0 0 \times \bar { 1 0 ^ { e + 1 } }$ the spacing is $1 0 ^ { e - \hat { 2 } }$ .

The careful design and analysis of a floating point computation requires an understanding of these inexactitudes and limitations. How are results rounded? How accurate is floating point arithmetic? What can we say about a sequence of floating point operations?

# 2.7.2 IEEE Floating Point Arithmetic

To build a solid, practical understanding of finite precision computation, we set aside our toy, motivational base-10 calculator and consider the key ideas behind the widely accepted IEEE floating point standard. The IEEE standard includes a 32-bit single format and a 64-bit double format. We will illustrate concepts using the latter as an example because typical accuracy requirements make it the format of choice.

The importance of having a standard for floating point arithmetic that is upheld by hardware manufacturers cannot be overstated. After all, floating point arithmetic is the foundation upon which all of scientific computing rests. The IEEE standard promotes software reliability and enables numerical analysts to make rigorous statements about computed results. Our discussion is based on the excellent book by Overton (2001).

The 64-bit double format allocates a single bit for the sign of the floating point number, 52 bits for the mantissa , and eleven bits for the exponent:

$$
x: \boxed {\pm} \left| a _ {1} a _ {2} \dots a _ {1 1} \right| b _ {1} b _ {2} \dots b _ {5 2}. \tag {2.7.1}
$$

The “formula” for the value of this representation depends upon the exponent bits:

If $a _ { 1 } \ldots a _ { 1 1 }$ is neither all 0’s nor all 1’s, then x is a normalized floating point number with value

$$
x = \pm (1. b _ {1} b _ {2} \dots b _ {5 2}) _ {2} \times 2 ^ {(a _ {1} a _ {2} \dots a _ {1 1}) _ {2} - 1 0 2 3}. \tag {2.7.2}
$$

The “1023 bias” in the exponent supports the graceful inclusion of various “unnormalized” floating numbers which we describe shortly. Several important quantities capture the finiteness of the representation. The machine epsilon is the gap between 1 and the next largest floating point number. Its value is $2 ^ { - 5 2 } \approx 1 0 ^ { - 1 6 }$ for the double format. Among the positive normalized floating point numbers, $N _ { \mathrm { m i n } } = 2 ^ { - 1 0 2 2 } \approx 1 0 ^ { - 3 0 8 }$ is the smallest and $N _ { \mathrm { { m a x } } } = ( 2 - 2 ^ { - 5 2 } ) 2 ^ { 1 0 2 3 } \approx 1 0 ^ { 3 0 8 }$ is the largest. A real number x is within the normalized range if $N _ { \mathrm { m i n } } \le | x | \le N _ { \mathrm { m a x } }$ .

If $a _ { 1 } \ldots a _ { 1 1 }$ is all $0 \mathrm { { s } }$ , then the value of the representation (2.7.1) is

$$
x = \pm \left(0. b _ {1} b _ {2} \dots b _ {5 2}\right) _ {2} \times 2 ^ {\left(a _ {1} a _ {2} \dots a _ {1 1}\right) _ {2} - 1 0 2 2} \tag {2.7.3}
$$

This includes 0 and the subnormal floating point numbers. This feature creates a uniform spacing of the floating point numbers between $- N _ { \mathrm { m i n } }$ and $+ N _ { \mathrm { m i n } }$ .

If $a _ { 1 } \ldots a _ { 1 1 }$ is all 1’s, then the encoding (2.7.1) represents inf for +∞, -inf for $- \infty$ , or NaN for “not-a-number.” The determining factor is the value of the $b _ { i }$ . (If the $b _ { i }$ are not all zero, then the value of x is NaN.) Quotients like $1 / 0 , - 1 / 0$ , and $0 / 0$ produce these special floating point numbers instead of prompting program termination.

There are four rounding modes: round down (toward $- \infty )$ , round up (toward $+ \infty )$ , round-toward-zero, and round-toward-nearest. We focus on round-toward-nearest since it is the mode almost always used in practice.

If a real number x is outside the range of the normalized floating point numbers then

$$
\operatorname{round} (x) = \left\{ \begin{array}{l l} - \infty & \text {if} x <   - N _ {\max} , \\ + \infty & \text {if} x > N _ {\max}. \end{array} \right.
$$

Otherwise, the rounding process depends upon its floating point “neighbors”:

$x _ { - }$ is the nearest floating point number to x that is $\leq x .$

$x _ { + }$ is the nearest floating point number to x that is $\geq x$

Define $d _ { - } = x - x _ { - }$ and $d _ { + } = x _ { + } - x$ and let “lsb” stand for “least significant bit.” If $N _ { \mathrm { m i n } } \le | x | \le N _ { \mathrm { m a x } }$ , then

$$
\operatorname{round} (x) = \left\{ \begin{array}{l l} x _ {-} & \text { if } d _ {-} <   d _ {+} \text { or } d _ {-} = d _ {+} \text { and } \operatorname{lsb} (x _ {-}) = 0, \\ x _ {+} & \text { if } d _ {+} <   d _ {-} \text { or } d _ {+} = d _ {-} \text { and } \operatorname{lsb} (x _ {+}) = 0. \end{array} \right.
$$

The tie-breaking criteria is well-defined because $x _ { - }$ and $x _ { + }$ are adjacent floating point numbers and so must differ in their least significant bit.

Regarding the accuracy of the round-to-nearest strategy, suppose x is a real number that satisfies $N _ { \mathrm { m i n } } \le | x | \le N _ { \mathrm { m a x } }$ . Thus,

$$
| \operatorname{round} (x) - x | \leq \frac {2 ^ {- 5 2}}{2} 2 ^ {e} \leq \frac {2 ^ {- 5 2}}{2} | x |
$$

which says that relative error is bounded by half of the machine epsilon:

$$
\frac {| \operatorname{round} (x) - x |}{| x |} \leq 2 ^ {- 5 3}.
$$

The IEEE standard stipulates that each arithmetic operation be correctly rounded, meaning that the computed result is the rounded version of the exact result. The implementation of correct rounding is far from trivial and requires registers that are equipped with several extra bits of precision.

We mention that the IEEE standard also requires correct rounding in the square root operation, the remainder operation, and various format conversion operations.

# 2.7.3 The “fl” Notation

With intuition gleaned from the toy calculator example and an understanding of IEEE arithmetic, we are ready to move on to the roundoff analysis of some basic algebraic calculations. The challenge when presenting the effects of finite precision arithmetic in this section and throughout the book is to communicate essential behavior without excessive detail. To that end we use the notation fl(·) to identify a floating point storage and/or computation. Unless exceptions are a critical part of the picture, we freely invoke the fl notation without mentioning “−∞,” “∞,” “NaN,” etc.

If $x \in \mathbb { R }$ , then fl(x) is its floating point representation and we assume that

$$
\mathsf {f l} (x) = x (1 + \delta), \quad | \delta | \leq \mathbf {u}, \tag {2.7.4}
$$

where u is the unit roundoff defined by

$$
\mathbf {u} = \frac {1}{2} \times (\text { gap   between   1   and   next   largest   floating   point   number }). \tag {2.7.5}
$$

The unit roundoff for IEEE single format is about $1 0 ^ { - 7 }$ and for double format it is about 10−16.

If x and y are floating point numbers and “op” is any of the four arithmetic operations, then fl(x op y) is the floating point result from the floating point op. Following Trefethen and Bau (NLA), the fundamental axiom of floating point arithmetic is that

$$
\operatorname{fl} (x \text {   op   } y) = (x \text {   op   } y) (1 + \delta), \quad | \delta | \leq \mathbf {u}, \tag {2.7.6}
$$

where x and y are floating point numbers and the “op” inside the fl operation means “floating point operation.” This shows that there is small relative error associated with individual arithmetic operations:

$$
\frac {| \mathsf {f l} (x \textsf {o p} y) - (x \textsf {o p} y) |}{| x \textsf {o p} y |} \leq \mathbf {u}, \qquad x \textsf {o p} y \neq 0.
$$

Again, unless it is particularly relevant to the discussion, it will be our habit not to bring up the possibilities of an exception arising during the floating point operation.

# 2.7.4 Become a Floating Point Thinker

It is a good idea to have a healthy respect for the subleties of floating point calculation. So before we proceed with our first serious roundoff error analysis we offer three maxims to keep in mind when designing a practical matrix computation. Each reinforces the distinction between computer arithmetic and exact arithmetic.

# Maxim 1. Order is Important.

Floating point arithmetic is not associative. For example, suppose

$$
x = 1. 2 4 \times 1 0 ^ {0}, \qquad y = - 1. 2 3 \times 1 0 ^ {0}, \qquad z = 1. 0 0 \times 1 0 ^ {- 3}.
$$

Using toy calculator arithmetic we have

$$
\mathsf {f l} (\mathsf {f l} (x + y) + z)) = 1. 1 0 \times 1 0 ^ {- 2}
$$

while

$$
\operatorname{fl} (x + \operatorname{fl} (y + z)) = 1. 0 0 \times 1 0 ^ {- 2}.
$$

A consequence of this is that mathematically equivalent algorithms may produce different results in floating point.

# Maxim 2. Larger May Mean Smaller.

Suppose we want to compute the derivative of $f ( x ) = \sin ( x )$ using a divided difference. Calculus tells us that $d = ( \sin ( x + h ) - \sin ( x ) ) / h$ satisfies $| d - \cos ( x ) | = O ( h )$ which argues for making h as small as possible. On the other hand, any roundoff error sustained in the sine evaluations is magnified by $1 / h$ . By setting $h = { \sqrt { \mathbf { u } } } .$ , the sum of the calculus error and roundoff error is approximately minimized. In other words, a value of h much greater than u renders a much smaller overall error. See Overton(2001, pp. 70–72).

# Maxim 3. A Math Book Is Not Enough.

The explicit coding of a textbook formula is not always the best way to design an effective computation. As an example, we consider the quadratic equation $x ^ { 2 } - 2 p x - q =$ 0 where both p and q are positive. Here are two methods for computing the smaller (necessarily real) root:

$\mathrm { M e t h o d ~ 1 : } \quad r _ { \mathrm { m i n } } = p - \sqrt { p ^ { 2 } + q } ,$

$\mathrm { M e t h o d ~ 2 } { : } \quad r _ { \mathrm { m i n } } = \frac { q } { p + \sqrt { p ^ { 2 } + q } } .$

The first method is based on the familiar quadratic formula while the second uses the fact that −q is the product of $r _ { \mathrm { m i n } }$ and the larger root. Using IEEE double format arithmetic with input $p = 1 2 3 4 5 6 7 8$ and $q = 1$ we obtain these results:

$\mathrm { M e t h o d ~ 1 } { : } \quad r _ { \mathrm { m i n } } = - 4 . 0 9 7 8 1 9 3 2 8 3 0 8 1 0 6 \times 1 0 ^ { - 8 } ,$

$\mathrm { M e t h o d ~ 2 } { : } \quad r _ { \mathrm { m i n } } = - 4 . 0 5 0 0 0 0 0 3 3 2 1 0 0 0 2 1 \times 1 0 ^ { - 8 } \quad ( \mathrm { c o r r e c t } ) .$

Method 1 produces an answer that has almost no correct significant digits. It attempts to compute a small number by subtracting a pair of nearly equal large numbers. Almost all correct significant digits in the input data are lost during the subtraction, a phenomenon known as catastrophic cancellation. In contrast, Method 2 produces an answer that is correct to full machine precision. It computes a small number as a division of one number by a much larger number. See Forsythe (1970).

Keeping these maxims in mind does not guarantee the production of accurate, reliable software, but it helps.

# 2.7.5 Application: Storing a Real Matrix

Suppose $A \in \mathbb { R } ^ { m \times n }$ and that we wish to quantify the errors associated with its floating point representation. Denoting the stored version of A by ${ \mathsf { f l } } ( A )$ , we see that

$$
[ \mathsf {f l} (A) ] _ {i j} = \mathsf {f l} (a _ {i j}) = a _ {i j} (1 + \epsilon_ {i j}), \quad | \epsilon_ {i j} | \leq \mathbf {u}, \tag {2.7.7}
$$

for all i and j, i.e.,

$$
| \mathrm{fl} (A) - A | \leq \mathbf {u} | A |.
$$

A relation such as this can be easily turned into a norm inequality, e.g.,

$$
\| \mathsf {f l} (A) - A \| _ {1} \leq \mathbf {u} \| A \| _ {1}.
$$

However, when quantifying the rounding errors in a matrix manipulation, the absolute value notation is sometimes more informative because it provides a comment on each entry.

# 2.7.6 Roundoff in Dot Products

We begin our study of finite precision matrix computations by considering the rounding errors that result in the standard dot product algorithm:

$$
s = 0
$$

for $k = 1 { : } n$

$$
s = s + x _ {k} y _ {k} \tag {2.7.8}
$$

end

Here, x and y are n-by-1 floating point vectors.

In trying to quantify the rounding errors in this algorithm, we are immediately confronted with a notational problem: the distinction between computed and exact quantities. If the underlying computations are clear, we shall use the fl(·) operator to signify computed quantities. Thus, $\mathsf { I } ( x ^ { T } y )$ denotes the computed output of (2.7.8). Let us bound $| \mathsf { f l } ( x ^ { \hat { T } } y ) - x ^ { T } y |$ . If

$$
s _ {p} = \mathsf {f l} \left(\sum_ {k = 1} ^ {p} x _ {k} y _ {k}\right),
$$

then $s _ { 1 } = x _ { 1 } y _ { 1 } ( 1 + \delta _ { 1 } )$ with $| \delta _ { 1 } | \leq \mathbf { u }$ and for $p = 2 { : } n$

$$
\begin{array}{l} s _ {p} = \mathsf {f l} (s _ {p - 1} + \mathsf {f l} (x _ {p} y _ {p})) \\ = (s _ {p - 1} + x _ {p} y _ {p} (1 + \delta_ {p})) (1 + \epsilon_ {p}) \quad | \delta_ {p} |, | \epsilon_ {p} | \leq \mathbf {u}. \tag {2.7.9} \\ \end{array}
$$

A little algebra shows that

$$
\mathsf {f l} (x ^ {T} y) = s _ {n} = \sum_ {k = 1} ^ {n} x _ {k} y _ {k} (1 + \gamma_ {k})
$$

where

$$
(1 + \gamma_ {k}) = (1 + \delta_ {k}) \prod_ {j = k} ^ {n} (1 + \epsilon_ {j})
$$

with the convention that $\epsilon _ { 1 } = 0$ . Thus,

$$
| \mathfrak {f l} (x ^ {T} y) - x ^ {T} y | \leq \sum_ {k = 1} ^ {n} | x _ {k} y _ {k} | | \gamma_ {k} |. \tag {2.7.10}
$$

To proceed further, we must bound the quantities $| \gamma _ { k } |$ in terms of u. The following result is useful for this purpose.

Lemma 2.7.1. $I f \left( 1 + \alpha \right) = \prod _ { k = 1 } ^ { n } ( 1 + \alpha _ { k } )$ where $| \alpha _ { k } | \le \mathbf { u }$ and $n \mathbf { u } \leq . 0 1$ , then $| \alpha | \leq$ 1.01nu. k=1

Proof. See Higham (ASNA, p. 75).

Application of this result to (2.7.10) under the “reasonable” assumption $n \mathbf { u } \leq . 0 1$ gives

$$
\left| \mathrm{fl} \left(x ^ {T} y\right) - x ^ {T} y \right| \leq 1. 0 1 n \mathbf {u} | x | ^ {T} | y |. \tag {2.7.11}
$$

Notice that if $| x ^ { T } y | \ll | x | ^ { T } | y |$ , then the relative error in $\mathsf { f l } ( x ^ { T } y )$ may not be small.

# 2.7.7 Alternative Ways to Quantify Roundoff Error

An easier but less rigorous way of bounding α in Lemma 2.7.1 is to say $| \alpha | \le n \mathbf { u } \mathbf { + } O ( \mathbf { u } ^ { 2 } )$ . With this convention we have

$$
\left| \mathrm{fl} \left(x ^ {T} y\right) - x ^ {T} y \right| \leq n \mathbf {u} \left| x \right| ^ {T} \left| y \right| + O \left(\mathbf {u} ^ {2}\right). \tag {2.7.12}
$$

Other ways of expressing the same result include

$$
\left| \mathbf {f l} \left(x ^ {T} y\right) - x ^ {T} y \right| \leq \phi (n) \mathbf {u} \left| x \right| ^ {T} | y | \tag {2.7.13}
$$

and

$$
\left| \mathbf {f l} \left(x ^ {T} y\right) - x ^ {T} y \right| \leq c n \mathbf {u} | x | ^ {T} | y |, \tag {2.7.14}
$$

where $\phi ( n )$ is a “modest” function of n and c is a constant of order unity.

We shall not express a preference for any of the error bounding styles shown in (2.7.11)–(2.7.14). This spares us the necessity of translating the roundoff results that appear in the literature into a fixed format. Moreover, paying overly close attention to the details of an error bound is inconsistent with the “philosophy” of roundoff analysis. As Wilkinson (1971, p. 567) says,

There is still a tendency to attach too much importance to the precise error bounds obtained by an a priori error analysis. In my opinion, the bound itself is usually the least important part of it. The main object of such an analysis is to expose the potential instabilities, if any, of an algorithm so that hopefully from the insight thus obtained one might be led to improved algorithms. Usually the bound itself is weaker than it might have been because of the necessity of restricting the mass of detail to a reasonable level and because of the limitations imposed by expressing the errors in terms of matrix norms. A priori bounds are not, in general, quantities that should be used in practice. Practical error bounds should usually be determined by some form of a posteriori error analysis, since this takes full advantage of the statistical distribution of rounding errors and of any special features, such as sparseness, in the matrix.

It is important to keep these perspectives in mind.

# 2.7.8 Roundoff in Other Basic Matrix Computations

It is easy to show that if A and B are floating point matrices and α is a floating point number, then

$$
\mathsf {f l} (\alpha A) = \alpha A + E, \quad | E | \leq \mathbf {u} | \alpha A |, \tag {2.7.15}
$$

and

$$
\mathsf {f l} (A + B) = (A + B) + E, \quad | E | \leq \mathbf {u} | A + B |. \tag {2.7.16}
$$

As a consequence of these two results, it is easy to verify that computed saxpy’s and outer product updates satisfy

$$
\operatorname{fl} (y + \alpha x) = y + \alpha x + z, \quad | z | \leq \mathbf {u} (| y | + 2 | \alpha x |) + O \left(\mathbf {u} ^ {2}\right), \tag {2.7.17}
$$

$$
\mathsf {f l} (C + u v ^ {T}) = C + u v ^ {T} + E, \quad | E | \leq \mathbf {u} \left(| C | + 2 | u v ^ {T} |\right) + O \left(\mathbf {u} ^ {2}\right). \tag {2.7.18}
$$

Using (2.7.11) it is easy to show that a dot-product-based multiplication of two floating point matrices A and B satisfies

$$
\mathsf {f l} (A B) = A B + E, \quad | E | \leq n \mathbf {u} | A | | B | + O \left(\mathbf {u} ^ {2}\right). \tag {2.7.19}
$$

The same result applies if a gaxpy or outer product based procedure is used. Notice that matrix multiplication does not necessarily give small relative error since $| A B |$ may be much smaller than $| A | | B |$ , e.g.,

$$
\left[ \begin{array}{c c} 1 & 1 \\ 0 & 0 \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ -. 9 9 & 0 \end{array} \right] = \left[ \begin{array}{c c}. 0 1 & 0 \\ 0 & 0 \end{array} \right].
$$

It is easy to obtain norm bounds from the roundoff results developed thus far. If we look at the 1-norm error in floating point matrix multiplication, then it is easy to show from (2.7.19) that

$$
\left\| \mathfrak {f l} (A B) - A B \right\| _ {1} \leq n \mathbf {u} \| A \| _ {1} \| B \| _ {1} + O \left(\mathbf {u} ^ {2}\right). \tag {2.7.20}
$$

# 2.7.9 Forward and Backward Error Analyses

Each roundoff bound given above is the consequence of a forward error analysis. An alternative style of characterizing the roundoff errors in an algorithm is accomplished through a technique known as backward error analysis. Here, the rounding errors are related to the input data rather than the answer. By way of illustration, consider the n = 2 version of triangular matrix multiplication. It can be shown that:

$$
\mathsf {f l} (A B) = \left[ \begin{array}{c c} a _ {1 1} b _ {1 1} (1 + \epsilon_ {1}) & (a _ {1 1} b _ {1 2} (1 + \epsilon_ {2}) + a _ {1 2} b _ {2 2} (1 + \epsilon_ {3})) (1 + \epsilon_ {4}) \\ 0 & a _ {2 2} b _ {2 2} (1 + \epsilon_ {5}) \end{array} \right]
$$

where $| \epsilon _ { i } | \le \mathbf { u } .$ , for $i = 1 { : } 5$ . However, if we define

$$
\hat {A} = \left[ \begin{array}{c c} a _ {1 1} & a _ {1 2} (1 + \epsilon_ {3}) (1 + \epsilon_ {4}) \\ 0 & a _ {2 2} (1 + \epsilon_ {5}) \end{array} \right]
$$

and

$$
\hat {B} = \left[ \begin{array}{c c} b _ {1 1} (1 + \epsilon_ {1}) & b _ {1 2} (1 + \epsilon_ {2}) (1 + \epsilon_ {4}) \\ 0 & b _ {2 2} \end{array} \right],
$$

then it is easily verified that $\mathsf { f l } ( A B ) = \hat { A } \hat { B }$ . Moreover,

$$
\hat {A} = A + E, \quad | E | \leq 2 \mathbf {u} | A | + O (\mathbf {u} ^ {2}),
$$

$$
\hat {B} = B + F, \quad | F | \leq 2 \mathbf {u} | B | + O (\mathbf {u} ^ {2}).
$$

which shows that the computed product is the exact product of slightly perturbed A and B.

# 2.7.10 Error in Strassen Multiplication

In §1.3.11 we outlined a recursive matrix multiplication procedure due to Strassen. It is instructive to compare the effect of roundoff in this method with the effect of roundoff in any of the conventional matrix multiplication methods of §1.1.

It can be shown that the Strassen approach (Algorithm 1.3.1) produces a $\hat { C } =$ fl(AB) that satisfies an inequality of the form (2.7.20). This is perfectly satisfactory in many applications. However, the $\hat { C }$ that Strassen’s method produces does not always satisfy an inequality of the form (2.7.19). To see this, suppose that

$$
A = B = \left[ \begin{array}{c c}. 9 9 & . 0 0 1 0 \\ . 0 0 1 0 & . 9 9 \end{array} \right]
$$

and that we execute Algorithm 1.3.1 using 2-digit floating point arithmetic. Among other things, the following quantities are computed:

$$
\hat {P} _ {3} = \mathrm{fl} (. 9 9 (. 0 0 1 -. 9 9)) = -. 9 8,
$$

$$
\hat {P} _ {5} = \mathsf {f l} ((. 9 9 +. 0 0 1). 9 9) = . 9 8,
$$

$$
\hat {c} _ {1 2} = \mathsf {f l} (\hat {P} _ {3} + \hat {P} _ {5}) = 0. 0.
$$

In exact arithmetic $c _ { 1 2 } = 2 ( . 0 0 1 ) ( . 9 9 ) = . 0 0 1 9 8$ and thus Algorithm 1.3.1 produces a $\hat { c } _ { 1 2 }$ with no correct significant digits. The Strassen approach gets into trouble in this example because small off-diagonal entries are combined with large diagonal entries. Note that in conventional matrix multiplication the sums $b _ { 1 2 } + b _ { 2 2 }$ and $a _ { 1 1 } + a _ { 1 2 }$ do not arise. For that reason, the contribution of the small off-diagonal elements is not lost in this example. Indeed, for the above A and B a conventional matrix multiplication gives $\hat { c } _ { 1 2 } = . 0 0 2 0$ .

Failure to produce a componentwise accurate $\hat { C }$ can be a serious shortcoming in some applications. For example, in Markov processes the $a _ { i j } , b _ { i j }$ , and $c _ { i j }$ are transition probabilities and are therefore nonnegative. It may be critical to compute $c _ { i j }$ accurately if it reflects a particularly important probability in the modeled phenomenon. Note that if $A \geq 0$ and $B \geq 0$ , then conventional matrix multiplication produces a product $\hat { C }$ that has small componentwise relative error:

$$
| \hat {C} - C | \leq n \mathbf {u} | A | | B | + O \left(\mathbf {u} ^ {2}\right) = n \mathbf {u} | C | + O \left(\mathbf {u} ^ {2}\right).
$$

This follows from (2.7.19). Because we cannot say the same for the Strassen approach, we conclude that Algorithm 1.3.1 is not attractive for certain nonnegative matrix multiplication problems if relatively accurate $\hat { c } _ { i j }$ are required.

Extrapolating from this discussion we reach two fairly obvious but important conclusions:

• Different methods for computing the same quantity can produce substantially different results.   
• Whether or not an algorithm produces satisfactory results depends upon the type of problem solved and the goals of the user.

These observations are clarified in subsequent chapters and are intimately related to the concepts of algorithm stability and problem condition. See §3.4.10.

# 2.7.11 Analysis of an Ideal Equation Solver

A nice way to conclude this chapter and to anticipate the next is to analyze the quality of a “make-believe” $A x = b$ solution process in which all floating point operations are performed exactly except the storage of the matrix A and the right-hand-side b. It follows that the computed solution ˆx satisfies

$$
(A + E) \hat {x} = (b + e), \quad \| E \| _ {\infty} \leq \mathbf {u} \| A \| _ {\infty}, \quad \| e \| _ {\infty} \leq \mathbf {u} \| b \| _ {\infty}. \tag {2.7.21}
$$

where

$$
\mathsf {f l} (b) = b + e, \quad \mathsf {f l} (A) = A + E.
$$

If u $\kappa _ { \infty } ( A ) \leq { \frac { 1 } { 2 } }$ (say), then by Theorem 2.6.2 it can be shown that

$$
\frac {\| x - \hat {x} \| _ {\infty}}{\| x \| _ {\infty}} \leq 4 \mathbf {u} \kappa_ {\infty} (A). \tag {2.7.22}
$$

The bounds (2.7.21) and (2.7.22) are “best possible” norm bounds. No general ∞- norm error analysis of a linear equation solver that requires the storage of A and b can render sharper bounds. As a consequence, we cannot justifiably criticize an algorithm for returning an inaccurate ˆx if A is ill-conditioned relative to the unit roundoff, e.g., u $\kappa _ { \infty } ( A ) \approx 1$ . On the other hand, we have every “right” to pursue the development of a linear equation solver that renders the exact solution to a nearby problem in the style of (2.7.21).

# Problems

P2.7.1 Show that if (2.7.8) is applied with $y = x$ , then $\mathsf { f l } ( x ^ { T } x ) = x ^ { T } x ( 1 + \alpha )$ where $| \alpha | \le n \mathbf { u } + O ( \mathbf { u } ^ { 2 } )$

P2.7.2 Prove (2.7.4) assuming that $\operatorname { f l } ( x )$ is the nearest floating point number to $x \in \mathbb { R }$

P2.7.3 Show that if $E \in \mathbb { R } ^ { m \times n }$ with $m \geq n$ , then $\| \mathbf { \epsilon } | \| _ { 2 } \leq { \sqrt { n } } \| E \| _ { 2 }$ . This result is useful when deriving norm bounds from absolute value bounds.

P2.7.4 Assume the existence of a square root function satisfying ${ \mathfrak { f l } } ( { \sqrt { x } } ) = { \sqrt { x } } ( 1 + \epsilon )$ with $| \epsilon | \le \mathbf { u }$ Give an algorithm for computing $\parallel x \parallel _ { 2 }$ and bound the rounding errors.

P2.7.5 Suppose A and B are n-by-n upper triangular floating point matrices. If $\hat { C } = \mathfrak { f l } ( A B )$ is computed using one of the conventional §1.1 algorithms, does it follow that $\hat { C } = \hat { A } \hat { B }$ where Aˆ and Bˆ are close to A and B?

P2.7.6 Suppose A and B are n-by-n floating point matrices and that $\parallel | A ^ { - 1 } | | A | \parallel _ { \infty } = \tau .$ . Show that if $\hat { C } = \mathsf { f l } ( A B )$ is obtained using any of the §1.1 algorithms, then there exists a Bˆ so that $\hat { C } = A \hat { B }$ and $\| \hat { B } - B \| _ { \infty } \leq n \mathbf { u } \tau \| B \| _ { \infty } + O ( \mathbf { u } ^ { 2 } )$ .

P2.7.7 Prove (2.7.19).

P2.7.8 For the IEEE double format, what is the largest power of 10 that can be represented exactly? What is the largest integer that can be represented exactly?

P2.7.9 For $k = 1 { : } 6 2$ , what is the largest power of 10 that can be stored exactly if k bits are are allocated for the mantissa and $6 3 - k$ are allocated for the exponent?

P2.7.10 Consider the quadratic equation

$$
q (\lambda) = \det \left(\left[ \begin{array}{c c} w - \lambda & x \\ x & z - \lambda \end{array} \right]\right).
$$

This quadratic has two real roots $r _ { 1 }$ and $r _ { 2 }$ . Assume that $| r _ { 1 } - z | \leq | r _ { 2 } - z |$ . Give an algorithm that computes $r _ { 1 }$ to full machine precision.

# Notes and References for §2.7

For an excellent, comprehensive treatment of IEEE arithmetic and its implications, see:

M.L. Overton (2001). Numerical Computing with IEEE Arithmetic, SIAM Publications, Philadelphia, PA.

The following basic references are notable for the floating point insights that they offer: Wilkinson (AEP), Stewart (IMC), Higham (ASNA), and Demmel (ANLA). For high-level perspectives we recommend:

J.H. Wilkinson (1963). Rounding Errors in Algebraic Processes, Prentice-Hall, Englewood Cliffs, NJ.

G.E. Forsythe (1970). “Pitfalls in Computation or Why a Math Book is Not Enough,” Amer. Math. Monthly 77, 931–956.

J.H. Wilkinson (1971). “Modern Error Analysis,” SIAM Review 13, 548–68.

U.W. Kulisch and W.L. Miranker (1986). “The Arithmetic of the Digital Computer,” SIAM Review 28, 1–40.

F. Chaitin-Chatelin and V. Frayse´e (1996). Lectures on Finite Precision Computations, SIAM Publications, Philadelphia, PA.

The design of production software for matrix computations requires a detailed understanding of finite precision arithmetic, see:

J.W. Demmel (1984). “Underflow and the Reliability of Numerical Software,” SIAM J. Sci. Stat. Comput. 5, 887–919.

W.J. Cody (1988). “ALGORITHM 665 MACHAR: A Subroutine to Dynamically Determine Machine Parameters,” ACM Trans. Math. Softw. 14, 303–311.

D. Goldberg (1991). “What Every Computer Scientist Should Know About Floating Point Arithmetic,” ACM Surveys 23, 5–48.

Other developments in error analysis involve interval analysis, the building of statistical models of roundoff error, and the automating of the analysis itself:

J. Larson and A. Sameh (1978). “Efficient Calculation of the Effects of Roundoff Errors,” ACM Trans. Math. Softw. 4, 228–36.

W. Miller and D. Spooner (1978). “Software for Roundoff Analysis, II,” ACM Trans. Math. Softw. 4, 369–90.

R.E. Moore (1979). Methods and Applications of Interval Analysis, SIAM Publications, Philadelphia, PA.

J.M. Yohe (1979). “Software for Interval Arithmetic: A Reasonable Portable Package,” ACM Trans. Math. Softw. 5, 50–63.

The accuracy of floating point summation is detailed in:

S.M. Rump, T. Ogita, and S. Oishi (2008). “Accurate Floating-Point Summation Part I: Faithful Rounding,” SIAM J. Sci. Comput. 31, 189–224.

S.M. Rump, T. Ogita, and S. Oishi (2008). “Accurate Floating-Point Summation Part II: Sign, K-fold Faithful and Rounding to Nearest,” SIAM J. Sci. Comput. 31, 1269–1302.   
For an analysis of the Strassen algorithm and other “fast” linear algebra procedures, see:   
R.P. Brent (1970). “Error Analysis of Algorithms for Matrix Multiplication and Triangular Decomposition Using Winograd’s Identity,” Numer. Math. 16, 145–156.   
W. Miller (1975). “Computational Complexity and Numerical Stability,” SIAM J. Comput. 4, 97–107.   
N.J. Higham (1992). “Stability of a Method for Multiplying Complex Matrices with Three Real Matrix Multiplications,” SIAM J. Matrix Anal. Applic. 13, 681–687.   
J.W. Demmel and N.J. Higham (1992). “Stability of Block Algorithms with Fast Level-3 BLAS,” ACM Trans. Math. Softw. 18, 274–291.   
B. Dumitrescu (1998). “Improving and Estimating the Accuracy of Strassen’s Algorithm,” Numer. Math. 79, 485–499.   
The issue of extended precision has received considerable attention. For example, a superaccurate dot product results if the summation can be accumulated in a register that is “twice as wide” as the floating representation of vector components. The overhead may be tolerable in a given algorithm if extended precision is needed in only a few critical steps. For insights into this topic, see:   
R.P. Brent (1978). “A Fortran Multiple Precision Arithmetic Package,” ACM Trans. Math. Softw. 4, 57–70.   
R.P. Brent (1978). “Algorithm 524 MP, a Fortran Multiple Precision Arithmetic Package,” ACM Trans. Math. Softw. 4, 71–81.   
D.H. Bailey (1993). “Algorithm 719: Multiprecision Translation and Execution of FORTRAN Programs,” ACM Trans. Math. Softw. 19, 288–319.   
X.S. Li, J.W. Demmel, D.H. Bailey, G. Henry, Y. Hida, J. Iskandar, W. Kahan, S.Y. Kang, A. Kapur, M.C. Martin, B.J. Thompson, T. Tung, and D.J. Yoo (2002). “Design, Implementation and Testing of Extended and Mixed Precision BLAS,” ACM Trans. Math. Softw. 28, 152–205.   
J.W. Demmel and Y. Hida (2004). “Accurate and Efficient Floating Point Summation,” SIAM J. Sci. Comput. 25, 1214–1248.   
M. Baboulin, A. Buttari, J. Dongarra, J. Kurzak, J. Langou, J. Langou, P. Luszczek, and S. Tomov (2009). “Accelerating Scientific Computations with Mixed Precision Algorithms,” Comput. Phys. Commun. 180, 2526–2533.
