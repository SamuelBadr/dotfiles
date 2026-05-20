# Chapter 1

# Matrix Multiplication

1.1 Basic Algorithms and Notation   
1.2 Structure and Efficiency   
1.3 Block Matrices and Algorithms   
1.4 Fast Matrix-Vector Products   
1.5 Vectorization and Locality   
1.6 Parallel Matrix Multiplication

The study of matrix computations properly begins with the study of various matrix multiplication problems. Although simple mathematically, these calculations are sufficiently rich to develop a wide range of essential algorithmic skills.

In §1.1 we examine several formulations of the matrix multiplication update problem $C = C + A B$ . Partitioned matrices are introduced and used to identify linear algebraic “levels” of computation.

If a matrix has special properties, then various economies are generally possible. For example, a symmetric matrix can be stored in half the space of a general matrix. A matrix-vector product may require much less time to execute if the matrix has many zero entries. These matters are considered in §1.2.

A block matrix is a matrix whose entries are themselves matrices. The “language” of block matrices is developed in §1.3. It supports the easy derivation of matrix factorizations by enabling us to spot patterns in a computation that are obscured at the scalar level. Algorithms phrased at the block level are typically rich in matrixmatrix multiplication, the operation of choice in many high-performance computing environments. Sometimes the block structure of a matrix is recursive, meaning that the block entries have an exploitable resemblance to the overall matrix. This type of connection is the foundation for “fast” matrix-vector product algorithms such as various fast Fourier transforms, trigonometric transforms, and wavelet transforms. These calculations are among the most important in all of scientific computing and are discussed in §1.4. They provide an excellent opportunity to develop a facility with block matrices and recursion.

The last two sections set the stage for effective, “large-n” matrix computations. In this context, data locality affects efficiency more than the volume of actual arithmetic. Having an ability to reason about memory hierarchies and multiprocessor computation is essential. Our goal in §1.5 and §1.6 is to build an appreciation for the attendant issues without getting into system-dependent details.

# Reading Notes

The sections within this chapter depend upon each other as follows:

$$
\begin{array}{c c c c c c c c} \S 1. 1 & \to & \S 1. 2 & \to & \S 1. 3 & \to & \S 1. 4 \\ & & & & \downarrow & & \\ & & & & \S 1. 5 & \to & \S 1. 6 \end{array}
$$

Before proceeding to later chapters, §1.1, §1.2, and §1.3 are essential. The fast transform ideas in §1.4 are utilized in §4.8 and parts of Chapters 11 and 12. The reading of §1.5 and §1.6 can be deferred until high-performance linear equation solving or eigenvalue computation becomes a topic of concern.
