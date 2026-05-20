# 11.1.9 Sparse LU

The first step in a pivoted LU procedure applied to $A \in \mathbb { R } ^ { n \times n }$ computes the factorization

$$
P A Q ^ {T} = \left[ \begin{array}{c c} \alpha & w ^ {T} \\ v & B \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ v / \alpha & I _ {n - 1} \end{array} \right] \left[ \begin{array}{c c} \alpha & w ^ {T} \\ 0 & A ^ {(1)} \end{array} \right] \tag {11.1.10}
$$

where P and $Q$ are permutation matrices and

$$
A ^ {(1)} = B - \frac {1}{\alpha} v w ^ {T}. \tag {11.1.11}
$$

In §3.4 we discussed various choices for P and Q. Stability was the primary issue and everything revolved around making the pivot element α sufficiently large. If A is sparse, then in addition to stability we have to be concerned about the sparsity of A(1). Balancing the tension between stability and sparsity defines the Sparse LU challenge:

# The Sparse LU Challenge

Given a matrix $A \in \mathbb { R } ^ { n \times n }$ , efficiently determine permutations $p$ and $q$ of 1:n so that if $P = I _ { n } ( : , p )$ and $Q = I _ { n } ( : , q )$ , then the factorization $A ( p , q ) = P A Q ^ { T } = L U$ is reasonably stable and the triangular factors $L$ and U are close to being optimally sparse.

To meet the challenge we must interpolate between a pair of extreme strategies:

• Maximize stability by choosing P and $Q$ so that $\left| \alpha \right| = \operatorname* { m a x } \left| a _ { i j } \right|$ .   
• Maximize sparsity by choosing P and Q so that $\mathsf { n n z } \big ( A ^ { ( 1 ) } \big )$ is minimized.

Markowitz pivoting provides a framework for doing this. Given a threshold parameter $\tau$ that satisfies $0 \leq \tau \leq 1$ , choose $P$ and $Q$ in each step of the form (11.1.10) so that $\mathsf { n n z } \big ( A ^ { ( 1 ) } \big )$ is minimized subject to the constraint that $| \alpha | \geq \tau | v _ { i } |$ for $i = 1 { : } n - 1$ . Small values of $\tau$ jeopardize stability but create more opportunities to control fill-in. A typical compromise value is $\tau = 1 / 1 0$ .

Sometimes there is an advantage to choosing the pivot from the diagonal, i.e., setting $P = Q$ . This is the case when the matrix A is structurally symmetric. A matrix A is structurally symmetric if $a _ { i j }$ and $a _ { j i }$ are either both zero or both nonzero. Symmetric matrices whose rows and/or columns are scaled have this property. It is easy so show from (11.1.10) and (11.1.11) that if A is structurally symmetric and $P = Q$ , then $A ^ { ( 1 ) }$ is structurally symmetric. The Markowitz strategy can be generalized to express a preference for diagonal pivoting if it is “safe”. If a diagonal element is sufficiently large compared to other entries in its column, then $P$ is chosen so that $( P A P ^ { T } ) _ { 1 1 }$ is that element and structural symmetry is preserved. Otherwise, a sufficiently large off-diagonal element is brought to the (1,1) position using a $P A Q ^ { T }$ update.

# Problems

P11.1.1 Give an algorithm that solves an upper triangular system $T x = b$ given that T is stored in the compressed-column format.   
P11.1.2 If both indexing and flops are taken into consideration, is the sparse outer-product update (11.1.2) an $O ( \mathsf { n n z } ( u ) \cdot \mathsf { n n z } ( v ) )$ computation?   
P11.1.3 For example (11.1.5), what is the resulting profile if $S _ { 0 } = \{ 9 \} \ ?$ What if $S _ { 0 } = \{ 4 \} \{$   
P11.1.4 Prove that the Cuthill-McKee ordering permutes A into a block tridiagonal form where the kth diagonal block is r-by-r where r is the cardinality of $S _ { k - 1 }$ .   
P11.1.5 (a) What is the resulting profile if the reverse Cuthill-McKee ordering is applied to the example in §11.1.5? (b) What is the elimination tree for the matrix in (11.1.5)?   
P11.1.6 Show that if G is the Cholesky factor of A and an element $g _ { i j } \neq 0$ , then $j \geq f _ { i }$ where $f _ { i }$ is defined by (11.1.6). Conclude that $\mathsf { n n z } ( G ) \leq \mathsf { p r o f i l e } ( A )$ .   
P11.1.7 Show how the method of seminormal equations can be used efficiently to minimize $\parallel M x - d \parallel _ { 2 }$ where

$$
M   =   \left[ \begin{array}{c c c c} A _ {1} & 0 & 0 & C _ {1} \\ 0 & A _ {2} & 0 & C _ {2} \\ 0 & 0 & A _ {3} & C _ {3} \end{array} \right], \qquad d   =   \left[ \begin{array}{c} b _ {1} \\ b _ {2} \\ b _ {3} \end{array} \right],
$$

and $A _ { i } \in \mathbb { R } ^ { m \times n } , C _ { i } \in \mathbb { R } ^ { m \times p } .$ , and $b _ { i } \in \mathbb { R } ^ { m }$ for $i = 1 { : } 3$ . Assume that M has full column rank and that $m > n + p$ . Hint: Compute the Q-less QR factorizations of $[ A _ { i } \thinspace C _ { i } ]$ for $i = 1 { : } 3$ .

# Notes and References for 11.1

Early references for direct sparse matrix computations include the following textbooks:

A. George and J.W.-H. Liu (1981). Computer Solution of Large Sparse Positive Definite Systems, Prentice-Hall, Englewood Cliffs, NJ.

O. Osterby and Z. Zlatev (1983). Direct Methods for Sparse Matrices, Springer-Verlag, New York.   
S. Pissanetzky (1984). Sparse Matrix Technology, Academic Press, New York.   
I.S. Duff, A.M. Erisman, and J.K. Reid (1986). Direct Methods for Sparse Matrices, Oxford University Press, London.   
A more recent treatment that targets practitioners, provides insight into a range of implementation issues, and has an excellent annotated bibliography is the following:   
T.A. Davis (2006). Direct Methods for Sparse Linear Systems, SIAM Publications, Philadelphia, PA.   
The interplay between graph theory and sparse matrix computations with emphasis on symbolic factorizations that predict fill is nicely set forth in:   
J.W.H. Liu (1990). “The Role of Elimination Trees in Sparse Factorizations,” SIAM J. Matrix Anal. Applic. 11, 134–172.   
J.R. Gilbert (1994). “Predicting Structure in Sparse Matrix Computations,” SIAM J. Matrix Anal. Applic. 15, 62–79.   
S.C. Eisenstat and J.W.H. Liu (2008). “Algorithmic Aspects of Elimination Trees for Sparse Unsymmetric Matrices,” SIAM J. Matrix Anal. Applic. 29, 1363–1381.   
Relatively recent papers on profile reduction include:   
W.W. Hager (2002). “Minimizing the Profile of a Symmetric Matrix,” SIAM J. Sci. Comput. 23, 1799–1816.   
J.K. Reid and J.A. Scott (2006). “Reducing the Total Bandwidth of a Sparse Unsymmetric Matrix,” SIAM J. Matrix Anal. Applic. 28, 805–821.   
Efficient implementations of the minimum degree idea are discussed in:   
P.R. Amestoy, T.A. Davis, and I.S. Duff (1996). “An Approximate Minimum Degree Ordering Algorithm,” SIAM J. Matrix Anal. Applic. 17, 886–905.   
T.A. Davis, J.R. Gilbert, S.I. Larimore, and E.G. Ng (2004). “A Column Approximate Minimum Degree Ordering Algorithm,” ACM Trans. Math. Softw. 30, 353–376.   
For an overview of sparse least squares, see Bj¨orck (NMLS, Chap. 6)) and also:   
J.A. George and M.T. Heath (1980). “Solution of Sparse Linear Least Squares Problems Using Givens Rotations,” Lin. Alg. Applic. 34, 69–83.   
˚A. Bj¨orck and I.S. Duff (1980). “A Direct Method for the Solution of Sparse Linear Least Squares Problems,” Lin. Alg. Applic. 34, 43–67.   
A. George and E. Ng (1983). “On Row and Column Orderings for Sparse Least Squares Problems,” SIAM J. Numer. Anal. 20, 326–344.   
M.T. Heath (1984). “Numerical Methods for Large Sparse Linear Least Squares Problems,” SIAM J. Sci. Stat. Comput. 5, 497–513.   
˚A. Bj¨orck (1987). “Stability Analysis of the Method of Seminormal Equations for Least Squares Problems,” Lin. Alg. Applic. 88/89, 31–48.   
The design of a sparse LU procedure that is also stable is discussed in:   
J.W. Demmel, S.C. Eisenstat, J.R. Gilbert, X.S. Li, and J.W.H. Liu (1999). “A Supernodal Approach to Sparse Partial Pivoting,” SIAM J. Matrix Anal. Applic. 20, 720–755.   
L. Grigori, J.W. Demmel, and X.S. Li (2007). “Parallel Symbolic Factorization for Sparse LU with Static Pivoting,” SIAM J. Sci. Comput. 3, 1289–1314.   
L. Grigori, J.R. Gilbert, and M. Cosnard (2008). “Symbolic and Exact Structure Prediction for Sparse Gaussian Elimination with Partial Pivoting,” SIAM J. Matrix Anal. Applic. 30, 1520–1545.   
Frontal methods are a way of organizing outer-product updates so that the resulting implementation is rich in dense matrix operations, a maneuver that is critical from the standpoint of performance, see:   
J.W.H. Liu (1992). “The Multifrontal Method for Sparse Matrix Solution: Theory and Practice,” SIAM Review 34, 82–109.   
D.J. Pierce and J.G. Lewis (1997). “Sparse Multifrontal Rank Revealing QR Factorization,” SIAM J. Matrix Anal. Applic. 18, 159–180.   
T.A. Davis and I.S. Duff (1999). “A Combined Unifrontal/Multifrontal Method for Unsymmetric Sparse Matrices,” ACM Trans. Math. Softw. 25, 1–20.

Another important reordering challenge involves permuting to block triangular form, see:

A. Pothen and C.-J. Fan (1990). “Computing the Block Triangular Form of a Sparse Matrix,” ACM Trans. Math. Softw. 16, 303–324.   
I.S. Duff and B. U¸car (2010). “On the Block Triangular Form of Symmetric Matrices,” SIAM Review 52, 455–470.   
Early papers on parallel sparse matrix computations that are filled with interesting ideas include:   
M.T. Heath, E. Ng, and B.W. Peyton (1991). “Parallel Algorithms for Sparse Linear Systems,” SIAM Review 33, 420–460.   
J.R. Gilbert and R. Schreiber (1992). “Highly Parallel Sparse Cholesky Factorization,” SIAM J. Sci. Stat. Comput. 13, 1151–1172.   
For a sparse-matrix discussion of condition estimation, error analysis, and related problems, see:   
R.G. Grimes and J.G. Lewis (1981). “Condition Number Estimation for Sparse Matrices,” SIAM J. Sci. Stat. Comput. 2, 384–388.   
M. Arioli, J.W. Demmel, and I.S. Duff (1989). “Solving Sparse Linear Systems with Sparse Backward error,” SIAM J. Matrix Anal. Applic. 10, 165–190.   
C.H. Bischof (1990). “Incremental Condition Estimation for Sparse Matrices,” SIAM J. Matrix Anal. Applic. 11, 312–322.   
M.W. Berry, S.A. Pulatova, and G.W. Stewart (2005). “Algorithm 844: Computing Sparse Reduced-Rank Approximations to Sparse Matrices,” ACM Trans. Math. Softw. 31, 252–269.
