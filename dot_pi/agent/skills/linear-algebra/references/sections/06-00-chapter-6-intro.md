# Chapter 6

# Modified Least Squares Problems and Methods

6.1 Weighting and Regularization   
6.2 Constrained Least Squares   
6.3 Total Least Squares   
6.4 Subspace Computations with the SVD   
6.5 Updating Matrix Factorizations

In this chapter we discuss an assortment of least square problems that can be solved using QR and SVD. We also introduce a generalization of the SVD that can be used to simultaneously diagonalize a pair of matrices, a maneuver that is useful in certain applications.

The first three sections deal with variations of the ordinary least squares problem that we treated in Chapter 5. The unconstrained minimization of $\parallel A x - b \parallel _ { 2 }$ does not always make a great deal of sense. How do we balance the importance of each equation in $\boldsymbol { A } \boldsymbol { x } = \boldsymbol { b } \boldsymbol { ? }$ How might we control the size of x if A is ill-conditioned? How might we minimize $\parallel A x - b \parallel _ { 2 }$ over a proper subspace of $\mathbb { R } ^ { n } ?$ What if there are errors in the “data matrix” A in addition to the usual errors in the “vector of observations” b?

In §6.4 we consider a number of multidimensional subspace computations including the problem of determining the principal angles between a pair of given subspaces. The SVD plays a prominent role.

The final section is concerned with the updating of matrix factorizations. In many applications, one is confronted with a succession of least squares (or linear equation) problems where the matrix associated with the current step is highly related to the matrix associated with the previous step. This opens the door to updating strategies that can reduce factorization overheads by an order of magnitude.

# Reading Notes

Knowledge of Chapter 5 is assumed. The sections in this chapter are independent of each other except that §6.1 should be read before §6.2. Excellent global references include Bj¨orck (NMLS) and Lawson and Hansen (SLS).
