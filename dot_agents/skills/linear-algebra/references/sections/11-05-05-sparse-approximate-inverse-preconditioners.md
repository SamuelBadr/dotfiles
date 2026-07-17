# 11.5.5 Sparse Approximate Inverse Preconditioners

Instead of determining M so $\| A - M \| _ { F }$ is small, we can address Criterion 1 above by choosing $M ^ { - 1 }$ so that $\left. \mathbf { \nabla } A M ^ { - 1 } - I \right. _ { F }$ is small. This is the idea behind sparse approximate inverse preconditioners. To be precise about the nature of the approximation, we define the $\mathsf { s p } ( \cdot )$ operator. For any $T \in \mathbb { R } ^ { n \times n }$ define $\mathsf { s p } ( T ) \in \mathbb { R } ^ { n \times n }$ by

$$
[ \mathfrak {s p} (T) ] _ {i j} = \left\{ \begin{array}{l l} 1 & \text {if} t _ {i j} \neq 0 \\ 0 & \text {otherwise} \end{array} \right..
$$

Suppose $Z \in \mathbb { R } ^ { n \times n }$ is a given n-by-n matrix of zeros and ones with a manageable sparsity pattern and that we solve the constrained least squares problem

$$
\min _ {\mathfrak {s p} (T) = Z} \left\| A T - I \right\| _ {F}.
$$

The constraint says that T is to have the same zero-nonzero structure as $Z .$ . Thus, the preconditioner M is specified through its inverse: $M ^ { - 1 } = T$ . A fringe benefit of this type of preconditioner design is that the $M z = r$ system is solved via matrix-vector multiplication: $z = T r$ . This is what makes this preconditioning approach attractive from the parallel computing point of view. Moreover, the actual columns of $T$ can be computed in parallel because they are independent of each other.

It is important to appreciate that $T ( : , k )$ is a constrained minimizer of $\parallel A \tau - e _ { k } \parallel _ { 2 }$ . Let cols be the subvector of 1:n that identifies the nonzero components of $T ( : , k )$ . (These indices are determined by $Z ( : , k ) . )$ Let rows be a subset of 1:n that identifies the nonzero rows in $A ( : , c o l s )$ . If τ solves the (generally very small) LS problem

$$
\min \left\| A (r o w s, c o l s) \tau - e _ {k} (r o w s) \right\| _ {2}
$$

then $T ( : , k )$ is zero except $T ( r o w s , k ) = \tau$ . We mention that the sparsity pattern Z can be determined dynamically. For example, after completing the above column-k calculation, it is possible to expand col cheaply to include more nonzeros in $T ( : , k )$ . See Grote and Huckle (1997). Updating QR factorizations is part of their method.
