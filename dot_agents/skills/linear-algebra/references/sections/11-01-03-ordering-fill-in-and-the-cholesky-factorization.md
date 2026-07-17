# 11.1.3 Ordering, Fill-In, and the Cholesky Factorization

The first step in the outer-product Cholesky process involves computation of the factorization

$$
A = \left[ \begin{array}{c c} \alpha & v ^ {T} \\ v & B \end{array} \right] = \left[ \begin{array}{c c} \sqrt {\alpha} & 0 \\ v / \sqrt {\alpha} & I \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ 0 & A ^ {(1)} \end{array} \right] \left[ \begin{array}{c c} \sqrt {\alpha} & v ^ {T} / \sqrt {\alpha} \\ 0 & I \end{array} \right] \tag {11.1.3}
$$

where

$$
A ^ {(1)} = B - \frac {v v ^ {T}}{\alpha}. \tag {11.1.4}
$$

Recall from §4.2 that this reduction is repeated on the matrix $A ^ { ( 1 ) }$ .

Now suppose A is a sparse matrix. From the standpoint of both arithmetic and memory requirements, we have a vested interest in the sparsity of $A ^ { ( 1 ) }$ . Since B is sparse, everything hinges on the sparsity of the vector v. Here are two examples that dramatize what is at stake:

![](images/golub_600_649__94ea9515e69522bb039a994ce491f253f696e8848e2f89f422b350c42d6e82ff.jpg)

In Example 1, the vector v associated with the first step is dense and that results in a full A(1). All sparsity is lost and the remaining steps essentially carry out a dense Cholesky factorization. Example 2 tells a happier story. The first v-vector is sparse and the update matrix $A ^ { ( 1 ) }$ has the same “arrow” structure as A. Note that Example 2 can be obtained from Example 1 by a reordering of the form $P A P ^ { T }$ where $P = I _ { n } ( : , n : - 1 : 1 ) )$ ). This motivates the Sparse Cholesky challenge:

<table><tr><td>The Sparse Cholesky Challenge</td></tr><tr><td>Given a symmetric positive definite matrix  $A \in \mathbb{R}^{n \times n}$ , efficiently determine a permutation  $p$  of 1:n so that if  $P = I_n(:, p)$ , then the Cholesky factor in  $A(p, p) = PAP^T = GG^T$  is close to being optimally sparse.</td></tr></table>

Choosing P to actually minimize nnz(G) is a formidable combinatorial problem and is therefore not a viable option. Fortunately, there are several practical procedures based on heuristics that can be used to determine a good reordering permutation P . These include (1) the Cuthill-McKee ordering, (2) the minimum degree ordering, and (3) the nested dissection ordering. However, before we discuss these strategies, we need to present a few concepts from graph theory.
