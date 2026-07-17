# Chapter 12

# Special Topics

12.1 Linear Systems with Displacement Structure   
12.2 Structured-Rank Problems   
12.3 Kronecker Product Computations   
12.4 Tensor Unfoldings and Contractions   
12.5 Tensor Decompositions and Iterations

Prominent themes in this final chapter include data sparsity, low-rank approximation, exploitation of structure, the importance of representation, and large-scale problems. We revisit (unsymmetric) Toeplitz systems in §12.1 and show how fast stable methods can be developed through a clever data-sparse representation. The ideas extend to other types of structured matrices. Representation is also central to the O(n) methods developed in §12.2 for matrices that have low-rank off-diagonal blocks.

The next three sections form a sequence. The Kronecker product section has general utility, but it is used very heavily in both §12.4 and §12.5 which together provide a brief introduction to the rapidly developing field of tensor computations.

Reading Path

Within this chapter, there are the following dependencies

$$
\begin{array}{c c c c c c c c c} \S 3. 1 \text {-} \S 3. 4, \S 4. 7 & \to & \S 1 2. 1 & & & & \S 5. 1 \text {-} \S 5. 3 \\ \S 3. 1 \text {-} \S 3. 4, \S 5. 1 \text {-} \S 5. 3 & \to & \S 1 2. 2 & & & & \downarrow \\ \S 1. 4 & \to & \S 1 2. 3 & \to & \S 1 2. 4 & \to & \S 1 2. 5 \end{array}
$$

The schematic also hints at the minimum “prerequisites” for each topic.
