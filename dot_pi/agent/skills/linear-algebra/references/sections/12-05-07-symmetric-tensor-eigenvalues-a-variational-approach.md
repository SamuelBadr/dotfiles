# 12.5.7 Symmetric Tensor Eigenvalues: A Variational Approach

If $C \in \mathbb { R } ^ { N \times N }$ is symmetric, then its eigenvalues are the stationary values of

$$
\phi_ {C} (x) = \frac {x ^ {T} C x}{x ^ {T} x} = \frac {\sum_ {i _ {1} = 1} ^ {N} \sum_ {i _ {2} = 1} ^ {N} C (i _ {1} , i _ {2}) x (i _ {1}) x (i _ {2})}{x ^ {T} x} \tag {12.5.23}
$$

and the corresponding stationary vectors are eigenvectors. This follows by setting the gradient of $\phi _ { C }$ to zero.

If we are to generalize this notion to tensors, then we need to define what we mean by a symmetric tensor. An order-d tensor $\mathcal { C } \in \mathbb { R } ^ { N \times \cdots \times N }$ is symmetric if for any permutation p of 1:d we have

$$
\mathcal {C} (\mathbf {i}) = \mathcal {C} (\mathbf {i} (\mathbf {p})), \quad \mathbf {1} \leq \mathbf {i} \leq_ {N}.
$$

For the case $d = 3$ this means $c _ { i j k } \ = \ c _ { i k j } \ = \ c _ { j i k } \ = \ c _ { j k i } \ = \ c _ { k i j } \ = \ c _ { k j i }$ for all $i , j ,$ , and k that satisfy $1 \leq i \leq N , 1 \leq j \leq N$ , and $1 \leq k \leq N$ .

It is easy to generalize (12.5.23) to the case of symmetric tensors. If $\mathcal { C } \in \mathbb { R } ^ { N \times N \times N }$ is symmetric and $\boldsymbol { x } \in \mathbb { R } ^ { N }$ then we define φ by

$$
\phi_ {\mathcal {C}} (x) = \frac {\sum_ {\mathbf {i} = \mathbf {1}} ^ {\mathbf {N}} \mathcal {C} (\mathbf {i}) \cdot x (i _ {1}) x (i _ {2}) x (i _ {3})}{\| x \| _ {2} ^ {3}} = \frac {x ^ {T} \mathcal {C} _ {(1)} (x \otimes x)}{\| x \| _ {2} ^ {3}}. \tag {12.5.24}
$$

Note that if $\mathcal { C }$ is a symmetric tensor, then all its modal unfoldings are the same. The equation $\nabla \phi _ { C } ( x ) = 0$ with $\parallel x \parallel _ { 2 } = 1$ has the form

$$
\nabla \phi_ {\mathcal {C}} (x) = \mathcal {C} _ {(1)} (x \otimes x) - \phi_ {\mathcal {C}} (x) \cdot x = 0.
$$

If this holds then we refer to $\phi _ { \mathcal { C } } ( x )$ as an eigenvalue of the tensor ${ \mathcal { C } } ,$ a concept introduced by Lim (2005) and Li (2005). An interesting framework for solving this nonlinear equation has been proposed by Kolda and Mayo (2012). It involves repetition of the operation sequence

$$
\tilde {x} = \mathcal {C} _ {(1)} (x \otimes x) + \alpha x, \quad \lambda = \| \tilde {x} \| _ {2}, \quad x = \tilde {x} / \lambda
$$

where the shift parameter α is determined to ensure convexity and eventual convergence of the iteration. For further discussion of the symmetric tensor eigenvalue problem and various power iterations that can be used to solve it, see Zhang and Golub (2001) and Kofidis and Regalia (2002).
