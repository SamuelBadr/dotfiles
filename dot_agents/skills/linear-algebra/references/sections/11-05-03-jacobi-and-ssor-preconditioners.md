# 11.5.3 Jacobi and SSOR Preconditioners

We now begin a tour of the major preconditioning strategies. Since some strategies help motivate others, the order of presentation is pedagogical. It does not indicate relative importance, nor does it reflect how the research on preconditioning evolved.

Suppose $A \in \mathbb { R } ^ { n \times n }$ is diagonally dominant or positive definite. For such a matrix, the diagonal tells much of the story and so it makes a certain amount of sense to consider perhaps the simplest preconditioner of all:

$$
M = \mathrm{diag} (a _ {1 1}, \ldots , a _ {n n}).
$$

Diagonal preconditioners are called Jacobi preconditioners. Recall from §11.2.2 that Jacobi’s method is based on the splitting $A = M - N$ where M is the diagonal of A. Indeed, for any iteration of the form $M x _ { + } = N x _ { c } + b $ , we can regard M as a preconditioner. The requirement that

$$
\rho (M ^ {- 1} N) = \rho (M ^ {- 1} (M - A)) = \rho (I - M ^ {- 1} A) <   1
$$

is a way of saying that $M ^ { - 1 }$ must “look like” $A ^ { - 1 }$ . In this context, the SSOR preconditioner

$$
M = (D - \omega L) D ^ {- 1} (D - \omega L) ^ {T}
$$

is attractive for certain symmetric positive definite systems. Note that in this case M is also symmetric positive definite and so it can be used with PCG.

If $A = ( A _ { i j } )$ is a p-by-p block matrix that is (block) diagonally dominant or positive definite, then the block Jacobi preconditioner $M = \mathrm { d i a g } ( A _ { 1 1 } , \ldots , A _ { p p } )$ is sometimes a natural choice.
