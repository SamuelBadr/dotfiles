# 11.5.4 Normwise-Near Preconditioners

Sometimes A is near a data-sparse matrix for which there is a fast solution procedure. Circulant preconditioners for symmetric Toeplitz systems are a nice example. For $a \in \mathbb { R } ^ { n }$ define the Toeplitz matrix $\ b { T } ( a ) \in \mathbb { R } ^ { n \times n }$ and the circulant matrix $C ( a ) \in \mathbb { R } ^ { n \times n }$ by

$$
T (a) = \left[ \begin{array}{c c c c} a _ {0} & a _ {1} & a _ {2} & a _ {3} \\ a _ {1} & a _ {0} & a _ {1} & a _ {2} \\ a _ {2} & a _ {1} & a _ {0} & a _ {1} \\ a _ {3} & a _ {2} & a _ {1} & a _ {0} \end{array} \right], \qquad C (a) = \left[ \begin{array}{c c c c} a _ {0} & a _ {1} & a _ {2} & a _ {3} \\ a _ {3} & a _ {0} & a _ {1} & a _ {2} \\ a _ {2} & a _ {3} & a _ {0} & a _ {1} \\ a _ {1} & a _ {2} & a _ {3} & a _ {0} \end{array} \right], \qquad (n = 4).
$$

Suppose we determine ˜a so that $\parallel T ( a ) - C ( \tilde { a } ) \parallel _ { F }$ is minimized. A case can be made that $M = C ( \tilde { \boldsymbol { a } } )$ captures the essence of $T ( a )$ and thus has potential as a preconditioner for the Toeplitz system $T ( a ) x = b$ . Recall from §4.8.2 that circulant linear systems can be solved in n log n time using the fast Fourier transform. This style of Toeplitz system preconditioning was proposed by Chan (1988).

Because of their importance, there is a large body of work concerned with preconditioners for Toeplitz systems. An idea due to Chan and Strang (1989) is to set $M = C ( \tilde { \boldsymbol { a } } )$ where

$$
\tilde {a} = \left[ \begin{array}{c} a (0: m) \\ a (m - 1: - 1: 0) \end{array} \right]
$$

assuming that $n = 2 m$ and $A = T ( a )$ is positive definite. Intuition tells us that A’s central diagonals carry most of the information and so it makes sense that they define the preconditioner $C ( \tilde { a } )$ .
