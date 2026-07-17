# 4.2 Positive Definite Systems

A matrix $A \in \mathbb { R } ^ { n \times n }$ is positive definite if $x ^ { T } A x > 0$ for all nonzero $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , positive semidefinite if $x ^ { T } A x \geq 0$ for all $\boldsymbol { x } \in \mathbb { R } ^ { n }$ , and indefinite if we can find $x , y \in \mathbb { R } ^ { n }$ so $\left( x ^ { T } A x \right) \left( y ^ { T } A y \right) < 0 . \mathrm { ~ S } _  \mathrm  ~ \scriptsize ~  ~$ ymmetric positive definite systems constitute one of the most important classes of special $A x = b$ problems. Consider the 2-by-2 symmetric case. If

$$
A = \left[ \begin{array}{l l} \alpha & \beta \\ \beta & \gamma \end{array} \right]
$$

is positive definite then

$$
x = [ 1, 0 ] ^ {T} \Rightarrow x ^ {T} A x = \alpha > 0,
$$

$$
x = [ 0, 1 ] ^ {T} \quad \Rightarrow x ^ {T} A x = \gamma > 0,
$$

$$
x = [ 1, 1 ] ^ {T} \Rightarrow x ^ {T} A x = \alpha + 2 \beta + \gamma > 0,
$$

$$
x = [ 1, - 1 ] ^ {T} \Rightarrow x ^ {T} A x = \alpha - 2 \beta + \gamma > 0.
$$

The last two equations imply $\left| \beta \right| \le ( \alpha { + } \gamma ) / 2$ . From these results we see that the largest entry in A is on the diagonal and that it is positive. This turns out to be true in general. (See Theorem 4.2.8 below.) A symmetric positive definite matrix has a diagonal that is sufficiently “weighty” to preclude the need for pivoting. A special factorization called the Cholesky factorization is available for such matrices. It exploits both symmetry and definiteness and its implementation is the main focus of this section. However, before those details are pursued we discuss unsymmetric positive definite matrices. This class of matrices is important in its own right and and presents interesting pivot-related issues.
