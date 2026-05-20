# 11.5.10 Saddle Point System Preconditioners

A nonsingular 2-by-2 block system of the form

$$
K = \left[ \begin{array}{c c} A & B _ {1} ^ {T} \\ B _ {2} & - C \end{array} \right] \left[ \begin{array}{c} x \\ y \end{array} \right] = \left[ \begin{array}{c} f \\ g \end{array} \right],
$$

where $A \in \mathbb { R } ^ { n \times n }$ is positive semidefinite and $C \in \mathbb { R } ^ { m \times m }$ is symmetric and positive semidefinite is an example of a saddle point problem. Equilibrium systems $( \ S 4 . 4 . 6 )$ are a special case.

Problems with saddle point structure arise in many applications and there is a host of solution frameworks. Various special cases create multiple possibilities for a preconditioner. For example, if A is nonsingular and $C = 0$ , then

$$
\left[ \begin{array}{c c} A & B _ {1} \\ B _ {2} ^ {T} & 0 \end{array} \right] = \left[ \begin{array}{c c} I & 0 \\ B _ {2} ^ {T} A ^ {- 1} & I \end{array} \right] \left[ \begin{array}{c c} A & 0 \\ 0 & S \end{array} \right] \left[ \begin{array}{c c} I & A ^ {- 1} B _ {1} \\ 0 & I \end{array} \right], \qquad S = - B _ {2} ^ {T} A ^ {- 1} B _ {1}.
$$

Possible preconditioners include

$$
M = \left[ \begin{array}{c c} \tilde {A} & 0 \\ 0 & \tilde {S} \end{array} \right] \mathrm{or} \left[ \begin{array}{c c} \tilde {A} & B _ {1} \\ 0 & S \end{array} \right] \mathrm{or} \left[ \begin{array}{c c} \tilde {A} & 0 \\ B _ {2} ^ {T} & \tilde {S} \end{array} \right] \left[ \begin{array}{c c} I & \tilde {A} ^ {- 1} B _ {1} \\ 0 & I \end{array} \right]
$$

where $\tilde { A } \approx A$ and $\tilde { S } \approx S$

If A and C are positive definite, $H _ { 1 } = ( A + A ^ { T } ) / 2 , H _ { 2 } = ( A - A ^ { T } ) / 2$ , and $B = B _ { 1 } = B _ { 2 }$ , then

$$
\left[ \begin{array}{c c} A & B \\ - B ^ {T} & C \end{array} \right] = \left[ \begin{array}{c c} H _ {1} & 0 \\ 0 & C \end{array} \right] + \left[ \begin{array}{c c} H _ {2} & B \\ - B ^ {T} & 0 \end{array} \right] \equiv K _ {1} + K _ {2}
$$

is a symmetric positive definite/skew-symmetric splitting. Preconditioners based on

$$
M = (\alpha I + K _ {2}) ^ {- 1} (\alpha I - K _ {1}) (\alpha I + K _ {1}) ^ {- 1} (\alpha I - K _ {2})
$$

where $\alpha > 0$ have been shown to be effective. See the saddle point problem survey by Benzi, Golub, and Liesen (2005) for more details. Note that the above strategies are specialized ILU strategies.
