# 12.2.2 Quasiseparable Matrices

Certain products of Givens rotations exhibit rank structure, but we frame the key fact in more general terms. If $\alpha , \beta , \gamma , \delta \in \mathbb { R } ^ { n - 1 }$ and

$$
M _ {k} = \operatorname{diag} (I _ {k - 1}, \tilde {M} _ {k}, I _ {n - k - 1}), \qquad \tilde {M} _ {k} = \left[ \begin{array}{c c} \alpha_ {k} & \beta_ {k} \\ \gamma_ {k} & \delta_ {k} \end{array} \right],
$$

for $k = 1 { : } n - 1$ , then the matrix $M = M _ { 1 } \cdot \cdot \cdot M _ { n - 1 }$ is fully illustrated by

$$
M = M _ {1} M _ {2} M _ {3} M _ {4} = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} \alpha_ {2} & \beta_ {1} \beta_ {2} \alpha_ {3} & \beta_ {1} \beta_ {2} \beta_ {3} \alpha_ {4} & \beta_ {1} \beta_ {2} \beta_ {3} \beta_ {4} \\ \gamma_ {1} & \delta_ {1} \alpha_ {2} & \delta_ {1} \beta_ {2} \alpha_ {3} & \delta_ {1} \beta_ {2} \beta_ {3} \alpha_ {4} & \delta_ {1} \beta_ {2} \beta_ {3} \beta_ {4} \\ 0 & \gamma_ {2} & \delta_ {2} \alpha_ {3} & \delta_ {2} \beta_ {3} \alpha_ {4} & \delta_ {2} \beta_ {3} \beta_ {4} \\ 0 & 0 & \gamma_ {3} & \delta_ {3} \alpha_ {4} & \delta_ {3} \beta_ {4} \\ 0 & 0 & 0 & \gamma_ {4} & \delta_ {4} \end{array} \right]. \tag {12.2.4}
$$

It has the property that off-diagonal blocks have unit rank or less provided they do not “intersect” the diagonal. Quasiseparable matrices have this property and if A is such a matrix, then

$$
j _ {2} <   i _ {1} \text {   or   } i _ {2} <   j _ {1} \Rightarrow \operatorname{rank} (A (i _ {1}: i _ {2}, j _ {1}: j _ {2})) \leq 1. \tag {12.2.5}
$$

By comparing this with (12.2.1), it is clear that the class of semiseparable matrices is a subset of the class of quasiseparable matrices.
