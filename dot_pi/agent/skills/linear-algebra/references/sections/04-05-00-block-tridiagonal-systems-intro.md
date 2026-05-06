# 4.5 Block Tridiagonal Systems

Block tridiagonal linear systems of the form

$$
\left[ \begin{array}{c c c c c} D _ {1} & F _ {1} & & \dots & 0 \\ E _ {1} & D _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & F _ {N - 1} \\ 0 & \dots & & E _ {N - 1} & D _ {N} \end{array} \right] \left[ \begin{array}{c} x _ {1} \\ x _ {2} \\ \vdots \\ \vdots \\ x _ {N} \end{array} \right] = \left[ \begin{array}{c} b _ {1} \\ b _ {2} \\ \vdots \\ \vdots \\ b _ {N} \end{array} \right]. \tag {4.5.1}
$$

frequently arise in practice. We assume for clarity that all blocks are $q { \mathrm { - b y - } } q$ . In this section we discuss both a block LU approach to this problem as well as a pair of divide-and-conquer schemes.
