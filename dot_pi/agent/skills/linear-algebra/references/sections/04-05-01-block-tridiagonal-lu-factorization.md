# 4.5.1 Block Tridiagonal LU Factorization

If

$$
A = \left[ \begin{array}{c c c c c} D _ {1} & F _ {1} & & \dots & 0 \\ E _ {1} & D _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & F _ {N - 1} \\ 0 & \dots & & E _ {N - 1} & D _ {N} \end{array} \right] \tag {4.5.2}
$$

then by comparing blocks in

$$
A = \left[ \begin{array}{c c c c c} I & & & \dots & 0 \\ L _ {1} & I & & & \vdots \\ & \ddots & \ddots & & \\ \vdots & & \ddots & & \\ 0 & \dots & & L _ {N - 1} & I \end{array} \right] \left[ \begin{array}{c c c c c} U _ {1} & F _ {1} & & \dots & 0 \\ 0 & U _ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & F _ {N - 1} \\ 0 & \dots & & 0 & U _ {N} \end{array} \right] \tag {4.5.3}
$$

we formally obtain the following algorithm for computing the $L _ { i }$ and $U _ { i }$ :

$$
U _ {1} = D _ {1}
$$

for i = 2:N

$$
\text { Solve } L _ {i - 1} U _ {i - 1} = E _ {i - 1} \text { for } L _ {i - 1}. \tag {4.5.4}
$$

$$
U _ {i} = D _ {i} - L _ {i - 1} F _ {i - 1}
$$

end

The procedure is defined as long as the $U _ { i }$ are nonsingular.

Having computed the factorization (4.5.3), the vector x in (4.5.1) can be obtained via block forward elimination and block back substitution:

$$
y _ {1} = b _ {1}
$$

for i = 2:N

$$
y _ {i} = b _ {i} - L _ {i - 1} y _ {i - 1}
$$

end (4.5.5)

Solve UN xN = yN for $x _ { N }$

for i = N − 1: − 1:1

Solve Uixi = yi − Fixi+1 for xi

end

To carry out both (4.5.4) and (4.5.5), each $U _ { i }$ must be factored since linear systems involving these submatrices are solved. This could be done using Gaussian elimination with pivoting. However, this does not guarantee the stability of the overall process.
