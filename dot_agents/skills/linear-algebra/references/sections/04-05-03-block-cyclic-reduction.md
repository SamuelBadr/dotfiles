# 4.5.3 Block-Cyclic Reduction

We next describe the method of block-cyclic reduction that can be used to solve some important special instances of the block tridiagonal system (4.5.1). For simplicity, we assume that A has the form

$$
A = \left[ \begin{array}{c c c c c} D & F & & \dots & 0 \\ F & D & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & F \\ 0 & \dots & & F & D \end{array} \right] \in \mathbb {R} ^ {N q \times N q} \tag {4.5.9}
$$

where $F$ and $D$ are $q { \mathrm { - b y - } } q$ matrices that satisfy ${ D F } = F { D }$ . We also assume that $N = 2 ^ { k } - 1$ . These conditions hold in certain important applications such as the discretization of Poisson’s equation on a rectangle. (See §4.8.4.)

The basic idea behind cyclic reduction is to halve repeatedly the dimension of the problem on hand repeatedly until we are left with a single q-by-q system for the unknown subvector $x _ { 2 } k - 1$ . This system is then solved by standard means. The previously eliminated $x _ { i }$ are found by a back-substitution process.

The general procedure is adequately illustrated by considering the case $N = 7 { : }$

$$
\begin{array}{l} b _ {1} = D x _ {1} + F x _ {2}, \\ b _ {2} = F x _ {1} + D x _ {2} + F x _ {3}, \\ b _ {3} = F x _ {2} + D x _ {3} + F x _ {4}, \\ b _ {4} = F x _ {3} + D x _ {4} + F x _ {5}, \\ b _ {5} = F x _ {4} + D x _ {5} + F x _ {6}, \\ b _ {6} = F x _ {5} + D x _ {6} + F x _ {7}, \\ b _ {7} = F x _ {6} + D x _ {7}. \\ \end{array}
$$

For i = 2, 4, and 6 we multiply equations $i - 1$ , i, and $i + 1$ by $F , \mathrm { ~ } - D _ { \mathrm { { i } } }$ , and $F _ { ; }$ , respectively, and add the resulting equations to obtain

$$
\begin{array}{l} (2 F ^ {2} - D ^ {2}) x _ {2} + F ^ {2} x _ {4} = F (b _ {1} + b _ {3}) - D b _ {2}, \\ F ^ {2} x _ {2} + (2 F ^ {2} - D ^ {2}) x _ {4} + F ^ {2} x _ {6} = F (b _ {3} + b _ {5}) - D b _ {4}, \\ F ^ {2} x _ {4} + (2 F ^ {2} - D ^ {2}) x _ {6} = F (b _ {5} + b _ {7}) - D b _ {6}. \\ \end{array}
$$

Thus, with this tactic we have removed the odd-indexed $x _ { i }$ and are left with a reduced block tridiagonal system of the form

$$
\begin{array}{l} D ^ {(1)} x _ {2} + F ^ {(1)} x _ {4} = b _ {2} ^ {(1)}, \\ F ^ {(1)} x _ {2} + D ^ {(1)} x _ {4} + F ^ {(1)} x _ {6} = b _ {4} ^ {(1)}, \\ F ^ {(1)} x _ {4} + D ^ {(1)} x _ {6} = b _ {6} ^ {(1)}, \\ \end{array}
$$

where $D ^ { ( 1 ) } = 2 F ^ { 2 } - D ^ { 2 }$ and $F ^ { ( 1 ) } = F ^ { 2 }$ commute. Applying the same elimination strategy as above, we multiply these three equations respectively by $F ^ { ( 1 ) } , - D ^ { ( 1 ) }$ , and $F ^ { ( 1 ) }$ . When these transformed equations are added together, we obtain the single equation

$$
\left(2 [ F ^ {(1)} ] ^ {2} - D ^ {(1) ^ {2}}\right) x _ {4} = F ^ {(1)} \left(b _ {2} ^ {(1)} + b _ {6} ^ {(1)}\right) - D ^ {(1)} b _ {4} ^ {(1)},
$$

which we write as

$$
D ^ {(2)} x _ {4} = b ^ {(2)}.
$$

This completes the cyclic reduction. We now solve this (small) q-by-q system for $x _ { 4 }$ . The vectors $x _ { 2 }$ and $x _ { 6 }$ are then found by solving the systems

$$
\begin{array}{l} D ^ {(1)} x _ {2} = b _ {2} ^ {(1)} - F ^ {(1)} x _ {4}, \\ D ^ {(1)} x _ {6} = b _ {6} ^ {(1)} - F ^ {(1)} x _ {4}. \\ \end{array}
$$

Finally, we use the first, third, fifth, and seventh equations in the original system to compute $x _ { 1 } , x _ { 3 } , x _ { 5 }$ , and $x _ { 7 }$ , respectively.

The amount of work required to perform these recursions for general N depends greatly upon the sparsity of the $D ^ { ( p ) }$ and $F ^ { ( p ) }$ . In the worst case when these matrices are full, the overall flop count has order $\mathrm { l o g } ( N ) q ^ { 3 }$ . Care must be exercised in order to ensure stability during the reduction. For further details, see Buneman (1969).
