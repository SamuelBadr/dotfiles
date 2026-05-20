# 7.7.6 The QZ Step

We are now in a position to describe a QZ step. The basic idea is to update A and B as follows

$$
(\bar {A} - \lambda \bar {B}) = \bar {Q} ^ {T} (A - \lambda B) \bar {Z},
$$

where $\bar { A }$ is upper Hessenberg, $\bar { B }$ is upper triangular, $\bar { Q }$ and $\bar { Z }$ are each orthogonal, and $\bar { A } \bar { B } ^ { - 1 }$ is essentially the same matrix that would result if a Francis QR step (Algorithm 7.5.1) were explicitly applied to $A B ^ { - 1 }$ . This can be done with some clever zero-chasing and an appeal to the implicit Q theorem.

Let $M = A B ^ { - 1 }$ (upper Hessenberg) and let v be the first column of the matrix $( M - a I ) ( M - b I )$ , where a and b are the eigenvalues of $M \mathrm { { s } }$ lower 2-by-2 submatrix. Note that v can be calculated in $O ( 1 )$ flops. If $P _ { 0 }$ is a Householder matrix such that $P _ { 0 } v$ is a multiple of $e _ { 1 }$ , then

$$
A \leftarrow P _ {0} A = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \leftarrow P _ {0} B = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

The idea now is to restore these matrices to Hessenberg-triangular form by chasing the unwanted nonzero elements down the diagonal.

To this end, we first determine a pair of Householder matrices $Z _ { 1 }$ and $Z _ { 2 }$ to zero $b _ { 3 1 } , \ b _ { 3 2 }$ , and $b _ { 2 1 }$ :

$$
A \leftarrow A Z _ {1} Z _ {2} = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \leftarrow B Z _ {1} Z _ {2} = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

Then a Householder matrix $P _ { 1 }$ is used to zero $a _ { 3 1 }$ and $a _ { 4 1 }$ :

$$
A \gets P _ {1} A = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \end{array} \right], \quad B \gets P _ {1} B = \left[ \begin{array}{c c c c c c} \times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times \end{array} \right].
$$

Notice that with this step the unwanted nonzero elements have been shifted down and to the right from their original position. This illustrates a typical step in the QZ iteration. Notice that $Q = Q _ { 0 } Q _ { 1 } \cdot \cdot \cdot Q _ { n - 2 }$ has the same first column as $Q _ { 0 }$ . By the way the initial Householder matrix was determined, we can apply the implicit Q theorem and assert that $A B ^ { - 1 } = Q ^ { T } ( A B ^ { - 1 } ) Q$ is indeed essentially the same matrix that we would obtain by applying the Francis iteration to $M = A B ^ { - 1 }$ directly. Overall we have the following algorithm.

Algorithm 7.7.2 (The QZ Step) Given an unreduced upper Hessenberg matrix $A \in \mathbb { R } ^ { n \times n }$ and a nonsingular upper triangular matrix $B \in \mathbb { R } ^ { n \times n }$ , the following algorithm overwrites A with the upper Hessenberg matrix $Q ^ { T } A Z$ and $B$ with the upper triangular matrix $Q ^ { T } B Z$ where $Q$ and $Z$ are orthogonal and $Q$ has the same first column as the orthogonal similarity transformation in Algorithm 7.5.1 when it is applied to $A B ^ { - 1 }$ .

Let $M = A B ^ { - 1 }$ and compute $( M - a I ) ( M - b I ) e _ { 1 } = [ x , y , z , 0 , \ldots , 0 ] ^ { T }$

where a and b are the eigenvalues of $M \mathrm { { s } }$ lower 2-by-2.

for $k = 1 { : } n - 2$

$\mathrm { F i n d ~ H o u s e h o l d e r } Q _ { k } \mathrm { \ s o } Q _ { k } \left[ \begin{array} { l } { x } \\ { y } \\ { z } \end{array} \right] = \left[ \begin{array} { l } { * } \\ { 0 } \\ { 0 } \end{array} \right] .$

$$
A = \operatorname{diag} (I _ {k - 1}, Q _ {k}, I _ {n - k - 2}) \cdot A
$$

$$
B = \operatorname{diag} (I _ {k - 1}, Q _ {k}, I _ {n - k - 2}) \cdot B
$$

$\mathrm { F i n d ~ H o u s e h o l d e r ~ } Z _ { k 1 } \mathrm { ~ s o ~ } \left[ \begin{array} { l } { b _ { k + 2 , k } \ \big \vert \ b _ { k + 2 , k + 1 } \ \big \vert \ b _ { k + 2 , k + 2 } } \end{array} \right] Z _ { k 1 } = \left[ \begin{array} { l }  0 \ \big \vert \ 0 \ \big \vert \ \ast \ \right] . \end{array}$

$$
A = A \cdot \operatorname{diag} \left(I _ {k - 1}, Z _ {k 1}, I _ {n - k - 2}\right)
$$

$$
B = B \cdot \operatorname{diag} \left(I _ {k - 1}, Z _ {k 1}, I _ {n - k - 2}\right)
$$

$\mathrm { F i n d ~ H o u s e h o l d e r ~ } Z _ { k 2 } \mathrm { ~ s o ~ \ } \left[ \begin{array} { l } { b _ { k + 1 , k } \ | \ b _ { k + 1 , k + 1 } } \end{array} \right] Z _ { k 2 } = \left[ \begin{array} { l } { 0 \ | \ * \ } \end{array} \right] .$

$$
A = A \cdot \operatorname{diag} \left(I _ {k - 1}, Z _ {k 2}, I _ {n - k - 1}\right)
$$

$$
B = B \cdot \operatorname{diag} \left(I _ {k - 1}, Z _ {k 2}, I _ {n - k - 1}\right)
$$

$$
x = a _ {k + 1, k}; y = a _ {k + 2, k}
$$

$\mathbf { i f } \ k < n - 2$

$$
z = a _ {k + 3, k}
$$

end

end

Find Householder Qn−1 so Qn−1  xy 	 =  ∗0 	 .

$$
A = \operatorname{diag} (I _ {n - 2}, Q _ {n - 1}) \cdot A
$$

$$
B = \operatorname{diag} (I _ {n - 2}, Q _ {n - 1}) \cdot B.
$$

$\mathrm { F i n d ~ H o u s e h o l d e r ~ } Z _ { n - 1 } \mathrm { ~ s o ~ } \left[ \textit { b } _ { n , n - 1 } \mid \textit { b } _ { n n } \right] Z _ { n - 1 } = \left[ \textit { 0 } \mid * \right] .$

$$
A = A \cdot \operatorname{diag} \left(I _ {n - 2}, Z _ {n - 1}\right)
$$

$$
B = B \cdot \mathrm{diag} (I _ {n - 2}, Z _ {n - 1})
$$

This algorithm requires $2 2 n ^ { 2 }$ flops. $Q$ and $Z$ can be accumulated for an additional $8 n ^ { 2 }$ flops and $1 3 n ^ { 2 }$ flops, respectively.
