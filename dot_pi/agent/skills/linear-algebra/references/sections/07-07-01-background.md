# 7.7.1 Background

The first thing to observe about the generalized eigenvalue problem is that there are n eigenvalues if and only if rank $( B ) = n$ . If B is rank deficient then $\lambda ( A , B )$ may be finite, empty, or infinite:

$$
A   =   \left[ \begin{array}{c c} 1 & 2 \\ 0 & 3 \end{array} \right],    B   =   \left[ \begin{array}{c c} 1 & 0 \\ 0 & 0 \end{array} \right]    \Rightarrow    \lambda (A, B) = \{1 \},
$$

$$
A   =   \left[ \begin{array}{c c} 1 & 2 \\ 0 & 3 \end{array} \right], \quad B   =   \left[ \begin{array}{c c} 0 & 1 \\ 0 & 0 \end{array} \right] \quad \Rightarrow \quad \lambda (A, B) = \emptyset ,
$$

$$
A   =   \left[ \begin{array}{c c} 1 & 2 \\ 0 & 0 \end{array} \right],    B   =   \left[ \begin{array}{c c} 1 & 0 \\ 0 & 0 \end{array} \right]    \Rightarrow    \lambda (A, B) = \mathbb {C}.
$$

Note that if $0 \neq \lambda \in \lambda ( A , B )$ , then $( 1 / \lambda ) \in \lambda ( B , A )$ . Moreover, if B is nonsingular, then $\lambda ( A , B ) = \lambda ( B ^ { - 1 } A , I ) = \lambda ( B ^ { - 1 } A )$ . This last observation suggests one method for solving the $A - \lambda B$ problem if B is nonsingular:

Step 1. Solve $B C = A$ for C using (say) Gaussian elimination with pivoting.

Step 2. Use the QR algorithm to compute the eigenvalues of C.

In this framework, C is affected by roundoff errors of order $\mathbf { u } \parallel A \parallel _ { 2 } \parallel B ^ { - 1 } \parallel _ { 2 }$ . If B is illconditioned, then this precludes the possibility of computing any generalized eigenvalue accurately—even those eigenvalues that may be regarded as well-conditioned. For example, if

$$
A = \left[ \begin{array}{l l} 1. 7 4 6 & . 9 4 0 \\ 1. 2 4 6 & 1. 8 9 8 \end{array} \right] \qquad \text {and} \qquad B = \left[ \begin{array}{l l}. 7 8 0 & . 5 6 3 \\ . 9 1 3 & . 6 5 9 \end{array} \right],
$$

then $\lambda ( A , B ) = \{ 2 , 1 . 0 7 \times 1 0 ^ { 6 } \}$ . With 7-digit floating point arithmetic, we find $\lambda ( \mathrm { f l } ( A B ^ { - 1 } ) ) = \{ 1 . 5 6 2 5 3 9 , 1 . 0 1 \times 1 0 ^ { 6 } \}$ . The poor quality of the small eigenvalue is because $\kappa _ { 2 } ( B ) \approx 2 \times 1 0 ^ { 6 }$ . On the other hand, we find that

$$
\lambda (I, \mathsf {f l} (A ^ {- 1} B)) \approx \{2. 0 0 0 0 0 1, 1. 0 6 \times 1 0 ^ {6} \}.
$$

The accuracy of the small eigenvalue is improved because $\kappa _ { 2 } ( A ) \approx 4$

The example suggests that we seek an alternative approach to the generalized eigenvalue problem. One idea is to compute well-conditioned Q and Z such that the matrices

$$
A _ {1} = Q ^ {- 1} A Z, \quad B _ {1} = Q ^ {- 1} B Z \tag {7.7.2}
$$

are each in canonical form. Note that $\lambda ( A , B ) { = } \lambda ( A _ { 1 } , B _ { 1 } )$ since

$$
A x = \lambda B x \Leftrightarrow A _ {1} y = \lambda B _ {1} y, x = Z y.
$$

We say that the pencils A − λB and $A _ { 1 } - \lambda B _ { 1 }$ are equivalent if (7.7.2) holds with nonsingular Q and Z.

As in the standard eigenproblem A − λI there is a choice between canonical forms. Corresponding to the Jordan form is a decomposition of Kronecker in which both $A _ { 1 }$ and $B _ { 1 }$ are block diagonal with blocks that are similar in structure to Jordan blocks. The Kronecker canonical form poses the same numerical challenges as the Jordan form, but it provides insight into the mathematical properties of the pencil $A - \lambda B$ . See Wilkinson (1978) and Demmel and K˚agstr¨om (1987) for details.
