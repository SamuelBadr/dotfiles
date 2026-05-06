# 4.4 Symmetric Indefinite Systems

Recall that a matrix whose quadratic form $x ^ { T }$ Ax takes on both positive and negative values is indefinite. In this section we are concerned with symmetric indefinite linear systems. The $\mathrm { L D L ^ { T } }$ factorization is not always advisable as the following 2-by-2 example illustrates:

$$
\left[ \begin{array}{c c} \epsilon & 1 \\ 1 & 0 \end{array} \right] = \left[ \begin{array}{c c} 1 & 0 \\ 1 / \epsilon & 1 \end{array} \right] \left[ \begin{array}{c c} \epsilon & 0 \\ 0 & - 1 / \epsilon \end{array} \right] \left[ \begin{array}{c c} 1 & 0 \\ 1 / \epsilon & 1 \end{array} \right] ^ {T}.
$$

Of course, any of the pivot strategies in §3.4 could be invoked. However, they destroy symmetry and, with it, the chance for a “Cholesky speed” symmetric indefinite system solver. Symmetric pivoting, i.e., data reshufflings of the form $A  P A P ^ { T }$ , must be used as we discussed in §4.2.8. Unfortunately, symmetric pivoting does not always stabilize the $\mathrm { L D L ^ { T } }$ computation. If $\epsilon _ { 1 }$ and $\epsilon _ { 2 }$ are small, then regardless of $P _ { \cdot }$ the matrix

$$
\tilde {A} = P \left[ \begin{array}{c c} \epsilon_ {1} & 1 \\ 1 & \epsilon_ {2} \end{array} \right] P ^ {T}
$$

has small diagonal entries and large numbers surface in the factorization. With symmetric pivoting, the pivots are always selected from the diagonal and trouble results if these numbers are small relative to what must be zeroed off the diagonal. Thus, $\mathrm { L D L ^ { T } }$ with symmetric pivoting cannot be recommended as a reliable approach to symmetric indefinite system solving. It seems that the challenge is to involve the off-diagonal entries in the pivoting process while at the same time maintaining symmetry.

In this section we discuss two ways to do this. The first method is due to Aasen (1971) and it computes the factorization

$$
P A P ^ {T} = L T L ^ {T}, \tag {4.4.1}
$$

where $L = ( \ell _ { i j } )$ is unit lower triangular and $T$ is tridiagonal. P is a permutation chosen such that $| \ell _ { i j } | \le 1$ . In contrast, the diagonal pivoting method due to Bunch and Parlett (1971) computes a permutation P such that

$$
P A P ^ {T} = L D L ^ {T}, \tag {4.4.2}
$$

where D is a direct sum of 1-by-1 and 2-by-2 pivot blocks. Again, $P$ is chosen so that the entries in the unit lower triangular L satisfy $| \ell _ { i j } | \le 1$ . Both factorizations involve $n ^ { 3 } / 3$ flops and once computed, can be used to solve $A x = b$ with $O ( n ^ { 2 } )$ work:

$$
P A P ^ {T} = L T L ^ {T}, L z = P b, T w = z, L ^ {T} y = w, x = P ^ {T} y \Rightarrow A x = b,
$$

$$
P A P ^ {T} = L D L ^ {T}, L z = P b, D w = z, L ^ {T} y = w, x = P ^ {T} y \Rightarrow A x = b.
$$

A few comments need to be made about the $T w = z$ and $D w = z$ systems that arise when these methods are invoked.

In Aasen’s method, the symmetric indefinite tridiagonal system $T w = z$ is solved in $O ( n )$ time using band Gaussian elimination with pivoting. Note that there is no serious price to pay for the disregard of symmetry at this level since the overall process is $O ( n ^ { 3 } )$ .

In the diagonal pivoting approach, the $D w = z$ system amounts to a set of 1-by-1 and $2 \mathrm { - b y { - } 2 }$ symmetric indefinite systems. The 2-by-2 problems can be handled via Gaussian elimination with pivoting. Again, there is no harm in disregarding symmetry during this $O ( n )$ phase of the calculation. Thus, the central issue in this section is the efficient computation of the factorizations (4.4.1) and (4.4.2).

# 4.4.1 The Parlett-Reid Algorithm

Parlett and Reid (1970) show how to compute (4.4.1) using Gauss transforms. Their algorithm is sufficiently illustrated by displaying the k = 2 step for the case $n = 5$ . At the beginning of this step the matrix A has been transformed to

$$
A ^ {(1)} = M _ {1} P _ {1} A P _ {1} ^ {T} M _ {1} ^ {T} = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & 0 & 0 & 0 \\ \beta_ {1} & \alpha_ {2} & v _ {3} & v _ {4} & v _ {5} \\ 0 & v _ {3} & \times & \times & \times \\ 0 & v _ {4} & \times & \times & \times \\ 0 & v _ {5} & \times & \times & \times \end{array} \right],
$$

where $P _ { 1 }$ is a permutation chosen so that the entries in the Gauss transformation $M _ { 1 }$ are bounded by unity in modulus. Scanning the vector $[ v _ { 3 } \ v _ { 4 } \ v _ { 5 } ] ^ { T }$ for its largest entry, we now determine a 3-by-3 permutation ${ \tilde { P } } _ { 2 }$ such that

$$
\tilde {P} _ {2} \left[ \begin{array}{l} v _ {3} \\ v _ {4} \\ v _ {5} \end{array} \right] = \left[ \begin{array}{l} \tilde {v} _ {3} \\ \tilde {v} _ {4} \\ \tilde {v} _ {5} \end{array} \right] \qquad \Rightarrow \qquad | \tilde {v} _ {3} | = \max \{| \tilde {v} _ {3} |, | \tilde {v} _ {4} |, | \tilde {v} _ {5} | \}.
$$

If this maximal element is zero, we set $M _ { 2 } = P _ { 2 } = I$ and proceed to the next step. Otherwise, we set $P _ { 2 } = \mathrm { d i a g } ( I _ { 2 } , \tilde { P } _ { 2 } )$ and $M _ { 2 } = I \ - \ \alpha ^ { ( 2 ) } e _ { 3 } ^ { T }$ with

$$
\alpha^ {(2)} = \left[ \begin{array}{l l l l l} 0 & 0 & 0 & \tilde {v} _ {4} / \tilde {v} _ {3} & \tilde {v} _ {5} / \tilde {v} _ {3} \end{array} \right] ^ {T}.
$$

Observe that

$$
A ^ {(2)} = M _ {2} P _ {2} A ^ {(1)} P _ {2} ^ {T} M _ {2} ^ {T} = \left[ \begin{array}{l l l l l} \alpha_ {1} & \beta_ {1} & 0 & 0 & 0 \\ \beta_ {1} & \alpha_ {2} & \tilde {v} _ {3} & 0 & 0 \\ 0 & \tilde {v} _ {3} & \times & \times & \times \\ 0 & 0 & \times & \times & \times \\ 0 & 0 & \times & \times & \times \end{array} \right].
$$

In general, the process continues for $n - 2$ steps leaving us with a tridiagonal matrix

$$
T = A ^ {(n - 2)} = (M _ {n - 2} P _ {n - 2} \dots M _ {1} P _ {1}) A (M _ {n - 2} P _ {n - 2} \dots M _ {1} P _ {1}) ^ {T}.
$$

It can be shown that (4.4.1) holds with $P = P _ { n - 2 } \cdot \cdot \cdot P _ { 1 }$ and

$$
L = \left(M _ {n - 2} P _ {n - 2} \dots M _ {1} P _ {1} P ^ {T}\right) ^ {- 1}.
$$

Analysis of L reveals that its first column is $e _ { 1 }$ and that its subdiagonal entries in column k with $k > 1$ are “made up” of the multipliers in $M _ { k - 1 }$ .

The efficient implementation of the Parlett-Reid method requires care when computing the update

$$
A ^ {(k)} = M _ {k} (P _ {k} A ^ {(k - 1)} P _ {k} ^ {T}) M _ {k} ^ {T}. \tag {4.4.3}
$$

To see what is involved with a minimum of notation, suppose B = BT ∈ IR(n−k)×(n−k) $B = B ^ { T } \in \mathbb { R } ^ { ( n - k ) \times ( n - k ) }$ has and that we wish to form

$$
B _ {+} = (I - w e _ {1} ^ {T}) B (I - w e _ {1} ^ {T}) ^ {T},
$$

where $w \in \mathbb { R } ^ { n - k }$ and $e _ { 1 }$ is the first column of $I _ { n - k }$ . Such a calculation is at the heart of (4.4.3). If we set

$$
u = B e _ {1} - \frac {b _ {1 1}}{2} w,
$$

then $\begin{array} { r } { B _ { + } = B - w u ^ { T } - u w ^ { T } } \end{array}$ and its lower triangular portion can be formed in $2 ( n - k ) ^ { 2 }$ flops. Summing this quantity as k ranges from 1 to $n - 2$ indicates that the Parlett-Reid procedure requires $2 n ^ { 3 } / 3$ flops—twice the volume of work associated with Cholesky.

# 4.4.2 The Method of Aasen

An $n ^ { 3 } / 3$ approach to computing (4.4.1) due to Aasen (1971) can be derived by reconsidering some of the computations in the Parlett-Reid approach. We examine the no-pivoting case first where the goal is to compute a unit lower triangular matrix L with $L ( : , 1 ) = e _ { 1 }$ and a tridiagonal matrix

$$
T = \left[ \begin{array}{c c c c c} \alpha_ {1} & \beta_ {1} & & \dots & 0 \\ \beta_ {1} & \alpha_ {2} & \ddots & & \vdots \\ & \ddots & \ddots & \ddots & \\ \vdots & & \ddots & \ddots & \beta_ {n - 1} \\ 0 & \dots & & \beta_ {n - 1} & \alpha_ {n} \end{array} \right].
$$

such that $A = L T L ^ { T }$ . The Aasen method is structured as follows:

for $j = 1:n$ $\{\alpha(1:j - 1), \beta(1:j - 1) \text{ and } L(:,1:j) \text{ are known}\}$ Compute $\alpha_j$ .  
    if $j \leq n - 1$ Compute $\beta_j$ .  
    end  
    if $j \leq n - 2$ Compute $L(j + 2:n, j + 1)$ .  
    end  
end

To develop recipes for $\alpha _ { j } , \beta _ { j }$ , and $L ( j + 2 ; n , j + 1 )$ , we compare the jth columns in the equation $A = L H$ where $\boldsymbol { H } ^ { \setminus } = \boldsymbol { T } \boldsymbol { L } ^ { T }$ . Noting that H is an upper Hessenberg matrix we obtain

$$
A (:, j) = L H (:, j) = \sum_ {k = 1} ^ {j + 1} L (:, k) \cdot h (k), \tag {4.4.5}
$$

where $h ( 1 { : } j + 1 ) = H ( 1 { : } j + 1 , j )$ and we assume that $j \le n - 1$ . It follows that

$$
h _ {j + 1} \cdot L (j + 1: n, j + 1) = v (j + 1: n), \tag {4.4.6}
$$

where

$$
v (j + 1: n) = A (j + 1: n, j) - L (j + 1: n, 1: j) \cdot h (1: j). \tag {4.4.7}
$$

Since L is unit lower triangular and $L ( : , 1 : j )$ is known, this gives us a working recipe for $L ( j + 2 { : } n , j + 1 )$ provided we know $h ( 1 { : } j )$ . Indeed, from (4.4.6) and (4.4.7) it is easy to show that

$$
L (j + 2: n, j + 1) = v (j + 2: n) / v (j + 1). \tag {4.4.8}
$$

To compute $h ( 1 { : } j )$ we turn to the equation $H = T L ^ { T }$ and examine its jth column. The case $j = 5$ amply displays what is going on:

$$
\left[ \begin{array}{l} h _ {1} \\ h _ {2} \\ h _ {3} \\ h _ {4} \\ h _ {5} \\ h _ {6} \end{array} \right] = \left[ \begin{array}{l l l l l} \alpha_ {1} & \beta_ {1} & 0 & 0 & 0 \\ \beta_ {1} & \alpha_ {2} & \beta_ {2} & 0 & 0 \\ 0 & \beta_ {2} & \alpha_ {3} & \beta_ {3} & 0 \\ 0 & 0 & \beta_ {3} & \alpha_ {4} & \beta_ {4} \\ 0 & 0 & 0 & \beta_ {4} & \alpha_ {5} \\ 0 & 0 & 0 & 0 & \beta_ {5} \end{array} \right] \left[ \begin{array}{c} 0 \\ \ell_ {5 2} \\ \ell_ {5 3} \\ \ell_ {5 4} \\ 1 \end{array} \right] = \left[ \begin{array}{c} \beta_ {1} \ell_ {5 2} \\ \alpha_ {2} \ell_ {5 2} + \beta_ {2} \ell_ {5 3} \\ \beta_ {2} \ell_ {5 2} + \alpha_ {3} \ell_ {5 3} + \beta_ {3} \ell_ {5 4} \\ \beta_ {3} \ell_ {5 3} + \alpha_ {4} \ell_ {5 4} + \beta_ {4} \\ \beta_ {4} \ell_ {5 4} + \alpha_ {5} \\ \beta_ {5} \end{array} \right] \tag {4.4.9}
$$

At the start of step $j ,$ we know $\alpha ( 1 { : } j - 1 ) , \beta ( 1 { : } j - 1 )$ and $L ( : , 1 : j )$ . Thus, we can determine $h ( 1 { : } j - 1 )$ as follows

$$
h _ {1} = \beta_ {1} \ell_ {j 2}
$$

for $k = 1 { : } j - 1$

$$
h _ {k} = \beta_ {k - 1} \ell_ {j, k - 1} + \alpha_ {k} \ell_ {j k} + \beta_ {k} \ell_ {j, k + 1} \tag {4.4.10}
$$

end

Equation (4.4.5) gives us a formula for $h _ { j }$ :

$$
h _ {j} = A (j, j) - \sum_ {k = 1} ^ {j - 1} L (j, k) h _ {k}. \tag {4.4.11}
$$

From (4.4.9) we infer that

$$
\alpha_ {j} = h _ {j} - \beta_ {j - 1} \ell_ {j, j - 1}, \tag {4.4.12}
$$

$$
\beta_ {j} = h _ {j + 1}. \tag {4.4.13}
$$

Combining these equations with (4.4.4), (4.4.7), (4.4.8), (4.4.10), and (4.4.11) we obtain the Aasen method without pivoting:

$$
L = I _ {n}
$$

for j = 1:n

$$
\text { if } j = 1
$$

$$
\alpha_ {1} = a _ {1 1}
$$

$$
v (2: n) = A (2: n, 1)
$$

else

$$
h _ {1} = \beta_ {1} \cdot \ell_ {j 2}
$$

for k = 2:j − 1

$$
h _ {k} = \beta_ {k - 1} \ell_ {j, k - 1} + \alpha_ {k} \ell_ {j k} + \beta_ {k} \ell_ {j, k + 1}
$$

end

$$
h _ {j} = a _ {j j} - L (j, 1: j - 1) \cdot h (1: j - 1)
$$

$$
\alpha_ {j} = h _ {j} - \beta_ {j - 1} \ell_ {j, j - 1} \tag {4.4.14}
$$

$$
v (j + 1: n) = A (j + 1: n, j) - L (j + 1: n, 1: j) \cdot h (1: j)
$$

end

$$
\text { if   } j <   = n - 1
$$

$$
\beta_ {j} = v (j + 1)
$$

end

$$
\text { if } j <   = n - 2
$$

$$
L (j + 2: n, j + 1) = v (j + 2: n) / v (j + 1)
$$

end

end

The dominant operation each pass through the j-loop is an (n−j)-by-j gaxpy operation. Accounting for the associated flops we see that the overall Aasen ccomputation involves $n ^ { 3 } / 3$ flops, the same as for the Cholesky factorization.

As it now stands, the columns of L are scalings of the v-vectors in (4.4.14). If any of these scalings are large, i.e., if any $v ( j + 1 )$ is small, then we are in trouble. To circumvent this problem, it is only necessary to permute the largest component of $v ( j + 1 { : } n )$ to the top position. Of course, this permutation must be suitably applied to the unreduced portion of A and the previously computed portion of L. With pivoting, Aasen’s method is stable in the same sense that Gaussian elimination with partial pivoting is stable.

In a practical implementation of the Aasen algorithm, the lower triangular portion of A would be overwritten with L and T , e.g.,

$$
A \leftarrow \left[ \begin{array}{c c c c c} \alpha_ {1} & & & & \\ \beta_ {1} & \alpha_ {2} & & & \\ \ell_ {3 2} & \beta_ {2} & \alpha_ {3} & & \\ \ell_ {4 2} & \ell_ {4 3} & \beta_ {3} & \alpha_ {4} & \\ \ell_ {5 2} & \ell_ {5 3} & \ell_ {5 4} & \beta_ {4} & \alpha_ {5} \end{array} \right].
$$

Notice that the columns of L are shifted left in this arrangement.

# 4.4.3 Diagonal Pivoting Methods

We next describe the computation of the block $L D L ^ { T }$ factorization (4.4.2). We follow the discussion in Bunch and Parlett (1971). Suppose

$$
P _ {1} A P _ {1} ^ {T} = \left[ \begin{array}{c c} E & C ^ {T} \\ C & B \end{array} \right] _ {n - s} ^ {s}
$$

where $P _ { 1 }$ is a permutation matrix and s = 1 or 2. If A is nonzero, then it is always possible to choose these quantities so that E is nonsingular, thereby enabling us to write

$$
P _ {1} A P _ {1} ^ {T} = \left[ \begin{array}{c c} I _ {s} & 0 \\ C E ^ {- 1} & I _ {n - s} \end{array} \right] \left[ \begin{array}{c c} E & 0 \\ 0 & B - C E ^ {- 1} C ^ {T} \end{array} \right] \left[ \begin{array}{c c} I _ {s} & E ^ {- 1} C ^ {T} \\ 0 & I _ {n - s} \end{array} \right].
$$

For the sake of stability, the s-by-s “pivot” E should be chosen so that the entries in

$$
\tilde {A} = \left(\tilde {a} _ {i j}\right) \equiv B - C E ^ {- 1} C ^ {T} \tag {4.4.15}
$$

are suitably bounded. To this end, let $\alpha \in ( 0 , 1 )$ be given and define the size measures

$$
\begin{array}{l} \mu_ {0} = \max _ {i, j} | a _ {i j} |, \\ \mu_ {1} = \max _ {i} | a _ {i i} |. \\ \end{array}
$$

The Bunch-Parlett pivot strategy is as follows:

$$
\begin{array}{l} \mathbf {i f} \mu_ {1} \geq \alpha \mu_ {0} \\ s = 1 \\ \text { Choose } P _ {1} \text { so } | e _ {1 1} | = \mu_ {1}. \\ s = 2 \\ \text { Choose } P _ {1} \text { so } | e _ {2 1} | = \mu_ {0}. \\ \end{array}
$$

It is easy to verify from (4.4.15) that if s = 1, then

$$
\left| \tilde {a} _ {i j} \right| \leq \left(1 + \alpha^ {- 1}\right) \mu_ {0}, \tag {4.4.16}
$$

while s = 2 implies

$$
\left| \tilde {a} _ {i j} \right| \leq \frac {3 - \alpha}{1 - \alpha} \mu_ {0}. \tag {4.4.17}
$$

By equating $( 1 + \alpha ^ { - 1 } ) ^ { 2 }$ , the growth factor that is associated with two $s = 1$ steps, and $( 3 - \alpha ) / ( 1 - \alpha )$ , the corresponding s = 2 factor, Bunch and Parlett conclude that $\alpha = ( 1 + { \sqrt { 1 7 } } ) / 8$ is optimum from the standpoint of minimizing the bound on element growth.

The reductions outlined above can be repeated on the order- $( n - s )$ symmetric matrix ${ \tilde { A } } .$ . A simple induction argument establishes that the factorization (4.4.2) exists and that $n ^ { 3 } / 3$ flops are required if the work associated with pivot determination is ignored.

# 4.4.4 Stability and Efficiency

Diagonal pivoting with the above strategy is shown by Bunch (1971) to be as stable as Gaussian elimination with complete pivoting. Unfortunately, the overall process requires between $n ^ { 3 } / 1 2$ and $n ^ { 3 } / 6$ comparisons, since $\mu _ { 0 }$ involves a two-dimensional search at each stage of the reduction. The actual number of comparisons depends on the total number of 2-by-2 pivots but in general the Bunch-Parlett method for computing (4.4.2) is considerably slower than the technique of Aasen. See Barwell and George (1976).

This is not the case with the diagonal pivoting method of Bunch and Kaufman (1977). In their scheme, it is only necessary to scan two columns at each stage of the reduction. The strategy is fully illustrated by considering the very first step in the reduction:

$\alpha = (1 + \sqrt{17}) / 8$ $\lambda = |a_{r1}| = \max \{|a_{21}|, \ldots, |a_{n1}|\}$ if $\lambda > 0$ if $|a_{11}| \geq \alpha \lambda$ Set $s = 1$ and $P_{1} = I$ .   
else $\sigma = |a_{pr}| = \max \{|a_{1r}, \ldots, |a_{r-1,r}|, |a_{r+1,r}|, \ldots, |a_{nr}|\}$ if $\sigma |a_{11}| \geq \alpha \lambda^2$ Set $s = 1$ and $P_{1} = I$ elseif $|a_{rr}| \geq \alpha \sigma$ Set $s = 1$ and choose $P_{1}$ so $(P_{1}^{T}AP_{1})_{11} = a_{rr}$ .   
else   
Set $s = 2$ and choose $P_{1}$ so $(P_{1}^{T}AP_{1})_{21} = a_{rp}$ .   
end   
end   
end

Overall, the Bunch-Kaufman algorithm requires $n ^ { 3 } / 3$ flops, $O ( n ^ { 2 } )$ comparisons, and, like all the methods of this section, $n ^ { 2 } / 2$ storage.

# 4.4.5 A Note on Equilibrium Systems

A very important class of symmetric indefinite matrices have the form

$$
A = \left[ \begin{array}{c c} C & B \\ B ^ {T} & 0 \\ n & p \end{array} \right] _ {p} ^ {n} \tag {4.4.18}
$$

where C is symmetric positive definite and B has full column rank. These conditions ensure that A is nonsingular.

Of course, the methods of this section apply to A. However, they do not exploit its structure because the pivot strategies “wipe out” the zero (2,2) block. On the other hand, here is a tempting approach that does exploit A’s block structure:

Step 1. Compute the Cholesky factorization $C = G G ^ { T }$ .

Step 2. Solve GK = B for $K \in \mathbb { R } ^ { n \times p }$ .

Step 3. Compute the Cholesky factorization $H H ^ { T } = K ^ { T } K = B ^ { T } C ^ { - 1 } B .$

From this it follows that

$$
A = \left[ \begin{array}{c c} G & 0 \\ K ^ {T} & H \end{array} \right] \left[ \begin{array}{c c} G ^ {T} & K \\ 0 & - H ^ {T} \end{array} \right].
$$

In principle, this triangular factorization can be used to solve the equilibrium system

$$
\left[ \begin{array}{c c} C & B \\ B ^ {T} & 0 \end{array} \right] \left[ \begin{array}{l} x \\ y \end{array} \right] = \left[ \begin{array}{l} f \\ g \end{array} \right]. \tag {4.4.19}
$$

However, it is clear by considering steps (b) and (c) above that the accuracy of the computed solution depends upon $\kappa ( C )$ and this quantity may be much greater than $\kappa ( A )$ . The situation has been carefully analyzed and various structure-exploiting algorithms have been proposed. A brief review of the literature is given at the end of the section.

It is interesting to consider a special case of (4.4.19) that clarifies what it means for an algorithm to be stable and illustrates how perturbation analysis can structure the search for better methods. In several important applications, $g = 0 , C$ is diagonal, and the solution subvector y is of primary importance. A manipulation shows that this vector is specified by

$$
y = (B ^ {T} C ^ {- 1} B) ^ {- 1} B ^ {T} C ^ {- 1} f. \tag {4.4.20}
$$

Looking at this we are again led to believe that $\kappa ( C )$ should have a bearing on the accuracy of the computed y. However, it can be shown that

$$
\left\| \left(B ^ {T} C ^ {- 1} B\right) ^ {- 1} B ^ {T} C ^ {- 1} \right\| \leq \psi_ {B} \tag {4.4.21}
$$

where the upper bound $\psi _ { B }$ is independent of C, a result that (correctly) suggests that y is not sensitive to perturbations in C. A stable method for computing this vector should respect this, meaning that the accuracy of the computed y should be independent of C. Vavasis (1994) has developed a method with this property. It involves the careful assembly of a matrix $V \in \mathbb { R } ^ { n \times ( n - p ) }$ whose columns are a basis for the nullspace of $B ^ { T } C ^ { - 1 }$ . The n-by-n linear system

$$
[ B \mid V ] \left[ \begin{array}{l} y \\ q \end{array} \right] = f
$$

is then solved implying $f = B y + V q$ . Thus, $B ^ { T } C ^ { - 1 } f = B ^ { T } C ^ { - 1 } B y$ and (4.4.20) holds.

# Problems

P4.4.1 Show that if all the 1-by-1 and 2-by-2 principal submatrices of an n-by-n symmetric matrix A are singular, then A is zero.

P4.4.2 Show that no 2-by-2 pivots can arise in the Bunch-Kaufman algorithm if A is positive definite.

P4.4.3 Arrange (4.4.14) so that only the lower triangular portion of A is referenced and so that $\alpha ( j )$ overwrites $A ( j , j )$ for $j = 1 { : } n , \beta ( j )$ overwrites $A ( j + 1 , j )$ for $j = 1 { : } n - 1$ , and $L ( i , j )$ overwrites $A ( i , j - 1 )$ for $j = 2 { : } n - 1$ and $i = j + 1 { : } n$ .

P4.4.4 Suppose $A \in \mathbb { R } ^ { n \times n }$ is symmetric and strictly diagonally dominant. Give an algorithm that computes the factorization

$$
\Pi A \Pi^ {T} = \left[ \begin{array}{c c} R & 0 \\ S & - M \end{array} \right] \left[ \begin{array}{c c} R ^ {T} & S ^ {T} \\ 0 & M ^ {T} \end{array} \right]
$$

where Π is a permuation and the diagonal blocks R and M are lower triangular.

P4.4.5 A symmetric matrix A is quasidefinite if it has the form

$$
A = \left[ \begin{array}{c c} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & - A _ {2 2} \end{array} \right] _ {p} ^ {n}
$$

with $A _ { 1 1 }$ and $A _ { 2 2 }$ positive definite. (a) Show that such a matrix has an $\mathrm { L D L } ^ { T }$ factorization with the property that

$$
D = \left[ \begin{array}{c c} D _ {1} & 0 \\ 0 & - D _ {2} \end{array} \right]
$$

where $D _ { 1 } \in \mathbb { R } ^ { n \times n }$ and $D _ { 2 } \in \mathbb { R } ^ { p \times p }$ have positive diagonal entries. (b) Show that if A is quasidefinite then all its principal submatrices are nonsingular. This means that ${ \bf \vec { P } } { \bf \vec { { A } } } { \cal { P } } ^ { T }$ has an $\mathrm { L D L } ^ { T }$ factorization for any permutation matrix P .

P4.4.6 Prove (4.4.16) and (4.4.17).

P4.4.7 Show that $- ( B ^ { T } C ^ { - 1 } B ) ^ { - 1 }$ is the (2,2) block of $A ^ { - 1 }$ where A is given by equation (4.4.18).

P4.4.8 The point of this problem is to consider a special case of (4.4.21). Define the matrix

$$
M (\alpha) = (B ^ {T} C ^ {- 1} B) ^ {- 1} B ^ {T} C ^ {- 1}
$$

where $C = ( I _ { n } + \alpha e _ { k } e _ { k } ^ { T } )$ , α > −1, and $e _ { k } = I _ { n } ( : , k )$ . (Note that C is just the identity with α added to the (k, k) entry.) Assume that $B \in \mathbb { R } ^ { n \times p }$ has rank p and show that

$$
M (\alpha) = (B ^ {T} B) ^ {- 1} B ^ {T} \left(I _ {n} - \frac {\alpha}{1 + \alpha w ^ {T} w} e _ {k} w ^ {T}\right)
$$

where

$$
w = (I _ {n} - B (B ^ {T} B) ^ {- 1} B ^ {T}) e _ {k}.
$$

Show that if $\parallel w \parallel _ { 2 } = 0 \mathrm { o r } \parallel w \parallel _ { 2 } = 1$ , then $\parallel M ( \alpha ) \parallel _ { 2 } = 1 / \sigma _ { \mathrm { m i n } } ( B )$ . Show that if $0 < \parallel w \parallel _ { 2 } < 1$ , then

$$
\| M (\alpha) \| _ {2} \leq \max \left\{\frac {1}{1 - \| w \| _ {2}}, 1 + \frac {1}{\| w \| _ {2}} \right\} \Bigg / \sigma_ {\min} (B).
$$

Thus,  $M ( \alpha ) \parallel _ { 2 }$ has an α-independent upper bound.

# Notes and References for §4.4

The basic references for computing (4.4.1) are as follows:

J.O. Aasen (1971). “On the Reduction of a Symmetric Matrix to Tridiagonal Form,” BIT 11, 233–242.

B.N. Parlett and J.K. Reid (1970). “On the Solution of a System of Linear Equations Whose Matrix Is Symmetric but not Definite,”BIT 10, 386–397.

The diagonal pivoting literature includes:

J.R. Bunch and B.N. Parlett (1971). “Direct Methods for Solving Symmetric Indefinite Systems of Linear Equations,” SIAM J. Numer. Anal. 8, 639–655.

J.R. Bunch (1971). “Analysis of the Diagonal Pivoting Method,” SIAM J. Numer. Anal. 8, 656–680.

J.R. Bunch (1974). “Partial Pivoting Strategies for Symmetric Matrices,” SIAM J. Numer. Anal. 11, 521–528.

J.R. Bunch, L. Kaufman, and B.N. Parlett (1976). “Decomposition of a Symmetric Matrix,” Numer. Math. 27, 95–109.   
J.R. Bunch and L. Kaufman (1977). “Some Stable Methods for Calculating Inertia and Solving Symmetric Linear Systems,” Math. Comput. 31, 162–79.   
M.T. Jones and M.L. Patrick (1993). “Bunch-Kaufman Factorization for Real Symmetric Indefinite Banded Matrices,” SIAM J. Matrix Anal. Applic. 14, 553–559.   
Because “future” columns must be scanned in the pivoting process, it is awkward (but possible) to obtain a gaxpy-rich diagonal pivoting algorithm. On the other hand, Aasen’s method is naturally rich in gaxpys. Block versions of both procedures are possible. Various performance issues are discussed in:   
V. Barwell and J.A. George (1976). “A Comparison of Algorithms for Solving Symmetric Indefinite Systems of Linear Equations,” ACM Trans. Math. Softw. 2, 242–251.   
M.T. Jones and M.L. Patrick (1994). “Factoring Symmetric Indefinite Matrices on High-Performance Architectures,” SIAM J. Matrix Anal. Applic. 15, 273–283.   
Another idea for a cheap pivoting strategy utilizes error bounds based on more liberal interchange criteria, an idea borrowed from some work done in the area of sparse elimination methods, see:   
R. Fletcher (1976). “Factorizing Symmetric Indefinite Matrices,” Lin. Alg. Applic. 14, 257–272.   
Before using any symmetric Ax = b solver, it may be advisable to equilibrate A. An $O ( n ^ { 2 } )$ algorithm for accomplishing this task is given in:   
J.R. Bunch (1971). “Equilibration of Symmetric Matrices in the Max-Norm,”J. ACM 18, 566–572.   
N.J. Higham (1997). “Stability of the Diagonal Pivoting Method with Partial Pivoting,” SIAM J. Matrix Anal. Applic. 18, 52–65.   
Procedures for skew-symmetric systems similar to the methods that we have presented in this section also exist:   
J.R. Bunch (1982). “A Note on the Stable Decomposition of Skew Symmetric Matrices,” Math. Comput. 158, 475–480.   
J. Bunch (1982). “Stable Decomposition of Skew-Symmetric Matrices,” Math. Comput. 38, 475–479.   
P. Benner, R. Byers, H. Fassbender, V. Mehrmann, and D. Watkins (2000). “Cholesky-like Factorizations of Skew-Symmetric Matrices,” ETNA 11, 85–93.   
For a discussion of symmetric indefinite system solvers that are also banded or sparse, see:   
C. Ashcraft, R.G. Grimes, and J.G. Lewis (1998). “Accurate Symmetric Indefinite Linear Equation Solvers,” SIAM J. Matrix Anal. Applic. 20, 513–561.   
S.H. Cheng and N.J. Higham (1998). “A Modified Cholesky Algorithm Based on a Symmetric Indefinite Factorization,” SIAM J. Matrix Anal. Applic. 19, 1097–1110.   
J. Zhao, W. Wang, and W. Ren (2004). “Stability of the Matrix Factorization for Solving Block Tridiagonal Symmetric Indefinite Linear Systems,” BIT 44, 181–188.   
H. Fang and D.P. O’Leary (2006). “Stable Factorizations of Symmetric Tridiagonal and Triadic Matrices,” SIAM J. Matrix Anal. Applic. 28, 576–595.   
D. Irony and S. Toledo (2006). “The Snap-Back Pivoting Method for Symmetric Banded Indefinite Matrices,” SIAM J. Matrix Anal. Applic. 28, 398–424.

The equilibrium system literature is scattered among the several application areas where it has an important role to play. Nice overviews with pointers to this literature include:

G. Strang (1988). “A Framework for Equilibrium Equations,” SIAM Review 30, 283–297.

S.A. Vavasis (1994). “Stable Numerical Algorithms for Equilibrium Systems,” SIAM J. Matrix Anal. Applic. 15, 1108–1131.

P.E. Gill, M.A. Saunders, and J.R. Shinnerl (1996). “On the Stability of Cholesky Factorization for Symmetric Quasidefinite Systems,” SIAM J. Matrix Anal. Applic. 17, 35–46.

G.H. Golub and C. Greif (2003). “On Solving Block-Structured Indefinite Linear Systems,” SIAM J. Sci. Comput. 24, 2076–2092.

For a discussion of (4.4.21), see:

G.W. Stewart (1989). “On Scaled Projections and Pseudoinverses,” Lin. Alg. Applic. 112, 189–193.

D.P. O’Leary (1990). “On Bounds for Scaled Projections and Pseudoinverses,” Lin. Alg. Applic. 132, 115–117.

M.J. Todd (1990). “A Dantzig-Wolfe-like Variant of Karmarkar’s Interior-Point Linear Programming Algorithm,” Oper. Res. 38, 1006–1018.

An equilibrium system is a special case of a saddle point system. See §11.5.10.
