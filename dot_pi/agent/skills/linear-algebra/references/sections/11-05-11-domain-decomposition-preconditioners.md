# 11.5.11 Domain Decomposition Preconditioners

Domain decomposition is a framework that can be used to design a preconditioner for an Ax = b problem that arises from a discretized boundary value problem (BVP). Here are the main ideas behind the strategy:

Step 1. Express the given “complicated” BVP domain Ω as a union of smaller, “simpler” subdomains $\Omega _ { 1 } , \ldots , \Omega _ { s }$ .

Step 2. Consider what the discretized BVP “looks like” on each subdomain. Presumably, these subproblems are easier to solve because they are smaller and have a computationally friendly geometry.

Step 3. Build the preconditioner M out of the subdomain matrix problems, paying attention to the ordering of the unknowns and how the subdomain solutions relate to one another and the overall solution.

We illustrate this strategy by considering the Poisson problem $\Delta u \ = \ f$ on an Lshaped domain Ω with Dirichlet boundary conditions. (For discretization strategies and solution procedures that are applicable to rectangular domains, see §4.8.4.)

Refer to Figure 11.5.1 where we have subdivided Ω into three non-overlapping rectangular subdomains $\Omega _ { 1 } , \Omega _ { 2 }$ , and $\Omega _ { 3 }$ . As a result of this subdivision, there are five

![](images/golub_650_699__be57e4fabe56fdb5a98c6a6e308d28cfb1abb2adfdcc39fc05457f699f59df80.jpg)

<details>
<summary>text_image</summary>

Grid diagram showing coordinate points and grid points on a boundary, with legend explaining interior and boundary grid points.
</details>

Figure 11.5.1. The Nonoverlapping subdomain framework

“types” of gridpoints (and unknowns). With proper ordering, this leads to a block linear system of the form

$$
A u = \left[ \begin{array}{c c c c c} A _ {1} & 0 & 0 & B & C \\ 0 & A _ {2} & 0 & D & 0 \\ 0 & 0 & A _ {3} & 0 & E \\ F & H & 0 & Q _ {4} & 0 \\ G & 0 & K & 0 & Q _ {5} \end{array} \right] \left[ \begin{array}{l} u _ {\circ^ {1}} \\ u _ {\circ^ {2}} \\ u _ {\circ^ {3}} \\ u _ {\bullet^ {1 2}} \\ u _ {\bullet^ {1 3}} \end{array} \right] = \left[ \begin{array}{l} f _ {\circ^ {1}} \\ f _ {\circ^ {2}} \\ f _ {\circ^ {3}} \\ f _ {\bullet^ {1 2}} \\ f _ {\bullet^ {1 3}} \end{array} \right] = f \tag {11.5.11}
$$

where $A _ { 1 } , A _ { 2 }$ , and $A _ { 3 }$ have the discrete Laplacian structure encountered in §4.8.4. Our notation is intuitive: $u _ { \bullet ^ { 1 2 } }$ is the vector of unknowns associated with the $\bullet ^ { 1 2 }$ grid points. Note that A can be factored as

$$
A = \left[ \begin{array}{c c c c c} I & 0 & 0 & 0 & 0 \\ 0 & I & 0 & 0 & 0 \\ 0 & 0 & I & 0 & 0 \\ F A _ {1} ^ {- 1} & H A _ {2} ^ {- 1} & 0 & I & 0 \\ G A _ {1} ^ {- 1} & 0 & K A _ {3} ^ {- 1} & 0 & I \end{array} \right] \left[ \begin{array}{c c c c c} A _ {1} & 0 & 0 & B & C \\ 0 & A _ {2} & 0 & D & 0 \\ 0 & 0 & A _ {3} & 0 & E \\ 0 & 0 & 0 & S _ {4} & 0 \\ 0 & 0 & 0 & 0 & S _ {5} \end{array} \right] = L U,
$$

where $S _ { 4 }$ and $S _ { 5 }$ are the Schur complements

$$
S _ {4} = Q _ {4} - F A _ {1} ^ {- 1} B - H A _ {2} ^ {- 1} D,
$$

$$
S _ {5} = Q _ {5} - G A _ {1} ^ {- 1} C - K A _ {3} ^ {- 1} E.
$$

If it were not for these typically expensive, dense blocks, the system $A u = f$ could be solved very efficiently via this LU factorization. Fortunately, there are many ways to manage problematic Schur complements. See Saad (IMSLE, pp. 456–465). With appropriate approximations

$$
\tilde {S} _ {4} \approx S _ {4}, \qquad \tilde {S} _ {5} \approx S _ {5},
$$

we are led to a block ILU preconditioner of the form $M = L U _ { M }$ where

$$
U _ {M} = \left[ \begin{array}{c c c c c} A _ {1} & 0 & 0 & B & C \\ 0 & A _ {2} & 0 & D & 0 \\ 0 & 0 & A _ {3} & 0 & E \\ 0 & 0 & 0 & \tilde {S} _ {4} & 0 \\ 0 & 0 & 0 & 0 & \tilde {S} _ {5} \end{array} \right].
$$

With sufficient structure, fast Poisson solvers can be used during the L-solves while the efficiency of the $U _ { M }$ solver would depend upon the nature of the Schur complement approximations.

Although the example is simple, it highlights one of the essential ideas behind nonoverlapping domain decomposition preconditioners like M. Bordered block diagonal systems must be solved where (a) each diagonal block is associated with a subdomain and (b) the border is relatively “thin” because in the partitioning of the overall domain, the number of domain-coupling unknowns is typically an order of magnitude less than the total number of unknowns. A consequence of (b) is that A − M has low rank and this translates into rapid convergence in a Krylov setting. There are also significant opportunities for parallel computation because of the nearly decoupled subdomain computations. See Bjorstad, Gropp, and Smith (1996).

A similar strategy involves overlapping subdomains and we continue with the same example to illustrate the main ideas. Figure 11.5.2 displays a partitioning of the same L-shaped domain into three overlapping subdomains. With proper ordering we obtain

![](images/golub_650_699__c2d8fc1ba01bc18ec6e7b877716e6a4f130ab54ec492b27be911bc72e6e51eca.jpg)

<details>
<summary>text_image</summary>

Grid diagram with labeled points and coordinate labels for interior and boundary points in a 2D grid
</details>

Figure 11.5.2. The overlapping Schwarz framework

a block linear system of the form

$$
A u = \left[ \begin{array}{c c c c c c c} A _ {1} & 0 & 0 & B _ {1} & 0 & C _ {1} & 0 \\ 0 & A _ {2} & 0 & 0 & B _ {2} & 0 & 0 \\ 0 & 0 & A _ {3} & 0 & 0 & 0 & C _ {2} \\ F _ {1} & 0 & 0 & Q _ {4} & D & 0 & 0 \\ 0 & F _ {2} & 0 & H & \tilde {Q} _ {4} & 0 & 0 \\ G _ {1} & 0 & 0 & 0 & 0 & Q _ {5} & E \\ 0 & 0 & G _ {2} & 0 & 0 & K & \tilde {Q} _ {5} \end{array} \right] \left[ \begin{array}{c} u _ {\circ^ {1}} \\ u _ {\circ^ {2}} \\ u _ {\circ^ {3}} \\ u _ {\bullet^ {1 2}} \\ u _ {\bullet^ {2 1}} \\ u _ {\bullet^ {1 3}} \\ u _ {\bullet^ {3 1}} \end{array} \right] = \left[ \begin{array}{c} f _ {\circ^ {1}} \\ f _ {\circ^ {2}} \\ f _ {\circ^ {3}} \\ f _ {\bullet^ {1 2}} \\ f _ {\bullet^ {2 1}} \\ f _ {\bullet^ {1 3}} \\ f _ {\bullet^ {3 1}} \end{array} \right] = f.
$$

In the multiplicative Schwarz approach we cycle through the subdomains improving the interior unknowns along the way. For example, fixing all but the interior $\Omega _ { 1 }$ unknowns, we solve

$$
\left[ \begin{array}{c c c} A _ {1} & B _ {1} & C _ {1} \\ F _ {1} & Q _ {4} & 0 \\ G _ {1} & 0 & Q _ {5} \end{array} \right] \left[ \begin{array}{c} u _ {\circ^ {1}} \\ u _ {\bullet^ {1 2}} \\ u _ {\bullet^ {1 3}} \end{array} \right] = \left[ \begin{array}{c} f _ {\circ^ {1}} \\ f _ {\bullet^ {1 2}} \\ f _ {\bullet^ {1 3}} \end{array} \right] - \left[ \begin{array}{c} 0 \\ D u _ {\bullet^ {2 1}} \\ E u _ {\bullet^ {3 1}} \end{array} \right].
$$

After updating $u _ { \circ ^ { 1 } } , u _ { \bullet ^ { 1 2 } }$ , and $u _ { \bullet ^ { 1 3 } }$ we proceed to fix all but the interior $\Omega _ { 2 }$ unknowns and solve

$$
\left[ \begin{array}{c c} A _ {2} & B _ {2} \\ F _ {2} & \tilde {Q} _ {4} \end{array} \right] \left[ \begin{array}{c} u _ {\circ^ {2}} \\ u _ {\bullet^ {2 1}} \end{array} \right] = \left[ \begin{array}{c} f _ {\circ^ {2}} \\ f _ {\bullet^ {2 1}} \end{array} \right] - \left[ \begin{array}{c} 0 \\ H u _ {\bullet^ {1 2}} \end{array} \right],
$$

and update $u _ { \circ ^ { 2 } }$ and $u _ { \bullet ^ { 2 1 } }$ . Finally, we fix all but the interior $\Omega _ { 3 }$ unknowns and obtain improved versions by solving

$$
\left[ \begin{array}{c c} A _ {3} & C _ {2} \\ G _ {2} & \tilde {Q} _ {5} \end{array} \right] \left[ \begin{array}{c} u _ {\circ^ {3}} \\ u _ {\bullet^ {3 1}} \end{array} \right] = \left[ \begin{array}{c} f _ {\circ^ {3}} \\ f _ {\bullet^ {3 1}} \end{array} \right] - \left[ \begin{array}{c} 0 \\ K u _ {\bullet^ {1 3}} \end{array} \right].
$$

This completes one cycle of multiplicative Schwarz. It is Gauss-Seidel-like in that the most recent values of the current solution are used in each of the three subdomain solves. In the additive Schwarz approach, no part if the solution vector u is updated until after the last subdomain solve. This Jacobi-like approach has certain advantages from the standpoint of parallel computing.

For either the multiplicative or additive approach, it is possible to relate $u ^ { \mathrm { ( n e w ) } }$ to $u ^ { \mathrm { ( o l d ) } }$ via an expression of the form

$$
u ^ {(\mathrm{new})} = u ^ {(\mathrm{old})} + M ^ {- 1} (f - A u ^ {(\mathrm{old})}),
$$

which opens the door to a new family of preconditioning techniques. The geometry of the subdomains and the extent of their overlap critically affects efficiency. Simple geometries can clear a path to fast subdomain solving. Overlap promotes the flow of information between the subdomains but leads to more complicated preconditioners. For an in-depth review of domain decomposition ideas, see Saad (IMSLE, pp. 451–493).

# Problems

P11.5.1 Verify (11.5.2).

P11.5.2 Suppose $H \in \mathbb { R } ^ { n \times n }$ is large sparse upper Hessenberg matrix and that we want to solve $H x = b $ . Note that $H ( [ 2 { : } n 1 ] , : )$ has the form $\boldsymbol { R } + \boldsymbol { e } _ { n } \boldsymbol { v } ^ { T }$ where R is upper triangular and $v \in \mathbb { R } ^ { n }$ . Show how GMRES with preconditioner R can (in principle) be used to solve the system in two iterations.

P11.5.3 Show that

$$
A = \left[ \begin{array}{c c c c} 1 & 1 & 3 & 0 \\ 1 & 2 & 0 & 3 \\ 3 & 0 & 1 9 & - 8 \\ 0 & 3 & - 8 & 1 1 \end{array} \right] = \left[ \begin{array}{c c c c} 1 & 0 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 3 & - 3 & 1 & 0 \\ 0 & 3 & 1 & 1 \end{array} \right] \left[ \begin{array}{c c c c} 1 & 1 & 3 & 0 \\ 0 & 1 & - 3 & 3 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{array} \right]
$$

does not have an incomplete Cholesky factorization if $P = \{ ( 4 , 1 ) , ( 3 , 2 ) \}$ .

P11.5.4 Prove that Equations (11.5.6)–(11.5.8) hold if incChol executes without breakdown.

P11.5.5 Suppose $T \in \mathbb { R } ^ { m \times m }$ is symmetric, tridiagonal, and positive definite. There exist $u , v \in \mathbb { R } ^ { m }$ so that

$$
[ T ^ {- 1} ] _ {i j} = u _ {i} v _ {j}
$$

for all i and j that satisfy $1 \leq j < i < m$ . Give an $O ( m )$ algorithm for computing u and v.

P11.5.6 Suppose $B \in \mathbb { R } ^ { m \times m }$ is a nonsingular, lower bidiagonal matrix. Give an $O ( m )$ algorithm for computing the lower bidiagonal portion of $B ^ { - 1 }$ .

P11.5.7 Consider the computation (11.5.10). Suppose $A _ { 1 } , \dotsc , A _ { p }$ are symmetric with bandwidth q and that $E _ { 1 } , \ldots , E _ { p - 1 }$ have upper bandwidth 0 and lower bandwidth r. What bandwidth constraints on $\Lambda _ { 1 } , \ldots , \Lambda _ { p }$ are necessary if $G _ { 1 } , \ldots , G _ { p }$ are to have lower bandwidth $q ?$

P11.5.8 This problem provides further insight into both the multiplicative Schwarz and additive Schwarz frameworks. Consider the block tridiagonal system

$$
A u = \left[ \begin{array}{c c c} A _ {1 1} & A _ {1 2} & 0 \\ A _ {2 1} & A _ {2 2} & A _ {2 3} \\ 0 & A _ {3 1} & A _ {3 3} \end{array} \right] \left[ \begin{array}{c} u _ {1} \\ u _ {2} \\ u _ {3} \end{array} \right] = \left[ \begin{array}{c} f _ {1} \\ f _ {2} \\ f _ {3} \end{array} \right] = f
$$

where we assume that $A _ { 2 2 }$ is much smaller than either $A _ { 1 1 }$ and $A _ { 3 3 }$ . Assume that an approximate solution $u ^ { ( k ) }$ is improved to $\boldsymbol u ^ { ( k + 1 ) }$ via the following multiplicative Schwarz update procedure:

$$
\left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right] \left[ \begin{array}{l} \Delta_ {1} ^ {(k)} \\ \widetilde {\Delta} _ {2} ^ {(k)} \end{array} \right] = \left[ \begin{array}{l} f _ {1} \\ f _ {2} \end{array} \right] - \left[ \begin{array}{l l l} A _ {1 1} & A _ {1 2} & 0 \\ A _ {2 1} & A _ {2 2} & A _ {2 3} \end{array} \right] \left[ \begin{array}{l} u _ {1} ^ {(k)} \\ u _ {2} ^ {(k)} \\ u _ {3} ^ {(k)} \end{array} \right],
$$

$$
\left[ \begin{array}{c c} A _ {2 2} & A _ {2 3} \\ A _ {3 2} & A _ {3 3} \end{array} \right] \left[ \begin{array}{c} \Delta_ {2} ^ {(k)} \\ \Delta_ {3} ^ {(k)} \end{array} \right] = \left[ \begin{array}{c} f _ {2} \\ f _ {3} \end{array} \right] - \left[ \begin{array}{c c c} A _ {2 1} & A _ {2 2} & A _ {2 3} \\ 0 & A _ {3 2} & A _ {3 3} \end{array} \right] \left[ \begin{array}{c} u _ {1} ^ {(k)} + \Delta_ {1} ^ {(k)} \\ u _ {2} ^ {(k)} + \widetilde {\Delta} _ {2} ^ {(k)} \\ u _ {3} ^ {(k)} \end{array} \right],
$$

$$
\left[ \begin{array}{c} u _ {1} ^ {(k + 1)} \\ u _ {2} ^ {(k + 1)} \\ u _ {3} ^ {(k + 1)} \end{array} \right] = \left[ \begin{array}{c} u _ {1} ^ {(k)} \\ u _ {2} ^ {(k)} \\ u _ {3} ^ {(k)} \end{array} \right] + \left[ \begin{array}{c} \Delta_ {1} ^ {(k)} \\ \Delta_ {2} ^ {(k)} \\ \Delta_ {3} ^ {(k)} \end{array} \right].
$$

(a) Determine a matrix M so that $u ^ { ( k + 1 ) } = u ^ { ( k ) } + M ^ { - 1 } ( f - A u ^ { ( k ) } )$ . (b) Repeat for the additive Schwarz update:

$$
\left[ \begin{array}{l l} A _ {1 1} & A _ {1 2} \\ A _ {2 1} & A _ {2 2} \end{array} \right] \left[ \begin{array}{l} \Delta_ {1} ^ {(k)} \\ \widetilde {\Delta} _ {2} ^ {(k)} \end{array} \right] = \left[ \begin{array}{l} f _ {1} \\ f _ {2} \end{array} \right] - \left[ \begin{array}{l l l} A _ {1 1} & A _ {1 2} & 0 \\ A _ {2 1} & A _ {2 2} & A _ {2 3} \end{array} \right] \left[ \begin{array}{l} u _ {1} ^ {(k)} \\ u _ {2} ^ {(k)} \\ u _ {3} ^ {(k)} \end{array} \right],
$$

$$
\left[ \begin{array}{c c} A _ {2 2} & A _ {2 3} \\ A _ {3 2} & A _ {3 3} \end{array} \right] \left[ \begin{array}{c} \Delta_ {2} ^ {(k)} \\ \Delta_ {3} ^ {(k)} \end{array} \right] = \left[ \begin{array}{c} f _ {2} \\ f _ {3} \end{array} \right] - \left[ \begin{array}{c c c} A _ {2 1} & A _ {2 2} & A _ {2 3} \\ 0 & A _ {3 2} & A _ {3 3} \end{array} \right] \left[ \begin{array}{c} u _ {1} ^ {(k)} \\ u _ {2} ^ {(k)} \\ u _ {3} ^ {(k)} \end{array} \right],
$$

$$
\left[ \begin{array}{c} u _ {1} ^ {(k + 1)} \\ u _ {2} ^ {(k + 1)} \\ u _ {3} ^ {(k + 1)} \end{array} \right] = \left[ \begin{array}{c} u _ {1} ^ {(k)} \\ u _ {2} ^ {(k)} \\ u _ {3} ^ {(k)} \end{array} \right] + \left[ \begin{array}{c} \Delta_ {1} ^ {(k)} \\ \widetilde {\Delta} _ {2} + \Delta_ {2} ^ {(k)} \\ \Delta_ {3} ^ {(k)} \end{array} \right].
$$

For further discussion, see Greenbaum (IMSL, pp. 198–201).

# Notes and References for §11.5

Early papers concerned with preconditioning include:

O. Axelsson (1972). “A Generalized SSOR Method,” BIT 12, 443–467.   
D.J. Evans (1973). “The Analysis and Application of Sparse Matrix Algorithms in the Finite Element Method,” in The Mathematics of Finite Elements and Applications, J.R. Whiteman (ed), Academic Press, New York, 427–447.   
R.H. Bartels and J.W. Daniel (1974). “A Conjugate Gradient Approach to Nonlinear Elliptic Boundary Value Problems,” in Conference on the Numerical Solution of Differential Equations, Dundee, 1973, G.A. Watson (ed), Springer Verlag, New York.   
R.S. Chandra, S.C. Eisenstat, and M.H. Shultz (1975). “Conjugate Gradient Methods for Partial Differential Equations,” in Advances in Computer Methods for Partial Differential Equations, R. Vichnevetsky (ed), Rutgers University, New Brunswick, NJ.   
O. Axelsson (1976). “A Class of Iterative Methods for Finite Element Equations,” Computer Methods in Applied Mechanics and Engineering 9, 123–137.   
P. Concus, G.H. Golub, and D.P. O’Leary (1976). “A Generalized Conjugate Gradient Method for the Numerical Solution of Elliptic Partial Differential Equations,” in Sparse Matrix Computations, J.R. Bunch and D.J. Rose (eds), Academic Press, New York, 309–332.   
J. Douglas Jr. and T. Dupont (1976). “Preconditioned Conjugate Gradient Iteration Applied to Galerkin Methods for a Mildly-Nonlinear Dirichlet Problem,” in Sparse Matrix Computations, J.R. Bunch and D.J. Rose (eds), Academic Press, New York, 333–348.

For an overview of preconditioning techniques, see Greenbaum (IMSL), Meurant (LCG), Saad (IS-PLA), van der Vorst (IMK), LIN TEMPLATES as well as the following surveys:

O. Axelsson (1985). “A Survey of Preconditioned Iterative Methods for Linear Systems of Equations,” BIT 25, 166–187.   
M. Benzi (2002). “Preconditioning for Large Linear Systems: A Survey,” J. Comp. Phys. 182, 418–477.

Papers concerned with sparse approximate inverse preconditioners include:

M. Benzi, C.D. Meyer, and M. Tuma (1996). “A Sparse Approximate Inverse Preconditioner for the Conjugate Gradient Method,” SIAM J. Sci. Comput. 17, 1135–1149.

E. Chow and Y. Saad (1997). “Approximate Inverse Techniques for Block–Partitioned Matrices,” SIAM J. Sci. Comput. 18, 1657–1675.

M.J. Grote and T. Huckle (1997). “Parallel Preconditioning with Sparse Approximate Inverses,” SIAM J. Sci. Comput. 18, 838–853.

N.I.M. Gould and J.A. Scott (1998). “Sparse Approximate-Inverse Preconditioners Using Norm-Minimization Techniques,” SIAM J. Sci. Comput. 19, 605–625.

M. Benzi and M. Tuma (1998). “A Sparse Approximate-Inverse Preconditioner for Nonsymmetric Linear Systems,” SIAM J. Sci. Comput. 19, 968–994.

Various aspects of polynomial preconditioners are discussed in:

O.G. Johnson, C.A. Micchelli, and G. Paul (1983). “Polynomial Preconditioners for Conjugate Gradient Calculations,” SIAM J. Numer. Anal. 20, 362–376.

L. Adams (1985). “m-step Preconditioned Congugate Gradient Methods,” SIAM J. Sci. Stat. Comput. 6, 452–463.

S. Ashby, T. Manteuffel, and P. Saylor (1989). “Adaptive Polynomial Preconditioning for Hermitian Indefinite Linear Systems,” BIT 29, 583–609.

R.W. Freund (1990). “On Conjugate Gradient Type Methods and Polynomial Preconditioners for a Class of Complex Non-Hermitian Matrices,” Numer. Math. 57, 285–312.

S. Ashby, T. Manteuffel, and J. Otto (1992). “A Comparison of Adaptive Chebyshev and Least Squares Polynomial Preconditioning for Hermitian Positive Definite Linear Systems,” SIAM J. Sci. Stat. Comput. 13, 1–29.

The incomplete Cholesky factorization idea is set forth and analyzed in:

J.A. Meijerink and H.A. van der Vorst (1977). “An Iterative Solution Method for Linear Equation Systems of Which the Coefficient Matrix is a Symmetric M-Matrix,” Math. Comput. 31, 148–162.

T.A. Manteuffel (1980). “An Incomplete Factorization Technique for Positive Definite Linear Systems,” Math. Comput. 34, 473–497.

C.-J. Lin and J.J. Mor´e (1999). “Incomplete Cholesky Factorizations with Limited Memory,” SIAM J. Sci. Comput. 21, 24–45.

Likewise, for the incomplete LU factorization strategy we have:

M. Bollhofer and Y. Saad (2006). “Multilevel Preconditioners Constructed From Inverse-Based ILUs,” SIAM J. Sci. Comput. 27, 1627–1650.

H. Elman (1986). “A Stability Analysis of Incomplete LU Factorization,” Math. Comput. 47, 191–218.

Incomplete QR factorizations have also been devised. See Bj¨orck (NMLS, pp. 297–299) as well as:

Z. Jia (1998). “On IOM(q): The Incomplete Orthogonalization Method for Large Unsymmetric Linear Systems,” Numer. Lin. Alg. 3, 491–512.

Z.-Z. Bai, I.S. Duff, and A.J. Wathen (2001). “A Class of Incomplete Orthogonal Factorization Methods. I: Methods and Theories,” BIT 41, 53–70.

Incomplete block factorizations are discussed in:

G. Roderigue and D. Wolitzer (1984). “Preconditioning by Incomplete Block Cyclic Reduction,” Math. Comput. 42, 549–566.

P. Concus, G.H. Golub, and G. Meurant (1985). “Block Preconditioning for the Conjugate Gradient Method,” SIAM J. Sci. Stat. Comput. 6, 220–252.

O. Axelsson (1985). “Incomplete Block Matrix Factorization Preconditioning Methods. The Ultimate Answer?”, J. Comput. Appl. Math. 12–13, 3–18.

O. Axelsson (1986). “A General Incomplete Block Matrix Factorization Method,” Lin. Alg. Applic. 74, 179–190.

The analysis of incomplete factorizations is both difficult and important, see:   
Y. Notay (1992). “On the Robustness of Modified Incomplete Factorization Methods,” J. Comput. Math. 40, 121–141.   
H. Lu and O. Axelsson (1997). “Conditioning Analysis of Block Incomplete Factorizations and Its Application to Elliptic Equations,” Numer. Math. 78, 189–209.   
M. Bollhofer and Y. Saad (2002). “On the Relations between ILUs and Factored Approximate Inverses,” SIAM J. Matrix Anal. Applic. 24, 219–237.   
Numerous vector/parallel implementations of the preconditioned CG method have been developed, see:   
G. Meurant (1984). “The Block Preconditioned Conjugate Gradient Method on Vector Computers,” BIT 24, 623–633.   
C.C. Ashcraft and R. Grimes (1988). “On Vectorizing Incomplete Factorization and SSOR Preconditioners,” SIAM J. Sci. Stat. Comp. 9, 122–151.   
U. Meier and A. Sameh (1988). “The Behavior of Conjugate Gradient Algorithms on a Multivector Processor with a Hierarchical Memory,” J. Comput. Appl. Math. 24, 13–32.   
H. van der Vorst (1989). “High Performance Preconditioning,” SIAM J. Sci. Stat. Comput. 10, 1174–1185.   
V. Eijkhout (1991). “Analysis of Parallel Incomplete Point Factorizations,” Lin. Alg. Applic. 154– 156, 723–740.   
Preconditioners for large Toeplitz systems are discussed in:   
T.F. Chan (1988). “An Optimal Circulant Preconditioner for Toeplitz Systems,” SIAM. J. Sci. Stat. Comput. 9, 766–771.   
R.H. Chan and G. Strang (1989). “Toeplitz Equations by Conjugate Gradients with Circulant Preconditioner,” SIAM J. Sci. Stat. Comput. 10, 104–119.   
T. Huckle (1992). “Circulant and Skew-circulant Matrices for Solving Toeplitz Matrix Problems,” SIAM J. Matrix Anal. Applic. 13, 767–777.   
R.H. Chan, J.G. Nagy, and R.J. Plemmons (1994). “Circulant Preconditioned Toeplitz Least Squares Iterations,” SIAM J. Matrix Anal. Applic. 15, 80–97.   
T.F. Chan and J.A. Olkin (1994). “Circulant Preconditioners for Toeplitz Block Matrices,” Numer. alg. 6, 89–101.   
R.H. Chan and M.K. Ng (1996). “Conjugate Gradient Methods for Toeplitz Systems,” SIAM Review 38, 427–482.   
R.H. Chan and X.-Q. Jin (2007). An Introduction to Iterative Toeplitz Solvers, SIAM Publications, Philadelphia, PA.   
Preconditioners based on the splitting of a matrix into the sum of its symmetric and skew-symmetric parts is covered in the following papers:   
Z.-Z. Bai, G.H. Golub, and M.K. Ng (2003). “Hermitian and Skew-Hermitian Splitting Methods for Non-Hermitian Positive Definite Linear Systems,” SIAM J. Matrix Anal. Applic. 24, 603–626.   
Z.-Z. Bai, G.H. Golub, and J.-Y. Pan (2004). “Preconditioned Hermitian and Skew-Hermitian Splitting Methods for Non-Hermitian Positive Semidefinite Linear Systems,” Numer. Math. 98, 1–32.   
Z.-Z. Bai, G.H. Golub, L.-Z. Lu, and J.-F. Yin (2005). “Block Triangular and Skew-Hermitian Splitting Methods for Positive-Definite Linear Systems,” SIAM J. Sci. Comput. 26, 844–863.   
For a discussion of saddle point systems and their preconditioning, see:   
M. Benzi, G.H. Golub, and J. Liesen (2005). “Numerical Solution of Saddle Point Problems,” Acta Numerica 14, 1–137.   
G.H. Golub, C. Greif, and J.M. Varah (2005). “An Algebraic Analysis of a Block Diagonal Preconditioner for Saddle Point Systems,” SIAM J. Matrix Anal. Applic. 27, 779–792.   
H.S. Dollar, N.I.M. Gould, W.H.A. Schilders, and A.J. Wathen (2006). “Implicit-Factorization Preconditioning and Iterative Solvers for Regularized Saddle-Point Systems,” SIAM J. Matrix Anal. Applic. 28, 170–189.   
C. Greif and D. Schtzau (2006). “Preconditioners for Saddle Point Linear Systems with Highly Singular (1,1) Blocks,” ETNA 22, 114–121.   
M.A. Botchev and G.H. Golub (2006). “A Class of Nonsymmetric Preconditioners for Saddle Point Problems,” SIAM J. Matrix Anal. Applic. 27, 1125–1149.

The handling of problematic Schur complements has attracted much attention. For an appreciation of the challenge and what to do about it, see:

H. Elman (1989). “Approximate Schur Complement Preconditioners on Serial and Parallel Computers,” SIAM J. Sci. Stat. Comput. 10, 581–605.   
S.C. Brenner (1999). “The Condition Number of the Schur Complement in Domain Decomposition,” Numer. Math. 83, 187–203.   
F. Zhang (2005). The Schur Complement and its Applications, Springer-Verlag, New York.   
Z. Li and Y. Saad (2006). “SchurRAS: A Restricted Version of the Overlapping Schur Complement Preconditioner,” SIAM J. Sci. Comput. 27, 1787–1801.

For an overview of the domain decomposition paradigm, see Demmel (ANLA, pp. 347–356) as well as:

T.F. Chan and T.P. Mathew (1994). “Domain Decomposition Algorithms,” Acta Numerica 3, 61–143.   
W.D. Gropp and D.E. Keyes (1992). “Domain Decomposition with Local Mesh Refinement,” SIAM J. Sci. Statist. Comput. 13, 967–993.   
D.E. Keyes, T.F. Chan, G. Meurant, J.S. Scroggs, and R.G. Voigt (eds) (1992). Domain Decomposition Methods for Partial Differential Equations, SIAM Publications, Philadelphia, PA.   
T.F. Chan and D. Goovaerts (1992). “On the Relationship Between Overlapping and Nonoverlapping Domain Decomposition Methods,” SIAM J. Matrix Anal. Applic. 13, 663–670.   
B. Smith, P. Bjorstad, and W. Gropp (1996). Domain Decomposition–Parallel Multilevel Methods for Elliptic Partial Differential Equations, Cambridge University Press, Cambridge, U.K.   
J. Xu and J. Xou (1998) “Some Nonoverlapping Domain Decomposition Methods,” SIAM Review 40, 857–914.   
A. Tosseli and O. Widlund (2010). Domain Decomposition Methods: Theory and Algorithms, Springer-Verlag, New York.

For insight into the role of preconditioning for least squares problems and more generally in numerical optimization, see:

P.E. Gill, W. Murray, D.B. Poncele´on, and M.A. Saunders (1992). “Preconditioners for Indefinite Systems Arising in Optimization,” SIAM J. Matrix Anal. Applic. 13, 292–311.   
A. Bj¨orck and J. Y. Yuan (1999). “Preconditioners for Least Squares Problems by LU Factorization,” ETNA 8, 26–35.   
M. Benzi and M. Tuma (2003). “A Robust Preconditioner with Low Memory Requirements for Large Sparse Least Squares Problems,” SIAM J. Sci. Comput. 25, 499–512.   
M. Jacobsen, P. C. Hansen, and M. A. Saunders (2003). “Subspace Preconditioned LSQR for Discrete Ill-Posed Problems,” BIT 43, 975–989.   
O. Axelsson and M. Neytcheva (2003). “Preconditioning Methods for Linear Systems Arising in Constrained Optimization Problems,” Numer. Lin. Alg. Applic. 10, 3–31.   
A.R.L. Oliveira and D.C. Sorensen (2004). “A New Class of Preconditioners for Large-Scale Linear Systems from Interior Point Methods for Linear Programming,” Lin. Alg. Applic. 394, 1–24.

Other ideas associated with preconditioning include inexact solution of the preconditioned system $M z = r$ and variation of M from iteration to iteration, see:

J. Baglama, D. Calvetti, G. H. Golub, and L. Reichel (1998). “Adaptively Preconditioned GMRES Algorithms,” SIAM J. Sci. Comput. 20, 243–269.   
G.H. Golub and Q.Ye (1999). “Inexact Preconditioned Conjugate Gradient Method with Inner-Outer Iteration,” SIAM J. Sci. Comput. 21, 1305–1320.   
Y. Notay (2000). “Flexible Conjugate Gradients,” SIAM J. Sci. Comput. 22, 1444–1460.

Error estimation in the preconditioned CG context is discussed in:

O. Axelsson and I. Kaporin (2001). “Error Norm Estimation and Stopping Criteria in Preconditioned Conjugate Gradient Iterations,” Numer. Lin. Alg. 8, 265–286.   
Z. Strakos and P. Tichy (2005). “Error Estimation in Preconditioned Conjugate Gradients,” BIT 45, 789–817.
