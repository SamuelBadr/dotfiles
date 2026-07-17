# MATRIX

#

Gene H.Golub

# Matrix Computations

Johns Hopkins Studies in the Mathematical Sciences in association with the Department of Mathematical Sciences, The Johns Hopkins University

# Matrix Computations

Fourth Edition

Gene H. Golub

Department of Computer Science

Stanford University

Charles F. Van Loan

Department of Computer Science

Cornell University

The Johns Hopkins University Press

Baltimore

c 1983, 1989, 1996, 2013 The Johns Hopkins University Press

All rights reserved. Published 2013

Printed in the United States of America on acid-free paper

9 8 7 6 5 4 3 2 1

First edition 1983

Second edition 1989

Third edition 1996

Fourth edition 2013

The Johns Hopkins University Press

2715 North Charles Street

Baltimore, Maryland 21218-4363

www.press.jhu.edu

Library of Congress Control Number: 2012943449

ISBN 13: 978-1-4214-0794-4 (hc)

ISBN 10: 1-4214-0794-9 (hc)

ISBN 13: 978-1-4214-0859-0 (eb)

ISBN 10: 1-4214-0859-7 (eb)

A catalog record for this book is available from the British Library.

MATLABr is a registered trademark of The Mathworks Inc.

Special discounts are available for bulk purchases of this book. For more information, please contact Special Sales at 410-516-6936 or specialsales@press.jhu.edu.

The Johns Hopkins University Press uses environmentally friendly book materials, including recycled text paper that is composed of at least 30 percent post-consumer waste, whenever possible.

To

# ALSTON S.HOUSEHOLDER

and

JAMES H. WILKINSON

This page intentionally left blank

# Contents

Preface xi

Global References xiii

Other Books xv

Useful URLs xix

Common Notation xxi

# 1 Matrix Multiplication 1

1.1 Basic Algorithms and Notation 2

1.2 Structure and Efficiency 14

1.3 Block Matrices and Algorithms 2 2

1.4 Fast Matrix-Vector Products 3 3

1.5 Vectorization and Locality 4 3

1.6 Parallel Matrix Multiplication 4 9

# 2 Matrix Analysis 63

2.1 Basic Ideas from Linear Algebra 6 4

2.2 Vector Norms 68

2.3 Matrix Norms 71

2.4 The Singular Value Decomposition 7 6

2.5 Subspace Metrics 8 1

2.6 The Sensitivity of Square Systems 8 7

2.7 Finite Precision Matrix Computations 9 3

# 3 General Linear Systems 105

3.1 Triangular Systems 106

3.2 The LU Factorization 111

3.3 Roundoff Error in Gaussian Elimination 122

3.4 Pivoting 125

3.5 Improving and Estimating Accuracy 137

3.6 Parallel LU 144

# 4 Special Linear Systems 153

4.1 Diagonal Dominance and Symmetry 154

4.2 Positive Definite Systems 159

4.3 Banded Systems 176   
4.4 Symmetric Indefinite Systems 186   
4.5 Block Tridiagonal Systems 196   
4.6 Vandermonde Systems 203   
4.7 Classical Methods for Toeplitz Systems 208   
4.8 Circulant and Discrete Poisson Systems 219

# 5 Orthogonalization and Least Squares 233

5.1 Householder and Givens Transformations 234   
5.2 The QR Factorization 246   
5.3 The Full-Rank Least Squares Problem 260   
5.4 Other Orthogonal Factorizations 274   
5.5 The Rank-Deficient Least Squares Problem 288   
5.6 Square and Underdetermined Systems 298

# 6 Modified Least Squares Problems and Methods 303

6.1 Weighting and Regularization 304   
6.2 Constrained Least Squares 313   
6.3 Total Least Squares 320   
6.4 Subspace Computations with the SVD 327   
6.5 Updating Matrix Factorizations 334

# 7 Unsymmetric Eigenvalue Problems 347

7.1 Properties and Decompositions 348   
7.2 Perturbation Theory 357   
7.3 Power Iterations 365   
7.4 The Hessenberg and Real Schur Forms 376   
7.5 The Practical QR Algorithm 385   
7.6 Invariant Subspace Computations 394   
7.7 The Generalized Eigenvalue Problem 405   
7.8 Hamiltonian and Product Eigenvalue Problems 420   
7.9 Pseudospectra 426

# 8 Symmetric Eigenvalue Problems 439

8.1 Properties and Decompositions 440   
8.2 Power Iterations 450   
8.3 The Symmetric QR Algorithm 458   
8.4 More Methods for Tridiagonal Problems 467   
8.5 Jacobi Methods 476   
8.6 Computing the SVD 486   
8.7 Generalized Eigenvalue Problems with Symmetry 497

# 9 Functions of Matrices 513

9.1 Eigenvalue Methods 514

9.2 Approximation Methods 522

9.3 The Matrix Exponential 530

9.4 The Sign, Square Root, and Log of a Matrix 536

# 10 Large Sparse Eigenvalue Problems 545

10.1 The Symmetric Lanczos Process 546

10.2 Lanczos, Quadrature, and Approximation 556

10.3 Practical Lanczos Procedures 562

10.4 Large Sparse SVD Frameworks 571

10.5 Krylov Methods for Unsymmetric Problems 579

10.6 Jacobi-Davidson and Related Methods 589

# 11 Large Sparse Linear System Problems 597

11.1 Direct Methods 598

11.2 The Classical Iterations 611

11.3 The Conjugate Gradient Method 625

11.4 Other Krylov Methods 639

11.5 Preconditioning 650

11.6 The Multigrid Framework 670

# 12 Special Topics 681

12.1 Linear Systems with Displacement Structure 681

12.2 Structured-Rank Problems 691

12.3 Kronecker Product Computations 707

12.4 Tensor Unfoldings and Contractions 719

12.5 Tensor Decompositions and Iterations 731

Index 747

This page intentionally left blank

# Preface

My thirty-year book collaboration with Gene Golub began in 1977 at a matrix computation workshop held at Johns Hopkins University. His interest in my work at the start of my academic career prompted the writing of GVL1. Sadly, Gene died on November 16, 2007. At the time we had only just begun to talk about GVL4. While writing these pages, I was reminded every day of his far-reaching impact and professional generosity. This edition is a way to thank Gene for our collaboration and the friendly research community that his unique personality helped create.

It has been sixteen years since the publication of the third edition—a power-of-two reminder that what we need to know about matrix computations is growing exponentially! Naturally, it is impossible to provide in-depth coverage of all the great new advances and research trends. However, with the relatively recent publication of so many excellent textbooks and specialized volumes, we are able to complement our brief treatments with useful pointers to the literature. That said, here are the new features of GVL4:

# Content

The book is about twenty-five percent longer. There are new sections on fast transforms (§1.4), parallel LU (§3.6), fast methods for circulant systems and discrete Poisson systems (§4.8), Hamiltonian and product eigenvalue problems (§7.8), pseudospectra (§7.9), the matrix sign, square root, and logarithm functions (§9.4), Lanczos and quadrature (§10.2), large-scale SVD (§10.4), Jacobi-Davidson (§10.6), sparse direct methods (§11.1), multigrid (§11.6), low displacement rank systems (§12.1), structuredrank systems (§12.2), Kronecker product problems (§12.3), tensor contractions (§12.4), and tensor decompositions (§12.5).

New topics at the subsection level include recursive block LU (§3.2.11), rook pivoting (§3.4.7), tournament pivoting (§3.6.3), diagonal dominance (§4.1.1), recursive block structures (§4.2.10), band matrix inverse properties (§4.3.8), divide-and-conquer strategies for block tridiagonal systems (§4.5.4), the cross product and various point/plane least squares problems (§5.3.9), the polynomial eigenvalue problem (§7.7.9), and the structured quadratic eigenvalue problem (§8.7.9).

Substantial upgrades include our treatment of floating-point arithmetic (§2.7), LU roundoff error analysis (§3.3.1), LS sensitivity analysis (§5.3.6), the generalized singular value decomposition (§6.1.6 and §8.7.4), and the CS decomposition (§8.7.6).

# References

The annotated bibliographies at the end of each section remain. Because of space limitations, the master bibliography that was included in previous editions is now available through the book website. References that are historically important have been retained because old ideas have a way of resurrecting themselves. Plus, we must never forget the 1950’s and 1960’s! As mentioned above, we have the luxury of being able to draw upon an expanding library of books on matrix computations. A mnemonic-based citation system has been incorporated that supports these connections to the literature.

# Examples

Non-illuminating, small-n numerical examples have been removed from the text. In their place is a modest suite of Matlab demo scripts that can be run to provide insight into critical theorems and algorithms. We believe that this is a much more effective way to build intuition. The scripts are available through the book website.

# Algorithmic Detail

It is important to have an algorithmic sense and an appreciation for high-performance matrix computations. After all, it is the clever exploitation of advanced architectures that account for much of the field’s soaring success. However, the algorithms that we “formally” present in the book must never be considered as even prototype implementations. Clarity and communication of the big picture are what determine the level of detail in our presentations. Even though specific strategies for specific machines are beyond the scope of the text, we hope that our style promotes an ability to reason about memory traffic overheads and the importance of data locality.

# Acknowledgements

I would like to thank everybody who has passed along typographical errors and suggestions over the years. Special kudos to the Cornell students in CS 4220, CS 6210, and CS 6220, where I used preliminary versions of GVL4. Harry Terkelson earned big bucks through through my ill-conceived \$5-per-typo program!

A number of colleagues and students provided feedback and encouragement during the writing process. Others provided inspiration through their research and books. Thank you all: Diego Accame, David Bindel, ˚Ake Bj¨orck, Laura Bolzano, Jim Demmel, Jack Dongarra, Mark Embree, John Gilbert, David Gleich, Joseph Grcar, Anne Greenbaum, Nick Higham, Ilse Ipsen, Bo K˚agstr¨om, Vel Kahan, Tammy Kolda, Amy Langville, Julian Langou, Lek-Heng Lim, Nicola Mastronardi, Steve McCormick, Mike McCourt, Volker Mehrmann, Cleve Moler, Dianne O’Leary, Michael Overton, Chris Paige, Beresford Parlett, Stefan Ragnarsson, Lothar Reichel, Yousef Saad, Mike Saunders, Rob Schreiber, Danny Sorensen, Pete Stewart, Gil Strang, Francoise Tisseur, Nick Trefethen, Raf Vandebril, and Jianlin Xia.

Chris Paige and Mike Saunders were especially helpful with the editing of Chapters 10 and 11.

Vincent Burke, Jennifer Mallet, and Juliana McCarthy at Johns Hopkins University Press provided excellent support during the production process. Jennifer Slater did a terrific job of copy-editing. Of course, I alone am responsible for all mistakes and oversights.

Finally, this book would have been impossible to produce without my great family and my 4AM writing companion: Henry the Cat!

Charles F. Van Loan

Ithaca, New York

July, 2012

# Global References

A number of books provide broad coverage of the field and are cited multiple times. We identify these global references using mnemonics. Bibliographic details are given in the Other Books section that follows.

AEP Wilkinson: Algebraic Eigenvalue Problem   
ANLA Demmel: Applied Numerical Linear Algebra   
ASNA Higham: Accuracy and Stability of Numerical Algorithms, second edition   
EOM Chatelin: Eigenvalues of Matrices   
FFT Van Loan: Computational Frameworks for the Fast Fourier Transform   
FOM Higham: Functions of Matrices   
FMC Watkins: Fundamentals of Matrix Computations   
IMC Stewart: Introduction to Matrix Computations   
IMK van der Vorst: Iterative Krylov Methods for Large Linear Systems   
IMSL Greenbaum: Iterative Methods for Solving Linear Systems   
ISM Axelsson: Iterative Solution Methods   
IMSLE Saad: Iterative Methods for Sparse Linear Systems, second edition   
LCG Meurant: The Lanczos and Conjugate Gradient Algorithms . . .   
MA Horn and Johnson: Matrix Analysis   
MABD Stewart: Matrix Algorithms: Basic Decompositions   
MAE Stewart: Matrix Algorithms Volume II: Eigensystems   
MEP Watkins: The Matrix Eigenvalue Problem: GR and Krylov Subspace Methods   
MPT Stewart and Sun: Matrix Perturbation Theory   
NLA Trefethen and Bau: Numerical Linear Algebra   
NMA Ipsen: Numerical Matrix Analysis: Linear Systems and Least Squares   
NMLE Saad: Numerical Methods for Large Eigenvalue Problems, revised edition   
NMLS Bj¨orck: Numerical Methods for Least Squares Problems   
NMSE Kressner: Numerical Methods for General and Structured Eigenvalue Problems   
SAP Trefethen and Embree: Spectra and Pseudospectra   
SEP Parlett: The Symmetric Eigenvalue Problem   
SLAS Forsythe and Moler: Computer Solution of Linear Algebraic Systems   
SLS Lawson and Hanson: Solving Least Squares Problems   
TMA Horn and Johnson: Topics in Matrix Analysis

LAPACK LAPACK Users’ Guide, third edition

E. Anderson, Z. Bai, C. Bischof, S. Blackford, J. Demmel, J. Dongarra,

J. Du Croz, A. Greenbaum, S. Hammarling, A. McKenney, and D. Sorensen.

scaLAPACK ScaLAPACK Users’ Guide

L.S. Blackford, J. Choi, A. Cleary, E. D’Azevedo, J. Demmel, I. Dhillon,

J. Dongarra, S. Hammarling, G. Henry, A. Petitet, K. Stanley, D. Walker, and R. C. Whaley.

LIN TEMPLATES Templates for the Solution of Linear Systems . . .

R. Barrett, M.W. Berry, T.F. Chan, J. Demmel, J. Donato, J. Dongarra, V. Eijkhout,

R. Pozo, C. Romine, and H. van der Vorst.

EIG TEMPLATES Templates for the Solution of Algebraic Eigenvalue Problems . . .

Z. Bai, J. Demmel, J. Dongarra, A. Ruhe, and H. van der Vorst.

This page intentionally left blank

# Other Books

The following volumes are a subset of a larger, ever-expanding library of textbooks and monographs that are concerned with matrix computations and supporting areas. The list of references below captures the evolution of the field and its breadth. Works that are more specialized are cited in the annotated bibliographies that appear at the end of each section in the chapters.

# Early Landmarks

V.N. Faddeeva (1959). Computational Methods of Linear Algebra, Dover, New York.   
E. Bodewig (1959). Matrix Calculus, North-Holland, Amsterdam.   
J.H. Wilkinson (1963). Rounding Errors in Algebraic Processes, Prentice-Hall, Englewood Cliffs, NJ.   
A.S. Householder (1964). Theory of Matrices in Numerical Analysis, Blaisdell, New York. Reprinted in 1974 by Dover, New York.   
L. Fox (1964). An Introduction to Numerical Linear Algebra, Oxford University Press, Oxford.   
J.H. Wilkinson (1965). The Algebraic Eigenvalue Problem, Clarendon Press, Oxford.

# General Textbooks on Matrix Computations

G.W. Stewart (1973). Introduction to Matrix Computations, Academic Press, New York.   
R.J. Goult, R.F. Hoskins, J.A. Milner, and M.J. Pratt (1974). Computational Methods in Linear Algebra, John Wiley and Sons, New York.   
W.W. Hager (1988). Applied Numerical Linear Algebra, Prentice-Hall, Englewood Cliffs, NJ.   
P.G. Ciarlet (1989). Introduction to Numerical Linear Algebra and Optimisation, Cambridge University Press, Cambridge.   
P.E. Gill, W. Murray, and M.H. Wright (1991). Numerical Linear Algebra and Optimization, Vol. 1, Addison-Wesley, Reading, MA.   
A. Jennings and J.J. McKeowen (1992). Matrix Computation,second edition, John Wiley and Sons, New York.   
L.N. Trefethen and D. Bau III (1997). Numerical Linear Algebra, SIAM Publications, Philadelphia, PA.   
J.W. Demmel (1997). Applied Numerical Linear Algebra, SIAM Publications, Philadelphia, PA.   
A.J. Laub (2005). Matrix Analysis for Scientists and Engineers, SIAM Publications, Philadelphia, PA.   
B.N. Datta (2010). Numerical Linear Algebra and Applications, second edition, SIAM Publications, Philadelphia, PA.   
D.S. Watkins (2010). Fundamentals of Matrix Computations, John Wiley and Sons, New York.   
A.J. Laub (2012). Computational Matrix Analysis, SIAM Publications, Philadelphia, PA.

# Linear Equation and Least Squares Problems

G.E. Forsythe and C.B. Moler (1967). Computer Solution of Linear Algebraic Systems, Prentice-Hall, Englewood Cliffs, NJ.   
A. George and J.W-H. Liu (1981). Computer Solution of Large Sparse Positive Definite Systems. Prentice-Hall, Englewood Cliffs, NJ.

I.S. Duff, A.M. Erisman, and J.K. Reid (1986). Direct Methods for Sparse Matrices, Oxford University Press, New York.   
R.W. Farebrother (1987). Linear Least Squares Computations, Marcel Dekker, New York.   
C.L. Lawson and R.J. Hanson (1995). Solving Least Squares Problems, SIAM Publications, Philadelphia, PA.   
˚A. Bj¨orck (1996). Numerical Methods for Least Squares Problems, SIAM Publications, Philadelphia, PA.   
G.W. Stewart (1998). Matrix Algorithms: Basic Decompositions, SIAM Publications, Philadelphia, PA.   
N.J. Higham (2002). Accuracy and Stability of Numerical Algorithms, second edition, SIAM Publications, Philadelphia, PA.   
T.A. Davis (2006). Direct Methods for Sparse Linear Systems, SIAM Publications, Philadelphia, PA.   
I.C.F. Ipsen (2009). Numerical Matrix Analysis: Linear Systems and Least Squares, SIAM Publications, Philadelphia, PA.

# Eigenvalue Problems

A.R. Gourlay and G.A. Watson (1973). Computational Methods for Matrix Eigenproblems, John Wiley & Sons, New York.   
F. Chatelin (1993). Eigenvalues of Matrices, John Wiley & Sons, New York.   
B.N. Parlett (1998). The Symmetric Eigenvalue Problem, SIAM Publications, Philadelphia, PA.   
G.W. Stewart (2001). Matrix Algorithms Volume II: Eigensystems, SIAM Publications, Philadelphia, PA.   
L. Komzsik (2003). The Lanczos Method: Evolution and Application, SIAM Publications, Philadelphia, PA.   
D. Kressner (2005). Numerical Methods for General and Structured Eigenvalue Problems, Springer, Berlin.   
D.S. Watkins (2007). The Matrix Eigenvalue Problem: GR and Krylov Subspace Methods, SIAM Publications, Philadelphia, PA.   
Y. Saad (2011). Numerical Methods for Large Eigenvalue Problems, revised edition, SIAM Publications, Philadelphia, PA.

# Iterative Methods

R.S. Varga (1962). Matrix Iterative Analysis, Prentice-Hall, Englewood Cliffs, NJ.   
D.M. Young (1971). Iterative Solution of Large Linear Systems, Academic Press, New York.   
L.A. Hageman and D.M. Young (1981). Applied Iterative Methods, Academic Press, New York.   
J. Cullum and R.A. Willoughby (1985). Lanczos Algorithms for Large Symmetric Eigenvalue Computations, Vol. I Theory, Birkha¨user, Boston.   
J. Cullum and R.A. Willoughby (1985). Lanczos Algorithms for Large Symmetric Eigenvalue Computations, Vol. II Programs, Birkha¨user, Boston.   
W. Hackbusch (1994). Iterative Solution of Large Sparse Systems of Equations, Springer-Verlag, New York.   
O. Axelsson (1994). Iterative Solution Methods, Cambridge University Press.   
A. Greenbaum (1997). Iterative Methods for Solving Linear Systems, SIAM Publications, Philadelphia, PA.   
Y. Saad (2003). Iterative Methods for Sparse Linear Systems, second edition, SIAM Publications, Philadelphia, PA.   
H. van der Vorst (2003). Iterative Krylov Methods for Large Linear Systems, Cambridge University Press, Cambridge, UK.

G. Meurant (2006). The Lanczos and Conjugate Gradient Algorithms: From Theory to Finite Precision Computations, SIAM Publications, Philadelphia, PA.

# Special Topics/Threads

L.N. Trefethen and M. Embree (2005). Spectra and Pseudospectra—The Behavior of Nonnormal Matrices and Operators, Princeton University Press, Princeton and Oxford.   
R. Vandebril, M. Van Barel, and N. Mastronardi (2007). Matrix Computations and Semiseparable Matrices I: Linear Systems, Johns Hopkins University Press, Baltimore, MD.   
R. Vandebril, M. Van Barel, and N. Mastronardi (2008). Matrix Computations and Semiseparable Matrices II: Eigenvalue and Singular Value Methods, Johns Hopkins University Press, Baltimore, MD.   
N.J. Higham (2008) Functions of Matrices, SIAM Publications, Philadelphia, PA.

# Collected Works

R.H. Chan, C. Greif, and D.P. O’Leary, eds. (2007). Milestones in Matrix Computation: Selected Works of G.H. Golub, with Commentaries, Oxford University Press, Oxford. M.E. Kilmer and D.P. O’Leary, eds. (2010). Selected Works of G.W. Stewart, Birkhauser, Boston, MA.

# Implementation

B.T. Smith, J.M. Boyle, Y. Ikebe, V.C. Klema, and C.B. Moler (1970). Matrix Eigensystem Routines: EISPACK Guide, second edition, Lecture Notes in Computer Science, Vol. 6, Springer-Verlag, New York.   
J.H. Wilkinson and C. Reinsch, eds. (1971). Handbook for Automatic Computation, Vol. 2, Linear Algebra, Springer-Verlag, New York.   
B.S. Garbow, J.M. Boyle, J.J. Dongarra, and C.B. Moler (1972). Matrix Eigensystem Routines: EISPACK Guide Extension, Lecture Notes in Computer Science, Vol. 51, Springer-Verlag, New York.   
J.J Dongarra, J.R. Bunch, C.B. Moler, and G.W. Stewart (1979). LINPACK Users’ Guide, SIAM Publications, Philadelphia, PA.   
K. Gallivan, M. Heath, E. Ng, B. Peyton, R. Plemmons, J. Ortega, C. Romine, A. Sameh, and R. Voigt (1990). Parallel Algorithms for Matrix Computations, SIAM Publications, Philadelphia, PA.   
R. Barrett, M.W. Berry, T.F. Chan, J. Demmel, J. Donato, J. Dongarra, V. Eijkhout, R. Pozo, C. Romine, and H. van der Vorst (1993). Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods, SIAM Publications, Philadelphia, PA.   
L.S. Blackford, J. Choi, A. Cleary, E. D’Azevedo, J. Demmel, I. Dhillon, J. Dongarra, S. Hammarling, G. Henry, A. Petitet, K. Stanley, D. Walker, and R.C. Whaley (1997). ScaLA-PACK Users’ Guide, SIAM Publications, Philadelphia, PA.   
J.J. Dongarra, I.S. Duff, D.C. Sorensen, and H.A. van der Vorst (1998). Numerical Linear Algebra on High-Performance Computers, SIAM Publications, Philadelphia, PA.   
E. Anderson, Z. Bai, C. Bischof, S. Blackford, J. Demmel, J. Dongarra, J. Du Croz, A. Greenbaum, S. Hammarling, A. McKenney, and D. Sorensen (1999). LAPACK Users’ Guide, third edition, SIAM Publications, Philadelphia, PA.   
Z. Bai, J. Demmel, J. Dongarra, A. Ruhe, and H. van der Vorst (2000). Templates for the Solution of Algebraic Eigenvalue Problems: A Practical Guide, SIAM Publications, Philadelphia, PA.   
V.A. Barker, L.S. Blackford, J. Dongarra, J. Du Croz, S. Hammarling, M. Marinova, J. Wasniewski, and P. Yalamov (2001). LAPACK95 Users’ Guide, SIAM Publications, Philadelphia.

# Matlab

D.J. Higham and N.J. Higham (2005). MATLAB Guide, second edition, SIAM Publications, Philadelphia, PA.

R. Pratap (2006). Getting Started with Matlab 7, Oxford University Press, New York.

C.F. Van Loan and D. Fan (2009). Insight Through Computing: A Matlab Introduction to Computational Science and Engineering, SIAM Publications, Philadelphia, PA.

# Matrix Algebra and Analysis

R. Horn and C. Johnson (1985). Matrix Analysis, Cambridge University Press, New York.

G.W. Stewart and J. Sun (1990). Matrix Perturbation Theory, Academic Press, San Diego.

R. Horn and C. Johnson (1991). Topics in Matrix Analysis, Cambridge University Press, New York.

D.S. Bernstein (2005). Matrix Mathematics, Theory, Facts, and Formulas with Application to Linear Systems Theory, Princeton University Press, Princeton, NJ.

L. Hogben (2006). Handbook of Linear Algebra, Chapman and Hall, Boca Raton, FL.

# Scientific Computing/Numerical Analysis

G.W. Stewart (1996). Afternotes on Numerical Analysis, SIAM Publications, Philadelphia, PA.

C.F. Van Loan (1997). Introduction to Scientific Computing: A Matrix-Vector Approach Using Matlab, Prentice Hall, Upper Saddle River, NJ.

G.W. Stewart (1998). Afternotes on Numerical Analysis: Afternotes Goes to Graduate School, SIAM Publications, Philadelphia, PA.

M.T. Heath (2002). Scientific Computing: An Introductory Survey, second edition), McGraw-Hill, New York.

C.B. Moler (2008) Numerical Computing with MATLAB, revised reprint, SIAM Publications, Philadelphia, PA.

G. Dahlquist and ˚A. Bj¨orck (2008). Numerical Methods in Scientific Computing, Vol. 1, SIAM Publications, Philadelphia, PA.

U. Ascher and C. Greif (2011). A First Course in Numerical Methods, SIAM Publications, Philadelphia, PA.

# Useful URLs

# GVL4

Matlab demo scripts and functions, master bibliography, list of errata.

http://www.cornell.edu/cv/GVL4

# Netlib

Huge repository of numerical software including LAPACK.

http://www.netlib.org/index.html

# Matrix Market

Test examples for matrix algorithms.

http://math.nist.gov/MatrixMarket/

# Matlab Central

Matlab functions, demos, classes, toolboxes, videos.

http://www.mathworks.com/matlabcentral/

# University of Florida Sparse Matrix Collections

Thousands of sparse matrix examples in several formats.

http://www.cise.ufl.edu/research/sparse/matrices/

# Pseudospectra Gateway

Grapical tools for pseudospectra.

http://www.cs.ox.ac.uk/projects/pseudospectra/

# ARPACK

Software for large sparse eigenvalue problems

http://www.caam.rice.edu/software/ARPACK/

# Innovative Computing Laboratory

State-of-the-art high performance matrix computations.

http://icl.cs.utk.edu/

This page intentionally left blank

# Common Notation

<table><tr><td> $\mathbb{R}, \mathbb{R}^{n}, \mathbb{R}^{m \times n}$ </td><td>set of real numbers, vectors, and matrices (p. 2)</td></tr><tr><td> $\mathbb{C}, \mathbb{C}^{n}, \mathbb{C}^{m \times n}$ </td><td>set of complex numbers, vectors, and matrices (p. 13)</td></tr><tr><td> $a_{ij}, A(i,j), [A]_{ij}$ </td><td> $(i,j)$  entry of a matrix (p. 2)</td></tr><tr><td> $\mathbf{u}$ </td><td>unit roundoff (p. 96)</td></tr><tr><td> $\text{fl}(\cdot)$ </td><td>floating point operator (p. 96)</td></tr><tr><td> $\| x \|_{p}$ </td><td> $p$ -norm of a vector (p. 68)</td></tr><tr><td> $\| A \|_{p}, \| A \|_{F}$ </td><td> $p$ -norm and Frobenius norm of a matrix (p. 71)</td></tr><tr><td> $\text{length}(x)$ </td><td>dimension of a vector (p. 236)</td></tr><tr><td> $\kappa_{p}(A)$ </td><td> $p$ -norm condition (p. 87)</td></tr><tr><td> $|A|$ </td><td>absolute value of a matrix (p. 91)</td></tr><tr><td> $A^{T}, A^{H}$ </td><td>transpose and conjugate transpose (p. 2, 13)</td></tr><tr><td> $\text{house}(x)$ </td><td>Householder vector (p. 236)</td></tr><tr><td> $\text{givens}(a,b)$ </td><td>cosine-sine pair (p. 240)</td></tr><tr><td> $x_{LS}$ </td><td>minimum-norm least squares solution (p. 260)</td></tr><tr><td> $\text{ran}(A)$ </td><td>range of a matrix (p. 64)</td></tr><tr><td> $\text{null}(A)$ </td><td>nullspace of a matrix (p. 64)</td></tr><tr><td> $\text{span}\{v_{1},\ldots,v_{n}\}$ </td><td>span defined by vectors (p. 64)</td></tr><tr><td> $\text{dim}(S)$ </td><td>dimension of a subspace (p. 64)</td></tr><tr><td> $\text{rank}(A)$ </td><td>rank of a matrix (p. 65)</td></tr><tr><td> $\text{det}(A)$ </td><td>determinant of a matrix (p. 66)</td></tr><tr><td> $\text{tr}(A)$ </td><td>trace of a matrix (p. 327)</td></tr><tr><td> $\text{vec}(A)$ </td><td>vectorization of a matrix (p. 28)</td></tr><tr><td> $\text{reshape}(A,p,q)$ </td><td>reshaping a matrix (p. 28)</td></tr><tr><td> $\text{Re}(A), \text{Im}(A)$ </td><td>real and imaginary parts of a matrix (p. 13)</td></tr><tr><td> $\text{diag}(d_{1},\ldots,d_{n})$ </td><td>diagonal matrix (p. 18)</td></tr><tr><td> $I_{n}$ </td><td> $n$ -by- $n$  identity matrix (p. 19)</td></tr><tr><td> $e_{i}$ </td><td> $i$ th column of the identity matrix (p. 19)</td></tr><tr><td> $\mathcal{E}_{n}, \mathcal{D}_{n}, \mathcal{P}_{p,q}$ </td><td>exchange, downshift, and perfect shuffle permutations (p. 20)</td></tr><tr><td> $\sigma_{i}(A)$ </td><td> $i$ th largest singular value (p. 77)</td></tr><tr><td> $\sigma_{\text{max}}(A), \sigma_{\text{min}}(A)$ </td><td>largest and smallest singular value (p. 77)</td></tr><tr><td> $\text{dist}(S_{1},S_{2})$ </td><td>distance between two subspaces (p. 82)</td></tr><tr><td> $\text{sep}(A_{1},A_{2})$ </td><td>separation between two matrices (p. 360)</td></tr><tr><td> $\lambda(A)$ </td><td>set of eigenvalues (p. 66)</td></tr><tr><td> $\lambda_{i}(A)$ </td><td> $i$ th largest eigenvalue of a symmetric matrix (p. 66)</td></tr><tr><td> $\lambda_{\text{max}}(A), \lambda_{\text{min}}(A)$ </td><td>largest and smallest eigenvalue of a symmetric matrix (p. 66)</td></tr><tr><td> $\rho(A)$ </td><td>spectral radius (p. 349)</td></tr><tr><td> $\mathcal{K}(A,q,j)$ </td><td>Krylov subspace (p. 548)</td></tr></table>

This page intentionally left blank

# Matrix Computations

This page intentionally left blank
