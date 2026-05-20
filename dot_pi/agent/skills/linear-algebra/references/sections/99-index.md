# Index

A-conjugate, 633

A-norm, 629

Aasen’s method, 188–90

Absolute value notation, 91

Additive Schwarz, 665

Adjacency set, 602

Affine space, 628

Algebraic multiplicity, 353

Algorithm, 135

Algorithmic detail, xii

Angles between subspaces, 329–31

Antidiagonal, 208

Approximate inverse preconditioner, 654

Approximate Newton method, 590–1

Approximation of a matrix function, 522–3

Arnoldi process, 579–83

implicit restarting, 581–3

k-step decomposition, 580

rational, 588

Arnoldi vectors, 580

Augmented system method, 316

Back substitution, 107

Backward error analysis, 100–1

Backward stable, 136

Backward successive over-relaxation, 619

Balancing, 392

Band algorithms, 176ff

Cholesky, 180

Gaussian elimination, 178–9

Hessenberg LU, 179

triangular systems, 177–8

Band matrix, 15

data structures and, 17

inverse, 182–3

LU factorization and, 176–7

pivoting and, 178–9

profile Cholesky, 184

Bandwidth, 176

barrier, 57

Bartels-Stewart algorithm, 398–400

Basic solution for least squares, 292

Basis, 64

eigenvector, 400

Bauer-Fike theorem, 357–8

BiCGstab method, 647

Biconjugate gradient method, 645

Bidiagonalization

Golub-Kahan, 572

Householder, 284–5

Paige-Saunders, 575

upper triangularizing first, 285

Bidiagonal matrix, 15

Big-Oh notation, 12

Binary powering, 527

Bisection methods

for Toeplitz eigenproblem, 216

for tridiagonal eigenproblem, 467

BLAS, 12ff

Block algorithms, 196ff

Cholesky, 168–70

cyclic reduction, 197–8

data reuse and, 47ff

diagonalization, 352–3, 397–400

Gaussian elimination, 144–6

Gauss-Seidel, 613

Jacobi method for eigenvalues, 481–2

Jacobi method for linear systems, 613

Lanczos, 566–9

LU, 118–20, 196–7

LU with pivoting, 144–6

multiple right-hand-side triangular, 109–10

QR factorization, 250–1

recursive QR factorization, 251

SPIKE, 197–8

Tridiagonal, 196ff

unsymmetric Lanczos, 586

Block-cyclic distribution layout, 50

and parallel LU, 146

Block diagonal dominance, 197

Block distribution layout, 50

Block Householder, 238–9

Block matrices, 22ff

data reuse and, 55

diagonal dominance of, 197

Block tridiagonal systems, 196ff

Bordered linear systems, 202

Bunch-Kaufman algorithm, 192

Bunch-Parlett pivoting, 191

Cache, 46

Cancellation, 97

Cannon’s algorithm, 60–1

Canonical correlation problem, 330

Cauchy-like matrix, 682

conversion to, 689–90

Cauchy-Schwarz inequality, 69

Cayley transform, 68, 245

CGNE, 636–7

CGNR, 636

Characteristic polynomial, 66, 348

generalized eigenproblem and, 405

Chebyshev polynomials, 621, 653

Chebyshev semi-iterative method, 621–2

Cholesky factorization, 163

band, 180

block, 168–70

downdating and, 338–41

gaxpy version, 164

matrix square root and, 163

profile, 184

recursive block, 172–3

stability of, 164–5

Cholesky reduction of A − λB, 500

Chordal metric, 407

Circulant systems, 220–2

Classical Gram-Schmidt, 254

Coarse grid role in multigrid, 673

Collatz-Wielandt formula, 373

Colon notation, 6, 16

Column

deletion or addition in QR, 235–8

major order, 45

ordering in QR factorization, 279–80

orientation, 5, 107–8

partitioning, 6

pivoting, 276–7

weighting in least squares, 306–7

Communication costs, 52ff

Compact WY transformation, 244

Companion matrix, 382–3

Complete orthogonal decomposition, 283

Complete pivoting, 131–3

Complex

Givens transformation, 243–4

Householder transformation, 243

matrices, 13

matrix multiplication, 29

QR factorization, 256

Complexity of matrix inversion, 174

Componentwise bounds, 92

Compressed column representation, 598–9

Computation/communication ratio, 53ff

Condition estimation, 140, 142–3, 436

Condition of

eigenvalues, 359–60

invariant subspaces, 360–1

least squares problem, 265–7

linear systems, 87–8

multiple eigenvalues, 360

rectangular matrix, 248

similarity transformation, 354

Confluent Vandermonde matrix, 206

Conformal partition, 23

Congruence transformation, 449

Conjugate

directions, 633

transpose, 13

Conjugate gradient method, 625ff

derivation and properties, 629–30, 633

Hestenes-Stiefel version, 634–5

Lanczos version, 632

practical, 635–6

pre-conditioned, 651–2

Conjugate gradient squared method, 646

Consistent norms, 71

Constrained least squares, 313–4

Contour integral and f(A), 528–9

Convergence. See under particular algorithm

Courant-Fischer minimax theorem, 441

CP approximation, 735–8

Craig’s method, 637

Crawford number, 499

Cross product, 70

Cross-validation, 308

CS decomposition, 84–5, 503–6

hyperbolic, 344

subset selection and, 294

thin version, 84

CUR decomposition, 576

Curse of dimensionality, 741

Cuthill-McKee ordering, 602–4

Cyclic Jacobi method, 480–1

Cyclic reduction, 197–8

Data least squares, 325

Data motion overhead, 53

Data reuse, 46–8

Data sparse, 154

Davidson method, 593–4

Decompositions and factorizations

Arnoldi, 580

bidiagonal, 5

block diagonal, 397–9

Cholesky, 163

companion matrix, 382

complete orthogonal, 283

CS (general), 85

CS (thin), 84

generalized real Schur, 407

generalized Schur, 406–7

Hessenberg, 378ff

Hessenberg-triangular, 408–9

Jordan, 354

LU, 114, 128

QR, 247

QR (thin version), 248

real Schur, 376

Schur, 351

singular value, 76

singular value (thin), 80

symmetric Schur, 440

tridiagonal, 458–9

Decoupling in eigenproblem, 349–50

Defective eigenvalue, 66, 353

Deflating subspace, 404

Deflation and

bidiagonal form, 490

Hessenberg-triangular form, 409–10

QR algorithm, 385

Denman-Beavers iteration, 539–40

Departure from normality, 351

Derogatory matrix, 383

Determinant, 66, 348

Gaussian elimination and, 114

and singularity, 89

Vandermonde matrix, 206

Diagonal dominance, 154–6, 615

block, 197

LU and, 156

Diagonal matrix, 18

Diagonal pivoting method, 191–2

Diagonal plus rank-1, 469–71

Diagonalizable, 67, 353

Differentiation of matrices, 67

Dimension, 64

Direct methods, 598ff

Dirichlet end condition, 222

Discrete cosine transform (DCT), 39

Discrete Fourier transform (DFT), 33-6

circulant matrices and, 221-2

factorizations and, 41

matrix, 34

Discrete Poisson problem

1-dimensional, 222-4

2-dimensional, 224-31

Discrete sine transform (DST), 39

Displacement rank, 682

Distance between subspaces, 82

Distributed memory model, 57

Divide-and-conquer algorithms

cyclic reduction, 197–8

Strassen, 30–1

tridiagonal eigenvalue, 471–3

Domain decomposition, 662–5

Dominant

eigenvalue, 366

invariant subspace, 368

Dot product, 4, 10

Dot product roundoff, 98

Double implicit shift, 388

Doubling formulae, 526

Downdating Cholesky, 338–41

Drazin inverse, 356

Durbin’s algorithm, 210

Eckhart-Young theorem, 79

Eigenproblem

diagonal plus rank-1, 469–71

generalized, 405ff, 497ff

inverse, 473–4

orthogonal Hessenberg matrix, 703–4

symmetric, 439ff

Toeplitz, 214–6

unsymmetric, 347ff

Eigensystem

fast, 219

Eigenvalue decompositions

Jordan, 354

Schur, 351

Eigenvalues

algebraic multiplicity, 353

characteristic polynomial and, 348

computing selected, 453

defective, 66

determinant and, 348

dominant, 366

generalized, 405

geometric multiplicity, 353

ordering in Schur form, 351, 396–7

orthogonal Hessenberg, 703–4

relative perturbation, 365

repeated, 360

sensitivity (symmetric case), 441–3

sensitivity (unsymmetric case), 359–60

singular values and, 355

Sturm sequence and, 468

symmetric tridiagonal, 467ff

trace, 348

unstable, 363

Eigenvector, 67

basis, 400

dominant, 366

left, 349

matrix and condition, 354

perturbation, 361–2

right, 349

Elementary Hermitian matrices.

See Householder matrix

Elementary transformations. See

Gauss transformations

Equality constained least squares, 315–7

Equilibration, 139

Equilibrium systems, 192–3

Equivalence of vector norms, 69

Error

absolute, 69

damping in multigrid, 622-3

relative, 70

roundoff, 96–102

Error analysis

backward, 100

forward, 100

Euclidean matrix norm. See

Frobenius matrix norm

Exchange permutation matrix, 20

Explicit shift in QR algorithm symmetric case, 461

unsymmetric case, 385–8

Exponential of matrix, 530–6

Factored form representation, 237–8

Factorization. See Decompositions and factorizations

Fast methods

cosine transform, 36ff

eigensystem, 219, 228–31

Fourier transform, 33ff

Givens QR, 245

Poisson solver, 226–7

sine transform, 36

Field of values, 349

Fine grid role in multigrid, 673

Floating point

fundamental axiom, 96

maxims, 96–7

normalized, 94

numbers, 93

storage of matrix, 97–8

Flopcounts, 12, 16

for square system methods, 298

F -norm, 71

Forward error analysis, 100

Forward substitution, 106

Francis QR step, 390

Frechet derivative, 521

Frobenius matrix norm, 71

Frontal methods, 610

Full multigrid, 678

Function of matrix, 513ff

eigenvectors and, 517–8

Schur decomposition and, 518–20

Taylor series and, 524–6

Gauss-Jordan transformations, 121

Gauss-Radau rule, 560–1

Gauss rules, 557–9

Gauss-Seidel iteration, 611-2

block, 613

Poisson equation and, 617

positive definite systems and, 615

Gauss transformations, 112-3

Gaussian elimination, 111ff

banded version, 176–9

block version, 144–5

complete pivoting and, 131–2

gaxpy version, 117

outer product version, 116

partial pivoting and, 127

rook pivoting and, 133

roundoff error and, 122–3

tournament pivoting and, 150

Gaxpy, 5

blocked, 25

Gaxpy-rich algorithms

Cholesky, 164

Gaussian elimination, 129–30

LDLT , 157–8

Gaxpy vs. outer product, 45

Generalized eigenproblem, 405ff

Generalized eigenvalues, 405

sensitivity, 407

Generalized least squares, 305–6

Generalized Schur decomposition, 406–7 computation of, 502–3

Generalized singular vectors, 502

Generalized SVD, 309–10, 501–2 constrained least squares and, 316–7

Generalized Sylvester equation, 417

Generator representation, 693

Geometric multiplicity, 353

Gershgorin theorem, 357, 442

Ghost eigenvalues, 566

givens, 240

Givens QR, 252–3 parallel, 257

Givens rotations, 239–42 complex, 243–4

fast, 245

rank-revealing decompositions and, 280–2 square-root free, 246

Global memory, 55

GMRES, 642–4 m-step, 644 preconditioned, 652–3

Golub-Kahan bidiagonalization, 571–3 SVD step, 491

Gram-Schmidt classical, 254 modified, 254–5

Graph, 602

Graphs and sparsity, 601–2

Growth in Gaussian elimination, 130–2

Haar wavelet transform, 40ff factorization, 41

Hadamard product, 710

Hamiltonian matrix, 29, 420 eigenvalue problem, 420–1

Hankel-like, 688–9

Hermitian matrix, 18

Hessenberg form, 15 Arnoldi process and, 579–80 Householder reduction to, 378–9 inverse iteration and, 395 properties, 381–2 QR factorization and, 253–4 QR iteration and, 385–6 unreduced, 381

Hessenberg QR step, 377-8

Hessenberg systems, 179 LU and, 179

Hessenberg-triangular form, 408–9

Hierarchical memory, 46

Hierarchical rank structure, 702

Higher-order SVD, 732–3 truncated, 734

Holder inequality, 69

Horner algorithm, 526–7

house, 236

Householder bidiagonalization, 284–5 tridiagonalization, 458–9

Householder matrix, 234–8 complex, 243 operations with, 235–7

Hyperbolic CS decomposition, 344 rotations, 339 transformations, 339

Identity matrix, 19

Ill-conditioned matrix, 88

IEEE arithmetic, 94

Im, 13

Implicit Q theorem symmetric matrix version, 460 unsymmetric matrix version, 381

Implicit symmetric QR step with Wilkinson Shift, 461–2

Implicitly restarted Arnoldi method, 581–3

Incomplete block preconditioners, 657–60

Incomplete Cholesky, 357–60

Indefinite least squares, 344

Indefinite symmetric matrix, 159

Indefinite systems, 639–41

Independence, 64

Inertia of symmetric matrix, 448

inf , 95

Integrating f(A), 527–8

Interchange permutation, 126

Interlacing property singular values, 487 symmetric eigenvalues, 443

Intersection nullspaces, 328–9 subspaces, 331

Invariant subspace approximate, 446–8 dominant, 378 perturbation of (sy perturbation of (un Schur vectors and, 3

Inverse, 19 band matrices and, 182–3

Inverse eigenvalue problems, 473–4

Inverse error analysis. See Backward error analysis

Inverse fast transforms cosine, 227–8 Fourier, 220 sine, 227–8

Inverse iteration generalized eigenproblem, 414 symmetric case, 453 unsymmetric case, 394–5

Inverse low-rank perturbation, 65

Inverse of matrix, perturbation of, 74 Toeplitz case, 212–3

Inverse orthogonal iteration, 374

Inverse power method, 374

Inverse scaling and squaring, 542

Irreducible, 373

Iteration matrix, 613

Iterative improvement fixed precision and, 140 least squares, 268–9, 272 linear systems, 139–40

Iterative methods, 611–50

Jacobi iteration for the SVD, 492–3

Jacobi iteration for symmetric eigenproblem, 476ff classical, 479–80 cyclic, 480 error, 480–1 parallel version, 482–3

Jacobi method for linear systems, block version, 613

diagonal dominance and, 615

preconditioning with, 653

Jacobi orthogonal correction method, 591–3

Jacobi rotations, 477

Jacobi-Davidson method, 594–5

Jordan blocks, 400-2

Jordan decomposition, 354

computation of, 400-2

matrix functions and, 514, 522-3

Kaniel-Paige-Saad theory, 552–4

Khatri-Rao product, 710

Kogbetiantz algorithm, 506

Kronecker product, 27

basic properties, 27, 707–8

multiple, 28, 716

SVD 712–4

Kronecker structure, 418

Krylov

matrix, 459

subspaces, 548

Krylov-Schur algorithm, 584

Krylov subspace methods

biconjugate gradients, 645

CG (conjugate gradients), 625ff

CGNE (conjugate gradient normal equation error), 637–8

CGNR (conjugate gradient normal equation residual), 637–8

CGS (conjugate gradient squared), 646

general linear systems and, 579ff

GMRES (general minimum residual), 642–5

MINRES (minimum residual), 639–40

QMR (quasi-minimum residual), 647

SYMMLQ, 640–1

Krylov subspace methods for

general linear systems, 636–7, 642–7

least squares, 641–2

singular values, 571–8

symmetric eigenproblem, 546–56, 562–71

symmetric indefinite systems, 639–41

symmetric positive definite systems, 625–39

unsymmetric eigenproblem, 579–89

Lagrange multipliers, 313

Lanczos tridiagonalization, 546ff

block version, 566-9

complete reorthogonalization and, 564–5

conjugate gradients and, 628–32

convergence of, 552–4

Gauss quadrature and, 560–1

interior eigenvalues and, 553-4

orthogonality loss, 564

power method and, 554–5

practical, 562ff

Ritz approximation and, 551–2

roundoff and, 563–4

selective orthogonalization and, 565–6

s-step, 569

termination of, 549

unsymmetric, 584–7

Lanczos vectors, 549

LDLT , 156–8

conjugate gradients and, 631

with pivoting, 165–6

Leading principal submatrix, 24

Least squares methods, flopcounts for, 293

Least squares problem

basic solution to, 292

cross-validation and, 308

equality constraints and, 315–7

full rank, 260ff

generalized, 305–6

indefinite, 344

iterative improvement, 268–9

Khatri-Rao product and, 737

minimum norm solution to, 288–9

quadratic inequality constraint, 313–5

rank deficient, 288ff

residual vs. column independence, 295–6

sensitivity of, 265–7

solution set of, 288

solution via Householder QR, 263–4

sparse, 607–8, 641–2

SVD and, 289

Least squares solution using

LSQR, 641–2

modified Gram-Schmidt, 264–5

normal equations, 262–3

QR factorization, 263–4

seminormal equations, 607

SVD, 289

Left eigenvector, 349

Left-looking, 117

Levels of linear algebra, 12

Level-3 fraction, 109

block Cholesky, 170

block LU, 120

Hessenberg reduction, 380

Levinson algorithm, 211

Linear equation sensitivity, 102, 137ff

Linear independence, 64

Linearization, 415–6

Load balancing, 50ff

Local memory, 50

Local program, 50

Log of a matrix, 541–2

Look-ahead, 217, 586–7

Loop reordering, 9

Loss of orthogonality

Gram-Schmidt, 254

Lanczos, 564

Low-rank approximation

randomized, 576–7

SVD, 79

LR iteration, 370

LSMR, 642

LSQR, 641–2

LU factorization, 111ff

band, 177

block, 196–7

Cauchy-like, 685–6

determinant and, 114

diagonal dominance and, 155

differentiation of, 120

existence of, 114

gaxpy version, 117

growth factor and, 130–1

Hessenberg, 179

mentality, 134

outer product version, 116

partial pivoting and, 128

rectangular matrices and, 118

roundoff and, 122-3

semiseparable, 695–7

sparse, 608–9

Machine precision, 95

Markov chain, 374

Markowitz pivoting, 609

Matlab, xix

Matrix functions, 513ff

integrating, 527–8

Jordan decomposition and, 514–5

polynomial evaluation and, 526–7

sensitivity of, 520–1

Matrix multiplication, 2, 8ff

blocked, 26

Cannon’s algorithm, 60–1

distributed memory, 50ff

dot product version, 10

memory hierarchy and, 47

outer product version, 11

parallel, 49ff

saxpy version, 11

Strassen, 30–1

tensor contractions and, 726–7

Matrix norms, 71–3

consistency, 71

Frobenius, 71

relations between, 72–3

subordinate, 72

Matrix-vector products, 33ff

blocked, 25

Memory hierarchy, 46

Minimax theorem for

singular values, 487

symmetric eigenvalues, 441

Minimum degree ordering, 604–5

Minimum singular value, 78

MINRES, 639–41

Mixed packed format, 171

Mixed precision, 140

Modal product, 727–8

Modal unfoldings, 723

Modified Gram-Schmidt, 254–5

and least squares, 264–5

Modified LR algorithm, 392

Moore-Penrose conditions, 290

Multigrid, 670ff

Multilinear product, 728–9

Multiple eigenvalues,

matrix functions and, 520

unreduced Hessenberg matrices and, 382

Multiple-right-hand-side problem, 108

Multiplicative Schwarz, 664

Multipliers in Gauss transformations, 112

NaN, 95

Nearness to

Kronecker product, 714–5

singularity, 88

skew-hermitian, 449

Nested-dissection ordering, 605–6

Netlib, xix

Neumann end condition, 222

Newton method for Toeplitz eigenvalue, 215

Newton-Schultz iteration, 538

nnz, 599

Node degree, 602

Nonderogatory matrices, 383

Nongeneric total least squares, 324

Nonsingular matrix, 65

Norm

matrix, 71–3

vector, 68

Normal equations, 262–3, 268

Normal matrix, 351

departure from, 351

Normwise-near preconditioners, 654

null, 64

Nullity theorem, 185

Nullspace, 64

intersection of, 328–9

Numerical radius, 349

Numerical range, 349

Numerical rank

least squares and, 291

QR with column pivoting and, 278–9

SVD and, 275–6

off, 477

Ordering eigenvalues, 396–7

Ordering for sparse matrices

Cuthill-McKee, 602–4

minimum degree, 604–6

nested dissection, 605–7

Orthogonal

complement, 65

invariance, 75

matrix, 66, 234

Procrustes problem, 327–8

projection, 82

symplectic matrix, 420

vectors, 65

Orthogonal iteration

symmetric, 454–5, 464–5

unsymmetric, 367–8, 370–3

Orthogonal matrix representations

factored form, 237–8

Givens rotations, 242

WY block form, 238–9

Orthogonality between subspaces, 65

Orthonormal basis computation, 247

Outer product, 7

Gaussian elimination and, 115

LDLT and, 166

sparse, 599–600

between tensors, 724

versus gaxpy, 45

Overdetermined system, 260

Packed format, 171

Pad´e approximation, 530–1

PageRank, 374

Parallel computation

divide and conquer eigensolver, 472–3

Givens QR, 257

Jacobi’s eigenvalue method, 482–3

LU, 144ff

matrix multiplication, 49ff

Parlett-Reid method, 187–8

Parlett-Schur method, 519

block version, 520

Partitioning

conformable, 23

matrix, 5–6

Pencils, equivalence of, 406

Perfect shuffle permutation, 20, 460, 711–2

Periodic end conditions, 222

Permutation matrices, 19ff

Perron-Frobenius theorem, 373

Perron’s theorem, 373

Persymmetric matrix, 208

Perturbation results

eigenvalues (symmetric case), 441–3

eigenvalues (unsymmetric case), 357–60

eigenvectors (symmetric case), 445–6

eigenvectors (unsymmetric case), 361–2

generalized eigenvalue, 407

invariant subspaces (symmetric case), 444–5

invariant subspaces (unsymmetric case), 361

least squares problem, 265–7

linear equation problem, 82–92

singular subspace pair, 488

singular values, 487

underdetermined systems, 301

Pipelining, 43

Pivoting

Aasen’s method and, 190

Bunch-Kaufman, 192

Bunch-Parlett, 191

Cauchy-like and, 686–7

column, 276–8

complete, 131–2

LU and, 125ff

Markowitz, 609

partial, 127

QR and, 279–80

rook, 133

symmetric matrices and, 165–6

tournament, 150

Plane rotations. See Givens rotations

p-norms, 71

minimization in, 260

Point, line, plane problems, 269–271

Pointwise operations, 3

Polar decomposition, 328, 540–1

Polynomial approximation and GMRES, 644

Polynomial eigenvalue problem, 414–7

Polynomial interpolation, Vandermonde systems and, 203–4

Polynomial preconditioner, 655-6

Positive definite systems, 159ff

Gauss-Seidel and, 615–6

LDLT and, 165ff

properties of, 159–61

unsymmetric, 161–3

Positive matrix, 373

Positive semidefinite matrix, 159

Post-smoothing in multigrid, 675

Power method, 365ff

error estimation in, 367

symmetric case, 451–2

Power series of matrix, 524

Powers of a matrix, 527

Preconditioned

conjugate gradient method, 651–2, 656ff

GMRES, 652–3

Preconditioners, 598

approximate inverse, 654-5

domain decomposition, 662–5

incomplete block, 660-1

incomplete Cholesky, 657–60

Jacobi and SSOR, 653

normwise-near, 654

polynomial, 655

saddle point, 661

Pre-smoothing role in multigrid, 675

Principal angles and vectors, 329–31

Principal square root, 539

Principal submatrix, 24

Probability vector, 373

Procrustes problem, 327–8

Product eigenvalue problem, 423–5

Product SVD problem, 507

Profile, 602

Cholesky, 184

indices, 184, 602

Projections, 82

Prolongation matrix, 673

Pseudo-eigenvalue, 428

Pseudoinverse, 290, 296

Pseudospectra, 426ff

computing plots, 433–4

matrix exponential and, 533–4

properties, 431–3

Pseudospectral abscissa, 434–5

Pseudospectral radius, 434–5

QMR, 647

QR algorithm for eigenvalues

Hessenberg form and 377–8

shifts and, 385ff

symmetric version, 456ff

tridiagonal form and, 460

unsymmetric version, 391ff

Wilkinson shift, 462–3

QR factorization, 246ff

block Householder, 250–1

block recursive, 251

classical Gram-Schmidt and, 254

column pivoting and, 276–8

complex, 256

Givens computation of, 252–3

Hessenberg matrices and, 253–4

Householder computation of, 248–9

least square problem and, 263–4

modified Gram-Schmidt and, 254–5

properties of, 246–7

range space and, 247

rank of matrix and, 274

sparse, 606–8

square systems and, 298–9

thin version, 248

tridiagonal matrix and, 460

underdetermined systems and, 300

updating, 335–8

Quadratic eigenvalue problem, 507–8

Quadratically constrained least squares, 314–5

Quasidefinite matrix, 194

Quasiseparable matrix, 693

Quotient SVD, 507

QZ algorithm, 412–3

step, 411-2

ran, 64

Randomization, 576–7

Range of a matrix,

orthonormal basis for, 247

Rank of matrix, 64

QR factorization and, 278–9

SVD and, 275–6

Rank-deficient LS problem, 288ff

breakdown of QR method, 264

Rank-revealing decomposition, 280–3

Rank-structured matrices, 691ff

Rayleigh quotient iteration, 453–4

QR algorithm and, 464

symmetric-definite pencils and, 501

R-bidiagonalization, 285

Re, 13

Real Schur decomposition, 376–7

generalized, 407

ordering in, 396–7

Rectangular LU, 118

Recursive algorithms

block Cholesky, 169

Strassen, 30–1

Reducible, 373

Regularized least squares, 307ff

Regularized total least squares, 324

Relative error, 69

Relaxation parameter, 619–20

Reorthogonalization

complete, 564

selective, 565

Representation, 681–2

generator, 693

Givens, 697–8

quasiseparable, 694

reshape, 28, 711

and Kronecker product, 28

Residuals vs. accuracy, 138

Restarting

Arnoldi method and, 581–2

block Lanczos and, 569

GMRES and, 644

Restricted generalized SVD, 507

Restricted total least squares, 324

Restriction matrix, 673

Ricatti equation, 422-3

Ridge regression, 307–8

Riemann-Stieltjes integral, 556–7

Right eigenvector, 349

Right-looking, 117

Ritz acceleration, orthogonal

iteration and, 464–5

Ritz approximation

eigenvalues, 551-2

singular values, 573

Rook pivoting, 133

Rotation of subspaces, 327–8

Rotation plus rank-1 (ROPR), 332

Rounding errors. See under particular algorithm

Roundoff error analysis, 100

dot product, 98–9

Wilkinson quote, 99

Row orientation, 5

Row partition, 6

Row scaling, 139

Row weighting in LS problem, 304–5

Saddle point preconditioners, 661

Saxpy, 4, 11

Scaling, linear systems and, 138–9

Scaling and squaring for exp(A), 531

Schur complement, 118–9, 663

Schur decomposition, 67, 350–1

generalized, 406–7

matrix functions and, 523–4

normal matrices and, 351

real matrices and, 376–7

symmetric matrices and, 440

2-by-2 symmetric, 478

Schur vectors, 351

Secular equations, 313–4

Selective reorthogonalizaton, 565–6

Semidefinite systems, 167–8

Semiseparable

eigenvalue problem, 703–4

LU factorization, 695–8

matrix, 682

plus diagonal, 694

QR factorization, 698–701

Sensitivity. See Perturbation results

sep

symmetric matrices and, 444

unsymmetric matrices and, 360

Shared-memory systems, 54–6

Shared-memory traffic, 55–6

Sherman-Morrison formula, 65

Sherman-Morrison-Woodbury formula, 65

Shifts in

QZ algorithm, 411

SVD algorithm, 489

symmetric QR algorithm, 461–2

unsymmetric QR algorithm, 385–90

Sign function, 536–8

Similar matrices, 67, 349

Similarity transformation, 349

condition of, 354

nonunitary, 352–4

Simpson’s rule, 528

Simultaneous diagonalization, 499

Simultaneous iteration. See orthogonal iteration

Sine of matrix, 526

Singular matrix, 65

Singular subspace pair, 488

Singular value decomposition (SVD), 76–80

algorithm for, 488–92

constrained least squares and, 313–4

generalized, 309–10

geometry of, 77

higher-order, 732–3

Jacobi algorithm for, 492–3

Lanczos method for, 571ff

linear systems and, 87–8

minimum-norm least squares solution, 288–9

nullspace and, 78

numerical rank and, 275–6

perturbation of, 487–8

principal angles and, 329–31

projections and, 82

pseudo-inverse and, 290

rank of matrix and, 78

ridge regression and, 307–8

subset selection and, 293–6

subspace intersection and, 331

subspace rotation and, 327–8

symmetric eigenproblem and, 486

total least squares and, 321–2

truncated, 291

Singular values, 76

condition and, 88

eigenvalues and, 355

interlacing property, 48

minimax characterization, 487

perturbation of, 487–8

range and nullspace, 78

rank and, 78

smallest, 279–80

Singular vectors, 76

Skeel condition number, 91

and iterative improvement, 140

Skew-Hamiltonian matrix 420

Skew-Hermitian matrix, 18

Skew-symmetric matrix, 18

span, 64

Sparse factorization challenges

Cholesky, 601

QR, 607

Sparsity, 154

graphs and, 601–2

Spectral abscissa, 349

Spectral radius, 349, 427, 614

Spectrum of matrix, 348

Speed-up, 53–4

SPIKE framework, 199–201

Splitting, 613

Square root of a matrix, 163

s-step Lanczos, 569

Stable algorithm, 136

Stable matrix, 436

Steepest descent method, 626–7

Stieltjes matrix, 658

Strassen method, 30–1

error analysis and, 101-2

Strictly diagonally dominant, 155

Stride, 45

Structured rank, 691ff types of, 702

Sturm sequence property, 468–9

Submatrix, 24

Subnormal floating point number, 95

Subordinate norm, 72

Subset selection, 293–5

using QR with column pivoting, 293

Subspace, 64

angles between, 329–31

deflating, 414

distance between, 82–3, 331

dominant, 368

intersection, 331

invariant, 349

nullspace intersection, 328–9

orthogonal projections onto, 82

rotation of, 327–8

Successive over-relaxation (SOR), 619

Sweep, 480

Sylvester equation, 398

generalized, 417

Sylvester law of inertia, 448

Sylvester map, 682

Symmetric-definite eigenproblem, 497–501

Symmetric eigenproblem, 439ff

sparse methods, 546ff

Symmetric indefinite methods

Aasen, 188–90

Diagonal pivoting, 191–2

Parlett-Reid, 187–8

Symmetric matrix, 18

Symmetric pivoting, 165

Symmetric positive definite systems, 163ff

Symmetric semidefinite properties, 167–8

Symmetric successive over-relaxation, (SSOR), 620

SYMMLQ, 641

Symplectic matrix, 29, 420

symSchur, 478

Taylor approximation of eA, 530

Taylor series, matrix functions and, 515–7

Tensor

contractions, 726ff

eigenvalues, 740–1

networks, 741

notation, 721

rank, 738–9

rank-1, 725

singular values, 739–40

train, 741–3

transpose, 722-3

unfoldings, 720

Thin CS decomposition, 84

Thin QR factorization, 248

Thin SVD, 80

Threshold Jacobi, 483

Tikhonov regularization, 309

Toeplitz-like matrix, 688

Toeplitz matrix methods, classical, 208ff

Toroidal network, 58

Total least squares, 320ff

geometry, 323–4

Tournament pivoting, 150

Trace, 348–9

tr, 348

Trace-min method, 595

Tracy-Singh product, 709

Transition probability matrix, 374

Transpose, 2, 711-2

Trench algorithm, 213

Treppeniteration, 369

Triangular matrices,

multiplication between, 15

unit, 110

Triangular systems, 106–11

band, 177-8

nonsquare, 109–10

roundoff and, 124–5

semiseparable, 694–5

Tridiagonalization,

connection to bidiagonalization, 574

Householder, 458–60

Krylov subspaces and, 459–60

Lanczos, 548–9

Tridiagonal matrices, 15, 223–4

QR algorithm and, 460–4

Tridiagonal systems, 180–1

Truncated

higher-order SVD, 734

SVD, 291

total least squares, 324

Tucker approximation problem, 734–5

ULV decomposition, 282–3

ULV updating, 341–3

Underdetermined systems, 134, 299-301

Undirected graph, 602

Unfolding, 723–4

Unit roundoff, 96

Unit stride, 45

Unit vector, 69

Unitary matrix, 80

Unreduced Hessenberg matrices, 381

Unreduced tridiagonal matrices, 459

Unstable eigenvalue, 363

Unsymmetric

eigenproblem, 347ff

Lanczos method, 584–7

positive definite systems, 161–3

Toeplitz systems, 216–7

Cholesky, 338–41

QR factorization, 334–8

ULV, 341–3

UTV, 282

Vandermonde systems, 203ff

confluent, 206

V-cycle, 677–8

vec, 28, 710–11

for tensors, 722

Vector

computing, 43ff

loads and stores, 43

norms, 68

operations, 3–4, 44

processing, 43

Vectorization, tridiagonal system solving and, 181

Weighted Jacobi iteration, 672–3

Weighting least squares problems column, 306–7

row, 304–5

See also Scaling

Well-conditioned matrix, 88

Wielandt-Hoffman theorem eigenvalues, 442

singular values, 487

Wilkinson shift, 462–3

Work

least squares methods and, 293

linear system methods and, 298

SVD and, 493

WY representation, 238–9

compact version, 244

Yule-Walker problem, 201–10
