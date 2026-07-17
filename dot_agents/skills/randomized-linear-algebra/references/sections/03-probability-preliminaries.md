# Probability preliminaries {#sec:prob}

This section summarizes the key definitions and background from
probability and high-dimensional probability. Later in the survey, we
will present more complete statements of foundational results, as they
are needed.

provide an accessible overview of applied probability. introduces the
field of high-dimensional probability. For more mathematical
presentations, see the classic book of or the lecture notes of .

## Basics

We work in a master probability space that is rich enough to support all
of the random variables that are defined. We will not comment further
about the underlying model.

In this paper, the unqualified term *random variable* encompasses random
scalars, vectors, and matrices. Scalar-valued random variables are
usually (but not always) denoted by uppercase italic Roman letters
($X, Y, Z$). A random vector is denoted by a lowercase bold letter
($\bm{\mathsf{x}}, \bm{\mathsf{\omega}}$). A random matrix is denoted by
an uppercase bold letter
($\bm{\mathsf{X}}, \bm{\mathsf{Y}}, \bm{\mathsf{\Gamma}}, \bm{\mathsf{\Omega}}$).
This notation works in concert with the notation for deterministic
vectors and matrices.

The map $\mathbb{P}(E)$ returns the probability of an event $E$. We
usually specify the event using the compact set builder notation that is
standard in probability. For example,
$\mathbb{P}\left\{  X > t  \right\}$ is the probability that the scalar
random variable $X$ exceeds a level $t$.

The operator $\operatorname{\mathbb{E}}$ returns the expectation of a
random variable. For vectors and matrices, the expectation can be
computed coordinate by coordinate. The expectation is linear, which
justifies relations like
$$\operatorname{\mathbb{E}}[ \bm{\mathsf{AX}} ] = \bm{\mathsf{A}} \operatorname{\mathbb{E}}[ \bm{\mathsf{X}} ]
\quad\text{when $\bm{\mathsf{A}}$ is deterministic and $\bm{\mathsf{X}}$ is random.}$$
We use the convention that nonlinear functions bind before the
expectation; for instance,
$\operatorname{\mathbb{E}}X^2 = \operatorname{\mathbb{E}}[ X^2 ]$. The
operator $\operatorname{Var}[\cdot]$ returns the variance of a scalar
random variable.

We say that a random variable is *centered* when its expectation equals
zero. A random vector $\bm{\mathsf{x}}$ is *isotropic* when
$\operatorname{\mathbb{E}}[ \bm{\mathsf{xx}}^* ] = \bm{\mathsf{I}}$. A
random vector is *standardized* when it is both centered and isotropic.
In particular, a scalar random variable is standardized when it has
expectation zero and variance one.

When referring to independent random variables, we often include the
qualification "statistically independent" to make a distinction with
"linearly independent." We abbreviate the term (statistically)
"independent and identically distributed" as i.i.d.

## Distributions

To refer to a named distribution, we use small capitals. In this
context, the symbol $\sim$ means "has the same distribution as."

We write $\textsc{unif}$ for the uniform distribution over a finite set
(with counting measure). In particular, a scalar Rademacher random
variable has the distribution $\textsc{unif}\{ \pm 1 \}$. A Rademacher
random vector has iid coordinates, each distributed as a scalar
Rademacher random variable. We sometimes require the uniform
distribution over a Borel subset of $\mathbb{F}^n$, equipped with
Lebesgue measure.

We write $\textsc{normal}(\bm{\mathsf{\mu}}, \bm{\mathsf{C}})$ for the
normal distribution on $\mathbb{F}^n$ with expectation
$\bm{\mathsf{\mu}} \in \mathbb{F}^n$ and psd covariance matrix
$\bm{\mathsf{C}} \in \mathbb{H}_n(\mathbb{F})$. A *standard normal*
random variable or random vector has expectation zero and covariance
matrix equal to the identity. We often use the term *Gaussian* to refer
to normal distributions.

## Concentration inequalities

Concentration inequalities provide bounds on the probability that a
random variable is close to its expectation. A good reference for the
scalar case is the book by . In the matrix setting, closeness is
measured in the spectral norm. For an introduction to matrix
concentration, see
[@Tro15:Introduction-Matrix; @Tro19:Matrix-Concentration-LN]. These
results play an important role in randomized linear algebra.

## Gaussian random matrix theory

On several occasions, we use comparison principles to study the action
of a random matrix with iid Gaussian entries. In particular, these
methods can be used to control the largest and smallest singular values.
The main classical comparison theorems are associated with the names
Slepian, Chevet, and Gordon. For accounts of this, see or or . More
recently, it has been observed that Gordon's inequality can be reversed
in certain settings [@TOH14:Gaussian-Minmax].

In several instances, we require more detailed information about
Gaussian random matrices. General resources include and . Most of the
specific results we need are presented in .

