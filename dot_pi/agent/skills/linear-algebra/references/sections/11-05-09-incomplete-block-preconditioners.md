# 11.5.9 Incomplete Block Preconditioners

The incomplete factorization idea can be applied at the block level. For example, an incomplete block Cholesky factor $H = \left( H _ { i j } \right)$ of a block symmetric positive definite matrix $A = \left( A _ { i j } \right)$ could be obtained by forcing $H _ { i j }$ to be zero if $A _ { i j }$ is zero. However, there is another level of opportunity if the individual $A _ { i j }$ are themselves sparse, for then it may be necessary to impose constraints on the sparsity structure of the $H _ { i j }$ .

To illustrate this in a simple familiar setting, let us build an incomplete Cholesky factorization for a block tridiagonal matrix whose diagonal blocks are tridiagonal and whose subdiagonal and superdiagonal blocks are diagonal. (The §4.8.3 model problem matrices have this structure.) With

$$
A = \left[ \begin{array}{c c c} A _ {1} & E _ {1} ^ {T} & 0 \\ E _ {1} & A _ {2} & E _ {2} ^ {T} \\ 0 & E _ {2} & A _ {3} \end{array} \right] = \left[ \begin{array}{c c c} G _ {1} & 0 & 0 \\ F _ {1} & G _ {2} & 0 \\ 0 & F _ {2} & G _ {3} \end{array} \right] \left[ \begin{array}{c c c} G _ {1} ^ {T} & F _ {1} ^ {T} & 0 \\ 0 & G _ {2} ^ {T} & F _ {2} ^ {T} \\ 0 & 0 & G _ {3} ^ {T} \end{array} \right],
$$

here are the recipes for the $G _ { k }$ and $F _ { k }$ if A is $p { \mathrm { - } } \mathrm { b y } { \mathrm { - } } p$ as a block matrix:

$$
G _ {1} G _ {1} ^ {T} = A _ {1}
$$

for $k = 1 { : } p - 1$

$$
F _ {k} = E _ {k} G _ {k} ^ {- T}
$$

$$
G _ {k + 1} G _ {k + 1} ^ {T} = A _ {k + 1} - E _ {k} (G _ {k} G _ {k} ^ {T}) ^ {- 1} E _ {k} ^ {T}
$$

end

Except for $G _ { 1 }$ , all the Cholesky factor blocks are dense. A way around this difficulty is to replace $( \dot { G } _ { k } G _ { k } ^ { T } ) ^ { - 1 }$ with a suitably chosen tridiagonal approximation $\Lambda _ { k }$ :

$$
\tilde {G} _ {1} \tilde {G} _ {1} ^ {T} = A _ {1}
$$

for $k = 1 { : } p - 1$

$$
\tilde {F} _ {k} = E _ {k} \tilde {G} _ {k} ^ {- T} \tag {11.5.10}
$$

$$
\tilde {G} _ {k + 1} \tilde {G} _ {k + 1} ^ {T} = A _ {k + 1} - E _ {k} \Lambda_ {k} E _ {k} ^ {T}
$$

end

Note that with this strategy, each $\tilde { G } _ { k }$ is lower bidiagonal. The $\tilde { F } _ { k }$ are full, but they do not have to actually be formed in order to solve systems that involve the incomplete factors. For example,

$$
\left[ \begin{array}{c c c} \tilde {G} _ {1} & 0 & 0 \\ \tilde {F} _ {1} & \tilde {G} _ {2} & 0 \\ 0 & \tilde {F} _ {2} & \tilde {G} _ {3} \end{array} \right] \left[ \begin{array}{l} w _ {1} \\ w _ {2} \\ w _ {3} \end{array} \right] = \left[ \begin{array}{l} r _ {1} \\ r _ {2} \\ r _ {3} \end{array} \right], \qquad \begin{array}{l} \tilde {G} _ {1} w _ {1} = r _ {1}, \\ \tilde {G} _ {2} w _ {2} = r _ {2} - E _ {1} \tilde {G} _ {1} ^ {- T} w _ {1}, \\ \tilde {G} _ {3} w _ {3} = r _ {3} - E _ {2} \tilde {G} _ {2} ^ {- T} w _ {2}. \end{array}
$$

Each $w _ { k }$ requires a $\tilde { G } _ { k } { \mathrm { - s y s t e m } }$ solution and a G˜Tk -system solution.

There remains the issue of choosing $\Lambda _ { 1 } , \ldots , \Lambda _ { p - 1 }$ . The central problem is how to determine a symmetric tridiagonal Λ so that if $\mathbf { \bar { \boldsymbol { T } } } \in \mathbb { R } ^ { m \times m }$ is symmetric positive definite and tridiagonal itself, then $\Lambda \approx T ^ { - 1 }$ . Possibilities include:

• Let $\Lambda = \mathrm { d i a g } ( 1 / t _ { 1 1 } , \dots , 1 / t _ { m m } )$ .   
• Let Λ be the tridiagonal part of $T ^ { - 1 }$ , an $O ( m )$ computation. See P11.5.5.   
• Let $\Lambda = U ^ { T } U$ where U is the lower bidiagonal portion of $K ^ { - 1 }$ where $T = K K ^ { T }$ is the Cholesky factorization. This is an $O ( m )$ computation. See P11.5.6.

For a discussion of these approximations and what they imply about the associated preconditioners, see Concus, Golub, and Meurant (1985).
