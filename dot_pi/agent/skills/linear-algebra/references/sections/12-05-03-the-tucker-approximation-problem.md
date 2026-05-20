# 12.5.3 The Tucker Approximation Problem

Suppose $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times n _ { 2 } \times n _ { 3 } }$ and assume that $\mathbf { r } \leq \mathbf { r a n k } _ { * } ( \mathcal { A } )$ with inequality in at least one component. Prompted by the optimality properties of the matrix SVD, let us consider the following optimization problem:

$$
\min _ {\mathcal {X}} \| \mathcal {A} - \mathcal {X} \| _ {F} \tag {12.5.12}
$$

such that

$$
\mathcal {X} = \sum_ {\mathbf {j} = \mathbf {1}} ^ {\mathbf {r}} \mathcal {S} (\mathbf {j}) \cdot U _ {1} (:, j _ {1}) \circ U _ {2} (:, j _ {2}) \circ U _ {3} (:, j _ {3}). \tag {12.5.13}
$$

We refer to this as the Tucker approximation problem. Unfortunately, the truncated HOSVD tensor $\mathcal { A } ^ { ( \mathbf { r } ) }$ does not solve the Tucker approximation problem, prompting us to develop an appropriate optimization strategy.

To be clear, we are given A and r and seek a core tensor S that is $r _ { 1 } – \mathrm { b y } – r _ { 2 } – \mathrm { b y } – r _ { 3 }$ and matrices $U _ { 1 } \in \mathbb { R } ^ { n _ { 1 } \times r _ { 1 } } , U _ { 2 } \in \mathbb { R } ^ { n _ { 2 } \times r _ { 2 } }$ , and $U _ { 3 } \in \mathbb { R } ^ { n _ { 3 } \times r _ { 3 } }$ with orthonormal columns so that the tensor X defined by (12.5.13) solves (12.5.12). Using Theorem 12.4.1 we know that

$$
\left\| \mathcal {A} - \mathcal {X} \right\| _ {F} = \left\| \operatorname{vec} (\mathcal {A}) - \left(U _ {3} \otimes U _ {2} \otimes U _ {1}\right) \cdot \operatorname{vec} (\mathcal {S}) \right\| _ {2}.
$$

Since $U _ { 3 } \otimes U _ { 2 } \otimes U _ { 1 }$ has orthonormal columns, it follows that the “best” S given any triplet $\{ U _ { 1 } , U _ { 2 } , U _ { 3 } \}$ is

$$
\mathcal {S} = \left(U _ {3} ^ {T} \otimes U _ {2} ^ {T} \otimes U _ {1} ^ {T}\right) \cdot \operatorname{vec} (\mathcal {A}).
$$

Thus, we can remove $s$ from the search space and simply look for $U = U _ { 3 } \otimes U _ { 2 } \otimes U _ { 1 }$ so that

$$
\left\| \left(I - U U ^ {T}\right) \cdot \operatorname{vec} (\mathcal {A}) \right\| _ {F} ^ {2} = \left\| \operatorname{vec} (\mathcal {A}) \right\| _ {F} ^ {2} - \left\| U ^ {T} \cdot \operatorname{vec} (\mathcal {A}) \right\| _ {F} ^ {2}
$$

is minimized. In other words, determine $U _ { 1 } , U _ { 2 }$ , and $U _ { 3 }$ so that

$$
\| \left(U _ {3} ^ {T} \otimes U _ {2} ^ {T} \otimes U _ {1} ^ {T}\right) \cdot \mathsf {v e c} (\mathcal {A}) \| _ {F} = \left\{ \begin{array}{l l} \| U _ {1} ^ {T} \cdot A _ {(1)} \cdot (U _ {3} \otimes U _ {2}) \| _ {F} \\ \| U _ {2} ^ {T} \cdot A _ {(2)} \cdot (U _ {3} \otimes U _ {1}) \| _ {F} \\ \| U _ {3} ^ {T} \cdot A _ {(3)} \cdot (U _ {2} \otimes U _ {1}) \| _ {F} \end{array} \right.
$$

is maximized. By freezing any two of the three matrices $\{ U _ { 1 } , U _ { 2 } , U _ { 3 } \}$ we can improve the third by solving an optimization problem of the form (12.5.3). This suggests the following strategy:

# Repeat:

Maximize $\parallel \boldsymbol { U } _ { 1 } ^ { T } \cdot \boldsymbol { A } _ { ( 1 ) } \cdot \left( \boldsymbol { U } _ { 3 } \otimes \boldsymbol { U } _ { 2 } \right) \parallel _ { F }$ with respect to $U _ { 1 }$ by computing the

$$
\text { SVD } \mathcal {A} _ {(1)} \cdot (U _ {3} \otimes U _ {2}) = \tilde {U} _ {1} \Sigma_ {1} V _ {1} ^ {T}. \text { Set } U _ {1} = \tilde {U} _ {1} (:, 1: r _ {1}).
$$

Maximize $\parallel U _ { 2 } ^ { T } \cdot A _ { ( 2 ) } \cdot ( U _ { 3 } \otimes U _ { 1 } ) \parallel _ { F }$ with respect to $U _ { 2 }$ by computing the

$$
\text { SVD } \mathcal {A} _ {(2)} \cdot (U _ {3} \otimes U _ {1}) = \tilde {U} _ {2} \Sigma_ {2} V _ {2} ^ {T}. \text { Set } U _ {2} = \tilde {U} _ {2} (:, 1: r _ {2}).
$$

Maximize  $U _ { 3 } ^ { T } \cdot A _ { ( 3 ) } \cdot ( U _ { 2 } \otimes U _ { 1 } ) \parallel _ { F }$ with respect to $U _ { 3 } { \mathrm { : } }$ : by computing the

$$
\text { SVD } \mathcal {A} _ {(3)} \cdot (U _ {2} \otimes U _ {1}) = \tilde {U} _ {3} \Sigma_ {3} V _ {3} ^ {T}. \text { Set } U _ {3} = \tilde {U} _ {3} (:, 1: r _ {3}).
$$

This is an example of the alternating least squares framework. For order-d tensors, there are d optimizations to perform each step:

# Repeat:

for k = 1:d

Compute the SVD:

$$
\mathcal {A} _ {(k)} \left(U _ {d} \otimes \dots \otimes U _ {k + 1} \otimes U _ {k - 1} \otimes \dots \otimes U _ {1}\right) = \tilde {U} _ {k} \Sigma_ {k} V _ {k} ^ {T}.
$$

$$
U _ {k} = \tilde {U} _ {k} (:, 1: r _ {k})
$$

end

This is essentially the Tucker framework. For implementation details concerning this nonlinear iteration, see De Lathauwer, De Moor, and Vandewalle (2000b), Smilde, Bro, and Geladi (2004, pp. 119–123), and Kolda and Bader (2009).
