# 12.5.6 Tensor Singular Values: A Variational Approach

The singular values of a matrix $A \in \mathbb { R } ^ { n _ { 1 } \times n _ { 2 } }$ are the stationary values of

$$
\psi_ {A} (u, v) = \frac {u ^ {T} A v}{\| u \| _ {2} \| v \| _ {2}} = \frac {\sum_ {i _ {1} = 1} ^ {n _ {1}} \sum_ {i _ {2} = 1} ^ {n _ {2}} A (i _ {1} , i _ {2}) u (i _ {1}) v (i _ {2})}{\| u \| _ {2} \| v \| _ {2}} \tag {12.5.22}
$$

and the associated stationary vectors are the corresponding singular vectors. This follows by looking at the gradient equation $\nabla \psi ( u , v ) = 0$ . Indeed, if u and v are unit vectors, then this equation has the form

$$
\nabla \psi_ {A} (u, v) = \left[ \begin{array}{c} A v - \psi_ {A} (u, v) u \\ A ^ {T} u - \psi_ {A} (u, v) v \end{array} \right] = 0.
$$

This variational characterization of matrix singular values and vectors extends to tensors; see Lim (2005). Suppose $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times n _ { 2 } \times n _ { 3 } }$ and define

$$
\psi_ {\mathcal {A}} (u _ {1}, u _ {2}, u _ {3}) = \frac {\sum_ {\mathbf {i} = \mathbf {1}} ^ {\mathbf {n}} \mathcal {A} (\mathbf {i}) \cdot u _ {1} (i _ {1})   u _ {2} (i _ {2})   u _ {3} (i _ {3})}{\parallel u _ {1} \parallel_ {2} \parallel u _ {2} \parallel_ {2} \parallel u _ {3} \parallel_ {2}}
$$

where $u _ { 1 } \in \mathbb { R } ^ { n _ { 1 } } , u _ { 2 } \in \mathbb { R } ^ { n _ { 2 } }$ , and $u _ { 3 } \in \mathbb { R } ^ { n _ { 3 } }$ . It is easy to show that

$$
\psi_ {\mathcal {A}} (u _ {1}, u _ {2}, u _ {3}) = \left\{ \begin{array}{l} u _ {1} ^ {T} \mathcal {A} _ {(1)} (u _ {3} \otimes u _ {2}) / (\| u _ {1} \| _ {2} \| u _ {2} \| _ {2} \| u _ {3} \| _ {2}), \\ u _ {2} ^ {T} \mathcal {A} _ {(2)} (u _ {3} \otimes u _ {1}) / (\| u _ {1} \| _ {2} \| u _ {2} \| _ {2} \| u _ {3} \| _ {2}), \\ u _ {3} ^ {T} \mathcal {A} _ {(3)} (u _ {2} \otimes u _ {1}) / (\| u _ {1} \| _ {2} \| u _ {2} \| _ {2} \| u _ {3} \| _ {2}). \end{array} \right.
$$

If $u _ { 1 } , u _ { 2 }$ , and $u _ { 3 }$ are unit vectors, then the equation $\nabla \psi _ { A } = 0$ i s

$$
\nabla \psi_ {\mathcal {A}} = \left[ \begin{array}{l} \mathcal {A} _ {(1)} (u _ {3} \otimes u _ {2}) \\ \mathcal {A} _ {(2)} (u _ {3} \otimes u _ {1}) \\ \mathcal {A} _ {(3)} (u _ {2} \otimes u _ {1}) \end{array} \right] - \psi_ {\mathcal {A}} (u _ {1}, u _ {2}, u _ {3}) \left[ \begin{array}{l} u _ {1} \\ u _ {2} \\ u _ {3} \end{array} \right] = 0.
$$

If we can satisfy this equation, then we will call $\psi _ { \mathcal { A } } ( u _ { 1 } , u _ { 2 } , u _ { 3 } )$ a singular value of the tensor A. If we take a componentwise approach to this this nonlinear system we are led to the following iteration

# Repeat:

$$
\begin{array}{l} \tilde {u} _ {1} = \mathcal {A} _ {(1)} (u _ {3} \otimes u _ {2}), \quad u _ {1} = \tilde {u} _ {1} / \| \tilde {u} _ {1} \| _ {2} \\ \tilde {u} _ {2} = \mathcal {A} _ {(2)} (u _ {3} \otimes u _ {1}), \quad u _ {2} = \tilde {u} _ {2} / \| \tilde {u} _ {2} \| _ {2} \\ \tilde {u} _ {3} = \mathcal {A} _ {(3)} (u _ {2} \otimes u _ {1}), \quad u _ {3} = \tilde {u} _ {3} / \| \tilde {u} _ {3} \| _ {2} \\ \sigma = \psi (u _ {1}, u _ {2}, u _ {3}) \\ \end{array}
$$

This can be thought of as a higher-order power iteration. Upon comparison with the Tucker approximation problem with $\mathbf { r } = [ 1 , 1 , \ldots , 1 ]$ , we see that it is a strategy for computing a nearest rank-1 tensor.
