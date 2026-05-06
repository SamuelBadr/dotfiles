# 12.5.4 The CP Approximation Problem

A nice attribute of the matrix SVD that is that the “core matrix” in the rank-1 expansion is diagonal. This is not true when we graduate to tensors and work with the

Tucker representation. However, there is an alternate way to extrapolate from the matrix SVD if we prefer “diagonalness” to orthogonality. Given $\mathcal { X } \in \mathbf { \bar { R } } ^ { n _ { 1 } \times n _ { 2 } \times n _ { 3 } }$ and an integer r, we consider the problem

$$
\min _ {\mathcal {X}} \| \mathcal {A} - \mathcal {X} \| _ {F} \tag {12.5.14}
$$

such that

$$
\mathcal {X} = \sum_ {j = 1} ^ {r} \lambda_ {j} \cdot F (:, j) \circ G (:, j) \circ H (:, j) \tag {12.5.15}
$$

where $F \in \mathbb { R } ^ { n _ { 1 } \times r } , G \in \mathbb { R } ^ { n _ { 2 } \times r }$ , and $H \in \mathbb { R } ^ { n _ { 3 } \times r }$ . This is an example of the CP approximation problem. We assume that the columns of F , G, and H have unit 2-norm.

The modal unfoldings of the tensor (12.5.15) are neatly characterized through the Khatri-Rao product that we defined in §12.3.3. If

$$
F = \left[ f _ {1} \mid \dots \mid f _ {r} \right], \qquad G = \left[ g _ {1} \mid \dots \mid g _ {r} \right], \qquad H = \left[ h _ {1} \mid \dots \mid h _ {r} \right],
$$

then

$$
\begin{array}{l} \mathcal {X} _ {(1)} = \sum_ {j = 1} ^ {r} \lambda_ {j} \cdot f _ {j} \otimes (h _ {j} \otimes g _ {j}) ^ {T} = F \cdot \operatorname{diag} (\lambda_ {j}) \cdot (H \odot G) ^ {T}, \\ \mathcal {X} _ {(2)} = \sum_ {j = 1} ^ {r} \lambda_ {j} \cdot g _ {j} \otimes (h _ {j} \otimes f _ {j}) ^ {T} = G \cdot \operatorname{diag} (\lambda_ {j}) \cdot (H \odot F) ^ {T}, \\ \mathcal {X} _ {(3)} = \sum_ {j = 1} ^ {r} \lambda_ {j} \cdot h _ {j} \otimes (g _ {j} \otimes f _ {j}) ^ {T} = H \cdot \operatorname{diag} (\lambda_ {j}) \cdot (G \odot F) ^ {T}. \\ \end{array}
$$

These results follow from the previous section. For example,

$$
\begin{array}{l} \mathcal {X} _ {(1)} = \sum_ {j = 1} ^ {r} \lambda_ {j} \left(f _ {j} \circ g _ {j} \circ h _ {j}\right) _ {(1)} = \sum_ {j = 1} ^ {r} \lambda_ {j} f _ {j} (h _ {j} \otimes g _ {j}) ^ {T} \\ = \left[ \lambda_ {1} f _ {1} \mid \dots \mid \lambda_ {r} f _ {r} \right] \left[ h _ {1} \otimes g _ {1} \mid \dots \mid h _ {r} \otimes g _ {r} \right] ^ {T} = F \cdot \operatorname{diag} (\lambda_ {j}) \cdot (H \odot G) ^ {T}. \\ \end{array}
$$

Noting that

$$
\left\| \mathcal {A} - \mathcal {X} \right\| _ {F} = \left\| \mathcal {A} _ {(1)} - \mathcal {X} _ {(1)} \right\| _ {F} = \left\| \mathcal {A} _ {(2)} - \mathcal {X} _ {(2)} \right\| _ {F} = \left\| \mathcal {A} _ {(3)} - \mathcal {X} _ {(3)} \right\| _ {F},
$$

we see that the CP approximation problem can be solved by minimizing any one of the following expressions:

$$
\left\| \mathcal {A} _ {(1)} - \mathcal {X} _ {(1)} \right\| _ {F} = \left\| \mathcal {A} _ {(1)} - F \cdot \mathrm{diag} (\lambda_ {j}) \cdot (H \odot G) ^ {T} \right\| _ {F}, \tag {12.5.16}
$$

$$
\left\| \mathcal {A} _ {(2)} - \mathcal {X} _ {(2)} \right\| _ {F} = \left\| \mathcal {A} _ {(2)} - G \cdot \mathrm{diag} (\lambda_ {j}) \cdot (H \odot F) ^ {T} \right\| _ {F}, \tag {12.5.17}
$$

$$
\left\| \mathcal {A} _ {(3)} - \mathcal {X} _ {(3)} \right\| _ {F} = \left\| \mathcal {A} _ {(3)} - H \cdot \mathrm{diag} (\lambda_ {j}) \cdot (G \odot F) ^ {T} \right\| _ {F}. \tag {12.5.18}
$$

This is a multilinear least squares problem. However, observe that if we fix λ, H, and G in (12.5.16), then $\| \mathcal { A } _ { ( 1 ) } - \mathcal { X } _ { ( 1 ) } \| _ { F }$ is linear in F . Similar comments apply to (12.5.17) and (12.5.18) and we are led to the following alternating least squares minimization strategy:

# Repeat:

$\mathrm { L e t } ~ \tilde { F } \mathrm { ~ m i n i m i z e } \parallel \mathcal { A } _ { ( 1 ) } - \tilde { F } \cdot ( H \odot G ) ^ { T } \parallel _ { _ { F } } \mathrm { ~ a n d ~ f o r } \quad j = 1 \colon r \mathrm { ~ s e t }$

$$
\lambda_ {j} = \parallel \tilde {F} (:, j) \parallel_ {2} \text { and } F (:, j) = \tilde {F} (:, j) / \lambda_ {j}.
$$

$\mathrm { L e t } \ \tilde { G } \ \mathrm { m i n i m i z e } \ \lVert \ A _ { ( 2 ) } - \tilde { G } \cdot ( H \odot F ) ^ { T } \ \rVert _ { F } \quad \mathrm { a n d ~ f o r } \quad j = 1 \colon r \quad \mathrm { s e t }$

$$
\lambda_ {j} = \parallel \tilde {G} (:, j) \parallel_ {2} \text { and } G (:, j) = \tilde {G} (:, j) / \lambda_ {j}.
$$

$\mathrm { L e t } ~ \tilde { H } \mathrm { ~ m i n i m i z e } ~ \lVert ~ \mathcal { A } _ { ( 3 ) } - \tilde { H } \cdot ( G \odot F ) ^ { T } ~ \rVert _ { \epsilon } \quad \mathrm { a n d ~ f o r } \quad j = 1 \colon r \mathrm { ~ s e t }$

$$
\lambda_ {j} = \left\| \tilde {H} (:, j) \right\| _ {2} \text { and } H (:, j) = \tilde {H} (:, j) / \lambda_ {j}.
$$

The update calculations for F , G, and H are highly structured linear least squares problems. The central calculations involve linear least square problems of the form

$$
\min \left\| (B \odot C) z - d \right\| _ {2} \tag {12.5.19}
$$

where $B \in \mathbb { R } ^ { p _ { B } \times q } , ~ C \in \mathbb { R } ^ { p _ { C } \times q }$ , and $d \in \mathbb { R } ^ { p _ { B } p _ { C } }$ . This is typically a “tall skinny” LS problem. If we form the Khatri-Rao product and use the QR factorization in the usual way, then $O ( p _ { B } p _ { C } q ^ { 2 } )$ flops are required to compute z. On the other hand, the normal equation system corresponding to (12.5.19) is

$$
\left((B ^ {T} B). * (C ^ {T} C)\right) z = (B \odot C) ^ {T} d \tag {12.5.20}
$$

which can be formed and solved via the Cholesky factorization in $O ( ( p _ { B } + p _ { C } ) q ^ { 2 } )$ flops. For general tensors $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ there are d least squares problems to solve per pass. In particular, given A and r, the CP approximation problem involves finding matrices

$$
F ^ {(k)} = [ f _ {1} ^ {(k)} \mid \dots \mid f _ {r} ^ {(k)} ] \in \mathbb {R} ^ {n _ {k} \times r}, \qquad k = 1: d,
$$

with unit 2-norm columns and a vector $\lambda \in \mathbb { R } ^ { r }$ so that if

$$
\mathcal {X} = \sum_ {j = 1} ^ {r} \lambda_ {j} f _ {j} ^ {(1)} \circ \dots \circ f _ {j} ^ {(d)}, \tag {12.5.21}
$$

then $\| \mathcal { A } - \mathcal { X } \| _ { F }$ is minimized. Noting that

$$
\mathcal {X} _ {(k)} = F ^ {(k)} \mathrm{diag} (\lambda) \left(F ^ {(d)} \odot \dots \odot F ^ {(k + 1)} \odot F ^ {(k - 1)} \odot \dots \odot F ^ {(1)}\right) ^ {T},
$$

we obtain the following iteration.

#

$$
\lambda_ {j} = \left\| \tilde {F} _ {(k)} (:, j) \right\| _ {2}
$$

$$
F ^ {(k)} (:, j) = \tilde {F} _ {k} (:, j) / \lambda_ {j}
$$

end

end

This is the CANDECOMP/PARAFAC framework. For implementation details about this nonlinear iteration, see Smilde, Bro, and Geladi (2004, pp. 113–119) and Kolda and Bader (2009).
