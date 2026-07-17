# 12.4 Tensor Unfoldings and Contractions

An order-d tensor $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ is a real d-dimensional array $\mathcal { A } ( 1 { : } n _ { 1 } , \ldots , 1 { : } n _ { d } )$ where the index range in the kth mode is from 1 to $n _ { k }$ . Low-order examples include scalars (order-0), vectors (order-1), and matrices (order-2). Order-3 tensors can be visualized as “Rubik cubes of data,” although the dimensions do not have to be equal along each mode. For example, $\mathcal { A } \in \mathbb { R } ^ { m \times n \times 3 }$ might house the red, green, and blue pixel data for an m-by-n image, $\mathrm { { a } \ { \tilde { \ s t a c k i n g } } ^ { ; 5 } }$ of three matrices. In many applications, a tensor is used to capture what a multivariate function looks like on a lattice of points, e.g., $\mathcal { A } ( i , j , k , \ell ) \approx f ( w _ { i } , x _ { j } , y _ { k } , z _ { \ell } )$ . The function f could be the solution to a complicated partial differential equation or a general mapping from some high-dimensional space of input values to a measurement that is acquired experimentally.

Because of their higher dimension, tensors are harder to reason about than matrices. Notation, which is always important, is critically important in tensor computations where vectors of subscripts and deeply nested summations are the rule. In this section we examine some basic tensor operations and develop a handy, matrix type of notation that can be used to describe them. Kronecker products are central.

Excellent background references include De Lathauwer (1997), Smilde, Bro, and Geladi (2004), and Kolda and Bader (2009).

# 12.4.1 Unfoldings and Contractions: A Preliminary Look

To unfold a tensor is to systematically arrange its entries into a matrix.3 Here is one possible unfolding of a 2-by-2-by-3-by-4 tensor:

$$
A = \left[ \begin{array}{c c c c c c c} a _ {1 1 1 1} & a _ {1 2 1 1} & a _ {1 1 1 2} & a _ {1 2 1 2} & a _ {1 1 1 3} & a _ {1 2 1 3} & a _ {1 1 1 4} & a _ {1 2 1 4} \\ a _ {2 1 1 1} & a _ {2 2 1 1} & a _ {2 1 1 2} & a _ {2 2 1 2} & a _ {2 1 1 3} & a _ {2 2 1 3} & a _ {2 1 1 4} & a _ {2 2 1 4} \\ \hline a _ {1 1 2 1} & a _ {1 2 2 1} & a _ {1 1 2 2} & a _ {1 2 2 2} & a _ {1 1 2 3} & a _ {1 2 2 3} & a _ {1 1 2 4} & a _ {1 2 2 4} \\ a _ {2 1 2 1} & a _ {2 2 2 1} & a _ {2 1 2 2} & a _ {2 2 2 2} & a _ {2 1 2 3} & a _ {2 2 2 3} & a _ {2 1 2 4} & a _ {2 2 2 4} \\ \hline a _ {1 1 3 1} & a _ {1 2 3 1} & a _ {1 1 3 2} & a _ {1 2 3 2} & a _ {1 1 3 3} & a _ {1 2 3 3} & a _ {1 1 3 4} & a _ {1 2 3 4} \\ a _ {2 1 3 1} & a _ {2 2 3 1} & a _ {2 1 3 2} & a _ {2 2 3 2} & a _ {2 1 3 3} & a _ {2 2 3 3} & a _ {2 1 3 4} & a _ {2 2 3 4} \end{array} \right]
$$

Order-4 tensors are interesting because of their connection to block matrices. Indeed, a block matrix $\boldsymbol { A } = \left( A _ { k \ell } \right)$ with equally sized blocks can be regarded as an order-4 tensor $\mathcal { A } = \left( a _ { i j k \ell } \right)$ where $[ A _ { k \ell } ] _ { i j } = a _ { i j k \ell }$ .

Unfoldings have an important role to play in tensor computations for three reasons. (1) Operations between tensors can often be reformulated as a matrix computation between unfoldings. (2) Iterative multilinear optimization strategies for tensor problems typically involve one or more unfoldings per step. (3) Hidden structures within a tensor dataset can sometimes be revealed by discovering patterns within its unfoldings. For these reasons, it is important to develop a facility with tensor unfoldings because they serve as a bridge between matrix computations and tensor computations

Operations between tensors typically involve vectors of indices and deeply nested loops. For example, here is a matrix-multiplication-like computation that combines two order-4 tensors to produce a third order-4 tensor:

$$
\begin{array}{l} \text {for} i _ {1} = 1: n \\ \text {for} i _ {2} = 1: n \\ \text {for} i _ {3} = 1: n \\ \text {for} i _ {4} = 1: n \\ \mathcal {C} (i _ {1}, i _ {2}, i _ {3}, i _ {4}) = \sum_ {p = 1} ^ {n} \sum_ {q = 1} ^ {n} \mathcal {A} (i _ {1}, p, i _ {3}, q) \mathcal {B} (p, i _ {2}, q, i _ {4}) \\ \text {end} \\ \text {end} \\ \text {end} \\ \text {end} \end{array} \tag {12.4.1}
$$

This is an example of a tensor contraction. Tensor contractions are essentially reshaped, multi-indexed matrix multiplications and can be very expensive to compute. (The above example involves $O ( n ^ { 6 } )$ flops.) It is increasingly common to have $O ( n ^ { d } )$ contraction bottlenecks in a simulation. In order to successfully tap into the “culture” of of high-performance matrix computations, it is important to have an intuition about tensor contractions and how they can be organized.

# 12.4.2 Notation and Definitions

If $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ and $\mathbf { i } = ( i _ { 1 } , \dots , i _ { d } )$ with $1 \leq i _ { k } \leq n _ { k }$ for k = 1:d, then

$$
\mathcal {A} (\mathbf {i}) \equiv \mathcal {A} (i _ {1}, \dots , i _ {k}).
$$

The vector i is a subscript vector. Bold font is used designate subscript vectors while calligraphic font is used for tensors. For low-order tensors we sometimes use matrixstyle subscripting, $\mathrm { e . g . } , \mathcal { A } \ = \ ( a _ { i j k \ell } )$ . It is sometimes instructive to write $\mathbf { \mathcal { A } } ( \mathbf { i } , \mathbf { j } )$ for $\mathcal { A } ( [ \mathbf { i j } ] )$ . Thus,

$$
\mathcal {A} ([ 2 5 3 4 7 ]) = \mathcal {A} (2, 5, 3, 4, 7) = a _ {2 5 3 4 7} = a _ {2 5 3, 4 7} = \mathcal {A} ([ 2, 5, 3 ], [ 4, 7 ])
$$

shows the several ways that we can refer to a tensor entry.

We extend the Matlab colon notation in order to identify subtensors. If L and R are subscript vectors with the same dimension, then $\mathbf { L } \leq \mathbf { R }$ means that $L _ { k } \le R _ { k }$ for all k. The length-d subscript vector of all 1’s is denoted by $\mathbf { 1 } _ { d }$ . If the dimension is clear from the context, then we just write 1. Suppose $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ with $\mathbf { n } = [ n _ { 1 } , \ldots , n _ { d } ]$ . If $\mathbf { 1 } \leq \mathbf { L } \leq \mathbf { R } \leq \mathbf { n }$ , then $\scriptstyle A ( \mathbf { L } : \mathbf { R } )$ denotes the subtensor

$$
B = \mathcal {A} (L _ {1}: R _ {1}, \ldots , L _ {d}: R _ {d}).
$$

Just as we can extract an order-1 tensor from an order-2 tensor, $\mathrm { e . g . } , A ( : , k )$ , so can we extract a lower-order tensor from a given tensor. Thus, if $\mathcal { A } \in \mathbb { R } ^ { \bar { 2 } \times 3 \times 4 \times 5 }$ , then

(i) $\begin{array} { r } { \mathcal { B } = \mathcal { A } ( 1 , : , 2 , 4 ) \in \mathbb { R } ^ { 3 } \qquad \Rightarrow \quad \mathcal { B } ( i _ { 2 } ) = \mathcal { A } ( 1 , i _ { 2 } , 2 , 4 ) , } \end{array}$   
(ii) $\begin{array} { r } { \mathcal { B } = \mathcal { A } ( 1 , : , 2 , : ) \in \mathbb { R } ^ { 3 \times 5 } \qquad \Rightarrow \quad \mathcal { B } ( i _ { 2 } , i _ { 4 } ) = \mathcal { A } ( 1 , i _ { 2 } , 2 , i _ { 4 } ) , } \end{array}$   
$\mathrm { ( i i i ) } \quad B = A ( : , : , 2 , : ) \in \mathbb { R } ^ { 2 \times 3 \times 5 } \quad \Rightarrow \quad \mathcal { B } ( i _ { 1 } , i _ { 2 } , i _ { 4 } ) = A ( i _ { 1 } , i _ { 2 } , 2 , i _ { 4 } ) .$

Order-1 extractions like (i) are called fibers. Order-2 extractions like (ii) are called slices. More general extractions like (iii) are called subtensors.

It is handy to have a multi-index summation notation. If n is a length-d index vector, then

$$
\sum_ {\mathbf {i} = \mathbf {1}} ^ {\mathbf {n}} \equiv \sum_ {i _ {1} = 1} ^ {n _ {1}} \dots \sum_ {i _ {d} = 1} ^ {n _ {d}}.
$$

Thus, if $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ , then its Frobenius norm is given by

$$
\| \mathcal {A} \| _ {F} = \sqrt {\sum_ {\mathbf {i} = 1} ^ {\mathbf {n}} \mathcal {A} (\mathbf {i}) ^ {2}}.
$$

# 12.4.3 The Vec Operation for Tensors

As with matrices, the vec(·) operator turns tensors into column vectors, e.g.,

$$
\mathcal {A} \in \mathbb {R} ^ {2 \times 3 \times 2} \qquad \Longrightarrow \qquad \operatorname{vec} (\mathcal {A}) = \left[ \begin{array}{c} \mathcal {A} (:, 1, 1) \\ \hline \mathcal {A} (:, 2, 1) \\ \hline \mathcal {A} (:, 3, 1) \\ \hline \mathcal {A} (:, 1, 2) \\ \hline \mathcal {A} (:, 2, 2) \\ \hline \mathcal {A} (:, 3, 2) \end{array} \right] = \left[ \begin{array}{c} a _ {1 1 1} \\ \hline a _ {2 1 1} \\ \hline a _ {1 2 1} \\ \hline a _ {2 2 1} \\ \hline a _ {1 3 1} \\ \hline a _ {2 3 1} \\ \hline a _ {1 1 2} \\ \hline a _ {2 1 2} \\ \hline a _ {1 2 2} \\ \hline a _ {2 2 2} \\ \hline a _ {1 3 2} \\ a _ {2 3 2} \end{array} \right].
$$

Formally, if $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ , then

$$
\operatorname{vec} (\mathcal {A}) = \left[ \begin{array}{c} \operatorname{vec} \left(\mathcal {A} ^ {(1)}\right) \\ \vdots \\ \operatorname{vec} \left(\mathcal {A} ^ {\left(n _ {d}\right)}\right) \end{array} \right] \tag {12.4.2}
$$

where $\mathcal { A } ^ { ( k ) } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d - } }$ 1 is defined by

$$
\mathcal {A} ^ {(k)} (i _ {1}, \dots , i _ {d - 1}) = \mathcal {A} (i _ {1}, \dots , i _ {d - 1}, k) \tag {12.4.3}
$$

for $k = 1 { : } n _ { d }$ . Alternatively, if we define the integer-valued function col by

$$
\mathbf {c o l} (\mathbf {i}, \mathbf {n}) = i _ {1} + (i _ {2} - 1) n _ {1} + (i _ {3} - 1) n _ {1} n _ {2} + \dots + (i _ {d} - 1) n _ {1} \dots n _ {d - 1}, \tag {12.4.4}
$$

then $a = { \mathsf { v e c } } ( { \mathcal { A } } )$ is specified by

$$
a (\operatorname{col} (\mathbf {i}, \mathbf {n})) = \mathcal {A} (\mathbf {i}), \quad \mathbf {1} \leq \mathbf {i} \leq \mathbf {n}. \tag {12.4.5}
$$

# 12.4.4 Tensor Transposition

If $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times n _ { 2 } \times n _ { 3 } }$ , then there are 6 = 3! possible transpositions identified by the notation $\mathcal { A } ^ { < [ i j k ] > }$ where $[ i j k ]$ is a permutation of [1 2 3]:

$$
\mathcal {B} = \left\{ \begin{array}{l} \mathcal {A} ^ {<   [ 1 2 3 ] >} \\ \mathcal {A} ^ {<   [ 1 3 2 ] >} \\ \mathcal {A} ^ {<   [ 2 1 3 ] >} \\ \mathcal {A} ^ {<   [ 2 3 1 ] >} \\ \mathcal {A} ^ {<   [ 3 1 2 ] >} \\ \mathcal {A} ^ {<   [ 3 2 1 ] >} \end{array} \right\} \quad \Longrightarrow \quad \left\{ \begin{array}{l} b _ {i j k} \\ b _ {i k j} \\ b _ {j i k} \\ b _ {j k i} \\ b _ {k i j} \\ b _ {k j i} \end{array} \right\} = a _ {i j k}.
$$

These transpositions can be defined using the perfect shuffle and the vec operator. For example, if $\stackrel { \cdot } { B } = \mathcal { A } ^ { < [ 3 2 1 ] } >$ , then vec $\mathbf { \langle } B \mathbf { \rangle } = ( \mathcal { P } _ { n _ { 1 } , n _ { 2 } } \otimes I _ { n _ { 3 } } ) \mathcal { P } _ { n _ { 1 } n _ { 2 } , n _ { 3 } } \cdot \mathbf { \mathcal { v } } \mathbf { e } \mathbf { c } ( \mathcal { A } )$ .

In general, if $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ and $\mathbf { p } = [ p _ { 1 } , \ldots , p _ { d } ]$ is a permutation of the index vector 1:d, then $\mathcal { A } ^ { < \mathbf { p } > } \in \mathbb { R } ^ { n _ { p _ { 1 } } \times \cdots \times n _ { p _ { d } } }$ is the p-transpose of A defined by

$$
\mathcal {A} ^ {<   \mathbf {p} >} (j _ {p _ {1}}, \ldots , j _ {p _ {d}}) = \mathcal {A} (j _ {1}, \ldots , j _ {d}), \qquad 1 \leq j _ {k} \leq n _ {k}, k = 1: d,
$$

i.e.,

$$
\mathcal {A} ^ {<   \mathbf {p} >} (\mathbf {j} (\mathbf {p})) = \mathcal {A} (\mathbf {j}), \quad \mathbf {1} \leq \mathbf {j} \leq \mathbf {n}.
$$

For additional tensor transposition discussion, see Ragnarsson and Van Loan (2012).

# 12.4.5 The Modal Unfoldings

Recall that a tensor unfolding is a matrix whose entries come from the tensor. Particularly important are the modal unfoldings. If $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ and $N = n _ { 1 } \cdot \cdot \cdot n _ { d }$ , then its mode-k unfolding is an $n _ { k } { \mathrm { - b y } } { \mathrm { - } } ( N / n _ { k } )$ matrix whose columns are the mode-k fibers. To illustrate, here are the three modal unfoldings for A ∈ IR4×3×2: $\mathcal { A } \in \mathbb { R } ^ { 4 \times 3 \times 2 }$

$$
\mathcal {A} _ {(1)} = \left[ \begin{array}{c c c c c c} a _ {1 1 1} & a _ {1 2 1} & a _ {1 3 1} & a _ {1 1 2} & a _ {1 2 2} & a _ {1 3 2} \\ a _ {2 1 1} & a _ {2 2 1} & a _ {2 3 1} & a _ {2 1 2} & a _ {2 2 2} & a _ {2 3 2} \\ a _ {3 1 1} & a _ {3 2 1} & a _ {3 3 1} & a _ {3 1 2} & a _ {3 2 2} & a _ {3 3 2} \\ a _ {4 1 1} & a _ {4 2 1} & a _ {4 3 1} & a _ {4 1 2} & a _ {4 2 2} & a _ {4 3 2} \end{array} \right],
$$

$$
\mathcal {A} _ {(2)} = \left[ \begin{array}{c c c c c c c c} a _ {1 1 1} & a _ {2 1 1} & a _ {3 1 1} & a _ {4 1 1} & a _ {1 1 2} & a _ {2 1 2} & a _ {3 1 2} & a _ {4 1 2} \\ a _ {1 2 1} & a _ {2 2 1} & a _ {3 2 1} & a _ {4 2 1} & a _ {1 2 2} & a _ {2 2 2} & a _ {3 2 2} & a _ {4 2 2} \\ a _ {1 3 1} & a _ {2 3 1} & a _ {3 3 1} & a _ {4 3 1} & a _ {1 3 2} & a _ {2 3 2} & a _ {3 3 2} & a _ {4 3 2} \end{array} \right],
$$

$$
\mathcal {A} _ {(3)} = \left[ \begin{array}{c c c c c c c c c c c c} a _ {1 1 1} & a _ {2 1 1} & a _ {3 1 1} & a _ {4 1 1} & a _ {1 2 1} & a _ {2 2 1} & a _ {3 2 1} & a _ {4 2 1} & a _ {1 3 1} & a _ {2 3 1} & a _ {3 3 1} & a _ {4 3 1} \\ a _ {1 1 2} & a _ {2 1 2} & a _ {3 1 2} & a _ {4 1 2} & a _ {1 2 2} & a _ {2 2 2} & a _ {3 2 2} & a _ {4 2 2} & a _ {1 3 2} & a _ {2 3 2} & a _ {3 3 2} & a _ {4 3 2} \end{array} \right].
$$

We choose to order the fibers left to right according to the “vec” ordering. To be precise, if $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ , then its mode-k unfolding $\boldsymbol { \mathcal { A } } _ { ( \boldsymbol { k } ) }$ is completely defined by

$$
\mathcal {A} _ {(k)} (i _ {k}, \text { col } (\tilde {\mathbf {i}} _ {\mathbf {k}}, \tilde {\mathbf {n}})) = \mathcal {A} (\mathbf {i}) \tag {12.4.6}
$$

where $\tilde { \mathbf { i } } _ { \mathbf { k } } = [ i _ { 1 } , \dots , i _ { k - 1 } , i _ { k + 1 } , \dots , i _ { d } ]$ and $\tilde { \mathbf { n } } _ { \mathbf { k } } = [ n _ { 1 } , \dots , n _ { k - 1 } , n _ { k + 1 } , \dots , n _ { d } ]$ . The rows of $\boldsymbol { \mathcal { A } } _ { ( \boldsymbol { k } ) }$ are associated with subtensors of A. In particular, we can identify $\mathcal { A } _ { ( k ) } ( q , : )$ with the order-(d − 1) tensor $\mathbf { \mathcal { A } } ^ { ( q ) }$ defined by $\mathbf { \nabla } _ { \mathcal { A } } ^ { ( q ) } ( \tilde { \mathbf { i } } _ { k } ) = \mathbf { \nabla } _ { \mathcal { A } _ { ( k ) } } ( q , \mathrm { c o l } ( \tilde { \mathbf { i } } _ { k } ) , \tilde { \mathbf { n } } _ { k } )$ .

# 12.4.6 More General Unfoldings

In general, an unfolding for $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ is defined by choosing a set of row modes and a set of column modes. For example, if $\mathcal { A } \in \mathbb { R } ^ { 2 \times 3 \times 2 \times 2 \times 3 }$ , r = 1:3 and ${ \bf c } = 4 { : } 5$ , then

$$
\mathcal {A} _ {\mathbf {r} \times \mathbf {c}} = \left[ \begin{array}{l l l l l l} a _ {1 1 1, 1 1} & a _ {1 1 1, 2 1} & a _ {1 1 1, 1 2} & a _ {1 1 1, 2 2} & a _ {1 1 1, 1 3} & a _ {1 1 1, 2 3} \\ a _ {2 1 1, 1 1} & a _ {2 1 1, 2 1} & a _ {2 1 1, 1 2} & a _ {2 1 1, 2 2} & a _ {2 1 1, 1 3} & a _ {2 1 1, 2 3} \\ a _ {1 2 1, 1 1} & a _ {1 2 1, 2 1} & a _ {1 2 1, 1 2} & a _ {1 2 1, 2 2} & a _ {1 2 1, 1 3} & a _ {1 2 1, 2 3} \\ a _ {2 2 1, 1 1} & a _ {2 2 1, 2 1} & a _ {2 2 1, 1 2} & a _ {2 2 1, 2 2} & a _ {2 2 1, 1 3} & a _ {2 2 1, 2 3} \\ a _ {1 3 1, 1 1} & a _ {1 3 1, 2 1} & a _ {1 3 1, 1 2} & a _ {1 3 1, 2 2} & a _ {1 3 1, 1 3} & a _ {1 3 1, 2 3} \\ a _ {2 3 1, 1 1} & a _ {2 3 1, 2 1} & a _ {2 3 1, 1 2} & a _ {2 3 1, 2 2} & a _ {2 3 1, 1 3} & a _ {2 3 1, 2 3} \\ a _ {1 1 2, 1 1} & a _ {1 1 2, 2 1} & a _ {1 1 2, 1 2} & a _ {1 1 2, 2 2} & a _ {1 1 2, 1 3} & a _ {1 1 2, 2 3} \\ a _ {2 1 2, 1 1} & a _ {2 1 2, 2 1} & a _ {2 1 2, 1 2} & a _ {2 1 2, 2 2} & a _ {2 1 2, 1 3} & a _ {2 1 2, 2 3} \\ a _ {1 2 2, 1 1} & a _ {1 2 2, 2 1} & a _ {1 2 2, 1 2} & a _ {1 2 2, 2 2} & a _ {1 2 2, 1 3} & a _ {1 2 2, 2 3} \\ a _ {2 2 2, 1 1} & a _ {2 2 2, 2 1} & a _ {2 2 2, 1 2} & a _ {2 2 2, 2 2} & a _ {2 2 2, 1 3} & a _ {2 2 2, 2 3} \\ a _ {1 3 2, 1 1} & a _ {1 3 2, 2 1} & a _ {1 3 2, 1 2} & a _ {1 3 2, 2 2} & a _ {1 3 2, 1 3} & a _ {1 3 2, 2 3} \\ a _ {2 3 2, 1 1} & a _ {2 3 2, 2 1} & a _ {2 3 2, 1 2} & a _ {2 3 2, 2 2} & a _ {2 3 2, 1 3} & a _ {2 3 2, 2 3} \end{array} \right] \begin{array}{l} (1, 1, 1) \\ (2, 1, 1) \\ (1, 2, 1) \\ (2, 2, 1) \\ (1, 3, 1) \\ (2, 3, 1) \\ (1, 1, 2) \\ (2, 1, 2) \\ (1, 2, 2) \\ (2, 2, 2) \\ (1, 3, 2) \\ (2, 3, 2) \end{array} . \tag {12.4.7}
$$

In general, let p be a permutation of 1:d and define the row and column modes by

$$
\mathbf {r} = \mathbf {p} (1: e), \qquad \mathbf {c} = \mathbf {p} (e + 1: d),
$$

where $0 \leq e \leq d .$ . This partitioning defines a matrix $\mathcal { A } _ { \bf r \times c }$ that has $n _ { p _ { 1 } } \cdot \cdot \cdot n _ { p _ { e } }$ rows and $n _ { p _ { e + 1 } } \cdot \cdot \cdot n _ { p _ { d } }$ columns and whose entries are defined by

$$
\mathcal {A} _ {\mathbf {r} \times \mathbf {c}} (\operatorname{col} (\mathbf {i}, \mathbf {n} (\mathbf {r}))  ,   \operatorname{col} (\mathbf {j}, \mathbf {n} (\mathbf {c}))) = \mathcal {A} (\mathbf {i}, \mathbf {j}). \tag {12.4.8}
$$

Important special cases include the modal unfoldings

$$
\mathbf {r} = [ k ], \mathbf {c} = [ 1, \dots , k - 1, k + 1, \dots , d ] \quad \Longrightarrow \quad \mathcal {A} _ {\mathbf {r} \times \mathbf {c}} = \mathcal {A} _ {(k)}
$$

and the vec operation

$$
\mathbf {r} = 1: d, \mathbf {c} = [ \emptyset ] \quad \Longrightarrow \quad \mathcal {A} _ {\mathbf {r} \times \mathbf {c}} = \operatorname{vec} (\mathcal {A}).
$$

# 12.4.7 Outer Products

The outer product of tensor $B \in \mathbb { R } ^ { m _ { 1 } \times \cdots \times m _ { f } }$ with tensor $\mathcal { C } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { g } }$ is the order-$( f + g )$ tensor A defined by

$$
\mathcal {A} (\mathbf {i}, \mathbf {j}) = \mathcal {B} (\mathbf {i}) \circ \mathcal {C} (\mathbf {j}), \qquad \mathbf {1} \leq \mathbf {i} \leq \mathbf {m}, \mathbf {1} \leq \mathbf {j} \leq \mathbf {n}.
$$

Multiple outer products are similarly defined, e.g.,

$$
\mathcal {A} = \mathcal {B} \circ \mathcal {C} \circ \mathcal {D} \quad \implies \quad \mathcal {A} (\mathbf {i}, \mathbf {j}, \mathbf {k}) = \mathcal {B} (\mathbf {i}) \cdot \mathcal {C} (\mathbf {j}) \cdot \mathcal {D} (\mathbf {k}).
$$

Note that if B and C are order-2 tensors (matrices), then

$$
\mathcal {A} = \mathcal {B} \circ \mathcal {C} \quad \Rightarrow \quad \mathcal {A} (i _ {1}, i _ {2}, j _ {1}, j _ {2}) = \mathcal {B} (i _ {1}, i _ {2}) \cdot \mathcal {C} (j _ {1}, j _ {2})
$$

and

$$
\mathcal {A} _ {[ 3 1 ] \times [ 4 2 ]} = B \otimes C.
$$

Thus, the Kronecker product of two matrices corresponds to their outer product as tensors.

# 12.4.8 Rank-1 Tensors

Outer products between order-1 tensors (vectors) are particularly important. We say that $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ is a rank-1 tensor if there exist vectors $z ^ { ( 1 ) } , \dotsc , z ^ { ( d ) } \in \mathbb { R } ^ { n _ { k } }$ such that

$$
\mathcal {A} (\mathbf {i}) = z ^ {(1)} (i _ {1}) \dots z ^ {(d)} (i _ {d}), \quad \mathbf {1} \leq \mathbf {i} \leq \mathbf {n}.
$$

A small example clarifies the definition and reveals a Kronecker product connection:

$$
\mathcal {A} = \left[ \begin{array}{l} u _ {1} \\ u _ {2} \end{array} \right] \circ \left[ \begin{array}{l} v _ {1} \\ v _ {2} \\ v _ {3} \end{array} \right] \circ \left[ \begin{array}{l} w _ {1} \\ w _ {2} \end{array} \right] \quad \Leftrightarrow \quad \left[ \begin{array}{l} a _ {1 1 1} \\ a _ {2 1 1} \\ a _ {1 2 1} \\ a _ {2 2 1} \\ a _ {1 3 1} \\ a _ {2 3 1} \\ a _ {1 1 2} \\ a _ {2 1 2} \\ a _ {1 2 2} \\ a _ {2 2 2} \\ a _ {1 3 2} \\ a _ {2 3 2} \end{array} \right] = \left[ \begin{array}{l} u _ {1} v _ {1} w _ {1} \\ u _ {2} v _ {1} w _ {1} \\ u _ {1} v _ {2} w _ {1} \\ u _ {2} v _ {2} w _ {1} \\ u _ {1} v _ {3} w _ {1} \\ u _ {2} v _ {3} w _ {1} \\ u _ {1} v _ {1} w _ {2} \\ u _ {2} v _ {1} w _ {2} \\ u _ {1} v _ {2} w _ {2} \\ u _ {2} v _ {2} w _ {2} \\ u _ {1} v _ {3} w _ {2} \\ u _ {2} v _ {3} w _ {2} \end{array} \right] = w \otimes v \otimes u.
$$

The modal unfoldings of a rank-1 tensor are highly structured. For the above example we have

$$
\mathcal {A} _ {(1)} = \left[ \begin{array}{l l l l l l} u _ {1} v _ {1} w _ {1} & u _ {1} v _ {2} w _ {1} & u _ {1} v _ {3} w _ {1} & u _ {1} v _ {1} w _ {2} & u _ {1} v _ {2} w _ {2} & u _ {1} v _ {3} w _ {2} \\ u _ {2} v _ {1} w _ {1} & u _ {2} v _ {2} w _ {1} & u _ {2} v _ {3} w _ {1} & u _ {2} v _ {1} w _ {2} & u _ {2} v _ {2} w _ {2} & u _ {2} v _ {3} w _ {2} \end{array} \right] = u \otimes (w \otimes v) ^ {T},
$$

$$
\mathcal {A} _ {(2)} = \left[ \begin{array}{l l l l} u _ {1} v _ {1} w _ {1} & u _ {2} v _ {1} w _ {1} & u _ {1} v _ {1} w _ {2} & u _ {2} v _ {1} w _ {2} \\ u _ {1} v _ {2} w _ {1} & u _ {2} v _ {2} w _ {1} & u _ {1} v _ {2} w _ {2} & u _ {2} v _ {2} w _ {2} \\ u _ {1} v _ {3} w _ {1} & u _ {2} v _ {3} w _ {1} & u _ {1} v _ {3} w _ {2} & u _ {2} v _ {3} w _ {2} \end{array} \right] = v \otimes (w \otimes u) ^ {T},
$$

$$
\mathcal {A} _ {(3)} = \left[ \begin{array}{l l l l l l} u _ {1} v _ {1} w _ {1} & u _ {2} v _ {1} w _ {1} & u _ {1} v _ {2} w _ {1} & u _ {2} v _ {2} w _ {1} & u _ {1} v _ {3} w _ {1} & u _ {2} v _ {3} w _ {1} \\ u _ {1} v _ {1} w _ {2} & u _ {2} v _ {1} w _ {2} & u _ {1} v _ {2} w _ {2} & u _ {2} v _ {2} w _ {2} & u _ {1} v _ {3} w _ {2} & u _ {2} v _ {3} w _ {2} \end{array} \right] = w \otimes (v \otimes u) ^ {T}.
$$

In general, if $z ^ { ( k ) } \in \mathbb { R } ^ { n _ { k } }$ for k = 1:d and

$$
\mathcal {A} = z ^ {(1)} \circ \dots \circ z ^ {(d)} \in \mathbb {R} ^ {n _ {1} \times \dots \times n _ {d}},
$$

then its modal unfoldings are rank-1 matrices:

$$
\mathcal {A} _ {(k)} = z ^ {(k)} \cdot \left(z ^ {(d)} \otimes \dots z ^ {(k + 1)} \otimes z ^ {(k - 1)} \otimes \dots z ^ {(1)}\right) ^ {T}. \tag {12.4.9}
$$


---

<!-- golub_750_799 -->

For general unfoldings of a rank-1 tensor, if p is a permutation of 1:d, ${ \bf r } = { \bf p } ( 1 { : } e )$ , and $\mathbf { c } = \mathbf { p } ( e + 1 { : } d )$ , then

$$
\mathcal {A} _ {\mathbf {r} \times \mathbf {c}} = \left(z ^ {(p _ {e})} \circ \dots \circ z ^ {(p _ {1})}\right) \left(z ^ {(p _ {d})} \circ \dots \circ z ^ {(p _ {e + 1})}\right) ^ {T}. \tag {12.4.10}
$$

Finally, we mention that any tensor can be expressed as a sum of rank-1 tensors

$$
\mathcal {A} \in \mathbb {R} ^ {n _ {1} \times \dots \times n _ {d}} \quad \Longrightarrow \quad \mathcal {A} = \sum_ {\mathbf {i} = \mathbf {1}} ^ {\mathbf {n}} \mathcal {A} (\mathbf {i}) I _ {n _ {1}} (:, i _ {1}) \circ \dots \circ I _ {n _ {d}} (:, i _ {d}).
$$

An important §12.5 theme is to find more informative rank-1 summations than this!

# 12.4.9 Tensor Contractions and Matrix Multiplication

Let us return to the notion of a tensor contraction introduced in §12.4.1. The first order of business is to show that a contraction between two tensors is essentially a matrix multiplication between a pair of suitably chosen unfoldings. This is a useful connection because it facilitates reasoning about high-performance implementation.

Consider the problem of computing

$$
\mathcal {A} (i, j, \alpha_ {3}, \alpha_ {4}, \beta_ {3}, \beta_ {4}, \beta_ {5}) = \sum_ {k = 1} ^ {n _ {2}} \mathcal {B} (i, k, \alpha_ {3}, \alpha_ {4}) \cdot \mathcal {C} (k, j, \beta_ {3}, \beta_ {4}, \beta_ {5}) \tag {12.4.11}
$$

where

$$
\mathcal {A} = \mathcal {A} (1: n _ {1}, 1: m _ {2}, 1: n _ {3}, 1: n _ {4}, 1: m _ {3}, 1: m _ {4}, 1: m _ {5}),
$$

$$
\mathcal {B} = \mathcal {B} (1: n _ {1}, 1: n _ {2}, 1: n _ {3}, 1: n _ {4}),
$$

$$
\mathcal {C} = \mathcal {C} (1: m _ {1}, 1: m _ {2}, 1: m _ {3}, 1: m _ {4}, 1: m _ {5}),
$$

and $n _ { 2 } ~ = ~ m _ { 1 }$ . The index k is a contraction index. The example shows that in a contraction, the order of the output tensor can be (much) larger than the order of either input tensor, a fact that can prompt storage concerns. For example, if $n _ { 1 } =$ $\cdot \cdot \cdot = n _ { 4 } = r$ and $m _ { 1 } = \cdot \cdot \cdot = m _ { 5 } = r$ in (12.4.11), then B and C are $O ( r ^ { 5 } )$ while the output tensor A is $O ( r ^ { 7 } )$ .

The contraction (12.4.11) is a collection of related matrix-matrix multiplications. Indeed, at the slice level we have

$$
\mathcal {A} (:,:, \alpha_ {3}, \alpha_ {4}, \beta_ {3}, \beta_ {4}, \beta_ {5}) = \mathcal {B} (:,:, \alpha_ {3}, \alpha_ {4}) \cdot C (:,:, \beta_ {3}, \beta_ {4}, \beta_ {5}).
$$

Each A-slice is an $n _ { 1 } – \mathrm { b y } – m _ { 2 }$ matrix obtained as a product of an $n _ { \mathrm { 1 } } \mathrm { - } \mathrm { b y } \mathrm { - } n _ { \mathrm { 2 } }$ B-slice and an $m _ { 1 } – \mathrm { b y } – m _ { 2 } \ C – \mathrm { s l i c e }$ .

The summation in a contraction can be over more than just a single mode. To illustrate, assume that

$$
\mathcal {B} = \mathcal {B} (1: m _ {1}, 1: m _ {2}, 1: t _ {1}, 1: t _ {2}),
$$

$$
\mathcal {C} = \mathcal {C} (1: t _ {1}, 1: t _ {2}, 1: n _ {1}, 1: n _ {2}, 1: n _ {3}),
$$

and define $\mathscr { A } = \mathscr { A } ( 1 { : } m _ { 1 } , 1 { : } m _ { 2 } , 1 { : } n _ { 1 } , 1 { : } n _ { 2 } , 1 { : } n _ { 3 } )$ by

$$
\mathcal {A} (i _ {1}, i _ {2}, j _ {1}, j _ {2}, j _ {3}) = \sum_ {k _ {1} = 1} ^ {t _ {1}} \sum_ {k _ {2} = 1} ^ {t _ {2}} \mathcal {B} (i _ {1}, i _ {2}, k _ {1}, k _ {2}) \cdot \mathcal {C} (k _ {1}, k _ {2}, j _ {1}, j _ {2}, j _ {3}). \tag {12.4.12}
$$

Note how “matrix like” this computation becomes with multiindex notation:

$$
\mathcal {A} (\mathbf {i}, \mathbf {j}) = \sum_ {\mathbf {k} = 1} ^ {\mathbf {t}} \mathcal {B} (\mathbf {i}, \mathbf {k}) \cdot \mathcal {C} (\mathbf {k}, \mathbf {j}), \quad \mathbf {1} \leq \mathbf {i} \leq \mathbf {m}, \mathbf {1} \leq \mathbf {j} \leq \mathbf {n}. \tag {12.4.13}
$$

A fringe benefit of this formulation is how nicely it connects to the following matrixmultiplication specification of A:

$$
\mathcal {A} _ {[ 1 2 ] \times [ 3 4 5 ]} = \mathcal {B} _ {[ 1 2 ] \times [ 3 4 ]} \cdot \mathcal {C} _ {[ 1 2 ] \times [ 3 4 5 ]}.
$$

The position of the contraction indices in the example (12.4.12) is convenient from the standpoint of framing the overall operation as a product of two unfoldings. However, it is not necessary to have the contraction indices “on the right” in B and “on the left” in C to formulate the operation as a matrix multiplication. For example, suppose

$$
\mathcal {B} = \mathcal {B} (1: t _ {2}, 1: m _ {1}, 1: t _ {1}, 1: m _ {2}),
$$

$$
\mathcal {C} = \mathcal {C} (1: n _ {2}, 1: t _ {2}, 1: n _ {3}, 1: t _ {1}, 1: n _ {1}),
$$

and that we want to compute the tensor $\begin{array} { r } { A = \ A ( 1 { : } m _ { 1 } , 1 { : } m _ { 2 } , 1 { : } n _ { 1 } , 1 { : } n _ { 2 } , 1 { : } n _ { 3 } ) } \end{array}$ defined by

$$
\mathcal {A} (i _ {2}, j _ {3}, j _ {1}, i _ {1}, j _ {2}) = \sum_ {k _ {1} = 1} ^ {t _ {1}} \sum_ {k _ {2} = 1} ^ {t _ {2}} \mathcal {B} (k _ {2}, i _ {1}, k _ {1}, i _ {2}) \cdot \mathcal {C} (j _ {2}, k _ {2}, j _ {3}, k _ {1}, j _ {1}).
$$

It can be shown that this calculation is equivalent to

$$
\mathcal {A} _ {[ 4 1 ] \times [ 3 5 2 ]} = \mathcal {B} _ {[ 2 4 ] \times [ 3 1 ]} \cdot \mathcal {C} _ {[ 4 2 ] \times [ 5 1 3 ]}.
$$

Hidden behind these formulations are important implementation choices that define the overheads associated with memory access. Are the unfoldings explicitly set up? Are there any particularly good data structures that moderate the cost of data transfer? Etc. Because of their higher dimension, there are typically many more ways to organize a tensor contraction than there are to organize a matrix multiplication.

# 12.4.10 The Modal Product

A very simple but important family of contractions are the modal products. These contractions involve a tensor, a matrix, and a mode. In particular, if $\boldsymbol { S } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ , $M \in \mathbb { R } ^ { m _ { k } \times n _ { k } }$ , and $1 \leq k \leq d .$ , then A is the mode-k product of S and M if

$$
\mathcal {A} _ {(k)} = M \cdot \mathcal {S} _ {(k)}. \tag {12.4.14}
$$

We denote this operation by

$$
\mathcal {A} = \mathcal {S} \times_ {k} M
$$

and remark that

$$
\mathcal {A} \left(\alpha_ {1}, \dots , \alpha_ {k - 1}, i, \alpha_ {k + 1}, \dots , \alpha_ {d}\right) = \sum_ {j = 1} ^ {n _ {k}} M (i, j) \cdot \mathcal {S} \left(\alpha_ {1}, \dots , \alpha_ {k - 1}, j, \alpha_ {k + 1}, \dots , \alpha_ {d}\right)
$$

and

$$
\operatorname{vec} (\mathcal {A}) = \left(I _ {n _ {k + 1} \dots n _ {d}} \otimes M \otimes I _ {n _ {1} \dots n _ {k - 1}}\right) \cdot \operatorname{vec} (\mathcal {S}) \tag {12.4.15}
$$

are equivalent formulations. Every mode-k fiber in S is multiplied by the matrix M.

Using (12.4.15) and elementary facts about the Kronecker product, it is easy to show that

$$
(\mathcal {S} \times_ {k} F) \times_ {j} G = (\mathcal {S} \times_ {j} G) \times_ {k} F, \tag {12.4.16}
$$

$$
(\mathcal {S} \times_ {k} F) \times_ {k} G = \mathcal {S} \times_ {k} (F G), \tag {12.4.17}
$$

assuming that all the dimensions match up.

# 12.4.11 The Multilinear Product

Suppose we are given an order-4 tensor $S \in \mathbb { R } ^ { n _ { 1 } \times n _ { 2 } \times n _ { 3 } \times n _ { 4 } }$ and four matrices

$$
M _ {1} \in \mathbb {R} ^ {m _ {1} \times n _ {1}}, \qquad M _ {2} \in \mathbb {R} ^ {m _ {2} \times n _ {2}}, \qquad M _ {3} \in \mathbb {R} ^ {m _ {3} \times n _ {3}}, \qquad M _ {4} \in \mathbb {R} ^ {m _ {4} \times n _ {4}}.
$$

The computation

$$
\mathcal {A} (\mathbf {i}) = \sum_ {\mathbf {j} = \mathbf {1}} ^ {\mathbf {n}} \mathcal {S} (\mathbf {j}) \cdot M _ {1} (i _ {1}, j _ {1}) \cdot M _ {2} (i _ {2}, j _ {2}) \cdot M _ {3} (i _ {3}, j _ {3}) \cdot M _ {4} (i _ {4}, j _ {4}) \tag {12.4.18}
$$

is equivalent to

$$
\operatorname{vec} (\mathcal {A}) = \left(M _ {4} \otimes M _ {3} \otimes M _ {2} \otimes M _ {1}\right) \operatorname{vec} (\mathcal {S}) \tag {12.4.19}
$$

and is an order-4 example of a multilinear product. As can be seen in the following table, a multilinear product is a sequence of contractions, each being a modal product:

$$
a ^ {(0)} = \operatorname{vec} (S)
$$

$$
a ^ {(1)} = \left(I _ {n _ {4}} \otimes I _ {n _ {3}} \otimes I _ {n _ {2}} \otimes M _ {1}\right) a ^ {(0)}
$$

$$
a ^ {(2)} = \left(I _ {n _ {4}} \otimes I _ {n _ {3}} \otimes M _ {2} \otimes I _ {n _ {1}}\right) a ^ {(1)}
$$

$$
a ^ {(3)} = \left(I _ {n _ {4}} \otimes M _ {3} \otimes I _ {n _ {2}} \otimes I _ {n _ {1}}\right) a ^ {(2)}
$$

$$
a ^ {(4)} = \left(M _ {4} \otimes I _ {n _ {3}} \otimes I _ {n _ {2}} \otimes I _ {n _ {1}}\right) a ^ {(3)}
$$

$$
\operatorname{vec} (\mathcal {A}) = a ^ {(4)}
$$

$$
\mathcal {A} ^ {(0)} = \mathcal {S}
$$

$$
\mathcal {A} _ {(1)} ^ {(1)} = M _ {1} \mathcal {A} _ {(1)} ^ {(0)} \quad (\text { Mode - 1   product })
$$

$$
\mathcal {A} _ {(2)} ^ {(2)} = M _ {2} \mathcal {A} _ {(2)} ^ {(1)} \quad (\text { Mode - 2   product })
$$

$$
\mathcal {A} _ {(3)} ^ {(3)} = M _ {3} \mathcal {A} _ {(3)} ^ {(2)} \quad (\text { Mode - 3   product })
$$

$$
\mathcal {A} _ {(4)} ^ {(4)} = M _ {4} \mathcal {A} _ {(4)} ^ {(3)} \quad (\text { Mode - 4   product })
$$

$$
\mathcal {A} = \mathcal {A} ^ {(4)}
$$

The left column specifies what is going on in Kronecker product terms while the right column displays the four required modal products. The example shows that mode-k operations can be sequenced,

$$
\mathcal {A} = \mathcal {S} \times_ {1} M _ {1} \times_ {2} M _ {2} \times_ {3} M _ {3} \times_ {4} M _ {4},
$$

and that their order is immaterial, e.g.,

$$
\mathcal {A} = \mathcal {S} \times_ {4} M _ {4} \times_ {1} M _ {1} \times_ {2} M _ {2} \times_ {3} M _ {3}.
$$

This follows from (12.4.16).

Because they are used in §12.5, we summarize two key properties of the multilinear product in the following theorem.

Theorem 12.4.1. Suppose $\boldsymbol { S } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ and $M _ { k } \in \mathbb { R } ^ { m _ { k } \times n _ { k } }$ for $k = 1 { : } d .$ . If the tensor $\mathcal { A } \in \mathbb { R } ^ { m _ { 1 } \times \cdots \times m _ { d } }$ is the multilinear product

$$
\mathcal {A} = \mathcal {S} \times_ {1} M _ {1} \times_ {2} M _ {2} \dots \times_ {d} M _ {d},
$$

then

$$
\mathcal {A} _ {(k)} = M _ {k} \cdot \mathcal {S} _ {(k)} \cdot \left(M _ {d} \otimes \dots \otimes M _ {k + 1} \otimes M _ {k - 1} \otimes \dots \otimes M _ {1}\right) ^ {T}.
$$

$I f M _ { 1 } , \ldots , M _ { d }$ are all nonsingular, then $\mathcal { S } = \mathcal { A } \times _ { 1 } M _ { 1 } ^ { - 1 } \times _ { 2 } M _ { 2 } ^ { - 1 } \cdot \cdot \cdot \times _ { d } M _ { d } ^ { - 1 }$ .

Proof. The proof involves equations (12.4.16) and (12.4.17) and the vec ordering of the mode-k fibers in $A _ { ( k ) }$ . □

# 12.4.12 Space versus Time

We close with an example from Baumgartner et al. (2005) that highlights the importance of order of operations and what the space-time trade-off can look like when a sequence of contractions is involved. Suppose that A, B, C and $\mathcal { D }$ are N-by-N-by-Nby-N tensors and that $s$ is defined as follows:

$$
\begin{array}{l} \text { for   } \mathbf {i} = \mathbf {1} _ {4}: \mathbf {N} \\ s = 0 \\ \text { for } \mathbf {k} = \mathbf {1} _ {6}: \mathbf {N} \\ s = s + \mathcal {A} (i _ {1}, k _ {1}, i _ {2}, k _ {2}) \cdot \mathcal {B} (i _ {2}, k _ {3}, k _ {4}, k _ {5}) \cdot \mathcal {C} (k _ {6}, k _ {4}, i _ {4}, k _ {2}) \cdot \mathcal {D} (k _ {1}, k _ {6}, k _ {3}, k _ {5}) \\ \end{array}
$$

$$
\mathbf {e n d}
$$

$$
\mathcal {S} (\mathbf {i}) = s
$$

end

Performed “as is,” this is an ${ \cal O } ( N ^ { 1 0 } )$ calculation. On the other hand, if we can afford an additional pair of $N { \mathrm { - b y } } - N { \mathrm { - b y } } - N { \mathrm { - b y } } - N$ arrays then work is reduced to $O ( N ^ { 6 } )$ . To see this, assume (for clarity) that we have a function $\mathcal { F } = \mathsf { C o n t r a c t 1 } ( \mathcal { G } , \mathcal { H } )$ that computes the contraction

$$
\mathcal {F} (\alpha_ {1}, \alpha_ {2}, \alpha_ {3}, \alpha_ {4}) = \sum_ {\beta_ {1} = 1} ^ {N} \sum_ {\beta_ {2} = 1} ^ {N} \mathcal {G} (\alpha_ {1}, \beta_ {1}, \alpha_ {2}, \beta_ {2}) \cdot \mathcal {H} (\alpha_ {3}, \alpha_ {4}, \beta_ {1}, \beta_ {2}),
$$

a function $\mathcal { F } = \mathsf { C o n t r a c t 2 } ( \mathcal { G } , \mathcal { H } )$ that computes the contraction

$$
\mathcal {F} (\alpha_ {1}, \alpha_ {2}, \alpha_ {3}, \alpha_ {4}) = \sum_ {\beta_ {1} = 1} ^ {N} \sum_ {\beta_ {2} = 1} ^ {N} \mathcal {G} (\alpha_ {1}, \beta_ {1}, \alpha_ {2}, \beta_ {2}) \cdot \mathcal {H} (\beta_ {2}, \beta_ {1}, \alpha_ {3}, \alpha_ {4}),
$$

and a function $\mathcal { F } = \mathsf { C o n t r a c t 3 } ( \mathcal { G } , \mathcal { H } )$ that computes the contraction

$$
\mathcal {F} (\alpha_ {1}, \alpha_ {2}, \alpha_ {3}, \alpha_ {4}) = \sum_ {\beta_ {1} = 1} ^ {N} \sum_ {\beta_ {2} = 1} ^ {N} \mathcal {G} (\alpha_ {2}, \beta_ {1}, \alpha_ {4}, \beta_ {2}) \cdot \mathcal {H} (\alpha_ {1} \beta_ {1}, \alpha_ {3}, \beta_ {2}).
$$

Each of these order-4 contractions requires $O ( N ^ { 6 } )$ flops. By exploiting common subexpressions suggested by the parentheses in

$$
((\mathcal {B} (i _ {2}, k _ {3}, k _ {4}, k _ {5}) \cdot \mathcal {D} (k _ {1}, k _ {6}, k _ {3}, k _ {5})) \cdot \mathcal {C} (k _ {6}, k _ {4}, i _ {4}, k _ {2})) \cdot \mathcal {A} (i _ {1}, k _ {1}, i _ {2}, k _ {2}),
$$

we arrive at the following $O ( N ^ { 6 } )$ specification of the tensor S:

$$
\mathcal {T} _ {1} = \text { Contract1 } (\mathcal {B}, \mathcal {D})
$$

$$
\mathcal {T} _ {2} = \text { Contract2 } (\mathcal {T} _ {1}, \mathcal {C})
$$

$$
\mathcal {S} = \text { Contract3 } (\mathcal {T} _ {2}, \mathcal {A})
$$

Of course, space-time trade-offs frequently arise in matrix computations. However, at the tensor level the stakes are typically higher and the number of options exponential. Systems that are able to chart automatically an optimal course of action subject to constraints that are imposed by the underlying computer system are therefore of interest. See Baumgartner et al. (2005).

# Problems

P12.4.1 Explain why (12.4.1) oversees a block matrix multiplication. Hint. Consider each of the three matrices as n-by-n block matrices with n-by-n blocks.

P12.4.2 Prove that the vec definition (12.4.2) and (12.4.3) is equivalent to the vec definition (12.4.4) and (12.4.5).

P12.4.3 How many fibers are there in the tensor $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d ? } }$ How many slices?

P12.4.5 Prove Theorem 12.4.1.

P12.4.6 Suppose $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ and that $B = A ^ { < } \mathbf { p } >$ where p is a permutation of 1:d. Specify a permutation matrix P so that $B _ { ( k ) } = \mathcal { A } _ { ( p ( k ) ) } P .$ .

P12.4.7 Suppose $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } } , \ N = n _ { 1 } \cdot \cdot \cdot n _ { d }$ , and that p is a permutation of 1:d that involves swapping a single pair of indices, $\mathrm { e . g . , [ 1 4 3 2 5 ] }$ . Determine a permutation matrix $P \in \mathbb { R } ^ { N \times N }$ so that if $B = { \mathcal { A } } ^ { < } \mathbf { p } >$ , then $\mathsf { v e c } ( \boldsymbol { B } ) \ : = \ : P \cdot \mathsf { v e c } ( \boldsymbol { A } )$ .

P12.4.8 Suppose $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ and that $\boldsymbol { \mathcal { A } } _ { ( k ) }$ has unit rank for some k. Does it follow that A is a rank-1 tensor?

P12.4.9 Refer to (12.4.18). Specify an unfolding S of S and an unfolding A of A so that $A =$ $( M _ { 1 } \otimes M _ { 3 } ) S ( M _ { 2 } \otimes M _ { 4 } )$ .

P12.4.10 Suppose $\mathcal { A } \in \mathbb { R } ^ { n _ { 1 } \times \cdots \times n _ { d } }$ and that both p and q are permutations of 1:d. Give a formula for r so that $( \mathcal { A } ^ { < \mathbf { p } > } ) ^ { < \mathbf { q } > } = \mathcal { A } ^ { < \mathbf { r } > }$ .

# Notes and References for §12.4

For an introduction to tensor computations, see:

L. De Lathauwer (1997). “Signal Processing Based on Multilinear Algebra,” PhD Thesis, K.U. Leuven. A. Smilde, R. Bro, and P. Geladi (2004). Multiway Analysis, John Wiley, Chichester, England. T.G. Kolda and B.W. Bader (2009). “Tensor Decompositions and Applications,” SIAM Review 51, 455–500.

For results that connect unfoldings, the vec operation, Kronecker products, contractions, and transposition, see:

S. Ragnarsson and C. Van Loan (2012). “Block Tensor Unfoldings,” SIAM J. Matrix Anal. Applic. 33, 149–169.

Matlab software that supports tensor computations as described in this section includes the Tensor Toolbox:

B.W. Bader and T.G. Kolda (2006). “Algorithm 862: MATLAB Tensor Classes for Fast Algorithm Prototyping,” ACM Trans. Math. Softw., 32, 635–653.

B.W. Bader and T.G. Kolda (2007). “Efficient MATLAB Computations with Sparse and Factored Tensors,” SIAM J. Sci. Comput. 30, 205–231.

The challenges associated with high-performance, large-scale tensor computations are discussed in:

W. Landry (2003). “Implementing a High Performance Tensor Library,” Scientific Programming 11, 273–290.

C. Lechner, D. Alic, and S. Husa (2004). “From Tensor Equations to Numerical Code,” Computer Algebra Tools for Numerical Relativity, Vol. 0411063.

G. Baumgartner, A. Auer, D. Bernholdt, A. Bibireata, V. Choppella, D. Cociorva, X. Gao, R. Harrison, S. Hirata, S. Krishnamoorthy, S. Krishnan, C. Lam, Q. Lu, M. Nooijen, R. Pitzer, J. Ramanujam, P. Sadayappan, and A. Sibiryakov (2005). “Synthesis of High-Performance Parallel Programs for a Class of Ab Initio Quantum Chemistry Models,” Proc. IEEE, 93, 276–292.

The multiway analysis community and the quantum chemistry/electronic structure community each have their own favored style of tensor notation and it is very different! See:

J.L. Synge and A. Schild (1978). Tensor Calculus, Dover Publications, New York. H.A.L. Kiers (2000). “Towards a Standardized Notation and Terminology in Multiway Analysis,” J. Chemometr. 14, 105–122.
