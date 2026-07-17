# 12.2.6 The Givens-Vector Representation

The QR factorization of a semiseparable matrix is also an $O ( n )$ computation. To motivate the algorithm we step through a simple special case that showcases the idea of a structured rank Givens update. Along the way we will discover yet another strategy that can be used to represent a semiseparable matrix.

Assume $A _ { L } \in \mathbb { R } ^ { n \times n }$ is a lower triangular semiseparable matrix and that $a \in \mathbb { R } ^ { n }$ is its first column. We can reduce this column to a multiple of $e _ { 1 }$ with a sequence of

n − 1 Givens rotations, $\mathrm { e . g . }$

$$
\left[ \begin{array}{c c c c} c _ {1} & s _ {1} & 0 & 0 \\ - s _ {1} & c _ {1} & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{array} \right] \left[ \begin{array}{c c c c} 1 & 0 & 0 & 0 \\ 0 & c _ {2} & s _ {2} & 0 \\ 0 & - s _ {2} & c _ {2} & 0 \\ 0 & 0 & 0 & 1 \end{array} \right] \left[ \begin{array}{c c c c} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & c _ {3} & s _ {3} \\ 0 & 0 & - s _ {3} & c _ {3} \end{array} \right] \left[ \begin{array}{c} a _ {1} \\ a _ {2} \\ a _ {3} \\ a _ {4} \end{array} \right] = \left[ \begin{array}{c} v _ {1} \\ 0 \\ 0 \\ 0 \end{array} \right].
$$

By moving the rotations to the right-hand side we see that

$$
A _ {L} (:, 1) = \left[ \begin{array}{c} a _ {1} \\ a _ {2} \\ a _ {3} \\ a _ {4} \end{array} \right] = v _ {1} \left[ \begin{array}{c} c _ {1} \\ c _ {2} s _ {1} \\ c _ {3} s _ {2} s _ {1} \\ s _ {3} s _ {2} s _ {1} \end{array} \right].
$$

Because this is the first column of a semiseparable matrix, it is not hard to show that there exist “weights” $v _ { 2 } , \ldots , v _ { n }$ so that

$$
A _ {L} = \left[ \begin{array}{c c c c} c _ {1} v _ {1} & 0 & 0 & 0 \\ c _ {2} s _ {1} v _ {1} & c _ {2} v _ {2} & 0 & 0 \\ c _ {3} s _ {2} s _ {1} v _ {1} & c _ {3} s _ {2} v _ {2} & c _ {3} v _ {3} & 0 \\ s _ {3} s _ {2} s _ {1} v _ {1} & s _ {3} s _ {2} v _ {2} & s _ {3} v _ {3} & v _ {4} \end{array} \right] = B (s) ^ {- T}. * \operatorname{tril} \left(c v ^ {T}\right) \tag {12.2.11}
$$

where

$$
v = \left[ \begin{array}{l} v _ {1} \\ v _ {2} \\ v _ {3} \\ v _ {4} \end{array} \right], \qquad c = \left[ \begin{array}{l} c _ {1} \\ c _ {2} \\ c _ {3} \\ 1 \end{array} \right], \qquad s = \left[ \begin{array}{l} s _ {1} \\ s _ {2} \\ s _ {3} \end{array} \right].
$$

The encoding (12.2.11) is an example of the Givens-vector representation for a triangular semiseparable matrix. It consists of a vector of cosines, a vector of sines, and a vector of weights. By “transposing” this idea, we can similarly represent an upper triangular semiseparable matrix. Thus, for a general semiseparable matrix A we may write

$$
A = A _ {L} + A _ {U},
$$

where

$$
A _ {L} = \operatorname{tril} (A) = B \left(s _ {L}\right) ^ {- T}. * \operatorname{tril} \left(c _ {L} v _ {L} ^ {T}\right),
$$

$$
A _ {U} = \operatorname{triu} (A, 1) = B \left(s _ {U}\right) ^ {- 1}. * \operatorname{triu} \left(v _ {U} c _ {U} ^ {T}, 1\right),
$$

where $c _ { L } , \ s _ { L } .$ , and $v _ { L }$ (resp. $c _ { U } , ~ s _ { U }$ , and $v _ { U } )$ are the cosine, sine, and weight vectors associated with the lower (resp. upper) triangular part. For more details on the properties and utility of this representation, see Vandebril and Van Barel (2005).
