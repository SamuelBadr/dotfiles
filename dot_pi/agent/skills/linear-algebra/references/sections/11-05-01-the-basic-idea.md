# 11.5.1 The Basic Idea

Suppose $M = M _ { 1 } M _ { 2 }$ is nonsingular and consider the linear system $\tilde { A } \tilde { x } = \tilde { b }$ where

$$
\tilde {A} = M _ {1} ^ {- 1} A M _ {2} ^ {- 1}, \qquad \tilde {b} = M _ {1} ^ {- 1} b.
$$

Note that if M looks like A, then $\tilde { A }$ looks like I. The proposal is to solve the “tilde problem” with a suitably chosen Krylov procedure and then determine x by solving $M _ { 2 } x \ = \ \tilde { x }$ . The matrix M is called a preconditioner and it must have two attributes for this solution framework to be of interest:

Criterion 1. M must capture the essence of A, for if $M \approx A$ , then we have $I \approx$ $M _ { 1 } ^ { - 1 } A M _ { 2 } ^ { - 1 } = \tilde { A }$ . (In settings where M is specified through its inverse, it is more appropriate to say that $M ^ { - 1 }$ captures the essence of $A ^ { - \bar { 1 } } . )$

Criterion 2. It must be easy to solve linear systems that involve the matrices $M _ { 1 }$ and $M _ { 2 }$ because the Krylov process involves the operation $( M _ { 1 } ^ { - 1 } A M _ { 2 } ^ { - 1 } )$ -times-vector.

Having a good preconditioner means fewer iterations. However, the cost of an iteration is an issue, as is the overhead associated with the construction of $M _ { 1 }$ and $M _ { 2 }$ . Thus, the enthusiasm for a preconditioner depends upon the strength of the inequality

$$
\left( \begin{array}{c} \text {Set up} \\ M \\ \text {cost} \end{array} \right) + \left( \begin{array}{c} \text {Single} \\ \tilde {A} \text {-iteration} \\ \text {cost} \end{array} \right) \cdot \left( \begin{array}{c} \text {Number} \\ \text {of} \tilde {A} \\ \text {iterations} \end{array} \right) <   \left( \begin{array}{c} \text {Single} \\ A \text {-iteration} \\ \text {cost} \end{array} \right) \cdot \left( \begin{array}{c} \text {Number} \\ \text {of} A \\ \text {iterations} \end{array} \right).
$$

There are several ways in which a preconditioner M can capture the essence of A. The difference $A - M$ could be small in norm or low in rank. More generally, if

$$
A = [ \text { friendly / important   part } ] + [ \text { troublesome / lesser   part } ],
$$

then the important part is an obvious candidate for a preconditioner subject to the constraint imposed by Criterion 2. For example, if A is symmetric positive definite, then its diagonal qualifies as an important part that is computationally friendly.
