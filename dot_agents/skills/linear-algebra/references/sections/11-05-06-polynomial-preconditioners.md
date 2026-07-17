# 11.5.6 Polynomial Preconditioners

Suppose $A = M _ { 1 } - N _ { 1 }$ is a splitting and that $\rho ( G ) < 1$ where $G = M _ { 1 } ^ { - 1 } N _ { 1 }$ . Since $A = M _ { 1 } ( I - G )$ , it follows that

$$
A ^ {- 1} = (I - G) ^ {- 1} M _ {1} ^ {- 1} = \left(\sum_ {k = 0} ^ {\infty} G ^ {k}\right) M _ {1} ^ {- 1}.
$$

This suggests another way to generate a preconditioner whose inverse resembles the inverse of A. We simply truncate the infinite series:

$$
M ^ {- 1} = \left(\sum_ {k = 0} ^ {m} G ^ {k}\right) M _ {1} ^ {- 1}.
$$

It follows that

$$
z = \left(I + G + G ^ {2} + \dots G ^ {m}\right) M _ {1} ^ {- 1} r
$$

solves $M z = r$ . Moreover, there is a very simple way to compute this vector:

$$
z _ {c} = 0
$$

for k = 1:m

$$
M _ {1} z _ {+} = N _ {1} z _ {c} + r
$$

$$
z _ {c} = z _ {+}
$$

end

$$
z = z _ {c}
$$

To see why this works, we note that $z _ { + } = G z _ { c } + d$ where $M _ { 1 } d = r$ , and apply induction:

$$
z _ {+} = G z _ {c} + d = G (I + G + \dots + G ^ {k - 1}) d + d = (I + G + \dots G ^ {k}) d.
$$

Thus, the $M z = r$ calculation requires m steps of the iteration $M _ { 1 } z _ { + } = N _ { 1 } z _ { c } + r$

In the polynomial preconditioner paradigm, the given system $A x = b$ is replaced by $M ^ { - 1 } A x = M ^ { - 1 } b$ where the preconditioner M is defined by

$$
M ^ {- 1} = p (M _ {1} ^ {- 1} A) M _ {1} ^ {- 1}. \tag {11.5.3}
$$

Here, p is a polynomial and $M _ { 1 }$ is itself a preconditioner, e.g., the diagonal of A. In the above example, p was determined by the parameter m and the chosen $M _ { 1 }$ .

We mention that there are more sophisticated ways to design a good polynomial preconditioner. With $M _ { 1 } = I$ for clarity in (11.5.3), the goal is for $p ( A )$ to look like $A ^ { - 1 }$ , i.e., we want $I \approx p ( A ) A$ . Note that $I - p ( A ) A = q ( A )$ where $q ( z ) = 1 - z p ( z )$ , so the challenge is to find $q \in \mathbb { P } _ { m + 1 }$ with the property that $q ( 0 ) = 1$ and $q ( A )$ is small. There are several ways to address this optimization problem in practice, see Ashby, Manteuffel, and Otto (1992) and Saad(1985).
