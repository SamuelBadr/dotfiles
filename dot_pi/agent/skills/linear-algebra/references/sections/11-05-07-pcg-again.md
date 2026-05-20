# 11.5.7 PCG—Again

The polynomial preconditioner discussion points to an important connection between the classical iterations and the preconditioned conjugate gradient algorithm. Many iterative methods have as their basic step

$$
x _ {k} = x _ {k - 2} + \omega_ {k} \left(\gamma_ {k - 1} z _ {k - 1} + x _ {k - 1} - x _ {k - 2}\right) \tag {11.5.4}
$$

where $M z _ { k - 1 } = r _ { k - 1 } = b - A x _ { k - 1 }$ . For example, if we set $\omega _ { k } = 1$ and $\gamma _ { k } = 1$ , then

$$
x _ {k} = M ^ {- 1} (b - A x _ {k - 1}) + x _ {k - 1},
$$

i.e., $M x _ { k } = N x _ { k - 1 } + b$ , where $A \ = \ M - N$ . Following Concus, Golub, and O’Leary (1976), it is also possible to organize the preconditioned CG method with a central step of the form (11.5.4):

$$
x _ {- 1} = 0; k = 0; r _ {0} = b - A x _ {0}
$$

while $r _ { k } \neq 0$

$$
k = k + 1
$$

$\mathrm { S o l v e ~ } M z _ { k - 1 } = r _ { k - 1 } \mathrm { ~ f o r ~ } z _ { k - 1 }$

$$
\gamma_ {k - 1} = z _ {k - 1} ^ {T} M z _ {k - 1} / z _ {k - 1} ^ {T} A z _ {k - 1}
$$

$\mathbf { i f } \ k = 1$

$$
\omega_ {1} = 1
$$

else

$$
\omega_ {k} = \left(1 - \frac {\gamma_ {k - 1}}{\gamma_ {k - 2}} \frac {z _ {k - 1} ^ {T} M z _ {k - 1}}{z _ {k - 2} ^ {T} M z _ {k - 2}} \frac {1}{\omega_ {k - 1}}\right) ^ {- 1}
$$

end

$$
x _ {k} = x _ {k - 2} + \omega_ {k} \big (\gamma_ {k - 1} z _ {k - 1} + x _ {k - 1} - x _ {k - 2} \big)
$$

$$
r _ {k} = b - A x _ {k}
$$

end

$$
x = x _ {k}
$$

Thus, we can think of the scalars $\omega _ { k }$ and $\gamma _ { k }$ in this iteration as acceleration parameters that can be chosen to speed the convergence of the iteration $M x _ { k } = N x _ { k - 1 } + b$ . Hence, any iterative method based on the splitting $A = M - N$ can be accelerated by the conjugate gradient algorithm as long as M (the preconditioner) is symmetric and positive definite.
