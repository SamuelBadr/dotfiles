# 7.7.3 Sensitivity Issues

The generalized Schur decomposition sheds light on the issue of eigenvalue sensitivity for the $A - \lambda B$ problem. Clearly, small changes in A and B can induce large changes in the eigenvalue $\lambda _ { i } = t _ { i i } / s _ { i i }$ if $s _ { i i }$ is small. However, as Stewart (1978) argues, it may not be appropriate to regard such an eigenvalue as “ill-conditioned.” The reason is that the reciprocal $\mu _ { i } = s _ { i i } / t _ { i i }$ might be a very well-behaved eigenvalue for the pencil $\mu A - B$ . In the Stewart analysis, A and B are treated symmetrically and the eigenvalues are regarded more as ordered pairs $( t _ { i i } , s _ { i i } )$ than as quotients. With this point of view it becomes appropriate to measure eigenvalue perturbations in the chordal metric chord(a, b) defined by

$$
\operatorname{chord} (a, b) = \frac {| a - b |}{\sqrt {1 + a ^ {2}} \sqrt {1 + b ^ {2}}}.
$$

Stewart shows that if λ is a distinct eigenvalue of A − λB and $\lambda _ { \epsilon }$ is the corresponding eigenvalue of the perturbed pencil $\tilde { A } - \lambda \tilde { B }$ with $\parallel A - \tilde { A } \parallel _ { 2 } \approx \parallel B - \tilde { B } \parallel _ { 2 } \approx \epsilon ,$ , then

$$
\operatorname{chord} (\lambda , \lambda_ {\epsilon}) \leq \frac {\epsilon}{\sqrt {(y ^ {H} A x) ^ {2} + (y ^ {H} B x) ^ {2}}} + O (\epsilon^ {2})
$$

where x and y have unit 2-norm and satisfy Ax = λBx and $y ^ { H } A { = } \lambda y ^ { H } B$ . Note that the denominator in the upper bound is symmetric in A and B. The “truly” ill-conditioned eigenvalues are those for which this denominator is small.

The extreme case when both $t _ { k k }$ and $s _ { k k }$ are zero for some k has been studied by Wilkinson (1979). In this case, the remaining quotients $t _ { i i } / s _ { i i }$ can take on arbitrary values.
