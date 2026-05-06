# Known errata

This file records local corrections or cautions for *Matrix Computations* as used by this skill.

## Algorithm 5.1.1 — Householder Vector

In the branch

- `σ = 0`
- `x(1) < 0`

the book gives

- `β = -2`

but the Householder reflector consistent with the stated contract and orthogonality requires

- `β = 2`

This matters for a 1-vector input like `x = [-5]`, where `β = -2` would produce a non-orthogonal update.

When implementing or reasoning from Algorithm 5.1.1, use `β = 2` in that branch.
