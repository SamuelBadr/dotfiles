---
name: julia
description: Julia 1.12 development guidance for implementing, reviewing, optimizing, testing, documenting, and packaging Julia code. Use when Codex needs modern Julia best practices for multiple dispatch, type stability, memory/performance work, environments, package setup, CI, docs, release automation, or package extensions.
---

# Julia

Assume Julia 1.12 unless the project declares tighter `julia` compat or an older CI matrix.

## Start Here

- Read `Project.toml`, `src/`, `test/`, `docs/`, and CI files before editing.
- Activate the project environment before running code: `julia --project=.` and `Pkg.instantiate()`.
- Match existing package structure and style instead of imposing a new layout.
- Keep `SKILL.md` advice high-level; load `references/performance.md` or `references/workflows.md` only when the task needs them.

## Write Julia Idiomatically

- Prefer small methods and multiple dispatch over large type switches.
- Prefer generic APIs over over-constrained signatures; add annotations when they encode a real interface or improve inference.
- Prefer immutable structs with concrete field types unless mutation is semantically necessary.
- Keep hot code inside functions and avoid untyped globals.
- Preserve type stability across branches and return values; use `zero`, `oneunit`, `oftype`, or explicit conversion when needed.
- Use function barriers when setup logic is dynamic but the inner kernel should specialize.
- Iterate with `eachindex`, `pairs`, or direct array traversal instead of assuming specific indices or storage.
- Respect Julia naming conventions already present in the codebase. Avoid style churn.

## Optimize Deliberately

- Measure before changing code. Use `@time` only for a first pass; use BenchmarkTools or Chairmarks for comparisons.
- Use `JET.@report_opt` and `@code_warntype` to investigate inference and type-instability issues.
- Preallocate and favor in-place `!` methods in repeated hot loops when allocations matter.
- Use broadcasting and dotted assignment deliberately; prefer clarity unless fusion or in-place execution is important.
- Use `@views` when repeated slice copies are wasteful, but copy irregular data if repeated contiguous work is faster.
- Keep dense-array inner loops aligned with column-major storage.
- Delay `@inbounds`, `@simd`, `@fastmath`, threading, or sysimage work until profiling shows they matter.

## Maintain Packages Cleanly

- Prefer the standard package layout: `Project.toml`, `src/`, `test/`, optionally `docs/`.
- Use project-local environments for code and keep personal tooling in the default environment.
- Prefer package extensions for optional dependencies when compat allows them.
- Add or update tests with each behavior change.
- Respect `.JuliaFormatter.toml` when present. If absent, avoid repository-wide formatting churn.
- Keep docstrings close to public APIs and integrate with Documenter when the package already uses it or needs hosted docs.

## Load References On Demand

- Read [`references/performance.md`](references/performance.md) for type stability, allocation control, benchmarking, profiling, latency, and threading guidance.
- Read [`references/workflows.md`](references/workflows.md) for environments, Revise-based development, package scaffolding, CI, formatting, Aqua, JET, Documenter, CompatHelper, TagBot, and extensions.
