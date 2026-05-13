# Julia Performance Notes

Use this file when the task is about performance, allocations, type instability, latency, or threading.

## Baseline First

- Put performance-critical work inside functions.
- Avoid untyped global variables in hot paths.
- Use `@time` for a quick first signal, then switch to BenchmarkTools or Chairmarks for real comparisons.
- Interpolate benchmark inputs with `$` so benchmark results are not distorted by global-variable lookups.

## Type Stability

- Keep return types stable across branches whenever practical.
- Prefer concrete field types in structs used on hot paths.
- Avoid containers with abstract element or field types in performance-sensitive code.
- Use function barriers when dynamic setup feeds a stable inner kernel.
- Use `@code_warntype` for manual inspection and `JET.@report_opt` for broader automated analysis.

## Memory and Arrays

- Preallocate outputs when the same computation runs repeatedly.
- Prefer dotted calls and `.=` when they remove temporary arrays without obscuring intent.
- Use `@views` or `view` when repeated slicing would otherwise allocate unnecessarily.
- Do not assume views are always faster. If repeated work on irregular access dominates, copying into contiguous storage can win.
- Traverse dense arrays in column-major order: make the first index vary fastest in inner loops.

## Profiling

- Use `Profile` for runtime hotspots and `Profile.Allocs` for memory hotspots.
- Prefer flame-graph tools such as VS Code profiling, ProfileView, or PProf when a call tree matters more than a single timing number.
- Optimize only after a benchmark or profile shows where time or allocations actually go.

## Latency and Compilation

- Separate startup latency from steady-state runtime when reporting performance.
- Reach for precompile work or PackageCompiler sysimages only if latency is a demonstrated problem.
- Treat package load time, first-call compilation, and steady-state execution as different costs.

## Parallel and BLAS Notes

- Watch for oversubscription when Julia threads call multithreaded BLAS kernels.
- When using multithreaded Julia code with BLAS-heavy kernels, `OPENBLAS_NUM_THREADS=1` is often a good baseline to test.
- Validate thread settings on the actual workload instead of assuming more threads help.

## Sources

- Julia manual performance tips: <https://docs.julialang.org/en/v1/manual/performance-tips/>
- Modern Julia Workflows optimizing guide: <https://modernjuliaworkflows.org/optimizing/>
