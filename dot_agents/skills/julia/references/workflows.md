# Julia Workflow Notes

Use this file when the task involves package setup, environments, testing, CI, documentation, or release workflow.

## Installation and Environments

- Prefer `juliaup` for installing and switching Julia versions.
- Activate a local environment for each project instead of working from the global default environment.
- Treat `Project.toml` as the source of direct dependencies and compat bounds.
- Treat `Manifest.toml` as the full resolved dependency graph for a concrete environment.
- Keep personal development tools such as Revise, BenchmarkTools, Aqua, JET, JuliaFormatter, or PackageCompiler in the default environment when that matches the user's setup.

## Development Loop

- Use Revise to shorten the edit-run loop during interactive development.
- Start Julia with the project environment active before testing package behavior.
- Prefer working from the package environment rather than ad hoc `include` patterns when editing packages.

## Package Scaffolding

- Prefer PkgTemplates for new packages instead of hand-assembling CI, docs, and release automation.
- Reuse the repository's existing template choices when contributing to an established package.
- Keep the standard Julia package shape unless the project already uses a deliberate alternative.

## Testing and Quality

- Put behavior changes behind tests in `test/`.
- Use GitHub Actions or the repository's existing CI provider for package checks.
- Use JuliaFormatter for consistent formatting when the repository already opted into it.
- Use Aqua to catch package hygiene issues such as ambiguous methods or dependency problems.
- Use JET as a complementary static and optimization-analysis tool.
- Use TestItemRunner when the project already uses test items or editor-driven test execution.

## Documentation and Release Automation

- Keep public APIs documented with docstrings close to the code.
- Use Documenter for package documentation when hosted docs are needed.
- Use CompatHelper to keep dependency compat bounds current when the project already uses GitHub automation.
- Use TagBot to automate release tagging when the package is registry-backed and already follows that workflow.

## Optional Dependencies

- Prefer package extensions over legacy conditional-loading patterns when optional functionality depends on extra packages.
- Reach for extension helper packages only if the repository already uses them or the setup becomes repetitive.

## Sources

- Modern Julia Workflows writing guide: <https://modernjuliaworkflows.org/writing/>
- Modern Julia Workflows sharing guide: <https://modernjuliaworkflows.org/sharing/>
