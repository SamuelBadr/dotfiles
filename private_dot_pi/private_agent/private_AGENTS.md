# Shared login node — computing rules

- This is a shared login node, not a compute node.
- Computation on the login node is fine as long as it doesn't hog the node.
- Use Slurm for anything plausibly exceeding ~8 CPU cores, ~64 GB RAM, or
  ~2 hours of compute.
- Always ask before submitting Slurm jobs, unless explicitly instructed to
  submit them.
- Do not disable or bypass the resource safeguards on this host.
