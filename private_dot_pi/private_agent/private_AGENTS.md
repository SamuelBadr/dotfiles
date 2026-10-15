# Global instructions

## Subagents

- You may spawn subagents on your own judgment with any model, except Astra or Fable models (model ID contains `astra` or `fable`): spawn those only when the user asks.
- Omit `model` when launching unless there's a specific reason; per-role models are configured in settings.

## Expensive commands

- Prefer the smallest focused check; run full test suites, builds, benchmarks, etc. only when necessary.
- Don't pipe expensive commands through `head`, `tail`, `grep`, etc. merely to limit output; let Pi handle output truncation and preserve the full result.
- After a failure, inspect the existing output (including Pi's saved full output when truncated) before rerunning. Rerun only when new execution is actually needed.
- To wait for running work (background jobs, CI, Slurm, remote runs, but not subagent runs, which wake you when they finish), poll in a loop: check every 15–30 s, exit as soon as it finishes or fails, and cap the loop at the expected remaining time.

## Compute fleet

Always-on machines, reachable from each other with `ssh <name>` (separate filesystems):

- `kauch` (`kauch01`): 12 cores, 122 GB RAM, GTX 1050 Ti. Shared with other people and a CI runner.
- `frank` (`ftuw`): 16 cores, 60 GB RAM. A desktop with a local user; also runs the user's GitHub Actions runners. Never use more than 8 cores there in total, counting everything of yours already running (the runners included), whatever `fleet` shows as free.
- `hclm` (login node `master.cluster`): Slurm cluster, 8 usable nodes with 96 cores and 1.1 TB RAM each. Partitions `short` (3 h), `medium` (1 day, default), `long` (no limit). `/home` and `/mnt/storage` are shared with the nodes; put big data in `/mnt/storage/$USER`. The login node itself (24 cores, 188 GB RAM) is shared by all cluster users: never use more than 8 cores and 64 GB RAM there in total, counting everything of yours already running. The 64 GB is also enforced for agents there (started via `~/bin/pi`, all combined): past it, processes die, possibly the agent itself.
- The Mac laptop sleeps: never leave long-running work on it.

`fleet` shows free cores and RAM on kauch, frank and the hclm login node, the user's jobs already running there, what is left of the caps on frank and the hclm login node, and idle cores per hclm partition. Place work like this:

- Before starting anything that may use more than 2 cores, 4 GB RAM or 10 minutes, run `fleet`, even to run it where you are, and stay within what it shows as free and left for you.
- Smaller work: where you are.
- Anything bigger: wherever it finishes soonest without queueing (kauch, frank, or a Slurm job on hclm that starts immediately). Split it across them if that helps. CPU is shared fairly between users on the desktops, so idle cores are yours; RAM is not, so check it first.
- Avoid queue time: before submitting, run the same `sbatch` command with `--test-only`, which prints the expected start time. If the job wouldn't start right away, shrink it to fit the idle nodes or use the desktops. Ask the user before submitting anything that would wait in the queue.
- Before submitting any job on hclm, check free space there (`df -h /home /mnt/storage /mnt/scratch`) and make sure the job's output fits with room to spare.
- In Slurm jobs, always pass `-c`, `--mem`, `--time` and `--hint=nomultithread`; on hclm a job without `--mem` claims the whole node's memory.
- Never use hyperthreads: at most one thread per physical core, BLAS threads included (set `OPENBLAS_NUM_THREADS` alongside `julia -t N`). Give every parallel command an explicit count, never `auto`.
- Before running on a machine, check its software environment: toolchains up to date and at the same version as where the code was developed (e.g. `juliaup update`, `julia --version`), the code at the same commit, and the project environment instantiated from the same lockfile (e.g. `Pkg.instantiate()`). Fix what you can; tell the user what you can't.
- Work runs where its files are: move code with git (push here, pull there) or `rsync`, and fetch results with `rsync`.
- Start remote runs detached (`nohup … > log 2>&1 &`) and report host, path, and PID or job ID.
- On `kauch`, herdr has `hclm` and `frank` saved: `herdr --machine <name> …` reaches the agents there.

### VSC: ask first

The user also has accounts on the Austrian Scientific Computing clusters: VSC-5 (`ssh VSC5B` as `samuelb`, `ssh VSC5` as `samuel`; Slurm, 128 cores and 512 GB–2 TB RAM per node, A40/A100 GPU nodes, 72 h job limit), VSC-4 (`VSC4B`/`VSC4`) and MUSICA (`ssh MUSICA`).

- Use them only when the desktops and hclm can't run the work without waiting (GPUs, more than hclm can give, or a long queue on hclm). Jobs are billed to the user's project.
- Logging in needs the user's password and 2FA, so you can't connect on your own. When you need access, ask the user to run `ssh -fN VSC5B` (or the host you need) in a terminal on this machine; `ssh -O check VSC5B` tells you when the connection is up.
