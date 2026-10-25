# Global instructions

## Next steps

- End every reply with a numbered list of the next steps I'm most likely to want, best first, so I can pick one by replying with just its number.

## Subagents

- Spawn subagents only when the user explicitly asks for them. Never on your own judgment, however large or risky the task: do the work yourself. This includes workflows, chains and scheduled runs.
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
- `hclm` (login node `master.cluster`): Slurm cluster, 8 usable nodes with 96 cores and 1.1 TB RAM each. Partitions `short` (3 h), `medium` (1 day, default), `long` (no limit); never use `cqs`. `/home` and `/mnt/storage` are shared with the nodes; put big data in `/mnt/storage/$USER`. The login node itself (24 cores, 188 GB RAM) is shared by all cluster users: never use more than 8 cores and 64 GB RAM there in total, counting everything of yours already running. The 64 GB is also enforced for agents there (started via `~/bin/pi`, all combined): past it, processes die, possibly the agent itself.
- The Mac laptop sleeps: never leave long-running work on it.

`fleet` shows free cores and RAM on kauch, frank and the hclm login node, the user's jobs already running there, what is left of the caps on frank and the hclm login node, idle cores per hclm partition, and idle nodes, your jobs and connection state on VSC and MUSICA (below). Place work like this:

- Before starting anything that may use more than 2 cores, 4 GB RAM or 10 minutes, run `fleet`, even to run it where you are, and stay within what it shows as free and left for you on the shared machines.
- Smaller work: where you are.
- A machine `fleet` shows as HELD is reserved for someone's quiet work (benchmarks): start nothing there, not even small work, unless you are its holder (the `host:dir` it names). To benchmark on kauch, frank or the hclm login node, first `fleet hold <host> <hours> <reason>` (refused if someone else holds it), wait until `fleet` shows it quiet, and `fleet release <host>` as soon as you are done; pick hours that cover the run, and re-run `hold` to extend. On Slurm, benchmark with `--exclusive` instead.
- Anything bigger: wherever it finishes soonest, queue wait included (kauch, frank, or a Slurm job on hclm, VSC or MUSICA). Joining a queue is fine: a job never has to start immediately; submit into a busy partition when that still gives the earliest finish. Split it across them if that helps. CPU is shared fairly between users on the desktops, so idle cores are yours; RAM is not, so check it first.
- Estimate queue wait: run the same `sbatch` command with `--test-only`, which prints the expected start time. Slurm assumes every running job uses its full time limit, so trust that time on VSC, MUSICA and hclm's `short` and `medium`. In hclm's `long`, many running jobs have no real limit (365 days or unlimited), so a start that isn't immediate often shows years away: treat it as an unknown wait, and prefer anywhere with a known finish. Shrinking a job to fit idle nodes often beats waiting.
- Before submitting any job on hclm, check free space there (`df -h /home /mnt/storage /mnt/scratch`) and make sure the job's output fits with room to spare.
- In Slurm jobs, always pass `-c`, `--mem`, `--time` and `--hint=nomultithread`; on hclm a job without `--mem` claims the whole node's memory.
- Never use hyperthreads: at most one thread per physical core, BLAS threads included (set `OPENBLAS_NUM_THREADS` alongside `julia -t N`). Give every parallel command an explicit count, never `auto`.
- Before running on a machine, check its software environment: toolchains up to date and at the same version as where the code was developed (e.g. `juliaup update`, `julia --version`), the code at the same commit, and the project environment instantiated from the same lockfile (e.g. `Pkg.instantiate()`). Fix what you can; tell the user what you can't.
- Work runs where its files are: move code with git (push here, pull there) or `rsync`, and fetch results with `rsync`.
- Start remote runs detached (`nohup … > log 2>&1 &`) and report host, path, and PID or job ID.
- On `kauch`, herdr has `hclm` and `frank` saved: `herdr --machine <name> …` reaches the agents there.

### VSC and MUSICA

The user also has accounts on the Austrian Scientific Computing (ASC) Slurm clusters, part of the fleet like hclm. Jobs are billed to the user's projects; that's fine, use whatever finishes first.

- VSC-5 (`VSC5B` as `samuelb`, project p72144; `VSC5` as `samuel`, project p71998): 128 cores and 512 GB–2 TB RAM per node; A40 (16 cores) and A100 (128 cores) nodes with 2 GPUs each; 72 h job limit.
- VSC-4 (`VSC4B`, `VSC4`, same accounts): 48 cores and 96–768 GB RAM per node. VSC-4 and VSC-5 share `/home` per account; big data goes in `/gpfs/data/fs72144/samuelb`.
- MUSICA (`MUSICA` Vienna, `MUSICA-LNZ` Linz, `MUSICA-INN` Innsbruck; user `sb69466`; a separate Slurm cluster per site): 192-core nodes with 768 GB, 4×H100 nodes, and in Vienna 8×B200 nodes with 2.3 TB. Slurm is on `PATH` only in login shells: `ssh MUSICA bash -lc '…'`.
- Login needs the user's 2FA, so connections exist only as SSH masters the user opened on hclm; they last for days. On hclm run `ssh VSC5B …` directly; elsewhere run `ssh hclm ssh VSC5B …`, or work from an agent on hclm.
- Never try to log in yourself: failed logins can get hclm blocked. Use only hosts `fleet` shows as connected, or check with `ssh -O check <host>` first. `fleet` (also run every 10 minutes by cron on hclm) reconnects dropped hosts whenever it can without 2FA. When a cluster you need is still unreachable, tell the user and give them the exact command `fleet` prints for it, to run in a terminal on hclm.
- Their login nodes are shared like hclm's: compile and copy files there, but run real work as Slurm jobs.
