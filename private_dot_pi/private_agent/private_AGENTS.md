# Global instructions

## Suggestions for next step

- End every reply with a numbered list of options for how to proceed, ranked by your recommendation (best first), so I can pick one by replying with just its number.
- Each item is an alternative to the others, complete on its own, not a step in a plan: picking 2 must not assume 1 was done.

## Math

- Write all math (no exceptions) as LaTeX (`$…$` inline, `$$…$$` display).

## Research wiki

- `~/.llm-wiki` is my research wiki (pi-llm-wiki; the same git repo on every machine, synced automatically). Before you create or edit wiki pages, read `~/.llm-wiki/WIKI_SCHEMA.md` and follow it: page types and format, math, links, and which pages you may change.
- As you work, record what is worth keeping (results and dead ends, derivations, settled conventions), as the schema describes.
- Never call `wiki_bootstrap`.

## Subagents

- Spawn subagents only when the user explicitly asks for them. Never on your own judgment, however large or risky the task: do the work yourself. This includes workflows, chains and scheduled runs.
- Omit `model` when launching unless there's a specific reason; per-role models are configured in settings.

## Readable tool calls

I read your tool calls as they run, so write them for a human reader:

- Prefer the dedicated tools (read, edit, write) over `cat`, `sed -i`, heredocs and the like.
- One step per line, not long `;`/`&&` chains. Clear variable names, long option names when the short ones are obscure, and a short `#` comment for any step whose purpose isn't obvious.
- No nested quoting or escaping puzzles. Anything longer than a few lines, or code in another language (`python -c`, `julia -e`), goes into a script file with a one-line header comment, which you then run.
- Label multi-part output (`echo "== tests =="`) so each part is easy to find.

## Expensive commands

- Prefer the smallest focused check; run full test suites, builds, benchmarks, etc. only when necessary.
- Don't pipe expensive commands through `head`, `tail`, `grep`, etc. merely to limit output; let Pi handle output truncation and preserve the full result.
- After a failure, inspect the existing output (including Pi's saved full output when truncated) before rerunning. Rerun only when new execution is actually needed.
- To wait for running work (background jobs, CI, Slurm, remote runs, but not subagent runs, which wake you when they finish), use one bash call whose shell loop checks every 15–30 s (Slurm queues: every few minutes) and exits as soon as the work finishes or fails, e.g. `timeout 1800 sh -c "while kill -0 $PID 2>/dev/null; do sleep 20; done"; tail log`. Set the `timeout` to the expected remaining time and the bash tool timeout a bit above it. A loop costs one turn however long it waits, so never use `sleep N; check` with N over 30 to save turns.

## Compute

- Before starting anything that may use more than 2 cores, 4 GB RAM or 10 minutes, run `fleet`, even to run it where you are. It shows the machines you can use, how to reach them, their rules and what is usable now. Stay within that: limits are for my user in total, shared by all agents and jobs.
- Smaller work: run it where you are.
- Bigger work: wherever it finishes soonest, counting setup, data transfer and queue wait; split it across machines if that helps. `sbatch --test-only` (same command) prints the expected start; shrinking a job to fit idle nodes often beats waiting.
- Slurm jobs: always pass `-c`, `--mem`, `--time` and `--hint=nomultithread`, and make sure the output fits with room to spare.
- Never use hyperthreads: at most one thread per physical core you use, counting Julia threads times BLAS threads (usually `OPENBLAS_NUM_THREADS=1` when Julia threads do the parallel work). Give every parallel command an explicit count, never `auto`.
- Before running on another machine: same toolchain version, same commit, environment instantiated from the same lockfile. Fix what you can; tell me what you can't. Copy code with `rsync` (git only for commits already pushed); fetch results with `rsync`.
- Start remote runs detached (`nohup … > log 2>&1 &`) and report host, path, and PID or job ID.
- Benchmarks: on a desktop or login node, `fleet hold <host> <hours> <reason>` first (rerun to extend), wait until `fleet` shows it quiet, and `fleet release <host>` when done. On Slurm, use `--exclusive`.
- Never attempt a login that may need my password or 2FA: failed logins can get machines blocked. Use only connections `fleet` shows as connected; for missing ones, give me the command it prints.
