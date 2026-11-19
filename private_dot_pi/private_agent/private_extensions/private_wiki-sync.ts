// Keeps the research wiki (~/.llm-wiki) in sync across machines: runs ~/bin/wiki-sync
// when a session starts, after every settled turn and at shutdown (wiki-sync only goes
// to GitHub when there is something to push or the last pull is a few minutes old, so
// these calls are cheap). When the wiki changed, it refreshes the qmd search index
// that the research-wiki skill searches.
import { execFile } from "node:child_process";
import { homedir } from "node:os";
import { join } from "node:path";
import { promisify } from "node:util";
import type { ExtensionAPI, ExtensionContext } from "@earendil-works/pi-coding-agent";

const run = promisify(execFile);
const repo = join(homedir(), ".llm-wiki");
const syncScript = join(homedir(), "bin", "wiki-sync");

// Syncs outlive turns and sessions, so they use node's execFile rather than pi.exec
// (which throws once the session is gone) and never let an error escape.
function warn(ctx: ExtensionContext, message: string) {
  try {
    ctx.ui.notify(message, "warning");
  } catch {
    // ctx is stale after shutdown; nowhere left to report
  }
}

function lastLine(error: unknown): string {
  const { stderr, stdout, message } = error as { stderr?: string; stdout?: string; message: string };
  return (stderr || stdout || message).trim().split("\n").at(-1) ?? "";
}

async function head(): Promise<string> {
  try {
    return (await run("git", ["-C", repo, "rev-parse", "HEAD"])).stdout.trim();
  } catch {
    return "";
  }
}

async function sync(ctx: ExtensionContext, when: "start" | "turn" | "shutdown") {
  const before = await head();
  if (!before) return; // wiki not cloned on this machine
  try {
    await run(syncScript, [], { timeout: 90_000 });
  } catch (error) {
    warn(ctx, `wiki-sync failed: ${lastLine(error)} (run ~/bin/wiki-sync for details)`);
  }
  // no reindex at shutdown, so pi exits promptly; the next session start catches up
  if (when === "shutdown") return;
  // wiki-sync commits local edits, so an unchanged HEAD means nothing to reindex;
  // at session start reindex anyway: the Mac's launchd job may have pulled since
  if (when === "turn" && (await head()) === before) return;
  try {
    await run("qmd", ["update"], { timeout: 60_000 });
    await run("qmd", ["embed"], { timeout: 600_000 });
  } catch (error) {
    // qmd not installed on this machine: search falls back to rg, nothing to reindex
    if ((error as { code?: string }).code === "ENOENT") return;
    warn(ctx, `wiki-sync: qmd reindex failed: ${lastLine(error)}`);
  }
}

export default function (pi: ExtensionAPI) {
  pi.on("session_start", (_event, ctx) => {
    void sync(ctx, "start");
  });
  pi.on("agent_settled", (_event, ctx) => {
    void sync(ctx, "turn");
  });
  pi.on("session_shutdown", async (_event, ctx) => {
    await sync(ctx, "shutdown");
  });
}
