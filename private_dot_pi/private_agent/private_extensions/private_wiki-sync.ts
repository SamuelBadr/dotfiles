// Keeps the research wiki (~/.llm-wiki, maintained by pi-llm-wiki) in sync across
// machines: runs ~/bin/wiki-sync when a session starts, after every settled turn
// and at shutdown. pi-llm-wiki searches a generated registry that only local
// writes refresh, so after a pull brings pages from another machine this
// rebuilds it, otherwise recall would miss them.
import { execFile } from "node:child_process";
import { homedir } from "node:os";
import { join } from "node:path";
import { promisify } from "node:util";
import type { ExtensionAPI, ExtensionContext } from "@earendil-works/pi-coding-agent";

const run = promisify(execFile);
const home = homedir();
const repo = join(home, ".llm-wiki");
const syncScript = join(home, "bin", "wiki-sync");
// shortcut: calls pi-llm-wiki internals (version pinned in settings.json); recheck after upgrading it
const llmWikiLib = join(home, ".pi/agent/npm/node_modules/@zosmaai/pi-llm-wiki/extensions/llm-wiki/lib");

// Syncs outlive turns and sessions, so they use node's execFile rather than pi.exec
// (which throws once the session is gone) and never let an error escape.
function warn(ctx: ExtensionContext, message: string) {
  try {
    ctx.ui.notify(message, "warning");
  } catch {
    // ctx is stale after shutdown; nowhere left to report
  }
}

async function head(): Promise<string> {
  try {
    return (await run("git", ["-C", repo, "rev-parse", "HEAD"])).stdout.trim();
  } catch {
    return "";
  }
}

async function rebuildIndex() {
  const { rebuildMetadata } = await import(join(llmWikiLib, "metadata.ts"));
  const { getVaultPaths } = await import(join(llmWikiLib, "utils.ts"));
  rebuildMetadata(getVaultPaths(home));
}

async function sync(ctx: ExtensionContext, alwaysRebuild = false) {
  const before = await head();
  if (!before) return; // wiki not cloned on this machine
  try {
    await run(syncScript, [], { timeout: 90_000 });
  } catch (error) {
    const { stderr, stdout, message } = error as { stderr?: string; stdout?: string; message: string };
    const reason = (stderr || stdout || message).trim().split("\n").at(-1);
    warn(ctx, `wiki-sync failed: ${reason} (run ~/bin/wiki-sync for details)`);
  }
  // at session start, rebuild even without new commits: launchd may have pulled since the last session
  if (!alwaysRebuild && (await head()) === before) return;
  try {
    await rebuildIndex();
  } catch (error) {
    warn(ctx, `wiki-sync: could not rebuild the wiki index: ${(error as Error).message}`);
  }
}

export default function (pi: ExtensionAPI) {
  pi.on("session_start", (_event, ctx) => {
    void sync(ctx, true);
  });
  pi.on("agent_settled", (_event, ctx) => {
    void sync(ctx);
  });
  pi.on("session_shutdown", async (_event, ctx) => {
    await sync(ctx);
  });
}
