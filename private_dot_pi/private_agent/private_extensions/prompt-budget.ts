// Trims fixed per-session prompt cost without losing features:
// - Rarely used wiki tools stay active and callable from codemode scripts (tools.wiki_lint({})),
//   but their declarations and rules leave the prompt. pi-llm-wiki registers them "direct", and
//   deactivating would make them uncallable, so an inert tool hides them via prepareLoadout,
//   the mechanism codemode itself uses (docs/extensions.md, Tool Exposure).
// - pi-emphasis's ~1.6k-token reply rules are off; a short version lives in ~/.pi/agent/AGENTS.md.
//   Set here rather than in common.sh, which pi launched outside an interactive shell misses.
import type { ExtensionAPI } from "@earendil-works/pi-coding-agent";

const SELF = "prompt_budget";
const SCRIPT_ONLY = [
	"wiki_capture_source", "wiki_ingest", "wiki_lint", "wiki_rebuild_meta", "wiki_status",
	"wiki_watch", "wiki_reindex", "wiki_reindex_embeddings", "wiki_log_event", "wiki_bootstrap",
];

process.env.PI_EMPHASIS_GUIDANCE ??= "0";

export default function promptBudget(pi: ExtensionAPI) {
	pi.registerTool({
		name: SELF,
		label: "prompt budget",
		description: "Internal; never call.",
		parameters: { type: "object", properties: {} } as any,
		exposure: "model-only",
		prepareLoadout: () => ({ hiddenDeclarations: [SELF, ...SCRIPT_ONLY] }),
		async execute() {
			return { content: [{ type: "text", text: "Internal tool; nothing to do." }], details: undefined };
		},
	});
}
