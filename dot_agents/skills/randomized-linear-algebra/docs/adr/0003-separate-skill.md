# Separate skill, not merged into linear-algebra

Randomized NLA is a distinct subfield (sketching, randomized SVD/rangefinder, trace estimation, Nyström) with its own methods and analysis tools (concentration inequalities, random matrix theory). It lives in a separate skill from `linear-algebra` for three reasons: (1) the addressing scheme differs (sequential + named labels vs. 3-part numeric), (2) the retrieval structure differs (single-document survey vs. 12-chapter textbook), and (3) the trigger boundary is clean (randomized vs. deterministic), so a merged skill would need a branching workflow.

Rejected: one skill with two references. Forcing both into one SKILL.md means the workflow says "if G&V, do X; if M&T, do Y" — a branching skill is a badly-factored skill. The two skills cross-link via bidirectional "see also" pointers instead.

See also: `linear-algebra` skill ADR-0001 for the shared retrieval architecture.
