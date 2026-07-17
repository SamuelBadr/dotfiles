# Addressing scheme: sequential number + named label (not 3-part renumbering)

The survey numbers theorems/lemmas/corollaries sequentially with a shared counter (so they jump: Theorem 6, 8, 9, 13…), and cross-references them two ways: by sequential number in the prose ("see Theorem 8") and by named label in the LaTeX (`thm:rand-power`). Algorithms and equations use named labels only.

We support both keys: `theorem 8` (sequential) and `theorem rand-power` (label) resolve to the same block. Equations resolve by label only (`eq rand-power-nogap`); the `(5.1.2)` paren form from the `linear-algebra` skill does not apply.

Rejected: auto-renumber everything into a fake `chapter.section.number` scheme to mirror the sibling skill. That would *lie* — printed numbers in the book wouldn't match, so every body-text cross-reference ("see Theorem 8") would point at the wrong thing. Also rejected: labels-only (loses "Theorem 8" prose resolution) and numbers-only (loses `thm:rand-power` label resolution).
