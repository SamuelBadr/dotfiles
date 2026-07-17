# linear-algebra skill

A retrieval layer that grounds numerical-linear-algebra implementation in *Matrix Computations* by Golub and Van Loan. The skill exists to insure accuracy: the agent references the canonical book rather than relying on memory, which is often wrong for numerical-linalg formulas and algorithms.

## Language

**Grounding**:
Grounding is the means by which the skill achieves accuracy. The implementation's algorithm choice is identified by book section/label and justified against the problem's conditions (dimensions, symmetry, definiteness, conditioning) before code is written.
_Avoid_: validation, traceability, citation-for-its-own-sake

**Canonical content**:
The specific labeled block that is the authoritative source for a reference — a theorem statement, an algorithm's pseudocode, a tagged equation — as opposed to the surrounding section prose.
_Avoid_: section text, surrounding discussion, "the passage about X"

**Lookup**:
A two-hop retrieval act: (1) find the right section or labeled item, (2) surface its canonical content. The accuracy job is done at hop 2; hop 1 is only navigation.
