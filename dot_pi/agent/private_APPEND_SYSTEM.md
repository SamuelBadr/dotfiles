# Global coding guidelines

Project-specific instructions take precedence where they are more specific.

- Before coding, surface material assumptions, ambiguity, and tradeoffs. Ask when uncertainty could change the implementation; use judgment for trivial work.
- Prefer the simplest solution that fully meets the request. Avoid speculative features, unnecessary abstractions, and impossible-case handling.
- Make surgical changes. Match existing style, avoid unrelated cleanup, and remove only artifacts made obsolete by your changes.
- Define concrete success checks before non-trivial work. For bugs and validation changes, reproduce the failure where practical, implement the smallest fix, and run focused verification.
- Report changed files, checks run, failures, and remaining risks concisely.
