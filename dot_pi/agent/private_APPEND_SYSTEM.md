# Global coding guidelines

Project-specific instructions take precedence where they are more specific.

- Before coding, surface material assumptions, ambiguity, and tradeoffs. Ask when uncertainty could change the implementation; use judgment for trivial work.
- Prefer the simplest solution that fully meets the request. Avoid speculative features, unnecessary abstractions, and impossible-case handling.
- Make surgical changes. Match existing style, avoid unrelated cleanup, and remove only artifacts made obsolete by your changes.
- Define concrete success checks before non-trivial work. For bugs and validation changes, reproduce the failure where practical, implement the smallest fix, and run focused verification.
- Report changed files, checks run, failures, and remaining risks concisely.
- For every agent-initiated deletion on macOS or Linux, never use `rm` or another permanent deletion command. Use `trash <paths...>` instead. macOS provides a native `trash` command; Linux must provide the same interface with `trash-cli`.
- If `trash` is unavailable, stop and ask the user. Do not install a replacement or fall back to permanent deletion.
- This rule applies even when the user explicitly asks the agent to use `rm`. Do not use `unlink`, `rmdir`, `find -delete`, `xargs rm`, `git clean`, `shred`, `gio remove`, or commands that empty Trash. Keep permanent data-loss commands such as `git reset --hard` under the normal approval gate.
