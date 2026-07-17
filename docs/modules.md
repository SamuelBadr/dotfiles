# Chezmoi Modules

This repo uses a module system to enable/disable optional components per machine.

## Modules

| Module | Description | Default |
|--------|-------------|---------|
| `pi_agent` | Pi coding agent configuration | `true` |
| `gui_desktop` | GUI apps: nvim, bat, btop, starship, fzf | `false` |
| `packages` | Brewfile and package installation | `false` |

## Configuration

### Repo Defaults (`.chezmoidata.toml`)

Safe defaults for all machines are in `.chezmoidata.toml` at the repo root.

### Machine-Local Overrides (`~/.config/chezmoi/chezmoi.toml`)

Enable modules per-machine in your local chezmoi config:

```toml
[data.modules]
pi_agent = true
gui_desktop = true
packages = true
```

## Example Configurations

### Cluster/Login Node (Minimal)

```toml
# ~/.config/chezmoi/chezmoi.toml
[data.modules]
pi_agent = true
gui_desktop = false
packages = false
```

### Workstation (Full)

```toml
# ~/.config/chezmoi/chezmoi.toml
[data.modules]
pi_agent = true
gui_desktop = true
packages = true
```

## Secrets

Secrets should NOT be stored in this repo. Use `chezmoi add --encrypt` with age/Boring
Crypto for any sensitive values, or keep them machine-local in `~/.config/chezmoi/chezmoi.toml`
(which is not managed by chezmoi).

## Setup Commands

### On a New Machine

1. Create local config:
   ```bash
   mkdir -p ~/.config/chezmoi
   nano ~/.config/chezmoi/chezmoi.toml
   ```

2. Apply dotfiles:
   ```bash
   chezmoi apply --verbose
   ```

### Testing Before Apply

```bash
# Check what modules are active
chezmoi execute-template '{{ .modules }}'

# Preview changes
chezmoi diff --verbose

# Dry run
chezmoi apply --dry-run --verbose
```

## Adding a New Module

1. Add default to `.chezmoidata.toml`:
   ```toml
   [modules]
     new_module = false
   ```

2. Gate files/scripts with:
   ```typst
   {{- if .modules.new_module -}}
   # file contents
   {{- end -}}
   ```

3. Gate scripts with:
   ```bash
   {{- if not .modules.new_module -}}
   #!/usr/bin/env bash
   exit 0
   {{- else -}}
   #!/usr/bin/env bash
   set -euo pipefail
   # script body
   {{- end -}}
   ```

4. Update this documentation.

## Troubleshooting

### Missing Key Errors

If you see `map has no entry for key "xxx"`, add the key to `.chezmoidata.toml` or enable the module that provides it.

### Disabled Module Still Prompts

Ensure the script starts with the proper conditional and `exit 0` for disabled state.
