# Chezmoi Modules

This repo uses a single module to gate workstation-only configuration (GUI tools
and Homebrew packages) so the same dotfiles repo can be applied on both
workstations and headless cluster/login nodes.

## Modules

| Module | Description | Default |
|--------|-------------|---------|
| `workstation` | GUI apps + Homebrew packages (nvim, bat, btop, starship, Brewfile) | `false` |

Everything else (shell, git, tmux, SSH, pi agent + skills) is always managed.

## Configuration

### Repo Defaults (`.chezmoidata.toml`)

Safe defaults for all machines are in `.chezmoidata.toml` at the repo root.

### Machine-Local Overrides (`~/.config/chezmoi/chezmoi.toml`)

Enable the workstation module per-machine:

```toml
[data.modules]
workstation = true
```

## Example Configurations

### Cluster/Login Node (Minimal)

No local config needed — `workstation = false` is the default.

### Workstation (Full)

```toml
# ~/.config/chezmoi/chezmoi.toml
[data.modules]
workstation = true
```

## Secrets

Secrets should NOT be stored in this repo. Use `chezmoi add --encrypt` with age
for any sensitive values, or keep them machine-local in
`~/.config/chezmoi/chezmoi.toml` (which is not managed by chezmoi).

## Setup Commands

### On a New Machine

1. Create local config (workstations only):
   ```bash
   mkdir -p ~/.config/chezmoi
   cat > ~/.config/chezmoi/chezmoi.toml <<'EOF'
   [data.modules]
   workstation = true
   EOF
   ```

2. Apply dotfiles:
   ```bash
   chezmoi apply
   ```

### Testing Before Apply

```bash
# Check what modules are active
chezmoi execute-template '{{ .modules }}'

# Preview changes
chezmoi diff

# Dry run
chezmoi apply --dry-run
```

## Adding a New Module

1. Add default to `.chezmoidata.toml`:
   ```toml
   [modules]
     new_module = false
   ```

2. Gate files/scripts with:
   ```
   {{- if .modules.new_module -}}
   # file contents
   {{- end -}}
   ```

3. Update this documentation.
