# Chezmoi Modules

Workstation configuration is optional so this repo works on Macs and headless
cluster/login nodes. Backup configuration requires a separate, explicit Mac opt-in.

## Modules

| Module | Description | Default |
|--------|-------------|---------|
| `workstation` | Workstation tool configuration, one-time macOS keyboard preferences | `false` |
| `backupMac` | This Mac's restic configuration and LaunchAgents (macOS only) | `false` |

Shell, Git, tmux, SSH, and Pi configuration are always managed. All skills
are installed and updated by `skills`; chezmoi does not manage `~/.agents/skills`
or its lock metadata.

## Homebrew packages

`Brewfile` in the repo root is the package manifest. It is repository-only
(`.chezmoiignore`d, never deployed to `~`); `HOMEBREW_BUNDLE_FILE` in
`~/.config/shell/common.sh` points Homebrew at it, so the native commands edit
the tracked file directly:

```bash
brew bundle add <pkg>       # install-and-record; use --cask for casks
brew bundle remove <pkg>    # drop the entry
brew bundle install         # install everything listed
brew bundle cleanup         # list installed packages missing from the manifest
```

The manifest is deliberately not templated, so Homebrew can parse it. It is
shared by all machines: the backup formulae (`restic`, `runitor`, `coreutils`)
are listed unconditionally even though only the `backupMac` machine runs them.

## Tool ownership

- **Shell shortcuts:** shared, shell-neutral setup lives in `~/.config/shell/common.sh`
  (sourced by both zsh and the bash fallback). By design it defines only `up` and
  `pluto`; treat additions as deliberate. zsh-specific config stays in `~/.zshrc`.
- **Editor:** Zed with `--wait`; nano fallback when Zed is unavailable.
- **Node:** Homebrew `node` and npm. NVM is a separate per-user Node version
  manager; no NVM installation or version file exists here, so its old shell
  initialization was dead code and has been removed.
- **Rust:** Homebrew `rustup` supplies the manager; rustup owns the stable
  toolchain in `~/.rustup` and the active Rust binaries in `~/.cargo/bin`.
  Do not install the Homebrew `rust` formula alongside it.

## Pi settings ownership

`.chezmoitemplates/pi-settings.json` is the shared source of truth for Pi
configuration (packages, subagent config, theme, tuiMode). A chezmoi
modify-template deploys it to `~/.pi/agent/settings.json` while preserving
machine-local runtime keys (`defaultProvider`, `defaultModel`,
`defaultThinkingLevel`, `modelThinkingLevels`, `enabledModels`,
`lastChangelogVersion`).

Package installs stay in sync automatically: `~/.config/shell/common.sh` wraps
`pi` so `pi install` / `pi remove` copy the live `packages` list back into the
manifest and commit+push it; every machine inherits the change on its next
`chezmoi update`. `pi update --all` only reconciles checkouts, never the list,
so no sync is needed there. Other stable config (theme, subagent models) is
edited in the repo manifest.

Edit the stable JSON in the source repo. Do not `chezmoi add
~/.pi/agent/settings.json`; that would recapture machine-local runtime choices.

## Mac backups

The current Mac opts in with `[data.modules] backupMac = true` in the untracked
local chezmoi config. `[data.backupPings]` supplies the `backup`, `check`, and
`prune` Healthchecks ping IDs. These are capability-bearing secrets; do not commit
them. Restic's repository password remains in Keychain, and SSH keys remain local.

Only the source list, exclusions, scripts, tests, README and three LaunchAgents
are managed. Logs, archives, repository caches, credentials and launchd runtime
state are not tracked. `chezmoi apply` deploys files but does not bootstrap/reload
jobs or start backups. On a replacement Mac, restore the Keychain/SSH credentials,
install the Brewfile dependencies, review paths, then load backup/check jobs manually.
Pruning is deliberately disabled in its plist until scope/restore/integrity checks
pass; see `~/.config/restic-hclm/README.md` before reenabling it.

## Skills

Install and update skills separately with `npx skills add ... -g` and
`npx skills update -g`. `chezmoi apply` does not install or update them.

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

# Offline regression checks (from the source repo)
python3 tests/test_dotfiles.py
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
