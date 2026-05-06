# Pi Setup Sync Guide

This document explains how to sync your pi coding agent setup to a new machine using chezmoi.

## What's Synced

The following pi configuration is now tracked in your chezmoi dotfiles repository:

### Core Configuration
- `~/.pi/agent/settings.json` - Packages, default models, theme preferences
- `~/.pi/agent/mcp.json` - MCP server configuration
- `~/.pi/agent/models.json.tmpl` - Custom model providers (API key stored in keyring)

### Customizations
- `~/.pi/agent/agents/` - Custom agent definitions
- `~/.pi/agent/extensions/` - Custom extensions
- `~/.pi/agent/skills/linear-algebra/` - Custom linear algebra skill
- `~/.agents/skills/` - All custom skills (touying-author, typst, julia, obsidian, etc.)

### Setup Script
- `run_onchange_install-pi-packages.sh` - Automatically installs pi packages after apply

## What's NOT Synced

These are intentionally excluded (cache/regenerable/sensitive):
- `~/.pi/cache/` - regenerable cache
- `~/.pi/agent/sessions/` - conversation history
- `~/.pi/agent/mcp-cache.json` - regenerable
- `~/.pi/agent/auth.json` - re-authenticate on new machine
- `~/.pi/agent/run-history.jsonl` - local history

## Setup on New Machine

### 1. Apply chezmoi
```bash
chezmoi apply
```

### 2. Retrieve the API Key Secret

The Aqueduct API key is stored in your system keyring. On the new machine, you need to add it:

**Option A: If you have keyring sync (recommended)**
- macOS: Ensure iCloud Keychain is enabled
- Linux: Use the same keyring backend (gnome-keyring, kwallet, etc.)

**Option B: Manual setup**
```bash
chezmoi secret keyring set --service=chezmoi-aqueduct --user=aqueduct_api_key --value="YOUR_API_KEY"
```

### 3. Re-authenticate with pi
```bash
# Pi will prompt for authentication on first run
pi
```

### 4. Verify packages are installed
The setup script should have automatically installed packages from settings.json. Verify:
```bash
pi list
```

If packages weren't auto-installed, manually install:
```bash
pi install npm:pi-markdown-preview
pi install git:https://github.com/badlogic/pi-diff-review
pi install https://github.com/davebcn87/pi-autoresearch
pi install npm:pi-subagents
pi install npm:@plannotator/pi-extension
pi install npm:pi-mcp-adapter
```

## Updating Configuration

When you modify pi settings on your main machine:

1. Changes to `settings.json`, `mcp.json`, or skills are auto-tracked
2. Commit and push happens automatically (configured in chezmoi.toml)
3. On other machines: `chezmoi apply`

## Troubleshooting

### Models.json not rendering correctly
Ensure the secret is set:
```bash
chezmoi secret keyring get --service=chezmoi-aqueduct --user=aqueduct_api_key
```

### Packages not installing
Check if jq is installed (required by the setup script):
```bash
brew install jq  # macOS
apt install jq   # Debian/Ubuntu
```

### Skills not loading
Verify symlinks are created correctly in `~/.pi/agent/skills/`. Some skills reference `~/.agents/skills/`.
