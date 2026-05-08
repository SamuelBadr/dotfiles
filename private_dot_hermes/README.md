# Hermes Agent Configuration

This directory contains the Hermes Agent configuration managed by chezmoi.

## Files

- `dot_config.yaml` - Main Hermes configuration (API keys templated)
- `dot_env.tmpl` - Environment variables template (secrets managed by chezmoi)

## Secrets

The following secrets are required and managed by chezmoi's keyring:

| Secret Name | Description | Example |
|-------------|-------------|---------|
| `hermes_api_key` | Custom LLM API key | `***REMOVED***...jWRU` |
| `firecrawl_api_key` | Firecrawl API for web search | `***REMOVED***...4f81` |
| `telegram_bot_token` | Telegram bot token | `***REMOVED***:***` |
| `sudo_password` | SUDO password for local commands | `***REMOVED***` |

## Setup on New Machine

### Option 1: Run the setup script (recommended)

```bash
chezmoi apply --init
./run_setup-hermes-secrets.sh
```

### Option 2: Manual setup with keyring

```bash
chezmoi secret keyring set --service=chezmoi --user=hermes_api_key
chezmoi secret keyring set --service=chezmoi --user=firecrawl_api_key
chezmoi secret keyring set --service=chezmoi --user=telegram_bot_token
chezmoi secret keyring set --service=chezmoi --user=sudo_password
chezmoi apply
```

### Option 3: Use .chezmoidata.toml (no keyring)

Create `~/.local/share/chezmoi/.chezmoidata.toml` (NOT committed to git):

```toml
hermes_api_key = "sk-your-key-here"
firecrawl_api_key = "fc-your-key-here"
telegram_bot_token = "***REMOVED***:your-token"
sudo_password = "your-password"
```

Then run:
```bash
chezmoi apply
```

## Configuration Highlights

- **Model**: Custom endpoint at `https://aqueduct.ai.datalab.tuwien.ac.at/v1`
- **Web Search**: Firecrawl backend
- **Terminal**: Local backend with 60s timeout
- **Reasoning Effort**: High
- **Browser**: Auto engine with 2min inactivity timeout

## Adding More Secrets

To add optional API keys (Exa, Parallel, FAL, etc.):

1. Uncomment the line in `dot_env.tmpl`
2. Add the secret: `chezmoi secret keyring set --service=chezmoi --user=exa_api_key`
3. Run `chezmoi apply`

## Troubleshooting

**Keyring not available:**
- Use `.chezmoidata.toml` instead (Option 3 above)
- Or install a keyring: `sudo apt install gnome-keyring`

**Secrets not being substituted:**
- Run `chezmoi doctor` to check keyring access
- Verify secrets exist: `chezmoi secret keyring get --service=chezmoi --user=firecrawl_api_key`

**Config changes not applying:**
- Run `chezmoi diff` to see what would change
- Run `chezmoi -n apply` for dry-run
