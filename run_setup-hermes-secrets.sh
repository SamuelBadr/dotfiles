#!/bin/bash
# Setup Hermes Agent secrets for chezmoi
# Run this on a new machine before `chezmoi apply`

set -e

echo "=== Hermes Agent Secrets Setup ==="
echo ""
echo "This script will help you set up the required secrets for Hermes Agent."
echo "Secrets are stored in your system keyring (not in git)."
echo ""

# Check if keyring is available
if ! chezmoi secret keyring get --service=chezmoi --user=test 2>&1 | grep -q "not found\|no entry"; then
    echo "✓ System keyring is available"
else
    echo "⚠ System keyring not available. Consider using a .chezmoidata.toml file instead."
    echo "  Create ~/.local/share/chezmoi/.chezmoidata.toml with:"
    echo "    firecrawl_api_key = \"your-key\""
    echo "    telegram_bot_token = \"your-token\""
    echo "    sudo_password = \"your-password\""
    echo "    hermes_api_key = \"your-api-key\""
    echo ""
    read -p "Continue with keyring anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo ""
echo "Enter your secrets (values will be hidden):"
echo ""

# Firecrawl API Key
read -sp "Firecrawl API Key (fc-...): " FIRECRAWL_KEY
echo ""
echo -n "$FIRECRAWL_KEY" | chezmoi secret keyring set --service=chezmoi --user=firecrawl_api_key --value="$(cat)"
echo "✓ Firecrawl API key stored"

# Telegram Bot Token
read -sp "Telegram Bot Token: " TELEGRAM_TOKEN
echo ""
echo -n "$TELEGRAM_TOKEN" | chezmoi secret keyring set --service=chezmoi --user=telegram_bot_token --value="$(cat)"
echo "✓ Telegram bot token stored"

# SUDO Password
read -sp "SUDO Password: " SUDO_PASS
echo ""
echo -n "$SUDO_PASS" | chezmoi secret keyring set --service=chezmoi --user=sudo_password --value="$(cat)"
echo "✓ SUDO password stored"

# Hermes API Key (custom LLM)
read -sp "Hermes API Key (sk-...): " HERMES_KEY
echo ""
echo -n "$HERMES_KEY" | chezmoi secret keyring set --service=chezmoi --user=hermes_api_key --value="$(cat)"
echo "✓ Hermes API key stored"

echo ""
echo "=== All secrets stored successfully ==="
echo ""
echo "Now run: chezmoi apply"
echo ""
echo "Optional secrets (run 'hermes tools' to configure):"
echo "  - EXA_API_KEY"
echo "  - PARALLEL_API_KEY"
echo "  - FAL_KEY"
echo "  - OPENROUTER_API_KEY"
echo "  - GOOGLE_API_KEY"
echo "  - GROQ_API_KEY"
echo "  - GITHUB_TOKEN"
