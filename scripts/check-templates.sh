#!/usr/bin/env bash
# Validate chezmoi templates render correctly with default data
# Run from the chezmoi source directory

set -euo pipefail

echo "==> Checking chezmoi template rendering..."

# Check modules are defined
echo -n "  - modules defined: "
if chezmoi execute-template '{{ .modules }}' >/dev/null 2>&1; then
    echo "✓"
else
    echo "✗"
    exit 1
fi

# Check for missing key errors
echo -n "  - no missing key errors: "
if chezmoi diff >/dev/null 2>&1; then
    echo "✓"
else
    echo "✗ (run 'chezmoi diff' for details)"
    exit 1
fi

# Check specific module conditionals
echo -n "  - pi_agent module conditional: "
if chezmoi execute-template '{{ if .modules.pi_agent }}PASS{{ else }}FAIL{{ end }}' | grep -q PASS; then
    echo "✓"
else
    echo "✗"
    exit 1
fi

echo ""
echo "==> All checks passed!"
echo ""
echo "Current modules:"
chezmoi execute-template '{{ range $k, $v := .modules }}{{ $k }}={{ $v }}
{{ end }}'
