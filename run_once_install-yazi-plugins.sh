#!/bin/sh
# Install Yazi packages pinned in package.toml (git rev + hash).
# Runs once per machine via chezmoi; idempotent.
if ! command -v ya >/dev/null 2>&1; then
	echo "ya (yazi package manager) not found; install yazi first, then: ya pkg install" >&2
	exit 0
fi
ya pkg install
