# ~/.config/shell/common.sh
# Shell-neutral setup shared by both zsh (~/.zshrc) and the bash fallback
# (~/.bashrc). Keep this POSIX-compatible so either shell can source it.

# ---------------------------------------------------------------------------
# Editor
# ---------------------------------------------------------------------------
# Prefer Zed, then nano.
if command -v zed >/dev/null 2>&1; then
  export EDITOR="zed --wait"
else
  export EDITOR="nano"
fi
export VISUAL="$EDITOR"
export GITHUB_USERNAME="SamuelBadr"
export BUN_INSTALL="$HOME/.bun"

# Point Homebrew at the manifest in the chezmoi source repo, so the native
# `brew bundle add/remove` commands edit the tracked file instead of a copy.
if [ -r "$HOME/.local/share/chezmoi/Brewfile" ]; then
  export HOMEBREW_BUNDLE_FILE="$HOME/.local/share/chezmoi/Brewfile"
fi

# ---------------------------------------------------------------------------
# PATH
# Prepend user-local tool dirs once each. Order (highest first) matches the
# original zsh array: bun, cargo, juliaup, julia, ~/bin, ~/.local/bin.
# Prepend in reverse so the desired first entry ends up at the front.
# ---------------------------------------------------------------------------
for _shell_path_dir in \
  "$HOME/.local/bin" \
  "$HOME/bin" \
  "$HOME/.julia/bin" \
  "$HOME/.juliaup/bin" \
  "$HOME/.cargo/bin" \
  "$BUN_INSTALL/bin"
do
  case ":$PATH:" in
    *":$_shell_path_dir:"*) ;;
    *) PATH="$_shell_path_dir:$PATH" ;;
  esac
done
export PATH
unset _shell_path_dir

# ---------------------------------------------------------------------------
# The only userland shortcuts that survive cleanup: update the system, and
# launch Pluto. Everything else was unmeasured dead weight and removed.
# ---------------------------------------------------------------------------
unalias up 2>/dev/null || true
up() {
  if command -v brew >/dev/null 2>&1; then
    brew update
    brew upgrade --greedy --overwrite -y
    brew cleanup
  fi

  if command -v pi >/dev/null 2>&1; then
    pi update --all
  fi

  if command -v skills >/dev/null 2>&1; then
    skills update -g -y
    skills update -p -y
  fi
}

alias pluto='julia -e "using Pluto; Pluto.run()"'

# ---------------------------------------------------------------------------
# pi: keep the repo manifest in sync with what is installed.
# After `pi install` / `pi remove`, copy the live `packages` list into
# .chezmoitemplates/pi-settings.json and commit+push, so every machine
# inherits the change on its next `chezmoi update`. The six runtime keys
# (defaultModel, provider, thinking...) stay machine-local via the
# modify-template. Override _PI_REPO for a relocated chezmoi source.
# ---------------------------------------------------------------------------
_pi_repo="${PI_CHEZMOI_REPO:-$HOME/.local/share/chezmoi}"
_pi_manifest="$_pi_repo/.chezmoitemplates/pi-settings.json"

_pi_sync_manifest() {
  [ -f "$_pi_manifest" ] && [ -f "$HOME/.pi/agent/settings.json" ] || return 0
  python3 - "$_pi_manifest" "$HOME/.pi/agent/settings.json" <<'PY' 2>/dev/null || return 0
import json, sys
man, live = sys.argv[1], sys.argv[2]
with open(live) as f: lj = json.load(f)
with open(man) as f: mj = json.load(f)
mj["packages"] = lj.get("packages", [])
with open(man, "w") as f:
    json.dump(mj, f, indent=2, ensure_ascii=False)
    f.write("\n")
PY
  if git -C "$_pi_repo" rev-parse --git-dir >/dev/null 2>&1; then
    git -C "$_pi_repo" add .chezmoitemplates/pi-settings.json \
      && git -C "$_pi_repo" commit -q -m "Sync pi packages" \
      && git -C "$_pi_repo" push -q
  fi
}

pi() {
  command pi "$@"; local _s=$?
  case "${1:-}" in install|remove) _pi_sync_manifest 2>/dev/null || true ;; esac
  return "$_s"
}
