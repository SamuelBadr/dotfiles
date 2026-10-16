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

# Point Homebrew at the manifest in the chezmoi source repo, so the native
# `brew bundle add/remove` commands edit the tracked file instead of a copy.
if [ -r "$HOME/.local/share/chezmoi/Brewfile" ]; then
  export HOMEBREW_BUNDLE_FILE="$HOME/.local/share/chezmoi/Brewfile"
fi

# ---------------------------------------------------------------------------
# PATH
# Move user-local tool dirs to the front, once each. Order (highest first)
# matches the original zsh array: cargo, juliaup, julia, ~/bin, ~/.local/bin.
# Prepend in reverse so the desired first entry ends up at the front.
# Node (newest nvm install, without sourcing the slow nvm.sh) has the lowest
# priority of these, so ~/bin wrappers still shadow npm-global binaries like pi.
# Existing entries are moved, not skipped: an inherited PATH (e.g. herdr panes)
# that already holds ~/bin but not node would otherwise put node in front.
# ---------------------------------------------------------------------------
export NVM_DIR="$HOME/.nvm"
# find, not a glob: zsh's nomatch aborts on an unmatched glob when nvm is absent.
_node_dir=$(find "$NVM_DIR/versions/node" -mindepth 2 -maxdepth 2 -path '*/v*/bin' 2>/dev/null | sort -V | tail -1)
[ -x "$_node_dir/node" ] || _node_dir=
for _shell_path_dir in \
  ${_node_dir:+"$_node_dir"} \
  "$HOME/.local/bin" \
  "$HOME/bin" \
  "$HOME/.julia/bin" \
  "$HOME/.juliaup/bin" \
  "$HOME/.cargo/bin"
do
  case ":$PATH:" in
    *":$_shell_path_dir:"*)
      _p=":$PATH:"; _p="${_p%%":$_shell_path_dir:"*}:${_p#*":$_shell_path_dir:"}"
      _p="${_p#:}"; PATH="${_p%:}" ;;
  esac
  PATH="$_shell_path_dir${PATH:+:$PATH}"
done
export PATH
unset _shell_path_dir _node_dir _p

# ---------------------------------------------------------------------------
# The only userland shortcut that survives cleanup: launch Pluto. (Updating
# lives in ~/bin/up, a script so `up --all` can run it over ssh.)
# ---------------------------------------------------------------------------
unalias up 2>/dev/null || true
alias pluto='julia -e "using Pluto; Pluto.run()"'

# ---------------------------------------------------------------------------
# pi: keep the repo in sync with what is installed.
# ~/.pi/agent/settings.json is tracked directly (a plain managed dotfile).
# After `pi install` / `pi remove`, re-add it so the change is captured and
# chezmoi's autoCommit/autoPush ships it; every machine inherits on its next
# `chezmoi update`. Other edits to the file behave like any dotfile: re-add to
# keep them. `pi update --all` never changes the packages list, so no sync.
# ---------------------------------------------------------------------------
pi() {
  command pi "$@"; local _s=$?
  case "${1:-}" in
    install|remove) chezmoi re-add "$HOME/.pi/agent/settings.json" >/dev/null 2>&1 || true ;;
  esac
  return "$_s"
}
