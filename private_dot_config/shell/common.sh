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
