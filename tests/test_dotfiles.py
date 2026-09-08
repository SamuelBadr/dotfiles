"""Offline checks. Run from anywhere with python3 tests/test_dotfiles.py."""
import json
import os
from pathlib import Path
import subprocess
import tempfile

source = Path(__file__).resolve().parents[1]

for shell, name in [('bash', 'dot_bashrc'), ('zsh', 'dot_zshrc'),
                    ('bash', 'dot_bash_profile'), ('bash', 'bin/executable_pi'),
                    ('bash', 'private_dot_config/shell/common.sh'),
                    ('zsh', 'private_dot_config/shell/common.sh')]:
    subprocess.run([shell, '-n', str(source / name)], check=True)
    assert 'nvim' not in (source / name).read_text()
assert 'defaults write' not in (source / 'dot_bashrc').read_text()

# Both shells must source the shared shell-neutral setup.
for rc in ('dot_zshrc', 'dot_bashrc'):
    assert '.config/shell/common.sh' in (source / rc).read_text(), f'{rc} must source common.sh'

# The shortcut layer is deliberately tiny: only `up` and `pluto` survive.
common = (source / 'private_dot_config/shell/common.sh').read_text()
bashrc = (source / 'dot_bashrc').read_text()
assert bashrc.count('alias ') == 0, 'dot_bashrc must not define aliases'
aliases = [l.split('alias ', 1)[1].split('=')[0] for l in common.splitlines() if l.startswith('alias ')]
assert aliases == ['pluto'], 'common.sh must define only the pluto alias'
assert 'up()' in common, 'common.sh must define the up function'

# Exercise the shared EDITOR selection (zed preferred, nano fallback).
editor = ('if command -v zed >/dev/null 2>&1; then\n' +
          common.split('command -v zed >/dev/null 2>&1; then', 1)[1].split('\nfi\n', 1)[0] +
          '\nfi\nprintf "%s" "$EDITOR"\n')
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    for command in ['zed', 'nano']:
        p = tmp / command
        p.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
        p.chmod(0o700)
    out = subprocess.check_output(['/bin/bash', '-c', editor],
                                  env=dict(os.environ, PATH=str(tmp)), text=True)
    assert out == 'zed --wait', out
    (tmp / 'zed').unlink()
    out = subprocess.check_output(['/bin/bash', '-c', editor],
                                  env=dict(os.environ, PATH=str(tmp)), text=True)
    assert out == 'nano', out

    config = tmp / 'chezmoi.toml'
    config.write_text('[data.modules]\nworkstation = false\nbackupMac = false\n')
    managed = subprocess.check_output(['chezmoi', '--config', str(config), '--source', str(source),
                                       'managed'], text=True).splitlines()
    assert 'restic-hclm' not in managed, 'Backup setup leaked to another machine'
    assert 'Library' not in managed, 'Mac-only parent directory leaked'
    assert 'tests/test_dotfiles.py' not in managed
    assert 'Brewfile' not in managed, 'Repo-only Homebrew manifest must not deploy'
    assert '.pi/agent/settings.json' in managed, 'Pi settings must be tracked as a plain file'

# Pi settings: tracked directly, not templated, must carry the current package set.
pi_settings = (source / 'private_dot_pi/private_agent/settings.json').read_text()
assert '{{' not in pi_settings, 'Pi settings must not be templated'
pi_json = json.loads(pi_settings)
assert 'https://github.com/ayghri/i-have-adhd' in pi_json['packages']
assert 'git:github.com/jonjonrankin/pi-caveman' not in pi_json['packages']

# Homebrew parses the manifest itself, so it must stay free of template syntax.
brewfile = (source / 'Brewfile').read_text()
assert '{{' not in brewfile, 'Brewfile must not be templated'
assert 'HOMEBREW_BUNDLE_FILE' in common, 'common.sh must point brew at the repo manifest'

# The `pi` wrapper re-adds settings.json after install/remove (and after no other verb).
def pi_dispatch_check():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        h = tmp / 'h'
        (h / '.pi' / 'agent').mkdir(parents=True)
        (h / '.pi' / 'agent' / 'settings.json').write_text('{}')
        fb = tmp / 'bin'
        fb.mkdir()
        trace = tmp / 'trace'
        (fb / 'pi').write_text('#!/bin/sh\necho "pi $*" >> "$TRACE"\nexit 0\n')
        (fb / 'chezmoi').write_text('#!/bin/sh\necho "chezmoi $*" >> "$TRACE"\nexit 0\n')
        for p in (fb / 'pi', fb / 'chezmoi'):
            p.chmod(0o700)
        env = dict(os.environ, HOME=str(h), TRACE=str(trace), PATH=f'{fb}:/usr/bin:/bin')
        common_path = source / 'private_dot_config/shell/common.sh'
        def run(v):
            return subprocess.run(['/bin/bash', '-c', f'. "{common_path}"; pi {v}'],
                                  env=env, capture_output=True, text=True).returncode
        assert run('install npm:x') == 0
        assert run('remove npm:x') == 0
        assert run('update --all') == 0
        calls = trace.read_text()
        assert 'pi install npm:x' in calls and 'pi remove npm:x' in calls
        assert calls.count('re-add') == 2, f'expected re-add only after install/remove:\n{calls}'
pi_dispatch_check()

print('PASS: shell syntax, common.sh wiring, editor fallbacks, machine gating, pi settings tracking, brew hygiene')
