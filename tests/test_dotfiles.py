"""Offline checks. Run from anywhere with python3 tests/test_dotfiles.py."""
import json
import os
from pathlib import Path
import subprocess
import tempfile

source = Path(__file__).resolve().parents[1]
base = json.loads((source / '.chezmoitemplates/pi-settings.json').read_text())
template = (source / 'private_dot_pi/private_agent/modify_private_settings.json').read_text()

def render(text):
    return subprocess.run(['chezmoi', '--source', str(source), 'execute-template',
                           '--with-stdin', template], input=text, text=True,
                          capture_output=True)

local = dict(base, defaultProvider='local', defaultModel='chosen-model',
             defaultThinkingLevel='high', enabledModels=['local/chosen-model'],
             modelThinkingLevels={'local/chosen-model': 'high'}, lastChangelogVersion='test')
text = json.dumps(local, indent=4) + '\n'
r = render(text)
assert r.returncode == 0, r.stderr
assert r.stdout == text, 'Runtime-only edits/formatting must not drift'
local['theme'] = 'deliberate drift'
r = render(json.dumps(local))
assert r.returncode == 0, r.stderr
fixed = json.loads(r.stdout)
assert fixed['theme'] == base['theme']
assert fixed['defaultModel'] == 'chosen-model'
assert fixed['enabledModels'] == ['local/chosen-model']
assert json.loads(render('').stdout) == base
assert render('{invalid').returncode != 0, 'Invalid input must fail, not overwrite'

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
                                       'managed'], text=True)
    assert 'restic-hclm' not in managed, 'Backup setup leaked to another machine'
    assert 'Library' not in managed.splitlines(), 'Mac-only parent directory leaked'
    assert 'tests/test_dotfiles.py' not in managed
print('PASS: Pi runtime preservation, stable preferences, invalid JSON, editor fallbacks, shell syntax, common.sh wiring, machine gating')
