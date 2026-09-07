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

for shell, name in [('bash', 'dot_bashrc'), ('zsh', 'dot_zshrc')]:
    subprocess.run([shell, '-n', str(source / name)], check=True)
    assert 'nvim' not in (source / name).read_text()
assert 'defaults write' not in (source / 'dot_bashrc').read_text()

# Exercise the actual file-picker editor helper without sourcing shell startup.
bashrc = (source / 'dot_bashrc').read_text()
helper = '_edit_file() {' + bashrc.split('_edit_file() {', 1)[1].split('\n}', 1)[0] + '\n}'
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    for command in ['zed', 'nano']:
        p = tmp / command
        p.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
        p.chmod(0o700)
    def edit():
        return subprocess.check_output(['/bin/bash', '-c', helper + '\n_edit_file "file with spaces" 12'],
                                       env=dict(os.environ, PATH=str(tmp)), text=True)
    assert edit() == '--wait\nfile with spaces:12\n'
    (tmp / 'zed').unlink()
    assert edit() == 'file with spaces\n'

    config = tmp / 'chezmoi.toml'
    config.write_text('[data.modules]\nworkstation = false\nbackupMac = false\n')
    managed = subprocess.check_output(['chezmoi', '--config', str(config), '--source', str(source),
                                       'managed'], text=True)
    assert 'restic-hclm' not in managed, 'Backup setup leaked to another machine'
    assert 'Library' not in managed.splitlines(), 'Mac-only parent directory leaked'
    assert '.config/nvim' not in managed
    assert 'tests/test_dotfiles.py' not in managed
print('PASS: Pi runtime preservation, stable preferences, invalid JSON, editor fallbacks, machine gating')
