"""Offline updater checks: python3 -m unittest discover -s tests -v."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

UP = Path(__file__).resolve().parents[1] / "bin/executable_up"

# Each fake updater records its lifetime. No real package manager or SSH runs.
STUB = r'''
import os
from pathlib import Path
import sys
import time
name = Path(sys.argv[0]).name
args = sys.argv[1:]
if "--version" in args:
    print("1.2.3")
    sys.exit(0)
if name == "pgrep":
    sys.exit(1)
if name == "cargo" and args == ["install", "--list"]:
    sys.exit(0)
label = name + (" " + args[0] if args else "")
def record(event):
    with open(os.environ["EVENTS"], "a") as f:
        f.write(event + " " + label + "\n")
record("start")
time.sleep(0.15)
record("end")
if os.environ.get("FAIL") == name:
    print(name + ": deliberate failure")
    sys.exit(1)
'''


class UpTests(unittest.TestCase):
    def run_up(self, fail="", sync=True):
        with tempfile.TemporaryDirectory() as tmp:
            home = Path(tmp)
            bin_dir = home / "bin"
            bin_dir.mkdir()
            config = home / ".config/shell/common.sh"
            config.parent.mkdir(parents=True)
            config.write_text("")
            # Restrict PATH: installed tools on the test machine must not run.
            for name in ("bash", "mktemp", "sed", "grep", "tr", "cat", "rm",
                         "awk", "head", "env"):
                (bin_dir / name).symlink_to(shutil.which(name))
            for name in ("brew", "pi", "skills-sync", "skills", "rustup",
                         "rustc", "cargo", "uv", "juliaup", "julia", "pgrep"):
                if name == "skills-sync" and not sync:
                    continue
                p = bin_dir / name
                p.write_text(f"#!{sys.executable}\n" + STUB)
                p.chmod(0o755)
            events = home / "events"
            result = subprocess.run(
                [shutil.which("bash"), str(UP)], capture_output=True, text=True,
                env=os.environ | {"HOME": tmp, "PATH": str(bin_dir),
                                  "EVENTS": str(events), "FAIL": fail}, timeout=10,
            )
            return result, events.read_text().splitlines()

    def test_parallel_updates_and_safe_brew_upgrade(self):
        result, events = self.run_up()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("versions: julia 1.2.3 | rustc 1.2.3 | pi 1.2.3", result.stdout)
        self.assertNotIn("start skills update", events)
        self.assertEqual(events.count("start skills-sync"), 1)
        # All lanes start before the foreground Julia update ends.
        julia_end = events.index("end julia --startup-file=no")
        for label in ("brew update", "pi update", "skills-sync", "rustup update", "uv tool"):
            self.assertLess(events.index("start " + label), julia_end)
        # Homebrew must not replace binaries while any updater is using them.
        upgrade = events.index("start brew upgrade")
        for i, event in enumerate(events):
            if event.startswith("end ") and event != "end brew upgrade":
                self.assertLess(i, upgrade)

    def test_failures_reach_parent_without_blocking_other_updates(self):
        for name in ("pi", "skills-sync", "rustup", "uv", "julia", "brew"):
            with self.subTest(name=name):
                result, events = self.run_up(fail=name)
                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                self.assertIn("up: failed:", result.stderr)
                self.assertIn(name + ": deliberate failure", result.stdout + result.stderr)
                self.assertIn("end uv tool", events)
                self.assertEqual(events.count("start pi update"), 2)  # self/packages + models

    def test_skills_fallback_updates_global_skills_once(self):
        result, events = self.run_up(sync=False)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(events.count("start skills update"), 1)

    def apply_dotfiles(self, diff, push):
        """Run apply_dotfiles on one file changed on both sides; return chezmoi's writes."""
        source = UP.read_text()
        function = source[source.index("apply_dotfiles() {"):source.index("\n}\n", source.index("apply_dotfiles() {")) + 3]
        stub = r'''
import os, sys
args = sys.argv[1:]
if args[0] == "status": print("MM /f")
elif args[0] == "diff": sys.stdout.write(os.environ["DIFF"])
else: open(os.environ["EVENTS"], "a").write(" ".join(args) + "\n")
'''
        with tempfile.TemporaryDirectory() as tmp:
            home = Path(tmp)
            (home / "bin").mkdir()
            chezmoi = home / "bin/chezmoi"
            chezmoi.write_text(f"#!{sys.executable}\n" + stub)
            chezmoi.chmod(0o755)
            if push:
                (home / ".config/chezmoi").mkdir(parents=True)
                (home / ".config/chezmoi/chezmoi.toml").write_text("[git]\n    autoPush = true\n")
            events = home / "events"
            events.touch()
            subprocess.run([shutil.which("bash"), "-c", function + "\napply_dotfiles"], check=True,
                           env=os.environ | {"HOME": tmp, "PATH": f"{home / 'bin'}:{os.environ['PATH']}",
                                             "EVENTS": str(events), "DIFF": diff}, timeout=10)
            return [e for e in events.read_text().splitlines() if not e.startswith("apply --no-tty")]

    def test_apply_dotfiles_resolves_pi_changelog_markers(self):
        head = "--- a/settings.json\n+++ b/settings.json\n@@ -1 +1 @@\n"
        marker = ('-    "{\\"source\\":\\"npm:x\\",\\"scope\\":\\"user\\"}": "0.48.0",\n'
                  '+    "{\\"source\\":\\"npm:x\\",\\"scope\\":\\"user\\"}": "0.47.0",\n')
        real = '-    "theme": "dark",\n'
        self.assertEqual(self.apply_dotfiles(head + marker, push=True), ["re-add /f"])
        self.assertEqual(self.apply_dotfiles(head + marker, push=False), ["apply --force /f"])
        self.assertEqual(self.apply_dotfiles(head + marker + real, push=True), [])
        self.assertEqual(self.apply_dotfiles(head + real, push=False), [])

    @unittest.skipUnless(shutil.which("julia"), "Julia not installed")
    def test_shared_environment_selection(self):
        source = UP.read_text()
        selection = source[source.index("    envs = "):source.index("    function quiet")]
        setup = '''
        using TOML, Dates
        empty!(DEPOT_PATH); push!(DEPOT_PATH, ARGS[1])
        root = joinpath(DEPOT_PATH[1], "environments")
        current = "v$(VERSION.major).$(VERSION.minor)"
        usage = Dict{String, Any}()
        for name in [current, "v0.1", "tool", "stale", "coveragetempenv", "__ci"]
            path = joinpath(root, name)
            mkpath(path); write(joinpath(path, "Project.toml"), "[deps]\\n")
            usage[joinpath(path, "Manifest.toml")] = [Dict("time" =>
                now(UTC) - Day(name == "stale" ? 31 : 1))]
        end
        log = joinpath(DEPOT_PATH[1], "logs", "manifest_usage.toml")
        mkpath(dirname(log)); open(io -> TOML.print(io, usage), log, "w")
        '''
        with tempfile.TemporaryDirectory() as depot:
            result = subprocess.run(
                [shutil.which("julia"), "--startup-file=no", "-e",
                 setup + selection + '\n@assert sort(shared) == sort([current, "tool"])', depot],
                capture_output=True, text=True, timeout=15,
                env=os.environ | {"OPENBLAS_NUM_THREADS": "1", "JULIA_NUM_THREADS": "1"},
            )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
