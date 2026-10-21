"""Prove that two builds of the same commit are byte-identical.

A reproducible build is one where the same sources produce the same wheel and
the same sdist: no timestamps, no paths, no ordering from the machine that built
them. That is what makes the hashes published on PyPI verifiable by a third
party who rebuilds from the tag, and it is the part of the 1.0 release checklist
that a script can check.

Each repetition exports the current ``HEAD`` with ``git archive``, unpacks it
into a fresh temporary directory (so the paths differ), builds a release wheel
and an sdist with vcraft, and hashes what came out. Any difference is reported
per file; the exit code is non-zero when two repetitions disagree or when an
artifact is missing.

It needs ``vcraft`` on ``PATH`` with its toolchain installed
(``vcraft toolchain install``); the CI job that runs it is
``.github/workflows/test.yml``.

Usage::

    python tools/check_reproducible.py
    python tools/check_reproducible.py --repetitions 3 --keep
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path


def _sha256(path: Path) -> str:
	return hashlib.sha256(path.read_bytes()).hexdigest()


def _artifacts(directory: Path) -> dict[str, str]:
	"""The hash of every wheel and sdist in ``directory``, keyed by name."""
	found: dict[str, str] = {}
	for path in sorted(directory.iterdir()):
		if path.suffix == ".whl" or path.name.endswith(".tar.gz"):
			found[path.name] = _sha256(path)
	return found


def _build(repetition: int, vcraft: str, root: Path) -> dict[str, str]:
	"""Export ``HEAD``, build in the export, and hash what the build wrote."""
	commit = subprocess.run(
		["git", "rev-parse", "--short", "HEAD"], cwd=root, capture_output=True, text=True, check=True
	).stdout.strip()
	archive = Path(tempfile.mkdtemp(prefix="vqubit-archive-")) / "source.tar"
	with archive.open("wb") as handle:
		subprocess.run(["git", "archive", "--format=tar", "HEAD"], cwd=root, stdout=handle, check=True)
	work = Path(tempfile.mkdtemp(prefix=f"vqubit-build-{repetition}-"))
	with tarfile.open(archive) as tar:
		tar.extractall(work, filter="data")
	for arguments in (["build", "--release", "--abi3", "3.11"], ["sdist"]):
		subprocess.run([vcraft, *arguments], cwd=work, capture_output=True, check=True)
	print(f"build {repetition}: {commit} exported to {work}")
	return _artifacts(work / "dist")


def main(argv: list[str]) -> int:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--repetitions", type=int, default=2, help="how many builds to compare (default 2)")
	parser.add_argument("--vcraft", default=None, help="the vcraft binary (default: the one on PATH)")
	parser.add_argument("--keep", action="store_true", help="leave the build directories on disk")
	options = parser.parse_args(argv)

	if options.repetitions < 2:
		print("--repetitions must be at least 2", file=sys.stderr)
		return 2
	vcraft = options.vcraft or shutil.which("vcraft")
	if not vcraft:
		print("vcraft is not on PATH; install it or pass --vcraft", file=sys.stderr)
		return 2

	root = Path(__file__).resolve().parent.parent
	results = [_build(index, vcraft, root) for index in range(options.repetitions)]
	reference = results[0]
	problems: list[str] = []
	same: dict[str, bool] = dict.fromkeys(reference, True)
	if not reference:
		problems.append("the build wrote no wheel and no sdist")
	for name, digest in sorted(reference.items()):
		for index, candidate in enumerate(results[1:], start=1):
			other = candidate.get(name)
			if other is None:
				same[name] = False
				problems.append(f"{name}: build {index} did not write it")
			elif other != digest:
				same[name] = False
				problems.append(f"{name}: build 0 has {digest[:16]}..., build {index} has {other[:16]}...")
	for candidate in results:
		for name in sorted(set(candidate) - set(reference)):
			problems.append(f"{name}: only some builds wrote it")

	for name, digest in sorted(reference.items()):
		print(f"{'identical' if same[name] else 'DIFFERENT'}  {digest[:16]}...  {name}")
	if problems:
		for problem in problems:
			print(f"::error title=Reproducible build::{problem}")
			print(problem)
		return 1
	print(f"{len(reference)} artifacts are byte-identical across {options.repetitions} builds")
	return 0


if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))
