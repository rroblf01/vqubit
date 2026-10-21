"""Verify the ``RECORD`` of a built wheel (or sdist).

``RECORD`` is the wheel's own manifest: one line per file with its ``sha256``
and size. ``pip`` does not check it, and nothing in the build path proves that
the archive a user downloads contains exactly the files that were recorded, so
a truncated or post-processed wheel could ship silently. This reads the archive
back, recomputes every hash and flags:

* a recorded file whose bytes or size do not match,
* a file in the archive that ``RECORD`` does not list,
* a recorded file that is missing from the archive.

That covers the "hash verification" part of the 1.0 reproducibility item. The
other part —two builds of the same commit producing the same bytes— needs
``SOURCE_DATE_EPOCH`` support in the packaging step and is tracked upstream
(``UPSTREAM.md`` §5); this script is what will prove it once it lands, by
comparing the digest of two builds.

Usage::

    python tools/check_wheel.py dist/*.whl
"""

from __future__ import annotations

import base64
import csv
import hashlib
import io
import sys
import tarfile
import zipfile
from collections.abc import Iterable
from pathlib import Path


def _digest(data: bytes) -> str:
	"""The ``RECORD`` spelling of a SHA-256 digest."""
	return base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode("ascii")


def _record_rows(text: str) -> list[tuple[str, str, str]]:
	"""Every ``RECORD`` row, including the one for ``RECORD`` itself."""
	rows: list[tuple[str, str, str]] = []
	for row in csv.reader(io.StringIO(text)):
		if not row:
			continue
		if len(row) != 3:
			raise ValueError(f"malformed RECORD row: {row!r}")
		rows.append((row[0], row[1], row[2]))
	return rows


def _check_wheel(path: Path) -> list[str]:
	problems: list[str] = []
	with zipfile.ZipFile(path) as archive:
		names = archive.namelist()
		records = [name for name in names if name.endswith(".dist-info/RECORD")]
		if len(records) != 1:
			return [f"{path.name}: expected one .dist-info/RECORD, found {len(records)}"]
		rows = _record_rows(archive.read(records[0]).decode("utf-8"))
		recorded = {name: (hash_value, size) for name, hash_value, size in rows}
		for name, hash_value, size in rows:
			if name not in names:
				problems.append(f"{path.name}: RECORD lists {name}, which is not in the archive")
				continue
			data = archive.read(name)
			# RECORD lists itself with no hash and no size; every other file has both.
			if hash_value and hash_value != f"sha256={_digest(data)}":
				problems.append(f"{path.name}: {name} does not match its RECORD hash")
			if size and size != str(len(data)):
				problems.append(f"{path.name}: {name} is {len(data)} bytes, RECORD says {size}")
		for name in names:
			if name not in recorded and not name.endswith("/"):
				problems.append(f"{path.name}: {name} is in the archive but not in RECORD")
	return problems


def _check_sdist(path: Path) -> list[str]:
	problems: list[str] = []
	with tarfile.open(path) as archive:
		members = [member for member in archive.getmembers() if member.isfile()]
		if not members:
			return [f"{path.name}: the archive has no files"]
		roots = {member.name.split("/", 1)[0] for member in members}
		if len(roots) != 1:
			problems.append(f"{path.name}: expected a single top-level directory, found {sorted(roots)}")
		for required in ("PKG-INFO", "pyproject.toml"):
			if not any(member.name.endswith("/" + required) for member in members):
				problems.append(f"{path.name}: {required} is missing")
	return problems


def check(paths: Iterable[Path]) -> list[str]:
	problems: list[str] = []
	for path in paths:
		if path.suffix == ".whl":
			problems.extend(_check_wheel(path))
		elif path.name.endswith(".tar.gz"):
			problems.extend(_check_sdist(path))
		else:
			problems.append(f"{path.name}: not a wheel or an sdist")
	return problems


def main(argv: list[str]) -> int:
	paths = [Path(argument) for argument in argv]
	if not paths:
		print("usage: python tools/check_wheel.py <dist/*.whl> [...]")
		return 2
	missing = [str(path) for path in paths if not path.exists()]
	if missing:
		print("missing: " + ", ".join(missing))
		return 1
	problems = check(paths)
	for problem in problems:
		print(f"::error title=Archive check::{problem}")
		print(problem)
	if problems:
		return 1
	for path in paths:
		print(f"{path.name}: RECORD verified")
	return 0


if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))
