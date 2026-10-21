"""The archive checker really catches a bad wheel.

``tools/check_wheel.py`` guards the RECORD of what gets published, so it is
tested against wheels that are wrong in each of the ways it claims to detect.
"""

from __future__ import annotations

import base64
import hashlib
import importlib.util
import io
import tarfile
import zipfile
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).resolve().parent.parent


def _load() -> ModuleType:
	spec = importlib.util.spec_from_file_location("check_wheel", ROOT / "tools" / "check_wheel.py")
	assert spec is not None and spec.loader is not None
	module = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(module)
	return module


CHECK = _load()


def _hash(data: bytes) -> str:
	return base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()


def _record(entries: list[tuple[str, bytes]]) -> str:
	rows = [f"{name},sha256={_hash(data)},{len(data)}" for name, data in entries]
	rows.append("demo.dist-info/RECORD,,")
	return "\n".join(rows) + "\n"


def _wheel(path: Path, entries: list[tuple[str, bytes]], record: str | None = None) -> Path:
	with zipfile.ZipFile(path, "w") as archive:
		for name, data in entries:
			archive.writestr(name, data)
		archive.writestr("demo.dist-info/RECORD", record if record is not None else _record(entries))
	return path


ENTRIES = [("vqubit/__init__.py", b"x = 1\n"), ("vqubit/_core.abi3.so", b"\x7fELF fake")]


def test_a_good_wheel_has_no_problems(tmp_path: Path) -> None:
	path = _wheel(tmp_path / "demo.whl", ENTRIES)
	assert CHECK.check([path]) == []
	assert CHECK.main([str(path)]) == 0


def test_a_tampered_file_is_caught(tmp_path: Path) -> None:
	path = _wheel(tmp_path / "demo.whl", ENTRIES, _record(ENTRIES).replace(_hash(ENTRIES[0][1]), _hash(b"other")))
	problems = CHECK.check([path])
	assert any("does not match its RECORD hash" in problem for problem in problems)
	assert CHECK.main([str(path)]) == 1


def test_a_file_outside_record_is_caught(tmp_path: Path) -> None:
	path = tmp_path / "demo.whl"
	_wheel(path, ENTRIES)
	with zipfile.ZipFile(path, "a") as archive:
		archive.writestr("vqubit/extra.py", b"y = 2\n")
	assert any("not in RECORD" in problem for problem in CHECK.check([path]))


def test_a_missing_file_is_caught(tmp_path: Path) -> None:
	path = tmp_path / "demo.whl"
	with zipfile.ZipFile(path, "w") as archive:
		for name, data in ENTRIES:
			archive.writestr(name, data)
		archive.writestr("demo.dist-info/RECORD", _record([*ENTRIES, ("vqubit/gone.py", b"")]))
	assert any("which is not in the archive" in problem for problem in CHECK.check([path]))


def test_a_wheel_without_one_record_is_rejected(tmp_path: Path) -> None:
	path = tmp_path / "demo.whl"
	with zipfile.ZipFile(path, "w") as archive:
		archive.writestr("vqubit/__init__.py", b"x = 1\n")
	problems = CHECK.check([path])
	assert problems and "expected one .dist-info/RECORD" in problems[0]


def test_a_malformed_record_row_is_rejected(tmp_path: Path) -> None:
	path = _wheel(tmp_path / "demo.whl", ENTRIES, "vqubit/__init__.py,sha256=abc\n")
	try:
		CHECK.check([path])
	except ValueError as error:
		assert "malformed RECORD row" in str(error)
	else:  # pragma: no cover - the assertion below is the point
		raise AssertionError("a two-field RECORD row should be rejected")


def _sdist(path: Path, members: dict[str, bytes]) -> Path:
	with tarfile.open(path, "w:gz") as archive:
		for name, data in members.items():
			info = tarfile.TarInfo(name)
			info.size = len(data)
			archive.addfile(info, io.BytesIO(data))
	return path


def test_an_sdist_needs_a_single_root_and_its_metadata(tmp_path: Path) -> None:
	good = _sdist(
		tmp_path / "demo.tar.gz",
		{"demo-1.0/PKG-INFO": b"Metadata-Version: 2.1\n", "demo-1.0/pyproject.toml": b"[build-system]\n"},
	)
	assert CHECK.check([good]) == []

	no_pyproject = _sdist(tmp_path / "missing.tar.gz", {"demo-1.0/PKG-INFO": b"x"})
	assert any("pyproject.toml is missing" in problem for problem in CHECK.check([no_pyproject]))

	two_roots = _sdist(
		tmp_path / "two.tar.gz",
		{"demo-1.0/PKG-INFO": b"x", "other-1.0/pyproject.toml": b"y"},
	)
	assert any("single top-level directory" in problem for problem in CHECK.check([two_roots]))


def test_empty_archives_and_unknown_suffixes_are_rejected(tmp_path: Path) -> None:
	empty = _sdist(tmp_path / "empty.tar.gz", {})
	assert any("no files" in problem for problem in CHECK.check([empty]))

	other = tmp_path / "demo.zip"
	other.write_bytes(b"")
	assert any("not a wheel or an sdist" in problem for problem in CHECK.check([other]))


def test_the_command_line_reports_usage_and_missing_files(tmp_path: Path) -> None:
	assert CHECK.main([]) == 2
	assert CHECK.main([str(tmp_path / "nope.whl")]) == 1
