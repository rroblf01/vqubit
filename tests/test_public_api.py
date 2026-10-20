"""Tests for the public API contract itself.

The API is the product: these tests check that every exported name exists, that
the export list is clean, that the version is the one the package metadata
declares, and that the API snapshot used by CI is in step.
"""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path

import pytest

import vqubit

ROOT = Path(__file__).resolve().parent.parent


def test_every_exported_name_exists() -> None:
	assert vqubit.__all__, "the package exports nothing"
	missing = [name for name in vqubit.__all__ if not hasattr(vqubit, name)]
	assert missing == []


def test_the_export_list_is_clean() -> None:
	# No duplicates, no private names (a dunder such as ``__version__`` is part
	# of the public surface), no whitespace.
	assert len(set(vqubit.__all__)) == len(vqubit.__all__)
	private = [
		name for name in vqubit.__all__ if name.startswith("_") and not (name.startswith("__") and name.endswith("__"))
	]
	assert private == []
	stripped = [name for name in vqubit.__all__ if name != name.strip()]
	assert stripped == []


def test_the_version_matches_the_package_metadata() -> None:
	# ``vcraft.toml`` is where the build reads the version from, and the native
	# ``version()`` string has to agree with it.
	build = (ROOT / "vcraft.toml").read_text(encoding="utf-8")
	declared = re.search(r'^version\s*=\s*"([^"]+)"', build, re.MULTILINE)
	assert declared is not None
	assert vqubit.__version__ == declared.group(1)
	# The changelog documents the same version.
	changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
	assert f"[{vqubit.__version__}]" in changelog or "[Unreleased]" in changelog


def test_the_api_snapshot_is_current() -> None:
	# The same check CI runs: a public API change has to be declared.
	result = subprocess.run(
		[sys.executable, str(ROOT / "tools" / "check_api_surface.py")],
		capture_output=True,
		text=True,
		check=False,
	)
	assert result.returncode == 0, result.stdout + result.stderr


def test_the_docstring_ratchet_is_current() -> None:
	result = subprocess.run(
		[sys.executable, str(ROOT / "tools" / "check_api_docs.py")],
		capture_output=True,
		text=True,
		check=False,
	)
	assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("submodule", ["ops", "sim", "linalg", "protocols", "transformers", "work"])
def test_the_cirq_submodules_import(submodule: str) -> None:
	module = importlib.import_module(f"vqubit.{submodule}")
	assert module is not None
