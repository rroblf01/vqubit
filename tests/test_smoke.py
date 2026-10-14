"""Smoke tests: the package imports and its metadata is self-consistent."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import vqubit

ROOT = Path(__file__).resolve().parent.parent


def test_imports() -> None:
	assert vqubit.__version__


def test_version_matches_vcraft_toml() -> None:
	"""The native version must match the version in `vcraft.toml`."""
	text = (ROOT / "vcraft.toml").read_text(encoding="utf-8")
	match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
	assert match is not None, "vcraft.toml has no version"
	assert vqubit.__version__ == match.group(1)


@pytest.mark.parametrize("name", ["vqubit_native"])
def test_native_module_available(name: str) -> None:
	module = pytest.importorskip(name)
	assert module.version() == vqubit.__version__
