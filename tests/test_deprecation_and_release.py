"""Tests for the deprecation cycle and the release-notes toolchain.

The deprecation policy is only real if the mechanism warns, and the changelog
workflow is only real if towncrier can read the fragments; both are checked
here, with towncrier skipped when it is not installed (it is a development
dependency).
"""

from __future__ import annotations

import inspect
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

import vqubit

ROOT = Path(__file__).resolve().parent.parent


def test_deprecated_function_warns_and_keeps_its_contract() -> None:
	@vqubit.deprecated("vqubit.old_add", since="0.9", replacement="vqubit.add", removed_in="1.0")
	def old_add(first: int, second: int) -> int:
		"""Add two numbers the old way."""
		return first + second

	with pytest.warns(DeprecationWarning, match="deprecated since vqubit 0.9") as record:
		assert old_add(1, 2) == 3
	message = str(record[0].message)
	# The warning says what to use instead and when the name goes away.
	assert "use vqubit.add instead" in message
	assert "removed in 1.0" in message
	# The wrapper keeps the name, the docstring and the signature.
	assert old_add.__name__ == "old_add"
	assert old_add.__doc__ == "Add two numbers the old way."
	assert str(inspect.signature(old_add)) == "(first: 'int', second: 'int') -> 'int'"


def test_deprecated_class_warns_on_each_instantiation() -> None:
	@vqubit.deprecated("vqubit.OldThing", since="0.9", replacement="vqubit.NewThing")
	class OldThing:
		def __init__(self, value: int = 0) -> None:
			self.value = value

	with pytest.warns(DeprecationWarning, match="OldThing is deprecated"):
		first = OldThing(1)
	with pytest.warns(DeprecationWarning, match="use vqubit.NewThing instead"):
		second = OldThing(2)
	assert (first.value, second.value) == (1, 2)


def test_warn_deprecated_points_at_the_caller() -> None:
	def uses_old_parameter() -> None:
		vqubit.warn_deprecated("the `p` argument", since="0.9", replacement="probability")

	with warnings.catch_warnings(record=True) as caught:
		warnings.simplefilter("always")
		uses_old_parameter()
	assert caught[0].category is DeprecationWarning
	assert "use probability instead" in str(caught[0].message)
	# The warning points at the caller's line, not at the helper.
	assert caught[0].filename == __file__


def test_the_changelog_fragments_feed_the_changelog() -> None:
	towncrier = pytest.importorskip("towncrier", exc_type=ImportError)
	del towncrier
	fragments = ROOT / "changelog.d"
	assert (fragments / "README.md").exists()
	# A released tree has no pending fragments, so the test brings its own: the
	# point is the mechanism, not that there is work waiting.
	fragment = fragments / "+test.added.md"
	fragment.write_text("A fragment the test wrote.\n", encoding="utf-8")
	try:
		result = subprocess.run(
			[sys.executable, "-m", "towncrier", "build", "--draft", "--version", "0.0.0"],
			cwd=ROOT,
			capture_output=True,
			text=True,
			check=False,
		)
	finally:
		fragment.unlink()
	assert result.returncode == 0, result.stderr
	# The draft carries the fragment text, which is what the release will say.
	assert "### Added" in result.stdout
	assert "A fragment the test wrote." in result.stdout
	# Drafting must not touch the changelog.
	assert "## [0.0.0]" not in (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
