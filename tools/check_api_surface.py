#!/usr/bin/env python3
"""Snapshot the public API, so a change to it is always deliberate.

The check writes (with ``--write``) or compares (by default) a canonical text
form of ``vqubit.__all__``: the kind of each symbol, the signature of every
public function and class constructor, and the names of the public methods.

Any difference fails the check. That is the point: adding, removing or changing
a public signature must be a conscious edit to ``tools/api_snapshot.txt``, which
a reviewer can read, and the pull request then has to say whether the change is
major, minor or a patch under SemVer.

Usage::

    python tools/check_api_surface.py          # compare with the snapshot
    python tools/check_api_surface.py --write  # regenerate the snapshot
"""

from __future__ import annotations

import difflib
import inspect
import sys
import types
from pathlib import Path
from typing import Any

import vqubit

SNAPSHOT = Path(__file__).with_name("api_snapshot.txt")


def describe_signature(value: Any) -> str:
	"""The signature without annotations, which render differently per version."""
	try:
		signature = inspect.signature(value)
	except (TypeError, ValueError):
		return "(...)"
	parameters = [parameter.replace(annotation=inspect.Parameter.empty) for parameter in signature.parameters.values()]
	normalised = signature.replace(parameters=parameters, return_annotation=inspect.Signature.empty)
	return str(normalised)


def interpreter() -> str:
	return f"CPython {sys.version_info.major}.{sys.version_info.minor}"


def describe(value: Any) -> list[str]:
	"""One line per public symbol, and one per public member of a class."""
	if inspect.isclass(value):
		lines = [f"class {value.__name__}{describe_signature(value)}"]
		members = []
		for name, member in vars(value).items():
			if name.startswith("_"):
				continue
			if isinstance(member, property):
				members.append(f"    property {name}")
			elif isinstance(member, (staticmethod, classmethod)):
				kind = "staticmethod" if isinstance(member, staticmethod) else "classmethod"
				members.append(f"    {kind} {name}{describe_signature(member.__func__)}")
			elif inspect.isfunction(member):
				members.append(f"    method {name}{describe_signature(member)}")
		return lines + sorted(members)
	if isinstance(value, types.FunctionType):
		return [f"function {value.__name__}{describe_signature(value)}"]
	if isinstance(value, types.ModuleType):
		return [f"module {value.__name__}"]
	return [f"value {getattr(value, '__name__', type(value).__name__)}"]


def surface() -> list[str]:
	lines: list[str] = []
	for name in sorted(vqubit.__all__):
		value = getattr(vqubit, name, None)
		if value is None:
			continue
		for line in describe(value):
			lines.append(line if line.startswith("    ") else f"{name}: {line}")
	return lines


def main(argv: list[str]) -> int:
	lines = [f"# generated with {interpreter()}", *surface()]
	text = "\n".join(lines) + "\n"
	if "--write" in argv:
		SNAPSHOT.write_text(text, encoding="utf-8")
		print(f"wrote {len(lines)} lines to {SNAPSHOT}")
		return 0
	if not SNAPSHOT.exists():
		print(f"{SNAPSHOT} is missing; run with --write")
		return 1
	previous = SNAPSHOT.read_text(encoding="utf-8").splitlines()
	current = text.splitlines()
	if previous and previous[0].startswith("# generated with"):
		generated = previous[0].removeprefix("# generated with ").strip()
		if generated != interpreter():
			# The rendering of ``inspect.signature`` is interpreter-specific; the
			# snapshot is the gate in the lint job, which runs the version it was
			# generated with.
			print(f"snapshot generated with {generated}, running {interpreter()}: comparison skipped")
			return 0
	if previous == current:
		print(f"public API unchanged ({len(current)} lines)")
		return 0
	difference = difflib.unified_diff(
		previous, current, fromfile="tools/api_snapshot.txt", tofile="current", lineterm=""
	)
	print("the public API changed:")
	for line in difference:
		print(line)
	print(
		"\nIf the change is intended, regenerate the snapshot with "
		"`python tools/check_api_surface.py --write` and say in the pull request "
		"whether it is major, minor or a patch under SemVer."
	)
	return 1


if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))
