#!/usr/bin/env python3
"""Check that the public API is documented, with a ratchet on what is not.

The check walks ``vqubit.__all__`` and the public members of every exported
class, and compares the result with ``tools/undocumented.txt``:

* a symbol without a docstring that is **not** in the baseline fails the check,
  so new public API has to be documented;
* a symbol in the baseline that now **has** a docstring is reported as stale,
  so the list only shrinks (run with ``--write`` to prune it).

The baseline is plain text, one symbol per line, so a pull request shows exactly
which names it is documenting.

Usage::

    python tools/check_api_docs.py          # check against the baseline
    python tools/check_api_docs.py --write  # rewrite the baseline
"""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import vqubit

BASELINE = Path(__file__).with_name("undocumented.txt")


def documented(value: object) -> bool:
	return bool((getattr(value, "__doc__", None) or "").strip())


def public_members() -> dict[str, bool]:
	"""Map every public symbol to whether it has a docstring."""
	result: dict[str, bool] = {}
	for name in vqubit.__all__:
		value = getattr(vqubit, name, None)
		if value is None:
			continue
		result[name] = documented(value)
		if not inspect.isclass(value):
			continue
		for attribute_name, member in vars(value).items():
			if attribute_name.startswith("_"):
				continue
			if isinstance(member, property):
				function = member.fget
			elif isinstance(member, (staticmethod, classmethod)):
				function = member.__func__
			elif inspect.isfunction(member):
				function = member
			else:
				continue
			if function is None:
				continue
			result[f"{name}.{attribute_name}"] = documented(function)
	return result


def read_baseline() -> set[str]:
	if not BASELINE.exists():
		return set()
	return {line.strip() for line in BASELINE.read_text(encoding="utf-8").splitlines() if line.strip()}


def main(argv: list[str]) -> int:
	members = public_members()
	missing = sorted(name for name, has_doc in members.items() if not has_doc)
	baseline = read_baseline()

	if "--write" in argv:
		BASELINE.write_text("\n".join(missing) + "\n", encoding="utf-8")
		print(f"wrote {len(missing)} undocumented symbols to {BASELINE}")
		return 0

	known = set(missing)
	new = sorted(known - baseline)
	fixed = sorted(baseline - known)

	total = len(members)
	documented_count = total - len(missing)
	print(f"{documented_count}/{total} public symbols have a docstring ({documented_count / total:.0%})")

	status = 0
	if new:
		print("\nThese public symbols have no docstring and are not in the baseline:")
		for name in new:
			print(f"  - {name}")
		print("Document them, or record them with `python tools/check_api_docs.py --write`.")
		status = 1
	if fixed:
		print("\nThese baseline entries now have a docstring: remove them from the baseline")
		print("with `python tools/check_api_docs.py --write`:")
		for name in fixed:
			print(f"  - {name}")
		status = 1
	return status


if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))
