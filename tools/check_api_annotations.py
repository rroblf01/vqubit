#!/usr/bin/env python3
"""Check that every public callable is fully annotated.

vqubit ships ``py.typed``, so the annotations *are* API: an unannotated
parameter forces every user onto ``Any``. This walks ``vqubit.__all__`` and the
public members of every exported class and fails if a parameter (other than
``self``/``cls``) or a return value has no annotation.

Usage::

    python tools/check_api_annotations.py
"""

from __future__ import annotations

import inspect
import sys
from typing import Any

import vqubit


def missing_annotations(function: Any) -> list[str]:
	"""The parameters (and ``return``) of ``function`` that have no annotation."""
	try:
		signature = inspect.signature(function)
	except (TypeError, ValueError):
		return []
	missing = [
		parameter.name
		for parameter in signature.parameters.values()
		if parameter.annotation is inspect.Parameter.empty
		and parameter.name not in ("self", "cls")
		and parameter.kind not in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD)
	]
	if signature.return_annotation is inspect.Signature.empty:
		missing.append("return")
	return missing


def public_functions() -> dict[str, Any]:
	"""Every public function, method and property, keyed by its dotted name."""
	functions: dict[str, Any] = {}
	for name in vqubit.__all__:
		value = getattr(vqubit, name, None)
		if value is None or inspect.ismodule(value):
			continue
		if inspect.isclass(value):
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
				if function is not None:
					functions[f"{name}.{attribute_name}"] = function
		elif inspect.isfunction(value):
			functions[name] = value
	return functions


def main() -> int:
	functions = public_functions()
	missing = {name: missing_annotations(function) for name, function in functions.items()}
	missing = {name: gaps for name, gaps in missing.items() if gaps}
	print(f"{len(functions) - len(missing)}/{len(functions)} public callables are fully annotated")
	for name, gaps in sorted(missing.items()):
		print(f"  - {name}: missing {', '.join(gaps)}")
	return 1 if missing else 0


if __name__ == "__main__":
	sys.exit(main())
