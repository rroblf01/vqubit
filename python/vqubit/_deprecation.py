"""Deprecation helpers: warn once per call, with the way forward.

vqubit promises not to break a public name without a deprecation cycle (see
``docs/es/guide/stability.md``). These helpers are what a cycle is written with:
a decorator for a function or class, and a plain warning for anything else.
"""

from __future__ import annotations

import functools
import warnings
from collections.abc import Callable
from typing import Any, TypeVar

T = TypeVar("T")


def _message(name: str, since: str, replacement: str | None, removed_in: str | None) -> str:
	message = f"{name} is deprecated since vqubit {since}"
	if replacement is not None:
		message += f"; use {replacement} instead"
	if removed_in is not None:
		message += f" and will be removed in {removed_in}"
	return message + "."


def warn_deprecated(
	name: str,
	*,
	since: str,
	replacement: str | None = None,
	removed_in: str | None = None,
	stacklevel: int = 2,
) -> None:
	"""Warn that ``name`` is deprecated, pointing at its replacement.

	Args:
		name: The deprecated name, as the user writes it (``vqubit.old_thing``).
		since: The version that deprecated it.
		replacement: What to use instead, if there is one.
		removed_in: The version that will remove it, if it is scheduled.
		stacklevel: How far up the stack the warning should point.
	"""
	warnings.warn(_message(name, since, replacement, removed_in), DeprecationWarning, stacklevel=stacklevel)


def deprecated(
	name: str,
	*,
	since: str,
	replacement: str | None = None,
	removed_in: str | None = None,
) -> Callable[[T], T]:
	"""Mark a function or class as deprecated, warning on every use.

	The wrapped object keeps its name, docstring and signature, so
	``inspect.signature`` and the API snapshot keep working; instantiating a
	decorated class warns once per instantiation.
	"""
	text = _message(name, since, replacement, removed_in)

	def decorate(value: Any) -> Any:
		if isinstance(value, type):
			original_init = getattr(value, "__init__")  # noqa: B009 - mypy rejects the attribute read

			@functools.wraps(original_init)
			def __init__(self: Any, *args: Any, **kwargs: Any) -> None:
				warnings.warn(text, DeprecationWarning, stacklevel=2)
				original_init(self, *args, **kwargs)

			setattr(value, "__init__", __init__)  # noqa: B010 - mypy rejects the attribute write
			return value

		@functools.wraps(value)
		def wrapper(*args: Any, **kwargs: Any) -> Any:
			warnings.warn(text, DeprecationWarning, stacklevel=2)
			return value(*args, **kwargs)

		return wrapper

	return decorate


__all__ = ["deprecated", "warn_deprecated"]
