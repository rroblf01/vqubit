"""Helpers for working with op trees (nested iterables of operations).

An *op tree* is anything Cirq accepts where a collection of operations is
expected: an :class:`~vqubit.Operation`, a :class:`~vqubit.Moment`, a
:class:`~vqubit.Circuit`, or arbitrarily nested iterables of those. These
helpers flatten, freeze and transform such trees.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from typing import Any

from ._circuits import Circuit
from ._ops import Moment, Operation

OP_TREE = Any


def _is_iterable_tree(value: Any) -> bool:
	return isinstance(value, Iterable) and not isinstance(value, (str, bytes))


def flatten_to_ops(op_tree: Any) -> list[Operation]:
	"""Flatten an op tree into a list of operations."""
	return list(_iter_ops(op_tree))


def flatten_to_ops_or_moments(op_tree: Any) -> list[Any]:
	"""Flatten an op tree into a list of operations and moments."""
	return list(_iter_ops_or_moments(op_tree))


def _iter_ops(op_tree: Any) -> Iterator[Operation]:
	if op_tree is None:
		return
	if isinstance(op_tree, Operation):
		yield op_tree
	elif isinstance(op_tree, Moment):
		yield from op_tree.operations
	elif isinstance(op_tree, Circuit):
		yield from op_tree.all_operations()
	elif _is_iterable_tree(op_tree):
		for item in op_tree:
			yield from _iter_ops(item)
	else:
		raise TypeError(f"Not an op tree: {op_tree!r}")


def _iter_ops_or_moments(op_tree: Any) -> Iterator[Any]:
	if op_tree is None:
		return
	if isinstance(op_tree, (Operation, Moment)):
		yield op_tree
	elif isinstance(op_tree, Circuit):
		yield from op_tree.moments
	elif _is_iterable_tree(op_tree):
		for item in op_tree:
			yield from _iter_ops_or_moments(item)
	else:
		raise TypeError(f"Not an op tree: {op_tree!r}")


def flatten_op_tree(op_tree: Any, preserve_moments: bool = False) -> list[Any]:
	"""Flatten an op tree, optionally pairing each operation with its moment."""
	if not preserve_moments:
		return flatten_to_ops(op_tree)
	result: list[Any] = []
	for item in _iter_ops_or_moments(op_tree):
		if isinstance(item, Moment):
			result.extend((operation, item) for operation in item.operations)
		else:
			result.append((item, Moment([item])))
	return result


def freeze_op_tree(op_tree: Any, preserve_moments: bool = False) -> tuple[Any, ...]:
	"""Return an immutable (tuple) version of an op tree."""
	return tuple(flatten_op_tree(op_tree, preserve_moments=preserve_moments))


def transform_op_tree(
	op_tree: Any,
	op_transformation: Callable[[Operation], Any],
	moment_transformation: Callable[[Moment], Any] | None = None,
	iter_transformation: Callable[[Iterator[Any]], Any] | None = None,
) -> Any:
	"""Apply ``op_transformation`` to every operation in an op tree."""
	if op_tree is None:
		return None
	if isinstance(op_tree, Operation):
		return op_transformation(op_tree)
	if isinstance(op_tree, Moment):
		transformed = [op_transformation(operation) for operation in op_tree.operations]
		if moment_transformation is not None:
			return moment_transformation(Moment(transformed))
		return Moment(transformed)
	if isinstance(op_tree, Circuit):
		transformed = [op_transformation(operation) for operation in op_tree.all_operations()]
		if moment_transformation is not None:
			return [moment_transformation(moment) for moment in op_tree.moments]
		return Circuit(transformed)
	if _is_iterable_tree(op_tree):
		children = (transform_op_tree(item, op_transformation, moment_transformation) for item in op_tree)
		if iter_transformation is not None:
			return iter_transformation(children)
		return list(children)
	raise TypeError(f"Not an op tree: {op_tree!r}")


__all__ = [
	"OP_TREE",
	"flatten_op_tree",
	"flatten_to_ops",
	"flatten_to_ops_or_moments",
	"freeze_op_tree",
	"transform_op_tree",
]
