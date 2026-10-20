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
	"""Flatten an op tree into a list of operations.

	>>> import vqubit
	>>> q0 = vqubit.LineQubit(0)
	>>> len(flatten_to_ops([vqubit.X(q0), [vqubit.Y(q0), vqubit.Z(q0)]]))
	3
	"""
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
	root: Any,
	op_transformation: Callable[[Operation], Any] = lambda operation: operation,
	iter_transformation: Callable[[Iterator[Any]], Any] = lambda items: list(items),
	preserve_moments: bool = False,
	*,
	moment_transformation: Callable[[Moment], Any] | None = None,
) -> Any:
	"""Apply the transformation functions to the nodes of an op tree.

	``op_transformation`` rewrites every operation (the leaves);
	``iter_transformation`` rewrites every iterable (the internal nodes), so
	passing ``tuple`` builds an immutable tree. ``preserve_moments`` leaves
	moments and their contents untouched, as Cirq's flag does.

	As an extension, ``moment_transformation`` (keyword-only) rewrites each
	moment after its operations have been transformed, which is what a caller
	that wants to inspect or rebuild the moment needs.
	"""
	if root is None:
		return None
	if isinstance(root, Operation):
		return op_transformation(root)
	if preserve_moments and isinstance(root, Moment):
		return root
	if isinstance(root, Moment):
		transformed = [op_transformation(operation) for operation in root.operations]
		moment = Moment(transformed)
		return moment_transformation(moment) if moment_transformation is not None else moment
	if isinstance(root, Circuit):
		moments: list[Moment] = [
			moment if preserve_moments else Moment([op_transformation(operation) for operation in moment.operations])
			for moment in root.moments
		]
		if moment_transformation is not None and not preserve_moments:
			# The caller asked to see every moment and decide what to build (a
			# list of lengths, a new circuit, ...), so hand back their results.
			return [moment_transformation(moment) for moment in moments]
		return Circuit.from_moments(moments)
	if _is_iterable_tree(root):
		children = (
			transform_op_tree(
				item,
				op_transformation,
				iter_transformation,
				preserve_moments,
				moment_transformation=moment_transformation,
			)
			for item in root
		)
		return iter_transformation(children)
	raise TypeError(f"Not an op tree: {root!r}")


__all__ = [
	"OP_TREE",
	"flatten_op_tree",
	"flatten_to_ops",
	"flatten_to_ops_or_moments",
	"freeze_op_tree",
	"transform_op_tree",
]
