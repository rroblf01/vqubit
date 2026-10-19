"""Noise channels.

A channel is a gate that is not unitary: it describes noise, so it maps a state
to a mixed state. Channels are defined by their Kraus operators, the matrices
``K_i`` such that the state evolves as ``rho -> sum_i K_i rho K_i^dagger``. The
unitary gates of the rest of vqubit are the special case with a single Kraus
operator.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from itertools import product
from typing import Any

from ._gates import Gate
from ._linalg import Matrix, kron
from ._parameter import ParamResolver, is_symbolic, resolve_value
from ._qubits import Qid
from ._value import value_equality
from ._value_extra import validate_probability

_I: Matrix = [[1 + 0j, 0j], [0j, 1 + 0j]]
_X: Matrix = [[0j, 1 + 0j], [1 + 0j, 0j]]
_Y: Matrix = [[0j, -1j], [1j, 0j]]
_Z: Matrix = [[1 + 0j, 0j], [0j, -1 + 0j]]
_PAULIS = {"I": _I, "X": _X, "Y": _Y, "Z": _Z}


def _sqrt(value: float) -> float:
	return math.sqrt(max(value, 0.0))


def _scale(factor: float, matrix: Matrix) -> Matrix:
	return [[factor * entry for entry in row] for row in matrix]


def _pauli_matrix(letters: str) -> Matrix:
	result: Matrix = [[1 + 0j]]
	for letter in letters:
		result = kron(result, _PAULIS[letter])
	return result


class _Channel(Gate):
	"""Shared behaviour for the noise channels."""

	_symbol = "C"

	def _num_qubits_(self) -> int:
		return 1

	def _unitary_(self) -> Matrix | None:
		return None

	def _has_unitary_(self) -> bool:
		return False

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return (self._symbol,) * self.num_qubits

	def _kraus_(self) -> list[Matrix] | None:
		raise NotImplementedError

	def __str__(self) -> str:
		return type(self).__name__


@value_equality
class BitFlipChannel(_Channel):
	"""Flips the qubit to the orthogonal state with probability ``p``."""

	_symbol = "BF"

	def __init__(self, probability: Any) -> None:
		self._probability = probability if is_symbolic(probability) else validate_probability(probability, "p")

	@property
	def probability(self) -> Any:
		return self._probability

	def _kraus_(self) -> list[Matrix] | None:
		if self._is_parameterized_():
			return None
		p = self._probability
		return [_scale(_sqrt(1 - p), _I), _scale(_sqrt(p), _X)]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._probability)

	def _resolve_parameters_(self, resolver: ParamResolver) -> BitFlipChannel:
		if not self._is_parameterized_():
			return self
		return BitFlipChannel(resolve_value(self._probability, resolver))

	def _value_equality_values_(self) -> Any:
		return (self._probability,)

	def __str__(self) -> str:
		return f"bit_flip({self._probability})"


@value_equality
class PhaseFlipChannel(_Channel):
	"""Applies a Z with probability ``p``."""

	_symbol = "PF"

	def __init__(self, probability: Any) -> None:
		self._probability = probability if is_symbolic(probability) else validate_probability(probability, "p")

	@property
	def probability(self) -> Any:
		return self._probability

	def _kraus_(self) -> list[Matrix] | None:
		if self._is_parameterized_():
			return None
		p = self._probability
		return [_scale(_sqrt(1 - p), _I), _scale(_sqrt(p), _Z)]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._probability)

	def _resolve_parameters_(self, resolver: ParamResolver) -> PhaseFlipChannel:
		if not self._is_parameterized_():
			return self
		return PhaseFlipChannel(resolve_value(self._probability, resolver))

	def _value_equality_values_(self) -> Any:
		return (self._probability,)

	def __str__(self) -> str:
		return f"phase_flip({self._probability})"


@value_equality
class BitPhaseFlipChannel(_Channel):
	"""Applies a Y with probability ``p``."""

	_symbol = "BPF"

	def __init__(self, probability: Any) -> None:
		self._probability = probability if is_symbolic(probability) else validate_probability(probability, "p")

	@property
	def probability(self) -> Any:
		return self._probability

	def _kraus_(self) -> list[Matrix] | None:
		if self._is_parameterized_():
			return None
		p = self._probability
		return [_scale(_sqrt(1 - p), _I), _scale(_sqrt(p), _Y)]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._probability)

	def _resolve_parameters_(self, resolver: ParamResolver) -> BitPhaseFlipChannel:
		if not self._is_parameterized_():
			return self
		return BitPhaseFlipChannel(resolve_value(self._probability, resolver))

	def _value_equality_values_(self) -> Any:
		return (self._probability,)

	def __str__(self) -> str:
		return f"bit_phase_flip({self._probability})"


@value_equality
class DepolarizingChannel(_Channel):
	"""Applies a uniformly random non-identity Pauli with probability ``p``."""

	_symbol = "D"

	def __init__(self, probability: Any, num_qubits: int = 1) -> None:
		self._probability = probability if is_symbolic(probability) else validate_probability(probability, "p")
		self._num_qubits = int(num_qubits)

	def _num_qubits_(self) -> int:
		return self._num_qubits

	@property
	def probability(self) -> Any:
		return self._probability

	def _kraus_(self) -> list[Matrix] | None:
		if self._is_parameterized_():
			return None
		p = self._probability
		others = 4**self._num_qubits - 1
		matrices = [_scale(_sqrt(1 - p), _pauli_matrix("I" * self._num_qubits))]
		for letters in product("IXYZ", repeat=self._num_qubits):
			if all(letter == "I" for letter in letters):
				continue
			matrices.append(_scale(_sqrt(p / others), _pauli_matrix("".join(letters))))
		return matrices

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._probability)

	def _resolve_parameters_(self, resolver: ParamResolver) -> DepolarizingChannel:
		if not self._is_parameterized_():
			return self
		return DepolarizingChannel(resolve_value(self._probability, resolver), self._num_qubits)

	def _value_equality_values_(self) -> Any:
		return (self._probability, self._num_qubits)

	def __str__(self) -> str:
		return f"depolarize({self._probability})"


@value_equality
class AmplitudeDampingChannel(_Channel):
	"""Models energy loss (``|1>`` decaying to ``|0>``) with rate ``gamma``."""

	_symbol = "AD"

	def __init__(self, gamma: Any) -> None:
		self._gamma = gamma if is_symbolic(gamma) else validate_probability(gamma, "gamma")

	@property
	def gamma(self) -> Any:
		return self._gamma

	def _kraus_(self) -> list[Matrix] | None:
		if self._is_parameterized_():
			return None
		g = self._gamma
		return [
			[[1 + 0j, 0j], [0j, _sqrt(1 - g)]],
			[[0j, _sqrt(g)], [0j, 0j]],
		]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._gamma)

	def _resolve_parameters_(self, resolver: ParamResolver) -> AmplitudeDampingChannel:
		if not self._is_parameterized_():
			return self
		return AmplitudeDampingChannel(resolve_value(self._gamma, resolver))

	def _value_equality_values_(self) -> Any:
		return (self._gamma,)

	def __str__(self) -> str:
		return f"amplitude_damp({self._gamma})"


@value_equality
class PhaseDampingChannel(_Channel):
	"""Models dephasing with rate ``gamma``."""

	_symbol = "PD"

	def __init__(self, gamma: Any) -> None:
		self._gamma = gamma if is_symbolic(gamma) else validate_probability(gamma, "gamma")

	@property
	def gamma(self) -> Any:
		return self._gamma

	def _kraus_(self) -> list[Matrix] | None:
		if self._is_parameterized_():
			return None
		g = self._gamma
		return [
			[[1 + 0j, 0j], [0j, _sqrt(1 - g)]],
			[[0j, 0j], [0j, _sqrt(g)]],
		]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._gamma)

	def _resolve_parameters_(self, resolver: ParamResolver) -> PhaseDampingChannel:
		if not self._is_parameterized_():
			return self
		return PhaseDampingChannel(resolve_value(self._gamma, resolver))

	def _value_equality_values_(self) -> Any:
		return (self._gamma,)

	def __str__(self) -> str:
		return f"phase_damp({self._gamma})"


@value_equality
class ResetChannel(_Channel):
	"""Resets a qubit to ``|0>``."""

	_symbol = "R"

	def __init__(self, num_qubits: int = 1) -> None:
		self._num_qubits = int(num_qubits)

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _kraus_(self) -> list[Matrix]:
		if self._num_qubits != 1:
			raise NotImplementedError("ResetChannel Kraus operators are implemented for one qubit")
		return [
			[[1 + 0j, 0j], [0j, 0j]],
			[[0j, 1 + 0j], [0j, 0j]],
		]

	def _value_equality_values_(self) -> Any:
		return (self._num_qubits,)

	def __str__(self) -> str:
		return "reset()"


def bit_flip(probability: Any) -> BitFlipChannel:
	"""Return a bit-flip channel with the given probability."""
	return BitFlipChannel(probability)


def phase_flip(probability: Any) -> PhaseFlipChannel:
	"""Return a phase-flip channel with the given probability."""
	return PhaseFlipChannel(probability)


def bit_phase_flip(probability: Any) -> BitPhaseFlipChannel:
	"""Return a bit-phase-flip channel with the given probability."""
	return BitPhaseFlipChannel(probability)


def depolarize(probability: Any, num_qubits: int = 1) -> DepolarizingChannel:
	"""Return a depolarizing channel with the given probability."""
	return DepolarizingChannel(probability, num_qubits)


def amplitude_damp(gamma: Any) -> AmplitudeDampingChannel:
	"""Return an amplitude-damping channel with the given rate."""
	return AmplitudeDampingChannel(gamma)


def phase_damp(gamma: Any) -> PhaseDampingChannel:
	"""Return a phase-damping channel with the given rate."""
	return PhaseDampingChannel(gamma)


def reset(qubit: Qid) -> Any:
	"""Return an operation that resets ``qubit`` to ``|0>``."""
	return ResetChannel(1).on(qubit)


R = reset


def reset_each(*qubits: Qid) -> tuple[Any, ...]:
	"""Return a reset operation on each of the given qubits."""
	return tuple(reset(qubit) for qubit in qubits)
