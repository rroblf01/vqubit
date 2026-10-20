"""Additional noise channels, matching the rest of Cirq's channel surface.

The elementary channels live in :mod:`vqubit._channels`; this module adds the
asymmetric depolarizing channel, the generalized amplitude damping channel, and
the generic ``KrausChannel``/``MixedUnitaryChannel``/``RandomGateChannel``/
``StatePreparationChannel``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from ._channels import _pauli_matrix, _scale, _sqrt
from ._gates import Gate, identity
from ._linalg import Matrix
from ._value import value_equality

_I2: Matrix = [[1 + 0j, 0j], [0j, 1 + 0j]]


def _validate_probability(value: float, name: str) -> float:
	value = float(value)
	if not 0.0 <= value <= 1.0:
		raise ValueError(f"{name} must be a probability, got {value}")
	return value


@value_equality
class AsymmetricDepolarizingChannel(Gate):
	"""Depolarizes asymmetrically: applies Pauli ``i`` with probability ``p_i``."""

	def __init__(
		self,
		p_x: float | None = None,
		p_y: float | None = None,
		p_z: float | None = None,
		error_probabilities: dict[str, float] | None = None,
		tol: float = 1e-8,
	) -> None:
		if error_probabilities:
			num_qubits = len(next(iter(error_probabilities)))
			for key in error_probabilities:
				if not set(key).issubset({"I", "X", "Y", "Z"}):
					raise ValueError(f"{key} is not made solely of I, X, Y, Z.")
				if len(key) != num_qubits:
					raise ValueError(f"{key} must have {num_qubits} Pauli gates.")
			for key, value in error_probabilities.items():
				_validate_probability(value, f"p({key})")
			probs = dict(error_probabilities)
			total = sum(probs.values())
			identity_key = "I" * num_qubits
			if total < 1.0 - tol and identity_key not in probs:
				probs[identity_key] = 1.0 - total
			elif abs(total - 1.0) > tol:
				raise ValueError(f"Probabilities do not add up to 1 but to {total}")
			self._num_qubits = num_qubits
			self._error_probabilities = probs
		else:
			px = _validate_probability(0.0 if p_x is None else p_x, "p_x")
			py = _validate_probability(0.0 if p_y is None else p_y, "p_y")
			pz = _validate_probability(0.0 if p_z is None else p_z, "p_z")
			pi = 1 - _validate_probability(px + py + pz, "p_x + p_y + p_z")
			self._num_qubits = 1
			self._error_probabilities = {"I": pi, "X": px, "Y": py, "Z": pz}

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _kraus_(self) -> list[Matrix]:
		return [
			_scale(_sqrt(probability), _pauli_matrix(pauli)) for pauli, probability in self._error_probabilities.items()
		]

	def _value_equality_values_(self) -> Any:
		return (self._num_qubits, tuple(sorted(self._error_probabilities.items())))

	@property
	def error_probabilities(self) -> dict[str, float]:
		"""A dictionary from Pauli strings to probabilities."""
		return dict(self._error_probabilities)

	@property
	def p_x(self) -> float:
		"""The probability of an X error."""
		return self._error_probabilities.get("X", 0.0)

	@property
	def p_y(self) -> float:
		"""The probability of a Y error."""
		return self._error_probabilities.get("Y", 0.0)

	@property
	def p_z(self) -> float:
		"""The probability of a Z error."""
		return self._error_probabilities.get("Z", 0.0)

	def __str__(self) -> str:
		return f"asymmetric_depolarize(error_probabilities={self._error_probabilities})"

	def __repr__(self) -> str:
		return f"vqubit.AsymmetricDepolarizingChannel(error_probabilities={self._error_probabilities!r})"


def asymmetric_depolarize(
	p_x: float | None = None,
	p_y: float | None = None,
	p_z: float | None = None,
	error_probabilities: dict[str, float] | None = None,
	tol: float = 1e-8,
) -> AsymmetricDepolarizingChannel:
	"""Return an :class:`AsymmetricDepolarizingChannel` with the given parameters."""
	return AsymmetricDepolarizingChannel(p_x, p_y, p_z, error_probabilities, tol)


@value_equality
class GeneralizedAmplitudeDampingChannel(Gate):
	"""Energy dissipation and excitation, with probability ``p`` and rate ``gamma``."""

	def __init__(self, p: float, gamma: float) -> None:
		self._p = _validate_probability(p, "p")
		self._gamma = _validate_probability(gamma, "gamma")

	def _num_qubits_(self) -> int:
		return 1

	def _kraus_(self) -> list[Matrix]:
		p0 = _sqrt(self._p)
		p1 = _sqrt(1.0 - self._p)
		sg = _sqrt(self._gamma)
		sg1 = _sqrt(1.0 - self._gamma)
		return [
			[[p0, 0j], [0j, p0 * sg1]],
			[[0j, p0 * sg], [0j, 0j]],
			[[p1 * sg1, 0j], [0j, p1]],
			[[0j, 0j], [p1 * sg, 0j]],
		]

	def _value_equality_values_(self) -> Any:
		return (self._p, self._gamma)

	@property
	def p(self) -> float:
		"""The excited-state population at equilibrium."""
		return self._p

	@property
	def gamma(self) -> float:
		"""The probability that the qubit decays."""
		return self._gamma

	def __str__(self) -> str:
		return f"generalized_amplitude_damp(p={self._p!r}, gamma={self._gamma!r})"

	def __repr__(self) -> str:
		return f"vqubit.GeneralizedAmplitudeDampingChannel({self._p!r}, {self._gamma!r})"


def generalized_amplitude_damp(p: float, gamma: float) -> GeneralizedAmplitudeDampingChannel:
	"""Return a :class:`GeneralizedAmplitudeDampingChannel` with ``p`` and ``gamma``."""
	return GeneralizedAmplitudeDampingChannel(p, gamma)


@value_equality
class KrausChannel(Gate):
	"""A channel defined by an explicit list of Kraus operators."""

	def __init__(
		self, kraus_ops: Sequence[Sequence[Sequence[complex]]], key: Any = None, validate: bool = False
	) -> None:
		ops = [[[complex(z) for z in row] for row in op] for op in kraus_ops]
		if not ops:
			raise ValueError("KrausChannel must have at least one operation.")
		dim = len(ops[0])
		if dim == 0 or dim & (dim - 1) or any(len(row) != dim for row in ops[0]):
			raise ValueError(f"Input Kraus ops of shape {dim} do not represent a square operator over qubits.")
		if any(len(op) != dim for op in ops):
			raise ValueError("Inconsistent Kraus operator shapes.")
		if validate:
			from ._linalg_public import is_cptp

			if not is_cptp(kraus_ops=ops):
				raise ValueError("Kraus operators do not describe a CPTP map.")
		self._kraus_ops = ops
		self._num_qubits = dim.bit_length() - 1
		self._key = str(key) if key is not None else None

	@staticmethod
	def from_channel(channel: Gate, key: Any = None) -> KrausChannel:
		"""Return a copy of a channel with the given measurement key."""
		from ._protocols import kraus

		return KrausChannel(list(kraus(channel)), key=key)

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _kraus_(self) -> list[Matrix]:
		return self._kraus_ops

	def _value_equality_values_(self) -> Any:
		return (tuple(tuple(tuple(row) for row in op) for op in self._kraus_ops), self._key)

	def __str__(self) -> str:
		return (
			f"KrausChannel({self._kraus_ops})"
			if self._key is None
			else f"KrausChannel({self._kraus_ops}, key={self._key})"
		)

	def __repr__(self) -> str:
		return f"vqubit.KrausChannel({self._kraus_ops!r})"


@value_equality
class MixedUnitaryChannel(Gate):
	"""A mixture of unitaries, ``rho -> sum_i p_i U_i rho U_i^dagger``."""

	def __init__(
		self, mixture: Sequence[tuple[float, Sequence[Sequence[complex]]]], key: Any = None, validate: bool = False
	) -> None:
		entries = [(float(p), [[complex(z) for z in row] for row in u]) for p, u in mixture]
		if not entries:
			raise ValueError("MixedUnitaryChannel must have at least one unitary.")
		if abs(sum(p for p, _ in entries) - 1.0) > 1e-8:
			raise ValueError("Unitary probabilities must sum to 1.")
		dim = len(entries[0][1])
		if dim == 0 or dim & (dim - 1):
			raise ValueError("Input mixture does not represent a square operator over qubits.")
		if validate:
			from ._linalg_public import is_unitary

			for _, u in entries:
				if not is_unitary(u):
					raise ValueError("Element of mixture is non-unitary.")
		self._mixture = entries
		self._num_qubits = dim.bit_length() - 1
		self._key = str(key) if key is not None else None

	@staticmethod
	def from_mixture(mixture: Any, key: Any = None) -> MixedUnitaryChannel:
		"""Return a copy of a mixture with the given measurement key."""
		return MixedUnitaryChannel(list(mixture), key=key)

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _mixture_(self) -> list[tuple[float, Matrix]]:
		return [(p, [list(row) for row in u]) for p, u in self._mixture]

	def _kraus_(self) -> list[Matrix]:
		return [_scale(_sqrt(p), u) for p, u in self._mixture]

	def _value_equality_values_(self) -> Any:
		return (tuple((p, tuple(tuple(row) for row in u)) for p, u in self._mixture), self._key)

	def __str__(self) -> str:
		return f"MixedUnitaryChannel({self._mixture})"

	def __repr__(self) -> str:
		return f"vqubit.MixedUnitaryChannel({self._mixture!r})"


@value_equality
class RandomGateChannel(Gate):
	"""Applies ``sub_gate`` with probability ``probability`` and nothing otherwise."""

	def __init__(self, *, sub_gate: Gate, probability: Any) -> None:
		if isinstance(probability, (int, float)) and not 0.0 <= float(probability) <= 1.0:
			raise ValueError("not 0 <= probability <= 1")
		self._sub_gate = sub_gate
		self._probability = probability
		if isinstance(sub_gate, RandomGateChannel):
			self._probability = self._probability * sub_gate.probability
			self._sub_gate = sub_gate.sub_gate

	@property
	def sub_gate(self) -> Gate:
		"""The gate applied when the channel fires."""
		return self._sub_gate

	@property
	def probability(self) -> Any:
		"""The probability that the gate is applied."""
		return self._probability

	def _num_qubits_(self) -> int:
		return self._sub_gate.num_qubits

	def _kraus_(self) -> list[Matrix] | None:
		from ._protocols import kraus

		sub = kraus(self._sub_gate)
		if sub is None:
			return None
		dim = 1 << self._sub_gate.num_qubits
		result = [_scale(_sqrt(self._probability), k) for k in sub]
		result.append(_scale(_sqrt(1 - self._probability), identity(dim)))
		return result

	def _value_equality_values_(self) -> Any:
		return (self._sub_gate, self._probability)

	def __str__(self) -> str:
		return f"{self._sub_gate}[prob={self._probability}]"

	def __repr__(self) -> str:
		return f"vqubit.RandomGateChannel(sub_gate={self._sub_gate!r}, probability={self._probability!r})"


@value_equality
class StatePreparationChannel(Gate):
	"""Prepares a target state vector on its qubits."""

	def __init__(self, target_state: Sequence[complex], *, name: str = "StatePreparation") -> None:
		state = [complex(z) for z in target_state]
		dim = len(state)
		if dim == 0 or dim & (dim - 1):
			raise ValueError(f"Matrix width ({dim}) is not a power of 2")
		norm = math.sqrt(sum(abs(z) ** 2 for z in state))
		self._state = [z / norm for z in state]
		self._num_qubits = dim.bit_length() - 1
		self._name = name

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _kraus_(self) -> list[Matrix]:
		dim = 1 << self._num_qubits
		return [[[self._state[j] if k == i else 0j for k in range(dim)] for j in range(dim)] for i in range(dim)]

	@property
	def state(self) -> list[complex]:
		"""The state the channel prepares."""
		return list(self._state)

	def _value_equality_values_(self) -> Any:
		return (tuple(self._state), self._name)

	def __str__(self) -> str:
		return f"StatePreparationChannel({self._state})"

	def __repr__(self) -> str:
		return f"vqubit.StatePreparationChannel({self._state!r}, name={self._name!r})"


__all__ = [
	"AsymmetricDepolarizingChannel",
	"GeneralizedAmplitudeDampingChannel",
	"KrausChannel",
	"MixedUnitaryChannel",
	"RandomGateChannel",
	"StatePreparationChannel",
	"asymmetric_depolarize",
	"generalized_amplitude_damp",
]
