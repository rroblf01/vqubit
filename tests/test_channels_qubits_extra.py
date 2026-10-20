"""Extra tests for noise channels and qubit identifiers.

``test_channels.py`` and ``test_qubits.py`` cover the common paths. Here we pin
down the Kraus matrices themselves, their Choi representations, parameter
resolution and the less-used qid variants.
"""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close
from vqubit._linalg import dagger, matmul

I2 = [[1 + 0j, 0j], [0j, 1 + 0j]]
X2 = [[0j, 1 + 0j], [1 + 0j, 0j]]
Y2 = [[0j, -1j], [1j, 0j]]
Z2 = [[1 + 0j, 0j], [0j, -1 + 0j]]
ZERO2 = [[0j, 0j], [0j, 0j]]


def _identity(dim: int) -> list[list[complex]]:
	return [[1 + 0j if i == j else 0j for j in range(dim)] for i in range(dim)]


def _assert_trace_preserving(operators: list[list[list[complex]]], dim: int) -> None:
	total = [[0j] * dim for _ in range(dim)]
	for operator in operators:
		product = matmul(dagger(operator), operator)
		for i in range(dim):
			for j in range(dim):
				total[i][j] += product[i][j]
	assert_matrix_close(total, _identity(dim))


def _apply_kraus(operators: list[list[list[complex]]], rho: list[list[complex]]) -> list[list[complex]]:
	dim = len(rho)
	result = [[0j] * dim for _ in range(dim)]
	for operator in operators:
		term = matmul(matmul(operator, rho), dagger(operator))
		for i in range(dim):
			for j in range(dim):
				result[i][j] += term[i][j]
	return result


def test_bit_flip_kraus_at_the_ends_of_the_probability_range() -> None:
	zero = vqubit.kraus(vqubit.bit_flip(0.0))
	assert_matrix_close(zero[0], I2)
	assert_matrix_close(zero[1], ZERO2)

	full = vqubit.kraus(vqubit.bit_flip(1.0))
	assert_matrix_close(full[0], ZERO2)
	assert_matrix_close(full[1], X2)


def test_bit_flip_at_half_probability_is_an_even_pauli_mixture() -> None:
	operators = vqubit.kraus(vqubit.bit_flip(0.5))
	assert_matrix_close(operators[0], [[math.sqrt(0.5) + 0j, 0j], [0j, math.sqrt(0.5) + 0j]])
	assert_matrix_close(operators[1], [[0j, math.sqrt(0.5) + 0j], [math.sqrt(0.5) + 0j, 0j]])


def test_phase_flip_and_bit_phase_flip_kraus_operators() -> None:
	phase = vqubit.kraus(vqubit.phase_flip(1.0))
	assert_matrix_close(phase[0], ZERO2)
	assert_matrix_close(phase[1], Z2)

	bit_phase = vqubit.kraus(vqubit.bit_phase_flip(1.0))
	assert_matrix_close(bit_phase[0], ZERO2)
	assert_matrix_close(bit_phase[1], Y2)

	half = vqubit.kraus(vqubit.phase_flip(0.5))
	assert_matrix_close(half[1], [[math.sqrt(0.5) + 0j, 0j], [0j, -math.sqrt(0.5) + 0j]])


def test_amplitude_damping_kraus_operators() -> None:
	gamma = 0.25
	operators = vqubit.kraus(vqubit.amplitude_damp(gamma))
	assert_matrix_close(
		operators[0],
		[[1 + 0j, 0j], [0j, math.sqrt(1 - gamma) + 0j]],
	)
	assert_matrix_close(operators[1], [[0j, math.sqrt(gamma) + 0j], [0j, 0j]])

	# At gamma = 1 the excited state decays completely.
	full = vqubit.kraus(vqubit.amplitude_damp(1.0))
	assert_matrix_close(full[0], [[1 + 0j, 0j], [0j, 0j]])
	assert_matrix_close(full[1], [[0j, 1 + 0j], [0j, 0j]])


def test_phase_damping_kraus_operators() -> None:
	gamma = 0.25
	operators = vqubit.kraus(vqubit.phase_damp(gamma))
	assert_matrix_close(operators[0], [[1 + 0j, 0j], [0j, math.sqrt(1 - gamma) + 0j]])
	assert_matrix_close(operators[1], [[0j, 0j], [0j, math.sqrt(gamma) + 0j]])

	# At gamma = 1 the channel is a fully dephasing measurement of Z.
	full = vqubit.kraus(vqubit.phase_damp(1.0))
	assert_matrix_close(full[0], [[1 + 0j, 0j], [0j, 0j]])
	assert_matrix_close(full[1], [[0j, 0j], [0j, 1 + 0j]])


def test_depolarizing_channel_kraus_operators() -> None:
	full = vqubit.kraus(vqubit.depolarize(1.0))
	assert len(full) == 4
	assert_matrix_close(full[0], ZERO2)
	scale = 1 / math.sqrt(3)
	assert_matrix_close(full[1], [[0j, scale + 0j], [scale + 0j, 0j]])
	assert_matrix_close(full[2], [[0j, -scale * 1j], [scale * 1j, 0j]])
	assert_matrix_close(full[3], [[scale + 0j, 0j], [0j, -scale + 0j]])

	# At p = 0 only the identity operator survives.
	none = vqubit.kraus(vqubit.depolarize(0.0))
	assert len(none) == 4
	assert_matrix_close(none[0], I2)
	for operator in none[1:]:
		assert_matrix_close(operator, ZERO2)


def test_depolarizing_channel_on_two_qubits_has_sixteen_kraus_operators() -> None:
	probability = 0.4
	operators = vqubit.kraus(vqubit.depolarize(probability, 2))
	assert len(operators) == 16
	assert_matrix_close(operators[0], [[math.sqrt(1 - probability) * entry for entry in row] for row in _identity(4)])
	# Operator 5 is XX (the first non-identity product after IX, IY, IZ, XI),
	# scaled by sqrt(p / 15).
	scaled_xx = [[math.sqrt(probability / 15) * entry for entry in row] for row in vqubit.kron(X2, X2)]
	assert_matrix_close(operators[5], scaled_xx)
	assert vqubit.is_cptp(kraus_ops=operators)
	_assert_trace_preserving(operators, 4)


def test_every_channel_is_trace_preserving() -> None:
	channels = [
		vqubit.bit_flip(0.2),
		vqubit.phase_flip(0.3),
		vqubit.bit_phase_flip(0.4),
		vqubit.depolarize(0.15),
		vqubit.amplitude_damp(0.25),
		vqubit.phase_damp(0.35),
	]
	for channel in channels:
		operators = vqubit.kraus(channel)
		assert vqubit.is_cptp(kraus_ops=operators), str(channel)
		_assert_trace_preserving(operators, 2)


def test_channels_map_density_matrices_through_their_kraus_operators() -> None:
	gamma = 0.3
	# An excited state loses its population at rate gamma.
	damped = _apply_kraus(vqubit.kraus(vqubit.amplitude_damp(gamma)), [[0j, 0j], [0j, 1 + 0j]])
	assert_matrix_close(damped, [[gamma + 0j, 0j], [0j, 1 - gamma + 0j]])

	# A bit flip moves population from |0> to |1> with probability p.
	flipped = _apply_kraus(vqubit.kraus(vqubit.bit_flip(0.25)), [[1 + 0j, 0j], [0j, 0j]])
	assert_matrix_close(flipped, [[0.75 + 0j, 0j], [0j, 0.25 + 0j]])

	# At p = 3/4 the depolarizing channel sends any state to the maximally
	# mixed state (p counts the probability of a non-identity Pauli).
	mixed = _apply_kraus(vqubit.kraus(vqubit.depolarize(0.75)), [[1 + 0j, 0j], [0j, 0j]])
	assert_matrix_close(mixed, [[0.5 + 0j, 0j], [0j, 0.5 + 0j]])
	# For p = 1 the identity term vanishes but the state is not yet maximally
	# mixed: it is the uniform mixture of X, Y and Z errors.
	fully_depolarized = _apply_kraus(vqubit.kraus(vqubit.depolarize(1.0)), [[1 + 0j, 0j], [0j, 0j]])
	assert_matrix_close(fully_depolarized, [[1 / 3 + 0j, 0j], [0j, 2 / 3 + 0j]])


def test_choi_matrix_is_trace_preserving_and_roundtrips() -> None:
	for channel in (vqubit.bit_flip(0.3), vqubit.amplitude_damp(0.4), vqubit.phase_flip(0.2)):
		choi = vqubit.kraus_to_choi(vqubit.kraus(channel))
		assert len(choi) == 4
		assert sum(choi[i][i] for i in range(4)) == pytest.approx(2)
		recovered = vqubit.kraus_to_choi(vqubit.choi_to_kraus(choi))
		assert_matrix_close(recovered, choi)

	two_qubit_choi = vqubit.kraus_to_choi(vqubit.kraus(vqubit.depolarize(0.2, 2)))
	assert len(two_qubit_choi) == 16
	assert sum(two_qubit_choi[i][i] for i in range(16)) == pytest.approx(4)


def test_channels_have_no_unitary_but_do_have_kraus_operators() -> None:
	channel = vqubit.bit_flip(0.1)
	assert channel.num_qubits == 1
	assert channel.qid_shape == (2,)
	assert not vqubit.has_unitary(channel)
	assert vqubit.has_kraus(channel)
	assert not vqubit.is_measurement(channel)
	with pytest.raises(TypeError):
		vqubit.unitary(channel)


def test_channel_operations_diagrams_and_strings() -> None:
	q0 = vqubit.LineQubit(0)
	operation = vqubit.bit_flip(0.1)(q0)
	assert operation.qubits == (q0,)
	assert operation.gate == vqubit.bit_flip(0.1)
	assert vqubit.num_qubits(operation) == 1

	expected = [
		(vqubit.bit_flip(0.1), ("BF",), "bit_flip(0.1)"),
		(vqubit.phase_flip(0.1), ("PF",), "phase_flip(0.1)"),
		(vqubit.bit_phase_flip(0.1), ("BPF",), "bit_phase_flip(0.1)"),
		(vqubit.depolarize(0.1), ("D",), "depolarize(0.1)"),
		(vqubit.amplitude_damp(0.1), ("AD",), "amplitude_damp(0.1)"),
		(vqubit.phase_damp(0.1), ("PD",), "phase_damp(0.1)"),
		(vqubit.ResetChannel(), ("R",), "reset()"),
	]
	for channel, symbols, text in expected:
		assert vqubit.circuit_diagram_info(channel).wire_symbols == symbols
		assert str(channel) == text
	assert vqubit.circuit_diagram_info(vqubit.depolarize(0.1, 2)).wire_symbols == ("D", "D")


def test_channel_equality_and_hashing_are_type_specific() -> None:
	assert vqubit.bit_flip(0.1) == vqubit.BitFlipChannel(0.1)
	assert vqubit.bit_flip(0.1) != vqubit.bit_flip(0.2)
	assert vqubit.bit_flip(0.1) != vqubit.phase_flip(0.1)
	assert hash(vqubit.bit_flip(0.1)) == hash(vqubit.bit_flip(0.1))
	assert len({vqubit.bit_flip(0.1), vqubit.bit_flip(0.1)}) == 1
	assert vqubit.depolarize(0.1, 2) != vqubit.depolarize(0.1, 1)


def test_channel_parameterization_and_resolution() -> None:
	parameter = vqubit.Parameter("t")
	resolver = vqubit.ParamResolver({parameter: 0.2})
	builders = [
		vqubit.bit_flip,
		vqubit.phase_flip,
		vqubit.bit_phase_flip,
		vqubit.amplitude_damp,
		vqubit.phase_damp,
	]
	for build in builders:
		channel = build(parameter)
		assert vqubit.is_parameterized(channel), build.__name__
		assert not vqubit.has_kraus(channel)
		with pytest.raises(TypeError):
			vqubit.kraus(channel)
		resolved = vqubit.resolve_parameters(channel, resolver)
		assert resolved == build(0.2)
		assert not vqubit.is_parameterized(resolved)

	two_qubit = vqubit.depolarize(parameter, 2)
	assert vqubit.is_parameterized(two_qubit)
	resolved = vqubit.resolve_parameters(two_qubit, resolver)
	assert resolved == vqubit.depolarize(0.2, 2)
	assert resolved.num_qubits == 2


def test_reset_channel_kraus_qid_shape_and_validation() -> None:
	channel = vqubit.ResetChannel()
	assert channel.num_qubits == 1
	assert channel.qid_shape == (2,)
	operators = vqubit.kraus(channel)
	assert len(operators) == 2
	assert_matrix_close(operators[0], [[1 + 0j, 0j], [0j, 0j]])
	assert_matrix_close(operators[1], [[0j, 1 + 0j], [0j, 0j]])
	_assert_trace_preserving(operators, 2)

	qutrit = vqubit.ResetChannel(3)
	assert qutrit.qid_shape == (3,)
	assert qutrit.num_qubits == 1
	qutrit_operators = vqubit.kraus(qutrit)
	assert len(qutrit_operators) == 3
	_assert_trace_preserving(qutrit_operators, 3)

	with pytest.raises(ValueError):
		vqubit.ResetChannel(1)
	assert vqubit.ResetChannel(2) == vqubit.ResetChannel()
	assert vqubit.ResetChannel(2) != vqubit.ResetChannel(3)


def test_reset_helpers_build_operations_on_qubits() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.reset(q0)
	assert operation.qubits == (q0,)
	assert isinstance(operation.gate, vqubit.ResetChannel)
	assert vqubit.R is vqubit.reset
	# reset_each is re-exported from both _channels and _gates_more.
	assert vqubit.reset_each(q0, q1) == (vqubit.reset(q0), vqubit.reset(q1))
	grid_operation = vqubit.reset(vqubit.GridQubit(1, 2))
	assert grid_operation.qubits == (vqubit.GridQubit(1, 2),)


def test_qid_dimension_and_hierarchy() -> None:
	assert vqubit.Qid().dimension == 2
	assert vqubit.LineQubit(0).dimension == 2
	assert vqubit.GridQubit(0, 0).dimension == 2
	assert vqubit.NamedQubit("a").dimension == 2
	assert isinstance(vqubit.LineQubit(0), vqubit.Qubit)
	assert isinstance(vqubit.GridQubit(0, 0), vqubit.Qubit)
	assert isinstance(vqubit.NamedQubit("a"), vqubit.Qubit)
	# The qid variants are qids, but they are not qubits.
	assert isinstance(vqubit.LineQid(0), vqubit.Qid)
	assert not isinstance(vqubit.LineQid(0), vqubit.Qubit)
	assert not isinstance(vqubit.GridQid(0, 0), vqubit.Qubit)


def test_line_qubit_range_and_validation_edges() -> None:
	assert vqubit.LineQubit.range(0) == ()
	assert vqubit.LineQubit.range(-3) == ()
	assert vqubit.LineQubit(0) == vqubit.LineQubit(0)
	with pytest.raises(TypeError):
		vqubit.LineQubit(1.0)  # ty: ignore[invalid-argument-type]
	with pytest.raises(TypeError):
		vqubit.LineQubit("0")  # ty: ignore[invalid-argument-type]
	with pytest.raises(TypeError):
		vqubit.LineQubit(False)
	with pytest.raises(ValueError):
		vqubit.LineQubit(-1)


def test_line_qubit_ordering_and_non_qid_comparisons() -> None:
	assert vqubit.LineQubit(0) < vqubit.LineQubit(1)
	assert vqubit.LineQubit(0) <= vqubit.LineQubit(0)
	assert vqubit.LineQubit(2) > vqubit.LineQubit(0)
	assert vqubit.LineQubit(2) >= vqubit.LineQubit(2)
	with pytest.raises(TypeError):
		_ = vqubit.LineQubit(0) < "1"
	assert vqubit.LineQubit(0) != 0
	assert vqubit.NamedQubit("0") != vqubit.LineQubit(0)


def test_grid_qubit_rect_and_square_layouts() -> None:
	rect = vqubit.GridQubit.rect(rows=range(2), cols=range(3))
	assert rect == [
		vqubit.GridQubit(0, 0),
		vqubit.GridQubit(0, 1),
		vqubit.GridQubit(0, 2),
		vqubit.GridQubit(1, 0),
		vqubit.GridQubit(1, 1),
		vqubit.GridQubit(1, 2),
	]
	square = vqubit.GridQubit.square(2, top=1, left=2)
	assert set(square) == {
		vqubit.GridQubit(1, 2),
		vqubit.GridQubit(1, 3),
		vqubit.GridQubit(2, 2),
		vqubit.GridQubit(2, 3),
	}
	assert vqubit.GridQubit.rect(rows=[0], cols=[]) == []


def test_grid_qubit_neighbors_and_adjacency() -> None:
	q = vqubit.GridQubit(0, 0)
	# Without a grid, the four geometrically adjacent positions are returned.
	assert set(q.neighbors()) == {
		vqubit.GridQubit(-1, 0),
		vqubit.GridQubit(1, 0),
		vqubit.GridQubit(0, -1),
		vqubit.GridQubit(0, 1),
	}
	assert q.neighbors([]) == []
	assert set(q.neighbors(vqubit.GridQubit.square(2))) == {vqubit.GridQubit(1, 0), vqubit.GridQubit(0, 1)}
	assert q.is_adjacent(vqubit.GridQubit(0, 1))
	assert q.is_adjacent(vqubit.GridQubit(1, 0))
	assert not q.is_adjacent(vqubit.GridQubit(1, 1))
	assert not q.is_adjacent(q)


def test_grid_qubit_validates_integer_coordinates() -> None:
	with pytest.raises(TypeError):
		vqubit.GridQubit(0.5, 0)  # ty: ignore[invalid-argument-type]
	with pytest.raises(TypeError):
		vqubit.GridQubit(0, "1")  # ty: ignore[invalid-argument-type]
	with pytest.raises(TypeError):
		vqubit.GridQubit(True, 0)
	# Negative coordinates are allowed: grids need no origin.
	assert vqubit.GridQubit(-2, -3) == vqubit.GridQubit(-2, -3)


def test_named_qubit_names_and_validation() -> None:
	assert vqubit.NamedQubit("").name == ""
	assert str(vqubit.NamedQubit("")) == ""
	assert vqubit.NamedQubit("a") < vqubit.NamedQubit("b")
	assert vqubit.NamedQubit("alice") == vqubit.NamedQubit("alice")
	assert vqubit.NamedQubit("alice") != vqubit.NamedQubit("bob")
	with pytest.raises(TypeError):
		vqubit.NamedQubit(3)  # ty: ignore[invalid-argument-type]


def test_line_qid_and_grid_qid_variants() -> None:
	assert vqubit.LineQid.range(3) == (vqubit.LineQid(0), vqubit.LineQid(1), vqubit.LineQid(2))
	assert vqubit.LineQid(2).x == 2
	assert str(vqubit.LineQid(2)) == "2"
	assert repr(vqubit.LineQid(2)) == "vqubit.LineQid(2)"
	assert vqubit.LineQid(0) != vqubit.LineQubit(0)

	grid = vqubit.GridQid(1, 2)
	assert (grid.row, grid.col) == (1, 2)
	assert str(grid) == "(1, 2)"
	assert repr(grid) == "vqubit.GridQid(1, 2)"
	assert grid == vqubit.GridQid(1, 2)
	assert grid != vqubit.GridQubit(1, 2)

	with pytest.raises(TypeError):
		vqubit.LineQid(1.5)  # ty: ignore[invalid-argument-type]
	with pytest.raises(ValueError):
		vqubit.LineQid(-1)
	with pytest.raises(TypeError):
		vqubit.LineQid(True)
	with pytest.raises(TypeError):
		vqubit.GridQid(0.5, 1)  # ty: ignore[invalid-argument-type]


def test_named_qid_and_variable_qids() -> None:
	assert vqubit.NamedQid("a").name == "a"
	assert str(vqubit.NamedQid("a")) == "a"
	assert repr(vqubit.NamedQid("a")) == "vqubit.NamedQid('a')"
	assert vqubit.NamedQid("a") != vqubit.NamedQubit("a")
	with pytest.raises(TypeError):
		vqubit.NamedQid(3)  # ty: ignore[invalid-argument-type]

	parameter = vqubit.Parameter("i")
	variable = vqubit.VariableQid(parameter)
	assert variable.name == parameter
	assert str(variable) == "i"
	assert repr(variable) == "vqubit.VariableQid(vqubit.Parameter('i'))"
	assert variable == vqubit.VariableQid(vqubit.Parameter("i"))

	line = vqubit.VariableLineQid(parameter)
	assert line.x == parameter
	assert str(line) == "i"
	assert repr(line) == "vqubit.VariableLineQid(vqubit.Parameter('i'))"

	grid = vqubit.VariableGridQid(parameter, 2)
	assert grid.row == parameter
	assert grid.col == 2
	assert str(grid) == "(i, 2)"
	assert repr(grid) == "vqubit.VariableGridQid(vqubit.Parameter('i'), 2)"


def test_hashes_and_sets_are_type_specific() -> None:
	values = [
		vqubit.LineQubit(0),
		vqubit.GridQubit(0, 0),
		vqubit.NamedQubit("0"),
		vqubit.LineQid(0),
		vqubit.GridQid(0, 0),
		vqubit.NamedQid("0"),
	]
	assert len(set(values)) == len(values)
	assert len({vqubit.LineQubit(0), vqubit.LineQubit(0)}) == 1
	assert hash(vqubit.GridQubit(1, 2)) == hash(vqubit.GridQubit(1, 2))
	assert hash(vqubit.NamedQubit("a")) == hash(vqubit.NamedQubit("a"))


def test_qubit_order_helpers() -> None:
	qubits = [vqubit.LineQubit(2), vqubit.LineQubit(0), vqubit.LineQubit(1)]
	assert vqubit.QubitOrder.DEFAULT.order_for(qubits) == (
		vqubit.LineQubit(0),
		vqubit.LineQubit(1),
		vqubit.LineQubit(2),
	)
	assert vqubit.QubitOrder.explicit([vqubit.LineQubit(2)]).order_for(qubits) == (
		vqubit.LineQubit(2),
		vqubit.LineQubit(0),
		vqubit.LineQubit(1),
	)
	named = [vqubit.NamedQubit("b"), vqubit.NamedQubit("a")]
	assert vqubit.QubitOrder.sorted_by(str).order_for(named) == (vqubit.NamedQubit("a"), vqubit.NamedQubit("b"))
	assert vqubit.QubitOrder.DEFAULT(qubits) == vqubit.QubitOrder.DEFAULT.order_for(qubits)
