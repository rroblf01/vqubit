"""Fuzzers for the three importer parsers: QASM, Quirk JSON and vqubit JSON.

A parser is only as good as its behaviour on garbage: it may return a circuit or
raise the documented error, but it must never leak an ``AttributeError``,
``KeyError`` or ``IndexError``, and it must terminate. These fuzzers mutate a
valid document with a fixed seed, so a failure is reproducible.

The seeds are deliberate: the suite must be deterministic in CI.
"""

from __future__ import annotations

import base64
import contextlib
import json
import random

import pytest

import vqubit

Q0, Q1 = vqubit.LineQubit.range(2)
CIRCUIT = vqubit.Circuit(vqubit.H(Q0), vqubit.CNOT(Q0, Q1), vqubit.measure(Q0, key="m"))
NOISE = "abcXYZ[]{}();,=@^0123456789 \n\"'"


def mutations(text: str, *, seed: int, count: int = 200) -> list[str]:
	"""Random single- and multi-character edits of ``text``."""
	random_source = random.Random(seed)
	result: list[str] = []
	for _ in range(count):
		characters = list(text)
		for _ in range(random_source.randint(1, 4)):
			if not characters:
				break
			position = random_source.randrange(len(characters))
			operation = random_source.choice(["replace", "delete", "insert", "truncate"])
			if operation == "replace":
				characters[position] = random_source.choice(NOISE)
			elif operation == "delete":
				del characters[position]
			elif operation == "insert":
				characters.insert(position, random_source.choice(NOISE))
			else:
				del characters[position:]
		result.append("".join(characters))
	return result


def quirk_json(circuit: vqubit.Circuit) -> str:
	"""The JSON payload behind ``quirk_url``."""
	payload = vqubit.quirk_url(circuit).split("circuit=", 1)[1]
	payload += "=" * (-len(payload) % 4)
	return base64.urlsafe_b64decode(payload).decode()


def test_fuzzed_qasm_only_raises_qasm_errors() -> None:
	valid = vqubit.qasm(CIRCUIT)
	parsed = 0
	for text in mutations(valid, seed=1):
		try:
			circuit = vqubit.circuit_from_qasm(text)
		except vqubit.QasmParseError:
			continue
		assert isinstance(circuit, vqubit.Circuit)
		parsed += 1
	# The mutations must not all be rejected: the fuzzer has to exercise the
	# parser, not just its error path.
	assert parsed > 0


def test_fuzzed_quirk_json_only_raises_value_errors() -> None:
	valid = quirk_json(CIRCUIT)
	parsed = 0
	for text in mutations(valid, seed=2):
		try:
			circuit = vqubit.quirk_json_to_circuit(text)
		except ValueError:  # includes json.JSONDecodeError
			continue
		assert isinstance(circuit, vqubit.Circuit)
		parsed += 1
	assert parsed > 0


def test_fuzzed_vqubit_json_only_raises_value_errors() -> None:
	valid = vqubit.to_json(CIRCUIT)
	parsed = 0
	for text in mutations(valid, seed=3):
		try:
			value = vqubit.read_json(text)
		except ValueError:  # includes json.JSONDecodeError
			continue
		parsed += 1
		assert value is not None
	assert parsed > 0


def test_fuzzed_quirk_url_only_raises_value_errors() -> None:
	valid = vqubit.quirk_url(CIRCUIT)
	for text in mutations(valid, seed=4, count=100):
		try:
			circuit = vqubit.quirk_url_to_circuit(text)
		except ValueError:
			continue
		assert isinstance(circuit, vqubit.Circuit)


@pytest.mark.parametrize("seed", range(5))
def test_random_json_objects_are_rejected_cleanly(seed: int) -> None:
	"""Arbitrary JSON shapes must not crash the readers either."""
	random_source = random.Random(seed)

	def random_value(depth: int = 0) -> object:
		choice = random_source.choice(["int", "float", "str", "bool", "null", "list", "dict"])
		if depth > 3:
			choice = random_source.choice(["int", "str", "bool", "null"])
		if choice == "int":
			return random_source.randint(-5, 5)
		if choice == "float":
			return random_source.random()
		if choice == "str":
			return "".join(random_source.choice("abc[]{}") for _ in range(random_source.randint(0, 4)))
		if choice == "bool":
			return random_source.choice([True, False])
		if choice == "null":
			return None
		if choice == "list":
			return [random_value(depth + 1) for _ in range(random_source.randint(0, 3))]
		return {random_source.choice(["cirq_type", "cols", "qubits", "x"]): random_value(depth + 1)}

	for _ in range(40):
		document = random_value()
		text = json.dumps(document)
		# Both readers may accept a document (`read_json` decodes any value, and
		# an empty `cols` is an empty circuit); what they may not do is leak an
		# unexpected exception type.
		for parser in (vqubit.read_json, vqubit.quirk_json_to_circuit):
			with contextlib.suppress(ValueError):
				parser(text)
