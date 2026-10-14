"""The Bell state: the smallest example of quantum entanglement.

Run it after building the extension::

    python examples/bell_state.py

It prints the circuit, the entangled state it produces and a histogram of
measurement outcomes. Only ``00`` and ``11`` ever appear, never ``01`` or
``10``, because the two qubits are entangled.
"""

from __future__ import annotations

import vqubit as vq


def main() -> None:
	q0, q1 = vq.LineQubit.range(2)

	circuit = vq.Circuit(
		vq.H(q0),
		vq.CNOT(q0, q1),
	)
	print("Circuit:")
	print(circuit)
	print()

	state = vq.Simulator().simulate(circuit)
	print("State after the circuit:")
	print(" ", state.dirac_notation())
	print()

	circuit.append(vq.measure(q0, q1, key="result"))
	result = vq.Simulator(seed=0).run(circuit, repetitions=1000)
	print("Measurement histogram over 1000 shots:")
	print(" ", result.histogram(key="result"))


if __name__ == "__main__":
	main()
