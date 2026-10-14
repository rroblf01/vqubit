"""Grover's search over two qubits, run on the local simulator.

    python examples/grover.py

One Grover iteration finds the marked state ``|11>`` with certainty. The script
prints the circuit, then a histogram over 1000 shots; ``3`` is the integer for
``11``.
"""

from __future__ import annotations

import vqubit as vq


def main() -> None:
	q0, q1 = vq.LineQubit.range(2)

	circuit = vq.Circuit(
		# Start in an equal superposition of all four states.
		vq.H(q0),
		vq.H(q1),
		# Oracle: flip the phase of |11>.
		vq.CZ(q0, q1),
		# Diffusion: reflect about the average.
		vq.H(q0),
		vq.H(q1),
		vq.X(q0),
		vq.X(q1),
		vq.CZ(q0, q1),
		vq.X(q0),
		vq.X(q1),
		vq.H(q0),
		vq.H(q1),
	)

	print("Circuit:")
	print(circuit)
	print()

	circuit.append(vq.measure(q0, q1, key="result"))
	result = vq.Simulator(seed=0).run(circuit, repetitions=1000)
	print("Measurement histogram:")
	print(" ", result.histogram(key="result"))


if __name__ == "__main__":
	main()
