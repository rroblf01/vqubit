# vqubit

**A dependency-free quantum computing framework for Python, written in [V](https://vlang.io).**

vqubit lets you build quantum circuits and run them on local simulators with
**zero third-party dependencies**. The numerically heavy work — linear algebra,
gate matrices and simulation — is written in V and compiled into a native
extension; the public API is dependency-free Python on top of it.

The API mirrors [Cirq](https://quantumai.google/cirq) closely, so if you know
Cirq you already know vqubit.

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.measure(q0, q1, key="result"))
print(circuit)

result = vq.Simulator().run(circuit, repetitions=100)
print(result.histogram(key="result"))
```

## Highlights

- **No dependencies.** `pip install vqubit` installs nothing else. Everything
  runs offline, on your machine.
- **Native speed.** Linear algebra, gate matrices and simulation are compiled
  from V, not interpreted.
- **One wheel, every Python.** A single abi3 wheel covers CPython 3.11–3.15 per
  platform.
- **Cross-platform.** Linux and macOS (x86-64 and ARM).
- **Cirq-compatible API.** Build, transform and simulate circuits the way Cirq
  documents them.

## Installation

```console
$ pip install vqubit
```

Requires Python 3.11 or newer. No other packages are installed.

## Documentation

The full documentation is published at **https://rroblf01.github.io/vqubit/**
in both [Spanish](https://rroblf01.github.io/vqubit/es/) and
[English](https://rroblf01.github.io/vqubit/en/). It starts from "what is a
qubit?" and builds up to circuits and simulation.

## Status

vqubit is **alpha**. The core object model (qubits, gates, operations, moments
and circuits), the protocols and the local simulators are being built out
phase by phase. See [CHANGELOG.md](CHANGELOG.md) for what is in place.

Windows support is currently blocked by the build tool [vcraft](https://github.com/rroblf01/vcraft),
which does not yet produce Windows wheels.

## Development

vqubit is built with [vcraft](https://github.com/rroblf01/vcraft). You need the
V compiler, a C compiler, and [uv](https://docs.astral.sh/uv/) for the
development tooling.

```console
$ python -m venv .venv && source .venv/bin/activate
$ pip install vcraft
$ uv pip install --group dev   # pytest, ruff, ty, build
$ vcraft develop               # build and install into the active virtualenv
$ python main.py               # run the example
$ pytest                       # run the test suite
$ vcraft build --release       # write an abi3 wheel into dist/
```

The dev and docs tooling are `[dependency-groups]` in `pyproject.toml`; use
`uv pip install --group docs` for the documentation toolchain. (`uv pip install`
adds to the virtualenv without pruning it — `uv sync` would remove the editable
extension that `vcraft develop` installs.)

The development setup is the only thing that needs the V toolchain — the
published wheels need nothing.

## License

MIT. See [LICENSE](LICENSE). © Ricardo Robles Fernández.
