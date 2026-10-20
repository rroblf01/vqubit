# Security policy

vqubit is a quantum-computing library: it does not execute untrusted code by
itself, but it does ship a **native extension** (compiled from V) and three
**parsers** (OpenQASM, Quirk, JSON). A bug in either can crash the interpreter or
use more memory than a program asked for, so both are in scope for this policy.

## Supported versions

| Version | Supported |
| ------- | --------- |
| 0.3.x   | ✅        |
| < 0.3   | ❌        |

The project is pre-1.0: security fixes land in a patch release of the latest
minor version. Once 1.0 is out, the latest minor gets fixes.

## Reporting a vulnerability

**Please do not open a public issue for a security problem.** Use either:

- GitHub's private vulnerability reporting: *Security → Report a vulnerability*
  on <https://github.com/rroblf01/vqubit>; or
- email to **ricardo.robles.fernandez@gmail.com** with the details.

A good report includes:

- what you were running (a minimal script is ideal) and the version
  (`python -c "import vqubit; print(vqubit.__version__)"`);
- the platform, the Python version and whether the wheel or a local build was in
  use;
- what happened (crash, memory growth, wrong result) and what you expected;
- if it is a parser problem, the exact input.

## What to expect

| Step | Target |
| ---- | ------ |
| Acknowledge the report | within 7 days |
| First assessment and severity | within 14 days |
| Fix or mitigation in a release | within 90 days |

You will be credited in the release notes unless you ask otherwise.

## Scope

In scope:

- **The native kernels** (`vqubit_native`): out-of-bounds access, integer
  overflow in the shape/size arithmetic, use-after-free of the state buffers.
  The extension is built with bounds checking on and the release wheels are
  reproducible from the tagged commit.
- **The parsers** (`circuit_from_qasm`, `quirk_json_to_circuit`, `read_json`):
  inputs that make them hang, grow memory without bound, or raise an exception
  other than the documented `QasmParseError` / `ValueError`.
- **The build and release pipeline**: a wheel that is not built from the tag, a
  dependency added at runtime, a signature that does not match the published
  hashes.

Out of scope:

- Denial of service from an intentionally huge but well-formed circuit (the
  library is a simulator; ask it to simulate a 30-qubit state vector and it will
  try);
- vulnerabilities in optional third-party packages (`numpy`, `sympy`,
  `matplotlib`, `cirq`) when the user installs the extras;
- the example URLs and links in the documentation.

## Hardening already in place

- No runtime dependencies: the wheel needs only the standard library and the
  compiled extension, so there is no supply chain beyond the build itself.
- The parsers are fuzzed in CI (`tests/test_fuzz_parsers.py`), and malformed
  input is required to raise `QasmParseError` or `ValueError`.
- The test matrix runs on three operating systems and five Python versions, and
  the numeric harness (`benchmarks/verify.py`) checks the kernels against closed
  forms.
