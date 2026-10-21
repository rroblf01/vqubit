# Installation

## Requirements

- **Python 3.11 or newer.** vqubit works on CPython 3.11, 3.12, 3.13, 3.14 and
  3.15.
- **Nothing else.** vqubit has no dependencies: it does not install numpy,
  sympy, or any other package.

## Install from PyPI

The recommended way is with `pip`:

```console
$ pip install vqubit
```

That is all. To check that it works:

```console
$ python -c "import vqubit; print(vqubit.__version__)"
{{ vqubit_version }}
```

!!! note "No dependencies"
    Installing vqubit does **not** download any other package. The heavy
    lifting is compiled into the native extension, so you do not even need a
    compiler to use it.

## Optional integrations

vqubit **declares no extras** (`vqubit[vis]`, `vqubit[interop]`): the optional
pieces are installed separately, under names you already know, and vqubit
notices them at run time. It is the same story as Cirq, which publishes no
extras either.

```console
$ pip install matplotlib   # vqubit.vis
$ pip install cirq-core    # to_cirq / from_cirq
$ pip install sympy        # ParamResolver keys that are sympy.Symbols
```

Everything else keeps working without them; the functions that need one say so
when called, with the missing `pip install` in the message. The `interop` and
`vis` groups in `pyproject.toml` are for developing vqubit and are not part of
the published API.

## Install from source

If you want to build vqubit yourself (for example, to contribute), you need:

- the **V** compiler,
- a **C** compiler (gcc or clang),
- [**vcraft**](https://github.com/rroblf01/vcraft), the build tool,
- [**uv**](https://docs.astral.sh/uv/) for the development tooling.

```console
$ git clone https://github.com/rroblf01/vqubit
$ cd vqubit
$ python -m venv .venv && source .venv/bin/activate
$ pip install vcraft
$ uv pip install --group dev
$ vcraft develop
```

`vcraft develop` builds the native extension and installs it into the active
virtual environment.

## Platforms

| Operating system | Architecture   | Status        |
| ---------------- | -------------- | ------------- |
| Linux            | x86-64         | Supported     |
| Linux            | ARM64 (aarch64) | Supported    |
| macOS            | Intel (x86-64) | From source   |
| macOS            | Apple Silicon (ARM64) | Supported |
| Windows          | x86-64         | Supported     |

macOS on Intel (x86-64) is not built; any platform without a wheel can build
vqubit from source as above.

## Faster builds

The published wheels target a generic x86-64 so they run on any machine. To
squeeze more speed out of a build you make yourself, add compiler flags to
`vcraft.toml`:

```toml
[build]
cflags = "-march=native"
```

`-march=native` (or a fixed baseline such as `-march=x86-64-v3`) lets the
compiler vectorise with the CPU's SIMD instructions, which makes the native
kernels faster whenever the state fits in cache. The resulting wheel is tied to
that CPU level, so this is for a local or self-hosted build, not for a wheel you
publish.
