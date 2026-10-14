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
0.1.0
```

!!! note "No dependencies"
    Installing vqubit does **not** download any other package. The heavy
    lifting is compiled into the native extension, so you do not even need a
    compiler to use it.

## Install from source

If you want to build vqubit yourself (for example, to contribute), you need:

- the **V** compiler,
- a **C** compiler (gcc or clang),
- [**vcraft**](https://github.com/rroblf01/vcraft), the build tool.

```console
$ git clone https://github.com/rroblf01/vqubit
$ cd vqubit
$ python -m venv .venv && source .venv/bin/activate
$ pip install vcraft
$ pip install -r requirements-dev.txt
$ vcraft develop
```

`vcraft develop` builds the native extension and installs it into the active
virtual environment.

## Platforms

| Operating system | Architecture   | Status        |
| ---------------- | -------------- | ------------- |
| Linux            | x86-64         | Supported     |
| Linux            | ARM64 (aarch64) | Supported    |
| macOS            | Intel (x86-64) | Supported     |
| macOS            | Apple Silicon (ARM64) | Supported |
| Windows          | x86-64         | Pending (blocked by vcraft) |

!!! warning "Windows"
    Windows is not supported yet because the build tool
    [vcraft](https://github.com/rroblf01/vcraft) does not produce Windows wheels
    at this time. It is on its roadmap.
