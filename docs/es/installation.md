# Instalación

## Requisitos

- **Python 3.11 o superior.** vqubit funciona con CPython 3.11, 3.12, 3.13,
  3.14 y 3.15.
- **Nada más.** vqubit no tiene dependencias: no instala numpy, ni sympy, ni
  ningún otro paquete.

## Instalar desde PyPI

La forma recomendada es con `pip`:

```console
$ pip install vqubit
```

Eso es todo. Para comprobar que funciona:

```console
$ python -c "import vqubit; print(vqubit.__version__)"
{{ vqubit_version }}
```

!!! note "Sin dependencias"
    Instalar vqubit **no** descarga ningún otro paquete. La parte pesada está
    compilada dentro de la extensión nativa, así que no hace falta ni siquiera un
    compilador para usarla.

## Instalar desde el código fuente

Si quieres compilar vqubit tú mismo (por ejemplo, para contribuir), necesitas:

- el compilador de **V**,
- un compilador de **C** (gcc o clang),
- [**vcraft**](https://github.com/rroblf01/vcraft), la herramienta de
  construcción,
- [**uv**](https://docs.astral.sh/uv/) para las herramientas de desarrollo.

```console
$ git clone https://github.com/rroblf01/vqubit
$ cd vqubit
$ python -m venv .venv && source .venv/bin/activate
$ pip install vcraft
$ uv pip install --group dev
$ vcraft develop
```

`vcraft develop` construye la extensión nativa y la instala en el entorno
virtual activo.

## Plataformas

| Sistema operativo | Arquitectura   | Estado          |
| ----------------- | -------------- | --------------- |
| Linux             | x86-64         | Soportado       |
| Linux             | ARM64 (aarch64) | Soportado      |
| macOS             | Intel (x86-64) | Desde fuente    |
| macOS             | Apple Silicon (ARM64) | Soportado |
| Windows           | x86-64         | Soportado       |

macOS en Intel (x86-64) no se compila; cualquier plataforma sin rueda puede
compilar vqubit desde el código fuente como arriba.

## Builds más rápidos

Las ruedas publicadas apuntan a un x86-64 genérico, para que corran en cualquier
máquina. Para exprimir más velocidad en una compilación tuya, añade flags de
compilador a `vcraft.toml`:

```toml
[build]
cflags = "-march=native"
```

`-march=native` (o una base fija como `-march=x86-64-v3`) deja que el compilador
vectorice con las instrucciones SIMD de tu CPU, lo que hace más rápidos los
kernels nativos siempre que el estado quepa en caché. La rueda resultante queda
atada a ese nivel de CPU, así que es para una compilación local o propia, no
para publicar.
