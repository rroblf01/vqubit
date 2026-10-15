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
0.1.0
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
| macOS             | Intel (x86-64) | Soportado       |
| macOS             | Apple Silicon (ARM64) | Soportado |
| Windows           | x86-64         | Pendiente (bloqueado por vcraft) |

!!! warning "Windows"
    Windows aún no está soportado porque la herramienta de construcción
    [vcraft](https://github.com/rroblf01/vcraft) todavía no genera ruedas para
    Windows. Está en su hoja de ruta.
