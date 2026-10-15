# Benchmarks

Esta página compara vqubit y [Cirq](https://quantumai.google/cirq) en las mismas
tareas, para que veas el equilibrio entre velocidad y memoria de una librería
sin dependencias y en Python puro, frente a una basada en NumPy.

El harness está en `benchmarks/compare.py` y solo usa la biblioteca estándar.
Ejecuta cada tarea en un **subproceso nuevo** (para que la memoria pico no se
contamine con tareas anteriores) e informa de:

- **tiempo**: el mejor tiempo de pared de varias ejecuciones;
- **RAM extra**: la memoria residente pico (`VmHWM`) durante la tarea menos la
  línea base tras importar el paquete.

## Cómo ejecutarlo

```bash
# Instala la dependencia de comparación (Cirq) en tu entorno.
uv pip install --group bench

# Imprime la tabla (la columna de Cirq queda vacía si no está instalado).
python benchmarks/compare.py

# O los números en crudo como JSON.
python benchmarks/compare.py --json
```

## Metodología

Cada tarea está implementada dos veces —una con vqubit y otra con Cirq— haciendo
el mismo trabajo. Cada tarea se ejecuta en su propio proceso:

1. importar el paquete y ejecutar la tarea una vez para calentar;
2. registrar la RSS de línea base;
3. ejecutar la tarea `N` veces y quedarse con el mejor tiempo;
4. leer `VmHWM` (RSS pico) para la cifra de memoria.

Los tiempos son de reloj de pared en un solo núcleo. La máquina usada para los
números de abajo es un equipo de desarrollo Linux x86-64 normal; **ejecuta el
harness en tu propia máquina** para obtener cifras que te importen.

## Resultados

Ejecución representativa (Python 3.15, Cirq 1.7.0, NumPy 2.x). La "RAM" es la
memoria residente pico durante la tarea (que, con la importación aislada, ya
incluye el coste de importar el paquete).

| Tarea | vqubit | Cirq | Aceleración | RAM vqubit | RAM Cirq |
|---|---:|---:|---:|---:|---:|
| Importar el paquete | — | — | — | 77 MB | 208 MB |
| Construir un circuito de 5 qubits y 300 momentos | 19 ms | 19 ms | ≈1× | 78 MB | 209 MB |
| Unitario de un circuito de 6 qubits | 0.5 ms | 0.8 ms | 1.5× | 77 MB | 208 MB |
| Simular vector de estado de 12 qubits | 0.3 ms | 1.9 ms | 6.3× | 77 MB | 208 MB |
| Muestrear 1000 disparos de 12 qubits | 2.7 ms | 6.7 ms | 2.5× | 82 MB | 209 MB |
| Simular matriz de densidad de 7 qubits | 3.2 ms | 7.8 ms | 2.5× | 78 MB | 209 MB |
| Simular vector de estado de 20 qubits | 41 ms | 42 ms | ≈1× | 93 MB | 246 MB |
| Simular matriz de densidad de 10 qubits | 293 ms | 370 ms | 1.3× | 109 MB | 258 MB |
| Descomposición KAK de un unitario 4×4 | 0.26 ms | 0.48 ms | 1.8× | 77 MB | 209 MB |
| Matriz de 1 qubit a PhasedXZ | <0.1 ms | <0.1 ms | 5.6× | 77 MB | 208 MB |
| Multiplicar cadenas de Pauli 1000× | 5.7 ms | 9.5 ms | 1.7× | 77 MB | 208 MB |
| Exportar OpenQASM de 50 puertas | 0.6 ms | 1.4 ms | 2.4× | 77 MB | 208 MB |
| Vector de Bloch de 12 qubits | 0.3 ms | 2.0 ms | 7.7× | 77 MB | 208 MB |

Una **aceleración por encima de 1×** significa que vqubit es más rápido; por
debajo de 1× significa que Cirq es más rápido.

## Cómo interpretarlo

- **vqubit gana o empata en todas las tareas**, y además importa unas **2.7×
  menos memoria** (77 MB frente a 208 MB) y usa menos memoria pico en las
  simulaciones grandes. Su wheel no tiene dependencias de ejecución, así que
  importarlo es mucho más barato que arrastrar NumPy, SciPy, SymPy, pandas y
  matplotlib.
- **Las dos tareas en empate están limitadas por el ancho de banda de memoria.**
  Construir un circuito y simular 20 qubits mueven decenas de megabytes por paso;
  una vez vectorizada la aritmética (lo está en ambos), el tiempo lo fija la
  velocidad a la que fluyen las amplitudes por memoria, así que las dos
  implementaciones convergen a los mismos ~40 ms. vqubit lo logra con los kernels
  nativos de V y vistas NumPy sin copia de su buffer de estado; Cirq lo logra con
  arrays de NumPy.
- **Las simulaciones pequeñas y medianas favorecen a vqubit.** El vector de
  estado y la matriz de densidad permanecen en buffers nativos y se devuelven como
  vistas NumPy sin copia, así que no hay ida y vuelta a Python por puerta; el
  unitario del circuito, el muestreo, el vector de Bloch y las descomposiciones
  también son más rápidos.
- Cuando NumPy está instalado, la API pública devuelve `numpy.ndarray`
  (ver [Instalación](../installation.md)); sin él todo sigue funcionando con
  listas, solo que sin las rutas aceleradas por NumPy.

## Reproducir los números

```bash
uv pip install --group bench
python benchmarks/compare.py --json > benchmarks.json
```

El grupo de dependencias `bench` es solo para la comparación; el wheel publicado
de vqubit nunca depende de Cirq ni de NumPy.
