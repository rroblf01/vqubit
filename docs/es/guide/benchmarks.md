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

Ejecución representativa (Python 3.15, Cirq 1.7.0, NumPy 2.x):

| Tarea | vqubit | Cirq | Aceleración | RAM extra vqubit | RAM extra Cirq |
|---|---:|---:|---:|---:|---:|
| Construir un circuito de 5 qubits y 300 momentos | 18.4 ms | 18.7 ms | 1.02× | 0 MB | 0 MB |
| Unitario de un circuito de 6 qubits | 45.7 ms | 0.7 ms | 0.02× | 1 MB | 0 MB |
| Simular vector de estado de 12 qubits | 4.0 ms | 1.7 ms | 0.41× | 0 MB | 0 MB |
| Muestrear 1000 disparos de 12 qubits | 96.7 ms | 6.5 ms | 0.07× | 0 MB | 0 MB |
| Simular matriz de densidad de 7 qubits | 22.7 ms | 8.3 ms | 0.37× | 7 MB | 0 MB |
| Simular vector de estado de 20 qubits | 1453 ms | 42 ms | 0.03× | 199 MB | 39 MB |
| Simular matriz de densidad de 10 qubits | 2366 ms | 387 ms | 0.16× | 593 MB | 49 MB |
| Descomposición KAK de un unitario 4×4 | 0.3 ms | 0.9 ms | 3.2× | 0 MB | 0 MB |
| Matriz de 1 qubit a PhasedXZ | <0.1 ms | <0.1 ms | 5.6× | 0 MB | 0 MB |
| Multiplicar cadenas de Pauli 1000× | 6.2 ms | 9.7 ms | 1.6× | 0 MB | 0 MB |
| Exportar OpenQASM de 50 puertas | 0.5 ms | 1.3 ms | 2.5× | 0 MB | 0 MB |
| Vector de Bloch de 12 qubits | 8.2 ms | 1.7 ms | 0.21× | 0 MB | 0 MB |

Una **aceleración por encima de 1×** significa que vqubit es más rápido; por
debajo de 1× significa que Cirq es más rápido.

## Cómo interpretarlo

- **Construcción de circuitos, descomposiciones, álgebra de Pauli y QASM** son
  donde vqubit es competitivo o más rápido. Son operaciones en Python puro con
  pocas reservas de memoria, y las representaciones compactas de vqubit salen a
  cuenta.
- **Las simulaciones grandes y el unitario de un circuito completo** favorecen a
  Cirq. Los simuladores de vqubit hacen el álgebra lineal en los kernels nativos
  de V, pero el resultado se materializa como objetos `complex` de Python, lo
  que cuesta tiempo y memoria a 20 qubits. Cirq lo mantiene todo en arrays de
  NumPy, así que sus conversiones son baratas. Cuando NumPy está instalado, los
  accesores `unitary`/`kraus`/estado de vqubit devuelven `numpy.ndarray` (ver
  [Instalación](../installation.md)), lo que reduce la brecha para quien los
  consume, pero no cambia el coste interno del simulador.
- **La memoria** está dominada por la importación en tareas pequeñas (Cirq
  arrastra NumPy, SciPy, SymPy, pandas y matplotlib) y por el búfer de
  simulación en tareas grandes. El wheel de vqubit no tiene dependencias de
  ejecución, así que importarlo es mucho más barato; el harness informa de la
  memoria *extra* durante la tarea para que el coste de importación no oculte el
  de la simulación.

## Reproducir los números

```bash
uv pip install --group bench
python benchmarks/compare.py --json > benchmarks.json
```

El grupo de dependencias `bench` es solo para la comparación; el wheel publicado
de vqubit nunca depende de Cirq ni de NumPy.
