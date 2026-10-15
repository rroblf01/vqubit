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

Ejecución representativa (Python 3.15, Cirq 1.7.0, NumPy 2.x). Las columnas de
RAM son la memoria residente pico durante la tarea, que incluye importar el
paquete (cada tarea corre en un subproceso nuevo que importa solo el paquete bajo
prueba).

| Tarea | vqubit | Cirq | Aceleración | RAM vqubit | RAM Cirq |
|---|---:|---:|---:|---:|---:|
| Importar el paquete | — | — | — | 77 MB | 208 MB |
| Construir un circuito de 5 qubits y 300 momentos | 10 ms | 19 ms | 1.9× | 78 MB | 209 MB |
| Unitario de un circuito de 6 qubits | 0.4 ms | 0.8 ms | 2.1× | 77 MB | 208 MB |
| Simular vector de estado de 12 qubits | 0.2 ms | 2.5 ms | 12.2× | 77 MB | 209 MB |
| Muestrear 1000 disparos de 12 qubits | 3.6 ms | 6.6 ms | 1.8× | 82 MB | 209 MB |
| Simular matriz de densidad de 7 qubits | 2.1 ms | 7.3 ms | 3.5× | 78 MB | 209 MB |
| Simular vector de estado de 20 qubits | ≈40 ms | ≈42 ms | ≈1× | 93 MB | 247 MB |
| Simular matriz de densidad de 10 qubits | 162 ms | 360 ms | 2.2× | 125 MB | 257 MB |
| Simular un circuito de 250 capas y 12 qubits | 31 ms | 193 ms | 6.3× | 82 MB | 212 MB |
| Descomposición KAK de un unitario 4×4 | 0.25 ms | 0.45 ms | 1.8× | 77 MB | 209 MB |
| Matriz de 1 qubit a PhasedXZ | <0.1 ms | <0.1 ms | 5.0× | 77 MB | 208 MB |
| Multiplicar cadenas de Pauli 1000× | 3.6 ms | 8.7 ms | 2.4× | 77 MB | 208 MB |
| Exportar OpenQASM de 50 puertas | 0.4 ms | 1.6 ms | 4.3× | 77 MB | 208 MB |
| Vector de Bloch de 12 qubits | 0.2 ms | 1.8 ms | 9.6× | 77 MB | 208 MB |

Una **aceleración por encima de 1×** significa que vqubit es más rápido; por
debajo de 1× significa que Cirq es más rápido.

## Cómo interpretarlo

- **vqubit gana en todas las tareas**, y además importa unas **2.7× menos
  memoria** (77 MB frente a 208 MB) y usa menos memoria pico en las simulaciones
  grandes. Su wheel no tiene dependencias de ejecución, así que importarlo es
  mucho más barato que arrastrar NumPy, SciPy, SymPy, pandas y matplotlib.
- **Los circuitos con muchas puertas son donde la diferencia es mayor.** Un
  circuito de 250 capas y 12 qubits corre ~6× más rápido que Cirq: vqubit
  construye la matriz de cada puerta de 1 qubit en forma cerrada y la mantiene
  (junto con su forma plana nativa) memoizada, así que el coste Python por puerta
  que suele dominar estas cargas casi desaparece.
- **Simular 20 qubits está limitado por el ancho de banda de memoria.** Mueve
  decenas de megabytes por puerta; una vez vectorizada la aritmética (lo está en
  ambos), el tiempo lo fija la velocidad a la que fluyen las amplitudes por
  memoria, así que las dos implementaciones convergen a los mismos ~40 ms. vqubit
  lo logra con los kernels nativos de V y vistas NumPy sin copia de su buffer de
  estado; Cirq lo logra con arrays de NumPy.
- **La simulación de matrices de densidad es ~2× más rápida.** Los kernels de
  canal escalan la diagonal (CZ/CCZ) directamente y factorizan el producto
  general por bloques en dos productos de matrices.
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
