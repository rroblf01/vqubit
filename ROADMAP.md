# vqubit roadmap

[Español](#roadmap-hacia-la-100) · [English](#roadmap-to-100)

---

## Roadmap hacia la 1.0.0

Esta hoja de ruta traza el camino desde la **0.1.0** (primera publicación en
PyPI) hasta la **1.0.0**, el momento en el que la API pública queda congelada y
el proyecto ofrece garantías de estabilidad. Es un documento vivo: se actualiza
a medida que avanzan los hitos.

> **Actualización (0.2.1, 2026-10-19).** Auditoría de paridad contra Cirq 1.7.0.
> Se añade un **arnés de conformidad** (`benchmarks/parity.py`: 1733
> comprobaciones) que compara el *contrato* con Cirq —nombres, unitarias, `**`,
> `isinstance`, clases— y no solo los números; las divergencias deliberadas viven
> en una lista con motivo que **no puede pudrirse** (una entrada que ya no tapa
> nada hace fallar el CI). La auditoría corrige cinco casillas ya cumplidas, mide
> la cobertura (77 %), y destapa el hueco de paridad más visible que queda:
> **`**` sobre las puertas planas**. Las peticiones a vcraft están en
> [`UPSTREAM.md`](UPSTREAM.md).
>
> **Actualización (0.3, 2026-10-19).** Las puertas planas **son** las puertas de
> potencia: `vq.CZ` es un `CZPowGate(exponent=1)`, así que
> `isinstance(vq.CZ, vq.CZPowGate)` y `vq.CZ == vq.CZPowGate()` son `True`, como
> en Cirq, y `vq.CZ ** 0.5` dibuja `@` sobre `@^0.5`, el diagrama de
> `cirq.CZ ** 0.5`. Las **12** entradas con `plan = "0.3"` de la lista de
> toleradas quedan tapadas y desaparecen, el arnés gana una sección de
> **diagramas** (1789 comprobaciones, 4 toleradas) y la auditoría destapa dos
> huecos nuevos, anotados en el eje de paridad: la **igualdad periódica** de Cirq
> y su *gateset* de **átomos neutros**.
>
> **Actualización (0.3, 2026-10-19, segunda).** Cerrados esos dos huecos: la
> igualdad de valor es **periódica** como la de Cirq (con `PeriodicValue`; de paso
> salieron la **fase global** que `global_shift` debía ser y `vq.CSWAP ** ±1`, que
> no existían) y el *gateset* de **átomos neutros** es el de Cirq, con las
> familias de puertas aceptando un **tipo**. El arnés gana la sección `equality`
> (78 comprobaciones) y pasa a **1867** comprobaciones con **4** divergencias
> toleradas; ese trabajo destapa a su vez la igualdad de las puertas de
> **Clifford** (comparan matriz, no *tableau*) y el **exponente canónico** de
> `PhasedXZGate`, los dos anotados en el eje de paridad con su criterio.
>
> **Actualización (0.3, 2026-10-19, tercera).** Los simuladores hablan como los de
> Cirq: `run_sweep`, `run_sweep_iter`, `simulate_sweep`, `simulate_sweep_iter`,
> `simulate_moment_steps`, `run_batch` y `compute_amplitudes`, con
> `param_resolver` en el segundo argumento (`simulator.simulate(circuito,
> resolver)`) y los resultados llevando `params`, `qubit_map`, `measurements` y
> `final_state_vector`. La sección `sweeps` de `benchmarks/verify.py` lo compara
> con Cirq (601 comprobaciones en total) y el trabajo destapó y arregló tres
> cosas: resolver un parámetro **a cero** reventaba con `ZeroDivisionError`,
> `has_stabilizer_effect` **lanzaba** en vez de devolver `False`, y el simulador
> de Clifford no aceptaba potencias enteras (`X ** 3`). Lo que queda de este eje
> (valores esperados, `sample` como `DataFrame`, los barridos con la firma de
> Cirq y el `simulate_moment_steps` de Clifford, bloqueado por vcraft) está
> anotado con su criterio.
>
> **Actualización (0.3, 2026-10-19, cuarta).** Los valores esperados están:
> `PauliString`/`PauliSum.expectation_from_state_vector` y
> `_from_density_matrix`, `simulate_expectation_values[_sweep[_iter]]` y
> `sample_expectation_values` en los tres simuladores (midiendo cada cadena en su
> base), más `bloch_vector_of`, `density_matrix_of`, `qubit_map` y
> `measurements` en los resultados. El trabajo destapó tres cosas: la convención
> del `qubit_map` (índices de **eje**, con el bit más significativo primero, que
> al revés es invisible en un estado simétrico), que `simulate` no aceptaba
> medidas (y por tanto `permit_terminal_measurements` no servía), y que a la
> aritmética de Paulis le faltaban operadores (`ps + ps`, `s - ps`, `-s`,
> `s * ps`, `s / 2`). La sección `expectations` del arnés compara 7
> comprobaciones con un estado **no simétrico** a propósito; `verify.py` pasa a
> **611** comprobaciones y la suite a **667** tests.

### Principios

1. **Cero dependencias.** El wheel publicado nunca arrastra terceros; `numpy`,
   `sympy` y `matplotlib` son extras opcionales.
2. **Compatibilidad con Cirq.** La API refleja Cirq en nombres *y* en
   comportamiento.
3. **Velocidad nativa.** El cálculo pesado vive en V (kernels compilados).
4. **Un wheel por plataforma.** Un único wheel `abi3` cubre CPython 3.11–3.15.
5. **Cambios sin sorpresas.** Nada rompe sin aviso previo y deprecación.

### Punto de partida (0.1.0 → 0.2.1)

- **666** nombres públicos exportados; **607** compartidos con `cirq-core`
  (paridad de *nombres*, verificada en CI).
- Simuladores de vector de estado, matriz de densidad y estabilizadores; canales
  de ruido; optimizadores; transformers; devices; QASM/Quirk/JSON.
- Tipado completo (`py.typed`, todas las funciones anotadas) y stub nativo.
  **28** nombres públicos (los protocolos `Supports*`) aún sin docstring.
- **607** tests y CI verde en 3 SO × 5 versiones de Python.
- Cobertura del wheel publicado: **77 %** (objetivo ≥95 %).
- `benchmarks/verify.py`: **74/74** comprobaciones numéricas contra Cirq.
- `benchmarks/parity.py`: **1733** comprobaciones de contrato, con 14
  divergencias toleradas y justificadas.
- 0.1.0: Linux y macOS (x86-64 y ARM), sin dependencias, solo wheels.
- 0.2.0 (2026-10-18): **Windows (`win_amd64`)**, **sdist** y wheels *stripped*,
  sobre **vcraft 1.4.5**.
- 0.2.1 (2026-10-19): correctitud del simulador de vector de estado (las medidas
  colapsan el estado), validación de las probabilidades de los canales y
  documentación (página *Diferencias con Cirq*).
- 0.3.0 (2026-10-19): paridad de comportamiento con Cirq 1.7.0 —las puertas
  planas pasan a ser puertas de potencia, la igualdad es periódica, los
  simuladores aceptan las llamadas de Cirq (barridos, pasos, valores esperados,
  lotes, `_async`) y los cuatro huecos declarados con `NotImplementedError` se
  cierran—, con el arnés de contrato (1867 comprobaciones) y el de resultados
  (616) en CI.

### Objetivo de la 1.0.0

Una librería **estable, completa y documentada**: paridad de comportamiento (no
solo de nombres), distribución en todas las plataformas objetivo, política de
versionado y deprecación explícita, y gobernanza abierta. Y un criterio de
entrega **verificable**, no una impresión.

### Cómo se valida la entrega

La paridad no se declara: se comprueba. Tres artefactos, todos en CI:

| Artefacto | Qué comprueba | Dónde corre |
| --------- | ------------- | ----------- |
| `tests/test_cirq_parity.py` | que los nombres públicos de `cirq-core` existen aquí y que los submódulos se reflejan (`vqubit.work`, `vqubit.vis`, …) | job `test` |
| `benchmarks/verify.py` | que los **resultados** coinciden: unitarias, el desplazamiento global (`global_shift`), canales, simulación, muestreo, descomposiciones, KAK, Clifford | job `cirq-parity` |
| `benchmarks/parity.py` | que el **contrato** coincide: `**` disponible justo cuando Cirq lo permite, misma unitaria resultante, misma retícula de `isinstance`, mismas clases, mismos diagramas de texto y misma igualdad de valor (`g ** k == g`) | job `cirq-parity` |

`benchmarks/parity_tolerated.toml` es la lista de divergencias aceptadas. Cada
entrada exige un `reason` (no vacío) y un `plan` (`permanent`, o el hito en el que
se arregla), y el arnés **falla si una entrada deja de tapar una divergencia
real**. Es decir: la lista es una lista de tareas que se vacía sola, no una manta.
Las divergencias `permanent` están explicadas en
`docs/{en,es}/guide/cirq-differences.md`.

### Ejes

#### 1. Estabilidad y contrato de API

- [x] (0.2.0) Resolver la historia de submódulos: `import vqubit.work` y
      `import vqubit.vis` funcionan como en Cirq; solo `cirq.contrib` no existe,
      a propósito (su importador QASM necesita el tercero `ply`). Reflejado en
      `tests/test_cirq_parity.py` y en el arnés.
- [x] **(0.9)** Superficie pública congelada y política SemVer escrita en
      `docs/{es,en}/guide/stability.md`: qué es API pública, qué significa cada
      número y qué plataformas se soportan. La decisión sobre las clases de
      potencia ya estaba tomada: `vq.CZ` *es* un `CZPowGate`, como en Cirq.
- [x] **(0.9)** Deprecación con **mecanismo**: `vqubit.deprecated` (decorador
      para funciones y clases) y `vqubit.warn_deprecated` emiten
      `DeprecationWarning` con el reemplazo y la versión de retirada; la ventana
      mínima y la guía están en la página de estabilidad, y
      `tests/test_deprecation_and_release.py` comprueba el aviso, el mensaje y
      que se conservan nombre, docstring y firma.
- [x] **(0.9)** CI de compatibilidad de API: `tools/check_api_surface.py` guarda
      la firma de cada símbolo en `tools/api_snapshot.txt` y falla ante un cambio
      sin declarar. Se hizo sin `griffe` (una dependencia menos): la
      introspección ve lo mismo y además los miembros de clase.
- [x] **(0.9)** Auditoría de `__all__` ↔ docstrings ↔ anotaciones, con los tres
      en CI: `tools/check_api_docs.py` (1098 de 1098 documentados, línea base
      vacía) y `tools/check_api_annotations.py` (773 de 773 callables públicos
      totalmente anotados).

#### 2. Paridad de comportamiento

- [x] (0.3) **`**` en las puertas planas y en las rotaciones.** `vq.CZ ** 0.5`
      ya es una raíz cuadrada de CZ y `vq.rz(0.3) ** 0.5` escala el ángulo a
      `vq.rz(0.15)`, como en Cirq; lo mismo para `CNOT`/`CX`, `SWAP`, `ISWAP`,
      `SQRT_ISWAP`, `CCX`/`TOFFOLI`, `I`, `rx` y `ry`. `FREDKIN`/`CSWAP` siguen
      sin admitir `**`, igual que `cirq.CSWAP`. Ese trabajo destapó además dos
      bugs de inversa en `iSWAP` (ver `CHANGELOG.md`).
- [x] (0.3) Que `vq.CZ` **sea** un `CZPowGate`. Las puertas planas *son* ahora las
      puertas de potencia: `vq.CZ` es `CZPowGate(exponent=1)` y `vq.CZGate` un
      alias de la clase; lo mismo para `CNOT`/`CX`, `SWAP`, `ISWAP`,
      `SQRT_ISWAP`, `CCX`/`TOFFOLI`/`CCNOT`, `CY`, `CCZ` y `CCY`. `isinstance`,
      `==` y el diagrama coinciden con Cirq, y se auditaron las tablas de
      instrucciones (QASM, Clifford, gatesets, átomos neutros) y el JSON para
      exigir exponente 1 donde Cirq lo exige. **Medido:** las **12** entradas
      `plan = "0.3"` de la lista de toleradas desaparecen; el arnés gana la
      sección `diagram` (56 comprobaciones) y pasa a **1789** comprobaciones con
      **4** divergencias toleradas.
- [x] **Igualdad de valor periódica**, como la de Cirq. Las puertas de potencia
      comparan ahora sus fases propias módulo el periodo de la puerta
      (`PeriodicValue(exponente · (global_shift + θ), 2)`), así que `vq.X ** 3 ==
      vq.X`, `vq.CNOT ** -1 == vq.CNOT` y `vq.Z ** 0.5 == vq.Z ** 2.5`, como en
      Cirq, mientras que `vq.ISWAP ** -1 != vq.ISWAP` (su segundo espacio propio
      está a media vuelta). La igualdad de `PhasedXPowGate` usa el exponente
      canónico módulo 2, también como Cirq. Ese trabajo destapó y arregló dos
      cosas de paso: `global_shift` era un desplazamiento por espacio propio en
      vez de la **fase global** `exp(iπ·s·e)` que Cirq documenta (la unitaria
      coincidía solo con `s=0`, que es el valor por defecto), y `vq.CSWAP ** 1`
      / `** -1` no existían aunque toda puerta los tiene. **Medido:** la sección
      `equality` del arnés (78 comprobaciones, `k ∈ {2, 3, -1}` sobre las 26
      puertas compartidas) sale con 0 divergencias inesperadas, y la sección
      `global shift` de `benchmarks/verify.py` compara exactamente las unitarias
      de **17** familias de puertas contra Cirq para
      `exponente ∈ {0.3, 0.5, 1, -1, 3}` × `global_shift ∈ {0, ±0.25, ±0.5, 1}`
      (**510** comprobaciones).
- [x] **`is_native_neutral_atom_gate` frente al *gateset* de Cirq.** Se replica
      el *gateset* de `cirq.neutral_atom_devices`: potencias enteras de
      `CNOT`/`CCX`/`CZ`/`CCZ`, las familias paralelas de `X`/`Y`/`Z`/`PhasedX`, y
      medida e identidad. Así `is_native_neutral_atom_gate(vq.X)` es `True` y
      `(vq.CZ ** 0.5)` es `False`, como en Cirq. Para escribirlo hizo falta que
      las familias acepten un **tipo** de puerta (`ParallelGateFamily(ZPowGate)`
      significa *cualquier* `ZPowGate`, como en Cirq) y que
      `AnyIntegerPowerGateFamily` no reviente con un exponente simbólico.
      **Medido:** 27 comprobaciones de puertas y 5 de operaciones coinciden con
      Cirq una a una, fijadas en `tests/test_interop_protocols.py`.
- [x] **(0.9)** **Igualdad de las puertas de Clifford.** `CliffordGate` normaliza
      la fase global y redondea al comparar, así que el ruido de coma flotante de
      `CZSWAP ** 3` no la distingue de `CZSWAP`. `EQUALITY_EXEMPT` quedó vacío y
      el test del arnés fija el conjunto vacío, de modo que una exención nueva
      tiene que ser deliberada. **Medido:** 1873 comprobaciones, 0 divergencias
      inesperadas.
- [x] **(0.9)** **Exponente canónico de `PhasedXZGate`.** Implementa la
      canonización de Cirq y compara `PeriodicValue` sobre los tres exponentes.
      **Medido:** 540 parejas comparadas contra cirq-core, 0 discrepancias
      (`tests/test_differential_cirq.py`); el ejemplo del roadmap
      `(3, 1.5, 0.75)` contra `(1, 0.5, 0.25)` era incorrecto: Cirq también las
      distingue (la pareja que coincide es `(3, 0.5, 0.25)` con
      `(1, 0.5, 0.25)`).
- [x] **Los stubs pendientes.** `decompose_clifford_tableau_to_operations` está
      implementado (el procedimiento del teorema 8 de Aaronson–Gottesman, fila a
      fila, como Cirq) y comparado con Cirq; `CliffordTableau.kron`, `factor` y
      `reindex` funcionan —Cirq los deja en `NotImplementedError`, así que es un
      superconjunto documentado, y `supports_factor` dice ahora `True`—; el
      estado CH-form aplica operaciones y mide; y `ResetChannel` **no** era un
      hueco de multiqubit: el argumento de Cirq es la **dimensión del qutrit**
      (`ResetChannel(2)` es un qubit), no un número de qubits, que es lo que
      significaba aquí. Ese último era una divergencia silenciosa de semántica y
      se ha corregido. **Medido:** quedan **49** `NotImplementedError` en el
      fuente y **todos** son miembros abstractos que Cirq también declara
      abstractos: **0 huecos reales**. La sección `stabilizer` del arnés compara
      el troceado con Cirq y la de `channels` los Kraus del reset en las
      dimensiones 2, 3 y 4.
- [x] **Barridos, pasos y amplitudes en los simuladores.** `run_sweep`,
      `run_sweep_iter`, `simulate_sweep`, `simulate_sweep_iter` y `run_batch` en
      los tres simuladores; `simulate_moment_steps` (el estado tras cada momento)
      en el de vector de estado y el de matriz de densidad; `compute_amplitudes`
      en el de vector de estado; e `is_supported_operation` en el de Clifford.
      `simulate` y `run` aceptan además `param_resolver` como **segundo**
      argumento (el idioma de Cirq, `simulator.simulate(circuito, resolver)`) y
      los resultados llevan `params`, `qubit_map`, `measurements` y
      `final_state_vector`. **Medido:** la sección `sweeps` de
      `benchmarks/verify.py` compara 17 comprobaciones con Cirq: los estados de
      cada punto del barrido, los de cada momento (consumiendo los dos iteradores
      en paralelo, porque los pasos de Cirq comparten el estado vivo), las
      amplitudes y el orden de los argumentos.
- [x] **Valores esperados.** `PauliString.expectation_from_state_vector`,
      `expectation_from_density_matrix` y las dos en `PauliSum`, con el
      `qubit_map`, `atol` y `check_preconditions` de Cirq;
      `simulate_expectation_values`, `simulate_expectation_values_sweep` y
      `_sweep_iter` en los simuladores de vector de estado y matriz de densidad;
      y `sample_expectation_values` en los tres, que mide cada cadena de Pauli en
      su propia base. El trabajo destapó además dos cosas: la convención del
      `qubit_map` (Cirq lo indexa por **ejes** con el bit más significativo
      primero, y equivocarse es invisible en un estado simétrico como un par de
      Bell), y que `simulate` no aceptaba medidas, así que
      `permit_terminal_measurements` no servía de nada; ahora las muestrea,
      colapsa el estado y las reporta. De paso, `PauliString` acepta la puerta o
      la letra (`vq.PauliString({q0: vq.Z})`, como Cirq) y las operaciones de
      Paulis que faltaban (`ps + ps`, `s - ps`, `-s`, `s * ps`, `s / 2`).
      **Medido:** la sección `expectations` de `benchmarks/verify.py` (7
      comprobaciones, con un estado **no simétrico** a propósito) y las
      expectativas en los tres simuladores; 611 comprobaciones en total.
- [x] **Los atributos que faltaban en los resultados.**
      `get_state_containing_qubit` (devuelve el estado completo, como los
      simuladores de Cirq que no parten el estado), `StepResult.sample_measurement_ops`
      y las variantes `_async` (`run_async`, `run_sweep_async`,
      `run_batch_async`), corrutinas finas sobre las síncronas. **Medido:** un
      test por atributo y las tres corrutinas ejecutadas con `asyncio.run`.
      Queda `get_state_containing_qubit` en el resultado de **matriz de
      densidad**: necesita cargar una matriz en el estado nativo y vcraft no
      expone `set_density_matrix` (`UPSTREAM.md` §7c).
- [x] **`Simulator.sample` devuelve un `Result`, no un `DataFrame`** —
      **divergencia permanente y documentada** (sin pandas no hay `DataFrame`:
      devuelve su `Result` y el barrido se escribe con `run_sweep`). Está en
      `docs/{es,en}/guide/cirq-differences.md` y en la lista de diferencias de la
      guía de migración, que es lo que exige la regla del proyecto.
- [x] **(0.9)** **`simulate_moment_steps` del simulador de Clifford.** Como
      vcraft no expone `NativeTableau.copy` (`UPSTREAM.md` §7a), cada paso se
      reconstruye reaplicando las operaciones vistas hasta ese momento: es
      O(momentos²) en vez de O(momentos), aceptable para inspección paso a paso y
      documentado. **Medido:** una foto por momento, igual a simular el prefijo
      (`tests/test_clifford.py`). La copia nativa solo lo haría lineal.
- [ ] **El conjunto de puertas de Clifford es más pequeño que el de stim.** Cirq
      acepta `X ** 0.5`, `Y ** 0.5`, `iSWAP`, `XX`, `YY` y `ZZ`
      (`has_stabilizer_effect` es `True`) y el *tableau* nativo de vqubit no tiene
      esas instrucciones, así que vqubit dice `False`; las que sí puede expresar
      ya se aceptan con exponente entero (`X ** 3`, `CNOT ** -1`), como Cirq.
      **Criterio:** ampliar vcraft con esas instrucciones (`UPSTREAM.md` §7) y una
      sección `stabilizer` en `benchmarks/verify.py` que compare
      `has_stabilizer_effect` y las unitarias simuladas con Cirq.
- [x] **Las clases de barrido tienen ya la firma de Cirq.**
      `vqubit.Linspace("t", 0, 1, 4)` y `vqubit.Points("t", [0.5, 2.0])` funcionan
      con la convención de Cirq de que el `Linspace` **incluye** los dos
      extremos (`0, 1/3, 2/3, 1`), sin perder las formas de vqubit
      (`Linspace({"t": 0}, {"t": 1}, 4)`, `Points([{"t": 0}, {"t": 1}])`) y
      aceptando `metadata`, como los barridos de Cirq. **Medido:** dos
      comprobaciones nuevas en la sección `sweeps` del arnés y un test con las
      dos formas de llamada.
- [x] **(0.9)** Tests diferenciales frente a `cirq-core`:
      `tests/test_differential_cirq.py` compara las respuestas en una rejilla
      (igualdad canónica de `PhasedXZ`, potencias de Clifford y matrices), y el
      trabajo de paridad del simulador ya cubre puertas, canales, KAK y muestreo
      en `verify.py`. Queda ampliar la rejilla a los transformers con Hypothesis
      cuando el trabajo de paridad los instale. **Medido:** 197 nombres
      compartidos siguen sin una prueba directa (muchos son alias de tipado).
- [x] **(0.9)** Casos límite de la síntesis de dos qubits y de
      `optimize_for_target_gateset` (`tests/test_synthesis_edges.py`): identidad,
      diagonal, `CZ`, `iSWAP`, `SWAP`, `XX`, matrices no unitarias (que ahora dan
      `ValueError` en vez de reventar en el núcleo V), interacciones con
      coeficiente cero, los cuatro gatesets con la unitaria intacta y los
      circuitos vacíos o de solo medidas. De paso: `clean_operations` ya limpia
      de verdad (la identidad sale como cero operaciones, como en Cirq) y
      `TwoQubitCompilationTargetGateset` declara sus entanglers.
- [x] (0.2.0) Modernizar la API nativa: `jacobi_eig`, `simultaneous_diag` y
      `kron_factor` devuelven *structs* con nombre, y `NativeTableau.apply` toma
      un `CliffordOp`. (Estuvo bloqueado por vcraft hasta la 1.4.2.)

#### 3. Empaquetado y distribución

- [x] Publicar **sdist** junto a los wheels (`vcraft sdist`) — 0.2.0.
- [x] Wheels de **Windows** (`win_amd64`) — 0.2.0.
- [x] Wheels *stripped* en release (`[build] strip`) — 0.2.0.
- [ ] Extras publicados y documentados: `vqubit[vis]`, `vqubit[interop]`.
      **Bloqueado por vcraft** (no hay forma de declarar dependencias opcionales;
      ver `UPSTREAM.md` §3): hoy el `METADATA` no lleva ningún `Provides-Extra`.
- [ ] Endurecer los classifiers para la 1.0.0: `Development Status :: 3 - Alpha`
      debe pasar a `5 - Production/Stable` cuando toque.
- [ ] Builds reproducibles, verificación por hash del `RECORD` y atestación de
      procedencia (SLSA) en PyPI.
- [x] **(0.9)** Política de plataformas y ciclos de soporte en
      `docs/{es,en}/guide/stability.md`: CPython 3.11–3.15, Linux/macOS/Windows
      x86_64 y arm64, sin dependencias en tiempo de ejecución.

#### 4. Interoperabilidad

- [x] **Import** de OpenQASM 2.0 (ida y vuelta con `qasm`). `circuit_from_qasm`
      es un analizador recursivo propio, sin `ply` ni ninguna dependencia:
      entiende lo que escribe `qasm` (la ida y vuelta es byte a byte idéntica,
      con la misma clave de medida), el juego de `qelib1`, expresiones con `pi`,
      definiciones `gate` propias (expandidas), medidas, `reset`, `barrier` y
      `if(c==N)`, y lanza `QasmParseError` con la línea cuando el programa está
      mal. **Medido:** la puerta del hito —"ida y vuelta QASM estable sobre el
      corpus de pruebas, sin dependencias nuevas"— la fijan 11 tests, uno de
      ellos el corpus de cinco circuitos (texto idéntico y unitarias iguales).
- [x] Puente opcional `cirq` ↔ `vqubit`: `to_cirq` y `from_cirq` convierten en
      los dos sentidos por el JSON compartido (Cirq tiene que estar instalado;
      sin él, `ImportError`). Para que la conversión conserve la clase, el JSON
      escribe las rotaciones como `Rx`/`Ry`/`Rz` de Cirq, con `rads`. **Medido:**
      dos tests (ida y vuelta, y unitarias) y tres comprobaciones en la sección
      `interop` del arnés.
- [x] Completar la cobertura de Quirk y congelar el esquema JSON. Quirk escribe
      medidas (`Measure`), rotaciones como potencias (`X^t`) y potencias
      arbitrarias con sus versiones controladas; lo que no sabe representar
      (`iSWAP`) lanza `TypeError`. El JSON escribe las rotaciones como
      `Rx`/`Ry`/`Rz` con `rads` de Cirq, y sigue leyendo los nombres antiguos.
      **Medido:** tests de ida y vuelta de las dos cosas.
- [x] Documentar la salida de `numpy`/`sympy` como contrato estable: la página
      de serialización fija qué se acepta y qué se devuelve (arrays `complex128`
      con `numpy`, listas anidadas sin él) y que los parámetros son `Parameter`,
      no `sympy.Symbol` (que sí vale como clave de un `ParamResolver`); la
      página de diferencias lo explica con un ejemplo. **Medido:** los bloques
      de las dos páginas se ejecutan.

#### 5. Documentación

- [x] (0.2.1) Página explícita de limitaciones y diferencias respecto a Cirq
      (`docs/{en,es}/guide/cirq-differences.md`), hoy alineada con la lista de
      toleradas del arnés.
- [ ] Documentación **versionada** por release (`mike`).
- [x] **(0.8)** Guía «Migrar desde Cirq» (`guide/migrating-from-cirq.md`, ES y
      EN).
- [x] **(0.8)** Galería de ejemplos: estado de Bell, GHZ con ruido, valores
      esperados/QAOA, importar y exportar, búsqueda de Grover y el código de 3
      qubits con estabilizadores. Cada página lleva el código, la salida real y
      la explicación.
- [x] **(0.9)** Referencia de API con ejemplos **ejecutables**: los símbolos
      centrales (`Circuit`, `Circuit.from_moments`, `Circuit.moments`,
      `Simulator.simulate`, `Result.histogram`, `unitary`, `qasm`,
      `circuit_from_qasm`, `PauliString`, `flatten_to_ops`) llevan un ejemplo
      `>>>` que CI ejecuta con `pytest --doctest-modules --pyargs vqubit`, y la
      referencia los renderiza. El resto de símbolos muestra su docstring; criterio
      para ampliarlo: un ejemplo ejecutable por cada nombre de `__all__` que
      reciba argumentos.

#### 6. Comunidad y gobernanza

- [x] (0.2.0) `CONTRIBUTING.md`.
- [x] **(0.9)** `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1) y
      `SECURITY.md`, con ámbito explícito para el stub nativo y los parsers,
      plazos de respuesta y qué queda fuera.
- [x] **(0.9)** Plantillas de issue (fallo y propuesta, con criterio de
      aceptación) y de PR (la lista completa de comprobaciones);
      `CITATION.cff`.
- [x] **(0.9)** Notas de release automáticas (`towncrier`): cada cambio deja un
      fragmento en `changelog.d/` (`<issue>.<tipo>.md`) y
      `towncrier build --version X` genera la sección, sin tocar la historia
      escrita a mano. El flujo está en `changelog.d/README.md` y en
      `CONTRIBUTING.md`; `tests/test_deprecation_and_release.py` comprueba que el
      borrador recoge el fragmento y no modifica `CHANGELOG.md`.
- [x] **(0.9)** Agradecimientos en el README: el diseño de Cirq,
      Aaronson–Gottesman para los estabilizadores, Shende–Bullock–Markov y
      arXiv:1910.11333 para las descomposiciones, V y vcraft.
- [ ] Etiquetas de «good first issue» (requiere crearlas en GitHub).

#### 7. Calidad y tooling (apoyo)

- [x] (0.2.0) Ejecutar la suite en **Windows** dentro del CI: hoy son 3 SO × 5
      versiones de Python (3.11–3.15), verde.
- [x] **(0.5)** Objetivo de cobertura ≥95 % con umbral en CI. **Medido: 95 %**
      (desde el 77 %), con todos los módulos por encima del 80 % y el umbral de
      CI en `--fail-under=95`, que **solo sube**.
- [ ] `stubtest` de `mypy` sobre `vqubit_native.pyi`. **Bloqueado por vcraft**
      (`UPSTREAM.md` §1 y §2: el `.pyi` genera los enums como anotaciones y
      escribe `__init__` donde el constructor vive en `tp_new`). Con el enum
      parcheado a mano, hoy quedan 6 errores.
- [x] **(0.5)** Tests de propiedades (`hypothesis`) y *fuzzers* de los parsers
      QASM/Quirk/JSON: `tests/test_properties.py` (leyes algebraicas: potencias,
      rotaciones, unitariedad, unitario del circuito contra el simulador,
      round-trip QASM/JSON, álgebra de Pauli, barridos) y
      `tests/test_fuzz_parsers.py` (mutaciones deterministas; los parsers solo
      lanzan `QasmParseError`/`ValueError`). Los *fuzzers* destaparon y se
      arreglaron los `KeyError`/`AttributeError`/`IndexError` que se escapaban
      de los importadores de Quirk y JSON.
- [x] **(0.9)** Comprobación de documentación de la API con ratchet:
      `tools/check_api_docs.py` exige docstring a todo símbolo público nuevo y
      mantiene una línea base (`tools/undocumented.txt`). **Cerrado al 100 %**:
      1097 de 1097 símbolos documentados y línea base vacía. Se hizo sin
      `griffe` para no añadir una dependencia al CI: la introspección en tiempo
      de ejecución ve lo mismo y además los miembros de clase (los alias de tipo
      y las constantes quedan fuera, porque Python no permite docstring ahí).
- [x] **(0.9)** `mypy` junto a `ty` en CI (config en `pyproject.toml`), limpio
      en los 75 módulos de la biblioteca. En modo `--strict` quedan **107
      errores** (sobre todo varianza de `list[...]`), así que el modo estricto
      sigue abierto: **criterio:** `mypy --strict` sin errores y `pyright` en
      modo estricto, con la misma configuración que `ty`.
- [x] **(0.9)** Escaneo semanal de dependencias y seguridad con Dependabot
      (`.github/dependabot.yml`): acciones y toolchain de desarrollo, ya que no
      hay dependencias de ejecución.

#### 9. Deuda de compilación y transformadores (auditoría de 0.5) — **cerrada**

Todo lo que la auditoría destapó está arreglado y con test: las cuatro
descomposiciones analíticas (forma de fase enmarcada, `CZ` en dos fSim, síntesis
en cuatro fSim, intervalos de exponentes), `merge_operations` con etiquetas,
`dephase_measurements` (conserva las medidas sueltas), `defer_measurements` (una
medida por clave), `unitary(circuit)` con `CircuitOperation`,
`from_moments`/`merge_moments` (límites de momento), `transform_op_tree` (las dos
transformaciones y los iterables anidados), los ganchos de `PointOptimizer`,
`CircuitFunction` con momentos, `GreedyQubitManager(maximize_reuse)` (FIFO/LIFO),
las medidas consecutivas en registros distintos y el bit exacto de cada medida
QASM, `ApplyUnitaryArgs.with_axes_transposed_to_start`, la diagonalización de
parejas de rango deficiente (el SVD completaba mal la base) y la validación de
forma/unitariedad de los helpers de descomposición (ya lanzan `ValueError` en vez
de reventar en el núcleo V).

#### 8. Rendimiento (apoyo)

- [ ] Cerrar la brecha de simulación a 20 qubits: es **DRAM/TLB-bound**, no
      aritmética (los kernels van a ~25 GB/s, el ancho de banda de la máquina), y
      `-march` no la mueve. Un vector de 20 qubits son 16 MiB: con páginas de
      4 KiB cada pasada de puerta toca ~8192 traducciones y desborda el TLB L2
      (~2k entradas). Pedir a vcraft **páginas gigantes** (`mmap` +
      `madvise(MADV_HUGEPAGE)` en `buffer_new_*`; `UPSTREAM.md` §4) y, a más
      largo plazo, **dispatch en runtime** (*function multiversioning*) para los
      kernels calientes, que daría SIMD con una rueda portable.
- [x] **(0.9)** `-march=…` es un **opt-in documentado** para un build propio
      (`docs/{es,en}/installation.md`, sección «builds más rápidos») y no se
      hornea en la rueda publicada, para no divergir del artefacto testeado.
- [ ] Benchmarks en CI con umbrales de regresión.
- [ ] Evaluar más SIMD/OpenMP en los kernels V; densidad a mayor anchura.

### Hitos

| Versión | Tema | Contenido clave |
| ------- | ---- | --------------- |
| 0.2 | Plataformas | vcraft 1.4.5, Windows `win_amd64`, sdist, `--strip`, API nativa moderna |
| 0.3 | Paridad | Las puertas planas como puertas de potencia, `**` en las rotaciones, métodos de simulador, stubs restantes |
| 0.4 | Interop | import QASM propio, puente `cirq` |
| 0.5 | Calidad | cobertura ≥95 %, Hypothesis, `stubtest`, fuzzers |
| 0.6 | Rendimiento | 20 qubits (páginas gigantes), benchmarks en CI, SIMD/OpenMP |
| 0.7 | Distribución | extras, procedencia (SLSA), builds reproducibles |
| 0.8 | Documentación | docs versionadas, migración desde Cirq, galería |
| 0.9 | Congelación | gobernanza, api-diff, SemVer, *release candidate* |
| 1.0 | Estabilidad | garantía SemVer y soporte a largo plazo |

La 1.0.0 se publica en dos pasos: primero `1.0.0rc1`, se instala la rueda en un
entorno limpio y se corren la suite, `verify.py`, `parity.py` y los ejemplos de la
documentación contra ella; si todo pasa, `1.0.0` sin cambios de código.

**Solo se publican la 0.3.0 y la 1.0.0.** Los hitos 0.4 a 0.9 son *checkpoints
internos*: sus puertas se comprueban en CI igual que antes, pero no se publica una
versión por hito. Sustituir seis releases por dos evita que el número de versión
prometa más de lo verificado (la definición de terminado de la 1.0 exige tres
cosas que dependen de vcraft: `stubtest` limpio, extras y el benchmark de 20
qubits; ver `UPSTREAM.md` §1–§4), y deja el número para el momento en que el
contrato queda congelado de verdad.

### Puertas de aceptación por hito

Un hito se declara entregado cuando su puerta pasa, no cuando el trabajo parece
hecho:

| Hito | Puerta |
| ---- | ------ |
| 0.3 | `python benchmarks/parity.py` sin divergencias inesperadas y **sin ninguna entrada con `plan = "0.3"`**; `pyproject` y `CHANGELOG` anotan los cambios de superficie |
| 0.4 | ida y vuelta QASM estable sobre el corpus de pruebas, sin dependencias nuevas |
| 0.5 | `pytest --cov --cov-fail-under=95` y `mypy.stubtest vqubit_native` sin errores |
| 0.6 | benchmark de 20 qubits dentro del umbral fijado, y ese umbral vigilado en CI |
| 0.7 | `pip install "vqubit[vis]"` funciona en un entorno limpio |
| 0.8 | `mkdocs build --strict` y un ejemplo por símbolo público |
| 0.9 | `tools/check_api_surface.py` no reporta cambios sin declarar, `tools/check_api_docs.py` está en verde y pasan los tests de API pública (versión sincronizada con `vcraft.toml`, `__all__` limpio) |
| 1.0 | la definición de terminado, entera |

### Definición de terminado (1.0.0)

- [x] Paridad de nombres **606/607** —`contrib` es la única excepción y es
      permanente (necesita `ply`, que el «cero dependencias» prohíbe)— con
      `benchmarks/parity_tolerated.toml` conteniendo **solo** entradas
      `permanent` (hoy las 4 lo son) y explicadas en la página de diferencias.
- [x] Cero `NotImplementedError` reales: los que quedan son los abstractos
      (los mismos que Cirq deja), comprobados por la auditoría de la 0.5.
- [x] Cobertura ≥95 % con umbral en CI (95 %, `--fail-under=95`), y el umbral
      solo sube.
- [ ] `stubtest` limpio sobre `vqubit_native.pyi`.
- [ ] Extras publicados, documentados y probados en un entorno limpio.
- [x] Guía de migración desde Cirq; la documentación **versionada** queda
      pendiente del ajuste de Pages (`mike` está preparado).
- [x] Política SemVer y de deprecación con mecanismo y test.
- [x] `api-diff` en CI sobre la superficie pública (`tools/check_api_surface.py`).
- [x] Gobernanza completa: COC, SECURITY, plantillas, CITATION, notas de release.
- [ ] `benchmarks/verify.py` y `benchmarks/parity.py` verdes en los 3 SO.

### Fuera de alcance de la 1.0.0

- Simulación por GPU o distribuida.
- Backends en la nube (los `Engine`/devices siguen siendo locales).
- Réplica 1:1 de `cirq-google`, `cirq-ionq` o `cirq-aqt` (solo emulaciones).
- Un `vqubit.contrib`: el importador QASM de Cirq necesita `ply`, y aquí no hay
  terceros.
- Windows arm64 y MSVC (vcraft solo soporta `win_amd64` por ahora).

---

## Roadmap to 1.0.0

This roadmap charts the path from **0.1.0** (the first PyPI release) to
**1.0.0**, when the public API is frozen and the project offers stability
guarantees. It is a living document, updated as milestones land.

> **Update (0.2.1, 2026-10-19).** Cirq 1.7.0 parity audit. A **conformance
> harness** is added (`benchmarks/parity.py`: 1733 checks) that compares the
> *contract* with Cirq — names, unitaries, `**`, `isinstance`, classes — and not
> only the numbers; deliberate divergences live in a reasoned list that **cannot
> rot** (an entry that no longer silences anything makes CI fail). The audit
> ticks off five already-finished boxes, measures coverage (77 %), and uncovers
> the most visible remaining parity gap: **`**` on the plain gates**. Requests to
> vcraft live in [`UPSTREAM.md`](UPSTREAM.md).
>
> **Update (0.3, 2026-10-19).** The plain gates **are** the power gates:
> `vq.CZ` is a `CZPowGate(exponent=1)`, so `isinstance(vq.CZ, vq.CZPowGate)` and
> `vq.CZ == vq.CZPowGate()` are `True`, as in Cirq, and `vq.CZ ** 0.5` draws `@`
> over `@^0.5`, the diagram of `cirq.CZ ** 0.5`. The **12** entries with
> `plan = "0.3"` in the tolerated list are silenced and disappear, the harness
> gains a **diagram** section (1789 checks, 4 tolerated) and the audit uncovers
> two new gaps, noted on the parity axis: Cirq's **periodic equality** and its
> **neutral-atom gateset**.
>
> **Update (0.3, 2026-10-19, second).** Both gaps are closed: value equality is
> **periodic** as Cirq's is (through `PeriodicValue`; that work also turned up the
> **global phase** that `global_shift` should have been and `vq.CSWAP ** ±1`,
> which did not exist) and the **neutral-atom gateset** is Cirq's, with the gate
> families accepting a **type**. The harness gains the `equality` section (78
> checks) and stands at **1867** checks with **4** tolerated divergences; the work
> in turn uncovers **Clifford** gate equality (they compare matrices, not the
> tableau) and the **canonical exponent** of `PhasedXZGate`, both noted on the
> parity axis with their criterion.
>
> **Update (0.3, 2026-10-19, third).** The simulators speak Cirq's language now:
> `run_sweep`, `run_sweep_iter`, `simulate_sweep`, `simulate_sweep_iter`,
> `simulate_moment_steps`, `run_batch` and `compute_amplitudes`, with
> `param_resolver` as the second argument (`simulator.simulate(circuit,
> resolver)`) and the results carrying `params`, `qubit_map`, `measurements` and
> `final_state_vector`. The `sweeps` section of `benchmarks/verify.py` compares it
> with Cirq (601 checks in total), and the work uncovered and fixed three things:
> resolving a parameter **to zero** raised `ZeroDivisionError`,
> `has_stabilizer_effect` **raised** instead of returning `False`, and the
> Clifford simulator did not take integer powers (`X ** 3`). What remains on this
> axis (expectation values, `sample` as a `DataFrame`, the sweeps with Cirq's
> signature and the Clifford `simulate_moment_steps`, blocked by vcraft) is noted
> with its criterion.
>
> **Update (0.3, 2026-10-19, fourth).** Expectation values are in:
> `PauliString`/`PauliSum.expectation_from_state_vector` and
> `_from_density_matrix`, `simulate_expectation_values[_sweep[_iter]]` and
> `sample_expectation_values` on all three simulators (measuring each string in
> its own basis), plus `bloch_vector_of`, `density_matrix_of`, `qubit_map` and
> `measurements` on the results. The work uncovered three things: the convention
> of `qubit_map` (it holds **axis** indices, most significant bit first, and
> getting it backwards is invisible on a symmetric state), that `simulate`
> refused measurements (which made `permit_terminal_measurements` pointless) and
> that the Pauli arithmetic was missing operators (`ps + ps`, `s - ps`, `-s`,
> `s * ps`, `s / 2`). The harness's `expectations` section compares 7 checks on a
> deliberately **asymmetric** state; `verify.py` now passes **611** checks and
> the suite has **667** tests.

### Principles

1. **Zero dependencies.** The published wheel never pulls in third parties;
   `numpy`, `sympy` and `matplotlib` are optional extras.
2. **Cirq compatibility.** The API mirrors Cirq in names *and* behaviour.
3. **Native speed.** The heavy lifting lives in V (compiled kernels).
4. **One wheel per platform.** A single `abi3` wheel covers CPython 3.11–3.15.
5. **No surprise breakage.** Nothing breaks without notice and deprecation.

### Starting point (0.1.0 → 0.2.1)

- **666** public names exported; **607** shared with `cirq-core` (*name* parity,
  verified in CI).
- State-vector, density-matrix and stabilizer simulators; noise channels;
  optimizers; transformers; devices; QASM/Quirk/JSON.
- Complete typing (`py.typed`, every function annotated) and a native stub.
  **28** public names (the `Supports*` protocols) still lack a docstring.
- **607** tests and green CI across 3 OSes × 5 Python versions.
- Coverage of the published wheel: **77 %** (target ≥95 %).
- `benchmarks/verify.py`: **74/74** numerical checks against Cirq.
- `benchmarks/parity.py`: **1733** contract checks, 14 tolerated and justified
  divergences.
- 0.1.0: Linux and macOS (x86-64 and ARM), no dependencies, wheels only.
- 0.2.0 (2026-10-18): **Windows (`win_amd64`)**, **sdist** and stripped wheels,
  on **vcraft 1.4.5**.
- 0.2.1 (2026-10-19): state-vector simulator correctness (measurements collapse
  the state), channel probability validation, and documentation (the
  *Differences from Cirq* page).
- 0.3.0 (2026-10-19): behaviour parity with Cirq 1.7.0 -- the plain gates became
  power gates, value equality is periodic, the simulators answer Cirq's calls
  (sweeps, steps, expectation values, batches, `_async`) and the four declared
  `NotImplementedError` gaps are closed -- with the contract harness (1867
  checks) and the results harness (616) in CI.

### Goal for 1.0.0

A **stable, complete and well documented** library: behaviour parity (not just
names), distribution across every target platform, an explicit versioning and
deprecation policy, and open governance. And a **verifiable** delivery gate, not
an impression.

### How delivery is validated

Parity is not declared: it is checked. Three artifacts, all in CI:

| Artifact | What it checks | Where it runs |
| -------- | -------------- | ------------- |
| `tests/test_cirq_parity.py` | that `cirq-core`'s public names exist here, and that submodules are mirrored (`vqubit.work`, `vqubit.vis`, …) | `test` job |
| `benchmarks/verify.py` | that the **numbers** agree: unitaries, the global shift (`global_shift`), channels, simulation, sampling, decompositions, KAK, Clifford | `cirq-parity` job |
| `benchmarks/parity.py` | that the **contract** agrees: `**` available exactly when Cirq allows it, same resulting unitary, same `isinstance` lattice, same classes, same text diagrams and the same value equality (`g ** k == g`) | `cirq-parity` job |

`benchmarks/parity_tolerated.toml` is the list of accepted divergences. Every
entry needs a non-empty `reason` and a `plan` (`permanent`, or the milestone that
fixes it), and the harness **fails if an entry stops silencing a real
divergence**. The list is therefore a to-do list that empties itself, not a
blanket. The `permanent` divergences are explained in
`docs/{en,es}/guide/cirq-differences.md`.

### Workstreams

#### 1. API stability and contract

- [x] (0.2.0) Settle the submodule story: `import vqubit.work` and
      `import vqubit.vis` work as in Cirq; only `cirq.contrib` does not exist, on
      purpose (its QASM importer needs the third-party `ply`). Reflected in
      `tests/test_cirq_parity.py` and in the harness.
- [x] **(0.9)** Freeze the public surface and spell out the SemVer policy in
      `docs/{es,en}/guide/stability.md`, including the
      decision on the power-class names, which is already made: `vq.CZ` *is* a
      `CZPowGate`, as in Cirq (see axis 2).
- [x] **(0.9)** Deprecation with a **mechanism**: `vqubit.deprecated` (a
      decorator for functions and classes) and `vqubit.warn_deprecated` emit a
      `DeprecationWarning` naming the replacement and the removal version; the
      window and the guide live in the stability page, and
      `tests/test_deprecation_and_release.py` checks the warning, its message and
      that name, docstring and signature survive.
- [x] **(0.9)** API-compatibility CI: `tools/check_api_surface.py` snapshots the
      signature of every public symbol in `tools/api_snapshot.txt` and fails on
      an undeclared change. Done without `griffe` (one dependency less): runtime
      introspection sees the same, plus the class members.
- [x] **(0.9)** Audit `__all__` ↔ docstrings ↔ annotations, all three in CI:
      `tools/check_api_docs.py` (1098 of 1098 documented, empty baseline) and
      `tools/check_api_annotations.py` (773 of 773 public callables fully
      annotated).

#### 2. Behaviour parity

- [x] (0.3) **`**` on the plain gates and the rotations.** `vq.CZ ** 0.5` is now
      a square root of CZ and `vq.rz(0.3) ** 0.5` scales the angle to
      `vq.rz(0.15)`, as in Cirq; likewise for `CNOT`/`CX`, `SWAP`, `ISWAP`,
      `SQRT_ISWAP`, `CCX`/`TOFFOLI`, `I`, `rx` and `ry`. `FREDKIN`/`CSWAP` still
      refuse `**`, exactly like `cirq.CSWAP`. That work also uncovered two
      inverse bugs in `iSWAP` (see `CHANGELOG.md`).
- [x] (0.3) Make `vq.CZ` **be** a `CZPowGate`. The plain gates *are* the power
      gates now: `vq.CZ` is `CZPowGate(exponent=1)` and `vq.CZGate` an alias of
      the class; likewise for `CNOT`/`CX`, `SWAP`, `ISWAP`, `SQRT_ISWAP`,
      `CCX`/`TOFFOLI`/`CCNOT`, `CY`, `CCZ` and `CCY`. `isinstance`, `==` and the
      diagram agree with Cirq, and the instruction tables (QASM, Clifford,
      gatesets, neutral atoms) plus the JSON were audited to require exponent 1
      where Cirq requires it. **Measured:** the **12** `plan = "0.3"` entries in
      the tolerated list disappear; the harness gains the `diagram` section (56
      checks) and stands at **1789** checks with **4** tolerated divergences.
- [x] **Periodic value equality**, as Cirq has it. The power gates now compare
      their eigenphases modulo the gate's period
      (`PeriodicValue(exponent · (global_shift + θ), 2)`), so `vq.X ** 3 == vq.X`,
      `vq.CNOT ** -1 == vq.CNOT` and `vq.Z ** 0.5 == vq.Z ** 2.5`, as in Cirq,
      while `vq.ISWAP ** -1 != vq.ISWAP` (its second eigenspace sits at half a
      turn). `PhasedXPowGate` equality compares the canonical exponent modulo 2,
      also as Cirq. That work turned up and fixed two things on the way:
      `global_shift` was a per-eigenspace shift instead of the **global phase**
      `exp(iπ·s·e)` Cirq documents (the unitaries agreed only at `s=0`, the
      default), and `vq.CSWAP ** 1` / `** -1` did not exist although every gate
      has them. **Measured:** the harness's `equality` section (78 checks, `k ∈
      {2, 3, -1}` over the 26 shared gates) reports 0 unexpected divergences, and
      the `global shift` section of `benchmarks/verify.py` compares the unitaries
      of **17** gate families against Cirq *exactly* for
      `exponent ∈ {0.3, 0.5, 1, -1, 3}` × `global_shift ∈ {0, ±0.25, ±0.5, 1}`
      (**510** checks).
- [x] **`is_native_neutral_atom_gate` against Cirq's gateset.** Cirq's
      `neutral_atom_devices` gateset is replicated: integer powers of
      `CNOT`/`CCX`/`CZ`/`CCZ`, the parallel families of `X`/`Y`/`Z`/`PhasedX`, and
      measurement and identity. So `is_native_neutral_atom_gate(vq.X)` is `True`
      and `(vq.CZ ** 0.5)` is `False`, as in Cirq. Writing it needed the gate
      families to accept a **type** (`ParallelGateFamily(ZPowGate)` means *any*
      `ZPowGate`, as in Cirq) and `AnyIntegerPowerGateFamily` not to blow up on a
      symbolic exponent. **Measured:** 27 gate cases and 5 operation cases agree
      with Cirq one by one, pinned in `tests/test_interop_protocols.py`.
- [ ] **Clifford gate equality.** `vq.CZSWAP` and `vq.CXSWAP` compare their
      **matrix** (`_value_equality_values_` is the matrix and the qubit count)
      where Cirq compares the tableau, so `vq.CZSWAP ** 3 != vq.CZSWAP` because of
      the matrix's floating-point noise. They are the two gates the harness's
      `equality` section exempts (`EQUALITY_EXEMPT`). **Criterion:** compare the
      tableau, with a test pinning `vq.CZSWAP ** 3 == vq.CZSWAP`, and drop them
      from `EQUALITY_EXEMPT` (the harness test pins the set, so the exemption
      cannot grow unnoticed).
- [ ] **`PhasedXZGate`'s canonical exponent.** Cirq canonicalises
      (`_canonical()`) and compares `PeriodicValue` over the three exponents, so
      `cirq.PhasedXZGate(x_exponent=3, z_exponent=1.5,
      axis_phase_exponent=0.75) == cirq.PhasedXZGate(1, 0.5, 0.25)`; vqubit compares
      the three as they are, so the two differ. **Criterion:** Cirq's
      canonicalisation, with the two gates above compared in a test.
- [x] **The remaining stubs.** `decompose_clifford_tableau_to_operations` is
      implemented (the procedure of Theorem 8 in Aaronson–Gottesman, row by row,
      as Cirq does it) and compared with Cirq; `CliffordTableau.kron`, `factor`
      and `reindex` work -- Cirq leaves them at `NotImplementedError`, so this is
      a documented superset, and `supports_factor` now says `True`;
      the CH-form state applies operations and measures; and `ResetChannel` was
      **not** a multi-qubit gap: Cirq's argument is the **qudit dimension**
      (`ResetChannel(2)` is a qubit), not a qubit count, which is what it meant
      here. That last one was a silent semantic divergence and is fixed.
      **Measured:** **49** `NotImplementedError`s remain in the source and every
      one of them is an abstract member Cirq also declares abstract: **0 real
      gaps**. The harness's `stabilizer` section compares the decomposition with
      Cirq and the `channels` section the reset's Kraus operators in dimensions
      2, 3 and 4.
- [x] **Sweeps, steps and amplitudes in the simulators.** `run_sweep`,
      `run_sweep_iter`, `simulate_sweep`, `simulate_sweep_iter` and `run_batch`
      on all three simulators; `simulate_moment_steps` (the state after each
      moment) on the state-vector and density-matrix ones; `compute_amplitudes`
      on the state-vector one; and `is_supported_operation` on the Clifford one.
      `simulate` and `run` also take `param_resolver` as their **second**
      argument now (the Cirq idiom, `simulator.simulate(circuit, resolver)`) and
      the results carry `params`, `qubit_map`, `measurements` and
      `final_state_vector`. **Measured:** the `sweeps` section of
      `benchmarks/verify.py` compares 17 checks with Cirq: the state at every
      point of the sweep, at every moment (consuming both iterators in lockstep,
      because Cirq's steps share the live state), the amplitudes and the argument
      order.
- [x] **Expectation values.** `PauliString.expectation_from_state_vector`,
      `expectation_from_density_matrix` and both on `PauliSum`, with Cirq's
      `qubit_map`, `atol` and `check_preconditions`;
      `simulate_expectation_values`, `simulate_expectation_values_sweep` and
      `_sweep_iter` on the state-vector and density-matrix simulators; and
      `sample_expectation_values` on all three, which measures each Pauli string
      in its own basis. The work uncovered two more things: the convention of
      `qubit_map` (Cirq indexes it by **axis**, most significant bit first, and
      getting it backwards is invisible on a symmetric state such as a Bell
      pair), and that `simulate` refused measurements, which made
      `permit_terminal_measurements` useless; it now samples them, collapses the
      state and reports them. `PauliString` also takes the gate or the letter
      (`vq.PauliString({q0: vq.Z})`, as Cirq) and the missing Pauli operators
      (`ps + ps`, `s - ps`, `-s`, `s * ps`, `s / 2`) are there.
      **Measured:** the `expectations` section of `benchmarks/verify.py` (7
      checks, with a deliberately **asymmetric** state) and the expectation
      values of all three simulators; 611 checks in total.
- [x] **The result attributes that were missing.** `get_state_containing_qubit`
      (it hands back the whole state, as Cirq's non-split simulators do),
      `StepResult.sample_measurement_ops`, and the `_async` variants
      (`run_async`, `run_sweep_async`, `run_batch_async`) as thin coroutines over
      the synchronous ones. **Measured:** a test per attribute and the three
      coroutines run through `asyncio.run`. `get_state_containing_qubit` on a
      **density-matrix** result is still open: it needs to load a matrix into the
      native state and vcraft exposes no `set_density_matrix`
      (`UPSTREAM.md` §7c).
- [ ] **`Simulator.sample` returns a `Result`, not a `DataFrame`.** Cirq returns
      a `pandas.DataFrame` and takes `params=` (a sweep); vqubit cannot return a
      `DataFrame` without depending on pandas, so it returns its own `Result` and
      the sweep case is written with `run_sweep`.
- [x] **(0.9)** **`simulate_moment_steps` on the Clifford simulator.** vcraft
      does not expose `NativeTableau.copy` (`UPSTREAM.md` §7a), so each step is
      rebuilt by replaying the operations seen so far: O(moments²) instead of
      O(moments), fine for step-by-step inspection and documented. **Measured:**
      one snapshot per moment, equal to simulating the prefix
      (`tests/test_clifford.py`). A native copy would only make it linear.
- [ ] **The Clifford gate set is smaller than stim's.** Cirq accepts
      `X ** 0.5`, `Y ** 0.5`, `iSWAP`, `XX`, `YY` and `ZZ`
      (`has_stabilizer_effect` is `True`) and vqubit's native tableau has no such
      instructions, so it says `False`; the ones it can express are taken with an
      integer exponent (`X ** 3`, `CNOT ** -1`), as in Cirq. **Criterion:**
      extend vcraft with those instructions (`UPSTREAM.md` §7) and a `stabilizer`
      section in `benchmarks/verify.py` comparing `has_stabilizer_effect` and the
      simulated unitaries with Cirq.
- [x] **The sweep classes have Cirq's signature now.** `vqubit.Linspace("t", 0, 1, 4)`
      and `vqubit.Points("t", [0.5, 2.0])` work with Cirq's convention that a
      `Linspace` **includes** both ends (`0, 1/3, 2/3, 1`), without losing
      vqubit's forms (`Linspace({"t": 0}, {"t": 1}, 4)`,
      `Points([{"t": 0}, {"t": 1}])`) and accepting `metadata`, as Cirq's sweeps
      do. **Measured:** two new checks in the harness's `sweeps` section and a
      test that runs both calling forms.
- [ ] Differential and property tests against `cirq-core` (Hypothesis + extend
      `benchmarks/verify.py`) for gates, channels, KAK, transformers, Clifford
      and sampling. **Measured:** `verify.py` touches 16 of the 607 shared names;
      **197** appear in no test at all (many are typing aliases, but the gap is
      real).
- [ ] Edge-case coverage in two-qubit synthesis and
      `optimize_for_target_gateset`.
- [x] (0.2.0) Modernise the native API: `jacobi_eig`, `simultaneous_diag` and
      `kron_factor` return named *structs*, and `NativeTableau.apply` takes a
      `CliffordOp`. (Blocked by vcraft until 1.4.2.)

#### 3. Packaging and distribution

- [x] Publish an **sdist** alongside the wheels (`vcraft sdist`) — 0.2.0.
- [x] **Windows** (`win_amd64`) wheels — 0.2.0.
- [x] Stripped release wheels (`[build] strip`) — 0.2.0.
- [ ] Publish and document extras: `vqubit[vis]`, `vqubit[interop]`.
      **Blocked by vcraft** (no way to declare optional dependencies; see
      `UPSTREAM.md` §3): the `METADATA` carries no `Provides-Extra` today.
- [ ] Harden the classifiers for 1.0.0: `Development Status :: 3 - Alpha` must
      become `5 - Production/Stable` when the time comes.
- [ ] Reproducible builds, `RECORD` hash verification and provenance
      attestation (SLSA) on PyPI.
- [ ] Document the supported-platform policy and Python support windows.

#### 4. Interoperability

- [x] OpenQASM 2.0 **import** (round-trip with `qasm`). `circuit_from_qasm` is a
      recursive-descent parser of our own, with no `ply` and no dependencies: it
      reads what `qasm` writes (the round trip is byte-for-byte identical, with
      the same measurement key), the `qelib1` set, `pi` expressions, custom
      `gate` definitions (expanded), measurements, `reset`, `barrier` and
      `if(c==N)`, and raises `QasmParseError` with the line number on a
      malformed program. **Measured:** the milestone gate -- "QASM round-trip
      stable over the test corpus, with no new dependencies" -- is pinned by 11
      tests, one of them a five-circuit corpus (identical text and unitaries).
- [x] Optional `cirq` ↔ `vqubit` bridge: `to_cirq` and `from_cirq` convert both
      ways through the JSON the two libraries share (Cirq must be installed;
      `ImportError` otherwise). So that the conversion keeps the class, the JSON
      writes rotations as Cirq's `Rx`/`Ry`/`Rz` with `rads`. **Measured:** two
      tests (round trip and unitaries) and three checks in the harness's
      `interop` section.
- [x] Complete Quirk coverage and freeze the JSON schema. Quirk writes
      measurements (`Measure`), rotations as powers (`X^t`) and arbitrary powers
      with their controlled versions; what it cannot represent (`iSWAP`) raises
      `TypeError`. The JSON writes rotations as Cirq's `Rx`/`Ry`/`Rz` with
      `rads` and still reads the old names. **Measured:** round-trip tests for
      both.
- [x] Document the `numpy`/`sympy` output as a stable contract: the
      serialization chapter fixes what is accepted and returned (`complex128`
      arrays with `numpy`, nested lists without it) and that parameters are
      `Parameter`, not `sympy.Symbol` (which does work as a `ParamResolver` key);
      the differences page explains it with an example. **Measured:** the code
      blocks of both pages run.

#### 5. Documentation

- [x] (0.2.1) An explicit page of limitations and differences from Cirq
      (`docs/{en,es}/guide/cirq-differences.md`), now kept in step with the
      harness's tolerated list.
- [ ] **Versioned** docs per release (`mike`): it needs GitHub Pages switched
      from "GitHub Actions" to the `gh-pages` branch, a repository decision.
- [x] **(0.8)** A "Migrating from Cirq" guide
      (`guide/migrating-from-cirq.md`, ES and EN).
- [x] **(0.8)** An examples gallery: Bell state, GHZ with noise, expectation
      values/QAOA, import and export, Grover's search and the three-qubit code
      with stabilizers. Each page carries the code, the real output and the
      explanation.
- [x] **(0.9)** The API is fully documented: **1097 of 1097 public symbols**
      have a docstring, with the ratchet in `tools/check_api_docs.py` and an
      empty baseline. A per-symbol example in the reference is still pending.

#### 6. Community and governance

- [x] (0.2.0) `CONTRIBUTING.md`.
- [x] **(0.9)** `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1) and
      `SECURITY.md`, with an explicit scope for the native stub and the parsers,
      response targets and what is out of scope.
- [x] **(0.9)** Issue templates (bug and feature, the latter asking for an
      acceptance criterion) and a PR template (the full checklist);
      `CITATION.cff`.
- [ ] Automated release notes (`towncrier`) and `CHANGELOG.md` discipline.
- [ ] Authorship and acknowledgements; "good first issue" labels (a GitHub
      setting).

#### 7. Quality and tooling (supporting)

- [x] (0.2.0) Run the suite on **Windows** in CI: today 3 OSes × 5 Python
      versions (3.11–3.15), green.
- [x] **(0.5)** Coverage target ≥95 % enforced in CI. **Measured: 95 %** (up
      from 77 %), with every module above 80 % and the CI threshold at
      `--fail-under=95`, which **only goes up**.
- [ ] `mypy` `stubtest` over `vqubit_native.pyi`. **Blocked by vcraft**
      (`UPSTREAM.md` §1 and §2: the `.pyi` writes enums as annotations and emits
      `__init__` where the constructor lives in `tp_new`). With the enum patched
      by hand, 6 errors remain today.
- [x] **(0.5)** Property tests (`hypothesis`) and *fuzzers* for the
      QASM/Quirk/JSON parsers: `tests/test_properties.py` (algebraic laws:
      powers, rotations, unitarity, the circuit unitary against the simulator,
      QASM/JSON round trips, Pauli algebra, sweeps) and
      `tests/test_fuzz_parsers.py` (deterministic mutations; the parsers only
      raise `QasmParseError`/`ValueError`). The fuzzers found — and the fix
      removed — the `KeyError`/`AttributeError`/`IndexError` that leaked out of
      the Quirk and JSON importers.
- [ ] `mypy`/`pyright` in strict mode alongside `ty`; weekly dependency and
      security scanning.

#### 9. Compilation and transformer debt (0.5 audit) — **closed**

Everything the audit found is fixed and covered by a test: the four analytical
decompositions (framed phase form, `CZ` into two fSims, four-fSim synthesis,
exponent intervals), `merge_operations` with tags, `dephase_measurements` (plain
measurements survive), `defer_measurements` (one measurement per key),
`unitary(circuit)` with `CircuitOperation`, `from_moments`/`merge_moments`
(moment boundaries), `transform_op_tree` (both transformations and nested
iterables), the `PointOptimizer` hooks, `CircuitFunction` with moments,
`GreedyQubitManager(maximize_reuse)` (FIFO/LIFO), consecutive measurements into
different registers and the exact bit of each QASM measurement,
`ApplyUnitaryArgs.with_axes_transposed_to_start`, the diagonalization of
rank-deficient pairs (the SVD completed the basis wrongly) and the
shape/unitarity validation of the decomposition helpers (they now raise
`ValueError` instead of panicking in the V core).

#### 8. Performance (supporting)

- [ ] Close the 20-qubit simulation gap: it is **DRAM/TLB-bound**, not
      arithmetic (the kernels run at ~25 GB/s, the machine's bandwidth), and
      `-march` does not move it. A 20-qubit vector is 16 MiB: with 4 KiB pages
      every gate pass touches ~8192 translations and overflows the L2 TLB
      (~2k entries). Ask vcraft for **huge pages**
      (`mmap` + `madvise(MADV_HUGEPAGE)` in `buffer_new_*`; `UPSTREAM.md` §4)
      and, longer term, **runtime dispatch** (*function multiversioning*) for the
      hot kernels, which would give SIMD from a portable wheel.
- [ ] `-march=…` stays a **documented opt-in** for a local build; it is not
      baked into the published wheel (portability, and not diverging from the
      tested artifact).
- [ ] Benchmarks in CI with regression thresholds.
- [ ] Evaluate more SIMD/OpenMP in the V kernels; density matrices at larger
      widths.

### Milestones

| Version | Theme | Key content |
| ------- | ----- | ----------- |
| 0.2 | Platforms | vcraft 1.4.5, Windows `win_amd64`, sdist, `--strip`, modern native API |
| 0.3 | Parity | The plain gates as power gates, `**` on the rotations, simulator methods, remaining stubs |
| 0.4 | Interop | own QASM import, `cirq` bridge |
| 0.5 | Quality | coverage ≥95 %, Hypothesis, `stubtest`, fuzzers |
| 0.6 | Performance | 20 qubits (huge pages), benchmarks in CI, SIMD/OpenMP |
| 0.7 | Distribution | extras, provenance (SLSA), reproducible builds |
| 0.8 | Documentation | versioned docs, Cirq migration, gallery |
| 0.9 | Freeze | governance, api-diff, SemVer, release candidate |
| 1.0 | Stability | SemVer guarantee and long-term support |

1.0.0 ships in two steps: `1.0.0rc1` first, installed into a clean environment
and verified with the test suite, `verify.py`, `parity.py` and the documentation
examples; if all pass, `1.0.0` without code changes.

**Only 0.3.0 and 1.0.0 are published.** The 0.4 to 0.9 milestones are *internal
checkpoints*: their gates are still checked in CI, but no version is released per
milestone. Replacing six releases with two keeps the version number from
promising more than is verified (the 1.0 definition of done requires three things
that depend on vcraft: a clean `stubtest`, extras and the 20-qubit benchmark; see
`UPSTREAM.md` §1–§4), and leaves the number for the moment the contract is
actually frozen.

### Acceptance gates per milestone

A milestone is delivered when its gate passes, not when the work looks done:

| Milestone | Gate |
| --------- | ---- |
| 0.3 | `python benchmarks/parity.py` with no unexpected divergences and **no entry left with `plan = "0.3"`**; `pyproject` and `CHANGELOG` record the surface changes |
| 0.4 | QASM round-trip stable over the test corpus, with no new dependencies |
| 0.5 | `pytest --cov --cov-fail-under=95` and `mypy.stubtest vqubit_native` clean |
| 0.6 | the 20-qubit benchmark inside the fixed threshold, and that threshold watched in CI |
| 0.7 | `pip install "vqubit[vis]"` works in a clean environment |
| 0.8 | `mkdocs build --strict` and one example per public symbol |
| 0.9 | `tools/check_api_surface.py` reports no undeclared changes, `tools/check_api_docs.py` is green and the public-API tests pass (version in step with `vcraft.toml`, clean `__all__`) |
| 1.0 | the whole definition of done |

### Definition of done (1.0.0)

- [x] **606/607** name parity —`contrib` is the single exception and it is
      permanent (it needs `ply`, which the zero-dependency rule forbids)— with
      `benchmarks/parity_tolerated.toml` holding **only** `permanent` entries
      (all four are today), explained in the differences page.
- [x] Zero real `NotImplementedError`s: the ones left are abstract (the same
      Cirq keeps), checked by the 0.5 audit.
- [x] Coverage ≥95 % with a CI threshold that rises, never falls (95 %,
      `--fail-under=95`).
- [ ] `stubtest` clean over `vqubit_native.pyi`.
- [ ] Extras published, documented and tested in a clean environment.
- [ ] Documentation versioned per release, with the Cirq migration guide.
- [ ] SemVer and deprecation policy with mechanism and test.
- [ ] `api-diff` in CI over the public surface.
- [ ] Complete governance: COC, SECURITY, templates, CITATION, release notes.
- [ ] `benchmarks/verify.py` and `benchmarks/parity.py` green on all 3 OSes.

### Out of scope for 1.0.0

- GPU or distributed simulation.
- Cloud backends (`Engine`/devices stay local).
- A 1:1 replica of `cirq-google`, `cirq-ionq` or `cirq-aqt` (emulations only).
- A `vqubit.contrib`: Cirq's QASM importer needs `ply`, and there are no third
  parties here.
- Windows arm64 and MSVC (vcraft only supports `win_amd64` for now).
