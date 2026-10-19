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
| `benchmarks/verify.py` | que los **resultados** coinciden: unitarias, canales, simulación, muestreo, descomposiciones, KAK, Clifford | job `cirq-parity` |
| `benchmarks/parity.py` | que el **contrato** coincide: `**` disponible justo cuando Cirq lo permite, misma unitaria resultante, misma retícula de `isinstance` y mismas clases | job `cirq-parity` |

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
- [ ] Congelar la superficie pública y declarar la política SemVer en detalle,
      incluida la decisión sobre los nombres de las clases de potencia (`vq.CZ`
      debería *ser* un `CZPowGate`, como en Cirq).
- [ ] Política de deprecación con **mecanismo**: un ayudante que emita
      `DeprecationWarning`, ventana mínima, guía de migración en `CHANGELOG.md` y
      un test que compruebe que el aviso se emite de verdad.
- [ ] CI de compatibilidad de API (`griffe`, ya disponible en el entorno de
      desarrollo) que falle ante cambios no intencionados.
- [ ] Auditoría de `__all__` ↔ docstrings ↔ anotaciones; el CI falla si falta
      cualquiera de los tres. **Medido:** 0 nombres sin anotación, **28** sin
      docstring (todos los `Supports*`).

#### 2. Paridad de comportamiento

- [x] (0.3) **`**` en las puertas planas y en las rotaciones.** `vq.CZ ** 0.5`
      ya es una raíz cuadrada de CZ y `vq.rz(0.3) ** 0.5` escala el ángulo a
      `vq.rz(0.15)`, como en Cirq; lo mismo para `CNOT`/`CX`, `SWAP`, `ISWAP`,
      `SQRT_ISWAP`, `CCX`/`TOFFOLI`, `I`, `rx` y `ry`. `FREDKIN`/`CSWAP` siguen
      sin admitir `**`, igual que `cirq.CSWAP`. Ese trabajo destapó además dos
      bugs de inversa en `iSWAP` (ver `CHANGELOG.md`).
- [ ] Que `vq.CZ` **sea** un `CZPowGate`. El comportamiento ya coincide, pero la
      clase no: `isinstance(vq.CZ, vq.CZPowGate)` es `False`, y también
      `vq.CZ == vq.CZPowGate()`, mientras que en Cirq ambas son `True`. Son las
      **12** entradas que quedan en la lista de toleradas con `plan = "0.3"`.
      Cerrarlo exige decidir entre convertir las puertas planas en las clases de
      potencia (rompe `type()` y el diagrama, y obliga a auditar los
      `isinstance` de QASM, Clifford, gatesets y átomos neutros) o dejarlas como
      clases fijas que saltan a la de potencia (documentado y sin tocar el
      resto).
- [ ] Implementar los stubs pendientes: `decompose_clifford_tableau_to_operations`,
      `ResetChannel` multiqubit, `CliffordTableau.kron/factor/reindex` y
      `StabilizerStateChForm.apply_operation` (hoy es una cáscara que no aplica
      puertas), más una revisión de todos los `NotImplementedError`. **Medido:**
      52 en el fuente, de los cuales 48 son métodos abstractos legítimos (Cirq
      los tiene igual) y **4 son huecos reales**.
- [ ] Paridad de simuladores: `simulate_expectation_values`, `simulate_sweep`,
      `simulate_moment_steps`, `sample_expectation_values`, simuladores
      *split/step*, `SimulatorBase`/`SimulatesIntermediateState`,
      `DensityMatrixTrialResult`.
- [ ] Tests diferenciales y de propiedades frente a `cirq-core` (Hypothesis +
      ampliación de `benchmarks/verify.py`) para puertas, canales, KAK,
      transformers, Clifford y muestreo. **Medido:** `verify.py` toca 16 de los
      607 nombres compartidos; **197** no aparecen en ninguna prueba (muchos son
      alias de tipado, pero el hueco es real).
- [ ] Cobertura de casos límite en síntesis de dos qubits y
      `optimize_for_target_gateset`.
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
- [ ] Documentar la política de plataformas soportadas y los ciclos de soporte
      de versiones de Python.

#### 4. Interoperabilidad

- [ ] **Import** de OpenQASM 2.0 (ida y vuelta con `qasm`). El importador de
      `cirq.contrib` depende de `ply`; el nuestro será propio, sin dependencias.
- [ ] Puente opcional `cirq` ↔ `vqubit` (conversión de circuitos) como extra.
- [ ] Completar la cobertura de Quirk (más puertas) y congelar el esquema JSON.
- [ ] Documentar la salida de `numpy`/`sympy` como contrato estable.

#### 5. Documentación

- [x] (0.2.1) Página explícita de limitaciones y diferencias respecto a Cirq
      (`docs/{en,es}/guide/cirq-differences.md`), hoy alineada con la lista de
      toleradas del arnés.
- [ ] Documentación **versionada** por release (`mike`).
- [ ] Guía «Migrar desde Cirq» (la página de diferencias es su materia prima).
- [ ] Tutoriales y galería de ejemplos (Grover, estabilizadores/QEC, variacional).
- [ ] Referencia de API completa con ejemplo por símbolo (mkdocstrings).

#### 6. Comunidad y gobernanza

- [x] (0.2.0) `CONTRIBUTING.md`.
- [ ] `CODE_OF_CONDUCT.md` y `SECURITY.md` (este último, con política de reporte
      de vulnerabilidades en el stub nativo).
- [ ] Plantillas de issue y PR; `CITATION.cff`.
- [ ] Notas de release automáticas (`towncrier`) y disciplina de `CHANGELOG.md`.
- [ ] Autoría y agradecimientos; etiquetas de «good first issue».

#### 7. Calidad y tooling (apoyo)

- [x] (0.2.0) Ejecutar la suite en **Windows** dentro del CI: hoy son 3 SO × 5
      versiones de Python (3.11–3.15), verde.
- [ ] Objetivo de cobertura ≥95 % con umbral en CI. **Medido hoy: 77 %**, con los
      peores módulos en `_vis.py` (13 %), `_experiments.py` (27 %), `_work.py`
      (31 %) y `_clifford.py` (41 %) — este último importa, porque el simulador
      de estabilizadores es una de las banderas del proyecto.
- [ ] `stubtest` de `mypy` sobre `vqubit_native.pyi`. **Bloqueado por vcraft**
      (`UPSTREAM.md` §1 y §2: el `.pyi` genera los enums como anotaciones y
      escribe `__init__` donde el constructor vive en `tp_new`). Con el enum
      parcheado a mano, hoy quedan 6 errores.
- [ ] Tests de propiedades (`hypothesis`) y *fuzzers* de los parsers
      QASM/Quirk/JSON.
- [ ] `mypy`/`pyright` en modo estricto junto a `ty`; escaneo semanal de
      dependencias y seguridad.

#### 8. Rendimiento (apoyo)

- [ ] Cerrar la brecha de simulación a 20 qubits: es **DRAM/TLB-bound**, no
      aritmética (los kernels van a ~25 GB/s, el ancho de banda de la máquina), y
      `-march` no la mueve. Un vector de 20 qubits son 16 MiB: con páginas de
      4 KiB cada pasada de puerta toca ~8192 traducciones y desborda el TLB L2
      (~2k entradas). Pedir a vcraft **páginas gigantes** (`mmap` +
      `madvise(MADV_HUGEPAGE)` en `buffer_new_*`; `UPSTREAM.md` §4) y, a más
      largo plazo, **dispatch en runtime** (*function multiversioning*) para los
      kernels calientes, que daría SIMD con una rueda portable.
- [ ] `-march=…` queda como **opt-in documentado** para un build propio; no se
      hornea en la rueda publicada (portabilidad y no divergir del artefacto
      testeado).
- [ ] Benchmarks en CI con umbrales de regresión.
- [ ] Evaluar más SIMD/OpenMP en los kernels V; densidad a mayor anchura.

### Hitos

| Versión | Tema | Contenido clave |
| ------- | ---- | --------------- |
| 0.2 | Plataformas | vcraft 1.4.5, Windows `win_amd64`, sdist, `--strip`, API nativa moderna |
| 0.3 | Paridad | `**` en las puertas planas y las rotaciones, métodos de simulador, stubs restantes |
| 0.4 | Interop | import QASM propio, puente `cirq` |
| 0.5 | Calidad | cobertura ≥95 %, Hypothesis, `stubtest`, fuzzers |
| 0.6 | Rendimiento | 20 qubits (páginas gigantes), benchmarks en CI, SIMD/OpenMP |
| 0.7 | Distribución | extras, procedencia (SLSA), builds reproducibles |
| 0.8 | Documentación | docs versionadas, migración desde Cirq, galería |
| 0.9 | Congelación | gobernanza, api-diff, SemVer, *release candidate* |
| 1.0 | Estabilidad | garantía SemVer y soporte a largo plazo |

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
| 0.9 | `griffe check` no reporta cambios sin declarar; el test de deprecación pasa |
| 1.0 | la definición de terminado, entera |

### Definición de terminado (1.0.0)

- [ ] Paridad de nombres 607/607, con `benchmarks/parity_tolerated.toml`
      conteniendo **solo** entradas `permanent` y explicadas en la página de
      diferencias.
- [ ] Cero `NotImplementedError` reales (los abstractos, como en Cirq, se quedan).
- [ ] Cobertura ≥95 % con umbral en CI, y el umbral sube, no baja.
- [ ] `stubtest` limpio sobre `vqubit_native.pyi`.
- [ ] Extras publicados, documentados y probados en un entorno limpio.
- [ ] Documentación versionada por release, con la guía de migración desde Cirq.
- [ ] Política SemVer y de deprecación con mecanismo y test.
- [ ] `api-diff` en CI sobre la superficie pública.
- [ ] Gobernanza completa: COC, SECURITY, plantillas, CITATION, notas de release.
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
| `benchmarks/verify.py` | that the **numbers** agree: unitaries, channels, simulation, sampling, decompositions, KAK, Clifford | `cirq-parity` job |
| `benchmarks/parity.py` | that the **contract** agrees: `**` available exactly when Cirq allows it, same resulting unitary, same `isinstance` lattice and same classes | `cirq-parity` job |

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
- [ ] Freeze the public surface and spell out the SemVer policy, including the
      decision on the power-class names (`vq.CZ` should *be* a `CZPowGate`, as in
      Cirq).
- [ ] Deprecation policy with a **mechanism**: a helper that emits
      `DeprecationWarning`, a minimum window, migration notes in `CHANGELOG.md`
      and a test that the warning really fires.
- [ ] API-compatibility CI (`griffe`, already available in the dev environment)
      failing on unintended changes.
- [ ] Audit `__all__` ↔ docstrings ↔ annotations; CI fails if any is missing.
      **Measured:** 0 names without annotations, **28** without a docstring (all
      the `Supports*`).

#### 2. Behaviour parity

- [x] (0.3) **`**` on the plain gates and the rotations.** `vq.CZ ** 0.5` is now
      a square root of CZ and `vq.rz(0.3) ** 0.5` scales the angle to
      `vq.rz(0.15)`, as in Cirq; likewise for `CNOT`/`CX`, `SWAP`, `ISWAP`,
      `SQRT_ISWAP`, `CCX`/`TOFFOLI`, `I`, `rx` and `ry`. `FREDKIN`/`CSWAP` still
      refuse `**`, exactly like `cirq.CSWAP`. That work also uncovered two
      inverse bugs in `iSWAP` (see `CHANGELOG.md`).
- [ ] Make `vq.CZ` **be** a `CZPowGate`. The behaviour already matches, the class
      does not: `isinstance(vq.CZ, vq.CZPowGate)` is `False`, and so is
      `vq.CZ == vq.CZPowGate()`, where Cirq has both `True`. These are the **12**
      entries still carrying `plan = "0.3"` in the tolerated list. Closing it
      means choosing between folding the plain gates into the power classes
      (breaks `type()` and the diagram, and forces an audit of the `isinstance`
      checks in QASM, Clifford, gatesets and neutral atoms) or keeping them as
      fixed classes that bridge to the power gate (documented, nothing else
      touched).
- [ ] Implement the remaining stubs: `decompose_clifford_tableau_to_operations`,
      multi-qubit `ResetChannel`, `CliffordTableau.kron/factor/reindex` and
      `StabilizerStateChForm.apply_operation` (a shell that applies no gates
      today), plus a review of every `NotImplementedError`. **Measured:** 52 in
      the source, of which 48 are legitimate abstract methods (Cirq has them too)
      and **4 are real gaps**.
- [ ] Simulator parity: `simulate_expectation_values`, `simulate_sweep`,
      `simulate_moment_steps`, `sample_expectation_values`, split/step
      simulators, `SimulatorBase`/`SimulatesIntermediateState`,
      `DensityMatrixTrialResult`.
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

- [ ] OpenQASM 2.0 **import** (round-trip with `qasm`). Cirq's importer in
      `cirq.contrib` depends on `ply`; ours will be its own, dependency-free.
- [ ] Optional `cirq` ↔ `vqubit` bridge (circuit conversion) as an extra.
- [ ] Complete Quirk coverage (more gates) and freeze the JSON schema.
- [ ] Document `numpy`/`sympy` output as a stable contract.

#### 5. Documentation

- [x] (0.2.1) An explicit page of limitations and differences from Cirq
      (`docs/{en,es}/guide/cirq-differences.md`), now kept in step with the
      harness's tolerated list.
- [ ] **Versioned** docs per release (`mike`).
- [ ] A "Migrating from Cirq" guide (the differences page is its raw material).
- [ ] Tutorials and an examples gallery (Grover, stabilizers/QEC, variational).
- [ ] Complete API reference with an example per symbol (mkdocstrings).

#### 6. Community and governance

- [x] (0.2.0) `CONTRIBUTING.md`.
- [ ] `CODE_OF_CONDUCT.md` and `SECURITY.md` (the latter with a vulnerability
      reporting policy for the native stub).
- [ ] Issue and PR templates; `CITATION.cff`.
- [ ] Automated release notes (`towncrier`) and `CHANGELOG.md` discipline.
- [ ] Authorship and acknowledgements; "good first issue" labels.

#### 7. Quality and tooling (supporting)

- [x] (0.2.0) Run the suite on **Windows** in CI: today 3 OSes × 5 Python
      versions (3.11–3.15), green.
- [ ] Coverage target ≥95 % enforced in CI. **Measured today: 77 %**, with the
      worst modules at `_vis.py` (13 %), `_experiments.py` (27 %), `_work.py`
      (31 %) and `_clifford.py` (41 %) — the last one matters, because the
      stabilizer simulator is a flagship feature.
- [ ] `mypy` `stubtest` over `vqubit_native.pyi`. **Blocked by vcraft**
      (`UPSTREAM.md` §1 and §2: the `.pyi` writes enums as annotations and emits
      `__init__` where the constructor lives in `tp_new`). With the enum patched
      by hand, 6 errors remain today.
- [ ] Property tests (`hypothesis`) and *fuzzers* for the QASM/Quirk/JSON
      parsers.
- [ ] `mypy`/`pyright` in strict mode alongside `ty`; weekly dependency and
      security scanning.

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
| 0.3 | Parity | `**` on the plain gates and the rotations, simulator methods, remaining stubs |
| 0.4 | Interop | own QASM import, `cirq` bridge |
| 0.5 | Quality | coverage ≥95 %, Hypothesis, `stubtest`, fuzzers |
| 0.6 | Performance | 20 qubits (huge pages), benchmarks in CI, SIMD/OpenMP |
| 0.7 | Distribution | extras, provenance (SLSA), reproducible builds |
| 0.8 | Documentation | versioned docs, Cirq migration, gallery |
| 0.9 | Freeze | governance, api-diff, SemVer, release candidate |
| 1.0 | Stability | SemVer guarantee and long-term support |

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
| 0.9 | `griffe check` reports no undeclared changes; the deprecation test passes |
| 1.0 | the whole definition of done |

### Definition of done (1.0.0)

- [ ] Name parity 607/607, with `benchmarks/parity_tolerated.toml` holding
      **only** `permanent`, explained entries.
- [ ] Zero real `NotImplementedError`s (the abstract ones, as in Cirq, stay).
- [ ] Coverage ≥95 % with a CI threshold that rises, never falls.
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
