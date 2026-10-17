# vqubit roadmap

[Español](#roadmap-hacia-la-100) · [English](#roadmap-to-100)

---

## Roadmap hacia la 1.0.0

Esta hoja de ruta traza el camino desde la **0.1.0** (primera publicación en
PyPI) hasta la **1.0.0**, el momento en el que la API pública queda congelada y
el proyecto ofrece garantías de estabilidad. Es un documento vivo: se actualiza
a medida que avanzan los hitos.

> **Actualización (vcraft 1.4.0, 2026-10-16).** La nueva versión de vcraft
> desbloquea **Windows (`win_amd64`)** —el hito que estaba bloqueado— y confirma
> que `vcraft sdist`, `vcraft publish` y `--strip` ya están disponibles. Además
> cambia el **commit de V fijado** (`6851aaf…`, ≠ el `36be926…` = V 0.5.2 que
> usábamos), así que la actualización exige `vcraft toolchain install`. La
> **0.2.0** recoge Windows, sdist y strip.

### Principios

1. **Cero dependencias.** El wheel publicado nunca arrastra terceros; `numpy`,
   `sympy` y `matplotlib` son extras opcionales.
2. **Compatibilidad con Cirq.** La API refleja Cirq en nombres *y* en
   comportamiento.
3. **Velocidad nativa.** El cálculo pesado vive en V (kernels compilados).
4. **Un wheel por plataforma.** Un único wheel `abi3` cubre CPython 3.11–3.15.
5. **Cambios sin sorpresas.** Nada rompe sin aviso previo y deprecación.

### Punto de partida (0.1.0 → 0.2.0)

- 650 nombres públicos exportados; paridad de **nombres** con `cirq-core`.
- Simuladores de vector de estado, matriz de densidad y estabilizadores;
  canales de ruido; optimizadores; transformers; devices; QASM/Quirk/JSON.
- Tipado completo (`py.typed`, todas las funciones anotadas) y stub nativo.
- 548 tests y CI verde; cobertura de nombres verificada contra Cirq.
- 0.1.0: Linux y macOS (x86-64 y ARM), sin dependencias, solo wheels.
- 0.2.0 (2026-10-16): **Windows (`win_amd64`)**, **sdist** y wheels *stripped*,
  sobre **vcraft 1.4.1**.

### Objetivo de la 1.0.0

Una librería **estable, completa y documentada**: paridad de comportamiento (no
solo de nombres), distribución en todas las plataformas objetivo, política de
versionado y deprecación explícita, y gobernanza abierta.

### Ejes

#### 1. Estabilidad y contrato de API

- [ ] Congelar la superficie pública y declarar la política SemVer en detalle.
- [ ] Política de deprecación: `DeprecationWarning`, ventana mínima y guía de
      migración en `CHANGELOG.md`.
- [ ] CI de compatibilidad de API (`api-diff`/griffe) que falle ante cambios no
      intencionados.
- [ ] Auditoría de `__all__` ↔ docstrings ↔ anotaciones; el CI falla si falta
      cualquiera de los tres para un nombre público.
- [ ] Resolver la historia de submódulos: que `import vqubit.work` y
      `import vqubit.vis` funcionen (como en Cirq) o documentar formalmente que
      la API es solo plana.

#### 2. Paridad de comportamiento

- [ ] Implementar los stubs pendientes:
      `decompose_clifford_tableau_to_operations`, `ResetChannel` multiqubit y
      revisión de todos los `NotImplementedError`.
- [ ] Paridad de simuladores: `simulate_expectation_values`, simuladores
      *split/step*, `SimulatorBase`/`SimulatesIntermediateState`,
      `DensityMatrixTrialResult`.
- [ ] Tests diferenciales y de propiedades frente a `cirq-core` (Hypothesis +
      ampliación de `benchmarks/verify.py`) para puertas, canales, KAK,
      transformers, Clifford y muestreo.
- [ ] Cobertura de casos límite en síntesis de dos qubits y
      `optimize_for_target_gateset`.
- [x] (0.2.0) Modernizar la API nativa: `jacobi_eig`, `simultaneous_diag` y
      `kron_factor` devuelven *structs* con nombre. Un `enum` para
      `NativeTableau.apply` sigue pendiente (comprobado también en vcraft 1.4.1):
      el glue pasa un `int` a un parámetro enum, así que no compila.

#### 3. Empaquetado y distribución

- [x] Publicar **sdist** junto a los wheels (`vcraft sdist`) — 0.2.0.
- [x] Wheels de **Windows** (`win_amd64`) — 0.2.0, ya soportado por vcraft 1.4.0.
- [x] Wheels *stripped* en release (workaround: `--strip` es no-op en vcraft 1.4.0) — 0.2.0.
- [ ] Corregir classifiers si algún día cambia el soporte (hoy `Microsoft ::
      Windows` ya es correcto).
- [ ] Extras publicados y documentados: `vqubit[vis]`, `vqubit[interop]`.
- [ ] Builds reproducibles, verificación por hash del `RECORD` y atestación de
      procedencia (SLSA) en PyPI.
- [ ] Documentar la política de plataformas soportadas y los ciclos de soporte
      de versiones de Python.

#### 4. Interoperabilidad

- [ ] **Import** de OpenQASM 2.0 (ida y vuelta con `qasm`).
- [ ] Puente opcional `cirq` ↔ `vqubit` (conversión de circuitos) como extra.
- [ ] Completar la cobertura de Quirk (más puertas) y congelar el esquema JSON.
- [ ] Documentar la salida de `numpy`/`sympy` como contrato estable.

#### 5. Documentación

- [ ] Documentación **versionada** por release (`mike`).
- [ ] Guía «Migrar desde Cirq».
- [ ] Tutoriales y galería de ejemplos (Grover, estabilizadores/QEC, variacional).
- [ ] Referencia de API completa con ejemplo por símbolo (mkdocstrings).
- [ ] Página explícita de limitaciones y diferencias respecto a Cirq.

#### 6. Comunidad y gobernanza

- [ ] `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` y `SECURITY.md`.
- [ ] Plantillas de issue y PR; `CITATION.cff`.
- [ ] Notas de release automáticas (`towncrier`) y disciplina de `CHANGELOG.md`.
- [ ] Autoría y agradecimientos; etiquetas de «good first issue».

#### 7. Calidad y tooling (apoyo)

- [ ] Objetivo de cobertura (p. ej. ≥95 %) con umbral en CI.
- [ ] Tests de propiedades (`hypothesis`) y *fuzzers* de los parsers
      QASM/Quirk/JSON.
- [ ] `stubtest` de `mypy` sobre `vqubit_native.pyi` para mantener el stub fiel.
- [ ] `mypy`/`pyright` en modo estricto junto a `ty`; escaneo semanal de
      dependencias y seguridad.
- [ ] Ejecutar la suite en **Windows** dentro del CI (0.2.0) y corregir
      supuestos POSIX.

#### 8. Rendimiento (apoyo)

- [ ] Cerrar la brecha de simulación a 20 qubits (actualmente a la par).
- [ ] Benchmarks en CI con umbrales de regresión.
- [ ] Evaluar SIMD/OpenMP en los kernels V; densidad a mayor anchura.

### Hitos

| Versión | Tema | Contenido clave |
| ------- | ---- | --------------- |
| 0.2 | Plataformas | vcraft 1.4.0, Windows `win_amd64`, sdist, `--strip`, API nativa moderna |
| 0.3 | Interop | import QASM, puente `cirq` |
| 0.4 | Calidad | cobertura, Hypothesis, `stubtest`, fuzzers |
| 0.5 | Rendimiento | 20 qubits, benchmarks en CI, SIMD/OpenMP |
| 0.6 | Distribución | extras, procedencia (SLSA), builds reproducibles |
| 0.7 | Documentación | docs versionadas, migración desde Cirq, galería |
| 0.8 | Gobernanza | CONTRIBUTING, plantillas, SECURITY, CITATION |
| 0.9 | Congelación | API congelada, *release candidate* |
| 1.0 | Estabilidad | garantía SemVer y soporte a largo plazo |

### Fuera de alcance de la 1.0.0

- Simulación por GPU o distribuida.
- Backends en la nube (los `Engine`/devices siguen siendo locales).
- Réplica 1:1 de `cirq-google`, `cirq-ionq` o `cirq-aqt` (solo emulaciones).
- Windows arm64 y MSVC (vcraft solo soporta `win_amd64` por ahora).

---

## Roadmap to 1.0.0

This roadmap charts the path from **0.1.0** (the first PyPI release) to
**1.0.0**, when the public API is frozen and the project offers stability
guarantees. It is a living document, updated as milestones land.

> **Update (vcraft 1.4.0, 2026-10-16).** The new vcraft release unblocks
> **Windows (`win_amd64`)** — the milestone that was blocked — and confirms that
> `vcraft sdist`, `vcraft publish` and `--strip` are already available. It also
> changes the **pinned V commit** (`6851aaf…`, ≠ the `36be926…` = V 0.5.2 we
> used), so the upgrade requires `vcraft toolchain install`. **0.2.0** delivers
> Windows, sdist and strip.

### Principles

1. **Zero dependencies.** The published wheel never pulls in third parties;
   `numpy`, `sympy` and `matplotlib` are optional extras.
2. **Cirq compatibility.** The API mirrors Cirq in names *and* behaviour.
3. **Native speed.** The heavy lifting lives in V (compiled kernels).
4. **One wheel per platform.** A single `abi3` wheel covers CPython 3.11–3.15.
5. **No surprise breakage.** Nothing breaks without notice and deprecation.

### Starting point (0.1.0 → 0.2.0)

- 650 public names exported; **name** parity with `cirq-core`.
- State-vector, density-matrix and stabilizer simulators; noise channels;
  optimizers; transformers; devices; QASM/Quirk/JSON.
- Complete typing (`py.typed`, every function annotated) and a native stub.
- 548 tests and green CI; name coverage verified against Cirq.
- 0.1.0: Linux and macOS (x86-64 and ARM), no dependencies, wheels only.
- 0.2.0 (2026-10-16): **Windows (`win_amd64`)**, **sdist** and stripped wheels,
  on **vcraft 1.4.1**.

### Goal for 1.0.0

A **stable, complete and well documented** library: behaviour parity (not just
names), distribution across every target platform, an explicit versioning and
deprecation policy, and open governance.

### Workstreams

#### 1. API stability and contract

- [ ] Freeze the public surface and spell out the SemVer policy.
- [ ] Deprecation policy: `DeprecationWarning`, minimum window and migration
      notes in `CHANGELOG.md`.
- [ ] API-compatibility CI (`api-diff`/griffe) failing on unintended changes.
- [ ] Audit `__all__` ↔ docstrings ↔ annotations; CI fails if any is missing
      for a public name.
- [ ] Settle the submodule story: make `import vqubit.work` and
      `import vqubit.vis` work (as in Cirq) or formally document that the API is
      flat-only.

#### 2. Behaviour parity

- [ ] Implement the remaining stubs: `decompose_clifford_tableau_to_operations`,
      multi-qubit `ResetChannel` and a review of every `NotImplementedError`.
- [ ] Simulator parity: `simulate_expectation_values`, split/step simulators,
      `SimulatorBase`/`SimulatesIntermediateState`, `DensityMatrixTrialResult`.
- [ ] Differential and property tests against `cirq-core` (Hypothesis + extend
      `benchmarks/verify.py`) for gates, channels, KAK, transformers, Clifford
      and sampling.
- [ ] Edge-case coverage in two-qubit synthesis and
      `optimize_for_target_gateset`.
- [x] (0.2.0) Modernise the native API: `jacobi_eig`, `simultaneous_diag` and
      `kron_factor` return named *structs*. An `enum` for `NativeTableau.apply`
      is still pending (checked in vcraft 1.4.1 too): the glue passes an `int` to
      an `enum` parameter, so it does not compile.

#### 3. Packaging and distribution

- [x] Publish an **sdist** alongside the wheels (`vcraft sdist`) — 0.2.0.
- [x] **Windows** (`win_amd64`) wheels — 0.2.0, now supported by vcraft 1.4.0.
- [x] Stripped release wheels (workaround: `--strip` is a no-op in vcraft 1.4.0) — 0.2.0.
- [ ] Keep classifiers accurate (`Microsoft :: Windows` is correct now).
- [ ] Publish and document extras: `vqubit[vis]`, `vqubit[interop]`.
- [ ] Reproducible builds, `RECORD` hash verification and provenance
      attestation (SLSA) on PyPI.
- [ ] Document the supported-platform policy and Python support windows.

#### 4. Interoperability

- [ ] OpenQASM 2.0 **import** (round-trip with `qasm`).
- [ ] Optional `cirq` ↔ `vqubit` bridge (circuit conversion) as an extra.
- [ ] Complete Quirk coverage (more gates) and freeze the JSON schema.
- [ ] Document `numpy`/`sympy` output as a stable contract.

#### 5. Documentation

- [ ] **Versioned** docs per release (`mike`).
- [ ] A "Migrating from Cirq" guide.
- [ ] Tutorials and an examples gallery (Grover, stabilizers/QEC, variational).
- [ ] Complete API reference with an example per symbol (mkdocstrings).
- [ ] An explicit page of limitations and differences from Cirq.

#### 6. Community and governance

- [ ] `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` and `SECURITY.md`.
- [ ] Issue and PR templates; `CITATION.cff`.
- [ ] Automated release notes (`towncrier`) and `CHANGELOG.md` discipline.
- [ ] Authorship and acknowledgements; "good first issue" labels.

#### 7. Quality and tooling (supporting)

- [ ] Coverage target (e.g. ≥95 %) enforced in CI.
- [ ] Property tests (`hypothesis`) and *fuzzers* for the QASM/Quirk/JSON
      parsers.
- [ ] `mypy` `stubtest` over `vqubit_native.pyi` to keep the stub honest.
- [ ] `mypy`/`pyright` in strict mode alongside `ty`; weekly dependency and
      security scanning.
- [ ] Run the suite on **Windows** in CI (0.2.0) and fix POSIX assumptions.

#### 8. Performance (supporting)

- [ ] Close the 20-qubit simulation gap (currently at parity).
- [ ] Benchmarks in CI with regression thresholds.
- [ ] Evaluate SIMD/OpenMP in the V kernels; density matrices at larger widths.

### Milestones

| Version | Theme | Key content |
| ------- | ----- | ----------- |
| 0.2 | Platforms | vcraft 1.4.0, Windows `win_amd64`, sdist, `--strip`, modern native API |
| 0.3 | Interop | QASM import, `cirq` bridge |
| 0.4 | Quality | coverage, Hypothesis, `stubtest`, fuzzers |
| 0.5 | Performance | 20 qubits, benchmarks in CI, SIMD/OpenMP |
| 0.6 | Distribution | extras, provenance (SLSA), reproducible builds |
| 0.7 | Documentation | versioned docs, Cirq migration, gallery |
| 0.8 | Governance | CONTRIBUTING, templates, SECURITY, CITATION |
| 0.9 | Freeze | API frozen, release candidate |
| 1.0 | Stability | SemVer guarantee and long-term support |

### Out of scope for 1.0.0

- GPU or distributed simulation.
- Cloud backends (`Engine`/devices stay local).
- A 1:1 replica of `cirq-google`, `cirq-ionq` or `cirq-aqt` (emulations only).
- Windows arm64 and MSVC (vcraft only supports `win_amd64` for now).
