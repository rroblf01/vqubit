# Estabilidad y versiones

vqubit sigue [SemVer 2.0](https://semver.org/lang/es/) y esta página dice
exactamente qué significa aquí, qué se considera API pública y cómo se retira
algo sin romper el código de nadie.

!!! info "Antes del 1.0"
    La versión publicada hoy es la **0.3.0**. La superficie se congela en la
    **1.0.0**; hasta entonces un cambio incompatible puede llegar en una versión
    *menor* (0.4, 0.5…), siempre anotado en el `CHANGELOG`. A partir de la 1.0.0
    rige la regla completa.

## Qué es API pública

- Todo nombre de `vqubit.__all__` (`vq.X`, `vq.Circuit`, `vq.simulate`,
  `vq.circuit_from_qasm`…) y los **miembros públicos** de sus clases
  (`Circuit.moments`, `Gate.num_qubits`…).
- Los submódulos que imitan a Cirq: `vqubit.ops`, `vqubit.sim`, `vqubit.linalg`,
  `vqubit.protocols`, `vqubit.transformers`, `vqubit.work`.
- Los formatos: el JSON de Cirq, el texto OpenQASM 2.0 y el enlace de Quirk.

Lo que empieza por `_` (incluidos `vqubit._gates`, `_linalg`…) es **interno** y
puede cambiar en cualquier momento sin aviso.

## Qué significa cada número

| Cambio | Versión |
| ------ | ------- |
| Quitar un nombre público, o cambiar su comportamiento de forma incompatible | **mayor** |
| Añadir un nombre público, o un comportamiento nuevo compatible | **menor** |
| Corregir un fallo, sin cambiar la superficie | **parche** |

La superficie es comprobable: `tools/api_snapshot.txt` guarda la firma de cada
símbolo público y el CI falla si cambia sin declararlo (`python
tools/check_api_surface.py`). Una corrección que cambia un resultado numérico
incorrecto es un **parche**, aunque el resultado sea distinto: el contrato es lo
que documenta la referencia, no el valor erróneo.

## Cómo se retira algo

1. **Se anuncia** en el `CHANGELOG` (sección `Deprecated`) y en la referencia.
2. **Avisa en tiempo de ejecución** con `DeprecationWarning`, diciendo desde qué
   versión, qué usar en su lugar y cuándo desaparece:

   ```python
   import vqubit as vq

   @vq.deprecated(
       "vqubit.old_helper",
       since="0.9",
       replacement="vqubit.new_helper",
       removed_in="1.0",
   )
   def old_helper(value):
       return new_helper(value)
   ```

   Para algo que no es una función ni una clase (un parámetro renombrado, por
   ejemplo), usa `vq.warn_deprecated(...)` en el punto donde lo detectas.
3. **Se mantiene al menos un ciclo menor** (o un ciclo de parche si la versión
   menor no cambia). Nunca se retira en el mismo ciclo en que se anuncia.
4. **Se quita** solo en una versión mayor, y el `CHANGELOG` lo dice.

Los avisos se lanzan con `stacklevel` correcto, así que apuntan a la línea del
usuario y no a la biblioteca. Los tests usan `pytest.warns(DeprecationWarning)`
para fijar el aviso, de modo que un cambio de mensaje se note.

## Plataformas soportadas

| | |
| --- | --- |
| Python | CPython **3.11 – 3.15** (rueda `abi3` para 3.11) |
| Sistemas | Linux, macOS y Windows, x86_64 y arm64 |
| Dependencias | **ninguna** en tiempo de ejecución: solo la extensión compilada |
| Opcionales | `numpy`, `sympy`, `matplotlib` y `cirq` (grupos `interop` y `vis`) |

Cada versión publicada se prueba en los tres sistemas y las cinco versiones de
Python (`test.yml`), y los números se comparan con Cirq en un trabajo aparte.
Una plataforma que no aparezca en esa matriz no está soportada, aunque pueda
funcionar.

## Cómo lo comprueba el proyecto

- `benchmarks/verify.py`: 619 comprobaciones numéricas contra formas cerradas.
- `benchmarks/parity.py`: 1873 comprobaciones de contrato contra cirq-core 1.7.0.
- `tools/check_api_surface.py`: la superficie pública no cambia sin declararlo.
- `tools/check_api_docs.py`: todo símbolo público tiene docstring.
- `CHANGELOG.md`: qué cambió y en qué versión.

La documentación **versionada** por release (`mike`) está preparada pero
pendiente de un ajuste del repositorio (publicar desde la rama `gh-pages`); hoy
la web publica la rama `main`.
