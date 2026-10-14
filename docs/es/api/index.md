# Referencia de la API

Esta sección se genera automáticamente a partir de los *docstrings* del código.

## Paquete `vqubit`

La API pública es plana: todo se accede como `vqubit.Algo`. La implementación
está organizada en módulos internos, pero solo necesitas recordar el nombre del
paquete.

::: vqubit
    options:
      members: false

## Módulo nativo `vqubit_native`

El trabajo pesado (álgebra lineal, matrices de puertas y simulación) vive en la
extensión nativa `vqubit_native`. No necesitas importarla directamente: el
paquete `vqubit` la envuelve.

::: vqubit_native
