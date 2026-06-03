Experiment exp03 — General Variant Evaluation Pipeline

## Objetivo

Desarrollar un pipeline general y automatizable para la evaluación exploratoria de variantes genéticas mediante AlphaGenome.

A diferencia del experimento exp02, centrado en la replicación de un caso concreto de exon skipping en DLG1, este experimento busca construir una arquitectura generalizada que permita analizar variantes novedosas sin asumir previamente el mecanismo biológico implicado.

De este modo, se pretende desarrollar un flujo de trabajo reproducible capaz de:

- recibir variantes genómicas como entrada;
- realizar un primer cribado global mediante score_variant();
- priorizar señales biológicas potencialmente relevantes;
- generar predicciones directas REF vs ALT mediante predict_variant();
- exportar resultados estructurados para su posterior interpretación.

Para asegurar la trazabilidad y reproducibilidad el script se desarrollará siguiendo unos criterios similares a los utilizados en los experimentos previos:
registro de parámetros usados, variante evaluada,outputs solicitados, archivos generados, etc. Esto quedará registrado en archivos runlog.txt de este proyecto, el manifes-data.tsv y Git.

## Estructura pipeline

El pipeline se desarrollará en diferentes fases que modularán lo siguiente:

1. Lectura de argumentos obligatorios y opcionales
2. Validación de variantes
3. Construcción de intervalos;
4. Generación automática de directorios y runlogs
5. Análisis mediante score_variant()
6. Análisis mediante predict_variant()
7. Exportación automática de tablas resumen
8. Priorización automática de señales biológicas (mejora a desarrollar en futuro).

### Inputs
- position  
- reference (REF)
- alternate 
Argumentos opcionales previstos
- interval size

### Outputs previstos

El pipeline generará:

- tablas completas de scores
- tablas priorizadas
- predicciones REF vs ALT

### Estado

Pipeline funcional en versión inicial.

Actualmente el script permite:

- recibir una variante mediante argumentos de línea de comandos;
- construir el objeto `Variant` de AlphaGenome;
- construir el intervalo de entrada usando `variant.reference_interval.resize()`;
- ejecutar `score_variant()`;
- exportar la tabla completa `score_variant_all.tsv`;
- ejecutar `predict_variant()` para una ontología concreta indicada mediante `--ontology-curie`;
- exportar predicciones REF, ALT y delta para salidas tipo TrackData;
- exportar predicciones REF, ALT y delta para splice junctions;
- gestionar contact maps, exportándolos solo cuando existan tracks disponibles;
- generar tablas top positivas y negativas por `output_type`;
- registrar parámetros de ejecución en `runlog.txt`.

### Actualización: optimización de exportación predict_variant

Se identificó un cuello de botella de memoria durante la exportación local de predicciones generadas por `predict_variant()`, especialmente al trabajar con ontologías con múltiples tracks y ventanas grandes.

Para reducir el uso de RAM, la exportación de salidas tipo TrackData se modificó para escribir los resultados progresivamente por track, en lugar de acumular todas las filas en memoria antes de generar el archivo final.

Además, se mejoró el manejo de ejecuciones parciales mediante `--output-types`, permitiendo solicitar únicamente modalidades concretas, como `CAGE` o `SPLICE_SITES`, sin que el exportador falle al encontrar outputs no solicitados.


