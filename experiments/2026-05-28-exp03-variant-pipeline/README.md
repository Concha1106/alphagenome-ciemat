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
- Argumentos opcionales previstos
- interval size

### Outputs previstos

El pipeline generará:

- tablas completas de scores
- tablas priorizadas
- predicciones REF vs ALT

### Estado

Pipeline funcional.

Actualmente el script permite:

- recibir una variante mediante argumentos de línea de comandos
- construir el objeto `Variant` de AlphaGenome
- construir el intervalo de entrada usando `variant.reference_interval.resize()`
- ejecutar `score_variant()`
- exportar la tabla completa `score_variant_all.tsv`
- ejecutar `predict_variant()` para una ontología concreta indicada mediante `--ontology-curie`
- exportar predicciones REF, ALT y delta para salidas tipo TrackData
- exportar predicciones REF, ALT y delta para splice junctions
- gestionar contact maps, exportándolos solo cuando existan tracks disponibles
- generar tablas top positivas y negativas por `output_type`
- registrar parámetros de ejecución en `runlog.txt`
- guardar el objeto completo predict_variant() como prediction.pkl, permitiendo reutilizar las predicciones sin repetir llamadas a la API
- registrar prediction.pkl en runlog.txt como archivo generado
- desarrollar un script independiente de visualización (visualize-prediction.py) basado en el objeto serializado
- generar una primera figura local de RNA-seq con: señal REF vs ALT, delta ALT-REF, anotación MANE Select de SEC23B, posición y coordenada del KI y región de visualización parametrizable

### Optimización de exportación predict_variant

Se identificó un cuello de botella de memoria durante la exportación local de predicciones generadas por `predict_variant()`, especialmente al trabajar con ontologías con múltiples tracks y ventanas grandes.

Para reducir el uso de RAM, la exportación de salidas tipo TrackData se modificó para escribir los resultados progresivamente por track, en lugar de acumular todas las filas en memoria antes de generar el archivo final.

Además, se mejoró el manejo de ejecuciones parciales mediante `--output-types`, permitiendo solicitar únicamente modalidades concretas, como `CAGE` o `SPLICE_SITES`, sin que el exportador falle al encontrar outputs no solicitados.

### Mejoras añadidas

- calcular y exportar `score_variant_merged_splicing.tsv`, con una métrica integrada de splicing basada en `SPLICE_SITES`, `SPLICE_SITE_USAGE` y `SPLICE_JUNCTIONS`;
- exportar `splice_sites` en formato wide para facilitar su visualización en Excel;
- generar vistas locales alrededor de la variante para facilitar la interpretación regional;
- generar `predict_variant_track_summary.tsv`, una tabla resumen por output type y track con máximos, medias y número de posiciones con delta relevante.

## Visualización de predicciones

Se desarrolló un módulo independiente de visualización (visualize-prediction.py) con el objetivo de separar el análisis computacional de la generación de figuras.

Mientras que el pipeline principal (run-variant-pipeline.py) ejecuta score_variant() y predict_variant(), exporta tablas estructuradas y serializa el objeto completo de predicción (prediction.pkl), el script de visualización reutiliza dicho objeto para generar figuras sin realizar nuevas llamadas a la API de AlphaGenome. Esta estrategia mejora la reproducibilidad del análisis, reduce el consumo de recursos y facilita la iteración sobre las representaciones gráficas.

Actualmente, el visualizador permite generar figuras para los principales outputs relacionados con splicing y expresión génica (RNA-seq, splice_sites, splice_site_usage y splice_junctions).

Como mejora respecto a versiones anteriores, la anotación génica deja de estar codificada específicamente para SEC23B y pasa a obtenerse automáticamente a partir de un archivo GTF compatible con GRCh38. El script selecciona por defecto el transcrito MANE Select, aunque permite indicar un transcrito concreto cuando sea necesario. De este modo, el visualizador puede reutilizarse para cualquier gen presente en el archivo de anotación.

Además, el proceso genera automáticamente un visualization_runlog.txt, donde quedan registrados los parámetros de visualización, el transcrito utilizado, la región representada y las figuras generadas, favoreciendo la trazabilidad del análisis. En el caso de splice_junctions, también se exporta una tabla con las junctions finalmente representadas en la figura, facilitando la revisión e interpretación posterior.
