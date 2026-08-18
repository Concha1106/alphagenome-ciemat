# Experiment exp03 — General Variant Evaluation Pipeline

## Objetivo

Desarrollar un pipeline generalizable y reproducible para la evaluación funcional de variantes genómicas mediante AlphaGenome y aplicarlo al análisis de una estrategia de knock-in terapéutico en SEC23B.

A diferencia de exp02, inicialmente centrado en la exploración y reproducción del caso de exon skipping en DLG1, exp03 consolida un flujo de trabajo reutilizable que integra la ejecución de AlphaGenome, el procesamiento y exportación de sus predicciones y su posterior visualización.

## Estructura del pipeline

El flujo de trabajo se divide en dos módulos principales:

1. `run-variant-pipeline.py`: ejecución de AlphaGenome, procesamiento y exportación estructurada de las predicciones.
2. `visualize-prediction.py`: generación de visualizaciones a partir de predicciones previamente almacenadas, sin realizar nuevas consultas a la API.

### Pipeline de análisis y exportación

El script `run-variant-pipeline.py` permite:

- recibir y validar una variante genómica mediante argumentos de línea de comandos;
- construir el objeto `Variant` y el intervalo de análisis;
- ejecutar `score_variant()` para obtener puntuaciones resumidas del efecto de la variante;
- ejecutar `predict_variant()` para obtener predicciones directas REF y ALT;
- restringir el análisis a una ontología concreta y seleccionar los output types solicitados;
- exportar las predicciones en tablas estructuradas;
- generar tablas resumen para facilitar la priorización e interpretación de las señales;
- serializar el objeto completo devuelto por `predict_variant()` en `prediction.pkl`;
- registrar automáticamente los parámetros de ejecución y archivos generados en `runlog.txt`.

Entre las salidas adicionales generadas por el pipeline se incluyen resúmenes de las predicciones por track y tablas específicas para facilitar la interpretación de outputs relacionados con splicing.

### Optimización de la exportación

Para reducir el uso de memoria durante el procesamiento de predicciones con múltiples tracks, las salidas de tipo TrackData se exportan progresivamente por track en lugar de acumular todas las filas en memoria.

El pipeline permite además seleccionar modalidades concretas mediante `--output-types`, evitando procesar outputs no solicitados.

## Visualización de predicciones

El script `visualize-prediction.py` reutiliza el objeto `prediction.pkl` generado durante el análisis para producir figuras sin repetir llamadas a la API de AlphaGenome.

El visualizador permite representar:

- RNA-seq;
- splice sites;
- splice site usage;
- splice junctions.

La anotación génica se obtiene automáticamente a partir de un archivo GTF compatible con GRCh38. Por defecto se selecciona el transcrito MANE Select del gen analizado, aunque puede indicarse un transcrito concreto.

Las visualizaciones pueden restringirse a regiones genómicas específicas y mantienen una representación común de las predicciones REF, ALT y, cuando corresponde, de la diferencia ALT- REF.

El proceso genera `visualization_runlog.txt`, que registra los parámetros de visualización, el transcrito utilizado, la región representada y las figuras generadas. Para splice junctions se exporta además una tabla con las junctions seleccionadas para su representación.

Las constantes gráficas compartidas se encuentran centralizadas en `config.py`, mientras que distintas funciones auxiliares reutilizables permiten mantener un comportamiento homogéneo entre las representaciones.

## Aplicación a SEC23B

El pipeline se aplicó al análisis de un knock-in terapéutico en `SEC23B` diseñado para el tratamiento de la anemia diseritropoyética congénita tipo II (CDAII).
La estrategia analizada consistió en una inserción mediante HDR en `chr20:18510827`, en la que la secuencia endógena fue reemplazada por un cassette terapéutico de 2.562 pb que incluye el cDNA funcional de `SEC23B` y una señal de poliadenilación bGH.

Debido a que la inserción sustituye 3 pb de la secuencia de referencia, la secuencia endógena situada aguas abajo queda desplazada 2.559 pb en el alelo ALT respecto a REF. Esta diferencia de coordenadas se tuvo en cuenta durante la interpretación de las predicciones generadas por `predict_variant()`.

El análisis se centró principalmente en outputs relacionados con expresión y splicing:

- `RNA_SEQ`
- `SPLICE_SITES`
- `SPLICE_SITE_USAGE`
- `SPLICE_JUNCTIONS`

El análisis se realizó principalmente en las ontologías `CL:0001059` (common myeloid progenitor, CD34-positive) y `CL:0000837` (hematopoietic multipotent progenitor cell), correspondientes a progenitores hematopoyéticos relevantes para el contexto terapéutico estudiado.

De forma complementaria, se evaluaron `CL:0000182` (hepatocyte) y `CL:0000100` (motor neuron) para explorar si el patrón funcional predicho se reproducía en otros contextos celulares. La selección de las ontologías se realizó a partir de los metadatos de tracks disponibles en AlphaGenome.

## Scripts

- `scripts/2026-05-28-exp03-variant-pipeline/run-variant-pipeline.py`
- `scripts/2026-05-28-exp03-variant-pipeline/visualize-prediction.py`

## Resultados y trazabilidad

Los resultados de cada ejecución se almacenan en subdirectorios específicos dentro de:

`results/2026-05-28-exp03-variant-pipeline/`

Cada ejecución del pipeline genera su correspondiente `runlog.txt`, mientras que el módulo de visualización genera `visualization_runlog.txt`.

Los datos externos utilizados en el proyecto se registran en `docs/manifest-data.tsv` y el entorno reproducible se encuentra definido en `environment/environment.yml`.

## Estado

Pipeline completado, validado mediante el caso DLG1 y aplicado al análisis del knock-in en SEC23B.
