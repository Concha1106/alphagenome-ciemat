# Experiment exp02 — Exon Skipping (DLG1)

## Objetivo

Validar el flujo de análisis mediante la reproducción del caso de exon skipping en DLG1 descrito por Avsec et al. (2026), utilizando tanto las métricas de efecto de `score_variant()` como las predicciones directas REF/ALT obtenidas mediante `predict_variant()`.

## Datos de entrada

- Variante: `chr3:197081044 TACTC>T`
- Gen: `DLG1`
- Genoma de referencia: GRCh38
- Ontología: `UBERON:0007610` (Artery Tibial)
- Transcrito MANE Select: `ENST00000667157`
- Anotación génica: Ensembl GRCh38, release 115

El caso de referencia corresponde al evento de exon skipping presentado en la Figura 3b de Avsec et al. (2026).

## Análisis

El experimento combina dos aproximaciones de AlphaGenome:

1. `score_variant()`, utilizado para obtener métricas resumidas del efecto de la variante sobre el splicing.
2. `predict_variant()`, utilizado para comparar directamente las predicciones REF y ALT de:
   - RNA-seq
   - splice junctions
   - splice sites
   - splice site usage

Las señales relacionadas con splicing se analizaron en la hebra negativa, correspondiente a DLG1.

La anotación del transcrito MANE Select se utilizó para contextualizar las splice junctions predichas respecto a la estructura exónica de DLG1.

## Scripts

- `scripts/2026-04-10-exp02-exon-skipping/exp02-score_exon_skipping.py`
- `scripts/2026-04-10-exp02-exon-skipping/exp02-ref_vs_alt_prediction.py`
- `scripts/2026-04-10-exp02-exon-skipping/plot-fig3b-junctions.py`
- `scripts/2026-04-10-exp02-exon-skipping/plot-fig3b-rna-seq.py`
- `scripts/2026-04-10-exp02-exon-skipping/plot-fig3b-splice-sites-and-usage.py`

## Resultados principales

El análisis mediante `score_variant()` mostró un efecto relevante de la variante sobre el splicing, con una puntuación combinada de aproximadamente 3.61.

La comparación REF/ALT mediante `predict_variant()` reprodujo el patrón esperado de exon skipping:

- pérdida de las junctions que conectan los exones 18–17 y 17–16;
- aparición en ALT de una junction que conecta directamente los exones 18–16;
- reducción de la señal de RNA-seq correspondiente al exón 17;
- cambios concordantes en splice sites y splice site usage.

En conjunto, las predicciones reproducen el evento de omisión del exón 17 descrito por Avsec et al. (2026).

## Visualización

Los principales outputs se representan mediante scripts específicos para RNA-seq, splice junctions, splice sites y splice site usage.

Las visualizaciones utilizan una región genómica común y la anotación del transcrito MANE Select para facilitar la comparación entre REF y ALT y la interpretación estructural del evento.

## Archivos asociados

**Datos intermedios:** `interm-data/2026-04-10-exp02-exon-skipping/`

**Resultados:** `results/2026-04-10-exp02-exon-skipping/`

El archivo `runlog.txt` registra los parámetros principales de la ejecución.

## Estado

Experimento completado y utilizado como caso de validación del flujo de análisis.
