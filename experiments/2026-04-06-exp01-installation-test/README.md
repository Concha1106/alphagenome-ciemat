# Experiment exp01 — Installation Test

## Objetivo

Verificar la correcta instalación y el funcionamiento básico de AlphaGenome, así como establecer la estructura inicial para la organización reproducible del proyecto.

## Datos de entrada

- GTF de Ensembl
  - Fuente: Ensembl
  - Archivo: `Homo_sapiens.GRCh38.115.abinitio.gtf.gz`

- GTF simulado
  - Fuente: archivo generado para pruebas
  - Archivo: `mini_clean2.gtf`

Los datos de entrada se encuentran en `raw-data/2026-04-06-exp01-installation-test/`.

## Scripts utilizados

- `scripts/2026-04-06-exp01-installation-test/process_gtf_test.py`
- `scripts/2026-04-06-exp01-installation-test/process_gtf.py`

## Resultados

- Ejecución correcta de las pruebas iniciales de procesamiento.
- Generación del archivo `mini_output.feather`.

## Trazabilidad

Los datos externos utilizados en el proyecto se registran en `docs/manifest-data.tsv`.

## Estado

Experimento completado.
