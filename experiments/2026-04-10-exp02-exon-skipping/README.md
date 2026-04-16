# Experiment exp02 — Exon Skipping (DLG1)

## Objetivo

Reproducir parcialmente un experimento descrito en el paper de AlphaGenome para una variante asociada a exon skipping en el gen DLG1.

El objetivo es obtener métricas (scores) que permitan observar cambios en splicing y apoyar la hipótesis de exon skipping.


## Datos de entrada

- Paper de referencia  
  - Fuente: Nature (AlphaGenome)  
  - Archivo: paper-alphagenome.pdf  

- Variante genética  
  - chromosome: chr3  
  - position: 197081044  
  - reference: TACTC  
  - alternate: T  

## Scripts utilizados

- scripts/2026-04-10-exp02-exon-skipping/exp02-score_exon_skipping.py  
- alphagenome_key.py  

## Ejecución

Desde la raíz del proyecto:

```bash
PYTHONPATH=. python scripts/2026-04-10-exp02-exon-skipping/exp02-score_exon_skipping.py

Nota: Se ejecuta con PYTHONPATH=. para que Python pueda encontrar alphagenome_key.py desde la raíz del proyecto.

## Resultados

- Creación correcta del cliente de AlphaGenome  
- Preparación del entorno para ejecutar score_variant()  
- Obtención de scores y análisis de resultados en formato .tsv. 

## Archivos asociados

- runlog.txt → registro reproducible de ejecución  
- manifest-data.tsv → inventario general de archivos del proyecto  

## Estado

Experimento en desarrollo.

Se ha completado:
- Definición de la variante
- Creación del cliente AlphaGenome
- Resolución de problemas de conectividad (DNS)
- Script para la obtención de scores para cada parámetro de interés en formato tabular. 

Pendiente:
- Tratar de replicar la fib 3b. para terminar añadir consistencia a la réplica del experimento y apoyar las interpretaciones.
