Experiment exp02 — Exon Skipping (DLG1)

## Objetivo

Reproducir un experimento descrito en el paper de AlphaGenome para una variante asociada a exon skipping en el gen DLG1.

El objetivo es obtener métricas (scores) que permitan observar cambios en splicing y apoyar la hipótesis de exon skipping.
Más adelante será necesario comparar predicciones referencia vs variante para interpretar con mayor precisión el mecanismo subyacente.


## Datos de entrada

- Paper de referencia  
  - Fuente: Nature (AlphaGenome)  
  - Archivo:docs/alphagenome.pdf  

- Variante genética  
  - chromosome: chr3  
  - position: 197081044  
  - reference: TACTC  
  - alternate: T  

## Scripts utilizados

- alphagenome_key.py 
- scripts/2026-04-10-exp02-exon-skipping/exp02-score_exon_skipping.py: calcula scores resumidos de efecto de variante sobre splicing y genera tablas de splice site usage, splice sites, splice junctions y un summary final con una junction candidata compatible con exon skipping.  


## Ejecución

Desde la raíz del proyecto:

```bash
PYTHONPATH=. python scripts/2026-04-10-exp02-exon-skipping/exp02-score_exon_skipping.py

´´´
Nota: Se ejecuta con PYTHONPATH=. para que Python pueda encontrar alphagenome_key.py desde la raíz del proyecto.

## Resultados

- Creación correcta del cliente de AlphaGenome  
- Preparación del entorno para ejecutar score_variant()  
- Obtención de scores y análisis de resultados en formato .tsv. 

### Archivos de salida principales:
    - results/2026-04-10-exp02-exon-skipping/dlg1_usage_tibial.tsv
    - results/2026-04-10-exp02-exon-skipping/dlg1_sites.tsv
    - results/2026-04-10-exp02-exon-skipping/dlg1_junctions_tibial.tsv
    - results/2026-04-10-exp02-exon-skipping/dlg1_summary.tsv

## Interpretación preliminar
Se calculó la métrica heuristica conjunta de splicing recomendada por AlphaGenome: 'max(splice_sites) + max(splice_site_usage) + max (splice_junction)/5'

El valor obtenido en la variante del paper fue de "3.61", lo cual sugiere un efecto fuerte de la varainte sobre el splicing.
Este resultado, junto con la señal detectada en "Artery_Tibial" y la presencia de una splice junction de mayor longitud, compatible con skipping, apoya la hipótesis de un evento de "exon skipping en DLG1".

Además, en la tabla de junctions se calcula la longitud de cada splice junction y una métrica adicional (max_tibial_score) que resume la mayor señal absoluta observada entre Artery_Tibial y Nerve_Tibial.

Esta interpretación es preliminar y se ajustará posteriormente mediante el análisis de predicciones directas de REF/ALT a través del desarrollo de un script con la función predict_variant(), con la que se tratará de recrear la figura 3b del paper. 

## Análisis exploratorios de variantes

- Sensibilidad posicional
Se introdujeron variantes en distintas posiciones respecto al locus original:
	- Variantes situadas a -42 nt y +10 nt mostraron valores bajos de combined_score, usage y splice sites, sin evidencia de un efecto relevante.
	- En cambio, una variante situada a +2 nt presentó valores elevados en todos los parámetros, incluyendo junctions con scores altos.

Esto sugiere que la respuesta del modelo depende fuertemente del contexto local, siendo la región inmediata a la variante original especialmente sensible a cambios.

- Efecto del tipo de mutación
Se evaluaron distintas modificaciones en la misma posición:
	- Sustitución completa (TACTC>ATGAG)
	- Deleción parcial (TACTC>TA)

Ambas produjeron scores elevados y aparición de junctions similares, lo que indica que el efecto no depende de una única mutación concreta, sino de la alteración del motivo de secuencia local.

- Comparación referencia vs mutación
La comparación entre la secuencia de referencia y las variantes muestra que:

La referencia no presenta junctions con señal relevante en comparación con las variantes evaluadas.
Las variantes inducen múltiples junctions con scores altos, destacando de forma recurrente:
197076685 → 197085579 (~8894 bp)

Esto sugiere una reorganización consistente del patrón de splicing predicho en presencia de la mutación.

### Interpretación conjunta
Estos resultados indican que el modelo responde de forma coherente al contexto de la variante, siendo la región analizada sensible a perturbaciones de secuencia.
Diferentes alteraciones en esa zona pueden inducir cambios relevantes en el splicing predicho.

Aunque no se demuestra de forma directa el mecanismo exacto, los resultados son compatibles con una alteración significativa del splicing, potencialmente asociada a eventos de exon skipping.

## Validación estructural mediante anotación (GTF)

Para apoyar la interpretación de las splice junctions predichas, se utilizó la anotación genómica de Ensembl (GRCh38, release 115).

Se filtraron los exones del gen DLG1 a partir del archivo GTF completo, generando un subconjunto intermedio (`dlg1_exons.gtf`) que permite contrastar las coordenadas de las junctions predichas con la estructura exónica real.

Este análisis fue consistente con la hipotesis de que la junction candidata conecta exones no consecutivos, lo cual apoya la idea de la existencia de un evento de exon skipping.

## Archivos asociados

- runlog.txt: registro reproducible de ejecución  
- manifest-data.tsv: inventario general de archivos del proyecto  

## Estado

Experimento en desarrollo.

Se ha completado:
- Definición de la variante
- Creación del cliente AlphaGenome
- Resolución de problemas de conectividad (DNS)
- Script para la obtención de scores para cada parámetro de interés en formato tabular. 
- Cálculo de longitud de splice junctions y resumen de señal tibial por junction

Pendiente:
- Desarrollo de predicciones directas de REF/ALT con función predict_variant().
- Tratar de replicar la fib 3b. para terminar añadir consistencia a la réplica del experimento y apoyar las interpretaciones.
