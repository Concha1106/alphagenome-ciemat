Experiment exp02 — Exon Skipping (DLG1)

## Objetivo

Reproducir un experimento descrito en el paper de AlphaGenome para una variante asociada a exon skipping en el gen DLG1.


El objetivo es obtener métricas resumidas (score_variant()) y predicciones directas REF vs ALT (predict_variant()) que permitan observar cambios en splicing y apoyar la hipótesis de exon skipping descrita en el paper.

El trabajo se desarrolla en dos fases principales:
	1. Análisis inicial mediante score_variant(), orientado a detectar señales globales de alteración del splicing.
	2. Comparación directa REF vs ALT mediante predict_variant(), orientada a validar estructuralmente el evento y aproximarse a una representación similar a la Figura 3b del paper.

## 1. ANÁLISIS INICIAL CON score_variant()

### Datos de entrada

- Paper de referencia  
  - Fuente: Nature (AlphaGenome)  
  - Archivo:docs/alphagenome.pdf  

- Variante genética  
  - chromosome: chr3  
  - position: 197081044  
  - reference: TACTC  
  - alternate: T  

### Scripts utilizados

- alphagenome_key.py 
- scripts/2026-04-10-exp02-exon-skipping/exp02-score_exon_skipping.py: este script calcula scores resumidos del efecto de lavariante sobre splicing y genera tablas resumen para los parámetros splice site usage, splice sites, splice junctions; y un summary final con una junction candidata compatible con exon skipping.  

### Ejecución

Desde la raíz del proyecto:

```bash
PYTHONPATH=. python scripts/2026-04-10-exp02-exon-skipping/exp02-score_exon_skipping.py
```
Nota: Se ejecuta con PYTHONPATH=. para que Python pueda encontrar alphagenome_key.py desde la raíz del proyecto.

### Resultados

- Creación correcta del cliente de AlphaGenome  
- Preparación del entorno para ejecutar score_variant()  
- Obtención de scores y análisis de resultados en formato .tsv. 

### Archivos de salida principales:
    - results/2026-04-10-exp02-exon-skipping/dlg1_usage_tibial.tsv
    - results/2026-04-10-exp02-exon-skipping/dlg1_sites.tsv
    - results/2026-04-10-exp02-exon-skipping/dlg1_junctions_tibial.tsv
    - results/2026-04-10-exp02-exon-skipping/dlg1_summary.tsv

### Interpretación preliminar
Se calculó la métrica heuristica conjunta de splicing recomendada por AlphaGenome: 'max(splice_sites) + max(splice_site_usage) + max (splice_junction)/5'

El valor obtenido en la variante del paper fue de "3.61", lo cual sugiere un efecto fuerte de la variante sobre el splicing.
Este resultado, junto con la señal detectada en "Artery_Tibial" y la presencia de una splice junction de mayor longitud, compatible con skipping, apoya la hipótesis de un evento de exon skipping en DLG1.

Además, en la tabla de junctions se calcula la longitud de cada splice junction y una métrica adicional (max_tibial_score) que resume la mayor señal absoluta observada entre Artery_Tibial y Nerve_Tibial.

Esta interpretación es preliminar y se ajustará posteriormente mediante el análisis de predicciones directas de REF/ALT a través del desarrollo de un script con la función predict_variant(), con la que se tratará de recrear la figura 3b del paper. 

## 1.1 ANÁLISIS EXPLORATORIOS DE VARIANTES 

- Sensibilidad posicional- se introdujeron variantes en distintas posiciones respecto al locus original:
	- Variantes situadas a -42 nt y +10 nt mostraron valores bajos de combined_score, usage y splice sites, sin evidencia de un efecto relevante.
	- En cambio, una variante situada a +2 nt presentó valores elevados en todos los parámetros, incluyendo junctions con scores altos.

Esto sugiere que la respuesta del modelo depende fuertemente del contexto local, siendo la región inmediata a la variante original especialmente sensible a cambios.

- Efecto del tipo de mutación- se evaluaron distintas modificaciones en la misma posición:
	- Sustitución completa (TACTC>ATGAG)
	- Deleción parcial (TACTC>TA)

Ambas produjeron scores elevados y aparición de junctions similares, lo que indica que el efecto no depende de una única mutación concreta, sino de la alteración del motivo de secuencia local.

- Comparación referencia vs mutación
La comparación entre la secuencia de referencia y las variantes muestra que:

La referencia no presenta junctions con señal relevante en comparación con las variantes evaluadas.
Las variantes inducen múltiples junctions con scores altos, destacando de forma recurrente:
197076685 - 197085579 (~8894 bp)

Esto sugiere una reorganización consistente del patrón de splicing predicho en presencia de la mutación.

### Interpretación conjunta
Estos resultados indican que el modelo responde de forma coherente al contexto de la variante, siendo la región analizada sensible a perturbaciones de secuencia.
Diferentes alteraciones en esa zona pueden inducir cambios relevantes en el splicing predicho.

Aunque no se demuestra de forma directa el mecanismo exacto, los resultados son compatibles con una alteración significativa del splicing, potencialmente asociada a eventos de exon skipping.

## 1.2 VALIDACIÓN ESTRUCTURAL MEDIANTE ANOTACIÓN (GTF)

Para apoyar la interpretación de las splice junctions predichas, se utilizó la anotación genómica de Ensembl (GRCh38, release 115).

Se filtraron los exones del gen DLG1 a partir del archivo GTF completo, generando un subconjunto intermedio ('dlg1_exons.gtf') que permite contrastar las coordenadas de las junctions predichas con la estructura exónica real.

Este análisis fue consistente con la hipotesis de que la junction candidata conecta exones no consecutivos, lo cual apoya la idea de la existencia de un evento de exon skipping.


## 2. COMPARACIÓN DIRECTA REF vs ALT CON 'predict_variant()'

### Datos de entrada
 
Se utilizan los mismos datos de entrada descritos en la sección anterior:
- variante chr3:197081044 TACTC>T
- paper de referencia de AlphaGenome

### Scripts utilizados:
- alphagenome_key.py
- scripts/2026-04-10-exp02-exon-skipping/exp02-ref_vs_alt_prediction.py: este segundo script genera predicciones directas de referencia (REF) frente a variante (ALT) utilizando la función 'predict_variant()' de AlphaGenome.

A diferencia del análisis anterior con 'score_variant()', que resume el efecto global de la variante mediante métricas agregadas, este enfoque permite inspeccionar directamente cómo cambian las señales de splicing y cobertura a nivel posicional, facilitando una futura representación visual similar a la Figura 3b del paper.

Se solicitaron los siguientes outputs:

- `SPLICE_JUNCTIONS`
- `SPLICE_SITE_USAGE`
- `SPLICE_SITES`
- `RNA_SEQ`

El tejido evaluado se restringió a 'UBERON:0007610', correspondiente a GTEx 'Artery_Tibial', ya que el paper describe el evento de exon skipping en tejido arterial tibial.

Dado que DLG1 se encuentra anotado en hebra negativa ('-'), las salidas relacionadas con splicing ('splice junctions', 'splice sites' y 'splice site usage') se filtraron específicamente a dicha hebra para evitar interpretar señales de la orientación opuesta.

Además, para facilitar la futura visualización, los datos tipo 'TrackData' se recortaron a una región local alrededor del evento de interés:

chr3:197075000–197087000

Esta región incluye los exones implicados en la hipótesis de exon skipping y permite comparar de forma más clara las diferencias entre REF y ALT.

### Ejecución

Desde la raíz del proyecto:

```bash
PYTHONPATH=. python scripts/2026-04-10-exp02-exon-skipping/exp02-ref_vs_alt_prediction.py 
```
Nota: Se ejecuta con PYTHONPATH=. para que Python pueda encontrar alphagenome_key.py desde la raíz del proyecto.

### Identificación estructural de exones implicados

Para interpretar correctamente las splice junctions predichas por predict_variant(), fue necesario localizar los exones reales implicados en la región de estudio utilizando la anotación oficial de Ensembl (GRCh38, release 115).

Aunque el paper proporciona directamente la variante chr3:197081044 TACTC>T fue necesario identificar qué transcrito concreto de DLG1 debía utilizarse como referencia estructural.

Para la selección del transcrito de referencia se filtraron los transcritos anotados para DLG1 en el archivo GTF y se seleccionó el transcrito marcado como MANE_Select:

grep 'gene_name "DLG1"' Homo_sapiens.GRCh38.115.gtf | \
grep $'\ttranscript\t' | \
grep 'MANE_Select'

El transcrito obtenido fue: ENST00000667157

Este transcrito corresponde a la isoforma canónica más representativa y consensuada entre Ensembl y RefSeq.

El uso de MANE_Select permite evitar ambigüedades entre múltiples isoformas y facilita una interpretación más robusta de los eventos de splicing.

Extracción de exones

A partir de este transcrito, se extrajeron sus exones anotados para visualizar con claridad la estructura local de la región afectada:
```
echo -e "chr\tstart\tend\tstrand\texon_number" > dlg1_mane_select_exons.tsv

grep 'transcript_id "ENST00000667157"' Homo_sapiens.GRCh38.115.gtf | \
grep $'\t'exon$'\t' | \
awk '{
    for(i=1; i<=NF; i++) {
        if($i == "exon_number") {
            exon = $(i+1)
            gsub(/"/, "", exon)
        }
    }
    print $1 "\t" $4 "\t" $5 "\t" $7 "\t" exon
}' >> dlg1_mane_select_exons.tsv
``
Esto permitió identificar los exones relevantes alrededor de la variante:

exon 16: 197085580–197085756
exon 17: 197081051–197081117
exon 18: 197076586–197076685

La relación con las splice junctions predichas por AlphaGenome:

197076685 - 197081050
197081117 - 197085579
197076685 - 197085579

coinciden con los límites exón–intrón anotados en el GTF.

La aparición de la junction:

197076685 - 197085579

indica una conexión directa entre los exones 18 y 16, omitiendo el exón intermedio (exón 17). Esto es consistente con un evento de exon skipping, exactamente el fenómeno descrito en la Figura 3b del paper.

### Resultados:
El script genera tablas en formato '.tsv' listas para graficar:

- 'dlg1_splice_junctions_ref_vs_alt.tsv'
- 'dlg1_splice_site_usage_ref_vs_alt.tsv'
- 'dlg1_splice_sites_ref_vs_alt.tsv'
- 'dlg1_rna_seq_ref_vs_alt.tsv'

Estas salidas permitirán analizar con mayor precisión qué junctions se pierden, cuáles aparecen de novo y cómo cambia la señal de uso de sitios de splicing, aportando una validación más directa de la hipótesis de exon skipping descrita en el artículo.

- Refinamiento de interpretación mediante tablas Top10: Para facilitar la interpretación biológica de los resultados obtenidos con predict_variant(), se añadieron tablas resumen que priorizan los eventoscon mayor cambio entre ALT y REF (delta_alt_ref)

Aunque las tablas completas permiten un análisis exhaustivo, su tamaño dificulta la identificación rápida de los eventos más relevantes. Por ello, se generaron tablas adicionales con los 10 eventos que más aumentan y los 10 que más disminuyen tras aplicar la variante.

Esto permite detectar de forma más directa:

- splice junctions ganadas o perdidas
- splice sites debilitados o reforzados
- cambios relevantes en splice site usage

Especialmente en el caso de splice_junctions, este enfoque ayuda a localizar rápidamente la aparición de junctions largas compatibles con exon skipping.
Se generaron los siguientes archivos: 

- dlg1_splice_junctions_top10_increased_alt_vs_ref.tsv
- dlg1_splice_junctions_top10_decreased_alt_vs_ref.tsv
- dlg1_splice_site_usage_top10_increased_alt_vs_ref.tsv
- dlg1_splice_site_usage_top10_decreased_alt_vs_ref.tsv
- dlg1_splice_sites_top10_increased_alt_vs_ref.tsv
- dlg1_splice_sites_top10_decreased_alt_vs_ref.tsv

Estas tablas resumen permiten priorizar las señales más relevantes antes de construir la visualización final tipo Figura 3b.

En particular, la aparición entre los valores más aumentados de una junction larga que conecta directamente los exones flanqueantes 197076685-197085579 refuerza la hipótesis de exon skipping del exón 17, ya que implica la omisión del exón intermedio y reproduce estructuralmente el fenómeno descrito en el paper

Este refinamiento mejora la trazabilidad biológica del análisis y facilita la justificación de la replicación incluso antes de generar la figura final.

## 3. SCRIPTS PARA LA VISUALIZACIÓN DE LA FIGURA 3b.

En esta fase del experimento se han desarrollado varios scripts para la visualización de los efectos de la variante sobre distintos parámetros de splicing:

- "plot-fig3b-junctions.py"  
  Representación de las splice junctions predichas por AlphaGenome.

- "plot-fig3b-rna-seq.py"  
  Visualización de la señal de RNA-seq asociada a la variante.

- "plot-fig3b-splice-sites-and-usage.py"  
  Representación de los splice sites y su uso (splice site usage).

- "plot-fig3b-complete.py"
  Integración de todos los parámetros anteriores en una única figura, aproximando la Figura 3b del artículo.

### ACTUALIZACIONES

- Se ha modificado el script `plot-fig3b-junctions.py` para mejorar la representación gráfica de las junctions y se ha añadido una breve descripción.
- Se han añadido nuevos scripts para representar de forma separada los distintos parámetros del modelo.
- Se ha comenzado la integración de todos los outputs en una visualización conjunta.


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
- Análisis exploratorio de sensibilidad posicional y tipo de mutación
- Validación estructural mediante anotación GTF
- Desarrollo del script predict_variant() para comparación REF vs AL y se ha añadido una breve descripción.
- Generación de tablas listas para visualización final
- Generación de visualización final tipo Figura 3b

Pendiente:
- Interpretación comparativa final de junctions perdidas/ganadas entre REF y ALT
