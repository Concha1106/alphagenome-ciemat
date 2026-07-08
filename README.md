#Repositorio para prácticas de bioinformática en CIEMAT usando AlphaGenome

Proyecto desarrollado durante prácticas en CIEMAT centrado en la exploración de AlphaGenome para el análisis de variantes genómicas y efectos en splicing.
Este proyecto integra procesamiento bioinformático local con modelos de deep learning accesibles vía API para la predicción de efectos de variantes.

## Objetivos

- Leer bibliografía acerca de AlphaGenome
- Instalar AlphaGenome
- Aprender uso de AlphaGenome
- Reproducir ejemplos del paper: Advancing regulatory variant effect prediction with AlphaGenome
- Explorar predicción de edición génica

## Requisitos

Para ejecutar los scripts es necesario:

- Tener Python y Conda instalados
- Crear el entorno a partir de:

	```bash
	conda env create -f environment/environment.yml
	conda activate alphagenome
	```
	El entorno incluye Alphagenome con sus dependencias necesarias.

- Disponer de una API key de AlphaGenome: 
	La API key no está incluida en el repositorio por seguridad.
	Cada usuario debe definirla como variable de entorno: export ALPHAGENOME_API_KEY="tu_api_key"
	Para hacerla persistente, se añade a al archivo de configuración con : vim ~/.bashrc
	Para hacer el cambio efectivo se ejecuta: source ~/.bashrc
	USO SIEMPRE DESDE LA RAÍZ DEL PROYECTO: python scripts/nombre_del_script.py

### Notas de uso de la API

No todos los scripts requieren una API key:

	- Scripts de procesamiento de datos (por ejemplo, manejo de archivos GTF) o de visualización de gráficos funcionan en local y no requieren API key.
	- Scripts que utilizan el modelo de AlphaGenome (por ejemplo, `score_variant`) sí requieren una API key configurada.

Esto se debe a que el modelo se ejecuta a través de la API, mientras que otras utilidades son locales.

## Arquitectura del proyecto: uso local vs remoto

Este proyecto combina dos tipos de funcionalidades dentro de AlphaGenome:

### Procesamiento local (sin API)

Incluye scripts que trabajan directamente con archivos en local, como:

	- Procesamiento de archivos GTF
	- Manipulación de datos de entrada

Ejemplos:
- `process_gtf.py`
- `process_gtf_test.py`

Estas operaciones no requieren conexión a la API, así no  necesitan API key. Se ejecutan completamente en el entorno local.

### Predicción mediante modelo (requiere API)

Incluye scripts que utilizan el modelo de AlphaGenome para predecir efectos de variantes genómicas.

Ejemplo de uso:
```python
dna_model = dna_client.create(API_KEY)
scores = dna_model.score_variant(...)

Estas operaciones requieren conexión a la API de AlphaGenome, necesitan una API key válida y ejecutan inferencia del modelo de forma remota. 

## Estructura del proyecto

alphagenome-ciemat/
├── alphagenome_key.py        # Cliente reutilizable para AlphaGenome
├── .gitignore                # Configuración de exclusión de archivos en Git
├── scripts/                  # Scripts de análisis
├── experiments/              # Experimentos organizados
├── raw-data/                 # Datos originales
├── interm-data/              # Datos intermedios
├── results/                  # Resultados finales
├── docs/                     # Documentación (manifest, paper, etc.)
├── config/                   # Configuración
├── environment/              # Entorno reproducible (environment.yml)
└── README.md                 # Documentación general del proyeto
-- common.py                  # Funciones reutilizables compartidas
-- config.py                  # Constantes compartidas para visualización

## Progreso

### Día 1
- Instalación de AlphaGenome y sus dependencias
- Primer procesamiento GTF con script process_gtf_test.py  (uso un archivo de prueba procedente de Emsembl)
- Prueba procesamiento con script process_gtf.py para formato GTF (archivo ficticio elaborado para pruebas)
- Registro de las tareas principales en git
### Día 2
- Enviroment set up: exporto el entorno Conda tras instalar y probar que funciona exitosamente la herramienta. Para ello:
	conda activate alphagenome
	conda env export > environment.yml
### Día 3
- Reestructuración del proyecto para reproducibilidad
   - Creación de carpetas `config` y `docs`
   - Eliminación de carpeta `notes`
   - Renombrado del experimento a formato reproducible
   - Reorganización de datos de entrada a `raw-data`
- Creación de `runlog` del experimento exp01
- Creación de `manifest` del experimento exp01
- Creación de `.gitignore`
- Mejora de documentación del proyecto

### Días 4, 5 y 6
- Lectura y familiarización con la documentación de alphagenome.
- Elaboración y exploración del código para desarrollar el script para el análisis de parámetros relevantes en la evaluación de exon skipping.

### Días 7, 8 y 9
- Pulimiento del script final exp02-score_exon_skipping.py
- Obtención de tablas de interés como primeros resultados útiles para la interpretación de exon skipping.

### Días 10 y 11 
- Reflexión e interpretación de los parámetros obtenidos a partir del script final exp02-score_exon_skipping.py
- Pruebas orientadas a comprender el análisis del parámetro combined_score en combinación con el splice_junctions, splice_sites y splice_sites_usage.
- Exploración de la posible junction causante de exon skipping.

### Día 12
- Búsqueda de los exones del gen DLG1 en el genoma de referencia a través de Ensembl.
- Confirmación de la coincidencia de las coordenada de inicio y final de la junction candidata con exones no consecutivos, dejando un exón intermedio fuera, dejando entrever la presencia de exon skipping.

### Día 13 
- Recuepración y reestructuración del git del proyecto a causa de un mensaje que mostraba que estaba corrupto.

### Días 14, 15 y 16
- Desarrollo de un script para ejecutar la función predict_variant() con el objetivo de obtener los valores de REF y ALT por separado y poder evaluar pérdidas o ganancias de junctions, así como los scores de parametros asociados que puedan servir de apoyo para confirmar la hipótesis y lograr representar la figura 3b del paper.

### Días 17, 18 y 19
- Desarrollo de tres scripts individuales que permiten graficar los parámetros splice_junctions, rna_seq y splice_sites/splice_site_usage por separado. 

### Día 20
- Desarrollo de scripts específicos para representar los principales outputs de `predict_variant()` (RNA-seq, splice junctions, splice sites y splice site usage) e integración preliminar en una figura conjunta inspirada en la Figura 3b del artículo de AlphaGenome.

### Dias 21-35

- Estudio de la documentación oficial de AlphaGenome.
- Profundización teórica y exploratoria de score_variant() y predict_variant(), así como familiarización con los métodos predict_sequence(), predict_interval() y scovre_variants(). 
- Exploración de los diferentes output types y variant scorers disponibles.
- Comprensión de señales biológicas relacionadas con splicing, expresión y regulación génica, estructura cromatínica y estructura 3D.
- Análisis detallado de la Figura 3b del artículo de AlphaGenome.

### Días 36 y 37
- Reflexión sobre posibles estrategias para automatizar la evaluación de variantes novedosas.
- Diseño conceptual de un pipeline general para análisis exploratorio de variantes.

### Día 38

- Inicio del experimento exp03.
- Creación de la estructura reproducible del experimento.
- Redacción del README específico del experimento.
- Implementación inicial del pipeline mediante argparse.
- Implementación de generación automática de runlogs.

### Días 39, 40 y 41

- Desarrollo de la fase de validación y representación interna de variantes para exp03.
- Construcción del objeto `Variant` de AlphaGenome a partir de cromosoma, posición, REF y ALT.
- Construcción del intervalo de entrada mediante `variant.reference_interval.resize()` usando longitudes compatibles con AlphaGenome.
- Configuración del cliente de AlphaGenome mediante `get_dna_model()`.
- Ejecución de `score_variant()` y exportación de `score_variant_all.tsv`.
- Ejecución de `predict_variant()` restringido a una ontología concreta mediante `--ontology-curie`, para evitar solicitar todas las ontologías disponibles debido a posible saturación del sistema.
- Exportación de predicciones REF, ALT y delta para:
  - TrackData;
  - splice junctions;
  - contact maps cuando existan tracks disponibles.
- Generación de tablas top positivas y negativas por `output_type` para facilitar la interpretación en Excel.

### Día 42

- Pruebas exploratorias del pipeline exp03 con una variante compleja en SEC23B.
- Selección inicial de ontologías candidatas (CL:0000558, UBERON:0002371 y CL:0001059) a partir de metadatos de AlphaGenome y relevancia biológica eritroide/hematopoyética.
- Identificación de cuello de botella de memoria en la exportación local de `predict_variant()`.
- Optimización inicial de exportación TrackData mediante escritura progresiva por track.
- Mejora del manejo de `--output-types` para permitir ejecuciones parciales por modalidad.

### Días 43, 44 y 45

- Mejora del pipeline exp03 para facilitar la interpretación de variantes complejas.
- Incorporación de métricas de longitud de variante y longitud de junctions.
- Cálculo de una métrica integrada de splicing basada en `score_variant()`.
- Exportación de `splice_sites` en formato wide, reduciendo la complejidad de visualización en Excel.
- Generación de vistas locales alrededor de la variante.
- Implementación de `predict_variant_track_summary.tsv` para resumir el impacto por track y orientar la interpretación de outputs con múltiples señales.

### Días 46-53
- Establecimiento de la estrategia de análisis e interpretación del KI en SEC23B para evaluar su impacto funcional.
- Exploración y elección de ontologías relevantes para este caso de estudio.
- Ejecución del pipeline del exp03 end diferentes ontologías.
- Interpretación por separado y conjunta de los distintos output types obtenidos en cada ontología para la evaluación del cambio en la expresión y del patrón de splicing debido al  KI simulado en AlphaGenome.

### Días 54, 55 y 56

- Serialización del objeto completo devuelto por predict_variant() mediante pickle, generando prediction.pkl junto a los TSV.
- Separación conceptual entre:
  --> pipeline de análisis/exportación (run-variant-pipeline.py)
  --> script de visualización (visualize-prediction.py)
- Desarrollo inicial de un script de visualización basado en objetos de AlphaGenome previamente guardados.
- Implementación de visualización local de RNA-seq REF vs ALT y delta ALT-REF.
- Incorporación de anotación MANE Select de SEC23B para contextualizar los cambios predichos en exones e intrones.
- Parametrización de la región visualizada mediante --region-start, --region-end y/o ventana local inferida desde runlog.txt.
- Inclusión de marca visual del punto de inserción KI y coordenada genómica real.

### Días 57 y 58

- Ampliación del script visualize-prediction.py para generar visualizaciones locales de splice_sites, splice_site_usage y splice_junctions a partir de prediction.pkl.
- Implementación de selección parametrizable de splice junctions según delta ALT-REF positivo y negativo.
- Desarrollo de una representación tipo arco para splice_junctions, codificando coordenadas genómicas, longitud de junction, magnitud del delta y dirección del cambio.
- Incorporación de visualization_runlog.txt para registrar parámetros de visualización, región usada, umbrales de junctions seleccionadas y figuras generadas.
- Queda pendiente unificar estilo visual entre figuras y refactorizar el script.

## Días 59-62
- Refactorización del script de visualización para reducir código duplicado y mejorar su mantenibilidad.
- Homogeneización del diseño de las figuras generadas y reorganización interna del código.


## Experimentos

### exp01 — Installation Test

Objetivo:
Verificar instalación y funcionamiento básico de AlphaGenome

Ubicación: experiments/2026-04-06-exp01-instalation-test/
Incluye:
- README del experimento
- runlog reproducible

### exp02 — Exon Skipping (DLG1)

Objetivo:
Reproducir parcialmente un experimento del paper de AlphaGenome para una variante asociada a exon skipping en el gen DLG1

Ubicación: experiments/2026-04-10-exp02-exon-skipping/
Incluye:
- README del experimento
- runlog reproducible
- uso de score_variant para análisis de splicing en script exp02-score_exon_skipping.py
- uso de predict_variant para completar el análisis de splicing en script exp02-ref_vs_alt_prediction.py
- scripts independientes para reproducir de forma aproximada la Figura 3b del artículo mediante representaciones independientes de RNA-seq, splice junctions, splice sites y splice site usage.

### exp03 - General Variant Evaluation Pipeline

Objetivo:
Desarrollar un pipeline general y automatizable para la evaluación exploratoria de variantes genómicas mediante AlphaGenome.

A diferencia de exp02, centrado en un caso concreto de exon skipping, este experimento busca construir una arquitectura reutilizable capaz de:

- recibir variantes genómicas como entrada;
- ejecutar análisis mediante score_variant();
- ejecutar análisis mediante predict_variant();
- exportar resultados estructurados para su interpretación;
- mantener trazabilidad mediante runlogs y control de versiones.

Ubicación:
experiments/2026-05-28-exp03-variant-pipeline/

Estado actual:
- estructura inicial del pipeline implementada;
- lectura de argumentos mediante argparse;
- validación básica de variante;
- construcción de intervalo mediante AlphaGenome;
- conexión con el modelo usando API key desde variable de entorno;
- ejecución de score_variant();
- ejecución de predict_variant() para una ontología concreta;
- exportación de tablas completas REF/ALT/delta;
- exportación de tablas top positivas y negativas por output type;
- generación automática de runlogs.
- generación de figuras locales para RNA-seq, splice_sites, splice_site_usage y splice_junctions;
- generación de visualization_runlog.txt para trazabilidad de la visualización;

Últimas mejoras

- separación completa entre pipeline de análisis y pipeline de visualización;
- visualización reutilizable a partir de prediction.pkl;
- carga automática de anotación génica desde GTF;
- selección automática del transcript MANE Select (o transcript indicado por el usuario);
- compatibilidad con cualquier gen anotado en GRCh38;
- validación previa de argumentos y mensajes de error orientados al usuario;
- generación automática de visualization_runlog.txt;
- exportación automática de la tabla de splice junctions representadas.
- refactorización del visualizador mediante funciones auxiliares compartidas para reducir duplicación de código;
- centralización de constantes de representación gráfica en `config.py`;
- reorganización del script de visualización por bloques funcionales para mejorar su legibilidad y mantenimiento;
- homogeneización del diseño de las figuras generadas.

## Entorno

Entorno reproducible disponible en: environment/environment.yml  
