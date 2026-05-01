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
└── README.md                 # Documentación general del proyecto

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
- Desarrollo de un scripot para graficar conjuntamente los parámetros de splice_junctions, rna_seq, splice_sites y splice_site_usage por separado.

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
- scripts para replicar la figura 3b.

## Entorno

Entorno reproducible disponible en: environment/environment.yml  
