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

	El entorno incluye Alphagenome con sus dependencias necesarias.

- Disponer de una API key de AlphaGenome: 
	La API key no está incluida en el repositorio por seguridad.
	Cada usuario debe definirla como variable de entorno: export ALPHAGENOME_API_KEY="tu_api_key"
	Para hacerla persistente, se añade a al archivo de configuración con : vim ~/.bashrc
	Para hacer el cambio efectivo se ejecuta: source ~/.bashrc
	USO SIEMPRE DESDE LA RAÍZ DEL PROYECTO: python scripts/nombre_del_script.py

### Notas de uso de la API

No todos los scripts requieren una API key:

	- Scripts de procesamiento de datos (por ejemplo, manejo de archivos GTF) funcionan en local y no requieren API key.
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

### Día 4, 5 y 6
- Lectura y familiarización con la documentación de alphagenome.
- Elaboración y exploración del código para desarrollar el script para el análisis de parámetros relevantes en la evaluación de exon skipping.

### Día 7 y 8
- Pulimiento del script final exp02-score_exon_skipping.py
- Obtención de tablas de interés como primeros resultados útiles para la interpretación de exon skipping.

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

## Entorno

Entorno reproducible disponible en: environment/environment.yml  
