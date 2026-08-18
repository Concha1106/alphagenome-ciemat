# Evaluación funcional de estrategias de edición génica mediante AlphaGenome

Proyecto desarrollado en el CIEMAT para evaluar AlphaGenome como herramienta de apoyo a la interpretación funcional de estrategias de edición génica. El trabajo incluye el desarrollo de un pipeline bioinformático reproducible, su validación mediante el caso de exon skipping en DLG1 y su aplicación al análisis de un knock-in terapéutico en SEC23B.

## Objetivos

- Explorar las capacidades de AlphaGenome para la predicción de efectos funcionales de variantes genómicas.
- Desarrollar un pipeline bioinformático reproducible para la ejecución, procesamiento, exportación y visualización de predicciones de AlphaGenome.
- Validar el flujo de trabajo mediante la reproducción del caso de exon skipping en DLG1 descrito por Avsec et al. (2026).
- Aplicar el pipeline al análisis funcional de un knock-in terapéutico en SEC23B.

## Requisitos

Para ejecutar los scripts es necesario:

- Tener Python y Conda instalados
- Crear el entorno a partir de:

	```bash
	conda env create -f environment/environment.yml
	conda activate alphagenome
	```
	El entorno incluye Alphagenome con sus dependencias necesarias.

- Disponer de una API key de AlphaGenome. La API key no está incluida en el repositorio por seguridad y debe definirse como variable de entorno:

    ```bash
    export ALPHAGENOME_API_KEY="tu_api_key"
    ```

    Para hacerla persistente, puede añadirse esta línea al archivo `~/.bashrc` y ejecutar:

    ```bash
    source ~/.bashrc
    ```

Los scripts deben ejecutarse desde la raíz del proyecto.

### Notas de uso de la API

No todos los scripts requieren una API key:

	- Scripts de procesamiento de datos (por ejemplo, manejo de archivos GTF) o de visualización de gráficos funcionan en local y no requieren API key.
	- Scripts que utilizan el modelo de AlphaGenome (por ejemplo, `score_variant`) sí requieren una API key configurada.

Esto se debe a que el modelo se ejecuta a través de la API, mientras que otras utilidades son locales.

## Estructura del proyecto

alphagenome-ciemat/
├── alphagenome_key.py        # Cliente reutilizable para AlphaGenome
├── common.py                 # Funciones reutilizables compartidas
├── config.py                 # Constantes compartidas para visualización
├── .gitignore                # Configuración de exclusión de archivos en Git
├── scripts/                  # Scripts organizados por experimento
├── experiments/              # Documentación específica de los experimentos
├── raw-data/                 # Datos originales
├── interm-data/              # Datos intermedios
├── results/                  # Resultados generados
├── docs/                     # Documentación y manifest de datos
├── environment/              # Entorno reproducible (environment.yml)
└── README.md                 # Documentación general del proyecto


## Experimentos

### exp01 — Installation Test

**Objetivo:** verificar la instalación y el funcionamiento básico de AlphaGenome y establecer la estructura reproducible inicial del proyecto.

**Documentación:** `experiments/2026-04-06-exp01-installation-test/`

**Scripts:** `scripts/2026-04-06-exp01-installation-test/`

**Resultados:** `results/2026-04-06-exp01-installation-test/`

### exp02 — Exon Skipping (DLG1)

**Objetivo:** explorar las capacidades de AlphaGenome para el análisis de variantes y validar el flujo de trabajo mediante la reproducción del caso de exon skipping en DLG1 descrito por Avsec et al. (2026).

**Documentación:** `experiments/2026-04-10-exp02-exon-skipping/`

**Scripts:** `scripts/2026-04-10-exp02-exon-skipping/`

**Resultados:** `results/2026-04-10-exp02-exon-skipping/`

El experimento incluye el análisis mediante `score_variant()` y `predict_variant()`, así como la representación de los principales outputs relacionados con el caso de exon skipping: RNA-seq, splice junctions, splice sites y splice site usage.

### exp03 - General Variant Evaluation Pipeline

**Objetivo:** desarrollar un pipeline generalizable y reproducible para la evaluación funcional de variantes genómicas mediante AlphaGenome y aplicarlo al análisis del knock-in terapéutico en SEC23B.

**Documentación:** `experiments/2026-05-28-exp03-variant-pipeline/`

**Scripts:** `scripts/2026-05-28-exp03-variant-pipeline/`

**Resultados:** `results/2026-05-28-exp03-variant-pipeline/`

El pipeline se organiza en dos módulos principales:

- `run-variant-pipeline.py`: ejecuta `score_variant()` y `predict_variant()`, procesa y exporta los resultados y serializa las predicciones completas en `prediction.pkl`.
- `visualize-prediction.py`: reutiliza `prediction.pkl` para generar visualizaciones de RNA-seq, splice sites, splice site usage y splice junctions sin realizar nuevas consultas a la API.

El flujo incorpora selección de ontología, anotación génica automática a partir de GTF, selección del transcrito MANE Select, exportación estructurada de resultados y generación de archivos de registro para favorecer la trazabilidad y reproducibilidad del análisis.


## Entorno

Entorno reproducible disponible en: environment/environment.yml  
