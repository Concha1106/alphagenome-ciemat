#Repositorio para prácticas de bioinformática en CIEMAT usando AlphaGenome

## Objetivos

- Leer bibliografía acerca de AlphaGenome
- Instalar AlphaGenome
- Aprender uso de AlphaGenome
- Reproducir ejemplos del paper: Advancing regulatory variant effect prediction with AlphaGenome
- Explorar predicción de edición génica

## Estructura del proyecto


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


## Experimentos

### exp01 — Installation Test

Objetivo:
Verificar instalación y funcionamiento básico de AlphaGenome

Ubicación: experiments/2026-04-06-exp01-instalation-test/
Incluye:
- README del experimento
- runlog reproducible
- manifest de archivos


## Entorno

Entorno reproducible disponible en: environment/environment.yml  
