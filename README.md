# Modelización del Hábitat y Vulnerabilidad Climática del *Abies pinsapo* Boiss

Este repositorio contiene el código fuente y el flujo de trabajo estadístico completo desarrollado para mi **Trabajo de Fin de Grado (TFG) en Estadística** por la **Universidad de Sevilla**.

## 📋 Descripción del Proyecto
El objetivo principal de esta investigación es evaluar y comparar la robustez de una aproximación paramétrica tradicional (Modelo Lineal Generalizado - GLM) frente a dos algoritmos no paramétricos de *Machine Learning* (*Random Forest* y *Gradient Boosting Machine* - GBM) en la tarea de modelizar el nicho ecológico del pinsapo (*Abies pinsapo* Boiss.) y proyectar su vulnerabilidad geográfica ante escenarios de cambio climático (CMIP6: SSP2-4.5 y SSP5-8.5) con horizonte al año 2050.

## 🗂️ Estructura del Repositorio
* `/datos`: Carpeta destinada a almacenar la caché local de registros de GBIF y las capas raster de WorldClim.
* `script_sdm_pinsapo.R`: Script único en lenguaje R que automatiza todo el proceso de modelización y proyección.

## 🛠️ Instalación y Uso Rápido

Para evitar tiempos de espera lentos con las APIs cartográficas, sigue estos tres sencillos pasos:

1. **Clona o descarga** este repositorio en tu ordenador.
2. **Descarga el dataset pre-procesado** (`datos.zip`) desde este enlace de Google Drive:  
   👉 **[PEGA AQUÍ TU ENLACE COPIADO DE GOOGLE DRIVE]**
3. **Descomprime** el archivo `datos.zip` en la raíz del proyecto. Asegúrate de que se quede una carpeta llamada `/datos` con todos los archivos `.rds` y `.tif` dentro.

### Librerías de R necesarias
Abre tu entorno de R y asegúrate de tener instalados los siguientes paquetes antes de ejecutar el script:

```r
install.packages(c("rgbif", "dplyr", "geodata", "terra", "corrplot", 
                   "usdm", "caTools", "caret", "randomForest", "pdp", 
                   "ggplot2", "gridExtra", "gbm", "pROC", "blockCV", 
                   "sf", "rnaturalearth"))
