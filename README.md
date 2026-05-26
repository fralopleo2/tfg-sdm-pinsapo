# Modelización del Hábitat y Vulnerabilidad Climática del *Abies pinsapo* Boiss

Este repositorio contiene el código fuente y el flujo de trabajo estadístico completo desarrollado para mi **Trabajo de Fin de Grado (TFG) en Estadística** por la **Universidad de Sevilla**.

## 📋 Descripción del Proyecto
El objetivo principal de esta investigación es evaluar y comparar la robustez de una aproximación paramétrica tradicional (Modelo Lineal Generalizado - GLM) frente a dos algoritmos no paramétricos de *Machine Learning* (*Random Forest* y *Gradient Boosting Machine* - GBM) en la tarea de modelizar el nicho ecológico del pinsapo (*Abies pinsapo* Boiss.) y proyectar su vulnerabilidad geográfica ante escenarios de cambio climático con horizonte al año 2050.

### 🔍 Hallazgos Clave
* **Rendimiento:** Los modelos de ensamble superaron al GLM, alcanzando un empate técnico matemático (AUC > 0.96).
* **Validación espacial:** La auditoría mediante *Spatial Block CV* (bloques de 15x15 km) demostró que *Random Forest* ofrece la proyección cartográfica más conservadora y realista del nicho realizado, evitando los falsos refugios del GBM y GLM.
* **Transferibilidad:** El modelo demostró una robustez biológica excepcional al predecir de manera "ciega" (sin datos previos de entrenamiento) el refugio rifeño en Marruecos de su especie vicariante (*Abies marocana*).
* **Horizonte 2050:** Las proyecciones alertan de un colapso severo en sus refugios béticos actuales y dibujan un "espejismo ecológico" no colonizable en el noroeste peninsular.

## 🗂️ Estructura del Repositorio
* `/datos`: Carpeta destinada a almacenar la caché local de registros de GBIF y las capas raster de WorldClim.
* `Abies Pinsapo VC.R`: Script único en lenguaje R que automatiza todo el proceso de modelización y proyección espacial.

## 🛠️ Instalación y Uso Rápido

Para evitar tiempos de espera lentos con las APIs cartográficas y problemas de rutas locales, sigue estos pasos:

1. **Clona o descarga** este repositorio en tu ordenador.
2. **Descarga el dataset pre-procesado** (`datos.zip`) desde este enlace de Google Drive:  
   👉 **[PEGA AQUÍ TU ENLACE COPIADO DE GOOGLE DRIVE]**
3. **Descomprime** el archivo `datos.zip` en la misma carpeta donde hayas guardado el script. Asegúrate de que se quede una subcarpeta llamada `datos/` justo al lado del archivo de R.
4. **Abre el archivo `Abies Pinsapo VC.R`** en RStudio.
5. ⚠️ **¡Paso crucial! Configura el directorio de trabajo:** En el menú superior de RStudio, ve a **Session** -> **Set Working Directory** -> **Choose Directory...** y selecciona la carpeta principal que contiene tanto el script como la subcarpeta `datos`.

### Librerías de R necesarias
Asegúrate de tener instalados los siguientes paquetes en tu entorno antes de ejecutar el script:

```r
install.packages(c("rgbif", "dplyr", "geodata", "terra", "corrplot", 
                   "usdm", "caTools", "caret", "randomForest", "pdp", 
                   "ggplot2", "gridExtra", "gbm", "pROC", "blockCV", 
                   "sf", "rnaturalearth"))
