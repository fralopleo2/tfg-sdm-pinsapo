## 🛠️ Instalación y Uso Rápido

Para evitar tiempos de espera lentos con las APIs cartográficas y problemas de rutas, sigue estos pasos:

1. **Clona o descarga** este repositorio en tu ordenador.
2. **Descarga el dataset pre-procesado** (`datos.zip`) desde este enlace de Google Drive:  
   👉 **[PEGA AQUÍ TU ENLACE COPIADO DE GOOGLE DRIVE]**
3. **Descomprime** el archivo `datos.zip` en la misma carpeta donde hayas guardado el script. Asegúrate de que se quede una subcarpeta llamada `datos/` justo al lado del archivo de R.
4. **Abre el archivo `Abies Pinsapo VC.R`** en RStudio.
5. ⚠️ **¡Paso crucial! Configura el directorio de trabajo:** En el menú superior de RStudio, ve a **Session** -> **Set Working Directory** -> **Choose Directory...** y selecciona la carpeta principal que contiene tanto el script como la subcarpeta `datos`.

### Librerías de R necesarias
Asegúrate de tener instalados los siguientes paquetes antes de ejecutar el script:

```r
install.packages(c("rgbif", "dplyr", "geodata", "terra", "corrplot", 
                   "usdm", "caTools", "caret", "randomForest", "pdp", 
                   "ggplot2", "gridExtra", "gbm", "pROC", "blockCV", 
                   "sf", "rnaturalearth"))
                   "ggplot2", "gridExtra", "gbm", "pROC", "blockCV", 
                   "sf", "rnaturalearth"))
