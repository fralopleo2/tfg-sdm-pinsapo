# 1. Cargar librerías necesarias
library(rgbif)
library(dplyr)
library(geodata)
library(terra)
library(corrplot)
library(usdm)
library(caTools)
library(caret)
library(randomForest)
library(pdp)
library(ggplot2)
library(gridExtra)
library(gbm)
library(pROC)
library(blockCV)
library(sf)
library(rnaturalearth)

# 2. Obtención del código identificador único del Pinsapo (Abies pinsapo)
especie <- name_backbone("Abies pinsapo")
id_especie <- especie$usageKey

# ==============================================================================
# EXTRACCIÓN DE DATOS BRUTOS
# ==============================================================================

# Definimos dónde se va a guardar nuestro archivo
archivo_gbif <- "datos/registros_pinsapo_gbif.rds"

if (file.exists(archivo_gbif)) {
  
  cat("\nCargando registros de GBIF desde el disco duro local...\n")
  df_bruto <- readRDS(archivo_gbif)
  
} else {
  
  cat("\nDescargando registros de GBIF de internet por primera vez...\n")
  # Descargamos los registros
  datos_gbif <- occ_data(taxonKey = id_especie, 
                         country = "ES", 
                         hasCoordinate = TRUE, 
                         limit = 20000)
  
  # Extraemos el dataframe
  df_bruto <- datos_gbif$data
  
  # Lo guardamos
  saveRDS(df_bruto, archivo_gbif)
  cat("Archivo guardado con éxito en:", archivo_gbif, "\n")
}

# Registros obtenidos:
nrow(df_bruto)

# ==============================================================================
# FASE 1: ANÁLISIS DEL DISEÑO Y CALIDAD DE LA BASE DE DATOS (METADATOS)
# ==============================================================================

analisis_metadatos <- df_bruto %>%
  # Agrupamos por la naturaleza del registro (Herbario vs. Observación humana)
  group_by(basisOfRecord) %>%
  
  # Calculamos métricas estadísticas de calidad para cada grupo
  summarise(
    Volumen_Total = n(),
    Falta_Año_pct = round(sum(is.na(year)) / n() * 100, 1),
    Incertidumbre_Media_m = round(mean(coordinateUncertaintyInMeters,
                                       na.rm = TRUE), 1),
    Año_Mas_Antiguo = min(year, na.rm = TRUE),
    Año_Mas_Reciente = max(year, na.rm = TRUE)
  ) %>%
  arrange(desc(Volumen_Total)) # Ordenamos de mayor a menor volumen

# Imprimimos la tabla comparativa en la consola
print(analisis_metadatos)

# ==============================================================================
# FASE 2: INGENIERÍA DE DATOS Y LIMPIEZA ESTADÍSTICA (PREPARACIÓN PARA MODELO)
# ==============================================================================

df_limpio <- df_bruto %>%
  # Seleccionamos solo las columnas útiles para el modelo y el análisis
  select(scientificName, decimalLongitude, decimalLatitude, year, 
         basisOfRecord, institutionCode, coordinateUncertaintyInMeters) %>%
  
  # Eliminamos registros históricos que no tengan año de muestreo
  filter(!is.na(year)) %>%
  
  # CONTROL DE CALIDAD ESPACIAL: 
  # Nos quedamos estrictamente con datos de altísima precisión (<= 500m)
  filter(coordinateUncertaintyInMeters <= 500) %>%
  
  # Filtramos la naturaleza del dato (evitamos material de jardinería...)
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN"))

# Registros obtenidos tras filtros de incertidumbre:
nrow(df_limpio)

# ==============================================================================
# AUDITORÍA DE LA MUESTRA FINAL
# ==============================================================================

composicion_final <- df_limpio %>%
  group_by(basisOfRecord) %>%
  summarise(
    Total = n(),
    Porcentaje = round((n() / nrow(df_limpio)) * 100, 1)
  ) %>%
  arrange(desc(Total))

# Tabla
print(composicion_final)

# ==============================================================================
# FASE 3: CRUCE DE COORDENADAS (GBIF) CON VARIABLES AMBIENTALES (WORLDCLIM)
# ==============================================================================

# 1. DESCARGA DEL CLIMA DE ESPAÑA
ruta_local = "datos"
clima_global <- worldclim_country(country = "Spain", 
                                  var = "bio", 
                                  res = 0.5, 
                                  path = ruta_local)
names(clima_global) <- paste0("BIO", 1:19)

# --- Variables de Temperatura (En °C) ---
# BIO1  : Temperatura Media Anual
# BIO2  : Oscilación Diurna Media (Media mensual de (Temp máx - Temp mín))
# BIO3  : Isotermalidad (BIO2 / BIO7 * 100) -> [Estabilidad térmica día
#                                               /noche vs año]
# BIO4  : Estacionalidad de la Temperatura (Desviación estándar * 100)
# BIO5  : Temperatura Máxima del Mes Más Cálido
# BIO6  : Temperatura Mínima del Mes Más Frío
# BIO7  : Oscilación Térmica Anual (BIO5 - BIO6)
# BIO8  : Temperatura Media del Trimestre Más Húmedo
# BIO9  : Temperatura Media del Trimestre Más Seco
# BIO10 : Temperatura Media del Trimestre Más Cálido
# BIO11 : Temperatura Media del Trimestre Más Frío

# --- Variables de Precipitación (En milímetros - mm) ---
# BIO12 : Precipitación Anual
# BIO13 : Precipitación del Mes Más Húmedo
# BIO14 : Precipitación del Mes Más Seco
# BIO15 : Estacionalidad de la Precipitación (Coeficiente de Variación)
#                                             [Ritmo de lluvias]
# BIO16 : Precipitación del Trimestre Más Húmedo
# BIO17 : Precipitación del Trimestre Más Seco
# BIO18 : Precipitación del Trimestre Más Cálido
# BIO19 : Precipitación del Trimestre Más Frío

# 2. DESCARGA DE ELEVACIÓN (ALTITUD)
elev_global <- elevation_30s(country = "Spain", path = ruta_local)
names(elev_global) <- "elevacion"

# Aseguramos que coincidan perfectamente en resolución y extensión
elev_global <- resample(elev_global, clima_global)
clima_global <- c(clima_global, elev_global)

# 3. CONVERSIÓN ESPACIAL
# Convertimos nuestro dataframe limpio en "puntos espaciales" 
puntos_pinsapo <- vect(df_limpio, 
                       geom = c("decimalLongitude", "decimalLatitude"), 
                       crs = "EPSG:4326") # Coordenadas GPS estándar

# 4. CRUCE (EXTRACCIÓN)
# Extraemos el clima exacto de cada una de las 467 coordenadas.
clima_pinsapo <- terra::extract(clima_global, puntos_pinsapo)

# Unimos los datos climáticos a nuestro dataframe original
df_modelo <- cbind(df_limpio, clima_pinsapo[, -1])
colnames(df_modelo)

# Comprobación rápida de que la extracción ha funcionado
head(df_modelo$BIO5)

# ==============================================================================
# FASE 4: EXTRACCIÓN TOTAL Y ANÁLISIS DE MULTICOLINEALIDAD (EDA)
# ==============================================================================

# 1. Limpieza de la matriz
# La función terra::extract() añade una primera columna llamada "ID" 
# que no es clima. La quitamos.
clima_vars <- clima_pinsapo[, -1]

# 2. Calculamos la Matriz de Correlación de Pearson
# Usamos use="complete.obs" por si se nos ha colado algún NA en el mar o bordes
matriz_correlacion <- cor(clima_vars, use = "complete.obs")

# 3. Visualización de la Matriz
corrplot(matriz_correlacion, 
         method = "color",       
         type = "upper",         
         tl.col = "black",       
         tl.srt = 45,            
         addCoef.col = "black",  
         number.cex = 0.5,       
         diag = FALSE,           
         title = "Matriz de Correlación de Pearson - Variables Bioclimáticas",
         mar = c(0,0,1,0))       

# ==============================================================================
# FASE 5: VIF AUTOMATIZADO SOBRE EL NICHO REALIZADO (FILTRADO)
# ==============================================================================

# 1. Recortamos Andalucía temporalmente con las 19 variables completas
limites_sur <- ext(-8.5, -1.5, 35.5, 39.0) 
clima_sur_completo <- crop(clima_global, limites_sur)

# 2. Filtro Espacial (Thinning): Máximo 1 punto por píxel de 1km
celdas_ocupadas <- terra::extract(clima_sur_completo[[1]], puntos_pinsapo, 
                                  cells = TRUE)$cell
celdas_unicas <- unique(na.omit(celdas_ocupadas))
coords_thin <- xyFromCell(clima_sur_completo[[1]], celdas_unicas)
puntos_presencia_limpios <- vect(coords_thin, crs = crs(clima_sur_completo))

# 3. Extraemos el clima EXACTO de los 101 pinsapos limpios
clima_pinsapos_limpios <- terra::extract(clima_sur_completo, 
                                         puntos_presencia_limpios)[, -1]
df_clima_pinsapo <- na.omit(as.data.frame(clima_pinsapos_limpios))

# 4. Cálculo del VIF dinámico
# Pasamos de th = 10 a th = 5 para máximo rigor ecológico
seleccion_vif <- vifstep(df_clima_pinsapo, th = 5)
clima_final_pinsapo <- exclude(df_clima_pinsapo, seleccion_vif)

# 5. Extraemos los nombres ganadores a la variable oro
variables_oro <- colnames(clima_final_pinsapo)

## VARIABLES SELECCIONADAS AUTOMÁTICAMENTE POR EL VIF
print(variables_oro)


# ==============================================================================
# FASE 6: PREPARACIÓN ESPACIAL Y RATIO 1:1
# ==============================================================================

# 1. Recortamos el mapa oficial usando solo las variables de oro
clima_sur <- subset(clima_sur_completo, variables_oro)

# ---- FIX PARA ELIMINAR EL MAR ---------------------
# Descargamos el polígono de España y nos quedamos solo con Andalucía
fronteras_espana <- gadm(country = "ESP", level = 1, path = ruta_local)
andalucia <- fronteras_espana[fronteras_espana$NAME_1 == "Andalucía", ]

# Enmascaramos el raster climático (todo lo que caiga fuera de Andalucía o
# en el mar será NA)
clima_sur <- mask(clima_sur, andalucia)
# ----------------------------------------------------

# 2. Pseudo-ausencias con Buffer de Exclusión y Ratio 1:1
buffer_presencias <- buffer(puntos_presencia_limpios, width = 2000) 
clima_para_ausencias <- mask(clima_sur, buffer_presencias, inverse = TRUE)

num_presencias <- nrow(coords_thin)

set.seed(123) 
puntos_ausencia <- spatSample(clima_para_ausencias[[1]], 
                              size = num_presencias, 
                              method = "random", 
                              na.rm = TRUE, 
                              as.points = TRUE)

## BALANCE DE CLASES
num_presencias

# 3. Construcción del Dataset Maestro Limpio
clima_ausencias <- terra::extract(clima_sur, puntos_ausencia)[, variables_oro] 
clima_ausencias$Presencia <- 0

clima_presencias <- df_clima_pinsapo[, variables_oro]
clima_presencias$Presencia <- 1

datos_modelo <- rbind(clima_presencias, clima_ausencias)
datos_modelo <- na.omit(datos_modelo)

# ==============================================================================
# FASE 7: ENTRENAMIENTO RANDOM FOREST CON VALIDACIÓN CRUZADA (K-FOLD CV)
# ==============================================================================

# 1. Preparamos los datos
datos_cv <- datos_modelo
datos_cv$Presencia <- factor(datos_cv$Presencia, 
                             levels = c(0, 1), 
                             labels = c("Ausencia", "Presencia"))

# 2. Configurar la Validación Cruzada (5 exámenes distintos)
control_cv <- trainControl(method = "cv", 
                           number = 5,                       
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,                
                           savePredictions = "final")

# 3. Entrenamiento
set.seed(123)
modelo_rf_cv <- train(Presencia ~ ., 
                      data = datos_cv, 
                      method = "rf", 
                      ntree = 500,
                      importance = TRUE,
                      trControl = control_cv, 
                      metric = "ROC")

# 4. Ver la nota final real
print(modelo_rf_cv)

# 5. Resultados
print(modelo_rf_cv$finalModel$confusion)
confusionMatrix(modelo_rf_cv)

# 6. Cálculo del TSS
mejor_fila_rf <- which.max(modelo_rf_cv$results$ROC)
sens_rf <- modelo_rf_cv$results$Sens[mejor_fila_rf]
spec_rf <- modelo_rf_cv$results$Spec[mejor_fila_rf]

tss_rf <- sens_rf + spec_rf - 1
round(tss_rf, 3)

## GRAFICOS PDP (Partial Dependence Plots)

# 1. Generar las curvas individuales guardándolas en variables (p1, p2, p3)
p1 <- partial(modelo_rf_cv, pred.var = "BIO3", prob = TRUE, plot = TRUE,
              which.class = 2, plot.engine = "ggplot2", ylab = "Probabilidad", 
              main = "BIO3 (Isotermalidad)")

p2 <- partial(modelo_rf_cv, pred.var = "BIO15", prob = TRUE, plot = TRUE, 
              which.class = 2, plot.engine = "ggplot2", ylab = "Probabilidad",
              main = "BIO15 (Estacionalidad)")

p3 <- partial(modelo_rf_cv, pred.var = "BIO8", prob = TRUE, plot = TRUE,
              which.class = 2, plot.engine = "ggplot2", ylab = "Probabilidad",
              main = "BIO8 (Temp. Húmedo)")

# 2. Mostrar los tres gráficos juntos
grid.arrange(p1, p2, p3, ncol = 3)

# ==============================================================================
# FASE 8: ENTRENAMIENTO GLM CON VALIDACIÓN CRUZADA (K-FOLD CV)
# ==============================================================================

set.seed(123) 

modelo_glm_cv <- train(Presencia ~ ., 
                       data = datos_cv,        
                       method = "glm",         
                       family = "binomial",    
                       trControl = control_cv, 
                       metric = "ROC")

# 1. Ver la nota final del GLM
print(modelo_glm_cv)

# 2. Ver la matriz de confusión promedio
confusionMatrix(modelo_glm_cv)

# 3. Comprobar la Devianza en consola
summary(modelo_glm_cv$finalModel)

# 4. Cálculo del TSS
mejor_fila_glm <- which.max(modelo_glm_cv$results$ROC)
sens_glm <- modelo_glm_cv$results$Sens[mejor_fila_glm]
spec_glm <- modelo_glm_cv$results$Spec[mejor_fila_glm]

tss_glm <- sens_glm + spec_glm - 1
round(tss_glm, 3)

# 5. Generar el gráfico con los 4 diagnósticos
par(mfrow = c(2, 2)) 
plot(modelo_glm_cv$finalModel)  
par(mfrow = c(1, 1)) 

# ==============================================================================
# FASE 9: ENTRENAMIENTO GRADIENT BOOSTING (GBM) CON K-FOLD CV
# ==============================================================================

# Usamos la misma semilla para que los exámenes sean idénticos al RF y GLM
set.seed(123) 

modelo_gbm_cv <- train(Presencia ~ ., 
                       data = datos_cv, 
                       method = "gbm", 
                       trControl = control_cv, 
                       metric = "ROC",
                       verbose = FALSE) 

# 1. Ver la nota final del GBM
print(modelo_gbm_cv)

# 2. Ver la matriz de confusión promedio
confusionMatrix(modelo_gbm_cv)

# 3. Cálculo del TSS
mejor_fila_gbm <- which.max(modelo_gbm_cv$results$ROC)
sens_gbm <- modelo_gbm_cv$results$Sens[mejor_fila_gbm]
spec_gbm <- modelo_gbm_cv$results$Spec[mejor_fila_gbm]

tss_gbm <- sens_gbm + spec_gbm - 1
round(tss_gbm, 3)

## GRAFICOS PDP

# 1. Generar las curvas individuales guardándolas en variables (p1, p2, p3)
p1_gbm <- partial(modelo_gbm_cv, pred.var = "BIO3", prob = TRUE, plot = TRUE, 
                  which.class = 2, plot.engine = "ggplot2", 
                  ylab = "Probabilidad", main = "BIO3 (Isotermalidad)")

p2_gbm <- partial(modelo_gbm_cv, pred.var = "BIO15", prob = TRUE, plot = TRUE, 
                  which.class = 2, plot.engine = "ggplot2", 
                  ylab = "Probabilidad", main = "BIO15 (Estacionalidad)")

p3_gbm <- partial(modelo_gbm_cv, pred.var = "BIO8", prob = TRUE, plot = TRUE, 
                  which.class = 2, plot.engine = "ggplot2", 
                  ylab = "Probabilidad", main = "BIO8 (Temp. Húmedo)")

# 2. Mostrar los tres gráficos juntos
grid.arrange(p1_gbm, p2_gbm, p3_gbm, ncol = 3)


# ==============================================================================
# FASE 10: VISUALIZACIÓN COMPARATIVA (RF vs GLM vs GBM)
# ==============================================================================

# 1. CÁLCULO DE CURVAS ROC
roc_rf <- roc(modelo_rf_cv$pred$obs, modelo_rf_cv$pred$Presencia, 
              quiet = TRUE)
roc_glm <- roc(modelo_glm_cv$pred$obs, modelo_glm_cv$pred$Presencia,
               quiet = TRUE)
roc_gbm <- roc(modelo_gbm_cv$pred$obs, modelo_gbm_cv$pred$Presencia, 
               quiet = TRUE)

# Usamos max() para coger el mejor resultado
auc_rf <- auc(roc_rf)
auc_glm <- auc(roc_glm)
auc_gbm <- auc(roc_gbm)

# 2. COMPARACIÓN DE CURVAS ROC (3 MODELOS)
plot(roc_rf, col = "forestgreen", lwd = 2, 
     main = "Comparativa de Modelos (5-fold CV)")
plot(roc_gbm, col = "darkorange", lwd = 2, add = TRUE)
plot(roc_glm, col = "steelblue", lwd = 2, add = TRUE)

# 2. Leyenda
legend("bottomright", 
       legend = c(
         paste0("GBM (AUC: ", round(auc_gbm,3), ")"), 
         paste0("GLM (AUC: ", round(auc_glm,3), ")"),
         paste0("Random Forest (AUC: ", round(auc_rf,3), ")")
       ), 
       col = c("darkorange", "steelblue", "forestgreen"), 
       lwd = 2, 
       cex = 0.8)

# 3. IMPORTANCIA DE VARIABLES (Basado en RF)
# Extraemos la importancia unificada usando caret
importancia_rf <- varImp(modelo_rf_cv, scale = FALSE)

valores_rf <- importancia_rf$importance[, 1]
nombres_rf <- rownames(importancia_rf$importance)

datos_barras_rf <- setNames(valores_rf, nombres_rf)
datos_barras_rf <- sort(datos_barras_rf, decreasing = FALSE)

par(mar = c(5, 5, 4, 2)) 
barplot(datos_barras_rf, 
        horiz = TRUE,          
        las = 1,               
        col = "forestgreen",    
        border = "black",      
        main = "Importancia de Variables (Random Forest)", 
        xlab = "Nivel de Importancia")

# 4. IMPORTANCIA DE VARIABLES (Basado en GLM)
importancia_glm <- varImp(modelo_glm_cv, scale = FALSE)

valores_glm <- importancia_glm$importance[, 1]
nombres_glm <- rownames(importancia_glm$importance)

datos_barras_glm <- setNames(valores_glm, nombres_glm)
datos_barras_glm <- sort(datos_barras_glm, decreasing = FALSE)

par(mar = c(5, 5, 4, 2)) 
barplot(datos_barras_glm, 
        horiz = TRUE,          
        las = 1,               
        col = "steelblue",      
        border = "black",      
        main = "Importancia de Variables (Regresión Logística - GLM)", 
        xlab = "Valor Z absoluto (Fuerza de la variable)")

# 5. IMPORTANCIA DE VARIABLES (Basado en GBM)
importancia_gbm <- varImp(modelo_gbm_cv, scale = FALSE)

valores_gbm <- importancia_gbm$importance[, 1]
nombres_gbm <- rownames(importancia_gbm$importance)

datos_barras_gbm <- setNames(valores_gbm, nombres_gbm)
datos_barras_gbm <- sort(datos_barras_gbm, decreasing = FALSE)

par(mar = c(5, 5, 4, 2)) 
barplot(datos_barras_gbm, 
        horiz = TRUE,          
        las = 1,               
        col = "darkorange",     
        border = "black",      
        main = "Importancia de Variables (Gradient Boosting - GBM)", 
        xlab = "Nivel de Importancia")

par(mar = c(5, 4, 4, 2) + 0.1)

# ==============================================================================
# FASE 11: PROYECCIÓN ESPACIAL DEL NICHO ECOLÓGICO (MAPA FINAL GLM)
# ==============================================================================

# Con caret, type = "prob" devuelve 2 valores (probabilidad de Ausencia 
# y de Presencia). Con index = 2 le decimos que mapee 
# la probabilidad de PRESENCIA.
mapa_prob_glm <- predict(clima_sur, modelo_glm_cv, type = "prob", index = 2, 
                         na.rm = TRUE)

# 1. DISEÑO VISUA
paleta_colores <- colorRampPalette(c("#e0e0e0", "#abdda4", "#ffffbf", 
                                     "#fdae61", "#d7191c"))(100)

# 2. DIBUJAMOS EL MAPA BASE
plot(mapa_prob_glm, 
     main = "Hábitat Potencial de Abies pinsapo (GLM)",
     col = paleta_colores,
     plg = list(title = "Idoneidad\n(0 a 1)"))

# 3. SUPERPONEMOS LA REALIDAD
plot(puntos_pinsapo, add = TRUE, col = "black", pch = 16, cex = 0.4)

# ==============================================================================
# FASE 12: PROYECCIÓN ESPACIAL DEL NICHO ECOLÓGICO (MAPA FINAL RANDOM FOREST)
# ==============================================================================

# Proyección idéntica a la del GLM pero con el modelo RF
mapa_prob_rf <- predict(clima_sur, modelo_rf_cv, type = "prob", index = 2,
                        na.rm = TRUE)

# 1. DIBUJAMOS EL MAPA BASE
plot(mapa_prob_rf, 
     main = "Hábitat Potencial de Abies pinsapo (Random Forest)",
     col = paleta_colores,
     plg = list(title = "Idoneidad\n(0 a 1)"))

# 2. SUPERPONEMOS LA REALIDAD
plot(puntos_pinsapo, add = TRUE, col = "black", pch = 16, cex = 0.4)

# ==============================================================================
# FASE 13: PROYECCIÓN ESPACIAL DEL NICHO ECOLÓGICO (MAPA FINAL GBM)
# ==============================================================================

# Proyección idéntica a la del GLM pero con el modelo RF
mapa_prob_gbm <- predict(clima_sur, modelo_gbm_cv, type = "prob", index = 2,
                         na.rm = TRUE)

# 1. DIBUJAMOS EL MAPA BASE
plot(mapa_prob_gbm, 
     main = "Hábitat Potencial de Abies pinsapo (GBM)",
     col = paleta_colores,
     plg = list(title = "Idoneidad\n(0 a 1)"))

# 2. SUPERPONEMOS LA REALIDAD
plot(puntos_pinsapo, add = TRUE, col = "black", pch = 16, cex = 0.4)

# ==============================================================================
# FASE 14: Spatial Block Cross-Validation
# ==============================================================================

# RECUPERAR COORDENADAS Y UNIRLAS A LOS DATOS
# Extraemos las coordenadas (lon, lat) de los SpatVector
coords_presencia <- crds(puntos_presencia_limpios)
coords_ausencia <- crds(puntos_ausencia)

# Las juntamos una encima de otra (primero 101 pres, luego 101 aus)
todas_las_coords <- rbind(coords_presencia, coords_ausencia)
coords_df <- as.data.frame(todas_las_coords)
colnames(coords_df) <- c("lon", "lat")

# Pegamos las coordenadas
datos_completos <- cbind(coords_df, datos_cv)

# 1. CONVERTIR DATOS A FORMATO ESPACIAL
datos_sf <- st_as_sf(datos_completos, coords = c("lon", "lat"), crs = 4326)

# 2. GENERAR LOS BLOQUES ESPACIALES (Cuadrículas de 15 km)
set.seed(123)
bloques <- cv_spatial(x = datos_sf,
                      column = "Presencia", 
                      k = 5,                
                      size = 15000,         
                      iteration = 50)       

# 1. Extraer SOLO los índices de entrenamiento (train) de los bloques
# espaciales
indices_espaciales <- lapply(bloques$folds_list, function(x) x[[1]])
names(indices_espaciales) <- paste0("Fold", 1:5)

# 2. Re-crear el control espacial con la lista limpia
control_espacial <- trainControl(method = "cv",
                                 index = indices_espaciales,
                                 classProbs = TRUE,
                                 summaryFunction = twoClassSummary)

# 3. Re-entrenar los modelos
formula_dinamica <- as.formula(paste("Presencia ~", paste(variables_oro, 
                                                          collapse = " + ")))
modelo_glm_esp <- train(formula_dinamica,
                        data = datos_completos,
                        method = "glm",
                        family = "binomial",
                        metric = "ROC",
                        trControl = control_espacial)

modelo_rf_esp <- train(formula_dinamica,
                       data = datos_completos,
                       method = "rf",
                       metric = "ROC",
                       trControl = control_espacial)

modelo_gbm_esp <- train(formula_dinamica,
                        data = datos_completos,
                        method = "gbm",
                        metric = "ROC",
                        trControl = control_espacial,
                        verbose = FALSE)

# 4. Ver los resultados

modelo_glm_esp$results[which.max(modelo_glm_esp$results$ROC), ]
modelo_rf_esp$results[which.max(modelo_rf_esp$results$ROC), ]
modelo_gbm_esp$results[which.max(modelo_gbm_esp$results$ROC), ]

## MAPA DE LOS BLOQUES ESPACIALES

# 1. Descargar el mapa base de España
espana <- ne_countries(scale = "medium", country = "Spain", returnclass = "sf")

# 2. Asignamos el número de bloque directamente a cada punto
datos_sf$Fold <- as.factor(bloques$folds_ids)

# 3. Creamos el mapa
mapa_ev <- ggplot() +
  # Fondo de Andalucía
  geom_sf(data = espana, fill = "gray98", color = "gray70") +
  
  # Dibujamos los puntos: Color = Pliegue (Fold), Forma = Presencia/Ausencia
  geom_sf(data = datos_sf, aes(color = Fold, shape = Presencia), size = 2.5,
          alpha = 0.9) +
  
  # Zoom a Andalucía
  coord_sf(xlim = c(-8.5, -1.5), ylim = c(35.5, 39.0)) +
  
  # Colores
  scale_color_brewer(palette = "Set1", name = "Bloque Espacial") +
  scale_shape_manual(values = c("Ausencia" = 4, "Presencia" = 16)) + 
  
  theme_minimal() +
  labs(title = "Validación Cruzada por Bloques Espaciales (Spatial CV)",
       subtitle = "Agrupación geográfica de las muestras para
       evitar autocorrelación",
       x = "Longitud", y = "Latitud") +
  theme(legend.position = "right",
        panel.background = element_rect(fill = "#eef7fa", color = NA), 
        panel.grid.major = element_line(color = "white")) 

# Mostrar el mapa
print(mapa_ev)

# ==============================================================================
# FASE 15: EXTRAPOLACIÓN PENINSULAR DEL NICHO (Refugios climáticos)
# ==============================================================================

# 1. DESCARGAR Y FUSIONAR ESPAÑA Y PORTUGAL

frontera_esp <- gadm(country = "ESP", level = 0, path = ruta_local)
frontera_prt <- gadm(country = "PRT", level = 0, path = ruta_local)

# Unimos los países para el mapa principal
macro_region <- rbind(frontera_esp, frontera_prt)

# Creamos una máscara de corte (extent) centrada en Península
limites_mapa_principal <- ext(-10, 5, 35.8, 44) 

# 2. Hacemos un corte para reducir el tamaño a nuestra zona.
clima_actual_corte <- crop(clima_global, limites_mapa_principal)

# 3. Enmascaramos con los países (para quitar el mar)
clima_actual_macro <- mask(clima_actual_corte, macro_region)

# 4. PREDICCIONES
mapa_actual_rf <- predict(subset(clima_actual_macro, variables_oro),
                              modelo_rf_esp, type = "prob", index = 2,
                          na.rm = TRUE)


# Mapa: Estado Actual
plot(mapa_actual_rf,main = "Extrapolación del Hábitat Potencial",
     col = paleta_colores,
     plg = list(title = "Idoneidad", x = "right"),
     axes = FALSE)

### MAPA PRESENCIAS GBIF

par(mar = c(1, 1, 3, 1))

plot(0, 0, type = "n", 
     xlim = c(xmin(limites_mapa_principal), xmax(limites_mapa_principal)), 
     ylim = c(ymin(limites_mapa_principal), ymax(limites_mapa_principal)), 
     axes = FALSE,
     xlab = "", ylab = "",
     main = "Distribución Observada (Registros de Muestreo)",
     asp = 1.3)

plot(macro_region, 
     col = "white",      
     border = "gray88",
     axes = FALSE,       
     add = TRUE)         

coords_puntos <- crds(puntos_pinsapo)
points(coords_puntos[,1], coords_puntos[,2], 
       col = "black", 
       pch = 16,
       cex = 0.1)

par(mar = c(5, 4, 4, 2) + 0.1)

# ==============================================================================
# FASE 16: TRANSFERIBILIDAD ESPACIAL (MARRUECOS Y LA CORDILLERA DEL RIF)
# ==============================================================================

# Recuperamos el polígono de España

clima_nacional <- subset(clima_global, variables_oro)

# 1. Descargamos las fronteras y ajustamos al Rif
fronteras_marruecos <- gadm(country = "MAR", level = 0, path = ruta_local)
clima_marruecos_norte <- crop(clima_nacional, ext(-6.5, -3.0, 34.0, 36.0))
clima_marruecos_norte <- mask(clima_marruecos_norte, fronteras_marruecos)

# 2. PREDECIR usando el modelo RF espacial
mapa_marruecos_rf <- predict(clima_marruecos_norte, modelo_rf_esp,
                             type = "prob", index = 2, na.rm = TRUE)

# 3. DIBUJAMOS EL MAPA
par(mar = c(3, 3, 4, 6)) 

plot(mapa_marruecos_rf, 
     main = "Transferibilidad: Nicho Climático en Marruecos",
     col = paleta_colores,
     plg = list(title = "Idoneidad\n(0 a 1)", x = "right"), 
     axes = TRUE)

par(mar = c(5, 4, 4, 2) + 0.1)

# ==============================================================================
# FASE 17: PROYECCIÓN CAMBIO CLIMÁTICO 2050 (MÁXIMA RESOLUCIÓN 0.5)
# ==============================================================================

## DESCARGA DE LOS ESCENARIOS CLIMÁTICOS
# ESCENARIO MODERADO (SSP2-4.5)
# clima_245_global <- cmip6_world(model = "MPI-ESM1-2-HR", ssp = "245", 
#                                time = "2041-2060", var = "bioc", 
#                                res = 0.5, path = ruta_local)

# ESCENARIO PESIMISTA (SSP5-8.5)
# clima_585_global <- cmip6_world(model = "MPI-ESM1-2-HR", ssp = "585", 
#                                time = "2041-2060", var = "bioc", 
#                                res = 0.5, path = ruta_local)

# Definimos las rutas a los archivos
ruta_245 <- "datos/cmip6/30s/MPI-ESM1-2-HR/ssp245/
              wc2.1_30s_bioc_MPI-ESM1-2-HR_ssp245_2041-2060.tif"
ruta_585 <- "datos/cmip6/30s/MPI-ESM1-2-HR/ssp585/
              wc2.1_30s_bioc_MPI-ESM1-2-HR_ssp585_2041-2060.tif"

# Usamos rast() para cargar los archivos
clima_245_global <- rast(ruta_245)
clima_585_global <- rast(ruta_585)

# Nombramos las capas como lo hace WorldClim (BIO1, BIO2...)
names(clima_245_global) <- paste0("BIO", 1:19)
names(clima_585_global) <- paste0("BIO", 1:19)

# 1. Recorte y Enmascarado
clima_245_sur <- mask(crop(clima_245_global, limites_sur), andalucia)
clima_585_sur <- mask(crop(clima_585_global, limites_sur), andalucia)

# 2. Predicción con Random Forest
mapa_245_rf <- predict(subset(clima_245_sur, variables_oro), modelo_rf_esp, 
                       type = "prob", index = 2, na.rm = TRUE)

mapa_585_rf <- predict(subset(clima_585_sur, variables_oro), modelo_rf_esp, 
                       type = "prob", index = 2, na.rm = TRUE)

# ==============================================================================
# FASE 18: PROYECCIÓN A NIVEL DE LA PENÍNSULA IBÉRICA (BIOGEOGRAFÍA REAL)
# ==============================================================================

frontera_mar <- gadm(country = "MAR", level = 0, path = ruta_local)
macro_region_2 <- rbind(frontera_esp, frontera_prt, frontera_mar)
limites_mapa_cambio <- ext(-10, 5, 34.0, 44) 

# 1. Hacemos un corte
clima_245_corte    <- crop(clima_245_global, limites_mapa_cambio)
clima_585_corte    <- crop(clima_585_global, limites_mapa_cambio)

# 2. Enmascaramos con los países
clima_245_macro    <- mask(clima_245_corte, macro_region_2)
clima_585_macro    <- mask(clima_585_corte, macro_region_2)

# 3. PREDICCIONES
mapa_245_mac_rf    <- predict(subset(clima_245_macro, variables_oro),
                              modelo_rf_esp, type = "prob", index = 2,
                              na.rm = TRUE)
mapa_585_mac_rf    <- predict(subset(clima_585_macro, variables_oro),
                              modelo_rf_esp, type = "prob", index = 2,
                              na.rm = TRUE)

# DIBUJADO DE MAPAS LIMPIO

par(mfrow = c(1, 1))

# Mapa 1: 2050 Moderado
plot(mapa_245_mac_rf, main = "2050 Moderado (SSP2-4.5)", col = paleta_colores,
     axes = FALSE, range = c(0, 1))

# Mapa 2: 2050 Pesimista
plot(mapa_585_mac_rf, main = "2050 Pesimista (SSP5-8.5)", col = paleta_colores,
     axes = FALSE, range = c(0, 1))

# ==============================================================================
# COMPARATIVA DE OVERFITTING (RF, GBM y GLM)
# ==============================================================================

# 1. Extraemos todas las notas (AUC) de los exámenes
# Random Forest
auc_rf_ran <- modelo_rf_cv$resample$ROC
auc_rf_esp <- modelo_rf_esp$resample$ROC
# GBM
auc_gbm_ran <- modelo_gbm_cv$resample$ROC
auc_gbm_esp <- modelo_gbm_esp$resample$ROC
# GLM
auc_glm_ran <- modelo_glm_cv$resample$ROC
auc_glm_esp <- modelo_glm_esp$resample$ROC

# 2. Construimos la tabla completa
datos_multi_overfitting <- data.frame(
  AUC = c(auc_rf_ran, auc_rf_esp, auc_gbm_ran, auc_gbm_esp, auc_glm_ran,
          auc_glm_esp),
  
  Modelo = factor(rep(c("Random Forest", "Gradient Boosting (GBM)",
                        "Regresión Logística (GLM)"), each = 10),
                  levels = c("Random Forest", "Gradient Boosting (GBM)",
                             "Regresión Logística (GLM)")),
  
  Validacion = rep(rep(c("1. Aleatoria", "2. Espacial"), each = 5), 3)
)

grafico_multi <- ggplot(datos_multi_overfitting, aes(x = Validacion,
                                                     y = AUC,
                                                     fill = Validacion)) +
  geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA,
               color = "darkgrey") +
  geom_jitter(width = 0.1, size = 3, color = "#2c3e50", alpha = 0.8) + 
  
  # Separamos en 3 paneles
  facet_wrap(~Modelo) +
  
  scale_fill_manual(values = c("1. Aleatoria" = "#e74c3c", 
                               "2. Espacial" = "#2ecc71")) +
  theme_bw() + # Fondo blanco con bordes limpios
  labs(title = "Impacto de la Autocorrelación Espacial por Algoritmo",
       subtitle = "Comparativa del AUC (5-fold) entre validación Aleatoria
                    y Bloques Espaciales de 15km",
       y = "Rendimiento (AUC)", 
       x = "") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 11, color = "white"),
        strip.background = element_rect(fill = "#2c3e50"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Mostrar el gráfico
print(grafico_multi)
