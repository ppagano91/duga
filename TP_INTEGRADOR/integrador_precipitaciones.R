## Integrador ADE-DUGA 2025 - Precipitaciones
## Flujo basado en scripts de clase M3U3 (geoR/gstat/sp/sf/raster)

## Librerías requeridas
suppressPackageStartupMessages({
  library(readxl)
  library(geoR)
  library(sp)
  library(sf)
  library(raster)
})

## Rutas y configuraciones
output_dir <- "/workspace/TP_INTEGRADOR"
data_path  <- file.path(output_dir, "Precipitaciones.xlsx")

if (!file.exists(data_path)) {
  stop(sprintf("No se encuentra el archivo de datos: %s", data_path))
}

## Lectura de datos
prec <- read_excel(data_path)

required_cols <- c("longitud", "latitud", "precipitacion")
missing_cols <- setdiff(required_cols, names(prec))
if (length(missing_cols) > 0) {
  stop(sprintf("Faltan columnas requeridas en el Excel: %s", paste(missing_cols, collapse = ", ")))
}

## Limpieza básica
prec <- prec[!is.na(prec$longitud) & !is.na(prec$latitud) & !is.na(prec$precipitacion), ]
if (nrow(prec) < 5) {
  stop("Muy pocos registros válidos tras filtrar NAs. Se requieren al menos 5 puntos.")
}

## Estadística descriptiva básica y gráficos exploratorios
summary_stats <- summary(prec$precipitacion)
print(summary_stats)

try({
  png(filename = file.path(output_dir, "01_descriptiva_hist_box_precipitacion.png"), width = 1600, height = 900, res = 150)
  par(mfrow = c(1, 2))
  hist(prec$precipitacion, main = "Histograma de Precipitación", xlab = "precipitacion", col = "lightblue", border = "gray40")
  boxplot(prec$precipitacion, main = "Boxplot de Precipitación", horizontal = TRUE, col = "lightgreen")
  dev.off()
}, silent = TRUE)

## Conversión a geoR::geodata
prec_geo <- as.geodata(prec, coords.col = c("longitud", "latitud"), data.col = "precipitacion")

## Variograma experimental
vario <- variog(prec_geo)

try({
  png(filename = file.path(output_dir, "02_variograma_experimental.png"), width = 1400, height = 900, res = 150)
  plot(vario, main = "Semivariograma experimental: precipitacion")
  dev.off()
}, silent = TRUE)

## Ajuste de modelo de variograma (REML) con iniciales heurísticas
## Heurística simple: pepita ~ min(2*gamma en primer lag), psill ~ var(y)*0.5, rango ~ 1/4 del max rango observado (escala relativa)
y_var <- var(prec$precipitacion, na.rm = TRUE)
ini_psill <- max(y_var * 0.5, 1e-6)
ini_range <- 0.25
ini_nug   <- max(1e-6, if (length(vario$u) > 0 && length(vario$v) > 0) {
  if (!is.na(vario$v[1])) 2 * vario$v[1] else 1e-6
} else 1e-6)

## Modelos candidatos (usar el de clase: "wave"). Si falla, intentar exponencial.
fit_reml <- NULL
fit_ml   <- NULL
fit_ok   <- FALSE

## Intento 1: "wave"
try({
  fit_reml <- likfit(prec_geo,
                      cov.model = "wave",
                      ini.cov.pars = c(ini_psill, ini_range),
                      lik.method = "REML",
                      nugget = ini_nug,
                      messages = FALSE)
  fit_ml <- likfit(prec_geo,
                   cov.model = "wave",
                   ini.cov.pars = c(ini_psill, ini_range),
                   lik.method = "ML",
                   nugget = ini_nug,
                   messages = FALSE)
  fit_ok <- TRUE
}, silent = TRUE)

## Intento 2: "exponential" si falla el anterior
if (!fit_ok) {
  try({
    fit_reml <- likfit(prec_geo,
                        cov.model = "exponential",
                        ini.cov.pars = c(ini_psill, ini_range),
                        lik.method = "REML",
                        nugget = ini_nug,
                        messages = FALSE)
    fit_ml <- likfit(prec_geo,
                     cov.model = "exponential",
                     ini.cov.pars = c(ini_psill, ini_range),
                     lik.method = "ML",
                     nugget = ini_nug,
                     messages = FALSE)
    fit_ok <- TRUE
  }, silent = TRUE)
}

if (!fit_ok) {
  stop("No fue posible ajustar un modelo de variograma (ni wave ni exponential). Revise los datos.")
}

print(fit_ml)
print(summary(fit_ml))
print(fit_reml)
print(summary(fit_reml))

## Graficar variograma con ajustes ML/REML
try({
  png(filename = file.path(output_dir, "03_variograma_ajustes.png"), width = 1400, height = 900, res = 150)
  plot(vario, pch = 20, col = "red", main = "Semivariograma experimental y ajustes")
  lines(fit_ml, col = "blue",  lwd = 2)
  lines(fit_reml, col = "darkgreen", lwd = 2)
  legend("topleft", legend = c("ML", "REML"), lty = c(1,1), col = c("blue", "darkgreen"), bty = "n")
  dev.off()
}, silent = TRUE)

## Grilla de predicción
pred_grilla <- expand.grid(
  x = seq(min(prec_geo$coords[, 1]), max(prec_geo$coords[, 1]), length.out = 100),
  y = seq(min(prec_geo$coords[, 2]), max(prec_geo$coords[, 2]), length.out = 100)
)

## Kriging con REML
ko_reml <- krige.conv(prec_geo, locations = pred_grilla, krige = krige.control(obj.model = fit_reml))

## A objetos espaciales
prec_sp <- prec
coordinates(prec_sp) <- ~ longitud + latitud

ko_reml_sp <- SpatialPixelsDataFrame(points = pred_grilla, data = data.frame(ko_reml[1:2]))
names(ko_reml_sp@data) <- c("predict", "krige.var")

## Mapas de predicción y de varianza (PNG)
try({
  png(filename = file.path(output_dir, "04_map_prediccion_precipitacion.png"), width = 1400, height = 1200, res = 150)
  spplot(ko_reml_sp, zcol = "predict",
         sp.layout = list("sp.points", prec_sp, col = "black", pch = 20),
         col.regions = heat.colors(100), contour = TRUE,
         main = "Predicciones de precipitación",
         xlab = "Longitud", ylab = "Latitud")
  dev.off()
}, silent = TRUE)

try({
  png(filename = file.path(output_dir, "05_map_varianza_precipitacion.png"), width = 1400, height = 1200, res = 150)
  spplot(ko_reml_sp, zcol = "krige.var",
         sp.layout = list("sp.points", prec_sp, col = "black", pch = 20),
         col.regions = rev(heat.colors(100)), contour = FALSE,
         main = "Varianza de predicción (precipitación)",
         xlab = "Longitud", ylab = "Latitud",
         key.space = "right")
  dev.off()
}, silent = TRUE)

## Exportación a archivos para QGIS
## 1) Polígono de área de predicción (bounding box)
min_lon_grid <- min(pred_grilla$x, na.rm = TRUE)
max_lon_grid <- max(pred_grilla$x, na.rm = TRUE)
min_lat_grid <- min(pred_grilla$y, na.rm = TRUE)
max_lat_grid <- max(pred_grilla$y, na.rm = TRUE)

bbox_coords_grid <- matrix(
  c(min_lon_grid, min_lat_grid,
    max_lon_grid, min_lat_grid,
    max_lon_grid, max_lat_grid,
    min_lon_grid, max_lat_grid,
    min_lon_grid, min_lat_grid),
  ncol = 2,
  byrow = TRUE
)

bbox_polygon_sfc_grid <- st_sfc(st_polygon(list(bbox_coords_grid)), crs = 4326)
bbox_sf_grid <- st_sf(geometry = bbox_polygon_sfc_grid)

st_write(bbox_sf_grid, dsn = file.path(output_dir, "area_prediccion_poligono.shp"), driver = "ESRI Shapefile", append = FALSE)
cat("area_prediccion_poligono.shp OK\n")

## 2) Puntos de muestreo
puntos_muestreo_sf <- st_as_sf(prec, coords = c("longitud", "latitud"), crs = 4326)
st_write(puntos_muestreo_sf, dsn = file.path(output_dir, "puntos_muestreo.shp"), driver = "ESRI Shapefile", append = FALSE)
cat("puntos_muestreo.shp OK\n")

## 3) Curvas de nivel (contornos) de la superficie de predicción
x_grid_unique <- unique(pred_grilla$x)
y_grid_unique <- unique(pred_grilla$y)
z_matrix_predict <- matrix(ko_reml$predict, nrow = length(x_grid_unique), ncol = length(y_grid_unique))

contour_lines_predict_data <- contourLines(x = x_grid_unique, y = y_grid_unique, z = z_matrix_predict)

if (length(contour_lines_predict_data) > 0) {
  lines_sf_predict_list <- vector("list", length(contour_lines_predict_data))
  levels_predict_list <- numeric(length(contour_lines_predict_data))
  for (i in seq_along(contour_lines_predict_data)) {
    coords <- cbind(contour_lines_predict_data[[i]]$x, contour_lines_predict_data[[i]]$y)
    lines_sf_predict_list[[i]] <- st_linestring(coords)
    levels_predict_list[i] <- contour_lines_predict_data[[i]]$level
  }
  all_lines_sfc_predict <- st_sfc(lines_sf_predict_list, crs = 4326)
  contours_predict_sf <- st_sf(level = levels_predict_list, geometry = all_lines_sfc_predict)
  st_write(contours_predict_sf, dsn = file.path(output_dir, "contornos_prediccion.shp"), driver = "ESRI Shapefile", append = FALSE)
  cat("contornos_prediccion.shp OK\n")
} else {
  warning("No se generaron contornos de predicción (posible superficie plana o datos insuficientes).")
}

## 4) Rasters GeoTIFF (predicción y varianza)
raster_predicciones_sp <- raster(ko_reml_sp["predict"])
crs(raster_predicciones_sp) <- CRS("+init=epsg:4326")
writeRaster(raster_predicciones_sp, filename = file.path(output_dir, "predicciones_precipitacion.tif"), format = "GTiff", overwrite = TRUE)
cat("predicciones_precipitacion.tif OK\n")

raster_varianzas_sp <- raster(ko_reml_sp["krige.var"])
crs(raster_varianzas_sp) <- CRS("+init=epsg:4326")
writeRaster(raster_varianzas_sp, filename = file.path(output_dir, "varianzas_precipitacion.tif"), format = "GTiff", overwrite = TRUE)
cat("varianzas_precipitacion.tif OK\n")

cat("Integrador precipitaciones finalizado. Salidas disponibles en:", output_dir, "\n")

