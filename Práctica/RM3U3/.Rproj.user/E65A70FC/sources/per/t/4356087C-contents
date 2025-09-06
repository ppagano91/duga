library(readxl)

suelo <- read_excel("data/suelo.xlsx")

library(geoR)

suelo1 <- as.geodata(suelo, coords.col =c("longitud","latitud"), data.col = "pH")
dup.coords(suelo1)

variograma_suelo <- variog(suelo1)
plot(variograma_suelo, main="Varigorama: valores de pH")

# windows()
# eyefit(variograma_suelo)

# dev.off()

plot(variograma_suelo, pch = 20, col = "red")
lines.variomodel(cov.model="wave", cov.pars = c(1.1,0.09), nug = 3.8, col = "blue", lwd = 2, lty = 1)

variograma_ml <- likfit(suelo1, cov.model = "wave", ini.cov.pars = c(1.1, 0.09), lik.method = "ML", nugget = 3.8)
variograma_ml
summary(variograma_ml)

variograma_reml <- likfit(suelo1, cov.model = "wave", ini.cov.pars = c(1.1, 0.09), lik.method = "REML", nugget = 3.8)
variograma_reml
summary(variograma_reml)

plot(variograma_suelo)
lines(variograma_ml, col = "red", lwd = 2)
lines(variograma_reml, col = "green", lwd = 2)

legend("topleft", legend = c("ML", "REML"), lty=c(1,1), col=c("red", "green"), bty="n")

# Grilla de Predicción
pred_grilla <- expand.grid(x = seq(min(suelo1$coords[,1]),
                                   max(suelo1$coords[,1]), l= 100),
                           y = seq(min(suelo1$coords[,2]),
                                   max(suelo1$coords[,2]), l= 100)
                           )

plot(suelo1$coords, pch = 20, asp = 1)

points(pred_grilla, pch = "+", cex = 0.2, col="green")


ko_reml <- krige.conv(suelo1, locations = pred_grilla, krige = krige.control(obj.model = variograma_reml))

summary(ko_reml)


# install.packages("sp")
library(sp)

suelo_sp <- suelo
coordinates(suelo_sp) <- ~ longitud + latitud

ko_reml_sp <- SpatialPixelsDataFrame(points=pred_grilla, data = data.frame(ko_reml[1:2]))

spplot(ko_reml_sp, zcol="predict", col.regions=heat.colors(100), main="Predicciones de pH", xlab="X Coord", ylab="Y Coord")

spplot(ko_reml_sp, zcol="predict", col.regions=heat.colors(100), contour=TRUE, main="Predicciones de pH", xlab="X Coord", ylab="Y Coord")


spplot(ko_reml_sp, zcol="predict",
       sp.layout = list("sp.points", suelo_sp, col="black", pch=20),
       col.regions=heat.colors(100),
       contour=TRUE,
       main="Predicciones de pH",
       xlab="X Coord",
       ylab="Y Coord")

ko_reml_var_sp <- SpatialPixelsDataFrame(points = pred_grilla, data = data.frame(ko_reml$krige.var))

spplot(ko_reml_var_sp, zcol="ko_reml.krige.var",
       sp.layout = list("sp.points", suelo_sp, col="black", pch=20),
       col.regions=rev(heat.colors(100)),
       contour=FALSE,
       main="Superficie de varianzas y puntos de muestreo (pH)",
       xlab="X Coord", ylab="Y Coord",
       key.space = "right"
       )

x_ticks <- round(seq(min(suelo1$coords[,1]), max(suelo1$coords[,1]), by = 0.5), digits = 2)
y_ticks <- round(seq(min(suelo1$coords[,2]), max(suelo1$coords[,2]), by = 0.5), digits = 2)

spplot(ko_reml_sp, zcol="predict",
       sp.layout = list("sp.points", suelo_sp, col="black", pch=20),
       col.regions=heat.colors(100),
       contour=FALSE,
       main="Predicciones de pH",
       xlab="Longitud", ylab="Latitud",
       scales = list(x = list(at = x_ticks),
                     y = list(at = y_ticks)
                     )
)

ko_reml_var_sp <- SpatialPixelsDataFrame(points = pred_grilla, data = data.frame(ko_reml$krige.var))

spplot(ko_reml_var_sp, zcol="ko_reml.krige.var",
       sp.layout = list("sp.points", suelo_sp, col="black", pch=20),
       col.regions=rev(heat.colors(100)),
       contour=FALSE,
       main="Superficie de varianzas y puntos de muestreo (pH)",
       xlab="Longitud", ylab="Latitud",
       key.space = "right",
       scales = list(x = list(at = x_ticks),
                     y = list(at = y_ticks)
       )
)

# Exportación
library(sf)
library(raster)

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

st_write(bbox_sf_grid, "area_prediccion_poligono.shp", driver="ESRI Shapefile",  append=FALSE)
cat("area_prediccion_poligono.shp OK")

puntos_muestreo_sf <- st_as_sf(suelo_sp)
if (is.na(st_crs(puntos_muestreo_sf))){
  st_crs(puntos_muestreo_sf) <- 4326
}

st_write(puntos_muestreo_sf, "puntos_muestreo.shp", driver="ESRI Shapefile",  append=FALSE)
cat("puntos_muestreo.shp OK")

x_grid_unique <- unique(pred_grilla$x)
y_grid_unique <- unique(pred_grilla$y)
z_matrix_predict <- matrix(ko_reml$predict, nrow = length(x_grid_unique), ncol=length(y_grid_unique))

contour_lines_predict_data <-  contourLines(x=x_grid_unique, y=y_grid_unique, z=z_matrix_predict)

lines_sf_predict_list <- list()
levels_predict_list <- numeric()

for (i in seq_along(contour_lines_predict_data)) {
  coords <- cbind(contour_lines_predict_data[[i]]$x, contour_lines_predict_data[[i]]$y)
  line_sf <- st_linestring(coords)
  lines_sf_predict_list[[i]] <- line_sf
  levels_predict_list <- c(levels_predict_list, contour_lines_predict_data[[i]]$level)
}

all_lines_sfc_predict <- st_sfc(lines_sf_predict_list, crs=4326)

contours_predict_sf <- st_sf(level = levels_predict_list, geometry = all_lines_sfc_predict)

st_write(contours_predict_sf, "contornoes_prediccion.shp", driver= "ESRI Shapefile", append=FALSE)
cat("contornoes_prediccion.shp OK")


raster_predicciones_sp <- raster(ko_reml_sp)

crs(raster_predicciones_sp) <-CRS("+init=epsg:4326")

writeRaster(raster_predicciones_sp, "predicciones_pH.tif", format="GTiff", overwrite=TRUE)
cat("predicciones_pH.tif OK")


raster_varianzas_sp <- raster(ko_reml_var_sp)

crs(raster_varianzas_sp) <-CRS("+init=epsg:4326")

writeRaster(raster_varianzas_sp, "varianzas_pH.tif", format="GTiff", overwrite=TRUE)
cat("varianzas_pH.tif OK")