library(sp)
data("meuse")

summary(meuse$zinc)

par(mfrow = c(1, 2))
hist(meuse$zinc, main="Histograma", xlab = "zinc")
boxplot(meuse$zinc, main="Boxplot")

par(mfrow = c(1, 2))
library(akima)
int_zinc <- interp(x = meuse$x, y = meuse$y, z =meuse$zinc)

filled.contour(int_zinc,
               plot.axes= {
                      axis(1)
                      axis(2)
                      contour(int_zinc, add = TRUE, lwd = 1)},
               asp = 1)


par(mfrow = c(2, 2))
persp(int_zinc$x, int_zinc$y, int_zinc$z,
      xlab ="Longitud", las = 3, font.axis = 8, cex.lab = 0.8,
      ylab ="Latitud", las = 3, font.axis = 8, cex.lab = 0.8,
      zlab ="Zinc", las = 3, font.axis = 8, cex.lab = 0.8,
      phi = 35, theta = 40, col = "lightgreen", expand = 0.5,
      ticktype = "simple")

persp(int_zinc$x, int_zinc$y, int_zinc$z,
      xlab ="Longitud", las = 3, font.axis = 8, cex.lab = 0.8,
      ylab ="Latitud", las = 3, font.axis = 8, cex.lab = 0.8,
      zlab ="Zinc", las = 3, font.axis = 8, cex.lab = 0.8,
      phi = 35, theta = 50, col = "lightgreen", expand = 0.5,
      ticktype = "simple")

persp(int_zinc$x, int_zinc$y, int_zinc$z,
      xlab ="Longitud", las = 3, font.axis = 8, cex.lab = 0.8,
      ylab ="Latitud", las = 3, font.axis = 8, cex.lab = 0.8,
      zlab ="Zinc", las = 3, font.axis = 8, cex.lab = 0.8,
      phi = 35, theta = 60, col = "lightgreen", expand = 0.5,
      ticktype = "simple")

persp(int_zinc$x, int_zinc$y, int_zinc$z,
      xlab ="Longitud", las = 3, font.axis = 8, cex.lab = 0.8,
      ylab ="Latitud", las = 3, font.axis = 8, cex.lab = 0.8,
      zlab ="Zinc", las = 3, font.axis = 8, cex.lab = 0.8,
      phi = 35, theta = 70, col = "lightgreen", expand = 0.5,
      ticktype = "simple")



par(mfcol=c(1,2))
hist(meuse$zinc, prob = TRUE,
     main="Histograma Zinc",
     ylab= "Frecuencia",
     xlab = "zinc")
x <- seq(min(meuse$zinc), max(meuse$zinc), length = 40)
f <- dnorm(x, mean = mean(meuse$zin), sd=sd(meuse$zinc))
lines(x, f, col="red", lwd=2)

qqnorm(meuse$zinc)
qqline(meuse$zinc)


shapiro.test(meuse$zinc)

library(car)
bc_trans <- powerTransform(meuse$zinc)
bc_trans$roundlam

meuse$zinc_t <- bcPower(meuse$zinc, lambda = bc_trans$roundlam)

par(mfcol = c(1,2))

hist(meuse$zinc_t, prob = TRUE,
     main = "Histograma con curva normal",
     ylab = "Frecuencia",
     xlab = "Zinc"
     )
x <- seq(min(meuse$zinc_t), max(meuse$zinc_t), length = 40)
f <- dnorm(x, mean = mean(meuse$zinc_t), sd=sd(meuse$zinc_t))
lines(x, f, col="red", lwd=2)

qqnorm(meuse$zinc_t)
qqline(meuse$zinc_t)

shapiro.test(meuse$zinc_t)



meuse$reci_zinc <- 1/(meuse$zinc)
summary(meuse$reci_zinc)

par(mfcol = c(1,2))

hist(meuse$reci_zinc, prob = TRUE,
     main = "Histograma reci_zinc",
     ylab = "Frecuencia",
     xlab = "Zinc"
)
x <- seq(min(meuse$reci_zinc), max(meuse$reci_zinc), length = 40)
f <- dnorm(x, mean = mean(meuse$reci_zinc), sd=sd(meuse$reci_zinc))
lines(x, f, col="red", lwd=2)

qqnorm(meuse$reci_zinc)
qqline(meuse$reci_zinc)

shapiro.test(meuse$reci_zinc)


par(mfcol = c(1,1))
install.packages("moments")
library(moments)
skewness(meuse$zinc)
skewness(meuse$zinc_t)
skewness(meuse$reci_zinc)


library(geoR)
zinc_geo <- as.geodata(meuse, coords.col = c("x","y"), data.col = "zinc_t")
duplicated(zinc_geo)

plot(zinc_geo, lowess= TRUE)

library(sf)

coordinates(meuse) = c("x", "y")

library(gstat)

variograma_zinc <- variogram(zinc_t ~ 1, meuse)
plot(variograma_zinc, lwd = 2)

variograma_zinc_2 <- variogram(zinc_t ~ 1, meuse, cutoff = 1550, width = 100)
plot(variograma_zinc_2, lwd = 2)

variograma_zinc$np
variograma_zinc_2$np

mod_sph <- fit.variogram(variograma_zinc_2, vgm("Sph"))
mod_exp <- fit.variogram(variograma_zinc_2, vgm("Exp"))
mod_gau <- fit.variogram(variograma_zinc_2, vgm("Gau"))


vg_line <- 
  rbind(cbind(variogramLine(mod_sph, maxdist = max(variograma_zinc_2$dist)),
                       id = "Esférico"),
                 cbind(variogramLine(mod_sph, maxdist = max(variograma_zinc_2$dist)),
                       id= "Exponencial"),
                 cbind(variogramLine(mod_sph, maxdist = max(variograma_zinc_2$dist)),
                       id = "Gaussiano"))

library(ggplot2)

graf <- 
  ggplot(variograma_zinc_2, aes(x = dist, y = gamma, color = id)) +
        geom_line(data = vg_line) +
        geom_point() +
        labs(
                title = "Semivariograma experimental y teóricas ajustados") +
                xlab("Distancia") +
                ylab("Semivarianza") +
                scale_color_discrete(name = "Modelo",
                               labels = c("Esférico", "Exponencial", "Gaussiano", "Experimental"))

graf

comp_mod <- data.frame(Modelo = c("Esférico", "Exponencial", "Gaussiano"))