library(readxl)

suelo <- read_excel("data/suelo.xlsx")

library(geoR)

suelo1 <- as.geodata(suelo, coords.col =c("longitud","latitud"), data.col = "pH")
dup.coords(suelo1)

variograma_suelo <- variog(suelo1)
plot(variograma_suelo, main="Varigorama: valores de pH")

windows()
eyefit(variograma_suelo)

dev.off()

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
legend("topleft", legend = c("ML", "REML"), lty=c(1,1), col=("red", "green"), bty = "n")

