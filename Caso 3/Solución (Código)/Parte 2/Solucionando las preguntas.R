# -------------------------------------------------------
# Caso de Estudio 3: Punto 2. Estadística Bayesiana (2023-I)
# Datos de diabetes: Desarrollo de las preguntas
# Fecha de Última Modificación: 26/Noviembre/2023
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchez
# -------------------------------------------------------

# ================================ BASE DE DATOS ===============================
testData  = source('http://www2.stat.duke.edu/~pdh10/FCBS/Inline/yX.diabetes.test')$value
str(testData)

# ---> Dividimos los datos entre variables dependientes e independientes:
y = testData[,1]
X = testData[,-1]

# ============================== PUNTO 1 =======================================
# ================= FITTED VS PREDICTED Y ERROR MEDIO ABSOLUTO =================

# ---> Cargando los modelos:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
Modelo1 = colMeans(read.table('Modelos/Modelo1'))
Modelo2 = colMeans(read.table('Modelos/Modelo2'))
Modelo3 = colMeans(read.table('Modelos/Modelo3'))
Modelo4 = colMeans(read.table('Modelos/Modelo4'))
# Nota: Para el primer punto sólo necesitamos las medias de los modelos.

# ---> Algunas constantes necesarias para el desarrollo de este punto:
n =  length(y)                                            # Número de datos de test
p = length(Modelo1) - 1                                   # Número de betas estimados

# ---> CALCULANDO LAS ESTIMACIONES DE Y:
Estimaciones1 = X %*% Modelo1[1:p]
Estimaciones2 = X %*% Modelo2[1:p]
Estimaciones3 = X %*% Modelo3[1:p]
Estimaciones4 = X %*% Modelo4[1:p]

# ---> GRAFICANDO VALORES PREDICHOS VS. ESTIMADOS:
yl = range(y)                                             # Las gráficas tienen el mismo rango en  
xl = range(Estimaciones1,Estimaciones2,                   # el eje y y en el eje x.
           Estimaciones3,Estimaciones4)

# Para dibujar las 4 gráficas simultáneamente:
par(oma = c(2,1,2,1),
    mfrow = c(2,2))

# MODELO 1
plot(x = Estimaciones1, y = y, xlab = expression(hat(y)[italic("test")]),
     ylab = expression(y[italic("test")]),
     main = paste('Modelo de previa unitaria\nEAM = ',round(sum(abs(Estimaciones1 - y)) * 1/n,3)),
     xlim = xl, ylim = yl, cex.main = 1.1,
     pch = 1, cex = 0.5)
abline(a = 0, b = 1, col = 'blue', lty = 1) 

# MODELO 2
plot(x = Estimaciones2, y = y, xlab = expression(hat(y)[italic("test")]),
     ylab = expression(y[italic("test")]),
     main = paste('Modelo de previa g\nEAM = ',round(sum(abs(Estimaciones2 - y)) * 1/n,3)), 
     xlim = xl, ylim = yl, cex.main = 1.1,
     pch = 1, cex =0.5)
abline(a = 0, b = 1, col = 'blue', lty = 1) 

# MODELO 3
plot(x = Estimaciones3, y = y, xlab = expression(hat(y)[italic("test")]),
     ylab = expression(y[italic("test")]),
     main = paste('Regresión rígida\nEAM = ', round(sum(abs(Estimaciones3 - y)) * 1/n,3)),
     xlim = xl, ylim = yl, cex.main = 1.1,
     pch = 1, cex =0.5)
abline(a = 0, b = 1, col = 'blue', lty = 1) 

# MODELO 4
plot(x = Estimaciones4, y = y, xlab = expression(hat(y)[italic("train")]),
     ylab = expression(y[italic("train")]),
     main = paste('Regresión con errores\ncorrelacionados\nEAM = ', round(sum(abs(Estimaciones4 - y)) * 1/n,3)), 
     xlim = xl, ylim = yl, cex.main = 1.1,
     pch = 1, cex =0.5)
abline(a = 0, b = 1, col = 'blue', lty = 1) 

mtext('Valores reales vs. valores predichos', side = 3, outer = TRUE, cex = 1.5, font = 2)

# ============================== PUNTO 2 =======================================
# ================ POSTERIOR DE LA MEDIA Y VALORES PPP =========================

rm(list = ls())

# ---> Cargando los datos de entrenamiento para revisar la bondad de ajuste:
trainData = source('http://www2.stat.duke.edu/~pdh10/FCBS/Inline/yX.diabetes.train')$value

# Separando las variables dependientes de las independientes:
y = trainData[,1]
yMean = mean(y)                                   # De y sólo necesitamos la media.
XTrain = as.matrix(trainData[, -1])
n = length(y)

# ---> Cargando los modelos:
Modelo1 = as.matrix(read.table('Modelos/Modelo1'))
Modelo2 = as.matrix(read.table('Modelos/Modelo2'))
Modelo3 = as.matrix(read.table('Modelos/Modelo3'))
Modelo4 = as.matrix(read.table('Modelos/Modelo4'))

# Para poder usar lapply:
Modelo1 = split(Modelo1, seq_len(nrow(Modelo1)))
Modelo2 = split(Modelo2, seq_len(nrow(Modelo2)))
Modelo3 = split(Modelo3, seq_len(nrow(Modelo3)))
Modelo4 = split(Modelo4, seq_len(nrow(Modelo4)))

# ---> CALCULANDO LA DISTRIBUCIÓN POSTERIOR DE LA MEDIA DE Y:

# MODELO 1:
set.seed(2012)
MediasModelo1 = unlist(lapply(Modelo1, 
                                FUN = function(x) 
                                  mean(rnorm(n = n, mean = XTrain %*% x[1:64],
                                             sd = sqrt(x[65])))))

# MODELO 2:
set.seed(2012)
MediasModelo2 = unlist(lapply(Modelo2, 
                                FUN = function(x) 
                                  mean(rnorm(n = n, mean = XTrain %*% x[1:64],
                                             sd = sqrt(x[65])))))

# MODELO 3:
set.seed(2012)
MediasModelo3 = unlist(lapply(Modelo3, 
                                FUN = function(x) 
                                  mean(rnorm(n = n, mean = XTrain %*% x[1:64],
                                             sd = sqrt(x[65])))))

# MODELO 4:
set.seed(2012)
MediasModelo4 = unlist(lapply(Modelo4, 
                                FUN = function(x) 
                                  mean(rnorm(n = n, mean = XTrain %*% x[1:64],
                                             sd = sqrt(x[65])))))

# ---> GRAFICANDO  LAS DISTRIBUCIONES POSTERIORES Y CALCULADO LOS
#      VALORES PPP:
par(oma = c(2,1,2,1),
    mfrow = c(2,2))

# ---> MODELO 1: 
hist(MediasModelo1, col = 'gray90', border = 'gray90', 
     main = paste0('Modelo 1 (Previa unitaria)\n(Valor ppp = ', round(mean(MediasModelo1 < yMean),3),')'),
     ylab = 'Densidad', xlab = 'Media posterior de y ',
     cex.main = 1.1,prob = T)
# ------> DENSIDAD Y VALOR REAL DE LA MEDIA DE Y:
lines(density(MediasModelo1, kernel = 'biweight'),col = 'black')
abline(v = yMean, col = 'red', lty = 2)

# ---> MODELO 2:
hist(MediasModelo2, col = 'gray90', border = 'gray90', 
     main = paste0('Modelo 2 (Previa g)\n(Valor ppp = ', round(mean(MediasModelo2 < yMean),3),')'),
     ylab = 'Densidad', xlab = 'Media posterior de y ', cex.main = 1.1,
     prob = T)
# ------> DENSIDAD Y VALOR REAL DE LA MEDIA DE Y:
lines(density(MediasModelo2, kernel = 'biweight'),col = 'black')
abline(v = yMean, col = 'red', lty = 2)

# ---> MODELO 3:
hist(MediasModelo3, col = 'gray90', border = 'gray90', 
     main = paste0('Modelo 3 (Regresión rígida)\n(Valor ppp = ', round(mean(MediasModelo3 < yMean),3),')'),
     ylab = 'Densidad', xlab = 'Media posterior de y ', cex.main = 1.1,
     prob = T)
# ------> DENSIDAD Y VALOR REAL DE LA MEDIA DE Y:
lines(density(MediasModelo3, kernel = 'biweight'),col = 'black')
abline(v = yMean, col = 'red', lty = 2)

# ---> MODELO 4:
hist(MediasModelo4, col = 'gray90', border = 'gray90', 
     main = paste0('Modelo 4 (Errores correlacionados)\n(Valor ppp = ', round(mean(MediasModelo4 < yMean),3),')'),
     ylab = 'Densidad', xlab = 'Media posterior de y ', cex.main = 1.1,
     prob = T)
# ------> DENSIDAD Y VALOR REAL DE LA MEDIA DE Y:
lines(density(MediasModelo4, kernel = 'biweight'),col = 'black')
abline(v = yMean, col = 'red', lty = 2)

# Título general del gráfico:
mtext('Bondad de ajuste', side = 3, outer = TRUE, cex = 1.5, font = 2)

# Leyenda general del gráfico:
par(fig = c(0, 1, 0, 1), oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = 'Media real de y',
       col = 'red' , xpd = TRUE, lty = 2,
       horiz = FALSE, cex = 1, seg.len=1, bty = 'n')

# Corra la siguiente línea para reestablecer los parámetros de los gráficos en R:
graphics.off()


# ============================== PUNTO 3 =======================================
# ========================== CÁLCULO DEL DIC ===================================

# Recuerde que DIC =  -2 log p(y | \hat\theta_{bayes}) + 2 p_{DIC}
# Con p_{DIC} = 2 (log p(y | \hat\theta_{bayes}) - 1/B sum log p(y | theta^{b}))
rm(list = ls())

# ---> Cargando los datos de entrenamiento para calcular el DIC:
trainData = source('http://www2.stat.duke.edu/~pdh10/FCBS/Inline/yX.diabetes.train')$value

# Separando las variables dependientes de las independientes:
y = trainData[,1]
XTrain = as.matrix(trainData[, -1])

# ---> Cargando los modelos:
Modelo1 = as.matrix(read.table('Modelos/Modelo1'))
hatModelo1 = colMeans(Modelo1)
Modelo2 = as.matrix(read.table('Modelos/Modelo2'))
hatModelo2 = colMeans(Modelo2)
Modelo3 = as.matrix(read.table('Modelos/Modelo3'))
hatModelo3 = colMeans(Modelo3)
Modelo4 = as.matrix(read.table('Modelos/Modelo4'))
hatModelo4 = colMeans(Modelo4)

# RESUMEN COEFICIENTES DE VARIACIÓN DE MONTECARLO PARA REVISAR EL TAMAÑO DE MUESTRA:
CV = matrix(ncol = 4, nrow = ncol(Modelo4)); colnames(CV) = c('Modelo 1', 'Modelo 2', 'Modelo 3', 'Modelo 4')
CV[,1] = c(100* abs((apply(Modelo1,MARGIN = 2, FUN = sd)/(sqrt(coda::effectiveSize(Modelo1)))) / hatModelo1), NA)
CV[,2] = c(100* abs((apply(Modelo2,MARGIN = 2, FUN = sd)/(sqrt(coda::effectiveSize(Modelo2)))) / hatModelo2), NA)
CV[,3] = 100* abs((apply(Modelo3,MARGIN = 2, FUN = sd)/(sqrt(coda::effectiveSize(Modelo3)))) / hatModelo3)
CV[,4] = 100* abs((apply(Modelo4,MARGIN = 2, FUN = sd)/(sqrt(coda::effectiveSize(Modelo4)))) / hatModelo4)
boxplot(CV[,1:3], ylab = 'Porcentaje', main = 'Coeficiente de variación de Monte Carlo')
boxplot(CV[,4], ylab = 'Porcentaje', main = 'CV\nMonte Carlo', xlab = 'Modelo 4')

CoefVar = apply(CV, MARGIN = 2, function(x) c(min(x, na.rm = T), 
                                    quantile(x, 0.25, na.rm = T),
                                    mean(x, na.rm = T),
                                    quantile(x, 0.75, na.rm = T),
                                    max(x,na.rm = T)))
rownames(CoefVar)[1] = 'min'
rownames(CoefVar)[3] = 'mean'
rownames(CoefVar)[5] = 'max'

# library(knitr)
# kable(round(CoefVar,3), format = 'latex')

Modelo1 = split(Modelo1, seq_len(nrow(Modelo1)))
Modelo2 = split(Modelo2, seq_len(nrow(Modelo2)))
Modelo3 = split(Modelo3, seq_len(nrow(Modelo3)))
Modelo4 = split(Modelo4, seq_len(nrow(Modelo4)))


# DIC MODELO 1:
Log1Est = sum(dnorm(y, mean = XTrain %*% hatModelo1[1:64],sd = sqrt(hatModelo1[65]), log = T))
LLM1 = unlist(lapply(Modelo1, FUN = function(x) sum(dnorm(y, mean = XTrain %*% x[1:64],
                                                          sd = sqrt(x[65]),
                                                          log = T))))
pDIC1 = 2 * (Log1Est - mean(LLM1))
DIC1 = -2*Log1Est + 2 * pDIC1

# DIC MODELO 2:
Log2Est = sum(dnorm(y, mean = XTrain %*% hatModelo2[1:64],sd = sqrt(hatModelo2[65]), log = T))
LLM2 = unlist(lapply(Modelo2, FUN = function(x) sum(dnorm(y,
                                                          mean = XTrain %*% x[1:64],
                                                          sd = sqrt(x[65]),
                                                          log = T))))
pDIC2 = 2 * (Log2Est - mean(LLM2))
DIC2 = - 2 * Log2Est + 2 * pDIC2

# DIC MODELO 3:
Log3Est = sum(dnorm(y, mean = XTrain %*% hatModelo3[1:64],sd = sqrt(hatModelo3[65]), log = T))
LLM3 = unlist(lapply(Modelo3, FUN = function(x) sum(dnorm(y,
                                                          mean = XTrain %*% x[1:64],
                                                          sd = sqrt(x[65]),
                                                          log = T))))
pDIC3 = 2 * (Log3Est - mean(LLM3))
DIC3 = - 2 * Log3Est + 2 * pDIC3

# DIC MODELO 4:
Log4Est = sum(dnorm(y, mean = XTrain %*% hatModelo4[1:64],sd = sqrt(hatModelo4[65]), log = T))
LLM4 = unlist(lapply(Modelo4, FUN = function(x) sum(dnorm(y,
                                                          mean = XTrain %*% x[1:64],
                                                          sd = sqrt(x[65]),
                                                          log = T))))
pDIC4 = 2 * (Log4Est - mean(LLM4))
DIC4 = - 2 * Log4Est + 2 * pDIC4

DICS = matrix(ncol = 4, nrow = 2)
colnames(DICS) = c('Previa unitaria', 'Previa g', 'Regresión rígida', 'Regresión con errores correlacionados')
rownames(DICS) = c('DIC', 'pDIC')
DICS[1,] = c(DIC1, DIC2, DIC3, DIC4)
DICS[2,] = c(pDIC1, pDIC2, pDIC3, pDIC4) 
DICS
# library(knitr)
# kable(t(round(DICS,3)), format = 'latex')

# PLOT DE LAS LOG VEROSIMILITUDES PARA REVISAR LA ESTABILIDAD DE LOS MODELOS:
par(mfrow = c(1,3))
yRange = range(LLM1,LLM3, LLM4)
plot(LLM1, col = 'darkgreen', lty = 1, type = 'l', ylim = yRange, xlab = 'Iteración',
     ylab = 'Log-verpsimilitud', main = 'Modelo 1 (Previa unitaria)')
plot(LLM3, col = 'darkred', lty = 1, type = 'l', ylim = yRange, xlab = 'Iteración',
     ylab = 'Log-verpsimilitud', main = 'Modelo 3 (Regresión rígida)')
plot(LLM4, col = 'darkblue', lty = 1, type = 'l', ylim = yRange, xlab = 'Iteración',
     ylab = 'Log-verpsimilitud', main = 'Modelo 4 (Errores correlacionados)')


