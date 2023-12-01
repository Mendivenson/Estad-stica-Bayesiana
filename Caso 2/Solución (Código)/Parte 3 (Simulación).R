# -------------------------------------------------------
# Caso de Estudio 2: Simulación. Estadística Bayesiana (2023-I)
# Prueba Saber 11 2022-2: Una perspectiva multinivel
# Fecha de Última Modificación: 30/Octubre/2023
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchéz
# -------------------------------------------------------

if (file.exists('Modelos') == FALSE){dir.create('Modelos')}

# En este script finalmente se realiza la simulación del muestreador de Gibbs.
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("datos.RData")
# ============================= SIMULACIÓN =====================================

source('Gibbs.R')

# install.packages('tictoc') # Para contabilizar el tiempo que le toma a los loops muy pesados generarse
# install.packages(''beepr)  # Para reproducir un sonido cada vez que un muestreador termine de correr
# ----------------------------- Modelo 1 ---------------------------------------
# --------------------------- Modelo Normal ------------------------------------
set.seed(2012)
tictoc::tic()
write.table(Gibbs1(y = puntajes$puntaje), file = 'Modelos/M1')
tictoc::toc()  

# ----------------------------- Modelo 2 ---------------------------------------
# ----------- Modelo Normal con medias específicas por departamento ------------
set.seed(2012)
tictoc::tic()
write.table(Gibbs2(medias_j =   dataDepartamentos$mean,
                   names = dataDepartamentos$CodDepartamento,
                   n_j = dataDepartamentos$n, 
                   y = puntajes$puntaje), file = 'Modelos/M2')
tictoc::toc()  

# ----------------------------- Modelo 3 ---------------------------------------
# ------ Modelo Normal con medias y varianzas específicas por departamento -----

set.seed(2012)
tictoc::tic()
write.table(Gibbs3(y = puntajes$puntaje,
           medias = dataDepartamentos$mean,
           varianzas = dataDepartamentos$var,
           n_j = dataDepartamentos$n,
           names = dataDepartamentos$CodDepartamento), file = 'Modelos/M3')
tictoc::toc() 

# ----------------------------- Modelo 4 ---------------------------------------

set.seed(2012)
tictoc::tic()
write.table(Gibbs4(puntajes$puntaje,
           mean_jk = dataMunicipios$mean, mean_k = dataDepartamentos$mean,
           n_jk = dataMunicipios$n, n_k = dataDepartamentos$n,
           NMnpioNDep = dataDepartamentos$nMunicipios, var_jk = dataMunicipios$var,
           var_k = dataDepartamentos$var, names_jk = dataMunicipios$CodMunicipio,
           names_k = dataDepartamentos$CodDepartamento),
           file = 'Modelos/M4')
tictoc::toc()
     
# ----------------------------- Modelo 5 ---------------------------------------

set.seed(2012)
tictoc::tic()
write.table(Gibbs5(puntajes$puntaje,
                   mean_jk = dataMunicipios$mean, mean_k = dataDepartamentos$mean,
                   n_jk = dataMunicipios$n, n_k = dataDepartamentos$n,
                   NMnpioNDep = dataDepartamentos$nMunicipios, var_jk = dataMunicipios$var,
                   var_k = dataDepartamentos$var, names_jk = dataMunicipios$CodMunicipio,
                   names_k = dataDepartamentos$CodDepartamento),
            file = 'Modelos/M5')
write.table(M6, file = 'Modelos/M5')
tictoc::toc()

# =============================== PUNTO 4 ======================================
# En un gráfico de cuatro paneles presente la cadena de la log-verosimilitud de 
# M2, M3, M4 y M5 (Los gráficos deben tener la misma escala para facilitar su
# comparación)
# ==============================================================================

# install.packages(data.table) # Para lectura más rápida de las tablas
library(data.table)

M1 = fread('Modelos/M1', drop = 1); M2 = fread('Modelos/M2', drop = 1);
M3 = fread('Modelos/M3', drop = 1);M4 = fread('Modelos/M4', drop = 1)
M5 = fread('Modelos/M5', drop = 1)

LLM1 = M1$LogVerosimilitud;LLM2 = M2$LogVerosimilitud
LLM3 = M3$LogVerosimilitud;LLM4 = M4$LogVerosimilitud
LLM5 = M5$LogVerosimilitud

# CADENAS DE LOG VEROSIMILITUD

par(mfrow = c(2, 2), mar = c(5, 4, 2, 1))

b = c(min(LLM2, LLM5),max(LLM2, LLM5)) # Límite en y común de los gráficos

# Plot Modelo 2
plot(LLM2, pch = 19, cex = 0.22, col = 'darkred',ylab = 'Log-verosimilitud', 
     xlab = 'Iteración',main = 'Modelo 2',ylim = b)

# Plot Modelo 3
plot(LLM3, pch = 19, cex = 0.22, col = 'darkorange',ylim = b, ylab = 'Log-verosimilitud',
     xlab = 'Iteración',main = 'Modelo 3')

# Plot Modelo 4
plot(LLM4, pch = 19, cex = 0.22, col = 'yellow4',ylim = b, ylab = 'Log-verosimilitud', 
     xlab = 'Iteración',main = 'Modelo 4')

# Plot Modelo 5
plot(LLM5, pch = 19, cex = 0.22, col = 'darkgreen',ylim = b, ylab = 'Log-verosimilitud',
     xlab = 'Iteración',main = 'Modelo 5')

# Plot Modelo 5 y 4
# par(mfrow = c(1,2))
# b = range(LLM4,LLM5)
# plot(LLM4, pch = 19, cex = 0.1, col = 'yellow4',ylim = b, ylab = 'Log-verosimilitud',
#      xlab = 'Iteración',main = 'Modelo 4')
# abline(h = mean(LLM4), lty = 'dashed')
# plot(LLM5, pch = 19, cex = 0.1, col = 'darkgreen',ylim = b, ylab = 'Log-verosimilitud',
#      xlab = 'Iteración',main = 'Modelo 5')
# abline(h = mean(LLM5), lty = 'dashed')

#  COEFICIENTE DE VARIACIÓN DE MONTECARLO

## Quitando la log verosimilitud de la matriz de los modelos
M1 = M1[,-ncol(M1),with = F]; M2 = M2[,-ncol(M2),with = F]
M3 = M3[,-ncol(M3),with = F]; M4 = M4[,-ncol(M4),with = F]
M5 = M5[,-ncol(M5),with = F]

## Calculando los tamaños efectivos de muestra:

Modelos = c('M1','M2','M3','M4','M5')

# install.packages('coda')
neff = list()
for (i in Modelos){
  neff[[i]] = coda::effectiveSize(get(i))
}

## Cálculando los errores estándar de MonteCarlo
EMC = list()
for (i in Modelos){
  EMC[[i]] = apply(X = get(i), MARGIN = 2, FUN = sd)/sqrt(neff[[i]])
}

## Finalmente, calculando los coeficientes de variación de Monte Carlo
CV = list()
for (i in Modelos){
  CV[[i]] = 100*abs(EMC[[i]]/colMeans(get(i)))
}

## Presentación de los resultados:

summaryCV = matrix(NA, nrow = 5, ncol = 6)
colnames(summaryCV) = c('Mín', '1er cuartil', 'Mediana', 'Media','3er cuartil', 'Máx')
rownames(summaryCV) = Modelos
for (i in Modelos){
  summaryCV[i,] = summary(CV[[i]])
}
rownames(summaryCV) = c('Modelo 1', 'Modelo 2', 'Modelo 3', 'Modelo 4', 'Modelo 5')
summaryCV


# install.packages('knitr')
# library(knitr)
# kable(round(summaryCV,4), format = 'latex')

# =============================== PUNTO 5 ======================================
# Calcular el DIC y el WAIC de los cinco modelos y presentar los resultados de
# de forma tabular.
# ==============================================================================

# DIC (Deviance Information Criterion)

## Modelo 1
Log1Est = sum(dnorm(puntajes$puntaje, mean = mean(M1$theta), sd = sqrt(mean(M1$sigma2)), log = T))
pDIC1 = 2 * (Log1Est - mean(LLM1))
DIC1 = -2*Log1Est + 2 *pDIC1

## Modelo 2
Log2Est = sum(dnorm(puntajes$puntaje,mean = rep(colMeans(M2[,1:32]),dataDepartamentos$n),sd = sqrt(mean(M2$sigma2)), log = T))
pDIC2 = 2 * (Log2Est - mean(LLM2))
DIC2 = -2 * Log2Est + 2 * pDIC2
  
## Modelo 3
Log3Est = sum(dnorm(puntajes$puntaje,mean = rep(colMeans(M3[,1:32]),dataDepartamentos$n),sd = rep(sqrt(colMeans(M3[,33:64])),dataDepartamentos$n), log = T))
pDIC3 = 2 * (Log3Est - mean(LLM3))
DIC3 = -2 * Log3Est + 2 * pDIC3

## Modelo 4
Log4Est = sum(dnorm(puntajes$puntaje,mean = rep(colMeans(M4[,1:1112]),dataMunicipios$n),sd = sqrt(mean(M4$kappaSq)), log = T))
pDIC4 = 2 * (Log4Est - mean(LLM4))
DIC4 = -2 * Log4Est + 2 * pDIC4
  
## Modelo 5
Log5Est = sum(dnorm(puntajes$puntaje,mean = rep(colMeans(M5[,1:1112]),dataMunicipios$n),sd = sqrt(mean(M5$kappaSq)), log = T))
pDIC5 = 2 * (Log5Est - mean(LLM5))
DIC5 = -2 * Log5Est + 2 * pDIC5

DIC = rbind('Modelo 1' = c(pDIC1,DIC1),'Modelo 2' = c(pDIC2,DIC2),
            'Modelo 3' = c(pDIC3,DIC3),'Modelo 4' = c(pDIC4,DIC4),
            'Modelo 5' = c(pDIC5,DIC5))

colnames(DIC) = c('pDIC','DIC')

# WAIC (Watanabe-Akaike Information Criterion)

## Modelo 1:
lppd1 =  0; pWAIC1 =  0
Progress = txtProgressBar(min = 1, max = nrow(puntajes), style = 3); tictoc::tic()
for (i in 1:nrow(puntajes)){
  dnorm1 = dnorm(puntajes$puntaje[i], 
                 mean = M1$theta, 
                 sd = sqrt(M1$sigma2))
  lppd1  = lppd1 + log(mean(dnorm1))
  dnorm2 = dnorm(puntajes$puntaje[i], 
                 mean = M1$theta, 
                 sd = sqrt(M1$sigma2), log = T)
  pWAIC1 = pWAIC1 + 2 * (log(mean(dnorm1)) - mean(dnorm2))
  setTxtProgressBar(Progress, i)
}
close(Progress); tictoc::toc()
WAIC1 = - 2 * lppd1 + 2 * pWAIC1

## Modelo 2:
lppd2 = 0
pWAIC2 = 0
Progress = txtProgressBar(min = 1, max = nrow(puntajes), style = 3)
tictoc::tic()
for (i in 1:nrow(puntajes)){
  dep = puntajes$CodDepartamento[i]
  dnorm1 = dnorm(puntajes$puntaje[i],
                 mean = M2[, ..dep][[1]],
                 sd = sqrt(M2$sigma2))
  lppd2  = lppd2 + log(mean(dnorm1))
  
  dnorm2 = dnorm(puntajes$puntaje[i],
                 mean = M2[, ..dep][[1]],
                 sd = sqrt(M2$sigma2),
                 log = T)
  pWAIC2 = pWAIC2 +  2 * (log(mean(dnorm1)) - mean(dnorm2))
  setTxtProgressBar(Progress,i)
}
close(Progress); tictoc::toc()
WAIC2 = - 2 * lppd2 + 2 * pWAIC2

## Modelo 3
lppd3 = 0
pWAIC3 = 0
Progress = txtProgressBar(min = 1, max = nrow(puntajes), style = 3)
tictoc:tic()
for (i in 1:nrow(puntajes)){
  dep = puntajes$CodDepartamento[i]
  thetadep = paste('theta',dep)
  sigmadep = paste('sigma',dep)
  dnorm1 = dnorm(puntajes$puntaje[i],
                 mean = M3[, ..thetadep][[1]],
                 sd = sqrt(M3[, ..sigmadep][[1]]))
  lppd3  = lppd3 + log(mean(dnorm1))
  
  dnorm2 = dnorm(puntajes$puntaje[i],
                 mean = M3[, ..thetadep][[1]],
                 sd = sqrt(M3[, ..sigmadep][[1]]),
                 log = T)
  pWAIC3 = pWAIC3 +  2 * (log(mean(dnorm1)) - mean(dnorm2))
  setTxtProgressBar(Progress,i)
}
close(Progress); tictoc::toc()
WAIC3 = - 2 * lppd3 + 2 * pWAIC3

## Modelo 4
lppd4 = 0
pWAIC4 = 0
Progress = txtProgressBar(min = 1, max = nrow(puntajes), style = 3)
tictoc::tic()
for (i in 1:nrow(puntajes)){
  mun = paste('M', puntajes$CodMunicipio[i])
  dnorm1 = dnorm(puntajes$puntaje[i],
                 mean = M4[, ..mun][[1]],
                 sd = sqrt(M4$kappaSq))
  lppd4  = lppd4 + log(mean(dnorm1))
  
  dnorm2 = dnorm(puntajes$puntaje[i],
                 mean = M4[, ..mun][[1]],
                 sd = sqrt(M4$kappaSq),
                 log = T)
  pWAIC4 = pWAIC4 +  2 * (log(mean(dnorm1)) - mean(dnorm2))
  setTxtProgressBar(Progress,i)
}
close(Progress); tictoc::toc()
WAIC4 = - 2 * lppd4 + 2 * pWAIC4

## Modelo 5
lppd5 = 0
pWAIC5 = 0
Progress = txtProgressBar(min = 1, max = nrow(puntajes), style = 3)
tictoc::tic()
for (i in 1:nrow(puntajes)){
  mun = paste('M', puntajes$CodMunicipio[i])
  dnorm1 = dnorm(puntajes$puntaje[i],
                 mean = M5[, ..mun][[1]],
                 sd = sqrt(M5$kappaSq))
  lppd5  = lppd5 + log(mean(dnorm1))
  
  dnorm2 = dnorm(puntajes$puntaje[i],
                 mean = M5[, ..mun][[1]],
                 sd = sqrt(M5$kappaSq),
                 log = T)
  pWAIC5 = pWAIC5 +  2 * (log(mean(dnorm1)) - mean(dnorm2))
  setTxtProgressBar(Progress,i)
}
close(Progress); tictoc::toc()
WAIC5 = - 2 * lppd5 + 2 * pWAIC5


WAIC = rbind('Modelo 1' = c(pWAIC1, WAIC1),
             'Modelo 2' = c(pWAIC2, WAIC2),
             'Modelo 3' = c(pWAIC3, WAIC3),
             'Modelo 4' = c(pWAIC4, WAIC4),
             'Modelo 5' = c(pWAIC5, WAIC5))
colnames(WAIC) = c('pWAIC', 'WAIC')

### Resultados
DICyWAIC = cbind(DIC,WAIC)
DICyWAIC

#                 pDIC     DIC      pWAIC    WAIC
# Modelo 1    1.986433 5636785   1.716754 5636785
# Modelo 2   32.856347 5595932  31.080767 5595930
# Modelo 3   63.423157 5594247  59.554916 5594243
# Modelo 4 1001.970860 5551814 790.618132 5551603
# Modelo 5  989.116552 5551802 780.864368 5551594

# write.table(DICyWAIC, 'DICyWAIC')
# library(knitr)
# kable(DICyWAIC, format = 'latex')


# =============================== PUNTO 6 ======================================
# Calcular la media posterior y el intervalo de credibilidad al 95% basado en 
# percentiles de \mu de cada modelo. Presentar los resultados tabularmente.
# ==============================================================================

# En el modelo 1 no hay mu propiamente dicho entonces usaremos la estimación de 
# theta pues en este caso representa lo mismo que en los demás modelos (Media global)

ICyMedia = rbind('Modelo 1' = c(quantile(M1$theta, c(0.025,0.975)), 'mu' = mean(M1$theta)),
                 'Modelo 2' = c(quantile(M2$mu, c(0.025,0.975)), 'mu' = mean(M2$mu)),
                 'Modelo 3' = c(quantile(M3$mu, c(0.025,0.975)), 'mu' = mean(M3$mu)),
                 'Modelo 4' = c(quantile(M4$mu, c(0.025,0.975)), 'mu' = mean(M4$mu)),
                 'Modelo 5' = c(quantile(M5$mu, c(0.025,0.975)), 'mu' = mean(M5$mu)))

print(ICyMedia)

# kable(ICyMedia, format = 'latex')
