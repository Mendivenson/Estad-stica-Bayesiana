# -------------------------------------------------------
# Caso de Estudio 1: Parte 1. Estadística Bayesiana (2023-I)
# Análisis de Tiempo de Reparación (Verizon)
# Fecha de Última Modificación: 2023-09-11
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchéz
# -------------------------------------------------------

# En este script se encuentran las líneas de código necesarias para ajustar un
# modelo bayesiano a los datos contenidos en Verizon.csv y dar solución a la 
# Parte 1 del caso de estudio presentado en clase. Los datos se separan en dos 
# poblaciones: ILEC y CLEC de las cuales se habla más en el informe entregado 
# junto con este script.


# ---------------- TRATAMIENTO PREVIO DE LA INFORMACIÓN ------------------------

rm(list = ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv('Verizon.csv', header = T, sep = ',')
dataILEC <- data[data$Group == 'ILEC',][['Time']]
dataCLEC <- data[data$Group == 'CLEC',][['Time']]
rm(data)

# Se usa la semilla 20122002 para que el código sea reproducible en todos los 
# lugares necesarios del siguiente script.

# ==================== PUNTO 1. DISTRIBUCIÓN POSTERIOR DE ETA ==================

# ------------------- Ajuste del modelo para Eta -------------------------------

# Función para hacer simulaciones de Montecarlo de la distribución posterior (Gamma Inversa)
# ajustada a unos datos que se suponen provenientes de una exponencial. (Modelo Gamma Inversa- Exponencial)

AjusteMCGIExp <- function(x,a,b,n = 10000){
  AlphaNew = a + length(x)
  BetaNew = b + sum(x)
  resultadoMC = 1/rgamma(n, shape = AlphaNew, rate = BetaNew)
  return(resultadoMC)
}


AlphaBefore = 3
BetaBefore = 17

# Obteniendo muestras de la media la población ILEC (Clientes de la empresa)
set.seed(12202002)
GenDataILEC = AjusteMCGIExp(dataILEC,AlphaBefore,BetaBefore, n = 50000)

# Obteniendo muestras de la media CLEC (Clientes por fuera de la empresa)
GenDataCLEC = AjusteMCGIExp(dataCLEC,AlphaBefore,BetaBefore, n = 50000)

# Finalmente las muestras de eta pueden ser obtenidas como:
GenDataEta = GenDataILEC - GenDataCLEC


# ----------------------------- Reporte Gráfico --------------------------------

# install.packages('latex2exp') # Paquete para etiquetas con ecuaciones LaTeX.
library(latex2exp)

hist(GenDataEta, main = TeX('Distribución posterior de $\\eta$'),
     ylab = TeX('$p(\\eta | \\, y)'),
     xlab = TeX('$\\eta$'),
     col = 'gray90',
     border = 'gray90',
     prob = T)
abline(v = mean(GenDataEta),col = "red", lty = 2)
abline(v = c(quantile(GenDataEta, 0.025),quantile(GenDataEta, 0.975)),
       col = "blue", lty = 2)
abline(v = c(quantile(GenDataEta, 0.01), quantile(GenDataEta, 0.99)),
       col = "green", lty = 2)
lines(density(GenDataEta, kernel = 'biweight'), col = 'darkgreen',lwd = 1)
legend("topleft", border = NULL,
       legend =  c('Media',
                   'IC 95%', 
                   'IC 98%',
                   'Distribución Posterior'),
       col = c('red','blue','green','darkgreen'),
       lwd = 1,lty = c(2,2,2,1), bty = 'n')

# ----------------------------- Reporte Tabular --------------------------------

InfoParte1 = matrix(0, nrow = 1,ncol=3)
colnames(InfoParte1) = c('Media', 
                         'Intervalo de credibilidad (95%)', 
                         'Coef. de variación') 
InfoParte1[1,'Media'] = round(mean(GenDataEta),3)
InfoParte1[1,'Intervalo de credibilidad (95%)'] = paste0('(',round(quantile(GenDataEta, 0.025),2),', '
         , round(quantile(GenDataEta, 0.975),2),')')
InfoParte1[1,'Coef. de variación'] = 
  paste(100*round(sd(GenDataEta)/abs(mean(GenDataEta)),3),'%')

print(InfoParte1)
 
print(paste0('La probabilidad de que eta sea mayor que cero ', 
             round(mean(GenDataEta < 0) * 100,2)))

# ----------------- ADICIONAL: Distribuciones posteriores ----------------------

# Un gráfico extra que pudiera ser de intéres sería el de las distribuciones posteriores
# del modelo gamma inversa exponencial ajustado con los datos de la población CLEC 
# y de la población ILEC por separado.

# install.packages('invgamma') # Una librería que agrega la distribución invgamma a las distribuciones disponibles en R.
library(invgamma)
curve(expr = dinvgamma(x,
                       shape = AlphaBefore + length(dataILEC), 
                       rate = BetaBefore + sum(dataILEC)),
      from = 1, to = 25,col = 'blue',
      n = 100000,
      xlab = expression(lambda),
      ylab = 'Densidad',
      main = 'Distribuciones posteriores')

curve(expr = dinvgamma(x,shape = AlphaBefore, rate = BetaBefore),
      n = 100000, add = T)
curve(expr = dinvgamma(x,
                       shape = AlphaBefore + length(dataCLEC), 
                       rate = BetaBefore + sum(dataCLEC)), col = 'green',
      add = T, n = 100000)
legend('topright',
       legend = c('Distribución previa', TeX('Distribución posterior ILEC $(\\lambda_1)$'),
                  TeX('Distribución posterior CLEC $(\\lambda_2)$')),
       col = c('black', 'blue', 'green'),
       lty = 1, cex = 0.6,bty = 'n')


# ================= PUNTO 2. ANÁLISIS DE SENSITIVIDAD PARA ETA =================

# Con la función ya creada debería ser relativamente sencillo hacer nuestro análisis de
# sensitividad:

# Valores dados para el análisis mediante comparación:
ListaDeValores = matrix(c(c(3,17,17),
                          c(2,8.5,8.5),
                          c(3,16.8,33),
                          c(2,8.4,16.5)),
                        byrow=T, nrow = 4)

# Una matriz para guardar la información que se presentará de forma tabular.
InfoParte1 <- matrix(0,ncol = 4,nrow = 6)
rownames(InfoParte1) = c('Media', 'Límite Inf', 'Límite Sup', 
                         'Coef. Variación', 'Media a priori', 'Coef. Variación a priori')

set.seed(20122002)
# Ahora hacemos el muestreo de Eta para cada configuración de distribución previa.
for (j in 1:length(ListaDeValores[,1])){
  AlphaBefore = ListaDeValores[j,1]
  BetaBefore1 = ListaDeValores[j,2]
  BetaBefore2 = ListaDeValores[j,3]
  GenDataILEC = AjusteMCGIExp(dataILEC,AlphaBefore,BetaBefore1, n = 10000)
  GenDataCLEC = AjusteMCGIExp(dataCLEC,AlphaBefore,BetaBefore2, n = 10000)
  GenDataEta = GenDataILEC - GenDataCLEC
  
  # Almacenamos la información solicitada.
  assign(paste('Distribución previa',j),GenDataEta)
  InfoParte1['Media',j] = round(mean(GenDataEta),3)
  InfoParte1['Límite Inf',j] = round(quantile(GenDataEta, 0.025),3)
  InfoParte1['Límite Sup',j] = round(quantile(GenDataEta,0.975),3)
  InfoParte1['Coef. Variación',j] = 100 * round(sd(GenDataEta)/
                                            abs(mean(GenDataEta)),3)
}

# Cálculo de los coeficientes y medias a priori para Eta como la resta de gammas inversas.
for (j in 1:length(ListaDeValores[,1])){
  AlphaBefore = ListaDeValores[j,1]
  BetaBefore1 = ListaDeValores[j,2]
  BetaBefore2 = ListaDeValores[j,3]
  
  Media1 = BetaBefore1/(AlphaBefore  - 1)
  Media2 = BetaBefore2/(AlphaBefore  - 1)
  DenVar =  ((AlphaBefore - 1)^2) * (AlphaBefore-2)
  Var1 = (BetaBefore1^2)/DenVar
  Var2 = (BetaBefore2^2)/DenVar

  InfoParte1['Media a priori', j] = 
    round(Media1 - Media2,3)
  InfoParte1['Coef. Variación a priori', j] = 
    round(100*sqrt(Var1+Var2)/ 
            abs(Media1 - Media2),3)
}

# --------------------- Estimación CV a priori por Monte Carlo -----------------

set.seed(20122002)

# Un plot vacío para poner la densidad previa de eta con los distintos hiperparámetros.
plot(NULL, ylim = c(0, 30), ylab = "Densidad",
     xlim = c(-2,2), xlab = TeX("$\\eta$"), 
     main = TeX("Distribuciones posteriores de $\\eta$"))
mtext(text = TeX('\\textbf{usando Monte Carlo con los hiperparámetros propuestos}'),
      side = 3, line = 0.5, cex = 1, font = 2)

colores = c('orange', 'blue', 'black','green')

for (j in 1:length(ListaDeValores[,1])){
  AlphaBefore = ListaDeValores[j,1]
  BetaBefore1 = ListaDeValores[j,2]
  BetaBefore2 = ListaDeValores[j,3]
  
  GenLambda1 = 1/rgamma(n = 50000, shape = AlphaBefore, rate = 1/BetaBefore1)
  GenLambda2 = 1/rgamma(n = 50000, shape = AlphaBefore, rate = 1/BetaBefore2)
  GenEta = GenLambda1 - GenLambda2
  lines(density(GenEta),col = colores[j])
  
  InfoParte1['Coef. Variación a priori', j] = 100 * sd(GenEta)/abs(mean(GenEta))
}

legend("topleft",legend = c('Distribución previa 1','Distribución previa 2',
                            'Distribución previa 3','Distribución previa 4'),
       col = colores,
       lty = 1, bty = 'n')

# ----------------------------- Reporte Gráfico --------------------------------

# Distribución posterior para la previa 1
plot(density(`Distribución previa 1`),
     main = TeX('Análisis de sensitividad $\\eta$'), col = 'orange',
     ylab = TeX('$p(\\eta | \\, y)'),
     xlab = TeX('$\\eta$'), bty = 'n', cex = 0.7)

# Distribución posterior para la previa 2 
lines(density(`Distribución previa 2`), col = 'blue')

# Distribución posterior para la previa 3
lines(density(`Distribución previa 3`), col = 'black')

# Distribución posterior para la previa 4
lines(density(`Distribución previa 4`), col = 'green')

# Explicación gráfica
legend("topleft",legend = c('Distribución previa 1','Distribución previa 2',
                            'Distribución previa 3','Distribución previa 4'),
       col = colores,
       lty = 1, bty = 'n')

abline(v = c(InfoParte1['Media',]), 
       col = colores, lty = 2)

abline(v = c(InfoParte1['Límite Inf',], InfoParte1['Límite Sup',]), 
       col = colores, lty = 'longdash')

text(x = c(-6.5,-14,-1.77), 0.08, 
     c('Medias', 'Límites Inferiores IC 95%', 'Límites superiores IC 95%'), 
     col = 'black',
     cex = 0.7,srt = -90,font = 2)

# ----------------------------- Reporte Tabular --------------------------------

InfoParte1 <- t(InfoParte1)
rownames(InfoParte1) = c('Distribución previa 1', 'Distribución previa 2',
                           'Distribución previa 3', 'Distribución previa 4')
print(InfoParte1)


# ========= PUNTO 3. BONDAD DE AJUSTE PARA LAS POBLACIONES SEPARADAS ===========

# Para revisar la bondad de ajuste lo que se hace es primero simular valores de la 
# distribución posterior (parámetro) y con cada uno de estos valores como parámetro
# de la exponencial se simulan datos sobre la población y se calculan estadísticos
# de prueba sobre estas. Si la distribución de estos estadísticos tiene valores usuales
# en los estadísticos de la población o muestra entonces diremos que nuestro modelo
# está bien ajustado a los datos.

AlphaBefore = 3 # Parámetros de la distribución previa 1.
BetaBefore = 17

set.seed(20122002)
GenDataCLEC = AjusteMCGIExp(dataCLEC,AlphaBefore,BetaBefore, n = 50000)
GenDataILEC = AjusteMCGIExp(dataILEC,AlphaBefore,BetaBefore, n = 50000)

# ================= CÓDIGO SOLO PARA LINUX, mcapply no funciona en windows ===== 
# # install.packages('parallel')
# library(parallel)         # De esta librería el comando mclapply es para ejecutar
# library(doParallel)       # el comando lapply de forma simultánea en varios elementos.
# cores <- detectCores()-1  # Usaremos paralelamente todos los núcleos menos uno.
# cl <- makeCluster(cores)
# registerDoParallel(cl)	
# 
# MediaMuestrasCLEC = unlist(mclapply(GenDataCLEC, 
#                                     FUN = function(x) 
#                                       mean(rexp(rate = 1/x, n = length(dataCLEC))),
#                                     mc.set.seed = 20122002, mc.cores = cores))
# DesvMuestrasCLEC = unlist(mclapply(GenDataCLEC, 
#                                    FUN = function(x) 
#                                      sd(rexp(rate = 1/x, n = length(dataCLEC))),
#                                    mc.set.seed = 20122002, mc.cores = cores))
# 
# MediaMuestrasILEC = unlist(mclapply(GenDataILEC, 
#                                     FUN = function(x) 
#                                       mean(rexp(rate = 1/x, n = length(dataILEC))),
#                                     mc.set.seed = 20122002, mc.cores = cores))
# DesvMuestrasILEC = unlist(mclapply(GenDataILEC, 
#                                    FUN = function(x) 
#                                      sd(rexp(rate = 1/x, n = length(dataILEC))),
#                                    mc.set.seed = 20122002, mc.cores = cores))
# ==============================================================================

set.seed(20122002)
MediaMuestrasCLEC = unlist(lapply(GenDataCLEC,
                                    FUN = function(x)
                                     mean(rexp(rate = 1/x, n = length(dataCLEC)))))

set.seed(20122002)
DesvMuestrasCLEC = unlist(lapply(GenDataCLEC,
                                   FUN = function(x)
                                     sd(rexp(rate = 1/x, n = length(dataCLEC)))))

set.seed(20122002)
MediaMuestrasILEC = unlist(lapply(GenDataILEC,
                                    FUN = function(x)
                                      mean(rexp(rate = 1/x, n = length(dataILEC)))))

set.seed(20122002)
DesvMuestrasILEC = unlist(lapply(GenDataILEC,
                                   FUN = function(x)
                                     sd(rexp(rate = 1/x, n = length(dataILEC)))))

# Tenga en cuenta que si bien las simulaciones para la media y la desviación estándar 
# se hacen por separado como definimos la semilla dentro de cada mclapply estamos usando
# las mismas muestras realmente.

# ----------------------------- Reporte Gráfico --------------------------------
#install.packages('ggplot2')
#install.packages('ggExtra')
library(ggplot2) # Librerías para las visualizaciones de los dispersogramas
library(ggExtra) # con los histogramas en los margenes.


# Dispersograma de el modelo ajustado para los clientes de otras empresas.
# ...........................................................................
captionMediaCLEC = paste0(' (Valor ppp = ', mean(MediaMuestrasCLEC < mean(dataCLEC)),')')
captionsdCLEC = paste0(' (Valor ppp = ', mean(DesvMuestrasCLEC < sd(dataCLEC)),')')

CLEC <- ggplot(as.data.frame(cbind('Media' = MediaMuestrasCLEC,'Desviación típica' = DesvMuestrasCLEC)), 
            aes(x = Media, y = `Desviación típica`)) +
  geom_point(shape = 18, color = "salmon2", size = 0.5) +
  geom_vline(xintercept = quantile(MediaMuestrasCLEC, 0.025),color = 'red', linetype = "dashed") + 
  geom_vline(xintercept = quantile(MediaMuestrasCLEC, 0.975),color = 'red', linetype = "dashed") + 
  geom_hline(yintercept = quantile(DesvMuestrasCLEC, 0.025),color = 'red', linetype = "dashed") + 
  geom_hline(yintercept = quantile(DesvMuestrasCLEC, 0.975),color = 'red', linetype = "dashed") + 
  labs(x = paste0("Media",captionMediaCLEC),y = paste0("Desviación típica",captionsdCLEC),title = 'Bondad de ajuste CLEC') +
  geom_point(data = data.frame(x = mean(dataCLEC), y = sd(dataCLEC)),aes(x = x, y = y),
             shape = 17, color = "#551A8B",size = 3)
ggMarginal(CLEC, type = "histogram", xparams = list(fill = 'salmon1', color = 'salmon'),
           yparams = list(fill = 'salmon1'), color = 'salmon')

# Dispersograma del modelo ajustado con los datos de los clientes de Verizon.
# ...........................................................................
captionMediaILEC = paste0(' (Valor ppp = ', mean(MediaMuestrasILEC < mean(dataILEC)),')')
captionsdILEC = paste0(' (Valor ppp = ', mean(DesvMuestrasILEC < sd(dataILEC)),')')

ILEC <- ggplot(as.data.frame(cbind('Media' = MediaMuestrasILEC,'Desviación típica' = DesvMuestrasILEC)), 
            aes(x = Media, y = `Desviación típica`)) +
  geom_point(shape = 18, color = "aquamarine3", size = 0.5) +
  geom_vline(xintercept = quantile(MediaMuestrasILEC, 0.025),color = 'blue', linetype = "dashed") + 
  geom_vline(xintercept = quantile(MediaMuestrasILEC, 0.975),color = 'blue', linetype = "dashed") + 
  geom_hline(yintercept = quantile(DesvMuestrasILEC, 0.025),color = 'blue', linetype = "dashed") + 
  geom_hline(yintercept = quantile(DesvMuestrasILEC, 0.975),color = 'blue', linetype = "dashed") + 
  labs(x = paste0("Media",captionMediaILEC),y = paste0("Desviación típica",captionsdILEC),title = 'Bondad de ajuste ILEC') +
  geom_point(data = data.frame(x = mean(dataILEC), y = sd(dataILEC)),aes(x = x, y = y),
             shape = 17, color = "#68838B",size = 3)
ggMarginal(ILEC, type = "histogram", xparams = list(fill = 'aquamarine2', color = 'aquamarine4'),
              yparams = list(fill = 'aquamarine2'), color = 'aquamarine4')

# ----------------------------- Reporte Tabular --------------------------------

InfoParte1 <- cbind(c(mean(MediaMuestrasILEC < mean(dataILEC)), 
                      mean(DesvMuestrasILEC < sd(dataILEC))),
                    c(mean(MediaMuestrasCLEC < mean(dataCLEC)),
                      mean(DesvMuestrasCLEC < sd(dataCLEC))))

colnames(InfoParte1) = c('ILEC','CLEC'); rownames(InfoParte1) = c('Media','Desviación estándar')

print('Los valores ppp se encuentran en la siguiente tabla: ')
print(InfoParte1)
