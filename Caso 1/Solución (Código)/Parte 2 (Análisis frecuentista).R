# -------------------------------------------------------
# Caso de Estudio 1: Parte 2. Estadística Bayesiana (2023-I)
# Análisis de Tiempo de Reparación (Verizon)
# Fecha de Última Modificación: 2023-09-11
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchéz
# -------------------------------------------------------

# En este script está el código necesario para hacer las mismas estimaciones sobre eta
# del punto anterior, pero mediante modelos frecuentistas: Bootstrap, Bootstrap 
# paramétrico y valiendonos de la normalidad asintótica del MLE.

# ---------------- TRATAMIENTO PREVIO DE LA INFORMACIÓN ------------------------

rm(list = ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv('Verizon.csv', header = T, sep = ',')
dataILEC <- data[data$Group == 'ILEC',][['Time']]
dataCLEC <- data[data$Group == 'CLEC',][['Time']]
rm(data)

# Se usa la semilla 20122002 para que el código sea reproducible en todos los 
# lugares necesarios del siguiente script.

# ==================== PUNTO 1. NORMALIDAD ASINTÓTICA DEL MLE ==================

# Como lo que queremos es estimar sobre Eta modelamos cada población exponencial
# por separado y restamos las normales resultantes sobre cada distribución.

# Definición de los parámetros para las normales asintóticas para cada población.
muCLEC = mean(dataCLEC); muILEC = mean(dataILEC)                                  # El estimador MLE para el caso exponencial es la media.
InvInfCLEC = (muCLEC^2)/length(dataCLEC);InvInfILEC = (muILEC^2)/length(dataILEC) # Cálculo de la información observada de Fisher.

# La media y la varianza de la resta de normales respectivamente.
mu = muILEC - muCLEC; var = InvInfCLEC + InvInfILEC
rm(muILEC,muCLEC,InvInfCLEC,InvInfILEC)

# ----------------------------- Reporte Tabular --------------------------------
print(cbind('Media' = round(mu,2), 
      'Coef. Variación' = paste0(round(100 * (sqrt(var)/abs(mu)),2),'%'), 
      'Lím. Inferior IC 95%' = round(qnorm(0.025,mean = mu, sd = sqrt(var)),2),
      'Lím. Superior IC 95%' = round(qnorm(0.975,mean = mu, sd = sqrt(var)),2)))

# ----------------------------- Reporte Gráfico --------------------------------
# install.packages('latex2exp')
library(latex2exp)
# Densidad
curve(expr = dnorm(x, mean= mu, sd = sqrt(var)),from = -17,to = 2,
      main = 'Análisis frecuentista',ylab = 'Densidad',xlab = TeX('$\\eta = \\lambda_1 - \\lambda_2$'))

# Intervalo de confianza 95% y media
abline(v = c(qnorm(0.025,mean = mu, sd = sqrt(var)),qnorm(0.9755,mean = mu, sd = sqrt(var))), col = 'red', lty = 2)
abline (v = mu, col = 'blue', lty = 3)

# Ayuda para la interpretación de la gráfica
mtext(text = TeX('\\textbf{usando la normalidad asintótica del MLE}'),side = 3, line = 0.5, cex = 1, font = 2)
legend('topright',legend = c('Media', 'IC 95%'),col = c('blue','red'),lty = 2, cex = 0.8)

# ====================== PUNTO 2. BOOTSTRAP NO PARAMÉTRICO =====================

# En este caso, hacemos bootstrap no paramétrico para cada población por separado
# y estimamos eta como la diferencia de las medias de estas estmaciones.

set.seed(20122002)
mBoots = 50000                        # Se hacen 5000 remuestras Bootstrap
MediaEta = NULL                       # Vector para guardar las estimaciones de Eta 
for (i in 1:mBoots){
  GenDataCLEC <- sample(dataCLEC, size = 1000,replace = T) # Remuestreo CLEC
  GenDataILEC <- sample(dataILEC, size = 1000,replace = T) # Remuestreo ILEC
  MediaEta[i] = mean(GenDataILEC) - mean(GenDataCLEC)
}

# ----------------------------- Reporte Gráfico --------------------------------
# Histograma
hist(MediaEta,col = 'gray90',border = 'gray90',main = 'Análisis frecuentista',
     prob = T,ylab = 'Densidad',xlab = TeX('$\\eta = \\lambda_1 - \\lambda_2$'))

# Densidad
lines(density(MediaEta), col = '#8B668B') 

# Media e intervalo de confianza 95%
abline(v = mean(MediaEta), col = 'red', lty = 2)
abline(v = c(quantile(MediaEta,0.025),quantile(MediaEta,0.975)), col = 'blue', lty = 2)

# Ayuda para la interpretación de la gráfica
mtext(text = TeX('\\textbf{usando Bootstrap no paramétrico}'),side = 3, line = 0.5, cex = 1, font = 2)
legend('topright',legend = c('Media', 'IC 95%'),col = c('blue','red'),lty = 2, cex = 0.8)

# ----------------------------- Reporte Tabular --------------------------------
q1 = quantile(MediaEta,0.025);names(q1) = NULL
q2 = quantile(MediaEta, 0.975);names(q2) = NULL
print(cbind('Media' = round(mean(MediaEta),2),
            'Coef. Variación' = paste0(round(100 * (sd(MediaEta)/abs(mean(MediaEta))),2),'%'),
            'Lím. Inferior IC 95%' = round(q1,2),'Lím. Superior IC 95%' = round(q2,2)))

# ======================== PUNTO 3. BOOTSTRAP PARAMÉTRICO ======================

# En este caso lo que se hace es muestrear directamente de las poblaciones suponiendo 
# una distribución exponencial de parámetros de razón como los estimadores MLE y para 
# estimar eta se hace la resta de las medias de cada una de las muestras nuevamente.

lambda1 = 1 / mean(dataCLEC)                    # Estimador MLE población CLEC
lambda2 = 1 / mean(dataILEC)                    # Estimador MLE población ILEC
B = 50000
MediaEta = NULL
set.seed(20122002)
for (i in 1:B){
  remuestraCLEC = rexp(1000, rate = lambda1)
  remuestraILEC = rexp(1000, rate = lambda2)
  remuestraEta = remuestraILEC - remuestraCLEC
  MediaEta[i] = mean(remuestraEta)
}

# ----------------------------- Reporte Gráfico --------------------------------
# Histograma
hist(MediaEta,col = 'gray90',border = 'gray90',main = 'Análisis frecuentista',prob = T,
     ylab = 'Densidad',ylim = c(0,0.7),xlab = TeX('$\\eta = \\lambda_1 - \\lambda_2$'))

# Densidad
lines(density(MediaEta), col = '#8B668B')

# Media e intervalo de confianza 95%
abline(v = mean(MediaEta), col = 'red', lty = 2)
abline(v = c(quantile(MediaEta,0.025),quantile(MediaEta,0.975)), col = 'blue', lty = 2)

# Ayuda para la interpretación gráfica
mtext(text = TeX('\\textbf{usando Bootstrap paramétrico}'),side = 3, line = 0.5, cex = 1, font = 2)
legend('topright',legend = c('Media', 'IC 95%'),col = c('blue','red'),lty = 2, cex = 0.8)

# ----------------------------- Reporte Tabular --------------------------------
q1 = quantile(MediaEta,0.025);names(q1) = NULL
q2 = quantile(MediaEta, 0.975);names(q2) = NULL
print(cbind('Media' = round(mean(MediaEta),2),
            'Coef. Variación' = paste0(round(100 * (sd(MediaEta)/abs(mean(MediaEta))),2),'%'),
            'Lím. Inferior IC 95%' = round(q1,2),
            'Lím. Superior IC 95%' = round(q2,2)))