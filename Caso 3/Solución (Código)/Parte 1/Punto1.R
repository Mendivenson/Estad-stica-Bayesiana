# -------------------------------------------------------
# Caso de Estudio 3: Punto 1. Estadística Bayesiana (2023-I)
# Alcaldía de Bogotá 2023 (Inferencia)
# Fecha de Última Modificación: 17/Noviembre/2023
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchez
# -------------------------------------------------------

# ============================= BASE DE DATOS ==================================
encuesta = matrix(NA, ncol = 2, nrow = 10)
rownames(encuesta) = c('C. F. Galán', 'G. Bolívar', 'J. D. Oviedo', 'D. Molano',
                       'R. Lara', 'J. L. Vargas', 'J. E. Robledo', 'N. Ramos',
                       'R. A. Quintero', 'Voto en\nBlanco')
colnames(encuesta) = c('Cantidad','Proporción')

encuesta[,1] = c(493,257,227,48,41,38,28,11,3,54); encuesta[,2] = encuesta[,1]/sum(encuesta[,1])

# =================== EXPLORACIÓN DEL ESPACIO DE PARÁMETROS ====================
# Para la distribución Dirichlet se usa el siguiente paquete:
# install.package('LaplaceDemons')
library(LaplacesDemon)

# == INFORMACIÓN PREVIA:
n = sum(encuesta[,1])                             # Tamaño de la muestra
n. = encuesta[,"Cantidad"]                        # Vector de proporciones (Información externa)
k = nrow(encuesta)                                # No. de grupos.

# == DATOS PARA LA SIMULACIÓN: 
S = 105000; AnchoMuestreo = 10; quemar = 5000; aceptance = 0
param = matrix(NA,ncol= k + 1, nrow = (S-quemar)/AnchoMuestreo) 
colnames(param) = c(rownames(encuesta), 'alpha'); set.seed(2012)
begin = sample(1:AnchoMuestreo, size = 1); contador = 0

# == HIPERPARÁMETROS (y parámetro de ajuste)
a = b = 1; SqAjuste = 0.3

# INICIALIZANDO LA CADENA:
alpha.b = 5; theta.b = encuesta[,"Proporción"]

# Barra de progreso:
Progress = txtProgressBar(min = 0, max = S, title = 'Avance del algoritmo',style = 3)

# CADENA EN EL ESPACIO DE PARÁMETROS:
for (i in 1:S){
  # MUESTREADOR DE GIBBS:
  theta.b = rdirichlet(1, alpha = n. + alpha.b)
  
  # METROPOLIS-HASTING
  # ====> Propuesta:
  alpha.star = rgamma(1, shape = alpha.b^2 / SqAjuste, rate = alpha.b / SqAjuste)
  
  # ====> Cálculo del parámetro de ajuste:
  ajuste = ddirichlet(theta.b, alpha = rep(alpha.star, k), log = T) + 
    dgamma(alpha.star, shape = 1, rate = 1, log = T) - 
    ddirichlet(theta.b, alpha = rep(alpha.b, k), log = T) -
    dgamma(alpha.b, shape = 1, rate = 1, log = T) +
    dgamma(alpha.b, shape = alpha.star^2 / SqAjuste, rate = alpha.star / SqAjuste, log = T) -
    dgamma(alpha.star, shape = alpha.b^2 / SqAjuste, rate = alpha.b / SqAjuste, log = T) 
  
  if (runif(1) < exp(ajuste)){
    alpha.b = alpha.star
    aceptance = aceptance + 1
  }
  
  if (i > quemar & ((i - begin) %% AnchoMuestreo == 0)){
    contador = contador + 1
    param[contador,] = c(theta.b, alpha.b)
  }
  setTxtProgressBar(Progress, i, title = 'Avance del algoritmo')
}
close(Progress)
plot(param[,11], pch = '19', type = 'l')

# Note que como los componentes del vector de proporciones están siendo obtenidos
# mediante un muestreador de Gibbs, no es necesario hacer los siguientes análisis:

# TASA DE ACEPTACIÓN
100 * aceptance/S

# AUTOCORRELACIÓN
AutoCorrelacion = acf(param[,"alpha"],main = 'Autocorrelación')
print(AutoCorrelacion)

# ESTACIONARIEDAD DE LA CADENA
# install.packages('latex2exp')
library(latex2exp)
plot(param[,"alpha"], type = 'l', xlab = 'Iteración', ylab = TeX('\\alpha'),
     main = TeX('Cadena en el espacio de parámetros de \\alpha'),
     )

# TAMAÑO EFECTIVO DE MUESTRA:
coda::effectiveSize(param)

# COEFICIENTES DE VARIACIÓN DE MONTECARLO:
100 * (apply(param,MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(param)))/colMeans(param)

# ESTACIONARIEDAD DE LA CADENA COMPLETA:
plot(apply(param[,1:10], MARGIN = 1, FUN = function(x) dmultinom(n., size = 1200, prob = x, log = T)),
     type = 'l', main = 'Log-verosimilitud de la cadena', cex.main = 1.2, ylab = 'Log-verosimilitud',
     xlab = 'Iteración')

# ============================= INFERENCIA =====================================
realizaciones = param[,1:10]

# Las columnas de inferencia son realizaciones del vector de conteo.
#ealizaciones = t(realizaciones)
colnames(realizaciones) = rownames(encuesta)


inferencia = matrix(nrow = k, ncol = 3); rownames(inferencia) = rownames(encuesta)
colnames(inferencia) = c('Lim. Inf','Media','Lim. Sup')

inferencia[,"Lim. Inf"] = apply(realizaciones, MARGIN = 2, FUN = function(x) quantile(x,probs =  0.025))
inferencia[,"Media"] = colMeans(realizaciones)
inferencia[,"Lim. Sup"] = apply(realizaciones, MARGIN = 2, FUN = function(x) quantile(x,probs =  0.975))
inferencia = inferencia * 100
inferencia = inferencia[order(inferencia[,"Media"]),]

# ===> Datos oficiales de la registraduría;
registraduria = matrix(ncol = 2, nrow = k)
rownames(registraduria) = c('C. F. Galán', 'G. Bolívar', 'J. D. Oviedo', 'D. Molano',
                            'R. Lara', 'J. L. Vargas', 'J. E. Robledo', 'N. Ramos',
                            'R. A. Quintero', 'Voto en\nBlanco')
registraduria[,1] = c(1499734, 571948,616902,
                      65681,69639,31153,34099,
                      14855,2730,148794)
registraduria[,2] = 100 * registraduria[,1]/sum(registraduria[,1])
registraduria = registraduria[rownames(inferencia),]

# ===> GRÁFICO
plot(NA, NA, xlab = 'Porcentaje de votos', ylab = '',
     main = 'Inferencia bayesiana', xlim = c(0,50),
     cex.main = 1.6,
     ylim = c(1,k),  cex.axis = 0.9, yaxt = 'n')
abline(h = 1:k, col = "lightgray", lty = 1)
abline(v = seq(0,50, by = 5), col = "lightgray", lty = 1)
for (j in 1:k){
  segments(x0 = inferencia[j,'Lim. Sup'], y0 = j, x1 = inferencia[j,'Lim. Inf'], y1 = j, col = 'red')
  lines(x = inferencia[j,'Media'], y = j, type = "p", pch = 1, cex = 0.6, col = 'red')
  points(x = registraduria[j,2], y = j, pch = 2, col = 'blue')
}
axis(side = 2, at = 1:k, labels = rownames(inferencia), las = 2, cex.axis = 0.7)  
legend("bottomright", legend = c("Votaciones oficiales", 'Inferencia bayesiana', 'IC del 95% y media'),
       pch = c(2,21,NA), col = c("blue",'red','red'), bty = "n", cex = 0.8)

# ===> TABLA
# install.packages('knitr')
# library(knitr)
# inferencia = inferencia[order(inferencia[,"Media"], decreasing = T),]
# registraduria = registraduria[rownames(inferencia),]
# kable(round(cbind(inferencia, 'oficial' = registraduria[,2]),2), format = 'latex')
print(cbind(inferencia, 'oficial' = registraduria[,2]))

      
