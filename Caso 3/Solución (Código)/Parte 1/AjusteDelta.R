# -------------------------------------------------------
# Caso de Estudio 3: Punto 1. Estadística Bayesiana (2023-I)
# Alcaldía de Bogotá 2023 (Elección parámetro de ajuste)
# Fecha de Última Modificación: 17/Noviembre/2023
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchez
# -------------------------------------------------------

# ============================= BASE DE DATOS ==================================
encuesta = matrix(NA, ncol = 2, nrow = 10)
rownames(encuesta) = c('C. F. Galán', 'G. Bolívar', 'J. D. Oviedo', 'D. Molano',
                       'R. Lara', 'J. L. Vargas', 'J. E. Robledo', 'N. Ramos',
                       'R. A. Quintero', 'Voto en Blanco')
colnames(encuesta) = c('Cantidad','Proporción')
encuesta[,1] = c(493,257,227,48,41,38,28,11,3,54);encuesta[,2] = encuesta[,1]/sum(encuesta[,1])

# ====================== Elección del parámetro de ajuste ======================
# Para la distribución Dirichlet se usa el siguiente paquete:
# install.package('LaplaceDemons')
library(LaplacesDemon)

# == DATOS NECESARIOS:
n = sum(encuesta[,1])                             # Tamaño de la muestra
n. = encuesta[,"Cantidad"]                        # Vector de proporciones (Información externa)
k = nrow(encuesta)                                # No. de grupos.

# == PARA GUARDAR LOS DATOS:
S = 1000; VectorAceptance = as.numeric()

# == HIPERPARÁMETROS:
a = b = 1; Ajuste = c(0.1,0.2,0.3,0.4,0.5,0.6)

# INICIALIZANDO EL ALGORITMO:
alpha.b = 5; theta.b = encuesta[,"Proporción"]

for (j in Ajuste){
  SqAjuste = j
  Aceptance = 0
  Progress = txtProgressBar(min = 0, max = S, title = 'Avance del algoritmo',
                            style = 3)
  set.seed(2012)
  # CADENA EN EL ESPACIO DE PARÁMETROS:
  for (i in 1:S){
    
    # MUESTREADOR DE GIBBS:
    theta.b = rdirichlet(1, alpha = n. + alpha.b)
    
    # METROPOLIS-HASTING
    # ====> Propuesta:
    alpha.star = rgamma(1, shape = alpha.b^2 / SqAjuste, rate = alpha.b / SqAjuste)
    
    # TASA DE ACEPTACIÓN
    ajuste = ddirichlet(theta.b, alpha = rep(alpha.star, k), log = T) + 
      dgamma(alpha.star, shape = 1, rate = 1, log = T) - 
      ddirichlet(theta.b, alpha = rep(alpha.b, k), log = T) -
      dgamma(alpha.b, shape = 1, rate = 1, log = T) +
      dgamma(alpha.b, shape = alpha.star^2 / SqAjuste, rate = alpha.star / SqAjuste, log = T) -
      dgamma(alpha.star, shape = alpha.b^2 / SqAjuste, rate = alpha.b / SqAjuste, log = T) 
    
    if (runif(1) < exp(ajuste)){
      alpha.b = alpha.star
      Aceptance = Aceptance + 1
    }
    setTxtProgressBar(Progress, i, title = 'Avance del algoritmo')
  }
  close(Progress)
  VectorAceptance = c(VectorAceptance, Aceptance)
}

VectorAceptance = 100 * VectorAceptance/S
names(VectorAceptance) = Ajuste
print(VectorAceptance)
# Por los resultados se puede escoger a delta como 0.2, 0.3, 0.4 o 0.5. En el caso
# de este grupo de trabajo se escoje e delta = 0.3