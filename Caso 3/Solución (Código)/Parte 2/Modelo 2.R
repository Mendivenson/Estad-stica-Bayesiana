# -------------------------------------------------------
# Caso de Estudio 3: Punto 2. Estadística Bayesiana (2023-I)
# Datos de diabetes (Modelo 2: Previa g)
# Fecha de Última Modificación: 17/Noviembre/2023
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchez
# -------------------------------------------------------

# ================================ BASE DE DATOS ===============================
trainData = source('http://www2.stat.duke.edu/~pdh10/FCBS/Inline/yX.diabetes.train')$value
str(trainData)

# ============================== INFORMACIÓN MCO ===============================
y = trainData[,"y"]
X = as.matrix(trainData[,-c(1)])
n = length(y)
p = ncol(X)

# ---> Previa g (usando a g = n, para una previa no informativa)
beta_mco = solve(t(X) %*% X) %*% t(X) %*% y       # Nos ayudará a calcular la media para generar beta
C = solve(t(X) %*% X)
g = n/(n+1)
SSRg = t(y) %*% y - g * t(y) %*% X %*% beta_mco

nu0 = 1; sigma2_0 = sum((y - X %*% beta_mco)^2)/(n-p)
# sigma2_0 = 100 ¡HACER ANÁLISIS DE SENSITIVIDAD!
# ==================== CONSTANTES DEL MUESTREADOR ==============================
S = 20000                                                           # Iteraciones
param = matrix(ncol = (p+1), nrow = S)                              # Matriz para guardar lo muestreado
colnames(param) = c( colnames(X), 'sigma')       
ProgressBar = txtProgressBar(min = 0, max = S, style = 3)           # Barra de progreso

# ============================ MONTE CARLO =====================================

# Inicialización del algoritmo:
sigma2.star = 1/rgamma(1, 0.5 * (nu0 + n), 0.5 * (nu0 * sigma2_0 + SSRg))
beta.star = beta_mco

set.seed(2012)
# install.packages('mvtnorm')
for (i in 1:S){
  sigma2.star = 1/rgamma(1, (nu0 + n)/2, (nu0 * sigma2_0 + SSRg)/2)
  beta.star = mvtnorm::rmvnorm(1, mean = g * beta_mco, 
                                sigma = g * sigma2.star * C)
  param[i,] = c(beta.star, sigma2.star)
  setTxtProgressBar(ProgressBar,value = i)
}
close(ProgressBar)

# Creando la carpeta modelos y guardando la información obtenida:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (file.exists('Modelos') == F){dir.create('Modelos')}
write.table(param, file = 'Modelos/Modelo2')
