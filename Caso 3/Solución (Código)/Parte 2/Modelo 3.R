# -------------------------------------------------------
# Caso de Estudio 3: Punto 2. Estadística Bayesiana (2023-I)
# Daatos de diabetes (Modelo 3: Regresión rígida)
# Fecha de Última Modificación: 17/Noviembre/2023
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchez
# -------------------------------------------------------

# ================================ BASE DE DATOS ===============================
trainData = source('http://www2.stat.duke.edu/~pdh10/FCBS/Inline/yX.diabetes.train')$value
str(trainData)

# =========================== CONSTANTES NECESARIAS ============================
y = trainData[,"y"]
X = as.matrix(trainData[,-c(1)])
n = length(y)
p = ncol(X)

beta_mco = solve(t(X) %*% X) %*% t(X) %*% y
  
# ---> Regresión rígida:

# Hiperparámetros sigma^2
nu_0 = 1
sigma2_0 = sum((y - X %*% beta_mco)^2)/(n-p)

# Hiperparámetros lambda
a = 1
b = 2

xTx = t(X) %*% X
xTy = t(X) %*% y

# ==================== CONSTANTES DEL MUESTREADOR ==============================
S = 201000                                                          # Iteraciones
AnchoMuestreo = 10                                                  # Ancho muestreo sistemático
quemar = 1000                                                       # Período de precalentamiento
begin = sample(1:10, size = 1)                                      # Inicio muestreo sistemático
param = matrix(ncol = (p+2), nrow = (S -quemar)/AnchoMuestreo)      # Matriz para guardar lo muestreado
colnames(param) = c( colnames(X), 'sigma', 'lambda')       
contador = 0                                                        # Para poder agregar a la matriz de muestras las muestras
ProgressBar = txtProgressBar(min = 0, max = S, style = 3)           # Barra de progreso

# ======================== MUESTREADOR DE GIBBS ================================

# install.packages('mvtnorm')
# Inicialización del algoritmo:
lambda.star =  rgamma(n = 1, a, b)
sigma2.star = 1/rgamma(n = 1, 0.5 * nu_0, nu_0 * sigma2_0 * 0.5)
beta.star = beta_mco

set.seed(2012)
for (i in 1:S){
  Sigma = solve(diag(x = lambda.star, nrow = p) + xTx)
  beta.star = as.numeric(mvtnorm::rmvnorm(n = 1, 
                                          mean = Sigma %*% t(X) %*% y,
                                          sigma = sigma2.star * Sigma))
  betaTbeta = sum(beta.star^2)
  SSR = sum((y - X %*% beta.star)^2)
  
  lambda.star = rgamma(n = 1, 0.5 * (p + 2 * a), b + 0.5 * (betaTbeta))
  
  sigma2.star = 1/rgamma(n = 1, 0.5 * (n + p + nu_0), 
                        0.5 * (SSR + lambda.star * betaTbeta + nu_0 * sigma2_0))
  if( i > quemar & ((i - begin) %% AnchoMuestreo == 0)){
    contador = contador + 1
    param[contador,] = c(beta.star, sigma2.star, lambda.star)
  }
  setTxtProgressBar(ProgressBar,value = i)
}
close(ProgressBar)

# Creando la carpeta modelos y guardando la información obtenida:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (file.exists('Modelos') == F){dir.create('Modelos')}
write.table(param, file = 'Modelos/Modelo3')
