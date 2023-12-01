# -------------------------------------------------------
# Caso de Estudio 3: Punto 2. Estadística Bayesiana (2023-I)
# Datos de diabetes (Modelo 1: previa unitaria)
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
# X = cbind(1,trainData[,-c(1)]) Las variables ya están centradas
X = trainData[, -c(1)]
n = length(y)
p = ncol(X)

# ---> Previa unitaria (Sólo se modifican los hiperparámetros)
beta0 = solve(t(X) %*% X) %*% t(X) %*% y                          # Beta por mínimos cuadrados
sigma2_0 = sum((y - X%*%beta0)^2)/(n-p)                           # MSE por mínimos cuadrados
Sigma0inv = solve(n * sigma2_0 * solve(t(X) %*% X))               # Varianza de los betas
nu0 = 1                                                           # nu0

# ==================== CONSTANTES DEL MUESTREADOR ==============================
S = 201000                                                          # Iteraciones
AnchoMuestreo = 10                                                  # Ancho muestreo sistemático
quemar = 1000                                                       # Período de precalentamiento
begin = sample(1:10, size = 1)                                      # Inicio muestreo sistemático
param = matrix(ncol = (p+1), nrow = (S -quemar)/AnchoMuestreo)      # Matriz para guardar lo muestreado
colnames(param) = c( colnames(X), 'sigma')       
contador = 0                                                        # Para poder agregar a la matriz de muestras las muestras
ProgressBar = txtProgressBar(min = 0, max = S, style = 3)           # Barra de progreso

# ======================== MUESTREADOR DE GIBBS ================================

# Inicialización del algoritmo:
beta.star = as.numeric(beta0)                                       # Inicializamos el vector beta en el estimador por MCO
sigma2.star = sigma2_0                                              # Inicializamos a sigma^2 con MSE (Proveniente de MCO)

# install.packages('mvtnorm')
xTx = t(X) %*% X
xTy = t(X) %*% y

set.seed(2012)
for (i in 1:S){
  VarCovar = solve(Sigma0inv + xTx * (sigma2.star)^-1)
  beta.star = as.numeric(mvtnorm::rmvnorm(n = 1,
                               mean = VarCovar %*% (Sigma0inv %*% beta0 + xTy * (sigma2.star)^-1),
                               sigma = VarCovar))
  sigma2.star = 1/rgamma(1, 0.5 * (n + nu0), 0.5 * (nu0* sigma2_0 + sum((y - X %*% beta.star)^2)))
  if( i > quemar & ((i - begin) %% AnchoMuestreo == 0)){
    contador = contador + 1
    param[contador,] = c(beta.star, sigma2.star)
  }
  setTxtProgressBar(ProgressBar,value = i)
}
close(ProgressBar)

# Creando la carpeta modelos y guardando la información obtenida:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (file.exists('Modelos') == F){dir.create('Modelos')}
write.table(param, file = 'Modelos/Modelo1')
