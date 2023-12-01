# -------------------------------------------------------
# Caso de Estudio 3: Punto 2. Estadística Bayesiana (2023-I)
# Daatos de diabetes (Modelo 4 Errores correlacionados)
# Fecha de Última Modificación: 24/Noviembre/2023
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

# Para construir la matriz de correlación:
DY <-abs(outer( (1:n),(1:n) ,"-"))

lmAjustado = lm(y ~ X - 1)

# S = 2000

# ==================== CONSTANTES DEL MUESTREADOR ==============================
S = 201000                                                          # Iteraciones
AnchoMuestreo = 10                                                  # Ancho muestreo sistemático
quemar = 1000                                                       # Período de precalentamiento
begin = sample(1:10, size = 1)                                      # Inicio muestreo sistemático
param = matrix(ncol = (p+2), nrow = (S - quemar)/AnchoMuestreo)     # Matriz para guardar lo muestreado
colnames(param) = c( colnames(X), 'sigma', 'lambda')
# contador = 0                                                        # Para poder agregar a la matriz de muestras las muestras
ProgressBar = txtProgressBar(min = 0, max = S, style = 3)           # Barra de progreso

# ========================== GIBBS-HASITING POR BLOQUES ========================

# install.packages('mvtnorm')
# Valores iniciales del algoritmo:
beta = lmAjustado$coefficients
sigma2 = summary(lmAjustado)$sigma^2
phi = acf(lmAjustado$res, plot = FALSE)$acf[2]

# ---> Regresión con errores correlacionados:

# Hiperparamétros:
nu_0 = 1
sigma2_0 = 1
T0 = diag(1/1000, nrow = p)

set.seed(2012)
ac = 0
OUT = NULL
for (i in 1:S){
  # Cor = phi^DY                                               # Matriz de correlación
  iCor = solve(phi^DY )                                          # Inversa de la matriz de correlación 
  V.beta = solve( t(X)%*%iCor%*%X/sigma2 + T0)
  # E.beta = V.beta%*%( t(X)%*%iCor%*%y/sigma2)
  beta   = as.numeric(mvtnorm::rmvnorm(1,V.beta%*%( t(X)%*%iCor%*%y/sigma2),V.beta))
  sigma2 = 1/rgamma(1,(nu_0+n)/2,(nu_0*sigma2_0+t(y-X%*%beta)%*%iCor%*%(y-X%*%beta)) /2 )
  phi.p <- abs(runif(1,phi-.28,phi+.28))
  phi.p <- min(phi.p, 2-phi.p)
  
  lr <- -.5*(determinant(phi.p^DY,log=TRUE)$mod - determinant(phi^DY,log=TRUE)$mod + 
                sum(diag((y-X%*%beta)%*%t(y-X%*%beta)%*%(solve(phi.p^DY) - solve(phi^DY)))/sigma2))
  
  if( log(runif(1)) < lr ) { 
    phi = phi.p
    ac = ac+1 
  }
  
  if( i > quemar & ((i - begin) %% AnchoMuestreo == 0)){
    OUT = rbind(OUT,c(beta,sigma2,phi))
  }
  setTxtProgressBar(ProgressBar, i)
  # OUT <- rbind(OUT,c(beta,sigma2,phi)) 
}
close(ProgressBar)

colnames(OUT) = c(colnames(X), 'sigma', 'phi')
# Creando la carpeta modelos y guardando la información obtenida:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (file.exists('Modelos') == F){dir.create('Modelos')}
write.table(OUT, file = 'Modelos/Modelo4')
