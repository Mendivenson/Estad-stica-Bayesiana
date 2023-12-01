# -------------------------------------------------------
# Caso de Estudio 1: Parte 3. Estadística Bayesiana (2023-I)
# Análisis de Tiempo de Reparación (Verizon) - Simulación
# Fecha de Última Modificación: 2023-09-15
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastian Gil Sanchéz
# -------------------------------------------------------

# En este script todo el código está contenido en una función 'metodo' la cual
# es necesaria para usar parLapply. 

setwd("C:/Users/gsgil/Documents/sebastiangils/Semestre_9/Bayesiana/CasosDeEstudio/Caso1/")

library(readr)
library(invgamma)
library(parallel)
library(doParallel)

Metodo <- function(n , m = 100000, mBoots= 10000){
  # n <- tamaño de muestras
  # nBoots <- multiplicador al tamaño de las muestras Bootstrap (size = nBoots*n)
  # m <- numero de muestras
  # mBoots <- numero de muestras Bootstrap
  library(readr)
  library(invgamma)
  library(parallel)
  library(doParallel)
  # Informacion original del problema
  # 1: ILEC - 2: CLEC
  verizon <- read.csv("Verizon.csv")
  dataILEC <- verizon[verizon$Group == 'ILEC',][['Time']]
  dataCLEC <- verizon[verizon$Group == 'CLEC',][['Time']]
  y1 = mean(dataILEC) 
  y2 = mean(dataCLEC)
  eta = y1 - y2
  # Variables procedimentales
  muestras = MediaEta = MediaEta1 = MediaEta2 = BNP = NULL
  InBootP = InBootNoP = InBayes = InAsintotic = 0
  
  
  # ------------------ Creación de muestras ---------------------------------
  set.seed(4321)
  for (i in 1:m){
    # Crear las nuevas muestras de tamaño 'n'
    Muestra1 <- rexp(n, 1/y1)
    Muestra2 <- rexp(n, 1/y2)
    # Bootstrap no parametricos
    for (i in 1:mBoots){
      MediaEta[i] = mean( sample(Muestra1, size = 0.5*n,replace = T) -
                          sample(Muestra2, size = 0.5*n,replace = T) )
    }
    # Limites de IC al 95% de analisis Frecuentista Bootstrap No Parametrico
    LimitFrecNoPara <- c(quantile(MediaEta, 0.025), quantile(MediaEta,0.975))
    if (eta >= LimitFrecNoPara[1] && eta <=  LimitFrecNoPara[2]) {
      InBootNoP = 1
    } 
    muestras <- rbind(muestras, c(sum(Muestra1), sum(Muestra2), InBootNoP))
  }
  intervalos <- NULL
  set.seed(4321)
  for (i in 1:m){
# --------------------- Modelo Bayesiano - Previa 1 ---------------------------
    # Hiperparametros  
    AlphaBefore = 3
    BetaBefore = 17
    # Distribución posterior de 'eta'
    PostMC1 <- rinvgamma(1000, shape = AlphaBefore + n ,
                         rate = BetaBefore + muestras[i,1])
    PostMC2 <- rinvgamma(1000, shape = AlphaBefore + n ,
                         rate = BetaBefore + muestras[i,2])
    PostEta <- PostMC1 - PostMC2
    # Limites de IC al 95% de analisis Bayes
    LimitBayes <- c(quantile(PostEta, 0.025), quantile(PostEta,0.975))
    
# ------------------- Normalidad Asintotica ----------------------------------
    ## Definición de los parámetros para las normales asintóticas para cada población.
    muCLEC = muestras[i,2]/n ; muILEC= muestras[i,1]/n
    InvInfCLEC = (muCLEC^2)/n;InvInfILEC = (muILEC^2)/n
    ## Definición nuevos parámetros
    mu = muILEC - muCLEC; var = InvInfCLEC + InvInfILEC
    rm(muILEC,muCLEC,InvInfCLEC,InvInfILEC)
    # Limites de IC al 95% de Analisis Frecuentista Asintotico
    LimitAsintotic <- c(round(qnorm(0.025,mean = mu, sd = sqrt(var)),2) ,
                        round(qnorm(0.975,mean = mu, sd = sqrt(var)),2))
    
# ------------------ Bootstrap parametrico -----------------------------
    lambda1 = 1 / (muestras[i,1]/n) 
    lambda2 = 1 / (muestras[i,2]/n) 
    MediaEta2 = NULL
    set.seed(4321)
    for (i in 1:mBoots){
      MediaEta2[i] = mean(rexp(n, rate = lambda1) - rexp(n, rate = lambda2))
    }
    # Limites de IC al 95% de analisis Frecuentista Bootstrap Parametrico
    LimitFrecPara <- c(quantile(MediaEta2, 0.025), quantile(MediaEta2,0.975))
    
#---------------------- Verificacion de etas en IC -----------------------------
    InBayes = 0
    InAsintotic = 0
    InBootP = 0
    # Verificar Eta al intervalo
    if (eta >= LimitBayes[1] && eta <= LimitBayes[2]) {
      InBayes = 1
    } 
    if (eta >= LimitAsintotic[1] && eta <= LimitAsintotic[2]) {
      InAsintotic = 1
    } 
    if (eta >= LimitFrecPara[1] && eta <= LimitFrecPara[2]) {
      InBootP = 1
    }
    intervalos <- rbind(intervalos, c(InBayes, InAsintotic, InBootP))
  }
  
  Escenario <- cbind(muestras, intervalos)
  as.data.frame(Escenario)
  colnames(Escenario) = c("M1","M2","BootNoP","Bayes","Asintotic","BootP")
  Proportion <- Escenario[,3:6]
  Proportion <- Proportion[,c("Bayes","Asintotic","BootP","BootNoP")]
  Proportion <- colSums(Proportion)/m
  return(Proportion)
}

# Prueba con un escenario
E1 <- Metodo(n = 10);E1
E4 <- Metodo(n = 100);E4

# ---------------------- Aplicación en paralelo de la función ------------------
# Lista de parametros, es decir, tamaño de muestras de cada escenario
parameters <- list(10, 20, 50, 100)  # Example parameters
# Incializar el cluster (en Windows)
cl <- makeCluster(detectCores())
# Ejecutar de manera paralela con parLapply 
results <- parLapply(cl, parameters, Metodo)
# Resultados
results
# Para el cluster 
stopCluster(cl)

