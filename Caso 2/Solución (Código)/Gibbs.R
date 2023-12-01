# -------------------------------------------------------
# Caso de Estudio 2: Muestreadores de Gibbs. Estadística Bayesiana (2023-I)
# Prueba Saber 11 2022-2: Una perspectiva multinivel
# Fecha de Última Modificación: 19/Octubre/2023
# Autores: 
#         - Michel Mendivenson Barragán Zabala
# -------------------------------------------------------


# ============================= SIMULACIÓN =====================================

## -----------------------  M1: MODELO NORMAL. ---------------------------------
Gibbs1 <- function(y,                                                          # Vector de datos
               mu_0 = 250,gammaSq_0 = 2500, nu_0 = 1, sigmaSq_0 = 2500,         # Hiperparámetros
               NoDeIteraciones = 101000, Quemar = 1000, AnchoMuestreo = 10){  
  
  # Antes de empezar con el muestreador necesitamos inicializar:
  Progress = txtProgressBar(min = 0, max = NoDeIteraciones,                     # Barra de progreso
                            style = 3)
  muestras = matrix(NA, nrow = (NoDeIteraciones - Quemar)/AnchoMuestreo,        # Matriz donde se guardarán las muestras
                    ncol = 3)
  
  colnames(muestras) = c('theta', 'sigma2', 'LogVerosimilitud')
  
  begin = sample(1:AnchoMuestreo, size = 1)                                     # Número de partida para el muestreo sistemático
  
  n = length(y)                                                                 # Tamaño del vector de datos
  SumaTotal = sum(y)                                                            # Suma total de los datos
  
  theta = mean(y)                                                               # Valor de theta (Estimación máxima verosimilitud)
  sigma2 = var(y)                                                               # Valor de sigma^2 (Estimación máxima verosimilitud)
  
  # El parámetro shape para la dist cond. compl. de sigma2 es siempre constante.
  
  contador = 0
  # ==== MUESTREADOR DE GIBBS
  for (i in 1:NoDeIteraciones){
    
    # Actualización de theta
    theta = rnorm(n = 1,
                  mean = (mu_0/gammaSq_0 + SumaTotal/sigma2)/
                    (1/gammaSq_0 + n/sigma2),
                  sd = 1/sqrt(1/gammaSq_0 + n/sigma2))
    
    # Actualización de sigma^2
    sigma2 = 1/rgamma(n = 1, shape = 0.5*(n + nu_0), rate = 0.5 * ((nu_0 * sigmaSq_0) + sum((y - theta)^2)))
    
    # Condición ciclo de quemado y muestreo sistemático
    if (i > Quemar & ((i - begin) %% AnchoMuestreo == 0)){
      contador = contador + 1
      muestras[contador,] = c(theta,sigma2, 
                              sum(dnorm(x = y,mean = theta,sd = sqrt(sigma2),
                                        log = T)))
    }
    setTxtProgressBar(Progress, i)
  }
  close(Progress)
  beepr::beep()
  return(muestras)
}

## ---------- M2: Modelo Normal con medias específicas por departamento. ----------

Gibbs2 = function(medias_j, n_j, y, names,                                     # media_j es el vector de medias y n_j es el vector de totales espcíficos
                   nu_0 = 1, sigmaSq_0 = 2500, mu_0 = 250,                      # Hiperparámetros
                   tauSq_0 = 2500, eta_0 = 1, gammaSq_0 = 2500,
                   NoDeIteraciones = 101000, Quemar = 1000, AnchoMuestreo = 10){
  # Antes de iniciar el muestreador:
  Progress = txtProgressBar(min = 0, max = NoDeIteraciones,
                            style = 3)
  m = length(medias_j)
  
  muestras = matrix(NA, ncol = m + 4, 
                    nrow = (NoDeIteraciones - Quemar)/AnchoMuestreo)
  colnames(muestras) = c(names, 'sigma2', 'mu', 'tauSq', 'LogVerosimilitud')
  
  begin = sample(1:AnchoMuestreo, size = 1)                                     # Número de partida para el muestreo sistemático
  
  SumaTotales = n_j * medias_j                                                  # La suma de los puntajes globales por departamentos
  # SumaCuadradosTotales = (n_j-1)*s2                                           # La suma de la diferencia
  Total = sum(n_j)
  Vmedias_j = medias_j                                                          # Vector de medias por departamento
  sigma2 = 2690.825                                                             # Varianza conjunta (Tomada de var(puntajes$puntaje))
  mu = mean(medias_j)                                                           # Media de las medias de los departamentos
  tauSq = var(medias_j)                                                         # Varianzas de las medias de los departamentos
  
  contador = 0
  # MUESTREADOR DE GIBBS:
  for (i in 1:NoDeIteraciones){
    
    # Actualizando el vector de medias
    Vmedias_j = rnorm(n = m, mean = 
                        (SumaTotales/sigma2 + mu/tauSq)/(1/tauSq + n_j/sigma2),
                         sd = 1/sqrt(1/tauSq + n_j/sigma2))
    # Actualizando  sigma^2
    sigma2 = 1/rgamma(n = 1, shape = 0.5 * (nu_0 + Total), 
                      rate = 0.5 * (nu_0 * sigmaSq_0 + 
                              sum((y - rep(Vmedias_j, n_j))^ 2)))
    # Actualizando mu
    mu = rnorm(n = 1, 
               mean = (mu_0/gammaSq_0 + sum(Vmedias_j)/tauSq)/
                 (1/gammaSq_0 + m/tauSq),
               sd = 1/sqrt(1/gammaSq_0 + m/tauSq))
    
    # Actualizando tau
    tauSq = 1/rgamma(n = 1, 
                     shape = 0.5 * (eta_0 + m),
                     rate = 0.5*(eta_0 * tauSq_0 + sum((Vmedias_j - mu)^2)))
    
    if (i > Quemar & ((i - begin) %% AnchoMuestreo == 0)){
      contador = contador + 1
      muestras[contador,] = c(Vmedias_j,sigma2, mu, tauSq, sum(dnorm(x = y,mean = rep(Vmedias_j, n_j),sd = sqrt(sigma2),log = T)))
    }
    setTxtProgressBar(Progress, i)
  }
  close(Progress)
  beepr::beep()
  
  return(muestras)
}

## --- M3: Modelo Normal con medias y varianzas específicas por departamento ---

Gibbs3 = function(y, medias, varianzas, n_j, names,
                  mu_0 = 250, gammaSq_0 = 2500, eta_0 = 1, tauSq_0 = 2500,
                  nu = 1, alpha_0 = 1, beta_0 = 1/2500,
                  NoDeIteraciones = 101000, Quemar = 1000, AnchoMuestreo = 10){
  
  Progress = txtProgressBar(min = 0, max = NoDeIteraciones,
                            style = 3)
  begin = sample(1:AnchoMuestreo, size = 1) 
  
  m = length(medias)
  
  muestras = matrix(NA, ncol = 2 * m + 4, 
                    nrow = (NoDeIteraciones - Quemar)/AnchoMuestreo)
  colnames(muestras) = c(paste('theta', names), paste('sigma', names),
                         'mu', 'tauSq', 'sigma 2', 'LogVerosimilitud')
  
  theta_j = medias
  sigma2_j = varianzas
  n = sum(n_j)
  mu = mean(theta_j)
  tauSq = var(theta_j)
  sigma2 = 1/mean(sigma2_j)
  contador = 0
  for (i in 1: NoDeIteraciones){
    # Actualizando el vector de medias
    theta_j = rnorm(n = m,
                    mean = ((n_j * medias)/sigma2_j + mu/tauSq) / 
                      (1/tauSq + n_j/sigma2_j),
                    sd = 1/sqrt((1/tauSq + n_j/sigma2_j)))
    
    # Actualizando el vector de varianzas
    sigma2_j = 1/rgamma(n = m, shape = 0.5* (nu + n_j),
                        rate = 0.5 * (nu * sigma2 + 
                                        ((n_j-1) * varianzas + 
                                           n_j*(medias - theta_j)^2)))
    
    # Actualizando mu
    mu = rnorm(n = 1,
               mean = (mu_0/gammaSq_0 + sum(theta_j)/tauSq) / (1/gammaSq_0 + m/tauSq),
               sd = 1/sqrt(1/gammaSq_0 + m/tauSq))
    
    # Actualizando tau^2
    tauSq = 1/rgamma(n = 1,
                   shape = 05 * (eta_0 + m),
                   rate = 0.5 * (eta_0 * tauSq_0 + sum((theta_j - mu)^2)))
    
    # Actualizando sigma^2
    sigma2 = rgamma(n = 1, 
                    shape = alpha_0 + (0.5 * m * nu),
                    rate = beta_0 + 0.5 * nu * sum(1/sigma2_j))
    
    if (i > Quemar & ((i - begin) %% AnchoMuestreo == 0)){
      contador = contador + 1
      muestras[contador,] = c(theta_j,sigma2_j, mu,
                              tauSq, sigma2,
                              sum(dnorm(y, 
                                        mean = rep(theta_j, n_j),
                                        sd = sqrt(rep(sigma2_j, n_j)),log = T)))
    }
    setTxtProgressBar(Progress, i)
  }
  close(Progress)
  beepr::beep()
  return(muestras)
}

## --- M4: Modelo Normal con medias y varianzas específicas por departamento y municipio ---

Gibbs4 = function(y,mean_k, mean_jk, var_k, var_jk,
                  # y = 'Todos los ptjs'. _jk = 'Valor mncpio'. _k = 'Valor dpto'
                  # n = 'No. Estd'. NMpioNdep = 'Municipios por depto'.
                  n_k,n_jk,NMnpioNDep,names_k, names_jk,
                  # Hiperparámetros del modelo
                  eps_0 = 1, kappaSq_0 = 2500,mu_0 = 250, gammaSq_0 = 2500, 
                  eta_0 = 1, tauSq_0 = 2500, nu_0 = 1, sigmaSq_0 = 2500,
                  # Otros argumentos importantes
                  NoDeIteraciones = 101000, Quemar = 1000, AnchoMuestreo = 10){
  # Barra de progreso
  Progress = txtProgressBar(min = 0, max = NoDeIteraciones,style = 3)
  
  # Valores auxiliares
  m = length(mean_k)
  j = length(mean_jk)
  n = length(y)
  
  # Almacenamiento del muestreador
  muestras = matrix(NA, ncol = m + j + 5, 
                    nrow = (NoDeIteraciones - Quemar)/AnchoMuestreo)
  colnames(muestras) = c(paste('M', names_jk), paste('D', names_k),
                         'mu', 'tauSq', 'kappaSq', 'sigma2', 'LogVerosimilitud')
  
  
  # Definición inicial de los parámetros
  zeta_jk = mean_jk; theta_k = mean_k; sigma2 = var(mean_jk); mu = mean(mean_k)
  tauSq = var(theta_k); kappaSq = var(y)
  
  # Muestreo sistemático 
  begin = sample(1:AnchoMuestreo, size = 1); contador = 0
  
  # Muestreador
  for (i in 1:NoDeIteraciones){
    # ACTUALIZANDO LAS MEDIAS POR MUNICIPIO [x]
    zeta_jk = rnorm(n = j,
                    mean = ((n_jk * mean_jk)/kappaSq + 
                              (rep(theta_k, NMnpioNDep)/sigma2)) / 
                      (n_jk/kappaSq + 1/sigma2),
                    sd = 1/sqrt(n_jk/kappaSq + 1/sigma2))
    
    SumaMediasMnpioPorDepto = tapply(zeta_jk, rep(names_k, NMnpioNDep), FUN = sum)
    
    # ACTUALIZANDO LAS MEDIAS POR DEPARTAMENTO [x]
    theta_k = rnorm(n = m, 
                    ## Revisar el cálculo de la suma de medias 
                    ## de municipios por departamento
                    mean = (SumaMediasMnpioPorDepto/sigma2 + mu/tauSq) / 
                      (NMnpioNDep/sigma2 + 1/tauSq),
                    ####
                    sd = 1/sqrt(NMnpioNDep/sigma2 + 1/tauSq))
    
    # ACTUALIZANDO SIGMA^2 [x]
    sigma2 = 1/rgamma(n = 1,
                      shape = 0.5 * (nu_0 + j),
                      rate  = 0.5 * (nu_0 * sigmaSq_0 + 
                                      sum((zeta_jk - rep(theta_k, NMnpioNDep))^2)))
    
    # ACTUALIZANDO KAPPA^2 [x]
    kappaSq = 1/rgamma(n = 1,
                       shape = 0.5 * (n + eps_0),
                       rate = 0.5 * (eps_0 * kappaSq_0 + 
                                       sum((y - rep(zeta_jk,n_jk))^2)))
    
    # ACTUALIZANDO MU [x]
    mu = rnorm(n = 1,
               mean = (sum(theta_k)/tauSq + mu_0/gammaSq_0) /
                 (m/tauSq + 1/gammaSq_0),
               sd = 1/sqrt(m/tauSq + 1/gammaSq_0))
    
    # ACTUALIZANDO TAU^2 [x]
    tauSq = 1/rgamma(n = 1,
                     shape = 0.5 * (eta_0 + m),
                     rate = 0.5 * (eta_0 * tauSq_0 +
                                     sum((theta_k - mu)^2)))
    
    # Condición ciclo de quemado y muestreo sistemático
    if (i > Quemar & ((i - begin) %% AnchoMuestreo == 0)){
      contador = contador + 1
      muestras[contador,] = c(zeta_jk,theta_k,mu, tauSq, kappaSq,sigma2,
                              sum(dnorm(y,mean = rep(zeta_jk,n_jk),
                                        sd = sqrt(kappaSq),log = T)))
    }
    setTxtProgressBar(Progress, i)
  }
  close(Progress);beepr::beep();return(muestras)
}

## --- M4: Modelo Normal con medias y varianzas específicas por departamento y municipio ---


Gibbs5 = function(y,mean_k, mean_jk, var_k, var_jk,
                  # y = 'Todos los ptjs'. _jk = 'Valor mncpio'. _k = 'Valor dpto'
                  # n = 'No. Estd'. NMpioNdep = 'Municipios por depto'.
                  n_k,n_jk,NMnpioNDep,names_k, names_jk,
                  # Hiperparámetros del modelo
                  eps_0 = 1, kappaSq_0 = 2500,mu_0 = 250, gammaSq_0 = 2500, 
                  eta_0 = 1, tauSq_0 = 2500, nu = 1, alpha_0 = 1, beta_0 = 1/2500,
                  # Otros argumentos importantes
                  NoDeIteraciones = 101000, Quemar = 1000, AnchoMuestreo = 10){
  # Barra de progreso
  Progress = txtProgressBar(min = 0, max = NoDeIteraciones,style = 3)
  
  # Valores auxiliares
  k = length(mean_k)
  jk = length(mean_jk)
  n = length(y)
  
  # Almacenamiento del muestreador
  muestras = matrix(NA, ncol = 2*k + jk + 5, 
                    nrow = (NoDeIteraciones - Quemar)/AnchoMuestreo)
  colnames(muestras) = c(paste('M', names_jk), paste('MD', names_k),paste('VD',names_k),
                         'mu', 'tauSq', 'kappaSq', 'sigma2', 'LogVerosimilitud')
  
  
  # Definición inicial de los parámetros
  zeta_jk = mean_jk; theta_k = mean_k; sigma2 = 1/mean(var_jk); mu = mean(mean_k)
  tauSq = var(theta_k); kappaSq = var(y); sigma2_k = var_k
  
  # Muestreo sistemático 
  begin = sample(1:AnchoMuestreo, size = 1); contador = 0
  
  for (i in 1:NoDeIteraciones){
    
    # ACTUALIZANDO zeta_jk
    zeta_jk = rnorm(n = jk,
                    mean = (n_jk * mean_jk/kappaSq + rep(theta_k/sigma2_k, NMnpioNDep)) / 
                      (n_jk/kappaSq + rep(1/sigma2_k, NMnpioNDep)), 
                    sd = 1/sqrt(n_jk/kappaSq + rep(1/sigma2_k, NMnpioNDep)))
    
    # ACTUALIZANDO theta_k
    SumaZeta = tapply(zeta_jk, rep(names_k, NMnpioNDep), sum)
    theta_k = rnorm(n = k, 
                    mean = (SumaZeta/sigma2_k + mu/tauSq) / 
                      (NMnpioNDep/sigma2_k + 1/tauSq),
                    sd = 1/sqrt(NMnpioNDep/sigma2_k + 1/tauSq))
    
    # ACTUALIZANDO sigma^2_k
    SumaZetaDifSq = (zeta_jk - rep(theta_k, NMnpioNDep))^2
    SumaZetaDifSq = tapply(SumaZetaDifSq, rep(names_k, NMnpioNDep), FUN = sum)
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # ¿Será este el problema?
    # SumaZetaDifSq = (NMnpioNDep - 1) * tapply(zeta_jk, rep(names_k, NMnpioNDep), var) + NMnpioNDep * (SumaZeta / NMnpioNDep - theta_k)^2
    # SumaZetaDifSq[is.na(SumaZetaDifSq)] = 0
    sigma2_k = 1/rgamma(n = k,
                        shape = 0.5 * (NMnpioNDep + nu),
                        rate  = 0.5 * (nu * sigma2 + SumaZetaDifSq))
    
    # ACTUALIZANDO sigma2
    sigma2 = rgamma(n = 1,
                    shape = 0.5 * (nu * k + alpha_0),
                    rate = 0.5 * (nu * sum(1/sigma2_k) + beta_0))
    
    # ACTUALIZANDO mu
    mu = rnorm(n = 1,
               mean = (sum(theta_k)/tauSq + mu_0/tauSq_0) / 
                 (k/tauSq + 1/gammaSq_0),
               sd = 1/sqrt(k/tauSq + 1/gammaSq_0))
    
    # ACTUALIZANDO tau^2
    tauSq = 1/rgamma(n = 1,
                     shape = 0.5 * (eta_0 + k),
                     rate = 0.5 * (sum((theta_k - mu)^2) + eta_0 * tauSq_0))
    
    # ACTUALIZANDO kappa^2
    kappaSq = 1/rgamma(n = 1,
                       shape = 0.5 * (n + eps_0),
                       rate = 0.5 * (sum((y - rep(zeta_jk, n_jk))^2) + eps_0 * kappaSq_0))
    
    # Condición ciclo de quemado y muestreo sistemático
    if (i > Quemar & ((i - begin) %% AnchoMuestreo == 0)){
      contador = contador + 1
      muestras[contador,] = c(zeta_jk,theta_k,sigma2_k,mu, tauSq, kappaSq,sigma2,
                              sum(dnorm(y,mean = rep(zeta_jk,n_jk),
                                        sd = sqrt(kappaSq),log = T)))
    }
    setTxtProgressBar(Progress, i)
  }
  close(Progress);beepr::beep();return(muestras)
}

