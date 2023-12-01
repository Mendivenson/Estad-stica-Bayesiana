# -------------------------------------------------------
# Caso de Estudio 2: Resultados. Estadística Bayesiana (2023-I)
# Prueba Saber 11 2022-2: Una perspectiva multinivel
# Fecha de Última Modificación: 03/Noviembre/2023
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchéz
# -------------------------------------------------------
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load('datos.RData')
# =============================== PUNTO 7 ======================================
# Usando M5 hacer un ranking de los departamentos basado en las medias específicas.
# Comparar este ranking con un ranking frecuentista.
# ==============================================================================
# Convenciones: 
#   - Rojo oscuro (Promedios significativamente menores al 250)
#   - Negro (Promedios que no difieren significativamente del 250)
#   - Verde Oscuro (Promedios significativamente superiores a 250)

library(data.table)
M5 = fread('Modelos/M5', drop = 1)

par(mfrow = c(1,2), mar = c(5,4,1,2))

# RANKING BAYESIANO:
RankingBayesiano = matrix(NA, ncol = 3, nrow = nrow(dataDepartamentos))
colnames(RankingBayesiano) = c('Inf', 'Medias', 'Sup')
rownames(RankingBayesiano) = colnames(M5)[grepl('MD',colnames(M5))]

## Llenado de la matriz 
RankingBayesiano[,"Medias"] = colMeans(M5[, .SD, .SDcols = rownames(RankingBayesiano)])
RankingBayesiano[,c("Inf",'Sup')] = t(apply(M5[, .SD, .SDcols = rownames(RankingBayesiano)],
                                          MARGIN = 2,
                                          FU = quantile, probs = c(0.025, 0.975)))
## Organizando la matriz
RankingBayesiano = RankingBayesiano[order(RankingBayesiano[,"Medias"]),]

nombres_organizados = departamentos[order(match(departamentos$CodDepartamento, 
                                  gsub('MD ','',rownames(RankingBayesiano)))),]$departamento
nombres_organizados[nombres_organizados == 'NORTE SANTANDER'] = 'N. SANTANDER'
## nombres_organizados[nombres_organizados == 'NORTE SANTANDER'] = 'N. SANTANDER'

m = nrow(RankingBayesiano)

plot(NA, NA, xlab = 'Media del puntaje global', ylab = "", main = "Ranking Bayesiano (Modelo 5)", 
     xlim = c(180,280), ylim = c(1,m),cex.axis = 0.75, yaxt = "n")
axis(side = 2, at = 1:m, labels = nombres_organizados, las = 2, cex.axis = 0.5)
abline(v = 250,  col = "gray", lwd = 3);abline(h = 1:m, col = "lightgray", lwd = 1)

colores = rep('black', m)
colores[which(RankingBayesiano[,"Sup"] < 250)] = 'darkred' 
colores[which(RankingBayesiano[,"Inf"] > 250) ] = 'darkgreen'

for (j in 1:m) {
  segments(x0 = RankingBayesiano[j,'Sup'], y0 = j, x1 = RankingBayesiano[j,'Inf'], y1 = j, col = colores[j])
  lines(x = RankingBayesiano[j,'Medias'], y = j, type = "p", pch = 16, cex = 0.6, col = colores[j])
}

# RANKING FRECUENTISTA:
dataDepartamentos = dataDepartamentos[order(dataDepartamentos$mean),]

plot(NA, NA, xlab = 'Media del puntaje global', ylab = "", main = "Ranking Frecuentista", 
     xlim = c(180,280), ylim = c(1,m),cex.axis = 0.75, yaxt = "n")
axis(side = 2, at = 1:m, labels = dataDepartamentos$departamento, las = 2, cex.axis = 0.5)
abline(v = 250,  col = "gray", lwd = 3);abline(h = 1:m, col = "lightgray", lwd = 1)

inf = dataDepartamentos$mean - 1.96 * 50/sqrt(dataDepartamentos$n)
sup = dataDepartamentos$mean + 1.96 * 50/sqrt(dataDepartamentos$n)

colores = rep('black', m)
colores[which(sup < 250)] = 'darkred' 
colores[which(inf > 250)] = 'darkgreen'

for (j in 1:m) {
  segments(x0 = sup[j], y0 = j, x1 = inf[j], y1 = j, col = colores[j])
  lines(x = dataDepartamentos$mean[j], y = j, type = "p", pch = 16, cex = 0.6, col = colores[j])
}

# =============================== PUNTO 8 ======================================
# Usando M5 hacer una segmentación con k means en cinco grupos de los departamentos
# respecto a sus medias específicas.
# ==============================================================================

# install.packages('corrplot')
library(corrplot)

# Calculando K medias
Progress = txtProgressBar(min = 1, max = 5000, style = 3)
matrizIncidencia = matrix(0, nrow = nrow(dataDepartamentos), ncol = nrow(dataDepartamentos))
for (i in 1:10000){
  a = kmeans(unlist(M5[i, 1113:1144]), 5)$cluster
  # Aquí se pasa del vector de clusters a la matriz de incidencia individual de cada
  # realización de medias de los departamentos.
  for (j in 1:5){
    b = (a == j)
    matrizIncidencia = matrizIncidencia + (b %*% t(b))
  }
  setTxtProgressBar(Progress,i)
}
close(Progress)
rownames(matrizIncidencia) = colnames(matrizIncidencia)
matrizIncidencia = matrizIncidencia/10000
matrizIncidencia = matrizIncidencia[rev(rownames(RankingBayesiano)),
                                    rev(rownames(RankingBayesiano))] 
rownames(matrizIncidencia) = colnames(matrizIncidencia) = rev(nombres_organizados)

corrplot::corrplot(corr = matrizIncidencia, 
                   is.corr = FALSE,
                   addgrid.col = NA, 
                   method = "color", 
                   tl.pos = "lt",
                   tl.cex = 0.8,
                   tl.col = 'black')
                   #mar = c(1,1,1,1),
                   #title = 'Matriz de incidencia por departamentos')
title(main = 'Matriz de incidencia por departamentos', line = 1.8, cex.main = 1.6)


### Mapa ----
RankingBayesiano = as.data.frame(RankingBayesiano)
colnames(RankingBayesiano) = c('Inf', 'mean', 'Sup')
RankingBayesiano = cbind(RankingBayesiano,'CodDepartamento' = gsub('MD ', '', rownames(RankingBayesiano)))
Grupos <- cbind(RankingBayesiano[,c("CodDepartamento","mean")] , kmeans(RankingBayesiano$mean, centers = 5, nstart = 25)$cluster)
names(Grupos) <- c("CodDepartamento","Media","Grupo")
dat_map = Grupos[,c("CodDepartamento","Grupo")]
colnames(dat_map) <- c("DPTO_CCDGO","Grupo")
dat_map$Grupo <- as.factor(dat_map$Grupo)

shp = sf::st_read("Mapas/MGN_DPTO_POLITICO.shp", quiet = T)
shp = shp[shp$DPTO_CCDGO != '88',]

colores <- c("#ff4848", "#ff8531", "#ffc131", "#8cff31", "#31aeff")

MapaMedia = left_join(x = shp, y = dat_map, by = c("DPTO_CCDGO")) %>% 
  select(DPTO_CCDGO, Grupo, geometry) %>%
  ggplot() +
  geom_sf(aes(fill = Grupo), size = 0.125, color = "#b2b2b2") +
  theme_bw() + 
  xlab("Longitud") + ylab("Latitud") +
  labs(title = "Conglomerados de departamentos por puntaje")+
  scale_fill_manual(values = colores, name = 'Grupo')

MapaMedia

# =============================== PUNTO 9 ======================================
# Calcular la media posterior y un intervalo de credibilidad al 95% del IPM para
# para los departamentos no medidos por el DANE.
# ==============================================================================

IPMreg = t(M5[,1113:1144])
rownames(IPMreg) = gsub('MD ','',rownames(IPMreg))
IPMreg = cbind('IPM' = Pobreza$`2018`, IPMreg)

IPMpred = IPMreg[is.na(IPMreg[,'IPM']),]
IPMreg = IPMreg[!is.na(IPMreg[,'IPM']),]

IPMpredicciones = matrix(NA, ncol = nrow(IPMpred), nrow = nrow(M5))
colnames(IPMpredicciones) = rownames(IPMpred)
Progress = txtProgressBar(min = 2, max = nrow(M5)+1, style = 3)
for (i in 2:(nrow(M5) + 1)){
  coeficientes = lm(IPMreg[,'IPM'] ~ IPMreg[,i])$coefficients
  IPMpredicciones[i-1,] = coeficientes[1] + coeficientes[2] * IPMpred[,i]
  setTxtProgressBar(Progress,i)
}
close(Progress)

ResultadosIPM = cbind('CodDepartamento' = colnames(IPMpredicciones),
                      'Departamento' = departamentos$departamento[departamentos$CodDepartamento %in% colnames(IPMpredicciones)],
                      'Media posterior' = colMeans(IPMpredicciones),
                      '2.5%' = apply(IPMpredicciones, 2, quantile, probs = 0.025),
                      '97.5%' = apply(IPMpredicciones, 2, quantile, probs = 0.975))

ResultadosIPM

# library(knitr)
# kable(ResultadosIPM, format = 'latex')

# =============================== PUNTO 10 =====================================
# Usando M5 hacer el ranking de los municipios y un agrupamiento de k means en 
# ocho grupos, presentar la matriz de incidencia organizada respecto al ranking
# bayesiano previo.
# ==============================================================================
par(mfrow = c(1,1))

RankingBayesiano = M5[,1:1112]
n = nrow(RankingBayesiano)

RankingBayesiano = rbind(RankingBayesiano, 'mean' = t(colMeans(RankingBayesiano)))
RankingBayesiano = rbind(RankingBayesiano, 
                         '2.5%' = t(apply(X = RankingBayesiano[1:n,],
                                          2,quantile, 0.25)),
                         '97.5%' = t(apply(X = RankingBayesiano[1:n,],
                                           2,quantile, 0.975)))
RankingBayesiano = as.data.frame(t(RankingBayesiano[(n+1):(n+3),]))
ncol(RankingBayesiano)
colnames(RankingBayesiano) <- c("mean","2.5%","97.5%")
RankingBayesiano = RankingBayesiano[order(RankingBayesiano[,'mean']),]
RankingBayesiano

## Matriz de incidencia ----
# install.packages('corrplot')
library(corrplot)

# Calculando K medias
Progress = txtProgressBar(min = 1, max = 5000, style = 3)
matriz = matrix(0, nrow = nrow(dataMunicipios), ncol = nrow(dataMunicipios))
for (i in 1:10000){
  a = kmeans(unlist(M5[i, 1:1112]), 8)$cluster
  # Aquí se pasa del vector de clusters a la matriz de incidencia individual de cada
  # realización de medias de los departamentos.
  for (j in 1:8){
    b = (a == j)
    matriz = matriz + (b %*% t(b))
  }
  setTxtProgressBar(Progress,i)
}
close(Progress)
rownames(matriz) = colnames(matriz)
matriz = matriz/10000
corrplot::corrplot(corr = matriz[rev(rownames(RankingBayesiano)),rev(rownames(RankingBayesiano))], 
                   is.corr = FALSE,
                   addgrid.col = NA, 
                   method = "color", 
                   tl.pos = "n",
                   title = 'Matriz de incidencia entre municipios')

## Mapa Colombia ----
library(tibble)
RankingBayesiano <- tibble::rownames_to_column(RankingBayesiano, var = "CodMunicipio")
RankingBayesiano$NewCode <- as.numeric(gsub("[^0-9]", "", RankingBayesiano$CodMunicipio))
RankingBayesiano <- RankingBayesiano[,c("NewCode","mean")]
RankingBayesiano

Grupos <- cbind(RankingBayesiano[,c("NewCode","mean")] , kmeans(RankingBayesiano$mean, centers = 8, nstart = 25)$cluster)
names(Grupos) <- c("CodMunicipio","Media","Grupo")
dat_map = Grupos[,c("CodMunicipio","Grupo")]
colnames(dat_map) <- c("CodMunicipio","Grupo")
dat_map$Grupo <- as.factor(dat_map$Grupo)


shp = sf::st_read("Mapas/MGN_MPIO_POLITICO.shp", quiet = T)

# Quitando a San Andŕes del mapa:
shp = shp[shp$DPTO_CCDGO != '88',]
shp$CodMunicipio = as.numeric(paste(shp$DPTO_CCDGO, shp$MPIO_CCDGO, sep = ''))

colores <- c("#ff4848", "#ff8531", "#ffc131", "#8cff31", "#31aeff","#33f08f","#ffe90a","#f577ff") 

MapaMediaMun = left_join(x = shp, y = dat_map, by = 'CodMunicipio') %>% 
  select(CodMunicipio, Grupo, geometry) %>%
  ggplot() +
  geom_sf(aes(fill = Grupo), size = 0.025, color = "#4F4F4F") +
  theme_bw() + 
  xlab("Longitud") + ylab("Latitud") +
  labs(title = "Conglomerado de puntajes por municipios")+
  scale_fill_manual(values = colores, name = 'Grupo')

MapaMediaMun


# =============================== PUNTO 11 =====================================
# Calcular la medias posterior y el intervalo de credibilidad del 95% de la cobertura
# neta secundaria en 2022 para todos los municipios no medidos por el MEN
# ==============================================================================

EDUreg = t(M5[,1:1112])
rownames(EDUreg) = gsub('M ','',rownames(EDUreg))
educacion = right_join(educacion,dataMunicipios)[,c('CodMunicipio', 'cobertura')]
educacion = educacion[order(match(educacion$CodMunicipio,rownames(EDUreg))),]
EDUreg = cbind('cobertura' = educacion$cobertura, EDUreg)

EDUpred = EDUreg[is.na(EDUreg[,'cobertura']),]
EDUreg = EDUreg[!is.na(EDUreg[,'cobertura']),]

EDUpredicciones = matrix(NA, ncol = nrow(EDUpred), nrow = nrow(M5))
colnames(EDUpredicciones) = rownames(EDUpred)
Progress = txtProgressBar(min = 2, max = nrow(M5)+1, style = 3)
for (i in 2:(nrow(M5) + 1)){
  coeficientes = lm(EDUreg[,'cobertura'] ~ EDUreg[,i])$coefficients
  EDUpredicciones[i-1,] = coeficientes[1] + coeficientes[2] * EDUpred[,i]
  setTxtProgressBar(Progress,i)
}
close(Progress)

ResultadosCobertura = cbind('CodMunicipio' = colnames(EDUpredicciones),
                            'Municipio' = municipios$municipio[municipios$CodMunicipio %in% colnames(EDUpredicciones)],
                            'Media posterior' = colMeans(EDUpredicciones),
                            '2.5%' = apply(EDUpredicciones, 2, quantile, probs = 0.025),
                            '97.5%' = apply(EDUpredicciones, 2, quantile, probs = 0.975))

ResultadosCobertura

# kable(ResultadosCobertura, format = 'latex')

# =============================== PUNTO 12 =====================================
# Validar la bondad de ajuste de M5 para cada municipio, usando como estadísticos
# de prueba:
#   - Mínimo        - Media       - Desviación estándar
#   - Máximo        - Mediana     - Rango intercuartílico
# ==============================================================================

library(progress)



pppValues = matrix(0, nrow = nrow(dataMunicipios), ncol = 6)
colnames(pppValues) = c('Mín', 'Media', 'Mediana', 'Máx', 'Desv. Est.', 'IQR')

desv = sqrt(M5$kappaSq)

Progress = txtProgressBar(min = 0, max = 10000, style = 3)
tictoc::tic()
for (i in 1:10000){
  desvi = desv[i]
  for (j in 1:nrow(dataMunicipios)){
    muestra = rnorm(n = dataMunicipios$n[j],
                    mean = M5[[i, paste('M', dataMunicipios$CodMunicipio[j])]],
                    sd = desvi)
    pppValues[j,"Mín"] =  pppValues[j,"Mín"] + (min(muestra) < dataMunicipios$min[j])
    pppValues[j,"Media"] =  pppValues[j,"Media"] + (mean(muestra) < dataMunicipios$mean[j])
    pppValues[j,"Mediana"] =  pppValues[j,"Mediana"] + (median(muestra) < dataMunicipios$mediana[j])
    pppValues[j,"Máx"] =  pppValues[j,"Máx"] + (max(muestra) < dataMunicipios$max[j])
    pppValues[j,"Desv. Est."] =  pppValues[j,"Desv. Est."] + (sd(muestra) < sqrt(dataMunicipios$var[j]))
    pppValues[j,"IQR"] =  pppValues[j,"IQR"] + (IQR(muestra) < dataMunicipios$IQR[j])
    setTxtProgressBar(Progress,i)
  }
}
close(Progress)
tictoc::toc()
beepr::beep()

pppValues = pppValues/10000
boxplot(pppValues, main = 'Bondad de ajuste', sub = 'Valores ppp para el modelo 5 por departamentos', pch = '*')
