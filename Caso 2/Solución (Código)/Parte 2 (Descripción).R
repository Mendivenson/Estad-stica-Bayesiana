# -------------------------------------------------------
# Caso de Estudio 2: Descripción. Estadística Bayesiana (2023-I)
# Prueba Saber 11 2022-2: Una perspectiva multinivel
# Fecha de Última Modificación: 12/Octubre/2023
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchéz
# -------------------------------------------------------

# En este script se presenta la descripción de los datos y de los métodos
# bayesianos utilizados
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("datos.RData")

# =============================== PUNTO 1 ======================================
# En un gráficos con dos paneles hacer un mapa de Colombia por departamentos donde
# se muestren los valores de la media muestral del puntaje global y la incidencia
# de la pobreza monetaria (Del año 2018).
# ==============================================================================

# install.packages('ggplot2')
# install.packages('sf')
# install.packages('dplyr')

library(ggplot2)
library(dplyr)


shp = sf::st_read("Mapas/MGN_DPTO_POLITICO.shp", quiet = T)

dat_map = dataDepartamentos[,c("CodDepartamento","mean")]
colnames(dat_map) <- c("DPTO_CCDGO","Media")

MapaMedia = left_join(x = shp, y = dat_map, by = c("DPTO_CCDGO")) %>% 
  select(DPTO_CCDGO, Media, geometry) %>%
  ggplot() +
  geom_sf(aes(fill = Media), size = 0.125, color = "#b2b2b2") +
  theme_bw() + 
  xlab("Longitud") + ylab("Latitud") +
  labs(title = "Media del puntaje global por departamentos")+
  scale_fill_gradient(name = 'Media',low = '#FF4500', high = "#FFB06A")

dat_map = Pobreza[,c('CodDepartamento', '2018')]
colnames(dat_map) = c("DPTO_CCDGO","Incidencia de la pobreza monetaria (2018)")

MapaIncidencia = inner_join(x = shp, y = dat_map, by = c("DPTO_CCDGO")) %>% 
  select(DPTO_CCDGO, `Incidencia de la pobreza monetaria (2018)`, geometry) %>%
  ggplot() +
  geom_sf(aes(fill = `Incidencia de la pobreza monetaria (2018)`), 
          size = 0.125, color = "#b2b2b2") +
  theme_bw() + 
  xlab("Longitud") + ylab("Latitud") +
  labs(title = "Incidencia de la pobreza monetaria (2018)") +
  scale_fill_gradient(name = 'Incidencia',low = '#FFB06A', high = "#FF4500") 

# install.packages('gridExtra')
library(gridExtra)
grid.arrange(MapaMedia, MapaIncidencia, ncol = 2)

# =============================== PUNTO 2 ======================================
# En un gráficos con dos paneles hacer un mapa de Colombia por municipios donde
# se muestren los valores de la media muestral del puntaje global y la cobertura
# neta secundaria en 2022.
# ==============================================================================

shp = sf::st_read("Mapas/MGN_MPIO_POLITICO.shp", quiet = T)

# Quitando a San Andŕes del mapa:
shp = shp[shp$DPTO_CCDGO != '88',]

dat_map = dataMunicipios[,c('CodMunicipio', 'mean')]
shp$CodMunicipio = as.numeric(paste(shp$DPTO_CCDGO, shp$MPIO_CCDGO, sep = ''))

MapaMediaMun = left_join(x = shp, y = dat_map, by = 'CodMunicipio') %>% 
  select(CodMunicipio, mean, geometry) %>%
  ggplot() +
  geom_sf(aes(fill = mean), size = 0.025, color = "#4F4F4F") +
  theme_bw() + 
  xlab("Longitud") + ylab("Latitud") +
  labs(title = "Media del puntaje global por municipios")+
  scale_fill_gradient(name = 'Media',low = '#000080', high = "#87CEEB")


MapaEducacion = left_join(x = shp, y = educacion, by = 'CodMunicipio') %>% 
  select(CodMunicipio, cobertura, geometry) %>%
  ggplot() +
  geom_sf(aes(fill = cobertura), size = 0.025, color = "#4F4F4F") +
  theme_bw() + 
  xlab("Longitud") + ylab("Latitud") +
  labs(title = "Cobertura neta secundaria (2022)")+
  scale_fill_gradient(name = paste('Cobertura','neta', sep = '\n')
                      ,low = '#000080', high = "#87CEEB")

grid.arrange(MapaMediaMun, MapaEducacion, ncol = 2)

# =============================== PUNTO 3 ======================================
# En un gráfico con cuatro paneles hacer el DAG de los modelos M2, M3, M4 y M5. 
# ==============================================================================

# install.packages('latex2exp')
# install.packages ('igraph')
library(igraph)
library(latex2exp)

# ================================ MODELO 2 ====================================
Modelo2 = graph(edges = c(1, 2, 
                          1, 3, 
                          2, 4,
                          2, 5,
                          3, 6,
                          3, 7,
                          4, 8,
                          4, 9,
                          5, 10,
                          5, 11), directed = TRUE)

# Atributos de los nodos:
V(Modelo2)$label = c(TeX("$y_{i,j}$"), 
                      TeX("$\\theta_{j}$"), TeX("$\\sigma^2$"),
                      TeX('$\\tau^2$'), TeX('$\\mu$'), TeX('$\\nu_0$'), TeX('$\\sigma^2_0$'),
                      TeX('$\\eta_0$$'), TeX('$\\tau^2_0$'), TeX('$\\mu_0$'), TeX('$\\gamma^2_0$'))

V(Modelo2)$label.color = "black"; V(Modelo2)$color = "white"

V(Modelo2)[1:5]$shape = 'circle'; V(Modelo2)[6:11]$shape = "square"

V(Modelo2)$size = 30

# Atributos de las flechas del DAG
E(Modelo2)$arrow.size = 0.05;E(Modelo2)$color = 'black'

# ================================ MODELO 3 ====================================

Modelo3 =  graph(edges = c(1, 2,                     # Primer nivel 
                        1, 3, 
                        2, 4,                        # Segundo nivel
                        2, 5,
                        3, 6,
                        3, 7,
                        4, 8,                        # Tercer nivel
                        4, 9,
                        5, 10,
                        5, 11, 
                        6, 12, 
                        6, 13), directed = TRUE)

# Atributos de los nodos:
V(Modelo3)$label = c(TeX("$y_{i,j}$"), 
                     TeX("$\\theta_{j}$"), TeX("$\\sigma^2_j$"),
                     TeX('$\\tau^2$'), TeX('$\\mu$'), TeX('$\\sigma^2$'), TeX('$\\nu$'),
                     TeX('$\\eta_0$$'), TeX('$\\tau^2_0$'), TeX('$\\gamma^2_0$'),TeX('$\\mu_0$'),
                     TeX('$\\alpha_0$'), TeX('$\\beta_0$') )

V(Modelo3)$label.color = "black"; V(Modelo3)$color = "white"

#V(Modelo3)[1:12]$shape = 'square';# V(Modelo2)[8:11]$shape = "square"
V(Modelo3)[1:6]$shape = 'circle'; V(Modelo3)[7:13]$shape = 'square'
V(Modelo3)$size = 30

# Atributos de las flechas del DAG
E(Modelo3)$arrow.size = 0.05;E(Modelo3)$color = 'black'

# ================================ MODELO 4 ====================================

Modelo4 = graph(edges = c(1, 2, 
                          1, 3, 
                          2, 4,
                          2, 5,
                          3, 6,
                          3, 7,
                          4, 8,
                          4, 9,
                          5, 10,
                          5, 11,
                          8, 12,
                          8, 13,
                          9, 14,
                          9, 15), directed = TRUE)

# Atributos de los nodos:
V(Modelo4)$label = c(TeX("$y_{_{i, j, k}}$"), 
                     TeX("$\\zeta_{j, k}$"), TeX("$\\kappa^2$"),
                     TeX('$\\theta_k$'), TeX('$\\sigma^2$'), TeX('$\\epsilon_0$'), TeX('$\\kappa^2_0$'),
                     TeX('$\\mu$$'), TeX('$\\tau^2$'), TeX('$\\nu_0$'), TeX('$\\sigma^2_0$'),
                     TeX('$\\mu_0$'), TeX('$\\gamma^2_0$'), TeX('$\\eta_0$'), TeX('$\\Tau^2_0$')  )

V(Modelo4)$label.color = "black"; V(Modelo4)$color = "white"

V(Modelo4)[1:9]$shape = 'circle'; V(Modelo4)[10:15]$shape = "square"; V(Modelo4)[6:7]$shape = 'square'

V(Modelo4)$size = 30
V(Modelo4)$size[1:2] = 40

# Atributos de las flechas del DAG
E(Modelo4)$arrow.size = 0.05;E(Modelo4)$color = 'black'

# ================================ MODELO 4 ====================================

Modelo5 = graph(edges = c(1, 2, 
                          1, 3, 
                          2, 4,
                          2, 5,
                          3, 6,
                          3, 7,
                          4, 8,
                          4, 9,
                          5, 10,
                          5, 11,
                          8, 12,
                          8, 13,
                          9, 14,
                          9, 15,
                          11, 16,
                          11, 17), directed = TRUE)

# Atributos de los nodos:
V(Modelo5)$label = c(TeX("$y_{_{i, j, k}}$"), 
                     TeX("$\\zeta_{j, k}$"), TeX("$\\kappa^2$"),
                     TeX('$\\theta_k$'), TeX('$\\sigma^2_k$'), TeX('$\\epsilon_0$'), TeX('$\\kappa^2_0$'),
                     TeX('$\\mu$$'), TeX('$\\tau^2$'), TeX('$\\nu$'), TeX('$\\sigma^2$'),
                     TeX('$\\mu_0$'), TeX('$\\gamma^2_0$'), TeX('$\\eta_0$'), TeX('$\\Tau^2_0$'), 
                     TeX('$\\alpha_0$'),TeX('$\\beta_0$') )

V(Modelo5)$label.color = "black"; V(Modelo5)$color = "white"

V(Modelo5)[1:11]$shape = 'circle'; V(Modelo5)[12:17]$shape = "square"; V(Modelo5)[6:7]$shape = 'square'
V(Modelo5)[10]$shape = 'square'

V(Modelo5)$size = 30
V(Modelo5)$size[1:2] = 40
# Atributos de las flechas del DAG
E(Modelo5)$arrow.size = 0.05;E(Modelo5)$color = 'black'

# =============================== LA GRILLA ====================================

# https://bugs.launchpad.net/igraph/+bug/1039350 (Cambiar el sentido de las flechas)
par(mfrow = c(2,2), tcl  = 0.5, 
    mar = c(1, 1, 2, 1), mgp = c(1.5, 0.5, 0), oma = c(4, 10, 0.2, 10))


plot(Modelo2, layout = layout.reingold.tilford,
     vertex.label.cex = 1.5, ylim = c(1.5, -1.5),
     edge.arrow.mode = T)
title("Modelo 2", line=-1.2)
mtext(paste('a) Modelo normal con medias específicas', 
            'por departamento', sep ='\n'), 
      side = 1, line = -0.5, at = 0.1, cex = 0.8)

plot(Modelo3, layout = layout.reingold.tilford,
     vertex.label.cex = 1.5, ylim = c(1.5, -1.5),
     edge.arrow.mode = T)
title("Modelo 3", line=-1.2)
mtext(paste('b) Modelo normal con medias y varianzas específicas',
            'por departamento', sep = '\n'), 
      side = 1, line = -0.5, at = 0.1, cex = 0.8)

plot(Modelo4, layout = layout.reingold.tilford,
     vertex.label.cex = 1.5, ylim = c(1.5, -1.5),
     edge.arrow.mode = T)
title("Modelo 4", line= -1.2)
mtext(paste('c) Modelo normal con medias específicas',
            'por municipio y departamento', sep = '\n'), 
      side = 1, line = -0.5, at = 0.1, cex = 0.8)

plot(Modelo5, layout = layout.reingold.tilford,
     vertex.label.cex = 1.5, ylim = c(1.5, -1.5),
     edge.arrow.mode = T)
title("Modelo 5", line=-1.2)
mtext(paste('d) Modelo normal con medias específicas'
            ,'por municipio y departamento', sep =  '\n'), 
      side = 1, line = -0.5, at = 0.1, cex = 0.8)

par(cex.main = 1.8)
title(main = "Gráficos acíclicos dirigidos", outer = TRUE, line = -1)

# Al hacer zoom al gráfico este se ve mejor.

# Correr la siguiente línea para volver a los valores predeterminados de par
# graphics.off()


