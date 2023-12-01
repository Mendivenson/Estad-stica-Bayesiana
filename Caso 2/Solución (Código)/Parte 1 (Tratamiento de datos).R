# -------------------------------------------------------
# Caso de Estudio 2: Tratamiento de datos. Estadística Bayesiana (2023-I)
# Prueba Saber 11 2022-2: Una perspectiva multinivel
# Fecha de Última Modificación: 12/Octubre/2023
# Autores: 
#         - Michel Mendivenson Barragán Zabala
#         - Gerardo Sebastián Gil Sanchéz
# -------------------------------------------------------

# **Nota:** Para correr todos los script basta con ponerlos en la misma carpeta
# y correr primero esta parte, de lo contrario en los demás scripts habrán
# problemas con encontrar los datos en los demás scripts.

# ------------------ TRATAMIENTO PREVIO DE LA INFORMACIÓN ----------------------
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Se define un 'diccionario' con el código de departamentos DANE para posterior presentación:
departamentos = cbind('CodDepartamento' = c('05', '08', '11', '13', '15', '17', '18', '19', '20',
                              '23', '25', '27', '41', '44', '47', '50', '52', '54', 
                              '63', '66', '68', '70', '73', '76', '81', '85','86',
                              '91', '94', '95', '97', '99'),
                      'departamento' = c("ANTIOQUIA", "ATLANTICO", "BOGOTA", "BOLIVAR",
                                      "BOYACA", "CALDAS", "CAQUETA", "CAUCA", "CESAR",
                                      "CORDOBA", "CUNDINAMARCA", "CHOCO", "HUILA", 
                                      "LA GUAJIRA", "MAGDALENA", "META", "NARIÑO", 
                                      "NORTE SANTANDER", "QUINDIO", "RISARALDA", "SANTANDER", 
                                      "SUCRE", "TOLIMA", "VALLE", "ARAUCA", 
                                      "CASANARE", "PUTUMAYO", "AMAZONAS",
                                      "GUAINIA", "GUAVIARE", "VAUPES", "VICHADA"))
departamentos = as.data.frame(departamentos)

# LEER EL ARCHIVO SEPARADO POR ;

data = read.csv('Datos/SB11_20222.TXT', sep = ';')


# QUITANDO LOS DATOS DESCRITOS EN EL CASO                                       # Únicamente los estudiantes:
data = subset(data, ESTU_NACIONALIDAD == 'COLOMBIA')                            # - Con nacionalidad colombiana
data = subset(data,ESTU_PAIS_RESIDE == 'COLOMBIA')                              # - Que reside en Colombia
data = subset(data, ESTU_ESTADOINVESTIGACION == 'PUBLICAR')                     # - Resultados para publicación
data = subset(data, COLE_DEPTO_UBICACION != 'SAN ANDRES')                       # - Colegios que no están en San Andrés
# colnames((head(FilteredData[, grepl("PUNT", names(FilteredData))])))

# Seleccionando las columnas de intéres de la base de datos:
data = data[, c('COLE_DEPTO_UBICACION', 'COLE_COD_DEPTO_UBICACION', 
                                'COLE_MCPIO_UBICACION','COLE_COD_MCPIO_UBICACION',
                                'PUNT_GLOBAL')]
colnames(data) = c('departamento','CodDepartamento','municipio','CodMunicipio','puntaje')

# Quitando los datos sin datos en algunas de las columnas seleccionadas
data = data[complete.cases(data),]                                              # - Sin datos faltantes en ninguna de las columnas seleccionadas.

if (nrow(data) == 525061){
  print('Los datos fueron filtrados correctamente')
} else {
  print('Existe algún error en el filtrado')
}

# Corrección código 5 a 05 y 8 a 08
data$CodDepartamento = as.character(data$CodDepartamento)
data$CodDepartamento[data$CodDepartamento == '5'] = '05'
data$CodDepartamento[data$CodDepartamento == '8'] = '08'

# CALCULANDO LAS ESTADÍSTICAS SUFICIENTES POR DEPARTAMENTOS Y POR MUNICIPIOS:
# install.packages('dplyr')
library(dplyr)
dataDepartamentos = data %>%
  group_by(CodDepartamento) %>%
  summarise(CodDepartamento = unique(CodDepartamento),
            departamento = unique(departamento),
            n = n(),
            mean = mean(puntaje),
            var = var(puntaje),
            nMunicipios = n_distinct(CodMunicipio),)
head(dataDepartamentos)

dataMunicipios = data %>%
  group_by(CodMunicipio) %>%
  summarise(CodMunicipio = unique(CodMunicipio),
            municipio = unique(municipio),
            n = n(),
            mean = mean(puntaje),
            var = var(puntaje),
            min = min(puntaje),
            max = max(puntaje),
            IQR = IQR(puntaje),
            mediana = median(puntaje))
head(dataMunicipios)

puntajes = data[,c('CodDepartamento', 'CodMunicipio','puntaje')]
puntajes = puntajes[order(puntajes$CodMunicipio, 
                          decreasing = F),]
head(puntajes)
rm(data)

# AGREGANDO LOS DATOS DE LA COBERTURA NETA SECUNDARIA:
educacion = read.csv(file = 'Datos/estadísticas educación.csv')
educacion = educacion[educacion$AÑO == 2022,]
educacion = educacion[,c('CÓDIGO_MUNICIPIO', 'COBERTURA_NETA_SECUNDARIA')]
colnames(educacion) = c('CodMunicipio', 'cobertura')
head(educacion)

# AGREGANDO LOS DATOS DE LA INCIDENCIA DE LA POBREZA MONETARIA:
# install.packages('readxl')
Pobreza = readxl::read_xls(path = 'Datos/pobreza monetaria.xls',
                           sheet = 'Pobreza Monetaria (%)',range = 'A16:P40', col_names = T)
names(Pobreza)[1] = 'departamento'
Pobreza = Pobreza[,c('departamento','2018')]

# Para poder poner códigos a la base de datos:
Pobreza$departamento = toupper(iconv(Pobreza$departamento, to = "ASCII//TRANSLIT"))
Pobreza$departamento = gsub("\\bde \\b", "", Pobreza$departamento, ignore.case = TRUE)
Pobreza$departamento = gsub("\\b del cauca\\b", "", Pobreza$departamento, ignore.case = TRUE)
Pobreza$departamento =gsub("\\bnarino\\b", "NARIÑO", Pobreza$departamento, ignore.case = TRUE)
Pobreza$departamento = gsub("\\b d.c.", "", Pobreza$departamento, ignore.case = TRUE)

Pobreza = left_join(x = departamentos, y = Pobreza)


# Diccionario de municipios:
municipios = dataMunicipios[c('CodMunicipio', 'municipio')]
dataMunicipios = dataMunicipios[, -which(names(dataMunicipios) == "municipio")]

save.image("datos.RData")
#load("datos.RData")
