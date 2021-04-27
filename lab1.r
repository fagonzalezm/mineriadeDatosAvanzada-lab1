require(mclust)
require(reshape2)
require(ggplot2)

######### Se carga la base de datos #########
names <- c("id", "clumpThickness", "uniformityOfCellSize", "uniformityOfCellShape",
           "marginalAdhesion", "singleEpithelialCellSize", "bareNuclei",
           "blandChromatin", "normalNucleoli", "mitoses", "class")
data <- read.table("./breast-cancer-wisconsin.data", sep=",", col.names = names)

set.seed(20)

################################
####### Preprocesamiento #######
################################


####### Missing Values ###########
# - Como existen 16 datos que presentan missing values en la variable barnuclei y el total de datos es 699,
#   se opta por eliminar estos datos.
data.original <- data

data.n = nrow(data)
data.m = ncol(data)

for (row in 1:data.n) {
  for (col in 1:data.m) {
    if (data[row, col] == "?") {
      data[row, col] <- NA
    }
  }
}
data$bareNuclei <- as.integer(data$bareNuclei)
data <- na.omit(data)

table(data$class)
table(data.original$class)

# - Se cuenta con 683 datos: 444 células benignas y 239 cancerosas.
# - Se descartaron 14 células benignas y 2 cancerosas.

####### Elección de variables de interés #######
# - La variable id no proporciona información útil para el estudio, así que no se considera
# - La variable class se usa como gold standard, por lo que su uso se limita a evaluar el modelo,
#   mas no a su confección.
# - Se estudia la distribución del resto de las variables para decidir si serán consideradas en
#   el modelo

# Se quita la variable id
features <- data[,2:11]

# Se transforman los 2 y 4 de la variable class a "B" y "M" según corresponda
features$class <- as.character(features$class)
features$class[features$class == "2"] <- "B"
features$class[features$class == "4"] <- "M"

# Se grafica un boxplot de las demas variables respecto a class
features.melt <- melt(features, id.var = "class")
p <- ggplot(data = features.melt, aes(x=variable, y=value))
p <- p +  geom_boxplot(aes(color=class))
p <- p + geom_jitter(aes(color=class), alpha=0.2, size=0.5)
p <- p + facet_wrap( ~ variable, scales="free")
p

# Se observa que para todas las variables existe una relación entre sus valores
# y si la célula corresponde a una cancerosa o no, excepto con mitoses.

#summary(features)
#var(features)
#by(features, features$class, summary)

####### Test de normalidad #######
# H0: La muestra sigue una distribución normal
# H1: La muestra no sigue una distribución normal

shapiro.test(features$clumpThickness)
shapiro.test(features$uniformityOfCellSize)
shapiro.test(features$uniformityOfCellShape)
shapiro.test(features$marginalAdhesion)
shapiro.test(features$singleEpithelialCellSize)
shapiro.test(features$bareNuclei)
shapiro.test(features$blandChromatin)
shapiro.test(features$normalNucleoli)
shapiro.test(features$mitoses)

#  - Todas con un p-value < 2.2e-16
#  - Por lo tanto se rechaza la hipótesis nula.
#  - Existe evidencia estadística suficiente para rechazar H0, entonces.
#    se puede afirmar que las varaibles no siguen una distribución normal.

# Esto se puede confirmar revisando sus respectivos histogramas

p <- ggplot(data = features.melt, aes(x=value, fill=variable))
p <- p + geom_bar()
p <- p + facet_wrap( ~ variable, scales="free")
p


features <- data[,1:9]
class <- data$class


################mclust######################

#primer intento

mod1 = Mclust(features) #DEFAULT
summary(mod1)

#segundo intento
mod2 = Mclust(features, G = 3)  #Numero de grupos = 3.
summary(mod2, parameters = TRUE)

#tercer intento
#mod6 = Mclust(features, prior = priorControl(functionName="defaultPrior", shrinkage=0.1), modelNames ="EII")  
#"EII" = spherical, equal volume # Using prior #The function priorControl is used to specify a conjugate prior for EM within MCLUST. 
#summary(mod6,parameter = TRUE)
#plot(mod6, what = "classification")

#BIC
BIC<-mclustBIC(features, prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
plot(BIC)  #se grafican los BIC por configuraci?n de par?metros
summary(BIC)  # se presentan los mejores valores BIC

#Usando VVV,3
mod11=Mclust(features,x=BIC) # en base al mejor valor BIC se realiza el mclust
summary(mod11)#se muestra resultado y tabla de clustering

#Graficando
plot(mod11, what = "classification")  #se grafica la configuraci?n de agrupamientos.
legend("bottomright", legend = 1:3,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(class, mod11$classification) #distribuci?n de clases por cada grupo.


#Usando otros valores
mod12 = Mclust(features, G=7, prior = priorControl(functionName="defaultPrior", shrinkage=0.1), modelNames ="VVI")  
plot(mod12, what = "classification")  #se grafica la configuraci?n de agrupamientos.
legend("bottomright", legend = 1:7,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")
table(class, mod12$classification) #distribuci?n de clases por cada grupo.