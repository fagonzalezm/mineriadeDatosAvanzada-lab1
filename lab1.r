require(mclust)
require(reshape2)
require(ggplot2)
require(corrplot)

library(cluster)
library(factoextra)

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

# Estadistica descriptiva
# meadias
means <- sapply(features[, 1:9], mean)

# medianas
medians <- sapply(features[, 1:9], median)

#varianzas
vars <- sapply(features[, 1:9], var)

#coeficientes de variacion
coefs <- sqrt(vars)/means

means
medians
vars
coefs


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

#Ninguna sigue una distribuici�n normal, se opta por un test no param�trico
data$bareNuclei <- as.numeric(data$bareNuclei)
data.cor = cor(data, method = 'spearman')

corrplot(data.cor, method = 'ellipse')
title("Matriz de Correlacion")


features <- data[,2:10]
class <- data$class


################mclust######################


#Quitando la variable shape.
features <- data[,2:10]
features <- features[, -3]
class <- data$class

#BIC
BIC<-mclustBIC(features, prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
plot(BIC)  #se grafican los BIC por configuraci?n de par?metros
summary(BIC)  # se presentan los mejores valores BIC

#Best BIC values:
#  VVI,9       VVI,8       VVI,7
#BIC      -12747.68 -13390.6688 -13409.6683
#BIC diff      0.00   -642.9912   -661.9907

#Quitando la variable size.
features <- data[,2:10]
features <- features[, -2]
class <- data$class

#BIC
BIC2<-mclustBIC(features, prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
plot(BIC2)  #se grafican los BIC por configuraci?n de par?metros
summary(BIC2)  # se presentan los mejores valores BIC

#Best BIC values:
#  VVI,8        VVI,9       VVV,5
#BIC      -13382.17 -13425.72389 -13574.5108
#BIC diff      0.00    -43.55636   -192.3433


# Se puede observar que se obtiene un mejor BIC sacando la variable shape,
# por lo que se opta por sacar esta variable.


mod=Mclust(features,x=BIC) # en base al mejor valor BIC se realiza el mclust
summary(mod)#se muestra resultado y tabla de clustering

#Graficando
plot(mod, what = "classification")  #se grafica la configuraci?n de agrupamientos.
legend("bottomright", legend = 1:9,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(class, mod$classification) #distribuci?n de clases por cada grupo.



diss.matrix = daisy(features, metric = "euclidean",stand = FALSE)

clusters = pam(diss.matrix,2,diss=TRUE, metric="euclidean")

summary(clusters)

clusplot(clusters)





