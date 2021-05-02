require(mclust)
require(reshape2)
require(ggplot2)
require(corrplot)

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

# - Se cuenta con 683 datos: 444 c茅lulas benignas y 239 cancerosas.
# - Se descartaron 14 c茅lulas benignas y 2 cancerosas.

####### Elecci贸n de variables de inter茅s #######
# - La variable id no proporciona informaci贸n 煤til para el estudio, as铆 que no se considera
# - La variable class se usa como gold standard, por lo que su uso se limita a evaluar el modelo,
#   mas no a su confecci贸n.
# - Se estudia la distribuci贸n del resto de las variables para decidir si ser谩n consideradas en
#   el modelo

# Se quita la variable id
features <- data[,2:11]

# Se transforman los 2 y 4 de la variable class a "B" y "M" seg煤n corresponda
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

# Se observa que para todas las variables existe una relaci贸n entre sus valores
# y si la c茅lula corresponde a una cancerosa o no, excepto con mitoses.

# Estadistica descriptiva
# meadias
means <- sapply(features[, 1:9], mean)

# medianas
medians <- sapply(features[, 1:9], median)

#varianzas
vars <- sapply(features[, 1:9], var)

#coeficientes de variacion
coefs <- sqrt(vars)/means



####### Test de normalidad #######
# H0: La muestra sigue una distribuci贸n normal
# H1: La muestra no sigue una distribuci贸n normal

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
#  - Por lo tanto se rechaza la hip贸tesis nula.
#  - Existe evidencia estad铆stica suficiente para rechazar H0, entonces.
#    se puede afirmar que las varaibles no siguen una distribuci贸n normal.

# Esto se puede confirmar revisando sus respectivos histogramas

p <- ggplot(data = features.melt, aes(x=value, fill=variable))
p <- p + geom_bar()
p <- p + facet_wrap( ~ variable, scales="free")
p

#Ninguna sigue una distribuici髇 normal, se opta por un test no param閠rico
data$bareNuclei <- as.numeric(data$bareNuclei)
data.cor = cor(data, method = 'spearman')

corrplot(data.cor, method = 'ellipse')
title("Matriz de Correlacion")


features <- data[,2:10]
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

#Usando VVV,4
mod11=Mclust(features,x=BIC) # en base al mejor valor BIC se realiza el mclust
summary(mod11)#se muestra resultado y tabla de clustering

#Graficando
plot(mod11, what = "classification")  #se grafica la configuraci?n de agrupamientos.
legend("bottomright", legend = 1:4,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(class, mod11$classification) #distribuci?n de clases por cada grupo.


#Usando segundo mejor BIC
mod12 = Mclust(features, G=5, prior = priorControl(functionName="defaultPrior", shrinkage=0.1), modelNames ="VVI")  
plot(mod12, what = "classification")  #se grafica la configuraci?n de agrupamientos.
legend("bottomright", legend = 1:5,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")
table(class, mod12$classification) #distribuci?n de clases por cada grupo.



#Quitando la variable mitosis.
features <- data[,2:9]
class <- data$class

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