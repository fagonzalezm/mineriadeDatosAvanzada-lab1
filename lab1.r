names <- c("id", "clumpThickness", "uniformityOfCellSize", "uniformityOfCellShape",
           "marginalAdhesion", "singleEpithelialCellSize", "bareNuclei",
           "blandChromatin", "normalNucleoli", "mitoses", "class")
data <- read.table("./breast-cancer-wisconsin.data", sep=",", col.names = names)

set.seed(20)


#######Quitar datos###########
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

################################


features <- data[,2:9]
class <- data$class

################mclust######################

#primer intento
require(mclust)
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
plot(BIC)  #se grafican los BIC por configuración de parámetros
summary(BIC)  # se presentan los mejores valores BIC

#Usando VVV,3
mod11=Mclust(features,x=BIC) # en base al mejor valor BIC se realiza el mclust
summary(mod11)#se muestra resultado y tabla de clustering

#Graficando
plot(mod11, what = "classification")  #se grafica la configuración de agrupamientos.
legend("bottomright", legend = 1:3,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(class, mod11$classification) #distribución de clases por cada grupo.


#Usando otros valores
mod12 = Mclust(features, G=7, prior = priorControl(functionName="defaultPrior", shrinkage=0.1), modelNames ="VVI")  
plot(mod12, what = "classification")  #se grafica la configuración de agrupamientos.
legend("bottomright", legend = 1:7,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")
table(class, mod12$classification) #distribución de clases por cada grupo.