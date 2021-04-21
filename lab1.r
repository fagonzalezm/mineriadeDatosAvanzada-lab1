names <- c("id", "clumpThickness", "uniformityOfCellSize", "uniformityOfCellShape",
           "marginalAdhesion", "singleEpithelialCellSize", "bareNuclei",
           "blandChromatin", "normalNucleoli", "mitoses", "class")
data <- read.table("./breast-cancer-wisconsin.data", sep=",", col.names = names)
