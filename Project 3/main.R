# startfil for å lese inn data
library(reshape2)
library(ggplot2)

# seismic:
seismic <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/seismic.dat")
seismic <- seismic[,1]
plot(seq(1,length(seismic)),seismic)

# complit:
complit <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/complit.dat")
com.arr <- as.matrix(complit)
com.arr <- melt(com.arr)
# flip y-axis:
com.arr <- data.frame(Var1 = (max(com.arr$Var1) - com.arr$Var1 + 1),
                      Var2 = com.arr$Var2, value = com.arr$value)

ggplot(data=com.arr, aes(x=Var2, y=Var1, fill = value)) + geom_tile() 
