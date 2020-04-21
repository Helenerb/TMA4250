# visualizing data in a)

library(reshape2)
library(ggplot2)
library(spatial)

seismic <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/seismic.dat")
seismic <- seismic[,1]
seismic.mat <- matrix(seismic, nrow = 75, ncol = 75)
seismic.melt <- melt(seismic.mat)

seismic.gg <- ggplot(data = seismic.melt, aes(x = Var1, y=Var2, fill = value)) + geom_tile()
seismic.gg <- seismic.gg + ggtitle(label="Seismic data") + xlab(label="") + ylab(label="") 
seismic.gg <- seismic.gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank())
seismic.gg
