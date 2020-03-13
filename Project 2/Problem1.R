# main
library(spatial)
library(MASS)
library(ggplot2)

cells <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/cells.dat")
redwood <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/redwood.dat")
#TODO: VOJIN! Endre filen her til der den ligger hos deg! 
pines <- read.table("C://Users//helen//Documents//GitHub//Spatial Statistics//TMA4250//Project 2//pines.txt")

# a)
# Visualize plots:

# Cells:
cells.plot <- ggplot(data=cells, aes(x=V1, y=V2)) + geom_point(color="darkolivegreen4") +
  xlab(label=" ") + ylab(label=" ") + ggtitle(label="Cells")
cells.plot

# Redwood:
redwood.plot <- ggplot(data=redwood, aes(x=V1, y=V2)) + geom_point(color="sienna") +
  xlab(label=" ") + ylab(label=" ") + ggtitle(label="Redwood")
redwood.plot

# Pines:
pines.plot <- ggplot(data=pines, aes(x=V1, y=V2)) + geom_point(color="springgreen4") +
  xlab(label=" ") + ylab(label=" ") + ggtitle(label="Pines")
pines.plot


## b)
pp.cells <- ppinit("cells.dat")
pp.pines <- ppinit("pines.dat")
pp.redwood <- ppinit("redwood.dat")

L.cells <- Kfn(pp.cells, 1)
L.pines <- Kfn(pp.pines, sqrt(2))
L.redwood <- Kfn(pp.redwood, sqrt(2))

as.data.frame(L.cells)

L.plot <- ggplot() + geom_point(data=as.data.frame(L.cells),  aes(x=x, y=y), color="darkolivegreen4")

read.table("towns.dat")


