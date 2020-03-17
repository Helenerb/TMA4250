# Problem 3

library(spatial)
library(MASS)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(datasets)

redwood <- read.table("redwood.dat.txt", col.names=c('x',  'y'))

redwood.plot <- ggplot(data=redwood, aes(x=x, y=y)) + geom_point(color="sienna") +
  xlab(label=" ") + ylab(label=" ") + ggtitle(label="Redwood")
redwood.plot

ppregion()
L.redwood <- Kfn(redwood, fs = sqrt(2))

# Simulation algorithm for Neumann-Scott 


NS <- function(lamb.M, sigma.c, p.mu, p.sigma){
  k <- 0
  k.M <- rpois(1,lamb.M)
}
