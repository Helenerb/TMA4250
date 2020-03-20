# problem 1a and c, considering the L-function instead of the J-function:
library(spatial)
library(MASS)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(datasets)

cells <- read.table("cells.dat.txt", col.names=c('x', 'y'))
redwood <- read.table("redwood.dat.txt", col.names=c('x',  'y'))
pines <- read.table(file = 'pines.dat.txt', skip=3, col.names = c('x','y'))


# a)
# Visualize plots:

# Cells:
cells.plot <- ggplot(data=cells, aes(x=x, y=y)) + geom_point(color="darkolivegreen4") +
  xlab(label="x") + ylab(label="y") #+ ggtitle(label="Cells")
cells.plot

# Redwood:
redwood.plot <- ggplot(data=redwood, aes(x=x, y=y)) + geom_point(color="sienna") +
  xlab(label="x") + ylab(label="y") #+ ggtitle(label="Redwood")
redwood.plot

# Pines:
pines.plot <- ggplot(data=pines, aes(x=x, y=y)) + geom_point(color="springgreen4") +
  xlab(label="x") + ylab(label="y") #+ ggtitle(label="Pines")
pines.plot

ppregion()
L.cells <- Kfn(cells, fs = sqrt(2))
L.pines <- Kfn(pines, fs = sqrt(2))
L.redwood <- Kfn(redwood, fs = sqrt(2))


sim.stat.poiss <- function(no){
  x <- runif(no)
  y <- runif(no)
  df<- data.frame(x, y)
  return(df)
}


L.sim <- function(no, S){
  # initate matrix to store simulated values
  L.x <- matrix(0L, nrow=0, ncol=50)        
  L.y <- matrix(0L, nrow=0, ncol=50)
  
  for(s in 1:S){
    real <- sim.stat.poiss(no)          # returns df with columns "x" and "y"
    L.real <- Kfn(real, fs = sqrt(2)) 
    L.x <- rbind(L.x, L.real$x)
    L.y <- rbind(L.y, L.real$y)
    
  }
  L.mean = colMeans(L.y)
  L.min = colMins(L.y)
  L.max = colMaxs(L.y)
  
  return(list("L.x" = L.x, "L.y" = L.y, "mean" = L.mean, "min" = L.min, "max" = L.max ))
}

# Redwood:
# redwood <- read.table("redwood.dat.txt", col.names=c('x',  'y'))
# ppregion()
# L.redwood <- Kfn(redwood, fs = sqrt(2))
# 
# Stat.sim <- L.sim(length(redwood$x), 100)
# 
# Stat.lower <- Stat.sim$mean + 0.95*(Stat.sim$max - Stat.sim$mean)
# Stat.upper <- Stat.sim$mean - 0.95*(Stat.sim$mean - Stat.sim$min)
# 
# Stat.gg.df <- data.frame("redwood" = L.redwood$y, "x" = L.redwood$x, "mean" = Stat.sim$mean,
#                     "lower" = Stat.lower, "upper" = Stat.upper)
# 
# Stat.plot <- ggplot(data=Stat.gg.df) + geom_point(aes(x=x, y=redwood), color="sienna") +
#   geom_line(aes(x=x, y=upper), color="red", size=0.75) +
#   geom_line(aes(x=x, y=lower), color="red", size=0.75) +
#   xlab(label="t") + ylab(label="L")
# Stat.plot

# Cells:

# cells <- read.table("cells.dat.txt", col.names=c('x', 'y'))
# ppregion()
# L.cells <- Kfn(cells, fs = sqrt(2))
# 
# Stat.sim <- L.sim(length(cells$x), 100)
# 
# Stat.lower <- Stat.sim$mean + 0.95*(Stat.sim$max - Stat.sim$mean)
# Stat.upper <- Stat.sim$mean - 0.95*(Stat.sim$mean - Stat.sim$min)
# 
# Stat.gg.df <- data.frame("cells" = L.cells$y, "x" = L.cells$x, "mean" = Stat.sim$mean,
#                     "lower" = Stat.lower, "upper" = Stat.upper)
# 
# Stat.plot <- ggplot(data=Stat.gg.df) + geom_point(aes(x=x, y=cells), color="darkolivegreen4") +
#   geom_line(aes(x=x, y=upper), color="red", size=0.75) +
#   geom_line(aes(x=x, y=lower), color="red", size=0.75) +
#   xlab(label="t") + ylab(label="L")
# Stat.plot

# pines:
pines <- read.table(file = 'pines.dat.txt', skip=3, col.names = c('x','y'))
ppregion()
L.pines <- Kfn(pines, fs = sqrt(2))

Stat.sim <- L.sim(length(pines$x), 100)

Stat.lower <- Stat.sim$mean + 0.95*(Stat.sim$max - Stat.sim$mean)
Stat.upper <- Stat.sim$mean - 0.95*(Stat.sim$mean - Stat.sim$min)

Stat.gg.df <- data.frame("pines" = L.pines$y, "x" = L.pines$x, "mean" = Stat.sim$mean,
                         "lower" = Stat.lower, "upper" = Stat.upper)

Stat.plot <- ggplot(data=Stat.gg.df) + geom_point(aes(x=x, y=pines), color="springgreen4") +
  geom_line(aes(x=x, y=upper), color="red", size=0.75) +
  geom_line(aes(x=x, y=lower), color="red", size=0.75) +
  xlab(label="t") + ylab(label="L")
Stat.plot

