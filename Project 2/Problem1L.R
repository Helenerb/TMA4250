# problem 1c, considering the L-function instead of the J-function:

L.sim <- function(no, S){
  # initate matrix to store simulated values
  L.x <- matrix(0L, nrow=0, ncol=50)                  # is 50 general?
  L.y <- matrix(0L, nrow=0, ncol=50)
  
  for(s in 1:S){
    # bytt ut med en simulering av stat poiss med k
    #NS <- NS(lamb.M, sigma.c, p.mu, p.sigma)
    real <- sim.stat.poiss(no)          # returns df with columns "x" and "y"
    
    #x <- NS$x.C.all[,"x"]
    #y <- NS$x.C.all[,"y"]
    
    #NS.df <- data.frame("x" = x, "y" = y)
    L.real <- Kfn(real, fs = sqrt(2))     # kanskje 1???
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

Stat.gg.df <- data.frame("redwood" = L.pines$y, "x" = L.pines$x, "mean" = Stat.sim$mean,
                         "lower" = Stat.lower, "upper" = Stat.upper)

Stat.plot <- ggplot(data=Stat.gg.df) + geom_point(aes(x=x, y=pines), color="springgreen4") +
  geom_line(aes(x=x, y=upper), color="red", size=0.75) +
  geom_line(aes(x=x, y=lower), color="red", size=0.75) +
  xlab(label="t") + ylab(label="L")
Stat.plot

