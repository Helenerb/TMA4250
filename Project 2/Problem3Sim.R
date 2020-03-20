# NS simulation
library(matrixStats)
library(spatial)
library(MASS)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(datasets)
library(MASS)

torus.2d <- function(x){ #transforms to torus domain
  while(x[1] < 0 | x[1] > 1 | x[2] < 0 | x[2] > 1){
    if(x[1] < 0){
      x[1] <- 1 + x[1]
    }
    if(x[1] > 1){
      x[1] <- x[1] - 1
    }
    if(x[2] < 0){
      x[2] <- 1 + x[2]
    }
    if(x[2] > 1){
      x[2] <- x[2] - 1
    }
  }
  return(x)
}

NS <- function(lamb.M, sigma.c, lamb.p){
  k <- 0
  k.M <- rpois(1,lamb.M)                    # sample no of mother points
  
  # initiate matrix to hold child events 
  x.C.all <- matrix(0L, nrow = 0, ncol=2)
  colnames(x.C.all) <- c("x", "y")
  
  #initiate matrix to hold mother events
  x.mother <- matrix(0L, nrow = 0, ncol=2)
  colnames(x.mother) <- c("x", "y")
  
  for(j in 1:k.M){
    # sample location of mother node
    x.M <- runif(1)
    y.M <- runif(1)
    x.mother <- rbind(x.mother, c(x.M, y.M)) # store mother event
    
    k.C <- rpois(1,lamb.p)                   # sample k.C child points              
    for(i in 1:k.C){
      x.C <- mvrnorm(1,c(x.M, y.M),sigma.c*diag(2))
      
      x.C <- torus.2d(x.C)                  # torus representation of D
      
      x.C.all <- rbind(x.C.all, x.C)        # store child events
    }
  }
  return(list("x.C.all" = x.C.all, "x.mother"=x.mother))
}

NS.sim <- function(lamb.M, sigma.c, lamb.p, S){ 
  # initate matrix to store simulated values
  L.x <- matrix(0L, nrow=0, ncol=50)                  
  L.y <- matrix(0L, nrow=0, ncol=50)
  
  for(s in 1:S){
    NS <- NS(lamb.M, sigma.c, lamb.p)
    
    x <- NS$x.C.all[,"x"]
    y <- NS$x.C.all[,"y"]
    
    NS.df <- data.frame("x" = x, "y" = y)
    L.NS <- Kfn(NS.df, fs = sqrt(2))
    L.x <- rbind(L.x, L.NS$x)
    L.y <- rbind(L.y, L.NS$y)
    
  }
  L.mean = colMeans(L.y)
  L.min = colMins(L.y)
  L.max = colMaxs(L.y)
  L.var = colVars(L.y)
  
  return(list("L.x" = L.x, "L.y" = L.y, "mean" = L.mean, "min" = L.min, "max" = L.max, "var" = L.var ))
}


redwood <- read.table("redwood.dat.txt", col.names=c('x',  'y'))
ppregion()
L.redwood <- Kfn(redwood, fs = sqrt(2))


#potensielt første, sånn som du tenkte tidligere: generelt for høy tetthet
#NS.sim.result <- NS.sim(lamb.M = 7, sigma.c = 0.125^2, lamb.p = 8, 100)

#potensielt andre, minke verdi for sigma: innenfor, men fortsatt på den høye siden
#NS.sim.result <- NS.sim(lamb.M = 10, sigma.c = 0.05^2, lamb.p = 6, 100)

# øker litt til:
#NS.sim.result <- NS.sim(lamb.M = 11, sigma.c = 0.05^2, lamb.p = 5, 100)

#prøver å minke enda litt mer: fortsatt for høyt!
#NS.sim.result <- NS.sim(lamb.M = 12, sigma.c = 0.05^2, lamb.p = 4, 100)

#prøver å senke variansen litt
#NS.sim.result <- NS.sim(lamb.M = 12, sigma.c = 0.03^2, lamb.p = 4, 100)

#pøver å øke variansen i stedet:
#NS.sim.result <- NS.sim(lamb.M = 12, sigma.c = 0.06^2, lamb.p = 4, 100)

#prøver nå med lamb.p 5 i stedet:
#NS.sim.result <- NS.sim(lamb.M = 12, sigma.c = 0.05^2, lamb.p = 5, 100)

#minker variansen bittelitt: så langt beste? Valgt som beste
NS.sim.result <- NS.sim(lamb.M = 12, sigma.c = 0.045^2, lamb.p = 5, 100)

#prøver med 4 og 12:
#NS.sim.result <- NS.sim(lamb.M = 12, sigma.c = 0.045^2, lamb.p = 4, 100)

# 95% confidence interval:
NS.upper <- NS.sim.result$mean + 0.9*(NS.sim.result$max - NS.sim.result$mean)
NS.lower <- NS.sim.result$mean - 0.9*(NS.sim.result$mean - NS.sim.result$min)

NS.lower.regular <- NS.sim.result$mean - 1.645*sqrt(NS.sim.result$var)/sqrt(100)
NS.upper.regular <- NS.sim.result$mean + 1.645*sqrt(NS.sim.result$var)/sqrt(100)

# plot confidence interval together with redwood:

#collect in dataframe for ggplot:

gg.df <- data.frame("redwood" = L.redwood$y, "x" = L.redwood$x, "mean" = NS.sim.result$mean,
                    "lower" = NS.lower, "upper" = NS.upper,
                    "upper.reg" = NS.upper.regular, "lower.reg" = NS.lower.regular)

NS.plot <- ggplot(data=gg.df) + geom_point(aes(x=x, y=redwood), color="sienna") + 
  geom_line(aes(x=x, y=upper), color="red", size=0.75) + 
  geom_line(aes(x=x, y=lower), color="red", size=0.75) + 
  #geom_line(aes(x=x, y=upper.reg), color="green", size=0.75) + 
  #geom_line(aes(x=x, y=lower.reg), color="green", size=0.75) + 
  xlab(label="t") + ylab(label="L")
NS.plot

# function forplotting one realization
plot.real.NS <- function(real){
  x <- real$x.C.all[,"x"]
  y <- real$x.C.all[,"y"]
  
  x.M <- real$x.mother[,1]
  y.M <- real$x.mother[,2]
  
  gg.df <- data.frame("x" = x, "y"=y)
  gg.mother.df <- data.frame("x.M" = x.M, "y.M" = y.M)
  
  real.plot <- ggplot(data=gg.df) + geom_point(aes(x=x, y=y), color="sienna3")
  real.plot <- real.plot +geom_point(data=gg.mother.df, aes(x=x.M, y=y.M), color="limegreen")
  real.plot
}

# plot final guestimate:
plot.real.NS(real = NS(12, 0.045^2, lamb.p = 5))






