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

torus.2d <- function(x){
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

NS <- function(lamb.M, sigma.c, p.mu, p.sigma, lamb.p){ #DEBUG
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
    
    #k.C <- rnorm(1,mean = p.mu, sd=sqrt(p.sigma))             # sample no of child points 
    k.C <- rpois(1,lamb.p)                      # DEBUG
    
    # sample k.C child points
    for(i in 1:k.C){
      x.C <- mvrnorm(1,c(x.M, y.M),sigma.c*diag(2))
      
      # transform to torus representation of D
      x.C <- torus.2d(x.C)
      
      x.C.all <- rbind(x.C.all, x.C)       # store child events
    }
  }
  return(list("x.C.all" = x.C.all, "x.mother"=x.mother))
}

NS.sim <- function(lamb.M, sigma.c, p.mu, p.sigma, lamb.p, S){ #DEBUG
  # initate matrix to store simulated values
  L.x <- matrix(0L, nrow=0, ncol=50)                  # is 50 general?
  L.y <- matrix(0L, nrow=0, ncol=50)
  
  for(s in 1:S){
    NS <- NS(lamb.M, sigma.c, p.mu, p.sigma, lamb.p)
    
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

# first guestimate:
#NS.sim <- NS.sim(7, 0.125^2, 8, 10, 100)

# second guestimate:
#NS.sim <- NS.sim(7, 0.2^2, 8, 10, 100)

# third guestimate:
#NS.sim <- NS.sim(7, 0.125^2, 6, 6, 100)

# fourth guestimate:
# NS.sim <- NS.sim(7, 0.1^2, 6, 6, 100)

# fifth guestimate:
#NS.sim <- NS.sim(7, 0.05^2, 6, 4, 100)

# sixth guestimate:
#NS.sim <- NS.sim(7, 0.05^2, 5, 4, 100) #gir for høye estimater nært

#seventh:
#dytter trærne lenger vekk
NS.sim <- NS.sim(lamb.M = 12, sigma.c = 0.05^2, 5, 4, lamb.p = 5, 100) #funker faktisk ganske bra!

# 95% confidence interval:
NS.upper <- NS.sim$mean + 0.9*(NS.sim$max - NS.sim$mean)
NS.lower <- NS.sim$mean - 0.9*(NS.sim$mean - NS.sim$min)

NS.lower.regular <- NS.sim$mean - 1.645*sqrt(NS.sim$var)/sqrt(100)
NS.upper.regular <- NS.sim$mean + 1.645*sqrt(NS.sim$var)/sqrt(100)

# plot confidence interval together with redwood:

#collect in dataframe for ggplot:

gg.df <- data.frame("redwood" = L.redwood$y, "x" = L.redwood$x, "mean" = NS.sim$mean,
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
  
  plot(x,y, col="red", xlim=c(0,1), ylim=c(0,1))
  points(x.M,y.M, col="green")
}

# fourth guestimate:
plot.real.NS(real = NS(12, 0.05^2, 6, 6, 5))

# fifth guestimate:
plot.real.NS(real = NS(7, 0.05^2, 6, 4, 5))

# sixth guestimate:
plot.real.NS(real = NS(7, 0.05^2, 4, 4))

# sevent guestimate:
plot.real.NS(real = NS(7, 0.1^2, 5, 4))






