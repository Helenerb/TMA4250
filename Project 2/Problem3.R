# Problem 3

library(spatial)
library(MASS)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(datasets)
library(MASS)

redwood <- read.table("redwood.dat.txt", col.names=c('x',  'y'))

redwood.plot <- ggplot(data=redwood, aes(x=x, y=y)) + geom_point(color="sienna") +
  xlab(label=" ") + ylab(label=" ") + ggtitle(label="Redwood")
redwood.plot

ppregion()
L.redwood <- Kfn(redwood, fs = sqrt(2))

# Simulation algorithm for Neumann-Scott 

# torus.2d help function:

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


NS <- function(lamb.M, sigma.c, p.mu, p.sigma){
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
    
    k.C <- rnorm(1,mean = p.mu, sd=sqrt(p.sigma))             # sample no of child points 
    
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

test <- NS(7, 0.125^2, 8, 10)
test$x.C.all
x <- test$x.C.all[,"x"]
y <- test$x.C.all[,"y"]

#making it into a dataframe, to be able to compute L! 
test.df <- data.frame("x" = x, "y" = y)

x.M <- test$x.mother[,1]
y.M <- test$x.mother[,2]

plot(x,y, col="red", xlim=c(0,1), ylim=c(0,1))
points(x.M,y.M, col="green")
test.df <- data.frame("x"=x,"y"=y)
