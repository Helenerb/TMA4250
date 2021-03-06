# main
library(spatial)
library(MASS)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(datasets)

cells <- read.table("cells.dat.txt", col.names=c('x', 'y'))
redwood <- read.table("redwood.dat.txt", col.names=c('x',  'y'))
#TODO: VOJIN! Endre filen her til der den ligger hos deg! 
pines <- read.table(file = 'pines.dat.txt', skip=3, col.names = c('x','y'))


# a)
# Visualize plots:

# Cells:
cells.plot <- ggplot(data=cells, aes(x=x, y=y)) + geom_point(color="darkolivegreen4") +
  xlab(label=" ") + ylab(label=" ") + ggtitle(label="Cells")
cells.plot

# Redwood:
redwood.plot <- ggplot(data=redwood, aes(x=x, y=y)) + geom_point(color="sienna") +
  xlab(label=" ") + ylab(label=" ") + ggtitle(label="Redwood")
redwood.plot

# Pines:
pines.plot <- ggplot(data=pines, aes(x=x, y=y)) + geom_point(color="springgreen4") +
  xlab(label=" ") + ylab(label=" ") + ggtitle(label="Pines")
pines.plot


## b)

ppregion()
L.cells <- Kfn(cells, fs = sqrt(2))
L.pines <- Kfn(pines, fs = sqrt(2))
L.redwood <- Kfn(redwood, fs = sqrt(2))

plot(L.redwood, type="b")
plot(L.pines,type="b")
plot(L.cells, type="b")

## c)

# implementerer simuleringsalgoritme:

sim.stat.poiss <- function(no){
  x <- runif(no)
  y <- runif(no)
  df<- data.frame(x, y)
  return(df)
}

J.hat <- function(df, t.end){
  no <- length(df$x)
  ts = seq(0.01,t.end,by=0.01)
  result = rep(0,length(ts))
  for(k in 1:length(ts)){
    t <- ts[k]
    norm <- no*(pi*t^2)          # area of Ball with radius t with arbitrary center x. 
    # ^Not accounted for boundary
    sum.outer <- 0
    for( i in 1:no){
      x0 <- df$x[i]
      y0 <- df$y[i]
      sum.inner <- 0
      for (j in 1:no){
        if( sqrt((x0 - df$x[j])^2 + (y0 - df$y[j])^2) <= t){
          sum.inner <- sum.inner + 1
        }
      }
      sum.outer = sum.outer + sum.inner - 1
    }
    result[k] <- sum.outer/norm
  }
  return(data.frame(result))
}

lambda.hat <- function(df, t){
  norm <- pi*t^2 # not taking boundary into account! 
  # discretize area:
  xs <- seq(0.1,0.9,by=0.1)
  ys <- seq(0.1,0.9, by=0.1)
  grid <- expand.grid(xs,ys)
  result <- matrix(nrow=length(ys),ncol=length(xs)) # den er fyllt med NAs nå
  for(i in 1:length(xs)){
    for(j in 1:length(ys)){
      # regn ut!
      sum.inner <- 0
      for(n in 1:length(df$x)){
        if( sqrt((xs[i] - df$x[n])^2 + (ys[j] - df$y[n])^2) <= t){
          sum.inner <- sum.inner + 1
        }
      }
      result[i,j] <- sum.inner/norm
    }
  }
  return(result)
}

test.df <- sim.stat.poiss(42)
plot(test.df$x, test.df$y)
result.test <- J.hat(test.df,0.7)
lambda.test <- lambda.hat(test.df,0.2)

# simulation algoritm:

MC.sim <- function(obs.df, S, t.max, t.lamb){
  no <- length(obs.df$x)
  sim.J <- rep(0,t.max/0.01)
  sim.lamb <- matrix(0L,nrow=9,ncol=9)
  
  # for variance, storing simulation results:
  lamb.arr <- array(data=NA, dim=c(9,9,S))
  J.arr <- matrix(0L,nrow=S,ncol=length(sim.J))
  
  for(s in 1:S){
    cat("s:", s, " ")
    sim.df <- sim.stat.poiss(no)
    J.hat <- J.hat(sim.df,t.max)
    sim.J = sim.J + J.hat
    cat("length: ", length(J.hat), "\n")
    J.arr[s,] <- J.hat$result
    lambda.hat <- lambda.hat(sim.df,t.lamb)
    lamb.arr[,,s] <- lambda.hat
    sim.lamb = sim.lamb + lambda.hat
  }
  J.mean = sim.J/S
  lamb.mean = sim.lamb/S
  
  #finding the variance:
  J.var <- (J.mean-J.arr[1,])^2
  for(i in 2:S){
    J.var = J.var + (J.mean-J.arr[i,])^2
  }
  J.var = J.var/(S-1)
  
  lamb.var <- (lamb.mean - lamb.arr[,,1])^2
  for (j in 2:S){
    lamb.var <- lamb.var + (lamb.mean-lamb.arr[,,j])^2
  }
  lamb.var = lamb.var/(S-1)
  return(list("J.mean"=J.mean, "lamb.mean"=lamb.mean, "J.var"=J.var, "lamb.var" = lamb.var))
}


# redwood:
test.MC <- MC.sim(redwood, 100, 0.7, 0.3)
J.mean <- test.MC$J.mean
J.var <- test.MC$J.var
lamb.mean <- test.MC$lamb.mean
lamb.var <- test.MC$lamb.var

# confidence interval:
J.upper <- J.mean + 1.645*sqrt(J.var)/sqrt(100)
J.lower <- J.mean - 1.645*sqrt(J.var)/sqrt(100)

#for comparison:
#J.redwood <- J.hat(redwood, 0.7)
J.redwood <- data.frame("result" = J.hat(redwood, 0.7)$result,"t"=seq(0.01,0.7,by=0.01), "upper" = J.upper$result, "lower" = J.lower$result)
#lambda.redwood <- lambda.hat(redwood, 0.3)

# plot confidence interval:
conf.redwood <- ggplot(data=J.redwood) + geom_line(aes(x=t, y=result), color="sienna", size=0.5) + 
  geom_line(aes(y=upper, x=t), color="red") + geom_line(aes(y=lower, x=t), color="red") + 
  ggtitle(label="J-function for Redwood") + xlab(label="t") + ylab(label="J")
conf.redwood
  
# cells:
test.MC <- MC.sim(cells, 100, 0.7, 0.3)
J.mean <- test.MC$J.mean
J.var <- test.MC$J.var
lamb.mean <- test.MC$lamb.mean
lamb.var <- test.MC$lamb.var

# confidence interval:
J.upper <- J.mean + 1.645*sqrt(J.var)/sqrt(100)
J.lower <- J.mean - 1.645*sqrt(J.var)/sqrt(100)

#for comparison:
J.cells <- data.frame("result" = J.hat(cells, 0.7)$result,"t"=seq(0.01,0.7,by=0.01), "upper" = J.upper$result, "lower" = J.lower$result)

# plot confidence interval:
conf.cells <- ggplot(data=J.cells) + geom_line(aes(x=t, y=result), color="darkolivegreen4", size=0.5) + 
  geom_line(aes(y=upper, x=t), color="red") + geom_line(aes(y=lower, x=t), color="red") + 
  ggtitle(label="J-function for Cells") + xlab(label="t") + ylab(label="J")
conf.cells

# pines:
test.MC <- MC.sim(pines, 100, 0.7, 0.3)
J.mean <- test.MC$J.mean
J.var <- test.MC$J.var
lamb.mean <- test.MC$lamb.mean
lamb.var <- test.MC$lamb.var

# confidence interval:
J.upper <- J.mean + 1.645*sqrt(J.var)/sqrt(100)
J.lower <- J.mean - 1.645*sqrt(J.var)/sqrt(100)

#for comparison:
J.pines <- data.frame("result" = J.hat(pines, 0.7)$result,"t"=seq(0.01,0.7,by=0.01), "upper" = J.upper$result, "lower" = J.lower$result)

# plot confidence interval:
conf.pines <- ggplot(data=J.pines) + geom_line(aes(x=t, y=result), color="springgreen4", size=0.5) + 
  geom_line(aes(y=upper, x=t), color="red") + geom_line(aes(y=lower, x=t), color="red") + 
  ggtitle(label="J-function for Pines") + xlab(label="t") + ylab(label="J")
conf.pines

