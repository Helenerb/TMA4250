# problem 4
library(matrixStats)
library(spatial)
library(MASS)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(datasets)
library(MASS)

Phi <- function(tau, tau.0, phi.0, phi.1){
  if (tau < 0){
    cat("uffda, du har fått negativ tau")
    return(0)
  }
  else if(tau < tau.0){
    return(phi.0)
  }
  else{
    return(phi.0*exp(-phi.1*(tau - tau.0)))
  }
}


Strauss <- function(tau.0, phi.0, phi.1, k, T, X.0, Y.0){
  X.s <- X.0                            # initialize variable to hold x values 
  Y.s <- Y.0                            # initialize variable to hold y values 
  X.D.x <- X.0
  X.D.y <- Y.0
  accprob = 0
  
  for(i in 1:T){                        # T: number of MCMC iterations
    u <- sample(seq(1:k),size=1)        # sample random event point
    #x.p <- runif(1,-1,1)
    #y.p <- runif(1,-1,1)
    x.p <- runif(1,-0.5,0.5)
    y.p <- runif(1,-0.5,0.5)
    
    # finding alpha:
    taus.p <- sqrt((x.p-X.D.x)^2 + (y.p - X.D.y)^2) # distance between proposed and other events
    taus.u <- sqrt((X.D.x[u]-X.D.x)^2 + (X.D.y[u] - X.D.y)^2) # distance between old and other events
    
    taus.p <- taus.p[-u]                
    taus.u <- taus.u[-u]
    if(i == 2){
      cat(length(taus.u), length(cells$x))
    }
    
    sum.phi <- sum(sapply(taus.p, FUN=Phi,tau.0=tau.0, phi.0=phi.0, phi.1=phi.1) -
                     sapply(taus.u, FUN=Phi, tau.0=tau.0, phi.0=phi.0, phi.1=phi.1))      
    alpha <- min(1,alpha)
    
    
    if(runif(1) < alpha){               # accept with propb alpha
      X.D.x[u] <- x.p
      X.D.y[u] <- y.p
      accprob  = accprob + 1
    }
    X.s <- rbind(X.s,X.D.x)
    Y.s <- rbind(Y.s,X.D.y)
  }
  accprob = accprob/T
  return(list("X" = X.D.x, "Y" = X.D.y, "X.s"=X.s, "Y.s"=Y.s, "accprob"=accprob))
}

cells <- read.table("cells.dat.txt", col.names=c('x', 'y'))
ppregion()
L.cells <- Kfn(cells, fs = sqrt(2))

# making init from cells to decrease time to convergence:
X.init <- append(append(append(cells$x,cells$x),cells$x-1),cells$x-1)
Y.init <- append(append(append(cells$y,cells$y-1),cells$y),cells$y-1)
init.df = data.frame("x" = X.init, "y" = Y.init)
init.df = init.df[which(X.init > 0.5 | X.init < -0.5 | Y.init > 0.5 | Y.init < -0.5),]

X.init = append(runif(length(cells$x),-0.5,0.5),init.df$x)
Y.init = append(runif(length(cells$x),-0.5,0.5),init.df$y)
plot(X.init, Y.init)

# plotting plotting realizations with different parameters:

#guess pushing too much into edges:
#test.Strss <- Strauss(tau.0 = 0.1, phi.0 = 10, phi.1 = 1, k = length(cells$x), T=500, X.0 = X.init, Y.0 = Y.init)

# guess with lower tau.0, but less rapidsly decreasing phi,
# test.Strss <- Strauss(tau.0 = 0.02, phi.0 = 100, phi.1 = 50, k = length(cells$x), T=2000, X.0 = X.init, Y.0 = Y.init)

#trying to improve last:
#test.Strss <- Strauss(tau.0 = 0.075, phi.0 = 1.5, phi.1 = 3, k = length(cells$x), T=2000, X.0 = X.init, Y.0 = Y.init)

# testing how our final result looks:
test.Strss <- Strauss(tau.0 = 0.02, phi.0 = 100, phi.1 = 50, k = length(cells$x), T=1000, X.0 = X.init, Y.0 = Y.init)

plot(test.Strss$X,test.Strss$Y)
test.Strss$accprob

Strauss.sim <- function(S, X.0, Y.0, tau.0, phi.0, phi.1, T){
  L.x <- matrix(0L, nrow=0, ncol=50)                        # matrix storing simulations
  L.y <- matrix(0L, nrow=0, ncol=50)
  for(s in 1:S){
    real <- Strauss(tau.0, phi.0, phi.1, k = length(cells$x), T, X.0, Y.0)
    
    D.plus.df = data.frame("x"=real$X, "y"=real$Y)          # extended domain
    D.df = D.plus.df[which(D.plus.df$x > -0.5 &             # original domain
                             D.plus.df$x < 0.5 & D.plus.df$y < 0.5 &
                             D.plus.df$y > -0.5),]
    D.df = data.frame("x"=D.df$x + 0.5, "y"=D.df$y + 0.5)   # shifting to [0,1]x[0,1]
    D.L <- Kfn(D.df, fs=sqrt(2))                            # finding L function
    L.x <- rbind(L.x, D.L$x)
    L.y <- rbind(L.y, D.L$y)
  }
  L.mean = colMeans(L.y)
  L.min = colMins(L.y)
  L.max = colMaxs(L.y)
  L.var = colVars(L.y)
  return(list("L.x" = L.x, "L.y" = L.y, "mean" = L.mean, "min" = L.min, "max" = L.max, "var" = L.var ))
}

compare.sims <- function(sim, L.cells){
  upper <- sim$mean + 0.9*(sim$max - sim$mean)              # upper 90% span
  lower <- sim$mean - 0.9*(sim$mean - sim$min)              # lower 90% span
  gg.df <- data.frame("cells" = L.cells$y, "x" = L.cells$x, "mean" = sim$mean,
                      "lower" = lower, "upper" = upper)
  
  NS.plot <- ggplot(data=gg.df) + geom_point(aes(x=x, y=cells), color="darkolivegreen4") + 
    geom_line(aes(x=x, y=upper), color="red", size=0.75) + 
    geom_line(aes(x=x, y=lower), color="red", size=0.75) + 
    xlab(label="t") + ylab(label="L")
  NS.plot
}

#TODO; Testing different parameters. Note: long run-time, only one at a time

#dette er guess 2; trenger noe som er mer repulsive!
#test.sim <- Strauss.sim(100, X.init, Y.init, tau.0=0.1, phi.0=2, phi.1=10, 500) 

#guess 3:
#test.sim <- Strauss.sim(100, X.init, Y.init, tau.0=0.1, phi.0=5, phi.1=10, 1000) 

#guess 4:
#test.sim <- Strauss.sim(20, X.init, Y.init, tau.0=0.125, phi.0=10, phi.1=20, 1000) 

#guess 5: for repulsive langt unna?
#test.sim <- Strauss.sim(20, X.init, Y.init, tau.0 = 0.05, phi.0 = 2, phi.1 = 3, 1000) 

#guess 6: kanskje forbedring av forrige: nei, denne var ikke spesielt mye bedre...
#test.sim <- Strauss.sim(20, X.init, Y.init, tau.0 = 0.075, phi.0 = 1.5, phi.1 = 3, 1000)

#guess 7: høyere verdier:
#test.sim <- Strauss.sim(20, X.init, Y.init, tau.0 = 0.02, phi.0 = 100, phi.1 = 50, 1000)

#guess 8: høyere verdier: bruk dette som second guess!!
#test.sim <- Strauss.sim(100, X.init, Y.init, tau.0 = 0.02, phi.0 = 10, phi.1 = 5, 1000)

#øk phi.0 og phi.1 
#test.sim <- Strauss.sim(100, X.init, Y.init, tau.0 = 0.02, phi.0 = 80, phi.1 = 40, 700)

#øk litt mer, lavere tau.0:
#test.sim <- Strauss.sim(100, X.init, Y.init, tau.0 = 0.015, phi.0 = 120, phi.1 = 60, 700)

#guess 12: øker tau.0 og phi.1 litt etter at du har tenkt deg om: fortsatt litt høy ved lave t
#test.sim <- Strauss.sim(50, X.init, Y.init, tau.0 = 0.025, phi.0 = 100, phi.1 = 70, 700)

# guess 13: øker phi.1 litt fra forrige: da ble det faktisk verrre
#test.sim <- Strauss.sim(50, X.init, Y.init, tau.0 = 0.025, phi.0 = 100, phi.1 = 80, 700)

# øker phi.0 litt
#test.sim <- Strauss.sim(50, X.init, Y.init, tau.0 = 0.025, phi.0 = 110, phi.1 = 80, 700)

#går tilbake på tau.0: her ligger de ganske godt innenfor langt oppe
#test.sim <- Strauss.sim(50, X.init, Y.init, tau.0 = 0.02, phi.0 = 110, phi.1 = 80, 700)

#prøver mer repulsive: det hadde motsatt effekt
#test.sim <- Strauss.sim(50, X.init, Y.init, tau.0 = 0.02, phi.0 = 110, phi.1 = 90, 700)

#går tilbake på phi.1, øker litt på phi.0 i stedet. Det ble litt bedre, fortsetter sånn
#test.sim <- Strauss.sim(50, X.init, Y.init, tau.0 = 0.02, phi.0 = 120, phi.1 = 80, 700)

#øker phi.0 enda mer. marginalt bedre
#test.sim <- Strauss.sim(50, X.init, Y.init, tau.0 = 0.02, phi.0 = 130, phi.1 = 80, 700)

#lowering the phi values a little:
#test.sim <- Strauss.sim(100, X.init, Y.init, tau.0 = 0.02, phi.0 = 110, phi.1 = 70, 700)

#even more: This is the best so far?
test.sim <- Strauss.sim(100, X.init, Y.init, tau.0 = 0.02, phi.0 = 100, phi.1 = 50, 700)

compare.sims(test.sim, L.cells)

plot.Strss <- function(real){
  D.plus.df <- data.frame("x"=real$X, "y"=real$Y)
  D.df = D.plus.df[which(D.plus.df$x > -0.5 &
                           D.plus.df$x < 0.5 & D.plus.df$y < 0.5 &
                           D.plus.df$y > -0.5),]
  D.df = data.frame("x"=D.df$x + 0.5, "y"=D.df$y + 0.5)
  real.plot <- ggplot(data=D.df) + geom_point(aes(x=x, y=y), color="darkolivegreen3")
  real.plot
}

# plot realizations:
real.Strss <- Strauss(tau.0 = 0.02, phi.0 = 100, phi.1 = 50, k = length(cells$x), T=1000, X.0 = X.init, Y.0 = Y.init)
plot.Strss(real.Strss)

