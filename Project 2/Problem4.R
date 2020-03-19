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

# hold litt igjen på denne; du trenger kanskje ikke faktisk bruke den 
# p.Strauss <- function(X,Y,tau.0, phi.0, phi.1){
#   taus <- c()
#   for(i in 1:length(X)){
#     for(j in (i+1):length(Y)){
#       taus <- append(taus)
#     }
#   }
#   return(tau.0)
# }

Strauss <- function(tau.0, phi.0, phi.1, k, T, X.0, Y.0){
  X.s <- X.0                            # initialize variable to hold x values 
  Y.s <- Y.0                            # initialize variable to hold y values 
  X.D.x <- X.0
  X.D.y <- Y.0
  
  for(i in 1:T){                        # T: number of MCMC iterations
    u <- sample(seq(1:k),size=1)        # sample random event point
    x.p <- runif(1,-1,1)
    y.p <- runif(1,-1,1)
    
    # finding alpha:
    taus.p <- sqrt((x.p-X.D.x)^2 + (y.p - X.D.y)^2) # distance between proposed and other events
    taus.u <- sqrt((X.D.x[u]-X.D.x)^2 + (X.D.y[u] - X.D.y)^2) # distance between old and other events
    
    sum.phi <- sum(sapply(taus.p, FUN=Phi,tau.0=tau.0, phi.0=phi.0, phi.1=phi.1) -
                     sapply(taus.u, FUN=Phi, tau.0=tau.0, phi.0=phi.0, phi.1=phi.1))      # bruk apply hvis dette ikke går
    alpha <- exp(-sum.phi)
    alpha <- min(1,alpha)
    
    
    if(runif(1) < alpha){               # accept with propb alpha
      X.D.x[u] <- x.p
      X.D.y[u] <- y.p
    }
    X.s <- rbind(X.s,X.D.x)
    Y.s <- rbind(Y.s,X.D.y)
  }
  return(list("X" = X.D.x, "Y" = X.D.y, "X.s"=X.s, "Y.s"=Y.s))
}

cells <- read.table("cells.dat.txt", col.names=c('x', 'y'))
ppregion()
L.cells <- Kfn(cells, fs = sqrt(2))

test.Strss <- Strauss(tau.0 = 0.1, phi.0 = 10, phi.1 = 1, k = 4*length(cells$x), T=3, X.0 = cells$x, Y.0 = cells$y)

test.Strss$X
test.Strss$Y
test.Strss$X.s
test.Strss$Y.s



