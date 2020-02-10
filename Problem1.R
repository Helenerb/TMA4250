library(fields); 
library(geoR);  
library(MASS)
library(akima)
library(mvtnorm)
library(matlib)

d <- 0:49    # All possible distances 
phi <- 10    # Usikker på denne, men tror det kommer av def av rho i oppgaven 

# Powered exponential:

nu.exp <- c(1,1.9)

corr.exp1 <- exp(-(d/phi)^nu.exp[1])
corr.exp1.9 <- exp(-(d/phi)^nu.exp[2])

plot(d/10,corr.exp1,type="l")
lines(d/10,corr.exp1.9,type="l",col="blue")

#Matern: 

nu.mat <- c(1,3)

corr.mat1 = Matern(d,range=10,nu=nu.mat[1])  # usikker på om det er riktig at range=10
corr.mat3 = Matern(d,range=10,nu=nu.mat[2],phi=1)

lines(d/10,corr.mat1,type="l",col="red")  # usikker på om vi skal plotte som funk av d eller d/10 = tau?
lines(d/10,corr.mat3,type="l",col="green")

#Related variograms:
sigma <- c(1,5)

var.exp1.1 <- sigma[1]*(1-corr.exp1)
var.exp1.5 <- sigma[2]*(1-corr.exp1)

var.exp1.9.1 <- sigma[1]*(1-corr.exp1.9)
var.exp1.9.5 <- sigma[2]*(1-corr.exp1.9)

var.mat1.1 <- sigma[1]*(1-corr.mat1)
var.mat1.5 <- sigma[2]*(1-corr.mat1)

var.mat3.1 <- sigma[1]*(1-corr.mat3)
var.mat3.5 <- sigma[2]*(1-corr.mat3)

#Plot correlation function with corresponding variograms, variance {1,5}

plot(d/10,var.exp1.5,type="l")
lines(d/10,var.exp1.1,type="l",col="orange")
lines(d/10,corr.exp.1,type="l",col="blue")

plot(d/10,var.exp1.9.5,type="l")
lines(d/10,var.exp1.9.1,type="l",col="orange")
lines(d/10,corr.exp1.9,type="l",col="blue")

plot(d/10,var.mat1.5,type="l")
lines(d/10,var.mat1.1,type="l",col="orange")
lines(d/10,corr.mat1,type="l",col="blue")

plot(d/10,var.mat3.5,type="l")
lines(d/10,var.mat3.1,type="l",col="orange")
lines(d/10,corr.mat3,type="l",col="blue")


# takes in vector of covariance between distances, and returns the corresponding covariance matrix. 
getCovMatrix <- function(cov.vec){
  cov.matrix <- matrix(data=0,nrow=length(cov.vec),ncol=length(cov.vec))
  for(i in 1:length(cov.vec)){
    for(j in 1:length(cov.vec)){
      cov.matrix[i,j] = cov.vec[abs(i-j)+1]
    }
  }
  return(cov.matrix)
}

get.prior.real <- function(model,sigma.sq,phi,nu,L=50){
  cov.vec <- cov.spatial(0:(L-1),cov.model=model,cov.pars=c(sigma.sq,phi),kappa=nu) # vector of covariances
  cov.matrix <- getCovMatrix(cov.vec)
  prior.real <- rmvnorm(4,rep(0,L),cov.matrix)
  return(prior.real)
}


plot.reals <- function(prior.real,title){
  plot(d,rep(0,length(prior.real[1,])),type="l",
       main=title,xlab="x",ylab="r(x)",ylim=c(min(prior.real),max(prior.real)))
  colors <- c("orange","red","green","blue")
  for(i in 1:4){
    lines(d,prior.real[i,],type="l",col=colors[i])
  }
}


# Powered exponential, nu = 1, var = 1:
prior.exp.1.1 <- get.prior.real("powered.exponential", sigma.sq=1, phi, nu=1)
plot.reals(prior.exp.1.1,"Powered Exponential, nu = 1, var = 1")

# Powered exponential, nu = 1, var = 5:
prior.exp.1.5 <- get.prior.real("powered.exponential", sigma.sq=5, phi, nu=1)
plot.reals(prior.exp.1.5,"Powered Exponential, nu = 1, var = 5")

# Powered exponential, nu = 1.9, var = 1:
prior.exp.19.1 <- get.prior.real("powered.exponential", sigma.sq=1, phi, nu=1.9)
plot.reals(prior.exp.19.1,"Powered Exponential, nu = 1.9, var = 1")

# Powered exponential, nu = 1.9, var = 5:
prior.exp.19.5 <- get.prior.real("powered.exponential", sigma.sq=5, phi, nu=1.9)
plot.reals(prior.exp.19.5,"Powered Exponential, nu = 1.9, var = 5")

# Matérn, nu = 1, var = 1:
prior.mat.1.1 <- get.prior.real("matern", sigma.sq=1, phi, nu=1)
plot.reals(prior.mat.1.1,"Matérn, nu = 1, var = 1")

# Matérn, nu = 1, var = 5:
prior.mat.1.5 <- get.prior.real("matern", sigma.sq=5, phi, nu=1)
plot.reals(prior.mat.1.5,"Matérn, nu = 1, var = 5")

# Matérn, nu = 3, var = 1:
prior.mat.3.1 <- get.prior.real("matern", sigma.sq=1, phi, nu=3)
plot.reals(prior.mat.3.1,"Matérn, nu = 3, var = 1")

# Matérn, nu = 3, var = 5: 
prior.mat.3.5 <- get.prior.real("matern", sigma.sq=5, phi, nu=3)
plot.reals(prior.mat.3.5,"Matérn, nu = 3, var = 5")

# Problem d) choosing one realization to use as our data points:
# we choose matérn with nu=1 and variance=5

observations <- c(prior.mat.1.5[1,10], prior.mat.1.5[1,25], prior.mat.1.5[1,30])
H <- matrix(data=0,nrow=3,ncol=50)
H[1,10] <- 1; H[2,25] <- 1; H[3,30] <- 1

covmatrix.r <- getCovMatrix(cov.spatial(0:49,cov.model="matern",cov.pars=c(5,10),kappa=1))

error.var.0 <- matrix(data=0, nrow=3, ncol=3)
error.var.025 <- 0.25*diag(3)

mu.error <- covmatrix.r%*%t(H)%*%inv(H%*%covmatrix.r%*%t(H) + error.var.025)%*%observations
mu.no.error <- covmatrix.r%*%t(H)%*%inv(H%*%covmatrix.r%*%t(H) + error.var.0)%*%observations

cov.rd.025 <- covmatrix.r - covmatrix.r%*%t(H)%*%inv(H%*%covmatrix.r%*%t(H) + error.var.025)%*%H%*%covmatrix.r
sd.025 <- sqrt(diag(cov.rd.025))

low.lim <- mu.error - qnorm(0.95)*sd.025
high.lim <- mu.error + qnorm(0.95)*sd.025


plot(0:49, mu.error, type="l",col="blue",ylab="^r(x|d)",xlab="x",ylim=c(min(low.lim),max(high.lim)))
lines(0:49, low.lim, lty="dashed",col="red")
lines(0:49, high.lim, lty="dashed",col="red")

cov.rd.0 <- covmatrix.r - covmatrix.r%*%t(H)%*%inv(H%*%covmatrix.r%*%t(H) + error.var.0)%*%H%*%covmatrix.r
sd.0 <- sqrt(abs(diag(cov.rd.0)))

low.lim.0 <- mu.no.error - qnorm(0.95)*sd.0
high.lim.0 <- mu.no.error + qnorm(0.95)*sd.0

plot(0:49, mu.no.error, type="l",col="blue",ylab="^r(x|d)",xlab="x",ylim=c(min(low.lim.0),max(high.lim.0)))
lines(0:49, low.lim.0, lty="dashed",col="red")
lines(0:49, high.lim.0, lty="dashed",col="red")







