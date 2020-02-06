library(fields); 
library(geoR);  
library(MASS)
library(akima)
library(mvtnorm)

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

cov.vec <- cov.spatial(0:49,cov.model="matern",cov.pars=c(sigma[1],phi),kappa=nu.mat[1])

getCovMatrix <- function(cov.vec,n){
  cov.matrix <- matrix(data=0,nrow=50,ncol=50)
  for(i in 1:50){
    for(j in 1:50){
      cov.matrix[i,j] = cov.vec[abs(i-j)+1]
    }
  }
  return(cov.matrix)
}

prior.real <- rmvnorm(4,rep(0,50),getCovMatrix(cov.vec,50))

par(mfrow=c(2,2))

plot(d,prior.real[4,],type="l",col="orange",main="Maternal, nu=1")
lines(d,prior.real[1,],type="l",col="red")
lines(d,prior.real[2,],col="green")
lines(d,prior.real[3,],col="blue")

