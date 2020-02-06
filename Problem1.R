library(fields); 
library(geoR);  
library(MASS)
library(akima)

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

