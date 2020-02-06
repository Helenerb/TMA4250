library(fields); 
library(geoR);  
library(MASS)
library(akima)

distance=0:0.01:20
CorrelationFunction=Matern(distance , range = 2,nu=0.5) 
plot(distance,CorrelationFunction,type="l")
lines(distance,CorrelationFunction,type="l",col="blue")


