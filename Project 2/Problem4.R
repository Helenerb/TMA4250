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

p.Strauss <- function(X,Y, tau.0, phi.0, phi.1){
  return(tau.0)
}

Strauss <- function(tau, tau.0, phi.0, phi.1){
  return(0)
}