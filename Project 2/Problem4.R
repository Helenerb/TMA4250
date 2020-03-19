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
  return(tau)
}

p.Strauss <- function(X,Y, tau.0, phi.0, phi.1){
  return(tau.0)
}

Strauss <- function(tau, tau.0, phi.0, phi.1){
  return(0)
}