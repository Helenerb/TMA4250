# code for b
library(reshape2)
library(ggplot2)
library(spatial)

seismic <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/seismic.dat")
seismic <- seismic[,1]
seismic.mat <- matrix(seismic, nrow = 75, ncol = 75)
seismic.melt <- melt(seismic.mat)

posterior.b <- function(data, nrow, ncol){
  l <- matrix(nrow=nrow, ncol=ncol)
  for(i in 1:nrow){
    for(j in 1:ncol){
      d <- data[i,j]
      p <- exp(-(d - 0.08)^2/0.0072)/(exp(-(d - 0.08)^2/0.0072) + exp(-(d - 0.02)^2/0.0072))
      if(p > runif(1)){
        l[i,j] <- 1
      }
      else{
        l[i,j] <- 0
      }
    }
  }
  return(l)
}

post <- posterior.b(seismic.mat, 75, 75)
post.melt <- melt(post)

# plot realizations:
plot.gg <- ggplot(data = post.melt, aes(x = Var1, y=Var2, fill = value)) + geom_tile()
plot.gg <- plot.gg + xlab(label="") + ylab(label="") 
plot.gg <- plot.gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank())
plot.gg

# function to find expectation and variance:
p.i1 <- function(data, nrow, ncol){
  p <- matrix(nrow=nrow, ncol=ncol)
  var <- matrix(nrow=nrow, ncol=ncol)
  for(i in 1:nrow){
    for(j in 1:ncol){
      d <- data[i,j]
      p[i,j] <- exp(-(d - 0.08)^2/0.0072)/(exp(-(d - 0.08)^2/0.0072) + exp(-(d - 0.02)^2/0.0072))
      var[i,j] <- p[i,j]*(1 - p[i,j])
    }
  }
  # reformatting for easier access:
  p.melt <- melt(p)
  var.melt <- melt(var)
  res <- data.frame(Var1 = p.melt$Var1, Var2=p.melt$Var2, mean = p.melt$value, var = var.melt$value)
  return(res)
}

ExpVar <- p.i1(seismic.mat, 75, 75)

# expextation:
exp.gg <- ggplot(data = ExpVar, aes(x = Var1, y=Var2, fill = mean)) + geom_tile()
exp.gg <- exp.gg + xlab(label="") + ylab(label="") 
exp.gg <- exp.gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank())
exp.gg

# variance: 
var.gg <- ggplot(data = ExpVar, aes(x = Var1, y=Var2, fill = var)) + geom_tile()
var.gg <- var.gg + xlab(label="") + ylab(label="") 
var.gg <- var.gg + scale_fill_distiller(palette="Greens")
var.gg <- var.gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank())
var.gg



