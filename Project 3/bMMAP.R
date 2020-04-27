library(reshape2)
library(ggplot2)
library(spatial)
source("b.R")

seismic <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/seismic.dat")
seismic <- seismic[,1]
seismic.mat <- matrix(seismic, nrow = 75, ncol = 75)
seismic.melt <- melt(seismic.mat)

pi1 <- p.i1(seismic.mat, 75, 75)

indicator <- function(x){
  if(x >= 0.5){
    return(1)
  }
  else{
    return(0)
  }
} 

predict <- sapply(pi1$mean, FUN = indicator)
predict <- data.frame(Var1 = seismic.melt$Var1, Var2 = seismic.melt$Var2, value = predict)

# plot analytical prediction 
pred.gg <- ggplot(data = predict, aes(x = Var1, y=Var2, fill = value)) + geom_tile()
pred.gg <- pred.gg + xlab(label="") + ylab(label="") 
pred.gg <- pred.gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank())
pred.gg

# by simulation:
S = 1000
sims <- matrix(nrow=S, ncol=length(predict$value))
# sett inn observasjon på sims[s,]
for(s in 1:S){
  obs <- posterior.b(seismic.mat, 75, 75)
  obs <- melt(obs)
  sims[s,] <- obs$value
}
sims.mean <- apply(sims,2,mean)
sims.pred <- sapply(sims.mean,FUN = indicator)
sims.pred <- data.frame(Var1 = seismic.melt$Var1, Var2 = seismic.melt$Var2, value = sims.pred)

sims.pred.gg <- ggplot(data = sims.pred, aes(x = Var1, y=Var2, fill = value)) + geom_tile()
sims.pred.gg <- sims.pred.gg + xlab(label="") + ylab(label="") 
sims.pred.gg <- sims.pred.gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank())
sims.pred.gg
