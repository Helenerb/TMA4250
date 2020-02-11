library(fields); 
library(geoR);  
library(MASS)
#install.packages("akima")
library(akima)
library(ggplot2)

topo <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise1/topo.dat")
topo.geodata <- list("coords"=cbind(topo$x,topo$y),"data"=cbind(topo$z))
pred.points <- expand.grid(seq(1,315,by=1),seq(1,315,by=1))

krige.params <- krige.control(type.krige="ok",trend.d="2nd",trend.l = "2nd",cov.model="powered.exponential",cov.pars=c(2500,100),kappa=1.5)
krige.pred <- krige.conv(coords = topo.geodata$coords,data=topo.geodata$data,locations=pred.points,krige=krige.params)

##Plot: 

# make krige.pred into dataframe for easier plotting
krige.df <- data.frame("x"=pred.points$Var1,"y"=pred.points$Var2,"z"=krige.pred$predict, "v"=krige.pred$krige.var)

# plotting the Krige prediction 
krige.plot = ggplot() + geom_tile(data=krige.df, aes(x=x, y=y, fill=z)) + scale_fill_gradientn(colors=rainbow(4))
krige.plot = krige.plot + theme_classic()
krige.plot

# plotting the variance from the Krige prediction
krige.plot.var = ggplot() + geom_tile(data=krige.df, aes(x=x, y=y, fill=v)) + scale_fill_gradientn(colors=c("blue", "lightblue","lightgreen","green","orange","red"),breaks=c(0,200,400,600,800,Inf))
krige.plot.var = krige.plot.var + theme_classic()
krige.plot.var = krige.plot.var + geom_point(data=topo, aes(x=topo$x,y=topo$y),shape=4)
krige.plot.var


# 2e:

#prediction at [100,100]:
krige.pred$predict[100*315+100]

# variance at [100,100]:
krige.pred$krige.var[100*315+100]
