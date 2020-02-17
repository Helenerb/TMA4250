library(fields); 
library(geoR);  
library(MASS)
#install.packages("akima")
library(akima)
library(ggplot2)
library(scales)

set.seed(123)

# 2a) 

#  loading data: 
sample <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise1/topo.dat")

#linear interpolation
plot(y ~ x, data = sample)
with(sample, text(x, y, formatC(z,dig=2), adj = -0.1))
sample.smooth <-
  with(sample, interp(x, y, z, nx=100, ny=100))
si.zmin <- min(sample.smooth$z,na.rm=TRUE)
si.zmax <- max(sample.smooth$z,na.rm=TRUE)
breaks <- pretty(c(si.zmin,si.zmax),10)
colors <- heat.colors(length(breaks)-1)

png(filename="krigdisp1.png")
image  (sample.smooth,
        breaks=breaks, col=colors)
contour(sample.smooth, add = TRUE, levels=breaks, col = "thistle")
points(sample, pch = 3, cex = 2, col = "blue")
dev.off()

## spline interpolation
sample.spl <- with(sample, interp(x, y, z, nx=100, ny=100, linear=FALSE))

png(filename="krigdisp2.png")  # save image
contour(sample.spl)
points(sample)
dev.off()

full.pal <- function(n) hcl.colors(n, "Oslo")

png(filename="krigdisp3.png")  # save image

filled.contour(sample.spl, color.palette = full.pal,
               plot.axes = { axis(1); axis(2);
                 points(sample, pch = 3, col= hcl(c=100, l = 20))})

dev.off()

# 2c) and 2d)

# argument trd either cte in c) or 2nd in d)
krige.predict <- function(trd, title){
  topo <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise1/topo.dat")
  topo.geodata <- list("coords"=cbind(topo$x,topo$y),"data"=cbind(topo$z))
  pred.points <- expand.grid(seq(1,315,by=1),seq(1,315,by=1))
  
  krige.params <- krige.control(type.krige="ok",trend.d=trd,trend.l = trd,cov.model="powered.exponential",cov.pars=c(2500,100),kappa=1.5)
  krige.pred <- krige.conv(coords = topo.geodata$coords,data=topo.geodata$data,locations=pred.points,krige=krige.params)
  cat(max(krige.pred$krige.var), "\n")
  cat(mean(krige.pred$krige.var))
  
  ##Plot: 
  
  # make krige.pred into dataframe for easier plotting
  krige.df <- data.frame("x"=pred.points$Var1,"y"=pred.points$Var2,"z"=krige.pred$predict, "v"=krige.pred$krige.var)
  
  # Comment in or out for plot of prediction or plot of variance
  
  # plotting the Krige prediction 
  krige.plot = ggplot() + geom_tile(data=krige.df, aes(x=x, y=y, fill=z)) + scale_fill_gradientn(colors=c("deepskyblue4", "darkturquoise","darkseagreen1","darkgoldenrod1","chocolate","brown4"),values=c(0,0.2,0.4,0.6,0.8,1), limits=c(650,1000))
  krige.plot = krige.plot + theme_classic()
  krige.plot = krige.plot + geom_contour(data=krige.df, aes(x=x, y=y, z=z), colour="black") + ggtitle(label=title)
  krige.plot
  
  # plotting the variance from the Krige prediction
  #krige.plot.var = ggplot() + geom_tile(data=krige.df, aes(x=x, y=y, fill=v)) + scale_fill_gradientn(colors=c("deepskyblue4", "darkturquoise","darkseagreen1","darkgoldenrod1","chocolate","brown4"),values=c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1000)) + ggtitle(label=title)
  #krige.plot.var = krige.plot.var + theme_classic()
  #krige.plot.var = krige.plot.var + geom_point(data=topo, aes(x=topo$x,y=topo$y),shape=4)
  #krige.plot.var
}

krige.predict("2nd","g(x) second degree")
krige.predict("cte","g(x) constant")
