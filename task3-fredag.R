library(RandomFields)
library(geoR)
library(mvtnorm)
library(ggplot2)
#a)

set.seed(124)
n=30  # grid size
L=seq(1:n)
L.grid=expand.grid(L,L)

var_r = 2; ksi_r = 3  # global variables for marameters
cov.mat=varcov.spatial(coords=L.grid,dists.lowertri=NULL,cov.model="powered.exponential", kappa=1, cov.pars=c(var_r,3))
realizations.RF=as.vector(rmvnorm(n=1,mean=numeric(n*n),cov.mat$varcov,method="chol"))

df=data.frame("x"=L.grid$Var1, "y"=L.grid$Var2, "z"=realizations.RF)
ggplot() + geom_tile(data=df, aes(x=x, y=y, fill=z)) + theme_classic() + scale_fill_gradientn(colors=heat.colors(10))

#b
vgram.sampled=variog(coords=L.grid,data=df$z,uvec=13)
vgram.analytic =var_r*(1-exp(-vgram.sampled$u/ksi_r))
vgram = data.frame("sampled"=vgram.sampled$v, "analytic"=vgram.analytic, "dist"=vgram.sampled$u)
vgram <- reshape2::melt(vgram, id.var="dist")
ggplot(data=vgram, aes(x=dist, y=value, col=variable)) + geom_line() + theme_classic() + ggtitle(label="seed = 124")

# c): change seed to 124 and 125 and run lines 20-24, from task b)

# d

# sample 36 locations
sample.idx <- sort(sample.int(900,36))
sample.locs <- L.grid[sample.idx,]
sample.reals <- df[sample.idx,]  # get realization at sampled points

vgram.sampled.36 = variog(coords=sample.locs, data=sample.reals$z, uvec=13)  # creating variogram from 36 samples
vgram.analytic.36 = var_r*(1-exp(-vgram.sampled.36$u/ksi_r))
vgram.36 <- data.frame("sampled"=vgram.sampled.36$v, "analytic"=vgram.analytic, "dist"=vgram.sampled.36$u)
vgram.36 <- reshape2::melt(vgram.36, id.var="dist")
ggplot(data=vgram.36, aes(x=dist, y=value, col=variable)) + geom_line() + theme_classic() + ggtitle("From 36 sampled values")

# considering sigma squared and xi to be unknown:

ML.variograms <- function(n,ini.cov.pars,title,L.grid,df,var_r,ksi_r){
  sample.idx <- sort(sample.int(900,n))  # sample n random locations
  sample.locs <- L.grid[sample.idx,]
  sample.reals <- df[sample.idx,]
  
  vgram.analytic.d =var_r*(1-exp(-seq(0:40)/ksi_r))  #defining the analytic variogram
  
  fitted.params = likfit(coords=sample.locs, data=sample.reals$z, ini.cov.pars = ini.cov.pars)
  fitted.params
  
  vgram.fitted = fitted.params$sigmasq*(1- exp(-seq(0:40)/fitted.params$phi))  # using the analytical expression with fitted parameters to make variogram
  
  vgram.d = data.frame("fitted"=vgram.fitted, "analytic"=vgram.analytic.d, "dist"=seq(0:40))
  vgram.d <- reshape2::melt(vgram.d, id.var="dist")
  plot <- ggplot(data=vgram.d, aes(x=dist, y=value, col=variable)) + geom_line() + theme_classic() + ggtitle(title) + geom_text(x = 25, y=1, label=expression(sigma^2), show.legend = FALSE) + geom_text(x=30,y=1,label=(round(fitted.params$sigmasq,2)),show.legend = FALSE) + geom_text(x = 25, y=0.7, label=expression(xi^2), show.legend=FALSE) + geom_text(x=30,y=0.7,label=(round(fitted.params$phi,2)),show.legend = FALSE)
  return(plot)
}

require(gridExtra)

# plots for 3d)
plot9 = ML.variograms(9,c(2,1),"Using 9 observations",L.grid,df,var_r,ksi_r)
plot64 = ML.variograms(64,c(2,1),"Using 64 observations",L.grid,df,var_r,ksi_r)
plot100 = ML.variograms(100,c(2,1),"Using 100 observations",L.grid,df,var_r,ksi_r)

grid.arrange(plot9,plot64,plot100,ncol=3)


# plots for 3e)
plotall = ML.variograms(900,c(2,1),"Using all observations",L.grid,df,var_r,ksi_r)
plot36 = ML.variograms(36,c(2,1),"Using 36 observations",L.grid,df,var_r,ksi_r)

grid.arrange(plotall,plot36,ncol=2)


