library(RandomFields)
library(geoR)
library(mvtnorm)
library(ggplot2)
#a)

set.seed(125)
n=30
L=seq(1:n)
L.grid=expand.grid(L,L)

var_r = 2; ksi_r = 3
cov.mat=varcov.spatial(coords=L.grid,dists.lowertri=NULL,cov.model="powered.exponential", kappa=1, cov.pars=c(var_r,3))
realizations.RF=as.vector(rmvnorm(n=1,mean=numeric(n*n),cov.mat$varcov,method="chol"))

df=data.frame("x"=L.grid$Var1, "y"=L.grid$Var2, "z"=realizations.RF)
ggplot() + geom_tile(data=df, aes(x=x, y=y, fill=z)) + theme_classic() + scale_fill_gradientn(colors=heat.colors(10))

#b
vgram.sampled=variog(coords=L.grid,data=df$z,uvec=13)
vgram.analytic =var_r*(1-exp(-vgram.sampled$u/ksi_r))
vgram = data.frame("sampled"=vgram.sampled$v, "analytic"=vgram.analytic, "dist"=vgram.sampled$u)
vgram <- reshape2::melt(vgram, id.var="dist")
ggplot(data=vgram, aes(x=dist, y=value, col=variable)) + geom_line() + theme_classic()

# d

#vgram.sampled=variog(coords=L.grid,data=df$z,uvec=13)
vgram.analytic =var_r*(1-exp(-vgram.sampled$u/ksi_r))
#vgram = data.frame("sampled"=vgram.sampled$v, "analytic"=vgram.analytic, "dist"=vgram.sampled$u)
#vgram <- reshape2::melt(vgram, id.var="dist")
#ggplot(data=vgram, aes(x=dist, y=value, col=variable)) + geom_line() + theme_classic() + ggtitle("seed = 125")

# sample 36 locations
sample.idx <- sort(as.integer(runif(36,1,900)))
sample.locs <- L.grid[sample.idx,]
sample.reals <- df[sample.idx,]

#vgram.sampled.36 = variog(coords=sample.locs, data=sample.reals$z, uvec=13)
#vgram.analytic.36 = var_r*(1-exp(-vgram.sampled.36$u/ksi_r))
#vgram.36 = data.frame("sampled"=vgram.sampled.36$v, "analytic"=vgram.analytic, "dist"=vgram.sampled.36$u)
#vgram.36 <- reshape2::melt(vgram.36, id.var="dist")
#ggplot(data=vgram.36, aes(x=dist, y=value, col=variable)) + geom_line() + theme_classic() + ggtitle("From 36 sampled values")

# considering sigma squared and xi to be unknown:

fitted.params = likfit(coords=L.grid, data=df$z, ini.cov.pars = c(2,1))
fitted.params
vgram.fitted = fitted.params$sigmasq*(1- exp(-seq(0:40)/fitted.params$phi))

vgram.analytic.d =var_r*(1-exp(-seq(0:40)/ksi_r))

vgram.d = data.frame("fitted"=vgram.fitted, "analytic"=vgram.analytic.d, "dist"=seq(0:40))
vgram.d <- reshape2::melt(vgram.d, id.var="dist")
ggplot(data=vgram.d, aes(x=dist, y=value, col=variable)) + geom_line() + theme_classic() + ggtitle("All observations used for ML estimation, estimated variance = 2.40 and estimated range = 3.91")

ML.variograms <- function(n,ini.cov.pars,title,L.grid,df,var_r,ksi_r){
  sample.idx <- sort(sample.int(900,n))
  sample.locs <- L.grid[sample.idx,]
  sample.reals <- df[sample.idx,]
  
  vgram.analytic.d =var_r*(1-exp(-seq(0:40)/ksi_r))
  
  fitted.params = likfit(coords=sample.locs, data=sample.reals$z, ini.cov.pars = ini.cov.pars)
  fitted.params
  
  vgram.fitted = fitted.params$sigmasq*(1- exp(-seq(0:40)/fitted.params$phi))
  
  vgram.d = data.frame("fitted"=vgram.fitted, "analytic"=vgram.analytic.d, "dist"=seq(0:40))
  vgram.d <- reshape2::melt(vgram.d, id.var="dist")
  ggplot(data=vgram.d, aes(x=dist, y=value, col=variable)) + geom_line() + theme_classic() + ggtitle(title) + geom_text(x=30,y=1,label=(fitted.params$sigmasq),show.legend = FALSE)
}

ML.variograms(36,c(2,1),"heii!!",L.grid,df,var_r,ksi_r)


