library(RandomFields)
library(geoR)
library(mvtnorm)
library(ggplot2)
#a)

set.seed(123)
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
