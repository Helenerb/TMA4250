rm(list=ls())

library(tidyverse)
library(gridExtra)
set.seed(123)

# Load data
probs = read.table("obsprob.txt", header=TRUE)
pines = read.table("obspines.txt", header=TRUE)
domain = 300^2
grid.unit = 100



# Plot number of observed pines and observation probabilities

ggplot(pines, aes(x, y, fill=N_obs)) + geom_tile() + 
  scale_fill_gradient(expression('count d'), low='white',high='springgreen4') + theme_minimal() + 
  geom_hline(yintercept=seq(0,300,by=10), size=0.1) + 
  geom_vline(xintercept=seq(0,300,by=10), size=0.1)

ggplot(probs, aes(x, y, fill=alpha)) + geom_tile() + 
  scale_fill_gradient(expression('prob '*alpha), low='white',high='springgreen4') + theme_minimal() +
  geom_hline(yintercept=seq(0,300,by=10), size=0.1) + 
  geom_vline(xintercept=seq(0,300,by=10), size=0.1)

#estimate total number of trees
N.p = sum(pines$N_obs/probs$alpha)
lambda = N.p/domain

#six realizations
prior.reals=pines
prior.plots = list()
set.seed(123)
for(i in 1:6){
  prior.reals["N_obs"]=rpois(900,lambda*grid.unit)
  prior.plots[[i]]=ggplot(prior.reals, aes(x, y, fill=N_obs)) + geom_tile() + 
    scale_fill_gradient(expression('count d'), low='white',high='springgreen4', lim=c(0,7)) + theme_minimal() + 
    geom_hline(yintercept=seq(0,300,by=10), size=0.1) + 
    geom_vline(xintercept=seq(0,300,by=10), size=0.1)
}
 #line 44 needs to be run separately, for some reason ... 
do.call("grid.arrange", c(prior.plots, ncol=2))



#d
obs = pines$N_obs
N.grids=900; N.sims=100
lambda.pri = rep(lambda*grid.unit,N.grids); lambda.post = grid.unit*lambda*(1-probs$alpha)

#simulate and store
sims.pri <- t(mapply(function(lambda){rpois(N.sims, lambda)}, lambda=lambda.pri))
sims.post <- t(mapply(function(lambda){rpois(N.sims, lambda)}, lambda=lambda.post))

#calculate means for each node
means.pri=rowMeans(sims.pri) 
means.post=rowMeans(sims.post)
data.d = prior.reals=pines[c("x","y")]
data.d[c("pri","post")] = c(means.pri,means.post)
maks.val=max(data.d[c("pri","post")]); min.val=min(data.d[c("pri","post")])
         
#plot prior)d
ggplot(data=data.d, aes(x, y, fill=pri)) + geom_tile() + 
  scale_fill_gradient(expression('count d'), low='white',high='springgreen4', lim=c(min.val,maks.val)) + theme_minimal() + 
  geom_hline(yintercept=seq(0,300,by=10), size=0.1) +
  geom_vline(xintercept=seq(0,300,by=10), size=0.1) +
  ggtitle(paste(" Prior model \n", "Total count = ",round(sum(means.pri),0)))
#plot posterior
ggplot(data=data.d, aes(x, y, fill=post)) + geom_tile() + 
  scale_fill_gradient(expression('count d'), low='white',high='springgreen4', lim=c(min.val,maks.val)) + theme_minimal() + 
  geom_hline(yintercept=seq(0,300,by=10), size=0.1) + 
  geom_vline(xintercept=seq(0,300,by=10), size=0.1) +
  ggtitle(paste(" Posterior model \n", "Total count = ",round(sum(means.post),0)))














