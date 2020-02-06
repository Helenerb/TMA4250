library(ggplot2)
# sizes
n=100

# define regular grid of locations
sites1v=array((1:n),c(n,1))

# Prior mean
m=0
# compute East and North distances on grid
ww=array(1,c(n,1))
H=abs(sites1v%*%t(ww)-ww %*% t(sites1v))
# Exponential covariance model
#range=40 # correlation range
#Sigma=exp(-(3/range)*H)
# Matern 3/2
phiM=0.19
Sigma=(1+phiM*H)*exp(-phiM*H)
image(Sigma)
# Compute the Cholesky factor
L=chol(Sigma)
Li=solve(L)
image(Li)

# sample zero mean part
r=t(L)%*%rnorm(n)
# sample by adding mean
r=r+m
plot(r,type="l", xlim=c(0,100), ylim=c(-3,3))

#
# Data
M=5
F=matrix(0,M,n)
des=ceiling(n*runif(M))
for(i in 1:M)
{
  F[i,des[i]]=1
}

# measurement standard deviation
tau=0.05
# sample data
y=F%*% r + tau *rnorm(M)

#plot(y)

# Prediction surface from data
C=F %*%Sigma%*%t(F)+diag(tau^2,nrow=M,ncol=M)
mx=matrix(m,nrow=n,ncol=1) 
rhat=mx+Sigma %*% t(F) %*% solve(C,(y-F %*% mx));
# Prediction variances
Vvhat=Sigma-Sigma %*% t(F)%*% solve(C,(F %*%Sigma));
vr=diag(Vvhat)

rlow=rhat-1.28*sqrt(vr)
rupp=rhat+1.28*sqrt(vr)
df=data.frame(r,rhat,rlow,rupp,sites1v)
head(df)

ggplot(data=df, aes(x=sites1v, y=rhat)) +
  geom_line(aes(x=sites1v, y=rhat)) +
  geom_line(aes(x=sites1v, y=r),linetype=2,color=2) +
  geom_line(aes(x=sites1v, y = rlow),linetype=3,colour=3) +
  geom_line(aes(x=sites1v, y = rupp),linetype=3, colour=3) 
