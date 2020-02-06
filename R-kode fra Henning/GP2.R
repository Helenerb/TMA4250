# sizes
n1=25
n=n1*n1

# define regular grid of locations
vv=array((1:n1),c(n1,1))
uu=array(1,c(n1,1))
sites1=vv %*% t(uu)
sites2=uu %*% t(vv)
#plot(sites1,sites2)
# vector lengths of sites in East and North

sites1v=matrix(sites1,nrow=n,ncol=1)
sites2v=matrix(sites2,nrow=n,ncol=1)
# Plot sites
plot(sites1v,sites2v)

# Prior mean
m=0
# compute East and North distances on grid
ww=array(1,c(n,1))
ddE=sites1v%*%t(ww)-ww %*% t(sites1v)
dd2E=ddE*ddE
ddN=sites2v%*%t(ww)-ww %*% t(sites2v)
dd2N=ddN*ddN
H=sqrt(dd2N+dd2E)
# Exponential covariance model
range=40 # correlation range
Sigma=exp(-(3/range)*H)
image(Sigma)
# Compute the Cholesky factor
L=chol(Sigma)

# sample zero mean part
x=t(L)%*%rnorm(n)
# sample by adding mean
mx=matrix(m,nrow=n,ncol=1)  
x=x+mx
xm=matrix(x,nrow=n1,ncol=n1)
image(xm)

# Set a random design of 100 locations (could be copies)
M=20
des=ceiling(n*runif(M))
F=matrix(0,M,n)
for(i in 1:M)
{
  F[i,des[i]]=1
}
# measurement standard deviation
tau=0.05
# sample data
y=F%*% x + tau *rnorm(M)

#plot(y)

# Prediction surface from data
C=F %*%Sigma%*%t(F)+diag(tau^2,nrow=M,ncol=M)
xhat=mx+Sigma %*% t(F) %*% solve(C,(y-F %*% mx));
xhatm=matrix(xhat,nrow=n1,ncol=n1)
image(xhatm)

# Prediction variances
Vvhat=Sigma-Sigma %*% t(F)%*% solve(C,(F %*%Sigma));
vhat=diag(Vvhat)
vhatm=matrix(vhat,nrow=n1,ncol=n1)
image(vhatm)