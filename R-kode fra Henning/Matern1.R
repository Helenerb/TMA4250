library(fields); 
library(geoR);  
library(MASS)
library(akima)

#HENNEKSEMPEL
distance=0:0.01:20
CorrelationFunction=Matern(distance , range = 2,nu=0.5) 
plot(distance,CorrelationFunction,type="l")
lines(distance,CorrelationFunction,type="l",col="blue")

#VÅRT
#skikkelig hofteskudd på parameteren phi
library(geoR)
CorrPexp=function(x,nu){
  return(exp(-x^nu))
}

curve(Matern(x,nu=1),from = 0, to = 5,
      xlab = expression(tau == (x_i - x_j)/10), ylab = expression(rho(tau)), lty = 2, col='blue')
curve(Matern(x,nu=3),from = 0, to = 5, add = TRUE, col='blue')
curve(CorrPexp(x,nu=1), from = 0, to = 5, add = TRUE, col ='red', lty=2)
curve(CorrPexp(x,nu=1.9), from = 0, to = 5, add = TRUE, col ='red')


legend("topright", expression("matern":nu==1, "matern":nu==3, "pexp":nu==1, "pexp":nu==1.9),
       lty=c(2,1,2,1), lwd=c(1,1,2,2), col=c("blue", "blue", "red", "red"))

### Eight different variograms
Vario=function(x,corrfunc,nu,var){
  return(var*(1-corrfunc(x,nu=nu)))
}
#sigma2 = 1
par(mfrow=c(1,2))
curve(Vario(x,Matern,nu=1, var = 1),from = 0, to = 5, main=expression(sigma^2==1),
      xlab = expression(tau == (x_i - x_j)/10), ylab = expression(gamma(tau)), lty = 2, col='blue')
curve(Vario(x,Matern,nu=3, var = 1),from = 0, to = 5, add=TRUE, col='blue')
curve(Vario(x,CorrPexp,nu=1, var = 1),from = 0, to = 5, add=TRUE, col='red')
curve(Vario(x,CorrPexp,nu=1.9, var = 1),from = 0, to = 5, add=TRUE, lty= 2, col='red')
legend("bottomright", expression(gamma=="mat(1)", gamma=="mat(3)", gamma=="pexp(1)", gamma=="pexp(1.9)"),
       lty=c(2,1,2,1), lwd=c(1,1,2,2), col=c("blue", "blue", "red", "red"), cex=0.85)
#sigma2=5
curve(Vario(x,Matern,nu=1, var = 5),from = 0, to = 5, main=expression(sigma^2==5),
      xlab = expression(tau == (x_i - x_j)/10), ylab = expression(gamma(tau)), lty = 2, col='blue')
curve(Vario(x,Matern,nu=3, var = 5),from = 0, to = 5, add=TRUE, col='blue')
curve(Vario(x,CorrPexp,nu=1, var = 5),from = 0, to = 5, add=TRUE, col='red')
curve(Vario(x,CorrPexp,nu=1.9, var = 5),from = 0, to = 5, add=TRUE, lty= 2, col='red')
legend("bottomright", expression(gamma=="mat(1)", gamma=="mat(3)", gamma=="pexp(1)", gamma=="pexp(1.9)"),
       lty=c(2,1,2,1), lwd=c(1,1,2,2), col=c("blue", "blue", "red", "red"), cex=0.85)




Mat.nu1.s5
Mat.nu3.s5





