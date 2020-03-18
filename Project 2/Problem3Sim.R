# NS simulation
library(matrixStats)

ppregion()
L.redwood <- Kfn(redwood, fs = sqrt(2))

NS.sim <- function(lamb.M, sigma.c, p.mu, p.sigma, S){
  # initate matrix to store simulated values
  L.x <- matrix(0L, nrow=0, ncol=50)                  # is 50 general?
  L.y <- matrix(0L, nrow=0, ncol=50)
  
  for(s in 1:S){
    NS <- NS(lamb.M, sigma.c, p.mu, p.sigma)
    
    x <- NS$x.C.all[,"x"]
    y <- NS$x.C.all[,"y"]
    
    NS.df <- data.frame("x" = x, "y" = y)
    L.NS <- Kfn(NS.df, fs = sqrt(2))
    L.x <- rbind(L.x, L.NS$x)
    L.y <- rbind(L.y, L.NS$y)
    
  }
  L.mean = colMeans(L.y)
  L.min = colMins(L.y)
  L.max = colMaxs(L.y)
  
  return(list("L.x" = L.x, "L.y" = L.y, "mean" = L.mean, "min" = L.min, "max" = L.max ))
}

# first guestimate:
#NS.sim <- NS.sim(7, 0.125^2, 8, 10, 100)

# second guestimate:
#NS.sim <- NS.sim(7, 0.2^2, 8, 10, 100)

# third guestimate:
#NS.sim <- NS.sim(7, 0.125^2, 6, 6, 100)

# fourth guestimate:
NS.sim <- NS.sim(7, 0.1^2, 6, 6, 100)


NS.sim$L.y
NS.sim$L.x
NS.sim$min
NS.sim$max
NS.sim$mean

# 95% confidence interval:
NS.lower <- NS.sim$mean + 0.95*(NS.sim$max - NS.sim$mean)
NS.upper <- NS.sim$mean - 0.95*(NS.sim$mean - NS.sim$min)

# plot confidence interval together with redwood:

#collect in dataframe for ggplot:

gg.df <- data.frame("redwood" = L.redwood$y, "x" = L.redwood$x, "mean" = NS.sim$mean,
                    "lower" = NS.lower, "upper" = NS.upper)

NS.plot <- ggplot(data=gg.df) + geom_point(aes(x=x, y=redwood), color="sienna") + 
  geom_line(aes(x=x, y=upper), color="red", size=0.75) + 
  geom_line(aes(x=x, y=lower), color="red", size=0.75) + 
  xlab(label="t") + ylab(label="L")
NS.plot










