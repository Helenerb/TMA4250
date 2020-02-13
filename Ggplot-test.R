library(ggplot2)

mtcars

basic <- ggplot(data, aes(wt, mpg, colour = factor(cyl), shape = factor(vs) )) +
  geom_point()
basic