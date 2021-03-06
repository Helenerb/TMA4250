#Load libraries
library(spatial)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)

#read data
cells   = read.table(file = 'cells.txt', col.names = c('x', 'y'))
pines   = read.table(file = 'pines.txt', col.names = c('x', 'y')) #first 3 lines not important?
redwood = read.table(file = 'redwood.txt', col.names = c('x', 'y'))

#initial plots
ggplot(data = cells, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Cells") +
  theme_bw()

ggplot(data = pines, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Pine trees") +
  theme_bw()

ggplot(data = redwood, aes(x=x, y=y)) +
  geom_point() +
  ggtitle("Redwood trees") +
  theme_bw()

#set up for Kfn
ppregion()
L_cells = Kfn(cells, fs = 1);      L_cells_df = data.frame(t = L_cells$x, L = L_cells$y)

L_redwood = Kfn(redwood, fs = 1) ; L_redwood_df = data.frame(t = L_redwood$x, L = L_redwood$y)

L_pines = Kfn(pines, fs = 1);      L_pines_df = data.frame(t = L_pines$x, L = L_pines$y)



#empirical and theoretical

ggplot() + 
  geom_line(data=L_cells_df, aes(x=t,y=L, colour='darkolivegreen4')) +
  geom_line(data=L_pines_df, aes(x=t,y=L, colour='springgreen4')) +
  geom_line(data=L_redwood_df, aes(x=t,y=L, colour='sienna')) +
  geom_abline(aes(slope = 1, intercept = 0, colour="black"),show.legend=FALSE)+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line()) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0))+
  geom_path() + 
  ggtitle("L-interaction functions")+
  scale_color_identity(name='',labels=c('darkolivegreen4'='cells','springgreen4'='pines','sienna'='redwood','black'='t'), guide="legend")

