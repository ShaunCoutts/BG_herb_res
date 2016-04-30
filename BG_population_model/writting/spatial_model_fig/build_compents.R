# script to make the components of a figure that explains how the spatial model works
library(colorspace)

#set up the pallets
G_pallet = rainbow_hcl(3, start = 180, end = 360)
G_pallet_fill = adjustcolor(G_pallet, alpha.f = 0.2)
names(G_pallet) <- c('RR', 'Rr', 'rr')
names(G_pallet_fill) <- c('RR', 'Rr', 'rr')

# set up the plotting areas, three plots along the bottom for seed bank distribution
# one wide one with three sets if three above ground distrubtions
# one wide dispersal kerenl plot with both seed and pollen
# then 9 mixed col dists at each of the three locations, possibly summarise with ellipsis.

layout_mat = rbind(c(), c(), c(),
  c(),
  c(),
  c(), c(), c(),
  c(), c(), c(),
  c(), c(), c(),
  c(), c(), c())

## make a distributionf for the seed bank at each location
x1_means = c(0, 4, -5)
RR_x1 = dnorm(seq(-10, 10, 0.1), x1_means[1], 1.41) * 50
Rr_x1 = dnorm(seq(-10, 10, 0.1), x1_means[2], 1.41) * 10
rr_x1 = dnorm(seq(-10, 10, 0.1), x1_means[3]5, 1.41) * 100

par(mar = c(4, 1, 1, 1))
plot(seq(-10, 10, 0.1), seq(-10, 10, 0.1), ylim = c(0, max(c(RR_x1, Rr_x1, rr_x1))), type = 'n',
  bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n')
polygon(seq(-10, 10, 0.1), rr_x1, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
polygon(seq(-10, 10, 0.1), RR_x1, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
polygon(seq(-10, 10, 0.1), Rr_x1, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
mtext('g', side = 1, line = 1, cex = 2)
text(x = x1_means, y = c(max(RR_x1), max(Rr_x1), max(rr_x1)) / 2, labels = c('RR', 'Rr', 'rr'), cex = 2, 
  col = G_pallet)

## above ground

## dispersal kernel

## mixed seed 


