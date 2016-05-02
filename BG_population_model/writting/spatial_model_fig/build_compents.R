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

layout_mat = rbind(c(0, 1, 0, 0.15), c(0, 1, 0.15, 0.3), c(0, 1, 0.3, 0.45),
  c(0, 1, 0.45, 0.5412), c(0, 1, 0.5412, 0.63333), c(0, 1, 0.63333, 0.752),
  c(0, 1, 0.752, 0.81667), c(0, 1, 0.81667, 0.90833), c(0, 1, 0.81667, 1))

## make a distributionf for the seed bank at each location
x_means = c(4, 0, -5, 28, 25, 20, 55, 50, 47)
RR_x1 = dnorm(seq(-10, 10, 0.1), x_means[1], 1.41) * 50
Rr_x1 = dnorm(seq(-10, 10, 0.1), x_means[2], 1.41) * 50
rr_x1 = dnorm(seq(-10, 10, 0.1), x_means[3], 1.41) * 100

RR_x2 = dnorm(seq(15, 35, 0.1), x_means[4], 1.41) * 55
Rr_x2 = dnorm(seq(15, 35, 0.1), x_means[5], 1.41) * 60
rr_x2 = dnorm(seq(15, 35, 0.1), x_means[6], 1.41) * 50

RR_x3 = dnorm(seq(40, 60, 0.1), x_means[7], 1.41) * 100
Rr_x3 = dnorm(seq(40, 60, 0.1), x_means[8], 1.41) * 50
rr_x3 = dnorm(seq(40, 60, 0.1), x_means[9], 1.41) * 50

# make the above ground distributions
RR_x1_ag = RR_x1 * 0.7
Rr_x1_ag = Rr_x1 * 0.7
rr_x1_ag = rr_x1 * (1 / (1 + exp(2 - 5 - pmin(5, seq(-10, 10, 0.1) * 1))))

RR_x2_ag = RR_x2 * 0.7
Rr_x2_ag = Rr_x2 * 0.7
rr_x2_ag = rr_x2 * (1 / (1 + exp(2 - 5 - pmin(5, seq(-10, 10, 0.1) * 1))))

RR_x3_ag = RR_x3 * 0.7
Rr_x3_ag = Rr_x3 * 0.7
rr_x3_ag = rr_x3 * (1 / (1 + exp(2 - 5 - pmin(5, seq(-10, 10, 0.1) * 1))))

# make the mixed distributions 
RR_RR_2_RR = dnorm(seq(15, 35, 0.1), 30, 1) * 100

RR_Rr_2_RR = dnorm(seq(15, 35, 0.1), 30, 1) * 50
RR_Rr_2_Rr = dnorm(seq(15, 35, 0.1), 25, 1) * 50

Rr_Rr_2_rr = dnorm(seq(15, 35, 0.1), 20, 1) * 25
Rr_Rr_2_Rr = dnorm(seq(15, 35, 0.1), 25, 1) * 50
Rr_Rr_2_RR = dnorm(seq(15, 35, 0.1), 30, 1) * 25

rr_RR_2_Rr = dnorm(seq(15, 35, 0.1), 25, 1) * 100

rr_Rr_2_rr = dnorm(seq(15, 35, 0.1), 20, 1) * 50
rr_Rr_2_Rr = dnorm(seq(15, 35, 0.1), 25, 1) * 50

rr_rr_2_rr = dnorm(seq(15, 35, 0.1), 20, 1) * 100

setwd(some place)
pdf(file = 'BG_pop_spatial_mod_schematic.pdf', height = 30, width = 10)
  split.screen(layout_mat)
  
  # single plot for all locations on one long axis
  screen(1)
    par(mar = c(4, 1, 1, 1))
    plot(-10:60, -10:60, ylim = c(0, max(c(RR_x1, Rr_x1, rr_x1, RR_x2, Rr_x2, rr_x2, RR_x3, Rr_x3, rr_x3))), type = 'n',
      bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n')
      
    polygon(seq(-10, 10, 0.1), rr_x1, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    polygon(seq(-10, 10, 0.1), RR_x1, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(-10, 10, 0.1), Rr_x1, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    text(x = x_means[1:3], y = c(max(RR_x1), max(Rr_x1), max(rr_x1)) / 2, labels = c('RR', 'Rr', 'rr'), cex = 2, 
      col = G_pallet)
      
    polygon(seq(15, 35, 0.1), rr_x2, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    polygon(seq(15, 35, 0.1), RR_x2, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(15, 35, 0.1), Rr_x2, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    text(x = x_means[4:6], y = c(max(RR_x2), max(Rr_x2), max(rr_x2)) / 2, labels = c('RR', 'Rr', 'rr'), cex = 2, 
      col = G_pallet)
      
    polygon(seq(40, 60, 0.1), rr_x3, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    polygon(seq(40, 60, 0.1), RR_x3, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(40, 60, 0.1), Rr_x3, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    mtext('g', side = 1, line = 1, cex = 2)
    text(x = x_means[7:9], y = c(max(RR_x3), max(Rr_x3), max(rr_x3)) / 2, labels = c('RR', 'Rr', 'rr'), cex = 2, 
      col = G_pallet)

  # above ground
  screen(2)
    par(mar = c(4, 1, 1, 1))
    plot(-10:60, -10:60, ylim = c(0, max(c(RR_x1, Rr_x1, rr_x1, RR_x2, Rr_x2, 
      rr_x2, RR_x3, Rr_x3, rr_x3))), type = 'n', bty = 'n', xlab = '', ylab = '', 
      yaxt = 'n', xaxt = 'n')
      
    polygon(seq(-10, 10, 0.1), rr_x1_ag, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    polygon(seq(-10, 10, 0.1), RR_x1_ag, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(-10, 10, 0.1), Rr_x1_ag, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
      
    polygon(seq(15, 35, 0.1), rr_x2_ag, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    polygon(seq(15, 35, 0.1), RR_x2_ag, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(15, 35, 0.1), Rr_x2_ag, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
      
    polygon(seq(40, 60, 0.1), rr_x3_ag, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    polygon(seq(40, 60, 0.1), RR_x3_ag, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(40, 60, 0.1), Rr_x3_ag, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)

    axis(side = 1, at = c(-20, 0, 25, 50, 100), labels = c('', expression(x[1]), expression(x[2]), bquote(x[3]), ''),
      cex.axis = 2, tck = 0.015)
    
    u <- par('usr')
    arrows(u[1], u[3], u[2], u[3], code = 2, xpd = TRUE)
    arrows(u[2], u[3], u[1], u[3], code = 2, xpd = TRUE)
 
  ## dispersal kernel
  screen(3)
    par(mar = c(4, 1, 1, 1))
    plot(-10:10, -10:10, type = 'n', ylim = c(0, 0.5), bty = 'n', axes = FALSE, xlab = 'distance', ylab = '', cex.lab = 2)
    axis(side = 1, pos = 0, at = c(-10, 10), label = c('', ''), tck = 0)
    axis(side = 2, pos = 0, at = c(0, 0.5), label = c('', ''), tck = 0)
    text(x = 0, y = 0.5, label = 'p(distance)', pos = 4, cex = 2)
    dists = c(rev(seq(0, 10, 0.1)), seq(0, 10, 0.1))
    seed_disp = 0.47 * exp(-0.8 * dists)
    ls_1D = c(seq(-10, 0, 0.1), seq(0, 10, 0.1))
    lines(ls_1D, seed_disp, lwd = 2, col = 'blue') 
    text(x = ls_1D[110], y = seed_disp[110], label = 'seed', col = 'blue', cex = 1.5, pos = 4)
    pollen_disp = 0.2 * exp(-0.3 * dists)
    lines(ls_1D, pollen_disp, lwd = 2, col = 'magenta') 
    text(x = ls_1D[110], y = pollen_disp[110] * 0.95, label = 'pollen', col = 'magenta', cex = 1.5, pos = 1)
    
  ## mixed seed 
  screen(4)
    par(mar = c(1, 1, 1, 1))
    plot(x = -10:60, y = -10:60, type = 'n', ylim = c(0, max(RR_RR_2_RR)), bty = 'n', axes = FALSE, xlab = '',
      ylab = '')
    polygon(seq(15, 35, 0.1), RR_RR_2_RR, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    text(x = 30, y = 0.5 * max(RR_RR_2_RR), labels = 'RR', cex = 2, col = G_pallet['RR'])
    text(x = c(15, 35), y = 0.5 * max(RR_RR_2_RR), labels = '...', cex = 4)
    text(x = 4, y = 0.5 * max(RR_RR_2_RR), labels = 'RR', cex = 3, col = G_pallet['RR'])
    text(x = 6, y = 0.5 * max(RR_RR_2_RR), labels = 'x', cex = 2)
    text(x = 8, y = 0.5 * max(RR_RR_2_RR), labels = 'RR', cex = 3, col = G_pallet['RR'])
    
  screen(5)
    par(mar = c(1, 1, 1, 1))
    plot(x = -10:60, y = -10:60, type = 'n', ylim = c(0, max(RR_RR_2_RR)), bty = 'n', axes = FALSE, xlab = '',
      ylab = '')
    polygon(seq(15, 35, 0.1), RR_Rr_2_RR, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    text(x = 30, y = 0.5 * max(RR_Rr_2_RR), labels = 'RR', cex = 2, col = G_pallet['RR'])
    polygon(seq(15, 35, 0.1), RR_Rr_2_Rr, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    text(x = 25, y = 0.5 * max(RR_Rr_2_Rr), labels = 'Rr', cex = 2, col = G_pallet['Rr'])
    
    text(x = c(15, 35), y = 0.5 * max(RR_RR_2_RR), labels = '...', cex = 4)
    text(x = 4, y = 0.5 * max(RR_RR_2_RR), labels = 'RR', cex = 3, col = G_pallet['RR'])
    text(x = 6, y = 0.5 * max(RR_RR_2_RR), labels = 'x', cex = 2)
    text(x = 8, y = 0.5 * max(RR_RR_2_RR), labels = 'Rr', cex = 3, col = G_pallet['Rr'])


  close.screen(all = TRUE)
dev.off()
