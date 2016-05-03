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

layout_mat = rbind(c(0, 1, 0, 0.2), c(0, 1, 0.20, 0.4), c(0, 1, 0.4, 0.55),
  c(0, 1, 0.55, 0.625), c(0, 1, 0.625, 0.7), c(0, 1, 0.7, 0.775),
c(0, 1, 0.775, 0.85), c(0, 1, 0.85, 0.925), c(0, 1, 0.925, 1))

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

setwd('/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/writting/spatial_model_fig')
pdf(file = 'BG_pop_spatial_mod_schematic.pdf', height = 15, width = 10)
  split.screen(layout_mat)
  
  # single plot for all locations on one long axis
  screen(1)
    par(mar = c(4, 1, 1.5, 1))
    plot(-10:60, -10:60, ylim = c(0, max(c(RR_x1, Rr_x1, rr_x1, RR_x2, Rr_x2, rr_x2, RR_x3, Rr_x3, rr_x3))), type = 'n',
      bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n')
      
    polygon(seq(-10, 10, 0.1), rr_x1, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    polygon(seq(-10, 10, 0.1), RR_x1, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(-10, 10, 0.1), Rr_x1, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    text(x = x_means[1:3], y = c(max(RR_x1), max(Rr_x1), max(rr_x1)) / 2, labels = c('RR', 'Rr', 'rr'), cex = 1.5, 
      col = G_pallet)
      
    polygon(seq(15, 35, 0.1), rr_x2, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    polygon(seq(15, 35, 0.1), RR_x2, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(15, 35, 0.1), Rr_x2, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    text(x = x_means[4:6], y = c(max(RR_x2), max(Rr_x2), max(rr_x2)) / 2, labels = c('RR', 'Rr', 'rr'), cex = 1.5, 
      col = G_pallet)
      
    polygon(seq(40, 60, 0.1), rr_x3, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    polygon(seq(40, 60, 0.1), RR_x3, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(40, 60, 0.1), Rr_x3, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    mtext('metabolic reistance score (g)', side = 1, line = 2.5, cex = 1.5)
    text(x = x_means[7:9], y = c(max(RR_x3), max(Rr_x3), max(rr_x3)) / 2, labels = c('RR', 'Rr', 'rr'), cex = 1.5, 
      col = G_pallet)
    axis(side = 1, at = c(-10, 0, 10), labels = c(-10, 0, 10), tck = 0.015)
    axis(side = 1, at = c(15, 25, 35), labels = c(-10, 0, 10), tck = 0.015)
    axis(side = 1, at = c(40, 50, 60), labels = c(-10, 0, 10), tck = 0.015)
    
    mtext(expression(paste('seed bank at each location: ', italic('b(g, G, x, t)'))), side = 3, adj = 0, cex = 1.5)

  # above ground
  screen(2)
    par(mar = c(4, 1, 1, 1))
    plot(-10:60, -10:60, ylim = c(0, max(c(RR_x1, Rr_x1, rr_x1, RR_x2, Rr_x2, 
      rr_x2, RR_x3, Rr_x3, rr_x3))), type = 'n', bty = 'n', xlab = '', ylab = '', 
      yaxt = 'n', xaxt = 'n')
    mtext('location', side = 1, line = 2, cex = 1.5)
    
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
      cex.axis = 1.5, tck = 0.015)
    
    u <- par('usr')
    arrows(u[1], u[3], u[2], u[3], code = 2, xpd = TRUE)
    arrows(u[2], u[3], u[1], u[3], code = 2, xpd = TRUE)
 
    mtext(expression(paste('above ground population: ', italic('n(g, G, x, t)s(g, G, x, h')[italic(x)], italic(')'))), 
      side = 3, adj = 0, cex = 1.5)   
  ## dispersal kernel
  screen(3)
    par(mar = c(4, 1, 2, 1))
    plot(-10:10, -10:10, type = 'n', ylim = c(0, 0.5), bty = 'n', axes = FALSE, xlab = '', ylab = '', cex.lab = 1.5)
    axis(side = 1, pos = 0, at = c(-10, 10), label = c('', ''), tck = 0)
    axis(side = 2, pos = 0, at = c(0, 0.5), label = c('', ''), tck = 0)
    dists = c(rev(seq(0, 10, 0.1)), seq(0, 10, 0.1))
    seed_disp = 0.47 * exp(-0.8 * dists)
    ls_1D = c(seq(-10, 0, 0.1), seq(0, 10, 0.1))
    lines(ls_1D, seed_disp, lwd = 2, col = 'blue') 
    text(x = ls_1D[110], y = seed_disp[110], label = 'seed', col = 'blue', cex = 1, pos = 4)
    pollen_disp = 0.2 * exp(-0.3 * dists)
    lines(ls_1D, pollen_disp, lwd = 2, col = 'magenta') 
    text(x = ls_1D[110], y = pollen_disp[110] * 0.90, label = 'pollen', col = 'magenta', cex = 1, pos = 1)
    mtext('distance', side = 1, line = 1, cex = 1.5)

    mtext(expression(paste('seed and pollen dispersal kernel: ', italic('d(i, j)')[italic('m')], ' and ', italic('d(i,j)')['p'])), 
      side = 3, adj = 0, cex = 1.5)   
  ## mixed seed 
  screen(4)
    par(mar = c(1, 1, 1, 1))
    plot(x = -10:60, y = -10:60, type = 'n', ylim = c(0, max(RR_RR_2_RR)), bty = 'n', axes = FALSE, xlab = '',
      ylab = '')
    polygon(seq(15, 35, 0.1), RR_RR_2_RR, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    text(x = 7, y = 0.5 * max(RR_RR_2_RR), labels = 'RR', cex = 2, col = G_pallet['RR'])
    text(x = 9.5, y = 0.5 * max(RR_RR_2_RR), labels = 'x', cex = 1.5)
    text(x = 12, y = 0.5 * max(RR_RR_2_RR), labels = 'RR', cex = 2, col = G_pallet['RR'])
    
  screen(5)
    par(mar = c(1, 1, 1, 1))
    plot(x = -10:60, y = -10:60, type = 'n', ylim = c(0, max(RR_RR_2_RR)), bty = 'n', axes = FALSE, xlab = '',
      ylab = '')
    polygon(seq(15, 35, 0.1), RR_Rr_2_RR, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(15, 35, 0.1), RR_Rr_2_Rr, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    text(x = 7, y = 0.5 * max(RR_RR_2_RR), labels = 'RR', cex = 2, col = G_pallet['RR'])
    text(x = 9.8, y = 0.5 * max(RR_RR_2_RR), labels = 'x', cex = 1.5)
    text(x = 12, y = 0.5 * max(RR_RR_2_RR), labels = 'Rr', cex = 2, col = G_pallet['Rr'])

  screen(6)
    par(mar = c(1, 1, 1, 1))
    plot(x = -10:60, y = -10:60, type = 'n', ylim = c(0, max(RR_RR_2_RR)), bty = 'n', axes = FALSE, xlab = '',
      ylab = '')
    polygon(seq(15, 35, 0.1), Rr_Rr_2_RR, col = G_pallet_fill['RR'], border = G_pallet['RR'], lwd = 3)
    polygon(seq(15, 35, 0.1), Rr_Rr_2_Rr, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    polygon(seq(15, 35, 0.1), Rr_Rr_2_rr, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    text(x = 8, y = 0.5 * max(RR_RR_2_RR), labels = 'Rr', cex = 2, col = G_pallet['Rr'])
    text(x = 10, y = 0.5 * max(RR_RR_2_RR), labels = 'x', cex = 1.5)
    text(x = 12, y = 0.5 * max(RR_RR_2_RR), labels = 'Rr', cex = 2, col = G_pallet['Rr'])
    
  screen(7)
    par(mar = c(1, 1, 1, 1))
    plot(x = -10:60, y = -10:60, type = 'n', ylim = c(0, max(RR_RR_2_RR)), bty = 'n', axes = FALSE, xlab = '',
      ylab = '')
    polygon(seq(15, 35, 0.1), rr_RR_2_Rr, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    text(x = 8, y = 0.5 * max(RR_RR_2_RR), labels = 'rr', cex = 2, col = G_pallet['rr'])
    text(x = 9.5, y = 0.5 * max(RR_RR_2_RR), labels = 'x', cex = 1.5)
    text(x = 12, y = 0.5 * max(RR_RR_2_RR), labels = 'RR', cex = 2, col = G_pallet['RR'])
    
  screen(8)
    par(mar = c(1, 1, 1, 1))
    plot(x = -10:60, y = -10:60, type = 'n', ylim = c(0, max(RR_RR_2_RR)), bty = 'n', axes = FALSE, xlab = '',
      ylab = '')
    polygon(seq(15, 35, 0.1), rr_Rr_2_Rr, col = G_pallet_fill['Rr'], border = G_pallet['Rr'], lwd = 3)
    polygon(seq(15, 35, 0.1), rr_Rr_2_rr, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    text(x = 8, y = 0.5 * max(RR_RR_2_RR), labels = 'rr', cex = 2, col = G_pallet['rr'])
    text(x = 10, y = 0.5 * max(RR_RR_2_RR), labels = 'x', cex = 1.5)
    text(x = 12, y = 0.5 * max(RR_RR_2_RR), labels = 'Rr', cex = 2, col = G_pallet['Rr'])
    
  screen(9)
    par(mar = c(1, 1, 1, 1))
    plot(x = -10:60, y = -10:60, type = 'n', ylim = c(0, max(RR_RR_2_RR)), bty = 'n', axes = FALSE, xlab = '',
      ylab = '')
    polygon(seq(15, 35, 0.1), rr_rr_2_rr, col = G_pallet_fill['rr'], border = G_pallet['rr'], lwd = 3)
    text(x = 8, y = 0.5 * max(RR_RR_2_RR), labels = 'rr', cex = 2, col = G_pallet['rr'])
    text(x = 10, y = 0.5 * max(RR_RR_2_RR), labels = 'x', cex = 1.5)
    text(x = 12, y = 0.5 * max(RR_RR_2_RR), labels = 'rr', cex = 2, col = G_pallet['rr'])
    
    mtext(expression(paste('offspring distribution: ', italic('f(g, G, h, x, t)'))), 
      side = 3, adj = 0, cex = 1.5, line = -1)   
    mtext('Offspring distributions of each\ngenotype G, produced by\neach combination of parent\ngenotypes G x G.',
      side = 3, adj = 1, cex = 1.3, line = -10)   
 
  close.screen(all = TRUE)
dev.off()
