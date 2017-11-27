# R plotting script to go with the output of GA_solve
# to make use of the more flexible ggplot2.
library(dplyr)
library(ggplot2)
library(RColorBrewer)

data_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/data_out'

plot_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/plot_out'

##########################################################################
# a function to recode actions so they can all share the same colour
# pallett, but with distinct hues for each action 
recode_act = function(sub_act, act_code){

	recoded = numeric(length(sub_act))
	
	# herbicide is simply itself 
	recoded[sub_act == 'herbicide'] = act_code[sub_act == 'herbicide']

	# for crop add 4 (to shift 1 past herbicide)
	recoded[sub_act == 'crop'] = act_code[sub_act == 'crop'] + 4

	# for crop add 7 (to shift 1 past herbicide + crop)
	recoded[sub_act == 'plow'] = act_code[sub_act == 'plow'] + 7
	
	# for crop add 10 (to shift 1 past herbicide + crop + plow)
	recoded[sub_act == 'spot'] = act_code[sub_act == 'spot'] + 10

	return(recoded)
	
}

# make a colour pallett for actions that is consistant across plots
set_act_col = function(){

	no_act = grey(1)
	herb_cols = c(no_act, brewer.pal(n = 9, "Blues")[c(3, 6, 9)])
	crop_cols = brewer.pal(n = 9, "Greens")[c(3, 6, 9)]
	plow_cols = c(no_act, grey(0))
	spot_cols = c(no_act, grey(0))

	return(c(herb_cols, crop_cols, plow_cols, spot_cols))
}

# make a color legend for the action sequences over time
make_act_key = function(){

	act_pal = set_act_col()
	
	# set up the dummy dataframe
       df  = data.frame(sub_act = rep(c('herbicide', 'crop', 'plow', 
			'spot'), each = 12),
		dum_x = rep(1:12, 4),
		act_val = c(rep(1:4, each = 3), rep(5:7, each = 4), 
			rep(8:9, each = 6), rep(10:11, each = 6)))	

       text_size = 5
	
       key_plt = ggplot(df, aes(x = dum_x, y = sub_act)) + 
		geom_tile(aes(fill = factor(act_val))) + 
		scale_fill_manual(values = setNames(act_pal, 
			1:length(act_pal))) + 
		labs(x = '', y = '', title = 'actions') + xlim(0.5, 12.5) +  
		annotate('segment', x = 0.5, xend = 12.5, 
			y = c(1.5, 2.5, 3.5, 4.5), 
			yend = c(1.5, 2.5, 3.5, 4.5), color = grey(0.7)) + 
		annotate('text', x = c(2.5, 6.5, 10.5), y = c(1, 1, 1), 
			label = c('WW', 'Alt', 'Lay'), size = text_size,
			colour = c(grey(0), grey(0), grey(1))) + 
		annotate('text', x = c(2, 5, 8, 11), y = c(2, 2, 2, 2), 
			colour = c(grey(0), grey(0),grey(0), grey(1)),  
			label = c('none', 'H1', 'H2', 'H1 & H2'), size = text_size) +
		annotate('text', x = c(3.5, 9.5, 3.5, 9.2), y = c(3, 3, 4, 4), 
			colour = c(grey(0), grey(1), grey(0), grey(1)),  
			label = c('no', 'yes', 'no', 'yes'), size = text_size) + 
		theme(legend.position = 'none', 
		      	panel.background = element_rect(fill = 'white'),
			axis.text.x = element_blank(), 
			axis.text.y = element_text(size = 15),
			axis.ticks = element_blank(),
			plot.title = element_text(size = 20, hjust = 0.5))

	return(key_plt)

}

##########################################################################

# plot of best action sequence found 
setwd(data_loc)
ga_test = read.csv('test_GA.csv', header = TRUE, stringsAsFactors = FALSE)

ga_test = mutate(ga_test, act_reco = recode_act(sub_act, best_act))

# pallett that links with the recoded action codes
act_pal = set_act_col()

hm = ggplot(ga_test, aes(x = ts, y = sub_act)) + 
	geom_tile(aes(fill = factor(act_reco)), color = grey(0.85)) +
	scale_fill_manual(values = setNames(act_pal, 1:length(act_pal))) +
	labs(x = 'time step', y = '') + 
	theme(legend.position = 'none', 
		panel.background = element_rect(fill = 'white'),
		axis.text.y = element_text(size = 15),
		axis.text.x = element_text(size = 10),
		axis.title = element_text(size = 15)) 

##########################################################################
# plot some parameter sweeps 
setwd(data_loc)
dis_sweep = read.csv('dis_rate_sweep.csv', header = TRUE, stringsAsFactors = FALSE)

my_blues = brewer.pal(9, 'Blues')[c(5, 7, 9)]

# grab the mesaures and parameters that I want 
df_dis = select(dis_sweep, int_g1, int_g2, off_cv, Y0, dis_rate, proWW:scen)
df_dis = mutate(df_dis, 
	int_scen = ifelse(int_g1 < 0.2 & int_g2 < 0.2, 'int state: susceptible h1 and h2',
	ifelse(int_g1 > 0.2 & int_g2 < 0.2, 'int state: resistant h1, susceptible h2',
	'int state: resistant h1 and h2')),
	Y0_pretty = paste0('Wheat yeild: ', Y0, ' £/ha'))

# split up the dataframes to seperate different measures
# total reward as prportion of maximum avilable
rew_res = select(df_dis, int_g1:dis_rate, pro_max, scen:Y0_pretty) %>%
	filter(scen == 'optimGA')

rew_fix = select(df_dis, int_g1:dis_rate, pro_max, scen:Y0_pretty) %>%
	filter(scen != 'optimGA') %>% filter(scen != 'allalt')

rew_plt = ggplot(rew_res, aes(x = dis_rate, y = pro_max, 
		group = as.factor(off_cv))) + geom_point() +
	geom_line(aes(linetype = as.factor(off_cv))) + 
	geom_line(data = rew_fix, aes(x = dis_rate, y = pro_max, 
		linetype = as.factor(off_cv), colour = scen), 
		 inherit.aes = FALSE) +  
	scale_colour_manual(values = my_blues) +
	facet_grid(int_scen ~ Y0_pretty)

rew_plt

# resistance of most effective herbicide after 20 years
min_res = select(df_dis, int_g1:dis_rate, fin_min_res, scen:Y0_pretty) %>%
	filter(scen == 'optimGA')

res_fix = select(df_dis, int_g1:dis_rate, fin_min_res, scen:Y0_pretty) %>%
	filter(scen != 'optimGA') %>% filter(scen != 'allalt')

res_plt = ggplot(min_res, aes(x = dis_rate, y = fin_min_res, 
		group = as.factor(off_cv))) + geom_point() +
	geom_line(aes(linetype = as.factor(off_cv))) + 
	geom_line(data = res_fix, aes(x = dis_rate, y = fin_min_res, 
		linetype = as.factor(off_cv), colour = scen), 
		inherit.aes = FALSE) +  
	scale_colour_manual(values = my_blues) +
	facet_grid(int_scen ~ Y0_pretty)

res_plt

# total number of herbicide applications 
herb_app = select(df_dis, int_g1:dis_rate, herb_apps, scen:Y0_pretty) %>%
	filter(scen == 'optimGA')

herb_fix = select(df_dis, int_g1:dis_rate, herb_apps, scen:Y0_pretty) %>%
	filter(scen != 'optimGA') %>% filter(scen != 'allalt')

herb_plt = ggplot(herb_app, aes(x = dis_rate, y = herb_apps, 
		group = as.factor(off_cv))) + geom_point() +
	geom_line(aes(linetype = as.factor(off_cv))) + 
	facet_grid(int_scen ~ Y0_pretty)

herb_plt

# total spent on spot control 
spot = select(df_dis, int_g1:dis_rate, spot_spend, scen:Y0_pretty) %>%
	filter(scen == 'optimGA')

spot_plt = ggplot(spot, aes(x = dis_rate, y = spot_spend, 
		group = as.factor(off_cv))) + geom_point() +
	geom_line(aes(linetype = as.factor(off_cv))) + 
	facet_grid(int_scen ~ Y0_pretty)

spot_plt

# proportion of winter wheat in the rotation
WW_sum = select(df_dis, int_g1:dis_rate, proWW, scen:Y0_pretty) %>%
	filter(scen == 'optimGA')

WW_plt = ggplot(WW_sum, aes(x = dis_rate, y = proWW, 
		group = as.factor(off_cv))) + geom_point() +
	geom_line(aes(linetype = as.factor(off_cv))) + 
	facet_grid(int_scen ~ Y0_pretty)

WW_plt

# measures of SB size, showing tolerance for high populations early on
SB5_opt = select(df_dis, int_g1:dis_rate, SB5, scen:Y0_pretty) %>%
	filter(scen == 'optimGA')

SB5_plt = ggplot(SB5_opt, aes(x = dis_rate, y = SB5, 
		group = as.factor(off_cv))) + geom_point() +
	geom_line(aes(linetype = as.factor(off_cv))) + 
	facet_grid(int_scen ~ Y0_pretty)

SB5_plt

SB10_opt = select(df_dis, int_g1:dis_rate, SB10, scen:Y0_pretty) %>%
	filter(scen == 'optimGA')

SB10_plt = ggplot(SB10_opt, aes(x = dis_rate, y = SB10, 
		group = as.factor(off_cv))) + geom_point() +
	geom_line(aes(linetype = as.factor(off_cv))) + 
	facet_grid(int_scen ~ Y0_pretty)
SB10_plt

##########################################################################
# make a more focused set of plots to tell the story for the poster 
rew_df = select(df_dis, int_g1:dis_rate, pro_max, scen:Y0_pretty) %>%
	filter(scen == 'optimGA', Y0 == 1022, int_scen != 'int state: resistant h1 and h2') %>% 
	mutate(cv_pretty = ifelse(off_cv == 0, 'none', 'high'))

rew_fix = select(df_dis, int_g1:dis_rate, pro_max, scen:Y0_pretty) %>%
	filter(scen != 'optimGA') %>% 
	filter(scen != 'allalt', Y0 == 1022, int_scen != 'int state: resistant h1 and h2') %>%
	mutate(cv_pretty = ifelse(off_cv == 0, 'none', 'high'),
		pretty_scen = ifelse(scen == 'allh12', 'cont. h1 and h2',
			ifelse(scen == 'cych12', 'h1->h2->h1...', 
			       'h1->h2 & rotate crop')))


cur_con = mean(rew_df$Y0) * 20 

rew_plt = ggplot(rew_df, aes(x = dis_rate, y = (pro_max * cur_con) / 1000, 
		group = as.factor(off_cv))) + geom_point(size = 4) +
	geom_line(aes(linetype = cv_pretty), size = 1.5) + 
	scale_linetype_discrete(name = 'cross resistance') +
	labs(y = 'profit over 20 years (£ x1000)', x = 'discount rate') +
	geom_line(data = rew_fix, aes(x = dis_rate, 
		y = (pro_max * cur_con) / 1000, 
		linetype = cv_pretty, colour = pretty_scen), 
		inherit.aes = FALSE) +  
	scale_colour_manual(name = 'fixed strategy', values = my_blues) +
	theme(axis.title = element_text(size = 20),
	      axis.text = element_text(size = 20),
	      legend.background = element_rect(color = 'white'),
	      legend.text = element_text(size = 15),
	      legend.title = element_text(size = 15),
	      legend.position = c(0.8, 0.82), 
	      strip.text = element_text(size = 20)) +
	annotate('text', y = 10, x = 0.75, label = c('a)', 'b)'), size = 10) + 
	facet_grid(int_scen ~ .)

setwd(plot_loc)
pdf('reward_poster.pdf', height = 12, width = 8)

	rew_plt

dev.off()

# seed bank plot
SB10_df = select(df_dis, int_g1:dis_rate, SB10, scen:Y0_pretty) %>%
	filter(scen == 'optimGA', Y0 == 1022, int_scen != 'int state: resistant h1 and h2') %>% 
	mutate(cv_pretty = ifelse(off_cv == 0, 'none', 'high'))

SB10_plt = ggplot(SB10_df, aes(x = dis_rate, y = SB10 / 1000, 
		group = cv_pretty)) + geom_point(size = 4) +
	geom_line(aes(linetype = cv_pretty), size = 1.5) + 
	scale_linetype_discrete(name = 'cross resistance') +
	labs(y = 'average seed bank over first 10 years (x1000 seeds)', x = 'discount rate') +
	theme(axis.title = element_text(size = 20),
	      axis.text = element_text(size = 20),
	      legend.background = element_rect(color = 'white'),
	      legend.text = element_text(size = 20),
	      legend.title = element_text(size = 20),
	      legend.position = c(0.8, 0.82), 
	      strip.text = element_text(size = 20)) +
	annotate('text', y = 220, x = 0.75, label = c('a)', 'b)'), size = 10) + 
	facet_grid(int_scen ~ .)

setwd(plot_loc)
pdf('early_seedbank.pdf', height = 12, width = 8)

	SB10_plt

dev.off()

# resistance plot
res_df = select(df_dis, int_g1:dis_rate, fin_min_res, scen:Y0_pretty) %>%
	filter(scen == 'optimGA', Y0 == 1022, int_scen != 'int state: resistant h1 and h2') %>% 
	mutate(cv_pretty = ifelse(off_cv == 0, 'none', 'high'))

res_fix = select(df_dis, int_g1:dis_rate, fin_min_res, scen:Y0_pretty) %>%
	filter(scen != 'optimGA') %>% 
	filter(scen != 'allalt', Y0 == 1022, int_scen != 'int state: resistant h1 and h2') %>%
	mutate(cv_pretty = ifelse(off_cv == 0, 'none', 'high'),
		pretty_scen = ifelse(scen == 'allh12', 'cont. h1 and h2',
			ifelse(scen == 'cych12', 'h1->h2->h1...', 
			       'h1->h2 & rotate crop')))

res_plt = ggplot(res_df, aes(x = dis_rate, y = fin_min_res, 
		group = as.factor(off_cv))) + geom_point(size = 4) +
	geom_line(aes(linetype = cv_pretty), size = 1.5) + 
	scale_linetype_discrete(name = 'cross resistance') +
	labs(y = 'resistance to most effective herbicide (% survival)', x = 'discount rate') +
	geom_line(data = res_fix, aes(x = dis_rate, y = fin_min_res, 
		linetype = cv_pretty, colour = pretty_scen), 
		inherit.aes = FALSE) +  
	scale_colour_manual(name = 'fixed strategy', values = my_blues) +
	theme(axis.title = element_text(size = 20),
	      axis.text = element_text(size = 20),
	      legend.background = element_rect(color = 'white'),
	      legend.text = element_text(size = 15),
	      legend.title = element_text(size = 15),
	      legend.position = c(0.8, 0.3), 
	      strip.text = element_text(size = 20)) +
	annotate('text', y = Inf, x = 0.75, label = c('a)', 'b)'), 
		 size = 10, vjust = 1) + 
	facet_grid(int_scen ~ .)

setwd(plot_loc)
pdf('resistance_min_poster', height = 12, width = 8)

	res_plt

dev.off()






