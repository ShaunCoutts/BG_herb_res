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






