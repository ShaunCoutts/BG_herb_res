## Plotting script to visulise the population and optimisation outputs
using Plots 

pyplot();



# to plot population over 2d
heatmap(matrix)


# to animate a plot over time
f0s = convert(Array{Float64}, 1:10);
fr = 0.5
an_test = zeros(10, size(g1_vals)[1])

for i = 1:10 
  an_test[i, :] = 1 ./ fec_cost_maker(fr, f0s[i], g1_vals, g2_vals);
end

anim = @animate for i = 1:10
  
 plot(an_test[i,:])
 
end

gif(anim, "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/test.gif", 
  fps = 1)

