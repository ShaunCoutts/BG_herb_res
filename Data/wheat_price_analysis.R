# simple script to take wheat prices over each year, get a yearly national average with 
#  So what does the average farmer in the average year get
# then the extrams of the best and worst a farmer could have gotten over the last 10 years 

#clean up data
setwd('/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/Data')
price_data = read.csv('wheat_price_2006_2016.csv', header = TRUE, stringsAsFactors = FALSE)
price_data = price_data[!is.na(price_data$Breadmaking_Wheat), ]
#look at market names
unique(price_data$Region)

#get mean per year along with min and max
price_mean = tapply(price_data$Breadmaking_Wheat, INDEX = price_data$Delv_Year, FUN = mean)
price_min = tapply(price_data$Breadmaking_Wheat, INDEX = price_data$Delv_Year, FUN = min)
price_max = tapply(price_data$Breadmaking_Wheat, INDEX = price_data$Delv_Year, FUN = max)

#should probably adjust for inflation as the low year that occurs in 2006 might not be the real low year and some of the variation seen is just inflation
# CPI data from https://www.ons.gov.uk/economy/inflationandpriceindices/timeseries/d7bt down loaded 23/3/2016
CPI_data = read.csv('CPI_data.csv', header = TRUE, stringsAsFactors = TRUE)
#take out 2016 as only a few months of data
price_data2 = price_data[price_data$Delv_Year != 2016, ]
price_year = data.frame(year = 2006:2015, 
  mean_price = tapply(price_data2$Breadmaking_Wheat, INDEX = price_data2$Delv_Year, FUN = mean),
  price_min = tapply(price_data2$Breadmaking_Wheat, INDEX = price_data2$Delv_Year, FUN = min),
  price_max = tapply(price_data2$Breadmaking_Wheat, INDEX = price_data2$Delv_Year, FUN = max),
  CPI = CPI_data$CPI)
#change in index realtive to base year 2014 which is year for which we have nix data
price_year$CPI_ratio = price_year$CPI / price_year$CPI[price_year$year == 2014]
price_year$mean_price_adj = price_year$mean_price / price_year$CPI_ratio
price_year$min_price_adj = price_year$price_min / price_year$CPI_ratio
price_year$max_price_adj = price_year$price_max / price_year$CPI_ratio

#Nix 2014 yeild data averages gives 
# milling wheat 
# low: 6.15 tonnes.ha
# average: 7.7
# high: 9.0

mean_est = mean(price_year$mean_price_adj) * 7.7
lower = min(price_year$min_price_adj) * 6.15
upper = max(price_year$max_price_adj) * 9

mean_est
lower
upper
