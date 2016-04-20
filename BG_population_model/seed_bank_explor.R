#seed bank overlap
sb <- read.csv('/home/shauncoutts/Dropbox/projects/seed_bank/SEEDBANK.CSV', header = T)
ob_name <- load("/home/shauncoutts/Dropbox/projects/COMPADRE_spatial_replication/Data/Demography/Rlease_version/COMPADRE Jul 21 2014.RData")

comp_spp <- unique(gsub('_', x = compadre$metadata$SpeciesAccepted, ' '))
sb_spp <- unique(as.character(sb$Species))

overlap <- comp_spp[comp_spp %in% sb_spp ] # 92 spp overlap

#dispersal data
disp <- read.csv('/home/shauncoutts/Dropbox/projects/seed_bank/DispersalDistanceData.csv', header = TRUE) 
d_spp <- as.character(disp$Species)
overlap_disp <- d_spp[d_spp %in%  overlap] #43 spp overlap


#Black grass seed bank
sb[grepl('Alopecurus myosuroides', sb$Species), ]
#Information extracted from seed bank data base to look at black grass seed banks 

#interpriting seed longevity codes
# 1: seeds survive <1 year
# 2: seeds surviv >1 and <5 years
# 3: seeds survive >5 years
# 4: not able to determine

#interpriting methodology code
# 1: seeds buried in garden plot, no additional disturbance
# 2: seeds buried in garden plot, additional disturbance
# 3: deliberate burial in field
# 4: field sample taken and seeds extracted and tested for germination or viability
# 5: soil sample taken from field, seeds germinated without extraction
# 6: as 5 but seeds germinated in the field
# 7: sequential sampling of seed bank

black_grass_sb = sb[grepl('Alopecurus myosuroides', sb$Species), ]
#look only at the seed buiral experiments 
black_grass_sb_bury = black_grass_sb[black_grass_sb$Method <= 3, ]
#This data suggests that 12 burial experiments put seed survival betwen 2 and 5 years, with one outlier at >11 years
#translating this to a constant per year survival find seed_sur such that M seeds surive y years 
M = 0.03
y = 5

seed_sur = exp(log(M) / y)
#go back the other way to chech we get M
seed_sur ^ y

#now look for the mean of all burial trials, will need the estimated max years, and class since one often provides a lower while the other the upper limit
bg_sb = black_grass_sb_bury[c('Seed.bank.type', 'Max.longevity..y.')] 
# don't know the size of class 3 which means that we cannot get an average of the burial trials, but we can get a median of the mid-point estimates per study
bg_sb$point_est = sapply(strsplit(as.character(bg_sb$Max.longevity..y. ), split = '>'), FUN = function(x){
  return(ifelse(length(x) == 1, TRUE, FALSE))
})
bg_sb$min_long = sapply(strsplit(as.character(bg_sb$Max.longevity..y. ), split = '>'), FUN = function(x){
  return(as.numeric(x[length(x)]))
})
bg_sb$est_long = sapply(seq_along(bg_sb$min_long), FUN = function(x){
  if(bg_sb$point_est[x]){ 
    return(bg_sb$min_long[x])
  }else{
    if(bg_sb$Seed.bank.type[x] == 1){
      return((1 - bg_sb$min_long[x]) / 2)
    }else{
      if(bg_sb$Seed.bank.type[x] == 2){
	return(((5 - bg_sb$min_long[x]) / 2) + bg_sb$min_long[x])
      }else{
	return(bg_sb$min_long[x] + 1) #This just says the max longevity is observed max + 1 as we can't know the true max, why we can only use a median not a mean 
      }
    }
  }
})

med_long = median(bg_sb$est_long)
med_seed_sur = exp(log(0.05) / med_long)# turns out to be 0.4498 


m = s^y















