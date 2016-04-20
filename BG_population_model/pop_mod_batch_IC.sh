#!/bin/bash

#request 7.5 hrs of time estimated time to run is 5hr  
#$ -l h_cpu=07:30:00

#$ -q short.q

#request 4G per core for 6 cores = 24G (48 virtual)
#$ -l mem=30G -l rmem=15G

#request 6 cores in an openMP environment
#$ -pe openmp 6

#email me the results
#$ -M shaun.coutts@gmail.com
#$ -m be

#set the directory to the one with the R script
cd /data/bo1src/BG_HR_pop_model

#load the R modual
module load apps/R

#call the batch script to run the script
R CMD BATCH --vanilla pop_model_filtering_sens_IB.R pop_mod_R_std_output2.Rout
