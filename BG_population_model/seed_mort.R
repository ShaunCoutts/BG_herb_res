# function for seed mortality per age from Colbach et al 2006
seed_mort = function(age, seed0){
  ((age / (age - 1))^-0.02444) * seed0
}

#age seems to be in days and survival gets better as age increases, but if we intergrate across a whole year
x = 2:365
S = numeric()
S[1] = 100
for(a in x) S[a] = seed_mort(a, S[a - 1])

#This function implies that seed survival is up around 86%, which seems much too high, but also exclues germination
#so maybe this is the actual mortality portion
