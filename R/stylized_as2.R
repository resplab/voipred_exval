library(MASS); data(birthwt)
library(predtools)

n <- dim(birthwt)[1]
z <- 0.2 #This is the risk threshold

#Step 1: calculate predicted probabilities. #Our model is logit(pi) = 2 -0.05*age -0.01*lwt
pi <- 1/(1+exp(-(2-0.05*birthwt$age-0.01*birthwt$lwt))) 

#Step 2: Estimate the moments* 
Y <- birthwt$low; 

moments <- predtools::calc_NB_moments(Y,pi,z)

#Step 3: Closed-form solution for E{max{0,NB_model,NB_all } }
A <- do.call(mu_max_trunc_bvn, as.list(moments))

# Step 4: EVPI calculation
EVPI <- A - max(0,moments[1],moments[2])
