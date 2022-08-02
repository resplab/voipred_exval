library(MASS); data(birthwt)
library(predtools)

n <- dim(birthwt)[1]
z <- 0.2 #This is the risk threshold

#Step 1: calculate predicted probabilities. #Our model is logit(pi) = 2 -0.05*age -0.01*lwt
pi <- 1/(1+exp(-(2-0.05*birthwt$age-0.01*birthwt$lwt))) 

#Step 2: Estimate the moments
Y <- birthwt$low; a <- I(pi>z); p0 <- mean(Y); TPR <- mean(Y*a)/p0; FPR <- mean(a*(1-Y)/(1-p0)); tz <- z/(1-z)

NB_model <- mean(a*(Y-(1-Y)*tz)); NB_all <- mean(Y-(1-Y)*tz)

var_Y <- p0*(1-p0)/n
var_TPR <- p0*TPR*(1-p0*TPR)/n
var_FPR <- (1-p0)*FPR*(1-(1-p0)*FPR)/n

sd_model <- sqrt(var_TPR+tz^2*var_FPR+2*tz*p0*(1-p0)*TPR*FPR/n)
sd_all <- sqrt((1+tz)^2 * var_Y)
rho <- var_Y * (1+tz)* (TPR + tz * FPR)/sd_model/sd_all

#Step 3: Calculate E{max⁡{0,NB_model,NB_all } }
A <- mu_max_trunc_bvn(NB_model,NB_all,sd_model,sd_all,rho)

# Step 4: EVPI calculation
EVPI <- A - max(0,NB_model,NB_all)
