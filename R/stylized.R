library(MASS)
remotes::install_github("resplab/predtools")
library("predtools")
data(birthwt)

n <- dim(birthwt)[1]
z <- 0.2 #This is the risk threshold

#Our model is logit(pi) = 2 -0.05*age -0.01*lwt

# set seed for reproducibility
set.seed(2022)

#Step 1: calculate predicted probabilities
birthwt$pi <- 1/(1+exp(-(2-0.05*birthwt$age-0.01*birthwt$lwt)))  #Predicted risks

# Bootstrap method

#Step 2: The bootstrap Monte Carlo simulation
N <- 10000
NBmodel <- NBall <- NBmax <- rep(0,N)
for(i in 1:N) 
{
  bsdata <-  birthwt[sample(1:n, n, replace = T),] 
  NBall[i] <- mean(bsdata$low-(1-bsdata$low)*z/(1-z)) 
  NBmodel[i] <- mean((bsdata$pi>z)*(bsdata$low-(1-bsdata$low)*z/(1-z))) #NB of using the model
}

#Step 3: EVPI calculation
EVPI <- mean(pmax(0,NBmodel,NBall))-max(0,mean(NBmodel),mean(NBall))
# 0.0008698413

# Asymptotic method

#Step 2: Estimate the moments
Y <- birthwt$low; pi <- birthwt$pi; a <- I(pi>z); rho <- mean(Y); TPR <- mean(Y*a)/rho; FPR <- mean(a*(1-Y)/(1-rho)); tz <- z/(1-z)
NBmodel <- mean(a*(Y-(1-Y)*tz))
NBall <- mean(Y-(1-Y)*tz)

Sigma <- matrix(0,nrow=2,ncol=2)

sig_Y <- rho*(1-rho)/n
sig_TPR <- rho*TPR*(1-rho*TPR)/n
sig_FPR <- (1-rho)*FPR*(1-(1-rho)*FPR)/n

Sigma[1,1] <-  sig_TPR + tz^2 * sig_FPR + 2 * tz * rho * (1-rho) * TPR * FPR / n
Sigma[2,2] <- (1+tz)^2 * sig_Y
Sigma[1,2] <- Sigma[2,1] <- sig_Y * (1+tz)* (TPR + tz * FPR)

#Step 3: Calculate E{max⁡{0,NB_model,NB_all } }
A <- mu_max_trunc_bvn(c(NBmodel,NBall),Sigma)

# Step 4: EVPI calculation
EVPI <- A - max(0,NBmodel,NBall)
# 0.0007555751