library(MASS); data(birthwt)

n <- dim(birthwt)[1]
z <- 0.2 #This is the risk threshold

#Step 1: calculate predicted probabilities. #Our model is logit(pi) = 2 -0.05*age -0.01*lwt
pi <- 1/(1+exp(-(2-0.05*birthwt$age-0.01*birthwt$lwt)))  #Predicted risks

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
