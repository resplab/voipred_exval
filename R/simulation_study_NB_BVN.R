source("NB_BVN.R")
source("mu_trunc_BVN.R")

# sim study ---------------------------------------------------------------
stylized_sim <- function(betas=c(-2,0,1), N=100, sample_size=5000, n_sim=1000,z=0.2)
{
 set.seed(1234)
  res <-c()
  for(i in 1:N){
      res[[i]] <- stylized_sim_internal(betas, sample_size, n_sim, z)
  }
  return(res)
}


stylized_sim_internal <- function(betas, sample_size, n_sim,z)
{
  formula=paste("Y~",paste0("X",1:(length(betas)-1),collapse = "+"))
  
  sample_data <<-generate_data(sample_size, betas)
  
  reg <- glm(formula=formula, data=sample_data, family=binomial(link="logit"))
  
  pi <- predict(reg, type="response") #Predicted risks
  
  exact <- NB_BVN(sample_data$Y,pi,z)
  
  bootstrap <- evpp.glm(reg, formula,n_sim,sample_data,pi,z)
  
  out <- data.frame(cbind(exact,bootstrap))
  rownames(out) <- c("mu_model","mu_all","sig_model","sig_all","cor")
  
  return(out)
}



generate_data <- function(n, betas)
{
  n_var <- length(betas)-1 #one is intercept
  X <- NULL
  
  for(i in 1:n_var)
  {
    X <- cbind(X, rnorm(n))
  }
  
  colnames(X)<-paste0("X",1:n_var)
  
  lin <- betas %*% t(cbind(1,X))
  
  Y <- rbinom(n, 1, 1/(1+exp(-lin)))
  
  # message(paste0("P(Y=1)=",mean(Y)))
  
  return(as.data.frame(cbind(X, Y=Y)))
}

evpp.glm <- function(reg, formula,n_sim,data_set,pi,z){
  n <- nrow(data_set)
  
  #We hold simulation results in memory for exploration (not necessary)
  # out <- data.frame("NB_model"=double(),"NB_all"=double(),"NB_max"=double())
  out <- data.frame("NB_model"=double(),"NB_all"=double())
  
  for(i in 1:n_sim){
  bs_data_set <- data_set[sample(1:n, n, replace = T),] #Create a bootstrap sample and fit the model again
  bs_model <- glm(formula=formula, data=bs_data_set, family=binomial(link="logit"))
  
  #Predict risks from this model are applied to the original data
  p <- predict(bs_model, newdata = data_set, type="response")
  
  #Bayesiuan NB calculations. p are taken as random draws from the distribution of correct risks
  NB_all <- mean(p-(1-p)*z/(1-z)) #NB of treating all
  NB_model <- mean((pi>z)*(p-(1-p)*z/(1-z))) #NB of using the model
  # NB_max <- mean((p>z)*(p-(1-p)*z/(1-z))) #NB of using the correct risks
  
  # out[i,] <- c(NB_model,NB_all,NB_max)
  out[i,] <- c(NB_model,NB_all)
  
  }
  var <- apply(out,2,var)
  
  c(colMeans(out),var,cov(out$NB_model,out$NB_all)/prod(sqrt(var)))
  
}


out <- stylized_sim(N=5)


# exemplary dataset -------------------------------------------------------
library(MASS)
#birthw is an exemlary dataset containing the low birthweight status of newborns and some covairates
data_set <- birthwt
n <- dim(birthwt)[1]

z <- 0.2 #This is the risk threshold

model <- glm(low ~ age + lwt, family=binomial(link="logit"), data=data_set) #Our risk prediciton model

pi <- predict(model, type="response") #Predicted risks

look <- NB_BVN(data_set$low,pi,z)

#We hold simulation results in memory for exploration (not necessary)
out <- data.frame("NB_model"=double(),"NB_all"=double(),"NB_max"=double())
results <- c()
#Looping over 1000 simulations
for(i in 1:1000){
  bs_data_set <- data_set[sample(1:n, n, replace = T),] #Create a bootstrap sample and fit the model again
  bs_model <- glm(low ~ age + lwt, family=binomial(link="logit"), data=bs_data_set)
  #Predict risks from this model are applied to the original data
  p <- predict(bs_model, newdata = data_set, type="response")
  
  #Bayesiuan NB calculations. p are taken as random draws from the distribution of correct risks
  NB_all <- mean(p-(1-p)*z/(1-z)) #NB of treating all
  NB_model <- mean((pi>z)*(p-(1-p)*z/(1-z))) #NB of using the model
  NB_max <- mean((p>z)*(p-(1-p)*z/(1-z))) #NB of using the correct risks
  
  out[i,] <- c(NB_model,NB_all,NB_max)
}
c(colMeans(out[,-3]),apply(out[,-3],2,var),cov(out$NB_model,out$NB_all))
look
