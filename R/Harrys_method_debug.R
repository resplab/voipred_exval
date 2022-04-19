z <- 0.02
n_sim <- 10^4

NB_hat_model <- mean((val_data$pi>z)*(val_data$Y-(1-val_data$Y)*z/(1-z))) 
NB_hat_all <- mean(1*(val_data$Y-(1-val_data$Y)*z/(1-z)))


NB_model <- NB_all <- X <- rep(NA,n_sim)
for(i in 1:n_sim) 
{
  bs_data <- val_data[sample(1:dim(val_data)[1],size = dim(val_data)[1], T),]
  NB_model[i] <- mean((bs_data$pi>z)*(bs_data$Y-(1-bs_data$Y)*z/(1-z))) 
  NB_all[i] <- mean(1*(bs_data$Y-(1-bs_data$Y)*z/(1-z)))
  
  X[i] <- max(0,NB_model[i],NB_all[i]) 
}

plot(NB_model,NB_all)
lines(c(0,1),c(0,1))

c(mean(NB_model),mean(NB_all),var(NB_model),var(NB_all),cor(NB_model,NB_all))
cat(res <- NB_BVN(val_data$Y,val_data$pi,z))

first_term <- do.call(mu_max_truncated_bvn,as.list(res))
mean(X)

EVPI <- c()
for(z in (10:99)/1000)
{
  cat(z,'\n')
  res <- NB_BVN(val_data$Y,val_data$pi,z)
  first_term <- do.call(mu_max_truncated_bvn,as.list(res))
  second_term <- max(0,res[1],res[2])
  EVPI <- c(EVPI,first_term-second_term)
}