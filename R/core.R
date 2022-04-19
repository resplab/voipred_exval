
voi_ex_glm <- function(model, val_data, method=c("bootstrap","model_based_ll","model_based_bs","asymptotic"), Bayesian_bootstrap=F, n_sim=1000, zs=(0:99)/100)
{
  n <- dim(val_data)[1]
  
  if(method=="asymptotic")
  {
    ENB_perfect <- ENB_current <- rep(0, length(zs))
    for(j in 1:length(zs))
    {
      NB_model <- sum((val_data$pi>zs[j])*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
      NB_all <- sum((val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
      
      parms <- NB_BVN(val_data$Y,val_data$pi,zs[j])
      
      ENB_perfect[j] <- do.call(mu_max_truncated_bvn,as.list(parms))
      ENB_current[j] <- max(0,NB_model,NB_all)
    }
    return(data.frame(z=zs, ENB_perfect=ENB_perfect, ENB_current=ENB_current, EVPIv=ENB_perfect-ENB_current))
  }
  
  NB_model <- NB_all <- matrix(0, n_sim, ncol=length(zs))
  
  if(method=="bootstrap")
  {
    for(i in 1:n_sim)
    {
      w_x <- voipred:::bootstrap(n, Bayesian_bootstrap)
      for(j in 1:length(zs))
      {
        NB_model[i,j] <- sum(w_x*(val_data$pi>zs[j])*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
        NB_all[i,j] <- sum(w_x*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
      }
      #if(i%%100 ==0) {plot(zs,NB_model[i,]); lines(zs,NB_all[i,])}
    }
  }
  else if(method=="model_based_ll")
  {
    local_model <-suppressWarnings(glm(formula=model$call$formula, data=val_data, family=binomial(link="logit")))
    mus <- as.vector(coefficients(local_model))
    covmat <- vcov(local_model)
    for(i in 1:n_sim)
    {
      w_x <- voipred:::bootstrap(n, Bayesian_bootstrap)
      betas <- rmvnorm(1,mus,covmat)
      p <- as.vector(1/(1+exp(-(model.matrix(model, data=val_data)%*%t(betas)))))
      for(j in 1:length(zs))
      {
        NB_model[i,j] <- sum(w_x*(val_data$pi>zs[j])*(p-(1-p)*zs[j]/(1-zs[j])))/n
        NB_all[i,j] <- sum(w_x*(p-(1-p)*zs[j]/(1-zs[j])))/n
      }
      #if(i%%100 ==0) {plot(zs,NB_model[i,]); lines(zs,NB_all[i,])}
    }
  }
  else if(method=="model_based_bs")
  {
    for(i in 1:n_sim)
    {
      w_x <- voipred:::bootstrap(n, Bayesian_bootstrap)
      val_data$w_y <- voipred:::bootstrap(n, Bayesian_bootstrap)
      local_model <-suppressWarnings(glm(formula=model$call$formula, data=val_data, family=binomial(link="logit"), weights = w_y))
      p <- predict(local_model, type="response", newdata = val_data)
      #betas <- as.vector(coefficients(local_model))
      #p <- as.vector(1/(1+exp(-(model.matrix(model, data=val_data)%*%(betas)))))
      for(j in 1:length(zs))
      {
        NB_model[i,j] <- sum(w_x*(val_data$pi>zs[j])*(p-(1-p)*zs[j]/(1-zs[j])))/n
        NB_all[i,j] <- sum(w_x*(p-(1-p)*zs[j]/(1-zs[j])))/n
      }
      #if(i%%100 ==0) {plot(zs,NB_model[i,]); lines(zs,NB_all[i,])}
    }
  }
  else
  {
    stop("Method ",method," is not recognized.")
  }
  
  ENB_perfect <- ENB_current <- EVPIv <- rep(NA,length(zs))
  for(i in 1:length(zs))
  {
    ENB_perfect[i] <- mean(pmax(NB_model[,i],NB_all[,i],0))
    ENB_current[i] <- max(mean(NB_model[,i]),mean(NB_all[,i]),0)
    EVPIv[i] <- ENB_perfect[i] - ENB_current[i]
  }
  
  data.frame(z=zs, ENB_perfect=ENB_perfect, ENB_current=ENB_current, EVPIv=EVPIv)
}







NB_BVN <- function(y,pi,z){
  # set up
  n <- length(y)
  a <- as.numeric(pi>z)
  rho <- mean(y)
  TPR <- mean(y*a)/rho
  FPR <- mean(a*(1-y))/(1-rho)
  tz <- z/(1-z)
  
  # mean
  mu <- c(rho*TPR-tz*(1-rho)*FPR,rho-tz*(1-rho))
  
  # var
  sig_Y <- rho*(1-rho)/n
  sig_TPR <- rho*TPR*(1-rho*TPR)/n
  sig_FPR <- (1-rho)*FPR*(1-(1-rho)*FPR)/n
  
  sig11 <-  sig_TPR + tz^2 * sig_FPR + 2 * tz * rho * (1-rho) * TPR * FPR / n
  sig22 <- (1+tz)^2 * sig_Y
  sig12 <- sig_Y * (1+tz)* (TPR + tz * FPR)
  cor_coefficient <- sig12/(sqrt(sig11)*sqrt(sig22))
  
  # mean of NB_model, mean of NB_all,
  # var of NB_model, var of NB_all,
  # correlation of NB_model and NB_all
  return(c(mu,sig11,sig22,cor_coefficient))
}


mu_max_truncated_bvn <- function(mu1,mu2,sig1,sig2,rho,exact=T){
  f1 <-  function(mu1,mu2,sig1,sig2,rho){
    tmp1 <- sig1-rho*sig2
    tmp2 <- (-sig1*mu2+rho*sig2*mu1)/(sig1*sig2*sqrt(1-rho^2))
    mu1 * ( as.numeric(tmp1>0) +  as.numeric(tmp1==0)*pnorm(tmp2) ) - pnorm(tmp2)*(-sig1*dnorm(-mu1/sig1) + mu1*pnorm(-mu1/sig1))
  }
  
  f2 <- function(mu1,mu2,sig1,sig2,rho,exact=T){
    if(exact){
      
      tmp1 <- (sig1-rho*sig2)
      alpha <- (sig1*mu2-rho*sig2*mu1)/tmp1
      beta <- sig1*sig2*sqrt(1-rho^2)/tmp1
      if(beta>0){
        a <- (mu1-alpha)/beta
        b <- sig1/beta
        T1_rho <- -1/sqrt(1+b^2)
        T1 <- mu1 * (pnorm((-a/b)/sqrt(1+(1/b)^2) )  - as.numeric(pmvnorm(lower = c(-Inf,-Inf),upper=c(-a/sqrt(1+b^2),-alpha/beta),
                                                                          mean=c(0,0),sigma=matrix(c(1,T1_rho,T1_rho,1),2))[[1]]))
        a <- (alpha-mu1)/sig1 
        b <- beta/sig1
        T2_t <- sqrt(1+b^2)
        T2 <- -sig1 / T2_t * dnorm(a/T2_t) * (1-pnorm(-T2_t*alpha/beta+a*b/T2_t))
        return(T1+T2)
      } else{
        beta <- abs(beta)
        a <- (mu1-alpha)/beta
        b <- sig1/beta
        T1_rho <- -1/sqrt(1+b^2)
        T1 <- -mu1 * (pnorm((-a/b)/sqrt(1+(1/b)^2) )  - as.numeric(pmvnorm(lower = c(-Inf,-Inf),upper=c(-a/sqrt(1+b^2),-alpha/beta),
                                                                           mean=c(0,0),sigma=matrix(c(1,T1_rho,T1_rho,1),2))[[1]]))
        a <- (alpha-mu1)/sig1 
        b <- beta/sig1
        T2_t <- sqrt(1+b^2)
        T2 <- sig1 / T2_t * dnorm(a/T2_t) * (1-pnorm(-T2_t*alpha/beta+a*b/T2_t))
        return(T1+T2)
        
      }
    } else{
      
      tmp1 <- (sig1-rho*sig2)/(sig1*sig2*sqrt(1-rho^2))
      tmp2 <- (-sig1*mu2+rho*sig2*mu1)/(sig1*sig2*sqrt(1-rho^2))
      
      int_LHS <- function(y){
        (-sig1*dnorm((y-mu1)/sig1)+mu1*pnorm((y-mu1)/sig1))*tmp1*dnorm(y*tmp1+tmp2)
      }
      return(integrate(int_LHS,0,Inf,subdivisions=5000L)[[1]])
    }
  }
  
  f1(mu1,mu2,sig1,sig2,rho) + f1(mu2,mu1,sig2,sig1,rho) - f2(mu1,mu2,sig1,sig2,rho,exact) - f2(mu2,mu1,sig2,sig1,rho,exact)
}


