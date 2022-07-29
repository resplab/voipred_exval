dca <- function(val_data, zs=(0:99)/100, weights=NULL)
{
  n <- dim(val_data)[1]
  
  NB_model <- NB_all <- rep(0, length(zs))
  
  if(is.null(weights)) weights <- rep(1,n)
    
  for(j in 1:length(zs))
  {
    NB_model[j] <- sum(weights*(val_data$pi>zs[j])*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
    NB_all[j] <- sum(weights*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
  }
  
  return(data.frame(z=zs, NB_model=NB_model, NB_all=NB_all))
}



voi_ex_glm <- function(model, val_data, method=c("bootstrap","model_based_ll","model_based_bs","asymptotic"), Bayesian_bootstrap=F, n_sim=1000, zs=(0:99)/100, weights=NULL)
{
  n <- dim(val_data)[1]
  
  if(method=="asymptotic")
  {
    if(is.null(weights)) weights <- rep(1,n)
        
    ENB_perfect <- ENB_current <- rep(0, length(zs))
    for(j in 1:length(zs))
    {
      NB_model <- sum(weights*(val_data$pi>zs[j])*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
      NB_all <- sum(weights*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
      parms <- NB_BVN(val_data$Y,val_data$pi,zs[j], weights)
      # print(paste0(n,"_",zs[j],"_",parms[5]))
      if(is.na(parms[5])){
        ENB_perfect[j] <- ENB_current[j] <- max(0,NB_model,NB_all)
        } else{
      if(parms[5]>0.999999) parms[5]<- 0.999999
      if(parms[5]< -0.999999) parms[5]<- -0.999999
  
    
      tryCatch(
        {ENB_perfect[j] <- do.call(mu_max_truncated_bvn,as.list(parms))}
      , error=function(cond) {
        return(NULL)
      })
      ENB_current[j] <- max(0,NB_model,NB_all)
        }
    }
    return(data.frame(z=zs, ENB_perfect=ENB_perfect, ENB_current=ENB_current, EVPIv=ENB_perfect-ENB_current))
  }
  
  NB_model <- NB_all <- matrix(0, n_sim, ncol=length(zs))
  
  if(method=="bootstrap")
  {
    for(i in 1:n_sim)
    {
      w_x <- bootstrap(n, Bayesian_bootstrap, weights = weights)
      for(j in 1:length(zs))
      {
        NB_model[i,j] <- sum(w_x*(val_data$pi>zs[j])*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
        NB_all[i,j] <- sum(w_x*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
      }
      #if(i%%100 ==0) {plot(zs,NB_model[i,]); lines(zs,NB_all[i,])}
    }
  }
  else
  {
    stop("Method ",method," is not recognized.")
  }
  
  ENB_model <- ENB_all <- ENB_perfect <- ENB_current <- EVPIv <- p_useful <- rep(NA,length(zs))
  for(i in 1:length(zs))
  {
    ENB_model[i] <- mean(NB_model[,i])
    ENB_all[i] <- mean(NB_all[,i])
                                              
    ENB_perfect[i] <- mean(pmax(NB_model[,i],NB_all[,i],0))
    ENB_current[i] <- max(ENB_model[i],ENB_all[i],0)
    EVPIv[i] <- ENB_perfect[i] - ENB_current[i]
    p_useful[i] <- mean((pmax(NB_model[,i],NB_all[,i],0)-NB_model[,i])==0)
  }
  
  data.frame(z=zs, ENB_model=ENB_model, ENB_all=ENB_all, ENB_current=ENB_current, ENB_perfect=ENB_perfect, EVPIv=EVPIv, p_useful=p_useful)
}







NB_BVN <- function(y,pi,z,weights=NULL){
  # set up
  n <- length(y)
  a <- as.numeric(pi>z)
  if(is.null(weights))
  {
    rho <- mean(y)
    TPR <- mean(y*a)/rho
    FPR <- mean(a*(1-y))/(1-rho)
  }
  else
  {
    rho <- weighted.mean(y,w = weights)
    TPR <- weighted.mean(y*a, w = weights)/rho
    FPR <- weighted.mean(a*(1-y), w=weights)/(1-rho) 
  }
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


mu_max_truncated_bvn <-
  function(mu1, mu2, sig1sq, sig2sq, rho) {
    sig1 <- sqrt(sig1sq)
    sig2 <- sqrt(sig2sq)
    f1 <-  function(mu1, mu2, sig1, sig2, rho) {
      tmp1 <- sig1 - rho * sig2
      tmp2 <-
        (-sig1 * mu2 + rho * sig2 * mu1) / (sig1 * sig2 * sqrt(1 - rho ^
                                                                 2))
      mu1 * (as.numeric(tmp1 > 0) +  0 * as.numeric(tmp1 == 0) * pnorm(tmp2)) -
        pnorm(tmp2) * (-sig1 * dnorm(-mu1 / sig1) + mu1 * pnorm(-mu1 / sig1))
    }
    
    f2 <- function(mu1, mu2, sig1, sig2, rho) {
      tmp1 <- (sig1 - rho * sig2)
      alpha_num <- (sig1 * mu2 - rho * sig2 * mu1)
      beta_num <- sig1 * sig2 * sqrt(1 - rho ^ 2)
      alpha <- alpha_num / tmp1
      beta <- beta_num / tmp1
      
      if (tmp1 > 0) {
        a <- (tmp1 * mu1 - alpha_num) / beta_num
        b <- tmp1 * sig1 / beta_num
        T1_rho <- -1 / sqrt(1 + b ^ 2)
        
        if (tmp1 < 1e-10) {
          T1_1 <- pnorm(alpha_num / beta_num)
        } else{
          T1_1 <-  pnorm((-a / b) / sqrt(1 + (1 / b) ^ 2))
        }
        
        T1 <-
          mu1 * (T1_1  - as.numeric(pmvnorm(
            lower = c(-Inf,-Inf),
            upper = c(-a / sqrt(1 + b ^ 2),-alpha_num / beta_num),
            mean =
              c(0, 0),
            sigma = matrix(c(1, T1_rho, T1_rho, 1), 2)
          )[[1]]))
        a <- (alpha - mu1) / sig1
        b <- beta / sig1
        T2_t <- sqrt(1 + b ^ 2)
        
        if (tmp1 < 1e-10) {
          T2 <-
            -sig1 / b * dnorm(alpha_num / beta_num) * (1 - pnorm(-mu1 / sig1))
        } else{
          T2 <-
            -sig1 / T2_t * dnorm(a / T2_t) *
            (1 - pnorm(-T2_t * alpha_num / beta_num + a * b / T2_t))
        }
        
        return(T1 + T2)
      }
      else{
        beta <- abs(beta)
        
        a <- (abs(tmp1) * mu1 + alpha_num) / beta_num
        b <- abs(tmp1) * sig1 / beta_num
        T1_rho <- -1 / sqrt(1 + b ^ 2)
        
        if (abs(tmp1) < 1e-10) {
          T1_1 <- pnorm(-alpha_num / beta_num)
        } else{
          T1_1 <-  pnorm((-a / b) / sqrt(1 + (1 / b) ^ 2))
        }
        
        T1 <-
          -mu1 * (T1_1  - as.numeric(pmvnorm(
            lower = c(-Inf,-Inf),
            upper = c(-a / sqrt(1 + b ^ 2), alpha_num / beta_num),
            mean =
              c(0, 0),
            sigma = matrix(c(1, T1_rho, T1_rho, 1), 2)
          )[[1]]))
        
        a <- (alpha - mu1) / sig1
        b <- beta / sig1
        T2_t <- sqrt(1 + b ^ 2)
        
        if (abs(tmp1) < 1e-10) {
          T2 <-
            sig1 / b * dnorm(-alpha_num / beta_num) * (1 - pnorm(-mu1 / sig1))
        } else{
          T2 <-
            sig1 / T2_t * dnorm(a / T2_t) * (1 - pnorm(T2_t * alpha_num / beta_num +
                                                         a * b / T2_t))
        }
        
        return(T1 + T2)
        
      }
    }
    
    f1(mu1, mu2, sig1, sig2, rho) + f1(mu2, mu1, sig2, sig1, rho) -
      f2(mu1, mu2, sig1, sig2, rho) - f2(mu2, mu1, sig2, sig1, rho)
  }


bootstrap <- function (n, Bayesian=F, weights=NULL)
{
  if(Bayesian)
  {
    u <- c(0,sort(runif(n-1)),1)
    u <- (u[-1] - u[-length(u)])*n
    if(!is.null(weights)) u <- u*weights*n/sum(u*weights)
  }
  else
  {
    if(is.null(weights)) weights <-rep(1/n,n)
    u <- rmultinom(1,n,weights)
  }
  as.vector(u)
}



