library(truncnorm)
library(tidyverse)
source("mu_trunc_BVN.R")

# normal approximation method ---------------------------------------------
# first and second moments of the max of BVN
mu_var_max_bvn <- function(mu1,mu2,sig1,sig2,rho){
  diff_mu <- mu1-mu2
  theta <- sqrt(sig1^2+sig2^2-2*rho*sig1*sig2)
  alpha <- diff_mu/theta
  mu <- mu1*pnorm(alpha) + mu2*pnorm(-alpha) + theta*dnorm(alpha)
  sigma_sq <- (sig1^2+mu1^2)*pnorm(alpha)+(sig2^2+mu2^2)*pnorm(-alpha)+(mu1+mu2)*theta*dnorm(alpha)-
    mu^2
  return(c(mu,sigma_sq))
}

# normal apprximation
# Y_normal_approx <- rtruncnorm(n,a=0,mean=normal_approx_moments[1],sd = sqrt(normal_approx_moments[2]))
mu_trunc_normal <- function(mu,sigma,a,b=Inf){
  alpha <- (a-mu)/sigma
  beta <- (b-mu)/sigma
  Z <- pnorm(beta)-pnorm(alpha)
  mu+(dnorm(alpha)-dnorm(beta))/Z*sigma
}



simulation_study <- function(num_sim,n=10000,std_case=F,exact=T,range_sigma=c(0.1,3)){
  results <- c()
  
  for(i in 1:num_sim){
    if(!std_case){
    mu1 <- rnorm(1)

    mu2 <- rnorm(1)

    sig1 <- runif(1,range_sigma[1],range_sigma[2])

    sig2 <- runif(1,range_sigma[1],range_sigma[2])

    rho <- runif(n = 1,min = -0.8, max = 0.8)

    sig <- matrix(c(sig1^2, sig1*sig2*rho, sig1*sig2*rho, sig2^2),
                  2)
    } else{
      mu1 <- 0
      
      mu2 <- 0
      
      sig1 <- 1
      
      sig2 <- 1
      
      rho <- runif(n = 1,min = -0.8, max = 0.8)
      
      sig <- matrix(c(sig1^2, sig1*sig2*rho, sig1*sig2*rho, sig2^2),
                    2)
      
    }
    
    
    dat <- rmvnorm(n,mean=c(mu1,mu2),sigma = sig)
    est <- mean(pmax(dat[,1],dat[,2],0))
    
    normal_approx_moments <- mu_var_max_bvn(mu1,mu2,sig1,sig2,rho)
    scalar <- (1-pmvnorm(upper=c(0,0),
                         mean=c(mu1,mu2),
                         sigma=sig))
    approx <- mu_trunc_normal(normal_approx_moments[1],sqrt(normal_approx_moments[2]),0)*scalar
    
    mine <- mu_max_truncated_bvn(mu1,mu2,sig1,sig2,rho,exact)
    err_approx <- est-approx
    err_harry <- est-mine
    results[[i]] <- c(est,approx,mine,err_approx,err_harry, mu1,mu2,sig1,sig2,rho,scalar)
  }
  
  results <- do.call(rbind,results) %>% 
    as.data.frame()
  
  colnames(results) <- c("sampling","approx","harry","diff_approx","diff_harry","mu1","mu2","sig1","sig2","rho","scalar")
  
  return(results)
}

n_sim <- 100

# standard BVN case
std_case <- simulation_study(n_sim,std_case = T)
# in comparsion to the MC estimator
boxplot(look1$diff_approx) # normal approximation
boxplot(look1$diff_harry) # exact method

# non-standard BVN case 
nonstd_case <- simulation_study(n_sim,std_case = F)
boxplot(look2$diff_approx) # normal approximation
boxplot(look2$diff_harry) #exact method

# sanity check ------------------------------------------------------------
# check integral of dv
check_v <- function(n_sim){
  tot <- 0 
  for(i in 1:n_sim){
    mu1 <- rnorm(1)
    sigma1 <- runif(1,0.5,3)
    
    LHS_int <- function(y){
      y*1/sigma1*dnorm((y-mu1)/sigma1)
    }
    
    RHS <- function(y){
      -sigma1*dnorm((y-mu1)/sigma1) + mu1*pnorm((y-mu1)/sigma1)
    }
    
    upper <- rnorm(1)
    
    tot <- tot+abs(integrate(LHS_int,lower=-Inf,upper,subdivisions = 5000L)[[1]]-RHS(upper))
    
  }
  
  tot/n_sim
}

check_v(1000)

# evaluation of uv
check_uv <- function(n_sim){
  tot <- 0
  for(i in 1:n_sim){
    mu1 <- rnorm(1)
    
    mu2 <- rnorm(1)
    
    sig1 <- runif(1,0.5,2)
    
    sig2 <- runif(1,0.5,2)
    
    rho <- runif(n = 1,min = -1, max = 1)
    
    tmp1 <- (sig1-rho*sig2)/(sig1*sig2*sqrt(1-rho^2))
    tmp2 <- (-sig1*mu2+rho*sig2*mu1)/(sig1*sig2*sqrt(1-rho^2))
    
    LHS <- function(y){
      pnorm(y*tmp1 + tmp2)*(-sig1*dnorm((y-mu1)/sig1) + mu1*pnorm((y-mu1)/sig1))
    }
    
    RHS <- function(y){
      mu1*(as.numeric((sig1-rho*sig2)>0) + as.numeric((sig1-rho*sig2)==0)*pnorm(tmp2)) -
        pnorm(tmp2)*(-sig1*dnorm(-mu1/sig1)+mu1*pnorm(-mu1/sig1))
    }
    
    tot <- tot + abs((LHS(Inf)-LHS(0))-RHS(0))
  }
  tot/n_sim
}

check_uv(1000)
