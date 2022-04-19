library(mvtnorm)

# exact method ------------------------------------------------------------
mu_max_truncated_bvn <- function(mu1,mu2,sig1,sig2,rho,exact=T){
  sig1 <- sqrt(sig1)
  sig2 <- sqrt(sig2)
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





