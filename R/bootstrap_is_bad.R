library(mvtnorm)
p <- (1:50)/100
lp <- log(p/(1-p))

sim <- function()
{
  y <- rbinom(length(p),1,1/(1+exp(-(1+1*lp))))

  reg <- glm(y~lp, family=binomial(link='logit'))
  mu <- unname(coefficients(reg))
  sigma <- vcov(reg)

  n_sim <- 100
  out <- rep(NA,n_sim)
  for(i in 1:n_sim)
  {
    bs <- rmvnorm(1,mu,sigma)
    out[i] <- mean(1/(1+exp(-(bs[1]+bs[2]*lp))))
  }

  return(c(mean(y)*(1-mean(y))/length(p), var(out)))
}


out <- matrix(NA,nrow=100,ncol=2)
for(i in 1:(dim(out)[1]))
{
  out[i,] <- sim()
}

colMeans(out)[2] / colMeans(out)[1]
