library(pROC)
library(mvtnorm)
library(sqldf)


perturbed_scenario <- function(sample_size,  true_model = c(-1.55,0.77), n_dev, max_n_inner_sim=1000, zs=0.1, seed=NULL)
{
  #cat("Doing ", "Sample_size:",sample_size, "n_dev:", n_dev, "\n")
  if(!is.null(seed))  set.seed(seed)
  
  sample <- as.data.frame(rmvnorm(sample_size, true_model*0, diag(length(true_model))))
  sample[,1] <- 1
  colnames(sample)=paste0("X",1:length(true_model))
  p <- as.vector(1/(1+exp(-as.matrix(sample)%*%as.matrix(true_model))))
  Y <- rbinom(sample_size, 1, p)
  sample$Y <- Y
  
  prev <- mean(Y)
  
  if(!is.na(n_dev))
  {
    while(T)
    {
      dev_sample <- as.data.frame(rmvnorm(n_dev, true_model*0, diag(length(true_model))))
      colnames(dev_sample)=paste0("X",1:length(true_model))
      dev_sample[,1] <- 1
      dev_p <- as.vector(1/(1+exp(-as.matrix(dev_sample)%*%as.matrix(true_model))))
      dev_Y <- rbinom(n_dev, 1, dev_p)
      if(sum(dev_Y)>2) break;
    }
    dev_sample$Y <- dev_Y
    model <- glm(Y ~ 0 + X1 + X2, data=dev_sample, family = binomial(link="logit"))
    
    logit_pi <- predict(model, newdata = sample)
    pi <- 1/(1+exp(-logit_pi))
    sample$logit_pi <- logit_pi
    sample$pi <- pi
 
    model2 <- glm(Y ~ logit_pi, data=sample, family = binomial(link="logit"))
    tmp <- coefficients(model2)
    intercept <- tmp[1]
    slope <- tmp[2]
  }
  
  auc <- roc(Y, pi, quiet=T)$auc
  #message("AUC is", auc)
  
  res <- voi_ex_glm(model=NULL, val_data=sample, method = "bootstrap", zs = zs)
  res2 <- dca(sample, zs=zs)
  
  c(EVPIv=res[,'EVPIv'], NB_model=res2[,'NB_model'], NB_all=res2[,'NB_all'],  auc=auc, prev=prev, intercept=intercept, slope=slope)
}




perturbed_scenarios <- function(sample_sizes = c(125, 250, 500, 1000),
                                 n_devs = c(50,100,200,400,800),
                                 zs = c(0.1,0.2,0.3),
                                 seed=NULL)
{
  out <- NULL #data.frame("sample_size"=integer(),"event_p"=double(), "noise_sd"=double(),"bias_OR"=double(),"voi_1"=double(), "voi_2"=double(),"voi_3"=double(),"voi_4"=double(),  "auc"=double())
  for(sample_size in sample_sizes)
  {
    for(n_dev in n_devs)
    {
      res <- perturbed_scenario(sample_size, n_dev = n_dev, zs = zs, seed=seed)
      out <- rbind(out,c(sample_size, n_dev, res))
    }
  }
  
  colnames(out) <- c("sample_size","n_dev", names(res))
  
  out
}




main <- function(n_sim=100)
{
  out <- NULL
  for(i in 1:n_sim)
  {
    cat(i,"\n")
    out <- rbind(out,perturbed_scenarios(seed=i))
  }
  
  out <- as.data.frame(out)
  
  res <<- out 
  
  out
}



process_results <- function()
{
  x<-sqldf("SELECT COUNT(*) AS N, sample_size, n_dev, 
        AVG(EVPIv1) AS EVPIv1, AVG(EVPIv2) AS EVPIv2, AVG(EVPIv3) AS EVPIv3, 
        SQRT(VARIANCE(EVPIv1)/COUNT(*)) AS seEVPIv1, SQRT(VARIANCE(EVPIv2)/COUNT(*)) AS seEVPIv2, SQRT(VARIANCE(EVPIv3)/COUNT(*)) AS seEVPIv3, 
        AVG(NB_model1) AS NB_model1, AVG(NB_model2) AS NB_model2, AVG(NB_model3) AS NB_model3, 
        AVG(NB_all1) AS NB_all1, AVG(NB_all2) AS NB_all2, AVG(NB_all3) AS NB_all3, 
        AVG(auc) AS auc, AVG(prev) AS prev,
        AVG([intercept.(Intercept)]) AS intercept, AVG([slope.logit_pi]) AS slope
        FROM res GROUP BY sample_size, n_dev")
  
  x
}


