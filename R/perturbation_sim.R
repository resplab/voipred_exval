library(pROC)
library(mvtnorm)
library(sqldf)


perturbed_scenario <- function(sample_size,  true_model = c(-2,1), event_p=NA, noise_sd=NA, bias_OR=NA, n_out_sim=200, max_n_inner_sim=1000, zs=0.1, seed=NULL)
{
  cat("Doing ", "Sample_size:",sample_size, "event_p:", event_p, "noise_sd:",  noise_sd, "bias_OR:", bias_OR)
  if(!is.null(seed))  set.seed(seed)
  
  sample <- as.data.frame(rmvnorm(sample_size, true_model*0, diag(length(true_model))))
  sample[,1] <- 1
  p <- as.vector(1/(1+exp(-as.matrix(sample)%*%as.matrix(true_model))))
  Y <- rbinom(sample_size, 1, p)
  
  if(!is.na(event_p))
  {
    cases <- which(Y==1)
    controls <- which(Y==0)
    P <- length(cases) / sample_size
    Q <- event_p
    weights <- rep(NA, sample_size)
    weights[cases] <- Q*(1-P)/P/(1-Q)
    weights[controls] <- 1
    weights <- weights / sum(weights) * sample_size
    prev <- weighted.mean(Y,w = weights)
    message("p(Y) in the sample is ", prev)
  }
  else
  {
    prev <- mean(p)
    weights <- NULL
    message("p(Y) in the sample is ", mean(Y))
  }
  
  pi <- p
  
  #bias and intercept of the proposed model
  if(!is.na(noise_sd))
  {
    log_lins <- log(pi/(1-pi)) + rnorm(length(pi),0, noise_sd)
    #Make sure the model has the correct intercept
    tmp <- glm(Y~offset(log_lins), family = binomial(link="logit"))
    log_lins <- log_lins + coefficients(tmp)
    pi <- as.vector(1/(1+exp(-log_lins)))
  }
  
  if(!is.na(bias_OR))
  {
    log_lins <- log(pi/(1-pi)) + log(bias_OR)
    pi <- as.vector(1/(1+exp(-log_lins)))
  }
  
  auc <- roc(Y, pi, quiet=T)$auc
  #message("AUC is", auc)
  
  sample$pi <- pi; sample$Y <- Y
  res <- voi_ex_glm(model=NULL, val_data=sample, method = "bootstrap", zs = zs, weights = weights)
  res2 <- dca(sample, zs=zs)
  
  c(EVPIv=res[,'EVPIv'], NB_model=res2[,'NB_model'], NB_all=res2[,'NB_all'],  auc=auc, prev=prev)
}




perturbed_scenarios <- function(sample_sizes = c(125, 250, 500, 1000),
                                 event_ps = c(0.05,0.1,0.25,0.5),
                                 noise_sds = c(0.25,0.5,0.75,1),
                                 bias_ORs = c(1,1.5,2.0,2.5,3,3.5,4),
                                 zs = c(0.05,0.1,0.15,0.20),
                                 seed=NULL)
{
  out <- NULL #data.frame("sample_size"=integer(),"event_p"=double(), "noise_sd"=double(),"bias_OR"=double(),"voi_1"=double(), "voi_2"=double(),"voi_3"=double(),"voi_4"=double(),  "auc"=double())
  for(sample_size in sample_sizes)
  {
    if(length(event_ps)>0)
    {
      for(event_p in event_ps)
      {
        res <- perturbed_scenario(sample_size, event_p = event_p, zs = zs, seed=seed)
        out <- rbind(out,c(sample_size, event_p, NA, NA, res))
      }
    }
    if(length(noise_sds)>0)
    {
      for(noise_sd in noise_sds)
      {
        res <- perturbed_scenario(sample_size, noise_sd = noise_sd, zs = zs, seed=seed)
        out <- rbind(out,c(sample_size, NA, noise_sd, NA, res))
      }
    }
    if(length(bias_ORs)>0)
    {
      for(bias_OR in bias_ORs)
      {
        res <- perturbed_scenario(sample_size, bias_OR = bias_OR, zs = zs, seed=seed)
        out <- rbind(out,c(sample_size, NA, NA, bias_OR, res))
      }
    }
  }
  
  colnames(out) <- c("sample_size","event_p", "noise_sd","bias_OR", names(res))
  
  out
}




main <- function(n_sim=10)
{
  out <- NULL
  for(i in 1:n_sim)
  {
    out <- rbind(out,perturbed_scenarios(seed=i))
  }
  
  out <- as.data.frame(out)
  
  res <<- out 
  
  out
}



process_results <- function()
{
  x<-sqldf("SELECT COUNT(*) AS N, sample_size, event_p, noise_sd, bias_OR, 
        AVG(EVPIv1) AS EVPIv1, AVG(EVPIv2) AS EVPIv2, AVG(EVPIv3) AS EVPIv3, AVG(EVPIv4) AS EVPIv4, 
        AVG(NB_model1) AS NB_model1, AVG(NB_model2) AS NB_model2, AVG(NB_model3) AS NB_model3, AVG(NB_model4) AS NB_model4,
        AVG(NB_all1) AS NB_all1, AVG(NB_all2) AS NB_all2, AVG(NB_all3) AS NB_all3, AVG(NB_all4) AS NB_all4,
        AVG(auc) AS auc, AVG(prev) AS prev
        FROM res GROUP BY sample_size, event_p, noise_sd, bias_OR")
  
  x
}


