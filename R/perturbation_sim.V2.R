#Simple change of intercept and slope, with threshold set at prevalence

library(pROC)
library(mvtnorm)
library(sqldf)


perturbed_scenario <- function(sample_size,  true_model = c(-2,1), event_p=NA, noise_sd=NA, bias_OR=NA, n_out_sim=200, max_n_inner_sim=1000, zs=0.1, seed=NULL)
{
  cat("Doing ", "Sample_size:", sample_size, "true model:", true_model, "event_p:", event_p, "noise_sd:",  noise_sd, "bias_OR:", bias_OR)
  if(!is.null(seed))  set.seed(seed)
  
  while(T)
  {
    sample <- as.data.frame(rmvnorm(sample_size, true_model*0, diag(length(true_model))))
    sample[,1] <- 1
    p <- as.vector(1/(1+exp(-as.matrix(sample)%*%as.matrix(true_model))))
    Y <- rbinom(sample_size, 1, p)
    if(sum(Y)>0) break;
  }
  
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
    weighted.mean(Y,w = weights)
    message("p(Y) in the sample is ", prev)
  }
  else
  {
    weights <- NULL
    prev <- mean(p)
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
  
  sample$pi <- pi
  
  auc <- roc(Y, pi, quiet=T)$auc
  #message("AUC is", auc)
  
  sample$pi <- pi; sample$Y <- Y
  res <- voi_ex_glm(model=NULL, val_data=sample, method = "bootstrap", zs = mean(Y), weights = weights)
  res2 <- dca(sample, zs=mean(Y))
  
  c(EVPIv=res[,'EVPIv'], NB_model=res2$NB_model, NB_all=res2$NB_all,  auc=auc, prev=prev)
}




perturbed_scenarios <- function(seed=NULL)
{
  out <- NULL #data.frame("sample_size"=integer(),"event_p"=double(), "noise_sd"=double(),"bias_OR"=double(),"voi_1"=double(), "voi_2"=double(),"voi_3"=double(),"voi_4"=double(),  "auc"=double())
  sample_sizes <- c(125,250,500,1000)
  #event_ps <- c(NA,0.075, 0.15, 0.3, 0.5)
  #noise_sds <- c(1/3, 2/3, 1, 3/2)
  #bias_ORs <- c(1/2, 3/4 , 4/3 , 2)
  intercepts <- c(-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0)
  slopes <- c(0.25,0.5,1,1.5,2)
  zs <- c(0.05,0.1,0.15,0.20)
  for(sample_size in sample_sizes)
  {
    # for(intercepts in event_ps)
    # {
    #   res <- perturbed_scenario(sample_size, event_p = event_p, zs = zs, seed=seed)
    #   out <- rbind(out,c(sample_size, event_p, NA, NA, res))
    # }
    # for(noise_sd in noise_sds)
    # {
    #   res <- perturbed_scenario(sample_size, noise_sd = noise_sd, zs = zs, seed=seed)
    #   out <- rbind(out,c(sample_size, NA, noise_sd, NA, res))
    # }
    # for(bias_OR in bias_ORs)
    # {
    #   res <- perturbed_scenario(sample_size, bias_OR = bias_OR, zs = zs, seed=seed)
    #   out <- rbind(out,c(sample_size, NA, NA, bias_OR, res))
    # }
    
    if(length(intercepts)>0)
    {
      for(intercept in intercepts)
      {
        res <- perturbed_scenario(sample_size, true_model=c(intercept,1), zs = zs, seed=seed)
        out <- rbind(out,c(sample_size, intercept, NA, res))
      }
    }
    if(length(slopes)>0)
    {
      for(slope in slopes)
      {
        res <- perturbed_scenario(sample_size, true_model=c(-2, slope), zs = zs, seed=seed)
        out <- rbind(out,c(sample_size, NA, slope, res))
      }
    }
  }
  
  colnames(out) <- c("sample_size","intercept", "slope", names(res))
  
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
  x <- sqldf("SELECT COUNT(*) AS N, sample_size, intercept, slope, 
        AVG(EVPIv) AS EVPIv, 
        AVG(NB_model) AS NB_model,
        AVG(NB_all) AS NB_all,
        AVG(auc) AS auc,
        AVG(prev) AS prev
        FROM res GROUP BY sample_size, intercept, slope")
}


