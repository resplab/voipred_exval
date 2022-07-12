library(pROC)
library(mvtnorm)
library(sqldf)
library(tidyverse)


perturbed_scenario <- function(sample_size,  true_model = c(-1.55,0.77), event_p=NA, noise_sds=NA, c_intercepts=NA, n_out_sim=200, max_n_inner_sim=1000, zs=0.1, seed=NULL)
{
  #cat("Doing ", "Sample_size:",sample_size, "event_p:", event_p, "noise_sd:",  noise_sd, "c_intercept:", c_intercept)
  if(!is.null(seed))  set.seed(seed)
  
  sample <- as.data.frame(rmvnorm(sample_size, true_model*0, diag(length(true_model))))
  sample[,1] <- 1
  p <- as.vector(1/(1+exp(-as.matrix(sample)%*%as.matrix(true_model))))
  Y <- rbinom(sample_size, 1, p)
  
  prev <- mean(p)
  weights <- NULL
  
  pi_original <- p
  
  if(is.na(noise_sds)){
    
    if(is.na(c_intercepts)){
      error("Input something")
    } else{
      df_test <- data.frame(value=c_intercepts,type='cintercept')
    } 

  } else{
    
    if(is.na(c_intercepts)){
      df_test <- data.frame(value=noise_sds,type='sd')
    } else{
      df_test <- rbind(data.frame(value=c_intercepts,type='cintercept'),
                       data.frame(value=noise_sds,type='sd'))
    }
  }
  
  results <- NULL
  
  for( i in 1:nrow(df_test)){
    #bias and intercept of the proposed model
    tmp_type <- df_test$type[i]
    
    if(tmp_type=="sd") {
      noise_sd <- df_test$value[i]
      log_lins <- log(pi_original/(1-pi_original)) + rnorm(length(pi_original),0, noise_sd)
      #Make sure the model has the correct intercept
      tmp <- glm(Y~offset(log_lins), family = binomial(link="logit"))
      log_lins <- log_lins + coefficients(tmp)
      pi <- as.vector(1/(1+exp(-log_lins)))

    }
    
    if(tmp_type=='cintercept'){
      c_intercept <- df_test$value[i]
      log_lins <- log(pi_original/(1-pi_original)) - c_intercept
      pi <- as.vector(1/(1+exp(-log_lins)))
    }
    
    predprev <- mean(pi)
    
    auc <- roc(Y, pi, quiet=T)$auc

    sample$pi <- pi; sample$Y <- Y
    sample$logit_pi <- log(pi/(1-pi))
    
    check <- unname(coefficients(glm(Y~logit_pi, data=sample, family=binomial())))
    
    res_ob <- voi_ex_glm(model=NULL, val_data=sample, method = "bootstrap", zs = zs, weights = weights)
    res_bb <- voi_ex_glm(model=NULL, val_data=sample, method = "bootstrap", Bayesian_bootstrap = T, zs = zs, weights = weights)
    res_as <- voi_ex_glm(model=NULL, val_data=sample, method = "asymptotic", zs = zs, weights = weights)
    
    res_dca <- dca(sample, zs=zs)
    
    res <-     c(EVPIv_ob=res_ob[,'EVPIv'], 
                 EVPIv_bb=res_bb[,'EVPIv'], 
                 EVPIv_as=res_as[,'EVPIv'], 
                 NB_model=res_dca[,'NB_model'], 
                 NB_all=res_dca[,'NB_all'], 
                 p_useful=res_ob$p_useful, 
                 auc=auc, 
                 prev=prev, 
                 predprev=predprev, 
                 check_intercept=check[1], 
                 check_slope=check[2])
    
    if(tmp_type == 'sd'){
      results[[i]] <- c(sample_size,NA,noise_sd,NA,res)
    }
    
    if (tmp_type =='cintercept'){
      results[[i]] <- c(sample_size, NA, NA, c_intercept, res)
      
    }

    
  }
  
  results %>% 
    do.call(rbind,.)
}




perturbed_scenarios_intercept_noise <- function(
                                sample_sizes = c(125, 250, 500, 1000,2000),
                                noise_sds = c(0,0.15,0.30,0.45,0.60),
                                c_intercepts = c(-3.2, -1.6, -0.8, -0.4, 0, 0.4, 0.8, 1.6, 3.2),
                                zs = c(0.1,0.2,0.3), #seq(0.25,0.75,length=21),
                                seed=NULL){
  out <- NULL #data.frame("sample_size"=integer(),"event_p"=double(), "noise_sd"=double(),"c_intercept"=double(),"voi_1"=double(), "voi_2"=double(),"voi_3"=double(),"voi_4"=double(),  "auc"=double())
  
  for(sample_size in sample_sizes) {
    res <- perturbed_scenario(sample_size, c_intercepts = c_intercepts, noise_sds = noise_sds, zs = zs, seed=seed)
    out <- rbind(out,res)
  }
  
  out
}

process_results <- function()
{
  x<-sqldf("SELECT COUNT(*) AS N, sample_size, event_p, noise_sd, c_intercept, 
        AVG(evpiv_ob1) AS EVPIv_ob1, AVG(evpiv_ob2) AS EVPIv_ob2, AVG(evpiv_ob3) AS EVPIv_ob3,
        AVG(evpiv_bb1) AS EVPIv_bb1, AVG(evpiv_bb2) AS EVPIv_bb2, AVG(evpiv_bb3) AS EVPIv_bb3,
        AVG(evpiv_as1) AS EVPIv_as1, AVG(evpiv_as2) AS EVPIv_as2, AVG(evpiv_as3) AS EVPIv_as3,
        SQRT(VARIANCE(evpiv_ob1)/count(*)) AS se1, SQRT(VARIANCE(evpiv_ob1)/count(*)) AS se2, SQRT(VARIANCE(evpiv_ob1)/count(*)) AS se3, 
        AVG(p_useful1) AS p_useful1, AVG(p_useful2) AS p_useful2, AVG(p_useful3) AS p_useful3, 
        SQRT(VARIANCE(NB_model1-MAX(NB_all1,0))) AS sddNB1, SQRT(VARIANCE(NB_model2-MAX(NB_all2,0))) AS sddNB2, SQRT(VARIANCE(NB_model3-MAX(NB_all3,0))) AS sddNB3, 
        AVG(NB_model1) AS NB_model1, AVG(NB_model2) AS NB_model2, AVG(NB_model3) AS NB_model3, 
        AVG(NB_all1) AS NB_all1, AVG(NB_all2) AS NB_all2, AVG(NB_all3) AS NB_all3, 
        AVG(auc) AS auc, AVG(prev) AS prev, AVG(predprev) AS predprev, 
        AVG(check_intercept) AS check_intercept, AVG(check_slope) AS check_slope
        FROM res GROUP BY sample_size, event_p, noise_sd, c_intercept")
  
  x
}

