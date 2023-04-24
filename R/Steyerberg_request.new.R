# title: "EVPI for external validation"
# author: "Mohsen Sadatsafavi"
# date: "2023.01.11"


source('core.R')

library(voipred)
library(sqldf)
set.seed(123) 

settings <- list(
  val_sample_size = 500,
  n_sim=10000,
  zs=c(0.01, 0.02, 0.05, 0.1)
)

data(gusto)
gusto$kill <- (as.numeric(gusto$Killip)>1)*1
gusto$Y <- gusto$day30


# data_us <- gusto[gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]
# data_other <- gusto[!gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]

tmp <- sqldf("SELECT regl, COUNT(*) AS N, AVG(Y) AS mean_Y FROM gusto GROUP BY regl ORDER BY mean_Y")
high_reg <- tmp$regl[9:16]
data_us <- gusto[gusto$regl %in% high_reg,]
data_other <- gusto[!gusto$regl %in% high_reg,]
# y1 <- sum(data_us$Y)
# y2 <- sum(data_other$Y)
# data_us <- data_us[1:min(which(cumsum(data_us$Y)==y2)),]

n_us <- dim(data_us)[1]
n_other <- dim(data_other)[1]


pb <- progress_bar$new(total=settings$n_sim)
out <- data.frame()
for(i in 1:settings$n_sim)
{
  pb$tick()
  conscripts_us <- sample(1:n_us, size=500)
  val_data_us <- data_us[conscripts_us,]
  conscripts_other <- sample(1:n_other, size=500)
  val_data_other <- data_other[conscripts_other,]
  
  dev_data_us <- data_us[-conscripts_us,]
  dev_data_other <- data_other#[-conscripts_other,]
  
  model_us <- glm(Y ~ age + miloc + pmi + kill + pmin(sysbp,100) + pulse, data=dev_data_us, family=binomial(link="logit"))
  model_other <- glm(Y ~ age + miloc + pmi + kill + pmin(sysbp,100) + pulse, data=dev_data_other, family=binomial(link="logit"))
  
  pi_us_us <- as.vector(predict(model_us, type="response", newdata=val_data_us))
  pi_us_other <- as.vector(predict(model_other, type="response", newdata=val_data_us))
  pi_other_us <- predict(model_us, type="response", newdata=val_data_other)
  pi_other_other <- predict(model_other, type="response", newdata=val_data_other)
  #
  # NB_model_us_us <- mean((pi_us_us>settings$z)*(val_data_us$Y - (1-val_data_us$Y)*settings$z/(1-settings$z)))
  # NB_model_us_other <- mean((pi_us_other>settings$z)*(val_data_us$Y - (1-val_data_us$Y)*settings$z/(1-settings$z)))
  # # NB_model_other_us[i] <- mean((val_data_other$pi_us>z)*(val_data_other$Y - (1-val_data_other$Y)*z/(1-z)))
  # NB_model_other_other[i] <- meansettings$((val_data_other$pi_other>z)*(val_data_other$Y - (1-val_data_other$Y)*z/(1-z)))
  # 
  #NB_all_us <- mean((val_data_us$Y - (1-val_data_us$Y)*settings$z/(1-settings$z)))
  #NB_all_other <- mean((val_data_other$Y - (1-val_data_other$Y)*z/(1-z)))
  
  #if (is.na(NB_model_us_other-NB_all_other)) browser()
  
  dt_us_us <- cbind(val_data_us, pi=pi_us_us)
  dt_us_other <- cbind(val_data_us, pi=pi_us_other)
  dt_other_us <- cbind(val_data_other, pi=pi_other_us)
  dt_other_other <- cbind(val_data_other, pi=pi_other_other)

  out <- data.frame(rbind(out,
          cbind(scenario="res_us_us", method="ob",
                #dNBh=NB_model_us_us-NB_all_us,
                voi_ex_glm(model_us,dt_us_us,method='bootstrap',Bayesian_bootstrap = F,z = settings$zs, n_sim=1000)),
          cbind(scenario="res_us_other", method="ob",
                #dNBh=NB_model_us_other-NB_all_us,
                voi_ex_glm(model_other,dt_us_other,method='bootstrap',Bayesian_bootstrap = F,z = settings$zs, n_sim=1000))
        ))
  out <- data.frame(rbind(out,
                          cbind(scenario="res_us_us", method="bb", 
                                #dNBh=NB_model_us_us-NB_all_us, 
                                voi_ex_glm(model_us,dt_us_us,method='bootstrap',Bayesian_bootstrap = T,z = settings$zs, n_sim=1000)),
                          cbind(scenario="res_us_other", method="bb",
                                #dNBh=NB_model_us_other-NB_all_us, 
                                voi_ex_glm(model_other,dt_us_other,method='bootstrap',Bayesian_bootstrap = T,z = settings$zs, n_sim=1000))
  ))
  out <- data.frame(rbind(out,
                                cbind(scenario="res_us_us", method="as",
                                      #dNBh=NB_model_us_us-NB_all_us,
                                      voi_ex_glm(model_us,dt_us_us,method='asymptotic',Bayesian_bootstrap = F,z = settings$zs, n_sim=1000)),
                                cbind(scenario="res_us_other", method="as",
                                      #dNBh=NB_model_us_other-NB_all_us,
                                      voi_ex_glm(model_other,dt_us_other,method='asymptotic',Bayesian_bootstrap = F,z = settings$zs, n_sim=1000))
        ))
    
  
  
  
  
  
  out <- data.frame(rbind(out,
                          cbind(scenario="res_other_us", method="ob",
                                #dNBh=NB_model_us_us-NB_all_us,
                                voi_ex_glm(model_us,dt_other_us,method='bootstrap',Bayesian_bootstrap = F,z = settings$zs, n_sim=1000)),
                          cbind(scenario="res_other_other", method="ob",
                                #dNBh=NB_model_us_other-NB_all_us,
                                voi_ex_glm(model_other,dt_other_other,method='bootstrap',Bayesian_bootstrap = F,z = settings$zs, n_sim=1000))
  ))
  out <- data.frame(rbind(out,
                          cbind(scenario="res_other_us", method="bb",
                                #dNBh=NB_model_us_us-NB_all_us,
                                voi_ex_glm(model_us,dt_other_us,method='bootstrap',Bayesian_bootstrap = T,z = settings$zs, n_sim=1000)),
                          cbind(scenario="res_other_other", method="bb",
                                #dNBh=NB_model_us_other-NB_all_us,
                                voi_ex_glm(model_other,dt_other_other,method='bootstrap',Bayesian_bootstrap = T,z = settings$zs, n_sim=1000))
  ))
  out <- data.frame(rbind(out,
                          cbind(scenario="res_other_us", method="as",
                                #dNBh=NB_model_us_us-NB_all_us,
                                voi_ex_glm(model_us,dt_other_us,method='asymptotic',Bayesian_bootstrap = F,z = settings$zs, n_sim=1000)),
                          cbind(scenario="res_other_other", method="as",
                                #dNBh=NB_model_us_other-NB_all_us,
                                voi_ex_glm(model_other,dt_other_other,method='asymptotic',Bayesian_bootstrap = F,z = settings$zs, n_sim=1000))
  ))
  
  
  # out <- data.frame(rbind(out,
  #   cbind(scenario="res_other_us", method="ob", 
  #         #dNBh=NB_model_us_us-NB_all_us, 
  #         voi_ex_glm(model_us,dt_other_us,method='bootstrap',Bayesian_bootstrap = F,z = settings$zs, n_sim=1000)),
  #   cbind(scenario="res_other_other", method="ob",
  #         #dNBh=NB_model_us_other-NB_all_us, 
  #         voi_ex_glm(model_other,dt_other_other,method='bootstrap',Bayesian_bootstrap = F,z = settings$zs, n_sim=1000))
  # ))
  #   
    
    # cbind(scenario="res_us_us", method="bb", 
    #       #dNBh=NB_model_us_us-NB_all_us, 
    #       voi_ex_glm(model_us,dt_us_us,method='bootstrap',Bayesian_bootstrap = T,z = settings$zs, n_sim=1000)),
    # cbind(scenario="res_us_other", method="bb",
    #       #dNBh=NB_model_us_other-NB_all_us, 
    #       voi_ex_glm(model_other,dt_us_other,method='bootstrap',Bayesian_bootstrap = T,z = settings$zs, n_sim=1000)),
    # cbind(scenario="res_us_us", method="as", 
    #       #dNBh=NB_model_us_us-NB_all_us, 
    #       voi_ex_glm(model_us,dt_us_us,method='asymptotic',z = settings$zs, n_sim=1000)),
    # cbind(scenario="res_us_other", method="as",
    #       #dNBh=NB_model_us_other-NB_all_us, 
    #       voi_ex_glm(model_other,dt_us_other,method='asymptotic',z = settings$zs, n_sim=1000))
    # 
    #,
    #cbind(scenario="res_other_us", voi_ex_glm(model_us,dt_other_us,method='bootstrap',Bayesian_bootstrap = T,z = settings$z, n_sim=1000)),
    #cbind(scenario="res_other_other", voi_ex_glm(model_other,dt_other_other,method='bootstrap',Bayesian_bootstrap = T,z = settings$z, n_sim=1000))
  # ))
  #
  
  if(is.null(sum(out$EVPIv))) browser()
  
  #cat('.')
}



#out %>% group_by(scenario) %>% summarise(n=n(), EVPI=mean(EVPIv), p_useful=mean(p_useful))


res <- sqldf("SELECT scenario, z, method,
      COUNT(*) AS n, 
      AVG(ENB_model-max(ENB_all,0)) AS dNB,
      SQRT(VARIANCE(ENB_model-ENB_all)/COUNT(*)) AS se_dNB, 
      AVG(EVPIv) AS EVPI, 
      SQRT(VARIANCE(EVPIv)/COUNT(*)) AS se_EVPI, 
      AVG(p_useful) AS p_useful, 
      SQRT(VARIANCE(p_useful)/COUNT(*)) AS se_p_useful
      FROM out GROUP BY scenario, z, method")

(res$EVPI[1:3]-res$EVPI[4:6])/sqrt(res$se_EVPI[1:3]^2+res$se_EVPI[4:6]^2)

