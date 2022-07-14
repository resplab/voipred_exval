library(tidyverse)
library(pROC)
library(mvtnorm)
library(sqldf)
library(tidyverse)
source("R/core.R")

main <- function(n_sim=10)
{
  library(foreach)
  library(doParallel)
  num_cores <- detectCores()-2
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  out <- foreach(i=1:n_sim) %dopar%{
    library(pROC)
    library(mvtnorm)
    source('R/core.R')
    source('R/perturbation_sim_HL.R')
    # cat(i,"\n")
    tmp <- perturbed_scenarios_intercept_noise(seed=i)
    write_rds(tmp,paste0("perturb_sim_files/perturb_",i,".rds"))
    tmp
  }
  stopCluster(cl)
  out <- do.call(rbind,out) %>% 
    as.data.frame()
  colnames(out)[1:4] <- c("sample_size","event_p","noise_sd","c_intercept")
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

# res <- main(n=100)
# write_rds(res,'perturb_sim_results.rds')

res <- read_rds("perturb_sim_results.rds") %>% 
  filter(sample_size>=200)

processed_res <- process_results()

types <- c("ob",'bb','as')

tidy_df <- function(df,types,index){
  tmp <- NULL
  for(i in 1:length(types)){
    type <- types[i]
    tmp[[i]] <-   df %>% 
      mutate(type=type) %>% 
      select(index,type,contains(type)) %>% 
      pivot_longer(cols=starts_with("EVPI"),names_to="threshold",values_to="EVPI") %>% 
      mutate(threshold=parse_number(threshold)/10)
  }
  do.call(rbind,tmp)

}

df_sample_size <- processed_res %>% 
  filter(c_intercept==0) %>% 
  select(sample_size,starts_with("EVPI"))

df_sample_size <- tidy_df(df_sample_size,types,index=1)

fig.size <- 1.3
alpha.size <- 0.7
ggplot(data=df_sample_size %>% 
         mutate(type = case_when(type=='as' ~'Asymptotic',
                                 type == 'bb' ~ "Bayesian",
                                 type == 'ob' ~ "Ordinary")),
       aes(y=EVPI,x=sample_size,colour=type))+
  
  geom_line(alpha=alpha.size,size=fig.size)+
  geom_point()+
  theme_classic() +
  facet_grid(.~threshold) +
  xlab("Sample size") +
  ylab("EVPI") +
  theme(legend.position = 'top',
        legend.title=element_blank(),
        text=element_text(size=18))

df_int <-  processed_res %>% 
  filter(!is.na(c_intercept)) %>% 
  select(sample_size,c_intercept,starts_with("EVPI"))

df_int <- tidy_df(df_int,types,index=c(1,2))

ggplot(data=df_int %>% 
         rename(intercept=c_intercept) %>% 
         mutate(type = case_when(type=='as' ~'Asymptotic',
                                 type == 'bb' ~ "Bayesian",
                                 type == 'ob' ~ "Ordinary")) %>% 
         mutate(sample_size=as.factor(sample_size)),
       aes(y=EVPI,x=intercept,color=sample_size))+
  geom_line(alpha=alpha.size,size=2)+
  geom_point()+
  theme_classic() +
  facet_grid(threshold~type) +
  xlab("Calibration intercept") +
  ylab("EVPI") +
  theme(legend.position = 'top',
        legend.title=element_blank(),
        text=element_text(size=18))


df_slope <-  processed_res %>% 
  filter(!is.na(noise_sd) & is.na(c_intercept)) %>% 
  select(sample_size,noise_sd,starts_with("EVPI"))

df_slope <- tidy_df(df_slope,types,index=c(1,2))

ggplot(data=df_slope %>% 
         rename(slope=noise_sd) %>% 
         mutate(type = case_when(type=='as' ~'Asymptotic',
                                 type == 'bb' ~ "Bayesian",
                                 type == 'ob' ~ "Ordinary")) %>% 
         mutate(sample_size=as.factor(sample_size)),
       aes(y=EVPI,x=slope,color=sample_size))+
  geom_line(alpha=alpha.size,size=2)+
  geom_point()+
  theme_classic() +
  facet_grid(threshold~type) +
  xlab("Slope") +
  ylab("EVPI") +
  theme(legend.position = 'top',
        legend.title=element_blank(),
        text=element_text(size=18))

