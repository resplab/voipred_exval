library(pROC)
library(mvtnorm)
library(sqldf)
source("R/core.R")


perturbed_scenario <- function(sample_size,  true_model = c(-1.55,0.77), cintercept=NA, cslope=NA, max_n_inner_sim=1000, zs=0.1, seed=NULL)
{
  #cat("Doing ", "Sample_size:",sample_size, "cintercept:", cintercept, "cslope:",  cslope,"\n")
  if(!is.null(seed))  set.seed(seed)
  
  sample <- as.data.frame(rmvnorm(sample_size, true_model*0, diag(length(true_model))))
  sample[,1] <- 1
  p <- as.vector(1/(1+exp(-as.matrix(sample)%*%as.matrix(true_model))))
  Y <- rbinom(sample_size, 1, p)
  sample$Y <- Y
  prev <- mean(Y)

  logit_p <- log(p/(1-p))
  
  if(is.na(cslope))
  {
    if(is.na(cintercept)) cintercept <- 0
    
    intercept <- true_model[1]-cintercept
    slope <- true_model[2]
  }
  else
  {
    slope <- true_model[2]/cslope
    if(is.na(cintercept)) #If not provided it means set it such that it is still mean-calibrated.
    {
      intercept <- my_optimize(o, interval=c(-3,3), b1=true_model[2]/cslope, prev=0.2)
    }
    else
    {
      intercept <- true_model[1]-cintercept
      stop("Error: this part is not implemented yet. You should not be here!")
    }
  }
  
  pred_model <- c(intercept, slope)
  pi <- as.vector(1/(1+exp(-as.matrix(sample[,1:2])%*%as.matrix(pred_model))))
  
  logit_pi <- log(pi/(1-pi))
  check <- unname(coefficients(glm(sample$Y~logit_pi, family=binomial())))
  check_mean <- mean(pi)
  
  auc <- roc(Y, pi, quiet=T)$auc
  #message("AUC is", auc)
  
  sample$pi <- pi
  
  res <- voi_ex_glm(model=NULL, val_data=sample, method = "bootstrap", zs = zs)
  res2 <- dca(sample, zs=zs)
  
  c(EVPIv=res[,'EVPIv'], NB_model=res2[,'NB_model'], NB_all=res2[,'NB_all'],  auc=auc, prev=prev, check_intercept=check[1], check_slope=check[2], check_mean=check_mean)
}




perturbed_scenarios <- function(sample_sizes = c(125, 250, 500, 1000,2000),
                                 cintercepts = seq(-0.8,0.8,by=0.1),
                                 cslopes = c(0.6,0.8,1,1.2,1.4),
                                 zs = seq(0.1,0.9,by=0.1),
                                 seed=NULL)
{
  out <- NULL #data.frame("sample_size"=integer(),"event_p"=double(), "noise_sd"=double(),"bias_OR"=double(),"voi_1"=double(), "voi_2"=double(),"voi_3"=double(),"voi_4"=double(),  "auc"=double())
  for(sample_size in sample_sizes)
  {
    if(length(cintercepts)>0)
    {
      for(cintercept in cintercepts)
      {
        res <- perturbed_scenario(sample_size, cintercept = cintercept, zs = zs, seed=seed)
        out <- rbind(out,c(sample_size, cintercept, NA, res))
      }
    }
    if(length(cslopes)>0)
    {
      for(cslope in cslopes)
      {
        res <- perturbed_scenario(sample_size, cslope = cslope, zs = zs, seed=seed)
        out <- rbind(out,c(sample_size, NA, cslope, res))
      }
    }
  }
  
  colnames(out) <- c("sample_size","cintercept", "cslope", names(res))
  
  out
}




main <- function(n_sim=10)
{
  optimize_memory <<- optimize_memory[0,]
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
  x<-sqldf("SELECT COUNT(*) AS N, sample_size, cintercept, cslope,  
        AVG(EVPIv1) AS EVPIv1, AVG(EVPIv2) AS EVPIv2, AVG(EVPIv3) AS EVPIv3,  
        SQRT(VARIANCE(EVPIv1)/COUNT(*)) AS seEVPIv1, SQRT(VARIANCE(EVPIv2)/COUNT(*)) AS seEVPIv2, SQRT(VARIANCE(EVPIv3)/COUNT(*)) AS seEVPIv3, 
        AVG(NB_model1) AS NB_model1, AVG(NB_model2) AS NB_model2, AVG(NB_model3) AS NB_model3, 
        AVG(NB_all1) AS NB_all1, AVG(NB_all2) AS NB_all2, AVG(NB_all3) AS NB_all3, 
        AVG(auc) AS auc, AVG(prev) AS prev,
        AVG(check_intercept) AS check_intercept, AVG(check_slope) AS check_slope, AVG(check_mean) AS check_mean
        FROM res GROUP BY sample_size, cintercept, cslope")
  
  x
}






f <- function(x, b0, b1)
{
  dnorm(x)/(1+exp(-(b0+b1*x)))
}

o <- function(b0, b1, prev)
{
  (integrate(f,-5,5, b0, b1)$value - prev)^2
}


optimize_memory <- data.frame(b1=double(), prev=double(), val=double())

my_optimize <- function (o, interval, b1, prev)
{
  i <- which(optimize_memory$b1==b1 & optimize_memory$prev==prev)
  if(length(i)>0) return(optimize_memory[i,]$val)
  res <- optimize(o, interval=interval, b1=b1, prev=prev)
  optimize_memory[nrow(optimize_memory) + 1,] <<- c(b1, prev, res$minimum)
  message("Added to memory:",c(b1, prev, res$minimum))
  res$minimum
}

# look <- main(n_sim=100)

library(tidyverse)

# saveRDS(look,"perturbation_sim_V5_HL.rds")
res <- read_rds("perturbation_sim_V5_HL.rds")

results <- res %>% 
  filter(!is.na(cintercept)) %>% 
  group_by(sample_size,cintercept) %>% 
  summarise(across(starts_with("EVPI"),mean))

df_int <- results %>% 
  rename(intercept=cintercept) %>% 
  pivot_longer(cols=starts_with("EVPI"),
               names_to="threshold",
               values_to="EVPI") %>% 
  mutate(threshold=parse_number(threshold)/10)
  
ggplot(data=df_int %>% 
         mutate(intercept=as.factor(intercept)),aes(y=EVPI,x=sample_size,color=intercept))+
  geom_line()+
  geom_point()+
  theme_classic() +
  facet_grid(.~threshold) +
  xlab("Sample size") +
  ylab("EVPI") +
  theme(legend.position = 'top')

ggplot(data=df_int %>% 
         filter(threshold <= 0.5) %>% 
         mutate(sample_size=as.factor(sample_size)),aes(y=EVPI,x=intercept,color=sample_size))+
  geom_line()+
  geom_point()+
  theme_classic() +
  facet_grid(.~threshold) +
  xlab("Intercept") +
  ylab("EVPI") +
  theme(legend.position = 'top')

library(plotly)
df_gg <- df_int %>% 
  group_by(sample_size) %>% 
  group_split(.) %>% 
  lapply(.,function(x){
    x %>% 
      select(-sample_size) %>% 
      filter(threshold <= 0.5) %>% 
      pivot_wider(names_from=intercept,values_from=EVPI)
  })

titles <- res$sample_size %>% 
  unique()

D3plotter <- function(index){
  look <- df_gg[[index]]
  look.matrix <- as.matrix(look)
  look.matrix <- look.matrix[,-1]
  plot_ly(z=look.matrix,
          y=look$threshold,
          x=seq(-0.8,0.8,by=0.1),
          ) %>% add_surface() %>% 
    layout(title =titles[index] )
  
}

D3plotter(2)
