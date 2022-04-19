source('./R/core.R')
library(voipred)
library(rms)
library(mvtnorm)
library(sqldf)

set.seed(1)

settings <- list(
  dev_sample_size = Inf,
  val_sample_size = 1000,
  zs=c(0.01,0.02,0.05,0.1),
  max_x_evpi_plot = 1   #For plots of evpi 
)



data(gusto)
gusto$kill <- (as.numeric(gusto$Killip)>1)*1
gusto$Y <- gusto$day30
data_us <- gusto[gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]
data_other <- gusto[!gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]

if(is.infinite(settings$dev_sample_size))
{
  dev_data <- data_other
}else
{
  dev_data <- data_other[sample(1:(dim(data_other)[1]),settings$dev_sample_size,F),]
}

if(is.infinite(settings$val_sample_size))
{
  val_data <- data_us
}else
{
  val_data <- data_us[sample(1:(dim(data_us)[1]),settings$val_sample_size,F),]  #data is for external validation
}

model <- glm(Y ~ age + miloc + pmi + kill + pmin(sysbp,100) + lsp(pulse,50), data=dev_data, family=binomial(link="logit"))
#model <- glm(Y ~ age + miloc + pmi + kill, data=dev_data, family=binomial(link="logit"))
pi <- predict(model, type="response", newdata=val_data)
val_data$pi <- pi
val_data$logit_pi <- log(pi/(1-pi))

zs <- sort(c(as.vector(quantile(pi,(0:99)/100)),settings$zs))





#res <- voi_ex_glm(model, val_data, method="asymptotic", zs = 0.1)
#lines(zs,res$EVPIv, type='l')



sim <- function(n_sim=10)
{
  for(i in 1:n_sim)
  {
    set.seed(i_sim)
    cat('.')
    this_data <- data_us[sample(1:(dim(data_us)[1]),settings$val_sample_size,F),]  #data is for external validation
    this_data$pi <- predict(model,newdata = this_data, type='response')
    saveRDS(this_data,"val_data.RDS")
    #this_model <- glm(Y ~ age + miloc + pmi + kill + pmin(sysbp,100) + lsp(pulse,50), data=this_data, family=binomial(link="logit"))
    #pi <- predict(this_model, type="response", newdata=this_data)
    #this_data$pi <- pi
    res_bs <- res_ll <- res_mb <- NULL
    tmp_bs <- voi_ex_glm(model, this_data, method="naive", zs = zs)
    tmp_ll <- voi_ex_glm(model, this_data, method="likelihood", zs = zs)
    tmp_mb <- voi_ex_glm(model, this_data, method="bootstrap", zs = zs)
    plot(zs,tmp_bs$EVPIv, type='l')
    lines(zs,tmp_ll$EVPIv, type='l', col='blue')
    lines(zs,tmp_mb$EVPIv, type='l', col='red')
    res_bs <- rbind(res_bs,tmp_bs)
    res_ll <- rbind(res_ll,tmp_ll)
    res_mb <- rbind(res_mb,tmp_mb)
  }
  list(res_bs=res_bs,res_ll=res_ll,res_mb=res_mb)
}









sim_by_size <- function(n_sim=10, sample_sizes=c(250,500,1000,4000),zs=c(0.01,0.02,0.05,0.1))
{
  out <- data.frame(method=character(), sample_size=integer())
  for(i in 1:length(zs))
  {
    out[paste0('val',i)] <- double()
  }
  
  index <- 1
  
  for(s in 1:length(sample_sizes))
  {
    sample_size <- sample_sizes[s]
    cat(sample_size)
    res_bb <- res_ob <- res_ll <- rep(0,length(zs))
    
    for(i in 1:n_sim)
    {
      cat('.')
      if(is.infinite(sample_size))
      {
        this_data <- data_us
        sample_size <- dim(this_data)[1]
      }
      else
      {
        this_data <- data_us[sample(1:(dim(data_us)[1]),sample_size,F),]  #data is for external validation
      }
      this_data$pi <- predict(model,newdata = this_data, type='response')
      
      tmp <- voi_ex_glm(model, this_data, method="bootstrap", Bayesian_bootstrap = T, zs = zs)
      out[index,'method'] <- "BB"; out[index,'sample_size']<-sample_size; out[index,c('val1','val2','val3','val4')] <- tmp$EVPIv
      index <- index+1
      tmp <- voi_ex_glm(model, this_data, method="bootstrap", Bayesian_bootstrap = F, zs = zs)
      out[index,'method'] <- "OB"; out[index,'sample_size']<-sample_size; out[index,c('val1','val2','val3','val4')] <- tmp$EVPIv
      index <- index+1
      tmp <- voi_ex_glm(model, this_data, method="model_based_bs", Bayesian_bootstrap = T, zs = zs)
      out[index,'method'] <- "likelihood"; out[index,'sample_size']<-sample_size; out[index,c('val1','val2','val3','val4')] <- tmp$EVPIv
      index <- index+1
    }
  }
  out
}
