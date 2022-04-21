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

zs <- sort(c(as.vector(quantile(pi,(1:99)/100)),settings$zs))





#res <- voi_ex_glm(model, val_data, method="asymptotic", zs = 0.1)
#lines(zs,res$EVPIv, type='l')



sim <- function(n_sim=10)
{
  set.seed(1)
  for(i in 1:n_sim)
  {
    cat('.')
    this_data <- data_us[sample(1:(dim(data_us)[1]),settings$val_sample_size,F),]  #data is for external validation
    this_data$pi <- predict(model,newdata = this_data, type='response')
    #this_model <- glm(Y ~ age + miloc + pmi + kill + pmin(sysbp,100) + lsp(pulse,50), data=this_data, family=binomial(link="logit"))
    #pi <- predict(this_model, type="response", newdata=this_data)
    #this_data$pi <- pi
    res_bs <- res_ll <- res_mb <- res_as <- NULL
    tmp_bs <- voi_ex_glm(model, this_data, method="bootstrap", zs = zs)
    #tmp_ll <- voi_ex_glm(model, this_data, method="model_based_ll", zs = zs)
    #tmp_mb <- voi_ex_glm(model, this_data, method="model_based_bs", zs = zs)
    tmp_as <- voi_ex_glm(model, this_data, method="asymptotic", zs = zs)
    
    plot(zs,tmp_bs$EVPIv, type='l', xlim=c(0,0.2))
    #lines(zs,tmp_ll$EVPIv, type='l', col='blue')
    #lines(zs,tmp_mb$EVPIv, type='l', col='red')
    lines(zs,tmp_as$EVPIv, type='l', col='orange')
    res_bs <- rbind(res_bs,tmp_bs)
    #res_ll <- rbind(res_ll,tmp_ll)
    #res_mb <- rbind(res_mb,tmp_mb)
    res_as <- rbind(res_as,tmp_as)
  }
  list(res_bs=res_bs,res_as)
}









sim_by_size <- function(n_sim=100, sample_sizes=c(500,1000,2000,4000,8000,Inf),zs=c(0.01,0.02,0.05,0.1))
{
  set.seed(1)
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
      #set.seed(i)
      cat('.')
      if(is.infinite(sample_size))
      {
        this_data <- data_us
        this_data$pi <- predict(model,newdata = this_data, type='response')
        sample_size <- dim(this_data)[1]
      }
      else
      {
        repeat
        {
          this_data <- data_us[sample(1:(dim(data_us)[1]),sample_size,F),]  #data is for external validation
          this_data$pi <- predict(model,newdata = this_data, type='response')
          if(min(this_data$pi)<min(zs))
          {
            break;
          }
          else
          {
            cat('bad')
          }
        }
      }
      
      bad <- F
      
      tmp <- voi_ex_glm(model, this_data, method="bootstrap", Bayesian_bootstrap = T, zs = zs)
      if(is.null(tmp)) bad <- T
      out[index,'method'] <- "BB"; out[index,'sample_size']<-sample_size; out[index,c('val1','val2','val3','val4')] <- tmp$EVPIv
      index <- index+1
      
      tmp <- voi_ex_glm(model, this_data, method="bootstrap", Bayesian_bootstrap = F, zs = zs)
      if(is.null(tmp)) bad <- T
      out[index,'method'] <- "OB"; out[index,'sample_size']<-sample_size; out[index,c('val1','val2','val3','val4')] <- tmp$EVPIv
      index <- index+1
      # tmp <- voi_ex_glm(model, this_data, method="model_based_bs", Bayesian_bootstrap = T, zs = zs)
      # out[index,'method'] <- "mbbs"; out[index,'sample_size']<-sample_size; out[index,c('val1','val2','val3','val4')] <- tmp$EVPIv
      # index <- index+1
      
      tmp <- voi_ex_glm(model, this_data, method="asymptotic", zs = zs)
      if(is.null(tmp)) bad <- T
      out[index,'method'] <- "asy"; out[index,'sample_size']<-sample_size; out[index,c('val1','val2','val3','val4')] <- tmp$EVPIv
      index <- index+1
      
      if(bad)
      {
        index <- index -3
        i <- i-1
        message("bad")
      }
    }
  }
  
  out
}

#pres <- sqldf("SELECT method, sample_size, AVG(val1) AS val1, AVG(val2) AS val2 , AVG(val3) AS val3 , AVG(val4) AS val4 from res GROUP BY method, sample_size")
#plot(pres[which(pres$method=="BB"),]$sample_size,pres[which(pres$method=="BB"),]$val1,type='l', xlab = "Sample size", ylab="EVPIv")
#lines(pres[which(pres$method=="OB"),]$sample_size,pres[which(pres$method=="OB"),]$val1,col='green')
#lines(pres[which(pres$method=="asy"),]$sample_size,pres[which(pres$method=="asy"),]$val1,col='orange')

#plot(pres[which(pres$method=="BB"),]$sample_size,pres[which(pres$method=="BB"),]$val2,type='l', xlab = "Sample size", ylab="EVPIv")
#lines(pres[which(pres$method=="OB"),]$sample_size,pres[which(pres$method=="OB"),]$val2,col='green')
#lines(pres[which(pres$method=="asy"),]$sample_size,pres[which(pres$method=="asy"),]$val2,col='orange')

#plot(pres[which(pres$method=="BB"),]$sample_size,pres[which(pres$method=="BB"),]$val3,type='l', xlab = "Sample size", ylab="EVPIv")
#lines(pres[which(pres$method=="OB"),]$sample_size,pres[which(pres$method=="OB"),]$val3,col='green')
#lines(pres[which(pres$method=="asy"),]$sample_size,pres[which(pres$method=="asy"),]$val3,col='orange')


#plot(pres[which(pres$method=="BB"),]$sample_size,pres[which(pres$method=="BB"),]$val4,type='l', xlab = "Sample size", ylab="EVPIv")
#lines(pres[which(pres$method=="OB"),]$sample_size,pres[which(pres$method=="OB"),]$val4,col='green')
#lines(pres[which(pres$method=="asy"),]$sample_size,pres[which(pres$method=="asy"),]$val4,col='orange')
