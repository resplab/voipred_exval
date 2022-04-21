library(mvtnorm)
library(progress)
library(sqldf)

if(exists("instanceId")) #Means GRcomp is loaded!
{
  source(GRdownload("core.R"))
}else
{
  source('./R/core.R')
}


settings <- list(
  true_model = c(-2,1),  #as.matrix(c(-2,1,0,1,0)),
  z = 0.1,
  sample_sizes= c(500,1000,2000,4000,8000),
  model_generator_sigma = diag(2)/4,
  model_generator_sample_size=50,
  universe_n = 10^7,
  n_sim_outer= 1000,
  max_n_sim_inner= 1000,
  verbose = F
)

universe <- new.env()


#Creates a model based on a sample of sample_size. Puts it in the universe (also calculates pi andNB_model)
generate_proposed_model <- function(sample_size=NULL, sigma=NULL)
{
  if(!is.null(sigma))
  {
    betas <- rmvnorm(1, settings$true_model, sigma)
    out <- as.vector(betas)
  }
  else if(!is.null(sample_size))
  {
    x <- generate_sample(sample_size, settings$true_model)
    p <- 1/(1+exp(-x%*%settings$true_model))
    y <- rbinom(sample_size,1,p)
    out <- as.vector(coefficients(glm(y ~ 0 + x), family=binomial(link='logit')))
  } else {stop("Neither sigma nor sample size provided. I dont know how to generate a model.")}

  out
}


generate_sample <- function(n, model)
{
  x <- rmvnorm(n, model*0,diag(length(model)))
  x[,1] <- 1
  #x[,4] <- (x[,4]>0)*1
  #x[,5] <- (x[,5]>0)*1

  x
}


#Creates a universe with covariates, true risks (from the true model), and will also contant predicted risks from a proposed model
generate_universe <- function(n, true_model)
{
  z <- settings$z
  universe <- new.env()

  x <- generate_sample(n, model=true_model)
  universe$x <- x
  p <- as.vector(1/(1+exp(-x%*%as.matrix(true_model))))
  universe$p <- p
  message("Average p in the universe is ",mean(p))
  universe$pi <- NULL #For now
  universe$NB_all <- mean((p-(1-p)*z/(1-z)))
  universe$NB_max <- mean((p>z)*(p-(1-p)*z/(1-z)))
  universe$NB_model <- NULL #For now
  universe$true_model <- true_model

  universe <<- universe
  universe
}



meta_sim <- function(n_sim_outer, max_n_sim_inner, sample_size, universe, method)
{
  #cat(method,'\n')
  z <- settings$z

  evpi_v <- winner_val <- rep(NA,n_sim_outer)
  for(i_sim in 1:n_sim_outer)
  {
    #plot(pi[1:1000],p[1:1000])
  
    repeat
    {
      x_ <- generate_sample(sample_size,universe$true_model)
      xx_ <- x_[,2]^2 #For fun
  
      p_ <- as.vector(1/(1+exp(-x_%*%as.matrix(universe$true_model))))
      y_ <- rbinom(sample_size, size = 1, prob=p_)
      pi_ <- as.vector(1/(1+exp(-x_%*%as.matrix(universe$proposed_model))))
      if(min(pi_)<z && max(pi_)>z) break;
    }
    
    bs_NB_model <- bs_NB_all <- bs_max_NB <- 0
    s2s <- c(0,0,0)
    j_sim <- 0

    if(method=="asymptotic")
    {
      parms <- NB_BVN(y_, pi_, z)
      #cat(parms,'\n')
      if(is.nan(parms[5])) {browser()}
      if(parms[5]>0.99999) 
      {
        #rbrowser();
        parms[5]<-0.999; 
      }
      tryCatch(
      {
        first_term <- do.call(mu_max_truncated_bvn,as.list(parms))
      }, 
      error=function(e) 
      {
        message("Error");  browser(); assign("first_term",max(0,parms[1],parms[2]),envir = parent.env(environment()))
      })
      
      evpi_v[i_sim] <- first_term-max(0,parms[1],parms[2])
      winner_val[i_sim] <- which.max(c(0,parms[1],parms[2]))
    }
    else
    {
      if(method=="model_based_ll")
      {
        local_model <- glm(formula = y_ ~ 0 + x_, family=binomial(link="logit"))
        mus <- as.vector(coefficients(local_model))
        covmat <- vcov(local_model)
      }
      
      repeat
      {
        j_sim <- j_sim+1
        w_x <- bootstrap(sample_size, Bayesian = T)
        w_y <- bootstrap(sample_size, Bayesian = T)
  
        if(method=="model_based_bs")
        {
          local_model <-suppressWarnings(glm(formula = y_ ~ 0 + x_, family=binomial(link="logit"), weights = w_y))
          p__ <- predict(local_model, type="response", newdata = as.data.frame(x_))
        }
        else if(method=="model_based_ll")
        {
          betas <- rmvnorm(1,mus,covmat)
          p__ <- as.vector(1/(1+exp(-x_%*%t(betas))))
        }
        else if(method=="BB")
        {
          p__ <- y_  #naive method with BB
        }
        else if(method=="OB")
        {
          w_y <- bootstrap(sample_size, Bayesian = F)
          p__ <- y_  #naive method with BB
        }
        else
        {
            stop("Method",method,"Not regonized")
        }
        
        bs_NB_model_this <- sum(w_x*(pi_>z)*(p__ - (1-p__)*z/(1-z)))/sample_size
        bs_NB_all_this <- sum(w_x*(p__ - (1-p__)*z/(1-z)))/sample_size
        bs_max_NB <- bs_max_NB + max(0, bs_NB_model_this, bs_NB_all_this)
  
        #tmp_model <- sum(1*(pi_>z)*(p__ - (1-p__)*z/(1-z)))/sample_size
        #tmp_all <- sum(1*(p__ - (1-p__)*z/(1-z)))/sample_size
  
        s2s <- s2s + c(bs_NB_model_this^2,bs_NB_all_this^2,(bs_NB_model_this-bs_NB_all_this)^2)
  
        bs_NB_model <- bs_NB_model + bs_NB_model_this
        bs_NB_all <- bs_NB_all + bs_NB_all_this
      
        if(j_sim%%100==0)
        {
          ses <- sqrt((s2s/j_sim-c((bs_NB_model/j_sim)^2,(bs_NB_all/j_sim)^2,((bs_NB_model-bs_NB_all)/j_sim)^2))/j_sim)
          cvs <- abs(ses/c(bs_NB_model/j_sim,bs_NB_all/j_sim,bs_NB_model/j_sim-bs_NB_all/j_sim))
          #if(sum(is.nan(cvs))>0) browser();
          if(max(cvs[which(ses>0)])<0.2 || j_sim==max_n_sim_inner) break;
        }
      } #repeat
      evpi_v[i_sim] <- bs_max_NB/j_sim - max(0,bs_NB_model,bs_NB_all)/j_sim
      winner_val[i_sim] <- which.max(c(0,bs_NB_model,bs_NB_all))
    }
  }

  if(settings$verbose) print(winner_val)

  c(
    evpi_v=mean(evpi_v),
    meta_evpi_v = max(0,universe$NB_model,universe$NB_all)-mean(c(0,universe$NB_model,universe$NB_all)[winner_val])
  )
}





main <- function()
{
  if(!exists("instanceId"))
  {
    instanceId <- 0
  }
  set.seed(instanceId)
    
  message("Instance id:",instanceId)
  
  generate_universe(n=settings$universe_n, true_model=settings$true_model)

  sample_sizes <- settings$sample_sizes

  pb <- progress_bar$new(total=settings$n_sim_outer)

  res <- data.frame()

  for(i in 1:settings$n_sim_outer)
  {
    pb$tick()

    repeat
    {
      proposed_model <- generate_proposed_model(sigma = settings$model_generator_sigma)
      universe$proposed_model <<- proposed_model
      universe$pi <<- as.vector(1/(1+exp(-universe$x%*%as.matrix(proposed_model))))
      qs <- quantile(universe$pi,c(0.01,0.99))
      if(qs[2]>settings$z && qs[1]<settings$z) 
      {
        break;
      }else
      {
        message("bad")
      }
    }
    universe$NB_model <<- mean((universe$pi>settings$z)*(universe$p-(1-universe$p)*settings$z/(1-settings$z)))
    print(proposed_model)
    
    for(sample_size in sample_sizes)
    {
      #  res <- rbind(res,c(method="model_based_bs",sample_size=sample_size,
      #                     as.list(meta_sim(n_sim_outer=1, max_n_sim_inner=settings$max_n_sim_inner, sample_size=sample_size, universe=universe, method="model_based_bs"))))
      #  res <- rbind(res,c(method="model_based_ll",sample_size=sample_size,
      #                     as.list(meta_sim(n_sim_outer=1, max_n_sim_inner=settings$max_n_sim_inner, sample_size=sample_size, universe=universe, method="model_based_ll"))))
        res <- rbind(res,c(method="BB",sample_size=sample_size,
                           as.list(meta_sim(n_sim_outer=1, max_n_sim_inner=settings$max_n_sim_inner, sample_size=sample_size, universe=universe, method="BB"))))
        res <- rbind(res,c(method="OB",sample_size=sample_size,
                           as.list(meta_sim(n_sim_outer=1, max_n_sim_inner=settings$max_n_sim_inner, sample_size=sample_size, universe=universe, method="OB"))))
        res <- rbind(res,c(method="asymptotic",sample_size=sample_size,
                         as.list(meta_sim(n_sim_outer=1, max_n_sim_inner=settings$max_n_sim_inner, sample_size=sample_size, universe=universe, method="asymptotic"))))
      
    }

    if(i%%50==0 && instanceId>0)
    {
      GRpush(cbind(instanceId=instanceId,res),TRUE)
    }
  }

  res <- as.data.frame(res)
  res <<- res
  res
}



process_results <- function()
{
  require("sqldf")
  #bsv <- sqldf("SELECT sample_size, COUNT(*) AS N, AVG(evpi_v) AS val, SQRT(VARIANCE(evpi_v)/COUNT(*)) AS val_se, AVG(meta_evpi_v) AS meta, SQRT(VARIANCE(meta_evpi_v)/COUNT(*)) AS meta_se,  STDEV(evpi_v)/SQRT(COUNT(*)) aS SE FROM res WHERE method='model_based_bs' GROUP by sample_size ORDER BY sample_size*1")
  #llv <- sqldf("SELECT sample_size, COUNT(*) AS N, AVG(evpi_v) AS val, SQRT(VARIANCE(evpi_v)/COUNT(*)) AS val_se, AVG(meta_evpi_v) AS meta, SQRT(VARIANCE(meta_evpi_v)/COUNT(*)) AS meta_se,  STDEV(evpi_v)/SQRT(COUNT(*)) aS SE FROM res WHERE method='model_based_ll' GROUP by sample_size ORDER BY sample_size*1")
  BB <- sqldf("SELECT sample_size, COUNT(*) AS N, AVG(evpi_v) AS val, SQRT(VARIANCE(evpi_v)/COUNT(*)) AS val_se, AVG(meta_evpi_v) AS meta, SQRT(VARIANCE(meta_evpi_v)/COUNT(*)) AS meta_se,  STDEV(evpi_v)/SQRT(COUNT(*)) aS SE FROM res WHERE method='BB' GROUP by sample_size ORDER BY sample_size*1")
  OB <- sqldf("SELECT sample_size, COUNT(*) AS N, AVG(evpi_v) AS val, SQRT(VARIANCE(evpi_v)/COUNT(*)) AS val_se, AVG(meta_evpi_v) AS meta, SQRT(VARIANCE(meta_evpi_v)/COUNT(*)) AS meta_se,  STDEV(evpi_v)/SQRT(COUNT(*)) aS SE FROM res WHERE method='OB' GROUP by sample_size ORDER BY sample_size*1")
  asy<- sqldf("SELECT sample_size, COUNT(*) AS N, AVG(evpi_v) AS val, SQRT(VARIANCE(evpi_v)/COUNT(*)) AS val_se, AVG(meta_evpi_v) AS meta, SQRT(VARIANCE(meta_evpi_v)/COUNT(*)) AS meta_se,  STDEV(evpi_v)/SQRT(COUNT(*)) aS SE FROM res WHERE method='asymptotic' GROUP by sample_size ORDER BY sample_size*1")
  
  #plot(bsv$sample_size, bsv$meta, type='l', ylim=c(0,max(bsv$val)), col='red', main="bsv")
  #lines(bsv$sample_size, bsv$val, type='l',col='blue')

  # plot(llv$sample_size, llv$meta, type='l', ylim=c(0,max(llv$val)), col='red', main="llv")
  # lines(llv$sample_size, llv$val, type='l',col='blue')
  # 
  plot(BB$sample_size, BB$meta, type='l',col='red', ylim=c(0,max(BB$val)), main="Bayesian bootstrap", xlab="Sample size", ylab="EVPIv")
  lines(BB$sample_size, BB$val, type='l',col='blue')
   
  plot(OB$sample_size, OB$meta, type='l',col='red', ylim=c(0,max(OB$val)), main="Ordinary bootstrap", xlab="Sample size", ylab="EVPIv")
  lines(OB$sample_size, OB$val, type='l',col='blue')
  
  plot(asy$sample_size, asy$meta, type='l',col='red', ylim=c(0,max(asy$val)), main="Asymptotic", xlab="Sample size", ylab="EVPIv")
  lines(asy$sample_size, asy$val, type='l',col='blue')
  

  # plot(bsd$sample_size, bsd$meta, type='l', ylim=c(0,max(bsd$val)), col='red', main="bsd")
  # lines(bsd$sample_size, bsd$val, type='l',col='blue')
  # 
  # plot(lld$sample_size, lld$meta, type='l', ylim=c(0,max(lld$val)), col='red', main="lld")
  # lines(lld$sample_size, lld$val, type='l',col='blue')
  # 
  # plot(nvd$sample_size, nvd$meta, type='l',col='red', ylim=c(0,max(nvd$val)), main="nvd")
  # lines(nvd$sample_size, nvd$val, type='l',col='blue')

  print(
    sqldf("SELECT method, sample_size, AVG(evpi_v-meta_evpi_v) AS bias, AVG(POWER(evpi_v-meta_evpi_v,2)) AS rmse FROM res GROUP BY method, sample_size")
  )

  return(list(
    #bsd = bsd,
    #model_based_bbs = bsv,
    #lld = lld,
    #model_based_ll = llv,
    #nvd = nvd,
    BB = BB,
    OB = OB,
    asy = asy
  ))
}


main()

