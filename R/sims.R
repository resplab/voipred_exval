library(voipred)
library(rms)

data(gusto)
gusto$kill <- (as.numeric(gusto$Killip)>1)*1
gusto$Y <- gusto$day30
data_us <- gusto[gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]
data_other <- gusto[!gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]


settings <- list(
  model_formula = Y ~ age + miloc + pmi + kill + pmin(sysbp,100) + lsp(pulse,50),
  z = 0.05,
  zs = (0:99)/100,
  Bayesian_bootstrap = T
)




voi_ex <- function(n_sim, n_dev, n_val, zs=settings$zs)
{
  dev_data <- gusto[sample(1:(dim(gusto)[1]),n_dev,F),]
  val_data <- gusto[sample(1:(dim(gusto)[1]),n_val,F),]

  model <- glm(formula = settings$model_formula, data=dev_data, family=binomial(link="logit"))

  pi <- predict(model, type="response", newdata=val_data)
  val_data$pi <- pi

  bs_NB_model <- bs_NB_all <- bs_max_NB <- bs_NB_max <- matrix(0, nrow=n_sim, ncol=length(zs))

  for(i_sim in 1:n_sim)
  {
    w_x <- voipred:::bootstrap(n,settings$Bayesian_bootstrap)
    val_data$w_pop <- voipred:::bootstrap(n_val,settings$Bayesian_bootstrap)

    pop_model <- glm(formula = settings$model_formula, family=binomial(link="logit"), data=val_data, weights = w_pop)

    p <- predict(pop_model, type="response", newdata = val_data)

    for(i in 1:length(z))
    {
      bs_NB_model[i_sim,i] <- sum(w_x*(val_data$pi>z[i])*(p - (1-p)*z[i]/(1-z[i])))/n
      bs_NB_all[i_sim,i] <- sum(w_x*(p - (1-p)*z[i]/(1-z[i])))/n
      bs_max_NB[i_sim,i] <- max(bs_NB_model[i_sim,i],bs_NB_all[i_sim,i],0)
      bs_NB_max[i_sim,i] <- sum(w_x*((p>z[i])*(p - (1-p)*z[i]/(1-z[i]))))/n
    }
  }

  ENB_current <- pmax(colMeans(bs_NB_all),colMeans(bs_NB_model),0)
  ENB_perfect_val <- colMeans(bs_max_NB)
  ENB_perfect_dev <- colMeans(bs_NB_max)
  EVPI_val <- ENB_perfect_val-ENB_current
  EVPI_dev <- ENB_perfect_dev-ENB_current

  # EINB_current <- (ENB_current-pmax(colMeans(bs_NB_all),0))
  # EINB_perfect <- (ENB_perfect_dev-pmax(colMeans(bs_NB_all),0))
  # EVPIr_dev <- EINB_perfect/EINB_current

  # plot(z,EINB_perfect,col='red', type='l')
  # lines(z,EINB_current)
  #
  # plot(z, EVPIr_dev, type='l', ylim=c(0,min(max(EVPIr_dev,na.rm=T),10)))
  #
  # max_y <- max(EVPI_val,EVPI_dev)
  # plot(z,EVPI_val,type='l', col='red', ylim=c(0,max_y))
  # lines(z,EVPI_dev,type='l', col='green')

  return(list(EVPI_val=EVPI_val,EVPI_dev=EVPI_dev))
}



main <- function(n_sim, n_dev, n_val)
{
  EVPI_val <- EVPI_dev <- rep(0,length(settings$zs))
  for(i in 1:n_sim)
  {
    cat('.')
    res <- voi_ex(200,n_dev,n_val)
    EVPI_val <- EVPI_val + res$EVPI_val
    EVPI_dev <- EVPI_dev + res$EVPI_dev
  }

  EVPI_val=EVPI_val/n_sim
  EVPI_dev=EVPI_dev/n_sim

  max_y <- max(EVPI_val,EVPI_dev)
  plot(z,EVPI_val,type='l', col='red', ylim=c(0,max_y))
  lines(z,EVPI_dev,type='l', col='green')

  return(list(EVPI_val, EVPI_dev))
}

res_5_5 <- main(100,500,500)
res_5_10 <- main(100,500,1000)
res_10_5 <- main(100,1000,500)
res_10_10 <- main(100,1000,1000)
