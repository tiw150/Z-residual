#input: survreg_fit is a survreg object
resid_survreg<-function(survreg_fit)
{
  distr<-survreg_fit$dist
  y<- survreg_fit$y
  m <- nrow (y)
  parms<- as.numeric(survreg_fit[["parms"]])
  alpha_hat<-1/survreg_fit$scale
  
  if (distr %in% c("weibull","exponential","logistic","lognormal",
                   "loglogistic","gaussian", "loggaussian","rayleigh"))
  {
    SP<-1-(psurvreg(as.data.frame(as.matrix(y))[,-2], 
                    survreg_fit[["linear.predictors"]], 
                    scale=1/alpha_hat, distribution=distr))
    haz <- dsurvreg(as.data.frame(as.matrix(y))[,-2],
                    survreg_fit[["linear.predictors"]],
                    scale=1/alpha_hat, distribution=distr)/SP
  }else if  (distr %in% "t") 
  {
    SP<-1-(psurvreg(as.data.frame(as.matrix(y))[,-2], 
                    survreg_fit[["linear.predictors"]], 
                    scale=1/alpha_hat, distribution=distr,
                    parms = parms))
    haz <- dsurvreg(as.data.frame(as.matrix(y))[,-2], 
                    survreg_fit[["linear.predictors"]],
                    scale=1/alpha_hat, distribution=distr,
                    parms = parms)/SP
  }else stop ("The distribution is not supported")
  censored <- which(as.data.frame(as.matrix(y))[,-1]==0)
  n.censored <- length(censored)
  # NRSP
  RSP <- SP
  RSP[censored] <- RSP[censored]*runif(n.censored)
  nrsp <- qnorm(RSP)
  #Normalized unmodified SPs (NUSP)
  USP<-SP
  USP[USP==1] <- .999999999
  nusp<- -qnorm(USP)
  # Unmodified CS residual
  ucs<- -log(SP)
  # Modified CS residual
  MSP<- SP
  MSP[censored] <- SP[censored]/exp(1)
  mcs <- -log(MSP)
  nmsp <- qnorm (MSP)
  #Martingale Residual
  martg<- as.data.frame(as.matrix(survreg_fit$y))[,-1]- ucs
  #Deviance Residual
  dev<- sign(martg)* sqrt((-2)*(martg+as.data.frame(as.matrix(survreg_fit$y))[,-1]*
                                  log(as.data.frame(as.matrix(survreg_fit$y))[,-1]-martg)))
  nrsp.sw.pvalue<-shapiro.test(nrsp)$p.value
  #hazard function
  haz_fn<- haz
  
  list(RSP=RSP,nrsp=nrsp,nrsp.sw.pvalue=nrsp.sw.pvalue,nusp=nusp,
       ucs=ucs, mcs=mcs,martg=martg, dev=dev, nmsp=nmsp,haz_fn=haz_fn)
}
# outputs:
## RSP --- Randomized Survival Probabilities
## nrsp --- Normalized RSP 
## nrsp.sw.pvalue --- GOF test p-values by applying SW test to NRSP
## ucs --- unmodified CS residuals
## mcs --- modified CS residuals
## nmsp --- normalized modified SP
## martg --- Martingale residuals
## dev --- Deviance residuals
## haz_fn --- hazard function of cs residuals

