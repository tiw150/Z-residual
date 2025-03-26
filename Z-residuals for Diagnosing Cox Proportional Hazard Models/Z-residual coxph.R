####### coxph package  ####################################################
resid_coxph<-function(coxfit_fit)
{
  y<- coxfit_fit$y
  m <- nrow (y)
  mre <- resid(coxfit_fit, type="martingale")
  dre <- resid(coxfit_fit, type="deviance")
  #Unmodified Cox-Snell residual
  ucs <- as.data.frame(as.matrix(y))[,-1] - mre
  #Survival Function
  SP<- exp(-ucs)
  censored <- which(as.data.frame(as.matrix(y))[,-1]==0)
  n.censored <- length(censored)
  #NRSP residual
  RSP <- SP
  RSP[censored] <- RSP[censored]*runif(n.censored)
  nrsp <- qnorm(RSP)
  #Normalized unmodified SPs (NUSP)
  USP<-SP
  USP[USP==1] <- .999999999
  nusp<- -qnorm(USP)
  
  list(nrsp=nrsp,nusp=nusp)
}
test.nl.aov <- function(Zresidual, fitted.value, k.anova=10)
{
  if(is.factor(fitted.value)){
    fitted.value<-as.numeric(fitted.value)-1
    lpred.bin <- fitted.value
    anova(lm(Zresidual ~ lpred.bin))$`Pr(>F)`[1]
  }
  
  if(!is.factor(fitted.value)){
    lpred.bin <- cut(fitted.value, k.anova)
    less2_factor<-which(tapply(lpred.bin,lpred.bin,length)<= 2)
    if(rlang::is_empty(names(less2_factor))){
      anova(lm(Zresidual ~ lpred.bin))$`Pr(>F)`[1]
    }else{
      list_less2_factor<-list()
      for(j in 1:length(less2_factor)){
        list_less2_factor[[j]]<-which(lpred.bin==names(less2_factor[j]))
      }
      vector_less2_factor<-unlist(list_less2_factor, use.names = FALSE)
      new.lpred.bin<- lpred.bin[-vector_less2_factor]
      new.Zresidual<-Zresidual[-vector_less2_factor]
      anova(lm(new.Zresidual ~ new.lpred.bin))$`Pr(>F)`[1]
    }
    
  }
  
}

test.var.bartl <- function(Zresidual, fitted.value, k.bl=10)
{
  if(is.factor(fitted.value)){
    lpred.bin <- fitted.value
    Z_group<- split(Zresidual, lpred.bin)
    bl.test<-bartlett.test(Z_group)[["p.value"]]
  }
  if(!is.factor(fitted.value)){
    lpred.bin <- cut(fitted.value, k.bl)
    Z_group<- split(Zresidual, lpred.bin)
    check_Z_group<-rep(k.bl)
    for(i in 1:k.bl)
    {
      fun<-function(x) x>2
      check_Z_group[i]<-fun(length(Z_group[[i]]))
    }
    if(all(check_Z_group!=0)){
      bl.test<-bartlett.test(Z_group)[["p.value"]]
    }else{
      Z_group<-Z_group[-which(check_Z_group==0)]
      bl.test<-bartlett.test(Z_group)[["p.value"]]
    }
  }
  bl.test
}





bounds_pvalues <- function (pv)
{
  pv <- pv[is.finite(pv)]
  n <- length (pv)
  if(n<=0) stop("There is no value in 'pv'")
  pv.sorted <- sort(pv)
  corrected.pv <- pmin(1, pv.sorted*n/1:n)
  pv.min <- min(corrected.pv)
  pv.min
}
