####### coxph package  ####################################################
zresidual.coxph <- function (fit_coxph, traindata, newdata)
{
  if (!requireNamespace("pacman")) {
    install.packages("pacman")
  }
  pacman::p_load(
    "stringr"
    # "data.table"
  )
  # data_name<- as.character(fit_coxph$call[[3]])
  # if(!exists(data_name)) stop(paste(data_name, "does not exist"))
  # traindata<-get(data_name)
  
  form<-(fit_coxph$formula)[[3]]
  group_id_name<-gsub(".*[(]([^.]+)[,].*", "\\1", form)[3]
  #  group_id_name<-noquote(str_replace_all(string=group_id_name, pattern=" ", repl=""))
  if(!is.factor(traindata[[group_id_name]])) stop("The group ID must be factor!")
  if(!is.factor(newdata[[group_id_name]])) stop("The group ID must be factor!")
  gpnumber<-length(levels(traindata[[group_id_name]]))
  
  mf <- model.frame(fit_coxph$formula, traindata)
  mf_nc<-ncol (mf)
  #gpnumber<-data.table::uniqueN(mf[,mf_nc])
  mm <- model.matrix(fit_coxph$formula, traindata)
  mm_nc <- ncol (mm)
  if(gpnumber>5){
    fix_var<-mm[,-c(1,mm_nc),drop=FALSE]
    explp<-exp(fix_var %*% fit_coxph$coefficients+fit_coxph$frail[mf[,mf_nc]])
    #z_hat is the estimate of random effect
    z_hat<-exp(fit_coxph$frail[mf[,mf_nc]])
  }
  if(gpnumber<=5){
    fix_var<-mm[,-c(1,mf_nc:mm_nc),drop=FALSE]
    coef_number<- ncol(fix_var)
    frailty<-fit_coxph$coefficients[coef_number+1:gpnumber]
    explp<-exp(fix_var %*% fit_coxph$coefficients[1:coef_number]+
                 frailty[mf[,mf_nc]])
    #z_hat is the estimate of random effect
    z_hat<-exp(frailty[mf[,mf_nc]])
  }
  
  Y <- mf[[1]]
  if(!inherits(Y, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y) != 3) {
    # making it all in (tstart, tstop) format
    Y <- Surv(rep(0, nrow(Y)), Y[,1], Y[,2])
  }
  
  # this one gives the baseline cumulative hazard at all the time points;
  getchz <- function(Y, explp) {
    death <- (Y[, ncol(Y)] == 1) # this is a TRUE FALSE for the status column
    dtime <- Y[, ncol(Y) - 1] # this is the tstop
    
    time <- sort(unique(dtime)) # unique tstops
    
    nevent <- as.vector(rowsum(1 * death, dtime))
    
    nrisk <- rev(cumsum(rev(rowsum(explp, dtime)))) # This gives the sum
    delta <- min(diff(time))/2
    etime <- c(sort(unique(Y[, 1])), max(Y[, 1]) + delta)  #unique entry times
    
    indx <- approx(etime, 1:length(etime), time, method = "constant", 
                   rule = 2, f = 1)$y
    
    esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))  #not yet entered
    nrisk <- nrisk - c(esum, 0)[indx]
    
    haz <- nevent/nrisk
    cumhaz <- cumsum(haz)
    
    out<-data.frame(time = time, haz = haz, cumhaz = cumhaz)
    return (out)
  }
  df<-getchz(Y=Y,explp=explp)
  t<-df$time
  haz<-df$haz
  H_0<-df$cumhaz
  f <- stepfun(t[-1], H_0)
  f_haz<-stepfun(t[-1], haz)
  
  full_data<-rbind(newdata,traindata)
  mf_new<-model.frame(fit_coxph$formula,full_data)[c(1:nrow(newdata)),]
  mf_nc_new<- ncol (mf_new)
  mm_new <-model.matrix(fit_coxph$formula,full_data)[c(1:nrow(newdata)),,drop=FALSE]
  mm_nc_new <- ncol (mm_new)
  
  if(gpnumber>5){
    fix_var_new<-mm_new[,-c(1,mm_nc_new),drop=FALSE]
    #z_hat is the estimate of random effect
    z_hat_new<-exp(fit_coxph$frail[mf_new[,mf_nc_new]])
    #explp is linear predictor with random effect
    explp_new<-exp(fix_var_new %*% fit_coxph$coefficients)
  }
  if(gpnumber<=5){
    fix_var_new<-mm_new[,-c(1,mf_nc_new:mm_nc_new),drop=FALSE]
    coef_number_new<- ncol(fix_var_new)
    frailty_new<-fit_coxph$coefficients[coef_number_new+1:gpnumber]
    explp_new<-exp(fix_var_new %*% fit_coxph$coefficients[1:coef_number_new])
    #z_hat is the estimate of random effect
    z_hat_new<-as.numeric(exp(frailty_new[mf_new[,mf_nc_new]]))
  }
  
  Y_new <- mf_new[[1]]
  if(!inherits(Y_new, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y_new) != 3) {
    # making it all in (tstart, tstop) format
    Y_new <- Surv(rep(0, nrow(Y_new)), Y_new[,1], Y_new[,2])
  }
  H0_new<-f(Y_new[,2])
  h0_new<-f_haz(Y_new[,2])
  #hazard function
  hazfn<- z_hat_new*as.vector(explp_new)*h0_new
  #Survival Function
  SP<- exp(-z_hat_new*as.vector(explp_new)*H0_new)
  censored <- which(Y_new[,3]==0)
  n.censored <- length(censored)
  #NRSP residual
  RSP <- SP
  RSP[censored] <- RSP[censored]*runif(n.censored)
  #R_H(t)
  log_Rcum_haz<- log(-log(RSP))
  Z_resid <- -qnorm(RSP)
  attr(Z_resid, "SP") <- SP
  attr(Z_resid, "log_Rcum_haz") <- log_Rcum_haz
  return(Z_resid)
}

#basehaz(fit_kidney1,centered = F)
test.nl.aov1 <- function(zresidual, fitted.values, k.anova=10)
{
  lpred.bin <- cut(fitted.values, k.anova)
  less2_factor<-which(tapply(lpred.bin,lpred.bin,length)<= 2)
  if(rlang::is_empty(names(less2_factor))){
    anova(lm(zresidual ~ lpred.bin))$`Pr(>F)`[1]
  }else{
    vector_less2_factor<-rep(length(less2_factor))
    for(j in 1:length(less2_factor)){
      vector_less2_factor[j]<-which(lpred.bin==names(less2_factor[j]))
    }
    new.lpred.bin<- lpred.bin[-vector_less2_factor]
    new.zresiduall<-zresidual[-vector_less2_factor]
    anova(lm(new.zresidual ~ new.lpred.bin))$`Pr(>F)`[1]
  }
}


test.nl.aov1.categ <- function(zresidual, fitted.values)
{
  lpred.bin <- fitted.values
  anova(lm(zresidual ~ lpred.bin))$`Pr(>F)`[1]
}


test.var.bartl <- function(zresidual, fitted.values, k=10)
{
  lpred.bin <- cut(fitted.values, k)
  Z_group<- split(zresidual, lpred.bin)
  check_Z_group<-rep(k)
  for(i in 1:k)
  {
    fun<-function(x) x>2
    check_Z_group[i]<-fun(length(Z_group[[i]]))
  }
  if(all(check_Z_group!=0))
  {
    bartlett.test(Z_group)[["p.value"]]
  }else{
    Z_group<-Z_group[-which(check_Z_group==0)]
    bartlett.test(Z_group)[["p.value"]]
  }
}

test.var.bartl.categ <- function(zresidual, fitted.values)
{
  lpred.bin <- fitted.values
  Z_group<- split(zresidual, lpred.bin)
  bartlett.test(Z_group)[["p.value"]]
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

allresidual.coxph <- function (fit_coxph, traindata, newdata)
{
  if (!requireNamespace("pacman")) {
    install.packages("pacman")
  }
  pacman::p_load(
    "stringr"
    # "data.table"
  )
  # data_name<- as.character(fit_coxph$call[[3]])
  # if(!exists(data_name)) stop(paste(data_name, "does not exist"))
  # traindata<-get(data_name)
  
  form<-(fit_coxph$formula)[[3]]
  group_id_name<-gsub(".*[(]([^.]+)[,].*", "\\1", form)[3]
  #  group_id_name<-noquote(str_replace_all(string=group_id_name, pattern=" ", repl=""))
  if(!is.factor(traindata[[group_id_name]])) stop("The group ID must be factor!")
  if(!is.factor(newdata[[group_id_name]])) stop("The group ID must be factor!")
  gpnumber<-length(levels(traindata[[group_id_name]]))
  
  mf <- model.frame(fit_coxph$formula, traindata)
  mf_nc<-ncol (mf)
  #gpnumber<-data.table::uniqueN(mf[,mf_nc])
  mm <- model.matrix(fit_coxph$formula, traindata)
  mm_nc <- ncol (mm)
  if(gpnumber>5){
    fix_var<-mm[,-c(1,mm_nc),drop=FALSE]
    explp<-exp(fix_var %*% fit_coxph$coefficients+fit_coxph$frail[mf[,mf_nc]])
    #z_hat is the estimate of random effect
    z_hat<-exp(fit_coxph$frail[mf[,mf_nc]])
  }
  if(gpnumber<=5){
    fix_var<-mm[,-c(1,mf_nc:mm_nc),drop=FALSE]
    coef_number<- ncol(fix_var)
    frailty<-fit_coxph$coefficients[coef_number+1:gpnumber]
    explp<-exp(fix_var %*% fit_coxph$coefficients[1:coef_number]+
                 frailty[mf[,mf_nc]])
    #z_hat is the estimate of random effect
    z_hat<-exp(frailty[mf[,mf_nc]])
  }
  
  Y <- mf[[1]]
  if(!inherits(Y, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y) != 3) {
    # making it all in (tstart, tstop) format
    Y <- Surv(rep(0, nrow(Y)), Y[,1], Y[,2])
  }
  
  # this one gives the baseline cumulative hazard at all the time points;
  getchz <- function(Y, explp) {
    death <- (Y[, ncol(Y)] == 1) # this is a TRUE FALSE for the status column
    dtime <- Y[, ncol(Y) - 1] # this is the tstop
    
    time <- sort(unique(dtime)) # unique tstops
    
    nevent <- as.vector(rowsum(1 * death, dtime))
    
    nrisk <- rev(cumsum(rev(rowsum(explp, dtime)))) # This gives the sum
    delta <- min(diff(time))/2
    etime <- c(sort(unique(Y[, 1])), max(Y[, 1]) + delta)  #unique entry times
    
    indx <- approx(etime, 1:length(etime), time, method = "constant", 
                   rule = 2, f = 1)$y
    
    esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))  #not yet entered
    nrisk <- nrisk - c(esum, 0)[indx]
    
    haz <- nevent/nrisk
    cumhaz <- cumsum(haz)
    
    out<-data.frame(time = time, haz = haz, cumhaz = cumhaz)
    return (out)
  }
  df<-getchz(Y=Y,explp=explp)
  t<-df$time
  haz<-df$haz
  H_0<-df$cumhaz
  f <- stepfun(t[-1], H_0)
  
  full_data<-rbind(newdata,traindata)
  mf_new<-model.frame(fit_coxph$formula,full_data)[c(1:nrow(newdata)),]
  mf_nc_new<- ncol (mf_new)
  mm_new <-model.matrix(fit_coxph$formula,full_data)[c(1:nrow(newdata)),,drop=FALSE]
  mm_nc_new <- ncol (mm_new)
  
  if(gpnumber>5){
    fix_var_new<-mm_new[,-c(1,mm_nc_new),drop=FALSE]
    #z_hat is the estimate of random effect
    z_hat_new<-exp(fit_coxph$frail[mf_new[,mf_nc_new]])
    #explp is linear predictor with random effect
    explp_new<-exp(fix_var_new %*% fit_coxph$coefficients)
  }
  if(gpnumber<=5){
    fix_var_new<-mm_new[,-c(1,mf_nc_new:mm_nc_new),drop=FALSE]
    coef_number_new<- ncol(fix_var_new)
    frailty_new<-fit_coxph$coefficients[coef_number_new+1:gpnumber]
    explp_new<-exp(fix_var_new %*% fit_coxph$coefficients[1:coef_number_new])
    #z_hat is the estimate of random effect
    z_hat_new<-as.numeric(exp(frailty_new[mf_new[,mf_nc_new]]))
  }
  
  Y_new <- mf_new[[1]]
  if(!inherits(Y_new, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y_new) != 3) {
    # making it all in (tstart, tstop) format
    Y_new <- Surv(rep(0, nrow(Y_new)), Y_new[,1], Y_new[,2])
  }
  H0_new<-f(Y_new[,2])
  #Survival Function
  SP<- exp(-z_hat_new*as.vector(explp_new)*H0_new)
  censored <- which(Y_new[,3]==0)
  n.censored <- length(censored)
  
  #NRSP residual
  RSP <- SP
  RSP[censored] <- RSP[censored]*runif(n.censored)
  nrsp <- -qnorm(RSP)
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
  #Normalized MSPs
  MSP<- SP
  MSP[censored] <- SP[censored]/exp(1)
  nmsp<- -qnorm(MSP)
  #Martingale Residual
  martg<- Y_new[,3] - ucs
  #Deviance Residual
  dev<- sign(martg)* sqrt((-2)*(martg+Y_new[,3]*log(Y_new[,3]-martg)))
  list(nrsp=nrsp,nusp=nusp,ucs=ucs,mcs=mcs,nmsp=nmsp,martg=martg,dev=dev)
}
