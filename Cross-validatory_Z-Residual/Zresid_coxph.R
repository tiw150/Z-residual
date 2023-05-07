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
  #Z-residual
  RSP <- SP
  RSP[censored] <- RSP[censored]*runif(n.censored)
  #R_H(t)
  Rcum_haz<- -log(RSP)
  log_Rcum_haz<- log(Rcum_haz)
  #RSP<- pmax(pmin(RSP, 1-10^(-6)),10^(-6))
  Z_resid <- -qnorm(RSP)
  attr(Z_resid, "SP") <- SP
  #  attr(nrsp, "z_hat_new") <- z_hat_new
  return(Z_resid)
}

