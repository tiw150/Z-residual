library("survival")
library("EnvStats")
library("nortest")
library("dplyr")
library("spBayesSurv")
library("dvmisc")
library("frailtySurv")
### Z-residual function can calculate 
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
  if(!is.factor(traindata[[group_id_name]])) stop("The group ID must be factor!")
  if(!is.factor(newdata[[group_id_name]])) stop("The group ID must be factor!")
  gpnumber<-length(levels(traindata[[group_id_name]]))
  
  mf <- model.frame(fit_coxph$formula, traindata)
  mf_nc<-ncol (mf)
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
  
  #Z-residual
  RSP <- SP
  RSP[censored] <- RSP[censored]*runif(n.censored)
  zresid <- -qnorm(RSP)
  #Normalized unmodified SPs (censored Z-residuals)
  USP<-SP
  USP[USP==1] <- .999999999
  czresid<- -qnorm(USP)
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
  list(zresid=zresid,czresid=czresid,ucs=ucs,mcs=mcs,nmsp=nmsp,martg=martg,dev=dev)
}

test.nl.aov <- function(qresidual, fitted.values, k.anova=10)
{
  lpred.bin <- cut(fitted.values, k.anova)
  less2_factor<-which(tapply(lpred.bin,lpred.bin,length)<= 2)
  if(rlang::is_empty(names(less2_factor))){
    anova(lm(qresidual ~ lpred.bin))$`Pr(>F)`[1]
  }else{
    vector_less2_factor<-rep(length(less2_factor))
    for(j in 1:length(less2_factor)){
      vector_less2_factor[j]<-which(lpred.bin==names(less2_factor[j]))
    }
    new.lpred.bin<- lpred.bin[-vector_less2_factor]
    new.qresidual<-qresidual[-vector_less2_factor]
    anova(lm(new.qresidual ~ new.lpred.bin))$`Pr(>F)`[1]
  }
}


### Simulate non-linear data ######
set.seed(1)
rexp2 <- function(n, rate){ if (rate==0) rep(Inf,n) else rexp(n = n, rate = rate)}
simulWeib <- function(n_clusters,n_individuals, lambda, alpha,
                      beta1,beta2,beta3, mean.censor,fv)
{
  # Make covariates 
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  grpid<- paste0("g",grpid)
  x1 = runif(n_clusters*n_individuals,0,1)
  x2 = pmax(abs(round(rnorm(n_clusters*n_individuals, 0,1),2)),0.01)
  x3 = rbinom(n_clusters*n_individuals, size = 1, p = 0.5)
  
  # Make frailty (shape = a and scale = s,E(X) = a*s and Var(X) = a*s^2)
  #frvec<-rep(rgamma(n_clusters,shape=0.5, rate=1), each = n_individuals)
  frvec<-rep(rgamma(n_clusters,shape = 1/fv,scale = fv), each = n_individuals)
  
  #Event times:baseline hazard is Weibull，lambda=scale parameter, alpha=shape parameter
  n<-n_individuals * n_clusters
  v <- runif(n)
  surv_time <-
    (- log(v)/(lambda*exp(x1*beta1+log(x2)*beta2+x3*beta3)*frvec))^(1 /alpha)
  
  # Censoring times
  t0<- rexp2(n, rate= 1/mean.censor)
  # follow-up times and event indicators
  t <- pmin(surv_time, t0)
  d <- as.numeric(surv_time <= t0)
  # Make data frame
  grpid<- as.factor(grpid)
  x3<-as.factor(x3)
  out<-data.frame(id,grpid,x1,x2,x3,t,d,frvec)
  return (out) 
}

para.frame <- data.frame(
  n_clusters = c(20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,
                 20,20,20,20,20,20,20,20,20,20,20,20,20,20,20),
  n_individuals = c(10,10,10,20,20,20,30,30,30,40,40,40,50,50,50,
                    60,60,60,70,70,70, 80,80,80,90,90,90,100,100,100),
  beta1 = rep(1,times=30),
  beta2 = rep(-2,times=30),
  beta3 = rep(0.5,times=30),
  fv =  rep(0.5,times=30),
  lambda=rep(0.007,times=30),
  alpha=rep(3,times=30),
  mean.censor=rep(c(17.8,5.5,2),times=10),
  c=rep(c(20,50,80),times=10)
)

i<-1 #### i is equal to 1 to 30
para<-para.frame[i,]
n_clusters<-para[1]
n_individuals<-para[2]
beta1<-para[3]
beta2<-para[4]
beta3<-para[5]
fv<-para[6]
lambda<-para[7]
alpha<-para[8]
mean.censor<-para[9]
fit_coxph<-NA
fit_coxph_w<-NA
while(any(is.na(fit_coxph)) || any(is.na(fit_coxph_w))) {
  check_table<-FALSE
  while(isFALSE(check_table)){
    simulated_shared_frailty_data<-
      simulWeib(n_clusters=n_clusters[1,],n_individuals=n_individuals[1,],
                lambda=lambda[1,],alpha=alpha[1,],beta1=beta1[1,],
                beta2=beta2[1,],beta3=beta3[1,], 
                mean.censor= mean.censor[1,],fv=fv[1,])
    check_table<-all(apply(table(simulated_shared_frailty_data$x3,
                                 simulated_shared_frailty_data$d),
                           2,function(x) x>1))
  }
  ####fit true non-linear model############################################
  fit_coxph <- tryCatch(
    coxph(Surv(t, d) ~ x1+log(x2)+x3 +frailty(grpid,distribution = "gamma"),
          data = simulated_shared_frailty_data),
    error = function(e) NA,
    warning = function(w) NA
  )
  ####fit wrong linear model################################################
  fit_coxph_w <- tryCatch(
    coxph(Surv(t, d) ~ x1 +x2 + x3 +frailty(grpid, distribution="gamma"),
          data= simulated_shared_frailty_data),
    error = function(e) NA,
    warning = function(w) NA ) 
}
censorship<-as.numeric(table(simulated_shared_frailty_data$d)[1]/
                         nrow(simulated_shared_frailty_data))

resids_t<-zresidual.coxph (fit_coxph = fit_coxph,
                           traindata = simulated_shared_frailty_data,
                           newdata = simulated_shared_frailty_data)
resids_w<-zresidual.coxph (fit_coxph = fit_coxph_w,
                           traindata = simulated_shared_frailty_data,
                           newdata = simulated_shared_frailty_data)

sw_nmsp_t<- shapiro.test(resids_t$nmsp)$p.value
sw_nmsp_w<- shapiro.test(resids_w$nmsp)$p.value
sw_dev_t<- shapiro.test(resids_t$dev)$p.value
sw_dev_w<- shapiro.test(resids_w$dev)$p.value

ks_zresid_t<-ks.test(resids_t$zresid,"pnorm")$p.value
ks_zresid_w<-ks.test(resids_w$zresid,"pnorm")$p.value

sf_zresid_t<-sf.test(resids_t$zresid)$p.value
sf_zresid_w<-sf.test(resids_w$zresid)$p.value
sw_zresid_t<- shapiro.test(resids_t$zresid)$p.value;sw_zresid_t
sw_zresid_w<- shapiro.test(resids_w$zresid)$p.value;sw_zresid_w

anov_zresid_t<- test.nl.aov(qresidual=resids_t$zresid, 
                          fitted.values=fit_coxph$linear.predictors,
                          k.anova=10);anov_zresid_t
anov_zresid_w<- test.nl.aov(qresidual=resids_w$zresid,
                          fitted.values=fit_coxph_w$linear.predictors,
                          k.anova=10);anov_zresid_w

anov_zresid_x2_t<- test.nl.aov(qresidual=resids_t$zresid, 
                             fitted.values=(simulated_shared_frailty_data$x2),
                             k.anova=10);anov_zresid_x2_t
anov_zresid_x2_w<- test.nl.aov(qresidual=resids_w$zresid,
                             fitted.values=(simulated_shared_frailty_data$x2),
                             k.anova=10);anov_zresid_x2_w

censored<- simulated_shared_frailty_data$d ==0
gof_censored_t<-gofTestCensored(resids_t$czresid,censored, test = "sf", 
                                censoring.side = "right",
                                distribution = "norm")$p.value
gof_censored_t
gof_censored_w<-gofTestCensored(resids_w$czresid,censored, test = "sf", 
                                censoring.side = "right",
                                distribution = "norm")$p.value
gof_censored_w


#the cumulative hazard function estimated by Kaplan-Meier method of CS residuals
km.ucs.t <- survfit(Surv(resids_t$ucs,simulated_shared_frailty_data$d)~1,type='fleming')
id.ucs.t<-order(resids_t$ucs)
km.ucs.w <- survfit(Surv(resids_w$ucs,simulated_shared_frailty_data$d)~1,type='fleming')
id.ucs.w<-order(resids_w$ucs)

####### real data analysis#####################################################
data("LeukSurv")
LeukSurv<-LeukSurv[LeukSurv$age<60,]
is.factor(LeukSurv$district)
is.factor(LeukSurv$sex)
LeukSurv$district<-as.factor(LeukSurv$district)
LeukSurv$sex<-as.factor(LeukSurv$sex)
LeukSurv$wbc_log<- log(LeukSurv$wbc+0.001)

fit_LeukSurv_logwbc  <- tryCatch(
  coxph(Surv(time, cens) ~ age +sex + wbc_log + tpi +
          frailty(district, distribution="gamma"), data= LeukSurv),
  error = function(e) NA,
  warning = function(w) NA
)

AIC(fit_LeukSurv_logwbc)
resid_LeukSurv_logwbc<-zresidual.coxph (fit_coxph = fit_LeukSurv_logwbc,
                                        traindata = LeukSurv,
                                        newdata = LeukSurv)
sw_zresid<-shapiro.test(resid_LeukSurv_logwbc$zresid)$p.value;sw_zresid
anov_zresid1<-test.nl.aov(qresidual=resid_LeukSurv_logwbc$zresid, 
                        fitted.values=fit_LeukSurv_logwbc$linear.predictors,
                        k.anova=10);anov_zresid1

anov_zresid_wbc_log1<-test.nl.aov(qresidual=resid_LeukSurv_logwbc$zresid, 
                                fitted.values=LeukSurv$wbc_log,
                                k.anova=10);anov_zresid_wbc_log1

anov_zresid_tpi<-test.nl.aov(qresidual=resid_LeukSurv_logwbc$zresid, 
                           fitted.values=LeukSurv$tpi,
                           k.anova=10);anov_zresid_tpi
anov_zresid_age<-test.nl.aov(qresidual=resid_LeukSurv_logwbc$zresid, 
                           fitted.values=LeukSurv$age,
                           k.anova=10);anov_zresid_age
censored<- LeukSurv$cens ==0
sf_czresid_logwbc<-as.numeric(gofTestCensored(resid_LeukSurv_logwbc$czresid,
                                           censored, test = "sf", 
                                           censoring.side = "right",
                                           distribution = "norm")$p.value)

km.ucs.LeukSurv_logwbc <- survfit(Surv(resid_LeukSurv_logwbc$ucs,
                                       LeukSurv$cens)~1,type='fleming')
id.ucs.LeukSurv_logwbc<-order(resid_LeukSurv_logwbc$ucs)


fit_LeukSurv_wbc  <- tryCatch(
  coxph(Surv(time, cens) ~ age  +sex+ wbc +tpi  +
          frailty(district, distribution="gamma"), data= LeukSurv),
  error = function(e) NA,
  warning = function(w) NA
)
AIC(fit_LeukSurv_wbc)
resid_LeukSurv_wbc<-zresidual.coxph (fit_coxph = fit_LeukSurv_wbc,
                                     traindata = LeukSurv,
                                     newdata = LeukSurv)
sw_zresid_t<-shapiro.test(resid_LeukSurv_wbc$zresid)$p.value;sw_zresid_t
anov_zresid_t1<-test.nl.aov(qresidual=resid_LeukSurv_wbc$zresid, 
                          fitted.values=fit_LeukSurv_wbc$linear.predictors,
                          k.anova=10);anov_zresid_t1


anov_zresid_wbc1<-test.nl.aov(qresidual=resid_LeukSurv_wbc$zresid, 
                            fitted.values=LeukSurv$wbc_log,
                            k.anova=10);anov_zresid_wbc1


anov_zresid_tpi<-test.nl.aov(qresidual=resid_LeukSurv_wbc$zresid, 
                           fitted.values=LeukSurv$tpi,
                           k.anova=10);anov_zresid_tpi
anov_zresid_age<-test.nl.aov(qresidual=resid_LeukSurv_wbc$zresid, 
                           fitted.values=LeukSurv$age,
                           k.anova=10);anov_zresid_age
censored<- LeukSurv$cens ==0
sf_czresid_wbc<-as.numeric(gofTestCensored(resid_LeukSurv_wbc$czresid,
                                        censored, test = "sf", 
                                        censoring.side = "right",
                                        distribution = "norm")$p.value)

km.ucs.LeukSurv_wbc <- survfit(Surv(resid_LeukSurv_wbc$ucs,LeukSurv$cens)~1,
                               type='fleming')
id.ucs.LeukSurv_wbc<-order(resid_LeukSurv_wbc$ucs)

