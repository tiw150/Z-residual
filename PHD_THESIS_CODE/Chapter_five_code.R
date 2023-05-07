##############Detection of Non-PH Due to Time-varying Covariate Effects#####
#set.seed(1)
rexp2 <- function(n, rate){ if (rate==0) rep(Inf,n) else rexp(n = n, rate = rate)}
########beta2 = beta_a if  0<surv_time< t0 , beta2 = beta_b if  t0 <surv_time 
simul_nonph1 <- function(n_clusters,n_individuals,fv,lambda, alpha, 
                         mean.censor,beta_a,beta_b,beta1,t0,corr)
{
  # Make covariates 
  n<-n_individuals * n_clusters
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  grpid<- paste0("g",grpid)
  sigma<-rbind(c(2,corr),c(corr,2))
  mu<-c(0,0) 
  covar<-mvrnorm(n=n_clusters*n_individuals, mu=mu, Sigma=sigma)
  x1 = round(covar[,1],1)
  x2 = round(covar[,2],1)
  # Make frailty (shape = a and scale = s,E(X) = a*s and Var(X) = a*s^2)
  frvec<-rep(rgamma(n_clusters,shape = 1/fv,scale = fv), each = n_individuals)
  #Event times:baseline hazard is Weibull，lambda=scale parameter, alpha=shape parameter
  v <- runif(n)
  a<- -log(v)
  b<- lambda*exp(x1*beta1+beta_a*x2)*frvec * (t0^alpha)
  surv_time<-rep(0,n)
  group1<-which(a<b)
  group2<-which(a>=b)
  surv_time[group1]<- (-log(v)[group1]/
                         (lambda*exp(x1[group1]*beta1+x2[group1]*beta_a)
                          *frvec[group1]))^(1 /alpha)
  
  surv_time[group2]<- ((-log(v)[group2]-
                          ((t0^alpha)*lambda*frvec[group2]*
                             exp(x1[group2]*beta1+x2[group2]*beta_a))+
                          ((t0^alpha)*lambda*frvec[group2]*
                             exp(x1[group2]*beta1+x2[group2]*beta_b)))
                       /(lambda*frvec[group2]*exp(x1[group2]*beta1)
                         *exp(x2[group2]*beta_b)))^(1 /alpha)
  # Censoring times
  c0<- rexp2(n, rate= 1/mean.censor)
  # follow-up times and event indicators
  t <- pmin(surv_time, c0)
  d <- as.numeric(surv_time <= c0)
  # Make data frame
  grpid<- as.factor(grpid)
  # x3<-as.factor(x3)
  out<-data.frame(id,grpid,x1,x2,t,d,frvec)
  return (out) 
}

####### coxph package  ####################################################
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


test.nl.aov1 <- function(qresidual, fitted.values, k.anova=10)
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


test.nl.aov1.categ <- function(qresidual, fitted.values)
{
  lpred.bin <- fitted.values
  anova(lm(qresidual ~ lpred.bin))$`Pr(>F)`[1]
}


test.var.bartl <- function(qresidual, fitted.values, k=10)
{
  lpred.bin <- cut(fitted.values, k)
  Z_group<- split(qresidual, lpred.bin)
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

test.var.bartl.categ <- function(qresidual, fitted.values)
{
  lpred.bin <- fitted.values
  Z_group<- split(qresidual, lpred.bin)
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

library("survival")
library("EnvStats")
library("nortest")
library("backports")
library("stringr")
library("timereg")
library("MASS")
library("dvmisc")

#########piecewise####################
####wrong model
nonph_data<-simul_nonph1(n_clusters=10,n_individuals=50,fv=0.5,
                         lambda=0.007, alpha=3, 
                         mean.censor=7.5,beta1=0.3,beta_a=1.35,
                         beta_b=-1.35,t0=2,corr=0)

table(nonph_data$d)
fit_coxph_nonph_data <- tryCatch(
  coxph(Surv(t, d) ~ x1+x2+frailty(grpid,distribution = "gamma"),
        data = nonph_data),
  error = function(e) NA,
  warning = function(w) NA
)
coxzph_test_w<-cox.zph(fit_coxph_nonph_data, transform="identity");coxzph_test_w
plot(coxzph_test_w, var = 1)
plot(coxzph_test_w, var = 2)
coxph_qr_w<-zresidual.coxph (fit_coxph = fit_coxph_nonph_data,
                             traindata = nonph_data,
                             newdata = nonph_data)
sw_coxph_w<-shapiro.test(coxph_qr_w)$p.value;sw_coxph_w

## the coefficient of x1 
anov_x1_w<- test.nl.aov1(qresidual=coxph_qr_w,
                         fitted.values=nonph_data$x1,
                         k.anova=10);anov_x1_w
bl_x1_w<-test.var.bartl(qresidual=coxph_qr_w,
                        fitted.values=nonph_data$x1,
                        k=10);bl_x1_w

## the coefficient of x2 
anov_x2_w<- test.nl.aov1(qresidual=coxph_qr_w,
                         fitted.values=nonph_data$x2,
                         k.anova=10);anov_x2_w
bl_x2_w<-test.var.bartl(qresidual=coxph_qr_w,
                        fitted.values=nonph_data$x2,
                        k=10);bl_x2_w

anov_lp_w<- test.nl.aov1(qresidual=coxph_qr_w,
                         fitted.values=fit_coxph_nonph_data$linear.predictors,
                         k.anova=10);anov_lp_w

bl_lp_w<-test.var.bartl(qresidual=coxph_qr_w,
                        fitted.values = fit_coxph_nonph_data$linear.predictors,
                        k=10);bl_lp_w

####true model
ph_data<-simul_nonph1(n_clusters=10,n_individuals=50,fv=0.5,
                      lambda=0.007, alpha=3,
                      mean.censor=7.5,beta1=0.3,beta_a=1.35,
                      beta_b=1.35,t0=pi/2,corr=0)
table(ph_data$d)
fit_coxph_ph_data <- tryCatch(
  coxph(Surv(t, d) ~ x1+x2 +frailty(grpid,distribution = "gamma"),
        data = ph_data),
  error = function(e) NA,
  warning = function(w) NA
)
coxzph_test_t<-cox.zph(fit_coxph_ph_data, transform="identity");coxzph_test_t
plot(coxzph_test_t, var = 1)
plot(coxzph_test_t, var = 2)


coxph_qr_t<-zresidual.coxph (fit_coxph = fit_coxph_ph_data,
                             traindata = ph_data,
                             newdata = ph_data)

sw_coxph_t<-shapiro.test(coxph_qr_t)$p.value;sw_coxph_t

## the coefficient of x1 
anov_x1_t<- test.nl.aov1(qresidual=coxph_qr_t,
                         fitted.values=ph_data$x1);anov_x1_t
bl_x1_t<-test.var.bartl(qresidual=coxph_qr_t,
                        fitted.values=ph_data$x1);bl_x1_t

## the coefficient of x2 
anov_x2_t<- test.nl.aov1(qresidual=coxph_qr_t,
                         fitted.values=ph_data$x2);anov_x2_t
bl_x2_t<-test.var.bartl(qresidual=coxph_qr_t,
                        fitted.values=ph_data$x2);bl_x2_t

anov_lp_t<- test.nl.aov1(qresidual=coxph_qr_t,
                         fitted.values=fit_coxph_ph_data$linear.predictors)
anov_lp_t
bl_lp_t<-test.var.bartl(qresidual=coxph_qr_t,
                        fitted.values = fit_coxph_ph_data$linear.predictors)
bl_lp_t


##########################################################################
##############Detection of Non-PH Due to Accelerated Failure#############
library("survival")
rexp2 <- function(n, rate){ if (rate==0) rep(Inf,n) else rexp(n=n, rate = rate)}
AFTLN_data <- function(n_clusters,n_individuals, mu, sigma,beta1,
                       mean.censor,fv)
{
  # Make covariates 
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  grpid<- paste0("g",grpid)
  x1 = round(rnorm(n_clusters*n_individuals, 0,1),2)
  
  #  x2 = runif(n_clusters*n_individuals,0,1)
  #  x3 = rbinom(n_clusters*n_individuals, size = 1, p = 0.5)
  
  # Make frailty (shape = a and scale = s,E(X) = a*s and Var(X) = a*s^2)
  #frvec<-rep(rgamma(n_clusters,shape=0.5, rate=1), each = n_individuals)
  frvec<-rep(rgamma(n_clusters,shape = 1/fv,scale = fv), each = n_individuals)
  
  #Event times:baseline hazard is AFT lognormal
  n<-n_individuals * n_clusters
  v <- runif(n)
  surv_time <- exp(sigma* qnorm(v)+mu)*(1/exp(x1*beta1+log(frvec)))
  # Censoring times
  t0<- rexp2(n, rate= 1/mean.censor)
  # follow-up times and event indicators
  t <- pmin(surv_time, t0)
  d <- as.numeric(surv_time <= t0)
  # Make data frame
  grpid<- as.factor(grpid)
  #  x3<-as.factor(x3)
  out<-data.frame(id,grpid,x1,t,d,frvec)
  return (out) 
}

AFTWB_data <- function(n_clusters,n_individuals, lambda, alpha,beta1,
                       mean.censor,fv)
{
  # Make covariates 
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  grpid<- paste0("g",grpid)
  x1 = round(rnorm(n_clusters*n_individuals, 0,1),2)
  #  x1 = runif(n_clusters*n_individuals,0,1)
  #  x3 = rbinom(n_clusters*n_individuals, size = 1, p = 0.5)
  
  # Make frailty (shape = a and scale = s,E(X) = a*s and Var(X) = a*s^2)
  #frvec<-rep(rgamma(n_clusters,shape=0.5, rate=1), each = n_individuals)
  frvec<-rep(rgamma(n_clusters,shape = 1/fv,scale = fv), each = n_individuals)
  
  #Event times:baseline hazard is AFT weibull
  #(baseline hazard is Weibull，lambda=scale parameter, alpha=shape parameter)
  n<-n_individuals * n_clusters
  v <- runif(n)
  surv_time <- (-log(1-v)/lambda)^(1/alpha)*(exp(-(x1*beta1+log(frvec))))
  # Censoring times
  t0<- rexp2(n, rate= 1/mean.censor)
  # follow-up times and event indicators
  t <- pmin(surv_time, t0)
  d <- as.numeric(surv_time <= t0)
  # Make data frame
  grpid<- as.factor(grpid)
  #  x3<-as.factor(x3)
  out<-data.frame(id,grpid,x1,t,d,frvec)
  return (out) 
}


aft_ln_data<-AFTLN_data(n_clusters=10,n_individuals=50,fv=0.5,
                        mu=0, sigma=1,mean.censor=3.8,beta1=1)
table(aft_ln_data$d)

fit_coxph_aftln_data <- tryCatch(
  coxph(Surv(t, d) ~ x1+frailty(grpid,distribution = "gamma"),
        data = aft_ln_data),
  error = function(e) NA,
  warning = function(w) NA
)

coxzph_test_w<-cox.zph(fit_coxph_aftln_data,transform="identity");coxzph_test_w
#plot(score_test)
coxph_qr_w<-zresidual.coxph (fit_coxph = fit_coxph_aftln_data,
                             traindata = aft_ln_data,
                             newdata = aft_ln_data)

sw_coxph_w<-shapiro.test(coxph_qr_w)$p.value;sw_coxph_w
anov_lp_w<- test.nl.aov1(qresidual=coxph_qr_w,
                         fitted.values=fit_coxph_aftln_data$linear.predictors);anov_lp_w
anov_x1_w<- test.nl.aov1(qresidual=coxph_qr_w,
                         fitted.values=aft_ln_data$x1);anov_x1_w

bl_lp_w<-test.var.bartl(qresidual=coxph_qr_w,
                        fitted.values = fit_coxph_aftln_data$linear.predictors);bl_lp_w
bl_x1_w<-test.var.bartl(qresidual=coxph_qr_w,
                        fitted.values = aft_ln_data$x1);bl_x1_w


aft_wb_data<-AFTWB_data(n_clusters=10,n_individuals=50,fv=0.5,
                        lambda=1.67 , alpha=0.48 ,mean.censor=0.2,beta1=1)
table(aft_wb_data$d)
fit_coxph_aftwb_data <- tryCatch(
  coxph(Surv(t, d) ~ x1+frailty(grpid,distribution = "gamma"),
        data = aft_wb_data),
  error = function(e) NA,
  warning = function(w) NA
)

coxzph_test_t<-cox.zph(fit_coxph_aftwb_data,transform="identity");coxzph_test_t

coxph_qr_t<-zresidual.coxph (fit_coxph = fit_coxph_aftwb_data,
                             traindata = aft_wb_data,
                             newdata = aft_wb_data)

sw_coxph_t<-shapiro.test(coxph_qr_t)$p.value;sw_coxph_t
anov_lp_t<- test.nl.aov1(qresidual=coxph_qr_t,
                         fitted.values=fit_coxph_aftwb_data$linear.predictors);anov_lp_t
anov_x1_t<- test.nl.aov1(qresidual=coxph_qr_t,
                         fitted.values=aft_wb_data$x1);anov_x1_t

bl_lp_t<-test.var.bartl(qresidual=coxph_qr_t,
                        fitted.values = fit_coxph_aftwb_data$linear.predictors);bl_lp_t
bl_x1_t<-test.var.bartl(qresidual=coxph_qr_t,
                        fitted.values = aft_wb_data$x1);bl_x1_t

##########################################################################
################real data example ########################################
library("survival")
library("EnvStats")
library("nortest")
library("backports")
library("stringr")
library("timereg")
library("MASS")
library("dvmisc")
library("frailtySurv")
############################################################
data("drs")
drs$subject_id<-as.factor(drs$subject_id)
drs$treated<-as.factor(drs$treated)
drs$laser_type<-as.factor(drs$laser_type)
drs$diabetes_type<-as.factor(drs$diabetes_type)
#colnames(drs)<-c("subject","eye","time","status","Treated","Age","Laser Type", "Diabetes Type")
############ fit Cox PH linear model
fit_drs1  <- tryCatch(
  coxph(Surv(time, status) ~ treated + age_at_onset +
          laser_type+ diabetes_type+
          frailty(subject_id, distribution="gamma"), data= drs),
  error = function(e) NA,
  warning = function(w) NA
);fit_drs1

AIC(fit_drs1)

coxzph_test_drs<-cox.zph(fit_drs1, transform="identity");coxzph_test_drs
qr_drs1<-qresidual.coxph (fit_coxph = fit_drs1,
                          traindata = drs,newdata = drs)
resids_drs1<-allresidual.coxph (fit_coxph = fit_drs1,
                                traindata = drs, newdata = drs)
censored_cox<- drs$status ==0
gof_fit_drs1<-gofTestCensored(resids_drs1$nusp,censored_cox, test = "sf", 
                              censoring.side = "right",
                              distribution = "norm")$p.value
gof_fit_drs1
sw_qr<-shapiro.test(qr_drs1)$p.value;sw_qr
anov_drs<-test.nl.aov1(qresidual=qr_drs1, 
                       fitted.values=fit_drs1$linear.predictors,
                       k.anova=10);anov_drs
bl_drs<-test.var.bartl(qresidual=qr_drs1,
                       fitted.values=fit_drs1$linear.predictors);bl_drs
anov_qr_age<-test.nl.aov1(qresidual=qr_drs1, 
                          fitted.values=drs$age_at_onset,
                          k.anova=10);anov_qr_age
bl_qr_age<-test.var.bartl(qresidual=qr_drs1, 
                          fitted.values=drs$age_at_onset,
                          k=10);bl_qr_age
anov_qr_laser<-test.nl.aov1.categ(qresidual=qr_drs1, 
                                  fitted.values=drs$laser_type);anov_qr_laser
bl_qr_laser<-test.var.bartl.categ(qresidual=qr_drs1, 
                                  fitted.values=drs$laser_type);bl_qr_laser
anov_qr_diabete<-test.nl.aov1.categ(qresidual=qr_drs1, 
                                    fitted.values=drs$diabetes_type);anov_qr_diabete
bl_qr_diabete<-test.var.bartl.categ(qresidual=qr_drs1, 
                                    fitted.values=drs$diabetes_type);bl_qr_diabete
anov_qr_trt<-test.nl.aov1.categ(qresidual=qr_drs1, 
                                fitted.values=drs$treated);anov_qr_trt
bl_qr_trt<-test.var.bartl.categ(qresidual=qr_drs1, 
                                fitted.values=drs$treated);bl_qr_trt
#####fit AFT Lognormal model
fit_drs_survreg<-
  survreg(Surv(time, status) ~ treated + age_at_onset +
            laser_type+ diabetes_type, data= drs,dist="lognormal")
AIC(fit_drs_survreg)
summary(fit_drs_survreg)
qr_drs_aft<-(resid_survreg (survreg_fit = fit_drs_survreg))$nrsp
nusp_drs_aft<-(resid_survreg (survreg_fit = fit_drs_survreg))$nusp
sw_qr_aft<-shapiro.test(qr_drs_aft)$p.value;sw_qr_aft
censored_aft<- drs$status ==0
gof_fit_drs1_aft<-gofTestCensored(nusp_drs_aft,censored_aft, test = "sf", 
                                  censoring.side = "right",
                                  distribution = "norm")$p.value
gof_fit_drs1_aft
anov_drs_aft<-test.nl.aov1(qresidual=qr_drs_aft, 
                           fitted.values=fit_drs_survreg$linear.predictors,
                           k.anova=10);anov_drs_aft
bl_drs_aft<-test.var.bartl(qresidual=qr_drs_aft,
                           fitted.values=fit_drs_survreg$linear.predictors);bl_drs_aft
anov_qr_age_aft<-test.nl.aov1(qresidual=qr_drs_aft, 
                              fitted.values=drs$age_at_onset,
                              k.anova=10);anov_qr_age_aft
bl_qr_age_aft<-test.var.bartl(qresidual=qr_drs_aft, 
                              fitted.values=drs$age_at_onset,
                              k=10);bl_qr_age_aft
anov_qr_laser_aft<-test.nl.aov1.categ(qresidual=qr_drs_aft, 
                                      fitted.values=drs$laser_type);anov_qr_laser_aft
bl_qr_laser_aft<-test.var.bartl.categ(qresidual=qr_drs_aft, 
                                      fitted.values=drs$laser_type);bl_qr_laser_aft
anov_qr_diabete_aft<-test.nl.aov1.categ(qresidual=qr_drs_aft, 
                                        fitted.values=drs$diabetes_type);anov_qr_diabete_aft
bl_qr_diabete_aft<-test.var.bartl.categ(qresidual=qr_drs_aft, 
                                        fitted.values=drs$diabetes_type);bl_qr_diabete_aft
anov_qr_trt_aft<-test.nl.aov1.categ(qresidual=qr_drs_aft, 
                                    fitted.values=drs$treated);anov_qr_trt_aft
bl_qr_trt_aft<-test.var.bartl.categ(qresidual=qr_drs_aft, 
                                    fitted.values=drs$treated);bl_qr_trt_aft





