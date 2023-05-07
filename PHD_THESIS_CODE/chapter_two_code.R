ffrailtyEM <- function(data, null_obj) {
  # fit model
  cc <- frailtyEM::emfrail_control(se_adj = FALSE, ca_test = FALSE)
  fit <- tryCatch(
    frailtyEM::emfrail(
      Surv(t, d) ~ x1+x2+x3 + cluster(grpid),
      distribution = emfrail_dist(dist = "gamma"),
      data = data,
      control = cc
    ),
    error = function(e)
      return(null_obj),
    warning = function(w)
      return(null_obj)
  )
  
  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }
  
  # make summary object
  smr <- tryCatch(
    summary(fit),
    error = function(e)
      return(null_obj)
  )
  
  # break out of the function if the fit converged to bad values (e.g. hessian not invertible)
  if ("null_obj" %in% class(smr)) {
    return(null_obj)
  }
  
  # if the object has NaN values, break out and return the null object
  if (any(is.nan(smr[["coefmat"]][, "coef"]))) {
    return(null_obj)
  }
  if (any(is.nan(smr[["coefmat"]][, "se(coef)"]))) {
    return(null_obj)
  }
  if (any(is.nan(smr[["fr_var"]]["fr_var"]))) {
    return(null_obj)
  }
  if (any(is.nan(smr[["fr_var"]]["se_fr_var"]))) {
    return(null_obj)
  }
  
  #confidence interval for regression coefficients and frailty variance
  x1_lb<-(coef(fit)["x1"][[1]])-1.96*(smr[["coefmat"]]["x1", "se(coef)"])
  x1_ub<-(coef(fit)["x1"][[1]])+1.96*(smr[["coefmat"]]["x1", "se(coef)"])
  x2_lb<-(coef(fit)["x2"][[1]])-1.96*(smr[["coefmat"]]["x2", "se(coef)"])
  x2_ub<-(coef(fit)["x2"][[1]])+1.96*(smr[["coefmat"]]["x2", "se(coef)"])
  x3_lb<-(coef(fit)["x3"][[1]])-1.96*(smr[["coefmat"]]["x3", "se(coef)"])
  x3_ub<-(coef(fit)["x3"][[1]])+1.96*(smr[["coefmat"]]["x3", "se(coef)"])
  fv_lb<-(smr[["fr_var"]]["fr_var"][[1]])-1.96*(smr[["fr_var"]]["se_fr_var"][[1]])
  fv_ub<-(smr[["fr_var"]]["fr_var"][[1]])+1.96*(smr[["fr_var"]]["se_fr_var"][[1]])
  # CI symmetric on log(theta) scale
  var_logtheta<-(smr[["fr_var"]]["se_fr_var"][[1]])^2/(smr[["fr_var"]]["fr_var"][[1]])^2
  fv_llb<-smr[["fr_var"]]["ci_frvar_low"][[1]]
  fv_lub<-smr[["fr_var"]]["ci_frvar_high"][[1]]  
  
  # make object to return
  out <- data.frame(
    par = c("x1","x2","x3", "fv","fv_log(theta)"),
    coef = c(coef(fit)["x1"][[1]],coef(fit)["x2"][[1]],coef(fit)["x3"][[1]], 
             smr[["fr_var"]]["fr_var"][[1]],log(smr[["fr_var"]]["fr_var"][[1]])),
    se = c(smr[["coefmat"]]["x1", "se(coef)"],smr[["coefmat"]]["x2", "se(coef)"],
           smr[["coefmat"]]["x3", "se(coef)"], smr[["fr_var"]]["se_fr_var"][[1]],
           sqrt(var_logtheta)),
    CI_LB=c(x1_lb,x2_lb,x3_lb,fv_lb,fv_llb),
    CI_UB=c(x1_ub,x2_ub,x3_ub,fv_ub,fv_lub)
  )
  
  return(out)
}

fparfm <- function(data, null_obj) {
  # fit model
  fit <- tryCatch(
    parfm::parfm(
      Surv(t, d) ~ x1+x2+x3,
      cluster = "grpid",
      data = data,
      dist="weibull",
      frailty="gamma"
    ),
    error = function(e)
      return(null_obj),
    warning = function(w)
      return(null_obj)
  )
  
  # there is some weird bug in parfm...
  # sink()
  
  
  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }
  
  x1_lb<-log(ci.parfm(fit, level=0.05)["x1","low"])
  x2_lb<-log(ci.parfm(fit, level=0.05)["x2","low"])
  x3_lb<-log(ci.parfm(fit, level=0.05)["x3","low"])
  x1_ub<-log(ci.parfm(fit, level=0.05)["x1","up"])
  x2_ub<-log(ci.parfm(fit, level=0.05)["x2","up"])
  x3_ub<-log(ci.parfm(fit, level=0.05)["x3","up"])
  fv_lb<-(fit["theta","ESTIMATE"])-1.96*(fit["theta", "SE"])
  fv_ub<-(fit["theta","ESTIMATE"])+1.96*(fit["theta", "SE"])
  # CI symmetric on log(theta) scale
  var_logtheta<-(fit["theta", "SE"])^2 / (fit["theta","ESTIMATE"])^2
  fv_llb<-exp(log(fit["theta","ESTIMATE"])-1.96*sqrt(var_logtheta))
  fv_lub<-exp(log(fit["theta","ESTIMATE"])+1.96*sqrt(var_logtheta))
  
  # make object to return
  out <- data.frame(
    par = c("x1","x2","x3", "fv","fv_log(theta)"),
    coef = c(coef(fit)[["x1"]],coef(fit)[["x2"]],coef(fit)[["x3"]], 
             fit["theta","ESTIMATE"],log(fit["theta","ESTIMATE"])),
    se = c(fit["x1", "SE"],fit["x2", "SE"],fit["x3", "SE"],
           fit["theta", "SE"],sqrt(var_logtheta)),
    CI_LB=c(x1_lb,x2_lb,x3_lb,fv_lb,fv_llb),
    CI_UB=c(x1_ub,x2_ub,x3_ub,fv_ub,fv_lub)
  )
  return(out)
}

ffrailtypack <- function(data, null_obj) {
  # fit model
  fit <- tryCatch(
    frailtypack::frailtyPenal(
      Surv(t, d) ~ x1+x2+x3 + cluster(grpid),
      data = data,
      RandDist = "Gamma",
      n.knots=15,
      kappa=1,
      cross.validation=TRUE,
      print.times = FALSE
    ),
    error = function(e)
      return(null_obj),
    warning = function(w)
      return(null_obj)
  )
  
  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }
  
  x1_lb<-(as.numeric(fit[["coef"]]["x1"]))-1.96*(sqrt(fit[["varH"]][1,1]))
  x2_lb<-(as.numeric(fit[["coef"]]["x2"]))-1.96*(sqrt(fit[["varH"]][2,2]))
  x3_lb<-(as.numeric(fit[["coef"]]["x3"]))-1.96*(sqrt(fit[["varH"]][3,3]))
  x1_ub<-(as.numeric(fit[["coef"]]["x1"]))+1.96*(sqrt(fit[["varH"]][1,1]))
  x2_ub<-(as.numeric(fit[["coef"]]["x2"]))+1.96*(sqrt(fit[["varH"]][2,2]))
  x3_ub<-(as.numeric(fit[["coef"]]["x3"]))+1.96*(sqrt(fit[["varH"]][3,3]))
  fv_lb<-(fit[["theta"]])-1.96*(sqrt(fit[["varTheta"]][1]))
  fv_ub<-(fit[["theta"]])+1.96*(sqrt(fit[["varTheta"]][1]))
  # CI symmetric on log(theta) scale
  var_logtheta<-fit[["varTheta"]][1] / (fit[["theta"]])^2
  fv_llb<-exp(log(fit[["theta"]])-1.96*sqrt(var_logtheta))
  fv_lub<-exp(log(fit[["theta"]])+1.96*sqrt(var_logtheta))
  
  
  # make object to return
  out <- data.frame(
    par = c("x1","x2","x3", "fv","fv_log(theta)"),
    coef = c(fit[["coef"]]["x1"], fit[["coef"]]["x2"],fit[["coef"]]["x3"], 
             fit[["theta"]],log(fit[["theta"]])),
    se = c(sqrt(fit[["varH"]][1,1]),sqrt(fit[["varH"]][2,2]),sqrt(fit[["varH"]][3,3]),
           sqrt(fit[["varTheta"]][1]),sqrt(var_logtheta)),
    CI_LB=c(x1_lb,x2_lb,x3_lb,fv_lb,fv_llb),
    CI_UB=c(x1_ub,x2_ub,x3_ub,fv_ub,fv_lub)
  )
  return(out)
}


fcoxph <- function(data, null_obj) {
  # fit model
  fit <- tryCatch(
    survival::coxph(Surv(t, d) ~ x1+x2+x3 +frailty(grpid,distribution = "gamma"),
                    data = data),
    error = function(e)
      return(null_obj),
    warning = function(w)
      return(null_obj)
  )
  
  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }
  
  # make summary object
  smr <- tryCatch(
    summary(fit),
    error = function(e)
      return(null_obj)
  )
  
  # break out of the function if the fit converged to bad values (e.g. hessian not invertible)
  if ("null_obj" %in% class(smr)) {
    return(null_obj)
  }
  
  x1_lb<-log(smr$conf.int[["x1","lower .95"]])
  x2_lb<-log(smr$conf.int[["x2","lower .95"]])
  x3_lb<-log(smr$conf.int[["x3","lower .95"]])
  x1_ub<-log(smr$conf.int[["x1","upper .95"]])
  x2_ub<-log(smr$conf.int[["x2","upper .95"]])
  x3_ub<-log(smr$conf.int[["x3","upper .95"]])
  # make object to return
  out <- data.frame(
    par = c("x1","x2","x3", "fv","fv_log(theta)"),
    coef = c(coef(fit)[["x1"]],coef(fit)[["x2"]],coef(fit)[["x3"]], 
             as.numeric(substring(summary(fit)$print2, first=28,last=36)),NA),
    se = c(smr[["coefficients"]]["x1","se(coef)"],smr[["coefficients"]]["x2","se(coef)"], smr[["coefficients"]]["x3","se(coef)"],NA,NA),
    CI_LB=c(x1_lb,x2_lb,x3_lb,NA,NA),
    CI_UB=c(x1_ub,x2_ub,x3_ub,NA,NA)
  )
  return(out)
}


ffrailtysurv <- function(data, null_obj) {
  # fit model
  fit <- tryCatch(
    frailtySurv::fitfrail(
      Surv(t, d) ~ x1+x2+x3 + cluster(grpid),
      data,
      frailty="gamma"
    ),
    error = function(e)
      return(null_obj),
    warning = function(w)
      return(null_obj)
  )
  
  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }
  
  x1_lb<-(fit$beta[["x1"]])-1.96*(sqrt(diag(vcov(fit)))[["x1"]])
  x1_ub<-(fit$beta[["x1"]])+1.96*(sqrt(diag(vcov(fit)))[["x1"]])
  x2_lb<-(fit$beta[["x2"]])-1.96*(sqrt(diag(vcov(fit)))[["x2"]])
  x2_ub<-(fit$beta[["x2"]])+1.96*(sqrt(diag(vcov(fit)))[["x2"]])
  x3_lb<-(fit$beta[["x3"]])-1.96*(sqrt(diag(vcov(fit)))[["x3"]])
  x3_ub<-(fit$beta[["x3"]])+1.96*(sqrt(diag(vcov(fit)))[["x3"]])
  fv_lb<-(fit$theta)-1.96*(sqrt(diag(vcov(fit)))[["theta.1"]])
  fv_ub<-(fit$theta)+1.96*(sqrt(diag(vcov(fit)))[["theta.1"]])
  # CI symmetric on log(theta) scale
  var_logtheta<-(sqrt(diag(vcov(fit)))[["theta.1"]])^2 / (fit$theta)^2
  fv_llb<-exp(log(fit$theta)-1.96*sqrt(var_logtheta))
  fv_lub<-exp(log(fit$theta)+1.96*sqrt(var_logtheta))
  
  # make object to return
  out <- data.frame(
    par = c("x1","x2","x3", "fv","fv_log(theta)"),
    coef = c(fit$beta[["x1"]],fit$beta[["x2"]],fit$beta[["x3"]], 
             fit$theta,log(fit$theta)),
    se = c(sqrt(diag(vcov(fit)))[["x1"]],sqrt(diag(vcov(fit)))[["x2"]],
           sqrt(diag(vcov(fit)))[["x3"]], sqrt(diag(vcov(fit)))[["theta.1"]],
           sqrt(var_logtheta)),
    CI_LB=c(x1_lb,x2_lb,x3_lb,fv_lb,fv_llb),
    CI_UB=c(x1_ub,x2_ub,x3_ub,fv_ub,fv_lub)
  )
  return(out)
}


ffrailtyHL <- function(data, null_obj) {
  # fit model
  fit <- tryCatch(
    frailtyHL::frailtyHL(
      Surv(t,d)~x1+x2+x3+(1|grpid),
      data=data, RandDist = "Gamma",
      mord = 1, dord = 2
    ),
    error = function(e)
      return(null_obj)
  )
  
  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }
  
  x1_lb<-(fit$FixCoef[["x1","Estimate"]])-1.96*(fit$FixCoef[["x1","Std. Error"]])
  x1_ub<-(fit$FixCoef[["x1","Estimate"]])+1.96*(fit$FixCoef[["x1","Std. Error"]])
  x2_lb<-(fit$FixCoef[["x2","Estimate"]])-1.96*(fit$FixCoef[["x2","Std. Error"]])
  x2_ub<-(fit$FixCoef[["x2","Estimate"]])+1.96*(fit$FixCoef[["x2","Std. Error"]])
  x3_lb<-(fit$FixCoef[["x3","Estimate"]])-1.96*(fit$FixCoef[["x3","Std. Error"]])
  x3_ub<-(fit$FixCoef[["x3","Estimate"]])+1.96*(fit$FixCoef[["x3","Std. Error"]])
  fv_lb<-(fit$RandCoef[["grpid","Estimate"]])-1.96*(fit$RandCoef[["grpid","Std. Error"]])
  fv_ub<-(fit$RandCoef[["grpid","Estimate"]])+1.96*(fit$RandCoef[["grpid","Std. Error"]])
  # CI symmetric on log(theta) scale
  var_logtheta<-(fit$RandCoef[["grpid","Std. Error"]])^2/(fit$RandCoef[["grpid","Estimate"]])^2
  fv_llb<-exp(log(fit$RandCoef[["grpid","Estimate"]])-1.96*sqrt(var_logtheta))
  fv_lub<-exp(log(fit$RandCoef[["grpid","Estimate"]])+1.96*sqrt(var_logtheta))
  
  # make object to return
  out <- data.frame(
    par = c("x1","x2","x3", "fv","fv_log(theta)"),
    coef = c(fit$FixCoef[["x1","Estimate"]],fit$FixCoef[["x2","Estimate"]],
             fit$FixCoef[["x3","Estimate"]], fit$RandCoef[["grpid","Estimate"]],
             log(fit$RandCoef[["grpid","Estimate"]])),
    se = c(fit$FixCoef[["x1","Std. Error"]],fit$FixCoef[["x2","Std. Error"]],
           fit$FixCoef[["x3","Std. Error"]],fit$RandCoef[["grpid","Std. Error"]],
           sqrt(var_logtheta)),
    CI_LB=c(x1_lb,x2_lb,x3_lb,fv_lb,fv_llb),
    CI_UB=c(x1_ub,x2_ub,x3_ub,fv_ub,fv_lub)
  )
  return(out)
}

make_models <- function(data) {
  # Null object to return in case of errors:
  null_obj <- data.frame(
    par = c("x1","x2","x3", "fv","fv_log(theta)"),
    coef = rep(NA, 5),
    se = rep(NA, 5),
    CI_LB = rep(NA, 5),
    CI_UP = rep(NA, 5)
  )
  class(null_obj) <- c(class(null_obj), "null_obj")
  
  # Required packages:
  if (!requireNamespace("pacman")) {
    install.packages("pacman")
  }
  pacman::p_load(
    "dplyr",
    "tidyr",
    "survival",
    "parfm",
    "frailtyEM",
    "frailtypack",
    "frailtySurv",
    "frailtyHL"
  )
  
  # Models:
  tstart_em <- Sys.time()
  EM_out <- ffrailtyEM(
    data = data,
    null_obj = null_obj
  )
  tstop_em <- Sys.time()
  
  tstart_par <- Sys.time()
  Par_out<-fparfm(
    data = data,
    null_obj = null_obj
  )
  tstop_par <- Sys.time()
  
  tstart_fpack <- Sys.time()
  fpack_out<-ffrailtypack(
    data = data,
    null_obj = null_obj
  )
  tstop_fpack <- Sys.time()
  
  tstart_coxph <- Sys.time()
  coxph_out<-fcoxph(
    data = data,
    null_obj = null_obj
  )
  tstop_coxph <- Sys.time()
  
  tstart_surv <- Sys.time()
  surv_out<-ffrailtysurv(
    data,
    null_obj = null_obj
  )
  tstop_surv <- Sys.time()
  
  tstart_HL <- Sys.time()
  HL_out<-ffrailtyHL(
    data=data,
    null_obj = null_obj
  )
  tstop_HL <- Sys.time()
  
  # Add time info to return
  EM_out[["time"]] <- difftime(tstop_em, tstart_em, units = "mins")
  Par_out[["time"]] <- difftime(tstop_par, tstart_par, units = "mins")
  fpack_out[["time"]] <- difftime(tstop_fpack, tstart_fpack, units = "mins")
  coxph_out[["time"]] <- difftime(tstop_coxph, tstart_coxph, units = "mins")
  surv_out[["time"]] <- difftime(tstop_surv, tstart_surv, units = "mins")
  HL_out[["time"]] <- difftime(tstop_HL, tstart_HL, units = "mins")
  
  # Return results
  out<-list(EM_out=EM_out,Par_out=Par_out,fpack_out=fpack_out,coxph_out=coxph_out,
            surv_out=surv_out,HL_out=HL_out)
  return(out)
}

############ Simulate data############################
set.seed(1)
rexp2 <- function(n, rate){ if (rate==0) rep(Inf,n) else rexp(n = n, rate = rate)}
simulWeib <- function(n_clusters,n_individuals, lambda, alpha,beta1,beta2,beta3,
                      mean.censor,fv)
{
  # Make covariates 
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  x1 = rbinom(n_clusters*n_individuals, size = 1, p = 0.5)
  x2 = round(runif(n_clusters*n_individuals,10,60),0)
  x3 = round(rnorm(n_clusters*n_individuals, 0,1),2)
  
  # Make frailty (shape = a and scale = s,E(X) = a*s and Var(X) = a*s^2)
  #frvec<-rep(rgamma(n_clusters,shape=0.5, rate=1), each = n_individuals)
  frvec<-rep(rgamma(n_clusters,shape = 1/fv,scale = fv), each = n_individuals)
  
  #Event times:baseline hazard is Weibull，lambda=scale parameter, alpha=shape parameter
  n<-n_individuals * n_clusters
  v <- runif(n)
  surv_time <-
    (- log(v)/(lambda*exp(x1*beta1+x2*beta2+x3*beta3+frvec)))^(1 /alpha)
  # Censoring times
  t0<- rexp2(n, rate= 1/mean.censor)
  # follow-up times and event indicators
  t <- pmin(surv_time, t0)
  d <- as.numeric(surv_time <= t0)
  
  # Make data frame
  out<-data.frame(id,grpid,x1,x2,x3,t,d,frvec)
  return (out) 
}

para.frame <- data.frame(
  n_individuals = c(10,10,10,10,10,10,10,10,10,40,40,40,80,80,80),
  n_clusters = c(10,10,10,40,40,40,80,80,80,10,10,10,10,10,10),
  beta1 = rep(0.2,times=15),
  beta2 = rep(0.1,times=15),
  beta3 = rep(0.3,times=15),
  fv =  rep(0.5,times=15),
  lambda=rep(1,times=15),
  alpha=rep(1.8,times=15),
  mean.censor=c(0.355,0.08,0.02, 0.45 ,0.09 ,0.023,0.45,0.095,
                0.024,0.34 ,0.076 ,0.016,0.43,0.1,0.024),
  c=rep(c(20,50,80),times=5)
)

i<-1 #### i is equal to 1 to 15.
para<-para.frame[i,]
n_individuals<-para[1]
n_clusters<-para[2]
beta1<-para[3]
beta2<-para[4]
beta3<-para[5]
fv<-para[6]
lambda<-para[7]
alpha<-para[8]
mean.censor<-para[9]

library("foreach")
library("doParallel")
n_sims<-1000
cur_time = proc.time()
out <- vector(mode = "list", length = 3) 
foreach (j = 1:n_sims, .combine = rbind) %do% 
  { 
    cat(paste('Simulation ',j,' out of ',n_sims,'\n'))
    if(j ==2){
      elapsed=as.numeric(proc.time()-cur_time)[3]
      cat(paste("Time for 1 simulation: ",elapsed/3600," hours \n"))
      cat(paste("Estimated time remaining: ",elapsed/3600*(n_sims-1)," hours \n"))
    }
    simulated_shared_frailty_data<-
      simulWeib(n_clusters=n_clusters[1,],n_individuals=n_individuals[1,],
                lambda=lambda[1,],alpha=alpha[1,],beta1=beta1[1,],
                beta2=beta2[1,],beta3=beta3[1,], 
                mean.censor= mean.censor[1,],fv=fv[1,])
    out[[j]]<-make_models(simulated_shared_frailty_data)
  }

###save each i's result 
saveRDS(object = out, file= "result.RData")

#####Make a data.frame (True value,Mean_coef,Bias,Mean_SE,Empirical_SE, 
####Median_coef,Median_se,CP,CP_logTheta,Converged,time) for each result 

out<-readRDS(file="result.RData")
n_sims<-1000
x1=1
x2=-1
x3=0.5
fv=0.5
Colnames=c("EM_out_par","coef","se","CI_LB","CI_UB","time")
EM_x1=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(EM_x1) = Colnames
EM_x2=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(EM_x2) = Colnames
EM_x3=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(EM_x3) = Colnames
EM_fv=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(EM_fv) = Colnames
EM_fvl=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(EM_fvl) = Colnames

Colnames=c("Par_out_par","coef","se","CI_LB","CI_UB","time")
Par_x1=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(Par_x1) = Colnames
Par_x2=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(Par_x2) = Colnames
Par_x3=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(Par_x3) = Colnames
Par_fv=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(Par_fv) = Colnames
Par_fvl=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(Par_fvl) = Colnames

Colnames=c("fpack_out_par","coef","se","CI_LB","CI_UB","time")
fpack_x1=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(fpack_x1) = Colnames
fpack_x2=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(fpack_x2) = Colnames
fpack_x3=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(fpack_x3) = Colnames
fpack_fv=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(fpack_fv) = Colnames
fpack_fvl=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(fpack_fvl) = Colnames

Colnames=c("coxph_out_par","coef","se","CI_LB","CI_UB","time")
coxph_x1=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(coxph_x1) = Colnames
coxph_x2=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(coxph_x2) = Colnames
coxph_x3=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(coxph_x3) = Colnames
coxph_fv=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(coxph_fv) = Colnames
coxph_fvl=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(coxph_fvl) = Colnames


Colnames=c("surv_out_par","coef","se","CI_LB","CI_UB","time")
surv_x1=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(surv_x1) = Colnames
surv_x2=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(surv_x2) = Colnames
surv_x3=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(surv_x3) = Colnames
surv_fv=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(surv_fv) = Colnames
surv_fvl=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(surv_fvl) = Colnames


Colnames=c("HL_out_par","coef","se","CI_LB","CI_UB","time")
HL_x1=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(HL_x1) = Colnames
HL_x2=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(HL_x2) = Colnames
HL_x3=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(HL_x3) = Colnames
HL_fv=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(HL_fv) = Colnames
HL_fvl=data.frame(matrix(0,n_sims,length(Colnames)))
colnames(HL_fvl) = Colnames

for(i in 1:n_sims){
  EM_x1[i,]<-out[[i]]$EM_out[1,]
  EM_x2[i,]<-out[[i]]$EM_out[2,]
  EM_x3[i,]<-out[[i]]$EM_out[3,]
  EM_fv[i,]<-out[[i]]$EM_out[4,]
  EM_fvl[i,]<-out[[i]]$EM_out[5,]
  
  Par_x1[i,]<-out[[i]]$Par_out[1,]
  Par_x2[i,]<-out[[i]]$Par_out[2,]
  Par_x3[i,]<-out[[i]]$Par_out[3,]
  Par_fv[i,]<-out[[i]]$Par_out[4,]
  Par_fvl[i,]<-out[[i]]$Par_out[5,]
  
  fpack_x1[i,]<-out[[i]]$fpack_out[1,]
  fpack_x2[i,]<-out[[i]]$fpack_out[2,]
  fpack_x3[i,]<-out[[i]]$fpack_out[3,]
  fpack_fv[i,]<-out[[i]]$fpack_out[4,]
  fpack_fvl[i,]<-out[[i]]$fpack_out[5,]
  
  coxph_x1[i,]<-out[[i]]$coxph_out[1,]
  coxph_x2[i,]<-out[[i]]$coxph_out[2,]
  coxph_x3[i,]<-out[[i]]$coxph_out[3,]
  coxph_fv[i,]<-out[[i]]$coxph_out[4,]
  coxph_fvl[i,]<-out[[i]]$coxph_out[5,]
  
  surv_x1[i,]<-out[[i]]$surv_out[1,]
  surv_x2[i,]<-out[[i]]$surv_out[2,]
  surv_x3[i,]<-out[[i]]$surv_out[3,]
  surv_fv[i,]<-out[[i]]$surv_out[4,]
  surv_fvl[i,]<-out[[i]]$surv_out[5,]
  
  HL_x1[i,]<-out[[i]]$HL_out[1,]
  HL_x2[i,]<-out[[i]]$HL_out[2,]
  HL_x3[i,]<-out[[i]]$HL_out[3,]
  HL_fv[i,]<-out[[i]]$HL_out[4,]
  HL_fvl[i,]<-out[[i]]$HL_out[5,]
  
}
coxph_x3[!complete.cases(coxph_x3), ] <- NA
coxph_x1[which(is.na(coxph_x3[[1]])),] <- NA
coxph_x2[which(is.na(coxph_x3[[1]])),] <- NA
coxph_x2[!complete.cases(coxph_x2), ] <- NA
coxph_x1[which(is.na(coxph_x2[[1]])),] <- NA
coxph_x3[which(is.na(coxph_x2[[1]])),] <- NA
coxph_x1[!complete.cases(coxph_x1), ] <- NA
coxph_x2[which(is.na(coxph_x1[[1]])),] <- NA
coxph_x3[which(is.na(coxph_x1[[1]])),] <- NA

EM_fv[!complete.cases(EM_fv),] <- NA
EM_x1[which(is.na(EM_fv[[1]])),] <- NA
EM_x2[which(is.na(EM_fv[[1]])),] <- NA
EM_x3[which(is.na(EM_fv[[1]])),] <- NA
EM_fvl[which(is.na(EM_fv[[1]])),] <- NA
EM_x1[!complete.cases(EM_x1),] <- NA
EM_x2[which(is.na(EM_x1[[1]])),] <- NA
EM_x3[which(is.na(EM_x1[[1]])),] <- NA
EM_fv[which(is.na(EM_x1[[1]])),] <- NA
EM_fvl[which(is.na(EM_x1[[1]])),] <- NA
EM_x2[!complete.cases(EM_x2),] <- NA
EM_x1[which(is.na(EM_x2[[1]])),] <- NA
EM_x3[which(is.na(EM_x2[[1]])),] <- NA
EM_fv[which(is.na(EM_x2[[1]])),] <- NA
EM_fvl[which(is.na(EM_x2[[1]])),] <- NA
EM_x3[!complete.cases(EM_x3),] <- NA
EM_x1[which(is.na(EM_x3[[1]])),] <- NA
EM_x2[which(is.na(EM_x3[[1]])),] <- NA
EM_fv[which(is.na(EM_x3[[1]])),] <- NA
EM_fvl[which(is.na(EM_x3[[1]])),] <- NA

HL_fv[!complete.cases(HL_fv), ] <- NA
HL_x1[which(is.na(HL_fv[[1]])),] <- NA
HL_x2[which(is.na(HL_fv[[1]])),] <- NA
HL_x3[which(is.na(HL_fv[[1]])),] <- NA
HL_fvl[which(is.na(HL_fv[[1]])),] <- NA
HL_x1[!complete.cases(HL_x1), ] <- NA
HL_x2[which(is.na(HL_x1[[1]])),] <- NA
HL_x3[which(is.na(HL_x1[[1]])),] <- NA
HL_fv[which(is.na(HL_x1[[1]])),] <- NA
HL_fvl[which(is.na(HL_x1[[1]])),] <- NA
HL_x2[!complete.cases(HL_x2), ] <- NA
HL_x1[which(is.na(HL_x2[[1]])),] <- NA
HL_x3[which(is.na(HL_x2[[1]])),] <- NA
HL_fv[which(is.na(HL_x2[[1]])),] <- NA
HL_fvl[which(is.na(HL_x2[[1]])),] <- NA
HL_x3[!complete.cases(HL_x3), ] <- NA
HL_x1[which(is.na(HL_x3[[1]])),] <- NA
HL_x2[which(is.na(HL_x3[[1]])),] <- NA
HL_fv[which(is.na(HL_x3[[1]])),] <- NA
HL_fvl[which(is.na(HL_x3[[1]])),] <- NA

surv_fv[!complete.cases(surv_fv), ] <- NA
surv_x1[which(is.na(surv_fv[[1]])),] <- NA
surv_x2[which(is.na(surv_fv[[1]])),] <- NA
surv_x3[which(is.na(surv_fv[[1]])),] <- NA
surv_fvl[which(is.na(surv_fv[[1]])),] <- NA
surv_x1[!complete.cases(surv_x1), ] <- NA
surv_x2[which(is.na(surv_x1[[1]])),] <- NA
surv_x3[which(is.na(surv_x1[[1]])),] <- NA
surv_fv[which(is.na(surv_x1[[1]])),] <- NA
surv_fvl[which(is.na(surv_x1[[1]])),] <- NA
surv_x2[!complete.cases(surv_x2), ] <- NA
surv_x1[which(is.na(surv_x2[[1]])),] <- NA
surv_x3[which(is.na(surv_x2[[1]])),] <- NA
surv_fv[which(is.na(surv_x2[[1]])),] <- NA
surv_fvl[which(is.na(surv_x2[[1]])),] <- NA
surv_x3[!complete.cases(surv_x3), ] <- NA
surv_x1[which(is.na(surv_x3[[1]])),] <- NA
surv_x2[which(is.na(surv_x3[[1]])),] <- NA
surv_fv[which(is.na(surv_x3[[1]])),] <- NA
surv_fvl[which(is.na(surv_x3[[1]])),] <- NA


Par_fv[!complete.cases(Par_fv), ] <- NA
Par_x1[which(is.na(Par_fv[[1]])),] <- NA
Par_x2[which(is.na(Par_fv[[1]])),] <- NA
Par_x3[which(is.na(Par_fv[[1]])),] <- NA
Par_fvl[which(is.na(Par_fv[[1]])),] <- NA
Par_x1[!complete.cases(Par_x1), ] <- NA
Par_x2[which(is.na(Par_x1[[1]])),] <- NA
Par_x3[which(is.na(Par_x1[[1]])),] <- NA
Par_fv[which(is.na(Par_x1[[1]])),] <- NA
Par_fvl[which(is.na(Par_x1[[1]])),] <- NA
Par_x2[!complete.cases(Par_x2), ] <- NA
Par_x1[which(is.na(Par_x2[[1]])),] <- NA
Par_x3[which(is.na(Par_x2[[1]])),] <- NA
Par_fv[which(is.na(Par_x2[[1]])),] <- NA
Par_fvl[which(is.na(Par_x2[[1]])),] <- NA
Par_x3[!complete.cases(Par_x3), ] <- NA
Par_x1[which(is.na(Par_x3[[1]])),] <- NA
Par_x2[which(is.na(Par_x3[[1]])),] <- NA
Par_fv[which(is.na(Par_x3[[1]])),] <- NA
Par_fvl[which(is.na(Par_x3[[1]])),] <- NA

fpack_fv[!complete.cases(fpack_fv), ] <- NA
fpack_x1[which(is.na(fpack_fv[[1]])),] <- NA
fpack_x2[which(is.na(fpack_fv[[1]])),] <- NA
fpack_x3[which(is.na(fpack_fv[[1]])),] <- NA
fpack_fvl[which(is.na(fpack_fv[[1]])),] <- NA
fpack_x1[!complete.cases(fpack_x1), ] <- NA
fpack_x2[which(is.na(fpack_x1[[1]])),] <- NA
fpack_x3[which(is.na(fpack_x1[[1]])),] <- NA
fpack_fv[which(is.na(fpack_x1[[1]])),] <- NA
fpack_fvl[which(is.na(fpack_x1[[1]])),] <- NA
fpack_x2[!complete.cases(fpack_x2), ] <- NA
fpack_x1[which(is.na(fpack_x2[[1]])),] <- NA
fpack_x3[which(is.na(fpack_x2[[1]])),] <- NA
fpack_fv[which(is.na(fpack_x2[[1]])),] <- NA
fpack_fvl[which(is.na(fpack_x2[[1]])),] <- NA
fpack_x3[!complete.cases(fpack_x3), ] <- NA
fpack_x1[which(is.na(fpack_x3[[1]])),] <- NA
fpack_x2[which(is.na(fpack_x3[[1]])),] <- NA
fpack_fv[which(is.na(fpack_x3[[1]])),] <- NA
fpack_fvl[which(is.na(fpack_x3[[1]])),] <- NA

converged.coxph_x1 = as.numeric(!is.na(coxph_x1$coef) & !is.na(coxph_x1$se))
(sum(converged.coxph_x1)/n_sims)*100
converged.coxph_x2 = as.numeric(!is.na(coxph_x2$coef) & !is.na(coxph_x2$se))
(sum(converged.coxph_x2)/n_sims)*100
converged.coxph_x3 = as.numeric(!is.na(coxph_x3$coef) & !is.na(coxph_x3$se))
(sum(converged.coxph_x3)/n_sims)*100
converged.coxph_fv = NA

converged.Par_x1 = as.numeric(!is.na(Par_x1$coef) & !is.na(Par_x1$se))
(sum(converged.Par_x1)/n_sims)*100
converged.Par_x2 = as.numeric(!is.na(Par_x2$coef) & !is.na(Par_x2$se))
(sum(converged.Par_x2)/n_sims)*100
converged.Par_x3 = as.numeric(!is.na(Par_x3$coef) & !is.na(Par_x3$se))
(sum(converged.Par_x3)/n_sims)*100
converged.Par_fv = as.numeric(!is.na(Par_fv$coef) & !is.na(Par_fv$se))
(sum(converged.Par_fv)/n_sims)*100

converged.EM_x1 = as.numeric(!is.na(EM_x1$coef) & !is.na(EM_x1$se))
(sum(converged.EM_x1)/n_sims)*100
converged.EM_x2 = as.numeric(!is.na(EM_x2$coef) & !is.na(EM_x2$se))
(sum(converged.EM_x2)/n_sims)*100
converged.EM_x3 = as.numeric(!is.na(EM_x3$coef) & !is.na(EM_x3$se))
(sum(converged.EM_x3)/n_sims)*100
converged.EM_fv = as.numeric(!is.na(EM_fv$coef) & !is.na(EM_fv$se))
(sum(converged.EM_fv)/n_sims)*100

converged.surv_x1 = as.numeric(!is.na(surv_x1$coef) & !is.na(surv_x1$se))
(sum(converged.surv_x1)/n_sims)*100
converged.surv_x2 = as.numeric(!is.na(surv_x2$coef) & !is.na(surv_x2$se))
(sum(converged.surv_x2)/n_sims)*100
converged.surv_x3 = as.numeric(!is.na(surv_x3$coef) & !is.na(surv_x3$se))
(sum(converged.surv_x3)/n_sims)*100
converged.surv_fv = as.numeric(!is.na(surv_fv$coef) & !is.na(surv_fv$se))
(sum(converged.surv_fv)/n_sims)*100

converged.HL_x1 = as.numeric(!is.na(HL_x1$coef) & !is.na(HL_x1$se))
(sum(converged.HL_x1)/n_sims)*100
converged.HL_x2 = as.numeric(!is.na(HL_x2$coef) & !is.na(HL_x2$se))
(sum(converged.HL_x2)/n_sims)*100
converged.HL_x3 = as.numeric(!is.na(HL_x3$coef) & !is.na(HL_x3$se))
(sum(converged.HL_x3)/n_sims)*100
converged.HL_fv = as.numeric(!is.na(HL_fv$coef) & !is.na(HL_fv$se))
(sum(converged.HL_fv)/n_sims)*100

converged.fpack_x1 = as.numeric(!is.na(fpack_x1$coef) & !is.na(fpack_x1$se))
(sum(converged.fpack_x1)/n_sims)*100
converged.fpack_x2 = as.numeric(!is.na(fpack_x2$coef) & !is.na(fpack_x2$se))
(sum(converged.fpack_x2)/n_sims)*100
converged.fpack_x3 = as.numeric(!is.na(fpack_x3$coef) & !is.na(fpack_x3$se))
(sum(converged.fpack_x3)/n_sims)*100
converged.fpack_fv = as.numeric(!is.na(fpack_fv$coef) & !is.na(fpack_fv$se))
(sum(converged.fpack_fv)/n_sims)*100

Converged<-c((sum(converged.coxph_x1)/n_sims)*100, (sum(converged.coxph_x2)/n_sims)*100,
             (sum(converged.coxph_x3)/n_sims)*100, NA,
             (sum(converged.Par_x1)/n_sims)*100,(sum(converged.Par_x2)/n_sims)*100, 
             (sum(converged.Par_x3)/n_sims)*100,(sum(converged.Par_fv)/n_sims)*100, 
             (sum(converged.EM_x1)/n_sims)*100,(sum(converged.EM_x2)/n_sims)*100, 
             (sum(converged.EM_x3)/n_sims)*100,(sum(converged.EM_fv)/n_sims)*100, 
             (sum(converged.surv_x1)/n_sims)*100,(sum(converged.surv_x2)/n_sims)*100,
             (sum(converged.surv_x3)/n_sims)*100,(sum(converged.surv_fv)/n_sims)*100, 
             (sum(converged.HL_x1)/n_sims)*100,(sum(converged.HL_x2)/n_sims)*100,
             (sum(converged.HL_x3)/n_sims)*100,(sum(converged.HL_fv)/n_sims)*100, 
             (sum(converged.fpack_x1)/n_sims)*100,(sum(converged.fpack_x2)/n_sims)*100,
             (sum(converged.fpack_x3)/n_sims)*100,(sum(converged.fpack_fv)/n_sims)*100)



empse1_coxph<-sqrt(1/(n_sims-1)*sum((coxph_x1$coef-
                                       mean(coxph_x1$coef, na.rm=TRUE))^2,na.rm=TRUE))
empse2_coxph<-sqrt(1/(n_sims-1)*sum((coxph_x2$coef-
                                       mean(coxph_x2$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse3_coxph<-sqrt(1/(n_sims-1)*sum((coxph_x3$coef-
                                       mean(coxph_x3$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse_fv_coxph<-sqrt(1/(n_sims-1)*sum((coxph_fv$coef-
                                         mean(coxph_fv$coef,na.rm=TRUE))^2,na.rm=TRUE))

empse1_par<-sqrt(1/(n_sims-1)*sum((Par_x1$coef-
                                     mean(Par_x1$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse2_par<-sqrt(1/(n_sims-1)*sum((Par_x2$coef-
                                     mean(Par_x2$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse3_par<-sqrt(1/(n_sims-1)*sum((Par_x3$coef-
                                     mean(Par_x3$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse_fv_par<-sqrt(1/(n_sims-1)*sum((Par_fv$coef-
                                       mean(Par_fv$coef,na.rm=TRUE))^2,na.rm=TRUE))

empse1_em<-sqrt(1/(n_sims-1)*sum((EM_x1$coef-mean(EM_x1$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse2_em<-sqrt(1/(n_sims-1)*sum((EM_x2$coef-mean(EM_x2$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse3_em<-sqrt(1/(n_sims-1)*sum((EM_x3$coef-mean(EM_x3$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse_fv_em<-sqrt(1/(n_sims-1)*sum((EM_fv$coef-mean(EM_fv$coef,na.rm=TRUE))^2,na.rm=TRUE))

empse1_surv<-sqrt(1/(n_sims-1)*sum((surv_x1$coef-
                                      mean(surv_x1$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse2_surv<-sqrt(1/(n_sims-1)*sum((surv_x2$coef-
                                      mean(surv_x2$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse3_surv<-sqrt(1/(n_sims-1)*sum((surv_x3$coef-
                                      mean(surv_x3$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse_fv_surv<-sqrt(1/(n_sims-1)*sum((surv_fv$coef-
                                        mean(surv_fv$coef,na.rm=TRUE))^2,na.rm=TRUE))

empse1_HL<-sqrt(1/(n_sims-1)*sum((HL_x1$coef-mean(HL_x1$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse2_HL<-sqrt(1/(n_sims-1)*sum((HL_x2$coef-mean(HL_x2$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse3_HL<-sqrt(1/(n_sims-1)*sum((HL_x3$coef-mean(HL_x3$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse_fv_HL<-sqrt(1/(n_sims-1)*sum((HL_fv$coef-mean(HL_fv$coef,na.rm=TRUE))^2,na.rm=TRUE))

empse1_fpack<-sqrt(1/(n_sims-1)*sum((fpack_x1$coef-
                                       mean(fpack_x1$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse2_fpack<-sqrt(1/(n_sims-1)*sum((fpack_x2$coef-
                                       mean(fpack_x2$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse3_fpack<-sqrt(1/(n_sims-1)*sum((fpack_x3$coef-
                                       mean(fpack_x3$coef,na.rm=TRUE))^2,na.rm=TRUE))
empse_fv_fpack<-sqrt(1/(n_sims-1)*sum((fpack_fv$coef-
                                         mean(fpack_fv$coef,na.rm=TRUE))^2,na.rm=TRUE))


Empirical_SE<-c(empse1_coxph,empse2_coxph,empse3_coxph,empse_fv_coxph,
                empse1_par,empse2_par,empse3_par,empse_fv_par,
                empse1_em,empse2_em,empse3_em,empse_fv_em,
                empse1_surv,empse2_surv,empse3_surv,empse_fv_surv,
                empse1_HL,empse2_HL,empse3_HL,empse_fv_HL,
                empse1_fpack,empse2_fpack,empse3_fpack,empse_fv_fpack)


cpx1_coxph<-mean(coxph_x1$CI_LB<= x1 & x1 <= coxph_x1$CI_UB,na.rm=TRUE)*100
cpx2_coxph<-mean(coxph_x2$CI_LB<= x2 & x2 <= coxph_x2$CI_UB,na.rm=TRUE)*100
cpx3_coxph<-mean(coxph_x3$CI_LB<= x3 & x3 <= coxph_x3$CI_UB,na.rm=TRUE)*100
cpfv_coxph<-NA

cpx1_par<-mean(Par_x1$CI_LB<= x1 & x1 <= Par_x1$CI_UB,na.rm=TRUE)*100
cpx2_par<-mean(Par_x2$CI_LB<= x2 & x2 <= Par_x2$CI_UB,na.rm=TRUE)*100
cpx3_par<-mean(Par_x3$CI_LB<= x3 & x3 <= Par_x3$CI_UB,na.rm=TRUE)*100
cpfv_par<-mean(Par_fv$CI_LB<= fv & fv <= Par_fv$CI_UB,na.rm=TRUE)*100

cpx1_em<-mean(EM_x1$CI_LB<= x1 & x1 <= EM_x1$CI_UB,na.rm=TRUE)*100
cpx2_em<-mean(EM_x2$CI_LB<= x2 & x2 <= EM_x2$CI_UB,na.rm=TRUE)*100
cpx3_em<-mean(EM_x3$CI_LB<= x3 & x3 <= EM_x3$CI_UB,na.rm=TRUE)*100
cpfv_em<-mean(EM_fv$CI_LB<= fv & fv <= EM_fv$CI_UB,na.rm=TRUE)*100

cpx1_surv<-mean(surv_x1$CI_LB<= x1 & x1 <= surv_x1$CI_UB,na.rm=TRUE)*100
cpx2_surv<-mean(surv_x2$CI_LB<= x2 & x2 <= surv_x2$CI_UB,na.rm=TRUE)*100
cpx3_surv<-mean(surv_x3$CI_LB<= x3 & x3 <= surv_x3$CI_UB,na.rm=TRUE)*100
cpfv_surv<-mean(surv_fv$CI_LB<= fv & fv <= surv_fv$CI_UB,na.rm=TRUE)*100

cpx1_HL<-mean(HL_x1$CI_LB<= x1 & x1 <= HL_x1$CI_UB,na.rm=TRUE)*100
cpx2_HL<-mean(HL_x2$CI_LB<= x2 & x2 <= HL_x2$CI_UB,na.rm=TRUE)*100
cpx3_HL<-mean(HL_x3$CI_LB<= x3 & x3 <= HL_x3$CI_UB,na.rm=TRUE)*100
cpfv_HL<-mean(HL_fv$CI_LB<= fv & fv <= HL_fv$CI_UB,na.rm=TRUE)*100

cpx1_fpack<-mean(fpack_x1$CI_LB<= x1 & x1 <= fpack_x1$CI_UB,na.rm=TRUE)*100
cpx2_fpack<-mean(fpack_x2$CI_LB<= x2 & x2 <= fpack_x2$CI_UB,na.rm=TRUE)*100
cpx3_fpack<-mean(fpack_x3$CI_LB<= x3 & x3 <= fpack_x3$CI_UB,na.rm=TRUE)*100
cpfv_fpack<-mean(fpack_fv$CI_LB<= fv & fv <= fpack_fv$CI_UB,na.rm=TRUE)*100

CP<-c(cpx1_coxph,cpx2_coxph,cpx3_coxph,cpfv_coxph,
      cpx1_par,cpx2_par,cpx3_par,cpfv_par,
      cpx1_em,cpx2_em,cpx3_em,cpfv_em,
      cpx1_surv,cpx2_surv,cpx3_surv,cpfv_surv,
      cpx1_HL,cpx2_HL,cpx3_HL,cpfv_HL,
      cpx1_fpack,cpx2_fpack,cpx3_fpack,cpfv_fpack)

cpx1_coxph_log<-NA
cpx2_coxph_log<-NA
cpx3_coxph_log<-NA
cpfv_coxph_log<-NA

cpx1_par_log<-NA
cpx2_par_log<-NA
cpx3_par_log<-NA
cpfv_par_log<-mean(Par_fvl$CI_LB<= fv & fv <= Par_fvl$CI_UB,na.rm=TRUE)*100

cpx1_em_log<-NA
cpx2_em_log<-NA
cpx3_em_log<-NA
cpfv_em_log<-mean(EM_fvl$CI_LB<= fv & fv <= EM_fvl$CI_UB,na.rm=TRUE)*100

cpx1_surv_log<-NA
cpx2_surv_log<-NA
cpx3_surv_log<-NA
cpfv_surv_log<-mean(surv_fvl$CI_LB<= fv & fv <= surv_fvl$CI_UB,na.rm=TRUE)*100

cpx1_HL_log<-NA
cpx2_HL_log<-NA
cpx3_HL_log<-NA
cpfv_HL_log<-mean(HL_fvl$CI_LB<= fv & fv <= HL_fvl$CI_UB,na.rm=TRUE)*100

cpx1_fpack_log<-NA
cpx2_fpack_log<-NA
cpx3_fpack_log<-NA
cpfv_fpack_log<-mean(fpack_fvl$CI_LB<= fv & fv <= fpack_fvl$CI_UB,na.rm=TRUE)*100

CP_logTheta<-c(cpx1_coxph_log,cpx2_coxph_log,cpx3_coxph_log,cpfv_coxph_log,
               cpx1_par_log,cpx2_par_log,cpx3_par_log,cpfv_par_log,
               cpx1_em_log,cpx2_em_log,cpx3_em_log,cpfv_em_log,
               cpx1_surv_log,cpx2_surv_log,cpx3_surv_log,cpfv_surv_log,
               cpx1_HL_log,cpx2_HL_log,cpx3_HL_log,cpfv_HL_log,
               cpx1_fpack_log,cpx2_fpack_log,cpx3_fpack_log,cpfv_fpack_log)


True<-rep(c(x1,x2,x3,fv),times=6)
Mean_coef<-c(mean(coxph_x1$coef, na.rm=TRUE),mean(coxph_x2$coef, na.rm=TRUE),
             mean(coxph_x3$coef, na.rm=TRUE), mean(coxph_fv$coef,na.rm=TRUE),
             mean(Par_x1$coef,na.rm=TRUE),mean(Par_x2$coef,na.rm=TRUE),
             mean(Par_x3$coef,na.rm=TRUE),mean(Par_fv$coef,na.rm=TRUE),
             mean(EM_x1$coef,na.rm=TRUE),mean(EM_x2$coef,na.rm=TRUE),
             mean(EM_x3$coef,na.rm=TRUE),mean(EM_fv$coef,na.rm=TRUE),
             mean(surv_x1$coef,na.rm=TRUE),mean(surv_x2$coef,na.rm=TRUE),
             mean(surv_x3$coef,na.rm=TRUE), mean(surv_fv$coef,na.rm=TRUE),
             mean(HL_x1$coef,na.rm=TRUE),mean(HL_x2$coef,na.rm=TRUE),
             mean(HL_x3$coef,na.rm=TRUE),mean(HL_fv$coef,na.rm=TRUE), 
             mean(fpack_x1$coef,na.rm=TRUE),mean(fpack_x2$coef,na.rm=TRUE),
             mean(fpack_x3$coef,na.rm=TRUE),mean(fpack_fv$coef,na.rm=TRUE))

Mean_se<-c(mean(coxph_x1$se, na.rm=TRUE),mean(coxph_x2$se, na.rm=TRUE),
           mean(coxph_x3$se, na.rm=TRUE), mean(coxph_fv$se,na.rm=TRUE),
           mean(Par_x1$se,na.rm=TRUE),mean(Par_x2$se,na.rm=TRUE),
           mean(Par_x3$se,na.rm=TRUE),mean(Par_fv$se,na.rm=TRUE),
           mean(EM_x1$se,na.rm=TRUE),mean(EM_x2$se,na.rm=TRUE),
           mean(EM_x3$se,na.rm=TRUE),mean(EM_fv$se,na.rm=TRUE),
           mean(surv_x1$se,na.rm=TRUE),mean(surv_x2$se,na.rm=TRUE),
           mean(surv_x3$se,na.rm=TRUE), mean(surv_fv$se,na.rm=TRUE),
           mean(HL_x1$se,na.rm=TRUE),mean(HL_x2$se,na.rm=TRUE),
           mean(HL_x3$se,na.rm=TRUE),mean(HL_fv$se,na.rm=TRUE),
           mean(fpack_x1$se,na.rm=TRUE),mean(fpack_x2$se,na.rm=TRUE),
           mean(fpack_x3$se,na.rm=TRUE),mean(fpack_fv$se,na.rm=TRUE))

Bias<- Mean_coef-True

Median_coef<-c(median(coxph_x1$coef, na.rm=TRUE),median(coxph_x2$coef, na.rm=TRUE),
               median(coxph_x3$coef, na.rm=TRUE), median(coxph_fv$coef,na.rm=TRUE),
               median(Par_x1$coef,na.rm=TRUE),median(Par_x2$coef,na.rm=TRUE),
               median(Par_x3$coef,na.rm=TRUE),median(Par_fv$coef,na.rm=TRUE),
               median(EM_x1$coef,na.rm=TRUE),median(EM_x2$coef,na.rm=TRUE),
               median(EM_x3$coef,na.rm=TRUE),median(EM_fv$coef,na.rm=TRUE),
               median(surv_x1$coef,na.rm=TRUE),median(surv_x2$coef,na.rm=TRUE),
               median(surv_x3$coef,na.rm=TRUE), median(surv_fv$coef,na.rm=TRUE),
               median(HL_x1$coef,na.rm=TRUE),median(HL_x2$coef,na.rm=TRUE),
               median(HL_x3$coef,na.rm=TRUE),median(HL_fv$coef,na.rm=TRUE),
               median(fpack_x1$coef,na.rm=TRUE),median(fpack_x2$coef,na.rm=TRUE),
               median(fpack_x3$coef,na.rm=TRUE),median(fpack_fv$coef,na.rm=TRUE))

Median_se<-c(median(coxph_x1$se, na.rm=TRUE),median(coxph_x2$se, na.rm=TRUE),
             median(coxph_x3$se, na.rm=TRUE), median(coxph_fv$se,na.rm=TRUE),
             median(Par_x1$se,na.rm=TRUE),median(Par_x2$se,na.rm=TRUE),
             median(Par_x3$se,na.rm=TRUE),median(Par_fv$se,na.rm=TRUE),
             median(EM_x1$se,na.rm=TRUE),median(EM_x2$se,na.rm=TRUE),
             median(EM_x3$se,na.rm=TRUE),median(EM_fv$se,na.rm=TRUE),
             median(surv_x1$se,na.rm=TRUE),median(surv_x2$se,na.rm=TRUE),
             median(surv_x3$se,na.rm=TRUE), median(surv_fv$se,na.rm=TRUE),
             median(HL_x1$se,na.rm=TRUE),median(HL_x2$se,na.rm=TRUE),
             median(HL_x3$se,na.rm=TRUE),median(HL_fv$se,na.rm=TRUE),
             median(fpack_x1$se,na.rm=TRUE),median(fpack_x2$se,na.rm=TRUE),
             median(fpack_x3$se,na.rm=TRUE),median(fpack_fv$se,na.rm=TRUE))

time<-c(mean(coxph_x1$time,),mean(coxph_x2$time),
        mean(coxph_x3$time), mean(coxph_fv$time),
        mean(Par_x1$time),mean(Par_x2$time),
        mean(Par_x3$time),mean(Par_fv$time),
        mean(EM_x1$time,na.rm=TRUE),mean(EM_x2$time,na.rm=TRUE),
        mean(EM_x3$time,na.rm=TRUE),mean(EM_fv$time,na.rm=TRUE),
        mean(surv_x1$time,na.rm=TRUE),mean(surv_x2$time,na.rm=TRUE),
        mean(surv_x3$time,na.rm=TRUE), mean(surv_fv$time,na.rm=TRUE),
        mean(HL_x1$time,na.rm=TRUE),mean(HL_x2$time,na.rm=TRUE),
        mean(HL_x3$time,na.rm=TRUE),mean(HL_fv$time,na.rm=TRUE), 
        mean(fpack_x1$time),mean(fpack_x2$time),
        mean(fpack_x3$time),mean(fpack_fv$time))

df<-data.frame(True,round(Mean_coef,3),round(Bias,3),round(Mean_se,3),
               round(Empirical_SE,3),round(Median_coef,3), 
               round(Median_se,3),round(CP,2),
               round(CP_logTheta,2),Converged,round(time,5))

row.names(df)<-c("coxph.x1","coxph.x2","coxph.x3","coxph.fv",
                 "Par.x1","Par.x2","Par.x3","Par.fv",
                 "EM.x1","EM.x2","EM.x3","EM.fv",
                 "surv.x1","surv.x2","surv.x3","surv.fv",
                 "HL.x1","HL.x2","HL.x3","HL.fv",
                 "fpack.x1","fpack.x2","fpack.x3","fpack.fv")
colnames(df)<-c("True","Mean_coef","Bias","Mean_SE","Empirical_SE",
                "Median_coef","Median_se","CP","CP_logTheta",
                "Converged","time")



