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

kfold_fn<-function (y, k, list = TRUE, returnTrain = FALSE) 
{
  if (class(y)[1] == "Surv") 
    y <- y[, "time"]
  if (is.numeric(y)) y<-as.factor(y)
  # {
  #   cuts <- floor(length(y)/k)
  #   if (cuts < 2) 
  #     cuts <- 2
  #   if (cuts > 5) 
  #     cuts <- 5
  #   breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
  #   y <- cut(y, breaks, include.lowest = TRUE)
  # }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0) 
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <-
          sample(1:k,size = numInClass[i])
      }
    }
    fold_num<-length(table(foldVector))
    if (fold_num !=k){
      subsets <- replicate(1,sample(length(y)))
      fold <- rep(seq_len(k), length.out = length(y))
      foldVector<-fold[order(subsets)]
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
                        sep = "")
    if (returnTrain) 
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}


make_fold<-function(fix_var,y,k,censor)
{
  ########Find out which column is factor data
  col<-as.logical(rep(0,ncol(fix_var)))
  for(i in 1:ncol(fix_var)){
    col[i]<-is.factor(fix_var[,i])
  }
  categ_col<-which(col==TRUE)
  
  #######combine fix variable and censor column#####
  fix_censor<-data.frame(fix_var,censor)
  fix_censor_col<-ncol (fix_censor)
  
  #####make a fold list
  if(length(categ_col)==0) fold_list<-kfold_fn(y,k)
  if(length(categ_col)!=0){
    counter <- 0
    allfoldisgood<-FALSE
    while(isFALSE(allfoldisgood)){
      fold_list<-kfold_fn(y,k)
      k2<-length(fold_list)
      if(k2 != k) warning(paste("It cannot get", k, "folds,the fold number is",k2))
      test_data <- as.list(1:k2)
      train_data<- as.list(1:k2)
      foldisgood<-matrix(FALSE,nrow =length(categ_col),ncol =k2)
      
      train_data_censor<-as.list(1:k2)
      check_table<-matrix(FALSE,nrow =length(categ_col),ncol =k2)
      
      for(i in 1:length(categ_col)){
        for(fid in 1:length (fold_list)) {
          #check factor column in test data within train data
          test_data[[fid]] <- fix_var[fold_list[[fid]], ,drop=FALSE]
          train_data[[fid]]<- fix_var[-fold_list[[fid]], ,drop=FALSE]
          foldisgood[i,fid]<- all(test_data[[fid]][,categ_col[i]]%in% 
                                    train_data[[fid]][,categ_col[i]])
          
          #check cross table contain zero or not in r (factor cross with censor)
          train_data_censor[[fid]] <- fix_censor[-fold_list[[fid]], ,drop=FALSE]
          check_table[i,fid]<-all(apply(table(train_data_censor[[fid]][,categ_col[i]],
                                              train_data_censor[[fid]][,fix_censor_col]),
                                        2,function(x) x!=0))
        }
      }
      
      allfoldisgood<-all(foldisgood,check_table)
      counter <- counter+1
      if (counter>500) stop("It cannot get the fold list ")
    }
  }
  fold_list
}

cv_zresidual.coxph<- function( fit, data, nfolds=NULL,foldlist=NULL)
{
  # Required packages:
  if (!requireNamespace("pacman")) {
    install.packages("pacman")
  }
  pacman::p_load(
    "foreach",
    "parallel",
    "doParallel",
    "stringr"
  )
  
  mf <- model.frame(fit$formula, data)
  nc <- ncol (mf)
  fix_var<-mf[,-c(1,nc),drop=FALSE]
  
  if(is.null(nfolds))
  {
    ncores <- detectCores() 
    cl <- makeForkCluster(ncores)
    registerDoParallel(cl)
    nfolds<-10 %/% ncores*ncores
  }
  
  if(is.null(foldlist)) foldlist<-make_fold(fix_var=fix_var,y=mf[,nc],
                                            k=nfolds,censor=mf[,1][,2])
  
  res_fold <- as.list(1:nfolds)
  for(fid in 1:length (foldlist)) 
  {
    testcases <- foldlist[[fid]]
    data.train <- data[-testcases, ]
    data.test <- data[testcases, ]
    fit_traindata <- tryCatch(
      coxph(fit$formula,data=data.train),
      error = function(e) NA,
      warning = function(w) NA)
    
    if(any(is.na(fit_traindata))) {
      res_fold[[fid]]<-rep(NA,length(foldlist[[fid]]))
      attr(res_fold[[fid]], "SP") <-rep(NA,length(foldlist[[fid]]))
      attr(res_fold[[fid]], "hazfn") <-rep(NA,length(foldlist[[fid]]))
    }
    
    if(any(!is.na(fit_traindata))){
      res_fold[[fid]]<-zresidual.coxph(fit_coxph=fit_traindata,
                                       traindata=data.train,
                                       newdata = data.test)
    }
    
  }
  
  res_cv <- rep (0, nrow(data))
  attr(res_cv, "SP") <- rep(0, nrow(data))
  for(fid in 1:length (foldlist)) {
    testcases <- foldlist[[fid]]
    res_cv[testcases]<-res_fold[[fid]]
    attr(res_cv, "SP")[testcases] <- attr(res_fold[[fid]],"SP")
  }
  res_cv
}



### Simulate nonlinear dataset ##################################
#set.seed(1)
rexp2 <- function(n, rate){ if (rate==0) rep(Inf,n) else rexp(n = n, rate = rate)}
simulWeib <- function(n_clusters,n_individuals, lambda, alpha,beta1,beta2,beta3,
                      mean.censor,fv)
{
  # Make covariates 
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  grpid<- paste0("g",grpid)
  x1 = runif(n_clusters*n_individuals,0,1)
  x2 = abs(round(rnorm(n_clusters*n_individuals, 0,1),2))
  x3 = rbinom(n_clusters*n_individuals, size = 1, p = 0.5)
  
  # Make frailty (shape = a and scale = s,E(X) = a*s and Var(X) = a*s^2)
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

### Simulate linear dataset#####################
#set.seed(1)
rexp2 <- function(n, rate){ if (rate==0) rep(Inf,n) else rexp(n = n, rate = rate)}
simulWeib <- function(n_clusters,n_individuals, lambda, alpha,beta1,beta2,beta3,
                      mean.censor,fv)
{
  # Make covariates 
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  grpid<- paste0("g",grpid)
  x1 = runif(n_clusters*n_individuals,0,1)
  x2 = round(rnorm(n_clusters*n_individuals, 0,1),2)
  x3 = rbinom(n_clusters*n_individuals, size = 1, p = 0.5)
  
  # Make frailty (shape = a and scale = s,E(X) = a*s and Var(X) = a*s^2)
  frvec<-rep(rgamma(n_clusters,shape = 1/fv,scale = fv), each = n_individuals)
  
  #Event times:baseline hazard is Weibull，lambda=scale parameter, alpha=shape parameter
  n<-n_individuals * n_clusters
  v <- runif(n)
  surv_time <-
    (- log(v)/(lambda*exp(x1*beta1+x2*beta2+x3*beta3)*frvec))^(1 /alpha)
  
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


### Simulate data with 10% outliers

simulWeib_outlier1 <- function(n_clusters,n_individuals, lambda, alpha,
                               beta1,beta2,beta3, mean.censor,fv)
{
  # Make covariates
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  grpid<- paste0("g",grpid)
  x1 = runif(n_clusters*n_individuals,0,1)
  x2 = round(rnorm(n_clusters*n_individuals, 0,1),2)
  x3 = rbinom(n_clusters*n_individuals, size = 1, p = 0.5)
  
  # Make frailty (shape = a and scale = s,E(X) = a*s and Var(X) = a*s^2)
  #frvec<-rep(rgamma(n_clusters,shape=0.5, rate=1), each = n_individuals)
  frvec<-rep(rgamma(n_clusters,shape = 1/fv,scale = fv), each = n_individuals)
  
  #Event times:baseline hazard is Weibull，lambda=scale parameter, alpha=shape parameter
  n<-n_individuals * n_clusters
  v <- runif(n)
  surv_time <-
    (- log(v)/(lambda*exp(x1*beta1+x2*beta2+x3*beta3)*frvec))^(1 /alpha)
  
  # Censoring times
  t0<- rexp2(n, rate= 1/mean.censor)
  # follow-up times and event indicators
  t <- pmin(surv_time, t0)
  d <- as.numeric(surv_time <= t0)
  #make outliers around 10%
  outlier<-sample(which(d==1),size=round(length(which(d==1))*0.1))
  t[outlier]<- t[outlier]+ pmax(rexp(length(outlier),rate=1/2),2)
  outlier_indicator<-rep(0,n)
  outlier_indicator[outlier]<-1
  # Make data frame
  grpid<- as.factor(grpid)
  x3<-as.factor(x3)
  out<-data.frame(id,grpid,x1,x2,x3,t,d,frvec,outlier_indicator)
  return (out)
}

### Simulate data with 10 outliers
simulWeib_outlier4 <- function(n_clusters,n_individuals, lambda, alpha,beta1,beta2,beta3,
                               mean.censor,fv)
{
  # Make covariates
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  grpid<- paste0("g",grpid)
  x1 = runif(n_clusters*n_individuals,0,1)
  x2 = round(rnorm(n_clusters*n_individuals, 0,1),2)
  x3 = rbinom(n_clusters*n_individuals, size = 1, p = 0.5)
  
  # Make frailty (shape = a and scale = s,E(X) = a*s and Var(X) = a*s^2)
  #frvec<-rep(rgamma(n_clusters,shape=0.5, rate=1), each = n_individuals)
  frvec<-rep(rgamma(n_clusters,shape = 1/fv,scale = fv), each = n_individuals)
  
  #Event times:baseline hazard is Weibull，lambda=scale parameter, alpha=shape parameter
  n<-n_individuals * n_clusters
  v <- runif(n)
  surv_time <-
    (- log(v)/(lambda*exp(x1*beta1+x2*beta2+x3*beta3)*frvec))^(1 /alpha)
  # Censoring times
  t0<- rexp2(n, rate= 1/mean.censor)
  # follow-up times and event indicators
  t <- pmin(surv_time, t0)
  d <- as.numeric(surv_time <= t0)
  #make 10 outliers
  outlier<-sample(which(d==1),size=10)
  t[outlier]<- t[outlier]+ pmax(rexp(length(outlier),rate=1/2),2)
  outlier_indicator<-rep(0,n)
  outlier_indicator[outlier]<-1
  # Make data frame
  grpid<- as.factor(grpid)
  x3<-as.factor(x3)
  out<-data.frame(id,grpid,x1,x2,x3,t,d,frvec,outlier_indicator)
  return (out)
}

library("survival")
library("EnvStats")
library("nortest")
library("backports")
library("stringr")
library("pROC")

#####simulation study for detecting the non-linear covariate effect
para.frame <- data.frame(
  n_clusters = c(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
                 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10),
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
i=1
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
  
  ####fit true model##################################################### 
  fit_coxph <- tryCatch(
    coxph(Surv(t, d) ~ x1+log(x2)+x3 +frailty(grpid,distribution = "gamma"),
          data = simulated_shared_frailty_data),
    error = function(e) NA,
    warning = function(w) NA
  )
  
  ####fit wrong model##################################################### 
  #simulated_shared_frailty_data$ey5<-exp(simulated_shared_frailty_data$x2)
  fit_coxph_w <- tryCatch(
    coxph(Surv(t, d) ~ x1 + x2 + x3 +frailty(grpid, distribution="gamma"),
          data= simulated_shared_frailty_data),
    error = function(e) NA,
    warning = function(w) NA )
  
  
}
censorship<-as.numeric(table(simulated_shared_frailty_data$d)[1]/
                         nrow(simulated_shared_frailty_data))

true_qr<- true.qresidual(data=simulated_shared_frailty_data,
                         lambda=0.007,alpha=3,beta1=1,beta2=-2,beta3=0.5)

coxph_qr<-zresidual.coxph (fit_coxph = fit_coxph,
                           traindata = simulated_shared_frailty_data,
                           newdata = simulated_shared_frailty_data)
coxph_cvqr<-cv_zresidual.coxph(fit = fit_coxph,
                              data=simulated_shared_frailty_data,
                              nfolds = 10)
coxph_loocvqr<-cv_zresidual.coxph(fit = fit_coxph,
                                 data=simulated_shared_frailty_data,
                                 nfolds =nrow(simulated_shared_frailty_data))

coxph_qr_w<-zresidual.coxph (fit_coxph = fit_coxph_w, 
                             traindata = simulated_shared_frailty_data,
                             newdata = simulated_shared_frailty_data)
coxph_cvqr_w<-cv_zresidual.coxph(fit = fit_coxph_w,
                                data=simulated_shared_frailty_data,
                                nfolds = 10)
coxph_loocvqr_w<-cv_zresidual.coxph(fit = fit_coxph_w,
                                   data=simulated_shared_frailty_data,
                                   nfolds = nrow(simulated_shared_frailty_data))

a1<-which(is.infinite(coxph_cvqr))
a2<-which(is.na(coxph_cvqr))
b1<-which(is.infinite(coxph_loocvqr))
b2<-which(is.na(coxph_loocvqr))
c1<-which(is.infinite(coxph_cvqr_w))
c2<-which(is.na(coxph_cvqr_w))
d1<-which(is.infinite(coxph_loocvqr_w))
d2<-which(is.na(coxph_loocvqr_w))

coxph_cvqr = coxph_cvqr[is.finite(coxph_cvqr)]
coxph_loocvqr = coxph_loocvqr[is.finite(coxph_loocvqr)]
coxph_cvqr = coxph_cvqr[!is.na(coxph_cvqr)]
coxph_loocvqr = coxph_loocvqr[!is.na(coxph_loocvqr)]

coxph_cvqr_w = coxph_cvqr_w[is.finite(coxph_cvqr_w)]
coxph_loocvqr_w = coxph_loocvqr_w[is.finite(coxph_loocvqr_w)]
coxph_cvqr_w = coxph_cvqr_w[!is.na(coxph_cvqr_w)]
coxph_loocvqr_w = coxph_loocvqr_w[!is.na(coxph_loocvqr_w)]

sw_coxph_qr_t<-shapiro.test(coxph_qr)$p.value;sw_coxph_qr_t
sw_coxph_cvqr_t<-shapiro.test(coxph_cvqr)$p.value;sw_coxph_cvqr_t
sw_coxph_loocvqr_t<-shapiro.test(coxph_loocvqr)$p.value;sw_coxph_loocvqr_t
sw_coxph_qr_w<-shapiro.test(coxph_qr_w)$p.value;sw_coxph_qr_w
sw_coxph_cvqr_w<-shapiro.test(coxph_cvqr_w)$p.value;sw_coxph_cvqr_w
sw_coxph_loocvqr_w<-shapiro.test(coxph_loocvqr_w)$p.value;sw_coxph_loocvqr_w

#####simulation study for detecting Outliers
para.frame <- data.frame(
  n_clusters = c(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
                 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10),
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

i=5
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
fit_coxph_w1<-NA
fit_coxph_w2<-NA
fit_coxph_w3<-NA
fit_coxph_w4<-NA
while(any(is.na(fit_coxph)) || any(is.na(fit_coxph_w1)) || 
      any(is.na(fit_coxph_w2))|| any(is.na(fit_coxph_w3))|| 
      any(is.na(fit_coxph_w4)) ) {
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
  
  check_table1<-FALSE
  while(isFALSE(check_table1)){
    simulated_shared_frailty_data1<-
      simulWeib_outlier1(n_clusters=n_clusters[1,],n_individuals=n_individuals[1,],
                         lambda=lambda[1,],alpha=alpha[1,],beta1=beta1[1,],
                         beta2=beta2[1,],beta3=beta3[1,], 
                         mean.censor= mean.censor[1,],fv=fv[1,])
    check_table1<-all(apply(table(simulated_shared_frailty_data1$x3,
                                  simulated_shared_frailty_data1$d),
                            2,function(x) x>1))
  }
  
  check_table2<-FALSE
  while(isFALSE(check_table2)){
    simulated_shared_frailty_data2<-
      simulWeib_outlier2(n_clusters=n_clusters[1,],
                         n_individuals=n_individuals[1,],
                         lambda=lambda[1,],alpha=alpha[1,],beta1=beta1[1,],
                         beta2=beta2[1,],beta3=beta3[1,], 
                         mean.censor= mean.censor[1,],fv=fv[1,])
    check_table2<-all(apply(table(simulated_shared_frailty_data2$x3,
                                  simulated_shared_frailty_data2$d),
                            2,function(x) x>1))
    
  }
  check_table3<-FALSE
  while(isFALSE(check_table3)){
    simulated_shared_frailty_data3<-
      simulWeib_outlier3(n_clusters=n_clusters[1,],
                         n_individuals=n_individuals[1,],
                         lambda=lambda[1,],alpha=alpha[1,],beta1=beta1[1,],
                         beta2=beta2[1,],beta3=beta3[1,], 
                         mean.censor= mean.censor[1,],fv=fv[1,])
    check_table3<-all(apply(table(simulated_shared_frailty_data3$x3,
                                  simulated_shared_frailty_data3$d),
                            2,function(x) x>1))
  }
  check_table4<-FALSE
  while(isFALSE(check_table4)){
    simulated_shared_frailty_data4<-
      simulWeib_outlier4(n_clusters=n_clusters[1,],
                         n_individuals=n_individuals[1,],
                         lambda=lambda[1,],alpha=alpha[1,],
                         beta1=beta1[1,],
                         beta2=beta2[1,],beta3=beta3[1,], 
                         mean.censor= mean.censor[1,],fv=fv[1,])
    check_table4<-all(apply(table(simulated_shared_frailty_data4$x3,
                                  simulated_shared_frailty_data4$d),
                            2,function(x) x>1))
  }
  
  
  ####fit true model##################################################### 
  fit_coxph <- tryCatch(
    coxph(Surv(t, d) ~ x1+x2+x3 +frailty(grpid,distribution = "gamma"),
          data = simulated_shared_frailty_data),
    error = function(e) NA,
    warning = function(w) NA
  )
  
  
  ####10% outliers dataset fit model ##################################
  fit_coxph_w1 <- tryCatch(
    coxph(Surv(t, d) ~ x1+x2+x3 +frailty(grpid,distribution = "gamma"),
          data = simulated_shared_frailty_data1),
    error = function(e) NA,
    warning = function(w) NA
  )
  fit_coxph_w2 <- tryCatch(
    coxph(Surv(t, d) ~ x1+x2+x3 +frailty(grpid,distribution = "gamma"),
          data = simulated_shared_frailty_data2),
    error = function(e) NA,
    warning = function(w) NA
  ) 
  
  
  ####10 outliers dataset fit model ##################################
  fit_coxph_w3 <- tryCatch(
    coxph(Surv(t, d) ~ x1+x2+x3 +frailty(grpid,distribution = "gamma"),
          data = simulated_shared_frailty_data3),
    error = function(e) NA,
    warning = function(w) NA
  )
  fit_coxph_w4 <- tryCatch(
    coxph(Surv(t, d) ~ x1+x2+x3 +frailty(grpid,distribution = "gamma"),
          data = simulated_shared_frailty_data4),
    error = function(e) NA,
    warning = function(w) NA
  )
}
censorship1<-as.numeric(table(simulated_shared_frailty_data$d)[1]/
                          nrow(simulated_shared_frailty_data))
censorship2<-as.numeric(table(simulated_shared_frailty_data2$d)[1]/
                          nrow(simulated_shared_frailty_data2))
censorship4<-as.numeric(table(simulated_shared_frailty_data4$d)[1]/
                          nrow(simulated_shared_frailty_data4))

coxph_qr<-zresidual.coxph (fit_coxph = fit_coxph,
                           traindata = simulated_shared_frailty_data,
                           newdata = simulated_shared_frailty_data)
coxph_cvqr<-cv_zresidual.coxph(fit = fit_coxph,
                              data=simulated_shared_frailty_data,
                              nfolds = 10)
coxph_loocvqr<-cv_zresidual.coxph(fit = fit_coxph,
                                 data=simulated_shared_frailty_data,
                                 nfolds =nrow(simulated_shared_frailty_data))

coxph_qr_w1<-zresidual.coxph (fit_coxph = fit_coxph_w1, 
                              traindata = simulated_shared_frailty_data1,
                              newdata = simulated_shared_frailty_data1)
coxph_cvqr_w1<-cv_zresidual.coxph(fit = fit_coxph_w1,
                                 data=simulated_shared_frailty_data1,
                                 nfolds = 10)
coxph_loocvqr_w1<-cv_zresidual.coxph(fit = fit_coxph_w1,
                                    data=simulated_shared_frailty_data1,
                                    nfolds = nrow(simulated_shared_frailty_data1))

coxph_qr_w2<-zresidual.coxph (fit_coxph = fit_coxph_w2, 
                              traindata = simulated_shared_frailty_data2,
                              newdata = simulated_shared_frailty_data2)
coxph_cvqr_w2<-cv_zresidual.coxph(fit = fit_coxph_w2,
                                 data=simulated_shared_frailty_data2,
                                 nfolds = 10)
coxph_loocvqr_w2<-cv_zresidual.coxph(fit = fit_coxph_w2,
                                    data=simulated_shared_frailty_data2,
                                    nfolds = nrow(simulated_shared_frailty_data2))

coxph_qr_w3<-zresidual.coxph (fit_coxph = fit_coxph_w3, 
                              traindata = simulated_shared_frailty_data3,
                              newdata = simulated_shared_frailty_data3)
coxph_cvqr_w3<-cv_zresidual.coxph(fit = fit_coxph_w3,
                                 data=simulated_shared_frailty_data3,
                                 nfolds = 10)
coxph_loocvqr_w3<-cv_zresidual.coxph(fit = fit_coxph_w3,
                                    data=simulated_shared_frailty_data3,
                                    nfolds = nrow(simulated_shared_frailty_data3))

coxph_qr_w4<-zresidual.coxph (fit_coxph = fit_coxph_w4, 
                              traindata = simulated_shared_frailty_data4,
                              newdata = simulated_shared_frailty_data4)
coxph_cvqr_w4<-cv_zresidual.coxph(fit = fit_coxph_w4,
                                 data=simulated_shared_frailty_data4,
                                 nfolds = 10)
coxph_loocvqr_w4<-cv_zresidual.coxph(fit = fit_coxph_w4,
                                    data=simulated_shared_frailty_data4,
                                    nfolds = nrow(simulated_shared_frailty_data4))

a1<-which(is.infinite(coxph_cvqr_w1))
a2<-which(is.na(coxph_cvqr_w1))

b1<-which(is.infinite(coxph_loocvqr_w1))
b2<-which(is.na(coxph_loocvqr_w1))

c1<-which(is.infinite(coxph_cvqr_w2))
c2<-which(is.na(coxph_cvqr_w2))

d1<-which(is.infinite(coxph_loocvqr_w2))
d2<-which(is.na(coxph_loocvqr_w2))

e1<-which(is.infinite(coxph_cvqr_w3))
e2<-which(is.na(coxph_cvqr_w3))

f1<-which(is.infinite(coxph_loocvqr_w3))
f2<-which(is.na(coxph_loocvqr_w3))

g1<-which(is.infinite(coxph_cvqr_w4))
g2<-which(is.na(coxph_cvqr_w4))

h1<-which(is.infinite(coxph_loocvqr_w4))
h2<-which(is.na(coxph_loocvqr_w4))


coxph_cvqr = coxph_cvqr[is.finite(coxph_cvqr)]
coxph_loocvqr = coxph_loocvqr[is.finite(coxph_loocvqr)]
coxph_cvqr = coxph_cvqr[!is.na(coxph_cvqr)]
coxph_loocvqr = coxph_loocvqr[!is.na(coxph_loocvqr)]

coxph_cvqr_w1 = coxph_cvqr_w1[is.finite(coxph_cvqr_w1)]
coxph_loocvqr_w1 = coxph_loocvqr_w1[is.finite(coxph_loocvqr_w1)]
coxph_cvqr_w1 = coxph_cvqr_w1[!is.na(coxph_cvqr_w1)]
coxph_loocvqr_w1 = coxph_loocvqr_w1[!is.na(coxph_loocvqr_w1)]

coxph_cvqr_w2 = coxph_cvqr_w2[is.finite(coxph_cvqr_w2)]
coxph_loocvqr_w2 = coxph_loocvqr_w2[is.finite(coxph_loocvqr_w2)]
coxph_cvqr_w2 = coxph_cvqr_w2[!is.na(coxph_cvqr_w2)]
coxph_loocvqr_w2 = coxph_loocvqr_w2[!is.na(coxph_loocvqr_w2)]

coxph_cvqr_w3 = coxph_cvqr_w3[is.finite(coxph_cvqr_w3)]
coxph_loocvqr_w3 = coxph_loocvqr_w3[is.finite(coxph_loocvqr_w3)]
coxph_cvqr_w3 = coxph_cvqr_w3[!is.na(coxph_cvqr_w3)]
coxph_loocvqr_w3 = coxph_loocvqr_w3[!is.na(coxph_loocvqr_w3)]

coxph_cvqr_w4 = coxph_cvqr_w4[is.finite(coxph_cvqr_w4)]
coxph_loocvqr_w4 = coxph_loocvqr_w4[is.finite(coxph_loocvqr_w4)]
coxph_cvqr_w4 = coxph_cvqr_w4[!is.na(coxph_cvqr_w4)]
coxph_loocvqr_w4 = coxph_loocvqr_w4[!is.na(coxph_loocvqr_w4)]



######real data analysis#######################################
###delete outlier can be solve problem

data(kidney)
kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
kidney$id<-as.factor(kidney$id)
fit_kidney1 <- tryCatch(
  coxph(Surv(time, status) ~ age + sex + disease+
          frailty(id, distribution="gamma"), data= kidney),
  error = function(e) NA,
  warning = function(w) NA
)

fit_kidney1_qr<-zresidual.coxph(fit_coxph = fit_kidney1,
                                traindata = kidney,
                                newdata = kidney)
fit_kidney1_loocvqr<-cv_zresidual.coxph(fit = fit_kidney1,
                                       data=kidney,
                                       nfolds = nrow(kidney))

a1<-which(is.infinite(fit_kidney1_loocvqr))
a2<-which(is.na(fit_kidney1_loocvqr))

fit_kidney1_loocvqr = fit_kidney1_loocvqr[is.finite(fit_kidney1_loocvqr)]
fit_kidney1_loocvqr = fit_kidney1_loocvqr[!is.na(fit_kidney1_loocvqr)]

sw_nrsp_qr<-shapiro.test(fit_kidney1_qr)$p.value;sw_nrsp_qr
sw_nrsp_loocvqr<-shapiro.test(fit_kidney1_loocvqr)$p.value;sw_nrsp_loocvqr
