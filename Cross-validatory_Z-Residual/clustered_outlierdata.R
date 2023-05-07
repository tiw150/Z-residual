### Simulate data with 10% outliers(moderate jitters)

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

### Simulate data with 10% outliers (strong jitters)
simulWeib_outlier2 <- function(n_clusters,n_individuals, lambda, alpha,
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
  t[outlier]<- t[outlier]+ pmax(rexp(length(outlier),rate=1/2),4)
  outlier_indicator<-rep(0,n)
  outlier_indicator[outlier]<-1
  # Make data frame
  grpid<- as.factor(grpid)
  x3<-as.factor(x3)
  out<-data.frame(id,grpid,x1,x2,x3,t,d,frvec,outlier_indicator)
  return (out)
}



### Simulate data with 10 outliers(moderate jitters)
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

### Simulate data with 10 outliers(strong jitters)
simulWeib_outlier3 <- function(n_clusters,n_individuals, lambda, alpha,beta1,beta2,beta3,
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
  t[outlier]<- t[outlier]+ pmax(rexp(length(outlier),rate=1/2),4)
  outlier_indicator<-rep(0,n)
  outlier_indicator[outlier]<-1
  # Make data frame
  grpid<- as.factor(grpid)
  x3<-as.factor(x3)
  out<-data.frame(id,grpid,x1,x2,x3,t,d,frvec,outlier_indicator)
  return (out)
}

