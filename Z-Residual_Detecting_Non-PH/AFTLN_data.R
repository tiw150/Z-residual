library("survival")
rexp2 <- function(n, rate){ if (rate==0) rep(Inf,n) else rexp(n=n, rate = rate)}
AFTLN_data <- function(n_clusters,n_individuals, mu, sigma,beta1,
                      mean.censor,fv)
{
  # Make covariates 
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  grpid<- paste0("g",grpid)
  X = round(rnorm(n_clusters*n_individuals, 0,1),2)
  frvec<-rep(rgamma(n_clusters,shape = 1/fv,scale = fv), each = n_individuals)
  #Event times:baseline hazard is AFT lognormal
  n<-n_individuals * n_clusters
  v <- runif(n)
  surv_time <- exp(sigma* qnorm(v)+mu)*(1/exp(X*beta1+log(frvec)))
  # Censoring times
  t0<- rexp2(n, rate= 1/mean.censor)
  # follow-up times and event indicators
  t <- pmin(surv_time, t0)
  d <- as.numeric(surv_time <= t0)
  # Make data frame
  grpid<- as.factor(grpid)
#  x3<-as.factor(x3)
  out<-data.frame(id,grpid,X,t,d,frvec)
  return (out) 
}

AFTWB_data <- function(n_clusters,n_individuals, lambda, alpha,beta1,
                       mean.censor,fv)
{
  # Make covariates 
  id = 1:(n_individuals * n_clusters)
  grpid = rep(1:n_clusters, each = n_individuals)
  grpid<- paste0("g",grpid)
  X = round(rnorm(n_clusters*n_individuals, 0,1),2)
  frvec<-rep(rgamma(n_clusters,shape = 1/fv,scale = fv), each = n_individuals)
  
  n<-n_individuals * n_clusters
  v <- runif(n)
  #Event times:baseline hazard is Cox weibull
  surv_time <- (1/lambda^(1/alpha))*(-log(1-v)/exp(X*beta1+log(frvec)))^(1/alpha)
  # Censoring times
  t0<- rexp2(n, rate= 1/mean.censor)
  # follow-up times and event indicators
  t <- pmin(surv_time, t0)
  d <- as.numeric(surv_time <= t0)
  # Make data frame
  grpid<- as.factor(grpid)
#  x3<-as.factor(x3)
  out<-data.frame(id,grpid,X,t,d,frvec)
  return (out) 
}
