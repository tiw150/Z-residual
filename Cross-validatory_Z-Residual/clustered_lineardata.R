### Simulate data
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
  # Make data frame
  grpid<- as.factor(grpid)
  x3<-as.factor(x3)
  out<-data.frame(id,grpid,x1,x2,x3,t,d,frvec)
  return (out) 
}


true.qresidual <- function (data,lambda=0.007,alpha=3,beta1=1,beta2=-2,beta3=0.5)
{
  explp<-exp(data$x1*beta1+log(data$x2)*beta2+as.numeric(as.character(data$x3))*beta3)
  z<- data$frvec
  H0<-lambda*(data$t)^alpha
  #Survival Function
  SP<- exp(-z*as.vector(explp)*H0)
  censored <- which(data$d==0)
  n.censored <- length(censored)
  #NRSP residual
  RSP <- SP
  RSP[censored] <- RSP[censored]*runif(n.censored)
  nrsp <- -qnorm(RSP)
  attr(nrsp, "SP") <- SP
  attr(nrsp, "z") <- z
  return(nrsp)
}