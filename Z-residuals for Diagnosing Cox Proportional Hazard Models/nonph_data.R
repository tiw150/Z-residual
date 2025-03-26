#set.seed(1)
rexp2 <- function(n, rate){ if (rate==0) rep(Inf,n) else rexp(n = n, rate = rate)}
########beta2 = beta_a if  0<surv_time< t0 , beta2 = beta_b if  t0 <surv_time 
simul_nonph <- function(n_clusters,n_individuals,lambda, alpha, 
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
  x1 = pmax(abs(round(covar[,1],1)),0.01)
  x2 = round(covar[,2],1)
  
  #Event times:baseline hazard is Weibull，lambda=scale parameter, alpha=shape parameter
  v <- runif(n)
  a<- -log(v)
  b<- lambda*exp((x1)*beta1+beta_a*x2) * (t0^alpha)
  surv_time<-rep(0,n)
  group1<-which(a<b)
  group2<-which(a>=b)
  surv_time[group1]<- (-log(v)[group1]/
                         (lambda*exp((x1)[group1]*beta1+x2[group1]*beta_a)))^(1 /alpha)
  
  surv_time[group2]<- ((-log(v)[group2]-
                          ((t0^alpha)*lambda*
                             exp((x1)[group2]*beta1+x2[group2]*beta_a))+
                          ((t0^alpha)*lambda*
                             exp((x1)[group2]*beta1+x2[group2]*beta_b)))
                       /(lambda*exp((x1)[group2]*beta1)
                         *exp(x2[group2]*beta_b)))^(1 /alpha)
  # Censoring times
  c0<- rexp2(n, rate= 1/mean.censor)
  # follow-up times and event indicators
  t <- pmin(surv_time, c0)
  d <- as.numeric(surv_time <= c0)
  # Make data frame
  grpid<- as.factor(grpid)
  # x3<-as.factor(x3)
  out<-data.frame(id,grpid,x1,x2,t,d)
  return (out) 
}