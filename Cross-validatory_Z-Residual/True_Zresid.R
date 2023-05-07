true.zresidual <- function (data,lambda=0.007,alpha=3,beta1=1,beta2=-2,beta3=0.5)
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