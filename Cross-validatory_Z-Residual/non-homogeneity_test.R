test.nl.aov <- function(zresidual, fitted.values, k.anova=10)
{
  lpred.bin <- cut(fitted.values, k.anova)
  less2_factor<-which(tapply(lpred.bin,lpred.bin,length)<= 2)
  if(rlang::is_empty(names(less2_factor))){
    anova(lm(zresidual ~ lpred.bin))$`Pr(>F)`[1]
  }else{
    vector_less2_factor<-rep(length(less2_factor))
    for(j in 1:length(less2_factor)){
      vector_less2_factor[j]<-which(lpred.bin==names(less2_factor[j]))
    }
    new.lpred.bin<- lpred.bin[-vector_less2_factor]
    new.zresidual<-zresidual[-vector_less2_factor]
    anova(lm(zresidual ~ new.lpred.bin))$`Pr(>F)`[1]
  }
}


test.var.bartl <- function(zresidual, fitted.values, k=10)
{
  lpred.bin <- cut(fitted.values, k)
  Z_group<- split(zresidual, lpred.bin)
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
