setwd("~/Cross-validatory_Z-Residual")
source("~/clustered_nonlineardata.R")
source("~/split_kfold.R")
source("~/Zresid_coxph.R")
source("~/cv_Zresid_coxph.R")
source("~/non-homogeneity_test.R")
source("~/pmin_value")

library("survival")
library("EnvStats")
library("nortest")
library("backports")
library("stringr")
library("pROC")

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
i=2
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

true_qr<- true.zresidual(data=simulated_shared_frailty_data,
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

compare_sp<-data.frame("ID"=simulated_shared_frailty_data$id,
                       "GroupID"=simulated_shared_frailty_data$grpid,
                       "t"=simulated_shared_frailty_data$t,
                       "d"=simulated_shared_frailty_data$d,
                       "true_qr"=true_qr, "true_sp"=attr(true_qr,"SP"),
                       "true_z"=attr(true_qr,"z"),
                       "nocv_qr"=coxph_qr,"nocv_sp"=attr(coxph_qr,"SP"),
                       "nocv_z"=attr(coxph_qr,"z_hat_new"),
                       "tenfcv_qr"=coxph_cvqr,
                       "tenfcv_sp"=attr(coxph_cvqr,"SP"),
                       "tenfcv_z"=attr(coxph_cvqr,"z_hat_new"),
                       "loocv_qr"=coxph_loocvqr,
                       "loocv_sp"=attr(coxph_loocvqr,"SP"),
                       "loocv_z"=attr(coxph_loocvqr,"z_hat_new"),
                       "nocv_qr_w"=coxph_qr_w,"nocv_sp_w"=attr(coxph_qr_w,"SP"),
                       "nocv_z_w"=attr(coxph_qr_w,"z_hat_new"),
                       "tenfcv_qr_w"=coxph_cvqr_w,
                       "tenfcv_sp_w"=attr(coxph_cvqr_w,"SP"),
                       "tenfcv_z_w"=attr(coxph_cvqr_w,"z_hat_new"),
                       "loocv_qr_w"=coxph_loocvqr_w,
                       "loocv_sp_w"=attr(coxph_loocvqr_w,"SP"),
                       "loocv_z_w"=attr(coxph_loocvqr_w,"z_hat_new"))

rsq_coxph_qr_t<-(cor(compare_sp$true_sp,compare_sp$nocv_sp))^2
rsq_coxph_cvqr_t<-(cor(compare_sp$true_sp,compare_sp$tenfcv_sp))^2
rsq_coxph_loocvqr_t<-(cor(compare_sp$true_sp,compare_sp$loocv_sp))^2
rsq_coxph_qr_w<-(cor(compare_sp$true_sp,compare_sp$nocv_sp_w))^2
rsq_coxph_cvqr_w<-(cor(compare_sp$true_sp,compare_sp$tenfcv_sp_w))^2
rsq_coxph_loocvqr_w<-(cor(compare_sp$true_sp,compare_sp$loocv_sp_w))^2


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

# save(coxph_qr,coxph_cvqr,coxph_loocvqr,
#      coxph_qr_w,coxph_cvqr_w,coxph_loocvqr_w,
#      a1,a2,b1,b2,c1,c2,d1,d2,simulated_shared_frailty_data,
#      sw_coxph_qr_t,sw_coxph_cvqr_t,sw_coxph_loocvqr_t,
#      sw_coxph_qr_w,sw_coxph_cvqr_w,sw_coxph_loocvqr_w,
#      fix_var,risk_score_t,risk_score_w,file="sigledata_n500c50.RData")


sw_coxph_qr_t<-shapiro.test(coxph_qr)$p.value;sw_coxph_qr_t
sw_coxph_cvqr_t<-shapiro.test(coxph_cvqr)$p.value;sw_coxph_cvqr_t
sw_coxph_loocvqr_t<-shapiro.test(coxph_loocvqr)$p.value;sw_coxph_loocvqr_t
sw_coxph_qr_w<-shapiro.test(coxph_qr_w)$p.value;sw_coxph_qr_w
sw_coxph_cvqr_w<-shapiro.test(coxph_cvqr_w)$p.value;sw_coxph_cvqr_w
sw_coxph_loocvqr_w<-shapiro.test(coxph_loocvqr_w)$p.value;sw_coxph_loocvqr_w

pdf(sprintf ("nonlinear_residualplot2.pdf"),width=12,height=8)
par(mfrow = c(2,3))
#par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0))
ymax=max(range(coxph_qr,coxph_cvqr,coxph_loocvqr,
               coxph_qr_w,coxph_cvqr_w,coxph_loocvqr_w))
ymin=min(range(coxph_qr,coxph_cvqr,coxph_loocvqr,
               coxph_qr_w,coxph_cvqr_w,coxph_loocvqr_w))
plot(simulated_shared_frailty_data$x2,coxph_qr,
     xlab=expression(paste(x[2])),ylab="Z Residual",
     col=c("blue","darkolivegreen4")[simulated_shared_frailty_data$d+1],
     pch=c(3,2)[simulated_shared_frailty_data$d+1],
     main="True Model, Z Residual plot, No CV",
     ylim=c(ymin,ymax))
abline(h=c(3,-3),col="grey")
legend(x = "bottom", inset = c(0, -0.17), 
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")
plot((simulated_shared_frailty_data$x2)[-c(a1,a2)],coxph_cvqr,
     xlab=expression(paste(x[2])),ylab="Z Residual",
     col=c("blue","darkolivegreen4")[simulated_shared_frailty_data$d+1],
     pch=c(3,2)[simulated_shared_frailty_data$d+1],
     main="True Model, Z Residual plot, 10-fold",
     ylim=c(ymin,ymax))
abline(h=c(3,-3),col="grey")
legend(x = "bottom", inset = c(0, -0.17), 
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")
plot((simulated_shared_frailty_data$x2)[-c(b1,b2)],coxph_loocvqr,
     xlab=expression(paste(x[2])),ylab="Z Residual",
     col=c("blue","darkolivegreen4")[simulated_shared_frailty_data$d+1],
     pch=c(3,2)[simulated_shared_frailty_data$d+1],
     main="True Model, Z Residual plot, LOOCV",
     ylim=c(ymin,ymax))
abline(h=c(3,-3),col="grey")
legend(x = "bottom", inset = c(0, -0.17), 
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")

plot(simulated_shared_frailty_data$x2,coxph_qr_w,
     xlab=expression(paste(x[2])),ylab="Z Residual",
     col=c("blue","darkolivegreen4")[simulated_shared_frailty_data$d+1],
     pch=c(3,2)[simulated_shared_frailty_data$d+1],
     main="Wrong Model, Z Residual plot, No CV",
     ylim=c(ymin,ymax))
abline(h=c(3,-3),col="grey")
legend(x = "bottom", inset = c(0, -0.17), 
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")
plot((simulated_shared_frailty_data$x2)[-c(c1,c2)],coxph_cvqr_w,
     xlab=expression(paste(x[2])),ylab="Z Residual",
     col=c("blue","darkolivegreen4")[simulated_shared_frailty_data$d+1],
     pch=c(3,2)[simulated_shared_frailty_data$d+1],
     main="Wrong Model, Z Residual plot, 10-fold",
     ylim=c(ymin,ymax))
abline(h=c(3,-3),col="grey")
legend(x = "bottom", inset = c(0, -0.17), 
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")
plot((simulated_shared_frailty_data$x2)[-c(d1,d2)],coxph_loocvqr_w,
     xlab=expression(paste(x[2])),ylab="Z Residual",
     col=c("blue","darkolivegreen4")[simulated_shared_frailty_data$d+1],
     pch=c(3,2)[simulated_shared_frailty_data$d+1],
     main="Wrong Model, Z Residual plot, LOOCV",
     ylim=c(ymin,ymax))
abline(h=c(3,-3),col="grey")
legend(x = "bottom", inset = c(0, -0.17), 
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")
dev.off()

pdf(sprintf ("nonlinear_qqplot2.pdf"),width=12,height=8)
par(mfrow = c(2,3))
#par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0))
ymax=max(range(coxph_qr,coxph_cvqr,coxph_loocvqr,
               coxph_qr_w,coxph_cvqr_w,coxph_loocvqr_w))
ymin=min(range(coxph_qr,coxph_cvqr,coxph_loocvqr,
               coxph_qr_w,coxph_cvqr_w,coxph_loocvqr_w))
qqnorm(coxph_qr,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("True Model, QQ plot, No CV, SW p-value = ", 
                   sprintf("%3.2f",sw_coxph_qr_t)),
       ylim=c(ymin,ymax))
qqline(coxph_qr)
qqnorm(coxph_cvqr,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("True Model, QQ plot, 10-fold, SW p-value = ", 
                   sprintf("%3.2f",sw_coxph_cvqr_t)),
       ylim=c(ymin,ymax))
qqline(coxph_cvqr)
qqnorm(coxph_loocvqr,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("True Model, QQ plot, LOOCV, SW p-value = ", 
                   sprintf("%3.2f",sw_coxph_loocvqr_t)),
       ylim=c(ymin,ymax))
qqline(coxph_loocvqr)

qqnorm(coxph_qr_w,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("Wrong Model, QQ plot, No CV, SW p-value = ", 
                   sprintf("%3.2f",sw_coxph_qr_w)),
       ylim=c(ymin,ymax))
qqline(coxph_qr_w)
qqnorm(coxph_cvqr_w,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("Wrong Model, QQ plot, 10-fold, SW p-value = ", 
                   sprintf("%3.2f",sw_coxph_cvqr_w)),
       ylim=c(ymin,ymax))
qqline(coxph_cvqr_w)
qqnorm(coxph_loocvqr_w,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("Wrong Model, QQ plot, LOOCV, SW p-value = ", 
                   sprintf("%3.2f",sw_coxph_loocvqr_w)),
       ylim=c(ymin,ymax))
qqline(coxph_loocvqr_w)
dev.off()


sw_trueqr<-shapiro.test(true_qr)$p.value
sw_coxph_qr_t<-shapiro.test(coxph_qr)$p.value
sw_coxph_cvqr_t<-shapiro.test(coxph_cvqr)$p.value
sw_coxph_loocvqr_t<-shapiro.test(coxph_loocvqr)$p.value

sf_trueqr<-sf.test(true_qr)$p.value
sf_coxph_qr_t<-sf.test(coxph_qr)$p.value
sf_coxph_cvqr_t<-sf.test(coxph_cvqr)$p.value
sf_coxph_loocvqr_t<-sf.test(coxph_loocvqr)$p.value

tprob_trueqr<-mean(abs(true_qr)>1.96)
tprob_coxph_qr_t<-mean(abs(coxph_qr)>1.96)
tprob_coxph_cvqr_t<-mean(abs(coxph_cvqr)>1.96)
tprob_coxph_loocvqr_t<-mean(abs(coxph_loocvqr)>1.96)

sd_trueqr<-sd(true_qr)
sd_coxph_qr_t<-sd(coxph_qr)
sd_coxph_cvqr_t<-sd(coxph_cvqr)
sd_coxph_loocvqr_t<-sd(coxph_loocvqr)

sw_coxph_qr_w<-shapiro.test(coxph_qr_w)$p.value
sw_coxph_cvqr_w<-shapiro.test(coxph_cvqr_w)$p.value
sw_coxph_loocvqr_w<-shapiro.test(coxph_loocvqr_w)$p.value

sf_coxph_qr_w<-sf.test(coxph_qr_w)$p.value
sf_coxph_cvqr_w<-sf.test(coxph_cvqr_w)$p.value
sf_coxph_loocvqr_w<-sf.test(coxph_loocvqr_w)$p.value

tprob_coxph_qr_w<-mean(abs(coxph_qr_w)>1.96)
tprob_coxph_cvqr_w<-mean(abs(coxph_cvqr_w)>1.96)
tprob_coxph_loocvqr_w<-mean(abs(coxph_loocvqr_w)>1.96)

sd_coxph_qr_w<-sd(coxph_qr_w)
sd_coxph_cvqr_w<-sd(coxph_cvqr_w)
sd_coxph_loocvqr_w<-sd(coxph_loocvqr_w)


