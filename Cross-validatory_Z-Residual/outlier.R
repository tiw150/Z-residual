setwd("~/Cross-validatory_Z-Residual")
source("~/clustered_lineardata.R")
source("~/split_kfold.R")
source("~/Zresid_coxph.R")
source("~/cv_Zresid_coxph.R")
source("~/non-homogeneity_test.R")
source("~/pmin_value")
source("~/clustered_outlierdata.R")

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



pdf( "outlier_residualplot.pdf",width=9,height=6)
par(mfrow=c(2,3))
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0))
ymax=max(max(range(coxph_qr,coxph_cvqr,coxph_loocvqr,
               coxph_qr_w2,coxph_cvqr_w2,coxph_loocvqr_w2)),3)
ymin=min(min(range(coxph_qr,coxph_cvqr,coxph_loocvqr,
               coxph_qr_w2,coxph_cvqr_w2,coxph_loocvqr_w2)),-3)
###first row
plot(coxph_qr, 
     ylim=c(ymin,ymax),main="Clean Data, No CV",
     ylab="Z Residual")
abline(h=c(3,-3),col="grey")
plot(coxph_cvqr, 
     ylim=c(ymin,ymax),main="Clean Data, 10-fold",
     ylab="Z Residual")
abline(h=c(3,-3),col="grey")
plot(coxph_loocvqr, 
     ylim=c(ymin,ymax),main="Clean Data, LOOCV",
     ylab="Z Residual")
abline(h=c(3,-3),col="grey")
#text(coxph_cvqr_w1)

##second row
plot(coxph_qr_w2, 
     col = c("black","red")
     [simulated_shared_frailty_data3$outlier_indicator+1],
     pch=c(1,16)[simulated_shared_frailty_data3$outlier_indicator+1],
     ylim=c(ymin,ymax),main="10 Outliers Data, No CV",
     ylab="Z Residual")
abline(h=c(3,-3),col="grey")
legend(x = "bottom", inset = c(0, -0.18), 
       legend="Outliers", col="red",
       pch=16,horiz=TRUE,cex=1,xpd = TRUE,bty="n")
plot(coxph_cvqr_w2, 
     col = c("black","red")
     [simulated_shared_frailty_data3$outlier_indicator+1],
     pch=c(1,16)[simulated_shared_frailty_data3$outlier_indicator+1],
     ylim=c(ymin,ymax),main="10 Outliers Data, 10-fold",
     ylab="Z Residual")
abline(h=c(3,-3),col="grey")
legend(x = "bottom", inset = c(0, -0.18), 
       legend="Outliers", col="red",
       pch=16,horiz=TRUE,cex=1,xpd = TRUE,bty="n")
plot(coxph_loocvqr_w2, 
     col = c("black","red")
     [simulated_shared_frailty_data3$outlier_indicator+1],
     pch=c(1,16)[simulated_shared_frailty_data3$outlier_indicator+1],
     ylim=c(ymin,ymax),main="10 Outliers Data, LOOCV",
     ylab="Z Residual")
abline(h=c(3,-3),col="grey")
legend(x = "bottom", inset = c(0, -0.18), 
       legend="Outliers", col="red",
       pch=16,horiz=TRUE,cex=1,xpd = TRUE,bty="n")
dev.off()

save(simulated_shared_frailty_data,simulated_shared_frailty_data3,
     coxph_qr,coxph_qr_w2,coxph_cvqr,coxph_cvqr_w2,coxph_loocvqr,
     coxph_loocvqr_w2, file="sigledata1.RData")

sw_coxph_qr_t<-shapiro.test(coxph_qr)$p.value
sw_coxph_cvqr_t<-shapiro.test(coxph_cvqr)$p.value
sw_coxph_loocvqr_t<-shapiro.test(coxph_loocvqr)$p.value

sf_coxph_qr_t<-sf.test(coxph_qr)$p.value
sf_coxph_cvqr_t<-sf.test(coxph_cvqr)$p.value
sf_coxph_loocvqr_t<-sf.test(coxph_loocvqr)$p.value

tprob_coxph_qr_t<-mean(abs(coxph_qr)>1.96)
tprob_coxph_cvqr_t<-mean(abs(coxph_cvqr)>1.96)
tprob_coxph_loocvqr_t<-mean(abs(coxph_loocvqr)>1.96)

sw_coxph_qr_w1<-shapiro.test(coxph_qr_w1)$p.value
sw_coxph_cvqr_w1<-shapiro.test(coxph_cvqr_w1)$p.value
sw_coxph_loocvqr_w1<-shapiro.test(coxph_loocvqr_w1)$p.value

sf_coxph_qr_w1<-sf.test(coxph_qr_w1)$p.value
sf_coxph_cvqr_w1<-sf.test(coxph_cvqr_w1)$p.value
sf_coxph_loocvqr_w1<-sf.test(coxph_loocvqr_w1)$p.value

tprob_coxph_qr_w1<-mean(abs(coxph_qr_w1)>1.96)
tprob_coxph_cvqr_w1<-mean(abs(coxph_cvqr_w1)>1.96)
tprob_coxph_loocvqr_w1<-mean(abs(coxph_loocvqr_w1)>1.96)

sw_coxph_qr_w2<-shapiro.test(coxph_qr_w2)$p.value
sw_coxph_cvqr_w2<-shapiro.test(coxph_cvqr_w2)$p.value
sw_coxph_loocvqr_w2<-shapiro.test(coxph_loocvqr_w2)$p.value

sf_coxph_qr_w2<-sf.test(coxph_qr_w2)$p.value
sf_coxph_cvqr_w2<-sf.test(coxph_cvqr_w2)$p.value
sf_coxph_loocvqr_w2<-sf.test(coxph_loocvqr_w2)$p.value

tprob_coxph_qr_w2<-mean(abs(coxph_qr_w2)>1.96)
tprob_coxph_cvqr_w2<-mean(abs(coxph_cvqr_w2)>1.96)
tprob_coxph_loocvqr_w2<-mean(abs(coxph_loocvqr_w2)>1.96)

model_qr_outlier1<-simulated_shared_frailty_data2$outlier_indicator
resid_qr_outlier1<-abs(coxph_qr_w1)
rcurve_qr_outlier1<-roc(model_qr_outlier1,resid_qr_outlier1)
sensitiv_qr_outlier1<-as.numeric(coords(rcurve_qr_outlier1, 3,ret=c("sensitivity",
                                                "1-specificity"),
                               transpose = FALSE)[1])
fpr_qr_outlier1<-as.numeric(coords(rcurve_qr_outlier1, 3, ret=c("sensitivity",
                                            "1-specificity"),
                          transpose = FALSE)[2])
auc_qr_outlier1<-as.numeric(auc(rcurve_qr_outlier1))

if(length(a1)==0 && length(a2)==0){
  model_cvqr_outlier1<-simulated_shared_frailty_data2$outlier_indicator
}else{
  model_cvqr_outlier1<-simulated_shared_frailty_data2$outlier_indicator[-c(a1,a2)]
}
resid_cvqr_outlier1<-abs(coxph_cvqr_w1)
rcurve_cvqr_outlier1<-roc(model_cvqr_outlier1,resid_cvqr_outlier1)
sensitiv_cvqr_outlier1<-as.numeric(coords(rcurve_cvqr_outlier1,3, ret=c("sensitivity",
                                                  "1-specificity"),
                                 transpose = FALSE)[1])
fpr_cvqr_outlier1<-as.numeric(coords(rcurve_cvqr_outlier1, 3, ret=c("sensitivity",
                                              "1-specificity"),
                            transpose = FALSE)[2])
auc_cvqr_outlier1<-as.numeric(auc(rcurve_cvqr_outlier1))

if(length(b1)==0 && length(b2)==0){
  model_loocvqr_outlier1<-simulated_shared_frailty_data2$outlier_indicator
}else{
  model_loocvqr_outlier1<-simulated_shared_frailty_data2$outlier_indicator[-c(b1,b2)]
}
resid_loocvqr_outlier1<-abs(coxph_loocvqr_w1)
rcurve_loocvqr_outlier1<-roc(model_loocvqr_outlier1,resid_loocvqr_outlier1)
sensitiv_loocvqr_outlier1<-as.numeric(coords(rcurve_loocvqr_outlier1,3,
                                             ret=c("sensitivity","1-specificity"),
                                    transpose = FALSE)[1])
fpr_loocvqr_outlier1<-as.numeric(coords(rcurve_loocvqr_outlier1, 3, 
                                        ret=c("sensitivity","1-specificity"),
                               transpose = FALSE)[2])
auc_loocvqr_outlier1<-as.numeric(auc(rcurve_loocvqr_outlier1))

model_qr_outlier2<-simulated_shared_frailty_data3$outlier_indicator
resid_qr_outlier2<-abs(coxph_qr_w2)
rcurve_qr_outlier2<-roc(model_qr_outlier2,resid_qr_outlier2)
sensitiv_qr_outlier2<-as.numeric(coords(rcurve_qr_outlier2, 3,
                                        ret=c("sensitivity","1-specificity"),
                                        transpose = FALSE)[1])
fpr_qr_outlier2<-as.numeric(coords(rcurve_qr_outlier2, 3,
                                   ret=c("sensitivity","1-specificity"),
                                   transpose = FALSE)[2])
auc_qr_outlier2<-as.numeric(auc(rcurve_qr_outlier2))

if(length(c1)==0 && length(c2)==0){
  model_cvqr_outlier2<-simulated_shared_frailty_data3$outlier_indicator
}else{
  model_cvqr_outlier2<-simulated_shared_frailty_data3$outlier_indicator[-c(c1,c2)]
}
resid_cvqr_outlier2<-abs(coxph_cvqr_w2)
rcurve_cvqr_outlier2<-roc(model_cvqr_outlier2,resid_cvqr_outlier2)
sensitiv_cvqr_outlier2<-as.numeric(coords(rcurve_cvqr_outlier2,3, 
                                          ret=c("sensitivity","1-specificity"),
                                          transpose = FALSE)[1])
fpr_cvqr_outlier2<-as.numeric(coords(rcurve_cvqr_outlier2, 3, 
                                     ret=c("sensitivity","1-specificity"),
                                     transpose = FALSE)[2])
auc_cvqr_outlier2<-as.numeric(auc(rcurve_cvqr_outlier2))

if(length(d1)==0 && length(d2)==0){
  model_loocvqr_outlier2<-simulated_shared_frailty_data3$outlier_indicator
}else{
  model_loocvqr_outlier2<-simulated_shared_frailty_data3$outlier_indicator[-c(d1,d2)]
}
resid_loocvqr_outlier2<-abs(coxph_loocvqr_w2)
rcurve_loocvqr_outlier2<-roc(model_loocvqr_outlier2,resid_loocvqr_outlier2)
sensitiv_loocvqr_outlier2<-as.numeric(coords(rcurve_loocvqr_outlier2,3, 
                                             ret=c("sensitivity","1-specificity"),
                                    transpose = FALSE)[1])
fpr_loocvqr_outlier2<-as.numeric(coords(rcurve_loocvqr_outlier2, 3,
                                        ret=c("sensitivity","1-specificity"),
                                        transpose = FALSE)[2])
auc_loocvqr_outlier2<-as.numeric(auc(rcurve_loocvqr_outlier2))







