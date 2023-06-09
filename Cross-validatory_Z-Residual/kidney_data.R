setwd("~/Cross-validatory_Z-Residual")
source("~/split_kfold.R")
source("~/Zresid_coxph.R")
source("~/cv_Zresid_coxph.R")
source("~/non-homogeneity_test.R")
source("~/pmin_value")
source("~/csresidual.R")
source("~/cvcsresidual.R")

library("survival")
library("EnvStats")
library("nortest")
library("dplyr")
library("frailtyHL")


#Data Description
#Data on the recurrence times to infection, at the point of insertion of the catheter, for kidney patients using portable dialysis equipment. Catheters may be removed for reasons other than infection, in which case the observation is censored. Each patient has exactly 2 observations. This data has often been used to illustrate the use of random effects (frailty) in a survival model. However, one of the males (id 21) is a large outlier, with much longer survival than his peers. If this observation is removed no evidence remains for a random subject effect.

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
AIC(fit_kidney1)

fix_var<- model.matrix(~.,data = kidney[,c(4,5,6),drop=FALSE])[,-1,drop=FALSE]
risk_score<-fix_var %*% fit_kidney1$coefficients+fit_kidney1$frail[kidney$id]


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

fit_kidney1_csr<-csresidual.coxph(fit = fit_kidney1,
                                  traindata = kidney,
                                  newdata = kidney)
fit_kidney1_loocvcsr<-cvcsresidual.coxph(fit = fit_kidney1,
                                      data=kidney,
                                      nfolds = nrow(kidney))

km.ucs.t.loocv <- survfit(Surv(fit_kidney1_loocvcsr,kidney$status)~1,
                          type='fleming')
id.ucs.t.loocv<-order(fit_kidney1_loocvcsr)

km.ucs.t <- survfit(Surv(fit_kidney1_csr,kidney$status)~1,type='fleming')
id.ucs.t<-order(fit_kidney1_csr)

#mean(sw.fit_kidney1.qr < 0.05) is 0.002
#mean(sw.fit_kidney1.cvqr<0.05) is 0.993

pdf("kidneyplot.pdf",width=12,height=6)
par(mfrow = c(2,4))
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0))
ymax1=max(max(range(fit_kidney1_qr,fit_kidney1_loocvqr)),3)
ymin1=min(min(range(fit_kidney1_qr,fit_kidney1_loocvqr)),-3)
ymax2=max(range(km.ucs.t$cumhaz,km.ucs.t.loocv$cumhaz))
ymin2=min(range(km.ucs.t$cumhaz,km.ucs.t.loocv$cumhaz))
xmax1=max(range(sw_kidney_qr,sw_kidney_loocvqr,
                sw_qr_no_outlier,sw_loocvqr_no_outlier))
xmin1=min(range(sw_kidney_qr,sw_kidney_loocvqr,
                sw_qr_no_outlier,sw_loocvqr_no_outlier))
xmax2=max(range(km.ucs.t$time,km.ucs.t.loocv$time))
xmin2=min(range(km.ucs.t$time,km.ucs.t.loocv$time))
####the first row: No CV method

plot(fit_kidney1_qr,ylab="Z Residual plot",
     col=c("blue","darkolivegreen4")[kidney$status+1],
     pch=c(3,2)[kidney$status+1],
     main="Z Residual plot, No CV",
     ylim=c(ymin1,ymax1))
abline(h=c(3,-3),col="grey")
symbols(which(abs(fit_kidney1_loocvqr)>3),
        (fit_kidney1_qr)[which(abs(fit_kidney1_loocvqr)>3)],
        circles=c(2, 2),fg=c('red', 'red'),
        add=T, inches=F)
text(which(abs(fit_kidney1_loocvqr)>3),
     (fit_kidney1_qr)[which(abs(fit_kidney1_loocvqr)>3)],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 1,col="red")
legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)
qqnorm(fit_kidney1_qr,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("QQ plot, No CV, SW p-value = ", 
                   sprintf("%3.2f",sw_nrsp_qr)),
       ylim=c(ymin1,ymax1))
qqline(fit_kidney1_qr)
symbols((qqnorm(fit_kidney1_qr,plot.it=FALSE)
         $x)[which(abs(fit_kidney1_loocvqr)>3)],
        (qqnorm(fit_kidney1_qr,plot.it=FALSE)
         $y)[which(abs(fit_kidney1_loocvqr)>3)],
        circles=c(0.1, 0.1),fg=c('red', 'red'),
        add=T, inches=F)
text((qqnorm(fit_kidney1_qr,plot.it=FALSE)
      $x)[which(abs(fit_kidney1_loocvqr)>3)],
     (qqnorm(fit_kidney1_qr,plot.it=FALSE)
      $y)[which(abs(fit_kidney1_loocvqr)>3)],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 0.8,col="red")
hist(sw_kidney_qr,main="Replicated SW p-values, No CV",
     xlab="SW p-values",xlim=c(xmin1,xmax1))
plot(km.ucs.t, fun="cumhaz", xlab="Unmodified Cox-Snell Residuals",
     ylab="Cumulative Hazard Function",
     ylim=c(ymin2,ymax2),xlim=c(xmin2,xmax2),
     main="Cum. Hazard of CS Residuals, No CV")
abline(0, 1, col="chocolate4", lty=2)
points(km.ucs.t$time, -log(km.ucs.t$surv), 
       col=c("blue","darkolivegreen4")[kidney$status[id.ucs.t]+1],
       pch=c(3,2)[kidney$status[id.ucs.t]+1])
symbols(km.ucs.t$time[c(which(id.ucs.t==20),
                        which(id.ucs.t==42))],
        (-log(km.ucs.t$surv))[c(which(id.ucs.t==20),
                                which(id.ucs.t==42))],
        circles=c(0.3,0.3),fg=c('red', 'red'),
        add=T, inches=F)
text(km.ucs.t$time[c(which(id.ucs.t==20),
                     which(id.ucs.t==42))],
     (-log(km.ucs.t$surv))[c(which(id.ucs.t==20),
                             which(id.ucs.t==42))],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 1,col="red")
legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)

####the second row: LOOCV method########
plot(fit_kidney1_loocvqr,ylab="Z Residual plot",
     col=c("blue","darkolivegreen4")[kidney$status+1],
     pch=c(3,2)[kidney$status+1],
     main="Z Residual plot, LOOCV",
     ylim=c(ymin1,ymax1))
abline(h=c(3,-3),col="grey")
symbols(which(abs(fit_kidney1_loocvqr)>3),
        fit_kidney1_loocvqr[which(abs(fit_kidney1_loocvqr)>3)],
        circles=c(2, 2),fg=c('red', 'red'),
        add=T, inches=F)
text(which(abs(fit_kidney1_loocvqr)>3),
     fit_kidney1_loocvqr[which(abs(fit_kidney1_loocvqr)>3)],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 1,col="red")
legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)
qqnorm(fit_kidney1_loocvqr,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("QQ plot, LOOCV, SW p-value = ", 
                   sprintf("%3.2f",sw_nrsp_loocvqr)),
       ylim=c(ymin1,ymax1))
qqline(fit_kidney1_loocvqr)
symbols((qqnorm(fit_kidney1_loocvqr,plot.it=FALSE)
         $x)[which(abs(fit_kidney1_loocvqr)>3)],
        (qqnorm(fit_kidney1_loocvqr,plot.it=FALSE)
         $y)[which(abs(fit_kidney1_loocvqr)>3)],
        circles=c(0.1, 0.1),fg=c('red', 'red'),
        add=T, inches=F)
text((qqnorm(fit_kidney1_loocvqr,plot.it=FALSE)
      $x)[which(abs(fit_kidney1_loocvqr)>3)],
     (qqnorm(fit_kidney1_loocvqr,plot.it=FALSE)
      $y)[which(abs(fit_kidney1_loocvqr)>3)],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 0.8,col="red")
hist(sw_kidney_loocvqr,main="Replicated SW p-values, LOOCV",
     xlab="SW p-values",xlim=c(xmin1,xmax1))
plot(km.ucs.t.loocv, fun="cumhaz", xlab="Unmodified Cox-Snell Residuals",
     ylab="Cumulative Hazard Function",
     ylim=c(ymin2,ymax2),xlim=c(xmin2,xmax2),
     main="Cum. Hazard of CS Residuals, LOOCV")
abline(0, 1, col="chocolate4", lty=2)
points(km.ucs.t.loocv$time, -log(km.ucs.t.loocv$surv), 
       col=c("blue","darkolivegreen4")[kidney$status[id.ucs.t.loocv]+1],
       pch=c(3,2)[kidney$status[id.ucs.t.loocv]+1])
symbols(km.ucs.t.loocv$time[c(which(id.ucs.t.loocv==20),
                        which(id.ucs.t.loocv==42))],
        (-log(km.ucs.t.loocv$surv))[c(which(id.ucs.t.loocv==20),
                                which(id.ucs.t.loocv==42))],
        circles=c(0.3,0.3),fg=c('red', 'red'),
        add=T, inches=F)
text(km.ucs.t.loocv$time[c(which(id.ucs.t.loocv==20),
                     which(id.ucs.t.loocv==42))],
     (-log(km.ucs.t.loocv$surv))[c(which(id.ucs.t.loocv==20),
                             which(id.ucs.t.loocv==42))],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 1,col="red")
legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)

dev.off()

data(kidney)
kidney<-kidney[-c(20,42),]
kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
kidney$id<-as.factor(kidney$id)
fit_kidney1 <- tryCatch(
  coxph(Surv(time, status) ~ age + sex + disease+
          frailty(id, distribution="gamma"), data= kidney),
  error = function(e) NA,
  warning = function(w) NA
)

fix_var<- model.matrix(~.,data = kidney[,c(4,5,6),drop=FALSE])[,-1,drop=FALSE]
risk_score<-fix_var %*% fit_kidney1$coefficients+fit_kidney1$frail[kidney$id]

fit_kidney1_qr<-qresidual.coxph(fit_coxph = fit_kidney1,
                                traindata = kidney,
                                newdata = kidney)
fit_kidney1_loocvqr<-cvqresidual.coxph(fit = fit_kidney1,
                                       data=kidney,
                                       nfolds = nrow(kidney))

a1<-which(is.infinite(fit_kidney1_loocvqr))
a2<-which(is.na(fit_kidney1_loocvqr))

fit_kidney1_loocvqr = fit_kidney1_loocvqr[is.finite(fit_kidney1_loocvqr)]
fit_kidney1_loocvqr = fit_kidney1_loocvqr[!is.na(fit_kidney1_loocvqr)]

sw_nrsp_qr<-shapiro.test(fit_kidney1_qr)$p.value;sw_nrsp_qr
sw_nrsp_loocvqr<-shapiro.test(fit_kidney1_loocvqr)$p.value;sw_nrsp_loocvqr

fit_kidney1_csr<-csresidual.coxph(fit = fit_kidney1,
                                  traindata = kidney,
                                  newdata = kidney)
fit_kidney1_loocvcsr<-cvcsresidual.coxph(fit = fit_kidney1,
                                         data=kidney,
                                         nfolds = nrow(kidney))

km.ucs.t.loocv <- survfit(Surv(fit_kidney1_loocvcsr,kidney$status)~1,
                          type='fleming')
id.ucs.t.loocv<-order(fit_kidney1_loocvcsr)

km.ucs.t <- survfit(Surv(fit_kidney1_csr,kidney$status)~1,type='fleming')
id.ucs.t<-order(fit_kidney1_csr)

pdf("kidneyplot_no_outlier.pdf",width=12,height=6)
par(mfrow = c(2,4))
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0))
ymax1=max(max(range(fit_kidney1_qr,fit_kidney1_loocvqr)),3)
ymin1=min(min(range(fit_kidney1_qr,fit_kidney1_loocvqr)),-3)
xmax1=max(range(sw_kidney_qr,sw_kidney_loocvqr,
                sw_qr_no_outlier,sw_loocvqr_no_outlier))
xmin1=min(range(sw_kidney_qr,sw_kidney_loocvqr,
                sw_qr_no_outlier,sw_loocvqr_no_outlier))
ymax2=max(range(km.ucs.t$cumhaz,km.ucs.t.loocv$cumhaz))
ymin2=min(range(km.ucs.t$cumhaz,km.ucs.t.loocv$cumhaz))
xmax2=max(range(km.ucs.t$time,km.ucs.t.loocv$time))
xmin2=min(range(km.ucs.t$time,km.ucs.t.loocv$time))
####the first row: No CV method
plot(fit_kidney1_qr,ylab="Z Residual plot",
     col=c("blue","darkolivegreen4")[kidney$status+1],
     pch=c(3,2)[kidney$status+1],
     main="Z Residual plot, No CV",
     ylim=c(ymin1,ymax1))
abline(h=c(3,-3),col="grey")
symbols(which(abs(fit_kidney1_loocvqr)>3),
        (fit_kidney1_qr)[which(abs(fit_kidney1_loocvqr)>3)],
        circles=2,fg='red',add=T, inches=F)
text(which(abs(fit_kidney1_loocvqr)>3),
     (fit_kidney1_qr)[which(abs(fit_kidney1_loocvqr)>3)],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 1,col="red")
legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)
qqnorm(fit_kidney1_qr,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("QQ plot, No CV, SW p-value = ", 
                   sprintf("%3.2f",sw_nrsp_qr)),
       ylim=c(ymin1,ymax1))
qqline(fit_kidney1_qr)
symbols((qqnorm(fit_kidney1_qr,plot.it=FALSE)
         $x)[which(abs(fit_kidney1_loocvqr)>3)],
        (qqnorm(fit_kidney1_qr,plot.it=FALSE)
         $y)[which(abs(fit_kidney1_loocvqr)>3)],
        circles=0.1,fg='red',add=T, inches=F)
text((qqnorm(fit_kidney1_qr,plot.it=FALSE)
      $x)[which(abs(fit_kidney1_loocvqr)>3)],
     (qqnorm(fit_kidney1_qr,plot.it=FALSE)
      $y)[which(abs(fit_kidney1_loocvqr)>3)],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 0.8,col="red")

hist(sw_qr_no_outlier,main="Replicated SW p-values, No CV",
     xlab="SW p-values",xlim=c(xmin1,xmax1))
plot(km.ucs.t, fun="cumhaz", xlab="Unmodified Cox-Snell Residuals",
     ylab="Cumulative Hazard Function",
     ylim=c(ymin2,ymax2),xlim=c(xmin2,xmax2),
     main="Cum. Hazard of CS Residuals, No CV")
abline(0, 1, col="chocolate4", lty=2)
points(km.ucs.t$time, -log(km.ucs.t$surv), 
       col=c("blue","darkolivegreen4")[kidney$status[id.ucs.t]+1],
       pch=c(3,2)[kidney$status[id.ucs.t]+1])
symbols(km.ucs.t$time[which(id.ucs.t.loocv==15)],
        (-log(km.ucs.t$surv))[which(id.ucs.t.loocv==15)],
        circles=0.2,fg='red',add=T, inches=F)
text(km.ucs.t$time[which(id.ucs.t.loocv==15)],
     (-log(km.ucs.t$surv))[which(id.ucs.t.loocv==15)],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 1,col="red")
legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)

####the second row: LOOCV method########
plot(fit_kidney1_loocvqr,ylab="Z Residual plot",
     col=c("blue","darkolivegreen4")[kidney$status+1],
     pch=c(3,2)[kidney$status+1],
     main="Z Residual plot, LOOCV",
     ylim=c(ymin1,ymax1))
abline(h=c(3,-3),col="grey")
symbols(which(abs(fit_kidney1_loocvqr)>3),
        fit_kidney1_loocvqr[which(abs(fit_kidney1_loocvqr)>3)],
        circles=2,fg='red',add=T, inches=F)
text(which(abs(fit_kidney1_loocvqr)>3),
     fit_kidney1_loocvqr[which(abs(fit_kidney1_loocvqr)>3)],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 1,col="red")

legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)
qqnorm(fit_kidney1_loocvqr,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("QQ plot, LOOCV, SW p-value = ", 
                   sprintf("%3.2f",sw_nrsp_loocvqr)),
       ylim=c(ymin1,ymax1))
qqline(fit_kidney1_loocvqr)
symbols((qqnorm(fit_kidney1_loocvqr,plot.it=FALSE)
         $x)[which(abs(fit_kidney1_loocvqr)>3)],
        (qqnorm(fit_kidney1_loocvqr,plot.it=FALSE)
         $y)[which(abs(fit_kidney1_loocvqr)>3)],
        circles=0.1,fg='red', add=T, inches=F)
text((qqnorm(fit_kidney1_loocvqr,plot.it=FALSE)
      $x)[which(abs(fit_kidney1_loocvqr)>3)],
     (qqnorm(fit_kidney1_loocvqr,plot.it=FALSE)
      $y)[which(abs(fit_kidney1_loocvqr)>3)],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 0.8,col="red")
hist(sw_loocvqr_no_outlier,main="Replicated SW p-values, LOOCV",
     xlab="SW p-values",xlim=c(xmin1,xmax1))
plot(km.ucs.t.loocv, fun="cumhaz", xlab="Unmodified Cox-Snell Residuals",
     ylab="Cumulative Hazard Function",
     ylim=c(ymin2,ymax2),xlim=c(xmin2,xmax2),
     main="Cum. Hazard of CS Residuals, LOOCV")
abline(0, 1, col="chocolate4", lty=2)
points(km.ucs.t.loocv$time, -log(km.ucs.t.loocv$surv), 
       col=c("blue","darkolivegreen4")[kidney$status[id.ucs.t.loocv]+1],
       pch=c(3,2)[kidney$status[id.ucs.t.loocv]+1])
symbols(km.ucs.t.loocv$time[which(id.ucs.t.loocv==15)],
        (-log(km.ucs.t.loocv$surv))[which(id.ucs.t.loocv==15)],
        circles=0.2,fg='red',add=T, inches=F)
text(km.ucs.t.loocv$time[which(id.ucs.t.loocv==15)],
     (-log(km.ucs.t.loocv$surv))[which(id.ucs.t.loocv==15)],
     pos=1,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 1,col="red")
legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)
dev.off()

data(kidney)
kidney<-kidney[-c(15,20,42),]
kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
kidney$id<-as.factor(kidney$id)
fit_kidney1 <- tryCatch(
  coxph(Surv(time, status) ~ age + sex + disease+
          frailty(id, distribution="gamma"), data= kidney),
  error = function(e) NA,
  warning = function(w) NA
)
AIC(fit_kidney1)
fix_var<- model.matrix(~.,data = kidney[,c(4,5,6),drop=FALSE])[,-1,drop=FALSE]
risk_score<-fix_var %*% fit_kidney1$coefficients+fit_kidney1$frail[kidney$id]

fit_kidney1_qr<-qresidual.coxph(fit_coxph = fit_kidney1,
                                traindata = kidney,
                                newdata = kidney)
fit_kidney1_loocvqr<-cvqresidual.coxph(fit = fit_kidney1,
                                       data=kidney,
                                       nfolds = nrow(kidney))

a1<-which(is.infinite(fit_kidney1_loocvqr))
a2<-which(is.na(fit_kidney1_loocvqr))

fit_kidney1_loocvqr = fit_kidney1_loocvqr[is.finite(fit_kidney1_loocvqr)]
fit_kidney1_loocvqr = fit_kidney1_loocvqr[!is.na(fit_kidney1_loocvqr)]

sw_nrsp_qr<-shapiro.test(fit_kidney1_qr)$p.value;sw_nrsp_qr
sw_nrsp_loocvqr<-shapiro.test(fit_kidney1_loocvqr)$p.value;sw_nrsp_loocvqr

fit_kidney1_csr<-csresidual.coxph(fit = fit_kidney1,
                                  traindata = kidney,
                                  newdata = kidney)
fit_kidney1_loocvcsr<-cvcsresidual.coxph(fit = fit_kidney1,
                                         data=kidney,
                                         nfolds = nrow(kidney))

km.ucs.t.loocv <- survfit(Surv(fit_kidney1_loocvcsr,kidney$status)~1,
                          type='fleming')
id.ucs.t.loocv<-order(fit_kidney1_loocvcsr)

km.ucs.t <- survfit(Surv(fit_kidney1_csr,kidney$status)~1,type='fleming')
id.ucs.t<-order(fit_kidney1_csr)

pdf("kidneyplot_no_outlier2.pdf",width=12,height=6)
par(mfrow = c(2,4))
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0))
ymax1=max(max(range(fit_kidney1_qr,fit_kidney1_loocvqr)),3)
ymin1=min(min(range(fit_kidney1_qr,fit_kidney1_loocvqr)),-3)
xmax1=max(range(sw_kidney_qr,sw_kidney_loocvqr,
                sw_qr_no_outlier,sw_loocvqr_no_outlier))
xmin1=min(range(sw_kidney_qr,sw_kidney_loocvqr,
                sw_qr_no_outlier,sw_loocvqr_no_outlier))
ymax2=max(range(km.ucs.t$cumhaz,km.ucs.t.loocv$cumhaz))
ymin2=min(range(km.ucs.t$cumhaz,km.ucs.t.loocv$cumhaz))
xmax2=max(range(km.ucs.t$time,km.ucs.t.loocv$time))
xmin2=min(range(km.ucs.t$time,km.ucs.t.loocv$time))
####the first row: No CV method
plot(fit_kidney1_qr,ylab="Z Residual plot",
     col=c("blue","darkolivegreen4")[kidney$status+1],
     pch=c(3,2)[kidney$status+1],
     main="Z Residual plot, No CV",
     ylim=c(ymin1,ymax1))
abline(h=c(3,-3),col="grey")
legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)
qqnorm(fit_kidney1_qr,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("QQ plot, No CV, SW p-value = ", 
                   sprintf("%3.2f",sw_nrsp_qr)),
       ylim=c(ymin1,ymax1))
qqline(fit_kidney1_qr)
hist(sw_qr_no_outlier1,main="Replicated SW p-values, No CV",
     xlab="SW p-values",xlim=c(xmin1,xmax1))
plot(km.ucs.t, fun="cumhaz", xlab="Unmodified Cox-Snell Residuals",
     ylab="Cumulative Hazard Function",
     ylim=c(ymin2,ymax2),xlim=c(xmin2,xmax2),
     main="Cum. Hazard of CS Residuals, No CV")
abline(0, 1, col="chocolate4", lty=2)
points(km.ucs.t$time, -log(km.ucs.t$surv), 
       col=c("blue","darkolivegreen4")[kidney$status[id.ucs.t]+1],
       pch=c(3,2)[kidney$status[id.ucs.t]+1])
legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)

####the second row: LOOCV method########
plot(fit_kidney1_loocvqr,ylab="Z Residual plot",
     col=c("blue","darkolivegreen4")[kidney$status+1],
     pch=c(3,2)[kidney$status+1],
     main="Z Residual plot, LOOCV",
     ylim=c(ymin1,ymax1))
abline(h=c(3,-3),col="grey")

legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)
qqnorm(fit_kidney1_loocvqr,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       main=paste0("QQ plot, LOOCV, SW p-value = ", 
                   sprintf("%3.2f",sw_nrsp_loocvqr)),
       ylim=c(ymin1,ymax1))
qqline(fit_kidney1_loocvqr)
hist(sw_loocvqr_no_outlier1,main="Replicated SW p-values, LOOCV",
     xlab="SW p-values",xlim=c(xmin1,xmax1))
plot(km.ucs.t.loocv, fun="cumhaz", xlab="Unmodified Cox-Snell Residuals",
     ylab="Cumulative Hazard Function",
     ylim=c(ymin2,ymax2),xlim=c(xmin2,xmax2),
     main="Cum. Hazard of CS Residuals, LOOCV")
abline(0, 1, col="chocolate4", lty=2)
points(km.ucs.t.loocv$time, -log(km.ucs.t.loocv$surv), 
       col=c("blue","darkolivegreen4")[kidney$status[id.ucs.t.loocv]+1],
       pch=c(3,2)[kidney$status[id.ucs.t.loocv]+1])
legend (x = "bottom", inset = c(0, -0.18),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)
dev.off()

data(kidney)
kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
kidney$id<-as.factor(kidney$id)
pdf("kidney_timeplot.pdf",width=5,height=5)
plot(kidney$time,ylab="Infection Time",
     col=c("blue","darkolivegreen4")[kidney$status+1],
     pch=c(3,2)[kidney$status+1],
     main="Infection time plot")
symbols(c(20,15,42),
        kidney$time[c(20,15,42)],
        circles=c(1.5, 1.5,1.5),fg=c('red', 'red', 'red'),
        add=T, inches=F)
text(c(20,15,42),
     kidney$time[c(20,15,42)],
     pos=1,label = c(20,15,42),
     cex = 1,col="red")
legend (x = "bottom", inset = c(0, -0.22),
        col=c("darkolivegreen4","blue"),
        pch=c(2,3),bty = "n",
        legend = c("Uncensored", "Censored"),
        horiz=TRUE,cex=1,xpd = TRUE)
dev.off()

plot(risk_score,fit_kidney1_qr,
     xlab="linear predictor",ylab="Z Residual plot",
     col=kidney$status+2,pch=kidney$status+1,
     main="Z Residual plot, No CV")
abline(h=c(3,-3))
plot(risk_score[-c(a1,a2)],fit_kidney1_loocvqr,
     xlab="linear predictor",ylab="Z Residual plot",
     col=kidney$status+2,pch=kidney$status+1,
     main="Z Residual plot, LOOCV")
abline(h=c(3,-3))
text(fit_kidney1_loocvqr[which(abs(fit_kidney1_loocvqr)>3)],
     pos=4,label = which(abs(fit_kidney1_loocvqr)>3),
     cex = 0.8,offset = 4.5)

plot(kidney$age,fit_kidney1_qr,
     xlab="age",ylab="Z Residual plot",
     col=kidney$status+2,pch=kidney$status+1,
     main="Z Residual plot, No CV")

plot(kidney$age[-c(a1,a2)],fit_kidney1_loocvqr,
     xlab="age",ylab="Z Residual plot",
     col=kidney$status+2,pch=kidney$status+1,
     main="Z Residual plot, LOOCV")

