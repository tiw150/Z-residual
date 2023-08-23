setwd("~/Z-Residual Diagnostics for Detecting Non-Proportional Hazards/")
source("~/qresidual_coxph.R")
source("~/allresiduals.R")
source("~/Z_residual_AFT.R")

library("survival")
library("EnvStats")
library("nortest")
library("backports")
library("stringr")
library("timereg")
library("MASS")
library("dvmisc")
library("frailtySurv")
############################################################
data("drs")
drs$subject_id<-as.factor(drs$subject_id)
drs$treated<-as.factor(drs$treated)
drs$laser_type<-as.factor(drs$laser_type)
drs$diabetes_type<-as.factor(drs$diabetes_type)
#colnames(drs)<-c("subject","eye","time","status","Treated","Age","Laser Type", "Diabetes Type")
############ fit Cox PH linear model
fit_drs1  <- tryCatch(
  coxph(Surv(time, status) ~ treated + age_at_onset +
          laser_type+ diabetes_type+
          frailty(subject_id, distribution="gamma"), data= drs),
  error = function(e) NA,
  warning = function(w) NA
);fit_drs1

AIC(fit_drs1)

coxzph_test_drs<-cox.zph(fit_drs1, transform="identity");coxzph_test_drs

pdf( "plot/drs_scaledsch_plot.pdf",width=8,height=8)
par(mfrow = c(2,2),mar=c(4,4,2,2))
plot(coxzph_test_drs, var = 1,main="Cox PH Model",xlab="Event Time",ylim=c(-4,3),
     ylab="Beta(t) for Treated")
text(24, 2.8, paste0("Score Tests Treated p-value = ", sprintf("%3.2f",coxzph_test_drs$table[1,3])))
plot(coxzph_test_drs, var = 2,main="Cox PH Model",xlab="Event Time",ylim=c(-0.3,0.4),
     ylab="Beta(t) for Age")
text(24, 0.38, paste0("Score Tests Age p-value = ", sprintf("%3.2f",coxzph_test_drs$table[2,3])))
plot(coxzph_test_drs, var = 3,main="Cox PH Model",xlab="Event Time",ylim=c(-3,3.5),
     ylab="Beta(t) for Laser Type")
text(24, 3.3, paste0("Score Tests Laser p-value = ", sprintf("%3.2f",coxzph_test_drs$table[3,3])))
plot(coxzph_test_drs, var = 4,main="Cox PH Model",xlab="Event Time",ylim=c(-6,9),
     ylab="Beta(t) for Diabetes Type")
text(24, 8.8, paste0("Score Tests Diabetes p-value = ", sprintf("%3.2f",coxzph_test_drs$table[4,3])))
dev.off()

qr_drs1<-qresidual.coxph (fit_coxph = fit_drs1,
                          traindata = drs,newdata = drs)
resids_drs1<-allresidual.coxph (fit_coxph = fit_drs1,
                                traindata = drs, newdata = drs)
censored_cox<- drs$status ==0
gof_fit_drs1<-gofTestCensored(resids_drs1$nusp,censored_cox, test = "sf", 
                                 censoring.side = "right",
                                 distribution = "norm")$p.value
gof_fit_drs1
sw_qr<-shapiro.test(qr_drs1)$p.value;sw_qr
anov_drs<-test.nl.aov1(qresidual=qr_drs1, 
                        fitted.values=fit_drs1$linear.predictors,
                        k.anova=10);anov_drs
bl_drs<-test.var.bartl(qresidual=qr_drs1,
                       fitted.values=fit_drs1$linear.predictors);bl_drs
anov_qr_age<-test.nl.aov1(qresidual=qr_drs1, 
                          fitted.values=drs$age_at_onset,
                          k.anova=10);anov_qr_age
bl_qr_age<-test.var.bartl(qresidual=qr_drs1, 
                          fitted.values=drs$age_at_onset,
                          k=10);bl_qr_age
anov_qr_laser<-test.nl.aov1.categ(qresidual=qr_drs1, 
                                  fitted.values=drs$laser_type);anov_qr_laser
bl_qr_laser<-test.var.bartl.categ(qresidual=qr_drs1, 
                                  fitted.values=drs$laser_type);bl_qr_laser
anov_qr_diabete<-test.nl.aov1.categ(qresidual=qr_drs1, 
                                  fitted.values=drs$diabetes_type);anov_qr_diabete
bl_qr_diabete<-test.var.bartl.categ(qresidual=qr_drs1, 
                                  fitted.values=drs$diabetes_type);bl_qr_diabete
anov_qr_trt<-test.nl.aov1.categ(qresidual=qr_drs1, 
                                    fitted.values=drs$treated);anov_qr_trt
bl_qr_trt<-test.var.bartl.categ(qresidual=qr_drs1, 
                                    fitted.values=drs$treated);bl_qr_trt


fit_timecox_drs1 <- tryCatch(
  timecox(Surv(time, status)~treated + age_at_onset +
            laser_type+ diabetes_type+ cluster(frailty(subject_id, distribution="gamma")),
          residuals=1,data=drs,max.time=7,n.sim=100),
  error = function(e) NA,
  warning = function(w) NA
)

cumresid_test<-cum.residuals(fit_timecox_drs1,data=drs,cum.resid=1,n.sim=100)
cum_resid<-cumresid_test[["pval.test"]];cum_resid

#####fit AFT Lognormal model###############
fit_drs_survreg<-
  survreg(Surv(time, status) ~ treated + age_at_onset +
            laser_type+ diabetes_type, data= drs,dist="lognormal")

AIC(fit_drs_survreg)

summary(fit_drs_survreg)
qr_drs_aft<-(resid_survreg (survreg_fit = fit_drs_survreg))$nrsp
nusp_drs_aft<-(resid_survreg (survreg_fit = fit_drs_survreg))$nusp
sw_qr_aft<-shapiro.test(qr_drs_aft)$p.value;sw_qr_aft
censored_aft<- drs$status ==0
gof_fit_drs1_aft<-gofTestCensored(nusp_drs_aft,censored_aft, test = "sf", 
                                     censoring.side = "right",
                                     distribution = "norm")$p.value
gof_fit_drs1_aft
anov_drs_aft<-test.nl.aov1(qresidual=qr_drs_aft, 
                           fitted.values=fit_drs_survreg$linear.predictors,
                           k.anova=10);anov_drs_aft
bl_drs_aft<-test.var.bartl(qresidual=qr_drs_aft,
                           fitted.values=fit_drs_survreg$linear.predictors);bl_drs_aft
anov_qr_age_aft<-test.nl.aov1(qresidual=qr_drs_aft, 
                          fitted.values=drs$age_at_onset,
                          k.anova=10);anov_qr_age_aft
bl_qr_age_aft<-test.var.bartl(qresidual=qr_drs_aft, 
                          fitted.values=drs$age_at_onset,
                          k=10);bl_qr_age_aft
anov_qr_laser_aft<-test.nl.aov1.categ(qresidual=qr_drs_aft, 
                                  fitted.values=drs$laser_type);anov_qr_laser_aft
bl_qr_laser_aft<-test.var.bartl.categ(qresidual=qr_drs_aft, 
                                  fitted.values=drs$laser_type);bl_qr_laser_aft
anov_qr_diabete_aft<-test.nl.aov1.categ(qresidual=qr_drs_aft, 
                                    fitted.values=drs$diabetes_type);anov_qr_diabete_aft
bl_qr_diabete_aft<-test.var.bartl.categ(qresidual=qr_drs_aft, 
                                    fitted.values=drs$diabetes_type);bl_qr_diabete_aft
anov_qr_trt_aft<-test.nl.aov1.categ(qresidual=qr_drs_aft, 
                                fitted.values=drs$treated);anov_qr_trt_aft
bl_qr_trt_aft<-test.var.bartl.categ(qresidual=qr_drs_aft, 
                                fitted.values=drs$treated);bl_qr_trt_aft

pdf("plot/drs_overall_qq_resid_plot.pdf",width=7,height=9)
par(mfrow = c(3,2),mar=c(4,4,2,2))
ymax=max(range(qr_drs1,qr_drs_aft),3)
ymin=min(range(qr_drs1,qr_drs_aft),-3)
qqnorm(qr_drs1,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",ylim = c(ymin,ymax),
       main="Cox PH Model, Z-Residual QQ Plot")
qqline(qr_drs1)
text(-1.2, 2.8, paste0("Z-SW p-value = ", sprintf("%3.2f",sw_qr)))
qqnorm(qr_drs_aft,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",ylim = c(ymin,ymax),
       main="Lognormal AFT Model, QQ plot of Z-Residuals")
qqline(qr_drs_aft)
text(-1.2, 2.8, paste0("Z-SW p-value = ", sprintf("%3.2f",sw_qr_aft)))

plot(fit_drs1$linear.predictors,qr_drs1,
     xlab="Linear Predictor",ylab="Z-Residual",
     col=c("blue","darkolivegreen4")[drs$status+1],
     pch=c(3,2)[drs$status+1],ylim = c(ymin,ymax),
     main="Cox PH Model, Z-Residual Scatterplot")
abline(h=c(3,-3),col="grey")
lines(lowess(qr_drs1 ~ fit_drs1$linear.predictors),
      col = "red",lwd = 3)
legend(x = "bottom", inset = c(0, -0.2),
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")

plot(fit_drs_survreg$linear.predictors,
     qr_drs_aft,ylab="Z-Residual",
     xlab="Linear Predictor",
     col=c("blue","darkolivegreen4")[drs$status+1],
     pch=c(3,2)[drs$status+1],ylim = c(ymin,ymax),
     main="Lognormal AFT Model, Z-Residual Scatterplot")
abline(h=c(3,-3),col="grey")
lines(lowess(qr_drs_aft ~ fit_drs_survreg$linear.predictors),
      col = "red",lwd = 3)
legend(x = "bottom", inset = c(0, -0.2),
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")

lpred.bin <- cut(fit_drs1$linear.predictors, 10)
plot(lpred.bin,qr_drs1,xlab="Linear Predictor",ylab="Z-Residual",
     main="Cox PH Model, Z-Residual Boxplot",ylim=c(-3,5)) 
text(3.5,4.8, paste0("Z-AOV-LP p-value = ", sprintf("%3.2f",anov_drs)))
text(3.5,4.3, paste0("Z-BL-LP p-value = ", sprintf("%3.2f",bl_drs)))

lpred.bin_AFT <- cut(fit_drs_survreg$linear.predictors, 10)
plot(lpred.bin_AFT,qr_drs_aft,xlab="Linear Predictor",ylab="Z-Residual",
     main="Lognormal AFT Model, Z-Residual Boxplot",ylim=c(-3,5)) 
text(3.5,4.8, paste0("Z-AOV-LP p-value = ", sprintf("%3.2f",anov_drs_aft)))
text(3.5,4.3, paste0("Z-BL-LP p-value = ", sprintf("%3.2f",bl_drs_aft)))

dev.off()

pdf("plot/drs_resid_plot.pdf",width=7,height=15)
par(mfrow = c(5,2),mar=c(4,4,2,2))

ymax=max(range(qr_drs1,qr_drs_aft),3)
ymin=min(range(qr_drs1,qr_drs_aft),-3)

plot(drs$age_at_onset ,qr_drs1,ylab="Z-Residual",
     col=c("blue","darkolivegreen4")[drs$status+1],
     pch=c(3,2)[drs$status+1],xlab="Age",ylim = c(ymin,ymax),
     main="Cox PH Model, Z-Residual Scatterplot")
abline(h=c(3,-3),col="grey")
lines(lowess(qr_drs1 ~ drs$age_at_onset),
      col = "red",lwd = 3)
legend(x = "bottom", inset = c(0, -0.2),
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")

plot(drs$age_at_onset ,qr_drs_aft,ylab="Z-Residual",xlab="Age",
     col=c("blue","darkolivegreen4")[drs$status+1],
     pch=c(3,2)[drs$status+1],ylim = c(ymin,ymax),
     main="Lognormal AFT Model, Z-Residual Scatterplot")
abline(h=c(3,-3),col="grey")
lines(lowess(qr_drs_aft ~ drs$age_at_onset),
      col = "red",lwd = 3)
legend(x = "bottom", inset = c(0, -0.2),
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")

lpred.bin <- cut(drs$age_at_onset, 10)
plot(lpred.bin,qr_drs1,xlab="Age",ylab="Z-Residual",
     main="Cox PH Model, Z-Residual Boxplot",ylim=c(-3,4.5)) 
text(4,4.3, paste0("Z-AOV-Age p-value = ", sprintf("%3.2f",anov_qr_age)))
text(4,3.8, paste0("Z-BL-Age p-value = ", sprintf("%3.2f",bl_qr_age)))

lpred.bin <- cut(drs$age_at_onset, 10)
plot(lpred.bin,qr_drs_aft,xlab="Age",ylab="Z-Residual",
     main="Lognormal AFT Model, Z-Residual Boxplot",ylim=c(-3,4.5)) 
text(4,4.3, paste0("Z-AOV-Age p-value = ", sprintf("%3.2f",anov_qr_age_aft)))
text(4,3.8, paste0("Z-BL-Age p-value = ", sprintf("%3.2f",bl_qr_age_aft)))


plot(drs$treated,qr_drs1,xlab="Treated",ylab="Z-Residual",
     main="Cox PH Model, Z-Residual Boxplot",ylim=c(-3,5)) 
text(1.2,4.8, paste0("Z-AOV-Treated p-value = ", sprintf("%3.2f",anov_qr_trt)))
text(1.2,4.3, paste0("Z-BL-Treated p-value = ", sprintf("%3.2f",bl_qr_trt)))

plot(drs$treated,qr_drs_aft,xlab="Treated",ylab="Z-Residual",
     main="Lognormal AFT Model, Z-Residual Boxplot",ylim=c(-3,5)) 
text(1.2,4.8, paste0("Z-AOV-Treated p-value = ", sprintf("%3.2f",anov_qr_trt_aft)))
text(1.2,4.3, paste0("Z-BL-Treated p-value = ", sprintf("%3.2f",bl_qr_trt_aft)))

plot(drs$laser_type,qr_drs1,xlab="Laser Type",ylab="Z-Residual",
     main="Cox PH Model, Z-Residual Boxplot",ylim=c(-3,5)) 
text(1.2,4.8, paste0("Z-AOV-Laser p-value = ", sprintf("%3.2f",anov_qr_laser)))
text(1.2,4.3, paste0("Z-BL-Laser p-value = ", sprintf("%3.2f",bl_qr_laser)))
plot(drs$laser_type,qr_drs_aft,xlab="Laser Type",ylab="Z-Residual",
     main="Lognormal AFT Model, Z-Residual Boxplot",ylim=c(-3,5)) 
text(1.2,4.8, paste0("Z-AOV-Laser p-value = ", sprintf("%3.2f",anov_qr_laser_aft)))
text(1.2,4.3, paste0("Z-BL-Laser p-value = ", sprintf("%3.2f",bl_qr_laser_aft)))

plot(drs$diabetes_type,qr_drs1,xlab="Diabetes Type",ylab="Z-Residual",
     main="Cox PH Model, Z-Residual Boxplot",ylim=c(-3,5)) 
text(1.2,4.8, paste0("Z-AOV-Diabetes p-value = ", sprintf("%3.2f",anov_qr_diabete)))
text(1.2,4.3, paste0("Z-BL-Diabetes p-value = ", sprintf("%3.2f",bl_qr_diabete)))

plot(drs$diabetes_type,qr_drs1,xlab="Diabetes Type",ylab="Z-Residual",
     main="Lognormal AFT Model, Z-Residual Boxplot",ylim=c(-3,5)) 
text(1.2,4.8, paste0("Z-AOV-Diabetes p-value = ", sprintf("%3.2f",anov_qr_diabete_aft)))
text(1.2,4.3, paste0("Z-BL-Diabetes p-value = ", sprintf("%3.2f",bl_qr_diabete_aft)))

dev.off()



n_sims<-1000
cur_time = proc.time()
sw_fit_drs1<- rep(0,n_sims)
anov_fit_drs1<- rep(0,n_sims)
bl_fit_drs1<- rep(0,n_sims)
anov_qr_age<- rep(0,n_sims)
bl_qr_age<- rep(0,n_sims)
anov_qr_laser<- rep(0,n_sims)
bl_qr_laser<- rep(0,n_sims)
anov_qr_diabete<- rep(0,n_sims)
bl_qr_diabete<- rep(0,n_sims)
anov_qr_treat<- rep(0,n_sims)
bl_qr_treat<- rep(0,n_sims)
sw_fit_drs1_aft<- rep(0,n_sims)
anov_fit_drs1_aft<- rep(0,n_sims)
bl_fit_drs1_aft<- rep(0,n_sims)
anov_qr_laser_aft<- rep(0,n_sims)
bl_qr_laser_aft<- rep(0,n_sims)
anov_qr_diabete_aft<- rep(0,n_sims)
bl_qr_diabete_aft<- rep(0,n_sims)
anov_qr_treat_aft<- rep(0,n_sims)
bl_qr_treat_aft<- rep(0,n_sims)
anov_qr_age_aft<- rep(0,n_sims)
bl_qr_age_aft<- rep(0,n_sims)
for(j in 1:n_sims ){
  cat(paste('Simulation ',j,' out of ',n_sims,'\n'))
  if(j ==2){
    elapsed=as.numeric(proc.time()-cur_time)[3]
    cat(paste("Time for 1 simulation: ",elapsed/3600," hours \n"))
    cat(paste("Estimated time remaining: ",elapsed/3600*(n_sims-1)," hours \n"))
  }
  qr_drs1<-qresidual.coxph (fit_coxph = fit_drs1,
                            traindata = drs,newdata = drs)
  sw_fit_drs1[j]<-shapiro.test(qr_drs1)$p.value
  anov_fit_drs1[j]<-test.nl.aov1(qresidual=qr_drs1, 
                                fitted.values=fit_drs1$linear.predictors)
  bl_fit_drs1[j]<-test.var.bartl(qresidual=qr_drs1,
                             fitted.values=fit_drs1$linear.predictors)
  
  anov_qr_age[j]<-test.nl.aov1(qresidual=qr_drs1, 
                                 fitted.values=drs$age_at_onset)
  bl_qr_age[j]<-test.var.bartl(qresidual=qr_drs1,
                                 fitted.values=drs$age_at_onset)
  anov_qr_laser[j]<-test.nl.aov1.categ(qresidual=qr_drs1, 
                                       fitted.values=drs$laser_type)
  bl_qr_laser[j]<-test.var.bartl.categ(qresidual=qr_drs1, 
                                    fitted.values=drs$laser_type)
  anov_qr_diabete[j]<-test.nl.aov1.categ(qresidual=qr_drs1, 
                                      fitted.values=drs$diabetes_type)
  bl_qr_diabete[j]<-test.var.bartl.categ(qresidual=qr_drs1, 
                                      fitted.values=drs$diabetes_type)
  anov_qr_treat[j]<-test.nl.aov1.categ(qresidual=qr_drs1, 
                                      fitted.values=drs$treated)
  bl_qr_treat[j]<-test.var.bartl.categ(qresidual=qr_drs1, 
                                      fitted.values=drs$treated)
  
  #####################################################################
  qr_drs_aft<-(resid_survreg (survreg_fit = fit_drs_survreg))$nrsp
  sw_fit_drs1_aft[j]<-shapiro.test(qr_drs_aft)$p.value
  anov_fit_drs1_aft[j]<-test.nl.aov1(qresidual=qr_drs_aft, 
                                 fitted.values=fit_drs_survreg$linear.predictors)
  bl_fit_drs1_aft[j]<-test.var.bartl(qresidual=qr_drs_aft,
                                 fitted.values=fit_drs_survreg$linear.predictors)
  
  anov_qr_age_aft[j]<-test.nl.aov1.categ(qresidual=qr_drs_aft, 
                                           fitted.values=drs$age_at_onset)
  bl_qr_age_aft[j]<-test.var.bartl.categ(qresidual=qr_drs_aft, 
                                           fitted.values=drs$age_at_onset)
  anov_qr_laser_aft[j]<-test.nl.aov1.categ(qresidual=qr_drs_aft, 
                                       fitted.values=drs$laser_type)
  bl_qr_laser_aft[j]<-test.var.bartl.categ(qresidual=qr_drs_aft, 
                                       fitted.values=drs$laser_type)
  anov_qr_diabete_aft[j]<-test.nl.aov1.categ(qresidual=qr_drs_aft, 
                                         fitted.values=drs$diabetes_type)
  bl_qr_diabete_aft[j]<-test.var.bartl.categ(qresidual=qr_drs_aft, 
                                         fitted.values=drs$diabetes_type)
  anov_qr_treat_aft[j]<-test.nl.aov1.categ(qresidual=qr_drs_aft, 
                                       fitted.values=drs$treated)
  bl_qr_treat_aft[j]<-test.var.bartl.categ(qresidual=qr_drs_aft, 
                                       fitted.values=drs$treated)
  
}

mean(sw_fit_drs1 < 0.05)

mean(anov_fit_drs1< 0.05)

mean(bl_fit_drs1< 0.05)

mean(anov_qr_laser< 0.05)

mean(bl_qr_laser< 0.05)

mean(anov_qr_diabete< 0.05)

mean(bl_qr_diabete< 0.05)

mean(anov_qr_treat< 0.05)

mean(bl_qr_treat< 0.05)

mean(anov_qr_age< 0.05)

mean(bl_qr_age< 0.05)

mean(sw_fit_drs1_aft < 0.05)

mean(anov_fit_drs1_aft< 0.05)

mean(bl_fit_drs1_aft< 0.05)

mean(anov_qr_laser_aft< 0.05)

mean(bl_qr_laser_aft< 0.05)

mean(anov_qr_diabete_aft< 0.05)

mean(bl_qr_diabete_aft< 0.05)

mean(anov_qr_treat_aft< 0.05)

mean(bl_qr_treat_aft< 0.05)

mean(anov_qr_age_aft< 0.05)

mean(bl_qr_age_aft< 0.05)


pmin_sw_drs<-bounds_pvalues(pv=sw_fit_drs1);pmin_sw_drs

pmin_aov_lp_drs<-bounds_pvalues(pv=anov_fit_drs1);pmin_aov_lp_drs

pmin_bl_lp_drs<-bounds_pvalues(pv=bl_fit_drs1);pmin_bl_lp_drs

pmin_aov_age_drs<-bounds_pvalues(pv=anov_qr_laser);pmin_aov_age_drs

pmin_bl_age_drs<-bounds_pvalues(pv=bl_qr_laser);pmin_bl_age_drs

pmin_aov_laser_drs<-bounds_pvalues(pv=anov_qr_laser);pmin_aov_laser_drs

pmin_bl_laser_drs<-bounds_pvalues(pv=bl_qr_laser);pmin_bl_laser_drs

pmin_aov_diabete_drs<-bounds_pvalues(pv=anov_qr_diabete);pmin_aov_diabete_drs

pmin_bl_diabete_drs<-bounds_pvalues(pv=bl_qr_diabete);pmin_bl_diabete_drs

pmin_aov_treat_drs<-bounds_pvalues(pv=anov_qr_treat);pmin_aov_treat_drs

pmin_bl_treat_drs<-bounds_pvalues(pv=bl_qr_treat);pmin_bl_treat_drs


pmin_sw_drs_aft<-bounds_pvalues(pv=sw_fit_drs1_aft);pmin_sw_drs_aft

pmin_aov_lp_drs_aft<-bounds_pvalues(pv=anov_fit_drs1_aft);pmin_aov_lp_drs_aft

pmin_bl_lp_drs_aft<-bounds_pvalues(pv=bl_fit_drs1_aft);pmin_bl_lp_drs_aft

pmin_aov_age_drs_aft<-bounds_pvalues(pv=anov_qr_laser_aft);pmin_aov_age_drs_aft

pmin_bl_age_drs_aft<-bounds_pvalues(pv=bl_qr_laser_aft);pmin_bl_age_drs_aft

pmin_aov_laser_drs_aft<-bounds_pvalues(pv=anov_qr_laser_aft);pmin_aov_laser_drs_aft

pmin_bl_laser_drs_aft<-bounds_pvalues(pv=bl_qr_laser_aft);pmin_bl_laser_drs_aft

pmin_aov_diabete_drs_aft<-bounds_pvalues(pv=anov_qr_diabete_aft);pmin_aov_diabete_drs_aft

pmin_bl_diabete_drs_aft<-bounds_pvalues(pv=bl_qr_diabete_aft);pmin_bl_diabete_drs_aft

pmin_aov_treat_drs_aft<-bounds_pvalues(pv=anov_qr_treat_aft);pmin_aov_treat_drs_aft

pmin_bl_treat_drs_aft<-bounds_pvalues(pv=bl_qr_treat_aft);pmin_bl_treat_drs_aft



pdf("plot/drs_hist_overall.pdf",width=8,height=9)
par(mfrow = c(3,2),mar=c(4,4,2,2))
xmax1=max(range(sw_fit_drs1,anov_fit_drs1,bl_fit_drs1,
                sw_fit_drs1_aft,anov_fit_drs1_aft,bl_fit_drs1_aft))
xmin1=min(range(sw_fit_drs1,anov_fit_drs1,bl_fit_drs1,
                sw_fit_drs1_aft,anov_fit_drs1_aft,bl_fit_drs1_aft))
hist(sw_fit_drs1,main="Replicated Z-SW P-values for Cox PH",
     xlab="Z-SW P-values for Cox PH Model",xlim=c(xmin1,xmax1))
abline(v=pmin_sw_drs,col="red")

hist(sw_fit_drs1_aft,main="Replicated Z-SW P-values for Lognormal AFT",breaks=20,
     xlab="Z-SW P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_sw_drs_aft,col="red")

hist(anov_fit_drs1,main="Replicated Z-AOV-LP P-values for Cox PH", breaks=20, 
     xlab="Z-AOV-LP P-values for Cox PH Model",xlim=c(0,0.0000005))
abline(v=pmin_aov_lp_drs,col="red")
hist(anov_fit_drs1_aft,main="Replicated Z-AOV-LP P-values for Lognormal AFT",breaks=20,
     xlab="Z-AOV-LP P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_aov_lp_drs_aft,col="red")

hist(bl_fit_drs1,main="Replicated Z-BL-LP P-values for Cox PH",breaks=20,
     xlab="Z-BL-LP P-values for Cox PH Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_lp_drs,col="red")
hist(bl_fit_drs1_aft,main="Replicated Z-BL-LP P-values for Lognormal AFT",breaks=20,
     xlab="Z-BL-LP P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_lp_drs_aft,col="red")
dev.off()


pdf("plot/drs_hist_covariate1.pdf",width=8,height=12)
par(mfrow = c(4,2),mar=c(4,4,2,2))

hist(anov_qr_treat,main="Replicated Z-AOV-Treat P-values for Cox PH",breaks=20,
     xlab="Z-AOV-Treat P-values for Cox PH Model",xlim=c(xmin1,xmax1))
abline(v=pmin_aov_treat_drs,col="red")
hist(anov_qr_treat_aft,main="Replicated Z-AOV-Treat P-values for Lognormal AFT",breaks=20,
     xlab="Z-AOV-Treat P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_aov_treat_drs_aft,col="red")

hist(bl_qr_treat,main="Replicated Z-BL-Treat P-values for Cox PH",breaks=20,
     xlab="Z-BL-Treat P-values for Cox PH Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_treat_drs,col="red")
hist(bl_qr_treat_aft,main="Replicated Z-BL-Treat P-values for Lognormal AFT",breaks=20,
     xlab="Z-BL-Treat P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_treat_drs_aft,col="red")


hist(anov_qr_age,main="Replicated Z-AOV-Age P-values for Cox PH",breaks=20,
     xlab="Z-AOV-Age P-values for Cox PH Model",xlim=c(xmin1,xmax1))
abline(v=pmin_aov_age_drs,col="red")
hist(anov_qr_age_aft,main="Replicated Z-AOV-Age P-values for Lognormal AFT",breaks=20,
     xlab="Z-AOV-Age P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_aov_age_drs_aft,col="red")

hist(bl_qr_age,main="Replicated Z-BL-Age P-values for Cox PH",breaks=20,
     xlab="Z-BL-Age P-values for Cox PH Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_age_drs,col="red")
hist(bl_qr_age_aft,main="Replicated Z-BL-Age P-values for Lognormal AFT",breaks=20,
     xlab="Z-BL-Age P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_age_drs_aft,col="red")
dev.off()

pdf("plot/drs_hist_covariate2.pdf",width=8,height=12)
par(mfrow = c(4,2),mar=c(4,4,2,2))

hist(anov_qr_laser,main="Replicated Z-AOV-Laser P-values for Cox PH",breaks=20,
     xlab="Z-AOV-Laser P-values for Cox PH Model",xlim=c(xmin1,xmax1))
abline(v=pmin_aov_laser_drs,col="red")
hist(anov_qr_laser_aft,main="Replicated Z-AOV-Laser P-values for Lognormal AFT",breaks=20,
     xlab="Z-AOV-Laser P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_aov_laser_drs_aft,col="red")

hist(bl_qr_laser,main="Replicated Z-BL-Laser P-values for Cox PH",breaks=20,
     xlab="Z-BL-Laser P-values for Cox PH Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_laser_drs,col="red")
hist(bl_qr_laser_aft,main="Replicated Z-BL-Laser P-values for Lognormal AFT",breaks=20,
     xlab="Z-BL-Laser P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_laser_drs_aft,col="red")

hist(anov_qr_diabete,main="Replicated Z-AOV-Diabete P-values for Cox PH",breaks=20,
     xlab="Z-AOV-Diabete P-values for Cox PH Model",xlim=c(xmin1,xmax1))
abline(v=pmin_aov_diabete_drs,col="red")
hist(anov_qr_diabete_aft,main="Replicated Z-AOV-Diabete P-values for Lognormal AFT",breaks=20,
     xlab="Z-AOV-Diabete P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_aov_diabete_drs_aft,col="red")

hist(bl_qr_diabete,main="Replicated Z-BL-Diabete P-values for Cox PH",breaks=20,
     xlab="Z-BL-Diabete P-values for Cox PH Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_diabete_drs,col="red")
hist(bl_qr_diabete_aft,main="Replicated Z-BL-Diabete P-values for Lognormal AFT",breaks=20,
     xlab="Z-BL-Diabete P-values for Lognormal AFT Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_diabete_drs_aft,col="red")

dev.off()

