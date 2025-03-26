setwd("~/Z-residual for Diagnosing Cox Proportional Hazard Models/")
source("~/Z-residual coxph.R")

library("survival")
library("EnvStats")
############################################################
data(pbc, package="survival")
pbc2<-pbc
pbc2<-subset(pbc, status != 1)
pbc2$status<-ifelse(pbc2$status==2,1,0)
pbc2$edema<-as.factor(pbc2$edema)

############ fit Cox PH linear and nonlinear models
fit_pbc2  <- tryCatch(
  coxph(Surv(time, status) ~ age + edema +
          (bili) + albumin, data= pbc2),
  error = function(e) NA,
  warning = function(w) NA
);
AIC(fit_pbc2)

fit_pbc2_nl  <- tryCatch(
  coxph(Surv(time, status) ~ age + edema +
          log(bili) + albumin, data= pbc2),
  error = function(e) NA,
  warning = function(w) NA
);
AIC(fit_pbc2_nl)

######score test
coxzph_test_pbc2<-cox.zph(fit_pbc2, transform="identity");coxzph_test_pbc2
coxzph_test_pbc2_NL<-cox.zph(fit_pbc2_nl, transform="identity");coxzph_test_pbc2_NL
##### Z-residuals test
coxph_pbc<-resid_coxph (coxfit_fit = fit_pbc2)
coxph_pbc_nl<-resid_coxph (coxfit_fit = fit_pbc2_nl)

#########sw test##############
sw_coxph_pbc<-shapiro.test(coxph_pbc$nrsp)$p.value;sw_coxph_pbc
sw_coxph_pbc_nl<-shapiro.test(coxph_pbc_nl$nrsp)$p.value;sw_coxph_pbc_nl

censored_pbc<- pbc2$status ==0
gof_censored_pbc<-gofTestCensored(coxph_pbc$nusp,censored_pbc, test = "sf", 
                                  censoring.side = "right",
                                  distribution = "norm")$p.value;gof_censored_pbc
censored_pbc_nl<- pbc2$status ==0
gof_censored_pbc_nl<-gofTestCensored(coxph_pbc_nl$nusp,censored_pbc_nl, test = "sf", 
                                     censoring.side = "right",
                                     distribution = "norm")$p.value;gof_censored_pbc_nl

#####anova test#########
anov_lp_pbc<- test.nl.aov(Zresidual=coxph_pbc$nrsp,
                          fitted.value=fit_pbc2$linear.predictors);anov_lp_pbc
anov_lp_pbc_nl<- test.nl.aov(Zresidual=coxph_pbc_nl$nrsp,
                             fitted.value=fit_pbc2_nl$linear.predictors);anov_lp_pbc_nl

bl_lp_pbc<-test.var.bartl(Zresidual=coxph_pbc$nrsp,
                          fitted.value = fit_pbc2$linear.predictors);bl_lp_pbc

bl_lp_pbc_nl<-test.var.bartl(Zresidual=coxph_pbc_nl$nrsp,
                             fitted.value = fit_pbc2_nl$linear.predictors);bl_lp_pbc_nl

###############for each covariate
anov_age<- test.nl.aov(Zresidual=coxph_pbc$nrsp,
                       fitted.value=pbc2$age);anov_age
anov_edema<- test.nl.aov(Zresidual=coxph_pbc$nrsp,
                         fitted.value=pbc2$edema);anov_edema
anov_albumin<- test.nl.aov(Zresidual=coxph_pbc$nrsp,
                           fitted.value=pbc2$albumin);anov_albumin
anov_bili<- test.nl.aov(Zresidual=coxph_pbc$nrsp,
                        fitted.value=pbc2$bili);anov_bili
bl_age<- test.var.bartl(Zresidual=coxph_pbc$nrsp,
                        fitted.value=pbc2$age);bl_age
bl_edema<- test.var.bartl(Zresidual=coxph_pbc$nrsp,
                          fitted.value=pbc2$edema);bl_edema
bl_albumin<- test.var.bartl(Zresidual=coxph_pbc$nrsp,
                            fitted.value=pbc2$albumin);bl_albumin
bl_bili<- test.var.bartl(Zresidual=coxph_pbc$nrsp,
                         fitted.value=pbc2$bili);bl_bili
anov_age_nl<- test.nl.aov(Zresidual=coxph_pbc_nl$nrsp,
                          fitted.value=pbc2$age);anov_age_nl
anov_albumin_nl<- test.nl.aov(Zresidual=coxph_pbc_nl$nrsp,
                              fitted.value=pbc2$albumin);anov_albumin_nl
anov_bili_nl<- test.nl.aov(Zresidual=coxph_pbc_nl$nrsp,
                           fitted.value=log(pbc2$bili));anov_bili_nl

bl_age_nl<- test.var.bartl(Zresidual=coxph_pbc_nl$nrsp,
                           fitted.value=pbc2$age);bl_age_nl
bl_albumin_nl<- test.var.bartl(Zresidual=coxph_pbc_nl$nrsp,
                               fitted.value=pbc2$albumin);bl_albumin_nl
bl_bili_nl<- test.var.bartl(Zresidual=coxph_pbc_nl$nrsp,
                            fitted.value=log(pbc2$bili));bl_bili_nl

pdf("plot/pbc_overall_zresid_schoenfresid_plot.pdf",width=16,height=7)
par(mfrow = c(2,4),mar=c(4,4,5,4))
ymax=max(range(coxph_pbc$nrsp,coxph_pbc_nl$nrsp),3.5)
ymin=min(range(coxph_pbc$nrsp,coxph_pbc_nl$nrsp),-3.5)
###the first row
qqnorm(coxph_pbc$nrsp,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       ylim = c(ymin,ymax),
       main="Cox PH Linear Model, Z-Residual QQ Plot")
qqline(coxph_pbc$nrsp)
text(-1.5, 3, paste0("Z-SW p-value = ", sprintf("%3.2f",sw_coxph_pbc)))
legend(x = "top", inset = c(0,-0.36), text.font=1,
       legend = c("(a)"),cex=2,xpd = TRUE,bty="n")

plot(pbc2$bili,coxph_pbc$nrsp,
     xlab="bilirunbin",ylab="Z-Residual",
     col=c("blue","darkolivegreen4")[pbc2$status+1],
     pch=c(3,2)[pbc2$status+1],
     ylim = c(ymin,ymax),
     main="Cox PH Linear Model, Z-Residual Scatterplot")
abline(h=c(3,-3),col="grey")
lines(lowess(coxph_pbc$nrsp ~ pbc2$bili),
      col = "red",lwd = 3)
legend(x = "bottom", inset = c(0, -0.38),
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")
legend(x = "top", inset = c(0,-0.36), text.font=1,
       legend = c("(b)"),cex=2,xpd = TRUE,bty="n")

bili.bin <- cut(pbc2$bili, 10)
plot(bili.bin,coxph_pbc$nrsp,xlab="bilirunbin",ylab="Z-Residual",
     main="Cox PH Linear Model, Z-Residual Boxplot",ylim=c(-4,5)) 
text(4.5,4.5, paste0("Z-AOV-bilirunbin p-value < 0.01 "))
#, sprintf("%3.2f",anov_bili)))
text(4.5,3.7, paste0("Z-BL-bilirunbin p-value = ", sprintf("%3.2f",bl_bili)))
legend(x = "top", inset = c(0,-0.36), text.font=1,
       legend = c("(c)"),cex=2,xpd = TRUE,bty="n")

plot(coxzph_test_pbc2, var = 3,main="Cox PH Linear Model, Schoenfeld Residual Plot",xlab="Event Time",
     ylim=c(-0.5,1))
text(2000, 0.95, paste0("Score Tests bilirunbin p-value < 0.01"))
#, sprintf("%3.2f",coxzph_test_pbc2$table[3,3])))
legend(x = "top", inset = c(0,-0.36), text.font=1,
       legend = c("(d)"),cex=2,xpd = TRUE,bty="n")

###the second row
qqnorm(coxph_pbc_nl$nrsp,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",
       ylim = c(ymin,ymax),
       main="Cox PH Log Model, Z-Residual QQ plot")
qqline(coxph_pbc_nl$nrsp)
text(-1.5, 3, paste0("Z-SW p-value = ", sprintf("%3.2f",sw_coxph_pbc_nl)))
legend(x = "top", inset = c(0,-0.36), text.font=1,
       legend = c("(e)"),cex=2,xpd = TRUE,bty="n")


plot(log(pbc2$bili),coxph_pbc_nl$nrsp,ylab="Z-Residual",
     xlab="log(bilirunbin)",
     col=c("blue","darkolivegreen4")[pbc2$status+1],
     pch=c(3,2)[pbc2$status+1],
     ylim = c(ymin,ymax),
     main="Cox PH Log Model, Z-Residual Scatterplot")
abline(h=c(3,-3),col="grey")
lines(lowess(coxph_pbc_nl$nrsp ~ log(pbc2$bili)),
      col = "red",lwd = 3)
legend(x = "bottom", inset = c(0, -0.36),
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),horiz=TRUE,cex=1,xpd = TRUE,bty="n")
legend(x = "top", inset = c(0,-0.33), text.font=1,
       legend = c("(f)"),cex=2,xpd = TRUE,bty="n")


logbili.bin <- cut(log(pbc2$bili), 10)
plot(logbili.bin,coxph_pbc_nl$nrsp,xlab="log(bilirunbin)",ylab="Z-Residual",
     main="Cox PH Log Model, Z-Residual Boxplot",ylim=c(-4,5)) 
text(5,4.5, paste0("Z-AOV-log(bilirunbin) p-value = ", sprintf("%3.2f",anov_bili_nl)))
text(5,3.7, paste0("Z-BL-log(bilirunbin) p-value = ", sprintf("%3.2f",bl_bili_nl)))
legend(x = "top", inset = c(0,-0.36), text.font=1,
       legend = c("(g)"),cex=2,xpd = TRUE,bty="n")

plot(coxzph_test_pbc2_NL, var = 3,main="Cox PH Log Model, Schoenfeld Residual Plot",xlab="Event Time",
     ylim=c(-2,4))
text(2000, 3.8, paste0("Score Tests log(bilirunbin) p-value = ", sprintf("%3.2f",coxzph_test_pbc2_NL$table[3,3])))
legend(x = "top", inset = c(0,-0.36), text.font=1,
       legend = c("(h)"),cex=2,xpd = TRUE,bty="n")

dev.off()

n_sims<-1000
cur_time = proc.time()
sw_coxph_rep<- rep(0,n_sims)
anov_bili_rep<- rep(0,n_sims)
bl_bili_rep<- rep(0,n_sims)
sw_coxph_nl_rep<- rep(0,n_sims)
anov_bili_nl_rep<- rep(0,n_sims)
bl_bili_nl_rep<- rep(0,n_sims)
for(j in 1:n_sims ){
  cat(paste('Simulation ',j,' out of ',n_sims,'\n'))
  if(j ==2){
    elapsed=as.numeric(proc.time()-cur_time)[3]
    cat(paste("Time for 1 simulation: ",elapsed/3600," hours \n"))
    cat(paste("Estimated time remaining: ",elapsed/3600*(n_sims-1)," hours \n"))
  }
  coxph_pbc<-resid_coxph (coxfit_fit = fit_pbc2)
  sw_coxph_rep[j]<-shapiro.test(coxph_pbc$nrsp)$p.value
  anov_bili_rep[j]<-test.nl.aov(Zresidual=coxph_pbc$nrsp,
                                fitted.value=pbc2$bili)
  bl_bili_rep[j]<-test.var.bartl(Zresidual=coxph_pbc$nrsp,
                                 fitted.value=pbc2$bili)
  #####################################################################
  coxph_pbc_nl<-resid_coxph (coxfit_fit = fit_pbc2_nl)
  
  sw_coxph_nl_rep[j]<-shapiro.test(coxph_pbc_nl$nrsp)$p.value
  
  anov_bili_nl_rep[j]<-test.nl.aov(Zresidual=coxph_pbc_nl$nrsp,
                                   fitted.value=log(pbc2$bili))
  bl_bili_nl_rep[j]<-test.var.bartl(Zresidual=coxph_pbc_nl$nrsp,
                                    fitted.value=log(pbc2$bili))
}
mean(sw_coxph_rep < 0.05)
mean(anov_bili_rep< 0.05)
mean(bl_bili_rep< 0.05)
mean(sw_coxph_nl_rep< 0.05)
mean(anov_bili_nl_rep< 0.05)
mean(bl_bili_nl_rep< 0.05)

bounds_pvalues <- function (pv)
{
  pv <- pv[is.finite(pv)]
  n <- length (pv)
  if(n<=0) stop("There is no value in 'pv'")
  pv.sorted <- sort(pv)
  corrected.pv <- pmin(1, pv.sorted*n/1:n)
  pv.min <- min(corrected.pv)
  pv.min
}

pmin_sw_pbc<-bounds_pvalues(pv=sw_coxph_rep);pmin_sw_pbc
pmin_anov_bili<-bounds_pvalues(pv=anov_bili_rep);pmin_anov_bili
pmin_bl_bili<-bounds_pvalues(pv=bl_bili_rep);pmin_bl_bili
pmin_sw_pbc_nl<-bounds_pvalues(pv=sw_coxph_nl_rep);pmin_sw_pbc_nl
pmin_anov_bili_nl<-bounds_pvalues(pv=anov_bili_nl_rep);pmin_anov_bili_nl
pmin_bl_bili_nl<-bounds_pvalues(pv=bl_bili_nl_rep);pmin_bl_bili_nl

pdf("plot/pbc_hist_overall.pdf",width=10,height=12)
par(mfrow = c(3,2),mar=c(4,4,5,2))
xmax1=max(range(sw_coxph_rep,anov_bili_rep,bl_bili_rep,
                #anov_protime_rep,bl_protime_rep,anov_protime_nl_rep,bl_protime_nl_rep,
                sw_coxph_nl_rep,anov_bili_nl_rep,bl_bili_nl_rep))
xmin1=min(range(sw_coxph_rep,anov_bili_rep,bl_bili_rep,
                #anov_protime_rep,bl_protime_rep,anov_protime_nl_rep,bl_protime_nl_rep,
                sw_coxph_nl_rep,anov_bili_nl_rep,bl_bili_nl_rep))
hist(sw_coxph_rep,main="Replicated Z-SW P-values for Cox PH Linear Model",
     xlab="Z-SW P-values for Cox PH Linear Model",xlim=c(xmin1,xmax1))
abline(v=pmin_sw_pbc,col="red")
legend(x = "top", inset = c(0,-0.29), text.font=1,
       legend = c("(a)"),cex=2,xpd = TRUE,bty="n")


hist(sw_coxph_nl_rep,main="Replicated Z-SW P-values for Cox PH Nonlinear Model",breaks=20,
     xlab="Z-SW P-values for Cox PH Nonlinear Model",xlim=c(xmin1,xmax1))
abline(v=pmin_sw_pbc_nl,col="red")
legend(x = "top", inset = c(0,-0.29), text.font=1,
       legend = c("(b)"),cex=2,xpd = TRUE,bty="n")


hist(anov_bili_rep,main="Replicated Z-AOV-bili P-values for Cox PH Linear Model", breaks=20, 
     xlab="Z-AOV-bili P-values for Cox PH Linear Model",xlim=c(xmin1,0.01))
abline(v=pmin_anov_bili,col="red")
legend(x = "top", inset = c(0,-0.29), text.font=1,
       legend = c("(c)"),cex=2,xpd = TRUE,bty="n")

hist(anov_bili_nl_rep,main="Replicated Z-AOV-log(bili) P-values for Cox PH Nonlinear Model",breaks=20,
     xlab="Z-AOV-log(bili) P-values for Cox PH Nonlinear Model",xlim=c(xmin1,xmax1))
abline(v=pmin_anov_bili_nl,col="red")
legend(x = "top", inset = c(0,-0.29), text.font=1,
       legend = c("(d)"),cex=2,xpd = TRUE,bty="n")


hist(bl_bili_rep,main="Replicated Z-BL-bili P-values for Cox PH Linear Model",breaks=20,
     xlab="Z-BL-bili P-values for Cox PH Linear Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_bili,col="red")
legend(x = "top", inset = c(0,-0.29), text.font=1,
       legend = c("(e)"),cex=2,xpd = TRUE,bty="n")

hist(bl_bili_nl_rep,main="Replicated Z-BL-log(bili) P-values for Cox PH Nonlinear Model",breaks=20,
     xlab="Z-BL-log(bili) P-values for Cox PH Nonlinear Model",xlim=c(xmin1,xmax1))
abline(v=pmin_bl_bili_nl,col="red")
legend(x = "top", inset = c(0,-0.29), text.font=1,
       legend = c("(f)"),cex=2,xpd = TRUE,bty="n")

dev.off()

