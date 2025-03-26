setwd("~/Z-residuals for Diagnosing Cox Proportional Hazard Models/")
source("~/nonlinear_data.R")
source("~/Z-residual_coxph.R")
library("survival")
library("EnvStats")
library("nortest")
library("backports")
library("stringr")
library("timereg")
library("MASS")
library("dvmisc")

###simulate data
nonph_data<-simul_nonph_linear(n_clusters=10,n_individuals=50,lambda=0.007, alpha=3, 
                               mean.censor=7.5,beta1=-0.5, beta_a=1.35,beta_b=-1.35, t0=pi/2,corr=0)
ph_data<-simul_ph_nonlinear(n_clusters=10,n_individuals=50,lambda=0.007, alpha=3, 
                            mean.censor=7.5,beta1=-0.5, beta_a=-2,beta_b=-2,t0=pi/2,corr=0)

#######fit coxph model################################################
fit_coxph_ph_data <- tryCatch(
  coxph(Surv(t, d) ~ x1+x2 ,
        data = ph_data),
  error = function(e) NA,
  warning = function(w) NA
)
fit_coxph_nonph_data <- tryCatch(
  coxph(Surv(t, d) ~ x1+x2,
        data = nonph_data),
  error = function(e) NA,
  warning = function(w) NA
)
#####Score test##########################################################
coxzph_test_ph<-cox.zph(fit_coxph_ph_data, transform="identity");coxzph_test_ph
coxzph_test_nonph<-cox.zph(fit_coxph_nonph_data, transform="identity");coxzph_test_nonph
plot(coxzph_test_ph, var = 1)
plot(coxzph_test_ph, var = 2)

plot(coxzph_test_nonph, var = 1)
plot(coxzph_test_nonph, var = 2)
##### Z-residuals test#########################################################
coxph_qr_ph<-resid_coxph (coxfit_fit = fit_coxph_ph_data)
coxph_qr_nonph<-resid_coxph (coxfit_fit = fit_coxph_nonph_data)

#########sw test##############
sw_coxph_ph<-shapiro.test(coxph_qr_ph$nrsp)$p.value;sw_coxph_ph
sw_coxph_nonph<-shapiro.test(coxph_qr_nonph$nrsp)$p.value;sw_coxph_nonph

censored_ph<- ph_data$d ==0
gof_censored_ph<-gofTestCensored(coxph_qr_ph$nusp,censored_ph, test = "sf", 
                                 censoring.side = "right",
                                 distribution = "norm")$p.value;gof_censored_ph
censored_nonph<- nonph_data$d ==0
gof_censored_nonph<-gofTestCensored(coxph_qr_nonph$nusp,censored_nonph, test = "sf", 
                                    censoring.side = "right",
                                    distribution = "norm")$p.value;gof_censored_nonph

#####anova test#########
anov_lp_ph<- test.nl.aov(Zresidual=coxph_qr_ph$nrsp,
                         fitted.value=fit_coxph_ph_data$linear.predictors);anov_lp_ph
anov_lp_nonph<- test.nl.aov(Zresidual=coxph_qr_nonph$nrsp,
                            fitted.value=fit_coxph_nonph_data$linear.predictors);anov_lp_nonph

bl_lp_ph<-test.var.bartl(Zresidual=coxph_qr_ph$nrsp,
                         fitted.value = fit_coxph_ph_data$linear.predictors);bl_lp_ph

bl_lp_nonph<-test.var.bartl(Zresidual=coxph_qr_nonph$nrsp,
                            fitted.value = fit_coxph_nonph_data$linear.predictors);bl_lp_nonph



anov_x1_ph<- test.nl.aov(Zresidual=coxph_qr_ph$nrsp,
                         fitted.value=ph_data$x1);anov_x1_ph
anov_x2_ph<- test.nl.aov(Zresidual=coxph_qr_ph$nrsp,
                         fitted.value=ph_data$x2);anov_x2_ph

bl_x1_ph<- test.var.bartl(Zresidual=coxph_qr_ph$nrsp,
                          fitted.value=ph_data$x1);bl_x1_ph
bl_x2_ph<- test.var.bartl(Zresidual=coxph_qr_ph$nrsp,
                          fitted.value=ph_data$x2);bl_x2_ph


anov_x1_nonph<- test.nl.aov(Zresidual=coxph_qr_nonph$nrsp,
                            fitted.value=nonph_data$x1);anov_x1_nonph
anov_x2_nonph<- test.nl.aov(Zresidual=coxph_qr_nonph$nrsp,
                            fitted.value=nonph_data$x2);anov_x2_nonph


bl_x1_nonph<- test.var.bartl(Zresidual=coxph_qr_nonph$nrsp,
                             fitted.value=nonph_data$x1);bl_x1_nonph
bl_x2_nonph<- test.var.bartl(Zresidual=coxph_qr_nonph$nrsp,
                             fitted.value=nonph_data$x2);bl_x2_nonph

pdf( "plot/qq_boxplot.pdf",width=10,height=7)
par(mfrow = c(2,3),mar=c(4,4,5,4))
ymax=max(range(coxph_qr_ph$nrsp,coxph_qr_nonph$nrsp),3)
ymin=min(range(coxph_qr_ph$nrsp,coxph_qr_nonph$nrsp),-3)

qqnorm(coxph_qr_ph$nrsp,xlab="Theoretical Quantiles",
       ylab="Sample Quantiles",ylim = c(ymin,ymax),
       main="Nonlinear Data, Z-Residual QQ Plot");qqline(coxph_qr_ph$nrsp)
text(-1, 3.3, paste0("Z-SW p-value = ", sprintf("%3.2f",sw_coxph_ph)))
legend(x = "top", inset = c(0,-0.35), text.font=1,
       legend = c("(a)"),cex=2,xpd = TRUE,bty="n")

x2_t <- cut(ph_data$x2, 10)
plot(x2_t,coxph_qr_ph$nrsp,xlab=expression(paste(X[2])),ylab="Z-Residual",
     main="Nonlinear Data, Z-Residual Boxplot",ylim=c(ymin,4.5)) 
text(3.8, 4.4, expression("Z-AOV-X"[2]*" p-value < 0.01"))
text(3.8, 3.8, expression("Z-BL-X"[2]*" p-value = "))
text(7, 3.8,paste0(sprintf("%3.2f",bl_x2_ph)))
legend(x = "top", inset = c(0,-0.35), text.font=1,
       legend = c("(b)"),cex=2,xpd = TRUE,bty="n")

plot(coxzph_test_ph, var = 2,main="Nonlinear Data, Schoenfeld Residual Plot",xlab="Event Time",ylim=c(-6,18))
text(6, 18, expression("Score Tests X"[2]*" p-value < 0.01"))
text(7.5, 16, paste0("Score Tests GLOBAL p-value < 0.01 "))
legend(x = "top", inset = c(0,-0.35), text.font=1,
       legend = c("(c)"),cex=2,xpd = TRUE,bty="n")

qqnorm(coxph_qr_nonph$nrsp,xlab="Theoretical Quantiles", 
       ylab="Sample Quantiles",ylim = c(ymin,ymax),
       main="Time-varying Data, Z-Residual QQ Plot");qqline(coxph_qr_nonph$nrsp)
text(-1, 3.3, paste0("Z-SW p-value < 0.01"))
legend(x = "top", inset = c(0,-0.35), text.font=1,
       legend = c("(d)"),cex=2,xpd = TRUE,bty="n")

x2_w <- cut(nonph_data$x2, 10)
plot(x2_w,coxph_qr_nonph$nrsp,xlab=expression(paste(X[2])),ylab="Z-Residual",
     main="Time-varying Data, Z-Residual Boxplot",ylim=c(ymin,4.5)) 
text(3.8, 4.4, expression("Z-AOV-X"[2]*" p-value < 0.01"))
text(3.8, 3.8, expression("Z-BL-X"[2]*" p-value < 0.01"))
legend(x = "top", inset = c(0,-0.35), text.font=1,
       legend = c("(e)"),cex=2,xpd = TRUE,bty="n")

plot(coxzph_test_nonph, var = 2,main="Time-varying Data, Schoenfeld Residual Plot",xlab="Event Time"
     ,ylim=c(-3,7))
text(7.5, 6.8, expression("Score Tests X"[2]*" p-value < 0.01 "))
text(9, 6, paste0("Score Tests GLOBAL p-value < 0.01 "))
legend(x = "top", inset = c(0,-0.35), text.font=1,
       legend = c("(f)"),cex=2,xpd = TRUE,bty="n")

dev.off()