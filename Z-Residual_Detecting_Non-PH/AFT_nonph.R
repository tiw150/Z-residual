setwd("~/Z-Residual Diagnostics for Detecting Non-Proportional Hazards/")
source("~/AFTLN_data.R")
source("~/qresidual_coxph.R")
library("survival")
library("EnvStats")
library("nortest")
library("backports")
library("stringr")
library("timereg")
library("MASS")
library("dvmisc")
library("timereg")

###########AFT Log-normal Model#################################
aft_ln_data<-AFTLN_data(n_clusters=10,n_individuals=50,fv=0.5,
                         mu=0, sigma=1,mean.censor=3.8,beta1=1)
table(aft_ln_data$d)

fit_coxph_aftln_data <- tryCatch(
  coxph(Surv(t, d) ~ X+frailty(grpid,distribution = "gamma"),
        data = aft_ln_data),
  error = function(e) NA,
  warning = function(w) NA
)

coxzph_test_w<-cox.zph(fit_coxph_aftln_data,transform="identity");coxzph_test_w
coxph_qr_ln<-zresidual.coxph (fit_coxph = fit_coxph_aftln_data,
                             traindata = aft_ln_data,
                             newdata = aft_ln_data)

sw_coxph_ln<-shapiro.test(coxph_qr_ln)$p.value;sw_coxph_ln
anov_lp_ln<- test.nl.aov1(qresidual=coxph_qr_ln,
                         fitted.values=fit_coxph_aftln_data$linear.predictors);anov_lp_ln
anov_x1_ln<- test.nl.aov1(qresidual=coxph_qr_ln,
                         fitted.values=aft_ln_data$x1);anov_x1_ln

bl_lp_ln<-test.var.bartl(qresidual=coxph_qr_ln,
                        fitted.values = fit_coxph_aftln_data$linear.predictors);bl_lp_ln
bl_x1_ln<-test.var.bartl(qresidual=coxph_qr_ln,
                        fitted.values = aft_ln_data$x1);bl_x1_ln

fit_timecox_aftln_data<- tryCatch(
  timecox(Surv(t,d)~X+ cluster(frailty(grpid,distribution = "gamma")),
          residuals=1,data=aft_ln_data,max.time=7,n.sim=100),
  error = function(e) NA,
  warning = function(w) NA
)
cumresid_test<-cum.residuals(fit_timecox_aftln_data,data=aft_ln_data,cum.resid=1,n.sim=100)
cum_x1_ln<-cumresid_test[["pval.test"]]



###########AFT Weibull Model##############################################
aft_wb_data<-AFTWB_data(n_clusters=10,n_individuals=50,fv=0.5,
                        lambda=1.67 , alpha=0.48 ,mean.censor=0.2,beta1=1)
table(aft_wb_data$d)
fit_coxph_aftwb_data <- tryCatch(
  coxph(Surv(t, d) ~ X+frailty(grpid,distribution = "gamma"),
        data = aft_wb_data),
  error = function(e) NA,
  warning = function(w) NA
)

coxzph_test_wb<-cox.zph(fit_coxph_aftwb_data,transform="identity");coxzph_test_wb
coxph_qr_wb<-zresidual.coxph (fit_coxph = fit_coxph_aftwb_data,
                             traindata = aft_wb_data,
                             newdata = aft_wb_data)

sw_coxph_wb<-shapiro.test(coxph_qr_wb)$p.value;sw_coxph_wb
anov_lp_wb<- test.nl.aov1(qresidual=coxph_qr_wb,
                         fitted.values=fit_coxph_aftwb_data$linear.predictors);anov_lp_wb
anov_x1_wb<- test.nl.aov1(qresidual=coxph_qr_wb,
                         fitted.values=aft_wb_data$x1);anov_x1_wb

bl_lp_wb<-test.var.bartl(qresidual=coxph_qr_wb,
                        fitted.values = fit_coxph_aftwb_data$linear.predictors);bl_lp_wb
bl_x1_wb<-test.var.bartl(qresidual=coxph_qr_wb,
                        fitted.values = aft_wb_data$x1);bl_x1_wb


fit_timecox_aftwb_data<- tryCatch(
  timecox(Surv(t,d)~X+ cluster(frailty(grpid,distribution = "gamma")),
          residuals=1,data=aft_wb_data,max.time=7,n.sim=100),
  error = function(e) NA,
  warning = function(w) NA
)
cumresid_test_wb<-cum.residuals(fit_timecox_aftwb_data,data=aft_wb_data,cum.resid=1,n.sim=100)
cum_x1_wb<-cumresid_test[["pval.test"]]
