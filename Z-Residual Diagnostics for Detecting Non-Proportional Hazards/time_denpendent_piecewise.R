setwd("~/cox_time_varying_nonph")
source("~/qresidual_coxph.R")
source("~/nonph_data.R")
library("survival")
library("EnvStats")
library("nortest")
library("backports")
library("stringr")
library("timereg")
library("MASS")
library("dvmisc")

#########piecewise###############################
####wrong model#################################################################
nonph_data<-simul_nonph1(n_clusters=10,n_individuals=50,fv=0.5,
                    lambda=0.007, alpha=3, 
                    mean.censor=7.5,beta1=0.3,beta_a=1.35,
                    beta_b=-1.35,t0=2,corr=0)

table(nonph_data$d)
fit_coxph_nonph_data <- tryCatch(
  coxph(Surv(t, d) ~ x1+x2+frailty(grpid,distribution = "gamma"),
        data = nonph_data),
  error = function(e) NA,
  warning = function(w) NA
)
coxzph_test_w<-cox.zph(fit_coxph_nonph_data, transform="identity");coxzph_test_w

coxph_qr_w<-zresidual.coxph (fit_coxph = fit_coxph_nonph_data,
                             traindata = nonph_data,
                             newdata = nonph_data)
sw_coxph_w<-shapiro.test(coxph_qr_w)$p.value;sw_coxph_w

## the coefficient of x1 is ph
anov_x1_w<- test.nl.aov1(qresidual=coxph_qr_w,
                         fitted.values=nonph_data$x1,
                         k.anova=10);anov_x1_w
bl_x1_w<-test.var.bartl(qresidual=coxph_qr_w,
                        fitted.values=nonph_data$x1,
                        k=10);bl_x1_w

## the coefficient of x2 is non-ph
anov_x2_w<- test.nl.aov1(qresidual=coxph_qr_w,
                      fitted.values=nonph_data$x2,
                      k.anova=10);anov_x2_w
bl_x2_w<-test.var.bartl(qresidual=coxph_qr_w,
                        fitted.values=nonph_data$x2,
                        k=10);bl_x2_w

anov_lp_w<- test.nl.aov1(qresidual=coxph_qr_w,
                         fitted.values=fit_coxph_nonph_data$linear.predictors,
                         k.anova=10);anov_lp_w

bl_lp_w<-test.var.bartl(qresidual=coxph_qr_w,
                        fitted.values = fit_coxph_nonph_data$linear.predictors,
                        k=10);bl_lp_w


################################################################################
##########################true model###########################################
ph_data<-simul_nonph1(n_clusters=10,n_individuals=50,fv=0.5,
                     lambda=0.007, alpha=3,
                     mean.censor=7.5,beta1=0.3,beta_a=1.35,
                     beta_b=1.35,t0=pi/2,corr=0)
table(ph_data$d)
fit_coxph_ph_data <- tryCatch(
  coxph(Surv(t, d) ~ x1+x2 +frailty(grpid,distribution = "gamma"),
        data = ph_data),
  error = function(e) NA,
  warning = function(w) NA
)
coxzph_test_t<-cox.zph(fit_coxph_ph_data, transform="identity");coxzph_test_t

coxph_qr_t<-zresidual.coxph (fit_coxph = fit_coxph_ph_data,
                             traindata = ph_data,
                             newdata = ph_data)
sw_coxph_t<-shapiro.test(coxph_qr_t)$p.value;sw_coxph_t

## the coefficient of x1 is ph
anov_x1_t<- test.nl.aov1(qresidual=coxph_qr_t,
                         fitted.values=ph_data$x1);anov_x1_t
bl_x1_t<-test.var.bartl(qresidual=coxph_qr_t,
                        fitted.values=ph_data$x1);bl_x1_t

## the coefficient of x2 is ph
anov_x2_t<- test.nl.aov1(qresidual=coxph_qr_t,
                         fitted.values=ph_data$x2);anov_x2_t
bl_x2_t<-test.var.bartl(qresidual=coxph_qr_t,
                        fitted.values=ph_data$x2);bl_x2_t

anov_lp_t<- test.nl.aov1(qresidual=coxph_qr_t,
                         fitted.values=fit_coxph_ph_data$linear.predictors)
anov_lp_t
bl_lp_t<-test.var.bartl(qresidual=coxph_qr_t,
                        fitted.values = fit_coxph_ph_data$linear.predictors)
bl_lp_t

