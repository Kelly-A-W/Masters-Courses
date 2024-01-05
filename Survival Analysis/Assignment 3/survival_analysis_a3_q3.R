# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#------------------------------------
# Load Libraries
#------------------------------------

library(flexsurv)
library(survival)
library(KMsurv)
library(ggplot2)
library(survminer)
library(haven)
library(cmprsk)
library(nlme)
library(JM)
library(JMbayes)



#------------------------------------
# Load Data
#------------------------------------

dat <- read.csv("epileptic.csv")
str(dat)

# Rescale Time Days to Years 
# divide by 365 (ignoring leap years)
dat$time <- dat$time/365
dat$with_time <- dat$with_time/365

# Set Factor Variables
dat$treat <- as.factor(dat$treat)
dat$gender <- as.factor(dat$gender)
dat$learn_dis <- as.factor(dat$learn_dis)


# Need only one observation per patient
dat_id <- dat[!duplicated(dat$id),-(13:18)]




#------------------------------------
# Question 3 (a)
#------------------------------------

# ISC
cox_isc <- coxph(Surv(with_time, with_status2 == 1) ~ treat + age + gender + learn_dis, data=dat_id, x=T)
lm      <- lme(dose ~ time, data = dat, random = ~time|id, na.action = na.omit)
jm_isc  <- jointModel(lm, cox_isc, timeVar = "time", method = "spline-PH-GH")
summary(jm_isc)
# Cox Snell Residuals
# "If the assumed model fits the data well, we expect the Cox-Snell residuals to have a unit exponential distribution;
# however, when Ti is censored, the residuals will be censored as well. To take censoring into account
# in checking the fit of the model, we can compare graphically the Kaplan-Meier estimate of
# the survival function of residuals with the survival function of the unit exponential distribution.
# JM: An R Package for the Joint Modelling of Longitudinal and Time-to-Event Data, 2010
resCST_isc <- residuals(jm_isc, process = "Event", type = "CoxSnell")
sfit_isc <- survfit(Surv(resCST_isc, with_status_isc) ~ 1, data = dat_id)
plot(sfit_isc, mark.time = FALSE, conf.int = TRUE, lty = c(1,2,2),
     xlab="Cox-Snell Residuals",ylab="Survival Probability")
curve(exp(-x), from = 0, to = max(dat$time), add = TRUE, col="red",lwd=2)
grid()



# UAE
cox_uae <- coxph(Surv(with_time, with_status2 == 2) ~ treat + age + gender + learn_dis, data=dat_id, x=T)
jm_uae <- jointModel(lm, cox_uae, timeVar = "time", method = "spline-PH-GH")
summary(jm_uae)
# Cox-Snell Residuals
resCST_uae <- residuals(jm_uae, process = "Event", type = "CoxSnell")
sfit_uae <- survfit(Surv(resCST_uae, with_status_uae) ~ 1, data = dat_id)
plot(sfit_uae, mark.time = FALSE, conf.int = TRUE, lty = c(1,2,2),
     xlab="Cox-Snell Residuals",ylab="Survival Probability")
curve(exp(-x), from = 0, to = max(dat$time), add = TRUE, col="red",lwd=2)
grid()



#------------------------------------
# Question 3 (b)
#------------------------------------

######## (i) #########
# Lagged Parameterisation 
# ISC
jm_isc_i <- jointModel(lm, cox_isc, timeVar = "time", lag = 1, method = "spline-PH-GH")  
summary(jm_isc_i)
# Cox Snell 
resCST_isc_i <- residuals(jm_isc_i, process = "Event", type = "CoxSnell")
sfit_isc_i   <- survfit(Surv(resCST_isc_i, with_status_isc) ~ 1, data = dat_id)
plot(sfit_isc_i, mark.time = FALSE, conf.int = TRUE, lty = c(1,2,2),
     xlab="Cox-Snell Residuals",ylab="Survival Probability")
curve(exp(-x), from = 0, to = max(dat$time), add = TRUE, col="red",lwd=2)
grid()
# UAE
jm_uae_i <- jointModel(lm, cox_uae, timeVar = "time", lag = 1, method = "spline-PH-GH")
summary(jm_uae_i)
# Cox-Snell Residuals
resCST_uae_i <- residuals(jm_uae_i, process = "Event", type = "CoxSnell")
sfit_uae_i <- survfit(Surv(resCST_uae_i, with_status_uae) ~ 1, data = dat_id)
plot(sfit_uae_i, mark.time = FALSE, conf.int = TRUE, lty = c(1,2,2),
     xlab="Cox-Snell Residuals",ylab="Survival Probability")
curve(exp(-x), from = 0, to = max(dat$time), add = TRUE, col="red",lwd=2)
grid()


######## (ii) #########
# Time Dependent Slope Parameterisation
lm_bii <- lme(dose ~ time + I(time^2), data = dat, random = list(id = pdDiag(form = ~ time + I(time)^2)))
dForm <- list(fixed = ~ I(2*time), 
              random = ~ I(2*time), 
              indFixed = c(2,3), 
              indRandom = c(2,3))
# ISC
jm_isc_ii <- jointModel(lm_bii, cox_isc, timeVar = "time", parameterization = "both", derivForm = dForm, method = "spline-PH-GH")  
summary(jm_isc_ii)
# Cox Snell 
resCST_isc_ii <- residuals(jm_isc_ii, process = "Event", type = "CoxSnell")
sfit_isc_ii   <- survfit(Surv(resCST_isc_ii, with_status_isc) ~ 1, data = dat_id)
plot(sfit_isc_ii, mark.time = FALSE, conf.int = TRUE, lty = c(1,2,2),
     xlab="Cox-Snell Residuals",ylab="Survival Probability")
curve(exp(-x), from = 0, to = max(dat$time), add = TRUE, col="red",lwd=2)
grid()
# UAE
jm_uae_ii <- jointModel(lm_bii, cox_uae, timeVar = "time", parameterization = "both", derivForm = dForm, method = "spline-PH-GH")  
summary(jm_uae_ii)
# Cox Snell 
resCST_uae_ii <- residuals(jm_uae_ii, process = "Event", type = "CoxSnell")
sfit_uae_ii   <- survfit(Surv(resCST_uae_ii, with_status_uae) ~ 1, data = dat_id)
plot(sfit_uae_ii, mark.time = FALSE, conf.int = TRUE, lty = c(1,2,2),
     xlab="Cox-Snell Residuals",ylab="Survival Probability")
curve(exp(-x), from = 0, to = max(dat$time), add = TRUE, col="red",lwd=2)
grid()


######## (iii) #########
# Random Slope Parameterisation
# ISC
jm_isc_iii <- jointModelBayes(lm, cox_isc, timeVar = "time", param = "shared-RE")
summary(jm_isc_iii)
# Cox Snell 
resCST_isc_iii <- dat$with_status_isc - residuals(jm_isc_iii, process="event", type = "Martingale")
sfit_isc_iii   <- survfit(Surv(resCST_isc_iii, with_status_isc) ~ 1, data = dat)
plot(sfit_isc_iii, mark.time = FALSE, conf.int = TRUE, lty = c(1,2,2),
     xlab="Cox-Snell Residuals",ylab="Survival Probability")
curve(exp(-x), from = 0, to = max(dat$time), add = TRUE, col="red",lwd=2)
grid()
# UAE
jm_uae_iii <- jointModelBayes(lm, cox_uae, timeVar = "time", param = "shared-RE")
summary(jm_uae_iii)
# Cox Snell 
resCST_uae_iii <- residuals(jm_uae_iii, process = "Event", type = "CoxSnell")
sfit_uae_iii   <- survfit(Surv(resCST_uae_ii, with_status_uae) ~ 1, data = dat_id)
plot(sfit_uae_iii, mark.time = FALSE, conf.int = TRUE, lty = c(1,2,2),
     xlab="Cox-Snell Residuals",ylab="Survival Probability")
curve(exp(-x), from = 0, to = max(dat$time), add = TRUE, col="red",lwd=2)
grid()



#------------------------------------
# Question 3 (c)
#------------------------------------

# Dichotomise Dose
dat$dose_bin <- dat$dose
dat$dose_bin[which(dat$dose <= 2)] <- 0
dat$dose_bin[which(dat$dose > 2)]  <- 1
#dat$dose_bin <- as.factor(dat$dose_bin)
# Check:
unique(dat$dose_bin)

# ISC
cox_isc <- coxph(Surv(with_time, with_status2 == 1) ~ treat + age + gender + learn_dis, data=dat_id, model=T)
lm_bin <- mvglmer(list(dose_bin ~ time + (time|id)), data = dat, families = list("binomial"))
set.seed(123)
jm_isc_bin <- mvJointModelBayes(lm_bin, cox_isc, timeVar = "time")
summary(jm_isc_bin)


# UAE
cox_uae <- coxph(Surv(with_time, with_status2 == 2) ~  treat + age + gender + learn_dis, data=dat_id, model=T)
set.seed(123)
jm_uae_bin <- mvJointModelBayes(lm_bin, cox_uae, timeVar = "time")
summary(jm_uae_bin)






