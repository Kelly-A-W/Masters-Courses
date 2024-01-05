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
library(JMbayes2)


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
dat_w <- dat[!duplicated(dat$id),-(13:18)]




#------------------------------------
# Question 2 (a)
#------------------------------------
# Fit separate Cox models for treatment failure due to ISC and for UAE
# Adjust for baseline dose, treatment, age, gender and learning disability
# Add interaction between dose and treatment

# ISC:
cox_isc <- coxph(Surv(with_time, with_status2 == 1) ~ dose*treat + age + gender + learn_dis, data=dat_w)
summary(cox_isc)
# exp(treatLTG) = 1.7414, thus the hazard of ISC is 1.7414  greater for patients receiving treatment LTG compared to those receiving treatment CBZ 
# (ignoring UAE)
# exp(age) = 0.9861, thus the hazard of ISC decreases slightly by less than 2% for each year increase in age
# (ignoring UAE)

# UAE:
cox_uae <- coxph(Surv(with_time, with_status2 == 2) ~ dose*treat + age + gender + learn_dis, data=dat_w)
summary(cox_uae)
# exp(treatLTG) = 0.5222, thus the hazard of UAE is 0.5222 greater for patients receiving treatment LTG compared
# to those receiving treatment CBZ  (ignoring ISC)



#------------------------------------
# Question 2 (b)
#------------------------------------
# Fit separated extended Cox models for treatment failure due to ISC and due to UAE
# Adjust for baseline dose, treatment, age, gender and learning disability
# Add interaction between dose and treatment
# Compare results to a

# Now we want data in log format
# so use dat

# ISC
# need to first remove redundant rows where dat$start1==dat$stop1 (seem to be one for each id)
# these cause problems for extended cox model
dat_isc_ec <- dat[-which(dat$start1==dat$stop1),]
ext_cox_isc <- coxph(Surv(start1, stop1, event1) ~ dose*treat + age + gender + learn_dis, data=dat_isc_ec)
summary(ext_cox_isc)
ext_cox_isc$loglik
AIC(ext_cox_isc)
BIC(ext_cox_isc)
# Cox-Snell Residuals:
csr_isc<- dat_isc_ec$event1 - residuals(ext_cox_isc, "martingale")
# Use NA method to estimate the cumulative hazard function for residuals
fit <- survfit(Surv(csr_isc,dat_isc_ec$event1)~1)
nacumhaz <- cumsum(fit$n.event/fit$n.risk)
# plot the results
par(mar = c(5.1, 5, 4.1, 2.1))
plot(fit$time,nacumhaz,type='s',col='blue', xlab = "Time", ylab = "Cumulative Hazard",  cex.axis = 1.5, 
     cex.lab = 2.5, lwd = 2)
abline(0,1,col='red',lty=2, lwd = 2)
grid()



# UAE
# need to first remove redundant rows where dat$start2==dat$stop2 (seem to be one for each id)
# these cause problems for extended cox model
dat_uae_ec <- dat[-which(dat$start2==dat$stop2),]
ext_cox_uae <- coxph(Surv(start2, stop2, event2) ~ dose*treat + age + gender + learn_dis, data=dat_uae_ec)
summary(ext_cox_uae)
ext_cox_uae$loglik
AIC(ext_cox_uae)
BIC(ext_cox_uae)
# Cox-Snell Residuals:
csr_uae<- dat_uae_ec$event2 - residuals(ext_cox_uae, "martingale")
# Use NA method to estimate the cumulative hazard function for residuals
fit <- survfit(Surv(csr_uae,dat_uae_ec$event2)~1)
nacumhaz <- cumsum(fit$n.event/fit$n.risk)
# plot the results
par(mar = c(5.1, 5, 4.1, 2.1))
plot(fit$time,nacumhaz,type='s',col='blue', xlab = "Time", ylab = "Cumulative Hazard",  cex.axis = 1.5, 
     cex.lab = 2.5, lwd = 2)
abline(0,1,col='red',lty=2, lwd = 2)
grid()




#------------------------------------
# Question 2 (c)
#------------------------------------
# Fit separate two-stage joint models for:
# (i) dose and treatment failure due to ISC
# (ii) dose and treatment failure due to UAE
# Adjust for baseline treatment, age, gender and learning disability
# Add interaction between dose and treatment
# Compare results to b

ggplot(data = dat, aes(x = time, y = dose, group = id))+geom_line()+
  theme_minimal(base_size = 22)

# ISC
lm_isc <- lme(dose ~ time, data = dat_isc_ec, random = ~time|id, na.action = na.omit)
fitval_isc <- fitted(lm_isc)
jm2stage_isc <- coxph(Surv(start1, stop1, event1) ~ treat*fitval_isc + age + gender + learn_dis , data = dat_isc_ec)
summary(jm2stage_isc)
# Cox-Snell Residuals:
csr_isc2<- dat_isc_ec$event1 - residuals(jm2stage_isc, "martingale")
# Use NA method to estimate the cumulative hazard function for residuals
fit <- survfit(Surv(csr_isc2,dat_isc_ec$event1)~1)
nacumhaz <- cumsum(fit$n.event/fit$n.risk)
# plot the results
par(mar = c(5.1, 5, 4.1, 2.1))
plot(fit$time,nacumhaz,type='s',col='blue', xlab = "Time", ylab = "Cumulative Hazard",  cex.axis = 1.5, 
     cex.lab = 2.5, lwd = 2)
abline(0,1,col='red',lty=2, lwd = 2)
grid()



# UAE
lm_uae <- lme(dose ~ time, data = dat_uae_ec, random = ~time|id, na.action = na.omit)
fitval_uae <- fitted(lm_uae)
jm2stage_uae<- coxph(Surv(start1, stop1, event2) ~ treat*fitval_uae + age + gender + learn_dis , data = dat_uae_ec)
summary(jm2stage_uae)
# Cox-Snell Residuals:
csr_uae2 <- dat_uae_ec$event2 - residuals(jm2stage_uae, "martingale")
# Use NA method to estimate the cumulative hazard function for residuals
fit <- survfit(Surv(csr_uae2,dat_uae_ec$event2)~1)
nacumhaz <- cumsum(fit$n.event/fit$n.risk)
# plot the results
par(mar = c(5.1, 5, 4.1, 2.1))
plot(fit$time,nacumhaz,type='s',col='blue', xlab = "Time", ylab = "Cumulative Hazard",  cex.axis = 1.5, 
     cex.lab = 2.5, lwd = 2)
abline(0,1,col='red',lty=2, lwd = 2)
grid()

# ISC
#lm_isc <- lme(dose ~ treat + age + gender + learn_dis + time, data = dat, random = ~time|id, na.action = na.omit)
#summary(lm_isc)
#jm_isc <- jm(cox_isc, lm_isc, time_var = "time")
#summary(jm_isc)


# UAE
#lm_uae <- lme(dose ~ treat + age + gender + learn_dis + time, data = dat, random = ~time|id, na.action = na.omit)
#summary(lm_uae)
#jm_uae <- jm(cox_uae, lm_uae, time_var = "time")
#summary(jm_uae)

