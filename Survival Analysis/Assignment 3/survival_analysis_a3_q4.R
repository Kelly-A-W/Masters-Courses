# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#------------------------------------
# Load Libraries
#------------------------------------


library(survival)
library(KMsurv)
library(ggplot2)
library(survminer)
library(cmprsk)
library(nlme)
library(JM)



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
# Question 4:
#------------------------------------
# Joint modelling of longitudinal data and time to treatment failure in a COMPETING RISK setting
# UAE is the competing event 


dat_idCR <- crLong(dat_id, statusVar = "with_status2", censLevel = 0, nameStrata = "CR")
head(dat_idCR)

lm  <- lme(dose ~ time, data = dat, random = ~time|id)
cox <- coxph(Surv(with_time, status2) ~ strata(CR)*treat + age + gender + learn_dis, data=dat_idCR, x=T)
jm  <- jointModel(lm, cox, timeVar = "time", CompRisk = T, 
                  method="spline-PH-GH", 
                  interFact = list(value = ~CR, data=dat_idCR))
summary(jm)

# residuals() is not currently implemented for competing risks joint models.




