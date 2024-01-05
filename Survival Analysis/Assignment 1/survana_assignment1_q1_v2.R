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


#------------------------------------
# Load Data
#------------------------------------

data("larynx")
head(larynx, 20)

# rename 
dat <- larynx
# set stage as factor, NB for cox model
dat$stage <- as.factor(dat$stage)
#levels(dat$stage) <- c(0, 1, 2, 3)

# Create an event indicator:
rcensored <- dat$delta == 1  # T if patient died, False otherwise (if a patient didn't died, then they are censored)
rcensored

# Create the Survival object:
surv.rel <- Surv(dat$time, rcensored) # Produces survival times with a '+' added to those that were censored
surv.rel
survdata <- data.frame(dat$stage, dat$time, dat$age, dat$diagyr, dat$delta, surv.rel) 
head(survdata, 10)




#------------------------------------
# 1 (a)
#------------------------------------
# Test the hypothesis that the distribution of the times to death are the same in the 4 stages of the disease.

# First lets look at each stage by themselves:
# Stage 1
# surv1 <- Surv(dat$time[dat$stage == 1], rcensored[dat$stage == 1])
# surv1  
# fit01 <- survfit(surv1 ~ 1, conf.type = "log-log") # creating a survival curve
# km01  <- summary(fit01, censor=TRUE)
# km01
# plot(fit01)
# legend(0,0.3,legend=c("Lower CI","Survival","Upper CI"),lty=c(2,1,2))
# title(main="", xlab="Time (months)", ylab = "Proportion")
# grid()
# fit01

# Stage 2
# surv2 <- Surv(dat$time[dat$stage == 2], rcensored[dat$stage == 2])
# surv2  
# fit02 <- survfit(surv2 ~ 1, conf.type = "log-log") # creating a survival curve
# km02  <- summary(fit02, censor=TRUE)
# km02
# plot(fit02)
# legend(0,0.3,legend=c("Lower CI","Survival","Upper CI"),lty=c(2,1,2))
# title(main="", xlab="Time (months)", ylab = "Proportion")
# grid()
# fit02

# Stage 3
# surv3 <- Surv(dat$time[dat$stage == 3], rcensored[dat$stage == 3])
# surv3  
# fit03 <- survfit(surv3 ~ 1, conf.type = "log-log") # creating a survival curve
# km03  <- summary(fit03, censor=TRUE)
# km03
# plot(fit03)
# legend(0,0.25,legend=c("Lower CI","Survival","Upper CI"),lty=c(2,1,2))
# title(main="", xlab="Time (months)", ylab = "Proportion")
# grid()
# fit03

# Stage 4
# surv4 <- Surv(dat$time[dat$stage == 4], rcensored[dat$stage == 4])
# surv4  
# fit04 <- survfit(surv4 ~ 1, conf.type = "log-log") # creating a survival curve
# km04  <- summary(fit04, censor=TRUE)
# km04
# plot(fit04)
# legend(2.5,1,legend=c("Lower CI","Survival","Upper CI"),lty=c(2,1,2))
# title(main="", xlab="Time (months)", ylab = "Proportion")
# grid()
# fit04





# survival curve fit using all stages, where stage is a covariate
fit1 <- survfit(surv.rel ~ stage, data = dat)
fit1
km1 <- summary(fit1,censor=TRUE)
km1
plot(fit1, col=c("red","blue", "green", "black"))  #, xaxp=c(0,160,16)
title(main="", xlab="Time (months)", ylab = "Proportion")
legend(8,1,legend=c(1,2,3,4), lty=c(1,1,1,1), col=c("red","blue", "green", "black"))
grid()
ggsurvplot(fit1,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=TRUE, # show a risk table below the plot
           legend.labs=c(1,2,3,4), # change group labels
           legend.title="Stage",  # add legend title
           #palette=c("dodgerblue4", "orchid2"), # change colors of the groups
           title="", # add title to plot
           risk.table.height=.2)
ggsurvplot(fit1,
           conf.int=F, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=TRUE, # show a risk table below the plot
           legend.labs=c(1,2,3,4), # change group labels
           legend.title="Stage",  # add legend title
           #palette=c("dodgerblue4", "orchid2"), # change colors of the groups
           title="", # add title to plot
           risk.table.height=.2)




# Use the survdiff function to compare the 2 distributions using the logrank test:
survdiff(surv.rel ~ dat$stage)
# p=5e-05, therefore reject H0
# evidence that the death rates in the four stages are different



#------------------------------------
# 1 (b)
#------------------------------------
# Fit a Cox proportional hazards model including the covariates age and stage of disease. 
# Interpret the model results.
# Check the overall fit of the model and check all model assumptions.

cox1 <- coxph(surv.rel ~ dat$stage + dat$age)
summary(cox1)
# The hazard of death for patients in stage 2 is 1.150 more then that of patients in stage 1
# The hazard of death for patients in stage 3 is 1.901 more then that of patients in stage 1
# The hazard of death for patients in stage 4 is 5.507 more then that of patients in stage 1
# positive beta coefficient indicates an increased hazard with increasing age
# The hazard of death increases by 1.9% for every one year increase in age 
exp(0.01903*10)  # 1.209612
# The hazard of death increases by almost 21% for every 10 year increase in age
# CI for age for 10yrs:
exp(0.01903*10-1.96*10*0.01426) 
exp(0.01903*10+1.96*10*0.01426) 


# CHECK FIT: 


# Cox-Snell Residuals
# Assess overall model fit


# Score Residuals  (dfbeta)
# Look for influential observation 
ggcoxdiagnostics(cox1, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
dfbeta <- residuals(cox1,type="dfbeta")
for (j in 1:4) {
  plot(dfbeta[,j], ylab = names(coef(cox1))[j])
  abline(h=0,lty=2)
}
# "comparing the magnitudes of the largest dfbeta values to the regression coefficients"
# source: https://shariq-mohammed.github.io/files/cbsa2019/1-intro-to-survival.html#7_paramteric_surival_models
# Dont see any very influential values for age, maybe some for stages


# Martingale Residuals
# Check form of covariates
df <- dat
df$res <- residuals(cox1, type = "martingale")
ggplot(data = df, mapping = aes(x = age, y = res)) +
  geom_point() +
  geom_smooth() +
  labs(xlab = "Age", ylab = "Martingale Residuals") +
  theme_bw() + theme(legend.key = element_blank())
# roughly linear but not horizontal, indicating age should be included as is (not transformed)
ggplot(data = df, mapping = aes(x = age^2, y = res)) +
  geom_point() +
  geom_smooth() +
  labs(xlab = "Age Squared", ylab = "Martingale Residuals") +
  theme_bw() + theme(legend.key = element_blank())
ggplot(data = df, mapping = aes(x = sqrt(age), y = res)) +
  geom_point() +
  geom_smooth() +
  labs(xlab = "Age Square Rooted", ylab = "Martingale Residuals") +
  theme_bw() + theme(legend.key = element_blank())





# Deviance Residuals 
ggcoxdiagnostics(cox1, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())
# "The deviance residual is a normalized transform of the martingale residual."
# "These residuals should be roughly symmetrically distributed about zero with a standard deviation of 1."
# "Positive values correspond to individuals that “died too soon” compared to expected survival times."
# "Negative values correspond to individual that “lived too long”.
# "Very large or small values are outliers, which are poorly predicted by the model."
# source: https://shariq-mohammed.github.io/files/cbsa2019/1-intro-to-survival.html#7_paramteric_surival_models
# Seem to be symmetrically distributed around zero, suggesting adequate fit
# No obvious outliers



# Test the Proportional Hazards Assumption of a Cox Regression
cox.zph(cox1)  
# P-values are not statistically significant for either covariates or for the overall model
# Thus the assumption of proportional hazards seems reasonable
cz <- cox.zph(cox1)
plot(cz)
# "produces, for each covariate, graphs of the scaled Schoenfeld residuals against the transformed time."
# "Systematic departures from a horizontal line are indicative of non-proportional hazards, since proportional hazards 
# assumes that estimates β1,β2,β3 do not vary much over time."
# source: https://shariq-mohammed.github.io/files/cbsa2019/1-intro-to-survival.html#7_paramteric_surival_models
# Age seems to show some variation with time, suggesting that the assumption of proportional hazard my in fact not hold




#------------------------------------
# 1 (c)
#------------------------------------
# Using the fitted Cox Model in (b) to plot the survival functions for an individual who is aged 40 and is in stage 1
# and for an individual who is aged 50 and is in stage 4


# (i) aged 40 and is in stage 1
# We need to generate new variables to redefine the baseline function
dat1 <- dat
# stage 1 is already baseline so no need to adjust
dat1$age <- dat1$age - (mean(dat$age)-50)    # re-centering around 40
#cox2 <- coxph(surv.rel ~ stage + age, data = dat1)
# Plot the baseline survival function
#basehaz(cox2, centered=F)
ggsurvplot(survfit(cox1, data = dat1), palette = "#2E9FDF",
           ggtheme = theme_minimal())


# or
dat_new <- data.frame(stage=1,age=50)
survfit(cox1, newdata = dat_new)
plot(survfit(cox1, newdata = dat_new))


# or
cox2 <- coxph(surv.rel ~ stage + age, data = dat1)
# Plot the baseline survival function
b1 <- basehaz(cox2, centered=F)
b1
plot(b1, type='l')





# (i) aged 50 and is in stage 4
# We need to generate new variables to redefine the baseline function
dat2 <- dat
levels(dat2$stage) <-  c(4,3,2,1)
dat2$age <- dat1$age - (mean(dat$age)-50) 
cox3 <- coxph(surv.rel ~ stage + age, data = dat2)
# Plot the baseline survival function
#basehaz(cox3, centered=F)
ggsurvplot(survfit(cox1, data = dat2), palette = "#2E9FDF",
           ggtheme = theme_minimal())



# or
# Plot the baseline survival function
b2 <- basehaz(cox3, centered=F)
b2
plot(b2, type='l')
ggsurvplot(b2)

pred <- predict(cox1,data.frame(stage=rep(1,nrow(dat)), time=dat$time, age=rep(50,nrow(dat))))
plot(pred)

#------------------------------------
# 1 (d)
#------------------------------------
# Based on your survival functions, give the survival function estimates at 3 months for an
# individaul who is aged 50 and in stage 1 and for an individaul who is aged 50 and in stage 4. Give
# 95% confidence intervals for your estimates.

df1 <- read_dta("age50stage1.dta")
as.data.frame(df1)



df2 <- read_dta("age50stage4.dta")
as.data.frame(df2)


