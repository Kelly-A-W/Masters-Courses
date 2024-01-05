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
# Question 1 (a)
#------------------------------------
# give summary statistics by treatment

# summary for continuous variables for treatment CBZ:
summary(dat_w[which(dat_w$treat=="CBZ"),c(2,4,10)])
# summary for factor variables for treatment CBZ:
round(table(dat_w[which(dat_w$treat=="CBZ"),5])/length(which(dat_w$treat=="CBZ")), 3) # with_status
round(table(dat_w[which(dat_w$treat=="CBZ"),6])/length(which(dat_w$treat=="CBZ")), 3) # with_status2   
round(table(dat_w[which(dat_w$treat=="CBZ"),7])/length(which(dat_w$treat=="CBZ")), 3) # with_status_uae
round(table(dat_w[which(dat_w$treat=="CBZ"),8])/length(which(dat_w$treat=="CBZ")), 3) # with_status_isc
round(table(dat_w[which(dat_w$treat=="CBZ"),11])/length(which(dat_w$treat=="CBZ")), 3) # gender
round(table(dat_w[which(dat_w$treat=="CBZ"),12])/length(which(dat_w$treat=="CBZ")), 3) # learn_dis

# summary for continuous variables for treatment LTG:
summary(dat_w[which(dat_w$treat=="LTG"),c(2,4,10)])
# summary for factor variables for treatment LTG:
round(table(dat_w[which(dat_w$treat=="LTG"),5])/length(which(dat_w$treat=="LTG")), 3) # with_status
round(table(dat_w[which(dat_w$treat=="LTG"),6])/length(which(dat_w$treat=="LTG")), 3) # with_status2
round(table(dat_w[which(dat_w$treat=="LTG"),7])/length(which(dat_w$treat=="LTG")), 3) # with_status_uae
round(table(dat_w[which(dat_w$treat=="LTG"),8])/length(which(dat_w$treat=="LTG")), 3) # with_status_isc
round(table(dat_w[which(dat_w$treat=="LTG"),11])/length(which(dat_w$treat=="LTG")), 3) # gender
round(table(dat_w[which(dat_w$treat=="LTG"),12])/length(which(dat_w$treat=="LTG")), 3) # learn_dis



#------------------------------------
# Question 1 (b)
#------------------------------------
# Fit cause-specific hazard models separate Cox models for treatment failure due to ISC and UAE
# Adjust for baseline dose, treatment, age, gender and learning disability
# Add interaction between dose and treatment

# ISC:
cox_isc <- coxph(Surv(with_time, with_status2 == 1) ~ dose*treat + age + gender + learn_dis, data=dat_w)
summary(cox_isc)
cox_isc$loglik
AIC(cox_isc)
BIC(cox_isc)
# Cox-Snell Residuals:
# Equal to the event indicator - the martingale residuals
# source:https://people.math.aau.dk/~rw/Undervisning/DurationAnalysis/Slides/lektion4.pdf
# and: https://longjp.github.io/survival/lectures/04coxph.html
csr_isc <- dat_w$with_status_isc - residuals(cox_isc, "martingale")
# Use NA method to estimate the cumulative hazard function for residuals
fit <- survfit(Surv(csr_isc,dat_w$with_status_isc)~1)
nacumhaz <- cumsum(fit$n.event/fit$n.risk)
# plot the results
par(mar = c(5.1, 5, 4.1, 2.1))
plot(fit$time,nacumhaz,type='s',col='blue', xlab = "Time", ylab = "Cumulative Hazard",  cex.axis = 1.5, 
     cex.lab = 2.5, lwd = 2)
abline(0,1,col='red',lty=2, lwd = 2)
grid()


# UAE:
cox_uae <- coxph(Surv(with_time, with_status2 == 2) ~ dose*treat + age + gender + learn_dis, data=dat_w)
summary(cox_uae)
print(summary(cox_uae),digits=3) 
cox_uae$loglik
AIC(cox_uae)
BIC(cox_uae)
# Cox-Snell Residuals:
csr_uae <- dat_w$with_status_uae - residuals(cox_uae, "martingale")
# Use NA method to estimate the cumulative hazard function for residuals
fit <- survfit(Surv(csr_uae,dat_w$with_status_uae)~1)
nacumhaz <- cumsum(fit$n.event/fit$n.risk)
# plot the results
par(mar = c(5.1, 5, 4.1, 2.1))
plot(fit$time,nacumhaz,type='s',col='blue', xlab = "Time", ylab = "Cumulative Hazard",  cex.axis = 1.5, 
     cex.lab = 2.5, lwd = 2)
abline(0,1,col='red',lty=2, lwd = 2)
grid()

#------------------------------------
# Question 1 (c)
#------------------------------------
# Fit a sub-distribution hazard model (Fine and Gray model) for ISC with UAE as a competing event
# Adjust for baseline dose, treatment, age, gender and learning disability
# Compare to model (b)

# Convert factors to numeric
dat_fg <- dat_w
dat_fg$treat <- gsub("CBZ", 0, dat_fg$treat)  # baseline (consistent with cox model!)
dat_fg$treat <- gsub("LTG", 1, dat_fg$treat)
dat_fg$treat <- as.numeric(dat_fg$treat)
dat_fg$gender <- gsub("F", 0, dat_fg$gender)  # baseline
dat_fg$gender <- gsub("M", 1, dat_fg$gender)
dat_fg$gender <- as.numeric(dat_fg$gender)
dat_fg$learn_dis <- gsub("No", 0, dat_fg$learn_dis)  # baseline
dat_fg$learn_dis <- gsub("Yes", 1, dat_fg$learn_dis)
dat_fg$learn_dis <- as.numeric(dat_fg$learn_dis)
str(dat_fg)

# Create matrix of covariates
x <- matrix(cbind(dat_fg$dose, dat_fg$treat, dat_fg$age, dat_fg$gender, dat_fg$learn_dis, dat_fg$dose*dat_fg$treat), 
            nrow = nrow(dat_fg), ncol = 6)

# Fit Fine and Gray Model
mod_fg <- crr(dat_fg$with_time, dat_fg$with_status2, x, failcode = 1) 
# failcode denotes which event in dat_fg$with_status2 is the event of interest
# censored observations are assumed to be = 0 by default
summary(mod_fg)
mod_fg$loglik
# "The residuals returned are analogous to the Schoenfeld residuals in ordinary survival models. 
# Plotting the jth column of res against the vector of unique failure times checks for lack of fit over time
# in the corresponding covariate (column of cov1)"

plot(mod_fg$uftime, mod_fg$res[,1])

# Maybe see https://cran.r-project.org/web/packages/mstate/vignettes/Tutorial.pdf page 43 onwards





