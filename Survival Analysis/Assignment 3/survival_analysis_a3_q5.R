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
library(lcmm)


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
# Question 5(a):
#------------------------------------

# ISC
pre_lcjm_isc <- Jointlcmm(fixed = dose ~ time, 
                          random = ~ time,  # random-effects in the linear mixed model.
                          subject = "id", 
                          ng = 1,  # optional number of latent classes considered
                          data = dat, 
                          survival = Surv(with_time, with_status_isc)~treat + age + gender + learn_dis, 
                          hazard = "6-quant-piecewise", 
                          hazardtype = "Specific")  #  class-specific baseline r
lcjm_isc <- Jointlcmm(fixed = dose ~ time, 
                      mixture = ~ time, # class-specific fixed effects in the linear mixed model 
                      random = ~ time,  # random-effects in the linear mixed model.
                      classmb = ~ treat, 
                      subject = "id", 
                      ng = 3,  # optional number of latent classes considered
                      data = dat, 
                      survival = Surv(with_time, with_status_isc)~mixture(treat) + age + gender + learn_dis, 
                      hazard = "6-quant-piecewise", 
                      hazardtype = "Specific",   #  class-specific baseline risk function
                      B = pre_lcjm_isc)  
summary(lcjm_isc)
# ng=2  => maximum log-likelihood: -3153.74, AIC: 6357.48, BIC: 6467.61 
# ng=3  => maximum log-likelihood: -3108.27, AIC: 6286.55, BIC: 6440.73  
# ng=4  => maximum log-likelihood: -3081.99, AIC: 6253.98, BIC: 6452.21  
# Choose ng=3
# Assess the goodness-of-fit of the model:
postprob(lcjm_isc)  # the posterior class-membership probability for subject i 
# The proportion of subjects classified in each latent class with a posterior probability
# above 0.7, 0.8 and 0.9. This indicates the proportion of subjects not ambiguously
# classified in each latent class.
# The posterior classification table as defined in Table 3 which computes the mean of the
# posterior probabilities of belonging to the latent class among the subjects classified a
# posteriori in each latent class. A perfect classification would provide ones in the diagonal
# and zeros elsewhere. In practice, high diagonal terms indicate a good discrimination of
# the population.
# Source: Estimation of Extended Mixed Models Using Latent Classes and Latent Processes
#plot(lcjm_isc, which = "fit", var.time = "age65", marg = FALSE,
#     break.times = 10, bty = "l", ylab = "normMMSE",
#     xlab = "Age in decades from 65 years")




# UAE
pre_lcjm_uae <- Jointlcmm(fixed = dose ~ time, 
                          random = ~ time,  # random-effects in the linear mixed model.
                          subject = "id", 
                          ng = 1,  # optional number of latent classes considered
                          data = dat, 
                          survival = Surv(with_time, with_status_uae)~treat + age + gender + learn_dis, 
                          hazard = "6-quant-piecewise", 
                          hazardtype = "Specific")  #  class-specific baseline r
lcjm_uae <- Jointlcmm(fixed = dose ~ time, 
                      mixture = ~ time, # class-specific fixed effects in the linear mixed model 
                      random = ~ time,  # random-effects in the linear mixed model.
                      classmb = ~ treat, 
                      subject = "id", 
                      ng = 2,  # optional number of latent classes considered
                      data = dat, 
                      survival = Surv(with_time, with_status_uae)~mixture(treat) + age + gender + learn_dis, 
                      hazard = "6-quant-piecewise", 
                      hazardtype = "Specific",   #  class-specific baseline risk function
                      B = pre_lcjm_uae) 
summary(lcjm_uae)
# ng=2  => maximum log-likelihood: -3101.88, AIC: 6253.77, BIC: 6363.9 
# ng=3  => maximum log-likelihood: -3087.96, AIC: 6245.92, BIC: 6400.11  
# ng=4  => maximum log-likelihood: -3045.71, AIC: 6181.42, BIC: 6379.65  
# Choose ng=2
postprob(lcjm_uae)





#------------------------------------
# Question 5(b):
#------------------------------------

# Dichotomise Dose
dat$dose_bin <- dat$dose
dat$dose_bin[which(dat$dose <= 2)] <- 0
dat$dose_bin[which(dat$dose > 2)] <- 1
#dat$dose_bin <- as.factor(dat$dose_bin)
# Check:
unique(dat$dose_bin)


# ISC
pre_lcjm_isc_d <- Jointlcmm(fixed = dose_bin ~ time, 
                            random = ~ time,  # random-effects in the linear mixed model.
                            subject = "id", 
                            ng = 1,  # optional number of latent classes considered
                            data = dat, 
                            survival = Surv(with_time, with_status_isc)~treat + age + gender + learn_dis, 
                            hazard = "6-quant-piecewise", 
                            hazardtype = "Specific")  #  class-specific baseline r
lcjm_isc_d <- Jointlcmm(fixed = dose_bin ~ time, 
                      mixture = ~ time, # class-specific fixed effects in the linear mixed model 
                      random = ~ time,  # random-effects in the linear mixed model.
                      classmb = ~ treat, 
                      subject = "id", 
                      ng = 3,  # optional number of latent classes considered
                      data = dat, 
                      survival = Surv(with_time, with_status_isc)~mixture(treat) + age + gender + learn_dis, 
                      hazard = "6-quant-piecewise", 
                      hazardtype = "Specific",   #  class-specific baseline risk function
                      B = pre_lcjm_isc_d)  
summary(lcjm_isc_d)
# ng=2  => maximum log-likelihood: -1099.92, AIC: 2249.84, BIC: 2359.97  
# ng=3  => maximum log-likelihood: -910.83, AIC: 1891.65, BIC: 2045.84  
# ng=4  =>  maximum log-likelihood: -810.23, AIC: 1710.47, BIC: 1908.7  
# ng=5  => maximum log-likelihood: -764.28, AIC: 1638.56, BIC: 1880.85 
# Choose ng= 5
postprob(lcjm_isc_d)


# UAE
pre_lcjm_uae_d <- Jointlcmm(fixed = dose_bin ~ time, 
                            random = ~ time,  # random-effects in the linear mixed model.
                            subject = "id", 
                            ng = 1,  # optional number of latent classes considered
                            data = dat, 
                            survival = Surv(with_time, with_status_uae)~treat + age + gender + learn_dis, 
                            hazard = "6-quant-piecewise", 
                            hazardtype = "Specific")  #  class-specific baseline r
lcjm_uae_d <- Jointlcmm(fixed = dose_bin ~ time, 
                        mixture = ~ time, # class-specific fixed effects in the linear mixed model 
                        random = ~ time,  # random-effects in the linear mixed model.
                        classmb = ~ treat, 
                        subject = "id", 
                        ng = 3,  # optional number of latent classes considered
                        data = dat, 
                        survival = Surv(with_time, with_status_uae)~mixture(treat) + age + gender + learn_dis, 
                        hazard = "6-quant-piecewise", 
                        hazardtype = "Specific",   #  class-specific baseline risk function
                        B = pre_lcjm_uae_d)  
summary(lcjm_uae_d)
# ng=2  => maximum log-likelihood: -1013.12, AIC: 2076.25, BIC: 2186.38 
# ng=3  => maximum log-likelihood: -831.4, AIC: 1732.81, BIC: 1886.99
# ng=4  =>  maximum log-likelihood: -826.59, AIC: 1743.18, BIC: 1941.42  
# Choose ng=3
postprob(lcjm_uae_d)





