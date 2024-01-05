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
library(eha)
library(ggfortify)


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

# Create the Survival object:
surv.rel <- Surv(dat$time, dat$delta)
surv.rel
survdata <- data.frame(stage=dat$stage, time=dat$time, age=dat$age, 
                       diagyr=dat$diagyr, delta=dat$delta, surv.rel) 
head(survdata, 10)

#write.csv(survdata, "larynx.csv", row.names = FALSE)  

#------------------------------------
# 2 
#------------------------------------
# Fit a Weibull PH model with stage as the only covariate.

fit1 <- flexsurvreg(surv.rel ~ stage, data = survdata, dist = "weibullph")
fit1

#------------------------------------
# 2 (a)
#------------------------------------
# Refit the model but now allow the shape parameter of the Weibull distributiuon to depend on stage

fit2 <- flexsurvreg(surv.rel~stage,
                    anc = list(shape = ~ stage),
                    data=survdata, dist = "weibullph")   
fit2




#------------------------------------
# 2 (a)(i)
#------------------------------------
# Write down the functional form of the fitted model





#------------------------------------
# 2 (a)(ii)
#------------------------------------
# Even though λ is parameterized similarly to that in a Weibull PH model, this model is not
# a PH model because the shape parameter p varies across treatment groups. 
# Show the PH assumption is violated in this model by estimating the hazard ratios for Stage 4 vs 
# Stage 1 after 1 month and after 3 months of follow-up.


# Note on weibull parameterisation used:
# The Weibull distribution with shape parameter a and scale parameter b has density given by
# f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a)
# for x > 0. The cumulative distribution function is F(x) = 1 - exp(- (x/b)^a) on x > 0, the mean is E(X) = b Γ(1 + 1/a), and the Var(X) = b^2 * (Γ(1 + 2/a) - (Γ(1 + 1/a))^2).

theta <- coef(fit2)
theta

# 1 months (see calculations in diary):
exp(-0.47929)*(1/0.02492)^(exp(1.68850-0.47929)-exp(1.68850))*exp(0.14444)
# therefor the hazard of death at 1 month for patients in stage 4 is 6.604521  that of patients 
# in stage 1 ?


# 3 months (see calculations in diary):
exp(theta[5]+theta[8])*3^(theta[8])




#------------------------------------
# 3 
#------------------------------------
# Fit a Weibull AFT model including the variables age and stage as covariates.
fit3 <- flexsurvreg(surv.rel~stage+age, data=survdata, dist = "weibull")
fit3

#------------------------------------
# 3 (a)
#------------------------------------
# From your fittted Weibull AFT model obtain the regression coefficients of the corresponding
# Webull PH model
beta <- -coef(fit3)
beta


#------------------------------------
# 3 (b)
#------------------------------------
# Give the value of the acceleration factor for an individual of age 50 in the fourth stage of the
# disease. Give a 95% confidence interval for your estimate.
# acceleration factor (af) = exp(beta*X)
af <- exp(beta[6]*50 + beta[5])
af

# 95 % confidence interval:
lage <- -(coef(fit3)[6]*50 + 0.01278*50*1.96)
uage <- -(coef(fit3)[6]*50 - 0.01278*50*1.96)

l4 <- -(coef(fit3)[5] + 0.36323*1.96)
u4 <- -(coef(fit3)[5] - 0.36323*1.96)

# lower bound:
exp(lage+l4)

# upper bound:
exp(uage+u4)




#------------------------------------
# 4
#------------------------------------

# Fit the exponential, lognormal and gamma AFT models including the variables age and stage as covariates. 
# Plot the hazard function by stage for each of these models. 
# Check the overall fits of the models. 
# Which model (including the Weibull AFT) is the most suitable for these data? 
# Compare the decceleration factors between these models


fit4a <- survreg(surv.rel~stage+age, data=survdata, dist = "exponential")
fit4b <- flexsurvreg(surv.rel~stage+age, data=survdata, dist = "lognormal")
fit4c <- flexsurvreg(surv.rel~stage+age, data=survdata, dist = "gamma")
fit4c
BIC(fit4c)


# NOT WORKING:
hazard_preds <- data.frame(age = c(65 , 65 , 65 , 65), stage = c("1","2","3","4")) # use the model to predict the hazards
haz_preds_exp <- summary(fit4a , newdata = hazard_preds ,
                         type = "hazard", tidy = T) # plot the predicted hazards
exp_haz_plot <- ggplot(haz_preds_exp, aes(x = time , y = est , col = stage))+ 
  geom_line(size = 1)+ 
  labs(x = "Months", y = "Exponential AFT n predicted  hazard", tag = "B")+ 
  scale_colour_manual(labels = c("Stage 1", "Stage 2", "Stage 3", "Stage 4"),
                      values = c("springgreen 3", "steelblue 3", "darkorange 2", "red2"),
                      name = "")+ 
  theme(legend.position = "top")
exp_haz_plot





hazard_preds <- data.frame(age = c(65, 65, 65, 65), stage = c("1","2","3","4")) # use the model to predict the hazards
haz_preds_exp <- summary(fit4b , newdata = hazard_preds ,
                         type = "hazard", tidy = T) # plot the predicted hazards
exp_haz_plot <- ggplot(haz_preds_exp, aes(x = time , y = est , col = stage))+ 
  geom_line(size = 1)+ 
  labs(x = "Months", y = "Exponential AFT n predicted  hazard", tag = "B")+ 
  scale_colour_manual(labels = c("Stage 1", "Stage 2", "Stage 3", "Stage 4"),
                      values = c("springgreen 3", "steelblue 3", "darkorange 2", "red2"),
                      name = "")+ 
  theme(legend.position = "top")
exp_haz_plot




hazard_preds <- data.frame(age = c(65, 65, 65, 65), stage = c("1","2","3","4")) # use the model to predict the hazards
haz_preds_exp <- summary(fit4c , newdata = hazard_preds ,
                         type = "hazard", tidy = T) # plot the predicted hazards
exp_haz_plot <- ggplot(haz_preds_exp, aes(x = time , y = est , col = stage))+ 
  geom_line(size = 1)+ 
  labs(x = "Months", y = "Gamma AFT Predicted  Hazard")+ 
  scale_colour_manual(labels = c("Stage 1", "Stage 2", "Stage 3", "Stage 4"),
                      values = c("springgreen 3", "steelblue 3", "darkorange 2", "red2"),
                      name = "")+ 
  theme(legend.position = "top", text = element_text(size = 30))
exp_haz_plot





cs <- coxsnell_flexsurvreg(fit4c)
qy <- qexp(ppoints(nrow(cs),0))
qqplot(cs$est, qy,  pch=19, xlab = "Cox-Snell Residuals", ylab = "")
abline(a=0,b=1, col="red", lwd=2)
grid()







#------------------------------------
# 5
#------------------------------------
# Fit the generalized Gamma AFT model ncluding the variables age and stage as covariates. 
# Check the overall fit of the model. 
# Use the fitted model to assess the suitability of the models fitted in questions 3 and 4.

fit5 <- flexsurvreg(surv.rel~stage+age, data=survdata, dist = "gengamma")
fit5a <- flexsurvreg(surv.rel~stage+age, data=survdata, dist = "gengamma", Q = 0, fixedpars = 3)  # lognormal
fit5b <- flexsurvreg(surv.rel~stage+age, data=survdata, dist = "gengamma", Q = 1, fixedpars = 3)  # weibull
flexsurvreg(surv.rel~stage+age, data=survdata, dist = "weibull", shape=1, fixedpars = 2)




# ok but now fit5c is ph ???
fit5c <- flexsurvreg(surv.rel~stage+age, data=survdata, dist = "gengamma", sigma = 1, fixedpars = 2) # gamma



# Not working:
fit5d <- flexsurvreg(surv.rel~stage+age, data=survdata, dist = "gengamma", Q = 1, sigma = 1, fixedpars = c(2,3)) #expo










