# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(nlme)
library(lmerTest)


#-------------------
# Load Data
#-------------------
head(Oxboys)
str(Oxboys) # Note that it seems actualy better to use data in this format, some functions don't seem to work using data.frame

dat <- as.data.frame(Oxboys)
head(dat)
str(dat)


#-------------------
# Question 2(a)
#-------------------
# plot data and verify that a simple linear regression model gives a suitable representation of the boys' growth patterns
# Do there appear to be significant differences in the individual growth patterns?

# spaghetti plots
p1 <- ggplot(data = dat, aes(x = age, y = height, group = Subject, colour = Subject))
p2 <- ggplot(data = dat, aes(x = Occasion, y = height, group = Subject, colour = Subject))
p1 + geom_line()
p2 + geom_line()
# two plots are pretty much the same, as expected
# can see a definite linear increase in height for all subjects as age (and thus occasion) increases - thus there seems to be a time effect
# no obvious difference in slope between subjects, lines fairly parallel - so no evidence of age*height interaction
# levels do seem to be different, different intercepts - evidence of subject random effect

# could also plot using:
interaction.plot(dat$Occasion,dat$Subject, dat$height, xlab="occasion", ylab="height", legend=F)

#-------------------
# Question 2(b)
#-------------------
# fit a simple linear regression model to height versus age using the lm function, ignoring the subject effects
# obtain the boxplots of the residuals by Subject and explain observed pattern

fit1 <- lm(height~age, data = dat)  # use Oxboys instead !!!!!
summary(fit1)
boxplot(resid(fit1)~dat$Subject,horizontal=TRUE, xlab = "Residuals", ylab = "Subject")
# plotting the residuals from fit1 grouped according to which subject
# can see that grouping by subject is evident in the residuals. 
# from this plot, we can see that subjects that are relatively short (eg: subject 10) tend to have their height overestimated, thus have high neagtive residuals
# subjects that are relatively tall (eg subject 4), tend to have their height underestimated, and thus have high positive residuals
# we want to have a relatively 50/50 overestimation/underestimation for all subjects 
# thus we want our residuals to ideally be centered around zero


#-------------------
# Question 2(c)
#-------------------
# use lmList function to fit separate simple linear regression models for each subject
# compare boxplots of the residuals by subject for the lmList fit to lm fit
# compare also the residual standard errors from both fits

fit2 <- lmList(height~age|Subject, data = dat) # use Oxboys instead !!!!!
summary(fit2)
boxplot(resid(fit2)~dat$Subject,horizontal=TRUE, xlab = "Residuals", ylab = "Subject")
# now all boxplots are centered around zero, because fitted to each subject so no longer have subject effect skewing residuals
# evidence that it was the subject effect that skewed the residuals

# fit1 Residual standard error: 8.081
# fit2 Residual standard error: 0.6598878
# residual standard errors is much lower for the second fit2, suggesting that fit2 fits the data better
# further evidence of significant subject effect


#-------------------
# Question 2(d)
#-------------------
# plot the individual confidence intervals on the parameters estimated by the lmList fit
# and verify that both the intercept and slope vary significantly with Subject

plot(intervals(fit2))
# can definitely see a difference between intercepts for diiferent subjects and the same for slope
# more evidence of a subject effect as well as a subject*age interaction



#-------------------
# Question 2(e)
#-------------------
# fit an LME model to the data with random effects for both the intercept and the slope
# examine the boxplots of the residuals by subject, comparing them to those obtained for the lm and lmList fits

fit3<-lme(fixed=height~age,random=~age|Subject,dat)
# adding slope random effect implies adding interaction effect between age (fixed effect) and subject (random effect)
summary(fit3)

boxplot(resid(fit3)~dat$Subject,horizontal=TRUE, xlab = "Residuals", ylab = "Subject")
# boxplots well centered around zero, since subject effect and subject*ge interaction effect have been accounted for in the model
# very similar to fit2 boxplots, which makes sense because fit2 had different slopes and intercepts for each subject, thus also accounting  for subject effect and subject*age interaction

# can also plot using
plot(fit3,Subject~resid(.))

#-------------------
# Question 2(f)
#-------------------
# produce the plot of the standardised residuals versus fitted values and the normal plot of standardizes residuals
# can you identify any departures from model's assumptions

qqnorm(fit3)
# some departure from normality at the tails, as to be expected
# right seems to show some departure from normality, but nothing serve
plot(fit3)
# mostly random scatter but there does seem to be a slight increase in spread as the fitted values increase, 
# suggesting that the assumption of common error variance could be violated
# and the model may not fit the data well
# although its really only slight


#-------------------
# Question 2(g)
#-------------------
# Plot the augemented predictions for the lme fit 
# Do the linear models for each subject appear adequate ?

plot(augPred(fit3)) # Doesn't work !

fit3<-lme(fixed=height~age,random=~age|Subject,Oxboys) 
plot(augPred(fit3)) # for some reason only works Oxboys data of class c("nfnGroupedData", "nfGroupedData" , "groupedData",    "data.frame"    )
# not a data frame 
# Even though both give the exact same output of the same class ?????

# Linear models for each subject seem to be very good, with each line passing through many points


#-------------------
# Question 2(h)
#-------------------
# Another way of assessing the linear models for each subject is to plot the residuals versus age by subject

plot(fit3,resid(.)~age|Subject)
# Several subjects have a noticeable "scooping" pattern in their residuals, indicating the need for a model with curvature



#-------------------
# Question 2(i)
#-------------------
# Use the lmList function to fit separate quadratic models for each subject
# The quadratic model is quadratic in age

fit4 <- lmList(height~age+I(age^2)|Subject,Oxboys)
summary(fit4)



#-------------------
# Question 2(j)
#-------------------
# examine a plot of the confidence intervals on coefficients from this second lmList fit
# are there indications that the coefficients differ between subjects
# are the quadratic coefficients significantly different from zero for some subject

plot(intervals(fit4))
# for all coefficients, especially the intercept and age coefficients, seem to differ from subject to subject
# coefficients for age^2 not so much, but still some evidence of difference
# a few of the quadratic coefficients are significantly different from zero, but most are not, most intervals containing zero, indicating that the age^2 term is probably unnecessary 


#-------------------
# Question 2(k)
#-------------------
# Fit the full mixed-effects model corresponding to the last lmList fit
# The model will have linear and quadratic terms for age in the fixed effects and the random effects
fit5<-lme(fit4)
summary(fit5)


#-------------------
# Question 2(l)
#-------------------
# check residual plots and numerical summaries for the lme model
# Do there appear to be deficiencies in the fit ?
# Do there appear to be terms in the morel that could be eliminated?

plot(fit5,Subject~resid(.))
# No obvious difference in residuals between subjects
# all centered around zero
# evidence of good model fit

qqnorm(fit5)
# again evidence of deviation from normality at the tales, as to be expected
# deviation is fairly pronounced though and over too large an area to just be in the tales
# so there seems to be some evidence of deviation from normaility

plot(fit5)
# scatter seems to be concentrated in the middle and becomes more dispersed further away
# so not really a random scatter
# thus possible deviation from the assumption of common error variance

plot(augPred(fit5))
# again fit seems pretty good for each subject, with no obvious discrepancies
# all lines are pretty close to straight lines, suggesting that the quadratic term may not be necessary and that linear is fine


plot(fit5,resid(.)~age|Subject)
# still presence of dips, despite addition of quadratic term

# so all in all, it looks like the quadratic term is unnecessary and a linear model is sufficient 

