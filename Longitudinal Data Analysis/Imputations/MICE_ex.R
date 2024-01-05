# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(mice)
library(VIM)
library(lattice)

# Load data
data("nhanes")

#----------------------------------------
# STEP 0: Assess Structure of Missingness
#----------------------------------------

head(nhanes)

md.pattern(nhanes)  # shows the different possible missingness patterns possible and the number of observations 
                    # with that pattern
                    # eg: first row shows that 13 observations have no missing measurements 

md.pairs(nhanes)    # rr for both present, rm for column variable missing, row variable present, mr is the transpose 
                    # of rm and mm is both missing

marginplot(nhanes[,c("chl", "bmi")], 
           col = mdc(1:2), cex=1.2, cex.lab=1.2, cex.numbers = 1.3, pch=19)
# do see a difference between the blue and pink distributions, suggesting that missingness may not be MCAR
# hyp seems to be categorical, so won't check plot with hyp



# HELP CHOOSING PREDICTORS

# Pearson correlations:
round(cor(nhanes, use="pair"), 3) # highly correlated variables will make better predictors

# Correlations between variables and missing data indicators:
round(cor(y=nhanes, x=!is.na(nhanes), use="pair"), 3)  # high correlation make beter predictors

# Proportion of usable cases:
p <- md.pairs(nhanes)
round(p$mr/(p$mr+p$mm), 3) # measures how many cases on the target variable have observations on the predictor
# high makes good predictors 

# Guessing predictor matrix based off correlations:
quickpred(nhanes)
quickpred(nhanes, minpuc = 0.25, include = "age") # override default parameters


#----------------------------------------
# STEP 1: MICE Imputation
#----------------------------------------

imp1 <- mice(nhanes, seed=23109)  # MICE imputation using DEFAULT IMPUTATION MODEL pmm
# imputation method is reflective of the (what you think are ) underlying distributions of the variables 
                                  # M = number of iterations = number of imputed datasets = 5

print(imp1) # shows what imputation methods were used for which variables (age had no missing values so was not imputed)
            # gives the PREDICTOR MATRIX which shows which variables were used to impute which
            # here every variable was used to impute every other variable (default)
imp1$imp$bmi # gives all M = 5 imputations for bmi variable
complete(imp1,2)  # the complete dataset where NAs have been fulled in using the second imputed data set
complete(imp1, "long") # complete datasets for all 5 imputations in long format
complete(imp1, "long", include = T) # will include the original missing data dataset as well at the top


# INCREASE NUMBER OF ITERATIONS M
imp2 <- mice(nhanes2, seed=23109,
             m=50)  # number of iterations to minimize simulation error


# CHOOSE IMPUTATION MODEL
imp3 <- mice(nhanes2, seed=23109,
             meth = c("", "pmm", "logreg", "norm")) # choice of imputation model depends on variable type (see slide 11)
                                                    # "" used to suppress imputation of age variable (which had no missing values)
                                                    # note that you can skip imputation of any variable using "", even ones with NAs

# SPECIFY PREDICTORS TO BE USED IN IMPUTATION
imp0 <- mice(nhanes, maxit = 0)   # create a "shell" (could just use imp1)
(pred <- imp0$predictorMatrix)  
pred[,"bmi"] <- 0  # editing matrix so that bmi is not used to predict any variable
pred
imp4a <- mice(nhanes, seed=23109,
             pred = pred) 
imp4a$imp$bmi
# Alternatively, 
(pred <- imp0$predictorMatrix)  
pred["bmi",] <- 0  # now no variables are used to predict bmi !
pred
imp4b <- mice(nhanes, seed=23109,
             pred = pred)   
imp4b$imp$bmi


# IMPUTATION FOR LONGITUDINAL DATA
head(popmis)  # longitudinal data: dependent variable = popular, random effects = const (intercept) and sex (slope), fixed effect = texp, class variable = school
              # only missing values for popular; pupil is nested within school
imp0 <- mice(popmis, maxit = 0)   # create shell
(pred <- imp0$predictorMatrix)
pred["popular",] <- c(0,-2,0,2,1,2,0)  # specifying predictor variables for popular
                                       # -2 for class variable, 2 for random effect
# We need to do all this because mice uses mixed effect models for longitudinal data
# we could also convert data to wide format and impute that, then we dont have do all this
imp5 <- mice(popmis, 
             meth = c("", "", "2l.norm", "", "", "", ""),  # only popular had missing values so only need to impute popular
             pred = pred,
             maxit = 1, seed = 71152)


# PASSIVE IMPUTATION 

# Passive imputation to PRESERVE TRANSFORMATIONS:
nhanes2.ext <- cbind(nhanes2, lchl = log(nhanes2$chl))
imp0 <- mice(nhanes2.ext, max=0)
meth <- imp0$method
meth["lchl"] <- "~log(chl)"    # specifying PASSIVE IMPUTATION of lchl (=log(chl)) to ensure that imputation of lchl
                               # is based on the latest imputation of chl
pred <- imp0$predictorMatrix
pred[c("hyp","chl"), "lchl"] <- 0  # don't use lchl to impute hyp or chl
pred["bmi", "chl"] <- 0    # don't use chl to impute bmi
pred
imp6 <- mice(nhanes2.ext, seed=38788, print = F,
             meth=meth, 
             pred=pred)
# We can then take it a step further by using the SQUEEZE FUNCTION prevent implausible lchl values:
meth["lchl"] <- "~log(squeeze(chl, bounds=c(100,300)))"  # keeping chl in bounds to allow for computation of logged values

# Using passive imputation for MORE THAN ONE VARIABLE:
head(boys)
imp0 <- mice(boys, max=0)
meth <- imp0$method
meth["bmi"] <- "~I(wgt/(hgt/100)^2)"  # calculate bmi from imputed values of weight and height 
# I just means that we're putting in a formula, not an indicator
pred <- imp0$predictorMatrix
pred[c("wgt", "hgt", "hc", "reg"), "bmi"] <- 0  # don't use bmi to impute variables wgt, hgt, hc or reg
pred[c("gen", "phb", "tv"), c("hgt", "wgt", "hc")] <- 0  # don't use hgt, wgt or hc to impute gen, phb or tv
imp6 <- mice(boys, maxit=20, seed = 9212, print=F, 
             pred=pred, 
             meth=meth)
head(complete(imp6)[is.na(boys$bmi),], 3)

# Using passive imputation for dealing with INTERACTIONS
# do this if you want to include interaction terms in you model
# note that you will NOT use this "bmi.chl" variable instead of a normal interaction term in your model
# including "bmi.chl" in imputation just helps effect convergence of your imputations and make sure you imputations
# are appropriate for including an interaction term in your model
nhanes2.ext <- cbind(nhanes2, bmi.chl=NA)
imp0 <- mice(nhanes2.ext, max=0, print=F)
meth <- imp0$method
meth["bmi.chl"] <- "~I((bmi-25)*(chl-200))"
pred <- imp0$predictorMatrix
pred[c("bmi", "chl"), "bmi.chl"] <- 0
imp7 <- mice(nhanes2.ext, seed = 9212, print=F, 
             pred=pred, 
             meth=meth)


# ORDER OF IMPUTATIONS
# set the order of your imputations to speed up convergence or to ensure meaningful imputations if variables depende on one another
imp8a <- mice(nhanes2.ext, meth=meth, pred=pred, 
              vis = c(2,4,5,3)) # order of which variables to impute
imp8b <- mice(nhanes2.ext, meth=meth, pred=pred, 
              vis = "monotone") # imputing variables increasing order of the number of missing data




#----------------------------------------
# STEP 2: Assess Imputations
#----------------------------------------

# scatter-plots of the observed and the imputed values for each variable
stripplot(imp1, pch=20,cex=1.2) 
# distribution of the imputed values should not be very different from that of the observed values

# scatter-plots of bmi vs chl for observed data and for the M imputed datasets
xyplot(imp1, bmi~chl|.imp, pch=20, cex=1.4)
# we want the all the scatter-plots to look alike
# (can do for other variables but it in our case chl and bmi are the only continuous variables with imputations)
# dont want our imputations to change the relationship between our variables

# CHECKING CONVERGENCE OF IMPUTATIONS
plot(imp6, c("hgt", "wgt", "bmi")) # good mixing of patterns indicates good convergence
densityplot(imp6)  # want densities to overlap and be similar to density of observed values
# hgt and tv imputations don't look great 

#-----------------------------------------------------------------
# STEP 3: Extract Parameter Estimates (for each imputed dataset)
#-----------------------------------------------------------------

fit1 <- with(imp1, lm(chl~age+bmi)) # in this case our model of choice is a simple linear regression model
fit1
# model fitted to each imputed dataset, resulting in M sets of estimates for each parameter
fit1$analyses[[3]] # analysis of third imputed dataset

# we can actually insert any calculation into with(), even one made by us!
expr <- expression(ov <- cut(bmi, c(10,25,50)), table(age,ov))
fit2 <- with(imp1, eval(expr))
fit2$an[c(2,5)] # contingency tables for 2nd and 5th imputation data sets

#--------------------------------------------------------
# STEP 4: Combine Parameter Estimates using Rubin's Rule
#--------------------------------------------------------

print(pool(fit1))
# fmi = fraction of missing information
# lambda = proportion of total variance attributable to missing data

summary(pool(fit1))


#--------------------------------------------------------
# STEP 5: Model Comparison
#--------------------------------------------------------
#number iterations in imputation is guided by convergence
# number of imputed data sets, 5 generally used, more might be better
imp <- mice(nhanes2, print=F, m=50, seed=219)
fit0 <- with(data=imp, expr=lm(bmi~age+hyp))
fit1 <- with(data=imp, expr=lm(bmi~age+hyp+chl))

# WALD TEST:
wt <- D1(fit1, fit0)  # D1 is a multivariate Wald Test to compare two nested models (pool.compare() function has been removed from r)
wt # gives a very different p-value to notes........

imp <- mice(boys, print=F, seed=60019)
fit0 <- with(data=imp, expr = glm(I(gen>levels(gen)[1])~hgt+hc, family = binomial))
fit1 <- with(data=imp, expr = glm(I(gen>levels(gen)[1])~hgt+hc+reg, family = binomial))

# LIKELIHOOD RATIO TEST
lrt <- D3(fit1, fit0)
lrt # also different to notes.... BUT same p-value that i got when using pool.compare



#--------------------------------------------------------
# STEP 6: Sensitivity Analysis
#--------------------------------------------------------
# investigate the sensitivity of the results to different values of a variable, in this case chl
# Suppose we suspect that our imputations under MAR assumption are too low
# We can check this by multiplying the imputations by a factor
ini  <- mice(nhanes2, maxit=0, print=F) # create shell
post <- ini$post # A vector of strings of length length(blocks) with commands for post-processing for each variable
k    <- seq(1, 1.5, 0.1)
est  <- vector("list", length(k))
for(i in 1:length(k)) {
  post["chl"] <- paste("imp[[j]][,i] <- ", k[i], "* imp[[j]][,i]")  
  imp <- mice(nhanes2, post=post, seed=10, print=F, maxit=20)
  fit <- with(imp, lm(bmi~age+chl))
  est[[i]] <- summary(pool(fit))
}
est    # scenarios involving increasing imputations for chl by 0% (MAR), 10%, 20%, 30%, 40% and 50% (MNAR)
# here we see that the intercept changes but the other effects less so
# and we dont mind that the intercept changes because we are more interested in the other estimates being consistent and accurate
# because we are more interested in the relationships between variables













