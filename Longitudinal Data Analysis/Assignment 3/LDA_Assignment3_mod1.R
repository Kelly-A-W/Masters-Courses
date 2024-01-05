# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(nlme)
library(readxl)
library(ggcorrplot)
library(RColorBrewer)
library(lmerTest)

#-------------------
# Load Data
#-------------------

df <- read_excel("malariadata.xls", sheet = 1)
df <- as.data.frame(df)
head(df)
str(df)

#---------------------------------
# Extract Variables of Interest
#---------------------------------

# remove, country, Studyyear, PIoutcomeday
mydatl <- df[,-c(13, 14, 16)]
str(mydatl)

#--------------------------------------
# Remove Day 1, 2 and 21 Observations 
#--------------------------------------
# gamedens was not collected on days 1, 2 and 21, so remove all observations for these days
day1 <- which(mydatl$pday==1)
day2 <- which(mydatl$pday==2)
day21 <- which(mydatl$pday==21)

mydatl <- mydatl[-c(day1,day2,day21),]
str(mydatl)

#-------------------
# Set Factors
#-------------------

mydatl$site <- as.factor(mydatl$site)  # 4 sites
mydatl$arm <- as.factor(mydatl$arm)    # 2 arms
mydatl$gender <- as.factor(mydatl$gender)  
str(mydatl)

#--------------------------------------
# Only include patients in SP arm
#--------------------------------------

spart <- which(mydatl$arm=="SP/ART")
mydatl <- mydatl[-spart,]
mydatl <- mydatl[,-2]  # don't need treatment variable anymore

#--------------------------
# Variable Transformations
#--------------------------

mydatl$pardens <- log10(1+mydatl$pardens)
colnames(mydatl)[4] <- "lpardens"
mydatl$Pyrconcentration <- log10(1+mydatl$Pyrconcentration)
colnames(mydatl)[7] <- "lpyrconcentration"
mydatl$Sulfconcentration <- log10(1+mydatl$Sulfconcentration)
colnames(mydatl)[8] <- "lsulfconcentration"
str(mydatl)

# Wide format:
mydatw <- reshape(mydatl,idvar= c("pid", "gender", "site", "age", "weight", "PIoutcome"), direction="wide",timevar="pday")
str(mydatw)


# need log transformation
mydatl$gamedens<- log10(1+mydatl$gamedens)
colnames(mydatl)[5] <- "lgamedens"
head(mydatl)

# update wide data sets
mydatw <- reshape(mydatl,idvar= c("pid", "gender", "site", "age", "weight", "PIoutcome"), direction="wide",timevar="pday")
str(mydatw)




#----------------------------------------------------------------
# Model 1: base model with random effects on intercept AND slope
#---------------------------------------------------------------
mydatl <- na.omit(mydatl)

mod1 <- lmer(lgamedens~pday + I(pday^2) + site  + (pday|pid), data = mydatl, REML = F) 
mod1 <- lmer(lgamedens~pday + I(pday^2) + site  + (pday|pid), data = mydatl, REML = F, 
             control = lmerControl(optimizer ="Nelder_Mead")) 
mod1 <- lmer(lgamedens~pday + I(pday^2) + site  + (pday|pid), data = mydatl, REML = F, 
             control = lmerControl(optimizer ="bobyqa")) 
mod1 <- lmer(lgamedens~pday + I(pday^2) + site  + (pday|pid), data = mydatl, REML = F, 
             control = lmerControl(optimizer ="nlminbwrap")) 
# tried all 4 optimisers but got a singular fit each time
# so the model is probably overfitted
summary(mod1)

# residuals versus fitted values
df1 <- data.frame(residuals = resid(mod1), 
                  fitted = predict(mod1))
ggplot(data = df1,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)

plot(mod1, resid(.,type="pearson")~fitted(.)|site, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")


#-------------------------------------------------------------
# Model 2: base model with random effects on intercept ONLY
#------------------------------------------------------------

mod2 <- lmer(lgamedens~pday + I(pday^2) + site  + (1|pid), data = mydatl, REML = F) 
summary(mod2)


# residuals versus fitted values
df2 <- data.frame(residuals = resid(mod2), 
                  fitted = predict(mod2))
ggplot(data = df2,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)

plot(mod2, resid(.,type="pearson")~fitted(.)|site, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")


anova(mod2, mod1)
# choose simpler model: mod2




#--------------------------
# Model 3: choosing fixed effects
#---------------------------
#options(scipen=999)

#options(scipen=0)

mod3 <- lmer(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + lpyrconcentration + 
               lsulfconcentration + gender + weight + age +(1|pid), data = mydatl, REML = F) 
summary(mod3)

mod4 <- lmer(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + lpyrconcentration + 
               lsulfconcentration + gender + weight  +(1|pid), data = mydatl, REML = F) 
summary(mod4)
anova(mod3, mod4) # choose mod4

mod5 <- lmer(lgamedens~pday + I(pday^2) + site  + lpardens + Hb  + 
               lsulfconcentration + gender + weight  +(1|pid), data = mydatl, REML = F) 
summary(mod5)
anova(mod4, mod5) # mod5

mod6 <- lmer(lgamedens~pday + I(pday^2) + site  + lpardens + Hb  + gender + weight  +(1|pid), 
             data = mydatl, REML = F) 
summary(mod6)
anova(mod5, mod6) # choose mod6

mod7 <- lmer(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + weight  +(1|pid), 
             data = mydatl, REML = F) 
summary(mod7)
anova(mod6, mod7) # choose model 7

mod8 <- lmer(lgamedens~pday + I(pday^2) + site  + lpardens + weight  +(1|pid), 
             data = mydatl, REML = F) 
summary(mod8)
anova(mod7, mod8) # choose mod8

mod9 <- lmer(lgamedens~pday + I(pday^2) + site   + weight  +(1|pid), 
             data = mydatl, REML = F) 
summary(mod9)
anova(mod8, mod9)  # CHOOSE MOD8



# ADD INTERACTION TERMS:

mod10 <- lmer(lgamedens ~ pday*site + I(pday^2)  + lpardens*weight  + (1|pid), data = mydatl, REML = F) 
summary(mod10)

mod11 <- lmer(lgamedens ~ pday + site + I(pday^2)  + lpardens*weight  + (1|pid), data = mydatl, REML = F) 
summary(mod11)
anova(mod10, mod11)  # choose mod 11
anova(mod11, mod8)   # choose mod 11



#--------------------------
# Choosing variance function
#---------------------------

# residuals versus fitted values
df11 <- data.frame(residuals = resid(mod11), 
                   fitted = predict(mod11))
ggplot(data = df11,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)

plot(mod11, resid(.,type="pearson")~fitted(.)|site, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")


# Parbase
dfa <- data.frame(residuals = resid(mod11,type="pearson"), 
                  fitted = mydatl$lpardens)
ggplot(data = dfa,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Logged baseline parasite density")+
  theme_minimal(base_size = 22)

# Weight
dfb <- data.frame(residuals = resid(mod11,type="pearson"), 
                  fitted = mydatl$weight)
ggplot(data = dfb,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Weight")+
  theme_minimal(base_size = 22)





# CANNOT use a variance function in lmer ! Have to use lme
mydatl <- groupedData(lgamedens~pday|pid, data = mydatl)
mod20 <- lme(lgamedens ~ pday + site + I(pday^2)  + lpardens*weight, 
             random = ~1, data = mydatl, method="ML")
mod21 <- update(mod20, weights = varPower(fixed=0.5))
anova(mod20,mod21)
mod22 <- update(mod20, weights = varFixed(~ 1/weight))
anova(mod20,mod22)
mod23 <- update(mod20, weights = varExp())
anova(mod20,mod21,mod22,mod23)




# residuals versus fitted values
df11 <- data.frame(residuals = resid(mod23), 
                   fitted = predict(mod23))
ggplot(data = df11,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)
plot(mod11, resid(.,type="pearson")~fitted(.)|site, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")
# Parbase
dfa <- data.frame(residuals = resid(mod23,type="pearson"), 
                  fitted = mydatl$lpardens)
ggplot(data = dfa,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Logged baseline parasite density")+
  theme_minimal(base_size = 22)
# Weight
dfb <- data.frame(residuals = resid(mod23,type="pearson"), 
                  fitted = mydatl$weight)
ggplot(data = dfb,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Weight")+
  theme_minimal(base_size = 22)



#--------------------------
# Choosing correlation structure
#---------------------------


ACF(mod23)

mod31 <- update(mod23,  correlation = corAR1())
anova(mod23,mod31)
mod32 <- update(mod23,  correlation = corCompSymm())
anova(mod23,mod32)
mod33 <- update(mod23,  correlation = corSymm())
anova(mod23,mod33)



#--------------------------
# FINAL MODEL
#---------------------------

final <- update(mod23, method = "REML")
summary(final)

set.seed(120)
round(intervals(final), 3)

qqnorm(final, col="black", pch=20)


qqnorm(final, ~ranef(.), col="black", pch=20)


plot(final, lgamedens ~fitted(.))


df <- data.frame(lgamedens  = mydatl$lgamedens , 
                 fitted = fitted(final))
ggplot(data = df,
       aes(x = fitted, y = lgamedens )) +
  geom_point(size=2) + 
  ylab("Observed values")+
  xlab("Fitted  values")+
  theme_minimal(base_size = 22)
  #geom_abline(intercept = 0, slope = 1, color="red", 
              #linetype="dashed", size=1.5)






library(pROC)
library(ROCR)
predict_value <- as.numeric(predict(final) )
actual_values <- mydatl$lgamedens

ROCR_pred_test <- prediction(predict_value,actual_values)
ROCR_perf_test <- performance(ROCR_pred_test,'tpr','fpr')
plot(ROCR_perf_test,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))

