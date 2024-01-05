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

#-------------------
# Sample
#-------------------

pids <- df$pid
head(pids)

uni_pids <- unique(pids)
str(uni_pids)

set.seed(1111996)
samp_pids <- sample(uni_pids, size=length(uni_pids)*0.75, replace = F)
str(samp_pids)

mydat <- df[df$pid %in% samp_pids, ]
str(mydat)


#---------------------------------
# Extract Variables of Interest
#---------------------------------

# remove gamedens, Pyrconcentration, Sulfconcentration, country, Studyyear,  PIoutcome, PIoutcomeday
mydatl <- mydat[,-c(6, 8, 9, 13, 14, 15, 16)]
str(mydatl)

# remove outlier
out <- which(mydatl$Hb==0)
mydatl <- mydatl[-out,]

# also want to include BASELINE parasite density
# first convert to wide format
mydatw <- reshape(mydatl,idvar= c("pid", "arm", "gender", "site", "age", "weight"), direction="wide",timevar="pday")
str(mydatw)
# log and store parasite density for day 0
mydatw$parbase <- log10(mydatw$pardens.0 + 1) 
str(mydatw)
# now remove unnecessary pardens variables 
mydatw <- mydatw[,-c(7,9,11,13,15,17,19,21,23)]
str(mydatw)

# now we need to transform back to long format
mydatl <- reshape(mydatw, idvar = "pid", varying = list(7:15), v.names = "Hb", direction = "long")
str(mydatl)
head(mydatl, 20)
# check: 
mydatl[mydatl$pid=="MOB2004_001",]$parbase   # same baseline value for same patient !

# remove NAs ( basically just getting rid of observations when Hg wasn't measured)
mydatl <- na.omit(mydatl)
str(mydatl)

# Note that there were certain days when Hg was NEVER measured
# looking at full data set mydat with no NAs removed
mydat[mydat$pday==0,]$Hb 
mydat[mydat$pday==1,]$Hb # All NAs
mydat[mydat$pday==2,]$Hb # All NAs
mydat[mydat$pday==3,]$Hb
mydat[mydat$pday==7,]$Hb
mydat[mydat$pday==14,]$Hb
mydat[mydat$pday==21,]$Hb # All NAs
mydat[mydat$pday==28,]$Hb
mydat[mydat$pday==42,]$Hb
# so actually ony have 6 days of measurements 

#-------------------
# Set Factors
#-------------------

str(mydatl)

mydatl$site <- as.factor(mydatl$site)  # 4 sites
levels(mydatl$site) <- c("South", "South", "North", "South")
levels(mydatl$site)
mydatl$arm <- as.factor(mydatl$arm)    # 2 arms
#mydatl$pid <- as.factor(mydatl$pid)    # 306 patients
mydatl$gender <- as.factor(mydatl$gender)  
mydatl$time <- as.factor(mydatl$time)  # 6 time points
levels(mydatl$time) <- c(0,3,7,14,28,42) # prefer to use actual days for time

str(mydatl)


# update wide data
mydatw <- reshape(mydatl,idvar= c("pid", "arm", "gender", "site", "age", "weight", "parbase"), direction="wide",timevar="time")
str(mydatw)


#--------------------------
# Change Response Variable
#---------------------------

# create haemoglobin baseline
mydatw$Hbbase <- mydatw$Hb.0
str(mydatw)
mydatw <- mydatw[,-8]

# add to long data
mydatl <- reshape(mydatw, idvar = "pid", varying = list(8:12), v.names = "Hb", direction = "long")
str(mydatl)

# fix time variable
mydatl$time <- as.factor(mydatl$time)  
levels(mydatl$time) <- c(3,7,14,28,42) - 3  # need to shift time to zero
mydatl <- na.omit(mydatl)
mydatl$time <- as.numeric(levels(mydatl$time))[mydatl$time]
str(mydatl)

# update wide data
mydatw <- reshape(mydatl,idvar= c("pid", "arm", "gender", "site", "age", "weight", "parbase", "Hbbase"), 
                  direction="wide",timevar="time")
str(mydatw)
# now only have 287 patients because removing any that only had measurements for day 0
# note that there will be NAs in mydatw because otherwise this dataframe would only contain  patients that were present every measurement day

#--------------------------
# Model 1: base model
#---------------------------
# first fit a model with only treatment and time as predictors

# Make a grouped data object
#groupdat1 <- groupedData(Hb~arm|pid, data=mydatl[,c(2,3,9,10)])
#lme(Hb~arm, data = groupdat1, method = "ML")
#lmer(Hb~(arm|pid), data = mydatl, REML = F,
#control = lmerControl(optimizer ="Nelder_Mead"))



mod1 <- lmer(Hb~time + arm  + (time|pid), data = mydatl, REML = F) 
summary(mod1)

mod1 <- lmer(Hb~time + arm +  (time|pid),data = mydatl,    # model with intercept and slope random effects
             REML = F,     # fit with ML to allow for comparison of models, will fit final model with REML               
             control = lmerControl(optimizer ="Nelder_Mead")) 
# using lmer to ensure correct p-values 
summary(mod1)

set.seed(120)
round(confint(mod1), 3)

# residuals versus fitted values
df1 <- data.frame(residuals = resid(mod1), 
                  fitted = predict(mod1))
ggplot(data = df1,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)

plot(mod1, resid(.,type="pearson")~fitted(.)|arm, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")


#--------------------------
# Model 2: base model - slope random effect
#---------------------------

mod2 <- lmer(Hb~time + arm +  (1|pid),data = mydatl,    # model with intercept and slope random effects
                   REML = F,     # fit with ML to allow for comparison of models, will fit final model with REML               
                   control = lmerControl(optimizer ="Nelder_Mead")) 
# using lmer to ensure correct p-values 
summary(mod2)


set.seed(120)
round(confint(mod2), 3)

# residuals versus fitted values
df2 <- data.frame(residuals = resid(mod2), 
                  fitted = predict(mod2))
ggplot(data = df2,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)

plot(mod2, resid(.,type="pearson")~fitted(.)|arm, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")


anova(mod2, mod1)
# keep random effect on slope !!!!




#--------------------------
# Model 3: choosing fixed effects
#---------------------------
#options(scipen=999)

#options(scipen=0)

mod3 <- lmer(Hb~time + arm + parbase  + Hbbase + gender + weight + age + site +  (time|pid), data = mydatl, 
             REML = F,  control = lmerControl(optimizer = "bobyqa")) 
summary(mod3)

mod4 <- lmer(Hb~time + arm + parbase  + Hbbase + gender + weight + age  +  (time|pid), data = mydatl, 
             REML = F, control = lmerControl(optimizer ="bobyqa")) 
summary(mod4)
anova(mod3, mod4) # choose mod4

mod5 <- lmer(Hb~time + arm + parbase  + Hbbase  + weight + age  +  (time|pid), data = mydatl, 
             REML = F, control = lmerControl(optimizer ="bobyqa")) 
summary(mod5)
anova(mod4, mod5) # AIC and likelihood-ratio test say mod4, BIC says mod5


# ADD INTERACTION TERMS:

mod6 <- lmer(Hb~time + arm*parbase  + arm*Hbbase + gender*weight + weight*age + parbase*Hbbase +  (time|pid), data = mydatl, 
             REML = F, control = lmerControl(optimizer ="Nelder_Mead")) 
summary(mod6)


mod7 <- lmer(Hb~time + arm*parbase + gender*weight + weight*age + parbase*Hbbase +  (time|pid), data = mydatl, 
             REML = F, control = lmerControl(optimizer ="Nelder_Mead")) 
summary(mod7)
anova(mod6, mod7) # choose mod7


mod8 <- lmer(Hb~time + arm*parbase + gender*weight + weight*age + Hbbase +  (time|pid), data = mydatl, 
             REML = F, control = lmerControl(optimizer = "bobyqa")) 
summary(mod8)
anova(mod8, mod7) # choose mod8


mod9 <- lmer(Hb~time + arm*parbase + gender*weight + age + Hbbase +  (time|pid), data = mydatl, 
             REML = F, control = lmerControl(optimizer = "Nelder_Mead")) 
summary(mod9)
anova(mod8, mod9) # choose mod9


mod10 <- lmer(Hb~time + arm + parbase + gender*weight + age + Hbbase +  (time|pid), data = mydatl, 
              REML = F, control = lmerControl(optimizer = "bobyqa")) 
summary(mod10)
anova(mod10, mod9) # choose mod10


mod11 <- lmer(Hb~time + arm + parbase + weight + gender + age + Hbbase +  (time|pid), data = mydatl, 
              REML = F, control = lmerControl(optimizer = "bobyqa")) 
summary(mod11)
anova(mod10, mod11) # choose mod10


mod13 <- lmer(Hb~time + arm + parbase + gender*weight + Hbbase +  (time|pid), data = mydatl, 
              REML = F, control = lmerControl(optimizer = "Nelder_Mead")) 
summary(mod13)
anova(mod13, mod10) # choose mod10


mod14 <- lmer(Hb~time + arm + parbase + weight + age + Hbbase +  (time|pid), data = mydatl, 
              REML = F, control = lmerControl(optimizer = "bobyqa")) 
summary(mod14)
anova(mod10, mod14) # choose mod10




set.seed(120)
con <- confint(mod10)
round(con, 3)



#--------------------------
# Choosing variance function
#---------------------------

# residuals versus fitted values
df10 <- data.frame(residuals = resid(mod10), 
                  fitted = predict(mod10))
ggplot(data = df10,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)

plot(mod10, resid(.,type="pearson")~fitted(.)|arm, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")

plot(mod10, resid(.,type="pearson")~fitted(.)|gender, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black", ylim=c(-4,4))


# Parbase
dfa <- data.frame(residuals = resid(mod10,type="pearson"), 
                   fitted = mydatl$parbase)
ggplot(data = dfa,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Logged baseline parasite density")+
  theme_minimal(base_size = 22)

# Weight
dfb <- data.frame(residuals = resid(mod10,type="pearson"), 
                  fitted = mydatl$weight)
ggplot(data = dfb,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Weight")+
  theme_minimal(base_size = 22)

# Age
dfc <- data.frame(residuals = resid(mod10,type="pearson"), 
                  fitted = mydatl$age)
ggplot(data = dfc,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Age")+
  theme_minimal(base_size = 22)

# Hbbase
dfd <- data.frame(residuals = resid(mod10,type="pearson"), 
                  fitted = mydatl$Hbbase)
ggplot(data = dfd,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Baseline haemoglobin")+
  theme_minimal(base_size = 22)



#--------------------------
# Choosing correlation structure
#---------------------------



# Re-fit Model 10 using nlme:
mod10_refit <- nlme::lme(Hb~time + arm + parbase + gender*weight + age + Hbbase,
                         data = mydatl, 
                         random = ~1 + time|pid)
summary(mod10_refit)  

ACF(mod10_refit)

mod_lme <- nlme::lme(Hb~time + arm + parbase + gender*weight + age + Hbbase,
                     data = mydatl, 
                     random = ~1 + time|pid,
                     correlation = corAR1())
summary(mod_lme)  


anova(mod10_refit, mod_lme) # mod10_refit is nested within mod_lme (corresponding to phi =0)
# insignificant p-value implies independent errors better (mod10_refit)



qqnorm(mod10_refit, col="black", pch=20)


mod13 <- nlme::lme(Hb~time + arm + parbase + gender*weight + age + Hbbase,
                   data = mydatl, 
                   random = list(pid=pdDiag(~time)))
summary(mod13)
anova(mod13,mod10_refit)


qqnorm(mod10_refit, ~ranef(.), col="black", pch=20)

pairs(mod10_refit, ~ranef(.)|arm, col="black", pch=20)
pairs(mod10_refit, ~ranef(.)|gender, col="black", pch=20)


#--------------------------
# FINAL MODEL
#---------------------------
# lmerTest fit with REML
final <- lmer(Hb~time + arm + parbase + gender*weight + age + Hbbase +  (time|pid), data = mydatl, 
              REML = T, control = lmerControl(optimizer = "bobyqa")) 
summary(final)

set.seed(120)
round(confint(final), 3)

plot(final, Hb~fitted(.), id=0.005)



df <- data.frame(Hb = mydatl$Hb, 
                  fitted = fitted(final))
ggplot(data = df,
       aes(x = fitted, y = Hb)) +
  geom_point(size=2) + 
  ylab("Observed haemoglobin values")+
  xlab("Fitted haemoglobin values")+
  theme_minimal(base_size = 22)+
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1.5)










