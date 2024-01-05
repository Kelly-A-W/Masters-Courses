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

#------------------------------------
# Deal with Zero-Inflation
#------------------------------------

pids <- mydatl$pid
uni_pids <- unique(pids)
str(uni_pids) # 198 patients


samp_pids <- c()
for(i in 1:length(uni_pids)){
  rows <- which(mydatl$pid==uni_pids[i])
  gam <- mydatl[rows,]$lgamedens
  gam <- na.omit(gam)
  add <- sum(gam != 0)
  
  if(add > 1){
    samp_pids <- append(samp_pids, uni_pids[i])
  }
}

str(samp_pids) # 372 patients
new_dfl <- mydatl[mydatl$pid %in% samp_pids, ]
str(mydatl)
nrow(new_dfl[new_dfl$site=="Boane",])
nrow(new_dfl[new_dfl$site=="Catuane",])
nrow(new_dfl[new_dfl$site=="Magude",])
nrow(new_dfl[new_dfl$site=="Namaacha",])
levels(new_dfl$site) <- c("Boane", "Boane", "Magude", "Namaacha")
levels(new_dfl$site)


ggplot(new_dfl,aes(x=lgamedens)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  xlab("Gametocyte Density")+
  theme_minimal(base_size = 27)


mydatl <- new_dfl
# Wide format:
mydatw <- reshape(mydatl,idvar= c("pid", "gender", "site", "age", "weight", "PIoutcome"), direction="wide",timevar="pday")
str(mydatw)




#----------------------------------------------------------------
# Model 1: base model with random effects on intercept AND slope
#---------------------------------------------------------------
mydatl <- na.omit(mydatl)
mydatl <- groupedData(lgamedens~pday|pid, data = mydatl)

mod1 <- lme(lgamedens~pday + I(pday^2) + site, random = list(~pday|pid), method="ML", data = mydatl)

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
mod2 <- lme(lgamedens~pday + I(pday^2) + site, random = list(~1|pid), method="ML", data = mydatl)
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




#-------------------------------------
# Model 3: choosing fixed effects
#-------------------------------------
#options(scipen=999)

#options(scipen=0)
mod3 <- lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + lpyrconcentration + 
              lsulfconcentration + gender + weight + age, random = list(~1|pid), method="ML", data = mydatl)
summary(mod3)


mod4 <- lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + lpyrconcentration +
              gender + weight + age, random = list(~1|pid), method="ML", data = mydatl)
summary(mod4)
anova(mod3, mod4) # choose mod4

mod5 <- lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + lpyrconcentration +
              gender + weight, random = list(~1|pid), method="ML", data = mydatl)
summary(mod5)
anova(mod4, mod5) # mod5

mod6 <- lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + lpyrconcentration +
              gender, random = list(~1|pid), method="ML", data = mydatl)
summary(mod6)
anova(mod5, mod6) # choose mod6

mod7 <- lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb +
              gender, random = list(~1|pid), method="ML", data = mydatl)
summary(mod7)
anova(mod6, mod7) # choose model 7


mod8 <- lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb, 
            random = list(~1|pid), method="ML", data = mydatl)
anova(mod7, mod8) 
mod8 <- lme(lgamedens~pday + I(pday^2) + site  + lpardens  +
              gender, random = list(~1|pid), method="ML", data = mydatl)
anova(mod7, mod8)
mod8 <- lme(lgamedens~pday + I(pday^2) + site   + Hb +
              gender, random = list(~1|pid), method="ML", data = mydatl)
anova(mod7, mod8)
# CHOOSE MOD7


# ADD INTERACTION TERMS:

mod9 <- lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + gender +
              pday*site + pday*gender + lpardens*Hb, random = list(~1|pid), method="ML", data = mydatl)
summary(mod9)

mod10 <- lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + gender +
              pday*site  + lpardens*Hb, random = list(~1|pid), method="ML", data = mydatl)
summary(mod10)
anova(mod9, mod10)  # choose mod10

mod11 <- lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + gender +
               pday*site, random = list(~1|pid), method="ML", data = mydatl)
summary(mod11)
anova(mod10, mod11)  # choose mod11

anova(mod11, mod7)   # choose mod7




#--------------------------
# Choosing variance function
#---------------------------

# residuals versus fitted values
df11 <- data.frame(residuals = resid(mod7), 
                   fitted = predict(mod7))
ggplot(data = df11,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)

plot(mod7, resid(.,type="pearson")~fitted(.)|site, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")
plot(mod7, resid(.,type="pearson")~fitted(.)|gender, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")
# Parbase
dfa <- data.frame(residuals = resid(mod7,type="pearson"), 
                  fitted = mydatl$lpardens)
ggplot(data = dfa,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Logged baseline parasite density")+
  theme_minimal(base_size = 22)
# HB
dfb <- data.frame(residuals = resid(mod7,type="pearson"), 
                  fitted = mydatl$Hb)
ggplot(data = dfb,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Haemoglobin")+
  theme_minimal(base_size = 22)




mod21 <- update(mod7, weights = varPower(fixed=0.5))
mod22 <- update(mod7, weights = varExp())
mod23 <- update(mod7, weights = varIdent(form = ~1|gender))
anova(mod7,mod21, mod22, mod23)





# residuals versus fitted values
df11 <- data.frame(residuals = resid(mod21), 
                   fitted = predict(mod21))
ggplot(data = df11,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)
plot(mod21, resid(.,type="pearson")~fitted(.)|site, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")
plot(mod21, resid(.,type="pearson")~fitted(.)|gender, id=0.0005, adj=-0.3, 
     xlab = "Fitted values", ylab = "Pearson residuals", pch = 20, col= "black")
# Parbase
dfa <- data.frame(residuals = resid(mod21,type="pearson"), 
                  fitted = mydatl$lpardens)
ggplot(data = dfa,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Logged baseline parasite density")+
  theme_minimal(base_size = 22)
# HB
dfb <- data.frame(residuals = resid(mod21,type="pearson"), 
                  fitted = mydatl$Hb)
ggplot(data = dfb,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Haemoglobin")+
  theme_minimal(base_size = 22)




#--------------------------
# Choosing correlation structure
#---------------------------


ACF(mod21)

mod31 <- update(mod21,  correlation = corAR1())
mod32 <- update(mod21,  correlation = corCompSymm())
mod33 <- update(mod21,  correlation = corSymm())
anova(mod21,mod31, mod32, mod33) 



#--------------------------
# FINAL MODEL
#---------------------------

final <- update(mod32, method = "REML", control=lmeControl(tolerance = 0.001))
summary(final)

set.seed(120)
round(intervals(final), 3)

df11 <- data.frame(residuals = resid(final), 
                   fitted = predict(final))
ggplot(data = df11,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)



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
  theme_minimal(base_size = 22)+
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5)










