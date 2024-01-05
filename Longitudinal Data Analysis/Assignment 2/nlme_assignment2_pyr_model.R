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
library(latex2exp)

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

# remove gamedens, Hb, country, Studyyear,  PIoutcome, PIoutcomeday
mydatl <- df
mydatl <- mydatl[,-c(6, 7, 13, 14, 15, 16)]
str(mydatl)

# also want to include BASELINE parasite density
# first convert to wide format
mydatw <- reshape(mydatl,idvar= c("pid", "arm", "gender", "site", "age", "weight"), direction="wide",timevar="pday")
str(mydatw)
# log and store parasite density for day 0
mydatw$parbase <- log10(mydatw$pardens.0 + 1) 
str(mydatw)
# now remove unnecessary pardens variables 
mydatw <- mydatw[,-c(7,10,13,16,19,22,25,28,31)]
str(mydatw)

# now we need to transform back to long format
mydatl <- reshape(mydatw, idvar = "pid", 
                  varying = list(Pyrconcentration=c(7,9,11,13,15,17,19,21,23), Sulfconcentration=c(8,10,12,14,16,18,20,22,24)), 
                  v.names = c("Pyrconcentration","Sulfconcentration"), 
                  direction = "long", 
                  times = c(0,1,2,3,7,14,21,28,42), # setting day numbers
                  timevar = "day") # naming the time column
str(mydatl)
head(mydatl, 20)
# check: 
mydatl[mydatl$pid=="MOB2004_001",]$parbase   # same baseline value for same patient !


#-------------------
# Set Factors
#-------------------

mydatl$site <- as.factor(mydatl$site)  # 4 sites
# Collapse sites into N and S
levels(mydatl$site) <- c("South", "South", "North", "South")
levels(mydatl$site)                    # 2 sites

mydatl$arm <- as.factor(mydatl$arm)    # 2 arms
mydatl$gender <- as.factor(mydatl$gender)  

str(mydatl)


#-------------------
# Sample
#-------------------
# only keep subjects with at least 4 repeated measurements

# need to create two datasets
dfl_sulf <- mydatl[,-9]
dfl_pyr <- mydatl[,-10]

# Remove rows with missing repeated measures
# check if covariates have NAs: (should be none!)
sum(is.na(mydatl$site))
sum(is.na(mydatl$arm))
sum(is.na(mydatl$pid))
sum(is.na(mydatl$gender))
sum(is.na(mydatl$age))
sum(is.na(mydatl$parbase))
sum(is.na(mydatl$day))
sum(is.na(mydatl$weight)) # some weights are missing !
# so need to be careful about how we remove NA's
miss1 <- which(is.na(dfl_sulf$Sulfconcentration))  # rows (observations) for which Sulfconcentration == NA
miss2 <- which(is.na(dfl_pyr$Pyrconcentration))
dfl_sulf <- dfl_sulf[-miss1,]  # removing those rows with missing values
dfl_pyr <- dfl_pyr[-miss2,]

# Now we only want to keep subjects with at least 4 repeated measurements
# i.e. only keep subject pids that are repeated more than 3 times in the dataset 
pids1 <- dfl_sulf$pid
pids2 <- dfl_pyr$pid
uni_pids1 <- unique(pids1)
uni_pids2 <- unique(pids2)
str(uni_pids1) # 408 patients
str(uni_pids2) # 408 patients

samp_pids1 <- c()
for(i in 1:length(uni_pids1)){
  if(length(pids1[pids1==uni_pids1[i]])>3){
    samp_pids1 <- append(samp_pids1, uni_pids1[i])
  }
}
str(samp_pids1) # 372 patients
dfl_sulf <- dfl_sulf[dfl_sulf$pid %in% samp_pids1, ]
str(dfl_sulf)

samp_pids2 <- c()
for(i in 1:length(uni_pids2)){
  if(length(pids2[pids2==uni_pids2[i]])>3){
    samp_pids2 <- append(samp_pids2, uni_pids2[i])
  }
}
str(samp_pids2) # 372
dfl_pyr <- dfl_pyr[dfl_pyr$pid %in% samp_pids2, ]
str(dfl_pyr)

# can look at the intersection of samp_pids to see how many patients have four observations for both variables
length(intersect(samp_pids1,samp_pids2))  # quite a lot!

# Create wide dataset
dfw_sulf <- reshape(dfl_sulf,idvar= c("pid", "arm", "gender", "site", "age", "weight", "parbase"), direction="wide",timevar="day")
dfw_pyr <- reshape(dfl_pyr,idvar= c("pid", "arm", "gender", "site", "age", "weight", "parbase"), direction="wide",timevar="day")
str(dfw_sulf)
str(dfw_pyr)


# create an joint dataset
un <- union(samp_pids1,samp_pids2)
mydatl <- mydatl[mydatl$pid %in% un, ]
str(mydatl)
mydatw <- reshape(mydatl,idvar= c("pid", "arm", "gender", "site", "age", "weight", "parbase"), direction="wide",timevar="day")
str(mydatw)


#----------------------------------
# LOG RESPONSES
#----------------------------------

dfl_pyr$Pyrconcentration <- log10(dfl_pyr$Pyrconcentration + 1) 
dfl_sulf$Sulfconcentration <- log10(dfl_sulf$Sulfconcentration + 1)
mydatl$Pyrconcentration <- log10(mydatl$Pyrconcentration + 1) 
mydatl$Sulfconcentration <- log10(mydatl$Sulfconcentration + 1) 
str(dfl_pyr)
str(dfl_sulf)
str(mydatl)
# update wide data sets
dfw_pyr <- reshape(dfl_pyr,idvar= c("pid", "arm", "gender", "site", "age", "weight", "parbase"), direction="wide",timevar="day")
str(dfw_pyr)
dfw_sulf <- reshape(dfl_sulf,idvar= c("pid", "arm", "gender", "site", "age", "weight", "parbase"), direction="wide",timevar="day")
str(dfw_sulf)
mydatw <- reshape(mydatl,idvar= c("pid", "arm", "gender", "site", "age", "weight", "parbase"), direction="wide",timevar="day")
str(mydatw)



#-----------------------------------------------------------------------
# STEP 1: Fit a single regression model, ignoring grouping by subject
#-----------------------------------------------------------------------
dfl_pyr <- na.omit(dfl_pyr)
fit_gen <- nls(Pyrconcentration~SSbiexp(day,A1,lrc1,A2,lrc2), data=dfl_pyr)
summary(fit_gen)
plot(fit_gen, pid~resid(.),abline = 0)

#-----------------------------------------------------------------------
# STEP 2: Fit individual models for each subject
#-----------------------------------------------------------------------
fit_each <- nlsList(Pyrconcentration~SSbiexp(day,A1,lrc1,A2,lrc2)|pid, data=dfl_pyr)
summary(fit_each)
plot(intervals(fit_each))

#-----------------------------------------------------------------------
# STEP 3: NLME with block diagonal
#-----------------------------------------------------------------------
start <- fixef(fit_each)
start
dfl_pyr1 <- groupedData(Pyrconcentration~day|pid, data = dfl_pyr)

fit1 <- nlme(Pyrconcentration~SSbiexp(day,A1,lrc1,A2,lrc2), data=dfl_pyr1,
             fixed = A1 + lrc1 + A2 + lrc2 ~ 1,
             start = c(start[1], start[2], start[3], start[4]),
             random = pdDiag(A1 + lrc1 + A2 + lrc2 ~ 1))
summary(fit1)
plot(fit1, adj=-1, pch=20, col="black")
intervals(fit1)


#-----------------------------------------------------------------------
# STEP 4: NLME with different random effect structures
#-----------------------------------------------------------------------
fit2 <- update(fit1, random=pdDiag(A1 + A2 + lrc2 ~ 1))
summary(fit2)
plot(fit2, adj=-1, pch=20, col="black")
intervals(fit2)
anova(fit2,fit1)

fit3 <- update(fit2, random=A1 + A2 + lrc2 ~ 1 ,
               control=nlmeControl(msMaxIter=100))
summary(fit3)
plot(fit3, adj=-1, pch=20, col="black")
anova(fit2,fit3)




#-----------------------------------------------------------------------
# STEP 5: NLME with different covariates
#-----------------------------------------------------------------------
re1 <- ranef(fit3, augFrame = T)

plot(re1, form=A1~arm, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4)) 
plot(re1, form=A2~arm, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4))
plot(re1, form=lrc2~arm, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4))

plot(re1, form=A1~gender, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4))
plot(re1, form=A2~gender, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4))
plot(re1, form=lrc2~gender, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4))

plot(re1, form=A1~site, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4))
plot(re1, form=A2~site, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4))
plot(re1, form=lrc2~site, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4))

plot(re1, form=A1~weight, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=A2~weight, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=lrc2~weight, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4), xlab = "", pch=20, col="black")

plot(re1, form=A1~age, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=A2~age, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=lrc2~age, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4), xlab = "", pch=20, col="black")

plot(re1, form=A1~parbase, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=A2~parbase, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=lrc2~parbase, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4), xlab = "", pch=20, col="black")


fit4 <- update(fit3, fixed = c(A1~gender+site, lrc1~1,  A2 ~ gender+site+weight+age, lrc2~arm+gender+site+weight+age),
               start=c(start[1],0,0,start[2],start[3],0,0,0,0,start[4],0,0,0,0,0))
summary(fit4)


fit5 <- update(fit3, fixed = c(A1~gender+site, lrc1~1,  A2 ~ site+weight+age, lrc2~arm+gender+site+weight+age),
               start=c(start[1],0,0,start[2],start[3],0,0,0,start[4],0,0,0,0,0))
summary(fit5)
anova(fit5,fit4)  # choose 5


fit6 <- update(fit3, fixed = c(A1~gender+site, lrc1~1,  A2 ~ weight+age, lrc2~arm+gender+site+weight+age),
               start=c(start[1],0,0,start[2],start[3],0,0,start[4],0,0,0,0,0))
summary(fit6)
anova(fit6,fit5)  # choose 6


fit7 <- update(fit3, fixed = c(A1~gender+site, lrc1~1,  A2 ~ weight, lrc2~arm+gender+site+weight+age),
               start=c(start[1],0,0,start[2],start[3],0,start[4],0,0,0,0,0))
summary(fit7)
anova(fit7,fit6)  # choose 7


fit8 <- update(fit3, fixed = c(A1~site, lrc1~1,  A2 ~ weight, lrc2~arm+gender+site+weight+age),
               start=c(start[1],0,start[2],start[3],0,start[4],0,0,0,0,0))
summary(fit8)
anova(fit8,fit7)  # MAYBE 8?


fit9 <- update(fit3, fixed = c(A1~site, lrc1~1,  A2 ~ weight, lrc2~arm+site+weight+age),
               start=c(start[1],0,start[2],start[3],0,start[4],0,0,0,0))
summary(fit9)
anova(fit9,fit8)
anova(fit9,fit7)


# CHOOSE MODEL 7

anova(fit3, fit7)
plot(fit7, adj=-1, pch=20, col="black")











#-----------------------------------------------------------------------
# STEP 5: variance function
#-----------------------------------------------------------------------

# residuals versus fitted values
df10 <- data.frame(residuals = resid(fit7), 
                   fitted = predict(fit7))
ggplot(data = df10,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)

plot(fit7, resid(.,type="pearson")~fitted(.)|arm, adj=-0.3, 
     pch = 20, col= "black",
     par=list(cex=2), ylab = list("Pearson residuals", cex=2), xlab = list("Fitted values", cex=2))

plot(fit7, resid(.,type="pearson")~fitted(.)|gender, adj=-0.3, 
     pch = 20, col= "black",
     par=list(cex=2), ylab = list("Pearson residuals", cex=2), xlab = list("Fitted values", cex=2))

plot(fit7, resid(.,type="pearson")~fitted(.)|site, adj=-0.3, 
     pch = 20, col= "black",
     par=list(cex=2), ylab = list("Pearson residuals", cex=2), xlab = list("Fitted values", cex=2))



# Weight
dfb <- data.frame(residuals = resid(fit7,type="pearson"), 
                  fitted = dfl_pyr$weight)
ggplot(data = dfb,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Weight")+
  theme_minimal(base_size = 22)

# Age
dfb <- data.frame(residuals = resid(fit7,type="pearson"), 
                  fitted = dfl_pyr$age)
ggplot(data = dfb,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Age")+
  theme_minimal(base_size = 22)




#fit20 <- update(fit7, weights=varExp(form= ~age))
#fit20 <- update(fit7, weights=varExp())
fit20 <- update(fit7, weights=varPower(form= ~age))
anova(fit20,fit7)
fit20 <- update(fit7, weights=varPower(fixed = ~weight))
anova(fit20,fit7)
fit20 <- update(fit7, weights=varFixed(~age))
fit20 <- update(fit7, weights=varFixed(~weight))
fit20 <- update(fit7, weights=varFixed(~age+weight))
summary(fit20)

plot(fit7, resid(.)~day, abline=0)
plot(Variogram(fit7), form=~day)

fit30 <- update(fit7, corr=corCAR1(form=~day), control=nlmeControl(msMaxIter = 200))
anova(fit7,fit30)

#-----------------------------------------------------------------------
# STEP 6: CORRELATION STRUCTURE
#-----------------------------------------------------------------------

ACF(fit15)
plot(ACF(fit15), alpha=0.05, ylab = list("Autocorrelation", cex=2), xlab = list("Lag", cex=2))
plot(ACF(fit15, resType = "n"), alpha=0.05)

fit30 <- update(fit7, correlation=corCompSymm())
fit31 <- update(fit7, correlation=corSymm())
fit32 <- update(fit7, correlation=corAR1())
anova(fit7, fit30)


#-----------------------------------------------------------------------
# STEP 7: Residual Checks
#-----------------------------------------------------------------------


qqnorm(fit7, col="black", pch=20, ylab = list("Quantiles of standard normal", cex=2), xlab = list("Standardised residuals", cex=2))
qqnorm(fit7, ~resid(.)|arm, col="black", pch=20, ylab = list("Quantiles of standard normal", cex=2), xlab = list("Standardised residuals", cex=2))
qqnorm(fit7, ~resid(.)|gender, col="black", pch=20, ylab = list("Quantiles of standard normal", cex=2), xlab = list("Standardised residuals", cex=2))
qqnorm(fit7, ~resid(.)|site, col="black", pch=20, ylab = list("Quantiles of standard normal", cex=2), xlab = list("Standardised residuals", cex=2))


qqnorm(fit7, ~ranef(.), col="black", pch=20, ylab = list("Quantiles of standard normal", cex=2), xlab = list("Standardised residuals", cex=2), par=list(cex=2))


df <- data.frame(sul = dfl_pyr$Pyrconcentration, 
                 fitted = fitted(fit7))
ggplot(data = df,
       aes(x = fitted, y = sul)) +
  geom_point(size=2) + 
  ylab("Observed Pyrimethamine values")+
  xlab("Fitted Pyrimethamine values")+
  theme_minimal(base_size = 22)+
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1.5)










qqnorm(fit7, col="black", pch=20)

qqnorm(fit7, ~ranef(.), col="black", pch=20)

pairs(fit7, ~ranef(.)|arm, col="black", pch=20)
pairs(fit7, ~ranef(.)|site, col="black", pch=20)
pairs(fit7, ~ranef(.)|gender, col="black", pch=20)

plot(fit7, Pyrconcentration~fitted(.), id=0.005)

df <- data.frame(pyr = dfl_pyr$Pyrconcentration, 
                 fitted = fitted(fit7))
ggplot(data = df,
       aes(x = fitted, y = pyr)) +
  geom_point(size=2) + 
  ylab("Observed Sulfadoxine values")+
  xlab("Fitted Sulfadoxine values")+
  theme_minimal(base_size = 22)+
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1.5)


fit <- update(fit3, fixed = c(A1~gender+site, lrc1~1,  A2 ~ weight, lrc2~arm+gender+site+weight+age),
                             start=c(start[1],0,0,start[2],start[3],0,start[4],0,0,0,0,0), method="REML")
summary(fit)

