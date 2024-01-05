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
dfl_sulf <- na.omit(dfl_sulf)
fit1 <- nls(Sulfconcentration~SSbiexp(day,A1,lrc1,A2,lrc2), data=dfl_sulf)
summary(fit1)
plot(fit1, pid~resid(.),abline = 0)

#-----------------------------------------------------------------------
# STEP 2: Fit individual models for each subject
#-----------------------------------------------------------------------
fit2 <- nlsList(Sulfconcentration~SSbiexp(day,A1,lrc1,A2,lrc2)|pid, data=dfl_sulf)
summary(fit2)
plot(intervals(fit2))

#-----------------------------------------------------------------------
# STEP 3: NLME with block diagonal
#-----------------------------------------------------------------------
start <- fixef(fit2)
start
fit3 <- nlme(fit2, random = pdDiag(A1+lrc1+A2+lrc2~1)) # random effects on each parameter, but still no covariates !
dfl_sulf1 <- groupedData(Sulfconcentration~day|pid, data = dfl_sulf)
fit3 <- nlme(Sulfconcentration~SSbiexp(day,A1,lrc1,A2,lrc2), data=dfl_sulf1,
             fixed = A1 + lrc1 + A2 + lrc2 ~ 1,
             start = c(start[1], start[2], start[3], start[4]),
             random = pdDiag(A1 + lrc1 + A2 + lrc2 ~ 1))
summary(fit3)
plot(fit3, adj=-1, pch=20, col="black")
intervals(fit3)


#-----------------------------------------------------------------------
# STEP 4: NLME with different random effect structures
#-----------------------------------------------------------------------
fit4 <- update(fit3, random=pdDiag(A1 + lrc1 + lrc2 ~ 1), 
               control=nlmeControl(pnlsMaxIter=20))
summary(fit4)
plot(fit4, adj=-1, pch=20, col="black")
intervals(fit4)
anova(fit4,fit3)

fit5 <- update(fit3, random=A1 + lrc1 + A2 + lrc2 ~ 1 ,
               control=nlmeControl(msMaxIter=100))

fit6 <- update(fit3, random=pdBlocked(list(A1 + lrc1 ~ 1, A2 + lrc2 ~ 1)))
summary(fit6)
plot(fit6, adj=-1, pch=20, col="black")
intervals(fit6)
anova(fit3,fit6)



#-----------------------------------------------------------------------
# STEP 5: NLME with different covariates
#-----------------------------------------------------------------------
re1 <- ranef(fit3, augFrame = T)

plot(re1, form=A1~arm, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4)) 
plot(re1, form=lrc1~arm, par=list(cex=4), ylab = list(TeX(r'($\b_{2i}$)'), cex=4))
plot(re1, form=A2~arm, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4))
plot(re1, form=lrc2~arm, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4))

plot(re1, form=A1~gender, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4))
plot(re1, form=lrc1~gender, par=list(cex=4), ylab = list(TeX(r'($\b_{2i}$)'), cex=4))
plot(re1, form=A2~gender, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4))
plot(re1, form=lrc2~gender, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4))

plot(re1, form=A1~site, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4))
plot(re1, form=lrc1~site, par=list(cex=4), ylab = list(TeX(r'($\b_{2i}$)'), cex=4))
plot(re1, form=A2~site, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4))
plot(re1, form=lrc2~site, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4))

plot(re1, form=A1~weight, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=lrc1~weight, par=list(cex=4), ylab = list(TeX(r'($\b_{2i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=A2~weight, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=lrc2~weight, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4), xlab = "", pch=20, col="black")

plot(re1, form=A1~age, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=lrc1~age, par=list(cex=4), ylab = list(TeX(r'($\b_{2i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=A2~age, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=lrc2~age, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4), xlab = "", pch=20, col="black")

plot(re1, form=A1~parbase, par=list(cex=4), ylab = list(TeX(r'($\b_{1i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=lrc1~parbase, par=list(cex=4), ylab = list(TeX(r'($\b_{2i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=A2~parbase, par=list(cex=4), ylab = list(TeX(r'($\b_{3i}$)'), cex=4), xlab = "", pch=20, col="black")
plot(re1, form=lrc2~parbase, par=list(cex=4), ylab = list(TeX(r'($\b_{4i}$)'), cex=4), xlab = "", pch=20, col="black")


fit7 <- update(fit3, fixed = c(A1~weight, lrc1~site + weight, A2 ~ arm +gender+site+weight+age, lrc2~arm+gender+site+weight+parbase),
               start=c(start[1],0,start[2],0,0,start[3],0,0,0,0,0,start[4],0,0,0,0,0))
summary(fit7)


fit8 <- update(fit3, fixed = c(A1~weight, lrc1~site + weight, A2 ~ arm +gender+site+weight+age, lrc2~arm+site+weight+parbase),
               start=c(start[1],0,start[2],0,0,start[3],0,0,0,0,0,start[4],0,0,0,0))
summary(fit8)
anova(fit8, fit7)


fit9 <- update(fit3, fixed = c(A1~weight, lrc1~site + weight, A2 ~ gender+site+weight+age, lrc2~arm+site+weight+parbase),
               start=c(start[1],0,start[2],0,0,start[3],0,0,0,0,start[4],0,0,0,0))
summary(fit9)
anova(fit9, fit8) #choose 9


fit10 <- update(fit3, fixed = c(A1~weight, lrc1~site + weight, A2 ~ gender+site+weight+age, lrc2~arm+weight+parbase),
               start=c(start[1],0,start[2],0,0,start[3],0,0,0,0,start[4],0,0,0))
summary(fit10)
anova(fit10, fit9) # choose 10


fit11 <- update(fit3, fixed = c(A1~weight, lrc1~weight, A2 ~ gender+site+weight+age, lrc2~arm+weight+parbase),
                start=c(start[1],0,start[2],0,start[3],0,0,0,0,start[4],0,0,0))
summary(fit11)
anova(fit11, fit10) # choose 11


fit12 <- update(fit3, fixed = c(A1~weight, lrc1~weight, A2 ~ gender+site+weight+age, lrc2~arm+weight),
                start=c(start[1],0,start[2],0,start[3],0,0,0,0,start[4],0,0))
# not converging for LME step: nlminb(), so increase iterations 
fit12 <- update(fit3, fixed = c(A1~weight, lrc1~weight, A2 ~ gender+site+weight+age, lrc2~arm+weight),
                start=c(start[1],0,start[2],0,start[3],0,0,0,0,start[4],0,0),
                control=nlmeControl(msMaxIter = 100, tolerance = 1e-4))
summary(fit12)
anova(fit12, fit11) #choose 12


fit13 <- update(fit3, fixed = c(A1~weight, lrc1~1, A2 ~ gender+site+weight+age, lrc2~arm+weight),
                start=c(start[1],0,start[2],start[3],0,0,0,0,start[4],0,0))
summary(fit13)
anova(fit13, fit12)  # choose 13


fit14 <- update(fit3, fixed = c(A1~weight, lrc1~1, A2 ~ site+weight+age, lrc2~arm+weight),
                start=c(start[1],0,start[2],start[3],0,0,0,start[4],0,0))
summary(fit14)
anova(fit14, fit13) #choose 14


fit15 <- update(fit3, fixed = c(A1~weight, lrc1~1, A2 ~ site+weight, lrc2~arm+weight),
                start=c(start[1],0,start[2],start[3],0,0,start[4],0,0))
summary(fit15)
anova(fit15, fit14) # fit 15
anova(fit3, fit15)  #fit15

fit16 <- update(fit3, fixed = c(A1~1, lrc1~1, A2 ~ site+weight, lrc2~arm+weight),
                start=c(start[1],start[2],start[3],0,0,start[4],0,0))
anova(fit16, fit15) # choose 15
fit17 <- update(fit3, fixed = c(A1~weight, lrc1~1, A2 ~ site, lrc2~arm+weight),
                start=c(start[1],0,start[2],start[3],0,start[4],0,0))
anova(fit17, fit15) # choose 15
fit18 <- update(fit3, fixed = c(A1~weight, lrc1~1, A2 ~ site+weight, lrc2~arm),
                start=c(start[1],0,start[2],start[3],0,0,start[4],0))
anova(fit18, fit15) # choose 15
fit18 <- update(fit3, fixed = c(A1~weight, lrc1~1, A2 ~ weight, lrc2~arm+weight),
                start=c(start[1],0,start[2],start[3],0,start[4],0,0))
anova(fit18, fit15) # choose 15

# CHOOSE MODEL 15 !!!!!!!!!!!!


#-----------------------------------------------------------------------
# STEP 5: variance function
#-----------------------------------------------------------------------

# residuals versus fitted values
df10 <- data.frame(residuals = resid(fit15), 
                   fitted = predict(fit15))
ggplot(data = df10,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)

plot(fit15, resid(.,type="pearson")~fitted(.)|arm, adj=-0.3, 
      pch = 20, col= "black",
     par=list(cex=2), ylab = list("Pearson residuals", cex=2), xlab = list("Fitted values", cex=2))

plot(fit15, resid(.,type="pearson")~fitted(.)|site, adj=-0.3, 
     pch = 20, col= "black",
     par=list(cex=2), ylab = list("Pearson residuals", cex=2), xlab = list("Fitted values", cex=2))


# Weight
dfb <- data.frame(residuals = resid(fit15,type="pearson"), 
                  fitted = dfl_sulf$weight)
ggplot(data = dfb,
       aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Pearson residuals")+
  xlab("Weight")+
  theme_minimal(base_size = 22)




fit20 <- update(fit15, weights=varConstPower(power = 0.1))
fit20 <- update(fit15, weights=varExp())
fit20 <- update(fit15, weights=varPower(form= ~weight))
fit20 <- update(fit15, weights=varFixed(~-weight))
summary(fit20)


#-----------------------------------------------------------------------
# STEP 6: CORRELATION STRUCTURE
#-----------------------------------------------------------------------

ACF(fit15)
plot(ACF(fit15), alpha=0.05, ylab = list("Autocorrelation", cex=2), xlab = list("Lag", cex=2))
plot(ACF(fit15, resType = "n"), alpha=0.05)

fit30 <- update(fit15, correlation=corCompSymm())
fit31 <- update(fit15, correlation=corSymm())
fit32 <- update(fit15, correlation=corAR1())
anova(fit15, fit30)


#-----------------------------------------------------------------------
# STEP 7: Residual Checks
#-----------------------------------------------------------------------


qqnorm(fit15, col="black", pch=20, ylab = list("Quantiles of standard normal", cex=2), xlab = list("Standardised residuals", cex=2))
qqnorm(fit15, ~resid(.)|arm, col="black", pch=20, ylab = list("Quantiles of standard normal", cex=2), xlab = list("Standardised residuals", cex=2))
qqnorm(fit15, ~resid(.)|site, col="black", pch=20, ylab = list("Quantiles of standard normal", cex=2), xlab = list("Standardised residuals", cex=2))


qqnorm(fit15, ~ranef(.), col="black", pch=20, ylab = list("Quantiles of standard normal", cex=2), xlab = list("Standardised residuals", cex=2), par=list(cex=2))




df <- data.frame(sul = dfl_sulf$Sulfconcentration, 
                 fitted = fitted(fit15))
ggplot(data = df,
       aes(x = fitted, y = sul)) +
  geom_point(size=2) + 
  ylab("Observed Sulfadoxine values")+
  xlab("Fitted Sulfadoxine values")+
  theme_minimal(base_size = 22)+
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1.5)





fit <- update(fit3, fixed = c(A1~weight, lrc1~1, A2 ~ site+weight, lrc2~arm+weight),
                start=c(start[1],0,start[2],start[3],0,0,start[4],0,0), method ="REML")
summary(fit)







