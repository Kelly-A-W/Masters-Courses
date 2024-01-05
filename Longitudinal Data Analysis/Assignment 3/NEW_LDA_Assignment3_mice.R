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
library(VIM)
library(mice)
library(miceadds)
library(broom.mixed)

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

ggplot(data = new_dfl,
       aes(x = 1:nrow(new_dfl), y = lgamedens)) +
  geom_point(size=2) + 
  ylab("Logged gametocyte density")+
  xlab("Index")+
  theme_minimal(base_size = 22)
ggplot(new_dfl,aes(x=lgamedens)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  xlab("Gametocyte Density")+
  theme_minimal(base_size = 27)



#------------------------------------
# Choose Predictorss 
#------------------------------------
dfl <- new_dfl[,-12]
dfl$pid <- as.integer(as.factor(dfl$pid))
dfl$pday_squared <- (dfl$pday)^2     # need to include this in predictor matrix

# predictor matrix
imp0 <- mice(dfl, maxit = 0)
pred <- imp0$predictorMatrix

# set variables that don't need imputation
pred["site",] <- 0
pred["pid",] <- 0
pred["pday",] <- 0
pred["gender",] <- 0
pred["weight",] <- 0
pred["age",] <- 0
pred["pday_squared",] <- 0

# set variables that do need imputation
pred["lpardens", ] <- c(1,1,1,0,0,0,0,0,1,1,1,1)
pred["lgamedens", ] <- c(1,1,1,1,0,1,0,0,1,1,1,1)
pred["Hb", ] <- c(1,1,1,1,0,0,1,1,1,1,1,1)
pred["lpyrconcentration", ] <- c(1,1,1,1,0,0,0,0,1,1,1,1)
pred["lsulfconcentration", ] <- c(1,1,1,1,0,0,1,0,1,1,1,1)

# set longitudinal data specific preds
pred[,"pid"] <- -2
pred[,"pday"] <- 2


# set imputation method for repeated measures
meth <- imp0$method
meth["lpardens"] <- "2l.norm"
meth["lgamedens"] <- "2l.norm"
meth["Hb"] <- "2l.norm"
meth["lpyrconcentration"] <- "2l.norm"
meth["lsulfconcentration"] <- "2l.norm"

# set passive imputation method to preserve pday^2 transformation
meth["pday_squared"] <- "~(pday)^2"    

# ensure imputations are all non-negative
post<-imp0$post
post[4:8] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 20))"

imp1 <- mice(dfl, seed=23109,
             pred = pred,
             meth = meth,
             maxit = 5,
             post=post)   




#----------------------------------------
# STEP 2: Assess Imputations
#----------------------------------------

# scatter-plots of the observed and the imputed values for each variable
# distribution of the imputed values should not be very different from that of the observed values
stripplot(imp1,lpardens + lgamedens + Hb + lpyrconcentration + lsulfconcentration ~ .imp, pch=20, 
          cex=2, cex.lab=2, cex.numbers = 1.5, cex.axis=1.5)  


# CHECKING CONVERGENCE OF IMPUTATIONS
plot(imp1, c("lpardens", "lgamedens", "Hb", "lpyrconcentration", "lsulfconcentration")) # good mixing of patterns indicates good convergence
densityplot(imp1)  # want densities to overlap and be similar to density of observed values
# hgt and tv imputations don't look great 


# scatter-plots of bmi vs chl for observed data and for the M imputed datasets
xyplot(imp1,  lgamedens~lpardens|.imp, pch=20, cex=2, cex.lab=2, cex.numbers = 1.5, cex.axis=1.5,
       ylab= "Logged gametocyte density", xlab="Logged parasite density")
xyplot(imp1,  lgamedens~Hb|.imp, pch=20, cex=2, cex.lab=2, cex.numbers = 1.5, cex.axis=1.5,
       ylab= "Logged gametocyte density", xlab="Haemoglobine")
# we want the all the scatter-plots to look alike
# dont want our imputations to change the relationship between our variables


#-----------------------------------------------------------------
# STEP 3: Extract Parameter Estimates (for each imputed dataset)
#-----------------------------------------------------------------


mod1 <- with(imp1, lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + gender, 
                       random = list(~1|pid), method="REML", 
                       weights = varPower(fixed=0.5), 
                       correlation = corCompSymm(), 
                       control=lmeControl(tolerance = 0.01)))
mod1$analyses



pool(mod1)

summary(pool(mod1))

sapply(mod1$analyses, AIC) 
sapply(mod1$analyses, BIC) 





#--------------------------------------------------------
# STEP 6: Sensitivity Analysis
#--------------------------------------------------------
# investigate the sensitivity of the results to different values of a variable, in this case chl
# Suppose we suspect that our imputations under MAR assumption are too low
# We can check this by multiplying the imputations by a factor
k    <- seq(1, 1.5, 0.1)
est  <- vector("list", length(k))
for(i in 1:length(k)) {
  post["lgamedens"] <- paste("imp[[j]][,i] <- ", k[i], "* imp[[j]][,i]")  
  imp <- mice(dfl, post=post, seed=23109, print=T, maxit=20, meth=meth, pred=pred)
  fit <- with(imp, lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + gender, 
                        random = list(~1|pid), method="REML", 
                        weights = varPower(fixed=0.5), 
                        correlation = corCompSymm(), 
                        control=lmeControl(tolerance = 0.01)))
  est[[i]] <- summary(pool(fit))
}
est    # scenarios involving increasing imputations for chl by 0% (MAR), 10%, 20%, 30%, 40% and 50% (MNAR)
# here we see that the intercept changes but the other effects less so
# and we dont mind that the intercept changes because we are more interested in the other estimates being consistent and accurate
# because we are more interested in the relationships between variables
est2 <- est


k    <- seq(1, 1.5, 0.1)
est  <- vector("list", length(k))
for(i in 1:length(k)) {
  post["lgamedens"] <- paste("imp[[j]][,i] <- ", k[i], "* imp[[j]][,i]") 
  post["lpardens"] <- paste("imp[[j]][,i] <- ", k[i], "* imp[[j]][,i]") 
  post["Hb"] <- paste("imp[[j]][,i] <- ", k[i], "* imp[[j]][,i]") 
  post["lpyrconcentration"] <- paste("imp[[j]][,i] <- ", k[i], "* imp[[j]][,i]") 
  post["lsulfconcentration"] <- paste("imp[[j]][,i] <- ", k[i], "* imp[[j]][,i]") 
  imp <- mice(dfl, post=post, seed=23109, print=T, maxit=20, meth=meth, pred=pred)
  fit <- with(imp, lme(lgamedens~pday + I(pday^2) + site  + lpardens + Hb + gender, 
                       random = list(~1|pid), method="REML", 
                       weights = varPower(fixed=0.5), 
                       correlation = corCompSymm(), 
                       control=lmeControl(tolerance = 0.01)))
  est[[i]] <- summary(pool(fit))
}
est   















