# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(ggcorrplot)
library(RColorBrewer)
library(readxl)
library(VIM)

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


#--------------------------------------------------
# EDA: 
#--------------------------------------------------

#--------------------------------------------------
# UNIVARIATE EDA
#--------------------------------------------------


# CONTINOUS VARIABLES:

#-------------------------------------------------- Box plots for gamedens
# we are expecting that a log transformation will probably be necessary
ggplot(mydatl,aes(x=gamedens)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  xlab("Gametocyte Density")+
  theme_minimal(base_size = 27)
# EXTREMELY skew data !!!!

# need log transformation
mydatl$gamedens<- log10(1+mydatl$gamedens)
colnames(mydatl)[5] <- "lgamedens"
head(mydatl)

# update wide data sets
mydatw <- reshape(mydatl,idvar= c("pid", "gender", "site", "age", "weight", "PIoutcome"), direction="wide",timevar="pday")
str(mydatw)

ggplot(mydatl,aes(x=lgamedens)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  xlab("Logged Gametocyte density")+
  theme_minimal(base_size = 27)
# still very skew !!!!!

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
  xlab("Logged Gametocyte density")+
  theme_minimal(base_size = 27)


mydatl <- new_dfl
# Wide format:
mydatw <- reshape(mydatl,idvar= c("pid", "gender", "site", "age", "weight", "PIoutcome"), direction="wide",timevar="pday")
str(mydatw)


#-------------------------------------------------------------





#--------------------------------------------------------------
# BIVARIATE EDA
#--------------------------------------------------------------

# FACTOR VARIABLES 

# need day as a factor variable to plot
dfl <- mydatl
dfl$pday <- as.factor(dfl$pday)  # 6 time points
levels(dfl$pday) <- c(0,3,7,14,28,42) # prefer to use actual days for time



#-------------------------------------------------- All Variables
df <- data.frame(Variable = c(rep("Parasite", length(mydatl$lpardens)), rep("Gametocyte", length(mydatl$lgamedens)), 
                              rep("Haemoglobin", length(mydatl$Hb)), rep("Pyrmethamine", length(mydatl$lpyrconcentration)), 
                              rep("Sulfadoxine", length(mydatl$lsulfconcentration))), 
                 Amount = c(mydatl$lpardens, mydatl$lgamedens, mydatl$Hb, mydatl$lpyrconcentration, mydatl$lsulfconcentration),
                 Day = c(dfl$pday, dfl$pday, dfl$pday, dfl$pday, dfl$pday))
df$Variable <- as.factor(df$Variable)
Variable <- df$Variable
ggplot(df,aes(x=Day,y=Amount,fill=Variable)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  facet_wrap(~Variable) +
  ylab("Concentration or Density")+
  xlab("Day")+
  theme_minimal(base_size = 27)
ggplot() +
  aes(x = df$Day, color = Variable, group = Variable, y = df$Amount) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  ylab("Mean Concentration or Density")+
  xlab("Day") +
  theme_minimal(base_size = 27) 


#-------------------------------------------------- lgamedens - site
# Box-plots
ggplot(dfl,aes(x=pday,y=lgamedens,fill=site)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  facet_wrap(~site) +
  ylab("Logged gametocyte density")+
  xlab("Day")+
  theme_minimal(base_size = 27) 
# Interaction Plot 
Site <- dfl$site
ggplot() +
  aes(x = dfl$pday, color = Site, group = Site, y = dfl$lgamedens) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  ylab("Mean logged gametocyte density")+
  xlab("Day") +
  theme_minimal(base_size = 27)


#-------------------------------------------------- lgamedens - gender
# box-plots
ggplot(dfl,aes(x=pday,y=lgamedens,fill=gender)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  facet_wrap(~gender) +
  ylab("Logged gametocyte density")+
  xlab("Day")+
  theme_minimal(base_size = 27) 
# Interaction Plot 
Sex <- dfl$gender
ggplot() +
  aes(x = dfl$pday, color = Sex, group = Sex, y = dfl$lgamedens) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  ylab("Mean logged gametocyte density")+
  xlab("Day") +
  theme_minimal(base_size = 27)

#-------------------------------------------------- CONTINOUS VARIABLES 

corr <- cor(mydatl[,-c(1,2,3,9,12)],  use = "pairwise.complete.obs")
ggcorrplot(corr, 
           hc.order = F, 
           type = "lower",
           show.legend = F,
           lab = T, 
           tl.cex = 20, 
           lab_size = 5) 




#--------------------------------------------------------------
# IMPUTATIONS EDA
#--------------------------------------------------------------

#-------------------------------------------------- EVALUATING MISSINGNESS

df <- mydatl[,-c(2,3,12)]
colnames(df) <- c("site", "lpar", "lgam", "Hb", "lpyr", "lsulf", "sex", "weight","age")
md.pattern(df)  

md.pairs(mydatl[,-c(2,3,12)])    


marginplot(mydatl[,c("lpardens", "lgamedens")], 
           col = mdc(1:2), cex=2, cex.lab=2, cex.numbers = 1.5, pch=19, cex.axis=1.5,
           xlab="Logged parasite density", ylab = "Logged gametocyte density")
marginplot(mydatl[,c("Hb", "lgamedens")],  
           col = mdc(1:2), cex=2, cex.lab=2, cex.numbers = 1.5, pch=19, cex.axis=1.5,
           xlab="Haemoglobin", ylab = "Logged gametocyte density")
marginplot(mydatl[,c("lpyrconcentration", "lgamedens")],  
           col = mdc(1:2), cex=2, cex.lab=2, cex.numbers = 1.5, pch=19, cex.axis=1.5,
           xlab="Logged Pyrmethamine concentration", ylab = "Logged gametocyte density")
marginplot(mydatl[,c("lsulfconcentration", "lgamedens")],  
           col = mdc(1:2), cex=2, cex.lab=2, cex.numbers = 1.5, pch=19, cex.axis=1.5,
           xlab="Logged Sulfadoxine concentration", ylab = "Logged gametocyte density")
marginplot(mydatl[,c("Hb", "lpardens")],  
           col = mdc(1:2), cex=2, cex.lab=2, cex.numbers = 1.5, pch=19, cex.axis=1.5,
           xlab="Haemoglobin", ylab = "Logged parasite density")
marginplot(mydatl[,c("lpyrconcentration", "lpardens")],  
           col = mdc(1:2), cex=2, cex.lab=2, cex.numbers = 1.5, pch=19, cex.axis=1.5,
           xlab="Logged Pyrmethamine concentration", ylab = "Logged parasite density")
marginplot(mydatl[,c("lsulfconcentration", "lpardens")],  
           col = mdc(1:2), cex=2, cex.lab=2, cex.numbers = 1.5, pch=19, cex.axis=1.5,
           xlab="Logged Sulfadoxine concentration", ylab = "Logged parasite density")
marginplot(mydatl[,c("lpyrconcentration", "Hb")],  
           col = mdc(1:2), cex=2, cex.lab=2, cex.numbers = 1.5, pch=19, cex.axis=1.5,
           xlab="Logged Pyrmethamine concentration", ylab = "Haemoglobin")
marginplot(mydatl[,c("lsulfconcentration", "Hb")],  
           col = mdc(1:2), cex=2, cex.lab=2, cex.numbers = 1.5, pch=19, cex.axis=1.5,
           xlab="Logged Sulfadoxine concentration", ylab = "Haemoglobin")
marginplot(mydatl[,c("lsulfconcentration", "lpyrconcentration")],  
           col = mdc(1:2), cex=2, cex.lab=2, cex.numbers = 1.5, pch=19, cex.axis=1.5,
           xlab="Logged Sulfadoxine concentration", ylab = "Logged Pyrmethamine concentration")



#-------------------------------------------------- CHOOSING PREDICTORS

# Pearson correlations:
round(cor(mydatl[,-c(1,2,3,9,12)], use="pair"), 3) # highly correlated variables will make better predictors

# Correlations between variables and missing data indicators:
round(cor(y=mydatl[,-c(1,2,3,9,12)], x=!is.na(mydatl[,-c(1,2,3,9,12)]), use="pair"), 3)  # high correlation make beter predictors

# Proportion of usable cases:
p <- md.pairs(mydatl[,-c(1,2,3,9,10,11,12)])
round(p$mr/(p$mr+p$mm), 3) # measures how many cases on the target variable have observations on the predictor
# high makes good predictors 


