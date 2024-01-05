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
library(Amelia)
library(mitools)

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


a.out <- amelia(dfl, m=5, ts="pday", cs="pid", noms = c("site", "gender"))
a.out<-transform(a.out, pday2 = pday^2)  

plot(a.out, which.vars=6:9)

mi_list <- imputationList((a.out$imputations))
#mod1 <- with(mi_list, lmer(lgamedens~pday + I(pday^2) + site  + (1|pid), REML = F))
mod1 <- with(mi_list, lme(lgamedens~pday + I(pday^2) + site, random = list(~pday|pid), method="ML"))

pool(mod1)
MIcombine(mod1)
summary(MIcombine(mod1))

mod2 <- with(mi_list, lme(lgamedens~pday + I(pday^2) + site, random = list(~1|pid), method="ML"))






library(merTools)
fastdisp(mod1)
modelFixedEff(mod1)



