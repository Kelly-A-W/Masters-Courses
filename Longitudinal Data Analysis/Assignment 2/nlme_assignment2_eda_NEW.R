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

#--------------------------------------------------
# UNIVARIATE EDA
#--------------------------------------------------

# FACTOR VARIABLES


#-------------------------------------------------- SITE:
# use joint dataset !
nrow(mydatw[mydatw$site=="South",])
nrow(mydatw[mydatw$site=="North",])
tab1 <- data.frame(
  Group=levels(mydatl$site),
  value=c(nrow(mydatw[mydatw$site=="South",]),nrow(mydatw[mydatw$site=="North",])))
tab1$value <- tab1$value/sum(tab1$value)*100
tab1
labels <- c(round(tab1$value[1]), round(tab1$value[2]))
labels <- c(paste(labels[1],"%"), paste(labels[2],"%"))
ggplot(tab1, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.5), size=20) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 45, face = "bold"), legend.text=element_text(size=40))+
  scale_fill_manual(values=c( "#CCECE6", "#99D8C9"))

#-------------------------------------------------- ARM:
# use joint dataset
nrow(mydatw[mydatw$arm=="SP",])
nrow(mydatw[mydatw$arm=="SP/ART",])
# virtually 50:50 !
tab2 <- data.frame(
  Group=levels(mydatl$arm),
  value=c(nrow(mydatw[mydatw$arm=="SP",]),nrow(mydatw[mydatw$arm=="SP/ART",])))
tab2$value <- tab2$value/sum(tab2$value)*100
labels <- c(round(tab2$value[1]), round(tab2$value[2]))
labels <- c(paste(labels[1],"%"), paste(labels[2],"%"))
ggplot(tab2, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.5), size=20) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 45, face = "bold"), legend.text=element_text(size=40))+
  scale_fill_manual(values=c( "#CCECE6", "#99D8C9"))

#--------------------------------------------------  GENDER:
# Use joint dataset
length(unique(mydatl[mydatl$gender=="F",3]))
length(unique(mydatl[mydatl$gender=="M",3]))
tab3 <- data.frame(
  Group=levels(mydatl$gender),
  value=c(nrow(mydatw[mydatw$gender=="F",]),nrow(mydatw[mydatw$gender=="M",])))
tab3$value <- tab3$value/sum(tab3$value)*100
labels <- c(round(tab3$value[1]), round(tab3$value[2]))
labels <- c(paste(labels[1],"%"), paste(labels[2],"%"))
ggplot(tab3, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.5), size=20) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 45, face = "bold"), legend.text=element_text(size=40))+
  scale_fill_manual(values=c( "#CCECE6", "#99D8C9"))

#--------------------------------------------------  OBSERVATIONS PER DAY:
# Pyrconcentration:
length(unique(dfl_pyr$pid)) # only 302 patients 
nrow(dfl_pyr)
df <- data.frame(Present = round(c(nrow(dfl_pyr[dfl_pyr$day==0,])/nrow(dfw_pyr),nrow(dfl_pyr[dfl_pyr$day==1,])/nrow(dfw_pyr),
                                   nrow(dfl_pyr[dfl_pyr$day==2,])/nrow(dfw_pyr),nrow(dfl_pyr[dfl_pyr$day==3,])/nrow(dfw_pyr),
                                   nrow(dfl_pyr[dfl_pyr$day==7,])/nrow(dfw_pyr),nrow(dfl_pyr[dfl_pyr$day==14,])/nrow(dfw_pyr),
                                   nrow(dfl_pyr[dfl_pyr$day==21,])/nrow(dfw_pyr),nrow(dfl_pyr[dfl_pyr$day==28,])/nrow(dfw_pyr),
                                   nrow(dfl_pyr[dfl_pyr$day==42,])/nrow(dfw_pyr)),2),
                 Day = c("0","1","2","3","7","14","21","28","42"))
ggplot(data=df, aes(x=reorder(Day,c(1,2,3,4,5,6,7,8,9)), y=Present)) + 
  xlab("Day")+
  geom_bar(stat="identity")+
  geom_text(aes(label=Present), vjust=1.7, color="white")+
  theme_minimal(base_size = 22)
# Sulfconcentration
length(unique(dfl_sulf$pid)) # only 304 patients 
nrow(dfl_sulf)
df <- data.frame(Present = round(c(nrow(dfl_sulf[dfl_sulf$day==0,])/nrow(dfw_sulf),nrow(dfl_sulf[dfl_sulf$day==1,])/nrow(dfw_sulf),
                                   nrow(dfl_sulf[dfl_sulf$day==2,])/nrow(dfw_sulf),nrow(dfl_sulf[dfl_sulf$day==3,])/nrow(dfw_sulf),
                                   nrow(dfl_sulf[dfl_sulf$day==7,])/nrow(dfw_sulf),nrow(dfl_sulf[dfl_sulf$day==14,])/nrow(dfw_sulf),
                                   nrow(dfl_sulf[dfl_sulf$day==21,])/nrow(dfw_sulf),nrow(dfl_sulf[dfl_sulf$day==28,])/nrow(dfw_sulf),
                                   nrow(dfl_sulf[dfl_sulf$day==42,])/nrow(dfw_sulf)),2),
                 Day = c("0","1","2","3","7","14","21","28","42"))
ggplot(data=df, aes(x=reorder(Day,c(1,2,3,4,5,6,7,8,9)), y=Present)) + 
  xlab("Day")+
  geom_bar(stat="identity")+
  geom_text(aes(label=Present), vjust=1.7, color="white")+
  theme_minimal(base_size = 22)



# CONTINOUS VARIABLES:

#-------------------------------------------------- Box plots for Pyrimethamine and Sulfadoxine:
ggplot(dfl_pyr,aes(x=Pyrconcentration)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  xlab("Pyrimethamine Concentration")+
  theme_minimal(base_size = 27)
ggplot(dfl_sulf,aes(x=Sulfconcentration)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  xlab("Sulfadoxine Concentration")+
  theme_minimal(base_size = 27)

# think should probably do a log transformation since these variables are very skew
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


#-------------------------------------------------- Logged Baseline Parasite Density
ggplot(data = mydatw,
       aes(x = 1:nrow(mydatw), y = parbase)) +
  geom_point(size=2) + 
  ylab("Logged baseline parasite density")+
  xlab("Index")+
  theme_minimal(base_size = 22)

#-------------------------------------------------- Weight
ggplot(data = mydatw,
       aes(x = 1:nrow(mydatw), y = weight)) +
  geom_point(size=2) + 
  ylab("Weight")+
  xlab("Index")+
  theme_minimal(base_size = 22)

#-------------------------------------------------- Age
ggplot(data = mydatw,
       aes(x = 1:nrow(mydatw), y = age)) +
  geom_point(size=2) + 
  ylab("Age")+
  xlab("Index")+
  theme_minimal(base_size = 22)

#-------------------------------------------------- Pyrimethamine Concentration
ggplot(data = dfl_pyr,
       aes(x = 1:nrow(dfl_pyr), y = Pyrconcentration)) +
  geom_point(size=2) + 
  ylab("Logged Pyrimethamine concentration")+
  xlab("Index")+
  theme_minimal(base_size = 22)

#-------------------------------------------------- Sulfadoxine Concentration
ggplot(data = dfl_sulf,
       aes(x = 1:nrow(dfl_sulf), y = Sulfconcentration)) +
  geom_point(size=2) + 
  ylab("Logged Sulfadoxine concentration")+
  xlab("Index")+
  theme_minimal(base_size = 22)


#-------------------------------------------------- Check normality 
# Logged Baseline Parasite Density
ggplot(mydatw, aes(x=parbase)) + 
  geom_density() +
  xlab("Logged Baseline Parasite Density")+
  ylab("Density")+ 
  theme_minimal(base_size = 22) 
# Weight
ggplot(mydatw, aes(x=weight)) + 
  geom_density() +
  xlab("Weight")+
  ylab("Density")+ 
  theme_minimal(base_size = 22) 
# Age 
ggplot(mydatw, aes(x=age)) + 
  geom_density() +
  xlab("Age")+
  ylab("Density")+ 
  theme_minimal(base_size = 22) 
# Pyrimethamine Concentration
ggplot(dfl_pyr, aes(x=Pyrconcentration)) + 
  geom_density() +
  xlab("logged Pyrimethamine Concentration")+
  ylab("Density")+ 
  theme_minimal(base_size = 22) 
# Sulfadoxine Concentration
ggplot(dfl_sulf, aes(x=Sulfconcentration)) + 
  geom_density() +
  xlab("Logged Sulfadoxine Concentration")+
  ylab("Density")+ 
  theme_minimal(base_size = 22) 




#--------------------------------------------------------------
# BIVARIATE EDA
#--------------------------------------------------------------

# FACTOR VARIABLES 

# need day as a factor variable to plot
dfl <- mydatl
dfl$day <- as.factor(dfl$day)  
levels(dfl$day) <- c(0,1,2,3,7,14,21,28,42) # prefer to use actual days for time
pyr_fac <- dfl_pyr
pyr_fac$day <- as.factor(dfl_pyr$day)
levels(pyr_fac$day) <- c(0,1,2,3,7,14,21,28,42)
sulf_fac <- dfl_sulf
sulf_fac$day <- as.factor(dfl_sulf$day)
levels(sulf_fac$day) <- c(0,1,2,3,7,14,21,28,42)


#-------------------------------------------------- PYRIMETHAMINE AND SULFADOXINE
df <- data.frame(Drug = c(rep("Pyrimethamine", nrow(mydatl)), rep("Sulfadoxine", nrow(mydatl))), 
                 Concentration = c(mydatl$Pyrconcentration, mydatl$Sulfconcentration),
                 Day = c(dfl$day, dfl$day))
df$Drug <- as.factor(df$Drug)
df$Day <- as.factor(df$Day)
Drug <- df$Drug
ggplot(df,aes(x=Day,y=Concentration,fill=Drug)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  facet_wrap(~Drug) +
  ylab("Logged Drug Concentration")+
  xlab("Day")+
  theme_minimal(base_size = 27)+
  coord_cartesian(ylim = c(0, 3)) 
ggplot() +
  aes(x = df$Day, color = Drug, group = Drug, y = df$Concentration) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  ylab("Mean Log Drug Concentration")+
  xlab("Day") +
  theme_minimal(base_size = 27) +
  coord_cartesian(ylim = c(0, 3))


#-------------------------------------------------- Pyrimethamine Concentration - treatment arm
# Box-plots
ggplot(pyr_fac,aes(x=day,y=Pyrconcentration,fill=arm)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  facet_wrap(~arm) +
  ylab("Logged Pyrimethamine Concentration")+
  xlab("Day")+
  theme_minimal(base_size = 27)+
  coord_cartesian(ylim = c(0, 3)) 
# Interaction Plot 
Arm <- pyr_fac$arm
ggplot() +
  aes(x = pyr_fac$day, color = Arm, group = Arm, y = pyr_fac$Pyrconcentration) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  ylab("Mean Log Pyrimethamine Concentration")+
  xlab("Day") +
  theme_minimal(base_size = 27)+
  coord_cartesian(ylim = c(0, 3))
#-------------------------------------------------- Sulfadoxine Concentration - treatment arm
# Box-plots
ggplot(sulf_fac,aes(x=day,y=Sulfconcentration,fill=arm)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  facet_wrap(~arm) +
  ylab("Logged Sulfadoxine Concentration")+
  xlab("Day")+
  theme_minimal(base_size = 27)+
  coord_cartesian(ylim = c(0, 3)) 
# Interaction Plot 
Arm <- sulf_fac$arm
ggplot() +
  aes(x = sulf_fac$day, color = Arm, group = Arm, y = sulf_fac$Sulfconcentration) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  ylab("Mean Log Sulfadoxine Concentration")+
  xlab("Day") +
  theme_minimal(base_size = 27) +
  coord_cartesian(ylim = c(0, 3))


#-------------------------------------------------- Pyrimethamine Concentration - site
# Box-plots
ggplot(pyr_fac,aes(x=day,y=Pyrconcentration,fill=site)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  facet_wrap(~site) +
  ylab("Logged Pyrimethamine Concentration")+
  xlab("Day")+
  theme_minimal(base_size = 30)+
  coord_cartesian(ylim = c(0, 3))
# Interaction Plot 
Site <- pyr_fac$site
ggplot() +
  aes(x = pyr_fac$day, color = Site, group = Site, y = pyr_fac$Pyrconcentration) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  ylab("Mean Log Pyrimethamine Concentration")+
  xlab("Day") +
  theme_minimal(base_size = 30) +
  coord_cartesian(ylim = c(0,3))
#-------------------------------------------------- Sulfadoxine Concentration - site
# box-plots
ggplot(sulf_fac,aes(x=day,y=Sulfconcentration,fill=site)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  facet_wrap(~site) +
  ylab("Logged Sulfadoxine Concentration")+
  xlab("Day")+
  theme_minimal(base_size = 30)+
  coord_cartesian(ylim = c(0, 3)) 
# Interaction Plot 
Site <- sulf_fac$site
ggplot() +
  aes(x = sulf_fac$day, color = Site, group = Site, y = sulf_fac$Sulfconcentration) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  ylab("Mean Log Sulfadoxine Concentration")+
  xlab("Day") +
  theme_minimal(base_size = 30) +
  coord_cartesian(ylim = c(0, 3))


#-------------------------------------------------- Pyrimethamine Concentration - gender
# box-plots
ggplot(pyr_fac,aes(x=day,y=Pyrconcentration,fill=gender)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  facet_wrap(~gender) +
  ylab("Logged Pyrimethamine Concentration")+
  xlab("Day")+
  theme_minimal(base_size = 27) +
  coord_cartesian(ylim = c(0, 3))
# Interaction Plot 
Sex <- pyr_fac$gender
ggplot() +
  aes(x = pyr_fac$day, color = Sex, group = Sex, y = pyr_fac$Pyrconcentration) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  ylab("Mean Log Pyrimethamine Concentration")+
  xlab("Day") +
  theme_minimal(base_size = 27) +
  coord_cartesian(ylim = c(0, 3))
#-------------------------------------------------- Sulfadoxine Concentration - gender
# box-plots
ggplot(sulf_fac,aes(x=day,y=Sulfconcentration,fill=gender)) + 
  geom_boxplot(alpha=0.7, show.legend = F) + 
  facet_wrap(~gender) +
  ylab("Logged Sulfadoxine Concentration")+
  xlab("Day")+
  theme_minimal(base_size = 27)+
  coord_cartesian(ylim = c(0, 3)) 
# Interaction Plot 
Sex <- sulf_fac$gender
ggplot() +
  aes(x = sulf_fac$day, color = Sex, group = Sex, y = sulf_fac$Sulfconcentration) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  ylab("Mean Log Sulfadoxine Concentration")+
  xlab("Day") +
  theme_minimal(base_size = 27) +
  coord_cartesian(ylim = c(0, 3))







#-------------------------------------------------- CONTINOUS VARIABLES 


corr <- cor(dfw_pyr[,-c(1,2,3,4)],  use = "pairwise.complete.obs")
ggcorrplot(corr, 
           hc.order = F, 
           type = "lower",
           show.legend = F,
           lab = T, 
           tl.cex = 20, 
           lab_size = 5) 

corr <- cor(dfw_sulf[,-c(1,2,3,4)],  use = "pairwise.complete.obs")
ggcorrplot(corr, 
           hc.order = F, 
           type = "lower",
           show.legend = F,
           lab = T, 
           tl.cex = 20, 
           lab_size = 5) 



# Box Plots for baseline parasite density (Not really necessary !!!)
# (can do the same for the other continuous variables)
ggplot(mydatw,aes(x=arm,y=parbase,fill=arm)) + 
  #geom_jitter() +     # adds data points
  geom_boxplot(alpha=0.7, show.legend = F) + 
  ylab("Logged Baseline Parasite Density")+
  xlab("Arm")+
  theme_minimal(base_size = 22) # this argument increases text sizes consistently!
# Not much difference in the baseline parasite densities, which is what we want to see!

# Quickly check other continous variables 
ggplot(mydatw,aes(x=arm,y=weight,fill=arm)) + 
  #geom_jitter() +     # adds data points
  geom_boxplot(alpha=0.7, show.legend = F) + 
  ylab("Weight")+
  xlab("Arm")+
  theme_minimal(base_size = 22) 
ggplot(mydatw,aes(x=arm,y=age,fill=arm)) + 
  #geom_jitter() +     # adds data points
  geom_boxplot(alpha=0.7, show.legend = F) + 
  ylab("Age")+
  xlab("Arm")+
  theme_minimal(base_size = 22) 
