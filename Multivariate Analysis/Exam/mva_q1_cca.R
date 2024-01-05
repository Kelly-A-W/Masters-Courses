# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(CCA)


#-------------------
# Load Data
#-------------------

dat <- read.csv("PrimatesData.csv", check.names = FALSE)
head(dat)
str(dat)

#-------------------
# Data Wrangling
#-------------------

# Set Factor Variables
dat$genus <- as.factor(dat$genus)
dat$class <- as.factor(dat$class)
dat$classdigit <- as.factor(dat$classdigit)
# and GENUSID is just a character variable (like pid)

str(dat) # as desired

# scale continuous variables:
dat_scal <- scale(dat[,3:10])
dat_scal <- as.data.frame(dat_scal)
str(dat_scal) # looks good !


#-------------------
# CCA
#-------------------
# AIM:
# The aim of CCA is to divide your data into two groups of variables and to look fo linear combinations of variables in each group such that these new variables will be maximally correlated.
# These linear combinations are known as vector variates and CCA is simply a method of studing the linear relationship between these two vector variates.
# We are essentially trying to extract as much of the correlation that exists between the two original groups and deposit it in decresing correlation order into pairs of new variables. 
# The two groups we will be looking at are the indices variables and the angle variables.

# separate the data:
y <- dat_scal[,1:5]
x <- dat_scal[,6:8]

cca <- cc(x,y)
cca$cor
round(cca$xcoef,3)
round(cca$ycoef,3)

# Figure ? shows the correlations between the three pairs of vector variates that were formed. 
par(mfrow=c(1,3))
plot(cca$scores$yscores[,1], cca$scores$xscores[,1])
plot(cca$scores$yscores[,2], cca$scores$xscores[,2])
plot(cca$scores$yscores[,3], cca$scores$xscores[,3])


cc2 <- comput(x, y, cca)

round(cc2$corr.X.xscores,3)
round(cc2$corr.Y.yscores,3)
