# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(PPCI)


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
# PPCI
#-------------------
X <- as.matrix(dat_scal)

mod1 <- mddc(X, K=5)
plot(mod1)
nhc_dat <- data.frame("Class" = dat[,11], 
                      "Cluster" = as.factor(mod1$cluster))
table(nhc_dat)

mod2 <- mcdc(X, K=5)
plot(mod2)
nhc_dat <- data.frame("Class" = dat[,11], 
                      "Cluster" = as.factor(mod2$cluster))
table(nhc_dat)

mod3 <- ncutdc(X, K=5)
plot(mod3)
nhc_dat <- data.frame("Class" = dat[,11], 
                      "Cluster" = as.factor(mod3$cluster))
table(nhc_dat)

