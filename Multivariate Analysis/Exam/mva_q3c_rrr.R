# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(rrr)


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
# RRR: PCA
#-------------------

X <- as.matrix(dat_scal, nrow = nrow(dat_scal), ncol = ncol(dat_scal))
pca_rrr <- rrr(X,
               X, 
               k = 0)
# extract the loading of each variable of the principal components
# rows correspond to variables and columns to principal components
# so the first column is the loadings of each variable on the first principal component
round(pca_rrr$A[,1:4], 3)

# Next we can calculate the cumulative proportion of variance explained by component i by taking the the sum of the the eigenvalues of the principal components 1 to i and 
# then dividing this over the sum of the eigenvalues of all the principal components. 
# This is equivalent to dividing variance explained by the i-dimensional approximation by the total variance.
(cp1 <- pca_rrr$eigen_values[1]/sum(pca_rrr$eigen_values))
(cp2 <- (pca_rrr$eigen_values[1]+pca_rrr$eigen_values[2])/sum(pca_rrr$eigen_values))
(cp3 <- (pca_rrr$eigen_values[1]+pca_rrr$eigen_values[2]+pca_rrr$eigen_values[3])/sum(pca_rrr$eigen_values))
# etc 


#-------------------
# RRR: CCA
#-------------------
indice <- dat_scal[,1:5] # X
angle <- dat_scal[,6:8] # Y
indice <- as.matrix(indice, nrow=nrow(indice), ncol=ncol(indice))
angle <- as.matrix(angle, nrow=nrow(angle), ncol=ncol(angle ))

cca_rrr <- rrr(angle, indice, k=0,
            type = "cva")
# correlations:
cca_rrr$canonical_corr
# correlations are slightly different


# coefficients
t(cca_rrr$G) # x coefficients 
t(cca_rrr$H) # y coefficients 
# coefficients different than before.
# This is due to the different scaling methods used by the two algorithms.




