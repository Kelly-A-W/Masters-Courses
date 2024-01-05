# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(Rtsne)


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


#-------------------
# RRR
#-------------------

# Split data into X (indices) and Y (angle measurements)
# X will be our input variables, and Y will be our response variables
(X <- as.matrix(dat[,3:7], nrow=nrow(dat), ncol=5, byrow=T))
(Y <- as.matrix(dat[,8:10], nrow=nrow(dat), ncol=2, byrow=T))

# covariance matrices 
(Signma_XX <- cov(X)) 
(Signma_YY <- cov(Y)) 

# transpose X and Y so that variables are rows and samples are columns
X <- t(X)
Y <- t(Y)

# create mean vectors mu for X and Y
(mux <- apply(X, 1, mean))  # row means (i.e. variable means)
(muy <- apply(Y, 1, mean))

# center X and Y (i.e. subtract variable means)
(xc <- apply(X, 1, function(X) (X-mean(X))))    # X_c
(yc <- apply(Y, 1, function(Y) (Y-mean(Y))))    # Y_c

# transpose X_c and Y_c back into variables are columns and samples are rows
xc <- t(xc)
yc <- t(yc)

# compute sample sigma matrices using equations 6.94-6.96 in Izenman textbook
(n <- nrow(dat))     # number of observations
(Sxx <- (1/(n-1))*tcrossprod(xc))   # estimates covariance matrix of X
(Syy <- (1/(n-1))*tcrossprod(yc))   # estimates covariance matrix of Y
(Syx <- (1/(n-1))*tcrossprod(yc,xc))   # estimates covariance matrix of YX
(Sxy <- t(Syx))  # estimates covariance matrix of XY

# now we want to find the eigenvectors and eigenvalues of \Gamma^{1/2} S_yx S_xx^{-1} S_xy \Gamma^{1/2}
# where \Gamma is some positive-definite symmetric matrix of weights
# however we will bet setting \Gamma is the identity matrix for regression
# so just need to calculate eigenvalues/vectors of S_yx S_xx^{-1} S_xy 
# first calculate S_yx S_xx^{-1} S_xy 
(SyxSxxSxy <- Syx %*% solve(Sxx) %*% Sxy)
# now find eigen-decomposition of S_yx S_xx^{-1} S_xy 
(eigen1 <- eigen(SyxSxxSxy))

# separate and name eigenvectors
(V <- eigen1$vectors)
(V1 <- V[,1])
(V2 <- V[,2])

# recall that the matrix of regression coefficients, with rank allowed to be deficient, is C
# and C = A%*%B

# Since \Gamma is the identity matrix, A = V
A <- V

# Since \Gamma is the identity matrix, B = V' S_yx S_xx^{-1}
B <- t(V)%*%Syx%*%solve(Sxx)

# Rank 1 solution (i.e. if we set C to have rank of 1, where full rank would be 3):
A1 <- matrix(V[,1], ncol = 1, byrow = T)
B1 <- matrix(B[1,], ncol = 5, byrow = T)
(C1 <- A1%*%B1)   # regression coefficient matrix with rank 1

# Rank 2 solution:
A2 <- matrix(V[,1:2], ncol = 2, byrow = F)
B2 <- matrix(B[1:2,], ncol = 5, byrow = F)
(C2 <- A2%*%B2)   # regression coefficient matrix with rank 1

# Rank 3 solution (i.e. FULL rank solution):
(C3 <- A%*%B)   # regression coefficient matrix with rank 1


# We can now calculate the estimate mean vectors for the regression equation
# this is just the sample mean of the Ys minus the regression coefficients %*% X sample means
mu1 <- muy-A1%*%B1%*%mux    # rank 1 mean
mu2 <- muy-A2%*%B2%*%mux    # rank 2 mean
mu1 <- muy-A%*%B%*%mux      # full rank mean

#----------------------------------------------------------------------------------------------------------------------------

# Will now begin algorithm for using the rank trace to assess the effective dimensionality of the multivariate regression
# Define C = 0 and estimated error covariance matrix See_0 <- Syy
C0 <- matrix(rep(0,15), nrow = 3, ncol = 5) 
See0 <- Syy
# now compute the estimated error covariance matrices See_t for ranks t = 1, 2 and 3
(See1 <- (1/n)*(yc-(C1%*%xc))%*%t(yc-(C1%*%xc))) 
(See2 <- (1/n)*(yc-(C2%*%xc))%*%t(yc-(C2%*%xc))) 
(See3 <- (1/n)*(yc-(C3%*%xc))%*%t(yc-(C3%*%xc)))  
# compute \Delta C_t for t=1,2,3
# formula: \Delta C_t = || \theta - C_t || / || \theta ||, where here \theta=C_3
(dc0 = (norm(C3-C0))/norm(C3))  
(dc1 = (norm(C3-C1))/norm(C3))
(dc2 = (norm(C3-C2))/norm(C3))
(dc3 = (norm(C3-C3))/norm(C3))
dc <- rbind(dc0, dc1, dc2, dc3)
# now compute \Delta See_t for ranks t = 1, 2 and 3
(ds0 <- (norm(See3-See0))/norm(See3-Syy))
(ds1 <- (norm(See3-See1))/norm(See3-Syy))
(ds2 <- (norm(See3-See2))/norm(See3-Syy))
(ds3 <- (norm(See3-See3))/norm(See3-Syy))
ds <- rbind(ds0, ds1, ds2, ds3)
# now make a scatter plot of the points (\Delta C_t, \Delta See_t)
plot_dat1 = data.frame("dc" = dc,
                       "ds" = ds,
                       "Rank" = as.factor(0:3))
ggplot(plot_dat1, aes(x=ds, y=ds, color = Rank)) +  # Rank Trace
  geom_point(size = 3) +
  theme(legend.position = "right")
# want the smallest rank for which both \Delta C_t and \Delta See_t are approximately zero
# these seems to be t=2

#----------------------------------------------------------------------------------------------------------------------------

# Now use lm function
X <- t(X)
Y <- t(Y)
mod1 <- lm(Y~X)
summary(mod1)

# lm coefficients match C3 coefficients
mod1
t(C3)
# the intercepts match mean for rank 3 regression
mu1






