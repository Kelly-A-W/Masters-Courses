# WEEK 1 CLASS EXERCISES

# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#----------------------------------
# Question 2
#----------------------------------
# One treatment with 4 levels
# Delta = 2*sigma

#--------------------------------------------- (a) 
# Calculate the power for various numbers of replicates r per treatment level
library(daewr)

rmin   <- 2
rmax   <- 12
alpha  <- rep(0.05, (rmax-rmin+1)) # choosing 5% significance level
sigma  <- 2                        # doesn't matter what this is  (cancels out)
nlev   <- 4
nreps  <- rmin:rmax
Delta  <- 2*sigma
(power <- Fpower1(alpha, nlev, nreps, Delta, sigma))   # CRD (only one treatment), so use Fpower1


#--------------------------------------------- (b) 
# Calculate the number of replicates needed to have 0.90 power of detecting a difference as large as Delta

# Can see from power output put that this would need to be at least 9 replicates


#--------------------------------------------- (c) 
# How would the results you got in (b) change if the number of treatment levels increased to 8 or 
# decreased to 2

# If the number of treatment levels increased, then more replicates would needed to produce the 
# same power, while if the number of treatment levels decreased, them fewer replicates would be needed
# to produce the same power

# We can show this explicitly.

# For 8 treatment levels:
Fpower1(alpha, 8, nreps, Delta, sigma) # 11 replicates are needed
# For 2 treatment levels:
Fpower1(alpha, 2, nreps, Delta, sigma) # 7 replicates are needed


#----------------------------------
# Question 3
#----------------------------------
# One treatment with 4 levels
# Delta = 2*sigma

#--------------------------------------------- (a) 
# If by blocking the experimental units into blocks it was believed that the variance of the 
# experimental error, σ², could be reduced by 50%, calculate the number of blocks that would be 
# required to have a power of 0.90 for detecting a maximum difference in treatment means as large as ∆.

re <- 2
(b <- 9/re) # where 5 is the number of replicates need for the CRD to have a power of 0.90

# Thus five blocks would be needed to achieve a power of 90%, given delta


#--------------------------------------------- (b) 
# If by using a Latin-square design the variance of the experimental error, σ², could be reduced 
# another 10%, determine the power for detecting treatment differences when using a 4 4 Latin-square
# design.

bmin <- 3
bmax <- 15
alpha <- rep(0.05, bmax - bmin + 1)
nblocks <- bmin:bmax 
t = 4
nu1 <- 4-1
nu2 <- (nblocks-2)*(t-1)
sigma2 <- 0.1*(0.5*(1))  # pretending variance = 1, it does not matter what the actual value is
css <- 2*sigma2
nc <- nblocks*t*css/sigma2  
Fpower(alpha, nu1, nu2, nc)



#----------------------------------
# Question 4
#----------------------------------

#--------------------------------------------- (a) 
# The experimental unit is the ball flip


#--------------------------------------------- (b) 
# Using the numbers 1, 2, and 3 to represent the levels of start angle and stop angle, and holding 
# the pivot height constant at its high level, make a randomized list of experiments for a 3 by 3 
# factorial experiment with r=2 replicates per cell.

D <- expand.grid(start = c(1,2,3), stop = c(1,2,3))  # 3x3
D <- rbind(D,D)  # r=2

# Randomise order of experiments
set.seed(123)
D <- D[order(sample(1:18)),]
D


#--------------------------------------------- (c) 
# If the variance of the experimental error in the measured distance was σ² =12 inches, calculate 
# the number of replicates you would need to have a power of 0.90 for detecting a difference in 10 
# inches in cell means.
rmin <- 2
rmax <- 10
sigma <- sqrt(12)
alpha <- 0.05
Delta <- 10
nlev <- 3*3
nreps <- c(rmin:rmax)
(power <- Fpower1(alpha, nlev, nreps, Delta, sigma))
# You would need r=6 replicates

#--------------------------------------------- (d) 
# Calculate the number of replicates you would need to have a power of 0.90 for detecting a
# difference of 24 inches in marginal means for either factor.

Delta = 24
nlev <- c(3,3)
(result <- Fpower2(alpha, nlev, nreps, Delta, sigma))
# only 2 replicates are needed 


#----------------------------------
# Question 5
#----------------------------------

#--------------------------------------------- (a) 
t = 10
k = 3
factorial(t)/(factorial(k)*factorial(t-k)) # 120 blocks

k = 4
factorial(t)/(factorial(k)*factorial(t-k)) # 210 blocks

k = 5
factorial(t)/(factorial(k)*factorial(t-k)) # 252 blocks

k = 6
factorial(t)/(factorial(k)*factorial(t-k)) # 120 blocks


#--------------------------------------------- (b) 

# See word Doc.


#--------------------------------------------- (c) 

BIBsize(10, 3)
BIBsize(10, 4)
BIBsize(10, 5)
BIBsize(10, 6)

library(AlgDesign)
BIB3 <- optBlock( ~ ., withinData = factor(1:10), blocksizes = rep(3, 30))
des3 <- BIB3$rows
dim(des3) <- NULL
des3 <- matrix(des3, nrow = 30, ncol = 3, byrow = TRUE)
des3

BIB4 <- optBlock( ~ ., withinData = factor(1:10), blocksizes = rep(4, 10))
des4 <- BIB4$rows
dim(des4) <- NULL
des4 <- matrix(des4, nrow = 10, ncol = 4, byrow = TRUE)
des4

BIB5 <- optBlock( ~ ., withinData = factor(1:10), blocksizes = rep(5, 18))
des5 <- BIB5$rows
dim(des5) <- NULL
des5 <- matrix(des5, nrow = 18, ncol = 5, byrow = TRUE)
des5

BIB6 <- optBlock( ~ ., withinData = factor(1:10), blocksizes = rep(6, 15))
des6 <- BIB6$rows
dim(des6) <- NULL
des6 <- matrix(des6, nrow = 15, ncol = 6, byrow = TRUE)
des6


#--------------------------------------------- (d) 

# The number of blocks will always be equal to the number of treatment levels
# So no matter what k and r are, b=10

library(agricolae)

treat <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
des <- design.cyclic(treat, k = 3, r = 3)
des$book
