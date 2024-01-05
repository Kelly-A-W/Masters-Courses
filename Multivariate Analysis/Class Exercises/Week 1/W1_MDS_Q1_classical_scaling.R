
# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load data
# install.packages('psychTools')
library(psychTools)
data(cities)
cities

#------------------------------
# Question 1.a:
#------------------------------
# Using the first 4 (four) cities, ATL, BOS, ORD and DCA. 
# Perform the classical scaling algorithm by first principals to obtain the
# principals coordinates. 
# Plot the principal coordinates in 2 dimensions. 

# extract first 4 cities
dissim <- cities[1:4,1:4]
dissim  # dissimilarity matrix of interpoint distances

# from matrix A
A <- -1/2*dissim^2
A <- as.matrix(A)
dimnames(A) <- NULL
A

# create H matrix
I <- diag(4)   # 4x4 identity matrix
I
J <- matrix(rep(1,16), nrow=4, ncol=4)  # 4x4 matrix of of 1's
J
H <- I-1/4*J
H

# obtaining the B matrix
B <- H%*%A%*%H
B

# obtain eigen values and vectors
EVD <- eigen(B)
eigen_val <- EVD$values[1:2] # only interested in the first 2 because our plot will be 2D
eigen_val <- sqrt(eigen_val)
eigen_val
eigen_val_mat <- matrix(c(eigen_val[1],0,0,eigen_val[2]), nrow=2, ncol=2)
eigen_vect <- EVD$vectors
eigen_vect <- eigen_vect[,1:2] # only want the first two

# obtaining principal components
principal_coord <- eigen_vect%*%eigen_val_mat
principal_coord

# plotting the principal coordinates in 2D
plot(principal_coord, pch=19, las=1, xlab="Dimension x", ylab="Dimension y", xlim=c(-500,450))
city_names <- c("ATL", "BOS", "ORD", "DCA")
text(principal_coord, pos=4, labels=city_names)
grid()


#------------------------------
# Question 1.b:
#------------------------------
# Use the function cmdscale in R to check check your answer in part (a).

check <- cmdscale(dissim)
check  # same as principal_coord

# plot cmdscale solution
require(graphics)
x <- check[,1]
y <- check[,2]
plot(x,y,type="n",xlab="Dimension x",ylab="Dimension y", asp=1, axes=T, main="cmdscale")
text(x,y,rownames(check),cex=0.6)


#------------------------------
# Question 1.c:
#------------------------------
# Perform the classical scaling algorithm to the entire dissimilarity matrix and 
# obtain the principal coordinates. 
# Plot the resulting coordinates in 2 dimensions. 

# Classical scaling using cmdscale
full <- cmdscale(cities)
x <- full[,1]
y <- full[,2]
plot(x,y,type="n",xlab="Dimension x",ylab="Dimension y", asp=1, axes=T, main="Classical Scaling on Full Data")
text(x,y,rownames(full),cex=0.6)



