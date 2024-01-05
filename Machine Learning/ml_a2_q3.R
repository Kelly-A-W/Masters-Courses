rm(list = ls(all = TRUE))

#-----------------
# Load Libraries
#-----------------
library(ggplot2)
library(pixmap)
library(dplyr)


#-----------------
# Load Data
#-----------------
x=read.pnm(file = "Face_Data/1.pgm")
#cols = colorRampPalette(c('grey90','grey50','grey10'))


#=================================
# Question 3(a)
#=================================


m  = prod(x@size) # number of pixels per image = number of variables
n  = 400          # number of observations 
xx = t(matrix(rev(x@grey),x@size[1],x@size[2]))  # x@grey = a 112x92 matrix/array with the brightness value for each pixel
# rev() converts this to a vector of length 112x92 = 10304 = m
# so xx is essentially converting x@grey from a matix/array to a matrix
# this is only for image 1 
image(xx)  # this is image 1
dims = x@size

# Create observation matrix
X = matrix(0, n, prod(x@size))   # each row is a different image, each column is a different pixel
for(i in 1:n)
{
  x     = read.pnm(file = paste0("Face_Data/",i,".pgm"))    # importing all our faces 
  X[i,] = rev(x@grey)
}



mean_x <- t(matrix(colMeans(X), x@size[1],x@size[2]))
image(mean_x)


sd_x <- t(matrix(apply(X, 2, sd), x@size[1],x@size[2]))
image(sd_x)


# Original 168
x_168 <- t(matrix(X[168,], x@size[1],x@size[2]))
image(x_168)


# Make Scaled Data
X_scal <- scale(X, center = FALSE, scale = apply(X, 2, sd, na.rm = TRUE))


# Scaled 168
x_scal_168 <- t(matrix(X_scal[168,], x@size[1],x@size[2]))
image(x_scal_168)




#=================================
# Question 3(b)
#=================================

# Calculate Eigenvalues:	
N   = dim(X)[1]
SIG = 1/N*(X_scal%*%t(X_scal))  # covariance matrix
sv  = eigen(SIG)      # eigenvectors

# Construct components:
temp  = t(X_scal)%*%sv$vectors 
temp  = temp*matrix(1/sqrt(sv$values),dim(X_scal)[2],dim(X_scal)[1], byrow = TRUE)

# Eigenfaces:
#par(mfrow = c(1,1))
for(i in 1:10){
  image(t(matrix(temp[,i],dims[1],dims[2])))
  }


#=================================
# Question 3(c)
#=================================


# plot scaled 115
x_scal_115 <- t(matrix(X_scal[115,], x@size[1],x@size[2]))
image(x_scal_115)


# Now, reconstruct 115 using 5, 50, 200 eigenfaces and compare construction 
# to full image:
ss = c(5,50,200)

ii = 115  
# par(mfrow = c(2,5))
for(i in ss){
  Xhat =  (X_scal[ii,]%*%temp[,1:i])%*%t(temp[,1:i])
  image(t(matrix(Xhat,dims[1],dims[2])))
  }






