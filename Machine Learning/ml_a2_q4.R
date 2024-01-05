
rm(list = ls())


#-----------------
# Load Libraries
#-----------------
library('e1071')


#-----------------
# Load Data
#-----------------
dat    = read.csv('Digits2020.csv')


#-----------------
# Template
#-----------------

Y  =  dat[,1]              # is it a 3 or 7
X  =  as.matrix(dat[,-1])  # Each column is a different pixel

image(matrix(X[1,],28,28)) # dimensions of the plot is 28 x 28 (total of 28X28=784 pixels)
# plotting the first observation/number

# A function for finding K nearest neighbours: 
k_nn = function(X1,X2,k = 10)
{
   # Initialise:
   N1 = dim(X1)[1]    # number of observations for data frame X1
   N2 = dim(X2)[1]    # number of observations for data frame X2
   d  = dim(X1)[2]    # number of variables in data frame X1 (in our case variables = pixels)
   ones  = matrix(1,N1,1)  # column vector of size N1 of 1's
   inds  = matrix(0,N2,k)  # matrix of zeros with N2 rows and k columns
   edges = inds            # edges set to zero 
   
   
   for(i in 1:N2)     # for each observation in X2
   {
       dists     = sqrt(rowSums((ones%*%X2[i,] - X1)^2))  # sum of euclidean distances between observation i in X2
                                                          # and ALL observations in X1
       wh        = order(dists)[2:(k+1)]                  # selecting the k shortest distances, excluding distance 
                                                          # between i and itself
       inds[i,]  = wh                                     # store which observations are neighbours
       edges[i,] = dists[wh]                              # store distances to neighbours
   }
   return(list(edges = edges, neighbours = inds, k = k))
}

# Calculate K nearest neighbours:
K   = 35
res = k_nn(X,X,K)   
# gives which observations are nieghbours for each observation and what the distances are



#=================================
# Question 4(a)
#=================================




# A function that calculates 
IsoMap=function(res,d)   # d is the components we want
{
  N    = dim(res$edges)[1] # Number of observations
  Dg   = matrix(Inf,N,N)   # initially assume that all observations have no neighbours
  
  # Initialize Dg
  for(i in 1:N){
    Dg[i,res$neighbours[i,]] = res$edges[i,] 
    # creating a matrix where element ij is the distance between observation i and observation j iff j is a nearest
    # neighbour of i - otherwise the value of element ij is inf
  }
  
  diag(Dg) = 0                # because the distance between each observation and itself is zero
  asp = allShortestPaths(Dg)  # find the shortest (gedesic) path between each observation using Floyd's algorithm
  Dg  = asp$length            # Shortest geodesic paths between any two observations
  
  # Find Embedding
  N = dim(Dg)[1]
  D = Dg
  S = D^2
  H = diag(N)-1/N*matrix(1,N,N)
  
  tau = -0.5*H%*%S%*%H
  sv = svd(tau)
  
  # Extract eigenvalues and eigenvectors
  lambda <- sv$d   # eigenvalues (only the first N=515 because N<p)
  v      <- sv$v   # eigenvectors (only the first N=515 because N<p)
  
  # Embedding
  U <- matrix(NA,N,d)  
  for (i in 1:N) {
    for(j in 1:d){
      U[i,j] <- lambda[j]*v[i,j]
    }
    
  }
  
  return(list(U = U, v=v, lambda=lambda, Dg = Dg, d = d))
}


im <- IsoMap(res, d = 2)


#=================================
# Question 4(b)
#=================================

plot(im$lambda[1:10]^2)   # plotting eigenvalues 


plot(x = im$U[,1], y= im$U[,2], xlab = "First Component", ylab = "Second Component",pch = 19)


df <- data.frame(im$U, Y)
df$Y <- as.factor(df$Y)
ggplot(df, aes(X1, X2, group = Y, color = Y)) +
  geom_point(size=3) +
  scale_x_continuous("First Component")+
  scale_y_continuous(name="Second Component")+
  theme_minimal(base_size = 22)
# Look very well separated !!!!

