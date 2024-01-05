# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#-----------------
# Load Libraries
#-----------------
library(ggplot2)



#=================================
# Question 2(a)
#=================================

# Simulate a dataset
N     <- 50
set.seed(123)
x     <- runif(N, min=-1, max=1)
noise <- rnorm(N)   # need a different epsilon for each y
y     <- sin(pi*x) + noise


domain <- seq(-1, 1, 0.01)
df <- data.frame(domain)
pointdata <- data.frame(x = x, y=y)
ggplot(df,aes(domain))+
  stat_function(fun= function(domain) sin(pi*domain))+ 
  geom_point(data = pointdata, 
             mapping = aes(x = x, y = y))+
  scale_x_continuous("x")+
  scale_y_continuous(name="y")+
  theme_minimal(base_size = 22)




#=================================
# Question 2(b)
#=================================

# Legendre polynomials
Legendre=function(x,q){     # q = order of the Legendre polynomial
  val=0
  for(k in 0:q){
    val=val+((x^k)*choose(q,k)*choose((q+k-1)/2,q))
  }
  return((2^q)*val)
}



# Transform input data:
z <- matrix(NA, nrow=length(x), ncol = 10)
for (i in 1:10) {
  z[,i] <- Legendre(x,i)
}


#------------------------------------
# Activation Functions
#------------------------------------
# Specify some appropriate activation functions:
# Hidden layers

# First Layer =  Output Layer 
sig1 = function(x)    
{
  x       # activation function
}



#------------------------------------
# Neural Net Function
#------------------------------------
# Write a function that evaluates a neural network. 
# Takes arguments:
# X     - Design matrix
# Y     - Response matrix
# theta - A vector of parameters: first the weights ordered by layer and index, then basis vectors
# nu    - A regularisation parameter 

neural_net = function(X,Y,theta,nu)
{
  # Infer dimensions:
  N = dim(X)[1]  # nxp design matrix
  p = dim(X)[2]  # number of parameters
  q = dim(Y)[2]  # dimension of response vector - here we have a single output per observation
  
  dims = c(p,q)
  
  # Populate weight and bias matrices:
  index = 1:(dims[1]*dims[2])             # indices for weight matrix 1:  1 to p*q
  W1    = matrix(theta[index],dims[1],dims[2])

  index = length(theta)           
  b1    = matrix(theta[index],dims[2],1)  # qx1 = 1x1
  
  ones = matrix(1,1,N)   # transpose of a column vector of 1s
  a0   = t(X)       
  
  # Evaluate the updating equation in matrix form
  a1 = sig1(t(W1)%*%a0+b1%*%ones)   # evaluating all observations SIMULTANEUOUSLY - the reason we're using matrix equations (so we can avoid loops)
  
  
  # Evaluate an appropriate objective function and return some predictions:
  # classification problem therefor use cross entropy objective
  # have to split calculation to avoid numerical underflow
  # NB this is not bulletpoof, numerical underflow can still happen!! We're just taking steps to try avoid it
  out        = t(a1)
  E1         <- 1/(2*N)*sum((Y-out)^2)
  E2         = E1 + nu/N*(sum(W1^2))
  
  # Return a list of relevant objects:
  return(list(out = out,a1 = a1, E1 = E1, E2 = E2))
}


# Tester 
X     = z
Y     = matrix(y, ncol = 1, nrow = length(y))
p     = dim(X)[2]
q     = dim(Y)[2]
(npars = p*q + q)
set.seed(123)
theta_rand = runif(npars, -1, 1)

res_1 = neural_net(X,Y,theta_rand,nu=0)
res_2 = neural_net(X,Y,theta_rand,nu=5)


#------------------------------------
# Objective Function
#------------------------------------
# Create an objective function and evaluate it at a random coordinate in the 
# parameter space:

obj_1 = function(pars)
{
  res = neural_net(X,Y,pars,nu=0)
  return(res$E2)
}


obj_2 = function(pars)
{
  res = neural_net(X,Y,pars,nu=5)
  return(res$E2)
}


obj_1(theta_rand)
obj_2(theta_rand)


#------------------------------------
# Fit Neural Net
#------------------------------------
# Fit the neural network using a standard optimizer in R:
betas_1 = nlm(obj_1,theta_rand)   # can  adjust using interlim 
betas_1   # theta hat 
# now put this in our neural net to get y hat
fit_1 = neural_net(X, Y, betas_1$estimate, nu=0)


# Fit the neural network using a standard optimizer in R:
betas_2 = nlm(obj_2,theta_rand)   # can  adjust using interlim 
betas_2   # theta hat 
# now put this in our neural net to get y hat
fit_2 = neural_net(X, Y, betas_2$estimate, nu=5)





#------------------------------------
# Plot
#------------------------------------

# Create sequence of input data
z_dom <- matrix(NA, nrow=length(domain), ncol = 10)
for (i in 1:10) {
  z_dom[,i] <- Legendre(domain,i)
}

# Create predictions
pred1 <- neural_net(z_dom, matrix(1, nrow=nrow(z_dom), ncol=1),  # putting in a random Y because we're not interested in error
                    betas_1$estimate, nu=0)$out
pred2 <- neural_net(z_dom, matrix(1, nrow=nrow(z_dom), ncol=1),  # putting in a random Y because we're not interested in error
                    betas_2$estimate, nu=5)$out


target = sin(pi*domain)
df <- data.frame(domain, pred1, pred2, target)
df <- pivot_longer(df, cols = 2:4, names_to = "Y", values_to = "Value")
df$Y <- as.factor(df$Y)
df
pointdata <- data.frame(x = x, y=y)
ggplot()+
  geom_line(df,mapping =aes(domain, Value, group=Y, color=Y))+
  geom_point(data = pointdata, mapping = aes(x = x, y = y))+
  scale_x_continuous("x")+
  scale_y_continuous("y")+
  theme_minimal(base_size = 22)









