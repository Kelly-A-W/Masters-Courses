# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#-----------------
# Load Libraries
#-----------------
library(ggplot2)
library(caret)


# Simulate a dataset
N     <- 50
set.seed(123)
x     <- runif(N, min=-1, max=1)
noise <- rnorm(N)   # need a different epsilon for each y
y     <- sin(pi*x) + noise





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

res = neural_net(X,Y,theta_rand,nu=0)






#------------------------------------
# 10-FOLD CROSS-VALIDATION
#------------------------------------
ind <- 1:nrow(z)
set.seed(123)
gr1 <- sample(ind, nrow(z)*0.1)
gr2 <- sample(ind[-gr1], nrow(z)*0.1)
gr3 <- sample(ind[-c(gr1, gr2)], nrow(z)*0.1)
gr4 <- sample(ind[-c(gr1, gr2, gr3)], nrow(z)*0.1)
gr5 <- sample(ind[-c(gr1, gr2, gr3, gr4)], nrow(z)*0.1)
gr6 <- sample(ind[-c(gr1, gr2, gr3, gr4, gr5)], nrow(z)*0.1)
gr7 <- sample(ind[-c(gr1, gr2, gr3, gr4, gr5, gr6)], nrow(z)*0.1)
gr8 <- sample(ind[-c(gr1, gr2, gr3, gr4, gr5, gr6, gr7)], nrow(z)*0.1)
gr9 <- sample(ind[-c(gr1, gr2, gr3, gr4, gr5, gr6, gr7, gr8)], nrow(z)*0.1)
gr10 <- sample(ind[-c(gr1, gr2, gr3, gr4, gr5, gr6, gr7, gr8, gr9)], nrow(z)*0.1)

fold <- matrix(NA, nrow = nrow(z), ncol = 1)
fold[gr1]  <- 1
fold[gr2]  <- 2
fold[gr3]  <- 3
fold[gr4]  <- 4
fold[gr5]  <- 5
fold[gr6]  <- 6
fold[gr7]  <- 7
fold[gr8]  <- 8
fold[gr9]  <- 9
fold[gr10] <- 10

cv_data <- data.frame(z, fold)
cv_data$fold <- as.factor(cv_data$fold)




# define 10-fold cv function - outputs the the sum of the cv errors
cv <- function(cv_X,Y,theta,nu){
  val_err <- 0
  for (i in 1:10) {
    train_x <- cv_X[-which(cv_X$fold == i), -11]
    train_y <- matrix(Y[-which(cv_X$fold == i),], nrow = 45, ncol = 1)
    test_x  <- cv_X[which(cv_X$fold == i), -11]
    test_y  <- matrix(Y[which(cv_X$fold == i),], nrow = 5, ncol = 1)
    res     <- neural_net(train_x,train_y,theta_rand,nu=nu)
    
    # define objective function
    obj = function(pars)
    {
      res = neural_net(train_x,train_y,pars,nu=nu)
      return(res$E2)
    }
    
    # fit
    fit     <- nlm(obj,theta_rand)
    pred    <- neural_net(test_x,test_y,fit$estimate,nu=nu)
    val_err <- val_err + pred$E1
  }
  return(val_err/10)
}

cv(cv_data, Y, theta_rand, nu=2)


#------------------------------------
# Iterate for different nu's
#------------------------------------

nu_seq <- seq(0.1, 10, 0.1)
cv_errs <- matrix(NA, nrow = length(nu_seq), ncol=1)
for (i in 1:length(nu_seq)) {
  nu = nu_seq[i]
  er <- cv(cv_data, Y, theta_rand, nu=nu)
  cv_errs[i,1] <- er
}



#------------------------------------
# Plot
#------------------------------------

df <- data.frame(cv_errs, nu_seq)
(opt_nu <- df$nu_seq[which.min(df$cv_errs)])   # Optimal Nu



ggplot()+
  geom_line(df, mapping = aes(nu_seq, cv_errs))+
  geom_vline(xintercept=opt_nu, linetype="dashed", color="red")+
  scale_x_continuous("lambda")+
  scale_y_continuous("Cross Valdation Errors")+
  theme_minimal(base_size = 22)









#------------------------------------
# Plot Optimal Function
#------------------------------------



X     = z
Y     = matrix(y, ncol = 1, nrow = length(y))

# define objective function
obj = function(pars)
{
  res = neural_net(X,Y,pars,nu=opt_nu)
  return(res$E2)
}

# fit
fit     <- nlm(obj,theta_rand)

# Create sequence of input data
domain <- seq(-1, 1, 0.01)
z_dom <- matrix(NA, nrow=length(domain), ncol = 10)
for (i in 1:10) {
  z_dom[,i] <- Legendre(domain,i)
}


# Generate Predictions
pred <- neural_net(z_dom, matrix(1, nrow=nrow(z_dom), ncol=1),  # putting in a random Y because we're not interested in error
                   fit$estimate, nu=opt_nu)$out




target = sin(pi*domain)
df <- data.frame(domain, pred, target)
df <- pivot_longer(df, cols = 2:3, names_to = "Y", values_to = "Value")
df$Y <- as.factor(df$Y)
df
pointdata <- data.frame(x = x, y=y)
ggplot()+
  geom_line(df,mapping =aes(domain, Value, group=Y, color=Y))+
  geom_point(data = pointdata, mapping = aes(x = x, y = y))+
  scale_x_continuous("x")+
  scale_y_continuous("y")+
  theme_minimal(base_size = 22)









