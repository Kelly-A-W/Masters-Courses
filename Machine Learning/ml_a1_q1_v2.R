rm(list = ls(all = TRUE))

#---------------------------------
# Simulate "Observed" Data
#---------------------------------

# Define f, the true Target Model:
f.mod=function(x)
{
  return(sin(pi*x) + cos(2*pi*x) + sin(3*pi*x) + cos(5*pi*x))
}

# Simulate Data from f
set.seed(4)
simul=function(n)
{
  x    <- runif(n, min=-1.5, max=1.5)  # Generating Input Data        
  y    <- f.mod(x)   # Generating Output Data
  data <- cbind(x, y)
  return(data)
}


# Fix the parameters:
n=3

# Simulated data (with noise)
(dat <- simul(n=n)) 


#------------------------------------------------------------
# Data produced by true target function ( no noise )
#------------------------------------------------------------
# Generate "True" Responses from Target Distribution

x <- seq(-1.5,1.5,0.01)    # domain
y <- f.mod(x)              # true output (for domain x)


#---------------------------------
# Final Hypotheses 
#---------------------------------

# H0
g_0.mod=function(x, data)
{
  g_0 <- lm(data[,2]~1)
  return(g_0$coefficients[1] + 0*x)
}


# H1
g_1.mod=function(x, data) 
{
  g_1 <- lm(data[,2]~data[,1])
  return(g_1$coefficients[1]+(g_1$coefficients[2]*x))
}




#---------------------------------
# Plot Data
#---------------------------------
plot(c(-1.5,1.5),c(-4,3),type="n",main="In-sample Errors",xlab=expression(x),ylab=expression(y))
lines(x,y,type="l",lwd=2)  # plotting the true responses (no error/noise)
points(dat[,1],dat[,2], pch=19);      # plotting generated data (with noise)
abline(h=g_0.mod(x, dat), col="green")    # H0 model 
lines(x,g_1.mod(x, dat), col="blue")      # H1 model 





#===============================================================================
# Let's see the trade-off for varying model complexity:
#===============================================================================

# Fit M constant models to random data:
M = 10000
curves = matrix(rep(NA,M*length(x)),ncol=length(x),nrow=M)
avg    = rep(NA,length(x))
vari   = rep(NA,length(x))

for(i in 1:M){
  data       <- simul(n=3)  # each hypothesis can only be fit with 3 data points 
  curves[i,] <- g_0.mod(x,data)
}
curves  # an element in row i, column j is the prediction for input x_j using g_0.mod fitted to dataset i
# thus since we have a constant function, all elements in the same row are equal
# thus each row is different but all columns are equal
for(i in 1:length(x)){
  avg[i]  <- mean(curves[,i])  # this will be the same for all i because all columns are equal because g_o.mod is constant function
  vari[i] <- var(curves[,i])   
}

# In-sample Errors
plot(c(-1.5,1.5),c(-4,3),type="n",main="",xlab=expression(x),ylab=expression(y))
lines(x,y,type="l",col="black",lwd=2)  # plotting the true responses (no error/noise)
lines(x,avg,type="l",lwd=2,lty=2,col="red")
bias=(1/3)*(t(f.mod(x)-avg)%*%(f.mod(x)-avg))*0.01
lines(x, avg+sqrt(vari),type="l",lwd=0.5,lty=2,col="red")
lines(x, avg-sqrt(vari),type="l",lwd=0.5,lty=2,col="red")
c(bias, mean(vari), bias + mean(vari))
grid()



# fit M linear models to random data
M = 10000
curves = matrix(rep(NA,M*length(x)),ncol=length(x),nrow=M)
avg    = rep(NA,length(x))
vari   = rep(NA,length(x))
for(i in 1:M){
  data       <- simul(n=3)  # each hypothesis can only be fit with 3 data points 
  curves[i,] <- g_1.mod(x,data)
}
for(i in 1:length(x)){
  avg[i]  <- mean(curves[,i])  # mean of each column - i.e. mean response
  vari[i] <- var(curves[,i])
}

# In-sample Errors
plot(c(-1.5,1.5),c(-5,5),type="n",main="",xlab=expression(x),ylab=expression(y))
lines(x,y,type="l",col="black",lwd=2)  # plotting the true responses (no error/noise)
lines(x,avg,type="l",lwd=2,lty=2,col="red")
#abline(h=0)
bias=(1/3)*(t(f.mod(x)-avg)%*%(f.mod(x)-avg))*0.01
lines(x,avg+sqrt(vari),type="l",lwd=0.5,lty=2,col="red")
lines(x,avg-sqrt(vari),type="l",lwd=0.5,lty=2,col="red")
c(bias,mean(vari),bias+mean(vari))
grid()



#bias=function(x){return((abs(f.mod(x)))^2)}
#integrate(bias,-1.5,1.5)
























