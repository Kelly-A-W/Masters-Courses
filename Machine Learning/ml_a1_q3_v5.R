rm(list = ls(all = TRUE))

#------------------------------------
# Load Libraries
#------------------------------------

library(ggplot2)
library(tidyr)
library(viridis)


#---------------------------------
# QUESTION 3(a)
#---------------------------------

# Legendre polynomials
Legendre=function(x,q){     # q = order of the Legendre polynomial
  val=0
  for(k in 0:q){
    val=val+((x^k)*choose(q,k)*choose((q+k-1)/2,q))
  }
  return((2^q)*val)
}



# Creates random polynomial target functions: Equation (3)
rLegfunc=function(x,Q_f){      # Q_f = order of the polynomial f(x)
  val=0
  beta=runif((Q_f+1),-1,1)  
  for(q in 0:Q_f){
    val = val+(Legendre(x,q)*beta[q+1])  
    val = scale(val, center = F, scale = T)  # scaling betas so that E[f^2]=1
  }
  return(val)
}



# simulates a dataset with target y=f(x) of size n with noise sig
generator=function(n,x,y,sig){
  l       <- length(y)       
  dat     <- matrix(rep(NA,2*n),ncol=n)  # n x 2 matrix - storage 
  dat_ind <- floor(runif(n)*l)+1         # randomly choosing elements to extract from our data [x,y]
  # runif(n) generates n random numbers between 0-1 - this can be seen as proportion
  # runif(n)*l then generates n proportions of l
  # floor() rounds these numbers down - because to select elements we need whole numbers
  # +1 ensures that the smallest number is 1 and not 0
  ydat    <- y[dat_ind]+rnorm(n,0,sig)  # picking out the responses from the target function and adding noise
  xdat    <- x[dat_ind]                 # picking out the corresponding x values
  Data    <- data.frame(xdat,ydat)      # combine
  return(Data)
}




# 10-th order target with noise (plot of what one of the datasets could look like)
x<-seq(-1,1,0.01);n=15;sig=0.5;Q_f=10;
plot(c(min(x),max(x)),c(-4,4),main="10-th order target + noise",type="n",xlab="x",ylab=expression(f(x)))
set.seed(123)
y=rLegfunc(x,Q_f)
lines(x,y,type="l",lwd=2)  # true target function
d<-generator(n,x,y,sig)
points(d$xdat,d$ydat)      # simulated data with noise
mean(y^2)   # scaled as desired 


# Create sequence of N's and sigma's
N   <- seq(20, 110, 1)       # 20 ≤ N ≤ 110
K   <- 0.01                    # INCREASE RESOLUTION LATER
sig <- seq(0.2, (1.1-K), K)   # 0.2 ≤ σ < 1.1 

# creat a vector of all possible combinations of N and sig
lat     <- expand.grid(N,sig)    # each row is a different N*sig combo
names(lat) <- c("N", "sig")



# Do M times:
# Need to generate a separate dataset for each row of lat (i.e. each combination of N and sig).
# Then for each dataset: 
# fit g_2 and g_10
# Then calculate Eout for g_2 and g_10 separately
# Then calculate  the difference Eout(g10)-Eout(g2)
# Take average to get overfit measure E[Eout(g10)-Eout(g2)]
# should have 2*dim(lat)[1] overfit estimates



# Initialise:
x <- seq(-1,1,0.01)       # domain
set.seed(123)
y <- rLegfunc(x,Q_f=10)   # target function / true responses




# Calculate Responses for fitted Model:
fit=function(x,model){
  v=0
  for(i in 1:length(model$coefficient)){
    v=v+(model$coefficient[i]*(x^(i-1)))
  }
  return(v)
}

# gives the bias for a given fitted model (DOES NOT ACCOUNT FOR VARIABILITY IN THE DATA)
fdiff=function(x,target,model){
  f=fit(x,model)
  return((t(f-target)%*%(f-target))*(x[2]-x[1])*1/(max(x)-min(x)))  #  CHECK 1/(max(x)-min(x)) TERM !!!!!!!!!!!
}




M = 1000
sum.err.overfit = 0
set.seed(123)
for(i in 1:M){
  
  # Intialise
  dat <- list()
  g_2 <- list()
  g_10 <- list()
  err.out.2 <- c()
  err.out.10 <- c()
  
  # For each combo of N and sig:
  for (j in 1:(dim(lat)[1])) {
    
    # Generate dataset
    dat[[j]] <- generator(lat[j,1], x, y, lat[j,2])
    
    # Fit g_2 
    g_2[[j]] <- lm(dat[[j]]$ydat~dat[[j]]$xdat+I(dat[[j]]$xdat^2))
    # Fit g_10
    g_10[[j]] <- lm(dat[[j]]$ydat~dat[[j]]$xdat+I(dat[[j]]$xdat^2)+I(dat[[j]]$xdat^3)+I(dat[[j]]$xdat^4)+I(dat[[j]]$xdat^5)+I(dat[[j]]$xdat^6)+I(dat[[j]]$xdat^7)+I(dat[[j]]$xdat^8)+I(dat[[j]]$xdat^9)+I(dat[[j]]$xdat^10))
    
    # Estimate Eout for g2
    err.out.2[j] <- fdiff(x,y,g_2[[j]]) + (lat[j,2])^2 
    # Estimate Eout for g10
    err.out.10[j] <- fdiff(x,y,g_10[[j]]) + (lat[j,2])^2
  }
  
  
  Err.overfit <- err.out.10 - err.out.2
  
  sum.err.overfit <- sum.err.overfit + Err.overfit
  
}


# Expected Overfit:
(exp.err.overfit <- sum.err.overfit/M)









# Need to convert err.out's to matrices for contour plot function

# Z
overfit_trans <- log(exp.err.overfit - min(exp.err.overfit)+0.01) - log(-min(exp.err.overfit)-0.01)
mat      <- cbind(lat, overfit_trans)
mat$N    <- as.factor(mat$N)
mat$sig  <- as.factor(mat$sig)
z        <- matrix(NA, nrow = length(levels(mat$sig)),  
                   ncol = length(levels(mat$N)))   # cols are N, rows are Sig 
# set x = N, y = sigma
z
for (i in 1:length(levels(mat$N))) {        # columns / N
  for(j in 1: length(levels(mat$sig))){     # rows / sig
    a      <- which(mat$N == levels(mat$N)[i])
    b      <- which(mat$sig == levels(mat$sig)[j])
    ind    <- intersect(a,b)
    z[j,i] <- mat$overfit_trans[ind]
  }
  
}

z







# Contour Plot:
filled.contour(x=N, y=sig,  
               t(z), 
               xlab="N", ylab="Sigma",
               color.palette=viridis) 




