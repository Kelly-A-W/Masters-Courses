#--------------------------------------------------------------
# Exercise 2: SIR MODEL CONTINUED
#--------------------------------------------------------------

# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Characteristics of a population of interest:
# At the start of the outbreak, there is 1 infectious person in the population of 1000
# The number of contacts per day between a susceptible and infectious individual is 0.8.
# The average duration of infectiousness is 10 days. 



#-----------------------------------------
# Part 1: SIR model with drug therapy
#-----------------------------------------

# Using your SIR model from exercise 1 and running the model for 1000 days,
  # a.	Incorporate a birth and death rate of 1/800 persons per day.
  # b.	Incorporate a loss of immunity rate of 1/50 days 
  # c.	Incorporate drug therapy into your model in the following ways
    # i.	as a flow from I to R
    # ii.	as a flow from I to T
    # iii.	as a flow from S to I and T

# Use the following parameters: 
  # Tau: rate of treatment seeking =  1/2
  # r: drug recovery rate = 1/5 
  # p: proportion of the population who are treated = 30% 
  # gamma: natural recovery rate = 1/20

# Assume that Treated individuals are not infectious. 
# Comment on the differences in the models in terms of R0, proportion of infections that are treated. 



#------------------------------------- Part (1.a):

# Incorporate a birth and death rate of 1/800 persons per day.

library(deSolve)


# birth + death rate SIR 
sir_br_dr <- function(t, x, parms)  { 
  with(as.list(c(parms, x)), { 
    N = S + I + R   # NEED TO INCLUDE THIS BECAUSE POPULATION IS NOT CONSTANT (although her our birth and death rate is the same so our population total is actually constant)
    
    dS <- b*N - beta*I/N*S - mu*S  # b = birth rate, mu = death rate     
    dI <- beta*I/N*S - gam*I - mu*I
    dR <- gam*I - mu*R
    
    output <- c(dS,dI,dR)       
    list(output)
  })
}


# Initial values
start <- c(S=999, I=1, R=0)
## The parameters 
parms <- c(beta=0.8, gam=1/20, b=1/800, mu=1/800)
## vector of timesteps
times <- seq(0, 1000, 1) 

run<-ode(times=times, y=start, func=sir_br_dr,parms=parms)
plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Birth and Death Rate Model", xlab = "Day", ylab = "Number of Individuals")
lines(run[,3], col="blue")
lines(run[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))


# Calculate R0:
as.numeric(parms["beta"]/(parms["gam"]+parms["mu"]))




#------------------------------------- Part (1.b):

# Incorporate a loss of immunity rate of 1/50 days

sir_loi <- function(t, x, parms)  { 
  with(as.list(c(parms, x)), { 
    N = S + I + R   # NEED TO INCLUDE THIS BECAUSE POPULATION IS NOT CONSTANT (although her our birth and death rate is the same so our population total is actually constant)
    
    dS <- b*N - beta*I/N*S - mu*S + rho*R  # rho = rate of loss of immunity   
    dI <- beta*I/N*S - gam*I - mu*I
    dR <- gam*I - mu*R - rho*R
    
    output <- c(dS,dI,dR)       
    list(output)
  })
}


# Initial values
start <- c(S=999, I=1, R=0)
## The parameters 
parms <- c(beta=0.8, gam=1/20, rho=1/50, b=1/800, mu=1/800)
## vector of timesteps
times <- seq(0, 1000, 1) 

run<-ode(times=times, y=start, func=sir_loi,parms=parms)

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Loss of Immunity Model", xlab = "Day", ylab = "Number of Individuals")
lines(run[,3], col="blue")
lines(run[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))





#------------------------------------- Part (1.c)(i):

# Incorporate drug therapy into your model as a flow from I to R
# r: drug recovery rate = 1/5 
# p: proportion of the population who are treated = 30% 

sir_treat1 <- function(t, x, parms)  { 
  with(as.list(c(parms, x)), { 
    N = S + I + R   # NEED TO INCLUDE THIS BECAUSE POPULATION IS NOT CONSTANT (although her our birth and death rate is the same so our population total is actually constant)
    
    dS <- b*N - beta*I/N*S - mu*S + rho*R
    dI <- beta*I/N*S - p*r*I - (1-p)*gam*I - mu*I
    dR <- p*r*I + (1-p)*gam*I - mu*R - rho*R
    
    output <- c(dS,dI,dR)       
    list(output)
  })
}


# Initial values
start <- c(S=999, I=1, R=0)
## The parameters 
parms <- c(beta=0.8, gam=1/20, r=1/(2+5), p=0.3, rho=1/50, b=1/800, mu=1/800) # NB check notes for explanation on r
## vector of timesteps
times <- seq(0, 1000, 1) 

run<-ode(times=times, y=start, func=sir_treat1,parms=parms)

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Treatment 1", xlab = "Day", ylab = "Number of Individuals")
lines(run[,3], col="blue")
lines(run[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))

# Calculate R0:
R0_i <- as.numeric(parms["beta"]/(parms["mu"] + (1-parms["p"])*parms["gam"] + parms["p"]*parms["r"]))
R0_i

# Proportion of infected treated:
prop_treat <- as.numeric(parms["p"]*parms["r"]/(parms["p"]*parms["r"] + (1-parms["p"])*parms["gam"]))
prop_treat

#------------------------------------- Part (1.c)(ii):

# Incorporate drug therapy into your model as a flow from I to T
# r: drug recovery rate = 1/5 
# p: proportion of the population who are treated = 30% 
# Tau: rate of treatment seeking =  1/2  (people getting treatment when sick)

sir_treat2 <- function(t, x, parms)  { 
  with(as.list(c(parms, x)), { 
    N = S + I + R   # NEED TO INCLUDE THIS BECAUSE POPULATION IS NOT CONSTANT (although her our birth and death rate is the same so our population total is actually constant)
    
    dS <- b*N - beta*I/N*S - mu*S + rho*R
    dI <- beta*I/N*S - (1-p)*gam*I - p*tau*I - mu*I 
    dT <- p*tau*I - r*Treat - mu*Treat
    dR <- (1-p)*gam*I + r*Treat - mu*R - rho*R
    
    output <- c(dS,dI,dT,dR)       
    list(output)
  })
}

# Initial values
start <- c(S=999, I=1, Treat=0, R=0)
## The parameters 
parms <- c(beta=0.8, gam=1/20, r=1/5, p=0.3, tau=1/2, rho=1/50, b=1/800, mu=1/800)
## vector of timesteps
times <- seq(0, 1000, 1) 

run<-ode(times=times, y=start, func=sir_treat2, parms=parms)

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Treatment 2", xlab = "Day", ylab = "Number of Individuals")
lines(run[,3], col="blue")
lines(run[,4], col="black")
lines(run[,5], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))


# Calculate R0:
R0_ii <- as.numeric(parms["beta"]/(parms["mu"] + (1-parms["p"])*parms["gam"] + parms["p"]*parms["tau"]))
R0_ii

# Proportion of infected treated:
prop_treat <- as.numeric(parms["p"]*parms["tau"]/(parms["p"]*parms["tau"] + (1-parms["p"])*parms["gam"]))
prop_treat

#------------------------------------- Part (1.c)(iii):

# Incorporate drug therapy into your model as a flow from S to I and T
# r: drug recovery rate = 1/5 
# p: proportion of the population who are treated = 30% 

sir_treat3 <- function(t, x, parms)  { 
  with(as.list(c(parms, x)), { 
    N = S + I + R   # NEED TO INCLUDE THIS BECAUSE POPULATION IS NOT CONSTANT (although her our birth and death rate is the same so our population total is actually constant)
    
    dS <-  b*N - beta*I/N*S - mu*S + rho*R 
    dI <- (1-p)*beta*I/N*S - gam*I - mu*I 
    dT <- p*beta*I/N*S - r*Treat - mu*Treat
    dR <- gam*I + r*Treat - mu*R - rho*R
    
    
    output <- c(dS,dI,dT,dR)       
    list(output)
  })
}

# Initial values
start <- c(S=999, I=1, Treat=0, R=0)
## The parameters 
parms <- c(beta=0.8, gam=1/20, r=1/5, p=0.3, rho=1/50, b=1/800, mu=1/800)
## vector of timesteps
times <- seq(0, 1000, 1) 

run<-ode(times=times, y=start, func=sir_treat3, parms=parms)

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Treatment 3", xlab = "Day", ylab = "Number of Individuals")
lines(run[,3], col="blue")
lines(run[,4], col="black")
lines(run[,5], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))


# Calculate R0:
R0_iii <- as.numeric((1-parms["p"])*parms["beta"]/(parms["mu"] + parms["gam"]))
R0_iii







