#--------------------------------------------------------------
# Exercise 1: Modelling Measles using Differential Equations
#--------------------------------------------------------------

# Measles is a highly contagious disease that commonly occurs in young children.  
# Initial signs and symptoms typically include fever, often greater than 40 °C (104.0 °F), cough, 
# runny nose, and inflamed eyes. 
# You are tasked with modelling a measles epidemic in a small town. 
# Getting measles confers permanent immunity. 

# Characteristics of a population of interest:
  # At the start of the outbreak, there is 1 infectious person in the population of 1000
  # The number of contacts per day between a susceptible and infectious individual is 0.8.
  # The average duration of infectiousness is 10 days. 




#------------------------
# Question 1.a: SIR MODEL
#------------------------
# Use the deSolve package in R to model the measles epidemic over 100 days. 

library(deSolve)


#SIR model for measles
sir <- function(t, x, parms)  { # t is time, parms is vector of parameters, x is states/compartments
  with(as.list(c(parms, x)), {  # refer to elements of parms and x by name in the function (like attach)
    dS <- -beta*I/1000*S        # I is the number of infectious people, thus I/1000 is proportion of infectiousness
                                # (since for this example thesize of the total population is 1000)
    dR <- 1/gam*I               # 1 person recovers per 10 days 
    dI <- -dS - dR
    
    output <- c(dS,dI,dR)       
    list(output)
  })
}

# Initial values
start <- c(S=999, I=1, R=0)

## The parameters 
parms <- c(beta=0.8, gam=10)

## vector of timesteps
times <- seq(0, 100, 1)  # how long we want model to run for
# parameters need to be in same unit as this
# unit of time must be consistent
# here we choose 100 DAYS

run<-ode(times=times, y=start, func=sir,parms=parms)

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model")
lines(run[,3], col="blue")
lines(run[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))




#------------------------
# Question 1.b: SIR MODEL
#------------------------
# Manually assess the impact of changing the number of contacts per day and the average duration of infectiousness. 


# increase number of contacts per day:
parms <- c(beta=5, gam=10)
run<-ode(times=times, y=start, func=sir,parms=parms)

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model")
lines(run[,3], col="blue")
lines(run[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))
# Shifts functions to the left -  the disease epidemic is much shorter (population becomes recovered and thus immune 
# much quicker).
# Curves much steeper - the disease spreads quicker.
# The blue curve reaches a higher maximum - the maximum number of infected individuals at a given time reached is 
# greater, i.e. a greater proportion of the population will need to be infected before the rate of infection 
# starts to decrease.


# increase duration of infectiousness
parms <- c(beta=0.8, gam=20)
run<-ode(times=times, y=start, func=sir,parms=parms)

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model")
lines(run[,3], col="blue")
lines(run[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))
# Blue curve reaches a higher maximum - the maximum number of infected individuals at a given time reached is greater.
# Decreasing side of the blue curve not as steep - duration of infectiousness is longer.
# Green curve not as steep - takes longer for recovered population to increase - disease epidemic takes longer.
# Red curve slightly steeper - susceptible population decreases slightly faster




#------------------------
# Question 1.c: SIR MODEL
#------------------------
# Using a “for” loop, create a graph/video to show the impact of changing the number of contacts per day on the epidemic

cpd   <- seq(0.4,10,0.2)
start <- c(S=999, I=1, R=0)
times <- seq(0, 100, 1) 

parms <- c(beta=0.2, gam=10)
run <- ode(times=times, y=start, func=sir,parms=parms)
plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model")
lines(run[,3], col="blue")
lines(run[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))


for (i in cpd) {
  parms <- c(beta=i, gam=10)
  
  run <- ode(times=times, y=start, func=sir,parms=parms)
  lines(run[,2], col="red")
  lines(run[,3], col="blue")
  lines(run[,4], col="green")
  Sys.sleep(0.5)
  
}






#------------------------
# Question 2.a: SEIR MODEL
#------------------------
# Measles is often characterized by an SEIR model as it has a latent period of 12 days before an 
# individual becomes infectious.
# Use the deSolve package in R to model the measles epidemic over 100 days. 



#SEIR model for measles
seir <- function(t, x, parms)  { 
  with(as.list(c(parms, x)), { 
    
    dS <- -beta*I/1000*S        
    dR <- 1/gam*I       
    dE <- -dS - 1/sigma*E  # 1/sigma because 1 person per 12 days leaves latent period = rate of latency
    dI <- 1/sigma*E - dR
    
    output <- c(dS,dE,dI,dR)       
    list(output)
  })
}

# Initial values
start <- c(S=999, E=0, I=1, R=0)

## The parameters 
parms <- c(beta=0.8, gam=10, sigma=12)

## vector of timesteps
times <- seq(0, 100, 1) 

run<-ode(times=times, y=start, func=seir,parms=parms)

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SEIR Model")
lines(run[,3], col="black")
lines(run[,4], col="blue")
lines(run[,5], col="green")
legend("topright",legend=c("S", "E", "I", "R"), col=c("red", "black", "blue", "green"), lty=c(1,1,1,1))


# compare to SIR model:
start_sir <- c(S=999, I=1, R=0)
parms_sir <- c(beta=0.8, gam=10)
run_sir  <- ode(times=times, y=start_sir, func=sir,parms=parms_sir)
lines(run_sir[,2], col="red", type = "l", lty = "dashed")
lines(run_sir[,3], col="blue", type = "l", lty = "dashed")
lines(run_sir[,4], col="green", type = "l", lty = "dashed")
legend("right",legend=c("SEIR", "SIR"), lty=c(1,2))
# all curves are much less steep for SEIR plot - disease does not progress as fast
# smaller maximum reached for blue curve - maximum number of infected people reached is less for SEIR


#------------------------
# Question 2.b: SEIR MODEL
#------------------------
# Manually assess the impact of changing the number of contacts per day and the average duration of infectiousness. 


# increase number of contacts per day:
parms_test <- c(beta=5, gam=10, sigma=12)
run_test   <- ode(times=times, y=start, func=seir,parms=parms_test)

plot(run_test[,2], col="red", ylim=c(0,1000), type="l", main="SEIR Model")
lines(run_test[,3], col="black")
lines(run_test[,4], col="blue")
lines(run_test[,5], col="green")
legend("topright",legend=c("S", "E", "I", "R"), col=c("red", "black", "blue", "green"), lty=c(1,1,1,1))

# compare to 0.8 contacts per day
lines(run[,2], col="red", type = "l", lty = "dashed")
lines(run[,3], col="black", type = "l", lty = "dashed")
lines(run[,4], col="blue", type = "l", lty = "dashed")
lines(run[,5], col="green", type = "l", lty = "dashed")
legend("right",legend=c("5 cpd", "0.8 cpd"), lty=c(1,2))

# Shifts functions to the left -  the disease epidemic is much shorter (population becomes recovered and thus immune 
# much quicker).
# Curves much steeper - the disease spreads quicker.
# The blue curve reaches a slightly higher maximum - the maximum number of infected individuals reached is greater
# Black curve reaches a much higher maximum - the maximum number of individuals with latent disease reached is greater


# increase duration of infectiousness
parms_test <- c(beta=0.8, gam=20, sigma=12)
run_test   <- ode(times=times, y=start, func=seir,parms=parms_test)

plot(run_test[,2], col="red", ylim=c(0,1000), type="l", main="SEIR Model")
lines(run_test[,3], col="black")
lines(run_test[,4], col="blue")
lines(run_test[,5], col="green")
legend("topright",legend=c("S", "E", "I", "R"), col=c("red", "black", "blue", "green"), lty=c(1,1,1,1))

# compare to 0.8 contacts per day
lines(run[,2], col="red", type = "l", lty = "dashed")
lines(run[,3], col="black", type = "l", lty = "dashed")
lines(run[,4], col="blue", type = "l", lty = "dashed")
lines(run[,5], col="green", type = "l", lty = "dashed")
legend("right",legend=c("20 days I", "10 days I"), lty=c(1,2))

# Red curve shifted left - susceptible population decreases and depletes much sooner
# Black curve shifted left - exposed population increases much sooner but gets depleted at roughly same time as before
# Blue curve shifted left and higher - infectious population increases quicker and decreases slower
# Green curve lower - recovered population takes longer to increase




#------------------------
# Question 2.c: SEIR MODEL
#------------------------
# Using a “for” loop, create a graph/video to show the impact of changing the number of contacts per day on the epidemic

cpd   <- seq(0.4,10,0.2)


parms <- c(beta=0.2, gam=10, sigma=12)
run<-ode(times=times, y=start, func=seir,parms=parms)
plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SEIR Model")
lines(run[,3], col="black")
lines(run[,4], col="blue")
lines(run[,5], col="green")
legend("topright",legend=c("S", "E", "I", "R"), col=c("red", "black", "blue", "green"), lty=c(1,1,1,1))


for (i in cpd) {
  parms <- c(beta=i, gam=10, sigma=12)
  
  run <- ode(times=times, y=start, func=seir,parms=parms)
  lines(run[,2], col="red")
  lines(run[,3], col="black")
  lines(run[,4], col="blue")
  lines(run[,5], col="green")
  Sys.sleep(0.5)
  
}

# more contacts, the quicker the disease spreads



