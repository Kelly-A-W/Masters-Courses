#--------------------------------------------------------------
# Exercise 3: Spatial heterogeneity 
#--------------------------------------------------------------

# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Using your SIRS model from exercise 2 (treatment model 2): 
  # a. Duplicate the 4 compartment model to create a second population/patch using initial values S1=900, I1=20, T1=10, R1=70, S2=900, I2=20,T2=10, R2=70. Make use of the following parameters: 
    # beta=0.4, effective contact rate
    # r=1/7, drug therapy cure rate
    # p=1/50, loss of immunity rate
    # mu=1/100, birth/death rate
    # ppi=0.3, probability of receiving treatment
    # a=1/50, natural cure rate. 

  # b.	Incorporate physical movement between the two patches at a rate of m=1/20 days. 

  # c.	Explore the impact of movement on the two patch model if 
    # i. 	beta=0.1 in Patch 2
    # ii.	ppi = 0.6 in Patch 1
    # iii.	I2=T2=0, S2=1000, in Patch 2


#---------------------------------------------(a)

sir_treat2_2pop <- function(t, x, parms)  { 
  with(as.list(c(parms, x)), { 
    N1=S1+I1+Treat1+R1
    N2=S2+I2+Treat2+R2
    
    # population 1
    dS1 <- b1*N1 - beta1*I1/N1*S1 - mu1*S1 + rho1*R1
    dI1 <- beta1*I1/N1*S1 - (1-p1)*gam1*I1 - p1*tau1*I1 - mu1*I1 
    dT1 <- p1*tau1*I1 - r1*Treat1 - mu1*Treat1
    dR1 <- (1-p1)*gam1*I1 + r1*Treat1 - mu1*R1 - rho1*R1
    
    # population 2
    dS2 <- b2*N2 - beta2*I2/N2*S2 - mu2*S2 + rho2*R2
    dI2 <- beta2*I2/N2*S2 - (1-p2)*gam2*I2 - p2*tau2*I2 - mu2*I2 
    dT2 <- p2*tau2*I2 - r2*Treat2 - mu2*Treat2
    dR2 <- (1-p2)*gam2*I2 + r2*Treat2 - mu2*R2 - rho2*R2
    
    output <- c(dS1,dI1,dT1,dR1,dS2,dI2,dT2,dR2)       
    list(output)
  })
}


start <- c(S1=900, I1=20, Treat1=10, R1=70, S2=900, I2=20, Treat2=10, R2=70)
parms <- c(beta1=0.4, gam1=1/50, r1=1/7, p1=0.3, tau1=1/2, rho1=1/50, b1=1/800, mu1=1/800,
           beta2=0.4, gam2=1/50, r2=1/7, p2=0.3, tau2=1/2, rho2=1/50, b2=1/800, mu2=1/800)
times <- seq(0, 1000, 1) 
run   <- ode(times=times, y=start, func=sir_treat2_2pop,parms=parms)

attr(run, "dimnames")

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Population 1", xlab = "Day", ylab = "Number of Individuals")
lines(run[,3], col="blue")
lines(run[,4], col="black")
lines(run[,5], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))


plot(run[,6], col="red", ylim=c(0,1000), type="l", main="SIR Population 2", xlab = "Day", ylab = "Number of Individuals")
lines(run[,7], col="blue")
lines(run[,8], col="black")
lines(run[,9], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))



#---------------------------------------------(b)
# Incorporate physical movement between the two patches at a rate of m=1/20 days.

sir_treat2_2pop_wm <- function(t, x, parms)  { 
  with(as.list(c(parms, x)), { 
    
    N1=S1+I1+Treat1+R1
    N2=S2+I2+Treat2+R2
    
    # population 1
    dS1 <- b1*N1 - beta1*I1/N1*S1 - mu1*S1 + rho1*R1 + m21*S2 - m12*S1
    dI1 <- beta1*I1/N1*S1 - (1-p1)*gam1*I1 - p1*tau1*I1 - mu1*I1 + m21*I2 - m12*I1
    dT1 <- p1*tau1*I1 - r1*Treat1 - mu1*Treat1 + m21*Treat2 - m12*Treat1
    dR1 <- (1-p1)*gam1*I1 + r1*Treat1 - mu1*R1 - rho1*R1 + m21*R2 - m12*R1
    
    # population 2
    dS2 <- b2*N2 - beta2*I2/N2*S2 - mu2*S2 + rho2*R2 - m21*S2 + m12*S1
    dI2 <- beta2*I2/N2*S2 - (1-p2)*gam2*I2 - p2*tau2*I2 - mu2*I2 - m21*I2 + m12*I1
    dT2 <- p2*tau2*I2 - r2*Treat2 - mu2*Treat2 - m21*Treat2 + m12*Treat1
    dR2 <- (1-p2)*gam2*I2 + r2*Treat2 - mu2*R2 - rho2*R2 - m21*R2 + m12*R1
    
    output <- c(dS1,dI1,dT1,dR1,dS2,dI2,dT2,dR2)       
    list(output)
  })
}


start <- c(S1=900, I1=20, Treat1=10, R1=70, S2=900, I2=20, Treat2=10, R2=70)
parms <- c(beta1=0.4, gam1=1/50, r1=1/7, p1=0.3, tau1=1/2, rho1=1/50, b1=1/800, mu1=1/800, m12 = 1/20,
           beta2=0.4, gam2=1/50, r2=1/7, p2=0.3, tau2=1/2, rho2=1/50, b2=1/800, mu2=1/800, m21 = 1/20)
times <- seq(0, 1000, 1) 
run   <- ode(times=times, y=start, func=sir_treat2_2pop_wm, parms=parms)

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Population 1 With Movement", xlab = "Day", ylab = "Number of Individuals")
lines(run[,3], col="blue")
lines(run[,4], col="black")
lines(run[,5], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))

plot(run[,6], col="red", ylim=c(0,1000), type="l", main="SIR Population 2 With Movement", xlab = "Day", ylab = "Number of Individuals")
lines(run[,7], col="blue")
lines(run[,8], col="black")
lines(run[,9], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))



#---------------------------------------------(c)(i)
# c.	Explore the impact of movement on the two patch model if beta=0.1 in Patch 2

parmsi <- c(beta1=0.4, gam1=1/50, r1=1/7, p1=0.3, tau1=1/2, rho1=1/50, b1=1/800, mu1=1/800, m12 = 1/20,
           beta2=0.1, gam2=1/50, r2=1/7, p2=0.3, tau2=1/2, rho2=1/50, b2=1/800, mu2=1/800, m21 = 1/20)
runi   <- ode(times=times, y=start, func=sir_treat2_2pop_wm, parms=parmsi)

plot(runi[,2], col="red", ylim=c(0,1000), type="l", main="SIR Population 1 With Movement", xlab = "Day", ylab = "Number of Individuals")
lines(runi[,3], col="blue")
lines(runi[,4], col="black")
lines(runi[,5], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))

plot(runi[,6], col="red", ylim=c(0,1000), type="l", main="SIR Population 2 With Movement", xlab = "Day", ylab = "Number of Individuals")
lines(runi[,7], col="blue")
lines(runi[,8], col="black")
lines(runi[,9], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))




#---------------------------------------------(c)(ii)
# c.	Explore the impact of movement on the two patch model if ppi = 0.6 in Patch 1


parmsii <- c(beta1=0.4, gam1=1/50, r1=1/7, p1=0.6, tau1=1/2, rho1=1/50, b1=1/800, mu1=1/800, m12 = 1/20,
            beta2=0.4, gam2=1/50, r2=1/7, p2=0.3, tau2=1/2, rho2=1/50, b2=1/800, mu2=1/800, m21 = 1/20)
runii   <- ode(times=times, y=start, func=sir_treat2_2pop_wm, parms=parmsii)

plot(runii[,2], col="red", ylim=c(0,1000), type="l", main="SIR Population 1 With Movement", xlab = "Day", ylab = "Number of Individuals")
lines(runii[,3], col="blue")
lines(runii[,4], col="black")
lines(runii[,5], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))

plot(runii[,6], col="red", ylim=c(0,1000), type="l", main="SIR Population 2 With Movement", xlab = "Day", ylab = "Number of Individuals")
lines(runii[,7], col="blue")
lines(runii[,8], col="black")
lines(runii[,9], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))





#---------------------------------------------(c)(iii)
# c.	Explore the impact of movement on the two patch model if I2=T2=0, S2=1000, in Patch 2

startiii <- c(S1=900, I1=20, Treat1=10, R1=70, S2=1000, I2=0, Treat2=0, R2=70)
runiii   <- ode(times=times, y=startiii, func=sir_treat2_2pop_wm, parms=parms)

plot(runiii[,2], col="red", ylim=c(0,1000), type="l", main="SIR Population 1 With Movement", xlab = "Day", ylab = "Number of Individuals")
lines(runiii[,3], col="blue")
lines(runiii[,4], col="black")
lines(runiii[,5], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))

plot(runiii[,6], col="red", ylim=c(0,1000), type="l", main="SIR Population 2 With Movement", xlab = "Day", ylab = "Number of Individuals")
lines(runiii[,7], col="blue")
lines(runiii[,8], col="black")
lines(runiii[,9], col="green")
legend("topright",legend=c("S", "I", "T", "R"), col=c("red", "blue", "black", "green"), lty=c(1,1,1))

# starting points differ but otherwise plots are identical









