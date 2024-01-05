#install.packages("pwr")


# Question 1(f)

library("pwr")

# Method 1
h= 2*asin(sqrt(0.02))-2*asin(sqrt(0.013))
pwr.2p.test(h=h, sig.level = 0.05, power = 0.80, alternative = "greater")


# Method 2 - not independent of sigma so wrong
library("daewr")
Fpower1(alpha=0.05, nlev=2, nreps=11000, Delta=0.65*0.02, sigma=1)

#################################################################
# Method 3 # correct !
(((1.96+0.845)/(asin(sqrt(0.013))-asin(sqrt(0.02))))^2)/2
#################################################################

# Method 4
2*(((1.645+0.845)/(0.35*0.02-0.035))^2)*0.02*(1-0.02)

# Method 5
h= 2*asin(sqrt(0.02))-2*asin(sqrt(0.35*0.02))
pwr.2p.test(h=0.35*0.02, sig.level = 0.05, power = 0.80, alternative = "greater")


# Method 6
power.prop.test(p1=0.02, p2=0.013, power=0.8, alternative = "two.sided")


# Method 7
pwr.2p.test(h=0.013, sig.level = 0.05, power = 0.80, alternative = "greater")


##########################
# Question 1(k)
##########################

sampsize <- function(rho, rho12, m){
  ((1.96 + 0.845)^2)*((0.02*(1-0.02) + 0.013*(1-0.013))/(0.02-0.013)^2)*(1+(m-1)*rho - m*rho12)/m
}

sampsize(0.015, 0.015, 200)
26*200*2


# Sensitivity Analysis
ICC = seq(0.015, (0.002+0.03), 0.0001)
IPC = rep(0.015, length(ICC))
M = rep(200, length(ICC))
plot(y=sampsize(ICC, IPC, M), ICC, pch=20, ylab = "Number of clusters (k)", cex.lab = 2, cex.axis = 1.5)
grid()


IPC = seq(0.0018, 0.015, 0.0001)
ICC = rep(0.015, length(IPC))
M = rep(200, length(IPC))
plot(y=sampsize(ICC, IPC, M), IPC, pch=20, ylab = "Number of clusters (k)", cex.lab = 2, cex.axis = 1.5)
grid()

M = seq(100, 300, 5)
ICC = rep(0.015, length(M))
IPC = rep(0.015, length(M))
plot(y=sampsize(ICC, IPC, M), M, pch=20, ylab = "Number of clusters (k)", 
     xlab="Average cluster size per arm (m)", cex.lab = 2, cex.axis = 1.5)
grid()




