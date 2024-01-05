
#-----------------
# Question 3a:
#-----------------
library("pwr")

power.prop.test(p1=0.30, p2=0.2, power=0.8, sig.level = 0.05, alternative = "two.sided")

h= 2*asin(sqrt(0.30))-2*asin(sqrt(0.2))
pwr.2p.test(h=h, sig.level = 0.05, power = 0.80, alternative = "two.sided")

pwr.2p.test(h=ES.h(p1=0.30, p2=0.2), sig.level = 0.05, power = 0.80, alternative = "two.sided")


#-----------------
# Question 3b:
#-----------------

# Plan a Factorial Design
D <- expand.grid(A = c("A1", "A2", "P"), B = c("A1", "P") )

# To replicate design
D <- rbind(D, D, D)   # 3 relicates

# Randomise order of runs
set.seed(123)
D <- D[order(sample(1:nrow(D))), ]
D


#-----------------
# Question 3c:
#-----------------

pwr.2p.test(h=ES.h(p1=0.30, p2=0.2), sig.level = 0.01, power = 0.80, alternative = "two.sided")




#------------------------------------------
# Question 3e:
#------------------------------------------
# We have a parallel design
# With 2 arms 
# and 600 patients

N = 600
trt <- c("A", "B")

#------------------------------------------ Completely Randomised Design:
# Initialise:
pids1 <- data.frame(PID = 1:N, Treat = NA)
D1 <- c(0)
guess <- c()
H  <-0 # hit
MT <-0 # miss or tie
Mi <-0 # miss
numA = 0
numB = 0

set.seed(12)
for(i in 1:N){
  
  # Simulate Clinicians Guess:
  if(D1[i]>0){
    guess[i] <- "B"
  }else if(D1[i]<0){
    guess[i] <- "A"
  }else{guess[i] <- "T"}   # T for Tie
  
  # Assign Treatment:
  pids1[i,2]  <- sample(trt, 1)
  
  # Update Imbalance:
  if(pids1[i,2] == "A"){numA=numA+1}else{numB=numB+1}
  D1[i+1]       <- numA - numB
  
  # Compare Guess to Actual Allocation:
  if(pids1[i,2] == guess[i]){
    H = H + 1
  }else{
      MT = MT + 1
      }
  
}
pids1
D1[601]
(M = MT - length(guess[guess=="T"]))
(biasF1 <- (H-M)/2)

# Plot Random Walk
x = seq(0,N,1)
plot(x, D1, type="l", xlab ="N", ylab="Dj", ylim = c(-23, 27), cex.lab = 2, cex.axis = 1.5)
abline(h=0, col = "red")
grid()



#------------------------------------------  Randomised Block Design:
# Initialise:
urn   <- rep(trt, N/2)
pids2 <- data.frame(PID = 1:N, Treat = NA)
D2    <- c(0)
guess <- c()
H  <-0 # hit
MT <-0 # miss or tie
Mi <-0 # miss
numA = 0
numB = 0

set.seed(10)
for(i in 1:N){
  
  # Simulate Clinicians Guess:
  if(D2[i]>0){
    guess[i] <- "B"
  }else if(D2[i]<0){
    guess[i] <- "A"
  }else{guess[i] <- "T"}   # T for Tie
  
  # Assign Treatment:
  pids2[i,2]  <- sample(urn, 1)
  
  # Update Imbalance:
  if(pids2[i,2] == "A"){numA=numA+1}else{numB=numB+1}
  D2[i+1]       <- numA - numB
  
  # Update Urn:
  urn <- c(rep("A",(N/2 - numA)), rep("B", (N/2 - numB)))
  
  # Compare Guess to Actual Allocation:
  if(pids2[i,2] == guess[i]){
    H = H + 1
  }else{
    MT = MT + 1
  }
}
pids2
D2[601]
(M = MT - length(guess[guess=="T"]))
(biasF2 <- (H-M)/2)

# Plot Random Walk
x = seq(0,N,1)
plot(x, D2, type="l", xlab ="N", ylab="Dj",ylim = c(-23, 27), cex.lab = 2, cex.axis = 1.5)
abline(h=0, col = "red")
grid()



#------------------------------------------ Efron's Biased Coin:

# Initialise:
pids3 <- data.frame(PID = 1:N, Treat = NA)
D3    <- c(0)
guess <- c()
H  <-0 # hit
MT <-0 # miss or tie
Mi <-0 # miss
numA = 0
numB = 0

set.seed(10)
for(i in 1:N){
  
  # Simulate Clinicians Guess and Allocation:
  if(D3[i]>0){
    guess[i] <- "B"
    pids3[i,2]  <- sample(trt, 1, prob = c(0.25, 0.75))  # Assign B with higher prob
    
  }else if(D3[i]<0){
    guess[i] <- "A"
    pids3[i,2]  <- sample(trt, 1, prob = c(0.75, 0.25))  # Assign A with higher prob
    
  }else{
    guess[i]    <- "T"    # T for Tie
    pids3[i,2]  <- sample(trt, 1)
    }   
  
  # Update Imbalance:
  if(pids3[i,2] == "A"){numA=numA+1}else{numB=numB+1}
  D3[i+1]       <- numA - numB
  
  # Compare Guess to Actual Allocation:
  if(pids3[i,2] == guess[i]){
    H = H + 1
  }else{
    MT = MT + 1
  }
}
pids3
D3[601]
(M = MT - length(guess[guess=="T"]))
(biasF3 <- (H-M)/2)

# Plot Random Walk
x = seq(0,N,1)
plot(x, D3, type="l", xlab ="N", ylab="Dj", ylim = c(-23, 27), cex.lab = 2, cex.axis = 1.5)
abline(h=0, col = "red")
grid()






#-----------------
# Question 3f:
#-----------------
(t = ((1/40 + 1/38)^{-1})/((2/300)^{-1}))
(phat = (22+20)/(40+38))
(z = (22/40 - 20/38)/sqrt(phat*(1-phat)*(1/40 + 1/38)))
(b=z*sqrt(t))
(theta = 1.96 + 0.845 )
(expb1 <- b + theta*(1-t))
(q = (1.96 - expb1)/sqrt(1-t))
(cpnull <- 1-pnorm(q=q))
(thetahat = b/t)
(qcur = (1.96 - thetahat)/sqrt(1-t))
(cpcur <- 1-pnorm(q=qcur))




#-----------------
# Question 3g:
#-----------------

library(gsDesign)

x1 <- gsDesign(k=3, test.type = 2, alpha = 0.05, beta = 0.2, sfu = "OF", n.fix = 600, n.I = c(80,300,600))
x1
gsBoundSummary(x1)

x2 <- gsDesign(k=3, test.type = 2, alpha = 0.05, beta = 0.2, sfu = "Pocock", n.fix = 600, n.I = c(80,300,600))
x2
gsBoundSummary(x2)


library(AGSDest)
x3 <- plan.GST(K=3, t=c(80,300,600)/600, SF=1, alpha = 0.05, pow=0.8)


