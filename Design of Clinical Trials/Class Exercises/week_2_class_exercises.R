# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#--------------------
# Load Libraries
#--------------------
library(iterpc)

#-----------------------------------------
# One Permuted Block Design Allocation
#-----------------------------------------
# Permuted block design with 4 blocks of size 6
# Two treatments: A and B
# Doing one allocation before simulating 10 000
# Use random allocation rule

# Need to generate all the possible combinations of 3 A's and 3 B's that we can have in each block
# then from these combinations, we randomly assign one (with replacement) to each block

# Possible Combinations:
x <- iterpc(table(rep(c("A", "B"), 3)), ordered = T, labels = c("A", "B"))
trt_combos <- getall(x)

# Now implement the algorithm underlying the randomization design
allocate_pbd <- function(m, M){     # PBD allocation function for treatments A and B
  
  # Initialize:
  D <- NULL   # recording imbalance
  D[1]<-0
  # Initially we have zero people allocated to group A and B
  numberA <- 0   
  numberB <- 0
  # Keep track of numbers of As and Bs for each round
  storeA <- NULL
  storeB <- NULL
  # Recording patient allocation (one for each block)
  allocation <- matrix(NA, nrow = m, ncol = M)  # each column is a different block
  
  # loop over blocks:
  for(i in 1:M){       
    
    # Choose a set of treatments from trt_combos
    row <-  sample(1:nrow(trt_combos), 1)
    allocation[,i] <- trt_combos[row,]
    
    # loop within blocks:
    for (j in 1:m) { 
      
      if(allocation[j,i] == "A"){             
        numberA = numberA + 1           
      }else{numberB = numberB + 1}
      
      # Store numbers
      storeA[j + m*(i-1)] <- numberA
      storeB[j + m*(i-1)] <- numberB
      
      # Record Imbalance
      D[j + m*(i-1) + 1]<- numberA-numberB    
      
    }
    
    
  }
  
  return(list(allocation = allocation, D= D, numberA = storeA, numberB = storeB))
}


set.seed(123)
(pbd <- allocate_pbd(m=6, M=4))

# Plot random walk
patients <- seq(0,24,1)
plot(patients, pbd$D, type="l", xaxt="n", ylim=c(-3,3), ylab="Dj")
axis(1,0:24)
abline(h=0,col='red')


#-----------------------------------------
# One Random Block Design
#-----------------------------------------
# Random block design using blocks of size 2, 4 or 6
# i.e. 12 blocks of size 2, 6 blocks of size 4, 4 blocks of size 6
# Two treatments: A and B
# Doing one allocation before simulating 10 000
# Use random allocation rule


# Using random allocation rule, we need an urn with m/2 A's and m/2 B's 
# and that we sample WITHOUT replacement for EACH BLOCK


# Implement the algorithm underlying the randomization design
allocate_rcbd <- function(m, M){     # RCBD allocation function for treatments A and B
  
  # Initialize:
  D <- NULL   # recording imbalance
  D[1]<-0
  # Initially we have zero people allocated to group A and B
  numberA <- 0   
  numberB <- 0
  # Make random allocation urn
  urn <- rep(c("A", "B"), m/2)
  # Keep track of numbers of As and Bs for each round
  storeA <- NULL
  storeB <- NULL
  # Recording patient allocation (one for each block)
  allocation <- matrix(NA, nrow = m, ncol = M)  # each column is a different block
  
  # loop over blocks:
  for(i in 1:M){       
    
    # sample without replacement for Block i
    allocation[,i] <- sample(urn, m, replace = FALSE) 
    
    # loop within blocks:
    for (j in 1:m) { 
      
      if(allocation[j,i] == "A"){             
        numberA = numberA + 1           
      }else{numberB = numberB + 1}
      
      # Store numbers
      storeA[j + m*(i-1)] <- numberA
      storeB[j + m*(i-1)] <- numberB
      
      # Record Imbalance
      D[j + m*(i-1) + 1]<- numberA-numberB    
      
    }
    
    
  }
  
  return(list(allocation = allocation, D= D, numberA = storeA, numberB = storeB))
}


set.seed(123)
# 12 blocks of size 2, 6 blocks of size 4, 4 blocks of size 6
rcbd12 <- allocate_rcbd(m=2, M=12)
rcbd6 <- allocate_rcbd(m=4, M=6)
rcbd4 <- allocate_rcbd(m=6, M=4)


# Plot random walk
patients <- seq(0,24,1)
plot(patients, rcbd12$D, type="l", xaxt="n", ylim=c(-3,3), ylab="Dj")
axis(1,0:24)
abline(h=0,col='red')

plot(patients, rcbd6$D, type="l", xaxt="n", ylim=c(-3,3), ylab="Dj")
axis(1,0:24)
abline(h=0,col='red')

plot(patients, rcbd4$D, type="l", xaxt="n", ylim=c(-3,3), ylab="Dj")
axis(1,0:24)
abline(h=0,col='red')

# Notice how as the size of the blocks increase, so does the maximum possible imbalance


#-----------------------------------------------------------
# 10000 Simulations of Permuted Block Design Allocation
#-----------------------------------------------------------

# Initialise:
exp_pbd    <- NULL # record expectation of T => the proportion of A  
biasF_pbd  <- NULL # record bias factor
D_pbd      <- NULL # Keep track of imbalance during trial
M = 4
m = 6

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = 10000, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
) 
for(k in 1:10000){
  
  # Initialize:
  
  # Keeping track of the clinicians guessing for selection bias
  guess <- NULL
  H  <-0 # hit
  MT <-0 # miss or tie
  Mi <-0 # miss
  expec <- NULL
  
  # ACTUALL ALLOCATION
  actalo <- allocate_pbd(m,M)
  
  
  for(j in 1:(m*M)){
    
    # GUESS 
    if(actalo$numberA[j] > actalo$numberB[j]){
      guess[j] <- "B"
    }else if(actalo$numberB[j] > actalo$numberA[j]){
      guess[j] <- "A"
    }else{guess[j]<-"T"}   # T for Tie
    
    # Update hit and miss scores
    if(j%%m == 0){t=m}else{t=j%%m }
    if(actalo$allocation[t, ceiling(j/m)]==guess[j]){
      H = H + 1
    }else{
      MT = MT + 1}
    
    # Update Expectation
    expec[j] <- actalo$numberA[j]/(actalo$numberA[j] + actalo$numberB[j])
    
  }
  
  Mi           <- MT - length(guess[guess=="T"])
  biasF_pbd[k] <- (H-Mi)/2
  D_pbd        <- append(D_pbd, actalo$D)
  exp_pbd      <- append(exp_pbd, expec)
  
  # Add progress bar because slow :(
  setTxtProgressBar(pb, k)
  
}


mean(biasF_pbd)
mean(exp_pbd)


# Plot
barplot(table(D_pbd),main="Distribution of the Imbalance",xlab="D",ylab="Frequency")
hist(exp_pbd)
hist(biasF_pbd)

# Plotting probabilities of imbalance.
probabilities <- NULL
for(i in 1:7){
  val <- -8+(i*2)
  probabilities[i] <- (length(D_pbd[D_pbd==val]))/length(D_pbd)
}
r <- seq(-6,6,2)
plot(r,probabilities,type="l",main="Probability of Imbalance r",ylab="Probability")



#-----------------------------------------------------------
# 10000 Simulations of RCBD Allocation: 12 blocks of size 2
#-----------------------------------------------------------

# Initialise:
exp_rcbd   <- NULL # record expectation of T => the proportion of A  
biasF_rcbd <- NULL # record bias factor
D_rcbd     <- NULL
M = 12
m = 2

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = 10000,  # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
) 
for(k in 1:10000){
  
  # Initialize:
  
  # Keeping track of the clinicians guessing for selection bias
  guess <- NULL
  H  <-0 # hit
  MT <-0 # miss or tie
  Mi  <-0 # miss
  expec <- NULL
  
  # ACTUALL ALLOCATION
  actalo <- allocate_rcbd(m,M)
  
  for(j in 1:(m*M)){
    
    # GUESS 
    if(actalo$numberA[j] > actalo$numberB[j]){
      guess[j]<-"B"
    }else if(actalo$numberB[j] > actalo$numberA[j]){
      guess[j]<-"A"
    }else{guess[j]<-"T"}   # T for Tie
    
    # Update hit and miss scores
    if(actalo$allocation[(j%%m + 1), ceiling(j/m)]==guess[j]){
      H = H + 1
    }else{
      MT = MT + 1}
    
    # Update Expectation
    expec[j] <- actalo$numberA[j]/(actalo$numberA[j] + actalo$numberB[j])
    
    
  }
  
  Mi        <- MT - length(guess[guess=="T"])
  biasF_rcbd[k] <- (H-Mi)/2
  D_rcbd        <- append(D_rcbd, actalo$D)
  exp_rcbd      <- append(exp_rcbd, expec)
  
  # Add progress bar because slow :(
  setTxtProgressBar(pb, k)
  
}

mean(biasF_rcbd)
mean(exp_rcbd)


# Plot
barplot(table(D_rcbd),main="Distribution of the Imbalance",xlab="D",ylab="Frequency")
hist(exp_rcbd)
hist(biasF_rcbd)


# Plotting probabilities of imbalance.
probabilities <- NULL
for(i in 1:7){
  val <- -8+(i*2)
  probabilities[i] <- (length(D_rcbd[D_rcbd==val]))/length(D_rcbd)
}
r <- seq(-6,6,2)
plot(r,probabilities,type="l",main="Probability of Imbalance r",ylab="Probability")



#-----------------------------------------------------------
# 10000 Simulations of RCBD Allocation: 6 blocks of size 4
#-----------------------------------------------------------

# Initialise:
exp_rcbd   <- NULL # record expectation of T => the proportion of A  
biasF_rcbd <- NULL # record bias factor
D_rcbd     <- NULL
M = 6
m = 4

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = 10000, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
) 
for(k in 1:10000){
  
  # Initialize:
  
  # Keeping track of the clinicians guessing for selection bias
  guess <- NULL
  H  <-0 # hit
  MT <-0 # miss or tie
  Mi  <-0 # miss
  expec <- NULL
  
  # ACTUALL ALLOCATION
  actalo <- allocate_rcbd(m,M)
  
  for(j in 1:(m*M)){
    
    # GUESS 
    if(actalo$numberA[j] > actalo$numberB[j]){
      guess[j]<-"B"
    }else if(actalo$numberB[j] > actalo$numberA[j]){
      guess[j]<-"A"
    }else{guess[j]<-"T"}   # T for Tie
    
    # Update hit and miss scores
    if(j%%m == 0){t=m}else{t=j%%m }
    if(actalo$allocation[t, ceiling(j/m)]==guess[j]){
      H = H + 1
    }else{
      MT = MT + 1}
    
    # Update Expectation
    expec[j] <- actalo$numberA[j]/(actalo$numberA[j] + actalo$numberB[j])
    
    
  }
  
  Mi        <- MT - length(guess[guess=="T"])
  biasF_rcbd[k] <- (H-Mi)/2
  D_rcbd        <- append(D_rcbd, actalo$D)
  exp_rcbd      <- append(exp_rcbd, expec)
  
  # Add progress bar because slow :(
  setTxtProgressBar(pb, k)
  
}

mean(biasF_rcbd)
mean(exp_rcbd)


# Plot
barplot(table(D_rcbd),main="Distribution of the Imbalance",xlab="D",ylab="Frequency")
hist(exp_rcbd)
hist(biasF_rcbd)


# Plotting probabilities of imbalance.
probabilities <- NULL
for(i in 1:7){
  val <- -8+(i*2)
  probabilities[i] <- (length(D_rcbd[D_rcbd==val]))/length(D_rcbd)
}
r <- seq(-6,6,2)
plot(r,probabilities,type="l",main="Probability of Imbalance r",ylab="Probability")


#-----------------------------------------------------------
# 10000 Simulations of RCBD Allocation: 4 blocks of size 6
#-----------------------------------------------------------

# Initialise:
exp_rcbd   <- NULL # record expectation of T => the proportion of A  
biasF_rcbd <- NULL # record bias factor
D_rcbd     <- NULL
M = 4
m = 6

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = 10000, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
) 
for(k in 1:10000){
  
  # Initialize:
  
  # Keeping track of the clinicians guessing for selection bias
  guess <- NULL
  H  <-0 # hit
  MT <-0 # miss or tie
  Mi  <-0 # miss
  expec <- NULL
  
  # ACTUALL ALLOCATION
  actalo <- allocate_rcbd(m,M)
  
  for(j in 1:(m*M)){
    
    # GUESS 
    if(actalo$numberA[j] > actalo$numberB[j]){
      guess[j]<-"B"
    }else if(actalo$numberB[j] > actalo$numberA[j]){
      guess[j]<-"A"
    }else{guess[j]<-"T"}   # T for Tie
    
    # Update hit and miss scores
    if(j%%m == 0){t=m}else{t=j%%m }
    if(actalo$allocation[t, ceiling(j/m)]==guess[j]){
      H = H + 1
    }else{
      MT = MT + 1}
    
    # Update Expectation
    expec[j] <- actalo$numberA[j]/(actalo$numberA[j] + actalo$numberB[j])
    
    
  }
  
  Mi        <- MT - length(guess[guess=="T"])
  biasF_rcbd[k] <- (H-Mi)/2
  D_rcbd        <- append(D_rcbd, actalo$D)
  exp_rcbd      <- append(exp_rcbd, expec)
  
  # Add progress bar because slow :(
  setTxtProgressBar(pb, k)
  
}

mean(biasF_rcbd)
mean(exp_rcbd)


# Plot
barplot(table(D_rcbd),main="Distribution of the Imbalance",xlab="D",ylab="Frequency")
hist(exp_rcbd)
hist(biasF_rcbd)


# Plotting probabilities of imbalance.
probabilities <- NULL
for(i in 1:7){
  val <- -8+(i*2)
  probabilities[i] <- (length(D_rcbd[D_rcbd==val]))/length(D_rcbd)
}
r <- seq(-6,6,2)
plot(r,probabilities,type="l",main="Probability of Imbalance r",ylab="Probability")






#-----------------------
# Use RandomiseR
#----------------------

devtools::install_github("cran/randomizeR")
library(randomizeR)

# PBD
pbd <- pbrPar(c(6,6,6,6), 2, groups = c("A", "B"))
getRandList(genSeq(pbd))

pbd_seq <- getAllSeq(pbd)
imbalance <- assess(pbd_seq, selBias("DS", method="sim"), # how do we get eta ?!
                    imbal("maxImb")) # no point in doing final imbalance, so do max
imbalance
