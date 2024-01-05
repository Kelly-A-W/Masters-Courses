
library(lme4)
library(simr)

df <- read.csv("LDAdat.csv")
head(df)


df$group <- as.factor(df$group)
df$pid <- as.factor(df$pid)

mod1 <- glmer(score ~ time*group + (1|pid), data = df, family = "poisson", 
              control = glmerControl(optimizer ="Nelder_Mead"))
summary(mod1)
length(unique(df$pid))/nrow(df) # on average less than one observation per pid !
# so the reason why our boundary has a singular fit and we have no values for mean and var of the 
# random effect is because we dont have enough observations per pid !

#------------------------------------------------------- (a)

set.seed(123)
powerSim(mod1, test = "time")
powerSim(mod1, test = "group")

#------------------------------------------------------- (b)

n = 200
(p1 <- length(which(df$group==1))/nrow(df)) # proportion allocated to group 1
grp <- c(rep(1,round(n*p1)), rep(0,round(n*(1-p1))))
grp <- rep(grp, 5)  # because 5 time points (adding one more)
pid <- rep(1:n, 5)
time <- c(rep(0.75, n), rep(1, n), rep(1.5, n), rep(2, n), rep(2.5, n))
dfsim <- data.frame(pid=pid, group=grp, time=time)
head(dfsim)


# Generate outcome variable:
fixed <- c(2.237043, 0.564608, -0.007293, -0.027194)
var(df[which(df$time==0.75),2]) # 3.294678
rand <- c(1.5)
res <- 2

# Create Model
mod2 <- makeLmer(score ~  time*group + (1|pid), fixef = fixed, VarCorr = rand, sigma = res, data=dfsim)
mod2

set.seed(12)
powerSim(mod2, nsim=100, test = fcompare(score~time))
powerSim(mod2, nsim=100, test = fcompare(score~group))

#------------------------------------------------------- (c)
mod3 <- mod2
fixef(mod3)["group"] <- 0.1
set.seed(123)
powerSim(mod3, nsim=100, test = fcompare(score~group))

mod3 <- mod2
fixef(mod3)["group"] <- 0.05
set.seed(123)
powerSim(mod3, nsim=100, test = fcompare(score~group))


#------------------------------------------------------- (d)


mod4 <- extend(mod2, along = "pid", n=300)

set.seed(123)
powerSim(mod4, nsim=100, test = fcompare(score~group))
summary(mod4)


