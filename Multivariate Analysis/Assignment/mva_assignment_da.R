# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(MASS)
library(lattice)

#-------------------
# Load Data
#-------------------

dat <- read.csv("mva_assignment_final_data.csv", check.names = FALSE)
head(dat)
str(dat)

#-------------------
# Data Wrangling
#-------------------

# Set Factor Variables
dat$Group <- as.factor(dat$Group) # level_1 = 5W, level_2 = 9W
dat$Status <- as.factor(dat$Status)  # level_1 = Unvaccinated, level_2 = Vaccinated
str(dat)  # 64 observations and 
head(dat)

# Remove Factor Variables and Scale
dat_scal <- scale(dat[, -c(1, 2, 3)])
dat_scal <- as.data.frame(dat_scal)
head(dat_scal)
dim(dat_scal)  


#-------------------
# Discriminant Analysis
#-------------------
mydat <- dat_scal
mydat$Status <- as.factor(dat$Status)
fit1 <- lda(Status ~. , data = mydat)
fit1
fit1.hat <- predict(fit1)
fit1.hat$class
tab1 <- table(dat$Status, fit1.hat$class)
conCV1 <- rbind(tab1[1,]/sum(tab1[1,]), tab1[2,]/sum(tab1[2,]))
dimnames(conCV1) <- list("Actual" = c("Unvaccinated", "Vaccinated"), "Predicted (cv)" = c("Unvaccinated", "Vaccinated"))
print(round(conCV1,3))


fit2 <- lda(Status ~. , data = mydat, CV=T)
fit2
tab2 <- table(mydat$Status, fit2$class)
conCV2 <- rbind(tab2[1,]/sum(tab2[1,]), tab2[2,]/sum(tab2[2,]))
dimnames(conCV2) <- list("Actual" = c("Unvaccinated", "Vaccinated"), "Predicted (cv)" = c("Unvaccinated", "Vaccinated"))
print(round(conCV2,3))

(E_AER <- (tab2[1,2] + tab2[2,1])/sum(tab2)) 


# use qda so when covariance matrices between the two groups are unequal
#fit3 <- qda(Status ~. , data = mydat, CV=T)   #some group is too small for 'qda'

#The test in MASS:::qda.default is if (any(counts < p + 1)) stop("some group is too small for 'qda'")
#where counts is the number of occurrences in each category and p is the number of columns in the 
#predictor matrix .

# I have very few observations, so its very likely that for each group I have less observations 
# than variables :(
df<-mydat
colnames(df) <- c("Status", "v1", "v2", "v3", "v4", "v5",
                 "v6", "v7", "v8", "v9", "v10",
                 "v11", "v12", "v13", "v14", "v15",
                 "v16", "v17", "v18", "v19", "v20",
                 "v21", "v22", "v23", "v24", "v25")
fit3 <- rda(Status ~. , data = df, crossval = T)
fit3
pred <- fit3 %>% predict(fit3)
tab3 <- table(mydat$Status, fit3$class)
conCV3 <- rbind(tab3[1,]/sum(tab3[1,]), tab3[2,]/sum(tab3[2,]))
dimnames(conCV3) <- list("Actual" = c("Unvaccinated", "Vaccinated"), "Predicted (cv)" = c("Unvaccinated", "Vaccinated"))
print(round(conCV3,3))

(E_AER3 <- (tab3[1,2] + tab3[2,1])/sum(tab3)) 


