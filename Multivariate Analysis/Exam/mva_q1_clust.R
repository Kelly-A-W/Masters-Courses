# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(NbClust)
library(cluster)


#-------------------
# Load Data
#-------------------

dat <- read.csv("PrimatesData.csv", check.names = FALSE)
head(dat)
str(dat)

#-------------------
# Data Wrangling
#-------------------

# Set Factor Variables
dat$genus <- as.factor(dat$genus)
dat$class <- as.factor(dat$class)
dat$classdigit <- as.factor(dat$classdigit)
# and GENUSID is just a character variable (like pid)

str(dat) # as desired

# scale continuous variables:
dat_scal <- scale(dat[,3:10])
dat_scal <- as.data.frame(dat_scal)
str(dat_scal) # looks good !

#-------------------
# PAM Cluster Analysis
#-------------------
# AIM:
# The aim of using PAM cluster analysis is to see if the the 8 skull measurement variables contain enough information on their own to differentiate observations (primates) according to their class.
# So, by performing PAM clustering on the 8 skull measurement variables, we hope to see the formation of five different clusters, each corresponding to a different primate class. 

# Although we know that we would like to see 5 clusters formed, we can look at the Gap statistic to see how many clusters are suggested, especially since there is evidence that the variables are not very good at diffrientiating between Gorilla and Pan classes.

set.seed(123)
fviz_nbclust(dat_scal, pam, method="gap_stat" )+theme_classic()
set.seed(123)
fviz_nbclust(dat_scal, pam, method="silhouette")+theme_classic()
set.seed(123)
fviz_nbclust(dat_scal, pam, method="wss")+theme_classic()
# From Figure 1 we can see that the suggested number of clusters is 4, while Figure 2 using the average silhouette width suggests only 2 clusters. 
# Figure 3 shows that the within cluster sum of squares suggests 4 or 5 clusters. 
# Since we would like to see how well the variables are able to differeniate between all 5 classes, we decided to choose 5 clusters. 

set.seed(123)
pam5 <- pam(dat_scal, 5)
summary(pam5)
plot(pam5, main ="")
fviz_cluster(pam5)
# average silhouette width is only 0.34, suggesting that the structure is weak and could be artificial 
# From the plot we can see that in total 76.02\% of the point variability is explained, which is not too bad. 
# The two clusters, clusters one and four, on the left of the plot are well separated, but the remaining three clusters are not, especially clusters two and three at the top left of the plot.
# Comparing these observations to the silhouette plot, we can see that cluster 4 has the highest average silhouette width of 0.45, while clusters three and five have very low silhoutte widths


nhc_dat <- data.frame("Class" = dat[,11], 
                      "Cluster" = as.factor(pam5$clustering))
table(nhc_dat)
# Pam has successfully clustered all the Hylobates observations together in cluster 1 and all of the Pongo observations together in cluster 2.
# However it is clear than PAM could not differentiate between the Gorilla observations and the Pan observations, clustering them both into cluster 4.
# Instead of separating the Gorilla and Pan observations into separate clusters, PAM choices to rather split the Homo observations into clusters 3 and 5 to meet the five cluster requirement.
# This provides further strong evidence that the 8 skull measurement variables do not contain enough information to differeniate between Pan and Gorilla primates. 
# However there is also evidence that the remaining groups can be classified using the 8 skull variables. 




