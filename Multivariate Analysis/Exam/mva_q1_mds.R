# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(MASS)


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
# MDS
#------------------
# AIM: 
# Metric MDS aims to obtain a map which approximates the given inter-sample dissimilarities as closely as possible
# So in our case using euclidean distance, we want to obtain a map to a lower-dimensional space that best represents the eucliden distances between observations  
# This will allow us to identify clusters of observations based on how similar their skull measurements are

# Which MDS to use?
# We would like to use euclidean distance to form our dissimilarity matrix since our variables are all continous and euclidean distance makes the most sense for the problem
# Thus we will not use classical MDS since doing classical MDS with euclidean distances is equivalent to preforming PCA, which we have already done.
# Therefore we chose to use metric MDS, since our data is continuous and we are not interested in scrictly ordering our dissimilarities. 
# We decided to use Sammon mapping instead of metric least-squares scaling because sammon mapping preserves small dissimilarities, making it better for trying to identif clusters \cit{izeman}


d <- dist(dat_scal, method = "euclidean", diag = T, upper = T)
d

set.seed(123)
sam_fit <- sammon(d, k=2, niter = 1000)
# chose 1000 iterations and a tolerance of 0.0004, and reached convergence after only 50 iterations.
# We chose the classical MDS as the initial configuration, since the PCA had been pretty successful.
# Again 2 dimensions were chosen for the configuration


# Before plotting our results, we need to check the fit and whether the map accurately represents the dissimilarities.
sam_fit$stress
# The stress value is very low, only 0.034, suggesting that the map does do a good job of representing the disimilarities in a lower dimension


plot_dat1 = data.frame("dim1" = sam_fit$points[,1],
                       "dim2" = sam_fit$points[,2],
                       "Class" = as.factor(dat$class))

ggplot(plot_dat1, aes(x=dim1, y=dim2, color = Class)) +
  geom_point(size=5) +
  #scale_color_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))+
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_text(size = 25, face = "bold"), 
        legend.text=element_text(size=23), axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) 

# The MDS plot is very similar to the PCA plot that was produced, thus the horizontal axis and the vertical axis will have very similar meanings to the first and second components respectively.
# Again four clusters are formed, one containing Hylobates observations, another containing Homo observations, another containing mostly Pongo observations, with a few Homo observations, and another containing a mixture of Gorilla and Pam observations.
# This indicates that observations in any of the five primate classes generally have skull measurements very similar to other primates in their same class, with regards to the 8 skull measurement variables.
# Generally primates tend to have skull measurements very different to primates from a different class, with the exception of primates from the Gorilla and Pan classes who have skull measurements very similar to each other. 
# So again there is evidence that these 8 skull measurement variables are able classify primates by their class, with the exception of Gorilla and Pan primates which the variables are not good a differentiating. 




