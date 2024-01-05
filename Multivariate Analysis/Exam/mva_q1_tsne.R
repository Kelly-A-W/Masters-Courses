# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(Rtsne)


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
# t-SNE
#-------------------
# AIM:
# The aim of t-SNE is to map high dimensional data that lies on or near a low-dimensional, non-linear manifold to this low-dimensional space while capturing \emph{both} the local and global structure, including the presence of clusters \cite{izeman; laurens van der Maaten, Geoffrey Hinton, Visualizing data using t-SNE, Journal of machine learning research, vol 9, 2008, 2579-2605}
# we aim to use this to represent the data as best as possible using a \emph{non-linear} transformation, in order to identify clusters in order to investigate if the skull measurement variables can differeniate between primate classes
# t-SNE has the advantage of kernel PCA of being able to capture both global and local structure


# The main tuning parameter for t-SNE is perplexity, which cannot be bigger than one third of the number of number of observations in the data minus 1
# Thus for our example the perplexity parameter must be less than (105-1)/3 = 34.67 
# cite t-SNE web resource for the following:
# The perplexity parameter essentially balances the importance of global structure and local structure
# It can be thought of roughly as the assumed amount of close neighbours that every point has.
# To choose the best perplexity value, we fitted t-SNE and evaluated the plots using perplexity values of 0 to 30, in increments of 5. 
#  Another parameter is the speed to accuracy trade-off parameter, which was set to 0 since the data set was relatively small. Increasing this parameter would result in a faster algorithm but a decrease in accuracy.
set.seed(123)
tsne1 <- Rtsne(dat_scal, 
               dims = 2, 
               intial_dims = ncol(dat_scal),
               perplexity = 5,  
               theta = 0,   
               verbose=TRUE, 
               max_iter = 100000)
plot_dat1 = data.frame("tSNE1" = tsne1$Y[,1],
                       "tSNE2" = tsne1$Y[,2],
                       "Class" = as.factor(dat$class))
ggplot(plot_dat1, aes(x=tSNE1, y=tSNE2, color = Class)) +
  geom_point() +
  theme(legend.position = "right")




set.seed(123)
tsne2 <- Rtsne(dat_scal, 
               dims = 2, 
               intial_dims = ncol(dat_scal),
               perplexity = 10,  
               theta = 0,   
               verbose=TRUE, 
               max_iter = 100000)
plot_dat2 = data.frame("tSNE1" = tsne2$Y[,1],
                       "tSNE2" = tsne2$Y[,2],
                       "Class" = as.factor(dat$class))
ggplot(plot_dat2, aes(x=tSNE1, y=tSNE2, color = Class)) +
  geom_point() +
  theme(legend.position = "right")


set.seed(123)
tsne3 <- Rtsne(dat_scal, 
               dims = 2, 
               intial_dims = ncol(dat_scal),
               perplexity = 15,  
               theta = 0,   
               verbose=TRUE, 
               max_iter = 100000)
plot_dat3 = data.frame("tSNE1" = tsne3$Y[,1],
                       "tSNE2" = tsne3$Y[,2],
                       "Class" = as.factor(dat$class))
ggplot(plot_dat3, aes(x=tSNE1, y=tSNE2, color = Class)) +
  geom_point() +
  theme(legend.position = "right")




set.seed(123)
tsne4 <- Rtsne(dat_scal, 
               dims = 2, 
               intial_dims = ncol(dat_scal),
               perplexity = 20,  
               theta = 0,   
               verbose=TRUE, 
               max_iter = 100000)
plot_dat4 = data.frame("tSNE1" = tsne4$Y[,1],
                       "tSNE2" = tsne4$Y[,2],
                       "Class" = as.factor(dat$class))
ggplot(plot_dat4, aes(x=tSNE1, y=tSNE2, color = Class)) +
  geom_point() +
  theme(legend.position = "right")




set.seed(123)
tsne5 <- Rtsne(dat_scal, 
               dims = 2, 
               intial_dims = ncol(dat_scal),
               perplexity = 25,  
               theta = 0,   
               verbose=TRUE, 
               max_iter = 100000)
plot_dat5 = data.frame("tSNE1" = tsne5$Y[,1],
                       "tSNE2" = tsne5$Y[,2],
                       "Class" = as.factor(dat$class))
ggplot(plot_dat5, aes(x=tSNE1, y=tSNE2, color = Class)) +
  geom_point() +
  theme(legend.position = "right")



set.seed(123)
tsne6 <- Rtsne(dat_scal, 
               dims = 2, 
               intial_dims = ncol(dat_scal),
               perplexity = 30,  
               theta = 0,   
               verbose=TRUE, 
               max_iter = 100000)
plot_dat6 = data.frame("tSNE1" = tsne6$Y[,1],
                       "tSNE2" = tsne6$Y[,2],
                       "Class" = as.factor(dat$class))
ggplot(plot_dat6, aes(x=tSNE1, y=tSNE2, color = Class)) +
  geom_point() +
  theme(legend.position = "right")



# Using a perplexity parameter of 15 resulted in the best separation of the clusters. 
# Looking at Figure, we can see that all the classes are fairly well separated, even the Pan and Gorilla classes.
# However the Pan and Gorilla classes are still in the same cluster, and even though the Pan observations are all on the right (except one) and the Gorilla observations on the left,  without knowing the classes of the observations it would be impossible to tell them apart because they are all well contained in the same cluster. 
# Thus there is again evidence that the skull measurements cannot be used to differeniate between the Pan and Gorilla primates.
# However all of the other clusters are well separated, suggesting again that the skull measurement variables can differeniate between the other classes.
# Note that when using t-SNE, we can't infer anything about the distances between clusters: so for instance the Homo cluster being close to the Pongo cluster does not mean that the observations are very similar to each other. 





