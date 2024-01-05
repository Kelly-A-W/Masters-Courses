# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(kernlab)

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
# Kernel PCA
#-------------------

# Gaussian radial basis function
kpc1 <- kpca(~., 
            data=dat_scal, 
            kernel="rbfdot",  
            kpar = list(sigma=1),
            features=2)   # want to visualise in 2D, thus want 2 PCs returned
kpc1
plot(rotated(kpc1), col=as.numeric(dat$class))            



#  Laplace radial basis kernel
kpc2 <- kpca(~., 
             data=dat_scal, 
             kernel="laplacedot",  
             kpar = list(sigma=1),
             features=2)   # want to visualise in 2D, thus want 2 PCs returned
plot(rotated(kpc2), col=as.numeric(dat$class))  


plot_dat1 = data.frame("PC1" = rotated(kpc2)[,1],
                       "PC2" = rotated(kpc2)[,2],
                       "Class" = as.factor(dat$class))

ggplot(plot_dat1, aes(x=PC1, y=PC2, color = Class)) +
  geom_point(size=5) +
  #scale_color_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))+
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_text(size = 25, face = "bold"), 
        legend.text=element_text(size=23), axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) 



