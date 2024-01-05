# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(RColorBrewer)
library(factoextra)


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
# PCA
#-------------------

prcfit <- prcomp(dat_scal,scale=TRUE) # PCA on correlation matrix
prcfit


# HOW MANY COMPONENTS:

summary(prcfit)   


plot(x=1:8, 
     y=prcfit$sdev^2,        # eigenvalues = variance explained by that component
     xlab = "Components",
     ylab = "Variance",
     pch = 16,
     type = "b",
     main = "" )
grid()
# can't really tell, maybe 4 or 5

plot(prcfit)

prcfit$sdev^2  # eigenvalues


#-------------------------
# VISUALISATION
#------------------------

plot_dat1 = data.frame("PC1" = prcfit$x[,1],
                       "PC2" = prcfit$x[,2],
                       "Class" = as.factor(dat$class))

ggplot(plot_dat1, aes(x=PC1, y=PC2, color = Class)) +
  geom_point(size=5) +
  #scale_color_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))+
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_text(size = 25, face = "bold"), 
        legend.text=element_text(size=23), axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) 

round(prcfit$rotation[,1:4], 3)

fviz_pca_biplot(prcfit, habillage=dat$class)



