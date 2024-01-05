# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(RColorBrewer)


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
#dat_scal <- log(dat[, -c(1, 2, 3)]+1)
#dat_scal <- scale(dat_scal)
dat_scal <- scale(dat[, -c(1, 2, 3)])
dat_scal <- as.data.frame(dat_scal)
head(dat_scal)
dim(dat_scal)  


#-------------------
# PCR
#-------------------

prcfit <- prcomp(dat_scal,scale=TRUE) # PCA on correlation matrix
prcfit
summary(prcfit)

plot(x=1:25, 
     y=prcfit$sdev^2,        # eigenvalues = variance explained by that component
     xlab = "Components",
     ylab = "Variance",
     pch = 16,
     type = "b",
     main = "" )
grid()


# prcfit$sdev^2 = eigenvalues
# recall that because the data is standardised, the total variation of all 25 
# standardised variables is 25
# so for the first component, prcfit$sdev[1]^2/25 = 2.9638^2/25 gives the proportion
# of variance explained by the first component

# quite hard to see an elbow - could be 5

# if we used the 90% criteria, we'd have to retain 10 components !

prcfit$sdev^2
# if we used the eigenvalues > 1, we'd have to retain 6 components - but the eigenvalue
# condition often underestimates the amount of pc needed



# So really not much point in plotting.....
# could very well be that a linear reduction doesnt exist



#-------------------------
# VISUALISATION
#------------------------

plot_dat1 = data.frame("PC1" = prcfit$x[,1],
                           "PC2" = prcfit$x[,2],
                           "GroupxStatus" = as.factor(interaction(dat$Group,dat$Status)))

ggplot(plot_dat1, aes(x=PC1, y=PC2, color = GroupxStatus)) +
  geom_point(size=5) +
  scale_color_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))+
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_text(size = 25, face = "bold"), 
        legend.text=element_text(size=23), axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) 

round(prcfit$rotation[,1:6], 3)

fviz_pca_biplot(prcfit)
