# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(sparsepca)


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
# SPCA
#-------------------
set.seed(123)
spca_fit <- sparsepca::spca(dat_scal,   
                 k = 5,           # using scree plot as reference
                 alpha = 1e-2,  #Sparsity controlling parameter. Higher values lead to sparser components
                 beta = 1e-4) # Amount of ridge shrinkage to apply in order to improve conditioning
# keep convergence tolerance and max iterations as default
# data already  centred and scaled so dont need to specify that
# beta was keep at default 
# alpha was adjusted to get desired sparsity without compromising too much on the percentage of variance explained
spca_fit
round(summary(spca_fit), 3)

length(which(spca_fit$loadings[,1]!=0))
length(which(spca_fit$loadings[,2]!=0))
length(which(spca_fit$loadings[,3]!=0))
length(which(spca_fit$loadings[,4]!=0))
length(which(spca_fit$loadings[,5]!=0))

plot_dat1 = data.frame("PC1" = spca_fit$scores[,1],
                       "PC2" = spca_fit$scores[,2],
                       "GroupxStatus" = as.factor(interaction(dat$Group,dat$Status)))

ggplot(plot_dat1, aes(x=PC1, y=PC2, color = GroupxStatus)) +
  geom_point(size=5) +
  scale_color_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))+
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_text(size = 25, face = "bold"), 
        legend.text=element_text(size=23), axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) 

round(spca_fit$loadings, 3)



# ROBUST SPCA:
set.seed(123)
robspca_fit <- sparsepca::robspca(dat_scal,
                       k=5,
                       alpha = 1e-2, 
                       beta = 1e-04, 
                       gamma = 2)  # Sparsity controlling parameter for the error matrix S. 
                                     # Smaller values lead to a larger amount of noise removal.

robspca_fit
summary(robspca_fit)

which(robspca_fit$sparse!=0,arr.ind = T) 

length(which(robspca_fit$loadings[,1]!=0))
length(which(robspca_fit$loadings[,2]!=0))
length(which(robspca_fit$loadings[,3]!=0))
length(which(robspca_fit$loadings[,4]!=0))
length(which(robspca_fit$loadings[,5]!=0))


plot_dat2 = data.frame("PC1" = robspca_fit$scores[,1],
                       "PC2" = robspca_fit$scores[,2],
                       "GroupxStatus" = as.factor(interaction(dat$Group,dat$Status)))

ggplot(plot_dat2, aes(x=PC1, y=PC2, color = GroupxStatus)) +
  geom_point(size=5) +
  scale_color_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))+
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_text(size = 25, face = "bold"), 
        legend.text=element_text(size=23), axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) 

round(robspca_fit$loadings, 3)




