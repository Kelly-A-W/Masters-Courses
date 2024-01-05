# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggplot2)
library(psych)
library(GPArotation)
library(fastICA)
library(ggcorrplot)
library(ggpubr)

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
str(dat)  
head(dat)

# Remove Factor Variables and Scale
dat_scal <- scale(dat[, -c(1, 2, 3)])
dat_scal <- as.data.frame(dat_scal)
head(dat_scal)
dim(dat_scal)   


#-------------------
# Factor Analysis
#-------------------
# CANNOT use mle estimation because haven't checked data normality
# so cant use factanal function
# don't want to use PCA method because our PCAs were useless

# CORRELATIONS
corr <- cor(dat_scal)
ggcorrplot(corr, type = "lower",
           outline.color = "white",
           hc.order = TRUE)
p.mat <- cor_pmat(dat_scal)
ggcorrplot(corr, type = "lower",
           outline.color = "white",
           p.mat = p.mat,
           insig = "blank",
           hc.order = TRUE) # ordering variables using hclust function
# hierarchical clustering 


# how many factors:
set.seed(123)
fa.parallel(dat_scal, 
            fm = "wls", 
            fa = "fa", 
            main = "")

# perform EFA
set.seed(123)
fa_fit <- fa(dat_scal, 
             nfactors = 4,
             rotate = "oblimin", # this and varimax are most popular
             fm="wls")
summary(fa_fit)
print(fa_fit$loadings,cutoff = 0.3)
fa_fit$EBIC

# lets compare EFA using different number of factors
set.seed(123)
fa_fit2 <- fa(dat_scal, 
             nfactors = 5,
             rotate = "oblimin", # this and varimax are most popular
             fm="wls")
summary(fa_fit2)
print(fa_fit2$loadings,cutoff = 0.3)
fa_fit2$EBIC

set.seed(123)
fa_fit3 <- fa(dat_scal, 
              nfactors = 3,
              rotate = "oblimin", # this and varimax are most popular
              fm="wls")
summary(fa_fit3)
print(fa_fit3$loadings,cutoff = 0.3)
fa_fit3$EBIC

# looking at TLI, m=4 produces the largest value
# so choose that model

fa.diagram(fa_fit)

# measures of fit for non-normal data
round(fa_fit$fit, 3)
print(fa_fit$loadings, cutoff=0.3)

 
 
plot_dat = data.frame("F1" = fa_fit$scores[,1],
                       "F2" = fa_fit$scores[,2],
                       "GroupxStatus" = as.factor(interaction(dat$Group,dat$Status)))

ggplot(plot_dat, aes(x=F1, y=F2, color = GroupxStatus)) +
  geom_point(size=5) +
  scale_color_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))+
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_text(size = 25, face = "bold"), 
        legend.text=element_text(size=23), axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) 

 
 
