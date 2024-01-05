# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(Rtsne)
library(ggplot2)

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
# tSNE
#-------------------
set.seed(123)
tsne1 <- Rtsne(dat_scal, 
              dims = 2, 
              intial_dims = ncol(dat_scal),
              perplexity = 5,  # Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
              theta = 0,   # Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
              verbose=TRUE, 
              max_iter = 100000)
plot_dat1 = data.frame("tSNE1" = tsne1$Y[,1],
                        "tSNE2" = tsne1$Y[,2],
                        "Status" = as.factor(dat$Status),
                        "Group" = as.factor(dat$Group))
ggplot(plot_dat1, aes(x=tSNE1, y=tSNE2, color = Status, shape = Group)) +
  geom_point() +
   scale_shape_manual(values = c(0, 8)) +
   scale_color_manual(values = c("blue", "red")) +
   theme(legend.position = "right")


set.seed(123)
tsne2 <- Rtsne(dat_scal, 
               dims = 2, 
               intial_dims = ncol(dat_scal),
               perplexity = 10,  # Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
               theta = 0,      # Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
               verbose=TRUE, 
               max_iter = 100000)
plot_dat2 = data.frame("tSNE1" = tsne2$Y[,1],
                       "tSNE2" = tsne2$Y[,2],
                       "Status" = as.factor(dat$Status),
                       "Group" = as.factor(dat$Group))
ggplot(plot_dat2, aes(x=tSNE1, y=tSNE2, color = Status, shape = Group)) +
  geom_point() +
  scale_shape_manual(values = c(0, 8)) +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "right")



set.seed(123)
tsne3 <- Rtsne(dat_scal, 
               dims = 2, 
               intial_dims = ncol(dat_scal),
               perplexity = 15,  # Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
               theta = 0,      # Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
               verbose=TRUE, 
               max_iter = 100000)
plot_dat3 = data.frame("tSNE1" = tsne3$Y[,1],
                       "tSNE2" = tsne3$Y[,2],
                       "Status" = as.factor(dat$Status),
                       "Group" = as.factor(dat$Group))
ggplot(plot_dat3, aes(x=tSNE1, y=tSNE2, color = Status, shape = Group)) +
  geom_point() +
  scale_shape_manual(values = c(0, 8)) +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "right")




set.seed(123)
tsne4 <- Rtsne(dat_scal, 
               dims = 2, 
               intial_dims = ncol(dat_scal),
               perplexity = 20,  # Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
               theta = 0,   # Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
               verbose=TRUE, 
               max_iter = 100000)
plot_dat4 = data.frame("tSNE1" = tsne4$Y[,1],
                       "tSNE2" = tsne4$Y[,2],
                       "Status" = as.factor(dat$Status),
                       "Group" = as.factor(dat$Group))
ggplot(plot_dat4, aes(x=tSNE1, y=tSNE2, color = Status, shape = Group)) +
  geom_point() +
  scale_shape_manual(values = c(0, 8)) +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "right")

# tsne4 with perplexity=20 looks like the best representation because it has the lowest error 
# co indication of clustering whatsoever 


