
# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load packages
library(ggcorrplot)
library(ggpubr)
library(mclust)
library(MASS)
library(cluster)
library(fpc)
library(dbscan)
library(NbClust)
library(factoextra)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(RColorBrewer)
library(biclust)

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
str(dat)  # 64 observations and 28 variables 
head(dat)

# Remove Factor Variables and Scale
dat_scal <- scale(dat[, -c(1, 2, 3)])
dat_scal <- as.data.frame(dat_scal)
head(dat_scal)
dim(dat_scal)   # 25 continuous variables


#-------------------
# EDA: CONTINOUS VARIABLES
#-------------------

# CORRELATIONS
corr <- cor(dat_scal)
ggcorrplot(corr, type = "lower",
                       outline.color = "white")
p.mat <- cor_pmat(dat_scal)
ggcorrplot(corr, type = "lower",
                       outline.color = "white",
                       p.mat = p.mat,
                       insig = "blank",
                       hc.order = TRUE) # ordering variables using hclust function
                                        # hierarchical clustering 


# NORMALITY
# If the sample size is large enough (n > 30), we can ignore the distribution of 
# the data and use parametric tests.
# The central limit theorem tells us that no matter what distribution things have,
# the sampling distribution tends to be normal if the sample is large enough (n > 30).
# http://www.sthda.com/english/wiki/normality-test-in-r
# We have n>30 observations, so should be fine ! 
# can't really check anyway because too many variables !!!





#-------------------
# EDA: DISCRETE VARIABLES
#-------------------
# We are primarily interested in two groups: vaccinated vs unvaccinated and 5W vs 9W



# Can't really do the usual plots such as box plots, density plots, histograms, etc
# because there are too many variables !!! Hence why we are doing a multivariate analysis
# so can for instance perform pca and look at then colour code observations by their
# vaccination status and/or age group

# One thing we can look at is the number of kids in each category
(tab <- table(dat[,c(2,3)]))
round(prop.table(tab),3)
rowSums(prop.table(tab)) # roughly equal split between 9W and 5W
colSums(prop.table(tab)) # almost double the proportion of vaccinated babies !


#-------------------
# CLUSTERING
#-------------------

# CLUSTERING OBSERVATIONS: 
# single and average linkage (agglomorative clustering) more likely to produce long stringy
# clusters and are thus more likely to produce outliers than aggolomorative complete linkage
# or divisive clustering

# 1. Hierarchical Clustering: 

# 1.1 Aggolomorative Complete Linkage Clustering
ahc_cl <- agnes(dat_scal,
               diss=F,
               metric="euclidean",
               method="complete")
plot(ahc_cl)
ahc_cl
dend <- dat_scal %>% dist %>% 
  hclust(method = "complete") %>% as.dendrogram %>%
  set("branches_k_color", k=3) %>% set("branches_lwd", 0.8) %>%
  set("labels_colors") %>% set("labels_cex", c(.9,1.2)) %>% 
  set("leaves_pch", 20) %>% set("leaves_cex", c(0.1,0.1))
ggd1 <- as.ggdend(dend)
ggplot(ggd1, labels = F) +
  scale_y_reverse(expand = c(0.2,0)) + 
  coord_polar(theta = "x")


# 1.2 Divisive Clustering
dhc_cl <- diana(dat_scal)
dhc_cl$dc   # Divisive coefficient: 0.7023754
# plot dendrogram
plot(dhc_cl)
dend2 <- dat_scal %>% dist %>% 
  diana(metric = "euclidean") %>% as.dendrogram %>%
  set("branches_k_color", k=3) %>% set("branches_lwd", 0.8) %>%
  set("labels_colors") %>% set("labels_cex", c(.9,1.2)) %>% 
  set("leaves_pch", 20) %>% set("leaves_cex", c(0.1,0.1))
ggd2 <- as.ggdend(dend2)
ggplot(ggd2, labels = F) +
  scale_y_reverse(expand = c(0.2,0)) + 
  coord_polar(theta = "x")





# 2. Non-Hierarchical Clustering: 
# It seems like we definitely have outliers
# so lets use k-means partitioning around medoids, since this is more robust to data anomolies
# such as outliers and missing values
# disadvantage of pam is that it is not good for large data sets

# K-Means Partitioning Around Mediods (PAM)
# we know that ideally there should be 2 or 4 groups, but lets look at the Gap stat anyway
set.seed(123)
fviz_nbclust(dat_scal, pam, method="gap_stat" )+theme_classic()
# suggests 8
# lets look at average silhoutte width
set.seed(123)
fviz_nbclust(dat_scal, pam, method="silhouette")+theme_classic()
# this suggests 2 !
# and within cluster sum of squares
set.seed(123)
fviz_nbclust(dat_scal, pam, method="wss")+theme_classic()
# very hard to see where the elbow is, so ignore
set.seed(123)
nhc_pam2 <- pam(dat_scal, 2)
summary(nhc_pam2)
plot(nhc_pam2, main ="")



# 3. DBSCAN
set.seed(123)
# note that the smallest statusXgroup has 9 observations
dbscan::kNNdistplot(dat_scal, k =  5)
# quite hard to tell, looks like either 5 or 8?
abline(h = 6, lty = 2)
set.seed(123)
db <- fpc::dbscan(dat_scal,
                  eps = 6,
                  MinPts = 5)
db
which(db$cluster==0)
# identifies only one cluster and classifies everything else as outliers/noise :(
# again 34 and 39 are picked out as outliers



# VISUALISATION:
# lets focus on results from complete linkage and pam
# lets see the breakdown of age group and status in each cluster

# Complete Linkage:
ahc_cut <- cutree(ahc_cl, k=3)
length(ahc_cut[ahc_cut==1]) # number of elements in cluster 1
length(ahc_cut[ahc_cut==2]) # number of elements in cluster 2
length(ahc_cut[ahc_cut==3]) # number of elements in cluster 3
ahc_dat <- dat[,c(2,3)]
ahc_dat$Cluster <- as.vector(ahc_cut)
ahc_dat1 <- ahc_dat[ahc_dat$Cluster==1,-3]
ahc_dat2 <- ahc_dat[ahc_dat$Cluster==2,-3]
ahc_dat3 <- ahc_dat[ahc_dat$Cluster==3,-3]
# Cluster 1
ahc_tab1 <- data.frame(
       Group=c("Unvaccinated 5W", "Unvaccinated 9W", "Vaccinated 5W", "Vaccinated 9W"),
       value=c(table(ahc_dat1)[1,1],table(ahc_dat1)[2,1],table(ahc_dat1)[1,2],table(ahc_dat1)[2,2]) )
ahc_tab1$value <- ahc_tab1$value/sum(ahc_tab1$value)*100
labels <- c(round(ahc_tab1$value[1]), round(ahc_tab1$value[2]), round(ahc_tab1$value[3]), round(ahc_tab1$value[4]))
labels <- c(paste(labels[1],"%"), paste(labels[2],"%"), paste(labels[3],"%"), paste(labels[4],"%"))
ggplot(ahc_tab1, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.5), size=10) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 25, face = "bold"), legend.text=element_text(size=23))+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))
# Cluster 2
ahc_tab2 <- data.frame(
  Group=c("Vaccinated 5W", "Vaccinated 9W"),
  value=c(table(ahc_dat2)[1,2],table(ahc_dat2)[2,2]) )
ahc_tab2$value <- ahc_tab2$value/sum(ahc_tab2$value)*100
labels <- c(round(ahc_tab2$value[1]), round(ahc_tab2$value[2]))
labels <- c(paste(labels[1],"%"), paste(labels[2],"%"))
ggplot(ahc_tab2, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.5), size=10) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 25, face = "bold"), legend.text=element_text(size=23))+
  scale_fill_manual(values=c("#CCECE6", "#99D8C9"))
# Cluster 3
ahc_tab3 <- data.frame(
  Group=c("Vaccinated 9W"),
  value=c(1))
ahc_tab3$value <- ahc_tab3$value/sum(ahc_tab3$value)*100
labels <- c(round(ahc_tab3$value[1]))
labels <- c(paste(labels[1],"%"))
ggplot(ahc_tab3, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0), size=10) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 25, face = "bold"), legend.text=element_text(size=23))+
  scale_fill_manual(values=c("#99D8C9"))



# Divisive:
dhc_cut <- cutree(dhc_cl, k=3)
length(dhc_cut[dhc_cut==1]) # number of elements in cluster 1
length(dhc_cut[dhc_cut==2]) # number of elements in cluster 2
length(dhc_cut[dhc_cut==3]) # number of elements in cluster 3
dhc_dat <- dat[,c(2,3)]
dhc_dat$Cluster <- as.vector(dhc_cut)
dhc_dat1 <- dhc_dat[dhc_dat$Cluster==1,-3]
dhc_dat2 <- dhc_dat[dhc_dat$Cluster==2,-3]
dhc_dat3 <- dhc_dat[dhc_dat$Cluster==3,-3]
# Cluster 1
dhc_tab1 <- data.frame(
  Group=c("Unvaccinated 5W", "Unvaccinated 9W", "Vaccinated 5W", "Vaccinated 9W"),
  value=c(table(dhc_dat1)[1,1],table(dhc_dat1)[2,1],table(dhc_dat1)[1,2],table(dhc_dat1)[2,2]) )
dhc_tab1$value <- dhc_tab1$value/sum(dhc_tab1$value)*100
labels <- c(round(dhc_tab1$value[1]), round(dhc_tab1$value[2]), round(dhc_tab1$value[3]), round(dhc_tab1$value[4]))
labels <- c(paste(labels[1],"%"), paste(labels[2],"%"), paste(labels[3],"%"), paste(labels[4],"%"))
ggplot(dhc_tab1, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.5), size=10) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 25, face = "bold"), legend.text=element_text(size=23))+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))
# Cluster 2
dhc_tab2 <- data.frame(
  Group=c("Vaccinated 5W", "Vaccinated 9W"),
  value=c(table(dhc_dat2)[1,2],table(dhc_dat2)[2,2]) )
dhc_tab2$value <- dhc_tab2$value/sum(dhc_tab2$value)*100
labels <- c(round(dhc_tab2$value[1]), round(dhc_tab2$value[2]))
labels <- c(paste(labels[1],"%"), paste(labels[2],"%"))
ggplot(dhc_tab2, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.5), size=10) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 25, face = "bold"), legend.text=element_text(size=23))+
  scale_fill_manual(values=c("#CCECE6", "#99D8C9"))
# Cluster 3
dhc_tab3 <- data.frame(
  Group=c("Vaccinated 9W"),
  value=c(1))
dhc_tab3$value <- dhc_tab3$value/sum(dhc_tab3$value)*100
labels <- c(round(dhc_tab3$value[1]))
labels <- c(paste(labels[1],"%"))
ggplot(dhc_tab3, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0), size=10) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 25, face = "bold"), legend.text=element_text(size=23))+
  scale_fill_manual(values=c("#99D8C9"))




# PAM
nhc_dat <- dat[,c(2,3)]
length(nhc_dat[nhc_dat==1]) # number of elements in cluster 1
length(nhc_dat[nhc_dat==2]) # number of elements in cluster 2
nhc_dat$Cluster <- as.vector(nhc_pam2$clustering)
nhc_dat1 <- nhc_dat[nhc_dat$Cluster==1,-3]
nhc_dat2 <- nhc_dat[nhc_dat$Cluster==2,-3]
# Cluster 1
nhc_tab1 <- data.frame(
  Group=c("Unvaccinated 5W", "Unvaccinated 9W", "Vaccinated 5W", "Vaccinated 9W"),
  value=c(table(nhc_dat1)[1,1],table(nhc_dat1)[2,1],table(nhc_dat1)[1,2],table(nhc_dat1)[2,2]) )
nhc_tab1$value <- nhc_tab1$value/sum(nhc_tab1$value)*100
labels <- c(round(nhc_tab1$value[1]), round(nhc_tab1$value[2]), round(nhc_tab1$value[3]), round(nhc_tab1$value[4]))
labels <- c(paste(labels[1],"%"), paste(labels[2],"%"), paste(labels[3],"%"), paste(labels[4],"%"))
ggplot(nhc_tab1, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.5), size=10) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 25, face = "bold"), legend.text=element_text(size=23))+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))
# Cluster 2
nhc_tab2 <- data.frame(
  Group=c("Vaccinated 5W", "Vaccinated 9W"),
  value=c(table(nhc_dat2)[1,2],table(nhc_dat2)[2,2]) )
nhc_tab2$value <- nhc_tab2$value/sum(nhc_tab2$value)*100
labels <- c(round(nhc_tab2$value[1]), round(nhc_tab2$value[2]))
labels <- c(paste(labels[1],"%"), paste(labels[2],"%"))
ggplot(nhc_tab2, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.5), size=10) +
  coord_polar("y", start=0)+ 
  theme_void()+ 
  theme(legend.title = element_text(size = 25, face = "bold"), legend.text=element_text(size=23))+
  scale_fill_manual(values=c("#CCECE6", "#99D8C9"))




# EXAMINE OBSERVATIONS IN EACH CLUSTER

# Complete Linkage
which(ahc_dat$Cluster==1)
which(ahc_dat$Cluster==2)
which(ahc_dat$Cluster==3)  # 39

# Divisive 
which(dhc_dat$Cluster==1)
which(dhc_dat$Cluster==2)
which(dhc_dat$Cluster==3)  # 39

# PAM
which(nhc_dat$Cluster==1)
which(nhc_dat$Cluster==2)

# seems like the second and third clusters all contain outliers






