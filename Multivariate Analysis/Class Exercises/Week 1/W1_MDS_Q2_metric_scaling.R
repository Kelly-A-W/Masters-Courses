
# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Install packages
# install.packages('psychTools')
library(psychTools)
# install.packages("smacof")
library(smacof)

# Load data
data(cities)
cities

#------------------------------
# Question 2.a:
#------------------------------
# (i) Use the mds() function – in the smacof package – to perform metric multidimensional
# scaling. 
# Locate the cities in t=2 dimensions, starting with a classical scaling solution 
# as the initial configuration. 

mds1 <- mds(cities,
            ndim=2,             # 2 dimensions
            type="ratio",       # metric MDS
            init="torgerson")   # use classical scaling as initial configuration

mds1           # stress value seems to be relatively low, suggesting good fit


# (ii) Which of the 11 cities has the poorest fit? Explain your answer.

summary(mds1)  # we can see that SEA contributes the largest percentage to the stress, 
# a massive 33%. Thus SEA contributes the most to the lack of fit and thus will have
# the poorest fit. LAX also contributes quite highly to the stress at 20.15%, and thus has
# the second poorest fit. The remaining cities all contribute relatively little to the
# overall stress. 


# (iii) Plot the 2D multidimensional scaling configuration and compare this the 
# locations of the cities on a map from the atlas.

plot(mds1)
plot(city.location, xlab="Dimension 1", ylab="Dimension 2",
     main ="Multidimensional scaling of US cities", type="n", xlim=c(65,130))
text(city.location,labels=names(cities), pos=4)
# mds1 locations need to be reflected along the horizontal axis
# i.e. y values must become -y
city.loc <- mds1$conf 
city.loc[,2] <- -city.loc[,2]  # flip
city.loc <- psych::rescale(city.loc,apply(city.location,2,mean),apply(city.location,2,sd))
points(city.loc,type="n") #add the date point to the map
text(city.loc,labels=names(cities), pos=4, col = "red")

# you can see that the 2D mds configuration (red) is very close to the actual locations of
# the airports (black), however they do not quite match up. While some are very close to the 
# real physical locations such as MSY and ATL, many are a bit off such as BOS and JFK.
# The reason for this discrepancy is because we used a 2D surface to display the 
# distances when in reality the airports are located on a globe, so in 3D, meaning that
# 2D euclidean distance will not give us the exact real life locations 
# of these airports and the results will be distorted.


#------------------------------
# Question 2.b:
#------------------------------
# Repeat the analysis in part (a) using a random starting solution. 
# How does the different initial configuration change the location of the cities 
# in the 2D configuration?

# (i) 
set.seed(123)
mds2 <- mds(cities,
            ndim=2,             # 2 dimensions
            type="ratio",       # metric MDS
            init="random")      # use random starting solution

mds2    # stress value is slightly higher than before but still relatively low, 
        # suggesting good fit


# (ii) 
summary(mds2)  
# As before SEA has the highest stress contribution at 38.52%, slightly more than before,
# and thus again has the poorest fit. LAX fits much better with its stress contribution
# going down to only 3.29%, but BOS and DEN fit much worse with stress contributions
# of 16.03% and 11.89% respectively.


# (iii) 
plot(mds2)
plot(city.location, xlab="Dimension 1", ylab="Dimension 2",
     main ="Multidimensional scaling of US cities", type="n", xlim=c(65,130))
text(city.location,labels=names(cities), pos=4)
# mds2 locations need to be reflected along the vertical axis
city.loc2 <- mds2$conf 
city.loc2[,1] <- -city.loc2[,1]  # flip
city.loc2 <- psych::rescale(city.loc2,apply(city.location,2,mean),apply(city.location,2,sd))
points(city.loc2,type="n") 
text(city.loc2,labels=names(cities), pos=4, col = "red")
# immediately we can see that the mds2 configuration (red) did not do a good job of 
# approximating that actual locations of the airports (black). A fee airports such as DEN and
# ATL are close to their actaul locations but most are not. 

# lets add mds1
points(city.loc,type="n") 
text(city.loc,labels=names(cities), pos=4, col = "blue")
# as expected from the previous plot, the mds1 configuration is quite different to the 
# mds2 configuration. 
# this highlights that these results are not unique, and that different starting 
# configurations can result in different results. A good starting configuration is
# very important when performing mds, and multiple starting configurations should
# always be tried.

