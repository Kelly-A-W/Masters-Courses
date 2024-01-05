
library(princurve)
x <- runif(100,-1,1) # generating 100 random uniform numbers from -1 to 1
x <- cbind(x, x ^ 2 + rnorm(100, sd = 0.1)) # one column is x, the other is x squared plus noise
View(x)
X<-as.data.frame(x)
plot(X$x,X$V2)
fit <- principal_curve(x)
plot(X$x,X$V2)
lines(fit)
whiskers(x, fit$s)



## load Podani and Miklos data-set-generating function from github
## Data set 1 - Fig1 in Podani & Miklos (2002) Ecology 83
podani1 <- function(species = 42, sites = 40, labels = TRUE) {
  mat <- matrix(0, ncol = species, nrow = sites)
  diag(mat) <- 1
  mat[seq(sites + 1, prod(sites, species), by = sites + 1)] <- 3
  mat[seq((2 * sites + 1), prod(sites, species),
          by = sites + 1)] <- 1
  if(labels) {
    colnames(mat) <- paste("Spp", seq_len(species), sep = "")
    rownames(mat) <- paste("", seq_len(sites), sep = "")
  }
  mat
}

## Data set 2 - Fig 2 in Podani & Miklos (2002) Ecology 83
podani2 <- function(species = 6, sites = 40, labels = TRUE) {
  matrix(c((g <- rep(1:40, species / 2)),
           rep(rev(g), species / 2)),
         ncol = species, nrow = sites)
}

## Data set 3 - Fig9 in Podani & Miklos (2002) Ecology 83
podani3 <- function(species = 11, sites = 19,
                    abund = c(1,2,3,7,8,7,4,2,1),
                    labels = TRUE) {
  mat <- matrix(0, ncol = species, nrow = sites)
  ind <- c(outer(seq_len(length(abund)),
                 seq(0, by = sites + 1, length = species), "+"))
  mat[ind] <- abund
  if(labels) {
    colnames(mat) <- paste("Spp", seq_len(species), sep = "")
    rownames(mat) <- paste("", seq_len(sites), sep = "")
  }
  mat
}

## generate data set 1 of Podani & Miklos
p1 <- podani1()
View(p1)
## ordinate
library("vegan")
pca1 <- rda(p1)
plot(pca1)

## plot data and PCA
layout(matrix(c(1,2,2), nrow = 3))
matplot(p1, type = "o", pch = 19, ylab = "Abundance", xlab = "Gradient")
ordipointlabel(pca1, display = "sites")
layout(1)

## Load analogue
library("analogue")
View(p1)
dim(p1)
prc1 <- prcurve(p1, trace = TRUE, plotit = TRUE, thresh = 0.0005, maxit = 50)
prc1
varExpl(pca1)
prc1$s
plot(prc1)

## Load the Abernethy Forest data
data(abernethy)
View(abernethy)
## Remove the Depth and Age variables
abernethy2 <- abernethy[, -(37:38)]
help(abernethy)
dim(abernethy2)
pca1<-prcomp(abernethy2)
biplot(pca1)
## Fit the principal curve varying complexity for each species
aber.pc <- prcurve(abernethy2, trace = TRUE, vary = TRUE, penalty = 1.4)
aber.pc
attributes(aber.pc)
aber.pc$s
aber.pc$dist
fit<-fitted(aber.pc)
fit
## Plot
op <- par(mar = c(5,4,2,2) + 0.1)
plot(aber.pc)

varExpl(rda(abernethy2), axes = 1:5, cumulative = TRUE)

ord <- rda(abernethy2)
plot(ord, scaling = 1)
lines(aber.pc, scaling = 1)