##############
## Starter code for fitting joint models:
## (1) A multivariate GLMM using the lme4 package
## (2) A LVM using the boral package
###############

## We will use the spider dataset available in mvabund, pitfall trap counts of 12 hunting spider species at 28 sites. First load it from mvabund:
library(mvabund)
data(spider)
head(spider$abund[,1:6])
head(spider$abund[,7:12])
head(spider$x)

require(graphics)

spiddat <- as.mvabund(spider$abund)

plot(spiddat)
## Cherry-pick a few variables (lme4 can't handle many responses, especially when we only have 28 sites)
spid.abund = spider$abund[,c(1:3,7,8,12)]
spid.abund = spider$abund
##########
## (1) Fitting a multivariate GLMM using the lme4 package ##
## Use lme4 to fit a multivariate GLMM like in equation 1 of main text, i.e. (random) site effect, and a quadratic effect of soil for each spp. 
## Using a Poisson model with a log link, other options are available via the family argument, as usual.
## Note this function is super-fussy about data - the model didn't converge even for these six species!
##########   

spid.vec = c(as.matrix(spid.abund))
n.site = dim(spid.abund)[1]
n.spp = dim(spid.abund)[2]
n.spp
spid.vec
## construct a data frame with soil (standardized), species labels, and site labels:
X = data.frame(soil = rep(scale(spider$x[,1]), n.spp), spp = rep(dimnames(spid.abund)[[2]], each=n.site), site = rep(dimnames(spid.abund)[[1]], n.spp) )
View(X)
## fit the GLMM using lme4 and look at results
library(lme4)
fit.glmm = glmer(spid.vec~0+spp+spp:soil+spp:I(soil^2)+(1|site)+(0+spp|site), data=X, family=poisson())
# Returns some warnings about non-convergence -- the best solution to this is more samples and less variables!!

## The key term in the above is (0+spp|site), which introduces a multivariate random intercept at each site.
## This technique will work for count data but for binary or ordinal data, the variance of the random effect
## needs to be fixed (e.g. at one) for identifiability, which lme4 doesn't do - instead something like MCMCglmm might work.

print(summary(fit.glmm),correlation=FALSE) # To look at estimated parameter values
confint(fit.glmm, method="Wald") # 95% confidence intervals for model parameters. Very approximate though, other more computationally intensive options are available.

## Use corrplot package to calculate and plot residual correlations between species, e.g. possibly due to species interaction etc...
library(corrplot)
vrcorrs=VarCorr(fit.glmm)
corrs=attr(vrcorrs$site.1,"corr")
corrplot(corrs, diag = F, type = "lower", title = "Residual correlations from GLMM", method = "color", tl.srt = 45)



##########
## (2) Fitting a Latent Variable Model using the boral package ##
## Use boral to fit a LVM like in equation 3 of main text, i.e. two latent variables, (random) site effect, and a quadratic effect of soil. 
## Also note that this can fit models to more data (e.g. try full dataset) but takes longer to fit it - for full spdier it will still be a few min, for full aravo dataset over half an hour.
##########   

## Covariates need to be stored as a matrix for boral:
covX <- cbind(scale(spider$x[,1]),scale(spider$x[,1])^2)
View(covX)
## fit the LVM using boral and look at results.  Need version 0.7 or later, available on CRAN.
## NOTE: boral uses JAGS so you need to first install JAGS-3.0.0.exe or higher from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

library(boral)
fit1.lvm<-boral(y=spid.abund,family="poisson",lv.control=list(num.lv=2,type="independent"),row.eff="fixed")

summary(fit1.lvm)
par(mfrow=c(2,2))
plot(fit1.lvm)
par(mfrow=c(1,1))
lvsplot(fit1.lvm,alpha=0.55,main="Unconstrained",jitter=TRUE,return.vals=TRUE,which.lvs=c(1,2))

fit2.lvm <- boral(y = spider$abund, X = spider$x ,lv.control=list(num.lv=2,type="independent"), family = "poisson", row.eff = "random", save.model = TRUE)
summary(fit2.lvm) # To look at estimated parameter values
fit.lvm$hpdintervals # 95% credible intervals for model parameters.
par(mfrow=c(2,2))
plot(fit2.lvm)
par(mfrow=c(1,1))
lvsplot(fit2.lvm,alpha=0.55,main="Residual biplot",jitter=TRUE,return.vals=TRUE,which.lvs=c(1,2))

envcors<-get.enviro.cor(fit2.lvm)
rescors<-get.residual.cor(fit2.lvm)
library(corrplot)
par(mfrow=c(1,2))
corrplot(envcors$sig.cor,type="lower",diag=F,title="Correlations due to covariates",mar=c(3,0.5,2,1),tl.srt=45)
corrplot(rescors$sig.cor,type="lower",diag=F,title="Residual correlations",mar=c(3,0.5,2,1),tl.srt=45)


data(antTraits)
help(antTraits)
head(antTraits$abund)
head(antTraits$env)
head(antTraits$traits)

y <- antTraits$abun
dim(y)
sel.spp <- colSums(y>0)>4
y <- y[,sel.spp]
dim(y)
X <- antTraits$env ## Scale covariates to ease interpretability
traits <-antTraits$traits[sel.spp,]
## An intercept column must be included if species-specific intercepts are to be treated as a function of species traits as well, i.e. species prevalence and well as their environmental response is driven by traits
which.traits <- vector("list",ncol(X)+1)
for(i in 1:length(which.traits)) which.traits[[i]] <- 1:ncol(traits)

## Tighter priors required on negative binomial dispersion parameter.
## This can take a while!!!
fit.traits <- boral(y, X = X, traits = traits, which.traits = which.traits, family = "negative.binomial", lv.control=list(num.lv=2,type="independent"), calc.ics = FALSE, save.model = TRUE, hypparams = c(100,20,100,50))

summary(fit.traits)