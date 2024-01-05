
men<-read.csv(file.choose())
head(men)
cov_men<-cov(men[,3:10])
cov_men
eigencov<-eigen(cov_men)
eigencov

X<-as.matrix(men[,3:10])

rrr_men <- rrr(X, X, k=0, type="Identity", rrr.plot=TRUE,
                missing="omit", title="")
rrr_men$A[8]

FA<-factanal(covmat=corr_men,correlation,factors=2,
             rotation="none")
FA

