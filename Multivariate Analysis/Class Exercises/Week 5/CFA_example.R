# Confirmatory Factor Analysis using Lavaan in R
library(lavaan)

CFAdata<-read.csv(file.choose())
CFAdata[1:10,]


Ch1.model1<- 'Conservative=~ x1 + x2 + x3+x4+x5+x6+x7+x8+x9'
fit1 <- cfa(Ch1.model1, data = CFAdata)
summary(fit1,fit.measures=TRUE)       

MI <- modificationIndices(fit1)    
MI    
subset(MI, mi > 20)      

Ch1.model2<- 'Conservative=~ x1 + x3+x4+x5+x6+x7+x9
                                          x3~~x4 '
fit2 <- cfa(Ch1.model2, data = CFAdata)
summary(fit2,fit.measures=TRUE) 
Est2<- parameterEstimates(fit2, ci = FALSE, standardized = TRUE)
Est2
#subset(Est, op == "=~")

Ch1.model3<- 'Depress=~ x11 + x12 + x13'
fit3 <- cfa(Ch1.model3, data = CFAdata)
summary(fit3,fit.measures=TRUE)  

Ch1.model4<- 'Conservative=~ x1+ x3+x4+x5+x6+x7+x9
                          x3~~x4
                         Depress =~ x11 + x12 + x13
                         Conservative~~Depress'

fit4 <- cfa(Ch1.model4, data = CFAdata)

summary(fit4,fit.measures=TRUE)
Est4<- parameterEstimates(fit4, ci = FALSE, standardized = TRUE)
Est4


