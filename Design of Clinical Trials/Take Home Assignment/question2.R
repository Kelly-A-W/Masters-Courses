library(daewr)


#---------------------------------------------------(b)
t = 10
k = 4
factorial(t)/(factorial(k)*factorial(t-k)) # 210


#---------------------------------------------------(c)

BIBsize(10, 4)  # 15


#---------------------------------------------------(d)
library(AlgDesign)

set.seed(123)

BIB <- optBlock(~., withinData = factor(1:10), blocksizes = rep(4,15))

des <- matrix(BIB$rows, nrow = 15, ncol = 4, byrow = TRUE,
              dimnames = list(c( "Block1", "Block2", "Block3", "Block4","Block5", "Block6", "Block7", "Block8", "Block9","Block10", "Block11", "Block12", "Block13", "Block14","Block15"), 
                              c("unit1", "unit2", "unit3", "unit4")))
des









