rm(list = ls(all = TRUE))

#------------------------------------
# Load Libraries
#------------------------------------

library(ggplot2)
library(tidyr)



#---------------------------------
# QUESTION 2(a)
#---------------------------------

# Legendre polynomials
Legendre=function(x,q){     # q = order of the Legendre polynomial
  val=0
  for(k in 0:q){
    val=val+((x^k)*choose(q,k)*choose((q+k-1)/2,q))
  }
  return((2^q)*val)
}




x<-seq(-1,1,0.01)
#par(mar = c(5.1, 4.1, 4.1, 2.1))
df <- data.frame(x = x, 
                 "0" = Legendre(x,0), 
                 "1" = Legendre(x,1), 
                 "2" = Legendre(x,2), 
                 "3" = Legendre(x,3), 
                 "4" = Legendre(x,4), 
                 "5" = Legendre(x,5))
df <- pivot_longer(df, cols = 2:7, names_to = "Order", values_to = "Value")
df$Order <- as.factor(df$Order)
levels(df$Order) <- c(0,1,2,3,4,5)
ggplot(df, aes(x, Value, group = Order, color = Order)) +
  geom_line()+ 
  scale_x_continuous("x")+
  scale_y_continuous(name=expression(L[q](x)))+
  theme_minimal(base_size = 22)



#---------------------------------
# QUESTION 2(b)
#---------------------------------


# Creates random polynomial target functions: Equation (3)
rLegfunc=function(x,Q_f){      # Q_f = order of the polynomial f(x)
  val=0
  beta=runif((Q_f+1),-1,1)
  for(q in 0:Q_f){
    val=val+(Legendre(x,q)*beta[q+1])    
  }
  return(val)
}


x<-seq(-1,1,0.01)
set.seed(123)
df <- data.frame(x = x,
                 "3" = rLegfunc(x,3), 
                 "4" = rLegfunc(x,4), 
                 "5" = rLegfunc(x,5))
df <- pivot_longer(df, cols = 2:4, names_to = "Order", values_to = "Value")
df$Order <- as.factor(df$Order)
levels(df$Order) <- c(3,4,5)
ggplot(df, aes(x, Value, group = Order, color = Order)) +
  geom_line()+ 
  scale_x_continuous("x")+
  scale_y_continuous(name=expression(f(x)))+
  theme_minimal(base_size = 22)





# Creates random polynomial target functions: Equation (2)
rPolyfunc=function(x,Q_f){
  val=0
  alpha=runif((Q_f+1),-1,1)
  for(q in 0:Q_f){
    val=val+(x^q*alpha[q+1])
  }
  return(val)
}


x<-seq(-1,1,0.01)
set.seed(123)
df <- data.frame(x = x,
                 "3" = rPolyfunc(x,3), 
                 "4" = rPolyfunc(x,4), 
                 "5" = rPolyfunc(x,5))
df <- pivot_longer(df, cols = 2:4, names_to = "Order", values_to = "Value")
df$Order <- as.factor(df$Order)
levels(df$Order) <- c(3,4,5)
ggplot(df, aes(x, Value, group = Order, color = Order)) +
  geom_line()+ 
  scale_x_continuous("x")+
  scale_y_continuous(name=expression(f(x)))+
  theme_minimal(base_size = 22)









