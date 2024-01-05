rm(list = ls(all = TRUE))

# INSTALL LIABRARIES
library(ggplot2)



#=================================
# Question 1(a)
#=================================

# Simulate a dataset
N     <- 30
set.seed(123)
x     <- runif(N, min=-1, max=1)
noise <- rnorm(N)   # need a different epsilon for each y
y     <- 0.8*x + noise

# Fit two linear models
g1    <- lm(y~0+x, offset = rep(0.5,N))
summary(g1)
# lets check that the intercept was in fact added
predict(g1) - (g1$coefficients*x +0.5)  # all zero, so g1 is doing what it is supposed to !
g2    <- lm(y~0+x, offset = rep(-0.5,N))
summary(g2)


# Plot data, the true underlying model as well g1 and g2
plot(c(min(x),max(x)),c(-4,4),main="",type="n",xlab="x",ylab="y")
points(x,y,pch=19)         # data
lines(0.8*x ~x,type="l",lwd=2)  # true underlying model
abline(g1, col = "blue")   # g1
abline(g2, col = "green")  # g2
grid()

domain <- seq(-1, 1, 0.01)
predg1 <- 0.5 + g1$coefficients*domain
predg2 <- -0.5 + g2$coefficients*domain
df <- data.frame(domain = domain, target = 0.8*domain, g1 = predg1, g2 = predg2)
df <- pivot_longer(df, cols = 2:4, names_to = "Y", values_to = "Value")
df
pointdata <- data.frame(x = x, y=y)
ggplot()+
  geom_line(df,mapping =aes(domain, Value, group=Y, color=Y))+
  geom_point(data = pointdata, mapping = aes(x = x, y = y))+
  scale_x_continuous("x")+
  scale_y_continuous("y")+
  theme_minimal(base_size = 22)




#=================================
# Question 1(b)
#=================================

M <- 10000
N <- 30
i <- seq(5,25,1)

# Set ups domain and target function
domain  <- seq(-1,1,0.01)
target  <- 0.8*domain

# Create Storage
set.seed(123)
xdat <- list()
ydat <- list()
min_val_err <- matrix(NA, nrow = M, ncol = length(i))
out_err     <- matrix(NA, nrow = M, ncol = length(i))

for (j in 1:M) {
  xdat[[j]] <- runif(N, min=-1, max=1)
  noise     <- rnorm(N) 
  ydat[[j]] <- 0.8*xdat[[j]] + noise
  dat       <- data.frame(x=xdat[[j]], y=ydat[[j]])
  
  for (k in 1:length(i)) {
    ind   <- sample(1:N, i[k], replace = F)
    val   <- dat[ind,]
    train <- dat[-ind,]
    
    # Fit models
    g1    <- lm(train$y~0+train$x, offset = rep(0.5,nrow(train)))
    g2    <- lm(train$y~0+train$x, offset = rep(-0.5,nrow(train)))
    
    # Calculate validation error
    val_err_g1 <- 1/(2*nrow(val))*sum((val$y - (g1$coefficients*val$x + 0.5))^2)  # MSE
    val_err_g2 <- 1/(2*nrow(val))*sum((val$y - (g2$coefficients*val$x - 0.5))^2) 
    
    # Choose best fitting function
    min_val_err[j,k]  <- list(val_err_g1, val_err_g2)[[which.min(c(val_err_g1,val_err_g2))]]
    gstar             <- list(g1,g2)[[which.min(c(val_err_g1,val_err_g2))]]
    
    # Calculate and store Eout 
    fit_gstar         <- gstar$coefficients*domain + gstar$offset[1]
    out_err[j,k]      <- (t(fit_gstar - target)%*%(fit_gstar - target))*(domain[2]-domain[1])*1/(max(domain)-min(domain)) + 1
    
    
  }
}



min_val_err
out_err


# Calculate expectations for each value of i
(E_min_val_err <- colMeans(min_val_err))
(E_out_err     <- colMeans(out_err))



# Plots
df <- data.frame(i = i, Eval = E_min_val_err, Eout=E_out_err)
df <- pivot_longer(df, cols = 2:3, names_to = "EE", values_to = "Value")
df
ggplot(df, aes(i, Value, group = EE, color = EE)) +
  geom_line(size=1)+ 
  geom_point(size=3) +
  scale_x_continuous("Valdation Set Size (i)")+
  scale_y_continuous(name="Expected Error")+
  theme_minimal(base_size = 22)






