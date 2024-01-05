# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#------------------------------------
# Load Libraries
#------------------------------------

library(deSolve)
library(gtools)
library(ggplot2)
library(tidyr)
library(gganimate)
library(transformr)


#------------------------------------
# Load Data
#------------------------------------

dat <- read.csv("TB_notifications_2022-08-17.csv")
str(dat)
dat$country <- as.factor(dat$country)

# Get South Africa Data
dat_za <- dat[which(dat$country=="South Africa"), ]
# Remove unnecessary variables (country indentifiers)
dat_za <- dat_za[-(1:5)]
str(dat_za)

# Extract variables of interest
mydat_za <- data.frame(year = dat_za$year,
                       child = dat_za$newrel_f014 + dat_za$newrel_m014,       # (reported) incidence of TB in children (regardless of HIV status)
                       adult = (dat_za$newrel_f15plus+dat_za$newrel_m15plus), # (reported) incidence of TB in adults (regardless of HIV status)
                       hiv = dat_za$newrel_hivpos)                            # (reported) incidence of TB in HIV positive individuals (regardless of age)
mydat_za
# Only have data from 2013 onwards
# subset:
mydat_za_2013 <- mydat_za[34:41,]
mydat_za_2013




df <- pivot_longer(mydat_za_2013, cols = 2:4, names_to = "Group", values_to = "Value")
df
ggplot(df, aes(year, Value, group = Group, color = Group)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2013,2020,1))+
  scale_y_continuous(name="Reported Incidence", labels = scales::comma, breaks=seq(0,275709,50000))+
  theme_minimal(base_size = 22)



#------------------------------------
# Create SIR Model
#------------------------------------


sir <- function(t, x, logkappa)  {
  with(as.list(c(logkappa, x)), {
    
    # Interaction Terms:
    kappa    = exp(logkappa)
    beta_cc  = 0.005*kappa   # not differentiating between contact between HIV neg and HIV pos people
    beta_ca  = 0.5*kappa
    beta_aa  = 1*kappa
    beta_ac  = 0.1*kappa
    
    # HIV-Negative Children
    p_nc     = 0.125      # proportion of children who develop active TB
    
    # HIV-Negative Adults
    p_na     = 0.145
    
    
    # HIV-Negative (children and adults)
    dtb_n   = 0.6      # death rate due to TB 
    d_n     = 11/1000  # death rate 
    x_n     = 0.33     # protection against development of reinfection or developing active TB due to being latently infected or recovered 
    r_n     = 0.01     # annual risk of relapse from recivered to active TB
    theta_n = 0.0002   # rate of progression to active disease due to existing infection 
    
    
    # HIV-Positive (children and adults)
    dtb_h   = 0.9       # death rate due to TB in HIV positive individuals 
    d_h     = (1-0.86) +  11/1000   # death rate for HIV positive individuals
    x_h     = 0.755     # protection against development of reinfection or developing active TB due to being latently infected or recovered fo HIV positive 
    r_h     = 0.425
    theta_h = 0.00236   # rate of progression to active disease due to existing infection 
    
    # HIV-Positive Adults
    p_ha    = 0.55
    
    # HIV-Positive Children
    p_hc    = 0.475
    mu      = 0.01    # proportion of hiv positive births
    
    # Whole Population
    t      = 0.58      # proportion of population who receive treatment 
    b      = 21.5/1000  # birth rate
    a      = 1/14       # aging rate 
    gammat = 0.764      # treatment recovery rate
    gamman = 0.32       # natural recovery rate
    s      = 60/365     # rate of treatment seeking 
    nu     = 0.185      # rate of HIV infection
    
    I_c = In_nc + It_nc + In_hc + It_hc
    I_a = In_na + It_na + In_ha + It_ha
    
    N_c = S_nc + L_nc + In_nc + It_nc + T_nc + R_nc + S_hc + L_hc + In_hc + It_hc + T_hc + R_hc
    N_a = S_na + L_na + In_na + It_na + T_na + R_na + S_ha + L_ha + In_ha + It_ha + T_ha + R_ha
    N = N_a + N_c
    
    lambda_c <- beta_cc*I_c/N_c + beta_ca*I_a/N_a
    lambda_a <- beta_aa*I_a/N_a + beta_ac*I_c/N_c
    
    
    # HIV Negative Children 
    dS_nc  = b*(1-mu)*N -S_nc*(a + d_n + lambda_c)
    dL_nc  = S_nc*lambda_c*(1-p_nc) + R_nc*x_n*lambda_c*(1-p_nc) - L_nc*(a + d_n + theta_n + lambda_c*p_nc*x_n)
    dIn_nc = S_nc*lambda_c*p_nc*(1-t) + L_nc*(theta_n + lambda_c*p_nc*x_n)*(1-t) + R_nc*(r_n + lambda_c*p_nc*x_n)*(1-t) - In_nc*(d_n + a + dtb_n + gamman)
    dIt_nc = S_nc*lambda_c*p_nc*t + L_nc*(theta_n + lambda_c*p_nc*x_n)*t + R_nc*(r_n + lambda_c*p_nc*x_n)*t - It_nc*(d_n + a + s)
    dT_nc  = It_nc*s - T_nc*(d_n + a + gammat)
    dR_nc  = T_nc*gammat + In_nc*gamman - R_nc*(d_n + a + r_n + x_n*lambda_c)
    
    
    # HIV Negative Adults
    dS_na  = a*S_nc - S_na*(nu + d_n + lambda_a)
    dL_na  = a*L_nc + S_na*lambda_a*(1-p_na) + R_na*x_n*lambda_a*(1-p_na) - L_na*(nu + d_n + theta_n + lambda_a*p_na*x_n)
    dIn_na = a*In_nc + S_na*lambda_a*p_na*(1-t) + L_na*(theta_n + lambda_a*p_na*x_n)*(1-t) + R_na*(r_n + lambda_a*p_na*x_n)*(1-t) - In_na*(d_n + nu + dtb_n + gamman)
    dIt_na = a*It_nc + S_na*lambda_a*p_na*t + L_na*(theta_n + lambda_a*p_na*x_n)*t + R_na*(r_n + lambda_a*p_na*x_n)*t - It_na*(d_n + nu + s)
    dT_na  = a*T_nc + It_na*s - T_na*(d_n + nu + gammat)
    dR_na  = a*R_nc + T_na*gammat + In_na*gamman - R_na*(d_n + nu + r_n + x_n*lambda_a)
    
    
    
    # HIV Positive Children 
    dS_hc  = b*mu*N - S_hc*(a + d_h + lambda_c)
    dL_hc  = S_hc*lambda_c*(1-p_hc) + R_hc*x_h*lambda_c*(1-p_hc) - L_hc*(a + d_h + theta_h + lambda_c*p_hc*x_h)
    dIn_hc = S_hc*lambda_c*p_hc*(1-t) + L_hc*(theta_h + lambda_c*p_hc*x_h)*(1-t) + R_hc*(r_h + lambda_c*p_hc*x_h)*(1-t) - In_hc*(d_h + a + dtb_h + gamman)
    dIt_hc = S_hc*lambda_c*p_hc*t + L_hc*(theta_h + lambda_c*p_hc*x_h)*t + R_hc*(r_h + lambda_c*p_hc*x_h)*t - It_hc*(d_h + a + s)
    dT_hc  = It_hc*s - T_hc*(d_h + a + gammat)
    dR_hc  = T_hc*gammat + In_hc*gamman - R_hc*(d_h + a + r_h + x_h*lambda_c)
    
    
    # HIV Positive Adults
    dS_ha  = a*S_hc + nu*S_na - S_ha*(d_h + lambda_a)
    dL_ha  = a*L_hc + nu*L_na + S_ha*lambda_a*(1-p_ha) + R_ha*x_h*lambda_a*(1-p_ha) - L_ha*(d_h + theta_h + lambda_a*p_ha*x_h)
    dIn_ha = a*In_hc + nu*In_na + S_ha*lambda_a*p_ha*(1-t) + L_ha*(theta_h + lambda_a*p_ha*x_h)*(1-t) + R_ha*(r_h + lambda_a*p_ha*x_h)*(1-t) -In_ha*(d_h + dtb_h + gamman)
    dIt_ha = a*It_hc + nu*It_na + S_ha*lambda_a*p_ha*t + L_ha*(theta_h + lambda_a*p_ha*x_n)*t + R_ha*(r_h + lambda_a*p_ha*x_h)*t -It_ha*(d_h + s)
    dT_ha  = a*T_hc + nu*T_na + It_ha*s - T_ha*(d_h + gammat)
    dR_ha  = a*R_hc + nu*R_na + T_ha*gammat + In_ha*gamman - R_ha*(d_h + r_h + x_h*lambda_a)
    
    
    
    output <- list(c(dS_nc, dL_nc, dIn_nc, dIt_nc, dT_nc, dR_nc, 
                     dS_na, dL_na, dIn_na, dIt_na, dT_na, dR_na,
                     dS_hc, dL_hc, dIn_hc, dIt_hc, dT_hc, dR_hc, 
                     dS_ha, dL_ha, dIn_ha, dIt_ha, dT_ha, dR_ha))
    output
  })
}










#-------------------------- Initial values (2008)
N = 48687000  
N_c <- 5139800 + 5254100 + 5278900  
N_a <- N-N_c



start<-c(S_nc = (1-0.167)*13156686,
         L_nc = (1-0.167)*2288229,
         In_nc = (1-0.167)*146697,   
         It_nc = (1-0.167)*35107,
         T_nc = (1-0.167)*37615,
         R_nc = (1-0.167)*102503,
         S_na = (1-0.08)*24163579, 
         L_na = (1-0.08)*8121493,
         In_na = (1-0.08)*194124, 
         It_na = (1-0.08)*129415,
         T_na = (1-0.08)*138660, 
         R_na = (1-0.08)*266929, 
         S_hc = 0.167*13156686,
         L_hc = 0.167*2288229,
         In_hc = 0.167*146697,   
         It_hc = 0.167*35107,
         T_hc = 0.167*37615,
         R_hc = 0.167*102503,
         S_ha = 0.08*24163579, 
         L_ha = 0.08*8121493,
         In_ha = 0.08*194124, 
         It_ha = 0.08*129415,
         T_ha = 0.08*138660, 
         R_ha = 0.08*266929)

# check
as.numeric(sum(start)) 
as.numeric(N)




#-------------------------- Run unfit model

logkappa_start <- log(46)


## vector of timesteps
times <- seq(0, 12, 1)  # 2008 to 2020

run_d<-ode(times=times, y=start, func=sir,parms=logkappa_start,maxsteps=100000)
run_d



#------------------------------------
# Fit SIR Model: Least Squares
#------------------------------------
rep = 0.5      # reporting parameter

#Standardise data to use for fit 
dat_fit <- data.frame(scale(mydat_za_2013))



sir.sse<-function(p){
  
  model <- ode(times=times, y=start, func=sir,parms=p)
  data = dat_fit
  
  a      = 1/14       
  s      = 60/365      
  nu     = 0.185 
  
  count_nc <- as.data.frame(model)$It_nc*s 
  count_na <- a*as.data.frame(model)$T_nc + as.data.frame(model)$It_na*s 
  count_hc <- as.data.frame(model)$It_hc*s 
  count_ha <- a*as.data.frame(model)$T_hc + nu*as.data.frame(model)$T_na + as.data.frame(model)$It_ha*s 
  
  
  inc_a <- diff((count_nc + count_hc))    # Children receiving treatment (regardless of HIV Status) at baseline (2008) 
  inc_b <- diff((count_na + count_ha))    # Adults receiving treatment (regardless of HIV Status) at baseline (2008)
  inc_c <- diff((count_ha + count_hc))    # All HIV positive individuals receiving treatment (regardless of age) at baseline (2008)
  
  
  
  # Standardise inc's
  inc_a_scal <- data.frame(scale(inc_a[5:12]))
  inc_b_scal <- data.frame(scale(inc_b[5:12]))
  inc_c_scal <- data.frame(scale(inc_c[7:12]))
  
  error_a<-rep*inc_a_scal-(data$child)
  error_b<-rep*inc_b_scal-(data$adult)
  error_c<-rep*inc_c_scal-data$hiv[3:8]
  
  sse <- c(sum(error_a^2), sum(error_b^2), sum(error_c^2))   
  return(sum(sse))
}

sir.sse(p = logkappa_start)  # check current error


# Find Optimal log(kappa) estimate
fit_sse<-optimise(sir.sse, c(0,10)) 
fit_sse
(parm_fit <- fit_sse$minimum)  # fitted log(kappa)
exp(parm_fit)

# Create fitted model
times <- seq(0, 12, 1)  # 2008 to 2020
run_fit <- ode(times=times, y=start, func=sir,parms=parm_fit,maxsteps=100000)
run_fit
sir.sse(p = parm_fit)



#--------------------------  Make a Counter for Incidence 
a      = 1/14       
s      = 60/365      
nu     = 0.185 

count_nc <- as.data.frame(run_fit)$It_nc*s 
count_na <-  a*as.data.frame(run_fit)$T_nc + as.data.frame(run_fit)$It_na*s 
count_hc <- as.data.frame(run_fit)$It_hc*s 
count_ha <- a*as.data.frame(run_fit)$T_hc + nu*as.data.frame(run_fit)$T_na + as.data.frame(run_fit)$It_ha*s 



# Plot Incidence:
inc <- diff((count_nc+count_hc))  # Children receiving treatment (regardless of HIV Status)
df <- data.frame(fitted = rep*inc[5:12], reported= mydat_za_2013$child, year=mydat_za_2013$year)
df <- pivot_longer(df, cols = 1:2, names_to = "Incidence", values_to = "Value")
df
ggplot(df, aes(year, Value, group = Incidence, color = Incidence)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2013,2020,1))+
  scale_y_continuous(name="Incidence of TB in Children")+
  theme_minimal(base_size = 22)

inc <- diff((count_na+count_ha))   # Adults receiving treatment (regardless of HIV Status)
df <- data.frame(fitted = rep*inc[5:12], reported= mydat_za_2013$adult, year=mydat_za_2013$year)
df <- pivot_longer(df, cols = 1:2, names_to = "Incidence", values_to = "Value")
df
ggplot(df, aes(year, Value, group = Incidence, color = Incidence)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2013,2020,1))+
  scale_y_continuous(name="Incidence of TB in Adults")+
  theme_minimal(base_size = 22)

inc <- diff((count_ha+count_hc))   # All HIV positive individuals receiving treatment (regardless of age)
df <- data.frame(fitted = rep*inc[7:12], reported= mydat_za_2013$hiv[3:8], year=mydat_za_2013$year[3:8])
df <- pivot_longer(df, cols = 1:2, names_to = "Incidence", values_to = "Value")
df
ggplot(df, aes(year, Value, group = Incidence, color = Incidence)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2013,2020,1))+
  scale_y_continuous(name="Incidence of TB in HIV Positive Individuals")+
  theme_minimal(base_size = 22)




# Child HIV Negative Population
df <- as.data.frame(run_fit[,2:7])
df$year <- 2008:2020
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2008,2020,1))+
  scale_y_continuous(name="Number of HIV Negative Children", labels = scales::comma, breaks=seq(0,22521578,4000000))+
  theme_minimal(base_size = 22)

# Adult HIV Negative Population
df <- as.data.frame(run_fit[,8:13])
df$year <- 2008:2020
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2008,2020,1))+
  scale_y_continuous(name="Number of HIV Negative Adults", labels = scales::comma, breaks=seq(0,22521578,4000000))+
  theme_minimal(base_size = 22)

# Child HIV Positive Population
df <- as.data.frame(run_fit[,14:19])
df$year <- 2008:2020
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2008,2020,1))+
  scale_y_continuous(name="Number of HIV Positive Children", labels = scales::comma, breaks=seq(0,22521578,4000000))+
  theme_minimal(base_size = 22)

# Adult HIV Positive Population
df <- as.data.frame(run_fit[,20:25])
df$year <- 2008:2020
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2008,2020,1))+
  scale_y_continuous(name="Number of HIV Positive Adults", labels = scales::comma, breaks=seq(0,22521578,4000000))+
  theme_minimal(base_size = 22)




#------------------------------------
# Sensitivity Analysis
#------------------------------------


logkappa_sa  <- seq(0,6,0.5)



DF <- data.frame()

for (i in 1:length(logkappa_sa)) {
  
  run <- ode(times=times, y=start, func=sir,parms=logkappa_sa[i],maxsteps=100000)
  
  df <- as.data.frame(run[,2:7]+run[,8:13]+run[,14:19]+run[,20:25])
  names(df)  <- c("S", "L", "In", "It", "T", "R")
  df$year <- 2008:2020
  df$kappa <- as.integer(rep(round(exp(logkappa_sa[i])), nrow(df)))
  df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
  DF <- rbind(DF, df)
}



gp <- ggplot(DF, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous(breaks=seq(2008,2020,1))+
  scale_y_continuous(labels = scales::comma)+
  theme_minimal(base_size = 22)+
  ## gganimate functionality starts here
  labs(x = "Year", y = "Number of People", title = "kappa: {frame_time}") +
  transition_time(kappa) +
  ease_aes("linear")

animate(gp, height = 800, width =1000)
anim_save("kappa.gif")











#----------------------------------------
# MAKE PREDICTIONS : Fitted Sir Model
#----------------------------------------


## vector of timesteps
times_pred <- seq(0, 30, 1)      # 2020 to 2050
start_pred <- run_fit[13,-1]     # starting values are now 2020


run_pred <- ode(times=times_pred, y=start_pred, func=sir,parms=parm_fit,maxsteps=100000)
run_pred





# Child HIV Negative Population
df <- as.data.frame(run_pred[,2:7])
df$year <- 2020:2050
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2020,2050,5))+
  scale_y_continuous(name="Number of HIV Negative Children", labels = scales::comma)+
  theme_minimal(base_size = 22)


# Adult HIV Negative Population
df <- as.data.frame(run_pred[,8:13])
df$year <- 2020:2050
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2020,2050,5))+
  scale_y_continuous(name="Number of HIV Negative Adults", labels = scales::comma)+
  theme_minimal(base_size = 22)



# Child HIV Positive Population
df <- as.data.frame(run_pred[,14:19])
df$year <- 2020:2050
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2020,2050,5))+
  scale_y_continuous(name="Number of HIV Positive Children", labels = scales::comma)+
  theme_minimal(base_size = 22)


# Adult HIV Positive Population
df <- as.data.frame(run_pred[,20:25])
df$year <- 2020:2050
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2020,2050,5))+
  scale_y_continuous(name="Number of HIV Positive Adults", labels = scales::comma)+
  theme_minimal(base_size = 22)






#-----------------------------------------------------
# MAKE PREDICTIONS : Fitted Sir Model + Intervention
#-----------------------------------------------------


sir_int <- function(t, x, logkappa)  {
  with(as.list(c(logkappa, x)), {
    
    # Interaction Terms:
    kappa    = exp(logkappa)
    beta_cc  = 0.005*kappa   # not differentiating between contact between HIV neg and HIV pos people
    beta_ca  = 0.5*kappa
    beta_aa  = 1*kappa
    beta_ac  = 0.1*kappa
    
    # HIV-Negative Children
    p_nc     = 0.125      # proportion of children who develop active TB
    
    # HIV-Negative Adults
    p_na     = 0.145
    
    
    # HIV-Negative (children and adults)
    dtb_n   = 0.6      # death rate due to TB 
    d_n     = 11/1000  # death rate 
    x_n     = 0.33     # protection against development of reinfection or developing active TB due to being latently infected or recovered 
    r_n     = 0.01     # annual risk of relapse from recivered to active TB
    theta_n = 0.0002   # rate of progression to active disease due to existing infection 
    
    
    # HIV-Positive (children and adults)
    dtb_h   = 0.9       # death rate due to TB in HIV positive individuals 
    d_h     = (1-0.86) +  11/1000   # death rate for HIV positive individuals
    x_h     = 0.755     # protection against development of reinfection or developing active TB due to being latently infected or recovered fo HIV positive 
    r_h     = 0.425
    theta_h = 0.00236   # rate of progression to active disease due to existing infection 
    
    # HIV-Positive Adults
    p_ha    = 0.55
    
    # HIV-Positive Children
    p_hc    = 0.475
    mu      = 0.01    # proportion of hiv positive births
    
    # Whole Population
    t      = 0.58      # proportion of population who receive treatment 
    b      = 21.5/1000  # birth rate
    a      = 1/14       # aging rate 
    gammat = 0.764      # treatment recovery rate
    gamman = 0.32       # natural recovery rate
    s      = 60/365     # rate of treatment seeking 
    nu     = 0.185      # rate of HIV infection
    
    # Intervention
    omega = 0.85
    
    I_c = In_nc + It_nc + In_hc + It_hc
    I_a = In_na + It_na + In_ha + It_ha
    
    N_c = S_nc + L_nc + In_nc + It_nc + T_nc + R_nc + S_hc + L_hc + In_hc + It_hc + T_hc + R_hc
    N_a = S_na + L_na + In_na + It_na + T_na + R_na + S_ha + L_ha + In_ha + It_ha + T_ha + R_ha
    N = N_a + N_c
    
    lambda_c <- beta_cc*I_c/N_c + beta_ca*I_a/N_a
    lambda_a <- beta_aa*I_a/N_a + beta_ac*I_c/N_c
    
    
    # HIV Negative Children 
    dS_nc  = b*(1-mu)*N -S_nc*(a + d_n + lambda_c)
    dL_nc  = S_nc*lambda_c*(1-p_nc) + R_nc*x_n*lambda_c*(1-p_nc) - L_nc*(a + d_n + theta_n + lambda_c*p_nc*x_n)
    dIn_nc = S_nc*lambda_c*p_nc*(1-t) + L_nc*(theta_n + lambda_c*p_nc*x_n)*(1-t) + R_nc*(r_n + lambda_c*p_nc*x_n)*(1-t) - In_nc*(d_n + a + dtb_n + gamman)
    dIt_nc = S_nc*lambda_c*p_nc*t + L_nc*(theta_n + lambda_c*p_nc*x_n)*t + R_nc*(r_n + lambda_c*p_nc*x_n)*t - It_nc*(d_n + a + s)
    dT_nc  = It_nc*s - T_nc*(d_n + a + gammat)
    dR_nc  = T_nc*gammat + In_nc*gamman - R_nc*(d_n + a + r_n + x_n*lambda_c)
    
    
    # HIV Negative Adults
    dS_na  = (1-omega)*a*S_nc - S_na*(nu + d_n + lambda_a)
    dL_na  = (1-omega)*a*L_nc + S_na*lambda_a*(1-p_na) + R_na*x_n*lambda_a*(1-p_na) - L_na*(nu + d_n + theta_n + lambda_a*p_na*x_n)
    dIn_na = (1-omega)*a*In_nc + S_na*lambda_a*p_na*(1-t) + L_na*(theta_n + lambda_a*p_na*x_n)*(1-t) + R_na*(r_n + lambda_a*p_na*x_n)*(1-t) - In_na*(d_n + nu + dtb_n + gamman)
    dIt_na = (1-omega)*a*It_nc + S_na*lambda_a*p_na*t + L_na*(theta_n + lambda_a*p_na*x_n)*t + R_na*(r_n + lambda_a*p_na*x_n)*t - It_na*(d_n + nu + s)
    dT_na  = (1-omega)*a*T_nc + It_na*s - T_na*(d_n + nu + gammat)
    dR_na  = (1-omega)*a*R_nc + T_na*gammat + In_na*gamman - R_na*(d_n + nu + r_n + x_n*lambda_a)
    
    
    
    # HIV Positive Children 
    dS_hc  = b*mu*N - S_hc*(a + d_h + lambda_c)
    dL_hc  = S_hc*lambda_c*(1-p_hc) + R_hc*x_h*lambda_c*(1-p_hc) - L_hc*(a + d_h + theta_h + lambda_c*p_hc*x_h)
    dIn_hc = S_hc*lambda_c*p_hc*(1-t) + L_hc*(theta_h + lambda_c*p_hc*x_h)*(1-t) + R_hc*(r_h + lambda_c*p_hc*x_h)*(1-t) - In_hc*(d_h + a + dtb_h + gamman)
    dIt_hc = S_hc*lambda_c*p_hc*t + L_hc*(theta_h + lambda_c*p_hc*x_h)*t + R_hc*(r_h + lambda_c*p_hc*x_h)*t - It_hc*(d_h + a + s)
    dT_hc  = It_hc*s - T_hc*(d_h + a + gammat)
    dR_hc  = T_hc*gammat + In_hc*gamman - R_hc*(d_h + a + r_h + x_h*lambda_c)
    
    
    # HIV Positive Adults
    dS_ha  = a*S_hc + nu*S_na - S_ha*(d_h + lambda_a)
    dL_ha  = a*L_hc + nu*L_na + S_ha*lambda_a*(1-p_ha) + R_ha*x_h*lambda_a*(1-p_ha) - L_ha*(d_h + theta_h + lambda_a*p_ha*x_h)
    dIn_ha = a*In_hc + nu*In_na + S_ha*lambda_a*p_ha*(1-t) + L_ha*(theta_h + lambda_a*p_ha*x_h)*(1-t) + R_ha*(r_h + lambda_a*p_ha*x_h)*(1-t) -In_ha*(d_h + dtb_h + gamman)
    dIt_ha = a*It_hc + nu*It_na + S_ha*lambda_a*p_ha*t + L_ha*(theta_h + lambda_a*p_ha*x_n)*t + R_ha*(r_h + lambda_a*p_ha*x_h)*t -It_ha*(d_h + s)
    dT_ha  = a*T_hc + nu*T_na + It_ha*s - T_ha*(d_h + gammat)
    dR_ha  = a*R_hc + nu*R_na + T_ha*gammat + In_ha*gamman - R_ha*(d_h + r_h + x_h*lambda_a)
    
    
    
    output <- list(c(dS_nc, dL_nc, dIn_nc, dIt_nc, dT_nc, dR_nc, 
                     dS_na, dL_na, dIn_na, dIt_na, dT_na, dR_na,
                     dS_hc, dL_hc, dIn_hc, dIt_hc, dT_hc, dR_hc, 
                     dS_ha, dL_ha, dIn_ha, dIt_ha, dT_ha, dR_ha))
    output
  })
}




run_int <- ode(times=times_pred, y=start_pred, func=sir_int,parms=parm_fit,maxsteps=100000)
run_int





# Child HIV Negative Population
df <- as.data.frame(run_int[,2:7])
df$year <- 2020:2050
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2020,2050,5))+
  scale_y_continuous(name="Number of HIV Negative Children", labels = scales::comma)+
  theme_minimal(base_size = 22)


# Adult HIV Negative Population
df <- as.data.frame(run_int[,8:13])
df$year <- 2020:2050
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2020,2050,5))+
  scale_y_continuous(name="Number of HIV Negative Adults", labels = scales::comma)+
  theme_minimal(base_size = 22)



# Child HIV Positive Population
df <- as.data.frame(run_int[,14:19])
df$year <- 2020:2050
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2020,2050,5))+
  scale_y_continuous(name="Number of HIV Positive Children", labels = scales::comma)+
  theme_minimal(base_size = 22)


# Adult HIV Positive Population
df <- as.data.frame(run_int[,20:25])
df$year <- 2020:2050
df <- pivot_longer(df, cols = 1:6, names_to = "State", values_to = "Value")
df
ggplot(df, aes(year, Value, group = State, color = State)) +
  geom_line()+ 
  geom_point() +
  scale_x_continuous("Year", breaks=seq(2020,2050,5))+
  scale_y_continuous(name="Number of HIV Positive Adults", labels = scales::comma)+
  theme_minimal(base_size = 22)






