clear
cd "/Users/kellywilliams/Documents/Uni Courses/Masters/Survival Analysis/Assignment 1"
import delimited "larynx.csv"


* look at the data 
describe
tab stage
summarize stage time age diagyr delta survrel


* create survival object
stset time, failure(delta)


* Fit the generalized Gamma AFT model including the variables age and stage as covariates. 
* Check the overall fit of the model. 
* Use the fitted model to assess the suitability of the models fitted in questions 3 and 4.

* Fit model:
streg i.stage age, dist(ggamma)
* store generalised gamma model
estimates store sat  

* Check the overall fit
predict double cs, csnell
quietly stset cs, failure(delta)
quietly sts generate km=s
quietly generate double H=-ln(km)
quietly line H cs cs, sort


* Assess the suitability of the models fitted in questions 3 and 4

* Wald Test: 

* Weibull
test [/kappa] = 1.           
* Ho (Weibull) not rejected

* log-normal
test [/kappa] = 0            
* Ho (log-normal) not rejected

* Gamma
test [/kappa] = 1.103853     
* Ho (Gamma) not rejected

* Exponential
test [/kappa] = 1, notest
test [/lnsigma] = 0, accum  
* Ho (exponential) not rejected


* LRT Test:
* Weibull
streg age i.stage, dist(weibull)  
lrtest sat, force
* log-normal
streg age i.stage, dist(lognormal)  
lrtest sat, force
* Gamma
streg age i.stage, dist(weibull)  
lrtest sat, force
* Exponential
streg age i.stage, dist(exponential) time  
lrtest sat, force








