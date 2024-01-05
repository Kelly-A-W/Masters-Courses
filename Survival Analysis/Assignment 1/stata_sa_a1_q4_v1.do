clear
cd "/Users/kellywilliams/Documents/Uni Courses/Masters/Survival Analysis/Assignment 1"
import delimited "larynx.csv"


* look at the data 
describe
tab stage
summarize stage time age diagyr delta survrel


* create survival object
stset time, failure(delta)


* Fit the exponential, lognormal and gamma AFT models including the variables age and stage as
* covariates. 
* Plot the hazard function by stage for each of these models. 
* Check the overall fits of the models. 
* Which model (including the Weibull AFT) is the most suitable for these data? 
* Compare the decceleration factors between these models. - see slide 31

* Exponential AFT
* Fit model
streg age i.stage, dist(exponential) time  
* Plot the hazard function by stage
stcurve, hazard at(stage = (1 2 3 4)) ylabel(, format(%9.1f)) title("")
* Check the overall fits 
predict double cs, csnell
quietly stset cs, failure(delta)
quietly sts generate km=s
quietly generate double H=-ln(km)
quietly line H cs cs, sort
estat ic

* Lognormal AFT
streg age i.stage, dist(lognormal)  
stcurve, hazard at(stage = (1 2 3 4)) ylabel(, format(%9.1f)) title("")
predict double cs, csnell
quietly stset cs, failure(delta)
quietly sts generate km=s
quietly generate double H=-ln(km)
quietly line H cs cs, sort
estat ic


* gamma AFT
streg age i.stage, dist(gamma) 
stcurve, hazard at(stage = (1 2 3 4)) ylabel(, format(%9.1f)) title("")
predict double cs, csnell
quietly stset cs, failure(delta)
quietly sts generate km=s
quietly generate double H=-ln(km)
quietly line H cs cs, sort

