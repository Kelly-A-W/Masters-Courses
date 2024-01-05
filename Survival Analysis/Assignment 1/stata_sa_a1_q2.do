clear
cd "/Users/kellywilliams/Documents/Uni Courses/Masters/Survival Analysis/Assignment 1"
import delimited "larynx.csv"



* create survival object
stset time, failure(delta)



* Weibul PH with constant shape parameter
streg i.stage, dist(weibull)

streg i.stage, dist(weibull) ancillary(i.stage)
