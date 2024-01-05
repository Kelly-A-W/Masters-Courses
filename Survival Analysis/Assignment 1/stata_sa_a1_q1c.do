clear
cd "/Users/kellywilliams/Documents/Uni Courses/Masters/Survival Analysis/Assignment 1"
import delimited "larynx.csv"


* look at the data 
describe
tab stage
summarize stage time age diagyr delta survrel


* create survival object
stset time, failure(delta)


stcox age i.stage
stcurve, survival at1(age=50 stage=1) at2(age=50 stage=4)
