clear
cd "/Users/kellywilliams/Documents/Uni Courses/Masters/Survival Analysis/Assignment 1"
import delimited "larynx.csv"


* look at the data 
describe
tab stage
summarize stage time age diagyr delta survrel


* create survival object
stset time, failure(delta)


* fit a cox model with age and stage
stcox age i.stage



gen age50=age-50


survci, at(age=50 stage=1) outfile(age50stage1)




survci, at(age=50 stage=4) outfile(age50stage4)



