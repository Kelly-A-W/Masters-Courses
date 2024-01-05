clear
cd "/Users/kellywilliams/Documents/Uni Courses/Masters/Survival Analysis/Assignment 1"
import delimited "larynx.csv"


* create survival object
stset time, failure(delta)

streg age i.stage, dist(weibull) time  

estat ic
