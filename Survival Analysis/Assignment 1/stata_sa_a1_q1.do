clear
cd "/Users/kellywilliams/Documents/Uni Courses/Masters/Survival Analysis/Assignment 1"
import delimited "larynx.csv"


* look at the data 
describe
tab stage
summarize stage time age diagyr delta survrel


* create survival object
stset time, failure(delta)


* ----------------------------- 1(a) --------------------------------------
sts list, by(stage) compare
sts test stage, logrank
sts graph, by(stage) risktable(0 1 2 3 4 5 6 7 8 9 10, failevents) tmax(11) xlabel(0 1 2 3 4 5 6 7 8 9 10)



* ----------------------------- 1(b) --------------------------------------
* fit a cox model with age and stage
stcox age i.stage

* Cox-Snell
predict cs, csnell
quietly stset cs, failure(delta)
quietly sts generate km=s
quietly generate H=-ln(km)
line H cs cs, sort

* Martingale
quietly  stcox age i.stage, mgale(mg1) nolog
lowess mg1 age, mean noweight title("") note("")

* Deviance Residuals
quietly  stcox age i.stage, mgale(mg) 
quietly predict xb,xb
generate obs = _n
twoway scatter mg xb, mlabel(obs) 

* Score Residuals
quietly stcox age i.stage
predict inf*, dfbeta
scatter inf* time, yline(0) mlabel(obs)

* Schoenfeld Residuals 
stcox age i.stage, scaledsch(sca*) schoenfeld(sch*) nolog
stphtest, rank detail


* ----------------------------- 1(c) --------------------------------------
stcox age i.stage
stcurve, survival at1(age=50 stage=1) at2(age=50 stage=4)


* ----------------------------- 1(d) --------------------------------------
gen age50=age-50


stcox age50 i.stage, basesurv(basesurv501)
stcox age50 b4.stage, basesurv(basesurv504)

sort time

list time basesurv501 basesurv504
