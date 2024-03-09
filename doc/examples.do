* end
* exit, clear
* stata14-mp

* Data from Fryer and Levitt (2013)
use test/example_fryer_levitt.dta, clear

* Regress IQ at 24 months on race indicators and baseline controls
multe std_iq_24 i.age_24 female [w=W2C0], treat(race) stratum(SES_quintile)

* Post ATE estimates
multe, est(ATE)

* Add covariates that fail overlap
local controls i.age_24 female i.siblings i.family_structure
multe std_iq_24 `controls' [w=W2C0], treat(race) stratum(SES_quintile)
tab race if siblings == 6, mi

* Show collected results for differences
matlist e(diffmatrix), format(%7.4g)
matlist e(overlapdiffmatrix), format(%7.4g)

* Post OWN estimates
multe, est(OWN) diff
multe, est(OWN) diff overlap

* Oracle SE
multe, est(ATE) oracle

* Cluster SE
local controls i.age_24 female
local options  treat(race) stratum(SES_quintile)
multe std_iq_24 `controls' [w=W2C0], `options' cluster(interviewer_ID_24)

* Example from help file
local nobs   1000
local ktreat 5
clear
set seed 1729
set obs `nobs'
gen T = ceil(runiform() * `ktreat')
gen W = mod(_n, 10)
gen Y = T + runiform()
multe Y, treat(T) strat(W)
ereturn list
