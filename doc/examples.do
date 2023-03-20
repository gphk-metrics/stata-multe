* Generated data
local nobs   1000
local ktreat 5
clear
set seed 1729
set obs `nobs'
gen T = ceil(runiform() * `ktreat')
gen W = mod(_n, 10)
gen Y = T + runiform()
multe Y T, control(W)
ereturn list

* Use cached results to compute decomposition, lambda, tau
multe, vce(oracle)
multe, decomposition
multe, decomposition minmax
multe, gen(lambda tau)

* Compute decomposition, lambda, tau from the onset
multe Y T, control(W) decomp gen(lambda(awesomeName) tau(coolerName))
desc, full

* Project STAR
use test/example_star.dta, clear
multe score treatment, control(school)
ereturn list
multe, vce(oracle)
multe, gen(lambda(M_) tau(tauhat_))
multe, decomp
desc, full
corr tauhat_? M_??

multe score treatment, control(school) matasave(matastructname)
mata mata desc
