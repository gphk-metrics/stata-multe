local nobs   50000
local ktreat 10
clear
mata mata clear
qui do src/ado/multe.ado
qui do src/mata/multe_helpers.mata
qui do src/mata/multe.mata
set seed 1729
set obs `nobs'
gen T = ceil(runiform() * `ktreat')
gen W = mod(_n, 100)
gen Y = T + runiform()
multe Y T, control(W)
