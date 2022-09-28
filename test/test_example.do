clear
set obs 2
gen W = mod(_n, 2)
expand 20
bys W: gen i = _n
gen X1 = 0
gen X2 = 0
bys W (i): replace X1 = 1 if i <= 1  & W == 0
bys W (i): replace X1 = 1 if i <= 9  & W == 1
bys W (i): replace X2 = 1 if i <= 10 & W == 0
bys W (i): replace X2 = 1 if i <= 18 & W == 1
replace X2 = 0 if X1 == 1
tab X1 if W == 0
tab X2 if W == 0
tab X1 if W == 1
tab X2 if W == 1
gen Y  = 0 * X1 + X2 * (1 - W)
reg Y X1 X2 W
gen byte X = 1 * X1 + 2 * X2
multe Y X, control(W)
reg Y X1 W
reg Y X2 W
