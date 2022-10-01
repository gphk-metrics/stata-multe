* ------------
* Misc example
* ------------
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

* ---------------
* Weights example
* ---------------
sysuse auto, clear
mata {
    Y    = st_data(., "price")
    X    = st_data(., "turn"), J(st_nobs(), 1, 1)
    w    = st_data(., "mpg")
    AW   = diag(st_nobs() * w / sum(w))
    FW   = diag(w)
    b    = invsym(X' * AW * X) * X' * AW * Y
    reldif(b, invsym(X' * FW * X) * X' * FW * Y)
    e    = Y - X * b
    Vaw  = invsym(X' * AW * X) * X' * AW * diag(e:^2) * AW * X * invsym(X' * AW * X)
    Vfw  = invsym(X' * FW * X) * X' * FW * diag(e:^2) * X * invsym(X' * FW * X)
    qcaw = st_nobs() / (st_nobs() - 2)
    qcfw = sum(w) / (sum(w) - 2)
    seaw = sqrt(diagonal(Vaw) * qcaw)
    sefw = sqrt(diagonal(Vfw) * qcfw)
}

reg price turn [aw = mpg], r
mata reldif(seaw, st_matrix("r(table)")[2, .]')
reg price turn [fw = mpg], r
mata reldif(sefw, st_matrix("r(table)")[2, .]')
mata reldif(sefw, seaw)

expand mpg
reg price turn, r
mata reldif(sefw, st_matrix("r(table)")[2, .]')
mata reldif(seaw, st_matrix("r(table)")[2, .]')
