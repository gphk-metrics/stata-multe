capture program drop multe_weight_tests
program multe_weight_tests
    multe_load_test_data
    gen `c(obs_t)' G = _n
    gen `c(obs_t)' _expand = ceil(runiform() * _N/50)
    gen double Y = Ydouble + rnormal() * _N / 50
    qui sum _expand
    local ratio = sqrt(r(sum) / _N)

    * --------------------
    * Internal Consistency
    * --------------------

    qui multe Y Cdouble             , treat(Tbyte) strat(Wbyte) mata(now)
    qui multe Y Cdouble [aw=_expand], treat(Tbyte) strat(Wbyte) mata(aw)
    qui multe Y Cdouble [pw=_expand], treat(Tbyte) strat(Wbyte) mata(pw)

    mata assert(max(reldif(aw.full.estA, now.full.estA)) > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.estB, now.full.estB)) > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.seP,  now.full.seP))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.seB,  now.full.seB))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.seO,  now.full.seO))  > (epsilon(1))^(1/4))

    mata assert(max(reldif(pw.full.estA, now.full.estA)) > (epsilon(1))^(1/4))
    mata assert(max(reldif(pw.full.estB, now.full.estB)) > (epsilon(1))^(1/4))
    mata assert(max(reldif(pw.full.seP,  now.full.seP))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(pw.full.seB,  now.full.seB))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(pw.full.seO,  now.full.seO))  > (epsilon(1))^(1/4))

    mata assert(max(reldif(pw.full.estA, aw.full.estA)) < (epsilon(1))^(3/4))
    mata assert(max(reldif(pw.full.estB, aw.full.estB)) < (epsilon(1))^(3/4))
    mata assert(max(reldif(pw.full.seP,  aw.full.seP))  < (epsilon(1))^(3/4))
    mata assert(max(reldif(pw.full.seB,  aw.full.seB))  < (epsilon(1))^(3/4))
    mata assert(max(reldif(pw.full.seO,  aw.full.seO))  < (epsilon(1))^(3/4))

    disp "(multe test success): multe basic weights consistency on simulated data"

    * NB: Fundamentally aw and fw have different variances. The first
    * gives V((X' W X)^-1 X' W e) = (X' W X)^-1 X' W V(e) W X (X' W X)^-1
    * and the second gives (X' W X)^-1 X' (W V(e)) X (X' W X)^-1

    * ----------------------
    * vs teffects (ate only)
    * ----------------------

    qui tab Wbyte, gen(_W)
    local k = r(r) + 1

    qui teffects ipw (Y) (Tbyte b0.Wbyte)
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(now.full.estA[.,3], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(now.full.seP[.,3],  table[2,.]')) < (epsilon(1)^(3/4)))

    qui teffects ipw (Y) (Tbyte b0.Wbyte) [pw=_expand]
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(pw.full.estA[.,3], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(pw.full.seP[.,3],  table[2,.]')) < (epsilon(1)^(3/4)))

    cap drop _aw
    qui sum _expand
    gen double _aw = _N * _expand / r(sum)
    qui teffects ipw (Y) (Tbyte b0.Wbyte) [pw=_aw]
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(aw.full.estA[.,3], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(aw.full.seP[.,3],  table[2,.]')) < (epsilon(1)^(3/4)))

    disp "(multe test success): multe ATE consistency on simulated data"

    * ----------------------
    * vs regress (oaat only)
    * ----------------------

    mata oaat = J(4, 2, .)
    forvalues i = 2 / 5 {
        qui reg Y i`i'.Tbyte _W* if inlist(Tbyte, 1, `i'), noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "Tbyte"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(now.full.estA[., 4], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (now.full.seP[., 4] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    mata oaat = J(4, 2, .)
    forvalues i = 2 / 5 {
        qui reg Y i`i'.Tbyte _W* if inlist(Tbyte, 1, `i') [pw=_expand], noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "Tbyte"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(pw.full.estA[., 4], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (pw.full.seP[., 4] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    mata oaat = J(4, 2, .)
    forvalues i = 2 / 5 {
        qui reg Y i`i'.Tbyte _W* if inlist(Tbyte, 1, `i') [aw=_expand], noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "Tbyte"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(aw.full.estA[., 4], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (aw.full.seP[., 4] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    disp "(multe test success): multe OAAT consistency on simulated data"

    * --------------------
    * Internal Consistency
    * --------------------

    qui expand _expand
    qui multe Y Cdouble, treat(Tbyte) strat(Wbyte) mata(ex)

    * mata assert(max(reldif(fw.full.estA, ex.full.estA)) < (epsilon(1))^(3/4))
    * mata assert(max(reldif(fw.full.estB, ex.full.estB)) < (epsilon(1))^(3/4))
    * mata assert(max(reldif(fw.full.seP,  ex.full.seP))  < (epsilon(1))^(3/4))
    * mata assert(max(reldif(fw.full.seB,  ex.full.seB))  < (epsilon(1))^(3/4))
    * mata assert(max(reldif(fw.full.seO,  ex.full.seO))  < (epsilon(1))^(3/4))

    mata assert(max(reldif(aw.full.estA, ex.full.estA)) < (epsilon(1))^(3/4))
    mata assert(max(reldif(aw.full.estB, ex.full.estB)) < (epsilon(1))^(3/4))
    mata assert(max(reldif(aw.full.seP,  ex.full.seP))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.seB,  ex.full.seB))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.seO,  ex.full.seO))  > (epsilon(1))^(1/4))

    * NB: Fundamentally aw and fw have different variances. The first
    * gives V((X' W X)^-1 X' W e) = (X' W X)^-1 X' W V(e) W X (X' W X)^-1
    * and the second gives (X' W X)^-1 X' (W V(e)) X (X' W X)^-1

    * --------
    * External
    * --------

    mata oaat = J(4, 2, .)
    forvalues i = 2 / 5 {
        qui reg Y i`i'.Tbyte _W* if inlist(Tbyte, 1, `i'), noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "Tbyte"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(ex.full.estA[., 4], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (ex.full.seP[., 4] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    qui teffects ipw (Y) (Tbyte b0.Wbyte)
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(ex.full.estA[.,3], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(ex.full.seP[.,3],  table[2,.]')) < (epsilon(1)^(3/4)))

    * -------------
    * More internal
    * -------------

    contract Y Cdouble Tbyte Wbyte G _expand, f(_nobs)
    assert _expand == _nobs
    gisid G
    qui multe Y Cdouble [aw=_expand], treat(Tbyte) strat(Wbyte) mata(con)

    * mata assert(max(reldif(fw.full.estA, con.full.estA)) < (epsilon(1))^(3/4))
    * mata assert(max(reldif(fw.full.estB, con.full.estB)) < (epsilon(1))^(3/4))
    * mata assert(max(reldif(fw.full.seP,  con.full.seP))  > (epsilon(1))^(1/4))
    * mata assert(max(reldif(fw.full.seB,  con.full.seB))  > (epsilon(1))^(1/4))
    * mata assert(max(reldif(fw.full.seO,  con.full.seO))  > (epsilon(1))^(1/4))

    mata assert(max(reldif(aw.full.estA, con.full.estA)) < (epsilon(1))^(3/4))
    mata assert(max(reldif(aw.full.estB, con.full.estB)) < (epsilon(1))^(3/4))
    mata assert(max(reldif(aw.full.seP,  con.full.seP))  < (epsilon(1))^(3/4))
    mata assert(max(reldif(aw.full.seB,  con.full.seB))  < (epsilon(1))^(3/4))
    mata assert(max(reldif(aw.full.seO,  con.full.seO))  < (epsilon(1))^(3/4))

    disp "(multe test success): multe internal expand/contract consistency on simulated data"
end

capture program drop multe_weight_startest
program multe_weight_startest
    set seed 1729
    use test/example_star.dta, clear
    gen `c(obs_t)' G = _n
    gen `c(obs_t)' _expand = ceil(runiform() * 10)

    qui multe score             , treat(treatment) strat(school) mata(now)
    qui multe score [aw=_expand], treat(treatment) strat(school) mata(aw)
    qui multe score [pw=_expand], treat(treatment) strat(school) mata(pw)

    * ---------------
    * Internal Checks
    * ---------------

    mata assert(max(reldif(aw.full.estA, now.full.estA)) > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.estB, now.full.estB)) > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.seP,  now.full.seP))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.seB,  now.full.seB))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.seO,  now.full.seO))  > (epsilon(1))^(1/4))

    mata assert(max(reldif(pw.full.estA, now.full.estA)) > (epsilon(1))^(1/4))
    mata assert(max(reldif(pw.full.estB, now.full.estB)) > (epsilon(1))^(1/4))
    mata assert(max(reldif(pw.full.seP,  now.full.seP))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(pw.full.seB,  now.full.seB))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(pw.full.seO,  now.full.seO))  > (epsilon(1))^(1/4))

    mata assert(max(reldif(pw.full.estA, aw.full.estA)) < (epsilon(1))^(3/4))
    mata assert(max(reldif(pw.full.estB, aw.full.estB)) < (epsilon(1))^(3/4))
    mata assert(max(reldif(pw.full.seP,  aw.full.seP))  < (epsilon(1))^(3/4))
    mata assert(max(reldif(pw.full.seB,  aw.full.seB))  < (epsilon(1))^(3/4))
    mata assert(max(reldif(pw.full.seO,  aw.full.seO))  < (epsilon(1))^(3/4))

    disp "(multe test success): multe basic weights consistency on simulated data"

    * ---------------------
    * vs teffects ipw (ATE)
    * ---------------------

    qui tab school, gen(_W)
    local k = r(r) + 1

    qui teffects ipw (score) (treatment i.school)
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(now.full.estA[.,3], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(now.full.seP[.,3],  table[2,.]')) < (epsilon(1)^(3/4)))

    qui teffects ipw (score) (treatment i.school) [pw=_expand]
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(pw.full.estA[.,3], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(pw.full.seP[.,3],  table[2,.]')) < (epsilon(1)^(3/4)))

    cap drop _aw
    qui sum _expand
    gen double _aw = _N * _expand / r(sum)
    qui teffects ipw (score) (treatment i.school) [pw=_aw]
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(aw.full.estA[.,3], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(aw.full.seP[.,3],  table[2,.]')) < (epsilon(1)^(3/4)))

    disp "(multe test success): multe ATE consistency on simulated data"

    * -----------------
    * vs regress (OAAT)
    * -----------------

    mata oaat = J(2, 2, .)
    forvalues i = 2 / 3 {
        qui reg score i`i'.treatment _W* if inlist(treatment, 1, `i'), noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "treatment"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(now.full.estA[., 4], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (now.full.seP[., 4] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    mata oaat = J(2, 2, .)
    forvalues i = 2 / 3 {
        qui reg score i`i'.treatment _W* if inlist(treatment, 1, `i') [pw=_expand], noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "treatment"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(pw.full.estA[., 4], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (pw.full.seP[., 4] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    mata oaat = J(2, 2, .)
    forvalues i = 2 / 3 {
        qui reg score i`i'.treatment _W* if inlist(treatment, 1, `i') [aw=_expand], noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "treatment"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(aw.full.estA[., 4], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (aw.full.seP[., 4] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    disp "(multe test success): multe OAAT consistency on simulated data"

    * ----------------------
    * Check in expanded data
    * ----------------------

    qui expand _expand
    qui multe score, treat(treatment) strat(school) mata(ex)

    * mata assert(max(reldif(fw.full.estA, ex.full.estA)) < (epsilon(1))^(3/4))
    * mata assert(max(reldif(fw.full.estB, ex.full.estB)) < (epsilon(1))^(3/4))
    * mata assert(max(reldif(fw.full.seP,  ex.full.seP))  < (epsilon(1))^(3/4))
    * mata assert(max(reldif(fw.full.seB,  ex.full.seB))  < (epsilon(1))^(3/4))
    * mata assert(max(reldif(fw.full.seO,  ex.full.seO))  < (epsilon(1))^(3/4))

    mata assert(max(reldif(aw.full.estA, ex.full.estA)) < (epsilon(1))^(3/4))
    mata assert(max(reldif(aw.full.estB, ex.full.estB)) < (epsilon(1))^(3/4))
    mata assert(max(reldif(aw.full.seP,  ex.full.seP))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.seB,  ex.full.seB))  > (epsilon(1))^(1/4))
    mata assert(max(reldif(aw.full.seO,  ex.full.seO))  > (epsilon(1))^(1/4))

    mata oaat = J(2, 2, .)
    forvalues i = 2 / 3 {
        qui reg score i`i'.treatment _W* if inlist(treatment, 1, `i'), noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "treatment"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(ex.full.estA[., 4], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (ex.full.seP[., 4] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    qui teffects ipw (score) (treatment i.school)
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(ex.full.estA[.,3], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(ex.full.seP[.,3],  table[2,.]')) < (epsilon(1)^(3/4)))

    disp "(multe test success): multe internal expand consistency on STAR data"
end
