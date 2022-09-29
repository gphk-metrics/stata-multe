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

    multe Y Tbyte             , control(Wbyte) mata(now)  decomp gen(lambda(lnow) tau(tnow))
    multe Y Tbyte [aw=_expand], control(Wbyte) mata(aw)   decomp gen(lambda(law)  tau(taw))
    multe Y Tbyte [fw=_expand], control(Wbyte) mata(fw)   decomp gen(lambda(lfw)  tau(tfw))

    unab taus: tnow?
    forvalues i = 1 / `:list sizeof taus' {
        assert reldif(taw`i', tfw`i')  < 1e-12
        assert reldif(taw`i', tnow`i') > 1e-6
        forvalues j = 1 / `:list sizeof taus' {
            assert reldif(law`i'`j', lfw`i'`j') < 1e-12
            assert reldif(law`i'`j', lnow`i'`j') > 1e-6
        }
    }

    mata assert(max(reldif(aw.estimates.est,   now.estimates.est))   > sqrt(epsilon(1)))
    mata assert(max(reldif(aw.estimates.se_po, now.estimates.se_po)) > sqrt(epsilon(1)))
    mata assert(max(reldif(aw.estimates.se_or, now.estimates.se_or)) > sqrt(epsilon(1)))
    mata assert(max(reldif(fw.estimates.se_po, now.estimates.se_po)) > sqrt(epsilon(1)))
    mata assert(max(reldif(fw.estimates.se_or, now.estimates.se_or)) > sqrt(epsilon(1)))

    mata assert(max(reldif(aw.estimates.est, fw.estimates.est)) < (epsilon(1)^(3/4)))
    * NB: Fundamentally aw and fw have different variances. The first
    * gives V((X' W X)^-1 X' W e) = (X' W X)^-1 X' W V(e) W X (X' W X)^-1
    * and the second gives (X' W X)^-1 X' (W V(e)) X (X' W X)^-1
    *
    * mata assert(max(abs(`ratio' :- aw.estimates.se_po[., (1, 3)] :/ fw.estimates.se_po[., (1, 3)])) < (epsilon(1)^(3/4)))
    * mata assert(max(abs(`ratio' :- aw.estimates.se_or[., (1, 3)] :/ fw.estimates.se_or[., (1, 3)])) < (epsilon(1)^(3/4)))

    * ----------------------
    * vs teffects (ate only)
    * ----------------------

    qui tab Wbyte, gen(_W)
    local k = r(r) + 1

    teffects ipw (Y) (Tbyte b0.Wbyte)
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(now.estimates.est[.,1], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(now.estimates.se_po[.,1], table[2,.]')) < (epsilon(1)^(3/4)))

    teffects ipw (Y) (Tbyte b0.Wbyte) [fw=_expand]
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(fw.estimates.est[.,1], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(fw.estimates.se_po[.,1], table[2,.]')) < (epsilon(1)^(3/4)))

    * NB: pw not allowed in multe so vce cannot be compared
    teffects ipw (Y) (Tbyte b0.Wbyte) [pw=_expand]
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(aw.estimates.est[.,1], table[1,.]')) < (epsilon(1)^(3/4)))

    cap drop _aw
    qui sum _expand
    gen double _aw = _N * _expand / r(sum)
    teffects ipw (Y) (Tbyte b0.Wbyte) [pw=_aw]
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(aw.estimates.est[.,1], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(aw.estimates.se_po[.,1], table[2,.]')) < (epsilon(1)^(3/4)))

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
    mata assert(max(reldif(now.estimates.est[., 2], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (now.estimates.se_po[., 2] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    mata oaat = J(4, 2, .)
    forvalues i = 2 / 5 {
        qui reg Y i`i'.Tbyte _W* if inlist(Tbyte, 1, `i') [aw=_expand], noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "Tbyte"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(aw.estimates.est[., 2], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (aw.estimates.se_po[., 2] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    mata oaat = J(4, 2, .)
    forvalues i = 2 / 5 {
        qui reg Y i`i'.Tbyte _W* if inlist(Tbyte, 1, `i') [fw=_expand], noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata X    = sort(st_data(., "Tbyte _expand"), 1)
    mata info = panelsetup(X[., 1], 1)
    mata nj   = panelsum(X[., 2], info)
    mata nj   = nj[1] :+ nj[|2 \ length(nj)|]
    mata assert(max(reldif(fw.estimates.est[., 2], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (fw.estimates.se_po[., 2] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    * --------------------
    * Internal Consistency
    * --------------------

    expand _expand
    multe Y Tbyte, control(Wbyte) mata(ex) decomp gen(lambda(lex) tau(tex))

    unab taus: tex?
    forvalues i = 1 / `:list sizeof taus' {
        assert reldif(tex`i', taw`i')  < 1e-12
        assert reldif(tex`i', tfw`i')  < 1e-12
        assert reldif(tex`i', tnow`i') > 1e-6
        forvalues j = 1 / `:list sizeof taus' {
            assert reldif(lex`i'`j', law`i'`j') < 1e-12
            assert reldif(lex`i'`j', lfw`i'`j') < 1e-12
            assert reldif(lex`i'`j', lnow`i'`j') > 1e-6
        }
    }

    mata assert(max(reldif(fw.estimates.est,   ex.estimates.est))   < (epsilon(1)^(3/4)))
    mata assert(max(reldif(fw.estimates.se_po, ex.estimates.se_po)) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(fw.estimates.se_or, ex.estimates.se_or)) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(aw.estimates.est,   ex.estimates.est))   < (epsilon(1)^(3/4)))
    * NB: Fundamentally aw and fw have different variances. The first
    * gives V((X' W X)^-1 X' W e) = (X' W X)^-1 X' W V(e) W X (X' W X)^-1
    * and the second gives (X' W X)^-1 X' (W V(e)) X (X' W X)^-1
    *
    * mata assert(max(abs(`ratio' :- aw.estimates.se_po[., (1, 3)] :/ ex.estimates.se_po[., (1, 3)])) < (epsilon(1)^(3/4)))
    * mata assert(max(abs(`ratio' :- aw.estimates.se_or[., (1, 3)] :/ ex.estimates.se_or[., (1, 3)])) < (epsilon(1)^(3/4)))

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
    mata assert(max(reldif(ex.estimates.est[., 2], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (ex.estimates.se_po[., 2] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    teffects ipw (Y) (Tbyte b0.Wbyte)
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(ex.estimates.est[.,1], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(ex.estimates.se_po[.,1], table[2,.]')) < (epsilon(1)^(3/4)))

    * -------------
    * More internal
    * -------------

    contract Y Tbyte Wbyte G _expand lex* tex*, f(_nobs)
    assert _expand == _nobs
    gisid G
    multe Y Tbyte [fw=_expand], control(Wbyte) mata(con) decomp gen(lambda(lcon) tau(tcon))

    unab taus: tcon?
    forvalues i = 1 / `:list sizeof taus' {
        assert reldif(tex`i', tcon`i') < 1e-12
        forvalues j = 1 / `:list sizeof taus' {
            assert reldif(lex`i'`j', lcon`i'`j') < 1e-12
        }
    }

    mata assert(max(reldif(fw.estimates.est,   con.estimates.est))   < (epsilon(1)^(3/4)))
    mata assert(max(reldif(fw.estimates.se_po, con.estimates.se_po)) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(fw.estimates.se_or, con.estimates.se_or)) < (epsilon(1)^(3/4)))
end

capture program drop multe_weight_startest
program multe_weight_startest
    set seed 1729
    use test/example_star.dta, clear
    gen `c(obs_t)' G = _n
    gen `c(obs_t)' _expand = ceil(runiform() * _N/50)

    multe score treatment             , control(school) mata(now)  decomp gen(lambda(lnow) tau(tnow))
    multe score treatment [aw=_expand], control(school) mata(aw)   decomp gen(lambda(law)  tau(taw))
    multe score treatment [fw=_expand], control(school) mata(fw)   decomp gen(lambda(lfw)  tau(tfw))

    * ---------------
    * Internal Checks
    * ---------------

    unab taus: tnow?
    forvalues i = 1 / `:list sizeof taus' {
        assert reldif(taw`i', tfw`i')  < 1e-12
        forvalues j = 1 / `:list sizeof taus' {
            assert reldif(law`i'`j', lfw`i'`j') < 1e-12
        }
    }
    mata assert(max(reldif(aw.estimates.est, fw.estimates.est)) < (epsilon(1)^(3/4)))

    * ---------------------
    * vs teffects ipw (ATE)
    * ---------------------

    qui tab school, gen(_W)
    local k = r(r) + 1

    teffects ipw (score) (treatment i.school)
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(now.estimates.est[.,1], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(now.estimates.se_po[.,1], table[2,.]')) < (epsilon(1)^(3/4)))

    teffects ipw (score) (treatment i.school) [fw=_expand]
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(fw.estimates.est[.,1], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(fw.estimates.se_po[.,1], table[2,.]')) < (epsilon(1)^(3/4)))

    * NB: pw not allowed in multe so vce cannot be compared
    teffects ipw (score) (treatment i.school) [pw=_expand]
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(aw.estimates.est[.,1], table[1,.]')) < (epsilon(1)^(3/4)))

    cap drop _aw
    qui sum _expand
    gen double _aw = _N * _expand / r(sum)
    teffects ipw (score) (treatment i.school) [pw=_aw]
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(aw.estimates.est[.,1], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(aw.estimates.se_po[.,1], table[2,.]')) < (epsilon(1)^(3/4)))

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
    mata assert(max(reldif(now.estimates.est[., 2], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (now.estimates.se_po[., 2] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    mata oaat = J(2, 2, .)
    forvalues i = 2 / 3 {
        qui reg score i`i'.treatment _W* if inlist(treatment, 1, `i') [aw=_expand], noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "treatment"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(aw.estimates.est[., 2], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (aw.estimates.se_po[., 2] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    mata oaat = J(2, 2, .)
    forvalues i = 2 / 3 {
        qui reg score i`i'.treatment _W* if inlist(treatment, 1, `i') [fw=_expand], noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata X    = sort(st_data(., "treatment _expand"), 1)
    mata info = panelsetup(X[., 1], 1)
    mata nj   = panelsum(X[., 2], info)
    mata nj   = nj[1] :+ nj[|2 \ length(nj)|]
    mata assert(max(reldif(fw.estimates.est[., 2], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (fw.estimates.se_po[., 2] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    * ----------------------
    * Check in expanded data
    * ----------------------

    expand _expand
    multe score treatment, control(school) mata(ex) decomp gen(lambda(lex) tau(tex))

    unab taus: tex?
    forvalues i = 1 / `:list sizeof taus' {
        assert reldif(tex`i', taw`i')  < 1e-12
        assert reldif(tex`i', tfw`i')  < 1e-12
        forvalues j = 1 / `:list sizeof taus' {
            assert reldif(lex`i'`j', law`i'`j') < 1e-12
            assert reldif(lex`i'`j', lfw`i'`j') < 1e-12
        }
    }
    mata assert(max(reldif(fw.estimates.est,   ex.estimates.est))   < (epsilon(1)^(3/4)))
    mata assert(max(reldif(fw.estimates.se_po, ex.estimates.se_po)) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(fw.estimates.se_or, ex.estimates.se_or)) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(aw.estimates.est,   ex.estimates.est))   < (epsilon(1)^(3/4)))

    mata oaat = J(2, 2, .)
    forvalues i = 2 / 3 {
        qui reg score i`i'.treatment _W* if inlist(treatment, 1, `i'), noc r
        mata oaat[`i'-1, .] = st_matrix("r(table)")[1::2, 1]'
    }
    mata info = panelsetup(sort(st_data(., "treatment"), 1), 1)
    mata nj   = info[2::rows(info), 2] :- info[2::rows(info), 1] :+ 1 :+ info[1, 2]
    mata assert(max(reldif(ex.estimates.est[., 2], oaat[., 1])) < (epsilon(1)^(3/4)))
    mata assert(max(reldif((nj :- `k') :/ nj, (ex.estimates.se_po[., 2] :/ oaat[., 2]):^2)) < (epsilon(1)^(3/4)))

    teffects ipw (score) (treatment i.school)
    matrix table = r(table)
    matrix table = table[1..2,"ATE:"]
    mata table = st_matrix("table")
    mata assert(max(reldif(ex.estimates.est[.,1], table[1,.]')) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(ex.estimates.se_po[.,1], table[2,.]')) < (epsilon(1)^(3/4)))

end
