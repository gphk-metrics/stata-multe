capture program drop multe_weight_tests
program multe_weight_tests
    multe_load_test_data
    gen `c(obs_t)' G = _n
    gen `c(obs_t)' _expand = ceil(runiform() * _N/50)
    gen double Y = Ydouble + rnormal() * _N / 50

    multe Y Tbyte             , control(Wbyte) mata(now)
    multe Y Tbyte [aw=_expand], control(Wbyte) mata(aw)
    multe Y Tbyte [fw=_expand], control(Wbyte) mata(fw)
    multe Y Tbyte [pw=_expand], control(Wbyte) mata(pw)

    mata assert(max(reldif(aw.estimates.est,   now.estimates.est))   > sqrt(epsilon(1)))
    mata assert(max(reldif(aw.estimates.se_po, now.estimates.se_po)) > sqrt(epsilon(1)))
    mata assert(max(reldif(aw.estimates.se_or, now.estimates.se_or)) > sqrt(epsilon(1)))

    mata assert(max(reldif(aw.estimates.est,   fw.estimates.est))   < epsilon(1))
    mata assert(max(reldif(aw.estimates.se_po, fw.estimates.se_po)) < epsilon(1))
    mata assert(max(reldif(aw.estimates.se_or, fw.estimates.se_or)) < epsilon(1))

    mata assert(max(reldif(aw.estimates.est,   pw.estimates.est))   < epsilon(1))
    mata assert(max(reldif(aw.estimates.se_po, pw.estimates.se_po)) < epsilon(1))
    mata assert(max(reldif(aw.estimates.se_or, pw.estimates.se_or)) < epsilon(1))

    expand _expand
    multe Y Tbyte, control(Wbyte) mata(ex)

    mata assert(max(reldif(aw.estimates.est,   ex.estimates.est))   < (epsilon(1)^(3/4)))
    mata assert(max(reldif(aw.estimates.se_po, ex.estimates.se_po)) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(aw.estimates.se_or, ex.estimates.se_or)) < (epsilon(1)^(3/4)))

    contract Y Tbyte Wbyte G _expand, f(_nobs)
    assert _expand == _nobs
    multe Y Tbyte [fw=_expand], control(Wbyte) mata(con)

    mata assert(max(reldif(aw.estimates.est,   con.estimates.est))   < (epsilon(1)^(3/4)))
    mata assert(max(reldif(aw.estimates.se_po, con.estimates.se_po)) < (epsilon(1)^(3/4)))
    mata assert(max(reldif(aw.estimates.se_or, con.estimates.se_or)) < (epsilon(1)^(3/4)))
end
