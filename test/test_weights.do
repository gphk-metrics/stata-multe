capture program drop multe_weight_tests
program multe_weight_tests
* TODO: xx test lambda and tau (i.e. decomposition)
    multe_load_test_data
    gen `c(obs_t)' G = _n
    gen `c(obs_t)' _expand = ceil(runiform() * _N/50)
    gen double Y = Ydouble + rnormal() * _N / 50
    qui sum _expand
    local ratio = sqrt(r(sum) / _N)

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

    mata assert(max(reldif(aw.estimates.est, fw.estimates.est)) < (epsilon(1)^(3/4)))
    mata assert(max(abs(`ratio' :- aw.estimates.se_po :/ fw.estimates.se_po)) < (epsilon(1)^(3/4)))
    mata assert(max(abs(`ratio' :- aw.estimates.se_or :/ fw.estimates.se_or)) < (epsilon(1)^(3/4)))

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

    mata assert(max(reldif(aw.estimates.est, ex.estimates.est)) < (epsilon(1)^(3/4)))
    mata assert(max(abs(`ratio' :- aw.estimates.se_po :/ ex.estimates.se_po)) < (epsilon(1)^(3/4)))
    mata assert(max(abs(`ratio' :- aw.estimates.se_or :/ ex.estimates.se_or)) < (epsilon(1)^(3/4)))

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
