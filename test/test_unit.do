capture program drop multe_unit_tests
program  multe_unit_tests
    syntax, [Verbose]
    multe_load_test_data
    multe_unit_fail
    multe_unit_test_Ytype, `verbose'
    multe_unit_test_Ttype, `verbose'
    multe_unit_test_Wtype, `verbose'
end

capture program drop multe_unit_fail
program multe_unit_fail
    * missing control|strat
    cap multe Ydouble, treat(Tbyte)
    local rc = ( _rc == 198 )
    * missing treatment
    cap multe Ydouble C*
    local rc = `rc' & ( _rc == 198 )
    if ( `rc' == 1 ) {
        disp "(multe test success): fail checks present"
    }
    else {
        disp "(multe test fail): fail checks missing; check stop logic"
        exit 9
    }
end

capture program drop multe_load_test_data
program multe_load_test_data
    syntax, [nobs(int 1000) ktreat(int 5)]
    clear
    qui set obs `nobs'
    gen byte   Tbyte   = ceil(runiform() * `ktreat')
    gen int    Tint    = Tbyte * 1234
    gen long   Tlong   = Tbyte * 123456
    gen float  Tfloat  = (Tbyte - 2.5) * `c(pi)'
    gen double Tdouble = (Tbyte - 2.5) * `c(pi)' / 123456
    gen str4   Tstr    = "str" + string(Tbyte)
    gen byte   Wbyte   = mod(_n, 10)
    gen int    Wint    = Wbyte * 1234
    gen long   Wlong   = Wbyte * 123456
    gen float  Wfloat  = (Wbyte - 2.5) * `c(pi)'
    gen double Wdouble = (Wbyte - 2.5) * `c(pi)' / 123456
    gen str4   Wstr    = "str" + string(Wbyte)
    gen byte   Cbyte   = mod(_n, 10)
    gen int    Cint    = Wbyte * 1234
    gen long   Clong   = Wbyte * 123456
    gen float  Cfloat  = (Wbyte - 2.5) * `c(pi)'
    gen double Cdouble = (Wbyte - 2.5) * `c(pi)' / 123456
    gen byte   Ybyte   = Tbyte - Wbyte + Cdouble + runiform()
    gen int    Yint    = Tbyte - Wbyte + Cdouble + runiform()
    gen long   Ylong   = Tbyte - Wbyte + Cdouble + runiform()
    gen float  Yfloat  = Tbyte - Wbyte + Cdouble + runiform()
    gen double Ydouble = Tbyte - Wbyte + Cdouble + runiform()
    gen str4   Ystr    = "str" + string(Tbyte - Wbyte + runiform())
end

capture program drop multe_unit_test_Ytype
program multe_unit_test_Ytype
    syntax, [Reload Verbose *]

    if "`reload'" != "" multe_load_test_data, `options'
    if ( "`verbose'" != "" ) {
        foreach Y in Ybyte Yint Ylong Yfloat Ydouble Ystr {
            tab `Y'
        }
    }

    tempname Yexpected Yresult
    local Ypass Ybyte Yint Ylong Yfloat Ydouble
    local Yfail Ystr
    foreach Y in `Ypass' `Yfail' {
        cap multe `Y', treat(Tbyte) strat(Wbyte)
        local rc1 = _rc
        cap multe `Y' C*, treat(Tbyte) strat(Wbyte)
        local rc2 = _rc
        cap multe `Y' Cdouble, treat(Tbyte)
        local rc3 = _rc
        if ( `:list Y in Ypass' ) {
            if ( (`rc1' != 0)  | (`rc2' != 0) | (`rc3' != 0) ) {
                disp "(multe test fail): depvar type `:subinstr local Y "Y" ""' failed with _rc = `rc1'"
                exit 9
            }
            else {
                disp "(multe test success): depvar type `:subinstr local Y "Y" ""' gave no errors"
            }
        }
        if ( `:list Y in Yfail' ) {
            if ( (`rc1' != 109) | (`rc2' != 109) | (`rc3' != 109) ) {
                disp "(multe test fail): depvar type `:subinstr local Y "Y" ""' did not fail with _rc =  109 as expected (_rc = `rc1')"
                exit 9
            }
            else {
                disp "(multe test success): depvar type `:subinstr local Y "Y" ""' failed with _rc = 109 as expected"
            }
        }
    }
end

* Double-check: Allow multe to take any multi-valued treatment T and
* automatically take the lowest value as the control
capture program drop multe_unit_test_Ttype
program multe_unit_test_Ttype
    syntax, [Reload Verbose *]

    if "`reload'" != "" multe_load_test_data, `options'
    if ( "`verbose'" != "" ) {
        foreach T in Tbyte Tint Tlong Tfloat Tdouble Tstr {
            tab `T'
        }
    }

    tempname Texpected Tresult
    local Tpass Tbyte Tint Tlong Tfloat Tdouble Tstr
    foreach T in `Tpass' `Tfail' {
        cap multe Ydouble Cdouble, treat(`T') strat(Wbyte)
        local rc1 = _rc
        if ( `:list T in Tpass' ) {
            if strpos("`:type `T''", "str") {
                cap mata assert(all(sort(uniqrows(st_sdata(., "`T'")), 1) :== `e(mata)'.Tlevels'))
            }
            else {
                cap mata assert(all(reldif(sort(uniqrows(st_data(., "`T'")), 1), strtoreal(`e(mata)'.Tlevels')) :< epsilon(1)^(0.75)))
            }
            local rc2 = _rc
            if ( `rc1' != 0 ) {
                disp "(multe test fail): treatment type `:subinstr local T "T" ""' failed with _rc = `rc1'"
                exit 9
            }
            else if ( `rc2' != 0 ) {
                disp "(multe test fail): treatment type `:subinstr local T "T" ""' levels computed internally do not match the data"
                exit 9
            }
            else {
                disp "(multe test success): treatment type `:subinstr local T "T" ""' gave no errors"
            }
        }
        if ( `:list T in Tfail' ) {
            if ( `rc1' != 109 ) {
                disp "(multe test fail): treatment type `:subinstr local T "T" ""' did not fail with _rc =  109 as expected (_rc = `rc1')"
                exit 9
            }
            else {
                disp "(multe test success): treatment type `:subinstr local T "T" ""' failed with _rc = 109 as expected"
            }
        }
    }
end

* Double-check: Allow multe to take any multi-valued control W
capture program drop multe_unit_test_Wtype
program multe_unit_test_Wtype
    syntax, [Verbose *]

    if "`reload'" != "" multe_load_test_data, `options'
    if ( "`verbose'" != "" ) {
        foreach W in Wbyte Wint Wlong Wfloat Wdouble Wstr {
            tab `W'
        }
    }

    tempname Wexpected Wresult
    local Wpass Wbyte Wint Wlong Wfloat Wdouble Wstr
    foreach W in `Wpass' {
        cap multe Ydouble Cdouble, treat(Tbyte) strat(`W')
        local rc1 = _rc
        if ( `rc1' != 0 ) {
            disp "(multe test fail): control type `:subinstr local W "W" ""' failed with _rc = `rc1'"
            exit 9
        }
        else {
            disp "(multe test success): control type `:subinstr local W "W" ""' gave no errors"
        }
    }
end
