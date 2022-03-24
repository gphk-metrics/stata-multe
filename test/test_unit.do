capture program drop multe_unit_tests
program  multe_unit_tests
    syntax, [Verbose]
    multe_load_test_data
    multe_unit_test_Ttype, `verbose'
    multe_unit_test_Wtype, `verbose'
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

    gen Y = Tbyte - Wbyte + runiform()
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
    local Tpass Tbyte Tint Tlong Tfloat Tdouble
    local Tfail Tstr
    foreach T in `Tpass' `Tfail' {
        cap multe Y `T', control(Wbyte)
        local rc1 = _rc
        if ( `:list T in Tpass' ) {
            mata st_numscalar("`Texpected'", sort(uniqrows(st_data(., "`T'")), 1)[1])
            mata st_numscalar("`Tresult'",   `e(mata)'.estimates.Tvalues[1])
            cap assert scalar(`Texpected') == scalar(`Tresult')
            local rc2 = _rc

            if ( `rc1' != 0 ) {
                disp "(multe test fail): treatment type `:subinstr local T "T" ""' failed with _rc = `rc1'"
                exit 9
            }
            else if ( `rc2' != 0 ) {
                disp "(multe test fail): treatment type `:subinstr local T "T" ""' levels not internally sorted"
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
        cap multe Y Tbyte, control(`W')
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
