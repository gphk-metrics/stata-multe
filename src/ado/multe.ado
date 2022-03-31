*! version 0.2.1 31Mar2022
*! Multiple Treatment Effects Regression
*! Based code and notes by Michal Kolesár <kolesarmi@googlemail dotcom>
*! Adapted for Stata by Mauricio Caceres Bravo <mauricio.caceres.bravo@gmail.com>

capture program drop multe
program multe, eclass
    version 14.1

    if replay() {
        Replay `0'
        exit 0
    }

    local 0bak: copy local 0
    syntax varlist(numeric fv ts min=2 max=2)  /// depvar treatment
           [if] [in] ,                         /// subset
        control(varlist)                       /// control variable
    [                                          ///
        vce(str)                               /// SEs to be displayed
        GENerate(str)                          /// save lambdas/taus
        MATAsave(str)                          /// Save resulting mata object
    ]

    local vce `vce'
    if !inlist("`vce'", "", "oracle") {
        disp as err "vce() option `vce' not known"
        exit 198
    }

    * Varlist must be depvar and indepvar
    gettoken depvar treatment: varlist
    local depvar    `depvar'
    local treatment `treatment'

    * Mark if, in, and any missing values by obs
    local varlist `varlist' `control'
    marksample touse, strok

    * Copy to mata for mata fun
    if "`matasave'" == "" local results MulTEResults
    else local results: copy local matasave

    * Check each control strata has all treatment levels
* TODO: xx this is very inefficient; optimize
    qui levelsof `treatment' if `touse', loc(Tlevels)
    local Tk: list sizeof Tlevels

    tempvar W Wtag Ttag Tuniq Tdrop
    egen long  `W'      = group(`control')
    egen byte  `Wtag'   = tag(`W') if `touse'
    egen byte  `Ttag'   = tag(`W' `treatment') if `touse'
    egen double `Tuniq' = sum(`Ttag') if `touse', by(`W')
    gen byte `Tdrop'    = (`Tuniq' < `Tk') if `touse'
    cap assert `Tuniq' <= `Tk' if `touse'
    if ( _rc ) {
        disp as err "found more levels within a control strata than levels overall"
        exit 9
    }

    qui count if (`Wtag' == 1) & (`Tdrop' == 1) & `touse'
    local nstrata = `r(N)'
    qui count if (`Tdrop' == 1) & `touse'
    local nobs = `r(N)'
    qui replace `touse' = 0 if (`Tdrop' == 1)
    if ( `nobs' | `nstrata' ) {
        disp as txt "`nstrata' control strata without sufficient overlap dropped (`nobs' obs)"
        * exit 498
    }

    mata Wm = designmatrix(st_data(., "`W'", "`touse'"))
    mata `results' = MulTE("`depvar'", "`treatment'", Wm, "`touse'")

    * Save lambda and tauhat, if requested
    tempname types names
    local 0, `generate'
    syntax, [lambda LAMBDAprefix(str) tau TAUprefix(str)]

    local savelambda = ("`lambdaprefix'" != "") | ("`lambda'" != "")
    local savetau    = ("`tauprefix'"    != "") | ("`tau'"    != "")

    if "`lambdaprefix'" == "" local lambdaprefix lambda
    if "`tauprefix'"    == "" local tauprefix    tau

    if ( `savelambda' ) {
        mata: `names' = `results'.decomposition.lambda_names
        mata: `names' = "`lambdaprefix'" :+ `names'
        mata: `types' = J(1, cols(`names'), "`:set type'")
        mata: (void) st_addvar(`types', `names')
        mata: (void) st_store(., `names', `results'.decomposition.lambda(Wm))
    }

    if ( `savetau' ) {
        mata: `names' = `results'.decomposition.tauhat_names
        mata: `names' = "`tauprefix'" :+ `names'
        mata: `types' = J(1, cols(`names'), "`:set type'")
        mata: (void) st_addvar(`types', `names')
        mata: (void) st_store(., `names', `results'.decomposition.tauhat(Wm))
    }

    mata mata drop Wm

    * Save estimates in r()
    tempname estmatrix decompmatrix
    mata `results'.estimates.save("`estmatrix'")
    mata `results'.decomposition.save("`decompmatrix'")

    * Display standard Stata table; save in e()
    Display `results', vce(`vce') touse(`touse')
    mata st_local("cmdline", "multe " + st_local("0bak"))
    ereturn local cmd     = "multe"
    ereturn local cmdline: copy local cmdline
    ereturn local depvar  = "`depvar'"
    ereturn local control = "`control'"
    ereturn local predict = "multe_p"
    ereturn local mata    = "`results'"

    ereturn matrix estimates     = `estmatrix'
    ereturn matrix decomposition = `decompmatrix'
end

capture program drop Replay
program Replay, eclass
    syntax, [vce(str) *]
    if (`"`e(cmd)'"' != "multe") error 301
    Display `e(mata)', vce(`vce') repost `options'
end

capture program drop Display
program Display, eclass
    syntax namelist(max = 1), [vce(str) touse(str) repost *]
    if "`post'" == "" local post post
    FreeMatrix b V
    mata `namelist'.estimates.post("`b'", "`V'", "`vce'")
    mata st_local("N", strofreal(`namelist'.estimates.n))
    if "`repost'" == "repost" {
        ereturn repost b = `b' V = `V'
    }
    else {
        ereturn post `b' `V', esample(`touse') obs(`N')
    }

    if ( "`vce'" == "oracle" ) ereturn local vcetype "Oracle"
    else ereturn local vcetype ""
    ereturn local vce `vce'

    _coef_table, noempty noomitted `options'
    disp as txt "({bf:warning:} off-diagonal entries of Vcov matrix hard-coded as 0)"
    //     level(95)
    //     bmatrix(`b')      // e(b)
    //     vmatrix(`V')      // e(V)
    //     dfmatrix(matname) // e(mi_df)
    //     ptitle(title)
    //     coeftitle(title)
    //     cititle(title)
    //     cformat(format)
    //     pformat(format)
    //     sformat(format)
end

capture program drop FreeMatrix
program FreeMatrix
    local FreeCounter 0
    local FreeMatrix
    foreach FM of local 0 {
        cap error 0
        while ( _rc == 0 ) {
            cap confirm matrix MulTE`++FreeCounter'
            c_local `FM' MulTE`FreeCounter'
        }
    }
end
