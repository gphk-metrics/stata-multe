*! version 0.3.1 31May2022
*! Multiple Treatment Effects Regression
*! Based code and notes by Michal Kolesár <kolesarmi@googlemail dotcom>
*! Adapted for Stata by Mauricio Caceres Bravo <mauricio.caceres.bravo@gmail.com> and Jerray Chang <jerray@bu.edu>

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
        DECOMPosition                          /// Compute and display decomposition
        minmax                                 /// Print full decomposition table
    ]

    local decomp1 = ("`decomposition'" != "")
    local decomp2 = (`"`generate'"'    != "")
    local decomp3 = (`"`minmax'"'      != "")
    local decomp  = `decomp1' | `decomp2'

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
    if "`matasave'" == "" local results multe_results
    else local results: copy local matasave

    * Check each control strata has all treatment levels
    * NB: This seems inefficient but doesn't in practice take very long
    qui levelsof `treatment' if `touse', loc(Tlevels)
    local Tk: list sizeof Tlevels

    tempvar W Wtag Ttag Tuniq Tdrop
    egen `c(obs_t)' `W'     = group(`control')
    egen byte       `Wtag'  = tag(`W') if `touse'
    egen byte       `Ttag'  = tag(`W' `treatment') if `touse'
    egen double     `Tuniq' = sum(`Ttag') if `touse', by(`W')
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

    qui drop `W'
    egen `c(obs_t)' `W' = group(`control') if `touse'
    mata `results' = MulTE()
    mata `results'.estimates("`depvar'", "`treatment'", "`W'", "`touse'")

    if ( `decomp' ) {
        mata `results'.decomposition("`depvar'", "`treatment'", "`W'", "`touse'")

        * Save lambda and tauhat, if requested
        LambdaTau, results(`results') `generate'
    }
    mata `results'.cache_drop()

    * Save estimates in r()
    tempname estmatrix decompmatrix
    mata `results'.estimates.save("`estmatrix'")
    mata `results'.decomposition.save("`decompmatrix'")

    * Display standard Stata table; save in e()
    Display `results', vce(`vce') touse(`touse')
    mata st_local("cmdline", "multe " + st_local("0bak"))
    ereturn local cmdline: copy local cmdline
    ereturn local cmd       = "multe"
    ereturn local depvar    = "`depvar'"
    ereturn local treatment = "`treatment'"
    ereturn local control   = "`control'"
    ereturn local mata      = "`results'"

    ereturn matrix estimates = `estmatrix'
    if ( `decomp' ) {
        ereturn matrix decomposition = `decompmatrix'
        if ( `decomp1' ) {
            mata `results'.decomposition.print(`decomp3')
        }
    }
end

capture program drop LambdaTau
program LambdaTau
    syntax, results(str) [lambda LAMBDAprefix(str) tau TAUprefix(str)]
    tempname types names

    local savelambda = ("`lambdaprefix'" != "") | ("`lambda'" != "")
    local savetau    = ("`tauprefix'"    != "") | ("`tau'"    != "")

    if "`lambdaprefix'" == "" local lambdaprefix lambda
    if "`tauprefix'"    == "" local tauprefix    tau

    if ( `savelambda' ) {
        mata: `names' = `results'.decomposition.lambda_names
        mata: `names' = "`lambdaprefix'" :+ `names'
        mata: `types' = J(1, cols(`names'), "`:set type'")
        mata: (void) st_addvar(`types', `names')
        mata: (void) st_store(., `names', `results'.decomposition.lambda(`results'.Wm))
    }

    if ( `savetau' ) {
        mata: `names' = `results'.decomposition.tauhat_names
        mata: `names' = "`tauprefix'" :+ `names'
        mata: `types' = J(1, cols(`names'), "`:set type'")
        mata: (void) st_addvar(`types', `names')
        mata: (void) st_store(., `names', `results'.decomposition.tauhat(`results'.Wm))
    }
end

capture program drop Replay
program Replay, eclass
    syntax, [vce(str) GENerate(str) DECOMPosition minmax *]
    local decomp1 = ("`decomposition'" != "")
    local decomp2 = (`"`generate'"'    != "")
    local decomp3 = (`"`minmax'"'      != "")
    local decomp  = `decomp1' | `decomp2'
    if (`"`e(cmd)'"' != "multe") error 301
    if ( `decomp' ) {
        Decomposition, `generate'
        if ( `decomp1' ) {
            mata `e(mata)'.decomposition.print(`decomp3')
        }
        tempname decompmatrix
        mata `e(mata)'.decomposition.save("`decompmatrix'")
        ereturn matrix decomposition = `decompmatrix'
    }
    else {
        Display `e(mata)', vce(`vce') repost `options'
    }
end

capture program drop Decomposition
program Decomposition
    syntax, [*]
    tempvar touse W
    gen byte `touse' = e(sample)
    egen `c(obs_t)' `W' = group(`e(control)') if `touse'
    mata `e(mata)'.cache_load("`e(depvar)'", "`e(treatment)'", "`W'", "`touse'")
    mata `e(mata)'.decomposition("`e(depvar)'", "`e(treatment)'", "`W'", "`touse'")
    LambdaTau, results(`e(mata)') `options'
    mata `e(mata)'.cache_drop()
end

capture program drop Display
program Display, eclass
    syntax namelist(max = 1), [vce(str) touse(str) repost *]
    mata printf("\nTreatment Effect Estimates\n")
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

    _coef_table, noempty `options'
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
