*! version 1.0.0 29Feb2024
*! Multiple Treatment Effects Regression
*! Based code and notes by Michal Kolesár <kolesarmi@googlemail dotcom>
*! Adapted for Stata by Mauricio Caceres Bravo <mauricio.caceres.bravo@gmail.com> and Jerray Chang <jerray@bu.edu>

capture program drop multe
program multe, eclass
    version 14.1

    * xx if replay() {
    * xx     Replay `0'
    * xx     exit 0
    * xx }

    syntax varlist(numeric fv ts) /// depvar controls (sans stratum)
           [if] [in] [aw pw],     /// subset, weights
        TREATment(varname)        /// treatment variable
    [                             ///
        STRATum(varname)          /// stratum name
        cluster(varname)          /// cluster by
        CW_uniform                ///
        GENerate(str)             /// save lambdas/taus
        MATAsave(str)             /// Save resulting mata object
        debug                     ///
    ]
    * absorb(str) /// absorb var(s)

    * NB: Sometimes Stata idiosyncrasies drive me to madness.  For example,
    * in this case to read the correct set of variables used in an interacted
    * regression into mata, I have to introduce a fake interaction with a
    * constant, lest Stata automagically exclude a base level for mata that it
    * actually included in the regression. Is Stata Basilio, keeping me locked
    * in the tower, punishing me for the mere crime of daring to code?

    tempvar cons cons2
    gen byte `cons'  = 1
    gen byte `cons2' = 1
    gettoken Y X: varlist
    local Xbak: copy local X
    local Y `Y'
    local X `cons' `X'

    if ( (`"`X'"' == "") & (`"`stratum'"' == "") ) {
        disp as err "no controls; add covariates or specify -stratum()-"
        exit 198
    }

    * Mark if, in, and any missing values by obs
    local varlist `varlist' `treatment' `stratum'
    local wbak: copy local weight
    if ( `"`weight'"' != "" ) {
        tempvar touse wgt
        qui gen double `wgt' `exp' `if' `in'
        local wcall [`weight'=`wgt']
        mark `touse' `if' `in' `wcall'
        markout `touse' `varlist', strok
    }
    else {
        local wcall
        local wgt
        marksample touse, strok
    }

    qui count if `touse'
    if ( `r(N)' == 0 ) {
        disp as txt "insufficient observations for estimation"
        exit 498
    }

    * Copy to mata for mata fun
    if `"`matasave'"' == "" local results multe_results
    else local results: copy local matasave

    * Encode variables for internal use
    qui levelsof `treatment' if `touse', loc(Tlevels)
    local Tk: list sizeof Tlevels

    tempvar T W C
    qui egen `c(obs_t)' `T' = group(`treatment') if `touse'
    if ( `"`stratum'"' != "" ) {
        qui egen `c(obs_t)' `W' = group(`stratum') if `touse'
        local Wi i.`W'

        qui levelsof `stratum' if `touse', loc(Slevels)
        local Sk: list sizeof Slevels
    }
    else {
        local Wi
        local Sk = 1
    }

    if ( `"`cluster'"' != "" ) {
        qui egen `c(obs_t)' `C' = group(`cluster') if `touse'
        local cluster: copy local C
    }

    * Get omitted varlist
    tempname zomit zomitbase
    helper_strip_omitted `X' `Wi' if `touse' `wcall', noconstant extra(b1.`T')
    local zfull  `r(expanded)'
    local zindep `r(varlist)'
    local zreg   `r(regress)'
    matrix `zomit'     = r(omit)
    matrix `zomitbase' = r(omitbase)
    foreach var of local zfull {
        if regexm("`var'", "([0-9]+).*o\.(.+)") {
            mata st_local("map", tokens(st_local("Slevels"))[`=regexs(1)'])
            local name = regexs(2)
            if ( `"`stratum'"' != "" ) {
                local name: subinstr local name "`W'" "`stratum'"
            }
            disp "note: `map'.`name' omitted because of collinearity"
        }
    }

    * 0. Base decomposition
    * ---------------------

    tempname multeworker
    mata `multeworker'  = MulTE()
    mata `results'      = MulTE_Return()
    mata `results'.full = `multeworker'.decomposition()

    * 1. Drop strata with no overlap
    * ------------------------------

    local rerun = 0
    if ( `Sk' > 1 ) {
        tempvar Wtag Ttag Tuniq Tdrop
        egen byte   `Wtag'  = tag(`W')         if `touse'
        egen byte   `Ttag'  = tag(`W' `T')     if `touse'
        egen double `Tuniq' = sum(`Ttag')      if `touse', by(`W')
        gen   byte  `Tdrop' = (`Tuniq' < `Tk') if `touse'

        qui count if (`Wtag' == 1) & (`Tdrop' == 1) & `touse'
        local nsdrop = `r(N)'

        qui count if (`Tdrop' == 1) & `touse'
        local nobs = `r(N)'

        qui replace `touse' = 0 if (`Tdrop' == 1)
        if ( `nobs' | `nsdrop' ) {
            disp as txt "dropped `nsdrop' control strata without sufficient overlap (`nobs' obs)"
            qui drop `W'
            qui egen `c(obs_t)' `W' = group(`control') if `touse'
            local rerun = 1
        }

        qui count if `touse'
        if ( `r(N)' == 0 ) {
            disp "overlap sample is empty; cannot run overlap"
        }
    }

    * 2. Drop controls that don't have within-treatment variation
    * -----------------------------------------------------------

    tempname zany zanybase zlevel zlevelbase
    mata `zany'     = J(1, `:list sizeof zfull', 0)
    mata `zanybase' = J(1, `:list sizeof zfull', 0)
    forvalues j = 1 / `Tk' {
        helper_strip_omitted `zfull' if `touse' & (`T' == `j') `wcall', noconstant
        mata `zlevel'     = st_matrix("r(omit)")
        mata `zlevelbase' = st_matrix("r(omitbase)")
        mata `zany'    [selectindex(`zlevel')]     = select(`zlevel',     `zlevel')
        mata `zanybase'[selectindex(`zlevelbase')] = select(`zlevelbase', `zlevelbase')
    }
    mata st_local("zreg2", invtokens(select(tokens(st_local("zfull")), !`zanybase')))
    mata st_matrix("`zomit'",     `zany')
    mata st_matrix("`zomitbase'", `zanybase')

    if ( "`zreg2'" != "`zreg'" ) {
        local zdrop: list zreg - zreg2
        local zreg: copy local zreg2
        disp "The following variables have no within-treatment variation and are dropped:"
        foreach var of local zdrop {
            disp "    `var'"
        }
        local rerun = 1
    }

    if ( `rerun' ) {
        mata `results'.overlap = `multeworker'.decomposition()
    }

    * Save estimates in r()
    * xx tempname estmatrix decompmatrix
    * xx mata `results'.estimates.save("`estmatrix'")
    * xx mata `results'.decomposition.save("`decompmatrix'")

    * Display standard Stata table; save in e()
    * xx Display `results', vce(`vce') touse(`touse')
    mata st_local("cmdline", "multe " + st_local("0"))
    ereturn local cmdline: copy local cmdline
    ereturn local wtype          = "`weight'"
    ereturn local wexp           = "`exp'"
    ereturn local cmd            = "multe"
    ereturn local depvar         = "`depvar'"
    ereturn local treatment      = "`treatment'"
    ereturn local controls       = "`X'"
    ereturn local stratum        = "`stratum'"
    ereturn local mata           = "`results'"

    tempname estmatrix
    mata st_matrix("`estmatrix'", `results'.full.estA#(1\0) :+ `results'.full.seP#(0\1))
    mata st_matrixcolstripe("`estmatrix'", (J(5, 1, ""), `results'.full.labels'))

    tempname sub labs
    mata `labs' = tokens(st_local("Tlevels"))[2..`Tk']
    if ( strpos("`:type `treatment''", "str") == 0 ) {
        mata `sub'  = st_varvaluelabel(st_local("treatment"))
        mata `labs' = (`sub' != "")? st_vlmap(`sub', strtoreal(`labs')): `labs'
    }
    mata st_matrixrowstripe("`estmatrix'", (J(2*(`Tk'-1), 1, ""), vec(`labs' \ J(1, `Tk'-1, "SE"))))
    matlist `estmatrix', format(%7.4g)

    ereturn matrix estmatrix = `estmatrix'
end

* xx capture program drop Replay
* xx program Replay, eclass
* xx     syntax, [vce(str) GENerate(str) DECOMPosition minmax *]
* xx     local decomp1 = ("`decomposition'" != "")
* xx     local decomp2 = (`"`generate'"'    != "")
* xx     local decomp3 = (`"`minmax'"'      != "")
* xx     local decomp  = `decomp1' | `decomp2'
* xx     if (`"`e(cmd)'"' != "multe") error 301
* xx     if ( `decomp' ) {
* xx         Decomposition, `generate'
* xx         if ( `decomp1' ) {
* xx             mata `e(mata)'.decomposition.print(`decomp3')
* xx         }
* xx         tempname decompmatrix
* xx         mata `e(mata)'.decomposition.save("`decompmatrix'")
* xx         ereturn matrix decomposition = `decompmatrix'
* xx     }
* xx     else {
* xx         Display `e(mata)', vce(`vce') repost `options'
* xx     }
* xx end
* xx 
* xx capture program drop Display
* xx program Display, eclass
* xx     syntax namelist(max = 1), [vce(str) touse(str) repost *]
* xx     * mata printf("\nTreatment Effect Estimates\n")
* xx     if "`post'" == "" local post post
* xx     FreeMatrix b V
* xx     mata `namelist'.estimates.post("`b'", "`V'", "`vce'")
* xx     mata st_local("N", strofreal(`namelist'.estimates.n))
* xx     if "`repost'" == "repost" {
* xx         ereturn repost b = `b' V = `V'
* xx     }
* xx     else {
* xx         ereturn post `b' `V', esample(`touse') obs(`N')
* xx     }
* xx 
* xx     if ( "`vce'" == "oracle" ) ereturn local vcetype "Oracle"
* xx     else ereturn local vcetype ""
* xx     ereturn local vce `vce'
* xx 
* xx     _coef_table_header, nomodeltest title(Treatment Effect Estimates)
* xx     disp ""
* xx     _coef_table, noempty `options'
* xx     //     level(95)
* xx     //     bmatrix(`b')      // e(b)
* xx     //     vmatrix(`V')      // e(V)
* xx     //     dfmatrix(matname) // e(mi_df)
* xx     //     ptitle(title)
* xx     //     coeftitle(title)
* xx     //     cititle(title)
* xx     //     cformat(format)
* xx     //     pformat(format)
* xx     //     sformat(format)
* xx end

* xx capture program drop FreeMatrix
* xx program FreeMatrix
* xx     local FreeCounter 0
* xx     local FreeMatrix
* xx     foreach FM of local 0 {
* xx         cap error 0
* xx         while ( _rc == 0 ) {
* xx             cap confirm matrix MulTE`++FreeCounter'
* xx             c_local `FM' MulTE`FreeCounter'
* xx         }
* xx     }
* xx end

capture program drop helper_strip_omitted
program helper_strip_omitted, rclass
    syntax anything(equalok) [if] [aw fw pw], [extra(str) *]
	if ( `"`weight'"' != "" ) local wcall [`weight' `exp']
    fvexpand `extra'
    local extra `r(varlist)'
    qui _rmcoll `extra' `anything' `if' `wcall', expand `options'
    local expanded `r(varlist)'

    tempname b omit omitbase final keep
    mata `keep' = (`:list sizeof extra'+1)..`:list sizeof expanded'
    matrix `b' = J(1, `:list sizeof expanded', .)
    matrix colnames `b' = `expanded'
    matrix `b' = `b'[1,(`:list sizeof extra'+1)..`:list sizeof expanded']

    _ms_omit_info `b'
    matrix `omit' = r(omit)
    mata `final' = select(st_matrixcolstripe("`b'")[., 2]', !st_matrix("r(omit)"))
    mata st_local("varlist",  invtokens(cols(`final')? `final': ""))
    mata st_local("expanded", invtokens(tokens(st_local("expanded"))[`keep']))
    local controls: list expanded - extra

    local i = 0
    matrix `omitbase' = J(1, `:list sizeof controls', 0)
    foreach var of local controls {
        local ++i
        if regexm("`var'", "([0-9]+).*o\.(.+)") {
            matrix `omitbase'[`i'] = 1
        }
    }
    mata st_local("regress", invtokens(select(tokens(st_local("controls")), !st_matrix("`omitbase'"))))

    return local expanded: copy local expanded
    return local varlist:  copy local varlist
    return local regress:  copy local regress
    return matrix omit     = `omit'
    return matrix omitbase = `omitbase'
end
