*! version 1.0.1 08Mar2024
*! Multiple Treatment Effects Regression
*! Based code and notes by Michal Kolesár <kolesarmi@googlemail dotcom>
*! Adapted for Stata by Mauricio Caceres Bravo <mauricio.caceres.bravo@gmail.com> and Jerray Chang <jerray@bu.edu>

capture program drop multe
program multe, eclass
    version 14.1

    qui syntax [anything(everything)] [if] [in] [aw pw], [ESTimates(str) full overlap diff oracle *]
    if `"`anything'"' == "" {
        if ( `"`e(cmd)'"' != "multe" ) error 301
        if ( "`estimates'`full'`overlap'`diff'`oracle'" == "" & replay() ) {
            Display `e(mata)', `e(displayopts)' repost
        }
        else {
            Display `e(mata)', est(`estimates') `full' `overlap' `diff' `oracle' `options' repost
        }
        exit 0
    }

    syntax varlist(numeric fv) /// depvar controls (sans stratum)
           [if] [in] [aw pw],  /// subset, weights
        TREATment(varname)     /// treatment variable
    [                          ///
        STRATum(varname)       /// stratum name
        cluster(varname)       /// cluster by
        CW_uniform             ///
        GENerate(str)          /// save lambdas/taus
        MATAsave(str)          /// Save resulting mata object
        debug                  ///
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
    mata `multeworker'        = MulTE()
    mata `results'            = MulTE_Return()
    mata `results'.Tk         = `Tk'
    mata `results'.Tvar       = st_local("treatment")
    mata `results'.Yvar       = st_local("Y")
    mata `results'.Tlevels    = tokens(st_local("Tlevels"))
    mata `results'.full       = `multeworker'.decomposition()
    mata `results'.full.touse = st_data(., st_local("touse"))
    mata `results'.full.N     = sum(`results'.full.touse)

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
        disp "The following variables have no within-treatment variation"
        disp "and are dropped from the overlap sample:"
        foreach var of local zdrop {
            disp "    `var'"
        }
        local rerun = 1
    }

    * Re-run in overlap sample
    if ( `rerun' ) {
        mata `results'.overlap = `multeworker'.decomposition()
        mata `results'.overlap.touse = st_data(., st_local("touse"))
        mata `results'.overlap.N     = sum(`results'.overlap.touse)
        mata `results'.has_overlap   = 1
    }

    * Display standard Stata table; save in e()
    * -----------------------------------------

    Display `results', est(PL) cluster(`cluster')
    mata st_local("cmdline", "multe " + st_local("0"))
    ereturn local cmdline: copy local cmdline
    ereturn local wtype          = "`weight'"
    ereturn local wexp           = "`exp'"
    ereturn local cmd            = "multe"
    ereturn local depvar         = "`Y'"
    ereturn local treatment      = "`treatment'"
    ereturn local controls       = "`X'"
    ereturn local cluster        = "`cluster'"
    ereturn local stratum        = "`stratum'"
    ereturn local mata           = "`results'"

    * Save full results
    tempname full_beta full_pop_se full_oracle_se full_diff full_diff_se levels

    mata st_matrix(st_local("full_beta"),      `results'.full.estA)
    mata st_matrix(st_local("full_diff"),      `results'.full.estB)
    mata st_matrix(st_local("full_pop_se"),    `results'.full.seP)
    mata st_matrix(st_local("full_diff_se"),   `results'.full.seB)
    mata st_matrix(st_local("full_oracle_se"), `results'.full.seO)

    mata st_matrixcolstripe(st_local("full_beta"),      (J(5, 1, ""), `results'.full.labels'))
    mata st_matrixcolstripe(st_local("full_diff"),      (J(4, 1, ""), `results'.full.labels[2..5]'))
    mata st_matrixcolstripe(st_local("full_pop_se"),    (J(5, 1, ""), `results'.full.labels'))
    mata st_matrixcolstripe(st_local("full_diff_se"),   (J(4, 1, ""), `results'.full.labels[2..5]'))
    mata st_matrixcolstripe(st_local("full_oracle_se"), (J(5, 1, ""), `results'.full.labels'))

    mata `levels' = `results'.Tlevels[2..`Tk']'
    mata st_matrixrowstripe(st_local("full_beta"),      (J(`Tk'-1, 1, ""), `levels'))
    mata st_matrixrowstripe(st_local("full_diff"),      (J(`Tk'-1, 1, ""), `levels'))
    mata st_matrixrowstripe(st_local("full_pop_se"),    (J(`Tk'-1, 1, ""), `levels'))
    mata st_matrixrowstripe(st_local("full_diff_se"),   (J(`Tk'-1, 1, ""), `levels'))
    mata st_matrixrowstripe(st_local("full_oracle_se"), (J(`Tk'-1, 1, ""), `levels'))

    ereturn matrix full_beta      = `full_beta'
    ereturn matrix full_pop_se    = `full_pop_se'
    ereturn matrix full_oracle_se = `full_oracle_se'
    ereturn matrix full_diff      = `full_diff'
    ereturn matrix full_diff_se   = `full_diff_se'

    * Results from re-run in overlap sample
    * -------------------------------------

    if ( `rerun' ) {
        tempname overlap_beta overlap_pop_se overlap_oracle_se overlap_diff overlap_diff_se

        mata st_matrix(st_local("overlap_beta"),      `results'.overlap.estA)
        mata st_matrix(st_local("overlap_diff"),      `results'.overlap.estB)
        mata st_matrix(st_local("overlap_pop_se"),    `results'.overlap.seP)
        mata st_matrix(st_local("overlap_diff_se"),   `results'.overlap.seB)
        mata st_matrix(st_local("overlap_oracle_se"), `results'.overlap.seO)

        mata st_matrixcolstripe(st_local("overlap_beta"),      (J(5, 1, ""), `results'.overlap.labels'))
        mata st_matrixcolstripe(st_local("overlap_diff"),      (J(4, 1, ""), `results'.overlap.labels[2..5]'))
        mata st_matrixcolstripe(st_local("overlap_pop_se"),    (J(5, 1, ""), `results'.overlap.labels'))
        mata st_matrixcolstripe(st_local("overlap_diff_se"),   (J(4, 1, ""), `results'.overlap.labels[2..5]'))
        mata st_matrixcolstripe(st_local("overlap_oracle_se"), (J(5, 1, ""), `results'.overlap.labels'))

        mata st_matrixrowstripe(st_local("overlap_beta"),      (J(`Tk'-1, 1, ""), `levels'))
        mata st_matrixrowstripe(st_local("overlap_diff"),      (J(`Tk'-1, 1, ""), `levels'))
        mata st_matrixrowstripe(st_local("overlap_pop_se"),    (J(`Tk'-1, 1, ""), `levels'))
        mata st_matrixrowstripe(st_local("overlap_diff_se"),   (J(`Tk'-1, 1, ""), `levels'))
        mata st_matrixrowstripe(st_local("overlap_oracle_se"), (J(`Tk'-1, 1, ""), `levels'))

        ereturn matrix overlap_beta      = `overlap_beta'
        ereturn matrix overlap_pop_se    = `overlap_pop_se'
        ereturn matrix overlap_oracle_se = `overlap_oracle_se'
        ereturn matrix overlap_diff      = `overlap_diff'
        ereturn matrix overlap_diff_se   = `overlap_diff_se'
    }
    else {
        mata `results'.has_overlap = 0
    }

    * Show estimates (including tests)
    tempname fullmatrix diffmatrix
    mata st_matrix("`fullmatrix'", `results'.full.estA#(1\0) :+ `results'.full.seP#(0\1))
    mata st_matrix("`diffmatrix'", `results'.full.estB#(1\0) :+ `results'.full.seB#(0\1))
    mata st_matrixcolstripe("`fullmatrix'", (J(5, 1, ""), `results'.full.labels'))
    mata st_matrixcolstripe("`diffmatrix'", (J(4, 1, ""), `results'.full.labels[2..5]'))

    tempname sub labs
    mata `labs' = tokens(st_local("Tlevels"))[2..`Tk']
    if ( strpos("`:type `treatment''", "str") == 0 ) {
        mata `sub'  = st_varvaluelabel(st_local("treatment"))
        mata `labs' = (`sub' != "")? st_vlmap(`sub', strtoreal(`labs')): `labs'
    }
    mata st_matrixrowstripe("`fullmatrix'", (J(2*(`Tk'-1), 1, ""), vec(`labs' \ J(1, `Tk'-1, "SE"))))
    mata st_matrixrowstripe("`diffmatrix'", (J(2*(`Tk'-1), 1, ""), vec(`labs' \ J(1, `Tk'-1, "SE"))))

    mata printf("\nAlternative Estimates on Full Sample:\n")
    matlist `fullmatrix', format(%7.4g)
    ereturn matrix fullmatrix = `fullmatrix'
    ereturn matrix diffmatrix = `diffmatrix'

    mata printf("\nP-values for null hypothesis of no propensity score variation:\n")
    mata printf("Wald test: %9.6g\n", `results'.full.Wa.pval)
    mata printf("  LM test: %9.6g\n", `results'.full.LM.pval)

    if ( `rerun' ) {
        tempname overlapmatrix overlapdiffmatrix
        mata st_matrix("`overlapmatrix'",     `results'.overlap.estA#(1\0) :+ `results'.overlap.seP#(0\1))
        mata st_matrix("`overlapdiffmatrix'", `results'.overlap.estB#(1\0) :+ `results'.overlap.seB#(0\1))

        mata st_matrixcolstripe("`overlapmatrix'",     (J(5, 1, ""), `results'.overlap.labels'))
        mata st_matrixcolstripe("`overlapdiffmatrix'", (J(4, 1, ""), `results'.overlap.labels[2..5]'))
        mata st_matrixrowstripe("`overlapmatrix'",     (J(2*(`Tk'-1), 1, ""), vec(`labs' \ J(1, `Tk'-1, "SE"))))
        mata st_matrixrowstripe("`overlapdiffmatrix'", (J(2*(`Tk'-1), 1, ""), vec(`labs' \ J(1, `Tk'-1, "SE"))))

        mata printf("\nAlternativ Estimates on Overlap Sample:")
        matlist `overlapmatrix', format(%7.4g)
        ereturn matrix overlapmatrix     = `overlapmatrix'
        ereturn matrix overlapdiffmatrix = `overlapdiffmatrix'
    }

    disp ""
    disp "Note: You can post any combination of results from the table to Stata:"
    disp ""
    disp "    multe, est(estimate) [{full|overlap} diff oracle]"
    disp ""
    disp "Examples:"
    disp ""
    disp "    {stata multe, est(ATE) full oracle}"
    disp "    {stata multe, est(CW)  overlap diff}"
end

capture program drop Display
program Display, eclass
    syntax namelist(max=1), ESTimates(str) [full overlap diff oracle repost cluster(str) *]

    if ( ("`full'" != "") & ("`overlap'" != "") ) {
        disp as err "unable to display both full and overlap sample results"
        exit 198
    }

    if ( ("`full'" == "") & ("`overlap'" == "") ) {
        local full full
    }

    local estimates = upper("`estimates'")
    if ( !inlist(upper("`estimates'"), "PL", "OWN", "ATE", "EW", "CW") ) {
        disp as err "estimates(`estimates') unknown; must be one of PL, OWN, ATE, EW, CW"
        exit 198
    }

    if ( ("`diff'" != "") & ("`estimates'" == "PL") ) {
        disp as err "diff estimates only available relative to PL"
        exit 198
    }

    if ( ("`diff'" != "") & ("`oracle'" == "oracle") ) {
        disp as err "oracle SE not available for diff estimates"
        exit 198
    }

    if ( ("`oracle'" == "oracle") & inlist("`estimates'", "PL", "OWN") ) {
        disp as err "oracle SE not available for PL, OWN"
        exit 198
    }

    mata st_local("has_overlap", strofreal(`namelist'.has_overlap))
    if ( ("`overlap'" != "") & (`has_overlap' == 0) ) {
        disp as err "requested overlap sample results but overlap was not computed"
        exit 198
    }

    local labPL  "Partly Linear Model"
    local labOWN "Own Treatment Effects"
    local labATE "ATE"
    local labEW  "Easiest-to-estimate Weighted ATE"
    local labCW  "Easiest-to-estimate Common Weighted ATE"

    FreeMatrix b V
    tempname index
    mata `index' = selectindex(`namelist'.`full'`overlap'.labels :== "`estimates'")
    if ( "`diff'" == "" ) {
        mata st_matrix(st_local("b"), `namelist'.`full'`overlap'.estA[., `index'])
        if ( "`oracle'" == "oracle" ) {
            mata st_matrix(st_local("V"), `namelist'.`full'`overlap'.Vo_`estimates')
        }
        else {
            mata st_matrix(st_local("V"), `namelist'.`full'`overlap'.Vpop_`estimates')
        }
    }
    else {
        mata st_matrix(st_local("b"), `namelist'.`full'`overlap'.estB[., `index'-1])
        mata st_matrix(st_local("V"), `namelist'.`full'`overlap'.Vdiff_`estimates')
    }
    mata multe_helper_post(st_local("b"), st_local("V"), `namelist')

    tempvar touse
    qui gen byte `touse' = .
    mata st_store(., st_local("touse"), `namelist'.`full'`overlap'.touse)
    qui count if `touse'

    if "`repost'" == "repost" {
        ereturn scalar N = r(N)
        ereturn repost b=`b' V=`V', esample(`touse')
    }
    else ereturn post `b' `V', esample(`touse') obs(`r(N)')

    if ( "`cluster'" == "" ) {
        if ( "`oracle'" == "oracle" ) ereturn local vcetype "Oracle"
        else ereturn local vcetype ""
    }
    else {
        if ( "`oracle'" == "oracle" ) ereturn local vcetype "Oracle Cluster"
        else ereturn local vcetype "Cluster"
    }
    ereturn local displayopts estimates(`estimates') `full' `overlap' `diff' `oracle'

    _coef_table_header, nomodeltest title(`lab`estimates'' Estimates (`full'`overlap' sample))
    disp ""
    _coef_table, noempty `options'
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
        if ( regexm("`var'", "([0-9]+|^).*o\.(.+)") ) {
            matrix `omitbase'[1, `i'] = 1
        }
    }
    mata st_local("regress", invtokens(select(tokens(st_local("controls")), !st_matrix("`omitbase'"))))

    return local expanded: copy local expanded
    return local varlist:  copy local varlist
    return local regress:  copy local regress
    return matrix omit     = `omit'
    return matrix omitbase = `omitbase'
end
