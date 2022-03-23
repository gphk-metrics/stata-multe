*! version 0.1.5 23Mar2022
*! Multiple Treatment Effects Regression
*! Based code and notes by Michal Kolesár <kolesarmi@googlemail dotcom>
*! Adapted for Stata by Mauricio Caceres Bravo <mauricio.caceres.bravo@gmail.com>

capture program drop multe
program multe, rclass
    version 14.1

    syntax varlist(numeric fv ts min=2 max=2)  /// depvar treatment
           [if] [in] ,                         /// subset
        control(varname)                       /// control variable
    [                                          ///
        GENerate(str)                          /// save lambdas/taus
        MATAsave(str)                          /// Save resulting mata object
    ]

    * Varlist must be depvar and indepvar
    gettoken depvar treatment: varlist
    local depvar    `depvar'
    local treatment `treatment'

    * Mark if, in, and any missing values by obs
    local varlist `varlist' `control'
    marksample touse, strok

    * Copy to mata for mata fun
    if "`matasave'" == "" tempname results
    else local results: copy local matasave

    * Force treatment, control into indices
    tempvar T W
    egen `T' = group(`treatment'), label
    egen `W' = group(`control')

    mata Wm = designmatrix(st_data(., "`W'", "`touse'"))
    mata `results' = MulTE("`depvar'", "`T'", Wm, "`touse'")
    mata mata drop Wm

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
        mata: (void) st_store(., `names', `results'.decomposition.lambda)
    }
    if ( `savetau' ) {
        mata: `names' = `results'.decomposition.tauhat_names
        mata: `names' = "`tauprefix'" :+ `names'
        mata: `types' = J(1, cols(`names'), "`:set type'")
        mata: (void) st_addvar(`types', `names')
        mata: (void) st_store(., `names', `results'.decomposition.tauhat)
    }

    * Save estimates in r()
    tempname estmatrix decompmatrix
    mata `results'.estimates.print()
    mata `results'.estimates.save("`estmatrix'")
    mata `results'.decomposition.save("`decompmatrix'")
    return matrix estimates     = `estmatrix'
    return matrix decomposition = `decompmatrix'

    * Cleanup mata object if no save
    if "`matasave'" == "" mata mata drop `results'
end
