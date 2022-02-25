*! version 0.1.2 25Feb2022
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

    mata Wm = designmatrix(st_data(., "`control'", "`touse'"))
    mata `results' = MulTE("`depvar'", "`treatment'", Wm, "`touse'")

    tempname estmatrix decompmatrix
    mata `results'.estimates.print()
    mata `results'.estimates.save("`estmatrix'")
    mata `results'.decomposition.save("`decompmatrix'")
    return matrix estimates     = `estmatrix'
    return matrix decomposition = `decompmatrix'

    * Cleanup mata object if no save
    if "`matasave'" == "" mata mata drop `results'
end
