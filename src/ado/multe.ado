*! version 0.1.0 16Feb2022
*! Multiple Treatment Effects Regression
*! Based code and notes by Michal Kolesár <kolesarmi@googlemail dotcom>
*! Adapted for Stata by Mauricio Caceres Bravo <mauricio.caceres.bravo@gmail.com>

capture program drop multe
program multe, eclass
    version 14.1

    syntax varlist(numeric fv ts min=2 max=2) /// depvar treatment
           [if] [in] ,                        /// subset
        control(varname)                      /// control variable
    [                                         ///
        *                                     ///
    ]

    * Varlist must be depvar and indepvar
    gettoken depvar treatment: varlist

    * Mark if, in, and any missing values by obs
    local varlist `varlist' `control'
    marksample touse, strok

    * Copy to mata for mata fun
    mata Wm = designmatrix(st_data(., "`control'", "`touse'"))
    mata Y  = st_data(., "`depvar'",    "`touse'")
    mata X  = st_data(., "`treatment'", "`touse'")
end
