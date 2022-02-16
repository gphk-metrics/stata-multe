global star "~/Dropbox/GPH_ExaminerDesign/Applications/STAR"
global data	"${star}/Data"

capture program drop main
program main

    * do ../src/ado/multe.ado
    * do ../src/mata/multe.mata
    loat_test_data
    multe score treatment, control(school)
    mata mata desc
end

capture program drop loat_test_data
program loat_test_data
    use `"${data}/STARgk_Lambdas.dta"', clear

    gen treatment = "regular"
    replace treatment = "small" if gksmall == 1
    replace treatment = "aide"  if gkaide  == 1

    keep gkaggregate2 treatment gkschid gktchid
    rename (gkaggregate2 treatment gkschid gktchid) (score treatmentlab schoollab teacherlab)
    gen treatment = 1 
    replace treatment = 2 if treatmentlab == "small"
    replace treatment = 3 if treatmentlab == "aide"
    label define treatmentlab 1 "regular" 2 "small" 3 "aide"
    factor school    = schoollab,  replace
    factor teacher   = teacherlab, replace
end

capture program drop factor
program factor
    syntax anything(equalok), [replace]
    gettoken gen varlist: anything, p(=)
    gettoken _ varlist: varlist
    confirm var `varlist'
    foreach var of varlist `varlist' {
        gettoken outvar gen: gen
        local i = 0
        local lbldef
        glevelsof `var', silent loc(`var')
        foreach level of local `var' {
            mata st_local("lbldef`++i'", `"`i' "`level'""')
            mata st_local("lbldef",  st_local("lbldef") + st_local("lbldef`i'"))
        }
        if "`replace'" != "" cap label drop `var'
        label define `var' `lbldef'
        gegen `outvar' = group(`var'), `replace'
        label values `outvar' `var'
    }
end

main
