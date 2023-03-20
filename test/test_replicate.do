global star "~/Dropbox/GPH_ExaminerDesign/Applications/STAR"
global data	"${star}/Data"
* global star "/Users/jchang42/Dropbox (Brown)/GPH_ExaminerDesign/Applications/STAR"
* global data "${star}/Data"

capture program drop multe_replicate_tests
program multe_replicate_tests
    syntax, [Verbose]
    if "`verbose'" != "" local noisily noisily
    cap findfile example_star.dta
    if ( _rc ) {
        qui multe_replicate_load
    }
    else {
        use `"`r(fn) '"', clear
    }
    cap `noisily' multe score treatment, control(school) decomp
    if ( _rc != 0 ) {
        disp "(multe test fail): multe run on STAR data failed with _rc = `rc1'"
        exit 9
    }
    else {
        disp "(multe test success): multe ran successfully on STAR data"
    }
	* output, treatment(treatment) matasave(multe_results) ///
    *     outpath("${star}/Output/tables/test") // pick your outpath
    * mata mata desc
end

capture program drop multe_replicate_load
program multe_replicate_load
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
    label values treatment treatmentlab
    multe_test_factor school    = schoollab,  replace
    multe_test_factor teacher   = teacherlab, replace
end

capture program drop multe_test_factor
program multe_test_factor
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
