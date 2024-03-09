{smcl}
{* *! version 1.0.1 08Mar2024}{...}
{viewerdialog multe "dialog multe"}{...}
{vieweralsosee "[R] multe" "mansection R multe"}{...}
{viewerjumpto "Syntax" "multe##syntax"}{...}
{viewerjumpto "Description" "multe##description"}{...}
{viewerjumpto "Options" "multe##options"}{...}
{viewerjumpto "Examples" "multe##examples"}{...}
{title:Title}

{p2colset 5 14 14 2}{...}
{p2col :{cmd:multe} {hline 2}}Multiple Treatment Effects for Contamination Bias Diagnostics{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}
This package computes multiple treatment effects regression for contamination bias diagnostics, using the procedures from Goldsmith-Pinkham, Hull, and Kolesár (2024).

{p 8 15 2}
{cmd:multe}
{depvar}
[{varlist}]
{ifin}
[{it:{help multe##weight:weight}}]
{cmd:,}
{opth treat(varname)}
[{opth stratum(varname)} {it:{help multe##table_options:options}}]

{pstd}
with {it:depvar} the dependent variable, {opt treat(varname)} the multi-valued treatment, and at least one of {it:varlist} with controls or {opt stratum(varname)} with a stratifying variable must be specified.

{synoptset 18 tabbed}{...}
{marker table_options}{...}
{synopthdr}
{synoptline}
{syntab :Options}
{synopt :{opt treat:ment(varname)}} Name of multi-valued treatment variable.
{p_end}
{synopt :{opt strat:um(varname)}} Name of stratifying covariate.
{p_end}
{synopt :{opt cluster(varname)}} Name of clustering covariate.
{p_end}
{synopt :{opt cw:_uniform}} Whether the target weighting scheme should give all comparisons equal weight or draw from the marginal empirical treatment distribution.
{p_end}
{synopt :{opt mata:save(str)}} Name of mata object with results (default: multe_results).
{p_end}

{p2colreset}{...}
{p 4 6 2}

{marker weight}{...}
{p 4 6 2}
{opt aweight}s and {opt pweight}s are allowed.

{marker description}{...}
{title:Description}

{pstd}
{cmd:multe} computes five estimates of treatment effects in settings with multiple treatments and provides estimates of contamination bias. See the package README at {browse "https://github.com/gphk-metrics/stata-multe/blob/main/README.md":github.com/gphk-metrics/stata-multe} for a detailed description and examples.

{pstd}
Heteroskedasticity-robust standard errors (default) and heteroskedasticity-robust standard errors that treat the treatment propensity scores as known (oracle) are also reported; standard errors can be clustered via {opt cluster()}. Each unique value of the {it:treatment} variable is taken as a distinct level of the treatment. The lowest value of {it:treatment} is picked to be the control group. The {cmd:stratum()} variables can be numeric or string, but should define a series of categories. Finally, the estimates are also reported for an "overlap" sample, which drops covariates that have no variation within a treatment level as well as observations from stratum that do not contain at least one observation from each treatment level. Unlike the full sample, all estimates are identified in the "overlap" sample.

{marker options}{...}
{title:Options}

{dlgtab:Command Options}

{phang}{opth treatment(varname)} Specifies the name of the treatment variable. Each value is taken as a different reatment level.

{phang}{opth stratum(varname)} Not optional if no covariates are specified. This can be any variable, and each value is taken to define a stratum. If any level does not have all values of the treatment, those observations are dropped in the overlap sample.

{phang}{opth cluster(varname)} Clustering covariate.

{phang}{opt cw_uniform} For the CW estimator, should the target weighting scheme give all comparisons equal weight (default), or should it draw from the marginal empirical treatment distribution (if specified)?

{phang}{opth mata:save(str)} supplies an alternative name for the mata structure which stores all estimates and variables in mata (default name is "multe_results"). Note this is in addition to results stored in {cmd:e()}; see {it:{help multe##results:stored results}} below for details.

{marker example}{...}
{title:Example 1: Generated data}

{phang2}{cmd:. local nobs   1000                   }{p_end}
{phang2}{cmd:. local ktreat 5                      }{p_end}
{phang2}{cmd:. clear                               }{p_end}
{phang2}{cmd:. set seed 1729                       }{p_end}
{phang2}{cmd:. set obs `nobs'                      }{p_end}
{phang2}{cmd:. gen T = ceil(runiform() * `ktreat') }{p_end}
{phang2}{cmd:. gen W = mod(_n, 10)                 }{p_end}
{phang2}{cmd:. gen Y = T + runiform()              }{p_end}
{phang2}{cmd:. multe Y, treat(T) strat(W)          }{p_end}
{phang2}{cmd:. ereturn list                        }{p_end}

{title:Example 2: Fryer and Levitt (2013)}

{pstd}The data for this example can be downloaded with the {cmd:multe} package by specifying the option {cmd:all} (i.e. {it:net install multe, all}).

{phang2}{cmd:. use example_fryer_levitt.dta, clear                                   }{p_end}
{phang2}{cmd:. local controls i.age_24 female i.siblings i.family_structure          }{p_end}
{phang2}{cmd:. multe std_iq_24 `controls' [w=W2C0], treat(race) stratum(SES_quintile)}{p_end}
{phang2}{cmd:. ereturn list    }{p_end}
{phang2}{cmd:. multe, est(ATE) }{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:multe} stores the following in {cmd:e()}:

{synoptset 23 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}

{p2col 5 23 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:multe}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(wtype)}}weight type{p_end}
{synopt:{cmd:e(wexp)}}weight expression{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(treatment)}}name of multi-valued treatment{p_end}
{synopt:{cmd:e(controls)}}name of control variables{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(mata)}}name of mata object where results are stored (see below){p_end}

{p2col 5 23 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector for requested estimates (default is partly linear model).{p_end}
{synopt:{cmd:e(V)}}cvariance matrix corresponding to the coefficient vector.{p_end}
{synopt:{cmd:e(fullmatrix)}}matrix of results with all five estimates and corresponding SEs.{p_end}
{synopt:{cmd:e(diffmatrix)}}matrix of differences between partly linear estimates and the rest of the estimates, with corresponding SEs.{p_end}
{synopt:{cmd:e(full_beta)}}coefficients for each of the five estimates.{p_end}
{synopt:{cmd:e(full_pop_se)}}SEs for each of the five estimates.{p_end}
{synopt:{cmd:e(full_oracle_se)}}Oracle SEs for each of the five estimates.{p_end}
{synopt:{cmd:e(full_diff)}}differences between partly linear estimates and the rest of the estimates.{p_end}
{synopt:{cmd:e(full_diff_se)}}SEs for each of the differences.{p_end}
{synopt:{cmd:e(overlapdiffmatrix)}}only if overlap sample is computed; analogous to {cmd:e(fullmatrix)}.{p_end}
{synopt:{cmd:e(overlapmatrix)}}only if overlap sample is computed; analogous to {cmd:e(diffmatrix)}.{p_end}
{synopt:{cmd:e(overlap_diff)}}only if overlap sample is computed; analogous to {cmd:e(full_beta)}{p_end}
{synopt:{cmd:e(overlap_diff_se)}}only if overlap sample is computed; analogous to {cmd:e(full_pop_se)}{p_end}
{synopt:{cmd:e(overlap_oracle_se)}}only if overlap sample is computed; analogous to {cmd:e(full_oracle_se}){p_end}
{synopt:{cmd:e(overlap_pop_se)}}only if overlap sample is computed; analogous to {cmd:e(full_diff)}{p_end}
{synopt:{cmd:e(overlap_beta}})only if overlap sample is computed; analogous to {cmd:e(full_diff_se)}{p_end}

{p2col 5 23 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker mata}{...}
{pstd}
In addition, the following data are available in {cmd:e(mata)} (default name: multe_results):

        real scalar has_overlap
            whether estimates for the overlap sample were computed

        real scalar Tk
            number of treatment levels

        string rowvector Tlevels
            treatment level names or labels

        string scalar Tvar
            name of treatment variable

        string scalar Yvar
            name of dependent variable

        real matrix full.N
            number of observations in full sample

        real matrix full.touse
            vector if 1s and 0s flagging used observations

        real matrix full.estA
            matrix with each of the give estimates

        real matrix full.estB
            matrix with the difference between PL and each of the other estimates

        real matrix full.seP
            SEs for estA

        real matrix full.seO
            Oracle SEs for estA

        real matrix full.seB
            SEs for estB

        real matrix full.Vpop_PL
            Covariance for PL

        real matrix full.Vpop_OWN
            Covariance for OWN

        real matrix full.Vpop_ATE
            Covariance for ATE

        real matrix full.Vpop_EW
            Covariance for EW

        real matrix full.Vpop_CW
            Covariance for CW

        real matrix full.Vo_OWN
            Oracle covariance OWN

        real matrix full.Vo_ATE
            Oracle covariance ATE

        real matrix full.Vo_EW
            Oracle covariance EW

        real matrix full.Vo_CW
            Oracle covariance CW

        real matrix full.Vdiff_OWN
            Covariance for differences between PL and OWN

        real matrix full.Vdiff_ATE
            Covariance for differences between PL and ATE

        real matrix full.Vdiff_EW
            Covariance for differences between PL and EW

        real matrix full.Vdiff_CW
            Covariance for differences between PL and CW

        real matrix full.LM.stat
            LM stat

        real matrix full.LM.df
            degrees of freedom

        real matrix full.LM.pval
            p-value

        struct full.Wa
            analogous to full.LM; contains stat, df, and pval for Wald test

        struct overlap
            analogous to full; contains all the same data for the overlap sample

{marker references}{...}
{title:References}

{pstd}
Goldsmith-Pinkham, Paul, Peter Hull, Michal Koles{c a'}r (2024): "Contamination Bias in Linear Regressions" Working Paper
