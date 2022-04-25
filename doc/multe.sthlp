{smcl}
{* *! version 0.2.3 22Apr2022}{...}
{viewerdialog multe "dialog multe"}{...}
{vieweralsosee "[R] multe" "mansection R multe"}{...}
{viewerjumpto "Syntax" "multe##syntax"}{...}
{viewerjumpto "Description" "multe##description"}{...}
{viewerjumpto "Options" "multe##options"}{...}
{viewerjumpto "Examples" "multe##examples"}{...}
{title:Title}

{p2colset 5 14 14 2}{...}
{p2col :{cmd:multe} {hline 2}}Multiple Treatment Effects regression with saturated group control{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}
Multiple treatment effects regression.
Given a multi-valued treatment, a saturated group variable (or a {varlist} which will be used to create a single, saturated group variable), and a dependent variable,
{cmd:multe} computes equal-weighted (ATE), variance-weighted, efficiently-weighted treatment effects estimates and contamination bias decomposition as in Goldsmith-Pinkham et al. (2022).

{p 8 15 2}
{cmd:multe}
{depvar}
{it:treatment}
{ifin}
{cmd:,}
{opth control(varlist)}
[{it:{help multe##table_options:options}}]

{synoptset 18 tabbed}{...}
{marker table_options}{...}
{synopthdr}
{synoptline}
{syntab :Options}
{synopt :{opt vce(str)}} Type of standard errors to print: "" (heteroskedasticity-robust) or "oracle" (heteroskedasticity-robust, treating propensity scores as known).
{p_end}
{synopt :{opt mata:save(str)}} Name of mata object with results (default: multe_results).
{p_end}
{synopt :{opt gen:erate(options)}} Optionally save tau (saturated group-specific treatment effects), lambda (implicit ATE regression weights). See {it:{help multe##gen_options:generate options}}.
{p_end}

{marker gen_options}{...}
{syntab :Generate Options}
{synopt :{cmd:lambda}[{cmd:(}str{cmd:)}]} Save lambdas in dataset; optionally specify prefix (default: {cmd:lambda}).
{p_end}
{synopt :{cmd:tau}[{cmd:(}str{cmd:)}]} Save taus in dataset; optionally specify prefix (default: {cmd:tau}).
{p_end}

{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:multe} computes three types of weighted-average treatment effects in settings with multiple treatments and a saturated group control, within which treatment is as-good-as-randomly
assigned. These are the equal-weighted averages (i.e. average treatment effects), treatment-specific variance-weighted averages (as in Angrist et al. (1998)), and comparable efficiently
weighted averages (as in Goldsmith-Pinkham et al. (2022)).

{pstd}
It also computes and saves a decomposition of regression estimates of treatment effects into an own-effect weighted average and a contamination bias term, following
Goldsmith-Pinkham et al. (2022). It also provides an option to save the implicit ATE regression weights (lambda) and/or the saturated group-specific treatment effects (tau) as variables
(see {it:{help multe##options:options}} for details). Results are saved in {cmd:e()} (see {it:{help multe##results:stored results}} below for details).

{pstd}
Heteroskedasticity-robust standard errors (default) and heteroskedasticity-robust standard errors that treat the treatment propensity scores as known (oracle) are also reported.

{pstd}
The {depvar} and {it:treatment} can be any numeric variables. However, each unique value of the {it:treatment} variable is taken as a distinct level of the treatment. The lowest value of
{it:treatment} is picked to be the control group. The {cmd:control} variables can be numeric or string, but should define a series of categories ({cmd:multe} will turn the controls into a
single, saturated group variable). Groups which do not satisfy overlap (i.e. the number of unique treatment levels in that group is less than the total number of unique treatment levels)
will be dropped.

{pstd}
For a detailed theoretical discussion of calculations done by {cmd:multe}, see Goldsmith-Pinkham et al. (2022). Examples are provided below, including the Project STAR
example used in Goldsmith-Pinkham et al. (2022).

{marker options}{...}
{title:Options}

{dlgtab:Command Options}

{phang}{opth vce(str)} specifies the type of standard errors to print. The default "" is heteroskedasticity-robust, and "oracle" specifies heteroskedasticty-robust standard errors that treat the
propensity score for each treatment level as known.

{phang}{opth mata:save(str)} supplies an alternative name for the mata structure which stores all estimates and variables in mata (default name is "multe_results").
Note this is in addition to results stored in {cmd:e()}; see {it:{help multe##results:stored results}} below for details.

{phang}{opt gen:erate(options)} specifies whether to save lambda and/or tau (as defined above) as variables. The user can optionally specify the names of these two sets of variables via the
options {cmd:lambda}[{cmd:(}str{cmd:)}] and {cmd:tau}[{cmd:(}str{cmd:)}]. For example, {cmd:gen(lambda tau)} would generate both with default names, while {cmd:gen(lambda(lname) tau(tname))}
would generate them with custom names.

{dlgtab:Generate Options}

{phang}{cmd:lambda}[{cmd:(}str{cmd:)}] saves the set of implicit ATE regression weights as variables and optionally specifies an alternative prefix. The default prefix is "lambda".

{phang}{cmd:tau}[{cmd:(}str{cmd:)}] saves the saturated group-specific treatment effects as variables and optionally specifies an alternative prefix. The default prefix is "tau".

{marker example}{...}
{title:Example 1: Generated data}

{phang2}{cmd:. local nobs   1000                                             }{p_end}
{phang2}{cmd:. local ktreat 5                                                }{p_end}
{phang2}{cmd:. clear                                                         }{p_end}
{phang2}{cmd:. set seed 1729                                                 }{p_end}
{phang2}{cmd:. set obs `nobs'                                                }{p_end}
{phang2}{cmd:. gen T = ceil(runiform() * `ktreat')                           }{p_end}
{phang2}{cmd:. gen W = mod(_n, 10)                                           }{p_end}
{phang2}{cmd:. gen Y = T + runiform()                                        }{p_end}
{phang2}{cmd:. multe Y T, control(W)                                         }{p_end}
{phang2}{cmd:. ereturn list                                                  }{p_end}
{phang2}{cmd:. multe, vce(oracle)                                            }{p_end}
{phang2}{cmd:. multe Y T, control(W) gen(lambda tau)                         }{p_end}
{phang2}{cmd:. multe Y T, control(W) gen(lambda(awesomeName) tau(coolerName))}{p_end}
{phang2}{cmd:. mata `e(mata)'.decomposition.print(1)                         }{p_end}
{phang2}{cmd:. desc, full                                                    }{p_end}

{title:Example 2: Project STAR}

{pstd}The data for this example can be downloaded with the {cmd:multe} package by specifying the option {cmd:all} (e.g. {it:ssc install multe, all}) or from our online repository {browse "https://raw.githubusercontent.com/gphk-metrics/stata-multe/ab353845e9cc4d3f30563c345342daff2ee1dec8/test/example_star.dta":here}.

{phang2}{cmd:. use example_star.dta, clear                                         }{p_end}
{phang2}{cmd:. multe score treatment, control(school)                              }{p_end}
{phang2}{cmd:. ereturn list                                                        }{p_end}
{phang2}{cmd:. multe, vce(oracle)                                                  }{p_end}
{phang2}{cmd:. multe score treatment, control(school) gen(lambda(M_) tau(tauhat_)) }{p_end}
{phang2}{cmd:. desc, full                                                          }{p_end}

{pstd}After obtaining the implicit equal-weighted regression weights (lambda) and group-specific treatment effects (tau), you can calculate
the correlations to get a sense of how much contamination bias might affect ATE estimates:

{phang2}{cmd:. forval i=1/2 {c -(}      }{p_end}
{phang2}{cmd:. forval j=1/2 {c -(}      }{p_end}
{phang2}{cmd:. corr tauhat_`j' M_`i'`j' }{p_end}
{phang2}{cmd:. {c )-}                   }{p_end}
{phang2}{cmd:. {c )-}                   }{p_end}

{pstd}You can also optionally specify an alternative name for the mata struct which contains store results (see {it:{help multe##mata:Stored mata results}}).

{phang2}{cmd:. multe score treatment, control(school) matasave(matastructname)}{p_end}
{phang2}{cmd:. mata mata desc}{p_end}

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
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(control)}}name of control variable{p_end}
{synopt:{cmd:e(vce)}}only populated when {opt vce(oracle)} is requested{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V} ({bf:NB}: V is a diagonal matrix of squared SEs.){p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}
{synopt:{cmd:e(mata)}}name of mata object where results are stored (see below){p_end}

{p2col 5 23 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}block diagonal matrix of three covariance matrices corresponding to the ATE estimates, one-at-a-time estimates, and efficiently-weighted
estimates. Note the covariance matrix for the one-at-a-time estimates only reports diagonal terms.{p_end}
{synopt:{cmd:e(estimates)}}matrix of coefficients, including SE and orable SE.{p_end}
{synopt:{cmd:e(decomposition)}}decomposition matrix (beta, own effect, contamination bias, minimum bias, maximum bias){p_end}

{p2col 5 23 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker mata}{...}
{pstd}
In addition, the following data are available in {cmd:e(mata)} (default name: multe_results):

        real scalar estimates.n
            number of observations

        real scalar estimates.k
            number of treatment levels

        real matrix estimates.est
            (k-1) by 3 matrix of coefficients

        real matrix estimates.se_po
            (k-1) by 3 matrix of heteroskedasticity-robust standard errors

        real matrix estimates.se_or
            (k-1) by 3 matrix of heteroskedasticity-robust standard errors treating treatment propensity scores as known

        real vector estimates.Tvalues
            k by 1 vector with treatment levels

        string vector estimates.Tlabels
            k by 1 vector with treatment labels

        string scalar estimates.Tvar
            treatment variable name

        string scalar estimates.Yvar
            outcome variable name

        string vector estimates.colnames
            column names for printing/saving estimates

        void estimates.print(| real scalar digits)
            print all estimates to console (default 6 significant digits)

        void estimates.save(string scalar outmatrix)
            save estimates to outmatrix

        real matrix decomposition.est
            (k-1) by 5 matrix with coefficient decomposition

        real matrix decomposition.se
            (k-1) by 5 matrix with decomposition SEs

        string vector decomposition.tauhat_names
            variable names for tauhat

        string vector decomposition.lambda_names
            variable names for lambda

        real matrix decomposition.tauhat(real matrix Wm)
            compute tauhat from matrix of dummies corresponding to control

        real matrix decomposition.lambda(real matrix Wm)
            compute lambda from matrix of dummies corresponding to control

{marker references}{...}
{title:References}

{pstd}
Goldsmith-Pinkham, Paul, Peter Hull, Michal Koles{c a'}r (2022): "Contamination Bias in Linear Regressions" Working Paper

