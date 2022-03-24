{smcl}
{* *! version 0.2.0 24Mar2022}{...}
{viewerdialog multe "dialog multe"}{...}
{vieweralsosee "[R] multe" "mansection R multe"}{...}
{viewerjumpto "Syntax" "multe##syntax"}{...}
{viewerjumpto "Description" "multe##description"}{...}
{viewerjumpto "Options" "multe##options"}{...}
{viewerjumpto "Examples" "multe##examples"}{...}
{title:Title}

{p2colset 5 14 14 2}{...}
{p2col :{cmd:multe} {hline 2}}Multiple Treatment Effects regression{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}
Multiple treatment effects regression. TODO: xx short description.

{p 8 15 2}
{cmd:multe}
{depvar}
{it:treatment}
{ifin}
{opth control(varname)}
[{cmd:,} {it:{help multe##table_options:options}}]

{synoptset 18 tabbed}{...}
{marker table_options}{...}
{synopthdr}
{synoptline}
{syntab :Options}
{synopt :{opt vce(str)}} Type of standard errors to print: "" (default) or "oracle" (both are computed internally).
{p_end}
{synopt :{opth mata:save(str)}} Name of mata object with results (default: MulTEResults).
{p_end}
{synopt :{opt gen:erate(options)}} Optionally save tau, lambda. See {it:{help multe##gen_options:generate options}}.
{p_end}

{marker gen_options}{...}
{syntab :Generate Options}
{synopt :{cmd:lambda}[{cmd:(}str{cmd:)}]} Save lambdas in dataset; optinally specity prefix (default: {cmd:lambda}).
{p_end}
{synopt :{cmd:tau}[{cmd:(}str{cmd:)}]} Save taus in dataset; optinally specity prefix (default: {cmd:tau}).
{p_end}

{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
Alpha package for multiple treatment effects regression. TODO: xx longer description.

{marker options}{...}
{title:Options}

{dlgtab:Command Options}

{phang}{opth vce(str)} TODO: xx

{phang}{opth mata:save(str)} TODO: xx

{phang}{opt gen:erate(options)} TODO: xx

{dlgtab:Generate Options}

{phang}{cmd:lambda}[{cmd:(}str{cmd:)}] TODO: xx

{phang}{cmd:tau}[{cmd:(}str{cmd:)}] TODO: xx

{marker example}{...}
{title:Examples}

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
{phang2}{cmd:. desc, full                                                    }{p_end}

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
{synopt:{cmd:e(V)}}diagonal matrix of squared standard errors{p_end}
{synopt:{cmd:e(estimates)}}matrix of coefficients, including SE and orable SE.{p_end}
{synopt:{cmd:e(decomposition)}}decomposition matrix (beta, own effect, contamination bias, minimum bias, maximum bias){p_end}

{p2col 5 23 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{pstd}
In addition, the following data are available in {cmd:e(mata)} (default name: MulTEResults):

        real scalar estimates.n
            number of observations

        real scalar estimates.k
            number of treatment levels

        real matrix estimates.est
            (k-1) by 3 matrix of coefficients

        real matrix estimates.se_po
            (k-1) by 3 matrix of standard errors

        real matrix estimates.se_or
            (k-1) by 3 matrix of standard errors

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
            print all estimates to console (default 6 dignificant digits)

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
TODO: xx Add reference?

