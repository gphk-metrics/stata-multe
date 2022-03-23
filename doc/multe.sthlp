{smcl}
{* *! version 0.1.5 23Mar2022}{...}
{viewerdialog multe "dialog multe"}{...}
{vieweralsosee "[R] multe" "mansection R multe"}{...}
{viewerjumpto "Syntax" "multe##syntax"}{...}
{viewerjumpto "Description" "multe##description"}{...}
{viewerjumpto "Options" "multe##options"}{...}
{viewerjumpto "Examples" "multe##examples"}{...}
{title:Title}

{p2colset 5 17 17 2}{...}
{p2col :{cmd:multe} {hline 2}}Multiple Treatment Effects regression{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}
Multiple treatment effects regression.

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
{synopt :{opth mata:save(str)}} Save mata object with results.
{p_end}
{synopt :{opt gen:erate(options)}} Optionally save tau, lambda. See {it:{help multe##gen_options:generate options}}.
{p_end}

{marker gen_options}{...}
{syntab :Generate Options}
{synopt :{cmd:lambda}[{cmd:(}str{cmd:)}]} Save lambdas in dataset; optinally specity prefix (default {cmd:lambda}).
{p_end}
{synopt :{cmd:tau}[{cmd:(}str{cmd:)}]} Save taus in dataset; optinally specity prefix (default {cmd:tau}).
{p_end}

{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
Alpha package for multiple treatment effects regression.

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
{phang2}{cmd:. return list                                                   }{p_end}
{phang2}{cmd:. multe Y T, control(W) gen(lambda tau)                         }{p_end}
{phang2}{cmd:. multe Y T, control(W) gen(lambda(awesomeName) tau(coolerName))}{p_end}
{phang2}{cmd:. desc, full                                                    }{p_end}
