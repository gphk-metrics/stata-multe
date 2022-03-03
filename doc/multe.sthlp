{smcl}
{* *! version 0.1.3 02Mar2022}{...}
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
Run multiple IV regressions:

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
{synopt :{opth mata:save(str)}} Save mata object with results.
{p_end}
{synopt :{opth xx(str)}} Option xx.
{p_end}

{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
Alpha package for multiple treatment effects regression.

{marker example}{...}
{title:Examples}

{phang2}{cmd:. clear                   }{p_end}
{phang2}{cmd:. set seed 1729           }{p_end}
{phang2}{cmd:. set obs 1000            }{p_end}
{phang2}{cmd:. gen T = runiform() > 0.5}{p_end}
{phang2}{cmd:. gen W = mod(_n, 10)     }{p_end}
{phang2}{cmd:. gen Y = T + runiform()  }{p_end}
{phang2}{cmd:. multe Y T, control(W)   }{p_end}
