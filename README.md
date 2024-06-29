MulTE
=====

Multiple Treatment Effects

`version 1.1.0 09Mar2024` | [Installation](#installation) | [Usage](#usage) | [Examples](#examples)

This package implements contamination bias diagnostics in Stata using procedures from
[Goldsmith-Pinkham, Hull, and Kolesár (2024)](https://arxiv.org/abs/2106.05024),
based on R's [multe](https://github.com/kolesarm/multe?tab=readme-ov-file) package.

### Installation

From Stata, latest version:

```stata
local github "https://raw.githubusercontent.com"
cap noi net uninstall multe
net install multe, from(`github'/gphk-metrics/stata-multe/main/)
```

You can also clone or download the code manually, e.g. to
`stata-multe-main`, and install from a local folder:

```stata
cap noi net uninstall multe
net install multe, from(`c(pwd)'/stata-multe-main)
```

We typically recommend installing the latest version; however, if you need a specific
tagged version (e.g. for replication purposes):

```stata
local github "https://raw.githubusercontent.com"
local multeversion "1.1.0"
cap noi net uninstall multe
net install multe, from(`github'/gphk-metrics/stata-multe/`mutleversion'/)
```

See the [tags](https://github.com/gphk-metrics/stata-multe/tags) page for all
available tagged versions.

### Usage

```
multe depvar [controls] [w=weight], treat(varname) [stratum(varname) options]
help multe
```

At least one of `controls` or `stratum` must be specified.

### Examples

The examples here are copied from the [vignette for the R package](https://github.com/kolesarm/multe/blob/master/doc/multe.pdf) and adapted for Stata. The vignette also includes an overview of the methods used in this package.

See [Fryer and Levitt (2013)](https://www.aeaweb.org/articles?id=10.1257/aer.103.2.981) for a description of the data. First, we fit a regression of test scores on a race dummy (the treatment of interest) and a few controls, using sampling weights.

```stata
use test/example_fryer_levitt.dta, clear
* Regress IQ at 24 months on race indicators and baseline controls
multe std_iq_24 i.age_24 female [w=W2C0], treat(race) stratum(SES_quintile)
```

The function posts the estimates from the partly linear model (OLS with robust
SE) and reports a table with alternative estimates free of contamination bias.

```stata
Partly Linear Model Estimates (full sample)     Number of obs     =      8,806

------------------------------------------------------------------------------
   std_iq_24 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
        race |
      Black  |  -.2574154   .0281244    -9.15   0.000    -.3125383   -.2022925
   Hispanic  |  -.2931465    .025962   -11.29   0.000    -.3440311    -.242262
      Asian  |   -.262108    .034262    -7.65   0.000    -.3292603   -.1949557
      Other  |  -.1563373   .0369127    -4.24   0.000    -.2286848   -.0839898
------------------------------------------------------------------------------

Alternative Estimates on Full Sample:

             |      PL      OWN      ATE       EW       CW
-------------+---------------------------------------------
       Black |  -.2574   -.2482   -.2655    -.255   -.2604
          SE |  .02812   .02906   .02983   .02888   .02925
    Hispanic |  -.2931   -.2829   -.2992   -.2862   -.2944
          SE |  .02596   .02673   .02988    .0268   .02792
       Asian |  -.2621   -.2609   -.2599   -.2611   -.2694
          SE |  .03426   .03432   .04177   .03433   .04751
       Other |  -.1563   -.1448   -.1503   -.1447   -.1522
          SE |  .03691   .03696   .03594   .03684   .03698

P-values for null hypothesis of no propensity score variation:
Wald test:  2.1e-188
  LM test:  8.8e-197
```

In particular, the package computes the following estimates in addition to the
partly linear model (PL):

- OWN: The own treatment effect component of the PL estimator that subtracts an estimate of the contamination bias
- ATE: The unweighted average treatment effect, implemented using regression that includes interactions of covariates with the treatment indicators
- EW: Weighted ATE estimator based on easiest-to-estimate weighting (EW) scheme, implemented by running one-treatment-at-a-time regressions.
- CW: Weighted ATE estimator using easiest-to-estimate common weighting (CW) scheme, implemented using weighted regression.

For more precise definitions, see the methods section of the [R vignette](https://github.com/kolesarm/multe/blob/master/doc/multe.pdf). Note any combination of estimates can be posted after `multe` has run:

```stata
. multe, est(ATE)

ATE Estimates (full sample)                     Number of obs     =      8,806

------------------------------------------------------------------------------
   std_iq_24 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
        race |
      Black  |  -.2655242   .0298285    -8.90   0.000    -.3239869   -.2070615
   Hispanic  |  -.2992432   .0298811   -10.01   0.000     -.357809   -.2406774
      Asian  |  -.2598644    .041769    -6.22   0.000    -.3417301   -.1779987
      Other  |  -.1502789   .0359405    -4.18   0.000     -.220721   -.0798368
------------------------------------------------------------------------------
```

In this example, the propensity score varies significantly with covariates, as indicated by the p-values of the Wald and LM tests. Including many controls may result in overlap failure, as the next example demonstrates:

```stata
. local controls i.age_24 female i.siblings i.family_structure

. multe std_iq_24 `controls' [w=W2C0], treat(race) stratum(SES_quintile)
(analytic weights assumed)
The following variables have no within-treatment variation
and are dropped from the overlap sample:
    6.siblings

Partly Linear Model Estimates (full sample)     Number of obs     =      8,806

------------------------------------------------------------------------------
   std_iq_24 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
        race |
      Black  |  -.2437653   .0307696    -7.92   0.000    -.3040726    -.183458
   Hispanic  |  -.2928037   .0259029   -11.30   0.000    -.3435725   -.2420348
      Asian  |   -.273888   .0341814    -8.01   0.000    -.3408823   -.2068938
      Other  |  -.1519798   .0368914    -4.12   0.000    -.2242857    -.079674
------------------------------------------------------------------------------

Alternative Estimates on Full Sample:

             |      PL      OWN      ATE       EW       CW
-------------+---------------------------------------------
       Black |  -.2438   -.2043   -.2482    -.218   -.2415
          SE |  .03077   .03321   .03553   .03276   .03853
    Hispanic |  -.2928   -.2801   -.2878    -.285   -.3001
          SE |   .0259   .02671      .03   .02671   .02984
       Asian |  -.2739   -.2836   -.2742   -.2839   -.2863
          SE |  .03418    .0343   .04195   .03426   .04549
       Other |   -.152   -.1277        .   -.1295   -.1452
          SE |  .03689   .03736        .   .03713   .03817

P-values for null hypothesis of no propensity score variation:
Wald test:  1.3e-275
  LM test:  2.2e-245

Alternativ Estimates on Overlap Sample:
             |      PL      OWN      ATE       EW       CW
-------------+---------------------------------------------
       Black |  -.2458   -.2062   -.2503   -.2199   -.2436
          SE |  .03074   .03318   .03532    .0327   .03824
    Hispanic |  -.2932   -.2809    -.288   -.2858   -.2987
          SE |  .02595   .02676   .03001   .02678   .02987
       Asian |  -.2741   -.2839    -.274   -.2841   -.2884
          SE |  .03417   .03429   .04187   .03426   .04563
       Other |  -.1511    -.127   -.1392   -.1289   -.1459
          SE |  .03688   .03735   .03618   .03711   .03853
```

The issue is that no observations with 6 siblings have race equal to "Other":

```stata
. tab race if siblings == 6, mi

       race |      Freq.     Percent        Cum.
------------+-----------------------------------
      White |         18       40.00       40.00
      Black |         10       22.22       62.22
   Hispanic |         12       26.67       88.89
      Asian |          5       11.11      100.00
------------+-----------------------------------
      Total |         45      100.00
```

Thus, the ATE estimator comparing other to white is not identified. The package drops observations with 6 siblings from the sample to form an "overlap sample," where the all estimators are identified. The overlap sample drops all observations from strata that do not have all levels of the treatment as well as all control variables that do not have all levels of the treatment. In this example `siblings` is a control so the dummy for 6 siblings is dropped; however, if this happened for a level of `SES_quintile` then all the observations associated with that level would be dropped.

For a researcher who wants to check whether there is a significant difference between the PL estimator and the other estimators, `e(diffmatrix)` and `e(overlapdiffmatrix)` report the differences between the estimates in the full sample and the overlap sample, respectively, with the corresponding standard errors.

```stata
. matlist e(diffmatrix), format(%7.4g)

             |     OWN      ATE       EW       CW
-------------+------------------------------------
       Black | -.03947    .0044  -.02573  -.00229
          SE |  .01204    .0173   .01027    .0205
    Hispanic | -.01269  -.00496  -.00783   .00725
          SE |  .00571   .01314   .00497   .01245
       Asian |  .00975   .00026   .01001   .01238
          SE |  .00487    .0246    .0048   .02814
       Other |  -.0243        .  -.02246  -.00679
          SE |  .00738        .   .00684   .01326

. matlist e(overlapdiffmatrix), format(%7.4g)

             |     OWN      ATE       EW       CW
-------------+------------------------------------
       Black | -.03958   .00453  -.02592  -.00218
          SE |  .01201   .01714   .01023   .02023
    Hispanic | -.01229  -.00523  -.00738   .00554
          SE |  .00569   .01319   .00494   .01247
       Asian |  .00978  -8.7e-05      .01   .01433
          SE |  .00486   .02454    .0048    .0282
       Other | -.02404  -.01189  -.02221   -.0052
          SE |  .00739   .01098   .00684   .01337
```

We see statistically significant difference between the OWN and PL estimate (i.e. significant contamination bias) for all races, both in the full sample and in the overlap sample. Differences can also be posted after the main estimation:

```
. multe, est(OWN) diff

Own Treatment Effects Estimates (full sample)   Number of obs     =      8,806

------------------------------------------------------------------------------
   std_iq_24 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
        race |
      Black  |  -.0394671   .0120418    -3.28   0.001    -.0630686   -.0158656
   Hispanic  |  -.0126872   .0057131    -2.22   0.026    -.0238847   -.0014898
      Asian  |    .009752   .0048651     2.00   0.045     .0002166    .0192874
      Other  |  -.0242985   .0073833    -3.29   0.001    -.0387696   -.0098274
------------------------------------------------------------------------------

. multe, est(OWN) diff overlap

Own Treatment Effects Estimates (overlap sample)

                                                Number of obs     =      8,806

------------------------------------------------------------------------------
   std_iq_24 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
        race |
      Black  |  -.0395843   .0120096    -3.30   0.001    -.0631227   -.0160458
   Hispanic  |    -.01229   .0056924    -2.16   0.031    -.0234468   -.0011332
      Asian  |    .009776   .0048602     2.01   0.044     .0002502    .0193019
      Other  |  -.0240377   .0073857    -3.25   0.001    -.0385134    -.009562
------------------------------------------------------------------------------
```

The package also computes "oracle" standard errors, in addition to the usual
standard errors reported above. These can be accessed by adding the `oracle`
option to `multe` after the main estimation for `ATE`, `EW`, `CW`:

```stata
. multe, est(ATE) oracle

ATE Estimates (full sample)                     Number of obs     =      8,806

------------------------------------------------------------------------------
             |               Oracle
   std_iq_24 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
        race |
      Black  |  -.2481627   .0354964    -6.99   0.000    -.3177343   -.1785911
   Hispanic  |  -.2878414   .0299246    -9.62   0.000    -.3464926   -.2291902
      Asian  |  -.2741527   .0418856    -6.55   0.000    -.3562469   -.1920584
      Other  |          0  (omitted)
------------------------------------------------------------------------------
```

These oracle standard errors don't account for estimation error in the
propensity score, in contrast to the default standard errors (note since ATE
was not identified for Other, Stata interprets it as an omitted category
for this estimator). Last, Specifying the `cluster()` argument allows for
computation of clustered standard errors:

```
. local controls i.age_24 female

. local options  treat(race) stratum(SES_quintile)

. multe std_iq_24 `controls' [w=W2C0], `options' cluster(interviewer_ID_24)
(analytic weights assumed)

Partly Linear Model Estimates (full sample)     Number of obs     =      8,806

------------------------------------------------------------------------------
             |               Cluster
   std_iq_24 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
        race |
      Black  |  -.2574154   .0411564    -6.25   0.000    -.3380804   -.1767504
   Hispanic  |  -.2931465   .0440827    -6.65   0.000     -.379547   -.2067461
      Asian  |   -.262108   .0521067    -5.03   0.000    -.3642353   -.1599807
      Other  |  -.1563373   .0403717    -3.87   0.000    -.2354645   -.0772102
------------------------------------------------------------------------------

Alternative Estimates on Full Sample:

             |      PL      OWN      ATE       EW       CW
-------------+---------------------------------------------
       Black |  -.2574   -.2482   -.2655    -.255   -.2604
          SE |  .04116   .04248   .04099    .0422   .04199
    Hispanic |  -.2931   -.2829   -.2992   -.2862   -.2944
          SE |  .04408   .04541    .0495   .04572   .04738
       Asian |  -.2621   -.2609   -.2599   -.2611   -.2694
          SE |  .05211   .05234   .06194   .05232   .06755
       Other |  -.1563   -.1448   -.1503   -.1447   -.1522
          SE |  .04037   .04161   .04099   .04145   .04278

P-values for null hypothesis of no propensity score variation:
Wald test:  9.4e-133
  LM test:  6.79e-07
```
