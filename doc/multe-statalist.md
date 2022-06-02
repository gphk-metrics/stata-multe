## multe: a new package to examine and provide solutions for contamination bias in your multiple treatment effects regressions with saturated group control
Thanks to Kit Baum, a new package -multe- is now available on SSC. Using Stata 14.1 or later, install via:

```
ssc install multe
```

If your research setting has multiple treatments and a saturated group control within which treatment is as-good-as-randomly assigned, treatment propensity scores which vary across saturated group, and heterogeneous treatment effects, you may be getting treatment effect estimates with contamination bias. Examples of such settings include multi-armed RCTs and two-way fixed effects.

**multe** can help you measure and avoid contamination bias in such settings.

**Avoids contamination bias:** Based on Goldsmith-Pinkham, Hull, and Kolesár (2022), the **multe** command computes three types of weighted-average treatment effects that avoid contamination bias. These are the equal-weighted averages (i.e. average treatment effect), treatment-specific variance-weighted averages (as in Angrist et al. (1998)), and comparably efficiently-weighted averages (as in Goldsmith-Pinkham et al. (2022)).

**Measures contamination bias:** Our command also measures the extent to which contamination bias affects treatment effect estimates based on a partially linear model (i.e. the coefficients on X when regressing Y on X and W, given an outcome Y, a vector of treatment indicators X, and a vector of saturated group indicators W). It decomposes such treatment effect estimates into own-effect and contamination bias terms.

For an example of usage, see the Stata helpfile. Example data from Project STAR are available for download directly with the package (via `ssc install multe, all`) or through our online repository [here](https://raw.githubusercontent.com/gphk-metrics/stata-multe/ab353845e9cc4d3f30563c345342daff2ee1dec8/test/example_star.dta).
