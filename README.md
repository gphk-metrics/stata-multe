MulTE
=====

Multiple Treatment Effects regression

`version 0.2.2 13Apr2022` | [Installation](#installation) | [Usage](#usage) | [Examples](#examples)

### Installation

From the command line:

```
git clone git@github.com:gphk-metrics/stata-multe
```

(or download the code manually and unzip). From Stata:

```
cap noi net uninstall multe
net install multe, from(`c(pwd)'/stata-multe)
```

(Change `stata-multe` if you download the package to a different
folder; e.g. `stata-multe-main`.) Note if the repo were public, this
could be installed directly from Stata:

```
local github "https://raw.githubusercontent.com"
net install manyiv, from(`github'/gphk-metrics/stata-multe/main/)
```

### Usage

```
multe depvar treatment, control(varname) [options]
help multe
```

### Examples

```stata
local nobs   1000
local ktreat 5
clear
set seed 1729
set obs `nobs'
gen T = ceil(runiform() * `ktreat')
gen W = mod(_n, 10)
gen Y = T + runiform()
multe Y T, control(W)
ereturn list
mata `e(mata)'.decomposition.print(0)
mata `e(mata)'.decomposition.print(1)
multe, vce(oracle)
multe Y T, control(W) gen(lambda tau)
multe Y T, control(W) gen(lambda(awesomeName) tau(coolerName))
desc, full
```
