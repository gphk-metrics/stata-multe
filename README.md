MulTE
=====

<!--

    README.md
    multe.pkg
    stata.toc
    doc/multe.sthlp
    src/ado/multe.ado
    src/mata/multe.mata

-->

Multiple Treatment Effects regression

`version 0.1.2 25Feb2022` | [Installation](#installation) | [Usage](#usage) | [Examples](#examples)

### Installation

From the command line:

```
git clone git@github.com:mcaceresb/stata-multe
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
net install manyiv, from(`github'/mcaceresb/stata-multe/main/)
```

### Usage

```
multe depvar treatment, control(varname) [options]
help multe
```

### Examples

```stata
clear
set seed 1729
set obs 1000
gen T = (runiform() > 0.5) + 1
gen W = mod(_n, 10) + 1
gen Y = T + runiform()
multe Y T, control(W)
return list
```
