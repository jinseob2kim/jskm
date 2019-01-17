jskm
================

Kaplan-Meier Plot with 'ggplot2': 'survfit' and 'svykm' objects from 'survival' and 'survey' packages.

[![Build
Status](https://travis-ci.org/jinseob2kim/jskm.svg?branch=master)](https://travis-ci.org/jinseob2kim/jskm)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/jinseob2kim/jskm?branch=master&svg=true)](https://ci.appveyor.com/project/jinseob2kim/jskm)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/jskm)](https://cran.r-project.org/package=jskm)
[![codecov](https://codecov.io/github/jinseob2kim/jskm/branch/master/graphs/badge.svg)](https://codecov.io/github/jinseob2kim/jskm)
[![GitHub
issues](https://img.shields.io/github/issues/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/issues)
[![GitHub
forks](https://img.shields.io/github/forks/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/network)
[![GitHub
stars](https://img.shields.io/github/stars/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/stargazers)
[![GitHub
license](https://img.shields.io/github/license/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/blob/master/LICENSE)
[![GitHub last
commit](https://img.shields.io/github/last-commit/google/skia.svg)](https://github.com/jinseob2kim/jskm)
[![GitHub
contributors](https://img.shields.io/github/contributors/jinseob2kim/jskm.svg?maxAge=2592000)](https://github.com/jinseob2kim/jskm/graphs/contributors)

## Install

``` r
install.packages("devtools")
library(devtools)
install_github("jinseob2kim/jskm")
library(jskm)
```

## Example

### Survival probability

``` r
#Load dataset
library(survival)
data(colon)
fit <- survfit(Surv(time,status)~rx, data=colon)

#Plot the data
jskm(fit)
```

``` r
jskm(fit, table = T, pval = T, label.nrisk = "No. at risk", size.label.nrisk = 8, 
     xlabs = "Time(Day)", ylabs = "Survival", ystratalabs = c("Obs", "Lev", "Lev + 5FU"), ystrataname = "rx",
     marks = F, timeby = 365, xlims = c(0, 3000), ylims = c(0.25, 1))
```


### Cumulative hazard: 1- Survival probability

``` r
jskm(fit, ci = T, cumhaz = T,  mark = F, ylab = "Cumulative hazard (%)", surv.scale = "percent", pval =T, pval.size = 6, pval.coord = c(300, 0.7))
```


### Weighted Kaplan-Meier plot - `svykm.object` in **survey** package

``` r
library(survey)
data(pbc, package="survival")
pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
biasmodel <- glm(randomized~age*edema,data=pbc)
pbc$randprob <- fitted(biasmodel)

dpbc<-svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))

s1 <-svykm(Surv(time,status>0) ~ 1, design = dpbc)
s2 <-svykm(Surv(time,status>0) ~ sex, design = dpbc)

svyjskm(s1)
```

``` r
svyjskm(s2, pval = T, design = dpbc)
```

``` r
svyjskm(s2, cumhaz = T, ylab = "Cumulative (%)", surv.scale = "percent", pval = T, design = dpbc, pval.coord = c(300, 0.7)) 
```

If you want to get **confidence interval**, you should apply `se = T`
option to `svykm` object.

``` r
s3 <- svykm(Surv(time,status>0) ~ sex, design=dpbc, se = T)
svyjskm(s3)
```


``` r
svyjskm(s3, ci = F)
```

