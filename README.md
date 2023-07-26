jskm
================

Kaplan-Meier Plot with ‘ggplot2’: ‘survfit’ and ‘svykm’ objects from
‘survival’ and ‘survey’ packages.

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/jinseob2kim/jskm?branch=master&svg=true)](https://ci.appveyor.com/project/jinseob2kim/jskm)
[![Github
action](https://github.com/jinseob2kim/jskm/workflows/R-CMD-check/badge.svg)](https://github.com/jinseob2kim/jskm/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/jskm)](https://cran.r-project.org/package=jskm)
[![CRAN_Download_Badge](https://cranlogs.r-pkg.org/badges/jskm)](https://CRAN.R-project.org/package=jskm)
[![codecov](https://codecov.io/github/jinseob2kim/jskm/branch/master/graphs/badge.svg)](https://app.codecov.io/github/jinseob2kim/jskm)
[![GitHub
issues](https://img.shields.io/github/issues/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/issues)
[![GitHub
stars](https://img.shields.io/github/stars/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/stargazers)
[![GitHub
license](https://img.shields.io/github/license/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/blob/master/LICENSE)

## Install

``` r
install.packages("jskm")

## From github: latest version
install.packages("remotes")
remotes::install_github("jinseob2kim/jskm")

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

![](man/figures/README-unnamed-chunk-1-1.png)<!-- -->

``` r
jskm(fit, table = T, pval = T, label.nrisk = "No. at risk", size.label.nrisk = 8, 
     xlabs = "Time(Day)", ylabs = "Survival", ystratalabs = c("Obs", "Lev", "Lev + 5FU"), ystrataname = "rx",
     marks = F, timeby = 365, xlims = c(0, 3000), ylims = c(0.25, 1), showpercent = T)
```

![](man/figures/README-unnamed-chunk-1-2.png)<!-- -->

### Cumulative hazard: 1- Survival probability

``` r
jskm(fit, ci = T, cumhaz = T,  mark = F, ylab = "Cumulative incidence (%)", surv.scale = "percent", pval =T, pval.size = 6, pval.coord = c(300, 0.7))
```

![](man/figures/README-unnamed-chunk-2-1.png)<!-- -->

### Landmark analysis

``` r
jskm(fit, mark = F,  surv.scale = "percent", pval =T, table = T, cut.landmark = 500)
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->

``` r
jskm(fit, mark = F,  surv.scale = "percent", pval =T, table = T, cut.landmark = 500, showpercent = T)
```

![](man/figures/README-unnamed-chunk-3-2.png)<!-- -->

### Competing risk analysis

`status2` variable: 0 - censoring, 1 - event, 2 - competing risk

``` r
## Make competing risk variable, Not real
colon$status2 <- colon$status
colon$status2[1:400] <- 2
colon$status2 <- factor(colon$status2)
fit2 <- survfit(Surv(time,status2)~rx, data=colon)
jskm(fit2, mark = F, surv.scale = "percent", table = T, status.cmprsk = "1")
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->

``` r
jskm(fit2, mark = F, surv.scale = "percent", table = T, status.cmprsk = "1", showpercent = T, cut.landmark = 500)
```

![](man/figures/README-unnamed-chunk-4-2.png)<!-- -->

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

![](man/figures/README-unnamed-chunk-5-1.png)<!-- -->

``` r
svyjskm(s2, pval = T,  table = T, design = dpbc)
```

![](man/figures/README-unnamed-chunk-5-2.png)<!-- -->

``` r
svyjskm(s2, cumhaz = T, ylab = "Cumulative incidence (%)", surv.scale = "percent", pval = T, design = dpbc, pval.coord = c(300, 0.7), showpercent = T) 
```

![](man/figures/README-unnamed-chunk-5-3.png)<!-- -->

If you want to get **confidence interval**, you should apply `se = T`
option to `svykm` object.

``` r
s3 <- svykm(Surv(time,status>0) ~ sex, design=dpbc, se = T)
svyjskm(s3)
```

![](man/figures/README-unnamed-chunk-6-1.png)<!-- -->

``` r
svyjskm(s3, ci = F)
```

![](man/figures/README-unnamed-chunk-6-2.png)<!-- -->

``` r
svyjskm(s3, ci = F,  surv.scale = "percent", pval =T, table = T, cut.landmark = 1000, showpercent = T)
```

![](man/figures/README-unnamed-chunk-6-3.png)<!-- -->
