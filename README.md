# jskm
R package for kaplan meier plot: Modified ggkm

[![Travis-CI Build Status](https://travis-ci.org/jinseob2kim/jskm.svg?branch=master)](https://travis-ci.org/jinseob2kim/jskm)
[![GitHub issues](https://img.shields.io/github/issues/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/issues)
[![GitHub forks](https://img.shields.io/github/forks/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/network)
[![GitHub stars](https://img.shields.io/github/stars/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/stargazers)
[![GitHub license](https://img.shields.io/github/license/jinseob2kim/jskm.svg)](https://github.com/jinseob2kim/jskm/blob/master/LICENSE)
[![GitHub last commit](https://img.shields.io/github/last-commit/google/skia.svg)](https://github.com/jinseob2kim/jskm)
[![GitHub contributors](https://img.shields.io/github/contributors/jinseob2kim/jskm.svg?maxAge=2592000)](https://github.com/jinseob2kim/jskm/graphs/contributors)


## Install
```r
install.packages("devtools")
library(devtools)
install_github("jinseob2kim/jskm")
library(jskm)
```


## Example

### Survival probability
```r
#Load dataset
library(survival)
data(colon)
fit <- survfit(Surv(time,status)~rx, data=colon)

#Plot the data
jskm(fit)
```

### Cumulative hazard: 1- Survival probability
```r
jskm(fit, cumhaz = T, ylab = "Cumulative hazard (%)")
```




