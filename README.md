# jskm
R package for kaplan meier plot: Modified ggkm


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




