# jskm 0.5.13

* Update: Add `nejm.surv.by` option to y-axis breaks in `jskm` and `svyjskm`. 

# jskm 0.5.12

* Update: Add censoring marks

# jskm 0.5.11

* Update: Add Hazard ratio display: **hr** option.


# jskm 0.5.10

* Update: Display p value for IPTW with adjusted log-rank test and delete "group =" in strata.

# jskm 0.5.9

* Update: Display p value for competing risk using Gray's test and also reorganize pval.testname option to display corresponding test names. 


# jskm 0.5.8

* Update: Add `table.censor` option to show censored number in `jskm`. 

# jskm 0.5.7

* Update: Add `surv.by` option to y-axis breaks in `jskm` and `svyjskm`. 

# jskm 0.5.6

* Update: Add median expression to the graph in `jskm` and `svyjskm`. 

# jskm 0.5.5

* Update: Align color of censoring marks with color of lines in `jskm`.

# jskm 0.5.5

* Update: Align color of censoring marks with color of lines in `jskm`.

# jskm 0.5.4

* Update: Add customization of line colors to `jskm` and `svyjskm`

# jskm 0.5.3

* Update: Add theme('JAMA','NEJM') to `jskm` and `svyjskm`

# jskm 0.5.2

* Update: Add `linewidth` option

# jskm 0.5.1

* Fix: line color problem when apply landmark analysis

# jskm 0.5

* Add competing risk analysis to `jskm`

# jskm 0.4.4

* Fix: align plot with table

# jskm 0.4.3

* Bugfix: `svyjskm` **ystratalab** problem with releveled factor variable.

# jskm 0.4.2

* Update: Add option to show percent(%)

# jskm 0.4.1

* Update: Align risktable with plot.

# jskm 0.4.0

* Add landmark analysis: **cut.landmark** option.

* Remove dependency with **plyr** package.

# jskm 0.3.9

* Bug fixes: show p value "< 0.001" when p < 0.001

# jskm 0.3.8

* Bug fixes: Incorrect p-value when applying cluster options.  

# jskm 0.3.7

* Match the color: Line and 95% CI.  

# jskm 0.3.6

* Add linecols option: **"black"** for black with dashed line.

# jskm 0.3.5

## Fix

* Apply appropriate p-value to `svyjskm`

# jskm 0.3.1

## Update

* The p-value is expressed as the value rounded to the 3rd decimal place.

* Add **pval.testname** option, p-value is expressed with **(Log-rank)** text if **pval.testname = T**.

# jskm 0.3.0

## New feature

* Add **Number at risk table** option to `svyjskm`. The number at risk is no-weighted values (same to `jskm`).  

## Bug fixes

* Appropriate default **pval.coord** in `svyjskm`.

# jskm 0.2.8

## Bug fixes

* Add namespace **survival::frailty, survival::cluster**.

# jskm 0.2.7

## Bug fixes

* Fix some spell for cran release.

## Update

* Update **travis-ci**.

* Add **appveyor** CI to test **window** environment. 


# jskm 0.2.6

## Bug fixes

* Change `||`, `&&` to `|`, `&` for **Debian** environment. 

* Set some variable's initial values to `NULL` for cran release.

# jskm 0.2.5

## Bug fixes

* Remove gray rectangle above the **Number at risk** table.

* Change **p-value** position according **ylims** option.

## Update

* Add **Numbers at risk** label option to `jskm`: `label.nrisk`, `size.label.nrisk`.

* Add **percent scale** option : `surv.scale`.

* Add **pvalue position & font size** option : `pval.size`, `pval.coord` 

# jskm 0.2.4

## Update

* Add **ci** option to `svyjskm`.

* Get **p-value** without **design** option in `svyjskm`

# jskm 0.2.2

* Apply **testthat**

# jskm 0.2.1

* Can run when reactive data

# jskm 0.2.0

## New function

* `svyjskm` : Weighted Kaplan-Meier plot - `svykm.object` in **survey** package

# jskm 0.1.0

* `frailty`, `cluster` options and their **p value**