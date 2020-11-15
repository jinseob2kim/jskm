context("Kaplan-meier plot")

library(survival)
data(colon);data(pbc)

test_that("Run jskm", {
  fit <- survfit(Surv(time,status)~rx, data=colon)
  expect_is(jskm(fit, timeby=500, table = T, pval = T), "gtable")
  expect_is(jskm(fit, table = T, pval = T, label.nrisk = "No. at risk", timeby = 365, xlabs = "Time(Day)", ylabs = "Survival", marks = F, xlims = c(0, 3500), ylims = c(0.25, 1),
       ystratalabs = c("Obs", "Lev", "Lev + 5FU"), ystrataname = "rx"), "gtable")
  expect_is(jskm(fit, timeby=500), "gg")
  expect_is(jskm(fit, timeby=500, main = "kaplan", xlabs = "time", ylabs = "Suvrival (%)", surv.scale = "percent"), "gg")
  expect_is(jskm(fit, pval.size = 7, pval.coord = c(100, 0.2),  pval.testname = T, showpercent = T), "gg")
  expect_is(jskm(fit, timeby=500, ci = T), "gg")
  expect_is(jskm(fit, timeby=500, legend = F, showpercent = T), "gg")
  expect_is(jskm(fit, timeby=500, cumhaz = T, cut.landmark = 500, showpercent = T), "gg")
  
  expect_is(jskm(fit, cluster.option = "cluster", cluster.var = "id", pval = T), "gg")
  expect_warning(jskm(fit, cluster.option = "frailty", cluster.var = "id", pval = T))
  
  fit3 <- survfit(Surv(time,status)~rx + frailty(id), data=colon)
  expect_is(jskm(fit3), "gg")
  expect_is(jskm(fit3, pval = T), "gg")
})




test_that("Run svyjskm", {
  pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
  biasmodel<-glm(randomized~age*edema,data=pbc)
  pbc$randprob<-fitted(biasmodel)
  dpbc<-survey::svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))
  s1 <- survey::svykm(Surv(time,status>0)~sex, design=dpbc, se = T)
  expect_is(svyjskm(s1, ci = T), "gg")
  expect_is(svyjskm(s1, ci = F, ylabs = "Suvrival (%)", surv.scale = "percent"), "gg")
  expect_is(svyjskm(s1, table = T, pval = T, design = dpbc), "gtable")
  expect_is(svyjskm(s1, ci = T, cumhaz = T), "gg")
  expect_is(svyjskm(s1, ci = F, cumhaz = T), "gg")
  s2 <- survey::svykm(Surv(time,status>0)~sex, design=dpbc, se = F)
  expect_is(svyjskm(s2, ci = F), "gg")
  expect_is(svyjskm(s2, design = dpbc, pval = T, pval.size = 7, pval.coord = c(100, 0.2), pval.testname = T), "gg")
  #pv <- svyjskm(s2, pval = T, design = dpbc)
  #expect_is(pv, "gg")
  expect_error(svyjskm(s2, ci = T))
  s3 <- survey::svykm(Surv(time,status>0)~1, design=dpbc, se = T)
  expect_is(svyjskm(s3, ci = T), "gg")
  expect_is(svyjskm(s3, ci = F, showpercent = T), "gg")
  s4 <- survey::svykm(Surv(time,status>0)~1, design=dpbc, se = F)
  expect_is(svyjskm(s4, ci = F), "gg")
  expect_error(svyjskm(s4, ci = T, showpercent = T))
  
})
