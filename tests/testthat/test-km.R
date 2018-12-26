context("Kaplan-meier plot")

library(survival)
data(colon);data(pbc)

test_that("Run jskm", {
  fit <- survfit(Surv(time,status)~rx, data=colon)
  jskm(fit, timeby=500)
  jskm(fit, timeby=500, table = T)
  jskm(fit, timeby=500, main = "kaplan", xlabs = "time", ylabs = "suvrival")
  expect_is(jskm(fit, timeby=500), "gg")
  expect_is(jskm(fit), "gg")
  expect_is(jskm(fit, timeby=500, ci = T), "gg")
  expect_is(jskm(fit, timeby=500, legend = F), "gg")
  expect_is(jskm(fit, timeby=500, cumhaz = T), "gg")
})



test_that("Run svyjskm", {
  pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
  biasmodel<-glm(randomized~age*edema,data=pbc)
  pbc$randprob<-fitted(biasmodel)
  dpbc<-survey::svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))
  s1 <- survey::svykm(Surv(time,status>0)~sex, design=dpbc, se = T)
  expect_is(svyjskm(s1, ci = T), "gg")
  expect_is(svyjskm(s1, ci = F), "gg")
  expect_is(svyjskm(s1, ci = T, cumhaz = T), "gg")
  expect_is(svyjskm(s1, ci = F, cumhaz = T), "gg")
  s2 <- survey::svykm(Surv(time,status>0)~sex, design=dpbc, se = F)
  expect_is(svyjskm(s2, ci = F), "gg")
  expect_error(svyjskm(s2, ci = T))
  s3 <- survey::svykm(Surv(time,status>0)~1, design=dpbc, se = T)
  expect_is(svyjskm(s3, ci = T), "gg")
  expect_is(svyjskm(s3, ci = F), "gg")
  s4 <- survey::svykm(Surv(time,status>0)~1, design=dpbc, se = F)
  expect_is(svyjskm(s4, ci = F), "gg")
  expect_error(svyjskm(s4, ci = T))
  
})
