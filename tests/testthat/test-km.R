context("Kaplan-meier plot")

library(survival)
data(colon);data(pbc)

test_that("Run jskm", {
  fit <- survfit(Surv(time,status)~rx, data=colon)
  jskm(fit, timeby=500)
  expect_is(jskm(fit, timeby=500), "gg")
})



test_that("Run svyjskm", {
  pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
  biasmodel<-glm(randomized~age*edema,data=pbc)
  pbc$randprob<-fitted(biasmodel)
  dpbc<-survey::svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))
  s1 <- survey::svykm(Surv(time,status>0)~sex, design=dpbc)
  expect_is(svyjskm(s1), "gg")
})
