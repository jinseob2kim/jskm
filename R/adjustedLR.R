# R/adjusted_LR.R

#' Adjusted Log-Rank Test
#'
#' Performs an adjusted log-rank test considering weights.
#'
#' @param times A numeric vector of survival times.
#' @param failures A binary vector indicating event occurrence.
#' @param variable A numeric binary variable (0 and 1) indicating group membership.
#' @param weights A numeric vector of weights.
#'
#' @return A list containing the test statistic and p-value.
#'
#' @importFrom survival survfit
#' @importFrom stats pnorm
#' @keywords internal



adjusted.LR <- function(times, failures, variable, weights = NULL) {
  if (sum(times < 0) > 0) {
    print("Error times must be positive")
  } else {
    if (sum(weights <= 0) > 0) {
      print("Error weights must be superior to 0")
    } else {
      if (sum(failures != 0 & failures != 1) > 0) {
        print("Error failures must be must be a vector of 0 or 1")
      } else {
        if (is.null(weights)) {
          .w <- rep(1, length(times))
        } else {
          .w <- weights
        }

        crosssum <- function(seq, point, value) {
          crosssum <- numeric(length(seq))
          for (i in 1:length(point)) {
            loc <- sum(seq <= point[i])
            crosssum[loc] <- crosssum[loc] + value[i]
          }
          crosssum
        }
        times <- round(times, 4)
        n.unit <- rep(1, length(times))
        d.time <- survival::survfit(Surv(times, n.unit) ~ 1)$time
        n.risk0 <- crosssum(d.time, times[variable == 0], n.unit[variable == 0])
        n.risk1 <- crosssum(d.time, times[variable == 1], n.unit[variable == 1])
        n.risk0 <- rev(cumsum(rev(n.risk0)))
        n.risk1 <- rev(cumsum(rev(n.risk1)))
        n.risk <- n.risk0 + n.risk1
        n.event <- crosssum(d.time, times[failures == 1], failures[failures == 1])
        mod.rate <- ifelse(n.risk == 1, 0, n.event * (n.risk - n.event) / (n.risk * (n.risk - 1)))
        mod.rate[is.na(mod.rate)] <- 0

        w.risk0 <- crosssum(d.time, times[variable == 0], .w[variable == 0])
        w.risk1 <- crosssum(d.time, times[variable == 1], .w[variable == 1])
        w.risk0 <- rev(cumsum(rev(w.risk0)))
        w.risk1 <- rev(cumsum(rev(w.risk1)))
        w.risk <- w.risk0 + w.risk1
        w.risk0[w.risk0 == 0] <- 0.0001
        w.risk1[w.risk1 == 0] <- 0.0001

        subset0 <- (failures == 1) * (variable == 0)
        w.event0 <- crosssum(d.time, times[subset0 == 1], .w[subset0 == 1])
        subset1 <- (failures == 1) * (variable == 1)
        w.event1 <- crosssum(d.time, times[subset1 == 1], .w[subset1 == 1])
        w.event0 <- w.event0 * n.risk0 / w.risk0
        w.event1 <- w.event1 * n.risk1 / w.risk1
        w.event <- w.event0 + w.event1

        ww.risk0 <- crosssum(d.time, times[variable == 0], (.w * .w)[variable == 0])
        ww.risk1 <- crosssum(d.time, times[variable == 1], (.w * .w)[variable == 1])
        ww.risk0 <- rev(cumsum(rev(ww.risk0)))
        ww.risk1 <- rev(cumsum(rev(ww.risk1)))
        ww.risk0 <- ww.risk0 * n.risk0 * n.risk0 / (w.risk0 * w.risk0)
        ww.risk1 <- ww.risk1 * n.risk1 * n.risk1 / (w.risk1 * w.risk1)

        log.mean <- sum(w.event1 - n.risk1 * w.event / n.risk)

        log.var <- sum(mod.rate * (n.risk0 * n.risk0 * ww.risk1 + n.risk1 * n.risk1 * ww.risk0) / (n.risk * n.risk))
        v.akm <- log.mean / sqrt(log.var)

        return(list(
          statistic = v.akm,
          p.value = 2 * (1 - stats::pnorm(abs(v.akm)))
        ))
      }
    }
  }
}
