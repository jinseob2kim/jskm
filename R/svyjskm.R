#' @title Creates a Weighted Kaplan-Meier plot - svykm.object in survey package
#' @description Creates a Weighted Kaplan-Meier plot - svykm.object in survey package
#' @param sfit a svykm object
#' @param xlabs x-axis label, Default: 'Time-to-event'
#' @param ylabs y-axis label.
#' @param xlims numeric: list of min and max for x-axis. Default: NULL
#' @param ylims numeric: list of min and max for y-axis. Default: c(0, 1)
#' @param surv.scale 	scale transformation of survival curves. Allowed values are "default" or "percent".
#' @param ystratalabs character list. A list of names for each strata. Default: NULL
#' @param ystrataname The legend name. Default: 'Strata'
#' @param timeby numeric: control the granularity along the time-axis; defaults to 7 time-points.
#' @param main plot title, Default: ''
#' @param pval logical: add the pvalue to the plot?, Default: FALSE
#' @param pval.size numeric value specifying the p-value text size. Default is 5.
#' @param pval.coord numeric vector, of length 2, specifying the x and y coordinates of the p-value. Default values are NULL
#' @param pval.testname logical: add '(Log-rank)' text to p-value. Default = F
#' @param legend logical. should a legend be added to the plot? Default: TRUE
#' @param ci logical. Should confidence intervals be plotted. Default = NULL
#' @param legendposition numeric. x, y position of the legend if plotted. Default: c(0.85, 0.8)
#' @param linecols Character or Character vector. Colour brewer pallettes too colour lines. Default ="Set1", "black" for black with dashed line, character vector for the customization of line colors.
#' @param dashed logical. Should a variety of linetypes be used to identify lines. Default: FALSE
#' @param cumhaz Show cumulaive incidence function, Default: F
#' @param design Data design for reactive design data , Default: NULL
#' @param subs = NULL,
#' @param table logical: Create a table graphic below the K-M plot, indicating at-risk numbers?
#' @param label.nrisk Numbers at risk label. Default = "Numbers at risk"
#' @param size.label.nrisk Font size of label.nrisk. Default = 10
#' @param cut.landmark cut-off for landmark analysis, Default = NULL
#' @param showpercent Shows the percentages on the right side.
#' @param linewidth Line witdh, Default = 0.75
#' @param theme Theme of the plot, Default = NULL, "nejm" for NEJMOA style, "jama" for JAMA style
#' @param nejm.infigure.ratiow Ratio of infigure width to total width, Default = 0.6
#' @param nejm.infigure.ratioh Ratio of infigure height to total height, Default = 0.5
#' @param nejm.infigure.ylim y-axis limit of infigure, Default = c(0,1)
#' @param ... PARAM_DESCRIPTION
#' @return plot
#' @details DETAILS
#' @examples
#' library(survey)
#' data(pbc, package = "survival")
#' pbc$randomized <- with(pbc, !is.na(trt) & trt > 0)
#' biasmodel <- glm(randomized ~ age * edema, data = pbc)
#' pbc$randprob <- fitted(biasmodel)
#' dpbc <- svydesign(id = ~1, prob = ~randprob, strata = ~edema, data = subset(pbc, randomized))
#' s1 <- svykm(Surv(time, status > 0) ~ sex, design = dpbc)
#' svyjskm(s1)
#' @rdname svyjskm
#' @import ggplot2
#' @importFrom stats formula
#' @importFrom survey svyranktest
#' @importFrom survival Surv
#' @importFrom ggpubr ggarrange
#' @importFrom patchwork inset_element
#' @export

svyjskm <- function(sfit,
                    theme = NULL,
                    xlabs = "Time-to-event",
                    ylabs = "Survival probability",
                    xlims = NULL,
                    ylims = c(0, 1),
                    ystratalabs = NULL,
                    ystrataname = NULL,
                    surv.scale = c("default", "percent"),
                    timeby = NULL,
                    main = "",
                    pval = FALSE,
                    pval.size = 5,
                    pval.coord = c(NULL, NULL),
                    pval.testname = F,
                    legend = TRUE,
                    legendposition = c(0.85, 0.8),
                    ci = NULL,
                    linecols = "Set1",
                    dashed = FALSE,
                    cumhaz = F,
                    design = NULL,
                    subs = NULL,
                    table = F,
                    label.nrisk = "Numbers at risk",
                    size.label.nrisk = 10,
                    cut.landmark = NULL,
                    showpercent = F,
                    linewidth = 0.75,
                    nejm.infigure.ratiow = 0.6,
                    nejm.infigure.ratioh = 0.5,
                    nejm.infigure.ylim = c(0, 1),
                    ...) {
  surv <- strata <- lower <- upper <- NULL

  if (!is.null(theme) && theme == "nejm") legendposition <- "right"
  if (is.null(timeby)) {
    if (inherits(sfit, "svykmlist")) {
      timeby <- signif(max(sapply(sfit, function(x) {
        max(x$time)
      })) / 7, 1)
    } else if (inherits(sfit, "svykm")) {
      timeby <- signif(max(sfit$time) / 7, 1)
    }
  }

  if (is.null(ci)) {
    if (inherits(sfit, "svykmlist")) {
      ci <- "varlog" %in% names(sfit[[1]])
    } else if (inherits(sfit, "svykm")) {
      ci <- "varlog" %in% names(sfit)
    }
  }


  if (ci & !is.null(cut.landmark)) {
    if (is.null(design)) {
      design <- tryCatch(get(as.character(attr(sfit, "call")$design)), error = function(e) e)
      if ("error" %in% class(design)) {
        stop("'pval' option requires design object. please input 'design' option")
      }
    }
    var.time <- as.character(formula(sfit)[[2]][[2]])
    var.event <- as.character(formula(sfit)[[2]][[3]])
    if (length(var.event) > 1) {
      var.event <- setdiff(var.event, as.character(as.symbol(var.event)))
      var.event <- var.event[sapply(var.event, function(x) {
        "warning" %in% class(tryCatch(as.numeric(x), warning = function(w) w))
      })]
    }
    design1 <- design
    design1$variables[[var.event]][design1$variables[[var.time]] >= cut.landmark] <- 0
    design1$variables[[var.time]][design1$variables[[var.time]] >= cut.landmark] <- cut.landmark

    sfit2 <- survey::svykm(formula(sfit), design = subset(design, get(var.time) >= cut.landmark), se = T)
  }




  if (inherits(sfit, "svykmlist")) {
    if (is.null(ystrataname)) ystrataname <- as.character(formula(sfit)[[3]])

    if (ci) {
      if ("varlog" %in% names(sfit[[1]])) {
        df <- do.call(rbind, lapply(names(sfit), function(x) {
          data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv, "lower" = pmax(0, exp(log(sfit[[x]]$surv) - 1.96 * sqrt(sfit[[x]]$varlog))), "upper" = pmin(1, exp(log(sfit[[x]]$surv) + 1.96 * sqrt(sfit[[x]]$varlog))))
        }))
        if (!is.null(cut.landmark)) {
          df2 <- do.call(rbind, lapply(names(sfit2), function(x) {
            data.frame("strata" = x, "time" = sfit2[[x]]$time, "surv" = sfit2[[x]]$surv, "lower" = pmax(0, exp(log(sfit2[[x]]$surv) - 1.96 * sqrt(sfit2[[x]]$varlog))), "upper" = pmin(1, exp(log(sfit2[[x]]$surv) + 1.96 * sqrt(sfit2[[x]]$varlog))))
          }))
          df <- rbind(df[df$time < cut.landmark, ], data.frame("strata" = unique(df$strata), "time" = cut.landmark, "surv" = 1, "lower" = 1, "upper" = 1), df2)
        }
      } else {
        stop("No CI information in svykmlist object. please run svykm with se = T option.")
      }
    } else {
      df <- do.call(rbind, lapply(names(sfit), function(x) {
        data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv)
      }))
      if (!is.null(cut.landmark)) {
        for (v in unique(df$strata)) {
          if (nrow(subset(df, time == cut.landmark & strata == v)) == 0) {
            df <- rbind(df, data.frame(strata = v, time = cut.landmark, surv = 1))
          } else {
            df[df$time == cut.landmark & df$strata == v, "surv"] <- 1
          }

          df[df$time > cut.landmark & df$strata == v, "surv"] <- df[df$time > cut.landmark & df$strata == v, "surv"] / min(df[df$time < cut.landmark & df$strata == v, "surv"])
        }
      }
    }

    df$strata <- factor(df$strata, levels = names(sfit))
    times <- seq(0, max(sapply(sfit, function(x) {
      max(x$time)
    })), by = timeby)
    if (is.null(ystratalabs)) {
      ystratalabs <- levels(df$strata)
    }
    if (is.null(xlims)) {
      xlims <- c(0, max(sapply(sfit, function(x) {
        max(x$time)
      })))
    }
  } else if (inherits(sfit, "svykm")) {
    if (is.null(ystrataname)) ystrataname <- "Strata"

    if (ci) {
      if ("varlog" %in% names(sfit)) {
        df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv, "lower" = pmax(0, exp(log(sfit$surv) - 1.96 * sqrt(sfit$varlog))), "upper" = pmax(0, exp(log(sfit$surv) + 1.96 * sqrt(sfit$varlog))))
        if (!is.null(cut.landmark)) {
          df2 <- data.frame("strata" = "All", "time" = sfit2$time, "surv" = sfit2$surv, "lower" = pmax(0, exp(log(sfit2$surv) - 1.96 * sqrt(sfit2$varlog))), "upper" = pmax(0, exp(log(sfit2$surv) + 1.96 * sqrt(sfit2$varlog))))
          df <- rbind(df[df$time < cut.landmark, ], data.frame("strata" = "All", "time" = cut.landmark, "surv" = 1, "lower" = 1, "upper" = 1), df2)
        }
      } else {
        stop("No CI information in svykm object. please run svykm with se = T option.")
      }
    } else {
      df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv)
      if (!is.null(cut.landmark)) {
        if (nrow(subset(df, time == cut.landmark)) == 0) {
          df <- rbind(df, data.frame(strata = "All", time = cut.landmark, surv = 1))
        } else {
          df[df$time == cut.landmark, "surv"] <- 1
        }

        df[df$time > cut.landmark, "surv"] <- df[df$time > cut.landmark, "surv"] / min(df[df$time < cut.landmark, "surv"])
      }
    }

    times <- seq(0, max(sfit$time), by = timeby)
    if (is.null(ystratalabs)) {
      ystratalabs <- "All"
    }
    if (is.null(xlims)) {
      xlims <- c(0, max(sfit$time))
    }
  }

  m <- max(nchar(ystratalabs))





  if (cumhaz) {
    df$surv <- 1 - df$surv
    if (ci) {
      upper.new <- 1 - df$lower
      lower.new <- 1 - df$upper
      df$lower <- lower.new
      df$upper <- upper.new
    }
  }

  # Final changes to data for survival plot
  levels(df$strata) <- ystratalabs
  zeros <- data.frame("strata" = factor(ystratalabs, levels = levels(df$strata)), "time" = 0, "surv" = 1)
  if (ci) {
    zeros$upper <- 1
    zeros$lower <- 1
  }

  if (cumhaz) {
    zeros$surv <- 0
    if (ci) {
      zeros$lower <- 0
      zeros$upper <- 0
    }
  }

  df <- rbind(zeros, df)
  d <- length(levels(df$strata))

  ###################################
  # specifying axis parameteres etc #
  ###################################

  if (dashed == TRUE | all(linecols == "black")) {
    linetype <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")
  } else {
    linetype <- c("solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")
  }

  # Scale transformation
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  surv.scale <- match.arg(surv.scale)
  scale_labels <- ggplot2::waiver()
  if (surv.scale == "percent") scale_labels <- scales::percent


  p <- ggplot2::ggplot(df, aes(x = time, y = surv, colour = strata, linetype = strata)) +
    ggtitle(main)
  linecols2 <- linecols
  if (all(linecols == "black")) {
    linecols <- "Set1"
    p <- ggplot2::ggplot(df, aes(x = time, y = surv, linetype = strata)) +
      ggtitle(main)
  }

  # Set up theme elements
  p <- p + theme_bw() +
    theme(
      axis.title.x = element_text(vjust = 0.7),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      legend.position = legendposition,
      legend.background = element_rect(fill = NULL),
      legend.key = element_rect(colour = NA),
      panel.border = element_blank(),
      plot.margin = unit(c(0, 1, .5, ifelse(m < 10, 1.5, 2.5)), "lines"),
      axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
      axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black")
    ) +
    scale_x_continuous(xlabs, breaks = times, limits = xlims) +
    scale_y_continuous(ylabs, limits = ylims, labels = scale_labels)

  if (!is.null(theme) && theme == "jama") {
    p <- p + theme(
      panel.grid.major.x = element_blank()
    )
  } else {
    p <- p + theme(
      panel.grid.major = element_blank()
    )
  }


  # Removes the legend:
  if (legend == FALSE) {
    p <- p + theme(legend.position = "none")
  }

  # Add lines too plot
  if (is.null(cut.landmark)) {
    p <- p + geom_step(linewidth = linewidth) +
      scale_linetype_manual(name = ystrataname, values = linetype)
  } else {
    p <- p + geom_step(data = subset(df, time < cut.landmark), linewidth = linewidth) + geom_step(data = subset(df, time >= cut.landmark), linewidth = linewidth) +
      scale_linetype_manual(name = ystrataname, values = linetype)
  }

  brewer.palette <- c(
    "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2",
    "Set1", "Set2", "Set3", "Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples",
    "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"
  )

  if (!is.null(theme) && theme == "jama") {
    col.pal <- c("#00AFBB", "#E7B800", "#FC4E07")
    col.pal <- rep(col.pal, ceiling(length(ystratalabs) / 3))
  } else if (all(linecols %in% brewer.palette)) {
    col.pal <- NULL
  } else {
    col.pal <- linecols
    col.pal <- rep(col.pal, ceiling(length(ystratalabs) / length(linecols)))
  }

  if (is.null(col.pal)) {
    p <- p + scale_colour_brewer(name = ystrataname, palette = linecols)
  } else {
    p <- p + scale_color_manual(name = ystrataname, values = col.pal)
  }

  # Add 95% CI to plot
  if (ci == TRUE) {
    if (all(linecols2 == "black")) {
      p <- p + geom_ribbon(data = df, aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA)
    } else if (is.null(col.pal)) {
      p <- p + geom_ribbon(data = df, aes(ymin = lower, ymax = upper, fill = strata), alpha = 0.25, colour = NA) + scale_fill_brewer(name = ystrataname, palette = linecols)
    } else {
      p <- p + geom_ribbon(data = df, aes(ymin = lower, ymax = upper, fill = strata), alpha = 0.25, colour = NA) + scale_fill_manual(name = ystrataname, values = col.pal)
    }
  }

  if (!is.null(cut.landmark)) {
    p <- p + geom_vline(xintercept = cut.landmark, lty = 2)
  }

  p1 <- p
  ## p-value
  if (inherits(sfit, "svykm")) pval <- FALSE
  # if(is.null(design)) pval <- FALSE
  if (showpercent == TRUE) {
    if (is.null(cut.landmark)) {
      y.percent <- df[df$time %in% tapply(df$time, df$strata, max), "surv"]
      p <- p + annotate(geom = "text", x = xlims[2], y = y.percent, label = paste0(round(100 * y.percent, 1), "%"), color = "black")
      if (!is.null(theme) && theme == "nejm") {
        p1 <- p1 + annotate(geom = "text", x = xlims[2], y = y.percent, label = paste0(round(100 * y.percent, 1), "%"), color = "black", size = nejm.infigure.ratiow * 5)
      }
    } else {
      df.cut <- df[df$time < cut.landmark, ]
      y.percent1 <- df.cut[df.cut$time %in% tapply(df.cut$time, df.cut$strata, max), "surv"]
      y.percent2 <- df[df$time %in% tapply(df$time, df$strata, max), "surv"]
      p <- p + annotate(geom = "text", x = cut.landmark, y = y.percent1, label = paste0(round(100 * y.percent1, 1), "%"), color = "black") +
        annotate(geom = "text", x = xlims[2], y = y.percent2, label = paste0(round(100 * y.percent2, 1), "%"), color = "black")
      if (!is.null(theme) && theme == "nejm") {
        p1 <- p1 + annotate(geom = "text", x = cut.landmark, y = y.percent1, label = paste0(round(100 * y.percent1, 1), "%"), color = "black", size = nejm.infigure.ratiow * 5) +
          annotate(geom = "text", x = xlims[2], y = y.percent2, label = paste0(round(100 * y.percent2, 1), "%"), color = "black", size = nejm.infigure.ratiow * 5)
      }
    }
  }
  if (pval) {
    if (is.null(design)) {
      design <- tryCatch(get(as.character(attr(sfit, "call")$design)), error = function(e) e)
      if ("error" %in% class(design)) {
        stop("'pval' option requires design object. please input 'design' option")
      }
    }

    if (is.null(cut.landmark)) {
      sdiff <- survey::svylogrank(formula(sfit), design = design)
      pvalue <- sdiff[[2]][2]

      pvaltxt <- ifelse(pvalue < 0.001, "p < 0.001", paste("p =", round(pvalue, 3)))
      if (pval.testname) pvaltxt <- paste0(pvaltxt, " (Log-rank)")

      # MOVE P-VALUE LEGEND HERE BELOW [set x and y]
      if (is.null(pval.coord)) {
        p <- p + annotate("text", x = as.integer(max(sapply(sfit, function(x) {
          max(x$time) / 10
        }))), y = 0.1 + ylims[1], label = pvaltxt, size = pval.size)
      } else {
        p <- p + annotate("text", x = pval.coord[1], y = pval.coord[2], label = pvaltxt, size = pval.size)
      }
    } else {
      if (is.null(design)) {
        design <- tryCatch(get(as.character(attr(sfit, "call")$design)), error = function(e) e)
        if ("error" %in% class(design)) {
          stop("'pval' option requires design object. please input 'design' option")
        }
      }
      var.time <- as.character(formula(sfit)[[2]][[2]])
      var.event <- as.character(formula(sfit)[[2]][[3]])
      if (length(var.event) > 1) {
        var.event <- setdiff(var.event, as.character(as.symbol(var.event)))
        var.event <- var.event[sapply(var.event, function(x) {
          "warning" %in% class(tryCatch(as.numeric(x), warning = function(w) w))
        })]
      }
      design1 <- design
      design1$variables[[var.event]][design1$variables[[var.time]] >= cut.landmark] <- 0
      design1$variables[[var.time]][design1$variables[[var.time]] >= cut.landmark] <- cut.landmark

      sdiff1 <- survey::svylogrank(formula(sfit), design = design1)
      sdiff2 <- survey::svylogrank(formula(sfit), design = subset(design, get(var.time) >= cut.landmark))
      pvalue <- sapply(list(sdiff1, sdiff2), function(x) {
        x[[2]][2]
      })

      pvaltxt <- ifelse(pvalue < 0.001, "p < 0.001", paste("p =", round(pvalue, 3)))
      if (pval.testname) pvaltxt <- paste0(pvaltxt, " (Log-rank)")

      if (is.null(pval.coord)) {
        p <- p + annotate("text", x = c(as.integer(max(sapply(sfit, function(x) {
          max(x$time) / 10
        }))), as.integer(max(sapply(sfit, function(x) {
          max(x$time) / 10
        }))) + cut.landmark), y = 0.1 + ylims[1], label = pvaltxt, size = pval.size)
      } else {
        p <- p + annotate("text", x = c(pval.coord[1], pval.coord[1] + cut.landmark), y = pval.coord[2], label = pvaltxt, size = pval.size)
      }
    }
  }





  ## Create a blank plot for place-holding
  blank.pic <- ggplot(df, aes(time, surv)) +
    geom_blank() +
    theme_void() + ## Remove gray color
    theme(
      axis.text.x = element_blank(), axis.text.y = element_blank(),
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(), panel.border = element_blank()
    )

  ###################################################
  # Create table graphic to include at-risk numbers #
  ###################################################

  n.risk <- NULL
  if (table == TRUE) {
    if (is.null(design)) {
      sfit2 <- survival::survfit(formula(sfit), data = get(as.character(attr(sfit, "call")$design))$variables)
    } else {
      sfit2 <- survival::survfit(formula(sfit), data = design$variables)
    }

    # times <- seq(0, max(sfit2$time), by = timeby)

    if (is.null(subs)) {
      if (length(levels(summary(sfit2)$strata)) == 0) {
        subs1 <- 1
        subs2 <- 1:length(summary(sfit2, censored = T)$time)
        subs3 <- 1:length(summary(sfit2, times = times, extend = TRUE)$time)
      } else {
        subs1 <- 1:length(levels(summary(sfit2)$strata))
        subs2 <- 1:length(summary(sfit2, censored = T)$strata)
        subs3 <- 1:length(summary(sfit2, times = times, extend = TRUE)$strata)
      }
    } else {
      for (i in 1:length(subs)) {
        if (i == 1) {
          ssvar <- paste("(?=.*\\b=", subs[i], sep = "")
        }
        if (i == length(subs)) {
          ssvar <- paste(ssvar, "\\b)(?=.*\\b=", subs[i], "\\b)", sep = "")
        }
        if (!i %in% c(1, length(subs))) {
          ssvar <- paste(ssvar, "\\b)(?=.*\\b=", subs[i], sep = "")
        }
        if (i == 1 & i == length(subs)) {
          ssvar <- paste("(?=.*\\b=", subs[i], "\\b)", sep = "")
        }
      }
      subs1 <- which(regexpr(ssvar, levels(summary(sfit2)$strata), perl = T) != -1)
      subs2 <- which(regexpr(ssvar, summary(sfit2, censored = T)$strata, perl = T) != -1)
      subs3 <- which(regexpr(ssvar, summary(sfit2, times = times, extend = TRUE)$strata, perl = T) != -1)
    }

    if (!is.null(subs)) pval <- FALSE



    if (length(levels(summary(sfit2)$strata)) == 0) {
      Factor <- factor(rep("All", length(subs3)))
    } else {
      Factor <- factor(summary(sfit2, times = times, extend = TRUE)$strata[subs3])
    }


    risk.data <- data.frame(
      strata = Factor,
      time = summary(sfit2, times = times, extend = TRUE)$time[subs3],
      n.risk = summary(sfit2, times = times, extend = TRUE)$n.risk[subs3]
    )


    risk.data$strata <- factor(risk.data$strata, levels = rev(levels(risk.data$strata)))

    data.table <- ggplot(risk.data, aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
      geom_text(size = 3.5) +
      theme_bw() +
      scale_y_discrete(
        breaks = as.character(levels(risk.data$strata)),
        labels = rev(ystratalabs)
      ) +
      scale_x_continuous(label.nrisk, limits = xlims) +
      theme(
        axis.title.x = element_text(size = size.label.nrisk, vjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_text(face = "bold", hjust = 1)
      )
    data.table <- data.table +
      theme(legend.position = "none") + xlab(NULL) + ylab(NULL)


    # ADJUST POSITION OF TABLE FOR AT RISK
    data.table <- data.table +
      theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 3.1, 4.3) - 0.38 * m), "lines"))
  }

  #######################
  # Plotting the graphs #
  #######################
  if (!is.null(theme) && theme == "nejm") {
    p2 <- p1 + coord_cartesian(ylim = nejm.infigure.ylim) + theme(
      legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
      axis.text = element_text(size = 10 * nejm.infigure.ratiow)
    )
    p <- p + patchwork::inset_element(p2, 1 - nejm.infigure.ratiow, 1 - nejm.infigure.ratioh, 1, 1, align_to = "panel")
  }

  if (table == TRUE) {
    ggpubr::ggarrange(p, blank.pic, data.table,
      nrow = 3,
      # align = "v",
      heights = c(2, .1, .25)
    )
  } else {
    p
  }
}
