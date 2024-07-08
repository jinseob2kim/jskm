#' @title Creates a Kaplan-Meier plot for survfit object.
#' @description Creates a Kaplan-Meier plot with at risk tables below for survfit object.
#' @param sfit a survfit object
#' @param table logical: Create a table graphic below the K-M plot, indicating at-risk numbers?
#' @param xlabs x-axis label
#' @param ylabs y-axis label
#' @param xlims numeric: list of min and max for x-axis. Default = c(0,max(sfit$time))
#' @param ylims numeric: list of min and max for y-axis. Default = c(0,1)
#' @param surv.scale 	scale transformation of survival curves. Allowed values are "default" or "percent".
#' @param ystratalabs character list. A list of names for each strata. Default = names(sfit$strata)
#' @param ystrataname The legend name. Default = "Strata"
#' @param timeby numeric: control the granularity along the time-axis; defaults to 7 time-points. Default = signif(max(sfit$time)/7, 1)
#' @param main plot title
#' @param pval logical: add the pvalue to the plot?
#' @param pval.size numeric value specifying the p-value text size. Default is 5.
#' @param pval.coord numeric vector, of length 2, specifying the x and y coordinates of the p-value. Default values are NULL
#' @param pval.testname logical: add '(Log-rank)' text to p-value. Default = F
#' @param marks logical: should censoring marks be added?
#' @param shape what shape should the censoring marks be, default is a vertical line
#' @param legend logical. should a legend be added to the plot?
#' @param legendposition numeric. x, y position of the legend if plotted. Default=c(0.85,0.8)
#' @param ci logical. Should confidence intervals be plotted. Default = FALSE
#' @param subs = NULL,
#' @param label.nrisk Numbers at risk label. Default = "Numbers at risk"
#' @param size.label.nrisk Font size of label.nrisk. Default = 10
#' @param linecols Character or Character vector. Colour brewer pallettes too colour lines. Default ="Set1", "black" for black with dashed line, character vector for the customization of line colors.
#' @param dashed logical. Should a variety of linetypes be used to identify lines. Default = FALSE
#' @param cumhaz Show cumulative incidence function, Default: F
#' @param cluster.option Cluster option for p value, Option: "None", "cluster", "frailty", Default: "None"
#' @param cluster.var Cluster variable
#' @param data select specific data - for reactive input, Default = NULL
#' @param cut.landmark cut-off for landmark analysis, Default = NULL
#' @param showpercent Shows the percentages on the right side.
#' @param status.cmprsk Status value when competing risk analysis, Default = 2nd level of status variable
#' @param linewidth Line witdh, Default = 0.75
#' @param theme Theme of the plot, Default = NULL, "nejm" for NEJMOA style, "jama" for JAMA style
#' @param nejm.infigure.ratiow Ratio of infigure width to total width, Default = 0.6
#' @param nejm.infigure.ratioh Ratio of infigure height to total height, Default = 0.5
#' @param nejm.infigure.ylim y-axis limit of infigure, Default = c(0,1)
#' @param ... PARAM_DESCRIPTION
#' @return Plot
#' @details DETAILS
#' @author Jinseob Kim, but heavily modified version of a script created by Michael Way.
#' \url{https://github.com/michaelway/ggkm/}
#' I have packaged this function, added functions to namespace and included a range of new parameters.
#' @examples
#' library(survival)
#' data(colon)
#' fit <- survfit(Surv(time, status) ~ rx, data = colon)
#' jskm(fit, timeby = 500)
#' @rdname jskm
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_step
#' @importFrom ggplot2 scale_linetype_manual
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_blank
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 scale_colour_brewer
#' @importFrom ggplot2 geom_ribbon
#' @importFrom grid unit
#' @importFrom ggpubr ggarrange
#' @importFrom stats pchisq time as.formula
#' @importFrom patchwork inset_element
#' @importFrom survival survfit survdiff coxph Surv cluster frailty
#' @export


jskm <- function(sfit,
                 table = FALSE,
                 xlabs = "Time-to-event",
                 ylabs = NULL,
                 xlims = c(0, max(sfit$time)),
                 ylims = c(0, 1),
                 surv.scale = c("default", "percent"),
                 ystratalabs = names(sfit$strata),
                 ystrataname = "Strata",
                 timeby = signif(max(sfit$time) / 7, 1),
                 main = "",
                 pval = FALSE,
                 pval.size = 5,
                 pval.coord = c(NULL, NULL),
                 pval.testname = F,
                 marks = TRUE,
                 shape = 3,
                 legend = TRUE,
                 legendposition = c(0.85, 0.8),
                 ci = FALSE,
                 subs = NULL,
                 label.nrisk = "Numbers at risk",
                 size.label.nrisk = 10,
                 linecols = "Set1",
                 dashed = FALSE,
                 cumhaz = F,
                 cluster.option = "None",
                 cluster.var = NULL,
                 data = NULL,
                 cut.landmark = NULL,
                 showpercent = F,
                 status.cmprsk = NULL,
                 linewidth = 0.75,
                 theme = NULL,
                 nejm.infigure.ratiow = 0.6,
                 nejm.infigure.ratioh = 0.5,
                 nejm.infigure.ylim = c(0, 1),
                 ...) {
  #################################
  # sorting the use of subsetting #
  #################################

  n.risk <- n.censor <- surv <- strata <- lower <- upper <- NULL

  times <- seq(0, max(sfit$time), by = timeby)
  if (!is.null(theme) && theme == "nejm") legendposition <- "right"
  if (is.null(subs)) {
    if (length(levels(summary(sfit)$strata)) == 0) {
      subs1 <- 1
      subs2 <- 1:length(summary(sfit, censored = T)$time)
      subs3 <- 1:length(summary(sfit, times = times, extend = TRUE)$time)
    } else {
      subs1 <- 1:length(levels(summary(sfit)$strata))
      subs2 <- 1:length(summary(sfit, censored = T)$strata)
      subs3 <- 1:length(summary(sfit, times = times, extend = TRUE)$strata)
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
    subs1 <- which(regexpr(ssvar, levels(summary(sfit)$strata), perl = T) != -1)
    subs2 <- which(regexpr(ssvar, summary(sfit, censored = T)$strata, perl = T) != -1)
    subs3 <- which(regexpr(ssvar, summary(sfit, times = times, extend = TRUE)$strata, perl = T) != -1)
  }

  if (!is.null(subs) | !is.null(sfit$states)) pval <- FALSE

  ##################################
  # data manipulation pre-plotting #
  ##################################

  if (is.null(ylabs)) {
    if (cumhaz | !is.null(sfit$states)) {
      ylabs <- "Cumulative incidence"
    } else {
      ylabs <- "Survival probability"
    }
  }


  if (length(levels(summary(sfit)$strata)) == 0) {
    # [subs1]
    if (is.null(ystratalabs)) ystratalabs <- as.character(sub("group=*", "", "All"))
  } else {
    # [subs1]
    if (is.null(ystratalabs)) ystratalabs <- as.character(sub("group=*", "", names(sfit$strata)))
  }

  if (is.null(ystrataname)) ystrataname <- "Strata"
  m <- max(nchar(ystratalabs))
  times <- seq(0, max(sfit$time), by = timeby)

  if (length(levels(summary(sfit)$strata)) == 0) {
    Factor <- factor(rep("All", length(subs2)))
  } else {
    Factor <- factor(summary(sfit, censored = T)$strata[subs2], levels = names(sfit$strata))
  }

  # Data to be used in the survival plot


  if (is.null(sfit$state)) { # no cmprsk
    df <- data.frame(
      time = sfit$time[subs2],
      n.risk = sfit$n.risk[subs2],
      n.event = sfit$n.event[subs2],
      n.censor = sfit$n.censor[subs2],
      surv = sfit$surv[subs2],
      strata = Factor,
      upper = sfit$upper[subs2],
      lower = sfit$lower[subs2]
    )
  } else { # cmprsk
    if (is.null(status.cmprsk)) {
      status.cmprsk <- sfit$states[2]
    }
    col.cmprsk <- which(sfit$state == status.cmprsk)
    df <- data.frame(
      time = sfit$time[subs2],
      n.risk = sfit$n.risk[, 1][subs2],
      n.event = sfit$n.event[, col.cmprsk][subs2],
      n.censor = sfit$n.censor[subs2],
      surv = sfit$pstate[, col.cmprsk][subs2],
      strata = Factor,
      upper = sfit$upper[, col.cmprsk][subs2],
      lower = sfit$lower[, col.cmprsk][subs2]
    )
  }

  form <- sfit$call$formula

  if (!is.null(cut.landmark)) {
    if (is.null(data)) {
      data <- tryCatch(eval(sfit$call$data), error = function(e) e)
      if ("error" %in% class(data)) {
        stop("Landmark analysis requires data object. please input 'data' option")
      }
    }

    var.time <- as.character(form[[2]][[2]])
    var.event <- as.character(form[[2]][[3]])
    if (length(var.event) > 1) {
      var.event <- setdiff(var.event, as.character(as.symbol(var.event)))
      var.event <- var.event[sapply(var.event, function(x) {
        "warning" %in% class(tryCatch(as.numeric(x), warning = function(w) w))
      })]
    }
    data1 <- data
    data1[[var.event]][data1[[var.time]] >= cut.landmark] <- 0
    data1[[var.time]][data1[[var.time]] >= cut.landmark] <- cut.landmark

    sfit1 <- survfit(as.formula(form), data1)
    sfit2 <- survfit(as.formula(form), data[data[[var.time]] >= cut.landmark, ])

    if (is.null(sfit$states)) {
      if (length(levels(Factor)) == 1) {
        df2 <- merge(subset(df, time >= cut.landmark)[, c("time", "n.risk", "n.event", "n.censor", "strata")],
          data.frame(time = sfit2$time, surv = sfit2$surv, strata = "All", upper = sfit2$upper, lower = sfit2$lower),
          by = c("time", "strata")
        )
      } else {
        df2 <- merge(subset(df, time >= cut.landmark)[, c("time", "n.risk", "n.event", "n.censor", "strata")],
          data.frame(time = sfit2$time, surv = sfit2$surv, strata = rep(names(sfit2$strata), sfit2$strata), upper = sfit2$upper, lower = sfit2$lower),
          by = c("time", "strata")
        )
      }

      df11 <- rbind(subset(df, time < cut.landmark), df2[, names(df)])
      df <- rbind(df11, data.frame(time = cut.landmark, n.risk = summary(sfit, times = cut.landmark)$n.risk[[1]], n.event = 0, n.censor = 0, surv = 1, strata = levels(df$strata), upper = 1, lower = 1))
    } else {
      if (is.null(status.cmprsk)) {
        status.cmprsk <- sfit$states[2]
      }
      col.cmprsk <- which(sfit$state == status.cmprsk)

      if (length(levels(Factor)) == 1) {
        df2 <- merge(subset(df, time >= cut.landmark)[, c("time", "n.risk", "n.event", "n.censor", "strata")],
          data.frame(time = sfit2$time, surv = sfit2$pstate[, col.cmprsk], strata = "All", upper = sfit2$upper[, col.cmprsk], lower = sfit2$lower[, col.cmprsk]),
          by = c("time", "strata")
        )
      } else {
        df2 <- merge(subset(df, time >= cut.landmark)[, c("time", "n.risk", "n.event", "n.censor", "strata")],
          data.frame(time = sfit2$time, surv = sfit2$pstate[, col.cmprsk], strata = rep(names(sfit2$strata), sfit2$strata), upper = sfit2$upper[, col.cmprsk], lower = sfit2$lower[, col.cmprsk]),
          by = c("time", "strata")
        )
      }
      df11 <- rbind(subset(df, time < cut.landmark), df2[, names(df)])
      df <- rbind(df11, data.frame(time = cut.landmark, n.risk = summary(sfit, times = cut.landmark)$n.risk[[1]], n.event = 0, n.censor = 0, surv = 0, strata = levels(df$strata), upper = 0, lower = 0))
    }
  }


  if (cumhaz & is.null(sfit$states)) {
    upper.new <- 1 - df$lower
    lower.new <- 1 - df$upper
    df$surv <- 1 - df$surv
    df$lower <- lower.new
    df$upper <- upper.new
  }

  # Final changes to data for survival plot
  levels(df$strata) <- ystratalabs
  zeros <- data.frame(
    time = 0, n.risk = NA, n.event = NA, n.censor = NA, surv = 1,
    strata = factor(ystratalabs, levels = levels(df$strata)),
    upper = 1, lower = 1
  )
  if (cumhaz | !is.null(sfit$states)) {
    zeros$surv <- 0
    zeros$lower <- 0
    zeros$upper <- 0
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
    p <- p +
      scale_linetype_manual(name = ystrataname, values = linetype) +
      geom_step(data = subset(df, time >= cut.landmark), linewidth = linewidth) + geom_step(data = subset(df, time < cut.landmark), linewidth = linewidth)
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

  # Add censoring marks to the line:
  if (marks == TRUE) {
    p <- p + geom_point(data = subset(df, n.censor >= 1), aes(x = time, y = surv, colour = strata), shape = shape)
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
  if (showpercent == T) {
    if (is.null(cut.landmark)) {
      y.percent <- summary(sfit, times = xlims[2], extend = T)$surv
      if (!is.null(sfit$states)) {
        y.percent <- summary(sfit, times = xlims[2], extend = T)$pstate[, col.cmprsk]
      }
      if (cumhaz == T & is.null(sfit$states)) y.percent <- 1 - y.percent
      p <- p + annotate(geom = "text", x = xlims[2], y = y.percent, label = paste0(round(100 * y.percent, 1), "%"), color = "black")
      if (!is.null(theme) && theme == "nejm") {
        p1 <- p1 + annotate(geom = "text", x = xlims[2], y = y.percent, label = paste0(round(100 * y.percent, 1), "%"), color = "black", size = nejm.infigure.ratiow * 5)
      }
    } else {
      y.percent1 <- summary(sfit, times = cut.landmark, extend = T)$surv
      y.percent2 <- summary(sfit2, times = xlims[2], extend = T)$surv
      if (!is.null(sfit$states)) {
        y.percent1 <- summary(sfit, times = cut.landmark, extend = T)$pstate[, col.cmprsk]
        y.percent2 <- summary(sfit2, times = xlims[2], extend = T)$pstate[, col.cmprsk]
      }
      if (cumhaz == T & is.null(sfit$states)) {
        y.percent1 <- 1 - y.percent1
        y.percent2 <- 1 - y.percent2
      }
      p <- p + annotate(geom = "text", x = cut.landmark, y = y.percent1, label = paste0(round(100 * y.percent1, 1), "%"), color = "black") +
        annotate(geom = "text", x = xlims[2], y = y.percent2, label = paste0(round(100 * y.percent2, 1), "%"), color = "black")
      if (!is.null(theme) && theme == "nejm") {
        p1 <- p1 + annotate(geom = "text", x = cut.landmark, y = y.percent1, label = paste0(round(100 * y.percent1, 1), "%"), color = "black", size = nejm.infigure.ratiow * 5) +
          annotate(geom = "text", x = xlims[2], y = y.percent2, label = paste0(round(100 * y.percent2, 1), "%"), color = "black", size = nejm.infigure.ratiow * 5)
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

  #####################
  # p-value placement #
  ##################### a

  if (length(levels(summary(sfit)$strata)) == 0) pval <- F
  # if(!is.null(cut.landmark)) pval <- F

  if (pval == TRUE) {
    if (is.null(data)) {
      data <- tryCatch(eval(sfit$call$data), error = function(e) e)
      if ("error" %in% class(data)) {
        stop("'pval' option requires data object. please input 'data' option")
      }
    }

    if (is.null(cut.landmark)) {
      sdiff <- survival::survdiff(as.formula(form), data = data)
      pvalue <- pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE)

      ## cluster option
      if (cluster.option == "cluster" & !is.null(cluster.var)) {
        form.old <- as.character(form)
        form.new <- paste(form.old[2], form.old[1], " + ", form.old[3], " + cluster(", cluster.var, ")", sep = "")
        sdiff <- survival::coxph(as.formula(form.new), data = data, model = T, robust = T)
        pvalue <- summary(sdiff)$robscore["pvalue"]
      } else if (cluster.option == "frailty" & !is.null(cluster.var)) {
        form.old <- as.character(form)
        form.new <- paste(form.old[2], form.old[1], " + ", form.old[3], " + frailty(", cluster.var, ")", sep = "")
        sdiff <- survival::coxph(as.formula(form.new), data = data, model = T)
        pvalue <- summary(sdiff)$logtest["pvalue"]
      }

      pvaltxt <- ifelse(pvalue < 0.001, "p < 0.001", paste("p =", round(pvalue, 3)))
      if (pval.testname) pvaltxt <- paste0(pvaltxt, " (Log-rank)")

      # MOVE P-VALUE LEGEND HERE BELOW [set x and y]
      if (is.null(pval.coord)) {
        p <- p + annotate("text", x = (as.integer(max(sfit$time) / 5)), y = 0.1 + ylims[1], label = pvaltxt, size = pval.size)
      } else {
        p <- p + annotate("text", x = pval.coord[1], y = pval.coord[2], label = pvaltxt, size = pval.size)
      }
    } else {
      sdiff1 <- survival::survdiff(as.formula(form), data1)
      sdiff2 <- survival::survdiff(as.formula(form), data[data[[var.time]] >= cut.landmark, ])
      pvalue <- sapply(list(sdiff1, sdiff2), function(x) {
        pchisq(x$chisq, length(x$n) - 1, lower.tail = FALSE)
      })

      ## cluster option
      if (cluster.option == "cluster" & !is.null(cluster.var)) {
        form.old <- as.character(form)
        form.new <- paste(form.old[2], form.old[1], " + ", form.old[3], sep = "")
        sdiff1 <- survival::coxph(as.formula(form.new), data = data1, model = T, cluster = get(cluster.var))
        sdiff2 <- survival::coxph(as.formula(form.new), data = data[data[[var.time]] >= cut.landmark, ], model = T, cluster = get(cluster.var))
        pvalue <- sapply(list(sdiff1, sdiff2), function(x) {
          summary(x)$robscore["pvalue"]
        })
      } else if (cluster.option == "frailty" & !is.null(cluster.var)) {
        form.old <- as.character(form)
        form.new <- paste(form.old[2], form.old[1], " + ", form.old[3], " + frailty(", cluster.var, ")", sep = "")
        sdiff1 <- survival::coxph(as.formula(form.new), data = data1, model = T)
        sdiff2 <- survival::coxph(as.formula(form.new), data = data[data[[var.time]] >= cut.landmark, ], model = T)
        pvalue <- sapply(list(sdiff1, sdiff2), function(x) {
          summary(x)$logtest["pvalue"]
        })
      }

      pvaltxt <- ifelse(pvalue < 0.001, "p < 0.001", paste("p =", round(pvalue, 3)))

      if (pval.testname) pvaltxt <- paste0(pvaltxt, " (Log-rank)")

      if (is.null(pval.coord)) {
        p <- p + annotate("text", x = c(as.integer(max(sfit$time) / 10), as.integer(max(sfit$time) / 10) + cut.landmark), y = 0.1 + ylims[1], label = pvaltxt, size = pval.size)
      } else {
        p <- p + annotate("text", x = c(pval.coord[1], pval.coord[1] + cut.landmark), y = pval.coord[2], label = pvaltxt, size = pval.size)
      }
    }
  }

  ###################################################
  # Create table graphic to include at-risk numbers #
  ###################################################

  n.risk <- NULL
  if (length(levels(summary(sfit)$strata)) == 0) {
    Factor <- factor(rep("All", length(subs3)))
  } else {
    Factor <- factor(summary(sfit, times = times, extend = TRUE)$strata[subs3])
  }

  if (table == TRUE) {
    risk.data <- data.frame(
      strata = Factor,
      time = summary(sfit, times = times, extend = TRUE)$time[subs3],
      n.risk = summary(sfit, times = times, extend = TRUE)$n.risk[subs3]
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
