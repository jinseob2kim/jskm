#' @title Creates a Kaplan-Meier plot for survfit object.
#' @description Creates a Kaplan-Meier plot with at risk tables below for survfit object.
#' @param sfit a survfit object
#' @param table logical: Create a table graphic below the K-M plot, indicating at-risk numbers?
#' @param table.censor logical: Add numbers of censored in table graphic
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
#' @param med should a median line be added to the plot? Default = F
#' @param legend logical. should a legend be added to the plot?
#' @param legendposition numeric. x, y position of the legend if plotted. Default=c(0.85,0.8)
#' @param ci logical. Should confidence intervals be plotted. Default = FALSE
#' @param subs = NULL,
#' @param label.nrisk Numbers at risk label. Default = "Numbers at risk"
#' @param size.label.nrisk Font size of label.nrisk. Default = 10
#' @param linecols Character or Character vector. Colour imported from ggsci. Default ="Set1", "black" for black with dashed line, character vector for the customization of line colors.
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
#' @param nejm.infigure.xlim x-axis limit of infigure, Default = NULL
#' @param surv.by breaks unit in y-axis, default = NULL(ggplot default)
#' @param nejm.surv.by breaks unit in y-axis in nejm figure, default = NULL(ggplot default)
#' @param hr logical: add the hazard ratio to the plot?
#' @param hr.size numeric value specifying the HR text size. Default is 5.
#' @param hr.coord numeric vector, of length 2, specifying the x and y coordinates of the p-value. Default values are NULL
#' @param hr.testname logical: add '(Log-rank)' text to p-value. Default = F
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
#' @importFrom ggplot2 geom_ribbon
#' @importFrom grid unit
#' @importFrom ggpubr ggarrange
#' @importFrom stats pchisq time as.formula
#' @importFrom patchwork inset_element
#' @importFrom survival survfit survdiff coxph Surv cluster frailty
#' @importFrom cmprsk cuminc crr
#' @importFrom ggsci scale_color_npg scale_fill_npg scale_color_aaas scale_fill_aaas scale_color_nejm scale_fill_nejm scale_color_lancet scale_fill_lancet scale_color_jama scale_fill_jama scale_color_jco scale_fill_jco scale_color_frontiers scale_fill_frontiers
#' @export


jskm <- function(sfit,
                 table = FALSE,
                 table.censor = FALSE,
                 xlabs = "Time-to-event",
                 ylabs = NULL,
                 xlims = c(0, max(sfit$time)),
                 ylims = c(0, 1),
                 surv.scale = c("default", "percent"),
                 ystratalabs = NULL,
                 ystrataname = "Strata",
                 timeby = signif(max(sfit$time) / 7, 1),
                 main = "",
                 pval = FALSE,
                 pval.size = 5,
                 pval.coord = c(NULL, NULL),
                 pval.testname = T,
                 marks = TRUE,
                 shape = 3,
                 med = FALSE,
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
                 nejm.infigure.xlim = NULL,
                 nejm.infigure.ylim = c(0, 1),
                 surv.by = NULL,
                 nejm.surv.by = NULL,
                 hr = FALSE,
                 hr.size = 5,
                 hr.coord = c(NULL, NULL),
                 hr.testname = F,
                 ...) {
  #################################
  # sorting the use of subsetting #
  #################################
  
  test_type <- n.risk <- n.censor <- surv <- strata <- lower <- upper <- NULL
  
  times <- seq(0, max(sfit$time), by = timeby)
  has_weights <- !is.null(sfit$call$weights)
  if (!is.null(theme) && theme == "nejm") legendposition <- legendposition
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
  
  if ((!is.null(subs) | !is.null(sfit$states)) & is.null(status.cmprsk)) pval <- FALSE
  
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
  
  if (!is.null(status.cmprsk)) {
    if (is.null(data)) {
      data <- tryCatch(eval(sfit$call$data), error = function(e) e)
      if ("error" %in% class(data)) {
        stop("Competing-risk analysis requires data object. please input 'data' option")
      }
    }
    if (length(levels(summary(sfit)$strata)) == 0) {
      # [subs1]
      if (is.null(ystratalabs)) ystratalabs <- as.character("All")
    } else {
      # [subs1]
      if (is.null(ystratalabs)) {
        ystratalabs <- as.character(names(sfit$strata))
        ystratalabs <- gsub("^group=*", "", ystratalabs)
      }
    }
  } else {
    if (length(levels(summary(sfit)$strata)) == 0) {
      # [subs1]
      if (is.null(ystratalabs)) {
        ystratalabs <- as.character("All")
      }
      nc <- length(summary(sfit)$table)
      L <- summary(sfit)$table[nc - 1][[1]]
      U <- summary(sfit)$table[nc][[1]]
      median_time <- summary(sfit)$table["median"][[1]]
      ystratalabs2 <- paste0(ystratalabs, " (median : ", median_time, ", ", sfit$conf.int * 100, "% CI : ", L, " - ", U, ")")
    } else {
      # [subs1]
      if (is.null(ystratalabs)) {
        ystratalabs <- as.character(names(sfit$strata))
        ystratalabs <- gsub("^group=*", "", ystratalabs)
      }
      ystratalabs2 <- NULL
      for (i in 1:length(levels(summary(sfit)$strata))) {
        nc <- ncol(summary(sfit)$table)
        L <- summary(sfit)$table[, nc - 1][[i]]
        U <- summary(sfit)$table[, nc][[i]]
        median_time <- summary(sfit)$table[, "median"][[i]]
        ystratalabs2 <- c(ystratalabs2, paste0(ystratalabs[[i]], " (median : ", median_time, ", ", sfit$conf.int * 100, "% CI : ", L, " - ", U, ")"))
      }
    }
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
  time_var <- all.vars(form[[2]])[1]
  event_var <- all.vars(form[[2]])[2]
  group_var <- all.vars(form)[3]
  
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
  #직접 surv.scale을 지정해주는게 맞지 않나?
  surv.scale <- match.arg(surv.scale)
  scale_labels <- ggplot2::waiver()
  if (surv.scale == "percent") scale_labels <- scales::percent
  
  p <- ggplot2::ggplot(df, aes(x = time, y = surv, colour = strata, linetype = strata)) +
    ggtitle(main)
  
  
  linecols2 <- linecols
  if (all(linecols == "black")) {
    # linecols <- "Set1"
    p <- ggplot2::ggplot(df, aes(x = time, y = surv, linetype = strata)) +
      ggtitle(main)
  }
  
  # Set up theme elements
  p <- p + theme_bw() +
    theme(
      axis.title.x = element_text(vjust = 0.7),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      legend.position = "inside",
      legend.position.inside = legendposition,
      legend.background = element_rect(fill = NULL),
      legend.key = element_rect(colour = NA),
      panel.border = element_blank(),
      # plot.margin = unit(c(0, 1, .5, ifelse(m < 10, 1.5, 2.5)), "lines"),
      axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
      axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black")
    ) +
    scale_x_continuous(xlabs, breaks = times, limits = xlims)
  
  if (!is.null(surv.by)) {
    p <- p + scale_y_continuous(ylabs, limits = ylims, labels = scale_labels, breaks = seq(ylims[1], ylims[2], by = surv.by))
  } else {
    p <- p + scale_y_continuous(ylabs, limits = ylims, labels = scale_labels)
  }
  
  
  
  
  
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
    p <- p + guides(colour = "none", linetype = "none")
  }
  
  # Add lines too plot
  if (is.null(cut.landmark)) {
    if (med == T & is.null(status.cmprsk)) {
      p <- p + geom_step(linewidth = linewidth) +
        scale_linetype_manual(name = ystrataname, values = linetype, labels = ystratalabs2)
    } else {
      p <- p + geom_step(linewidth = linewidth) +
        scale_linetype_manual(name = ystrataname, values = linetype, labels = ystratalabs)
    }
  } else {
    if (med == T & is.null(status.cmprsk)) {
      p <- p +
        scale_linetype_manual(name = ystrataname, values = linetype, labels = ystratalabs2) +
        geom_step(data = subset(df, time >= cut.landmark), linewidth = linewidth) + geom_step(data = subset(df, time < cut.landmark), linewidth = linewidth)
    } else {
      p <- p +
        scale_linetype_manual(name = ystrataname, values = linetype, labels = ystratalabs) +
        geom_step(data = subset(df, time >= cut.landmark), linewidth = linewidth) + geom_step(data = subset(df, time < cut.landmark), linewidth = linewidth)
    }
  }
  
  
  Set1 <- c("#E41A1CFF", "#377EB8FF", "#4DAF4AFF", "#984EA3FF", "#FF7F00FF", "#FFFF33FF", "#A65628FF", "#F781BFFF", "#999999FF")
  ggsci_palettes <- c("Set1", "npg", "aaas", "nejm", "lancet", "jama", "jco", "frontiers")
  
  if (!is.null(theme) && theme == "jama") {
    col.pal <- c("#00AFBB", "#E7B800", "#FC4E07")
    col.pal <- rep(col.pal, ceiling(length(ystratalabs) / 3))
  } else if (linecols[1] %in% ggsci_palettes) {
    col.pal <- NULL
  } else {
    col.pal <- linecols
    col.pal <- rep(col.pal, ceiling(length(ystratalabs) / length(linecols)))
  }
  
  
  if (is.null(col.pal)) {
    p <- p + switch(linecols[1],
                    "Set1" = scale_color_manual(name = ystrataname, values = Set1, labels = if(med == T & is.null(status.cmprsk)) ystratalabs2 else ystratalabs),
                    "npg" = scale_color_npg(name = ystrataname, labels = if(med == T & is.null(status.cmprsk)) ystratalabs2 else ystratalabs),
                    "aaas" = scale_color_aaas(name = ystrataname, labels = if(med == T & is.null(status.cmprsk)) ystratalabs2 else ystratalabs),
                    "nejm" = scale_color_nejm(name = ystrataname, labels = if(med == T & is.null(status.cmprsk)) ystratalabs2 else ystratalabs),
                    "lancet" = scale_color_lancet(name = ystrataname, labels = if(med == T & is.null(status.cmprsk)) ystratalabs2 else ystratalabs),
                    "jama" = scale_color_jama(name = ystrataname, labels = if(med == T & is.null(status.cmprsk)) ystratalabs2 else ystratalabs),
                    "jco" = scale_color_jco(name = ystrataname, labels = if(med == T & is.null(status.cmprsk)) ystratalabs2 else ystratalabs),
                    "frontiers" = scale_color_frontiers(name = ystrataname, labels = if(med == T & is.null(status.cmprsk)) ystratalabs2 else ystratalabs))
    } else {
      p <- p + scale_color_manual(name = ystrataname, values = col.pal, labels = if(med == T & is.null(status.cmprsk)) ystratalabs2 else ystratalabs)
      }
  
  # Add censoring marks to the line:
  if (marks == TRUE) {
    p <- p + geom_point(data = subset(df, n.censor >= 1), aes(x = time, y = surv, colour = strata), shape = shape)
  }
  
  # Add median value
  
  
  if (med == TRUE & is.null(cut.landmark) & is.null(status.cmprsk)) {
    if (length(levels(summary(sfit)$strata)) == 0) {
      median_time <- summary(sfit)$table["median"][[1]]
      if (!is.na(median_time)) {
        p <- p + annotate("segment", x = xlims[1], xend = median_time, y = 0.5, yend = 0.5, linewidth = 0.3, linetype = "dashed") +
          annotate("segment", x = median_time, xend = median_time, y = ylims[1], yend = 0.5, linewidth = 0.3, linetype = "dashed")
      }
    } else {
      for (i in 1:length(levels(summary(sfit)$strata))) {
        median_time <- summary(sfit)$table[, "median"][[i]]
        if (!is.na(median_time)) {
          p <- p +
            annotate("segment", x = xlims[1], xend = median_time, y = 0.5, yend = 0.5, linewidth = 0.3, linetype = "dashed") +
            annotate("segment", x = median_time, xend = median_time, y = ylims[1], yend = 0.5, linewidth = 0.3, linetype = "dashed")
        }
      }
    }
  }
  
  if (med == TRUE & !is.null(cut.landmark) & is.null(status.cmprsk)) {
    if (length(levels(summary(sfit)$strata)) == 0) {
      median_time <- summary(sfit1)$table[, "median"][[1]]
      if (!is.na(median_time)) {
        p <- p + annotate("segment", x = xlims[1], xend = median_time, y = 0.5, yend = 0.5, linewidth = 0.3, linetype = "dashed") + annotate("segment", x = median_time, xend = median_time, y = ylims[1], yend = 0.5, linewidth = 0.3, linetype = "dashed")
      }
      median_time <- summary(sfit2)$table[, "median"][[1]]
      if (!is.na(median_time)) {
        p <- p + annotate("segment", x = xlims[1], xend = median_time, y = 0.5, yend = 0.5, linewidth = 0.3, linetype = "dashed") + annotate("segment", x = median_time, xend = median_time, y = ylims[1], yend = 0.5, linewidth = 0.3, linetype = "dashed")
      }
    } else {
      for (i in 1:length(levels(summary(sfit)$strata))) {
        median_time <- summary(sfit1)$table[, "median"][[i]]
        if (!is.na(median_time)) {
          p <- p + annotate("segment", x = xlims[1], xend = median_time, y = 0.5, yend = 0.5, linewidth = 0.3, linetype = "dashed") + annotate("segment", x = median_time, xend = median_time, y = ylims[1], yend = 0.5, linewidth = 0.3, linetype = "dashed")
        }
        median_time <- summary(sfit2)$table[, "median"][[i]]
        if (!is.na(median_time)) {
          p <- p + annotate("segment", x = xlims[1], xend = median_time, y = 0.5, yend = 0.5, linewidth = 0.3, linetype = "dashed") + annotate("segment", x = median_time, xend = median_time, y = ylims[1], yend = 0.5, linewidth = 0.3, linetype = "dashed")
        }
      }
    }
  }
  
  
  
  # Add 95% CI to plot
  if (ci == TRUE) {
    if (all(linecols2 == "black")) {
      p <- p + geom_ribbon(data = df, aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA)
    } else {
      p <- p + geom_ribbon(data = df, aes(ymin = lower, ymax = upper, fill = strata), alpha = 0.25, colour = NA)
      
      if (is.null(col.pal)) {
        p <- p + switch(linecols[1],
                        "Set1"      = scale_fill_manual(name = ystrataname, values = Set1, labels = if (isTRUE(med) && is.null(status.cmprsk) && (is.null(theme) || theme != "nejm")) ystratalabs2 else ystratalabs),
                        "npg"       = scale_fill_npg(name = ystrataname, labels = if (isTRUE(med) && is.null(status.cmprsk) && (is.null(theme) || theme != "nejm")) ystratalabs2 else ystratalabs),
                        "aaas"      = scale_fill_aaas(name = ystrataname, labels = if (isTRUE(med) && is.null(status.cmprsk) && (is.null(theme) || theme != "nejm")) ystratalabs2 else ystratalabs),
                        "nejm"      = scale_fill_nejm(name = ystrataname, labels = if (isTRUE(med) && is.null(status.cmprsk) && (is.null(theme) || theme != "nejm")) ystratalabs2 else ystratalabs),
                        "lancet"    = scale_fill_lancet(name = ystrataname, labels = if (isTRUE(med) && is.null(status.cmprsk) && (is.null(theme) || theme != "nejm")) ystratalabs2 else ystratalabs),
                        "jama"      = scale_fill_jama(name = ystrataname, labels = if (isTRUE(med) && is.null(status.cmprsk) && (is.null(theme) || theme != "nejm")) ystratalabs2 else ystratalabs),
                        "jco"       = scale_fill_jco(name = ystrataname, labels = if (isTRUE(med) && is.null(status.cmprsk) && (is.null(theme) || theme != "nejm")) ystratalabs2 else ystratalabs),
                        "frontiers" = scale_fill_frontiers(name = ystrataname, labels = if (isTRUE(med) && is.null(status.cmprsk) && (is.null(theme) || theme != "nejm")) ystratalabs2 else ystratalabs),
                        scale_fill_brewer(name = ystrataname, palette = linecols, labels = if (isTRUE(med) && is.null(status.cmprsk) && (is.null(theme) || theme != "nejm")) ystratalabs2 else ystratalabs)
        )
      } else {
        p <- p + scale_fill_manual(
          name   = ystrataname,
          values = col.pal,
          labels = if (isTRUE(med) && is.null(status.cmprsk) && (is.null(theme) || theme != "nejm")) ystratalabs2 else ystratalabs
        )
      }
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
  ##################### 
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
      if (!is.null(status.cmprsk)) {
        ci_obj <- cmprsk::cuminc(ftime = data[[time_var]], fstatus = data[[event_var]], group = data[[group_var]])
        pvalue <- ci_obj$Tests[, "pv"][1]
        test_type <- "Gray's Test"
      } else if (has_weights) {
        vv <- data[[group_var]]
        unique_groups <- unique(vv)
        n_groups <- length(unique_groups)
        if (n_groups != 2) {
          warning("P-value calculation is only available for binary group variables (2 groups). Number of groups found: ", n_groups)
          pval <- FALSE 
        } else {
          if (is.factor(vv) || is.character(vv)) {
            vv <- as.character(vv)
            unique_groups_sorted <- sort(unique_groups)
            vv <- ifelse(vv == unique_groups_sorted[1], 0, 1)
          } else if (is.numeric(vv)) {
            unique_values_sorted <- sort(unique_groups)
            vv <- ifelse(vv == unique_values_sorted[1], 0, 1)
          } else {
            warning("Unsupported group_var data type for p-value calculation.")
            pval <- FALSE
          }
          tt <- data[[time_var]]
          ff <- data[[event_var]]
          weight_var <- as.character(sfit$call$weights)
          weights <- data[[weight_var]]
          adj_lr_result <- adjusted.LR(tt, ff, vv, weights)
          pvalue <- adj_lr_result$p.value
          test_type <- "Adjusted Log-Rank Test"
        }
      } else {
        sdiff <- survival::survdiff(as.formula(form), data = data)
        pvalue <- pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE)
        test_type <- "Log-rank Test"
        ## cluster option
        if (cluster.option == "cluster" & !is.null(cluster.var)) {
          form.old <- as.character(form)
          form.new <- paste(form.old[2], form.old[1], " + ", form.old[3], " + cluster(", cluster.var, ")", sep = "")
          sdiff <- survival::coxph(as.formula(form.new), data = data, model = T, robust = T)
          pvalue <- summary(sdiff)$robscore["pvalue"]
          test_type <- "Cox (Cluster Robust)"
        } else if (cluster.option == "frailty" & !is.null(cluster.var)) {
          form.old <- as.character(form)
          form.new <- paste(form.old[2], form.old[1], " + ", form.old[3], " + frailty(", cluster.var, ")", sep = "")
          sdiff <- survival::coxph(as.formula(form.new), data = data, model = T)
          pvalue <- summary(sdiff)$logtest["pvalue"]
          test_type <- "Cox (Frailty)"
        }
      }
      pvaltxt <- ifelse(pvalue < 0.001, "p < 0.001", paste("p =", round(pvalue, 3)))
      if (pval.testname & !is.null(test_type)) {
        pvaltxt <- paste0(pvaltxt, " (", test_type, ")")
      }
      
      # MOVE P-VALUE LEGEND HERE BELOW [set x and y]
      if (is.null(pval.coord)) {
        p <- p + annotate("text", x = (as.integer(max(sfit$time) / 5)), y = 0.1 + ylims[1], label = pvaltxt, size = pval.size)
      } else {
        p <- p + annotate("text", x = pval.coord[1], y = pval.coord[2], label = pvaltxt, size = pval.size)
      }
    } else {
      if (!is.null(status.cmprsk)) {
        ci_obj1 <- cmprsk::cuminc(ftime = data1[[time_var]], fstatus = data1[[event_var]], group = data1[[group_var]])
        data2 <- data[data[[var.time]] >= cut.landmark, ]
        data2[[time_var]] <- data2[[time_var]] - cut.landmark
        ci_obj2 <- cmprsk::cuminc(ftime = data2[[time_var]], fstatus = data2[[event_var]], group = data2[[group_var]])
        pvalue1 <- ci_obj1$Tests[, "pv"][1]
        pvalue2 <- ci_obj2$Tests[, "pv"][1]
        pvalue <- c(pvalue1, pvalue2)
        test_type <- "Gray's Test"
      } else if (has_weights) {
        compute_pval_weighted <- function(sub_data, sfit, group_var, time_var, event_var) {
          vv_sub <- sub_data[[group_var]]
          unique_groups_sub <- unique(vv_sub)
          n_groups_sub <- length(unique_groups_sub)
          
          if (n_groups_sub != 2) {
            warning("P-value calculation is only available for binary group variables (2 groups) in landmark subset. Number of groups found: ", n_groups_sub)
            return(NA)
          } else {
            if (is.factor(vv_sub) || is.character(vv_sub)) {
              vv_sub <- as.character(vv_sub)
              unique_groups_sorted_sub <- sort(unique_groups_sub)
              vv_sub <- ifelse(vv_sub == unique_groups_sorted_sub[1], 0, 1)
            } else if (is.numeric(vv_sub)) {
              unique_values_sorted_sub <- sort(unique_groups_sub)
              vv_sub <- ifelse(vv_sub == unique_values_sorted_sub[1], 0, 1)
            } else {
              warning("Unsupported group_var data type for p-value calculation in landmark subset.")
              return(NA)
            }
            tt_sub <- sub_data[[time_var]]
            ff_sub <- sub_data[[event_var]]
            weight_var_sub <- as.character(sfit$call$weights)
            weights_sub <- sub_data[[weight_var_sub]]
            
            # Adjusted Log-Rank Test 
            adj_lr_result_sub <- adjusted.LR(tt_sub, ff_sub, vv_sub, weights_sub)
            return(adj_lr_result_sub$p.value)
          }
        }
        data2 <- data[data[[var.time]] >= cut.landmark, ]
        data2[[time_var]] <- data2[[time_var]] - cut.landmark
        pvalue_1 <- compute_pval_weighted(data1, sfit, group_var, time_var, event_var)
        pvalue_2 <- compute_pval_weighted(data2, sfit, group_var, time_var, event_var)
        pvalue <- c(pvalue_1, pvalue_2)
        test_type <- "Adjusted Log-Rank Test"
      } else {
        sdiff1 <- survival::survdiff(as.formula(form), data1)
        sdiff2 <- survival::survdiff(as.formula(form), data[data[[var.time]] >= cut.landmark, ])
        pvalue <- sapply(list(sdiff1, sdiff2), function(x) {
          pchisq(x$chisq, length(x$n) - 1, lower.tail = FALSE)
        })
        test_type <- "Log-rank Test"
        ## cluster option
        if (cluster.option == "cluster" & !is.null(cluster.var)) {
          form.old <- as.character(form)
          form.new <- paste(form.old[2], form.old[1], " + ", form.old[3], sep = "")
          sdiff1 <- survival::coxph(as.formula(form.new), data = data1, model = T, cluster = get(cluster.var))
          sdiff2 <- survival::coxph(as.formula(form.new), data = data[data[[var.time]] >= cut.landmark, ], model = T, cluster = get(cluster.var))
          pvalue <- sapply(list(sdiff1, sdiff2), function(x) {
            summary(x)$robscore["pvalue"]
          })
          test_type <- "Cox (Cluster Robust)"
        } else if (cluster.option == "frailty" & !is.null(cluster.var)) {
          form.old <- as.character(form)
          form.new <- paste(form.old[2], form.old[1], " + ", form.old[3], " + frailty(", cluster.var, ")", sep = "")
          sdiff1 <- survival::coxph(as.formula(form.new), data = data1, model = T)
          sdiff2 <- survival::coxph(as.formula(form.new), data = data[data[[var.time]] >= cut.landmark, ], model = T)
          pvalue <- sapply(list(sdiff1, sdiff2), function(x) {
            summary(x)$logtest["pvalue"]
          })
          test_type <- "Cox (Frailty)"
        }
      }
      pvaltxt <- ifelse(pvalue < 0.001, "p < 0.001", paste("p =", round(pvalue, 3)))
      
      if (pval.testname & !is.null(test_type)) {
        pvaltxt <- paste0(pvaltxt, " (", test_type, ")")
      }
      
      if (is.null(pval.coord)) {
        p <- p + annotate("text", x = c(as.integer(max(sfit$time) / 10), as.integer(max(sfit$time) / 10) + cut.landmark), y = 0.1 + ylims[1], label = pvaltxt, size = pval.size)
      } else {
        p <- p + annotate("text", x = c(pval.coord[1], pval.coord[1] + cut.landmark), y = pval.coord[2], label = pvaltxt, size = pval.size)
      }
    }
  }
  
  ##########################
  # Hazard Ratio placement #
  ##########################
  
  if (hr == TRUE) {
    if (is.null(data)) {
      data <- tryCatch(eval(sfit$call$data), error = function(e) e)
      if ("error" %in% class(data)) {
        stop("'HR' option requires data object. Please input 'data' option")
      }
    }
    
    # binary check.
    group_values <- data[[group_var]]
    unique_groups <- unique(group_values)
    n_groups <- length(unique_groups)
    if (n_groups != 2) {
      stop("Currently, HR calculation is only available for binary group variables. Number of groups found: ", n_groups)
    }
    
    # w/o Landmark
    if (is.null(cut.landmark)) {
      
      # 1) competing risk: Fine-Gray 
      if (!is.null(status.cmprsk)) {
        fg_model <- cmprsk::crr(ftime = data[[time_var]],
                                fstatus = data[[event_var]],
                                cov1 = as.matrix(data[[group_var]]))
        HR_value    <- exp(fg_model$coef[1])
        HR_ci_lower <- exp(fg_model$coef[1] - 1.96 * sqrt(fg_model$var[1,1]))
        HR_ci_upper <- exp(fg_model$coef[1] + 1.96 * sqrt(fg_model$var[1,1]))
        test_type <- "Fine-Gray Model"
        pval <- 2 * (1 - pnorm(abs(fg_model$coef[1] / sqrt(fg_model$var[1,1]))))
        
        
        # 2) weights: Weighted cox
      } else if (has_weights){
        weight_var <- as.character(sfit$call$weights)
        cox_model <- survival::coxph(as.formula(form), data = data, weights = data[[weight_var]])
        HR_value    <- summary(cox_model)$coefficients[,"exp(coef)"][1]
        HR_ci_lower <- summary(cox_model)$conf.int[1, "lower .95"]
        HR_ci_upper <- summary(cox_model)$conf.int[1, "upper .95"]
        pval <- summary(cox_model)$coefficients[, "Pr(>|z|)"][1]
        test_type <- "Weighted Cox Model"
        
        # 3) else, Cox__ w/ cluster(3-1), w/o(3-2)
      } else {
        cox_model <- survival::coxph(as.formula(form), data = data)
        HR_value    <- summary(cox_model)$coefficients[,"exp(coef)"][1]
        HR_ci_lower <- summary(cox_model)$conf.int[1, "lower .95"]
        HR_ci_upper <- summary(cox_model)$conf.int[1, "upper .95"]
        pval <- summary(cox_model)$coefficients[, "Pr(>|z|)"][1]
        test_type <- "Cox Model"
        # w/ cluster(3-1)
        if (cluster.option == "cluster" & !is.null(cluster.var)) {
          form.old <- as.character(form)
          form.new <- paste(form.old[2], form.old[1], " + ", form.old[3],
                            " + cluster(", cluster.var, ")", sep = "")
          cox_model <- survival::coxph(as.formula(form.new), data = data,
                                       model = TRUE, robust = TRUE)
          HR_value    <- summary(cox_model)$coefficients[,"exp(coef)"][1]
          HR_ci_lower <- summary(cox_model)$conf.int[1, "lower .95"]
          HR_ci_upper <- summary(cox_model)$conf.int[1, "upper .95"]
          test_type <- "Cox (Cluster Robust)"
          pval <- summary(cox_model)$coefficients[, "Pr(>|z|)"][1]
          # w/o cluster (3-2)
        } else if (cluster.option == "frailty" & !is.null(cluster.var)) {
          form.old <- as.character(form)
          form.new <- paste(form.old[2], form.old[1], " + ", form.old[3],
                            " + frailty(", cluster.var, ")", sep = "")
          cox_model <- survival::coxph(as.formula(form.new), data = data, model = TRUE)
          HR_value    <- summary(cox_model)$coefficients[,"exp(coef)"][1]
          HR_ci_lower <- summary(cox_model)$conf.int[1, "lower .95"]
          HR_ci_upper <- summary(cox_model)$conf.int[1, "upper .95"]
          test_type <- "Cox (Frailty)"
          pval <- summary(cox_model)$coefficients[, "Pr(>|z|)"][1]
        }
      }
      # HR text
      hr_txt <- ifelse(HR_value < 0.001, "HR < 0.001", paste("HR =", round(HR_value, 2)))
      hr_txt <- paste0(hr_txt, " (95% CI: ", round(HR_ci_lower, 2), " ", round(HR_ci_upper, 2), "; P = ", round(pval, 3), ")")
      
      if ((hr.testname == T) & !is.null(test_type)) {
        hr_txt <- paste0(hr_txt, " (", test_type, ")")
      }
      
      if (is.null(hr.coord)) {
        p <- p + annotate("text", x = (as.integer(max(sfit$time) / 5)),
                          y = 0.2 + ylims[1], label = hr_txt, size = hr.size)
      } else {
        p <- p + annotate("text", x = hr.coord[1], y = hr.coord[2],
                          label = hr_txt, size = hr.size)
      }
      
      #  w Landmark(2 HRs)
    } else {
      data1 <- data[data[[var.time]] < cut.landmark, ]
      cox_model1 <- survival::coxph(as.formula(form), data = data1)
      HR1    <- summary(cox_model1)$coefficients[,"exp(coef)"][1]
      HR1_ci_lower <- summary(cox_model1)$conf.int[1, "lower .95"]
      HR1_ci_upper <- summary(cox_model1)$conf.int[1, "upper .95"]
      pval1 <- summary(cox_model1)$coefficients[, "Pr(>|z|)"][1]
      
      
      data2 <- data[data[[var.time]] >= cut.landmark, ]
      data2[[time_var]] <- data2[[time_var]] - cut.landmark
      cox_model2 <- survival::coxph(as.formula(form), data = data2)
      HR2    <- summary(cox_model2)$coefficients[,"exp(coef)"][1]
      HR2_ci_lower <- summary(cox_model2)$conf.int[1, "lower .95"]
      HR2_ci_upper <- summary(cox_model2)$conf.int[1, "upper .95"]
      pval2 <- summary(cox_model2)$coefficients[, "Pr(>|z|)"][1]
      
      test_type <- "Cox Model"
      hr_txt <- paste0("HR1 = ", round(HR1, 2), 
                       " (95% CI: ", round(HR1_ci_lower, 2)," ", round(HR1_ci_upper, 2), ");", 
                       " P = ", round(pval1, 3),
                       "\n",
                       "HR2 = ", round(HR2, 2), 
                       " (95% CI: ", round(HR2_ci_lower, 2)," ", round(HR2_ci_upper, 2), ");",
                       " P = ", round(pval2, 3)
      )
      if ((hr.testname==T) & !is.null(test_type)) {
        hr_txt <- paste0(hr_txt, "\n(", test_type, ")")
      }
      
      if (is.null(hr.coord)) {
        p <- p + annotate("text",
                          x = as.integer(max(sfit$time) / 10),
                          y = 0.2 + ylims[1], label = hr_txt, size = hr.size)
      } else {
        p <- p + annotate("text",
                          x = hr.coord[1],
                          y = hr.coord[2], label = hr_txt, size = hr.size)
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
    sfit_unweighted <- survfit(as.formula(form), data = data)
    summary_unweighted <- summary(sfit_unweighted, times = times, extend = TRUE)
    
    risk.data <- data.frame(
      strata = Factor,
      time = summary_unweighted$time[subs3],
      n.risk = summary_unweighted$n.risk[subs3]
    )
    if (table.censor) {
      risk.data <- data.frame(
        strata = Factor,
        time = summary_unweighted$time[subs3],
        n.risk = summary_unweighted$n.risk[subs3],
        n.censor = summary_unweighted$n.censor[subs3]
      )
      risk.data$n.risk <- paste0(risk.data$n.risk, " (", risk.data$n.censor, ")")
      risk.data$n.censor <- NULL
      
      if (identical(label.nrisk, "Numbers at risk")) {
        label.nrisk <- "Numbers at risk (number censored)"
      }
    }
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
      guides(colour = "none", linetype = "none") + xlab(NULL) + ylab(NULL)
    
    
    # ADJUST POSITION OF TABLE FOR AT RISK
    data.table <- data.table +
      theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 3.1, 4.3) - 0.38 * m), "lines"))
  }
  
  
  #######################
  # Plotting the graphs #
  #######################
  
  if (!is.null(theme) && theme == "nejm") {
    ## both are NULL
    if (is.null(nejm.infigure.xlim) && is.null(nejm.surv.by)) {
      p2 <- p1 +
        coord_cartesian(ylim = nejm.infigure.ylim) +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text    = element_text(size = 10 * nejm.infigure.ratiow)
        ) +
        guides(colour = "none", linetype = "none") +
        scale_y_continuous(
          limits = nejm.infigure.ylim,
          breaks = waiver(),
          labels = scale_labels
        )
      
      ## nejm.infigure.xlim: NOT NULL, nejm.surv.by: NULL
    } else if (!is.null(nejm.infigure.xlim) && is.null(nejm.surv.by)) {
      p2 <- p1 +
        coord_cartesian(
          xlim = nejm.infigure.xlim,
          ylim = nejm.infigure.ylim
        ) +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text    = element_text(size = 10 * nejm.infigure.ratiow)
        ) +
        guides(colour = "none", linetype = "none") +
        scale_x_continuous(
          limits = nejm.infigure.xlim,
          breaks = signif(
            seq(
              nejm.infigure.xlim[1],
              nejm.infigure.xlim[2],
              length.out = 7
            ),
            2
          ),
          labels = waiver()
        ) +
        scale_y_continuous(
          limits = nejm.infigure.ylim,
          breaks = waiver(),
          labels = scale_labels
        )
      
      ## nejm.infigure.xlim: NULL, nejm.surv.by: NOT NULL
    } else if (is.null(nejm.infigure.xlim) && !is.null(nejm.surv.by)) {
      p2 <- p1 +
        coord_cartesian(ylim = nejm.infigure.ylim) +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text    = element_text(size = 10 * nejm.infigure.ratiow)
        ) +
        guides(colour = "none", linetype = "none") +
        scale_y_continuous(
          limits = nejm.infigure.ylim,
          breaks = seq(
            nejm.infigure.ylim[1],
            nejm.infigure.ylim[2],
            by = nejm.surv.by
          ),
          labels = scale_labels
        )
      
      ## both of them are NOT NULL
    } else {
      p2 <- p1 +
        coord_cartesian(
          xlim = nejm.infigure.xlim,
          ylim = nejm.infigure.ylim
        ) +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text    = element_text(size = 10 * nejm.infigure.ratiow)
        ) +
        guides(colour = "none", linetype = "none") +
        scale_x_continuous(
          limits = nejm.infigure.xlim,
          breaks = signif(
            seq(
              nejm.infigure.xlim[1],
              nejm.infigure.xlim[2],
              length.out = 7
            ),
            2
          ),
          labels = waiver()
        ) +
        scale_y_continuous(
          limits = nejm.infigure.ylim,
          breaks = seq(
            nejm.infigure.ylim[1],
            nejm.infigure.ylim[2],
            by = nejm.surv.by
          ),
          labels = scale_labels
        )
    }
    
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