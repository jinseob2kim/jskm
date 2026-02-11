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
#' @param pval.size numeric value specifying the p-value text size. Default is 4.
#' @param pval.coord numeric vector, of length 2, specifying the x and y coordinates of the p-value. Default values are NULL
#' @param pval.testname logical: add '(Log-rank)' text to p-value. Default = F
#' @param marks logical: should censoring marks be added?
#' @param hr logical: add the Hazard Ratio to the plot?, Default: FALSE
#' @param hr.size numeric value specifying the Hazard Ratio text size. Default is 2.
#' @param hr.coord numeric vector, of length 2, specifying the x and y coordinates of the Hazard Ratio. Default values are NULL
#' @param med should a median line be added to the plot? Default = F
#' @param med.decimal Select the decimal to round to. Default = 2
#' @param legend logical. should a legend be added to the plot?
#' @param legendposition numeric. x, y position of the legend if plotted. Default=c(0.85,0.8)
#' @param ci logical. Should confidence intervals be plotted. Default = NULL
#' @param linecols Character or Character vector. Colour imported from ggsci to colour lines. Default ="Set1", "black" for black with dashed line, character vector for the customization of line colors.
#' @param dashed logical. Should a variety of linetypes be used to identify lines. Default: FALSE
#' @param cumhaz Show cumulaive incidence function, Default: F
#' @param design Data design for reactive design data , Default: NULL
#' @param subs = NULL,
#' @param table logical: Create a table graphic below the K-M plot, indicating at-risk numbers?
#' @param table.censor logical: Add numbers of censored in table graphic
#' @param label.nrisk Numbers at risk label. Default = "Numbers at risk"
#' @param left.nrisk  Should risk label be at top left? Default = "False"
#' @param size.label.nrisk Font size of label.nrisk. Default = 10
#' @param cut.landmark cut-off for landmark analysis, Default = NULL
#' @param showpercent Shows the percentages on the right side.
#' @param linewidth Line witdh, Default = 0.75
#' @param theme Theme of the plot, Default = NULL, "nejm" for NEJMOA style, "jama" for JAMA style
#' @param nejm.infigure.ratiow Ratio of infigure width to total width, Default = 0.6
#' @param nejm.infigure.ratioh Ratio of infigure height to total height, Default = 0.5
#' @param nejm.infigure.ylim y-axis limit of infigure, Default = c(0,1)
#' @param nejm.infigure.xlim x-axis limit of infigure, Default = NULL
#' @param surv.by breaks unit in y-axis, default = NULL(ggplot default)
#' @param nejm.surv.by breaks unit in y-axis in nejm figure, default = NULL(ggplot default)
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
#' @importFrom survey svyranktest svycoxph
#' @importFrom survival Surv
#' @importFrom ggpubr ggarrange
#' @importFrom patchwork inset_element
#' @importFrom ggsci scale_color_npg scale_fill_npg scale_color_aaas scale_fill_aaas scale_color_nejm scale_fill_nejm scale_color_lancet scale_fill_lancet scale_color_jama scale_fill_jama scale_color_jco scale_fill_jco scale_color_frontiers scale_fill_frontiers
#' @export

svyjskm <- function(sfit,
                    theme = NULL,
                    xlabs = "Time-to-event",
                    ylabs = NULL,
                    xlims = NULL,
                    ylims = c(0, 1),
                    ystratalabs = NULL,
                    ystrataname = "Group",           # 'Strata'에서 'Group'으로 변경
                    surv.scale = c("percent", "default"),
                    timeby = NULL,
                    main = "",
                    pval = FALSE,
                    pval.size = 5,                  # 4에서 5로 상향 (jskm 기준)
                    pval.coord = c(NULL, NULL),
                    pval.testname = F,
                    marks = FALSE, 
                    hr = FALSE,
                    hr.size = 5,                    # 2에서 5로 상향
                    hr.coord = c(NULL, NULL),
                    med = FALSE,
                    med.decimal = 2,
                    legend = TRUE,
                    legendposition = c(0.85, 0.8),
                    ci = NULL,
                    linecols = "Set1",
                    dashed = FALSE,
                    cumhaz = F,
                    design = NULL,
                    subs = NULL,
                    table = TRUE,                   # F에서 TRUE로 변경 (표 기본 출력)
                    table.censor = F,
                    label.nrisk = "No. at Risk",    # jskm과 동일하게 통일
                    left.nrisk = TRUE,              # FALSE에서 TRUE로 변경 (왼쪽 정렬 기본값)
                    size.label.nrisk = 11,          # 10에서 11로 상향
                    cut.landmark = NULL,
                    showpercent = F,
                    linewidth = 1.2,                # 0.75에서 1.2로 상향 (굵은 선)
                    nejm.infigure.ratiow = 0.6,
                    nejm.infigure.ratioh = 0.5,
                    nejm.infigure.xlim = NULL,
                    nejm.infigure.ylim = c(0, 1),
                    surv.by = 0.1,                  # NULL에서 0.1로 변경
                    nejm.surv.by = NULL,
                    ...) {
  n.censor <- surv <- strata <- lower <- upper <- NULL

  if (!is.null(theme) && theme == "nejm") legendposition <- legendposition
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
  
  
  if (inherits(sfit, "svykmlist")) {
    temp <- "varlog" %in% names(sfit[[1]])
    if (temp & marks){
      warning("Cannot mark when se=T")
      marks <- F
    }
  } else if (inherits(sfit, "svykm")) {
    temp <- "varlog" %in% names(sfit)
    if (temp & marks){
      warning("Cannot mark when se=T")
      marks <- F
    }
  }
  
  surv.scale <- match.arg(surv.scale)
  
  # 2. 축 숫자 레이블 설정 (scales::percent 대신 직접 계산)
  scale_labels <- ggplot2::waiver()
  if (surv.scale == "percent") {
    # 0.25 -> 25 로 표시해주는 익명 함수 사용
    scale_labels <- function(x) x * 100
  }
  
  
  # 1. y축 제목 설정 로직
  if (is.null(ylabs)) {
    if (cumhaz | (!is.null(sfit$states) && inherits(sfit, "svykm"))) {
      ylabs <- "Cumulative incidence"
    } else {
      ylabs <- "Survival probability"
    }
    
    # surv.scale이 percent인 경우 제목 뒤에 (%) 자동 추가
    if (surv.scale == "percent") {
      ylabs <- paste0(ylabs, " (%)")
    }
  }

  
  ##여기서 n.censor 시작할꺼야
  #if (marks == T){
    if (is.null(design)){
      design <- tryCatch(get(as.character(attr(sfit, "call")$design)), error = function(e) e)
      if ("error" %in% class(design)) {
        warning("Design not found, can't mark it")
        marks <- F
      }
    }
  if(!is.null(design)){
    var.time <- as.character(formula(sfit)[[2]][[2]])
    var.event <- as.character(formula(sfit)[[2]][[3]])
    if (length(var.event) > 1) {
      var.event <- setdiff(var.event, as.character(as.symbol(var.event)))
      var.event <- var.event[sapply(var.event, function(x) {
        "warning" %in% class(tryCatch(as.numeric(x), warning = function(w) w))
      })]
    }
    design1 <- design
    form <- formula(sfit)
    
    fit2 <- survfit(form, data = design1$variables)
  }
    
    
  #}
  
  
  

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
  
  
  if (!is.numeric(med.decimal) || length(med.decimal) != 1 || med.decimal < 0 || med.decimal %% 1 != 0) {
    stop("med.decimal must be a non-negative integer (e.g., 0, 1, 2).")
  }

  if (inherits(sfit, "svykmlist")) {
    if (is.null(ystrataname)) ystrataname <- as.character(formula(sfit)[[3]])

    if (ci) {
      if ("varlog" %in% names(sfit[[1]])) {
        if (!is.null(cut.landmark)) {
          df <- do.call(rbind, lapply(names(sfit), function(x) {
            data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv, "lower" = pmax(0, exp(log(sfit[[x]]$surv) - 1.96 * sqrt(sfit[[x]]$varlog))), "upper" = pmin(1, exp(log(sfit[[x]]$surv) + 1.96 * sqrt(sfit[[x]]$varlog))))
          }))
          df2 <- do.call(rbind, lapply(names(sfit2), function(x) {
            data.frame("strata" = x, "time" = sfit2[[x]]$time, "surv" = sfit2[[x]]$surv, "lower" = pmax(0, exp(log(sfit2[[x]]$surv) - 1.96 * sqrt(sfit2[[x]]$varlog))), "upper" = pmin(1, exp(log(sfit2[[x]]$surv) + 1.96 * sqrt(sfit2[[x]]$varlog))))
          }))
          df <- rbind(df[df$time < cut.landmark, ], data.frame("strata" = unique(df$strata), "time" = cut.landmark, "surv" = 1, "lower" = 1, "upper" = 1), df2)
        } else {
          if (med == TRUE) {
            df <- do.call(rbind, lapply(names(sfit), function(x) {
              df2 <- data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv, "lower" = pmax(0, exp(log(sfit[[x]]$surv) - 1.96 * sqrt(sfit[[x]]$varlog))), "upper" = pmin(1, exp(log(sfit[[x]]$surv) + 1.96 * sqrt(sfit[[x]]$varlog))))
              valid_indices <- which(df2$surv <= 0.5)
              closest_index <- valid_indices[which.min(abs(df2$surv[valid_indices] - 0.5))]
              df2$med <- df2[closest_index, "time"]
              return(df2)
            }))
          } else {
            df <- do.call(rbind, lapply(names(sfit), function(x) {
              data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv, "lower" = pmax(0, exp(log(sfit[[x]]$surv) - 1.96 * sqrt(sfit[[x]]$varlog))), "upper" = pmin(1, exp(log(sfit[[x]]$surv) + 1.96 * sqrt(sfit[[x]]$varlog))))
            }))
          }
        }
      } else {
        stop("No CI information in svykmlist object. please run svykm with se = T option.")
      }
    } else {
      if (!is.null(cut.landmark)) {
        df <- do.call(rbind, lapply(names(sfit), function(x) {
          data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv)
        }))
        for (v in unique(df$strata)) {
          if (nrow(subset(df, time == cut.landmark & strata == v)) == 0) {
            df <- rbind(df, data.frame(strata = v, time = cut.landmark, surv = 1))
          } else {
            df[df$time == cut.landmark & df$strata == v, "surv"] <- 1
          }

          df[df$time > cut.landmark & df$strata == v, "surv"] <- df[df$time > cut.landmark & df$strata == v, "surv"] / min(df[df$time < cut.landmark & df$strata == v, "surv"])
        }
      } else {
        if (med == TRUE) {
          df <- do.call(rbind, lapply(names(sfit), function(x) {
            df2 <- data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv)
            valid_indices <- which(df2$surv <= 0.5)
            closest_index <- valid_indices[which.min(abs(df2$surv[valid_indices] - 0.5))]
            df2$med <- df2[closest_index, "time"]
            return(df2)
          }))
          # mat <- cbind(fit2$time, fit2$n.censor)
          # colnames(mat) <- c("time", "n.censor")
          # mat_df <- as.data.frame(mat)
          # df <- merge(x = df, y = mat_df, by = "time", all.x = TRUE)
        } else {
          df <- do.call(rbind, lapply(names(sfit), function(x) {
            data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv)
          }))
         # if (marks){
            # mat <- cbind(fit2$time, fit2$n.censor)
            # colnames(mat) <- c("time", "n.censor")
            # mat_df <- as.data.frame(mat)
            # df <- merge(x = df, y = mat_df, by = "time", all.x = TRUE)
            #print(df2)
          #}
        }
      }
    }

    df$strata <- factor(df$strata, levels = names(sfit))
    times <- seq(0, max(sapply(sfit, function(x) {
      max(x$time)
    })), by = timeby)
    
    if (is.null(ystratalabs)) {
      df3 <- df[-c(1, 2), ]
      ystratalabs <- levels(df$strata)
      if (med == TRUE) {
        ystratalabs2 <- NULL
        for (i in 1:length(names(sfit))) {
          median_time <- unique(df3[df3$strata == names(sfit)[[i]], "med"])
          ystratalabs2 <- c(ystratalabs2, paste0(ystratalabs[[i]], " (median : ", round(median_time, med.decimal), ")"))
        }
      }
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
        if (!is.null(cut.landmark)) {
          df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv, "lower" = pmax(0, exp(log(sfit$surv) - 1.96 * sqrt(sfit$varlog))), "upper" = pmax(0, exp(log(sfit$surv) + 1.96 * sqrt(sfit$varlog))))
          df2 <- data.frame("strata" = "All", "time" = sfit2$time, "surv" = sfit2$surv, "lower" = pmax(0, exp(log(sfit2$surv) - 1.96 * sqrt(sfit2$varlog))), "upper" = pmax(0, exp(log(sfit2$surv) + 1.96 * sqrt(sfit2$varlog))))
          df <- rbind(df[df$time < cut.landmark, ], data.frame("strata" = "All", "time" = cut.landmark, "surv" = 1, "lower" = 1, "upper" = 1), df2)
        } else {
          if (med == T) {
            df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv, "lower" = pmax(0, exp(log(sfit$surv) - 1.96 * sqrt(sfit$varlog))), "upper" = pmax(0, exp(log(sfit$surv) + 1.96 * sqrt(sfit$varlog))))
            valid_indices <- which(df$surv <= 0.5)
            closest_index <- valid_indices[which.min(abs(df$surv[valid_indices] - 0.5))]
            df$med <- df[closest_index, "time"]
          } else {
            df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv, "lower" = pmax(0, exp(log(sfit$surv) - 1.96 * sqrt(sfit$varlog))), "upper" = pmax(0, exp(log(sfit$surv) + 1.96 * sqrt(sfit$varlog))))
          }
        }
      } else {
        stop("No CI information in svykm object. please run svykm with se = T option.")
      }
    } else {
      if (!is.null(cut.landmark)) {
        df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv)
        if (nrow(subset(df, time == cut.landmark)) == 0) {
          df <- rbind(df, data.frame(strata = "All", time = cut.landmark, surv = 1))
        } else {
          df[df$time == cut.landmark, "surv"] <- 1
        }

        df[df$time > cut.landmark, "surv"] <- df[df$time > cut.landmark, "surv"] / min(df[df$time < cut.landmark, "surv"])
      } else {
        if (med == T) {
          df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv)
          valid_indices <- which(df$surv <= 0.5)
          closest_index <- valid_indices[which.min(abs(df$surv[valid_indices] - 0.5))]
          df$med <- df[closest_index, "time"]
        } else {
          df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv)
        }
      }
    }

    times <- seq(0, max(sfit$time), by = timeby)
    if (is.null(ystratalabs)) {
      ystratalabs <- "All"
      if (med == TRUE) {
        ystratalabs2 <- paste0(ystratalabs, " (median : ", round(unique(df$med), med.decimal), ")")
      } else {
        ystratalabs2 <- ystratalabs
      }    
    }
    if (is.null(xlims)) {
      xlims <- c(0, max(sfit$time))
    }
  }

  
  if (!is.null(design)){
    mat <- cbind(fit2$time, fit2$n.censor)
    colnames(mat) <- c("time", "n.censor")
    mat_df <- as.data.frame(mat)
    df <- merge(x = df, y = mat_df, by = "time", all.x = TRUE)
  }else{
    df$n.censor <- NA
  }
  
  # [수정] 변수명(rx= 등) 제거 및 0점 가려짐 방지 공백 추가
  if (!is.null(ystratalabs)) {
    ystratalabs <- gsub("^.*=", "", ystratalabs) # 모든 변수명 대응
    ystratalabs <- paste0(ystratalabs, "    ")    # 물리적 여백
  }
  if (exists("ystratalabs2") && !is.null(ystratalabs2)) {
    ystratalabs2 <- gsub("^.*=", "", ystratalabs2)
    ystratalabs2 <- paste0(ystratalabs2, "    ")
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
  zeros <- if (med == T & is.null(cut.landmark)) {
    data.frame("strata" = factor(ystratalabs, levels = levels(df$strata)), "time" = 0, "surv" = 1, "med" = 0.5, "n.censor"=NA)
  } else {
    data.frame("strata" = factor(ystratalabs, levels = levels(df$strata)), "time" = 0, "surv" = 1, "n.censor"=NA)
  }

  if (ci) {
    if (med == T & is.null(cut.landmark)) {
      zeros$med <- NULL
      zeros$upper <- 1
      zeros$lower <- 1
      zeros$med <- 0.5
    } else {
      zeros$upper <- 1
      zeros$lower <- 1
    }
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
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      # x축 제목 간격 확보 (t=15)
      axis.title.x = element_text(size = 11, colour = "black", face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 11, colour = "black", face = "bold"),
      axis.text.x = element_text(size = 10, colour = "black"),
      axis.text.y = element_text(size = 10, colour = "black"),
      # 눈금 선 길게 (0.25cm)
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      legend.position = "inside",
      legend.position.inside = legendposition,
      legend.background = element_rect(fill = NULL),
      legend.key = element_rect(colour = NA),
      panel.border = element_blank()
    ) +
    # [수정] x축 여백 변경
    scale_x_continuous(xlabs, breaks = times, limits = xlims, expand = expansion(mult = c(0.001, 0.001)))

  # y축 설정 부분
  if (!is.null(surv.by)) {
    p <- p + scale_y_continuous(
      name = ylabs,                  
      limits = ylims, 
      labels = scale_labels, 
      breaks = seq(ylims[1], ylims[2], by = surv.by),
      expand = expansion(mult = c(0.001, 0.001)) # 상하 여백 변경
    )
  } else {
    p <- p + scale_y_continuous(
      name = ylabs,                  
      limits = ylims, 
      labels = scale_labels,
      expand = expansion(mult = c(0.001, 0.001)) # 상하 여백 변경
    )
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
    if (med == T) {
      p <- p + geom_step(linewidth = linewidth) +
        scale_linetype_manual(name = ystrataname, values = linetype, labels = ystratalabs2)
    } else {
      p <- p + geom_step(linewidth = linewidth) +
        scale_linetype_manual(name = ystrataname, values = linetype, labels = ystratalabs)
    }
  } else {
    p <- p + geom_step(data = subset(df, time < cut.landmark), linewidth = linewidth) + geom_step(data = subset(df, time >= cut.landmark), linewidth = linewidth) +
      scale_linetype_manual(name = ystrataname, values = linetype, labels = ystratalabs)
  }

  if (marks == TRUE) {
    p <- p + geom_point(data = subset(df, n.censor >= 1), aes(x = time, y = surv, colour = strata), shape = 3)
  }
  
  
  
  # Add median value
  if (med == TRUE & is.null(cut.landmark)) {
    df3 <- df[-c(1, 2), ]
    # [수정 후] 길이 체크 추가
    if (inherits(sfit, "svykm")) {
      median_time <- unique(df3$med)
      # 길이가 있고(>0), NA가 아닐 때만 실행
      if (length(median_time) > 0 && !is.na(median_time)) { 
        p <- p + annotate("segment", x = xlims[1], xend = median_time, y = 0.5, yend = 0.5, linewidth = 0.3, linetype = "dashed") + 
          annotate("segment", x = median_time, xend = median_time, y = ylims[1], yend = 0.5, linewidth = 0.3, linetype = "dashed")
      }
    } else {
      for (i in 1:length(names(sfit))) {
        median_time <- unique(df3[df3$strata == names(sfit)[[i]], "med"])
        # 길이가 있고(>0), NA가 아닐 때만 실행
        if (length(median_time) > 0 && !is.na(median_time)) {
          p <- p + annotate("segment", x = xlims[1], xend = median_time, y = 0.5, yend = 0.5, linewidth = 0.3, linetype = "dashed") + 
            annotate("segment", x = median_time, xend = median_time, y = ylims[1], yend = 0.5, linewidth = 0.3, linetype = "dashed")
        }
      }
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
  #
  #  if (is.null(cut.landmark)) {
  #    if (is.null(col.pal)) {
  #      p <- p + scale_colour_brewer(name = ystrataname, palette = linecols, labels = ystratalabs2)
  #    } else {
  #      p <- p + scale_color_manual(name = ystrataname, values = col.pal, labels = ystratalabs2)
  #    }} else {if (is.null(col.pal)) {
  #      p <- p + scale_colour_brewer(name = ystrataname, palette = linecols)
  #    } else {
  #      p <- p + scale_color_manual(name = ystrataname, values = col.pal)
  #    } }


  if (is.null(col.pal)) {
    p <- p + switch(linecols[1],
                    "Set1"      = scale_color_manual(name = ystrataname, values = Set1, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                    "npg"       = scale_color_npg(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                    "aaas"      = scale_color_aaas(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                    "nejm"      = scale_color_nejm(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                    "lancet"    = scale_color_lancet(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                    "jama"      = scale_color_jama(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                    "jco"       = scale_color_jco(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                    "frontiers" = scale_color_frontiers(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                    scale_colour_brewer(name = ystrataname, palette = linecols, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs)
    )
  } else {
    p <- p + scale_color_manual(name   = ystrataname, values = col.pal, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs)
  }
  # Add 95% CI to plot


  if (ci == TRUE) {
    if (all(linecols2 == "black")) {
      p <- p + geom_ribbon(data = df, aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA)
    } else {
      p <- p + geom_ribbon(data = df, aes(ymin = lower, ymax = upper, fill = strata), alpha = 0.25, colour = NA)
      
      if (is.null(col.pal)) {
        p <- p + switch(linecols[1],
                        "Set1"      = scale_fill_manual(name = ystrataname, values = Set1, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                        "npg"       = scale_fill_npg(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                        "aaas"      = scale_fill_aaas(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                        "nejm"      = scale_fill_nejm(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                        "lancet"    = scale_fill_lancet(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                        "jama"      = scale_fill_jama(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                        "jco"       = scale_fill_jco(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                        "frontiers" = scale_fill_frontiers(name = ystrataname, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs),
                        scale_fill_brewer(name = ystrataname, palette = linecols, labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs)
        )
      } else {
        p <- p + scale_fill_manual(
          name   = ystrataname,
          values = col.pal,
          labels = if (isTRUE(med) && is.null(cut.landmark)) ystratalabs2 else ystratalabs
        )
      }
    }
  }

  if (!is.null(cut.landmark)) {
    p <- p + geom_vline(xintercept = cut.landmark, lty = 2)
  }

  p1 <- p
  
  ###############
  # p-value     # 
  ###############
  if (inherits(sfit, "svykm")) pval <- FALSE
  # if(is.null(design)) pval <- FALSE
  if (showpercent == TRUE) {
  # is.null(landmark)
      if (is.null(cut.landmark)) {
      y.percent <- df[df$time %in% tapply(df$time, df$strata, max), "surv"]
      p <- p + annotate(geom = "text", x = xlims[2], y = y.percent, label = paste0(round(100 * y.percent, 1), "%"), color = "black")
      if (!is.null(theme) && theme == "nejm") {
        p1 <- p1 + annotate(geom = "text", x = xlims[2], y = y.percent, label = paste0(round(100 * y.percent, 1), "%"), color = "black", size = nejm.infigure.ratiow * 5)
      }
    } else {
      # landmark yes.
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

  ##################
  # modify. HR. #
  ###############
  
  if (hr == TRUE) {
    # design check
    
    # if (is.null(design)) {
    #   design <- tryCatch(get(as.character(attr(sfit, "call")$design)), error = function(e) e)
    #   if ("error" %in% class(design)) {
    #     stop("'HR' option requires design object. please input 'design' option")
    #   }
    # }
      if (is.null(design)) {
      design <- tryCatch(get(as.character(attr(sfit, "call")$design)), error = function(e) e)
      if ("error" %in% class(design) || !inherits(design, "survey.design")) {
        stop("'HR' option requires a valid survey design object. Please input a valid `design`.")
      }
    }
    
    # independent variable #n check
    independent_var <- all.vars(formula(sfit))[-c(1, 2)]  # 첫 번째는 Surv(), 두 번째는 time, 이후가 독립변수
    if (length(independent_var) != 1) {
      stop("HR calculation requires exactly one independent variable. Found: ", paste(independent_var, collapse = ", "))
    }
    
    # binary check
    independent_var_name <- independent_var[1]
    independent_var_data <- design$variables[[independent_var_name]]
    
    if (is.factor(independent_var_data)) {
      n_levels <- length(levels(independent_var_data))
    } else {
      n_levels <- length(unique(independent_var_data))
    }
    
    if (n_levels > 2) {
      stop(paste0("HR calculation only supports binary independent variables. Found ", 
                  n_levels, " levels in variable '", independent_var_name, "'."))
    }
    
    # landmark check
    if (is.null(cut.landmark)) {
      # no landmark - cox 비례위험.
      cox_model <- survey::svycoxph(formula(sfit), design = design)
      cox_summary <- summary(cox_model)
      
        # HR CI
        hr_value <- round(cox_summary$conf.int[1, "exp(coef)"], 2)  
        hr_ci_lower <- round(cox_summary$conf.int[1, "lower .95"], 2)  
        hr_ci_upper <- round(cox_summary$conf.int[1, "upper .95"], 2) 
        hr_pval <- round(cox_summary$coefficients[1, "Pr(>|z|)"], 3)  
        
        # HR text
        hr_text <- paste0("HR = ", hr_value, " (95% CI: ", hr_ci_lower, " - ", hr_ci_upper, 
                          ifelse(hr_pval < 0.001, "; P < 0.001", paste("; P =", hr_pval)), ")")
        
        # HR placement
        if (is.null(hr.coord)) {
          hr_x <- max(sapply(sfit, function(x) { max(x$time) })) / 4
          hr_y <- 0.15 + ylims[1]
        } else {
          hr_x <- hr.coord[1]
          hr_y <- hr.coord[2]
        }
        
        # annotate 
        p <- p + annotate("text", x = hr_x, y = hr_y, label = hr_text, size = hr.size)
      }else{
      
      # yes landmark
        # # modify
        var.time <- as.character(attr(sfit, "formula")[[2]][[2]])
        var.event_all <- as.character(attr(sfit, "formula")[[2]][[3]])
         var.event <- setdiff(var.event_all, c(">", "0", "1"))
         var.event <- if (length(var.event) > 1) var.event[1] else var.event

         if (!var.event %in% names(design$variables)) {
          stop(paste("Error: Variable", var.event, "not found in design$variables"))
        }
        # 
        
      # design1 <- subset(design, get(attr(sfit, "formula")[[2]]) < cut.landmark)
      # design2 <- subset(design, get(attr(sfit, "formula")[[2]]) >= cut.landmark)
      # design2$variables[[as.character(attr(sfit, "formula")[[2]])]] <- 
      #   design2$variables[[as.character(attr(sfit, "formula")[[2]])]] - cut.landmark
      design1 <- subset(design, design$variables[[var.time]] < cut.landmark)
      design2 <- subset(design, design$variables[[var.time]] >= cut.landmark)
      design2$variables[[var.time]] <- design2$variables[[var.time]] - cut.landmark
      
      # before landmark_HR
        cox_model1 <- survey::svycoxph(formula(sfit), design = design1)
        cox_summary1 <- summary(cox_model1)
        
        hr1 <- round(cox_summary1$conf.int[1, "exp(coef)"], 2)
        hr1_ci_lower <- round(cox_summary1$conf.int[1, "lower .95"], 2)
        hr1_ci_upper <- round(cox_summary1$conf.int[1, "upper .95"], 2)
        hr1_pval <- round(cox_summary1$coefficients[1, "Pr(>|z|)"], 3)
      
      # after landmark_HR
        cox_model2 <- survey::svycoxph(formula(sfit), design = design2)
        cox_summary2 <- summary(cox_model2)
        
        hr2 <- round(cox_summary2$conf.int[1, "exp(coef)"], 2)
        hr2_ci_lower <- round(cox_summary2$conf.int[1, "lower .95"], 2)
        hr2_ci_upper <- round(cox_summary2$conf.int[1, "upper .95"], 2)
        hr2_pval <- round(cox_summary2$coefficients[1, "Pr(>|z|)"], 3)
      
      # HR text
      hr_text1 <- paste0("HR1 = ", hr1, " (95% CI: ", hr1_ci_lower, " - ", hr1_ci_upper, 
                         ifelse(hr1_pval < 0.001, "; p < 0.001", paste("; p =", hr1_pval)), ")")
      hr_text2 <- paste0("HR2 = ", hr2, " (95% CI: ", hr2_ci_lower, " - ", hr2_ci_upper, 
                         ifelse(hr2_pval < 0.001, "; p < 0.001", paste("; p =", hr2_pval)), ")")
      
      # HR placement
      if (is.null(pval.coord)) {
        hr_x1 <- max(sapply(sfit, function(x) max(x$time))) / 8  # HR1
        hr_x2 <- hr_x1 + cut.landmark  # HR2
        hr_y <- 0.15 + ylims[1] #same
      } else {
         hr_x1 <- hr.coord[1]
         hr_x2 <- hr_x1 + cut.landmark  # HR2
         hr_y <- rep(hr.coord[2], 2)
      }
      
      # plot 
      p <- p + annotate("text", x = hr_x1, y = hr_y, label = hr_text1, size = hr.size) +
        annotate("text", x = hr_x2, y = hr_y, label = hr_text2, size = hr.size)
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
  # [수정] No. at risk 라벨을 굵은(Bold) 11폰트로 정확하게 설정
  header.pic <- ggplot(df, aes(x = time, y = surv)) +
    geom_blank() +
    theme_void() +
    annotate("text", 
             x = -Inf, y = Inf, 
             label = label.nrisk, 
             hjust = 0, vjust = 1, 
             size = size.label.nrisk / .pt,  # [수정] .pt로 나누어야 정확한 포인트 크기가 됨
             fontface = "bold") +            # [수정] 굵게 설정
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 18, unit = "pt"))

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
    if (table.censor) {
      risk.data <- data.frame(
        strata = Factor,
        time = summary(sfit2, times = times, extend = TRUE)$time[subs3],
        n.risk = summary(sfit2, times = times, extend = TRUE)$n.risk[subs3],
        n.censor = summary(sfit2, times = times, extend = TRUE)$n.censor[subs3]
      )
      risk.data$n.risk <- paste0(risk.data$n.risk, " (", risk.data$n.censor, ")")
      risk.data$n.censor <- NULL
      
      if (identical(label.nrisk, "Numbers at risk")) {
        label.nrisk <- "Numbers at risk (number censored)"
      }
    }

    risk.data$strata <- factor(risk.data$strata, levels = rev(levels(risk.data$strata)))
    
    # 1. data.table 생성 및 테마 설정 (여기를 수정하세요)
    data.table <- ggplot(risk.data, aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
      geom_text(size = 3.5) +
      theme_bw() +
      scale_y_discrete(breaks = as.character(levels(risk.data$strata)), labels = rev(ystratalabs)) +
      coord_cartesian(clip = "off") + 
      theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(face = "plain", colour = "black", size = 10, hjust = 0, margin = margin(r = 10)),
        
        # [핵심 수정] 왼쪽/오른쪽 불필요한 텍스트 및 잔상 완벽 제거
        axis.title.y = element_blank(),         # 왼쪽 세로 글씨 제거
        axis.text.y.right = element_blank(),    # 오른쪽 잔상 제거
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_blank(),
        
        panel.background = element_blank(),
        plot.background = element_blank()
      )
    
    # [수정] 표 x축도 그래프와 똑같이 expand = c(0, 0) 적용해야 정렬이 맞음
    if (left.nrisk == TRUE) {
      data.table <- data.table + scale_x_continuous(NULL, limits = xlims, expand = c(0, 0))
    } else {
      data.table <- data.table + scale_x_continuous(label.nrisk, limits = xlims, expand = c(0, 0))
    }
  }

  #   # ADJUST POSITION OF TABLE FOR AT RISK
  #   data.table <- data.table +
  #     theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 3.1, 4.3) - 0.38 * m), "lines"))
  # }

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
    g_plot <- ggplot2::ggplotGrob(p)
    g_table <- ggplot2::ggplotGrob(data.table)
    
    # 왼쪽 레이블 너비 일치화 (핵심)
    plot_left <- g_plot$widths[g_plot$layout[g_plot$layout$name == "axis-l", ]$l]
    table_left <- g_table$widths[g_table$layout[g_table$layout$name == "axis-l", ]$l]
    max_left <- grid::unit.pmax(plot_left, table_left)
    
    g_plot$widths[g_plot$layout[g_plot$layout$name == "axis-l", ]$l] <- max_left
    g_table$widths[g_table$layout[g_table$layout$name == "axis-l", ]$l] <- max_left
    
    # 전체 너비 맞춤
    maxWidth <- grid::unit.pmax(g_plot$widths, g_table$widths)
    g_plot$widths <- maxWidth
    g_table$widths <- maxWidth
    
    # 정렬 배치
    if (left.nrisk == TRUE) {
      ggpubr::ggarrange(
        plotlist = list(g_plot, header.pic, g_table),
        nrow = 3,
        heights = c(2, 0.07, 0.5)  # [수정] 0.08 -> 0.07 (jskm과 통일)
      )
      
    } else {
      ggpubr::ggarrange(
        plotlist = list(g_plot, blank.pic, g_table),
        nrow = 3,
        heights = c(2, 0.05, 0.5)
      )
    }
  } else {
    p
  }
}