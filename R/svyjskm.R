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
#' @param linecols Character. Colour brewer pallettes too colour lines. Default: 'Set1', "black" for black with dashed line.
#' @param dashed logical. Should a variety of linetypes be used to identify lines. Default: FALSE
#' @param cumhaz Show cumulaive hazard function, Default: F
#' @param design Data design for reactive design data , Default: NULL
#' @param subs = NULL,
#' @param table logical: Create a table graphic below the K-M plot, indicating at-risk numbers?
#' @param label.nrisk Numbers at risk label. Default = "Numbers at risk"
#' @param size.label.nrisk Font size of label.nrisk. Default = 10
#' @param ... PARAM_DESCRIPTION
#' @return plot
#' @details DETAILS
#' @examples 
#'  library(survey)
#'  data(pbc, package="survival")
#'  pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
#'  biasmodel <- glm(randomized~age*edema,data=pbc)
#'  pbc$randprob <- fitted(biasmodel)
#'  dpbc <- svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))
#'  s1 <- svykm(Surv(time,status>0)~sex, design=dpbc)
#'  svyjskm(s1)
#' @rdname svyjskm
#' @import ggplot2
#' @importFrom stats formula
#' @importFrom survey svyranktest
#' @importFrom survival Surv
#' @export 

svyjskm <- function(sfit,
                    xlabs = "Time-to-event",
                    ylabs = "Survival probability",
                    xlims = NULL,
                    ylims = c(0,1),
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
                    legendposition=c(0.85,0.8),
                    ci = NULL,
                    linecols="Set1",
                    dashed= FALSE,
                    cumhaz = F,
                    design = NULL,
                    subs = NULL,
                    table = F,
                    label.nrisk = "Numbers at risk",
                    size.label.nrisk = 10,
                    ...) {
  
  surv <- strata <- lower <- upper <- NULL
  
  
  if (is.null(timeby)){
    if (class(sfit) == "svykmlist"){
      timeby <- signif(max(sapply(sfit, function(x){max(x$time)}))/7, 1)
    } else if(class(sfit) == "svykm"){
      
      timeby <- signif(max(sfit$time)/7, 1)
      }
  }
  
  if (class(sfit) == "svykmlist"){
    if(is.null(ystrataname)) ystrataname <- as.character(formula(sfit)[[3]])
    if (is.null(ci)){
      ci <- "varlog" %in% names(sfit[[1]])
    }
    if (ci){
      if ("varlog" %in% names(sfit[[1]])){
        df <- Reduce(rbind, lapply(names(sfit), function(x){data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv, "lower" = pmax(0, exp(log(sfit[[x]]$surv) - 1.96 * sqrt(sfit[[x]]$varlog))), "upper" = pmin(1, exp(log(sfit[[x]]$surv) + 1.96 * sqrt(sfit[[x]]$varlog))))}))
      } else{
        stop("No CI information in svykmlist object. please run svykm with se = T option.")
      }
    } else{
      df <- Reduce(rbind, lapply(names(sfit), function(x){data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv)}))
    }
    
    df$strata <- as.factor(df$strata)
    times <- seq(0, max(sapply(sfit, function(x){max(x$time)})), by = timeby)
    if (is.null(ystratalabs)){
      ystratalabs <- names(sfit)
    }
    if (is.null(xlims)){
      xlims <- c(0,max(sapply(sfit, function(x){max(x$time)})))
    }
    
  } else if(class(sfit) == "svykm"){
    if(is.null(ystrataname)) ystrataname <- "Strata"
    if (is.null(ci)){
      ci <- "varlog" %in% names(sfit)
    }
    if (ci){
      if ("varlog" %in% names(sfit)){
        df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv,  "lower" = pmax(0, exp(log(sfit$surv) - 1.96 * sqrt(sfit$varlog))), "upper" = pmax(0, exp(log(sfit$surv) + 1.96 * sqrt(sfit$varlog))))
      } else{
        stop("No CI information in svykm object. please run svykm with se = T option.")
      }
    } else{
      df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv)
    }
    
    times <- seq(0, max(sfit$time), by = timeby)
    if (is.null(ystratalabs)){
      ystratalabs <- "All"
    }
    if (is.null(xlims)){
      xlims <- c(0,max(sfit$time))
    }
  }
  
  m <- max(nchar(ystratalabs))
  
  
  
  
  
  if (cumhaz){
    df$surv <- 1 - df$surv
    if (ci){
      df$lower <- 1 - df$upper
      df$upper <- 1 - df$lower
    }
    }
  
  #Final changes to data for survival plot
  levels(df$strata) <- ystratalabs
  zeros <- data.frame("strata" = factor(ystratalabs, levels=levels(df$strata)), "time" = 0, "surv" = 1)
  if (ci){
    zeros$upper <- 1
    zeros$lower <- 1
  }
  
  if (cumhaz){
    zeros$surv <- 0
    if (ci){
      zeros$lower <- 0
      zeros$upper <- 0
    }
  }
  
  df <- rbind(zeros, df)
  d <- length(levels(df$strata))
  
  ###################################
  # specifying axis parameteres etc #
  ###################################
  
  if(dashed == TRUE | linecols == "black"){
    linetype=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")
  } else {
    linetype=c("solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")
  }
  
  # Scale transformation
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  surv.scale <- match.arg(surv.scale)
  scale_labels <-  ggplot2::waiver()
  if (surv.scale == "percent") scale_labels <- scales::percent
  
  
  p <- ggplot2::ggplot( df, aes(x=time, y=surv, colour=strata, linetype=strata)) + ggtitle(main)
  if (linecols == "black"){
    linecols <- "Set1"
    p <- ggplot2::ggplot( df, aes(x=time, y=surv, linetype=strata)) + ggtitle(main)
  }
  
  #Set up theme elements
  p <- p + theme_bw() +
    theme(axis.title.x = element_text(vjust = 0.7),
          panel.grid.minor = element_blank(),
          axis.line = element_line(size =0.5, colour = "black"),
          legend.position = legendposition,
          legend.background = element_rect(fill = NULL),
          legend.key = element_rect(colour = NA),
          panel.border = element_blank(),
          plot.margin = unit(c(0, 1, .5,ifelse(m < 10, 1.5, 2.5)),"lines"),
          panel.grid.major = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) +
    scale_x_continuous(xlabs, breaks = times, limits = xlims) +
    scale_y_continuous(ylabs, limits = ylims, labels = scale_labels)
  
  
  #Add 95% CI to plot
  if(ci == TRUE)
    p <- p +  geom_ribbon(data=df, aes(ymin = lower, ymax = upper), fill = "grey", alpha=0.25, colour=NA)
  
  #Removes the legend:
  if(legend == FALSE)
    p <- p + theme(legend.position="none")
  
  #Add lines too plot
  p <- p + geom_step(size = 0.75) +
    scale_linetype_manual(name = ystrataname, values=linetype) +
    scale_colour_brewer(name = ystrataname, palette=linecols)
  
  
  ## p-value
  if(class(sfit) == "svykm") pval <- FALSE
  #if(is.null(design)) pval <- FALSE
  
  if(pval == TRUE) {
    if(is.null(design)){
      sdiff <- survey::svylogrank(formula(sfit), design = get(as.character(attr(sfit, "call")$design)))
    } else{
      sdiff <- survey::svylogrank(formula(sfit), design = design)
    }
    pvalue <- sdiff[[2]][2]
    
    pvaltxt <- ifelse(pvalue < 0.0001,"p < 0.0001",paste("p =", round(pvalue, 3)))
    if (pval.testname) pvaltxt <- paste0(pvaltxt, " (Log-rank)")
    
    # MOVE P-VALUE LEGEND HERE BELOW [set x and y]
    if (is.null(pval.coord)){
      p <- p + annotate("text",x = (as.integer(max(sapply(sfit, function(x){max(x$time)/5})))), y = 0.1 + ylims[1],label = pvaltxt, size  = pval.size)
    } else{
      p <- p + annotate("text",x = pval.coord[1], y = pval.coord[2], label = pvaltxt, size  = pval.size)
    }
  }
  
  ## Create a blank plot for place-holding
  blank.pic <- ggplot(df, aes(time, surv)) +
    geom_blank() + theme_void() +             ## Remove gray color
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),panel.border = element_blank())
  
  ###################################################
  # Create table graphic to include at-risk numbers #
  ###################################################
  
  n.risk <- NULL
  if(table == TRUE) {
    
    if(is.null(design)){
      sfit2 <- survival::survfit(formula(sfit), data = get(as.character(attr(sfit, "call")$design))$variables)
    } else{
      sfit2 <- survival::survfit(formula(sfit), data = design$variables)
    }
    
    #times <- seq(0, max(sfit2$time), by = timeby)
    
    if(is.null(subs)){
      if(length(levels(summary(sfit2)$strata)) == 0) {
        subs1 <- 1
        subs2 <- 1:length(summary(sfit2,censored=T)$time)
        subs3 <- 1:length(summary(sfit2,times = times,extend = TRUE)$time)
      } else {
        subs1 <- 1:length(levels(summary(sfit2)$strata))
        subs2 <- 1:length(summary(sfit2,censored=T)$strata)
        subs3 <- 1:length(summary(sfit2,times = times,extend = TRUE)$strata)
      }
    } else{
      for(i in 1:length(subs)){
        if(i==1){
          ssvar <- paste("(?=.*\\b=",subs[i],sep="")
        }
        if(i==length(subs)){
          ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],"\\b)",sep="")
        }
        if(!i %in% c(1, length(subs))){
          ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],sep="")
        }
        if(i==1 & i==length(subs)){
          ssvar <- paste("(?=.*\\b=",subs[i],"\\b)",sep="")
        }
      }
      subs1 <- which(regexpr(ssvar,levels(summary(sfit2)$strata), perl=T)!=-1)
      subs2 <- which(regexpr(ssvar,summary(sfit2,censored=T)$strata, perl=T)!=-1)
      subs3 <- which(regexpr(ssvar,summary(sfit2,times = times,extend = TRUE)$strata, perl=T)!=-1)
    }
    
    if(!is.null(subs)) pval <- FALSE
    

    
    if(length(levels(summary(sfit2)$strata)) == 0) {
      Factor <- factor(rep("All",length(subs3)))
    } else {
      Factor <- factor(summary(sfit2,times = times,extend = TRUE)$strata[subs3])
    }
    
    
    risk.data <- data.frame(
      strata = Factor,
      time = summary(sfit2,times = times,extend = TRUE)$time[subs3],
      n.risk = summary(sfit2,times = times,extend = TRUE)$n.risk[subs3]
    )
    
    
    risk.data$strata <- factor(risk.data$strata, levels=rev(levels(risk.data$strata)))
    
    data.table <- ggplot(risk.data,aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) + 
      geom_text(size = 3.5) + theme_bw() +
      scale_y_discrete(breaks = as.character(levels(risk.data$strata)),
                       labels = rev(ystratalabs)) +
      scale_x_continuous(label.nrisk, limits = xlims) +
      theme(axis.title.x = element_text(size = size.label.nrisk, vjust = 1),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),axis.text.x = element_blank(),
            axis.ticks = element_blank(),axis.text.y = element_text(face = "bold",hjust = 1)) 
    data.table <- data.table +
      theme(legend.position = "none") + xlab(NULL) + ylab(NULL)
    
    
    # ADJUST POSITION OF TABLE FOR AT RISK
    data.table <- data.table +
      theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 2.5, 3.5) - 0.15 * m), "lines"))
  }
  
  #######################
  # Plotting the graphs #
  #######################
  
  if(table == TRUE){
    grid.arrange(p, blank.pic, data.table, clip = FALSE, nrow = 3,
                 ncol = 1, heights = unit(c(2, .1, .25),c("null", "null", "null")))
  } else {
    p
  }

  
}



