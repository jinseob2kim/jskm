
#' @title Creates a Weighted Kaplan-Meier plot - svykm.object in survey package
#' @description Creates a Weighted Kaplan-Meier plot - svykm.object in survey package
#' @param sfit a svykm object
#' @param xlabs x-axis label, Default: 'Time-to-event'
#' @param ylabs y-axis label.
#' @param xlims numeric: list of min and max for x-axis. Default: NULL
#' @param ylims numeric: list of min and max for y-axis. Default: c(0, 1)
#' @param ystratalabs character list. A list of names for each strata. Default: NULL
#' @param ystrataname The legend name. Default: 'Strata'
#' @param timeby numeric: control the granularity along the time-axis; defaults to 7 time-points. 
#' @param main plot title, Default: ''
#' @param pval logical: add the pvalue to the plot?, Default: FALSE
#' @param legend logical. should a legend be added to the plot? Default: TRUE
#' @param legendposition numeric. x, y position of the legend if plotted. Default: c(0.85, 0.8)
#' @param linecols Character. Colour brewer pallettes too colour lines. Default: 'Set1'
#' @param dashed logical. Should a variety of linetypes be used to identify lines. Default: FALSE
#' @param cumhaz Show cumulaive hazard function, Default: F
#' @param design Data design for reactive design data , Default: NULL
#' @param ... PARAM_DESCRIPTION
#' @return plot
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  data(pbc, package="survival")
#'  pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
#'  biasmodel<-glm(randomized~age*edema,data=pbc)
#'  pbc$randprob<-fitted(biasmodel)
#'  
#'  dpbc<-svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))
#'  
#'  s1 <- svykm(Surv(time,status>0)~sex, design=dpbc)
#'  svyjskm(s1)
#'  }
#' }
#' @rdname svyjskm
#' @import ggplot2
#' @importFrom stats formula
#' @importFrom survey svyranktest
#' @export 

svyjskm <- function(sfit,
                    xlabs = "Time-to-event",
                    ylabs = "Survival (%)",
                    xlims = NULL,
                    ylims = c(0,1),
                    ystratalabs = NULL,
                    ystrataname = "Strata",
                    timeby = NULL,
                    main = "",
                    pval = FALSE,
                    legend = TRUE,
                    legendposition=c(0.85,0.8),
                    linecols="Set1",
                    dashed= FALSE,
                    cumhaz = F,
                    design = NULL,
                    ...) {

  if(is.null(ystrataname)) ystrataname <- "Strata"
  
  
  if (is.null(timeby)){
    if (class(sfit) == "svykmlist"){
      timeby <- signif(max(sapply(sfit, function(x){max(x$time)}))/7, 1)
    } else if(class(sfit) == "svykm"){
      
      timeby <- signif(max(sfit$time)/7, 1)
      }
  }
  
  if (class(sfit) == "svykmlist"){
    df <- Reduce(rbind, lapply(names(sfit), function(x){data.frame("strata" = x, "time" = sfit[[x]]$time, "surv" = sfit[[x]]$surv)}))
    df$strata <- as.factor(df$strata)
    times <- seq(0, max(sapply(sfit, function(x){max(x$time)})), by = timeby)
    if (is.null(ystratalabs)){
      ystratalabs <- names(sfit)
    }
    if (is.null(xlims)){
      xlims <- c(0,max(sapply(sfit, function(x){max(x$time)})))
    }
    
  } else if(class(sfit) == "svykm"){
    df <- data.frame("strata" = "All", "time" = sfit$time, "surv" = sfit$surv)
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
    df$surv = 1 - df$surv
    }
  
  #Final changes to data for survival plot
  levels(df$strata) <- ystratalabs
  zeros <- data.frame("strata" = factor(ystratalabs, levels=levels(df$strata)), "time" = 0, "surv" = 1)
  
  if (cumhaz){
    zeros$surv = 0
  }
  
  df <- rbind(zeros, df)
  d <- length(levels(df$strata))
  
  ###################################
  # specifying axis parameteres etc #
  ###################################
  
  if(dashed == TRUE){
    linetype=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")
  } else {
    linetype=c("solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")
  }
  
  
  p <- ggplot( df, aes(x=time, y=surv, colour=strata, linetype=strata)) +
    ggtitle(main)
  
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
    scale_y_continuous(ylabs, limits = ylims)
  
  
  
  #Removes the legend:
  if(legend == FALSE)
    p <- p + theme(legend.position="none")
  
  #Add lines too plot
  p <- p + geom_step(size = 0.75) +
    scale_linetype_manual(name = ystrataname, values=linetype) +
    scale_colour_brewer(name = ystrataname, palette=linecols)
  
  
  ## p-value
  if(class(sfit) == "svykm") pval <- FALSE
  if(is.null(design)) pval <- FALSE
  
  if(pval == TRUE) {
    
    sdiff <- survey::svyranktest(formula(sfit), design = design)
    pvalue <- sdiff$p.value
    
    pvaltxt <- ifelse(pvalue < 0.0001,"p < 0.0001",paste("p =", signif(pvalue, 3)))
    # MOVE P-VALUE LEGEND HERE BELOW [set x and y]
    p <- p + annotate("text",x = (as.integer(max(sapply(sfit, function(x){max(x$time)}))/5)), y = 0.1,label = pvaltxt)
  }
  
  p

  
}



