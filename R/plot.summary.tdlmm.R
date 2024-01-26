#' plot.summary.tdlmm
#'
#' @description Method for plotting DLMMs from summary of tdlmm. Includes plots
#' for marginal exposure effects as well as interactions between two exposures.
#'
#' @param object object of type 'summary.tdlmm' from summary.tdlmm() output
#' @param type plot type, 'marginal' (default)
#' @param exposure1 exposure for plotting DLM
#' @param exposure2 exposure paired with 'exposure1' for plotting interaction
#' @param time1 plot a cross section from an interaction plot at specific time
#' for 'exposure1'
#' @param time2 plot a cross section from an interaction plot at specific time
#' for 'exposure2'
#' @param show.cw indicate location of critical windows in interaction plot
#' with red points
#' @param cw.plots.only show only plots with critical windows
#' @param ... additional plotting parameters for manipulating title and labels
#'
#' @return ggplot
#' @export plot.summary.tdlmm
#' @export
#'
plot.summary.tdlmm <- function(object,
                               type = "marginal",
                               exposure1 = NULL,
                               exposure2 = NULL,
                               time1 = c(),
                               time2 = c(),
                               show.cw = TRUE,
                               cw.plots.only = TRUE,
                               trueDLM = NULL,
                               scale = NULL,
                               ...)
{
  # cycle through plotting all exposures and interactions if none specified
  if (is.null(exposure1)) {
    if (type == "marginal") {
      cat("Plotting DLM marginal effects:\n")
    } #else if (type == "nonlinear") {
    #cat("Plotting DLM marginal nonlinear effects:\n")
    #}
    for (ex.name in object$expNames) {
      if (!cw.plots.only | any(object$TreeStructs[[ex.name]]$marg.cw)) {
        plot(plot.summary.tdlmm(object, type, ex.name, NULL, time1, time2, ...))
        readline(prompt = "Press [enter] to continue")
      }
    }

    if (length(object$mixNames) > 0 & object$interaction > 0) {
      cat("Plotting interaction effects:\n")
      for (ex.name1 in object$expNames) {
        for (ex.name2 in object$expNames) {
          if (paste0(ex.name1, "-", ex.name2) %in% names(object$MIX)) {
            if (!cw.plots.only | any(object$MIX[[paste0(ex.name1, "-", ex.name2)]]$cw)) {
              plot(plot.summary.tdlmm(object, type, ex.name1,
                                      ex.name2, time1, time2, ...))
              readline(prompt = "Press [enter] to continue")
            }
          }
        }
      }
    }
    return(invisible())
  }


  # Plot setup
  args <- list(...)
  start.time <- ifelse(!is.null(args$start.time), args$start.time, 1)
  Lags <- start.time:(start.time + object$nLags - 1)
  base_size <- ifelse(!is.null(args$base_size), args$base_size, 11)

  if (is.numeric(exposure1)) {
    if (exposure1 > length(object$expNames)) {
      stop("exposure1 incorrectly specified")
    }
      
    exposure1 <- object$expNames[exposure1]
  }

  # Plot DLM
  if (is.null(exposure2)) {
    xlab <- ifelse(!is.null(args$xlab), args$xlab, "Time")
    ylab <- ifelse(!is.null(args$ylab), args$ylab, "Effect")

    if (type == "marginal") {
      main <- ifelse(!is.null(args$main), args$main,
                     paste0("Marginal effect: ", exposure1))

      dat <- data.frame("Est" = object$TreeStructs[[exposure1]]$marg.matfit,
                        "CIMin" = object$TreeStructs[[exposure1]]$marg.cilower,
                        "CIMax" = object$TreeStructs[[exposure1]]$marg.ciupper,
                        "X" = Lags)

      if (!is.null(scale)) {  # Scaling for dlm
        dat[, c("Est", "CIMin", "CIMax")] <- exp(dat[, c("Est", "CIMin", "CIMax")])
      }

      if (!is.null(trueDLM)) {  # SI: df for a plot returning trueDLM
        if (is.null(scale)) {
          dat$trueDLM <- trueDLM
        } else {
          dat$trueDLM <- exp(trueDLM)
        }
      }
    }
    


    if (!is.null(trueDLM)) { # draws an additional line of true DLM effect
      p <- ggplot(dat) +
        geom_hline(yintercept = ifelse(is.null(scale), 0, 1), color = "red", linetype = "dashed") +
        geom_ribbon(aes(x = `X`, ymin = `CIMin`, ymax = `CIMax`), fill = "grey") +
        geom_line(aes(x = `X`, y = `Est`)) +
        geom_line(aes(x = `X`, y = `trueDLM`), col = "blue", linetype = "dashed") + # SI
        theme_bw(base_size = base_size) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        labs(x = xlab, y = ylab, title = main)
    } else { # No true DLM (original code)
      p <- ggplot(dat) +
        geom_hline(yintercept = ifelse(is.null(scale), 0, 1), color = "red", linetype = "dashed") +
        geom_ribbon(aes(x = `X`, ymin = `CIMin`, ymax = `CIMax`), fill = "grey") +
        geom_line(aes(x = `X`, y = `Est`)) +
        theme_bw(base_size = base_size) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        labs(x = xlab, y = ylab, title = main)
    }

    return(p)


    # Plot Mixture
  } else {
    if (is.numeric(exposure2)) {
      if (exposure2 > length(object$expNames)){
        stop("exposure2 incorrectly specified")
      }
        
      exposure2 <- object$expNames[exposure2]
    }

    plotDat <- data.frame(x = rep(Lags, length(Lags)),
                          y = rep(Lags, each = length(Lags)))

    if (paste0(exposure1, "-", exposure2) %in% names(object$MIX)){
      mix.name <- paste0(exposure1, "-", exposure2)
    } else if (paste0(exposure2, "-", exposure1) %in% names(object$MIX)){
      mix.name <- paste0(exposure2, "-", exposure1)
    } else {
      stop("mixture not found, check exposures names or numbers")
    }

    plotDat <- cbind.data.frame(plotDat,
                                Effect = c(object$MIX[[mix.name]]$matfit),
                                CW = c(object$MIX[[mix.name]]$cw.plot))

    main <- ifelse(!is.null(args$main), args$main, "Interaction effect")
    xlab <- ifelse(!is.null(args$xlab), args$xlab, paste0("Time: ",exposure1))
    ylab <- ifelse(!is.null(args$ylab), args$ylab, paste0("Time: ",exposure2))

    if (!is.null(args$cw.point.range)) {
      cw.point.range <- args$cw.point.range
    } else {
      cw.point.range <- c(0.1, 3)
    }
    cw.lab <- ifelse(!is.null(args$cw.lab), args$cw.lab, "Conf. Level")
    cw.lab.show <- ifelse(!is.null(args$cw.lab.show), args$cw.lab.show, FALSE)

    p <- ggplot(plotDat, aes(x = `x`, y = `y`, z = `Effect`, fill = `Effect`)) +
      geom_tile() +
      scale_fill_viridis()

    if (show.cw){
      p <- p + geom_point(data = plotDat[which(plotDat$CW != 0),],
                          aes(x = `x`, y = `y`, size = `CW`),
                          color = "red", show.legend = cw.lab.show) +
      scale_size_continuous(name = cw.lab, range = cw.point.range)
    }
      
    p <- p +
      theme_bw(base_size = base_size) +
      coord_equal() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = xlab, y = ylab, title = main)

    return(p)
  }
}
