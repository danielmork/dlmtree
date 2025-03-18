#' plot.summary.tdlm
#'
#' @title Plots a distributed lag function for model summary of 'tdlm'
#' @description Method for plotting a distributed lag function for model summary of 'tdlm'
#' 
#' @param x object of class 'summary.tdlm', output of summary of 'tdlm'
#' @param ... additional plotting parameters for title and labels
#' 'start.time' which sets the first time value
#'
#' @returns A plot of distributed lag effect estimated with tdlm
#' @export
#'
plot.summary.tdlm <- function(x, ...) {

  args <- list(...)
  main <- ifelse(!is.null(args$main), args$main, "DLM")
  xlab <- ifelse(!is.null(args$xlab), args$xlab, "Time")
  ylab <- ifelse(!is.null(args$ylab), args$ylab, "Effect")
  start.time  <- ifelse(!is.null(args$start.time), args$start.time, 1)

  d <- data.frame("Est" = x$matfit, 
                  "CIMin" = x$cilower,
                  "CIMax" = x$ciupper,
                  "X" = start.time:(start.time + length(x$matfit) - 1))

  p <- ggplot(d) +
        geom_hline(yintercept = 0, color = "red") +
        geom_ribbon(aes(x = X, ymin = CIMin, ymax = CIMax), fill = "grey", alpha = 0.7) +
        geom_line(aes(x = X, y = Est)) +
        theme_bw() +
        scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
        labs(x = xlab, y = ylab, title = main)
  
  return(p)
}
