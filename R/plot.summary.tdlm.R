#' plot.summary.tdlnm
#'
#' @param x object of class 'summary.tdlnm', output of summary of 'tdlnm'
#' @param trueDLM A vector of true effects that can be obtained from the simulated data. Only applicable for simulation studies
#' @param ... additional plotting parameters for title and labels
#' 'flab' which sets the effect label for surface plots,
#' 'start.time' which sets the first time value
#'
#' @import ggplot2
#' @export plot.summary.tdlm
#' @export
#'
plot.summary.tdlm <- function(x,  trueDLM = NULL, ...) {

  args        <- list(...)
  main <- ifelse(!is.null(args$main), args$main, "DLM")
  xlab <- ifelse(!is.null(args$xlab), args$xlab, "Time")
  ylab <- ifelse(!is.null(args$ylab), args$ylab, "Effect")
  start.time  <- ifelse(!is.null(args$start.time), args$start.time, 1)

  if (!is.null(trueDLM)) {
    d <- data.frame("Est" = x$matfit, 
                    "CIMin" = x$cilower,
                    "CIMax" = x$ciupper,
                    "X" = start.time:(start.time + length(x$matfit) - 1),
                    "trueDLM" = trueDLM)
  } else {
    d <- data.frame("Est" = x$matfit, 
                    "CIMin" = x$cilower,
                    "CIMax" = x$ciupper,
                    "X" = start.time:(start.time + length(x$matfit) - 1))
  }
  
  if (!is.null(trueDLM)) {
    p <- ggplot(d) +
          geom_hline(yintercept = 0, color = "red") +
          geom_ribbon(aes(x = X, ymin = CIMin, ymax = CIMax), fill = "grey") +
          geom_line(aes(x = X, y = Est)) +
          geom_line(aes(x = X, y = trueDLM), col = "blue", linetype = "dashed") + # SI
          theme_bw() +
          scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
          labs(x = xlab, y = ylab, title = main)
  } else {
    p <- ggplot(d) +
          geom_hline(yintercept = 0, color = "red") +
          geom_ribbon(aes(x = X, ymin = CIMin, ymax = CIMax), fill = "grey") +
          geom_line(aes(x = X, y = Est)) +
          theme_bw() +
          scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
          labs(x = xlab, y = ylab, title = main)
  }
  
  
  return(p)
}