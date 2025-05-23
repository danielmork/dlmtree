#' plot.summary.monotone
#'
#' @title Returns variety of plots for model summary of class 'monotone'
#' @description Method for returning variety of plots for model summary of class 'monotone'
#'
#' @param x object of class 'summary.monotone', output of summary of 'monotone'
#' @param plot.type string indicating plot type, options are 'mean' (default)
#' which shows mean exposure-time response surface, 'se', 'ci-min', 'ci-max',
#' 'slice' which takes a slice of the plot at a given 'val' or 'time',
#' 'animate' which creates a animation of slices of the surface plot across
#' exposure values (requires package gganimate)
#' @param val exposure value for slice plot
#' @param time time value for slice plot
#' @param ... additional parameters to alter plots: 'main', 'xlab', 'ylab',
#' 'flab' which sets the effect label for surface plots,
#' 'start.time' which sets the first time value
#' 
#' @returns A plot of distributed lag effect estimated with monotone-TDLNM
#' @export
#'
plot.summary.monotone <- function(x, plot.type = "mean", val = c(), time = c(), ...)
{
  args <- list(...)
  main <- ifelse(!is.null(args$main), args$main, "")
  xlab <- ifelse(!is.null(args$xlab), args$xlab, "Time")
  ylab <- ifelse(!is.null(args$ylab), args$ylab, "Exposure-Concentration")
  start.time <- ifelse(!is.null(args$start.time), args$start.time, 1)
  if (plot.type == "mean") {
    flab <- ifelse(!is.null(args$flab), args$flab, "Est Effect")
    p <- ggplot(x$plot.dat, aes(xmin = Tmin + start.time,
                                      xmax = Tmax + start.time,
                                      ymin = Xmin,
                                      ymax = Xmax, fill = Est)) +
      geom_rect() +
      scale_fill_viridis_c() +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (plot.type == "se") {
    flab <- ifelse(!is.null(args$flab), args$flab, "SE Effect")
    p <- ggplot(x$plot.dat, aes(xmin = Tmin + start.time,
                                      xmax = Tmax + start.time,
                                      ymin = Xmin,
                                      ymax = Xmax, fill = SD)) +
      geom_rect() +
      scale_fill_viridis_c() +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (plot.type == "ci-min") {
    flab <- ifelse(!is.null(args$flab), args$flab, "CI Min")
    p <- ggplot(x$plot.dat, aes(xmin = Tmin, xmax = Tmax,
                                      ymin = Xmin,
                                      ymax = Xmax, fill = CIMin)) +
      geom_rect() +
      scale_fill_viridis_c() +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (plot.type == "ci-max") {
    flab <- ifelse(!is.null(args$flab), args$flab, "CI Max")
    p <- ggplot(x$plot.dat, aes(xmin = Tmin + start.time,
                                      xmax = Tmax + start.time,
                                      ymin = Xmin,
                                      ymax = Xmax, fill = CIMax)) +
      geom_rect() +
      scale_fill_viridis_c() +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (plot.type == "cumulative") {
    xlab <- ifelse(!is.null(args$xlab), args$xlab, "Exposure-concentration")
    ylab <- ifelse(!is.null(args$ylab), args$ylab, "Cumulative effect")
    p <- ggplot(x$cumulative.effect, aes(x = vals, y = mean, ymin = lower, ymax = upper)) +
      geom_hline(yintercept = 0, color = "red") +
      geom_ribbon(fill = "grey", alpha = 0.7) +
      geom_line() +
      theme_bw() +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      labs(x = xlab, y = ylab, title = main)
  } else if (plot.type == "effect") {
    flab <- ifelse(!is.null(args$flab), args$flab, "Effect")
    p <- ggplot(x$plot.dat, aes(xmin = Tmin + start.time,
                                      xmax = Tmax + start.time,
                                      ymin = Xmin,
                                      ymax = Xmax, fill = Effect)) +
      geom_rect() +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      scale_fill_manual(breaks = c("+", "-", " "),
                        values = c("+" = "dodgerblue1", "-" = "tomato1", " " = "white")) +
      labs(x = xlab, y = ylab, fill = flab, title = main)

  } else if (plot.type == "slice") {
    if (length(val) != 0) {
      main <- ifelse(!is.null(args$main), args$main, paste0("Exposure = ", val))
      xlab <- ifelse(!is.null(args$xlab), args$xlab, "Time")
      ylab <- ifelse(!is.null(args$ylab), args$ylab, "Est Effect")
      idx <- which(x$plot.dat$Xmin <= val &
                     x$plot.dat$Xmax > val)
      p <- ggplot(x$plot.dat[idx,]) +
        geom_hline(yintercept = 0, color = "red") +
        geom_ribbon(aes(x = Tmin + start.time, ymin = CIMin, ymax = CIMax), fill = "grey", alpha = 0.7) +
        geom_line(aes(x = Tmin + start.time, y = Est)) +
        theme_bw() +
        scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
        labs(x = xlab, y = ylab, title = main)
    } else if (length(time) != 0) {
      main <- ifelse(!is.null(args$main), args$main, paste0("Time = ", time))
      xlab <- ifelse(!is.null(args$xlab), args$xlab, "Exposure")
      ylab <- ifelse(!is.null(args$ylab), args$ylab, "Est Effect")
      idx <- which(x$plot.dat$Tmin + start.time <= time &
                     x$plot.dat$Tmax + start.time > time)
      p <- ggplot(x$plot.dat[idx,]) +
        geom_hline(yintercept = 0, color = "red") +
        geom_ribbon(aes(x = PredVal, ymin = CIMin, ymax = CIMax), fill = "grey", alpha = 0.7) +
        geom_line(aes(x = PredVal, y = Est)) +
        theme_bw() +
        scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
        labs(x = xlab, y = ylab, title = main)
    }
  } 

  return(p)
}