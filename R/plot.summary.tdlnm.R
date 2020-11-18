#' plot.summary.tdlnm
#'
#' @param summary object of class 'summary.tdlnm', output of summary of 'tdlnm'
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
#' @return
#' @export
#' @import ggplot2
#' @import viridis
#'
plot.summary.tdlnm <- function(summary, plot.type = "mean", val = c(), time = c(), ...)
{
  if (summary$ctr$dl.function == "tdlm")
    plot.type = "dlm"
  args <- list(...)
  main <- ifelse(!is.null(args$main), args$main, "")
  xlab <- ifelse(!is.null(args$xlab), args$xlab, "Time")
  ylab <- ifelse(!is.null(args$ylab), args$ylab, "Exposure-Concentration")
  start.time <- ifelse(!is.null(args$start.time), args$start.time, 1)
  if (plot.type == "mean") {
    flab <- ifelse(!is.null(args$flab), args$flab, "Est Effect")
    p <- ggplot(summary$plot.dat, aes(xmin = `Tmin` + start.time,
                                      xmax = `Tmax` + start.time,
                                      ymin = `Xmin`,
                                      ymax = `Xmax`, fill = `Est`)) +
      geom_rect() +
      scale_color_viridis(aesthetics = "fill", option = "D") +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (plot.type == "se") {
    flab <- ifelse(!is.null(args$flab), args$flab, "SE Effect")
    p <- ggplot(summary$plot.dat, aes(xmin = `Tmin` + start.time,
                                      xmax = `Tmax` + start.time,
                                      ymin = `Xmin`,
                                      ymax = `Xmax`, fill = `SD`)) +
      geom_rect() +
      scale_color_viridis(aesthetics = "fill", option = "B") +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (plot.type == "ci-min") {
    flab <- ifelse(!is.null(args$flab), args$flab, "CI Min")
    p <- ggplot(summary$plot.dat, aes(xmin = `Tmin`, xmax = `Tmax`,
                                      ymin = `Xmin`,
                                      ymax = `Xmax`, fill = `CIMin`)) +
      geom_rect() +
      scale_color_viridis(aesthetics = "fill", option = "A") +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (plot.type == "ci-max") {
    flab <- ifelse(!is.null(args$flab), args$flab, "CI Max")
    p <- ggplot(summary$plot.dat, aes(xmin = `Tmin` + start.time,
                                      xmax = `Tmax` + start.time,
                                      ymin = `Xmin`,
                                      ymax = `Xmax`, fill = `CIMax`)) +
      geom_rect() +
      scale_color_viridis(aesthetics = "fill", option = "A") +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  # } else if (plot.type == "animate") {
  #   if (!require(gganimate))
  #     stop("Package gganimate required.")
  #   if (!require(transformr))
  #     stop("Package transformr required.")
  #   summary$plot.dat$Xmin <- round(summary$plot.dat$Xmin, 2)
  #   p <- ggplot(summary$plot.dat) +
  #     geom_hline(yintercept = 0, color = "red") +
  #     geom_ribbon(aes(x = Tmin + start.time, ymin = CIMin, ymax = CIMax), fill = "grey") +
  #     geom_line(aes(x = Tmin + start.time, y = Est)) +
  #     transition_states(Xmin, transition_length = 1, state_length = 0) +
  #     theme_bw() +
  #     scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  #     labs(x = xlab, y = ylab) +
  #     ggtitle('Exposure {closest_state}')
  #   return(animate(p, duration = 20, fps = 10))
  } else if (plot.type == "effect") {
    flab <- ifelse(!is.null(args$flab), args$flab, "Effect")
    p <- ggplot(summary$plot.dat, aes(xmin = `Tmin` + start.time,
                                      xmax = `Tmax` + start.time,
                                      ymin = `Xmin`,
                                      ymax = `Xmax`, fill = `Effect`)) +
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
      idx <- which(summary$plot.dat$Xmin <= val &
                     summary$plot.dat$Xmax > val)
      p <- ggplot(summary$plot.dat[idx,]) +
        geom_hline(yintercept = 0, color = "red") +
        geom_ribbon(aes(x = `Tmin` + start.time, ymin = `CIMin`, ymax = `CIMax`), fill = "grey") +
        geom_line(aes(x = `Tmin` + start.time, y = `Est`)) +
        theme_bw() +
        scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
        labs(x = xlab, y = ylab, title = main)
    } else if (length(time) != 0) {
      main <- ifelse(!is.null(args$main), args$main, paste0("Time = ", time))
      xlab <- ifelse(!is.null(args$xlab), args$xlab, "Exposure")
      ylab <- ifelse(!is.null(args$ylab), args$ylab, "Est Effect")
      idx <- which(summary$plot.dat$Tmin + start.time <= time &
                     summary$plot.dat$Tmax + start.time > time)
      p <- ggplot(summary$plot.dat[idx,]) +
        geom_hline(yintercept = 0, color = "red") +
        geom_ribbon(aes(x = `PredVal`, ymin = `CIMin`, ymax = `CIMax`), fill = "grey") +
        geom_line(aes(x = `PredVal`, y = `Est`)) +
        theme_bw() +
        scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
        labs(x = xlab, y = ylab, title = main)
    }
  } else if (plot.type == "dlm") {
    main <- ifelse(!is.null(args$main), args$main, "DLM")
    xlab <- ifelse(!is.null(args$xlab), args$xlab, "Time")
    ylab <- ifelse(!is.null(args$ylab), args$ylab, "Effect")
    d <- data.frame("Est" = summary$matfit, "CIMin" = summary$cilower,
                    "CIMax" = summary$ciupper,
                    "X" = start.time:(start.time + length(summary$matfit) - 1))
    p <- ggplot(d) +
      geom_hline(yintercept = 0, color = "red") +
      geom_ribbon(aes(x = `X`, ymin = `CIMin`, ymax = `CIMax`), fill = "grey") +
      geom_line(aes(x = `X`, y = `Est`)) +
      theme_bw() +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      labs(x = xlab, y = ylab, title = main)
  }
  return(p)
}
