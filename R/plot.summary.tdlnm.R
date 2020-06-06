#' plot.summary.tdlnm
#'
#' @param summary
#' @param ptype
#' @param var
#' @param lag
#' @param ...
#'
#' @return
#' @export
#' @import ggplot2
#' @import viridis
#'
#' @examples
plot.summary.tdlnm <- function(summary, ptype = "mean", var = c(), lag = c(), ...)
{
  args <- list(...)
  main <- ifelse(!is.null(args$main), args$main, "")
  xlab <- ifelse(!is.null(args$xlab), args$xlab, "Time")
  ylab <- ifelse(!is.null(args$ylab), args$ylab, "Exposure-Concentration")
  start.time <- ifelse(!is.null(args$start.time), args$start.time, 1)
  if (ptype == "mean") {
    flab <- ifelse(!is.null(args$flab), args$flab, "Est Effect")
    p <- ggplot(summary$dlnm.estimates, aes(xmin = Tmin + start.time,
                                            xmax = Tmax + start.time,
                                           ymin = Xmin,
                                           ymax = Xmax, fill = Est)) +
      geom_rect() +
      scale_color_viridis(aesthetics = "fill", option = "D") +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (ptype == "se") {
    flab <- ifelse(!is.null(args$flab), args$flab, "SE Effect")
    p <- ggplot(summary$dlnm.estimates, aes(xmin = Tmin + start.time,
                                            xmax = Tmax + start.time,
                                           ymin = Xmin,
                                           ymax = Xmax, fill = SD)) +
      geom_rect() +
      scale_color_viridis(aesthetics = "fill", option = "B") +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (ptype == "ci-min") {
    flab <- ifelse(!is.null(args$flab), args$flab, "CI Min")
    p <- ggplot(summary$dlnm.estimates, aes(xmin = Tmin, xmax = Tmax,
                                           ymin = Xmin,
                                           ymax = Xmax, fill = CIMin)) +
      geom_rect() +
      scale_color_viridis(aesthetics = "fill", option = "A") +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (ptype == "ci-max") {
    flab <- ifelse(!is.null(args$flab), args$flab, "CI Max")
    p <- ggplot(summary$dlnm.estimates, aes(xmin = Tmin + start.time,
                                            xmax = Tmax + start.time,
                                           ymin = Xmin,
                                           ymax = Xmax, fill = CIMax)) +
      geom_rect() +
      scale_color_viridis(aesthetics = "fill", option = "A") +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      labs(x = xlab, y = ylab, fill = flab, title = main)
  } else if (ptype == "animate") {
    if (!require(gganimate))
     stop("Package gganimate required.")
    summary$dlnm.estimates$Xmin <- round(summary$dlnm.estimates$Xmin, 2)
    p <- ggplot(summary$dlnm.estimates) +
      geom_hline(yintercept = 0, color = "red") +
      geom_ribbon(aes(x = Tmin + start.time, ymin = CIMin, ymax = CIMax), fill = "grey") +
      geom_line(aes(x = Tmin + start.time, y = Est)) +
      transition_states(Xmin, transition_length = 1, state_length = 0) +
      theme_bw() +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      labs(x = xlab, y = ylab) +
      ggtitle('Exposure {closest_state}')
    return(animate(p, duration = 20, fps = 10))
  } else if (ptype == "effect") {
    flab <- ifelse(!is.null(args$flab), args$flab, "Effect")
    p <- ggplot(summary$dlnm.estimates, aes(xmin = Tmin + start.time,
                                            xmax = Tmax + start.time,
                                       ymin = Xmin,
                                       ymax = Xmax, fill = Effect)) +
      geom_rect() +
      scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
      theme_bw() +
      scale_fill_manual(breaks = c("+", "-", " "),
                        values = c("+" = "dodgerblue1", "-" = "tomato1", " " = "white")) +
      labs(x = xlab, y = ylab, fill = flab, title = main)

  } else if (ptype == "slice") {
    if (length(var) != 0) {
      main <- ifelse(!is.null(args$main), args$main, paste0("Exposure = ", var))
      xlab <- ifelse(!is.null(args$xlab), args$xlab, "Time")
      ylab <- ifelse(!is.null(args$ylab), args$ylab, "Est Effect")
      idx <- which(summary$dlnm.estimates$Xmin <= var &
                     summary$dlnm.estimates$Xmax > var)
      p <- ggplot(summary$dlnm.estimates[idx,]) +
        geom_hline(yintercept = 0, color = "red") +
        geom_ribbon(aes(x = Tmin + start.time, ymin = CIMin, ymax = CIMax), fill = "grey") +
        geom_line(aes(x = Tmin + start.time, y = Est)) +
        theme_bw() +
        scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
        labs(x = xlab, y = ylab, title = main)
    } else if (length(time) != 0) {
      main <- ifelse(!is.null(args$main), args$main, paste0("Time = ", lag))
      xlab <- ifelse(!is.null(args$xlab), args$xlab, "Exposure")
      ylab <- ifelse(!is.null(args$ylab), args$ylab, "Est Effect")
      idx <- which(summary$dlnm.estimates$Tmin + start.time <= lag &
                     summary$dlnm.estimates$Tmax + start.time > lag)
      p <- ggplot(summary$dlnm.estimates[idx,]) +
        geom_hline(yintercept = 0, color = "red") +
        geom_ribbon(aes(x = Xmin, ymin = CIMin, ymax = CIMax), fill = "grey") +
        geom_line(aes(x = Xmin, y = Est)) +
        theme_bw() +
        scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
        labs(x = xlab, y = ylab, title = main)
    }
  }
  return(p)
}
