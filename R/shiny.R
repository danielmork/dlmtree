#' shiny
#'
#' @description shiny generic function for S3method
#'
#' @param x an object to which S3method is applied
#'
#' @export shiny
shiny <- function(x, ...) {
  UseMethod("shiny")
}