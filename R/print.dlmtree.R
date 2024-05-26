#' Print a tdlnm Object
#'
#' @param x An object of class tdlnm.
#' @param ... Not used.
#'
#' @return Assorted model output.
#' @export
#'
print.tdlnm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))

}

#' Print a tdlm Object
#'
#' @param x An object of class tdlm.
#' @param ... Not used.
#'
#' @return Assorted model output.
#' @export
#'
print.tdlm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}


#' Print a tdlmm Object
#'
#' @param x An object of class tdlmm.
#' @param ... Not used.
#'
#' @return Assorted model output.
#' @export
#'
print.tdlmm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}

#' Print a hdlm Object
#'
#' @param x An object of class hdlm.
#' @param ... Not used.
#'
#' @return Assorted model output.
#' @export
#'
print.hdlm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}


#' Print a hdlmm Object
#'
#' @param x An object of class hdlmm.
#' @param ... Not used.
#'
#' @return Assorted model output.
#' @export
#'
print.hdlmm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}


#' Print a monotone Object
#'
#' @param x An object of class monotone
#' @param ... Not used.
#'
#' @return Assorted model output.
#' @export
#'
print.monotone <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}