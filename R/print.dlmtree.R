#' @method print tdlnm
#' @rdname print
#'
#' @export
print.tdlnm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}

#' @method print tdlm
#' @rdname print
#'
#' @export
print.tdlm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}


#' @method print tdlmm
#' @rdname print
#'
#' @export
print.tdlmm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}

#' @method print hdlm
#' @rdname print
#'
#' @export
print.hdlm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}



#' @method print hdlmm
#' @rdname print
#'
#' @export
print.hdlmm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}


#' @method print monotone
#' @rdname print
#'
#' @export
print.monotone <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}