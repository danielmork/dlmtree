#' @method print tdlnm
#' @rdname print
#'
print.tdlnm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}

#' @method print tdlm
#' @rdname print
#'
print.tdlm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}


#' @method print tdlmm
#' @rdname print
#'
print.tdlmm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}

#' @method print hdlm
#' @rdname print
#'
print.hdlm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}



#' @method print hdlmm
#' @rdname print
#'
print.hdlmm <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}


#' @method print monotone
#' @rdname print
#'
print.monotone <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}