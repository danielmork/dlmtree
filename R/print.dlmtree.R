#' Print method for model class 'tdlm', 'tdlmm', 'tdlnm', 'hdlm', 'hdlmm', 'monotone'
#' @method print tdlnm
#' @rdname print
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

#' @method print tdlm
#' @rdname print
#' 
#' @param x An object of class tdlm.
#' @param ... Not used.
#'
#' @export
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
#' @param x An object of class tdlmm.
#' @param ... Not used.
#'
#' @export
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
#' @param x An object of class hdlm.
#' @param ... Not used.
#'
#' @export
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
#' @param x An object of class hdlmm.
#' @param ... Not used.
#'
#' @export
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
#' @param x An object of class monotone
#' @param ... Not used.
#'
#' @export
#'
print.monotone <- function(x, ...){
  
  cat("Object of class",x$class)
  cat("\n\nCall:\n")
  print(x$call)
  cat("\nAvailable methods:",paste(methods(class=x$class), sep=", "))
  
}