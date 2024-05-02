#' ppRange
#'
#' @title Makes a 'pretty' output of a group of numbers
#' @description Method for making a 'pretty' output of a group of numbers. For
#' example: 2,3,4,5,8,9,12,15,16 becomes 2-5,8-9,12,15-16
#'
#' @param r set of integers to make 'pretty'
#'
#' @returns character string of values representing 'r'
#'
ppRange <- function(r) {
  if (length(r) <= 1){
    return(paste0(r))
  }

  rMin  <- min(r)
  rMax  <- max(r)
  rAll  <- rMin:rMax
  rMiss <- which(!(rAll %in% r))

  if (length(rMiss) == 0){
    return(paste0(rMin, "-", rMax))
  }
    
  rOut  <- ""
  rMiss <- c(0, rMiss, length(rAll) + 1)
  for (s in 1:(length(rMiss) - 1)) {
    if (rMiss[s + 1] == rMiss[s] + 1){
      next;
    }
      
    if (rMiss[s + 1] - rMiss[s] == 2) {
      rOut <- paste0(rOut, rAll[rMiss[s] + 1])
    } else {
      rOut <- paste0(rOut, rAll[rMiss[s] + 1], "-", rAll[rMiss[s + 1] - 1])
    }
      
    if (s < length(rMiss) - 1) {
      rOut <- paste0(rOut, ",")
    }
  }
  
  return(rOut)
}
