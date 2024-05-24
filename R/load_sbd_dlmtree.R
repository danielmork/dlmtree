
load_sbd_dlmtree <- function() {
  # Construct the raw file URL

  load(url("https://github.com/danielmork/dlmtree/blob/master/data/sbd_dlmtree.rda"))
    #   # Download the file
    #   download.file(url(), destfile = "sbd_dlmtree.rda")
    
    #   # Load the data into R
    #   data <- load("sbd_dlmtree.rda")
    
    #   # Delete the temporary file
    #   file.remove("sbd_dlmtree.rda")
    
    #   return(data)
}