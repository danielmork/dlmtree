###
### R routines for the R package dlmtree (c)
#


.onAttach <-
  function(lib, pkg) {
    #
    ################################################################################
    #
    meta <- packageDescription("dlmtree")
    attachmsg <- paste("This is dlmtree ",meta$Version,
                       ". For details visit https://danielmork.github.io/dlmtree/.",
                       sep="")
    packageStartupMessage(attachmsg, domain = NULL, appendLF = TRUE)
  }
