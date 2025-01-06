#' Retrieve Package Path
#'
#' @description
#' @param subdir
#' @return
#' @export
#' @author Mohammed H. Cherkaoui (MMV, \email{cherkaouim@@mmv.org})
#' @family General Functions
get_MMVmalariaPath <- function(subdir = NULL){

  # Depending on how MMVmalaria was installed:
  #   Sources for development:
  if(exists("MMVmalariaPathDev")){
    if(is.null(subdir)){
      out <- file.path(MMVmalariaPathDev)
    }else{
      out <- file.path(MMVmalariaPathDev, subdir)
    }

    #   Installed as a package
  }else if("MMVmalaria" %in% .packages()){
    out <- system.file("",package = "MMVmalaria")

  }else {
    stop("Neither 'MMVmalaria' is installed as a package or 'MMVmalariaPathDev' exist")
  }

  # Output:
  return(out)
}

