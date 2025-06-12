#' Retrieve Package Path
#'
#' @description This function retrieves the path to the MMVmalaria package directory
#' @param subdir A subdirectory within the MMVmalaria package directory. If NULL, the function returns the main package path.
#' @return A character string representing the path to the MMVmalaria package directory or a specified subdirectory.
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

