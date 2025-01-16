# Functions for manipulation of IQR General Parameter Format (GPF) objects.

#' Remove parameters with NA in VALUE column from the GPF object
#' @param gpf a GPF object or filepath to an gpf file.
#' @return a GPF object
#'
#' @export
#' @family GPF
removeParametersWithNAValue_GPF <- function(gpf) {
  if(is.character(gpf)) {
    gpf <- suppressWarnings(GPF(gpf))
  } else if(!is_GPF(gpf)) {
    stop("removeParametersWithNAValue_GPF: argument gpf should be a gpf object or a file path to a GPF excel file.")
  }
  paramsWithNAValue <- gpf$estimates$PARAMETER[is.na(gpf$estimates$VALUE)]
  paramsToKeep <- setdiff(gpf$estimates$PARAMETER, paramsWithNAValue)
  if(length(paramsToKeep) == 0) {
    stop("removeParametersWithNAValue_GPF: the provided GPF object does not have non-NA values.")
  }
  gpf$estimates <- gpf$estimates[gpf$estimates$PARAMETER %in% paramsToKeep,]
  gpf$uncertainty_correlation <- gpf$uncertainty_correlation[gpf$uncertainty_correlation$PARAMETER %in% paramsToKeep,
                                                             colnames(gpf$uncertainty_correlation) %in% c("PARAMETER", paramsToKeep),
                                                             with = FALSE]
  gpf
}

#' Check if an object inherits from S3 class GPF
#' @param x any R object.
#' @return a logical
#' @author Venelin Mitov (IntiQuan)
#' @family GPF
is_GPF <- function (x) {
  inherits(x, "GPF")
}