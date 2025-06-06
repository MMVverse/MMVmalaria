#' Sample NLME parameters from their uncertainty distribution
#'
#' @inheritParams sampleIndParamValues_MMVmalariaXLS
#'
#' @return a data.frame of `Npop` rows with columns `ID.POP` and parameter
#' value columns as follows:
#' \describe{
#' \item{X,Y...: }{reference population values for model parameters X,Y...}
#' \item{omega(X),omega(Y)...: }{IIV population values for model parameters X,Y,...}
#' \item{corr(X,Y),...: }{IIV correlation population values for model parameter pairs.}
#' \item{beta_X(COVRT1),beta_Y(COVRT2),...: }{Covariate coefficient population values.}
#' \item{error_PROP1,...: }{Residual error population parameters.}
#' }
#' @examples
#' \dontrun{
#' #See the examples for sampleIndParamValues_MMVmalariaXLS:
#' ?sampleIndParamValues_MMVmalariaXLS
#' }
#' @family functions in the NLME parameter sampling API
#'
#' @export
#'
#' @author Venelin Mitov (IntiQuan)
sampleParamFromUncertainty_MMVmalariaXLS <- function(spec,
                                                     Npop = if(is_IQRnlmeParamSpec(spec)) spec$Npop else 1L) {

  IQRtools:::sampleParamFromUncertainty(spec, Npop)
}


#' Sample NLME random effects given population parameter values
#' describing the inter-individual variability distribution
#'
#' @param spec A filename (character string) denoting the path to a GPF excel file, or a GPF object, or an IQRnlmeParamSpec object.
#' @param iivPopParamValues a data.frame with population parameter values
#' describing the IIV distribution. The number of rows in this data.frame
#' defines the number of populations for which random effects are to be
#' sampled.
#'
#' @return a data.frame of `Nsamples*nrow(iivPopParamValues)` rows and
#' columns `ID.POP`, `ID` and a column for each model parameter
#' defined in `spec`.
#' @examples
#' \dontrun{
#' #See the examples for sampleIndParamValues_MMVmalariaXLS:
#' ?sampleIndParamValues_MMVmalariaXLS
#' }
#' @family functions in the NLME parameter sampling API
#'
#' @export
#'
#' @author Venelin Mitov (IntiQuan)
sampleRandomEffects_MMVmalariaXLS <- function(spec,
                                              iivPopParamValues,
                                              Nsamples = if(is_IQRnlmeParamSpec(spec)) spec$Nsamples else 1) {

  IQRtools:::sampleRandomEffects(spec, iivPopParamValues, Nsamples)
}


#' Sample or subset records from a matrix or a data.frame
#'
#' @param data A data.frame describing patient or other subject covariate and regressor data. If `nrow(data)==Nsamples`, then all rows in data are returned (in their original order) in the resulting data.frame; if `nrow(data)>Nsamples`, sampling without replacement from the rows in data is done; if `nrow(data)<Nsamples` then sampling with replacement is done.
#' @param Nsamples Integer denoting the number of individuals in a population. Note: in the case of `FLAG_SAMPLE=2`, this value is used to specify the number of populations (see `Npop`). Default: `if(is.IQRnlmeParamSpec(obj)) obj$Nsamples else 1L`.
#'
#' @return a data.frame. If nrow(data)==Nsamples, then all rows in data
#' are included (in their original order); if nrow(data)>Nsamples, sampling
#' without replacement from the rows in data is done; if nrow(data)<Nsamples
#' then sampling with replacement is done. The returned data.frame has a
#' column ID set to the integers 1:Nsamples. If the original data object
#' did have a column ID, this column is copied as a column ID0 in the returned
#' data.frame.
#'
#' @examples
#' \dontrun{
#' #See the examples for sampleIndParamValues_MMVmalariaXLS:
#' ?sampleIndParamValues_MMVmalariaXLS
#' }
#' @family functions in the NLME parameter sampling API
#'
#' @export
#'
#' @author Venelin Mitov (IntiQuan)
sampleIDs_MMVmalariaXLS <- function(data, Nsamples) {
  IQRtools:::sampleIDs(data, Nsamples)
}

#' Calculate typical individual parameter values by applying covariate
#' formulae to reference population parameter values
#'
#' @param spec A filename (character string) denoting the path to a GPF excel file, or a GPF object, or an IQRnlmeParamSpec object.
#' @param referencePopParamValues a data.frame with a column ID.POP and columns
#' including at least the model parameters and the betas.
#' @param data a data.frame describing patient or other subject covariate and
#' regressor data. If nrow(data) == nrow(referencePopParamValues), then
#' the function supports two possible treatments depending on doCartesian (see
#' description for doCartesian argument below).
#' @param doCartesian a logical. When TRUE (the DEFAULT), each row in referencePopParamValues
#' is bound to data, so that the returned data.frame contains
#' nrow(referencePopParamValues)*nrow(data) rows. If FALSE
#' (only allowed when `nrow(referencePopParamValues) == nrow(data)`), then
#' the two tables are joined using [cbind()].
#' @return a data.frame with number of rows equal either to
#' `nrow(referencePopParamValues)*nrow(data)` (when `doCartesian == TRUE) (default)`),
#' or `nrow(referencePopParamValues)` (when
#' `doCartesian == FALSE && nrow(referencePopParamValues) == nrow(data)`).
#' The columns in this data.frame are ID.POP, ID and other columns matching
#' the model parameters according to spec.
#' @examples
#' \dontrun{
#' #See the examples for IQRtools::sampleIndParamValues
#' }
#' @import data.table
#'
#' @family functions in the NLME parameter sampling API
#'
#' @export
#'
#' @author Venelin Mitov (IntiQuan)
calcTypicalIndParamValues_MMVmalariaXLS <- function(spec,
                                                    referencePopParamValues,
                                                    data,
                                                    doCartesian = TRUE) {

  IQRtools:::calcTypicalIndParamValues(spec, referencePopParamValues, data, doCartesian)
}


#' Calculate individual parameter values by applying sampled random effects
#' to typical individual parameter values
#'
#' @param spec A filename (character string) denoting the path to a GPF excel file, or a GPF object, or an IQRnlmeParamSpec object.
#' @param typicalIndParamValues a data.frame specifying typical individual
#' parameter values on an original (not transformed to normal) scale for
#' the model parameters in `spec`. Apart from the
#' parameter columns, this data.frame should have an `ID.POP` and an
#' `ID` column corresponding to the columns `ID.POP` and `ID`
#' in `randomEffects`. This data.frame
#' should have the same number of rows as `randomEffects`. See
#' [calcTypicalIndParamValues()].
#' @param randomEffects a data.frame specifying random effect values on a
#' normal (transformed) scale for the model parameters in `spec`.
#' Apart from the parameter columns, this data.frame should have an `ID.POP`
#' and an `ID` integer columns corresponding to population id and individual id
#' within a population respectively. See [sampleRandomEffects()].
#'
#' @return a data.frame with the same column names as the arguments `typicalIndParamValues`
#' and `randomEffects`. The individual parameter values in this data.frame are on
#' the original scale (not transformed to a normal scale).
#' @examples
#' \dontrun{
#' #See the examples for IQRtools::sampleIndParamValues
#' }
#' @family functions in the NLME parameter sampling API
#'
#' @export
#' @examples
#' # See examples in [sampleIndParamValues_MMVmalariaXLS()].
#'
#' @author Venelin Mitov (IntiQuan)
calcIndParamValues_MMVmalariaXLS <- function(spec,
                                             typicalIndParamValues,
                                             randomEffects) {

  IQRtools:::calcIndParamValues(spec, typicalIndParamValues, randomEffects)
}


#' Sample individual parameter values, based on a NLME parameter
#' specification
#'
#' @param spec A filename (character string) denoting the path to a GPF excel file, or a GPF object, or an IQRnlmeParamSpec object.
#' @param data a data.frame describing patient or other subject covariate and regressor data.
#' If nrow(data)==Nsamples, then all rows in data
#' are returned (in their original order) in the resulting data.frame; if
#' nrow(data)>Nsamples, sampling without replacement from the rows in data
#' is done; if nrow(data)<Nsamples then sampling with replacement is done.
#' @param Nsamples integer denoting the number of individuals in a
#' population. Note: in the case of `FLAG_SAMPLE=2`, this value is used to specify
#' the number of populations (see `Npop`). Default:
#' `if(is.IQRnlmeParamSpec(obj)) obj$Nsamples else 1L`.
#' @param Npop integer denoting the number of populations to sample. Default:
#'  `if(is.IQRnlmeParamSpec(obj)) obj$Npop else if(FLAG_SAMPLE==2) Nsamples else 1L`.
#' @param FLAG_SAMPLE  Flag indicating the type of sampling:
#' \tabular{ll}{
#' FLAG_SAMPLE       \tab Meaning \cr
#' ===========       \tab =================================================== \cr
#' 0                 \tab Use point estimates of population parameters (do not consider uncertainty) and sample Nsample
#'                        individual patients based on these. Covariates considered if defined by user and used in model.
#'                        Please note: population parameters do not take covariates into account!\cr
#' 1                 \tab Sample single set of population parameters from uncertainty distribution and sample Nsample
#'                        individual patient parameters based on these. Covariates considered if defined by user and
#'                        used in model. Please note: population parameters do not take covariates into account!\cr
#' 2                 \tab Sample Nsample sets of population parameters from uncertainty distribution.
#'                        Do not sample from variability distribution and do not take into account covariates (even
#'                        if user specified).\cr
#' 3                 \tab Use point estimates of population parameters (do not consider uncertainty)
#'                        Return Nsamples sets of population parameters with covariates taken into account.\cr
#' 4                 \tab Sample single set of population parameters from uncertainty distribution
#'                        Return Nsamples sets of population parameters with covariates taken into account.\cr
#' 5                 \tab Use point estimates of population parameters (do not consider uncertainty) and sample Nsample
#'                        individual patients based on these. Covariates considered if defined by user and used in model.
#'                        Population parameters, typical individual parameters (no IIV, but considering covariates),
#'                        and individual parameters are returned.\cr
#' 6                 \tab Sample single set of population parameters from uncertainty distribution and sample Nsample
#'                        individual patient parameters based on these. Covariates considered if defined by user and used in model.
#'                        Sampled population parameters, typical individual parameters (no IIV, but considering covariates),
#'                        and individual parameters are returned.\cr
#' 7                 \tab Same as FLAG_SAMPLE=0 but obtaining the empirical random effect covariance matrix from the post-hoc ETAs.
#'                        Only used for IQRsysProject. If IQRnlmeProject then fallback to flag 0.\cr
#' 8                 \tab Same as FLAG_SAMPLE=1 but obtaining the empirical random effect covariance matrix from the post-hoc ETAs.\cr
#'                        Only used for IQRsysProject. If IQRnlmeProject then fallback to flag 1.\cr
#' }
#'
#' @details If `spec` is an IQRnlmeParamSpec object and `Npop` is specified,
#' the specified argument Npop overwrites the `spec$Npop`.
#' Default: `if(is.IQRnlmeParamSpec(spec)) spec$Npop else 1L`.
#'
#' @family functions in the NLME parameter sampling API
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' covrtData <- data.frame(
#'   ID = seq_len(4),
#'   WT0 = rlnorm(4, 3, sdlog = 0.2),
#'   SEX = sample(c(0, 1), 4, replace = TRUE))
#'
#' set.seed(1)
#' listResultsSampling <- sampleIndParamValues_MMVmalariaXLS(
#'   fileXLS,
#'   Nsamples = 10,
#'   FLAG_SAMPLE = 1,
#'   data = covrtData)
#'
#' # create a IQRnlmeParamSpec object for future use (not suitable for reading by humans):
#' spec <- specifyParamSampling(fileXLS, Nsamples = 10, FLAG_SAMPLE = 1)
#' is.IQRnlmeParamSpec(spec)
#'
#' set.seed(1)
#'
#' popParamValues<- sampleParamFromUncertainty_MMVmalariaXLS(spec)
#' # should be all 0s
#' listResultsSampling$popParamValues - popParamValues
#'
#' sampledData <- sampleIDs_MMVmalariaXLS(covrtData, spec$Nsamples)
#' # should be all 0s
#' listResultsSampling$sampledData - sampledData
#'
#' typicalIndParamValues <- calcTypicalIndParamValues_MMVmalariaXLS(spec, popParamValues, sampledData)
#' # should be all 0s
#' listResultsSampling$typicalIndParamValues - typicalIndParamValues
#'
#' randomEffects <- sampleRandomEffects_MMVmalariaXLS(spec, popParamValues)
#' # should be all 0s
#' listResultsSampling$randomEffects - randomEffects
#'
#' indParamValues <- calcIndParamValues_MMVmalariaXLS(spec, typicalIndParamValues, randomEffects)
#' # should be all 0s
#' listResultsSampling$indParamValues - indParamValues
#' }
#'
#' @export
#'
#' @author Venelin Mitov (IntiQuan)
sampleIndParamValues_MMVmalariaXLS <- function(spec,
                                               data        = NULL,
                                               Nsamples    = 1,
                                               Npop        = if(FLAG_SAMPLE==2) Nsamples else 1L,
                                               FLAG_SAMPLE = 0) {

  IQRtools::sampleIndParamValues(spec, data, Nsamples, Npop, FLAG_SAMPLE)
}

