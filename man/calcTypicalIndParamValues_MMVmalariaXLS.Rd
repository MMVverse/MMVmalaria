% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/14-MMVnlmeParamSamplingXLS.R
\name{calcTypicalIndParamValues_MMVmalariaXLS}
\alias{calcTypicalIndParamValues_MMVmalariaXLS}
\title{Calculate typical individual parameter values by applying covariate
formulae to reference population parameter values}
\usage{
calcTypicalIndParamValues_MMVmalariaXLS(
  spec,
  referencePopParamValues,
  data,
  doCartesian = TRUE
)
}
\arguments{
\item{spec}{A filename (character string) denoting the path to a GPF excel file, or a GPF object, or an IQRnlmeParamSpec object.}

\item{referencePopParamValues}{a data.frame with a column ID.POP and columns
including at least the model parameters and the betas.}

\item{data}{a data.frame describing patient or other subject covariate and
regressor data. If nrow(data) == nrow(referencePopParamValues), then
the function supports two possible treatments depending on doCartesian (see
description for doCartesian argument below).}

\item{doCartesian}{a logical. When TRUE (the DEFAULT), each row in referencePopParamValues
is bound to data, so that the returned data.frame contains
nrow(referencePopParamValues)*nrow(data) rows. If FALSE
(only allowed when \code{nrow(referencePopParamValues) == nrow(data)}), then
the two tables are joined using \code{\link{cbind}}.}
}
\value{
a data.frame with number of rows equal either to
\code{nrow(referencePopParamValues)*nrow(data)} (when \code{doCartesian == TRUE) (default)}),
or \code{nrow(referencePopParamValues)} (when
\code{doCartesian == FALSE && nrow(referencePopParamValues) == nrow(data)}).
The columns in this data.frame are ID.POP, ID and other columns matching
the model parameters according to spec.
}
\description{
Calculate typical individual parameter values by applying covariate
formulae to reference population parameter values
}
\examples{
\dontrun{
#See the examples for IQRtools::sampleIndParamValues
}
}
\seealso{
Other functions in the NLME parameter sampling API: 
\code{\link{calcIndParamValues_MMVmalariaXLS}()},
\code{\link{sampleIDs_MMVmalariaXLS}()},
\code{\link{sampleIndParamValues_MMVmalariaXLS}()},
\code{\link{sampleParamFromUncertainty_MMVmalariaXLS}()},
\code{\link{sampleRandomEffects_MMVmalariaXLS}()}
}
\author{
Venelin Mitov (IntiQuan)
}
\concept{functions in the NLME parameter sampling API}
