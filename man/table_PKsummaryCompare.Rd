% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/13-Reporting.R
\name{table_PKsummaryCompare}
\alias{table_PKsummaryCompare}
\title{table_PKsummaryCompare}
\usage{
table_PKsummaryCompare(
  tableInputPath1,
  tableInputPath2,
  tableOutputPath = NULL,
  renameMODEL = NULL,
  addDescription = NULL,
  Metric = "BIC",
  xTitle = NULL,
  xFooter = NULL,
  report = TRUE
)
}
\arguments{
\item{tableOutputPath}{(Default: \code{NULL}).}

\item{renameMODEL}{(Default: \code{NULL}).}

\item{addDescription}{(Default: \code{NULL}).}

\item{Metric}{Default: 'BIC'}

\item{xTitle}{(Default: \code{NULL}).}

\item{xFooter}{(Default: \code{NULL}).}

\item{report}{Default: TRUE}
}
\description{
table_PKsummaryCompare
}
\seealso{
Other Reporting: 
\code{\link{generate_EstimateTableForReport}()},
\code{\link{parameter_strToLaTeX}()},
\code{\link{reportText_PKmodelHypos}()},
\code{\link{reportText_convertListToText}()},
\code{\link{table_PKsummaryResults}()},
\code{\link{table_covariateSummaryResults}()},
\code{\link{table_strToLaTeX}()}
}
\author{
Mohammed H. Cherkaoui (MMV, \email{cherkaouim@mmv.org})
}
\concept{Reporting}
