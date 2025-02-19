% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07-ModelAssessment.R
\name{likelihoodRatioTestIQR}
\alias{likelihoodRatioTestIQR}
\title{likelihoodRatioTestIQR}
\usage{
likelihoodRatioTestIQR(
  ModelPaths,
  ModelName = NULL,
  outputFolder = "../04-Output/Sxxx-ModelComparison",
  parameterPNG = NULL
)
}
\arguments{
\item{ModelName}{Default: \code{NULL}}

\item{outputFolder}{Default: '../04-Output/Sxxx-ModelComparison'}

\item{parameterPNG}{Default: \code{NULL}}
}
\description{
likelihoodRatioTestIQR
}
\seealso{
Other Model Assessment: 
\code{\link{add_CovariateToModelSpec}()},
\code{\link{add_OneCovariateToModelSpec}()},
\code{\link{assess_PDcomboModel}()},
\code{\link{backwardElimination}()},
\code{\link{compare_DataModelMMV}()},
\code{\link{compare_DataPopPred}()},
\code{\link{compare_IndPredFits}()},
\code{\link{compare_ModelEstimateRobustness}()},
\code{\link{compare_PredFitsSYS}()},
\code{\link{forwardSelection}()},
\code{\link{plotSummaryComparisonBinary}()},
\code{\link{plotSummaryComparisonContinuous}()},
\code{\link{plotVPC_MMVdataVPC}()},
\code{\link{remove_OneCovariateToModelSpec}()},
\code{\link{tablePopPDmalaria}()},
\code{\link{table_EstimatesComboNLME}()}
}
\author{
Aline Fuchs (MMV), Anne Kümmel (IntiQuan), Mohammed H. Cherkaoui (MMV)
}
\concept{Model Assessment}
