% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07-ModelAssessment.R
\name{compare_ModelEstimateRobustness}
\alias{compare_ModelEstimateRobustness}
\title{compare_ModelEstimateRobustness}
\usage{
compare_ModelEstimateRobustness(
  FitList,
  filename,
  ModelAlignment = "Vertical",
  PDddiPara = c("Alpha", "Beta", "Alpha12", "Alpha21", "Beta12", "Beta21", "Gamma")
)
}
\arguments{
\item{ModelAlignment}{Default: 'Vertical'}

\item{PDddiPara}{Default: c("Alpha", "Beta", "Alpha12", "Alpha21", "Beta12", "Beta21",
"Gamma")}
}
\description{
compare_ModelEstimateRobustness
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
\code{\link{compare_PredFitsSYS}()},
\code{\link{forwardSelection}()},
\code{\link{likelihoodRatioTestIQR}()},
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
